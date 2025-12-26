# src/tools/prot/node.py
# Author: Dmitri Kosenkov
# Created: 2025-09-20
# Updated: 2025-12-08
"""
node.py
=======

Main entry point for the protein visualization and summarization tool.
Refactored from standalone tool_prot.py into the src.tools.prot package.

Responsibilities:
-----------------
1. Resolve gene symbols to UniProt IDs and local protein records
2. Fetch per-residue pLDDT, fpocket, SASA, polarity index, disorder consensus, MoRF propensity, and BioLiP2 evidence
3. Assemble summaries, warnings, and interpretation notes
4. Build a standalone interactive HTML viewer (3Dmol.js)
5. Return updated State and write HTML/TXT outputs to OUTPUT_DIR
"""

import sys
import warnings
from pathlib import Path
from typing import Dict, List

PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import src.tools.prot as prot  # loads constants from __init__.py

# --- Imports ---
from src.state import State
from src.alvessa.domain.gene_class import Gene
from src.alvessa.domain.gene_components import GeneIdentifiers
from src.tools.base import Node

# --- Local storage Layout ---
from src.tools.prot import (
    TOOL_PROT_MAX_GENES,
    OUTPUT_DIR,
    HTML_TEMPLATE,
    CSS_TEMPLATE,
    JS_TEMPLATE,
)

from src.tools.prot.utils import (
    ensure_fresh_db,
    get_connection,
    load_pdb_inline,
    log,
    make_full_summary,
    interpretation_notes,
    derive_warning_flags,
    append_flags_to_summary_text,
    inject_frontend_assets,
)

from src.tools.prot.plddt import fetch_plddt
from src.tools.prot.fpocket import fetch_fpocket
from src.tools.prot.sasa_pi import fetch_sasa_pi
from src.tools.prot.disorder import fetch_disorder
from src.tools.prot.biolip2 import fetch_biolip2
from src.tools.prot.cysdb import fetch_cysdb


# ----------------------------------------------------------------------
# Gene Preparation
# ----------------------------------------------------------------------
def _prepare_genes(state: "State") -> List[str]:
    """Resolve, deduplicate, and prepare gene symbols for downstream processing."""
    genes = state.get("genes") or []
    if not genes and state.get("gene_symbol"):
        genes = [state["gene_symbol"]]

    if not genes:
        warnings.warn("No gene symbols provided.")
        state["used_tools"] = state.get("used_tools", []) + ["prot"]
        return []

    # Deduplicate while preserving input order
    seen = set()
    unique_genes = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            unique_genes.append(g)
    genes = unique_genes

    if len(genes) > TOOL_PROT_MAX_GENES:
        warnings.warn(f"Truncating gene list from {len(genes)} to {TOOL_PROT_MAX_GENES}")
        genes = genes[:TOOL_PROT_MAX_GENES]

    # Ensure gene_entities dict exists and contains Gene objects (corrected)
    gene_entities = state.get("gene_entities") or {}
    for g in genes:
        if not isinstance(gene_entities.get(g), Gene):
            gene_entities[g] = Gene(identifiers=GeneIdentifiers(symbol=g))
    state["gene_entities"] = gene_entities

    return genes


# ----------------------------------------------------------------------
# Main Protein Agent
# ----------------------------------------------------------------------
def prot_agent(state: "State") -> "State":
    """Main agent pipeline: resolves UniProt IDs, fetches per-residue features,
       generates summaries, and builds interactive HTML for the frontend.
    """
    ensure_fresh_db()

    # --- Gene handling via helper ---
    genes = _prepare_genes(state)
    if not genes:
        return state

    conn = get_connection()
    prot_data_all: Dict[str, Dict] = {}
    include_plddt_any = include_fpocket_any = include_sasa_any = include_pi_any = False
    include_disorder_any = include_morf_any = False
    include_biolip2_any = False
    include_cysdb_any = False

    for gene_symbol in genes:
        log(f"Processing gene: {gene_symbol}")

        uniprot_id, protein_id, pdb_file = _resolve_gene_to_protein(conn, gene_symbol)
        if not (uniprot_id and protein_id and pdb_file):
            log(f"Skipping {gene_symbol}: could not resolve UniProt/Protein/PDB mapping")
            continue

        plddt_stats, plddt_norm = fetch_plddt(conn, protein_id)
        fpocket_stats, fpocket_norm = fetch_fpocket(conn, protein_id)
        sasa_stats, sasa_norm, pi_stats, pi_norm, residue_labels = fetch_sasa_pi(conn, protein_id)
        disorder_stats, disorder_norm, morf_stats, morf_norm = fetch_disorder(conn, protein_id, uniprot_id)
        biolip2_summary, biolip2_norm = fetch_biolip2(conn, uniprot_id)
        cysdb_stats, cysdb_tracks = fetch_cysdb(conn, protein_id, uniprot_id)

        pdb_data = load_pdb_inline(pdb_file)
        if not pdb_data:
            log(f"Skipping {gene_symbol}: missing PDB data ({pdb_file})")
            continue

        include_plddt_any |= bool(plddt_stats)
        include_fpocket_any |= bool(fpocket_stats)
        include_sasa_any |= bool(sasa_stats)
        include_pi_any |= bool(pi_stats)
        include_disorder_any |= bool(disorder_stats)
        include_morf_any |= bool(morf_stats)
        include_biolip2_any |= bool(biolip2_summary)
        if cysdb_stats and cysdb_stats.get("include_cysdb"):
            include_cysdb_any = True

        # --- Text summary ---
        summary = make_full_summary(
            gene_symbol,
            uniprot_id,
            protein_id,
            pdb_file,
            plddt_stats,
            fpocket_stats,
            sasa_stats,
            pi_stats,
            disorder_stats,
            morf_stats,
            biolip2_summary,
            cysdb_stats,
        )

        flags = derive_warning_flags(plddt_stats)
        if flags:
            summary = append_flags_to_summary_text(summary, flags)

        state["gene_entities"][gene_symbol].update_text_summaries(
            "*Protein: " + summary.replace("\n", " ").strip()
        )

        # --- Store all features for frontend ---
        prot_data_all[gene_symbol] = {
            "plddt": plddt_norm or [],
            "fpocket": fpocket_norm or [],
            "sasa": sasa_norm or [],
            "pi": pi_norm or [],
            "disorder": disorder_norm or [],
            "morf": morf_norm or [],
            "pdb": pdb_data,
            "biolip2": biolip2_norm or [],
            "cysdb_hyperreactive": cysdb_tracks.get("cysdb_hyperreactive", []) if cysdb_tracks else [],
            "cysdb_ligandable": cysdb_tracks.get("cysdb_ligandable", []) if cysdb_tracks else [],
            "cysdb_is_act_site": cysdb_tracks.get("cysdb_is_act_site", []) if cysdb_tracks else [],
            "cysdb_near_act_site": cysdb_tracks.get("cysdb_near_act_site", []) if cysdb_tracks else [],
            "cysdb_is_bind_site": cysdb_tracks.get("cysdb_is_bind_site", []) if cysdb_tracks else [],
            "cysdb_near_bind_site": cysdb_tracks.get("cysdb_near_bind_site", []) if cysdb_tracks else [],
            "stats": {
                "plddt": plddt_stats,
                "fpocket": fpocket_stats,
                "sasa": sasa_stats,
                "pi": pi_stats,
                "disorder": disorder_stats,
                "morf": morf_stats,
                "biolip2": biolip2_summary,
                "cysdb": cysdb_stats,
            },
            "labels": residue_labels or {},
        }

    conn.close()

    if prot_data_all:
        notes = interpretation_notes(
            include_fpocket_any,
            include_sasa_any,
            include_pi_any,
            include_plddt_any,
            include_disorder_any,
            include_morf_any,
            include_biolip2_any,
            include_cysdb_any,
        )
        last_gene = genes[-1]
        if state["gene_entities"].get(last_gene):
            state["gene_entities"][last_gene].update_text_summaries("*Protein: " + notes.replace("\n", " ").strip())
    else:
        log("No proteins resolved for any of the requested genes.")

    with open(HTML_TEMPLATE, "r", encoding="utf-8") as f:
        html_template = f.read()
    with open(CSS_TEMPLATE, "r", encoding="utf-8") as f:
        css_template = f.read()
    with open(JS_TEMPLATE, "r", encoding="utf-8") as f:
        js_template = f.read()

    html = inject_frontend_assets(
        html_template, css_template, js_template, prot_data_all, ", ".join(genes)
    )
    state["used_tools"] = state.get("used_tools", []) + ["prot"]
    if prot_data_all:
        state["prot_html"] = html
        state["prot_data_all"] = prot_data_all
    return state


# ----------------------------------------------------------------------
# Helper: Resolve Gene -> Protein Mapping
# ----------------------------------------------------------------------
def _resolve_gene_to_protein(conn, gene_symbol: str):
    """Resolve a gene symbol to (uniprot_id, protein_id, pdb_file)."""
    import requests
    try:
        resp = requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={
                "query": f"gene:{gene_symbol} AND organism_id:9606 AND reviewed:true",
                "format": "json",
                "size": 1,
            },
            timeout=12,
        )
        resp.raise_for_status()
        results = resp.json().get("results", [])
        if not results:
            return None, None, None
        uniprot_id = results[0].get("primaryAccession")
    except Exception as e:
        warnings.warn(f"Error fetching UniProt entry for {gene_symbol}: {e}")
        return None, None, None

    row = conn.execute(
        "SELECT protein_id, pdb_file FROM data_proteins WHERE uniprot_id=? LIMIT 1",
        (uniprot_id,),
    ).fetchone()
    if not row:
        return uniprot_id, None, None
    protein_id, pdb_file = row
    return uniprot_id, protein_id, pdb_file


# ----------------------------------------------------------------------
# CLI (testing)
# ----------------------------------------------------------------------
if __name__ == "__main__":
    genes = sys.argv[1:] if len(sys.argv) > 1 else ["TP53", "EGFR", "AFM", "GAPDH"]
    base_name = genes[0] if len(genes) == 1 else f"{genes[0]}_plus{len(genes)-1}"

    state = State({
        "genes": genes,
        "gene_entities": {g: Gene(identifiers=GeneIdentifiers(symbol=g)) for g in genes},
    })
    result = prot_agent(state)

    # --- Outputs ---
    html_out = OUTPUT_DIR / f"{base_name}_prot.html"
    with html_out.open("w", encoding="utf-8") as f:
        f.write(result.get("prot_html", "<p>No HTML produced.</p>"))

    txt_out = OUTPUT_DIR / f"{base_name}_prot.txt"
    lines: List[str] = []
    for g in genes[:TOOL_PROT_MAX_GENES]:
        obj = result.get("gene_entities", {}).get(g)
        if obj and getattr(obj, "text_summaries_from_tools", None):
            lines.append(g)
            lines.extend(obj.text_summaries_from_tools)
            lines.append("")
    with txt_out.open("w", encoding="utf-8") as f:
        f.write("\n".join(lines) if lines else "No summaries produced.")

    # Also dump prot_data_all for debugging
    json_out = OUTPUT_DIR / f"{base_name}_prot.json"
    with json_out.open("w", encoding="utf-8") as f:
        import json
        f.write(json.dumps(result.get("prot_data_all", {}), indent=2, ensure_ascii=True))

    print("[OK] Generated outputs:")
    print(f"  * HTML file: {html_out.resolve()}")
    print(f"  * TXT file:  {txt_out.resolve()}")
    print(f"  * JSON file: {json_out.resolve()}")
    print(f"[INFO] Genes processed: {', '.join(genes[:TOOL_PROT_MAX_GENES])}")


# ----------------------------------------------------------------------
# Node Registration
# ----------------------------------------------------------------------
NODES: tuple[Node, ...] = (
    Node(
        name="prot",
        entry_point=prot_agent,
        description=(
            "Visualize AlphaFold-predicted protein structures for one or more genes. "
            "Overlays include: per-residue confidence scores (pLDDT, AlphaFold); "
            "pocket druggability and geometric features (FPocket); "
            "solvent-accessible surface area and polarity index (FreeSASA); "
            "intrinsic disorder consensus (IUPred3 + DisProt); "
            "MoRF propensity predictions (IUPred3/ANCHOR2); "
            "ligand-binding evidence (BioLiP2 + ChEMBL); "
            "and cysteine chemoproteomics annotations from CysDB, including detected, "
            "ligandable, hyperreactive, and active-site or binding-site–proximal cysteines. "
            "Generates interactive 3Dmol.js views with color-coded surfaces, binding-site "
            "context, summary statistics, and interpretation notes."
        ),
    ),
)
