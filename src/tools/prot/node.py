# src/tools/prot/node.py
# Author: Dmitri Kosenkov
# Created: 2025-08-25
# Updated: 2025-09-18
"""
node.py
=======

Main entry point for the protein visualization and summarization tool.
Refactored from standalone tool_prot.py into the src.tools.prot package.

Responsibilities:
-----------------
1. Resolve gene symbols to UniProt IDs and local protein records
2. Fetch per-residue pLDDT, fpocket, SASA, and polarity index features
3. Assemble summaries, warnings, and interpretation notes
4. Build a standalone interactive HTML viewer (3Dmol.js)
5. Return updated State and write HTML/TXT outputs to OUTPUT_DIR
"""

import sys
from pathlib import Path
import warnings
from typing import Dict, List

PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import src.tools.prot as prot  # loads constants from  __init__.py

# --- Imports ---
from src.state import State
from src.config import DEBUG, TOOL_PROT_MAX_GENES
from src.alvessa.domain.gene_class import Gene
from src.alvessa.domain.gene_components import GeneIdentifiers
from src.tools.base import Node

# --- Local storage Layout ---
from src.tools.prot import (
    REPO_ROOT,
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
    make_summary_text,
    interpretation_notes,
    derive_warning_flags,
    append_flags_to_summary_text,
    inject_frontend_assets,
)

from src.tools.prot.plddt import fetch_plddt
from src.tools.prot.fpocket import fetch_fpocket
from src.tools.prot.sasa_pi import fetch_sasa_pi

def prot_agent(state: "State") -> "State":
    """Main agent pipeline: resolves UniProt IDs, fetches per-residue features,
    generates summaries, and builds interactive HTML for the frontend.
    """
    ensure_fresh_db()
    genes = state.get("genes") or []
    if not genes and state.get("gene_symbol"):
        genes = [state["gene_symbol"]]
    if not genes:
        warnings.warn("No gene symbols provided.")
        state["used_tools"] = state.get("used_tools", []) + ["prot"]
        return state
    if len(genes) > TOOL_PROT_MAX_GENES:
        warnings.warn(f"Truncating gene list from {len(genes)} to {TOOL_PROT_MAX_GENES}")
        genes = genes[:TOOL_PROT_MAX_GENES]

    gene_entities = state.get("gene_entities") or {}
    state["gene_entities"] = gene_entities
    for g in genes:
        if not isinstance(gene_entities.get(g), Gene):
            gene_entities[g] = Gene(symbol=g)

    conn = get_connection()
    prot_data_all: Dict[str, Dict] = {}
    include_plddt_any = include_fpocket_any = include_sasa_any = include_pi_any = False

    for gene_symbol in genes:
        log(f"Processing gene: {gene_symbol}")

        uniprot_id, protein_id, pdb_file = _resolve_gene_to_protein(conn, gene_symbol)
        if not (uniprot_id and protein_id and pdb_file):
            log(f"Skipping {gene_symbol}: could not resolve UniProt/Protein/PDB mapping")
            continue

        plddt_stats, plddt_norm = fetch_plddt(conn, protein_id)
        fpocket_stats, fpocket_norm = fetch_fpocket(conn, protein_id)
        sasa_stats, sasa_norm, pi_stats, pi_norm, residue_labels = fetch_sasa_pi(conn, protein_id)
        pdb_data = load_pdb_inline(pdb_file)
        if not pdb_data:
            log(f"Skipping {gene_symbol}: missing PDB data ({pdb_file})")
            continue

        include_plddt_any |= bool(plddt_stats)
        include_fpocket_any |= bool(fpocket_stats)
        include_sasa_any |= bool(sasa_stats)
        include_pi_any |= bool(pi_stats)

        summary = make_summary_text(
            gene_symbol, uniprot_id, protein_id, pdb_file,
            plddt_stats, fpocket_stats, sasa_stats, pi_stats
        )
        flags = derive_warning_flags(plddt_stats)
        if flags:
            summary = append_flags_to_summary_text(summary, flags)
        gene_entities[gene_symbol].update_text_summaries(
            "Protein structure and druggability:\n" + summary
        )

        prot_data_all[gene_symbol] = {
            "plddt": plddt_norm or [],
            "fpocket": fpocket_norm or [],
            "sasa": sasa_norm or [],
            "pi": pi_norm or [],
            "pdb": pdb_data,
            "stats": {
                "plddt": plddt_stats,
                "fpocket": fpocket_stats,
                "sasa": sasa_stats,
                "pi": pi_stats,
            },
            "labels": residue_labels or {},
        }

    conn.close()

    if prot_data_all:
        notes = interpretation_notes(
            include_fpocket_any, include_sasa_any, include_pi_any, include_plddt_any
        )
        last_gene = genes[-1]
        if gene_entities.get(last_gene):
            gene_entities[last_gene].update_text_summaries(notes)
    else:
        log("No proteins resolved for any of the requested genes.")

    with open(HTML_TEMPLATE, "r", encoding="utf-8") as f:
        html_template = f.read()
    with open(CSS_TEMPLATE, "r", encoding="utf-8") as f:
        css_template = f.read()
    with open(JS_TEMPLATE, "r", encoding="utf-8") as f:
        js_template = f.read()

    html = inject_frontend_assets(html_template, css_template, js_template, prot_data_all, ", ".join(genes))
    state["used_tools"] = state.get("used_tools", []) + ["prot"]
    if prot_data_all:
        state["prot_html"] = html
    return state


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

    cur = conn.cursor()
    cur.execute(
        "SELECT protein_id, pdb_file FROM data_proteins WHERE uniprot_id=? LIMIT 1",
        (uniprot_id,),
    )
    row = cur.fetchone()
    if not row:
        return uniprot_id, None, None
    protein_id, pdb_file = row
    return uniprot_id, protein_id, pdb_file


# ----------------------------------------------------------------------
# CLI (testing)
# ----------------------------------------------------------------------
if __name__ == "__main__":
    genes = sys.argv[1:] if len(sys.argv) > 1 else ["TP53", "EGFR", "AFM"]
    base_name = genes[0] if len(genes) == 1 else f"{genes[0]}_plus{len(genes)-1}"

    state = State({"genes": genes, "gene_entities": {g: Gene(GeneIdentifiers(symbol=g)) for g in genes}})
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

    print("[OK] Generated outputs:")
    print(f"  * HTML file: {html_out.resolve()}")
    print(f"  * TXT file:  {txt_out.resolve()}")
    print(f"[INFO] Genes processed: {', '.join(genes[:TOOL_PROT_MAX_GENES])}")


# --- Node Registration ---
NODES: tuple[Node, ...] = (
    Node(
        name="prot",
        entry_point=prot_agent,
        description = (
        "Visualize AlphaFold-predicted protein structures for one or more genes, overlaying "
        "per-residue pLDDT confidence, FPocket-derived druggability scores, solvent-accessible "
        "surface area (SASA), and polarity index (hydrophilic vs. hydrophobic regions). "
        "Generates interactive 3D views with color-coded surfaces, summary statistics "
        "(min/max/average per protein), and interpretation notes for druggability assessment.")
    ),
)
