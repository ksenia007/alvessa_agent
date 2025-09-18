"""
Author: Dmitri Kosenkov
Created: 2025-08-25
Updated: 2025-09-10

Description:
Agentic-style protein visualization and summarization tool supporting multiple
human genes. For each provided gene symbol, the tool resolves the corresponding
UniProt entry and retrieves structural confidence (pLDDT) and fpocket-based
druggability data from the local SQLite database alvessa_proteins.db, along
with the corresponding AlphaFold PDB structure.

Features:
  - Generates an HTML page with an interactive 3Dmol.js viewer:
    * Dropdown to select which protein to visualize
    * Style options: cartoon, surface colored by pLDDT, or surface colored
      by fpocket druggability scores
    * Legends display raw pLDDT (0-100 percent) or fpocket min-max values
      before normalization, while visualization uses normalized scores
  - Updates each Gene entity with a concise text summary reporting UniProt ID,
    protein ID, and raw statistics (min, max, avg) for pLDDT and fpocket
  - Appends reusable interpretation notes (pLDDT reliability categories,
    fpocket druggability score interpretation) once at the end of the output

Returns:
  - Updated State object with:
    * HTML string (prot_html) for interactive visualization
    * Updated Gene entities containing protein summaries and interpretation notes
"""

import os
import sys
import json
import sqlite3
import warnings
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, Tuple

import requests

from src.state import State
from src.config import DEBUG, TOOL_PROT_MAX_GENES
from src.alvessa.domain.gene_class import Gene
from src.alvessa.domain.gene_components import GeneIdentifiers
from src.tools.base import Node

# --- Local storage layout ---
PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]
LOCAL_DBS_DIR = REPO_ROOT / "local_dbs"
STATIC_DIR = REPO_ROOT / "static"
DB_PATH = LOCAL_DBS_DIR / "alvessa_proteins.db"
PDB_DIR = LOCAL_DBS_DIR / "pdb"
HTML_TEMPLATE = STATIC_DIR / "tool_prot.html"
CSS_TEMPLATE = STATIC_DIR / "tool_prot.css"
JS_TEMPLATE = STATIC_DIR / "tool_prot.js"

# ----------------------------------------------------------------------
# Utilities
# ----------------------------------------------------------------------
def _log(msg: str):
    if DEBUG:
        print(f"[Protein] {msg} @ {datetime.now()}")


def _get_connection():
    return sqlite3.connect(str(DB_PATH), check_same_thread=False)


def _mark_used(state, name="prot"):
    used = list(state.get("used_tools", []))
    if name not in used:
        used.append(name)
    return used


# ----------------------------------------------------------------------
# UniProt resolver
# ----------------------------------------------------------------------
def get_uniprot_entry_for_gene(gene_symbol: str) -> Optional[Dict]:
    """Retrieve the reviewed UniProtKB record for a human gene symbol."""
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene:{gene_symbol} AND organism_id:9606 AND reviewed:true",
        "format": "json",
        "size": 1,
    }
    try:
        resp = requests.get(base_url, params=params, timeout=12)
        resp.raise_for_status()
        results = resp.json().get("results", [])
        return results[0] if results else None
    except Exception as exc:
        warnings.warn(f"Error fetching UniProt entry for {gene_symbol}: {exc}")
        return None


# ----------------------------------------------------------------------
# Database accessors
# ----------------------------------------------------------------------
def _resolve_protein_id(conn, uniprot_id: str) -> Optional[Tuple[str, str]]:
    cur = conn.cursor()
    cur.execute(
        """
        SELECT protein_id, pdb_file
        FROM data_proteins
        WHERE uniprot_id = ?
        LIMIT 1
        """,
        (uniprot_id,),
    )
    return cur.fetchone()


def _fetch_plddt(conn, protein_id: str):
    cur = conn.cursor()
    cur.execute(
        """
        SELECT residue_no, plddt
        FROM data_plldt_residues
        WHERE protein_id = ?
        ORDER BY residue_no ASC
        """,
        (protein_id,),
    )
    rows = cur.fetchall()
    if not rows:
        return None, None

    values = [r[1] for r in rows if r[1] is not None]
    if not values:
        return None, None

    stats = {"min": min(values), "max": max(values), "avg": sum(values) / len(values)}
    normalized = [{"residue_no": r[0], "score": r[1] / 100.0} for r in rows if r[1] is not None]
    return stats, normalized


def _fetch_fpocket(conn, protein_id: str):
    cur = conn.cursor()
    cur.execute(
        """
        SELECT residue_no, druggability_score
        FROM data_fpocket_residues_raw
        WHERE protein_id = ?
        ORDER BY residue_no ASC
        """,
        (protein_id,),
    )
    rows = cur.fetchall()
    if not rows:
        return None, None

    values = [r[1] for r in rows if r[1] is not None]
    if not values:
        return None, None

    stats = {
        "min": min(values),
        "max": max(values),
        "avg": sum(values) / len(values),
    }

    min_val, max_val = stats["min"], stats["max"]
    rng = max_val - min_val if max_val > min_val else 1.0

    # --- Visualization normalization: min : 0.0, max : 1.0
    normalized = [
        {"residue_no": r[0], "score": (r[1] - min_val) / rng}
        for r in rows if r[1] is not None
    ]

    # Return both raw stats (for summaries) and normalized (for visualization)
    return stats, normalized


def _load_pdb_inline(protein_id: str, pdb_file: str):
    pdb_path = PDB_DIR / pdb_file
    if not pdb_path.exists():
        warnings.warn(f"PDB not found for {protein_id}: {pdb_file}")
        return None
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        return f.read()


# ----------------------------------------------------------------------
# Summary preparation
# ----------------------------------------------------------------------
def _make_summary(
    gene: str,
    uniprot_id: str,
    protein_id: str,
    pdb_file: str,
    plddt_stats: Optional[Dict],
    fpocket_stats: Optional[Dict],
) -> str:
    """Return concise per-protein summary for reasoning (ASCII only)."""
    lines = [
        f"Gene: {gene}",
        f"  UniProt ID: {uniprot_id}",
        f"  Protein ID: {protein_id}",
    ]
    if plddt_stats:
        lines.append(
            f"  pLDDT: min={plddt_stats['min']:.2f}, max={plddt_stats['max']:.2f}, avg={plddt_stats['avg']:.2f}"
        )
    else:
        lines.append("  pLDDT: not available")
    if fpocket_stats:
        lines.append(
            f"  FPocket: min={fpocket_stats['min']:.3f}, max={fpocket_stats['max']:.3f}, avg={fpocket_stats['avg']:.3f}"
        )
    else:
        lines.append("  FPocket: not available")
    return "\n".join(lines)


def _interpretation_notes() -> str:
    """Reusable block of interpretation notes for pLDDT and FPocket."""
    return "\n".join([
        "",
        "Interpretation Notes",
        "",
        "pLDDT (Predicted Local Distance Difference Test)",
        "Per-residue confidence score from AlphaFold, scaled 0-100:",
        "  - >90 : very high reliability (backbone + side chains accurate)",
        "  - 70-90 : backbone usually correct; side chains less certain",
        "  - <70 : lower confidence, often flexible/disordered regions",
        "Refs: Mariani V. et al., 2013; Guo C. et al., 2022; EBI AlphaFold course",
        "",
        "FPocket Druggability Score",
        "Numerical score (0-1) estimating drug-likeness of a pocket:",
        "  - 0   : very unlikely to bind",
        "  - ~0.5: borderline",
        "  - 1   : highly druggable",
        "Refs: Schmidtke P. and Barril X., J. Med. Chem. 2010",
    ])


# ----------------------------------------------------------------------
# Agentic Node
# ----------------------------------------------------------------------
def prot_agent(state: "State") -> "State":
    genes = state.get("genes") or []
    if not genes and state.get("gene_symbol"):
        genes = [state["gene_symbol"]]

    if not genes:
        warnings.warn("No gene symbols provided.")
        return {**state, "used_tools": _mark_used(state)}

    omitted = 0
    if len(genes) > TOOL_PROT_MAX_GENES:
        warnings.warn(
            f"Truncating gene list from {len(genes)} to {TOOL_PROT_MAX_GENES}"
        )
        omitted = len(genes) - TOOL_PROT_MAX_GENES

    genes = genes[:TOOL_PROT_MAX_GENES]

    # Ensure gene_entities exists
    gene_entities = state.get("gene_entities")
    if not isinstance(gene_entities, dict):
        gene_entities = {}
        state["gene_entities"] = gene_entities
    for g in genes:
        if not isinstance(gene_entities.get(g), Gene):
            gene_entities[g] = Gene(GeneIdentifiers(symbol=g))

    conn = _get_connection()
    prot_data_all: Dict[str, Dict[str, object]] = {}

    summaries_for_all_genes = []

    for gene_symbol in genes:
        entry = get_uniprot_entry_for_gene(gene_symbol)
        if not entry:
            continue

        uniprot_id = entry.get("primaryAccession")
        row = _resolve_protein_id(conn, uniprot_id) if uniprot_id else None
        if not row:
            continue

        protein_id, pdb_file = row
        plddt_stats, plddt_norm = _fetch_plddt(conn, protein_id)
        fpocket_stats, fpocket_norm = _fetch_fpocket(conn, protein_id)
        pdb_data = _load_pdb_inline(protein_id, pdb_file)
        if not pdb_data:
            continue

        # Per-gene summary
        summary = _make_summary(
            gene_symbol, uniprot_id, protein_id, pdb_file, plddt_stats, fpocket_stats
        )
        summaries_for_all_genes.append(summary)

        gene_entities[gene_symbol].update_text_summaries(
            f"Protein structure and druggability: {summary}"
        )

        prot_data_all[gene_symbol] = {
            "plddt": plddt_norm or [],
            "fpocket": fpocket_norm or [],
            "pdb": pdb_data,
            "stats": {
                "plddt": plddt_stats,
                "fpocket": fpocket_stats,
            },
        }

    conn.close()

    # Append interpretation notes ONCE, to the last gene only
    if summaries_for_all_genes:
        notes = _interpretation_notes()
        last_gene = genes[-1]
        gene_obj = gene_entities.get(last_gene)
        if gene_obj:
            gene_obj.update_text_summaries(notes)

    data_script = f"""
    <script>
      const protData = {json.dumps(prot_data_all)};
    </script>
    """

    with open(HTML_TEMPLATE, "r", encoding="utf-8") as f:
        html_template = f.read()
    with open(CSS_TEMPLATE, "r", encoding="utf-8") as f:
        css_template = f.read()
    with open(JS_TEMPLATE, "r", encoding="utf-8") as f:
        js_template = f.read()

    html = (
        html_template.replace("{{GENE_SYMBOL}}", ", ".join(genes))
        .replace("{{CSS_INLINE}}", f"<style>\n{css_template}\n</style>")
        .replace("{{JS_INLINE}}", f"<script>\n{js_template}\n</script>")
        .replace("</body>", data_script + "\n</body>")
    )

    result = {**state, "used_tools": _mark_used(state)}
    if prot_data_all:  # only attach HTML if we actually built data
        result["prot_html"] = html
    return result


# ----------------------------------------------------------------------
# CLI (testing)
# ----------------------------------------------------------------------
if __name__ == "__main__":
    genes = sys.argv[1:] if len(sys.argv) > 1 else ["TP53", "EGFR"]
    base_name = genes[0] if len(genes) == 1 else f"{genes[0]}_plus{len(genes)-1}"

    state = State({
        "genes": genes,
        "gene_entities": {g: Gene(GeneIdentifiers(symbol=g)) for g in genes}
    })

    result = prot_agent(state)

    html_out = f"{base_name}_prot.html"
    with open(html_out, "w", encoding="utf-8") as f:
        f.write(result.get("prot_html", "<p>No HTML produced.</p>"))

    txt_out = f"{base_name}_prot.txt"
    lines = []

    processed_genes = genes[:TOOL_PROT_MAX_GENES]
    gene_entities = result.get("gene_entities", {})

    for g in processed_genes:
        obj = gene_entities.get(g)
        if not obj:
            continue
        if getattr(obj, "text_summaries_from_tools", None):
            lines.append(g)
            lines.extend(obj.text_summaries_from_tools)
            lines.append("")

    with open(txt_out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) if lines else "No summaries produced.")

    print(f"[OK] Generated outputs for {', '.join(processed_genes)}: {html_out}, {txt_out}")


NODES: tuple[Node, ...] = (
    Node(
        name="prot",
        entry_point=prot_agent,
        description="Retrieves structural data for a single or several proteins given one or several gene symbols. Resolves UniProt ID and AlphaFold Protein Database structure, and provides per-residue metrics: pLDDT confidence scores (structural reliability) and FPocket druggability scores. Outputs both an interactive 3Dmol.js visualization and a text summary of min, max and averaged metrics.",
    ),
)
