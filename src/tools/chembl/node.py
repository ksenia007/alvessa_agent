# src/tools/chembl/node.py
# ===========================================================
# ChEMBL Drug–Target Summarization and Visualization Tool
# ===========================================================
#
# Author: Dmitri Kosenkov
# Created: 2025-09-25
# Updated: 2025-10-07

"""
Overview
--------
This module provides a complete, self-contained system for extracting,
summarizing, and interactively visualizing drug–target relationships
from a local ChEMBL database snapshot (schema v35+). It mirrors the
design and architecture of the `prot` tool and integrates seamlessly
into the Alvessa agentic workflow.

The tool supports both textual and interactive outputs, offering a
concise summary of all ChEMBL evidence related to one or more human
genes, including FDA-approved drugs, clinical/preclinical candidates,
and bioactivity assay data. An interactive HTML viewer renders 2D
molecular structures using RDKit.js and supports collapsible sections
for each evidence type.
"""

import sys
import warnings
from pathlib import Path
from typing import List, Dict

PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.state import State
from src.alvessa.domain.gene_class import Gene
from src.alvessa.domain.gene_components import GeneIdentifiers
from src.tools.base import Node

from src.tools.chembl.utils import (
    get_connection,
    fetch_target_info,        # now includes InChIKeys and molecule types
    make_summary_text,
    interpretation_notes,
    log,
    inject_frontend_assets,   # chembl-local version (mirrors prot tool)
)

from src.tools.chembl import OUTPUT_DIR, HTML_TEMPLATE, CSS_TEMPLATE, JS_TEMPLATE


# ------------------------------
# Gene preparation helper
# ------------------------------
def _prepare_genes_chembl(state: "State") -> List[str]:
    """Resolve, deduplicate, and prepare gene symbols for ChEMBL agent."""
    genes = state.get("genes") or []
    if not genes and state.get("gene_symbol"):
        genes = [state["gene_symbol"]]

    if not genes:
        warnings.warn("No gene symbols provided.")
        state["used_tools"] = state.get("used_tools", []) + ["chembl"]
        return []

    # Deduplicate and preserve order
    unique_genes = []
    seen = set()
    for g in genes:
        if g not in seen:
            seen.add(g)
            unique_genes.append(g)

    # Ensure Gene objects exist (corrected instantiation)
    gene_entities = state.get("gene_entities") or {}
    for g in unique_genes:
        if not isinstance(gene_entities.get(g), Gene):
            gene_entities[g] = Gene(identifiers=GeneIdentifiers(symbol=g))
    state["gene_entities"] = gene_entities

    return unique_genes


# ------------------------------
# Main ChEMBL agent
# ------------------------------
def chembl_agent(state: "State") -> "State":
    """Resolve UniProt IDs, query ChEMBL for drug-target data,
    and return combined text + interactive HTML summary."""
    genes = _prepare_genes_chembl(state)
    if not genes:
        return state

    conn = get_connection()
    any_bioactivity = False
    chembl_data_all: Dict[str, Dict[str, List]] = {}

    for gene_symbol in genes:
        log(f"Processing gene: {gene_symbol}")

        # Resolve UniProt via UniProt API
        import requests
        uniprot_id = None
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
            if results:
                uniprot_id = results[0].get("primaryAccession")
        except Exception as e:
            warnings.warn(f"Error fetching UniProt for {gene_symbol}: {e}")

        if not uniprot_id:
            log(f"Skipping {gene_symbol}: could not resolve UniProt ID")
            continue

        # Fetch ChEMBL data (includes InChIKeys for 2D rendering)
        target_data = fetch_target_info(conn, uniprot_id)
        if target_data.get("bioactivity"):
            any_bioactivity = True

        # Store structured text summary
        summary = make_summary_text(gene_symbol, uniprot_id, target_data)
        gene_obj = state["gene_entities"][gene_symbol]
        gene_obj.update_text_summaries("ChEMBL drug-target evidence:\n" + summary)

        # Collect for frontend viewer
        chembl_data_all[gene_symbol] = target_data

    conn.close()

    # Add interpretation notes (only once, for the last gene)
    if any_bioactivity and genes:
        last_gene = genes[-1]
        if state["gene_entities"].get(last_gene):
            state["gene_entities"][last_gene].update_text_summaries(
                interpretation_notes(include_bioactivity=True)
            )

    # --- Build HTML frontend ---
    try:
        with open(HTML_TEMPLATE, "r", encoding="utf-8") as f:
            html_template = f.read()
        with open(CSS_TEMPLATE, "r", encoding="utf-8") as f:
            css_template = f.read()
        with open(JS_TEMPLATE, "r", encoding="utf-8") as f:
            js_template = f.read()
    except FileNotFoundError as e:
        raise RuntimeError(f"Missing frontend template: {e.filename}")

    html = inject_frontend_assets(
        html_template, css_template, js_template, chembl_data_all, ", ".join(genes)
    )

    state["used_tools"] = state.get("used_tools", []) + ["chembl"]
    if chembl_data_all:
        state["chembl_html"] = html

    return state


# ------------------------------
# CLI (testing mode)
# ------------------------------
if __name__ == "__main__":
    #genes = sys.argv[1:] if len(sys.argv) > 1 else ["TP53", "EGFR", "DRD2"]
    genes = sys.argv[1:] if len(sys.argv) > 1 else ["TMEM179", "C20orf85", "POTEM", "ZNF533", "WDR90"]
    
    base_name = genes[0] if len(genes) == 1 else f"{genes[0]}_plus{len(genes)-1}"

    state = State({
        "genes": genes,
        "gene_entities": {g: Gene(identifiers=GeneIdentifiers(symbol=g)) for g in genes},
    })

    result = chembl_agent(state)

    # --- Write outputs ---
    html_out = OUTPUT_DIR / f"{base_name}_chembl.html"
    html_out.write_text(result.get("chembl_html", "<p>No HTML produced.</p>"), encoding="utf-8")

    txt_out = OUTPUT_DIR / f"{base_name}_chembl.txt"
    lines: List[str] = []
    for g in genes:
        obj = result.get("gene_entities", {}).get(g)
        if obj and getattr(obj, "text_summaries_from_tools", None):
            lines.append(g)
            lines.extend(obj.text_summaries_from_tools)
            lines.append("")
    txt_out.write_text("\n".join(lines) if lines else "No summaries produced.", encoding="utf-8")

    print("[OK] Generated outputs:")
    print(f"  * HTML file: {html_out.resolve()}")
    print(f"  * TXT file:  {txt_out.resolve()}")
    print(f"[INFO] Genes processed: {', '.join(genes)}")


# ------------------------------
# Node Registration
# ------------------------------
NODES: tuple[Node, ...] = (
    Node(
        name="chembl",
        entry_point=chembl_agent,
        description=(
            "Query ChEMBL for drug-target information about one or more genes. "
            "For each gene: UniProt ID, summarizes FDA-approved drugs, "
            "clinical and preclinical trials, and assay bioactivity data. "
            "Generates text summaries and an interactive HTML viewer with 2D molecular renderings."
        ),
    ),
)
