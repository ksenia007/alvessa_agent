# src/tools/chembl/node.py
# Author: Dmitri Kosenkov
# Created: 2025-09-25
# Updated: 2025-09-26
#
# Main entry point for ChEMBL drug-target summarization tool.

import sys
import warnings
from pathlib import Path
from typing import List

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
    fetch_target_info,
    make_summary_text,
    interpretation_notes,
    log,
)
from src.tools.chembl import OUTPUT_DIR


# ------------------------------
# Main agent
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

    # Deduplicate while preserving input order
    seen = set()
    unique_genes = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            unique_genes.append(g)
    genes = unique_genes

    # Ensure gene_entities dict exists and contains Gene objects
    gene_entities = state.get("gene_entities") or {}
    for g in genes:
        if not isinstance(gene_entities.get(g), Gene):
            gene_entities[g] = Gene(symbol=g)
    state["gene_entities"] = gene_entities

    return genes

def chembl_agent(state: "State") -> "State":
    """Resolve UniProt IDs, query ChEMBL for drug-target data, and return text summaries."""
    # --- Gene handling via helper ---
    genes = _prepare_genes_chembl(state)
    if not genes:
        return state

    conn = get_connection()
    any_bioactivity = False

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

        target_data = fetch_target_info(conn, uniprot_id)
        if target_data.get("bioactivity"):
            any_bioactivity = True

        summary = make_summary_text(gene_symbol, uniprot_id, target_data)
        state["gene_entities"][gene_symbol].update_text_summaries(
            "ChEMBL drug-target evidence:\n" + summary
        )

    conn.close()

    # Add interpretation notes if bioactivity data was found
    if any_bioactivity and genes:
        last_gene = genes[-1]
        if state["gene_entities"].get(last_gene):
            state["gene_entities"][last_gene].update_text_summaries(
                interpretation_notes(include_bioactivity=True)
            )

    state["used_tools"] = state.get("used_tools", []) + ["chembl"]
    return state

# ------------------------------
# CLI (testing)
# ------------------------------
if __name__ == "__main__":
    genes = sys.argv[1:] if len(sys.argv) > 1 else ["TP53", "EGFR", "DRD2"]
    base_name = genes[0] if len(genes) == 1 else f"{genes[0]}_plus{len(genes)-1}"

    state = State({
        "genes": genes,
        "gene_entities": {g: Gene(GeneIdentifiers(symbol=g)) for g in genes},
    })
    result = chembl_agent(state)

    txt_out = OUTPUT_DIR / f"{base_name}_chembl.txt"
    lines: List[str] = []
    for g in genes:
        obj = result.get("gene_entities", {}).get(g)
        if obj and getattr(obj, "text_summaries_from_tools", None):
            lines.append(g)
            lines.extend(obj.text_summaries_from_tools)
            lines.append("")
    with txt_out.open("w", encoding="utf-8") as f:
        f.write("\n".join(lines) if lines else "No summaries produced.")

    print("[OK] Generated outputs:")
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
            "For each gene -> UniProt ID, summarizes FDA-approved drugs, "
            "clinical and preclinical trials, and assay bioactivity data. "
            "Produces text summaries only, with interpretation notes for assays."
        ),
    ),
)
