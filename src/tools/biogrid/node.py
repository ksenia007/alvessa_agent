"""
Author: Keerthana Nallamotu <kn6412@princeton.edu>
Contributors: 
Created: 2024-06-27
Updated: 2025-09-12


Description: 

Tool to query BioGRID interactions with a list of gene symbols"""

from __future__ import annotations

import time

from src.state import State
from src.tools.base import Node
from .utils import _fetch_predictions_BioGRID

DEBUG = True

def bioGRID_predictions_agent(state: "State") -> "State":
    """
    LangGraph node that annotates each gene with bioGRID interactions.

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state with the `"bioGRID_predictions"` field filled.
        
    # TODO: Split the functionality, BioGRID should only fetch interactions, and HB called outside
    """
    gene_entities = state.get("gene_entities") or {}
    for gene_name, gene in gene_entities.items():
        if not gene_name:
            continue
        # if BioGrid was run before, skip genes we already have predictions for
        if gene.has_interactions_collected():
            if DEBUG:
                print(f"[BioGRID] Skipping {gene.symbol}, already has interactions collected.")
            continue

        try:
            _, human_set, nonhuman_set = _fetch_predictions_BioGRID(gene_name)
        except Exception as exc:
            print(f"[BioGRID] {gene}: {exc}")
        else:
            # populate Gene object
            for exp, all_genes_interacting in human_set.items():
                # key is the experimental system
                gene.add_many_interactions(exp, all_genes_interacting)
            for exp, all_genes_interacting in nonhuman_set.items():
                gene.add_many_nonhuman_interactions(exp, all_genes_interacting)
            gene.add_tool("BioGRID")

        time.sleep(0.3)  # courteous pause

    if DEBUG:
        print(f"[BioGRID] Predictions fetched")

    return 
NODES: tuple[Node, ...] = (
    Node(
        name="BioGRID",
        entry_point=bioGRID_predictions_agent,
        description="Fetches gene interactions from BioGRID for the input genes. Provides a curated context-specific list of protein-protein, genetic and chemical interactions.",
    ),
)
