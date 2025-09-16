"""
Author: Keerthana Nallamotu <kn6412@princeton.edu>
Contributors: 
Created: 2024-06-27
Updated: 2025-09-12


Description: 

Tool to query BioGRID interactions with a list of gene symbols"""

from __future__ import annotations
import requests
import time
from typing import Any, Dict, List, Optional
from src.config import BioGRID_API_KEY
from src.state import State 
from src.tools.base import Node
from src.tools.shared.word2vec import fps_word2vec
from collections import defaultdict
from tools.biogrid.utils import _fetch_predictions_BioGRID

DEBUG=True

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
    gene_obj = state.get("gene_entities", [])
    for gene_name in gene_obj.keys():
        if not gene_name:
            continue
        gene = gene_obj[gene_name]
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

        time.sleep(0.3)  # courteous pause

    if DEBUG:
        print(f"[BioGRID] Predictions fetched")

    return 


NODES: tuple[Node, ...] = (
    Node(
        name="BioGRID",
        entry_point=bioGRID_predictions_agent,
        description="Fetch BioGRID interactions for each gene and record human/non-human evidence.",
    ),
)
