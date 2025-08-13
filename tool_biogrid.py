"""
Author: Keerthana Nallamotu <kn6412@princeton.edu>
Contributors: 
Created: 2024-06-27
Updated: 2025-06-30


Description: 

Tool to query BioGRID interactions with a list of gene symbols"""

from __future__ import annotations
import requests
import time
from typing import Any, Dict, List, Optional
from config import BioGRID_API_KEY
from state import State 
from tools.word2vec import fps_word2vec
from tool_go_summarization import make_go_summarization_node

DEBUG=True


def _fetch_predictions_BioGRID(gene_of_interest: str) -> List[Dict[str, float]]:
    """Download BioGRID functional predictions for a given gene symbol."""
    if DEBUG:
        print(f"[BioGRID] Fetching predictions for gene symbol: {gene_of_interest}")

    accesskey = BioGRID_API_KEY
    
    # API endpoint
    url = "https://webservice.thebiogrid.org/interactions/"
    params = {
        "accesskey": accesskey,
        "geneList": gene_of_interest,
        "format": "json",
        "includeInteractors": "true"
    }
    
    response = requests.get(url, params=params)
    data = response.json()
    
    output_list = set()
    
    for interaction_id, interaction_info in data.items():
        
        if gene_of_interest == interaction_info["OFFICIAL_SYMBOL_A"]:
            interactor = interaction_info["OFFICIAL_SYMBOL_B"]
        else:
            interactor = interaction_info["OFFICIAL_SYMBOL_A"]
    
        output_list.add(interactor)

    return output_list


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
    preds = state.get("bioGRID_predictions", {}).copy()

    for gene in state.get("genes", []):
        if gene in preds:
            continue

        try:
            interactions = _fetch_predictions_BioGRID(gene)
        except Exception as exc:
            print(f"[BioGRID] {gene}: {exc}")
        else:
            preds[gene] = list(interactions)[:100]

        time.sleep(0.3)  # courteous pause

    if DEBUG:
        print(f"[BioGRID] Predictions fetched for {len(preds)} genes.")
        print("[BioGRID] Example predictions:", list(preds.items())[:2])

    return {
        "biogrid_predictions": preds}
