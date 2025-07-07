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
from tool_humanbase import _fetch_predictions_HB, _filter_predictions_HB, _symbol_to_entrez
from state import State 

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
    """
    preds = state.get("bioGRID_predictions", {}).copy()

    for gene in state.get("genes", []):
        if gene in preds:
            continue

        preds[gene] = {}

        try:
            interactions = _fetch_predictions_BioGRID(gene)
        except Exception as exc:
            print(f"[BioGRID] {gene}: {exc}")
        else:

            for interactor in list(interactions)[:30]:
                entrez = _symbol_to_entrez(interactor)
                if not entrez:
                    functions_list = []
                    continue
        
                try:
                    tmp = _fetch_predictions_HB(entrez)
                except Exception as exc:
                    functions_list = []
                else:
                    functions_list = _filter_predictions_HB(tmp, threshold=0.95)

                if functions_list:
                    terms = [hit["term"] for hit in functions_list if "term" in hit]
                    if terms:
                        preds[gene][interactor] = {'functions': terms[:30]}
        
        time.sleep(0.3)  # courteous pause

    return {
        **state,
        "bioGRID_predictions": preds}
