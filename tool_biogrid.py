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
from collections import defaultdict

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

    human_set = defaultdict(list)
    nonhuman_set = defaultdict(list)
    allowed_nonhuman_experiment_systems = ['Dosage Growth Defect', 'Dosage Lethality', 'Dosage Rescue', 'Negative Genetic', 'Phenotypic Enhancement', 'Phenotypic Suppression', 'Positive Genetic', 'Synthetic Growth Defect', 'Synthetic Haploinsufficiency', 'Synthetic Lethality', 'Synthetic Rescue', 'Synthetic Growth Defect']
    
    for interaction_id, interaction_info in data.items():
        # print(interaction_info)
        
        if gene_of_interest == interaction_info["OFFICIAL_SYMBOL_A"]:
            interactor = interaction_info["OFFICIAL_SYMBOL_B"]
            interactor_organism = interaction_info["ORGANISM_B"]

        else:
            interactor = interaction_info["OFFICIAL_SYMBOL_A"]
            interactor_organism = interaction_info["ORGANISM_A"]
        
        experimental_system = interaction_info['EXPERIMENTAL_SYSTEM']
        if interactor_organism == 9606:
            human_set[experimental_system].append(interactor)
        elif experimental_system in allowed_nonhuman_experiment_systems:
            # print(interactor_organism)
            nonhuman_set[experimental_system].append(interactor)
        output_list.add(interactor)

    return output_list, human_set, nonhuman_set


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
    human_interactions = state.get("biogrid_interaction_groups", {}).copy()
    nonhuman_select_interactions = state.get("biogrid_interactions_select_nonhuman", {}).copy()

    # human_interactions = defaultdict(list, state.get("biogrid_interaction_groups", {}))
    # nonhuman_select_interactions = defaultdict(list, state.get("biogrid_interactions_select_nonhuman", {}))

    for gene in state.get("genes", []):
        if gene in preds:
            continue

        try:
            interactions, human_set, nonhuman_set = _fetch_predictions_BioGRID(gene)

            # print(human_set)
            # print(nonhuman_set)
        except Exception as exc:
            print(f"[BioGRID] {gene}: {exc}")
        else:
            human_interactions[gene] = defaultdict(list)
            nonhuman_select_interactions[gene] = defaultdict(list)

            for key, val in human_set.items():
                human_interactions[gene][key].extend(list(set(val)))
            for key, val in nonhuman_set.items():
                nonhuman_select_interactions[gene][key].extend(list(set(val)))

            preds[gene] = list(set(interactions))

        time.sleep(0.3)  # courteous pause

    if DEBUG:
        print(f"[BioGRID] Predictions fetched for {len(preds)} genes.")
        print("[BioGRID] Example predictions:", list(preds.items())[:2])

    return {
        "biogrid_predictions": preds,
        "biogrid_interaction_groups": human_interactions,
        "biogrid_interactions_select_nonhuman": nonhuman_select_interactions}
