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
