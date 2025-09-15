from __future__ import annotations
import requests
import time
from typing import Any, Dict, List, Optional
from config import BioGRID_API_KEY
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