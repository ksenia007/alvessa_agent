from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Set, Tuple

import requests

from src.config import BioGRID_API_KEY, DEBUG

# ALLOWED_NONHUMAN_EXPERIMENTS = {
#     "Dosage Growth Defect",
#     "Dosage Lethality",
#     "Dosage Rescue",
#     "Negative Genetic",
#     "Phenotypic Enhancement",
#     "Phenotypic Suppression",
#     "Positive Genetic",
#     "Synthetic Growth Defect",
#     "Synthetic Haploinsufficiency",
#     "Synthetic Lethality",
#     "Synthetic Rescue",
# }


def _fetch_predictions_BioGRID(
    gene_symbol: str,
) -> Tuple[Set[str], Dict[str, List[str]], Dict[str, List[str]]]:
    """Return BioGRID interaction partners partitioned by organism.

    Parameters
    ----------
    gene_symbol:
        HGNC gene symbol to query.

    Returns
    -------
    tuple
        A triple of (all unique partners, human-only partners by experiment,
        non-human partners by experiment).
    """

    if DEBUG:
        print(f"[BioGRID] Fetching predictions for gene symbol: {gene_symbol}")

    params = {
        "accesskey": BioGRID_API_KEY,
        "geneList": gene_symbol,
        "format": "json",
        "includeInteractors": "true",
    }

    response = requests.get("https://webservice.thebiogrid.org/interactions/", params=params)
    response.raise_for_status()
    data = response.json() or {}

    unique_partners: Set[str] = set()
    human_partners: Dict[str, List[str]] = defaultdict(list)
    nonhuman_partners: Dict[str, List[str]] = defaultdict(list)

    for interaction in data.values():
        symbol_a = interaction.get("OFFICIAL_SYMBOL_A")
        symbol_b = interaction.get("OFFICIAL_SYMBOL_B")
        organism_a = interaction.get("ORGANISM_A")
        organism_b = interaction.get("ORGANISM_B")

        if gene_symbol == symbol_a:
            interactor = symbol_b
            interactor_organism = organism_b
        else:
            interactor = symbol_a
            interactor_organism = organism_a

        if not interactor:
            continue

        experimental_system = interaction.get("EXPERIMENTAL_SYSTEM")
        if interactor_organism == 9606:
            human_partners[experimental_system].append(interactor)
        else:
            nonhuman_partners[experimental_system].append(interactor)

        unique_partners.add(interactor)

    return unique_partners, human_partners, nonhuman_partners
