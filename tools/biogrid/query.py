from typing import List, Dict
from tools.biogrid.consts import BASE_URL, API_KEY
import requests
import pandas as pd
import io

INTERACTIONS_ENDPOINT = "interactions"

def get_interactions(gene_list: List[str]) -> List[Dict]:
    """
    Get interactions for a list of genes.
    """
    request_url = f"{BASE_URL}{INTERACTIONS_ENDPOINT}"
    gene_list = list(set(gene_list))
    gene_list_str = "|".join(gene_list) if len(gene_list) > 1 else gene_list[0]
    evidence_list = ["POSITIVE GENETIC", "PHENOTYPIC ENHANCEMENT"]

    params = {
        "accesskey": API_KEY,
        "format": "tab2",  # Return results in TAB2 format
        "geneList": gene_list_str,  # Must be | separated
        "searchNames": "true",  # Search against official names
        "includeInteractors": "true",  # Set to true to get any interaction involving EITHER gene, set to false to get interactions between genes
        "taxId": 559292,  # Limit to Saccharomyces cerevisiae
        "evidenceList": "|".join(evidence_list),  # Exclude these two evidence types
        "includeEvidence": "false",  # If false "evidenceList" is evidence to exclude, if true "evidenceList" is evidence to show
        "includeHeader": "true",
    }
    r = requests.get(request_url, params=params)
    print(r.text)
    interactions = pd.read_csv(io.StringIO(r.text), sep="\t")
    return interactions


if __name__ == "__main__":
    interactions = get_interactions(["STE11"])
    print(interactions.columns)
    print(interactions)