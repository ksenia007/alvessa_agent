"""
Author: Ksenia 
Crated: 2025-11-24


Description: 

Fetch alliance genome information"""

from __future__ import annotations

import time
import requests

from src.state import State
from src.tools.base import Node

DEBUG = True

def get_alliancegenome_ids(gene_name: str) -> list[str]:
    """
    Fetch Alliance of Genome IDs for a given gene name.

    Parameters
    ----------
    gene_name : str
        The gene name to query.

    Returns
    -------
    list[str]
        List of Alliance of Genome IDs associated with the gene.
    """
    url = f"https://www.alliancegenome.org/api/search_autocomplete?q={gene_name}"
    r = requests.get(url, timeout=12)
    if r.status_code == 404:  # error
        return []
    r.raise_for_status()
    r = r.json()  # we get a dict {results: [...]}
    
    r = r.get("results", [])
    
    if not r:
        return []
    
    valid_ids = []
    for entry in r:
        if entry.get("symbol", "").upper() == gene_name.upper():
            valid_ids.append(entry.get("primaryKey"))
            
    return valid_ids

def get_information_per_id(alliance_id: str) -> dict:
    """
    Fetch detailed information for a given Alliance of Genome ID.

    Parameters
    ----------
    alliance_id : str
        The Alliance of Genome ID to query.

    Returns
    -------
    dict
        Detailed information about the gene associated with the ID.
    """
    url = f"https://www.alliancegenome.org/api/gene/{alliance_id}"
    r = requests.get(url, timeout=12)
    if r.status_code == 404:  # error
        return {}
    r.raise_for_status()
    return r.json()
    
    
def alliancegenome_predictions_agent(state: "State") -> "State":
    """
    LangGraph node that annotates each gene with Alliance of Genomes info.

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state 
    """
    gene_entities = state.get("gene_entities") or {}
    for gene_name, gene in gene_entities.items():
        if not gene_name:
            continue
        
        valid_ids = get_alliancegenome_ids(gene_name)
        if DEBUG: 
            print(f"[AllianceOfGenomes] Found IDs for {gene_name}: {valid_ids}")
        summary = 'Collected Alliance of Genome information across species. '
        for alliance_id in valid_ids:
            info = get_information_per_id(alliance_id)
            if not info:
                if DEBUG:
                    print(f"[AllianceOfGenomes] No info found for ID {alliance_id}")
                continue
            gene.add_alliancegenome_info(info)  # need to implement this method in Gene class, structured to future UI consumption
            
            # save summary text info. Create a string:
            summary += f"Matched to {alliance_id}, {info.get('species', {}).get('name', 'unknown species')}. Summary: {info.get('geneSynopsis', 'No summary available.')}. Auto-summary: {info.get('automatedGeneSynopsis', 'No auto-summary available.')} |"
        
        summary += " End of Alliance of Genome information."
        gene.update_text_summaries(summary.strip())
        gene.add_tool("AllianceOfGenomes")
        time.sleep(0.3)

    return 
NODES: tuple[Node, ...] = (
    Node(
        name="AllianceOfGenomes",
        entry_point=alliancegenome_predictions_agent,
        description="Fetches Alliance of Genome - consortium of 7 model organism databases - information for genes across all available model organisms. Pulls species specific summaries and gene descriptions.",
    ),
)
