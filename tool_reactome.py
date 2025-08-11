"""
Author: Keerthana Nallamotu <kn6412@princeton.edu>
Contributors: 
Created: 2024-08-05
Updated: 2025-08-05


Description: 

Tool to query Reactome Pathways with a list of gene symbols"""

from __future__ import annotations
from state import State 
import pandas as pd
from tool_humanbase import _symbol_to_entrez

DEBUG=True


def reactome_pathways_agent(state: "State") -> "State":
    """
    LangGraph node that annotates each gene with reactome pathways.

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state with the `"reactome_pathways"` field filled.
        
    """
    preds = state.get("reactome_pathways", {}).copy()

    for gene in state.get("genes", []):
        if gene in preds:
            continue

        entrez = _symbol_to_entrez(gene)
        if not entrez:
            preds[gene] = []
            continue
        
        df = pd.read_csv('local_dbs/NCBI2Reactome_All_Levels.txt', names = ['geneID', 'pathwayID', 'url', 'pathway_name', 'evidence_code', 'species'], sep = '\t')
        match = df[df['geneID']==int(entrez)]

        preds[gene] = list(set(match['pathway_name'].values))

    return {
        **state,
        "reactome_pathways": preds}
