"""
Author: Keerthana Nallamotu <kn6412@princeton.edu>
Contributors: 
Created: 2024-08-05
Updated: 2025-08-05


Description: 

Tool to query Reactome Pathways with a list of gene symbols"""

from __future__ import annotations
from src.state import State 
import pandas as pd

from src.tools.base import Node

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
    
    gene_list = state.get("gene_entities", [])
    print(gene_list)
    
    
    try:
        df = pd.read_csv('local_dbs/NCBI2Reactome_All_Levels.txt', 
                         names = ['geneID', 'pathwayID', 'url', 'pathway_name', 'evidence_code', 'species'], sep = '\t')
    except:
        print("[Reactome] Could not load local Reactome database file.")
        return 


    for gene in gene_list.keys():
        # if reactome already ran
        if gene_list[gene].has_tool("reactome_pathways"):
            continue

        entrez = gene_list[gene].entrez_id
        if not entrez:
            if DEBUG:
                print(f"[Reactome] Could not find Entrez ID for gene symbol: {gene}")
            gene_list[gene].add_tool("reactome_pathways")
            continue
        
        match = df[df['geneID']==int(entrez)]
        gene_list[gene].add_many_pathways(list(set(match['pathway_name'].values)))
        gene_list[gene].add_tool("reactome_pathways")
    
    return 


NODES: tuple[Node, ...] = (
    Node(
        name="reactome",
        entry_point=reactome_pathways_agent,
        description="Fetches Reactome pathways associated with the input genes, providing curated, expert-reviewed pathway information that describe molecular interactions, signaling cascades, and biological processes in which the genes participate.",
    ),
)
