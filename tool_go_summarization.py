"""
Author: Keerthana Nallamotu <kn6412@princeton.edu>
Contributors: 
Created: 2024-07-11
Updated: 2025-07-14


Description: 

Tool to summarize the GO terms associated with a list of genes"""

from __future__ import annotations
import requests
import time
from typing import Any, Dict, List, Optional
from config import BioGRID_API_KEY
from state import State 
from tools.word2vec import fps_word2vec
from tool_uniprot import get_uniprot_entry_for_gene, extract_GO_from_uniprot_entry

DEBUG=True

def make_go_summarization_node(source, embedding_method):
    def node(state: State):
        genes = state.get(f"{source}_predictions", [])
        return go_summarization_agent(state, genes, source, embedding_method)
    return node


def go_summarization_agent(state: "State", genes: List, source: str, embedding_method: str) -> "State":
    """
    LangGraph node that summarizes the GO terms associated with a list of genes.

    Parameters
    ----------
    state
        Current graph state.

    genes
        List of genes to summarize

    source
        Source of the gene list of interest

    embedding_method
        Use 'tf_idf' or 'word2vec' to specify embedding method

    Returns
    -------
    State
        Updated state with the `"{source}_predictions"` field filled.
        
    """
    preds = state.get(f"{source}_summarized_go", {}).copy()
    uniprot_entries: Dict[str, Dict] = {}

    all_go_terms = []
    for gene in list(genes):
        if gene:
            interactor_entry = get_uniprot_entry_for_gene(gene)

            if interactor_entry:
                traits = extract_GO_from_uniprot_entry(interactor_entry)
                all_go_terms.extend(traits)
        
    if all_go_terms:
        fps_indices = fps_word2vec(list(set(all_go_terms)), 6, separate_sampling=True)
        selected_go_terms = [all_go_terms[i] for i in fps_indices]
        
        preds[gene] = selected_go_terms

    time.sleep(0.3)  # courteous pause

    return {
        **state,
        f"{source}_summarized_go": preds}
