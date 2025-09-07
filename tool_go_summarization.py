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

# def make_go_summarization_node(source, embedding_method='word2vec'):
#     def node(state: State):
#         genes = state.get(f"{source}_predictions", [])
#         return go_summarization_agent(state, genes, source, embedding_method)
#     return node


def go_summarization_agent(state: "State", embedding_method: str = 'word2vec') -> "State":
    """
    Summarize GO terms per gene and append a short line to each Gene's text summaries.
    Prefers the Gene.go_annotations field; falls back to UniProt if empty.
    """
    gene_objs = (state.get("gene_entities") or {}).copy()
    label = f"Summarized GO terms for: "

    for gene in gene_objs.keys():
        gobj = gene_objs.get(gene)
        if not gobj:
            continue

        # 1) Prefer existing GO annotations on the Gene objects
        terms = getattr(gobj, "go_annotations", []) or []
    
        # 2) If none, fetch from UniProt
        if not terms:
            try:
                entry = get_uniprot_entry_for_gene(gene)
                if entry:
                    terms = extract_GO_from_uniprot_entry(entry) or []
            except Exception:
                terms = []

        if not terms:
            continue

        # Dedup and select representative subset
        uniq = list(dict.fromkeys(terms))
        if embedding_method.lower() == "word2vec":
            try:
                idxs = fps_word2vec(uniq, 6, separate_sampling=True)
                picked = [uniq[i] for i in idxs]
            except Exception:
                picked = uniq[:6]
        else:  # fallback (tf-idf or anything else)
            raise NotImplementedError(f"Unknown embedding method: {embedding_method}")

        gobj.update_text_summaries(f"{gene} - " + ", ".join(picked)+". ")

    time.sleep(0.1)
    return
