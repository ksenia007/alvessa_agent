"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-08-11
Updated: 


Description: 

Agent wrapper to annotate genes and variants with GENCODE features.
"""

from __future__ import annotations
import requests
import re
from typing import Dict, List, Optional
from state import State
from tools.gencode.query import summarize_gene_structure, annotate_variant

DEBUG = True

def gencode_gene_node(state: "State") -> "State":
    """
    LangGraph node that annotates genes with GENCODE features.

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state with the `gene_structure` field filled.
    """
    gene_objects = state.get("gene_entities", [])
    gene_list = list(gene_objects.keys())
    
    if DEBUG:
        print("[gencode_gene_node] Annotating genes with GENCODE features:", gene_list)

    gene_structure = {}
    
    for gene in gene_list:
        if gene:
            gene_obj = gene_objects[gene]
            try:
                structure = summarize_gene_structure(gene)
                gene_obj.set_gene_type(structure.get("gene_type", "unknown"))
                gene_obj.set_gene_ids(ensemble_id = structure.get("gene_id", ""))
                gene_obj.set_chrom_location(chrom = structure.get("chromosome", ""),
                                            gene_span=structure.get("gene_span_bp", (0,0)))
                print('structure.get("transcript_ids", [])', structure.get("transcript_ids", []))
                for i in structure.get("transcript_ids", []):
                    gene_obj.add_transcript(i, n_exons=structure.get("exons_per_transcript", {}).get(i, 0))
                    
                gene_obj.add_tool("gencode_gene_node")
            except Exception as e:
                print(f"[gencode_gene_node] Error annotating gene {gene}: {e}")
                    
    gene_structure['assembly'] = "GRCh38"  # Gencode uses GRCh38 as the default assembly
    if DEBUG:
        print("[gencode_gene_node] Gene structure summary:", gene_structure)
    
    return {'gene_level_gencode': gene_structure}
