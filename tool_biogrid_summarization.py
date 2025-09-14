"""
Author: Keerthana Nallamotu <kn6412@princeton.edu>
Contributors: 
Created: 2024-08-29
Updated: 2025-09-12


Description: 

Tool to run GO enrichment on biogrid interactors"""

from __future__ import annotations
from tools.enrichment.enrichment_test import run_go_enrichment
from state import State 
import pandas as pd
import requests
from datetime import datetime
import re
from tools.biogrid.utils import _fetch_predictions_BioGRID

DEBUG=True

def Summarize_bioGRID_GO_agent(state: "State") -> "State":
    """
    LangGraph node that runs GO enrichment on all human BIOGRID interactors of a gene.

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state with the `"interactors_summary_*"` fields filled.
        
    """

    gene_objs = state.get('gene_entities', {})
    go_categories = ['BP', 'MF', 'CC']

    full_category_name_mapping = {
        'BP': 'Biological Process (larger sequences of molecular activities)',
        'MF': 'Molecular Function (biochemical activities of gene products)', 
        'CC': 'Cellular Component (locations within the cell where gene products are active)'
    }

    old_go_folder = 'local_dbs/old_go_data/'
    pan_go_folder = 'local_dbs/pan_go_data/'
    old_go_obo_file = 'local_dbs/old_go_data/go.obo'
    pan_go_obo_file = 'local_dbs/pan_go_data/go.obo'

    for gene_name in gene_objs.keys():

        if not gene_name:
            continue
        gene = gene_objs[gene_name]

        summary = ""
        all_human_interactors_set = set()

        try:
            _, human_set, nonhuman_set = _fetch_predictions_BioGRID(gene_name)
        except Exception as exc:
            print(f"[BioGRID] {gene}: {exc}")
        else:
            for exp, all_genes_interacting in human_set.items():
                all_human_interactors_set.update(all_genes_interacting)

        for category in go_categories:
            old_go_enrichment_results = run_go_enrichment(list(all_human_interactors_set), f'{old_go_folder}/goa_human_{category}.gmt', old_go_obo_file)
            pan_go_enrichment_results = run_go_enrichment(list(all_human_interactors_set), f'{pan_go_folder}/functionome_release_{category}.gmt', pan_go_obo_file)

            old_go_enriched_go_terms = list(filter(None, old_go_enrichment_results.GO_Name.values))
            pan_go_enriched_go_terms = list(filter(None, pan_go_enrichment_results.GO_Name.values))
            
            gene.add_interactors_enriched_GO_terms('old_go', category, old_go_enriched_go_terms)
            gene.add_interactors_enriched_GO_terms('pan_go', category, pan_go_enriched_go_terms)

            # create a text summary for enriched GO Terms
            summary += f"All enriched PAN-GO terms (functional annotations curated using phylogenetic trees that are more informative and less redundant than standard GO terms) in the category of {full_category_name_mapping[category]} for all human interacting genes (fetched from BioGRID) of {gene}, include: "
            summary += f"{', '.join(pan_go_enriched_go_terms)}; "

        summary += f" End of record for {gene} |"
        gene.update_text_summaries(summary)

        gene.add_tool("Summarize_bioGRID_GO_agent")



    return 
