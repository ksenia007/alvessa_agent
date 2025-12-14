"""
Author: Keerthana Nallamotu <kn6412@princeton.edu>
Contributors: 
Created: 2024-9-27
Updated: 2025-10-10


Description: 

Tool to query OMIM annotations with a list of gene symbols"""

from __future__ import annotations

import time

from src.state import State
from src.tools.base import Node
from pathlib import Path
from typing import Dict, Iterable, List, Sequence
import json
import pandas as pd

DEBUG = True

REPO_ROOT = Path(__file__).resolve().parents[3]
LOCAL_DBS = REPO_ROOT / "local_dbs"
OMIM = LOCAL_DBS / "omim_genemap2.txt"

def omim_agent(state: "State") -> "State":
    """
    LangGraph node that annotates each gene with OMIM data.

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state with the omim fields filled.
        
    """
    gene_entities = state.get("gene_entities") or {}

    omim_df = pd.read_csv(OMIM, sep = '\t', comment='#', header=None, names = ['Chromosome','Genomic Position Start','Genomic Position End','Cyto Location','Computed Cyto Location','MIM Number','Gene/Locus And Other Related Symbols','Gene Name','Approved Gene Symbol','Entrez Gene ID','Ensembl Gene ID','Comments','Phenotypes','Mouse Gene Symbol/ID'])

    # Fixing weird formatting they have for phenotypes
    pattern = r'(, \d+ \(\d+\))|(\s\(\d+\))'
    omim_df['extracted'] = omim_df.loc[:,'Phenotypes'].str.replace(pattern, '', regex=True).str.replace('[\{\}\[\]\?]', '', regex=True)

    for gene_name, gene in gene_entities.items():
        if not gene_name:
            continue

        summary_lines: List[str] = []

        try:
            relevant_phenotypes = omim_df[omim_df['Approved Gene Symbol'] == gene.symbol]['extracted'].values[0].split(';')

            if len(relevant_phenotypes) > 0:
                summary_lines.append(f"*OMIM: Phenotypes annotated to {gene.symbol}: " + "; ".join(relevant_phenotypes) + ".")
                gene.add_many_omim_annotated_terms(relevant_phenotypes)

                if summary_lines:
                    gene.update_text_summaries(" ".join(summary_lines))
                
                    gene.add_tool("OMIM")

        except:
            print(gene)


        time.sleep(0.3)  # courteous pause

    if DEBUG:
        print(f"[OMIM] Predictions fetched")

    return 
NODES: tuple[Node, ...] = (
    Node(
        name="OMIM",
        entry_point=omim_agent,
        description=" Fetches curated OMIM entries describing Mendelian disease phenotypes caused by pathogenic variants in the input genes, including clinical features, inheritance patterns, and molecular genetic evidence.",
    ),
)
