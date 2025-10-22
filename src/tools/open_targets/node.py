"""
Author: Keerthana Nallamotu <kn6412@princeton.edu>
Contributors: 
Created: 2024-10-4
Updated: 2025-10-6


Description: 

Tool to query OpenTarget annotations with a list of gene symbols"""

from __future__ import annotations

import time

from src.state import State
from src.tools.base import Node
from .utils import read_all_parquet_in_folder
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

DEBUG = True

REPO_ROOT = Path(__file__).resolve().parents[3]
LOCAL_DBS = REPO_ROOT / "local_dbs"
OPEN_TARGETS = LOCAL_DBS / "open_targets"
TARGET_DISEASE_DATA = OPEN_TARGETS / "final_association_overall_direct"
EXPRESSION_DATA = OPEN_TARGETS / "final_expression"
ESSENTIALITY_DATA = OPEN_TARGETS / "final_target_essentiality"
CONSTRAINT_DATA = OPEN_TARGETS / "target"

CONSTRAINT_TYPES = {
    'syn': "Synonymous (silent) variants - don't change amino acid sequence of a protein",
    'mis': "Missense variants - lead to amino acid substitution in the protein",
    'lof': "Loss-of-function variants - include nonsense, frameshift, and canonical splice site variants"
}

def opentargets_agent(state: "State") -> "State":
    """
    LangGraph node that annotates each gene with Open Target data.

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state with the open target fields filled.
        
    """
    gene_entities = state.get("gene_entities") or {}

    target_disease_df = read_all_parquet_in_folder(TARGET_DISEASE_DATA)
    expression_df = read_all_parquet_in_folder(EXPRESSION_DATA)
    essentiality_df = read_all_parquet_in_folder(ESSENTIALITY_DATA)
    constraint_df = read_all_parquet_in_folder(CONSTRAINT_DATA)

    for gene_name, gene in gene_entities.items():
        if not gene_name:
            continue

        summary_lines: List[str] = []

        try:
            # Disease annotations
            associated_disease_list = list(set(target_disease_df[target_disease_df['target_symbol'] == gene.symbol]['disease_name'].tolist()))
            gene.add_many_direct_disease_associations(associated_disease_list)

            summary_lines.append(f"All direct disease associations for {gene.symbol} (from the Open Targets database): " 
                    + ', '.join(associated_disease_list) + ".")

            print(len(associated_disease_list))
        
        except Exception as e:
            print(f"Failed to pull disease annotations: {gene}, {e}")

        try:
            # Tissue-specific expression
            tissue_zscore = {}
            for tissues_list in expression_df[expression_df['target_symbol'] == gene.symbol]['tissues']:
                if not tissues_list:
                    continue
                for tissue_dict in tissues_list:
                    label = tissue_dict.get('label', None)
                    zscore = tissue_dict.get('rna', {}).get('zscore', None)
                    if label is not None and zscore is not None:
                        tissue_zscore[label] = zscore
                    
        gene.add_many_tissues_expression(tissue_zscore)

            summary_lines.append(f"Tissue-specific expression z-scores for {gene.symbol} (from the Open Targets database). A gene is considered to be tissue specific if the z-score for that tissue is greater than 0.674 (or the 75th percentile of a perfect normal distribution): " 
                        + str(tissue_zscore) + ".")
            
        except Exception as e:
            print(f"Failed to pull tissue expression: {gene}, {e}")
        
        try:
            # DepMap Essentiality
            gene_is_essential = essentiality_df[essentiality_df['target_symbol'] == gene.symbol]['geneEssentiality'].values[0][0]['isEssential']
            gene.add_essentiality(gene_is_essential)

            if gene_is_essential:
                summary_lines.append(f"{gene.symbol} is a core essential gene, meaning that it is unlikely to tolerate inhibition and is susceptible to causing adverse events if modulated. The gene is crucial for basic cellular function in many tissue types.")
            else:
                summary_lines.append(f"{gene.symbol} is not a core essential gene, meaning that inhibition or knockout of this gene does not consistently result in cell death across the majority of cell lines tested. The gene's function is not universally vital for cell survival in diverse tissues, as determined by large-scale cell fitness and depletion assays.")
        
        except Exception as e:
            print(f"Failed to pull essentiality: {gene}, {e}")

        try:       
            # Genetic Constraint
            gene_all_constraint = constraint_df[constraint_df['approvedSymbol'] == gene.symbol]['constraint'].values[0]
            if not gene_all_constraint:
                if DEBUG:
                    print(f"[Open Targets] Empty constraint data for {gene.symbol}")
                    print(constraint_df[constraint_df['approvedSymbol'] == gene.symbol])
            else:
                for row in gene_all_constraint:
                    constraint_type_name = row['constraintType']
                    score = row['score']

                    gene.add_constraint(score, constraint_type_name)
                    summary_lines.append(f"Genetic constraint score for {constraint_type_name} variants ({CONSTRAINT_TYPES[constraint_type_name]}) in {gene.symbol} from the Open Targets database. A constraint score from -1 to 1 is given to genes depending on their LOEUF (loss-of-function observed/expected upper bound fraction) metric rank, with -1 being the least tolerant to LoF variation and 1 being the most tolerant.: " + str(score))

        gene.update_text_summaries(
            " ".join(summary_lines)
        )
    
        gene.add_tool("OpenTargets")

        except Exception as e:
            print(f"Failed to pull constraint: {gene}, {e}")

        time.sleep(0.3)  # courteous pause

    if DEBUG:
        print(f"[Open Targets] Predictions fetched")

    return 
NODES: tuple[Node, ...] = (
    Node(
        name="OpenTargets",
        entry_point=opentargets_agent,
        description="Tool that fetches detailed information on the target (gene)-disease associations, tissue-specific expression of the gene of interest (experimental), essentiality of the gene, and genetic constraint for synonymous, missense, and loss-of-function variants from Open Targets. This is a reliable and verified source for gene-disease associations and target annotations."
    ),
)