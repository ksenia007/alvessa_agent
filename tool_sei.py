"""
Author: Keerthana Nallamotu <kn6412@princeton.edu>
Contributors: 
Created: 2024-07-08
Updated: 2025-07-08


Description: 

Tool to run Sei model for variants of interest"""

import os
import sys

sys.path.append('../../sei')
import sei_argo as sei
import pandas as pd
import time
import shutil
from state import State 

DEBUG=True

def sei_predictions_agent(state: "State") -> "State":
    """
    LangGraph node that runs sei predictions for 

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state with the `"sei_predictions"` field filled.
        
    """
    preds = state.get("sei_predictions", {}).copy()
    variants = state.get("dbsnp_variants", {}).copy()

    seq_class_df = pd.read_csv('local_dbs/sorted.hg38.tiling.bed.ipca_randomized_300.labels.merged.bed', names=['chr', 'start_pos', 'end_pos', 'seq_class'], sep = '\t')

    with open('local_dbs/seqclass.names', 'r') as f:
        seq_class_names = [line.strip() for line in f]

    state_all_snps = {}

    for gene, gene_vars in variants.items():
        state_all_snps[gene] = {}
        for var, var_data in gene_vars.items():
            all_snps = var_data['coordinates']
            for snp in all_snps:
                chrom = snp['chrom']
                pos = snp['pos']
                ref_base = snp['ref']
                alt_base = snp['alt']
                state_all_snps[gene][var] = [chrom, pos, ref_base, alt_base]

    for gene, variants in state_all_snps.items():
        preds[gene] = {}
        for var_id, (chrom, variant_pos, ref_base, alt_base) in variants.items():
            chr_str = f'chr{chrom}'

            match = seq_class_df[
                        (seq_class_df['chr'] == chr_str) &
                        (seq_class_df['start_pos'] <= variant_pos) &
                        (seq_class_df['end_pos'] >= variant_pos)
                    ]
            
            if not match.empty:
                seq_class_num = match['seq_class'].iloc[0]
                seq_class_name = seq_class_names[int(seq_class_num)] 
            else:
                print(chrom, variant_pos, ref_base, alt_base)
                seq_class_name = None
                
            preds[gene][var_id] = seq_class_name

    print(preds)
    time.sleep(0.3)  # courteous pause

    return {
        **state,
        "sei_predictions": preds}


