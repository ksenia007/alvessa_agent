"""
Author: Keerthana Nallamotu <kn6412@princeton.edu>
Contributors: 
Created: 2024-07-08
Updated: 2025-07-08


Description: 

Tool to run Sei model for variants of interest"""

import os
import sys
import pandas as pd
import time
import shutil
from state import State 
import warnings

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
 
    # Gracefully handle file reading errors
    try:
        seq_class_df_hg38 = pd.read_csv('local_dbs/sorted.hg38.tiling.bed.ipca_randomized_300.labels.merged.bed', names=['chr', 'start_pos', 'end_pos', 'seq_class'], sep = '\t')
        with open('local_dbs/seqclass.names', 'r') as f:
            seq_class_names = [line.strip() for line in f]
    except Exception as e:
        warnings.warn(f"Failed to load required files: {e}. Cannot run SEI predictions.")
        return {**state, "sei_predictions": preds}

    state_all_snps = {}

    for gene, gene_vars in variants.items():
        state_all_snps[gene] = {}
        
        try:
            for var, var_data in gene_vars.items():
                try:
                    all_snps = var_data['coordinates']

                    if 'assembly_filter' in var_data:
                        assembly = var_data['assembly_filter']
                    else:
                        continue

                    for snp in all_snps:
                        chrom = snp.get('chrom')
                        pos = snp.get('pos')
                        ref_base = snp.get('ref')
                        alt_base = snp.get('alt')

                        if all(x is not None for x in [chrom, pos, ref_base, alt_base, assembly]):
                            state_all_snps[gene][var] = [chrom, pos, ref_base, alt_base, assembly]
                        else:
                            warnings.warn(f"Missing coordinate data for {gene} variant {var}")
                except Exception as e:
                    warnings.warn(f"Failed to process variant {var} for gene {gene}: {e}")
                    continue
        except Exception as e:
            warnings.warn(f"Sei prediction unavailable for gene {gene}: {e}")
            continue

    for gene, variants in state_all_snps.items():
        preds[gene] = {}
        for var_id, (chrom, variant_pos, ref_base, alt_base, assembly) in variants.items():
            try:
                chr_str = f'chr{chrom}'

                if 'GRCh38' in assembly:
                    seq_class_df = seq_class_df_hg38
                else:
                    preds[gene][var_id] = None
                    continue

                match = seq_class_df[
                            (seq_class_df['chr'] == chr_str) &
                            (seq_class_df['start_pos'] <= variant_pos) &
                            (seq_class_df['end_pos'] >= variant_pos)
                        ]
                
                if not match.empty:
                    try:
                        seq_class_num = match['seq_class'].iloc[0]
                        seq_class_name = seq_class_names[int(seq_class_num)]
                    except (IndexError, ValueError) as e:
                        warnings.warn(f"Failed to get seq class name for {gene} variant {var_id}: {e}")
                        seq_class_name = None
                else:
                    if DEBUG:
                        print(f"No sequence class match found for: {chrom}, {variant_pos}, {ref_base}, {alt_base}")
                    seq_class_name = None
                    
                preds[gene][var_id] = seq_class_name

            except Exception as e:
                warnings.warn(f"Failed to process prediction for {gene} variant {var_id}: {e}")
                preds[gene][var_id] = None

    print(preds)
    time.sleep(0.3)  # courteous pause

    return {
        **state,
        "sei_predictions": preds}


