"""
Author: Keerthana Nallamotu <kn6412@princeton.edu>
Contributors: 
Created: 2024-08-05
Updated: 2025-08-10


Description: 

Tool to extract AlphaMissense prediction for variants of interest"""

import pandas as pd
import time
from state import State 
import warnings
from datetime import datetime

DEBUG=True

def alphamissense_predictions_agent(state: "State") -> "State":
    """
    LangGraph node that runs alphamissense predictions for 

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state with the `"alphamissense_predictions"` field filled.
        
    """
    preds = state.get("alphamissense_predictions", {}).copy()
    variants = state.get("dbsnp_variants", {}).copy()

    # Gracefully handle file reading errors
    try:
        pathogenicity_class_df_hg38 = pd.read_parquet('local_dbs/AlphaMissense_hg38.parquet')
    except Exception as e:
        warnings.warn(f"Failed to load required files: {e}. Cannot run AlphaMissense predictions.")
        return {**state, "alphamissense_predictions": preds}
    
    state_all_snps = {}

    for gene, gene_vars in variants.items():
        state_all_snps[gene] = {}

        try:
            for var, var_data in gene_vars.items():
                try:
                    all_snps = var_data['coordinates']

                    state_all_snps[gene][var] = []

                    for snp in all_snps:
                        chrom = snp['chrom']
                        pos = snp['pos']
                        ref_base = snp['ref']
                        alt_base = snp['alt']
                        assembly = snp['assembly']

                        if 'GRCh38' in assembly:
                            if all(x is not None for x in [chrom, pos, ref_base, alt_base]):
                                state_all_snps[gene][var].append([chrom, pos, ref_base, alt_base])
                            else:
                                warnings.warn(f"Missing coordinate data for {gene} variant {var} (SNP {ref_base} -> {alt_base})")
                        else:
                            continue
                        
                except Exception as e:
                    warnings.warn(f"Failed to process variant {var} for gene {gene}: {e}")
                    continue

        except Exception as e:
            warnings.warn(f"AlphaMissense prediction unavailable for gene {gene}: {e}")
            continue

    for gene, variants in state_all_snps.items():
        preds[gene] = {}
        for var_id, snp_list in variants.items():
            preds[gene][var_id] = {}
            for (chrom, variant_pos, ref_base, alt_base) in snp_list:

                snp_key = f"SNP:{ref_base}->{alt_base}"
                try:
                    chr_str = f'chr{chrom}'

                    match = pathogenicity_class_df_hg38[
                                (pathogenicity_class_df_hg38['chrom'] == chr_str) &
                                (pathogenicity_class_df_hg38['pos'] == variant_pos) &
                                (pathogenicity_class_df_hg38['ref'] == ref_base) &
                                (pathogenicity_class_df_hg38['alt'] == alt_base)
                            ]
                    
                    if not match.empty:
                        try:
                            pathogenicity_class = match.iloc[0]['am_class'] 
                        except (IndexError, ValueError) as e:
                            warnings.warn(f"Failed to get pathogenicity class for {gene} variant {var_id} (SNP {ref_base} -> {alt_base}): {e}")
                            pathogenicity_class = None
                    else:
                        if DEBUG:
                            print(f"No pathogenicity class found for: {chrom}, {variant_pos}, {ref_base}, {alt_base}")
                        pathogenicity_class = None

                    preds[gene][var_id][snp_key] = pathogenicity_class

                except Exception as e:
                    warnings.warn(f"Failed to process prediction for {gene} variant {var_id} (SNP {ref_base} -> {alt_base}): {e}")
                    preds[gene][var_id][snp_key] = None

    print(preds)
    time.sleep(0.3)  # courteous pause

    return {
        "alphamissense_predictions": preds}



