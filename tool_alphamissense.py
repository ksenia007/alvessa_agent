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

    pathogenicity_class_df_hg38 = pd.read_csv('local_dbs/AlphaMissense_hg38.tsv', skiprows=4, names = ['chrom', 'pos', 'ref', 'alt', 'genome', 'uniprot_id', 'transcript_id', 'protein_variant', 'am_pathogenicity', 'am_class'], sep = '\t')

    state_all_snps = {}

    for gene, gene_vars in variants.items():
        state_all_snps[gene] = {}
        for var, var_data in gene_vars.items():
            all_snps = var_data['coordinates']
            if 'assembly_filter' in var_data:
                assembly = var_data['assembly_filter']
            else:
                continue

            for snp in all_snps:
                chrom = snp['chrom']
                pos = snp['pos']
                ref_base = snp['ref']
                alt_base = snp['alt']
                state_all_snps[gene][var] = [chrom, pos, ref_base, alt_base, assembly]

    for gene, variants in state_all_snps.items():
        preds[gene] = {}
        for var_id, (chrom, variant_pos, ref_base, alt_base, assembly) in variants.items():

            chr_str = f'chr{chrom}'

            if 'GRCh38' in assembly:
                pathogenicity_class_df = pathogenicity_class_df_hg38
            else:
                preds[gene][var_id] = None
                continue

            match = pathogenicity_class_df[
                        (pathogenicity_class_df['chrom'] == chr_str) &
                        (pathogenicity_class_df['pos'] == variant_pos) &
                        (pathogenicity_class_df['ref'] == ref_base) &
                        (pathogenicity_class_df['alt'] == alt_base)
                    ]
            
            if not match.empty:
                pathogenicity_class = match.iloc[0]['am_class'] 
            else:
                pathogenicity_class = None
                
            preds[gene][var_id] = pathogenicity_class

    print(preds)
    time.sleep(0.3)  # courteous pause

    return {
        **state,
        "alphamissense_predictions": preds}


