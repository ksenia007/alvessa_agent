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

    print(f"started... {datetime.now()}")
    preds = state.get("alphamissense_predictions", {}).copy()
    variants = state.get("dbsnp_variants", {}).copy()

    # Gracefully handle file reading errors
    try:
        pathogenicity_class_df_hg38 = pd.read_parquet('local_dbs/AlphaMissense_hg38.parquet')
    except Exception as e:
        warnings.warn(f"Failed to load required files: {e}. Cannot run AlphaMissense predictions.")
        return {**state, "alphamissense_predictions": preds}
    
    print(f"loaded am file... {datetime.now()}")
    
    snp_records = []
    for gene, gene_vars in variants.items():
        for var_id, var_data in gene_vars.items():
            all_snps = var_data.get("coordinates", [])
            for snp in all_snps:
                chrom, pos, ref, alt, assembly = snp.get("chrom"), snp.get("pos"), snp.get("ref"), snp.get("alt"), snp.get("assembly")

                if assembly and "GRCh38" in assembly:
                    if None not in (chrom, pos, ref, alt):
                        snp_records.append({
                            "gene": gene,
                            "var_id": var_id,
                            "snp_key": f"SNP:{ref}->{alt}",
                            "chrom": f"chr{chrom}",
                            "pos": pos,
                            "ref": ref,
                            "alt": alt
                        })
                    else:
                        warnings.warn(f"Missing coordinate data for {gene} variant {var_id} (SNP {ref}->{alt})")

    if not snp_records:
        return {"alphamissense_predictions": preds}

    snps_df = pd.DataFrame(snp_records)

    print(f"finished aggregating snps... {datetime.now()}")

    merged = snps_df.merge(
    pathogenicity_class_df_hg38[['chrom', 'pos', 'ref', 'alt', 'am_class']],
    on=['chrom', 'pos', 'ref', 'alt'],
    how='left'
    )

    print(f"finished merging... {datetime.now()}")

    for (gene, var_id), group in merged.groupby(["gene", "var_id"]):
        preds.setdefault(gene, {})
        preds[gene][var_id] = {
            row.snp_key: row.am_class if not pd.isna(row.am_class) else None
            for row in group.itertuples()
        }

    print(f"done... {datetime.now()}")

    print(preds)
    time.sleep(0.3)  # courteous pause

    return {
        "alphamissense_predictions": preds}


