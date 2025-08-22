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
import requests

DEBUG=True

def _symbol_to_uniprot(gene):

    all_symbols = []
    try:
        r = requests.get(f'https://mygene.info/v3/query?q={gene}&fields=uniprot')
        r.raise_for_status()
        hits = r.json().get("hits", [])
        for hit in hits:
            try:
                uniprot = hit.get("uniprot", {})
                if "Swiss-Prot" in uniprot:
                    all_symbols.append(uniprot["Swiss-Prot"])
            except Exception as inner_e:
                warnings.warn(
                    f"Skipping malformed UniProt entry for gene {gene}: {inner_e}"
                )    
    except Exception as e:
        print(e)

    return all_symbols

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

    print(f"[AlphaMissense] started... {datetime.now()}")
    preds = state.get("alphamissense_predictions", {}).copy()
    variants = state.get("dbsnp_variants", {}).copy()

    # Gracefully handle file reading errors
    try:
        pathogenicity_class_df_hg38 = pd.read_parquet('local_dbs/AlphaMissense_hg38.parquet')
    except Exception as e:
        warnings.warn(f"Failed to load required files: {e}. Cannot run AlphaMissense predictions.")
        return {**state, "alphamissense_predictions": preds}
    
    print(f"[AlphaMissense] loaded am file... {datetime.now()}")
    
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
                            "uniprot_IDs": _symbol_to_uniprot(gene),
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

    print(f"[AlphaMissense] finished aggregating snps... {datetime.now()}")

    if "uniprot_IDs" in snps_df.columns:
        snps_df["uniprot_IDs"] = snps_df["uniprot_IDs"].apply(lambda x: x if x else [None])
        snps_exploded = snps_df.explode("uniprot_IDs")
    else:
        warnings.warn("[AlphaMissense] No UniProt column generated.")
        return {"alphamissense_predictions": preds}

    try:
        merged = snps_exploded.merge(
            pathogenicity_class_df_hg38[["chrom", "pos", "ref", "alt", "uniprot_id", "am_class"]],
            left_on=["chrom", "pos", "ref", "alt", "uniprot_IDs"],
            right_on=["chrom", "pos", "ref", "alt", "uniprot_id"],
            how="left",
        )
    except Exception as e:
        warnings.warn(f"[AlphaMissense] Merge failed: {e}")
        return {"alphamissense_predictions": preds}


    print(f"[AlphaMissense] finished merging... {datetime.now()}")

    grouped = merged.groupby(['gene', 'var_id', 'snp_key'], as_index=False).agg({
        'am_class': lambda x: next(iter(set(filter(pd.notna, x))), None)
    })

    for row in grouped.itertuples(index=False):
        gene, var_id, snp_key, am_class = row
        preds.setdefault(gene, {})
        preds[gene].setdefault(var_id, {})
        preds[gene][var_id][snp_key] = am_class

    print(f"[AlphaMissense] done... {datetime.now()}")

    print("[AlphaMissense] Predictions: ", preds)
    time.sleep(0.3)  # courteous pause

    return {
        "alphamissense_predictions": preds}


