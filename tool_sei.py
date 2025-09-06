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
import datetime
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
    print('[SEI] started...')
    variants = state.get("variant_entities", {}).copy()

    # Load resources
    try:
        seq_class_df_hg38 = pd.read_csv(
            "local_dbs/sorted.hg38.tiling.bed.ipca_randomized_300.labels.merged.bed",
            sep="\t",
            names=["chr", "start_pos", "end_pos", "seq_class"],
        )
        with open("local_dbs/seqclass.names", "r") as f:
            seq_class_names = [line.strip() for line in f]
    except Exception as e:
        warnings.warn(f"[SEI] Failed to load required files: {e}. Cannot run SEI predictions.")
        return

    # Group BED rows by chromosome for faster filtering
    bed_by_chr = {c: df for c, df in seq_class_df_hg38.groupby("chr", sort=False)}

    assigned = 0
    for var_id, var_obj in variants.items():
        try:
            locs = var_obj.get_location("GrCh38")  
            genes = var_obj.get_related_genes()
            if not locs or not genes:
                continue

            chrom = locs.get("chrom")
            pos   = locs.get("pos")
            if chrom is None or pos is None:
                continue

            chr_str = f"chr{chrom}"
            bed_chr = bed_by_chr.get(chr_str)
            if bed_chr is None:
                # no tiles for this chromosome
                continue

            # Interval lookup 
            hit = bed_chr[(bed_chr["start_pos"] <= pos) & (bed_chr["end_pos"] >= pos)]
            if hit.empty:
                seq_class_name = None
            else:
                try:
                    cls_idx = int(hit.iloc[0]["seq_class"])
                    seq_class_name = seq_class_names[cls_idx] if 0 <= cls_idx < len(seq_class_names) else None
                except Exception as e:
                    warnings.warn(f"[SEI] seq_class index issue for {var_id} @ {chr_str}:{pos}: {e}")
                    seq_class_name = None

            if seq_class_name:
                for g in genes:
                    variants[var_id].add_functional_prediction(g, "SEI", seq_class_name)
                    assigned += 1

        except Exception as e:
            warnings.warn(f"[SEI] Failed to process {var_id}: {e}")
            continue

    time.sleep(0.3)  # courteous pause
    return


