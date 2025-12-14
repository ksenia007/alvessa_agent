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
from src.state import State 
import warnings
import datetime

from src.tools.base import Node
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
    gene_objs = state.get("gene_entities", {}).copy()

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
                    if g in gene_objs.keys():
                        gene_objs[g].add_tool("SEI")
                        gene_objs[g].update_text_summaries(
                            "*Sei: Values predicted by Sei: the Sei model clusters predicted effects on ChIP-seq, DNase, and histone mark profiles into 40 sequence classes. Each Sei class represents the inferred regulatory role of the local sequence where a variant occurs. These classes are annotation-agnostic and may be assigned to variants in both non-coding and coding (exonic) regions."
                        )
                    assigned += 1

        except Exception as e:
            warnings.warn(f"[SEI] Failed to process {var_id}: {e}")
            continue

    time.sleep(0.3)  # courteous pause
    return


NODES: tuple[Node, ...] = (
    Node(
        name="sei",
        entry_point=sei_predictions_agent,
        description="Fetches predictions of the sequence regulatory activity for given variants. This requires variant_annotations to be run first.",
        dependencies=("variant_annotations",),
    ),
)
