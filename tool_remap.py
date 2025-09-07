
"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-09-07
Updated: 


Description: 

ReMap tool to get bindings in front of TSS """

from __future__ import annotations
import re
from typing import Any, Dict, List, Optional, Tuple
import pandas as pd
import numpy as np
# Agent-compatible node
from state import State 

DEBUG=True



def remap_crm_agent(state: "State", window_width: int =1000) -> "State":
    """
    Annotate each Gene with nearby ReMap 2022 CRM peaks (hg38).
    For each gene:
      - get (chrom, start, end, strand) on GrCh38
      - TSS := start if is_positive_strand(strand) else end
      - find peaks within window width bp of TSS
      - from BED col4 collect list of 'genes that bind there' (CRM members)
      - write a concise summary to gene.update_text_summaries(...)
    
    NOTE: We also add list to the Gene object, but the main output is the text summary.
    """

    genes = (state.get("gene_entities") or {}).copy() # get all the Gene objects
    if not genes:
        return

    #  --- related to this specific tool - loading info
    # 1) Load ReMap CRM BED (expect ≥ 4 columns)
    try:
        bed = pd.read_csv(
            "local_dbs/remap2022_crm_macs2_hg38_v1_0.bed",
            sep="\t",
            header=None,
            comment="#",
            dtype={0: str, 1: int, 2: int, 3: str},
            usecols=[0, 1, 2, 3],  # chrom, start, end, name(col4)
            names=["chrom", "start", "end", "crm_genes"],
        )
    except Exception as e:
        if DEBUG:
            print(f"[ReMap] Could not load local ReMap CRM BED file.")
        return

    # 2) Group by chromosome for quick slicing
    bed_by_chr = {c: df.reset_index(drop=True) for c, df in bed.groupby("chrom", sort=False)}
    
    # Done loading and pre-processing

    # 3) Walk genes and annotate
    annotated = 0
    for gsym, gobj in genes.items():
        # Get hg38 coordinates
        try:
            chrom, gstart, gend, strand = gobj.get_location()
        except Exception:
            chrom, gstart, gend, strand = None, None, None, None

        if chrom is None or gstart is None or gend is None or strand is None:
            print(f"[ReMap] Skipping {gsym}: missing chrom/start/end/strand info.")
            continue

        # Normalize chromosome label to match BED (with 'chr' prefix)
        chr_str = f"chr{chrom}" if not str(chrom).startswith("chr") else str(chrom)
        chr_bed = bed_by_chr.get(chr_str)
        if chr_bed is None:
            continue

        # TSS based on strand
        if strand == "+":
            tss = gstart
        elif strand == "-":
            tss = gend
        else:
            if DEBUG:
                print(f"[ReMap] Skipping {gsym}: unrecognized strand '{strand}'.")
            continue
        
        # +/- window_width bp window
        win_start = tss - window_width
        win_end   = tss + window_width

        # Overlap condition: peak_end >= win_start AND peak_start <= win_end
        hits = chr_bed[(chr_bed["end"] >= win_start) & (chr_bed["start"] <= win_end)]
        if hits.empty:
            continue

        # Parse CRM “genes that bind there” from col4 (comma/semicolon/pipe)
        binders = []
        for name in hits["crm_genes"]:
            if not isinstance(name, str):
                continue
            # split on common delimiters for CRM member lists
            parts = re.split(r"[;,|]", name)
            for p in parts:
                p = p.strip()
                if p:
                    binders.append(p)

        # Dedup and keep it short
        uniq_binders = list(dict.fromkeys(binders))
        n_peaks = len(hits)
        n_bind  = len(uniq_binders)
        list_bind = ", ".join(uniq_binders)

        # Write one concise line to the Gene
        # IMPORTANT: this is the main output of this tool
        gobj.update_text_summaries(
            f"[ReMap Cis Regulatory Modules] Which TFs are predicted to bind for {gsym} {chr_str}:{tss} (+/-{window_width}bp): {n_peaks} peaks; "
            f"Total unique binders {n_bind}: {list_bind}"
        )
        annotated += 1
        
        # Also add to Gene object, to field binding_peaks
        dict_info = {
            'n_peaks': n_peaks,
            'n_unique_binders': n_bind,
            'unique_binders': uniq_binders,
        }
        gobj.add_binding_peaks(dict_info)
        
        # don't forget to note that we used this tool
        gobj.add_tool("remap_crm_agent")

    if DEBUG:
        print(f"[ReMap] Annotated {annotated} genes with nearby CRM peaks.")

    return
