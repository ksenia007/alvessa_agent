"""
Author: Keerthana Nallamotu <kn6412@princeton.edu>
Contributors: 
Created: 2024-08-05
Updated: 2025-08-10


Description: 

Tool to extract AlphaMissense prediction for variants of interest"""

import pandas as pd
import time
from src.state import State 
import warnings
from datetime import datetime
import requests

from src.tools.base import Node

DEBUG=True

def _uniprot_to_symbol(uniprot_id: str) -> str:
    """
    Convert UniProt accession(s) to gene symbol(s) using mygene.info.

    Parameters
    ----------
    uniprot_id : str
        UniProt accession (e.g. "P04637").

    Returns
    -------
    list[str]
        List of gene symbols (may be empty if none found).
    """
    symbols = []
    try:
        r = requests.get(f"https://mygene.info/v3/query?q=uniprot:{uniprot_id}&fields=symbol")
        r.raise_for_status()
        hits = r.json().get("hits", [])
        for hit in hits:
            sym = hit.get("symbol")
            if sym:
                symbols.append(sym)
    except Exception as e:
        warnings.warn(f"Failed UniProt→symbol lookup for {uniprot_id}: {e}")
    return ','.join(symbols)

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

    variants = state.get("variant_entities", {}).copy()
    gene_objs = state.get("gene_entities", {}).copy()

    # Early return if no variants to process
    if not variants:
        if DEBUG:
            print("[AlphaMissense] No variants to process, skipping")
        return

    if DEBUG:
        print(f"[AlphaMissense] Started with {len(variants)} variants... {datetime.now()}")

    # Build SNP records FIRST, before loading the large parquet file
    snp_records = []
    for var_id, var_obj in variants.items():
        locs = var_obj.get_location('GrCh38')
        if not locs:
            continue

        chrom, pos, alts = locs.get('chrom'), locs.get('pos'), locs.get('alt')

        # Skip if missing essential data
        if not alts or not chrom or not pos:
            if DEBUG:
                warnings.warn(f"Skipping variant {var_id} - missing coordinate data")
            continue

        # Add SNP records for each alt allele
        for alt in alts:
            snp_records.append({
                "var_id": var_id,
                "snp_key": f"SNP:REF->{alt}",
                "chrom": f"chr{chrom}",
                "pos": pos,
                "alt": alt
            })

    # Early return if no processable SNPs
    if not snp_records:
        if DEBUG:
            print("[AlphaMissense] No SNP records found to process")
        return

    if DEBUG:
        print(f"[AlphaMissense] Built {len(snp_records)} SNP records from {len(variants)} variants")

    # NOW load the large parquet file - only if we have SNPs to look up
    try:
        pathogenicity_class_df_hg38 = pd.read_parquet('local_dbs/AlphaMissense_hg38.parquet')
        if DEBUG:
            print(f"[AlphaMissense] Loaded AlphaMissense parquet... {datetime.now()}")
    except Exception as e:
        warnings.warn(f"Failed to load AlphaMissense parquet: {e}. Cannot run predictions.")
        return

    snps_df = pd.DataFrame(snp_records)

    try:
        merged = snps_df.merge(
            pathogenicity_class_df_hg38[["chrom", "pos", "ref", "alt", "uniprot_id", "am_class"]],
            left_on=["chrom", "pos",  "alt"],
            right_on=["chrom", "pos", "alt"],
            how="left",
        )
    except Exception as e:
        warnings.warn(f"[AlphaMissense] Merge failed: {e}")
        return 
    # Drop lines with NaN in am_class
    merged = merged.dropna(subset=['am_class'])

    if DEBUG:
        print(f"[AlphaMissense] Finished merging, found {len(merged)} predictions... {datetime.now()}")

    # Early return if no predictions found
    if merged.empty:
        if DEBUG:
            print("[AlphaMissense] No predictions found after merge")
        return

    # Add gene column based on uniprot_IDs
    if 'gene' not in merged.columns:
        needed_conversions = merged['uniprot_id'].dropna().unique()
        if DEBUG:
            print(f"[AlphaMissense] Converting {len(needed_conversions)} UniProt IDs to gene symbols...")

        uniprot_to_symbol_map = {}
        for i, uid in enumerate(needed_conversions):
            uniprot_to_symbol_map[uid] = _uniprot_to_symbol(uid)
            # Rate limit: 1 request per 0.2s to be courteous to mygene.info API
            if i < len(needed_conversions) - 1:
                time.sleep(0.2)

        merged['gene'] = merged['uniprot_id'].map(uniprot_to_symbol_map)

    grouped = merged.groupby(['gene', 'var_id', 'snp_key'], as_index=False).agg({
        'am_class': lambda x: next(iter(set(filter(pd.notna, x))), None)
    })

    predictions_added = 0
    for row in grouped.itertuples(index=False):
        gene, var_id, snp_key, am_class = row
        if am_class is None:
            continue
        variants[var_id].add_functional_prediction(gene, 'AlphaMissense', am_class)
        predictions_added += 1

    if DEBUG:
        print(f"[AlphaMissense] Done! Added {predictions_added} predictions... {datetime.now()}")

    return 


NODES: tuple[Node, ...] = (
    Node(
        name="alphamissense",
        entry_point=alphamissense_predictions_agent,
        description="Fetches Alphamissense predicted pathogenicity classes for given variants. This requires variant_annotations to be run first.",
        dependencies=("variant_annotations",),
    ),
)
