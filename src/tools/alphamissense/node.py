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

    print(f"[AlphaMissense] started... {datetime.now()}")
    # preds = state.get("alphamissense_predictions", {}).copy()
    variants = state.get("variant_entities", {}).copy()
    gene_objs = state.get("gene_entities", {}).copy()

    # Gracefully handle file reading errors
    try:
        pathogenicity_class_df_hg38 = pd.read_parquet('local_dbs/AlphaMissense_hg38.parquet')
    except Exception as e:
        warnings.warn(f"Failed to load required files: {e}. Cannot run AlphaMissense predictions.")
        return 
    
    
    print(f"[AlphaMissense] loaded am file... {datetime.now()}")
    
    snp_records = []
    for var_id, var_obj in variants.items():
        locs = var_obj.get_location('GrCh38')
        if not locs:
            continue
        related_genes = var_obj.get_related_genes()

        chrom, pos,  alts = locs.get('chrom'), locs.get('pos'),  locs.get('alt')
        if not alts:
            warnings.warn(f"Skipping variant {var_id} ")
            continue
        if locs:
            if not alts: 
                continue
            for alt in alts:
                snp_records.append({
                    "var_id": var_id,
                    "snp_key": f"SNP:REF->{alt}",
                    "chrom": f"chr{chrom}",
                    "pos": pos,
                    "alt": alt
                })
        else:
            warnings.warn(f"Missing coordinate data")
            raise ValueError(f"Missing coordinate data for variant")


    if not snp_records:
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
    # drop lines with NaN in am_class
    merged = merged.dropna(subset=['am_class'])
    print(f"[AlphaMissense] finished merging... {datetime.now()}")
    # add gene column based on uniprot_IDs
    if 'gene' not in merged.columns:
        needed_conversions = merged['uniprot_id'].dropna().unique()
        uniprot_to_symbol_map = {uid: _uniprot_to_symbol(uid) for uid in needed_conversions}
        merged['gene'] = merged['uniprot_id'].map(uniprot_to_symbol_map)

    grouped = merged.groupby(['gene', 'var_id', 'snp_key'], as_index=False).agg({
        'am_class': lambda x: next(iter(set(filter(pd.notna, x))), None)
    })

    for row in grouped.itertuples(index=False):
        gene, var_id, snp_key, am_class = row
        if am_class is None:
            continue
        variants[var_id].add_functional_prediction(gene, 'AlphaMissense', am_class)

    print(f"[AlphaMissense] done... {datetime.now()}")

    time.sleep(0.3)  # courteous pause

    return 


NODES: tuple[Node, ...] = (
    Node(
        name="alphamissense",
        entry_point=alphamissense_predictions_agent,
        description="Fetches Alphamissense predicted pathogenicity classes for given variants. This requires variant_annotations to be run first.",
        dependencies=("variant_annotations",),
    ),
)
