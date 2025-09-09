"""
Author: Keerthana Nallamotu <kn6412@princeton.edu>
Contributors: 
Created: 2024-08-21
Updated: 2025-09-09


Description: 

Tool to query miRDB with a list of miRNAs"""

from __future__ import annotations
from state import State 
import pandas as pd
import requests
from datetime import datetime
import re

DEBUG=True

organism_codes = {
    'cfa': 'dog',
    'gga': 'chicken',
    'hsa': 'human',
    'mmu': 'mouse',
    'rno': 'rat'
}

def _standardize_miRNA_name(symbol: str): 
    symbol = symbol.lower().replace('_', '-')
    organism = None

    for code in organism_codes.keys():
        if symbol.startswith(code + '-'):
            organism = code
            symbol = symbol[len(code)+1:] 
            break

    pattern = re.compile(r'mir-?(\d+)(-?(5p|3p))?')

    match = pattern.match(symbol)
    if not match:
        pattern_simple = re.compile(r'(\d+)(-?(5p|3p))?')
        match = pattern_simple.match(symbol)

    if not match:
        standardized_name = symbol
    else:
        number = match.group(1)
        arm = match.group(3)  
        standardized_name = f"miR-{number}"
        if arm:
            standardized_name += f"-{arm}"

    if organism:
        standardized_name = f"{organism}-{standardized_name}"

    return organism, standardized_name


import requests

def _refseq_to_symbol(refseq_ids):
    url = "http://mygene.info/v3/query"
    all_symbols = []

    for i in range(0, len(refseq_ids), 300):
        chunk = refseq_ids[i:i+300]
        query_str = " OR ".join([f"refseq:{ID}" for ID in chunk])
        params = {
            'q': query_str,
            'fields': 'symbol',
            'size': len(chunk)
        }

        try:
            r = requests.get(url, params=params)
            r.raise_for_status()
            hits = r.json().get("hits", [])
            all_symbols.extend(hit.get('symbol') for hit in hits if 'symbol' in hit)
        except Exception as e:
            print(e)

    return all_symbols


def miRDB_agent(state: "State") -> "State":
    """
    LangGraph node that annotates each miRNA with target genes.

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state with the `"mirDB_targets"` field filled.
        
    """

    #preds = state.get("mirDB_targets", {}).copy()
    gene_objs = state.get('gene_entities', {})

    df = pd.read_csv(
        'local_dbs/miRDB_v6.0_prediction_result.txt',
        sep='\t',
        header=None,
        names=['miRNAID', 'geneID', 'confidence']
    )

    for gene in gene_objs.keys():
        if 'mir' not in gene.lower():
            continue

        organism, standardized_miRNA = _standardize_miRNA_name(gene)
        if DEBUG:
            print(f"Processing {gene} as {standardized_miRNA} (organism={organism})")

        preds_temp = {}

        organisms_to_check = [organism] if organism else list(organism_codes.keys())

        for org_code in organisms_to_check:
            organism_name = organism_codes[org_code]
            
            if standardized_miRNA.startswith(org_code + '-'):
                miRNA_to_match = standardized_miRNA[len(org_code)+1:]
            else:
                miRNA_to_match = standardized_miRNA

            subset = df[
                df['miRNAID'].str.startswith(org_code) &
                df['miRNAID'].str.contains(miRNA_to_match)
            ]

            if subset.empty:
                if DEBUG:
                    print(f"No targets found for {standardized_miRNA} in {organism_name}")
                continue

            targets_entrez = subset['geneID'].values

            if DEBUG:
                print(f"started symbol querying ({organism_name})... {datetime.now()}")
            targets_symbols = _refseq_to_symbol(targets_entrez)
            if DEBUG:
                print(f"finished symbol querying ({organism_name})... {datetime.now()}")

            preds_temp[organism_name] = list(set(targets_symbols))
            
        gene_objs[gene].add_miRNA_targets(preds_temp)
        gene_objs[gene].add_tool("miRDB_agent")

        # create a text summary about miRNA targets
        summary = f"All computationally predicted gene targets for {gene} (from the miRDB database), separated by organism, include: "
        for org, targets in preds_temp.items():
            summary += f"{org} - {', '.join(targets)}; "
        summary += f" End of record for {gene} |"
        gene_objs[gene].update_text_summaries(summary)

    return 
