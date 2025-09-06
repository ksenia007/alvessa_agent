"""
Author: Keerthana Nallamotu <kn6412@princeton.edu>
Contributors: 
Created: 2024-08-21
Updated: 2025-08-21


Description: 

Tool to query miRDB with a list of miRNAs"""

from __future__ import annotations
from state import State 
import pandas as pd
import requests
from datetime import datetime

DEBUG=True

def _convert_to_miRBASE(symbol):
    components = symbol.split("_")

    name = components[0]
    if len(components)>1:
        arm = components[1]
    else:
        arm = None

    output = f'miR-{name[3:]}'
    if arm:
        output+=f'-{arm[0]}p'

    return output

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
    
    organism_codes = {
        'cfa': 'dog',
        'gga': 'chicken',
        'hsa': 'human',
        'mmu': 'mouse',
        'rno': 'rat'
    }

    df = pd.read_csv(
        'local_dbs/miRDB_v6.0_prediction_result.txt',
        sep='\t',
        header=None,
        names=['miRNAID', 'geneID', 'confidence']
    )

    for gene in gene_objs.keys():
        if 'mir' not in gene.lower():
            continue

        mirbase_ID = _convert_to_miRBASE(gene)
        print(f"Processing {gene} as {mirbase_ID}")
        preds_temp = {}

        for prefix, organism_name in organism_codes.items():
            subset = df[df['miRNAID'].str.startswith(prefix) & df['miRNAID'].str.contains(mirbase_ID)]
            if subset.empty:
                print(f"No targets found for {mirbase_ID} in {organism_name}")
                continue

            targets_entrez = subset['geneID'].values

            print(f"started symbol querying ({organism_name})... {datetime.now()}")
            targets_symbols = _refseq_to_symbol(targets_entrez)
            print(f"finished symbol querying ({organism_name})... {datetime.now()}")

            preds_temp[organism_name] = list(set((targets_symbols)))
            
        gene_objs[gene].add_miRNA_targets(preds_temp)
        gene_objs[gene].add_tool("miRDB_agent")
        # create a text summary about miRNA targets
        summary = f" Using mirDB, predicted targets for {gene}, separated by organism, include: "
        for org, targets in preds_temp.items():
            summary += f"{org} - {', '.join(targets)}; "
        summary += f' End of record for {gene} |'
        gene_objs[gene].update_text_summaries(summary)


    return 
