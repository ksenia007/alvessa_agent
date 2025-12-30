"""
Description:

Tool to query DisGeNet annotations with a list of gene symbols"""

from __future__ import annotations

import time

from src.state import State
from src.tools.base import Node
from typing import List
from .utils import _fetch_predictions_DisGeNet

DEBUG = True


def disgenet_agent(state: "State") -> "State":
    """
    LangGraph node that annotates each gene with DisGeNet data.

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state with the disgenet fields filled.
        
    """
    gene_entities = state.get("gene_entities") or {}

    for gene_name, gene in gene_entities.items():
        if not gene_name:
            continue

        summary_lines: List[str] = []

        try:
            disease_list = _fetch_predictions_DisGeNet(gene)

            if len(disease_list) > 0:
                summary_lines.append(f"*DisGeNet: Diseases annotated to {gene.symbol}: " + "; ".join(disease_list) + ".")
                gene.add_many_disgenet_diseases(disease_list)

            if summary_lines:
                gene.update_text_summaries(" ".join(summary_lines))
            
                gene.add_tool("DisGeNet")

        except:
            print(gene)


        time.sleep(0.3)  # courteous pause

    if DEBUG:
        print(f"[DisGeNet] Predictions fetched")

    return 

NODES: tuple[Node, ...] = (
    Node(
        name="DisGeNet",
        entry_point=disgenet_agent,
        description="Fetches information about disease annotations for the input genes from DisGeNet, which provides curated data linking human genes to a wide range of diseases, including Mendelian, complex, environmental, and rare diseases.",
    ),
)
