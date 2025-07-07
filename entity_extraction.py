"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-06-25
Updated: 2025-06-26


Description: 

Nodes to extract gene symbols from user input."""

from __future__ import annotations
from typing import List

from claude_client import claude_call
from config import DEBUG, GENE_EXTRACT_MODEL
from state import State


def gene_extraction_node(state: "State") -> "State":
    """
    Identify HGNC gene symbols mentioned in the last user message.

    Returns
    -------
    State
        Same state with a new `"genes"` list.
    """
    user_input: str = state["messages"][-1]["content"]
    system_message: str = (
        "Extract gene symbols from the message. Extract gene names **only if they appear verbatim in the input**. "
        "Reply with a comma-separated list of gene names only, no extra words."
    )
    response = claude_call(
        model=GENE_EXTRACT_MODEL,
        max_tokens=50,
        temperature=0,
        system=system_message,
        messages=[{"role": "user", "content": user_input}],
    )
    genes: List[str] = [g.strip() for g in response.content[0].text.split(",") if g.strip()]
    if DEBUG:
        print("[gene_extraction_node] extracted:", genes)
    return {**state, "genes": genes}


def has_genes(state: "State") -> bool:
    """
    Edge-condition helper for LangGraph: returns `True` if any genes were found.
    """
    if DEBUG:
        print("[has_genes] genes present?", bool(state.get("genes")))
    return bool(state.get("genes"))
