""" 
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-06-25
Updated: 2025-06-26


Description: 

Compile the LangGraph workflow and save the diagram as PNG."""

from __future__ import annotations
from langgraph.graph import StateGraph, END
from typing import Callable

from state import State
from entity_extraction import gene_extraction_node, has_genes
from tool_biogrid import bioGRID_predictions_agent
from tool_humanbase import humanbase_predictions_agent
from tool_uniprot import (
    uniprot_node,
    trait_disease_extraction_node,
    trait_function_extraction_node,
    trait_GO_extraction_node,
    has_uniprot_entries,
)
from conditioned_claude import conditioned_claude_node
from verify import verify_evidence_node

MAX_ATTEMPTS: int = 3  # verification retries


def build_graph() -> Callable[[State], State]:
    """Return a compiled LangGraph ready for invocation."""
    g = StateGraph(State)

    # Entry + gene extraction
    g.add_node("extract_genes", gene_extraction_node)
    g.set_entry_point("extract_genes")
    g.add_conditional_edges(
        "extract_genes",
        has_genes,
        {True: "biogrid", False: "claude"},
    )

    # Tool nodes
    g.add_node("biogrid", bioGRID_predictions_agent)
    g.add_node("humanbase", humanbase_predictions_agent)
    g.add_node("uniprot", uniprot_node)
    g.add_node("trait_disease", trait_disease_extraction_node)
    g.add_node("trait_function", trait_function_extraction_node)
    g.add_node("trait_go", trait_GO_extraction_node)

    g.add_edge("biogrid", "humanbase")
    g.add_edge("humanbase", "uniprot")
    g.add_edge("uniprot", "trait_disease")
    g.add_edge("trait_disease", "trait_function")
    g.add_edge("trait_function", "trait_go")

    # Main LLM
    g.add_node("claude", conditioned_claude_node)
    g.add_edge("trait_go", "claude")  # normal path

    # Verification loop
    g.add_node("verify", verify_evidence_node)
    g.add_edge("claude", "verify")

    g.add_conditional_edges(
        "verify",
        lambda s: (
            "retry"
            if s.get("verification") == "fail" and s.get("verify_attempts", 0) < MAX_ATTEMPTS
            else "done"
        ),
        {"retry": "claude", "done": END},
    )

    g = g.compile()
    
    # Save the graph diagram as PNG
    png_bytes = g.get_graph().draw_mermaid_png()
    with open("graph_diagram.png", "wb") as f:
        f.write(png_bytes)
    return g
