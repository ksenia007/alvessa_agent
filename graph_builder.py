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
from tool_go_summarization import go_summarization_agent, make_go_summarization_node
from tool_humanbase import humanbase_predictions_agent
from tool_uniprot import (
    uniprot_node,
    trait_disease_extraction_node,
    trait_function_extraction_node,
    trait_GO_extraction_node,
    has_uniprot_entries,
)
from tool_gwas import gwas_associations_agent
from conditioned_claude import conditioned_claude_node
from verify import verify_evidence_node

MAX_ATTEMPTS: int = 3  # verification retries

import asyncio
from typing import Callable, Any

def run_async_sync(fn: Callable[..., Any]) -> Callable[..., Any]:
    """
    Take an async function or coroutine-returning function and return a sync version
    that runs it via asyncio.run().
    """
    def wrapper(*args, **kwargs):
        result = fn(*args, **kwargs)
        if asyncio.iscoroutine(result):
            return asyncio.run(result)
        return result
    return wrapper


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
    g.add_node("biogrid_go_summarization", make_go_summarization_node("biogrid", "word2vec"))
    g.add_node("humanbase", humanbase_predictions_agent)
    g.add_node("uniprot_base", uniprot_node)
    g.add_node("uniprot_gwas", uniprot_node)
    g.add_node("trait_disease", trait_disease_extraction_node)
    g.add_node("trait_function", trait_function_extraction_node)
    g.add_node("trait_go", trait_GO_extraction_node)
    g.add_node("gwas", gwas_associations_agent)


    g.add_edge("biogrid", "biogrid_go_summarization")    
    g.add_edge("biogrid_go_summarization", "humanbase")
    g.add_edge("humanbase", "uniprot_base")
    g.add_edge("uniprot_base", "trait_disease")
    g.add_edge("trait_disease", "trait_function")
    g.add_edge("trait_function", "trait_go")
    g.add_edge("trait_go", "gwas")
    g.add_edge("gwas", "uniprot_gwas")
    g.add_edge("uniprot_gwas", "claude")

    # Main LLM
    g.add_node("claude", conditioned_claude_node)

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
