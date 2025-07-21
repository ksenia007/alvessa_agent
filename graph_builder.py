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
from tool_agent_node import select_tools_and_run_dynamic, run_async_sync
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
    g.add_node("select_tools", run_async_sync(select_tools_and_run_dynamic))
    g.set_entry_point("extract_genes")
    g.add_conditional_edges(
        "extract_genes",
        has_genes,
        {True: "select_tools", False: "claude"},
    )

    g.add_edge("select_tools", "claude")
    
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
