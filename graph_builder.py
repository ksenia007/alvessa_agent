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
from entity_extraction import entity_extraction_node, has_genes
from tool_agent_node import run_async_sync, select_tools, tool_invoke
from conditioned_claude import conditioned_claude_node
from verify import verify_evidence_node

MAX_ATTEMPTS: int = 3  # verification retries
MAX_TOOL_RESELECT: int = 1 # tool re-selection attempts before final output

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


def build_graph(mc_setup: bool = False) -> Callable[[State], State]:
    """Return a compiled LangGraph ready for invocation (diagram export is best‑effort)."""
    g = StateGraph(State)

    # Nodes
    g.add_node("select_tools", select_tools)
    g.add_node("tool_invoke", run_async_sync(tool_invoke))
    g.add_node("claude", conditioned_claude_node)
    g.add_node("verify", verify_evidence_node) 

    # Edges
    g.set_entry_point("select_tools")
    g.add_edge("select_tools", "tool_invoke")
    

    
    if not mc_setup:
        # in non-MC we let the model to re-select tools before final output  
        g.add_conditional_edges(
            "tool_invoke",
            lambda s: (
                "update_tools"
                if s.get("tool_updates", 0) < MAX_TOOL_RESELECT
                else "finalized_state"
            ),
            {"update_tools": "select_tools", "finalized_state": "claude"},
        )
        # non-MC setup: add verifier
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
    else:
        # If verification is not enabled, skip it
        g.add_edge("tool_invoke", "claude")
        g.add_edge("claude", END)

    
    compiled = g.compile()

    # Best‑effort diagram printing
    try:
        # Newer pattern (sometimes on the compiled object)
        png = compiled.draw_mermaid_png()
        with open("graph_diagram.png", "wb") as f:
            f.write(png)
    except AttributeError:
        try:
            # Older pattern via an internal graph on compiled
            png = compiled.get_graph().draw_mermaid_png()  # may exist in some versions
            with open("graph_diagram.png", "wb") as f:
                f.write(png)
        except Exception:
            try:
                # Fall back to Mermaid text (you can render externally)
                mermaid = getattr(compiled, "draw_mermaid", None)
                if callable(mermaid):
                    mm = compiled.draw_mermaid()
                else:
                    mm = compiled.get_graph().draw_mermaid()  # older fallback
                with open("graph_diagram.mmd", "w") as f:
                    f.write(mm)
            except Exception:
                print('UNABLE TO UPDATE THE DIAGRAM')
                pass

    return compiled


