""" 
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-06-25


Description: 

Compile the LangGraph workflow and save the diagram as PNG."""

from __future__ import annotations
from pathlib import Path

from langgraph.graph import StateGraph, END

from src.state import State
from src.alvessa.agents.entity_extraction import entity_extraction_node, has_genes
from src.alvessa.agents.tool_agent import select_tools, tool_invoke
from src.alvessa.agents.conditioned_claude import conditioned_claude_node
from src.alvessa.agents.verify import verify_evidence_node
from src.alvessa.agents.adversarial_agent import adversarial_node

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


def build_graph(mc_setup: bool = False, diagram_dir: Path | str | None = None) -> Callable[[State], State]:
    """Return a compiled LangGraph ready for invocation (diagram export is best‑effort)."""
    g = StateGraph(State)

    # Nodes
    g.add_node("select_tools", select_tools)
    g.add_node("tool_invoke", run_async_sync(tool_invoke))
    g.add_node("claude", conditioned_claude_node)
    g.add_node("verify", verify_evidence_node) 
    g.add_node("adversarial_agent", adversarial_node)

    # Edges
    g.set_entry_point("select_tools")
    g.add_edge("select_tools", "tool_invoke")
    

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
    g.add_edge("claude", "adversarial_agent")
    g.add_edge("adversarial_agent", "verify")
    g.add_edge("verify", END)
    
    compiled = g.compile()

    # Best‑effort diagram printing
    output_dir = Path(diagram_dir) if diagram_dir is not None else Path(__file__).resolve().parents[3]
    output_dir.mkdir(parents=True, exist_ok=True)

    png_path = output_dir / "graph_diagram.png"
    mmd_path = output_dir / "graph_diagram.mmd"

    try:
        # Newer pattern (sometimes on the compiled object)
        png = compiled.draw_mermaid_png()
        with png_path.open("wb") as f:
            f.write(png)
    except AttributeError:
        try:
            # Older pattern via an internal graph on compiled
            png = compiled.get_graph().draw_mermaid_png()  # may exist in some versions
            with png_path.open("wb") as f:
                f.write(png)
        except Exception:
            try:
                # Fall back to Mermaid text (you can render externally)
                mermaid = getattr(compiled, "draw_mermaid", None)
                if callable(mermaid):
                    mm = compiled.draw_mermaid()
                else:
                    mm = compiled.get_graph().draw_mermaid()  # older fallback
                with mmd_path.open("w", encoding="utf-8") as f:
                    f.write(mm)
            except Exception:
                print('UNABLE TO UPDATE THE DIAGRAM')
                pass

    return compiled
