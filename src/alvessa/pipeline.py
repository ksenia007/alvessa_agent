"""Top-level helpers to run the Alvessa LangGraph workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Dict

from src.alvessa.workflow.graph_builder import build_graph


def run_pipeline(
    user_message: str,
    prompt: str = "",
    mc_setup: bool = False,
    *,
    output_dir: Path | str | None = None,
) -> Dict:
    """Run the LangGraph workflow for a single user message.

    Parameters
    ----------
    user_message:
        The natural-language question to answer.
    prompt:
        Optional system prompt override.
    mc_setup:
        If ``True``, run the workflow without the verifier node.
    output_dir:
        Optional directory used for graph artifacts (diagram, log files, etc.).
    """
    graph = build_graph(mc_setup=mc_setup, diagram_dir=output_dir)
    state = graph.invoke(
        {
            "messages": [{"role": "user", "content": user_message}],
            "prompt": prompt,
            "mc_setup": mc_setup,
        }
    )
    return state
