"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-06-25
Updated: 2025-06-26


Description: 

LLM node that answers the user question using compiled context: main reasoning LLM"""

from __future__ import annotations
import json
from typing import Dict, List, Any

from claude_client import claude_call
from config import CONDITIONED_MODEL, DEBUG, N_CHARS
from state import State


def conditioned_claude_node(state: "State") -> "State":
    """
    Build a compact JSON context and ask Claude for an answer.

    Returns
    -------
    State
        Updated state with assistant message, context block and parsed JSON.
    """
    if DEBUG:
        print("[conditioned_claude_node] preparing context block...")
    
    # Build CONTEXT payload
    gene_payload: List[Dict[str, Any]] = []
    for g in state.get("genes", []):
        gene_info: Dict[str, Any] = {"gene": g}

        diseases = state.get("gene_disease_traits", {}).get(g, [])
        if diseases:
            gene_info["diseases"] = diseases

        humanbase_hits = state.get("humanbase_predictions", {}).get(g, [])
        if humanbase_hits:
            terms = [hit["term"] for hit in humanbase_hits if "term" in hit]
            if terms:
                gene_info["functions"] = terms[:30]

        gene_payload.append(gene_info)

    
    context_block: str = json.dumps(gene_payload, separators=(",", ":"))
    if DEBUG:
        print("[conditioned_claude_node] context length:", len(context_block))
    if len(context_block) > N_CHARS:
        context_block = context_block[:N_CHARS] + "...<truncated>"
    print(context_block)
    
    # Call Claude
    system_msg: str = (
        "You are a biology data analyst. Answer strictly with facts you can "
        "point to inside CONTEXT. Respond only with JSON with keys answer and evidence. Ensure proper JSON format. "
        "The 'evidence' field must always be a list of short strings."
    )
    user_question: str = state["messages"][-1]["content"]

    raw = claude_call(
        model=CONDITIONED_MODEL,
        temperature=0,
        max_tokens=5000,
        system=system_msg,
        messages=[{"role": "user", "content": f"User asked: {user_question}\n\nCONTEXT:\n{context_block}"}],
    )

    if hasattr(raw.content[0], "text"):
        llm_resp: str | Dict[str, Any] = raw.content[0].text.strip()
    else:
        llm_resp = raw.content[0]

    try:
        parsed_resp: Dict[str, Any]
        if isinstance(llm_resp, dict):
            parsed_resp = llm_resp
        else:
            parsed_resp = json.loads(llm_resp)
    except Exception as exc:  # Fall back to raw JSON string
        raise ValueError(
            f"Failed to parse LLM response as JSON. Response was:\n{llm_resp}\nError: {exc}"
        ) from exc

    return {
        **state,
        "messages": state["messages"] + [{"role": "assistant", "content": llm_resp}],
        "context_block": context_block,
        "llm_json": parsed_resp,
    }
