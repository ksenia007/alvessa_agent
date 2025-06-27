"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-06-25
Updated: 2025-06-26


Description: 

Self-verification node that checks if the answer matches context."""

from __future__ import annotations
import json
from typing import Any, Dict

from claude_client import claude_call
from config import VERIFY_MODEL, DEBUG
from state import State


def verify_evidence_node(state: "State") -> "State":
    """
    Ask model (Claude-Haiku) to judge whether the answer’s evidence is supported.

    Returns
    -------
    State
        Messages possibly annotated with failure-info and a `"verification"` key.
    """
    if DEBUG:
        print("[verify_evidence_node] preparing verification...")

    reply: Dict[str, Any] = state["llm_json"]
    context: str = state["context_block"]
    question: str = state["messages"][-1]["content"]
    answer: str = reply.get("answer", "")
    evidence_list = reply.get("evidence", [])

    system_msg: str = (
        "You are a meticulous fact-checker. Decide whether ANSWER is fully "
        "supported by the listed EVIDENCE within CONTEXT. Reply only with either\n"
        '{"verdict":"pass"}\n  or\n{"verdict":"fail","reason":"<brief>"}'
    )

    verify_prompt = (
        f"QUESTION:\n{question}\n\n"
        f"ANSWER:\n{answer}\n\n"
        f"EVIDENCE:\n{json.dumps(evidence_list, indent=2)}\n\n"
        f"CONTEXT:\n{context}"
    )

    raw = claude_call(
        model=VERIFY_MODEL,
        temperature=0,
        max_tokens=300,
        system=system_msg,
        messages=[{"role": "user", "content": verify_prompt}],
    ).content[0].text

    try:
        parsed: Dict[str, Any] = raw if isinstance(raw, dict) else json.loads(raw)
        verdict: str = parsed.get("verdict", "fail")
    except Exception:
        verdict = "fail"
        parsed = {"reason": "LLM output not valid JSON.", "raw_output": raw}

    if DEBUG:
        print("[verify_evidence_node]", verdict, parsed.get("reason", ""))

    msgs = list(state["messages"])
    if verdict == "fail":
        last = dict(msgs[-1])
        last["content"] += (
            "\n\nVERIFICATION FAILED: LLM check failed.\n"
            f"Reason: {parsed.get('reason', 'unknown')}"
        )
        msgs[-1] = last

    return {"messages": msgs, "verification": verdict}
