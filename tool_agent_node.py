"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-07-15
Updated: 2025-07-21


Description: 

Selecting tools"""

from state import State

from tool_humanbase import humanbase_predictions_agent
from tool_biogrid import bioGRID_predictions_agent
from tool_uniprot import uniprot_node
from tool_gwas import gwas_associations_agent
from tool_descriptions import TOOL_CATALOG, TOOL_FN_MAP, EXAMPLE_TOOL_SELECTION
from config import TOOL_SELECTOR_MODEL, N_CHARS
from claude_client import claude_call
from typing import Any, Dict, List
import inspect
import asyncio
import ast
from conditioned_claude import create_context_block

DEBUG = True


def format_state_for_prompt(state: State) -> str:
    """Format the state for the tool selection prompt."""

    question = state["messages"][-1]["content"]
    catalog = "\n".join(f"- {name}: {desc}" for name, desc in TOOL_CATALOG.items())

    return f"""You are an assistant deciding which tools to use to answer a biomedical question. User question: \"\"\"{question}\"\"\" \n\n Available tools: {catalog}.\n\n Examples workflows: {EXAMPLE_TOOL_SELECTION} \n Which tools should be called, and in what order? Respond *ONLY* with a Python list of tool names. Example: ["humanbase_functions", "uniprot_base"] or ["humanbase_functions", "uniprot_base", "query_gwas_by_gene"] or ["query_gwas_by_gene", "BioGRID"]"""

def _safe_merge(acc: dict, out: dict) -> None:
    """Merge tool output into acc without mutating caller state."""
    for k, v in list(out.items()):
        if k == "messages":
            continue  # never propagate messages

        # union semantics for genes
        if k == "genes":
            v = v or []
            prev = acc.get("genes", [])
            acc["genes"] = list(dict.fromkeys(list(prev) + list(v)))
            continue

        # shallow-merge dicts (last write wins)
        if isinstance(v, dict) and isinstance(acc.get(k), dict):
            merged = dict(acc[k])
            merged.update(v)
            acc[k] = merged
            continue

        # de-duplicate lists if both sides are lists
        if isinstance(v, list) and isinstance(acc.get(k), list):
            acc[k] = list(dict.fromkeys(acc[k] + v))
            continue

        # default: overwrite
        acc[k] = v


def tool_invoke(state: "State") -> "State":
    """Invoke selected tools and agents, returning ONLY updates (no in-place mutation)."""
    selected_tools = list(state.get("use_tools", []))
    used = set(state.get("used_tools", []))

    # filter out tools already used
    selected_tools = [t for t in selected_tools if t not in used]
    if DEBUG:
        print(f"[TOOL INVOKE] Selected tools: {selected_tools}")
        print(f"[TOOL INVOKE] Already used tools: {sorted(used)}")

    # no new tools - return a delta that clears use_tools (optional) and nothing else
    if not selected_tools:
        return {"use_tools": []}

    acc: dict = {}  # all updates accumulated here

    # run tools sequentially against a read-only view
    for name in selected_tools:
        fn = TOOL_FN_MAP.get(name)
        if not fn:
            continue

        print(f"[TOOL RUN] → {name}", flush=True)

        # compose a view for the tool (original state + updates so far), no mutation
        view = dict(state)
        view.update(acc)

        out = fn(view)
        if not isinstance(out, dict):
            used.add(name)
            continue

        # sanitize & merge
        out = dict(out)  # shallow copy
        out.pop("messages", None)
        _safe_merge(acc, out)

        used.add(name)

    # book-keeping updates
    acc["used_tools"] = sorted(used)
    acc["use_tools"] = []  # clear the plan so it isn't re-applied

    return acc


def select_tools(state: "State") -> "State":
    """Intent recognition and tool selection."""
    used_tools = set(state.get("used_tools", []))

    if used_tools:
        if DEBUG:
            print("[TOOL SELECTION] Selecting additional tools based on current state.")
        system_msg_base = format_state_for_prompt(state)
        current_context_block = create_context_block(state)
        used_tools_str = ", ".join(sorted(used_tools)) if used_tools else "None"
        system_msg = (
            f"{system_msg_base}\n\n"
            f"Already used tools: {used_tools_str}\n"
            f"Do not repeat tools already used.\n"
            f"IMPORTANT: if no additional tools are needed, return [] with no explanation.\n"
            f"Overall, NO explanation, downatream parser can only work with a list of tools.\n\n"
            f"Current context block:\n{current_context_block}\n\n"
        )
        if len(system_msg) > N_CHARS:
            system_msg = system_msg[:N_CHARS] + "...<truncated>"
        tool_updates = state.get("tool_updates", 0) + 1
    else:
        system_msg = format_state_for_prompt(state)
        tool_updates = state.get("tool_updates", 0)

    if DEBUG:
        print(f"[TOOL SELECTION PROMPT] {system_msg}")

    try:
        completion = claude_call(
            model=TOOL_SELECTOR_MODEL,
            max_tokens=256,
            temperature=0,
            messages=[{"role": "user", "content": system_msg}],
        )
        if DEBUG:
            print(f"[TOOL SELECTION RESPONSE] {completion}")
        tool_response = completion.content[0].text.strip()
    except Exception as e:
        print(f"[CLAUDE TOOL SELECTION ERROR] {e}")
        tool_response = "[]"

    try:
        selected_tools = ast.literal_eval(tool_response)
        assert isinstance(selected_tools, list)
    except Exception as e:
        print(f"[TOOL SELECTION ERROR] Could not parse: {tool_response}\n{e}")
        selected_tools = []

    selected_tools = [t for t in selected_tools if t in TOOL_FN_MAP and t not in used_tools]

    updates = {
        "use_tools": selected_tools,
        "tool_updates": tool_updates + (1 if used_tools else 0),
    }
    return updates

def run_async_sync(fn):
    def wrapper(*args, **kwargs):
        async def run_coroutine():
            return await fn(*args, **kwargs)
        return asyncio.run(run_coroutine())
    return wrapper