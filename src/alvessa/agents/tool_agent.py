"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-07-15
Updated: 2025-07-21


Description: 

Selecting tools"""

from src.state import State

from src.tools.humanbase.node import humanbase_predictions_agent
from src.tools.biogrid.node import bioGRID_predictions_agent
from src.tools.uniprot.node import uniprot_node
from src.tools.gwas.node import gwas_associations_agent
from src.tools.registry import TOOL_CATALOG, TOOL_FN_MAP, EXAMPLE_TOOL_SELECTION
from src.config import TOOL_SELECTOR_MODEL, N_CHARS
from src.alvessa.clients.claude import claude_call
from typing import Any, Dict, List
import inspect
import asyncio
import ast
from src.alvessa.agents.conditioned_claude import create_context_block

DEBUG = True


def format_state_for_prompt(state: State) -> str:
    """Format the state for the tool selection prompt."""

    question = state["messages"][-1]["content"]
    catalog = "\n".join(f"- {name}: {desc}" for name, desc in TOOL_CATALOG.items())

    return f"""You are an assistant deciding which tools to use to answer a biomedical question. User question: \"\"\"{question}\"\"\" \n\n If the question mentioned a database or tool and we do not have it, use the most similar available one. Available tools: {catalog}.\n\n Examples outputs look like this: {EXAMPLE_TOOL_SELECTION}"""

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
            merged = acc[k] + v
            try:
                # Works for hashable elements (str, int, tuples, etc.)
                acc[k] = list(dict.fromkeys(merged))
            except TypeError:
                if DEBUG:
                    print("[SAFE MERGE] list-of-dicts (or other unhashables) for key:", k)
                    print("  acc[k] type:", type(acc[k][0]) if acc[k] else None)
                    print("  v type:", type(v[0]) if v else None)
                # Elements are unhashable (e.g. dicts) – just concatenate
                acc[k] = merged
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

   
    view = dict(state)
    acc: dict = {}  # all updates accumulated here

    # run tools sequentially against a read-only view
    for name in selected_tools:
        print('*NAME*', name)
        view.update(acc)
        print(view)
        fn = TOOL_FN_MAP.get(name)
        if not fn:
            # check twosample_mr_agent special case
            if name.startswith("twosample_mr_agent"):
                print('[TOOL RUN] Special case: twosample_mr_agent with parameters')
                fn = TOOL_FN_MAP.get("twosample_mr_agent")
                if fn:
                    # parse parameters from the name, e.g. twosample_mr_agent-EXPOSURE-OUTCOME
                    parts = name.split("-")
                    if len(parts) == 3:
                        exposure, outcome = parts[1], parts[2]
                        if DEBUG:
                            print(f"[TOOL RUN] → {name} (exposure={exposure}, outcome={outcome})")
                        print('view in twosample_mr_agent:', view)
                        out = fn(view, exposure_gwas=exposure, outcome_gwas=outcome)
                        if out is None:
                            used.add(name)
                            continue
                        out = dict(out)
                        out.pop("messages", None)
                        _safe_merge(acc, out)

                        used.add(name)
                        continue 
                    if len(parts) == 4:
                        # also pass in optional gene name
                        exposure, outcome, gene = parts[1], parts[2], parts[3]
                        if DEBUG:
                            print(f"[TOOL RUN] → {name} (exposure={exposure}, outcome={outcome}, gene={gene})")
                        out = fn(view, exposure_gwas=exposure, outcome_gwas=outcome, gene_name=gene)
                        if out is None:
                            used.add(name)
                            continue
                        out = dict(out) # shallow copy
                        out.pop("messages", None)
                        _safe_merge(acc, out)
                        used.add(name)
                        view.update(acc)

                        continue    
            else:
                print(f"[TOOL RUN] Skipping unknown tool: {name}", flush=True)
                used.add(name)

                continue       
        else:
            print(f"[TOOL RUN] → {name}", flush=True)

            # compose a view for the tool (original state + updates so far), no mutation
        
            out = fn(view)
            if not isinstance(out, dict):
                used.add(name)
                continue

            # sanitize & merge
            out = dict(out)  # shallow copy
            out.pop("messages", None)
            _safe_merge(acc, out)

            used.add(name)
            view.update(acc)

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
            f"Based on the current context block, select ADDITIONAL tools needed to answer the question.\n"
            f"Be inclusive in tool selection. Include tools with possible relevance, not just certain relevance."
            f" For example in variant questions be thorough and unless specified use all variant models available.\n"
            f"IMPORTANT: if no additional tools are needed, return [] with no explanation.\n"
            f"Overall, NO explanation, downstream parser can only work with a list of tools.\n\n"
            f"Current context block:\n{current_context_block}\n\n No other words should be included in the response."
            f"{system_msg_base}"
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
    #Always call extract_entities() to extract genes from AA seq's and etc.
    selected_tools.insert(0, "extract_entities")

    selected_tools_main = [t for t in selected_tools if t in TOOL_FN_MAP and t not in used_tools]
    # special case for 'twosample_mr_agent' which allows more complicated semi match, as fomatted twosample_mr_agent-value1-value2
    found_2sample = [
        t for t in selected_tools if t.startswith("twosample_mr_agent") 
    ]
    if found_2sample:
        selected_tools_main.append(found_2sample)
    if DEBUG:
        print(f"[TOOL SELECTION] Selected tools: {selected_tools}")

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
