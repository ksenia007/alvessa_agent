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
from tool_descriptions import TOOL_CATALOG, TOOL_FN_MAP
from config import TOOL_SELECTOR_MODEL
from claude_client import claude_call
from typing import Any, Dict, List
import inspect
import asyncio

def select_tools_and_run_ALL(state: State) -> State:
    """Dynamic tool runner without LLM — just calls all tools in order. 
    Primarily used for debugging, but also could be viewed as data collection
    """
    # You already have genes from extract_genes
    state = humanbase_predictions_agent(state)
    state = uniprot_node(state)  # for base genes
    state = gwas_associations_agent(state)
    state = uniprot_node(state)  # again for gwas genes (you may want to scope this)
    state = bioGRID_predictions_agent(state)

    return state


def format_state_for_prompt(state: State) -> str:
    
    genes = state.get("genes", [])
    genes_str = ", ".join(genes) if genes else "None"
    question = state["messages"][-1]["content"]

    catalog = "\n".join(f"- {name}: {desc}" for name, desc in TOOL_CATALOG.items())
    
    return f"""You are an assistant deciding which tools to use to answer a biomedical question. User question: \"\"\"{question}\"\"\" Extracted gene symbols: {genes_str} Available tools: {catalog}.\n Which tools should be called, and in what order? Respond ONLY with a Python list of tool names. Example: ["query_by_trait", "humanbase_functions", "uniprot_base"] or ["humanbase_functions", "uniprot_base", "gwas"] or ["query_by_trait", "gwas", "BioGRID"]"""

def select_tools_and_run_dynamic(state: State) -> State:
    
    updates: State = dict(state)  # start with all existing keys
    updates.setdefault("used_tools", [])

    # prepare the prompt for tool selection
    prompt = format_state_for_prompt(state)
    
    try:
        completion = claude_call(
            model=TOOL_SELECTOR_MODEL,
            max_tokens=256,
            temperature=0,
            messages=[
                {"role": "user", "content": prompt}
            ]
        )
        tool_response = completion.content[0].text.strip()
        
    # if there is an error in the LLM call, we catch it and return an empty list
    except Exception as e:
        print(f"[CLAUDE ERROR] {e}")
        tool_response = "[]"

    try:
        selected_tools: List[str] = eval(tool_response)
        assert isinstance(selected_tools, list)
    except Exception as e:
        print(f"[TOOL SELECTION ERROR] Could not parse: {tool_response}\n{e}")
        selected_tools = []

    print(f"[TOOL SELECTION] Claude (Haiku) selected: {selected_tools}")

    state.update({'used_tools': selected_tools})  # record used tools
    
    # update state
    for name in selected_tools:
        fn = TOOL_FN_MAP.get(name)
        if not fn:
            continue

        print(f"[TOOL RUN] → {name}", flush=True)
        out = fn(state)
        if not isinstance(out, dict):
            continue

        out.pop("messages", None)

        # Keep existing genes if tool returns nothing useful
        if "genes" in out:
            # Add genes to the current genes list
            state["genes"] = list(set(state["genes"] + out["genes"]))
            if not isinstance(out["genes"], list) or not out["genes"]:
                out.pop("genes", None)

        # Merge tool output safely
        state.update(out)

    return state


def run_async_sync(fn):
    def wrapper(*args, **kwargs):
        async def run_coroutine():
            return await fn(*args, **kwargs)
        return asyncio.run(run_coroutine())
    return wrapper