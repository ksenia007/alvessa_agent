from functools import partial
from state import State
from typing import Callable, Any, Optional
import asyncio


def make_tool_wrapper(tool_fn: Callable[[Any], Any], *, state_key: Optional[str] = None):
    return 
    # """
    # Wrap a LangGraph-style (state→state) function so the wrapper:
    # 1. Accepts **text** (what the LLM agent passes in)
    # 2. Builds a minimal State (msg + optional genes)
    # 3. Runs the real tool (sync or async) **to completion**
    # 4. Returns **{state_key: value}**  (or the whole state)  –> _dict_, never str
    # """
    # async def _run(text: str) -> dict:
    #     from state import State

    #     st = State()
    #     st["messages"] = [{"role": "user", "content": text}]

    #     # call the real tool (it may be sync or async)
    #     out = tool_fn(st)
    #     if asyncio.iscoroutine(out):
    #         out = await out

    #     return {state_key: out[state_key]} if state_key else out