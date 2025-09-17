"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-06-25
Updated: 2025-06-26


Description: 

Thread-safe, rate-limited wrapper around Anthropic’s SDK."""

from __future__ import annotations
import collections
import threading
import time
from typing import Any

from anthropic import Anthropic
from src.config import RATE_LIMIT, ANTHROPIC_API_KEY

# Anthropic client (re-used across threads)
_model = Anthropic(api_key=ANTHROPIC_API_KEY)

# simple token-bucket state
_call_log: "collections.deque[float]" = collections.deque()
_log_lock: threading.Lock = threading.Lock()


def claude_call(**kwargs: Any):
    """
    Rate-limit calls to Anthropic.

    Parameters
    ----------
    **kwargs
        The exact keyword args forwarded to `Anthropic.messages.create`.

    Returns
    -------
    Any
        The raw SDK response object.
    """
    with _log_lock:
        now = time.time()

        # Drop timestamps older than 60 s
        while _call_log and now - _call_log[0] > 60:
            _call_log.popleft()

        # Throttle if limit exceeded
        if len(_call_log) >= RATE_LIMIT:
            sleep_for = 60 - (now - _call_log[0]) + 0.1
            time.sleep(sleep_for)

        _call_log.append(now)

    return _model.messages.create(**kwargs)
