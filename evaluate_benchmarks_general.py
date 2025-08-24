"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-08-21
Updated: 2025-08-21


Description: 

Main code to run MC evaluation; supports Alvessa and Claude. To addd more baselines reference `method_baseline`."""

from __future__ import annotations
import json
import pandas as pd
from typing import Dict, Any
import os
import re
import time
from claude_client import claude_call
from config import CONDITIONED_MODEL
import os, time, tempfile, shutil, re
from typing import Callable, Dict, Any, Optional, Tuple
import pandas as pd

from run import run_pipeline



def normalize_answer(x: Any) -> str:
    """Lowercase, strip; if text contains A-D inside backticks or JSON, try to extract."""
    if x is None:
        return ""
    s = str(x).strip()
    return s

def compute_is_correct(correct: Any, model_answer: Any) -> bool:
    return normalize_answer(correct).lower() == normalize_answer(model_answer).lower()

def atomic_to_csv(df: pd.DataFrame, path: str) -> None:
    """Write CSV atomically to avoid partial writes."""
    tmp_dir = os.path.dirname(os.path.abspath(path)) or "."
    fd, tmp = tempfile.mkstemp(prefix=".tmp_", dir=tmp_dir, suffix=".csv")
    os.close(fd)
    try:
        df.to_csv(tmp, index=False)
        shutil.move(tmp, path)
    finally:
        try:
            if os.path.exists(tmp):
                os.remove(tmp)
        except Exception:
            pass

def load_existing(file_save: Optional[str]) -> Tuple[pd.DataFrame, set[Tuple[str,str]]]:
    """Load existing results and return df + set of (question, method) already done."""
    if not file_save or not os.path.exists(file_save):
        return pd.DataFrame(), set()
    try:
        existing = pd.read_csv(file_save)
    except Exception as e:
        print(f"[resume] Warning: failed to load existing results at {file_save}: {e}")
        return pd.DataFrame(), set()

    # Ensure required columns exist
    if "question" not in existing.columns:
        existing["question"] = ""
    # Normalize method column (fallback to model when method missing)
    if "method" not in existing.columns:
        if "model" in existing.columns:
            existing["method"] = existing["model"].astype(str)
        else:
            existing["method"] = "unknown"

    done = set(existing["question"].astype(str)) # Set of unique questions
    print(f"[resume] Loaded {len(existing)} rows, {len(done)} unique (question, method) pairs.")
    return existing, done

def merge_and_dedupe(existing: pd.DataFrame, new: pd.DataFrame) -> pd.DataFrame:
    """Merge on rows and drop duplicate (question, method) keeping the first (existing wins)."""
    if existing is None or existing.empty:
        combined = new.copy()
    else:
        combined = pd.concat([existing, new], ignore_index=True)

    # Guarantee minimal columns
    for col in ["question", "correct_answer", "model_answer", "is_correct", "method", "model", "used_models"]:
        if col not in combined.columns:
            combined[col] = None

    # Deduplicate on (question, method)
    combined = combined.drop_duplicates(subset=["question", "method"], keep="first")

    # Recompute is_correct where missing
    if "is_correct" in combined.columns:
        mask = combined["is_correct"].isna()
        if mask.any():
            combined.loc[mask, "is_correct"] = combined.loc[mask].apply(
                lambda x: compute_is_correct(x.get("correct_answer", ""), x.get("model_answer", "")), axis=1
            )
    else:
        combined["is_correct"] = combined.apply(
            lambda x: compute_is_correct(x.get("correct_answer", ""), x.get("model_answer", "")), axis=1
        )
    return combined


def method_alvessa(question: str, system_msg: str) -> Dict[str, Any]:
    """
    Returns a unified record dict for the question using Alvessa (run_pipeline).
    """
    system_msg += "\n Try to use the following information (if available) to answer the question (while still only replying with A / B / C / D - \n"
    result = run_pipeline(question, prompt=system_msg, mc_setup=True)
    answer = result.get("llm_json", {}).get("answer", "")
    used = result.get("used_tools", None)

    return {
        "model_answer": answer,
        "model": "Alvessa",
        "used_models": used,
        "method": "alvessa",
    }

def method_baseline(question: str, system_msg: str) -> Dict[str, Any]:
    """
    Returns a unified record dict for the question using baseline claude_call.
    """
    print('Running baseline method with Claude with system message:', system_msg)
    raw = claude_call(
        model=CONDITIONED_MODEL,
        temperature=0,
        max_tokens=10,
        system=system_msg,
        messages=[{"role": "user", "content": f"User asked: {question}"}],
    )
    # Handle SDK shapes defensively
    first = raw.content[0]
    text = getattr(first, "text", first)
    llm_resp = str(text).strip()

    return {
        "model_answer": llm_resp,
        "model": CONDITIONED_MODEL,
        "used_models": None,
        "method": "Claude",
    }



METHODS: dict[str, Callable[[str, str], Dict[str, Any]]] = {
    "alvessa": method_alvessa,
    "claude": method_baseline,
}
def run_benchmark(
    file_path: str,
    system_msg: str,
    *,
    method: str = "alvessa",
    max_rows: int = -1,
    file_save: Optional[str] = None,
    throttle_s: float = 0,
    method_fn: Optional[Callable[[str, str], Dict[str, Any]]] = None,
) -> pd.DataFrame:
    """
    Unified runner supporting:
      - Resume/skip previously answered (question, method)
      - Incremental atomic saves
      - Pluggable method: 'alvessa', 'baseline', or custom method_fn(question, system_msg) -> dict
    Output columns: question, correct_answer, model_answer, is_correct, method, model, used_models
    """
    df = pd.read_csv(file_path)

    if max_rows > 0:
        df = df.sample(frac=1, random_state=42).reset_index(drop=True).head(max_rows)

    # Select method implementation
    if method_fn is None:
        if method not in METHODS:
            raise ValueError(f"Unknown method '{method}'. Available: {list(METHODS)}")
        method_impl = METHODS[method]
    else:
        method_impl = method_fn

    existing_df, done_questions = load_existing(file_save)
    results_rows = []

    for i, row in df.iterrows():
        q = str(row.get("question", ""))
        if not q:
            print("[skip] empty question row")
            continue

        # Skip if already done for this method
        if q in done_questions:
            print(f"[skip] Already answered: ({q[:100]!r}, method={method})")
            continue

        print("\n" + "=" * 80)
        print(f"[{method}] Q:", q)

        try:
            rec = method_impl(q, system_msg)
            model_answer = rec.get("model_answer", "")
        except Exception as e:
            print(f"[error] {method} failed: {e}")
            rec = {"model_answer": ""}
            model_answer = ""

        # Build normalized row
        correct = row.get("answer", "")
        out = {
            "question": q,
            "correct_answer": correct,
            "model_answer": model_answer,
            "is_correct": compute_is_correct(correct, model_answer),
            "method": rec.get("method", method),
            "model": rec.get("model", method),
            "used_models": rec.get("used_models", None),
        }
        results_rows.append(out)

        print(f"Model answer: {model_answer}")
        print(f"Correct answer: {correct}")

        # Incremental save
        if file_save:
            # check that the directory exists
            os.makedirs(os.path.dirname(file_save), exist_ok=True)
            
            new_df_so_far = pd.DataFrame(results_rows)
            combined = merge_and_dedupe(existing_df, new_df_so_far)
            atomic_to_csv(combined, file_save)

        if throttle_s and throttle_s > 0:
            time.sleep(throttle_s)

    new_df = pd.DataFrame(results_rows)
    final_df = merge_and_dedupe(existing_df, new_df)

    if not final_df.empty and "is_correct" in final_df.columns:
        total_correct = int(final_df["is_correct"].sum())
        print(f"Correct answers (all in file): {total_correct} / {len(final_df)}")
    else:
        print("No results to summarize.")

    return final_df




