"""
Evaluate entity extraction against generated question components.

Uses the same nodes and State objects as the main Alvessa pipeline to score:
  - simple_genes.csv (single gene)
  - multi_genes.csv (5–6 genes)
  - variants.csv (5–6 rsIDs)
  - drugs.csv (5–6 drug names)
  - mirnas.csv (miRNA mentions)
from benchmarks_generation/questions/components/entity_recognition.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Sequence, Set, Tuple

import pandas as pd

from src.alvessa.agents.entity_extraction import (
    entity_extraction_node,
    claude_entity_extraction_node,
    flair_entity_extraction_node,
    gliner_entity_extraction_node,
    claude_flair_entity_extraction_node,
    gliner_flair_entity_extraction_node,
    gliner_claude_entity_extraction_node,
)
from src.state import State


REPO_ROOT = Path(__file__).resolve().parents[1]
QUESTIONS_DIR = REPO_ROOT / "benchmarks_generation" / "questions" / "components" / "entity_recognition"


# -------------------------
# Normalization helpers
# -------------------------

def _norm_gene(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "", (value or "").upper())


def _norm_variant(value: str) -> str:
    return (value or "").strip().lower()


def _norm_drug(value: str) -> str:
    return (value or "").strip().lower()


def _split_entities(raw: str) -> List[str]:
    if isinstance(raw, list):
        return [str(x).strip() for x in raw if str(x).strip()]

    text = str(raw).strip()
    if not text:
        return []

    # Try JSON array first
    try:
        parsed = json.loads(text)
        if isinstance(parsed, list):
            return [str(x).strip() for x in parsed if str(x).strip()]
    except Exception:
        pass

    if "|" in text:
        return [part.strip() for part in text.split("|") if part.strip()]

    return [part.strip() for part in text.split(",") if part.strip()]


# -------------------------
# Extraction helpers
# -------------------------

def _extract_predictions(result: Dict, set_type: str) -> Set[str]:
    if set_type in {"simple_genes", "multi_genes"}:
        return {_norm_gene(g) for g in result.get("genes", [])}

    if set_type == "variants":
        variants = set(result.get("variants", {}).keys()) | set(result.get("chr_pos_variants", {}).keys())
        return {_norm_variant(v) for v in variants}

    if set_type == "drugs":
        return {_norm_drug(d) for d in result.get("drugs", [])}

    if set_type == "mirnas":
        # miRNAs are carried in genes; normalize similarly to genes
        return {_norm_gene(g) for g in result.get("genes", [])}

    return set()


def _normalize_expected(values: Sequence[str], set_type: str) -> Set[str]:
    if set_type in {"simple_genes", "multi_genes", "mirnas"}:
        return {_norm_gene(v) for v in values}
    if set_type == "variants":
        return {_norm_variant(v) for v in values}
    if set_type == "drugs":
        return {_norm_drug(v) for v in values}
    return set()


# -------------------------
# Metrics
# -------------------------

def _evaluate_row(query: str, expected: Set[str], set_type: str, node_fn: Callable[[State], Dict]) -> Dict:
    state = State({"messages": [{"role": "user", "content": query}]})
    result = node_fn(state)
    predicted = _extract_predictions(result, set_type)

    tp = len(expected & predicted)
    fp = len(predicted - expected)
    fn = len(expected - predicted)

    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

    return {
        "query": query,
        "expected": sorted(expected),
        "predicted": sorted(predicted),
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "perfect_match": expected == predicted,
    }


def _aggregate(rows: List[Dict]) -> Dict:
    if not rows:
        return {k: 0 for k in ["macro_precision", "macro_recall", "macro_f1", "micro_precision", "micro_recall", "micro_f1", "perfect_matches", "total", "error_rate"]}

    macro_precision = sum(r["precision"] for r in rows) / len(rows)
    macro_recall = sum(r["recall"] for r in rows) / len(rows)
    macro_f1 = sum(r["f1"] for r in rows) / len(rows)

    total_tp = sum(r["tp"] for r in rows)
    total_fp = sum(r["fp"] for r in rows)
    total_fn = sum(r["fn"] for r in rows)

    micro_precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) else 0.0
    micro_recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) else 0.0
    micro_f1 = 2 * micro_precision * micro_recall / (micro_precision + micro_recall) if (micro_precision + micro_recall) else 0.0

    perfect_matches = sum(1 for r in rows if r["perfect_match"])

    return {
        "macro_precision": macro_precision,
        "macro_recall": macro_recall,
        "macro_f1": macro_f1,
        "micro_precision": micro_precision,
        "micro_recall": micro_recall,
        "micro_f1": micro_f1,
        "perfect_matches": perfect_matches,
        "total": len(rows),
        "error_rate": 0.0,
    }


# -------------------------
# Data loading
# -------------------------

def _load_dataset(path: Path, set_type: str, limit: int | None) -> List[Tuple[str, Set[str]]]:
    df = pd.read_csv(path)
    rows = []
    for _, row in df.iterrows():
        query = str(row["query"])
        expected_raw = _split_entities(row.get("recognized_entities", ""))
        expected = _normalize_expected(expected_raw, set_type)
        rows.append((query, expected))
        if limit and len(rows) >= limit:
            break
    return rows


# -------------------------
# Runner
# -------------------------

def run_for_model(model_name: str, node_fn: Callable[[State], Dict], datasets: Dict[str, List[Tuple[str, Set[str]]]]) -> Dict:
    per_set_results = {}
    for set_type, rows in datasets.items():
        eval_rows: List[Dict] = []
        for query, expected in rows:
            eval_rows.append(_evaluate_row(query, expected, set_type, node_fn))
        per_set_results[set_type] = {
            "metrics": _aggregate(eval_rows),
            "rows": eval_rows,
        }
    return per_set_results


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate entity extraction nodes on generated question components.")
    parser.add_argument("--limit", type=int, default=None, help="Limit rows per dataset.")
    parser.add_argument("--output", type=Path, default=REPO_ROOT / "evals" / "results" / "entity_extraction_results.json", help="Path to save JSON results.")
    args = parser.parse_args()

    datasets = {
        "simple_genes": _load_dataset(QUESTIONS_DIR / "simple_genes.csv", "simple_genes", args.limit),
        "multi_genes": _load_dataset(QUESTIONS_DIR / "multi_genes.csv", "multi_genes", args.limit),
        "variants": _load_dataset(QUESTIONS_DIR / "variants.csv", "variants", args.limit),
        "drugs": _load_dataset(QUESTIONS_DIR / "drugs.csv", "drugs", args.limit),
        "mirnas": _load_dataset(QUESTIONS_DIR / "mirnas.csv", "mirnas", args.limit),
    }

    models: Dict[str, Callable[[State], Dict]] = {
        "merged": entity_extraction_node,
        "claude_only": claude_entity_extraction_node,
        "flair_only": flair_entity_extraction_node,
        "gliner_only": gliner_entity_extraction_node,
        "claude_flair": claude_flair_entity_extraction_node,
        "gliner_flair": gliner_flair_entity_extraction_node,
        "gliner_claude": gliner_claude_entity_extraction_node,
    }

    all_results: Dict[str, Dict] = {}
    for name, fn in models.items():
        print(f"Evaluating {name}...")
        all_results[name] = run_for_model(name, fn, datasets)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    serializable = {
        model: {
            set_type: {
                "metrics": metrics_rows["metrics"],
                # Keep only a lightweight per-row view for readability
                "rows": [
                    {
                        "query": r["query"],
                        "expected": r["expected"],
                        "predicted": r["predicted"],
                        "tp": r["tp"],
                        "fp": r["fp"],
                        "fn": r["fn"],
                        "precision": r["precision"],
                        "recall": r["recall"],
                        "f1": r["f1"],
                        "perfect_match": r["perfect_match"],
                    }
                    for r in metrics_rows["rows"]
                ],
            }
            for set_type, metrics_rows in model_results.items()
        }
        for model, model_results in all_results.items()
    }

    with args.output.open("w") as f:
        json.dump(serializable, f, indent=2)

    print("\nSummary (micro F1 per set):")
    for model_name, model_results in all_results.items():
        per_set = {s: v["metrics"]["micro_f1"] for s, v in model_results.items()}
        print(f"{model_name}: {per_set}")


if __name__ == "__main__":
    main()
