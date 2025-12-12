"""
Summarize entity extraction evaluation results (merged model only).

Outputs:
  1) Per-question CSV: rows per question with recall/precision.
  2) Summary CSV: mean and std of recall/precision grouped by entity set (simple_genes, multi_genes, variants, drugs, mirnas).
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from pathlib import Path
from typing import Dict, Any, List

import matplotlib.pyplot as plt
import pandas as pd

ORANGE = "#D95F02"
DEFAULT_INPUT = Path("evals/results/entity_extraction_results.json")
DEFAULT_PER_QUESTION = Path("evals/results/entity_extraction_per_question.csv")
DEFAULT_SUMMARY = Path("evals/results/entity_extraction_summary.csv")


def load_results(path: Path) -> Dict[str, Any]:
    with path.open() as f:
        return json.load(f)


def compute_precision(rec: Dict[str, Any]) -> float:
    tp = rec.get("tp", 0)
    fp = rec.get("fp", 0)
    denom = tp + fp
    return tp / denom if denom else 0.0


def write_per_question(data: Dict[str, Any], output_path: Path, model_name: str = "merged") -> None:
    rows: List[Dict[str, Any]] = []
    model_data = data.get(model_name, {})
    for set_name, payload in model_data.items():
        for row in payload.get("rows", []):
            recall = row.get("recall", 0)
            precision = compute_precision(row)
            rows.append(
                {
                    "model": model_name,
                    "set": set_name,
                    "query": row.get("query", ""),
                    "recall": recall,
                    "precision": precision,
                    "tp": row.get("tp", 0),
                    "fp": row.get("fp", 0),
                    "fn": row.get("fn", 0),
                }
            )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["model", "set", "query", "recall", "precision", "tp", "fp", "fn"]
    with output_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_summary(data: Dict[str, Any], output_path: Path, model_name: str = "merged") -> None:
    rows: List[Dict[str, Any]] = []
    model_data = data.get(model_name, {})
    for set_name, payload in model_data.items():
        recalls = []
        precisions = []
        for row in payload.get("rows", []):
            recalls.append(row.get("recall", 0))
            precisions.append(compute_precision(row))

        recall_mean = statistics.mean(recalls) if recalls else 0.0
        recall_std = statistics.pstdev(recalls) if len(recalls) > 1 else 0.0
        precision_mean = statistics.mean(precisions) if precisions else 0.0
        precision_std = statistics.pstdev(precisions) if len(precisions) > 1 else 0.0

        rows.append(
            {
                "model": model_name,
                "set": set_name,
                "recall_mean": recall_mean,
                "recall_std": recall_std,
                "precision_mean": precision_mean,
                "precision_std": precision_std,
                "recall_err": recall_std,
                "precision_err": precision_std,
                "n": len(recalls),
            }
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "model",
        "set",
        "recall_mean",
        "recall_std",
        "recall_err",
        "precision_mean",
        "precision_std",
        "precision_err",
        "n",
    ]
    with output_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _format_set_name(name: str) -> str:
    if name == "simple_genes":
        return "One gene"
    pretty = name.replace("_", " ").title()
    return pretty.replace("Mirnas", "miRNAs")


def _plot_metric(df: pd.DataFrame, metric: str, err: str, ylabel: str, title: str, out_path: Path) -> None:
    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 14,
            "axes.labelsize": 13,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
        }
    )

    fig, ax = plt.subplots(figsize=(8, 4), dpi=200)
    labels = [_format_set_name(s) for s in df["set"]]
    values = df[metric]
    errors = df[err].fillna(0)

    bars = ax.bar(labels, values, yerr=errors, color=ORANGE, edgecolor="none", linewidth=0, capsize=5)
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", linestyle="--", linewidth=0.7, alpha=0.5, color="#444444")
    plt.xticks(rotation=0, ha="center")

    # Match box styling: hide top/right, keep left/bottom subtle
    for spine_name, spine in ax.spines.items():
        if spine_name in ("top", "right"):
            spine.set_visible(False)
        else:
            spine.set_visible(True)
            spine.set_color("#444444")
            spine.set_linewidth(0.8)
    ax.set_facecolor("white")
    ax.figure.set_facecolor("white")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Report entity extraction evaluation results (merged model only).")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Path to JSON results file.")
    parser.add_argument("--per-question", type=Path, default=DEFAULT_PER_QUESTION, help="Path to write per-question CSV.")
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY, help="Path to write summary CSV.")
    parser.add_argument("--plot", action="store_true", help="Also save bar plots for recall and precision by set (merged).")
    args = parser.parse_args()

    data = load_results(args.input)
    write_per_question(data, args.per_question, model_name="merged")
    write_summary(data, args.summary, model_name="merged")
    print(f"Wrote per-question to {args.per_question}")
    print(f"Wrote summary to {args.summary}")

    if args.plot:
        df = pd.read_csv(args.summary)

        recall_path = args.summary.with_name(args.summary.stem + "_recall.png")
        _plot_metric(
            df,
            metric="recall_mean",
            err="recall_err",
            ylabel="Recall (mean ± std)",
            title="Merged model: Recall by set",
            out_path=recall_path,
        )

        precision_path = args.summary.with_name(args.summary.stem + "_precision.png")
        _plot_metric(
            df,
            metric="precision_mean",
            err="precision_err",
            ylabel="Precision (mean ± std)",
            title="Merged model: Precision by set",
            out_path=precision_path,
        )

        print(f"Wrote plots to {recall_path} and {precision_path}")


if __name__ == "__main__":
    main()
