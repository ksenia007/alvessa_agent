# %%
"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-08-21
Updated: 2025-08-21


Description: 

File to compare performance by the LabBench dataset. Uses direct match. Uses DEFAULT_INPUT for quick check. 
Not required for any other scripts."""

import re
import os
import argparse
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DEFAULT_INPUT = "benchmarks_generation/results/labbench/claude/dbqa_mc_results_subset.csv"

SOURCE_PATTERNS = {
    "miRDB": re.compile(r"\bmiRDB\b", re.IGNORECASE),
    "OMIM": re.compile(r"\bOMIM\b", re.IGNORECASE),
    "DisGeNet": re.compile(r"\bDisGeN?et\b", re.IGNORECASE),  
    "Gene Transcription Regulation Database": re.compile(
        r"\bGene\s+Transcription\s+Regulation\s+Database\b|\bGTRD\b", re.IGNORECASE
    ),
    "ClinVar": re.compile(r"\bClinVar\b", re.IGNORECASE),
    "P-HIPSter": re.compile(r"\bP[-\s]?HIPS?TER\b", re.IGNORECASE),
    "C7": re.compile(r"\b(MSigDB\s+)?C7(\s+immunologic\s+signatures)?\b", re.IGNORECASE),
}

def tag_sources(text: str, multi: bool = True):
    if not isinstance(text, str):
        return []
    hits = [name for name, pat in SOURCE_PATTERNS.items() if pat.search(text)]
    if not hits:
        return []
    return hits if multi else [hits[0]]

def ensure_is_correct(df: pd.DataFrame) -> pd.Series:
    """Always return boolean is_correct column."""
    if "is_correct" in df.columns:
        return df["is_correct"].astype(bool)

    ca = df.get("correct_answer", "").astype(str).str.strip().str.lower()
    ma = df.get("model_answer", "").astype(str).str.strip().str.lower()
    return (ca == ma)

def accuracy_table(df_typed: pd.DataFrame) -> pd.DataFrame:
    base = df_typed.copy()
    base["is_correct"] = ensure_is_correct(base)
    exp = base.copy()
    exp["source"] = exp["source_tag"].fillna("other")

    exp["source"] = exp["source"].replace("", "other").fillna("other")

    grp = (exp.groupby("source", dropna=False)
             .agg(n=("is_correct", "size"),
                  accuracy=("is_correct", "mean"))
             .reset_index()
             .sort_values(["accuracy", "n"], ascending=[False, False]))
    return grp

def plot_accuracy_by_source(table: pd.DataFrame, title: str, save_path = None):
    if table.empty:
        print("No data to plot."); return

    # Aesthetics
    plt.rcParams.update({
        "axes.titlesize": 22,
        "axes.labelsize": 18,
        "xtick.labelsize": 18,  
        "ytick.labelsize": 16,
        "legend.fontsize": 16,
    })

    labels = table["source"].tolist()
    x = np.arange(len(labels))
    vals = table["accuracy"].to_numpy()
    ns   = table["n"].to_numpy()

    fig, ax = plt.subplots(figsize=(12, 6))
    bars = ax.bar(x, vals, color="grey")  # Claude = grey

    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Accuracy")
    ax.set_xlabel("")  # drop x-axis label per your request
    ax.set_title(title, pad=14)
    ax.set_xticks(x)
    ax.set_xticklabels([textwrap.fill(s, width=28) for s in labels], rotation=0, ha="center")
    ax.yaxis.grid(True, linestyle="--", alpha=0.35)
    ax.set_axisbelow(True)

    # Annotate bars with accuracy and n
    for rect, acc, n in zip(bars, vals, ns):
        if not np.isnan(acc):
            ax.annotate(f"{acc:.2f}\n(n={int(n)})",
                        (rect.get_x() + rect.get_width()/2, acc),
                        textcoords="offset points", xytext=(0, 6),
                        ha="center", va="bottom", fontsize=13)

    fig.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches="tight")
        print(f"Saved plot to {save_path}")
    else:
        plt.show()

def main():
    ap = argparse.ArgumentParser(description="Tag LabBench questions by source and plot accuracies.")
    ap.add_argument("-i", "--input", default=DEFAULT_INPUT, help="Input CSV")
    ap.add_argument("-o", "--output", default=None, help="Output typed CSV (default: *_typed.csv)")
    ap.add_argument("--question-col", default="question", help="Question text column")
    ap.add_argument("--plot", default=None,
                    help="Where to save the plot (default: next to output, *_accuracy.png)")
    args = ap.parse_args()

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input not found: {args.input}")

    df = pd.read_csv(args.input)
    if args.question_col not in df.columns:
        raise ValueError(f"Column '{args.question_col}' not in CSV. Columns: {list(df.columns)}")

    # Tag
    tags = df[args.question_col].apply(lambda q: tag_sources(q))
    df["source_tag_list"] = tags
    df["source_tag"] = tags.apply(lambda lst: "; ".join(lst) if lst else "other")

    # Table & print
    tbl = accuracy_table(df)
    print("\nAccuracy by source:")
    print(tbl.to_string(index=False, formatters={"accuracy": "{:.3f}".format}))

    plot_title = "Claude Accuracy by db for LabBench"
    plot_accuracy_by_source(tbl, plot_title)

if __name__ == "__main__":
    main()
# %%