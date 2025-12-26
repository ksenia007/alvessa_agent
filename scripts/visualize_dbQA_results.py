from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import pandas as pd


GROUP_RULES = [
    ("DisGeNet+\nOMIM", lambda q: re.search(r"disgenet", q) and re.search(r"omim", q)),
    ("Ensembl", lambda q: "ensembl" in q),
    ("miRDB", lambda q: "mirdb" in q),
    ("MouseMine", lambda q: "mousemine" in q),
    ("GTRD", lambda q: "gene transcription regulation database" in q),
    ("ClinVar", lambda q: "clinvar" in q),
    ("P-HIPSter", lambda q: "p-hipster" in q),
    ("MSigDB", lambda q: re.search(r"c\d+\s+collection", q) or re.search(r"c\d+\s+subcollection", q)),
]


def _assign_group(question: str) -> str:
    q = (question or "").lower()
    for name, fn in GROUP_RULES:
        try:
            if fn(q):
                return name
        except Exception:
            continue
    return "other"


def _load_data(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df["is_correct"] = pd.to_numeric(df.get("is_correct", 0), errors="coerce").fillna(0).astype(int)
    df["question"] = df.get("question", "").fillna("").astype(str)
    df["question_group"] = df["question"].apply(_assign_group)
    return df


def _compute_accuracy(df: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        df.groupby("question_group")["is_correct"]
        .agg(["mean", "count"])
        .rename(columns={"mean": "accuracy", "count": "n"})
        .reset_index()
    )
    grouped = grouped.sort_values("accuracy", ascending=True)
    return grouped


def _style_axes(ax, *, theme: str) -> None:
    if theme == "white":
        ax.set_facecolor("white")
        ax.figure.set_facecolor("white")
        color = "#111111"
        spine_color = "#444444"
    else:  # black theme
        ax.set_facecolor("none")
        ax.figure.set_facecolor("none")
        color = "#f5f5f5"
        spine_color = "#888888"

    ax.tick_params(axis="both", labelsize=12, colors=color)
    ax.set_ylabel("Accuracy", fontsize=14, color=color)
    ax.set_xlabel("")
    ax.set_ylim(0, 1)
    for spine in ax.spines.values():
        spine.set_color(spine_color)
    ax.yaxis.label.set_color(color)
    ax.xaxis.label.set_color(color)
    ax.tick_params(colors=color)
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.5, color=spine_color)


def _plot_groups(df: pd.DataFrame, *, out_dir: Path, dpi: int = 300) -> Tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)

    labels = df["question_group"].tolist()
    accuracies = df["accuracy"].tolist()
    counts = df["n"].tolist()

    def _make_fig(theme: str, fname: str, bar_color: str) -> Path:
        fig, ax = plt.subplots(figsize=(max(9, len(labels) * 0.7), 4.5))
        ax.bar(labels, accuracies, color=bar_color, edgecolor="none")
        ax.set_xticklabels(labels, rotation=0, ha="center", color=("#111111" if theme == "white" else "#f5f5f5"))
        _style_axes(ax, theme=theme)
        outfile = out_dir / fname
        if theme == "white":
            fig.savefig(outfile, dpi=dpi, bbox_inches="tight", transparent=False, facecolor="white")
        else:
            fig.savefig(outfile, dpi=dpi, bbox_inches="tight", transparent=True)
        plt.close(fig)
        return outfile

    white_path = _make_fig("white", "dbqa_accuracy_by_group_white.png", bar_color="#D95F02")
    black_path = _make_fig("black", "dbqa_accuracy_by_group_black.png", bar_color="#FFB570")
    return white_path, black_path


def _save_group_assignments(df: pd.DataFrame, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "question_group_assignments.csv"
    df[["question", "question_group", "source_name", "source_folder"]].to_csv(out_path, index=False)
    return out_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Visualize dbQA accuracy by inferred question group.")
    parser.add_argument("file_path", help="Path to benchmark_summary.csv")
    parser.add_argument("--exclude-disgenet", action="store_true", help="Also plot a version excluding DisGeNet questions.")
    args = parser.parse_args()

    csv_path = Path(args.file_path).expanduser().resolve()
    if not csv_path.exists():
        print(f"Benchmark CSV not found: {csv_path}", file=sys.stderr)
        return 1

    try:
        df = _load_data(csv_path)
    except Exception as exc:
        print(f"Failed to load CSV: {exc}", file=sys.stderr)
        return 1

    if df.empty:
        print("No data to plot (empty CSV).", file=sys.stderr)
        return 1

    def plot_and_save(data: pd.DataFrame, suffix: str = ""):
        overall_acc = data["is_correct"].mean()
        overall_row = pd.DataFrame({"question_group": ["All"], "accuracy": [overall_acc], "n": [len(data)]})
        group_acc = _compute_accuracy(data)
        full = pd.concat([overall_row, group_acc], ignore_index=True)
        figures_dir = csv_path.parent / "figures"
        white_path, black_path = _plot_groups(full, out_dir=figures_dir, dpi=300 if not suffix else 300)
        assign_path = _save_group_assignments(data, figures_dir)
        print(f"Saved white plot{suffix} to: {white_path}")
        print(f"Saved black plot{suffix} to: {black_path}")
        print(f"Saved question->group assignments{suffix} to: {assign_path}")

    # Plot full data
    plot_and_save(df, suffix="")

    # Optionally plot excluding DisGeNet
    if args.exclude_disgenet:
        df_no_disgenet = df[df["question_group"] != "DisGeNet+\nOMIM"].copy()
        if df_no_disgenet.empty:
            print("No data left after excluding DisGeNet; skipping filtered plots.", file=sys.stderr)
        else:
            plot_and_save(df_no_disgenet, suffix=" (no DisGeNet)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
