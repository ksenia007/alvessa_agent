from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import pandas as pd


def _load_data(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    required_cols = {"source_folder", "is_correct"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {csv_path}: {', '.join(sorted(missing))}")

    # Coerce correctness to numeric 0/1
    df["is_correct"] = pd.to_numeric(df["is_correct"], errors="coerce").fillna(0).astype(int)
    df["source_folder"] = df["source_folder"].fillna("UNKNOWN").astype(str)
    return df


def _compute_accuracy(df: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        df.groupby("source_folder")["is_correct"]
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


def _plot_accuracy(df: pd.DataFrame, *, out_dir: Path, dpi: int = 300) -> Tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)

    labels = df["source_folder"].tolist()
    accuracies = df["accuracy"].tolist()
    counts = df["n"].tolist()

    def _make_fig(theme: str, fname: str, bar_color: str) -> Path:
        fig, ax = plt.subplots(figsize=(max(6, len(labels) * 0.5), 4))
        ax.bar(labels, accuracies, color=bar_color, edgecolor="none")
        for x, acc, n in zip(range(len(labels)), accuracies, counts):
            ax.text(x, acc + 0.02, f"{acc:.2f}\n(n={n})", ha="center", va="bottom", fontsize=10, color=("#111111" if theme == "white" else "#f5f5f5"))
        _style_axes(ax, theme=theme)
        outfile = out_dir / fname
        if theme == "white":
            fig.savefig(outfile, dpi=dpi, bbox_inches="tight", transparent=False, facecolor="white")
        else:
            fig.savefig(outfile, dpi=dpi, bbox_inches="tight", transparent=True)
        plt.close(fig)
        return outfile

    white_path = _make_fig("white", "accuracy_by_source_white.png", bar_color="#D95F02")
    black_path = _make_fig("black", "accuracy_by_source_black.png", bar_color="#FFB570")
    return white_path, black_path


def _compute_accuracy_by_folder_name(df: pd.DataFrame, folder_order: list[str] | None = None) -> pd.DataFrame:
    grouped = (
        df.groupby(["source_folder", "source_name"])["is_correct"]
        .agg(["mean", "count"])
        .rename(columns={"mean": "accuracy", "count": "n"})
        .reset_index()
    )
    if folder_order:
        grouped["_folder_rank"] = grouped["source_folder"].apply(lambda x: folder_order.index(x) if x in folder_order else len(folder_order))
        grouped = grouped.sort_values(["_folder_rank", "source_name"])
        grouped = grouped.drop(columns=["_folder_rank"])
    return grouped


FOLDER_NAME_MAP = {
    "BioGRID": "BioGRID",
    "chembl": "ChEMBL",
    "gencode": "GENCODE",
    "GWAS": "GWAS",
    "gwas+alphamissense": "GWAS+\nAlphaMissense",
    "miRDB": "miRDB",
    "MSigDB": "MSigDB",
    "OMIM": "OMIM",
    "OpenTargets": "OpenTargets",
    "prot": "Protein\nstructure",
    "reactome": "Reactome",
}


def _format_folder(name: str) -> str:
    return FOLDER_NAME_MAP.get(name, name)


def _plot_grouped(df: pd.DataFrame, *, out_dir: Path, folder_order: list[str], dpi: int = 300) -> Tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)

    positions = []
    heights = []
    labels = []
    folder_bounds = []
    x = 0.0
    gap_between_sets = 0.35
    gap_between_folders = gap_between_sets * 3.0

    for f in folder_order:
        sub = df[df["source_folder"] == f].copy()
        sub = sub.sort_values("source_name")  # preserve natural order (set1, set2,...)
        start = x
        for _, row in sub.iterrows():
            positions.append(x)
            heights.append(row["accuracy"])
            labels.append(str(row["source_name"]).removesuffix(".csv"))
            x += gap_between_sets
        end = x - gap_between_sets if sub.shape[0] else x
        folder_bounds.append((f, start, end))
        x += gap_between_folders

    def _make_fig(theme: str, fname: str, bar_color: str) -> Path:
        fig, ax = plt.subplots(figsize=(max(13, x * 0.7), 5.0))
        ax.bar(positions, heights, color=bar_color, width=0.32, edgecolor="none")

        # Set labels below axis
        for xpos, lbl in zip(positions, labels):
            ax.text(xpos, -0.04, lbl, ha="center", va="top", fontsize=9, color=("#111111" if theme == "white" else "#f5f5f5"))

        # Folder labels as bracket-like connectors under set labels
        for name, start, end in folder_bounds:
            mid = (start + end) / 2 if end >= start else start
            y_bracket = -0.08
            color = "#888888" if theme == "white" else "#aaaaaa"
            ax.hlines(y_bracket, start, end, colors=color, linewidth=1.2)
            ax.vlines([start, end], y_bracket, y_bracket - 0.015, colors=color, linewidth=1.2)
            ax.text(mid, y_bracket - 0.035, _format_folder(name), ha="center", va="top", fontsize=10.5, color=color)

        ax.set_xticks([])
        ax.set_xlim(min(positions) - 0.6, max(positions) + 0.6)
        ax.set_ylim(0, 1.08)
        _style_axes(ax, theme=theme)
        ax.set_xlabel("")
        ax.set_ylabel("Accuracy", fontsize=14, color=("#111111" if theme == "white" else "#f5f5f5"))

        outfile = out_dir / fname
        if theme == "white":
            fig.savefig(outfile, dpi=dpi, bbox_inches="tight", transparent=False, facecolor="white")
        else:
            fig.savefig(outfile, dpi=dpi, bbox_inches="tight", transparent=True)
        plt.close(fig)
        return outfile

    white_path = _make_fig("white", "accuracy_by_source_and_set_white.png", bar_color="#D95F02")
    black_path = _make_fig("black", "accuracy_by_source_and_set_black.png", bar_color="#FFB570")
    return white_path, black_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Visualize benchmark accuracy by source folder.")
    parser.add_argument("file_path", help="Path to benchmark_summary.csv (from benchmark or benchmark_all).")
    args = parser.parse_args()

    csv_path = Path(args.file_path).expanduser().resolve()
    if not csv_path.exists():
        print(f"Benchmark CSV not found: {csv_path}", file=sys.stderr)
        return 1

    try:
        df = _load_data(csv_path)
        acc_folder = _compute_accuracy(df)
        acc_folder_set = _compute_accuracy_by_folder_name(df, acc_folder["source_folder"].tolist())
    except Exception as exc:
        print(f"Failed to load/compute accuracy: {exc}", file=sys.stderr)
        return 1

    if acc_folder.empty:
        print("No data to plot (empty accuracy table).", file=sys.stderr)
        return 1

    figures_dir = csv_path.parent / "figures"
    white_path, black_path = _plot_accuracy(acc_folder, out_dir=figures_dir)
    print(f"Saved white-background-friendly plot to: {white_path}")
    print(f"Saved black-background-friendly plot to: {black_path}")

    if not acc_folder_set.empty:
        white2, black2 = _plot_grouped(acc_folder_set, out_dir=figures_dir, folder_order=acc_folder["source_folder"].tolist())
        print(f"Saved grouped white plot to: {white2}")
        print(f"Saved grouped black plot to: {black2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
