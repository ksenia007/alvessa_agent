from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch
import argparse

# Configure benchmark CSVs to compare (model name -> CSV path)
MODEL_FILES = {
    "Alvessa": "/Users/sokolova/Documents/research/alvessa_agent/out/FINAL_DBQA_20251216-125635_cli/benchmark_summary.csv",
    "Claude\nSonnet 4.5": "/Users/sokolova/Documents/research/alvessa_agent/chat_claude_baselines/FINAL_claude_baseline_dbQA_20251215-2308.csv",
    "ChatGPT 5.1": "/Users/sokolova/Documents/research/alvessa_agent/chat_claude_baselines/FINAL_chatgpt_baseline_dbQA_20251215-2330.csv"
}

ALVESSA_COLOR = "#D95F02"
# Palette/hatches for non-Alvessa models (cycled in order of appearance)
OTHER_COLORS = [ "#727272", "#555555","#8C8C8C", "#A6A6A6", "#BEBEBE"]
OTHER_HATCHES = ["..", "xx", "//", "\\\\",  "++", "--"]

WIDTH_PLOT = 3.5

# Same grouping logic as visualize_dbQA_results.py
GROUP_RULES: List[Tuple[str, callable]] = [
    ("DisGeNet+\nOMIM", lambda q: ("disgenet" in q) and ("omim" in q)),
    ("Ensembl", lambda q: "ensembl" in q),
    ("miRDB", lambda q: "mirdb" in q),
    ("MouseMine", lambda q: "mousemine" in q),
    ("GTRD", lambda q: "gene transcription regulation database" in q),
    ("ClinVar", lambda q: "clinvar" in q),
    ("P-HIPSter", lambda q: "p-hipster" in q),
    ("MSigDB", lambda q: (" c" in q and "collection" in q) or (" c" in q and "subcollection" in q)),
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


def _load_benchmark(csv_path: Path, *, drop_disgenet: bool = False) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df = df.drop_duplicates('question')
    df["is_correct"] = pd.to_numeric(df.get("is_correct", 0), errors="coerce").fillna(0).astype(int)
    df["question"] = df.get("question", "").fillna("").astype(str)
    if drop_disgenet:
        mask = ~df["question"].str.lower().str.contains("disgenet", na=False)
        df = df[mask].reset_index(drop=True)
    df["source_folder"] = df.get("source_folder", "").fillna("").astype(str).str.strip()
    df["source_name"] = df.get("source_name", "").fillna("").astype(str).str.strip()
    df["question_group"] = df["question"].apply(_assign_group)
    return df


def _compute_by_group(df: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        df.groupby("question_group")["is_correct"]
        .agg(["mean", "count"])
        .rename(columns={"mean": "accuracy", "count": "n"})
        .reset_index()
    )
    overall = pd.DataFrame({"question_group": ["All"], "accuracy": [df["is_correct"].mean()], "n": [len(df)]})
    full = pd.concat([overall, grouped], ignore_index=True)
    return full


def _style_axes(ax, theme: str) -> None:
    if theme == "white":
        ax.set_facecolor("white")
        ax.figure.set_facecolor("white")
        color = "#111111"
        spine_color = "#444444"
    else:
        ax.set_facecolor("none")
        ax.figure.set_facecolor("none")
        color = "#f5f5f5"
        spine_color = "#888888"
    ax.tick_params(axis="both", labelsize=12, colors=color)
    ax.set_ylabel("Accuracy", fontsize=14, color=color)
    for spine in ax.spines.values():
        spine.set_color(spine_color)
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.5, color=spine_color)
    ax.yaxis.label.set_color(color)
    ax.tick_params(colors=color)
    ax.set_ylim(0, 1.05)


def _compute_styles(models_order: List[str]) -> Dict[str, Tuple[str, str | None]]:
    styles: Dict[str, Tuple[str, str | None]] = {}
    other_idx = 0
    for model in models_order:
        if model.lower().startswith("alvessa"):
            styles[model] = (ALVESSA_COLOR, None)
        else:
            color = OTHER_COLORS[other_idx % len(OTHER_COLORS)]
            hatch = OTHER_HATCHES[other_idx % len(OTHER_HATCHES)]
            styles[model] = (color, hatch)
            other_idx += 1
    return styles


def _plot_by_group(data: Dict[str, pd.DataFrame], out_dir: Path, theme: str) -> Path:
    # Collect all groups and define order (prefer Alvessa ordering if present)
    all_groups = set()
    for df in data.values():
        all_groups.update(df["question_group"].unique())
    if "Alvessa" in data:
        order_df = data["Alvessa"].copy()
        order_df = order_df.sort_values("accuracy")
        ordered = [g for g in order_df["question_group"].tolist() if g in all_groups]
        for g in all_groups:
            if g not in ordered:
                ordered.append(g)
        groups = ordered
    else:
        groups = sorted(all_groups)

    # Ensure "All" stays first if present
    if "All" in groups:
        groups = ["All"] + [g for g in groups if g != "All"]

    x = np.arange(len(groups))
    fig_width = max(10, len(groups) * 0.6)
    fig, ax = plt.subplots(figsize=(fig_width, 5), dpi=300)

    models_order = list(data.keys())
    styles = _compute_styles(models_order)
    n_models = len(models_order)
    width = min(0.8 / max(1, n_models), 0.15)

    text_color = "#111111" if theme == "white" else "#F5F5F5"
    outline_color = "white" if theme == "white" else "#FFFFFF"

    rounded_patches = []
    for idx, model in enumerate(models_order):
        df = data[model]
        acc_map = dict(zip(df["question_group"], df["accuracy"]))
        vals = [acc_map.get(g, 0.0) for g in groups]
        color, hatch = styles.get(model, ("#888888", "//"))
        positions = x + (idx - (n_models - 1) / 5) * width * 1.3

        base_bars = ax.bar(
            positions,
            vals,
            width=width,
            color=color,
            edgecolor="none",
            alpha=0.95,
            label=model,
        )

        for bar in base_bars:
            height = bar.get_height()
            if height is None or height < 1e-6:
                continue
            bbox = bar.get_bbox()
            p = FancyBboxPatch(
                (bbox.xmin, bbox.ymin),
                bbox.width,
                bbox.height,
                boxstyle="round,pad=0.02",
                linewidth=1.0,
                facecolor=color,
                edgecolor=outline_color,
                hatch=hatch or "",
            )
            bar.remove()
            rounded_patches.append(p)

    for p in rounded_patches:
        ax.add_patch(p)

    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=45, ha="right", fontsize=11, color=text_color)
    _style_axes(ax, theme)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#A0A0A0")
    ax.spines["bottom"].set_color("#A0A0A0")
    ax.legend(fontsize=11, loc="upper right", bbox_to_anchor=(1.15, 1.0))

    fname = "benchmark_dbqa_by_group_white.png" if theme == "white" else "benchmark_dbqa_by_group_black.png"
    out_path = out_dir / fname
    if theme == "white":
        fig.savefig(out_path, dpi=300, bbox_inches="tight", transparent=False, facecolor="white")
    else:
        fig.savefig(out_path, dpi=300, bbox_inches="tight", transparent=True)
    plt.close(fig)
    return out_path


def _plot_overall(data: Dict[str, pd.DataFrame], out_dir: Path, theme: str) -> Path:
    labels = list(data.keys())
    styles = _compute_styles(labels)

    accuracies = []
    for name, df in data.items():
        # Prefer explicit "All" row; fallback to mean of is_correct
        all_row = df[df["question_group"] == "All"]
        if not all_row.empty:
            accuracies.append(all_row["accuracy"].iloc[0])
        else:
            accuracies.append(df["accuracy"].mean() if not df.empty else 0.0)

    text_color = "#111111" if theme == "white" else "#F5F5F5"
    outline_color = "white" if theme == "white" else "#FFFFFF"

    fig, ax = plt.subplots(figsize=(WIDTH_PLOT, 5.0), dpi=300) # (max(WIDTH_PLOT, len(labels) * 0.8)
    width = 0.8
    colors = []
    hatches = []
    for name in labels:
        color, hatch = styles.get(name, ("#888888", "//"))
        colors.append(color)
        hatches.append(hatch)
    x_pos = np.arange(len(labels))
    base_bars = ax.bar(x_pos, accuracies, color=colors, edgecolor="none", width=width)

    new_patches = []
    for bar, color, hatch in zip(base_bars, colors, hatches):
        height = bar.get_height()
        if height is None or height < 1e-6:
            continue
        bbox = bar.get_bbox()
        p = FancyBboxPatch(
            (bbox.xmin, bbox.ymin),
            bbox.width,
            bbox.height,
            boxstyle="round,pad=0.00",
            linewidth=1.0,
            facecolor=color,
            edgecolor=outline_color,
            hatch=hatch or "",
        )
        bar.remove()
        new_patches.append(p)

    for p in new_patches:
        ax.add_patch(p)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=0, ha="center", fontsize=11, color=text_color)

    for x, acc, patch in zip(x_pos, accuracies, new_patches):
        height = patch.get_height()
        ax.text(
            x,
            height + 0.005,
            f"{acc:.2f}",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
            color="black" if theme == "white" else "white",
        )

    _style_axes(ax, theme)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#A0A0A0")
    ax.spines["bottom"].set_color("#A0A0A0")

    fname = "benchmark_dbqa_overall_white.png" if theme == "white" else "benchmark_dbqa_overall_black.png"
    out_path = out_dir / fname
    if theme == "white":
        fig.savefig(out_path, dpi=300, bbox_inches="tight", transparent=False, facecolor="white")
    else:
        fig.savefig(out_path, dpi=300, bbox_inches="tight", transparent=True)
    plt.close(fig)
    return out_path


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Visualize dbQA benchmark accuracy.")
    parser.add_argument(
        "--remove_disgenet",
        action="store_true",
        help="If set, drop any questions mentioning DisGeNet before computing accuracies.",
    )
    args = parser.parse_args(argv)

    figures_dir = Path("results/benchmark_figures").resolve()
    figures_dir.mkdir(parents=True, exist_ok=True)

    model_dfs: Dict[str, pd.DataFrame] = {}
    
    N_total = 0

    for name, path_str in MODEL_FILES.items():
        p = Path(path_str).expanduser()
        if not p.exists():
            print(f"[visualize] Missing file for {name}: {p}")
            continue
        try:
            df = _load_benchmark(p, drop_disgenet=args.remove_disgenet)
            model_dfs[name] = _compute_by_group(df)
            N_total += len(df)
        except Exception as exc:
            print(f"[visualize] Failed to load {p} for {name}: {exc}")
            continue

    if not model_dfs:
        print("[visualize] No data loaded; aborting.")
        return 1

    white = _plot_by_group(model_dfs, figures_dir, theme="white")
    black = _plot_by_group(model_dfs, figures_dir, theme="black")
    print(f"Saved dbQA group plots: {white}, {black}")

    white_overall = _plot_overall(model_dfs, figures_dir, theme="white")
    black_overall = _plot_overall(model_dfs, figures_dir, theme="black")
    print(f"Saved dbQA overall plots: {white_overall}, {black_overall}")
    
    print(f"Total questions across all models: {N_total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
