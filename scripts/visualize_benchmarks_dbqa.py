from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd

# Configure benchmark CSVs to compare (model name -> CSV path)
MODEL_FILES = {
    "Alvessa": "/Users/sokolova/Documents/research/alvessa_agent/out/20251215-235658_cli/benchmark_summary.csv",
    "Claude\nSonnet 4.5": "/Users/sokolova/Documents/research/alvessa_agent/claude_baseline_dbQA_20251215-2308.csv",
    "ChatGPT 5.1": "/Users/sokolova/Documents/research/alvessa_agent/chatgpt_baseline_dbQA_20251215-2330.csv"
}

ALVESSA_COLOR = "#D95F02"
# Palette/hatches for non-Alvessa models (cycled in order of appearance)
OTHER_COLORS = ["#555555", "#727272", "#8C8C8C", "#A6A6A6", "#BEBEBE"]
OTHER_HATCHES = ["//", "\\\\", "xx", "..", "++", "--"]

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


def _load_benchmark(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df["is_correct"] = pd.to_numeric(df.get("is_correct", 0), errors="coerce").fillna(0).astype(int)
    df["question"] = df.get("question", "").fillna("").astype(str)
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

    x = range(len(groups))
    fig, ax = plt.subplots(figsize=(max(15, len(groups) * 0.7), 4.8))

    models_order = list(data.keys())
    styles = _compute_styles(models_order)
    n_models = len(models_order)
    width = min(0.8 / max(1, n_models), 0.28)
    for idx, model in enumerate(models_order):
        df = data[model]
        acc_map = dict(zip(df["question_group"], df["accuracy"]))
        vals = [acc_map.get(g, 0.0) for g in groups]
        color, hatch = styles.get(model, ("#888888", "//"))
        positions = [p + (idx - (n_models - 1) / 2) * width for p in x]
        ax.bar(positions, vals, width=width, color=color, edgecolor="black", hatch=hatch, label=model)

    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=45, ha="right", fontsize=11, color=("#111111" if theme == "white" else "#f5f5f5"))
    _style_axes(ax, theme)
    ax.legend(fontsize=11)

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

    fig, ax = plt.subplots(figsize=(max(5, len(labels) * 1.2), 4.0))
    width = 0.6
    colors = []
    hatches = []
    for name in labels:
        color, hatch = styles.get(name, ("#888888", "//"))
        colors.append(color)
        hatches.append(hatch)
    bars = ax.bar(labels, accuracies, color=colors, edgecolor="black", width=width)
    for bar, hatch in zip(bars, hatches):
        if hatch:
            bar.set_hatch(hatch)
    ax.set_xticklabels(labels, rotation=0, ha="center", fontsize=12, color=("#111111" if theme == "white" else "#f5f5f5"))
    _style_axes(ax, theme)

    fname = "benchmark_dbqa_overall_white.png" if theme == "white" else "benchmark_dbqa_overall_black.png"
    out_path = out_dir / fname
    if theme == "white":
        fig.savefig(out_path, dpi=300, bbox_inches="tight", transparent=False, facecolor="white")
    else:
        fig.savefig(out_path, dpi=300, bbox_inches="tight", transparent=True)
    plt.close(fig)
    return out_path


def main() -> int:
    figures_dir = Path("results/benchmark_figures").resolve()
    figures_dir.mkdir(parents=True, exist_ok=True)

    model_dfs: Dict[str, pd.DataFrame] = {}

    for name, path_str in MODEL_FILES.items():
        p = Path(path_str).expanduser()
        if not p.exists():
            print(f"[visualize] Missing file for {name}: {p}")
            continue
        try:
            df = _load_benchmark(p)
            model_dfs[name] = _compute_by_group(df)
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
