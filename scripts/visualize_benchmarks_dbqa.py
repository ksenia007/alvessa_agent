from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd

# Configure benchmark CSVs to compare (model name -> CSV path)
MODEL_FILES = {
    "Alvessa": "/Users/sokolova/Documents/research/alvessa_agent/out/sample_dbqa20251123-214810_cli/benchmark_summary.csv",
    "Claude Sonnet 4.5": "/Users/sokolova/Documents/research/alvessa_agent/results/benchmark_runs_done/claude_baseline_dbQA_20251124-1604.csv",
}

ALVESSA_COLOR = "#D95F02"
BASELINE_COLOR = "#888888"

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
    fig, ax = plt.subplots(figsize=(max(10, len(groups) * 0.7), 4.8))

    width = 0.35
    models_order = list(data.keys())
    n_models = len(models_order)
    for idx, model in enumerate(models_order):
        df = data[model]
        acc_map = dict(zip(df["question_group"], df["accuracy"]))
        vals = [acc_map.get(g, 0.0) for g in groups]
        if model.lower().startswith("alvessa"):
            color = ALVESSA_COLOR
            hatch = None
        else:
            color = BASELINE_COLOR
            hatch = "//"
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
