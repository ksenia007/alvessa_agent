from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


# Hardwired benchmark CSVs 
MODEL_FILES = {
    "Alvessa*": 'results/benchmark_results/FINAL_GA_20251216-162600_cli/benchmark_summary.csv', 
    "Biomni*": 'results/benchmark_results/biomni_baseline_GA_10_subset_20251218-2235.csv'
}
figure_ext = '_GA'


MODEL_FILES = {
    "Alvessa*": "results/benchmark_results/FINAL_DBQA_20251216-125635_cli/benchmark_summary.csv",
    "Biomni*": "results/benchmark_results/biomni_baseline_dbQA_subset_20251219-1701.csv",
}
figure_ext = '_dbQA'

ALVESSA_COLOR = "#D95F02"
OTHER_COLORS = [ "#727272", "#555555","#8C8C8C", "#A6A6A6", "#BEBEBE"]
OTHER_HATCHES = ["..", "xx", "//", "\\\\",  "++", "--"]
WIDTH_PLOT = (4/3)*2
MATCH_Q = True  # whether to only use questions that were attempted by all models
BAR_WIDTH = 0.9  # unified bar width

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

    ax.tick_params(axis="x", labelsize=14, colors=color)
    ax.tick_params(axis="y", labelsize=13, colors=color)
    ax.set_ylabel("Runtime (sec)", fontsize=15, color=color)
    for spine in ax.spines.values():
        spine.set_color(spine_color)
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.5, color=spine_color)
    ax.yaxis.label.set_color(color)
    ax.tick_params(colors=color)
    ax.set_ylim(bottom=0)

def _load_benchmark(csv_path: Path, *, drop_disgenet: bool = False) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df = df.drop_duplicates('question')
    df["runtime_seconds"] = pd.to_numeric(df.get("runtime_seconds", 0.0), errors="coerce")
    if drop_disgenet:
        mask = ~df["question"].str.lower().str.contains("disgenet", na=False)
        df = df[mask].reset_index(drop=True)
    return df


def _plot_runtime_overall(data: Dict[str, pd.DataFrame], out_dir: Path, theme: str) -> Path:
    labels = list(data.keys())
    styles = _compute_styles(labels)

    means = []
    stds = []
    counts = []
    for name, df in data.items():
        runtimes = df["runtime_seconds"].dropna()
        means.append(runtimes.mean() if not runtimes.empty else 0.0)
        stds.append(runtimes.std() if len(runtimes) > 1 else 0.0)
        counts.append(len(runtimes))

    text_color = "#111111" if theme == "white" else "#F5F5F5"
    outline_color = "white" if theme == "white" else "#FFFFFF"

    fig, ax = plt.subplots(figsize=(WIDTH_PLOT , 5.0), dpi=300) # standardize height with other plots
    x_pos = np.arange(len(labels))
    width = BAR_WIDTH

    colors = []
    hatches = []
    for name in labels:
        color, hatch = styles.get(name, ("#888888", "//"))
        colors.append(color)
        hatches.append(hatch)

    base_bars = ax.bar(
        x_pos,
        means,
        yerr=stds,
        capsize=5,
        color=colors,
        edgecolor="none",
        width=width,
        alpha=0.95,
        error_kw=dict(lw=2.3, capsize=5, capthick=3)
    )

    new_patches = []
    for bar, color, hatch in zip(base_bars, colors, hatches):
        height = bar.get_height()
        bbox = bar.get_bbox()
        p = FancyBboxPatch(
            (bbox.xmin, bbox.ymin),
            bbox.width,
            bbox.height,
            boxstyle="round,pad=0.0",
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
    ax.set_xticklabels(labels, rotation=0, ha="center", fontsize=14, color=text_color)

    # for x, mean, patch in zip(x_pos, means, new_patches):
    #     ax.text(
    #         x,
    #         mean + (0.02 * max(means) if max(means) else 0.05),
    #         f"{mean:.2f}s",
    #         ha="center",
    #         va="bottom",
    #         fontsize=10,
    #         fontweight="bold",
    #         color="black" if theme == "white" else "white",
    #     )

    _style_axes(ax, theme)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#A0A0A0")
    ax.spines["bottom"].set_color("#A0A0A0")

    plt.tight_layout()

    fname = f"benchmark_runtime_white{figure_ext}.png" if theme == "white" else f"benchmark_runtime_black{figure_ext}.png"
    out_path = out_dir / fname
    if theme == "white":
        fig.savefig(out_path, dpi=300, bbox_inches="tight", transparent=False, facecolor="white")
    else:
        fig.savefig(out_path, dpi=300, bbox_inches="tight", transparent=True)
    plt.close(fig)
    return out_path


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Visualize mean runtime per question from benchmark summaries.")
    parser.add_argument(
        "--theme",
        choices=["white", "black"],
        default="white",
        help="Plot theme background (default: white).",
    )
    parser.add_argument(
        "--out_dir",
        default="results/benchmark_figures",
        help="Directory to write figures (default: results/benchmark_figures).",
    )
    parser.add_argument(
        "--remove_disgenet",
        action="store_true",
        help="If set, drop any questions mentioning DisGeNet before computing accuracies.",
    )
    args = parser.parse_args(argv)

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not MODEL_FILES:
        print("[visualize] MODEL_FILES is empty; please populate it with name->csv_path.")
        return 1
    
    if MATCH_Q:
        print("[visualize] MATCH_Q is set to True")
        # preload all dataframes to find common questions
        all_dfs: Dict[str, pd.DataFrame] = {}
        for name, path_str in MODEL_FILES.items():
            p = Path(path_str).expanduser()
            if not p.exists():
                print(f"[visualize] Missing file for {name}: {p}")
                continue
            try:
                df = _load_benchmark(p)
                all_dfs[name] = df
            except Exception as exc:
                print(f"[visualize] Failed to load {p} for {name}: {exc}")
                continue
        if not all_dfs:
            print("[visualize] No data loaded; aborting.")
            return 1
        # find common questions
        common_questions = set.intersection(*(set(df["question"].tolist()) for df in all_dfs.values()))
        print(f"[visualize] Found {len(common_questions)} common questions across all models.")
    

    model_dfs: Dict[str, pd.DataFrame] = {}
    for name, path_str in MODEL_FILES.items():
        p = Path(path_str).expanduser()
        if not p.exists():
            print(f"[visualize] Missing file for {name}: {p}")
            continue
        try:
            df = _load_benchmark(p, drop_disgenet=args.remove_disgenet).drop_duplicates(subset=["question"])
            if MATCH_Q:
                df = df[df["question"].isin(common_questions)].reset_index(drop=True)
                print(f"[visualize] After filtering, {name} has {len(df)} questions.")
            model_dfs[name] = df
        except Exception as exc:
            print(f"[visualize] Failed to load {p} for {name}: {exc}")
            continue

    if not model_dfs:
        print("[visualize] No data loaded; aborting.")
        return 1

    plot_path = _plot_runtime_overall(model_dfs, out_dir, theme=args.theme)
    print(f"Saved runtime plot: {plot_path}")

    print("\nModel runtime stats:")
    for name, df in model_dfs.items():
        rt = pd.to_numeric(df.get("runtime_seconds", []), errors="coerce").dropna()
        mean = rt.mean() if not rt.empty else 0.0
        std = rt.std() if len(rt) > 1 else 0.0
        print(f"  - {name}: mean={mean:.3f}s, std={std:.3f}s, n={len(rt)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
