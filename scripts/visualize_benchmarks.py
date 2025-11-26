from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd

# Hardwired benchmark CSVs (extendable)
MODEL_FILES = {
    "Alvessa": "/Users/sokolova/Documents/research/alvessa_agent/out/alvessa_half_GA_20251124-012201_cli/benchmark_summary.csv",
    "Claude\nSonnet 4.5": "/Users/sokolova/Documents/research/alvessa_agent/results/benchmark_runs_done/claude_baseline_GenomeArena_20251124-1607.csv",
    "ChatGPT 5.1": "/Users/sokolova/Documents/research/alvessa_agent/results/benchmark_runs_done/chatgpt_baseline_GenomeArena_20251124-1750.csv"
}

ALVESSA_COLOR = "#D95F02"
# Palette/hatches for non-Alvessa models (cycled in order of appearance)
OTHER_COLORS = ["#555555", "#727272", "#8C8C8C", "#A6A6A6", "#BEBEBE"]
OTHER_HATCHES = ["//", "\\\\", "xx", "..", "++", "--"]

FOLDER_NAME_MAP = {
    "BioGRID": "BioGRID",
    "chembl": "ChEMBL",
    "gencode": "GENCODE",
    "GWAS": "GWAS",
    "gwas+alphamissense": "GWAS + AlphaMissense",
    "miRDB": "miRDB",
    "MSigDB": "MSigDB",
    "OMIM": "OMIM",
    "OpenTargets": "Open Targets",
    "prot": "Protein features",
    "reactome": "Reactome",
}


def _format_folder(name: str) -> str:
    return FOLDER_NAME_MAP.get(name, name)


def _load_benchmark(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    # normalize correctness to int
    df["is_correct"] = pd.to_numeric(df.get("is_correct", 0), errors="coerce").fillna(0).astype(int)
    df["source_folder"] = df.get("source_folder", "").fillna("").astype(str).str.strip()
    df["source_name"] = df.get("source_name", "").fillna("").astype(str).str.strip()
    df["question"] = df.get("question", "").fillna("").astype(str)
    return df


def _compute_by_folder(df: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        df.groupby("source_folder")["is_correct"]
        .agg(["mean", "count"])
        .rename(columns={"mean": "accuracy", "count": "n"})
        .reset_index()
    )
    grouped["display"] = grouped["source_folder"].apply(_format_folder)
    return grouped


def _compute_by_folder_set(df: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        df.groupby(["source_folder", "source_name"])["is_correct"]
        .agg(["mean", "count"])
        .rename(columns={"mean": "accuracy", "count": "n"})
        .reset_index()
    )
    grouped["display_folder"] = grouped["source_folder"].apply(_format_folder)
    return grouped


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
    ax.tick_params(colors=color)
    ax.yaxis.label.set_color(color)


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


def _plot_by_folder(data: Dict[str, pd.DataFrame], out_dir: Path, theme: str) -> Path:
    # Align folders across models; order by Alvessa accuracy low->high if available, else alphabetical
    all_folders = set()
    for df in data.values():
        all_folders.update(df["source_folder"].unique())
    folders = list(all_folders)

    if "Alvessa" in data:
        order_df = data["Alvessa"].sort_values("accuracy")
        ordered = [f for f in order_df["source_folder"].tolist() if f in all_folders]
        for f in folders:
            if f not in ordered:
                ordered.append(f)
        folders = ordered
    else:
        folders = sorted(folders)

    x = range(len(folders))
    fig, ax = plt.subplots(figsize=(max(15, len(folders) * 0.6), 4.5))

    models_order = list(data.keys())
    styles = _compute_styles(models_order)
    n_models = len(models_order)
    width = min(0.8 / max(1, n_models), 0.28)
    for idx, model in enumerate(models_order):
        df = data[model]
        acc_map = dict(zip(df["source_folder"], df["accuracy"]))
        vals = [acc_map.get(f, 0.0) for f in folders]
        color, hatch = styles.get(model, ("#888888", "//"))
        positions = [p + (idx - (n_models - 1) / 2) * width for p in x]
        ax.bar(positions, vals, width=width, color=color, edgecolor="black", hatch=hatch, label=model)

    ax.set_xticks(x)
    ax.set_xticklabels([_format_folder(f) for f in folders], rotation=90, ha="center", fontsize=11, color=("#111111" if theme == "white" else "#f5f5f5"))
    _style_axes(ax, theme)
    ax.legend(fontsize=11)

    fname = "benchmark_by_folder_white.png" if theme == "white" else "benchmark_by_folder_black.png"
    out_path = out_dir / fname
    if theme == "white":
        fig.savefig(out_path, dpi=300, bbox_inches="tight", transparent=False, facecolor="white")
    else:
        fig.savefig(out_path, dpi=300, bbox_inches="tight", transparent=True)
    plt.close(fig)
    return out_path


def _plot_by_folder_set(data: Dict[str, pd.DataFrame], out_dir: Path, theme: str) -> Path:
    # Build combined sorted positions by folder (aligned ordering as above)
    all_folders = set()
    for df in data.values():
        all_folders.update(df["source_folder"].unique())
    if "Alvessa" in data:
        order_df = data["Alvessa"].sort_values("accuracy")
        folders = [f for f in order_df["source_folder"].tolist() if f in all_folders]
        for f in all_folders:
            if f not in folders:
                folders.append(f)
    else:
        folders = sorted(all_folders)

    # Deduplicate folders by display label to avoid repeated buckets (e.g., trailing-space variants)
    seen_disp = set()
    deduped_folders: List[str] = []
    for f in folders:
        disp = _format_folder(f)
        if disp in seen_disp:
            continue
        seen_disp.add(disp)
        deduped_folders.append(f)
    folders = deduped_folders

    def _build_positions(folders_subset: List[str]) -> Tuple[List[float], List[str], List[Tuple[str, float, float]], List[str]]:
        positions: List[float] = []
        labels: List[str] = []
        folder_bounds: List[Tuple[str, float, float]] = []
        folder_order_flat: List[str] = []
        current_x = 0.0
        gap_sets = 0.5
        gap_folders = 1.2
        for f in folders_subset:
            names: List[str] = []
            for df in data.values():
                names.extend(df[df["source_folder"] == f]["source_name"].tolist())
            seen: List[str] = []
            for n in names:
                if n not in seen:
                    seen.append(n)
            if not seen:
                continue
            start = current_x
            for name in seen:
                positions.append(current_x)
                labels.append(name.replace(".csv", ""))
                folder_order_flat.append(f)
                current_x += gap_sets
            end = current_x - gap_sets if seen else current_x
            folder_bounds.append((f, start, end))
            current_x += gap_folders
        return positions, labels, folder_bounds, folder_order_flat

    styles = _compute_styles(list(data.keys()))
    width = min(0.8 / max(1, len(data)), 0.28)
    n_models = max(1, len(data))

    # split folders into two halves to reduce horizontal sprawl
    mid = (len(folders) + 1) // 2
    halves = [folders[:mid], folders[mid:]]

    # determine overall figure width based on the widest half
    half_lengths = []
    for half in halves:
        positions, _, _, _ = _build_positions(half)
        half_lengths.append(len(positions))
    fig_width = max(14, max(half_lengths) * 0.9 if half_lengths else 0)
    fig, axes = plt.subplots(2, 1, figsize=(fig_width, 9.0), constrained_layout=True)

    for ax, folders_subset in zip(axes, halves):
        positions_base, labels, folder_bounds, folder_order_flat = _build_positions(folders_subset)
        if not positions_base:
            ax.axis("off")
            continue
        for idx, (model, df) in enumerate(data.items()):
            acc_map = {(r["source_folder"], r["source_name"]): r["accuracy"] for _, r in df.iterrows()}
            vals = []
            for f, name in zip(folder_order_flat, labels):
                key = (f, name if name.endswith(".csv") else name + ".csv")
                val = acc_map.get(key)
                if val is None:
                    val = acc_map.get((f, name))
                vals.append(val if val is not None else 0.0)
            color, hatch = styles.get(model, ("#888888", "//"))
            offset_positions = [p + (idx - (n_models - 1) / 2) * width for p in positions_base]
            ax.bar(
                offset_positions,
                vals,
                width=width,
                color=color,
                edgecolor="black",
                hatch=hatch,
                alpha=0.9,
                label=model,
            )

        for xpos, lbl in zip(positions_base, labels):
            ax.text(
                xpos,
                -0.03,
                lbl,
                ha="center",
                va="top",
                fontsize=9.5,
                rotation=90,
                color=("#111111" if theme == "white" else "#f5f5f5"),
            )

        for name, start, end in folder_bounds:
            span_pad = width * (n_models - 1) / 2
            adj_start = start - span_pad
            adj_end = end + span_pad
            midpt = (adj_start + adj_end) / 2 if adj_end >= adj_start else adj_start
            y_bracket = -0.07
            color = "#000000" if theme == "white" else "#aaaaaa"
            ax.hlines(y_bracket, adj_start, adj_end, colors=color, linewidth=1.0)
            ax.vlines([adj_start, adj_end], y_bracket, y_bracket - 0.015, colors=color, linewidth=1.0)
            ax.text(midpt, y_bracket - 0.035, _format_folder(name), ha="center", va="top", fontsize=10, color=color)

        ax.set_xticks([])
        ax.set_ylim(0, 1.05)
        _style_axes(ax, theme)

    # single legend centered on top axes
    handles, labels_legend = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].legend(handles, labels_legend, fontsize=11, loc="upper right")

    fname = "benchmark_by_folder_set_white.png" if theme == "white" else "benchmark_by_folder_set_black.png"
    out_path = out_dir / fname
    if theme == "white":
        fig.savefig(out_path, dpi=300, bbox_inches="tight", transparent=False, facecolor="white")
    else:
        fig.savefig(out_path, dpi=300, bbox_inches="tight", transparent=True)
    plt.close(fig)
    return out_path


def _plot_overall(data: Dict[str, pd.DataFrame], out_dir: Path, theme: str) -> Path:
    labels = list(data.keys())
    accuracies = [df["accuracy"].mean() if not df.empty else 0.0 for df in data.values()]

    fig, ax = plt.subplots(figsize=(max(3, len(labels) * 1.1), 4.0))
    width = min(0.95 / max(1, len(labels)), 0.6)
    styles = _compute_styles(labels)
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
    ax.set_ylim(0, 1.05)

    fname = "benchmark_overall_white.png" if theme == "white" else "benchmark_overall_black.png"
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

    model_dfs_folder: Dict[str, pd.DataFrame] = {}
    model_dfs_set: Dict[str, pd.DataFrame] = {}

    for name, path_str in MODEL_FILES.items():
        p = Path(path_str).expanduser()
        if not p.exists():
            print(f"[visualize] Missing file for {name}: {p}")
            continue
        try:
            df = _load_benchmark(p)
            model_dfs_folder[name] = _compute_by_folder(df)
            model_dfs_set[name] = _compute_by_folder_set(df)
        except Exception as exc:
            print(f"[visualize] Failed to load {p} for {name}: {exc}")
            continue

    if not model_dfs_folder:
        print("[visualize] No data loaded; aborting.")
        return 1

    white1 = _plot_by_folder(model_dfs_folder, figures_dir, theme="white")
    black1 = _plot_by_folder(model_dfs_folder, figures_dir, theme="black")
    print(f"Saved folder plots: {white1}, {black1}")

    white2 = _plot_by_folder_set(model_dfs_set, figures_dir, theme="white")
    black2 = _plot_by_folder_set(model_dfs_set, figures_dir, theme="black")
    print(f"Saved set plots: {white2}, {black2}")

    white3 = _plot_overall(model_dfs_folder, figures_dir, theme="white")
    black3 = _plot_overall(model_dfs_folder, figures_dir, theme="black")
    print(f"Saved overall plots: {white3}, {black3}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
