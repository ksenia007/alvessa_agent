from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch
from numpy.random import default_rng
import re
import unicodedata


# Hardwired benchmark CSVs (extendable)
MODEL_FILES = {
    "Alvessa": 'results/benchmark_results/FINAL_GA_20251216-162600_cli/benchmark_summary.csv',
    "\nClaude\nSonnet 4.5": "results/benchmark_results/FINAL_claude_baseline_GenomeArena_20251216-0122.csv",
    "ChatGPT 5.1": "results/benchmark_results/FINAL_chatgpt_baseline_GenomeArena_20251216-0122.csv", 
    '\nClaude\nSonnet 4.5+\nsearch': "results/benchmark_results/claude_baseline_web_serch_N5_GenomeArena_20251216-2008.csv",
}

OTHER_HATCHES = [ "\\\\",  "xx", "//", "++", "--"]
WIDTH_PLOT = (4/3)*4
fig_ext = '_LLM'

MODEL_FILES = {
    "Alvessa*": 'results/benchmark_results/FINAL_GA_20251216-162600_cli/benchmark_summary.csv',
    'Biomni*\n\n\n': "results/benchmark_results/biomni_baseline_GA_10_subset_20251218-2235.csv", 
}
OTHER_HATCHES = ["..", "xx", "//", "\\\\",  "++", "--"]
WIDTH_PLOT = (4/3)*2
fig_ext = '_B'

BAR_WIDTH = 0.85  # unified bar width across all plots

MATCH_Q = True  # whether to only use questions that were attempted by all models
ALVESSA_COLOR = "#D95F02"
# Palette/hatches for non-Alvessa models (cycled in order of appearance)
OTHER_COLORS = [ "#727272", "#555555","#8C8C8C", "#A6A6A6", "#BEBEBE"]
# OTHER_HATCHES = ["..", "xx", "//", "\\\\",  "++", "--"]
FOLDER_NAME_MAP = {
    "BioGRID": "BioGRID",
    "chembl": "ChEMBL",
    "gencode": "GENCODE",
    "gencode_gene_node": "GENCODE",
    "GWAS": "GWAS",
    "gwas+alphamissense": "GWAS+\nAlphaMissense",
    "query_gwas_extensive,alphamissense": "GWAS+\nAlphaMissense",
    "miRDB": "miRDB",
    "MSigDB": "MSigDB",
    "OMIM": "OMIM",
    "OpenTargets": "Open Targets",
    "prot": "Protein features",
    "reactome": "Reactome",
    'aa_seq': 'Amino acid \n sequences',
}

BOOTSTRAP_RESAMPLES = 1000
RNG = default_rng(123)


def _bootstrap_ci(values: pd.Series, n_resamples: int = BOOTSTRAP_RESAMPLES, alpha: float = 0.05) -> tuple[float, float]:
    vals = pd.to_numeric(values, errors="coerce").dropna().to_numpy()
    if len(vals) < 2:
        m = float(vals.mean()) if len(vals) else 0.0
        return m, m
    resamples = RNG.choice(vals, size=(n_resamples, len(vals)), replace=True)
    means = resamples.mean(axis=1)
    low = float(np.percentile(means, 100 * (alpha / 2)))
    high = float(np.percentile(means, 100 * (1 - alpha / 2)))
    return low, high


def _format_folder(name: str) -> str:
    return FOLDER_NAME_MAP.get(name, name)

def canon_question(s: str) -> str:
    s = "" if pd.isna(s) else str(s)

    # Unicode normalize (handles curly quotes etc.)
    s = unicodedata.normalize("NFKC", s)

    # normalize whitespace
    s = s.replace("\u00A0", " ")
    s = re.sub(r"\s+", " ", s).strip()

    # strip ONLY trailing quote characters (and whitespace again)
    s = re.sub(r'[\s"\']+$', "", s).strip()
    return s


def _load_benchmark(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.drop_duplicates('question')
    # normalize correctness to int
    df["is_correct"] = pd.to_numeric(df.get("is_correct", 0), errors="coerce").fillna(0).astype(int)
    df["tool_tag"] = df.get("tool_tag", "").fillna("").astype(str).str.strip()
    df["tool_tag"] = df.get("tool_tag", "").fillna("").astype(str).str.strip()
    df["question"] = df.get("question", "").fillna("").astype(str)
    return df


def _compute_by_folder(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for tool, sub in df.groupby("tool_tag"):
        low, high = _bootstrap_ci(sub["is_correct"])
        rows.append(
            {
                "tool_tag": tool,
                "accuracy": sub["is_correct"].mean(),
                "n": len(sub),
                "ci_low": low,
                "ci_high": high,
            }
        )
    grouped = pd.DataFrame(rows)
    grouped["display"] = grouped["tool_tag"].apply(_format_folder)
    return grouped


def _compute_by_folder_set(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for tool, sub in df.groupby("tool_tag"):
        low, high = _bootstrap_ci(sub["is_correct"])
        rows.append(
            {
                "tool_tag": tool,
                "accuracy": sub["is_correct"].mean(),
                "n": len(sub),
                "ci_low": low,
                "ci_high": high,
            }
        )
    grouped = pd.DataFrame(rows)
    grouped["display_folder"] = grouped["tool_tag"].apply(_format_folder)
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
    ax.tick_params(axis="both", labelsize=13, colors=color)
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
        all_folders.update(df["tool_tag"].unique())
    folders = list(all_folders)

    if "Alvessa" in data:
        order_df = data["Alvessa"].sort_values("accuracy")
        ordered = [f for f in order_df["tool_tag"].tolist() if f in all_folders]
        for f in folders:
            if f not in ordered:
                ordered.append(f)
        folders = ordered
    else:
        folders = sorted(folders)

    out_dir.mkdir(parents=True, exist_ok=True)

    x = np.arange(len(folders))
    fig_width = max(12, len(folders) * 0.6)
    fig, ax = plt.subplots(figsize=(fig_width, 6), dpi=300)

    models_order = list(data.keys())
    styles = _compute_styles(models_order)
    n_models = len(models_order)
    width = BAR_WIDTH

    text_color = "#111111" if theme == "white" else "#F5F5F5"
    # match overall: white outline for hatching
    outline_color = "white" if theme == "white" else "#FFFFFF"

    # plot per model; then replace bars with FancyBboxPatch
    all_patches = []  # if you ever want labels later
    for idx, model in enumerate(models_order):
        df = data[model]
        acc_map = dict(zip(df["tool_tag"], df["accuracy"]))
        low_map = dict(zip(df["tool_tag"], df.get("ci_low", [0.0] * len(df))))
        high_map = dict(zip(df["tool_tag"], df.get("ci_high", [0.0] * len(df))))
        vals = [acc_map.get(f, 0.0) for f in folders]
        err_lower = [max(0.0, v - low_map.get(f, v)) for f, v in zip(folders, vals)]
        err_upper = [max(0.0, high_map.get(f, v) - v) for f, v in zip(folders, vals)]

        color, hatch = styles.get(model, ("#888888", "//"))
        positions = x + (idx - (n_models - 1) / 2) * width

        # base bars to get geometry
        base_bars = ax.bar(
            positions,
            vals,
            width=width,
            color=color,
            edgecolor="none",
            alpha=0.95,
            yerr=[err_lower, err_upper],
            capsize=5,
            error_kw=dict(lw=1.3, capsize=5, capthick=3),
            label=model,
        )

        rounded_for_model = []
        for bar, v in zip(base_bars, vals):
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
            rounded_for_model.append(p)

        for p in rounded_for_model:
            ax.add_patch(p)
        all_patches.extend(rounded_for_model)

    # x tick labels (folders)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [_format_folder(f) for f in folders],
        rotation=90,
        ha="center",
        fontsize=13,
        color=text_color,
    )

    # y axis + grid
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("Accuracy", fontsize=12, color=text_color)
    ax.set_xlabel("", labelpad=6)

    ax.tick_params(axis="y", labelsize=10, colors=text_color)
    ax.tick_params(axis="x", labelsize=13, colors=text_color)

    ax.grid(True, axis="y", linestyle="--", linewidth=0.6, alpha=0.6)
    ax.grid(False, axis="x")

    _style_axes(ax, theme)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#A0A0A0")
    ax.spines["bottom"].set_color("#A0A0A0")

    ax.legend(fontsize=11, loc="upper right", bbox_to_anchor=(1.2, 1.0))

    plt.tight_layout()

    fname = f"benchmark_by_folder_white{fig_ext}.png" if theme == "white" else f"benchmark_by_folder_black{fig_ext}.png"
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
        all_folders.update(df["tool_tag"].unique())
    if "Alvessa" in data:
        order_df = data["Alvessa"].sort_values("accuracy")
        folders = [f for f in order_df["tool_tag"].tolist() if f in all_folders]
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
        gap_sets = max(BAR_WIDTH + 0.1, 0.7)
        gap_folders = 1.2
        for f in folders_subset:
            names: List[str] = []
            for df in data.values():
                names.extend(df[df["tool_tag"] == f]["tool_tag"].tolist())
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
    width = BAR_WIDTH
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
    fig, axes = plt.subplots(2, 1, figsize=(fig_width, 10.0), constrained_layout=True)

    for ax, folders_subset in zip(axes, halves):
        positions_base, labels, folder_bounds, folder_order_flat = _build_positions(folders_subset)
        if not positions_base:
            ax.axis("off")
            continue
        for idx, (model, df) in enumerate(data.items()):
            acc_map = {r["tool_tag"]: r["accuracy"] for _, r in df.iterrows()}
            low_map = {r["tool_tag"]: r.get("ci_low", r["accuracy"]) for _, r in df.iterrows()}
            high_map = {r["tool_tag"]: r.get("ci_high", r["accuracy"]) for _, r in df.iterrows()}
            vals = [acc_map.get(f, 0.0) for f in folder_order_flat]
            err_lower = [max(0.0, v - low_map.get(f, v)) for f, v in zip(folder_order_flat, vals)]
            err_upper = [max(0.0, high_map.get(f, v) - v) for f, v in zip(folder_order_flat, vals)]
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
                yerr=[err_lower, err_upper],
                capsize=5,
                error_kw=dict(lw=1.3, capsize=5, capthick=3),
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

    fname = f"benchmark_by_folder_set_white{fig_ext}.png" if theme == "white" else f"benchmark_by_folder_set_black{fig_ext}.png"
    out_path = out_dir / fname
    if theme == "white":
        fig.savefig(out_path, dpi=300, bbox_inches="tight", transparent=False, facecolor="white")
    else:
        fig.savefig(out_path, dpi=300, bbox_inches="tight", transparent=True)
    plt.close(fig)
    return out_path



def _plot_overall(data: Dict[str, pd.DataFrame], raw_data: Dict[str, pd.DataFrame], out_dir: Path, theme: str) -> Path:
    labels = list(data.keys())
    accuracies = []
    ci_lows = []
    ci_highs = []
    stds = []
    for name in labels:
        raw = raw_data.get(name, pd.DataFrame())
        vals = pd.to_numeric(raw.get("is_correct", []), errors="coerce").dropna()
        mean = vals.mean() if not vals.empty else 0.0
        low, high = _bootstrap_ci(vals)
        accuracies.append(mean)
        ci_lows.append(low)
        ci_highs.append(high)
        stds.append(vals.std() if len(vals) > 1 else 0.0)

    out_dir.mkdir(parents=True, exist_ok=True)

    text_color = "#111111" if theme == "white" else "#F5F5F5"
    outline_color = "white" if theme == "white" else "#FFFFFF"

    fig_width = WIDTH_PLOT 
    fig, ax = plt.subplots(figsize=(fig_width, 6.0), dpi=300)

    x_pos = np.arange(len(labels))

    styles = _compute_styles(labels)
    colors = []
    hatches = []
    for name in labels:
        color, hatch = styles.get(name, ("#888888", ""))  # default color + no hatch
        colors.append(color)
        hatches.append(hatch)

    # Base bars just to get positions; we'll replace them
    base_bars = ax.bar(
        x_pos,
        accuracies,
        color=colors,
        edgecolor="none",
        width=BAR_WIDTH,
        yerr=[np.maximum(0, np.array(accuracies) - np.array(ci_lows)), np.maximum(0, np.array(ci_highs) - np.array(accuracies))],
        capsize=5,
        error_kw=dict(lw=2.3, capsize=5, capthick=3),
    )

    # X tick labels
    ax.set_xticks(x_pos)
    ax.set_xticklabels(
        labels,
        rotation=0,
        ha="center",
        fontsize=15,
        color=text_color,
    )

    # ---- Rounded corners + hatching + outline ----
    new_patches = []
    for bar, color, hatch in zip(base_bars, colors, hatches):
        height = bar.get_height()
        if height is None or height < 1e-6:
            continue
        bbox = bar.get_bbox()

        p = FancyBboxPatch(
            (bbox.xmin, bbox.ymin),
            BAR_WIDTH, #bbox.width,
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

    # Y axis setup
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("Accuracy", fontsize=15, color=text_color)
    # ax.set_xlabel("", labelpad=6)

    ax.tick_params(axis="y", labelsize=13, colors=text_color)
    ax.tick_params(axis="x", labelsize=14, colors=text_color)

    ax.grid(True, axis="y", linestyle="--", linewidth=0.6, alpha=0.6)
    ax.grid(False, axis="x")

    # External theming hook
    # _style_axes(ax, theme)

    # Spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#A0A0A0")
    ax.spines["bottom"].set_color("#A0A0A0")

    plt.tight_layout()

    fname = f"benchmark_overall_white{fig_ext}.png" if theme == "white" else f"benchmark_overall_black{fig_ext}.png"
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
    model_raw: Dict[str, pd.DataFrame] = {}
    
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
        questions = df["question"].tolist()
        questions = [canon_question(q) for q in questions]
        common_questions = set.intersection(*(set(questions) for df in all_dfs.values()))
        print(f"[visualize] Found {len(common_questions)} common questions across all models.")
        
        
    N_total = 0

    for name, path_str in MODEL_FILES.items():
        p = Path(path_str).expanduser()
        if not p.exists():
            print(f"[visualize] Missing file for {name}: {p}")
            continue
        try:
            df = _load_benchmark(p).drop_duplicates(subset=["question"])
            if MATCH_Q:
                df['formatted_question'] = df['question'].map(canon_question)
                df = df[df["formatted_question"].isin(common_questions)].reset_index(drop=True)
            N_total += len(df)
            model_raw[name] = df
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

    white3 = _plot_overall(model_dfs_folder, model_raw, figures_dir, theme="white")
    black3 = _plot_overall(model_dfs_folder, model_raw, figures_dir, theme="black")
    print(f"Saved overall plots: {white3}, {black3}")
    
    print(f"[visualize] Processed total of {N_total} questions across {len(model_dfs_folder)} models.")

    print("\n[visualize] Overall accuracy summary (mean and 95% CI):")
    for name, df in model_raw.items():
        vals = pd.to_numeric(df.get("is_correct", []), errors="coerce").dropna()
        mean = vals.mean() if not vals.empty else 0.0
        low, high = _bootstrap_ci(vals)
        print(f"  - {name}: {mean:.3f} (95% CI: {low:.3f}–{high:.3f}), n={len(vals)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
