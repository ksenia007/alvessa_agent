# %%
"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-08-21
Updated: 2025-08-21


Description: 

Simple plot of the results tables"""

import os
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

RESULTS_ROOT = "benchmarks_generation/results"

def collect_results(root=RESULTS_ROOT):
    rows = []
    for dataset in os.listdir(root):
        dpath = os.path.join(root, dataset)
        if not os.path.isdir(dpath):
            continue
        for model in os.listdir(dpath):
            mpath = os.path.join(dpath, model)
            if not os.path.isdir(mpath):
                continue
            for fname in os.listdir(mpath):
                if not fname.endswith(".csv"):
                    continue
                fpath = os.path.join(mpath, fname)
                try:
                    df = pd.read_csv(fpath)
                    if "is_correct" not in df.columns:
                        continue
                    rows.append({
                        "dataset": dataset,
                        "test_set": fname.replace(".csv",""),
                        "model": model,
                        "n": len(df),
                        "accuracy": float(df["is_correct"].mean())
                    })
                except Exception as e:
                    print(f"[warn] {fpath}: {e}")
    return pd.DataFrame(rows)

import numpy as np
import matplotlib.pyplot as plt
import textwrap
from matplotlib.transforms import blended_transform_factory

def plot_by_testset(df, save_path=None):
    if df.empty:
        print("No results found."); return

    # --- Normalize model names & fix plotting order
    name_map = {
        "alvessa": "Alvessa",
        "alvessa_pipeline": "Alvessa",
        "claude": "Claude",
        "anthropic": "Claude",
    }
    df = df.copy()
    df.loc[:, "model_norm"] = df["model"].str.lower().map(name_map).fillna(df["model"])
    model_order = ["Alvessa", "Claude"]

    # --- Aggregate in case multiple files exist per (dataset, test_set, model)
    g = (df.groupby(["dataset", "test_set", "model_norm"])
            .agg(accuracy=("accuracy", "mean"), n=("n", "sum"))
            .reset_index())

    # --- Clean names: dataset = "name", set = "set"
    # Special-case labbench (treat as single-set dataset with empty set label)
    g["dataset_clean"] = g["dataset"].fillna("").astype(str)
    g["set_clean"] = g["test_set"].fillna("").astype(str)
    is_labbench = g["dataset_clean"].str.lower().eq("labbench")
    g.loc[is_labbench, "set_clean"] = np.where(
        g.loc[is_labbench, "set_clean"].str.strip().eq(""),
        "overall",  # or keep "" if you truly want blank tick
        g.loc[is_labbench, "set_clean"]
    )

    # --- Build accuracy and N tables at (dataset,set) granularity
    g["pair"] = list(zip(g["dataset_clean"], g["set_clean"]))
    acc = g.pivot(index="pair", columns="model_norm", values="accuracy")
    Ns  = g.pivot(index="pair", columns="model_norm", values="n")

    # Ensure consistent column order
    for m in model_order:
        if m not in acc.columns:
            acc[m] = np.nan
            Ns[m]  = np.nan
    acc = acc[model_order]
    Ns  = Ns[model_order]

    # --- Compute per-(dataset,set) gain and per-dataset MAX gain
    acc["gain"] = acc["Alvessa"] - acc["Claude"]
    # Handle missing Claude gracefully
    acc["gain"] = acc["gain"].fillna(0.0)
    # Dataset-level max gain
    idx_df = acc.reset_index()
    idx_df["dataset"] = idx_df["pair"].str[0]
    idx_df["set"]     = idx_df["pair"].str[1]
    max_gain_by_dataset = (idx_df.groupby("dataset")["gain"].max().sort_values(ascending=False))

    # --- Order: datasets by descending max gain; within a dataset keep sets by:
    #  1) descending gain, 2) then alphabetical as tiebreaker (stable & readable)
    ordered_pairs = []
    for ds in max_gain_by_dataset.index:
        sub = idx_df[idx_df["dataset"] == ds].copy()
        sub = sub.sort_values(by=["set"], ascending=[True])
        ordered_pairs.extend(list(zip(sub["dataset"], sub["set"])))

    # Reindex acc/Ns to the new order
    acc = acc.loc[ordered_pairs]
    Ns  = Ns.loc[ordered_pairs]

    # --- Labels: show only SET on the tick; dataset name will be drawn as a bracket label
    set_labels = [p[1] for p in acc.index]
    set_labels_wrapped = [textwrap.fill(lbl.replace("_"," "), width=16) for lbl in set_labels]

    # --- Prepare plot
    plt.rcParams.update({
        "axes.titlesize": 22,
        "axes.labelsize": 18,
        "xtick.labelsize": 13,
        "ytick.labelsize": 16,
        "legend.fontsize": 15,
    })

    x = np.arange(len(acc))
    width = 0.36
    fig, ax = plt.subplots(figsize=(16, 7.5))

    # Colors
    colors = {"Alvessa": "darkorange", "Claude": "grey"}

    # Bars
    bars = []
    for i, m in enumerate(model_order):
        vals = acc[m].values
        b = ax.bar(x + (i - 0.5)*width, vals, width,
                   label=m, color=colors.get(m, "steelblue"))
        bars.append((m, b))

    # Axes & grid
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("Accuracy")
    ax.set_xlabel("")
    ax.set_title("Accuracy by Model and Set (grouped by Dataset)", pad=16)
    ax.set_xticks(x)
    ax.set_xticklabels(set_labels_wrapped, rotation=0, ha="center")
    ax.yaxis.grid(True, linestyle="--", alpha=0.35)
    ax.set_axisbelow(True)
    ax.legend(title="Model", title_fontsize=16, frameon=False, ncols=2)

    # --- Annotations on bars (accuracy and n)
    for model_name, bar_container in bars:
        for idx, rect in enumerate(bar_container):
            h = rect.get_height()
            if np.isnan(h):
                continue
            n_val = Ns.iloc[idx][model_name]
            txt = f"{h:.2f}" if (isinstance(n_val, float) and np.isnan(n_val)) else f"{h:.2f}\n(n={int(n_val)})"
            ax.annotate(txt,
                        xy=(rect.get_x() + rect.get_width()/2, h),
                        xytext=(0, 6),
                        textcoords="offset points",
                        ha="center", va="bottom", fontsize=12)

    # --- Draw dataset brackets underneath the x-axis
    # Build contiguous runs of x positions for each dataset (in the ordered_pairs sequence)
    dataset_runs = []
    current_ds, start_i = None, 0
    for i, (ds, st) in enumerate(ordered_pairs):
        if current_ds is None:
            current_ds, start_i = ds, i
        elif ds != current_ds:
            dataset_runs.append((current_ds, start_i, i-1))
            current_ds, start_i = ds, i
    if current_ds is not None:
        dataset_runs.append((current_ds, start_i, len(ordered_pairs)-1))

    # Place brackets a bit below the axis using blended transform (x in data, y in axes fraction)
    trans = ax.get_xaxis_transform()  # x in data coords; y in axes [0..1]
    y_line  = -0.07
    y_text  = -0.13
    # Expand bottom to make space for brackets
    fig.subplots_adjust(bottom=0.28)
    
    map_names = {
        'GWAS': 'GWAS',
        'labbench': 'LabBench',
        'reactome' : 'Reactome',
        'gencode': 'Gencode',
        'biogrid': 'BioGRID',
    }

    for ds, i0, i1 in dataset_runs:
        # horizontal line
        ax.plot([i0-0.4, i1+0.4], [y_line, y_line], transform=trans, clip_on=False, linewidth=1.5, color='grey')
        # short vertical "ticks" at ends (optional: subtle bracket look)
        ax.plot([i0-0.4, i0-0.4], [y_line, y_line+0.02], transform=trans, clip_on=False, linewidth=1.5, color='grey')
        ax.plot([i1+0.4, i1+0.4], [y_line, y_line+0.02], transform=trans, clip_on=False, linewidth=1.5, color='grey')
        # centered dataset name
        x_mid = (i0 + i1) / 2.0
        ax.text(x_mid, y_text, map_names[ds], #ds.replace("_"," "),
                transform=trans, ha="center", va="top", fontsize=13)

    fig.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=220, bbox_inches="tight")
        print(f"Plot saved to {save_path}")
    else:
        plt.show()


def plot_split(results, save_prefix="accuracy"):
    """
    Make two plots: one for labbench only, one for everything else.
    """
    if results.empty:
        print("No results to plot.")
        return
    
    if not os.path.exists("figures"):
        os.makedirs("figures")
    
    # Plot all together
    print("Plotting all results...")
    plot_by_testset(results, save_path=f"figures/{save_prefix}_all.png")
    
    
    # Split
    lab_df = results[results["dataset"].str.lower().str.contains("labbench")]
    other_df = results[~results["dataset"].str.lower().str.contains("labbench")]

    if not other_df.empty:
        print("Plotting non-LabBench...")
        plot_by_testset(other_df, save_path=f"figures/{save_prefix}_nonlabbench.png")

    if not lab_df.empty:
        print("Plotting LabBench...")
        plot_by_testset(lab_df, save_path=f"figures/{save_prefix}_labbench.png")
        
if __name__ == "__main__":
    results = collect_results()
    plot_split(results, save_prefix="accuracy_by_testset")


# %%
