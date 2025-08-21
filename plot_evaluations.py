# %%
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

def plot_by_testset(df, save_path=None):
    if df.empty:
        print("No results found."); return

    # Normalize model names & lock order
    name_map = {
        "alvessa": "Alvessa",
        "alvessa_pipeline": "Alvessa",
        "claude": "Claude",
        "anthropic": "Claude"
    }
    df["model_norm"] = df["model"].str.lower().map(name_map).fillna(df["model"])
    model_order = ["Alvessa", "Claude"]

    # Aggregate in case multiple files exist per (dataset, test_set, model)
    g = (df.groupby(["dataset","test_set","model_norm"])
            .agg(accuracy=("accuracy","mean"), n=("n","sum"))
            .reset_index())

    # Label = dataset – test_set (shortened)
    def short_label(r):
        ds = r["dataset"].replace("_"," ")
        ts = r["test_set"].replace("_"," ")
        return f"{ds} – {ts}"
    g["label"] = g.apply(short_label, axis=1)

    # Pivot to [labels x models] accuracy, and parallel pivot for n
    acc = g.pivot(index="label", columns="model_norm", values="accuracy")
    Ns  = g.pivot(index="label", columns="model_norm", values="n")

    # Ensure consistent column order
    for m in model_order:
        if m not in acc.columns:
            acc[m] = np.nan
            Ns[m]  = np.nan
    acc = acc[model_order]
    Ns  = Ns[model_order]

    # Sort labels for stability (by dataset then test_set alphabetically)
    acc = acc.sort_index()
    Ns  = Ns.loc[acc.index]

    plt.rcParams.update({
        "axes.titlesize": 22,
        "axes.labelsize": 18,
        "xtick.labelsize": 18,   
        "ytick.labelsize": 16,   
        "legend.fontsize": 16
    })

    labels = acc.index.tolist()
    labels_wrapped = [textwrap.fill(lbl, width=28) for lbl in labels]

    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(15, 7))

    # Colors
    colors = {
        "Alvessa": "darkorange",
        "Claude":  "grey"
    }

    # Bars
    bars = []
    for i, m in enumerate(model_order):
        vals = acc[m].values
        b = ax.bar(x + (i - 0.5)*width, vals, width,
                   label=m, color=colors.get(m, "steelblue"))
        bars.append((m, b))

    # Grid & axes
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("Accuracy")
    ax.set_xlabel("")   
    ax.set_title("Accuracy by Model and Test Set", pad=16)
    ax.set_xticks(x)
    ax.set_xticklabels(labels_wrapped, rotation=0, ha="center")
    ax.yaxis.grid(True, linestyle="--", alpha=0.35)
    ax.set_axisbelow(True)
    ax.legend(title="Model", title_fontsize=16, frameon=False, ncols=2)

    # Annotations (accuracy and n)
    for model_name, bar_container in bars:
        for idx, rect in enumerate(bar_container):
            h = rect.get_height()
            if np.isnan(h):
                continue
            n_val = Ns.iloc[idx][model_name]
            ax.annotate(f"{h:.2f}\n(n={int(n_val)})" if not np.isnan(n_val) else f"{h:.2f}",
                        xy=(rect.get_x() + rect.get_width()/2, h),
                        xytext=(0, 6),
                        textcoords="offset points",
                        ha="center", va="bottom", fontsize=13)

    fig.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches="tight")
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
    
    # Split
    lab_df = results[results["dataset"].str.lower().str.contains("labbench")]
    other_df = results[~results["dataset"].str.lower().str.contains("labbench")]
    
    

    if not other_df.empty:
        print("Plotting non-LabBench...")
        plot_by_testset(other_df, save_path=f"{save_prefix}_nonlabbench.png")

    if not lab_df.empty:
        print("Plotting LabBench...")
        plot_by_testset(lab_df, save_path=f"{save_prefix}_labbench.png")
        
if __name__ == "__main__":
    results = collect_results()
    plot_split(results, save_prefix="accuracy_by_testset")


# %%
