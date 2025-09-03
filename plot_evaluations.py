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
import seaborn as sns

RESULTS_ROOT = "benchmarks_generation"

def collect_results(root_base=RESULTS_ROOT):
    rows = []
    full_df = pd.DataFrame()
    root = root_base+'/results'
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
                    df['model'] = model
                    df['dataset'] = dataset
                    df['test_set'] = fname.replace(".csv","")
                    
                    # reference main file w/ models needed
                    ref = pd.read_csv(os.path.join(root_base, dataset, fname))
                    df['need_tool'] = ref.iloc[0]['tool']
                    
                    full_df = pd.concat([full_df, df], ignore_index=True)
                except Exception as e:
                    print(f"[warn] {fpath}: {e}")
    return pd.DataFrame(rows), full_df

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
    g = df
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
        'GWAS_AM': 'Variant pathogenicity',
        'labbench': 'LabBench',
        'reactome' : 'Reactome',
        'gencode': 'Gencode',
        'biogrid': 'BioGRID',
        'mirDB': 'miRNA targets',
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


def plot_by_dataset(df, save_path=None, figure_size=(9,6), max_y=1):
    if df.empty:
        print("No results found.")
        return

    # --- Normalize model names
    name_map = {
        "alvessa": "Alvessa",
        "alvessa_pipeline": "Alvessa",
        "claude": "Claude",
        "anthropic": "Claude",
    }
    df = df.copy()
    df["model_norm"] = df["model"].str.lower().map(name_map).fillna(df["model"])
    model_order = ["Alvessa", "Claude"]

    # --- Aggregate mean accuracy per dataset/model
    g = (df.groupby(["dataset", "model_norm"])
            .agg(accuracy=("accuracy", "mean"), n=("n", "sum"))
            .reset_index())

    # --- Pivot to dataset x model table
    acc = g.pivot(index="dataset", columns="model_norm", values="accuracy")
    for m in model_order:
        if m not in acc.columns:
            acc[m] = np.nan
    acc = acc[model_order]

    # --- Compute gain and sort datasets
    acc["gain"] = acc["Alvessa"] - acc["Claude"]
    acc = acc.sort_values("gain", ascending=False)

    # --- Map dataset names to nicer labels
    map_names = {
        'GWAS': 'Trait\nassociation',
        'GWAS_AM': 'Variant\nannotation',
        'labbench': 'dbQA LabBench',
        'reactome': 'Pathways',
        'gencode': 'Gene\nannotations',
        'biogrid': 'Interactions',
        'mirDB': 'miRNA\ntargets',
    }
    acc.index = acc.index.map(lambda d: map_names.get(d, d))

    # --- Plot
    plt.rcParams.update({
        "axes.titlesize": 18,
        "axes.labelsize": 15,
        "xtick.labelsize": 12,
        "ytick.labelsize": 13,
        "legend.fontsize": 12,
    })

    x = np.arange(len(acc))
    width = 0.35
    fig, ax = plt.subplots(figsize=figure_size)

    colors = {"Alvessa": "darkorange", "Claude": "grey"}

    for i, m in enumerate(model_order):
        ax.bar(x + (i - 0.5) * width, acc[m], width,
               label=m, color=colors.get(m, "steelblue"))

    ax.set_ylim(0, max_y+0.1)
    ax.set_ylabel("Accuracy")
    ax.set_title("", pad=14)
    ax.set_xticks(x)
    ax.set_xticklabels(acc.index, rotation=0, ha="center")
    # ax.legend(title="Model", frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=220, bbox_inches="tight")
        print(f"Plot saved to {save_path}")
    else:
        plt.show()

def convert_full_accuracy(full_df, group_by_set=False):
    
    if full_df.empty or "is_correct" not in full_df.columns:
        print("No valid results to convert.")
        return pd.DataFrame()

    df = full_df.copy()
    df["model_norm"] = df["model"].str.lower().map({
        "alvessa": "Alvessa",
        "alvessa_pipeline": "Alvessa",
        "claude": "Claude",
        "anthropic": "Claude",
    }).fillna(df["model"])


    if group_by_set:
        group_cols = ["dataset", "test_set", "model"]
    else:
        group_cols = ["dataset", "model"]
        
    agg_df = (
        df.groupby(group_cols, dropna=False)
          .agg(
              n=("is_correct", "size"),          # N questions
              accuracy=("is_correct", "mean")    # mean correctness
          )
          .reset_index()
          .sort_values(group_cols)
          .reset_index(drop=True)
    )
    return agg_df
      
def plot_heatmap_by_set(df):
    # --- Manual descriptionss
    desc_map = {
        ("GWAS","set1"): "GWAS set",
        ("GWAS","set2"): "GWAS set2",
        ("GWAS","set3"): "GWAS set3",
        ("GWAS_AM","set1"): "Variant pathogenicity, set1",
    }

    df["desc"] = df.apply(lambda r: desc_map.get((r["dataset"], r["test_set"]), 
                                                f"{r['dataset']}:{r['test_set']}"), axis=1)

    # --- Pivot to desc × model
    pivot = df.pivot(index="desc", columns="model", values="accuracy")

    # --- Heatmap
    plt.figure(figsize=(6,4))
    sns.heatmap(pivot, annot=True, fmt=".2f", cmap="YlOrRd", cbar_kws={"label": "Accuracy"})
    plt.title("Accuracy Heatmap by Dataset/Set and Model")
    plt.xlabel("Model")
    plt.ylabel("Dataset + Set")
    plt.tight_layout()
    plt.show()

def tool_selection(df):
    # --- Check if the needed tool was called
    print(df[['need_tool', 'used_models']])
    
    tool_selection_results = {}
    
    mapping_tool_names_acceptable = {
        'GWAS': ['query_gwas_extensive', 'query_gwas_by_gene'],
        'query_gwas_extensive': ['query_gwas_by_gene']
    }
    
    for i, row in df[['need_tool', 'used_models']].iterrows():
        
        # convert need_tool and used_models into lists, [' and splitting
        need_list = [x.strip() for x in str(row['need_tool']).replace('[','').replace(']','').replace("'",'').split(',')]
        used_list = [x.strip() for x in str(row['used_models']).replace('[','').replace(']','').replace("'",'').split(',')]
        
        missing = False
        for n in need_list:
            if n not in used_list:
                if n in mapping_tool_names_acceptable:
                    acc_names = mapping_tool_names_acceptable[n]
                    if not any(a in used_list for a in acc_names):
                        missing = True
                        break
                else:
                    missing = True
                    break
        got_all_needed = not missing
                            
                    
        print(f"Row {i}: Need {need_list}, Used {used_list} => Got all needed: {got_all_needed}")
        
        if not missing:
            num_extra = len(used_list) - len(need_list)
            if 'extract_entities' in used_list:
                num_extra -= 1
            if 'variant_annotations' in used_list and ('alphamissense' in used_list or 'query_gwas_by_gene' in used_list):
                num_extra -= 1
            if 'query_gwas_extensive' in used_list and 'query_gwas_by_gene' in used_list:
                num_extra -= 1
        else:
            num_extra = 0
            
        print(f"Row {i}: Number of extra tools called: {num_extra}")
        
        tool_selection_results[i] = {
            'need_tool': need_list,
            'need_tool_str': str(need_list),
            'used_models': used_list,
            'got_all_needed': got_all_needed,
            'num_extra': num_extra
        }
        
    tool_df = pd.DataFrame.from_dict(tool_selection_results, orient='index')    
    # print(tool_df)
    
    # accuracy for got_all_needed 
    print("Accuracy for calling the right tool", tool_df.got_all_needed.mean())
    print("Mean # extra tools called", tool_df.num_extra.mean())
    print(tool_df.groupby(['need_tool_str']).got_all_needed.mean())
    print(tool_df.groupby(['need_tool_str']).num_extra.mean())
    
    # simple barplot, grouped by need_tool_str, showing got_all_needed and num_extra
    overall_acc   = tool_df.got_all_needed.mean()
    overall_extra = tool_df.num_extra.mean()

    # --- Per-tool
    per_tool = (
        tool_df.groupby("need_tool_str")
            .agg(acc=("got_all_needed","mean"),
                    extra=("num_extra","mean"))
    )

    # --- Plot
    fig, axes = plt.subplots(1, 2, figsize=(12,5), gridspec_kw={'width_ratios':[1,2]})

    # Left: overall
    axes[0].bar(["Got all needed", "Mean # extras"], [overall_acc, overall_extra],
                color=["steelblue","olive"])
    axes[0].set_ylim(0, max(1, overall_extra+0.5))
    axes[0].set_title("Overall")
    axes[0].set_ylabel("Score")
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)

    # Right: per tool grouped bars
    x = range(len(per_tool))
    width = 0.35
    axes[1].bar([i - width/2 for i in x], per_tool["acc"], width, label="Got all needed", color="steelblue")
    axes[1].bar([i + width/2 for i in x], per_tool["extra"], width, label="Mean # extras", color="olive")

    # Define mapping from raw tool names to nicer labels
    label_map = {
        "['BioGRID']": "BioGRID",
        "['GWAS']": "GWAS",
        "['gencode_gene_node']": "Gencode",
        "['miRDB']": "miRNA\ntargets",
        "['query_gwas_extensive', 'alphamissense']": "GWAS+\nAlphaMissense",
        "['reactome']": "Reactome",
    }

    # After plotting per-tool bars:
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([label_map.get(lbl, lbl) for lbl in per_tool.index],
                            rotation=0, ha="center")

    # axes[1].set_xticks(x)
    # axes[1].set_xticklabels(per_tool.index, rotation=20, ha="right")
    axes[1].set_ylim(0, max(1, per_tool["extra"].max()+0.5))
    axes[1].set_title("By needed tool")
    axes[1].legend(frameon=False)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)


    plt.tight_layout()
    plt.savefig('figures/tool_selection.png', dpi=220, bbox_inches="tight")

if __name__ == "__main__":
    results, full_df = collect_results()
    
    results = convert_full_accuracy(full_df, group_by_set=True)
    plot_split(results, save_prefix="accuracy_by_testset")
    
    non_set = convert_full_accuracy(full_df, group_by_set=False)
    results_nonlabbench = non_set[~non_set["dataset"].str.lower().str.contains("labbench")]
    plot_by_dataset(results_nonlabbench, save_path="figures/accuracy_by_dataset.png", 
                    figure_size=(9.5,5))
    
    results_labbench = non_set[non_set["dataset"].str.lower().str.contains("labbench")]
    plot_by_dataset(results_labbench, save_path="figures/accuracy_dbqa.png", 
                    figure_size=(3,5), max_y=0.5)
    
    # print(convert_full_accuracy(full_df))
    
    # plot_heatmap_by_set(results)
    print(full_df)
    

    tool_selection(full_df[full_df.method=='alvessa'])



# %%
