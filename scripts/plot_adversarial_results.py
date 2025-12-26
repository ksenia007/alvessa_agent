#%%
import json
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

#%%

file_locations = [
    '/Users/sokolova/Documents/research/alvessa_agent/out/20251225-165728_adversarial_eval/baseline-20251225-181505', 
    '/Users/sokolova/Documents/research/alvessa_agent/out/20251225-165722_adversarial_eval/baseline-20251225-181449', 
    '/Users/sokolova/Documents/research/alvessa_agent/out/20251225-151850_adversarial_eval/baseline-20251225-181857', 
    '/Users/sokolova/Documents/research/alvessa_agent/out/20251225-151846_adversarial_eval/baseline-20251225-181942', 
]

# load all jsons from the file locations; 
# create 2 tables - per statement and overall
# overall
# first overall LLM- baseline1 - [baseline_overall_llm][verdict], 
# then LLM w/ per statement verification - baseline2 - [baseline_overall_per_statement][verdict],
# then alvessa overall - [alvessa_overall][verdict]

# then in per statement
# in "statements" key find "is_adversarial", 
# in which case record baseline_label , alvessa_label & modification_type + question ID (from filename) + "original_statement" + "adversarial_statement" + "original_proof"

# overall 
overall_data = []
for file_location in file_locations:
    path = Path(file_location)
    for json_file in path.glob('*.json'):
        with json_file.open('r', encoding='utf-8') as f:
            content = json.load(f)
            baseline1 = content.get('baseline_overall_llm', {}).get('baseline_full', {})
            baseline2 = content.get('baseline_overall_per_statement', {})
            alvessa = content.get('alvessa_overall', {})
            overall_data.append({
                'question_id': json_file.stem,
                'baseline1_verdict': baseline1.get('verdict'),
                'baseline2_verdict': baseline2.get('verdict'),
                'alvessa_verdict': alvessa.get('verdict'),
            })
            
overall_df = pd.DataFrame(overall_data)
overall_df.head()
# %%
# per statement data
data = []
for file_location in file_locations:
    path = Path(file_location)
    for json_file in path.glob('*.json'):
        with json_file.open('r', encoding='utf-8') as f:
            content = json.load(f)
            for statement in content.get('statements', []):
                if statement.get('is_adversarial'):
                    data.append({
                        'question_id': json_file.stem,
                        'baseline_label': statement.get('baseline_label'),
                        'alvessa_label': statement.get('alvessa_label'),
                        'modification_type': statement.get('modification_type'),
                        'original_statement': statement.get('original_statement'),
                        'adversarial_statement': statement.get('adversarial_statement'), 
                        'original_proof': statement.get('original_proof'),
                    })
                    
# convert data to DataFrame
df = pd.DataFrame(data)
df.head()
# %%
# fraction of correct baseline and alvessa labels - if label is false then correct for baseline and 'unsupported' or 'partial' for alvessa
df['baseline_correct'] = df['baseline_label'] == 'fail'
df['alvessa_correct'] = df['alvessa_label'].isin(['unsupported', 'partial', 'speculation-overreach'])

baseline_accuracy = df['baseline_correct'].mean()
alvessa_accuracy = df['alvessa_correct'].mean()
print(f'Baseline Accuracy: {baseline_accuracy:.2%}')
print(f'Alvessa Accuracy: {alvessa_accuracy:.2%}')
# %%
# Plot nicer bar chart
# ------------------------------------------------------------------
# Pretty labels for y-axis
label_map = {
    # "contradiction": "Contradiction\n(plain)",
    "contradiction_with_proofs": "Contradiction",
    "overstatement": "Overstatement",
    "hallucinated_alphanumeric_entities": "Wrong alphanumeric\nvalues",
    "hallucinated_numbers": "Wrong numerical\nvalues",
}

modification_summary = (
    df.groupby("modification_type")
      .agg(
          baseline_correct=("baseline_correct", "mean"),
          alvessa_correct=("alvessa_correct", "mean"),
      )
      .reset_index()
)
modification_summary["pretty_type"] = modification_summary["modification_type"].map(label_map).fillna(modification_summary["modification_type"])
modification_summary = modification_summary.melt(
    id_vars=["modification_type", "pretty_type"],
    var_name="model",
    value_name="accuracy",
)
# convert accuracy to percent
modification_summary['accuracy'] = modification_summary['accuracy']*100

# Styling
sns.set_theme(style="whitegrid")
plt.rcParams.update({
    "font.size": 22,
    "axes.titlesize": 24,
    "axes.labelsize": 22,
    "xtick.labelsize": 20,
    "ytick.labelsize": 22,
})

palette = {
    "baseline_correct": 'lightgrey', #"tan",  
    "alvessa_correct": "#FF7C07",   # vibrant orange
}

fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
hue_order = ["baseline_correct", "alvessa_correct"]
cat_order = modification_summary.sort_values("accuracy", ascending=False)["pretty_type"].unique().tolist()
model_sequence = modification_summary["model"].tolist()  # preserves melted row order
bar = sns.barplot(
    data=modification_summary,
    x="accuracy",
    y="pretty_type",
    hue="model",
    hue_order=hue_order,
    order=cat_order,
    palette=palette,
    ax=ax,
    dodge=True,
    edgecolor="none",
)
# Rounded corners for bars + hatching on baseline
from matplotlib.patches import FancyBboxPatch
new_patches = []
for patch in bar.patches:
    width = patch.get_width()
    if width is None or width < 1e-3:   
        continue
    bbox = patch.get_bbox()
    p = FancyBboxPatch(
        (bbox.xmin, bbox.ymin),
        bbox.width,
        bbox.height,
        boxstyle="round,pad=0.0",
        linewidth=0,
        facecolor=patch.get_facecolor(),
        edgecolor="none",
    )
    patch.remove()
    new_patches.append(p)
for p in new_patches:
    ax.add_patch(p)

# Apply hatching only to baseline bars using known hue order (no outline)
n_hues = len(hue_order)
for idx, patch in enumerate(new_patches):
    # Map by data row order to avoid surprises
    model_key = model_sequence[idx] if idx < len(model_sequence) else hue_order[idx % n_hues]
    # enforce facecolor by model to avoid bleed-through
    patch.set_facecolor(palette.get(model_key, patch.get_facecolor()))
    width = patch.get_width()
    if model_key == "baseline_correct":
        patch.set_hatch("//")
        patch.set_edgecolor("grey")
        patch.set_linewidth(0.8)
        y_center = patch.get_y() + patch.get_height() / 2.0
        print('y_center', y_center)
        ax.text(
            width - 1,              
            y_center,
            f"{width:.1f}%",      
            va="center",
            ha="right",
            fontsize=19,
            fontweight="bold",
            color="black",
        )

    elif model_key == "alvessa_correct":
        patch.set_hatch("")
        patch.set_edgecolor("none")
        y_center = patch.get_y() + patch.get_height() / 2.0
        ax.text(
            width - 1,              
            y_center,
            f"{width:.1f}%",      
            va="center",
            ha="right",
            fontsize=19,
            fontweight="bold",
            color="black",
        )
    

# Axes formatting
ax.set_xlim(0.0, 105)
ax.set_xlabel("% adversarial statements detected", labelpad=8)
ax.set_ylabel("")  # no y-axis label per request
ax.set_title("")   # no title
ax.grid(True, axis="x", linestyle="--", linewidth=0.6, alpha=0.6)
ax.grid(False, axis="y")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_color("#A0A0A0")
ax.spines["bottom"].set_color("#A0A0A0")


# Legend cleanup (manual handles so hatching shows)
legend_handles = [
    Patch(facecolor=palette["baseline_correct"], edgecolor="grey", linewidth=0.8, hatch="///", label="Ablated verifier"),
    Patch(facecolor=palette["alvessa_correct"], edgecolor="none", label="Context-aware verifier"),
]
ax.legend(legend_handles, [h.get_label() for h in legend_handles], title="", 
          frameon=False, loc="upper right",bbox_to_anchor=(1.1, 1.1), fontsize=18)

plt.tight_layout()
plt.show()

# %%
# write into .txt: original, adversarial,  baseline label, alvessa label
label_map = {
    "contradiction_with_proofs": "Contradiction",
    "overstatement": "Overstatement",
    "hallucinated_alphanumeric_entities": "Wrong alphanumeric values",
    "hallucinated_numbers": "Wrong numerical values",
}
txt_output_path = Path('/Users/sokolova/Documents/research/alvessa_agent/results/benchmark_results/adversarial_statements_summary.txt')
with txt_output_path.open('w', encoding='utf-8') as f:
    for idx, row in df.iterrows():
        f.write(f"Question ID: {row['question_id']}\n")
        f.write(f"Modification Type: {label_map.get(row['modification_type'], None)} \n")
        f.write(f"Original Statement: {row['original_statement']}\n")
        f.write(f"Adversarial Statement: {row['adversarial_statement']}\n")
        f.write(f"Baseline Label: {row['baseline_label']}\n")
        f.write(f"Alvessa Label: {row['alvessa_label']}\n")
        f.write("-" * 80 + "\n")

# %%