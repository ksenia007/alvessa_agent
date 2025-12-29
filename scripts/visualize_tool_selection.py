# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
# %%
loc = 'results/benchmark_results/FINAL2_GA_20251224-200223_cli/benchmark_summary.csv'
preds = pd.read_csv(loc)

# --- drop rows where either tool_tag or used_tools is NaN ---
df = preds.dropna(subset=['tool_tag', 'used_tools']).copy()


# df = df[df['tool_tag']!='query_gwas_extensive,alphamissense']
# replace query_gwas_extensive,alphamissense substring in tool_tag with  query_gwas_extensive+alphamissense
df['tool_tag'] = df['tool_tag'].str.replace(
    'query_gwas_extensive,alphamissense',
    'query_gwas_extensive+alphamissense'
)
# count query_gwas_extensive+alphamissense
print("Rows with query_gwas_extensive+alphamissense in tool_tag:",)
print((df['tool_tag'] == 'query_gwas_extensive+alphamissense').sum())
#%%

# --- parse lists ---
df['needed'] = df['tool_tag'].str.split(',').apply(
    lambda lst: [x.strip() for x in lst if x and x.strip()]
)
df['used'] = df['used_tools'].str.split(';').apply(
    lambda lst: [x.strip() for x in lst if x and x.strip()]
)

# unify GWAS naming
mapping = {
    'query_gwas_by_gene': 'GWAS',
    'query_gwas_extensive': 'GWAS',
    'uniprot_gwas': 'UniProt',
    'uniprot_base': 'UniProt',
    'query_gwas_extensive+alphamissense': 'Annot+AlphaMissense',
}

df['needed'] = df['needed'].apply(lambda lst: [mapping.get(x, x) for x in lst])
df['used']   = df['used'].apply(lambda lst: [mapping.get(x, x) for x in lst])

# drop extract_entities
df['needed'] = df['needed'].apply(lambda lst: [x for x in lst if x != 'extract_entities'])
df['used']   = df['used'].apply(lambda lst: [x for x in lst if x != 'extract_entities'])

# drop rows where needed list is empty after cleaning
df = df[df['needed'].map(len) > 0].copy()

needed_tools = sorted(set(sum(df['needed'], [])))   # tools ever REQUIRED
used_tools   = sorted(set(sum(df['used'], [])))     # tools ever USED
used_tools
# %%

# --- indicator matrices ---
needed_mat = pd.DataFrame(0, index=df.index, columns=needed_tools, dtype=int)
used_mat   = pd.DataFrame(0, index=df.index, columns=used_tools, dtype=int)

for t in needed_tools:
    needed_mat[t] = df['needed'].apply(lambda lst: t in lst).astype(int)

for t in used_tools:
    used_mat[t] = df['used'].apply(lambda lst: t in lst).astype(int)

# --- conditional probabilities P(tool used | tool needed) ---
crosstab_prop = pd.DataFrame(index=needed_tools, columns=used_tools, dtype=float)

for t_need in needed_tools:
    mask = needed_mat[t_need] == 1
    denom = mask.sum()
    if denom == 0:
        crosstab_prop.loc[t_need, :] = 0.0
    else:
        crosstab_prop.loc[t_need, :] = used_mat[mask].mean()
# %%
# --- ORDER rows / columns explicitly for the figure ---
row_order = needed_tools 
# move "'GWAS+AlphaMissense'" to the bottom
if 'Annot+AlphaMissense' in row_order:
    row_order.remove('Annot+AlphaMissense')
    row_order.append('Annot+AlphaMissense')
# align beginning of col_order with row_order
col_order = [i for i in row_order if i in used_tools] 
col_order.append('alphamissense')
# add remaining used_tools at the end
for t in used_tools:
    if t not in col_order:
        col_order.append(t)


row_order = [t for t in row_order if t in crosstab_prop.index]
col_order = [t for t in col_order if t in crosstab_prop.columns]
crosstab_prop = crosstab_prop.loc[row_order, col_order]

# --- pretty labels for paper ---
pretty_labels = {
    'BioGRID': 'BioGRID',
    'GWAS': 'GWAS Catalog',
    'MSigDB': 'MSigDB',
    'OMIM': 'OMIM',
    'OpenTargets': 'Open Targets',
    'gencode_gene_node': 'GENCODE',
    'miRDB': 'miRDB',
    'reactome': 'Reactome',
    'clinvar_gene_node': 'ClinVar (gene)',
    'uniprot_base': 'UniProt',
    'variant_annotations': 'dbSNP',
    'prot': 'Protein structure', 
    'alphamissense': 'AlphaMissense',
    'chembl': 'ChEMBL',
    'aa_seq': 'AA sequence',
    'Summarize_bioGRID_GO': 'BioGRID summ.',
    'drug_central': 'DrugCentral',
    'humanbase_functions': 'HumanBase',
    'intact_viral': 'IntAct Viral',
    
}

crosstab_pretty = crosstab_prop.copy()
crosstab_pretty.index = [pretty_labels.get(t, t) for t in crosstab_prop.index]
crosstab_pretty.columns = [pretty_labels.get(t, t) for t in crosstab_prop.columns]

# --- styling & plotting (fixed tick_params) ---
sns.set_theme(style="white")

fig, ax = plt.subplots(figsize=(9, 5), dpi=300)

heatmap = sns.heatmap(
    crosstab_pretty,
    ax=ax,
    annot=True,
    fmt=".2f",
    cmap="Blues",
    annot_kws={"size": 8.5},
    vmin=0,
    vmax=1,
    cbar_kws={"shrink": 0.9, "label": ""},
)

ax.set_title("", fontsize=14, pad=10)
ax.set_xlabel("Tool used", fontsize=12)
ax.set_ylabel("Tool needed", fontsize=12)

# tick label sizes
ax.tick_params(axis='x', labelrotation=0, labelsize=10)
ax.tick_params(axis='y', labelrotation=0, labelsize=10)

# remove colorbar
cbar = heatmap.collections[0].colorbar
cbar.remove()

# now set rotation + alignment on Text objects (this is where 'ha' belongs)
for label in ax.get_xticklabels():
    label.set_rotation(90)
    label.set_horizontalalignment('right')

plt.tight_layout()

# plt.savefig("tool_usage_heatmap.png", dpi=300)

plt.show()

# %%
import numpy as np
# --- ORDER rows / columns explicitly for the figure ---
row_order = needed_tools

# move "'GWAS+AlphaMissense'" to the bottom
if 'Annot+AlphaMissense' in row_order:
    row_order.remove('Annot+AlphaMissense')
    row_order.append('Annot+AlphaMissense')

# align beginning of col_order with row_order  (this defines the "square")
square_cols = [i for i in row_order if i in used_tools]

# ensure AlphaMissense is part of the square (if it exists at all)
if 'alphamissense' in used_tools and 'alphamissense' not in square_cols:
    square_cols.append('alphamissense')

# build final col_order: square first, then remaining tools
col_order = list(square_cols)
for t in used_tools:
    if t not in col_order:
        col_order.append(t)

# filter to actual table
row_order = [t for t in row_order if t in crosstab_prop.index]
col_order = [t for t in col_order if t in crosstab_prop.columns]

# --- insert a small visual gap after the square block (before the "remainder") ---
gap_key = "__GAP__"
insert_pos = len([t for t in square_cols if t in crosstab_prop.columns])  # robust to filtering

crosstab_prop[gap_key] = np.nan  # blank separator column
col_order = col_order[:insert_pos] + [gap_key] + col_order[insert_pos:]

crosstab_prop = crosstab_prop.loc[row_order, col_order]

# --- pretty labels for paper ---
pretty_labels = {
    'BioGRID': 'BioGRID',
    'GWAS': 'GWAS Catalog',
    'MSigDB': 'MSigDB',
    'OMIM': 'OMIM',
    'OpenTargets': 'Open Targets',
    'gencode_gene_node': 'GENCODE',
    'miRDB': 'miRDB',
    'reactome': 'Reactome',
    'clinvar_gene_node': 'ClinVar (gene)',
    'uniprot_base': 'UniProt',
    'variant_annotations': 'dbSNP',
    'prot': 'Protein structure',
    'alphamissense': 'AlphaMissense',
    'chembl': 'ChEMBL',
    'aa_seq': 'AA sequence',
    'Summarize_bioGRID_GO': 'BioGRID summ.',
    'drug_central': 'DrugCentral',
    'humanbase_functions': 'HumanBase',
    'intact_viral': 'IntAct Viral',
    gap_key: '',  # blank label for the gap column
}

crosstab_pretty = crosstab_prop.copy()
crosstab_pretty.index = [pretty_labels.get(t, t) for t in crosstab_prop.index]
crosstab_pretty.columns = [pretty_labels.get(t, t) for t in crosstab_prop.columns]

# --- styling & plotting (fixed tick_params) ---
sns.set_theme(style="white")

fig, ax = plt.subplots(figsize=(9, 5), dpi=300)

heatmap = sns.heatmap(
    crosstab_pretty,
    ax=ax,
    annot=True,
    fmt=".2f",
    cmap="Blues",
    annot_kws={"size": 8.9},
    vmin=0,
    vmax=1,
    cbar_kws={"shrink": 0.9, "label": ""},
    mask=crosstab_pretty.isna(),  # keeps the gap column blank + avoids "nan" annotations
)

ax.set_title("", fontsize=14, pad=10)
ax.set_xlabel("", fontsize=12)
ax.set_ylabel("", fontsize=12)

# tick label sizes
ax.tick_params(axis='x', labelrotation=0, labelsize=13.4)
ax.tick_params(axis='y', labelrotation=0, labelsize=13.4)

# remove colorbar
cbar = heatmap.collections[0].colorbar
cbar.remove()

# now set rotation + alignment on Text objects (this is where 'ha' belongs)
for label in ax.get_xticklabels():
    label.set_rotation(90)
    label.set_horizontalalignment('right')

plt.tight_layout()

# plt.savefig("tool_usage_heatmap.png", dpi=300)

plt.show()

# %%
