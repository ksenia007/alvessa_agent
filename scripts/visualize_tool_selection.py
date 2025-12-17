# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
# %%
loc = '/Users/sokolova/Documents/research/alvessa_agent/out/FINAL_GA_20251216-162600_cli/benchmark_summary.csv'
preds = pd.read_csv(loc)

# --- drop rows where either tool_tag or used_tools is NaN ---
df = preds.dropna(subset=['tool_tag', 'used_tools']).copy()
df = df[df['tool_tag']!='query_gwas_extensive,alphamissense']

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

# --- ORDER rows / columns explicitly for the figure ---
row_order = [
    'BioGRID',
    'GWAS',
    'OMIM',
    'OpenTargets',
    'gencode_gene_node',
    'miRDB',
    'MSigDB',
    'reactome',
    'aa_seq',
       'chembl', 'prot', 
]

col_order = [
    'BioGRID',
    'GWAS',
    'OMIM',
    'OpenTargets',
    'gencode_gene_node',
    'miRDB','MSigDB',
    'reactome',
    'aa_seq', 'chembl', 'prot', 
'AllianceOfGenomes','DisGeNet',
       'Summarize_bioGRID_GO', 
       'clinvar_node', 'drug_central',
       'humanbase_functions', 'intact_viral', 
       'uniprot_base', 'uniprot_gwas', 'variant_annotations'
]

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
    'clinvar_node': 'ClinVar',
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
    annot_kws={"size": 8},
    vmin=0,
    vmax=1,
    cbar_kws={"shrink": 0.9, "label": ""},
)

ax.set_title("", fontsize=14, pad=10)
ax.set_xlabel("Tool used", fontsize=12)
ax.set_ylabel("Tool needed", fontsize=12)

# tick label sizes
ax.tick_params(axis='x', labelrotation=0, labelsize=9)
ax.tick_params(axis='y', labelrotation=0, labelsize=9)

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
