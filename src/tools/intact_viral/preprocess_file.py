# %%
# load in the IntAct Viral interaction file and preprocess it for quick querying
import pandas as pd
file = '/Users/sokolova/Documents/research/alvessa_agent/local_dbs/intact_viral.txt'
intact_df = pd.read_csv(file, sep='\t')
intact_df.head()
# %%
intact_df['Alias(es) interactor A'].iloc[0]

# %%
# find "uniprotkb:NAME(gene name)""
import re
pattern = r'uniprotkb:([^(\s]+)\(([^)]+)\)'
matches = re.findall(pattern, intact_df['Alias(es) interactor A'].iloc[0])
matches

# now add columns with gene names for interactor A and B, use gene name only (not synonym)
def extract_gene_names(alias_str: str) -> list[str]:
    matches = re.findall(pattern, alias_str or "")
    # only keep the gene name
    matches = [i for i in matches if i[1]=='gene name']  # filter out entries without gene name
    gene_names = [match[0].upper() for match in matches]
    return gene_names
intact_df['Gene Names A'] = intact_df['Alias(es) interactor A'].apply(extract_gene_names)
intact_df['Gene Names B'] = intact_df['Alias(es) interactor B'].apply(extract_gene_names)
# save the dataframe with new columns
# intact_df.to_csv('/Users/sokolova/Documents/research/alvessa_agent/local_dbs/intact_viral_with_genes.txt', sep='\t', index=False)
intact_df.head()
# %%
# check how many A columns haeve AСE2 
ace2_count = intact_df['Gene Names A'].apply(lambda x: 'ACE2' in x).sum()
ace2_count
# %%
# now convert to a dictionary gene A: [gene B list]
gene_interactions = {}
for idx, row in intact_df.iterrows():
    aliases_a = row['Alias(es) interactor A']
    aliases_b = row['Alias(es) interactor B']
    
    matches_a = re.findall(pattern, aliases_a or "")
    matches_b = re.findall(pattern, aliases_b or "")
    
    gene_names_a = [match[0].upper() for match in matches_a]
    gene_names_b = [match[0].upper() for match in matches_b]
    
    for gene_a in gene_names_a:
        if gene_a not in gene_interactions:
            gene_interactions[gene_a] = set()
        for gene_b in gene_names_b:
            gene_interactions[gene_a].add(gene_b)

# %%
gene_interactions
# %%
