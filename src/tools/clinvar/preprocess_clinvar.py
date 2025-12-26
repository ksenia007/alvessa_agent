# %%
import pandas as pd

clinvar_path = "../../../local_dbs/clinvar/variant_summary.txt"
clinvar_gene_cond_path = "../../../local_dbs/clinvar/gene_condition_source_id.txt"
clinvar_allele_path = "../../../local_dbs/clinvar/allele_gene.txt"
# %%
# read in all the files as dataframes
clinvar_df = pd.read_csv(clinvar_path, sep="\t", dtype=str, low_memory=False)
# %%
# filter clinvardf to only 'GRCh38' in 'Assembly'
clinvar_df_grch38 = clinvar_df[clinvar_df['Assembly'] == 'GRCh38'].copy()
# %%
# rename 'RS# (dbSNP)' to 'rsid'
clinvar_df_grch38 = clinvar_df_grch38.rename(columns={'RS# (dbSNP)': 'rsid'})
# remove rsID=-1 
clinvar_df_grch38 = clinvar_df_grch38[clinvar_df_grch38['rsid'] != '-1'].copy()
clinvar_df_grch38
# %%
# read in allele
clinvar_allele_df = pd.read_csv(clinvar_gene_cond_path, sep="\t", dtype=str, low_memory=False)
clinvar_allele_df.head()
# %%
# save allele and clinvar_df_grch38 as parquet files
clinvar_df_grch38.to_parquet("clinvar_genes_variants_grch38.parquet", index=False)
clinvar_allele_df.to_parquet("clinvar_gene_condition_source_id.parquet", index=False)
# %%
# clinvar_df_grch38[['pathogen' in i for i in clinvar_df_grch38.ClinicalSignificance]]
# filter to pathogenic or likely pathogenic
keep = clinvar_df_grch38['ClinicalSignificance'].str.contains('pathogenic', case=False, na=False)
clinvar_pathogenic_df = clinvar_df_grch38[keep].copy()
clinvar_pathogenic_df.to_parquet("clinvar_pathogenic_variants_grch38.parquet", index=False)
# %%
