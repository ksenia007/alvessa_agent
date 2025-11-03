import pandas as pd
import os
from pathlib import Path
import pprint
import pickle

REPO_ROOT = Path(__file__).resolve().parents[3]
LOCAL_DBS = REPO_ROOT / "local_dbs"
OPEN_TARGETS = LOCAL_DBS / "open_targets"
TARGET_DISEASE_DATA = OPEN_TARGETS / "final_association_overall_direct"
EXPRESSION_DATA = OPEN_TARGETS / "final_expression"
ESSENTIALITY_DATA = OPEN_TARGETS / "final_target_essentiality"
CONSTRAINT_DATA = OPEN_TARGETS / "target"

def read_all_parquet_in_folder(folder_path):
    all_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.parquet')]

    if not all_files:
        print(f"No Parquet files found in '{folder_path}'")
        return pd.DataFrame()

    try:
        df = pd.read_parquet(all_files[0])
    except Exception as e:
        print(f"Error reading initial Parquet file '{all_files[0]}': {e}")
        return pd.DataFrame()

    for file_path in all_files[1:]:
        try:
            temp_df = pd.read_parquet(file_path)
            df = pd.concat([df, temp_df], ignore_index=True)
        except Exception as e:
            print(f"Error reading and concatenating Parquet file '{file_path}': {e}")
            continue 

    return df

def add_disease_name_to_target_disease_data():
    associations = read_all_parquet_in_folder('../../../local_dbs/open_targets/association_overall_direct')
    diseases = read_all_parquet_in_folder('../../../local_dbs/open_targets/disease')

    merged = pd.merge(associations, diseases[['id', 'name']], left_on='diseaseId', right_on='id', how='left')
    merged = merged.loc[:,['diseaseId', 'name', 'targetId', 'score', 'evidenceCount']]
    merged.rename(columns={'name': 'disease_name'}, inplace=True)

    return merged
    
def add_gene_symbol_to_all_data():
    associations = add_disease_name_to_target_disease_data()
    expression = read_all_parquet_in_folder('../../../local_dbs/open_targets/expression')
    essentiality = read_all_parquet_in_folder('../../../local_dbs/open_targets/target_essentiality')
    targets = read_all_parquet_in_folder('../../../local_dbs/open_targets/target')

    target_disease_merged = pd.merge(associations, targets[['id', 'approvedSymbol']], left_on='targetId', right_on='id', how='left')
    target_disease_merged = target_disease_merged.loc[:,['diseaseId', 'disease_name', 'targetId', 'approvedSymbol', 'score', 'evidenceCount']]
    target_disease_merged.rename(columns={'approvedSymbol': 'target_symbol'}, inplace=True)

    target_disease_merged.to_parquet('../../../local_dbs/open_targets/final_association_overall_direct/target_disease_direct_associations.parquet')

    expression_merged = pd.merge(expression, targets[['id', 'approvedSymbol']], on='id', how='left')
    expression_merged = expression_merged.loc[:,['id', 'approvedSymbol', 'tissues']]
    expression_merged.rename(columns={'approvedSymbol': 'target_symbol'}, inplace=True)

    expression_merged.to_parquet('../../../local_dbs/open_targets/final_expression/expression.parquet')

    essentiality_merged = pd.merge(essentiality, targets[['id', 'approvedSymbol']], on='id', how='left')
    essentiality_merged = essentiality_merged.loc[:,['id', 'approvedSymbol', 'geneEssentiality']]
    essentiality_merged.rename(columns={'approvedSymbol': 'target_symbol'}, inplace=True)

    essentiality_merged.to_parquet('../../../local_dbs/open_targets/final_target_essentiality/target_essentiality.parquet')

def convert_final_dfs_to_dicts():

    target_disease_df = read_all_parquet_in_folder(TARGET_DISEASE_DATA)

    gene_disease_dict = (target_disease_df.groupby('target_symbol')['disease_name'].unique().apply(list).to_dict())

    with open('../../../local_dbs/open_targets/final_association_overall_direct/target_disease_direct_associations.pkl', 'wb') as file:
        pickle.dump(gene_disease_dict, file)


    expression_df = read_all_parquet_in_folder(EXPRESSION_DATA)

    gene_tissue_zscore = {}
    for symbol, group_df in expression_df.groupby('target_symbol'):
        tissue_zscore = {}
        for tissues_list in group_df['tissues']:
            for tissue_dict in tissues_list:
                label = tissue_dict.get('label', None)
                zscore = tissue_dict.get('rna', {}).get('zscore', None)
                if label is not None and zscore is not None:
                    tissue_zscore[label] = zscore
        gene_tissue_zscore[symbol] = tissue_zscore

    with open('../../../local_dbs/open_targets/final_expression/expression.pkl', 'wb') as file:
        pickle.dump(gene_tissue_zscore, file)


    essentiality_df = read_all_parquet_in_folder(ESSENTIALITY_DATA)

    gene_essentiality = {}
    for _, row in essentiality_df.iterrows():
        symbol = row['target_symbol']
        is_essential = row['geneEssentiality'][0]['isEssential']
        gene_essentiality[symbol] = is_essential

    with open('../../../local_dbs/open_targets/final_target_essentiality/target_essentiality.pkl', 'wb') as file:
        pickle.dump(gene_essentiality, file)

    
    constraint_df = read_all_parquet_in_folder(CONSTRAINT_DATA)

    gene_constraint = {}
    for _, row in constraint_df.iterrows():
        symbol = row['approvedSymbol']
        constraints = {}
        all_constraints = row['constraint']
        if all_constraints is not None:
            for constraint_row in all_constraints:
                constraint_type_name = constraint_row['constraintType']
                score = constraint_row['score']
                constraints[constraint_type_name] = score
        gene_constraint[symbol] = constraints

    with open('../../../local_dbs/open_targets/final_constraint/constraint.pkl', 'wb') as file:
        pickle.dump(gene_constraint, file)

