import pandas as pd
import os

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
