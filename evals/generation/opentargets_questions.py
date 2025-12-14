import pandas as pd
import random
import os
import sys
from pathlib import Path
import json

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


ROOT = Path(__file__).resolve().parents[2]
LOCAL_DBS = ROOT / "local_dbs"
OPEN_TARGETS = LOCAL_DBS / "open_targets"
ESSENTIALITY_DATA = OPEN_TARGETS / "final_target_essentiality"
EXPRESSION_DATA = OPEN_TARGETS / "final_expression"

OUTPUT = ROOT / 'benchmarks_generation/benchmark_questions/OpenTargets'

essentiality_df = read_all_parquet_in_folder(ESSENTIALITY_DATA)
expression_df = read_all_parquet_in_folder(EXPRESSION_DATA)

essential_genes = essentiality_df[essentiality_df['geneEssentiality'].apply(lambda x: x[0]['isEssential'])]['target_symbol'].tolist()
nonessential_genes = essentiality_df[~essentiality_df['geneEssentiality'].apply(lambda x: x[0]['isEssential'])]['target_symbol'].tolist()

all_genes = essential_genes + nonessential_genes

# Which of the following is a core essential gene?
def opentargets_set1(count):
    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:

        random_gene_idx = random.sample(range(len(essential_genes)), 1)[0]
        correct_choice = essential_genes[random_gene_idx]

        false_choice_idxs = random.sample(range(len(nonessential_genes)), 3)
        false_choices = [list(nonessential_genes)[i] for i in false_choice_idxs]

        all_choices = false_choices + [correct_choice]
        random.shuffle(all_choices)

        question = f'Which of the following is a core essential gene? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'

        correct_choice_idx_in_question = all_choices.index(correct_choice)
        all_letters = ['A', 'B', 'C', 'D']
        correct_letter = all_letters[correct_choice_idx_in_question]
        final_correct_answer = f'{correct_letter}'
        new_questions.append({'question': question, 'answer': final_correct_answer, 'tool': 'OpenTargets'})

        chosen_gene_idxs.append(random_gene_idx)

    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    os.makedirs(OUTPUT, exist_ok=True)
    output_df.to_csv(f'{OUTPUT}/set1.csv', index = False)
    return

# Which of the following tissues is gene X most highly expressed in?
def opentargets_set2(count):
    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:
        available_values = list(set(range(len(all_genes))) - set(chosen_gene_idxs))
        random_gene_idx = random.sample(available_values, 1)[0]
        curr_gene = all_genes[random_gene_idx]

        highly_expressed = set()
        low_expressed = set()

        for tissues_list in expression_df[expression_df['target_symbol'] == curr_gene]['tissues']:
            for tissue_dict in tissues_list:
                label = tissue_dict.get('label', None)
                zscore = tissue_dict.get('rna', {}).get('zscore', None)
                if label is not None and zscore is not None:
                    if zscore > .674:
                        highly_expressed.add(label)
                    elif zscore <= -1:
                        low_expressed.add(label)

        highly_expressed = list(highly_expressed)
        low_expressed = list(low_expressed)

        if len(highly_expressed)==0 or len(low_expressed)<3:
            continue

        random_tissue_idx = random.sample(range(len(highly_expressed)), 1)[0]
        correct_choice = list(highly_expressed)[random_tissue_idx]

        false_choice_idxs = random.sample(range(len(low_expressed)), 3)
        false_choices = [list(low_expressed)[i] for i in false_choice_idxs]

        all_choices = false_choices + [correct_choice]
        random.shuffle(all_choices)

        question = f'Which of the following tissues or cell types is gene {curr_gene} most highly expressed in? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'

        correct_choice_idx_in_question = all_choices.index(correct_choice)
        all_letters = ['A', 'B', 'C', 'D']
        correct_letter = all_letters[correct_choice_idx_in_question]
        final_correct_answer = f'{correct_letter}'
        new_questions.append({'question': question, 'answer': final_correct_answer, 'tool': 'OpenTargets'})

        chosen_gene_idxs.append(random_gene_idx)

    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    os.makedirs(OUTPUT, exist_ok=True)
    output_df.to_csv(f'{OUTPUT}/set2.csv', index = False)
    return


if __name__ == "__main__":
    num_questions = 20
    # print("=== Generating OpenTargets Question Set 1 ===")
    # opentargets_set1(num_questions)
    print("=== Generating OpenTargets Question Set 2 ===")
    opentargets_set2(num_questions)
