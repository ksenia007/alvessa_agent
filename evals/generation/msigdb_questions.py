import pandas as pd
import random
import os
import sys
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]
LOCAL_DBS = ROOT / "local_dbs"
MSIGDB = LOCAL_DBS / "gene_centered_msigdb.v2025.1.Hs.json.txt"

OUTPUT = ROOT / 'benchmarks_generation/benchmark_questions/MSigDB_test'

with open(MSIGDB, "r") as file:
    msigdata = json.load(file)

all_genes = list(set(msigdata.keys()))

gene_set_names_c1 = set()
gene_set_names_c6 = set()

for gene, gene_sets in msigdata.items():
    for gene_set_name, collection in gene_sets.items():
        if "C1" in collection:
            gene_set_names_c1.add(gene_set_name)
        elif "C6" in collection:
            gene_set_names_c6.add(gene_set_name)

gene_set_names_c1 = list(gene_set_names_c1)
gene_set_names_c6 = list(gene_set_names_c6)


# Which of the following human chromosome cytogenetic bands is gene X located at?
def msigdb_set1(count):
    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:
        available_values = list(set(range(len(all_genes))) - set(chosen_gene_idxs))
        random_gene_idx = random.sample(available_values, 1)[0]
        curr_gene = all_genes[random_gene_idx]

        gene_annotations = set()

        for gene_set, collection in msigdata[curr_gene].items():
            if collection == 'C1':
                gene_annotations.add(gene_set)
        gene_annotations = list(gene_annotations)

        if len(gene_annotations) == 0:
            continue

        random_annotation_idx = random.sample(range(len(gene_annotations)), 1)[0]

        correct_choice = list(gene_annotations)[random_annotation_idx]

        gene_annotation_idxs = []
        for annotation in gene_annotations:
            gene_annotation_idxs.append(gene_set_names_c1.index(annotation))

        available_values = list(set(range(len(gene_set_names_c1))) - set(gene_annotation_idxs))
        false_choice_idxs = random.sample(available_values, 3)

        false_choices = [list(gene_set_names_c1)[i] for i in false_choice_idxs]

        all_choices = false_choices + [correct_choice]
        random.shuffle(all_choices)

        question = f'Which of the following human chromosome cytogenetic bands is gene {curr_gene} located at? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'

        correct_choice_idx_in_question = all_choices.index(correct_choice)
        all_letters = ['A', 'B', 'C', 'D']
        correct_letter = all_letters[correct_choice_idx_in_question]
        final_correct_answer = f'{correct_letter}'
        new_questions.append({'question': question, 'answer': final_correct_answer, 'tool': 'MSigDB'})

        chosen_gene_idxs.append(random_gene_idx)

    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    os.makedirs(OUTPUT, exist_ok=True)
    output_df.to_csv(f'{OUTPUT}/set1.csv', index = False)
    return

# Which of the following oncogenic signatures is gene X associated with?
def msigdb_set2(count):
    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:
        available_values = list(set(range(len(all_genes))) - set(chosen_gene_idxs))
        random_gene_idx = random.sample(available_values, 1)[0]
        curr_gene = all_genes[random_gene_idx]

        gene_annotations = set()

        for gene_set, collection in msigdata[curr_gene].items():
            if collection == 'C6':
                gene_annotations.add(gene_set)
        gene_annotations = list(gene_annotations)

        if len(gene_annotations) == 0:
            continue

        random_annotation_idx = random.sample(range(len(gene_annotations)), 1)[0]

        correct_choice = list(gene_annotations)[random_annotation_idx]

        gene_annotation_idxs = []
        for annotation in gene_annotations:
            gene_annotation_idxs.append(gene_set_names_c6.index(annotation))

        available_values = list(set(range(len(gene_set_names_c6))) - set(gene_annotation_idxs))
        false_choice_idxs = random.sample(available_values, 3)

        false_choices = [list(gene_set_names_c6)[i] for i in false_choice_idxs]

        all_choices = false_choices + [correct_choice]
        random.shuffle(all_choices)

        question = f'Which of the following oncogenic signatures is gene {curr_gene} associated with? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'

        correct_choice_idx_in_question = all_choices.index(correct_choice)
        all_letters = ['A', 'B', 'C', 'D']
        correct_letter = all_letters[correct_choice_idx_in_question]
        final_correct_answer = f'{correct_letter}'
        new_questions.append({'question': question, 'answer': final_correct_answer, 'tool': 'MSigDB'})

        chosen_gene_idxs.append(random_gene_idx)

    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    os.makedirs(OUTPUT, exist_ok=True)
    output_df.to_csv(f'{OUTPUT}/set2.csv', index = False)
    return


if __name__ == "__main__":
    num_questions = 20
    print("=== Generating MSigDB Question Set 1 ===")
    msigdb_set1(num_questions)
    print("=== Generating MSigDB Question Set 2 ===")
    msigdb_set2(num_questions)
