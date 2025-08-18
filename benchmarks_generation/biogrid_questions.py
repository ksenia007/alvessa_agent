import pandas as pd
import random
import os
import sys
sys.path.append('../')
from tool_biogrid import _fetch_predictions_BioGRID

with open('../local_dbs/gene_names_list.txt', 'r') as f:
    gene_list = [line.strip() for line in f]
all_genes_count = len(gene_list)

def biogrid_set1(count):
    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(chosen_gene_idxs)<count:
        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))
        random_gene_idx = random.sample(available_values, 1)[0]
        curr_gene = gene_list[random_gene_idx]
        try:
            interactors = _fetch_predictions_BioGRID(curr_gene)

            max_interactor_idx = min(100, len(interactors))

            random_interactor_idxs = random.sample(range(max_interactor_idx), 3)

            false_choices = [list(interactors)[i] for i in random_interactor_idxs]

            interactors_idxs_in_gene_list = []
            for interactor in list(interactors):
                try:
                    interactors_idxs_in_gene_list.append(gene_list.index(interactor))
                except:
                    continue
            interactors_idxs_in_gene_list.append(random_gene_idx)

            available_values = list(set(range(all_genes_count)) - set(interactors_idxs_in_gene_list))
            correct_choice_idx = random.sample(available_values, 1)[0]

            correct_choice = gene_list[correct_choice_idx]

            all_choices = false_choices + [correct_choice]
            random.shuffle(all_choices)

            question = f'Which of the following genes is NOT interacting with gene {curr_gene}? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'

            correct_choice_idx_in_question = all_choices.index(correct_choice)
            all_letters = ['A', 'B', 'C', 'D']
            correct_letter = all_letters[correct_choice_idx_in_question]
            final_correct_answer = f'{correct_letter}'
            new_questions.append({'question': question, 'answer': final_correct_answer, 'tool': 'BioGRID'})

            chosen_gene_idxs.append(random_gene_idx)
        except Exception as e:
            # continue
            print(curr_gene, e)

    
    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    output_path = 'benchmark_questions/BioGRID'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set1.csv', index = False)

def biogrid_set2(count):
    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(chosen_gene_idxs)<count:
        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))
        random_gene_idx = random.sample(available_values, 1)[0]
        curr_gene = gene_list[random_gene_idx]
        try:
            interactors = _fetch_predictions_BioGRID(curr_gene)

            max_interactor_idx = min(100, len(interactors))
            random_interactor_idxs = random.sample(range(max_interactor_idx), 1)[0]
            correct_choice = list(interactors)[random_interactor_idxs]

            interactors_idxs_in_gene_list = []
            for interactor in list(interactors):
                try:
                    interactors_idxs_in_gene_list.append(gene_list.index(interactor))
                except:
                    continue
            interactors_idxs_in_gene_list.append(random_gene_idx)

            available_values = list(set(range(all_genes_count)) - set(interactors_idxs_in_gene_list))
            false_choice_idxs = random.sample(available_values, 3)

            false_choices = [list(gene_list)[i] for i in false_choice_idxs]

            all_choices = false_choices + [correct_choice]
            random.shuffle(all_choices)

            question = f'Which of the following genes is interacting with gene {curr_gene}? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'

            correct_choice_idx_in_question = all_choices.index(correct_choice)
            all_letters = ['A', 'B', 'C', 'D']
            correct_letter = all_letters[correct_choice_idx_in_question]
            final_correct_answer = f'{correct_letter}'
            new_questions.append({'question': question, 'answer': final_correct_answer, 'tool': 'BioGRID'})

            chosen_gene_idxs.append(random_gene_idx)
        except Exception as e:
            # continue
            print(curr_gene, e)

    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)
    
    output_path = 'benchmark_questions/BioGRID'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set2.csv', index = False)


if __name__ == "__main__":
    num_questions = 20
    print("=== Generating BioGRID Question Set 1 ===")
    biogrid_set1(num_questions)
    print("=== Generating BioGRID Question Set 2 ===")
    biogrid_set2(num_questions)
    
    

    