import pandas as pd
import random
import os
import sys
sys.path.append('../')
from tool_biogrid import _fetch_predictions_BioGRID

with open('../local_dbs/gene_names_list.txt', 'r') as f:
    gene_list = [line.strip() for line in f]
all_genes_count = len(gene_list)

## Which of the following genes is interacting with gene X?
def biogrid_set1(count):
    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:
        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))
        random_gene_idx = random.sample(available_values, 1)[0]
        curr_gene = gene_list[random_gene_idx]
        try:
            all_interactors, human_interactors, nonhuman_interactors = _fetch_predictions_BioGRID(curr_gene)

            all_human_interactors = []
            for value_list in human_interactors.values():
                all_human_interactors.extend(value_list)

            # print("all_genes: ", all_interactors)
            # print("human_genes: ", all_human_interactors)  

            max_interactor_idx = len(all_human_interactors)
            random_interactor_idxs = random.sample(range(max_interactor_idx), 1)[0]
            correct_choice = list(all_human_interactors)[random_interactor_idxs]

            interactors_idxs_in_gene_list = []
            for interactor in list(all_interactors):
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
    output_df.to_csv(f'{output_path}/set1.csv', index = False)

## Which of the following genes interacts both with gene X and gene Y?
def biogrid_set2(count):
    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:
        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))

        chosen_gene_idxs = []
        chosen_genes = []
        all_interactors_union = []
        human_interactors_union = []

        while len(chosen_genes)<2:
            random_gene_idx = random.sample(available_values, 1)[0]
            curr_gene = gene_list[random_gene_idx]
            try:
                all_interactors, human_interactors, nonhuman_interactors = _fetch_predictions_BioGRID(curr_gene)
                all_interactors_union.extend(all_interactors)

                all_human_interactors = []
                for value_list in human_interactors.values():
                    all_human_interactors.extend(value_list)

                # print("all_genes: ", all_interactors_union)
                # print("human_genes: ", all_human_interactors)        

                human_interactors_union.extend(all_human_interactors)
                chosen_gene_idxs.append(random_gene_idx)
                chosen_genes.append(curr_gene)
            except Exception as e:
                # continue
                print(curr_gene, e)
            
        max_interactor_idx = len(human_interactors_union)
        random_interactor_idxs = random.sample(range(max_interactor_idx), 1)[0]
        correct_choice = list(human_interactors_union)[random_interactor_idxs]

        interactors_idxs_in_gene_list = []
        for interactor in list(all_interactors_union):
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

        question = f'Which of the following genes interacts both with gene {chosen_genes[0]} and gene {chosen_genes[1]}? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'

        correct_choice_idx_in_question = all_choices.index(correct_choice)
        all_letters = ['A', 'B', 'C', 'D']
        correct_letter = all_letters[correct_choice_idx_in_question]
        final_correct_answer = f'{correct_letter}'
        new_questions.append({'question': question, 'answer': final_correct_answer, 'tool': 'BioGRID'})

        chosen_gene_idxs.append(random_gene_idx)
        

    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)
    
    output_path = 'benchmark_questions/BioGRID'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set2.csv', index = False)

## Which of the following genes interacts with gene X through gene Y?
def biogrid_set3(count):
    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:
        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))
        random_gene_idx = random.sample(available_values, 1)[0]
        gene_x = gene_list[random_gene_idx]
        try:
            all_gene_x_interactors, human_gene_x_interactors, nonhuman_gene_x_interactors = _fetch_predictions_BioGRID(gene_x)
            # print("all_x: ", all_gene_x_interactors)
            # print("human_x: ", human_gene_x_interactors)
            all_human_gene_x_interactors = []
            for value_list in human_gene_x_interactors.values():
                all_human_gene_x_interactors.extend(value_list)

            max_interactor_idx = len(all_human_gene_x_interactors)
            gene_y_idx = random.sample(range(max_interactor_idx), 1)[0]
            gene_y = list(all_human_gene_x_interactors)[gene_y_idx]
            while True:
                try:
                    all_gene_y_interactors, human_gene_y_interactors, nonhuman_gene_y_interactors = _fetch_predictions_BioGRID(gene_y)
                    break
                except Exception as e:
                    print(gene_y, e)

            all_human_gene_y_interactors = []
            for value_list in human_gene_y_interactors.values():
                all_human_gene_y_interactors.extend(value_list)

            # print("all_y: ", all_gene_y_interactors)
            # print("human_y: ", human_gene_y_interactors)
            in_y_not_in_x = [gene for gene in all_human_gene_y_interactors if gene not in all_gene_x_interactors]

            max_interactor_idx = len(in_y_not_in_x)
            random_interactor_idxs = random.sample(range(max_interactor_idx), 1)[0]
            correct_choice = in_y_not_in_x[random_interactor_idxs]

            interactors_idxs_in_gene_list = []
            for interactor in list(all_gene_y_interactors):
                try:
                    interactors_idxs_in_gene_list.append(gene_list.index(interactor))
                except:
                    continue
            interactors_idxs_in_gene_list.append(random_gene_idx)
            interactors_idxs_in_gene_list.append(gene_y_idx)

            available_values = list(set(range(all_genes_count)) - set(interactors_idxs_in_gene_list))
            false_choice_idxs = random.sample(available_values, 3)

            false_choices = [list(gene_list)[i] for i in false_choice_idxs]

            all_choices = false_choices + [correct_choice]
            random.shuffle(all_choices)

            question = f'Which of the following genes interacts with gene {gene_x} through gene {gene_y}? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'

            correct_choice_idx_in_question = all_choices.index(correct_choice)
            all_letters = ['A', 'B', 'C', 'D']
            correct_letter = all_letters[correct_choice_idx_in_question]
            final_correct_answer = f'{correct_letter}'
            new_questions.append({'question': question, 'answer': final_correct_answer, 'tool': 'BioGRID'})

            chosen_gene_idxs.append(random_gene_idx)
        except Exception as e:
            # continue
            print(gene_x, e)

    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)
    
    output_path = 'benchmark_questions/BioGRID'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set3.csv', index = False)


if __name__ == "__main__":
    num_questions = 20
    print("=== Generating BioGRID Question Set 1 ===")
    biogrid_set1(num_questions)
    print("=== Generating BioGRID Question Set 2 ===")
    biogrid_set2(num_questions)
    print("=== Generating BioGRID Question Set 3 ===")
    biogrid_set3(num_questions)
    

    