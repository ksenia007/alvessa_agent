import pandas as pd
import random 
import os
import requests
import csv
from itertools import combinations

mirdb_df = pd.read_csv(
        '../local_dbs/miRDB_v6.0_prediction_result.txt',
        sep='\t',
        header=None,
        names=['miRNAID', 'geneID', 'confidence']
    )
mirna_list = [miRNAID[4:] for miRNAID in set(mirdb_df['miRNAID'].values)]
all_mirna_count = len(mirna_list)

def _convert_from_miRBASE(symbol):
    components = symbol.split("-")

    output = f'MIR{components[1]}'
    if len(components)==3:
        output+=f'_{components[2][0]}P'

    return output

def _refseq_to_symbol(refseq_id):

    try:
        r = requests.get(url = f"http://mygene.info/v3/query?q={refseq_id}&fields=symbol")
        r.raise_for_status()
        hits = r.json().get("hits", [])
        symbol = [hit.get('symbol') for hit in hits if 'symbol' in hit][0]
        return symbol
    except Exception as e:
        print(e)

    return None

# Which of the following is a predicted gene target of the miRNA X?
def mirdb_set1(count):
    chosen_mirna_idxs = []

    output_df = pd.DataFrame()
    new_questions = []
    all_genes_list = list(set(mirdb_df['geneID'].values))
    all_genes_list_len = len(all_genes_list)

    while len(new_questions)<count:
        available_values = list(set(range(all_mirna_count)) - set(chosen_mirna_idxs))
        random_mirna_idx = random.sample(available_values, 1)[0]
        curr_mirna = mirna_list[random_mirna_idx]

        try:
            match = mirdb_df[mirdb_df['miRNAID'].str.contains(curr_mirna)]
            if len(match)>0:
                gene_targets = list(set(match['geneID'].values))
                # print(gene_targets)

                random_pathway_idxs = random.sample(range(len(gene_targets)), 1)[0]

                correct_choice = _refseq_to_symbol(gene_targets[random_pathway_idxs])
                print(correct_choice)

                gene_pathways_idxs = []
                for pathway in gene_targets:
                    gene_pathways_idxs.append(all_genes_list.index(pathway))

                available_values = list(set(range(all_genes_list_len)) - set(gene_pathways_idxs))
                false_choice_idxs = random.sample(available_values, 3)

                false_choices = [_refseq_to_symbol(list(all_genes_list)[i]) for i in false_choice_idxs]

                all_choices = false_choices + [correct_choice]
                random.shuffle(all_choices)
                if any(x is None for x in all_choices):
                    continue
                
                curr_mirna = _convert_from_miRBASE(curr_mirna)
                question = f'Which of the following is a predicted gene target of the miRNA {curr_mirna}? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'

                correct_choice_idx_in_question = all_choices.index(correct_choice)
                all_letters = ['A', 'B', 'C', 'D']
                correct_letter = all_letters[correct_choice_idx_in_question]
                final_correct_answer = f'{correct_letter}'
                new_questions.append({'question': question, 'answer': final_correct_answer, 'tool': 'miRDB'})

                chosen_mirna_idxs.append(random_mirna_idx)
        except Exception as e:
            continue

    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    output_path = 'benchmark_questions/miRDB'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set1.csv', index = False)
    return

# Which of the following is a predicted gene target of both miRNA X and miRNA Y?
def mirdb_set2(count):

    output_df = pd.DataFrame()
    new_questions = []

    ### Code to create overlapping miRNAID_pairs.csv
    # pairs = set()
    # grouped = mirdb_df.groupby('geneID')['miRNAID'].apply(list)

    # for miRNAs in grouped:
    #     miRNA_str = [m[4:] for m in miRNAs]
    #     for pair in combinations(sorted(set(miRNA_str)), 2):
    #         pairs.add(pair)

    # pairs_list = list(pairs)

    # print(len(pairs_list))

    # with open("../local_dbs/overlapping_gene_targets_miRNAID_pairs.csv", "w", newline="") as f:
    #     writer = csv.writer(f)
    #     writer.writerow(["miRNAID1", "miRNAID2"])
    #     writer.writerows(pairs_list)

    pairs_list = []
    with open("../local_dbs/overlapping_gene_targets_miRNAID_pairs.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)  
        for row in reader:
            pairs_list.append((row[0], row[1]))


    chosen_pair_idxs = []
    all_genes_list = list(set(mirdb_df['geneID'].values))
    all_genes_list_len = len(all_genes_list)

    while len(new_questions)<count:
        available_values = list(set(range(len(pairs_list))) - set(chosen_pair_idxs))
        random_pair_idx = random.sample(available_values, 1)[0]
        curr_mirna_pair = pairs_list[random_pair_idx]

        all_targets_union = []
        should_skip = False
        for curr_mirna in curr_mirna_pair:
            try:
                match = mirdb_df[mirdb_df['miRNAID'].str.contains(curr_mirna)]
                if len(match)>0:
                    gene_targets = list(set(match['geneID'].values))
                    # print(gene_targets)
                    all_targets_union.append(gene_targets)
                else:
                    should_skip = True
                    break
            except:
                should_skip = True
                break

            if should_skip:
                continue

            # print(all_targets_union)

            try:
                mirna_x = _convert_from_miRBASE(curr_mirna_pair[0])
                mirna_y = _convert_from_miRBASE(curr_mirna_pair[1])
                target_intersection_list = list(set(all_targets_union[0]) & set(all_targets_union[1]))
            except:
                continue

            if len(target_intersection_list)==0:
                continue

            random_pathway_idxs = random.sample(range(len(target_intersection_list)), 1)[0]
            correct_choice = _refseq_to_symbol(target_intersection_list[random_pathway_idxs])

            gene_target_idxs = []
            for target in target_intersection_list:
                gene_target_idxs.append(all_genes_list.index(target))

            available_values = list(set(range(all_genes_list_len)) - set(gene_target_idxs))
            false_choice_idxs = random.sample(available_values, 3)

            false_choices = [_refseq_to_symbol(list(all_genes_list)[i]) for i in false_choice_idxs]

            all_choices = false_choices + [correct_choice]
            random.shuffle(all_choices)
            if any(x is None for x in all_choices):
                continue
            
            curr_mirna = _convert_from_miRBASE(curr_mirna)
            question = f'Which of the following is a predicted gene target of both miRNA {mirna_x} and miRNA {mirna_y}? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'

            correct_choice_idx_in_question = all_choices.index(correct_choice)
            all_letters = ['A', 'B', 'C', 'D']
            correct_letter = all_letters[correct_choice_idx_in_question]
            final_correct_answer = f'{correct_letter}'
            new_questions.append({'question': question, 'answer': final_correct_answer, 'tool': 'miRDB'})

            chosen_pair_idxs.append(random_pair_idx)

    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    output_path = 'benchmark_questions/miRDB'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set2.csv', index = False)
    return

if __name__ == "__main__":
    num_questions = 20
    # print("=== Generating MiRDB Question Set 1 ===")
    # mirdb_set1(num_questions)
    print("=== Generating MiRDB Question Set 2 ===")
    mirdb_set2(num_questions)

