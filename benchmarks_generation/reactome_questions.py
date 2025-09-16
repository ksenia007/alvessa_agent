import pandas as pd
import random
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

from src.tools.humanbase.node import _symbol_to_entrez
from itertools import combinations
import csv
import requests

with open('../local_dbs/gene_names_list.txt', 'r') as f:
    gene_list = [line.strip() for line in f]
all_genes_count = len(gene_list)

def entrez_to_symbol(gene):
    try:
        r = requests.get(f'https://mygene.info/v3/query?q=entrezgene:{gene}&fields=symbol')
        r.raise_for_status()
        hits = r.json().get("hits", [])
        for hit in hits:
            try:
                symbol = hit.get("symbol", {})
                if symbol:
                    return symbol
            except Exception as inner_e:
                print(
                    f"Skipping no symbol entry for gene {gene}: {inner_e}"
                )    
    except Exception as e:
        print(e)

    return None

def reactome_set1(count):
    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []
    reactome_df = pd.read_csv('../local_dbs/NCBI2Reactome_All_Levels.txt', names = ['geneID', 'pathwayID', 'url', 'pathway_name', 'evidence_code', 'species'], sep = '\t')
    all_pathways_list = list(set(reactome_df['pathway_name'].values))
    all_pathways_count = len(all_pathways_list)

    while len(new_questions)<count:
        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))
        random_gene_idx = random.sample(available_values, 1)[0]
        curr_gene = gene_list[random_gene_idx]

        try:
            entrez = _symbol_to_entrez(curr_gene)
            if not entrez:
                continue
        
            match = reactome_df[reactome_df['geneID']==int(entrez)]
            if len(match)>0:
                gene_pathways = list(set(match['pathway_name'].values))
                # print(gene_pathways)

                random_pathway_idxs = random.sample(range(len(gene_pathways)), 1)[0]

                correct_choice = list(gene_pathways)[random_pathway_idxs]

                gene_pathways_idxs = []
                for pathway in gene_pathways:
                    gene_pathways_idxs.append(all_pathways_list.index(pathway))

                available_values = list(set(range(all_pathways_count)) - set(gene_pathways_idxs))
                false_choice_idxs = random.sample(available_values, 3)

                false_choices = [list(all_pathways_list)[i] for i in false_choice_idxs]

                all_choices = false_choices + [correct_choice]
                random.shuffle(all_choices)

                question = f'Which of the following pathways is associated with gene {curr_gene}? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'

                correct_choice_idx_in_question = all_choices.index(correct_choice)
                all_letters = ['A', 'B', 'C', 'D']
                correct_letter = all_letters[correct_choice_idx_in_question]
                final_correct_answer = f'{correct_letter}'
                new_questions.append({'question': question, 'answer': final_correct_answer, 'tool': 'reactome'})

                chosen_gene_idxs.append(random_gene_idx)
        except Exception as e:
            continue

    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    output_path = 'benchmark_questions/reactome'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set1.csv', index = False)
    return

def reactome_set2(count):

    output_df = pd.DataFrame()
    new_questions = []
    reactome_df = pd.read_csv('../local_dbs/NCBI2Reactome_All_Levels.txt', names = ['geneID', 'pathwayID', 'url', 'pathway_name', 'evidence_code', 'species'], sep = '\t')
    
    ### Code to create overlapping geneID_pairs.csv
    # pairs = set()
    # grouped = reactome_df.groupby('pathwayID')['geneID'].apply(list)

    # for genes in grouped:
    #     genes_str = [str(g) for g in genes]
    #     for pair in combinations(sorted(set(genes_str)), 2):
    #         pairs.add(pair)

    # pairs_list = list(pairs)

    # print(len(pairs_list))

    # with open("../local_dbs/overlapping_pathways_geneID_pairs.csv", "w", newline="") as f:
    #     writer = csv.writer(f)
    #     writer.writerow(["geneID1", "geneID2"])
    #     writer.writerows(pairs_list)

    pairs_list = []
    with open("../local_dbs/overlapping_pathways_geneID_pairs.csv", "r") as f:
        reader = csv.reader(f)
        next(reader)  
        for row in reader:
            pairs_list.append((row[0], row[1]))
    
    chosen_pair_idxs = []
    all_pathways_list = list(set(reactome_df['pathway_name'].values))
    all_pathways_count = len(all_pathways_list)

    while len(new_questions)<count:
        available_values = list(set(range(len(pairs_list))) - set(chosen_pair_idxs))
        random_pair_idx = random.sample(available_values, 1)[0]
        curr_gene_pair = pairs_list[random_pair_idx]

        all_pathways_union = []

        should_skip = False
        for curr_gene in curr_gene_pair:

            try:
                match = reactome_df[reactome_df['geneID']==int(curr_gene)]
                if len(match)>0:
                    gene_pathways = list(set(match['pathway_name'].values))
                    all_pathways_union.append(gene_pathways)
                    # print(gene_pathways)
                else:
                    should_skip = True
                    break
            except:
                should_skip = True
                break

        if should_skip:
            continue
        
        try:
            gene_x = entrez_to_symbol(curr_gene_pair[0])
            gene_y = entrez_to_symbol(curr_gene_pair[1])
        except:
            continue
        pathway_intersection_list = list(set(all_pathways_union[0]) & set(all_pathways_union[1]))

        if len(pathway_intersection_list)==0:
            continue

        random_pathway_idxs = random.sample(range(len(pathway_intersection_list)), 1)[0]
        correct_choice = pathway_intersection_list[random_pathway_idxs]

        gene_pathways_idxs = []
        for pathway in pathway_intersection_list:
            gene_pathways_idxs.append(all_pathways_list.index(pathway))

        available_values = list(set(range(all_pathways_count)) - set(gene_pathways_idxs))
        false_choice_idxs = random.sample(available_values, 3)

        false_choices = [list(all_pathways_list)[i] for i in false_choice_idxs]

        all_choices = false_choices + [correct_choice]
        random.shuffle(all_choices)

        question = f'Which of the following pathways is associated with both gene {gene_x} and gene {gene_y}? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'

        correct_choice_idx_in_question = all_choices.index(correct_choice)
        all_letters = ['A', 'B', 'C', 'D']
        correct_letter = all_letters[correct_choice_idx_in_question]
        final_correct_answer = f'{correct_letter}'
        new_questions.append({'question': question, 'answer': final_correct_answer, 'tool': 'reactome'})

        chosen_pair_idxs.append(random_pair_idx)

    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    output_path = 'benchmark_questions/reactome'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set2.csv', index = False)
    return



if __name__ == "__main__":
    num_questions = 20
    print("=== Generating Reactome Question Set 1 ===")
    reactome_set1(num_questions)
    print("=== Generating Reactome Question Set 2 ===")
    reactome_set2(num_questions)
    
    

    
