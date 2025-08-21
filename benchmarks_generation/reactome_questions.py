import pandas as pd
import random
import os
import sys
sys.path.append('../')
from tool_humanbase import _symbol_to_entrez

with open('../local_dbs/gene_names_list.txt', 'r') as f:
    gene_list = [line.strip() for line in f]
all_genes_count = len(gene_list)

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

if __name__ == "__main__":
    num_questions = 20
    print("=== Generating Reactome Question Set 1 ===")
    reactome_set1(num_questions)
    
    

    