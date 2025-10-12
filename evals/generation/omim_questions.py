import pandas as pd
import random
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

LOCAL_DBS_DIR = ROOT / 'local_dbs'
gene_to_phen = LOCAL_DBS_DIR / "omim_genemap2.txt"
OUTPUT = ROOT / 'benchmarks_generation/benchmark_questions/OMIM'

gene_to_phen_df = pd.read_csv(gene_to_phen, sep = '\t', comment='#', header=None, names = ['Chromosome','Genomic Position Start','Genomic Position End','Cyto Location','Computed Cyto Location','MIM Number','Gene/Locus And Other Related Symbols','Gene Name','Approved Gene Symbol','Entrez Gene ID','Ensembl Gene ID','Comments','Phenotypes','Mouse Gene Symbol/ID'])

# Fixing weird formatting they have for phenotypes
pattern = r'(, \d+ \(\d+\))|(\s\(\d+\))'
gene_to_phen_df['extracted'] = gene_to_phen_df.loc[:,'Phenotypes'].str.replace(pattern, '', regex=True).str.replace('[\{\}\[\]\?]', '', regex=True)
    

genedf_exploded = gene_to_phen_df.assign(phenotype=gene_to_phen_df['extracted'].str.split(';')).explode('extracted')

# Creating dataframe of phenotypes with at least 2 associated genes
phen_to_gene_df = (
    genedf_exploded.groupby('extracted')['Approved Gene Symbol']
    .agg(lambda x: ','.join(sorted(set(map(str, x.dropna())))))
    .reset_index()
)

phen_to_gene_df.columns = ['Phenotype', 'Approved Gene Symbols']
gene_counts = phen_to_gene_df['Approved Gene Symbols'].str.count(',') + 1
phen_to_gene_df = phen_to_gene_df[gene_counts >= 2]

all_phenotypes = list(set(phen_to_gene_df['Phenotype'].values))
all_genes = list(set(gene_to_phen_df['Approved Gene Symbol'].values))


# Which of the following phenotypes is associated with gene X?
def omim_set1(count):
    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:
        available_values = list(set(range(len(all_genes))) - set(chosen_gene_idxs))
        random_gene_idx = random.sample(available_values, 1)[0]
        curr_gene = all_genes[random_gene_idx]

        match = gene_to_phen_df[gene_to_phen_df['Approved Gene Symbol']==curr_gene]
        
        try:
            gene_phenotypes = match['extracted'].values[0].split(';')

            random_phenotype_idx = random.sample(range(len(gene_phenotypes)), 1)[0]

            correct_choice = list(gene_phenotypes)[random_phenotype_idx]

            gene_phenotypes_idxs = []
            for phenotype in gene_phenotypes:
                gene_phenotypes_idxs.append(all_phenotypes.index(phenotype))

            available_values = list(set(range(len(all_phenotypes))) - set(gene_phenotypes_idxs))
            false_choice_idxs = random.sample(available_values, 3)

            false_choices = [list(all_phenotypes)[i] for i in false_choice_idxs]

            all_choices = false_choices + [correct_choice]
            random.shuffle(all_choices)

            question = f'Which of the following phenotypes is associated with gene {curr_gene}? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'

            correct_choice_idx_in_question = all_choices.index(correct_choice)
            all_letters = ['A', 'B', 'C', 'D']
            correct_letter = all_letters[correct_choice_idx_in_question]
            final_correct_answer = f'{correct_letter}'
            new_questions.append({'question': question, 'answer': final_correct_answer, 'tool': 'OMIM'})

            chosen_gene_idxs.append(random_gene_idx)

        except:
            continue 

    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    os.makedirs(OUTPUT, exist_ok=True)
    output_df.to_csv(f'{OUTPUT}/set1.csv', index = False)
    return

# Which of the following phenotypes is associated with both gene X and gene Y?
def omim_set2(count):
    chosen_phenotype_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:
        available_values = list(set(range(len(all_phenotypes))) - set(chosen_phenotype_idxs))
        random_phenotype_idx = random.sample(available_values, 1)[0]
        curr_phenotype = all_phenotypes[random_phenotype_idx]

        match = phen_to_gene_df[phen_to_gene_df['Phenotype']==curr_phenotype]

        
        try:
            phenotype_genes = match['Approved Gene Symbols'].values[0].split(',')

            if len(phenotype_genes) >= 2:
                random_genes_idxs = random.sample(range(len(phenotype_genes)), 2)

                gene_x = phenotype_genes[random_genes_idxs[0]]
                gene_y = phenotype_genes[random_genes_idxs[1]]

                correct_choice = curr_phenotype

                if gene_x in all_genes and gene_y in all_genes:
                    try: 
                        match_x = gene_to_phen_df[gene_to_phen_df['Approved Gene Symbol']==gene_x]['extracted'].values[0].split(';')
                    except:
                        match_x = []

                    try:
                        match_y = gene_to_phen_df[gene_to_phen_df['Approved Gene Symbol']==gene_y]['extracted'].values[0].split(';')
                    except:
                        match_y = []

                    intersection = [x for x in match_x if x in match_y]

                    if curr_phenotype in intersection:

                        intersection_idxs = []
                        for phenotype in intersection:
                            intersection_idxs.append(all_phenotypes.index(phenotype))

                        available_values = list(set(range(len(all_phenotypes))) - set(intersection_idxs))
                        false_choice_idxs = random.sample(available_values, 3)

                        false_choices = [list(all_phenotypes)[i] for i in false_choice_idxs]

                        all_choices = false_choices + [correct_choice]
                        random.shuffle(all_choices)

                        question = f'Which of the following phenotypes is associated with both gene {gene_x} and gene {gene_y}? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'

                        correct_choice_idx_in_question = all_choices.index(correct_choice)
                        all_letters = ['A', 'B', 'C', 'D']
                        correct_letter = all_letters[correct_choice_idx_in_question]
                        final_correct_answer = f'{correct_letter}'
                        new_questions.append({'question': question, 'answer': final_correct_answer, 'tool': 'OMIM'})

                        chosen_phenotype_idxs.append(random_phenotype_idx)

        except Exception as e:
            continue 

    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    os.makedirs(OUTPUT, exist_ok=True)
    output_df.to_csv(f'{OUTPUT}/set2.csv', index = False)
    return

if __name__ == "__main__":
    num_questions = 20
    print("=== Generating OMIM Question Set 1 ===")
    omim_set1(num_questions)
    print("=== Generating OMIM Question Set 2 ===")
    omim_set2(num_questions)
