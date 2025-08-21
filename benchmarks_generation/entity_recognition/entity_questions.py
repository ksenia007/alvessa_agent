import pandas as pd
import random
import os
import sys
sys.path.append('../')
import random

with open('../../local_dbs/gene_names_list.txt', 'r') as f:
    gene_list = [line.strip() for line in f]
all_genes_count = len(gene_list)

def random_case_string(gene):
    result_chars = []
    for char in gene:
        if char.isalpha():  
            if random.choice([True, False]):  
                result_chars.append(char.upper())
            else:
                result_chars.append(char.lower())
        else:
            result_chars.append(char) 
    return "".join(result_chars)

def set_case(gene, case_type):
    if case_type == 'original':
        return gene
    if case_type == 'lower':
        return gene.lower()
    if case_type == 'upper':
        return gene.upper()
    if case_type == 'mixed':
        return random_case_string(gene)
    return None

## What are variants for gene X?
def entity_set1(count):

    case_types = ['original', 'lower', 'upper', 'mixed']

    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:
        curr_case_type_idx = len(chosen_gene_idxs)%4

        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))
        random_gene_idx = random.sample(available_values, 1)[0]
        curr_gene = set_case(gene_list[random_gene_idx], case_types[curr_case_type_idx])

        query = f'What are variants for gene {curr_gene}?'

        genes = curr_gene

        new_questions.append({'query': query, 'recognized_genes': genes})

        chosen_gene_idxs.append(random_gene_idx)

    random.shuffle(new_questions)
    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    output_path = 'entity_recognition'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set1.csv', index = False)

## List common interactions for genes X, Y, and Z.
def entity_set2(count):

    case_types = ['original', 'lower', 'upper', 'mixed']

    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:
        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))
        random_gene_idxs = random.sample(available_values, 3)

        chosen_genes = [set_case(gene_list[i], random.sample(case_types, 1)[0]) for i in random_gene_idxs]

        query = f'List common interactions for genes {chosen_genes[0]}, {chosen_genes[1]}, and {chosen_genes[2]}.'
        new_questions.append({'query': query, 'recognized_genes': ",".join(chosen_genes)})
        chosen_gene_idxs.extend(random_gene_idxs)

    random.shuffle(new_questions)
    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    output_path = 'entity_recognition'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set2.csv', index = False)

## List all variants in X.
def entity_set3(count):

    case_types = ['original', 'lower', 'upper', 'mixed']

    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:
        curr_case_type_idx = len(chosen_gene_idxs)%4

        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))
        random_gene_idx = random.sample(available_values, 1)[0]
        curr_gene = set_case(gene_list[random_gene_idx], case_types[curr_case_type_idx])

        query = f'List all of the variants in {curr_gene}.'

        genes = curr_gene

        new_questions.append({'query': query, 'recognized_genes': genes})

        chosen_gene_idxs.append(random_gene_idx)

    random.shuffle(new_questions)
    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    output_path = 'entity_recognition'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set3.csv', index = False)

## What is common between X and Y?
def entity_set4(count):

    case_types = ['original', 'lower', 'upper', 'mixed']

    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:
        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))
        random_gene_idxs = random.sample(available_values, 2)

        chosen_genes = [set_case(gene_list[i], random.sample(case_types, 1)[0]) for i in random_gene_idxs]

        query = f'What is common between {chosen_genes[0]} and {chosen_genes[1]}?'
        new_questions.append({'query': query, 'recognized_genes': ",".join(chosen_genes)})
        chosen_gene_idxs.extend(random_gene_idxs)

    random.shuffle(new_questions)
    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    output_path = 'entity_recognition'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set4.csv', index = False)

## X is important for?
def entity_set5(count):

    case_types = ['original', 'lower', 'upper', 'mixed']

    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:
        curr_case_type_idx = len(chosen_gene_idxs)%4

        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))
        random_gene_idx = random.sample(available_values, 1)[0]
        curr_gene = set_case(gene_list[random_gene_idx], case_types[curr_case_type_idx])

        query = f'{curr_gene} is important for?'

        genes = curr_gene

        new_questions.append({'query': query, 'recognized_genes': genes})

        chosen_gene_idxs.append(random_gene_idx)

    random.shuffle(new_questions)
    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    output_path = 'entity_recognition'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set5.csv', index = False)

## What pathways are associated with X as well as Y but not with Z?
def entity_set6(count):

    case_types = ['original', 'lower', 'upper', 'mixed']

    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:
        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))
        random_gene_idxs = random.sample(available_values, 3)

        chosen_genes = [set_case(gene_list[i], random.sample(case_types, 1)[0]) for i in random_gene_idxs]

        query = f'What pathways are associated with {chosen_genes[0]} as well as {chosen_genes[1]} but not with {chosen_genes[2]}?'
        new_questions.append({'query': query, 'recognized_genes': ",".join(chosen_genes)})
        chosen_gene_idxs.extend(random_gene_idxs)

    random.shuffle(new_questions)
    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    output_path = 'entity_recognition'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set6.csv', index = False)

## Given gene X, which of the following are interacting? Y, Z.
def entity_set7(count):

    case_types = ['original', 'lower', 'upper', 'mixed']

    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:
        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))
        random_gene_idxs = random.sample(available_values, 3)

        chosen_genes = [set_case(gene_list[i], random.sample(case_types, 1)[0]) for i in random_gene_idxs]

        query = f'Given gene {chosen_genes[0]}, which of the following are interacting? {chosen_genes[1]}, {chosen_genes[2]}.'
        new_questions.append({'query': query, 'recognized_genes': ",".join(chosen_genes)})
        chosen_gene_idxs.extend(random_gene_idxs)

    random.shuffle(new_questions)
    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    output_path = 'entity_recognition'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set7.csv', index = False)

## What are the pathways associated with X, Y, Z… (random number of genes)
def entity_set8(count):

    case_types = ['original', 'lower', 'upper', 'mixed']

    chosen_gene_idxs = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:
        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))
        random_gene_idxs = random.sample(available_values, random.randint(5, 15))

        chosen_genes = [set_case(gene_list[i], random.sample(case_types, 1)[0]) for i in random_gene_idxs]

        query = f'What are the pathways associated with {",".join(chosen_genes)}'
        new_questions.append({'query': query, 'recognized_genes': ",".join(chosen_genes)})
        chosen_gene_idxs.extend(random_gene_idxs)

    random.shuffle(new_questions)
    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)

    output_path = 'entity_recognition'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set8.csv', index = False)

if __name__ == "__main__":
    num_questions = 50
    print("=== Generating Entity Questions Set 1 ===")
    entity_set1(num_questions)

    print("=== Generating Entity Questions Set 2 ===")
    entity_set2(num_questions)

    print("=== Generating Entity Questions Set 3 ===")
    entity_set3(num_questions)

    print("=== Generating Entity Questions Set 4 ===")
    entity_set4(num_questions)

    print("=== Generating Entity Questions Set 5 ===")
    entity_set5(num_questions)

    print("=== Generating Entity Questions Set 6 ===")
    entity_set6(num_questions)

    print("=== Generating Entity Questions Set 7 ===")
    entity_set7(num_questions)

    print("=== Generating Entity Questions Set 8 ===")
    entity_set8(num_questions)