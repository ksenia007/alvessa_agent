import pandas as pd
import random
import os
from collections import defaultdict
from utils import gwas_associations, get_variant_coords, check_key_in_nested_dicts, _symbol_to_uniprot

with open('../../../local_dbs/gene_names_list.txt', 'r') as f:
    gene_list = [line.strip() for line in f]
all_genes_count = len(gene_list)

pathogenicity_class_df_hg38 = pd.read_parquet('../../../local_dbs/AlphaMissense_hg38.parquet')

with open('../../../local_dbs/pathogenic_genes_list.txt', 'r') as f:
    pathogenic_genes = [line.strip() for line in f]

def gwas_alphamissense_set1(count):
    chosen_genes = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:

        curr_genes = {}

        selected_from_pathogenic = random.sample(list(pathogenic_genes), 1)[0]
        remaining_candidates = set(pathogenic_genes) - {selected_from_pathogenic}
        selected_remaining = random.sample(list(remaining_candidates), 1)[0]
        final_selection = {selected_from_pathogenic:True, selected_remaining:False}
        correct_choice = ""
        false_choices = []
        main_correct_choices = set()
        main_false_choices = set()
        alt_false_choices = set()

        should_continue = False

        main_var_ids = []

        for random_gene, main_gene in final_selection.items():

            variants = {}
            associations = {}

            succeeded, variants, associations = gwas_associations(random_gene, variants, associations)

            if succeeded:
                variants = get_variant_coords(variants, associations)
                if len(variants)>0 and check_key_in_nested_dicts(variants, 'coordinates'):
                    curr_genes[(random_gene, main_gene)] = variants[random_gene]
                else:
                    should_continue = True
                    break
            else:
                should_continue = True
                break

        if should_continue:
            continue
        
        snp_records = []
        for (gene, main_gene), gene_vars in curr_genes.items():
            for var_id, var_data in gene_vars.items():
                all_snps = var_data.get("coordinates", [])
                for snp in all_snps:
                    chrom, pos, ref, alt, assembly = snp.get("chrom"), snp.get("pos"), snp.get("ref"), snp.get("alt"), snp.get("assembly")

                    if assembly and "GRCh38" in assembly:
                        if None not in (chrom, pos, ref, alt):
                            snp_records.append({
                                "gene": f'{gene},{int(main_gene)}',
                                "uniprot_IDs": _symbol_to_uniprot(gene),
                                "var_id": var_id,
                                "snp_key": f"SNP:{ref}->{alt}",
                                "chrom": f"chr{chrom}",
                                "pos": pos,
                                "ref": ref,
                                "alt": alt
                            })
                            if main_gene:
                                main_var_ids.append(var_id)
                        else:
                            print(f"Missing coordinate data for {gene} variant {var_id} (SNP {ref}->{alt})")

        snps_df = pd.DataFrame(snp_records)
        
        snps_exploded = snps_df.explode("uniprot_IDs")
        merged = snps_exploded.merge(
            pathogenicity_class_df_hg38[['chrom', 'pos', 'ref', 'alt', 'uniprot_id', 'am_class']],
            left_on=['chrom', 'pos', 'ref', 'alt', 'uniprot_IDs'],
            right_on=['chrom', 'pos', 'ref', 'alt', 'uniprot_id'],
            how='left'
        )

        grouped = merged.groupby(['gene', 'var_id', 'snp_key'], as_index=False).agg({
        'am_class': lambda x: next(iter(set(filter(pd.notna, x))), None)
        })

        logs = defaultdict(set)
        for row in grouped.itertuples(index=False):
            gene_tag, var_id, snp_key, am_class = row

            gene, main_gene = gene_tag.split(",")
            main_gene = int(main_gene)

            if main_gene:
                logs[var_id].add(am_class)

            if am_class=='likely_pathogenic':
                if main_gene:
                    main_correct_choices.add(var_id)
                else:   
                    if var_id not in main_var_ids:
                        alt_false_choices.add(var_id)
            else:
                if main_gene:
                    main_false_choices.add(var_id)

        if len(main_correct_choices)<1 or len(main_false_choices)<2 or len(alt_false_choices)<1:
            continue

        correct_choice = random.sample(list(main_correct_choices), 1)[0]
        false_choices = random.sample(list(main_false_choices), 2) + random.sample(list(alt_false_choices), 1)
        all_choices = false_choices + [correct_choice]
        random.shuffle(all_choices)

        question = f'Which of the following variants is associated with gene {selected_from_pathogenic} and has the worst possible predicted coding downstream effect according to AlphaMissense? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'

        correct_choice_idx_in_question = all_choices.index(correct_choice)
        all_letters = ['A', 'B', 'C', 'D']
        correct_letter = all_letters[correct_choice_idx_in_question]
        final_correct_answer = f'{correct_letter}'
        new_questions.append({'question': question, 'answer': final_correct_answer, 'tool': 'query_gwas_extensive,alphamissense', 'tool_specified_in_question': True})

        chosen_genes.append(selected_from_pathogenic)

    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)
    
    output_path = '../../../benchmarks_generation/benchmark_questions/gwas_alphamissense'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set1.csv', index = False)


def gwas_alphamissense_set2(count):
    chosen_genes = []

    output_df = pd.DataFrame()
    new_questions = []

    while len(new_questions)<count:

        curr_genes = {}

        selected_from_pathogenic = random.sample(list(pathogenic_genes), 1)
        remaining_candidates = (set(gene_list) - set(pathogenic_genes)) - set(chosen_genes)
        selected_remaining = random.sample(list(remaining_candidates), 3)
        final_selection = selected_from_pathogenic + selected_remaining
        random.shuffle(final_selection)

        should_continue = False

        for random_gene in final_selection:
            variants = {}
            associations = {}

            succeeded, variants, associations = gwas_associations(random_gene, variants, associations)
            if succeeded:
                variants = get_variant_coords(variants, associations)
                if len(variants)>0 and check_key_in_nested_dicts(variants, 'coordinates'):
                    curr_genes[random_gene] = variants[random_gene]
                else:
                    should_continue = True
                    break
            else:
                should_continue = True
                break

        if should_continue:
            continue

        snp_records = []
        for gene, gene_vars in curr_genes.items():
            for var_id, var_data in gene_vars.items():
                all_snps = var_data.get("coordinates", [])
                for snp in all_snps:
                    chrom, pos, ref, alt, assembly = snp.get("chrom"), snp.get("pos"), snp.get("ref"), snp.get("alt"), snp.get("assembly")

                    if assembly and "GRCh38" in assembly:
                        if None not in (chrom, pos, ref, alt):
                            snp_records.append({
                                "gene": gene,
                                "uniprot_IDs": _symbol_to_uniprot(gene),
                                "var_id": var_id,
                                "snp_key": f"SNP:{ref}->{alt}",
                                "chrom": f"chr{chrom}",
                                "pos": pos,
                                "ref": ref,
                                "alt": alt
                            })
                        else:
                            print(f"Missing coordinate data for {gene} variant {var_id} (SNP {ref}->{alt})")

        snps_df = pd.DataFrame(snp_records)
        
        snps_exploded = snps_df.explode("uniprot_IDs")
        merged = snps_exploded.merge(
            pathogenicity_class_df_hg38[['chrom', 'pos', 'ref', 'alt', 'uniprot_id', 'am_class']],
            left_on=['chrom', 'pos', 'ref', 'alt', 'uniprot_IDs'],
            right_on=['chrom', 'pos', 'ref', 'alt', 'uniprot_id'],
            how='left'
        )

        grouped = merged.groupby(['gene', 'var_id', 'snp_key'], as_index=False).agg({
        'am_class': lambda x: next(iter(set(filter(pd.notna, x))), None)
        })

        counted = defaultdict(list)

        gene_pathogenic_counts = {}
        for row in grouped.itertuples(index=False):
            gene, var_id, snp_key, am_class = row

            if gene in counted and var_id in counted[gene]:
                continue

            if gene not in gene_pathogenic_counts:
                gene_pathogenic_counts[gene] = 0

            if am_class=='likely_pathogenic':
                gene_pathogenic_counts[gene] += 1
                counted[gene].append(var_id) 

        gene_pathogenic_counts_desc = sorted(gene_pathogenic_counts.items(), key=lambda item: item[1], reverse=True)

        if gene_pathogenic_counts_desc[0][1]==0 or (gene_pathogenic_counts_desc[0][1]==gene_pathogenic_counts_desc[1][1]):
            continue

        correct_choice = gene_pathogenic_counts_desc[0][0]
        false_choices = [gene_pathogenic_counts_desc[1][0], gene_pathogenic_counts_desc[2][0], gene_pathogenic_counts_desc[3][0]]

        all_choices = false_choices + [correct_choice]
        random.shuffle(all_choices)

        question = f'Which of the following genes has the most number of coding variants associated with any trait that are predicted by AlphaMissense to be pathogenic? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'

        correct_choice_idx_in_question = all_choices.index(correct_choice)
        all_letters = ['A', 'B', 'C', 'D']
        correct_letter = all_letters[correct_choice_idx_in_question]
        final_correct_answer = f'{correct_letter}'
        new_questions.append({'question': question, 'answer': final_correct_answer, 'tool': 'query_gwas_extensive,alphamissense', 'tool_specified_in_question': True})

        chosen_genes.append(selected_from_pathogenic[0])

    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)
    
    output_path = '../../../benchmarks_generation/benchmark_questions/gwas_alphamissense'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set2.csv', index = False)

 

if __name__ == "__main__":
    num_questions = 20
    print("=== Generating GWAS+AlphaMissense Question Set 1 ===")
    gwas_alphamissense_set1(num_questions)
    print("=== Generating GWAS+AlphaMissense Question Set 2 ===")
    gwas_alphamissense_set2(num_questions)
    
    

    