import pandas as pd
import random
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

LOCAL_DBS_DIR = ROOT / "local_dbs"

def load_gwas_data(file_path: str) -> pd.DataFrame:
    """Load GWAS catalogue data from TSV file."""
    print(f"Loading GWAS data from {file_path}...")
    
    # Read columns needed for variant-based questions
    columns_of_interest = ['DISEASE/TRAIT', 'MAPPED_GENE', 'P-VALUE', 'OR or BETA', 
                          'SNPS', 'STRONGEST SNP-RISK ALLELE', 'CONTEXT', 'CHR_ID', 'CHR_POS']
    
    try:
        # Use chunking to handle large file
        chunk_list = []
        chunk_size = 10000
        
        for chunk in pd.read_csv(file_path, sep='\t', chunksize=chunk_size, 
                                usecols=columns_of_interest, low_memory=False):
            chunk_list.append(chunk)
            
        df = pd.concat(chunk_list, ignore_index=True)
        print(f"Loaded {len(df)} rows")
        
        return df
    except Exception as e:
        print(f"Error loading data: {e}")
        sys.exit(1)

def clean_and_prepare_variant_data(df: pd.DataFrame, p_value_threshold: float = 5e-8) -> pd.DataFrame:
    """Clean and prepare the data for variant-based question generation."""
    print("Cleaning and preparing variant data...")
    
    # Remove rows with missing values in critical columns
    df = df.dropna(subset=['DISEASE/TRAIT', 'MAPPED_GENE', 'SNPS'])
    
    # Clean gene names - some entries have multiple genes separated by commas or semicolons
    # For simplicity, take the first gene mentioned
    df = df.copy()  # Avoid pandas warning
    df['MAPPED_GENE'] = df['MAPPED_GENE'].astype(str).str.split(',').str[0].str.strip()
    df['MAPPED_GENE'] = df['MAPPED_GENE'].str.split(';').str[0].str.strip()
    
    # Clean trait names
    df['DISEASE/TRAIT'] = df['DISEASE/TRAIT'].astype(str).str.strip()
    
    # Clean variant names - extract primary SNP ID
    df['SNPS'] = df['SNPS'].astype(str).str.strip()
    
    # Clean p-values and risk scores
    df['P-VALUE'] = pd.to_numeric(df['P-VALUE'], errors='coerce')
    df['OR or BETA'] = pd.to_numeric(df['OR or BETA'], errors='coerce')
    
    # Apply genome-wide significance filter (p-value < threshold)
    print(f"Before p-value filter: {len(df)} rows")
    df = df[df['P-VALUE'] < p_value_threshold]
    print(f"After p-value filter (< {p_value_threshold}): {len(df)} rows")
    
    # Remove entries with empty or very short names
    df = df[df['MAPPED_GENE'].str.len() > 1]
    df = df[df['DISEASE/TRAIT'].str.len() > 3]
    df = df[df['SNPS'].str.len() > 3]
    
    # Remove entries with excluded gene names
    excluded_genes = {'NR', 'INTERGENIC', '', 'nan'}
    df = df[~df['MAPPED_GENE'].isin(excluded_genes)]
    
    # Remove duplicate variant-gene-trait combinations, keeping the one with the smallest p-value
    df = df.sort_values('P-VALUE').drop_duplicates(subset=['SNPS', 'MAPPED_GENE', 'DISEASE/TRAIT'], keep='first')
    
    print(f"After cleaning: {len(df)} rows")
    return df

def get_variant_gene_trait_mapping(df: pd.DataFrame) -> dict:
    """Create a mapping from (gene, trait) pairs to their associated variants with statistical data."""
    gene_trait_variants = {}
    
    for _, row in df.iterrows():
        gene = row['MAPPED_GENE']
        trait = row['DISEASE/TRAIT']
        variant = row['SNPS']
        pvalue = row['P-VALUE']
        risk_score = row['OR or BETA']
        context = row.get('CONTEXT', 'Unknown')
        
        # Create gene-trait key
        gene_trait_key = (gene, trait)
        
        if gene_trait_key not in gene_trait_variants:
            gene_trait_variants[gene_trait_key] = {}
        
        # Store variant with its statistical information
        gene_trait_variants[gene_trait_key][variant] = {
            'p_value': pvalue,
            'risk_score': risk_score,
            'context': context
        }
    
    # Filter gene-trait pairs that have at least one variant
    gene_trait_variants = {key: variants for key, variants in gene_trait_variants.items() if variants}
    
    print(f"Found {len(gene_trait_variants)} gene-trait pairs with associated variants")
    return gene_trait_variants

def get_all_variants(df: pd.DataFrame):
    """Get all unique variants for generating random options."""
    variants = df['SNPS'].unique().tolist()
    # Filter out variants that are too short or invalid
    variants = [v for v in variants if len(str(v)) > 3 and str(v) != 'nan']
    print(f"Found {len(variants)} unique variants")
    return variants

def generate_variant_question(gene_trait_key: tuple, correct_variants: dict, all_variants: list):
    """Generate a single multiple-choice question for a gene-trait pair asking about variants."""
    gene, trait = gene_trait_key
    
    # Pick one correct variant randomly
    correct_variant = random.choice(list(correct_variants.keys()))
    variant_data = correct_variants[correct_variant]
    
    # Generate 3 random incorrect variants
    incorrect_variants = []
    attempts = 0
    max_attempts = 100
    
    while len(incorrect_variants) < 3 and attempts < max_attempts:
        random_variant = random.choice(all_variants)
        if random_variant not in correct_variants and random_variant not in incorrect_variants:
            incorrect_variants.append(random_variant)
        attempts += 1
    
    if len(incorrect_variants) < 3:
        # If we couldn't find 3 unique incorrect variants, pad with generic SNP IDs
        generic_variants = [
            "rs1234567",
            "rs2345678", 
            "rs3456789",
            "rs4567890",
            "rs5678901"
        ]
        for variant in generic_variants:
            if len(incorrect_variants) < 3 and variant not in correct_variants:
                incorrect_variants.append(variant)
    
    # Combine all options and shuffle
    all_options = [correct_variant] + incorrect_variants[:3]
    random.shuffle(all_options)
    
    # Find the correct answer index
    correct_index = all_options.index(correct_variant)
    correct_letter = chr(ord('A') + correct_index)
    
    # Format statistical values for display
    p_value = variant_data['p_value']
    risk_score = variant_data['risk_score']
    context = variant_data.get('context', 'Unknown')
    
    # Format p-value for display
    if pd.isna(p_value):
        p_value_str = "Not available"
    elif p_value < 0.001:
        p_value_str = f"{p_value:.2e}"
    else:
        p_value_str = f"{p_value:.3f}"
    
    # Format risk score for display  
    if pd.isna(risk_score):
        risk_score_str = "Not available"
    else:
        risk_score_str = f"{risk_score:.2f}"
    
    return {
        'gene': gene,
        'trait': trait,
        'question': f"Out of these variants, which one is associated with trait '{trait}' for gene {gene}?",
        'options': all_options,
        'correct_answer': correct_letter,
        'correct_variant': correct_variant,
        'p_value': p_value_str,
        'risk_score': risk_score_str,
        'context': context
    }

def gwas_variant_set1(count):
    """Generate Set 1: Out of these variants, which one is associated with trait X for gene Y? (one real, 3 random)"""
    
    # Load and process GWAS data
    gwas_file = LOCAL_DBS_DIR / 'gwas_catalogue_association.tsv'
    df = load_gwas_data(str(gwas_file))
    df = clean_and_prepare_variant_data(df, p_value_threshold=5e-8)
    
    # Get variant-gene-trait mapping and all variants
    gene_trait_variants = get_variant_gene_trait_mapping(df)
    all_variants = get_all_variants(df)
    
    # Filter gene-trait pairs that have at least 1 variant
    eligible_gene_trait_pairs = [key for key, variants in gene_trait_variants.items() if len(variants) >= 1]
    all_pairs_count = len(eligible_gene_trait_pairs)
    
    chosen_pair_idxs = []
    output_df = pd.DataFrame()
    new_questions = []
    
    while len(chosen_pair_idxs) < count:
        available_values = list(set(range(all_pairs_count)) - set(chosen_pair_idxs))
        if not available_values:
            print(f"Warning: Only found {len(chosen_pair_idxs)} eligible gene-trait pairs, requested {count}")
            break
            
        random_pair_idx = random.sample(available_values, 1)[0]
        curr_gene_trait_key = eligible_gene_trait_pairs[random_pair_idx]
        
        try:
            # Generate question
            question_data = generate_variant_question(curr_gene_trait_key, gene_trait_variants[curr_gene_trait_key], all_variants)
            
            # Format question exactly like the other benchmark files
            all_choices = question_data['options']
            question = f"Out of these variants, which one is associated with trait '{question_data['trait']}' for gene {question_data['gene']}? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}"
            
            new_questions.append({
                'question': question, 
                'answer': question_data['correct_answer'], 
                'tool': 'GWAS',
                'gene': question_data['gene'],
                'trait': question_data['trait'],
                'correct_variant': question_data['correct_variant'],
                'p_value': question_data['p_value'],
                'risk_score': question_data['risk_score'],
                'context': question_data['context']
            })
            
            chosen_pair_idxs.append(random_pair_idx)
        except Exception as e:
            print(f"Error with gene-trait pair {curr_gene_trait_key}: {e}")
            continue
    
    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)
    
    output_path = 'benchmark_questions/GWAS-variants'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set1.csv', index=False)

def gwas_variant_set2(count):
    """Generate Set 2: Same format but different question instances"""
    
    # Load and process GWAS data
    gwas_file = LOCAL_DBS_DIR / 'gwas_catalogue_association.tsv'
    df = load_gwas_data(str(gwas_file))
    df = clean_and_prepare_variant_data(df, p_value_threshold=5e-8)
    
    # Get variant-gene-trait mapping and all variants
    gene_trait_variants = get_variant_gene_trait_mapping(df)
    all_variants = get_all_variants(df)
    
    # Filter gene-trait pairs that have at least 1 variant
    eligible_gene_trait_pairs = [key for key, variants in gene_trait_variants.items() if len(variants) >= 1]
    all_pairs_count = len(eligible_gene_trait_pairs)
    
    chosen_pair_idxs = []
    output_df = pd.DataFrame()
    new_questions = []
    
    while len(chosen_pair_idxs) < count:
        available_values = list(set(range(all_pairs_count)) - set(chosen_pair_idxs))
        if not available_values:
            print(f"Warning: Only found {len(chosen_pair_idxs)} eligible gene-trait pairs, requested {count}")
            break
            
        random_pair_idx = random.sample(available_values, 1)[0]
        curr_gene_trait_key = eligible_gene_trait_pairs[random_pair_idx]
        
        try:
            # Generate question
            question_data = generate_variant_question(curr_gene_trait_key, gene_trait_variants[curr_gene_trait_key], all_variants)
            
            # Format question exactly like the other benchmark files
            all_choices = question_data['options']
            question = f"Out of these variants, which one is associated with trait '{question_data['trait']}' for gene {question_data['gene']}? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}"
            
            new_questions.append({
                'question': question, 
                'answer': question_data['correct_answer'], 
                'tool': 'GWAS',
                'gene': question_data['gene'],
                'trait': question_data['trait'],
                'correct_variant': question_data['correct_variant'],
                'p_value': question_data['p_value'],
                'risk_score': question_data['risk_score'],
                'context': question_data['context']
            })
            
            chosen_pair_idxs.append(random_pair_idx)
        except Exception as e:
            print(f"Error with gene-trait pair {curr_gene_trait_key}: {e}")
            continue
    
    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)
    
    output_path = 'benchmark_questions/GWAS-variants'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set2.csv', index=False)

def gwas_variant_set3(count):
    """Generate Set 3: Same format but different question instances"""
    
    # Load and process GWAS data
    gwas_file = LOCAL_DBS_DIR / 'gwas_catalogue_association.tsv'
    df = load_gwas_data(str(gwas_file))
    df = clean_and_prepare_variant_data(df, p_value_threshold=5e-8)
    
    # Get variant-gene-trait mapping and all variants
    gene_trait_variants = get_variant_gene_trait_mapping(df)
    all_variants = get_all_variants(df)
    
    # Filter gene-trait pairs that have at least 1 variant
    eligible_gene_trait_pairs = [key for key, variants in gene_trait_variants.items() if len(variants) >= 1]
    all_pairs_count = len(eligible_gene_trait_pairs)
    
    chosen_pair_idxs = []
    output_df = pd.DataFrame()
    new_questions = []
    
    while len(chosen_pair_idxs) < count:
        available_values = list(set(range(all_pairs_count)) - set(chosen_pair_idxs))
        if not available_values:
            print(f"Warning: Only found {len(chosen_pair_idxs)} eligible gene-trait pairs, requested {count}")
            break
            
        random_pair_idx = random.sample(available_values, 1)[0]
        curr_gene_trait_key = eligible_gene_trait_pairs[random_pair_idx]
        
        try:
            # Generate question
            question_data = generate_variant_question(curr_gene_trait_key, gene_trait_variants[curr_gene_trait_key], all_variants)
            
            # Format question exactly like the other benchmark files
            all_choices = question_data['options']
            question = f"Out of these variants, which one is associated with trait '{question_data['trait']}' for gene {question_data['gene']}? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}"
            
            new_questions.append({
                'question': question, 
                'answer': question_data['correct_answer'], 
                'tool': 'GWAS',
                'gene': question_data['gene'],
                'trait': question_data['trait'],
                'correct_variant': question_data['correct_variant'],
                'p_value': question_data['p_value'],
                'risk_score': question_data['risk_score'],
                'context': question_data['context']
            })
            
            chosen_pair_idxs.append(random_pair_idx)
        except Exception as e:
            print(f"Error with gene-trait pair {curr_gene_trait_key}: {e}")
            continue
    
    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)
    
    output_path = 'benchmark_questions/GWAS-variants'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set3.csv', index=False)

if __name__ == "__main__":
    # Set random seed for reproducibility (same as generate_mcq_questions.py)
    random.seed(42)
    
    num_questions = 20
    print("=== Generating GWAS Variant Question Set 1 ===")
    gwas_variant_set1(num_questions)
    print("=== Generating GWAS Variant Question Set 2 ===") 
    gwas_variant_set2(num_questions)
    print("=== Generating GWAS Variant Question Set 3 ===")
    gwas_variant_set3(num_questions)
