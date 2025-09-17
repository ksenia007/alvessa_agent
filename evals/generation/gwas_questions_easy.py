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
    
    # Read only the columns we need to reduce memory usage
    columns_of_interest = ['DISEASE/TRAIT', 'MAPPED_GENE', 'P-VALUE', 'OR or BETA']
    
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

def clean_and_prepare_data(df: pd.DataFrame, p_value_threshold: float = 5e-8) -> pd.DataFrame:
    """Clean and prepare the data for question generation."""
    print("Cleaning and preparing data...")
    
    # Remove rows with missing values in critical columns
    df = df.dropna(subset=['DISEASE/TRAIT', 'MAPPED_GENE'])
    
    # Clean gene names - some entries have multiple genes separated by commas or semicolons
    # For simplicity, take the first gene mentioned
    df = df.copy()  # Avoid pandas warning
    df['MAPPED_GENE'] = df['MAPPED_GENE'].astype(str).str.split(',').str[0].str.strip()
    df['MAPPED_GENE'] = df['MAPPED_GENE'].str.split(';').str[0].str.strip()
    
    # Clean trait names
    df['DISEASE/TRAIT'] = df['DISEASE/TRAIT'].astype(str).str.strip()
    
    # Clean p-values and risk scores
    df['P-VALUE'] = pd.to_numeric(df['P-VALUE'], errors='coerce')
    df['OR or BETA'] = pd.to_numeric(df['OR or BETA'], errors='coerce')
    
    # Apply genome-wide significance filter (p-value < threshold)
    print(f"Before p-value filter: {len(df)} rows")
    df = df[df['P-VALUE'] < p_value_threshold]
    print(f"After p-value filter (< {p_value_threshold}): {len(df)} rows")
    
    # Remove entries with empty or very short gene names
    df = df[df['MAPPED_GENE'].str.len() > 1]
    df = df[df['DISEASE/TRAIT'].str.len() > 3]
    
    # Remove duplicate gene-trait pairs, keeping the one with the smallest p-value
    df = df.sort_values('P-VALUE').drop_duplicates(subset=['MAPPED_GENE', 'DISEASE/TRAIT'], keep='first')
    
    print(f"After cleaning: {len(df)} rows")
    return df

def get_gene_trait_mapping(df: pd.DataFrame) -> dict:
    """Create a mapping from genes to their associated traits with statistical data."""
    gene_traits = {}
    
    for _, row in df.iterrows():
        gene = row['MAPPED_GENE']
        trait = row['DISEASE/TRAIT']
        pvalue = row['P-VALUE']
        risk_score = row['OR or BETA']
        
        if gene not in gene_traits:
            gene_traits[gene] = {}
        
        # Store trait with its statistical information
        gene_traits[gene][trait] = {
            'p_value': pvalue,
            'risk_score': risk_score
        }
    
    # Filter genes that have at least one trait
    gene_traits = {gene: traits for gene, traits in gene_traits.items() if traits}
    
    print(f"Found {len(gene_traits)} genes with associated traits")
    return gene_traits

def get_all_traits(df: pd.DataFrame):
    """Get all unique traits for generating random options."""
    traits = df['DISEASE/TRAIT'].unique().tolist()
    print(f"Found {len(traits)} unique traits")
    return traits

def generate_question(gene: str, correct_traits: dict, all_traits):
    """Generate a single multiple-choice question for a gene."""
    # Pick one correct trait randomly
    correct_trait = random.choice(list(correct_traits.keys()))
    trait_data = correct_traits[correct_trait]
    
    # Generate 3 random incorrect traits
    incorrect_traits = []
    attempts = 0
    max_attempts = 100
    
    while len(incorrect_traits) < 3 and attempts < max_attempts:
        random_trait = random.choice(all_traits)
        if random_trait not in correct_traits and random_trait not in incorrect_traits:
            incorrect_traits.append(random_trait)
        attempts += 1
    
    if len(incorrect_traits) < 3:
        # If we couldn't find 3 unique incorrect traits, pad with generic ones
        generic_traits = [
            "Height measurement",
            "Body mass index",
            "Blood pressure measurement",
            "Cholesterol levels",
            "Blood glucose levels"
        ]
        for trait in generic_traits:
            if len(incorrect_traits) < 3 and trait not in correct_traits:
                incorrect_traits.append(trait)
    
    # Combine all options and shuffle
    all_options = [correct_trait] + incorrect_traits[:3]
    random.shuffle(all_options)
    
    # Find the correct answer index
    correct_index = all_options.index(correct_trait)
    correct_letter = chr(ord('A') + correct_index)
    
    # Format statistical values for display
    p_value = trait_data['p_value']
    risk_score = trait_data['risk_score']
    
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
        'question': f"Which of the following traits are associated with gene {gene}?",
        'options': all_options,
        'correct_answer': correct_letter,
        'correct_trait': correct_trait,
        'p_value': p_value_str,
        'risk_score': risk_score_str
    }

def gwas_set1(count):
    """Generate Set 1: Which of the following traits are associated with gene X? (one real, 3 random)"""
    
    # Load and process GWAS data
    gwas_file = LOCAL_DBS_DIR / 'gwas_catalogue_association.tsv'
    df = load_gwas_data(str(gwas_file))
    df = clean_and_prepare_data(df, p_value_threshold=5e-8)
    
    # Get gene-trait mapping and all traits
    gene_traits = get_gene_trait_mapping(df)
    all_traits = get_all_traits(df)
    
    # Filter genes that have at least 1 trait
    eligible_genes = [gene for gene, traits in gene_traits.items() if len(traits) >= 1]
    all_genes_count = len(eligible_genes)
    
    chosen_gene_idxs = []
    output_df = pd.DataFrame()
    new_questions = []
    
    while len(chosen_gene_idxs) < count:
        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))
        random_gene_idx = random.sample(available_values, 1)[0]
        curr_gene = eligible_genes[random_gene_idx]
        
        try:
            # Generate question using same logic as generate_mcq_questions.py
            question_data = generate_question(curr_gene, gene_traits[curr_gene], all_traits)
            
            # Format question exactly like the other benchmark files
            all_choices = question_data['options']
            question = f'Which of the following traits are associated with gene {curr_gene}? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'
            
            new_questions.append({
                'question': question, 
                'answer': question_data['correct_answer'], 
                'tool': 'GWAS',
                'gene': curr_gene,
                'correct_trait': question_data['correct_trait'],
                'p_value': question_data['p_value'],
                'risk_score': question_data['risk_score']
            })
            
            chosen_gene_idxs.append(random_gene_idx)
        except Exception as e:
            print(f"Error with gene {curr_gene}: {e}")
            continue
    
    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)
    
    output_path = 'benchmark_questions/GWAS'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set1.csv', index=False)

def gwas_set2(count):
    """Generate Set 2: Same format but different question instances"""
    
    # Load and process GWAS data
    gwas_file = LOCAL_DBS_DIR / 'gwas_catalogue_association.tsv'
    df = load_gwas_data(str(gwas_file))
    df = clean_and_prepare_data(df, p_value_threshold=5e-8)
    
    # Get gene-trait mapping and all traits
    gene_traits = get_gene_trait_mapping(df)
    all_traits = get_all_traits(df)
    
    # Filter genes that have at least 1 trait
    eligible_genes = [gene for gene, traits in gene_traits.items() if len(traits) >= 1]
    all_genes_count = len(eligible_genes)
    
    chosen_gene_idxs = []
    output_df = pd.DataFrame()
    new_questions = []
    
    while len(chosen_gene_idxs) < count:
        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))
        random_gene_idx = random.sample(available_values, 1)[0]
        curr_gene = eligible_genes[random_gene_idx]
        
        try:
            # Generate question using same logic as generate_mcq_questions.py
            question_data = generate_question(curr_gene, gene_traits[curr_gene], all_traits)
            
            # Format question exactly like the other benchmark files
            all_choices = question_data['options']
            question = f'Which of the following traits are associated with gene {curr_gene}? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'
            
            new_questions.append({
                'question': question, 
                'answer': question_data['correct_answer'], 
                'tool': 'GWAS',
                'gene': curr_gene,
                'correct_trait': question_data['correct_trait'],
                'p_value': question_data['p_value'],
                'risk_score': question_data['risk_score']
            })
            
            chosen_gene_idxs.append(random_gene_idx)
        except Exception as e:
            print(f"Error with gene {curr_gene}: {e}")
            continue
    
    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)
    
    output_path = 'benchmark_questions/GWAS'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set2.csv', index=False)

def gwas_set3(count):
    """Generate Set 3: Same format but different question instances"""
    
    # Load and process GWAS data
    gwas_file = LOCAL_DBS_DIR / 'gwas_catalogue_association.tsv'
    df = load_gwas_data(str(gwas_file))
    df = clean_and_prepare_data(df, p_value_threshold=5e-8)
    
    # Get gene-trait mapping and all traits
    gene_traits = get_gene_trait_mapping(df)
    all_traits = get_all_traits(df)
    
    # Filter genes that have at least 1 trait
    eligible_genes = [gene for gene, traits in gene_traits.items() if len(traits) >= 1]
    all_genes_count = len(eligible_genes)
    
    chosen_gene_idxs = []
    output_df = pd.DataFrame()
    new_questions = []
    
    while len(chosen_gene_idxs) < count:
        available_values = list(set(range(all_genes_count)) - set(chosen_gene_idxs))
        random_gene_idx = random.sample(available_values, 1)[0]
        curr_gene = eligible_genes[random_gene_idx]
        
        try:
            # Generate question using same logic as generate_mcq_questions.py
            question_data = generate_question(curr_gene, gene_traits[curr_gene], all_traits)
            
            # Format question exactly like the other benchmark files
            all_choices = question_data['options']
            question = f'Which of the following traits are associated with gene {curr_gene}? [A] {all_choices[0]} [B] {all_choices[1]} [C] {all_choices[2]} [D] {all_choices[3]}'
            
            new_questions.append({
                'question': question, 
                'answer': question_data['correct_answer'], 
                'tool': 'GWAS',
                'gene': curr_gene,
                'correct_trait': question_data['correct_trait'],
                'p_value': question_data['p_value'],
                'risk_score': question_data['risk_score']
            })
            
            chosen_gene_idxs.append(random_gene_idx)
        except Exception as e:
            print(f"Error with gene {curr_gene}: {e}")
            continue
    
    output_df = pd.concat([output_df, pd.DataFrame(new_questions)], ignore_index=True)
    
    output_path = 'benchmark_questions/GWAS'
    os.makedirs(output_path, exist_ok=True)
    output_df.to_csv(f'{output_path}/set3.csv', index=False)

if __name__ == "__main__":
    # Set random seed for reproducibility (same as generate_mcq_questions.py)
    random.seed(42)
    
    num_questions = 20
    print("=== Generating GWAS Question Set 1 ===")
    gwas_set1(num_questions)
    print("=== Generating GWAS Question Set 2 ===") 
    gwas_set2(num_questions)
    print("=== Generating GWAS Question Set 3 ===")
    gwas_set3(num_questions)
