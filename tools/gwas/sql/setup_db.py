#!/usr/bin/env python3
"""
Script to convert GWAS TSV files to SQLite database for faster querying.
Run this once to set up the database, then use the optimized query engine.
"""

import sqlite3
import csv
import os
import sys
from pathlib import Path

def create_gwas_database(db_path: str = "local_dbs/gwas.db", tsv_dir: str = "local_dbs"):
    """
    Convert GWAS TSV files to SQLite database with proper indexing.
    
    Args:
        db_path: Path where to create the SQLite database
        tsv_dir: Directory containing the TSV files
    """
    
    associations_file = os.path.join(tsv_dir, "gwas_catalogue_association.tsv")
    studies_file = os.path.join(tsv_dir, "gwas_catalogue_studies.tsv")
    
    if not os.path.exists(associations_file) or not os.path.exists(studies_file):
        print(f"Error: TSV files not found in {tsv_dir}")
        return False
    
    # Remove existing database
    if os.path.exists(db_path):
        os.remove(db_path)
        print(f"Removed existing database: {db_path}")
    
    # Create database connection
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        # Create associations table
        print("Creating associations table...")
        cursor.execute('''
            CREATE TABLE associations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                date_added TEXT,
                pubmed_id TEXT,
                first_author TEXT,
                date TEXT,
                journal TEXT,
                link TEXT,
                study TEXT,
                disease_trait TEXT,
                initial_sample_size TEXT,
                replication_sample_size TEXT,
                region TEXT,
                chr_id TEXT,
                chr_pos TEXT,
                reported_genes TEXT,
                mapped_gene TEXT,
                upstream_gene_id TEXT,
                downstream_gene_id TEXT,
                snp_gene_ids TEXT,
                upstream_gene_distance TEXT,
                downstream_gene_distance TEXT,
                strongest_snp_risk_allele TEXT,
                snps TEXT,
                merged TEXT,
                snp_id_current TEXT,
                context TEXT,
                intergenic TEXT,
                risk_allele_frequency TEXT,
                p_value TEXT,
                p_value_mlog TEXT,
                p_value_text TEXT,
                or_beta TEXT,
                ci_95_text TEXT,
                platform TEXT,
                cnv TEXT
            )
        ''')
        
        # Load associations data
        print("Loading associations data...")
        with open(associations_file, 'r', encoding='utf-8', errors='ignore') as f:
            # Read header to get column names
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            
            # Create placeholders for all columns
            placeholders = ','.join(['?' for _ in header])
            insert_query = f"INSERT INTO associations ({','.join(['date_added', 'pubmed_id', 'first_author', 'date', 'journal', 'link', 'study', 'disease_trait', 'initial_sample_size', 'replication_sample_size', 'region', 'chr_id', 'chr_pos', 'reported_genes', 'mapped_gene', 'upstream_gene_id', 'downstream_gene_id', 'snp_gene_ids', 'upstream_gene_distance', 'downstream_gene_distance', 'strongest_snp_risk_allele', 'snps', 'merged', 'snp_id_current', 'context', 'intergenic', 'risk_allele_frequency', 'p_value', 'p_value_mlog', 'p_value_text', 'or_beta', 'ci_95_text', 'platform', 'cnv'])}) VALUES ({placeholders})"
            
            # Insert data in batches
            batch_size = 1000
            batch = []
            count = 0
            
            for row in reader:
                # Pad row to match expected columns
                while len(row) < 34:
                    row.append('')
                
                batch.append(row[:34])  # Take only first 34 columns
                count += 1
                
                if len(batch) >= batch_size:
                    cursor.executemany(insert_query, batch)
                    batch = []
                    if count % 10000 == 0:
                        print(f"  Processed {count} associations...")
            
            # Insert remaining batch
            if batch:
                cursor.executemany(insert_query, batch)
            
            print(f"  Loaded {count} associations total")
        
        # Create studies table
        print("Creating studies table...")
        cursor.execute('''
            CREATE TABLE studies (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                study_accession TEXT,
                pubmed_id TEXT,
                first_author TEXT,
                date TEXT,
                initial_sample_description TEXT,
                replication_sample_description TEXT,
                stage TEXT,
                number_of_individuals TEXT,
                broad_ancestral_category TEXT,
                country_of_origin TEXT,
                country_of_recruitment TEXT,
                additional_ancestry_description TEXT,
                ancestry_descriptor TEXT,
                founder_genetically_isolated_population TEXT,
                number_of_cases TEXT,
                number_of_controls TEXT,
                sample_description TEXT
            )
        ''')
        
        # Load studies data
        print("Loading studies data...")
        with open(studies_file, 'r', encoding='utf-8', errors='ignore') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            
            placeholders = ','.join(['?' for _ in range(17)])
            insert_query = f"INSERT INTO studies (study_accession, pubmed_id, first_author, date, initial_sample_description, replication_sample_description, stage, number_of_individuals, broad_ancestral_category, country_of_origin, country_of_recruitment, additional_ancestry_description, ancestry_descriptor, founder_genetically_isolated_population, number_of_cases, number_of_controls, sample_description) VALUES ({placeholders})"
            
            batch = []
            count = 0
            
            for row in reader:
                # Pad row to match expected columns
                while len(row) < 17:
                    row.append('')
                
                batch.append(row[:17])
                count += 1
                
                if len(batch) >= batch_size:
                    cursor.executemany(insert_query, batch)
                    batch = []
                    if count % 1000 == 0:
                        print(f"  Processed {count} studies...")
            
            if batch:
                cursor.executemany(insert_query, batch)
            
            print(f"  Loaded {count} studies total")
        
        # Create indexes for fast querying
        print("Creating indexes...")
        
        # Indexes for gene queries
        cursor.execute('CREATE INDEX idx_reported_genes ON associations(reported_genes)')
        cursor.execute('CREATE INDEX idx_mapped_gene ON associations(mapped_gene)')
        cursor.execute('CREATE INDEX idx_disease_trait ON associations(disease_trait)')
        cursor.execute('CREATE INDEX idx_p_value ON associations(p_value)')
        cursor.execute('CREATE INDEX idx_chr_id ON associations(chr_id)')
        cursor.execute('CREATE INDEX idx_study ON associations(study)')
        cursor.execute('CREATE INDEX idx_pubmed_id ON associations(pubmed_id)')
        cursor.execute('CREATE INDEX idx_snps ON associations(snps)')
        
        # Index for studies
        cursor.execute('CREATE INDEX idx_study_accession ON studies(study_accession)')
        cursor.execute('CREATE INDEX idx_study_pubmed ON studies(pubmed_id)')
        
        # Commit changes
        conn.commit()
        print(f"✅ Database created successfully: {db_path}")
        
        # Show database stats
        cursor.execute("SELECT COUNT(*) FROM associations")
        assoc_count = cursor.fetchone()[0]
        cursor.execute("SELECT COUNT(*) FROM studies")
        study_count = cursor.fetchone()[0]
        
        print(f"📊 Database stats:")
        print(f"   Associations: {assoc_count:,}")
        print(f"   Studies: {study_count:,}")
        print(f"   Database size: {os.path.getsize(db_path) / (1024*1024):.1f} MB")
        
        return True
        
    except Exception as e:
        print(f"❌ Error creating database: {e}")
        return False
    
    finally:
        conn.close()

def verify_database(db_path: str = "local_dbs/gwas.db"):
    """Verify the database was created correctly."""
    if not os.path.exists(db_path):
        print(f"❌ Database not found: {db_path}")
        return False
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        # Test basic queries
        cursor.execute("SELECT COUNT(*) FROM associations")
        assoc_count = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(*) FROM studies") 
        study_count = cursor.fetchone()[0]
        
        # Test gene search
        cursor.execute("SELECT COUNT(*) FROM associations WHERE reported_genes LIKE '%APOE%' OR mapped_gene LIKE '%APOE%'")
        apoe_count = cursor.fetchone()[0]
        
        print(f"✅ Database verification passed:")
        print(f"   Total associations: {assoc_count:,}")
        print(f"   Total studies: {study_count:,}")
        print(f"   APOE associations: {apoe_count:,}")
        
        return True
        
    except Exception as e:
        print(f"❌ Database verification failed: {e}")
        return False
    
    finally:
        conn.close()

if __name__ == "__main__":
    print("🔄 Setting up GWAS SQLite database...")
    print("This may take a few minutes for large datasets...")
    
    # Create database
    success = create_gwas_database()
    
    if success:
        # Verify database
        verify_database()
        print("\n🎉 Setup complete! You can now use the optimized query engine.")
    else:
        print("\n❌ Setup failed!")
        sys.exit(1) 