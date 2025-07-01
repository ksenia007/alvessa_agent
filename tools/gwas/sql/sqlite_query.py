#!/usr/bin/env python3
"""
Optimized GWAS query engine using SQLite database.
Much faster and more memory-efficient than loading pandas DataFrames.
"""

import sqlite3
import json
import os
from typing import Dict, List, Optional, Union

class GWASSQLiteQueryEngine:
    """
    Fast GWAS query engine using SQLite database with proper indexing.
    """
    
    def __init__(self, db_path: str = "local_dbs/gwas.db"):
        """
        Initialize the GWAS SQLite query engine.
        
        Args:
            db_path: Path to the SQLite database file
        """
        self.db_path = db_path
        
        if not os.path.exists(db_path):
            raise FileNotFoundError(
                f"SQLite database not found at {db_path}. "
                f"Please run 'python tools/gwas/setup_db.py' first to create the database."
            )
    
    def _get_connection(self):
        """Get a database connection."""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row  # Enable column access by name
        return conn
    
    def query_gene(self, gene_symbol: str, include_study_details: bool = True, limit: int = None) -> Dict:
        """
        Query GWAS associations for a specific gene using fast SQL queries.
        
        Args:
            gene_symbol: Gene symbol to search for (e.g., 'APOE', 'BRCA1')
            include_study_details: Whether to include detailed study information
            limit: Maximum number of associations to return (None for all)
            
        Returns:
            Dictionary containing structured results
        """
        gene_symbol = gene_symbol.upper().strip()
        
        conn = self._get_connection()
        cursor = conn.cursor()
        
        try:
            # Build the query
            base_query = """
                SELECT * FROM associations 
                WHERE (reported_genes LIKE ? OR mapped_gene LIKE ?)
                AND (reported_genes IS NOT NULL OR mapped_gene IS NOT NULL)
            """
            
            params = [f'%{gene_symbol}%', f'%{gene_symbol}%']
            
            if limit:
                base_query += " LIMIT ?"
                params.append(limit)
            
            # Execute query
            cursor.execute(base_query, params)
            rows = cursor.fetchall()
            
            if not rows:
                return {
                    "gene": gene_symbol,
                    "found": False,
                    "total_associations": 0,
                    "associations": [],
                    "summary": {
                        "unique_traits": [],
                        "unique_studies": [],
                        "significant_associations": 0,
                        "trait_count": 0,
                        "study_count": 0
                    }
                }
            
            # Process associations
            associations = []
            traits = set()
            studies = set()
            significant_count = 0
            
            for row in rows:
                # Convert row to dictionary
                association = {
                    "study_id": row['study'],
                    "pubmed_id": row['pubmed_id'],
                    "first_author": row['first_author'],
                    "date": row['date'],
                    "journal": row['journal'],
                    "link": row['link'],
                    "disease_trait": row['disease_trait'],
                    "reported_genes": row['reported_genes'],
                    "mapped_gene": row['mapped_gene'],
                    "chromosome": row['chr_id'],
                    "position": row['chr_pos'],
                    "snp": row['snps'],
                    "risk_allele": row['strongest_snp_risk_allele'],
                    "p_value": row['p_value'],
                    "p_value_text": row['p_value_text'],
                    "or_beta": row['or_beta'],
                    "confidence_interval": row['ci_95_text'],
                    "risk_allele_frequency": row['risk_allele_frequency'],
                    "sample_size_initial": row['initial_sample_size'],
                    "sample_size_replication": row['replication_sample_size'],
                    "platform": row['platform']
                }
                
                # Add study details if requested
                if include_study_details and row['study']:
                    study_details = self._get_study_details(row['study'])
                    association["study_details"] = study_details
                
                associations.append(association)
                
                # Collect summary data
                if row['disease_trait']:
                    traits.add(row['disease_trait'])
                if row['study']:
                    studies.add(row['study'])
                if self._is_significant_p_value(row['p_value']):
                    significant_count += 1
            
            # Build result
            result = {
                "gene": gene_symbol,
                "found": True,
                "total_associations": len(associations),
                "associations": associations,
                "summary": {
                    "unique_traits": list(traits),
                    "unique_studies": list(studies),
                    "significant_associations": significant_count,
                    "trait_count": len(traits),
                    "study_count": len(studies)
                }
            }
            
            return result
            
        finally:
            conn.close()
    
    def _get_study_details(self, study_id: str) -> Optional[Dict]:
        """Get detailed information about a study using SQL."""
        if not study_id:
            return None
        
        conn = self._get_connection()
        cursor = conn.cursor()
        
        try:
            cursor.execute(
                "SELECT * FROM studies WHERE study_accession = ?",
                (study_id,)
            )
            row = cursor.fetchone()
            
            if not row:
                return None
            
            return {
                "study_accession": row['study_accession'],
                "pubmed_id": row['pubmed_id'],
                "first_author": row['first_author'],
                "date": row['date'],
                "initial_sample_description": row['initial_sample_description'],
                "replication_sample_description": row['replication_sample_description'],
                "number_of_individuals": row['number_of_individuals'],
                "ancestral_category": row['broad_ancestral_category'],
                "country_of_origin": row['country_of_origin'],
                "number_of_cases": row['number_of_cases'],
                "number_of_controls": row['number_of_controls']
            }
            
        finally:
            conn.close()
    
    def search_by_trait(self, trait: str, limit: int = 100) -> Dict:
        """
        Search GWAS associations by disease/trait using fast SQL.
        
        Args:
            trait: Disease or trait to search for
            limit: Maximum number of results to return
            
        Returns:
            Dictionary containing structured results
        """
        conn = self._get_connection()
        cursor = conn.cursor()
        
        try:
            cursor.execute("""
                SELECT reported_genes, mapped_gene, disease_trait, p_value, 
                       snps, strongest_snp_risk_allele, pubmed_id, study
                FROM associations 
                WHERE disease_trait LIKE ? 
                AND disease_trait IS NOT NULL
                LIMIT ?
            """, (f'%{trait}%', limit))
            
            rows = cursor.fetchall()
            
            if not rows:
                return {
                    "trait": trait,
                    "found": False,
                    "total_associations": 0,
                    "associations": []
                }
            
            associations = []
            for row in rows:
                association = {
                    "reported_genes": row['reported_genes'],
                    "mapped_gene": row['mapped_gene'],
                    "disease_trait": row['disease_trait'],
                    "p_value": row['p_value'],
                    "snp": row['snps'],
                    "risk_allele": row['strongest_snp_risk_allele'],
                    "pubmed_id": row['pubmed_id'],
                    "study_id": row['study']
                }
                associations.append(association)
            
            return {
                "trait": trait,
                "found": True,
                "total_associations": len(associations),
                "associations": associations
            }
            
        finally:
            conn.close()
    
    def search_by_chromosome(self, chromosome: str, start_pos: int = None, 
                           end_pos: int = None, limit: int = 1000) -> Dict:
        """
        Search GWAS associations by chromosomal region.
        
        Args:
            chromosome: Chromosome (e.g., '1', '2', 'X')
            start_pos: Start position (optional)
            end_pos: End position (optional) 
            limit: Maximum results to return
        """
        conn = self._get_connection()
        cursor = conn.cursor()
        
        try:
            base_query = "SELECT * FROM associations WHERE chr_id = ?"
            params = [str(chromosome)]
            
            if start_pos is not None and end_pos is not None:
                base_query += " AND CAST(chr_pos AS INTEGER) BETWEEN ? AND ?"
                params.extend([start_pos, end_pos])
            elif start_pos is not None:
                base_query += " AND CAST(chr_pos AS INTEGER) >= ?"
                params.append(start_pos)
            elif end_pos is not None:
                base_query += " AND CAST(chr_pos AS INTEGER) <= ?"
                params.append(end_pos)
            
            base_query += " LIMIT ?"
            params.append(limit)
            
            cursor.execute(base_query, params)
            rows = cursor.fetchall()
            
            associations = []
            for row in rows:
                association = {
                    "reported_genes": row['reported_genes'],
                    "mapped_gene": row['mapped_gene'],
                    "disease_trait": row['disease_trait'],
                    "chromosome": row['chr_id'],
                    "position": row['chr_pos'],
                    "p_value": row['p_value'],
                    "snp": row['snps'],
                    "risk_allele": row['strongest_snp_risk_allele']
                }
                associations.append(association)
            
            return {
                "chromosome": chromosome,
                "start_position": start_pos,
                "end_position": end_pos,
                "found": len(associations) > 0,
                "total_associations": len(associations),
                "associations": associations
            }
            
        finally:
            conn.close()
    
    def get_top_genes_by_associations(self, limit: int = 50) -> Dict:
        """Get genes with the most associations."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        try:
            cursor.execute("""
                SELECT mapped_gene, COUNT(*) as association_count
                FROM associations 
                WHERE mapped_gene IS NOT NULL AND mapped_gene != ''
                GROUP BY mapped_gene 
                ORDER BY association_count DESC 
                LIMIT ?
            """, (limit,))
            
            rows = cursor.fetchall()
            
            genes = []
            for row in rows:
                genes.append({
                    "gene": row['mapped_gene'],
                    "association_count": row['association_count']
                })
            
            return {
                "top_genes": genes,
                "total_genes": len(genes)
            }
            
        finally:
            conn.close()
    
    def get_database_stats(self) -> Dict:
        """Get database statistics."""
        conn = self._get_connection()
        cursor = conn.cursor()
        
        try:
            # Basic counts
            cursor.execute("SELECT COUNT(*) FROM associations")
            total_associations = cursor.fetchone()[0]
            
            cursor.execute("SELECT COUNT(*) FROM studies")
            total_studies = cursor.fetchone()[0]
            
            # Unique genes
            cursor.execute("SELECT COUNT(DISTINCT mapped_gene) FROM associations WHERE mapped_gene IS NOT NULL")
            unique_genes = cursor.fetchone()[0]
            
            # Unique traits
            cursor.execute("SELECT COUNT(DISTINCT disease_trait) FROM associations WHERE disease_trait IS NOT NULL")
            unique_traits = cursor.fetchone()[0]
            
            # Significant associations
            cursor.execute("SELECT COUNT(*) FROM associations WHERE CAST(p_value AS REAL) < 5e-8")
            significant_associations = cursor.fetchone()[0]
            
            return {
                "total_associations": total_associations,
                "total_studies": total_studies,
                "unique_genes": unique_genes,
                "unique_traits": unique_traits,
                "significant_associations": significant_associations,
                "database_size_mb": os.path.getsize(self.db_path) / (1024 * 1024)
            }
            
        finally:
            conn.close()
    
    def _is_significant_p_value(self, p_value) -> bool:
        """Check if a p-value is genome-wide significant (< 5e-8)."""
        if not p_value:
            return False
        
        try:
            if isinstance(p_value, str):
                if 'e-' in p_value.lower() or 'E-' in p_value:
                    p_val = float(p_value)
                else:
                    p_val = float(p_value)
            else:
                p_val = float(p_value)
            
            return p_val < 5e-8
        except (ValueError, TypeError):
            return False


def query_gene_associations_sqlite(gene_symbol: str, db_path: str = "local_dbs/gwas.db", 
                                 include_study_details: bool = True, 
                                 output_format: str = "dict",
                                 limit: int = None) -> Union[Dict, str]:
    """
    Convenience function to query gene associations using SQLite.
    
    Args:
        gene_symbol: Gene symbol to search for
        db_path: Path to SQLite database
        include_study_details: Whether to include study details
        output_format: 'dict' or 'json'
        limit: Maximum number of associations to return
        
    Returns:
        Query results as dictionary or JSON string
    """
    engine = GWASSQLiteQueryEngine(db_path)
    results = engine.query_gene(gene_symbol, include_study_details, limit)
    
    if output_format.lower() == 'json':
        return json.dumps(results, indent=2, default=str)
    else:
        return results


# Example usage
if __name__ == "__main__":
    # Example queries using SQLite
    gene_to_query = "TP53"  # Change this to any gene of interest
    
    try:
        # Query for a specific gene
        results = query_gene_associations_sqlite(gene_to_query, limit=10)
        
        print(f"🧬 Results for gene {gene_to_query}:")
        print(f"   Found: {results['found']}")
        print(f"   Total associations: {results['total_associations']}")
        print(f"   Unique traits: {len(results['summary']['unique_traits'])}")
        print(f"   Significant associations: {results['summary']['significant_associations']}")
        
        # Print first few associations
        print(f"\n📊 First few associations:")
        for i, assoc in enumerate(results['associations'][:5]):
            print(f"   {i+1}. Trait: {assoc['disease_trait']}")
            print(f"      P-value: {assoc['p_value']}")
            print(f"      SNP: {assoc['snp']}")
            print(f"      Risk allele: {assoc['risk_allele']}")
            print()
        
        # Test database stats
        engine = GWASSQLiteQueryEngine()
        stats = engine.get_database_stats()
        print(f"📈 Database statistics:")
        for key, value in stats.items():
            if isinstance(value, float):
                print(f"   {key}: {value:.1f}")
            else:
                print(f"   {key}: {value:,}")
        
    except FileNotFoundError as e:
        print(f"❌ {e}")
        print("Run: python tools/gwas/setup_db.py")
    except Exception as e:
        print(f"❌ Error: {e}") 