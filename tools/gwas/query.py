import pandas as pd
import json
import os
from typing import Dict, List, Optional, Union


class GWASQueryEngine:
    """
    A class to query GWAS catalog data for gene associations and return structured results.
    """
    
    def __init__(self, db_path: str = "local_dbs"):
        """
        Initialize the GWAS query engine.
        
        Args:
            db_path: Path to the directory containing GWAS database files
        """
        self.db_path = db_path
        self.associations_file = os.path.join(db_path, "gwas_catalogue_association.tsv")
        self.studies_file = os.path.join(db_path, "gwas_catalogue_studies.tsv")
        
        # Load data lazily
        self._associations_df = None
        self._studies_df = None
    
    @property
    def associations_df(self):
        """Lazy load associations dataframe."""
        if self._associations_df is None:
            self._associations_df = pd.read_csv(self.associations_file, sep='\t', low_memory=False)
        return self._associations_df
    
    @property
    def studies_df(self):
        """Lazy load studies dataframe."""
        if self._studies_df is None:
            self._studies_df = pd.read_csv(self.studies_file, sep='\t', low_memory=False)
        return self._studies_df
    
    def query_gene(self, gene_symbol: str, include_study_details: bool = True) -> Dict:
        """
        Query GWAS associations for a specific gene.
        
        Args:
            gene_symbol: Gene symbol to search for (e.g., 'APOE', 'BRCA1')
            include_study_details: Whether to include detailed study information
            
        Returns:
            Dictionary containing structured results
        """
        gene_symbol = gene_symbol.upper().strip()
        
        # Search in both REPORTED GENE(S) and MAPPED_GENE columns
        reported_matches = self.associations_df[
            self.associations_df['REPORTED GENE(S)'].str.contains(
                gene_symbol, case=False, na=False, regex=False
            )
        ]
        
        mapped_matches = self.associations_df[
            self.associations_df['MAPPED_GENE'].str.contains(
                gene_symbol, case=False, na=False, regex=False
            )
        ]
        
        # Combine and remove duplicates
        all_matches = pd.concat([reported_matches, mapped_matches]).drop_duplicates()
        
        if len(all_matches) == 0:
            return {
                "gene": gene_symbol,
                "found": False,
                "total_associations": 0,
                "associations": [],
                "summary": {
                    "unique_traits": [],
                    "unique_studies": [],
                    "significant_associations": 0
                }
            }
        
        # Process associations
        associations = []
        for _, row in all_matches.iterrows():
            association = {
                "study_id": self._safe_get(row, 'STUDY'),
                "pubmed_id": self._safe_get(row, 'PUBMEDID'),
                "first_author": self._safe_get(row, 'FIRST AUTHOR'),
                "date": self._safe_get(row, 'DATE'),
                "journal": self._safe_get(row, 'JOURNAL'),
                "disease_trait": self._safe_get(row, 'DISEASE/TRAIT'),
                "reported_genes": self._safe_get(row, 'REPORTED GENE(S)'),
                "mapped_gene": self._safe_get(row, 'MAPPED_GENE'),
                "chromosome": self._safe_get(row, 'CHR_ID'),
                "position": self._safe_get(row, 'CHR_POS'),
                "snp": self._safe_get(row, 'SNPS'),
                "risk_allele": self._safe_get(row, 'STRONGEST SNP-RISK ALLELE'),
                "p_value": self._safe_get(row, 'P-VALUE'),
                "p_value_text": self._safe_get(row, 'P-VALUE (TEXT)'),
                "or_beta": self._safe_get(row, 'OR or BETA'),
                "confidence_interval": self._safe_get(row, '95% CI (TEXT)'),
                "risk_allele_frequency": self._safe_get(row, 'RISK ALLELE FREQUENCY'),
                "sample_size_initial": self._safe_get(row, 'INITIAL SAMPLE SIZE'),
                "sample_size_replication": self._safe_get(row, 'REPLICATION SAMPLE SIZE'),
                "platform": self._safe_get(row, 'PLATFORM [SNPS PASSING QC]')
            }
            
            # Add study details if requested
            if include_study_details:
                study_details = self._get_study_details(association["study_id"])
                association["study_details"] = study_details

            associations.append(association)
        
        import ipdb; ipdb.set_trace()
        
        # Generate summary statistics
        unique_traits = all_matches['DISEASE/TRAIT'].dropna().unique().tolist()
        unique_studies = all_matches['STUDY'].dropna().unique().tolist()
        
        # Count significant associations (p < 5e-8, standard GWAS threshold)
        significant_count = 0
        for _, row in all_matches.iterrows():
            p_val = self._safe_get(row, 'P-VALUE')
            if p_val and self._is_significant_p_value(p_val):
                significant_count += 1
        
        result = {
            "gene": gene_symbol,
            "found": True,
            "total_associations": len(all_matches),
            "associations": associations,
            "summary": {
                "unique_traits": unique_traits,
                "unique_studies": unique_studies,
                "significant_associations": significant_count,
                "trait_count": len(unique_traits),
                "study_count": len(unique_studies)
            }
        }
        
        return result
    
    def _get_study_details(self, study_id: str) -> Optional[Dict]:
        """Get detailed information about a study."""
        if not study_id or pd.isna(study_id):
            return None
            
        study_match = self.studies_df[self.studies_df['STUDY ACCESSION'] == study_id]
        
        if len(study_match) == 0:
            return None
        
        study_row = study_match.iloc[0]
        return {
            "study_accession": self._safe_get(study_row, 'STUDY ACCESSION'),
            "pubmed_id": self._safe_get(study_row, 'PUBMED ID'),
            "first_author": self._safe_get(study_row, 'FIRST AUTHOR'),
            "date": self._safe_get(study_row, 'DATE'),
            "initial_sample_description": self._safe_get(study_row, 'INITIAL SAMPLE DESCRIPTION'),
            "replication_sample_description": self._safe_get(study_row, 'REPLICATION SAMPLE DESCRIPTION'),
            "number_of_individuals": self._safe_get(study_row, 'NUMBER OF INDIVIDUALS'),
            "ancestral_category": self._safe_get(study_row, 'BROAD ANCESTRAL CATEGORY'),
            "country_of_origin": self._safe_get(study_row, 'COUNTRY OF ORIGIN'),
            "number_of_cases": self._safe_get(study_row, 'NUMBER OF CASES'),
            "number_of_controls": self._safe_get(study_row, 'NUMBER OF CONTROLS')
        }
    
    def _safe_get(self, row, column: str):
        """Safely get a value from a row, handling NaN values."""
        try:
            value = row[column]
            if pd.isna(value) or value == '' or str(value).lower() == 'nan':
                return None
            return value
        except (KeyError, IndexError):
            return None
    
    def _is_significant_p_value(self, p_value) -> bool:
        """Check if a p-value is genome-wide significant (< 5e-8)."""
        try:
            if isinstance(p_value, str):
                # Handle scientific notation in text
                if 'E-' in p_value or 'e-' in p_value:
                    p_val = float(p_value)
                else:
                    return False
            else:
                p_val = float(p_value)
            
            return p_val < 5e-8
        except (ValueError, TypeError):
            return False
    
    def search_by_trait(self, trait: str, limit: int = 100) -> Dict:
        """
        Search GWAS associations by disease/trait.
        
        Args:
            trait: Disease or trait to search for
            limit: Maximum number of results to return
            
        Returns:
            Dictionary containing structured results
        """
        trait_matches = self.associations_df[
            self.associations_df['DISEASE/TRAIT'].str.contains(
                trait, case=False, na=False, regex=False
            )
        ].head(limit)
        
        if len(trait_matches) == 0:
            return {
                "trait": trait,
                "found": False,
                "total_associations": 0,
                "associations": []
            }
        
        associations = []
        for _, row in trait_matches.iterrows():
            association = {
                "reported_genes": self._safe_get(row, 'REPORTED GENE(S)'),
                "mapped_gene": self._safe_get(row, 'MAPPED_GENE'),
                "disease_trait": self._safe_get(row, 'DISEASE/TRAIT'),
                "p_value": self._safe_get(row, 'P-VALUE'),
                "snp": self._safe_get(row, 'SNPS'),
                "risk_allele": self._safe_get(row, 'STRONGEST SNP-RISK ALLELE'),
                "pubmed_id": self._safe_get(row, 'PUBMEDID'),
                "study_id": self._safe_get(row, 'STUDY')
            }
            associations.append(association)
        
        return {
            "trait": trait,
            "found": True,
            "total_associations": len(trait_matches),
            "associations": associations
        }


def query_gene_associations(gene_symbol: str, db_path: str = "local_dbs", 
                          include_study_details: bool = True, 
                          output_format: str = "dict") -> Union[Dict, str]:
    """
    Convenience function to query gene associations.
    
    Args:
        gene_symbol: Gene symbol to search for
        db_path: Path to database files
        include_study_details: Whether to include study details
        output_format: 'dict' or 'json'
        
    Returns:
        Query results as dictionary or JSON string
    """
    engine = GWASQueryEngine(db_path)
    results = engine.query_gene(gene_symbol, include_study_details)
    
    if output_format.lower() == 'json':
        return json.dumps(results, indent=2, default=str)
    else:
        return results


# Example usage
if __name__ == "__main__":
    # Example queries
    gene_to_query = "TP53"  # Change this to any gene of interest
    
    # Query for a specific gene
    results = query_gene_associations(gene_to_query, output_format="dict")
    
    print(f"Results for gene {gene_to_query}:")
    print(f"Found: {results['found']}")
    print(f"Total associations: {results['total_associations']}")
    print(f"Unique traits: {len(results['summary']['unique_traits'])}")
    print(f"Significant associations: {results['summary']['significant_associations']}")
    
    # Print first few associations
    for i, assoc in enumerate(results['associations'][:10]):
        print(f"\nAssociation {i+1}:")
        print(f"  Trait: {assoc['disease_trait']}")
        print(f"  P-value: {assoc['p_value']}")
        print(f"  SNP: {assoc['snp']}")
        print(f"  Risk allele: {assoc['risk_allele']}")
