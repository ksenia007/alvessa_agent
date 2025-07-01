import pandas as pd
import os
import numpy as np
from typing import Dict, List, Optional

from tools.tf_idf import fps_tfidf

import warnings
warnings.filterwarnings("ignore")


def query_gene_associations(
    gene_symbol: str, 
    db_path: str = "local_dbs",
    p_value_threshold: float = 5e-8,
    min_associations_per_study: int = 1,
    top_studies_by_risk: int = 10,
    top_studies_by_significance: int = 10,
    fps_disease_traits: Optional[int] = None,
) -> Dict:
    """Query gene associations and return summary results with FPS applied to disease traits."""
    engine = GWASQueryEngine(db_path)
    results = engine.query_gene(
        gene_symbol, p_value_threshold, min_associations_per_study,
        top_studies_by_risk, top_studies_by_significance, fps_disease_traits
    )
    
    all_results = {
        k: results[k] for k in [
            "gene", 
            "found", 
            "total_associations", 
            "total_significant_associations",
            "total_studies_analyzed", 
            "p_value_threshold", 
            "support_types",
            "support_values",
            "summary_by_high_risk_alleles", 
            "summary_by_significance"
        ] if k in results
    }
    return all_results


class GWASQueryEngine:
    """Query GWAS catalog data for gene associations."""
    
    def __init__(self, db_path: str = "local_dbs"):
        self.associations_file = os.path.join(db_path, "gwas_catalogue_association.tsv")
        self._associations_df = None
    
    @property
    def associations_df(self):
        if self._associations_df is None:
            self._associations_df = pd.read_csv(self.associations_file, sep='\t', low_memory=False)
        return self._associations_df
    
    def query_gene(self, gene_symbol: str, p_threshold: float = 5e-8, min_assoc: int = 1,
                   top_risk: int = 10, top_sig: int = 10, fps_traits: Optional[int] = None) -> Dict:
        """Query GWAS associations for a gene."""
        gene_symbol = gene_symbol.upper().strip()
        
        # Find gene matches
        matches = pd.concat([
            self.associations_df[self.associations_df['REPORTED GENE(S)'].str.contains(
                gene_symbol, case=False, na=False, regex=False)],
            self.associations_df[self.associations_df['MAPPED_GENE'].str.contains(
                gene_symbol, case=False, na=False, regex=False)]
        ]).drop_duplicates()
        
        empty_summary = {"related_genes": [], "high_risk_snps": [], "proteins": [], "disease_traits": []}
        
        if len(matches) == 0:
            return self._create_empty_result(gene_symbol, p_threshold, empty_summary)
        
        # Filter by p-value
        significant = self._filter_by_pvalue(matches, p_threshold)
        if len(significant) == 0:
            return self._create_empty_result(gene_symbol, p_threshold, empty_summary, 
                                           len(matches), has_matches=True)
        
        # Group by study and extract summaries
        studies = []
        for pubmed_id, group in significant.groupby('PUBMEDID'):
            if not pd.isna(pubmed_id) and len(group) >= min_assoc:
                studies.append(self._extract_study_summary(pubmed_id, group))
        
        # Sort studies by different criteria
        risk_studies = sorted(studies, key=lambda x: (len(x['risk_alleles']), x['max_risk']), reverse=True)[:top_risk]
        sig_studies = sorted(studies, key=lambda x: (x['sig_count'], -x['best_log_p']), reverse=False)[:top_sig]
        
        return {
            "gene": gene_symbol,
            "found": True,
            "total_associations": len(matches),
            "total_significant_associations": len(significant),
            "total_studies_analyzed": len(studies),
            "p_value_threshold": p_threshold,
            "summary_by_high_risk_alleles": self._generate_summary(risk_studies, True, fps_traits),
            "summary_by_significance": self._generate_summary(sig_studies, False, fps_traits),
            "studies_by_high_risk_alleles": risk_studies,
            "studies_by_significance": sig_studies,
        }
    
    def _create_empty_result(self, gene: str, p_thresh: float, summary: Dict, 
                           total: int = 0, has_matches: bool = False) -> Dict:
        """Create empty result structure."""
        result = {
            "gene": gene, "found": has_matches, "total_associations": total,
            "p_value_threshold": p_thresh, "summary_by_high_risk_alleles": summary.copy(),
            "summary_by_significance": summary.copy()
        }
        if has_matches:
            result["total_significant_associations"] = 0
        return result
    
    def _filter_by_pvalue(self, df: pd.DataFrame, threshold: float) -> pd.DataFrame:
        """Filter by p-value threshold."""
        def is_significant(p_val):
            if pd.isna(p_val):
                return False
            try:
                p_str = str(p_val)
                return float(p_str) < threshold if ('E-' in p_str or 'e-' in p_str or 
                                                  isinstance(p_val, (int, float))) else False
            except (ValueError, TypeError):
                return False
        
        return df[df['P-VALUE'].apply(is_significant)].copy().sort_values('P-VALUE')
    
    def _extract_study_summary(self, pubmed_id: str, group_df: pd.DataFrame) -> Dict:
        """Extract study summary."""
        # Extract genes
        genes = set()
        for _, row in group_df.iterrows():
            for col in ['REPORTED GENE(S)', 'MAPPED_GENE']:
                val = self._safe_get(row, col)
                if val:
                    for gene in str(val).split(','):
                        gene = gene.strip().upper()
                        if gene and gene not in ['NR', 'INTERGENIC', '']:
                            genes.add(gene)
        
        # Extract traits
        traits = {str(self._safe_get(row, 'DISEASE/TRAIT')).strip() 
                 for _, row in group_df.iterrows() 
                 if self._safe_get(row, 'DISEASE/TRAIT')}
        
        # Extract risk alleles
        risk_alleles = []
        max_risk = 0
        for _, row in group_df.iterrows():
            risk_allele = self._safe_get(row, 'STRONGEST SNP-RISK ALLELE')
            or_beta = self._safe_get(row, 'OR or BETA')
            if risk_allele and or_beta:
                try:
                    beta_val = float(str(or_beta).replace('>', '').replace('<', ''))
                    risk_score = abs(beta_val) if beta_val <= 1 else beta_val
                    max_risk = max(max_risk, risk_score)
                    risk_alleles.append({
                        "risk_allele": risk_allele, "snp": self._safe_get(row, 'SNPS'),
                        "p_value": self._safe_get(row, 'P-VALUE'), "or_beta": or_beta,
                        "risk_score": risk_score, "disease_trait": self._safe_get(row, 'DISEASE/TRAIT'),
                        "mapped_gene": self._safe_get(row, 'MAPPED_GENE')
                    })
                except (ValueError, TypeError):
                    continue
        
        # Calculate significance metrics
        p_values = []
        for _, row in group_df.iterrows():
            p_val = self._safe_get(row, 'P-VALUE')
            if p_val:
                try:
                    p_float = float(p_val) if isinstance(p_val, str) and ('E-' in p_val or 'e-' in p_val) else float(p_val)
                    p_values.append(p_float)
                except (ValueError, TypeError):
                    continue
        
        best_p = min(p_values) if p_values else 1.0
        sig_count = sum(1 for p in p_values if p < 5e-8)
        
        return {
            "pubmed_id": str(pubmed_id), 
            "related_genes": sorted(genes),
            "disease_traits": sorted(traits), 
            "risk_alleles": sorted(risk_alleles, key=lambda x: x['risk_score'], reverse=True),
            "max_risk": max_risk, 
            "best_pvalue": best_p, 
            "best_log_p": -np.log10(best_p) if best_p > 0 else 100,
            "sig_count": sig_count
        }
    
    def _generate_summary(self, studies: List[Dict], sort_by_risk: bool = True, fps_traits: Optional[int] = None) -> Dict:
        """Generate summary from studies."""
        gene_scores = {}
        trait_scores = {}
        snp_scores = {}
        proteins = set()
        all_risk_alleles = []
        
        for study in studies:
            study_risk = study.get('max_risk', 0)
            study_pval = study.get('best_pvalue', 1.0)
            
            # Collect genes
            for gene in study['related_genes']:
                if gene not in gene_scores:
                    gene_scores[gene] = {'risk_score': 0, 'best_pvalue': 1.0}
                gene_scores[gene]['risk_score'] = max(gene_scores[gene]['risk_score'], study_risk)
                gene_scores[gene]['best_pvalue'] = min(gene_scores[gene]['best_pvalue'], study_pval)
            
            # Collect traits and extract proteins
            for trait in study['disease_traits']:
                if trait not in trait_scores:
                    trait_scores[trait] = {'risk_score': 0, 'best_pvalue': 1.0}
                trait_scores[trait]['risk_score'] = max(trait_scores[trait]['risk_score'], study_risk)
                trait_scores[trait]['best_pvalue'] = min(trait_scores[trait]['best_pvalue'], study_pval)
                
                # Extract proteins
                trait_lower = trait.lower()
                protein_keywords = ["protein level", "protein measurement", "protein concentration", "serum protein"]
                if any(kw in trait_lower for kw in protein_keywords):
                    protein_name = trait
                    for replacement in [" protein levels", " protein level", " protein measurement", 
                                      " protein concentration", " serum protein", "serum "]:
                        protein_name = protein_name.replace(replacement, "")
                    if protein_name.strip():
                        proteins.add(protein_name.strip())
            
            # Collect SNPs
            for allele in study['risk_alleles']:
                all_risk_alleles.append(allele)
                snp = allele.get('snp')
                if snp:
                    risk_score = allele.get('risk_score', 0)
                    p_val = self._parse_pvalue(allele.get('p_value', 1.0))
                    
                    if snp not in snp_scores:
                        snp_scores[snp] = {'risk_score': 0, 'best_pvalue': 1.0}
                    snp_scores[snp]['risk_score'] = max(snp_scores[snp]['risk_score'], risk_score)
                    snp_scores[snp]['best_pvalue'] = min(snp_scores[snp]['best_pvalue'], p_val)
        
        # Sort by appropriate criteria
        sort_key = (lambda x: gene_scores[x]['risk_score']) if sort_by_risk else (lambda x: gene_scores[x]['best_pvalue'])
        sorted_genes = sorted(gene_scores.keys(), key=sort_key, reverse=sort_by_risk)
        
        sort_key = (lambda x: trait_scores[x]['risk_score']) if sort_by_risk else (lambda x: trait_scores[x]['best_pvalue'])
        sorted_traits = sorted(trait_scores.keys(), key=sort_key, reverse=sort_by_risk)
        
        sort_key = (lambda x: snp_scores[x]['risk_score']) if sort_by_risk else (lambda x: snp_scores[x]['best_pvalue'])
        sorted_snps = sorted(snp_scores.keys(), key=sort_key, reverse=sort_by_risk)
        
        # Apply FPS to traits if requested
        if fps_traits and len(sorted_traits) > fps_traits:
            fps_indices = fps_tfidf(sorted_traits, fps_traits)
            sorted_traits = [sorted_traits[i] for i in fps_indices]
        
        return {
            "gwas_related_genes": sorted_genes,
            "gwas_high_risk_snps": sorted_snps,
            "gwas_affected_protein_levels": sorted(proteins),
            "gwas_disease_traits": sorted_traits
        }
    
    def _parse_pvalue(self, p_val) -> float:
        """Parse p-value to float."""
        try:
            if isinstance(p_val, str) and ('E-' in p_val or 'e-' in p_val):
                return float(p_val)
            return float(p_val) if p_val else 1.0
        except (ValueError, TypeError):
            return 1.0
    
    def _safe_get(self, row, column: str):
        """Safely get value from row."""
        try:
            value = row[column]
            return None if pd.isna(value) or value == '' or str(value).lower() == 'nan' else value
        except (KeyError, IndexError):
            return None


if __name__ == "__main__":
    gene = "NBR2"
    results = query_gene_associations(gene, fps_disease_traits=20)
    
    print(f"Results for gene {gene}:")
    print(f"Found: {results['found']}")
    print(f"Total associations: {results['total_associations']}")
    print(f"Significant: {results['total_significant_associations']}")
    print(f"Studies: {results['total_studies_analyzed']}")
    print(f"Studies by high-risk alleles: ")
    for study in results['studies_by_high_risk_alleles']:
        print(f"  {study['pubmed_id']}: {study['related_genes']}, {study['max_risk']}, {study['best_pvalue']}, {study['risk_alleles']}")
    
    print(f"--------------------------------")
    print(f"Studies by significance: ")
    for study in results['studies_by_significance']:
        print(f"  {study['pubmed_id']}: {study['related_genes']}, {study['best_pvalue']}, {study['risk_alleles']}")
    
    for summary_type in ["high_risk_alleles", "significance"]:
        summary = results[f'summary_by_{summary_type}']
        print(f"\nSummary by {summary_type.replace('_', ' ').title()}:")
        print(f"  Genes ({len(summary['related_genes'])}): {summary['related_genes']}")
        print(f"  SNPs ({len(summary['high_risk_snps'])}): {summary['high_risk_snps']}")
        print(f"  Proteins ({len(summary['proteins'])}): {summary['proteins']}")
        print(f"  Traits: {summary['disease_traits']}")
    