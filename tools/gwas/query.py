import pandas as pd
import os
import numpy as np
from typing import Dict, List, Optional, Tuple, Union
from tools.tf_idf import fps_tfidf
import warnings
warnings.filterwarnings("ignore")


def query_gene_associations(gene_symbol: str, db_path: str = "local_dbs", 
                          p_value_threshold: float = 5e-8, min_associations_per_study: int = 1,
                          top_studies_by_risk: int = 10, top_studies_by_significance: int = 10,
                          fps_disease_traits: Optional[int] = None) -> Dict:
    engine = GWASQueryEngine(db_path)
    results = engine.query_gene(gene_symbol, p_value_threshold, min_associations_per_study,
                               top_studies_by_risk, top_studies_by_significance, fps_disease_traits)
    return engine._format_gene_results(results)


def query_trait_associations(trait_term: str, db_path: str = "local_dbs",
                           p_value_threshold: float = 5e-8, min_associations_per_study: int = 1,
                           top_studies_by_risk: int = 10, top_studies_by_significance: int = 10,
                           fps_genes: Optional[int] = None, exact_match: bool = False) -> Dict:
    engine = GWASQueryEngine(db_path)
    results = engine.query_trait(trait_term, p_value_threshold, min_associations_per_study,
                               top_studies_by_risk, top_studies_by_significance, fps_genes, exact_match)
    return engine._format_trait_results(results)


class GWASQueryEngine:
    GENE_COLUMNS = ['REPORTED GENE(S)', 'MAPPED_GENE']
    PROTEIN_KEYWORDS = ["protein level", "protein measurement", "protein concentration", "serum protein"]
    PROTEIN_REPLACEMENTS = [" protein levels", " protein level", " protein measurement", 
                           " protein concentration", " serum protein", "serum "]
    EXCLUDED_GENES = {'NR', 'INTERGENIC', ''}
    VARIANT_COLUMNS = ['SNPS', 'CHR_ID', 'CHR_POS', 'REGION', 'CONTEXT', 'INTERGENIC', 
                      'STRONGEST SNP-RISK ALLELE', 'RISK ALLELE FREQUENCY', 'UPSTREAM_GENE_DISTANCE', 
                      'DOWNSTREAM_GENE_DISTANCE', 'UPSTREAM_GENE_ID', 'DOWNSTREAM_GENE_ID']
    
    # Variant context categories for analysis
    CODING_CONTEXTS = {'missense_variant', 'synonymous_variant', 'stop_gained', 'frameshift_variant'}
    UTR_CONTEXTS = {'3_prime_UTR_variant', '5_prime_UTR_variant'}
    SPLICE_CONTEXTS = {'splice_region_variant', 'splice_donor_variant', 'splice_acceptor_variant'}
    REGULATORY_CONTEXTS = {'regulatory_region_variant', 'TF_binding_site_variant'}
    
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
        gene_symbol = gene_symbol.upper().strip()
        matches = self._find_gene_matches(gene_symbol)
        return self._process_matches(matches, gene_symbol, p_threshold, min_assoc, 
                                   top_risk, top_sig, fps_traits, is_gene_query=True)
    
    def query_trait(self, trait_term: str, p_threshold: float = 5e-8, min_assoc: int = 1,
                   top_risk: int = 10, top_sig: int = 10, fps_genes: Optional[int] = None, 
                   exact_match: bool = False) -> Dict:
        trait_term = trait_term.strip()
        matches = self._find_trait_matches(trait_term, exact_match)
        return self._process_matches(matches, trait_term, p_threshold, min_assoc,
                                   top_risk, top_sig, fps_genes, is_gene_query=False)
    
    def _find_gene_matches(self, gene_symbol: str) -> pd.DataFrame:
        matches = []
        for col in self.GENE_COLUMNS:
            mask = self.associations_df[col].str.upper().str.strip() == gene_symbol
            matches.append(self.associations_df[mask])
        return pd.concat(matches).drop_duplicates() if matches else pd.DataFrame()
    
    def _find_trait_matches(self, trait_term: str, exact_match: bool) -> pd.DataFrame:
        if exact_match:
            return self.associations_df[self.associations_df['DISEASE/TRAIT'].str.strip() == trait_term]
        return self.associations_df[
            self.associations_df['DISEASE/TRAIT'].str.contains(trait_term, case=False, na=False, regex=False)
        ]
    
    def _process_matches(self, matches: pd.DataFrame, search_term: str, p_threshold: float,
                        min_assoc: int, top_risk: int, top_sig: int, fps_param: Optional[int],
                        is_gene_query: bool) -> Dict:
        if matches.empty:
            return self._create_empty_result(search_term, p_threshold, is_gene_query)
        
        significant = self._filter_by_pvalue(matches, p_threshold)
        if significant.empty:
            return self._create_empty_result(search_term, p_threshold, is_gene_query, 
                                           len(matches), has_matches=True)
        
        studies = self._extract_studies(significant, min_assoc)
        risk_studies, sig_studies = self._sort_studies(studies, top_risk, top_sig)
        
        return {
            "gene" if is_gene_query else "trait": search_term,
            "found": True,
            "total_associations": len(matches),
            "total_significant_associations": len(significant),
            "total_studies_analyzed": len(studies),
            "p_value_threshold": p_threshold,
            "summary_by_high_risk_alleles": self._generate_summary(risk_studies, True, fps_param, is_gene_query),
            "summary_by_significance": self._generate_summary(sig_studies, False, fps_param, is_gene_query),
            "studies_by_high_risk_alleles": risk_studies,
            "studies_by_significance": sig_studies,
        }
    
    def _create_empty_result(self, search_term: str, p_threshold: float, is_gene_query: bool,
                           total: int = 0, has_matches: bool = False) -> Dict:
        key_prefix = "" if is_gene_query else "gwas_"
        empty_summary = {
            f"{key_prefix}related_genes" if is_gene_query else f"{key_prefix}related_genes": [],
            f"{key_prefix}high_risk_snps": [],
            f"proteins" if is_gene_query else f"{key_prefix}affected_protein_levels": [],
            f"{key_prefix}disease_traits": []
        }
        
        result = {
            "gene" if is_gene_query else "trait": search_term,
            "found": has_matches,
            "total_associations": total,
            "p_value_threshold": p_threshold,
            "summary_by_high_risk_alleles": empty_summary.copy(),
            "summary_by_significance": empty_summary.copy()
        }
        
        if has_matches:
            result["total_significant_associations"] = 0
        return result
    
    def _filter_by_pvalue(self, df: pd.DataFrame, threshold: float) -> pd.DataFrame:
        def is_significant(p_val):
            if pd.isna(p_val):
                return False
            try:
                p_str = str(p_val)
                return (float(p_str) < threshold if ('E-' in p_str or 'e-' in p_str or 
                       isinstance(p_val, (int, float))) else False)
            except (ValueError, TypeError):
                return False
        
        return df[df['P-VALUE'].apply(is_significant)].copy().sort_values('P-VALUE')
    
    def _extract_studies(self, significant: pd.DataFrame, min_assoc: int) -> List[Dict]:
        studies = []
        for pubmed_id, group in significant.groupby('PUBMEDID'):
            if not pd.isna(pubmed_id) and len(group) >= min_assoc:
                studies.append(self._extract_study_summary(pubmed_id, group))
        return studies
    
    def _sort_studies(self, studies: List[Dict], top_risk: int, top_sig: int) -> Tuple[List[Dict], List[Dict]]:
        risk_studies = sorted(studies, key=lambda x: (len(x['risk_alleles']), x['max_risk']), reverse=True)[:top_risk]
        sig_studies = sorted(studies, key=lambda x: (x['sig_count'], -x['best_log_p']), reverse=False)[:top_sig]
        return risk_studies, sig_studies
    
    def _extract_study_summary(self, pubmed_id: str, group_df: pd.DataFrame) -> Dict:
        genes = self._extract_genes(group_df)
        traits = self._extract_traits(group_df)
        risk_alleles, max_risk = self._extract_risk_alleles(group_df)
        p_values = self._extract_p_values(group_df)
        
        best_p = min(p_values) if p_values else 1.0
        return {
            "pubmed_id": str(pubmed_id),
            "related_genes": sorted(genes),
            "disease_traits": sorted(traits),
            "risk_alleles": sorted(risk_alleles, key=lambda x: x['risk_score'], reverse=True),
            "max_risk": max_risk,
            "best_pvalue": best_p,
            "best_log_p": -np.log10(best_p) if best_p > 0 else 100,
            "sig_count": sum(1 for p in p_values if p < 5e-8)
        }
    
    def _extract_genes(self, group_df: pd.DataFrame) -> set:
        genes = set()
        for _, row in group_df.iterrows():
            for col in self.GENE_COLUMNS:
                val = self._safe_get(row, col)
                if val:
                    for gene in str(val).split(','):
                        gene = gene.strip().upper()
                        if gene and gene not in self.EXCLUDED_GENES:
                            genes.add(gene)
        return genes
    
    def _extract_traits(self, group_df: pd.DataFrame) -> set:
        return {str(self._safe_get(row, 'DISEASE/TRAIT')).strip() 
                for _, row in group_df.iterrows() 
                if self._safe_get(row, 'DISEASE/TRAIT')}
    
    def _extract_risk_alleles(self, group_df: pd.DataFrame) -> Tuple[List[Dict], float]:
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
                    
                    # Extract detailed variant information
                    variant_info = self._extract_variant_details(row)
                    
                    risk_alleles.append({
                        "risk_allele": risk_allele,
                        "snp": self._safe_get(row, 'SNPS'),
                        "p_value": self._safe_get(row, 'P-VALUE'),
                        "or_beta": or_beta,
                        "risk_score": risk_score,
                        "disease_trait": self._safe_get(row, 'DISEASE/TRAIT'),
                        "mapped_gene": self._safe_get(row, 'MAPPED_GENE'),
                        **variant_info  # Include all variant details
                    })
                except (ValueError, TypeError):
                    continue
        
        return risk_alleles, max_risk
    
    def _extract_p_values(self, group_df: pd.DataFrame) -> List[float]:
        p_values = []
        for _, row in group_df.iterrows():
            p_val = self._safe_get(row, 'P-VALUE')
            if p_val:
                try:
                    p_float = (float(p_val) if isinstance(p_val, str) and ('E-' in p_val or 'e-' in p_val) 
                              else float(p_val))
                    p_values.append(p_float)
                except (ValueError, TypeError):
                    continue
        return p_values
    
    def _extract_variant_details(self, row) -> Dict:
        """Extract comprehensive variant information from a row."""
        context = self._safe_get(row, 'CONTEXT') or ''
        chr_id = self._safe_get(row, 'CHR_ID')
        chr_pos = self._safe_get(row, 'CHR_POS')
        
        # Parse chromosome position
        try:
            chr_pos_int = int(chr_pos) if chr_pos else None
        except (ValueError, TypeError):
            chr_pos_int = None
        
        # Categorize variant context
        variant_category = self._categorize_variant_context(context)
        
        # Extract distance information
        upstream_dist = self._safe_get(row, 'UPSTREAM_GENE_DISTANCE')
        downstream_dist = self._safe_get(row, 'DOWNSTREAM_GENE_DISTANCE')
        
        try:
            upstream_dist = int(upstream_dist) if upstream_dist else None
        except (ValueError, TypeError):
            upstream_dist = None
            
        try:
            downstream_dist = int(downstream_dist) if downstream_dist else None
        except (ValueError, TypeError):
            downstream_dist = None
        
        return {
            "chromosome": chr_id,
            "position": chr_pos_int,
            "region": self._safe_get(row, 'REGION'),
            "context": context,
            "variant_category": variant_category,
            "is_intergenic": bool(self._safe_get(row, 'INTERGENIC')),
            "risk_allele_frequency": self._safe_get(row, 'RISK ALLELE FREQUENCY'),
            "upstream_gene_id": self._safe_get(row, 'UPSTREAM_GENE_ID'),
            "downstream_gene_id": self._safe_get(row, 'DOWNSTREAM_GENE_ID'),
            "upstream_distance": upstream_dist,
            "downstream_distance": downstream_dist
        }
    
    def _categorize_variant_context(self, context: str) -> str:
        """Categorize variant context into functional groups."""
        if not context:
            return "unknown"
        
        context_lower = context.lower()
        
        if any(ctx in context_lower for ctx in self.CODING_CONTEXTS):
            return "coding"
        elif any(ctx in context_lower for ctx in self.UTR_CONTEXTS):
            return "utr"
        elif any(ctx in context_lower for ctx in self.SPLICE_CONTEXTS):
            return "splice"
        elif any(ctx in context_lower for ctx in self.REGULATORY_CONTEXTS):
            return "regulatory"
        elif "intron" in context_lower:
            return "intronic"
        elif "intergenic" in context_lower:
            return "intergenic"
        else:
            return "other"
    
    def _generate_summary(self, studies: List[Dict], sort_by_risk: bool, fps_param: Optional[int], 
                         is_gene_query: bool) -> Dict:
        gene_scores, trait_scores, snp_scores = {}, {}, {}
        proteins = set()
        variant_annotations = {}  # Dictionary to store raw variant information
        
        for study in studies:
            study_risk = study.get('max_risk', 0)
            study_pval = study.get('best_pvalue', 1.0)
            
            self._collect_scores(study['related_genes'], gene_scores, study_risk, study_pval)
            self._collect_scores(study['disease_traits'], trait_scores, study_risk, study_pval)
            self._extract_proteins_from_traits(study['disease_traits'], proteins)
            self._collect_snp_scores(study['risk_alleles'], snp_scores)
            self._collect_variant_annotations(study['risk_alleles'], variant_annotations)
        
        sorted_genes = self._sort_by_scores(gene_scores, sort_by_risk)
        sorted_traits = self._sort_by_scores(trait_scores, sort_by_risk)
        sorted_snps = self._sort_by_scores(snp_scores, sort_by_risk)
        
        # Apply FPS
        if fps_param:
            if is_gene_query and len(sorted_traits) > fps_param:
                fps_indices = fps_tfidf(sorted_traits, fps_param)
                sorted_traits = [sorted_traits[i] for i in fps_indices]
            elif not is_gene_query and len(sorted_genes) > fps_param:
                fps_indices = fps_tfidf(sorted_genes, fps_param)
                sorted_genes = [sorted_genes[i] for i in fps_indices]
        
        summary = self._format_summary(sorted_genes, sorted_snps, sorted(proteins), sorted_traits, is_gene_query)
        summary["variant_annotations"] = variant_annotations  # Include raw variant data
        return summary
    
    def _collect_scores(self, items: List[str], scores_dict: Dict, study_risk: float, study_pval: float):
        for item in items:
            if item not in scores_dict:
                scores_dict[item] = {'risk_score': 0, 'best_pvalue': 1.0}
            scores_dict[item]['risk_score'] = max(scores_dict[item]['risk_score'], study_risk)
            scores_dict[item]['best_pvalue'] = min(scores_dict[item]['best_pvalue'], study_pval)
    
    def _extract_proteins_from_traits(self, traits: List[str], proteins: set):
        for trait in traits:
            trait_lower = trait.lower()
            if any(kw in trait_lower for kw in self.PROTEIN_KEYWORDS):
                protein_name = trait
                for replacement in self.PROTEIN_REPLACEMENTS:
                    protein_name = protein_name.replace(replacement, "")
                
                if "ratio" in protein_name:
                    for protein in protein_name.replace("ratio", "").split("/"):
                        if protein.strip():
                            proteins.add(protein.strip())
                elif protein_name.strip():
                    proteins.add(protein_name.strip())
    
    def _collect_snp_scores(self, risk_alleles: List[Dict], snp_scores: Dict):
        for allele in risk_alleles:
            snp = allele.get('snp')
            if snp:
                risk_score = allele.get('risk_score', 0)
                p_val = self._parse_pvalue(allele.get('p_value', 1.0))
                
                if snp not in snp_scores:
                    snp_scores[snp] = {'risk_score': 0, 'best_pvalue': 1.0}
                snp_scores[snp]['risk_score'] = max(snp_scores[snp]['risk_score'], risk_score)
                snp_scores[snp]['best_pvalue'] = min(snp_scores[snp]['best_pvalue'], p_val)
    
    def _sort_by_scores(self, scores_dict: Dict, sort_by_risk: bool) -> List[str]:
        sort_key = (lambda x: scores_dict[x]['risk_score']) if sort_by_risk else (lambda x: scores_dict[x]['best_pvalue'])
        return sorted(scores_dict.keys(), key=sort_key, reverse=sort_by_risk)
    
    def _collect_variant_annotations(self, risk_alleles: List[Dict], variant_annotations: Dict):
        """Collect raw variant information keyed by variant ID for downstream processing."""
        for variant in risk_alleles:
            variant_id = variant.get('snp')
            if not variant_id:
                continue
                
            disease_trait = variant.get('disease_trait')
            if not disease_trait:
                continue
                
            # Create comprehensive variant record with new nested format
            variant_record = {
                'mapped_gene': variant.get('mapped_gene'),
                'context': variant.get('context'),
                'variant_category': variant.get('variant_category'),
                'associated_disease_trait': {
                    disease_trait: {
                        'p_value': variant.get('p_value'),
                        'risk_score': variant.get('risk_score')
                    }
                }
            }
            
            # If variant already exists, merge all disease traits
            if variant_id in variant_annotations:
                new_trait_data = variant_record['associated_disease_trait']
                variant_annotations[variant_id]['associated_disease_trait'].update(new_trait_data)
            else:
                variant_annotations[variant_id] = variant_record
    
    def _format_summary(self, genes: List[str], snps: List[str], proteins: List[str], 
                       traits: List[str], is_gene_query: bool) -> Dict:
        if is_gene_query:
            return {
                "related_genes": genes,
                "high_risk_snps": snps,
                "affected_protein_levels": proteins,
                "associated_disease_traits": traits
            }
        return {
            "related_genes": genes,
            "high_risk_snps": snps,
            "affected_protein_levels": proteins,
            "associated_disease_traits": traits
        }
    
    def _format_gene_results(self, results: Dict) -> Dict:
        if not results.get("found", False):
            return results
        
        filtered_results = {k: results[k] for k in [
            "gene", "found", "total_associations", "total_significant_associations",
            "total_studies_analyzed", "p_value_threshold", "support_types", "support_values",
            "summary_by_high_risk_alleles", "summary_by_significance",
            "studies_by_high_risk_alleles", "studies_by_significance"
        ] if k in results}
        
        filtered_results["gwas_linked_genes"] = set([
            *results["summary_by_high_risk_alleles"]["related_genes"],
            *results["summary_by_significance"]["related_genes"]
        ])
        return filtered_results
    
    def _format_trait_results(self, results: Dict) -> Dict:
        if not results.get("found", False):
            return results
        
        filtered_results = {k: results[k] for k in [
            "trait", "found", "total_associations", "total_significant_associations",
            "total_studies_analyzed", "p_value_threshold", 
            "summary_by_high_risk_alleles", "summary_by_significance",
            "studies_by_high_risk_alleles", "studies_by_significance"
        ] if k in results}
        
        filtered_results["gwas_linked_genes"] = set([
            *results["summary_by_high_risk_alleles"]["related_genes"],
            *results["summary_by_significance"]["related_genes"]
        ])
        return filtered_results
    
    def _parse_pvalue(self, p_val) -> float:
        try:
            if isinstance(p_val, str) and ('E-' in p_val or 'e-' in p_val):
                return float(p_val)
            return float(p_val) if p_val else 1.0
        except (ValueError, TypeError):
            return 1.0
    
    def _safe_get(self, row, column: str):
        try:
            value = row[column]
            return None if pd.isna(value) or value == '' or str(value).lower() == 'nan' else value
        except (KeyError, IndexError):
            return None


if __name__ == "__main__":
    # Test gene search with raw variants
    gene = "APOE"
    results = query_gene_associations(gene, fps_disease_traits=3)
    print(f"Gene {gene}: Found={results['found']}, Associations={results['total_associations']}")
    
    if results['found']:
        raw_variants = results['summary_by_high_risk_alleles']['variant_annotations']
        print(f"  Raw variants collected: {len(raw_variants)}")
        
        # Show sample raw variant records
        sample_variants = list(raw_variants.items())[:3]
        for variant_id, variant_data in sample_variants:
            print(f"  Variant {variant_id}:")
            print(f"    Location: chr{variant_data['chromosome']}:{variant_data['position']}")
            print(f"    Context: {variant_data['context']} ({variant_data['variant_category']})")
            print(f"    Risk allele: {variant_data['risk_allele']}")
            trait_count = len(variant_data['associated_disease_trait'])
            print(f"    Associated disease traits ({trait_count} total):")
            for trait, trait_data in variant_data['associated_disease_trait'].items():
                print(f"      {trait}: p_value={trait_data['p_value']}, risk_score={trait_data['risk_score']}")
    
    print("\n" + "="*50)
    
    # Test trait search with raw variants
    trait = "Type 2 diabetes"
    trait_results = query_trait_associations(trait, fps_genes=3, exact_match=True)
    print(f"Trait '{trait}': Found={trait_results['found']}, Associations={trait_results['total_associations']}")
    
    if trait_results['found']:
        raw_variants = trait_results['summary_by_high_risk_alleles']['variant_annotations']
        print(f"  Raw variants collected: {len(raw_variants)}")
        
        # Show variant categories distribution
        categories = {}
        for variant_data in raw_variants.values():
            cat = variant_data['variant_category']
            categories[cat] = categories.get(cat, 0) + 1
        print(f"  Category distribution: {categories}")
        
        # Show sample high-risk variants - need to get max risk score from disease traits
        def get_max_risk_score(variant_data):
            max_risk = 0
            for trait_data in variant_data.get('associated_disease_trait', {}).values():
                max_risk = max(max_risk, trait_data.get('risk_score', 0) or 0)
            return max_risk
        
        high_risk_variants = sorted(raw_variants.items(), 
                                  key=lambda x: get_max_risk_score(x[1]), reverse=True)[:3]
        print(f"  Top risk variants:")
        for variant_id, variant_data in high_risk_variants:
            trait_count = len(variant_data['associated_disease_trait'])
            print(f"    {variant_id}: chr{variant_data['chromosome']}:{variant_data['position']}")
            print(f"      {variant_data['context']}, Max Risk: {get_max_risk_score(variant_data):.2f}")
            print(f"      Associated with {trait_count} disease traits:")
            for trait, trait_data in variant_data['associated_disease_trait'].items():
                print(f"        {trait}: risk_score={trait_data['risk_score']:.2f}, p_value={trait_data['p_value']}")
    
    print(f"\n  Raw variant structure (keys): {list(raw_variants[list(raw_variants.keys())[0]].keys()) if raw_variants else 'No variants'}")