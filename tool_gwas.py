from __future__ import annotations
import time
from typing import List, Dict, Any

from config import DEBUG
from tools.gwas.query import query_gene_associations, query_trait_associations
from state import State


def _create_empty_association_record(identifier: str, is_gene: bool = True) -> Dict[str, Any]:
    """Create a standardized empty association record."""
    base_record = {
        "found": False,
        "total_associations": 0,
        "total_significant_associations": 0,
        "total_studies_analyzed": 0,
        "summary_by_high_risk_alleles": {
            "related_genes": [],
            "high_risk_snps": [],
            "affected_protein_levels": [],
            "associated_disease_traits": []
        },
        "summary_by_significance": {
            "related_genes": [],
            "high_risk_snps": [],
            "affected_protein_levels": [],
            "associated_disease_traits": []
        },
        "variant_annotations": {}
    }
    
    if is_gene:
        base_record["gene"] = identifier
    else:
        base_record["trait_term"] = identifier
    
    return base_record


def _extract_genes_from_trait_result(result: Dict[str, Any]) -> List[str]:
    """
    Extract top 20 genes from trait query results based on multiple criteria:
    - Top 20 by number of studies
    - Top 20 by number of significant hits  
    - Top 20 by maximum risk scores
    """
    if not result.get("found", False):
        return []
    
    # Collect gene statistics from studies
    gene_stats = {}
    
    # Process high-risk studies
    for study in result.get("studies_by_high_risk_alleles", []):
        max_risk = study.get("max_risk", 0)
        sig_count = study.get("sig_count", 0)
        for gene in study.get("related_genes", []):
            if gene not in gene_stats:
                gene_stats[gene] = {"studies": 0, "total_sig_hits": 0, "max_risk": 0}
            gene_stats[gene]["studies"] += 1
            gene_stats[gene]["total_sig_hits"] += sig_count
            gene_stats[gene]["max_risk"] = max(gene_stats[gene]["max_risk"], max_risk)
    
    # Process significance studies
    for study in result.get("studies_by_significance", []):
        max_risk = study.get("max_risk", 0)
        sig_count = study.get("sig_count", 0)
        for gene in study.get("related_genes", []):
            if gene not in gene_stats:
                gene_stats[gene] = {"studies": 0, "total_sig_hits": 0, "max_risk": 0}
            gene_stats[gene]["studies"] += 1
            gene_stats[gene]["total_sig_hits"] += sig_count
            gene_stats[gene]["max_risk"] = max(gene_stats[gene]["max_risk"], max_risk)
    
    if not gene_stats:
        return []
    
    # Get top 20 genes by each criterion
    top_by_studies = sorted(gene_stats.items(), key=lambda x: x[1]["studies"], reverse=True)[:20]
    top_by_sig_hits = sorted(gene_stats.items(), key=lambda x: x[1]["total_sig_hits"], reverse=True)[:20]
    top_by_max_risk = sorted(gene_stats.items(), key=lambda x: x[1]["max_risk"], reverse=True)[:20]
    
    # Combine all top genes and deduplicate while preserving priority order
    selected_genes = []
    for gene, _ in top_by_studies + top_by_sig_hits + top_by_max_risk:
        if gene not in selected_genes:
            selected_genes.append(gene)
    
    return selected_genes


def _log_query_result(query_type: str, identifier: str, result: Dict[str, Any], 
                     related_genes: List[str] = None) -> None:
    """Centralized logging for query results."""
    if not DEBUG:
        return
    
    if result.get("found", False):
        total_assoc = result.get("total_associations", 0)
        significant_assoc = result.get("total_significant_associations", 0)
        studies = result.get("total_studies_analyzed", 0)
        
        if query_type == "gene":
            variant_count = len(result.get("summary_by_high_risk_alleles", {}).get("variant_annotations", {}))
            print(f"[GWAS] {identifier}: found, {total_assoc} total associations, "
                  f"{significant_assoc} significant, {studies} studies, {variant_count} variant annotations")
        else:  # trait
            genes_found = len(related_genes) if related_genes else 0
            print(f"[query_by_trait_agent] Found: {total_assoc} total associations, "
                  f"{significant_assoc} significant, {studies} studies, {genes_found} selected genes")
            if related_genes:
                print(f"[query_by_trait_agent] Top genes selected by studies/significance/risk: {related_genes[:10]}...")
    else:
        print(f"[{query_type.upper()}] {identifier}: not found")


def gwas_associations_agent(state: "State", mode: str = "summary") -> "State":
    """
    Query GWAS associations for each gene in the state.
    
    Annotates genes with comprehensive GWAS association data including:
    - Association counts and significance statistics
    - High-risk allele summaries
    - Variant annotations
    
    Parameters
    ----------
    state : State
        Current graph state containing genes to query
    mode : str
        Mode of operation: "summary" (default) or "extensive"
        
    Returns
    -------
    State
        Updated state with "gwas_associations" field
    """
    associations = state.get("gwas_associations", {}).copy()
    variants = state.get("dbsnp_variants", {}).copy()
    gene_list = list(set(state.get("genes", [])))
    
    for gene in gene_list:
        if gene in associations:
            if DEBUG:
                print(f"[GWAS] Skipping {gene} - already processed")
            continue

        if DEBUG:
            print(f"[GWAS] Querying associations for gene: {gene}")

        try:
            if mode == "extensive":
                print(f"[GWAS] Querying associations for gene: {gene} in extensive mode")
                result = query_gene_associations(
                    gene_symbol=gene,
                        db_path="local_dbs",
                        fps_disease_traits=200,
                        top_studies_by_risk=100,
                        top_studies_by_significance=100,
                    )
            else:
                print(f"[GWAS] Querying associations for gene: {gene} in summary mode")
                result = query_gene_associations(
                    gene_symbol=gene,
                    db_path="local_dbs",
                    fps_disease_traits=20,
                )
            
            if result.get("found", True):
                associations[gene] = {
                    "gene": gene,
                    "found": True,
                    "total_associations": result.get("total_associations", 0),
                    "total_significant_associations": result.get("total_significant_associations", 0),
                    "total_studies_analyzed": result.get("total_studies_analyzed", 0),
                    "summary_by_high_risk_alleles": result.get("summary_by_high_risk_alleles", {}),
                    "summary_by_significance": result.get("summary_by_significance", {}),

                }
                # Extract variant_annotations from both summary sections
                risk_variants = result.get("summary_by_high_risk_alleles", {}).get("variant_annotations", {})
                sig_variants = result.get("summary_by_significance", {}).get("variant_annotations", {})
                
                # Combine variant annotations from both sections
                all_variants = {}
                all_variants.update(risk_variants)
                all_variants.update(sig_variants)
                variants[gene] = all_variants
                
                # Pop out the variant annotations from summary_by_high_risk_alleles, and from summary_by_significance
                associations[gene]["summary_by_high_risk_alleles"].pop("variant_annotations", None)
                associations[gene]["summary_by_significance"].pop("variant_annotations", None)
                
                _log_query_result("gene", gene, result)
            else:
                associations[gene] = _create_empty_association_record(gene, is_gene=True)
                _log_query_result("gene", gene, result)
                
        except Exception as exc:
            if DEBUG:
                print(f"[GWAS] Error querying {gene}: {exc}")
            
            error_record = _create_empty_association_record(gene, is_gene=True)
            error_record["error"] = str(exc)
            associations[gene] = error_record

        # Rate limiting
        time.sleep(0.1)

    return {"gwas_associations": associations, "dbsnp_variants": variants}


def query_by_trait_agent(state: "State") -> "State":
    """
    Query GWAS associations by disease/trait terms.
    
    Uses extracted traits from state or falls back to full user input.
    Populates genes field with discovered genes for downstream processing.
    
    Parameters
    ----------
    state : State
        Current graph state
        
    Returns
    -------
    State
        Updated state with "trait_associations" and "genes" fields
    """
    user_input = state["messages"][-1]["content"]
    
    # Determine trait term to query
    extracted_traits = state.get("traits", [])
    if extracted_traits:
        trait_term = extracted_traits[0]  # Use most prominent trait
        if DEBUG:
            print(f"[query_by_trait_agent] Using extracted trait: {trait_term} "
                  f"(from {len(extracted_traits)} total traits)")
    else:
        trait_term = user_input.strip()
        if DEBUG:
            print(f"[query_by_trait_agent] No extracted traits found, using full query: {trait_term}")
    
    if DEBUG:
        print(f"[query_by_trait_agent] Querying trait associations for: {trait_term}")

    try:
        result = query_trait_associations(
            trait_term=trait_term,
            db_path="local_dbs",
            fps_genes=60,  # Here I think we should be generous?
            exact_match=False
        )
        
        trait_associations = {
            "trait_term": trait_term,
            "found": result.get("found", False),
            "total_associations": result.get("total_associations", 0),
            "total_significant_associations": result.get("total_significant_associations", 0),
            "total_studies_analyzed": result.get("total_studies_analyzed", 0),
            "summary_by_high_risk_alleles": result.get("summary_by_high_risk_alleles", {}),
            "summary_by_significance": result.get("summary_by_significance", {}),
        }
        
        # Extract variant_annotations from both summary sections
        risk_variants = result.get("summary_by_high_risk_alleles", {}).get("variant_annotations", {})
        sig_variants = result.get("summary_by_significance", {}).get("variant_annotations", {})
        
        # Combine variant annotations from both sections
        all_variants = {}
        all_variants.update(risk_variants)
        all_variants.update(sig_variants)
        trait_associations["variant_annotations"] = all_variants
        
        # Extract genes for downstream processing
        related_genes = _extract_genes_from_trait_result(result)
        _log_query_result("trait", trait_term, result, related_genes)
                
    except Exception as exc:
        if DEBUG:
            print(f"[query_by_trait_agent] Error querying trait '{trait_term}': {exc}")
        
        trait_associations = _create_empty_association_record(trait_term, is_gene=False)
        trait_associations["error"] = str(exc)
        related_genes = []

    return {"trait_associations": trait_associations, "genes": related_genes}


def has_gwas_associations(state: "State") -> bool:
    """Check if any GWAS associations were found for genes."""
    associations = state.get("gwas_associations", {})
    has_associations = any(
        assoc.get("found", False) and assoc.get("total_associations", 0) > 0
        for assoc in associations.values()
    )
    
    if DEBUG:
        print(f"[has_gwas_associations] associations found: {has_associations}")
    
    return has_associations


def has_significant_gwas_associations(state: "State") -> bool:
    """Check if any significant GWAS associations were found for genes."""
    associations = state.get("gwas_associations", {})
    has_significant = any(
        assoc.get("total_significant_associations", 0) > 0
        for assoc in associations.values()
    )
    
    if DEBUG:
        print(f"[has_significant_gwas_associations] significant associations found: {has_significant}")
    
    return has_significant


def has_trait_associations(state: "State") -> bool:
    """Check if any trait associations were found."""
    trait_associations = state.get("trait_associations", {})
    has_associations = (
        trait_associations.get("found", False) and 
        trait_associations.get("total_associations", 0) > 0
    )
    
    if DEBUG:
        print(f"[has_trait_associations] trait associations found: {has_associations}")
    
    return has_associations