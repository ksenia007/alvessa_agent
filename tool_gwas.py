from __future__ import annotations
import time

from config import DEBUG
from tools.gwas.query import query_gene_associations
from state import State


def gwas_associations_agent(state: "State") -> "State":
    """
    LangGraph node that annotates each gene with GWAS association summaries.

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state with the `"gwas_associations"` field filled with summaries only.
    """
    associations = state.get("gwas_associations", {}).copy()
    gene_list = list(set(state.get("genes", [])))
    for gene in gene_list:
        if gene in associations:
            if DEBUG:
                print(f"[GWAS] Skipping {gene} - already processed")
            continue

        if DEBUG:
            print(f"[GWAS] Querying associations for gene: {gene}")

        try:
            # Query GWAS associations using the simplified function
            result = query_gene_associations(
                gene_symbol=gene,
                db_path="local_dbs",
                fps_disease_traits=20  # Apply FPS to disease traits for manageable output
            )
            
            # Extract only the summaries for downstream processing
            if result.get("found", True):
                associations[gene] = {
                    "gene": gene,
                    "found": True,
                    "total_associations": result.get("total_associations", 0),
                    "total_significant_associations": result.get("total_significant_associations", 0),
                    "total_studies_analyzed": result.get("total_studies_analyzed", 0),
                    "summary_by_high_risk_alleles": result.get("summary_by_high_risk_alleles", {}),
                    "summary_by_significance": result.get("summary_by_significance", {})
                }
                
                if DEBUG:
                    total_assoc = result.get("total_associations", 0)
                    significant_assoc = result.get("total_significant_associations", 0)
                    studies = result.get("total_studies_analyzed", 0)
                    print(f"[GWAS] {gene}: found, {total_assoc} total associations, {significant_assoc} significant, {studies} studies")
            else:
                associations[gene] = {
                    "gene": gene,
                    "found": False,
                    "total_associations": 0,
                    "total_significant_associations": 0,
                    "total_studies_analyzed": 0,
                    "summary_by_high_risk_alleles": {"related_genes": [], "high_risk_snps": [], "proteins": [], "disease_traits": []},
                    "summary_by_significance": {"related_genes": [], "high_risk_snps": [], "proteins": [], "disease_traits": []}
                }
                
                if DEBUG:
                    print(f"[GWAS] {gene}: not found")
                
        except Exception as exc:
            if DEBUG:
                print(f"[GWAS] Error querying {gene}: {exc}")
            associations[gene] = {
                "gene": gene,
                "found": False,
                "total_associations": 0,
                "total_significant_associations": 0,
                "total_studies_analyzed": 0,
                "summary_by_high_risk_alleles": {"related_genes": [], "high_risk_snps": [], "proteins": [], "disease_traits": []},
                "summary_by_significance": {"related_genes": [], "high_risk_snps": [], "proteins": [], "disease_traits": []},
                "error": str(exc)
            }

        # Courteous pause to avoid overwhelming the system
        time.sleep(0.1)

    return {**state, "gwas_associations": associations}


def has_gwas_associations(state: "State") -> bool:
    """
    Edge-condition helper for LangGraph: returns `True` if any GWAS associations were found.
    """
    associations = state.get("gwas_associations", {})
    has_associations = any(
        assoc.get("found", False) and assoc.get("total_associations", 0) > 0
        for assoc in associations.values()
    )
    
    if DEBUG:
        print(f"[has_gwas_associations] associations found: {has_associations}")
    
    return has_associations


def has_significant_gwas_associations(state: "State") -> bool:
    """
    Edge-condition helper for LangGraph: returns `True` if any significant GWAS associations were found.
    """
    associations = state.get("gwas_associations", {})
    has_significant = any(
        assoc.get("total_significant_associations", 0) > 0
        for assoc in associations.values()
    )
    
    if DEBUG:
        print(f"[has_significant_gwas_associations] significant associations found: {has_significant}")
    
    return has_significant 