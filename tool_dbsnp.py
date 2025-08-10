from __future__ import annotations
import time
import re
from typing import Set

from config import DEBUG, DBSNP_DEFAULT_ASSEMBLY
from tools.dbsnp.query import get_variant_info
from state import State


def dbsnp_variants_agent(state: "State", assembly: str = None) -> "State":
    """
    LangGraph node that annotates SNPs with dbSNP variant information.
    Extracts rsIDs from GWAS associations and queries dbSNP for each variant.

    Args:
        state: Current workflow state containing gene symbols and GWAS associations
        assembly: Genome assembly to use (default: config.DBSNP_DEFAULT_ASSEMBLY)
                 Options: "GRCh38", "GRCh37", "all"

    Returns:
        Updated state with dbSNP variant information
    """
    if assembly is None:
        assembly = DBSNP_DEFAULT_ASSEMBLY
    
    if DEBUG:
        print(f"[dbSNP] Using genome assembly: {assembly}")

    variants = state.get("dbsnp_variants", {}).copy()
    
    # Extract rsIDs from GWAS associations
    rsids = _extract_rsids_from_gwas(state)
    
    if DEBUG:
        print(f"[dbSNP] Found {len(rsids)} unique rsIDs to query")
    
    for gene, gene_rsids in rsids.items():
        if gene in variants:
            if DEBUG:
                print(f"[dbSNP] Skipping {gene} - already processed")
            continue

        variants[gene] = {}
        
        for i, rsid in enumerate(gene_rsids):

            if rsid in variants[gene]:
                if DEBUG:
                    print(f"[dbSNP] Skipping {rsid} - already processed")
                continue

            if DEBUG:
                print(f"[dbSNP] Querying variant information for: {rsid}")

            try:
                # Query dbSNP variant information
                result = get_variant_info(rsid, assembly)
                
                variants[gene][rsid] = {
                    "rsid": result.get("rsid", rsid),
                    "found": True,
                    "coordinates": result.get("coordinates", []),
                    "coordinate_count": len(result.get("coordinates", [])),
                    "chromosomes": list(set(coord.get("chrom", "Unknown") for coord in result.get("coordinates", []))),
                    "assembly_filter": result.get("assembly_filter", assembly)
                }
                
                if DEBUG:
                    coord_count = len(result.get("coordinates", []))
                    print(f"[dbSNP] {rsid}: found, {coord_count} coordinates")
                    print(variants[gene][rsid])

                    
            except Exception as exc:
                if DEBUG:
                    print(f"[dbSNP] Error querying {rsid}: {exc}")
                variants[gene][rsid] = {
                    "rsid": rsid,
                    "found": False,
                    "coordinates": [],
                    "coordinate_count": 0,
                    "chromosomes": [],
                    "error": str(exc)
                }

            # Courteous pause
            time.sleep(0.1)

    # import ipdb; ipdb.set_trace()
    return {**state, "dbsnp_variants": variants}


def _extract_rsids_from_gwas(state: "State") -> Set[str]:
    """
    Extract unique rsIDs from GWAS associations in the state.
    
    Parameters
    ----------
    state
        Current graph state containing GWAS associations.
        
    Returns
    -------
    Set[str]
        Set of unique rsIDs found in GWAS data.
    """
    rsids = dict()
    associations = state.get("gwas_associations", {})
    
    for gene, data in associations.items():
        gene_rsids = set()
        if not data.get("found", False):
            continue
            
        # Extract from high risk alleles summary
        high_risk_summary = data.get("summary_by_high_risk_alleles", {})
        high_risk_snps = high_risk_summary.get("gwas_high_risk_snps", [])
        for snp in high_risk_snps:
            if isinstance(snp, str):
                extracted = _extract_rsids_from_text(snp)
                gene_rsids.update(extracted)
        
        # Extract from significance summary
        significance_summary = data.get("summary_by_significance", {})
        sig_snps = significance_summary.get("gwas_high_risk_snps", [])
        for snp in sig_snps:
            if isinstance(snp, str):
                extracted = _extract_rsids_from_text(snp)
                gene_rsids.update(extracted)

        rsids[gene] = gene_rsids
    
    return rsids


def _extract_rsids_from_text(text: str) -> Set[str]:
    """
    Extract rsIDs from text using regex pattern matching.
    
    Parameters
    ----------
    text
        Text that may contain rsIDs.
        
    Returns
    -------
    Set[str]
        Set of rsIDs found in the text.
    """
    # Pattern to match rsIDs (rs followed by digits)
    pattern = r'\brs\d+\b'
    matches = re.findall(pattern, text, re.IGNORECASE)
    return set(match.lower() for match in matches)


def has_dbsnp_variants(state: "State") -> bool:
    """
    Edge-condition helper for LangGraph: returns `True` if any dbSNP variants were found.
    """
    variants = state.get("dbsnp_variants", {})
    has_variants = any(
        variant.get("found", False) and variant.get("coordinate_count", 0) > 0
        for gene_variants in variants.values()
        for variant in gene_variants.values()
    )
    
    if DEBUG:
        print(f"[has_dbsnp_variants] variants found: {has_variants}")
    
    return has_variants


def has_dbsnp_coordinates(state: "State") -> bool:
    """
    Edge-condition helper for LangGraph: returns `True` if any variants have coordinate information.
    """
    variants = state.get("dbsnp_variants", {})
    has_coords = any(
        len(variant.get("coordinates", [])) > 0
        for gene_variants in variants.values()
        for variant in gene_variants.values()
    )
    
    if DEBUG:
        print(f"[has_dbsnp_coordinates] coordinates found: {has_coords}")
    
    return has_coords 