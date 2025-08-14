from __future__ import annotations
import time
import re
from typing import Set

from config import DEBUG, DBSNP_DEFAULT_ASSEMBLY
from tools.dbsnp.query import get_variant_info
from state import State
from tools.gencode.query import annotate_variant


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
                print(f"[dbSNP] Querying variant informations for: {rsid}")

            try:
                # Query dbSNP variant information
                result = get_variant_info(rsid, assembly)
                
                allele_frequencies = result.get("allele_frequencies", [])
                coordinates = result.get("coordinates", [])
                
                variants[gene][rsid] = {
                    "rsid": result.get("rsid", rsid),
                    "found": True,
                    "coordinates": coordinates,
                    "coordinate_count": len(coordinates),
                    "chromosomes": list(set(coord.get("chrom", "Unknown") for coord in coordinates)),
                    "allele_frequencies": allele_frequencies,
                    "frequency_study_count": len(allele_frequencies),
                    "frequency_studies": [freq.get("study_name", "Unknown") for freq in allele_frequencies],
                    "assembly_filter": result.get("assembly_filter", assembly)
                }
                
                if DEBUG:
                    print(f"[dbSNP] {rsid}: {len(coordinates)} coords, {len(allele_frequencies)} studies")
                    
                # go through coordinates, find build 38 chrom and position and annotate
                # set deafult annots to empty list
                variants[gene][rsid].setdefault("annotations", [])
                for coord in result.get("coordinates", []):
                    build = coord.get("assembly", "Unknown")
                    if 'GRCh37' in build or 'hg19' in build:
                        continue
                    chrom = coord.get("chrom", "Unknown")
                    pos = coord.get("pos", 0)
                    annots = annotate_variant(chrom, pos)
                    if annots:
                        variants[gene][rsid]["annotations"].extend(annots)
                        if DEBUG:
                            print(f"[dbSNP] Annotated {rsid} at {chrom}:{pos} with {annots}")
                        break
                    
            except Exception as exc:
                if DEBUG:
                    print(f"[dbSNP] Error querying {rsid}: {exc}")
                variants[gene][rsid] = {
                    "rsid": rsid,
                    "found": False,
                    "coordinates": [],
                    "coordinate_count": 0,
                    "chromosomes": [],
                    "allele_frequencies": [],
                    "frequency_study_count": 0,
                    "frequency_studies": [],
                    "error": str(exc)
                }

            # Courteous pause
            time.sleep(0.1)

    # Generate summary statistics for each gene
    variant_summaries = {}
    for gene, gene_variants in variants.items():
        if not gene_variants:
            continue
            
        # Calculate summary statistics
        found_variants = [v for v in gene_variants.values() if v.get("found", False)]
        
        if found_variants:
            chromosomes = set()
            for v in found_variants:
                chromosomes.update(v.get("chromosomes", []))
            
            # Categorize by frequency: rare (<1%), common (1-15%), very common (>15%)
            rare_variants, common_variants, very_common_variants = [], [], []
            
            for v in found_variants:
                afs = [freq.get('allele_frequency', 0) for freq in v.get("allele_frequencies", [])]
                min_af, _ = (min(afs), max(afs)) if afs else (0, 0)
                
                if min_af < 0.01:
                    rare_variants.append(v.get("rsid"))
                elif min_af < 0.15:
                    common_variants.append(v.get("rsid"))
                else:
                    very_common_variants.append(v.get("rsid"))
            
            variant_summaries[gene] = {
                "total_variants": len(found_variants),
                "chromosomes": sorted(list(chromosomes)),
                "assembly": assembly,
                "rare_variants": len(rare_variants),
                "common_variants": len(common_variants),
                "very_common_variants": len(very_common_variants),
                "rare_rsids (<1%)": rare_variants,
                "common_rsids (1-15%)": common_variants,
                "very_common_rsids (>15%)": very_common_variants,
            }

    return {"dbsnp_variants": variants, "dbsnp_summaries": variant_summaries}


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