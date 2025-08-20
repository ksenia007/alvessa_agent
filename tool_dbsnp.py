from __future__ import annotations
import time
import re
from typing import Set, Dict, Any, List

from config import DEBUG, DBSNP_DEFAULT_ASSEMBLY
from tools.dbsnp.query import get_variant_info
from state import State
from tools.gencode.query import annotate_variant


def _create_coordinate_summary(rsid: str, coordinates: List[Dict[str, Any]]) -> str:
    """Generates a concise, human-readable summary of variant coordinates."""
    if not coordinates:
        return f"{rsid} has no coordinate data available"

    primary_coord = next((c for c in coordinates if 'GRCh38' in c.get("assembly", "")), coordinates[0])
    chrom = primary_coord.get("chrom", "Unknown")
    pos38 = primary_coord.get("pos", "Unknown")
    
    return f"{rsid} is located at chr:{chrom}, pos: {pos38} (GRCh38)"

def _create_frequency_summary(allele_frequencies: List[Dict[str, Any]]) -> str:
    """
    Generates a concise, human-readable summary of allele frequency data,
    highlighting the top 2 and bottom 2 representative studies.
    """
    if not allele_frequencies:
        return "No allele frequency data is available."

    # Filter for valid entries that have the required keys
    valid_freqs = [f for f in allele_frequencies if 'allele_frequency' in f and 'study_name' in f]
    
    if not valid_freqs:
        num_studies = len(set(f.get("study_name") for f in allele_frequencies))
        return f"Frequency data from {num_studies} studies was found, but no valid frequencies were calculated."
        
    # Sort all valid observations by allele frequency
    sorted_freqs = sorted(valid_freqs, key=lambda x: x['allele_frequency'])
    count = len(sorted_freqs)
    num_unique_studies = len(set(f.get("study_name") for f in sorted_freqs))
    
    # Helper function to format a single observation string
    def format_obs(entry):
        return f"{entry['allele_frequency']:.4f} in {entry.get('study_name', 'N/A')}"

    # Handle cases with few observations
    if count == 1:
        return f"The only available frequency observation is {format_obs(sorted_freqs[0])}."
    if count == 2:
        return f"The two available frequency observations are {format_obs(sorted_freqs[0])} and {format_obs(sorted_freqs[1])}."
    if count == 3:
        mean_af = sum(f['allele_frequency'] for f in sorted_freqs) / count
        return (f"Frequencies from 3 observations range from {format_obs(sorted_freqs[0])} to "
                f"{format_obs(sorted_freqs[2])}, with a mean of {mean_af:.4f}.")

    # Handle the main case for 4 or more observations
    bottom1 = format_obs(sorted_freqs[0])
    bottom2 = format_obs(sorted_freqs[1])
    top1 = format_obs(sorted_freqs[-2])
    top2 = format_obs(sorted_freqs[-1])
    
    return (f"Allele frequencies from {count} observations across {num_unique_studies} studies range from as low as "
            f"{bottom1} and {bottom2}, to as high as {top1} and {top2}.")



def _create_af_summary_sorted_by_counts(allele_frequencies: List[Dict[str, Any]]) -> str:
    """
    Generates a summary of allele frequencies from studies that are first sorted
    by their sample size (total_count). This highlights the frequencies observed
    in the smallest and largest available cohorts.
    """
    if not allele_frequencies:
        return "No allele frequency data is available."

    # Filter for valid entries that have all required keys for sorting and reporting
    valid_data = [
        f for f in allele_frequencies
        if ('total_count' in f) and ('study_name' in f) and ('allele_frequency' in f)
    ]
    
    if not valid_data:
        return "No studies with complete allele frequency and total count data were found."
        
    # Sort all valid observations by total_count (sample size)
    sorted_by_counts = sorted(valid_data, key=lambda x: x['total_count'])
    count = len(sorted_by_counts)
    
    # Helper function to format the output string
    # It reports the allele frequency but includes the sample size for context
    def format_obs(entry):
        af = entry['allele_frequency']
        study = entry.get('study_name', 'N/A')
        n = entry['total_count']
        return f"{af:.4f} in {study} (n={n:,})"

    # Handle cases with few observations
    if count == 1:
        return f"The only available study reports a frequency of {format_obs(sorted_by_counts[0])}."
    if count == 2:
        # sorted_by_counts[0] is the smaller cohort, [1] is the larger
        return (f"The allele frequency is {format_obs(sorted_by_counts[0])} in the smaller cohort and "
                f"{format_obs(sorted_by_counts[1])} in the larger one.")
    if count == 3:
        return (f"From three cohorts sorted by size, the reported allele frequencies are "
                f"{format_obs(sorted_by_counts[0])} (smallest), "
                f"{format_obs(sorted_by_counts[1])} (medium), and "
                f"{format_obs(sorted_by_counts[2])} (largest).")

    # Handle the main case for 4 or more observations
    # Get the studies with the two smallest sample sizes
    bottom1 = format_obs(sorted_by_counts[0])
    bottom2 = format_obs(sorted_by_counts[1])
    # Get the studies with the two largest sample sizes
    top1 = format_obs(sorted_by_counts[-2])
    top2 = format_obs(sorted_by_counts[-1])
    
    num_unique_studies = len(set(f.get("study_name") for f in sorted_by_counts))
    
    return (f"From {count} observations across {num_unique_studies} studies sorted by sample size, "
            f"the reported allele frequencies range from {bottom1} and {bottom2} in the smallest cohorts, "
            f"to {top1} and {top2} in the largest cohorts.")
    

def dbsnp_variants_agent(state: "State", assembly: str = None) -> "State":
    """
    LangGraph node that annotates SNPs with dbSNP variant information.
    Extracts rsIDs from GWAS associations and queries dbSNP for each variant.
    """
    assembly = assembly or DBSNP_DEFAULT_ASSEMBLY
    if DEBUG:
        print(f"[dbSNP] Using genome assembly: {assembly}")

    variants = state.get("dbsnp_variants", {}).copy()
    
    rsids = _extract_rsids_from_gwas(state)
    
    if DEBUG:
        print(f"[dbSNP] Found {sum(len(v) for v in rsids.values())} total rsIDs to query across genes.")
    
    for gene, gene_rsids in rsids.items():
        if gene in variants:
            if DEBUG:
                print(f"[dbSNP] Skipping {gene} - already processed")
            continue

        variants[gene] = {}
        
        for rsid in gene_rsids:
            if rsid in variants[gene]:
                if DEBUG:
                    print(f"[dbSNP] Skipping {rsid} - already processed")
                continue

            if DEBUG:
                print(f"[dbSNP] Querying variant information for: {rsid}")

            try:
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
                    
                variants[gene][rsid].setdefault("annotations", [])
                for coord in result.get("coordinates", []):
                    build = coord.get("assembly", "Unknown")
                    if 'GRCh37' in build or 'hg19' in build:
                        continue
                    
                    chrom = coord.get("chrom", "Unknown")
                    pos = coord.get("pos", 0)
                    annots = annotate_variant(chrom, pos)
                    if annots:
                        # Correctly append the entire annotation dictionary as one item in the list
                        variants[gene][rsid]["annotations"].append(annots)
                        if DEBUG:
                            print(f"[dbSNP] Annotated {rsid} at {chrom}:{pos} with features.")
                        break
                    
            except Exception as exc:
                if DEBUG:
                    print(f"[dbSNP] Error querying {rsid}: {exc}")
                variants[gene][rsid] = {
                    "rsid": rsid, "found": False, "coordinates": [], "coordinate_count": 0,
                    "chromosomes": [], "allele_frequencies": [], "frequency_study_count": 0,
                    "frequency_studies": [], "error": str(exc)
                }
            time.sleep(0.1)

    # Generate new, improved per-variant summaries from the detailed variant data
    dbsnp_summaries = {}
    for gene, gene_variants in variants.items():
        if not gene_variants:
            continue
        
        dbsnp_summaries[gene] = {}
        for rsid, data in gene_variants.items():
            if not data.get("found"):
                dbsnp_summaries[gene][rsid] = {"summary": f"Variant {rsid} was not found in dbSNP.", "stats": {}}
                continue
            
            coord_summary = _create_coordinate_summary(rsid, data.get("coordinates", []))
            freq_summary = _create_frequency_summary(data.get("allele_frequencies", []))
            freq_summary_sorted_by_counts = _create_af_summary_sorted_by_counts(data.get("allele_frequencies", []))
            
            # Combine into a final summary object for this variant
            dbsnp_summaries[gene][rsid] = {
                "rsid": rsid,
                "summary": f"{coord_summary}. {freq_summary} {freq_summary_sorted_by_counts}"
            }
    # Pop allele frequencies from variants
    for gene, gene_variants in variants.items():
        for rsid, data in gene_variants.items():
            data.pop("allele_frequencies", None)
            data.pop("frequency_study_count", None)
            data.pop("frequency_studies", None)

    return {"dbsnp_variants": variants, "dbsnp_summaries": dbsnp_summaries}


def _extract_rsids_from_gwas(state: "State") -> Dict[str, Set[str]]:
    """Extracts unique rsIDs from GWAS associations in the state, grouped by gene."""
    rsids_by_gene = {}
    associations = state.get("gwas_associations", {})
    
    for gene, data in associations.items():
        if not data.get("found", False):
            continue
        
        gene_rsids = set()
        for summary_key in ["summary_by_high_risk_alleles", "summary_by_significance"]:
            summary = data.get(summary_key, {})
            snps = summary.get("high_risk_snps", [])
            for snp in snps:
                if isinstance(snp, str):
                    gene_rsids.update(_extract_rsids_from_text(snp))
        
        if gene_rsids:
            rsids_by_gene[gene] = gene_rsids
            
    return rsids_by_gene


def _extract_rsids_from_text(text: str) -> Set[str]:
    """Extracts rsIDs from text using regex pattern matching."""
    pattern = r'rs\d+'
    matches = re.findall(pattern, text, re.IGNORECASE)
    return set(match.lower() for match in matches)


def has_dbsnp_variants(state: "State") -> bool:
    """Edge-condition helper: returns `True` if any dbSNP variants were found."""
    variants = state.get("dbsnp_variants", {})
    has_variants = any(
        variant.get("found", False)
        for gene_variants in variants.values()
        for variant in gene_variants.values()
    )
    if DEBUG:
        print(f"[has_dbsnp_variants] Variants found: {has_variants}")
    return has_variants


def has_dbsnp_coordinates(state: "State") -> bool:
    """Edge-condition helper: returns `True` if any variants have coordinate information."""
    variants = state.get("dbsnp_variants", {})
    has_coords = any(
        len(variant.get("coordinates", [])) > 0
        for gene_variants in variants.values()
        for variant in gene_variants.values()
    )
    if DEBUG:
        print(f"[has_dbsnp_coordinates] Coordinates found: {has_coords}")
    return has_coords
