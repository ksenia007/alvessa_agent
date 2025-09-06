import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import re
from typing import List, Dict

# Import configuration
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from config import DBSNP_DEFAULT_ASSEMBLY, DBSNP_ASSEMBLY_PATTERNS


def get_with_retries(url, max_retries=5, backoff_factor=1.0):
    """Get URL with retry logic for robustness."""
    session = requests.Session()
    retries = Retry(
        total=max_retries,
        backoff_factor=backoff_factor,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"]
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session.get(url, timeout=10)


def fetch_snp_json(rsid: str) -> dict:
    """
    Fetch the dbSNP record for a given rsID via NCBI's Variation API.
    
    Args:
        rsid: The rsID to query (e.g., "rs1333049" or "1333049")
    
    Returns:
        dict: Raw JSON response from dbSNP API
    """
    # Strip "rs" prefix if present
    numeric_id = rsid.lower().lstrip("rs")
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{numeric_id}"
    
    try:
        resp = get_with_retries(url)
        resp.raise_for_status()
    except requests.exceptions.HTTPError as e:
        # If 500 error, retry once
        if resp.status_code == 500:
            resp = requests.get(url, timeout=10)
            resp.raise_for_status()
        else:
            raise
    return resp.json()


def _accession_to_chrom(seq_id: str) -> str:
    """
    Convert chromosome accession to standard format.
    NC_000023.11 → 'X', NC_000019.10 → '19'
    """
    m = re.match(r"NC_0*(\d+)\.", seq_id)
    if not m:
        return seq_id  # Non-chromosomal accession, return as-is
    
    num = int(m.group(1))
    if num == 23:
        return "X"
    if num == 24:
        return "Y"
    return str(num)  # 1-22


def _matches_assembly(assembly_name: str, target_assembly: str) -> bool:
    """
    Check if an assembly name matches the target assembly pattern.
    
    Args:
        assembly_name: Assembly name from dbSNP (e.g., "GRCh38.p14")
        target_assembly: Target assembly pattern (e.g., "GRCh38", "GRCh37", "all")
    
    Returns:
        bool: True if the assembly matches the target pattern
    """
    if target_assembly == "all":
        return True
    
    patterns = DBSNP_ASSEMBLY_PATTERNS.get(target_assembly, [])
    if not patterns:
        # If no patterns defined, do exact match
        return assembly_name.startswith(target_assembly)
    
    return any(pattern in assembly_name for pattern in patterns)


def get_variant_coordinates(rsid: str, assembly: str = None) -> List[Dict[str, str]]:
    """
    Get genomic coordinates for an rsID.
    
    Args:
        rsid: The rsID to query (e.g., "rs1333049" or "1333049")
        assembly: Genome assembly version (default: config.DBSNP_DEFAULT_ASSEMBLY)
                 Options: "GRCh38", "GRCh37", "all"
    
    Returns:
        List of dictionaries with variant coordinates:
        [{"chrom": str, "pos": int, "ref": str, "alt": str, "assembly": str}]
    """
    if assembly is None:
        assembly = DBSNP_DEFAULT_ASSEMBLY

    data = fetch_snp_json(rsid)
    try:
        placements = data["primary_snapshot_data"]["placements_with_allele"]
    except:
        return None

    bucket = []
    for plc in placements:
        chrom = _accession_to_chrom(plc["seq_id"])
        # Check if chrom is a valid chromosome (1-22, X, or Y)
        if chrom in ["X", "Y"]:
            pass  # Valid sex chromosome
        else:
            try:
                chrom_num = int(chrom)
                if not (1 <= chrom_num <= 22):
                    continue  # Skip if not in range 1-22
            except ValueError:
                continue  # Skip if chrom is not a number or X/Y
        
        # Get assembly information for this placement
        placement_assembly = None
        if 'placement_annot' in plc and 'seq_id_traits_by_assembly' in plc['placement_annot']:
            for assembly_info in plc['placement_annot']['seq_id_traits_by_assembly']:
                placement_assembly = assembly_info.get('assembly_name', '')
                break  # Use the first assembly name found
        
        # Skip if assembly doesn't match filter
        if placement_assembly and not _matches_assembly(placement_assembly, assembly):
            continue
        
        for al in plc["alleles"]:
            spdi = al["allele"]["spdi"]
            ref, alt = spdi["deleted_sequence"], spdi["inserted_sequence"]
            if ref == alt:  # Skip reference allele
                continue  
            bucket.append({
                "chrom": chrom,
                "pos": spdi["position"] + 1,  # Convert from 0-based to 1-based
                "ref": ref,
                "alt": alt,
                "assembly": placement_assembly or "Unknown"
            })

    return bucket


def extract_allele_frequencies(data: dict) -> List[Dict]:
    """
    Extract allele frequency data from dbSNP JSON response.
    
    Args:
        data: Raw JSON response from dbSNP API
    
    Returns:
        List of dictionaries containing frequency information from different studies
    """
    frequencies = []
    
    if 'primary_snapshot_data' not in data:
        return frequencies
    
    # Extract from all allele annotations
    allele_annotations = data['primary_snapshot_data'].get('allele_annotations', [])
    
    for annotation in allele_annotations:
        if 'frequency' not in annotation:
            continue
            
        freq_data = annotation['frequency']
        
        for study in freq_data:
            study_info = {
                "study_name": study.get("study_name", "Unknown"),
                "study_version": study.get("study_version", ""),
                "allele_count": study.get("allele_count", 0),
                "total_count": study.get("total_count", 0),
                "allele_frequency": 0.0
            }
            
            # Calculate allele frequency
            if study_info["total_count"] > 0:
                study_info["allele_frequency"] = study_info["allele_count"] / study_info["total_count"]
            
            # Add observation details if available
            if "observation" in study:
                obs = study["observation"]
                study_info["observation"] = {
                    "seq_id": obs.get("seq_id", ""),
                    "position": obs.get("position", 0),
                    "deleted_sequence": obs.get("deleted_sequence", ""),
                    "inserted_sequence": obs.get("inserted_sequence", "")
                }
            
            frequencies.append(study_info)
    
    return frequencies


def get_variant_info(rsid: str, assembly: str = None) -> Dict:
    """
    Get basic variant information for an rsID.
    
    Args:
        rsid: The rsID to query (e.g., "rs1333049" or "1333049")
        assembly: Genome assembly version (default: config.DBSNP_DEFAULT_ASSEMBLY)
    
    Returns:
        Dictionary with variant information including coordinates and allele frequencies
    """
    if assembly is None:
        assembly = DBSNP_DEFAULT_ASSEMBLY
        
    rsid_with_prefix = rsid if rsid.lower().startswith("rs") else f"rs{rsid}"
    raw_data = fetch_snp_json(rsid)
    coords = get_variant_coordinates(rsid, assembly)
    frequencies = extract_allele_frequencies(raw_data)
    return {
        "rsid": rsid_with_prefix,
        "coordinates": coords,
        "allele_frequencies": frequencies,
        "assembly_filter": assembly,
        "raw_data": raw_data
    }

if __name__ == "__main__":
    # Test with different assemblies
    print("=== Testing GRCh38 (default) ===")
    info_38 = get_variant_info("rs138188004", "GRCh38")
    print(f"Coordinates: {len(info_38['coordinates'])}")
    for i, c in enumerate(info_38['coordinates']):
        print(f"{i}: {c}")
    
    print(f"\nAllele frequencies: {len(info_38['allele_frequencies'])}")
    for i, freq in enumerate(info_38['allele_frequencies'][:5]):  # Show first 5
        print(f"{i}: {freq['study_name']} - AF: {freq['allele_frequency']:.4f} ({freq['allele_count']}/{freq['total_count']})")
    
    print("\n=== Testing GRCh37 ===")
    info_37 = get_variant_info("rs138188004", "GRCh37")
    print(f"Coordinates: {len(info_37['coordinates'])}")
    print(f"Allele frequencies: {len(info_37['allele_frequencies'])}")
    
    print("\n=== Testing ALL assemblies ===")
    info_all = get_variant_info("rs138188004", "all")
    print(f"Coordinates: {len(info_all['coordinates'])}")
    print(f"Allele frequencies: {len(info_all['allele_frequencies'])}")