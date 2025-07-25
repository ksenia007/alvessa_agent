import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import re
from typing import List, Dict


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


def get_variant_coordinates(rsid: str, assembly: str = "GRCh38") -> List[Dict[str, str]]:
    """
    Get genomic coordinates for an rsID.
    
    Args:
        rsid: The rsID to query (e.g., "rs1333049" or "1333049")
        assembly: Genome assembly version (default: "GRCh38")
    
    Returns:
        List of dictionaries with variant coordinates:
        [{"chrom": str, "pos": int, "ref": str, "alt": str, "assembly": str}]
    """
    def _same_major(a, b) -> bool:
        """Compare assembly names ignoring patch version."""
        if a is None or b is None:
            return True
        return a.split(".")[0] == b.split(".")[0]

    data = fetch_snp_json(rsid)
    placements = data["primary_snapshot_data"]["placements_with_allele"]

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
            })

    # Return the first non-empty bucket
    return bucket


def get_variant_info(rsid: str) -> Dict:
    """
    Get basic variant information for an rsID.
    
    Args:
        rsid: The rsID to query (e.g., "rs1333049" or "1333049")
    
    Returns:
        Dictionary with variant information including coordinates
    """
    rsid_with_prefix = rsid if rsid.lower().startswith("rs") else f"rs{rsid}"
    coords = get_variant_coordinates(rsid)

    return {
        "rsid": rsid_with_prefix,
        "coordinates": coords,
        "raw_data": fetch_snp_json(rsid)
    }

if __name__ == "__main__":
    # Get all basic info
    info = get_variant_info("rs138188004")

    # Get just coordinates  
    coords = get_variant_coordinates("rs138188004")

    # Get raw API data
    raw_data = fetch_snp_json("rs138188004")
    
    print(info['raw_data'].keys())

    print("-"*100)
    for i, c in enumerate(coords):
        print(f"{i}: {c}")
    # print(raw_data)