"""Global settings and constants shared across the package."""

import operator
import os

# -----------------------------------------------------------------------------
# Debug & runtime limits
# -----------------------------------------------------------------------------
DEBUG: bool = True
N_CHARS: int = 100_000
RATE_LIMIT: int = 60  # Claude calls per rolling minute

# -----------------------------------------------------------------------------
# dbSNP Configuration
# -----------------------------------------------------------------------------
# Default genome assembly for dbSNP queries
# Options: "GRCh38", "GRCh37", "all" (returns coordinates from all assemblies)
DBSNP_DEFAULT_ASSEMBLY: str = "GRCh38"

# Assembly mapping for filtering dbSNP results
DBSNP_ASSEMBLY_PATTERNS = {
    "GRCh38": ["GRCh38", "hg38"],
    "GRCh37": ["GRCh37", "hg19", "b37"],
    "all": []  # Empty list means no filtering
}

# -----------------------------------------------------------------------------
# Anthropic model versions
# -----------------------------------------------------------------------------
GENE_EXTRACT_MODEL: str = "claude-3-haiku-20240307"
CONDITIONED_MODEL: str = "claude-sonnet-4-20250514"
VERIFY_MODEL: str = "claude-3-haiku-20240307"
TOOL_SELECTOR_MODEL: str = "claude-3-haiku-20240307"

ANTHROPIC_API_KEY: str = os.environ["ANTHROPIC_API_KEY"]
BioGRID_API_KEY: str = os.environ["BIOGRID_API_KEY"]
