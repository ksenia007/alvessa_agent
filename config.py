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
# Anthropic model versions
# -----------------------------------------------------------------------------
GENE_EXTRACT_MODEL: str = "claude-3-haiku-20240307"
CONDITIONED_MODEL: str = "claude-sonnet-4-20250514"
VERIFY_MODEL: str = "claude-3-haiku-20240307"
TOOL_SELECTOR_MODEL: str = "claude-3-haiku-20240307"
GLINER_CONFIG: dict = {
    "model_name": "urchade/gliner_medium-v2.1",
    "threshold": 0.5,
    "labels": ["Gene", "Disease", "Trait", "Protein Function", "GO terms", "Protein", "Genetic Variant"]
}

ANTHROPIC_API_KEY: str = os.environ["ANTHROPIC_API_KEY"]
BioGRID_API_KEY: str = os.environ["BioGRID_API_KEY"]
