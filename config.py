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

ANTHROPIC_API_KEY: str = os.environ["ANTHROPIC_API_KEY"]
