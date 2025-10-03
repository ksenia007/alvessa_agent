# src/tools/prot/__init__.py
"""
Protein Visualization & Summarization Tool Package
Centralized constants, thresholds, and shared paths for the `prot` tool.
"""
from pathlib import Path

# --- Package/Repo Roots ---
PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]

# --- Local Storage Layout ---
LOCAL_DBS_DIR = REPO_ROOT / "local_dbs"
OUTPUT_DIR = REPO_ROOT / "out"
WEB_DIR = REPO_ROOT / "web"
STATIC_DIR = WEB_DIR / "static"

# --- Files ---
DB_PATH = LOCAL_DBS_DIR / "alvessa_proteins.db"
PDB_DIR = LOCAL_DBS_DIR / "pdb"
HTML_TEMPLATE = STATIC_DIR / "tool_prot.html"
CSS_TEMPLATE = STATIC_DIR / "tool_prot.css"
JS_TEMPLATE = STATIC_DIR / "tool_prot.js"

# --- External Datasets ---
DISPROT_JSON = LOCAL_DBS_DIR / "DisProt release_2025_06 with_ambiguous_evidences.json"

# --- Database validation metadata ---
REQUIRED_META = {
    "schema_version": "1.0",
    "build_time": "2025-09-26T18:53:44Z",
}
UPDATE_MSG = f"Please update to the latest {DB_PATH.name} file."

# --- Thresholds & Cutoffs ---
PLDDT_HIGH_CUTOFF: float = 90.0
PLDDT_MEDIUM_CUTOFF: float = 70.0
PLDDT_LOW_WARNING: float = 70.0   # threshold for warnings

# --- Formatting limits ---
BIO_RESIDUE_PREVIEW_MAX: int = 20     # residues shown in preview before truncation
BIO_EXCLUDED_IDS_MAX: int = 10        # excluded ligand IDs shown before truncation
FMT_FLOAT_NDIGITS: int = 2
FMT_FLOAT_DEFAULT: str = "NA"

# Ensure output dir exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

__all__ = [
    # Roots & storage
    "PACKAGE_ROOT", "REPO_ROOT",
    "LOCAL_DBS_DIR", "OUTPUT_DIR",
    "WEB_DIR", "STATIC_DIR",
    # Files
    "DB_PATH", "PDB_DIR",
    "HTML_TEMPLATE", "CSS_TEMPLATE", "JS_TEMPLATE",
    "DISPROT_JSON",
    # Metadata
    "REQUIRED_META", "UPDATE_MSG",
    # Thresholds
    "PLDDT_HIGH_CUTOFF", "PLDDT_MEDIUM_CUTOFF", "PLDDT_LOW_WARNING",
    # Formatting
    "BIO_RESIDUE_PREVIEW_MAX", "BIO_EXCLUDED_IDS_MAX",
    "FMT_FLOAT_NDIGITS", "FMT_FLOAT_DEFAULT",
]
