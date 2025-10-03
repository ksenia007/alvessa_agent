# src/tools/chembl/__init__.py
"""
ChEMBL Drug-Target Summarization Tool Package
Centralized constants and shared paths for the `chembl` tool.
"""

from pathlib import Path

# --- Package/Repo Roots ---
PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]

# --- Local Storage Layout ---
LOCAL_DBS_DIR = REPO_ROOT / "local_dbs"
OUTPUT_DIR = REPO_ROOT / "out"

# --- Files ---
CHEMBL_DB_PATH = LOCAL_DBS_DIR / "chembl_35.db"

# --- Database validation metadata ---
REQUIRED_META = {
    "schema_version": "35.0",
}
UPDATE_MSG = f"Please update to the latest {CHEMBL_DB_PATH.name} file."

# Ensure output dir exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

__all__ = [
    "PACKAGE_ROOT", "REPO_ROOT",
    "LOCAL_DBS_DIR", "OUTPUT_DIR",
    "CHEMBL_DB_PATH",
]
