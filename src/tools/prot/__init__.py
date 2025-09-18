# src/tools/prot/__init__.py
"""
Protein Visualization & Summarization Tool Package
Centralized constants and shared paths for the `prot` tool.
"""
from pathlib import Path
from datetime import datetime

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

# --- Database freshness requirement ---
DB_MIN_MTIME = datetime(2025, 9, 10)

# Ensure output dir exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

__all__ = [
    "PACKAGE_ROOT", "REPO_ROOT",
    "LOCAL_DBS_DIR", "OUTPUT_DIR",
    "WEB_DIR", "STATIC_DIR",
    "DB_PATH", "PDB_DIR",
    "HTML_TEMPLATE", "CSS_TEMPLATE", "JS_TEMPLATE",
    "DB_MIN_MTIME"
]
