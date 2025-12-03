"""
DrugCentral Drug-Centric Tool Package
Centralized constants and shared paths for the `drug_central` tool.
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

# --- Database File ---
DRUG_CENTRAL_DB_PATH = LOCAL_DBS_DIR / "drug_central.dump.11012023.db"

# --- Frontend Templates (for standalone Drug Central viewer) ---
HTML_TEMPLATE = STATIC_DIR / "tool_drug_central.html"
CSS_TEMPLATE = STATIC_DIR / "tool_drug_central.css"
JS_TEMPLATE = STATIC_DIR / "tool_drug_central.js"

# --- Database validation metadata (minimal for now) ---
REQUIRED_META = {
    # Drug Central dbversion.version should be >= this if you want a sanity check.
    # Left as a loose requirement for now; you can tighten later.
    "min_version": 20231101,
}
UPDATE_MSG = (
    "Please update to the latest Drug Central SQLite dump file: "
    f"{DRUG_CENTRAL_DB_PATH.name}"
)

# Ensure output dir exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Expansion behavior is controlled by src.tools.drug_central.DRUG_TARGET_RETURN_MODE:
#
#      - "none" -> do not resolve or return any gene targets for drugs
#                  (ChemBL and DrugCentral both skipped).
#      - "top"  -> resolve ChemBL + DrugCentral, but keep only the top
#                  DRUG_TARGET_TOP_N DrugCentral HUMAN targets (by score).
#      - "all"  -> resolve ChemBL + DrugCentral and return all targets.
DRUG_TARGET_RETURN_MODE = "none"

# Sort DrugCentral targets by descending score, keep top N
DRUG_TARGET_TOP_N = 5

__all__ = [
    # Paths
    "PACKAGE_ROOT",
    "REPO_ROOT",
    "LOCAL_DBS_DIR",
    "OUTPUT_DIR",
    "WEB_DIR",
    "STATIC_DIR",
    # Database + templates
    "DRUG_CENTRAL_DB_PATH",
    "HTML_TEMPLATE",
    "CSS_TEMPLATE",
    "JS_TEMPLATE",
    # Metadata
    "REQUIRED_META",
    "UPDATE_MSG",
    #Gene extraction mode
    "DRUG_TARGET_RETURN_MODE",
    "DRUG_TARGET_TOP_N",
]
