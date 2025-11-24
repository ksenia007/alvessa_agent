# src/tools/aa_seq/__init__.py
# ===========================================================
#  AA Sequence/ Gene Tool Package
# ===========================================================
#
# Author: Dmitri Kosenkov
# Created: 2025-11-20
#
# This package enabels AA sequence-to-gene resolution
# using a local copy of UniProtKB (fast).
#
#
# Expected SQLite schema (built outside of Alvessa):
#
#   Table: uniprot
#     - acc            (TEXT PRIMARY KEY)  -- full accession (e.g., P00533-2)
#     - canonical_acc  (TEXT NOT NULL)     -- root accession (e.g., P00533)
#     - sequence       (TEXT NOT NULL)     -- amino acid sequence
#     - gene_name      (TEXT)              -- primary gene symbol
#     - entrez_gene_id (TEXT)              -- Entrez Gene ID (optional)
#
#   Table: meta
#     - key   (TEXT PRIMARY KEY)
#     - value (TEXT NOT NULL)
#       * db_version    : integer-like string, e.g. 20251118
#       * tsv_source    : TSV filename
#       * fasta_sources : comma-separated FASTA filenames
#       * n_records     : row count in `uniprot`
#

from __future__ import annotations

from pathlib import Path
import sqlite3

# --- Package / Repo Roots ----------------------------------------------------

PACKAGE_ROOT = Path(__file__).resolve().parent
# src/tools/uniprot/__init__.py -> REPO_ROOT = project root (alvessa_agent)
REPO_ROOT = PACKAGE_ROOT.parents[2]

# --- Local Storage Layout ----------------------------------------------------

# Matches your actual tree:
#   <REPO_ROOT>/local_dbs
#   <REPO_ROOT>/out
#   <REPO_ROOT>/web/static
LOCAL_DBS_DIR = REPO_ROOT / "local_dbs"
OUTPUT_DIR = REPO_ROOT / "out"

# Ensure output dir exists (for any future reports / logs)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# --- UniProt Database File ---------------------------------------------------
#
# This should correspond to the DB built by uniprot_db_build.py
# (DB_PATH_DEFAULT = data/uniprot_hs_9606_reviewed.db), but
# relocated into <REPO_ROOT>/local_dbs.
#
UNIPROT_DB_PATH = LOCAL_DBS_DIR / "uniprot_hs_9606_reviewed.db"

# --- Database validation metadata (aligned with uniprot_db_build.py) ---------

REQUIRED_META = {
    # Must match or exceed the min_version used during DB build
    # (see uniprot_db_build.REQUIRED_META).
    "min_version": 20251101,
}

UPDATE_MSG = (
    "Please ensure the UniProt SQLite database exists at:\n"
    f"  {UNIPROT_DB_PATH}\n"
    "and that it was built with the expected schema (tables `uniprot` and `meta`)."
)


def validate_uniprot_db(db_path: Path | None = None) -> None:
    """
    Sanity check for the UniProt SQLite database.

    Validates:
      - DB file exists
      - Table `uniprot` is present
      - Table `meta` is present
      - `meta.db_version` exists and is >= REQUIRED_META['min_version']

    Raises:
        FileNotFoundError if the DB file is missing.
        RuntimeError if required tables/metadata are missing or too old.
    """
    path = db_path or UNIPROT_DB_PATH

    if not path.is_file():
        raise FileNotFoundError(
            f"UniProt SQLite DB not found at: {path}\n{UPDATE_MSG}"
        )

    conn = sqlite3.connect(path)
    try:
        cur = conn.cursor()

        # Check `uniprot` table
        cur.execute(
            "SELECT name FROM sqlite_master "
            "WHERE type='table' AND name='uniprot'"
        )
        row = cur.fetchone()
        if row is None:
            raise RuntimeError(
                f"Table 'uniprot' not found in DB: {path}\n"
                "Please rebuild the DB using your UniProt build pipeline."
            )

        # Check `meta` table
        cur.execute(
            "SELECT name FROM sqlite_master "
            "WHERE type='table' AND name='meta'"
        )
        row = cur.fetchone()
        if row is None:
            raise RuntimeError(
                f"Table 'meta' not found in DB: {path}\n"
                "The DB should be built by uniprot_db_build.py with provenance metadata."
            )

        # Check db_version in meta
        cur.execute("SELECT value FROM meta WHERE key='db_version'")
        row = cur.fetchone()
        if row is None or row[0] is None:
            raise RuntimeError(
                f"'meta.db_version' missing in DB: {path}\n"
                "Please rebuild the DB using uniprot_db_build.py."
            )

        db_version = int(row[0])
        min_version = int(REQUIRED_META.get("min_version", 0))

        if db_version < min_version:
            raise RuntimeError(
                f"UniProt DB version {db_version} is older than the required "
                f"min_version {min_version}. "
                "Please rebuild the UniProt DB from the latest UniProtKB files."
            )

    finally:
        conn.close()


__all__ = [
    # Paths
    "PACKAGE_ROOT",
    "REPO_ROOT",
    "LOCAL_DBS_DIR",
    "OUTPUT_DIR",
    # Database
    "UNIPROT_DB_PATH",
    "REQUIRED_META",
    "UPDATE_MSG",
    # Helpers
    "validate_uniprot_db",
]
