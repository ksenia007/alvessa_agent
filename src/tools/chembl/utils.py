# src/tools/chembl/utils.py
# Author: Dmitri Kosenkov
# Created: 2025-09-25
# Updated: 2025-10-01
#
# Optimized utilities for the agentic ChEMBL drug-target summarization tool.
# Revised to include molecule_type and mark non-small molecules explicitly.

import sqlite3
from typing import List, Dict, Tuple
from pathlib import Path
from datetime import datetime
from collections import defaultdict

from . import CHEMBL_DB_PATH, UPDATE_MSG


# ------------------------------
# Logging
# ------------------------------
def log(msg: str) -> None:
    print(f"[ChEMBL] {msg} @ {datetime.now()}")


# ------------------------------
# DB connection
# ------------------------------
def get_connection() -> sqlite3.Connection:
    if not CHEMBL_DB_PATH.exists():
        raise RuntimeError(f"Database missing: {CHEMBL_DB_PATH}. {UPDATE_MSG}")
    conn = sqlite3.connect(CHEMBL_DB_PATH, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    return conn


# ------------------------------
# Unit normalization
# ------------------------------
def _to_nM(value: float, unit: str) -> float:
    """
    Convert assay value into nM for consistency.
    Returns None if the unit is unsupported.
    """
    if value is None or value <= 0:
        return None
    unit = (unit or "").lower()
    if unit == "nm":
        return value
    elif unit == "um":
        return value * 1000.0
    elif unit == "pm":
        return value / 1000.0
    elif unit == "m":
        return value * 1e9
    else:
        return None


# ------------------------------
# Queries
# ------------------------------
def fetch_target_info(
    conn: sqlite3.Connection,
    uniprot_id: str,
    limit: int = 25
) -> Dict[str, List]:
    """
    Query ChEMBL for a UniProt target: approved drugs, clinical/preclinical trials,
    and assay bioactivities (normalized to nM).
    Returns dict with keys: approved_drugs, clinical_trials, bioactivity.
    Each entry includes molecule_type; if not 'Small molecule', it will be labeled.
    """
    data = {
        "approved_drugs": [],      # list of (chembl_id, molecule_type)
        "clinical_trials": [],     # list of (chembl_id, phase, molecule_type)
        "bioactivity": [],         # list of (chembl_id, molecule_type, summary_str)
    }
    cur = conn.cursor()

    # --- 1. Resolve UniProt ID -> ChEMBL target ---
    cur.execute("""
        SELECT td.chembl_id AS target_chembl_id, td.tid
        FROM target_dictionary td
        JOIN target_components tc ON td.tid = tc.tid
        JOIN component_sequences cs ON tc.component_id = cs.component_id
        WHERE cs.accession = ?
        LIMIT 1
    """, (uniprot_id,))
    row = cur.fetchone()
    if not row:
        return data

    tid = row["tid"]

    # --- 2. FDA-approved drugs (max_phase=4) ---
    cur.execute("""
        SELECT DISTINCT m.chembl_id, m.molecule_type
        FROM drug_mechanism dm
        JOIN molecule_dictionary m ON dm.molregno = m.molregno
        WHERE dm.tid=? AND m.max_phase=4
    """, (tid,))
    data["approved_drugs"] = [(r["chembl_id"], r["molecule_type"]) for r in cur.fetchall()]

    # --- 3. Clinical trial drugs (phases 1-3) ---
    cur.execute("""
        SELECT DISTINCT m.chembl_id, m.max_phase, m.molecule_type
        FROM drug_mechanism dm
        JOIN molecule_dictionary m ON dm.molregno = m.molregno
        WHERE dm.tid=? AND m.max_phase BETWEEN 1 AND 3
    """, (tid,))
    data["clinical_trials"] = [(r["chembl_id"], r["max_phase"], r["molecule_type"])
                               for r in cur.fetchall()]

    # --- 4. Preclinical compounds (phase=0) ---
    cur.execute("""
        SELECT DISTINCT m.chembl_id, m.molecule_type
        FROM drug_mechanism dm
        JOIN molecule_dictionary m ON dm.molregno = m.molregno
        WHERE dm.tid=? AND m.max_phase=0
    """, (tid,))
    preclinical = [(r["chembl_id"], r["molecule_type"]) for r in cur.fetchall()]
    if preclinical:
        data["clinical_trials"].extend([(d, 0, t) for d, t in preclinical])

    # --- 5. Bioactivity assays (normalized to nM) ---
    cur.execute("""
        SELECT m.chembl_id,
               m.molecule_type,
               a.standard_type,
               a.standard_value,
               a.standard_units,
               a.pchembl_value
        FROM activities a
        JOIN molecule_dictionary m ON a.molregno = m.molregno
        JOIN assays ass ON a.assay_id = ass.assay_id
        WHERE ass.tid=?
          AND a.standard_type IN ('IC50','EC50','Ki')
          AND a.standard_value > 0
          AND a.standard_units IN ('nM','uM','pM','M')
    """, (tid,))

    assay_data = defaultdict(lambda: {"values": [], "pchembl": [], "mtype": None})
    for r in cur.fetchall():
        nm_value = _to_nM(r["standard_value"], r["standard_units"])
        if nm_value is None:
            continue
        chembl_id = r["chembl_id"]
        stype = r["standard_type"]
        key = (chembl_id, stype)
        assay_data[key]["values"].append(nm_value)
        assay_data[key]["mtype"] = r["molecule_type"]
        if r["pchembl_value"] is not None:
            assay_data[key]["pchembl"].append(r["pchembl_value"])

    for (chembl_id, stype), vdict in assay_data.items():
        best_value = min(vdict["values"])
        n_assays = len(vdict["values"])
        line = f"{stype}={best_value:.2f} nM (n={n_assays})"
        if vdict["pchembl"]:
            line += f", pChEMBL={max(vdict['pchembl']):.2f}"
        data["bioactivity"].append((chembl_id, vdict["mtype"], line))

    # Sort strongest first (lowest potency value)
    data["bioactivity"].sort(key=lambda x: float(x[2].split("=")[1].split()[0]))
    if limit and len(data["bioactivity"]) > limit:
        data["bioactivity"] = data["bioactivity"][:limit]

    return data


# ------------------------------
# Summary formatting
# ------------------------------
def _fmt_entry(chembl_id: str, mtype: str) -> str:
    """Format a ChEMBL ID with molecule type if not a small molecule."""
    return chembl_id if mtype == "Small molecule" else f"{chembl_id} [Molecule type: {mtype}]"


def make_summary_text(gene: str, uniprot_id: str, target_data: Dict[str, List]) -> str:
    lines = [f"Gene: {gene}"]
    lines.append(f"  UniProt ID: {uniprot_id}")

    if target_data["approved_drugs"]:
        lines.append("  FDA-approved drugs:")
        lines.extend([f"    - {_fmt_entry(d, t)}" for d, t in target_data["approved_drugs"]])

    if target_data["clinical_trials"]:
        lines.append("  Clinical/preclinical trials:")
        for d, p, t in target_data["clinical_trials"]:
            phase_str = "preclinical" if p == 0 else f"phase {p}"
            lines.append(f"    - {_fmt_entry(d, t)} ({phase_str})")

    if target_data["bioactivity"]:
        lines.append("  Bioactivity evidence (normalized to nM, best potency per compound):")
        for d, t, rest in target_data["bioactivity"]:
            lines.append(f"    - {_fmt_entry(d, t)} {rest}")

    if len(lines) == 2:
        lines.append("  No ChEMBL druggability data found.")

    return "\n".join(lines)


# ------------------------------
# Interpretation notes
# ------------------------------
def interpretation_notes(include_bioactivity: bool = True) -> str:
    sections: List[str] = ["\nInterpretation Notes\n"]
    if include_bioactivity:
        sections.append(
            "- Bioactivity values are normalized to nM and aggregated per compound/assay type.\n"
            "- Only positive numeric IC50 / Ki / EC50 values with molar units were considered.\n"
            "- Reported values are the lowest (strongest) potency observed, with count of assays (n).\n"
            "- pChEMBL is -log10(molar potency): higher values = stronger activity.\n"
            "    * >=9  (~1 nM): very high potency\n"
            "    * 7-9 (~100 nM - 10 nM): strong potency\n"
            "    * 5-7 (~10 uM - 100 nM): moderate potency\n"
            "    * <5  : weak potency\n"
            "- Larger assay counts (n) suggest stronger evidence; single assay results may be less reliable.\n"
            "- Phase annotation comes from ChEMBL molecule max_phase field (0=preclinical, 1-3=clinical, 4=approved).\n"
            "- Molecule type is taken from ChEMBL (Small molecule, Protein, Antibody, Oligosaccharide, etc.).\n"
            "- Data were extracted from the local ChEMBL database snapshot.\n"
            "  Schema version=35.0, build_time=NA.\n"
            "- Assumptions: only molar-based assays included; units like ug/mL excluded.\n"
            "- These summaries are heuristic and not a substitute for detailed pharmacological evaluation.\n"
            "References:\n"
            "  - Davies M. et al., Nucleic Acids Res. 2015;43(Database issue):D1109-D1116.\n"
            "  - Gaulton A. et al., Nucleic Acids Res. 2017;45(Database issue):D945-D954.\n"
            "  - Bento A.P. et al., Nucleic Acids Res. 2014;42(Database issue):D1083-D1090.\n"
        )
    return "".join(sections)
