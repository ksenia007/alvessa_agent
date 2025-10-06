# src/tools/chembl/utils.py
# Author: Dmitri Kosenkov
# Created: 2025-09-25
# Updated: 2025-10-06
#
# Optimized utilities for the agentic ChEMBL drug-target summarization tool.
# Includes canonical SMILES for frontend RDKit.js 2D rendering.
# Also retains molecule_type and InChIKey for PubChem enrichment.
# Provides standalone HTML injection utilities (mirrors prot tool).

import sqlite3
from typing import List, Dict
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import json

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

    Returns dict with keys:
      - approved_drugs: list of (chembl_id, molecule_type, inchikey, smiles, first_approval)
      - clinical_trials: list of (chembl_id, phase, molecule_type, inchikey, smiles)
      - bioactivity: list of (chembl_id, molecule_type, summary_str, inchikey, smiles)
    """
    data = {"approved_drugs": [], "clinical_trials": [], "bioactivity": []}
    cur = conn.cursor()

    # --- 1. Resolve UniProt ID -> target ID ---
    cur.execute("""
        SELECT td.tid
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

    # --- 2. FDA-approved drugs ---
    cur.execute("""
        SELECT DISTINCT
            m.chembl_id,
            m.molecule_type,
            cs.standard_inchi_key,
            cs.canonical_smiles,
            m.first_approval
        FROM drug_mechanism dm
        JOIN molecule_dictionary m ON dm.molregno = m.molregno
        LEFT JOIN compound_structures cs ON m.molregno = cs.molregno
        WHERE dm.tid = ? AND m.max_phase = 4
        ORDER BY
            CASE
                WHEN m.molecule_type = 'Small molecule' THEN 1
                WHEN m.molecule_type = 'Antibody' THEN 2
                ELSE 3
            END,
            COALESCE(m.first_approval, 0) DESC,
            m.chembl_id
    """, (tid,))
    data["approved_drugs"] = [
        (
            r["chembl_id"],
            r["molecule_type"],
            r["standard_inchi_key"] or "",
            r["canonical_smiles"] or "",
            r["first_approval"] or ""
        )
        for r in cur.fetchall()
    ]

    # --- 3. Clinical trials (phases 1–3) ---
    cur.execute("""
        SELECT DISTINCT
            m.chembl_id,
            m.max_phase,
            m.molecule_type,
            cs.standard_inchi_key,
            cs.canonical_smiles
        FROM drug_mechanism dm
        JOIN molecule_dictionary m ON dm.molregno = m.molregno
        LEFT JOIN compound_structures cs ON m.molregno = cs.molregno
        WHERE dm.tid = ? AND m.max_phase BETWEEN 1 AND 3
        ORDER BY
            m.max_phase DESC,
            CASE
                WHEN m.molecule_type = 'Small molecule' THEN 1
                WHEN m.molecule_type = 'Antibody' THEN 2
                ELSE 3
            END,
            m.chembl_id
    """, (tid,))
    data["clinical_trials"] = [
        (
            r["chembl_id"],
            r["max_phase"],
            r["molecule_type"],
            r["standard_inchi_key"] or "",
            r["canonical_smiles"] or "",
        )
        for r in cur.fetchall()
    ]

    # --- 4. Preclinical (max_phase=0) ---
    cur.execute("""
        SELECT DISTINCT
            m.chembl_id,
            m.molecule_type,
            cs.standard_inchi_key,
            cs.canonical_smiles
        FROM drug_mechanism dm
        JOIN molecule_dictionary m ON dm.molregno = m.molregno
        LEFT JOIN compound_structures cs ON m.molregno = cs.molregno
        WHERE dm.tid = ? AND m.max_phase = 0
    """, (tid,))
    preclinical = [
        (
            r["chembl_id"],
            r["molecule_type"],
            r["standard_inchi_key"] or "",
            r["canonical_smiles"] or "",
        )
        for r in cur.fetchall()
    ]
    if preclinical:
        data["clinical_trials"].extend([(d, 0, t, k, s) for d, t, k, s in preclinical])

    # --- 5. Bioactivity assays ---
    cur.execute("""
        SELECT
            m.chembl_id,
            m.molecule_type,
            cs.standard_inchi_key,
            cs.canonical_smiles,
            a.standard_type,
            a.standard_value,
            a.standard_units,
            a.pchembl_value
        FROM activities a
        JOIN molecule_dictionary m ON a.molregno = m.molregno
        JOIN assays ass ON a.assay_id = ass.assay_id
        LEFT JOIN compound_structures cs ON m.molregno = cs.molregno
        WHERE ass.tid = ?
          AND a.standard_type IN ('IC50','EC50','Ki')
          AND a.standard_value > 0
          AND a.standard_units IN ('nM','uM','pM','M')
    """, (tid,))

    assay_data = defaultdict(
        lambda: {"values": [], "pchembl": [], "mtype": None, "inchikey": "", "smiles": ""}
    )

    for r in cur.fetchall():
        nm_value = _to_nM(r["standard_value"], r["standard_units"])
        if nm_value is None:
            continue
        chembl_id = r["chembl_id"]
        stype = r["standard_type"]
        key = (chembl_id, stype)
        assay_data[key]["values"].append(nm_value)
        assay_data[key]["mtype"] = r["molecule_type"]
        assay_data[key]["inchikey"] = r["standard_inchi_key"] or ""
        assay_data[key]["smiles"] = r["canonical_smiles"] or ""
        if r["pchembl_value"] is not None:
            assay_data[key]["pchembl"].append(r["pchembl_value"])

    for (chembl_id, stype), vdict in assay_data.items():
        best_value = min(vdict["values"])
        n_assays = len(vdict["values"])
        line = f"{stype}={best_value:.2f} nM (n={n_assays})"
        if vdict["pchembl"]:
            line += f", pChEMBL={max(vdict['pchembl']):.2f}"
        data["bioactivity"].append(
            (
                chembl_id,
                vdict["mtype"],
                line,
                vdict["inchikey"],
                vdict["smiles"],
            )
        )

    try:
        data["bioactivity"].sort(key=lambda x: float(x[2].split("=")[1].split()[0]))
    except Exception:
        log("Warning: could not sort bioactivity safely")

    if limit and len(data["bioactivity"]) > limit:
        data["bioactivity"] = data["bioactivity"][:limit]

    return data


# ------------------------------
# Summary formatting
# ------------------------------
def _fmt_entry(chembl_id: str, mtype: str) -> str:
    """Format ChEMBL ID, adding molecule type if non-small molecule."""
    return chembl_id if mtype == "Small molecule" else f"{chembl_id} [Molecule type: {mtype}]"


def make_summary_text(gene: str, uniprot_id: str, target_data: Dict[str, List]) -> str:
    lines = [f"Gene: {gene}", f"  UniProt ID: {uniprot_id}"]

    if target_data["approved_drugs"]:
        lines.append("  FDA-approved drugs:")
        for d, t, *_ in target_data["approved_drugs"]:
            lines.append(f"    - {_fmt_entry(d, t)}")

    if target_data["clinical_trials"]:
        lines.append("  Clinical/preclinical trials:")
        for d, p, t, *_ in target_data["clinical_trials"]:
            phase_str = "preclinical" if p == 0 else f"phase {p}"
            lines.append(f"    - {_fmt_entry(d, t)} ({phase_str})")

    if target_data["bioactivity"]:
        lines.append("  Bioactivity evidence (normalized to nM):")
        for d, t, rest, *_ in target_data["bioactivity"]:
            lines.append(f"    - {_fmt_entry(d, t)} {rest}")

    if len(lines) == 2:
        lines.append("  No ChEMBL druggability data found.")

    return "\n".join(lines)


# ------------------------------
# Interpretation notes
# ------------------------------
def interpretation_notes(include_bioactivity: bool = True) -> str:
    notes = ["\nInterpretation Notes\n"]
    if include_bioactivity:
        notes.append(
            "- Bioactivity values are normalized to nM and aggregated per compound.\n"
            "- Only positive numeric IC50 / Ki / EC50 values with molar units were considered.\n"
            "- Reported values are the lowest (strongest) potency observed per assay type.\n"
            "- pChEMBL is -log10(molar potency): higher values = stronger activity.\n"
            "- Molecule type is taken directly from ChEMBL (Small molecule, Protein, Antibody, etc.).\n"
            "- Canonical SMILES are provided for RDKit.js 2D rendering.\n"
            "- InChIKeys are provided for PubChem enrichment.\n"
            "- Schema version = 35.0 (local snapshot).\n"
            "- These summaries are heuristic and not a substitute for experimental validation.\n"
        )
    return "".join(notes)


# ------------------------------
# Frontend asset injection
# ------------------------------
def inject_frontend_assets(
    html_template: str,
    css_template: str,
    js_template: str,
    chembl_data_all: Dict[str, Dict],
    genes_str: str
) -> str:
    """
    Inject inline CSS, JS, and JSON ChemBL data into HTML template.
    Produces a standalone HTML viewer (parallel to prot tool).
    """
    css_inline = f"<style>\n{css_template}\n</style>"
    js_inline = (
        "<script>\n"
        f"var chemblData = {json.dumps(chembl_data_all, indent=2)};\n"
        f"{js_template}\n</script>"
    )

    return (
        html_template
        .replace("{{CSS_INLINE}}", css_inline)
        .replace("{{JS_INLINE}}", js_inline)
        .replace("{{GENE_SYMBOL}}", genes_str)
    )
