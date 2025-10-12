# src/tools/chembl/utils.py
# ===========================================================
# ChEMBL Drug-Target Evidence Utilities with OpenFDA Integration
# ===========================================================
#
# Author: Dmitri Kosenkov
# Created: 2025-09-25
# Updated: 2025-10-09
#
# Extends core ChEMBL utilities with OpenFDA black-box text retrieval.
# Adds black_box_text and black_box_url fields for frontend integration.
#

# --- Standard library imports ---
import sqlite3
import json
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from typing import List, Dict

# --- Package/Repo Roots ---
PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]

# --- Imports ---
from src.config import DEBUG
from src.tools.chembl import CHEMBL_DB_PATH, UPDATE_MSG
from src.tools.chembl.openfda import fetch_black_box_text_by_name


MAX_BB_TEXT = 800  # Max FDA black box warning text

# ------------------------------
# Logging
# ------------------------------
def log(msg: str) -> None:
    """Log ChEMBL messages with timestamp."""
    if DEBUG:
        print(f"[ChEMBL] {msg} @ {datetime.now()}")


# ------------------------------
# DB connection
# ------------------------------
def get_connection() -> sqlite3.Connection:
    """Return SQLite connection to the ChEMBL database."""
    if not CHEMBL_DB_PATH.exists():
        raise RuntimeError(f"Database missing: {CHEMBL_DB_PATH}. {UPDATE_MSG}")
    conn = sqlite3.connect(CHEMBL_DB_PATH, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    return conn


# ------------------------------
# Unit normalization
# ------------------------------
def _to_nM(value: float, unit: str) -> float:
    """Convert assay values into nanomolar units for consistency."""
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
def fetch_target_info(conn: sqlite3.Connection, uniprot_id: str, limit: int = 25) -> Dict[str, List]:
    """
    Query ChEMBL for a UniProt target: approved drugs, clinical/preclinical trials,
    and assay bioactivities (normalized to nM). Adds OpenFDA black-box text and URL.
    """
    data = {"approved_drugs": [], "clinical_trials": [], "bioactivity": []}
    cur = conn.cursor()

    # --- Resolve UniProt ID -> target ID ---
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

    # --- FDA-approved drugs ---
    cur.execute("""
        SELECT DISTINCT
            m.chembl_id,
            m.pref_name,
            m.molecule_type,
            cs.standard_inchi_key,
            cs.canonical_smiles,
            m.first_approval,
            m.withdrawn_flag,
            m.black_box_warning,
            di.efo_term,
            di.mesh_heading,
            dm.mechanism_of_action
        FROM drug_mechanism dm
        JOIN molecule_dictionary m ON dm.molregno = m.molregno
        LEFT JOIN compound_structures cs ON m.molregno = cs.molregno
        LEFT JOIN drug_indication di ON di.molregno = m.molregno
        WHERE dm.tid = ? AND m.max_phase = 4
        ORDER BY
            CASE
                WHEN m.withdrawn_flag = 1 THEN 0
                WHEN m.black_box_warning = 1 THEN 1
                WHEN m.first_approval IS NOT NULL THEN 2
                ELSE 3
            END,
            COALESCE(m.first_approval, 0) DESC,
            m.chembl_id
    """, (tid,))

    openfda_cache: Dict[str, tuple | None] = {}

    for r in cur.fetchall():
        indication = r["efo_term"] or r["mesh_heading"] or ""
        mechanism = r["mechanism_of_action"] or ""
        withdrawn = "Yes" if r["withdrawn_flag"] == 1 else ""
        black_box = "Yes" if r["black_box_warning"] == 1 else ""

        pref_name = r["pref_name"] or ""
        bb_text = None
        bb_url = None

        if black_box:
            if pref_name in openfda_cache:
                bb_text, bb_url = openfda_cache[pref_name]  # type: ignore[misc]
            else:
                res = fetch_black_box_text_by_name(pref_name)
                if res:
                    txt, link = res
                    txt = (txt or "").strip().replace("\n", " ")
                    if len(txt) > MAX_BB_TEXT:
                        txt = txt[:MAX_BB_TEXT].rstrip() + "... [truncated, see FDA label]"
                    bb_text, bb_url = txt, link
                openfda_cache[pref_name] = (bb_text, bb_url)

        data["approved_drugs"].append({
            "chembl_id": r["chembl_id"],
            "pref_name": pref_name,  # <-- included for frontend + summaries
            "molecule_type": r["molecule_type"],
            "inchi_key": r["standard_inchi_key"] or "",
            "smiles": r["canonical_smiles"] or "",
            "first_approval": r["first_approval"] or "",
            "indication": indication,
            "mechanism": mechanism,
            "withdrawn": withdrawn,
            "black_box": black_box,
            "black_box_text": bb_text,
            "black_box_url": bb_url,
        })

    # --- Clinical trials (include pref_name) ---
    cur.execute("""
        SELECT DISTINCT
            m.chembl_id,
            m.max_phase,
            m.molecule_type,
            cs.standard_inchi_key,
            cs.canonical_smiles,
            m.pref_name
        FROM drug_mechanism dm
        JOIN molecule_dictionary m ON dm.molregno = m.molregno
        LEFT JOIN compound_structures cs ON m.molregno = cs.molregno
        WHERE dm.tid = ? AND m.max_phase BETWEEN 1 AND 3
        ORDER BY m.max_phase DESC, m.chembl_id
    """, (tid,))
    data["clinical_trials"] = [
        (
            r["chembl_id"],
            r["max_phase"],
            r["molecule_type"],
            r["standard_inchi_key"] or "",
            r["canonical_smiles"] or "",
            r["pref_name"] or "",  # <-- added
        )
        for r in cur.fetchall()
    ]

    # --- Preclinical (include pref_name) ---
    cur.execute("""
        SELECT DISTINCT
            m.chembl_id,
            m.molecule_type,
            cs.standard_inchi_key,
            cs.canonical_smiles,
            m.pref_name
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
            r["pref_name"] or "",  # <-- added
        )
        for r in cur.fetchall()
    ]
    if preclinical:
        # Extend as phase 0 and preserve pref_name as last element
        data["clinical_trials"].extend([(d, 0, t, k, s, n) for d, t, k, s, n in preclinical])

    # --- Bioactivity assays (include pref_name) ---
    cur.execute("""
        SELECT
            m.chembl_id,
            m.molecule_type,
            cs.standard_inchi_key,
            cs.canonical_smiles,
            a.standard_type,
            a.standard_value,
            a.standard_units,
            a.pchembl_value,
            m.pref_name
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
        lambda: {
            "values": [],
            "pchembl": [],
            "mtype": None,
            "inchikey": "",
            "smiles": "",
            "pref_name": "",
        }
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
        if r["pref_name"] and not assay_data[key]["pref_name"]:
            assay_data[key]["pref_name"] = r["pref_name"]
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
                vdict["pref_name"],  # <-- added as last element
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
    """Create compact, human-readable druggability summary.
    Groups all indications for the same drug and prints the black-box warning once.
    Includes pref_name after the ChEMBL ID when available (approved, clinical/preclinical, bioactivity).
    """
    lines = [f"Gene: {gene}", f"  UniProt ID: {uniprot_id}"]

    if target_data["approved_drugs"]:
        lines.append("  FDA-approved drugs:")

        # Group all approved drugs by chembl_id
        grouped = {}
        for drug in target_data["approved_drugs"]:
            cid = drug["chembl_id"]
            if cid not in grouped:
                grouped[cid] = {
                    "chembl_id": cid,
                    "pref_name": drug.get("pref_name") or "",
                    "molecule_type": drug["molecule_type"],
                    "first_approval": drug.get("first_approval"),
                    "mechanism": drug.get("mechanism"),
                    "withdrawn": drug.get("withdrawn"),
                    "black_box": drug.get("black_box"),
                    "black_box_text": drug.get("black_box_text"),
                    "black_box_url": drug.get("black_box_url"),
                    "indications": set(),
                }
            else:
                if not grouped[cid]["pref_name"] and drug.get("pref_name"):
                    grouped[cid]["pref_name"] = drug["pref_name"]

            if drug.get("indication"):
                grouped[cid]["indications"].add(drug["indication"])

        # Format grouped entries
        for cid, info in grouped.items():
            # base: CHEMBL ID (+ molecule type if non-small molecule)
            entry_head = _fmt_entry(cid, info['molecule_type'])
            # add pref name right after the ID if present
            if info.get("pref_name"):
                entry_head += f" ({info['pref_name']})"

            entry = f"    - {entry_head}"
            details = []
            if info["first_approval"]:
                details.append(f"approved {info['first_approval']}")
            if info["withdrawn"]:
                details.append("withdrawn")
            if info["indications"]:
                ind_list = ", ".join(sorted(info["indications"]))
                details.append(f"for {ind_list}")
            if info["mechanism"]:
                details.append(f"MoA: {info['mechanism']}")
            if info["black_box"]:
                details.append("black box warning")
            if details:
                entry += " – " + ", ".join(details)
            lines.append(entry)

            # Single warning section per unique drug (indented for clarity)
            if info.get("black_box_text"):
                warning_text = info["black_box_text"].strip().replace("\n", " ").strip()
                lines.append(f"          {warning_text}")
            # if info.get("black_box_url"):
            #     lines.append(f"          FDA Label: {info['black_box_url']}")

    if target_data["clinical_trials"]:
        lines.append("  Clinical/preclinical trials:")
        for row in target_data["clinical_trials"]:
            # Supports old (5 fields) and new (6 fields with pref_name) shapes
            chembl_id, phase, mtype, *rest = row
            inchikey = rest[0] if len(rest) > 0 else ""
            smiles = rest[1] if len(rest) > 1 else ""
            pref_name = rest[2] if len(rest) > 2 else ""

            phase_str = "preclinical" if phase == 0 else f"phase {phase}"
            entry_head = _fmt_entry(chembl_id, mtype)
            if pref_name:
                entry_head += f" ({pref_name})"
            lines.append(f"    - {entry_head} ({phase_str})")

    if target_data["bioactivity"]:
        lines.append("  Bioactivity evidence (normalized to nM):")
        for row in target_data["bioactivity"]:
            # Supports old (5 fields) and new (6 fields with pref_name) shapes
            chembl_id, mtype, evidence_line, *rest = row
            inchikey = rest[0] if len(rest) > 0 else ""
            smiles = rest[1] if len(rest) > 1 else ""
            pref_name = rest[2] if len(rest) > 2 else ""

            entry_head = _fmt_entry(chembl_id, mtype)
            if pref_name:
                entry_head += f" ({pref_name})"
            lines.append(f"    - {entry_head} {evidence_line}")

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
            "- Includes drug indication, mechanism, approval, withdrawal flag, and black-box info (if available).\n"
            "- Black-box text and FDA label links are retrieved from the openFDA API.\n"
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
    Returns empty HTML if all genes have empty evidence sets.
    """
    has_data = False
    for gene, data in (chembl_data_all or {}).items():
        if not isinstance(data, dict):
            continue
        if (
            (data.get("approved_drugs") and len(data["approved_drugs"]) > 0)
            or (data.get("clinical_trials") and len(data["clinical_trials"]) > 0)
            or (data.get("bioactivity") and len(data["bioactivity"]) > 0)
        ):
            has_data = True
            break

    if not has_data:
        return ""  # disables card automatically

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