# src/tools/chembl/utils.py
# ===========================================================
# ChEMBL Drug-Target Evidence Utilities with OpenFDA Integration
# ===========================================================
#
# Author: Dmitri Kosenkov
# Created: 2025-09-25
# Updated: 2025-12-01
#
# Extends core ChEMBL utilities with OpenFDA black-box text retrieval.
# Adds black_box_text and black_box_url fields for frontend integration.
#
# Also provides drug-centric helpers for:
#   - Normalizing a drug identifier to a canonical ChEMBL record
#   - Retrieving targets for a given ChEMBL compound ID
#

# --- Standard library imports ---
import sqlite3
import json
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from typing import List, Dict, Any, Optional, Tuple

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
    """Log ChEMBL messages with timestamp (DEBUG-gated)."""
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
def _to_nM(value: float, unit: str) -> Optional[float]:
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
# Drug-centric helpers (new)
# ------------------------------
def _looks_like_chembl_id(value: str) -> bool:
    """Heuristic check for ChEMBL compound IDs like 'CHEMBL1234'."""
    v = (value or "").strip().upper()
    return v.startswith("CHEMBL") and len(v) > 6


def _maybe_extract_cas(value: str) -> Optional[str]:
    """
    Very simple CAS pattern detector: NNNNN-NN-N.
    If the string matches this pattern, return it; otherwise None.
    """
    s = (value or "").strip()
    parts = s.split("-")
    if len(parts) != 3:
        return None
    if not all(p.isdigit() for p in parts):
        return None
    # basic length sanity; CAS left block is usually 2-7 digits
    if not (2 <= len(parts[0]) <= 7 and len(parts[1]) == 2 and len(parts[2]) == 1):
        return None
    return s


def normalize_chembl_drug_identifier(name_or_id: str) -> Optional[Dict[str, Any]]:
    """
    Resolve a free-text or ID-like input (ChEMBL ID, preferred name,
    synonym, CAS-like string) to a canonical ChEMBL drug record.

    Returns a dict with:
        {
          "chembl_id": str,
          "pref_name": str,
          "molecule_type": str,
          "inchi_key": str,
          "smiles": str,
          "cas_number": Optional[str],
          "synonyms": List[str],
        }

    Returns None if resolution fails.
    """
    if not name_or_id or not str(name_or_id).strip():
        return None

    value = str(name_or_id).strip()
    conn = get_connection()
    try:
        cur = conn.cursor()

        row: Optional[sqlite3.Row] = None
        cas_candidate: Optional[str] = None

        # --------------------------------------------------
        # 1) Direct ChEMBL ID (e.g., CHEMBL1234)
        # --------------------------------------------------
        if _looks_like_chembl_id(value):
            cur.execute(
                """
                SELECT m.molregno,
                       m.chembl_id,
                       m.pref_name,
                       m.molecule_type,
                       cs.standard_inchi_key,
                       cs.canonical_smiles
                FROM molecule_dictionary m
                LEFT JOIN compound_structures cs
                  ON m.molregno = cs.molregno
                WHERE UPPER(m.chembl_id) = ?
                LIMIT 1
                """,
                (value.upper(),),
            )
            row = cur.fetchone()

        # --------------------------------------------------
        # 2) Exact preferred name match (case-insensitive)
        # --------------------------------------------------
        if row is None:
            cur.execute(
                """
                SELECT m.molregno,
                       m.chembl_id,
                       m.pref_name,
                       m.molecule_type,
                       cs.standard_inchi_key,
                       cs.canonical_smiles
                FROM molecule_dictionary m
                LEFT JOIN compound_structures cs
                  ON m.molregno = cs.molregno
                WHERE LOWER(m.pref_name) = ?
                LIMIT 1
                """,
                (value.lower(),),
            )
            row = cur.fetchone()

        # --------------------------------------------------
        # 3) Synonym match (molecule_synonyms), capturing CAS
        # --------------------------------------------------
        if row is None:
            cur.execute(
                """
                SELECT m.molregno,
                       m.chembl_id,
                       m.pref_name,
                       m.molecule_type,
                       cs.standard_inchi_key,
                       cs.canonical_smiles,
                       ms.synonyms
                FROM molecule_dictionary m
                JOIN molecule_synonyms ms
                  ON m.molregno = ms.molregno
                LEFT JOIN compound_structures cs
                  ON m.molregno = cs.molregno
                WHERE LOWER(ms.synonyms) = ?
                LIMIT 1
                """,
                (value.lower(),),
            )
            row = cur.fetchone()
            if row:
                cas_candidate = _maybe_extract_cas(row["synonyms"] or "")

        if row is None:
            return None

        molregno = row["molregno"]
        chembl_id = row["chembl_id"]
        pref_name = row["pref_name"] or ""
        molecule_type = row["molecule_type"] or ""
        inchi_key = row["standard_inchi_key"] or ""
        smiles = row["canonical_smiles"] or ""

        # --------------------------------------------------
        # All synonyms for this molregno
        # --------------------------------------------------
        cur.execute(
            """
            SELECT synonyms
            FROM molecule_synonyms
            WHERE molregno = ?
            """,
            (molregno,),
        )
        syns_raw = [r["synonyms"] for r in cur.fetchall() if r["synonyms"]]

        # Simple dedup (case-insensitive), preserving first occurrence
        seen: set[str] = set()
        synonyms: List[str] = []
        for s in syns_raw:
            s_norm = (s or "").strip()
            if not s_norm:
                continue
            key = s_norm.lower()
            if key in seen:
                continue
            seen.add(key)
            synonyms.append(s_norm)

        # Detect CAS from synonyms if not already captured
        cas_number = cas_candidate
        if cas_number is None:
            for s in synonyms:
                c = _maybe_extract_cas(s)
                if c:
                    cas_number = c
                    break

        return {
            "chembl_id": chembl_id,
            "pref_name": pref_name,
            "molecule_type": molecule_type,
            "inchi_key": inchi_key,
            "smiles": smiles,
            "cas_number": cas_number,
            "synonyms": synonyms,
        }
    finally:
        conn.close()

def get_chembl_drug_targets(
    chembl_id: str,
    only_human: bool = True,
) -> List[Dict[str, Any]]:
    """
    Given a ChEMBL compound ID, return its known targets from ChEMBL
    drug_mechanism / target_dictionary / component_sequences.

    Each record includes (when available):
        {
          "chembl_id": str,          # compound CHEMBL ID
          "molecule_type": str,
          "target_chembl_id": int,   # TID
          "target_name": str,        # target_dictionary.pref_name
          "target_type": str,        # target_dictionary.target_type
          "organism": str,           # component_sequences.organism or target_dictionary.organism
          "uniprot_id": str,         # component_sequences.accession
          "gene_symbol": str,        # component_synonyms.component_synonym (GENE_SYMBOL)
          "mechanism_of_action": str,
          "action_type": str,        # agonist, antagonist, inhibitor, etc.
        }
    """
    cid = (chembl_id or "").strip()
    if not cid:
        return []

    conn = get_connection()
    try:
        cur = conn.cursor()

        # ----------------------------------------------------------
        # 1) Resolve CHEMBL ID -> molregno + active_molregno
        #    and use the UNION of both when querying mechanisms.
        # ----------------------------------------------------------
        cur.execute(
            """
            SELECT
                md.molregno,
                mh.active_molregno
            FROM molecule_dictionary md
            LEFT JOIN molecule_hierarchy mh
              ON md.molregno = mh.molregno
            WHERE UPPER(md.chembl_id) = ?
            """,
            (cid.upper(),),
        )
        rows = cur.fetchall()
        if not rows:
            if DEBUG:
                log(f"[get_chembl_drug_targets] No molecule_dictionary row for {cid}")
            return []

        molregnos: set[int] = set()
        for r in rows:
            if r["molregno"] is not None:
                molregnos.add(r["molregno"])
            if r["active_molregno"] is not None:
                molregnos.add(r["active_molregno"])

        if not molregnos:
            if DEBUG:
                log(f"[get_chembl_drug_targets] No molregno/active_molregno for {cid}")
            return []

        if DEBUG:
            log(
                f"[get_chembl_drug_targets] {cid} -> molregno set "
                f"{sorted(molregnos)}"
            )

        def _run_main_query(molregno_ids: List[int]) -> List[sqlite3.Row]:
            if not molregno_ids:
                return []
            placeholders = ",".join(["?"] * len(molregno_ids))
            cur.execute(
                f"""
                SELECT
                    m.chembl_id,
                    m.molecule_type,
                    dm.tid AS target_chembl_id,
                    td.pref_name AS target_name,
                    td.target_type,
                    td.organism AS target_organism,
                    cs.accession AS uniprot_id,
                    cs.organism AS component_organism,
                    csy.component_synonym AS gene_symbol,
                    dm.mechanism_of_action,
                    dm.action_type
                FROM drug_mechanism dm
                JOIN molecule_dictionary m
                  ON dm.molregno = m.molregno
                LEFT JOIN target_dictionary td
                  ON dm.tid = td.tid
                LEFT JOIN target_components tc
                  ON td.tid = tc.tid
                LEFT JOIN component_sequences cs
                  ON tc.component_id = cs.component_id
                LEFT JOIN component_synonyms csy
                  ON tc.component_id = csy.component_id
                 AND csy.syn_type = 'GENE_SYMBOL'
                WHERE dm.molregno IN ({placeholders})
                """,
                molregno_ids,
            )
            return cur.fetchall()

        # First attempt: mechanisms on any of the molregnos/active_molregnos
        rows = _run_main_query(sorted(molregnos))

        # ----------------------------------------------------------
        # 2) Fallback: sometimes mechanisms are linked in odd ways.
        #    As a safety net, also try a direct join on chembl_id.
        # ----------------------------------------------------------
        if not rows:
            if DEBUG:
                log(
                    f"[get_chembl_drug_targets] No mechanisms via molregno set "
                    f"for {cid}, trying fallback join on chembl_id"
                )
            cur.execute(
                """
                SELECT
                    m.chembl_id,
                    m.molecule_type,
                    dm.tid AS target_chembl_id,
                    td.pref_name AS target_name,
                    td.target_type,
                    td.organism AS target_organism,
                    cs.accession AS uniprot_id,
                    cs.organism AS component_organism,
                    csy.component_synonym AS gene_symbol,
                    dm.mechanism_of_action,
                    dm.action_type
                FROM drug_mechanism dm
                JOIN molecule_dictionary m
                  ON dm.molregno = m.molregno
                LEFT JOIN target_dictionary td
                  ON dm.tid = td.tid
                LEFT JOIN target_components tc
                  ON td.tid = tc.tid
                LEFT JOIN component_sequences cs
                  ON tc.component_id = cs.component_id
                LEFT JOIN component_synonyms csy
                  ON tc.component_id = csy.component_id
                 AND csy.syn_type = 'GENE_SYMBOL'
                WHERE UPPER(m.chembl_id) = ?
                """,
                (cid.upper(),),
            )
            rows = cur.fetchall()

        if DEBUG:
            log(
                f"[get_chembl_drug_targets] Raw rows for {cid} from drug_mechanism: "
                f"{len(rows)}"
            )

        results: List[Dict[str, Any]] = []
        seen: set[tuple] = set()

        for r in rows:
            target_chembl_id = r["target_chembl_id"]
            target_name = r["target_name"] or ""
            target_type = r["target_type"] or ""

            # Prefer component organism when present; fallback to target organism
            organism = (r["component_organism"] or r["target_organism"] or "").strip()
            uniprot_id = (r["uniprot_id"] or "").strip()
            gene_symbol = (r["gene_symbol"] or "").strip()
            moa = (r["mechanism_of_action"] or "").strip()
            action_type = (r["action_type"] or "").strip()
            molecule_type = (r["molecule_type"] or "").strip()

            # Optional organism filter (focus on human targets if desired)
            if only_human:
                org_low = organism.lower()
                if org_low and ("homo sapiens" not in org_low and "human" not in org_low):
                    continue

            # Deduplicate by (tid, uniprot, gene_symbol)
            key = (target_chembl_id, uniprot_id, gene_symbol)
            if key in seen:
                continue
            seen.add(key)

            results.append(
                {
                    "chembl_id": r["chembl_id"],
                    "molecule_type": molecule_type,
                    "target_chembl_id": target_chembl_id,
                    "target_name": target_name,
                    "target_type": target_type,
                    "organism": organism,
                    "uniprot_id": uniprot_id,
                    "gene_symbol": gene_symbol,
                    "mechanism_of_action": moa,
                    "action_type": action_type,
                }
            )

        if DEBUG:
            log(
                f"[get_chembl_drug_targets] Resolved {len(results)} targets for {cid} "
                f"(only_human={only_human})"
            )
        return results
    finally:
        conn.close()

# ------------------------------
# Queries (existing, target-centric)
# ------------------------------
def fetch_target_info(conn: sqlite3.Connection, uniprot_id: str, limit: int = 25) -> Dict[str, List]:
    """
    Query ChEMBL for a UniProt target: approved drugs, clinical/preclinical trials,
    and assay bioactivities (normalized to nM). Adds OpenFDA black-box text and URL.
    """
    data: Dict[str, List] = {"approved_drugs": [], "clinical_trials": [], "bioactivity": []}
    cur = conn.cursor()

    # --- Resolve UniProt ID -> target ID ---
    cur.execute(
        """
        SELECT td.tid
        FROM target_dictionary td
        JOIN target_components tc ON td.tid = tc.tid
        JOIN component_sequences cs ON tc.component_id = cs.component_id
        WHERE cs.accession = ?
        LIMIT 1
        """,
        (uniprot_id,),
    )
    row = cur.fetchone()
    if not row:
        return data
    tid = row["tid"]

    # --- FDA-approved drugs ---
    cur.execute(
        """
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
        """,
        (tid,),
    )

    openfda_cache: Dict[str, Tuple[Optional[str], Optional[str]]] = {}

    for r in cur.fetchall():
        indication = r["efo_term"] or r["mesh_heading"] or ""
        mechanism = r["mechanism_of_action"] or ""
        withdrawn = "Yes" if r["withdrawn_flag"] == 1 else ""
        black_box = "Yes" if r["black_box_warning"] == 1 else ""

        pref_name = r["pref_name"] or ""
        bb_text: Optional[str] = None
        bb_url: Optional[str] = None

        if black_box:
            if pref_name in openfda_cache:
                bb_text, bb_url = openfda_cache[pref_name]
            else:
                res = fetch_black_box_text_by_name(pref_name)
                if res:
                    txt, link = res
                    txt = (txt or "").strip().replace("\n", " ")
                    if len(txt) > MAX_BB_TEXT:
                        txt = (
                            txt[:MAX_BB_TEXT].rstrip()
                            + "... [truncated, see FDA label]"
                        )
                    bb_text, bb_url = txt, link
                openfda_cache[pref_name] = (bb_text, bb_url)

        data["approved_drugs"].append(
            {
                "chembl_id": r["chembl_id"],
                "pref_name": pref_name,  # for frontend + summaries
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
            }
        )

    # --- Clinical trials (include pref_name) ---
    cur.execute(
        """
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
        """,
        (tid,),
    )
    data["clinical_trials"] = [
        (
            r["chembl_id"],
            r["max_phase"],
            r["molecule_type"],
            r["standard_inchi_key"] or "",
            r["canonical_smiles"] or "",
            r["pref_name"] or "",  # added
        )
        for r in cur.fetchall()
    ]

    # --- Preclinical (include pref_name) ---
    cur.execute(
        """
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
        """,
        (tid,),
    )
    preclinical = [
        (
            r["chembl_id"],
            r["molecule_type"],
            r["standard_inchi_key"] or "",
            r["canonical_smiles"] or "",
            r["pref_name"] or "",  # added
        )
        for r in cur.fetchall()
    ]
    if preclinical:
        # Extend as phase 0 and preserve pref_name as last element
        data["clinical_trials"].extend(
            [(d, 0, t, k, s, n) for d, t, k, s, n in preclinical]
        )

    # --- Bioactivity assays (include pref_name) ---
    cur.execute(
        """
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
        """,
        (tid,),
    )

    assay_data: Dict[
        Tuple[str, str],
        Dict[str, Any],
    ] = defaultdict(
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
                vdict["pref_name"],  # added as last element
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
    """
    Create compact, human-readable druggability summary.
    Groups all indications for the same drug and prints the black-box warning once.
    Includes pref_name after the ChEMBL ID when available (approved, clinical/preclinical, bioactivity).
    """
    lines = [f"Gene: {gene}", f"  UniProt ID: {uniprot_id}"]

    if target_data["approved_drugs"]:
        lines.append("  FDA-approved drugs:")

        # Group all approved drugs by chembl_id
        grouped: Dict[str, Dict[str, Any]] = {}
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

        # Format grouped entries & count approved, withdrawn, black-box
        counter_approved = 0
        counter_withdrawn = 0
        counter_black_box = 0
        for cid, info in grouped.items():
            entry_head = _fmt_entry(cid, info["molecule_type"])
            if info.get("pref_name"):
                entry_head += f" ({info['pref_name']})"

            entry = f"    - {entry_head}"
            details = []
            if info["first_approval"]:
                details.append(f"approved {info['first_approval']}")
                counter_approved += 1
            if info["withdrawn"]:
                details.append("withdrawn")
                counter_withdrawn += 1
            if info["indications"]:
                ind_list = ", ".join(sorted(info["indications"]))
                details.append(f"for {ind_list}")
            if info["mechanism"]:
                details.append(f"MoA: {info['mechanism']}")
            if info["black_box"]:
                details.append("black box warning")
                counter_black_box += 1
            if details:
                entry += " - " + ", ".join(details)
            

            # Single warning section per unique drug (indented for clarity)
            if info.get("black_box_text"):
                warning_text = info["black_box_text"].strip().replace("\n", " ").strip()
                lines.append(f"          {warning_text}")
            # if info.get("black_box_url"):
            #     lines.append(f"          FDA Label: {info['black_box_url']}")
            lines.append(entry+'\n')
        lines.append(
            f"    Summary: in total in ChEMBL there are {counter_approved} approved (ever) drug entities, "
            f"{counter_withdrawn} withdrawn, "
            f"{counter_black_box} with black-box warnings."
        )
        

    if target_data["clinical_trials"]:
        lines.append("  Clinical/preclinical trials:")
        # counter
        counter_trials ={
            'preclinical':0,
            'phase 1':0,
            'phase 2':0,
            'phase 3':0
        }
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
            if phase == 0:
                counter_trials['preclinical'] += 1
            elif phase == 1:
                counter_trials['phase 1'] += 1
            elif phase == 2:
                counter_trials['phase 2'] += 1
            elif phase == 3:
                counter_trials['phase 3'] += 1
        lines.append(
            f"    Summary: in total in the current snapshot of ChEMBL there are {counter_trials['preclinical']} preclinical, "
            f"{counter_trials['phase 1']} phase 1, "
            f"{counter_trials['phase 2']} phase 2, "
            f"{counter_trials['phase 3']} phase 3 trial drug entities."
        )

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
    genes_str: str,
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
        html_template.replace("{{CSS_INLINE}}", css_inline)
        .replace("{{JS_INLINE}}", js_inline)
        .replace("{{GENE_SYMBOL}}", genes_str)
    )
