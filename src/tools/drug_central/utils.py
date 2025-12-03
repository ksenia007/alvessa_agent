# src/tools/drug_central/utils.py
# ===========================================================
# Drug Central Utilities: Drug Normalization, Targets,
# Indications, Safety, and UI Summaries
# ===========================================================
#
# Author: Dmitri Kosenkov (patterned after chembl/utils.py)
# Created: 2025-11-17
#
# This module exposes a DRY, high-level API over the local
# Drug Central SQLite dump for use by the `drug_central` tool.
#
# NOTE: In this module, struct_id ALWAYS refers to structures.id
# (the stable Drug Central structure identifier).
#
# Only ASCII characters are used in code.
#

import sqlite3
import json
from datetime import datetime
from collections import Counter, defaultdict
from typing import Any, Dict, List, Optional, Tuple, TypedDict

from src.config import DEBUG
from src.tools.drug_central import DRUG_CENTRAL_DB_PATH, UPDATE_MSG


# ------------------------------
# Logging
# ------------------------------
def log(msg: str) -> None:
    """Log DrugCentral messages with timestamp when DEBUG is enabled."""
    if DEBUG:
        print(f"[Drug Central] {msg} @ {datetime.now()}")


# ------------------------------
# DB connection
# ------------------------------
def get_connection() -> sqlite3.Connection:
    """Return SQLite connection to the local Drug Central database."""
    if not DRUG_CENTRAL_DB_PATH.exists():
        raise RuntimeError(f"Database missing: {DRUG_CENTRAL_DB_PATH}. {UPDATE_MSG}")
    conn = sqlite3.connect(DRUG_CENTRAL_DB_PATH, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    return conn


# ------------------------------
# Data model (TypedDicts)
# ------------------------------
class DrugRecord(TypedDict, total=False):
    # struct_id is structures.id (stable DrugCentral drug ID)
    struct_id: int
    name: str
    cas_number: str
    inchikey: str
    smiles: str
    status: str
    synonyms: List[str]
    chembl_id: Optional[str]  # mapped from identifier.id_type = 'ChEMBL_ID'


class StructuralInfo(TypedDict, total=False):
    struct_id: int         # structures.id
    name: str
    cas_number: str
    inchikey: str
    smiles: str
    inchi: str
    mol_weight: float
    clogp: float
    tpsa: float
    lipinski_violations: int
    rotatable_bonds: int
    aromatic_carbons: int
    halogen_count: int
    hydrogen_bond_donors: int
    hydrogen_bond_acceptors: int


class RegulatoryInfo(TypedDict, total=False):
    struct_id: int         # structures.id
    status: str
    approvals: List[Dict[str, Any]]
    orange_book_products: List[Dict[str, Any]]
    exclusivity: List[Dict[str, Any]]
    patents: List[Dict[str, Any]]


class TargetRecord(TypedDict, total=False):
    struct_id: int           # Drug Central struct_id = structures.id
    target_id: int
    target_name: str
    target_class: str
    accession: str
    gene: str
    swissprot: str
    tdl: str
    action_type: str
    is_moa: bool
    act_type: str
    act_value: float
    act_unit: str
    relation: str
    organism: str


class DrugTargetLink(TypedDict, total=False):
    struct_id: int           # Drug Central struct_id = structures.id
    gene: str
    accession: str
    swissprot: str
    target_id: int
    action_type: str
    is_moa: bool
    act_summary: str


class GeneMapping(TypedDict, total=False):
    target_id: int
    name: str
    gene: str
    geneid: int
    accession: str
    swissprot: str
    tdl: str
    target_class: str
    organism: str


class IndicationRecord(TypedDict, total=False):
    struct_id: int           # structures.id
    concept_id: int
    concept_name: str
    relationship_name: str
    umls_cui: str
    snomed_conceptid: int
    snomed_full_name: str
    cui_semantic_type: str


class ContraindicationRecord(TypedDict, total=False):
    struct_id: int           # structures.id
    concept_name: str
    relationship_name: str
    umls_cui: str
    snomed_conceptid: int


class PharmActionRecord(TypedDict, total=False):
    struct_id: int           # structures.id
    type: str
    name: str
    class_code: str
    source: str


class SafetyProfile(TypedDict, total=False):
    struct_id: int           # structures.id
    faers_signals: List[Dict[str, Any]]
    faers_female: List[Dict[str, Any]]
    faers_male: List[Dict[str, Any]]
    faers_pediatric: List[Dict[str, Any]]
    faers_geriatrics: List[Dict[str, Any]]
    ddi_risks: List[Dict[str, Any]]


class CrossRefInfo(TypedDict, total=False):
    struct_id: int           # structures.id
    identifiers: List[Dict[str, str]]
    atc_codes: List[str]
    drug_classes: List[str]


class TargetClassSummary(TypedDict, total=False):
    struct_id: int           # structures.id
    total_targets: int
    by_class: Dict[str, int]
    by_tdl: Dict[str, int]
    moa_targets: int


class HealthStatus(TypedDict, total=False):
    ok: bool
    message: str
    db_path: str
    version: Optional[int]
    tables_checked: Dict[str, int]


class HumanGeneTarget(TypedDict, total=False):
    """Normalized human gene target projection for a drug."""

    symbol: str
    entrez_id: Optional[str]
    uniprot_id: Optional[str]
    ensembl_id: Optional[str]
    source: str

    # Optional scoring / summary fields used for target expansion
    score: float            # heuristic score: higher = more DrugCentral evidence
    n_activities: int       # number of act_table_full rows for this gene
    n_moa: int              # how many of those rows have is_moa = True


# ------------------------------
# Helper: struct_id mapping
# ------------------------------
def _parse_struct_id(value: Any) -> Optional[int]:
    """Try to parse a value as struct_id (structures.id) integer."""
    if value is None:
        return None
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


# ------------------------------
# Core lookup utilities
# ------------------------------
def _fetch_chembl_id(conn: sqlite3.Connection, struct_id: int) -> Optional[str]:
    """
    Return the ChEMBL ID for a given Drug Central struct_id (structures.id),
    if present in the identifier table as id_type = 'ChEMBL_ID'.
    """
    cur = conn.cursor()
    cur.execute(
        """
        SELECT identifier
        FROM identifier
        WHERE struct_id = ?
          AND upper(id_type) = 'CHEMBL_ID'
        ORDER BY id
        LIMIT 1
        """,
        (struct_id,),
    )
    row = cur.fetchone()
    if row:
        ident = row["identifier"]
        if ident:
            return str(ident)
    return None


def _fetch_struct_by_struct_id(conn: sqlite3.Connection, struct_id: int) -> Optional[DrugRecord]:
    """
    Fetch a DrugRecord by structures.id (struct_id).
    """
    cur = conn.cursor()
    cur.execute(
        """
        SELECT
            id AS struct_id,
            name,
            cas_reg_no AS cas_number,
            inchikey,
            smiles,
            status
        FROM structures
        WHERE id = ?
        """,
        (struct_id,),
    )
    row = cur.fetchone()
    if not row:
        return None

    # synonyms: synonyms.id = structures.id
    syn_cur = conn.cursor()
    syn_cur.execute(
        """
        SELECT s.name
        FROM synonyms s
        JOIN structures st ON s.id = st.id
        WHERE st.id = ?
        """,
        (struct_id,),
    )
    syns = [r["name"] for r in syn_cur.fetchall()]

    # ChEMBL ID from identifier table (id_type = 'ChEMBL_ID')
    chembl_id = _fetch_chembl_id(conn, row["struct_id"])

    return DrugRecord(
        struct_id=row["struct_id"],          # structures.id
        name=row["name"] or "",
        cas_number=row["cas_number"] or "",
        inchikey=row["inchikey"] or "",
        smiles=row["smiles"] or "",
        status=row["status"] or "",
        synonyms=syns,
        chembl_id=chembl_id,
    )


def _fetch_struct_by_name_or_synonym(conn: sqlite3.Connection, name: str) -> Optional[DrugRecord]:
    """Resolve a name or synonym to a single struct_id (structures.id) if possible."""
    cur = conn.cursor()
    name_l = name.lower().strip()

    # Try exact match on structures.name
    cur.execute(
        """
        SELECT id
        FROM structures
        WHERE lower(name) = ?
        LIMIT 1
        """,
        (name_l,),
    )
    row = cur.fetchone()
    if row:
        return _fetch_struct_by_struct_id(conn, row["id"])

    # Try synonyms
    cur.execute(
        """
        SELECT DISTINCT st.id
        FROM synonyms s
        JOIN structures st ON s.id = st.id
        WHERE lower(s.name) = ?
        LIMIT 1
        """,
        (name_l,),
    )
    row = cur.fetchone()
    if row:
        return _fetch_struct_by_struct_id(conn, row["id"])

    # Try product names -> active_ingredient.struct_id -> structures.id
    cur.execute(
        """
        SELECT DISTINCT st.id
        FROM product p
        JOIN active_ingredient ai
          ON p.ndc_product_code = ai.ndc_product_code
        JOIN structures st
          ON ai.struct_id = st.id
        WHERE lower(p.product_name) = ? OR lower(p.generic_name) = ?
        LIMIT 1
        """,
        (name_l, name_l),
    )
    row = cur.fetchone()
    if row:
        return _fetch_struct_by_struct_id(conn, row["id"])

    return None


def _fetch_struct_by_cas(conn: sqlite3.Connection, cas_number: str) -> Optional[DrugRecord]:
    cur = conn.cursor()
    cur.execute(
        """
        SELECT id
        FROM structures
        WHERE cas_reg_no = ?
        LIMIT 1
        """,
        (cas_number.strip(),),
    )
    row = cur.fetchone()
    if row:
        return _fetch_struct_by_struct_id(conn, row["id"])
    return None


def _fetch_struct_by_identifier(conn: sqlite3.Connection, identifier: str) -> Optional[DrugRecord]:
    """
    Try to resolve by identifier table (ChEMBL, PubChem CID, RxNorm, etc.).
    Case-insensitive over identifier to support 'chembl521' vs 'CHEMBL521'.
    """
    cur = conn.cursor()
    cur.execute(
        """
        SELECT DISTINCT struct_id
        FROM identifier
        WHERE lower(identifier) = lower(?)
        LIMIT 1
        """,
        (identifier.strip(),),
    )
    row = cur.fetchone()
    if row:
        return _fetch_struct_by_struct_id(conn, row["struct_id"])
    return None


# ------------------------------
# Public: Drug normalization
# ------------------------------
def normalize_drug_identifier(name_or_id: str) -> Optional[DrugRecord]:
    """
    Resolve a free-text or ID-like input (generic name, brand name, synonym,
    struct_id (structures.id), Drug Central identifier, ChEMBL ID if mapped,
    InChIKey, CAS, etc.) to a canonical DrugRecord.

    struct_id in the returned DrugRecord is ALWAYS structures.id.

    Returns None if resolution fails.
    """
    if not name_or_id or not str(name_or_id).strip():
        return None

    value = str(name_or_id).strip()
    conn = get_connection()
    try:
        # 1) Direct struct_id (structures.id)
        sid = _parse_struct_id(value)
        if sid is not None:
            rec = _fetch_struct_by_struct_id(conn, sid)
            if rec:
                return rec

        # 2) CAS-like pattern
        if value.count("-") == 2 and all(part.isdigit() for part in value.split("-")):
            rec = _fetch_struct_by_cas(conn, value)
            if rec:
                return rec

        # 3) Identifier table (ChEMBL, PubChem, RxNorm, etc.)
        rec = _fetch_struct_by_identifier(conn, value)
        if rec:
            return rec

        # 4) Name or synonym or product name
        rec = _fetch_struct_by_name_or_synonym(conn, value)
        if rec:
            return rec

        # 5) Try InChIKey match via structures table
        cur = conn.cursor()
        cur.execute(
            """
            SELECT id
            FROM structures
            WHERE inchikey = ?
            LIMIT 1
            """,
            (value,),
        )
        row = cur.fetchone()
        if row:
            return _fetch_struct_by_struct_id(conn, row["id"])

        # 6) Try SMILES exact match (rare but useful in some contexts)
        cur.execute(
            """
            SELECT id
            FROM structures
            WHERE smiles = ?
            LIMIT 1
            """,
            (value,),
        )
        row = cur.fetchone()
        if row:
            return _fetch_struct_by_struct_id(conn, row["id"])

        return None
    finally:
        conn.close()


def get_drug_by_id(drug_central_id: str) -> Optional[DrugRecord]:
    """
    Return the core Drug Central record for a known Drug Central struct_id
    (no fuzzy logic, direct lookup). Here struct_id = structures.id.
    """
    sid = _parse_struct_id(drug_central_id)
    if sid is None:
        return None
    conn = get_connection()
    try:
        return _fetch_struct_by_struct_id(conn, sid)
    finally:
        conn.close()


# ------------------------------
# Structural information
# ------------------------------
def get_drug_structural_info(drug_central_id: str) -> Optional[StructuralInfo]:
    """Return structural info and basic physchem properties for a struct_id (structures.id)."""
    sid = _parse_struct_id(drug_central_id)
    if sid is None:
        return None
    conn = get_connection()
    try:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT
                id AS struct_id,
                name,
                cas_reg_no AS cas_number,
                inchikey,
                smiles,
                inchi,
                cd_molweight AS mol_weight,
                clogp,
                tpsa,
                lipinski,
                rotb,
                arom_c,
                halogen,
                o_n,
                oh_nh
            FROM structures
            WHERE id = ?
            """,
            (sid,),
        )
        row = cur.fetchone()
        if not row:
            return None

        return StructuralInfo(
            struct_id=row["struct_id"],
            name=row["name"] or "",
            cas_number=row["cas_number"] or "",
            inchikey=row["inchikey"] or "",
            smiles=row["smiles"] or "",
            inchi=row["inchi"] or "",
            mol_weight=row["mol_weight"] or 0.0,
            clogp=row["clogp"] or 0.0,
            tpsa=row["tpsa"] or 0.0,
            lipinski_violations=row["lipinski"] or 0,
            rotatable_bonds=row["rotb"] or 0,
            aromatic_carbons=row["arom_c"] or 0,
            halogen_count=row["halogen"] or 0,
            hydrogen_bond_acceptors=row["o_n"] or 0,
            hydrogen_bond_donors=row["oh_nh"] or 0,
        )
    finally:
        conn.close()


# ------------------------------
# Regulatory information
# ------------------------------
def get_drug_regulatory_status(drug_central_id: str) -> Optional[RegulatoryInfo]:
    """Return regulatory status and metadata for a struct_id (structures.id)."""
    sid = _parse_struct_id(drug_central_id)
    if sid is None:
        return None

    conn = get_connection()
    try:
        cur = conn.cursor()

        # Base status from structures
        cur.execute(
            """
            SELECT status
            FROM structures
            WHERE id = ?
            """,
            (sid,),
        )
        row = cur.fetchone()
        status = row["status"] if row else ""

        # Approval table (FDA, EMA, etc.) - approval.struct_id -> structures.id
        cur.execute(
            """
            SELECT id,
                   approval,
                   type,
                   applicant,
                   orphan
            FROM approval
            WHERE struct_id = ?
            ORDER BY approval
            """,
            (sid,),
        )
        approvals = []
        for r in cur.fetchall():
            approvals.append(
                {
                    "id": r["id"],
                    "approval": r["approval"] or "",
                    "type": r["type"] or "",
                    "applicant": r["applicant"] or "",
                    "orphan": bool(r["orphan"]),
                }
            )

        # Orange Book products (ob_product) via struct2obprod.struct_id -> structures.id
        cur.execute(
            """
            SELECT op.id,
                   op.ingredient,
                   op.trade_name,
                   op.applicant_full_name,
                   op.strength,
                   op.dose_form,
                   op.route,
                   op.appl_type,
                   op.appl_no,
                   op.product_no,
                   op.te_code,
                   op.approval_date,
                   op.type,
                   op.rld
            FROM struct2obprod s2p
            JOIN ob_product op ON s2p.prod_id = op.id
            WHERE s2p.struct_id = ?
            """,
            (sid,),
        )
        ob_products = []
        for r in cur.fetchall():
            ob_products.append(
                {
                    "id": r["id"],
                    "ingredient": r["ingredient"],
                    "trade_name": r["trade_name"],
                    "applicant": r["applicant_full_name"],
                    "strength": r["strength"],
                    "dose_form": r["dose_form"],
                    "route": r["route"],
                    "appl_type": r["appl_type"],
                    "appl_no": r["appl_no"],
                    "product_no": r["product_no"],
                    "te_code": r["te_code"],
                    "approval_date": r["approval_date"],
                    "type": r["type"],
                    "rld": bool(r["rld"]),
                }
            )

        # Exclusivity
        cur.execute(
            """
            SELECT e.id,
                   e.appl_type,
                   e.appl_no,
                   e.product_no,
                   e.exclusivity_code,
                   e.exclusivity_date
            FROM struct2obprod s2p
            JOIN ob_product op
              ON s2p.prod_id = op.id
            JOIN ob_exclusivity e
              ON e.appl_type = op.appl_type
             AND e.appl_no = op.appl_no
             AND e.product_no = op.product_no
            WHERE s2p.struct_id = ?
            """,
            (sid,),
        )
        exclusivity = []
        for r in cur.fetchall():
            exclusivity.append(
                {
                    "id": r["id"],
                    "appl_type": r["appl_type"],
                    "appl_no": r["appl_no"],
                    "product_no": r["product_no"],
                    "exclusivity_code": r["exclusivity_code"],
                    "exclusivity_date": r["exclusivity_date"],
                }
            )

        # Patents
        cur.execute(
            """
            SELECT p.id,
                   p.appl_type,
                   p.appl_no,
                   p.product_no,
                   p.patent_no,
                   p.patent_expire_date,
                   p.drug_substance_flag,
                   p.drug_product_flag,
                   p.patent_use_code,
                   p.delist_flag
            FROM struct2obprod s2p
            JOIN ob_product op
              ON s2p.prod_id = op.id
            JOIN ob_patent p
              ON p.appl_type = op.appl_type
             AND p.appl_no = op.appl_no
             AND p.product_no = op.product_no
            WHERE s2p.struct_id = ?
            """,
            (sid,),
        )
        patents = []
        for r in cur.fetchall():
            patents.append(
                {
                    "id": r["id"],
                    "appl_type": r["appl_type"],
                    "appl_no": r["appl_no"],
                    "product_no": r["product_no"],
                    "patent_no": r["patent_no"],
                    "patent_expire_date": r["patent_expire_date"],
                    "drug_substance_flag": r["drug_substance_flag"],
                    "drug_product_flag": r["drug_product_flag"],
                    "patent_use_code": r["patent_use_code"],
                    "delist_flag": r["delist_flag"],
                }
            )

        return RegulatoryInfo(
            struct_id=sid,
            status=status or "",
            approvals=approvals,
            orange_book_products=ob_products,
            exclusivity=exclusivity,
            patents=patents,
        )
    finally:
        conn.close()


# ------------------------------
# Targets and gene-centric views
# ------------------------------
def get_drug_targets(
    drug_central_id: str,
    include_off_target: bool = True,
) -> List[TargetRecord]:
    """
    Return drug target associations for a drug, including MoA and
    optionally secondary/off-target interactions derived from act_table_full.

    The input `drug_central_id` is treated as Drug Central struct_id = structures.id.
    act_table_full.struct_id also points to structures.id.
    """
    sid = _parse_struct_id(drug_central_id)
    if sid is None:
        return []
    conn = get_connection()
    try:
        cur = conn.cursor()
        params: List[Any] = [sid]
        sql = """
            SELECT
                a.struct_id AS struct_id,
                a.target_id,
                a.target_name,
                a.target_class,
                a.accession,
                a.gene,
                a.swissprot,
                a.tdl,
                a.action_type,
                a.moa,
                a.act_type,
                a.act_value,
                a.act_unit,
                a.relation,
                a.organism
            FROM act_table_full a
            WHERE a.struct_id = ?
        """
        if not include_off_target:
            sql += " AND a.moa = 1"
        cur.execute(sql, tuple(params))

        results: List[TargetRecord] = []
        for r in cur.fetchall():
            results.append(
                TargetRecord(
                    struct_id=r["struct_id"],            # structures.id
                    target_id=r["target_id"],
                    target_name=r["target_name"] or "",
                    target_class=r["target_class"] or "",
                    accession=r["accession"] or "",
                    gene=r["gene"] or "",
                    swissprot=r["swissprot"] or "",
                    tdl=r["tdl"] or "",
                    action_type=r["action_type"] or "",
                    is_moa=bool(r["moa"]),
                    act_type=r["act_type"] or "",
                    act_value=r["act_value"] or 0.0,
                    act_unit=r["act_unit"] or "",
                    relation=r["relation"] or "",
                    organism=r["organism"] or "",
                )
            )
        return results
    finally:
        conn.close()


def get_drugs_for_gene(
    gene_symbol: str,
    role_filter: Optional[str] = None,
) -> List[DrugTargetLink]:
    """
    Given a gene symbol or UniProt accession, return drugs acting on that
    gene's product. Role filter can be "MOA", "off_target", or "binding_only".

    Returns DrugTargetLink records with struct_id = Drug Central struct_id
    (structures.id), so they can be used directly as tool/viewer IDs.
    """
    if not gene_symbol:
        return []

    g = gene_symbol.strip().upper()
    conn = get_connection()
    try:
        cur = conn.cursor()
        sql = """
            SELECT
                a.struct_id AS struct_id,
                a.target_id,
                a.gene,
                a.accession,
                a.swissprot,
                a.action_type,
                a.moa,
                a.act_type,
                a.act_value,
                a.act_unit,
                a.relation
            FROM act_table_full a
            WHERE upper(a.gene) = ?
               OR upper(a.accession) = ?
               OR upper(a.swissprot) = ?
        """
        params: List[Any] = [g, g, g]
        cur.execute(sql, tuple(params))

        links: List[DrugTargetLink] = []
        for r in cur.fetchall():
            is_moa = bool(r["moa"])
            if role_filter == "MOA" and not is_moa:
                continue
            if role_filter == "off_target" and is_moa:
                continue
            # binding_only is heuristic; no extra filter here.

            act_value = r["act_value"]
            act_unit = r["act_unit"] or ""
            act_type = r["act_type"] or ""
            relation = r["relation"] or ""
            if act_value:
                summary = f"{act_type} {relation} {act_value} {act_unit}".strip()
            else:
                summary = act_type or ""

            links.append(
                DrugTargetLink(
                    struct_id=r["struct_id"],          # structures.id
                    gene=r["gene"] or "",
                    accession=r["accession"] or "",
                    swissprot=r["swissprot"] or "",
                    target_id=r["target_id"],
                    action_type=r["action_type"] or "",
                    is_moa=is_moa,
                    act_summary=summary,
                )
            )
        return links
    finally:
        conn.close()


def map_target_identifier_to_gene(target_id: str) -> Optional[GeneMapping]:
    """
    Map a Drug Central target identifier (numeric target_id) to gene-centric
    identifiers using target_dictionary and target_component.
    """
    if not target_id:
        return None
    try:
        tid = int(target_id)
    except ValueError:
        return None

    conn = get_connection()
    try:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT td.id AS target_id,
                   td.name,
                   td.target_class,
                   td.tdl,
                   tc.accession,
                   tc.swissprot,
                   tc.organism,
                   tc.gene,
                   tc.geneid
            FROM target_dictionary td
            LEFT JOIN td2tc ON td.id = td2tc.target_id
            LEFT JOIN target_component tc ON td2tc.component_id = tc.id
            WHERE td.id = ?
            """,
            (tid,),
        )
        row = cur.fetchone()
        if not row:
            return None
        return GeneMapping(
            target_id=row["target_id"],
            name=row["name"] or "",
            gene=row["gene"] or "",
            geneid=row["geneid"] or 0,
            accession=row["accession"] or "",
            swissprot=row["swissprot"] or "",
            tdl=row["tdl"] or "",
            target_class=row["target_class"] or "",
            organism=row["organism"] or "",
        )
    finally:
        conn.close()


def project_targets_to_human_genes(
    targets: List[TargetRecord],
) -> List[HumanGeneTarget]:
    """
    Project raw Drug Central TargetRecord list onto a normalized, human-only
    gene-level view.

    This mirrors the logic used in the entity-extraction / drug_extraction
    pipeline for Drug Central targets, so that:
      - Genes added via drug->gene expansion in entity_extraction.py
      - and gene targets shown in DrugCentral text summaries
    are consistent.

    Rules:
      - Use target.gene, swissprot, accession, and map_target_identifier_to_gene
        to fill missing gene symbol and accession.
      - Keep only HUMAN entries when swissprot/accession ends with '_HUMAN'
        (e.g., COMT_HUMAN).
      - Deduplicate by (symbol, UniProt), keeping:
          * first non-empty Entrez ID,
          * activity counts (n_activities, n_moa),
          * a simple heuristic score:
              score = n_activities + n_moa
            (MoA-labeled interactions are weighted slightly higher).
    """
    by_key: Dict[Tuple[str, str], HumanGeneTarget] = {}
    gene_map_cache: Dict[int, Optional[GeneMapping]] = {}

    for t in targets or []:
        gene_symbol = (t.get("gene") or "").strip()
        accession = (t.get("swissprot") or t.get("accession") or "").strip()
        target_id_val = t.get("target_id")
        entrez_id: Optional[str] = None

        # Map via target_dictionary / target_component when possible
        if target_id_val is not None:
            try:
                tid = int(target_id_val)
            except (TypeError, ValueError):
                tid = None
            if tid is not None:
                if tid not in gene_map_cache:
                    gene_map_cache[tid] = map_target_identifier_to_gene(str(tid))
                gm = gene_map_cache.get(tid) or None
                if gm:
                    if not gene_symbol:
                        gene_symbol = (gm.get("gene") or "").strip()
                    if not accession:
                        accession = (gm.get("swissprot") or gm.get("accession") or "").strip()
                    geneid_val = gm.get("geneid")
                    if geneid_val:
                        entrez_id = str(geneid_val)

        # Skip completely unmapped rows
        if not gene_symbol and not accession:
            continue

        acc_upper = accession.upper()
        # HUMAN-only filter when accession contains species suffix
        if acc_upper and not acc_upper.endswith("_HUMAN"):
            continue

        symbol_norm = gene_symbol.upper() if gene_symbol else ""
        uniprot_norm = acc_upper if acc_upper else ""
        key = (symbol_norm, uniprot_norm)

        symbol_out = symbol_norm or uniprot_norm or ""

        existing = by_key.get(key)
        if existing is None:
            # Initialize entry with counts and optional IDs
            entry: HumanGeneTarget = HumanGeneTarget(
                symbol=symbol_out,
                entrez_id=entrez_id,
                uniprot_id=acc_upper or None,
                ensembl_id=None,
                source="drugcentral",
                n_activities=0,
                n_moa=0,
            )
            by_key[key] = entry
            existing = entry
        else:
            # Fill Entrez if missing
            if not existing.get("entrez_id") and entrez_id:
                existing["entrez_id"] = entrez_id
            # Ensure UniProt is set if we gained it later
            if not existing.get("uniprot_id") and acc_upper:
                existing["uniprot_id"] = acc_upper

        # Update counters
        existing["n_activities"] = int(existing.get("n_activities") or 0) + 1
        if t.get("is_moa"):
            existing["n_moa"] = int(existing.get("n_moa") or 0) + 1

    # Finalize scores
    for entry in by_key.values():
        n_act = int(entry.get("n_activities") or 0)
        n_moa = int(entry.get("n_moa") or 0)
        # Simple heuristic: MoA-labeled interactions get extra weight
        entry["score"] = float(n_act + n_moa)

    return list(by_key.values())


def compute_human_gene_target_scores(
    targets: List[TargetRecord],
) -> Dict[Tuple[str, str], float]:
    """
    Backwards-compatible wrapper for older code paths that expect a dict of
    numeric scores and iterate via `.items()`.

    Input:
      - `targets`: list of act_table_full-style TargetRecord rows.

    Output:
      - dict keyed by (symbol, uniprot_id) with a float score per key.

    Internally:
      1) Call project_targets_to_human_genes(targets) to aggregate counts and
         compute per-gene scores.
      2) Build a dict {(symbol, uniprot_id): score} where:
           - symbol and uniprot_id are uppercased,
           - score is taken from the HumanGeneTarget["score"] field,
           - if the same key appears multiple times we sum scores to remain
             compatible with existing aggregation code that does:
                 score_by_key[key] = score_by_key.get(key, 0) + val
    """
    human_targets = project_targets_to_human_genes(targets)

    scores: Dict[Tuple[str, str], float] = {}
    for entry in human_targets:
        symbol = (entry.get("symbol") or "").upper()
        uniprot = (entry.get("uniprot_id") or "").upper()
        key = (symbol, uniprot)
        score = float(entry.get("score") or 0.0)

        # Same semantics as existing aggregation in drug_extraction:
        # accumulate scores for the same (symbol, uniprot) key.
        scores[key] = scores.get(key, 0.0) + score

    return scores


# ------------------------------
# Indications and Contraindications
# ------------------------------
def get_drug_indications(drug_central_id: str) -> List[IndicationRecord]:
    """
    Return indications for a drug using OMOP relationship table.
    relationship_name is typically "indication" or similar.

    omop_relationship.struct_id references structures.id.
    """
    sid = _parse_struct_id(drug_central_id)
    if sid is None:
        return []
    conn = get_connection()
    try:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT
                struct_id,
                concept_id,
                concept_name,
                relationship_name,
                umls_cui,
                snomed_conceptid,
                snomed_full_name,
                cui_semantic_type
            FROM omop_relationship
            WHERE struct_id = ?
              AND lower(relationship_name) LIKE 'indication%%'
            """,
            (sid,),
        )
        res: List[IndicationRecord] = []
        for r in cur.fetchall():
            res.append(
                IndicationRecord(
                    struct_id=r["struct_id"],
                    concept_id=r["concept_id"],
                    concept_name=r["concept_name"] or "",
                    relationship_name=r["relationship_name"] or "",
                    umls_cui=r["umls_cui"] or "",
                    snomed_conceptid=r["snomed_conceptid"] or 0,
                    snomed_full_name=r["snomed_full_name"] or "",
                    cui_semantic_type=r["cui_semantic_type"] or "",
                )
            )
        return res
    finally:
        conn.close()


def get_drugs_for_indication(disease_term_or_id: str) -> List[DrugRecord]:
    """
    Given a disease/phenotype term or ID (concept_id, UMLS CUI, SNOMED),
    return drugs with matching indications.

    All struct_id values here are structures.id.
    """
    if not disease_term_or_id:
        return []
    val = disease_term_or_id.strip()
    conn = get_connection()
    try:
        cur = conn.cursor()
        # First try numeric concept_id
        concept_id = None
        try:
            concept_id = int(val)
        except ValueError:
            pass

        if concept_id is not None:
            cur.execute(
                """
                SELECT DISTINCT struct_id
                FROM omop_relationship
                WHERE concept_id = ?
                  AND lower(relationship_name) LIKE 'indication%%'
                """,
                (concept_id,),
            )
        else:
            cur.execute(
                """
                SELECT DISTINCT struct_id
                FROM omop_relationship
                WHERE (
                    lower(concept_name) LIKE ?
                    OR umls_cui = ?
                    OR CAST(snomed_conceptid AS TEXT) = ?
                )
                  AND lower(relationship_name) LIKE 'indication%%'
                """,
                (f"%{val.lower()}%", val, val),
            )

        struct_ids = [r["struct_id"] for r in cur.fetchall()]
        records: List[DrugRecord] = []
        for sid in struct_ids:
            rec = _fetch_struct_by_struct_id(conn, sid)
            if rec:
                records.append(rec)
        return records
    finally:
        conn.close()


def get_drug_contraindications(drug_central_id: str) -> List[ContraindicationRecord]:
    """
    Return contraindications for a drug using OMOP relationship table.
    We use relationship_name patterns like '%contraindication%'.

    omop_relationship.struct_id references structures.id.
    """
    sid = _parse_struct_id(drug_central_id)
    if sid is None:
        return []
    conn = get_connection()
    try:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT
                struct_id,
                concept_id,
                concept_name,
                relationship_name,
                umls_cui,
                snomed_conceptid
            FROM omop_relationship
            WHERE struct_id = ?
              AND lower(relationship_name) LIKE '%contraindication%'
            """,
            (sid,),
        )
        res: List[ContraindicationRecord] = []
        for r in cur.fetchall():
            res.append(
                ContraindicationRecord(
                    struct_id=r["struct_id"],
                    concept_name=r["concept_name"] or "",
                    relationship_name=r["relationship_name"] or "",
                    umls_cui=r["umls_cui"] or "",
                    snomed_conceptid=r["snomed_conceptid"] or 0,
                )
            )
        return res
    finally:
        conn.close()


# ------------------------------
# Pharmacologic actions / mechanisms
# ------------------------------
def get_drug_pharmacologic_actions(drug_central_id: str) -> List[PharmActionRecord]:
    """
    Return pharmacologic action annotations for a drug using pharma_class.
    type can be MOA, PE, EPC, etc.

    pharma_class.struct_id references structures.id.
    """
    sid = _parse_struct_id(drug_central_id)
    if sid is None:
        return []
    conn = get_connection()
    try:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT
                struct_id,
                type,
                name,
                class_code,
                source
            FROM pharma_class
            WHERE struct_id = ?
            """,
            (sid,),
        )
        res: List[PharmActionRecord] = []
        for r in cur.fetchall():
            res.append(
                PharmActionRecord(
                    struct_id=r["struct_id"],
                    type=r["type"] or "",
                    name=r["name"] or "",
                    class_code=r["class_code"] or "",
                    source=r["source"] or "",
                )
            )
        return res
    finally:
        conn.close()


def get_drugs_by_mechanism(mechanism_term_or_id: str) -> List[DrugRecord]:
    """
    Return all drugs annotated with a given pharmacologic action term
    (name or class_code). struct_id values are structures.id.
    """
    if not mechanism_term_or_id:
        return []
    val = mechanism_term_or_id.strip()
    conn = get_connection()
    try:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT DISTINCT struct_id
            FROM pharma_class
            WHERE lower(name) LIKE ?
               OR class_code = ?
            """,
            (f"%{val.lower()}%", val),
        )
        struct_ids = [r["struct_id"] for r in cur.fetchall()]
        res: List[DrugRecord] = []
        for sid in struct_ids:
            rec = _fetch_struct_by_struct_id(conn, sid)
            if rec:
                res.append(rec)
        return res
    finally:
        conn.close()


# ------------------------------
# Safety profile
# ------------------------------
def _fetch_top_faers(
    conn: sqlite3.Connection,
    table: str,
    struct_id: int,
    limit: int = 15,
) -> List[Dict[str, Any]]:
    allowed_tables = {
        "faers",
        "faers_female",
        "faers_male",
        "faers_ped",
        "faers_ger",
    }
    if table not in allowed_tables:
        raise ValueError(f"Unsupported FAERS table: {table!r}")

    cur = conn.cursor()
    cur.execute(
        f"""
        SELECT meddra_name,
               meddra_code,
               level,
               llr,
               llr_threshold,
               drug_ae,
               drug_no_ae,
               no_drug_ae,
               no_drug_no_ae
        FROM {table}
        WHERE struct_id = ?
        ORDER BY llr DESC
        LIMIT ?
        """,
        (struct_id, limit),
    )
    res: List[Dict[str, Any]] = []
    for r in cur.fetchall():
        res.append(
            {
                "meddra_name": r["meddra_name"] or "",
                "meddra_code": r["meddra_code"] or 0,
                "level": r["level"] or "",
                "llr": r["llr"] or 0.0,
                "llr_threshold": r["llr_threshold"] or 0.0,
                "drug_ae": r["drug_ae"] or 0,
                "drug_no_ae": r["drug_no_ae"] or 0,
                "no_drug_ae": r["no_drug_ae"] or 0,
                "no_drug_no_ae": r["no_drug_no_ae"] or 0,
            }
        )
    return res


def get_drug_safety_profile(drug_central_id: str) -> Optional[SafetyProfile]:
    """
    Aggregate safety-related information for a drug:
    FAERS signals (overall, sex/age strata) and simple DDI risks.

    All struct_id values correspond to structures.id.
    """
    sid = _parse_struct_id(drug_central_id)
    if sid is None:
        return None
    conn = get_connection()
    try:
        faers_all = _fetch_top_faers(conn, "faers", sid)
        faers_female = _fetch_top_faers(conn, "faers_female", sid)
        faers_male = _fetch_top_faers(conn, "faers_male", sid)
        faers_ped = _fetch_top_faers(conn, "faers_ped", sid)
        faers_ger = _fetch_top_faers(conn, "faers_ger", sid)

        # DDI data is mostly at drug_class level; here we just surface the raw entries
        ddi_cur = conn.cursor()
        ddi_cur.execute(
            """
            SELECT d.drug_class1,
                   d.drug_class2,
                   d.ddi_risk,
                   d.description
            FROM ddi d
            """,
        )
        ddi_risks = []
        for r in ddi_cur.fetchall():
            ddi_risks.append(
                {
                    "drug_class1": r["drug_class1"] or "",
                    "drug_class2": r["drug_class2"] or "",
                    "ddi_risk": r["ddi_risk"] or "",
                    "description": r["description"] or "",
                }
            )

        return SafetyProfile(
            struct_id=sid,
            faers_signals=faers_all,
            faers_female=faers_female,
            faers_male=faers_male,
            faers_pediatric=faers_ped,
            faers_geriatrics=faers_ger,
            ddi_risks=ddi_risks,
        )
    finally:
        conn.close()


# ------------------------------
# Cross references
# ------------------------------
def get_cross_references(drug_central_id: str) -> Optional[CrossRefInfo]:
    """
    Return external identifiers and ATC / drug class info for a drug.

    identifier.struct_id, struct2atc.struct_id, struct2drgclass.struct_id
    all reference structures.id.
    """
    sid = _parse_struct_id(drug_central_id)
    if sid is None:
        return None
    conn = get_connection()
    try:
        cur = conn.cursor()
        # identifier table
        cur.execute(
            """
            SELECT id_type, identifier
            FROM identifier
            WHERE struct_id = ?
            """,
            (sid,),
        )
        ids = []
        for r in cur.fetchall():
            ids.append(
                {
                    "id_type": r["id_type"] or "",
                    "identifier": r["identifier"] or "",
                }
            )

        # ATC codes
        cur.execute(
            """
            SELECT atc_code
            FROM struct2atc
            WHERE struct_id = ?
            """,
            (sid,),
        )
        atc_codes = [r["atc_code"] for r in cur.fetchall()]

        # Drug classes
        cur.execute(
            """
            SELECT dc.name
            FROM struct2drgclass s2d
            JOIN drug_class dc ON s2d.drug_class_id = dc.id
            WHERE s2d.struct_id = ?
            """,
            (sid,),
        )
        drug_classes = [r["name"] for r in cur.fetchall()]

        return CrossRefInfo(
            struct_id=sid,
            identifiers=ids,
            atc_codes=atc_codes,
            drug_classes=drug_classes,
        )
    finally:
        conn.close()


# ------------------------------
# Target class summary
# ------------------------------
def get_target_class_summary(drug_central_id: str) -> Optional[TargetClassSummary]:
    """
    Summarize MoA and secondary targets by target class and TDL for a drug.

    act_table_full.struct_id references structures.id.
    """
    sid = _parse_struct_id(drug_central_id)
    if sid is None:
        return None
    conn = get_connection()
    try:
        cur = conn.cursor()
        cur.execute(
            """
            SELECT target_class, tdl, moa
            FROM act_table_full
            WHERE struct_id = ?
            """,
            (sid,),
        )
        by_class: Counter[str] = Counter()
        by_tdl: Counter[str] = Counter()
        total_targets = 0
        moa_targets = 0
        for r in cur.fetchall():
            total_targets += 1
            tclass = r["target_class"] or "Unknown"
            tdl = r["tdl"] or "Unknown"
            by_class[tclass] += 1
            by_tdl[tdl] += 1
            if r["moa"]:
                moa_targets += 1

        if total_targets == 0:
            return TargetClassSummary(
                struct_id=sid,
                total_targets=0,
                by_class={},
                by_tdl={},
                moa_targets=0,
            )

        return TargetClassSummary(
            struct_id=sid,
            total_targets=total_targets,
            by_class=dict(by_class),
            by_tdl=dict(by_tdl),
            moa_targets=moa_targets,
        )
    finally:
        conn.close()


# ------------------------------
# Search and list helpers
# ------------------------------
def search_drugs(query: str, limit: int = 20) -> List[DrugRecord]:
    """
    General search over names, synonyms, identifiers, and basic indications.
    Intended for UI autocomplete and broad discovery.

    All returned struct_id values are structures.id.
    """
    if not query:
        return []
    q = query.strip().lower()
    conn = get_connection()
    try:
        cur = conn.cursor()
        # Search structures.name (structures.id) and CAS
        cur.execute(
            """
            SELECT DISTINCT id
            FROM structures
            WHERE lower(name) LIKE ?
               OR lower(cas_reg_no) = ?
            LIMIT ?
            """,
            (f"%{q}%", query.strip().lower(), limit),
        )
        struct_ids = [r["id"] for r in cur.fetchall()]

        # If not enough hits, add synonyms
        if len(struct_ids) < limit:
            needed = limit - len(struct_ids)
            cur.execute(
                """
                SELECT DISTINCT st.id
                FROM synonyms s
                JOIN structures st ON s.id = st.id
                WHERE lower(s.name) LIKE ?
                LIMIT ?
                """,
                (f"%{q}%", needed),
            )
            struct_ids.extend(r["id"] for r in cur.fetchall())

        # If still not enough hits, add product names
        if len(struct_ids) < limit:
            needed = limit - len(struct_ids)
            cur.execute(
                """
                SELECT DISTINCT st.id
                FROM product p
                JOIN active_ingredient ai
                  ON p.ndc_product_code = ai.ndc_product_code
                JOIN structures st
                  ON ai.struct_id = st.id
                WHERE lower(p.product_name) LIKE ?
                   OR lower(p.generic_name) LIKE ?
                LIMIT ?
                """,
                (f"%{q}%", f"%{q}%", needed),
            )
            struct_ids.extend(r["id"] for r in cur.fetchall())

        # If still not enough hits, add identifiers (ChEMBL, PubChem CID, RxNorm, etc.)
        if len(struct_ids) < limit:
            needed = limit - len(struct_ids)
            cur.execute(
                """
                SELECT DISTINCT struct_id
                FROM identifier
                WHERE lower(identifier) LIKE ?
                LIMIT ?
                """,
                (f"%{q}%", needed),
            )
            struct_ids.extend(r["struct_id"] for r in cur.fetchall())

        # Deduplicate
        struct_ids = list(dict.fromkeys(struct_ids))[:limit]

        results: List[DrugRecord] = []
        for sid in struct_ids:
            rec = _fetch_struct_by_struct_id(conn, sid)
            if rec:
                results.append(rec)
        return results
    finally:
        conn.close()


def list_approved_drugs(
    agency_filter: Optional[List[str]] = None,
    year_range: Optional[Tuple[int, int]] = None,
) -> List[DrugRecord]:
    """
    Return approved drugs, optionally filtered by regulatory agency
    (approval.type) and approval year range.

    approval.struct_id references structures.id.
    """
    conn = get_connection()
    try:
        cur = conn.cursor()
        sql = """
            SELECT DISTINCT struct_id, approval
            FROM approval
        """
        params: List[Any] = []
        filters: List[str] = []

        if agency_filter:
            placeholders = ",".join("?" for _ in agency_filter)
            filters.append(f"type IN ({placeholders})")
            params.extend(agency_filter)

        if year_range:
            start, end = year_range
            filters.append("CAST(substr(approval, 1, 4) AS INTEGER) BETWEEN ? AND ?")
            params.extend([start, end])

        if filters:
            sql += " WHERE " + " AND ".join(filters)

        cur.execute(sql, tuple(params))
        struct_ids = [r["struct_id"] for r in cur.fetchall()]
        results: List[DrugRecord] = []
        for sid in struct_ids:
            rec = _fetch_struct_by_struct_id(conn, sid)
            if rec:
                results.append(rec)
        return results
    finally:
        conn.close()


# ------------------------------
# Health check
# ------------------------------
def health_check() -> HealthStatus:
    """
    DB-level health check. Confirms presence of the file, connectability,
    version, and simple row counts from a few key tables.
    """
    if not DRUG_CENTRAL_DB_PATH.exists():
        return HealthStatus(
            ok=False,
            message=f"Drug Central DB file not found at {DRUG_CENTRAL_DB_PATH}",
            db_path=str(DRUG_CENTRAL_DB_PATH),
            version=None,
            tables_checked={},
        )

    conn = get_connection()
    try:
        cur = conn.cursor()
        # Version
        version = None
        try:
            cur.execute("SELECT version FROM dbversion ORDER BY version DESC LIMIT 1")
            row = cur.fetchone()
            if row:
                version = row["version"]
        except Exception:
            version = None

        tables_checked: Dict[str, int] = {}
        for table in ["structures", "act_table_full", "omop_relationship", "pharma_class"]:
            try:
                cur.execute(f"SELECT COUNT(1) AS cnt FROM {table}")
                row = cur.fetchone()
                tables_checked[table] = row["cnt"] if row else 0
            except Exception:
                tables_checked[table] = -1

        return HealthStatus(
            ok=True,
            message="Drug Central DB is present and basic tables are queryable.",
            db_path=str(DRUG_CENTRAL_DB_PATH),
            version=version,
            tables_checked=tables_checked,
        )
    finally:
        conn.close()


# ------------------------------
# Summary formatting for UI
# ------------------------------
def make_drug_summary_text(
    drug: DrugRecord,
    regulatory: Optional[RegulatoryInfo],
    targets: List[TargetRecord],
    indications: List[IndicationRecord],
    safety: Optional[SafetyProfile],
    pharm_actions: Optional[List[PharmActionRecord]] = None,
    human_gene_targets: Optional[List[HumanGeneTarget]] = None,
) -> str:
    """
    Create a compact, human-readable summary text for a single drug.

    Sections:
      - Core identifiers (name, Drug ID = structures.id, ChEMBL ID, CAS, InChIKey, SMILES)
      - Regulatory status, approvals, and Orange Book products
      - Pharmacologic actions (pharma_class)
      - Target / bioactivity overview
      - Human gene targets (Drug Central, HUMAN only)
      - Drug use (indications via OMOP)
      - Safety and DDI summary
    """
    lines: List[str] = []

    name = drug.get("name") or "Unknown drug"
    struct_id = drug.get("struct_id", "")
    lines.append(f"Drug Central record for: {name}")
    lines.append(f"  Drug ID (structures.id): {struct_id}")

    chembl_id = drug.get("chembl_id") or ""
    if chembl_id:
        lines.append(f"  ChEMBL ID: {chembl_id}")

    cas = drug.get("cas_number") or ""
    if cas:
        lines.append(f"  CAS: {cas}")
    inchikey = drug.get("inchikey") or ""
    if inchikey:
        lines.append(f"  InChIKey: {inchikey}")
    smiles = drug.get("smiles") or ""
    if smiles:
        lines.append(f"  SMILES: {smiles}")

    # Regulatory
    if regulatory:
        status = regulatory.get("status") or ""
        if status:
            lines.append(f"  Regulatory status: {status}")
        if regulatory.get("approvals"):
            lines.append("  Approvals:")
            for app in regulatory["approvals"]:
                date = app.get("approval") or ""
                agency = app.get("type") or ""
                company = app.get("applicant") or ""
                orphan = " (orphan)" if app.get("orphan") else ""
                detail = ", ".join(p for p in [date, agency, company] if p)
                if detail:
                    lines.append(f"    - {detail}{orphan}")
        if regulatory.get("orange_book_products"):
            lines.append("  Pharmaceutical products (Orange Book):")
            for prod in regulatory["orange_book_products"][:5]:
                trade = prod.get("trade_name") or ""
                ingr = prod.get("ingredient") or ""
                route = prod.get("route") or ""
                form = prod.get("dose_form") or ""
                strength = prod.get("strength") or ""
                appl_type = prod.get("appl_type") or ""
                appl_no = prod.get("appl_no") or ""
                product_no = prod.get("product_no") or ""
                category = prod.get("type") or ""
                company = prod.get("applicant") or ""
                rld_flag = " (RLD)" if prod.get("rld") else ""

                label = trade or ingr or "Product"
                meta_parts: List[str] = []
                if strength:
                    meta_parts.append(strength)
                if form:
                    meta_parts.append(form)
                if route:
                    meta_parts.append(route)
                if company:
                    meta_parts.append(company)
                appl_id = "-".join(x for x in [appl_type, appl_no, product_no] if x)
                if appl_id:
                    meta_parts.append(f"Application: {appl_id}")
                if category:
                    meta_parts.append(f"Category: {category}")
                meta = "; ".join(meta_parts)

                lines.append(f"    - {label}{rld_flag}")
                if meta:
                    lines.append(f"         {meta}")

    # Pharmacologic actions
    if pharm_actions:
        lines.append("  Pharmacologic actions (pharma_class):")
        for pa in pharm_actions:
            src = pa.get("source") or ""
            typ = pa.get("type") or ""
            code = pa.get("class_code") or ""
            pname = pa.get("name") or ""
            left = " ".join([p for p in [src, typ] if p])
            right = " - ".join([c for c in [code, pname] if c])
            if left and right:
                lines.append(f"    - {left}: {right}")
            elif right:
                lines.append(f"    - {right}")
            elif left:
                lines.append(f"    - {left}")

    # Targets / bioactivity (all organisms, as-is)
    if targets:
        lines.append("  Targets and bioactivity:")
        # Group by (gene, target_name)
        grouped: Dict[Tuple[str, str], List[TargetRecord]] = defaultdict(list)
        for t in targets:
            key = (t.get("gene") or "", t.get("target_name") or "")
            grouped[key].append(t)
        for (gene, tname), recs in grouped.items():
            gene_part = gene or "Unknown gene"
            tname_part = tname or ""
            head = f"    - {gene_part}"
            if tname_part and tname_part != gene_part:
                head += f" ({tname_part})"
            moa_count = sum(1 for r in recs if r.get("is_moa"))
            total = len(recs)
            head += f" [{moa_count} MoA / {total} total interactions]"
            lines.append(head)
            # up to a few representative activities
            for r in recs[:3]:
                act_type = r.get("act_type") or ""
                act_val = r.get("act_value")
                act_unit = r.get("act_unit") or ""
                rel = r.get("relation") or ""
                atype = r.get("action_type") or ""
                tclass = r.get("target_class") or ""
                tdl = r.get("tdl") or ""
                acc = r.get("accession") or ""
                sw = r.get("swissprot") or ""

                desc_parts: List[str] = []
                if atype:
                    desc_parts.append(atype)
                if act_type:
                    desc_parts.append(act_type)
                if act_val:
                    desc_parts.append(f"{rel} {act_val} {act_unit}".strip())
                if tclass:
                    desc_parts.append(f"class={tclass}")
                if tdl:
                    desc_parts.append(f"TDL={tdl}")
                if acc or sw:
                    acc_part = "/".join(p for p in [acc, sw] if p)
                    desc_parts.append(acc_part)

                desc = " ".join(p for p in desc_parts if p).strip()
                if desc:
                    lines.append(f"         {desc}")

    # Human gene targets (Drug Central, HUMAN-only, consistent with entity extraction)
    if human_gene_targets:
        lines.append("  Human gene targets (Drug Central, HUMAN only):")
        for tgt in human_gene_targets:
            symbol = tgt.get("symbol") or ""
            uniprot = tgt.get("uniprot_id") or ""
            entrez = tgt.get("entrez_id") or ""
            source = tgt.get("source") or ""
            score = tgt.get("score")
            n_act = tgt.get("n_activities")
            n_moa = tgt.get("n_moa")

            label = symbol or uniprot or "Unknown gene"
            meta_parts: List[str] = []
            if uniprot:
                meta_parts.append(f"UniProt: {uniprot}")
            if entrez:
                meta_parts.append(f"Entrez: {entrez}")
            if isinstance(score, (int, float)):
                meta_parts.append(f"score={score:.1f}")
            if isinstance(n_act, int) and n_act > 0:
                if isinstance(n_moa, int) and n_moa >= 0:
                    meta_parts.append(f"MoA/total: {n_moa}/{n_act}")
                else:
                    meta_parts.append(f"n_activities={n_act}")
            if source:
                meta_parts.append(f"source={source}")

            if meta_parts:
                lines.append(f"    - {label} ({'; '.join(meta_parts)})")
            else:
                lines.append(f"    - {label}")

    # Drug use / indications
    if indications:
        lines.append("  Drug use (indications via OMOP concepts):")
        seen = set()
        for ind in indications:
            cname = ind.get("concept_name") or ""
            if not cname or cname in seen:
                continue
            seen.add(cname)
            rel = ind.get("relationship_name") or ""
            cui = ind.get("umls_cui") or ""
            snomed = ind.get("snomed_conceptid") or 0
            extra: List[str] = []
            if rel:
                extra.append(rel)
            if cui:
                extra.append(f"UMLS CUI: {cui}")
            if snomed:
                extra.append(f"SNOMED CT: {snomed}")
            if extra:
                lines.append(f"    - {cname} ({'; '.join(extra)})")
            else:
                lines.append(f"    - {cname}")

    # Safety
    if safety:
        if safety.get("faers_signals"):
            lines.append("  Safety signals (FAERS, top by LLR):")
            for s in safety["faers_signals"][:10]:
                name_ = s.get("meddra_name") or ""
                llr = s.get("llr") or 0.0
                lvl = s.get("level") or ""
                lines.append(f"    - {name_} (LLR={llr:.2f}, level={lvl})")
        if safety.get("ddi_risks"):
            lines.append("  Drug-drug interaction risk classes (summary):")
            risk_counts: Counter[str] = Counter()
            for d in safety["ddi_risks"]:
                risk = d.get("ddi_risk") or ""
                if risk:
                    risk_counts[risk] += 1
            for risk, cnt in risk_counts.most_common(5):
                lines.append(f"    - {risk} (n={cnt} class pairs)")

    if len(lines) <= 2:
        lines.append("  No Drug Central annotations found for this record.")

    return "\n".join(lines)


def interpretation_notes(include_safety: bool = True) -> str:
    """
    Interpretation notes for Drug Central-based evidence.

    Acronyms are expanded for clarity and to emphasize that these
    summaries are heuristic and not a substitute for full label
    and primary data review.
    """
    lines: List[str] = []
    lines.append("")
    lines.append("Drug Central Interpretation Notes")
    lines.append("")
    lines.append(
        "- Drug Central struct_id in this tool corresponds to structures.id, "
        "the stable per-drug identifier used across Drug Central tables. "
        "The cd_id field in structures is an internal compound dictionary index."
    )
    lines.append(
        "- Structural properties such as logP (octanol/water partition coefficient), "
        "TPSA (topological polar surface area), Lipinski rule-of-five violations, "
        "rotatable bond count, and related features are read from the structures table."
    )
    lines.append(
        "- Targets, bioactivity measurements, and mechanism-of-action (MoA) versus "
        "off-target assignments are summarized from the act_table_full table."
    )
    lines.append(
        "- Gene symbols, UniProt accessions, and target development level (TDL, as "
        "used in the Pharos framework) are taken from target_dictionary and "
        "target_component mappings when available."
    )
    lines.append(
        "- Indications and contraindications use the OMOP (Observational Medical "
        "Outcomes Partnership) relationship tables, which map drugs to UMLS concept "
        "identifiers (CUI) and SNOMED CT concept identifiers."
    )
    if include_safety:
        lines.append(
            "- Safety profiles use FAERS (FDA Adverse Event Reporting System) "
            "disproportionality statistics such as LLR (log-likelihood ratio). "
            "These signals are hypothesis-generating and do not represent incidence "
            "or risk in the general population."
        )
        lines.append(
            "- Drug-drug interaction (DDI) information is largely captured at the level "
            "of drug classes and should be interpreted cautiously, in conjunction with "
            "formal labeling and clinical judgment."
        )
    lines.append(
        "- Pharmacologic actions (for example MOA = mechanism of action, EPC = "
        "established pharmacologic class, PE = pharmacologic effect) are derived "
        "from the pharma_class annotations."
    )
    lines.append(
        "- External cross-references such as ChEMBL identifiers, PubChem compound IDs, "
        "RxNorm codes, and ATC (Anatomical Therapeutic Chemical) classes are obtained "
        "from the identifier, struct2atc, and related mapping tables."
    )
    lines.append(
        "- Disease Ontology identifiers (DOID) are not present in the current Drug "
        "Central schema and would require an external mapping layer if needed."
    )
    lines.append(
        "- All tool-generated summaries are heuristic views over the underlying "
        "Drug Central tables and are not a substitute for detailed label review, "
        "regulatory documents, or experimental validation."
    )
    return "\n".join(lines)


# ------------------------------
# Frontend asset injection
# ------------------------------
def inject_frontend_assets(
    html_template: str,
    css_template: str,
    js_template: str,
    drug_central_data_all: Dict[str, Dict[str, Any]],
    drugs_str: str,
) -> str:
    """
    Inject inline CSS, JS, and JSON Drug Central data into HTML template.
    Produces a standalone HTML viewer (parallel to chembl tool).
    Returns empty HTML if all drugs have empty evidence sets.

    Keys in drug_central_data_all are typically Drug Central struct_id
    (structures.id) rendered as strings.
    """
    has_data = False
    for _, data in (drug_central_data_all or {}).items():
        if not isinstance(data, dict):
            continue
        if any(
            data.get(key)
            for key in [
                "drug",
                "regulatory",
                "targets",
                "indications",
                "safety",
                "pharm_actions",
                "gene_targets_human",
            ]
        ):
            has_data = True
            break

    if not has_data:
        return ""

    css_inline = f"<style>\n{css_template}\n</style>"
    js_inline = (
        "<script>\n"
        f"var drug_centralData = {json.dumps(drug_central_data_all, indent=2)};\n"
        f"{js_template}\n</script>"
    )

    return (
        html_template.replace("{{CSS_INLINE}}", css_inline)
        .replace("{{JS_INLINE}}", js_inline)
        .replace("{{DRUG_LABEL}}", drugs_str)
    )
