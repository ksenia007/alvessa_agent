# src/tools/prot/cysdb.py
# Author: Dmitri Kosenkov
# Created: 2025-12-08
#
# CysDB cysteine chemoproteomics and functional context loader.
#
# Table: data_cys_db_residues
#   protein_id TEXT,
#   uniprot_id TEXT,
#   residue_no INTEGER,
#   identified INTEGER,
#   hyperreactive INTEGER,
#   ligandable INTEGER,
#   is_act_site INTEGER,
#   near_act_site INTEGER,
#   near_active_site_neighbors TEXT,
#   is_bind_site INTEGER,
#   near_bind_site INTEGER,
#   near_bind_site_neighbors TEXT
#
# Returns:
#   - stats: per-protein summary including counts, residue labels, neighbor
#            summaries, and include_cysdb flag.
#   - tracks: per-residue binary tracks for frontend integration:
#       "cysdb_hyperreactive"
#       "cysdb_ligandable"
#       "cysdb_is_act_site"
#       "cysdb_near_act_site"
#       "cysdb_is_bind_site"
#       "cysdb_near_bind_site"
#
# Notes:
#   - Residue numbering is assumed to be 1-based and aligned with other
#     per-residue tables for the same protein_id.
#   - Neighbor lists (near_act_site_neighbors, near_bind_site_neighbors) are
#     global, deduplicated summaries, not per-residue maps.

from typing import Dict, List, Tuple, Any, Optional

from src.tools.prot.utils import uniq_preserve_order

CysdbStats = Dict[str, Any]
CysdbTrack = List[Dict[str, int]]
CysdbTracks = Dict[str, CysdbTrack]


def _flag_value(raw: Optional[int]) -> Optional[int]:
    """Normalize raw integer flags (0/1 or None) to 0/1/None."""
    if raw is None:
        return None
    try:
        iv = int(raw)
    except Exception:
        return None
    if iv == 0:
        return 0
    if iv == 1:
        return 1
    return None


def _label_cysteine(uniprot_id: str, residue_no: int) -> str:
    """Create labels like P12345_C11 using UniProt ID and residue number."""
    return f"{uniprot_id}_C{int(residue_no)}"

def fetch_cysdb(
    conn,
    protein_id: str,
    uniprot_id: str,
) -> Tuple[Optional[CysdbStats], CysdbTracks]:
    """
    Fetch CysDB cysteine annotations for a protein.

    Returns:
      stats: dict with counts, residue lists, neighbor summaries and
             per-residue neighbor maps, plus include_cysdb flag,
             or None if no rows found.
      tracks: dict of per-residue binary tracks for frontend use.

    Tracks (keys):
      "cysdb_hyperreactive"
      "cysdb_ligandable"
      "cysdb_is_act_site"
      "cysdb_near_act_site"
      "cysdb_is_bind_site"
      "cysdb_near_bind_site"
    """
    cur = conn.cursor()
    cur.execute(
        """
        SELECT
            residue_no,
            identified,
            hyperreactive,
            ligandable,
            is_act_site,
            near_act_site,
            near_active_site_neighbors,
            is_bind_site,
            near_bind_site,
            near_bind_site_neighbors
        FROM data_cys_db_residues
        WHERE protein_id=?
        ORDER BY residue_no
        """,
        (protein_id,),
    )
    rows = cur.fetchall()
    if not rows:
        return None, {}

    # Counts per flag
    n_identified = 0
    n_hyperreactive = 0
    n_ligandable = 0
    n_is_act_site = 0
    n_near_act_site = 0
    n_is_bind_site = 0
    n_near_bind_site = 0

    # Residue labels per flag for text output
    identified_residues: List[str] = []
    hyperreactive_residues: List[str] = []
    ligandable_residues: List[str] = []
    is_act_site_residues: List[str] = []
    near_act_site_residues: List[str] = []
    is_bind_site_residues: List[str] = []
    near_bind_site_residues: List[str] = []

    # Neighbor residue labels for global summaries
    near_act_neighbors_all: List[str] = []
    near_bind_neighbors_all: List[str] = []

    # Per-residue neighbor raw strings (for per-residue printing in summary)
    per_residue_near_act_site_neighbors: Dict[int, str] = {}
    per_residue_near_bind_site_neighbors: Dict[int, str] = {}

    # Per-residue binary tracks for frontend
    tracks: CysdbTracks = {
        "cysdb_hyperreactive": [],
        "cysdb_ligandable": [],
        "cysdb_is_act_site": [],
        "cysdb_near_act_site": [],
        "cysdb_is_bind_site": [],
        "cysdb_near_bind_site": [],
    }

    for (
        residue_no,
        identified_raw,
        hyperreactive_raw,
        ligandable_raw,
        is_act_site_raw,
        near_act_site_raw,
        near_act_neighbors_raw,
        is_bind_site_raw,
        near_bind_site_raw,
        near_bind_neighbors_raw,
    ) in rows:
        if residue_no is None:
            continue

        rn = int(residue_no)
        label = _label_cysteine(uniprot_id, rn)

        identified = _flag_value(identified_raw)
        hyperreactive = _flag_value(hyperreactive_raw)
        ligandable = _flag_value(ligandable_raw)
        is_act_site = _flag_value(is_act_site_raw)
        near_act_site = _flag_value(near_act_site_raw)
        is_bind_site = _flag_value(is_bind_site_raw)
        near_bind_site = _flag_value(near_bind_site_raw)

        # Counts and residue label lists
        if identified == 1:
            n_identified += 1
            identified_residues.append(label)
        if hyperreactive == 1:
            n_hyperreactive += 1
            hyperreactive_residues.append(label)
        if ligandable == 1:
            n_ligandable += 1
            ligandable_residues.append(label)
        if is_act_site == 1:
            n_is_act_site += 1
            is_act_site_residues.append(label)
        if near_act_site == 1:
            n_near_act_site += 1
            near_act_site_residues.append(label)
        if is_bind_site == 1:
            n_is_bind_site += 1
            is_bind_site_residues.append(label)
        if near_bind_site == 1:
            n_near_bind_site += 1
            near_bind_site_residues.append(label)

        # Tracks: store 0 or 1 explicitly for residues where data is present
        if hyperreactive is not None:
            tracks["cysdb_hyperreactive"].append(
                {"residue_no": rn, "score": int(hyperreactive)}
            )
        if ligandable is not None:
            tracks["cysdb_ligandable"].append(
                {"residue_no": rn, "score": int(ligandable)}
            )
        if is_act_site is not None:
            tracks["cysdb_is_act_site"].append(
                {"residue_no": rn, "score": int(is_act_site)}
            )
        if near_act_site is not None:
            tracks["cysdb_near_act_site"].append(
                {"residue_no": rn, "score": int(near_act_site)}
            )
        if is_bind_site is not None:
            tracks["cysdb_is_bind_site"].append(
                {"residue_no": rn, "score": int(is_bind_site)}
            )
        if near_bind_site is not None:
            tracks["cysdb_near_bind_site"].append(
                {"residue_no": rn, "score": int(near_bind_site)}
            )

        # Neighbors for global summary lists and per-residue maps
        if near_act_neighbors_raw:
            raw_str = str(near_act_neighbors_raw)
            tokens = [
                t.strip()
                for t in raw_str.split(";")
                if t and t.strip()
            ]
            near_act_neighbors_all.extend(tokens)
            per_residue_near_act_site_neighbors[rn] = raw_str

        if near_bind_neighbors_raw:
            raw_str = str(near_bind_neighbors_raw)
            tokens = [
                t.strip()
                for t in raw_str.split(";")
                if t and t.strip()
            ]
            near_bind_neighbors_all.extend(tokens)
            per_residue_near_bind_site_neighbors[rn] = raw_str

    # Deduplicate neighbor lists while preserving order
    near_act_neighbors_all = uniq_preserve_order(near_act_neighbors_all)
    near_bind_neighbors_all = uniq_preserve_order(near_bind_neighbors_all)

    has_any_flag = any(
        c > 0
        for c in (
            n_identified,
            n_hyperreactive,
            n_ligandable,
            n_is_act_site,
            n_near_act_site,
            n_is_bind_site,
            n_near_bind_site,
        )
    )

    # Per-residue neighbor maps keyed by residue label (e.g., P04406_C152)
    near_act_site_neighbors_per_residue: Dict[str, str] = {}
    near_bind_site_neighbors_per_residue: Dict[str, str] = {}
    for rn, neigh_raw in per_residue_near_act_site_neighbors.items():
        if neigh_raw:
            lab = _label_cysteine(uniprot_id, rn)
            near_act_site_neighbors_per_residue[lab] = neigh_raw
    for rn, neigh_raw in per_residue_near_bind_site_neighbors.items():
        if neigh_raw:
            lab = _label_cysteine(uniprot_id, rn)
            near_bind_site_neighbors_per_residue[lab] = neigh_raw

    stats: CysdbStats = {
        "n_identified": float(n_identified),
        "n_hyperreactive": float(n_hyperreactive),
        "n_ligandable": float(n_ligandable),
        "n_is_act_site": float(n_is_act_site),
        "n_near_act_site": float(n_near_act_site),
        "n_is_bind_site": float(n_is_bind_site),
        "n_near_bind_site": float(n_near_bind_site),
        "identified_residues": identified_residues,
        "hyperreactive_residues": hyperreactive_residues,
        "ligandable_residues": ligandable_residues,
        "is_act_site_residues": is_act_site_residues,
        "near_act_site_residues": near_act_site_residues,
        "is_bind_site_residues": is_bind_site_residues,
        "near_bind_site_residues": near_bind_site_residues,
        "near_act_site_neighbors": near_act_neighbors_all,
        "near_bind_site_neighbors": near_bind_neighbors_all,
        "near_act_site_neighbors_per_residue": near_act_site_neighbors_per_residue,
        "near_bind_site_neighbors_per_residue": near_bind_site_neighbors_per_residue,
        "has_any_flag": bool(has_any_flag),
        "include_cysdb": bool(has_any_flag),
    }

    return stats, tracks
