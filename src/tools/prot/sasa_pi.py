# src/tools/prot/tool_prot_sasa_pi.py
# Author: Dmitri Kosenkov
# Created: 2025-08-25
# Updated: 2025-09-18
#
# SASA, Polarity Index (PI), and residue label retrieval and normalization.
# ASCII-only. No omissions.

from typing import Optional, Dict, List, Tuple

# Import utilities from same package
from src.tools.prot.utils import minmax_normalize, sanitize_residue_label


def fetch_sasa_pi(
    conn,
    protein_id: str
) -> Tuple[
    Optional[Dict[str, float]],  # sasa_stats
    List[Dict[str, float]],      # sasa_norm
    Optional[Dict[str, float]],  # pi_stats
    List[Dict[str, float]],      # pi_norm
    Dict[int, str]               # residue_labels
]:
    """
    Fetch SASA (area_total), Polarity Index (PI), and residue labels in one query.

    Returns:
      sasa_stats: {"min": float, "max": float, "avg": float, "total_area": float}
      sasa_norm:  [{"residue_no": int, "score": 0..1}, ...]
      pi_stats:   {"min": float, "max": float, "avg": float}
      pi_norm:    [{"residue_no": int, "score": 0..1}, ...]
      residue_labels: {residue_no: "AAA<no>"}

    Notes:
      - SASA values are min-max normalized for visualization.
      - PI values are min-max normalized for visualization but raw stats are preserved.
      - Residue labels are built as ASCII-safe strings (e.g. "MET1").
    """
    cur = conn.cursor()
    cur.execute(
        """
        SELECT residue_no, residue_name, area_total, pol_index_raw
        FROM data_sasa_pi_residues
        WHERE protein_id=?
        ORDER BY residue_no
        """,
        (protein_id,),
    )
    rows = cur.fetchall()
    if not rows:
        return None, [], None, [], {}

    sasa_rows: List[Tuple[int, float]] = []
    pi_rows: List[Tuple[int, float]] = []
    residue_labels: Dict[int, str] = {}

    for rn, resname, area_total, pi_val in rows:
        if rn is None:
            continue
        # Build label safely
        lab = sanitize_residue_label(resname, rn)
        if lab:
            residue_labels[int(rn)] = lab
        if area_total is not None:
            sasa_rows.append((rn, area_total))
        if pi_val is not None:
            pi_rows.append((rn, pi_val))

    sasa_stats, sasa_norm = minmax_normalize(sasa_rows) if sasa_rows else (None, [])
    if sasa_stats:
        sasa_stats["total_area"] = sum(v for _, v in sasa_rows)

    pi_stats, pi_norm = minmax_normalize(pi_rows) if pi_rows else (None, [])

    return sasa_stats, sasa_norm, pi_stats, pi_norm, residue_labels
