# src/tools/prot/tool_prot_plddt.py
# Author: Dmitri Kosenkov
# Created: 2025-08-25
# Updated: 2025-09-18
#
# pLDDT-specific data retrieval and normalization.

from typing import Optional, Dict, List, Tuple

# Import from same package
from src.tools.prot.utils import compute_stats


def fetch_plddt(conn, protein_id: str) -> Tuple[Optional[Dict[str, float]], List[Dict[str, float]]]:
    """
    Fetch per-residue pLDDT values for a given protein_id.

    Returns:
      (stats, normalized_list)

      stats: {"min": float, "max": float, "avg": float}
      normalized_list: [{"residue_no": int, "score": 0..1}, ...]

    Normalization:
      - Raw pLDDT is scaled from [0..100] to [0..1] for visualization.
    """
    cur = conn.cursor()
    cur.execute(
        """
        SELECT residue_no, plddt
        FROM data_plldt_residues
        WHERE protein_id=?
        ORDER BY residue_no
        """,
        (protein_id,),
    )
    rows = cur.fetchall()
    if not rows:
        return None, []

    vals = [r[1] for r in rows if r[1] is not None]
    if not vals:
        return None, []

    stats = compute_stats(vals)
    normalized = [{"residue_no": int(r[0]), "score": float(r[1]) / 100.0} for r in rows if r[1] is not None]
    return stats, normalized
