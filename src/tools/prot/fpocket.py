# src/tools/prot/tool_prot_fpocket.py
# Author: Dmitri Kosenkov
# Created: 2025-08-25
# Updated: 2025-09-18
#
# FPocket-specific data retrieval and normalization.

from typing import Optional, Dict, List, Tuple

# Import from same package
from src.tools.prot.utils import minmax_normalize


def fetch_fpocket(conn, protein_id: str) -> Tuple[Optional[Dict[str, float]], List[Dict[str, float]]]:
    """
    Fetch per-residue mean_druggability_score values for a given protein_id.

    Returns:
      (stats, normalized_list)

      stats: {"min": float, "max": float, "avg": float}
      normalized_list: [{"residue_no": int, "score": 0..1}, ...]

    Normalization:
      - Min-max normalization is applied to map values into [0..1].
    """
    cur = conn.cursor()
    cur.execute(
        """
        SELECT residue_no, mean_druggability_score
        FROM data_fpocket_residues
        WHERE protein_id=?
        ORDER BY residue_no
        """,
        (protein_id,),
    )
    rows = [(rn, v) for rn, v in cur.fetchall() if v is not None]
    if not rows:
        return None, []
    return minmax_normalize(rows)
