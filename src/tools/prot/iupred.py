# src/tools/prot/iupred.py
# Updated: 2025-09-27
#
# IUPred3 / ANCHOR2 loader + lightweight metrics per track.

from typing import Optional, Dict, List, Tuple
from src.tools.prot.utils import compute_stats

TrackStats = Dict[str, float]
TrackList = List[Dict[str, float]]
IuPredOut = Dict[str, Tuple[Optional[TrackStats], TrackList]]


def _segments_above_threshold(values: List[float], thr: float = 0.5, min_len: int = 5):
    count = 0
    longest = 0
    cur = 0
    for v in values:
        if v is not None and v >= thr:
            cur += 1
        else:
            if cur >= min_len:
                count += 1
                if cur > longest:
                    longest = cur
            cur = 0
    if cur >= min_len:
        count += 1
        if cur > longest:
            longest = cur
    return count, longest


def _augmented_stats(vals: List[float]) -> TrackStats:
    """Return min/max/avg plus practical binary metrics at 0.5 threshold."""
    base = compute_stats(vals) or {"min": None, "max": None, "avg": None}
    n = len(vals)
    ge = sum(1 for v in vals if v is not None and v >= 0.5)
    frac_ge = (ge / n) if n else 0.0
    num_seg, longest = _segments_above_threshold(vals, thr=0.5, min_len=5)
    base.update({
        "frac_ge_0_5": frac_ge,
        "num_segments_ge_0_5": float(num_seg),
        "longest_segment_ge_0_5": float(longest),
    })
    return base


def fetch_iupred3(conn, protein_id: str) -> IuPredOut:
    """
    Fetch IUPred3 and ANCHOR2 per-residue scores.
    Table: data_iupred3_residues
      - iupred_short, iupred_long in [0..1]
      - anchor_short, anchor_long in [0..1]

    Returns dict with keys:
      "iupred_short", "iupred_long", "anchor_short", "anchor_long"
    Each maps to (stats, [{"residue_no": int, "score": float}, ...])

    Notes:
      - 'stats' now includes min/max/avg as well as:
          frac_ge_0_5, num_segments_ge_0_5, longest_segment_ge_0_5
        (Backwards compatible: existing consumers reading min/max/avg still work.)
    """
    cur = conn.cursor()
    cur.execute(
        """
        SELECT residue_no, iupred_short, anchor_short, iupred_long, anchor_long
        FROM data_iupred3_residues
        WHERE protein_id=?
        ORDER BY residue_no
        """,
        (protein_id,),
    )
    rows = cur.fetchall()
    if not rows:
        return {}

    buckets = {
        "iupred_short": [],
        "anchor_short": [],
        "iupred_long": [],
        "anchor_long": [],
    }
    norm = {k: [] for k in buckets.keys()}

    for rn, s_short, a_short, s_long, a_long in rows:
        if s_short is not None:
            f = float(s_short)
            buckets["iupred_short"].append(f)
            norm["iupred_short"].append({"residue_no": int(rn), "score": f})
        if a_short is not None:
            f = float(a_short)
            buckets["anchor_short"].append(f)
            norm["anchor_short"].append({"residue_no": int(rn), "score": f})
        if s_long is not None:
            f = float(s_long)
            buckets["iupred_long"].append(f)
            norm["iupred_long"].append({"residue_no": int(rn), "score": f})
        if a_long is not None:
            f = float(a_long)
            buckets["anchor_long"].append(f)
            norm["anchor_long"].append({"residue_no": int(rn), "score": f})

    out: IuPredOut = {}
    for k, vals in buckets.items():
        out[k] = (_augmented_stats(vals) if vals else None, norm[k] if vals else [])
    return out
