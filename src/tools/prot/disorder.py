# src/tools/prot/disorder.py
# Author: Dmitri Kosenkov
# Created: 2025-09-25
# Updated: 2025-09-27
#
# Consensus disorder analysis manager.
# Uses IUPred3 (short/long) + DisProt regions.
# MoRF propensity from IUPred3 ANCHOR2 predictions.

import sys
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# --- Package/Repo Roots ---
PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.tools.prot.utils import compute_stats, get_connection, log
from src.tools.prot.iupred import fetch_iupred3
from src.tools.prot.disprot import (
    fetch_disprot_for_uniprot,
    merge_disprot_regions,
    residue_in_disprot,
)
from src.tools.prot import OUTPUT_DIR


def fetch_disorder(
    conn,
    protein_id: str,
    uniprot_id: str,
    *,
    disorder_threshold: float = 0.5,
    morf_threshold: float = 0.5,
    min_segment_len: int = 5,
) -> Tuple[
    Optional[Dict[str, float]],  # disorder_stats (rich)
    List[Dict[str, float]],      # disorder_norm
    Optional[Dict[str, float]],  # morf_stats (rich)
    List[Dict[str, float]],      # morf_norm
]:
    """
    Build disorder and MoRF tracks and summarize with informative metrics.

    Returns:
      (disorder_stats, disorder_norm, morf_stats, morf_norm)

    disorder_stats keys:
      - min, max                             (for legend range)
      - pct_disordered                       (fraction of residues >= threshold, %)
      - num_regions                          (count of contiguous segments >= threshold with length >= min_segment_len)
      - longest_region                       (length of the longest such segment)
      - avg_short, avg_long                  (mean IUPred3 short/long scores where available)
      - n_disprot_regions                    (# of merged DisProt disorder regions used)
      - threshold                            (echo of disorder_threshold)
      - min_segment_len                      (echo of min_segment_len)

    morf_stats keys:
      - min, max                             (for legend range)
      - pct_high                             (% residues with MoRF >= threshold)
      - num_segments                         (contiguous segments at/above threshold with length >= min_segment_len)
      - longest_segment                      (length of the longest such segment)
      - threshold                            (echo of morf_threshold)
      - min_segment_len                      (echo of min_segment_len)
    """
    iupred = fetch_iupred3(conn, protein_id)
    if not iupred:
        return None, [], None, []

    # Collect all residue numbers present in any track
    resnos = {row["residue_no"] for _, lst in iupred.values() for row in lst}
    if not resnos:
        return None, [], None, []

    # DisProt consensus regions (optional, may be empty)
    dis_data = fetch_disprot_for_uniprot(uniprot_id)
    dis_regs = merge_disprot_regions(dis_data["consensus"]["regions"]) if dis_data else []

    # Fast lookup tables for IUPred/ANCHOR scores
    short_map = _to_map(iupred.get("iupred_short", (None, []))[1])
    long_map = _to_map(iupred.get("iupred_long", (None, []))[1])
    anch_s_map = _to_map(iupred.get("anchor_short", (None, []))[1])
    anch_l_map = _to_map(iupred.get("anchor_long", (None, []))[1])

    disorder_norm: List[Dict[str, float]] = []
    morf_norm: List[Dict[str, float]] = []
    disorder_scores: List[float] = []
    morf_scores: List[float] = []

    # Build per-residue consensus disorder and MoRF
    for res in sorted(resnos):
        s_short = short_map.get(res)
        s_long = long_map.get(res)
        a_short = anch_s_map.get(res)
        a_long = anch_l_map.get(res)

        in_disprot, amb = residue_in_disprot(res, dis_regs)

        # Consensus disorder: average available (short/long) + DisProt vote (1.0 or 0.5 if ambiguous)
        d_score = _fuse_disorder(s_short, s_long, in_disprot, amb)
        disorder_norm.append({"residue_no": res, "score": d_score})
        disorder_scores.append(d_score)

        # MoRF propensity: take maximum of available ANCHOR2 tracks
        m_score = max(a for a in [a_short, a_long, 0.0] if a is not None)
        morf_norm.append({"residue_no": res, "score": m_score})
        morf_scores.append(m_score)

    # Summaries
    disorder_stats = _make_disorder_stats(
        disorder_scores=disorder_scores,
        resnos_sorted=sorted(resnos),
        short_map=short_map,
        long_map=long_map,
        disprot_regions=dis_regs,
        threshold=disorder_threshold,
        min_segment_len=min_segment_len,
    )
    morf_stats = _make_morf_stats(
        morf_scores=morf_scores,
        threshold=morf_threshold,
        min_segment_len=min_segment_len,
    )

    return disorder_stats, disorder_norm, morf_stats, morf_norm


# ------------------------------
# Internals / helpers
# ------------------------------
def _to_map(rows: List[Dict[str, float]]) -> Dict[int, float]:
    """List[{'residue_no', 'score'}] -> {residue_no: score}"""
    return {int(r["residue_no"]): float(r["score"]) for r in rows or []}


def _lookup(iupred: Dict, key: str, resno: int) -> Optional[float]:
    """(Kept for compatibility; not used in the new path)"""
    if key not in iupred:
        return None
    _, lst = iupred[key]
    for row in lst:
        if row["residue_no"] == resno:
            return row["score"]
    return None


def _fuse_disorder(
    s_short: Optional[float],
    s_long: Optional[float],
    in_disprot: bool,
    ambiguous: bool,
) -> float:
    """Fuse predictors into a consensus score (0..1)."""
    vals: List[float] = []
    if s_short is not None:
        vals.append(s_short)
    if s_long is not None:
        vals.append(s_long)
    if in_disprot:
        vals.append(1.0 if not ambiguous else 0.5)
    if not vals:
        return 0.0
    avg = sum(vals) / float(len(vals))
    return max(0.0, min(1.0, avg))


def _segment_metrics_from_scores(
    scores: List[float],
    *,
    threshold: float,
    min_len: int,
) -> Tuple[int, int]:
    """
    From a list of per-residue scores aligned along the sequence, compute:
      - number of contiguous segments with score >= threshold and length >= min_len
      - longest such segment length
    """
    count = 0
    longest = 0
    cur_len = 0

    for s in scores:
        if s is not None and s >= threshold:
            cur_len += 1
        else:
            if cur_len >= min_len:
                count += 1
                longest = max(longest, cur_len)
            cur_len = 0

    # tail segment
    if cur_len >= min_len:
        count += 1
        longest = max(longest, cur_len)

    return count, longest


def _make_disorder_stats(
    disorder_scores: List[float],
    resnos_sorted: List[int],
    short_map: Dict[int, float],
    long_map: Dict[int, float],
    disprot_regions: List[Dict],
    *,
    threshold: float,
    min_segment_len: int,
) -> Optional[Dict[str, float]]:
    """Summarize consensus disorder with richer insight."""
    if not disorder_scores:
        return None

    n_res = len(disorder_scores)
    min_v = min(disorder_scores)
    max_v = max(disorder_scores)

    # fraction disordered at threshold
    n_dis = sum(1 for s in disorder_scores if s is not None and s >= threshold)
    pct_dis = (100.0 * n_dis / n_res) if n_res else 0.0

    # segments / longest on consensus
    num_regions, longest_region = _segment_metrics_from_scores(
        disorder_scores, threshold=threshold, min_len=min_segment_len
    )

    # means of IUPred short/long where available
    if short_map:
        avg_short = sum(short_map.values()) / float(len(short_map))
    else:
        avg_short = None
    if long_map:
        avg_long = sum(long_map.values()) / float(len(long_map))
    else:
        avg_long = None

    return {
        "min": min_v,
        "max": max_v,
        "pct_disordered": pct_dis,
        "num_regions": num_regions,
        "longest_region": longest_region,
        "avg_short": avg_short,
        "avg_long": avg_long,
        "n_disprot_regions": len(disprot_regions),
        "threshold": threshold,
        "min_segment_len": float(min_segment_len),
    }


def _make_morf_stats(
    morf_scores: List[float],
    *,
    threshold: float,
    min_segment_len: int,
) -> Optional[Dict[str, float]]:
    if not morf_scores:
        return None

    min_v = min(morf_scores)
    max_v = max(morf_scores)
    n_res = len(morf_scores)
    n_high = sum(1 for s in morf_scores if s is not None and s >= threshold)
    pct_high = (100.0 * n_high / n_res) if n_res else 0.0

    num_seg, longest_seg = _segment_metrics_from_scores(
        morf_scores, threshold=threshold, min_len=min_segment_len
    )

    return {
        "min": min_v,
        "max": max_v,
        "pct_high": pct_high,
        "num_segments": num_seg,
        "longest_segment": longest_seg,
        "threshold": threshold,
        "min_segment_len": float(min_segment_len),
    }


# ----------------------------------------------------------------------
# CLI (testing)
# ----------------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 3:
        protein_id = "P04637_F1"   # TP53 example
        uniprot_id = "P04637"
        print(f"[INFO] Using defaults: {protein_id} {uniprot_id}")
    else:
        protein_id = sys.argv[1]
        uniprot_id = sys.argv[2]

    conn = get_connection()
    disorder_stats, disorder_norm, morf_stats, morf_norm = fetch_disorder(conn, protein_id, uniprot_id)

    base_name = f"{protein_id}_{uniprot_id}_disorder"
    txt_out = OUTPUT_DIR / f"{base_name}.txt"
    json_out = OUTPUT_DIR / f"{base_name}.json"

    with txt_out.open("w", encoding="utf-8") as f:
        f.write(f"Disorder summary for {protein_id} ({uniprot_id})\n")
        f.write(f"Disorder stats: {json.dumps(disorder_stats, ensure_ascii=True)}\n")
        f.write(f"MoRF stats: {json.dumps(morf_stats, ensure_ascii=True)}\n\n")
        f.write("First 10 consensus scores:\n")
        for row in disorder_norm[:10]:
            f.write(f"  Residue {row['residue_no']}: {row['score']:.3f}\n")
        f.write("\nFirst 10 MoRF scores:\n")
        for row in morf_norm[:10]:
            f.write(f"  Residue {row['residue_no']}: {row['score']:.3f}\n")

    with json_out.open("w", encoding="utf-8") as f:
        json.dump(
            {
                "stats": disorder_stats,
                "consensus": disorder_norm,
                "morf_stats": morf_stats,
                "morf": morf_norm,
            },
            f,
            ensure_ascii=True,
            indent=2,
            sort_keys=True,
        )

    print("[OK] Disorder analysis complete.")
    print(f"  * TXT: {txt_out.resolve()}")
    print(f"  * JSON: {json_out.resolve()}")
    print(f"[INFO] Residues processed: {len(disorder_norm)}")
