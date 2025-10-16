# alvessa_agent/evals/generation/prot_questions.py
# Author: Dmitri Kosenkov
# Updated: 2025-10-16
#
# Description
# -----------
# Generates multiple-choice validation questions for the `prot` tool,
# computing correct answers using structural data from the local protein DB.
#
# Input (TSV, tab-separated):
#   alvessa_agent/local_dbs/known_tool_prot_genes_list.tsv
#   Columns: gene_symbol  uniprot_id
#
# Output (default, strict 3-column TSV):
#   alvessa_agent/out/prot_validation_questions.tsv
#   Columns (tab-separated, with header): question    answer  tool
#

import sys
import time
import math
import json
import random
import warnings
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple, Callable

import pandas as pd

# =============================================================================
# Bootstrap: set paths, load .env (optional), then import prot utils
# =============================================================================
ROOT = Path(__file__).resolve().parents[2]  # alvessa_agent/
SRC = ROOT / "src"
OUT_DIR = ROOT / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))
if str(SRC) not in sys.path:
    sys.path.append(str(SRC))

try:
    from dotenv import load_dotenv  # optional dependency
    ENV = ROOT / ".env"
    if ENV.exists():
        load_dotenv(dotenv_path=ENV, override=False)
except Exception:
    pass

DEBUG = False

# prot helpers
from src.tools.prot.utils import get_connection
from src.tools.prot.plddt import fetch_plddt
from src.tools.prot.fpocket import fetch_fpocket
from src.tools.prot.sasa_pi import fetch_sasa_pi
from src.tools.prot.disorder import fetch_disorder
from src.tools.prot.biolip2 import fetch_biolip2

# Hard-pinned input TSV path
KNOWN_GENE_LIST_PATH = ROOT / "local_dbs" / "known_tool_prot_genes_list.tsv"

# =============================================================================
# Script-local config (simplified)
# =============================================================================
SEED = 42
MAX_ATTEMPTS_PER_Q = 60
TEMPLATE_MAX_PER_ID = 3  # exactly 3 per template (best effort)

# Margin constants (percent) for float metrics; do not apply to plddt_mean or biolip_n_sites
DISORDER_PCT_LOWEST_PER_CENT = 30.0
PI_MEAN_SURFACE_HIGHEST_PER_CENT = 30.0
SASA_AVG_HIGHEST_PER_CENT = 30.0
MORF_MAX_HIGHEST_PER_CENT = 30.0

# Threshold constant for pLDDT "above X"
PLDDT_MEAN_GT_THRESHOLD = 70.0

# BioLiP2 threshold constant for ">= N sites" condition
BIOLIP_SITES_AT_LEAST = 1

# =============================================================================
# Logging
# =============================================================================
def _now_str() -> str:
    import datetime as _dt
    return _dt.datetime.now().strftime("%H:%M:%S")

class Progress:
    """
    Minimal progress tracker. Reports each attempt only when DEBUG is True.
    """
    def __init__(self, name: str, target: int):
        self.name = name
        self.target = target
        self.attempts = 0
        self.success = 0
        self.skipped = 0
        self.errors = 0
        self.last_genes: List[str] = []
        self.start = time.time()

    def bump_attempt(self, genes: Optional[List[str]] = None):
        self.attempts += 1
        if genes:
            self.last_genes = genes[:4]
        self.report(intermediate=True)

    def bump_success(self):
        self.success += 1

    def bump_skip(self):
        self.skipped += 1

    def bump_error(self):
        self.errors += 1

    def _rate(self) -> float:
        elapsed = max(1e-6, time.time() - self.start)
        return self.success / elapsed

    def report(self, intermediate: bool = False):
        if not DEBUG and intermediate:
            return
        elapsed = time.time() - self.start
        tag = self.name if intermediate else f"{self.name} DONE"
        last = ",".join(self.last_genes) if self.last_genes else "-"
        if DEBUG:
            print(
                f"[{_now_str()}] {tag}: rows={self.success}/{self.target} | "
                f"att={self.attempts} skip={self.skipped} err={self.errors} | "
                f"r/s={self._rate():.2f} | last={last} | {elapsed:.1f}s",
                flush=True,
            )

# =============================================================================
# Helpers (no caching)
# =============================================================================
def get_num(m: Optional[Dict[str, Any]], key: str) -> Optional[float]:
    if not m:
        return None
    v = m.get(key)
    if v is None:
        return None
    try:
        if isinstance(v, float) and math.isnan(v):
            return None
        return float(v)
    except Exception:
        return None

def resolve_by_uniprot(conn, gene_symbol: str, uniprot_id: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    protein_id = None
    pdb_file = None
    row = conn.execute(
        "SELECT protein_id, pdb_file FROM data_proteins WHERE uniprot_id=? LIMIT 1",
        (uniprot_id,),
    ).fetchone()
    if row:
        protein_id, pdb_file = row
    return uniprot_id, protein_id, pdb_file

# =============================================================================
# Metrics extractor (always retrieve fresh; no cache)
# =============================================================================
def metrics_for_gene(conn, gene_symbol: str, uniprot_id: str) -> Optional[Dict[str, Any]]:
    """
    Pull feature summaries for a gene. No caching.
    """
    uniprot_id, protein_id, _ = resolve_by_uniprot(conn, gene_symbol, uniprot_id)
    if not (uniprot_id and protein_id):
        return None

    try:
        plddt_stats, _ = fetch_plddt(conn, protein_id)
        fpocket_stats, _ = fetch_fpocket(conn, protein_id)
        sasa_stats, _, pi_stats, _, _ = fetch_sasa_pi(conn, protein_id)
        disorder_stats, _, morf_stats, _ = fetch_disorder(conn, protein_id, uniprot_id)
        biolip_summary, _ = fetch_biolip2(conn, uniprot_id)
    except Exception as e:
        warnings.warn(f"[prot] feature extraction failed for {gene_symbol}: {e}")
        return None

    def s(d, *keys):
        if not isinstance(d, dict):
            return None
        for k in keys:
            if k in d and d[k] is not None:
                return d[k]
        return None

    # BioLiP2 sites
    biolip_n_sites = None
    if isinstance(biolip_summary, dict):
        n_b = biolip_summary.get("n_binding_sites")
        if isinstance(n_b, (int, float)):
            try:
                biolip_n_sites = int(n_b)
            except Exception:
                biolip_n_sites = None
        if biolip_n_sites is None and "sites" in biolip_summary:
            scont = biolip_summary["sites"]
            if isinstance(scont, dict):
                biolip_n_sites = len(scont)
            elif isinstance(scont, list):
                biolip_n_sites = len(scont)

    return {
        "plddt_mean": s(plddt_stats, "avg", "mean"),
        "fpocket_max": s(fpocket_stats, "max", "max_score"),
        "sasa_avg": s(sasa_stats, "avg", "total_sasa"),
        "pi_mean_surface": s(pi_stats, "avg", "mean_surface_pi"),
        "disorder_pct": s(disorder_stats, "pct_disordered"),
        "morf_max": s(morf_stats, "max"),
        "morf_pct_high": s(morf_stats, "pct_high"),
        "biolip_n_sites": biolip_n_sites,
    }

# =============================================================================
# Evaluators and templates
# =============================================================================
def _idx_unique_extreme(values: List[Optional[float]], mode: str) -> Optional[int]:
    nums = [(i, v) for i, v in enumerate(values) if v is not None]
    if not nums:
        return None
    target = (max if mode == "max" else min)(v for _, v in nums)
    cand = [i for i, v in nums if v == target]
    return cand[0] if len(cand) == 1 else None

def _idx_first_true_unique(mask: List[bool]) -> Optional[int]:
    cand = [i for i, t in enumerate(mask) if t]
    return cand[0] if len(cand) == 1 else None

def _diagnose_failure(metric: str, vals: List[Optional[float]], mask: Optional[List[bool]] = None) -> str:
    finite = sum(v is not None for v in vals)
    if finite == 0:
        return "no_finite_values"
    if mask is not None:
        ntrue = sum(bool(x) for x in mask)
        return f"threshold_unique_fail ntrue={ntrue}"
    uniq = len({v for v in vals if v is not None})
    size = len([v for v in vals if v is not None])
    return "not_unique_extremum" if uniq != size else "internal_eval_fail"

Template = Tuple[str, str, Optional[str], Callable[[List[Optional[Dict[str, Any]]]], Tuple[Optional[int], Dict[str, Any]]]]

def _ev_stat(ms, key, mode):
    vals = [get_num(m, key) for m in ms]
    idx = _idx_unique_extreme(vals, mode)
    diag = {"metric": key, "vals": vals, "reason": None}
    if idx is None:
        diag["reason"] = _diagnose_failure(key, vals)
    return idx, diag

def _ev_threshold_unique(ms, key, thr, op: str):
    vals = [get_num(m, key) for m in ms]
    mask = [(v is not None and (v > thr if op == ">" else v >= thr)) for v in vals]
    idx = _idx_first_true_unique(mask)
    diag = {"metric": key, "vals": vals, "reason": None, "mask": mask, "thr": thr, "op": op}
    if idx is None:
        diag["reason"] = _diagnose_failure(key, vals, mask)
    return idx, diag

# Evaluator with margin requirement (percent) for unique max/min
def _ev_stat_with_margin(ms, key, mode, margin_pct):
    vals = [get_num(m, key) for m in ms]
    idx = _idx_unique_extreme(vals, mode)
    diag = {"metric": key, "vals": vals, "reason": None, "margin_pct": margin_pct, "mode": mode}
    if idx is None:
        diag["reason"] = _diagnose_failure(key, vals)
        return None, diag
    best = vals[idx]
    if best is None:
        diag["reason"] = "no_finite_values"
        return None, diag
    comp_vals = [v for i, v in enumerate(vals) if i != idx and v is not None]
    if not comp_vals:
        diag["reason"] = "no_comparable_values"
        return None, diag
    if mode == "max":
        factor = 1.0 + (margin_pct / 100.0)
        ok = all(best >= factor * v for v in comp_vals)
    else:
        factor = 1.0 - (margin_pct / 100.0)
        if factor < 0.0:
            factor = 0.0
        ok = all(best <= factor * v for v in comp_vals)
    if not ok:
        diag["reason"] = "margin_not_satisfied"
        return None, diag
    return idx, diag

# BioLiP2 evaluator using a parameterized threshold constant
def _ev_biolip_sites(ms, min_sites: int = BIOLIP_SITES_AT_LEAST):
    vals = [get_num(m, "biolip_n_sites") for m in ms]
    mask = [(v is not None and v >= float(min_sites)) for v in vals]
    idx = _idx_first_true_unique(mask)
    reason = None
    if idx is None:
        idx = _idx_unique_extreme(vals, "max")
        if idx is None:
            reason = _diagnose_failure("biolip_n_sites", vals, mask)
    diag = {"metric": "biolip_n_sites", "vals": vals, "mask": mask, "reason": reason, "min_sites": min_sites}
    return idx, diag

PROT_TEMPLATES: List[Template] = [
    ("According to AlphaFold pLDDT, which gene's protein has the highest average pLDDT?",
     "plddt_mean_highest", "plddt_mean", lambda ms: _ev_stat(ms, "plddt_mean", "max")),
    ("According to AlphaFold pLDDT, which gene's protein has an average pLDDT above 70?",
     "plddt_mean_gt_70", "plddt_mean", lambda ms: _ev_threshold_unique(ms, "plddt_mean", PLDDT_MEAN_GT_THRESHOLD, ">")),
    ("According to FPocket, which gene's protein has the most druggable pocket (highest pocket score across all pockets)?",
     "fpocket_max_highest", "fpocket_max", lambda ms: _ev_stat(ms, "fpocket_max", "max")),
    ("According to FreeSASA, which gene's protein has the highest average solvent-accessible surface area (SASA)?",
     "sasa_avg_highest", "sasa_avg", lambda ms: _ev_stat_with_margin(ms, "sasa_avg", "max", SASA_AVG_HIGHEST_PER_CENT)),
    ("According to the surface Polarity Index (PI), which gene's protein has the highest average surface polarity (PI)?",
     "pi_mean_surface_highest", "pi_mean_surface",
     lambda ms: _ev_stat_with_margin(ms, "pi_mean_surface", "max", PI_MEAN_SURFACE_HIGHEST_PER_CENT)),
    ("According to IUPred3, which gene's protein has the lowest percentage of disordered residues?",
     "disorder_pct_lowest", "disorder_pct", lambda ms: _ev_stat_with_margin(ms, "disorder_pct", "min", DISORDER_PCT_LOWEST_PER_CENT)),
    ("According to IUPred3 MoRF propensity, which gene's protein has the highest maximum MoRF score?",
     "morf_max_highest", "morf_max", lambda ms: _ev_stat_with_margin(ms, "morf_max", "max", MORF_MAX_HIGHEST_PER_CENT)),
    (f"According to BioLiP2, which gene's protein has the highest number of experimentally supported binding sites (≥{int(BIOLIP_SITES_AT_LEAST)}) ?",
     "biolip_sites_pref_ge1_else_max", "biolip_n_sites",
     lambda ms: _ev_biolip_sites(ms, BIOLIP_SITES_AT_LEAST)),
]

# =============================================================================
# Sampling
# =============================================================================
class QuartetSampler:
    def __init__(self, pairs: List[Tuple[str, str]], seed: Optional[int] = None):
        self.pairs = pairs[:]  # (gene_symbol, uniprot_id)
        random.Random(seed).shuffle(self.pairs)
        self._recent: List[str] = []

    def sample(self) -> List[Tuple[str, str]]:
        # Avoid very recent repeats; fall back to full pool if needed.
        candidates = [p for p in self.pairs if p[0] not in self._recent[-32:]]
        src = candidates if len(candidates) >= 4 else self.pairs
        pool = random.sample(src, 4)
        self._recent.extend([g for g, _ in pool])
        if len(self._recent) > 256:
            self._recent = self._recent[-256:]
        return pool

# =============================================================================
# Generation
# =============================================================================
def _letters():
    return ["A", "B", "C", "D"]

def _format_option(gene: str) -> str:
    return gene  # no metric text in options

def generate_prot_questions(
    out_path: Optional[Path] = None,
    seed: Optional[int] = SEED,
):
    if seed is not None:
        random.seed(seed)

    # Load TSV with columns: gene_symbol, uniprot_id
    tsv_path = KNOWN_GENE_LIST_PATH
    if not tsv_path.exists():
        raise FileNotFoundError(f"TSV not found: {tsv_path}")

    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    df = df.rename(columns={c: c.strip() for c in df.columns})
    if "gene_symbol" not in df.columns or "uniprot_id" not in df.columns:
        raise ValueError("Input TSV must contain columns: gene_symbol, uniprot_id")

    pairs = [(r["gene_symbol"].strip(), r["uniprot_id"].strip()) for _, r in df.iterrows()
             if str(r["gene_symbol"]).strip() and str(r["uniprot_id"]).strip()]
    if len(pairs) < 4:
        raise ValueError("Need at least 4 (gene_symbol, uniprot_id) pairs in known_tool_prot_genes_list.tsv")

    sampler = QuartetSampler(pairs, seed=seed)
    conn = get_connection()
    rows: List[Dict[str, Any]] = []

    # Build a worklist with each template repeated TEMPLATE_MAX_PER_ID times
    tmpl_ids = [t[1] for t in PROT_TEMPLATES]
    worklist: List[str] = []
    for tid in tmpl_ids:
        worklist.extend([tid] * TEMPLATE_MAX_PER_ID)
    random.shuffle(worklist)

    # Map id -> template tuple
    tmpl_by_id = {tid: t for (txt, tid, key, ev) in PROT_TEMPLATES for t in [(txt, tid, key, ev)]}

    p = Progress("ProtQGen", target=len(worklist))
    if DEBUG:
        print(f"[{_now_str()}] Starting generation: target={len(worklist)} (exactly {TEMPLATE_MAX_PER_ID} per template), "
              f"attempts per item={MAX_ATTEMPTS_PER_Q}, seed={seed}")

    try:
        for idx_wl, template_id in enumerate(worklist, start=1):
            template, _, metric_key, evaluator = tmpl_by_id[template_id]
            success = False
            last_diag: Optional[Dict[str, Any]] = None

            for attempt in range(1, MAX_ATTEMPTS_PER_Q + 1):
                quartet = sampler.sample()
                genes = [g for g, _ in quartet]
                p.bump_attempt(genes)

                # Collect metrics (fresh each time)
                metrics_list: List[Optional[Dict[str, Any]]] = []
                for g, u in quartet:
                    try:
                        metrics_list.append(metrics_for_gene(conn, g, u))
                    except Exception as e:
                        warnings.warn(f"[prot] metrics failed for {g}: {e}")
                        metrics_list.append(None)

                idx, diag = evaluator(metrics_list)
                last_diag = diag

                if idx is not None:
                    letters = _letters()
                    ans_letter = letters[idx]
                    correct_gene = genes[idx]

                    # Build options (no metrics in text)
                    opt_texts = [_format_option(genes[i]) for i in range(4)]
                    q_text = f"{template} [A] {opt_texts[0]} [B] {opt_texts[1]} [C] {opt_texts[2]} [D] {opt_texts[3]}"

                    # Append ONLY the 3 required fields
                    rows.append({
                        "question": q_text,
                        "answer": ans_letter,
                        "tool": "prot",
                    })

                    # ---- Console metrics print for selected questions ----
                    metric_name = diag.get("metric")
                    vals = diag.get("vals", [None, None, None, None])
                    vals_s = "[" + ",".join(("NA" if v is None else f"{float(v):.4g}") for v in vals) + "]"
                    print(
                        f"[{_now_str()}] WL{idx_wl:02d} attempt#{attempt:02d} ok tmpl={template_id} "
                        f"genes={genes} metric={metric_name} vals={vals_s} "
                        f"answer={ans_letter} {correct_gene}",
                        flush=True,
                    )
                    # ----------------------------------------------------

                    p.bump_success()
                    success = True
                    break
                else:
                    if DEBUG:
                        print(f"[{_now_str()}] WL{idx_wl:02d} attempt#{attempt:02d} fail tmpl={template_id} genes={genes}", flush=True)

            if not success:
                p.bump_skip()
                print(f"[{_now_str()}] WL{idx_wl:02d}: SKIPPED after {MAX_ATTEMPTS_PER_Q} attempts (tmpl={template_id})")
                if last_diag is not None:
                    vals = last_diag.get("vals", [])
                    vals_s = "[" + ",".join(("NA" if v is None else f"{float(v):.4g}") for v in vals) + "]"
                    print(f"[{_now_str()}] WL{idx_wl:02d} last diag metric={last_diag.get('metric')} "
                          f"vals={vals_s} reason={last_diag.get('reason')}")

    finally:
        conn.close()

    # ---------- Write ONLY question, answer, tool as TSV ----------
    out_path = out_path or (OUT_DIR / "prot_validation_questions.tsv")
    df_out = pd.DataFrame(rows, columns=["question", "answer", "tool"])
    df_out.to_csv(out_path, sep="\t", index=False, header=True)
    print(f"[{_now_str()}] OK Wrote {len(df_out)} -> {out_path}")
    # --------------------------------------------------------------

    p.report(intermediate=False)

# =============================================================================
# Main
# =============================================================================
if __name__ == "__main__":
    random.seed(SEED)
    generate_prot_questions()
