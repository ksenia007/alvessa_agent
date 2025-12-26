# alvessa_agent/evals/generation/chembl_questions.py
# Author: Dmitri Kosenkov
# Updated: 2025-12-10
#
# Description
# -----------
# Generates multiple-choice validation questions for the `chembl` tool
# using local ChEMBL (schema v35) evidence.
#
# Input (TSV, tab-separated):
#   alvessa_agent/local_dbs/known_tool_prot_genes_list.tsv
#   Columns: gene_symbol  uniprot_id
#
# Output (CSV):
#   alvessa_agent/results/set1.csv
#   alvessa_agent/results/set2.csv
#   ...
#   Each CSV has exactly 20 questions from a single template.
#   Columns: question, answer, tool, tool_specified_in_question
#

import sys
import time
import math
import json
import random
import sqlite3
import warnings
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple

import pandas as pd

# =============================================================================
# Bootstrap
# =============================================================================
ROOT = Path(__file__).resolve().parents[2]  # alvessa_agent/
SRC = ROOT / "src"
OUT_DIR = ROOT / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))
if str(SRC) not in sys.path:
    sys.path.append(str(SRC))

DEBUG = False

# Hard-pinned input TSV path
KNOWN_GENE_LIST_PATH = ROOT / "local_dbs" / "known_tool_prot_genes_list.tsv"

# Direct ChEMBL DB path
CHEMBL_DB_PATH = ROOT / "local_dbs" / "chembl_35.db"

# =============================================================================
# Config
# =============================================================================
SEED = 42
MAX_ATTEMPTS_PER_Q = 8000
QUESTIONS_PER_SET = 20  # exactly 20 per set

ASSAY_TYPES = ("IC50", "Ki", "EC50")
ASSAY_MARGIN_PERCENT = 30.0  # require winner to be at least 30% better (lower) than all others

RAND = random.Random(SEED)

# Keywords for detecting when a question explicitly specifies a tool or DB
_TOOL_NAME_KEYWORDS = [
    "chembl",
    "chembl tool",
    "drugcentral",
    "alphamissense",
    "clinvar",
    "cysdb",
    "protein tool",
]

# =============================================================================
# Logging / Progress
# =============================================================================
def _now_str() -> str:
    import datetime as _dt
    return _dt.datetime.now().strftime("%H:%M:%S")


class Progress:
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
# DB helpers (direct SQLite)
# =============================================================================
def get_conn() -> sqlite3.Connection:
    if not CHEMBL_DB_PATH.exists():
        raise FileNotFoundError(f"ChEMBL DB not found: {CHEMBL_DB_PATH}")
    conn = sqlite3.connect(str(CHEMBL_DB_PATH))
    conn.row_factory = sqlite3.Row
    return conn

def resolve_tid(conn: sqlite3.Connection, uniprot_id: str) -> Optional[int]:
    cur = conn.cursor()
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
    return int(row["tid"]) if row else None

def _to_nM(value: float, unit: str) -> Optional[float]:
    if value is None or value <= 0:
        return None
    u = (unit or "").lower()
    if u == "nm":
        return value
    if u == "um":
        return value * 1000.0
    if u == "pm":
        return value / 1000.0
    if u == "m":
        return value * 1e9
    return None

# Caches within a run (memory only)
_tid_cache: Dict[str, Optional[int]] = {}
_mech_cache: Dict[int, List[str]] = {}
_ind_cache: Dict[int, List[str]] = {}
_bb_cache: Dict[int, bool] = {}
_ph3_cache: Dict[int, bool] = {}
_assay_cache: Dict[int, Dict[str, Optional[float]]] = {}

def get_tid_cached(conn: sqlite3.Connection, uniprot_id: str) -> Optional[int]:
    if uniprot_id not in _tid_cache:
        _tid_cache[uniprot_id] = resolve_tid(conn, uniprot_id)
    return _tid_cache[uniprot_id]

def fetch_approved_mechanisms(conn: sqlite3.Connection, tid: int) -> List[str]:
    if tid in _mech_cache:
        return _mech_cache[tid]
    cur = conn.cursor()
    cur.execute(
        """
        SELECT DISTINCT dm.mechanism_of_action AS moa
        FROM drug_mechanism dm
        JOIN molecule_dictionary m ON dm.molregno = m.molregno
        WHERE dm.tid = ? AND m.max_phase = 4 AND dm.mechanism_of_action IS NOT NULL
        """,
        (tid,),
    )
    mechs = [str(r["moa"]).strip() for r in cur.fetchall() if str(r["moa"]).strip()]
    _mech_cache[tid] = mechs
    return mechs

def fetch_approved_indications(conn: sqlite3.Connection, tid: int) -> List[str]:
    if tid in _ind_cache:
        return _ind_cache[tid]
    cur = conn.cursor()
    cur.execute(
        """
        SELECT DISTINCT COALESCE(di.efo_term, di.mesh_heading) AS ind
        FROM drug_mechanism dm
        JOIN molecule_dictionary m ON dm.molregno = m.molregno
        JOIN drug_indication di ON di.molregno = m.molregno
        WHERE dm.tid = ? AND m.max_phase = 4
        """,
        (tid,),
    )
    inds = [str(r["ind"]).strip() for r in cur.fetchall() if r["ind"] and str(r["ind"]).strip()]
    _ind_cache[tid] = inds
    return inds

def has_black_box_approved(conn: sqlite3.Connection, tid: int) -> bool:
    if tid in _bb_cache:
        return _bb_cache[tid]
    cur = conn.cursor()
    cur.execute(
        """
        SELECT 1
        FROM drug_mechanism dm
        JOIN molecule_dictionary m ON dm.molregno = m.molregno
        WHERE dm.tid = ? AND m.max_phase = 4 AND m.black_box_warning = 1
        LIMIT 1
        """,
        (tid,),
    )
    _bb_cache[tid] = cur.fetchone() is not None
    return _bb_cache[tid]

def has_phase3(conn: sqlite3.Connection, tid: int) -> bool:
    if tid in _ph3_cache:
        return _ph3_cache[tid]
    cur = conn.cursor()
    cur.execute(
        """
        SELECT 1
        FROM drug_mechanism dm
        JOIN molecule_dictionary m ON dm.molregno = m.molregno
        WHERE dm.tid = ? AND m.max_phase = 3
        LIMIT 1
        """,
        (tid,),
    )
    _ph3_cache[tid] = cur.fetchone() is not None
    return _ph3_cache[tid]

def best_assay_values_by_type(conn: sqlite3.Connection, tid: int) -> Dict[str, Optional[float]]:
    if tid in _assay_cache:
        return _assay_cache[tid]
    cur = conn.cursor()
    cur.execute(
        """
        SELECT a.standard_type AS stype, a.standard_value AS sval, a.standard_units AS sunits
        FROM activities a
        JOIN assays ass ON a.assay_id = ass.assay_id
        WHERE ass.tid = ?
          AND a.standard_type IN ('IC50','EC50','Ki')
          AND a.standard_value > 0
          AND a.standard_units IN ('nM','uM','pM','M')
        """,
        (tid,),
    )
    buckets: Dict[str, List[float]] = {"IC50": [], "Ki": [], "EC50": []}
    for r in cur.fetchall():
        stype = r["stype"]
        sval = r["sval"]
        sunits = r["sunits"]
        nm = _to_nM(sval, sunits) if sval is not None else None
        if nm is not None and stype in buckets:
            buckets[stype].append(nm)
    best = {k: (min(v) if v else None) for k, v in buckets.items()}
    _assay_cache[tid] = best
    return best

# =============================================================================
# Templates and evaluators
# =============================================================================
TemplateId = str
TemplateDef = Tuple[TemplateId, str]  # (id, question_text_with_placeholders)

TEMPLATES: List[TemplateDef] = [
    ("q1_moa", "Which gene has an FDA-approved drug with the mechanism '{moa}'?"),
    ("q2_black_box", "Which gene has at least one FDA-approved drug with a black box warning?"),
    ("q3_phase3", "Which gene has a compound in phase 3 clinical trials?"),
    ("q4_indication", "Which gene has an FDA-approved drug for {indication}?"),
    ("q5_assay", "Which gene shows stronger binding (lower {assay}) among the options?"),
]

def _letters() -> List[str]:
    return ["A", "B", "C", "D"]

def _fmt_options(genes: List[str]) -> str:
    return f" [A] {genes[0]} [B] {genes[1]} [C] {genes[2]} [D] {genes[3]}"

def _unique_token_owner(token_sets: List[set], token: str) -> Optional[int]:
    owners = [i for i, s in enumerate(token_sets) if token in s]
    return owners[0] if len(owners) == 1 else None

def _pick_unique_token(token_sets: List[set]) -> Tuple[Optional[int], Optional[str]]:
    all_tokens: List[str] = []
    for s in token_sets:
        all_tokens.extend(list(s))
    for t in all_tokens:
        idx = _unique_token_owner(token_sets, t)
        if idx is not None:
            return idx, t
    RAND.shuffle(all_tokens)
    for t in all_tokens:
        idx = _unique_token_owner(token_sets, t)
        if idx is not None:
            return idx, t
    return None, None

def _idx_unique_min_with_margin(vals: List[Optional[float]], margin_pct: float) -> Optional[int]:
    finite = [(i, v) for i, v in enumerate(vals) if v is not None]
    if not finite:
        return None
    best_val = min(v for _, v in finite)
    winners = [i for i, v in finite if v == best_val]
    if len(winners) != 1:
        return None
    idx = winners[0]
    others = [v for j, v in enumerate(vals) if j != idx and v is not None]
    if not others:
        return None
    factor = 1.0 - (margin_pct / 100.0)
    if factor < 0.0:
        factor = 0.0
    ok = all(best_val <= factor * v for v in others)
    return idx if ok else None

# =============================================================================
# Sampler
# =============================================================================
class QuartetSampler:
    def __init__(self, pairs: List[Tuple[str, str]], seed: Optional[int] = None):
        self.pairs = pairs[:]
        random.Random(seed).shuffle(self.pairs)
        self._recent: List[str] = []

    def sample(self) -> List[Tuple[str, str]]:
        candidates = [p for p in self.pairs if p[0] not in self._recent[-32:]]
        src = candidates if len(candidates) >= 4 else self.pairs
        pool = random.sample(src, 4)
        self._recent.extend([g for g, _ in pool])
        if len(self._recent) > 256:
            self._recent = self._recent[-256:]
        return pool

# =============================================================================
# Tool name detection
# =============================================================================
def _tool_specified_in_question(q_text: str) -> str:
    """
    Return "True" if the question explicitly names a tool or database
    (AlphaMissense, ClinVar, DrugCentral, CysDB, ChEMBL, protein tool, etc.), else "False".
    Comparison is case-insensitive.
    """
    q_lower = q_text.lower()
    for kw in _TOOL_NAME_KEYWORDS:
        if kw.lower() in q_lower:
            return "True"
    return "False"

# =============================================================================
# Per-template generation
# =============================================================================
def _generate_set_for_template(
    conn: sqlite3.Connection,
    sampler: QuartetSampler,
    template: TemplateDef,
    set_index: int,
    out_dir: Path,
    seed: Optional[int] = SEED,
) -> Path:
    """
    Generate exactly QUESTIONS_PER_SET questions for a single template and write to CSV.
    Returns the output path.
    """
    template_id, qtext_base = template
    rows: List[Dict[str, Any]] = []

    p = Progress(f"ChemblQGen_set{set_index}_{template_id}", target=QUESTIONS_PER_SET)
    if DEBUG:
        print(
            f"[{_now_str()}] Starting set {set_index} for template={template_id}: "
            f"target={QUESTIONS_PER_SET}, attempts per item={MAX_ATTEMPTS_PER_Q}, seed={seed}",
            flush=True,
        )

    for q_idx in range(1, QUESTIONS_PER_SET + 1):
        success = False

        for attempt in range(1, MAX_ATTEMPTS_PER_Q + 1):
            quartet = sampler.sample()
            genes = [g for g, _ in quartet]
            unis = [u for _, u in quartet]
            p.bump_attempt(genes)

            # Resolve TIDs
            tids: List[Optional[int]] = []
            for u in unis:
                try:
                    tids.append(get_tid_cached(conn, u))
                except Exception:
                    tids.append(None)
            if any(t is None for t in tids):
                continue

            ans_idx: Optional[int] = None
            rendered: Optional[str] = None
            diag_info: Dict[str, Any] = {}

            try:
                if template_id == "q1_moa":
                    s_mechs = [set(fetch_approved_mechanisms(conn, int(t))) for t in tids]  # type: ignore[arg-type]
                    idx, moa = _pick_unique_token(s_mechs)
                    if idx is not None and moa:
                        rendered = f"{qtext_base.format(moa=moa)}{_fmt_options(genes)}"
                        ans_idx = idx
                        diag_info = {"metric": "moa_unique", "moa": moa}

                elif template_id == "q2_black_box":
                    flags = [has_black_box_approved(conn, int(t)) for t in tids]  # type: ignore[arg-type]
                    idxs = [i for i, v in enumerate(flags) if v]
                    if len(idxs) == 1:
                        rendered = f"{qtext_base}{_fmt_options(genes)}"
                        ans_idx = idxs[0]
                        diag_info = {"metric": "black_box_approved", "mask": flags}

                elif template_id == "q3_phase3":
                    flags = [has_phase3(conn, int(t)) for t in tids]  # type: ignore[arg-type]
                    idxs = [i for i, v in enumerate(flags) if v]
                    if len(idxs) == 1:
                        rendered = f"{qtext_base}{_fmt_options(genes)}"
                        ans_idx = idxs[0]
                        diag_info = {"metric": "phase3_presence", "mask": flags}

                elif template_id == "q4_indication":
                    s_inds = [set(fetch_approved_indications(conn, int(t))) for t in tids]  # type: ignore[arg-type]
                    idx, ind = _pick_unique_token(s_inds)
                    if idx is not None and ind:
                        rendered = f"{qtext_base.format(indication=ind)}{_fmt_options(genes)}"
                        ans_idx = idx
                        diag_info = {"metric": "indication_unique", "indication": ind}

                elif template_id == "q5_assay":
                    by_gene = [best_assay_values_by_type(conn, int(t)) for t in tids]  # type: ignore[arg-type]
                    chosen: Optional[Tuple[str, int, List[Optional[float]]]] = None
                    for assay in ASSAY_TYPES:
                        vals = [gd.get(assay) for gd in by_gene]
                        if any(v is None for v in vals):
                            continue
                        idx = _idx_unique_min_with_margin(vals, ASSAY_MARGIN_PERCENT)
                        if idx is not None:
                            chosen = (assay, idx, vals)
                            break
                    if chosen:
                        assay, idx, vals = chosen
                        rendered = f"{qtext_base.format(assay=assay)}{_fmt_options(genes)}"
                        ans_idx = idx
                        diag_info = {"metric": f"{assay}_min_nM_with_margin", "vals": vals}
            except Exception as e:
                warnings.warn(f"[chembl] evaluation failed: {e}")
                rendered = None
                ans_idx = None

            if rendered is not None and ans_idx is not None:
                letters = _letters()
                ans_letter = letters[ans_idx]
                tool_flag = _tool_specified_in_question(rendered)

                rows.append(
                    {
                        "question": rendered,
                        "answer": ans_letter,
                        "tool": "chembl",
                        "tool_specified_in_question": tool_flag,
                    }
                )

                metric = diag_info.get("metric")
                if "vals" in diag_info:
                    vals = diag_info["vals"]
                    vals_s = "[" + ",".join(("NA" if v is None else f"{float(v):.4g}") for v in vals) + "]"
                    extra = f" vals={vals_s}"
                elif "mask" in diag_info:
                    extra = f" mask={diag_info['mask']}"
                elif "moa" in diag_info:
                    extra = f" moa={diag_info['moa']}"
                elif "indication" in diag_info:
                    extra = f" indication={diag_info['indication']}"
                else:
                    extra = ""
                print(
                    f"[{_now_str()}] SET{set_index:02d} Q{q_idx:02d} attempt#{attempt:02d} "
                    f"ok tmpl={template_id} genes={genes} metric={metric}{extra} "
                    f"answer={ans_letter} {genes[ans_idx]}",
                    flush=True,
                )

                p.bump_success()
                success = True
                break  # next question in this set

        if not success:
            p.bump_skip()
            print(
                f"[{_now_str()}] SET{set_index:02d} Q{q_idx:02d}: "
                f"SKIPPED after {MAX_ATTEMPTS_PER_Q} attempts (tmpl={template_id})"
            )
            # Hard requirement: 20 questions per set. Fail fast if we cannot generate.
            raise RuntimeError(
                f"Failed to generate question {q_idx} for template {template_id} "
                f"in set {set_index} after {MAX_ATTEMPTS_PER_Q} attempts."
            )

    p.report(intermediate=False)

    out_path = out_dir / f"set{set_index}.csv"
    df_out = pd.DataFrame(
        rows,
        columns=["question", "answer", "tool", "tool_specified_in_question"],
    )
    df_out.to_csv(out_path, index=False, header=True)
    print(f"[{_now_str()}] OK Wrote {len(df_out)} questions -> {out_path}")
    return out_path

# =============================================================================
# Generation core
# =============================================================================
def generate_chembl_questions(
    out_path: Optional[Path] = None,
    seed: Optional[int] = SEED,
) -> List[Path]:
    """
    Generate question sets for all TEMPLATES.

    Each set:
      - Contains exactly QUESTIONS_PER_SET questions.
      - Uses a single question stem from TEMPLATES.
      - Is written to setN.csv (N starts at 1) in the target directory.

    Returns:
      List of Paths to the generated CSV files.
    """
    if seed is not None:
        random.seed(seed)

    # Decide output directory based on out_path
    if out_path is not None and out_path.is_dir():
        out_dir = out_path
    elif out_path is not None:
        out_dir = out_path.parent
    else:
        out_dir = OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load (gene_symbol, uniprot_id)
    if not KNOWN_GENE_LIST_PATH.exists():
        raise FileNotFoundError(f"TSV not found: {KNOWN_GENE_LIST_PATH}")
    df = pd.read_csv(KNOWN_GENE_LIST_PATH, sep="\t", dtype=str)
    df = df.rename(columns={c: c.strip() for c in df.columns})
    if "gene_symbol" not in df.columns or "uniprot_id" not in df.columns:
        raise ValueError("Input TSV must contain columns: gene_symbol, uniprot_id")

    pairs = [
        (r["gene_symbol"].strip(), r["uniprot_id"].strip())
        for _, r in df.iterrows()
        if str(r["gene_symbol"]).strip() and str(r["uniprot_id"]).strip()
    ]
    if len(pairs) < 4:
        raise ValueError(
            "Need at least 4 (gene_symbol, uniprot_id) pairs in known_tool_prot_genes_list.tsv"
        )

    sampler = QuartetSampler(pairs, seed=seed)
    conn = get_conn()
    output_files: List[Path] = []

    try:
        for set_index, template in enumerate(TEMPLATES, start=1):
            out_file = _generate_set_for_template(
                conn=conn,
                sampler=sampler,
                template=template,
                set_index=set_index,
                out_dir=out_dir,
                seed=seed,
            )
            output_files.append(out_file)
    finally:
        try:
            conn.close()
        except Exception:
            pass

    return output_files

# =============================================================================
# Main
# =============================================================================
if __name__ == "__main__":
    random.seed(SEED)
    generate_chembl_questions()
