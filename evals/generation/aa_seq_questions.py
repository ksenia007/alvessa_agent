# alvessa_agent/evals/generation/aa_seq_questions.py
# Author: Dmitri Kosenkov
# Updated: 2025-12-10
#
# Description
# -----------
# Generates multiple-choice validation questions for sequence-to-gene
# resolution using the `aa_seq` tool and the local UniProt DB.
#
# Input (SQLite):
#   local_dbs/uniprot_hs_9606_reviewed.db
#   Table: uniprot_gene_index
#     acc            TEXT PRIMARY KEY
#     canonical_acc  TEXT NOT NULL
#     sequence       TEXT NOT NULL
#     gene_name      TEXT NOT NULL
#     entrez_gene_id TEXT
#
# Output (CSV):
#   alvessa_agent/results/set1.csv
#   alvessa_agent/results/set2.csv
#   alvessa_agent/results/set3.csv
#   alvessa_agent/results/set4.csv
#
# Each CSV has exactly 20 questions from a single template.
# Columns: question, answer, tool, tool_specified_in_question
#

import sys
import time
import random
import sqlite3
from pathlib import Path
from dataclasses import dataclass
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

from src.tools.aa_seq import UNIPROT_DB_PATH, validate_uniprot_db

# =============================================================================
# Config
# =============================================================================
SEED = 42
QUESTIONS_PER_SET = 20  # exactly 20 per set
SEQ_FRAGMENT_MAX_LEN = 120
APPROX_MUTATIONS_MIN = 3
APPROX_MUTATIONS_MAX = 6
MIN_SEQUENCE_LENGTH = 80  # minimum full sequence length to sample from

RAND = random.Random(SEED)

# Keywords for detecting when a question explicitly specifies a tool or DB
_TOOL_NAME_KEYWORDS = [
    "aa_seq",
    "sequence-to-gene resolver",
    "sequence resolver",
    "sequence search",
    "uniprot",
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
# DB helpers
# =============================================================================
def _get_uniprot_conn() -> sqlite3.Connection:
    """
    Open validated UniProt DB connection.
    """
    validate_uniprot_db(UNIPROT_DB_PATH)
    if not UNIPROT_DB_PATH.exists():
        raise FileNotFoundError(f"UniProt DB not found: {UNIPROT_DB_PATH}")
    conn = sqlite3.connect(str(UNIPROT_DB_PATH))
    conn.row_factory = sqlite3.Row
    return conn


@dataclass
class SeqEntry:
    acc: str
    canonical_acc: str
    sequence: str
    gene_name: str


def _load_seq_entries(conn: sqlite3.Connection, min_len: int = MIN_SEQUENCE_LENGTH) -> List[SeqEntry]:
    """
    Load candidate entries from uniprot_gene_index with a minimum sequence length.
    """
    cur = conn.cursor()
    cur.execute(
        """
        SELECT acc, canonical_acc, sequence, gene_name
        FROM uniprot_gene_index
        WHERE sequence IS NOT NULL
          AND gene_name IS NOT NULL
          AND gene_name <> ''
          AND LENGTH(sequence) >= ?
        """,
        (min_len,),
    )
    rows = cur.fetchall()
    entries: List[SeqEntry] = []
    for r in rows:
        seq = str(r["sequence"]).strip().upper()
        gene = str(r["gene_name"]).strip()
        acc = str(r["acc"]).strip()
        canon = str(r["canonical_acc"]).strip()
        if not seq or not gene or not acc:
            continue
        entries.append(
            SeqEntry(
                acc=acc,
                canonical_acc=canon,
                sequence=seq,
                gene_name=gene,
            )
        )
    if not entries:
        raise RuntimeError("No suitable entries found in uniprot_gene_index.")
    return entries

# =============================================================================
# Sampling
# =============================================================================
class EntrySampler:
    def __init__(self, entries: List[SeqEntry], seed: Optional[int] = None):
        self.entries = entries[:]
        rand_local = random.Random(seed)
        rand_local.shuffle(self.entries)
        self._recent_genes: List[str] = []

    def sample_four(self, distinct_on_gene: bool = False) -> List[SeqEntry]:
        """
        Sample four entries. If distinct_on_gene is True, enforce
        distinct gene_name among the four options.
        """
        for _ in range(64):
            pool = RAND.sample(self.entries, 4)
            if not distinct_on_gene:
                self._recent_genes.extend([e.gene_name for e in pool])
                if len(self._recent_genes) > 256:
                    self._recent_genes = self._recent_genes[-256:]
                return pool
            genes = {e.gene_name for e in pool}
            if len(genes) == 4:
                self._recent_genes.extend([e.gene_name for e in pool])
                if len(self._recent_genes) > 256:
                    self._recent_genes = self._recent_genes[-256:]
                return pool
        # Fallback: return any four
        pool = RAND.sample(self.entries, 4)
        self._recent_genes.extend([e.gene_name for e in pool])
        if len(self._recent_genes) > 256:
            self._recent_genes = self._recent_genes[-256:]
        return pool

# =============================================================================
# Sequence utilities
# =============================================================================
_AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")


def _extract_fragment(seq: str, max_len: int = SEQ_FRAGMENT_MAX_LEN) -> str:
    """
    Extract a contiguous fragment (prefix or random window) up to max_len.
    """
    seq = seq.strip().upper()
    if len(seq) <= max_len:
        return seq
    start = RAND.randint(0, len(seq) - max_len)
    return seq[start : start + max_len]


def _mutate_sequence(seq: str) -> str:
    """
    Introduce a small number of point mutations into seq.
    """
    if not seq:
        return seq
    n = len(seq)
    if n == 1:
        choices = [aa for aa in _AMINO_ACIDS if aa != seq]
        return RAND.choice(choices) if choices else seq

    max_allowed = max(1, int(0.2 * n))
    n_mut = min(APPROX_MUTATIONS_MAX, max_allowed)
    n_mut = max(APPROX_MUTATIONS_MIN, n_mut)
    n_mut = min(n_mut, n)

    positions = RAND.sample(range(n), n_mut)
    seq_list = list(seq)
    for pos in positions:
        current = seq_list[pos]
        choices = [aa for aa in _AMINO_ACIDS if aa != current]
        if not choices:
            continue
        seq_list[pos] = RAND.choice(choices)
    mutated = "".join(seq_list)
    if mutated == seq:
        return _mutate_sequence(seq)
    return mutated

# =============================================================================
# Tool name detection
# =============================================================================
def _tool_specified_in_question(q_text: str) -> str:
    """
    Return "True" if the question explicitly names a tool or database
    (aa_seq, UniProt, sequence resolver, protein tool, etc.), else "False".
    Comparison is case-insensitive.
    """
    q_lower = q_text.lower()
    for kw in _TOOL_NAME_KEYWORDS:
        if kw.lower() in q_lower:
            return "True"
    return "False"

# =============================================================================
# Templates
# =============================================================================
# Each template:
#   (template_id, question_stem, is_approximate, answer_mode)
# answer_mode: "gene" or "acc"
SeqTemplate = Tuple[str, str, bool, str]

SEQ_TEMPLATES: List[SeqTemplate] = [
    (
        "exact_gene",
        "Given the amino acid sequence fragment below, which gene does it best match exactly in UniProt?",
        False,
        "gene",
    ),
    (
        "approx_gene",
        "Given the amino acid sequence fragment below (with a few residues mutated), which gene does it best match in UniProt?",
        True,
        "gene",
    ),
]


def _letters() -> List[str]:
    return ["A", "B", "C", "D"]


def _fmt_options_gene(entries: List[SeqEntry]) -> List[str]:
    return [e.gene_name for e in entries]


def _fmt_options_acc(entries: List[SeqEntry]) -> List[str]:
    return [e.acc for e in entries]

# =============================================================================
# Per-template generation
# =============================================================================
def _generate_set_for_template(
    sampler: EntrySampler,
    template: SeqTemplate,
    set_index: int,
    out_dir: Path,
) -> Path:
    """
    Generate exactly QUESTIONS_PER_SET questions for a single template and write to CSV.
    Returns the output path.
    """
    template_id, stem, is_approx, answer_mode = template
    rows: List[Dict[str, Any]] = []

    p = Progress(f"SeqQGen_set{set_index}_{template_id}", target=QUESTIONS_PER_SET)
    if DEBUG:
        print(
            f"[{_now_str()}] Starting set {set_index} for template={template_id}: "
            f"target={QUESTIONS_PER_SET}",
            flush=True,
        )

    for q_idx in range(1, QUESTIONS_PER_SET + 1):
        distinct_on_gene = (answer_mode == "gene")
        entries = sampler.sample_four(distinct_on_gene=distinct_on_gene)
        genes = [e.gene_name for e in entries]
        p.bump_attempt(genes)

        correct_idx = RAND.randint(0, 3)
        target = entries[correct_idx]

        base_fragment = _extract_fragment(target.sequence, max_len=SEQ_FRAGMENT_MAX_LEN)
        seq_for_question = _mutate_sequence(base_fragment) if is_approx else base_fragment

        if answer_mode == "gene":
            opt_texts = _fmt_options_gene(entries)
        else:
            opt_texts = _fmt_options_acc(entries)

        # Question text includes AA sequence and options, formatted like other tools.
        q_text = (
            f"{stem}\n\n"
            f"Sequence fragment:\n"
            f"{seq_for_question}\n\n"
            f"[A] {opt_texts[0]} [B] {opt_texts[1]} "
            f"[C] {opt_texts[2]} [D] {opt_texts[3]}"
        )

        letters = _letters()
        ans_letter = letters[correct_idx]
        tool_flag = _tool_specified_in_question(q_text)

        rows.append(
            {
                "question": q_text,
                "answer": ans_letter,
                "tool": "aa_seq",
                "tool_specified_in_question": tool_flag,
            }
        )

        if DEBUG:
            print(
                f"[{_now_str()}] SET{set_index:02d} Q{q_idx:02d} "
                f"ok tmpl={template_id} correct_idx={correct_idx} "
                f"gene={target.gene_name} acc={target.acc}",
                flush=True,
            )

        p.bump_success()

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
def generate_aa_seq_questions(
    out_path: Optional[Path] = None,
    seed: Optional[int] = SEED,
) -> List[Path]:
    """
    Generate four question sets for AA-sequence-based resolution.

    Set 1: exact sequence fragment -> gene_name
    Set 2: exact sequence fragment -> UniProt accession
    Set 3: approximate sequence fragment (mutated) -> gene_name
    Set 4: approximate sequence fragment (mutated) -> UniProt accession

    Each set:
      - Contains exactly QUESTIONS_PER_SET questions.
      - Uses a single question stem from SEQ_TEMPLATES.
      - Is written to setN.csv (N starts at 1) in the target directory.

    Returns:
      List of Paths to the generated CSV files.
    """
    if seed is not None:
        random.seed(seed)
        RAND.seed(seed)

    # Decide output directory based on out_path
    if out_path is not None and out_path.is_dir():
        out_dir = out_path
    elif out_path is not None:
        out_dir = out_path.parent
    else:
        out_dir = OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    conn = _get_uniprot_conn()
    try:
        entries = _load_seq_entries(conn, min_len=MIN_SEQUENCE_LENGTH)
    finally:
        try:
            conn.close()
        except Exception:
            pass

    sampler = EntrySampler(entries, seed=seed)
    output_files: List[Path] = []

    for set_index, template in enumerate(SEQ_TEMPLATES, start=1):
        out_file = _generate_set_for_template(
            sampler=sampler,
            template=template,
            set_index=set_index,
            out_dir=out_dir,
        )
        output_files.append(out_file)

    return output_files

# =============================================================================
# Main
# =============================================================================
if __name__ == "__main__":
    random.seed(SEED)
    RAND.seed(SEED)
    generate_aa_seq_questions()
