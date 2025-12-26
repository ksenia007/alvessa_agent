# src/tools/aa_seq/seq_search.py
# ===========================================================
#  Sequence-to-Gene Resolver
# ===========================================================
#
# This module resolves amino acid sequences to gene-associated UniProt
# entries using a local SQLite UniProt database with:
#   - uniprot            (substring search)
#   - uniprot_gene_index (k-mer approximate search)
#   - kmer_index         (k-mer approximate search)
#
# Search strategy (per query sequence):
#   1) Substring-based search in table "uniprot" using a local, gapless
#      similarity score (percent identity over a window of length equal
#      to the query).
#   2) K-mer approximate search over "kmer_index" joined with
#      "uniprot_gene_index" (shared k-mers based similarity).
#   3) Merge all hits from both stages by (gene_name, acc), keeping the
#      highest scoring Hit for each (gene_name, acc) pair.
#
# Similarity is always:
#       score = similarity_percent in range 0..100
#
# Public API:
#
#   from src.tools.aa_seq.seq_search import resolve_sequences_to_gene_records
#
#   result = resolve_sequences_to_gene_records(
#       aa_sequences,
#       top_n=5,          # max hits per input sequence
#       min_score=80.0,   # similarity threshold in percent
#       debug=DEBUG,
#   )
#
# Output:
#   Dict[str, Any]
#   {
#     "genes": ["ABC", "DEF", ...],  # unique gene symbols (strings)
#     "records": [
#       {
#         "sequence": "<original_sequence>",
#         "gene_name": ...,
#         "entrez_gene_id": ...,
#         "acc": ...,
#         "canonical_acc": ...,
#         "score": ...,
#         "coverage": ...,
#         "alignment_length": ...,
#         "uniprot_sequence": "<full_uniprot_sequence_for_acc>",
#       },
#       ...
#     ]
#   }

from __future__ import annotations

import sys
import json
import sqlite3
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Any, Dict, List, Optional, Tuple

# Make sure repo root is on sys.path when running this file directly
PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.tools.aa_seq import UNIPROT_DB_PATH, validate_uniprot_db

# --------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------

# k-mer size used to build kmer_index (must match DB builder)
KMER_K: int = 5

# Maximum number of raw rows returned from SQL before filtering
DEFAULT_SQL_LIMIT: int = 200


# --------------------------------------------------------------------
# DATA STRUCTURES
# --------------------------------------------------------------------


@dataclass
class Hit:
    acc: str                  # full accession, e.g. P00533 or P00533-2
    canonical_acc: str        # root accession, e.g. P00533
    gene_name: Optional[str]  # primary gene symbol (can be None)
    entrez_gene_id: Optional[str]
    score: float              # similarity in percent (0..100)
    coverage: float           # fraction of query covered (0..1), approx
    alignment_length: int     # approximate aligned length
    uniprot_sequence: Optional[str] = None  # full UniProt AA sequence for this acc


# --------------------------------------------------------------------
# INTERNAL HELPERS: DB RESOLUTION
# --------------------------------------------------------------------


def _resolve_db_path(db_path: Path | str | None) -> Path:
    """
    Resolve DB path, defaulting to UNIPROT_DB_PATH and validating it.
    """
    raw_path = db_path if db_path is not None else UNIPROT_DB_PATH
    path = raw_path if isinstance(raw_path, Path) else Path(raw_path)
    validate_uniprot_db(path)
    return path


# --------------------------------------------------------------------
# INTERNAL HELPERS: SUBSTRING SEARCH (TABLE "uniprot")
# --------------------------------------------------------------------


def _best_substring_similarity(
    full_seq: str,
    query: str,
) -> Tuple[float, int, int]:
    """
    Compute the best local, gapless similarity between `query` and any
    same-length window in `full_seq`.

    Returns:
        similarity_percent (float): 0..100
        best_start (int): 0-based index where best match starts (-1 if none)
        best_matches (int): number of identical residues in best window

    Similarity definition:
        similarity_percent = 100 * matches / len(query)
    """
    full_seq = full_seq.strip().upper()
    query = query.strip().upper()

    qlen = len(query)
    slen = len(full_seq)

    if qlen == 0 or slen < qlen:
        return 0.0, -1, 0

    best_matches = -1
    best_start = -1

    for start in range(slen - qlen + 1):
        window = full_seq[start : start + qlen]
        matches = sum(1 for a, b in zip(window, query) if a == b)
        if matches > best_matches:
            best_matches = matches
            best_start = start

    if best_matches <= 0:
        return 0.0, -1, 0

    similarity = 100.0 * best_matches / float(qlen)
    return similarity, best_start, best_matches


def _substring_search_hits(
    query_seq: str,
    db_path: Path | str | None = None,
    min_similarity_percent: float = 0.0,
    debug: bool = False,
) -> List[Hit]:
    """
    Run substring-based search in table "uniprot" and compute a local
    similarity score for each hit.

    Returns:
        List[Hit]

    Null gene_name values are allowed and kept.
    """
    query_seq = query_seq.strip().upper()
    if not query_seq:
        if debug:
            print("[UniProt seq_search] Substring search: empty query; no hits.")
        return []

    path = _resolve_db_path(db_path)
    conn = sqlite3.connect(str(path))
    try:
        cur = conn.cursor()

        sql = """
        SELECT
            acc,
            canonical_acc,
            gene_name,
            entrez_gene_id,
            sequence
        FROM uniprot
        WHERE sequence LIKE '%' || ? || '%';
        """

        cur.execute(sql, (query_seq,))
        rows = cur.fetchall()
    finally:
        conn.close()

    if debug:
        print(
            "[UniProt seq_search] Substring search raw SQL hits (LIKE): "
            f"{len(rows)}"
        )

    hits: List[Hit] = []

    for acc, canonical_acc, gene_name, entrez_gene_id, sequence in rows:
        similarity, best_start, best_matches = _best_substring_similarity(
            full_seq=sequence,
            query=query_seq,
        )

        if similarity < min_similarity_percent:
            if debug:
                print(
                    "[UniProt seq_search] Substring hit below threshold: "
                    f"acc={acc}, similarity={similarity:.2f} "
                    f"< min_similarity={min_similarity_percent:.2f}"
                )
            continue

        query_len = len(query_seq)
        if best_start < 0 or query_len == 0:
            coverage = 0.0
            align_len = 0
        else:
            coverage = 1.0
            align_len = query_len

        h = Hit(
            acc=str(acc),
            canonical_acc=str(canonical_acc),
            gene_name=gene_name,
            entrez_gene_id=entrez_gene_id,
            score=similarity,
            coverage=coverage,
            alignment_length=align_len,
            uniprot_sequence=sequence,
        )
        hits.append(h)

    if debug and hits:
        sims = [h.score for h in hits]
        print(
            "[UniProt seq_search] Substring hits kept after similarity filter: "
            f"{len(hits)}, score range min={min(sims):.2f}, max={max(sims):.2f}"
        )

    return hits


# --------------------------------------------------------------------
# INTERNAL HELPERS: K-MER APPROX SEARCH (TABLES "kmer_index" + "uniprot_gene_index")
# --------------------------------------------------------------------


def _generate_kmers(seq: str, k: int) -> List[str]:
    """
    Generate all overlapping k-mers from the query sequence.
    """
    seq = seq.strip().upper()
    n = len(seq)
    if n < k:
        return []
    return [seq[i : i + k] for i in range(n - k + 1)]


def _kmer_approx_search_hits(
    query_seq: str,
    db_path: Path | str | None = None,
    k: int = KMER_K,
    min_similarity_percent: float = 0.0,
    limit: int = DEFAULT_SQL_LIMIT,
    debug: bool = False,
) -> List[Hit]:
    """
    Run approximate k-mer search against the local SQLite DB and
    return Hit instances sorted by similarity (descending).

    Similarity definition:
        similarity_raw = distinct_shared_kmers / distinct_query_kmers
        similarity_percent = 100 * similarity_raw, clamped to [0, 100].

    Hit.score stores similarity_percent (0..100).
    """
    query_seq = query_seq.strip().upper()
    if not query_seq:
        if debug:
            print("[UniProt seq_search] K-mer search: empty query; no hits.")
        return []

    # All k-mers, then dedupe for counting
    query_kmers_all = _generate_kmers(query_seq, k)
    if not query_kmers_all:
        if debug:
            print(
                f"[UniProt seq_search] K-mer search: query too short for k={k}; "
                "no hits."
            )
        return []

    query_kmers = sorted(set(query_kmers_all))
    total_query_kmers = len(query_kmers)
    if total_query_kmers == 0:
        if debug:
            print("[UniProt seq_search] K-mer search: no distinct k-mers; no hits.")
        return []

    if debug:
        print(
            f"[UniProt seq_search] K-mer search: "
            f"distinct_query_kmers={total_query_kmers}, k={k}"
        )

    path = _resolve_db_path(db_path)
    conn = sqlite3.connect(str(path))
    try:
        cur = conn.cursor()
        placeholders = ", ".join("?" for _ in query_kmers)

        sql = f"""
        SELECT
            u.acc,
            u.canonical_acc,
            u.sequence,
            u.gene_name,
            u.entrez_gene_id,
            COUNT(DISTINCT ki.kmer) AS shared_kmers
        FROM kmer_index AS ki
        JOIN uniprot_gene_index AS u
          ON u.acc = ki.acc
        WHERE ki.kmer IN ({placeholders})
        GROUP BY
            u.acc,
            u.canonical_acc,
            u.sequence,
            u.gene_name,
            u.entrez_gene_id
        HAVING shared_kmers > 0
        ORDER BY shared_kmers DESC
        LIMIT ?;
        """

        cur.execute(sql, (*query_kmers, limit))
        rows = cur.fetchall()
    finally:
        conn.close()

    if debug:
        print(
            "[UniProt seq_search] K-mer search: raw rows from SQL "
            f"(shared_kmers > 0): {len(rows)}"
        )

    hits: List[Hit] = []
    if not rows:
        if debug:
            print("[UniProt seq_search] K-mer search: no sequences share k-mers.")
        return hits

    query_len = len(query_seq)

    for acc, canonical_acc, sequence, gene_name, entrez_gene_id, shared_kmers in rows:
        shared = int(shared_kmers)

        similarity_raw = shared / float(total_query_kmers)
        # Guard against any unexpected >1 due to quirks
        if similarity_raw > 1.0:
            similarity_raw = 1.0
        if similarity_raw < 0.0:
            similarity_raw = 0.0

        similarity_percent = 100.0 * similarity_raw
        coverage = similarity_raw
        align_len = int(round(coverage * query_len))

        h = Hit(
            acc=str(acc),
            canonical_acc=str(canonical_acc),
            gene_name=gene_name,
            entrez_gene_id=entrez_gene_id,
            score=similarity_percent,
            coverage=coverage,
            alignment_length=align_len,
            uniprot_sequence=sequence,
        )
        hits.append(h)

    if not hits:
        if debug:
            print("[UniProt seq_search] K-mer search: no raw hits after rows.")
        return hits

    hits.sort(key=lambda h: h.score, reverse=True)

    if debug:
        sims = [h.score for h in hits]
        print(
            "[UniProt seq_search] K-mer search: similarity_percent range "
            f"(raw hits): min={min(sims):.2f}, max={max(sims):.2f}"
        )
        print(
            "[UniProt seq_search] K-mer search: applying min_similarity_percent "
            f"threshold = {min_similarity_percent:.2f}"
        )

    if min_similarity_percent <= 0.0:
        return hits

    filtered = [h for h in hits if h.score >= min_similarity_percent]
    if not filtered and debug:
        print(
            "[UniProt seq_search] K-mer search: no hits above threshold; "
            "returning all raw hits."
        )
        return hits

    return filtered


# --------------------------------------------------------------------
# PUBLIC API
# --------------------------------------------------------------------


def resolve_sequences_to_gene_records(
    sequences: List[str],
    db_path: Path | str | None = None,
    top_n: int = 5,
    min_score: float = 80.0,
    debug: bool = False,
) -> Dict[str, Any]:
    """
    Resolve a list of amino acid sequences to gene-associated records.

    For each input sequence:
      - Perform substring search in table "uniprot".
      - Perform k-mer approximate search in "kmer_index" + "uniprot_gene_index".
      - Combine all hits from both stages.
      - Merge by (gene_name, acc), keeping the highest scoring Hit for each.
      - Sort hits so that entries with non-null gene_name come first,
        then by score (descending), then by canonical_acc, then by acc.
      - Apply top_n after merging and sorting.

    Parameters
    ----------
    sequences : List[str]
        List of amino acid sequences (one-letter codes).
    db_path : Path | str | None
        Path to the UniProt SQLite database. If None, UNIPROT_DB_PATH is used.
    top_n : int
        Maximum number of hits to keep per input sequence.
    min_score : float
        Minimum similarity score (percent, 0..100) to keep a hit.
    debug : bool
        If True, print detailed debug information.

    Returns
    -------
    Dict[str, Any]
        Dictionary with:
          - "genes": list of unique gene symbols (strings), suitable for
                     result.get("genes", []) -> ["ABC", "DEF", ...]
          - "records": flat list of per-hit records with the originating
                      sequence preserved.

        {
          "genes": ["ABC", "DEF", ...],
          "records": [
            {
              "sequence": "<original_sequence>",
              "gene_name": ...,
              "entrez_gene_id": ...,
              "acc": ...,
              "canonical_acc": ...,
              "score": ...,
              "coverage": ...,
              "alignment_length": ...,
              "uniprot_sequence": "<full_uniprot_sequence_for_acc>",
            },
            ...
          ]
        }
    """
    records: List[Dict[str, Any]] = []

    for idx, seq in enumerate(sequences):
        raw_seq = seq
        seq = (seq or "").strip().upper()

        if not seq:
            if debug:
                print(
                    f"[UniProt seq_search] Query {idx}: empty sequence, "
                    "skipping."
                )
            continue

        if debug:
            print(
                f"\n[UniProt seq_search] Processing query {idx} "
                f"(len={len(seq)})"
            )

        # 1) Substring-based hits (table "uniprot")
        substring_hits = _substring_search_hits(
            query_seq=seq,
            db_path=db_path,
            min_similarity_percent=min_score,
            debug=debug,
        )

        # 2) K-mer approximate hits (tables "kmer_index" + "uniprot_gene_index")
        kmer_hits = _kmer_approx_search_hits(
            query_seq=seq,
            db_path=db_path,
            min_similarity_percent=min_score,
            limit=DEFAULT_SQL_LIMIT,
            debug=debug,
        )

        # 3) Merge all hits and keep best per (gene_name, acc)
        best_by_key: Dict[Tuple[str, str], Hit] = {}

        def _update_best(h: Hit) -> None:
            gene_key = h.gene_name or ""
            key = (gene_key, h.acc)
            if h.score < min_score:
                return
            existing = best_by_key.get(key)
            if existing is None or h.score > existing.score:
                best_by_key[key] = h

        for h in substring_hits:
            _update_best(h)

        for h in kmer_hits:
            _update_best(h)

        hits = list(best_by_key.values())

        # 4) Sort:
        #    - prefer non-null gene_name
        #    - then by score descending
        #    - then by canonical_acc
        #    - then by acc
        hits.sort(
            key=lambda h: (
                h.gene_name is None,
                -h.score,
                h.canonical_acc,
                h.acc,
            )
        )

        # 5) Apply top_n
        if top_n > 0:
            hits = hits[:top_n]

        # 6) Convert to plain dicts and append to global records list
        for h in hits:
            d = asdict(h)

            # Defensive clamping to keep invariants:
            #   score in [0, 100], coverage in [0, 1],
            #   alignment_length <= query length.
            score = float(d.get("score") or 0.0)
            if score < 0.0:
                score = 0.0
            if score > 100.0:
                score = 100.0

            coverage = float(d.get("coverage") or 0.0)
            if coverage < 0.0:
                coverage = 0.0
            if coverage > 1.0:
                coverage = 1.0

            align_len_raw = int(d.get("alignment_length") or 0)
            if align_len_raw < 0:
                align_len_raw = 0
            alignment_length = min(align_len_raw, len(seq))

            records.append(
                {
                    "sequence": raw_seq,  # keep original sequence as a field
                    "gene_name": d["gene_name"],
                    "entrez_gene_id": d["entrez_gene_id"],
                    "acc": d["acc"],
                    "canonical_acc": d["canonical_acc"],
                    "score": score,
                    "coverage": coverage,
                    "alignment_length": alignment_length,
                    "uniprot_sequence": d.get("uniprot_sequence") or "",
                }
            )

    # 7) Build unique list of gene names (strings) for convenience
    genes: List[str] = []
    seen: set[str] = set()
    for rec in records:
        gene_name = rec.get("gene_name")
        if not gene_name:
            continue
        if gene_name not in seen:
            seen.add(gene_name)
            genes.append(gene_name)

    return {
        "genes": genes,       # what result.get("genes", []) will return
        "records": records,   # full records with sequence and scores
    }


__all__ = [
    "Hit",
    "resolve_sequences_to_gene_records",
]


# --------------------------------------------------------------------
# SELF-TEST (simple demo)
# --------------------------------------------------------------------

if __name__ == "__main__":
    DEBUG = True

    aa_sequences: List[str] = [
        "DIGDSFGHPACPLVSRSRNSPVEVDDDEDDVVFTEIIQPPSASWPKIADQRNFIFASSKNEKHKGNYSVIPPSSRDLASQKGNISETIVIDDEEDIETNGGARKKSSCWIEWTLPGTKNK"
    ]

    result = resolve_sequences_to_gene_records(
        aa_sequences,
        top_n=100,
        min_score=50.0,
        debug=DEBUG,
    )

    print("\n[UniProt seq_search DEMO] result:")
    print(json.dumps(result, indent=2))

    print("\n[UniProt seq_search DEMO] genes only:")
    print(json.dumps(result.get("genes", []), indent=2))
