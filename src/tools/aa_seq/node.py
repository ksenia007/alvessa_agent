# src/tools/aa_seq/node.py
# ===========================================================
# AA Sequence-to-Gene Tool Node
# ===========================================================
#
# Description:
#
# This node integrates the AA sequence resolver with the Alvessa
# agent framework.
#
# INPUT:
#   - amino acid sequences (NOT genes).
#
# BEHAVIOR:
#   - Reads AA sequences from the State (aa_sequences / sequences / aaseq / aa_sequence).
#   - Calls resolve_sequences_to_gene_records(...) to get UniProt mapping.
#   - Populates:
#       * state["aa_seq_result"]   : structured resolver output
#       * state["aa_seq_summary"]  : global human-readable text summary
#       * state["aa_seq_html"]     : interactive HTML viewer for AA-seq mappings
#   - Per-gene text summaries:
#       * For each gene_name found in the result, attaches a per-gene summary
#         to Gene objects via Gene.update_text_summaries(...).
#       * If a gene is not present in state["gene_entities"], creates it and
#         appends its symbol to state["genes"].
#   - Records without gene_name (unmapped accessions):
#       * Included at the end of aa_seq_summary (before interpretation notes).
#       * If there are no genes in the state at all and only unmapped records
#         exist, a synthetic gene 'UNKNOWN_GENE' is created and those records
#         are attached to it as text summaries (last resort).
#

from __future__ import annotations

import sys
import json
from typing import Any, Dict, List

from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.config import DEBUG
from src.state import State
from src.tools.base import Node
from src.alvessa.domain.gene_class import Gene, canon_gene_key
from src.alvessa.domain.gene_components import GeneIdentifiers

from src.tools.aa_seq import OUTPUT_DIR, HTML_TEMPLATE, CSS_TEMPLATE, JS_TEMPLATE
from src.tools.aa_seq.seq_search import resolve_sequences_to_gene_records
from src.tools.aa_seq.utils import (
    log,
    normalize_sequences,
    summarize_aa_seq_result,
    build_per_gene_summaries,
    inject_frontend_assets_aa_seq,
)


# --------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------

# Default parameters for the underlying resolver
AA_SEQ_TOP_N: int = 5
AA_SEQ_MIN_SCORE: float = 60.0


# --------------------------------------------------------------------
# INTERNAL HELPERS
# --------------------------------------------------------------------


def _prepare_sequences_aa_seq(state: "State") -> List[str]:
    """
    Extract amino acid sequences from the State and normalize them.

    Input is expected to be one of:
      - state["aa_sequences"]
      - state["sequences"]
      - state["aaseq"]
      - state["aa_sequence"]

    The canonical normalized list is written back into:
      - state["aa_sequences"]
    """
    raw = (
        state.get("aa_sequences")
        or state.get("sequences")
        or state.get("aaseq")
        or state.get("aa_sequence")
    )

    seqs = normalize_sequences(raw, min_len=10)
    state["aa_sequences"] = seqs

    used = state.get("used_tools", []) or []
    if "aa_seq" not in used:
        used.append("aa_seq")
    state["used_tools"] = used

    return seqs


def _collect_unmapped_records(records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Return records where gene_name is missing or empty."""
    unmapped: List[Dict[str, Any]] = []
    for rec in records or []:
        gene = (rec.get("gene_name") or "").strip()
        if not gene:
            unmapped.append(rec)
    return unmapped


def _build_unmapped_block_for_unknown_gene(unmapped: List[Dict[str, Any]]) -> str:
    """
    Build a text block suitable for attaching to the UNKNOWN_GENE object,
    based on unmapped records.
    """
    if not unmapped:
        return ""

    lines: List[str] = []
    lines.append("AA sequence to UniProt mapping (no gene symbol available).")
    lines.append("These UniProt entries were matched by sequence but have no gene_name.")
    lines.append("")

    # Reuse internal formatting helpers from utils
    from src.tools.aa_seq.utils import _fmt_score, _fmt_float_sig

    for i, rec in enumerate(unmapped, start=1):
        acc = rec.get("acc", "") or ""
        canonical_acc = rec.get("canonical_acc", "") or ""
        score = _fmt_score(rec.get("score"))
        coverage = _fmt_float_sig(rec.get("coverage"))
        align_len = rec.get("alignment_length")
        query_seq = (rec.get("sequence") or "").strip()
        uniprot_seq = (rec.get("uniprot_sequence") or "").strip()

        lines.append(f"  Unmapped Hit {i}:")
        lines.append(f"    Accession: {acc}")
        lines.append(f"    Canonical Accession: {canonical_acc}")
        lines.append(f"    Score Percent: {score}")
        lines.append(f"    Coverage: {coverage}")
        lines.append(f"    Alignment Length: {align_len if align_len is not None else 'N/A'}")

        if query_seq:
            lines.append("    Query Sequence:")
            lines.append(f"      {query_seq}")
        if uniprot_seq:
            lines.append("    Uniprot Reference Sequence:")
            lines.append(f"      {uniprot_seq}")

    return "\n".join(lines)


# --------------------------------------------------------------------
# MAIN AGENT
# --------------------------------------------------------------------


def aa_seq_agent(state: "State") -> "State":
    """
    Main AA sequence resolver agent.

    INPUT:
      - Amino acid sequences in the State (not genes), using keys:
          * "aa_sequences" (preferred)
          * "sequences"
          * "aaseq"
          * "aa_sequence"

    STEPS:
      1) Collect and normalize AA sequences from the State.
      2) Call resolve_sequences_to_gene_records(...) from seq_search.py.
      3) Populate:
           - state["aa_seq_result"]   : raw resolver output
           - state["aa_seq_summary"]  : human-readable text summary
           - state["aa_seq_html"]     : interactive HTML viewer
      4) For each gene_name in the result:
           - Ensure a Gene object exists in state["gene_entities"].
           - Attach a per-gene text summary via update_text_summaries().
      5) For records without gene_name:
           - They appear in the global aa_seq_summary.
           - If there were no genes in the state at all and only unmapped
             records exist, create an UNKNOWN_GENE and attach them there.

    OUTPUT:
      - Updated State with:
          * per-gene text summaries for sequence-derived evidence
          * aa_seq_result / aa_seq_summary / aa_seq_html for UI and downstream tools.
    """
    had_genes_before = bool((state or {}).get("gene_entities"))

    sequences = _prepare_sequences_aa_seq(state)
    if not sequences:
        log("No valid amino acid sequences found in State.")
        state["aa_seq_result"] = {"genes": [], "records": []}
        state["aa_seq_summary"] = (
            "AA sequence UniProt mapping: no valid input sequences were provided."
        )
        # In this case we do not generate HTML
        return state

    log(
        f"Resolving {len(sequences)} AA sequences "
        f"(top_n={AA_SEQ_TOP_N}, min_score={AA_SEQ_MIN_SCORE:.1f})"
    )

    result = resolve_sequences_to_gene_records(
        sequences=sequences,
        db_path=None,
        top_n=AA_SEQ_TOP_N,
        min_score=AA_SEQ_MIN_SCORE,
        debug=DEBUG,
    )

    # Ensure both "records" and "gene_records" keys are present for downstream users
    if "records" not in result and "gene_records" in result:
        result["records"] = result["gene_records"]
    if "gene_records" not in result:
        result["gene_records"] = result.get("records", []) or []

    # Structured data for future UI tables
    state["aa_seq_result"] = result

    # Human-readable text summary (includes query and UniProt sequences
    # and the interpretation note at the end).
    state["aa_seq_summary"] = summarize_aa_seq_result(result)

    # --------------------------------------------------------------
    # Per-gene text summaries
    # --------------------------------------------------------------
    per_gene_summaries = build_per_gene_summaries(result)
    records: List[Dict[str, Any]] = result.get("records", []) or []
    unmapped_records = _collect_unmapped_records(records)

    gene_entities: Dict[str, Gene] = state.get("gene_entities", {}) or {}
    genes_list: List[str] = state.get("genes", []) or []

    # Attach mapping evidence to each gene that has hits
    for gene_name, block in per_gene_summaries.items():
        key = canon_gene_key(gene_name)
        gene_obj = gene_entities.get(key)

        # representative record for IDs
        rep_rec = None
        for r in records:
            g = (r.get("gene_name") or "").strip()
            if g == gene_name:
                rep_rec = r
                break

        if gene_obj is None:
            ids = GeneIdentifiers(symbol=gene_name)
            if rep_rec is not None:
                entrez = rep_rec.get("entrez_gene_id")
                if entrez:
                    ids.entrez_id = str(entrez)
                uniprot_id = rep_rec.get("canonical_acc") or rep_rec.get("acc")
                if uniprot_id:
                    ids.uniprot_id = uniprot_id

            gene_obj = Gene(identifiers=ids)
            gene_obj.add_tool("AASeq")
            gene_entities[key] = gene_obj
            if gene_name not in genes_list:
                genes_list.append(gene_name)
        else:
            gene_obj.add_tool("AASeq")

        summary_line_parts = []
        if rep_rec:
            acc = (rep_rec.get("canonical_acc") or rep_rec.get("acc") or "").strip()
            coverage = rep_rec.get("coverage")
            align_len = rep_rec.get("alignment_length")
            query_seq = (rep_rec.get("sequence") or "").strip()
            uniprot_seq = (rep_rec.get("uniprot_sequence") or "").strip()
            q_len = len(query_seq) if query_seq else None
            u_len = len(uniprot_seq) if uniprot_seq else None
            diff_len = (q_len - u_len) if (q_len is not None and u_len is not None) else None
            coverage_str = f"{coverage:.2f}" if isinstance(coverage, (int, float)) else str(coverage) if coverage is not None else "N/A"
            summary_line = f"*AASeq: Matched UniProt {acc or 'unknown'}"
            summary_bits = []
            if coverage is not None:
                summary_bits.append(f"coverage {coverage_str}")
            if align_len is not None:
                summary_bits.append(f"alignment_len {align_len}")
            if q_len is not None:
                summary_bits.append(f"query_len {q_len}")
            if u_len is not None:
                summary_bits.append(f"uniprot_len {u_len}")
            if diff_len is not None:
                summary_bits.append(f"length_diff {diff_len}")
            if summary_bits:
                summary_line += " (" + "; ".join(summary_bits) + ")"
            summary_line += "."
            summary_line_parts.append(summary_line)

        # Attach per-gene block with evidence
        summary_line_parts.append("*AASeq: Evidence: " + block.replace("\n", " ").strip())
        gene_obj.update_text_summaries(" ".join(summary_line_parts))

    # --------------------------------------------------------------
    # Handle unmapped records with UNKNOWN_GENE (last resort)
    # --------------------------------------------------------------
    if unmapped_records and not had_genes_before and not per_gene_summaries:
        # No genes in the state and only unmapped records: create UNKNOWN_GENE
        key = "UNKNOWN_GENE"
        if key in gene_entities and isinstance(gene_entities[key], Gene):
            unknown_gene = gene_entities[key]
        else:
            ids = GeneIdentifiers(symbol=key)
            unknown_gene = Gene(identifiers=ids)
            unknown_gene.add_tool("AASeq")
            gene_entities[key] = unknown_gene
            if key not in genes_list:
                genes_list.append(key)

        unknown_block = _build_unmapped_block_for_unknown_gene(unmapped_records)
        if unknown_block:
            unknown_gene.update_text_summaries(
                "*AASeq: AA sequence to UniProt mapping evidence (no gene symbol): "
                + unknown_block.replace("\n", " ").strip()
            )

    # Write back gene_entities and genes list
    state["gene_entities"] = gene_entities
    state["genes"] = genes_list

    # --------------------------------------------------------------
    # Build HTML frontend (interactive AA-seq viewer)
    # --------------------------------------------------------------
    try:
        with open(HTML_TEMPLATE, "r", encoding="utf-8") as f:
            html_template = f.read()
        with open(CSS_TEMPLATE, "r", encoding="utf-8") as f:
            css_template = f.read()
        with open(JS_TEMPLATE, "r", encoding="utf-8") as f:
            js_template = f.read()
    except FileNotFoundError as e:
        # Mirror chembl behavior: escalate missing template as a runtime error
        raise RuntimeError(f"Missing AA-seq frontend template: {e.filename}")

    # inject_frontend_assets_aa_seq takes the raw result + sequences and
    # builds the aaSeqData object inside the HTML/JS template.
    aa_sequences = state.get("aa_sequences", []) or []
    aa_seq_html = inject_frontend_assets_aa_seq(
        html_template=html_template,
        css_template=css_template,
        js_template=js_template,
        aa_seq_result=result,
        sequences=aa_sequences,
    )

    if aa_seq_html:
        state["aa_seq_html"] = aa_seq_html

    return state


# --------------------------------------------------------------------
# CLI (testing mode)
# --------------------------------------------------------------------
if __name__ == "__main__":
    """
    CLI testing mode for the AA sequence node.

    USAGE EXAMPLES (INPUT IS AA SEQUENCES, NOT GENES):

        python node.py "MTEYKLVVVGAGGVGKSALTIQLIQNHFVD"

        python node.py "
            MTEYKLVVVGAGGVGKSALTIQLIQNHFVD
            MAAAGITSLILV...
        "

    If no sequences are provided, a short HRAS fragment is used as a demo.
    """
    raw_args = sys.argv[1:]
    if not raw_args:
        # fallback demo sequence (HRAS fragment, or any other test sequence)
        # raw_args = ["MTEYKLVVVGAGGVGKSALTIQLIQNHFVD"]
        #raw_args = ["CKCWHGQLRCFPQAFLPGCDGLVMDEHLVA"]
        #raw_args = ["LENLQIIRGN", "LENLQIIRGG", "LENLQIIGGG", "VDRVDYDRQSGSAVITFVEI","DDRVDYDRQSGSAVITFVEI"]
        raw_args = ["SWPSRSLPAQEAVVKLRVEGMTCQSCVSSIEGKVRKLQGVVRVKVSLSNQEAVITYQPYLIQPEDLRDHVNDMGFEAAIKSKVAPLSLGPIDIERLQSTNPKRPLSSANQNFNNSETLGH"]
        

    sequences = raw_args

    # Base name for output files: derived from first AA sequence
    first_seq = "".join(sequences[0].split()).upper()
    first_id = first_seq[:10] if len(first_seq) >= 10 else (first_seq or "seq")

    if len(sequences) == 1:
        base_name = f"{first_id}_aa"
    else:
        base_name = f"{first_id}_plus{len(sequences)-1}_aa"

    state: State = {
        "aa_sequences": sequences,
        "gene_entities": {},
        "genes": [],
    }

    result_state = aa_seq_agent(state)

    aa_seq_result = result_state.get("aa_seq_result", {}) or {}
    aa_seq_summary = result_state.get("aa_seq_summary", "") or ""
    aa_seq_html = result_state.get("aa_seq_html", "")

    # TXT summary (single block, ready to paste into logs or UI)
    txt_out = OUTPUT_DIR / f"{base_name}_aa_seq.txt"
    txt_out.write_text(
        aa_seq_summary if aa_seq_summary else "No AA-seq summary produced.",
        encoding="utf-8",
    )

    # JSON with raw resolver output (genes + records) for future UI tables
    json_out = OUTPUT_DIR / f"{base_name}_aa_seq_result.json"
    json_out.write_text(
        json.dumps(aa_seq_result, indent=2),
        encoding="utf-8",
    )

    # HTML viewer (if produced)
    html_out = OUTPUT_DIR / f"{base_name}_aa_seq.html"
    html_out.write_text(
        aa_seq_html if aa_seq_html else "<p>No AA-seq HTML produced.</p>",
        encoding="utf-8",
    )

    print("[OK] Generated outputs:")
    print(f"  * TXT file:   {txt_out.resolve()}")
    print(f"  * JSON file:  {json_out.resolve()}")
    print(f"  * HTML file:  {html_out.resolve()}")
    print(f"[INFO] Sequences processed: {len(sequences)}")


# --------------------------------------------------------------------
# Node registration
# --------------------------------------------------------------------

NODES: tuple[Node, ...] = (
    Node(
        name="aa_seq",
        entry_point=aa_seq_agent,
        description=(
            "Resolve amino acid sequences (not gene symbols) to UniProt and gene-level "
            "information using a local copy of the UniProt database. The tool updates "
            "per-gene text summaries and returns both a structured result, a global "
            "text summary, and an interactive HTML viewer for AA-sequence mappings."
        ),
    ),
)
