# src/tools/aa_seq/utils.py
# ===========================================================
# AA Sequence Tool Utilities
# ===========================================================
#
# Shared helpers for the aa_seq tool:
#   - logging
#   - sequence normalization
#   - text summarization of UniProt mapping results
#   - building frontend data and injecting HTML/CSS/JS for UI
#
# These utilities are intentionally UI-agnostic at the data level.
# They produce:
#   - plain text summaries (for Gene objects and logs)
#   - a structured JSON payload (aaSeqData) for the frontend
#

from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, Iterable, List

import json

from src.config import DEBUG


# --------------------------------------------------------------------
# Logging
# --------------------------------------------------------------------


def log(msg: str) -> None:
    """Simple timestamped logger for aa_seq tool."""
    if DEBUG:
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[aa_seq] {msg} @ {now}")


# --------------------------------------------------------------------
# Sequence normalization
# --------------------------------------------------------------------


def _as_iterable(value: Any) -> Iterable[str]:
    """Turn scalars or lists/tuples into a flat iterable of strings."""
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, (list, tuple, set)):
        return [str(v) for v in value if v is not None]
    return [str(value)]


def normalize_sequences(raw: Any, min_len: int = 10) -> List[str]:
    """
    Normalize arbitrary user input into a list of amino acid sequences.

    Rules:
      - Accepts str, list[str], tuple[str], or other iterables.
      - Strips whitespace.
      - Keeps only alphabetical characters and uppercases them.
      - Drops sequences shorter than min_len.
    """
    if raw is None:
        return []

    sequences: List[str] = []
    seen = set()

    for item in _as_iterable(raw):
        if not item:
            continue
        letters = "".join(ch for ch in item if ch.isalpha()).upper()
        if len(letters) < min_len:
            continue
        if letters in seen:
            continue
        seen.add(letters)
        sequences.append(letters)

    return sequences


# --------------------------------------------------------------------
# Text summarization helpers
# --------------------------------------------------------------------


def _fmt_score(value: Any) -> str:
    """Format similarity percent with two decimal places."""
    try:
        return f"{float(value):.2f}"
    except Exception:
        return "N/A"


def _fmt_float_sig(value: Any, digits: int = 3) -> str:
    """Format with a given number of significant digits (general format)."""
    try:
        return f"{float(value):.{digits}g}"
    except Exception:
        return "N/A"


def _group_records_by_gene(records: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    """
    Group flat hit records by gene_name.

    Returns:
        dict[gene_name, list[record]]
        Records with falsy or missing gene_name are not included here.
    """
    by_gene: Dict[str, List[Dict[str, Any]]] = {}
    for rec in records or []:
        gene = (rec.get("gene_name") or "").strip()
        if not gene:
            continue
        by_gene.setdefault(gene, []).append(rec)
    return by_gene


def _collect_unmapped_records(records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Return records that have no gene_name (unmapped UniProt accessions).
    """
    unmapped: List[Dict[str, Any]] = []
    for rec in records or []:
        gene = (rec.get("gene_name") or "").strip()
        if not gene:
            unmapped.append(rec)
    return unmapped


def _build_gene_block(gene: str, recs: List[Dict[str, Any]]) -> str:
    """
    Build a per-gene text block used both for:
      - per-gene summaries (attached to Gene objects)
      - global summary sections.
    """
    lines: List[str] = []

    # Pick a representative Entrez Gene ID if present
    entrez = None
    for r in recs:
        eid = r.get("entrez_gene_id")
        if eid:
            entrez = str(eid)
            break

    if entrez:
        lines.append(f"Gene: {gene}, Entrez Gene ID: {entrez}")
    else:
        lines.append(f"Gene: {gene}")

    for i, rec in enumerate(recs, start=1):
        acc = rec.get("acc", "") or ""
        canonical_acc = rec.get("canonical_acc", "") or ""
        score = _fmt_score(rec.get("score"))
        coverage = _fmt_float_sig(rec.get("coverage"))
        align_len = rec.get("alignment_length")

        query_seq = (rec.get("sequence") or "").strip()
        uniprot_seq = (rec.get("uniprot_sequence") or "").strip()

        lines.append(f"  Hit {i}:")
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


def build_per_gene_summaries(result: Dict[str, Any]) -> Dict[str, str]:
    """
    Build per-gene text summaries from the aa_seq resolver result.

    Returns:
        dict[gene_name, summary_text]
    """
    records: List[Dict[str, Any]] = result.get("records", []) or []
    by_gene = _group_records_by_gene(records)

    summaries: Dict[str, str] = {}
    for gene in sorted(by_gene.keys()):
        summaries[gene] = _build_gene_block(gene, by_gene[gene])

    return summaries


def _build_unmapped_block(unmapped: List[Dict[str, Any]]) -> str:
    """
    Build a text block for records that have no associated gene_name.
    These are UniProt accessions without an obvious gene symbol.
    """
    if not unmapped:
        return ""

    lines: List[str] = []
    lines.append("Records without associated gene_name (UniProt accessions only):")
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


def interpretation_notes() -> str:
    """
    Interpretation notes appended at the end of the global summary to help
    both human users and agent-level reasoning.
    """
    lines: List[str] = []
    lines.append("Interpretation Notes")
    lines.append("")
    lines.append("- This mapping is based on amino acid sequence similarity against a local UniProtKB snapshot.")
    lines.append("- Score Percent is a local, gapless similarity between the query sequence and the best window")
    lines.append("  in the UniProt reference sequence (0 to 100; higher is more similar).")
    lines.append("- Coverage reflects the fraction of the query sequence that is effectively aligned or supported")
    lines.append("  by shared k-mers. It is a heuristic, not a full gapped alignment.")
    lines.append("- Canonical Accession corresponds to the root UniProt ID; isoform accessions may have the form")
    lines.append("  ACC-1, ACC-2, etc.")
    lines.append("- Entrez Gene IDs are taken from the UniProt record when available and may not be present for")
    lines.append("  all accessions.")
    lines.append("- Records listed under 'Records without associated gene_name' are UniProt entries that lack a")
    lines.append("  clear gene symbol mapping in the local database.")
    lines.append("- This tool does not replace a full multiple-alignment or domain-level homology analysis;")
    lines.append("  it is intended as a fast, approximate sequence-to-gene resolver.")
    return "\n".join(lines)


def summarize_aa_seq_result(result: Dict[str, Any]) -> str:
    """
    Build a global, human-readable summary for the aa_seq result.

    Structure:
      - Header with total hit statistics.
      - Per-gene sections.
      - Unmapped UniProt entries (no gene_name).
      - Interpretation notes.
    """
    records: List[Dict[str, Any]] = result.get("records", []) or []
    by_gene = _group_records_by_gene(records)
    unmapped = _collect_unmapped_records(records)

    lines: List[str] = []
    lines.append("AA sequence UniProt mapping summary")

    lines.append(f"  Total hits: {len(records)}")
    lines.append(f"  Distinct genes with hits: {len(by_gene)}")
    lines.append("")

    # Per-gene blocks
    for gene in sorted(by_gene.keys()):
        lines.append(_build_gene_block(gene, by_gene[gene]))
        lines.append("")

    # Unmapped entries, if any
    if unmapped:
        lines.append(_build_unmapped_block(unmapped))
        lines.append("")

    # Interpretation notes at the end
    lines.append(interpretation_notes())

    return "\n".join(lines)


# --------------------------------------------------------------------
# Frontend data builder and HTML injectors
# --------------------------------------------------------------------


def build_frontend_data(result: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert the raw resolver result into a frontend-friendly aaSeqData payload.

    Input (result) is expected to have:
      - "genes": list of gene symbols (optional; as in your sample JSON)
      - "records" or "gene_records": flat list of hit dictionaries

    Output (aaSeqData) is shaped as:

      {
        "genes": {
          "CRIPTO": {
            "entrez_gene_id": "6997",
            "hits": [ {..hit1..}, {..hit2..}, ... ]
          },
          "CRIPTO3": { ... }
        },
        "unmapped": [ {..hit with no gene_name..}, ... ],
        "meta": {
          "total_hits": int,
          "n_genes": int
        }
      }
    """
    # Prefer "records", fall back to "gene_records" for compatibility
    records: List[Dict[str, Any]] = (
        result.get("records")
        or result.get("gene_records")
        or []
    )

    # Group by gene
    by_gene = _group_records_by_gene(records)
    unmapped = _collect_unmapped_records(records)

    genes_payload: Dict[str, Any] = {}
    for gene, recs in sorted(by_gene.items()):
        # Representative Entrez ID
        entrez = None
        for r in recs:
            eid = r.get("entrez_gene_id")
            if eid:
                entrez = str(eid)
                break

        # Normalize hits: keep only relevant fields for the UI
        hits: List[Dict[str, Any]] = []
        for r in recs:
            hits.append(
                {
                    "sequence": (r.get("sequence") or "").strip(),
                    "gene_name": (r.get("gene_name") or "").strip(),
                    "entrez_gene_id": r.get("entrez_gene_id") or "",
                    "acc": r.get("acc") or "",
                    "canonical_acc": r.get("canonical_acc") or "",
                    "score": r.get("score"),
                    "coverage": r.get("coverage"),
                    "alignment_length": r.get("alignment_length"),
                    "uniprot_sequence": (r.get("uniprot_sequence") or "").strip(),
                }
            )

        genes_payload[gene] = {
            "entrez_gene_id": entrez or "",
            "n_hits": len(hits),
            "hits": hits,
        }

    aa_seq_data: Dict[str, Any] = {
        "genes": genes_payload,
        "unmapped": unmapped,
        "meta": {
            "total_hits": len(records),
            "n_genes": len(genes_payload),
        },
    }

    return aa_seq_data


def inject_frontend_assets(
    *,
    html_template: str,
    css_template: str,
    js_template: str,
    aa_seq_data: Dict[str, Any],
    title: str,
) -> str:
    """
    Inject CSS, JS, and aaSeqData into the HTML template.

    Template placeholders (mirroring ChEMBL and prot tools):

      {{CSS_INLINE}}  - replaced with <style>...</style>
      {{JS_INLINE}}   - replaced with <script>const aaSeqData = ...; ...</script>
      {{GENE_SYMBOL}} - replaced with a string title used in the <title> tag
                        and possibly in the header.
    """
    css_block = "<style>\n" + css_template + "\n</style>"
    js_block = (
        "<script>\n"
        "const aaSeqData = "
        + json.dumps(aa_seq_data, ensure_ascii=False)
        + ";\n"
        + js_template
        + "\n</script>"
    )

    html = html_template.replace("{{CSS_INLINE}}", css_block)
    html = html.replace("{{JS_INLINE}}", js_block)
    html = html.replace("{{GENE_SYMBOL}}", title)

    return html


def inject_frontend_assets_aa_seq(
    *,
    html_template: str,
    css_template: str,
    js_template: str,
    aa_seq_result: Dict[str, Any],
    sequences: List[str] | None = None,
) -> str:
    """
    Convenience wrapper used by aa_seq.node:

      - Converts the raw resolver result into a frontend payload
        (build_frontend_data).
      - Adds the original AA sequences into the meta block.
      - Derives a title from the gene list, falling back to
        "AA sequence viewer".
      - Delegates to inject_frontend_assets(...) to actually inject
        CSS/JS + data into the HTML template.
    """
    # Build core payload from resolver result
    aa_seq_data = build_frontend_data(aa_seq_result)

    # Attach sequences into meta (for UI display)
    meta = aa_seq_data.setdefault("meta", {})
    meta["sequences"] = list(sequences or [])

    # Title: join distinct gene symbols if available
    gene_symbols = sorted(aa_seq_data.get("genes", {}).keys())
    if gene_symbols:
        title = ", ".join(gene_symbols)
    else:
        title = "AA sequence viewer"

    return inject_frontend_assets(
        html_template=html_template,
        css_template=css_template,
        js_template=js_template,
        aa_seq_data=aa_seq_data,
        title=title,
    )
