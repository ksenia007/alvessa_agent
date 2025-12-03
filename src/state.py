"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors:
Created: 2024-06-25
Updated: 2025-06-26

Description:

TypedDict schema describing the mutable LangGraph state + functions to extract info from it.
"""

from __future__ import annotations

import operator
import csv
import json
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Annotated, Any, Dict, List, Optional, Tuple

from typing_extensions import TypedDict

from src.alvessa.domain.gene_class import Gene
from src.alvessa.domain.variant_class import Variant
from src.alvessa.domain.drug_class import Drug

def merge_html(old: Optional[str], new: Optional[str]) -> Optional[str]:
    """
    Reducer for HTML blobs (prot_html, chembl_html, aa_seq_html, drug_central_html).
    Latest non-empty value wins instead of concatenating strings.
    """
    if new is None or new == "":
        return old
    return new

class State(TypedDict, total=False):
    # Conversation
    messages: Annotated[List[Dict[str, str]], operator.add]
    genes: Annotated[List[str], operator.add]
    traits: Annotated[List[str], operator.add]
    proteins: Annotated[List[str], operator.add]
    drugs: Annotated[List[str], operator.add]
    transcripts: Annotated[List[str], operator.add]
    variants: Annotated[Dict[str, Dict[str, Dict[str, Any]]], operator.or_]
    chr_pos_variants: Annotated[Dict[str, Dict[str, Dict[str, Any]]], operator.or_]
    gene_level_gencode: Annotated[Dict[str, Dict[str, Any]], operator.or_]
    prompt: Annotated[str, operator.add]
    mc_setup: Annotated[bool, operator.or_]
    aa_sequences: Annotated[List[str], operator.add]

    # AA Sequence-to-Gene data
    aa_seq_result: Annotated[Dict[str, Any], operator.or_]
    aa_seq_summary: Annotated[str, operator.add]

    # Object-level entities
    gene_entities: Annotated[Dict[str, "Gene"], operator.or_]
    variant_entities: Annotated[Dict[str, "Variant"], operator.or_]
    drug_entities: Annotated[Dict[str, "Drug"], operator.or_]

    # LLM bookkeeping
    context_block: Annotated[str, operator.add]
    llm_json: Annotated[Dict[str, Any], operator.or_]
    verification: str
    verify_attempts: int
    tool_updates: int
    used_tools: Annotated[List[str], operator.add]
    use_tools: Annotated[List[str], operator.add]  # tools to use in the current run

    # Interactive view
    ui: Annotated[Dict[str, Any], operator.or_]  # e.g., {"panels": [...]}

    # Protein structure and druggablility visualization
    #prot_html: Annotated[str, operator.add]
    prot_html: Annotated[str, merge_html]

    # For ChemBL drug-target interactive viewer
    #chembl_html: Annotated[str, operator.add]
    chembl_html: Annotated[str, merge_html]

    # For DrugCentral drug-target interactive viewer
    #drug_central_html: Annotated[str, operator.add]
    drug_central_html: Annotated[str, merge_html]

    # For AA Sequence-to-Gene interactive viewer
    #aa_seq_html: Annotated[str, operator.add]
    aa_seq_html: Annotated[str, merge_html]

    # General free-text annotations (not tied to a specific gene)
    text_notes: Annotated[List[str], operator.add]


# =========================
# SAVE FILES
# =========================
def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _safe_name(s: Optional[str], fallback: str = "UNKNOWN") -> str:
    if not s:
        return fallback
    s = s.strip()
    if not s:
        return fallback
    return "".join(c for c in s if c.isalnum() or c in ("_", "-", ".", ":"))


def _write_csv(path: Path, rows: List[Dict[str, Any]], field_order: List[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=field_order)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in field_order})


def _write_tsv(path: Path, rows: List[Dict[str, Any]], field_order: Optional[List[str]] = None) -> None:
    if not rows and not field_order:
        path.touch()
        return
    keys = field_order or list({k for r in rows for k in r.keys()})
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in keys})


# =========================
# Gene extractors
# =========================
def _gene_index_row(g: Gene) -> Dict[str, Any]:
    loc = g.get_location() or (None, None, None, None)
    chrom, start, end, strand = loc
    return {
        "symbol": g.symbol or "",
        "gene_type": g.gene_type or "",
        "entrez_id": g.entrez_id or "",
        "uniprot_id": g.uniprot_id or "",
        "ensembl_id": g.ensembl_id or "",
        "chrom": chrom or "",
        "start": start or "",
        "end": end or "",
        "strand": strand or "",
        "n_transcripts": getattr(g.transcriptome, "transcript_count", None) or g.get_n_transcripts(),
        "has_gwas": int(g.has_gwas_associations()),
        "gwas_total": g.get_all_gwas_associations() or 0,
        "n_variants": len(g.variants or {}),
        "n_interaction_partners": len(g.all_interaction_partners()) if g.has_interactions_collected() else 0,
        "n_binding_peaks": g.get_n_unique_binding_peaks() or 0,
        "n_mirna_targets": sum(len(v or []) for v in (g.mirna_targets or {}).values()),
        "n_traits": len(g.get_all_traits() or []),
        "tools_run": ";".join(g.tools_run or []),
    }


def _gene_transcripts_rows(g: Gene) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for tx_id, tx_info in (g.transcriptome.transcripts or {}).items():
        rows.append({"transcript_id": tx_id, "n_exons": tx_info.get("n_exons", ""), "source": "transcripts"})
    for iso_id, rec in (g.transcriptome.isoforms or {}).items():
        rows.append(
            {
                "transcript_id": iso_id,
                "name": rec.get("name", ""),
                "status": rec.get("status", ""),
                "aliases": ",".join(rec.get("aliases", []) or []),
                "locations": ",".join(rec.get("locations", []) or []),
                "notes": " | ".join(rec.get("notes", []) or []),
                "source": "isoforms",
            }
        )
    return rows


def _gene_interactions_rows(g: Gene, *, nonhuman: bool = False) -> List[Dict[str, Any]]:
    table = g.interactions.nonhuman_interactions if nonhuman else g.interactions.human_interactions
    rows: List[Dict[str, Any]] = []
    for exp_type, partners in (table or {}).items():
        for partner in partners or []:
            rows.append({"experiment_type": exp_type, "partner_symbol": partner})
    return rows


# =========================
# Drug extractors
# =========================
def _drug_index_row(d: Drug) -> Dict[str, Any]:
    """
    Compact index row for drugs.

    Focuses strictly on identifiers + basic counts, consistent with the
    lean Drug/DrugIdentifiers design (no derived properties).
    """
    ids = d.identifiers
    return {
        "name": ids.name or "",
        "chembl_id": ids.chembl_id or "",
        "drugcentral_id": ids.drugcentral_id or "",
        "cas_number": ids.cas_number or "",
        "catalog_number": ids.catalog_number or "",
        "n_synonyms": len(ids.synonyms or []),
        "n_mentions": len(d.mentions or []),
        "tools_run": ";".join(d.tools_run or []),
    }


# =========================
# Variant extractors (CSV-only, no per-variant dirs)
# =========================
def _best_variant_id(v: Variant) -> str:
    if v.rsID:
        return _safe_name(v.rsID)
    if v.loc_by_build:
        first_build = next(iter(v.loc_by_build))
        chrom = str(v.loc_by_build[first_build].get("chrom", "")).replace("chr", "").upper()
        pos = str(v.loc_by_build[first_build].get("pos", ""))
        if chrom and pos:
            return _safe_name(f"{chrom}:{pos}")
    return "VARIANT"


def _variant_index_row(v: Variant) -> Dict[str, Any]:
    build, chrom, pos = "", "", ""
    if v.loc_by_build:
        first_build = next(iter(v.loc_by_build))
        build = first_build
        chrom = str(v.loc_by_build[first_build].get("chrom", ""))
        pos = str(v.loc_by_build[first_build].get("pos", ""))
    return {
        "variant_id": v.rsID or _best_variant_id(v),
        "organism": v.organism or "",
        "primary_build": build,
        "primary_chrom": chrom,
        "primary_pos": pos,
        "n_related_genes": len(v.genes_related_to or []),
        "n_traits": len(v.traits or []),
        "n_af_freqs": len(v.af_freqs or []),
        "n_func_pred_genes": len(v.functional_predictions or {}),
        "tools_run": ";".join(v.tools_run or []),
        "has_text_summary": int(bool(v.text_summary or v.variant_summaries)),
    }


def _variant_locations_rows(v_id: str, v: Variant) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for build, rec in (v.loc_by_build or {}).items():
        refs = rec.get("ref", [])
        alts = rec.get("alt", [])
        ref = ",".join(refs) if isinstance(refs, list) else (refs or "")
        alt = ",".join(alts) if isinstance(alts, list) else (alts or "")
        rows.append(
            {
                "variant_id": v_id,
                "build": build,
                "chrom": rec.get("chrom", ""),
                "pos": rec.get("pos", ""),
                "ref": ref,
                "alt": alt,
            }
        )
    return rows


def _variant_per_gene_traits_rows(v_id: str, v: Variant) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for gene, traits in (v.per_gene_traits or {}).items():
        for trait, metrics in (traits or {}).items():
            row = {"variant_id": v_id, "gene": gene, "trait": trait}
            if isinstance(metrics, dict):
                # project simple metrics; keep everything as raw_json for fidelity
                for k, val in metrics.items():
                    if isinstance(val, (int, float, str, bool)) or val is None:
                        row[k] = val
                row["raw_json"] = json.dumps(metrics, ensure_ascii=False)
            rows.append(row)
    return rows


def _variant_per_gene_context_rows(v_id: str, v: Variant) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for gene, ctx in (v.per_gene_context or {}).items():
        rows.append(
            {
                "variant_id": v_id,
                "gene": gene,
                "context": ctx.get("context", ""),
                "variant_category": ctx.get("variant_category", ""),
            }
        )
    return rows


def _variant_func_pred_rows(v_id: str, v: Variant) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for gene, tools in (v.functional_predictions or {}).items():
        for tool, scores in (tools or {}).items():
            sc_list = scores if isinstance(scores, list) else [scores]
            for sc in sc_list:
                rows.append({"variant_id": v_id, "gene": gene, "tool": tool, "score": sc})
    return rows


def _variant_af_rows(v_id: str, v: Variant) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for entry in (v.af_freqs or []):
        row = {
            "variant_id": v_id,
            "source": entry.get("source", ""),
            "population": entry.get("population", entry.get("pop", "")),
            "af": entry.get("af", entry.get("allele_freq", "")),
            "ac": entry.get("ac", ""),
            "an": entry.get("an", ""),
            "raw_json": json.dumps(entry, ensure_ascii=False),
        }
        rows.append(row)
    return rows


def _variant_summary_text(v: Variant) -> str:
    # Prefer the class’s own builder, fall back gracefully if it raises
    try:
        return v.return_full_summary()
    except Exception:
        bits: List[str] = []
        if v.rsID:
            bits.append(f"Variant: {v.rsID}")
        if v.organism:
            bits.append(f"Organism: {v.organism}")
        for build, rec in (v.loc_by_build or {}).items():
            chrom = rec.get("chrom", "")
            pos = rec.get("pos", "")
            bits.append(f"{build}: chr{chrom}:{pos}")
        if v.text_summary:
            bits.append(v.text_summary)
        return " | ".join(bits)


# =========================
# Main export
# =========================
def create_files_from_state(state: "State", output_dir: str) -> None:
    """
    Emit UI-friendly artifacts.
    Genes: keep rich per-gene files (as before).
    Variants: CSV-only, centralized tables (no per-variant directories).
    """
    out_root = Path(output_dir).resolve()
    _ensure_dir(out_root)

    # ------------------ Genes (same as before) ------------------
    genes_dir = out_root / "genes"
    _ensure_dir(genes_dir)

    gene_entities: Dict[str, Gene] = (state or {}).get("gene_entities", {}) or {}
    gene_index_rows: List[Dict[str, Any]] = []

    for _, g in sorted(gene_entities.items(), key=lambda kv: (kv[0] or "")):
        if not isinstance(g, Gene):
            continue
        symbol = _safe_name(g.symbol, "UNKNOWN_GENE")
        gdir = genes_dir / symbol
        _ensure_dir(gdir)

        # summary.txt (bulleted for the card body)
        (gdir / "summary.txt").write_text(
            g.summarize_text(include_go=True, include_pathways=True),
            encoding="utf-8",
        )

        # lossless JSON
        (gdir / "gene.json").write_text(
            json.dumps(asdict(g), ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

        # transcripts + isoforms
        _write_tsv(
            gdir / "transcripts.tsv",
            _gene_transcripts_rows(g),
            field_order=[
                "transcript_id",
                "n_exons",
                "name",
                "status",
                "aliases",
                "locations",
                "notes",
                "source",
            ],
        )

        # interactions
        _write_tsv(
            gdir / "interactions_human.tsv",
            _gene_interactions_rows(g, nonhuman=False),
            field_order=["experiment_type", "partner_symbol"],
        )
        _write_tsv(
            gdir / "interactions_nonhuman.tsv",
            _gene_interactions_rows(g, nonhuman=True),
            field_order=["experiment_type", "partner_symbol"],
        )

        # compact index row
        idx = _gene_index_row(g)
        idx.update(
            {
                "summary_relpath": f"genes/{symbol}/summary.txt",
                "json_relpath": f"genes/{symbol}/gene.json",
                "transcripts_relpath": f"genes/{symbol}/transcripts.tsv",
                "interactions_human_relpath": f"genes/{symbol}/interactions_human.tsv",
                "interactions_nonhuman_relpath": f"genes/{symbol}/interactions_nonhuman.tsv",
            }
        )
        gene_index_rows.append(idx)

    _write_csv(
        genes_dir / "genes.index.csv",
        gene_index_rows,
        field_order=[
            "symbol",
            "gene_type",
            "entrez_id",
            "uniprot_id",
            "ensembl_id",
            "chrom",
            "start",
            "end",
            "strand",
            "n_transcripts",
            "has_gwas",
            "gwas_total",
            "n_variants",
            "n_interaction_partners",
            "n_binding_peaks",
            "n_mirna_targets",
            "n_traits",
            "tools_run",
            "summary_relpath",
            "json_relpath",
            "transcripts_relpath",
            "interactions_human_relpath",
            "interactions_nonhuman_relpath",
        ],
    )

    # ------------------ Drugs ------------------
    drugs_dir = out_root / "drugs"
    _ensure_dir(drugs_dir)

    drug_entities: Dict[str, Drug] = (state or {}).get("drug_entities", {}) or {}
    drug_index_rows: List[Dict[str, Any]] = []

    for _, d in sorted(drug_entities.items(), key=lambda kv: (kv[0] or "")):
        if not isinstance(d, Drug):
            continue
        name = _safe_name(d.identifiers.name, "UNKNOWN_DRUG")
        ddir = drugs_dir / name
        _ensure_dir(ddir)

        (ddir / "summary.txt").write_text(d.summarize_text(), encoding="utf-8")
        (ddir / "drug.json").write_text(
            json.dumps(asdict(d), ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

        idx = _drug_index_row(d)
        idx.update(
            {
                "summary_relpath": f"drugs/{name}/summary.txt",
                "json_relpath": f"drugs/{name}/drug.json",
            }
        )
        drug_index_rows.append(idx)

    _write_csv(
        drugs_dir / "drugs.index.csv",
        drug_index_rows,
        field_order=[
            "name",
            "chembl_id",
            "drugcentral_id",
            "cas_number",
            "catalog_number",
            "n_synonyms",
            "n_mentions",
            "tools_run",
            "summary_relpath",
            "json_relpath",
        ],
    )

    # ------------------ Variants (CSV-only) ------------------
    variants_dir = out_root / "variants"
    _ensure_dir(variants_dir)

    variants_csv = variants_dir / "variants.csv"
    locations_csv = variants_dir / "locations.csv"
    per_gene_traits_csv = variants_dir / "per_gene_traits.csv"
    per_gene_context_csv = variants_dir / "per_gene_context.csv"
    functional_predictions_csv = variants_dir / "functional_predictions.csv"
    allele_frequencies_csv = variants_dir / "allele_frequencies.csv"
    summaries_csv = variants_dir / "summaries.csv"

    variant_entities: Dict[str, Variant] = (state or {}).get("variant_entities", {}) or {}

    # Collect rows
    idx_rows: List[Dict[str, Any]] = []
    loc_rows: List[Dict[str, Any]] = []
    pgt_rows: List[Dict[str, Any]] = []
    pgc_rows: List[Dict[str, Any]] = []
    fp_rows: List[Dict[str, Any]] = []
    af_rows: List[Dict[str, Any]] = []
    sm_rows: List[Dict[str, Any]] = []

    for _, v in sorted(variant_entities.items(), key=lambda kv: (kv[0] or "")):
        if not isinstance(v, Variant):
            continue
        var_id = v.rsID or _best_variant_id(v)

        # master row
        idx_rows.append(_variant_index_row(v))

        # satellites
        loc_rows.extend(_variant_locations_rows(var_id, v))
        pgt_rows.extend(_variant_per_gene_traits_rows(var_id, v))
        pgc_rows.extend(_variant_per_gene_context_rows(var_id, v))
        fp_rows.extend(_variant_func_pred_rows(var_id, v))
        af_rows.extend(_variant_af_rows(var_id, v))

        # summary (single consolidated string per variant)
        sm_rows.append({"variant_id": var_id, "summary_text": _variant_summary_text(v)})

    # Write CSVs
    _write_csv(
        variants_csv,
        idx_rows,
        field_order=[
            "variant_id",
            "organism",
            "primary_build",
            "primary_chrom",
            "primary_pos",
            "n_related_genes",
            "n_traits",
            "n_af_freqs",
            "n_func_pred_genes",
            "tools_run",
            "has_text_summary",
        ],
    )
    _write_csv(
        locations_csv,
        loc_rows,
        field_order=["variant_id", "build", "chrom", "pos", "ref", "alt"],
    )
    # For traits, columns vary; we keep flexible plus raw_json
    trait_fields = sorted({k for r in pgt_rows for k in r.keys()} - {"raw_json"})
    _write_csv(
        per_gene_traits_csv,
        pgt_rows,
        field_order=["variant_id", "gene", "trait"] + trait_fields + ["raw_json"],
    )
    _write_csv(
        per_gene_context_csv,
        pgc_rows,
        field_order=["variant_id", "gene", "context", "variant_category"],
    )
    _write_csv(
        functional_predictions_csv,
        fp_rows,
        field_order=["variant_id", "gene", "tool", "score"],
    )
    _write_csv(
        allele_frequencies_csv,
        af_rows,
        field_order=["variant_id", "source", "population", "af", "ac", "an", "raw_json"],
    )
    _write_csv(
        summaries_csv,
        sm_rows,
        field_order=["variant_id", "summary_text"],
    )
