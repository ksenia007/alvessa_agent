"""LangGraph node that surfaces gene-centric CIGS drug perturbations."""

from __future__ import annotations

import csv
import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

try:  # Avoid hard failure if config imports require missing env vars
    from src.config import DEBUG
except Exception:  # pragma: no cover - fallback for sandboxed tests
    DEBUG = False

from src.state import State
from src.tools.base import Node
from src.alvessa.domain.gene_class import Gene

MODULE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = MODULE_ROOT.parents[2]
GENE_CENTRIC_ROOT = REPO_ROOT / "local_dbs" / "deseq" / "cigs" / "gene_centric"
MEDCHEM_LIBRARY = REPO_ROOT / "local_dbs" / "compound_library_medchemexpress.csv"

CELL_LINES: Tuple[str, ...] = ("HEK293T", "MDA_MB_231")
TOOL_NAME = "gene_drug_interactions"
SANITIZE_PATTERN = re.compile(r"[^A-Z0-9_.-]")


def _normalize_catalog_number(value: str) -> str:
    return re.sub(r"[^A-Z0-9]+", "", (value or "").upper())


def _sanitize_gene_symbol(symbol: str) -> str:
    cleaned = (symbol or "").strip().upper()
    return SANITIZE_PATTERN.sub("_", cleaned)


def _safe_float(value: Any) -> Optional[float]:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    if parsed != parsed:  # NaN
        return None
    return parsed


def _collect_gene_queries(state: State) -> List[Tuple[str, Optional[Gene]]]:
    gene_entities = state.get("gene_entities") or {}
    entity_index: Dict[str, Gene] = {}
    for key, gene in gene_entities.items():
        if isinstance(gene, Gene):
            entity_index[key.upper()] = gene

    queries: List[Tuple[str, Optional[Gene]]] = []
    seen: set[str] = set()

    for gene_name in state.get("genes") or []:
        symbol = (gene_name or "").strip()
        if not symbol:
            continue
        key = symbol.upper()
        if key in seen:
            continue
        seen.add(key)
        queries.append((symbol, entity_index.get(key)))

    if queries:
        return queries

    for key, gene in entity_index.items():
        if key in seen:
            continue
        seen.add(key)
        display = gene.symbol or key
        queries.append((display, gene))

    return queries


def _load_gene_cell_line_data(cell_line: str, sanitized_gene: str) -> List[Dict[str, Any]]:
    gene_file = GENE_CENTRIC_ROOT / cell_line / "genes" / f"{sanitized_gene}.csv"
    if not gene_file.exists():
        return []

    records: List[Dict[str, Any]] = []
    with gene_file.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            records.append(
                {
                    "drug_id": row.get("drug_id", ""),
                    "direction": row.get("direction", "unknown"),
                    "log2_fold_change": _safe_float(row.get("log2_fold_change")),
                    "pvalue": _safe_float(row.get("pvalue")),
                    "padj": _safe_float(row.get("padj")),
                    "source_file": row.get("source_file", ""),
                }
            )
    records.sort(key=lambda rec: rec["padj"] if rec["padj"] is not None else 1.0)
    return records


@lru_cache(maxsize=1)
def _load_medchem_lookup() -> Dict[str, str]:
    lookup: Dict[str, str] = {}
    if not MEDCHEM_LIBRARY.exists():
        if DEBUG:
            print(f"[gene_drug_interactions] Missing MedChemExpress library at {MEDCHEM_LIBRARY}")
        return lookup

    with MEDCHEM_LIBRARY.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            catalog_number = (row.get("Catalog Number") or "").strip()
            product_name = (row.get("Product Name") or "").strip()
            if not catalog_number or not product_name:
                continue
            lookup[_normalize_catalog_number(catalog_number)] = product_name
    return lookup


def _attach_drug_names(records: List[Dict[str, Any]]) -> None:
    medchem = _load_medchem_lookup()
    for rec in records:
        raw_id = (rec.get("drug_id") or "").strip()
        normalized = _normalize_catalog_number(raw_id.replace("_", "-"))
        rec["drug_name"] = medchem.get(normalized, raw_id)


def _summarize_hits(gene_symbol: str, cell_map: Dict[str, Dict[str, Any]]) -> str:
    pieces = []
    for cell_line, payload in cell_map.items():
        count = len(payload.get("drugs") or [])
        if count == 0:
            continue
        pieces.append(f"{cell_line}: {count} drugs")
    if not pieces:
        return f"No significant CIGS perturbations found for {gene_symbol}."
    joined = ", ".join(pieces)
    return f"CIGS gene-centric hits for {gene_symbol}: {joined}."


def gene_drug_interactions_node(state: State) -> State:
    if not GENE_CENTRIC_ROOT.exists():
        note = "CIGS gene-centric directory is missing; cannot look up gene→drug interactions."
        notes = state.get("text_notes") or []
        notes.append(note)
        state["text_notes"] = notes
        return state

    queries = _collect_gene_queries(state)
    if not queries:
        if DEBUG:
            print("[gene_drug_interactions] No genes available in state.")
        return state

    results: Dict[str, Any] = dict(state.get("gene_drug_interactions") or {})
    text_notes = state.get("text_notes") or []

    for raw_query, gene_obj in queries:
        canonical = (raw_query or "").strip()
        if not canonical:
            continue
        gene_key = canonical.upper()
        sanitized = _sanitize_gene_symbol(gene_key)
        cell_payload: Dict[str, Dict[str, Any]] = {}
        total_hits = 0

        for cell_line in CELL_LINES:
            drugs = _load_gene_cell_line_data(cell_line, sanitized)
            if drugs:
                _attach_drug_names(drugs)
                total_hits += len(drugs)
                status = f"{len(drugs)} drugs with padj<0.05"
            else:
                status = "No significant perturbations"
            cell_payload[cell_line] = {"drugs": drugs, "status": status}

        results[gene_key] = {
            "query": raw_query,
            "gene": gene_key,
            "cell_lines": cell_payload,
            "source": "CIGS gene-centric DESeq2",
            "total_hits": total_hits,
        }

        if gene_obj and isinstance(gene_obj, Gene):
            gene_obj.add_tool("CIGSGeneDrugInteractions")
            summary = _summarize_hits(gene_obj.symbol or gene_key, cell_payload)
            gene_obj.update_text_summaries(summary)

        text_notes.append(_summarize_hits(gene_key, cell_payload))

    state["gene_drug_interactions"] = results
    state["text_notes"] = text_notes

    used = state.get("used_tools") or []
    used.append(TOOL_NAME)
    state["used_tools"] = used
    return state


NODE = Node(
    name=TOOL_NAME,
    entry_point=gene_drug_interactions_node,
    description=(
        "For each gene, returns the list of drugs that significantly perturb it (padj<0.05) in the "
        "CIGS HEK293T and MDA-MB-231 datasets, including direction and fold-change."
    ),
)

NODES: Tuple[Node, ...] = (NODE,)
