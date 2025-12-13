"""LangGraph node that surfaces drug-induced gene expression changes from CIGS."""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

try:  # Import DEBUG lazily so missing API keys do not break scripting.
    from src.config import DEBUG
except Exception:  # pragma: no cover - dev environments without env vars
    DEBUG = False

from src.state import State
from src.tools.base import Node
from src.alvessa.domain.drug_class import Drug

MODULE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = MODULE_ROOT.parents[2]
LOCAL_DBS = REPO_ROOT / "local_dbs"
MEDCHEM_LIBRARY = LOCAL_DBS / "compound_library_medchemexpress.csv"
CIGS_ROOT = LOCAL_DBS / "deseq" / "cigs" / "drug_centric"

PADJ_THRESHOLD: float = 0.05
TOOL_NAME = "drug_gene_interactions"

CELL_LINES: Dict[str, Dict[str, str]] = {
    "HEK293T": {"folder": "HEK293T", "prefix": "MCE1_293T_24H"},
    "MDA_MB_231": {"folder": "MDA_MB_231", "prefix": "MCE1_231_24H"},
}


@dataclass(frozen=True)
class MedChemEntry:
    catalog_number: str
    product_name: str
    synonyms: Tuple[str, ...]

    def iter_terms(self) -> Iterable[str]:
        yield self.product_name
        yield self.catalog_number
        for syn in self.synonyms:
            yield syn


@dataclass
class DrugEntityView:
    key: str
    name: Optional[str]
    catalog_number: Optional[str]
    synonyms: Tuple[str, ...]
    obj: Optional[Drug]

    def aliases(self) -> Iterable[str]:
        if self.name:
            yield self.name
        for synonym in self.synonyms:
            yield synonym


_NORMALIZE_RE = re.compile(r"[^a-z0-9]+")


def _normalize_term(value: Optional[str]) -> str:
    if not value:
        return ""
    cleaned = value.strip().lower()
    if not cleaned:
        return ""
    return _NORMALIZE_RE.sub("", cleaned)


def _split_synonyms(raw: Optional[str]) -> Tuple[str, ...]:
    if not raw:
        return ()
    pieces = re.split(r"[;\n]", raw)
    uniq: List[str] = []
    for piece in pieces:
        token = piece.strip()
        if token and token not in uniq:
            uniq.append(token)
    return tuple(uniq)


@lru_cache(maxsize=1)
def _load_medchem_library() -> Tuple[Dict[str, MedChemEntry], Dict[str, MedChemEntry]]:
    """Return (name_index, catalog_index) for the MedChemExpress catalog."""
    name_index: Dict[str, MedChemEntry] = {}
    catalog_index: Dict[str, MedChemEntry] = {}

    if not MEDCHEM_LIBRARY.exists():
        if DEBUG:
            print(f"[drug_gene_interactions] Missing MedChemExpress library at {MEDCHEM_LIBRARY}")
        return name_index, catalog_index

    with MEDCHEM_LIBRARY.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            product_name = (row.get("Product Name") or "").strip()
            catalog_number = (row.get("Catalog Number") or "").strip()
            if not product_name or not catalog_number:
                continue
            synonyms = _split_synonyms(row.get("Synonyms"))
            entry = MedChemEntry(
                catalog_number=catalog_number,
                product_name=product_name,
                synonyms=synonyms,
            )
            catalog_key = _normalize_term(catalog_number)
            if catalog_key:
                catalog_index[catalog_key] = entry

            for term in entry.iter_terms():
                norm = _normalize_term(term)
                if not norm or norm in name_index:
                    continue
                name_index[norm] = entry

    return name_index, catalog_index


def _coerce_drug_entity_view(key: str, value: Any) -> DrugEntityView:
    if isinstance(value, Drug):
        ids = value.identifiers
        synonyms = tuple(ids.synonyms or [])
        return DrugEntityView(
            key=key,
            name=ids.name,
            catalog_number=ids.catalog_number,
            synonyms=synonyms,
            obj=value,
        )

    if isinstance(value, dict):
        payload = value.get("identifiers") if isinstance(value.get("identifiers"), dict) else value
        name = payload.get("name") or payload.get("Product Name")
        catalog_number = payload.get("catalog_number") or payload.get("Catalog Number")
        raw_synonyms = payload.get("synonyms")
        if isinstance(raw_synonyms, list):
            synonyms = tuple(x for x in raw_synonyms if isinstance(x, str))
        elif isinstance(raw_synonyms, str):
            synonyms = _split_synonyms(raw_synonyms)
        else:
            synonyms = ()
        return DrugEntityView(
            key=key,
            name=name,
            catalog_number=catalog_number,
            synonyms=synonyms,
            obj=None,
        )

    return DrugEntityView(key=key, name=str(value), catalog_number=None, synonyms=(), obj=None)


def _build_entity_index(drug_entities: Dict[str, Any]) -> Tuple[Dict[str, DrugEntityView], Dict[str, DrugEntityView]]:
    index: Dict[str, DrugEntityView] = {}
    normalized_views: Dict[str, DrugEntityView] = {}

    for key, value in (drug_entities or {}).items():
        view = _coerce_drug_entity_view(key, value)
        normalized_views[key] = view
        for alias in view.aliases():
            norm = _normalize_term(alias)
            if norm and norm not in index:
                index[norm] = view
        if view.catalog_number:
            cat_norm = _normalize_term(view.catalog_number)
            if cat_norm and cat_norm not in index:
                index[cat_norm] = view

    return index, normalized_views


def _collect_drug_queries(
    state: State,
    entity_index: Dict[str, DrugEntityView],
    normalized_views: Dict[str, DrugEntityView],
) -> List[Tuple[str, Optional[DrugEntityView]]]:
    seen: set[str] = set()
    queries: List[Tuple[str, Optional[DrugEntityView]]] = []

    for raw in state.get("drugs") or []:
        label = (raw or "").strip()
        if not label:
            continue
        sentinel = label.lower()
        if sentinel in seen:
            continue
        seen.add(sentinel)
        view = entity_index.get(_normalize_term(label))
        queries.append((label, view))

    if queries:
        return queries

    # Fallback: process known drug entities even if the top-level list is empty
    for view in normalized_views.values():
        label = view.name or (view.synonyms[0] if view.synonyms else None)
        if not label:
            continue
        sentinel = label.lower()
        if sentinel in seen:
            continue
        seen.add(sentinel)
        queries.append((label, view))

    return queries


def _lookup_medchem_entry(
    term: Optional[str],
    name_index: Dict[str, MedChemEntry],
    catalog_index: Dict[str, MedChemEntry],
) -> Optional[MedChemEntry]:
    norm = _normalize_term(term)
    if not norm:
        return None
    return catalog_index.get(norm) or name_index.get(norm)


def _resolve_entry_for_drug(
    query: str,
    view: Optional[DrugEntityView],
    name_index: Dict[str, MedChemEntry],
    catalog_index: Dict[str, MedChemEntry],
) -> Tuple[Optional[MedChemEntry], Optional[str]]:
    """Return (library entry, fallback_catalog_number)."""
    tried: set[str] = set()
    candidates: List[str] = []

    if view and view.catalog_number:
        candidates.append(view.catalog_number)
    candidates.append(query)
    if view:
        if view.name:
            candidates.append(view.name)
        candidates.extend(view.synonyms)

    for candidate in candidates:
        norm = _normalize_term(candidate)
        if not norm or norm in tried:
            continue
        tried.add(norm)
        entry = _lookup_medchem_entry(candidate, name_index, catalog_index)
        if entry:
            return entry, None

    fallback = view.catalog_number if view and view.catalog_number else None
    return None, fallback


def _catalog_to_filename_fragment(catalog_number: str) -> str:
    return catalog_number.strip().upper().replace("-", "_")


def _resolve_deg_file(cell_line: str, fragment: str) -> Optional[Path]:
    cfg = CELL_LINES[cell_line]
    folder = CIGS_ROOT / cfg["folder"]
    expected = folder / f"{cfg['prefix']}_{fragment}_DESeq2_DEG.csv"
    if expected.exists():
        return expected

    # Fallback: glob for close matches (case variations, extra underscores)
    pattern = f"{cfg['prefix']}_{fragment}*_DESeq2_DEG.csv"
    for candidate in sorted(folder.glob(pattern)):
        if candidate.name.endswith("_DESeq2_DEG.csv"):
            return candidate
    return None


def _safe_float(value: Any) -> Optional[float]:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(parsed):
        return None
    return parsed


def _load_significant_genes(csv_path: Path) -> List[Dict[str, Any]]:
    hits: List[Dict[str, Any]] = []
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            padj = _safe_float(row.get("padj"))
            if padj is None or padj >= PADJ_THRESHOLD:
                continue
            gene = row.get("gene") or row.get("Gene")
            if not gene:
                continue
            hits.append(
                {
                    "gene": gene,
                    "padj": padj,
                    "pvalue": _safe_float(row.get("pvalue")),
                    "log2_fold_change": _safe_float(row.get("log2FoldChange")),
                    "base_mean": _safe_float(row.get("baseMean")),
                }
            )
    hits.sort(key=lambda rec: rec["padj"])
    return hits


def _relative_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _update_drug_metadata(view: Optional[DrugEntityView], payload: Dict[str, Any]) -> None:
    if not view or not view.obj:
        return
    view.obj.add_tool("CIGSDrugGeneInteractions")
    metadata = view.obj.metadata.setdefault("cigs_drug_gene_interactions", {})
    metadata.update({k: v for k, v in payload.items() if k in {"medchem_id", "drug_name"}})
    cells = metadata.setdefault("cell_lines", {})
    for cell_line, info in payload.get("cell_lines", {}).items():
        cells[cell_line] = {
            "significant_gene_count": info["significant_gene_count"],
            "file": info["file"],
        }


def drug_gene_interactions_node(state: "State") -> "State":
    name_index, catalog_index = _load_medchem_library()
    if not name_index and not catalog_index:
        note = "MedChemExpress library not available; cannot resolve drug IDs for CIGS lookup."
        notes = state.get("text_notes") or []
        notes.append(note)
        state["text_notes"] = notes
        return state

    drug_entities = state.get("drug_entities") or {}
    entity_index, normalized_views = _build_entity_index(drug_entities)
    queries = _collect_drug_queries(state, entity_index, normalized_views)
    if not queries:
        if DEBUG:
            print("[drug_gene_interactions] No drug names available in state.")
        return state

    results: Dict[str, Any] = dict(state.get("drug_gene_interactions") or {})
    processed_keys: set[str] = set()
    narrative: List[str] = []

    for query, view in queries:
        entry, fallback_catalog = _resolve_entry_for_drug(query, view, name_index, catalog_index)
        medchem_id = entry.catalog_number if entry else fallback_catalog
        if not medchem_id:
            message = f"Drug '{query}' not found in MedChemExpress library; skipping CIGS lookup."
            narrative.append(message)
            continue

        dedup_key = medchem_id.upper()
        if dedup_key in processed_keys:
            continue
        processed_keys.add(dedup_key)

        fragment = _catalog_to_filename_fragment(medchem_id)
        drug_name = entry.product_name if entry else (view.name if view and view.name else query)
        cell_payload: Dict[str, Any] = {}
        cell_summaries: List[str] = []

        for cell_line in CELL_LINES:
            deg_path = _resolve_deg_file(cell_line, fragment)
            if not deg_path:
                cell_payload[cell_line] = {
                    "file": None,
                    "significant_gene_count": 0,
                    "genes": [],
                    "status": f"No DESeq2 file for {cell_line} and {fragment}.",
                }
                continue

            genes = _load_significant_genes(deg_path)
            summary = (
                f"{len(genes)} genes with padj<{PADJ_THRESHOLD}"
                if genes
                else "No significantly perturbed genes (padj<0.05)"
            )
            if genes:
                cell_summaries.append(f"{cell_line}: {len(genes)} hits")
            else:
                cell_summaries.append(f"{cell_line}: no significant hits")

            cell_payload[cell_line] = {
                "file": _relative_path(deg_path),
                "significant_gene_count": len(genes),
                "genes": genes,
                "status": summary,
            }

        result_key = dedup_key
        results[result_key] = {
            "query": query,
            "drug_name": drug_name,
            "medchem_id": medchem_id,
            "cell_lines": cell_payload,
            "source": "CIGS drug-centric DESeq2",
        }

        _update_drug_metadata(view, results[result_key])

        joined = ", ".join(cell_summaries)
        narrative.append(f"CIGS: {drug_name} ({medchem_id}) → {joined}.")

    if narrative:
        notes = state.get("text_notes") or []
        notes.extend(narrative)
        state["text_notes"] = notes

    if results:
        state["drug_gene_interactions"] = results

    used = state.get("used_tools") or []
    used.append(TOOL_NAME)
    state["used_tools"] = used
    return state


NODE = Node(
    name=TOOL_NAME,
    entry_point=drug_gene_interactions_node,
    description=(
        "Given drug names, looks up MedChemExpress IDs and reports genes perturbed in the data collected from Chemically Induced Gene Sets (CIGS; Nature Methods, 2025) for two cell lines: "
        "HEK293T and MDA-MB-231, using DESeq2 analysis (padj<0.05)."
    ),
)

NODES: Tuple[Node, ...] = (NODE,)
