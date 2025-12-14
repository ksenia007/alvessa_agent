from __future__ import annotations

from functools import partial
from typing import Any, Dict, List, Optional

from src.config import DEBUG
from src.state import State
from src.alvessa.domain.gene_class import Gene, GeneGWASProfile, GeneGWASTraitHit
from src.alvessa.domain.gene_components import GeneIdentifiers
from src.alvessa.domain.variant_class import Variant
from .query import query_gene_associations, query_trait_associations
from src.tools.base import Node

GWAS_TOOL_NAME = "gwas_associations_agent"
DB_PATH = "local_dbs"

MODE_CONFIG: Dict[str, Dict[str, Any]] = {
    "summary": {
        "fps_disease_traits": 50,
        "top_studies_by_risk": 100,
        "top_studies_by_significance": 100,
    },
    "extensive": {
        "fps_disease_traits": 200,
        "top_studies_by_risk": 300,
        "top_studies_by_significance": 300,
    },
}


def _debug(message: str) -> None:
    if DEBUG:
        print(f"[GWAS] {message}")


def _safe_float(value: Any, default: float = 1.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _unique(items, *, limit: Optional[int] = None) -> List[str]:
    seen: List[str] = []
    for raw in items or []:
        if raw is None:
            continue
        text = str(raw).strip()
        if not text or text in seen:
            continue
        seen.append(text)
        if limit is not None and len(seen) >= limit:
            break
    return seen


def _normalise_gwas_result(result: Dict[str, Any]) -> Dict[str, Any]:
    """Return a copy of the GWAS result with JSON-serialisable types."""
    cleaned = dict(result)
    linked = cleaned.get("gwas_linked_genes")
    if isinstance(linked, set):
        cleaned["gwas_linked_genes"] = sorted(linked)
    return cleaned


def _collect_variant_annotations(result: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """Merge variant annotations from all summary blocks."""
    annotations: Dict[str, Dict[str, Any]] = {}
    for section in ("summary_by_high_risk_alleles", "summary_by_significance"):
        summary = result.get(section) or {}
        for variant_id, payload in (summary.get("variant_annotations") or {}).items():
            record = annotations.setdefault(
                variant_id,
                {
                    "mapped_gene": payload.get("mapped_gene"),
                    "context": payload.get("context"),
                    "variant_category": payload.get("variant_category"),
                    "associated_disease_trait": {},
                    "total_trait_count": 0,
                    "traits_truncated": False,
                },
            )
            if payload.get("mapped_gene") and not record.get("mapped_gene"):
                record["mapped_gene"] = payload["mapped_gene"]
            if payload.get("context"):
                record["context"] = payload["context"]
            if payload.get("variant_category"):
                record["variant_category"] = payload["variant_category"]
            record["associated_disease_trait"].update(payload.get("associated_disease_trait") or {})
            record["total_trait_count"] = max(
                record.get("total_trait_count", 0),
                int(payload.get("total_trait_count", 0) or len(record["associated_disease_trait"]))
            )
            record["traits_truncated"] = bool(record.get("traits_truncated")) or bool(payload.get("traits_truncated"))
    return annotations


def _build_gene_summary(gene_symbol: str, profile: GeneGWASProfile) -> str:
    if not profile.found or profile.total_associations == 0:
        if profile.p_value_threshold:
            return (
                f"No genome-wide significant GWAS associations found for {gene_symbol} "
                f"(p<{profile.p_value_threshold})."
            )
        return f"No genome-wide significant GWAS associations found for {gene_symbol}."

    header = (
        f"*GWAS: Gene {gene_symbol} has {profile.total_associations} GWAS associations "
        f"({profile.significant_associations} significant) across {profile.total_studies} studies."
    )

    if not profile.top_traits:
        return header

    highlights = []
    for hit in profile.top_traits[:5]:
        p_str = hit.p_value if hit.p_value is not None else "N/A"
        risk_str = f", risk={hit.risk_score}" if hit.risk_score else ""
        highlights.append(f"{hit.trait} (rsID={hit.rsid}, p={p_str}{risk_str})")

    return header + " Top associations: " + "; ".join(highlights) + "."


def _apply_gene_result(
    gene: Gene,
    result: Dict[str, Any],
    variant_store: Dict[str, Variant],
) -> GeneGWASProfile:
    """Populate gene/variant objects and return the compact GWAS profile."""
    gene.add_tool(GWAS_TOOL_NAME)

    profile = GeneGWASProfile(
        found=bool(result.get("found", False)),
        total_associations=int(result.get("total_associations", 0) or 0),
        significant_associations=int(result.get("total_significant_associations", 0) or 0),
        total_studies=int(result.get("total_studies_analyzed", 0) or 0),
        p_value_threshold=(str(result.get("p_value_threshold")) if result.get("p_value_threshold") is not None else None),
    )

    if not profile.found or profile.total_associations == 0:
        gene.set_gwas_profile(profile)
        summary = _build_gene_summary(gene.symbol, profile)
        if summary:
            gene.update_text_summaries(summary)
        return profile

    annotations = _collect_variant_annotations(result)
    variant_count = 0
    trait_link_count = 0
    best_trait_hits: Dict[str, tuple] = {}
    variant_ranks: List[tuple] = []

    for variant_id, payload in annotations.items():
        trait_map = payload.get("associated_disease_trait") or {}
        if not trait_map:
            continue

        variant = variant_store.get(variant_id)
        if variant is None:
            variant = Variant(rsID=variant_id, organism="human")
            variant_store[variant_id] = variant

        variant.add_per_gene_traits(gene.symbol, trait_map)
        variant.add_per_gene_context(
            gene.symbol,
            payload.get("context", ""),
            payload.get("variant_category", ""),
        )

        if payload.get("traits_truncated"):
            total = payload.get("total_trait_count", len(trait_map))
            variant.update_text_summaries(
                f"*GWAS: Associations truncated to top 20 of {total} traits for {gene.symbol}."
            )

        gene.link_variant(variant)
        for trait_name in trait_map.keys():
            gene.link_trait(trait_name)

        variant_count += 1
        trait_link_count += len(trait_map)

        min_p_for_variant = float("inf")
        for trait, stats in trait_map.items():
            p_float = _safe_float(stats.get("p_value"), default=float("inf"))
            display_p = None if stats.get("p_value") is None else str(stats.get("p_value"))
            risk_value = None if stats.get("risk_score") is None else str(stats.get("risk_score"))
            hit = GeneGWASTraitHit(
                trait=trait,
                rsid=variant_id,
                p_value=display_p,
                risk_score=risk_value,
            )
            existing = best_trait_hits.get(trait)
            if existing is None or p_float < existing[0]:
                best_trait_hits[trait] = (p_float, hit)
            if p_float < min_p_for_variant:
                min_p_for_variant = p_float

        variant_ranks.append((min_p_for_variant, variant_id))

    top_trait_hits = [item[1] for item in sorted(best_trait_hits.values(), key=lambda x: x[0])[:10]]
    variant_ranks.sort(key=lambda x: x[0])
    top_variants = [vid for _, vid in variant_ranks[:10] if vid]

    profile.variant_count = variant_count
    profile.trait_link_count = trait_link_count
    profile.top_traits = top_trait_hits
    profile.top_variants = top_variants
    high_related = _unique((result.get("summary_by_high_risk_alleles", {}) or {}).get("related_genes") or [])
    sig_related = _unique((result.get("summary_by_significance", {}) or {}).get("related_genes") or [])
    high_set = set(high_related)
    profile.related_genes = _unique(high_related + [g for g in sig_related if g not in high_set], limit=10)
    profile.affected_proteins = _unique(
        ((result.get("summary_by_high_risk_alleles", {}) or {}).get("affected_protein_levels") or [])
        + ((result.get("summary_by_significance", {}) or {}).get("affected_protein_levels") or []),
        limit=10,
    )

    gene.set_gwas_profile(profile)
    summary = _build_gene_summary(gene.symbol, profile)
    if summary:
        gene.update_text_summaries(summary)

    return profile


def _ensure_gene_entities(state: State) -> Dict[str, Gene]:
    gene_objs: Dict[str, Gene] = dict(state.get("gene_entities", {}) or {})
    # Do NOT add new genes here; rely on existing gene_entities
    return gene_objs


def gwas_associations_agent(state: "State", mode: str = "summary") -> "State":
    """Query GWAS catalogue for each gene and enrich the state object model."""
    mode_key = mode.lower().strip()
    query_kwargs = MODE_CONFIG.get(mode_key, MODE_CONFIG["summary"])

    gene_objs = _ensure_gene_entities(state)
    variant_objs: Dict[str, Variant] = dict(state.get("variant_entities", {}) or {})
    results: Dict[str, Any] = dict(state.get("gwas_associations", {}) or {})

    for gene_symbol, gene_obj in gene_objs.items():
        if gene_obj.has_tool(GWAS_TOOL_NAME):
            _debug(f"Skipping {gene_symbol}: GWAS data already collected")
            continue

        _debug(f"Querying GWAS associations for {gene_symbol} (mode={mode_key})")
        try:
            result = query_gene_associations(
                gene_symbol=gene_symbol,
                db_path=DB_PATH,
                **query_kwargs,
            )
        except Exception as exc:
            _debug(f"Error querying {gene_symbol}: {exc}")
            continue

        if not result:
            continue

        results[gene_symbol] = _normalise_gwas_result(result)
        profile = _apply_gene_result(gene_obj, result, variant_objs)
        _debug(
            f"{gene_symbol}: {profile.total_associations} associations, "
            f"{profile.variant_count} variants, {profile.trait_link_count} trait links"
        )

    updates: Dict[str, Any] = {
        "gene_entities": gene_objs,
        "variant_entities": variant_objs,
    }

    if results:
        updates["gwas_associations"] = results
    return updates


def _create_empty_association_record(identifier: str, is_gene: bool = True) -> Dict[str, Any]:
    record: Dict[str, Any] = {
        "found": False,
        "total_associations": 0,
        "total_significant_associations": 0,
        "total_studies_analyzed": 0,
        "summary_by_high_risk_alleles": {},
        "summary_by_significance": {},
        "variant_annotations": {},
    }
    record["gene" if is_gene else "trait_term"] = identifier
    return record


def _select_trait_term(state: State) -> str:
    extracted = state.get("traits", []) or []
    if extracted:
        return extracted[0]
    messages = state.get("messages", []) or []
    return messages[-1]["content"] if messages else ""


def query_by_trait_agent(state: "State") -> "State":
    """Query GWAS catalogue for a trait/disease term."""
    trait_term = _select_trait_term(state)
    if not trait_term:
        return {"trait_associations": _create_empty_association_record("", is_gene=False), "genes": []}

    _debug(f"Querying trait associations for '{trait_term}'")
    try:
        result = query_trait_associations(
            trait_term=trait_term,
            db_path=DB_PATH,
            fps_genes=60,
            exact_match=False,
        )
    except Exception as exc:
        _debug(f"Trait query failed for '{trait_term}': {exc}")
        payload = _create_empty_association_record(trait_term, is_gene=False)
        payload["error"] = str(exc)
        return {"trait_associations": payload, "genes": []}

    if not result:
        return {"trait_associations": _create_empty_association_record(trait_term, is_gene=False), "genes": []}

    normalised = _normalise_gwas_result(result)
    payload = {
        "trait_term": trait_term,
        "found": normalised.get("found", False),
        "total_associations": normalised.get("total_associations", 0),
        "total_significant_associations": normalised.get("total_significant_associations", 0),
        "total_studies_analyzed": normalised.get("total_studies_analyzed", 0),
        "summary_by_high_risk_alleles": normalised.get("summary_by_high_risk_alleles", {}),
        "summary_by_significance": normalised.get("summary_by_significance", {}),
        "variant_annotations": _collect_variant_annotations(normalised),
    }

    # Do not introduce new genes into state; only return associations
    return {"trait_associations": payload, "genes": []}


def _extract_genes_from_trait_result(result: Dict[str, Any]) -> List[str]:
    if not result.get("found", False):
        return []

    gene_stats: Dict[str, Dict[str, Any]] = {}
    for key in ("studies_by_high_risk_alleles", "studies_by_significance"):
        for study in result.get(key, []) or []:
            max_risk = study.get("max_risk", 0)
            sig_count = study.get("sig_count", 0)
            for gene in study.get("related_genes", []) or []:
                stats = gene_stats.setdefault(gene, {"studies": 0, "sig_hits": 0, "max_risk": 0.0})
                stats["studies"] += 1
                stats["sig_hits"] += sig_count
                stats["max_risk"] = max(stats["max_risk"], max_risk)

    if not gene_stats:
        return []

    # prioritise by study count, then significant hits, then risk
    ranked = sorted(
        gene_stats.items(),
        key=lambda item: (item[1]["studies"], item[1]["sig_hits"], item[1]["max_risk"]),
        reverse=True,
    )
    return [gene for gene, _ in ranked[:20]]


NODES: tuple[Node, ...] = (
    Node(
        name="query_gwas_by_gene",
        entry_point=gwas_associations_agent,
        description=(
            "Retrieves genome-wide association study (GWAS) results for a given gene. It collects traits and diseases associated with genetic variants linked to that gene, along with the specific variants"
        ),
    ),
    Node(
        name="query_gwas_extensive",
        entry_point=partial(gwas_associations_agent, mode="extensive"),
        description=(
            "This is a more comprehensive version of the query_gwas_by_gene tool, and it is used to retrieve more detailed information about the GWAS results. It collects an extensive list of traits/diseases associated with an extensive list of genetic variants linked to that gene. Use this tool *ONLY* if the question is very specific that requires what is equivalent to an extensive database search, not to general characterisation of the gene."
        ),
    ),
    Node(
        name="expand_gene_set_by_trait",
        entry_point=query_by_trait_agent,
        description=(
            "Query GWAS traits by disease/phenotype keyword to expand the working gene list. This could be used to discover more relevant genes underlying a trait, if a broader genetic context would be valuable for the analysis."
        ),
    ),
)
