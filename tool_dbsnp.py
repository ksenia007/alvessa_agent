from __future__ import annotations

import statistics
from typing import Any, Dict, Iterable, List, Set

from config import DEBUG, DBSNP_DEFAULT_ASSEMBLY
from state import State
from tools.dbsnp.query import get_variant_info


def _debug(message: str) -> None:
    if DEBUG:
        print(f"[dbSNP] {message}")


def _create_coordinate_summary(rsid: str, coordinates: List[Dict[str, Any]]) -> str:
    """Generates a concise, human-readable summary of variant coordinates."""
    if not coordinates:
        return f"{rsid} has no coordinate data available"

    primary_coord = next(
        (c for c in coordinates if "GRCh38" in str(c.get("assembly", ""))),
        coordinates[0],
    )
    chrom = primary_coord.get("chrom", "Unknown")
    pos38 = primary_coord.get("pos", "Unknown")

    return f"{rsid} is located at chr{chrom}:{pos38} (GRCh38)"


def _create_af_summary_sorted_by_counts(rows, target_alt=None, min_people=10_000):
    """
    ALT-only AF summary across studies (3-decimal AFs, no percents).
    Adds heterogeneity (I²) and shows extreme 2 by AF (low/high) + largest-N examples.
    """
    # 1) Aggregate ALT rows per study
    per = {}  # study -> (K_alt, N)
    for r in rows or []:
        study = (r.get("study_name") or r.get("study") or "").strip()
        obs = r.get("observation") or {}
        ref = str(obs.get("deleted_sequence") or "").strip().upper()
        alt = str(obs.get("inserted_sequence") or "").strip().upper()
        if not study or not ref or not alt or alt == ref:
            continue  # drop REF→REF and incomplete
        if target_alt and alt != target_alt.upper():
            continue

        n_raw = r.get("total_count") or r.get("allele_number") or r.get("allele_denominator")
        k_raw = r.get("allele_count")
        if k_raw < min_people*2: 
            continue
        af_raw = r.get("allele_frequency")
        try:
            n = int(float(n_raw))
            k = int(round(float(k_raw))) if k_raw is not None else int(round(float(af_raw) * n))
        except (TypeError, ValueError):
            continue
        if n <= 0 or not (0 <= k <= n):
            continue

        K, N = per.get(study, (0, 0))
        per[study] = (K + k, N + n)

    items = [{"study": s, "k": K, "n": N, "p": K / N} for s, (K, N) in per.items() if N > 0]
    if not items:
        return "No ALT-allele observations match the criteria."

    # 2) Core stats
    items.sort(key=lambda x: x["n"])
    tot_k = sum(x["k"] for x in items)
    tot_n = sum(x["n"] for x in items)
    pooled = tot_k / tot_n
    med_all = statistics.median(x["p"] for x in items)

    # high-powered: ≈ 10k people -> ~20k alleles
    allele_cut = 2 * int(min_people or 0)
    hi = [x for x in items if x["n"] >= allele_cut]
    med_hi = statistics.median(x["p"] for x in hi) if hi else None

    # 3) Agreement (fixed-effect I² with plus-four variance)
    def _var_plus4(p, n):
        n4 = n + 4
        p4 = (p * n + 2) / n4
        return p4 * (1 - p4) / n4

    ws = [1.0 / max(_var_plus4(x["p"], x["n"]), 1e-12) for x in items]
    p_fe = sum(w * x["p"] for w, x in zip(ws, items)) / sum(ws)
    Q = sum(w * (x["p"] - p_fe) ** 2 for w, x in zip(ws, items))
    df = max(1, len(items) - 1)
    I2 = max(0.0, (Q - df) / max(Q, 1e-12)) * 100.0  # percent (conventional)

    # 4) Extremes by AF (low/high) and by N (largest two)
    by_af = sorted(items, key=lambda x: x["p"])
    low2 = by_af[:2]
    high2 = by_af[-2:] if len(by_af) >= 2 else by_af[-1:]

    by_n = sorted(items, key=lambda x: x["n"])
    topN2 = by_n[-2:] if len(by_n) >= 2 else by_n[-1:]

    f = lambda x: f"{x:.3f}"  # 3-decimal AFs
    label = f"ALT {target_alt.upper()}" if target_alt else "ALT (any substitution)"

    parts = [
        f"{label}: weighted AF {f(pooled)} across {len(items)} studies, median {f(med_all)}",
    ]
    if med_hi is not None:
        parts.append(f"median among large cohorts (n>={min_people:,} people) {f(med_hi)}")
    parts.append(f"heterogeneity I^2={I2:.1f}%")
    parts.append(
        "Low AF examples: "
        + ", ".join(f"{x['study']} {f(x['p'])} (n={x['n']:,})" for x in low2)
    )
    parts.append(
        "High AF examples: "
        + ", ".join(f"{x['study']} {f(x['p'])} (n={x['n']:,})" for x in reversed(high2))
    )
    parts.append(
        "Largest cohorts: "
        + ", ".join(f"{x['study']} {f(x['p'])} (n={x['n']:,})" for x in reversed(topN2))
    )
    return "; ".join(parts) + "."


def _format_population_summary(rows: Iterable[Dict[str, Any]]) -> str:
    if not rows:
        return "No population frequency data available."
    return _create_af_summary_sorted_by_counts(rows)


def _to_list(value: Any) -> List[str]:
    if not value:
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(v) for v in value if v]
    return [str(value)]


def dbsnp_variants_agent(
    state: "State",
    assembly: str = None,
    include_population_summaries: bool = False,
) -> "State":
    """Annotate existing Variant entities with dbSNP coordinates and population data."""

    assembly = assembly or DBSNP_DEFAULT_ASSEMBLY
    variant_objs: Dict[str, Any] = dict(state.get("variant_entities", {}) or {})

    if not variant_objs:
        _debug("No variants found in state; skipping dbSNP lookup.")
        return {}

    _debug(f"Using genome assembly: {assembly}")
    _debug(f"Processing {len(variant_objs)} rsIDs from prior GWAS calls.")

    for rsid, variant in variant_objs.items():
        try:
            result = get_variant_info(rsid, assembly)
        except Exception as exc:
            _debug(f"Error fetching {rsid}: {exc}")
            continue

        if not result:
            continue

        coordinates = result.get("coordinates") or []
        if coordinates:
            for coord in coordinates:
                variant.add_location(
                    build=coord.get("assembly", "Unknown"),
                    chrom=coord.get("chrom", "Unknown"),
                    pos=coord.get("pos", 0),
                    ref=_to_list(coord.get("ref")),
                    alt=_to_list(coord.get("alt")),
                )
        else:
            _debug(f"No coordinates returned for {rsid}; keeping GWAS-only data.")

        allele_freqs = result.get("allele_frequencies") or []
        if allele_freqs:
            variant.add_af_freq(allele_freqs)

        variant.add_tool("dbsnp_variants_agent")

        summaries: List[str] = []
        if coordinates:
            summaries.append(_create_coordinate_summary(rsid, coordinates))
        if include_population_summaries and allele_freqs:
            summaries.append(_format_population_summary(allele_freqs))

        for summary in summaries:
            variant.update_text_summaries(summary)

    return {"variant_entities": variant_objs}


def _extract_rsids_from_text(text: str) -> Set[str]:
    """Extracts rsIDs from text using regex pattern matching."""
    import re

    pattern = r"rs\d+"
    matches = re.findall(pattern, text, re.IGNORECASE)
    return set(match.lower() for match in matches)


def has_dbsnp_variants(state: "State") -> bool:
    variants = state.get("dbsnp_variants", {})
    has_variants = any(
        variant.get("found", False)
        for gene_variants in variants.values()
        for variant in gene_variants.values()
    )
    if DEBUG:
        _debug(f"Variants found: {has_variants}")
    return has_variants


def has_dbsnp_coordinates(state: "State") -> bool:
    variants = state.get("dbsnp_variants", {})
    has_coords = any(
        len(variant.get("coordinates", [])) > 0
        for gene_variants in variants.values()
        for variant in gene_variants.values()
    )
    if DEBUG:
        _debug(f"Coordinates found: {has_coords}")
    return has_coords
