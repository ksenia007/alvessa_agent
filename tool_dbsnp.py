from __future__ import annotations
import time
import re
from typing import Set, Dict, Any, List

from config import DEBUG, DBSNP_DEFAULT_ASSEMBLY
from tools.dbsnp.query import get_variant_info
from state import State
from tools.gencode.query import annotate_variant
import statistics


def _create_coordinate_summary(rsid: str, coordinates: List[Dict[str, Any]]) -> str:
    """Generates a concise, human-readable summary of variant coordinates."""
    if not coordinates:
        return f"{rsid} has no coordinate data available"

    primary_coord = next((c for c in coordinates if 'GRCh38' in c.get("assembly", "")), coordinates[0])
    chrom = primary_coord.get("chrom", "Unknown")
    pos38 = primary_coord.get("pos", "Unknown")
    
    return f"{rsid} is located at chr:{chrom}, pos: {pos38} (GRCh38)"
    
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
        p4 = (p*n + 2) / n4
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
        f"{label}: weighted AF {f(pooled)} (over {len(items)} studies), "
        f"median AF {f(med_all)}"
        + (f", median for studies with n>={min_people:,} ppl = {f(med_hi)}" if med_hi is not None else ""),
        f"agreement I^2 heterogeneity statistic={I2:.1f}%",
        "Extreme AFs: low in " + ", ".join(f"{x['study']} {f(x['p'])} (n={x['n']:,})" for x in low2),
        "high in " + ", ".join(f"{x['study']} {f(x['p'])} (n={x['n']:,})" for x in reversed(high2)),
        "largest-N studies: " + ", ".join(f"{x['study']} {f(x['p'])} (n={x['n']:,})" for x in reversed(topN2)),
    ]
    return "; ".join(parts)+". "
    

def dbsnp_variants_agent(state: "State", assembly: str = None, include_population_summaries: bool = False) -> "State":
    """
    LangGraph node that annotates SNPs with dbSNP variant information.
    Extracts rsIDs from GWAS associations and queries dbSNP for each variant.
    """
    assembly = assembly or DBSNP_DEFAULT_ASSEMBLY
    if DEBUG:
        print(f"[dbSNP] Using genome assembly: {assembly}")

    variants = state.get("variant_entities", {}).copy()
    
    rsids = variants.keys()
    
    if DEBUG:
        print(f"[dbSNP] Found {len(rsids)} total rsIDs to query across genes.")
    
    for rsid in rsids:
        result = get_variant_info(rsid, assembly)
        if not result: continue
        coordinates = result.get("coordinates", [])
        if not coordinates:
            if DEBUG:
                print(f"[dbSNP] No coordinates found for {rsid}. Skipping annotation.")
            continue

        for coord in coordinates:
            variants[rsid].add_location(
                build=coord.get("assembly", "Unknown"),
                chrom=coord.get("chrom", "Unknown"),
                pos=coord.get("pos", 0),
                ref=coord.get("ref", ""),
                alt=coord.get("alt", "")
            )
        variants[rsid].add_tool("dbsnp_variants_agent")
        variants[rsid].add_af_freq(result.get("allele_frequencies", []))
        
        variant_summary = f" {rsid} was found in dbSNP. "
        variant_summary += _create_af_summary_sorted_by_counts(result.get("allele_frequencies", []))
        
        variants[rsid].update_text_summaries(variant_summary)
    return



def _extract_rsids_from_text(text: str) -> Set[str]:
    """Extracts rsIDs from text using regex pattern matching."""
    pattern = r'rs\d+'
    matches = re.findall(pattern, text, re.IGNORECASE)
    return set(match.lower() for match in matches)


def has_dbsnp_variants(state: "State") -> bool:
    """Edge-condition helper: returns `True` if any dbSNP variants were found."""
    variants = state.get("dbsnp_variants", {})
    has_variants = any(
        variant.get("found", False)
        for gene_variants in variants.values()
        for variant in gene_variants.values()
    )
    if DEBUG:
        print(f"[has_dbsnp_variants] Variants found: {has_variants}")
    return has_variants


def has_dbsnp_coordinates(state: "State") -> bool:
    """Edge-condition helper: returns `True` if any variants have coordinate information."""
    variants = state.get("dbsnp_variants", {})
    has_coords = any(
        len(variant.get("coordinates", [])) > 0
        for gene_variants in variants.values()
        for variant in gene_variants.values()
    )
    if DEBUG:
        print(f"[has_dbsnp_coordinates] Coordinates found: {has_coords}")
    return has_coords
