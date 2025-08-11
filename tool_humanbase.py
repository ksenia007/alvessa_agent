"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-06-25
Updated: 2025-08-10


Description: 

Helpers to query HumanBase and MyGene.info. The latter is used to convert to entrez IDs, needed in HB"""

from __future__ import annotations
import requests
import time
from typing import Any, Dict, List, Optional
import pandas as pd
import numpy as np

DEBUG=True

def _symbol_to_entrez(symbol: str) -> Optional[str]:
    """Convert an HGNC symbol to an Entrez ID via MyGene.info"""
    if DEBUG:
        print(f"[HumanBase] Resolving symbol: {symbol}")
    try:
        r = requests.get(
            "https://mygene.info/v3/query",
            params={"q": symbol, "species": "human", "fields": "entrezgene", "size": 1},
            timeout=8,
        )
        r.raise_for_status()
        hits = r.json()["hits"]
        return None if not hits else str(hits[0]["entrezgene"])
    except Exception:
        return None


def _fetch_predictions_HB(entrez: str) -> List[Dict[str, float]]:
    """Download HumanBase functional predictions for a given Entrez ID."""
    if DEBUG:
        print(f"[HumanBase] Fetching predictions for Entrez ID: {entrez}")
    url = f"https://humanbase.io/api/genes/{entrez}/predictions/"
    r = requests.get(url, timeout=12)
    if r.status_code == 404:  # no gene found in HumanBase 
        return []
    r.raise_for_status()
    return r.json()  # list[{function, probability, tissue}]


def _filter_predictions_HB(
    preds: List[Dict[str, float]], *, threshold: float = 0.9
) -> List[Dict[str, float]]:
    """Keep only high-confidence predictions."""
    if DEBUG:
        print(f"[HumanBase] Filtering predictions with threshold: {threshold}")
    pared: List[Dict[str, str]] = []
    for hit in preds:
        if hit["score"] < threshold:
            continue
        t = hit["term"]
        pared.append(
            {
                "score": round(hit["score"], 3),
                "term": t["title"],
                "description": t["description"],
                "category": t["database"]["name"],
                "go_id": t["identifier"],
                "annotation_count": t.get("annotation_count", 0),
            }
        )
    return pared

def _fetch_tissue_variants_HB(entrez: str) -> Optional[pd.DataFrame]:
    """
    Download HumanBase tissue-specific variant effect predictions for a given Entrez ID.
    Note that this is only pre-computed scores, with hg19 build.
    
    Returns:
        pd.DataFrame: columns = chr, position, ref, alt, tissue, score
        None: if gene not found or no variants available.
    """
    if DEBUG:
        print(f"[HumanBase] Fetching tissue variant predictions for Entrez ID: {entrez}")

    url = f"https://humanbase.io/api/genes/{entrez}/tissue_variants"
    r = requests.get(url, timeout=12)

    if r.status_code == 404:  # no gene found in HumanBase
        if DEBUG:
            print(f"[HumanBase] No entry found for Entrez ID {entrez}")
        return None

    r.raise_for_status()
    data = r.json()

    # Validate presence of variants/tissues
    if not data.get("variants") or not data.get("tissues"):
        if DEBUG:
            print(f"[HumanBase] No variants or tissues found for Entrez ID {entrez}")
        return None

    records = []
    for variant_idx, var_meta in enumerate(data["variants"]):
        for tissue in data["tissues"]:
            tissue_name = tissue["title"]
            alt_scores = tissue["variants"][variant_idx]
            if alt_scores:  # dict of alt -> score
                for alt, score in alt_scores.items():
                    records.append({
                        "chr": var_meta["chr"],
                        "position": var_meta["position"],
                        "ref": var_meta["ref"],
                        "alt": alt,
                        "tissue": tissue_name,
                        "score": score
                    })

    if not records:
        if DEBUG:
            print(f"[HumanBase] No scored variants for Entrez ID {entrez}")
        return None

    return pd.DataFrame(records)

def summarize_tissue_variants_text_HB(df: Optional[pd.DataFrame]) -> Optional[str]:
    if df is None or df.empty:
        return None

    # Coverage
    n_variants = df[["chr", "position", "ref"]].drop_duplicates().shape[0]
    n_alts = df[["chr", "position", "ref", "alt"]].drop_duplicates().shape[0]
    n_tissues = df["tissue"].nunique()
    n_scores = len(df)

    # Score stats
    s = df["score"].describe(percentiles=[0.25, 0.5, 0.75])
    score_min, score_q25, score_median, score_q75, score_max, score_mean, score_std = [
        s[k] for k in ["min", "25%", "50%", "75%", "max", "mean", "std"]
    ]

    # Tissue means
    tm = df.groupby("tissue")["score"].mean().sort_values(ascending=False)
    top_tissues = ", ".join([f"{t} ({m:.4f})" for t, m in tm.head(3).items()])
    bottom_tissues = ", ".join([f"{t} ({m:.4f})" for t, m in tm.tail(3).items()])

    # Extremes
    max_row = df.loc[df["score"].idxmax()]
    min_row = df.loc[df["score"].idxmin()]
    pos_tissue = max_row["tissue"]
    neg_tissue = min_row["tissue"]

    def tissue_stats(tname):
        sub = df[df["tissue"] == tname]["score"]
        return f"mean={sub.mean():.4f}, median={sub.median():.4f}, min={sub.min():.4f}, max={sub.max():.4f}"

    def variant_other_context(row):
        sub = df[
            (df["chr"] == row["chr"]) &
            (df["position"] == row["position"]) &
            (df["ref"] == row["ref"]) &
            (df["alt"] == row["alt"])
        ]
        other = sub[sub["tissue"] != row["tissue"]]["score"]
        if other.empty:
            return "no other tissue data"
        mean_other = other.mean()
        sign_consistent = all((o >= 0) == (row["score"] >= 0) for o in other.dropna())
        if sign_consistent:
            sign_note = "all other tissues have the same direction of effect"
        else:
            sign_note = "effects vary in sign across tissues"
        return f"mean across other tissues={mean_other:.4f} ({sign_note})"

    pos_tissue_context = tissue_stats(pos_tissue)
    neg_tissue_context = tissue_stats(neg_tissue)
    pos_variant_context = variant_other_context(max_row)
    neg_variant_context = variant_other_context(min_row)

    return (
        "HumanBase tissue_variant summary "
        "(Expecto; positive=increased expression, negative=decreased; "
        "scores comparable across tissues/variants; subset of precomputed variants incl. common+selected).\n"
        f"Variants={n_variants}, Alts={n_alts}, Tissues={n_tissues}, Scores={n_scores}.\n"
        f"Score distribution: min={score_min:.4f}, Q25={score_q25:.4f}, median={score_median:.4f}, "
        f"Q75={score_q75:.4f}, max={score_max:.4f}, mean={score_mean:.4f}, std={score_std:.4f}.\n"
        f"Tissues with highest mean score: {top_tissues}.\n"
        f"Tissues with lowest mean score: {bottom_tissues}.\n"
        f"Strongest positive effect variant: {max_row['chr']}:{int(max_row['position'])} "
        f"{max_row['ref']}>{max_row['alt']} in {pos_tissue} ({max_row['score']:.4f}).\n"
        f"  Distribution of scores for all variants in {pos_tissue}: {pos_tissue_context}.\n"
        f"  In other tissues, this variant has {pos_variant_context}.\n"
        f"Strongest negative effect variant: {min_row['chr']}:{int(min_row['position'])} "
        f"{min_row['ref']}>{min_row['alt']} in {neg_tissue} ({min_row['score']:.4f}).\n"
        f"  Distribution of scores for all variants in {neg_tissue}: {neg_tissue_context}.\n"
        f"  In other tissues, this variant has {neg_variant_context}."
    )

from typing import Optional, Tuple, Dict

def describe_variant_in_tissues_HB(df: Optional[pd.DataFrame],
                                 chrom: str,
                                 position: int,
                                 ref: str,
                                 alt: str) -> Tuple[Optional[Dict], Optional[str]]:
    """
    Describe a specific variant's effects across tissues in the HumanBase dataset.
    Parameters:
        df: DataFrame with tissue variant predictions (columns: chr, position, ref, alt, tissue, score).
        chrom: Chromosome of the variant.
        position: Position of the variant.
        ref: Reference allele.
        alt: Alternate allele.
    Note: build should be hg19.
    """
    if df is None or df.empty:
        return None, None
    
    # enforce position to be integer in df and in parameters
    if not isinstance(position, int):
        position = int(position)
    df["position"] = df["position"].astype(int)
    
    # ensure ref and alt are strings and upper case
    ref = str(ref).upper()
    alt = str(alt).upper()

    # Subset for this variant
    sub = df[
        (df["chr"] == chrom) &
        (df["position"] == position) &
        (df["ref"] == ref) &
        (df["alt"] == alt)
    ]
    if sub.empty:
        return None, ""

    # Variant-wide stats
    score_desc = sub["score"].describe()
    mean_score = score_desc["mean"]
    median_score = score_desc["50%"]
    min_score = score_desc["min"]
    max_score = score_desc["max"]
    std_score = score_desc["std"]

    # Helper to get percentile rank of variant score within a tissue
    def score_percentile(tissue_name, score_value):
        scores_all = df[df["tissue"] == tissue_name]["score"]
        if scores_all.empty:
            return None
        # percentile of score_value among all variant scores in that tissue
        percentile = (np.sum(scores_all <= score_value) / len(scores_all)) * 100
        return round(percentile, 1)

    # Top 2 and bottom 2 tissues for this variant
    top_tissues_df = sub.sort_values("score", ascending=False).head(2)
    bottom_tissues_df = sub.sort_values("score", ascending=True).head(2)

    def tissue_records(df_slice):
        recs = []
        for _, row in df_slice.iterrows():
            pct = score_percentile(row["tissue"], row["score"])
            recs.append({
                "tissue": row["tissue"],
                "score": row["score"],
                "percentile_in_tissue": pct
            })
        return recs

    top_tissues = tissue_records(top_tissues_df)
    bottom_tissues = tissue_records(bottom_tissues_df)

    # Compare mean to dataset mean
    dataset_mean = df["score"].mean()
    relative_effect = "higher than average" if mean_score > dataset_mean else "lower than average"

    # Structured stats output
    stats_dict = {
        "variant": {"chr": chrom, "position": position, "ref": ref, "alt": alt},
        "n_tissues": sub["tissue"].nunique(),
        "mean_score": mean_score,
        "median_score": median_score,
        "min_score": min_score,
        "max_score": max_score,
        "std_score": std_score,
        "top_tissues": top_tissues,
        "bottom_tissues": bottom_tissues,
        "relative_effect": relative_effect
    }

    # Text description
    def tissue_str(tissues):
        parts = []
        for t in tissues:
            pct_str = f", percentile_in_tissue={t['percentile_in_tissue']}" if t['percentile_in_tissue'] is not None else ""
            parts.append(f"{t['tissue']} ({t['score']:.4f}{pct_str})")
        return ", ".join(parts)

    text_description = (
        f"Variant {chrom}:{position} {ref}>{alt} appears in {sub['tissue'].nunique()} tissues.\n"
        f"Scores: min={min_score:.4f}, median={median_score:.4f}, max={max_score:.4f}, "
        f"mean={mean_score:.4f}, std={std_score:.4f}.\n"
        f"Top tissues: {tissue_str(top_tissues)}.\n"
        f"Bottom tissues: {tissue_str(bottom_tissues)}.\n"
        f"Overall, this variant's mean effect is {relative_effect} compared to all variants in the dataset."
    )

    return stats_dict, text_description

# Agent-compatible node
from state import State  # noqa: E402


def humanbase_predictions_agent(state: "State") -> "State":
    """
    LangGraph node that annotates each gene with HumanBase predictions.

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state with the `"humanbase_predictions"` field filled.
    """
    preds = state.get("humanbase_predictions", {}).copy()

    for gene in state.get("genes", []):
        if gene in preds:
            continue

        entrez = _symbol_to_entrez(gene)
        if not entrez:
            preds[gene] = []
            continue

        try:
            tmp = _fetch_predictions_HB(entrez)
        except Exception as exc:
            print(f"[HumanBase] {gene}: {exc}")
            preds[gene] = []
        else:
            preds[gene] = _filter_predictions_HB(tmp, threshold=0.95)


    return {**state, "humanbase_predictions": preds}

def humanbase_expecto_agent(state: "State") -> "State":
    """
    LangGraph node that annotates each gene with Expecto predictions from HumanBase.

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state with the `"humanbase_predictions"` field filled.
    """
    if DEBUG:
        print("[HumanBase Expecto] Fetching tissue-specific predictions for genes")
    
    preds = state.get("humanbase_expecto", {}).copy()

    for gene in state.get("genes", []):
        if gene in preds:
            continue

        entrez = _symbol_to_entrez(gene)
        if not entrez:
            preds[gene] = {}
            continue

        try:
            tmp = _fetch_tissue_variants_HB(entrez)
        except Exception as exc:
            print(f"[HumanBase] {gene}: {exc}")
            preds[gene] = {}
        else:
            summary_text = summarize_tissue_variants_text_HB(tmp)
            if summary_text is None:
                preds[gene] = {}
            else:
                preds[gene] = {
                    "summary_text": summary_text,
                    "dataframe": tmp
                }

    return {**state, "humanbase_expecto": preds}

def humanbase_tissue_expecto_annotate_variants(state: "State") -> "State":
    """
    LangGraph node that annotates each variant with HumanBase tissue-specific predictions.

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state with the `"humanbase_tissue_expecto"` field filled.
    """
    print('***'*20)
    print(state)
    print('***'*20)
    if DEBUG:
        print("[HumanBase Expecto Variant] Annotating variants with tissue-specific predictions")
    variant_descr = {}
    preds = state.get("humanbase_expecto", {}).copy()
    variants = state.get("dbsnp_variants", {}).copy()
    if not variants:
        return {**state, "tissue_expression_preds_variant": []}

    # annotate each variant with tissue-specific predictions
    for gene, gene_vars in variants.items():
        if DEBUG:
            print('[HumanBase Expecto Variant] Processing gene:', gene)
        if gene not in preds:
            continue
        gene_descr = f"Gene: {gene}: " 
        if DEBUG:
            print(preds[gene])
        for variant_id, var_data in gene_vars.items():
            # use describe_variant_in_tissues_HB
            all_snps = var_data['coordinates']
            for snp in all_snps:
                if 'GRCh38' in snp['assembly']:
                    continue
                chrom = snp['chrom']
                if 'chr' not in chrom:
                    chrom = 'chr' + str(chrom)
                pos = snp['pos']
                ref = snp['ref']
                alt = snp['alt']
                if DEBUG:
                    print(f"[HumanBase Expecto per variant] Processing variant {chrom}:{pos} {ref}>{alt} for gene {gene}")
                descr = describe_variant_in_tissues_HB(
                    preds[gene]["dataframe"],
                    chrom,
                    pos,
                    ref,
                    alt
                )
                if descr[0] is None:
                    continue
                
                if len(descr[1])>1:
                    gene_descr += f" Variant {chrom}:{pos} {ref}>{alt}: {descr[1]};"
            
        if gene_descr:
            variant_descr[gene] = gene_descr
        
    return {
        **state, 
        "tissue_expression_preds_variant_text_description": variant_descr
    }