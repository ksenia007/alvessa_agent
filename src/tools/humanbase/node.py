"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-06-25
Updated: 2025-08-10


Description: 

Helpers to query HumanBase and MyGene.info. The latter is used to convert to entrez IDs, needed in HB"""

from __future__ import annotations
import json
import requests
import time
from typing import Any, Dict, List, Optional, Tuple
import pandas as pd
import numpy as np
# Agent-compatible node
from src.state import State 
from src.tools.base import Node

DEBUG=True



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

def _fetch_expectosc_variants_HB(entrez: str) -> Optional[pd.DataFrame]:
    """
    Download HumanBase tissue-specific variant effect predictions for a given Entrez ID.
    Note that this is only pre-computed scores, with hg19 build.
    
    Returns:
        pd.DataFrame: columns = chr, position, ref, alt, tissue, score
        None: if gene not found or no variants available.
    """
    if DEBUG:
        print(f"[HumanBase] Fetching tissue variant predictions for Entrez ID: {entrez}")

    url = f"https://humanbase.io/api/genes/{entrez}/tissue_variants/?database=clever-tissues&collapse=true"
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
                        # "ref": "", var_meta["ref"],
                        "alt": alt,
                        "tissue": tissue_name,
                        "score": score
                    })

    if not records:
        if DEBUG:
            print(f"[HumanBase] No scored variants for Entrez ID {entrez}")
        return None

    return pd.DataFrame(records)


def describe_variant_in_celltypes_HB(df: Optional[pd.DataFrame],
                                 chrom: str,
                                 position: int,
                                 alt: str, 
                                 variant_id: str) -> Tuple[Optional[Dict], Optional[str]]:
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
    num_label = {}
    if df is None or df.empty:
        return None, None
    
    # enforce position to be integer in df and in parameters
    if not isinstance(position, int):
        position = int(position)
    df["position"] = df["position"].astype(int)
    
    # ensure ref and alt are strings and upper case
    # ref = str(ref).upper()
    alt = str(alt).upper()
    
    if 'chr' not in chrom:  
        chrom = 'chr' + str(chrom)
    if DEBUG:
        print(f"[HumanBase Expecto Variant] Describing variant {chrom}:{position} >{alt}")

    # Subset for this variant
    sub = df[
        (df["chr"] == chrom) &
        (df["position"] == position) &
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
        "variant": {"chr": chrom, "position": position, "alt": alt},
        "variant_id": variant_id,
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
            # pct_str = f", percentile_in_tissue={t['percentile_in_tissue']}" if t['percentile_in_tissue'] is not None else ""
            # parts.append(f"{t['tissue']} ({t['score']:.4f}{pct_str})")
            if np.abs(t['score']) >= 1 :
                parts.append(f"{t['tissue']} ({t['score']:.4f})")
        return ", ".join(parts)

    top_tissues = f"Top: {tissue_str(top_tissues)}." if np.abs(max_score) >= 1 else ""
    bottom_tissues = f"Bottom: {tissue_str(bottom_tissues)}." if np.abs(min_score) >= 1 else ""
    

    text_description = (
        f"Scores across cell types: min={min_score:.4f}, median={median_score:.4f}, max={max_score:.4f}, "
        f"{top_tissues}"
        f"{bottom_tissues}"
        # f"Overall, this variant's mean effect is {relative_effect} compared to other variants for this gene."
    )

    return stats_dict, text_description, sub[["tissue", "score"]].to_dict("records")



def expectosc_predictions_agent(state: "State") -> "State":
    """
    Annotate variants with HumanBase ExpectoSC (cell-type/tissue) scores (GrCh37).
      1) Fetch per gene in gene_entities, annotate any variants that appear.
      2) For variants whose related genes weren't fetched, fetch-on-demand and annotate.
    Writes to each Variant:
      - add_functional_prediction(gene, "ExpectoSC", text)
      - update_text_summaries(text)
    """
    if DEBUG:
        print(f"[ExpectoSC] started...")

    variants = (state.get("variant_entities") or {}).copy()
    genes    = (state.get("gene_entities")    or {}).copy()
    if not variants or not genes:
        if DEBUG: print("[ExpectoSC] nothing to do (no variants or genes).")
        return

    def _get_entrez(gsym, gobj):
        eid = getattr(gobj, "entrez_id", None)
        if not eid and hasattr(gobj, "get_entrez_id"):
            try: eid = gobj.get_entrez_id()
            except Exception: eid = None
        if not eid:
            eid = None
        return str(eid) if eid else None

    def _annotate_if_match(var_id, var_obj, gene_sym, df):
        # GrCh37 coordinates from Variant
        locs = var_obj.get_location("GrCh37")
        if not locs: return 0
        chrom, pos = locs.get("chrom"), locs.get("pos")
        refs, alts = (locs.get("ref") or []), (locs.get("alt") or [])
        if chrom is None or pos is None or not refs or not alts: return 0
        if len(refs) != 1: return 0
        ref = refs[0]
        chr_str = f"chr{chrom}"

        hits = 0
        for alt in alts:
            try:
                dict_summary, text, variant_rows = describe_variant_in_celltypes_HB(df, chr_str, pos, alt, var_id)
            except Exception:
                print(f"[ExpectoSC] describe failed for var {var_id} {chr_str}:{pos}>{alt}")
                continue
            if not text:
                print(f"[ExpectoSC] no data for var {var_id} {chr_str}:{pos}>{alt}")
                continue
            var_obj.update_text_summaries(f"*ExpectoSC: Cell-type specific predicted GE disruption (target = {gene_sym}): {text}.")
            if not dict_summary:
                continue

            compact_top = [
                {"tissue": rec.get("tissue"), "score": float(rec.get("score", 0))}
                for rec in (dict_summary.get("top_tissues") or [])
                if rec and rec.get("tissue") is not None and rec.get("score") is not None
            ]
            compact_bottom = [
                {"tissue": rec.get("tissue"), "score": float(rec.get("score", 0))}
                for rec in (dict_summary.get("bottom_tissues") or [])
                if rec and rec.get("tissue") is not None and rec.get("score") is not None
            ]
            numeric_summary = {
                "top_tissues": compact_top,
                "bottom_tissues": compact_bottom,
                "scores": [
                    {
                        "tissue": rec.get("tissue"),
                        "score": float(rec.get("score", 0)),
                    }
                    for rec in (variant_rows or [])
                    if rec and rec.get("tissue") is not None and rec.get("score") is not None
                ],
            }
            try:
                payload = json.dumps(numeric_summary, ensure_ascii=False)
            except Exception:
                payload = json.dumps({"top_tissues": [], "bottom_tissues": [], "scores": []})
            var_obj.add_functional_prediction(gene_sym, "ExpectoSC", payload)
            hits += 1
        return hits

    fetched: dict[str, Optional[pd.DataFrame]] = {}
    fetched_genes: set[str] = set()
    total_ann = 0

    # --- Pass 1: fetch for genes in gene_entities, annotate matching variants ---
    for gsym, gobj in genes.items():
        entrez = _get_entrez(gsym, gobj)
        if not entrez:
            continue
        try:
            df = _fetch_expectosc_variants_HB(entrez)  # returns DataFrame or None
        except Exception as e:
            if DEBUG: print(f"[ExpectoSC] fetch failed {gsym}/{entrez}: {e}")
            df = None
        fetched[gsym] = df
        fetched_genes.add(gsym)
        if df is None or df.empty:
            continue
        
        gobj.update_text_summaries("*ExpectoSC: Cell-type specific GE disruption scores (scaled z vs 1000 Genomes; >3 high confidence; abs(score)<1 hides top tissues; predicts variants +/-20kbp of TSS; positive=increased expression, negative=decreased).")
        
        # annotate any variant that lists this gene
        for var_id, var_obj in variants.items():
            rel = (var_obj.get_related_genes() or [])
            if gsym in rel:
                total_ann += _annotate_if_match(var_id, var_obj, gsym, df)

    # --- Pass 2: for variants referencing genes we didn't fetch, fetch-on-demand & annotate ---
    for var_id, var_obj in variants.items():
        rel = (var_obj.get_related_genes() or [])
        for gsym in rel:
            if gsym in fetched_genes:
                continue
            gobj = genes.get(gsym)
            if not gobj:
                continue
            entrez = _get_entrez(gsym, gobj)
            if not entrez:
                continue
            try:
                df = _fetch_expectosc_variants_HB(entrez)
            except Exception as e:
                if DEBUG: print(f"[ExpectoSC] fetch (late) failed {gsym}/{entrez}: {e}")
                df = None
            fetched[gsym] = df
            fetched_genes.add(gsym)
            if df is None or df.empty:
                continue
            total_ann += _annotate_if_match(var_id, var_obj, gsym, df)

    if DEBUG:
        print(f"[ExpectoSC] done — annotated {total_ann}")

    # nothing to put back into state — annotations live on Variant objects
    return

    
    
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
    gene_objs = state.get("gene_entities", {})
    
    print(gene_objs)

    for gene in gene_objs.keys():
        if gene_objs[gene].has_tool("humanbase_predictions"):
            continue
        
        entrez = gene_objs[gene].entrez_id
        if not entrez:
            if DEBUG:
                print(f"[HumanBase] Could not find Entrez ID for gene symbol: {gene}")
            gene_objs[gene].add_tool("humanbase_predictions")
            continue
        try:
            tmp = _fetch_predictions_HB(entrez)
        except Exception as exc:
            gene_objs[gene].add_tool("humanbase_predictions")
            continue
        
        low_thr = _filter_predictions_HB(tmp, threshold=0.7)
        high_thr = _filter_predictions_HB(tmp, threshold=0.95)
        # split by category - disease ontology and GO terms
        disease_terms = []
        for d in low_thr:
            if 'category' not in d:
                continue
            if d['category'].lower() in ['disease ontology', 'disease']:
                disease_terms.append(d)
        
        go_terms = []
        for d in high_thr:
            if 'category' not in d:
                continue
            if d['category'].lower() in ['gene ontology (bp)', 'go']:
                go_terms.append(d)
        
        # convert to list of strings - term + description for disease and term for GO
        go_terms_list = [list(set([d['term'] for d in go_terms]))][0]
        disease_list = [f"{d['term']}: {d['description']} (score={d['score']})" for d in disease_terms]
        
        # collect predicted Disease Ontology (finter by category) and separately GO terms
        gene_objs[gene].add_many_predicted_diseases(disease_list)
        gene_objs[gene].add_many_predicted_go(go_terms_list)
        gene_objs[gene].add_tool("humanbase_predictions")
        
        # add text description of the disease terms
        if disease_list:
            disease_text = f"*HumanBase: Predicted disease associations for {gene} (score>0.7): " + "; ".join(disease_list) + "."
            gene_objs[gene].update_text_summaries(disease_text)

    return 


NODES: tuple[Node, ...] = (
    Node(
        name="humanbase_functions",
        entry_point=humanbase_predictions_agent,
        description="Fetch per-gene functional predictions from HumanBase tissue-specific networks. Provides expanded list of functions.",
    ),
    Node(
        name="expectosc_predictions_agent",
        entry_point=expectosc_predictions_agent,
        description="Annotates variants with predicted (from sequence) cell type-specific expression disruption predictions. Requires variant_annotations to be run first",
        dependencies=("variant_annotations",),
    ),
)
