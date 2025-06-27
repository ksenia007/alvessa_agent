"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-06-25
Updated: 2025-06-26


Description: 

Helpers to query HumanBase and MyGene.info. The latter is used to convert to entrez IDs, needed in HB"""

from __future__ import annotations
import requests
import time
from typing import Any, Dict, List, Optional

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

        time.sleep(0.3)  # courteous pause

    return {**state, "humanbase_predictions": preds}
