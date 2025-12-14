"""
Author: Keerthana Nallamotu <kn6412@princeton.edu>
Contributors: 
Created: 2024-07-11
Updated: 2025-07-14


Description: 
Tool to summarize the GO terms associated with a list of genes"""

from __future__ import annotations

import time
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

from src.state import State
from src.tools.base import Node
from src.tools.shared.word2vec import fps_word2vec
from src.tools.uniprot.node import get_uniprot_entry_for_gene, extract_GO_from_uniprot_entry
from src.tools.shared.go_enrichment import run_go_enrichment
from src.tools.biogrid.utils import _fetch_predictions_BioGRID

DEBUG=True

REPO_ROOT = Path(__file__).resolve().parents[3]
LOCAL_DBS = REPO_ROOT / "local_dbs"
OLD_GO_DIR = LOCAL_DBS / "old_go_data"
PAN_GO_DIR = LOCAL_DBS / "pan_go_data"

GO_CATEGORIES: Sequence[str] = ("BP", "MF", "CC")
CATEGORY_LABELS: Dict[str, str] = {
    "BP": "Biological Process (larger sequences of molecular activities)",
    "MF": "Molecular Function (biochemical activities of gene products)",
    "CC": "Cellular Component (locations within the cell where gene products are active)",
}
GO_SOURCES = {
    "old_go": {
        "dir": OLD_GO_DIR,
        "pattern": "goa_human_{category}.gmt",
        "obo": "go.obo",
    },
    "pan_go": {
        "dir": PAN_GO_DIR,
        "pattern": "functionome_release_{category}.gmt",
        "obo": "go.obo",
    },
}

def go_summarization_agent(state: "State", embedding_method: str = 'word2vec') -> "State":
    """
    Summarize GO terms per gene and append a short line to each Gene's text summaries.
    Prefers the Gene.go_annotations field; falls back to UniProt if empty.
    """
    gene_objs = (state.get("gene_entities") or {}).copy()

    for gene in gene_objs.keys():
        gobj = gene_objs.get(gene)
        if not gobj:
            continue

        # 1) Prefer existing GO annotations on the Gene objects
        terms = getattr(gobj, "go_annotations", []) or []
    
        # 2) If none, fetch from UniProt
        if not terms:
            try:
                entry = get_uniprot_entry_for_gene(gene)
                if entry:
                    terms = extract_GO_from_uniprot_entry(entry) or []
            except Exception:
                terms = []

        if not terms:
            continue

        # Dedup and select representative subset
        uniq = list(dict.fromkeys(terms))
        if embedding_method.lower() == "word2vec":
            try:
                idxs = fps_word2vec(uniq, 6, separate_sampling=True)
                picked = [uniq[i] for i in idxs]
            except Exception:
                picked = uniq[:6]
        else:  # fallback (tf-idf or anything else)
            raise NotImplementedError(f"Unknown embedding method: {embedding_method}")

        gobj.update_text_summaries(f"*GO: {gene} representative GO terms: " + ", ".join(picked) + ".")

    time.sleep(0.1)
    return


def _load_enrichment(
    partners: Iterable[str],
    *,
    source: str,
    category: str,
) -> List[str]:
    config = GO_SOURCES[source]
    gmt_path = config["dir"] / config["pattern"].format(category=category)
    obo_path = config["dir"] / config["obo"]

    if not gmt_path.exists() or not obo_path.exists():
        warnings.warn(
            f"[BioGRID] Missing GO enrichment resource for {source}/{category}: {gmt_path}"
        )
        return []

    df = run_go_enrichment(list(partners), str(gmt_path), str(obo_path))
    if df is None or df.empty:
        return []

    return [term for term in df["GO_Name"].tolist() if isinstance(term, str) and term]


def summarize_biogrid_go_agent(state: "State") -> "State":
    """Run GO enrichment on BioGRID interactors for each gene."""

    gene_entities = state.get("gene_entities") or {}
    if not gene_entities:
        return

    for gene_name, gene in gene_entities.items():
        if not gene_name:
            continue

        try:
            _, human_set, _ = _fetch_predictions_BioGRID(gene_name)
        except Exception as exc:
            print(f"[BioGRID] {gene_name}: {exc}")
            continue

        partners = {
            (partner or "").strip()
            for partner_list in human_set.values()
            for partner in partner_list
            if partner
        }

        if not partners:
            continue

        summary_lines: List[str] = []

        for category in GO_CATEGORIES:
            for source in GO_SOURCES.keys():
                enriched_terms = _load_enrichment(partners, source=source, category=category)
                if not enriched_terms:
                    continue

                gene.add_interactors_enriched_go_terms(source, category, enriched_terms)

                if source == "pan_go":
                    summary_lines.append(
                        (
                            f"*BioGRID GO: Enriched PAN-GO {CATEGORY_LABELS.get(category, category)} terms "
                            f"for interactors of {gene.symbol}: {', '.join(enriched_terms)}."
                        )
                    )

        if summary_lines:
            gene.update_text_summaries(" ".join(summary_lines))
            gene.add_tool("Summarize_bioGRID_GO")

    return


NODES: tuple[Node, ...] = (
    Node(
        name="Summarize_GO",
        entry_point=go_summarization_agent,
        description="Summarizes GO terms for the input genes, to condense long list into representative terms.",
    ),
    Node(
        name="Summarize_bioGRID_GO",
        entry_point=summarize_biogrid_go_agent,
        description=(
            "This is a summarized version of BIOGRID that should be run instead of BioGRID when individual interactor names are not needed to answer the query. Provides a summarized list of GO terms significantly enriched for the interacting genes of each input gene fetched from BioGRID."
        ),
        dependencies=("BioGRID",),
    ),
)
