"""Shared helpers for Gene Ontology enrichment calculations."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping, Sequence

import pandas as pd
from scipy.stats import fisher_exact, hypergeom

try:  # pragma: no cover - optional dependency
    import statsmodels.stats.multitest as smm
except ImportError:  # pragma: no cover - handled at runtime
    smm = None


def parse_obo_names(obo_file: str | Path) -> Mapping[str, str]:
    """Return a mapping from GO term ID to human-readable name."""
    go_names: dict[str, str] = {}
    current_id, current_name = None, None
    with open(obo_file, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("id: GO:"):
                current_id = line.split("id: ")[1]
            elif line.startswith("name:") and current_id:
                current_name = line.split("name: ")[1]
                go_names[current_id] = current_name
                current_id, current_name = None, None
    return go_names


def run_go_enrichment(
    gene_set: Iterable[str],
    gmt_file: str | Path,
    obo_file: str | Path,
    *,
    method: str = "fisher",
    fdr_threshold: float = 0.05,
    min_overlap: int = 2,
):
    """Perform GO enrichment using either Fisher's exact test or hypergeometric test."""

    go_terms: dict[str, set[str]] = {}
    all_genes: set[str] = set()
    with open(gmt_file, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            term = parts[0]
            genes = set(parts[2:])
            go_terms[term] = genes
            all_genes.update(genes)

    go_names = parse_obo_names(obo_file) if obo_file else {}

    gene_list = set(gene_set) & all_genes
    if not gene_list:
        return pd.DataFrame()

    population_size = len(all_genes)
    input_count = len(gene_list)

    results = []
    for term, term_genes in go_terms.items():
        term_gene_count = len(term_genes)
        overlap = len(gene_list & term_genes)

        if method == "fisher":
            if overlap == 0:
                continue

            a = overlap
            b = input_count - overlap
            c = term_gene_count - overlap
            d = population_size - term_gene_count - b
            _, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
        elif method == "hypergeom":
            if overlap < min_overlap:
                continue
            p_value = hypergeom.pmf(overlap, population_size, term_gene_count, input_count)
        else:
            raise ValueError("Method must be 'fisher' or 'hypergeom'")

        results.append(
            {
                "GO_Term": term,
                "GO_Name": go_names.get(term, "NA"),
                "Overlap_Count": overlap,
                "GO_Term_Gene_Count": term_gene_count,
                "Input_Gene_Count": input_count,
                "P-value": p_value,
                "Genes": ", ".join(sorted(gene_list & term_genes)),
            }
        )

    df = pd.DataFrame(results)
    if df.empty:
        return df

    if smm is None:
        raise ImportError("statsmodels is required to compute FDR for GO enrichment results.")

    df["FDR"] = smm.multipletests(df["P-value"], method="fdr_bh")[1]
    return df[df["FDR"] < fdr_threshold].sort_values("FDR")
