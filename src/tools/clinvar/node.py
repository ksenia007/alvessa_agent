"""
Description:

ClinVar tools (gene-level diseases and variant summaries)

"""
from __future__ import annotations

import time

from src.state import State
from src.tools.base import Node
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import json
import pandas as pd

DEBUG = True

CLINVAR_GENE_DISEASES = Path(__file__).resolve().parents[3] / "local_dbs" / "clinvar" / "clinvar_gene_condition_source_id.parquet"
CLINVAR_DF = Path(__file__).resolve().parents[3] / "local_dbs" / "clinvar" / "clinvar_pathogenic_variants_grch38.parquet"

def _load_gene_disease_df():
    return pd.read_parquet(CLINVAR_GENE_DISEASES)


def _load_variant_df():
    return pd.read_parquet(CLINVAR_DF)


def clinvar_gene_node(state: "State") -> "State":
    """Annotate genes with ClinVar gene–disease associations (no variant objects)."""
    gene_entities = state.get("gene_entities") or {}
    try:
        df = _load_gene_disease_df()
    except Exception as e:
        if DEBUG:
            print(f"[ClinVar-Gene] Error loading ClinVar gene-disease data: {e}")
        return state

    genes_with_disease = set(df.AssociatedGenes.dropna().unique())

    for gene_name, gene in gene_entities.items():
        if not gene_name:
            continue
        try:
            if gene_name in genes_with_disease:
                subset = df[df["AssociatedGenes"] == gene_name]
                diseases = subset[["DiseaseName", "SourceName", "LastUpdated"]].drop_duplicates()
                lines = []
                for _, row in diseases.iterrows():
                    lines.append(
                        f"{row['DiseaseName']} (source: {row['SourceName']}, updated: {row['LastUpdated']})"
                    )
                if lines:
                    gene.update_text_summaries(f"*ClinVar: Diseases for {gene_name}: " + "; ".join(lines) + ".")
                    gene.add_tool("clinvar_gene_node")
        except Exception as exc:
            if DEBUG:
                print(f"[ClinVar-Gene] Error processing {gene_name}: {exc}")

    return state


def clinvar_variants_node(state: "State") -> "State":
    """Summarize ClinVar pathogenic variants per gene as text (no Variant objects created)."""
    gene_entities = state.get("gene_entities") or {}
    try:
        df = _load_variant_df()
    except Exception as e:
        if DEBUG:
            print(f"[ClinVar-Variants] Error loading ClinVar variant data: {e}")
        return state

    genes_with_variants = set(df.GeneSymbol.dropna().unique())

    for gene_name, gene in gene_entities.items():
        if not gene_name:
            continue
        try:
            if gene_name in genes_with_variants:
                subset = df[df["GeneSymbol"] == gene_name]
                variant_lines = []
                for _, row in subset.iterrows():
                    parts = []
                    rsid = str(row.get("rsid", "")).strip()
                    if rsid:
                        rsid = rsid if rsid.startswith("rs") else f"rs{rsid}"
                        parts.append(f"{rsid}:")
                    for col in ["Type", "Name", "ClinicalSignificance"]:
                        val = row.get(col)
                        if pd.isna(val) or val in ["-", "na", ""]:
                            continue
                        parts.append(str(val))
                    if parts:
                        variant_lines.append(" ".join(parts))
                if variant_lines:
                    gene.update_text_summaries(
                        f"*ClinVar variants: Pathogenic/likely pathogenic records for {gene_name}: "
                        + "; ".join(variant_lines)
                        + "."
                    )
                    gene.add_tool("clinvar_variants_node")
        except Exception as exc:
            if DEBUG:
                print(f"[ClinVar-Variants] Error processing {gene_name}: {exc}")

    return state

NODES: tuple[Node, ...] = (
    Node(
        name="clinvar_gene_node",
        entry_point=clinvar_gene_node,
        description=(
            "Annotates genes with ClinVar gene–disease associations (disease name, source, last updated). "
        ),
    ),
    Node(
        name="clinvar_variants_node",
        entry_point=clinvar_variants_node,
        description=(
            "Retrieves pathogenic and likely pathogenic variants for input genes (rsID, variant type, name, clinical significance) from the ClinVar database. "
            "Should only be called if exact variant information such as amino acid substitution or variant ID is required to answer the question, otherwise use ClinVar gene-level."
            " Might create excessively long context for commonly implicated genes."
        ),
    ),
)
