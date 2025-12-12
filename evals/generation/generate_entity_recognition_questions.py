"""
Generate sample entity-recognition questions for benchmarks.

Creates five CSVs under benchmarks_generation/questions/components/entity_recognition:
  - simple_genes.csv     : single-gene queries
  - multi_genes.csv      : 5–6 gene queries
  - variants.csv         : 5–6 rsID queries (ClinVar-derived)
  - drugs.csv            : 5–6 drug queries (MedChemExpress library)
  - mirnas.csv           : miRNA queries (miRDB-derived)

All outputs use a fixed random seed for reproducibility and rely on local
resources already shipped with the repo.
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Iterable, List, Sequence
import json

import pandas as pd

# -------------------------
# Config
# -------------------------
SEED = 42
OUTPUT_DIR = Path("benchmarks_generation/questions/components/entity_recognition")

# Local data sources
GENCODE_PARQUET = Path("local_dbs/gencode.v48.expanded.parquet")
CLINVAR_PARQUET = Path("local_dbs/clinvar/clinvar_pathogenic_variants_grch38.parquet")
MEDCHEMEXPRESS_CSV = Path("local_dbs/compound_library_medchemexpress.csv")
MIRDB_PREDICTIONS = Path("local_dbs/miRDB_v6.0_prediction_result.txt")


def _seed_everything(seed: int) -> None:
    random.seed(seed)
    try:
        import numpy as np

        np.random.seed(seed)
    except Exception:
        pass


def _as_list(values: Iterable[str], k: int) -> List[str]:
    """Sample k unique items from the pool."""
    pool = list(values)
    if len(pool) < k:
        raise ValueError(f"Pool too small: need {k}, have {len(pool)}")
    return random.sample(pool, k)


def _load_gene_pool() -> List[str]:
    df = pd.read_parquet(GENCODE_PARQUET, columns=["gene_name"])
    genes = sorted(set(df["gene_name"].dropna().astype(str)))
    return [g for g in genes if g and g.lower() != "nan"]


def _load_variant_pool() -> List[str]:
    df = pd.read_parquet(CLINVAR_PARQUET, columns=["rsid"])
    rsids = []
    for val in df["rsid"].dropna().astype(str):
        rsid = val if val.startswith("rs") else f"rs{val}"
        rsids.append(rsid)
    return sorted(set(rsids))


def _load_drug_pool() -> List[str]:
    df = pd.read_csv(MEDCHEMEXPRESS_CSV)
    if "Product Name" in df.columns:
        names = df["Product Name"]
    elif "ProductName" in df.columns:
        names = df["ProductName"]
    else:
        raise ValueError("Expected 'Product Name' or 'ProductName' column in MedChemExpress CSV.")
    drugs = sorted({str(n).strip() for n in names.dropna() if str(n).strip()})
    return drugs


def _load_mirna_pool() -> List[str]:
    df = pd.read_csv(
        MIRDB_PREDICTIONS,
        sep="\t",
        header=None,
        names=["miRNAID", "geneID", "confidence"],
        usecols=["miRNAID"],
    )
    mirnas = sorted(set(df["miRNAID"].dropna().astype(str)))
    return mirnas


def _save_csv(name: str, rows: Sequence[dict]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_DIR / name, index=False)


def _join_entities(entities: Sequence[str]) -> str:
    return ", ".join(entities)


def build_simple_gene_questions(genes: List[str], n: int = 30) -> List[dict]:
    rows = []
    for gene in _as_list(genes, n):
        query = f"What is known about {gene} in oncology?"
        rows.append({"query": query, "recognized_entities": json.dumps([gene])})
    return rows


def build_multi_gene_questions(genes: List[str], n: int = 20, group_size: int = 5) -> List[dict]:
    rows = []
    pool = _as_list(genes, n * group_size)
    for i in range(n):
        chunk = pool[i * group_size : (i + 1) * group_size]
        query = f"Compare the roles of {', '.join(chunk[:-1])}, and {chunk[-1]} in immune signaling."
        rows.append({"query": query, "recognized_entities": json.dumps(chunk)})
    return rows


def build_variant_questions(rsids: List[str], n: int = 20, group_size: int = 5) -> List[dict]:
    rows = []
    pool = _as_list(rsids, n * group_size)
    for i in range(n):
        chunk = pool[i * group_size : (i + 1) * group_size]
        query = f"Interpret pathogenicity for {', '.join(chunk)}."
        rows.append({"query": query, "recognized_entities": json.dumps(chunk)})
    return rows


def build_drug_questions(drugs: List[str], n: int = 20, group_size: int = 5) -> List[dict]:
    rows = []
    pool = _as_list(drugs, n * group_size)
    for i in range(n):
        chunk = pool[i * group_size : (i + 1) * group_size]
        query = f"Which targets are modulated by {', '.join(chunk[:-1])}, and {chunk[-1]}?"
        rows.append({"query": query, "recognized_entities": json.dumps(chunk)})
    return rows


def build_mirna_questions(mirnas: List[str], n: int = 25) -> List[dict]:
    rows = []
    for mirna in _as_list(mirnas, n):
        query = f"List predicted targets for {mirna}."
        rows.append({"query": query, "recognized_entities": json.dumps([mirna])})
    return rows


def main() -> None:
    _seed_everything(SEED)

    genes = _load_gene_pool()
    rsids = _load_variant_pool()
    drugs = _load_drug_pool()
    mirnas = _load_mirna_pool()

    simple_gene_rows = build_simple_gene_questions(genes)
    multi_gene_rows = build_multi_gene_questions(genes)
    variant_rows = build_variant_questions(rsids)
    drug_rows = build_drug_questions(drugs)
    mirna_rows = build_mirna_questions(mirnas)

    _save_csv("simple_genes.csv", simple_gene_rows)
    _save_csv("multi_genes.csv", multi_gene_rows)
    _save_csv("variants.csv", variant_rows)
    _save_csv("drugs.csv", drug_rows)
    _save_csv("mirnas.csv", mirna_rows)

    print(f"Wrote {len(simple_gene_rows)} simple gene queries")
    print(f"Wrote {len(multi_gene_rows)} multi-gene queries")
    print(f"Wrote {len(variant_rows)} variant queries")
    print(f"Wrote {len(drug_rows)} drug queries")
    print(f"Wrote {len(mirna_rows)} miRNA queries")


if __name__ == "__main__":
    main()
