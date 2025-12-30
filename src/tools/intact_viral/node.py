"""
Annotate genes with host–virus interactions from IntAct (viral dataset with gene names).

Reads local_dbs/intact_viral_with_genes.txt and, for each gene in the state,
finds rows where the gene symbol appears in Gene Names A or Gene Names B. For matching rows,
emits a text summary listing the viral partner gene names (the opposite interactor) and taxid A/B.
No structured data is stored; only text summaries with prefix *IntActViral: ...
"""

from __future__ import annotations

import pandas as pd
from pathlib import Path
from typing import Dict

from src.state import State
from src.tools.base import Node

DEBUG = True

REPO_ROOT = Path(__file__).resolve().parents[3]
LOCAL_DBS = REPO_ROOT / "local_dbs"
INTACT_VIRAL_PATH = LOCAL_DBS / "intact_viral_with_genes.txt"

COL_GENE_A = "Gene Names A"
COL_GENE_B = "Gene Names B"
COL_TAX_A = "Taxid interactor A"
COL_TAX_B = "Taxid interactor B"


def intact_viral_agent(state: "State") -> "State":
    gene_entities = state.get("gene_entities") or {}
    if not gene_entities:
        return state

    if not INTACT_VIRAL_PATH.exists():
        if DEBUG:
            print(f"[IntActViral] Missing file: {INTACT_VIRAL_PATH}")
        return state

    try:
        df = pd.read_csv(INTACT_VIRAL_PATH, sep="\t", dtype=str, low_memory=False)
        df = df.rename(columns=lambda c: str(c).strip().lstrip("\ufeff"))
    except Exception as exc:
        if DEBUG:
            print(f"[IntActViral] Failed to load: {exc}")
        return state

    for sym, g in gene_entities.items():
        if not sym:
            continue
        sym_up = sym.upper()
        print(f"[IntActViral] Searching interactions for gene: {sym}")
        # Match gene name in either column
        print(COL_GENE_A, COL_GENE_B)
        # print columns
        print(df.columns.tolist())
        mask_a = df[COL_GENE_A].fillna("").str.contains(sym_up, case=False) 
        mask_b = df[COL_GENE_B].fillna("").str.contains(sym_up, case=False, regex=True) 
        print(f"[IntActViral] Processing {sym}: found {mask_a.sum() if len(mask_a) else 0} in A, {mask_b.sum() if len(mask_b) else 0} in B")
        matches = df[mask_a | mask_b] if len(mask_a) or len(mask_b) else pd.DataFrame()
        if matches.empty:
            continue

        rows = []
        for _, row in matches.iterrows():
            taxA = row.get(COL_TAX_A, "").split('|')[0]
            taxB = row.get(COL_TAX_B, "").split('|')[0]
            genesA = row.get(COL_GENE_A, "")
            genesB = row.get(COL_GENE_B, "")
            # also add Alias(es) interactor B and Alias(es) interactor A
            aliasA = row.get("Alias(es) interactor A", "")
            aliasB = row.get("Alias(es) interactor B", "")
            if sym_up in str(genesA).upper():
                partner = str(genesB)
                alias_partner = str(aliasB)
                taxa_partner = taxB
                taxa_main = taxA
            else:
                partner = str(genesA)
                alias_partner = str(aliasA)
                taxa_partner = taxA
                taxa_main = taxB
            rows.append(f"{partner}, full alias list {alias_partner} (taxa_this={taxa_main}, tax_parter={taxa_partner})")

        if not rows:
            continue

        summary = f"*IntActViral: Viral interactions for {sym} (listing parter/interactor gene synonyms, alias info for this partner and then taxa for input gene, followed by taxa for the partner): " + "; ".join(rows) + "."
        g.update_text_summaries(summary)

    return state


NODES: tuple[Node, ...] = (
    Node(
        name="intact_viral",
        entry_point=intact_viral_agent,
        description="Fetches curated host–virus molecular interactions from the IntAct Virus dataset, reporting experimentally supported protein–protein interactions between host genes/proteins and viral proteins, with organism and taxon context for each interactor.",
    ),
)


if __name__ == "__main__":
    from src.alvessa.domain.gene_class import Gene
    from src.alvessa.domain.gene_components import GeneIdentifiers

    test_genes = ["ACE2", "BRCA2"]
    state = {"gene_entities": {g: Gene(GeneIdentifiers(symbol=g)) for g in test_genes}}
    out = intact_viral_agent(state)
    for g in test_genes:
        obj = out["gene_entities"][g]
        print(f"{g}: summaries:")
        for s in obj.text_summaries_from_tools:
            if "IntActViral" in s:
                print("  ", s)
        print("-" * 20)
