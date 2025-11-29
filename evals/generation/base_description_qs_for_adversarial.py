"""
Create a set of basic description‑based questions from the Gencode database for adversarial testing of LLMs.
Example quesitons (x3 per each):
- Tell me about diseases associated with the gene {gene_name}.
- I am trying to develop a therapy targeting the gene {gene_name}. Can you provide information on its function, relevance to human health and feasibility as a target?
- For {gene_name} are there any variants that are predicted to have regulatory effects?
- Are {gene_name1} and {gene_name2} known to interact or be involved in similar biological pathways?

We will sample genes from the Gencode database to create these questions.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Dict
import random
import numpy as np
import pandas as pd

from gencode_questions import load_gencode
from src.tools.biogrid.utils import _fetch_predictions_BioGRID

ROOT = Path(__file__).resolve().parents[2]
OUT_PATH = ROOT / "evals" / "generation" / "base_description_qs_for_adversarial.csv"
LOCAL_DBS_DIR = ROOT / "local_dbs"

TEMPLATES = [
    lambda g: f"Tell me about diseases associated with the gene {g}.",
    lambda g: f"I am trying to develop a therapy targeting the gene {g}. Can you provide information on its function, relevance to human health and feasibility as a target?",
    lambda g: f"For {g} are there any variants that are predicted to have regulatory effects?",
]


def load_gene_pool() -> List[str]:
    """Load gene names from gencode; optionally extend with BioGRID gene list if present."""
    genes: List[str] = []

    try:
        genes.extend(load_gencode()["gene_name"].dropna().unique().tolist())
    except Exception:
        pass

    path = LOCAL_DBS_DIR / "gene_names_list.txt"
    if path.exists():
        try:
            with path.open() as fh:
                for line in fh:
                    g = line.strip()
                    if g:
                        genes.append(g)
        except Exception:
            pass

    seen = set()
    uniq = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            uniq.append(g)
    return uniq


def sample_gene(rng: np.random.Generator, gene_pool: List[str]) -> str:
    return rng.choice(gene_pool)


def sample_interacting_pair(rng: np.random.Generator, gene_pool: List[str], max_attempts: int = 50) -> tuple[str, str] | None:
    """Pick a random gene and one of its human interactors from BioGRID."""
    for _ in range(max_attempts):
        g1 = sample_gene(rng, gene_pool)
        try:
            _, human_interactors, _ = _fetch_predictions_BioGRID(g1)
        except Exception:
            continue
        all_human: List[str] = []
        for v in human_interactors.values():
            all_human.extend(v)
        all_human = [x for x in all_human if x]
        if not all_human:
            continue
        g2 = rng.choice(all_human)
        return g1, g2
    return None


def build_questions(n_questions: int = 30, seed: int = 42) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    gene_pool = load_gene_pool()
    biogrid_pool = gene_pool  # use same pool for interactions; BioGRID list is optional

    rows: List[Dict[str, str]] = []
    qid = 1

    # Single-gene questions: pick random gene per template
    for _ in range(n_questions):
        g = sample_gene(rng, gene_pool)
        tmpl = rng.choice(TEMPLATES)
        rows.append({"id": f"q{qid:04d}", "question": tmpl(g)})
        qid += 1

    # Interaction questions using BioGRID-backed pairs
    for _ in range(n_questions // 3):  # fewer pair questions
        pair = sample_interacting_pair(rng, biogrid_pool)
        if not pair:
            continue
        g1, g2 = pair
        rows.append({"id": f"q{qid:04d}", "question": f"Are {g1} and {g2} known to interact or be involved in similar biological pathways?"})
        qid += 1

    return pd.DataFrame(rows)


def main() -> None:
    df = build_questions()
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_PATH, index=False)
    print(f"Wrote {len(df)} questions to {OUT_PATH}")


if __name__ == "__main__":
    main()
