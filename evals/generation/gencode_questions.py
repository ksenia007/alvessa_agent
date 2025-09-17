"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-08-22
Updated: 2025-08-22


Description: 

Create an evaluation using Gencode database.


"""
#%%
from pathlib import Path

import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
LOCAL_DBS_DIR = ROOT / "local_dbs"


def load_gencode(file_path: Path = LOCAL_DBS_DIR / "gencode.v48.expanded.parquet") -> pd.DataFrame:
    """Load Gencode database from a CSV file."""
    df = pd.read_parquet(file_path)
    df = df[df.gene_type.isin(['protein_coding', 'lncRNA'])]
    return df



def get_counts(df: pd.DataFrame) -> pd.DataFrame:
    """
    Get counts of exons per gene name as a table to generate questions
    
    Args:
        df: DataFrame containing Gencode data.
        
    Returns:
        DataFrame with counts of the specified feature type.
    """
    ## Work only with exon features
    exons = df[df["feature"] == "exon"].copy()

    # Per‑gene quick summary
    gene_summary = (
        exons.groupby("gene_id")
            .agg(
                gene_name=("gene_name", "first"),
                num_transcripts=("transcript_id", "nunique"),
                num_exons=("exon_number", "count"),  # total exon rows across all transcripts
            )
            .reset_index()
    )

    # Per‑transcript exon counts (within each gene)
    transcript_exon_counts = (
        exons.groupby(["gene_id", "transcript_id"])
            .agg(exon_count=("exon_number", "count"))
            .reset_index()
    )

    # Gene‑level stats across transcripts
    gene_exon_stats = (
        transcript_exon_counts.groupby("gene_id")
            .agg(
                max_exons=("exon_count", "max"),
                median_exons=("exon_count", "median"),
                mean_exons=("exon_count", "mean"),
            )
            .reset_index()
    )

    # Merge everything
    summary = gene_summary.merge(gene_exon_stats, on="gene_id", how="left")
    return summary


def make_mc_questions_max_exons(df: pd.DataFrame, n: int = 20, seed: int = 42) -> pd.DataFrame:
    """ MCQ: 'What is the maximum # exons <gene_name> can have across all transcripts?'
    - Sample n genes, build 4 options: 1 correct (max_exons) + 3 random wrongs from 1..15.
    - Shuffle A–D; return columns: question, answer, tool.
    """
    rng = np.random.default_rng(seed)

    sample_df = df.sample(n=n, random_state=seed)

    out = []
    for _, row in sample_df.iterrows():
        gene_id = row["gene_name"]
        correct = int(row["max_exons"])

        # 3 wrong answers from 1..15, != correct
        pool = [x for x in range(1, 16) if x != correct]
        wrongs = rng.choice(pool, size=3, replace=False).tolist()

        # shuffle options
        options = wrongs + [correct]
        rng.shuffle(options)

        letters = ["A", "B", "C", "D"]
        correct_letter = letters[options.index(correct)]

        # build question text WITH options inline
        opts_text = "  ".join(f"[{letters[i]}] {options[i]}" for i in range(4))
        q = f"What is the maximum # exons {gene_id} can have across all transcripts? {opts_text}"

        out.append({
            "question": q,
            "answer": f"{correct_letter}",
            "tool": "gencode_gene_node",
        })

    return pd.DataFrame(out)


def make_mc_questions_most_exons_in_list(df: pd.DataFrame, n_questions: int = 20, 
                                 seed: int = 42, max_attempts_per_q: int = 200) -> pd.DataFrame:
    """
    MCQ: 'Out of the following genes, which one has the most number of exons?'
    - Sample 6 genes, require a unique highest max_exons among those 6.
    - Build 4 options from: 1 correct (the unique max) + 3 random others from the remaining 5.
    - Shuffle A–D; return columns: question, answer, tool.
    """
    rng = np.random.default_rng(seed)

    # label to display (prefer gene_name; fallback to gene_id)
    show = df["gene_name"].dropna()
    pool = df.assign(__label__=show)

    # need at least 6 and some diversity in max_exons
    if len(pool) < 6 or pool["max_exons"].nunique() < 2:
        raise ValueError("Dataset too small or lacks diversity in max_exons.")

    letters = ["A", "B", "C", "D"]
    out = []

    for _ in range(n_questions):
        for _attempt in range(max_attempts_per_q):
            six = pool.sample(n=6, random_state=rng.integers(0, 1_000_000)).reset_index(drop=True)

            # check uniqueness of the maximum within these 6
            vmax = six["max_exons"].max()
            if (six["max_exons"] == vmax).sum() != 1:
                continue  # resample if top is tied

            # identify the correct row (unique max)
            correct_idx = int(six["max_exons"].idxmax())
            correct_label = six.loc[correct_idx, "__label__"]

            # choose 3 distractors from the remaining 5
            remaining = six.drop(index=correct_idx)
            distractors = remaining.sample(n=3, random_state=rng.integers(0, 1_000_000))["__label__"].tolist()

            # assemble and shuffle options
            options = distractors + [correct_label]
            rng.shuffle(options)

            # correct letter
            correct_letter = letters[options.index(correct_label)]

            # build question text
            opts_text = "  ".join(f"[{letters[i]}] {options[i]}" for i in range(4))
            q = f"Out of the following genes, which one has the most number of exons? {opts_text}"

            out.append({
                "question": q,
                "answer": f"{correct_letter}",
                "tool": "gencode_gene_node",
            })
            break
        else:
            raise RuntimeError("Failed to build a tie-free question after many attempts. Check data diversity or increase max_attempts_per_q.")

    return pd.DataFrame(out)


if __name__ == "__main__":
    df = load_gencode()
    df = get_counts(df)
    mc_max_exons = make_mc_questions_max_exons(df)
    mc_max_exons.to_csv('gencode/set1.csv', index=False)
    mc_select_max_list = make_mc_questions_most_exons_in_list(df)
    mc_select_max_list.to_csv('gencode/set2.csv', index=False)
