from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import pandas as pd


def subset_dbqa(input_path: Path, output_path: Path, n: int = 100, seed: int = 42) -> None:
    df = pd.read_csv(input_path)
    if df.empty:
        raise ValueError("Input CSV is empty.")

    rng = random.Random(seed)
    if len(df) <= n:
        sampled = df.copy()
    else:
        sampled = df.sample(n=n, random_state=seed)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    sampled.to_csv(output_path, index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description="Create a random subset of dbQA questions.")
    parser.add_argument(
        "--input",
        default="/Users/sokolova/Documents/research/alvessa_agent/benchmarks_generation/questions/dbQA/dbqa_mc.csv",
        help="Path to the full dbQA CSV.",
    )
    parser.add_argument(
        "--output",
        default="/Users/sokolova/Documents/research/alvessa_agent/benchmarks_generation/questions/dbQA/dbqa_mc_SUBSET.csv",
        help="Path to write the subset CSV.",
    )
    parser.add_argument("--n", type=int, default=100, help="Number of questions to sample (default: 100).")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42).")
    args = parser.parse_args()

    try:
        subset_dbqa(Path(args.input), Path(args.output), n=args.n, seed=args.seed)
    except Exception as exc:
        print(f"Failed to subset dbQA: {exc}", file=sys.stderr)
        return 1
    print(f"Saved subset to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

