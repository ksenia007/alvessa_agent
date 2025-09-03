"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-08-21
Updated: 2025-08-21


Description: 

Run labBench & save results."""

import argparse
import os
from pathlib import Path

from evaluate_benchmarks_general import run_benchmark
SUBFOLDERS = ['labbench']
MAIN_PATH = Path("benchmarks_generation")
MAX_ROWS = 100 # -1 means all

SYSTEM_MSG = (
    """You are a multiple-choice answering system.
You must reply with exactly one of the following letters: A, B, C, or D.
Do not include any explanation, reasoning, or extra text.
Your response will be parsed by a program that will fail if you output anything other than a single capital letter.
Example valid output: C
Example invalid outputs: "Answer: C", "C.", "Option C", "B because...", "A is correct because..." """
)

def main():
    parser = argparse.ArgumentParser(description="Run benchmarks with a single switch.")
    parser.add_argument("runner", choices=["claude", "alvessa"],
                        help="Pick which pipeline to run.")
    args = parser.parse_args()

    for subfolder in SUBFOLDERS:
        sub_path = MAIN_PATH / subfolder
        if not sub_path.is_dir():
            print(f"[warn] Skipping missing subfolder: {sub_path}")
            continue

        out_dir = MAIN_PATH / "results" / subfolder / args.runner
        out_dir.mkdir(parents=True, exist_ok=True)

        for file_name in sorted(os.listdir(sub_path)):
            if not file_name.endswith(".csv"):
                continue

            if 'dbqa_' not in file_name:
                continue
            
            in_path = sub_path / file_name
            out_path = out_dir / file_name

            if args.runner == "claude":
                print(f"Running Claude on {in_path}")
                run_benchmark(
                    str(in_path),
                    system_msg=SYSTEM_MSG,
                    max_rows=MAX_ROWS,
                    file_save=str(out_path),
                    method='claude'
                )
            else:  # alvessa
                print(f"Running Alvessa on {in_path}")
                run_benchmark(
                    str(in_path),
                    system_msg=SYSTEM_MSG,
                    max_rows=MAX_ROWS,
                    file_save=str(out_path),
                    method='alvessa'
                )

if __name__ == "__main__":
    main()
