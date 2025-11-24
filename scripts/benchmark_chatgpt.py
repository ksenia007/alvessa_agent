from __future__ import annotations

import argparse
import csv
import os
import random
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import List

from openai import OpenAI

from benchmark_config import MC_BENCHMARK_SYS_MSG, MC_BENCHMARK_PROMPT_ADD, extract_answer_letter


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark ChatGPT on multiple-choice CSVs.")
    parser.add_argument("folder", help="Folder to search recursively for benchmark CSV files.")
    parser.add_argument("--N", type=int, default=0, help="Number of questions to run per CSV (default: all).")
    parser.add_argument("--shuffle", action="store_true", help="Shuffle rows per CSV (seed=42).")
    parser.add_argument(
        "--restart",
        help="Path to an existing chatgpt_baseline_*.csv; rows will be copied and questions skipped.",
    )
    return parser.parse_args()


def iter_csv_files(root: Path) -> List[Path]:
    return sorted(root.rglob("*.csv"))


def load_rows(csv_path: Path) -> List[dict]:
    with csv_path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        return [row for row in reader if row.get("question")]


def call_chatgpt(client: OpenAI, question: str) -> str:
    try:
        resp = client.chat.completions.create(
            model="gpt-5.1-2025-11-13",
            max_tokens=500,
            temperature=0,
            messages=[
                {"role": "system", "content": MC_BENCHMARK_SYS_MSG},
                {"role": "user", "content": MC_BENCHMARK_PROMPT_ADD + question},
            ],
        )
        choice = resp.choices[0] if resp.choices else None
        message = choice.message if choice else None
        text = message.content if message else ""
        return text.strip() if text else ""
    except Exception as exc:
        print(f"[chatgpt] error: {exc}", file=sys.stderr)
        return "ERROR"


def main() -> int:
    args = parse_args()
    root = Path(args.folder).expanduser().resolve()
    if not root.exists():
        print(f"Folder not found: {root}", file=sys.stderr)
        return 1

    api_key = os.environ.get("CHAT_API_KEY")
    if not api_key:
        print("CHAT_API_KEY not set", file=sys.stderr)
        return 1

    client = OpenAI(api_key=api_key)

    csv_files = iter_csv_files(root)
    if not csv_files:
        print(f"No CSV files under {root}", file=sys.stderr)
        return 1

    timestamp = datetime.now().strftime("%Y%m%d-%H%M")
    out_name = f"chatgpt_baseline_{root.name}_{timestamp}.csv"
    out_path = Path(out_name).resolve()

    header = [
        "index",
        "question",
        "expected_answer",
        "predicted_answer",
        "is_correct",
        "runtime_seconds",
        "tool_tag",
        "source_csv",
        "source_name",
        "source_folder",
    ]

    total_questions = 0
    total_correct = 0
    total_runtime = 0.0
    next_index = 1

    print(f"[chatgpt baseline] Starting benchmark at {datetime.now().isoformat()}")
    print(f"[chatgpt baseline] Output will be saved to: {out_path}")

    with out_path.open("w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)

        for csv_file in csv_files:
            print(f"[benchmark] processing file: {csv_file}", file=sys.stderr)
            rows = load_rows(csv_file)
            if not rows:
                print("[benchmark] no valid rows found, skipping", file=sys.stderr)
                raise ValueError
                continue

            if args.shuffle:
                random.Random(42).shuffle(rows)

            limit = args.N if args.N and args.N > 0 else len(rows)

            source_folder = csv_file.parent.name
            source_name = csv_file.name

            for local_idx, row in enumerate(rows):
                if local_idx >= limit:
                    break
                question = (row.get("question") or "").strip()
                if not question:
                    continue
                expected = (row.get("answer") or "").strip()
                tool_tag = row.get("tool") or ""

                try:
                    q_start = time.perf_counter()
                    answer_text = call_chatgpt(client, question)
                    q_elapsed = time.perf_counter() - q_start
                except Exception as exc:
                    print(f"[benchmark] error on question: {exc}", file=sys.stderr)
                    answer_text = "ERROR"
                    q_elapsed = 0.0
                total_runtime += q_elapsed

                pred_letter = extract_answer_letter(answer_text)
                correct_letter = expected[:1].upper() if expected else ""
                is_correct = int(bool(pred_letter and correct_letter and pred_letter == correct_letter))
                total_correct += is_correct
                total_questions += 1

                writer.writerow(
                    [
                        next_index,
                        question,
                        expected,
                        answer_text,
                        str(is_correct),
                        f"{q_elapsed:.3f}",
                        tool_tag,
                        str(csv_file),
                        source_name,
                        source_folder,
                    ]
                )
                csvfile.flush()
                next_index += 1

    accuracy = total_correct / total_questions if total_questions else 0.0
    print(f"[chatgpt baseline] Completed {total_questions} questions; accuracy={accuracy:.2%}; output={out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
