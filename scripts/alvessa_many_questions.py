"""
Batch utility: run Alvessa question mode over a CSV of {id, question} rows.

Outputs a JSON map {id: answer}. Optionally writes per-question artifacts
under out/<timestamp>_manyq/case_0001/, mirroring the CLI layout, so long runs
can be resumed/inspected without rerunning everything.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
import pandas as pd
from typing import Dict, List, Tuple

# Ensure repo root on path so `src.*` imports resolve when run from anywhere.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.alvessa.pipeline import run_pipeline  # noqa: E402
from src.alvessa.workflow.output_paths import create_run_directory, build_output_paths  # noqa: E402
from src.state import create_files_from_state  # noqa: E402


def _load_questions(csv_path: Path) -> List[Tuple[str, str]]:
    """Return list of (id, question) tuples; skip blank entries."""
    df = pd.read_csv(csv_path) 
    rows: List[Tuple[str, str]] = []
    for _, row in df.iterrows():
        qid = str(row.get("id", "")).strip()
        question = str(row.get("question", "")).strip()
        if qid and question:
            rows.append((qid, question))
    return rows


def _extract_answer(state: Dict) -> str:
    """Best-effort pull of the answer text from a pipeline state."""
    if not isinstance(state, dict):
        return ""
    llm_json = state.get("llm_json") or {}
    answer = (llm_json.get("answer") if isinstance(llm_json, dict) else None) or state.get("answer", "")
    if answer is None:
        return ""
    return str(answer).strip()


def _write_artifacts(case_dir: Path, question: str, result: Dict) -> None:
    """Persist pipeline outputs similar to CLI: demo.json/txt + rich state files."""
    case_dir.mkdir(parents=True, exist_ok=True)
    outputs = build_output_paths(case_dir)
    outputs["json"].write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")

    llm_json = result.get("llm_json", {}) or {}
    answer = llm_json.get("answer", "")
    evidence = llm_json.get("evidence", []) or []
    lines: List[str] = ["=" * 80, f"Q: {question}", "", "Answer:", str(answer), "", "Evidence:"]
    lines.extend(f"  - {item}" for item in evidence)
    outputs["txt"].write_text("\n".join(lines) + "\n", encoding="utf-8")

    # Gene/variant/drug artifacts for UI reuse
    create_files_from_state(result, case_dir)


def run_batch(csv_path: Path, *, limit: int | None = None, save_intermediate: bool = False, run_dir: Path | None = None, out_json: Path | None = None) -> Dict[str, str]:
    """Run pipeline for each question; return {id: answer} and optionally write artifacts."""
    questions = _load_questions(csv_path)
    if limit is not None and limit > 0:
        questions = questions[:limit]
        
    print(f"[alvessa] Loaded {len(questions)} question(s) from {csv_path}")

    results: Dict[str, str] = {}

    if save_intermediate:
        if run_dir is None:
            run_dir, _ = create_run_directory("manyq")
        run_dir.mkdir(parents=True, exist_ok=True)
        print(f"[alvessa] Saving per-question artifacts under {run_dir}")

    for idx, (qid, question) in enumerate(questions, start=1):
        print(f"[alvessa] ({idx}/{len(questions)}) {qid}: {question[:80]}{'...' if len(question) > 80 else ''}")
        try:
            state = run_pipeline(question)
            answer = _extract_answer(state)
        except Exception as exc:  # keep batch going
            print(f"[alvessa] ERROR for {qid}: {exc}", file=sys.stderr)
            state = {}
            answer = ""

        results[qid] = answer

        if save_intermediate and run_dir:
            case_dir = run_dir / f"case_{idx:04d}"
            _write_artifacts(case_dir, question, state)

        if out_json:
            out_json.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")

    return results


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Batch runner for Alvessa question mode.")
    parser.add_argument("csv_path", help="Input CSV with columns: id, question.")
    parser.add_argument(
        "-o",
        "--output",
        help="Output JSON path (default: <csv_stem>_answers.json in same folder).",
    )
    parser.add_argument(
        "-n",
        "--limit",
        type=int,
        default=0,
        help="Optional cap on number of questions (0 = all).",
    )
    parser.add_argument(
        "--save-intermediate",
        action="store_true",
        help="If set, write per-question artifacts under out/<timestamp>_manyq/case_xxxx/.",
    )
    parser.add_argument(
        "--run-dir",
        help="Optional run directory root (used only if --save-intermediate).",
    )
    args = parser.parse_args(argv)

    csv_path = Path(args.csv_path).expanduser().resolve()
    if not csv_path.exists():
        print(f"[alvessa] CSV not found: {csv_path}", file=sys.stderr)
        return 1

    limit = args.limit if args.limit and args.limit > 0 else None
    out_path = (
        Path(args.output).expanduser().resolve()
        if args.output
        else csv_path.with_name(f"{csv_path.stem}_answers.json")
    )
    print(f"[alvessa] Writing answers to: {out_path}")
    run_dir = Path(args.run_dir).expanduser().resolve() if args.run_dir else None

    results = run_batch(
        csv_path,
        limit=limit,
        save_intermediate=bool(args.save_intermediate),
        run_dir=run_dir,
        out_json=out_path,
    )

    out_path.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[alvessa] Wrote answers for {len(results)} question(s) to {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
