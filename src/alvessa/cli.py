"""Command-line entrypoint for running Alvessa workflows."""

from __future__ import annotations

import argparse
import csv
import json
import random
import sys
import textwrap
import time
from datetime import datetime
from pathlib import Path
from typing import Iterable, List

from src.alvessa.pipeline import run_pipeline
from src.alvessa.ui.server import serve as serve_ui
from src.alvessa.workflow.output_paths import build_output_paths, create_run_directory
from src.tools.base import Node
from src.tools.registry import NODES
from src.state import create_files_from_state

MC_BENCHMARK_PROMPT = (
    "You are answering multiple-choice questions. Each question lists answer choices "
    "labeled [A], [B], [C], [D]. Respond with exactly one capital letter (A, B, C, or D)."
    " Do not provide any explanations or additional text."
    " IMPORTANT: Response length should always be 1 character, independent of the question or provided information."
)


def _run_benchmark_rows(
    rows: List[dict],
    *,
    csv_path: Path,
    limit: int,
    run_dir: Path,
    writer: csv.writer,
    csv_handle,
    start_index: int,
    save_intermediate: bool,
) -> tuple[int, int, float, int]:
    """Execute up to ``limit`` question rows, appending metadata to the CSV writer."""
    processed = 0
    correct = 0
    total_runtime = 0.0
    next_index = start_index
    max_questions = len(rows) if limit <= 0 else min(limit, len(rows))

    for local_idx, row in enumerate(rows):
        if processed >= max_questions:
            break
        question = (row.get("question") or "").strip()
        if not question:
            continue
        expected = (row.get("answer") or "").strip()
        tool_tag = row.get("tool") or ""

        label = f"{csv_path}:{local_idx + 1}"
        print(f"[alvessa][{processed + 1}/{max_questions}] {label} — {question[:90]}{'...' if len(question) > 90 else ''}")

        case_dir = run_dir / f"case_{next_index:04d}" if save_intermediate else None
        if case_dir:
            case_dir.mkdir(parents=True, exist_ok=True)

        q_start = time.perf_counter()
        result = run_pipeline(
            question,
            prompt=MC_BENCHMARK_PROMPT,
            mc_setup=True,
            output_dir=case_dir if case_dir else run_dir,
        )
        q_elapsed = time.perf_counter() - q_start
        total_runtime += q_elapsed

        if save_intermediate and case_dir:
            _write_artifacts(case_dir, question, result)

        llm_json = result.get("llm_json", {}) or {}
        answer = (llm_json.get("answer") or "").strip()
        
        pred_letter = answer[-1:].upper() if answer else "" # use last letter
        correct_letter = expected[:1].upper() if expected else ""
        is_correct = bool(pred_letter and correct_letter and pred_letter == correct_letter)
        if is_correct:
            correct += 1

        used_tools = (
            result.get("used_tools")
            or llm_json.get("used_tools")
            or []
        )

        source_folder = csv_path.parent.name if csv_path.parent.name else ""
        writer.writerow(
            [
                next_index,
                question,
                expected,
                answer,
                "1" if is_correct else "0",
                f"{q_elapsed:.3f}",
                tool_tag,
                ";".join(used_tools) if used_tools else "",
                case_dir.name if case_dir else "",
                str((case_dir / "demo.json").relative_to(run_dir)) if case_dir else "",
                str(csv_path),
                csv_path.name,
                source_folder,
            ]
        )
        csv_handle.flush()

        processed += 1
        next_index += 1

    return processed, correct, total_runtime, next_index

def _write_artifacts(run_dir: Path, question: str, result: dict) -> None:
    """Persist canonical JSON/text artifacts for a pipeline run + the TABLES"""
    outputs = build_output_paths(run_dir)

    outputs["json"].write_text(
        json.dumps(result, indent=2, default=str),
        encoding="utf-8",
    )

    llm_json = result.get("llm_json", {}) or {}
    answer = llm_json.get("answer", "")
    evidence = llm_json.get("evidence", []) or []

    lines: List[str] = ["=" * 80, f"Q: {question}", "", "Answer:", str(answer), "", "Evidence:"]
    lines.extend(f"  - {item}" for item in evidence)
    outputs["txt"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    
    create_files_from_state(result, run_dir)


def _render_table(nodes: Iterable[Node]) -> str:
    """Return a formatted table describing the registered tool nodes."""
    rows = []
    for node in sorted(nodes, key=lambda n: n.name.lower()):
        rows.append(
            {
                "Name": node.name,
                "Description": (node.description or "").strip(),
                "Dependencies": ", ".join(node.dependencies) if node.dependencies else "—",
                "Aliases": ", ".join(node.aliases) if node.aliases else "—",
                "Tags": ", ".join(node.tags) if node.tags else "—",
            }
        )

    wrap_width = 72
    if not rows:
        return "No tools registered."
    col_widths = {
        "Name": max(len("Name"), *(len(row["Name"]) for row in rows)) + 2,
        "Description": wrap_width + 2,
        "Dependencies": max(len("Dependencies"), *(len(row["Dependencies"]) for row in rows)) + 2,
        "Aliases": max(len("Aliases"), *(len(row["Aliases"]) for row in rows)) + 2,
        "Tags": max(len("Tags"), *(len(row["Tags"]) for row in rows)) + 2,
    }

    def wrap(value: str, width: int) -> List[str]:
        text = value.strip()
        if not text:
            return [""]
        return textwrap.wrap(text, width=width) or [""]

    header = "".join(
        (
            "Name".ljust(col_widths["Name"]),
            "Description".ljust(col_widths["Description"]),
            "Dependencies".ljust(col_widths["Dependencies"]),
            "Aliases".ljust(col_widths["Aliases"]),
            "Tags".ljust(col_widths["Tags"]),
        )
    )
    divider = "-" * len(header)

    lines: List[str] = [header, divider]

    for row in rows:
        name_lines = [row["Name"]]
        desc_lines = wrap(row["Description"], wrap_width)
        dep_lines = wrap(row["Dependencies"], col_widths["Dependencies"] - 2)
        alias_lines = wrap(row["Aliases"], col_widths["Aliases"] - 2)
        tag_lines = wrap(row["Tags"], col_widths["Tags"] - 2)
        max_height = max(len(name_lines), len(desc_lines), len(dep_lines), len(alias_lines), len(tag_lines))

        for idx in range(max_height):
            line = (name_lines[idx] if idx < len(name_lines) else "")
            line = line.ljust(col_widths["Name"])
            line += (desc_lines[idx] if idx < len(desc_lines) else "").ljust(col_widths["Description"])
            line += (dep_lines[idx] if idx < len(dep_lines) else "").ljust(col_widths["Dependencies"])
            line += (alias_lines[idx] if idx < len(alias_lines) else "").ljust(col_widths["Aliases"])
            line += (tag_lines[idx] if idx < len(tag_lines) else "").ljust(col_widths["Tags"])
            lines.append(line)
        lines.append("")

    return "\n".join(lines).rstrip()


def _handle_question(args: argparse.Namespace) -> int:
    """Execute the LangGraph workflow for a single question."""
    query = args.query.strip()
    if not query:
        print("[alvessa] Please provide a non-empty question.", file=sys.stderr)
        return 1

    run_dir, _ = create_run_directory("cli")
    print(f"[alvessa] Running pipeline for: {query}\n[alvessa] Writing artifacts under {run_dir}")
    result = run_pipeline(query, output_dir=run_dir)

    _write_artifacts(run_dir, query, result)

    llm_json = result.get("llm_json", {}) or {}
    answer = llm_json.get("answer", "")

    print("\n=== Answer ===")
    print(answer if answer else "(empty response)")
    
    verification = result.get("verification", "")
    print("\n=== Evidence ===")
    print(verification)

    return 0


def _handle_ui(args: argparse.Namespace) -> int:
    """Launch the local FastAPI UI using uvicorn."""
    port = args.port if args.port is not None else 8000
    serve_ui(port=port)
    return 0


def _handle_tools(_: argparse.Namespace) -> int:
    """Print the available tool catalog in a human-readable table."""
    table = _render_table(NODES)
    print(table)
    return 0


def _handle_benchmark(args: argparse.Namespace) -> int:
    """Run a multiple-choice benchmark from a CSV file."""
    csv_path = Path(args.csv_path)
    if not csv_path.exists():
        print(f"[alvessa] Benchmark file not found: {csv_path}", file=sys.stderr)
        return 1

    with csv_path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        rows = [row for row in reader if row.get("question")]

    if not rows:
        print(f"[alvessa] No questions found in {csv_path}", file=sys.stderr)
        return 1

    restart_questions: List[str] = []
    restart_arg = getattr(args, "restart", None)
    if restart_arg:
        restart_path = Path(restart_arg)
        if not restart_path.exists():
            print(f"[alvessa] Restart file not found: {restart_path}", file=sys.stderr)
            return 1
        with restart_path.open("r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            restart_questions = [
                (row.get("question") or "").strip()
                for row in reader
                if row.get("question")
            ]

    if restart_questions:
        trimmed_rows: List[dict] = []
        skip_count = 0
        restart_idx = 0
        resume = False
        for row in rows:
            question = (row.get("question") or "").strip()
            if not resume and restart_idx < len(restart_questions):
                if question == restart_questions[restart_idx]:
                    restart_idx += 1
                    skip_count += 1
                    continue
                else:
                    resume = True
            trimmed_rows.append(row)
        rows = trimmed_rows
        print(f"[alvessa] Restart mode: skipped {skip_count} previously completed question(s).")

    restart_questions: List[str] = []
    if getattr(args, "restart", None):
        restart_path = Path(args.restart)
        if not restart_path.exists():
            print(f"[alvessa] Restart file not found: {restart_path}", file=sys.stderr)
            return 1
        with restart_path.open("r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            restart_questions = [
                (row.get("question") or "").strip()
                for row in reader
                if row.get("question")
            ]

    if restart_questions:
        trimmed = []
        skip_idx = 0
        skipped = 0
        resume = False
        for row in rows:
            question = (row.get("question") or "").strip()
            if not resume and skip_idx < len(restart_questions):
                if question == restart_questions[skip_idx]:
                    skip_idx += 1
                    skipped += 1
                    continue
                else:
                    resume = True
            trimmed.append(row)
        rows = trimmed
        print(f"[alvessa] Restart: skipped {skipped} completed question(s).")

    limit = args.N if args.N and args.N > 0 else len(rows)

    run_dir, _ = create_run_directory("cli")
    preview_count = min(limit, len(rows))
    print(f"[alvessa] Benchmarking up to {preview_count} question(s) from {csv_path}")
    print(f"[alvessa] Writing benchmark artifacts under {run_dir}")

    summary_csv = run_dir / "benchmark_summary.csv"
    header = [
        "index",
        "question",
        "expected_answer",
        "predicted_answer",
        "is_correct",
        "runtime_seconds",
        "tool_tag",
        "used_tools",
        "case_dir",
        "raw_output",
        "source_csv",
        "source_name",
        "source_folder",
    ]

    metadata_txt = run_dir / "metadata.txt"
    metadata_txt.write_text(
        "Benchmark metadata\n"
        f"Prompt:\n{MC_BENCHMARK_PROMPT}\n"
        f"Source file: {csv_path}\n"
        f"N per file: {limit}\n"
        f"Save intermediate: {bool(getattr(args, 'save_intermediate', False))}\n",
        encoding="utf-8",
    )

    with summary_csv.open("w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        processed, n_correct, total_runtime, _ = _run_benchmark_rows(
            rows,
            csv_path=csv_path,
            limit=limit,
            run_dir=run_dir,
            writer=writer,
            csv_handle=csvfile,
            start_index=1,
            save_intermediate=getattr(args, "save_intermediate", False),
        )

    outputs = build_output_paths(run_dir)
    accuracy = n_correct / processed if processed else 0.0
    summary = {
        "benchmark_file": str(csv_path),
        "num_questions": processed,
        "num_correct": n_correct,
        "accuracy": accuracy,
        "prompt": MC_BENCHMARK_PROMPT,
        "total_runtime_seconds": total_runtime,
        "summary_csv": str(summary_csv.relative_to(run_dir)),
        "timestamp": datetime.utcnow().isoformat() + "Z",
    }
    outputs["json"].write_text(json.dumps(summary, indent=2), encoding="utf-8")
    outputs["txt"].write_text(
        "\n".join(
            [
                f"Benchmark file: {csv_path}",
                f"Questions: {processed}",
                f"Correct: {n_correct}",
                f"Accuracy: {accuracy:.2%}",
                f"Metadata file: {metadata_txt.name}",
                f"Summary CSV: {summary_csv.name}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"[alvessa] Accuracy: {accuracy:.2%} ({n_correct}/{processed})")
    return 0


def _handle_benchmark_all(args: argparse.Namespace) -> int:
    """Run the benchmark across every CSV inside a folder."""
    folder = Path(args.folder)
    if not folder.exists():
        print(f"[alvessa] Folder not found: {folder}", file=sys.stderr)
        return 1

    csv_files = sorted(folder.rglob("*.csv"))
    if not csv_files:
        print(f"[alvessa] No CSV files under {folder}", file=sys.stderr)
        return 1

    run_dir, _ = create_run_directory("cli")
    print(f"[alvessa] Benchmark-all across {len(csv_files)} CSV files under {folder}")
    print(f"[alvessa] Writing benchmark artifacts under {run_dir}")

    summary_csv = run_dir / "benchmark_summary.csv"
    metadata_txt = run_dir / "metadata.txt"
    header = [
        "index",
        "question",
        "expected_answer",
        "predicted_answer",
        "is_correct",
        "runtime_seconds",
        "tool_tag",
        "used_tools",
        "case_dir",
        "raw_output",
        "source_csv",
        "source_name",
        "source_folder",
    ]

    metadata_txt.write_text(
        "Benchmark-all metadata\n"
        f"Prompt:\n{MC_BENCHMARK_PROMPT}\n"
        f"Root folder: {folder}\n"
        f"N per file: {args.N}\n"
        f"Save intermediate: {bool(getattr(args, 'save_intermediate', False))}\n",
        encoding="utf-8",
    )

    total_questions = 0
    total_correct = 0
    total_runtime = 0.0
    next_index = 1

    with summary_csv.open("w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)

        for csv_file in csv_files:
            with csv_file.open("r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh)
                rows = [row for row in reader if row.get("question")]

            if getattr(args, "shuffle", False):
                random.Random(42).shuffle(rows)

            limit = args.N if args.N and args.N > 0 else len(rows)
            processed, correct, runtime, next_index = _run_benchmark_rows(
                rows,
                csv_path=csv_file,
                limit=limit,
                run_dir=run_dir,
                writer=writer,
                csv_handle=csvfile,
                start_index=next_index,
                save_intermediate=getattr(args, "save_intermediate", False),
            )
            total_questions += processed
            total_correct += correct
            total_runtime += runtime

    outputs = build_output_paths(run_dir)
    accuracy = total_correct / total_questions if total_questions else 0.0
    summary = {
        "root_folder": str(folder),
        "benchmark_files": len(csv_files),
        "num_questions": total_questions,
        "num_correct": total_correct,
        "accuracy": accuracy,
        "prompt": MC_BENCHMARK_PROMPT,
        "total_runtime_seconds": total_runtime,
        "summary_csv": str(summary_csv.relative_to(run_dir)),
        "timestamp": datetime.utcnow().isoformat() + "Z",
    }
    outputs["json"].write_text(json.dumps(summary, indent=2), encoding="utf-8")
    outputs["txt"].write_text(
        "\n".join(
            [
                f"Root folder: {folder}",
                f"CSV files: {len(csv_files)}",
                f"Questions: {total_questions}",
                f"Correct: {total_correct}",
                f"Accuracy: {accuracy:.2%}",
                f"Metadata file: {metadata_txt.name}",
                f"Summary CSV: {summary_csv.name}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"[alvessa] Accuracy: {accuracy:.2%} ({total_correct}/{total_questions})")
    return 0


def _build_parser() -> argparse.ArgumentParser:
    """Construct the top-level argument parser."""
    parser = argparse.ArgumentParser(prog="alvessa", description="CLI for the Alvessa agent.")
    subparsers = parser.add_subparsers(dest="command")

    question_parser = subparsers.add_parser("question", help="Run the pipeline on a single question.")
    question_parser.add_argument("query", help="Question to submit to the agent.")
    question_parser.set_defaults(func=_handle_question)

    ui_parser = subparsers.add_parser("ui", help="Serve the local FastAPI dashboard.")
    ui_parser.add_argument(
        "port",
        nargs="?",
        type=int,
        default=8000,
        help="Port to bind (default: 8000).",
    )
    ui_parser.set_defaults(func=_handle_ui)

    tools_parser = subparsers.add_parser("tools", help="List the registered tool catalog.")
    tools_parser.set_defaults(func=_handle_tools)

    bench_parser = subparsers.add_parser("benchmark", help="Run a CSV multiple-choice benchmark.")
    bench_parser.add_argument("csv_path", help="Path to the benchmark CSV file.")
    bench_parser.add_argument(
        "--N",
        type=int,
        default=0,
        help="Number of questions to run (default: all remaining questions).",
    )
    bench_parser.add_argument("--restart", help="Path to an existing benchmark_summary.csv to resume from.")
    bench_parser.add_argument(
        "--save_intermediate",
        action="store_true",
        help="If set, save per-question case_* artifacts (default: off).",
    )
    bench_parser.set_defaults(func=_handle_benchmark)

    bench_all_parser = subparsers.add_parser("benchmark_all", help="Run benchmarks for every CSV under a folder.")
    bench_all_parser.add_argument("folder", help="Folder containing benchmark CSV files (searched recursively).")
    bench_all_parser.add_argument(
        "--N",
        type=int,
        default=0,
        help="Number of questions to run per CSV (default: all).",
    )
    bench_all_parser.add_argument(
        "--save_intermediate",
        action="store_true",
        help="If set, save per-question case_* artifacts for every CSV.",
    )
    bench_all_parser.add_argument(
        "--shuffle",
        action="store_true",
        help="Shuffle each CSV's rows (seed=42) before selecting N.",
    )
    bench_all_parser.set_defaults(func=_handle_benchmark_all)

    return parser


def main(argv: Iterable[str] | None = None) -> int:
    """Entry point used by `python -m alvessa.cli` and the `alvessa` script."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 1
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover - manual entry point
    sys.exit(main())
