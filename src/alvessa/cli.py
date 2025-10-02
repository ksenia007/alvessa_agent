"""Command-line entrypoint for running Alvessa workflows."""

from __future__ import annotations

import argparse
import json
import sys
import textwrap
from pathlib import Path
from typing import Iterable, List

from src.alvessa.pipeline import run_pipeline
from src.alvessa.ui.server import serve as serve_ui
from src.alvessa.workflow.output_paths import build_output_paths, create_run_directory
from src.tools.base import Node
from src.tools.registry import NODES
from src.state import create_files_from_state

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
    evidence = llm_json.get("evidence", []) or []

    print("\n=== Answer ===")
    print(answer if answer else "(empty response)")
    if evidence:
        print("\n=== Evidence ===")
        for item in evidence:
            print(f"- {item}")

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
