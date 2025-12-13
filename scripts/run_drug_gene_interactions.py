"""CLI helper to exercise the drug_gene_interactions node against local CIGS data."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.tools.drug_gene_interactions.node import drug_gene_interactions_node


def build_state(drugs: List[str]) -> Dict[str, Any]:
    """Construct the minimal LangGraph state required by the node."""
    return {
        "messages": [{"role": "user", "content": f"Querying drugs: {', '.join(drugs)}"}],
        "drugs": drugs,
        "drug_entities": {},
        "text_notes": [],
        "used_tools": [],
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Look up CIGS drug-gene interactions for the provided drug names."
    )
    parser.add_argument(
        "drugs",
        nargs="+",
        help="Drug names or MedChemExpress catalog numbers (e.g., 'SKF-96365 (hydrochloride)' or 'HY-100001').",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print JSON output for easier reading.",
    )
    args = parser.parse_args()

    state = build_state(args.drugs)
    result = drug_gene_interactions_node(state)
    payload = result.get("drug_gene_interactions")

    if not payload:
        notes = result.get("text_notes") or []
        print("No significant genes found or drugs could not be resolved.")
        if notes:
            print("Notes:")
            for note in notes:
                print(f"  - {note}")
        return

    indent = 2 if args.pretty else None
    print(json.dumps(payload, indent=indent, sort_keys=True))


if __name__ == "__main__":
    main()
