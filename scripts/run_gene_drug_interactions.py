"""CLI helper to inspect gene-centric CIGS drug perturbations."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.tools.gene_drug_interactions.node import gene_drug_interactions_node


def build_state(genes: List[str]) -> Dict[str, Any]:
    return {
        "messages": [{"role": "user", "content": f"Querying genes: {', '.join(genes)}"}],
        "genes": genes,
        "gene_entities": {},
        "text_notes": [],
        "used_tools": [],
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Look up drugs that perturb the specified genes in CIGS (padj<0.05)."
    )
    parser.add_argument("genes", nargs="+", help="Gene symbols (e.g., TP53 BRCA1)")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output")
    args = parser.parse_args()

    state = build_state(args.genes)
    result = gene_drug_interactions_node(state)
    payload = result.get("gene_drug_interactions")

    if not payload:
        notes = result.get("text_notes") or []
        print("No gene-centric drug interactions were found.")
        if notes:
            print("Notes:")
            for note in notes:
                print(f"  - {note}")
        return

    subset = {gene: payload.get(gene) for gene in (gene.upper() for gene in args.genes)}
    indent = 2 if args.pretty else None
    print(json.dumps(subset, indent=indent, sort_keys=True))


if __name__ == "__main__":
    main()
