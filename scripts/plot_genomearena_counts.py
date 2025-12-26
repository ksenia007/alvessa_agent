"""
Plot a donut chart of question counts per subfolder under benchmarks_generation/questions/GenomeArena.
Requires pandas and matplotlib.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main() -> int:
    base = Path("/Users/sokolova/Documents/research/alvessa_agent/benchmarks_generation/questions/GenomeArena")
    if not base.exists():
        print(f"Base path not found: {base}", file=sys.stderr)
        return 1

    counts = {}
    try:
        for csv in base.rglob("*.csv"):
            try:
                n = len(pd.read_csv(csv))
            except Exception as exc:
                print(f"Failed on {csv}: {exc}", file=sys.stderr)
                continue
            folder = csv.parent.name
            counts[folder] = counts.get(folder, 0) + n
    except Exception as exc:
        print(f"Error scanning CSVs: {exc}", file=sys.stderr)
        return 1

    if not counts:
        print("No CSVs found or all failed to load.", file=sys.stderr)
        return 1

    labels = list(counts.keys())
    sizes = [counts[k] for k in labels]
    total = sum(sizes)

    # Donut plot
    fig, ax = plt.subplots(figsize=(6, 6), dpi=300)
    wedges, texts = ax.pie(sizes, labels=None, startangle=90)
    # Draw center circle
    centre_circle = plt.Circle((0, 0), 0.60, fc="white")
    fig.gca().add_artist(centre_circle)
    ax.axis("equal")
    # Legend with counts and percentages
    legend_labels = [
        f"{lbl}: {cnt} ({cnt/total:.1%})" for lbl, cnt in zip(labels, sizes)
    ]
    ax.legend(wedges, legend_labels, title="Subfolders", loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_title(f"GenomeArena question counts (total {total})")

    out_path = base / "genomearena_counts_donut.png"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)

    print("Counts per folder:")
    for lbl, cnt in sorted(counts.items()):
        print(f"  {lbl}: {cnt}")
    print(f"Total: {total}")
    print(f"Saved plot to: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
