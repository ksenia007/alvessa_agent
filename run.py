""" 
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-06-25
Updated: 2025-06-26


Description: 

Entry point for LangGraph-based gene function reasoning demo.
At runtime it goes through the example questions listed below. 

Outputs:
Saves clean user-friendly output to demo.txt, log (if in debug mode in config) as demo.log,
and also a PNG visual of the graph as graph_diagram.png

"""

from __future__ import annotations
import pprint
from typing import Dict, List

from graph_builder import build_graph

import sys

logfile = open("demo.log", "w")
sys.stdout = logfile


def run_pipeline(user_message: str) -> Dict:
    """
    Execute the LangGraph workflow on a single user prompt.

    Parameters
    ----------
    user_message
        The natural-language question.

    Returns
    -------
    dict
        Final LangGraph state for inspection.
    """
    graph = build_graph()
    state = graph.invoke({"messages": [{"role": "user", "content": user_message}]})
    return state


if __name__ == "__main__":
    EXAMPLE_QUESTIONS: List[str] = [
        # "Which diseases and traits are associated with the genes TP53 and KRAS?",
        # "Which gene is the best drug target for virally induced cancers, KRAS or TP53?",
        # "Describe distinct roles of TP53 and BRCA1 in cancer biology.",
        # "How does NUCKS1 play a role in cancer and in viral infections, and what is the overlap of these roles?",
        # "Why is TP53 important for all cancers but BRCA1 only in breast and ovarian cancers?",
        # "Through what pathways or protein interactions does viral gene E6 (or E7) modulate cellular metabolism?"
        # "Is there evidence of protein-protein interaction or functional overlaps between TP53 and KRAS?"
        "How many variants are identified for TP53 and KRAS? What does it imply about the genes in terms of regulatory activity??",
    ]

    with open("demo.txt", "w") as f:
        for q in EXAMPLE_QUESTIONS:
            print("\n" + "=" * 80)
            print("Q:", q)
            result = run_pipeline(q)
            last_msg = result["messages"][-1]["content"]
            answer = result["llm_json"].get("answer", "")
            evidence = result["llm_json"].get("evidence", [])

            # Print to console
            pprint.pp(last_msg)
            pprint.pp(result["llm_json"])

            # Save to file
            f.write("=" * 80 + "\n")
            f.write(f"Q: {q}\n\n")
            f.write(f"Answer:\n{answer}\n\n")
            f.write("Evidence:\n")
            for ev in evidence:
                f.write(f"  - {ev}\n")
            f.write("\n\n")

