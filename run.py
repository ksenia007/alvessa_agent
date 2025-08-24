""" 
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-06-25
Updated: 2025-08-12


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
import json



def run_pipeline(user_message: str, prompt: str = '', mc_setup: bool = False) -> Dict:
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
    graph = build_graph(mc_setup = mc_setup)
    state = graph.invoke({"messages": [{"role": "user", "content": user_message}], 
                          "prompt": prompt})
    return state


if __name__ == "__main__":
    logfile = open("demo.log", "w")
    sys.stdout = logfile

    EXAMPLE_QUESTIONS: List[str] = [
        # "What microRNAs regulate the cancer-related functions of TP53, and how are these connected to pathways and protein interactions?"
        #"Do not run any of the tools, just proceed directly to the LLM.",
        #"Tell me about exons in TP53",
        #"Which diseases and traits are associated with the genes TP53 and KRAS?",
        # "Which gene is the best drug target for virally induced cancers, KRAS or TP53?",
        #"How does NUCKS1 play a role in cancer and in viral infections, and what is the overlap of these roles?",
        # "Why is TP53 important for all cancers but BRCA1 only in breast and ovarian cancers?",
        # "Through what pathways or protein interactions does viral gene E6 (or E7) modulate cellular metabolism?"
        # "Is there evidence of protein-protein interaction or functional overlaps between TP53 and KRAS?"
        # "Are there common patterns in regulatory activity of variants found for TP53 and KRAS in GWAS studies?"
        """Which of the following genes is a computationally predicted human gene target of the miRNA MIR3140_5P according to miRDB v6.0? 

A. TMCC3
B. EDEM3
C. ATP6V1B1
D. GBP2"""
    ]

    with open("demo.txt", "w") as f:
        for q in EXAMPLE_QUESTIONS:
            print("\n" + "=" * 80)
            print("Q:", q)
            result = run_pipeline(q)
            # save for the html
            with open("demo.json", "w") as jf:
                json.dump(result, jf, indent=2, default=str)
                
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
            

            with open("demo.json", "w") as jf:
                json.dump(result, jf, indent=2, default=str)


