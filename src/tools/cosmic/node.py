
from __future__ import annotations

import time

from src.state import State
from src.tools.base import Node
from typing import List
from pathlib import Path
import pickle

DEBUG = True

REPO_ROOT = Path(__file__).resolve().parents[3]
LOCAL_DBS = REPO_ROOT / "local_dbs"   
HALLMARK_DATA_PATH = LOCAL_DBS / "cosmic_gene_hallmark_dict.pkl"
GENE_CENSUS_DATA_PATH = LOCAL_DBS / "cosmic_gene_census_data_dict.pkl"


def cosmic_agent(state: "State") -> "State":
    gene_entities = state.get("gene_entities") or {}

    for gene_name, gene in gene_entities.items():
        if not gene_name:
            continue

        summary_lines: List[str] = []

        try:

            with open(HALLMARK_DATA_PATH, "rb") as f:  
                hallmark_data = pickle.load(f)

            with open(GENE_CENSUS_DATA_PATH, "rb") as f:  
                gene_census_data = pickle.load(f)

            if gene_name in gene_census_data:
                associated_data = gene_census_data[gene_name]

                for row in associated_data:
                    tumour_types_somatic = row['tumour_types_somatic']
                    tumour_types_germline = row['tumour_types_germline']
                    cancer_syndrome = row['cancer_syndrome']
                    tissue_type = row['tissue_type']
                    role_in_cancer = row['role_in_cancer']
                    tier = row['tier']

                    summary = ""
                    if tumour_types_somatic:
                        summary += f"Somatic mutations in {gene.symbol} are associated with the following diseases: {tumour_types_somatic}. "
                    if tumour_types_germline:
                        summary += f"Germline mutations in {gene.symbol} are associated with the following diseases: {tumour_types_germline}. "
                    if cancer_syndrome:
                        summary += f"Cancer syndrome associated with germline mutation: {cancer_syndrome}. "
                    if tissue_type:
                        summary += f"Types of tissue: {tissue_type}. "
                    if role_in_cancer:
                        summary += f"Role in cancer (oncogene: hyperactivity of the gene drives the transformation; TSG: loss of gene function drives the transformation. Some genes can play either of these roles depending on cancer type. Fusion: the gene is known to be involved in oncogenic fusions) associated with {gene.symbol}: {role_in_cancer}. "
                    if tier:
                        summary += f"Which tier of the Cancer Gene Census {gene.symbol} belongs to (Tier 1 indicates the highest confidence and strongest evidence of driving cancer while Tier 2 indicates strong but less extensive evidence): Tier {tier}. "                                                

                    summary_lines.append(f"*Cosmic: {summary}")

            if gene_name in hallmark_data:
                associated_hallmarks = hallmark_data[gene_name]

                for h in associated_hallmarks:
                    hallmark_name = h['hallmark']
                    cell_type = h['cell_type']
                    impact = h['impact']
                    description = h['description']

                    summary = f"Hallmark (Name of the biological process that when dysregulated, may promote cancer or other data category describing the role of a gene in cancer) associated with {gene.symbol}: {hallmark_name}. "

                    if cell_type:
                        summary += f"Cell type of hallmark (tissue or cancer for which it is described): {cell_type}. "
                    if impact:
                        summary += f"Impact of hallmark (describes how the gene activity impacts the hallmarks of cancer): {impact}. "
                    if description:
                        summary += f"Description of hallmark (describes how the gene activity impacts the hallmarks of cancer): {description}. "

                    summary_lines.append(f"*Cosmic: {summary}")

            print(summary_lines)

            if summary_lines:
                gene.update_text_summaries(" ".join(summary_lines))

                gene.add_tool("Cosmic")

        except Exception as e:
            print(gene_name, e)


        time.sleep(0.3)  # courteous pause

    if DEBUG:
        print(f"[Cosmic] Predictions fetched")

    return 

NODES: tuple[Node, ...] = (
    Node(
        name="Cosmic",
        entry_point=cosmic_agent,
        description="Fetches information about genes linked to cancer through the COSMIC Cancer Gene Census, including diseases associated with somatic/germline mutations, role in cancer, biological mechanisms (hallmarks) that underlie cancer, and evidence tiers.",
    ),
)
