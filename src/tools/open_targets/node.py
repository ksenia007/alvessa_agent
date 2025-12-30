"""
Description:

Tool to query OpenTarget annotations with a list of gene symbols"""

from __future__ import annotations

from src.state import State
from src.tools.base import Node
from pathlib import Path
from typing import List
import pickle

DEBUG = True

REPO_ROOT = Path(__file__).resolve().parents[3]
LOCAL_DBS = REPO_ROOT / "local_dbs"      
OPEN_TARGETS = LOCAL_DBS / "open_targets"
TARGET_DISEASE_DATA = OPEN_TARGETS / "final_association_overall_direct/target_disease_direct_associations.pkl"
EXPRESSION_DATA = OPEN_TARGETS / "final_expression/expression.pkl"
ESSENTIALITY_DATA = OPEN_TARGETS / "final_target_essentiality/target_essentiality.pkl"
CONSTRAINT_DATA = OPEN_TARGETS / "final_constraint/constraint.pkl"
PHARMACOVIGILANCE_DATA = OPEN_TARGETS / "final_target_pharmacovigilance/target_adrs.pkl"
PHARMACOGENOMICS_DATA = OPEN_TARGETS /  "final_variant_pharmacogenomics/variant_pharmacogenomics.pkl"

CONSTRAINT_TYPES = {
    'syn': "Synonymous (silent) variants - don't change amino acid sequence of a protein",
    'mis': "Missense variants - lead to amino acid substitution in the protein",
    'lof': "Loss-of-function variants - include nonsense, frameshift, and canonical splice site variants"
}

def opentargets_agent(state: "State"):
    """
    LangGraph node that annotates each gene with Open Target data.

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state with the open target fields filled.

    """
    gene_entities = state.get("gene_entities") or {}
    variant_entities = state.get("variant_entities") or {}

    # Early return if nothing to process
    if not gene_entities and not variant_entities:
        if DEBUG:
            print("[Open Targets] No entities to process, skipping")
        return

    # Load gene-related data only if we have genes
    target_disease_dict = {}
    expression_dict = {}
    essentiality_dict = {}
    constraint_dict = {}
    pharmocovigilance_dict = {}

    if gene_entities:
        try:
            if TARGET_DISEASE_DATA.exists():
                with open(TARGET_DISEASE_DATA, 'rb') as file:
                    target_disease_dict = pickle.load(file)
        except Exception as e:
            print(f"Failed to load target disease data: {e}")

        try:
            if EXPRESSION_DATA.exists():
                with open(EXPRESSION_DATA, 'rb') as file:
                    expression_dict = pickle.load(file)
        except Exception as e:
            print(f"Failed to load expression data: {e}")

        try:
            if ESSENTIALITY_DATA.exists():
                with open(ESSENTIALITY_DATA, 'rb') as file:
                    essentiality_dict = pickle.load(file)
        except Exception as e:
            print(f"Failed to load essentiality data: {e}")

        try:
            if CONSTRAINT_DATA.exists():
                with open(CONSTRAINT_DATA, 'rb') as file:
                    constraint_dict = pickle.load(file)
        except Exception as e:
            print(f"Failed to load constraint data: {e}")

        try:
            if PHARMACOVIGILANCE_DATA.exists():
                with open(PHARMACOVIGILANCE_DATA, 'rb') as file:
                    pharmocovigilance_dict = pickle.load(file)
        except Exception as e:
            print(f"Failed to load pharmacovigilance data: {e}")

    # Load variant-related data only if we have variants
    pharmacogenomics_dict = {}
    if variant_entities:
        try:
            if PHARMACOGENOMICS_DATA.exists():
                with open(PHARMACOGENOMICS_DATA, 'rb') as file:
                    pharmacogenomics_dict = pickle.load(file)
        except Exception as e:
            print(f"Failed to load pharmacogenomics data: {e}")

    # Process genes
    for gene_name, gene in gene_entities.items():
        if not gene_name:
            continue

        summary_lines: List[str] = []

        try:
            # Disease annotations
            associated_disease_list = target_disease_dict[gene.symbol]
            gene.add_many_direct_disease_associations(associated_disease_list)

            summary_lines.append(f"*OpenTargets: Direct disease associations for {gene.symbol}: " 
                        + ', '.join(associated_disease_list) + ".")

            # print(len(associated_disease_list))
        
        except Exception as e:
            print(f"Failed to pull disease annotations: {gene}, {e}")

        try:
            # Tissue-specific expression
            tissue_zscore = expression_dict[gene.symbol]
                        
            gene.add_many_tissues_expression(tissue_zscore)

            summary_lines.append(
                "*OpenTargets: Tissue-specific expression binned z-scores for "
                f"{gene.symbol} (>=2 indicates 75th percentile and is tissue-specific; "
                "z-scores per gene per tissue, binned into 10 quantile-based bins): "
                + str(tissue_zscore) + "."
            )
            
        except Exception as e:
            print(f"Failed to pull tissue expression: {gene}, {e}")
        
        try:
            # DepMap Essentiality
            gene_is_essential = essentiality_dict[gene.symbol]
            gene.add_essentiality(gene_is_essential)

            if gene_is_essential:
                summary_lines.append(
                    f"*OpenTargets: {gene.symbol} is a core essential gene (unlikely to tolerate inhibition; "
                    "crucial for basic cellular function across many tissues)."
                )
            else:
                summary_lines.append(
                    f"*OpenTargets: {gene.symbol} is not a core essential gene (inhibition/knockout not consistently lethal "
                    "across most cell lines; not universally vital for survival)."
                )
        
        except Exception as e:
            print(f"Failed to pull essentiality: {gene}, {e}")

        try:       
            # Genetic Constraint
            gene_all_constraint = constraint_dict[gene.symbol]
            
            for constraint_type_name, score in gene_all_constraint.items():
                gene.add_constraint(score, constraint_type_name)
                summary_lines.append(
                    f"*OpenTargets: Genetic constraint score for {constraint_type_name} variants "
                    f"({CONSTRAINT_TYPES[constraint_type_name]}) in {gene.symbol}; "
                    "scores -1 to 1 map to LOEUF rank (lower = less tolerant): "
                    + str(score)
                )

        except Exception as e:
            print(f"Failed to pull constraint: {gene}, {e}")

        try:       
            # Pharmacovigilance
            gene_adrs = pharmocovigilance_dict[gene.symbol]
            gene.add_many_adverse_reactions(gene_adrs)
            summary_lines.append(
                f"*OpenTargets: Significant adverse drug reactions for drugs targeting {gene.symbol}: "
                + ', '.join(gene_adrs) + "."
            )
        
        except Exception as e:
            print(f"Failed to pull pharmacovigilance: {gene}, {e}")

        if summary_lines:
            gene.update_text_summaries(" ".join(summary_lines))
            gene.add_tool("OpenTargets")

    # Process variants - ONLY loop over variants that exist in the data
    if pharmacogenomics_dict and variant_entities:
        # Find intersection: only variants we have AND have data for
        variant_rsids_with_data = set(pharmacogenomics_dict.keys())
        variants_to_process = {
            variant_name: variant
            for variant_name, variant in variant_entities.items()
            if variant_name and hasattr(variant, 'rsID') and variant.rsID in variant_rsids_with_data
        }

        if DEBUG and variants_to_process:
            print(f"[Open Targets] Processing {len(variants_to_process)}/{len(variant_entities)} variants with pharmacogenomics data")

        for variant_name, variant in variants_to_process.items():
            summary_lines: List[str] = []

            try:
                # Pharmacogenomics - we KNOW this variant is in the dict
                variant_effects = pharmacogenomics_dict[variant.rsID]
                variant.add_many_drug_response_effects(variant_effects)
                summary_lines.append(
                    f"*OpenTargets: Drug response annotations for {variant.rsID}: "
                    + ', '.join(variant_effects) + "."
                )

            except Exception as e:
                print(f"Failed to pull pharmacogenomics: {variant}, {e}")

            if summary_lines:
                variant.update_text_summaries(" ".join(summary_lines))
                variant.add_tool("OpenTargets")

    if DEBUG:
        print(f"[Open Targets] Predictions fetched")

    return 

NODES: tuple[Node, ...] = (
    Node(
        name="OpenTargets",
        entry_point=opentargets_agent,
        description="Fetches target-disease associations, tissue-specific expression, essentiality, genetic constraint for synonymous, missense, and loss-of-function variants, associations between target genes and drug-induced adverse effects, and variant to drug response associations from Open Targets for the input genes.",
    ),
)
