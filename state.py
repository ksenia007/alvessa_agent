"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-06-25
Updated: 2025-06-26


Description: 

TypedDict schema describing the mutable LangGraph state."""

from __future__ import annotations
import operator
from typing import Annotated, Any, Dict, List

from typing_extensions import TypedDict


class State(TypedDict, total=False):
    # Conversation
    messages: Annotated[List[Dict[str, str]], operator.add]
    genes: Annotated[List[str], operator.add]
    traits: Annotated[List[str], operator.add]
    proteins: Annotated[List[str], operator.add]
    transcripts: Annotated[List[str], operator.add]
    variants: Annotated[Dict[str, Dict[str, Dict[str, Any]]], operator.or_]
    chr_pos_variants: Annotated[Dict[str, Dict[str, Dict[str, Any]]], operator.or_]
    gene_level_gencode: Annotated[Dict[str, Dict[str, Any]], operator.or_]
    prompt: Annotated[str, operator.add]
    mc_setup: Annotated[bool, operator.and_]

    # UniProt / HumanBase look-ups
    uniprot_entries_base: Annotated[Dict[str, Dict], operator.or_]
    uniprot_entries_gwas: Annotated[Dict[str, Dict], operator.or_]
    gene_disease_traits: Annotated[Dict[str, List[str]], operator.or_]
    gene_function_traits: Annotated[Dict[str, List[str]], operator.or_]
    gene_GO_traits: Annotated[Dict[str, List[str]], operator.or_]
    humanbase_predictions: Annotated[Dict[str, List[Dict[str, Any]]], operator.or_]
    humanbase_expecto: Annotated[Dict[str, List[Dict[str, Any]]], operator.or_]
    tissue_expression_preds_variant_text_description: Annotated[Dict[str, Any], operator.or_]
    expression_preds_variant_table: Annotated[Dict[str, Dict[str, Any]], operator.or_]
    biogrid_predictions: Annotated[Dict[str, List[str]], operator.or_]
    biogrid_interaction_groups: Annotated[Dict[str, List[str]], operator.or_]
    biogrid_interactions_select_nonhuman: Annotated[Dict[str, List[str]], operator.or_]
    reactome_pathways: Annotated[Dict[str, List[str]], operator.or_]
    biogrid_summarized_go: Annotated[Dict[str, List[str]], operator.or_]
    dbsnp_variants: Annotated[Dict[str, Dict[str, Dict[str, Any]]], operator.or_]
    dbsnp_summaries: Annotated[Dict[str, Dict[str, Any]], operator.or_]
    sei_predictions: Annotated[Dict[str, Dict[str, Any]], operator.or_]
    alphamissense_predictions: Annotated[Dict[str, Dict[str, Any]], operator.or_]
    mirDB_targets: Annotated[Dict[str, List[str]], operator.or_]
    
    
    # GWAS associations
    gwas_associations: Annotated[Dict[str, Dict[str, Any]], operator.or_]

    # LLM bookkeeping
    context_block: Annotated[str, operator.add]
    llm_json: Annotated[Dict[str, Any], operator.or_]
    verification: str
    verify_attempts: int
    tool_updates: int
    used_tools: Annotated[List[str], operator.add]
    use_tools: Annotated[List[str], operator.add]  # tools to use in the current run
    
    # interactive view
    ui: Annotated[Dict[str, Any], operator.or_]  # e.g., {"panels": [...]}
