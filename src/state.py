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
from src.alvessa.domain.gene_class import Gene
from src.alvessa.domain.variant_class import Variant
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
    
    # to replace the keys commented out below w/ proper gene objects
    gene_entities: Annotated[Dict[str, "Gene"], operator.or_]
    variant_entities: Annotated[Dict[str, "Variant"], operator.or_]

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

    # Protein structure and druggablility visualization
    prot_html: Annotated[str, operator.add]
    
    # General free-text annotations (not tied to a specific gene)
    text_notes: Annotated[List[str], operator.add]
    two_sample_mr_ui: Annotated[Dict[str, Any], operator.or_]
