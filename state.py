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

    # UniProt / HumanBase look-ups
    uniprot_entries_base: Annotated[Dict[str, Dict], operator.or_]
    uniprot_entries_gwas: Annotated[Dict[str, Dict], operator.or_]
    gene_disease_traits: Annotated[Dict[str, List[str]], operator.or_]
    gene_function_traits: Annotated[Dict[str, List[str]], operator.or_]
    gene_GO_traits: Annotated[Dict[str, List[str]], operator.or_]
    humanbase_predictions: Annotated[Dict[str, List[Dict[str, Any]]], operator.or_]
    bioGRID_predictions: Annotated[Dict[str, List[str]], operator.or_]
    gwas_associations: Annotated[Dict[str, Dict[str, Any]], operator.or_]

    # LLM bookkeeping
    context_block: Annotated[str, operator.add]
    llm_json: Annotated[Dict[str, Any], operator.or_]
    verification: str
    verify_attempts: int
    used_tools: Annotated[List[str], operator.add]
