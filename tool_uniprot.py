"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-06-25
Updated: 2025-06-26


Description: 

UniProt tool: fetch entry + extract disease / function / GO traits.
"""

from __future__ import annotations
import requests
from typing import Dict, List, Optional
from state import State

DEBUG = True

def get_uniprot_entry_for_gene(gene_symbol: str) -> Optional[Dict]:
    """
    Retrieve the *reviewed* UniProtKB record for a human gene symbol.

    Returns
    -------
    dict | None
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene:{gene_symbol} AND organism_id:9606 AND reviewed:true",
        "format": "json",
        "size": 1,
    }

    if DEBUG:
        print(f"Fetching UniProt entry for gene: {gene_symbol}")

    try:
        resp = requests.get(base_url, params=params, timeout=12)
        resp.raise_for_status()
        results = resp.json().get("results", [])
        return results[0] if results else None
    except Exception as exc:
        print(f"Error fetching UniProt entry: {exc}")
        return None


def extract_disease_from_uniprot_entry(entry: Dict) -> List[str]:
    """Pull disease names from a UniProt entry."""
    if DEBUG:
        print("Extracting disease traits")
    traits: set[str] = set()
    for item in entry.get("comments", []):
        if item.get("commentType") == "DISEASE":
            disease = item.get("disease", {})
            name = disease.get("diseaseId") or disease.get("acronym") or disease.get(
                "name", {}
            ).get("value")
            if name:
                traits.add(name)
    return sorted(traits)


def extract_function_from_uniprot_entry(entry: Dict) -> List[str]:
    """Pull free-text *Function* comment blocks."""
    if DEBUG:
        print("Extracting function traits")
    traits: set[str] = set()
    for item in entry.get("comments", []):
        if item.get("commentType") == "FUNCTION":
            for text_block in item.get("texts", []):
                traits.add(text_block.get("value", ""))
    return sorted(traits)


def extract_GO_from_uniprot_entry(entry: Dict) -> List[str]:
    """Extract GO terms from cross-references."""
    if DEBUG:
        print("Extracting GO terms")
    traits: set[str] = set()
    for ref in entry.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "GO":
            for prop in ref.get("properties", []):
                if prop.get("key") == "GoTerm":
                    traits.add(prop["value"])
    return sorted(traits)



def uniprot_node(state: "State") -> "State":
    """Download UniProt entries for every gene symbol."""

    base_genes = list(set(state.get("genes", [])))
    gwas_linked_genes = list(set(state.get("gwas_linked_genes", [])))
    genes = list(set(base_genes + gwas_linked_genes))

    uniprot_entries_base: Dict[str, Dict] = {}
    uniprot_entries_gwas: Dict[str, Dict] = {}

    for gene in genes:
        entry_base = get_uniprot_entry_for_gene(gene)
        if entry_base:
            uniprot_entries_base[gene] = entry_base
        entry_gwas = get_uniprot_entry_for_gene(gene)
        if entry_gwas:
            uniprot_entries_gwas[gene] = entry_gwas

    return {**state, "uniprot_entries_base": uniprot_entries_base, "uniprot_entries_gwas": uniprot_entries_gwas}


def trait_disease_extraction_node(state: "State") -> "State":
    """Attach disease traits pulled from UniProt entries."""
    entries = state.get("uniprot_entries", {})
    gene_traits: Dict[str, List[str]] = {
        g: extract_disease_from_uniprot_entry(e) for g, e in entries.items()
    }
    # prune empties
    gene_traits = {g: t for g, t in gene_traits.items() if t}
    return {**state, "gene_disease_traits": gene_traits}


def trait_function_extraction_node(state: "State") -> "State":
    """Attach free-text functional annotations."""
    entries = state.get("uniprot_entries", {})
    gene_traits: Dict[str, List[str]] = {
        g: extract_function_from_uniprot_entry(e) for g, e in entries.items()
    }
    gene_traits = {g: t for g, t in gene_traits.items() if t}
    return {**state, "gene_function_traits": gene_traits}


def trait_GO_extraction_node(state: "State") -> "State":
    """Attach GO term lists."""
    entries = state.get("uniprot_entries", {})
    gene_traits: Dict[str, List[str]] = {
        g: extract_GO_from_uniprot_entry(e) for g, e in entries.items()
    }
    gene_traits = {g: t for g, t in gene_traits.items() if t}
    return {**state, "gene_GO_traits": gene_traits}


def has_uniprot_entries(state: "State") -> bool:
    """Whether any UniProt entries have been fetched (LangGraph edge helper)."""
    return bool(state.get("uniprot_entries"))
