"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-06-25
Updated: 2025-08-11


Description: 

Nodes to extract gene symbols from user input."""

from __future__ import annotations
from typing import List, Dict, Any
import re

from claude_client import claude_call
from config import DEBUG, GENE_EXTRACT_MODEL, GLINER_MODEL, GLINER_THRESHOLD, GLINER_ENTITY_LABELS
from state import State

# Global variable to cache the GLiNER model
_gliner_model = None


def _get_gliner_model():
    """
    Get or initialize the GLiNER model (cached for efficiency).
    """
    global _gliner_model
    if _gliner_model is None:
        try:
            from gliner import GLiNER
            _gliner_model = GLiNER.from_pretrained(GLINER_MODEL)
            if DEBUG:
                print(f"[_get_gliner_model] Loaded GLiNER model: {GLINER_MODEL}")
        except ImportError:
            raise ImportError("GLiNER is not installed. Please install it with: pip install gliner")
        except Exception as e:
            raise RuntimeError(f"Failed to load GLiNER model {GLINER_MODEL}: {e}")
    return _gliner_model


def _extract_genes_with_claude(text: str) -> List[str]:
    """
    Extract gene symbols using Claude model.
    
    Parameters
    ----------
    text : str
        Input text to extract genes from
        
    Returns
    -------
    List[str]
        List of unique gene symbols found
    """
    system_message: str = (
        "Extract gene symbols from the message. Extract gene names **only if they appear verbatim in the input**. "
        "Reply with a comma-separated list of gene names only, no extra words. If no gene names are found, reply with an empty string."
    )
    response = claude_call(
        model=GENE_EXTRACT_MODEL,
        max_tokens=50,
        temperature=0,
        system=system_message,
        messages=[{"role": "user", "content": text}],
    )
    try:
        genes = [g.strip() for g in response.content[0].text.split(",") if g.strip()]
    except:
        print("[_extract_genes_with_claude] Error parsing genes")
        genes = []
    
    genes = list(set(genes))  # Ensure unique genes
    if DEBUG:
        print(f"[_extract_genes_with_claude] Found genes: {genes}")
    
    return genes


def _extract_entities_merged(text: str) -> Dict[str, Any]:
    """
    Extract entities using both Claude and GLiNER models, then merge results.
    Returns both genes specifically and all entity types from GLiNER.
    
    Parameters
    ----------
    text : str
        Input text to extract entities from
        
    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - "genes": List of unique gene symbols found by both methods combined
        - "all_entities": Dict of all entity types found by GLiNER
        - "traits": List of disease/trait entities from GLiNER
    """
    # Extract genes using Claude
    claude_genes = _extract_genes_with_claude(text)
    
    # Extract all entities using GLiNER
    model = _get_gliner_model()
    gliner_entities = model.predict_entities(text, GLINER_ENTITY_LABELS, threshold=GLINER_THRESHOLD)
    # import ipdb; ipdb.set_trace()
    
    # Organize GLiNER entities by type
    all_entities = {}
    gliner_genes = []
    traits = []
    
    for entity in gliner_entities:
        label = entity["label"]
        text_entity = entity["text"].strip()
        
        # Populate all_entities dict
        if label not in all_entities:
            all_entities[label] = []
        if text_entity not in all_entities[label]:
            all_entities[label].append(text_entity)
        
        # Extract genes from GLiNER results
        if label.lower() in ["gene", "protein"]:
            # Gene validation
            if re.match(r'^[A-Z0-9-]+$', text_entity) and 2 <= len(text_entity) <= 15:
                if text_entity not in gliner_genes:
                    gliner_genes.append(text_entity)
        
        # Extract diseases/traits
        elif label.lower() in ["disease", "trait", "phenotype", "disorder", "syndrome", "condition"]:
            if len(text_entity) > 2 and text_entity not in traits:
                traits.append(text_entity)
    
    # Merge gene results (Claude + GLiNER)
    merged_genes = claude_genes + [gene for gene in gliner_genes if gene not in claude_genes]
    
    # Remove duplicates while preserving order
    merged_genes = list(dict.fromkeys(merged_genes))
    traits = list(dict.fromkeys(traits))
    
    if DEBUG:
        print(f"[_extract_entities_merged] Claude genes: {claude_genes}")
        print(f"[_extract_entities_merged] GLiNER genes: {gliner_genes}")
        print(f"[_extract_entities_merged] Merged genes: {merged_genes}")
        print(f"[_extract_entities_merged] GLiNER traits: {traits}")
        print(f"[_extract_entities_merged] All entity types: {list(all_entities.keys())}")
    
    return {
        "genes": merged_genes,
        "all_entities": all_entities,
        "traits": traits
    }


def entity_extraction_node(state: "State") -> "State":
    """
    Extract entities using both Claude and GLiNER methods for comprehensive query understanding.
    Provides extraction of genes, traits, and all biomedical entity types mentioned in user queries.
    Uses Claude for gene extraction and GLiNER for comprehensive biomedical entity extraction.

    Returns
    -------
    State
        Updated state with:
        - "genes": List of gene symbols from both Claude and GLiNER
        - "all_entities": Dict of all entity types found by GLiNER  
        - "traits": List of disease/trait entities from GLiNER
    """
    user_input: str = state["messages"][-1]["content"]
    
    # Use comprehensive merged extraction method
    extraction_result = _extract_entities_merged(user_input)
    
    if DEBUG:
        print(f"[entity_extraction_node] extracted genes: {extraction_result['genes']}")
        print(f"[entity_extraction_node] extracted traits: {extraction_result['traits']}")
        print(f"[entity_extraction_node] all entity types: {list(extraction_result['all_entities'].keys())}")
    
    return {
        "genes": extraction_result["genes"],
        "all_entities": extraction_result["all_entities"], 
        "traits": extraction_result["traits"]
    }


def gliner_entity_extraction_node(state: "State") -> "State":
    """
    Extract all types of entities using GLiNER model.
    This provides broader entity extraction beyond just genes.

    Returns
    -------
    State
        Same state with a new `"entities"` dict containing all extracted entities.
    """
    user_input: str = state["messages"][-1]["content"]
    model = _get_gliner_model()
    
    # Predict entities using GLiNER
    entities = model.predict_entities(user_input, GLINER_ENTITY_LABELS, threshold=GLINER_THRESHOLD)
    
    # Organize entities by type
    entities_by_type = {}
    for entity in entities:
        label = entity["label"]
        text = entity["text"].strip()
        
        if label not in entities_by_type:
            entities_by_type[label] = []
        
        # Avoid duplicates
        if text not in entities_by_type[label]:
            entities_by_type[label].append(text)
    
    if DEBUG:
        print(f"[gliner_entity_extraction_node] extracted entities: {entities_by_type}")
    
    return {"entities": entities_by_type}


def has_genes(state: "State") -> bool:
    """
    Edge-condition helper for LangGraph: returns `True` if any genes were found.
    """
    if DEBUG:
        print("[has_genes] genes present?", bool(state.get("genes")))
    return bool(state.get("genes"))


def has_traits(state: "State") -> bool:
    """
    Edge-condition helper for LangGraph: returns `True` if any traits/diseases were found.
    """
    traits_found = bool(state.get("traits"))
    if DEBUG:
        print(f"[has_traits] traits present? {traits_found}")
    return traits_found
