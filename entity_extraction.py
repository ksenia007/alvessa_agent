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
from config import DEBUG, GENE_EXTRACT_MODEL, ENTITY_EXTRACTION_METHOD, GLINER_MODEL, GLINER_THRESHOLD, GLINER_ENTITY_LABELS
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


def _extract_genes_with_gliner(text: str) -> List[str]:
    """
    Extract gene symbols using GLiNER model.
    
    Parameters
    ----------
    text : str
        Input text to extract genes from
        
    Returns
    -------
    List[str]
        List of unique gene symbols found
    """
    model = _get_gliner_model()
    
    # Predict entities using GLiNER
    entities = model.predict_entities(text, GLINER_ENTITY_LABELS, threshold=GLINER_THRESHOLD)
    
    # Filter for gene entities and extract text
    genes = []
    for entity in entities:
        if entity["label"].lower() in ["gene", "protein"]:  # Include both genes and proteins
            gene_text = entity["text"].strip()
            # Basic validation: gene symbols are typically uppercase and 3-10 characters
            if re.match(r'^[A-Z0-9-]+$', gene_text) and 2 <= len(gene_text) <= 15:
                genes.append(gene_text)
    
    # Remove duplicates and return
    genes = list(set(genes))
    if DEBUG:
        print(f"[_extract_genes_with_gliner] Found genes: {genes}")
    
    return genes


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


def gene_extraction_node(state: "State") -> "State":
    """
    Identify HGNC gene symbols mentioned in the last user message.
    Uses either Claude or GLiNER based on ENTITY_EXTRACTION_METHOD configuration.

    Returns
    -------
    State
        Same state with a new `"genes"` list.
    """
    user_input: str = state["messages"][-1]["content"]
    
    # Choose extraction method based on configuration
    if ENTITY_EXTRACTION_METHOD.lower() == "gliner":
        genes = _extract_genes_with_gliner(user_input)
    else:
        genes = _extract_genes_with_claude(user_input)
    
    if DEBUG:
        print(f"[gene_extraction_node] extracted using {ENTITY_EXTRACTION_METHOD}: {genes}")
    
    return {"genes": genes}


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


def trait_extraction_node(state: "State") -> "State":
    """
    Extract disease/trait entities from the user message using GLiNER.
    This is used to identify specific diseases or traits for trait-based queries.

    Returns
    -------
    State
        Same state with a new `"traits"` list containing extracted diseases/traits.
    """
    user_input: str = state["messages"][-1]["content"]
    
    try:
        model = _get_gliner_model()
        
        # Predict entities using GLiNER with focus on diseases and traits
        disease_trait_labels = ["Disease", "Trait", "Phenotype", "Disorder", "Syndrome", "Condition"]
        entities = model.predict_entities(user_input, disease_trait_labels, threshold=GLINER_THRESHOLD)
        
        # Extract and clean disease/trait entities
        traits = []
        for entity in entities:
            trait_text = entity["text"].strip()
            # Basic validation and cleaning
            if len(trait_text) > 2 and trait_text not in traits:  # Avoid very short or duplicate entries
                traits.append(trait_text)
        
        # Remove duplicates while preserving order
        traits = list(dict.fromkeys(traits))
        
        if DEBUG:
            print(f"[trait_extraction_node] extracted traits/diseases: {traits}")
        
        return {"traits": traits}
        
    except Exception as e:
        if DEBUG:
            print(f"[trait_extraction_node] Error extracting traits: {e}")
        return {"traits": []}


def comprehensive_entity_extraction_node(state: "State") -> "State":
    """
    Extract all types of entities (genes, diseases, traits, etc.) using GLiNER.
    This provides a complete entity extraction covering genes, diseases, and other biomedical entities.

    Returns
    -------
    State
        Updated state with `"genes"`, `"traits"`, and `"all_entities"` fields.
    """
    user_input: str = state["messages"][-1]["content"]
    
    try:
        model = _get_gliner_model()
        
        # Predict all entity types using GLiNER
        entities = model.predict_entities(user_input, GLINER_ENTITY_LABELS, threshold=GLINER_THRESHOLD)
        
        # Organize entities by type
        genes = []
        traits = []
        all_entities = {}
        
        for entity in entities:
            label = entity["label"]
            text = entity["text"].strip()
            
            # Populate all_entities dict
            if label not in all_entities:
                all_entities[label] = []
            if text not in all_entities[label]:
                all_entities[label].append(text)
            
            # Extract genes
            if label.lower() in ["gene", "protein"]:
                # Gene validation
                if re.match(r'^[A-Z0-9-]+$', text) and 2 <= len(text) <= 15:
                    if text not in genes:
                        genes.append(text)
            
            # Extract diseases/traits
            elif label.lower() in ["disease", "trait", "phenotype", "disorder", "syndrome", "condition"]:
                if len(text) > 2 and text not in traits:
                    traits.append(text)
        
        # Remove duplicates while preserving order
        genes = list(dict.fromkeys(genes))
        traits = list(dict.fromkeys(traits))
        
        if DEBUG:
            print(f"[comprehensive_entity_extraction_node] extracted genes: {genes}")
            print(f"[comprehensive_entity_extraction_node] extracted traits: {traits}")
            print(f"[comprehensive_entity_extraction_node] all entities: {list(all_entities.keys())}")
        
        return {
            "genes": genes,
            "traits": traits,
            "all_entities": all_entities
        }
        
    except Exception as e:
        if DEBUG:
            print(f"[comprehensive_entity_extraction_node] Error extracting entities: {e}")
        return {
            "genes": [],
            "traits": [],
            "all_entities": {}
        }


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
