from __future__ import annotations
from typing import List, Dict, Any, Optional
import re

from claude_client import claude_call
from config import DEBUG, GENE_EXTRACT_MODEL, GLINER_MODEL, GLINER_THRESHOLD, GLINER_ENTITY_LABELS
from state import State
from flair.data import Sentence

# Global variables to cache the models
_gliner_model = None
_flair_model = None


# --- Model Loading Helpers ---

def _get_gliner_model():
    """Initializes and caches the GLiNER model."""
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


def _get_flair_model():
    """Initializes and caches the Flair model."""
    global _flair_model
    if _flair_model is None:
        try:
            from flair.nn import Classifier
            _flair_model = Classifier.load("hunflair2")
            if DEBUG:
                print(f"[_get_flair_model] Loaded FLAIR model: hunflair2")
        except ImportError:
            raise ImportError("Flair is not installed. Please install it with: pip install flair")
        except Exception as e:
            raise RuntimeError(f"Failed to load FLAIR model: {e}")
    return _flair_model

# --- New Variant Extraction Helpers ---

def _extract_rsids(text: str) -> List[str]:
    """Extracts genetic variants in the rs<id> format (e.g., rs12345)."""
    rs_id_pattern = r"rs\d+"
    found_rs_ids = re.findall(rs_id_pattern, text, re.IGNORECASE)
    # Normalize to lowercase for consistency
    found_rs_ids = [rsid.lower() for rsid in found_rs_ids]
    if DEBUG:
        print(f"[_extract_rsids] Found: {found_rs_ids}")
    return found_rs_ids

def _extract_chr_allele_pos(text: str) -> List[str]:
    """Extracts variants in chr<num>:<pos>:<ref>><alt> format (e.g., chr7:55249071:C>T)."""
    # This regex captures the specific allele change format.
    # It looks for one of A, T, C, or G for the reference and alternate alleles.
    chr_pattern = r"chr\w+:\d+:[ATCG]>[ATCG]"
    # findall returns the full match for patterns without capture groups
    matches = re.findall(chr_pattern, text, re.IGNORECASE)
    if DEBUG:
        print(f"[_extract_chr_allele_pos] Found: {matches}")
    return matches

def get_variant_coordinates(variant_string: str) -> Optional[Dict[str, Any]]:
    """
    Parses a variant string and returns its components in a structured dictionary.
    
    Args:
        variant_string: A string in chr<num>:<pos>:<ref>><alt> format.
        
    Returns:
        A dictionary with variant coordinates or None if the format is invalid.
    """
    # Regex to capture the components of the variant string
    pattern = r"chr(\w+):(\d+):([ATCG])>([ATCG])"
    match = re.match(pattern, variant_string, re.IGNORECASE)
    
    if not match:
        if DEBUG:
            print(f"[get_variant_coordinates] Invalid format for: {variant_string}")
        return None
        
    # Extract components from the matched groups
    chrom, pos, ref, alt = match.groups()
    
    # Construct the result dictionary
    coordinates = {
        "chrom": chrom,
        "pos": int(pos),
        "ref": ref.upper(),
        "alt": alt.upper(),
        "assembly": "GRCh38.p14"  # Hard-coded as requested
    }
    return coordinates


def _extract_variants_with_regex(text: str) -> Dict[str, List[str]]:
    """
    Extracts all supported genetic variant formats using modular regex functions.
    """
    # Call each specific extraction function
    rsids = _extract_rsids(text)
    allele_pos_variants = _extract_chr_allele_pos(text)

    # Combine all found variants into a single list and remove duplicates
    variants = list(set(rsids))
    chr_pos_variants = list(set(allele_pos_variants))
    
    result = {
        "variants": {variant: {"rsid": variant} for variant in variants}, 
        "chr_pos_variants": {
            chr_pos_variant: {
                "coordinates": get_variant_coordinates(chr_pos_variant) 
            }
            for chr_pos_variant in chr_pos_variants
        }
    }
    if DEBUG:
        print(f"[_extract_variants_with_regex] Combined variants: {result}")
        
    return result


# --- Core Extraction Logic ---

def _extract_entities_with_claude(text: str) -> Dict[str, List[str]]:
    """Extracts gene symbols using Claude and returns them in a standard format."""
    system_message = (
        "Extract gene symbols from the message. Extract gene names **only if they appear verbatim in the input**. "
        "Reply with a comma-separated list of gene names only, no extra words. If no gene names are found, reply with an empty string."
    )
    try:
        response = claude_call(
            model=GENE_EXTRACT_MODEL,
            max_tokens=50,
            temperature=0,
            system=system_message,
            messages=[{"role": "user", "content": text}],
        )
        genes = [g.strip() for g in response.content[0].text.split(",") if g.strip()]
    except Exception as e:
        print(f"[_extract_entities_with_claude] Error parsing genes: {e}")
        genes = []

    result = {"genes": list(set(genes)), "traits": []}
    if DEBUG:
        print(f"[_extract_entities_with_claude] Found: {result}")
    return result


def _extract_entities_with_flair(text: str) -> Dict[str, List[str]]:
    """Extracts gene and trait entities using Flair."""
    flair_model = _get_flair_model()
    sentence = Sentence(text)
    flair_model.predict(sentence)

    genes, traits = [], []
    for label in sentence.get_labels():
        entity_text = label.data_point.text.strip()
        entity_type = label.value.lower()

        if entity_type in ["gene", "protein"]:
            genes.append(entity_text)
        elif entity_type in ["disease", "trait", "phenotype", "disorder", "syndrome", "condition"]:
            traits.append(entity_text)
            
    result = {"genes": list(set(genes)), "traits": list(set(traits))}
    if DEBUG:
        print(f"[_extract_entities_with_flair] Found: {result}")
    return result


def _extract_entities_with_gliner(text: str) -> Dict[str, Any]:
    """Extracts all entities using GLiNER, returning genes and traits."""
    gliner_model = _get_gliner_model()
    entities = gliner_model.predict_entities(text, GLINER_ENTITY_LABELS, threshold=GLINER_THRESHOLD)
    
    genes, traits = [], []
    for entity in entities:
        label = entity["label"]
        text_entity = entity["text"].strip()

        if label.lower() in ["gene", "protein"]:
            genes.append(text_entity)
        elif label.lower() in ["disease", "trait", "phenotype", "disorder", "syndrome", "condition"]:
            traits.append(text_entity)

    result = {
        "genes": list(set(genes)), 
        "traits": list(set(traits))
    }
    if DEBUG:
        print(f"[_extract_entities_with_gliner] Found genes: {result['genes']}, traits: {result['traits']}")
    return result

# --- Merging and Post-processing Logic ---

def _post_process_entities(genes: List[str], traits: List[str]) -> Dict[str, List[str]]:
    """Cleans and filters lists of genes and traits."""
    
    def clean_list(items: List[str]) -> List[str]:
        cleaned_items = []
        for item in items:
            cleaned_items.extend(sub_item.strip() for sub_item in item.split(","))
        return cleaned_items

    processed_genes = clean_list(genes)
    processed_traits = clean_list(traits)

    processed_genes = [g for g in processed_genes if "no gene" not in g.lower()]
    processed_genes = [g for g in processed_genes if 2 <= len(g) <= 20]
    processed_genes = [g for g in processed_genes if g.lower() != "mirna"]
    
    processed_traits = [t for t in processed_traits if len(t) > 2]
    
    return {
        "genes": list(dict.fromkeys(processed_genes)),
        "traits": list(dict.fromkeys(processed_traits))
    }


def _extract_entities_merged(text: str) -> Dict[str, Any]:
    """Extracts entities using all models, merges the results, and post-processes them."""
    claude_result = _extract_entities_with_claude(text)
    flair_result = _extract_entities_with_flair(text)
    gliner_result = _extract_entities_with_gliner(text)
    variant_result = _extract_variants_with_regex(text)
    
    raw_genes = claude_result["genes"] + flair_result["genes"] + gliner_result["genes"]
    raw_traits = flair_result["traits"] + gliner_result["traits"]
    
    processed_entities = _post_process_entities(raw_genes, raw_traits)

    final_result = {
        "genes": processed_entities["genes"],
        "traits": processed_entities["traits"],
        "variants": variant_result["variants"],
        "chr_pos_variants": variant_result["chr_pos_variants"]
    }
    
    if DEBUG:
        print(f"[_extract_entities_merged] Final Merged Genes: {final_result['genes']}")
        print(f"[_extract_entities_merged] Final Merged Traits: {final_result['traits']}")
        print(f"[_extract_entities_merged] Final Merged Variants: {final_result['variants']}")

    return final_result


# --- Agent Nodes ---

def entity_extraction_node(state: "State") -> "State":
    """Runs the comprehensive, merged entity extraction process."""
    user_input: str = state["messages"][-1]["content"]
    extraction_result = _extract_entities_merged(user_input)
    if DEBUG:
        print(f"[entity_extraction_node] Extracted: {extraction_result}")
    return extraction_result


def claude_entity_extraction_node(state: "State") -> "State":
    """Extracts gene entities using only the Claude model."""
    user_input: str = state["messages"][-1]["content"]
    claude_result = _extract_entities_with_claude(user_input)
    processed_result = _post_process_entities(claude_result["genes"], claude_result["traits"])
    if DEBUG:
        print(f"[claude_entity_extraction_node] Extracted: {processed_result}")
    return processed_result


def gliner_entity_extraction_node(state: "State") -> "State":
    """Extracts all entity types using only the GLiNER model."""
    user_input: str = state["messages"][-1]["content"]
    gliner_result = _extract_entities_with_gliner(user_input)
    processed_result = _post_process_entities(gliner_result["genes"], gliner_result["traits"])
    
    final_output = {
        "genes": processed_result["genes"],
        "traits": processed_result["traits"],
    }
    if DEBUG:
        print(f"[gliner_entity_extraction_node] Extracted: {final_output}")
    return final_output


def flair_entity_extraction_node(state: "State") -> "State":
    """Extracts gene and trait entities using only the Flair model."""
    user_input: str = state["messages"][-1]["content"]
    flair_result = _extract_entities_with_flair(user_input)
    processed_result = _post_process_entities(flair_result["genes"], flair_result["traits"])
    if DEBUG:
        print(f"[flair_entity_extraction_node] Extracted: {processed_result}")
    return processed_result

def variant_extraction_node(state: "State") -> "State":
    """Extracts variant entities using only regex."""
    user_input: str = state["messages"][-1]["content"]
    variant_result = _extract_variants_with_regex(user_input)
    if DEBUG:
        print(f"[variant_extraction_node] Extracted: {variant_result}")
    return variant_result
    

# --- Graph Edges ---

def has_genes(state: "State") -> bool:
    """Edge condition helper: returns `True` if any genes were found."""
    found = bool(state.get("genes"))
    if DEBUG:
        print(f"[has_genes] Genes present? {found}")
    return found


def has_traits(state: "State") -> bool:
    """Edge condition helper: returns `True` if any traits were found."""
    found = bool(state.get("traits"))
    if DEBUG:
        print(f"[has_traits] Traits present? {found}")
    return found

def has_variants(state: "State") -> bool:
    """Edge condition helper: returns `True` if any variants were found."""
    found = bool(state.get("variants"))
    if DEBUG:
        print(f"[has_variants] Variants present? {found}")
    return found
