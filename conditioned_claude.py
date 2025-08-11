"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-06-25
Updated: 2025-06-26


Description: 

LLM node that answers the user question using compiled context: main reasoning LLM"""

from __future__ import annotations
import json
from typing import Dict, List, Any

from claude_client import claude_call
from config import CONDITIONED_MODEL, DEBUG, N_CHARS
from state import State


def conditioned_claude_node(state: "State") -> "State":
    """
    Build a compact JSON context and ask Claude for an answer.

    Returns
    -------
    State
        Updated state with assistant message, context block and parsed JSON.
    """
    if DEBUG:
        print("[conditioned_claude_node] preparing context block...")
    
    # Build CONTEXT payload
    gene_payload: List[Dict[str, Any]] = []
    gene_list = list(set(state.get("genes", [])))
    
    if DEBUG:
        print("[conditioned_claude_node] genes to summarize:", gene_list)
    
    gene_list = list(set(gene_list))  # Ensure unique genes
    
    for g in gene_list:
        gene_info: Dict[str, Any] = {"gene": g}

        diseases = state.get("gene_disease_traits", {}).get(g, [])
        if diseases:
            gene_info["diseases"] = diseases
            
        humanbase_hits = state.get("humanbase_predictions", {}).get(g, [])
        terms = []
        if humanbase_hits:
            terms = [hit["term"] for hit in humanbase_hits if "term" in hit]
            if terms:
                gene_info["functions"] = terms[:30] # TODO: update this to a more flexible limit

        biogrid_hits = state.get("biogrid_summarized_go", {}).get(g, [])
        if biogrid_hits:
            if DEBUG:
                print('[conditioned_claude_node] Found BioGRID GO terms for', g, biogrid_hits)
            gene_info["gene_ontology_terms_of_interacting_genes"] = biogrid_hits

        # Add associations to the gene info
        associations = state["gwas_associations"].get(g, [])
        if associations:
            gene_info["gwas_associations"] = associations

        # Add UniProt entries
        uniprot_entries_base = state.get("uniprot_entries_base", {}).get(g, [])
        if uniprot_entries_base:
            uniprot_entries_base.pop('go_terms', None)
            gene_info["uniprot_entries_base"] = uniprot_entries_base
        uniprot_entries_gwas = state.get("uniprot_entries_gwas", {}).get(g, [])
        if uniprot_entries_gwas:
            uniprot_entries_gwas.pop('go_terms', None)
            gene_info["uniprot_entries_gwas"] = uniprot_entries_gwas

        # Add sei predictions
        sei_effect_predictions = state.get("sei_predictions", {}).get(g, [])
        if sei_effect_predictions:
            #gene_info["sei_sequence_class_(regulatory_activity_classification_for_variant)_predictions"] = sei_effect_predictions
            gene_info["Regulatory activity role for the regions where variants were found. Defined computationally through Sei, a deep learning model that predicts transcription factors, histone marks and dnase, though clustering prediction over the genome and then assinging values. Reperesents role of the region"] = sei_effect_predictions
            
        # Add summary text from Expecto from HB
        humanbase_expecto = state.get("humanbase_expecto", {}).get(g, [])
        if humanbase_expecto:
            gene_info["Tissue-specific expression disruption predictions from sequence"] = humanbase_expecto['summary_text']
        
        # Add per-variant info
        tissue_expression_preds_variant_text_description = state.get("tissue_expression_preds_variant_text_description", {}).get(g, {})
        if tissue_expression_preds_variant_text_description:
            gene_info["Per variant GE modulation predictions:"] = tissue_expression_preds_variant_text_description
            
        gene_payload.append(gene_info)

    context_block: str = json.dumps(gene_payload, separators=(",", ":"))
    if DEBUG:
        print("[conditioned_claude_node] context length:", len(context_block))
    if len(context_block) > N_CHARS:
        context_block = context_block[:N_CHARS] + "...<truncated>"
    else:
        print("[conditioned_claude_node] context", context_block)
    print(context_block)
    
    # Call Claude
    system_msg: str = (
        "You are a biology data analyst. Answer strictly with facts you can "
        "point to inside CONTEXT. Respond only with JSON with keys answer and evidence. Ensure proper JSON format. "
        "Produce raw json output. I don't want markdown."
        "The 'evidence' field must always be a list of short strings, and always reference the entity to which you are referring."
    )
    user_question: str = state["messages"][-1]["content"]

    raw = claude_call(
        model=CONDITIONED_MODEL,
        temperature=0,
        max_tokens=20000,
        system=system_msg,
        messages=[{"role": "user", "content": f"User asked: {user_question}\n\nCONTEXT:\n{context_block}"}],
    )

    if hasattr(raw.content[0], "text"):
        llm_resp: str | Dict[str, Any] = raw.content[0].text.strip()
    else:
        llm_resp = raw.content[0]

    try:
        parsed_resp: Dict[str, Any]
        if isinstance(llm_resp, dict):
            parsed_resp = llm_resp
        else:
            parsed_resp = json.loads(llm_resp)
    except Exception as exc:  # Fall back to raw JSON string
        raise ValueError(
            f"Failed to parse LLM response as JSON. Response was:\n{llm_resp}\nError: {exc}"
        ) from exc

    return {
        **state,
        "messages": state["messages"] + [{"role": "assistant", "content": llm_resp}],
        "context_block": context_block,
        "llm_json": parsed_resp,
    }
