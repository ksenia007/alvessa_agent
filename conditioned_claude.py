"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-06-25
Updated: 2025-06-26


Description: 

LLM node that answers the user question using compiled context: main reasoning LLM"""

from __future__ import annotations
import json
from typing import Dict, List, Any, Optional

from claude_client import claude_call
from config import CONDITIONED_MODEL, DEBUG, N_CHARS
from state import State
import re


# Data aggregation helpers
def _extract_gene_data(state: "State", gene: str) -> Dict[str, Any]:
    """Extract all data for a single gene from state."""
    gene_info = {"gene": gene}
    
    # Define data sources with their state keys and optional processing
    data_sources = [
        ("diseases", "gene_disease_traits"),
        ("functions", "humanbase_predictions", lambda hits: [hit["term"] for hit in hits if "term" in hit][:30]),
        ("gene_ontology_terms_of_interacting_genes", "biogrid_summarized_go"),
        ("Interacting genes based on BioGRID curated database", "biogrid_predictions"),
        ("Associated Reactome pathways (curated biological pathways which describe how molecules interact within a cell to carry out different biological processes)", "reactome_pathways"),
        ("gwas_associations", "gwas_associations"),
        ("uniprot_entries_base", "uniprot_entries_base", lambda data: {k: v for k, v in data.items() if k != 'go_terms'}),
        ("uniprot_entries_gwas", "uniprot_entries_gwas", lambda data: {k: v for k, v in data.items() if k != 'go_terms'}),
        ("Regulatory activity role for the regions where variants were found. Defined computationally through Sei, a deep learning model that predicts transcription factors, histone marks and dnase, though clustering prediction over the genome and then assinging values. Reperesents role of the region", "sei_predictions"),
        ("Tissue-specific expression disruption predictions from sequence", "humanbase_expecto", lambda data: data.get('summary_text')),
        ("Per variant GE modulation predictions:", "tissue_expression_preds_variant_text_description"),
        ("Pathogenicity predictions for each missense variant of interest. Computed through AlphaMissense, which predicts the likelihood that missense variants (genetic mutations where a single amino acid in a protein is changed) can cause disease", "alphamissense_predictions"),
        ("dbSNP variant annotations (genomic coordinates and allele frequencies from population studies)", "dbsnp_variants"),
        ("dbSNP variant summary (rare vs common variants, chromosomes, assembly info)", "dbsnp_summaries")
    ]
    
    for source in data_sources:
        field_name, state_key = source[:2]
        processor = source[2] if len(source) > 2 else None
        
        data = state.get(state_key, {}).get(gene)
        if data:
            if processor:
                data = processor(data)
            if data:  # Only add if data exists after processing
                gene_info[field_name] = data
                if DEBUG and state_key in ['biogrid_summarized_go', 'biogrid_predictions', 'reactome_pathways']:
                    print(f'[conditioned_claude_node] Found {state_key} for {gene}: {data}')
    
    return gene_info


def _build_trait_context(trait_associations: Dict[str, Any]) -> Dict[str, Any]:
    """Build trait-based context from trait associations."""
    trait_info = {
        "query_type": "trait_based",
        "trait_term": trait_associations.get("trait_term", ""),
        "total_associations": trait_associations.get("total_associations", 0),
        "total_significant_associations": trait_associations.get("total_significant_associations", 0),
        "total_studies_analyzed": trait_associations.get("total_studies_analyzed", 0)
    }
    
    # Add summary data with consistent structure
    for summary_type, prefix in [("summary_by_high_risk_alleles", ""), ("summary_by_significance", "significant_")]:
        summary = trait_associations.get(summary_type, {})
        if summary:
            for key in ["related_genes", "high_risk_snps", "proteins", "disease_traits"]:
                if summary.get(key):
                    trait_info[f"{prefix}{key}"] = summary[key]
    
    # Add variant annotations
    variant_annotations = trait_associations.get("variant_annotations", {})
    if variant_annotations:
        trait_info["variant_annotations"] = variant_annotations
    
    return trait_info


def conditioned_claude_node(state: "State") -> "State":
    """
    Build a compact JSON context and ask Claude for an answer.
    
    This is the main reasoning engine that synthesizes data from all biomedical tools
    into a coherent, evidence-based response.

    Returns
    -------
    State
        Updated state with assistant message, context block and parsed JSON.
    """
    if DEBUG:
        print("[conditioned_claude_node] preparing context block...")
    
    # Build gene-based context
    gene_list = list(set(state.get("genes", [])))
    context_payload = [_extract_gene_data(state, gene) for gene in gene_list]
    
    if DEBUG:
        print(f"[conditioned_claude_node] genes to summarize: {gene_list}")
    
    # Add trait-based context if no genes but trait data exists
    trait_associations = state.get("trait_associations", {})
    if not gene_list and trait_associations.get("found", False):
        if DEBUG:
            print("[conditioned_claude_node] No genes found, adding trait-based associations")
        context_payload.append(_build_trait_context(trait_associations))
    
    # Context building is now handled by helper functions above

    # Serialize context and handle truncation
    context_block = json.dumps(context_payload, separators=(",", ":"))
    if DEBUG:
        print(f"[conditioned_claude_node] context length: {len(context_block)}")
        if len(context_block) <= N_CHARS:
            print(f"[conditioned_claude_node] context: {context_block}")
    
    if len(context_block) > N_CHARS:
        context_block = context_block[:N_CHARS] + "...<truncated>"
    
    # Generate Claude response
    user_question = state["messages"][-1]["content"]
    system_msg = (
        "You are a biology data analyst. Answer strictly with facts you can "
        "point to inside CONTEXT. Respond only with JSON with keys answer and evidence. Ensure proper JSON format. "
        "Produce raw json output. I don't want markdown. "
        "The 'evidence' field must always be a list of short strings, and always reference the entity to which you are referring. "
        "If the CONTEXT contains trait-based associations (query_type: 'trait_based'), focus on the genetic associations "
        "with the queried trait/disease, including related genes, variants, and their biological significance."
    )
    
    raw = claude_call(
        model=CONDITIONED_MODEL,
        temperature=0,
        max_tokens=20000,
        system=system_msg,
        messages=[{"role": "user", "content": f"User asked: {user_question}\n\nCONTEXT:\n{context_block}"}],
    )
    
    # Parse Claude response
    llm_resp = raw.content[0].text.strip() if hasattr(raw.content[0], "text") else raw.content[0]
    
    if llm_resp.startswith("```"):
        llm_resp = re.sub(r"^```(?:json)?\s*|\s*```$", "", llm_resp.strip(), flags=re.DOTALL).strip()

    
    try:
        parsed_resp = json.loads(llm_resp) if isinstance(llm_resp, str) else llm_resp
    except Exception as exc:
        raise ValueError(
            f"Failed to parse LLM response as JSON. Response was:\n{llm_resp}\nError: {exc}"
        ) from exc
    
    return {
        "messages": [{"role": "assistant", "content": llm_resp}],
        "context_block": context_block,
        "llm_json": parsed_resp,
    }
