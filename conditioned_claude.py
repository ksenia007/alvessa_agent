"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-06-25
Updated: 2025-08-26


Description: 

LLM node that answers the user question using compiled context: main reasoning LLM"""

from __future__ import annotations
import json
from typing import Dict, List, Any, Optional

from claude_client import claude_call
from config import CONDITIONED_MODEL, DEBUG, N_CHARS
from state import State
import re
import numpy as np

def ensure_json_safe(x):
    if isinstance(x, dict):
        return {ensure_json_safe(k): ensure_json_safe(v) for k, v in x.items()}
    if isinstance(x, list):
        return [ensure_json_safe(v) for v in x]
    if isinstance(x, tuple):
        return [ensure_json_safe(v) for v in x]  # JSON has no tuple
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if hasattr(x, "item"):  # NumPy scalar fallback
        return x.item()
    return x

def _process_dbsnp_variants(state: "State") -> "State":
    """Process dbSNP variants to remove the `matches` field that is populated by gencode which is too verbose."""
    for variant_id in state.keys():
        if "annotations" in state[variant_id]:
            for annotation_i in state[variant_id]["annotations"]:
                try:
                    annotation_i.pop("matches")
                except KeyError:
                    pass
    return state


# Data aggregation helpers
def _extract_gene_data(state: "State", gene: str) -> Dict[str, Any]:
    """Extract all data for a single gene from state."""

    gene_info = {"gene": gene}
    
    # Define data sources with their state keys and optional processing
    data_sources = [
        ('Summary about the gene structure, transcripts and complexity from GENCODE:', 'gene_level_gencode'),
        ("Diseases associated with this gene:", "gene_disease_traits"),
        ("Top 30 predicted gene functions (from HumanBase database):", "humanbase_predictions", lambda hits: [hit["term"] for hit in hits if "term" in hit][:30]),
        ("Gene ontology (GO) summarized terms of interacting genes: ", "biogrid_summarized_go"),
        ("Interacting human genes based on BioGRID curated database, organized by experimental system", "biogrid_interaction_groups"),
        ("Interacting non-human genes based on BioGRID curated database, organized by experimental system", "biogrid_interactions_select_nonhuman"),
        ("Associated Reactome pathways (curated biological pathways which describe how molecules interact within a cell to carry out different biological processes)", "reactome_pathways"),
        ("GWAS-based associations and their statistical significance for the gene and its variants, for example connection to traits, protein levels etc. Returns associations that are high-risk or high-signifiance or both.", "gwas_associations"),
        ("Uniprot-derived annotations for the given gene, diseases and functions:", "uniprot_entries_base", lambda data: {k: v for k, v in data.items() if k != 'go_terms'}),
        # ("uniprot_entries_gwas", "uniprot_entries_gwas", lambda data: {k: v for k, v in data.items() if k != 'go_terms'}),
        ("Regulatory activity role for the regions where variants were found. Defined computationally through Sei, a deep learning model that predicts transcription factors, histone marks and dnase, though clustering prediction over the genome and then assinging values. Reperesents role of the region", "sei_predictions"),
        ("Gene expression context of the genes, gives context for the ranges of variants that *could* be observed in this gene:", "humanbase_expecto", lambda data: data.get('summary_text')),
        ("Per variant gene expression modulation predictions, linked to variants of interest. Note that it is z-scored to 1000Genomes, so values below 1 are a fairly small effect:", "tissue_expression_preds_variant_text_description", lambda txt: txt if isinstance(txt, str) and txt.strip() else None),
        ("Pathogenicity predictions for each missense variant of interest. Computed through AlphaMissense, which predicts the likelihood that missense variants (genetic mutations where a single amino acid in a protein is changed) can cause disease", "alphamissense_predictions"),
        ("dbSNP variant annotations (genomic coordinates and allele frequencies from population studies)", "dbsnp_variants", _process_dbsnp_variants),
        ("dbSNP variant summary (rare vs common variants, chromosomes, assembly info)", "dbsnp_summaries", _process_dbsnp_variants),
        ("All computationally predicted gene targets for this miRNA (from the miRDB database):", "mirDB_targets"),
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
                if DEBUG and state_key in ['biogrid_summarized_go', 'biogrid_predictions','reactome_pathways']:
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


def _clean_str(s: str) -> str:
    """Remove embedded newlines so they don't break the one-line-per-entry rule."""
    s = s.replace("\r", "")
    s = s.replace("\n", " ")   # or ", " if you want them comma-separated
    return re.sub(r"\s{2,}", " ", s).strip()

def _to_unquoted_inner(obj) -> str:
    if obj is None: return "null"
    if isinstance(obj, bool): return "true" if obj else "false"
    if isinstance(obj, (int, float)): return str(obj)
    if isinstance(obj, str): return _clean_str(obj)
    if isinstance(obj, (list, tuple, set)):
        return "[" + ", ".join(_to_unquoted_inner(x) for x in obj) + "]"
    if isinstance(obj, dict):
        return "{" + ", ".join(f"{k}: {_to_unquoted_inner(v)}" for k, v in obj.items()) + "}"
    return _clean_str(str(obj))

def to_unquoted_top(payload_list) -> str:
    """Exactly one newline between top-level entries, no extra newlines inside."""
    return "\n".join(_to_unquoted_inner(entry) for entry in payload_list)


def create_context_block(state: "State") -> str:
    """Build a compact JSON context. """
    
    if DEBUG:
        print("[create_context_block function] preparing context block...")
    
    # Build gene-based context
    gene_list = list(set(state.get("genes", [])))
    context_payload = [_extract_gene_data(state, gene) for gene in gene_list]

    if DEBUG:
        print(f"[create_context_block] genes to summarize: {gene_list}")
    
    # Add trait-based context if no genes but trait data exists
    trait_associations = state.get("trait_associations", {})
    if not gene_list and trait_associations.get("found", False):
        if DEBUG:
            print("[create_context_block] No genes found, adding trait-based associations")
        context_payload.append(_build_trait_context(trait_associations))
    
    # Serialize context and handle truncation
    marked_payload = []
    for entry in context_payload:
        marked_payload.append(entry)  # keep the entry
        if isinstance(entry, dict):
            if "gene" in entry:
                marker = f"end of entity {entry['gene']} description"
            elif "trait" in entry:
                marker = f"end of entity {entry['trait']} description"
            else:
                marker = "end of entity description"
        else:
            marker = "end of entity description"
        marked_payload.append({"_marker": marker})

    context_payload = ensure_json_safe(marked_payload)
    return to_unquoted_top(context_payload)


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
    context_block = create_context_block(state)

    if DEBUG:
        print(f"[conditioned_claude_node] context length before truncation: {len(context_block)}")
    
    if len(context_block) > N_CHARS:
        context_block = context_block[:N_CHARS] + "...<truncated>"
    
    # Generate Claude response
    user_question = state["messages"][-1]["content"]
    system_msg = state.get('prompt', '')    
    
    if len(system_msg)<2:
        system_msg = (
            "You are a biology data analyst. Answer strictly with facts you can "
            "point to inside CONTEXT. Respond only with JSON with keys answer and evidence. Ensure proper JSON format. "
            "Produce raw json output. I don't want markdown. "
            "The 'evidence' field must always be a list of short strings, and always reference the entity to which you are referring. "
            "If the CONTEXT contains trait-based associations (query_type: 'trait_based'), focus on the genetic associations "
            "with the queried trait/disease, including related genes, variants, and their biological significance."
        )
    
    if DEBUG:
        print(f"[conditioned_claude_node] system message: {system_msg}")
        print(f"[conditioned_claude_node] user question: {user_question}")
        print(f"[conditioned_claude_node] context block length: {len(context_block)}")
        print(f"[conditioned_claude_node] full context block: {context_block}")
    
    raw = claude_call(
        model=CONDITIONED_MODEL,
        temperature=0.1,
        max_tokens=20000,
        system=system_msg,
        messages=[{"role": "user", "content": f"User asked: {user_question}\n\nCONTEXT:\n{context_block}"}],
    )
    
    if DEBUG:
        print("[conditioned_claude_node] raw response from Claude:", raw)

    # Parse Claude response
    llm_resp = raw.content[0].text.strip() if hasattr(raw.content[0], "text") else raw.content[0]
    
    if DEBUG:
        print("[conditioned_claude_node] processed response text:", llm_resp)
    
    if llm_resp.startswith("```"):
        llm_resp = re.sub(r"^```(?:json)?\s*|\s*```$", "", llm_resp.strip(), flags=re.DOTALL).strip()

    try:
        parsed_resp = json.loads(llm_resp) if isinstance(llm_resp, str) else llm_resp
    except Exception as exc:
        parsed_resp = {
            "answer": llm_resp, 
            "evidence": '',
        }
    return {
        "messages": [{"role": "assistant", "content": llm_resp}],
        "context_block": context_block,
        "llm_json": parsed_resp,
    }
