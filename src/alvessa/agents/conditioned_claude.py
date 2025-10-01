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

from src.alvessa.clients.claude import claude_call
from src.config import CONDITIONED_MODEL, DEBUG, N_CHARS
from src.state import State
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


# Data aggregation helpers
def _extract_gene_data(state: "State", gene: str) -> Dict[str, Any]:
    """Extract all data for a single gene from state."""
    
    gene_objs = state.get("gene_entities", {})
    gene_obj = gene_objs.get(gene)
    print('****'*10)
    return gene_obj.summarize_text()


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


# --- New helpers for deterministic titles and pure text docs ---
def canon_gene_key(symbol: str) -> str:
    """Normalize gene keys to uppercase with trimmed whitespace (local fallback)."""
    return (symbol or "").strip().upper()

def _removeprefix(s: str, prefix: str) -> str:
    # Py<3.9 compat
    return s[len(prefix):] if s.startswith(prefix) else s

def _stable_gene_title(symbol: str) -> str:
    return f"GENE: {canon_gene_key(symbol)}" if symbol else "GENE: UNKNOWN"

def _stable_variant_title(rsid: str | None) -> str:
    return f"VAR: {rsid}" if rsid else "VAR: UNKNOWN"

def variant_text_document(v: "Variant") -> str:
    """
    Pure text (no in-place mutation) for sentence-chunked citation.
    Keep it simple and sentence-y so Claude's auto sentence chunking works well.
    """
    lines: list[str] = []

    if v.rsID:
        lines.append(f"Variant ID: {v.rsID}.")
    # locations (one sentence per build)
    for build, loc in (v.loc_by_build or {}).items():
        chrom = loc.get("chrom", "")
        pos = loc.get("pos", "")
        ref = ",".join(loc.get("ref", []) or []) if isinstance(loc.get("ref"), list) else (loc.get("ref") or "")
        alt = ",".join(loc.get("alt", []) or []) if isinstance(loc.get("alt"), list) else (loc.get("alt") or "")
        bits = [f"Location {build}: chr{chrom}:{pos}"]
        if ref: bits.append(f"Ref={ref}")
        if alt: bits.append(f"Alt={alt}")
        lines.append((", ".join(bits)) + ".")
    if v.organism:
        lines.append(f"Organism: {v.organism}.")

    # per-gene context (one sentence each)
    for gene, ctx in (v.per_gene_context or {}).items():
        category = ctx.get("variant_category", "")
        context = ctx.get("context", "")
        pieces = [f"Gene {gene}"]
        if category: pieces.append(f"category {category}")
        if context: pieces.append(f"context {context}")
        lines.append((", ".join(pieces)) + ".")

    # traits (short)
    if v.traits:
        lines.append(f"Associated traits: {', '.join(sorted(set(map(str, v.traits))))}.")

    # functional predictions: compress to short sentences by tool over genes
    if v.functional_predictions:
        for gene, tools in v.functional_predictions.items():
            tool_bits = []
            for tool, scores in (tools or {}).items():
                # shorten large lists conservatively
                s = scores if isinstance(scores, list) else [scores]
                s_str = ", ".join(map(lambda x: str(x)[:24], s[:5]))
                extra = "" if len(s) <= 5 else f" (+{len(s)-5} more)"
                tool_bits.append(f"{tool}: {s_str}{extra}")
            if tool_bits:
                lines.append(f"Functional predictions for {gene}: " + "; ".join(tool_bits) + ".")

    # variant text summary lines (if present)
    for t in (v.variant_summaries or []):
        t = str(t).strip()
        if t and not t.endswith("."):
            t += "."
        if t:
            lines.append(t)

    # fall back if nothing
    if not lines:
        lines.append("No additional summary available.")
    return " ".join(lines)

def _collect_notes_doc(state: "State") -> list[dict]:
    text_notes = state.get("text_notes", []) or []
    # Simple prose blob; one sentence per note helps chunking
    if not text_notes:
        return []
    joined = " ".join([n.strip().rstrip(".") + "." for n in text_notes if str(n).strip()])
    if not joined:
        return []
    return [{
        "type": "document",
        "source": {"type": "text", "media_type": "text/plain", "data": joined},
        "title": "NOTES",
        "citations": {"enabled": True},
    }]

def _collect_gene_docs(state: "State") -> tuple[list[dict], list[tuple[str, int]]]:
    """
    Returns (document_blocks, manifest) where manifest holds (stable_id, document_index)
    so you can map citations back if you ever need to.
    """
    docs = []
    manifest = []

    gene_objs = state.get("gene_entities", {}) or {}
    # If user stored as list, normalize to dict-like view
    if isinstance(gene_objs, list):
        # assume list of symbols; you can adjust if needed
        gene_objs = {g: state.get("gene_lookup", {}).get(g) for g in gene_objs}

    # Dedup + deterministic order
    items = []
    for k, g in (gene_objs or {}).items():
        if not g:
            continue
        symbol = getattr(g, "symbol", None) or str(k)
        items.append((canon_gene_key(symbol), g))

    items.sort(key=lambda x: x[0])
    seen = set()
    for symbol, g in items:
        if symbol in seen:
            continue
        seen.add(symbol)
        title = _stable_gene_title(symbol)
        # keep your presenter, sentence-y bullets work fine
        text = g.summarize_text(include_go=False, include_pathways=True)
        # strip leading bullets to keep sentences crisp
        text = text.replace("\n• ", " ").replace("• ", "")
        docs.append({
            "type": "document",
            "source": {"type": "text", "media_type": "text/plain", "data": text},
            "title": title,
            "citations": {"enabled": True},
        })
        manifest.append((title, len(manifest)))  # index is positional after NOTES (will adjust later)
    return docs, manifest

def _collect_variant_docs(state: "State") -> tuple[list[dict], list[tuple[str, int]]]:
    docs = []
    manifest = []

    var_objs = state.get("variant_entities", {}) or {}
    # Normalize to dict of rsid->Variant
    items = []
    for k, v in (var_objs or {}).items():
        if not v:
            continue
        rsid = getattr(v, "rsID", None) or str(k)
        items.append((rsid, v))
    # deterministic order, unknowns last
    items.sort(key=lambda x: (x[0] is None, str(x[0]).upper()))

    seen = set()
    for rsid, v in items:
        key = rsid or f"NOID_{id(v)}"
        if key in seen:
            continue
        seen.add(key)
        title = _stable_variant_title(rsid)
        text = variant_text_document(v)
        docs.append({
            "type": "document",
            "source": {"type": "text", "media_type": "text/plain", "data": text},
            "title": title,
            "citations": {"enabled": True},
        })
        manifest.append((title, len(manifest)))
    return docs, manifest
def _anthropic_text_blocks(msg) -> list[str]:
    parts = []
    content = getattr(msg, "content", []) or []
    for blk in content:
        t = getattr(blk, "type", None) if not isinstance(blk, dict) else blk.get("type")
        if t == "text":
            txt = getattr(blk, "text", None) if not isinstance(blk, dict) else blk.get("text")
            if txt:
                parts.append(txt)
    return parts

def _anthropic_join_text(msg) -> str:
    return "".join(_anthropic_text_blocks(msg))

def _anthropic_collect_citations(msg) -> list[dict]:
    """
    Return a flat list of normalized citations.
    Each item: {
        'type': 'char_location' | 'page_location' | 'content_block_location',
        'document_index': int,
        'document_title': str | None,
        'cited_text': str | None,
        # char:
        'start_char_index': int | None, 'end_char_index': int | None,
        # page:
        'start_page_number': int | None, 'end_page_number': int | None,
        # block:
        'start_block_index': int | None, 'end_block_index': int | None
    }
    """
    out = []
    content = getattr(msg, "content", []) or []
    for blk in content:
        citations = getattr(blk, "citations", None)
        if isinstance(blk, dict):
            citations = blk.get("citations")
        if not citations:
            continue
        for c in citations:
            isdict = isinstance(c, dict)
            typ = (c.get("type") if isdict else getattr(c, "type", None)) or "unknown"
            rec = {
                "type": typ,
                "document_index": c.get("document_index") if isdict else getattr(c, "document_index", None),
                "document_title": c.get("document_title") if isdict else getattr(c, "document_title", None),
                "cited_text": c.get("cited_text") if isdict else getattr(c, "cited_text", None),
                "start_char_index": c.get("start_char_index") if isdict else getattr(c, "start_char_index", None),
                "end_char_index": c.get("end_char_index") if isdict else getattr(c, "end_char_index", None),
                "start_page_number": c.get("start_page_number") if isdict else getattr(c, "start_page_number", None),
                "end_page_number": c.get("end_page_number") if isdict else getattr(c, "end_page_number", None),
                "start_block_index": c.get("start_block_index") if isdict else getattr(c, "start_block_index", None),
                "end_block_index": c.get("end_block_index") if isdict else getattr(c, "end_block_index", None),
            }
            out.append(rec)
    return out

def _map_doc_index_to_title(manifest: list[dict]) -> dict[int, str]:
    # manifest items look like {'title': ..., 'index': ...}
    return {int(m["index"]): str(m["title"]) for m in manifest if "index" in m and "title" in m}

def _shorten(s: str, n: int = 140) -> str:
    s = (s or "").strip().replace("\n", " ")
    return s if len(s) <= n else s[: n - 1].rstrip() + "…"

def _plaintext_evidence_from_citations(citations: list[dict], idx2title: dict[int, str]) -> list[str]:
    """
    Make compact human-readable evidence lines from citations.
    Examples:
      "GENE: FTO — “Symbol: FTO | Entrez: …” [chars 0–114]"
      "VAR: rs123 — “Location GRCh38: chr1:12345, Ref=A, Alt=G”"
      "GENE: TP53 — [pages 5–6]"
    """
    lines = []
    seen = set()
    for c in citations:
        idx = c.get("document_index")
        title = idx2title.get(idx, f"doc#{idx}") if idx is not None else "doc"
        snippet = _shorten(c.get("cited_text") or "")
        typ = c.get("type")

        where = ""
        if typ == "char_location":
            a, b = c.get("start_char_index"), c.get("end_char_index")
            if a is not None and b is not None:
                where = f"[chars {a}–{b}]"
        elif typ == "page_location":
            a, b = c.get("start_page_number"), c.get("end_page_number")
            if a is not None and b is not None:
                where = f"[pages {a}–{b}]"
        elif typ == "content_block_location":
            a, b = c.get("start_block_index"), c.get("end_block_index")
            if a is not None and b is not None:
                where = f"[blocks {a}–{b}]"

        # Compose and dedupe
        if snippet:
            line = f"{title} — “{snippet}” {where}".strip()
        else:
            line = f"{title} {where}".strip()

        if line not in seen:
            seen.add(line)
            lines.append(line)
    return lines


def create_context_block(state: "State") -> str:
    """Build a compact JSON context. """
    
    if DEBUG:
        print("[create_context_block function] preparing context block...")
        
    # Fiest, use text_notes in state if available
    text_notes = state.get("text_notes", [])
    context_payload = [{"general_notes": text_notes}] if text_notes else []
    
    # Build gene-based context
    gene_list = list(set(state.get("gene_entities", [])))
    context_payload += [_extract_gene_data(state, gene) for gene in gene_list]

    if DEBUG:
        print(f"[create_context_block] genes to summarize: {gene_list}")
        
    # Add variant-based context if available
    variant_objs = state.get("variant_entities", {})
    variant_summaries = ""
    for id in variant_objs:
        variant_summary = variant_objs[id].return_full_summary()
        if variant_summary:
            variant_summaries += f"{variant_summary}\n"

    context_payload += [{"variant_summaries": variant_summaries.strip()}] if variant_summaries else []
    context_payload = ensure_json_safe(context_payload)
    return to_unquoted_top(context_payload)

def create_context_block_with_citations(state: "State") -> tuple[list[dict], list[dict]]:
    """
    Assemble: NOTES → GENE docs → VAR docs, with citations enabled.
    Returns (document_blocks, manifest_list) where manifest_list has
    [{'title': ..., 'kind': 'NOTES'|'GENE'|'VAR', 'key': ... , 'index': ...}]
    """
    blocks: list[dict] = []
    manifest: list[dict] = []

    # NOTES first
    notes_blocks = _collect_notes_doc(state)
    if notes_blocks:
        # Single NOTES doc
        blocks.extend(notes_blocks)
        manifest.append({"title": "NOTES", "kind": "NOTES", "key": "NOTES", "index": 0})

    # GENES
    gene_blocks, gene_manifest = _collect_gene_docs(state)
    start_idx = len(blocks)
    blocks.extend(gene_blocks)
    for i, (title, _) in enumerate(gene_manifest):
        manifest.append({"title": title, "kind": "GENE", "key": title.removeprefix("GENE: ").strip(), "index": start_idx + i})

    # VARIANTS
    var_blocks, var_manifest = _collect_variant_docs(state)
    start_idx = len(blocks)
    blocks.extend(var_blocks)
    for i, (title, _) in enumerate(var_manifest):
        manifest.append({"title": title, "kind": "VAR", "key": title.removeprefix("VAR: ").strip(), "index": start_idx + i})

    return blocks, manifest



def conditioned_claude_node(state: "State") -> "State":
    """
    Replacement for conditioned_claude_node that uses Claude's citations feature.
    - Builds documents (NOTES, GENE:*, VAR:*) with sentence-level auto chunking.
    - Sends the user question as a separate text block with an explicit instruction to use citations.
    - Does NOT require 'evidence' JSON anymore.
    """
    if DEBUG:
        print("[conditioned_claude_node_with_citations] building documents...")

    # Build docs
    doc_blocks, manifest = create_context_block_with_citations(state)

    # Instruction: keep it simple; ask Claude to use citations
    system_msg = state.get("prompt", "").strip()
    if not system_msg:
        system_msg = (
            "You are a research scientist. Answer strictly based on the provided documents. "
            "Be accurate and reason through the provided information step-by-step. "
            "Make sure you only use citations to back up your answer. "
            "No additional information should be added beyond what is in the documents. "
        )

    # Pull latest user question
    user_question = state["messages"][-1]["content"]

    # User turn content: documents + question text
    user_content = list(doc_blocks) + [
        {
            "type": "text",
            "text": f"User asked: {user_question}\n\nUse citations to back up your answer."
        }
    ]
    
    print(user_content)
    print('***')

    # Call Claude
    raw = claude_call(
        model=CONDITIONED_MODEL,
        temperature=0.1,
        max_tokens=20000,
        system=system_msg,
        messages=[{"role": "user", "content": user_content}],
    )
    print(raw)

        
    # 1) Plain text answer: join all text blocks
    answer_text = _anthropic_join_text(raw)

    # 2) Collect normalized citations
    citations = _anthropic_collect_citations(raw)
    

    # 3) Map doc indices -> titles and make human evidence lines
    idx2title = _map_doc_index_to_title(manifest)
    evidence_lines = _plaintext_evidence_from_citations(citations, idx2title)

    # 4) UI-friendly raw artifacts
    #    - 'llm_raw_repr' keeps the full object repr for debugging
    #    - 'llm_citations' is a compact, JSON-safe projection for your hyperlink UI
    llm_raw_repr = repr(raw)  # safe, always string
    llm_citations = ensure_json_safe(citations)  # JSON-safe list[dict]

    parsed_resp = {
    "answer": answer_text,
    "evidence": evidence_lines,        
    "doc_manifest": manifest,          
    "citations": llm_raw_repr,        
    "citations_enabled": True,        
    }
    
    return {'llm_json': parsed_resp}

# def conditioned_claude_node(state: "State") -> "State":
#     """
#     Build a compact JSON context and ask Claude for an answer.
    
#     This is the main reasoning engine that synthesizes data from all biomedical tools
#     into a coherent, evidence-based response.

#     Returns
#     -------
#     State
#         Updated state with assistant message, context block and parsed JSON.
#     """
#     context_block = create_context_block(state)

#     if DEBUG:
#         print(f"[conditioned_claude_node] context length before truncation: {len(context_block)}")
    
#     if len(context_block) > N_CHARS:
#         context_block = context_block[:N_CHARS] + "...<truncated>"
    
#     # Generate Claude response
#     user_question = state["messages"][-1]["content"]
#     system_msg = state.get('prompt', '')    
    
#     if len(system_msg)<2:
#         system_msg = (
#             "You are a research scientist. Answer the question strictly with facts you can "
#             "point to inside CONTEXT. Respond only with JSON with keys answer and evidence." 
#             " Be descriptive and detailed, think in steps and outline your process. Ensure proper JSON format. "
#             "The 'evidence' field must always be a list of short strings "
#             # "If the CONTEXT contains trait-based associations (query_type: 'trait_based'), focus on the genetic associations "
#             # "with the queried trait/disease, including related genes, variants, and their biological significance."
#         )
    
#     if DEBUG:
#         print(f"[conditioned_claude_node] system message: {system_msg}")
#         print(f"[conditioned_claude_node] user question: {user_question}")
#         print(f"[conditioned_claude_node] context block length: {len(context_block)}")
#         print(f"[conditioned_claude_node] full context block: {context_block}")
    
#     raw = claude_call(
#         model=CONDITIONED_MODEL,
#         temperature=0.1,
#         max_tokens=20000,
#         system=system_msg,
#         messages=[{"role": "user", "content": f"User asked: {user_question}\n\nCONTEXT:\n{context_block}"}],
#     )
    
#     if DEBUG:
#         print("[conditioned_claude_node] raw response from Claude:", raw)

#     # Parse Claude response
#     llm_resp = raw.content[0].text.strip() if hasattr(raw.content[0], "text") else raw.content[0]
    
#     if DEBUG:
#         print("[conditioned_claude_node] processed response text:", llm_resp)
    
#     if llm_resp.startswith("```"):
#         llm_resp = re.sub(r"^```(?:json)?\s*|\s*```$", "", llm_resp.strip(), flags=re.DOTALL).strip()

#     try:
#         parsed_resp = json.loads(llm_resp) if isinstance(llm_resp, str) else llm_resp
#     except Exception as exc:
#         parsed_resp = {
#             "answer": llm_resp, 
#             "evidence": '',
#         }
#     return {
#         "messages": [{"role": "assistant", "content": llm_resp}],
#         "context_block": context_block,
#         "llm_json": parsed_resp,
#     }
