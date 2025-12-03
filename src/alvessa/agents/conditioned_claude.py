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
    s = s.replace("\n", " ")   
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

def _stable_drug_title(label: str | None) -> str:
    return f"DRUG: {label}" if label else "DRUG: UNKNOWN"

def _clean_bullets(text: str) -> str:
    return text.replace("\n• ", " ").replace("• ", " ")

def variant_text_document(v: "Variant") -> str:
    """
    Pure text (no in-place mutation) for sentence-chunked citation.
    Keeping it simple and sentence-y so Claude's auto sentence chunking works well.
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
                if tool=='expectosc_predictions_agent': continue  # skip long Expecto outputs, they are added in summary
                # shorten large lists conservatively
                s = scores if isinstance(scores, list) else [scores]
                s_str = ", ".join(map(lambda x: str(x), s))
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
    so we can map citations back
    """
    docs = []
    manifest = []

    gene_objs = state.get("gene_entities", {}) or {}
    # If user stored as list, normalize to dict-like view
    if isinstance(gene_objs, list):
        # assume list of symbols
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
        text = g.summarize_text(include_go=False, include_pathways=True, include_interactions=True)
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


def _collect_drug_docs(state: "State") -> tuple[list[dict], list[tuple[str, int]]]:
    docs = []
    manifest = []

    drug_objs = state.get("drug_entities", {}) or {}

    items = []
    for k, d in (drug_objs or {}).items():
        if not d:
            continue
        ids = getattr(d, "identifiers", None)
        name = getattr(ids, "name", None) if ids else None
        canon_key = getattr(ids, "canon_key", None) if ids else None
        label = name or canon_key or str(k)
        items.append((label, d))

    items.sort(key=lambda x: (x[0] is None, str(x[0]).lower()))

    seen = set()
    for label, d in items:
        key = label or f"UNKNOWN_{id(d)}"
        if key in seen:
            # ensure uniqueness even if label is missing/duplicated
            key = f"{key}_{len(seen)}"
        seen.add(key)

        title = _stable_drug_title(key)
        text = _clean_bullets(d.summarize_text())
        docs.append({
            "type": "document",
            "source": {"type": "text", "media_type": "text/plain", "data": text},
            "title": title,
            "citations": {"enabled": True},
        })
        manifest.append((title, len(manifest)))

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


def _slice_span(doc_text_by_index: dict[int, str], idx, a, b, max_len: int = 10000) -> str | None:
    if idx is None or a is None or b is None:
        return None
    try:
        txt = doc_text_by_index.get(int(idx), "")
        if not txt:
            return None
        a = max(0, int(a))
        b = max(a, int(b))
        snippet = txt[a:b].replace("\n", " ").strip()
        if not snippet:
            return None
        return snippet if len(snippet) <= max_len else snippet[: max_len - 1].rstrip() + "…"
    except Exception:
        return None

def _anthropic_join_text(
    msg,
    doc_text_by_index: dict[int, str] | None = None,
    idx2title: dict[int, str] | None = None,
) -> str:
    def _as_dict(obj):
        if isinstance(obj, dict):
            return obj
        return {
            "type": getattr(obj, "type", None),
            "text": getattr(obj, "text", None),
            "citations": getattr(obj, "citations", None),
        }

    def _safe_str(x):
        return "" if x is None else str(x)

    def _short_proof(s: str, n: int = 240) -> str:
        s = s.strip().replace("\n", " ")
        return s if len(s) <= n else s[: n - 1].rstrip() + "…"

    parts = []
    content = getattr(msg, "content", []) or []
    for blk in content:
        b = _as_dict(blk)
        if b.get("type") != "text":
            continue

        text = _safe_str(b.get("text"))
        proofs = []
        citations = b.get("citations") or []

        for c in citations:
            isdict = isinstance(c, dict)
            typ   = (c.get("type") if isdict else getattr(c, "type", None)) or "unknown"
            d_idx = c.get("document_index") if isdict else getattr(c, "document_index", None)
            a     = c.get("start_char_index") if isdict else getattr(c, "start_char_index", None)
            b_    = c.get("end_char_index") if isdict else getattr(c, "end_char_index", None)
            title = _safe_str(c.get("document_title") if isdict else getattr(c, "document_title", None)).strip()
            cited = _safe_str(c.get("cited_text") if isdict else getattr(c, "cited_text", None)).strip()

            # Prefer exact span when available
            snippet = None
            if typ == "char_location" and doc_text_by_index is not None:
                snippet = _slice_span(doc_text_by_index, d_idx, a, b_, max_len=20000)

            # Fallback to model-provided cited_text
            if not snippet and cited:
                snippet = _short_proof(cited, n=240)

            # Resolve title via manifest if missing
            if (not title) and idx2title is not None and d_idx is not None:
                title = idx2title.get(int(d_idx), "")

            if snippet:
                # Include where-span so the join text is self-sufficient
                where = ""
                if typ == "char_location" and a is not None and b_ is not None:
                    where = f" [chars {a}–{b_}]"
                if title:
                    proofs.append(f"\"{snippet}\"@{title}{where}")
                else:
                    proofs.append(f"\"{snippet}\"{where}")

        proofs_repr = "NONE" if not proofs else "; ".join(proofs)
        parts.append((proofs_repr, text))

    return parts


def _map_doc_index_to_title(manifest: list[dict]) -> dict[int, str]:
    # manifest items look like {'title': ..., 'index': ...}
    return {int(m["index"]): str(m["title"]) for m in manifest if "index" in m and "title" in m}

def _shorten(s: str, n: int = 140) -> str:
    s = (s or "").strip().replace("\n", " ")
    return s if len(s) <= n else s[: n - 1].rstrip() + "…"


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

    # DRUGS
    drug_blocks, drug_manifest = _collect_drug_docs(state)
    start_idx = len(blocks)
    blocks.extend(drug_blocks)
    for i, (title, _) in enumerate(drug_manifest):
        manifest.append({"title": title, "kind": "DRUG", "key": title.removeprefix("DRUG: ").strip(), "index": start_idx + i})

    # VARIANTS
    var_blocks, var_manifest = _collect_variant_docs(state)
    start_idx = len(blocks)
    blocks.extend(var_blocks)
    for i, (title, _) in enumerate(var_manifest):
        manifest.append({"title": title, "kind": "VAR", "key": title.removeprefix("VAR: ").strip(), "index": start_idx + i})

    return blocks, manifest


def _strip_citations(joined: List[tuple[str, str]]) -> str:
    """
    Input: text produced by _anthropic_join_text, i.e. tuples of (proofs, text)
    We will join only the text parts, stripping out proof markup.
    """
    out_lines = ""
    for proofs, text in joined:
        out_lines += text
    return out_lines



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

    # Map indices -> titles and texts (for span reconstruction inside joiner)
    idx2title = _map_doc_index_to_title(manifest)
    doc_text_by_index = {}
    for m in manifest:
        i = int(m["index"])
        try:
            src = doc_blocks[i]["source"]
            if isinstance(src, dict) and src.get("type") == "text":
                doc_text_by_index[i] = str(src.get("data", ""))
            else:
                doc_text_by_index[i] = ""
        except Exception:
            doc_text_by_index[i] = ""

    # Instruction
    system_msg = state.get("prompt", "").strip()
    if not system_msg:
        system_msg = (
            "You are a research assistant. Answer strictly from the provided documents.\n"
            "But be thorough and point out all relevant and potentially interesting information."
            "Requirements:\n"
            "- Every factual claim must be grounded in the documents and include citations.\n"
            "- Do not add outside knowledge. If key information is missing, say so briefly.\n"
            "- Write clear, smooth prose (avoid filler words like 'Based on the documents' or long repetitive lists).\n"
            "- Lead with the direct answer to the question. Then provide supporting details in decreasing order of importance, with concise synthesis where useful (also cited).\n"
            "- Length: write as much as needed for completeness and clarity.\n"
            "- Factoid questions (Who/What/When/Where/Which/How many): if a single sentence fully answers the question, provide that one sentence. Avoid lengthy explanations unless necessary for clarity or context.\n"
            "- If the question is ambiguous, briefly note the ambiguity and address the most common interpretations.\n"
            "- If the question is multi-part, address each part clearly and separately.\n"
            "- When citing, prioritize attaching citations to complete sentences or clear factual units. Avoid breaking sentences unnaturally unless absolutely necessary (e.g., when two distinct factual claims appear within one sentence).\n"
            "- At the end, you may include a short section labeled **Possible speculation: if—and only if—it is clearly marked as such and explicitly reasoned from the cited evidence."
        )
        
    # check if there is verification in State alwready that is non-empty and if it is, we get the answer_text and retry_feedback
    prior_verification = state.get("verification", {})
    added_feedback = ""
    if prior_verification:
        answer_text = prior_verification.get("answer_text", "")
        retry_feedback = prior_verification.get("retry_feedback", "")
        if answer_text and retry_feedback:
            added_feedback = f" Previous answer: {answer_text}. It was declined as low quality, with feedback: {retry_feedback}."

    # Pull latest user question
    user_question = state["messages"][-1]["content"]
        
    print('USER QUESTION')
    print(user_question)
    
    if added_feedback=="":
        
        # User turn content: documents + question text
        user_content = list(doc_blocks) + [
            {
                "type": "text",
                "text": f"{user_question}\n Use provided information to answer."
            }
        ]
    else:
        # User turn content: documents + question text + feedback
        user_content = list(doc_blocks) + [
            {
                "type": "text",
                "text": f"{user_question}\n Use provided information to answer. Additionally, consider the following feedback as you tried to answer this quesiton before:\n {added_feedback}"
            }
        ]
    
    print('USER CONTENT')
    print(user_content)
    print('***')

    # Call Claude
    raw = claude_call(
        model=CONDITIONED_MODEL,
        temperature=0.0,
        max_tokens=20000,
        system=system_msg,
        messages=[{"role": "user", "content": user_content}],
    )
    print(raw)

    answer_with_proofs = _anthropic_join_text(
            raw,
            doc_text_by_index=doc_text_by_index,
            idx2title=idx2title,
        )
    # answer_with_proofs is a list of tuples (proofs, text)
    answer_plain = _strip_citations(answer_with_proofs)

    parsed_resp = {
        "answer_with_proofs": answer_with_proofs,
        "answer": answer_plain,
        "manifest": manifest,          
        "raw_answer": raw,
    }
    
    return {"llm_json": parsed_resp}
