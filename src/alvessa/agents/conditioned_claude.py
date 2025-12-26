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
from src.config import CONDITIONED_MODEL, TOOL_SELECTOR_MODEL_BACKUP, DEBUG, N_CHARS
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


# === Token limit handling helpers ===

def _parse_token_error(exc_str: str) -> tuple[int, int] | None:
    """
    Parse error message like: 'prompt is too long: 217552 tokens > 200000 maximum'
    Returns: (current_tokens, max_tokens) or None if not a token limit error
    """
    match = re.search(r'(\d+) tokens > (\d+) maximum', str(exc_str))
    if match:
        return int(match.group(1)), int(match.group(2))
    return None


def _estimate_document_tokens(doc: dict) -> int:
    """
    Rough token estimate for a document block.
    Uses chars / 3 as conservative estimate 
    """
    try:
        src = doc.get("source", {})
        if isinstance(src, dict) and src.get("type") == "text":
            text = str(src.get("data", ""))
            return len(text) // 3  # Conservative estimate
    except:
        pass
    return 0


def _find_longest_documents(doc_blocks: list[dict], n: int = 1) -> list[int]:
    """
    Return indices of the N longest documents by token estimate.
    """
    doc_sizes = [(i, _estimate_document_tokens(doc)) for i, doc in enumerate(doc_blocks)]
    doc_sizes.sort(key=lambda x: x[1], reverse=True)  # Largest first
    return [idx for idx, _ in doc_sizes[:n] if idx is not None]


def _shorten_document(doc: dict, reduction_pct: float) -> dict:
    """
    Shorten a document by reduction_pct (0.0 to 1.0).
    Strategy: Keep first 70% of target length + last 30% of target length.
    This preserves critical info at start (names, IDs) and summaries at end.

    Returns a new document dict with shortened text.
    """
    try:
        src = doc.get("source", {})

        if DEBUG:
            print(f"    [_shorten_document] Source type: {type(src)}, is dict: {isinstance(src, dict)}")
            if isinstance(src, dict):
                print(f"    [_shorten_document] Source keys: {src.keys()}")
                print(f"    [_shorten_document] Source type field: {src.get('type')}")

        if not isinstance(src, dict) or src.get("type") != "text":
            if DEBUG:
                print(f"    [_shorten_document] Skipping: not a text document")
            return doc

        text = str(src.get("data", ""))
        current_len = len(text)
        target_len = int(current_len * (1.0 - reduction_pct))

        if DEBUG:
            print(f"    [_shorten_document] Current length: {current_len:,} chars")
            print(f"    [_shorten_document] Reduction pct: {reduction_pct:.1%}")
            print(f"    [_shorten_document] Target length: {target_len:,} chars")

        if target_len >= current_len:
            if DEBUG:
                print(f"    [_shorten_document] No reduction needed (target >= current)")
            return doc  # No reduction needed

        if target_len <= 0:
            if DEBUG:
                print(f"    [_shorten_document] ERROR: Target length is <= 0!")
            return doc

        # Split: 70% from start, 30% from end
        keep_start = int(target_len * 0.7)
        keep_end = int(target_len * 0.3)

        if DEBUG:
            print(f"    [_shorten_document] Keeping: first {keep_start:,} + last {keep_end:,} chars")

        shortened = (
            text[:keep_start] +
            "\n\n[... middle section truncated to fit token limit ...]\n\n" +
            text[-keep_end:]
        )

        shortened_len = len(shortened)
        if DEBUG:
            print(f"    [_shorten_document] Result length: {shortened_len:,} chars")
            print(f"    [_shorten_document] Actual reduction: {(current_len - shortened_len) / current_len:.1%}")

        # Create new doc (don't mutate original)
        new_doc = dict(doc)
        new_src = dict(src)
        new_src["data"] = shortened
        new_doc["source"] = new_src
        return new_doc

    except Exception as e:
        if DEBUG:
            print(f"    [_shorten_document] Exception: {e}")
            import traceback
            traceback.print_exc()
        return doc


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


def _slice_span(doc_text_by_index: dict[int, str], idx, a, b, max_len: int = 100000) -> str | None:
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
                snippet = _short_proof(cited, n=24000)

            # Resolve title via manifest if missing
            if (not title) and idx2title is not None and d_idx is not None:
                title = idx2title.get(int(d_idx), "")

            if snippet:
                where = ""
                if typ == "char_location" and a is not None and b_ is not None:
                    where = f" [chars {a}–{b_}]"
                idx_marker = f" [doc {d_idx}]" if d_idx is not None else ""
                if title:
                    proofs.append(f"\"{snippet}\"@{title}{where}{idx_marker}")
                else:
                    proofs.append(f"\"{snippet}\"{where}{idx_marker}")

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
            "Be thorough and point out all relevant and potentially interesting information."
            "Requirements:\n"
            "- Every factual claim must be grounded in the documents and include associated citations linking it to the document.\n"
            "- When a claim is supported by multiple facts from the documents, cite each supporting "
            "passage separately as a list of citations rather than citing only the first occurrence or a large chunk of text.\n"
            " - When citing, attach citations to complete sentences or clearly delimited factual units. Do not place citations on section headers or allow a single citation to apply across multiple sections. Avoid breaking sentences unnaturally unless needed for citations.\n"
            "- Do not add outside knowledge. If key information is missing, say so briefly.\n"
            "- Lead with the direct answer to the question. Then provide supporting details in decreasing order of importance, with concise synthesis where useful (also cited).\n"
            "- Length: write as much as needed for completeness and clarity.\n"
            "- Factoid questions (Who/What/When/Where/Which/How many): if a single sentence fully answers the question, provide that one sentence. Avoid lengthy explanations unless necessary for clarity or context.\n"
            "- If the question is ambiguous, briefly note the ambiguity and address the most common interpretations.\n"
            "- At the end, you may include a short section labeled **Possible speculation**: if and only if it is clearly marked as such and explicitly reasoned from the cited evidence."
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

    if DEBUG:
        print('[conditioned_claude_node] USER QUESTION:')
        print(user_question)
        print('---')

    # Call Claude with retry logic for token limit errors
    # NOTE: user_content is built inside the retry loop with doc_blocks_working
    MAX_RETRIES = 2
    doc_blocks_working = list(doc_blocks)  # Mutable copy for shortening
    shortened_docs = []  # Track what we shortened
    model_used = CONDITIONED_MODEL
    used_backup = False
    raw = None
    if DEBUG:
        doc_lengths = [_estimate_document_tokens(doc) for doc in doc_blocks_working]
        doc_count = len(doc_lengths)
        mean_len = (sum(doc_lengths) / doc_count) if doc_count else 0
        max_len = max(doc_lengths) if doc_lengths else 0
        print("[conditioned_claude_node] Doc stats:")
        print(f"  Count: {doc_count}")
        print(f"  Mean est tokens: {mean_len:.1f}")
        print(f"  Max est tokens: {max_len}")
        if doc_lengths:
            longest_idx = max(range(len(doc_lengths)), key=lambda i: doc_lengths[i])
            longest_doc = doc_blocks_working[longest_idx]
            src = longest_doc.get("source", {})
            if isinstance(src, dict) and src.get("type") == "text":
                text = str(src.get("data", ""))
                sections = text.split("*")
                db_char_counts = {}
                for section in sections:
                    section = section.strip()
                    if not section:
                        continue
                    first_word = section.split(None, 1)[0]
                    db_char_counts[first_word] = db_char_counts.get(first_word, 0) + len(section)
                print(f"  Longest doc index: {longest_idx}")
                print(f"  Longest doc sections: {sum(1 for s in sections if s.strip())}")
                if db_char_counts:
                    top_db, top_chars = max(db_char_counts.items(), key=lambda x: x[1])
                    print(f"  Top db in longest doc by chars: {top_db} ({top_chars} chars)")

    for attempt in range(MAX_RETRIES + 1):
        # Rebuild user_content with current doc_blocks_working
        if added_feedback == "":
            user_content = list(doc_blocks_working) + [
                {
                    "type": "text",
                    "text": f"{user_question}\n Use provided information to answer."
                }
            ]
        else:
            user_content = list(doc_blocks_working) + [
                {
                    "type": "text",
                    "text": f"{user_question}\n Use provided information to answer. Additionally, consider the following feedback as you tried to answer this question before:\n {added_feedback}"
                }
            ]

        try:
            # Try primary model
            raw = claude_call(
                model=CONDITIONED_MODEL,
                temperature=0.0,
                max_tokens=20000,
                system=system_msg,
                messages=[{"role": "user", "content": user_content}],
            )
            model_used = CONDITIONED_MODEL
            used_backup = False
            break  # Success!

        except Exception as exc:
            exc_str = str(exc)
            token_info = _parse_token_error(exc_str)

            # Check if it's a "prompt too long" error and we have retries left
            if token_info and attempt < MAX_RETRIES:
                current_tokens, max_tokens = token_info
                overage = current_tokens - max_tokens
                overage_pct = overage / current_tokens
                # Add buffer to ensure we get under the limit
                reduction_needed = overage_pct + 0.40*(attempt+1)

                if DEBUG:
                    print(f"\n[conditioned_claude_node] PROMPT TOO LONG")
                    print(f"  Current: {current_tokens:,} tokens")
                    print(f"  Maximum: {max_tokens:,} tokens")
                    print(f"  Overage: {overage:,} tokens ({overage_pct*100:.1f}%)")
                    print(f"  Retry {attempt + 1}/{MAX_RETRIES}: reducing by {reduction_needed*100:.1f}%")

                # Find longest documents (on second retry, shorten 2)
                num_to_shorten = 1 if attempt == 0 else 2
                longest_indices = _find_longest_documents(doc_blocks_working, n=num_to_shorten)
                if not longest_indices:
                    if DEBUG:
                        print("ERROR: No documents to shorten!")
                    raise

                for longest_idx in longest_indices:
                    old_doc = doc_blocks_working[longest_idx]
                    old_tokens = _estimate_document_tokens(old_doc)
                    title = old_doc.get("title", f"Doc #{longest_idx}")

                    # Shorten it
                    new_doc = _shorten_document(old_doc, reduction_needed)
                    new_tokens = _estimate_document_tokens(new_doc)
                    doc_blocks_working[longest_idx] = new_doc

                    shortened_docs.append({
                        "title": title,
                        "attempt": attempt + 1,
                        "old_tokens": old_tokens,
                        "new_tokens": new_tokens,
                        "reduction_pct": (old_tokens - new_tokens) / old_tokens if old_tokens > 0 else 0
                    })

                    if DEBUG:
                        print(f"Shortened '{title}':")
                        print(f"    {old_tokens:,} -> {new_tokens:,} tokens ({(old_tokens-new_tokens)/old_tokens*100:.1f}% reduction)")

                continue  # Retry with shortened docs

            # Not a token error, or out of retries - try backup model
            if DEBUG:
                print(f"[conditioned_claude_node] Primary model failed: {exc}")

            if TOOL_SELECTOR_MODEL_BACKUP and attempt == 0:  # Only try backup on first failure
                try:
                    raw = claude_call(
                        model=TOOL_SELECTOR_MODEL_BACKUP,
                        temperature=0.0,
                        max_tokens=20000,
                        system=system_msg,
                        messages=[{"role": "user", "content": user_content}],
                    )
                    model_used = TOOL_SELECTOR_MODEL_BACKUP
                    used_backup = True
                    break  # Success with backup!
                except Exception as exc2:
                    if DEBUG:
                        print(f"[conditioned_claude_node] Backup model also failed: {exc2}")
                    raise
            else:
                raise  # Out of retries or no backup

    if shortened_docs and DEBUG:
        print(f"\n[conditioned_claude_node] Document shortening summary:")
        for s in shortened_docs:
            print(f"  Attempt {s['attempt']}: '{s['title']}' reduced by {s['reduction_pct']*100:.1f}%")

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
        # pass raw documents so verifier can slice by char_range
        "documents": {i: t for i, t in doc_text_by_index.items()},
        "raw_answer": raw,
        "model_used": model_used,
    }
    
    return {"llm_json": parsed_resp, "llm_backup_used": used_backup}
