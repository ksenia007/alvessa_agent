"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-06-25
Updated: 2025-10-20


Description: 
Self-verification node that checks if the answer matches context."""

from __future__ import annotations
import json
import re
import unicodedata
import math
from typing import Any, Dict, List, Tuple, Iterable, Set, Optional

from src.alvessa.clients.claude import claude_call
from src.config import VERIFY_MODEL, DEBUG, CONDITIONED_MODEL, MAX_TOKENS
from src.state import State

# -----------------------------
# Lightweight text utilities
# -----------------------------
_WORD_RE = re.compile(r"[A-Za-z0-9_\-]+")
_NUM_RE = re.compile(
    r"(?<![A-Za-z])"                               # no letter before
    r"(?:(?:\d{1,3}(?:,\d{3})+)|\d+)(?:\.\d+)?%?"  # number w/ optional decimals/percent
    r"(?![A-Za-z])"                                # no letter after
)

def _numbers(s: str) -> Set[str]:
    """Extract standalone numeric tokens; never digits that are part of gene-like tokens."""
    return set(_NUM_RE.findall(s or ""))


_CHARS_RE = re.compile(r"\[chars\s+(\d+)\s*[–-]\s*(\d+)\]", re.I)


# split a semicolon list but ignore semicolons inside quotes
def _split_proofs_safely(proofs_raw: str) -> list[str]:
    parts, buf, in_q = [], [], False
    for ch in proofs_raw or "":
        if ch == '"' and (not buf or buf[-1] != "\\"):
            in_q = not in_q
        if ch == ';' and not in_q:
            parts.append("".join(buf).strip()); buf = []
        else:
            buf.append(ch)
    if buf:
        parts.append("".join(buf).strip())
    return [p for p in parts if p]

_CHARS_RE = re.compile(r"\[chars\s+(\d+)\s*[–-]\s*(\d+)\]", re.I | re.U)


def _parse_one_proof(p: str) -> Dict[str, Any]:
    p = (p or "").strip()
    snippet = ""
    title = None
    char_range = None

    first_quote = p.find('"')
    closing_marker = '"@'
    last_quote_before_title = p.rfind(closing_marker)
    if first_quote != -1 and last_quote_before_title != -1 and last_quote_before_title > first_quote:
        snippet = p[first_quote + 1:last_quote_before_title].strip()
    else:
        q = re.search(r'"(.*?)"', p, flags=re.DOTALL)
        if q:
            snippet = q.group(1).strip()

    t = re.search(r'@([^\[]+)$', p)  # everything after @ up to before '['
    if t:
        tval = t.group(1).strip()
        title = tval if tval else None

    cr = _CHARS_RE.search(p)
    if cr:
        a, b = int(cr.group(1)), int(cr.group(2))
        if a <= b:
            char_range = [a, b]

    return {
        "cited_text": snippet,
        "title": title,
        "document_index": None,  # will be resolved via manifest
        "char_range": char_range
    }

# turn [(proofs_str, text_str), ...] into [{"text":..., "proofs":[...]}...]
def _parse_answer_tuples(answer_with_proofs: list[tuple[str, str]]) -> List[Dict[str, Any]]:
    statements: List[Dict[str, Any]] = []
    for item in (answer_with_proofs or []):
        # Be defensive: accept lists, tuples, etc.
        try:
            proofs_raw, text = item
        except Exception:
            # if it’s a dict or malformed, skip safely
            continue
        proofs_raw = (proofs_raw or "").strip()
        text = (text or "").strip()

        proofs: List[Dict[str, Any]] = []
        if proofs_raw and proofs_raw.upper() != "NONE":
            for token in _split_proofs_safely(proofs_raw):
                proofs.append(_parse_one_proof(token))

        statements.append({"text": text, "proofs": proofs})
    return statements

def _link_titles_to_indices(proofs: List[Dict[str, Any]], manifest: List[Dict[str, Any]]) -> None:
    title2idx = {str(m["title"]).strip(): int(m["index"])
                 for m in manifest if "title" in m and "index" in m}
    for p in proofs:
        title = (p.get("title") or "").strip()
        if title and title in title2idx:
            p["document_index"] = title2idx[title]



def _link_titles_to_indices(proofs: List[Dict[str, Any]], manifest: List[Dict[str, Any]]) -> None:
    title2idx = {str(m["title"]).strip(): int(m["index"])
                 for m in manifest if "title" in m and "index" in m}
    for p in proofs:
        title = (p.get("title") or "").strip()
        if title and title in title2idx:
            p["document_index"] = title2idx[title]

# -----------------------------------
# Deterministic verdict
# -----------------------------------
def _deterministic_verdict(text: str, proofs: List[Dict[str, Any]], idx2title: Dict[int, str]) -> Tuple[str, List[str]]:
    reasons: List[str] = []

    for p in proofs:
        di = p.get("document_index")
        if di is not None and di not in idx2title:
            reasons.append(f"bad-doc-index:{di}")

    # Union proof text + titles
    union_proof = " ".join(p.get("cited_text") or "" for p in proofs)
    union_title = " ".join((p.get("title") or "") for p in proofs)
    union_all = (union_proof + " " + union_title).strip()
    
    # for each statement, check if there are any rsid* or ids starting with ENS* or r-hsa-*  or HGNC:* MGI:* RGD:* mentioned that are missing in proofs, where * is the numeric part
    # we will use regex to find the # patterns after thense IDs to the proof, removing all punctuation - substring match, all on .lower()
    id_patterns = [
        r"\brs\d+\b",          # rs + digits
        r"\bENS[A-Z]\d+\b",    # ENS + one letter + digits (e.g. ENSG00000123456)
        r"\br-hsa-\d+\b",      # r-hsa- + digits
        r"\bHGNC:\d+\b",       # HGNC: + digits
        r"\bMGI:\d+\b",        # MGI: + digits
        r"\bRGD:\d+\b",        # RGD: + digits
    ]

    text_lc = (text or "")
    union_lc = (union_all or "").lower()

    for pattern in id_patterns:
        for match in re.finditer(pattern, text_lc, flags=re.IGNORECASE):
            id_str = match.group(0)          # exactly the ID, nothing extra
            if id_str.lower() not in union_lc:
                reasons.append(f"missing-id:{id_str}")


    # Verdict (categorical only)
    if "no-citations" in reasons:
        verdict = "unsupported"
    elif any(r.startswith(( "bad-doc-index")) for r in reasons):
        verdict = "partial"
    elif any(r.startswith("missing-id") for r in reasons):
        verdict = "partial"
    else:
        verdict = "supported"

    return verdict, reasons

# -----------------------------------
# Main feedback
# -----------------------------------
def _llm_feedback(statements: List[Dict[str, Any]], question: str) -> Optional[Dict[str, Any]]:
    """
    Returns a dict like:
    {
      "overall": {
        "support_quality": "high|medium|low",
        "summary": "...",
        "concerns": ["..."],
        "suggestions": ["..."]
      },
      "per_statement": [
        {"idx": 0, "feedback": "…", "labels": ["overclaiming","needs-citation"]},
        ...
      ]
    }
    """
    pack = []
    for i, s in enumerate(statements):
        pack.append({
            "idx": i,
            "text": s["text"],
            "proofs": [
                {"title": p.get("title"), "snippet": (p.get("cited_text") or "")}
                for p in (s.get("proofs") or [])
            ]
        })

    system_msg = (
        "You are a meticulous research assitant. Provide qualitative feedback on whether each statement of the following answer is well supported "
        "by the quoted snippets (proofs). Do NOT invent facts or add citations, and double check the numbers.\n"
        "Special case: statements that begin with 'Possible speculation:' are allowed as cautious synthesis IF they are "
        "clearly grounded in the cited evidence. Flag only if they overreach or contradict the proofs.\n"
        "Return STRICT JSON with this schema:\n"
        "{\n"
        "  \"overall\": {\"support_quality\": \"high|medium|low\", \"summary\": \"...\", \"concerns\": [\"...\"], \"suggestions\": [\"...\"]},\n"
        "  \"per_statement\": {\n"
        "    \"0\": {\"label\": \"supported|partial|unsupported|speculation-ok|speculation-overreach\", \"feedback\": \"...\"},\n"
        "    \"1\": {\"label\": \"partial\", \"feedback\": \"...\"}\n ..."
        "  }\n"
        "}\n"
        "The keys in per_statement correspond to the statement indices. Ensure every statement is addressed.\n"
        "The output will be parsed with json.loads(), so it MUST be valid JSON.\n"
        "Note: citations consisting of long lists or large chunks may be truncated. If you encounter truncated proofs, advise the user to manually check the specific field.\n"
        "Non-problematic transitional or connective text may be labeled as 'supported' with a brief explanation.\n"
        "Label guidance:\n"
        "- supported: fully backed by the provided proofs.\n"
        "- partial: generally correct but missing a minor referenced entity/number, or only loosely aligned with the proofs, or contains minor exaggeration.\n"
        "- unsupported: not supported or contradicted by the proofs; introduces key entities or numbers not present; clear overreach."
        " If a statement contains multiple factual claims and even one key claim is unsupported, label the entire statement as 'unsupported' and explain in the feedback which part is unsupported.\n"
        "- speculation-ok: starts with 'Possible speculation:' AND stays within what the proofs plausibly support.\n"
        "- speculation-overreach: starts with 'Possible speculation:' BUT extends beyond or contradicts the proofs.\n"
        "Overall support_quality heuristic:\n"
        "- high: most statements are 'supported', at most a few 'partial', and no major claims are 'unsupported'.\n"
        "- medium: a mix of 'supported' and 'partial' statements, and/or a small number of 'unsupported' statements that do not affect the main conclusions.\n"
        "- low: several 'unsupported' statements, especially if they concern core claims, or pervasive partial/weak support.\n"
        "Avoid any numeric scoring."
    )
    payload = {
        "question": question,
        "statements": pack,
    }
    user_msg = json.dumps(payload, ensure_ascii=True)
    
    # check the length of user_msg, and if it exceeds MAX_TOKENS, we will split and make multiple calls
    # approx tokens as # words * 1.33
    approx_tokens = int(len(user_msg.split()) * 1.5)
    if approx_tokens > MAX_TOKENS:
        if DEBUG:
            print(f"[_llm_feedback] Warning: user_msg approx tokens {approx_tokens} exceeds MAX_TOKENS {MAX_TOKENS}. Truncating statements.")
        # split into 2 or 3 chunks, and if still too long, truncate each statement text; call LLM on each chunk and aggregate results
        n_chunks = max(1, math.ceil(approx_tokens / MAX_TOKENS))
        chunk_size = max(1, math.ceil(len(pack) / n_chunks))
        all_parsed: List[Dict[str, Any]] = []
        for start in range(0, len(pack), chunk_size):
            chunk_pack = pack[start:start + chunk_size]
            chunk_payload = {
                "question": question,
                "statements": chunk_pack,
            }
            chunk_user_msg = json.dumps(chunk_payload, ensure_ascii=True)
            raw_chunk = claude_call(
                model=CONDITIONED_MODEL,
                temperature=0,
                max_tokens=13000,
                system=system_msg,
                messages=[{"role": "user", "content": chunk_user_msg}],
            )
            txt_chunk = ""
            content_chunk = getattr(raw_chunk, "content", None)
            if content_chunk and isinstance(content_chunk, list) and getattr(content_chunk[0], "type", None) == "text":
                txt_chunk = getattr(content_chunk[0], "text", "") or ""
            else:
                txt_chunk = str(raw_chunk)
            cleaned_chunk = re.sub(r"```json\s*", "", txt_chunk)
            cleaned_chunk = re.sub(r"```", "", cleaned_chunk)
            try:
                parsed_chunk = json.loads(cleaned_chunk)
                all_parsed.append(parsed_chunk)
            except Exception:
                if DEBUG:
                    print("[_llm_feedback] JSON parse error in chunk verifier")
        if all_parsed:
            return _merge_chunked_feedback(all_parsed)
        return {}

    raw = claude_call(
        model=CONDITIONED_MODEL,
        temperature=0,
        max_tokens=13000,
        system=system_msg,
        messages=[{"role": "user", "content": user_msg}],
    )
    
    # Extract plain text from Anthropic object
    txt = ""
    content = getattr(raw, "content", None)
    if content and isinstance(content, list) and getattr(content[0], "type", None) == "text":
        txt = getattr(content[0], "text", "") or ""
    else:
        txt = str(raw)
    
    # if there is ```json and ``` markers, remove them
    cleaned = re.sub(r"```json\s*", "", txt)
    cleaned = re.sub(r"```", "", cleaned)
    try:
        parsed = json.loads(cleaned)
    except:
        if DEBUG:
            print("[_llm_feedback] JSON parse errorin verifier")
        parsed = {}    
    return parsed


def _merge_chunked_feedback(chunks: List[Dict[str, Any]] | None) -> Dict[str, Any]:
    if DEBUG:
        print(f"[_merge_chunked_feedback] Merging {len(chunks) if chunks else 0} chunks")
    if not chunks:
        return {}
    if len(chunks) == 1:
        return chunks[0]
    combined: Dict[str, Any] = {
        "overall": {},
        "per_statement": {},
    }
    for chunk in chunks:
        if not isinstance(chunk, dict):
            continue
        per_stmt = chunk.get("per_statement")
        if isinstance(per_stmt, dict):
            combined["per_statement"].update(per_stmt)
        overall = chunk.get("overall")
        if not isinstance(overall, dict):
            continue
        base = combined["overall"]
        if not base:
            combined["overall"] = dict(overall)
            continue
        # concatenate lists where applicable
        summary = overall.get("summary")
        if summary:
            prev = base.get("summary", "")
            base["summary"] = (prev + " " + str(summary)).strip()
        if overall.get("support_quality") and not base.get("support_quality"):
            base["support_quality"] = overall["support_quality"]
        base.setdefault("concerns", [])
        base.setdefault("suggestions", [])
        base["concerns"].extend(overall.get("concerns") or [])
        base["suggestions"].extend(overall.get("suggestions") or [])
    return combined



def _is_speculation_text(text: str) -> bool:
    return (text or "").lstrip().lower().startswith("possible speculation:")


def verify_evidence_node(state: "State") -> "State":
    if DEBUG:
        print("[verify_evidence_node] start")

    reply = state.get("llm_json", {}) or {}
    manifest = reply.get("manifest", []) or []

    answer_with_proofs = reply.get("answer_with_proofs", []) or []
    answer_plain = reply.get("answer", "") or ""

    question = (state.get("messages") or [{}])[-1].get("content", "")

    # 1) Parse statements from [(proofs, text), ...]
    statements = _parse_answer_tuples(answer_with_proofs)

    # 2) Link @TITLE -> document_index using manifest
    for s in statements:
        _link_titles_to_indices(s.get("proofs", []), manifest)

    idx2title = {int(m["index"]): str(m["title"]) for m in manifest if "index" in m and "title" in m}

    # 3) Deterministic categorical verdicts
    verified: List[Dict[str, Any]] = []
    for s in statements:
        verdict, reasons = _deterministic_verdict(s["text"], s.get("proofs", []), idx2title)
        verified.append({
            "text": s["text"],
            "proofs": s.get("proofs", []),   # includes optional char_range + title + document_index
            "verdict": verdict,
            "reasons": reasons,
            "is_speculation": _is_speculation_text(s["text"]),
        })

    # 4) LLM PASS
    overall_feedback = _llm_feedback(verified, question)
    
    support_quality = None
    retry_feedback = None
    if overall_feedback and isinstance(overall_feedback.get("overall"), dict):
        ovr = overall_feedback["overall"]
        support_quality = str(ovr.get("support_quality", "")).lower() or None
        retry_feedback = {
            "summary": ovr.get("summary") or "",
            "concerns": ovr.get("concerns") or [],
            "suggestions": ovr.get("suggestions") or []
        }
    verdict = "fail" if support_quality == "low" else "pass"
    
    # add per-statement LLM feedback in addition to deterministic verdicts
    per_stmt_feedback = overall_feedback.get("per_statement") or {}
    for i, s in enumerate(verified):
        fdbk = per_stmt_feedback.get(str(i)) or {}
        # if current verdict is 'unsupported' or 'partial' while in LLM feedback it is 'supported' or 'speculation-ok', we keep the more conservative label
        llm_label = fdbk.get("label")
        if llm_label:
            llm_label = llm_label.strip().lower()
            current_verdict = s.get("verdict")
            if current_verdict in {"unsupported", "partial"} and llm_label in {"supported", "speculation-ok"}:
                # keep current verdict & append feedback
                feedback_extra = fdbk.get("feedback")
                if feedback_extra:
                    existing_feedback = s.get("llm_feedback") or ""
                    s["llm_feedback"] = ('Deterministic feedbacl:'+ existing_feedback + "\n Additional feedback:" + feedback_extra).strip()
            else:
                # upgrade to LLM label
                s["verdict"] = llm_label
                s["llm_feedback"] = fdbk.get("feedback") or ""
        else:
            s["llm_feedback"] = ""
            s["verdict"] = s.get("verdict")  # keep deterministic verdict
            
    # recompute summary based on all verdicts when available
    def _norm_label(lbl):
        return (lbl or "").strip().lower()

    summary = {
        "supported":   sum(1 for s in verified if _norm_label(s.get("verdict")) == "supported"),
        "partial":     sum(1 for s in verified if _norm_label(s.get("verdict")) == "partial"),
        "unsupported": sum(1 for s in verified if _norm_label(s.get("verdict")) in {"unsupported", "speculation-overreach"}),
    }

    # Upgrade fail if any unsupported/speculation-overreach label
    total_statements = max(1, sum(summary.values()))
    unsupported_ratio = summary["unsupported"] / total_statements
    if verdict != "fail" and (summary["unsupported"] > 0 or unsupported_ratio >= 0.2):
        verdict = "fail"

    # recompute summary based on LLM verdicts when available
    def _norm_label(lbl):
        return (lbl or "").strip().lower()

    summary = {
        "supported":   sum(1 for s in verified if _norm_label(s.get("verdict")) == "supported"),
        "partial":     sum(1 for s in verified if _norm_label(s.get("verdict")) == "partial"),
        "unsupported": sum(1 for s in verified if _norm_label(s.get("verdict")) in {"unsupported", "speculation-overreach"}),
    }

    # Upgrade fail if any unsupported/speculation-overreach label
    total_statements = max(1, sum(summary.values()))
    unsupported_ratio = summary["unsupported"] / total_statements
    if verdict != "fail" and (summary["unsupported"] > 0 or unsupported_ratio >= 0.2):
        verdict = "fail"

    verification = {
        "verdict": verdict,
        "statements": verified,
        "summary": summary,
        "evidence": overall_feedback,
        "retry_feedback": retry_feedback,
        "answer": answer_plain,  # plain version
    }

    if DEBUG:
        print("[verify_evidence_node] verdict:", verdict, "summary:", summary)

    return {"verification": verification}
