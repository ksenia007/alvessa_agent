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
from src.config import VERIFY_MODEL, DEBUG, CONDITIONED_MODEL
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

    # 1) citation gate
    # if not proofs:
    #     if _WORD_RE.search(text) is None:
    #         return "unsupported", ["no-citations", "no-words"]
    #     reasons.append("no-citations")

    # 2) doc integrity
    for p in proofs:
        di = p.get("document_index")
        if di is not None and di not in idx2title:
            reasons.append(f"bad-doc-index:{di}")

    # Union proof text + titles
    union_proof = " ".join(p.get("cited_text") or "" for p in proofs)
    union_title = " ".join((p.get("title") or "") for p in proofs)
    union_all = (union_proof + " " + union_title).strip()


    # # numeric coverage
    # nums = _numbers(text)
    # missing_nums = set()
    # if nums and union_all:
    #     for n in nums:
    #         if n not in union_all:
    #             missing_nums.add(n)
    #     if missing_nums:
    #         reasons.append(f"missing-numbers:{','.join(sorted(missing_nums))}")
    
    # Verdict (categorical only)
    if "no-citations" in reasons:
        verdict = "unsupported"
    elif any(r.startswith(( "bad-doc-index")) for r in reasons):
        verdict = "partial"
    else:
        verdict = "supported"

    return verdict, reasons

# -----------------------------------
# Optional: get LLM feedback on the quality of the answer
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
        "You are a meticulous research assitant. Provide qualitative feedback on whether each sentence of the following statement is well supported "
        "by the quoted snippets (proofs). Do NOT invent facts or add citations.\n"
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
        "Where the keys in per_statement correspond to the statement indices.\n"
        "Note that the output will be parsed with json.loads(), so ensure it is valid JSON.\n"
        "Note that it is possible that citations consiting of long lists or long chunks got cutoff due to length"
        "If you see this case, in your feedback direct user to check specific field manually to verify. \n"
        "For transitional text that is non-problematic or connective, you can label as 'partial'.\n"
        "Label guidance:\n"
        "- supported: fully backed by the provided proofs.\n"
        "- partial: generally correct but missing a key entity/number or the proofs are only loosely aligned.\n"
        "- unsupported: lacks adequate proofs or conflicts with them (for non-speculation claims).\n"
        "- speculation-ok: starts with 'Possible speculation:' AND stays within what the proofs plausibly support.\n"
        "- speculation-overreach: starts with 'Possible speculation:' BUT extends beyond or contradicts the proofs.\n"
        "Avoid any numeric scoring."
    )
    user_msg = json.dumps({
        "question": question,
        "statements": pack
    }, ensure_ascii=True)

    raw = claude_call(
        model=CONDITIONED_MODEL,
        temperature=0,
        max_tokens=10000,
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

    summary = {
        "supported":   sum(1 for s in verified if s["verdict"] == "supported"),
        "partial":     sum(1 for s in verified if s["verdict"] == "partial"),
        "unsupported": sum(1 for s in verified if s["verdict"] == "unsupported"),
    }

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
    
    # add per-statement LLM feedback
    per_stmt_feedback = overall_feedback.get("per_statement") or {}
    for i, s in enumerate(verified):
        fdbk = per_stmt_feedback.get(str(i)) or {}
        s["verdict"] = fdbk.get("label")
        s["llm_feedback"] = fdbk.get("feedback")

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
