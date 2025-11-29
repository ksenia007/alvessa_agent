"""
Baseline adversarial evaluation using base Claude to flag statements as true/false.

Usage:
    python scripts/eval_adversarial.py /path/to/results.json --model claude-3-5-sonnet-20241022
    python scripts/eval_adversarial.py /path/to/results_dir --N 100
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

import anthropic
CLAUDE_MODEL = "claude-sonnet-4-5-20250929"

def call_claude(client: anthropic.Anthropic, model: str, system_msg: str, user_msg: str, max_tokens: int = 500) -> str:
    """Call Claude and return the concatenated text content."""
    try:
        resp = client.messages.create(
            model=model,
            max_tokens=max_tokens,
            temperature=0,
            system=system_msg,
            messages=[{"role": "user", "content": user_msg}],
        )
        parts = resp.content or []
        text = "".join([getattr(p, "text", "") for p in parts if getattr(p, "type", "") == "text"])
        return text.strip()
    except Exception as exc:
        print(f"[claude] error: {exc}", file=sys.stderr)
        return "ERROR"

EVAL_SYSTEM_OVERALL = (
    "You are a meticulous research assistant. Given this answer, find out if there are false statements in it. Return strict JSON with this shape:\n"
    "{\n"
    '  "verdict": "pass|fail", "concerns": ["..."]\n'
    "}. Where you set verdict to 'fail' if there are any likely false statements, else 'pass'.\n"
)

EVAL_SYSTEM = (
    "You are a meticulous research assistant. Given statements, "
    "flag each as likely true (pass) or false (fail). "
    "Return strict JSON with this shape:\n"
    "{\n"
    '  "overall": {"contains_false_statements": "true|false", "summary": "...", "concerns": ["..."], "suggestions": ["..."]},\n'
    '  "per_statement": {\n'
    '    "0": {"label": "pass|fail", "reason": "brief rationale"},\n'
    '    "1": {"label": "pass|fail", "reason": "brief rationale"}\n'
    '    "2": {"label": "pass|fail", "reason": "brief rationale"}\n ...'
    "  }\n"
    "}\n"
    "Where key is the statement index. Make sure you return same # of statements and they are aligned to original index\n"
)


def build_eval_prompt(statements: List[str]) -> str:
    payload = {"statements": [{"idx": i, "text": s} for i, s in enumerate(statements)]}
    return json.dumps(payload, ensure_ascii=True)


def parse_eval_response(raw: str) -> Dict[str, Any]:
    cleaned = raw.strip()
    if cleaned.startswith("```"):
        cleaned = cleaned.strip("`")
        cleaned = cleaned.replace("json", "", 1).strip()
    try:
        return json.loads(cleaned)
    except Exception:
        # try to take text inside ```json ... ``` in case there was also text
        start_idx = cleaned.find("```json")
        end_idx = cleaned.rfind("```")
        if start_idx != -1 and end_idx != -1 and end_idx > start_idx:
            json_part = cleaned[start_idx + len("```json") : end_idx].strip()
            try:
                return json.loads(json_part)
            except Exception:
                pass
        
    return {}


def evaluate_file(client: anthropic.Anthropic, model: str, path: Path, out_dir: Path) -> None:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        print(f"[eval] Failed to read {path}: {exc}", file=sys.stderr)
        return

    statements = data.get("verification", {}).get("statements", []) or []
    texts = [str(s.get("text", "")).strip() for s in statements if s.get("text")]
    if not texts:
        print(f"[eval] No statements found in {path}")
        return
    
    # first, glue them all together and call claude with EVAL_SYSTEM_OVERALL
    full_text = "".join(texts)
    print(f"[eval] Evaluating overall correctness of {len(texts)} statements in {path}...")
    resp_overall_one = call_claude(client, model, EVAL_SYSTEM_OVERALL, full_text, max_tokens=5000)
    print(resp_overall_one)
    print(f"[eval] Claude overall response recieved")
    overall_parsed_one = parse_eval_response(resp_overall_one)
    print("[eval] Overall parsed response (full text run):", overall_parsed_one)            
    
    # alvessa style per-statement evaluation
    user_msg = build_eval_prompt(texts)
    print(f"[eval] Evaluating {len(texts)} statements in {path}...")
    resp = call_claude(client, model, EVAL_SYSTEM, user_msg, max_tokens=10000)
    print(f"[eval] Claude response recieved")
    parsed = parse_eval_response(resp)
    per_stmt = parsed.get("per_statement") or {}
    overall = parsed.get("overall") or {}
    
    # find also adversarial_indices in data
    adv_indices = data.get("llm_json", {}).get("adversarial_indices", [])
    # but they are [{"index", "original_text", "modified_text", "modification"}], we need to know which index are there, flip to {index: Dict} structure
    adv_indices = {entry["index"]: entry for entry in adv_indices if "index" in entry}
    adv_indices_list = set(adv_indices.keys())
    
    print(f"[eval] Found {len(adv_indices_list)} adversarial statements to evaluate")
    
    results = []
    false_count = 0
    is_adversarial = False
    for i, text in enumerate(texts):
        
        if i in adv_indices_list:
            # verify that the text matches the new_statement
            adv_text = adv_indices[i]['new_statement']
            if adv_text is None:
                continue
            # strip of all extra spaces for comparison
            if text.strip() != adv_text.strip():
                print(f"[eval] Warning: statement index {i} text does not match adversarial new_statement text.", file=sys.stderr)
                print(f"[eval] Statement text: {text}", file=sys.stderr)
                print(f"[eval] Adversarial new_statement text: {adv_text}", file=sys.stderr)
                continue
            is_adversarial = True  
        else:
            is_adversarial = False

        entry = per_stmt.get(str(i), {}) or {}
        verdict = entry.get("label")
        reason = entry.get("reason")
        if verdict and verdict.lower() == "false":
            false_count += 1
        alvessa_label = statements[i].get("verdict")
        alvessa_reason = statements[i].get("llm_feedback")
        
        if is_adversarial:
            original_entry = adv_indices[i]['original_statement']
            adversarial_entry = adv_indices[i]['new_statement']
            original_proof = adv_indices[i].get('original_proofs', '')
            modification = adv_indices[i].get('modification_type', '')
        else:
            original_entry = ""
            adversarial_entry = ""
            original_proof = ""
            modification = ""
        
            
        results.append(
            {
                "text": text,
                "baseline_label": verdict,
                "baseline_reason": reason,
                "baseline_raw": entry,
                "alvessa_label": alvessa_label,
                "alvessa_reason": alvessa_reason,
                "is_adversarial": is_adversarial,
                "original_statement": original_entry,
                "adversarial_statement": adversarial_entry,
                "original_proof": original_proof,
                "modification_type": modification,
            }
        )

    baseline_overall = {
        "total": len(texts),
        "false": false_count,
        "verdict": "fail" if false_count > 0 else "pass",
        "llm_overall": overall,
    }
    alvessa_overall = {
        "verdict": data.get("verification", {}).get("verdict"),
        "summary": data.get("verification", {}).get("summary"),
    }
    
    baseline_overall_one = {
        "baseline_full": overall_parsed_one
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / path.name
    out_path.write_text(
        json.dumps(
            {
                "statements": results,
                "baseline_overall_per_statement": baseline_overall,
                "baseline_overall_llm": baseline_overall_one,
                "alvessa_overall": alvessa_overall,
                "raw_response": resp,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"[eval] Saved baseline verdicts to {out_path}")


def iter_json_files(root: Path) -> List[Path]:
    if root.is_file():
        return [root]
    return sorted([p for p in root.glob("*.json") if p.is_file()])


def main() -> int:
    parser = argparse.ArgumentParser(description="Baseline adversarial evaluation with Claude.")
    parser.add_argument("path", help="Path to a pipeline JSON file or a directory of such files.")
    parser.add_argument("--model", default=CLAUDE_MODEL, help="Claude model to use.")
    parser.add_argument("--N", type=int, default=0, help="Limit number of files to process (default all).")
    args = parser.parse_args()

    root = Path(args.path)
    files = iter_json_files(root)
    if args.N and args.N > 0:
        files = files[: args.N]
    if not files:
        print(f"[eval] No JSON files found under {root}", file=sys.stderr)
        return 1

    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    base_out = (root.parent if root.is_file() else root) / f"baseline-{ts}"

    client = anthropic.Anthropic()

    for path in files:
        evaluate_file(client, args.model, path, base_out)

    return 0


if __name__ == "__main__":
    sys.exit(main())
