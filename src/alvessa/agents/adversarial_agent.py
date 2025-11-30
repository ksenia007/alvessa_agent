"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-11-25


Description: 

Adversarial agent that indjects adversarial statements into the output to test verifier robustness
Assume that verification step is downstream 

Supports multiple modes of adversarial injection:

(1) Mode 1 - create contradicting statements using the given statememt but not it the provided proofs

(2) Mode 2 - create contradicting statements using the provided proofs & statement as a base

(3) Mode 3 - create overstatements that exaggerate the claims made in the statement using the provided proofs as a base

(4) Mode 4 - create subtle, numberical hallucinations by changing numbers in the statement to different numbers not supported by the proofs

(5) Mode 5 - create subtle, alphanumeric hallucinations by changing entities with numbers in them (e.g., gene names with numbers) to similar entities with different numbers not supported by the proofs


"""

from __future__ import annotations
import json
import re
import random
import unicodedata
import math
from typing import Any, Dict, List, Tuple, Iterable, Set, Optional

from src.alvessa.clients.claude import claude_call
from src.config import VERIFY_MODEL, DEBUG, CONDITIONED_MODEL, MAX_TOKENS
from src.state import State

from src.alvessa.agents.verify import _parse_answer_tuples, _link_titles_to_indices, _deterministic_verdict, _is_speculation_text, _llm_feedback

N_CHANGE = 1  # number of adversarial statements to inject
MODE = 4  # adversarial mode

def _create_contradiction_statement(statement: str, proofs: List[str]) -> str:
    """ Create a contradiction statement that contradicts the original statement without using the proofs. Use the LLM CONDITIONED_MODEL to generate the contradiction. """
    
    if DEBUG:
        print("[_create_contradiction_statement] Creating contradiction statement")
        
    content = f"""Original statement: {statement}"""
   
    prompt = f"""Your task is to generate a CONTRADICTION of the original scientific statement that still sounds plausible.

Goals:
- Reverse, negate, or meaningfully conflict with the original claim or its implications.
- The new statement should sound scientific and plausible, and could fool a casual reader.

Constraints:
1. Stay on-topic: keep the same main entities (genes, variants, pathways, diseases, tissues) and overall context.
2. Do NOT introduce completely new genes, variants, pathways, diseases, or datasets.
3. Do NOT copy the original statement verbatim; the meaning must conflict with it.
4. Keep it concise: 1–2 sentences.
5. Do NOT add explanations, reasoning steps, or meta-commentary.

Return ONLY the contradictory statement inside XML tags:

<answer>CONTRADICTORY_STATEMENT</answer>"""
  
    response = claude_call(
        model=CONDITIONED_MODEL,
        max_tokens=len(statement),
        temperature=0.0,
        system=prompt,
        messages=[{"role": "user", "content": content}],
    )
    response = getattr(response, "content", None)
    if response and isinstance(response, list) and getattr(response[0], "type", None) == "text":
        response = getattr(response[0], "text", "") or ""
    else:
        response = str(response)
    
    print("[_create_counterfactual_statement] LLM response:", response)
    
    match = re.search(r"<answer>(.*?)</answer>", response, re.DOTALL)
    if match:
        counterfactual_statement = match.group(1).strip()
        return counterfactual_statement
    else:
        if DEBUG:
            print("[_create_counterfactual_statement] No <answer> tags found in LLM response, returning original statement")
        return None 


def _format_proofs_for_prompt(proofs: List[str]) -> str:
    """Compact proofs for prompt readability."""
    if not proofs:
        return "None provided."
    joined = "\n".join(f"- {p}" for p in proofs if p)
    return joined or "None provided."


def _extract_text_response(response: Any) -> str:
    """Normalize Claude response object to plain text."""
    content = getattr(response, "content", None)
    if content and isinstance(content, list) and getattr(content[0], "type", None) == "text":
        return getattr(content[0], "text", "") or ""
    return str(response)


def _create_contradiction_with_proofs(statement: str, proofs: List[str]) -> str:
    """Create a contradiction using the provided proofs as a base."""
    if DEBUG:
        print("[_create_contradiction_with_proofs] Creating contradiction with proofs")
    prompt = f"""Your task is to generate a CONTRADICTION of the original scientific statement that is NOT supported by the provided proofs, 
but still sounds plausible and superficially consistent with them.

Goals:
- Reverse, negate, or meaningfully conflict with the original claim or its implications.
- The new statement should sound scientific and plausible, and could fool a casual reader into thinking it is based on the proofs.

Constraints:
1. Stay on-topic: keep the same main entities (genes, variants, pathways, diseases, tissues) and overall context.
2. Do NOT introduce completely new genes, variants, pathways, diseases, or datasets.
3. Do NOT copy the original statement verbatim; the meaning must conflict with it.
4. Keep it concise: 1–2 sentences.
5. Do NOT add explanations, reasoning steps, or meta-commentary.

Return ONLY the contradictory statement inside XML tags:

<answer>CONTRADICTORY_STATEMENT</answer>"""

    response = claude_call(
        model=CONDITIONED_MODEL,
        max_tokens=max(128, min(MAX_TOKENS // 4, len(statement) * 2)),
        temperature=0.0,
        system=prompt,
        messages=[{"role": "user", "content": f"""Original statement: "{statement}". Supporting proofs: {_format_proofs_for_prompt(proofs)}"""}],
    )
    response_text = _extract_text_response(response)
    if DEBUG:
        print("[_create_counterfactual_with_proofs] LLM response:", response_text)

    match = re.search(r"<answer>(.*?)</answer>", response_text, re.DOTALL)
    if match:
        return match.group(1).strip()
    if DEBUG:
        print("[_create_counterfactual_with_proofs] No <answer> tags found, returning None")
    return None


def _create_overstatement(statement: str, proofs: List[str]) -> str:
    """Create an exaggerated version of the claim using the provided proofs as a base."""
    
    if DEBUG:
        print("[_create_overstatement] Creating overstatement")
    
    prompt = """Your task is to rewrite the original scientific statement as an OVERSTATEMENT that goes clearly beyond what is warranted by the original wording and supporting proofs.

Goals:
- Make the claim stronger, more absolute, or more general than the original.
- Exaggerate causality, certainty, effect size, or scope 
  (e.g., from association to causation, from one cohort to all humans, from "may" to "is" or "clearly").
- The new statement should sound plausible and scientific, but more confident, assertive, or sweeping than the original.

Constraints:
1. Stay on-topic: keep the same main entities (genes, variants, pathways, diseases, tissues) and the same general context.
2. Do NOT introduce completely new genes, variants, pathways, diseases, or datasets.
3. You may:
   - Strengthen hedging language (e.g., "may be associated" → "is associated").
   - Broaden scope (e.g., "in specific cohort" → "across all individuals").
   - Upgrade correlation/association phrasing to causal phrasing.
4. Do NOT change numerical values unless it clearly strengthens the claim (e.g., "modest effect" → "strong effect" without touching exact numbers).
5. Keep it concise: 1–2 sentences maximum.
6. Do NOT add explanations or commentary about what you changed.
7. Do NOT produce an implausible or extreme claim; it must remain within the bounds of realistic scientific language.

Return ONLY the overstated claim inside XML tags:

<answer>OVERSTATED_STATEMENT</answer>"""
    
    content = f"""Original statement: "{statement}". Supporting proofs: {_format_proofs_for_prompt(proofs)}"""

    response = claude_call(
        model=CONDITIONED_MODEL,
        max_tokens=max(128, min(MAX_TOKENS // 4, len(statement) * 2)),
        temperature=0.0,
        system=prompt,
        messages=[{"role": "user", "content": content}],
    )
    response_text = _extract_text_response(response)
    if DEBUG:
        print("[_create_overstatement] LLM response:", response_text)

    match = re.search(r"<answer>(.*?)</answer>", response_text, re.DOTALL)
    if match:
        return match.group(1).strip()
    if DEBUG:
        print("[_create_overstatement] No <answer> tags found, returning None")
    return None


def _create_wrong_numbers(statement: str, proofs: List[str]) -> str:
    """
    Replace any number appearing in the statement with a different number to create a wrong statement. 
    """
    if DEBUG:
        print("[_create_hard_example] Creating hard adversarial example")
    prompt = f"""Your task is ONLY to alter numeric values, without changing any wording, structure, claims, or named entities.
The modified version must contradict the supporting proofs while still sounding plausible.

Rules:
1. Do NOT add explanations. Do NOT rewrite or rephrase any part of the statement.
2. Modify ONLY numbers. Every other character (words, punctuation, ordering) must remain identical.
   - Numbers appearing inside identifiers or words (e.g., BRCA1, rs12345, ENSG00000141510) should NOT be changed.
3. For every pure number in the original statement, replace it with a DIFFERENT plausible number.
4. "Plausible" means:
   - same type (integer stays integer, decimal stays decimal, percent stays percent)
   - similar scale (where possible, e.g., 0.05 → 0.03, 1000 → 1200)
5. If the original contains NO numbers, insert EXACTLY ONE plausible number that does NOT appear in the supporting proofs.
6. NEVER introduce new facts, entities, gene names, IDs, variants, datasets, or qualitative changes. Numeric substitution ONLY.
7. Keep spacing, hyphens, underscores, and formatting unchanged.
8. Do NOT alter the ordering or wording of the statement.
9. Return ONLY the modified statement inside the tags.

Return ONLY:
<answer>MODIFIED_STATEMENT</answer>"""

    content = f"""Original statement: "{statement}". Supporting proofs: {_format_proofs_for_prompt(proofs)}"""


    response = claude_call(
        model=CONDITIONED_MODEL,
        max_tokens=max(128, min(MAX_TOKENS // 4, len(statement) * 2)),
        temperature=0.0,
        system=prompt,
        messages=[{"role": "user", "content": content}],
    )
    response_text = _extract_text_response(response)
    if DEBUG:
        print("[_create_wrong_numbers] LLM response:", response_text)

    match = re.search(r"<answer>(.*?)</answer>", response_text, re.DOTALL)
    if match:
        return match.group(1).strip()
    if DEBUG:
        print("[_create_wrong_numbers] No <answer> tags found, returning None")
    return None

def _perturb_numerical_entities(statement: str, proofs: List[str]) -> str:
    """
    Replace entities w/ numbers in them (e.g., gene names with numbers) with similar entities but different numbers to create a wrong statement.
    """
    if DEBUG:
        print("[_create_hard_example] Creating hard adversarial example")
    prompt = f"""Your task is to ONLY alter identifiers that contain both letters and digits, without changing any other wording, structure, or numeric quantities.
The modified version must contradict the supporting proofs while still sounding plausible.

Definitions:
- An "alphanumeric identifier" is any token that contains at least one letter and at least one digit, possibly with separators like '-', '_', ':', or '.' inside the same token.
  Examples: BRCA1, IL6, TP53, rs12345, ENSG00000141510, MIR21-5p, R-HSA-123456.
- A "pure number" is a token that consists only of digits or standard numeric notation (e.g., 100, 3.14, 0.05, 2023, 1e-5).

Rules:
1. Do NOT add explanations. Do NOT rewrite or rephrase any part of the statement.
2. Modify ONLY the digits inside alphanumeric identifiers.
   - Do NOT change any letters.
   - Do NOT change punctuation, separators, or surrounding text.
3. Do NOT modify pure numbers (sample sizes, p-values, years, coordinates, etc.).
4. For every alphanumeric identifier in the original statement, change its numeric part to a DIFFERENT plausible pattern of digits.
5. "Plausible" means:
   - Keep the same general length or structure of the digits when possible (e.g., rs12345 → rs92347).
   - Preserve prefixes, suffixes, and formatting (e.g., keep 'rs', 'ENSG', 'R-HSA-' identical).
6. NEVER introduce new database prefixes, entity types, or change gene symbols.
7. Keep spacing, hyphens, underscores, colons, and formatting unchanged.
8. Do NOT alter the ordering or wording of the statement.
9. Return ONLY the modified statement inside the tags.

Return ONLY:
<answer>MODIFIED_STATEMENT</answer>"""

    content = f"""Original statement: "{statement}". Supporting proofs: {_format_proofs_for_prompt(proofs)}"""


    response = claude_call(
        model=CONDITIONED_MODEL,
        max_tokens=max(128, min(MAX_TOKENS // 4, len(statement) * 2)),
        temperature=0.0,
        system=prompt,
        messages=[{"role": "user", "content": content}],
    )
    print("[_create_wrong_numbers] LLM raw response:", response)
    response_text = _extract_text_response(response)
    if DEBUG:
        print("[_create_wrong_numbers] LLM response:", response_text)

    match = re.search(r"<answer>(.*?)</answer>", response_text, re.DOTALL)
    if match:
        return match.group(1).strip()
    if DEBUG:
        print("[_create_wrong_numbers] No <answer> tags found, returning None")
    return None


def _create_adversarial_statement(statement_with_proofs: List[Any], mode: int) -> Tuple[str, str]:
    """ Select which mode function to call to create adversarial statement """
    statement = statement_with_proofs[1]
    proofs = statement_with_proofs[0]
    
    if mode == 1:
        return _create_contradiction_statement(statement, proofs), "contradiction"
    if mode == 2:
        return _create_contradiction_with_proofs(statement, proofs), "contradiction_with_proofs"
    if mode == 3:
        return _create_overstatement(statement, proofs), "overstatement"
    if mode == 4:
        return _create_wrong_numbers(statement, proofs), "hallucinated_numbers"
    if mode == 5:
        return _perturb_numerical_entities(statement, proofs), "hallucinated_alphanumeric_entities"
    
    raise ValueError(f"Invalid adversarial mode: {mode}")

def adversarial_node(state: "State") -> "State":
    if DEBUG:
        print("[adversarial_node] Start adversarial eval")

    reply = state.get("llm_json", {}) or {}
    answer_with_proofs = reply.get("answer_with_proofs", []) or []

    # our goal is to modify the statements to create adversarial ones directly in the answer_with_proofs structure & record which ones were modified 
    # into state["adversarial_indices"] together with metadata about the modification: original, new, type of modification
    adversarial_indices = []
    
    # select random statements to modify, the ones with proofs that are not NONE & "possible speculation" is not in the text
    candidate_indices = [i for i, stmt in enumerate(answer_with_proofs) if stmt[0] != "NONE" and 'possible speculation' not in stmt[1].lower() and len(stmt[1].strip()) > 20]
    
    # if MODE==4, we want to select only statements that contain numbers w/ letters in it
    if MODE == 4:
        candidate_indices = [i for i in candidate_indices if re.search(r'\d', answer_with_proofs[i][1])]
    if MODE == 5:
        # we want to select only statements that contain numbers mergemd with letters in it (e.g., gene names with numbers or special characters)
        candidate_indices = [i for i in candidate_indices if re.search(r'(?i)[A-Za-z]\d|\d[A-Za-z]', answer_with_proofs[i][1])]
    
    selected_indices = random.sample(candidate_indices, min(N_CHANGE*3, len(candidate_indices)))
    if DEBUG:
        print(f"[adversarial_node] Selected indices for modification: {selected_indices}")
        print(f"[adversarial_node] Original statement: {[answer_with_proofs[i][1] for i in selected_indices]}")
        
    created_N = 0
    
    while created_N < N_CHANGE and selected_indices:
        i = selected_indices.pop(0)
        new_statement, modification_label = _create_adversarial_statement(answer_with_proofs[i], mode=MODE)
        if new_statement and new_statement != answer_with_proofs[i][1] and new_statement.lower()!='none':
            adversarial_indices.append({
                "index": i,
                "original_statement": answer_with_proofs[i][1],
                "original_proofs": answer_with_proofs[i][0],
                "new_statement": new_statement,
                "modification_type": modification_label
            })
                        
            # modify answer_with_proofs in place
            answer_with_proofs[i] = (answer_with_proofs[i][0], new_statement)
            created_N += 1
            if DEBUG:
                print(f"[adversarial_node] Created adversarial statement for index {i}: {new_statement}")
        else:
            if DEBUG:
                print(f"[adversarial_node] No adversarial statement created for index {i}, trying next candidate.")
        
    # for i in selected_indices:
    #     new_statement, modification_label = _create_adversarial_statement(answer_with_proofs[i], mode=MODE)
    #     if DEBUG:
    #         print(f"[adversarial_node] Original statement: {answer_with_proofs[i][1]}")
    #         print(f"[adversarial_node] New adversarial statement: {new_statement}")
            
    #     adversarial_indices.append({
    #         "index": i,
    #         "original_statement": answer_with_proofs[i][1],
    #         "original_proofs": answer_with_proofs[i][0],
    #         "new_statement": new_statement,
    #         "modification_type": modification_label
    #     })
                    
    #     # modify answer_with_proofs in place
    #     if new_statement:
    #         answer_with_proofs[i] = (answer_with_proofs[i][0], new_statement)
    #     else:
    #         if DEBUG:
    #             print(f"[adversarial_node] Skipping overwrite for index {i} because new_statement is None")
            
    # reassign modified structures back to state
    # state['llm_json']['answer_with_proofs'] = answer_with_proofs
    state['llm_json']['adversarial_indices'] = adversarial_indices
    
    return state
    
        
    
    
    
