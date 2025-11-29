"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-11-25


Description: 

Adversarial agent that indjects adversarial statements into the output to test verifier robustness
Assume that verification step is downstream 

Supports multiple modes of adversarial injection:

(1) Level 1 - create counterfactual statements using the given statememt but not it the provided proofs

(2) Level 2 - create counterfactual statements using the provided proofs & statement as a base

(3) Level 3 - create overstatements that exaggerate the claims made in the statement using the provided proofs as a base

(4) Level 4 - create subtle, plausible-sounding statements that include a small unsupported change while appearing well-grounded


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

def _create_counterfactual_statement(statement: str, proofs: List[str]) -> str:
    """ Create a counterfactual statement that contradicts the original statement without using the proofs. Use the LLM CONDITIONED_MODEL to generate the counterfactual. """
    
    if DEBUG:
        print("[_create_counterfactual_statement] Creating counterfactual statement")
   
    prompt = f"""Given the following statement, create a counterfactual statement that contradicts it. 
    Reply only with the new counterfactual statement without any additional explanations. Put it inside <answer></answer> tags, and everything inside them will be used."""
    
    response = claude_call(
        model=CONDITIONED_MODEL,
        max_tokens=len(statement),
        temperature=0.0,
        system=prompt,
        messages=[{"role": "user", "content": statement}],
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


def _create_counterfactual_with_proofs(statement: str, proofs: List[str]) -> str:
    """Create a counterfactual using the provided proofs as a base."""
    if DEBUG:
        print("[_create_counterfactual_with_proofs] Creating counterfactual with proofs")
    prompt = f"""Write a counterfactual statement that contradicts the original statement provided and evidence but sounds like it is based on the evidence provided. 
    Keep it concise (1–2 sentences). Put the result inside <answer></answer> tags and do not add extra text."""

    response = claude_call(
        model=CONDITIONED_MODEL,
        max_tokens=max(128, min(MAX_TOKENS // 4, len(statement) * 2)),
        temperature=0.0,
        system=prompt,
        messages=[{"role": "user", "content": f"Generate a counterfactual with proof: Original statement: {statement} Supporting proofs:{_format_proofs_for_prompt(proofs)}"}],
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
    
    prompt = """Write an overstatement that significantly exaggerates the claims (e.g., stronger causality, broader scope, bigger effects) 
    while staying on-topic and sounding plausible and scientific. Keep it concise (1–2 sentences). Put the result inside <answer></answer> tags."""
    
    content = f"""Original statement: "{statement}"""

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
    prompt = f"""Create a new statement by changing any numbers in the original statement to different numbers, while keeping the rest of the statement the same.
    If there are no numbers in text, introduce a reasonable number that can plausibly be there. Use plausible numbers that fit the context and that are not in the proofs. 
    Return only the statement inside <answer></answer> tags."""

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


def _create_adversarial_statement(statement_with_proofs: List[Any], mode: int) -> Tuple[str, str]:
    """ Select which mode function to call to create adversarial statement """
    statement = statement_with_proofs[1]
    proofs = statement_with_proofs[0]
    
    if mode == 1:
        return _create_counterfactual_statement(statement, proofs), "counterfactual"
    if mode == 2:
        return _create_counterfactual_with_proofs(statement, proofs), "counterfactual_with_proofs"
    if mode == 3:
        return _create_overstatement(statement, proofs), "overstatement"
    if mode == 4:
        return _create_wrong_numbers(statement, proofs), "hallucinated_numbers"
    
    raise ValueError(f"Invalid adversarial mode: {mode}")

def adversarial_node(state: "State") -> "State":
    if DEBUG:
        print("[adversarial_node] Start adversarial eval")

    reply = state.get("llm_json", {}) or {}
    answer_with_proofs = reply.get("answer_with_proofs", []) or []

    # our goal is to modify the statements to create adversarial ones directly in the answer_with_proofs structure & record which ones were modified 
    # into state["adversarial_indices"] together with metadata about the modification: original, new, type of modification
    adversarial_indices = []
    modified_statements = []
    
    # select random statements to modify, the ones with proofs that are not NONE & "possible speculation" is not in the text
    candidate_indices = [i for i, stmt in enumerate(answer_with_proofs) if stmt[0] != "NONE" and 'possible speculation' not in stmt[1].lower() and len(stmt[1].strip()) > 20]
    selected_indices = random.sample(candidate_indices, min(N_CHANGE, len(candidate_indices)))
    if DEBUG:
        print(f"[adversarial_node] Selected indices for modification: {selected_indices}")
        print(f"[adversarial_node] Original statement: {[answer_with_proofs[i][1] for i in selected_indices]}")
        
    for i in selected_indices:
        new_statement, modification_label = _create_adversarial_statement(answer_with_proofs[i], mode=MODE)
        if DEBUG:
            print(f"[adversarial_node] Original statement: {answer_with_proofs[i][1]}")
            print(f"[adversarial_node] New adversarial statement: {new_statement}")
            
        adversarial_indices.append({
            "index": i,
            "original_statement": answer_with_proofs[i][1],
            "original_proofs": answer_with_proofs[i][0],
            "new_statement": new_statement,
            "modification_type": modification_label
        })
                    
        # modify answer_with_proofs in place
        if new_statement:
            answer_with_proofs[i] = (answer_with_proofs[i][0], new_statement)
        else:
            if DEBUG:
                print(f"[adversarial_node] Skipping overwrite for index {i} because new_statement is None")
            
    # reassign modified structures back to state
    # state['llm_json']['answer_with_proofs'] = answer_with_proofs
    state['llm_json']['adversarial_indices'] = adversarial_indices
    
    return state
    
        
    
    
    
