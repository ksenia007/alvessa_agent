""" Functions shared across bechmarking scripts. """
import re

MC_BENCHMARK_SYS_MSG = (
    "You are answering multiple-choice questions. Each question lists answer choices "
    "labeled [A], [B], [C], [D]. Respond with exactly one capital letter (A, B, C, or D) inside <answer> tags."
    " Do not provide any explanations or additional text."
)

MC_BENCHMARK_PROMPT_ADD = (
    "Given the question, think silently."
    "Then output the final choice as a single letter only, inside <answer></answer> tags. \n"
)


def extract_answer_letter(answer: str) -> str:
    """Extract a single-letter MC answer from <answer> tags or best-effort tokens."""
    if not answer:
        return ""
    s = str(answer)
    # Prefer any <answer> tag that contains exactly one A-D character (ignore longer content)
    tag_matches = re.findall(r"<answer>\s*([A-D])\s*</answer>", s, re.IGNORECASE)
    if tag_matches:
        # take the last valid single-letter match
        return tag_matches[-1].upper()
    # Fallback: last standalone A-D token
    tokens = re.findall(r"\\b([A-D])\\b", s, re.IGNORECASE)
    if tokens:
        return tokens[-1].upper()
    ch = s.strip()[-1:].upper()
    return ch if ch in {"A", "B", "C", "D"} else ""
