from __future__ import annotations

from typing import List, Dict, Any, Optional, Tuple, Set, Iterable
import json
import re
import requests

from flair.data import Sentence

from src.alvessa.clients.claude import claude_call
from src.config import (
    DEBUG,
    GENE_EXTRACT_MODEL,
    GLINER_MODEL,
    GLINER_THRESHOLD,
    GLINER_ENTITY_LABELS,
)
from src.state import State
from src.alvessa.domain.gene_class import Gene, canon_gene_key
from src.alvessa.domain.gene_components import GeneIdentifiers
from src.alvessa.domain.variant_class import Variant

from src.tools.aa_seq.seq_search import resolve_sequences_to_gene_records
from src.tools.drug_central.drug_extraction import (
    extract_drug_matches,
    build_drug_entities,
    resolve_drug_gene_targets,
)

# ============================================================
# Amino acid recognition constants
# ============================================================

MIN_AA_LEN = 10
AA_CHARS = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA_CHARS + AA_CHARS.lower())

# Global variables to cache the models
_gliner_model = None
_flair_model = None


# ============================================================
# HGNC → Entrez helper
# ============================================================

def _symbol_to_entrez(symbol: str) -> Optional[str]:
    """Convert an HGNC symbol to an Entrez ID via MyGene.info."""
    if DEBUG:
        print(f"[_symbol_to_entrez] Resolving symbol: {symbol}")
    try:
        r = requests.get(
            "https://mygene.info/v3/query",
            params={"q": symbol, "species": "human", "fields": "entrezgene", "size": 1},
            timeout=8,
        )
        r.raise_for_status()
        hits = r.json().get("hits", [])
        return None if not hits else str(hits[0].get("entrezgene"))
    except Exception:
        return None


# ============================================================
# Model Loading Helpers
# ============================================================

def _get_gliner_model():
    """Initializes and caches the GLiNER model."""
    global _gliner_model
    if _gliner_model is None:
        try:
            from gliner import GLiNER

            _gliner_model = GLiNER.from_pretrained(GLINER_MODEL)
            if DEBUG:
                print(f"[_get_gliner_model] Loaded GLiNER model: {GLINER_MODEL}")
        except ImportError:
            raise ImportError("GLiNER is not installed. Please install it with: pip install gliner")
        except Exception as e:
            raise RuntimeError(f"Failed to load GLiNER model {GLINER_MODEL}: {e}")
    return _gliner_model


def _get_flair_model():
    """Initializes and caches the Flair model."""
    global _flair_model
    if _flair_model is None:
        try:
            from flair.nn import Classifier

            _flair_model = Classifier.load("hunflair2")
            if DEBUG:
                print(f"[_get_flair_model] Loaded FLAIR model: hunflair2")
        except ImportError:
            raise ImportError("Flair is not installed. Please install it with: pip install flair")
        except Exception as e:
            raise RuntimeError(f"Failed to load FLAIR model: {e}")
    return _flair_model


# ============================================================
# Variant Extraction Helpers
# ============================================================

def _extract_rsids(text: str) -> List[str]:
    """Extracts genetic variants in the rs<id> format (e.g., rs12345)."""
    rs_id_pattern = r"rs\d+"
    found_rs_ids = re.findall(rs_id_pattern, text, re.IGNORECASE)
    found_rs_ids = [rsid.lower() for rsid in found_rs_ids]
    if DEBUG:
        print(f"[_extract_rsids] Found: {found_rs_ids}")
    return found_rs_ids


def _extract_chr_allele_pos(text: str) -> List[str]:
    """Extracts variants in chr<num>:<pos>:<ref>><alt> format (e.g., chr7:55249071:C>T)."""
    chr_pattern = r"chr\w+:\d+:[ATCG]>[ATCG]"
    matches = re.findall(chr_pattern, text, re.IGNORECASE)
    if DEBUG:
        print(f"[_extract_chr_allele_pos] Found: {matches}")
    return matches


def _extract_ensembl_ids(text: str) -> Dict[str, List[str]]:
    """Extracts Ensembl IDs and categorizes them by type (gene, protein, transcript)."""
    gene_pattern = r"ENSG\d+"
    protein_pattern = r"ENSP\d+"
    transcript_pattern = r"ENST\d+"

    genes = re.findall(gene_pattern, text, re.IGNORECASE)
    proteins = re.findall(protein_pattern, text, re.IGNORECASE)
    transcripts = re.findall(transcript_pattern, text, re.IGNORECASE)

    result = {
        "genes": sorted({g.upper() for g in genes}),
        "proteins": sorted({p.upper() for p in proteins}),
        "transcripts": sorted({t.upper() for t in transcripts}),
    }
    if DEBUG:
        print(f"[_extract_ensembl_ids] Found: {result}")
    return result


def get_variant_coordinates(variant_string: str) -> Optional[Dict[str, Any]]:
    """
    Parses a variant string and returns its components in a structured dictionary.

    Args:
        variant_string: A string in chr<num>:<pos>:<ref>><alt> format.

    Returns:
        A dictionary with variant coordinates or None if the format is invalid.
    """
    pattern = r"chr(\w+):(\d+):([ATCG])>([ATCG])"
    match = re.match(pattern, variant_string, re.IGNORECASE)

    if not match:
        if DEBUG:
            print(f"[get_variant_coordinates] Invalid format for: {variant_string}")
        return None

    chrom, pos, ref, alt = match.groups()
    return {
        "chrom": chrom,
        "pos": int(pos),
        "ref": ref.upper(),
        "alt": alt.upper(),
        "assembly": "GRCh38.p14",
    }


def _extract_entities_with_regex(text: str) -> Dict[str, Any]:
    """Extracts genetic variants and Ensembl IDs using modular regex functions."""
    rsids = _extract_rsids(text)
    allele_pos_variants = _extract_chr_allele_pos(text)
    ensembl_ids = _extract_ensembl_ids(text)

    result = {
        "variants": {variant: {"rsid": variant} for variant in sorted(set(rsids))},
        "chr_pos_variants": {
            v: {"coordinates": get_variant_coordinates(v)}
            for v in sorted(set(allele_pos_variants))
        },
        "ensembl_genes": ensembl_ids.get("genes", []),
        "ensembl_proteins": ensembl_ids.get("proteins", []),
        "ensembl_transcripts": ensembl_ids.get("transcripts", []),
    }
    if DEBUG:
        print(f"[_extract_entities_with_regex] Found entities: {result}")
    return result


# ============================================================
# Amino-acid sequence extraction helpers
# ============================================================

def _looks_like_aa_sequence(token: str, min_len: int = MIN_AA_LEN) -> bool:
    """
    Heuristic: return True if `token` is essentially an amino acid sequence.

    - Strip non-letters.
    - Must be at least `min_len` long.
    - All letters must be in the AA alphabet (ACDEFGHIKLMNPQRSTVWY).
    """
    if not token:
        return False
    letters = re.sub(r"[^A-Za-z]", "", token)
    if len(letters) < min_len:
        return False
    return all(ch.upper() in AA_CHARS for ch in letters)


def _extract_candidate_sequences(text: str) -> List[str]:
    """
    Extract candidate amino acid sequences from free text.

    Rules:
      - AA sequences are contiguous runs of letters from the set:
          ACDEFGHIKLMNPQRSTVWY (case-insensitive).
      - Minimum length: 10 characters.
      - Handles:
          * Inline sequences (e.g., <sequence>MARLGN...</sequence>).
          * FASTA-style blocks (lines after '>' headers are joined).
      - Returns unique sequences in the order they are found.
    """
    sequences: List[str] = []
    seen: Set[str] = set()

    def _add_sequence(seq: str) -> None:
        """Normalize, validate, deduplicate, and store a sequence."""
        letters = "".join(ch for ch in seq if ch.isalpha()).upper()
        if len(letters) < MIN_AA_LEN:
            return
        if not all(ch in AA_CHARS for ch in letters):
            return
        if letters not in seen:
            seen.add(letters)
            sequences.append(letters)

    current_block: List[str] = []

    def _flush_block() -> None:
        nonlocal current_block
        if current_block:
            _add_sequence("".join(current_block))
            current_block = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            _flush_block()
            continue

        if line.startswith(">"):
            _flush_block()
            continue

        letters_only = "".join(ch for ch in line if ch.isalpha())
        if letters_only and all(ch in AA_SET for ch in letters_only):
            current_block.append(letters_only.upper())
        else:
            _flush_block()

    _flush_block()

    pattern = re.compile(rf"[{AA_CHARS}{AA_CHARS.lower()}]{{{MIN_AA_LEN},}}")
    for match in pattern.finditer(text):
        _add_sequence(match.group(0))

    if DEBUG and sequences:
        print(f"[_extract_candidate_sequences] Found {len(sequences)} candidate sequences.")

    return sequences


def _extract_sequence_based_genes(text: str) -> Dict[str, Any]:
    """Integration for sequence-based gene discovery."""
    sequences = _extract_candidate_sequences(text)
    if not sequences:
        if DEBUG:
            print("[_extract_sequence_based_genes] No candidate AA sequences found.")
        return {"genes": [], "sequences": []}

    records = resolve_sequences_to_gene_records(sequences)
    genes = records.get("genes", []) or []

    if DEBUG:
        print(
            f"[_extract_sequence_based_genes] Resolved "
            f"{len(genes)} genes from {len(sequences)} sequences."
        )

    return {
        "genes": genes,
        "sequences": sequences,
    }


def _verify_genes_in_query(text: str, candidate_genes: List[str]) -> List[str]:
    """Lightweight verification: keep candidates that visibly appear in the query."""
    unique_candidates = list(dict.fromkeys(g for g in candidate_genes if g))
    if not unique_candidates:
        return []

    lowered_text = text.lower()
    normalized_text = re.sub(r"[^a-z0-9]+", "", lowered_text)

    verified: List[str] = []
    for gene in unique_candidates:
        g_clean = gene.strip()
        if not g_clean:
            continue

        if re.search(r"\b" + re.escape(g_clean) + r"\b", text, re.IGNORECASE):
            verified.append(gene)
            continue

        norm_gene = re.sub(r"[^a-z0-9]+", "", g_clean.lower())
        if norm_gene and norm_gene in normalized_text:
            verified.append(gene)

    if DEBUG:
        missing = set(unique_candidates) - set(verified)
        if missing:
            print(f"[_verify_genes_in_query] Filtered out (no match): {missing}")
    return verified


# ============================================================
# Core Extraction Logic: Claude, Flair, GLiNER
# ============================================================

def _extract_entities_with_claude(text: str) -> Dict[str, List[str]]:
    """Extracts gene symbols and drug names using Claude and returns them in a standard format."""
    system_message = (
        "Extract gene symbols and drug names from the message, but only if they appear verbatim in the input. "
        "Reply with a dictionary of entities, without any additional text. It must be a valid JSON object "
        "with two keys: drugs and genes. "
        "Example question: Is HER2 or PTEN a drug target of neratinib in breast cancer? "
        "Example answer: ```json{\"genes\": [\"HER2\", \"PTEN\"], \"drugs\": [\"neratinib\"]} ``` "
        "Example invalid answers: 'Found these { \"genes\"...}', 'The drugs are neratinib...', "
        "'I do not see any genes here.'"
    )
    try:
        response = claude_call(
            model=GENE_EXTRACT_MODEL,
            max_tokens=50,
            temperature=0,
            system=system_message,
            messages=[{"role": "user", "content": text}],
        )
        raw_text = response.content[0].text.strip()
        raw_text = re.sub(r"^```(?:json)?\s*", "", raw_text).strip()
        raw_text = re.sub(r"```$", "", raw_text).strip()
        parsed: Any = None
        try:
            parsed = json.loads(raw_text)
        except json.JSONDecodeError:
            brace_block = re.search(r"\{[\s\S]*\}", raw_text)
            if brace_block:
                try:
                    parsed = json.loads(brace_block.group(0))
                except json.JSONDecodeError:
                    print("[_extract_entities_with_claude] JSON decode error after fence strip.")
                    print("Raw response:", raw_text)
            else:
                print("[_extract_entities_with_claude] JSON decode error, no braces found.")
                print("Raw response:", raw_text)
        if not isinstance(parsed, dict):
            parsed = {}

        genes = parsed.get("genes", []) if isinstance(parsed.get("genes", []), list) else []
        drugs = parsed.get("drugs", []) if isinstance(parsed.get("drugs", []), list) else []
    except Exception as e:
        print(f"[_extract_entities_with_claude] Error parsing: {e}")
        genes = []
        drugs = []

    result = {"genes": list(set(genes)), "traits": [], "drugs": list(set(drugs))}
    if DEBUG:
        print(f"[_extract_entities_with_claude] Found: {result}")
    return result


def _extract_entities_with_flair(text: str) -> Dict[str, List[str]]:
    """Extracts gene and protein entities using Flair."""
    flair_model = _get_flair_model()
    sentence = Sentence(text)
    flair_model.predict(sentence)

    genes: List[str] = []
    proteins: List[str] = []
    for label in sentence.get_labels():
        entity_text = label.data_point.text.strip()
        entity_type = label.value.lower()

        if entity_type == "gene":
            genes.append(entity_text)
        elif entity_type == "protein":
            proteins.append(entity_text)

    result = {
        "genes": list(set(genes)),
        "proteins": list(set(proteins)),
    }
    if DEBUG:
        print(f"[_extract_entities_with_flair] Found: {result}")
    return result


def _extract_entities_with_gliner(text: str) -> Dict[str, Any]:
    """Extracts all entities using GLiNER, returning genes and proteins."""
    gliner_model = _get_gliner_model()
    entities = gliner_model.predict_entities(text, GLINER_ENTITY_LABELS, threshold=GLINER_THRESHOLD)

    genes: List[str] = []
    proteins: List[str] = []
    for entity in entities:
        label = entity["label"].lower()
        text_entity = entity["text"].strip()

        if label == "gene":
            genes.append(text_entity)
        elif label == "protein":
            proteins.append(text_entity)

    result = {
        "genes": list(set(genes)),
        "proteins": list(set(proteins)),
    }
    if DEBUG:
        print(f"[_extract_entities_with_gliner] Found: {result}")
    return result


# ============================================================
# Gene post-processing and miRNA handling
# ============================================================

def _dedupe_mirna_variants(genes: Iterable[str]) -> List[str]:
    """
    Collapse miRNA-like names so that for MIR* we prefer specific 5p/3p
    (or longer) variants over bare MIR*, and we keep the species (hsa-, mmu-, etc.) prefix if present.
    """
    genes = list(genes)
    mir_like: List[str] = []
    non_mir: List[str] = []

    for g in genes:
        if re.match(r"^(?:hsa-|mmu-|rno-|ptr-)?miR|^MIR|^mir|^let-7", g, re.IGNORECASE):
            mir_like.append(g)
        else:
            non_mir.append(g)

    groups: Dict[str, List[str]] = {}
    for g in mir_like:
        # Extract species prefix if present; group by core (without prefix) so we dedupe prefixed/unprefixed together
        m = re.match(r"^((?:hsa|mmu|rno|ptr)[_-])?(.*)$", g, re.IGNORECASE)
        prefix = m.group(1) if m else ""
        body = m.group(2) if m else g

        core = re.sub(r"([_-]?[35][pP])$", "", body)
        core_key = core.lower()
        groups.setdefault(core_key, []).append(g)

    deduped_mir: List[str] = []

    def sort_key(name: str) -> Tuple[int, int, int]:
        has_arm = bool(re.search(r"[_-]?[35][pP]$", name))
        has_prefix = bool(re.match(r"^(?:hsa|mmu|rno|ptr)[_-]", name, re.IGNORECASE))
        # Prefer prefixed, then arm-specific, then longer names
        return (0 if has_prefix else 1, 0 if has_arm else 1, -len(name))

    for _, variants in groups.items():
        best = sorted(variants, key=sort_key)[0]
        deduped_mir.append(best)

    return non_mir + deduped_mir


def _expand_and_refine_gene_names(base_genes: List[str], text: str) -> List[str]:
    """
    Expands gene symbols by finding more complete patterns in the original text,
    and adds miRNA patterns directly from the text.
    """
    expanded_genes: Set[str] = {g for g in base_genes if g}

    for gene in base_genes:
        try:
            pattern = re.escape(gene) + r"[\w-]*"
            matches = re.findall(pattern, text, re.IGNORECASE)
            for match in matches:
                expanded_genes.add(match)
        except re.error:
            if DEBUG:
                print(f"[_expand_and_refine_gene_names] Regex error with gene: {gene}")
            continue

    mir_patterns = [
        r"\b(?:hsa-|mmu-|rno-|ptr-)?miR[-_]?\d+[a-zA-Z]*?(?:[-_]\d+)?(?:[-_][35]p)?\b",
        r"\bMIR\d+[\w-]*[35]?[pP]?\b",
        r"\bmir[-_]?\d+[a-zA-Z]*?(?:[-_]\d+)?(?:[-_][35]p)?\b",
        r"\blet-7[\w-]*\b",
    ]
    for pat in mir_patterns:
        try:
            matches = re.findall(pat, text, re.IGNORECASE)
        except re.error:
            continue
        for m in matches:
            if isinstance(m, tuple):
                for part in m:
                    if part:
                        expanded_genes.add(part)
            else:
                expanded_genes.add(m)

    expanded_genes_list = _dedupe_mirna_variants(expanded_genes)

    if DEBUG:
        print(f"[_expand_and_refine_gene_names] Expanded {len(base_genes)} base genes to {len(expanded_genes)}.")
    return expanded_genes_list


def _post_process_entities(genes: List[str], traits: List[str]) -> Dict[str, List[str]]:
    """Cleans and filters lists of genes and traits."""

    def clean_list(items: List[str]) -> List[str]:
        cleaned_items: List[str] = []
        for item in items:
            cleaned_items.extend(sub_item.strip() for sub_item in item.split(","))
        return cleaned_items

    processed_genes = clean_list(genes)
    processed_traits: List[str] = []  # traits currently unused but kept for API symmetry

    processed_genes = [g for g in processed_genes if "no gene" not in g.lower()]
    processed_genes = [g for g in processed_genes if 2 <= len(g) <= 30]

    filtered_genes: List[str] = []
    for g in processed_genes:
        if _looks_like_aa_sequence(g):
            if DEBUG:
                print(f"[_post_process_entities] Dropping AA-like token from genes: {g}")
            continue
        filtered_genes.append(g)

    return {
        "genes": sorted(dict.fromkeys(filtered_genes)),
        "traits": processed_traits,
    }


# ============================================================
# Core merged extraction (genes, variants, drugs, AA, etc.)
# ============================================================

def _extract_entities_merged(text: str) -> Dict[str, Any]:
    """Extracts entities using all models, merges the results, and post-processes them."""
    claude_result = _extract_entities_with_claude(text)
    flair_result = _extract_entities_with_flair(text)
    gliner_result = _extract_entities_with_gliner(text)
    regex_result = _extract_entities_with_regex(text)
    drug_matches = extract_drug_matches(text)
    seq_result = _extract_sequence_based_genes(text)

    raw_genes = (
        claude_result.get("genes", [])
        + flair_result.get("genes", [])
        + gliner_result.get("genes", [])
        + regex_result.get("ensembl_genes", [])
        + seq_result.get("genes", [])
    )
    raw_traits: List[str] = []

    raw_proteins = (
        flair_result.get("proteins", [])
        + gliner_result.get("proteins", [])
        + regex_result.get("ensembl_proteins", [])
    )

    # Claude + library-backed drug names (flat list)
    full_drugs = claude_result.get("drugs", []) + sorted(
        {info["name"] for info in drug_matches.values()}
    )
    full_drugs = sorted({name.lower() for name in full_drugs})

    expanded_genes = _expand_and_refine_gene_names(raw_genes, text)
    processed_entities = _post_process_entities(expanded_genes, raw_traits)

    unique_genes = sorted(set(processed_entities.get("genes", [])))
    aa_sequences = seq_result.get("sequences", []) or []

    final_result: Dict[str, Any] = {
        "genes": unique_genes,
        "traits": [],
        "proteins": sorted(set(raw_proteins)),
        "transcripts": sorted(set(regex_result.get("ensembl_transcripts", []))),
        "variants": regex_result.get("variants", {}),
        "chr_pos_variants": regex_result.get("chr_pos_variants", {}),
        "drugs": full_drugs,
        "drug_matches": drug_matches,
        "aa_sequences": aa_sequences,
    }

    if DEBUG:
        print(f"[_extract_entities_merged] Final Merged Genes: {final_result['genes']}")
        print(f"[_extract_entities_merged] Final Merged Traits: {final_result['traits']}")
        print(f"[_extract_entities_merged] Final Merged Proteins: {final_result['proteins']}")
        print(f"[_extract_entities_merged] Final Merged Transcripts: {final_result['transcripts']}")
        print(f"[_extract_entities_merged] Final Merged Variants: {final_result['variants']}")
        print(f"[_extract_entities_merged] Final Merged Drugs: {final_result['drugs']}")
        print(f"[_extract_entities_merged] AA sequences extracted: {len(final_result['aa_sequences'])}")

    return final_result


# ============================================================
# Agent Nodes
# ============================================================

def entity_extraction_node(state: "State") -> "State":
    """
    Runs the comprehensive, merged entity extraction process.
    """
    user_input: str = state["messages"][-1]["content"]
    extraction_result = _extract_entities_merged(user_input)
    drug_matches = extraction_result.pop("drug_matches", {})

    # ------------------------------------------------------------------
    # AA-sequence-derived genes as a state-level text note
    # ------------------------------------------------------------------
    seq_genes = sorted(
        {
            canon_gene_key(g)
            for g in (extraction_result.get("genes") or [])
            if g
        }
    )
    if extraction_result.get("aa_sequences") and seq_genes:
        note = (
            "Genes identified from amino-acid sequence analysis: "
            + ", ".join(seq_genes)
        )
        state.setdefault("text_notes", [])
        if note not in state["text_notes"]:
            state["text_notes"].append(note)

    if DEBUG:
        print(f"[entity_extraction_node] Extracted (pre-entities): {extraction_result}")

    # --------------------------------------------------------
    # Initialize / preserve entity dicts from state
    # --------------------------------------------------------
    state_gene_entities: Dict[str, Gene] = state.get("gene_entities") or {}
    state_variant_entities: Dict[str, Variant] = state.get("variant_entities") or {}
    state_drug_entities: Dict[str, Any] = state.get("drug_entities") or {}

    gene_entities: Dict[str, Gene] = dict(state_gene_entities)
    variant_entities: Dict[str, Variant] = dict(state_variant_entities)

    # --------------------------------------------------------
    # Gene entities from text/sequence
    # --------------------------------------------------------
    raw_genes = extraction_result.get("genes") or []
    for g_symbol in raw_genes:
        sym = (g_symbol or "").strip()
        if not sym:
            continue
        key = canon_gene_key(sym)
        if not key:
            continue
        if key in gene_entities:
            continue

        entrez = _symbol_to_entrez(sym)
        gid = GeneIdentifiers(symbol=sym, entrez_id=entrez)
        g = Gene(gid)
        g.add_tool("EntityExtraction")
        g.normalize()
        gene_entities[key] = g

        if DEBUG:
            print(
                f"[entity_extraction_node] Created Gene object from text: "
                f"{sym} (key={key}, entrez={entrez})"
            )

    # --------------------------------------------------------
    # Variant entities
    # --------------------------------------------------------
    extraction_result["variant_entities"] = variant_entities
    raw_variants = list((extraction_result.get("variants") or {}).keys())
    for var in raw_variants:
        if var in variant_entities:
            continue
        if "rs" not in var.lower():
            if DEBUG:
                print("[entity_extraction_node] Variant not in rsID format, skipping:", var)
            continue
        v = Variant(rsID=var)
        v.add_tool("EntityExtraction")
        variant_entities[var] = v

    # --------------------------------------------------------
    # Drug entities (MedChemExpress + Claude + ID-based)
    # --------------------------------------------------------
    claude_drugs = extraction_result.get("drugs") or []

    new_drug_entities = build_drug_entities(
        user_input=user_input,
        state=state,
        claude_drugs=claude_drugs,
        drug_matches=drug_matches,
    )

    merged_drug_entities: Dict[str, Any] = dict(state_drug_entities)
    merged_drug_entities.update(new_drug_entities)
    state["drug_entities"] = merged_drug_entities
    extraction_result["drug_entities"] = merged_drug_entities

    # --------------------------------------------------------
    # Drug → gene targets
    # --------------------------------------------------------
    dc_target_map = resolve_drug_gene_targets(list(new_drug_entities.values()))

    drug_targets: Dict[str, List[Dict[str, Any]]] = {}
    for struct_id, targets in (dc_target_map or {}).items():
        if targets:
            drug_targets[str(struct_id)] = list(targets)

    extraction_result["drug_targets"] = drug_targets

    target_gene_info: Dict[str, Dict[str, Any]] = {}
    for targets in drug_targets.values():
        for t in targets:
            sym_raw = (t.get("symbol") or "").strip()
            if not sym_raw:
                continue
            key = canon_gene_key(sym_raw)
            if not key:
                continue
            existing = target_gene_info.get(key)
            if existing is None:
                target_gene_info[key] = dict(t)
            else:
                for field in ("entrez_id", "uniprot_id", "ensembl_id"):
                    if not existing.get(field) and t.get(field):
                        existing[field] = t[field]

    for key, info in target_gene_info.items():
        if key in gene_entities:
            continue
        entrez = info.get("entrez_id")
        symbol_for_display = info.get("symbol") or key
        gid = GeneIdentifiers(symbol=symbol_for_display, entrez_id=entrez)
        g = Gene(gid)
        g.add_tool("DrugTargetResolution")
        g.normalize()
        gene_entities[key] = g

        if DEBUG:
            print(
                f"[entity_extraction_node] Created Gene object from drug targets: "
                f"{symbol_for_display} (key={key}, entrez={entrez})"
            )

    # --------------------------------------------------------
    # Sync top-level drug list with state
    # --------------------------------------------------------
    extraction_result["drugs"] = list(state.get("drugs") or [])

    # --------------------------------------------------------
    # Final gene list
    # --------------------------------------------------------
    all_gene_keys: Set[str] = set(gene_entities.keys())
    extraction_result["genes"] = sorted(all_gene_keys)
    extraction_result["gene_entities"] = gene_entities
    extraction_result["variant_entities"] = variant_entities

    if DEBUG:
        print(f"[entity_extraction_node] Final genes after target expansion: {extraction_result['genes']}")

    return extraction_result


def claude_entity_extraction_node(state: "State") -> "State":
    """Extracts gene entities using only the Claude model."""
    user_input: str = state["messages"][-1]["content"]
    claude_result = _extract_entities_with_claude(user_input)
    processed_result = _post_process_entities(
        claude_result.get("genes", []),
        claude_result.get("traits", []),
    )
    if DEBUG:
        print(f"[claude_entity_extraction_node] Extracted: {processed_result}")
    return processed_result


def gliner_entity_extraction_node(state: "State") -> "State":
    """Extracts all entity types using only the GLiNER model."""
    user_input: str = state["messages"][-1]["content"]
    gliner_result = _extract_entities_with_gliner(user_input)
    processed_result = _post_process_entities(
        gliner_result.get("genes", []),
        gliner_result.get("traits", []),
    )
    final_output = {
        "genes": processed_result["genes"],
        "traits": processed_result["traits"],
        "proteins": sorted(set(gliner_result.get("proteins", []))),
    }
    if DEBUG:
        print(f"[gliner_entity_extraction_node] Extracted: {final_output}")
    return final_output


def flair_entity_extraction_node(state: "State") -> "State":
    """Extracts gene, protein, and trait entities using only the Flair model."""
    user_input: str = state["messages"][-1]["content"]
    flair_result = _extract_entities_with_flair(user_input)
    processed_result = _post_process_entities(
        flair_result.get("genes", []),
        flair_result.get("traits", []),
    )
    final_output = {
        "genes": processed_result["genes"],
        "traits": processed_result["traits"],
        "proteins": sorted(set(flair_result.get("proteins", []))),
    }
    if DEBUG:
        print(f"[flair_entity_extraction_node] Extracted: {final_output}")
    return final_output


def gliner_claude_entity_extraction_node(state: "State") -> "State":
    """Extracts entities using a combination of GLiNER and Claude models."""
    user_input: str = state["messages"][-1]["content"]

    gliner_result = _extract_entities_with_gliner(user_input)
    claude_result = _extract_entities_with_claude(user_input)

    raw_genes = gliner_result.get("genes", []) + claude_result.get("genes", [])
    raw_traits = gliner_result.get("traits", []) + claude_result.get("traits", [])
    raw_proteins = gliner_result.get("proteins", [])

    expanded_genes = _expand_and_refine_gene_names(raw_genes, user_input)
    processed_entities = _post_process_entities(expanded_genes, raw_traits)

    final_output = {
        "genes": processed_entities["genes"],
        "traits": processed_entities["traits"],
        "proteins": sorted(set(raw_proteins)),
    }

    if DEBUG:
        print(f"[gliner_claude_entity_extraction_node] Extracted: {final_output}")

    return final_output


def gliner_flair_entity_extraction_node(state: "State") -> "State":
    """Extracts entities using a combination of GLiNER and Flair models."""
    user_input: str = state["messages"][-1]["content"]

    gliner_result = _extract_entities_with_gliner(user_input)
    flair_result = _extract_entities_with_flair(user_input)

    raw_genes = gliner_result.get("genes", []) + flair_result.get("genes", [])
    raw_traits = gliner_result.get("traits", []) + flair_result.get("traits", [])
    raw_proteins = gliner_result.get("proteins", []) + flair_result.get("proteins", [])

    expanded_genes = _expand_and_refine_gene_names(raw_genes, user_input)
    processed_entities = _post_process_entities(expanded_genes, raw_traits)

    final_output = {
        "genes": processed_entities["genes"],
        "traits": processed_entities["traits"],
        "proteins": sorted(set(raw_proteins)),
    }

    if DEBUG:
        print(f"[gliner_flair_entity_extraction_node] Extracted: {final_output}")

    return final_output


def claude_flair_entity_extraction_node(state: "State") -> "State":
    """Extracts entities using a combination of Claude and Flair models."""
    user_input: str = state["messages"][-1]["content"]

    claude_result = _extract_entities_with_claude(user_input)
    flair_result = _extract_entities_with_flair(user_input)

    raw_genes = claude_result.get("genes", []) + flair_result.get("genes", [])
    raw_traits = claude_result.get("traits", []) + flair_result.get("traits", [])
    raw_proteins = flair_result.get("proteins", [])

    expanded_genes = _expand_and_refine_gene_names(raw_genes, user_input)
    processed_entities = _post_process_entities(expanded_genes, raw_traits)

    final_output = {
        "genes": processed_entities["genes"],
        "traits": processed_entities["traits"],
        "proteins": sorted(set(raw_proteins)),
    }

    if DEBUG:
        print(f"[claude_flair_entity_extraction_node] Extracted: {final_output}")

    return final_output


def variant_extraction_node(state: "State") -> "State":
    """Extracts variant and Ensembl entities using only regex."""
    user_input: str = state["messages"][-1]["content"]
    regex_result = _extract_entities_with_regex(user_input)
    if DEBUG:
        print(f"[variant_extraction_node] Extracted: {regex_result}")
    return regex_result


# ============================================================
# Graph Edge Helpers
# ============================================================

def has_genes(state: "State") -> bool:
    """Edge condition helper: returns True if any genes were found."""
    found = bool(state.get("genes"))
    if DEBUG:
        print(f"[has_genes] Genes present? {found}")
    return found


def has_traits(state: "State") -> bool:
    """Edge condition helper: returns True if any traits were found."""
    found = bool(state.get("traits"))
    if DEBUG:
        print(f"[has_traits] Traits present? {found}")
    return found


def has_variants(state: "State") -> bool:
    """Edge condition helper: returns True if any variants were found."""
    found = bool(state.get("variants"))
    if DEBUG:
        print(f"[has_variants] Variants present? {found}")
    return found


def has_proteins(state: "State") -> bool:
    """Edge condition helper: returns True if any proteins were found."""
    found = bool(state.get("proteins"))
    if DEBUG:
        print(f"[has_proteins] Proteins present? {found}")
    return found


def has_transcripts(state: "State") -> bool:
    """Edge condition helper: returns True if any transcripts were found."""
    found = bool(state.get("transcripts"))
    if DEBUG:
        print(f"[has_transcripts] Transcripts present? {found}")
    return found


def has_drugs(state: "State") -> bool:
    """Edge condition helper: returns True if any small molecules were found."""
    found = bool(state.get("drugs"))
    if DEBUG:
        print(f"[has_drugs] Drugs present? {found}")
    return found
