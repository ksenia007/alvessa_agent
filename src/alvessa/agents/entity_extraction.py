from __future__ import annotations
from typing import List, Dict, Any, Optional, Tuple, Set, Iterable
import csv
import difflib
import json
import re
import requests
from pathlib import Path
from src.alvessa.clients.claude import claude_call
from src.config import DEBUG, GENE_EXTRACT_MODEL, GLINER_MODEL, GLINER_THRESHOLD, GLINER_ENTITY_LABELS, VERIFY_MODEL
from src.state import State
from flair.data import Sentence
from src.alvessa.domain.gene_class import Gene, canon_gene_key
from src.alvessa.domain.drug_class import Drug, DrugIdentifiers
from src.alvessa.domain.gene_components import GeneIdentifiers
from src.alvessa.domain.variant_class import Variant

from src.tools.aa_seq.seq_search import resolve_sequences_to_gene_records

#Amino acid recognition constants
MIN_AA_LEN = 10
AA_CHARS = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA_CHARS + AA_CHARS.lower())

# Global variables to cache the models
_gliner_model = None
_flair_model = None

# Drug lookup caches
_drug_library: Dict[str, Dict[str, Any]] = {}
_drug_term_lookup: Dict[str, Set[str]] = {}
_drug_lookup_keys: List[str] = []
_drug_max_tokens: int = 1


def _symbol_to_entrez(symbol: str) -> Optional[str]:
    """Convert an HGNC symbol to an Entrez ID via MyGene.info"""
    print(f"Resolving symbol: {symbol}")
    try:
        r = requests.get(
            "https://mygene.info/v3/query",
            params={"q": symbol, "species": "human", "fields": "entrezgene", "size": 1},
            timeout=8,
        )
        r.raise_for_status()
        hits = r.json()["hits"]
        return None if not hits else str(hits[0]["entrezgene"])
    except Exception:
        return None

# --- Model Loading Helpers ---

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

# --- New Variant Extraction Helpers ---

def _extract_rsids(text: str) -> List[str]:
    """Extracts genetic variants in the rs<id> format (e.g., rs12345)."""
    rs_id_pattern = r"rs\d+"
    found_rs_ids = re.findall(rs_id_pattern, text, re.IGNORECASE)
    # Normalize to lowercase for consistency
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
    # Define patterns for each Ensembl ID type
    gene_pattern = r"ENSG\d+"
    protein_pattern = r"ENSP\d+"
    transcript_pattern = r"ENST\d+"

    # Find all matches, ignoring case
    genes = re.findall(gene_pattern, text, re.IGNORECASE)
    proteins = re.findall(protein_pattern, text, re.IGNORECASE)
    transcripts = re.findall(transcript_pattern, text, re.IGNORECASE)

    # Normalize to uppercase and remove duplicates
    result = {
        "genes": sorted(list(set(g.upper() for g in genes))),
        "proteins": sorted(list(set(p.upper() for p in proteins))),
        "transcripts": sorted(list(set(t.upper() for t in transcripts))),
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
    
    coordinates = {
        "chrom": chrom,
        "pos": int(pos),
        "ref": ref.upper(),
        "alt": alt.upper(),
        "assembly": "GRCh38.p14"
    }
    return coordinates


def _extract_entities_with_regex(text: str) -> Dict[str, Any]:
    """
    Extracts genetic variants and Ensembl IDs using modular regex functions.
    """
    rsids = _extract_rsids(text)
    allele_pos_variants = _extract_chr_allele_pos(text)
    ensembl_ids = _extract_ensembl_ids(text)
    
    result = {
        "variants": {variant: {"rsid": variant} for variant in list(set(rsids))},
        "chr_pos_variants": {
            v: {"coordinates": get_variant_coordinates(v)} for v in list(set(allele_pos_variants))
        },
        "ensembl_genes": ensembl_ids.get("genes", []),
        "ensembl_proteins": ensembl_ids.get("proteins", []),
        "ensembl_transcripts": ensembl_ids.get("transcripts", [])
    }
    if DEBUG:
        print(f"[_extract_entities_with_regex] Found entities: {result}")
    return result


# --- Drug Extraction Helpers ---

_DRUG_TOKEN_STRIP = ".,;:!?\"'()[]{}<>"


def _canon_drug_key(name: str) -> str:
    """Normalizes a drug name for dictionary keys."""
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def _normalize_drug_phrase(value: str) -> str:
    """Creates a comparable representation for fuzzy drug matching."""
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def _split_drug_synonyms(raw_synonyms: Optional[str]) -> List[str]:
    if not raw_synonyms:
        return []
    parts = re.split(r"[;\n]", raw_synonyms)
    return [part.strip() for part in parts if part.strip()]


def _ensure_drug_library_loaded() -> None:
    """Loads the compound library once and prepares lookup tables."""
    global _drug_library, _drug_term_lookup, _drug_lookup_keys, _drug_max_tokens
    if _drug_library:
        return

    csv_path = Path(__file__).resolve().parents[3] / "local_dbs" / "compound_library_medchemexpress.csv"
    if not csv_path.exists():
        if DEBUG:
            print(f"[_ensure_drug_library_loaded] Missing compound library at {csv_path}")
        return

    try:
        with csv_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                product_name = (row.get("Product Name") or "").strip()
                if not product_name:
                    continue

                canon_key = _canon_drug_key(product_name)
                synonyms = _split_drug_synonyms(row.get("Synonyms"))
                entry = {
                    "name": product_name,
                    "catalog_number": (row.get("Catalog Number") or "").strip(),
                    "cas_number": (row.get("CAS Number") or "").strip(),
                    "synonyms": synonyms,
                }
                _drug_library[canon_key] = entry

                for term in [product_name] + synonyms:
                    term = term.strip()
                    if not term:
                        continue
                    normalized = _normalize_drug_phrase(term)
                    if not normalized:
                        continue
                    _drug_term_lookup.setdefault(normalized, set()).add(canon_key)
                    _drug_max_tokens = max(_drug_max_tokens, len(term.split()))
        _drug_lookup_keys = list(_drug_term_lookup.keys())
    except Exception as exc:
        if DEBUG:
            print(f"[_ensure_drug_library_loaded] Failed to load compound library: {exc}")
        _drug_library = {}
        _drug_term_lookup = {}
        _drug_lookup_keys = []
        _drug_max_tokens = 1


def _generate_drug_candidate_phrases(text: str) -> List[Tuple[str, str]]:
    """Returns n-gram phrases (and their normalized form) from the input text."""
    tokens: List[str] = []
    for raw in text.split():
        cleaned = raw.strip(_DRUG_TOKEN_STRIP)
        if cleaned:
            tokens.append(cleaned)

    if not tokens:
        return []

    max_window = max(1, min(_drug_max_tokens or 1, 8))  # cap to avoid combinatorial explosion
    phrases: List[Tuple[str, str]] = []
    n_tokens = len(tokens)
    for start in range(n_tokens):
        phrase_tokens: List[str] = []
        for end in range(start, min(n_tokens, start + max_window)):
            phrase_tokens.append(tokens[end])
            phrase = " ".join(phrase_tokens)
            normalized = _normalize_drug_phrase(phrase)
            if len(normalized) < 4:
                continue
            phrases.append((phrase, normalized))
    return phrases


def _register_drug_match(
    matches: Dict[str, Dict[str, Any]],
    canon_key: str,
    entry: Dict[str, Any],
    phrase: str,
    match_type: str,
) -> None:
    record = matches.setdefault(
        canon_key,
        {
            "name": entry.get("name", canon_key),
            "catalog_number": entry.get("catalog_number", ""),
            "cas_number": entry.get("cas_number", ""),
            "synonyms": entry.get("synonyms", []),
            "mentions": [],
        },
    )
    mention = {"text": phrase, "match_type": match_type}
    if mention not in record["mentions"]:
        record["mentions"].append(mention)


def _extract_drug_entities(text: str) -> Dict[str, Dict[str, Any]]:
    """Extracts small-molecule mentions by scanning against the compound library."""
    _ensure_drug_library_loaded()
    if not _drug_term_lookup:
        return {}

    matches: Dict[str, Dict[str, Any]] = {}
    phrases = _generate_drug_candidate_phrases(text)
    if not phrases:
        return matches

    for phrase, normalized in phrases:
        canon_keys = _drug_term_lookup.get(normalized)
        if canon_keys:
            for key in canon_keys:
                entry = _drug_library.get(key)
                if entry:
                    _register_drug_match(matches, key, entry, phrase, "exact")
            continue

        if len(normalized) < 6 or not _drug_lookup_keys:
            continue

        close_matches = difflib.get_close_matches(normalized, _drug_lookup_keys, n=2, cutoff=0.93)
        for close_norm in close_matches:
            for key in _drug_term_lookup.get(close_norm, set()):
                entry = _drug_library.get(key)
                if entry:
                    _register_drug_match(matches, key, entry, phrase, "approximate")

    return matches

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
    seen = set()

    def _add_sequence(seq: str) -> None:
        """Normalize, validate, deduplicate, and store a sequence."""
        # Keep only letters, normalize to upper
        letters = "".join(ch for ch in seq if ch.isalpha()).upper()
        if len(letters) < MIN_AA_LEN:
            return
        if not all(ch in AA_CHARS for ch in letters):
            return
        if letters not in seen:
            seen.add(letters)
            sequences.append(letters)

    # ------------------------------------------------------------------
    # 1) Handle FASTA-style blocks: join consecutive AA-only lines
    # ------------------------------------------------------------------
    current_block: List[str] = []

    def _flush_block():
        nonlocal current_block
        if current_block:
            _add_sequence("".join(current_block))
            current_block = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            # Blank line: end of a block
            _flush_block()
            continue

        if line.startswith(">"):
            # FASTA header: flush previous block, skip header
            _flush_block()
            continue

        # Check if this line is mostly an AA sequence
        # We allow tags/whitespace around, but require that the
        # alphabetical part is valid AAs.
        letters_only = "".join(ch for ch in line if ch.isalpha())
        if letters_only and all(ch in AA_SET for ch in letters_only):
            # Part of a FASTA-like block
            current_block.append(letters_only.upper())
        else:
            # Not a pure AA line: end of any ongoing block
            _flush_block()

    # Flush trailing block at EOF
    _flush_block()

    # ------------------------------------------------------------------
    # 2) Inline AA sequences in arbitrary prose / tags
    # ------------------------------------------------------------------
    # Regex to catch contiguous runs of AA letters, length >= MIN_AA_LEN
    pattern = re.compile(rf"[{AA_CHARS}{AA_CHARS.lower()}]{{{MIN_AA_LEN},}}")

    for match in pattern.finditer(text):
        _add_sequence(match.group(0))

    if DEBUG and sequences:
        print(f"[_extract_candidate_sequences] Found {len(sequences)} candidate sequences.")

    return sequences

def _extract_sequence_based_genes(text: str) -> Dict[str, Any]:
    """
    Integration for sequence-based gene discovery.

    This:
      1) extracts candidate protein sequences from the user text;
      2) calls resolve_sequences_to_gene_records(sequences);
      3) normalizes whatever structure comes back into:
         - a flat list of gene symbols;
         - the raw gene records (for downstream tools / debugging);
         - the list of AA sequences used for the search.
    """
    sequences = _extract_candidate_sequences(text)
    if not sequences:
        if DEBUG:
            print("[_extract_sequence_based_genes] No candidate AA sequences found.")
        return {"genes": [], "gene_records": [], "sequences": []}

    records = resolve_sequences_to_gene_records(sequences)

    genes = records.get("genes", []) or []

    if DEBUG:
        print(
            f"[_extract_sequence_based_genes] Resolved "
            f"{len(genes)} genes from {len(sequences)} sequences; "
            #f"{len(gene_records)} gene records."
        )

    # Normalized shape for the rest of the pipeline
    return {
        "genes": genes,
        #"gene_records": gene_records,
        "sequences": sequences,
    }

def _verify_genes_in_query(text: str, candidate_genes: List[str]) -> List[str]:
    """Lightweight verification: keep candidates that visibly appear in the query."""
    unique_candidates = list(dict.fromkeys(g for g in candidate_genes if g))
    if not unique_candidates:
        return []

    # Precompute normalized text variants
    lowered_text = text.lower()
    normalized_text = re.sub(r"[^a-z0-9]+", "", lowered_text)

    verified: List[str] = []
    for gene in unique_candidates:
        g_clean = gene.strip()
        if not g_clean:
            continue

        # 1) strict word boundary check
        if re.search(r"\b" + re.escape(g_clean) + r"\b", text, re.IGNORECASE):
            verified.append(gene)
            continue

        # 2) normalized substring check (drop non-alnum, lowercase)
        norm_gene = re.sub(r"[^a-z0-9]+", "", g_clean.lower())
        if norm_gene and norm_gene in normalized_text:
            verified.append(gene)
            continue

    if DEBUG:
        missing = set(unique_candidates) - set(verified)
        if missing:
            print(f"[_verify_genes_in_query] Filtered out (no match): {missing}")
    return verified


# --- Core Extraction Logic ---

def _extract_entities_with_claude(text: str) -> Dict[str, List[str]]:
    """Extracts gene symbols using Claude and returns them in a standard format."""
    system_message = (
        "Extract gene symbols and drug names from the message, but **only if they appear verbatim in the input**. "
        "Reply with a dictionary of entities, without any additional text. It needs to be a correct JSON object with two keys: drugs and genes. "
        "Example question: Is HER2 or PTEN a drug target of neratinib in breast cancer? "
        "Example answer: ```json{\"genes\": [\"HER2\", \"PTEN\"], \"drugs\": [\"neratinib\"]} ```"
        "Example **invalid** answers: \"Found these { \"genes\"...}\", \"The drugs are neratinib...\", \"I do not see any genes here.\""
    )
    try:
        response = claude_call(
            model=GENE_EXTRACT_MODEL,
            max_tokens=50,
            temperature=0,
            system=system_message,
            messages=[{"role": "user", "content": text}],
        )
        # Parse the response content dict, accounting for possible ```json fences or prefixes
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
    """Extracts gene, protein, and trait entities using Flair."""
    flair_model = _get_flair_model()
    sentence = Sentence(text)
    flair_model.predict(sentence)

    genes, proteins, traits = [], [], []
    for label in sentence.get_labels():
        entity_text = label.data_point.text.strip()
        entity_type = label.value.lower()

        if entity_type == "gene":
            genes.append(entity_text)
        elif entity_type == "protein":
            proteins.append(entity_text)
        # Traits disabled for now
        # elif entity_type in ["disease", "trait", "phenotype", "disorder", "syndrome", "condition"]:
        #     traits.append(entity_text)
            
    result = {
        "genes": list(set(genes)),
        "proteins": list(set(proteins)),
        # "traits": list(set(traits))
    }
    if DEBUG:
        print(f"[_extract_entities_with_flair] Found: {result}")
    return result


def _extract_entities_with_gliner(text: str) -> Dict[str, Any]:
    """Extracts all entities using GLiNER, returning genes, proteins, and traits."""
    gliner_model = _get_gliner_model()
    entities = gliner_model.predict_entities(text, GLINER_ENTITY_LABELS, threshold=GLINER_THRESHOLD)
    
    genes, proteins, traits = [], [], []
    for entity in entities:
        label = entity["label"].lower()
        text_entity = entity["text"].strip()

        if label == "gene":
            genes.append(text_entity)
        elif label == "protein":
            proteins.append(text_entity)
        # Traits disabled for now
        # elif label in ["disease", "trait", "phenotype", "disorder", "syndrome", "condition"]:
        #     traits.append(text_entity)

    result = {
        "genes": list(set(genes)),
        "proteins": list(set(proteins)),
        # "traits": list(set(traits))
    }
    if DEBUG:
        print(f"[_extract_entities_with_gliner] Found: {result}")
    return result

# --- Merging and Post-processing Logic ---

def _dedupe_mirna_variants(genes: Iterable[str]) -> List[str]:
    """
    Collapse miRNA-like names so that for MIR*
    we prefer specific 5p/3p (or longer) variants over bare MIR*.
    """
    genes = list(genes)

    mir_like = []
    non_mir = []
    for g in genes:
        if re.match(r"^(?:hsa-|mmu-|rno-|ptr-)?miR|^MIR|^mir|^let-7", g, re.IGNORECASE):
            mir_like.append(g)
        else:
            non_mir.append(g)

    groups = {}
    for g in mir_like:
        # core: strip 5p/3p if present at the end, and ignore trailing case diff
        core = re.sub(r"([_-]?[35][pP])$", "", g)
        core_key = core.lower()
        groups.setdefault(core_key, []).append(g)

    deduped_mir = []

    def sort_key(name: str):
        # Prefer: has 5p/3p > longer name > starts with miR > MIR > mir
        has_arm = bool(re.search(r"[_-]?[35][pP]$", name))
        return (
            1 if has_arm else 0,       # 5p/3p gets 1, bare MIR* gets 0
            len(name),                 # longer is usually more specific
            1 if name.lower().startswith("mir") else 0,
            2 if name.lower().startswith("mir") else
            3 if name.lower().startswith("mir") else 0,
        )

    for core, variants in groups.items():
        # pick the “best” representative per core
        best = sorted(variants, key=sort_key, reverse=True)[0]
        deduped_mir.append(best)

    return non_mir + deduped_mir

def _expand_and_refine_gene_names(base_genes: List[str], text: str) -> List[str]:
    """
    Expands gene symbols by finding more complete patterns in the original text.

    This function takes a list of candidate genes (e.g., "MIR") and searches the
    original text for more complete versions (e.g., "MIR_5a", "hsa-miR-21-5p").
    This is particularly useful for capturing full miRNA names.

    **miRNA Nomenclature Handling:**
    -   `[species]-miR-N`: e.g., `hsa-miR-N` (human), `mmu-miR-N` (mouse).
    -   `miR-`: Mature microRNA (functional form).
    -   `mir-`: Precursor hairpin gene (pri-/pre-miRNA).
    -   `MIR`: The gene locus encoding the precursor.
    -   `-5p` / `-3p`: Suffix indicating derivation from the 5' or 3' arm of the precursor.
    """
    expanded_genes = {g for g in base_genes if g}

    # Expand each seed by allowing common suffixes (e.g., _5p, -3p)
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

    # Direct miRNA capture, even if seeds are missing
    mir_patterns = [
        r"\b(?:hsa-|mmu-|rno-|ptr-)?miR[-_]?\d+[a-zA-Z]*?(?:[-_]\d+)?(?:[-_][35]p)?\b",  # miR-21-5p, miR520b_5P
        r"\bMIR\d+[\w-]*[35]?[pP]?\b",                                                    # MIR520b_5P
        r"\bmir[-_]?\d+[a-zA-Z]*?(?:[-_]\d+)?(?:[-_][35]p)?\b",                            # mir-21-3p
        r"\blet-7[\w-]*\b",                                                                # let-7 family
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
        cleaned_items = []
        for item in items:
            cleaned_items.extend(sub_item.strip() for sub_item in item.split(","))
        return cleaned_items

    processed_genes = clean_list(genes)
    processed_traits: List[str] = []  # traits disabled

    # Basic filters
    processed_genes = [g for g in processed_genes if "no gene" not in g.lower()]
    processed_genes = [g for g in processed_genes if 2 <= len(g) <= 30]

    # Drop things that look like amino acid sequences
    filtered_genes = []
    for g in processed_genes:
        if _looks_like_aa_sequence(g):
            if DEBUG:
                print(f"[_post_process_entities] Dropping AA-like token from genes: {g}")
            continue
        filtered_genes.append(g)
    
    # Return unique, sorted lists
    return {
        "genes": sorted(list(dict.fromkeys(filtered_genes))),
        "traits": []
    }

def _extract_entities_merged(text: str) -> Dict[str, Any]:
    """Extracts entities using all models, merges the results, and post-processes them."""
    claude_result = _extract_entities_with_claude(text)
    flair_result = _extract_entities_with_flair(text)
    gliner_result = _extract_entities_with_gliner(text)
    regex_result = _extract_entities_with_regex(text)
    drug_matches = _extract_drug_entities(text)
    seq_result = _extract_sequence_based_genes(text)

    # Merge gene symbols and Ensembl Gene IDs (ENSG)
    raw_genes = (
        claude_result.get("genes", []) +
        flair_result.get("genes", []) +
        gliner_result.get("genes", []) +
        regex_result.get("ensembl_genes", []) +
        seq_result.get("genes", [])
    )
    raw_traits = []  # traits disabled

    # Merge protein symbols and Ensembl Protein IDs (ENSP)
    raw_proteins = (
        flair_result.get("proteins", []) +
        gliner_result.get("proteins", []) +
        regex_result.get("ensembl_proteins", [])
    )

    full_drugs = claude_result.get("drugs", []) + sorted(
        list({info["name"] for info in drug_matches.values()})
    )
    # remove duplicates due to the upper / lower case differences and other minor variations
    full_drugs = [i.lower() for i in full_drugs]
    full_drugs = sorted(list(set(full_drugs)))

    # Does the microRNA expansion
    expanded_genes = _expand_and_refine_gene_names(raw_genes, text)
    processed_entities = _post_process_entities(expanded_genes, raw_traits)

    # De-duplicate genes
    unique_genes = sorted(set(processed_entities.get("genes", [])))

    # Pull sequence-derived records & sequences
    aa_sequences = seq_result.get("sequences", []) or []

    final_result = {
        # NOTE: verification is intentionally disabled to allow adding new genes
        "genes": unique_genes,
        "traits": [],
        "proteins": sorted(list(set(raw_proteins))),
        "transcripts": sorted(list(set(regex_result.get("ensembl_transcripts", [])))),
        "variants": regex_result.get("variants", {}),
        "chr_pos_variants": regex_result.get("chr_pos_variants", {}),
        "drugs": full_drugs,
        "drug_matches": drug_matches,
        # Sequence-derived info for downstream tools
        #"sequence_gene_records": sequence_gene_records,
        "aa_sequences": aa_sequences,
    }

    if DEBUG:
        print(f"[_extract_entities_merged] Final Merged Genes: {final_result['genes']}")
        print(f"[_extract_entities_merged] Final Merged Traits: {final_result['traits']}")
        print(f"[_extract_entities_merged] Final Merged Proteins: {final_result['proteins']}")
        print(f"[_extract_entities_merged] Final Merged Transcripts: {final_result['transcripts']}")
        print(f"[_extract_entities_merged] Final Merged Variants: {final_result['variants']}")
        print(f"[_extract_entities_merged] Final Merged Drugs: {final_result['drugs']}")
        print(
            f"[_extract_entities_merged] AA sequences extracted: "
            f"{len(final_result['aa_sequences'])}"
        )

    return final_result

# --- Agent Nodes ---

def entity_extraction_node(state: "State") -> "State":
    """Runs the comprehensive, merged entity extraction process."""
    user_input: str = state["messages"][-1]["content"]
    extraction_result = _extract_entities_merged(user_input)
    drug_matches = extraction_result.pop("drug_matches", {})
    if DEBUG:
        print(f"[entity_extraction_node] Extracted: {extraction_result}")
    # create gene objects in the state if any genes were found
    extraction_result['gene_entities'] = {}
    raw_genes = extraction_result.get("genes") or []
    gene_keys = [canon_gene_key(g) for g in raw_genes if g]
    for key in gene_keys:
        print(f"[entity_extraction_node] Processing gene key: {key}")
        if key in state.get("gene_entities", {}):
            continue
        g = Gene(GeneIdentifiers(symbol=key))
        g.add_tool("EntityExtraction")
        g.set_gene_ids(_symbol_to_entrez(key))
        print(_symbol_to_entrez(key))
        print(f"[entity_extraction_node] Created Gene object for: {key}")
        print(g.to_json())
        g.normalize()
        extraction_result['gene_entities'][key] = g 
        
        
    # pull in Variants
    extraction_result['variant_entities'] = {}
    raw_variants = extraction_result.get("variants") or []
    for var in raw_variants:
        if var in state.get("variant_entities", {}):
            continue
        if 'rs' not in var:
            print('Variant detected not in rsID format, skipping:', var) #TODO: add support for chr:pos:ref>alt
            continue
        v = Variant(rsID=var)
        v.add_tool("EntityExtraction")
        extraction_result['variant_entities'][var] = v

    # pull in Drugs
    extraction_result['drug_entities'] = {}
    existing_drugs = state.get("drug_entities", {})

    for canon_key, info in drug_matches.items():
        if canon_key in existing_drugs:
            continue

        identifiers = DrugIdentifiers(
            name=info.get("name", canon_key),
            canon_key=canon_key,
            synonyms=info.get("synonyms", []),
            catalog_number=info.get("catalog_number"),
            cas_number=info.get("cas_number"),
        )
        drug_obj = Drug(identifiers=identifiers)
        for mention in info.get("mentions", []):
            if isinstance(mention, dict):
                drug_obj.add_mention(mention)
        drug_obj.add_tool("EntityExtraction")
        extraction_result['drug_entities'][canon_key] = drug_obj

    # add drugs that were in "drugs" too, in case some were missing from drug_entities. These will not have the extra info
    existing_drug_names = {str(info.get("name", "")).lower() for info in drug_matches.values()}
    for info in drug_matches.values():
        for syn in info.get("synonyms", []):
            existing_drug_names.add(str(syn).lower())
        for mention in info.get("mentions", []):
            if isinstance(mention, dict):
                existing_drug_names.add(str(mention.get("text", "")).lower())

    for d in existing_drugs.values():
        if isinstance(d, Drug):
            if d.identifiers.name:
                existing_drug_names.add(str(d.identifiers.name).lower())
            for syn in d.identifiers.synonyms or []:
                existing_drug_names.add(str(syn).lower())

    claude_drugs = extraction_result.get("drugs") or []
    for drug in claude_drugs:
        if not drug:
            continue
        # ensure the drug actually appears in the user text to avoid LLM hallucinations/duplicates
        if not re.search(r"\b" + re.escape(drug) + r"\b", user_input, re.IGNORECASE):
            continue
        normalized = drug.lower()
        if normalized in existing_drug_names:
            continue
        if normalized in extraction_result['drug_entities'] or normalized in existing_drugs:
            continue

        identifiers = DrugIdentifiers(name=drug, canon_key=None, synonyms=[])
        drug_obj = Drug(identifiers=identifiers)
        drug_obj.add_mention({"text": drug, "match_type": "claude"})
        drug_obj.add_tool("EntityExtraction")
        extraction_result['drug_entities'][normalized] = drug_obj
        existing_drug_names.add(normalized)
    return extraction_result


def claude_entity_extraction_node(state: "State") -> "State":
    """Extracts gene entities using only the Claude model."""
    user_input: str = state["messages"][-1]["content"]
    claude_result = _extract_entities_with_claude(user_input)
    # TODO: handle microRNA expansion in isolated nodes
    processed_result = _post_process_entities(claude_result.get("genes", []), claude_result.get("traits", []))
    if DEBUG:
        print(f"[claude_entity_extraction_node] Extracted: {processed_result}")
    return processed_result


def gliner_entity_extraction_node(state: "State") -> "State":
    """Extracts all entity types using only the GLiNER model."""
    user_input: str = state["messages"][-1]["content"]
    gliner_result = _extract_entities_with_gliner(user_input)
    processed_result = _post_process_entities(gliner_result.get("genes", []), gliner_result.get("traits", []))
    # TODO: handle microRNA expansion in isolated nodes
    final_output = {
        "genes": processed_result["genes"],
        "traits": processed_result["traits"],
        "proteins": sorted(list(set(gliner_result.get("proteins", []))))
    }
    if DEBUG:
        print(f"[gliner_entity_extraction_node] Extracted: {final_output}")
    return final_output


def flair_entity_extraction_node(state: "State") -> "State":
    """Extracts gene, protein, and trait entities using only the Flair model."""
    user_input: str = state["messages"][-1]["content"]
    flair_result = _extract_entities_with_flair(user_input)
    processed_result = _post_process_entities(flair_result.get("genes", []), flair_result.get("traits", []))
    # TODO: handle microRNA expansion in isolated nodes
    final_output = {
        "genes": processed_result["genes"],
        "traits": processed_result["traits"],
        "proteins": sorted(list(set(flair_result.get("proteins", []))))
    }
    if DEBUG:
        print(f"[flair_entity_extraction_node] Extracted: {final_output}")
    return final_output

def gliner_claude_entity_extraction_node(state: "State") -> "State":
    """Extracts entities using a combination of GLiNER and Claude models."""
    user_input: str = state["messages"][-1]["content"]
    
    # Run both models
    gliner_result = _extract_entities_with_gliner(user_input)
    claude_result = _extract_entities_with_claude(user_input)
    
    # Merge raw results
    raw_genes = gliner_result.get("genes", []) + claude_result.get("genes", [])
    raw_traits = gliner_result.get("traits", []) + claude_result.get("traits", [])
    raw_proteins = gliner_result.get("proteins", [])  # Claude does not extract proteins
    
    # Expand and post-process
    expanded_genes = _expand_and_refine_gene_names(raw_genes, user_input)
    processed_entities = _post_process_entities(expanded_genes, raw_traits)
    
    # Compile final result
    final_output = {
        "genes": processed_entities["genes"],
        "traits": processed_entities["traits"],
        "proteins": sorted(list(set(raw_proteins)))
    }
    
    if DEBUG:
        print(f"[gliner_claude_entity_extraction_node] Extracted: {final_output}")
        
    return final_output


def gliner_flair_entity_extraction_node(state: "State") -> "State":
    """Extracts entities using a combination of GLiNER and Flair models."""
    user_input: str = state["messages"][-1]["content"]
    
    # Run both models
    gliner_result = _extract_entities_with_gliner(user_input)
    flair_result = _extract_entities_with_flair(user_input)
    
    # Merge raw results
    raw_genes = gliner_result.get("genes", []) + flair_result.get("genes", [])
    raw_traits = gliner_result.get("traits", []) + flair_result.get("traits", [])
    raw_proteins = gliner_result.get("proteins", []) + flair_result.get("proteins", [])
    
    # Expand and post-process
    expanded_genes = _expand_and_refine_gene_names(raw_genes, user_input)
    processed_entities = _post_process_entities(expanded_genes, raw_traits)
    
    # Compile final result
    final_output = {
        "genes": processed_entities["genes"],
        "traits": processed_entities["traits"],
        "proteins": sorted(list(set(raw_proteins)))
    }
    
    if DEBUG:
        print(f"[gliner_flair_entity_extraction_node] Extracted: {final_output}")
        
    return final_output


def claude_flair_entity_extraction_node(state: "State") -> "State":
    """Extracts entities using a combination of Claude and Flair models."""
    user_input: str = state["messages"][-1]["content"]
    
    # Run both models
    claude_result = _extract_entities_with_claude(user_input)
    flair_result = _extract_entities_with_flair(user_input)
    
    # Merge raw results
    raw_genes = claude_result.get("genes", []) + flair_result.get("genes", [])
    raw_traits = claude_result.get("traits", []) + flair_result.get("traits", [])
    raw_proteins = flair_result.get("proteins", []) # Claude does not extract proteins
    
    # Expand and post-process
    expanded_genes = _expand_and_refine_gene_names(raw_genes, user_input)
    processed_entities = _post_process_entities(expanded_genes, raw_traits)
    
    # Compile final result
    final_output = {
        "genes": processed_entities["genes"],
        "traits": processed_entities["traits"],
        "proteins": sorted(list(set(raw_proteins)))
    }
    
    if DEBUG:
        print(f"[claude_flair_entity_extraction_node] Extracted: {final_output}")
        
    return final_output


def variant_extraction_node(state: "State") -> "State":
    """Extracts variant and Ensembl entities using only regex."""
    user_input: str = state["messages"][-1]["content"]
    regex_result = _extract_entities_with_regex(user_input)
    # TODO: handle microRNA expansion in isolated nodes
    if DEBUG:
        print(f"[variant_extraction_node] Extracted: {regex_result}")
    return regex_result
    

# --- Graph Edges ---

def has_genes(state: "State") -> bool:
    """Edge condition helper: returns `True` if any genes were found."""
    found = bool(state.get("genes"))
    if DEBUG:
        print(f"[has_genes] Genes present? {found}")
    return found


def has_traits(state: "State") -> bool:
    """Edge condition helper: returns `True` if any traits were found."""
    found = bool(state.get("traits"))
    if DEBUG:
        print(f"[has_traits] Traits present? {found}")
    return found

def has_variants(state: "State") -> bool:
    """Edge condition helper: returns `True` if any variants were found."""
    found = bool(state.get("variants"))
    if DEBUG:
        print(f"[has_variants] Variants present? {found}")
    return found

def has_proteins(state: "State") -> bool:
    """Edge condition helper: returns `True` if any proteins were found."""
    found = bool(state.get("proteins"))
    if DEBUG:
        print(f"[has_proteins] Proteins present? {found}")
    return found

def has_transcripts(state: "State") -> bool:
    """Edge condition helper: returns `True` if any transcripts were found."""
    found = bool(state.get("transcripts"))
    if DEBUG:
        print(f"[has_transcripts] Transcripts present? {found}")
    return found


def has_drugs(state: "State") -> bool:
    """Edge condition helper: returns `True` if any small molecules were found."""
    found = bool(state.get("drugs"))
    if DEBUG:
        print(f"[has_drugs] Drugs present? {found}")
    return found
