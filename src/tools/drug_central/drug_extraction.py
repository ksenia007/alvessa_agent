# src/tools/drug_central/drug_extraction.py
from __future__ import annotations

from typing import Dict, Any, List, Set, Tuple, Optional
import csv
import difflib
import re
import sys
from pathlib import Path

# --- Package/Repo Roots ---
PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# --- Imports ---
from src.config import DEBUG
from src.alvessa.domain.drug_class import Drug, DrugIdentifiers
from src.state import State

import src.tools.drug_central as drug_central  # loads constants from __init__.py

# Import package-level paths and metadata
from src.tools.drug_central import LOCAL_DBS_DIR, DRUG_CENTRAL_DB_PATH


# --------------------------------------------------------------------
# Internal caches and constants
# --------------------------------------------------------------------

_DRUG_TOKEN_STRIP = ".,;:!?\"'()[]{}<>"

# Internal library:
#   internal_key -> { "name", "catalog_number", "cas_number", "synonyms" }
_drug_library: Dict[str, Dict[str, Any]] = {}

# Normalized phrase -> set of internal_keys
_drug_term_lookup: Dict[str, Set[str]] = {}

# All normalized phrases (for difflib)
_drug_lookup_keys: List[str] = []

# Longest token length (in words) across all library names/synonyms
_drug_max_tokens: int = 1


# --------------------------------------------------------------------
# Helper functions for library and matching
# --------------------------------------------------------------------

def _canon_drug_key(name: str) -> str:
    """
    Internal helper: normalize a name into an ASCII, lowercase, alnum-only key.

    NOTE:
    This is *not* stored in DrugIdentifiers. It is only used as an
    in-memory key for the MedChemExpress library.
    """
    return re.sub(r"[^a-z0-9]+", "", (name or "").lower())


def _normalize_drug_phrase(value: str) -> str:
    """Creates a comparable representation for fuzzy drug matching."""
    return re.sub(r"[^a-z0-9]+", "", (value or "").lower())


def _split_drug_synonyms(raw_synonyms: str | None) -> List[str]:
    if not raw_synonyms:
        return []
    parts = re.split(r"[;\n]", raw_synonyms)
    return [part.strip() for part in parts if part.strip()]


def _ensure_drug_library_loaded() -> None:
    """
    Load the MedChemExpress compound library (CSV) once and build lookup tables.

    Populates:
      - _drug_library: internal_key -> entry
      - _drug_term_lookup: normalized_phrase -> {internal_key, ...}
      - _drug_lookup_keys: list of all normalized phrases (for difflib)
      - _drug_max_tokens: max phrase length in tokens (for n-gram search)
    """
    global _drug_library, _drug_term_lookup, _drug_lookup_keys, _drug_max_tokens

    if _drug_library:
        return

    csv_path = LOCAL_DBS_DIR / "compound_library_medchemexpress.csv"
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

                internal_key = _canon_drug_key(product_name)
                synonyms = _split_drug_synonyms(row.get("Synonyms"))

                entry = {
                    "name": product_name,
                    "catalog_number": (row.get("Catalog Number") or "").strip(),
                    "cas_number": (row.get("CAS Number") or "").strip(),
                    "synonyms": synonyms,
                }
                _drug_library[internal_key] = entry

                for term in [product_name] + synonyms:
                    term = (term or "").strip()
                    if not term:
                        continue
                    normalized = _normalize_drug_phrase(term)
                    if not normalized:
                        continue
                    _drug_term_lookup.setdefault(normalized, set()).add(internal_key)
                    _drug_max_tokens = max(_drug_max_tokens, len(term.split()))

        _drug_lookup_keys = list(_drug_term_lookup.keys())
        if DEBUG:
            print(
                f"[_ensure_drug_library_loaded] Loaded {len(_drug_library)} compounds from {csv_path}; "
                f"{len(_drug_term_lookup)} normalized terms; max_tokens={_drug_max_tokens}"
            )
    except Exception as exc:
        if DEBUG:
            print(f"[_ensure_drug_library_loaded] Failed to load compound library: {exc}")
        _drug_library = {}
        _drug_term_lookup = {}
        _drug_lookup_keys = []
        _drug_max_tokens = 1


def _generate_drug_candidate_phrases(text: str) -> List[Tuple[str, str]]:
    """
    Generate candidate n-gram phrases from free text.

    Returns
    -------
    List[(phrase, normalized_phrase)]
    """
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
    internal_key: str,
    entry: Dict[str, Any],
    phrase: str,
    match_type: str,
) -> None:
    """
    Register a hit from the library into the matches dict.

    Parameters
    ----------
    matches : dict
        Accumulator mapping internal_key -> info dict.
    internal_key : str
        MedChemExpress internal name key (normalized).
    entry : dict
        Library entry with "name", "catalog_number", "cas_number", "synonyms".
    phrase : str
        Phrase from the input text that matched.
    match_type : str
        "exact" or "approximate".
    """
    record = matches.setdefault(
        internal_key,
        {
            "name": entry.get("name", internal_key),
            "catalog_number": entry.get("catalog_number", ""),
            "cas_number": entry.get("cas_number", ""),
            "synonyms": entry.get("synonyms", []),
            "mentions": [],
        },
    )
    mention = {"text": phrase, "match_type": match_type, "source": "medchemexpress"}
    if mention not in record["mentions"]:
        record["mentions"].append(mention)


# --------------------------------------------------------------------
# Public API: text -> library-backed drug_matches
# --------------------------------------------------------------------

def extract_drug_matches(text: str) -> Dict[str, Dict[str, Any]]:
    """
    Extract small-molecule mentions by scanning against the MedChemExpress library.

    Parameters
    ----------
    text : str
        Free text input.

    Returns
    -------
    dict
        Mapping:
          internal_key -> {
              "name": str,
              "catalog_number": str,
              "cas_number": str,
              "synonyms": List[str],
              "mentions": List[{"text": str, "match_type": "exact"|"approximate", "source": "medchemexpress"}]
          }
    """
    _ensure_drug_library_loaded()
    if not _drug_term_lookup:
        return {}

    matches: Dict[str, Dict[str, Any]] = {}
    phrases = _generate_drug_candidate_phrases(text)
    if not phrases:
        return matches

    for phrase, normalized in phrases:
        # Direct lookup (exact normalized phrase)
        internal_keys = _drug_term_lookup.get(normalized)
        if internal_keys:
            for key in internal_keys:
                entry = _drug_library.get(key)
                if entry:
                    _register_drug_match(matches, key, entry, phrase, "exact")
            continue

        # Fuzzy match
        if len(normalized) < 6 or not _drug_lookup_keys:
            continue

        close_matches = difflib.get_close_matches(normalized, _drug_lookup_keys, n=2, cutoff=0.93)
        for close_norm in close_matches:
            for key in _drug_term_lookup.get(close_norm, set()):
                entry = _drug_library.get(key)
                if entry:
                    _register_drug_match(matches, key, entry, phrase, "approximate")

    if DEBUG:
        print(f"[extract_drug_matches] Found {len(matches)} MedChemExpress-backed drugs.")
    return matches


# --------------------------------------------------------------------
# Public API: build Drug objects for the state
# --------------------------------------------------------------------

def build_drug_entities(
    user_input: str,
    state: State,
    claude_drugs: List[str],
    drug_matches: Dict[str, Dict[str, Any]],
) -> Dict[str, Drug]:
    """
    Turn raw `drug_matches` + Claude's flat drug list into Drug objects.

    Strategy
    --------
    0) Normalize any existing state["drug_entities"] entries so they are
       actual Drug objects, not plain dicts.

    1) Start from MedChemExpress-backed matches (`drug_matches`).
       - Each internal_key becomes one Drug instance.
       - Identifiers initially include: name, CAS, catalog_number, synonyms.
       - No canonical key is stored in the identifiers.

    2) Build a case-insensitive set of known drug names/synonyms
       from:
         - library-backed matches,
         - existing state drug entities,
         - newly created Drug objects in this function.

    3) For each Claude-only drug name:
       - Keep it only if it literally appears in the user text.
       - Skip if already known (by name or synonym).
       - Create a minimal Drug with name only.

    The keys in the returned dict are:
      - MedChemExpress internal keys for library-backed drugs.
      - lowercased name strings for Claude-only drugs.
    """

    # --------------------------------------------------------
    # 0) Normalize existing state["drug_entities"] to Drug objs
    # --------------------------------------------------------
    existing_raw = state.get("drug_entities", {}) or {}
    normalized_existing: Dict[str, Drug] = {}

    for key, val in existing_raw.items():
        if isinstance(val, Drug):
            normalized_existing[key] = val
            continue

        if isinstance(val, dict):
            # Try to read identifiers from nested "identifiers" dict if present
            id_dict = val.get("identifiers") if isinstance(val.get("identifiers"), dict) else val

            name = id_dict.get("name")
            chembl_id = id_dict.get("chembl_id")
            drugcentral_id = id_dict.get("drugcentral_id")
            cas_number = id_dict.get("cas_number")
            catalog_number = id_dict.get("catalog_number")
            synonyms = id_dict.get("synonyms") or []

            if not isinstance(synonyms, list):
                synonyms = [str(synonyms)]

            identifiers = DrugIdentifiers(
                name=name,
                chembl_id=chembl_id,
                drugcentral_id=drugcentral_id,
                cas_number=cas_number,
                catalog_number=catalog_number,
                synonyms=synonyms,
            )
            d_obj = Drug(identifiers=identifiers)

            # Optionally restore mentions if they existed
            mentions = val.get("mentions") if isinstance(val.get("mentions"), list) else []
            for m in mentions:
                if isinstance(m, dict):
                    d_obj.add_mention(m)

            # Optionally restore tools if they existed
            tools = val.get("tools") if isinstance(val.get("tools"), list) else []
            for tool_name in tools:
                try:
                    d_obj.add_tool(str(tool_name))
                except Exception:
                    pass

            normalized_existing[key] = d_obj
        else:
            # Fallback: opaque object, wrap as Drug with name=str(val)
            identifiers = DrugIdentifiers(
                name=str(val),
                chembl_id=None,
                drugcentral_id=None,
                cas_number=None,
                catalog_number=None,
                synonyms=[],
            )
            normalized_existing[key] = Drug(identifiers=identifiers)

    # Write normalized mapping back into state
    state["drug_entities"] = normalized_existing  # type: ignore[index]
    existing_drugs: Dict[str, Drug] = normalized_existing

    # ----------------------------
    # 1) Library-backed drugs
    # ----------------------------
    drug_entities: Dict[str, Drug] = {}

    for internal_key, info in (drug_matches or {}).items():
        if internal_key in existing_drugs:
            # State already has a Drug object keyed by this internal key
            continue

        identifiers = DrugIdentifiers(
            name=info.get("name", internal_key),
            chembl_id=None,
            drugcentral_id=None,
            cas_number=info.get("cas_number"),
            catalog_number=info.get("catalog_number"),
            synonyms=info.get("synonyms", []),
        )
        drug_obj = Drug(identifiers=identifiers)

        for mention in info.get("mentions", []):
            if isinstance(mention, dict):
                drug_obj.add_mention(mention)

        drug_obj.add_tool("EntityExtraction")
        drug_entities[internal_key] = drug_obj

    # ----------------------------
    # 2) Build a set of known names
    # ----------------------------
    existing_drug_names: Set[str] = set()

    # From library-backed matches
    for info in (drug_matches or {}).values():
        if info.get("name"):
            existing_drug_names.add(str(info["name"]).lower())
        for syn in info.get("synonyms", []):
            existing_drug_names.add(str(syn).lower())
        for mention in info.get("mentions", []):
            if isinstance(mention, dict):
                existing_drug_names.add(str(mention.get("text", "")).lower())

    # From pre-existing state drug entities (now guaranteed Drug)
    for d in existing_drugs.values():
        ids = d.identifiers
        if ids.name:
            existing_drug_names.add(str(ids.name).lower())
        for syn in ids.synonyms or []:
            existing_drug_names.add(str(syn).lower())

    # From newly created library-backed Drug objects in this call
    for d in drug_entities.values():
        ids = d.identifiers
        if ids.name:
            existing_drug_names.add(str(ids.name).lower())
        for syn in ids.synonyms or []:
            existing_drug_names.add(str(syn).lower())

    # ----------------------------
    # 3) Claude-only drugs
    # ----------------------------
    for drug_name in claude_drugs or []:
        if not drug_name:
            continue

        # Ensure the token actually appears in the text to reduce hallucinations
        if not re.search(r"\b" + re.escape(drug_name) + r"\b", user_input, re.IGNORECASE):
            continue

        normalized = drug_name.lower()
        if normalized in existing_drug_names:
            continue

        if normalized in drug_entities or normalized in existing_drugs:
            continue

        identifiers = DrugIdentifiers(
            name=drug_name,
            chembl_id=None,
            drugcentral_id=None,
            cas_number=None,
            catalog_number=None,
            synonyms=[],
        )
        drug_obj = Drug(identifiers=identifiers)
        drug_obj.add_mention({"text": drug_name, "match_type": "claude", "source": "claude"})
        drug_obj.add_tool("EntityExtraction")

        drug_entities[normalized] = drug_obj
        existing_drug_names.add(normalized)

    if DEBUG:
        print(f"[build_drug_entities] Built {len(drug_entities)} Drug entities.")

    return drug_entities


# --------------------------------------------------------------------
# API for drug -> gene target resolution
# --------------------------------------------------------------------

def _get_targets_from_chembl(ids: DrugIdentifiers) -> List[Dict[str, Any]]:
    """
    Look up gene targets from ChEMBL for a given drug.

    Uses:
      - src.tools.chembl.utils.normalize_chembl_drug_identifier
      - src.tools.chembl.utils.get_chembl_drug_targets

    Behavior
    --------
    - Tries to resolve a canonical ChEMBL compound (chembl_id) from:
        chembl_id (if already set), name, CAS, synonyms.
    - If resolution succeeds, it *enriches* the DrugIdentifiers in-place:
        - fills chembl_id if missing,
        - fills cas_number if missing,
        - extends synonyms with ChEMBL synonyms (dedup is handled upstream).
    - Returns a list of target records of the form:
        {
          "symbol": "EGFR",
          "entrez_id": None,
          "uniprot_id": "P00533",
          "ensembl_id": None,
          "source": "chembl",
        }
    """
    try:
        from src.tools.chembl.utils import (
            normalize_chembl_drug_identifier,
            get_chembl_drug_targets,
        )
    except Exception as exc:
        if DEBUG:
            print(f"[_get_targets_from_chembl] ChemBL utils unavailable: {exc}")
        return []

    # Build candidate identifiers in priority order
    candidates: List[str] = []

    if ids.chembl_id:
        candidates.append(str(ids.chembl_id).strip())

    if ids.name:
        candidates.append(ids.name.strip())

    if ids.cas_number:
        candidates.append(str(ids.cas_number).strip())

    for syn in ids.synonyms or []:
        s = (syn or "").strip()
        if s:
            candidates.append(s)

    # Deduplicate candidates while preserving order
    seen_cand: Set[str] = set()
    ordered_candidates: List[str] = []
    for c in candidates:
        key = c.lower()
        if key in seen_cand:
            continue
        seen_cand.add(key)
        ordered_candidates.append(c)

    chembl_record: Dict[str, Any] | None = None

    for token in ordered_candidates:
        rec = normalize_chembl_drug_identifier(token)
        if not rec:
            continue
        chembl_record = rec
        break

    if not chembl_record:
        if DEBUG:
            print(f"[_get_targets_from_chembl] No ChEMBL record resolved for {ids.name!r}")
        return []

    chembl_id = chembl_record.get("chembl_id")
    if not chembl_id:
        if DEBUG:
            print(f"[_get_targets_from_chembl] ChEMBL record without chembl_id for {ids.name!r}")
        return []

    # Enrich DrugIdentifiers in-place (use ChemBL to fill attributes)
    if not ids.chembl_id:
        ids.chembl_id = str(chembl_id).strip()

    pref_name = chembl_record.get("pref_name") or ""
    if not ids.name and pref_name:
        ids.name = pref_name

    cas_from_chembl = chembl_record.get("cas_number")
    if not ids.cas_number and cas_from_chembl:
        ids.cas_number = str(cas_from_chembl).strip()

    syns_from_chembl = chembl_record.get("synonyms") or []
    if syns_from_chembl:
        # Simple extend; dedup will be handled by DrugIdentifiers.normalize()
        ids.synonyms.extend([s for s in syns_from_chembl if s])

    # Now fetch targets for this ChEMBL compound
    try:
        targets_raw = get_chembl_drug_targets(str(chembl_id), only_human=True)
    except Exception as exc:
        if DEBUG:
            print(f"[_get_targets_from_chembl] Error querying ChemBL targets for {chembl_id}: {exc}")
        return []

    by_key: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for t in targets_raw or []:
        gene_symbol = (t.get("gene_symbol") or "").strip()
        uniprot_id = (t.get("uniprot_id") or "").strip()

        # Skip completely unmapped rows
        if not gene_symbol and not uniprot_id:
            continue

        symbol_norm = gene_symbol.upper() if gene_symbol else ""
        uniprot_norm = uniprot_id.upper() if uniprot_id else ""
        key = (symbol_norm, uniprot_norm)
        if key in by_key:
            continue

        symbol = gene_symbol or uniprot_id or symbol_norm or uniprot_norm
        by_key[key] = {
            "symbol": symbol,
            "entrez_id": None,      # can be filled later via gene mapping tools
            "uniprot_id": uniprot_id or None,
            "ensembl_id": None,
            "source": "chembl",
        }

    out = list(by_key.values())
    if DEBUG:
        print(
            f"[_get_targets_from_chembl] Resolved {len(out)} ChemBL targets "
            f"for drug {ids.name!r} (chembl_id={chembl_id})"
        )
    return out


def _get_targets_from_drugcentral(
    ids: DrugIdentifiers,
) -> Tuple[Optional[int], List[Dict[str, Any]], Dict[Tuple[str, str], int]]:
    """
    Look up HUMAN gene targets from Drug Central for a given drug using the
    local SQLite dump and utils.py.

    Returns
    -------
    (struct_id, targets, score_by_key)

    Where:
      - struct_id is Drug Central structures.id for the resolved drug
      - targets is a list of HUMAN-only targets of the form:
          {
            "symbol": "EGFR",
            "entrez_id": "1956" | None,
            "uniprot_id": "P00533" | None,
            "ensembl_id": None,
            "source": "drugcentral",
          }
      - score_by_key maps (symbol, UniProt) -> integer score,
        used for ranking when DRUG_TARGET_RETURN_MODE == "top".
    """
    try:
        from src.tools.drug_central.utils import (
            normalize_drug_identifier,
            get_drug_targets,
            project_targets_to_human_genes,
            compute_human_gene_target_scores,
        )
    except Exception as exc:
        if DEBUG:
            print(f"[_get_targets_from_drugcentral] Drug Central utils unavailable: {exc}")
        _ = DRUG_CENTRAL_DB_PATH  # keep linter happy
        return None, [], {}

    # Build candidate identifiers in a reasonable order
    candidates: List[str] = []

    if ids.drugcentral_id:
        candidates.append(str(ids.drugcentral_id).strip())

    if ids.name:
        candidates.append(ids.name.strip())

    if ids.cas_number:
        candidates.append(str(ids.cas_number).strip())

    for syn in ids.synonyms or []:
        s = (syn or "").strip()
        if s:
            candidates.append(s)

    # Deduplicate candidates while preserving order
    seen_cand: Set[str] = set()
    ordered_candidates: List[str] = []
    for c in candidates:
        key = c.lower()
        if key in seen_cand:
            continue
        seen_cand.add(key)
        ordered_candidates.append(c)

    # Resolve to one or more Drug Central struct_ids (structures.id)
    struct_ids: List[int] = []
    for token in ordered_candidates:
        rec = normalize_drug_identifier(token)
        if not rec:
            continue
        sid_raw = rec.get("struct_id")
        try:
            sid = int(sid_raw) if sid_raw is not None else None
        except (TypeError, ValueError):
            sid = None
        if sid is None:
            continue
        if sid not in struct_ids:
            struct_ids.append(sid)
        # If the user explicitly set a Drug Central id, stop after first hit
        if ids.drugcentral_id:
            break

    if not struct_ids:
        if DEBUG:
            print(f"[_get_targets_from_drugcentral] No Drug Central struct_id for drug: {ids.name!r}")
        return None, [], {}

    # Set the primary Drug Central id back onto the identifiers if missing
    primary_struct_id = struct_ids[0]
    if not ids.drugcentral_id:
        ids.drugcentral_id = str(primary_struct_id)

    # Collect raw targets and per-gene scores
    all_dc_raw_targets: List[Dict[str, Any]] = []
    score_by_key: Dict[Tuple[str, str], int] = {}

    for sid in struct_ids:
        try:
            dc_targets_raw = get_drug_targets(str(sid), include_off_target=True)
        except Exception as exc:
            if DEBUG:
                print(f"[_get_targets_from_drugcentral] Error querying targets for {sid}: {exc}")
            continue

        if not dc_targets_raw:
            continue

        all_dc_raw_targets.extend(dc_targets_raw)

        local_scores = compute_human_gene_target_scores(dc_targets_raw)
        for key, val in local_scores.items():
            score_by_key[key] = score_by_key.get(key, 0) + val

    if not all_dc_raw_targets:
        if DEBUG:
            print(
                f"[_get_targets_from_drugcentral] No act_table_full targets found "
                f"for struct_ids={struct_ids} ({ids.name!r})"
            )
        return None, [], {}

    human_targets = project_targets_to_human_genes(all_dc_raw_targets)
    out_targets = [dict(t) for t in human_targets]

    if DEBUG:
        print(
            f"[_get_targets_from_drugcentral] Resolved {len(out_targets)} HUMAN targets "
            f"for drug {ids.name!r} from Drug Central (struct_id={primary_struct_id})."
        )

    return primary_struct_id, out_targets, score_by_key


def _merge_target_records(
    records_a: List[Dict[str, Any]],
    records_b: List[Dict[str, Any]],
    context_label: str,
) -> List[Dict[str, Any]]:
    """
    Merge two lists of target records, preferring existing values and
    warning on hard-ID conflicts.

    Rules
    -----
    - Deduplicate by gene symbol (case-insensitive) + UniProt.
    - For each symbol, we keep the first non-empty hard IDs
      (entrez_id, uniprot_id, ensembl_id).
    - If a later record suggests a *different* hard ID for the same field,
      we print a warning and ignore that conflicting value.
    - Additional fields (for example Drug Central "score" etc.) are left as-is
      on the first record for that symbol.
    """
    by_key: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for rec in (records_a or []) + (records_b or []):
        if not isinstance(rec, dict):
            continue
        symbol_raw = str(rec.get("symbol", "")).strip()
        if not symbol_raw:
            continue
        uniprot_raw = str(rec.get("uniprot_id", "") or "").strip()
        symbol_key = symbol_raw.upper()
        uniprot_key = uniprot_raw.upper()
        key = (symbol_key, uniprot_key)

        existing = by_key.get(key)
        if existing is None:
            new_rec = dict(rec)
            new_rec["symbol"] = symbol_key
            if uniprot_raw:
                new_rec["uniprot_id"] = uniprot_key
            by_key[key] = new_rec
            continue

        # Merge fields, respecting "do not overwrite hard IDs on conflict"
        for field in ("entrez_id", "uniprot_id", "ensembl_id"):
            new_val = rec.get(field)
            if not new_val:
                continue
            if existing.get(field) in (None, "", 0):
                existing[field] = new_val
            elif str(existing[field]) != str(new_val):
                if DEBUG:
                    print(
                        f"[resolve_drug_gene_targets] WARNING: conflicting {field} for symbol "
                        f"{symbol_raw} while merging {context_label}: "
                        f"existing={existing[field]!r}, new={new_val!r}. Keeping existing."
                    )

        # Merge sources
        src_new = rec.get("source")
        if src_new:
            src_existing = existing.get("source")
            if not src_existing:
                existing["source"] = src_new
            else:
                parts = {s.strip() for s in str(src_existing).split(",") if s.strip()}
                parts.update(s.strip() for s in str(src_new).split(",") if s.strip())
                existing["source"] = ",".join(sorted(parts))

    return list(by_key.values())


def resolve_drug_gene_targets(drugs: List[Drug]) -> Dict[int, List[Dict[str, Any]]]:
    """
    High-level helper: resolve gene targets for one or more Drug objects
    using ChemBL and Drug Central.

    Parameters
    ----------
    drugs : List[Drug]
        Drug objects with identifiers (name, chembl_id, drugcentral_id, CAS, etc.)

    Returns
    -------
    Dict[int, List[Dict[str, Any]]]
        Mapping:
          Drug Central struct_id (structures.id) -> list of normalized HUMAN
          target records, each of the form:
            {
              "symbol": "EGFR",
              "entrez_id": "1956" | None,
              "uniprot_id": "P00533" | None,
              "ensembl_id": None,
              "source": "chembl,drugcentral" | "chembl" | "drugcentral",
            }

    Notes
    -----
    - ChemBL targets are resolved via _get_targets_from_chembl and merged with
      Drug Central targets via _merge_target_records.
    - Target expansion behavior is controlled by:
        src.tools.drug_central.DRUG_TARGET_RETURN_MODE = "none" | "top" | "all"
        src.tools.drug_central.DRUG_TARGET_TOP_N      = int (used when mode == "top")
      where "top" keeps at most DRUG_TARGET_TOP_N Drug Central HUMAN targets
      per drug, ranked by a simple per-gene score derived from act_table_full.
    """
    mode = getattr(drug_central, "DRUG_TARGET_RETURN_MODE", "top")
    top_n = getattr(drug_central, "DRUG_TARGET_TOP_N", 25)

    mode_str = str(mode or "").lower()
    if mode_str not in ("none", "top", "all"):
        if DEBUG:
            print(
                f"[resolve_drug_gene_targets] Unknown DRUG_TARGET_RETURN_MODE={mode!r}, "
                f"falling back to 'top'."
            )
        mode_str = "top"

    if not drugs or mode_str == "none":
        if DEBUG:
            print(
                f"[resolve_drug_gene_targets] Target expansion disabled "
                f"(mode={mode_str!r}) or no drugs provided."
            )
        return {}

    try:
        top_n_int = int(top_n)
    except (TypeError, ValueError):
        top_n_int = 25
    if top_n_int < 1:
        top_n_int = 1

    results: Dict[int, List[Dict[str, Any]]] = {}

    for drug in drugs:
        if not isinstance(drug, Drug):
            continue

        ids = drug.identifiers

        # Stage 1: ChemBL targets (may be empty)
        chembl_targets = _get_targets_from_chembl(ids)

        # Stage 2: Drug Central targets + struct_id + per-gene scores
        struct_id, dc_targets, score_by_key = _get_targets_from_drugcentral(ids)

        if struct_id is None:
            if DEBUG:
                print(
                    f"[resolve_drug_gene_targets] Skipping drug {ids.name!r}: "
                    "no Drug Central struct_id resolved."
                )
            continue

        merged = _merge_target_records(
            chembl_targets,
            dc_targets,
            context_label=f"targets for drug {ids.name} (struct_id={struct_id})",
        )

        if not merged:
            if DEBUG:
                print(
                    f"[resolve_drug_gene_targets] No targets after merge for "
                    f"drug {ids.name!r} (struct_id={struct_id})."
                )
            continue

        if mode_str == "top":
            # Rank primarily by Drug Central per-gene score, then by symbol
            def score_key(rec: Dict[str, Any]) -> Tuple[int, str]:
                symbol_raw = (rec.get("symbol") or "").strip().upper()
                uniprot_raw = (rec.get("uniprot_id") or "").strip().upper()
                key = (symbol_raw, uniprot_raw)
                score = score_by_key.get(key, 0)
                label = symbol_raw or uniprot_raw
                return (-score, label)

            ranked = sorted(merged, key=score_key)
            final_targets = ranked[:top_n_int]
        else:  # "all"
            final_targets = merged

        results[struct_id] = final_targets

        if DEBUG:
            print(
                f"[resolve_drug_gene_targets] Drug {ids.name!r} (struct_id={struct_id}): "
                f"{len(final_targets)} targets returned "
                f"(mode={mode_str}, ChemBL={len(chembl_targets)}, DrugCentral={len(dc_targets)})."
            )

    return results


# ------------------------------
# CLI (testing mode)
# ------------------------------

def _run_cli_self_test(drugs: List[str]) -> None:
    """
    CLI-style self-test.

    Example:
      python drug_extraction.py "Imatinib" "Atorvastatin"

    For each provided drug name, constructs a minimal Drug object,
    resolves ChemBL + DrugCentral targets, and prints a short
    text summary to stdout.
    """
    state: State = State({"drug_entities": {}})  # type: ignore[arg-type]

    mode = getattr(drug_central, "DRUG_TARGET_RETURN_MODE", "top")
    top_n = getattr(drug_central, "DRUG_TARGET_TOP_N", 25)

    print(f"[INFO] DRUG_TARGET_RETURN_MODE: {mode}")
    print(f"[INFO] DRUG_TARGET_TOP_N: {top_n}")

    print("[INFO] Running drug_extraction self-test")
    print(f"[INFO] Drugs requested: {', '.join(drugs)}")

    drug_objs: List[Drug] = []
    for name in drugs:
        identifiers = DrugIdentifiers(
            name=name,
            chembl_id=None,
            drugcentral_id=None,
            cas_number=None,
            catalog_number=None,
            synonyms=[],
        )
        drug = Drug(identifiers=identifiers)
        state["drug_entities"][name] = drug  # type: ignore[index]
        drug_objs.append(drug)

    targets_map = resolve_drug_gene_targets(drug_objs)

    for name in drugs:
        ent: Drug = state["drug_entities"][name]  # type: ignore[index]
        ent.normalize()

        print("")
        print(f"[DRUG] {name}")
        print(ent.summarize_text() or "• (no identifier details)")

        # Find any struct_id entries corresponding to this drugcentral_id if set
        ids = ent.identifiers
        printed_any = False

        candidate_struct_ids: List[int] = []
        if ids.drugcentral_id:
            try:
                candidate_struct_ids.append(int(str(ids.drugcentral_id)))
            except (TypeError, ValueError):
                pass

        # Fallback: print all entries if no explicit struct_id hint
        struct_ids_to_print = candidate_struct_ids or list(targets_map.keys())

        for sid in struct_ids_to_print:
            tlist = targets_map.get(sid, [])
            if not tlist:
                continue
            printed_any = True
            print(f"  DrugCentral struct_id={sid}: {len(tlist)} targets")
            for t in tlist[:10]:
                symbol = t.get("symbol") or "N/A"
                uniprot = t.get("uniprot_id") or "N/A"
                source = t.get("source") or "N/A"
                print(f"    - {symbol} (UniProt={uniprot}, source={source})")
            if len(tlist) > 10:
                print("    ... (truncated)")

        if not printed_any:
            print("  Targets: none resolved.")

    print("")
    print("[OK] drug_extraction self-test completed.")
    print(f"[INFO] Drugs processed: {', '.join(drugs)}")


if __name__ == "__main__":
    # Example:
    #   python drug_extraction.py "Imatinib" "Atorvastatin"
    cli_drugs = sys.argv[1:] if len(sys.argv) > 1 else ["Imatinib", "Atorvastatin"]
    _run_cli_self_test(cli_drugs)
