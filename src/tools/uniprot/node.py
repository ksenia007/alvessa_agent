"""
Description:

UniProt tool: fetch entry + extract disease / function / GO traits.
"""

from __future__ import annotations
import requests
import re
from typing import Dict, List, Optional, Any, Iterable

from src.state import State
from src.tools.base import Node

DEBUG = True

def get_uniprot_entry_for_gene(gene_symbol: str) -> Optional[Dict]:
    """
    Retrieve the *reviewed* UniProtKB record for a human gene symbol.

    Returns
    -------
    dict | None
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene:{gene_symbol} AND organism_id:9606 AND reviewed:true",
        "format": "json",
        "size": 1,
    }

    if DEBUG:
        print(f"Fetching UniProt entry for gene: {gene_symbol}")

    try:
        resp = requests.get(base_url, params=params, timeout=12)
        resp.raise_for_status()
        results = resp.json().get("results", [])
        print('UNIPROT RESULTS', results)
        return results[0] if results else None
    except Exception as exc:
        print(f"Error fetching UniProt entry: {exc}")
        return None


def extract_disease_from_uniprot_entry(entry: Dict) -> List[str]:
    """Pull disease names from a UniProt entry."""
    if DEBUG:
        print("Extracting disease traits")
    traits: set[str] = set()
    for item in entry.get("comments", []):
        if item.get("commentType") == "DISEASE":
            disease = item.get("disease", {})
            name = disease.get("diseaseId") or disease.get("acronym") or disease.get(
                "name", {}
            ).get("value")
            if name:
                traits.add(name)
    return sorted(traits)


def extract_function_from_uniprot_entry(entry: Dict) -> List[str]:
    """Pull free-text *Function* comment blocks."""
    if DEBUG:
        print("Extracting function traits")
    traits: set[str] = set()
    for item in entry.get("comments", []):
        if item.get("commentType") == "FUNCTION":
            for text_block in item.get("texts", []):
                traits.add(text_block.get("value", ""))
    return sorted(traits)


def extract_GO_from_uniprot_entry(entry: Dict) -> List[str]:
    """Extract GO terms from cross-references."""
    if DEBUG:
        print("Extracting GO terms")
    traits: set[str] = set()
    for ref in entry.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "GO":
            for prop in ref.get("properties", []):
                if prop.get("key") == "GoTerm":
                    traits.add(prop["value"])
    return sorted(traits)

def remove_pubmed_text(description: str) -> str:
    cleaned = re.sub(r'\s*PubMed:\d+\s*,?', '', description)
    cleaned = re.sub(r'\(\s*,*\s*\)', '()', cleaned)
    return re.sub(r'\(\s*\)', '', cleaned).strip()

def dedup_keep(seq):
    seen = set(); out = []
    for x in seq:
        if x not in seen:
            seen.add(x); out.append(x)
    return out


def _isoform_sort_key(iso_id: str, name: Optional[str]) -> tuple:
    """
    Sort by numeric isoform number if available, else by name, then by id.
    Examples:
      "P04637-2" -> (0, 2, "", "P04637-2")
      "Q9XYZ1-11" -> (0, 11, "", "Q9XYZ1-11")
      fallback -> (1, inf, name_or_empty, iso_id)
    """
    m = re.search(r"-(\d+)$", iso_id or "")
    if m:
        return (0, int(m.group(1)), "", iso_id)
    return (1, float("inf"), (name or ""), iso_id)


def summarize_isoform_localization(
    isoforms: Dict[str, Dict[str, Any]],
    *,
    line_separator: str = " |-| ",
    joiner_meta: str = "; ",
    joiner_locs: str = ", ",
    include_general: bool = True,
    max_notes_per_line: Optional[int] = None,
) -> str:
    """
    Build a compact, readable paragraph summarizing isoform localizations.
    Expects a dict shaped like the output of `parse_uniprot_isoforms_localization`.

    Returns a single string (paragraph). If nothing to report, returns "".

    Example fragment per isoform:
      "Isoform 2 [P04637-2] | aka: p53beta; status: Described; Locations: Nucleus, Cytoplasm; <notes>"

    Parameters:
      line_separator: placed between isoform summaries (default " |-| ").
      joiner_meta: joins metadata pieces after the header (default "; ").
      joiner_locs: joins multiple location names (default ", ").
      include_general: append a final “general” bucket if present.
      max_notes_per_line: cap number of notes appended per isoform/general (None = no cap).
    """
    if not isoforms:
        return ""

    lines: List[str] = []

    # --- Specific isoforms first (skip the general bucket) ---
    specific_items = [
        (iso_id, rec) for iso_id, rec in isoforms.items()
        if iso_id != "general_localization"
    ]
    # Sort for stable output: numeric isoform order when possible
    specific_items.sort(key=lambda kv: _isoform_sort_key(kv[0], kv[1].get("name")))

    for iso_id, rec in specific_items:
        name: str = rec.get("name") or iso_id
        aliases: List[str] = rec.get("aliases") or []
        status: Optional[str] = rec.get("status")
        locs: List[str] = dedup_keep(rec.get("locations") or [])
        notes: List[str] = rec.get("notes") or []

        header = f"{name} [{iso_id}]"
        meta_bits: List[str] = []
        if aliases:
            meta_bits.append("aka: " + ", ".join(aliases))
        if status:
            meta_bits.append(f"status: {status}")

        parts: List[str] = [header + ((" | " + joiner_meta.join(meta_bits)) if meta_bits else "")]
        if locs:
            parts.append("Locations: " + joiner_locs.join(locs))

        if max_notes_per_line is not None:
            notes = notes[:max_notes_per_line]
        for n in notes:
            if n:
                parts.append(n)

        line = "; ".join(p for p in parts if p).strip()
        if line:
            lines.append(line)

    # --- General bucket last (optional) ---
    if include_general and "general_localization" in isoforms:
        gen = isoforms["general_localization"]
        gen_locs: List[str] = dedup_keep(gen.get("locations") or [])
        gen_notes: List[str] = gen.get("notes") or []

        parts: List[str] = []
        if gen_locs:
            parts.append("General locations: " + joiner_locs.join(gen_locs))
        if max_notes_per_line is not None:
            gen_notes = gen_notes[:max_notes_per_line]
        for n in gen_notes:
            if n:
                parts.append(n)

        if parts:
            lines.append(" — ".join(parts))

    return line_separator.join(lines)




def parse_uniprot_isoforms_localization(entry: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """
    Return:
      {
        "P04637-2": {"name":"Isoform 2","aliases":[...],"status":"Described","locations":[...],"notes":[...]},
        ...,
        "general_localization": {"name":"general_localization","aliases":[],"status":None,"locations":[...],"notes":[...]}
      }
    """
    import re

    iso_map_by_num: Dict[str, Dict[str, Any]] = {}
    out: Dict[str, Dict[str, Any]] = {}

    # 1) ALTERNATIVE PRODUCTS -> map isoform number -> (ids, aliases, status, name)
    for item in entry.get("comments", []):
        if item.get("commentType") != "ALTERNATIVE PRODUCTS":
            continue
        for iso in item.get("isoforms", []) or []:
            nmv = iso.get("name", {})
            num = nmv.get("value") if isinstance(nmv, dict) else (nmv if isinstance(nmv, str) else None)
            if not num:
                continue
            aliases = []
            for syn in iso.get("synonyms", []) or []:
                if isinstance(syn, dict) and isinstance(syn.get("value"), str):
                    aliases.append(syn["value"])
                elif isinstance(syn, str):
                    aliases.append(syn)
            ids = list(iso.get("isoformIds", []) or [])
            status = iso.get("isoformSequenceStatus")
            num_map_entry = {
                "ids": ids,
                "aliases": dedup_keep(aliases),
                "status": status,
                "name": f"Isoform {str(num).strip()}",
            }
            iso_map_by_num[str(num).strip()] = num_map_entry

    # 2) Collect isoform-specific and general notes from TISSUE SPECIFICITY / INDUCTION
    iso_notes: Dict[str, List[str]] = {}
    general_notes: List[str] = []
    for item in entry.get("comments", []):
        ctype = item.get("commentType")
        if ctype not in ("TISSUE SPECIFICITY", "INDUCTION"):
            continue
        for t in item.get("texts") or []:
            val = t.get("value") if isinstance(t, dict) else (t if isinstance(t, str) else "")
            if not val:
                continue
            try:
                val = remove_pubmed_text(val)  # clean refs
            except NameError:
                pass
            tagged = False
            for m in re.finditer(r"\bIsoform\s+(\d+)\b", val, flags=re.IGNORECASE):
                n = m.group(1)
                iso_notes.setdefault(n, []).append(val)
                tagged = True
            if not tagged:
                general_notes.append(f"{ctype}: {val}")

    # 3) SUBCELLULAR LOCATION to build records (prefer isoform IDs; else general)
    def rec_for_general():
        return out.setdefault("general_localization", {"name": "general_localization", "aliases": [], "status": None, "locations": [], "notes": []})

    for item in entry.get("comments", []):
        if item.get("commentType") != "SUBCELLULAR LOCATION":
            continue

        molecule = item.get("molecule") or ""
        m = re.search(r"\bIsoform\s+(\d+)\b", molecule, flags=re.IGNORECASE)
        iso_id = None
        num = m.group(1) if m else None
        name = molecule or ""
        aliases: List[str] = []
        status = None

        if num and num in iso_map_by_num and iso_map_by_num[num]["ids"]:
            iso_id = iso_map_by_num[num]["ids"][0]
            name = iso_map_by_num[num]["name"] or name
            aliases = iso_map_by_num[num]["aliases"] or []
            status = iso_map_by_num[num]["status"]

        # Note text
        note_text = ""
        note = item.get("note")
        if isinstance(note, dict):
            if isinstance(note.get("value"), str):
                note_text = note["value"]
            elif isinstance(note.get("texts"), list):
                parts = []
                for tx in note["texts"]:
                    if isinstance(tx, dict) and isinstance(tx.get("value"), str):
                        parts.append(tx["value"])
                note_text = " ".join(parts).strip()
        elif isinstance(note, str):
            note_text = note
        try:
            if note_text:
                note_text = remove_pubmed_text(note_text)
        except NameError:
            pass

        # Locations
        locs_raw = []
        for sl in item.get("subcellularLocations", []) or []:
            loc = sl.get("location")
            v = loc.get("value") if isinstance(loc, dict) else (loc if isinstance(loc, str) else None)
            if v:
                locs_raw.append(v)
        locations = dedup_keep(locs_raw)

        # Upsert target record
        if iso_id:
            rec = out.setdefault(iso_id, {"name": name, "aliases": [], "status": status, "locations": [], "notes": []})
            for a in aliases:
                if a not in rec["aliases"]:
                    rec["aliases"].append(a)
            for v in locations:
                if v not in rec["locations"]:
                    rec["locations"].append(v)
            if note_text and note_text not in rec["notes"]:
                rec["notes"].append(note_text)
            # attach iso-specific notes
            if num and num in iso_notes:
                for txt in iso_notes[num]:
                    if txt not in rec["notes"]:
                        rec["notes"].append(txt)
        else:
            rec = rec_for_general()
            for v in locations:
                if v not in rec["locations"]:
                    rec["locations"].append(v)
            if note_text and note_text not in rec["notes"]:
                rec["notes"].append(note_text)

    # 4) Add general notes (no isoform mention)
    if general_notes:
        rec = out.setdefault("general_localization", {"name": "general_localization", "aliases": [], "status": None, "locations": [], "notes": []})
        for txt in general_notes:
            if txt not in rec["notes"]:
                rec["notes"].append(txt)

    # 5) Ensure any isoforms with notes but no SUBCELLULAR LOCATION still appear
    for num, meta in iso_map_by_num.items():
        if not meta["ids"]:
            continue
        iso_id = meta["ids"][0]
        if iso_id not in out and num in iso_notes:
            out[iso_id] = {
                "name": meta["name"],
                "aliases": meta["aliases"],
                "status": meta["status"],
                "locations": [],
                "notes": dedup_keep(iso_notes[num]),
            }
            
    summary = summarize_isoform_localization(out)

    return out, summary

                

def extract_summaries(entry: Dict) -> Dict[str, str]:
    """Extract summaries from a UniProt entry."""
    if DEBUG:
        print("Extracting summaries")
    summaries: Dict[str, str] = {}
    diseases = extract_disease_from_uniprot_entry(entry)
    function = extract_function_from_uniprot_entry(entry)
    go_terms = extract_GO_from_uniprot_entry(entry)
    uniprotID = entry.get("primaryAccession", "")
    locations, text_summary_loc = parse_uniprot_isoforms_localization(entry)
    print('parse_uniprot_isoforms_localization', locations)
    if locations:
        summaries["locations"] = locations
    if uniprotID:
        summaries["uniprotID"] = uniprotID
    if diseases:
        summaries["diseases"] = diseases #remove_pubmed_text(", ".join(diseases))
    if function:
        summaries["function"] = remove_pubmed_text(", ".join(function))
    if go_terms:
        summaries["go_terms"] = go_terms #remove_pubmed_text(", ".join(go_terms))
    if DEBUG:
        print(f"Extracted summaries: {summaries}")
        
    # create disease text summary:
    summary_full = "Information collected from UniProt: \n"
    if locations:
        summary_full += f"Isoform localizations: {text_summary_loc}\n"
    if diseases:
        summary_full += "\n Associated diseases: " + ", ".join(diseases) + "."
    
    return summaries, f"*UniProt: {summary_full.strip().replace(chr(10), ' ')}"

    
def uniprot_node(state: "State") -> "State":
    """Download UniProt entries for every gene symbol."""

    # gwas_linked_genes = list(set(state.get("gwas_linked_genes", [])))
    gene_objs = state.get("gene_entities", {})
    genes = list(gene_objs.keys())

    uniprot_entries_base: Dict[str, Dict] = {}
    # uniprot_entries_gwas: Dict[str, Dict] = {}  removed GWAS genes for now at least
    
    if DEBUG:
        print(f"Fetching UniProt entries for {len(genes)} genes: {genes}")

    for gene in genes:
        entry_base = get_uniprot_entry_for_gene(gene)
        print('ENTRY BASE', entry_base)
        if entry_base:
            info, summary_full = extract_summaries(entry_base)
            # raise Exception('STOP')
            # add to the gene object 
            gene_objs[gene].add_function_label(info.get('function', ""))
            gene_objs[gene].add_many_diseases(info.get('diseases', []))
            gene_objs[gene].add_many_go_terms(info.get('go_terms', []))
            gene_objs[gene].set_isoform_localizations(info.get('locations', {}))
            if "uniprotID" in info:
                gene_objs[gene].set_gene_ids(uniprot_id = info["uniprotID"])
            gene_objs[gene].add_tool("uniprot_entry")
            gene_objs[gene].update_text_summaries(summary_full)

    return 


NODES: tuple[Node, ...] = (
    Node(
        name="uniprot_base",
        entry_point=uniprot_node,
        description="Fetch UniProt entries for genes, extracting diseases, GO terms, and isoform annotations. Helps to expand the base annotations.",
        aliases=("uniprot_gwas",),
    ),
)

# not sure what these were used for, leaving commented out for now
# def trait_disease_extraction_node(state: "State") -> "State":
#     """Attach disease traits pulled from UniProt entries."""
#     entries = state.get("uniprot_entries", {})
#     gene_traits: Dict[str, List[str]] = {
#         g: extract_disease_from_uniprot_entry(e) for g, e in entries.items()
#     }
#     # prune empties
#     gene_traits = {g: t for g, t in gene_traits.items() if t}
#     return {"gene_disease_traits": gene_traits}


# def trait_function_extraction_node(state: "State") -> "State":
#     """Attach free-text functional annotations."""
#     entries = state.get("uniprot_entries", {})
#     gene_traits: Dict[str, List[str]] = {
#         g: extract_function_from_uniprot_entry(e) for g, e in entries.items()
#     }
#     gene_traits = {g: t for g, t in gene_traits.items() if t}
#     return {"gene_function_traits": gene_traits}


# def trait_GO_extraction_node(state: "State") -> "State":
#     """Attach GO term lists."""
#     entries = state.get("uniprot_entries", {})
#     gene_traits: Dict[str, List[str]] = {
#         g: extract_GO_from_uniprot_entry(e) for g, e in entries.items()
#     }
#     gene_traits = {g: t for g, t in gene_traits.items() if t}
#     return {"gene_GO_traits": gene_traits}


# def has_uniprot_entries(state: "State") -> bool:
#     """Whether any UniProt entries have been fetched (LangGraph edge helper)."""
#     return bool(state.get("uniprot_entries"))
