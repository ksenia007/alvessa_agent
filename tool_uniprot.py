"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-06-25
Updated: 2025-06-26


Description: 

UniProt tool: fetch entry + extract disease / function / GO traits.
"""

from __future__ import annotations
import requests
import re
from typing import Dict, List, Optional
from state import State

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


def extract_locations_ISOFORM(entry: Dict) -> List[str]:
    """Extract subcellular locations, enriched with isoform metadata and tissue/induction notes.
    Keeps molecule, note, location values; maps 'Isoform N' to ALTERNATIVE PRODUCTS isoformIds/synonyms.
    Also parses TISSUE SPECIFICITY / INDUCTION texts and attaches them to matching isoforms.
    Uses remove_pubmed_text if available to strip PubMed refs."""
    if DEBUG:
        print("Extracting subcellular locations (smart)")

    import re

    summaries: List[str] = []

    # -------- 1) Build isoform map from ALTERNATIVE PRODUCTS --------
    # iso_map['1'] = {'ids': ['P04637-1'], 'synonyms': ['p53','p53alpha'], 'status': 'Displayed'}
    iso_map: Dict[str, Dict[str, Any]] = {}
    for item in entry.get("comments", []):
        if item.get("commentType") == "ALTERNATIVE PRODUCTS":
            for iso in item.get("isoforms", []) or []:
                num = None
                name_obj = iso.get("name")
                if isinstance(name_obj, dict):
                    num = name_obj.get("value")
                elif isinstance(name_obj, str):
                    num = name_obj
                if not num:
                    continue
                syns = []
                for s in iso.get("synonyms", []) or []:
                    if isinstance(s, dict) and isinstance(s.get("value"), str):
                        syns.append(s["value"])
                    elif isinstance(s, str):
                        syns.append(s)
                ids = list(iso.get("isoformIds", []) or [])
                status = iso.get("isoformSequenceStatus")
                iso_map[str(num).strip()] = {
                    "ids": ids,
                    "synonyms": syns,
                    "status": status,
                }

    # -------- 2) Pull tissue specificity / induction notes and attach to isoforms if text mentions them --------
    # iso_notes['2'] = ["Isoform 2 is expressed ...", "Isoform 2 is not induced ..."]
    iso_notes: Dict[str, List[str]] = {}
    general_notes: List[str] = []
    for item in entry.get("comments", []):
        ctype = item.get("commentType")
        if ctype not in ("TISSUE SPECIFICITY", "INDUCTION"):
            continue
        texts = item.get("texts") or []
        for t in texts:
            val = ""
            if isinstance(t, dict) and isinstance(t.get("value"), str):
                val = t["value"]
            elif isinstance(t, str):
                val = t
            if not val:
                continue
            # Clean PubMed refs if utility exists
            try:
                val = remove_pubmed_text(val)
            except NameError:
                pass

            # Attach to specific isoforms if "Isoform N" is mentioned; else keep as general
            mentioned = False
            for m in re.finditer(r"\bIsoform\s+(\d+)\b", val, flags=re.IGNORECASE):
                n = m.group(1)
                iso_notes.setdefault(n, []).append(val)
                mentioned = True
            if not mentioned:
                general_notes.append(f"{ctype}: {val}")

    # -------- 3) Build enriched SUBCELLULAR LOCATION summaries --------
    for item in entry.get("comments", []):
        if item.get("commentType") != "SUBCELLULAR LOCATION":
            continue

        molecule = item.get("molecule") or ""  # e.g., "Isoform 1"
        # Try to extract the number part to link to iso_map/iso_notes
        mnum = None
        m = re.search(r"\bIsoform\s+(\d+)\b", molecule, flags=re.IGNORECASE)
        if m:
            mnum = m.group(1)

        # Note text normalization (supports note.value or note.texts[].value)
        note_text = ""
        note = item.get("note")
        if isinstance(note, dict):
            if isinstance(note.get("value"), str):
                note_text = note.get("value", "")
            elif isinstance(note.get("texts"), list):
                vals = []
                for t in note["texts"]:
                    if isinstance(t, dict) and isinstance(t.get("value"), str):
                        vals.append(t["value"])
                note_text = " ".join(vals).strip()
        elif isinstance(note, str):
            note_text = note

        try:
            if note_text:
                note_text = remove_pubmed_text(note_text)
        except NameError:
            pass

        # Collect unique location values
        locs_raw = []
        for sl in item.get("subcellularLocations", []) or []:
            loc = sl.get("location")
            val = None
            if isinstance(loc, dict):
                val = loc.get("value")
            elif isinstance(loc, str):
                val = loc
            if val:
                locs_raw.append(val)
        # de-dup keep order
        seen = set(); locations = []
        for v in locs_raw:
            if v not in seen:
                seen.add(v); locations.append(v)

        # Enrich with isoform metadata
        meta_parts: List[str] = []
        if molecule:
            meta_parts.append(molecule)
        if mnum and mnum in iso_map:
            ids = ", ".join(iso_map[mnum].get("ids") or [])
            syns = ", ".join(iso_map[mnum].get("synonyms") or [])
            status = iso_map[mnum].get("status")
            if ids:
                meta_parts.append(f"IDs: {ids}")
            if syns:
                meta_parts.append(f"aka: {syns}")
            if status:
                meta_parts.append(f"status: {status}")

        # Attach tissue/induction notes that reference this isoform
        extra_notes = []
        if mnum and mnum in iso_notes:
            extra_notes = iso_notes[mnum]

        # Build summary string
        parts: List[str] = []
        if meta_parts:
            parts.append(" | ".join(meta_parts))
        if note_text:
            parts.append(note_text)
        if locations:
            parts.append("Locations: " + ", ".join(locations))
        for en in extra_notes:
            parts.append(en)

        summary = " — ".join([p for p in parts if p]).strip()
        if summary:
            summaries.append(summary)

    # -------- 4) Optionally append general tissue/induction notes (not isoform-specific) --------
    for gnote in general_notes:
        if gnote and gnote not in summaries:
            summaries.append(gnote)

    return summaries

                

def extract_summaries(entry: Dict) -> Dict[str, str]:
    """Extract summaries from a UniProt entry."""
    if DEBUG:
        print("Extracting summaries")
    summaries: Dict[str, str] = {}
    diseases = extract_disease_from_uniprot_entry(entry)
    function = extract_function_from_uniprot_entry(entry)
    go_terms = extract_GO_from_uniprot_entry(entry)
    uniprotID = entry.get("primaryAccession", "")
    locations = extract_locations(entry)
    print('LOCATIONS', locations)
    test = extract_locations_ISOFORM(entry)
    print('TEST', test)
    raise Exception('STOP')
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
    return summaries

def additional_info_extraction(entry: Dict) -> Dict[str, str]:
    """Extract additional info from a UniProt entry."""
    if DEBUG:
        print("Extracting additional info")
    info: Dict[str, str] = {}
    # extract organism: scientificName
    # fullName': {'value'
    # 'genes': [{'geneName': {'value': 'TP53'}, 'synonyms': [{'value': 'P53'}]}]...
    
    
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
            info = extract_summaries(entry_base)
            print('INFO', info)
            # raise Exception('STOP')

            # add to the gene object 
            gene_objs[gene].add_function_label(info.get('function', ""))
            gene_objs[gene].add_many_diseases(info.get('diseases', []))
            gene_objs[gene].add_many_go_terms(info.get('go_terms', []))
            if "uniprotID" in info:
                gene_objs[gene].set_gene_ids(uniprot_id = info["uniprotID"])
        
    return 


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
