# src/tools/drug_central/node.py
# ===========================================================
# DrugCentral Drug-Target, Indication, and Safety Tool
# ===========================================================
#
# Author: Dmitri Kosenkov
# Created: 2025-11-17
# Updated: 2025-12-02
#
# This module queries a local DrugCentral SQLite database for
# drug-centric evidence: structure, regulatory status, targets,
# indications, pharmacologic actions, and safety profile.
#
# Inputs:
#   - state["drugs"]: list of free-text or ID-like drug inputs
#   - state["drug_entities"]: dict keyed by canonical drug key,
#       each entry being a Drug object (src.alvessa.domain.drug_class.Drug)
#
# Outputs:
#   - Enriched Drug objects in state["drug_entities"]
#   - Optional "drug_central_html" with inline viewer
#

import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.state import State
from src.tools.base import Node
from src.config import DEBUG

from src.alvessa.domain.drug_class import Drug, DrugIdentifiers

from src.tools.drug_central.utils import (
    normalize_drug_identifier,
    get_drug_regulatory_status,
    get_drug_targets,
    get_drug_indications,
    get_drug_safety_profile,
    get_cross_references,
    get_target_class_summary,
    get_drug_pharmacologic_actions,
    make_drug_summary_text,
    interpretation_notes,
    log,
    inject_frontend_assets,
    project_targets_to_human_genes,
)

from src.tools.drug_central import OUTPUT_DIR, HTML_TEMPLATE, CSS_TEMPLATE, JS_TEMPLATE


# ------------------------------
# Helpers for Drug objects
# ------------------------------
def _coerce_dict_to_drug(key: str, payload: Dict[str, Any]) -> Drug:
    """
    Backwards-compatible adapter: turn a legacy dict entry from
    state["drug_entities"][key] into a Drug object.
    """
    # Support both top-level dict and {"identifiers": {...}} layout
    id_src = payload.get("identifiers") if isinstance(payload.get("identifiers"), dict) else payload

    name = (id_src.get("name") or "").strip() or key
    chembl_id = id_src.get("chembl_id")
    drugcentral_id = id_src.get("drugcentral_id")
    cas_number = id_src.get("cas_number")
    catalog_number = id_src.get("catalog_number")
    synonyms = id_src.get("synonyms") or []
    if not isinstance(synonyms, list):
        synonyms = [str(synonyms)]

    ids = DrugIdentifiers(
        name=name,
        chembl_id=chembl_id,
        drugcentral_id=drugcentral_id,
        cas_number=cas_number,
        catalog_number=catalog_number,
        synonyms=synonyms,
    )
    d = Drug(identifiers=ids)

    # Mentions (if present)
    mentions = payload.get("mentions") if isinstance(payload.get("mentions"), list) else []
    for m in mentions:
        if isinstance(m, dict):
            d.add_mention(m)

    # Tools (if present)
    tools = payload.get("tools") if isinstance(payload.get("tools"), list) else []
    for t in tools:
        try:
            d.add_tool(str(t))
        except Exception:
            pass

    # Text summaries from tools (if present)
    summaries = payload.get("text_summaries_from_tools")
    if isinstance(summaries, list):
        setattr(d, "text_summaries_from_tools", list(map(str, summaries)))

    # Any other useful payload fields we want to keep can be attached
    # as attributes for now:
    for extra_key in ("cross_references", "target_class_summary", "pharmacologic_actions"):
        if extra_key in payload:
            setattr(d, extra_key, payload[extra_key])

    return d


def _normalize_drug_entities_in_state(state: "State") -> Dict[str, Drug]:
    """
    Ensure state["drug_entities"] is a dict[str, Drug].
    Converts legacy dict entries into Drug objects.
    """
    raw = state.get("drug_entities") or {}
    norm: Dict[str, Drug] = {}

    if not isinstance(raw, dict):
        if DEBUG:
            print(f"[drug_central] Warning: state['drug_entities'] is not a dict: {type(raw)}")
        state["drug_entities"] = {}
        return norm

    for k, v in raw.items():
        if isinstance(v, Drug):
            norm[k] = v
        elif isinstance(v, dict):
            if DEBUG:
                print(f"[drug_central] Coercing legacy dict to Drug for key {k!r}")
            norm[k] = _coerce_dict_to_drug(str(k), v)
        else:
            if DEBUG:
                print(f"[drug_central] Skipping unexpected drug_entities[{k!r}] type={type(v)}")

    state["drug_entities"] = norm  # write back normalized mapping
    return norm


def _prepare_drugs_drug_central(state: "State") -> List[str]:
    """
    Resolve, deduplicate, and prepare drug inputs for the DrugCentral agent.
    Uses state["drugs"] and state["drug_name"] (if present).
    """
    drugs = state.get("drugs") or []
    if not drugs and state.get("drug_name"):
        drugs = [state["drug_name"]]

    if not drugs:
        log("No drug inputs provided.")
        state["used_tools"] = state.get("used_tools", []) + ["drug_central"]
        return []

    # Deduplicate while preserving order
    unique_drugs: List[str] = []
    seen = set()
    for d in drugs:
        if d not in seen:
            seen.add(d)
            unique_drugs.append(d)

    # Ensure drug_entities is normalized (Drug objects)
    _normalize_drug_entities_in_state(state)

    return unique_drugs


def _find_or_create_drug_entity(
    state: "State",
    struct_id: int,
    drug_rec: Dict[str, Any],
    input_token: str,
) -> Tuple[str, Drug]:
    """
    Find an existing Drug object in state['drug_entities'] that matches
    this DrugCentral struct_id or name; otherwise create a new one.

    Returns (key, drug_obj).
    """
    drug_entities: Dict[str, Drug] = _normalize_drug_entities_in_state(state)

    struct_id_str = str(struct_id)
    name_dc = (drug_rec.get("name") or input_token).strip()

    # 1) Try to find by existing drugcentral_id
    for key, d in drug_entities.items():
        ids = getattr(d, "identifiers", None)
        if ids is None:
            continue
        if getattr(ids, "drugcentral_id", None) == struct_id:
            return key, d

    # 2) Try to find by name (case-insensitive)
    for key, d in drug_entities.items():
        ids = getattr(d, "identifiers", None)
        if ids is None:
            continue
        if ids.name and ids.name.strip().lower() == name_dc.lower():
            # Fill in drugcentral_id if missing
            if not getattr(ids, "drugcentral_id", None):
                ids.drugcentral_id = struct_id
            return key, d

    # 3) Create a new Drug object with canonical key
    canon_key = f"DC:{struct_id_str}"
    if canon_key in drug_entities:
        d = drug_entities[canon_key]
    else:
        ids = DrugIdentifiers(
            name=name_dc,
            chembl_id=None,
            drugcentral_id=struct_id,
            cas_number=drug_rec.get("cas_number"),
            catalog_number=None,
            synonyms=drug_rec.get("synonyms") or [],
        )
        d = Drug(identifiers=ids)
        d.add_tool("DrugCentral")
        drug_entities[canon_key] = d

    state["drug_entities"] = drug_entities
    return canon_key, d


def _append_tool_summary(drug_obj: Drug, text: str) -> None:
    """
    Attach a text summary string to a Drug object under
    .text_summaries_from_tools (list of str).
    """
    summaries = getattr(drug_obj, "text_summaries_from_tools", None)
    if summaries is None:
        summaries = []
        setattr(drug_obj, "text_summaries_from_tools", summaries)
    if isinstance(summaries, list):
        summaries.append(text)


# ------------------------------
# Main DrugCentral agent
# ------------------------------
def drug_central_agent(state: "State") -> "State":
    """
    Query DrugCentral for each drug in the state.
    Retrieves structural, regulatory, target, indication,
    pharmacologic action, and safety data.

    Also computes a normalized, human-only gene-target projection
    using project_targets_to_human_genes so that:

      - Genes added via the drug->gene expansion in the
        entity_extraction / drug_extraction pipeline, and
      - Human gene targets shown in DrugCentral text summaries

    are consistent.
    """
    drugs = _prepare_drugs_drug_central(state)
    if not drugs:
        return state

    drug_central_data_all: Dict[str, Dict[str, Any]] = {}
    any_data = False

    for input_token in drugs:
        log(f"Processing drug input: {input_token!r}")

        drug_rec = normalize_drug_identifier(input_token)
        if not drug_rec:
            log(f"Could not resolve drug: {input_token!r}")
            continue

        struct_id = drug_rec.get("struct_id")
        if struct_id is None:
            log(f"Resolved drug has no struct_id: {drug_rec}")
            continue

        # Find or create a proper Drug object in state["drug_entities"]
        canon_key, drug_obj = _find_or_create_drug_entity(
            state=state,
            struct_id=struct_id,
            drug_rec=drug_rec,
            input_token=input_token,
        )

        # Ensure identifiers reflect DrugCentral info
        ids = drug_obj.identifiers
        if not getattr(ids, "drugcentral_id", None):
            ids.drugcentral_id = struct_id
        if not ids.name and drug_rec.get("name"):
            ids.name = drug_rec["name"]
        if not ids.cas_number and drug_rec.get("cas_number"):
            ids.cas_number = drug_rec["cas_number"]
        # Merge synonyms
        syns_dc = drug_rec.get("synonyms") or []
        if syns_dc:
            existing_syns = set((ids.synonyms or []))
            for s in syns_dc:
                s_norm = (s or "").strip()
                if s_norm and s_norm not in existing_syns:
                    ids.synonyms.append(s_norm)
                    existing_syns.add(s_norm)

        # Fetch layered information
        struct_id_str = str(struct_id)
        regulatory = get_drug_regulatory_status(struct_id_str)
        targets = get_drug_targets(struct_id_str, include_off_target=True)
        indications = get_drug_indications(struct_id_str)
        safety = get_drug_safety_profile(struct_id_str)
        xrefs = get_cross_references(struct_id_str)
        tclass_summary = get_target_class_summary(struct_id_str)
        pharm_actions = get_drug_pharmacologic_actions(struct_id_str)

        # Normalized human gene targets, consistent with drug_extraction
        human_gene_targets = project_targets_to_human_genes(targets)

        if regulatory or targets or indications or safety or pharm_actions or human_gene_targets:
            any_data = True

        # Text summary and interpretation notes
        summary = make_drug_summary_text(
            drug=drug_rec,
            regulatory=regulatory,
            targets=targets,
            indications=indications,
            safety=safety,
            pharm_actions=pharm_actions,
            human_gene_targets=human_gene_targets,
        )
        _append_tool_summary(drug_obj, "DrugCentral evidence:\n" + summary)

        # Attach cross references, target class summary, pharmacologic actions, and gene targets
        if xrefs:
            setattr(drug_obj, "cross_references", xrefs)
        if tclass_summary:
            setattr(drug_obj, "target_class_summary", tclass_summary)
        if pharm_actions:
            setattr(drug_obj, "pharmacologic_actions", pharm_actions)
        if human_gene_targets:
            setattr(drug_obj, "gene_targets_human", human_gene_targets)

        # JSON payload for HTML viewer
        drug_central_data_all[canon_key] = {
            "drug": drug_rec,
            "regulatory": regulatory,
            "targets": targets,
            "indications": indications,
            "safety": safety,
            "cross_refs": xrefs,
            "target_class_summary": tclass_summary,
            "pharm_actions": pharm_actions,
            # Expose normalized human gene targets to the frontend
            "gene_targets_human": human_gene_targets,
        }

    # Add interpretation notes to the last processed Drug entity
    if any_data and state.get("drug_entities"):
        # Normalize again in case keys were updated
        drug_entities: Dict[str, Drug] = _normalize_drug_entities_in_state(state)
        if drug_entities:
            last_key = list(drug_entities.keys())[-1]
            ent = drug_entities.get(last_key)
            if ent is not None:
                _append_tool_summary(
                    ent,
                    interpretation_notes(include_safety=True),
                )

    # Build HTML frontend, if templates are present
    try:
        with open(HTML_TEMPLATE, "r", encoding="utf-8") as f:
            html_template = f.read()
        with open(CSS_TEMPLATE, "r", encoding="utf-8") as f:
            css_template = f.read()
        with open(JS_TEMPLATE, "r", encoding="utf-8") as f:
            js_template = f.read()
    except FileNotFoundError as e:
        # Non-fatal: we simply skip HTML generation
        log(f"Missing DrugCentral frontend template: {e.filename}")
        html = ""
    else:
        html = inject_frontend_assets(
            html_template,
            css_template,
            js_template,
            drug_central_data_all,
            ", ".join(drugs),
        )

    state["used_tools"] = state.get("used_tools", []) + ["drug_central"]
    if html:
        state["drug_central_html"] = html

    return state


# ------------------------------
# CLI (testing mode)
# ------------------------------
if __name__ == "__main__":
    # Example:
    #   python node.py "Imatinib" "Atorvastatin"
    drugs = sys.argv[1:] if len(sys.argv) > 1 else ["Imatinib", "Atorvastatin", "CHEMBL521"]
    base_name = drugs[0] if len(drugs) == 1 else f"{drugs[0]}_plus{len(drugs)-1}"

    state = State(
        {
            "drugs": drugs,
            "drug_entities": {},
        }
    )

    result = drug_central_agent(state)
    drug_entities: Dict[str, Any] = result.get("drug_entities", {})

    html_out = OUTPUT_DIR / f"{base_name}_drug_central.html"
    html_out.write_text(
        result.get("drug_central_html", "<p>No HTML produced.</p>"),
        encoding="utf-8",
    )

    txt_out = OUTPUT_DIR / f"{base_name}_drug_central.txt"
    lines: List[str] = []
    for key, obj in drug_entities.items():
        if isinstance(obj, Drug):
            summaries = getattr(obj, "text_summaries_from_tools", []) or []
        elif isinstance(obj, dict):
            summaries = obj.get("text_summaries_from_tools") or []
        else:
            continue

        if summaries:
            lines.append(str(key))
            lines.extend(str(s) for s in summaries)
            lines.append("")

    txt_out.write_text(
        "\n".join(lines) if lines else "No summaries produced.",
        encoding="utf-8",
    )

    print("[OK] DrugCentral outputs generated:")
    print(f"  * HTML file: {html_out.resolve()}")
    print(f"  * TXT file:  {txt_out.resolve()}")
    print(f"[INFO] Drugs processed: {', '.join(drugs)}")


# ------------------------------
# Node Registration
# ------------------------------
NODES: tuple[Node, ...] = (
    Node(
        name="drug_central",
        entry_point=drug_central_agent,
        description=(
            "Query Drug Central for drug-centric information given one or more drugs "
            "identified by name, synonym, ChEMBL ID, CAS number, struct_id, or external IDs. "
            "Summarizes drug description, information, structural properties, regulatory status, pharmacologic "
            "actions, HUMAN gene-level targets (MoA and off-target), indications, "
            "contraindications, FAERS safety signals, and cross references, "
            "with text summaries and an optional interactive HTML viewer."
        ),
    ),
)
