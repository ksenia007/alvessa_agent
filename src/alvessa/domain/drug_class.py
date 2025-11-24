from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional


def canon_drug_key(name: str) -> str:
    """Normalize drug keys to lowercase alnum-only."""
    return "".join(ch for ch in (name or "").lower() if ch.isalnum())


def _unique_preserve(items: List[str]) -> List[str]:
    seen = set()
    cleaned: List[str] = []
    for item in items or []:
        norm = (item or "").strip()
        if not norm or norm.lower() in seen:
            continue
        cleaned.append(norm)
        seen.add(norm.lower())
    return cleaned


@dataclass
class DrugIdentifiers:
    """Core identifiers for a drug/small molecule."""

    name: str
    canon_key: Optional[str] = None
    synonyms: List[str] = field(default_factory=list)
    catalog_number: Optional[str] = None
    cas_number: Optional[str] = None
    chembl_id: Optional[str] = None


@dataclass
class Drug:
    """Aggregates drug-related data, mirroring the Gene class pattern."""

    identifiers: DrugIdentifiers
    smiles: Optional[str] = None
    research_areas: List[str] = field(default_factory=list)
    clinical_status: Optional[str] = None
    targets: List[str] = field(default_factory=list)  # gene/protein targets
    notes: List[str] = field(default_factory=list)
    text_summaries_from_tools: List[str] = field(default_factory=list)
    mentions: List[Dict[str, Any]] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    tools_run: List[str] = field(default_factory=list)

    # ------------------------------------------------------------------
    # Normalization
    # ------------------------------------------------------------------
    def normalize(self) -> None:
        if self.identifiers.name:
            self.identifiers.name = self.identifiers.name.strip()
        if self.identifiers.canon_key:
            self.identifiers.canon_key = canon_drug_key(self.identifiers.canon_key)
        self.identifiers.synonyms = _unique_preserve(self.identifiers.synonyms)
        self.targets = [t.upper().strip() for t in _unique_preserve(self.targets)]
        self.research_areas = _unique_preserve(self.research_areas)
        self.notes = _unique_preserve(self.notes)
        self.text_summaries_from_tools = _unique_preserve(self.text_summaries_from_tools)
        self.tools_run = _unique_preserve(self.tools_run)

    # ------------------------------------------------------------------
    # Setters / adders
    # ------------------------------------------------------------------
    def set_smiles(self, smiles: str) -> None:
        smiles = (smiles or "").strip()
        if smiles:
            self.smiles = smiles

    def add_research_area(self, area: str) -> None:
        area = (area or "").strip()
        if area and area not in self.research_areas:
            self.research_areas.append(area)

    def set_clinical_status(self, status: str) -> None:
        status = (status or "").strip()
        if status:
            self.clinical_status = status

    def add_target(self, target: str) -> None:
        target = (target or "").strip()
        if target and target.upper() not in [t.upper() for t in self.targets]:
            self.targets.append(target)

    def add_note(self, note: str) -> None:
        note = (note or "").strip()
        if note and note not in self.notes:
            self.notes.append(note)

    def update_text_summaries(self, summary: str) -> None:
        summary = (summary or "").strip()
        if summary and summary not in self.text_summaries_from_tools:
            self.text_summaries_from_tools.append(summary)

    def add_mention(self, mention: Dict[str, Any]) -> None:
        if mention and mention not in self.mentions:
            self.mentions.append(mention)

    def add_tool(self, tool_name: str) -> None:
        name = (tool_name or "").strip()
        if name and name not in self.tools_run:
            self.tools_run.append(name)

    def has_tool(self, tool_name: str) -> bool:
        return tool_name in self.tools_run

    def clear_tools(self) -> None:
        self.tools_run.clear()

    # ------------------------------------------------------------------
    # Serialization & presentation
    # ------------------------------------------------------------------
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), ensure_ascii=False)

    def summarize_text(self) -> str:
        self.normalize()
        lines: List[str] = []

        ids = self.identifiers
        id_bits = []
        if ids.name:
            id_bits.append(f"Name: {ids.name}")
        if ids.canon_key:
            id_bits.append(f"Key: {ids.canon_key}")
        if ids.cas_number:
            id_bits.append(f"CAS: {ids.cas_number}")
        if ids.catalog_number:
            id_bits.append(f"Catalog: {ids.catalog_number}")
        if ids.chembl_id:
            id_bits.append(f"ChEMBL: {ids.chembl_id}")
        if ids.synonyms:
            id_bits.append(f"Synonyms: {', '.join(ids.synonyms)}")
        if id_bits:
            lines.append(" | ".join(id_bits))

        if self.clinical_status:
            lines.append(f"Clinical status: {self.clinical_status}")
        if self.research_areas:
            lines.append(f"Research areas: {', '.join(self.research_areas)}")
        if self.targets:
            lines.append(f"Targets: {', '.join(self.targets)}")
        if self.smiles:
            lines.append(f"SMILES: {self.smiles}")
        if self.notes:
            lines.append(f"Notes: {' | '.join(self.notes)}")
        if self.text_summaries_from_tools:
            lines.extend(self.text_summaries_from_tools)

        return "\n".join([f"• {line}" for line in lines if line])

