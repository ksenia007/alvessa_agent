# src/alvessa/domain/drug_class.py
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional


def _unique_preserve(items: List[str]) -> List[str]:
    """
    Return a list with duplicates removed, preserving the first occurrence
    (case-insensitive comparison).
    """
    seen = set()
    out: List[str] = []
    for item in items or []:
        s = (item or "").strip()
        if not s:
            continue
        key = s.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append(s)
    return out


@dataclass
class DrugIdentifiers:
    """
    Minimal identifier bundle for a small molecule / drug.

    This class is focused on *identification and lookup*, not on
    downstream properties like activity or clinical status.

    All IDs are "hard" identifiers coming from external databases.
    No canonical key is stored here; merging logic lives in higher layers.
    """

    # Primary human-readable name (e.g., "neratinib").
    name: str

    # Hard identifiers from external DBs (optional).
    chembl_id: Optional[str] = None
    drugcentral_id: Optional[str] = None
    cas_number: Optional[str] = None
    catalog_number: Optional[str] = None

    # Any additional aliases or trade names.
    synonyms: List[str] = field(default_factory=list)

    def normalize(self) -> None:
        """Normalize identifiers in-place (trim, deduplicate synonyms)."""
        if self.name:
            self.name = self.name.strip()

        if self.chembl_id:
            self.chembl_id = str(self.chembl_id).strip()

        if self.drugcentral_id:
            self.drugcentral_id = str(self.drugcentral_id).strip()

        if self.cas_number:
            self.cas_number = str(self.cas_number).strip()

        if self.catalog_number:
            self.catalog_number = str(self.catalog_number).strip()

        self.synonyms = _unique_preserve(self.synonyms)


@dataclass
class Drug:
    """
    Compact representation of a drug / small molecule.

    Responsibilities:
      - Hold *identifiers* and *how/where* it was mentioned.
      - Be easy to match and merge across data sources.
      - Stay agnostic to downstream properties (activities, PK, etc.).
    """

    identifiers: DrugIdentifiers

    # Where/how this drug was mentioned in the conversation / text.
    # Example entry: {"text": "neratinib", "match_type": "exact", "source": "medchemexpress"}
    mentions: List[Dict[str, Any]] = field(default_factory=list)

    # Free-form metadata; reserved for future use (e.g., provenance flags).
    metadata: Dict[str, Any] = field(default_factory=dict)

    # Bookkeeping for which tools have enriched this object.
    tools_run: List[str] = field(default_factory=list)

    # ------------------------------------------------------------------
    # Normalization
    # ------------------------------------------------------------------
    def normalize(self) -> None:
        """Normalize identifiers, mentions, and tool list."""
        self.identifiers.normalize()
        self.tools_run = _unique_preserve(self.tools_run)

        # Normalize mention texts a bit (optional, non-destructive)
        normalized_mentions: List[Dict[str, Any]] = []
        for m in self.mentions or []:
            if not isinstance(m, dict):
                continue
            m_copy = dict(m)
            if "text" in m_copy and isinstance(m_copy["text"], str):
                m_copy["text"] = m_copy["text"].strip()
            normalized_mentions.append(m_copy)
        self.mentions = normalized_mentions

    # ------------------------------------------------------------------
    # Mutators / helpers
    # ------------------------------------------------------------------
    def add_mention(self, mention: Dict[str, Any]) -> None:
        """
        Add a mention record for this drug.

        A "mention" is a small dict capturing how this drug appeared
        in the text (e.g., substring, match type, source).
        """
        if not mention or not isinstance(mention, dict):
            return
        if mention not in self.mentions:
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
        """
        Produce a short bullet-style textual summary suitable for
        UI cards or export. Focuses only on identifiers and mentions.
        """
        self.normalize()
        ids = self.identifiers

        lines: List[str] = []

        id_bits: List[str] = []
        if ids.name:
            id_bits.append(f"Name: {ids.name}")
        if ids.chembl_id:
            id_bits.append(f"ChEMBL: {ids.chembl_id}")
        if ids.drugcentral_id:
            id_bits.append(f"DrugCentral: {ids.drugcentral_id}")
        if ids.cas_number:
            id_bits.append(f"CAS: {ids.cas_number}")
        if ids.catalog_number:
            id_bits.append(f"Catalog: {ids.catalog_number}")
        if ids.synonyms:
            id_bits.append(f"Synonyms: {', '.join(ids.synonyms)}")

        if id_bits:
            lines.append(" | ".join(id_bits))

        if self.mentions:
            texts = _unique_preserve(
                [str(m.get("text", "")).strip() for m in self.mentions if isinstance(m, dict)]
            )
            if texts:
                lines.append(f"Mentions in text: {', '.join(texts)}")

        if self.tools_run:
            lines.append(f"Tools: {', '.join(self.tools_run)}")

        return "\n".join([f"• {line}" for line in lines if line])
