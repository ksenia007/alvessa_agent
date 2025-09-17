"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-09-01
Updated: 2025-09-15

Description: 
Defining specialized data-holding components for the Gene class to promote
a cleaner, more modular design using composition.
"""

from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple, Any

# --- Helper Functions ---

def _upper_unique(items: List[str]) -> List[str]:
    """Returns a sorted list of unique, uppercase strings from the input list."""
    return sorted({s.upper() for s in items if s})

# --- Component Dataclasses ---

@dataclass
class GeneIdentifiers:
    """A container for core gene identifiers."""
    symbol: str
    entrez_id: Optional[str] = None
    uniprot_id: Optional[str] = None
    ensembl_id: Optional[str] = None
    aliases: List[str] = field(default_factory=list)
    gene_type: Optional[str] = None

@dataclass
class GenomicLocation:
    """A container for genomic coordinates."""
    chrom: str
    start: int
    end: int
    strand: Optional[str] = None
    build: str = "GRCh38"

@dataclass
class FunctionalAnnotations:
    """A container for functional annotations like GO terms, pathways, and diseases."""
    go_annotations: List[str] = field(default_factory=list)
    functions: List[str] = field(default_factory=list)
    diseases: List[str] = field(default_factory=list)
    pathways: List[str] = field(default_factory=list)
    cell_localizations: List[str] = field(default_factory=list)
    predicted_go: List[str] = field(default_factory=list)
    predicted_disease: List[str] = field(default_factory=list)

@dataclass
class InteractionProfile:
    """A container for gene-gene interaction data."""
    # Key: experiment type, Val: list of partner gene symbols
    human_interactions: Dict[str, List[str]] = field(default_factory=dict)
    nonhuman_interactions: Dict[str, List[str]] = field(default_factory=dict)
    go_enrichment: Dict[str, Dict[str, List[str]]] = field(default_factory=dict)

    def all_interaction_partners(self) -> List[str]:
        """Flattens and returns a unique, sorted list of human interaction partners."""
        partners = {p for partners_list in self.human_interactions.values() for p in partners_list}
        return sorted(list(partners))

    def add_enrichment(self, source: str, category: str, terms: List[str]) -> None:
        """Record GO enrichment terms by source (e.g. old_go, pan_go) and category."""
        source = (source or "").strip()
        category = (category or "").strip()
        if not source or not category:
            return

        cleaned = []
        seen = set()
        for term in terms or []:
            term_norm = (term or "").strip()
            if not term_norm or term_norm in seen:
                continue
            cleaned.append(term_norm)
            seen.add(term_norm)

        bucket = self.go_enrichment.setdefault(source, {})
        bucket[category] = cleaned

    def get_enrichment(self, source: str | None = None) -> Dict[str, Dict[str, List[str]]]:
        """Return recorded GO enrichment grouped by source/category."""
        if source is None:
            return self.go_enrichment
        return {source: dict(self.go_enrichment.get(source, {}))}

@dataclass
class TranscriptomeProfile:
    """A container for transcript/isoform data and their derived stats."""
    # transcripts[transcript_id] = { "chrom": str, "start": int, "end": int, "strand": str, "exons": List[Tuple[int,int]] }
    transcripts: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    isoforms: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    # Derived stats
    avg_exons_per_transcript: Optional[float] = None
    transcript_count: Optional[int] = None
    total_exon_count: Optional[int] = None
    median_transcript_span_bp: Optional[float] = None
    max_transcript_span_bp: Optional[int] = None

    def compute_stats(self) -> None:
        """Calculates transcript-based statistics and updates instance attributes."""
        n_tx = len(self.transcripts)
        self.transcript_count = n_tx
        if n_tx == 0:
            self.avg_exons_per_transcript = None
            self.total_exon_count = None
            self.median_transcript_span_bp = None
            self.max_transcript_span_bp = None
            return

        total_exons = 0
        spans: List[int] = []

        for tx in self.transcripts.values():
            if tx.get("exons") is not None:
                total_exons += len(tx["exons"])
            elif tx.get("n_exons") is not None:
                total_exons += int(tx["n_exons"])

            if tx.get("start") is not None and tx.get("end") is not None:
                spans.append(int(tx["end"]) - int(tx["start"]) + 1)

        self.total_exon_count = total_exons
        self.avg_exons_per_transcript = (total_exons / n_tx) if n_tx else None

        if spans:
            spans.sort()
            mid = len(spans) // 2
            self.median_transcript_span_bp = float(spans[mid]) if len(spans) % 2 else (spans[mid-1] + spans[mid]) / 2.0
            self.max_transcript_span_bp = spans[-1]
        else:
            self.median_transcript_span_bp = None
            self.max_transcript_span_bp = None

@dataclass
class GeneGWASTraitHit:
    """Represents a single significant GWAS hit linking a trait and a variant."""
    trait: str
    rsid: str
    p_value: Optional[str] = None
    risk_score: Optional[str] = None

@dataclass
class GeneGWASProfile:
    """Aggregates GWAS summary data for a gene."""
    found: bool = False
    total_associations: int = 0
    significant_associations: int = 0
    total_studies: int = 0
    variant_count: int = 0
    trait_link_count: int = 0
    p_value_threshold: Optional[str] = None
    top_traits: List[GeneGWASTraitHit] = field(default_factory=list)
    top_variants: List[str] = field(default_factory=list)
    related_genes: List[str] = field(default_factory=list)
    affected_proteins: List[str] = field(default_factory=list)
