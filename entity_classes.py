"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-09-01
Updated: 2025-09-02


Description: 

Definint entity classes for Gene, Variant, Trait"""

from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple, Any
import json

def _upper_unique(items: List[str]) -> List[str]:
    return sorted({s.upper() for s in items if s})

def _dedup(items: List[str]) -> List[str]:
    seen = set(); out = []
    for x in items:
        if x not in seen:
            seen.add(x); out.append(x)
    return out

def canon_gene_key(x: str) -> str:
    # make sure that gene keys are in a consistent format
    return (x or "").strip().upper()

@dataclass
class Gene:
    # ---- Core identifiers ----
    symbol: str
    entrez_id: Optional[str] = None
    uniprot_id: Optional[str] = None
    aliases: List[str] = field(default_factory=list)
    gene_type: Optional[str] = None
    build: Optional[str] = "GRCh38"  

    # ---- Gene-level genomic location ----
    chrom: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None
    strand: Optional[str] = None

    # ---- Links ----
    variants: List[str] = field(default_factory=list)
    traits: List[str] = field(default_factory=list)

    # ---- Interactions (BioGRID later) ----
    # Key: experiment type (as in BioGRID "Experimental System", e.g. "Two-hybrid", "Co-IP",
    #      "Affinity Capture-MS", "PCA", "Dosage Rescue", etc.)
    # Val: list of partner gene IDs/symbols
    interactions_by_exp: Dict[str, List[str]] = field(default_factory=dict)
    nonhuman_interaction_by_exp: Dict[str, List[str]] = field(default_factory=dict)


    # ---- Functional annotations ----
    go_annotations: List[str] = field(default_factory=list)
    functions: List[str] = field(default_factory=list)
    diseases: List[str] = field(default_factory=list)
    pathways: List[str] = field(default_factory=list)
    cell_localizations: Optional[str] = None

    # ---- Transcripts ----
    # transcripts[transcript_id] = { "chrom": str, "start": int, "end": int, "strand": str,
    #                                "exons": List[Tuple[int,int]] }
    transcripts: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    # ---- Derived structure stats ----
    avg_exons_per_transcript: Optional[float] = None
    transcript_count: Optional[int] = None
    total_exon_count: Optional[int] = None
    median_transcript_span_bp: Optional[float] = None
    max_transcript_span_bp: Optional[int] = None
    
    # ---- Tools ----
    tools_run: List[str] = field(default_factory=list)


    # ------------ Normalization ------------
    # Make sure chrom, aliases etc are in a consistent format, variants and traits are unique
    def normalize(self) -> None:
        print('IN NORMALIZE', self.symbol)
        if self.symbol: self.symbol = self.symbol.upper()
        self.aliases = _upper_unique(self.aliases)
        if self.chrom: self.chrom = self.chrom.replace("chr", "").upper()
        if self.strand: self.strand = self.strand.strip()
        self.variants = _dedup(self.variants)
        self.traits   = _dedup(self.traits)
        # normalize interaction partner casing; keep experiment type as-is (BioGRID labels)
        self.interactions_by_exp = {
            exp: _upper_unique(partners)
            for exp, partners in sorted(self.interactions_by_exp.items())
        }
        # normalize transcript records
        self.transcripts = {k: self._norm_tx(v) for k, v in sorted(self.transcripts.items())}

    @staticmethod
    def _norm_tx(tx: Dict[str, Any]) -> Dict[str, Any]:
        # normalize a single transcript record
        d = dict(tx)
        if d.get("chrom"):  d["chrom"] = d["chrom"].replace("chr", "").upper()
        if d.get("strand"): d["strand"] = d["strand"].strip()
        exons = d.get("exons", [])
        d["exons"] = sorted([(int(a), int(b)) for a, b in exons], key=lambda x: (x[0], x[1]))
        for k in ("start", "end"):
            if d.get(k) is not None: d[k] = int(d[k])
        return d

    # ------------ Linking to other classes ------------
    def link_variant(self, variant_id: str) -> None:
        if variant_id not in self.variants:
            self.variants.append(variant_id)

    def link_trait(self, trait_id: str) -> None:
        if trait_id not in self.traits:
            self.traits.append(trait_id)
            
    def set_gene_type(self, gene_type: str) -> None:
        gene_type = (gene_type or "").strip()
        if gene_type:
            self.gene_type = gene_type
            
    def set_gene_ids(self, entrez_id: Optional[str] = None, uniprot_id: Optional[str] = None) -> None:
        if entrez_id:
            self.entrez_id = str(entrez_id).strip()
        if uniprot_id:
            self.uniprot_id = str(uniprot_id).strip()
            
    def set_chrom_location(self, chrom: str, gene_span: Tuple[Optional[int], Optional[int]], strand: Optional[str] = None) -> None:
        chrom = (chrom or "").replace("chr", "").upper().strip()
        if chrom:
            self.chrom = chrom
        start, end = gene_span
        if start is not None:
            self.start = int(start)
        if end is not None:
            self.end = int(end)
        if strand:
            self.strand = strand.strip()

    # ------------ Interaction functions ------------
    def add_interaction(self, experiment_type: str, partner_gene: str) -> None:
        """Add a partner gene under a given experimental system."""
        if not experiment_type: return
        partner_gene = (partner_gene or "").upper()
        if not partner_gene: return
        self.interactions_by_exp.setdefault(experiment_type, [])
        if partner_gene not in self.interactions_by_exp[experiment_type]:
            self.interactions_by_exp[experiment_type].append(partner_gene)
    
    def add_many_interactions(self, experiment_type: str, partner_genes: List[str]) -> None:
        """Add multiple partner genes under a given experimental system."""
        for pg in partner_genes:
            self.add_interaction(experiment_type, pg)

    def get_interactions(self, experiment_type: Optional[str] = None) -> Dict[str, List[str]] | List[str]:
        """Get all interactions or only a specific experiment type."""
        if experiment_type is None:
            return self.interactions_by_exp
        return self.interactions_by_exp.get(experiment_type, [])

    def all_interaction_partners(self) -> List[str]:
        """Flatten list of unique partners across all experiment types."""
        partners: List[str] = []
        for lst in self.interactions_by_exp.values():
            partners.extend(lst)
        return _upper_unique(partners)
    
    def has_interactions_collected(self) -> bool:
        """Check if any interactions have been recorded."""
        return any(self.interactions_by_exp.values())
    
    # we also have non-human interactions, we save them separately and always with experiment type
    def add_nonhuman_interaction(self, experiment_type: str, partner_gene: str) -> None:
        """Add a non-human partner gene under a given experimental system in special nonhuman field."""
        if not experiment_type: return
        partner_gene = (partner_gene or "").upper()
        if not partner_gene: return
        self.nonhuman_interaction_by_exp.setdefault(experiment_type, [])
        if partner_gene not in self.nonhuman_interaction_by_exp[experiment_type]:
            self.nonhuman_interaction_by_exp[experiment_type].append(partner_gene)
    
    def add_many_nonhuman_interactions(self, experiment_type: str, partner_genes: List[str]) -> None:
        """Add multiple non-human partner genes under a given experimental system."""
        for pg in partner_genes:
            self.add_nonhuman_interaction(experiment_type, pg)

    # ------------ GO & other functions ------------
    def add_go_terms(self, label: str) -> None:
        label = (label or "").strip()
        if label and label not in self.go_annotations:
            self.go_annotations.append(label)
            
    def add_many_go_terms(self, labels: List[str]) -> None:
        for label in labels:
            self.add_go_terms(label)
    
    def add_disease_label(self, label: str) -> None:
        label = (label or "").strip()
        if label and label not in self.diseases:
            self.diseases.append(label)
            
    def add_many_diseases(self, labels: List[str]) -> None:
        for label in labels:
            self.add_disease_label(label)

    def add_function_label(self, label: str) -> None:
        label = (label or "").strip()
        if label and label not in self.functions:
            self.functions.append(label)
    
    def add_many_functions(self, labels: List[str]) -> None:
        for label in labels:
            self.add_function_label(label)
            
    def add_pathway_label(self, label: str) -> None:
        label = (label or "").strip()
        if label and label not in self.pathways:
            self.pathways.append(label)
        
    def add_many_pathways(self, labels: List[str]) -> None:
        for label in labels:
            self.add_pathway_label(label)
    
    def already_called_pathways(self) -> bool:
        return bool(len(self.pathways)>1)
    
    def set_cell_localization(self, label: str) -> None:
        label = (label or "").strip()
        if label:
            self.cell_localizations = label
            

    # ------------ Transcripts ------------
    def add_transcript(self, transcript_id: str, n_exons: int) -> None:
        self.transcripts[transcript_id] = {
            "n_exons": n_exons,
        }

    def get_transcript(self, transcript_id: str) -> Optional[Dict[str, Any]]:
        return self.transcripts.get(transcript_id)
    
    def get_n_transcripts(self) -> int:
        return len(self.transcripts)

    # ------------ Structure stats ------------
    def compute_structure_stats(self) -> None:
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
            exons = tx.get("exons", [])
            total_exons += len(exons)
            if tx.get("start") is not None and tx.get("end") is not None:
                spans.append(int(tx["end"]) - int(tx["start"]) + 1)

        self.total_exon_count = total_exons
        self.avg_exons_per_transcript = total_exons / n_tx if n_tx else None
        spans.sort()
        if spans:
            mid = len(spans) // 2
            self.median_transcript_span_bp = float(spans[mid]) if len(spans) % 2 else (spans[mid-1] + spans[mid]) / 2.0
            self.max_transcript_span_bp = spans[-1]
        else:
            self.median_transcript_span_bp = None
            self.max_transcript_span_bp = None

    # ------------ Record tools that were run ------------
    def add_tool(self, tool_name: str) -> None:
        """Record that a tool was run for this gene."""
        name = tool_name.strip()
        if name and name not in self.tools_run:
            self.tools_run.append(name)

    def has_tool(self, tool_name: str) -> bool:
        """Check if a tool was already run."""
        return tool_name in self.tools_run

    def clear_tools(self) -> None:
        """Reset tool history"""
        self.tools_run.clear()


    # ------------ Serialization ------------
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "Gene":
        g = cls(**d)
        g.normalize()
        g.compute_structure_stats()
        return g

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), ensure_ascii=False)

    @classmethod
    def from_json(cls, s: str) -> "Gene":
        return cls.from_dict(json.loads(s))

