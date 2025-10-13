from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional, Tuple

from .variant_class import Variant
from .gene_components import (
    GeneIdentifiers,
    GenomicLocation,
    FunctionalAnnotations,
    InteractionProfile,
    MSigDBAnnotations,
    OMIMAnnotations,
    OpenTargetAnnotations,
    TranscriptomeProfile,
    GeneGWASProfile,
    GeneGWASTraitHit,
    _upper_unique,
)


def canon_gene_key(symbol: str) -> str:
    """Normalize gene keys to uppercase with trimmed whitespace."""
    return (symbol or "").strip().upper()


@dataclass
class Gene:
    """Aggregates gene-related data using composable component classes."""

    identifiers: GeneIdentifiers
    location: Optional[GenomicLocation] = None
    annotations: FunctionalAnnotations = field(default_factory=FunctionalAnnotations)
    interactions: InteractionProfile = field(default_factory=InteractionProfile)
    msigDB: MSigDBAnnotations = field(default_factory=MSigDBAnnotations)
    omim: OMIMAnnotations = field(default_factory=OMIMAnnotations)
    open_targets: OpenTargetAnnotations = field(default_factory=OpenTargetAnnotations)
    transcriptome: TranscriptomeProfile = field(default_factory=TranscriptomeProfile)
    gwas_profile: Optional[GeneGWASProfile] = None

    variants: Dict[str, Variant] = field(default_factory=dict)
    traits: List[str] = field(default_factory=list)

    binding_peaks: Dict[str, Any] = field(default_factory=dict)
    mirna_targets: Dict[str, List[str]] = field(default_factory=dict)
    text_summaries_from_tools: List[str] = field(default_factory=list)
    tools_run: List[str] = field(default_factory=list)

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------
    @property
    def symbol(self) -> str:
        return self.identifiers.symbol

    @property
    def gene_type(self) -> Optional[str]:
        return self.identifiers.gene_type

    @property
    def aliases(self) -> List[str]:
        return self.identifiers.aliases

    @property
    def entrez_id(self) -> Optional[str]:
        return self.identifiers.entrez_id

    @property
    def uniprot_id(self) -> Optional[str]:
        return self.identifiers.uniprot_id

    @property
    def ensembl_id(self) -> Optional[str]:
        return self.identifiers.ensembl_id

    @property
    def chrom(self) -> Optional[str]:
        return self.location.chrom if self.location else None

    @property
    def start(self) -> Optional[int]:
        return self.location.start if self.location else None

    @property
    def end(self) -> Optional[int]:
        return self.location.end if self.location else None

    @property
    def strand(self) -> Optional[str]:
        return self.location.strand if self.location else None

    # ------------------------------------------------------------------
    # Normalization
    # ------------------------------------------------------------------
    def normalize(self) -> None:
        if self.identifiers.symbol:
            self.identifiers.symbol = self.identifiers.symbol.upper()
        self.identifiers.aliases = _upper_unique(self.identifiers.aliases)

        if self.location and self.location.chrom:
            self.location.chrom = self.location.chrom.replace("chr", "").upper()

        self.interactions.human_interactions = {
            exp: _upper_unique(partners)
            for exp, partners in sorted(self.interactions.human_interactions.items())
        }
        self.interactions.nonhuman_interactions = {
            exp: _upper_unique(partners)
            for exp, partners in sorted(self.interactions.nonhuman_interactions.items())
        }
        self.traits = sorted(set(self.traits))

    # ------------------------------------------------------------------
    # Identifier helpers
    # ------------------------------------------------------------------
    def set_gene_type(self, gene_type: str) -> None:
        gene_type = (gene_type or "").strip()
        if gene_type:
            self.identifiers.gene_type = gene_type

    def set_gene_ids(
        self,
        entrez_id: Optional[str] = None,
        uniprot_id: Optional[str] = None,
        ensemble_id: Optional[str] = None,
    ) -> None:
        if entrez_id:
            self.identifiers.entrez_id = str(entrez_id).strip()
        if uniprot_id:
            self.identifiers.uniprot_id = str(uniprot_id).strip()
        if ensemble_id:
            self.identifiers.ensembl_id = str(ensemble_id).strip()

    def set_chrom_location(
        self,
        chrom: str,
        gene_span: Tuple[Optional[int], Optional[int]],
        strand: Optional[str] = None,
    ) -> None:
        chrom_norm = (chrom or "").replace("chr", "").upper().strip()
        start, end = gene_span
        if not any([chrom_norm, start, end, strand]):
            return
        build = self.location.build if self.location else "GRCh38"
        self.location = GenomicLocation(
            chrom=chrom_norm or (self.location.chrom if self.location else ""),
            start=int(start) if start is not None else (self.location.start if self.location else 0),
            end=int(end) if end is not None else (self.location.end if self.location else 0),
            strand=strand or (self.location.strand if self.location else None),
            build=build,
        )

    def get_location(self) -> Optional[Tuple[str, Optional[int], Optional[int], Optional[str]]]:
        if not self.location:
            return None
        return (
            self.location.chrom,
            self.location.start,
            self.location.end,
            self.location.strand,
        )

    def is_positive_strand(self) -> Optional[bool]:
        if self.strand == "+":
            return True
        if self.strand == "-":
            return False
        return None

    # ------------------------------------------------------------------
    # Functional annotations
    # ------------------------------------------------------------------
    def add_go_terms(self, label: str) -> None:
        label = (label or "").strip()
        if label and label not in self.annotations.go_annotations:
            self.annotations.go_annotations.append(label)

    def add_many_go_terms(self, labels: List[str]) -> None:
        for label in labels:
            self.add_go_terms(label)

    def add_function_label(self, label: str) -> None:
        label = (label or "").strip()
        if label and label not in self.annotations.functions:
            self.annotations.functions.append(label)

    def add_many_functions(self, labels: List[str]) -> None:
        for label in labels:
            self.add_function_label(label)

    def add_disease_label(self, label: str) -> None:
        label = (label or "").strip()
        if label and label not in self.annotations.diseases:
            self.annotations.diseases.append(label)

    def add_many_diseases(self, labels: List[str]) -> None:
        for label in labels:
            self.add_disease_label(label)

    def add_predicted_go(self, label: str) -> None:
        label = (label or "").strip()
        if label and label not in self.annotations.predicted_go:
            self.annotations.predicted_go.append(label)

    def add_many_predicted_go(self, labels: List[str]) -> None:
        for label in labels:
            self.add_predicted_go(label)

    def add_predicted_disease(self, label: str) -> None:
        label = (label or "").strip()
        if label and label not in self.annotations.predicted_disease:
            self.annotations.predicted_disease.append(label)

    def add_many_predicted_diseases(self, labels: List[str]) -> None:
        for label in labels:
            self.add_predicted_disease(label)

    def add_pathway(self, label: str) -> None:
        label = (label or "").strip()
        if label and label not in self.annotations.pathways:
            self.annotations.pathways.append(label)

    def add_many_pathways(self, labels: List[str]) -> None:
        for label in labels:
            self.add_pathway(label)

    def already_called_pathways(self) -> bool:
        return len(self.annotations.pathways) > 1

    # ------------------------------------------------------------------
    # Interactions
    # ------------------------------------------------------------------
    def add_interaction(self, experiment_type: str, partner_gene: str) -> None:
        partner_gene = (partner_gene or "").upper().strip()
        if not experiment_type or not partner_gene:
            return
        partners = self.interactions.human_interactions.setdefault(experiment_type, [])
        if partner_gene not in partners:
            partners.append(partner_gene)

    def add_many_interactions(self, experiment_type: str, partner_genes: List[str]) -> None:
        for partner in partner_genes:
            self.add_interaction(experiment_type, partner)

    def add_nonhuman_interaction(self, experiment_type: str, partner_gene: str) -> None:
        partner_gene = (partner_gene or "").upper().strip()
        if not experiment_type or not partner_gene:
            return
        partners = self.interactions.nonhuman_interactions.setdefault(experiment_type, [])
        if partner_gene not in partners:
            partners.append(partner_gene)

    def add_many_nonhuman_interactions(self, experiment_type: str, partner_genes: List[str]) -> None:
        for partner in partner_genes:
            self.add_nonhuman_interaction(experiment_type, partner)

    def add_interactors_enriched_go_terms(
        self,
        source: str,
        category: str,
        enriched_terms: List[str],
    ) -> None:
        """Record GO term enrichment calculated over interaction partners."""
        self.interactions.add_enrichment(source, category, enriched_terms)

    def get_interactors_enriched_go_terms(self) -> Dict[str, Dict[str, List[str]]]:
        """Return GO enrichment summaries grouped by source and category."""
        return self.interactions.get_enrichment()

    def get_interactions(self, experiment_type: Optional[str] = None) -> Dict[str, List[str]] | List[str]:
        if experiment_type is None:
            return self.interactions.human_interactions
        return self.interactions.human_interactions.get(experiment_type, [])

    def all_interaction_partners(self) -> List[str]:
        return self.interactions.all_interaction_partners()

    def has_interactions_collected(self) -> bool:
        return any(self.interactions.human_interactions.values())
        
    # ------------------------------------------------------------------
    # MSigDB
    # ------------------------------------------------------------------

    def add_msigdb_annotated_terms(self, category: str, annotation_terms: List[str]) -> None:
        self.msigDB.add_msigdb_annotation(category, annotation_terms)
        
    # ------------------------------------------------------------------
    # OMIM
    # ------------------------------------------------------------------

    def add_omim_annotated_term(self, phenotype: str) -> None:
        phenotype = (phenotype or "").strip()
        if phenotype and phenotype not in self.omim.phenotype_annotations:
            self.omim.phenotype_annotations.append(phenotype)
            
    def add_many_omim_annotated_terms(self, phenotypes: str) -> None:
        for phenotype in phenotypes:
            self.add_omim_annotated_term(phenotype)

    # ------------------------------------------------------------------
    # Open Targets
    # ------------------------------------------------------------------

    def add_direct_disease_association(self, disease: str) -> None:
        disease = (disease or "").strip()
        if disease and disease not in self.open_targets.disease_annotations:
            self.open_targets.disease_annotations.append(disease)

    def add_many_direct_disease_associations(self, diseases: List[str]) -> None:
        for disease in diseases:
            self.add_direct_disease_association(disease)

    def add_tissue_expression(self, tissue: str, zscore: int) -> None:
        tissue = (tissue or "").strip()
        if tissue and tissue not in self.open_targets.tissue_specific_expression:
            self.open_targets.tissue_specific_expression[tissue] = zscore

    def add_many_tissues_expression(self, tissues: Dict[str, Any]) -> None:
        for tissue, zscore in tissues.items():
            self.add_tissue_expression(tissue, zscore)

    def add_essentiality(self, is_essential) -> None:
        self.open_targets.is_essential = is_essential

    def add_constraint(self, score, constraint_type) -> None:
        if constraint_type not in self.open_targets.genetic_constraint:
            self.open_targets.genetic_constraint[constraint_type] = score

    # ------------------------------------------------------------------
    # Text summaries & traits
    # ------------------------------------------------------------------
    def update_text_summaries(self, summary: str) -> None:
        summary = (summary or "").strip()
        if summary and summary not in self.text_summaries_from_tools:
            self.text_summaries_from_tools.append(summary)

    def link_trait(self, trait_id: str) -> None:
        trait_id = (trait_id or "").strip()
        if trait_id and trait_id not in self.traits:
            self.traits.append(trait_id)

    def get_all_traits(self) -> List[str]:
        return list(self.traits)

    # ------------------------------------------------------------------
    # Variants & tools
    # ------------------------------------------------------------------
    def link_variant(self, variant: Variant) -> None:
        if variant.rsID not in self.variants:
            self.variants[variant.rsID] = variant

    def add_tool(self, tool_name: str) -> None:
        name = (tool_name or "").strip()
        if name and name not in self.tools_run:
            self.tools_run.append(name)

    def has_tool(self, tool_name: str) -> bool:
        return tool_name in self.tools_run

    def clear_tools(self) -> None:
        self.tools_run.clear()

    # ------------------------------------------------------------------
    # Transcriptome helpers
    # ------------------------------------------------------------------
    def add_transcript(self, transcript_id: str, n_exons: int) -> None:
        if not transcript_id:
            return
        self.transcriptome.transcripts[transcript_id] = {"n_exons": n_exons}

    def get_transcript(self, transcript_id: str) -> Optional[Dict[str, Any]]:
        return self.transcriptome.transcripts.get(transcript_id)

    def get_n_transcripts(self) -> int:
        return len(self.transcriptome.transcripts)

    def set_isoform_localizations(
        self,
        isoform_map: Dict[str, Dict[str, Any]],
        *,
        replace: bool = False,
    ) -> None:
        if replace:
            self.transcriptome.isoforms = {}

        for iso_id, rec in (isoform_map or {}).items():
            if not iso_id:
                continue
            dst = self.transcriptome.isoforms.setdefault(
                iso_id,
                {"name": "", "aliases": [], "status": None, "locations": [], "notes": []},
            )
            name = rec.get("name") or ""
            if name and not dst.get("name"):
                dst["name"] = name
            status = rec.get("status")
            if status and not dst.get("status"):
                dst["status"] = status
            for alias in rec.get("aliases", []) or []:
                if alias and alias not in dst["aliases"]:
                    dst["aliases"].append(alias)
            for loc in rec.get("locations", []) or []:
                if loc and loc not in dst["locations"]:
                    dst["locations"].append(loc)
            for note in rec.get("notes", []) or []:
                if note and note not in dst["notes"]:
                    dst["notes"].append(note)

    # ------------------------------------------------------------------
    # GWAS helpers
    # ------------------------------------------------------------------
    def set_gwas_profile(self, profile: Optional[GeneGWASProfile]) -> None:
        self.gwas_profile = profile

    def clear_gwas_profile(self) -> None:
        self.gwas_profile = None

    def has_gwas_associations(self) -> bool:
        return bool(self.gwas_profile and self.gwas_profile.total_associations)

    def get_all_gwas_associations(self) -> int:
        if not self.gwas_profile:
            return 0
        return self.gwas_profile.total_associations

    def get_all_gwas_info(self) -> Optional[Dict[str, Any]]:
        if not self.gwas_profile:
            return None
        return asdict(self.gwas_profile)

    # ------------------------------------------------------------------
    # Binding peaks & miRNA targets
    # ------------------------------------------------------------------
    def add_binding_peaks(self, peak_info: Dict[str, Any]) -> None:
        self.binding_peaks = peak_info or {}

    def get_n_unique_binding_peaks(self) -> int:
        return self.binding_peaks.get("n_peaks", 0)

    def get_n_binders(self) -> Any:
        return self.binding_peaks.get("n_unique_binders", 0)

    def get_all_binding_tfs(self) -> List[str]:
        return self.binding_peaks.get("unique_binders", [])

    def add_miRNA_targets(self, miRNA: Dict[str, Any]) -> None:
        self.mirna_targets = miRNA or {}

    # ------------------------------------------------------------------
    # Serialization & presentation
    # ------------------------------------------------------------------
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), ensure_ascii=False)

    def summarize_text(
        self,
        *,
        include_go: bool = False,
        include_pathways: bool = True,
        include_interactions: bool =True
    ) -> str:
        presenter = GeneTextPresenter(self)
        return presenter.summarize(include_go=include_go, include_pathways=include_pathways, include_interactions = include_interactions)


class GeneTextPresenter:
    """Formats Gene objects into concise textual summaries."""

    def __init__(self, gene: Gene):
        self.gene = gene

    def _kv(self, key: str, value: Any) -> Optional[str]:
        if value is None or (isinstance(value, str) and not value.strip()):
            return None
        return f"{key}: {value}"

    def _join(self, items: List[Optional[str]], sep: str = "; ") -> str:
        return sep.join([x for x in items if x])

    def _clean(self, values: List[Any]) -> List[str]:
        return [str(v).strip() for v in values if str(v).strip()]

    def summarize(self, include_go: bool = False, include_pathways: bool = True, include_interactions: bool = True) -> str:
        self.gene.normalize()
        self.gene.transcriptome.compute_stats()

        lines: List[str] = []

        ids = self.gene.identifiers
        id_line = self._join(
            [
                self._kv("Symbol", ids.symbol),
                self._kv("Entrez", ids.entrez_id),
                self._kv("UniProt", ids.uniprot_id),
                self._kv("Gene type", ids.gene_type),
            ],
            sep=" | ",
        )
        if id_line:
            lines.append(id_line)

        if self.gene.location:
            loc = self.gene.location
            span = None
            if loc.start is not None and loc.end is not None:
                span = f"{loc.start:,}-{loc.end:,}"
            loc_line = self._join([f"chr{loc.chrom}", span, loc.strand], sep=" ")
            if loc_line:
                lines.append(loc_line)

        txp = self.gene.transcriptome
        if txp.transcript_count is not None:
            tx_bits = [f"n_transcripts={txp.transcript_count}:"]
            # add transcript name + exon count
            for tx_id, tx_info in txp.transcripts.items():
                n_exons = tx_info.get("n_exons")
                if n_exons is not None:
                    tx_bits.append(f"{tx_id}({n_exons} exons)")
                else:
                    tx_bits.append(tx_id)
            if txp.max_transcript_span_bp is not None:
                tx_bits.append(f"max_span={txp.max_transcript_span_bp:,} bp")
            lines.append(self._kv("Transcripts", self._join(tx_bits)))

        if include_pathways and self.gene.annotations.pathways:
            lines.append(
                self._kv(
                    "Pathways",
                    self._join(self._clean(self.gene.annotations.pathways)),
                )
            )

        if include_go and self.gene.annotations.go_annotations:
            lines.append(
                self._kv(
                    "GO",
                    self._join(self._clean(self.gene.annotations.go_annotations)),
                )
            )
        
        if include_interactions and self.gene.has_interactions_collected():
            # include summary + LIST of the partners, comma-separated by human and non-human
            n_partners = len(self.gene.all_interaction_partners())
            n_exp_types = len(self.gene.interactions.human_interactions)
            n_nonhuman_exp_types = len(self.gene.interactions.nonhuman_interactions)
            # add string for human first, grouped by experiment type
            for exp_type, partners in self.gene.interactions.human_interactions.items():
                if partners:
                    lines.append(
                        self._kv(
                            f"Interactions from human data ({exp_type})",
                            self._join(self._clean(partners), sep=", "),
                        )
                    )
            # add string for non-human, grouped by experiment type
            for exp_type, partners in self.gene.interactions.nonhuman_interactions.items():
                if partners:
                    lines.append(
                        self._kv(
                            f"Interactions from additional (non-human) data: ({exp_type})",
                            self._join(self._clean(partners), sep=", "),
                        )
                    )
                    

        if self.gene.text_summaries_from_tools:
            lines.extend(self.gene.text_summaries_from_tools)

        return "\n".join([f"• {line}" for line in lines if line])
