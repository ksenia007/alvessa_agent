"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2024-09-01
Updated: 2025-09-02


Description: 

Defining entity classes for Gene, Variant, Trait"""

from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple, Any
import json
import requests





@dataclass
class Variant:
    rsID: Optional[str]
    loc_by_build: Optional[Dict[str, Any]] = field(default_factory=dict) # e.g., {"GrCh38": chrom, pos, ref, alt}
    organism: Optional[str] = None
    text_summary: Optional[str] = None   # Short text description of the variant, to be propagated to the LLM
    traits: List[str] = field(default_factory=list)  # List of associated traits/diseases across genes
    genes_related_to: List[str] = field(default_factory=list)
    per_gene_traits: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    per_gene_context: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    tools_run: List[str] = field(default_factory=list)  # List of tools that have been run for this variant
    af_freqs: Optional[List[Dict[str, Any]]] = field(default_factory=list)  # List of allele frequencies from various populations/databases
    variant_summaries: Optional[List[str]] = field(default_factory=list)  # List of textual summaries from various tools
    functional_predictions: Optional[Dict[str, Dict[str, Any]]] = field(default_factory=dict)  # e.g., {gene: {'AM': [score], 'CADD': score, 'expecto': score}}
    drug_response_effects: List[str] = field(default_factory=list) # variant effect on drug responses from OpenTargets
    
    def update_text_summaries(self, new_summary: str):
        """Add a new text summary to the variant"""
        if new_summary and new_summary not in self.variant_summaries:
            self.variant_summaries.append(new_summary)
            
    def return_full_summary(self) -> str:
        # add functional predictions summary first
        print(self.variant_summaries)
        func_summary = ""
        if self.functional_predictions:
            func_preds = []
            for gene, tools in self.functional_predictions.items():
                for tool, scores in tools.items():
                    score_str = ', '.join(map(str, scores))
                    func_preds.append(f"{tool}: {score_str}")
                func_summary += f"Functional predictions for {gene}: " + ", ".join(func_preds)
        
        if func_summary:
            self.variant_summaries.insert(0, func_summary)
            
        
        # add per_gene_context info
        for gene, context in self.per_gene_context.items():
            ctx = context.get('context', '')
            var_cat = context.get('variant_category', '')
            context_str = f"Gene: {gene}, Category: {var_cat}"
            self.variant_summaries.insert(0, context_str)  
            
        # add info about organism, location and per_gene_context (by gene) info
        if self.organism:
            self.variant_summaries.insert(0, f"Organism: {self.organism}")
        
        for build, loc in self.loc_by_build.items():

            chrom = loc.get('chrom', '')
            pos = loc.get('pos', '')
            ref = ','.join(loc.get('ref', [])) if loc.get('ref') else ''
            alt = ','.join(loc.get('alt', [])) if loc.get('alt') else ''
            loc_str = f"Location ({build}): chr{chrom}:{pos}"
            if ref:
                loc_str += f", Ref: {ref}"
            if alt:
                loc_str += f", Alt: {alt}"
            self.variant_summaries.insert(0, loc_str)
        
        # insert rsID at the very front
        if self.rsID:
            self.variant_summaries.insert(0, f"Begin info for Variant ID: {self.rsID}")
        self.variant_summaries.insert(0, "*")
        print(self.variant_summaries)
        return ' | '.join(self.variant_summaries)+'*.'
    
    def add_functional_prediction(self, gene: str, tool: str, score: Any):
        """Add a functional prediction score for a specific gene from a specific tool.
        If already exists, add new one to it
        """
        if gene not in self.functional_predictions:
            self.functional_predictions[gene] = {}
        if tool not in self.functional_predictions[gene]:
            self.functional_predictions[gene][tool] = []
        if isinstance(score, list):
            self.functional_predictions[gene][tool].extend(score)
        else:
            self.functional_predictions[gene][tool].append(score)
        # deduplicate
        self.functional_predictions[gene][tool] = list(set(self.functional_predictions[gene][tool]))
    
    def add_location(self, build: str, chrom: str, pos: int, ref: Optional[List[str]] = None, alt: Optional[List[str]] = None):
        """Add genomic location for a specific genome build, and reference/alternate alleles if available.
        note that there could be multiple ref/alt alleles for a given position; if so we update current list it is exists
        """
        if build not in self.loc_by_build:
            self.loc_by_build[build] = {}
        self.loc_by_build[build]['chrom'] = chrom
        self.loc_by_build[build]['pos'] = pos
        if ref:
            if 'ref' not in self.loc_by_build[build]:
                self.loc_by_build[build]['ref'] = []
            self.loc_by_build[build]['ref'] = list(set(self.loc_by_build[build]['ref']) | set(ref))
        if alt:
            if 'alt' not in self.loc_by_build[build]:
                self.loc_by_build[build]['alt'] = []
            self.loc_by_build[build]['alt'] = list(set(self.loc_by_build[build]['alt']) | set(alt))
            
    def get_location(self, build: str) -> Optional[Dict[str, Any]]:
        """Retrieve genomic location for a specific genome build."""
        # builds are tricky, will allow for partial matches if exact not found
        if build in self.loc_by_build:
            return self.loc_by_build[build]
        for b in self.loc_by_build.keys():
            if build.lower() in b.lower() or b.lower() in build.lower():
                return self.loc_by_build[b]
        return None
    
    def add_drug_response_effect(self, variant_effect: str) -> None:
        variant_effect = (variant_effect or "").strip()
        if variant_effect and variant_effect not in self.drug_response_effects:
            self.drug_response_effects.append(variant_effect)

    def add_many_drug_response_effects(self, variant_effects: List[str]) -> None:
        for variant_effect in variant_effects:
            self.add_drug_response_effect(variant_effect)
    
    def get_related_genes(self) -> List[str]:
        """Return a list of genes related to this variant."""   
        return self.genes_related_to
        
    def add_af_freq(self, info: List[Dict[Any, Any]]):
        """Add allele frequency information from a specific source and population."""
        self.af_freqs.extend(info)
    
    
    def add_per_gene_traits(self, gene: str, traits: Dict[str, Any]):
        """Add traits associated with a specific gene.
        Expect traits like "{trait: {'p_value': float, 'risk_score': float}"
        """
        if gene not in self.per_gene_traits.keys():
            self.per_gene_traits[gene] = {}
        for trait in traits.keys():
            if trait not in self.per_gene_traits[gene].keys():
                self.per_gene_traits[gene][trait] = traits[trait]
        self.traits = list(set(self.traits) | set(traits.keys()))
        if gene not in self.genes_related_to:
            self.genes_related_to.append(gene)
        
        
    def add_per_gene_context(self, gene: str, context: str, variant_category: str):
        """Add context information associated with a specific gene.
        Expect context like  {'context': 'intron_variant', 'variant_category': 'intronic'}
        """
        if gene not in self.per_gene_context:
            self.per_gene_context[gene] = {}
        self.per_gene_context[gene]['context'] = context
        self.per_gene_context[gene]['variant_category'] = variant_category
        if gene not in self.genes_related_to:
            self.genes_related_to.append(context['gene_name'])
            
    def add_tool(self, tool_name: str) -> None:
        """Record that a tool was run for this gene."""
        name = tool_name.strip()
        if name and name not in self.tools_run:
            self.tools_run.append(name)
        
    

    