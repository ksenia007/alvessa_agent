"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-11-24
Updated: 


Description: 

ClinVar tool

"""
from __future__ import annotations

import time

from src.state import State
from src.tools.base import Node
from pathlib import Path
from typing import Dict, Iterable, List, Sequence
from src.alvessa.domain.variant_class import Variant

import json
import pandas as pd

DEBUG = True

CLINVAR_GENE_DISEASES = Path(__file__).resolve().parents[3] / "local_dbs" / "clinvar" / "clinvar_gene_condition_source_id.parquet"
CLINVAR_DF = Path(__file__).resolve().parents[3] / "local_dbs" / "clinvar" / "clinvar_pathogenic_variants_grch38.parquet"

def clinvar_node(state: "State") -> "State":
    """
    LangGraph node that annotates genes and variants with ClinVar data.

    Parameters
    ----------
    state
        Current graph state.

    Returns
    -------
    State
        Updated state with the omim fields filled.
        
    """
    gene_entities = state.get("gene_entities") or {}
    variant_entities = state.get("variant_entities") or {}
    
    try: 
        clinvar_df = pd.read_parquet(CLINVAR_DF)
        clinvar_gene_diseases_df = pd.read_parquet(CLINVAR_GENE_DISEASES)
    except Exception as e:
        if DEBUG:
            print(f"[ClinVar] Error loading ClinVar data: {e}")
        return 
    
    genes_with_variants = set(clinvar_df.GeneSymbol.dropna().unique())
    genes_with_disease = set(clinvar_gene_diseases_df.AssociatedGenes.dropna().unique())

    # add ClinVar data to gene entities
    for gene_name, gene in gene_entities.items():
        if not gene_name:
            continue
        
        try:
            if gene_name in genes_with_disease:
                subset = clinvar_gene_diseases_df[
                    clinvar_gene_diseases_df['AssociatedGenes'] == gene_name
                ]
                # there is DiseaseName+SourceName+LastUpdated (the last is date)
                diseases = subset[['DiseaseName', 'SourceName', "LastUpdated"]].drop_duplicates()
                disease_lines = []
                for _, row in diseases.iterrows():
                    disease_lines.append(f"{row['DiseaseName']} (Source: {row['SourceName']}, Last Updated: {row['LastUpdated']})")
                gene.update_text_summaries(
                f"ClinVar-sourced diseases for {gene_name}: " + ", ".join(disease_lines) + f"|"
                )
        except:
            if DEBUG:
                print(f"[ClinVar] Error processing diseases for gene {gene_name}")
        try:         
            if gene_name in genes_with_variants:
                subset = clinvar_df[
                    clinvar_df['GeneSymbol'] == gene_name
                ]
                # for each line we want ALL columns that are not "-" and not NaN and not "na", and keep column name as explanaiton
                variant_lines = []
                for _, row in subset.iterrows():
                    variant_info_short = []
                    variant_info_short.append(f"rs{row['rsid']}:")
                    for col in ['Type', 'Name', 'ClinicalSignificance']:
                        val = row[col]
                        if pd.isna(val) or val in ["-", "na", ""]:
                            continue
                        variant_info_short.append(f"{val},")
                                            
                    # is rsID is not none and "rs"+str(rsID) is not in variant's existing identifiers, add it as new object
                    if 'rsid' in row and pd.notna(row['rsid']):
                        rsid = str(row['rsid'])
                        if rsid.startswith("rs"):
                            rsid_clean = rsid
                        else:
                            rsid_clean = "rs" + rsid
                        if rsid_clean not in variant_entities.keys():
                            # create new variant object
                            new_variant = Variant(
                                rsID = rsid_clean,
                                tools_run=['clinvar_node'],
                            )
                            # add_location for GRCh38 if available: Chromosome, Start, ReferenceAllele, AlternateAllele
                            if 'Chromosome' in row and pd.notna(row['Chromosome']) and \
                                'Start' in row and pd.notna(row['Start']):
                                chrom = str(row['Chromosome'])
                                pos = int(row['Start'])
                                ref = str(row['ReferenceAllele']) if 'ReferenceAllele' in row and pd.notna(row['ReferenceAllele']) else None
                                alt = str(row['AlternateAllele']) if 'AlternateAllele' in row and pd.notna(row['AlternateAllele']) else None
                                new_variant.add_location(build="GrCh38", chrom=chrom, pos=pos, ref=[ref] if ref else None, alt=[alt] if alt else None)
                            variant_info_long = []
                            for col in ['Name', 'ClinicalSignificance', 'SomaticClinicalImpact', 'Oncogenicity', 'OriginSimple', 'PhenotypeList']:
                                val = row[col]
                                if pd.isna(val) or val in ["-", "na", ""]:
                                    continue
                                variant_info_long.append(f"{val}")
                            # add variant_info as text summary with update_text_summaries
                            new_variant.update_text_summaries("ClinVar variant info: " + ", ".join(variant_info_long))
                            variant_entities[rsid_clean] = new_variant
                    variant_lines.append(", ".join(variant_info_short))
                gene.update_text_summaries(f"ClinVar pathogenic variants for {gene_name}: " + "; ".join(variant_lines)+ f"|")
        except:
            if DEBUG:
                print(f"[ClinVar] Error processing variants for gene {gene_name}")
    if DEBUG:
        print(f"[ClinVar] ClinVar fetched")

    return 

NODES: tuple[Node, ...] = (
    Node(
        name="clinvar_node",
        entry_point=clinvar_node,
        description="Add information for genes and pathogenic variants from the ClinVar database. For genes it pulls in associated diseases as well as variants. For variants it adds clinical significance and associated conditions. Should be run after GWAS (if used) and before other variant annotation tools.",
    ),
)