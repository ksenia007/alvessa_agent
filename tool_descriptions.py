from langchain.agents import Tool
from entity_extraction import gene_extraction_node
from tool_humanbase import humanbase_predictions_agent
from tool_biogrid import bioGRID_predictions_agent
from tool_go_summarization import make_go_summarization_node
from tool_uniprot import (
    uniprot_node,
)
from tool_gwas import gwas_associations_agent
from conditioned_claude import conditioned_claude_node
from verify import verify_evidence_node
from tool_dbsnp import dbsnp_variants_agent
from tool_sei import sei_predictions_agent

TOOL_CATALOG = {
    "humanbase": "Fetch functional predictions from HumanBase tissue-specific networks. Provides expanded list of functions.",
    "uniprot_base":  "Queries UniProt for functional annotations and disease links",
    "gwas":           "Fetches GWAS associations for genes, including diseases, related variants and genes. Population-level variant associations from common variants.",
    "uniprot_gwas":   "Runs UniProt query again on genes identified via GWAS associations. Helps to expand the base annotations with related genes.",
    "BioGRID": "Fetches BioGRID interactions and their functional annotations for the input genes. Provides a curated context-specific list of protein-protein, genetic and chemical interactions.",
    "Summarize_bioGRID_GO": "Required for BioGrid. Summarizes BioGRID GO terms for the input genes. Provides a compact list of GO terms for the input genes.",
    "dbsnp": "Fetches dbSNP data about the identified variants. This requires gwas to be run first. Returns genomic coordinates and allele/variant frequencies across different studies/populations.",
    "sei": "Fetches predictions of the sequence regulatory activity for given variants. This requires dbsnp to be run first."
}


TOOL_FN_MAP = {
    "humanbase":      humanbase_predictions_agent,
    "uniprot_base":   uniprot_node,
    "gwas":           gwas_associations_agent,
    "uniprot_gwas":   uniprot_node,       
    "BioGRID":        bioGRID_predictions_agent,
    "Summarize_bioGRID_GO": make_go_summarization_node, 
    "dbsnp":          dbsnp_variants_agent,
    "sei":          sei_predictions_agent,
}

