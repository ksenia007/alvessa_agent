from langchain.agents import Tool
from entity_extraction import gene_extraction_node
from tool_humanbase import humanbase_predictions_agent
from tool_biogrid import bioGRID_predictions_agent
from tool_uniprot import (
    uniprot_node,
    trait_disease_extraction_node,
    trait_function_extraction_node,
    trait_GO_extraction_node,
)
from tool_gwas import gwas_associations_agent
from conditioned_claude import conditioned_claude_node
from verify import verify_evidence_node

TOOL_CATALOG = {
    "humanbase": "Fetch functional predictions from HumanBase tissue-specific networks. Provides expanded list of functions.",
    "uniprot_base":  "Queries UniProt for annotations, disease links, and GO terms from the input genes. Provides base annotations.",
    "trait_disease":  "Extracts disease-related traits from UniProt annotations. Provides well known disease associations.",
    "trait_function": "Extracts functional traits from UniProt annotations. Provides base (non-enriched) functional annotations.",
    "trait_go":       "Extracts Gene Ontology terms from UniProt annotations. Provides a list of Gene Ontology terms associated with the genes.",
    "gwas":           "Fetches GWAS associations for genes, including diseases, related variants and genes. Population-level variant associations from common variants.",
    "uniprot_gwas":   "Runs UniProt query again on genes identified via GWAS associations. Helps to expand the base annotations with related genes.",
    "BioGRID": "Fetches BioGRID interactions and their functional annotations for the input genes. Provides a curated context-specific list of protein-protein, genetic and chemical interactions.",
}


TOOL_FN_MAP = {
    "humanbase":      humanbase_predictions_agent,
    "uniprot_base":   uniprot_node,
    "trait_disease":  trait_disease_extraction_node,
    "trait_function": trait_function_extraction_node,
    "trait_go":       trait_GO_extraction_node,
    "gwas":           gwas_associations_agent,
    "uniprot_gwas":   uniprot_node,       
    "BioGRID":        bioGRID_predictions_agent,
}

