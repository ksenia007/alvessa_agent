from langchain.agents import Tool
from entity_extraction import entity_extraction_node
from tool_humanbase import humanbase_predictions_agent, expectosc_predictions_agent
from tool_biogrid import bioGRID_predictions_agent
from tool_reactome import reactome_pathways_agent
from tool_go_summarization import go_summarization_agent
from tool_uniprot import (
    uniprot_node,
)
from tool_gwas import gwas_associations_agent, query_by_trait_agent
from conditioned_claude import conditioned_claude_node
from verify import verify_evidence_node
from tool_dbsnp import dbsnp_variants_agent
from tool_sei import sei_predictions_agent
from tool_alphamissense import alphamissense_predictions_agent
from tool_annotate_gencode import gencode_gene_node
from tool_miRDB import miRDB_agent
from tool_remap import remap_crm_agent
from tool_prot import prot_agent

EXAMPLE_TOOL_SELECTION = """EXAMPLE PIPELINES (pay attention to dependencies):

1. Variant regulatory activity (e.g. SEI):
   ["extract_entities", "query_gwas_by_gene", "variant_annotations", "sei"]

2. Variant pathogenicity (e.g. AlphaMissense):
   ["extract_entities", "query_gwas_by_gene", "variant_annotations", "alphamissense"] 

3. HumanBase Expecto variant annotation:
   ["extract_entities", "humanbase_expecto", "query_gwas_by_gene", "variant_annotations", "humanbase_tissue_expecto_annotate_variants"]

4. Gene-level functional annotation:
   ["extract_entities", "gencode_gene_node", "humanbase_functions", "uniprot_base", "reactome", "BioGRID", "Summarize_bioGRID_GO", "uniprot_gwas"]

5. Protein structure, visualization and druggability:
   ["prot"]


Note these are only examples, and in real life you may need to run combinations of these tools **depending on the user intent and the entities extracted**.

Respond *ONLY* with a Python list of tool names. Example: ["humanbase_functions", "uniprot_base"] or ["humanbase_functions", "uniprot_base", "query_gwas_by_gene"] or ["query_gwas_by_gene", "BioGRID"]

IF YOU INCLUDE ANYTHING OTHER THAN A LIST OF TOOLS, YOU WILL BE CONSIDERED INVALID.
"""

TOOL_CATALOG = {
    "extract_entities": "Extracts genes and all biomedical entities from the user question for query understanding using both Claude and GLiNER models. Returns genes (Claude + GLiNER), traits, and all entity types found by GLiNER. This tool MUST be called first before running any other tool as it provides essential entity extraction for pipeline execution.",
    # TODO: need to add support w/ classes "expand_gene_set_by_trait": "Use this tool ONLY IF the expanded gene set by a trait is REQUIRED to answer the question. This could be used to discover more relevant genes underlying a trait, if a broader genetic context would be valuable for the analysis. Searches for genetic associations with diseases, traits, or phenotypes mentioned in the user query and expands the gene list with additional related genes found through GWAS associations. ",
    "gencode_gene_node": "Annotates genes with GENCODE gene annotations such as transcripts, number exons, genomic span. Essential for many downstream analyses.",
    "humanbase_functions": "Fetch per-gene functional predictions from HumanBase tissue-specific networks. Provides expanded list of functions.",
    "uniprot_base":  "Queries UniProt for functional annotations and disease links.",
    "query_gwas_by_gene":  "Retrieves genome-wide association study (GWAS) results for a given gene. It collects traits and diseases associated with genetic variants linked to that gene, along with the specific variants",
    "query_gwas_extensive": "This is a more comprehensive version of the query_gwas_by_gene tool, and it is used to retrieve more detailed information about the GWAS results. It collects an extensive list of traits/diseases associated with an extensive list of genetic variants linked to that gene. Use this tool *ONLY* if the question is very specific that requires what is equivalent to an extensive database search, not to general characterisation of the gene.",
    "uniprot_gwas":   "Runs UniProt query again on genes identified via GWAS associations. Helps to expand the base annotations with related genes.",
    "BioGRID": "Fetches gene interactions from BioGRID and their functional annotations for the input genes. Provides a curated context-specific list of protein-protein, genetic and chemical interactions.",
    "reactome": "Fetches Reactome pathways associated with the input genes. Provides a curated collection of biological pathways which describe how molecules interact within a cell to carry out different biological processes.",
   #  "Summarize_bioGRID_GO": "Required for BioGrid. Summarizes BioGRID GO terms for the input genes. Provides a compact list of GO terms for the input genes.",
   "Summarize_GO": "Summarizes GO terms for the input genes. Can only ber run after UniProt, needed to answer questions about GO terms or functions of genes.",

    "variant_annotations": "Fetches dbSNP data about the identified variants. This requires gwas to be run first. This is to be used only when the variant needs to be annotated with its **coordinates** and **chromosome number** and similar details. This tool is also a prerequisite to run humanbase, sei, alphamissense, etc. Reason based on the user's question if sequence models need to be run, and if so, this tool is a prerequisite to run them.",
    "variant_population_summaries": "Fetches population-wide summaries of allele frequencies for the identified variants across studies. Useful for characterizing the frequency of variants in the population. This requires gwas to be run first.",
    "sei": "Fetches predictions of the sequence regulatory activity for given variants. This requires variant_annotations to be run first.", 
    "expectosc_predictions_agent": "Annotates variants with predicted (from sequence) cell type-specific expression disruption predictions. Requires variant_annotations to be run first.",
    "alphamissense": "Fetches Alphamissense predicted pathogenicity classes for given variants. This requires variant_annotations to be run first.",
    "miRDB": "Fetches miRDB computationally predicted gene targets of miRNA.",
    "remap_crm_agent": "Fetches cis-regulatory modules (CRMs) from the ReMap 2022 database for each gene, reporting transcription factors with ChIP-seq binding peaks near gene’s TSS. Requires Gencode (gencode_gene_node) run first. Useful for exploring transcriptional regulation.",
    "prot": "Retrieves structural data for a single or several proteins given one or several gene symbols. Resolves UniProt ID and AlphaFold Protein Database structure, and provides per-residue metrics: pLDDT confidence scores (structural reliability) and FPocket druggability scores. Outputs both an interactive 3Dmol.js visualization and a text summary of min, max and averaged metrics."
}


TOOL_FN_MAP = {
    "extract_entities":   entity_extraction_node,
    "expand_gene_set_by_trait":     query_by_trait_agent,
    "humanbase_functions":      humanbase_predictions_agent,
    "expectosc_predictions_agent": expectosc_predictions_agent,
    "uniprot_base":   uniprot_node,
    "query_gwas_by_gene":           gwas_associations_agent,
    "query_gwas_extensive":          lambda x: gwas_associations_agent(x, mode="extensive"),
    "uniprot_gwas":   uniprot_node,       
    "BioGRID":        bioGRID_predictions_agent,
    "reactome":       reactome_pathways_agent,
    "Summarize_GO": go_summarization_agent, 
    "variant_annotations":          dbsnp_variants_agent,
    "variant_population_summaries": lambda x: dbsnp_variants_agent(x, include_population_summaries=True),
    "sei":            sei_predictions_agent,
    "alphamissense":  alphamissense_predictions_agent,
    "gencode_gene_node": gencode_gene_node,
    "miRDB": miRDB_agent,
    "remap_crm_agent": remap_crm_agent,
    "prot": prot_agent
}

