"""
Description:
Gencode GTF query tool: annotate variants with GTF features or get locations of genes
"""
# %%
import pandas as pd
from typing import Dict, Any, Optional, List
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))


# Load expanded GTF
try:
    gtf_df = pd.read_parquet("local_dbs/gencode.v48.expanded.parquet")
except Exception as e:
    print(f"[Gencode] Failed to load Gencode data: {e}")
    gtf_df = pd.DataFrame()

FEATURES_TO_FLAG = ["gene", "transcript", "exon", "CDS", "UTR", "start_codon", "stop_codon"]

def annotate_variant(chrom: str, pos: int) -> Optional[Dict[str, Any]]:
    """
    Annotate a genomic variant by checking overlaps with GENCODE features.

    Args:
        chrom: Chromosome name (e.g., "1" or "chr1").
        pos: 1-based genomic position.

    Returns:
        Dictionary with:
            - chrom, pos
            - feature_flags: Boolean flags for each feature type
            - feature_gene_lists: List of genes for each feature type
            - matches: DataFrame with matching features
            - genes_related_to: List of unique gene names related to the variant
   
    """
    
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    chrom_data = gtf_df[gtf_df["seqname"] == chrom]
    matches = chrom_data[(chrom_data["start"] <= pos) & (chrom_data["end"] >= pos)]
    if matches.empty:
        return None

    matches = matches[matches["feature"] != "Selenocysteine"]

    # Create flags for features
    feature_flags = {feat: False for feat in FEATURES_TO_FLAG}
    for feat in FEATURES_TO_FLAG:
        if feat in matches["feature"].values:
            feature_flags[feat] = True
            
    genes_related_to = matches["gene_name"].unique()
    
    # Feature-specific gene lists
    feature_gene_lists: Dict[str, List[str]] = {}
    for feat in FEATURES_TO_FLAG:
        genes_for_feat = sorted(
            matches.loc[matches["feature"] == feat, "gene_name"].dropna().unique()
        )
        feature_gene_lists[f"{feat}_genes"] = genes_for_feat

    return {
        "chrom": chrom,
        "pos": pos,
        'feature_flags': feature_flags,
        'feature_gene_lists': feature_gene_lists,
        "matches": matches[["feature", "gene_name", "transcript_id", "start", "end", "strand"]]
                    .drop_duplicates()
                    .to_dict(orient="records"),
        "genes_related_to": genes_related_to.tolist()
    }


def summarize_gene_structure(gene_name: str) -> Optional[Dict[str, Any]]:
    """
    Summarize transcript and exon structure for a given gene.

    Args:
        gene_name: Gene symbol (e.g., "TP53")

    Returns:
        Dict with:
            - n_transcripts
            - exons_per_transcript
            - min_transcript_span_bp
            - max_transcript_span_bp
            - transcript_ids
            - gene_span_bp (tuple: start, end)
        Returns None if gene not found.

    """
    gene_df = gtf_df[gtf_df["gene_name"] == gene_name]
    if gene_df.empty:
        return None

    transcripts = gene_df["transcript_id"].dropna().unique().tolist()
    exon_counts = {}
    spans = []

    for tx in transcripts:
        tx_rows = gtf_df[gtf_df["transcript_id"] == tx]
        exon_count = (tx_rows["feature"] == "exon").sum()
        exon_counts[tx] = int(exon_count)
        spans.append(int(tx_rows["end"].max() - tx_rows["start"].min() + 1))

    gene_span = (int(gene_df["start"].dropna().min()), int(gene_df["end"].dropna().max()))
    gene_type = gene_df["gene_type"].unique()
    if len(gene_type) != 1:
        gene_type = ",".join(gene_type)  # Multiple types, join as string
    else:
        gene_type = gene_type[0]
        
    gene_id = gene_df["gene_id"].unique()
    if len(gene_id) != 1:
        gene_id = ",".join(gene_id)  # Multiple IDs, join as string
    else:
        gene_id = gene_id[0]
        
    chromosome = gene_df["seqname"].unique()
    if len(chromosome) != 1:
        chromosome = ",".join(chromosome)  # Multiple chromosomes, join as string
    else:
        chromosome = chromosome[0]
        
    strand = gene_df["strand"].dropna().unique()
    if len(strand) != 1:
        strand = None
    else:
        strand = strand[0]
    
    return {
        "gene_name": gene_name,
        "gene_type": gene_type,
        'gene_id': gene_id,
        "chromosome": chromosome,
        "exons_per_transcript": exon_counts,
        "transcript_ids": transcripts,
        "gene_span_bp": gene_span, 
        "strand": strand,
    }

# %%
