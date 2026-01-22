# Amino Acid Sequence–to–Gene Resolution Tool Tutorial

This tutorial explains how to use and interpret the **Amino Acid (AA) Sequence–to–Gene Resolution Tool** in **Alvessa**, using a concrete example to illustrate how short peptide sequences are resolved to genes, UniProt accessions, and isoforms.

The tool maps short peptide fragments or full-length amino acid sequences to **gene-associated UniProt records** using a **local SQLite snapshot of UniProtKB**, without requiring internet access.

---
## When the AA Sequence Tool Is Used 

The **AA Sequence Tool is invoked automatically** when Alvessa detects that the user input contains **amino acid sequences** rather than gene symbols. Users do not need to explicitly select the tool. 

## Underlying Database 
The resolver uses a **local SQLite snapshot of UniProtKB** that includes: 
- Canonical and non-canonical UniProt accessions
- Full protein amino acid sequences 
- Primary gene symbols
- Optional Entrez Gene identifiers 
 
 The database is constructed offline and reused across all Alvessa sessions, ensuring reproducibility and eliminating runtime dependency on external services.

### Example queries
 - “Identify the gene for this peptide: MEEPQSDPSV" 
 - “Which gene does this amino acid sequence correspond to? FGEVAKQEEFFNLSHCQLVTLISRDDLNVR” 
 - “Identify the gene for this peptide: CAQYWPQKEEKEMIFEDTNL” 
 - “Map this protein sequence to UniProt and gene information PQYRLEKQSEPNVAVDLDSTLESQSAWEFC"
 - “Resolve the following amino acid sequence to a human gene AAGGYDGQDQLNSVERYDVETETWTFVAPMKHRRSALGIT” 
 
 If multiple sequences are provided in a single query, each sequence is processed independently.

## Example question with output

> **Identify the gene for this peptide:**  
> `MEEPQSDPSV`

This is a short (10–amino acid) peptide corresponding to the extreme N-terminus of a well-known human protein. Because short peptides can be shared across multiple isoforms, this example illustrates how the tool reports **gene-level identity**, **canonical UniProt accessions**, and **isoform-level matches** together.

---

## Question and Text Summary Output

When the tool is triggered, Alvessa first displays a **text-based answer summary** that interprets the sequence resolution results at a high level.

<img src="img/aa_seq_text_summary_overview.png" alt="Question and text summary overview" width="80%">

*Question and text summary overview*

For this query, the summary reports:

- The peptide `MEEPQSDPSV` maps to the **TP53** gene  
- **Gene:** TP53  
- **Entrez Gene ID:** 7157  
- **Canonical UniProt accession:** P04637  

The summary further explains that:

- The query sequence matches **UniProt P04637** with **100% identity**
- The match corresponds to the **N-terminal sequence** of TP53
- TP53 has a total length of **393 amino acids**
- The same N-terminal sequence is present in multiple TP53 isoforms

Specifically, the following isoforms contain an identical match to this peptide:

- **Isoform 1:** P04637-1 (p53 / p53α)  
- **Isoform 2:** P04637-2 (p53β)  
- **Isoform 3:** P04637-3 (p53γ)  

All three isoforms show a **100% identity match** to the 10–amino acid query sequence at their N-terminus.

This summary is intended to answer the user’s question directly while flagging important biological context—namely, that **isoform-level ambiguity is expected for short peptides**.

---

## AA Sequence Mapping (UniProt / Genes) Panel

To inspect the matches in detail, the **AA Sequence Mapping (UniProt / Genes) panel** provides a structured, tabular view of all resolved UniProt hits.

<img src="img/aa_sequence_mapping_panel_overview.png" alt="AA Sequence Mapping panel overview" width="80%">

*AA Sequence Mapping panel overview*

This panel separates results into:
- **Gene-associated UniProt hits**
- **UniProt accessions without gene annotation**

---

## Mapped UniProt Hits for TP53

When **TP53 (Entrez Gene ID: 7157)** is selected, the panel displays all UniProt accessions associated with this gene that match the query sequence.

<img src="img/aa_sequence_gene_hits.png" alt="Mapped UniProt hits for TP53" width="80%">

*Mapped UniProt hits for TP53*

In this example, the following entries are shown:

- **TP53 — P04637**  
  - Canonical accession  

- **TP53 — P04637-2**  
  - Canonical accession: P04637  

- **TP53 — P04637-3**  
  - Canonical accession: P04637  

Each row reports:
- Gene symbol (TP53)
- UniProt accession (canonical and isoform-specific)
- Entrez Gene ID (7157)
- Coverage
- Alignment length
- Query length
- Query sequence
- UniProt reference sequence

These results indicate that the peptide sequence is **not isoform-specific**: it occurs identically in multiple TP53 isoforms.  
Accordingly, the tool resolves the query confidently at the **gene level**, while retaining **isoform-level transparency**.

---

## UniProt Accessions Without Gene Annotation

In addition to gene-associated hits, the panel may also display **UniProt accessions without gene annotation**.

<img src="img/aa_sequence_no_gene_hits.png" alt="UniProt accessions without gene annotation" width="80%">

*UniProt accessions without gene annotation*

In this example, one such entry is shown:

- **Accession:** Q9TTA1  
- **Canonical accession:** Q9TTA1  

This table includes:
- Accession and canonical accession
- Similarity score
- Coverage
- Alignment length
- Query length
- Query sequence
- UniProt reference sequence

Entries in this section typically correspond to:
- UniProt records lacking an assigned gene symbol
- Historical or poorly annotated entries
- Non-standard or unreviewed sequences

These hits are preserved for **completeness and auditability**, but they are **not promoted to gene-level objects** unless additional annotation becomes available.

---

## Interpretation for This Example

For the peptide `MEEPQSDPSV`:

- The **gene identity is unambiguous**: **TP53**
- The peptide matches the **canonical TP53 protein** (P04637) with 100% identity
- Multiple TP53 isoforms share the same N-terminal sequence
- Isoform-level resolution is **not possible or necessary** for this peptide alone
- The tool therefore:
  - Resolves the query to **TP53 at the gene level**
  - Retains isoform matches for transparency
  - Flags unannotated UniProt accessions separately

This behavior is expected and correct for short peptides derived from conserved regions.

---

## Key Takeaways

- Short peptides often map to **multiple isoforms of the same gene**
- The AA Sequence Tool prioritizes **gene-level resolution**
- Canonical UniProt accessions are preferred when available
- Isoforms are retained to reflect biological reality
- Unannotated UniProt accessions are shown but not promoted

---

## Summary

Using the example peptide `MEEPQSDPSV`, the Amino Acid Sequence–to–Gene Resolution Tool demonstrates how Alvessa:

- Resolves short amino acid sequences to genes
- Handles canonical UniProt accessions and isoforms transparently
- Distinguishes gene-associated hits from unannotated UniProt records
- Provides clear, interpretable results suitable for downstream analysis

**Note:** This tool is intended for research and hypothesis generation only and does not provide clinical or diagnostic interpretations.
