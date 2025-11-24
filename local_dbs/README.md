# Local Database Files

## Required Downloads

The following files must be downloaded before using the GWAS and Sei tools:

### Required Files:
- `gwas_catalogue_association.tsv`
**Download Location:** https://drive.google.com/drive/folders/1a0zqcyt84Iy2D0DcCr-1VGFaTM96opMS?usp=share_link
- `gwas_catalogue_studies.tsv`
- `seqclass.names`
**Download Location:** https://github.com/FunctionLab/sei-framework/blob/main/model/seqclass.names
- `sorted.hg38.tiling.bed.ipca_randomized_300.labels.merged.bed`
**Download Location:** https://zenodo.org/records/7113989
- `sorted.hg19.tiling.bed.ipca_randomized_300.labels.merged.bed`
**Download Location:** https://zenodo.org/records/7113989
- `AlphaMissense_hg38.tsv`
**Download Location:** https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
- `AlphaMissense_hg19.tsv`
**Download Location:** https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg19.tsv.gz
- `NCBI2Reactome_All_Levels.txt`
**Download Location:** https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt
- `old_go_data: Version 2025-07-22`
**Download Location:** https://current.geneontology.org/products/pages/downloads.html

**Donwload location** https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.primary_assembly.basic.annotation.gtf.gz
-  Gencode annotaitons

- `alvessa_proteins.db`   
**Download Location:** https://drive.google.com/file/d/1p9gg-iogIIu7LAmwPBmoSVQtJbalDrTj/view?usp=drive_link  
Database supporting protein structure visualization and druggability analysis (AlphaFold pLDDT, FPocket, SASA, IUPPred3, and BioLiP2 binding sites).

- `DisProt_release_2025_06_with_ambiguous_evidences.json`  
**Download Location:** https://drive.google.com/file/d/1sLM0mw3bu8rqOAX47Xurx32UujeLOVu4/view?usp=drive_link  
DisProt v2025-06 dataset with ambiguous evidence annotations for protein disorder, used in disorder consensus overlays, used by `tools.prot`

- `chembl_35.db`  
**Download Location:** https://drive.google.com/file/d/1P5WMBrbpmW6aOIY9hP2a74-zmIEwM7eC/view?usp=drive_link  
Local ChEMBL v35 database for querying drug–target, drug approval, trial, and assay bioactivity data.

- `pdb.zip` (unzip into the `pdb/` folder)  
**Download Location:** https://drive.google.com/file/d/17uZse6xB8T_U_XcUeVWdCntqfmHIsfv0/view?usp=drive_link  
Contains PDB protein structures (from AlphaFold DB) corresponding to proteins in `alvessa_proteins.db`, used by `tool_prot` for structural visualization and druggability.

- `uniprot_hs_9606_reviewed.db`
**Download Location:** https://drive.google.com/file/d/1PHNgkvjzYa0sfZ3pFi-fhA74xwVKp93j/view?usp=drive_link
Locally built from UniProtKB reviewed human proteome (Homo sapiens, tax_id 9606).
SQLite database of reviewed human UniProt entries for sequence-to-gene resolution and protein annotation (accessions, gene names, isoforms, and core metadata), including a custom index for quick searching genes by amino acid sequences, used by `aa_seq` tool.

- `drugcentral.dump.11012023.db`
**Download Location:** https://drive.google.com/file/d/1FrovOP2Uab4pyLsgc0XKc5qwNnTYEnQ9/view?usp=drive_link
Source: Locally built from DrugCentral PostgreSQL dump (2023-11-01).
SQLite snapshot of DrugCentral for querying approved and investigational drugs, indications, mechanisms of action, targets, and key safety/regulatory annotations, used by `drug_cntrl` tool.