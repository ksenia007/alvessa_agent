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
**Download Location:** https://drive.google.com/file/d/1YmZTGjNjmgWCbM2yKr_Ajp3YYeaHeFiE/view?usp=drive_link
Contains pLDDT and FPocket data for proteins, used by `tool_prot` to visualize protein structures and druggability.

- `pdb.zip` (unzip into the `pdb/` folder)  
**Download Location:** https://drive.google.com/file/d/17uZse6xB8T_U_XcUeVWdCntqfmHIsfv0/view?usp=drive_link  
Contains PDB protein structures (from AlphaFold DB) corresponding to proteins in `alvessa_proteins.db`, used by `tool_prot` for structural visualization and druggability.
