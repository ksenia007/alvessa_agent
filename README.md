<p align="center">
  <img src="https://github.com/user-attachments/assets/2bdd5321-8ffe-44ae-aaff-be11a715bdc5" alt="Alvessa" width="1100"/>
</p>
<h1 align="center">
  <strong>An Evidence-Grounded Research Assistant for Functional Genomics and Drug Target Assessment</strong>
</h1>

<p align="center">
  <a href="https://alvessa.ai"><img src="https://img.shields.io/badge/🌐_Live_Demo-alvessa.ai-blue?style=for-the-badge" alt="Live Demo"/></a>
  <a href="https://github.com/ksenia007/alvessa_agent/tree/main/tutorials"><img src="https://img.shields.io/badge/📚_Tutorials-GitHub-yellow?style=for-the-badge" alt="Tutorials"/></a>
  <a href="https://www.biorxiv.org/content/10.64898/2025.12.30.697073v1"><img src="https://img.shields.io/badge/📄_Manuscript-BioRxiv-red?style=for-the-badge" alt="Manuscript"/></a>
</p>

<p align="center">
  <a href="#-key-features">Features</a> •
  <a href="#-quick-start">Quick Start</a> •
  <a href="#-usage">Usage</a> •
  <a href="#-contribute-tools">Contribute Tools</a> •
  <a href="#-genomearena">Benchmarks</a> •
  <a href="#-contribute-questions">Contribute Questions</a> •
  <a href="#-citation">Citation</a>
</p>

---

Alvessa is a multi-agent framework that provides **verifiable, evidence-grounded answers** to genomics and drug target assessment questions. Unlike general-purpose LLMs that can hallucinate identifiers or fabricate specifics, Alvessa enforces statement-level verification against retrieved database records, with feedback and explanations for each claim.


## ✨ Key Features

| Capability | Description |
|------------|-------------|
| 🔍 **Entity Recognition** | Ensemble NER for genes, variants, drugs, miRNAs, and protein sequences |
| 🛠️ **Tool Orchestration** | Context-aware orchestration of validated tools (dbSNP, ClinVar, UniProt, ChEMBL, AlphaFold, etc.) |
| 📎 **Evidence Grounding** | Answers constrained to retrieved database records with inline citations |
| ✅ **Verification Loop** | Per-statement evaluation flags unsupported claims and triggers revision |
| 🖥️ **Interactive UI** | Web interface with citation links and verification feedback |

## 🚀 Quick Start

```bash
# Clone and install
git clone https://github.com/ksenia007/alvessa_agent.git
cd alvessa_agent
conda create -n agents python=3.10 -y
conda activate agents
pip install -e .

# download and unzip local_dbs (45GB)
cd local_dbs
curl -O https://alvessa-public-access-bucket.s3.us-east-1.amazonaws.com/local_dbs.zip
unzip local_dbs.zip
cd ..

# Set API keys
export ANTHROPIC_API_KEY="your-key"
export BIOGRID_API_KEY="your-key"        # optional
export DISGENET_API_KEY="your-key"       # optional
```
Some tools rely on databases that require separate registration and local downloads due to licensing restrictions. These resources are optional and only needed if you plan to use the corresponding tools.

**OMIM.**
Access to OMIM downloadable data requires registration at [https://www.omim.org/downloads](https://www.omim.org/downloads). Once access is granted, download the file `genemap2.txt` and place it directly into the `local_dbs/` folder. No additional preprocessing is required.

**MSigDB.**
Access to MSigDB gene set data requires registration at [https://www.gsea-msigdb.org/gsea/login.jsp](https://www.gsea-msigdb.org/gsea/login.jsp). Download the “Human Gene Set JSON file set (ZIP)”, extract the archive, and locate the full MSigDB annotation JSON (e.g. `msigdb.v2026.1.Hs.json`). Copy this file into the `local_dbs/` folder, then run the preprocessing step below to enable MSigDB-based tools:

```bash
cd src/tools/msigdb
python process_msigdb.py file_name
```


## 💻 Usage

To run Alvessa in UI mode:
```
# Launch the UI
alvessa ui
# Open http://127.0.0.1:8000
```

### CLI Commands

In addition to `alvessa ui`, the following commands are available:

```bash
# Ask a question
alvessa question "What pathways involve BRCA1 and what drugs target it?"

# List available tools
alvessa tools

# Run benchmarks
alvessa benchmark_all <folder> --shuffle --N <int> --save_intermediate --shuffle --restart <path.to.csv>
```

## 🔧 Contribute Tools

We welcome new tools! All tools are self-documenting with the structure `src/tools/<n>/node.py`. Register your tool in the tool catalog and Alvessa will automatically consider it during orchestration. See the [tutorials](tutorials/) for detailed information and requirements.

## 🏆 GenomeArena

A curated benchmark of **720 multiple-choice questions** spanning:

- Variant annotation
- Gene annotation  
- Pathways & interactions
- miRNA targets
- Drug-target relationships
- Protein structure
- Gene-phenotype associations

GenomeArena can be downloaded from [here](https://alvessa-public-access-bucket.s3.us-east-1.amazonaws.com/GenomeArena.zip) and scripts used to generate the questions are available [here](https://github.com/ksenia007/alvessa_agent/tree/main/evals/generation)

## 💡 Contribute Questions

Did you discover a good (or bad) example case? [Let us know!](https://github.com/ksenia007/alvessa_agent/issues)

We are also crowdsourcing a dataset of questions where the answer is known but requires multi-step reasoning. If you'd like to contribute evaluation cases, please contact us.

## 🗄️ Integrated Data Sources

<details>
<summary><strong>Click to expand full tool list</strong></summary>

**Variant Annotation:** dbSNP, ClinVar, AlphaMissense, Sei, ExpectoSC

**Gene Annotation:** GENCODE, UniProt, Alliance of Genome Resources, OMIM

**Pathways & Interactions:** Reactome, MSigDB, BioGRID, Gene Ontology, IntAct (viral)

**Regulatory:** ReMap, miRDB

**Drug & Druggability:** ChEMBL, DrugCentral, CysDB, OpenTargets

**Protein Structure:** AlphaFold, FPocket, FreeSASA, IUPred3, DisProt, BioLiP2

**Disease Associations:** GWAS Catalog, OpenTargets, ClinVar

</details>

## 📁 Outputs

Each run generates:

| File | Description |
|------|-------------|
| `demo.json` | Serialized pipeline state for UI |
| `demo.txt` | Human-readable answers with citations |
| `demo.log` | Full execution trace including tool calls |

As well as structured entity artifacts:

- 🧬 **Genes:** `genes/` — collected information per gene (`summary.txt`, `gene.json`, `transcripts.tsv`, `interactions_human.tsv`, `interactions_nonhuman.tsv`)
- 💊 **Drugs:** `drugs/` — per drug (`summary.txt`, `drug.json`)
- 🧬 **Variants:** `variants/` — CSVs for batch processing: `variants.csv` (index), `locations.csv`, `per_gene_traits.csv`, `per_gene_context.csv`, `functional_predictions.csv`, `allele_frequencies.csv`, `summaries.csv`

## 📄 Citation

If you use Alvessa in your research, please cite:

```bibtex
@article{sokolova2025alvessa,
  title   = {An Evidence-Grounded Research Assistant for Functional Genomics and Drug Target Assessment},
  author  = {Sokolova, Ksenia and Kosenkov, Dmitri and Nallamotu, Keerthana and Vedula, Sanketh and Sokolov, Daniil and Sapiro, Guillermo and Troyanskaya, Olga G},
  journal = {bioRxiv},
  year    = {2025},
  month   = dec,
  doi     = {10.64898/2025.12.30.697073},
  url     = {https://www.biorxiv.org/content/10.64898/2025.12.30.697073v1},
}
```
When results derived from specific resources are used in publications or downstream analyses, please also cite those primary databases directly, in accordance with their recommended citation policies.


### License
This project is released under the terms described in `LICENSE.pdf`.  
Please note that Alvessa uses data from third-party databases, each of which is governed by its own license and usage terms. Users are responsible for complying with the licenses and citation requirements of all external resources used in their analyses.


## 📧 Contact us
If you have any questions, please email us at `info@alvessa.ai`


---

<p align="center">
  <a href="https://github.com/ksenia007/alvessa_agent/issues">🐛 Report Bug</a> •
  <a href="https://github.com/ksenia007/alvessa_agent/issues">✨ Request Feature</a>
</p>
