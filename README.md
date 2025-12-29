<p align="center">
  <img src="https://github.com/user-attachments/assets/2bdd5321-8ffe-44ae-aaff-be11a715bdc5" alt="Alvessa" width="1100"/>
</p>
<h1 align="center">
  <strong>An Evidence-Grounded Research Assistant for Functional Genomics and Drug Target Assessment</strong>
</h1>

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

Alvessa is a multi-agent framework that provides **verifiable, evidence-grounded answers** to genomics questions. Unlike general-purpose LLMs that can hallucinate identifiers or fabricate specifics, Alvessa enforces statement-level verification against retrieved database records, with feedback and explanations for each claim.

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
conda activate agents
pip install -e .

# Set API keys
export ANTHROPIC_API_KEY="your-key"
export BIOGRID_API_KEY="your-key"        # optional
export DISGENET_API_KEY="your-key"       # optional

# Launch the UI
alvessa ui 8000
# Open http://127.0.0.1:8000
```

**Requirements:** Python 3.10+, local database files in `local_dbs/`

## 💻 Usage

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
@article{TBD}
```

## 📜 License

[License details here]

---

<p align="center">
  <a href="https://github.com/ksenia007/alvessa_agent/issues">🐛 Report Bug</a> •
  <a href="https://github.com/ksenia007/alvessa_agent/issues">✨ Request Feature</a>
</p>
