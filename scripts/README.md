# Scripts

Scripts for reproducing benchmark evaluations and figures in the paper.

## Benchmark Evaluation

### `entity_extraction.py`
Evaluates entity extraction on genes, variants, drugs, and miRNAs.

```bash
python entity_extraction.py --output evals/results/entity_extraction_results.json
```

### `report_entity_extraction_results.py`
Generates summary statistics and plots from entity extraction results.

```bash
python report_entity_extraction_results.py --input evals/results/entity_extraction_results.json --plot
```

## Visualization

### GenomeArena Benchmark

`visualize_benchmarks.py` - Overall accuracy and per-database breakdowns
```bash
python visualize_benchmarks.py
```

`visualize_tool_selection.py` - Tool selection heatmap

### dbQA Benchmark

`visualize_benchmarks_dbqa.py` - Accuracy by question database
```bash
python visualize_benchmarks_dbqa.py
```

### Runtime Analysis

`visualize_benchmarks_runtime.py` - Runtime comparison across models
```bash
python visualize_benchmarks_runtime.py --theme white
```

## Data Preparation

`subset_dbQA.py` - Create random subsets of dbQA questions
```bash
python subset_dbQA.py --n --seed 
```

`subset_ga.py` - Create 10% subset of GenomeArena questions

`plot_genomearena_counts.py` - Question distribution across databases

## Configuration

Edit `MODEL_FILES` dictionaries in visualization scripts to specify benchmark CSV paths for comparison.

Most scripts generate both white and black background variants for figures, saved to `results/benchmark_figures/`.
