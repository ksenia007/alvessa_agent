#!/usr/bin/env python3
"""
Entity Extraction Evaluation Script

This script evaluates the performance of multiple entity extraction methods
by running them against benchmark datasets and comparing their results.

Author: Generated for evaluation
Created: 2025-01-27
Modified: 2025-08-21
"""

import csv
import os
import time
import json
from pathlib import Path
from typing import Dict, List, Tuple, Callable
from collections import defaultdict

# Import the entity extraction nodes and the State class
from src.alvessa.agents.entity_extraction import (
    entity_extraction_node,
    claude_entity_extraction_node,
    flair_entity_extraction_node,
    gliner_entity_extraction_node,
    claude_flair_entity_extraction_node,
    gliner_flair_entity_extraction_node,
    gliner_claude_entity_extraction_node,
)
from src.state import State


def load_benchmark_data(benchmark_dir: str) -> Dict[str, List[Tuple[str, List[str]]]]:
    """
    Load all benchmark CSV files from the 'entity_recognition' subdirectory.
    """
    benchmark_data = {}
    benchmark_path = Path(benchmark_dir)
    entity_recognition_path = benchmark_path / "entity_recognition"

    if not entity_recognition_path.exists():
        raise FileNotFoundError(f"Benchmark subdirectory not found: {entity_recognition_path}")

    csv_files = sorted(list(entity_recognition_path.glob("*.csv")))
    if not csv_files:
        print(f"Warning: No benchmark CSV files found in {entity_recognition_path}")

    for csv_file in csv_files:
        set_name = csv_file.stem
        queries_and_genes = []
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                query = row['query'].strip()
                expected_genes_str = row.get('recognized_genes', '').strip()
                expected_genes = [g.strip() for g in expected_genes_str.split(',') if g.strip()]
                queries_and_genes.append((query, expected_genes))
        benchmark_data[set_name] = queries_and_genes
        print(f"Loaded {len(queries_and_genes)} queries from {set_name}")
    return benchmark_data


def evaluate_single_query(query: str, expected_genes: List[str], extraction_fn: Callable) -> Dict:
    """
    Evaluate entity extraction for a single query using a specified extraction function.
    """
    state = State({"messages": [{"role": "user", "content": query}]})
    try:
        result = extraction_fn(state)
        extracted_genes = result.get("genes", [])

        # Case-insensitive comparison
        expected_set = {gene.upper() for gene in expected_genes}
        extracted_set = {gene.upper() for gene in extracted_genes}

        true_positives = len(expected_set & extracted_set)
        false_positives = len(extracted_set - expected_set)
        false_negatives = len(expected_set - extracted_set)

        precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
        recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
        f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        
        return {
            "query": query, "expected_genes": expected_genes, "extracted_genes": extracted_genes,
            "true_positives": true_positives, "false_positives": false_positives, "false_negatives": false_negatives,
            "precision": precision, "recall": recall, "f1_score": f1_score,
            "perfect_match": expected_set == extracted_set, "success": True, "error": None,
        }
    except Exception as e:
        print(f"ERROR processing query '{query[:50]}...': {e}")
        return {
            "query": query, "expected_genes": expected_genes, "extracted_genes": [],
            "true_positives": 0, "false_positives": 0, "false_negatives": len(expected_genes),
            "precision": 0, "recall": 0, "f1_score": 0,
            "perfect_match": False, "success": False, "error": str(e),
        }


def calculate_aggregate_metrics(results: List[Dict]) -> Dict:
    """Calculate aggregate metrics (macro and micro averages) across all results."""
    successful_results = [r for r in results if r["success"]]
    if not successful_results:
        return {
            "total_queries": len(results), "successful_queries": 0, "error_rate": 1.0, "perfect_match_rate": 0.0,
            "macro_precision": 0.0, "macro_recall": 0.0, "macro_f1": 0.0,
            "micro_precision": 0.0, "micro_recall": 0.0, "micro_f1": 0.0, "perfect_matches": 0
        }

    total_tp = sum(r["true_positives"] for r in successful_results)
    total_fp = sum(r["false_positives"] for r in successful_results)
    total_fn = sum(r["false_negatives"] for r in successful_results)

    micro_precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    micro_recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
    micro_f1 = 2 * (micro_precision * micro_recall) / (micro_precision + micro_recall) if (micro_precision + micro_recall) > 0 else 0

    perfect_matches = sum(1 for r in successful_results if r["perfect_match"])
    
    return {
        "total_queries": len(results), "successful_queries": len(successful_results),
        "error_rate": (len(results) - len(successful_results)) / len(results),
        "perfect_match_rate": perfect_matches / len(successful_results),
        "perfect_matches": perfect_matches,
        "macro_precision": sum(r["precision"] for r in successful_results) / len(successful_results),
        "macro_recall": sum(r["recall"] for r in successful_results) / len(successful_results),
        "macro_f1": sum(r["f1_score"] for r in successful_results) / len(successful_results),
        "micro_precision": micro_precision, "micro_recall": micro_recall, "micro_f1": micro_f1,
    }


def run_evaluation(benchmark_dir: str, extraction_fn: Callable, output_file: str, max_queries: int = None) -> Dict:
    """Run the complete evaluation pipeline for a given extraction function."""
    print("Loading benchmark data...")
    benchmark_data = load_benchmark_data(benchmark_dir)
    all_results = {}
    total_queries_processed = 0

    for set_name, queries_and_genes in benchmark_data.items():
        print(f"\nEvaluating {set_name}...")
        set_results = []
        queries_to_run = queries_and_genes[:max_queries - total_queries_processed] if max_queries else queries_and_genes

        for i, (query, expected_genes) in enumerate(queries_to_run):
            if i % 20 == 0:
                print(f"  Processing query {i+1}/{len(queries_to_run)}...")
            result = evaluate_single_query(query, expected_genes, extraction_fn)
            set_results.append(result)
        
        total_queries_processed += len(queries_to_run)
        set_metrics = calculate_aggregate_metrics(set_results)
        all_results[set_name] = {"queries": set_results, "metrics": set_metrics}
        print(f"  {set_name} completed. Micro F1: {set_metrics['micro_f1']:.3f}")

        if max_queries and total_queries_processed >= max_queries:
            break

    overall_metrics = calculate_aggregate_metrics([q for r in all_results.values() for q in r["queries"]])
    final_results = {
        "evaluation_summary": {"overall_metrics": overall_metrics, "evaluation_time": time.strftime("%Y-%m-%d %H:%M:%S")},
        "results_by_set": all_results
    }

    with open(output_file, 'w') as f:
        json.dump(final_results, f, indent=2)
    print(f"\nDetailed results saved to: {output_file}")
    return final_results


def print_evaluation_summary(results: Dict, model_name: str):
    """Print a formatted summary of evaluation results for a single model."""
    summary = results["evaluation_summary"]
    overall = summary["overall_metrics"]
    
    print("\n" + "="*60)
    print(f"EVALUATION SUMMARY: {model_name.upper()}")
    print("="*60)
    print(f"Total queries processed: {overall['total_queries']}")
    print(f"Perfect Match Rate: {overall['perfect_match_rate']:.2%} ({overall['perfect_matches']}/{overall['successful_queries']})")
    print(f"Micro-averaged F1-Score: {overall['micro_f1']:.3f}")
    print(f"Micro-averaged Precision: {overall['micro_precision']:.3f}")
    print(f"Micro-averaged Recall: {overall['micro_recall']:.3f}")
    print("-" * 60)


def print_comparative_summary(all_metrics: Dict):
    """Prints a final comparative table of all evaluated models."""
    print("\n" + "="*70)
    print("COMPARATIVE EVALUATION SUMMARY")
    print("="*70)
    # Header
    print(f"{'Model':<30} | {'Micro F1':>10} | {'Precision':>10} | {'Recall':>10}")
    print("-" * 70)
    
    # Rows
    for model_name, metrics in all_metrics.items():
        f1 = f"{metrics['micro_f1']:.3f}"
        precision = f"{metrics['micro_precision']:.3f}"
        recall = f"{metrics['micro_recall']:.3f}"
        print(f"{model_name:<30} | {f1:>10} | {precision:>10} | {recall:>10}")
    print("="*70)


def main():
    """Main function to run evaluations for all specified models."""
    base_dir = Path(__file__).parent
    benchmark_dir = base_dir / "benchmark"
    output_dir = base_dir / "evaluation_results"
    os.makedirs(output_dir, exist_ok=True)

    # Use a small number for quick testing, set to None for the full benchmark
    MAX_QUERIES_PER_MODEL = None 

    nodes_to_evaluate = {
        "Merged (Claude+Flair+GLiNER)": entity_extraction_node,
        "Merged (Claude+Flair)": claude_flair_entity_extraction_node,
        "Merged (Claude+GLiNER)": gliner_claude_entity_extraction_node,
        "Merged (Flair+GLiNER)": gliner_flair_entity_extraction_node,
        "Claude Only": claude_entity_extraction_node,
        "Flair Only": flair_entity_extraction_node,
        "GLiNER Only": gliner_entity_extraction_node,
    }

    all_run_metrics = {}

    for name, node_fn in nodes_to_evaluate.items():
        print(f"\n\n{'#'*25} EVALUATING: {name.upper()} {'#'*25}")
        output_file = output_dir / f"evaluation_results_{name.replace(' ', '_').lower()}.json"
        
        results = run_evaluation(
            benchmark_dir=str(benchmark_dir),
            extraction_fn=node_fn,
            output_file=str(output_file),
            max_queries=MAX_QUERIES_PER_MODEL
        )
        
        print_evaluation_summary(results, name)
        all_run_metrics[name] = results["evaluation_summary"]["overall_metrics"]
    
    print_comparative_summary(all_run_metrics)
    print("\nEvaluation complete for all models.")


if __name__ == "__main__":
    main()
