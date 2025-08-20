"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-08-18
Updated: 2025-08-18


Description: 

Code to evaluate the base model performance."""

from __future__ import annotations
import json
import pandas as pd
from typing import Dict, Any
import os
import re
import time
from claude_client import claude_call
from config import CONDITIONED_MODEL

from run import run_pipeline

BASE_DIR = "benchmarks_generation"
RESULTS_DIR = os.path.join(BASE_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

system_msg = (
    "You are a multiple-choice answering system. "
    "You must only reply with one of the following single capital letters: A, B, C, or D. "
    "Do not add any words, punctuation, or explanation. "
    "Example valid output: 'C'. "
    "Example invalid outputs: 'Answer: C', 'C because...', 'Option C'."
)


def run_pipeline_alvessa(file_path: str, system_msg: str, max_rows: int = -1, file_save: str | None = None) -> pd.DataFrame:
    """
    Process a single CSV file and return results DataFrame from Alvessa.
    If `file_save` exists, previously answered questions are skipped.
    """
    df = pd.read_csv(file_path)

    # Optional cap & shuffle
    if max_rows > 0:
        df = df.sample(frac=1, random_state=42).reset_index(drop=True).head(max_rows)

    # Load existing results (resume)
    existing_df = None
    done_questions = set()
    if file_save and os.path.exists(file_save):
        try:
            existing_df = pd.read_csv(file_save)
            if "question" in existing_df.columns:
                done_questions = set(existing_df["question"].dropna().astype(str).unique())
            print(f"[resume] Loaded existing results: {len(existing_df)} rows, {len(done_questions)} unique questions.")
        except Exception as e:
            print(f"[resume] Warning: failed to load existing results at {file_save}: {e}")

    results = {}  # new results this run only

    # Helper for correctness
    def _is_correct(correct, model):
        return str(correct).strip().lower() == str(model).strip().lower()

    for i, row in df.iterrows():
        user_question = str(row["question"])
        if user_question in done_questions:
            # Already answered; skip
            print(f"[skip] Question already answered: {user_question}")
            continue

        print("\n" + "=" * 80)
        print("Q:", user_question)
        try:
            result = run_pipeline(user_question, prompt=system_msg, run_verifier=False)  
            answer = result["llm_json"].get("answer", "")
        except:
            print("[error] Exception during run_pipeline; skipping this question.")
            answer = ""
        

        results[i] = {
            "question": user_question,
            "correct_answer": row["answer"],
            "model_answer": answer,
            "is_correct": _is_correct(row["answer"], answer),
            "model": "Alvessa",
            "used_models": result["used_tools"]
        }

        print(f"Model answer: {answer}")
        print(f"Correct answer: {row['answer']}")

        # Save incremental progress
        if file_save:
            # Combine existing + new so far, de-duplicate by question (keep existing to avoid clobbering)
            new_df_so_far = pd.DataFrame.from_dict(results, orient="index")
            if existing_df is not None and not existing_df.empty:
                combined = pd.concat([existing_df, new_df_so_far], ignore_index=True)
                combined = combined.drop_duplicates(subset=["question"], keep="first")
            else:
                combined = new_df_so_far

            # Ensure is_correct present & consistent
            if "is_correct" not in combined.columns:
                combined["is_correct"] = combined.apply(
                    lambda x: _is_correct(x.get("correct_answer", ""), x.get("model_answer", "")), axis=1
                )

            combined.to_csv(file_save, index=False)

        time.sleep(10)  # throttle

    # Final assembly
    new_df = pd.DataFrame.from_dict(results, orient="index")

    # Ensure is_correct on new_df
    if not new_df.empty and "is_correct" not in new_df.columns:
        new_df["is_correct"] = new_df.apply(
            lambda x: _is_correct(x.get("correct_answer", ""), x.get("model_answer", "")),
            axis=1
        )

    if existing_df is not None and not existing_df.empty:
        final_df = pd.concat([existing_df, new_df], ignore_index=True)
        final_df = final_df.drop_duplicates(subset=["question"], keep="first")
    else:
        final_df = new_df

    # Summary
    if "is_correct" in final_df.columns and not final_df.empty:
        total_correct = int(final_df["is_correct"].sum())
        print(f"Correct answers (all in file): {total_correct} / {len(final_df)}")
    else:
        print("No results to summarize.")

    return final_df



def run_one_file(file_path: str, system_msg: str, max_rows: int = -1) -> pd.DataFrame:
    """
    Process a single CSV file and return results DataFrame.
    """
    df = pd.read_csv(file_path)
    results = {}
    
    if max_rows > 0:
        # shuffle and then head
        df = df.sample(frac=1, random_state=42).reset_index(drop=True)
        df = df.head(max_rows)

    for i, row in df.iterrows():
        user_question = row.question

        raw = claude_call(
            model=CONDITIONED_MODEL,
            temperature=0,
            max_tokens=2000,
            system=system_msg,
            messages=[{"role": "user", "content": f"User asked: {user_question}"}],
        )

        llm_resp = raw.content[0].text.strip() if hasattr(raw.content[0], "text") else raw.content[0]

        results[i] = {
            "question": user_question,
            "correct_answer": row.answer,
            "model_answer": llm_resp, 
            "model": CONDITIONED_MODEL
        }
        

    results_df = pd.DataFrame.from_dict(results, orient="index")
    results_df["is_correct"] = results_df.apply(
        lambda x: str(x["correct_answer"]).strip().lower() == str(x["model_answer"]).strip().lower(),
        axis=1
    )
    
    print(f"Correct answers: {results_df['is_correct'].sum()} / {len(results_df)}")
    print(f"Number of correct answers: {results_df['is_correct'].sum()} out of {len(results_df)}")
    
    return results_df



def run_eval_crawler(base_dir: str, results_dir: str, system_msg: str, function_name: str) -> None:
    """
    Evaluate all CSV files in the given base directory and save results.
    
    Parameters
    ----------
    base_dir : str
        Path to the base directory containing folders with CSV files.
    results_dir : str
        Path to the directory where results will be saved.
    system_msg : str
        System message for the Claude model.
    """
    all_results = []  # for combined results

    for folder in os.listdir(base_dir):
        print('***' * 10)
        print('Processing folder:', folder)
        folder_path = os.path.join(base_dir, folder)
        if not os.path.isdir(folder_path):
            continue
        if folder == "results":  # skip results folder itself
            continue

        # make a results subfolder for each dataset
        folder_results_dir = os.path.join(results_dir, folder)
        os.makedirs(folder_results_dir, exist_ok=True)

        for file in os.listdir(folder_path):
            print('Processing file:', file)
            if not file.endswith(".csv"):
                continue

            file_path = os.path.join(folder_path, file)
            results_df = run_one_file(file_path, system_msg)

            # save per-file
            out_path = os.path.join(folder_results_dir, f"{file.replace('.csv', '')}_results.csv")
            results_df.to_csv(out_path, index=False)

            print(f"{folder}/{file}: {results_df['is_correct'].sum()} / {len(results_df)} correct")

            # add to combined results
            all_results.append(results_df)

    # save combined
    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        combined_df.to_csv(os.path.join(results_dir, "all_results.csv"), index=False)
        print("Combined results saved.")

if __name__ == "__main__":
    # run_eval_crawler(BASE_DIR, RESULTS_DIR, system_msg)
    # print("Evaluation completed.")
    
    # One-off run for a specific file
    # # One-off run for a specific file
    # LOC_PATH = "benchmarks_generation/labbench/"
    # file_path = LOC_PATH + 'litqa_mc.csv' #"dbqa_mc.csv"
    # LOC_SAVE_PATH = "benchmarks_generation/results/labbench/alvessa"
    # # results_df = run_pipeline_alvessa(file_path, system_msg, max_rows = 50, 
    # #                                   file_save = os.path.join(LOC_SAVE_PATH, "dbqa_mc_results_subset.csv"))
    # results_df = run_pipeline_alvessa(file_path, system_msg, max_rows = 50, 
    #                                   file_save = os.path.join(LOC_SAVE_PATH, "litqa_mc_results_subset.csv"))
    # os.makedirs(LOC_SAVE_PATH, exist_ok=True)
    # #results_df.to_csv(os.path.join(LOC_SAVE_PATH, "dbqa_mc_results_subset.csv"), index=False)
    # results_df.to_csv(os.path.join(LOC_SAVE_PATH, "litqa_mc_results_subset.csv"), index=False)
    
    
    
    
    # ALVESSA ON BIOGRID
    # file_path =  "benchmarks_generation/biogrid/set2.csv"
    # LOC_SAVE_PATH = "benchmarks_generation/results/biogrid/alvessa"
    # os.makedirs(LOC_SAVE_PATH, exist_ok=True)
    # results_df = run_pipeline_alvessa(file_path, system_msg, max_rows = 50, 
    #                                   file_save = os.path.join(LOC_SAVE_PATH, "set2.csv"))
    # results_df.to_csv(os.path.join(LOC_SAVE_PATH, "set2.csv"), index=False)
    
    
    # # CLAUDE on BIOGRID
    # file_path =  "benchmarks_generation/biogrid/set2.csv"
    # LOC_SAVE_PATH = "benchmarks_generation/results/biogrid/claude"
    # os.makedirs(LOC_SAVE_PATH, exist_ok=True)
    # results_df = run_one_file(file_path, system_msg, max_rows = 50)
    # results_df.to_csv(os.path.join(LOC_SAVE_PATH, "set2.csv"), index=False)
    
    # CLAUDE on Reactome
    # file_path =  "benchmarks_generation/reactome/set1.csv"
    # LOC_SAVE_PATH = "benchmarks_generation/results/reactome/claude"
    # os.makedirs(LOC_SAVE_PATH, exist_ok=True)
    # results_df = run_one_file(file_path, system_msg, max_rows = 50)
    # results_df.to_csv(os.path.join(LOC_SAVE_PATH, "set1.csv"), index=False)
    
    # Alvessa on Reactome
    file_path =  "benchmarks_generation/reactome/set1.csv"
    LOC_SAVE_PATH = "benchmarks_generation/results/reactome/alvessa"
    os.makedirs(LOC_SAVE_PATH, exist_ok=True)
    results_df = run_pipeline_alvessa(file_path, system_msg, max_rows = 50, 
                                       file_save = os.path.join(LOC_SAVE_PATH, "set1.csv"))
    results_df.to_csv(os.path.join(LOC_SAVE_PATH, "set1.csv"), index=False)
    
    
    

    
    
    # file_path = LOC_PATH + "litqa_mc.csv"
    # results_df = run_one_file(file_path, system_msg, max_rows = 50)
    # # save files into new results folder: benchmarks_generation/results/labbench/claude
    # LOC_SAVE_PATH = "benchmarks_generation/results/labbench/claude"
    # os.makedirs(LOC_SAVE_PATH, exist_ok=True)
    # results_df.to_csv(os.path.join(LOC_SAVE_PATH, "litqa_mc_results_subset.csv"), index=False)
    
    # file_path = LOC_PATH + "dbqa_mc.csv"
    # results_df = run_one_file(file_path, system_msg, max_rows = 50)
    # # save files into new results folder: benchmarks_generation/results/labbench/claude
    # LOC_SAVE_PATH = "benchmarks_generation/results/labbench/claude"
    # os.makedirs(LOC_SAVE_PATH, exist_ok=True)
    # results_df.to_csv(os.path.join(LOC_SAVE_PATH, "dbqa_mc_results_subset.csv"), index=False)
   
    

