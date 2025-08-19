"""
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-08-18
Updated: 2025-08-18


Description: 

Code to evaluate the base model performance."""
# %%
from __future__ import annotations
import json
import pandas as pd
from typing import Dict, Any
import os
import re

from claude_client import claude_call
from config import CONDITIONED_MODEL

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

all_results = []  # for combined results

for folder in os.listdir(BASE_DIR):
    print('***' * 10)
    print('Processing folder:', folder)
    folder_path = os.path.join(BASE_DIR, folder)
    if not os.path.isdir(folder_path):
        continue
    if folder == "results":  # skip results folder itself
        continue

    # make a results subfolder for each dataset
    folder_results_dir = os.path.join(RESULTS_DIR, folder)
    os.makedirs(folder_results_dir, exist_ok=True)

    for file in os.listdir(folder_path):
        print('Processing file:', file)
        if not file.endswith(".csv"):
            continue

        file_path = os.path.join(folder_path, file)
        df = pd.read_csv(file_path)

        results = {}
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
                "folder": folder,
                "file": file
            }

        results_df = pd.DataFrame.from_dict(results, orient="index")
        results_df["is_correct"] = results_df.apply(
            lambda x: x["correct_answer"].strip().lower() == x["model_answer"].strip().lower(),
            axis=1
        )

        # save per-file
        out_path = os.path.join(folder_results_dir, f"{file.replace('.csv', '')}_results.csv")
        results_df.to_csv(out_path, index=False)

        print(f"{folder}/{file}: {results_df['is_correct'].sum()} / {len(results_df)} correct")

        # add to combined results
        all_results.append(results_df)

# save combined
if all_results:
    combined_df = pd.concat(all_results, ignore_index=True)
    combined_df.to_csv(os.path.join(RESULTS_DIR, "all_results.csv"), index=False)
    print("Combined results saved.")
