"""Take 10% of the questions in a GA dataset to create a smaller subset for benchmarking."""
# %%
import pandas as pd
import random
import os
from pathlib import Path
LOC = Path("/Users/sokolova/Documents/research/alvessa_agent/benchmarks_generation/questions/GenomeArena")
# loop though subfolders and folders, find all csv files, shuffle and take 10% of rows, save to a new file
subset_ratio = 0.1
FILE_SAVE = Path("/Users/sokolova/Documents/research/alvessa_agent/benchmarks_generation/questions/subset_ga_questions.csv")

all_files = []
for root, dirs, files in os.walk(LOC):
    for file in files:
        if file.endswith(".csv"):
            all_files.append(Path(root) / file)

total_quesitons = 0
subset_dfs = []
for file in all_files:
    df = pd.read_csv(file)
    total_quesitons += len(df)
    subset_df = df.sample(frac=subset_ratio, random_state=42, replace=False)
    print(f"Original questions: {len(df)}, Subset questions: {len(subset_df)}")
    subset_dfs.append(subset_df)
print(f"Total questions in original GA dataset: {total_quesitons}")
final_subset_df = pd.concat(subset_dfs, ignore_index=True)
print(f"Total questions in subset GA dataset: {len(final_subset_df)}")
os.makedirs(FILE_SAVE.parent, exist_ok=True)
final_subset_df.to_csv(FILE_SAVE, index=False)

# %%
final_subset_df.question.value_counts()
# %%
