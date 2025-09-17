# %%
import json
import random
from typing import List
import pandas as pd
# %%
def construct_question_mc(df: pd.Series):
    """
    Constructs a multiple-choice question from a DataFrame row.
    Answer is in "ideal" column, and distractors are in "distractors" column.
    The order of answers is shuffled.

    Returns
    -------
    tuple[str, str]
        (question_text, correct_letter)
    """
    q = df["question"]
    ideal = df["ideal"]
    distractors = df["distractors"]

    # If distractors come in as a string like "[MNX1, MCTP1, DMXL1]"
    if isinstance(distractors, str):
        distractors = distractors.strip("[]").split(",")
    distractors = [d.strip() for d in distractors if d.strip()]

    if not distractors:
        raise ValueError("No distractors found in the input data.")

    # Combine and shuffle
    options = distractors + [ideal]
    random.shuffle(options)

    labels = ["A", "B", "C", "D"]
    lines = [q, ""]
    correct_letter = None
    for label, opt in zip(labels, options):
        lines.append(f"{label}. {opt}")
        if opt == ideal:
            correct_letter = label

    return "\n".join(lines), correct_letter

def convert_to_csv(df: pd.DataFrame, outfile: str):
    """
    Given a DataFrame with columns ['question', 'ideal', 'distractors'],
    build multiple-choice questions and save them to a CSV.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain 'question', 'ideal', 'distractors'.
    outfile : str
        Path to output CSV file.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['question', 'answer'] that was written to CSV.
    """
    records = []
    for _, row in df.iterrows():
        q_text, correct_letter = construct_question_mc(row)
        records.append({"question": q_text, "answer": correct_letter})

    result = pd.DataFrame(records)
    result.to_csv(outfile, index=False)
    return result


# file = 'dbqa.parquet' 
# questions = pd.read_parquet(file) 
# convert_to_csv(questions, 'dbqa_mc.csv')
file = 'litqa.parquet' 
questions = pd.read_parquet(file) 
convert_to_csv(questions, 'litqa_mc.csv')


