import numpy as np
from sklearn.metrics import accuracy_score
from scipy.stats import bootstrap

def _compute_by_group(df: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        df.groupby("question_group")["is_correct"]
        .agg(["mean", "count"])
        .rename(columns={"mean": "accuracy", "count": "n"})
        .reset_index()
    )
    overall = pd.DataFrame({"question_group": ["All"], "accuracy": [df["is_correct"].mean()], "n": [len(df)]})
    full = pd.concat([overall, grouped], ignore_index=True)
    return full


# Combine true and predicted labels into a single dataset for consistent resampling
data = (y_true, y_pred)

def calculate_accuracy(data, axis):
    true_labels = data[0, axis]
    predicted_labels = data[1, axis]
    return accuracy_score(true_labels, predicted_labels)

def compute_accuracy_ci(y_true, y_pred, n_resamples=1000, confidence_level=0.95):
    data = (y_true, y_pred)

    results = bootstrap(
        data,
        calculate_accuracy,
        vectorized=True,
        paired=True,
        axis=-1,
        n_resamples=n_resamples,
        confidence_level=confidence_level,
        method='percentile'
    )

    lower_bound = results.confidence_interval.low
    upper_bound = results.confidence_interval.high
    return lower_bound, upper_bound 

def compute_ci(df: pd.DataFrame) -> pd.DataFrame:
    

# Perform bootstrapping
# n_resamples is typically a large number (e.g., 1000 or 5000)
# confidence_level is 0.95 for a 95% CI
results = bootstrap(
    data,
    calculate_accuracy,
    vectorized=True, # Set to True as our function handles vectorized input
    paired=True,     # Use paired=True since we are resampling corresponding true and pred pairs
    axis=-1,         # Resample along the last axis (the individual samples)
    n_resamples=1000,
    confidence_level=0.95,
    method='percentile' # Use the percentile method for CI calculation
)

# Extract the confidence interval bounds
lower_bound = results.confidence_interval.low
upper_bound = results.confidence_interval.high