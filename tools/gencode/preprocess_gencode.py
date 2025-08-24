# %%
import pandas as pd

gtf_path = "../../local_dbs/gencode.v48.primary_assembly.basic.annotation.gtf.gz"

# Load GTF into DataFrame
gtf_df = pd.read_csv(
    gtf_path,
    sep="\t", comment="#", header=None,
    names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
)

# Parse attributes into columns
def parse_attributes(attr_str):
    attrs = {}
    for item in attr_str.strip().split(";"):
        if item.strip():
            key, value = item.strip().replace('"', "").split(" ", 1)
            attrs[key] = value
    return attrs

attr_expanded = gtf_df["attributes"].apply(parse_attributes)
attr_df = pd.DataFrame(attr_expanded.tolist())

# Merge expanded attributes with original GTF data
gtf_df_expanded = pd.concat([gtf_df.drop(columns=["attributes"]), attr_df], axis=1)

# Save preprocessed version
out_path = "../../local_dbs/gencode.v48.expanded.parquet"
# gtf_df_expanded.to_parquet(out_path, index=False)

print(f"Expanded GTF saved to {out_path} with {gtf_df_expanded.shape[0]} rows.")
print("Available attribute columns:", list(attr_df.columns))

# %%
