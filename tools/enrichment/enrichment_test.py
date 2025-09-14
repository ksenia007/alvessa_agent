import pandas as pd
from scipy.stats import hypergeom
import statsmodels.stats.multitest as smm
from collections import defaultdict
from scipy.stats import hypergeom, fisher_exact

def parse_obo_names(obo_file):
    go_names = {}
    current_id, current_name = None, None
    with open(obo_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("id: GO:"):
                current_id = line.split("id: ")[1]
            elif line.startswith("name:") and current_id:
                current_name = line.split("name: ")[1]
                go_names[current_id] = current_name
                current_id, current_name = None, None
    return go_names


def run_go_enrichment(gene_set, gmt_file, obo_file,
                      method="fisher", fdr_threshold=0.05, min_overlap=2):

    go_terms = {}
    all_genes = set()
    with open(gmt_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            term = parts[0] 
            genes = set(parts[2:]) 
            go_terms[term] = genes
            all_genes.update(genes)

    go_names = {}
    if obo_file:
        go_names = parse_obo_names(obo_file)

    gene_list = set(gene_set) & all_genes
    M = len(all_genes)   
    n = len(gene_list)  

    results = []
    for term, term_genes in go_terms.items():
        N = len(term_genes) 
        k = len(gene_list & term_genes) 

        if method == "fisher":
            if k == 0: 
                continue 

            a = k
            b = n - k
            c = N - k
            d = M - N - b 
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")

        elif method == "hypergeom":
            if k < min_overlap:
                continue

            p_value = hypergeom.pmf(k, M, N, n)

        else:
            raise ValueError("Method must be 'fisher' or 'hypergeom'")

        results.append({
            "GO_Term": term,
            "GO_Name": go_names.get(term, "NA"),
            "Overlap_Count": k,
            "GO_Term_Gene_Count": N,
            "Input_Gene_Count": n,
            "P-value": p_value,
            "Genes": ", ".join(sorted(gene_list & term_genes)),
        })

    df = pd.DataFrame(results)
    if df.empty:
        print("No enriched terms found.")
        return df

    df["FDR"] = smm.multipletests(df["P-value"], method="fdr_bh")[1]
    df = df.sort_values("FDR")

    sig_df = df[df["FDR"] < fdr_threshold]

    # print(f"Significant GO terms (FDR < {fdr_threshold}): {len(sig_df)}")

    return sig_df
