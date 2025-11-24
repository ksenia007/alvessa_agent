import argparse
import random
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LOCAL_DBS = ROOT / "local_dbs"


def load_gwas_data(file_path: Path) -> pd.DataFrame:
    """Load GWAS catalogue data with variant columns."""
    print(f"Loading GWAS data from {file_path}...")

    columns_of_interest = [
        "SNPS",
        "MAPPED_GENE",
        "CHR_ID",
        "UPSTREAM_GENE_ID",
        "DOWNSTREAM_GENE_ID",
    ]

    try:
        chunks = []
        chunk_size = 10_000
        for chunk in pd.read_csv(
            file_path,
            sep="\t",
            chunksize=chunk_size,
            usecols=columns_of_interest,
            low_memory=False,
        ):
            chunks.append(chunk)

        df = pd.concat(chunks, ignore_index=True)
        print(f"Loaded {len(df)} rows")
        return df
    except Exception as exc:
        print(f"Error loading data: {exc}", file=sys.stderr)
        sys.exit(1)


def clean_and_expand(df: pd.DataFrame) -> pd.DataFrame:
    """Clean and expand SNP/gene mappings; keep one SNP-gene per row."""

    df = df.dropna(subset=["SNPS", "MAPPED_GENE", "CHR_ID"])

    # Normalize to strings
    df = df.copy()
    df["SNPS"] = df["SNPS"].astype(str)
    df["MAPPED_GENE"] = df["MAPPED_GENE"].astype(str)
    df["CHR_ID"] = df["CHR_ID"].astype(str)

    # Expand SNPs if multiple per row (split on comma/semicolon/plus)
    def split_field(val: str) -> list[str]:
        parts = []
        for token in str(val).replace(";", ",").replace("+", ",").split(","):
            token = token.strip()
            if token:
                parts.append(token)
        return parts

    rows = []
    for _, row in df.iterrows():
        snps = split_field(row["SNPS"])
        genes = [g.strip() for g in row["MAPPED_GENE"].replace(";", ",").split(",") if g.strip()]
        chr_tokens = [c.strip() for c in row["CHR_ID"].replace(";", ",").split(",") if c.strip()]
        if not genes or not chr_tokens:
            continue

        for snp in snps:
            for idx, gene in enumerate(genes):
                if " - " in gene or " x " in gene or "-.-" in gene:
                    continue
                chr_id = chr_tokens[idx] if idx < len(chr_tokens) else chr_tokens[0]
                rows.append({"SNP": snp, "GENE": gene, "CHR_ID": chr_id})

    expanded = pd.DataFrame(rows)
    if expanded.empty:
        print("No SNP-gene mappings after cleaning", file=sys.stderr)
        return expanded

    # Drop duplicates
    expanded = expanded.drop_duplicates()
    print(f"After cleaning/expansion: {len(expanded)} SNP-gene pairs")
    return expanded


def build_indices(expanded: pd.DataFrame):
    """Build lookup structures for SNP->genes and chromosome gene lists."""
    snp_to_genes: dict[str, set[str]] = {}
    chr_to_genes: dict[str, set[str]] = {}

    for _, row in expanded.iterrows():
        snp = row["SNP"]
        gene = row["GENE"]
        chr_id = row["CHR_ID"]

        snp_to_genes.setdefault(snp, set()).add(gene)
        chr_to_genes.setdefault(chr_id, set()).add(gene)

    return snp_to_genes, chr_to_genes


def generate_variant_questions(
    count: int,
    gwas_file: Path,
    out_dir: Path,
    seed: int = 42,
):
    rng = random.Random(seed)

    df = load_gwas_data(gwas_file)
    expanded = clean_and_expand(df)
    if expanded.empty:
        print("Nothing to generate for set3", file=sys.stderr)
        return

    snp_to_genes, chr_to_genes = build_indices(expanded)
    snps = [s for s, genes in snp_to_genes.items() if genes]
    rng.shuffle(snps)

    questions = []

    for snp in snps:
        if len(questions) >= count:
            break

        mapped_genes = list(snp_to_genes[snp])
        chr_id = expanded.loc[expanded["SNP"] == snp, "CHR_ID"].iloc[0]

        # Correct
        correct_gene = rng.choice(mapped_genes)

        # Distractors
        same_chr_pool = [g for g in chr_to_genes.get(chr_id, set()) if g not in snp_to_genes[snp]]
        other_genes_pool = [g for g in set(expanded["GENE"]) if g not in snp_to_genes[snp]]

        if not same_chr_pool or len(other_genes_pool) < 2:
            continue

        rng.shuffle(same_chr_pool)
        rng.shuffle(other_genes_pool)

        distractors = [same_chr_pool[0], other_genes_pool[0], other_genes_pool[1]]

        options = [correct_gene] + distractors
        rng.shuffle(options)
        correct_index = options.index(correct_gene)
        correct_letter = chr(ord("A") + correct_index)

        question_text = (
            f"Which gene is variant {snp} mapped to according to the GWAS catalog? "
            f"[A] {options[0]} [B] {options[1]} [C] {options[2]} [D] {options[3]}"
        )

        questions.append(
            {
                "question": question_text,
                "answer": correct_letter,
                "tool": "GWAS",
                "snp": snp,
                "chr_id": chr_id,
                "correct_gene": correct_gene,
            }
        )

    if not questions:
        print("No questions generated for set3; insufficient distractors", file=sys.stderr)
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "set3.csv"
    pd.DataFrame(questions).to_csv(out_path, index=False)
    print(f"Wrote {len(questions)} set3 questions to {out_path}")


def parse_args():
    parser = argparse.ArgumentParser(description="Generate GWAS variant->gene MCQ set (set3)")
    parser.add_argument("--count", type=int, default=20, help="Number of questions to generate")
    parser.add_argument(
        "--db-dir",
        type=Path,
        default=DEFAULT_LOCAL_DBS,
        help="Directory containing gwas_catalogue_association.tsv",
    )
    parser.add_argument(
        "--gwas-file",
        type=Path,
        default=None,
        help="Path to GWAS TSV (overrides --db-dir)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "benchmark_questions" / "GWAS",
        help="Directory to write set3.csv",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    gwas_file = args.gwas_file
    if gwas_file is None:
        gwas_file = Path(args.db_dir) / "gwas_catalogue_association.tsv"

    if not gwas_file.exists():
        print(f"GWAS file not found at {gwas_file}", file=sys.stderr)
        raise SystemExit(1)

    generate_variant_questions(
        count=args.count,
        gwas_file=gwas_file,
        out_dir=Path(args.out_dir),
        seed=args.seed,
    )
