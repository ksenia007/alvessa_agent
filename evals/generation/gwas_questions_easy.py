import argparse
import random
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

DEFAULT_LOCAL_DBS = ROOT / "local_dbs"


def load_gwas_data(file_path: Path) -> pd.DataFrame:
    """Load GWAS catalogue data from TSV file."""
    print(f"Loading GWAS data from {file_path}...")

    columns_of_interest = ["DISEASE/TRAIT", "MAPPED_GENE", "P-VALUE", "OR or BETA"]

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


def clean_data(df: pd.DataFrame) -> pd.DataFrame:
    """Clean strings and numeric columns; keep all associations (significant + limbo)."""
    print("Cleaning and preparing data...")

    df = df.dropna(subset=["DISEASE/TRAIT", "MAPPED_GENE"])

    df = df.copy()
    df["MAPPED_GENE"] = df["MAPPED_GENE"].astype(str).str.split(",").str[0].str.strip()
    df["MAPPED_GENE"] = df["MAPPED_GENE"].str.split(";").str[0].str.strip()
    df = df[~df["MAPPED_GENE"].str.contains(" - ")]  # drop composite/hyphenated gene listings
    df = df[~df["MAPPED_GENE"].str.contains(" x ")]  # drop compound listings with x
    df = df[~df["MAPPED_GENE"].str.contains("-.-")]  # drop genes with double dashes

    df["DISEASE/TRAIT"] = df["DISEASE/TRAIT"].astype(str).str.strip()

    df["P-VALUE"] = pd.to_numeric(df["P-VALUE"], errors="coerce")
    df["OR or BETA"] = pd.to_numeric(df["OR or BETA"], errors="coerce")

    df = df[df["MAPPED_GENE"].str.len() > 1]
    df = df[df["DISEASE/TRAIT"].str.len() > 3]

    df = df.sort_values("P-VALUE").drop_duplicates(subset=["MAPPED_GENE", "DISEASE/TRAIT"], keep="first")

    print(f"After cleaning: {len(df)} rows")
    return df


def partition_gene_traits(df: pd.DataFrame, p_value_threshold: float = 5e-8) -> tuple[dict, list]:
    """Split traits per gene into significant vs limbo; return mapping and all traits list."""

    gene_traits: dict[str, dict[str, dict]] = {}
    for _, row in df.iterrows():
        gene = row["MAPPED_GENE"]
        trait = row["DISEASE/TRAIT"]
        pvalue = row["P-VALUE"]
        risk_score = row["OR or BETA"]

        bucket = "significant" if (pd.notna(pvalue) and pvalue < p_value_threshold) else "limbo"

        if gene not in gene_traits:
            gene_traits[gene] = {"significant": {}, "limbo": {}}

        # Keep the best (lowest p) entry per gene/trait within the bucket
        current = gene_traits[gene][bucket].get(trait)
        if current is None or (pd.notna(pvalue) and pd.notna(current.get("p_value")) and pvalue < current["p_value"]):
            gene_traits[gene][bucket][trait] = {"p_value": pvalue, "risk_score": risk_score}

    all_traits = df["DISEASE/TRAIT"].unique().tolist()
    print(f"Found {len(gene_traits)} genes with any association; {len(all_traits)} unique traits")
    return gene_traits, all_traits


def format_metric(value: float, small_fmt: str = "{:.2e}", regular_fmt: str = "{:.3f}") -> str:
    if pd.isna(value):
        return "Not available"
    if value < 0.001:
        return small_fmt.format(value)
    return regular_fmt.format(value)


def generate_question(gene: str, traits: dict, all_traits: list, rng: random.Random):
    """Generate a single multiple-choice question for a gene; return None when impossible."""

    significant = traits.get("significant", {})
    limbo = traits.get("limbo", {})
    if not significant:
        return None

    correct_trait = rng.choice(list(significant.keys()))
    trait_data = significant[correct_trait]

    excluded = set(significant) | set(limbo)
    distractor_pool = [t for t in all_traits if t not in excluded]
    if len(distractor_pool) < 3:
        return None

    incorrect_traits = rng.sample(distractor_pool, 3)

    all_options = [correct_trait] + incorrect_traits
    rng.shuffle(all_options)

    correct_index = all_options.index(correct_trait)
    correct_letter = chr(ord("A") + correct_index)

    p_value_str = format_metric(trait_data.get("p_value"))
    risk_score_str = format_metric(trait_data.get("risk_score"), regular_fmt="{:.2f}")

    return {
        "gene": gene,
        "question": f"Which of the following traits are associated with gene {gene}?",
        "options": all_options,
        "correct_answer": correct_letter,
        "correct_trait": correct_trait,
        "p_value": p_value_str,
        "risk_score": risk_score_str,
    }

def generate_gwas_set(
    count: int,
    gwas_file: Path,
    out_dir: Path,
    seed: int = 42,
    p_value_threshold: float = 5e-8,
):
    """Generate one GWAS MCQ set using only genome-wide significant positives and clean negatives."""

    rng = random.Random(seed)

    df = load_gwas_data(gwas_file)
    df = clean_data(df)
    gene_traits, all_traits = partition_gene_traits(df, p_value_threshold=p_value_threshold)

    eligible_genes = [g for g, t in gene_traits.items() if t.get("significant")]
    if not eligible_genes:
        print("No genes with significant associations found; nothing to generate", file=sys.stderr)
        return

    rng.shuffle(eligible_genes)

    questions = []
    for gene in eligible_genes:
        if len(questions) >= count:
            break

        qdata = generate_question(gene, gene_traits[gene], all_traits, rng)
        if qdata is None:
            print(f"Skipping gene {gene}: insufficient distractors or no significant trait")
            continue

        options = qdata["options"]
        question = (
            f"Which of the following traits are associated with gene {gene}? "
            f"[A] {options[0]} [B] {options[1]} [C] {options[2]} [D] {options[3]}"
        )

        questions.append(
            {
                "question": question,
                "answer": qdata["correct_answer"],
                "tool": "GWAS",
                "gene": gene,
                "correct_trait": qdata["correct_trait"],
                "p_value": qdata["p_value"],
                "risk_score": qdata["risk_score"],
            }
        )

    if not questions:
        print("No questions generated; check data and thresholds", file=sys.stderr)
        return

    output_df = pd.DataFrame(questions)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "set.csv"
    output_df.to_csv(out_path, index=False)
    if len(questions) < count:
        print(
            f"Only generated {len(questions)} questions (requested {count}); "
            "not enough eligible genes with clean distractors",
            file=sys.stderr,
        )
    print(f"Wrote {len(questions)} questions to {out_path}")


def generate_gwas_set2(
    count: int,
    gwas_file: Path,
    out_dir: Path,
    seed: int = 42,
    p_value_threshold: float = 5e-8,
    max_traits_per_gene: int = 500,
    max_questions_per_gene: int = 1,
):
    """Generate Set 2: trait associated with exactly one of two genes (not both)."""

    rng = random.Random(seed)

    df = load_gwas_data(gwas_file)
    df = clean_data(df)
    gene_traits, all_traits = partition_gene_traits(df, p_value_threshold=p_value_threshold)

    gene_sig = {
        g: set(traits.get("significant", {}).keys())
        for g, traits in gene_traits.items()
        if traits.get("significant") and len(traits.get("significant", {})) <= max_traits_per_gene
    }
    gene_limbo = {g: set(traits.get("limbo", {}).keys()) for g, traits in gene_traits.items()}
    if not gene_sig:
        print("No genes with significant associations found; nothing to generate", file=sys.stderr)
        return

    trait_to_genes: dict[str, list[str]] = {}
    for gene, traits in gene_sig.items():
        for trait in traits:
            trait_to_genes.setdefault(trait, []).append(gene)

    traits_with_genes = list(trait_to_genes.items())
    rng.shuffle(traits_with_genes)

    questions = []
    used_pairs: set[frozenset[str]] = set()
    gene_usage: dict[str, int] = {}

    def eligible_exclusive(pos_gene: str, neg_gene: str) -> list[str]:
        return [
            t
            for t in gene_sig.get(pos_gene, set())
            if t not in gene_sig.get(neg_gene, set()) and t not in gene_limbo.get(neg_gene, set())
        ]

    for trait, genes_with_trait in traits_with_genes:
        if len(questions) >= count:
            break
        if len(genes_with_trait) < 2:
            continue

        genes = list(genes_with_trait)
        rng.shuffle(genes)

        for i in range(len(genes)):
            if len(questions) >= count:
                break
            for j in range(i + 1, len(genes)):
                if len(questions) >= count:
                    break
                gene_a, gene_b = genes[i], genes[j]
                pair_key = frozenset({gene_a, gene_b})
                if pair_key in used_pairs:
                    continue
                if gene_usage.get(gene_a, 0) >= max_questions_per_gene:
                    continue
                if gene_usage.get(gene_b, 0) >= max_questions_per_gene:
                    continue

                exclusives_a = eligible_exclusive(gene_a, gene_b)
                exclusives_b = eligible_exclusive(gene_b, gene_a)

                directions = []
                if exclusives_a:
                    directions.append((gene_a, gene_b, exclusives_a))
                if exclusives_b:
                    directions.append((gene_b, gene_a, exclusives_b))
                if not directions:
                    continue

                directions.sort(key=lambda d: gene_usage.get(d[0], 0))
                pos_gene, neg_gene, pos_exclusives = directions[0]
                correct_trait = rng.choice(pos_exclusives)

                pos_sig = gene_sig.get(pos_gene, set())
                pos_limbo = gene_limbo.get(pos_gene, set())
                neg_sig = gene_sig.get(neg_gene, set())
                neg_limbo = gene_limbo.get(neg_gene, set())

                # (1) significant for both
                shared = [t for t in pos_sig & neg_sig if t != correct_trait]

                # (2-3) significant for neither (exclude limbo for both)
                excluded_union = pos_sig | pos_limbo | neg_sig | neg_limbo | {correct_trait}
                neither = [t for t in all_traits if t not in excluded_union]

                if not shared or len(neither) < 2:
                    continue

                rng.shuffle(shared)
                rng.shuffle(neither)

                distractors = [shared[0], neither[0], neither[1]]

                options = [correct_trait] + distractors
                rng.shuffle(options)
                correct_index = options.index(correct_trait)
                correct_letter = chr(ord("A") + correct_index)

                display_genes = [gene_a, gene_b]
                rng.shuffle(display_genes)

                question_text = (
                    f"Which of the following traits is associated with exactly one of these genes (and not the other)? "
                    f"Genes: {display_genes[0]}, {display_genes[1]}. "
                    f"[A] {options[0]} [B] {options[1]} [C] {options[2]} [D] {options[3]}"
                )

                questions.append(
                    {
                        "question": question_text,
                        "answer": correct_letter,
                        "tool": "GWAS",
                        "gene_pair": f"{gene_a}|{gene_b}",
                        "associated_gene": pos_gene,
                        "correct_trait": correct_trait,
                    }
                )

                used_pairs.add(pair_key)
                gene_usage[pos_gene] = gene_usage.get(pos_gene, 0) + 1
                gene_usage[neg_gene] = gene_usage.get(neg_gene, 0) + 1

                if len(questions) >= count:
                    break

    if not questions:
        print("No questions generated for set2; check data and thresholds", file=sys.stderr)
        return

    rng.shuffle(questions)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "set2.csv"
    pd.DataFrame(questions).to_csv(out_path, index=False)
    if len(questions) < count:
        print(
            f"Only generated {len(questions)} set2 questions (requested {count}); "
            "insufficient eligible gene pairs",
            file=sys.stderr,
        )
    print(f"Wrote {len(questions)} set2 questions to {out_path}")


def parse_args():
    parser = argparse.ArgumentParser(description="Generate GWAS MCQ benchmark set")
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
        help="Directory to write the generated set.csv",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=5e-8,
        help="P-value threshold for genome-wide significance",
    )
    parser.add_argument(
        "--variant",
        choices=["set1", "set2"],
        default="set1",
        help="Which set to generate (set1: single gene; set2: one-of-two genes)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    gwas_file = args.gwas_file
    if gwas_file is None:
        gwas_file = Path(args.db_dir) / "gwas_catalogue_association.tsv"

    if not gwas_file.exists():
        print(f"GWAS file not found at {gwas_file}", file=sys.stderr)
        raise SystemExit(1)

    if args.variant == "set1":
        generate_gwas_set(
            count=args.count,
            gwas_file=gwas_file,
            out_dir=Path(args.out_dir),
            seed=args.seed,
            p_value_threshold=args.p_threshold,
        )
    else:
        generate_gwas_set2(
            count=args.count,
            gwas_file=gwas_file,
            out_dir=Path(args.out_dir),
            seed=args.seed,
            p_value_threshold=args.p_threshold,
        )
