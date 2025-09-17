from collections import defaultdict
from pathlib import Path
from typing import Dict


def parse_obo(obo_file: str | Path) -> Dict[str, str]:
    go_namespace: Dict[str, str] = {}
    with open(obo_file, "r", encoding="utf-8") as handle:
        current_id = None
        for line in handle:
            line = line.strip()
            if line.startswith("id: GO:"):
                current_id = line.split("id: ")[1]
            elif line.startswith("namespace:") and current_id:
                namespace = line.split("namespace: ")[1]
                go_namespace[current_id] = namespace
                current_id = None
    return go_namespace


def split_gaf_by_category(gaf_file: str | Path, go_namespace: Dict[str, str], output_file: str | Path) -> None:
    output_path = Path(output_file)
    out_files = {
        "biological_process": open(output_path.with_name(f"{output_path.name}_BP.gaf"), "w", encoding="utf-8"),
        "molecular_function": open(output_path.with_name(f"{output_path.name}_MF.gaf"), "w", encoding="utf-8"),
        "cellular_component": open(output_path.with_name(f"{output_path.name}_CC.gaf"), "w", encoding="utf-8"),
    }

    with open(gaf_file, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("!"):
                for out in out_files.values():
                    out.write(line)
                continue

            parts = line.strip().split("\t")
            if len(parts) > 4:
                go_id = parts[4]
                if go_id in go_namespace:
                    ns = go_namespace[go_id]
                    if ns in out_files:
                        out_files[ns].write(line)

    for out in out_files.values():
        out.close()


def convert_gaf_to_gmt(gaf_file: str | Path, output_file: str | Path) -> None:
    go_terms = defaultdict(set)

    with open(gaf_file, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 15:
                continue
            gene = cols[2]
            go_id = cols[4]
            go_terms[go_id].add(gene)

    with open(output_file, "w", encoding="utf-8") as out:
        for go_id, genes in go_terms.items():
            out.write(f"{go_id}\tNA\t" + "\t".join(sorted(genes)) + "\n")


def _process_dataset(base_dir: Path, gaf_name: str, output_stem: str) -> None:
    obo_file = base_dir / "go.obo"
    gaf_file = base_dir / gaf_name

    go_namespace = parse_obo(obo_file)
    output_file = base_dir / output_stem
    go_categories = ["BP", "MF", "CC"]

    split_gaf_by_category(gaf_file, go_namespace, output_file)
    print(
        "Files written to:",
        ", ".join(
            str(base_dir / f"{output_stem}_{category}.gaf") for category in go_categories
        ),
        flush=True,
    )

    for category in go_categories:
        input_gaf = base_dir / f"{output_stem}_{category}.gaf"
        output_gmt = base_dir / f"{output_stem}_{category}.gmt"
        convert_gaf_to_gmt(input_gaf, output_gmt)
        print(f"Converted {input_gaf} to {output_gmt}", flush=True)

    convert_gaf_to_gmt(gaf_file, base_dir / f"{output_stem}.gmt")


def main() -> None:
    repo_root = Path(__file__).resolve().parents[3]

    old_go_dir = repo_root / "local_dbs" / "old_go_data"
    _process_dataset(old_go_dir, "goa_human.gaf", "goa_human")

    pan_go_dir = repo_root / "local_dbs" / "pan_go_data"
    _process_dataset(pan_go_dir, "functionome_release.gaf", "functionome_release")


if __name__ == "__main__":
    main()
