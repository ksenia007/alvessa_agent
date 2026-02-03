import os
import json
import sys

def create_gene_based_annotation_file(orig_file_name):
    orig_annotation_file = os.path.join('../../../local_dbs/', orig_file_name)
    with open(orig_annotation_file, "r") as file:
        data = json.load(file)

    gene_dict = {}

    for gene_set_name, info in data.items():
        genes = info.get("geneSymbols", [])
        collection_name = info.get("collection", "Unknown")

        for gene in genes:
            if gene not in gene_dict:
                gene_dict[gene] = {}
            gene_dict[gene][gene_set_name] = collection_name

    prefix = os.path.dirname(orig_annotation_file)

    with open(f"{prefix}/gene_centered_msigdb.json", "w") as file:
        json.dump(gene_dict, file, indent=4)

def main():
    if len(sys.argv) != 2:
        print("Usage: python process_msigdb.py name_of_msigdb_file")
        sys.exit(1)

    orig_annotation_file = sys.argv[1]
    create_gene_based_annotation_file(orig_annotation_file)

if __name__ == "__main__":
    main()