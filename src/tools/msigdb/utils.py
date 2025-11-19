import json
import os

def create_gene_based_annotation_file(orig_annotation_file):
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

    orig_file_name = orig_annotation_file.split('/')[-1]
    prefix = os.path.dirname(orig_annotation_file)

    with open(f"{prefix}/gene_centered_{orig_file_name}", "w") as file:
        json.dump(gene_dict, file, indent=4)


def extract_geneSetNames(gene_dict, gene, collection_substring):
    geneSetNames = set()
    if gene in gene_dict:
        for geneSetName, collection in gene_dict[gene].items():
            if collection_substring == collection.split(":")[0]:
                geneSetNames.add(geneSetName)
    return list(geneSetNames)


