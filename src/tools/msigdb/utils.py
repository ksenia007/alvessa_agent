import json

def create_gene_based_annotation_file():
    with open("../../../local_dbs/msigdb.v2025.1.Hs.json.txt", "r") as file:
        data = json.load(file)

    gene_dict = {}

    for gene_set_name, info in data.items():
        genes = info["geneSymbols"]
        new_info = {k: v for k, v in info.items() if k != "geneSymbols"}

        collection_name = info.get("collection", "Unknown")

        for gene in genes:
            if gene not in gene_dict:
                gene_dict[gene] = {}
            if collection_name not in gene_dict[gene]:
                gene_dict[gene][collection_name] = []
            
            annotation_entry = {"geneSetName": gene_set_name}
            annotation_entry.update(new_info)
            gene_dict[gene][collection_name].append(annotation_entry)

    with open("../../../local_dbs/gene_centered_msigdb.v2025.1.Hs.json.txt", "w") as file:
        json.dump(gene_dict, file, indent=4)


def extract_geneSetNames(gene_dict, gene, collection_substring):
    geneSetNames = set()
    if gene in gene_dict:
        for collection, annotations in gene_dict[gene].items():
            if collection_substring in collection:
                for annotation in annotations:
                    geneSetNames.add(annotation["geneSetName"])
    return list(geneSetNames)
