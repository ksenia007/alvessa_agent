from collections import defaultdict

def parse_obo(obo_file):
    go_namespace = {}
    with open(obo_file, "r") as f:
        current_id = None
        for line in f:
            line = line.strip()
            if line.startswith("id: GO:"):
                current_id = line.split("id: ")[1]
            elif line.startswith("namespace:") and current_id:
                namespace = line.split("namespace: ")[1]
                go_namespace[current_id] = namespace
                current_id = None
    return go_namespace

def split_gaf_by_category(gaf_file, go_namespace, output_file):
    out_files = {
        "biological_process": open(f"{output_file}_BP.gaf", "w"),
        "molecular_function": open(f"{output_file}_MF.gaf", "w"),
        "cellular_component": open(f"{output_file}_CC.gaf", "w"),
    }

    with open(gaf_file, "r") as f:
        for line in f:
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

def convert_gaf_to_gmt(gaf_file, output_file):

    go_terms = defaultdict(set)

    with open(gaf_file, 'r') as f:
        for line in f:
            if line.startswith('!'):  
                continue
            cols = line.strip().split('\t')
            if len(cols) < 15:
                continue
            gene = cols[2] 
            go_id = cols[4] 
            go_terms[go_id].add(gene)

    with open(output_file, 'w') as out:
        for go_id, genes in go_terms.items():
            out.write(f"{go_id}\tNA\t" + "\t".join(genes) + "\n")

### OLD GO FILES
obo_file = "../../local_dbs/old_go_data/go.obo"
gaf_file = "../../local_dbs/old_go_data/goa_human.gaf"

go_namespace = parse_obo(obo_file)
output_file = '../../local_dbs/old_go_data/goa_human'
go_categories = ['BP', 'MF', 'CC']

split_gaf_by_category(gaf_file, go_namespace, output_file)
print(f"Files written to: {output_file}_BP.gaf, {output_file}_MF.gaf, {output_file}_CC.gaf")

for category in go_categories:
    input_gaf = f'{output_file}_{category}.gaf'
    output_gmt = f'{output_file}_{category}.gmt'
    convert_gaf_to_gmt(input_gaf, output_gmt)
    print(f"Converted {input_gaf} to {output_gmt}")

convert_gaf_to_gmt(gaf_file, "../../local_dbs/old_go_data/goa_human.gmt")

### PAN-GO FILES
obo_file = "../../local_dbs/pan_go_data/go.obo"
gaf_file = "../../local_dbs/pan_go_data/functionome_release.gaf"

go_namespace = parse_obo(obo_file)
output_file = '../../local_dbs/pan_go_data/functionome_release'
go_categories = ['BP', 'MF', 'CC']

split_gaf_by_category(gaf_file, go_namespace, output_file)
print(f"Files written to: {output_file}_BP.gaf, {output_file}_MF.gaf, {output_file}_CC.gaf")

for category in go_categories:
    input_gaf = f'{output_file}_{category}.gaf'
    output_gmt = f'{output_file}_{category}.gmt'
    convert_gaf_to_gmt(input_gaf, output_gmt)
    print(f"Converted {input_gaf} to {output_gmt}")

convert_gaf_to_gmt(gaf_file, "../../local_dbs/pan_go_data/functionome_release.gmt")