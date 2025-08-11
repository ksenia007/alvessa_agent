from Bio import Entrez
import xml.etree.ElementTree as ET

Entrez.email = "svedula@princeton.edu"

def fetch_clinvar_entrez(rsid: str):
    # Step 1: find relevant ClinVar UIDs
    handle = Entrez.esearch(db="clinvar", term=f"{rsid}[VariantSet]", retmax=10)
    uids = Entrez.read(handle)["IdList"]
    handle.close()
    results = []
    print(uids)
    
    # Step 2: fetch each record’s XML
    for uid in uids:
        print(uid)
        h = Entrez.efetch(db="clinvar", id=uid, rettype="xml", retmode="text")
        print(h)
        xml = ET.fromstring(h.read()); h.close()
        print(h)
        # pull out assertion summary
        summary = xml.find(".//TraitSet")
        results.append({
            "uid": uid,
            "name": summary.get("trait_name") if summary is not None else None,
            "clinical_significance": xml.findtext(".//ClinicalSignificance/Description")
        })
    return results

if __name__ == "__main__":
    for rec in fetch_clinvar_entrez("rs121909229"):
        print(rec)
