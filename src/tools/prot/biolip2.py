# src/tools/prot/biolip2.py
# Author: Dmitri Kosenkov
# Created: 2025-09-28
# Updated: 2025-10-01
#
# BioLiP2 + ChEMBL evidence manager.
# Exports ligand-binding evidence per UniProt ID.

import sys
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# --- Package/Repo Roots ---
PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.tools.prot.utils import (
    compute_stats,
    get_connection,
    merge_residue_fields,
    uniq_preserve_order,
    uniq_dicts_preserve_order,
    make_biolip2_summary,
)
from src.tools.prot import OUTPUT_DIR

# Common ions and solvents from PDB CCD (excluded)
EXCLUDED_IONS = {
    "NA", "K", "CL", "CA", "MG", "ZN", "FE", "MN", "CU", "CO", "CD", "NI", "HG",
    "IOD", "BR", "F", "SR", "BA", "AL",
    "SO4", "PO4", "GOL", "HOH", "ACT", "NO3", "CAC", "DMS",
}

# --- Index map for query results ---
IDX = {
    "pdb_id": 0, "chain_id": 1, "resolution": 2, "bs_code": 3,
    "ligand_comp_id": 4, "ligand_chain_id": 5, "ligand_count": 6,
    "residues_chain1": 7, "residues_chain2": 8,
    "residues_extra1": 9, "residues_extra2": 10,
    "ec_numbers": 11, "go_terms": 12,
    "annotation1": 13, "annotation2": 14, "annotation3": 15, "annotation4": 16,
    "uniprot_id": 17, "pubmed_id": 18,
    "protein_length": 19, "sequence": 20,
    "chembl_id": 21, "smiles": 22, "inchi": 23,
    "inchikey": 24, "pref_name": 25, "rn": 26,
}


# --- Helpers ------------------------------------------------------------
def _make_ligand_record(r) -> Dict[str, str]:
    """Build a clean ligand dictionary from a query row."""
    return {
        "pdb_id": r[IDX["pdb_id"]],
        "bs_code": r[IDX["bs_code"]],
        "chain_id": r[IDX["chain_id"]],
        "ligand_chain_id": r[IDX["ligand_chain_id"]],
        "ligand_comp_id": r[IDX["ligand_comp_id"]],
        "chembl_id": r[IDX["chembl_id"]],
        "pref_name": r[IDX["pref_name"]],
        "inchikey": r[IDX["inchikey"]],
    }


def _merge_site(site: Dict, r, ligand: Dict):
    """Add residues/ligands for a binding site."""
    residues = merge_residue_fields(
        r[IDX["residues_chain1"]],
        r[IDX["residues_chain2"]],
        r[IDX["residues_extra1"]],
        r[IDX["residues_extra2"]],
    )
    site["pdb_ids"].append(r[IDX["pdb_id"]])
    site["chains"].append(r[IDX["chain_id"]])
    site["residues"].extend(residues)
    site["ligands"].append(ligand)


# --- Main ---------------------------------------------------------------
def fetch_biolip2(conn, uniprot_id: str) -> Tuple[Optional[Dict[str, any]], List[Dict[str, str]]]:
    """
    Fetch BioLiP2 ligand-binding evidence for a UniProt ID.
    Returns:
      - biolip2_summary: dict with stats, sites, ligands, excluded info
      - biolip2_norm: list of ligand dicts (normalized) for JS frontend
    """
    
    excluded_ions_sql = ",".join([f"'{ion}'" for ion in EXCLUDED_IONS])

    query = f"""
        WITH filtered AS (
            SELECT *,
                   ROW_NUMBER() OVER (
                       PARTITION BY bs_code, ligand_chain_id, ligand_comp_id, chembl_id
                       ORDER BY resolution ASC
                   ) AS rn
            FROM biolip2_chembl
            WHERE uniprot_id = ?
              AND chembl_id IS NOT NULL
              AND TRIM(chembl_id) != ''
              AND UPPER(ligand_comp_id) NOT IN ({excluded_ions_sql})
        )
        SELECT * FROM filtered WHERE rn = 1
    """

    rows = conn.execute(query, (uniprot_id,)).fetchall()
    if not rows:
        return None, []

    pdb_ids, resolutions, go_terms, pubmed_ids = [], [], [], []
    ligands, sites = [], {}

    for r in rows:
        pdb_ids.append(r[IDX["pdb_id"]])
        if r[IDX["resolution"]] and r[IDX["resolution"]] > 0:
            resolutions.append(r[IDX["resolution"]])

        if r[IDX["go_terms"]]:
            go_terms.extend(r[IDX["go_terms"]].split(","))
        if r[IDX["pubmed_id"]]:
            pubmed_ids.extend(r[IDX["pubmed_id"]].split(","))

        ligand = _make_ligand_record(r)
        ligands.append(ligand)

        site = sites.setdefault(
            r[IDX["bs_code"]],
            {"pdb_ids": [], "chains": [], "residues": [], "ligands": []}
        )
        _merge_site(site, r, ligand)

    # Summary
    summary = {
        "n_structures": len(set(pdb_ids)),
        "pdb_ids": uniq_preserve_order(pdb_ids),
        "resolution": compute_stats(resolutions),
        "n_binding_sites": len(sites),
        "go_terms": uniq_preserve_order(go_terms),
        "pubmed_ids": uniq_preserve_order(pubmed_ids),
        "ligands": uniq_dicts_preserve_order(ligands),
        "sites": sites,
        "excluded": {"count": 0, "ids": []},  # reserved for future
    }

    # Normalize sites
    for site in sites.values():
        site["pdb_ids"] = uniq_preserve_order(site["pdb_ids"])
        site["chains"] = uniq_preserve_order(site["chains"])
        site["residues"] = uniq_preserve_order(site["residues"])
        site["ligands"] = uniq_dicts_preserve_order(site["ligands"])
        site["n_residues"] = len(site["residues"])

    # Normalized ligands for JS
    biolip2_norm = uniq_dicts_preserve_order(ligands)

    return summary, biolip2_norm


# ----------------------------------------------------------------------
# CLI (testing)
# ----------------------------------------------------------------------
if __name__ == "__main__":
    uniprot_id = sys.argv[1] if len(sys.argv) > 1 else "P04637"  # TP53 example
    print(f"[INFO] Using UniProt ID: {uniprot_id}")

    conn = get_connection()
    biolip2_summary, biolip2_norm = fetch_biolip2(conn, uniprot_id)

    base_name = f"{uniprot_id}_biolip2"
    txt_out = OUTPUT_DIR / f"{base_name}.txt"
    json_out = OUTPUT_DIR / f"{base_name}.json"

    with txt_out.open("w", encoding="utf-8") as f:
        if not biolip2_summary:
            f.write(f"No BioLiP2 evidence for {uniprot_id}\n")
        else:
            f.write(make_biolip2_summary(uniprot_id, biolip2_summary))

    with json_out.open("w", encoding="utf-8") as f:
        json.dump(
            {
                "summary": biolip2_summary,
                "biolip2": biolip2_norm,
            },
            f,
            ensure_ascii=True,
            indent=2,
            sort_keys=True,
        )

    print("[OK] BioLiP2 evidence exported.")
    if biolip2_summary:
        print(f"[INFO] Structures: {biolip2_summary['n_structures']}, "
              f"Sites: {biolip2_summary['n_binding_sites']}, "
              f"Ligands kept: {len(biolip2_summary['ligands'])}")
