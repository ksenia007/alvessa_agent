# src/tools/prot/disprot.py
# DisProt data parser, query helpers, and CLI testing harness.

import sys
from pathlib import Path
import json
from typing import List, Dict, Any, Optional, Tuple

# --- Package/Repo Roots ---
PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# --- Imports after sys.path fix ---
from src.tools.prot import OUTPUT_DIR, DISPROT_JSON
from src.tools.prot.utils import log


def extract_disprot_regions(json_path: str) -> List[Dict[str, Any]]:
    """Extract disorder region annotations from a DisProt JSON release."""
    path = Path(json_path)
    if not path.exists():
        raise FileNotFoundError(f"DisProt JSON file not found: {json_path}")

    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)

    entries = raw.get("data", [])
    regions_out: List[Dict[str, Any]] = []

    for entry in entries:
        uniprot_id = entry.get("acc") or entry.get("uniprot_id")
        disprot_id = entry.get("disprot_id")
        prot_name = entry.get("name")
        organism = entry.get("organism")

        gene_symbol = None
        if entry.get("genes"):
            try:
                gene_symbol = entry["genes"][0]["name"]["value"]
            except Exception:
                pass

        for region in entry.get("regions", []):
            region_dict = {
                "uniprot_id": uniprot_id,
                "disprot_id": disprot_id,
                "gene_symbol": gene_symbol,
                "protein_name": prot_name,
                "organism": organism,
                "region_id": region.get("region_id"),
                "start_res": region.get("start"),
                "end_res": region.get("end"),
                "evidence_type": _infer_evidence_type(region),
                "term_name": region.get("term_name"),
                "term_namespace": region.get("term_namespace"),
                "function_tags": _extract_function_tags(region),
                "notes": region.get("term_comment") or None,
            }
            regions_out.append(region_dict)

    return regions_out


def fetch_disprot_for_uniprot(uniprot_id: str, json_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Fetch DisProt annotations for a given UniProt ID.

    Returns:
      {
        "consensus": {
            "residues": [{"residue_no": int, "score": 1.0, "tag": str}, ...],
            "regions": [{"start": int, "end": int, "term": str, "evidence": str}, ...],
            "metrics": {
                "n_regions": int,
                "n_function_tags": int
            }
        },
        "summary": str
      }
    """
    path = Path(json_path) if json_path else DISPROT_JSON
    if not path.exists():
        return {"consensus": {"residues": [], "regions": [], "metrics": {"n_regions": 0, "n_function_tags": 0}}, "summary": ""}

    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)

    entries = raw.get("data", [])
    entry = next((e for e in entries if e.get("acc") == uniprot_id), None)
    if not entry:
        return {"consensus": {"residues": [], "regions": [], "metrics": {"n_regions": 0, "n_function_tags": 0}}, "summary": ""}

    prot_name = entry.get("name")
    organism = entry.get("organism")
    gene_symbol = None
    if entry.get("genes"):
        try:
            gene_symbol = entry["genes"][0]["name"]["value"]
        except Exception:
            pass

    residues: List[Dict[str, Any]] = []
    regions: List[Dict[str, Any]] = []
    fn_tags: List[str] = []

    for region in entry.get("regions", []):
        start, end = region.get("start"), region.get("end")
        term = region.get("term_name")
        evidence = _infer_evidence_type(region)
        if start and end:
            regions.append({"start": int(start), "end": int(end), "term": term, "evidence": evidence})
            if term:
                fn_tags.append(term)
            for res in range(int(start), int(end) + 1):
                residues.append({"residue_no": res, "score": 1.0, "tag": term})

    metrics = {
        "n_regions": len(regions),
        "n_function_tags": len(set(fn_tags)),
    }

    summary_lines = [
        f"DisProt: {gene_symbol or prot_name} ({uniprot_id}, {organism})",
        f"  {len(regions)} annotated regions"
    ]
    for r in regions[:3]:
        summary_lines.append(f"  - {r['start']}-{r['end']}: {r['term']} ({r['evidence']})")
    if len(regions) > 3:
        summary_lines.append(f"  ... and {len(regions)-3} more regions")

    return {
        "consensus": {
            "residues": residues,
            "regions": regions,
            "metrics": metrics,
        },
        "summary": "\n".join(summary_lines),
    }


def merge_disprot_regions(regs: List[Dict]) -> List[Dict]:
    """Merge overlapping DisProt regions, preserving ambiguous evidence flags."""
    if not regs:
        return []
    regs = sorted(regs, key=lambda r: (int(r["start"]), int(r["end"])))
    merged: List[Dict[str, Any]] = []
    cur = dict(regs[0])
    cur_amb = 1 if cur.get("evidence") == "ambiguous" else 0
    for r in regs[1:]:
        if int(r["start"]) <= int(cur["end"]) + 1:
            cur["end"] = max(int(cur["end"]), int(r["end"]))
            cur_amb = cur_amb or (1 if r.get("evidence") == "ambiguous" else 0)
        else:
            merged.append({"start": int(cur["start"]), "end": int(cur["end"]), "ambiguous": bool(cur_amb)})
            cur = dict(r)
            cur_amb = 1 if r.get("evidence") == "ambiguous" else 0
    merged.append({"start": int(cur["start"]), "end": int(cur["end"]), "ambiguous": bool(cur_amb)})
    return merged


def residue_in_disprot(resno: int, regs: List[Dict]) -> Tuple[bool, bool]:
    """
    Check if a residue falls in any merged DisProt disorder region.
    Returns (in_region, is_ambiguous).
    """
    in_region, amb = False, False
    for r in regs:
        if int(r["start"]) <= int(resno) <= int(r["end"]):
            in_region = True
            amb = amb or r.get("ambiguous", False)
    return in_region, amb


def _infer_evidence_type(region: Dict[str, Any]) -> str:
    if region.get("unpublished"):
        return "unpublished"
    if "ambiguous" in (region.get("term_name") or "").lower():
        return "ambiguous"
    ec_name = region.get("ec_name", "")
    if "manual" in str(ec_name).lower():
        return "manual"
    return "other"


def _extract_function_tags(region: Dict[str, Any]) -> Optional[List[str]]:
    tags: List[str] = []
    if region.get("term_name"):
        tags.append(region["term_name"])
    if region.get("term_namespace"):
        tags.append(region["term_namespace"])
    return tags or None


# ----------------------------------------------------------------------
# CLI (testing)
# ----------------------------------------------------------------------
if __name__ == "__main__":
    disprot_file = Path(sys.argv[1]) if len(sys.argv) > 1 else DISPROT_JSON
    log(f"Loading DisProt file: {disprot_file}")
    regions = extract_disprot_regions(str(disprot_file))

    base_name = disprot_file.stem.replace(" ", "_")
    txt_out = OUTPUT_DIR / f"{base_name}_summary.txt"
    json_out = OUTPUT_DIR / f"{base_name}_regions.json"

    with txt_out.open("w", encoding="utf-8") as f:
        f.write(f"Total regions extracted: {len(regions)}\n\n")
        for r in regions[:10]:
            line = (
                f"{r.get('gene_symbol') or r.get('protein_name','')} "
                f"({r['uniprot_id']}, {r.get('organism') or ''}) "
                f"{r['region_id']} [{r['start_res']}-{r['end_res']}] "
                f"{r['evidence_type']} {r.get('term_name','')}\n"
            )
            f.write(line)

    with json_out.open("w", encoding="utf-8") as f:
        json.dump(regions, f, ensure_ascii=True, indent=2, sort_keys=True)

    if regions:
        test_uniprot = regions[0]["uniprot_id"]
        log(f"Testing fetch_disprot_for_uniprot on {test_uniprot}")
        out = fetch_disprot_for_uniprot(test_uniprot, disprot_file)
        summary_file = OUTPUT_DIR / f"{test_uniprot}_disprot_summary.txt"
        with summary_file.open("w", encoding="utf-8") as f:
            f.write(out["summary"])
        print(out["summary"])
        print(f"[INFO] Consensus residues: {len(out['consensus']['residues'])}")
        print(f"[INFO] Summary written to {summary_file.resolve()}")

    print("[OK] DisProt parsing complete.")
    print(f"  * TXT preview: {txt_out.resolve()}")
    print(f"  * JSON output: {json_out.resolve()}")
    print(f"[INFO] Regions extracted: {len(regions)}")
