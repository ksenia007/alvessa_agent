# src/tools/prot/utils.py
# Author: Dmitri Kosenkov
# Created: 2025-09-20
# Updated: 2025-12-08
#
# Shared utilities for the agentic protein visualization tool.

import json
import sqlite3
import warnings
from typing import List, Dict, Optional, Tuple, Any
from collections import OrderedDict
from pathlib import Path
from datetime import datetime

# --- Package/Repo Roots ---
PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]

# --- Imports ---
from src.config import DEBUG

# --- Local storage layout and constants ---
from src.tools.prot import (
    DB_PATH,
    PDB_DIR,
    REQUIRED_META,
    UPDATE_MSG,
    PLDDT_HIGH_CUTOFF,
    PLDDT_MEDIUM_CUTOFF,
    PLDDT_LOW_WARNING,
    BIO_RESIDUE_PREVIEW_MAX,
    BIO_EXCLUDED_IDS_MAX,
    FMT_FLOAT_NDIGITS,
    FMT_FLOAT_DEFAULT,
)

# ------------------------------
# Logging and small utilities
# ------------------------------
def log(msg: str) -> None:
    if DEBUG:
        print(f"[Protein] {msg} @ {datetime.now()}")


def fmt_float(x: Optional[float], ndigits: int = FMT_FLOAT_NDIGITS, default: str = FMT_FLOAT_DEFAULT) -> str:
    if x is None:
        return default
    try:
        return f"{float(x):.{ndigits}f}"
    except Exception:
        return default


def json_dumps_compact(obj) -> str:
    return json.dumps(obj, ensure_ascii=True, separators=(",", ":"), sort_keys=True)


def sanitize_residue_label(res3: Optional[str], residue_no: Optional[int]) -> Optional[str]:
    """
    Create labels like ALA1, ARG9 using 3-letter residue name and residue_no.
    Ensures uppercase letters and ASCII-only output.
    """
    if not res3 or residue_no is None:
        return None
    lab = f"{str(res3).strip().upper()}{int(residue_no)}"
    try:
        lab.encode("ascii", "strict")
    except Exception:
        lab = f"RES{int(residue_no)}"
    return lab


def uniq_preserve_order(seq: List[Any]) -> List[Any]:
    """Deduplicate a list while preserving order (hashable values)."""
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def uniq_dicts_preserve_order(seq: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Deduplicate list of dicts by value while preserving order."""
    seen = set()
    out = []
    for d in seq:
        try:
            key = tuple(sorted(d.items()))
        except Exception:
            key = None
        if key not in seen:
            seen.add(key)
            out.append(d)
    return out


def merge_residue_fields(*fields: Optional[str]) -> List[str]:
    """
    Merge residue fields from BioLiP2 (chain1, chain2, extra1, extra2).
    Example values: 'Q185 D187 W189 E198 N210'
    Returns flat list of residue labels like ['Q185','D187',...].
    """
    residues: List[str] = []
    for f in fields:
        if not f:
            continue
        residues.extend([tok.strip() for tok in f.split() if tok.strip()])
    return residues


# ------------------------------
# Database helpers
# ------------------------------
def get_connection() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH, check_same_thread=False)
    try:
        conn.execute("PRAGMA journal_mode=WAL;")
        conn.execute("PRAGMA synchronous=NORMAL;")
        conn.execute("PRAGMA temp_store=MEMORY;")
        conn.execute("PRAGMA cache_size=-131072;")
    except Exception as e:
        log(f"PRAGMA setup failed: {e}")
    return conn


def ensure_fresh_db() -> None:
    """
    Ensure that the local database exists and matches the required metadata.
    """
    if not DB_PATH.exists():
        raise RuntimeError(f"Database missing: {DB_PATH}. {UPDATE_MSG}")

    try:
        conn = sqlite3.connect(DB_PATH)
        cur = conn.cursor()
    except sqlite3.Error as e:
        raise RuntimeError(f"Could not open {DB_PATH}: {e}. {UPDATE_MSG}")

    cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='db_metadata';")
    if cur.fetchone() is None:
        print(f"try: {DB_PATH}. {UPDATE_MSG}")
        raise RuntimeError(f"Invalid database: missing db_metadata table. {UPDATE_MSG}")

    cur.execute("SELECT key, value FROM db_metadata;")
    rows = dict(cur.fetchall())

    for key, expected in REQUIRED_META.items():
        actual = rows.get(key)
        if actual != expected:
            problem = (
                f"missing required key '{key}'"
                if actual is None
                else f"{key} mismatch (expected {expected}, found {actual})"
            )
            raise RuntimeError(f"Invalid database: {problem}. {UPDATE_MSG}")

    conn.close()
    print(
        f"Database {DB_PATH} passed validation "
        f"(schema_version={rows['schema_version']}, build_time={rows['build_time']})."
    )


def load_pdb_inline(pdb_file: str) -> Optional[str]:
    pdb_path = PDB_DIR / pdb_file
    if not pdb_path.exists():
        warnings.warn(f"Missing PDB: {pdb_file}")
        return None
    return pdb_path.read_text(encoding="utf-8", errors="ignore")


# ------------------------------
# Stats and normalization
# ------------------------------
def compute_stats(values: List[float]) -> Optional[Dict[str, float]]:
    if not values:
        return None
    return {
        "min": min(values),
        "max": max(values),
        "avg": sum(values) / float(len(values)),
    }


def minmax_normalize(rows: List[Tuple[int, float]]) -> Tuple[Optional[Dict[str, float]], List[Dict[str, float]]]:
    if not rows:
        return None, []
    vals = [v for (_, v) in rows if v is not None]
    if not vals:
        return None, []
    stats = compute_stats(vals)
    vmin, vmax = stats["min"], stats["max"]
    rng = (vmax - vmin) if vmax > vmin else 1.0
    norm = [{"residue_no": int(rn), "score": (float(v) - vmin) / rng} for rn, v in rows if v is not None]
    return stats, norm


# ------------------------------
# Interpretation notes
# ------------------------------
INTERPRETATION_TEXTS = OrderedDict([
    ("plddt", f"""
AlphaFold pLDDT (Predicted Local Distance Difference Test)
Per-residue confidence score from AlphaFold (0-100):
  - >{PLDDT_HIGH_CUTOFF:.0f} : very high reliability
  - {PLDDT_MEDIUM_CUTOFF:.0f}-{PLDDT_HIGH_CUTOFF:.0f} : backbone usually correct
  - <{PLDDT_MEDIUM_CUTOFF:.0f} : lower confidence, often flexible
"""),
    ("fpocket", """
FPocket Druggability Score (0-1)
Estimate of drug-likeness of a pocket.
"""),
    ("sasa", """
Solvent Accessible Surface Area (SASA)
Computed per residue using FreeSASA (A^2).
"""),
    ("pi", """
Polarity Index (PI)
Derived from polar vs. apolar surface area decomposition.
"""),
    ("disorder", """
Structural Disorder Consensus Score
Integrates IUPred3 and DisProt.
"""),
    ("morf", """
MoRF Propensity (ANCHOR2)
Regions disordered in isolation but fold on binding.
"""),
    ("biolip2", """
BioLiP2 ligand-binding evidence from experimental PDB structures.
Curated protein-ligand interactions with cross-references to ChEMBL.
For small-molecule ligands, ChEMBL IDs are retrieved.
Common ions and solvents from PDB are excluded.
Summarizes unique ligands, binding sites, and supporting PDB entries.
"""),
    ("cysdb", """
CysDB cysteine chemoproteomics and functional context.
Flags are derived from curated chemoproteomics datasets and structural annotations:
  - detected: cysteine was detected in at least one chemoproteomics experiment.
  - ligandable: competition with a covalent probe indicates ligandable behavior.
  - hyperreactive: elevated intrinsic reactivity relative to baseline.
  - is_act_site: annotated as an active-site (catalytic) cysteine.
  - near_act_site: spatially close to an active-site residue within a distance threshold.
  - is_bind_site: lies in a ligand-binding site based on structural or UniProt annotations.
  - near_bind_site: spatially close to a ligand-binding site residue.
Neighbor lists indicate nearby active-site or binding-site residues (for example: G324, G326, S329).
For details on CysDB, see: Boatner, L. M. et al. Cell Chem. Biol. 2023, 30(6), 683-698.e3.
DOI: 10.1016/j.chembiol.2023.04.004.
"""),
])


def interpretation_notes(
    include_fpocket: bool,
    include_sasa: bool,
    include_pi: bool,
    include_plddt: bool = True,
    include_disorder: bool = True,
    include_morf: bool = True,
    include_biolip2: bool = False,
    include_cysdb: bool = False,
) -> str:
    sections: List[str] = ["\nInterpretation Notes\n"]
    if include_plddt:
        sections.append(INTERPRETATION_TEXTS["plddt"])
    if include_fpocket:
        sections.append(INTERPRETATION_TEXTS["fpocket"])
    if include_sasa:
        sections.append(INTERPRETATION_TEXTS["sasa"])
    if include_pi:
        sections.append(INTERPRETATION_TEXTS["pi"])
    if include_disorder:
        sections.append(INTERPRETATION_TEXTS["disorder"])
    if include_morf:
        sections.append(INTERPRETATION_TEXTS["morf"])
    if include_biolip2:
        sections.append(INTERPRETATION_TEXTS["biolip2"])
    if include_cysdb:
        sections.append(INTERPRETATION_TEXTS["cysdb"])
    return "".join(sections)


# ------------------------------
# Summary formatting helpers
# ------------------------------
def _format_block(title: str, stats: Dict[str, Any], indent: int = 2) -> str:
    if not stats:
        return ""
    pad = " " * indent
    lines = [f"{pad}{title}:"]
    for k, v in stats.items():
        lines.append(f"{pad*2}{k}={v}")
    return "\n".join(lines)


# ------------------------------
# Summaries (text only)
# ------------------------------
def make_summary_text(
    gene: str,
    uniprot_id: Optional[str],
    protein_id: Optional[str],
    pdb_file: Optional[str],
    plddt_stats: Optional[Dict[str, float]],
    fpocket_stats: Optional[Dict[str, float]],
    sasa_stats: Optional[Dict[str, float]],
    pi_stats: Optional[Dict[str, float]],
) -> str:
    lines = [f"Gene: {gene}"]
    if uniprot_id:
        lines.append(f"  UniProt ID: {uniprot_id}")
    if protein_id:
        lines.append(f"  Protein ID: {protein_id}")
    if plddt_stats:
        lines.append(_format_block("pLDDT", plddt_stats))
    if fpocket_stats:
        lines.append(_format_block("FPocket", fpocket_stats))
    if sasa_stats:
        lines.append(_format_block("SASA", sasa_stats))
    if pi_stats:
        lines.append(_format_block("Polarity Index", pi_stats))
    return "\n".join(lines)


def make_full_summary(
    gene_symbol: str,
    uniprot_id: str,
    protein_id: str,
    pdb_file: str,
    plddt_stats: Optional[Dict[str, float]],
    fpocket_stats: Optional[Dict[str, float]],
    sasa_stats: Optional[Dict[str, float]],
    pi_stats: Optional[Dict[str, float]],
    disorder_stats: Optional[Dict[str, float]],
    morf_stats: Optional[Dict[str, float]] = None,
    biolip2_summary: Optional[Dict[str, Any]] = None,
    cysdb_stats: Optional[Dict[str, Any]] = None,
) -> str:
    parts: List[str] = []
    parts.append(
        make_summary_text(
            gene_symbol,
            uniprot_id,
            protein_id,
            pdb_file,
            plddt_stats,
            fpocket_stats,
            sasa_stats,
            pi_stats,
        )
    )
    if disorder_stats:
        parts.append(_format_block("Disorder consensus", disorder_stats))
    if morf_stats:
        parts.append(_format_block("MoRF propensity", morf_stats))
    if cysdb_stats and cysdb_stats.get("include_cysdb"):
        parts.append("CysDB cysteine chemoproteomics summary:")
        parts.append(make_cysdb_summary(uniprot_id, cysdb_stats))
    if biolip2_summary:
        parts.append("BioLiP2 ligand-binding evidence:")
        parts.append(make_biolip2_summary(uniprot_id, biolip2_summary))
    return "\n".join(p for p in parts if p)


# ------------------------------
# Flags and warnings
# ------------------------------
def derive_warning_flags(plddt_stats: Optional[Dict[str, float]]) -> List[str]:
    flags: List[str] = []
    try:
        if plddt_stats and plddt_stats.get("avg") is not None and plddt_stats["avg"] < PLDDT_LOW_WARNING:
            flags.append(f"Low mean pLDDT (<{PLDDT_LOW_WARNING}).")
    except Exception:
        pass
    return flags


def append_flags_to_summary_text(summary_text: str, flags: List[str]) -> str:
    if not flags:
        return summary_text
    lines = [summary_text, "", "Flags:"]
    lines.extend([f"- {f}" for f in flags])
    return "\n".join(lines)


# ------------------------------
# HTML template injection helper
# ------------------------------
def inject_frontend_assets(
    html_template_text: str,
    css_text: str,
    js_text: str,
    data_obj: Dict,
    genes_display: str,
) -> str:
    """
    Inject inline CSS, JS, and JSON protein data into HTML template.
    Returns empty string if all provided genes have no usable feature data.
    """

    # --- Check for meaningful content in any gene ---
    has_data = False
    for gene, d in (data_obj or {}).items():
        if not isinstance(d, dict):
            continue

        # skip completely empty dicts
        if not d:
            continue

        # check arrays that carry actual per-residue data
        feature_keys = [
            "plddt",
            "fpocket",
            "sasa",
            "pi",
            "disorder",
            "morf",
            "biolip2",
            "cysdb_hyperreactive",
            "cysdb_ligandable",
            "cysdb_is_act_site",
            "cysdb_near_act_site",
            "cysdb_is_bind_site",
            "cysdb_near_bind_site",
        ]
        if any(d.get(k) and len(d[k]) > 0 for k in feature_keys):
            has_data = True
            break

        # if pdb text exists and is not empty, still count as data
        if d.get("pdb"):
            has_data = True
            break

    # --- No usable data: skip viewer entirely ---
    if not has_data:
        return ""  # returning empty HTML disables the Prot card automatically

    # --- Build full viewer HTML ---
    data_script = (
        "<script>\nconst protData = "
        + json_dumps_compact(data_obj)
        + ";\n</script>"
    )

    html = (
        html_template_text.replace("{{GENE_SYMBOL}}", genes_display)
        .replace("{{CSS_INLINE}}", "<style>\n" + css_text + "\n</style>")
        .replace(
            "{{JS_INLINE}}",
            data_script + "\n<script>\n" + js_text + "\n</script>",
        )
    )

    return html


# ------------------------------
# BioLiP2 summary text
# ------------------------------
def make_biolip2_summary(
    uniprot_id: str,
    summary: Dict[str, Any],
) -> str:
    """
    Format BioLiP2 ligand-binding evidence for text output.
    Matches style of make_full_summary and biolip2.py pipeline.
    """
    if not summary:
        return f"UniProt {uniprot_id}: No BioLiP2/ChEMBL evidence found."

    lines: List[str] = []
    lines.append(f"UniProt: {uniprot_id}")
    lines.append(f"  Structures ({summary['n_structures']}): {', '.join(summary['pdb_ids'])}")

    rs = summary.get("resolution") or {}
    if rs:
        lines.append(
            f"  Resolution (min/max/avg): "
            f"{fmt_float(rs.get('min'))}, {fmt_float(rs.get('max'))}, {fmt_float(rs.get('avg'))}"
        )

    lines.append(f"  Binding sites: {summary['n_binding_sites']}")

    if summary.get("go_terms"):
        lines.append(f"  GO terms: {', '.join(summary['go_terms'])}")
    if summary.get("pubmed_ids"):
        lines.append(f"  PubMed IDs: {', '.join(summary['pubmed_ids'])}")

    if summary.get("ligands"):
        lines.append("  Ligands (unique):")
        for lg in summary["ligands"]:
            lines.append(
                f"    {lg.get('chembl_id','')}  {lg.get('pref_name','')}  "
                f"{lg.get('inchikey','')}  comp:{lg.get('ligand_comp_id','')}"
            )

    if summary.get("sites"):
        lines.append("  Experimental Binding site details:")
        for bs_code, site in summary["sites"].items():
            lines.append(f"    {bs_code}:")
            if site.get("pdb_ids"):
                lines.append(f"      PDBs: {', '.join(site['pdb_ids'])}")
            if site.get("chains"):
                lines.append(f"      Chains: {', '.join(site['chains'])}")
            if site.get("residues"):
                nres = site.get("n_residues", len(site['residues']))
                res_preview = " ".join(site["residues"][:BIO_RESIDUE_PREVIEW_MAX])
                suffix = " ..." if nres > BIO_RESIDUE_PREVIEW_MAX else ""
                lines.append(f"      Residues ({nres}): {res_preview}{suffix}")
            if site.get("ligands"):
                lines.append("      Ligands:")
                for lg in site["ligands"]:
                    lines.append(
                        f"        {lg.get('chembl_id','')}  {lg.get('pref_name','')}  "
                        f"{lg.get('inchikey','')}  comp:{lg.get('ligand_comp_id','')}"
                    )

    excluded = summary.get("excluded", {})
    if excluded and excluded.get("count", 0) > 0:
        lines.append(
            f"  [INFO] {excluded['count']} ligands excluded "
            f"(ions/solvents or missing ChEMBL IDs)."
        )
        if excluded.get("ids"):
            ids_preview = ", ".join(excluded["ids"][:BIO_EXCLUDED_IDS_MAX])
            suffix = " ..." if len(excluded["ids"]) > BIO_EXCLUDED_IDS_MAX else ""
            lines.append(f"  Excluded ligand IDs: {ids_preview}{suffix}")

    return "\n".join(lines)


# ------------------------------
# CysDB summary text
# ------------------------------
def make_cysdb_summary(
    uniprot_id: str,
    stats: Dict[str, Any],
) -> str:
    """
    Format CysDB cysteine chemoproteomics and functional context for text output.

    Revised formatting per user specification:
      - Uses "CysDB Cysteines Entries:" header line
      - Removes example neighbor section
      - For each cysteine in near_* categories, append its respective neighbor list
    """
    if not stats or not stats.get("has_any_flag"):
        return f"UniProt {uniprot_id}: No CysDB evidence found."

    lines: List[str] = []
    lines.append(f"UniProt: {uniprot_id}")
    lines.append(
        "  CysDB Cysteines Entries: "
        f"detected={int(stats.get('n_identified', 0))}, "
        f"hyperreactive={int(stats.get('n_hyperreactive', 0))}, "
        f"ligandable={int(stats.get('n_ligandable', 0))}, "
        f"is_act_site={int(stats.get('n_is_act_site', 0))}, "
        f"near_act_site={int(stats.get('n_near_act_site', 0))}, "
        f"is_bind_site={int(stats.get('n_is_bind_site', 0))}, "
        f"near_bind_site={int(stats.get('n_near_bind_site', 0))}"
    )

    def _emit_block(title: str, key: str, include_neighbors: bool = False,
                    neighbor_map: Optional[Dict[str, str]] = None):
        """
        Emit:
            Title:
              P04406_C152
              P04406_C156 <neighbor list, if any>
        """
        residues = stats.get(key) or []
        if not residues:
            return
        lines.append(f"{title}:")
        for lab in residues:
            if include_neighbors and neighbor_map:
                neigh = neighbor_map.get(lab)
                if neigh:
                    lines.append(f"  {lab}  {neigh}")
                else:
                    lines.append(f"  {lab}")
            else:
                lines.append(f"  {lab}")

    # Neighbor maps: residue label → neighbor text
    near_act_map: Dict[str, str] = stats.get("near_act_site_neighbors_per_residue", {})
    near_bind_map: Dict[str, str] = stats.get("near_bind_site_neighbors_per_residue", {})

    # Emit blocks in revised order and style
    _emit_block("CysDB Detected", "detected_residues")
    _emit_block("Hyperreactive", "hyperreactive_residues")
    _emit_block("Ligandable", "ligandable_residues")
    _emit_block("Active-site cysteines", "is_act_site_residues")

    # Per-residue neighbor printing for near-active-site cysteines
    _emit_block(
        "Near active-site cysteines",
        "near_act_site_residues",
        include_neighbors=True,
        neighbor_map=near_act_map,
    )

    # Per-residue neighbor printing for near-binding-site cysteines
    _emit_block(
        "Near binding-site cysteines",
        "near_bind_site_residues",
        include_neighbors=True,
        neighbor_map=near_bind_map,
    )

    return "\n".join(lines)
