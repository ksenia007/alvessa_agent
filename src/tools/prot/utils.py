# src/tools/prot/tool_prot_utils.py
# Author: Dmitri Kosenkov
# Created: 2025-08-25
# Updated: 2025-09-18
#
# Shared utilities for the agentic protein visualization tool.

import json
import sqlite3
import warnings
from typing import List, Dict, Optional, Tuple
from collections import OrderedDict
from pathlib import Path
from datetime import datetime

# --- Package/Repo Roots ---
PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]

# --- Imports ---
from src.config import DEBUG

# --- Local storage layout ---
from src.tools.prot import (
    REPO_ROOT,
    DB_PATH,
    PDB_DIR,
    STATIC_DIR,
    DB_MIN_MTIME
)

# ------------------------------
# Logging and small utilities
# ------------------------------
def log(msg: str) -> None:
    if DEBUG:
        print(f"[Protein] {msg} @ {datetime.now()}")


def fmt_float(x: Optional[float], ndigits: int = 2, default: str = "NA") -> str:
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
    if not DB_PATH.exists():
        raise RuntimeError(f"Database missing: {DB_PATH}")
    mtime = datetime.fromtimestamp(DB_PATH.stat().st_mtime)
    log(f"DB last modified: {mtime}")
    if mtime < DB_MIN_MTIME:
        raise RuntimeError(
            f"alvessa_proteins.db is outdated (last modified {mtime}), please update."
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
    ("plddt",
     "\n"
     "pLDDT (Predicted Local Distance Difference Test)\n"
     "Per-residue confidence score from AlphaFold (0-100):\n"
     "  - >90 : very high reliability (backbone + side chains accurate)\n"
     "  - 70-90 : backbone usually correct; side chains less certain\n"
     "  - <70 : lower confidence, often flexible or disordered regions\n"
     "Refs.: Mariani V. et al., 2013; Guo C. et al., 2022; EBI AlphaFold course\n"
    ),
    ("fpocket",
     "\n"
     "FPocket Druggability Score (mean_druggability_score)\n"
     "Numerical score (0-1) estimating drug-likeness of a pocket:\n"
     "  - 0.0 : very unlikely to bind\n"
     "  - ~0.5: borderline\n"
     "  - 1.0 : highly druggable\n"
     "Ref.: Schmidtke P. and Barril X., J. Med. Chem. 2010\n"
    ),
    ("sasa",
     "\n"
     "Total Solvent Accessible Surface Area (SASA)\n"
     "Computed per residue (Å²). High SASA implies solvent exposure; low SASA implies burial.\n"
     "Reported statistics include total SASA (sum across residues) and per-residue min, max, and mean values.\n"
     "Ref.: Mitternacht S. FreeSASA: An open source C library for solvent accessible surface area calculations  F1000Research 2016, 5:189 https://doi.org/10.12688/f1000research.7931.1\n"
    ),
    ("pi",
     "\n"
     "SASA Per-Residue Polarity Index (PI)\n"
     "Definition: PI = (P - A) / (P + A), where P = polar SASA and A = apolar SASA (Å²).\n"
     "Range: -1.0 <= PI <= +1.0\n"
     "  - +1.0 : completely polar (A = 0)\n"
     "  -  0.0 : equal polar and apolar contributions (P = A)\n"
     "  - -1.0 : completely apolar (P = 0)\n"
     "Meaning: Positive PI indicates predominantly polar surface; negative PI indicates predominantly hydrophobic surface.\n"
    ),
])


def interpretation_notes(include_fpocket: bool,
                         include_sasa: bool,
                         include_pi: bool,
                         include_plddt: bool = True) -> str:
    sections: List[str] = ["\nInterpretation Notes\n"]
    if include_plddt:
        sections.append(INTERPRETATION_TEXTS["plddt"])
    if include_fpocket:
        sections.append(INTERPRETATION_TEXTS["fpocket"])
    if include_sasa:
        sections.append(INTERPRETATION_TEXTS["sasa"])
    if include_pi:
        sections.append(INTERPRETATION_TEXTS["pi"])
    return "".join(sections)


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
        lines.append(
            "  pLDDT: min=" + fmt_float(plddt_stats.get("min"), 2) +
            ", max=" + fmt_float(plddt_stats.get("max"), 2) +
            ", avg=" + fmt_float(plddt_stats.get("avg"), 2)
        )
    if fpocket_stats:
        lines.append(
            "  FPocket (mean): min=" + fmt_float(fpocket_stats.get("min"), 3) +
            ", max=" + fmt_float(fpocket_stats.get("max"), 3) +
            ", avg=" + fmt_float(fpocket_stats.get("avg"), 3)
        )
    if sasa_stats:
        total_area = fmt_float(sasa_stats.get("total_area"), 1)
        lines.append(
            "  SASA total: " + total_area + " Å²"
            "(per-residue min=" + fmt_float(sasa_stats.get("min"), 2) +
            ", max=" + fmt_float(sasa_stats.get("max"), 2) +
            ", avg=" + fmt_float(sasa_stats.get("avg"), 2) + ")"
        )
    if pi_stats:
        lines.append(
            "  Polarity Index: min=" + fmt_float(pi_stats.get("min"), 2) +
            ", max=" + fmt_float(pi_stats.get("max"), 2) +
            ", avg=" + fmt_float(pi_stats.get("avg"), 2)
        )
    return "\n".join(lines)


# ------------------------------
# Flags and warnings (pLDDT only)
# ------------------------------
def derive_warning_flags(plddt_stats: Optional[Dict[str, float]]) -> List[str]:
    flags: List[str] = []
    try:
        if plddt_stats and plddt_stats.get("avg") is not None and plddt_stats["avg"] < 70.0:
            flags.append("Low mean pLDDT (<70).")
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
    data_script = "<script>\nconst protData = " + json_dumps_compact(data_obj) + ";\n</script>"
    html = (
        html_template_text.replace("{{GENE_SYMBOL}}", genes_display)
        .replace("{{CSS_INLINE}}", "<style>\n" + css_text + "\n</style>")
        .replace("{{JS_INLINE}}", "<script>\n" + js_text + "\n</script>")
        .replace("</body>", data_script + "\n</body>")
    )
    return html
