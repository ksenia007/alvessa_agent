"""Helpers for organising demo outputs under timestamped folders."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Dict, Tuple

REPO_ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = REPO_ROOT / "out"
LATEST_FILE = OUT_DIR / "latest_run.txt"


def _ensure_out_dir() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)


def _unique_run_dir(base_name: str) -> Path:
    candidate = OUT_DIR / base_name
    counter = 1
    while candidate.exists():
        candidate = OUT_DIR / f"{base_name}_{counter:02d}"
        counter += 1
    return candidate


def create_run_directory(prefix: str | None = None) -> Tuple[Path, str]:
    """Create (and mark) a fresh run directory under ``out/``.

    Parameters
    ----------
    prefix:
        Optional suffix to append after the timestamp (e.g. ``"ui"``).

    Returns
    -------
    tuple
        The created directory path and the base timestamp string.
    """
    _ensure_out_dir()
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    base_name = timestamp if not prefix else f"{timestamp}_{prefix}"
    run_dir = _unique_run_dir(base_name)
    run_dir.mkdir(parents=True, exist_ok=True)
    mark_latest_run(run_dir)
    return run_dir, timestamp


def build_output_paths(run_dir: Path) -> Dict[str, Path]:
    """Return standardised output file locations for a run directory."""
    return {
        "log": run_dir / "demo.log",
        "json": run_dir / "demo.json",
        "txt": run_dir / "demo.txt",
    }


def mark_latest_run(run_dir: Path) -> None:
    """Record the most recent run directory for other processes to reuse."""
    _ensure_out_dir()
    LATEST_FILE.write_text(run_dir.name, encoding="utf-8")


def get_latest_run_directory() -> Path | None:
    """Return the most recently marked run directory if it still exists."""
    try:
        name = LATEST_FILE.read_text(encoding="utf-8").strip()
    except FileNotFoundError:
        return None
    if not name:
        return None
    path = OUT_DIR / name
    return path if path.exists() else None
