"""Compatibility wrapper for launching the packaged FastAPI UI."""

from __future__ import annotations

from src.alvessa.ui.server import APP as app
from src.alvessa.ui.server import serve


if __name__ == "__main__":  # pragma: no cover - manual invocation
    serve()
