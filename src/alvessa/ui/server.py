"""FastAPI server utilities for the Alvessa UI."""

from __future__ import annotations

import io
import json
import os
import sys
import threading
from contextlib import contextmanager
from pathlib import Path

import numpy as np
from fastapi import FastAPI, Query
from fastapi.encoders import jsonable_encoder
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from src.alvessa.pipeline import run_pipeline
from src.alvessa.workflow.output_paths import (
    build_output_paths,
    create_run_directory,
    get_latest_run_directory,
)
from src.config import DEBUG

try:  # Optional dependency (installed by default in this project)
    import pandas as pd
    HAS_PD = True
except Exception:  # pragma: no cover - defensive
    HAS_PD = False


REPO_ROOT = Path(__file__).resolve().parents[3]
WEB_DIR = REPO_ROOT / "web"

if not WEB_DIR.exists():  # pragma: no cover - early failure for misconfigured installs
    raise RuntimeError(
        f"Expected UI assets under {WEB_DIR}. Did you move the web directory?"
    )

APP = FastAPI()
APP.mount("/images", StaticFiles(directory=WEB_DIR / "images"), name="images")
APP.mount("/static", StaticFiles(directory=WEB_DIR / "static"), name="static")

CURRENT_RUN_DIR = get_latest_run_directory()
CURRENT_OUTPUTS = build_output_paths(CURRENT_RUN_DIR) if CURRENT_RUN_DIR else None


class Tee(io.TextIOBase):
    """Duplicate writes to multiple streams (used for log teeing)."""

    def __init__(self, *streams: io.TextIOBase):
        self.streams = streams
        self._lock = threading.Lock()

    def write(self, data: str) -> int:  # pragma: no cover - simple passthrough
        with self._lock:
            for stream in self.streams:
                try:
                    stream.write(data)
                    stream.flush()
                except Exception:
                    pass
        return len(data)

    def flush(self) -> None:  # pragma: no cover - simple passthrough
        with self._lock:
            for stream in self.streams:
                try:
                    stream.flush()
                except Exception:
                    pass


@contextmanager
def redirect_stdio_to_log(log_path: Path):
    """Context manager that tees stdout/stderr into a log file."""

    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = log_path.open("a", encoding="utf-8", buffering=1)
    old_out, old_err = sys.stdout, sys.stderr

    for stream in (sys.stdout, sys.stderr):  # pragma: no cover - best effort
        try:
            stream.reconfigure(line_buffering=True)
        except Exception:
            pass

    sys.stdout = Tee(old_out, fh)
    sys.stderr = Tee(old_err, fh)

    try:
        yield
    finally:
        sys.stdout, sys.stderr = old_out, old_err
        try:
            fh.flush()
        except Exception:
            pass
        fh.close()


RUN_LOCK = threading.Lock()


def _set_current_run_dir(path: Path) -> None:
    global CURRENT_RUN_DIR, CURRENT_OUTPUTS
    CURRENT_RUN_DIR = path
    CURRENT_OUTPUTS = build_output_paths(path)


def _to_jsonable(obj):
    custom = {
        np.integer: int,
        np.floating: float,
        np.ndarray: lambda a: a.tolist(),
    }
    if HAS_PD:
        custom.update(
            {
                pd.Series: lambda s: s.tolist(),
                pd.DataFrame: lambda df: df.to_dict(orient="records"),
            }
        )
    return jsonable_encoder(obj, custom_encoder=custom)


@APP.get("/", response_class=HTMLResponse)
def index():
    return FileResponse(WEB_DIR / "ui.html")


@APP.get("/favicon.ico", include_in_schema=False)
def favicon():
    ico = WEB_DIR / "images" / "favicon.ico"
    if ico.exists():
        return FileResponse(ico, media_type="image/x-icon", headers={"Cache-Control": "no-store"})
    return FileResponse(
        WEB_DIR / "images" / "rocket_favicon_32.png",
        media_type="image/png",
        headers={"Cache-Control": "no-store"},
    )


@APP.get("/apple-touch-icon.png", include_in_schema=False)
def apple_touch_icon():
    return FileResponse(
        WEB_DIR / "images" / "rocket_favicon_180.png",
        media_type="image/png",
        headers={"Cache-Control": "no-store"},
    )


@APP.get("/state")
def get_state():
    if not CURRENT_OUTPUTS:
        return JSONResponse({"error": "No run yet. Click Run."}, status_code=404)

    state_path = CURRENT_OUTPUTS["json"]
    if not state_path.exists():
        return JSONResponse({"error": "No run yet. Click Run."}, status_code=404)
    return JSONResponse(json.loads(state_path.read_text()))


@APP.get("/run")
def run(q: str = Query(..., description="User question")):
    if not RUN_LOCK.acquire(blocking=False):
        return JSONResponse({"error": "A run is already in progress."}, status_code=409)

    try:
        run_dir, _ = create_run_directory("ui")
        _set_current_run_dir(run_dir)
        outputs = CURRENT_OUTPUTS
        assert outputs is not None
        log_path = outputs["log"]
        log_path.write_text("", encoding="utf-8")

        with redirect_stdio_to_log(log_path):
            print(f"=== START: {q} ===")
            state = run_pipeline(q, output_dir=run_dir)
            print("=== END ===")

        state_json = _to_jsonable(state)
        outputs["json"].write_text(json.dumps(state_json, indent=2), encoding="utf-8")

        answer = state_json.get("llm_json", {}).get("answer", "")
        evidence = state_json.get("llm_json", {}).get("evidence", [])
        with outputs["txt"].open("w", encoding="utf-8") as txt:
            txt.write("=" * 80 + "\n")
            txt.write(f"Q: {q}\n\n")
            txt.write(f"Answer:\n{answer}\n\n")
            txt.write("Evidence:\n")
            for ev in evidence or []:
                txt.write(f"  - {ev}\n")

        return JSONResponse(state_json)

    finally:
        RUN_LOCK.release()


@APP.get("/log")
def get_log(pos: int = 0):
    if not CURRENT_OUTPUTS:
        return JSONResponse({"error": "No run yet. Click Run."}, status_code=404)

    log_path = CURRENT_OUTPUTS["log"]
    log_path.touch(exist_ok=True)
    with log_path.open("r", encoding="utf-8", errors="replace") as fh:
        fh.seek(pos)
        data = fh.read()
        new_pos = fh.tell()
    return JSONResponse({"pos": new_pos, "data": data}, headers={"Cache-Control": "no-store"})


if sys.platform.startswith("win"):  # pragma: no cover - platform-specific
    os.environ.setdefault("HF_HUB_DISABLE_SYMLINKS", "1")
    from gliner import GLiNER

    m = GLiNER.from_pretrained("urchade/gliner_medium-v2.1")
    if DEBUG:
        print("Loaded:", type(m))


def get_app() -> FastAPI:
    """Return the FastAPI application instance."""

    return APP


def serve(*, host: str = "127.0.0.1", port: int = 8000) -> None:
    """Run the FastAPI app with uvicorn on the requested address."""

    if port <= 0 or port >= 65536:
        raise ValueError("Port must be between 1 and 65535")

    print(f"[alvessa] Serving UI at http://{host}:{port}")

    import uvicorn

    uvicorn.run(APP, host=host, port=port, log_level="info")


__all__ = ["get_app", "serve", "APP"]
