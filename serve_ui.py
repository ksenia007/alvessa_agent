""" 
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-08-12
Updated: 2025-08-12


Description: 

UI for the pipeline
"""
from __future__ import annotations

import io
import json
import sys
import os
import threading
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import uvicorn
from fastapi import FastAPI, Query
from fastapi.encoders import jsonable_encoder
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from src.config import DEBUG
from src.alvessa.workflow.output_paths import (
    build_output_paths,
    create_run_directory,
    get_latest_run_directory,
)
#------------------------------------------------
# For Windows Platform:

if sys.platform.startswith("win"):   # True only on Windows
    os.environ["HF_HUB_DISABLE_SYMLINKS"] = "1"
    from gliner import GLiNER

    m = GLiNER.from_pretrained("urchade/gliner_medium-v2.1")

    if DEBUG:
       print("Loaded:", type(m))
#------------------------------------------------

try:
    import pandas as pd
    HAS_PD = True
except Exception:
    HAS_PD = False

from run import run_pipeline  


# -----------------------------------------------------------------------------
# App & paths & images
# -----------------------------------------------------------------------------
app = FastAPI()
BASE = Path(__file__).parent.resolve()
WEB_DIR = BASE / "web"

if not WEB_DIR.exists():
    raise RuntimeError(f"Expected UI assets under {WEB_DIR}. Did you move the web directory?")


CURRENT_RUN_DIR = get_latest_run_directory()
CURRENT_OUTPUTS = build_output_paths(CURRENT_RUN_DIR) if CURRENT_RUN_DIR else None


def _set_current_run_dir(path: Path) -> None:
    global CURRENT_RUN_DIR, CURRENT_OUTPUTS
    CURRENT_RUN_DIR = path
    CURRENT_OUTPUTS = build_output_paths(path)

@app.get("/favicon.ico", include_in_schema=False)
def favicon():
    ico = WEB_DIR / "images" / "favicon.ico"
    if ico.exists():
        return FileResponse(ico, media_type="image/x-icon", headers={"Cache-Control":"no-store"})
    # fallback to your 32x32 PNG if you don't have an .ico yet
    return FileResponse(WEB_DIR / "images" / "rocket_favicon_32.png",
                        media_type="image/png",
                        headers={"Cache-Control":"no-store"})

@app.get("/apple-touch-icon.png", include_in_schema=False)
def apple_touch_icon():
    return FileResponse(WEB_DIR / "images" / "rocket_favicon_180.png",
                        media_type="image/png",
                        headers={"Cache-Control":"no-store"})

app.mount("/images", StaticFiles(directory=WEB_DIR / "images"), name="images")
app.mount("/static", StaticFiles(directory=WEB_DIR / "static"), name="static")

# -----------------------------------------------------------------------------
# JSON sanitizer for NumPy/Pandas types
# -----------------------------------------------------------------------------
def to_jsonable(obj):
    custom = {
        np.integer: int,
        np.floating: float,
        np.ndarray: lambda a: a.tolist(),
    }
    if HAS_PD:
        custom.update({
            pd.Series: lambda s: s.tolist(),
            pd.DataFrame: lambda df: df.to_dict(orient="records"),
        })
    return jsonable_encoder(obj, custom_encoder=custom)

# -----------------------------------------------------------------------------
# Tee stdout/stderr to file (line-buffered) while the pipeline runs
# -----------------------------------------------------------------------------
class Tee(io.TextIOBase):
    def __init__(self, *streams):
        self.streams = streams
        self._lock = threading.Lock()

    def write(self, data):
        with self._lock:
            for s in self.streams:
                try:
                    s.write(data)
                    s.flush()  # force line-by-line visibility
                except Exception:
                    pass
        return len(data)

    def flush(self):
        with self._lock:
            for s in self.streams:
                try:
                    s.flush()
                except Exception:
                    pass

@contextmanager
def redirect_stdio_to_log(log_path: Path):
    # open log file in text, line-buffered mode
    f = log_path.open("a", encoding="utf-8", buffering=1)

    # keep originals
    old_out, old_err = sys.stdout, sys.stderr

    # best-effort: enable line buffering on existing streams
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(line_buffering=True)  # py>=3.7
        except Exception:
            pass

    # tee both to console and file
    sys.stdout = Tee(old_out, f)
    sys.stderr = Tee(old_err, f)

    try:
        yield
    finally:
        sys.stdout, sys.stderr = old_out, old_err
        try:
            f.flush()
        except Exception:
            pass
        f.close()

# Only one run at a time (global stdio redirection)
RUN_LOCK = threading.Lock()

# -----------------------------------------------------------------------------
# Routes
# -----------------------------------------------------------------------------
@app.get("/", response_class=HTMLResponse)
def index():
    return FileResponse(WEB_DIR / "ui.html")

@app.get("/state")
def get_state():
    if not CURRENT_OUTPUTS:
        return JSONResponse({"error": "No run yet. Click Run."}, status_code=404)

    state_path = CURRENT_OUTPUTS["json"]
    if not state_path.exists():
        return JSONResponse({"error": "No run yet. Click Run."}, status_code=404)
    return JSONResponse(json.loads(state_path.read_text()))

@app.get("/run")
def run(q: str = Query(..., description="User question")):
    # prevent overlapping runs that would fight over sys.stdout/sys.stderr
    if not RUN_LOCK.acquire(blocking=False):
        return JSONResponse({"error": "A run is already in progress."}, status_code=409)

    try:
        run_dir, _ = create_run_directory("ui")
        _set_current_run_dir(run_dir)
        outputs = CURRENT_OUTPUTS
        log_path = outputs["log"]
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("")

        # redirect stdio for the duration of the pipeline
        with redirect_stdio_to_log(log_path):
            print(f"=== START: {q} ===")
            state = run_pipeline(q, output_dir=run_dir)  # all prints inside stream to demo.log in real time
            print("=== END ===")

        # write sanitized state and return it
        state_json = to_jsonable(state)
        outputs["json"].write_text(json.dumps(state_json, indent=2))

        # simple text summary akin to CLI run
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

@app.get("/log")
def get_log(pos: int = 0):
    if not CURRENT_OUTPUTS:
        return JSONResponse({"error": "No run yet. Click Run."}, status_code=404)

    log_path = CURRENT_OUTPUTS["log"]
    log_path.touch(exist_ok=True)
    with log_path.open("r", encoding="utf-8", errors="replace") as f:
        f.seek(pos)
        data = f.read()
        new_pos = f.tell()
    return JSONResponse(
        {"pos": new_pos, "data": data},
        headers={"Cache-Control": "no-store"},
    )

# -----------------------------------------------------------------------------
# Entrypoint
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    # Run with a single worker; multiple workers would have distinct processes/files.
    uvicorn.run(app, host="127.0.0.1", port=8000)
