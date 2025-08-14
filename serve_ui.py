""" 
Author: Ksenia Sokolova <sokolova@princeton.edu>
Contributors: 
Created: 2025-08-12
Updated: 2025-08-12


Description: 

UI for the pipeline
"""
# serve_ui.py
# serve_ui.py
from __future__ import annotations

import io
import json
import sys
import threading
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import uvicorn
from fastapi import FastAPI, Query
from fastapi.encoders import jsonable_encoder
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles


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
STATE_FILE = BASE / "demo.json"
LOG_FILE = BASE / "demo.log"   # single, absolute path used everywhere

@app.get("/favicon.ico", include_in_schema=False)
def favicon():
    ico = BASE / "images" / "favicon.ico"
    if ico.exists():
        return FileResponse(ico, media_type="image/x-icon", headers={"Cache-Control":"no-store"})
    # fallback to your 32x32 PNG if you don't have an .ico yet
    return FileResponse(BASE / "images" / "rocket_favicon_32.png",
                        media_type="image/png",
                        headers={"Cache-Control":"no-store"})

@app.get("/apple-touch-icon.png", include_in_schema=False)
def apple_touch_icon():
    return FileResponse(BASE / "images" / "rocket_favicon_180.png",
                        media_type="image/png",
                        headers={"Cache-Control":"no-store"})

app.mount("/images", StaticFiles(directory=BASE / "images"), name="images")

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
    return FileResponse(BASE / "ui.html")

@app.get("/state")
def get_state():
    if not STATE_FILE.exists():
        return JSONResponse({"error": "No run yet. Click Run."}, status_code=404)
    return JSONResponse(json.loads(STATE_FILE.read_text()))

@app.get("/run")
def run(q: str = Query(..., description="User question")):
    # prevent overlapping runs that would fight over sys.stdout/sys.stderr
    if not RUN_LOCK.acquire(blocking=False):
        return JSONResponse({"error": "A run is already in progress."}, status_code=409)

    try:
        # start fresh log
        LOG_FILE.write_text("")

        # redirect stdio for the duration of the pipeline
        with redirect_stdio_to_log(LOG_FILE):
            print(f"=== START: {q} ===")
            state = run_pipeline(q)  # all prints inside stream to demo.log in real time
            print("=== END ===")

        # write sanitized state and return it
        state_json = to_jsonable(state)
        STATE_FILE.write_text(json.dumps(state_json, indent=2))
        return JSONResponse(state_json)

    finally:
        RUN_LOCK.release()

@app.get("/log")
def get_log(pos: int = 0):
    LOG_FILE.touch(exist_ok=True)
    with LOG_FILE.open("r", encoding="utf-8", errors="replace") as f:
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
