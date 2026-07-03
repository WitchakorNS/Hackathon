"""
FastAPI app — serves the pipeline output to index.html.

Endpoints:
    GET /                     -> serves index.html (the storefront)
    GET /api/health           -> liveness + gate/QC status
    GET /api/pipeline         -> full PipelineOutput (metabolites, modules, ...)
    GET /api/metabolites      -> just the metData map (what the frontend renders)
    GET /api/modules          -> module cards
    GET /api/pathways         -> pathway enrichment table
    POST /api/pipeline/rerun  -> recompute (e.g. after swapping the data source)

The result is cached after first run; POST /api/pipeline/rerun refreshes it.
Run:  uvicorn api.main:app --reload  (from the backend/ directory)
"""
from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse

from pipeline import run_pipeline, GateError

app = FastAPI(title="Metabolic Signature Discovery Platform — Backend", version="1.0.0")

# Frontend is a standalone file served from the repo root (one level up).
_FRONTEND = Path(__file__).resolve().parent.parent.parent / "index.html"

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],          # dev only; lock this down for production
    allow_methods=["*"],
    allow_headers=["*"],
)

_CACHE: dict | None = None


def _get_output() -> dict:
    global _CACHE
    if _CACHE is None:
        _CACHE = run_pipeline().to_json()
    return _CACHE


@app.get("/", include_in_schema=False)
def index():
    if _FRONTEND.exists():
        return FileResponse(_FRONTEND)
    return JSONResponse({"detail": "index.html not found"}, status_code=404)


@app.get("/api/health")
def health():
    try:
        out = _get_output()
    except GateError as e:
        return JSONResponse({"status": "gate_failed", "detail": str(e)}, status_code=503)
    return {
        "status": "ok",
        "generated_at": out["generated_at"],
        "gate_passed": out["resolution"]["gate"]["passed"],
        "qc_passed": out["qc"]["passed"],
        "n_metabolites": len(out["metabolites"]),
    }


@app.get("/api/pipeline")
def pipeline():
    try:
        return _get_output()
    except GateError as e:
        return JSONResponse({"detail": str(e)}, status_code=503)


@app.get("/api/metabolites")
def metabolites():
    return _get_output()["metabolites"]


@app.get("/api/modules")
def modules():
    return _get_output()["modules"]


@app.get("/api/pathways")
def pathways():
    return _get_output()["pathways"]


@app.post("/api/pipeline/rerun")
def rerun():
    global _CACHE
    _CACHE = None
    try:
        return _get_output()
    except GateError as e:
        return JSONResponse({"detail": str(e)}, status_code=503)
