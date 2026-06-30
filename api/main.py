from __future__ import annotations

import os
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse

from api.routers import dashboard, lightcurves, lunar, plans

app = FastAPI(
    title="EP Follow-up API",
    description="Backend API for the EP follow-up observation dashboard",
    version="1.0.0",
)

# CORS: allow localhost for dev + optional CORS_ORIGIN for Render → Vercel
_cors_origins = ["http://localhost:3000", "http://localhost:5173"]
if cors_origin := os.environ.get("CORS_ORIGIN"):
    _cors_origins.append(cors_origin)

app.add_middleware(
    CORSMiddleware,
    allow_origins=_cors_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(dashboard.router, prefix="/api", tags=["dashboard"])
app.include_router(lunar.router, prefix="/api", tags=["lunar"])
app.include_router(plans.router, prefix="/api", tags=["plans"])
app.include_router(lightcurves.router, prefix="/api", tags=["lightcurves"])

RESULTS_DIR = Path(__file__).resolve().parents[1] / "results"


@app.get("/api/health")
def health():
    return {"status": "ok"}


@app.get("/api/results/{filename:path}")
def get_result_file(filename: str):
    """Serve static result files (PNGs, files) from the results directory."""
    file_path = RESULTS_DIR / filename
    if file_path.exists() and file_path.is_file():
        return FileResponse(str(file_path))
    return {"error": "File not found"}, 404
