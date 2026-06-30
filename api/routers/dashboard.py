from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from fastapi import APIRouter

from api.config import OPTICAL_DIR, RESULTS_DIR
from api.services.data_loader import (
    invalidate_cache,
    load_candidates,
    load_lunar,
    load_meta,
    load_timeline,
)

router = APIRouter()


def _nan_to_none(val):
    """Convert NaN/NaT to None for JSON serialization."""
    if isinstance(val, float) and np.isnan(val):
        return None
    if isinstance(val, (np.floating,)) and np.isnan(val):
        return None
    if val is pd.NA or val is pd.NaT:
        return None
    if isinstance(val, np.generic):
        return val.item()
    return val


def _df_to_records(df: pd.DataFrame) -> list[dict]:
    """Convert DataFrame to JSON-safe list of dicts, handling NaN/NaT."""
    records = df.to_dict(orient="records")
    for rec in records:
        for k, v in rec.items():
            rec[k] = _nan_to_none(v)
    return records


@router.get("/dashboard-data")
def get_dashboard_data():
    """Return all candidates, timeline, meta, and lunar data as one JSON payload."""
    cand = load_candidates()
    timeline = load_timeline()
    meta = load_meta()
    lunar_df = load_lunar()

    # Build target index (FITS file stats)
    target_index_rows = []
    if OPTICAL_DIR.exists():
        for d in sorted(OPTICAL_DIR.iterdir()):
            if not d.is_dir():
                continue
            fits_count = 0
            for f in d.rglob("*"):
                name = f.name.lower()
                if f.is_file() and name.endswith((".fits", ".fits.fz")):
                    fits_count += 1
            target_index_rows.append({
                "target": d.name,
                "num_fits": fits_count,
                "has_ps_csv": (d / "ps.csv").exists(),
            })

    # Summary stats
    summary = {
        "total_candidates": len(cand),
        "observed_targets": int(meta[meta["nwatch"] > 0]["target"].nunique()) if "target" in meta.columns and "nwatch" in meta.columns else 0,
        "total_observations": len(timeline),
        "alive_targets": int((cand.get("Priority", 0) > 2).sum()) if len(cand) else 0,
    }

    return {
        "summary": summary,
        "candidates": _df_to_records(cand) if not cand.empty else [],
        "timeline": _df_to_records(timeline) if not timeline.empty else [],
        "meta": _df_to_records(meta) if not meta.empty else [],
        "lunar": _df_to_records(lunar_df) if not lunar_df.empty else [],
        "target_index": target_index_rows,
    }


@router.post("/dashboard-data/refresh")
def refresh_dashboard_data():
    """Invalidate cached data so next request re-reads from CSVs."""
    invalidate_cache()
    return {"status": "ok"}
