from __future__ import annotations

import json

from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse

from api.config import LIGHTCURVE_DIR, OPTICAL_DIR
from api.services.data_loader import get_target_photometry_csv

router = APIRouter()


@router.get("/lightcurve/{target}")
def get_lightcurve(target: str):
    """Return lightcurve data for a target as JSON.

    The frontend renders it client-side with Plotly, replacing the old
    approach of pre-generating self-contained Plotly HTML files.
    """
    # First try the pre-generated lightcurve HTML (extract embedded data)
    lc_html = LIGHTCURVE_DIR / f"{target}_lc.html"
    if lc_html.exists():
        # Return the raw data from photometry.csv if available
        photometry_csv = get_target_photometry_csv(target)
        if photometry_csv is not None:
            import pandas as pd
            from astropy.table import Table

            data = Table.read(str(photometry_csv))
            valid = data[data["message"] == "ok"]
            if len(valid) == 0:
                return {"target": target, "data": [], "status": "no_valid_data"}

            bands = list(set(valid["band"]))
            telescopes = list(set(valid["telescope"]))
            records = _table_to_records(valid)
            return {
                "target": target,
                "data": records,
                "bands": bands,
                "telescopes": telescopes,
                "status": "ok",
            }

    # Fallback: try the optical_data pipeline photometry.csv
    photometry_csv = get_target_photometry_csv(target)
    if photometry_csv is not None:
        from astropy.table import Table

        data = Table.read(str(photometry_csv))
        valid = data[data["message"] == "ok"]
        records = _table_to_records(valid) if len(valid) > 0 else []
        bands = list(set(valid["band"])) if len(valid) > 0 else []
        telescopes = list(set(valid["telescope"])) if len(valid) > 0 else []
        return {
            "target": target,
            "data": records,
            "bands": bands,
            "telescopes": telescopes,
            "status": "ok" if records else "no_valid_data",
        }

    raise HTTPException(status_code=404, detail=f"No lightcurve data found for target '{target}'")


@router.get("/lightcurve")
def list_lightcurve_targets():
    """List all targets with available lightcurve data."""
    targets = []

    # Check lcs directory
    if LIGHTCURVE_DIR.exists():
        for f in LIGHTCURVE_DIR.glob("*_lc.html"):
            target = f.name.replace("_lc.html", "")
            targets.append({"target": target, "source": "html"})

    # Also check pipeline directories
    if OPTICAL_DIR.exists():
        for d in sorted(OPTICAL_DIR.iterdir()):
            if d.is_dir() and d.name not in {t["target"] for t in targets}:
                csv_path = d / "pipeline" / "photometry.csv"
                if csv_path.exists():
                    targets.append({"target": d.name, "source": "pipeline"})

    return {"targets": sorted(targets, key=lambda t: t["target"])}


def _table_to_records(table) -> list[dict]:
    """Convert an astropy Table rows to JSON-safe dicts."""
    import numpy as np

    records = []
    for row in table:
        rec = {}
        for col in table.colnames:
            val = row[col]
            if isinstance(val, (np.floating, float)) and np.isnan(val):
                rec[col] = None
            elif isinstance(val, (np.integer, int)):
                rec[col] = int(val)
            elif isinstance(val, (np.floating, float)):
                rec[col] = float(val)
            elif isinstance(val, np.ma.core.MaskedConstant):
                rec[col] = None
            else:
                try:
                    rec[col] = str(val)
                except Exception:
                    rec[col] = None
        records.append(rec)
    return records
