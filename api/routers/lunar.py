from __future__ import annotations

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from api.services.data_loader import load_candidates
from dashboard.services.lunar_service import compute_lunar_curve

router = APIRouter()


class LunarCurveRequest(BaseModel):
    target: str = ""
    ra: float | None = None
    dec: float | None = None
    start_time: str | None = None
    ndays: int = Field(default=3, ge=1, le=30)
    step_hours: int = Field(default=6, ge=1, le=24)
    threshold: float = Field(default=30.0, ge=0, le=180)


@router.post("/lunar-curve")
def get_lunar_curve(req: LunarCurveRequest):
    """Compute lunar distance curve for a target or RA/Dec."""
    ra, dec = req.ra, req.dec

    # If target name given but no RA/Dec, look up from candidates
    if (ra is None or dec is None) and req.target:
        cand = load_candidates()
        match = cand[cand["target"] == req.target]
        if not match.empty:
            ra = float(match.iloc[0]["RA"])
            dec = float(match.iloc[0]["Dec"])

    if ra is None or dec is None:
        raise HTTPException(
            status_code=400,
            detail="RA and Dec are required (provide them directly or via a valid target name)",
        )

    curve = compute_lunar_curve(
        ra=ra,
        dec=dec,
        start_time=req.start_time,
        ndays=req.ndays,
        step_hours=req.step_hours,
        threshold=req.threshold,
    )

    return {
        "target": req.target,
        "ra": ra,
        "dec": dec,
        "curve": curve.to_dict(orient="records") if not curve.empty else [],
    }
