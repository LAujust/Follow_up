from __future__ import annotations

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

from api.config import GENERATED_PLAN_DIR
from api.services.data_loader import load_candidates
from dashboard.services.plan_service import generate_sitian_plan_text, generate_tnot_plan_text, save_plan_outputs

router = APIRouter()


class PlanRequest(BaseModel):
    target_names: list[str] = Field(default_factory=list)
    count: int = Field(default=6, ge=1, le=100)
    interval: int = Field(default=300, ge=1)
    exptime: int = Field(default=300, ge=1)
    expcount: int = Field(default=6, ge=1, le=100)
    priority: int = Field(default=6, ge=1, le=10)
    save_to_disk: bool = True


@router.post("/plans")
def get_plans(req: PlanRequest):
    """Generate TNOT and Sitian observation plans."""
    cand = load_candidates()
    if cand.empty:
        raise HTTPException(status_code=400, detail="No candidates data available")

    # Filter by requested targets, or use all
    if req.target_names:
        filtered = cand[cand["target"].isin(req.target_names)]
        if filtered.empty:
            raise HTTPException(status_code=400, detail="No matching targets found")
        targets_df = filtered
    else:
        targets_df = cand

    tnot_text = generate_tnot_plan_text(targets_df, count=req.count, interval=req.interval)
    sitian_text = generate_sitian_plan_text(
        targets_df,
        exptime=req.exptime,
        expcount=req.expcount,
        priority=req.priority,
    )

    result = {
        "tnot": tnot_text,
        "sitian": sitian_text,
        "target_count": len(targets_df),
    }

    if req.save_to_disk:
        paths = save_plan_outputs(tnot_text, sitian_text, GENERATED_PLAN_DIR)
        result["saved_paths"] = {k: str(v) for k, v in paths.items()}

    return result
