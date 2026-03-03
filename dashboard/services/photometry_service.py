from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from config import OPTICAL_DIR
from photometry import Photometry


@dataclass
class PhotometryRequest:
    target: str
    fits_file: str
    ra: float
    dec: float
    method: str = "psf"
    fwhm: float = 3.0
    sigma: float = 5.0
    match_radius: float = 1.0
    forced: bool = True


def _target_output_dir(base_dir: Path, target: str) -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = base_dir / target / stamp
    out.mkdir(parents=True, exist_ok=True)
    return out


def run_photometry(req: PhotometryRequest, output_dir: Path) -> dict[str, Any]:
    out = _target_output_dir(output_dir, req.target)
    p = Photometry(path=str(out))
    p.read_image(req.fits_file)

    cat_csv = OPTICAL_DIR / req.target / "ps.csv"
    cat_dir = str(cat_csv) if cat_csv.exists() else None

    result: dict[str, Any] = {
        "status": "ok",
        "mag": None,
        "mag_err": None,
        "flux": None,
        "flux_err": None,
        "upper_limit": None,
        "zp": None,
        "zp_std": None,
        "output_dir": str(out),
        "diagnostic_paths": [],
        "raw_table_preview": None,
    }

    if req.method == "aperture":
        table = p.aperture_photometry(
            ra=req.ra,
            dec=req.dec,
            fwhm=req.fwhm,
            sigma=req.sigma,
            match_radius=req.match_radius,
            cat_dir=cat_dir,
            path=str(out),
            show=False,
        )
        result["raw_table_preview"] = table.to_pandas().head(20)
        try:
            ap_sum = float(table["aperture_sum_bkgsub"][0])
            ap_err = float(table["aperture_sum_err"][0])
            if ap_sum > 0 and p.zp is not None:
                result["mag"] = float(p.zp - 2.5 * __import__("numpy").log10(ap_sum))
                result["mag_err"] = float(
                    ((1.0857 * ap_err / ap_sum) ** 2 + (p.zp_std or 0.0) ** 2) ** 0.5
                )
                result["flux"] = ap_sum
                result["flux_err"] = ap_err
        except Exception:
            pass
    else:
        table_or_dict = p.psf_photometry(
            ra=req.ra,
            dec=req.dec,
            fwhm=req.fwhm,
            sigma=req.sigma,
            match_radius=req.match_radius,
            forced=req.forced,
            cat_dir=cat_dir,
            path=str(out),
            show=False,
        )
        if isinstance(table_or_dict, dict) and "upper_limit" in table_or_dict:
            result["upper_limit"] = float(table_or_dict["upper_limit"])
        else:
            df = table_or_dict.to_pandas().head(20)
            result["raw_table_preview"] = df
            if "mag" in df.columns:
                result["mag"] = float(df.iloc[0]["mag"])
            if "mag_err" in df.columns:
                result["mag_err"] = float(df.iloc[0]["mag_err"])
            if "flux" in df.columns:
                result["flux"] = float(df.iloc[0]["flux"])
            if "flux_err" in df.columns:
                result["flux_err"] = float(df.iloc[0]["flux_err"])

    result["zp"] = float(p.zp) if p.zp is not None else None
    result["zp_std"] = float(p.zp_std) if p.zp_std is not None else None

    for name in ("zeropoint_diagnostic.png", "psf_residual.png"):
        f = out / name
        if f.exists():
            result["diagnostic_paths"].append(str(f))

    if isinstance(result.get("raw_table_preview"), pd.DataFrame):
        pass

    return result
