from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd


def _iter_targets(targets_df: pd.DataFrame) -> Iterable[tuple[str, float, float]]:
    for _, row in targets_df.iterrows():
        yield str(row["target"]), float(row["RA"]), float(row["Dec"])


def generate_tnot_plan_text(targets_df: pd.DataFrame, count: int = 6, interval: int = 300) -> str:
    lines = [
        "#dir EP2025",
        "#Filter rp",
        "#Binning 1",
        f"#Count {count}",
        f"#Interval {interval}",
        "",
    ]
    for name, ra_deg, dec_deg in _iter_targets(targets_df):
        c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="fk5")
        ra_hms = c.ra.to_string(unit=u.hour, sep=":", precision=2, pad=True)
        dec_dms = c.dec.to_string(unit=u.deg, sep=":", precision=2, pad=True, alwayssign=True)
        lines.append(f"{name}\t{ra_hms}\t{dec_dms}")
    return "\n".join(lines) + "\n"


def generate_sitian_plan_text(
    targets_df: pd.DataFrame,
    exptime: int = 300,
    expcount: int = 6,
    priority: int = 6,
) -> str:
    lines = ["objName\tra\tdec\texpTime\texpCount\ttype\tPriority"]
    for name, ra_deg, dec_deg in _iter_targets(targets_df):
        c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="fk5")
        ra_hms = c.ra.to_string(unit=u.hour, sep=":", precision=2, pad=True)
        dec_dms = c.dec.to_string(unit=u.deg, sep=":", precision=2, pad=True, alwayssign=True)
        lines.append(f"{name}\t{ra_hms}\t{dec_dms}\t{exptime}\t{expcount}\tOBJECT\t{priority}")
    return "\n".join(lines) + "\n"


def save_plan_outputs(tnot_text: str, sitian_text: str, out_dir: Path) -> dict[str, Path]:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d")
    dated = out_dir / stamp
    dated.mkdir(parents=True, exist_ok=True)

    tnot_path = dated / f"plan_{stamp}.txt"
    sitian_path = dated / f"sitian_plan_{stamp}.txt"

    tnot_path.write_text(tnot_text, encoding="utf-8")
    sitian_path.write_text(sitian_text, encoding="utf-8")

    return {"tnot": tnot_path, "sitian": sitian_path}
