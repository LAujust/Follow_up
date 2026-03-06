from __future__ import annotations

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord, get_body
from astropy.time import Time
import sys


def compute_lunar_curve(
    ra: float,
    dec: float,
    start_time: str | None = None,
    ndays: int = 3,
    step_hours: int = 6,
    threshold: float = 30.0,
) -> pd.DataFrame:
    t0 = Time(start_time) if start_time else Time.now()
    times = t0 + np.arange(0, ndays * 24, step_hours) * u.hour
    target = SkyCoord(ra=ra, dec=dec, unit=u.deg)
    moon = get_body("moon", times)
    separation = moon.separation(target).deg
    if not isinstance(times, Time):
        times = Time(times)

        jd = times.jd

        synodic_month = 29.53058867

        # reference new moon (2000 Jan 6)
        jd_ref = 2451550.1

        age = (jd - jd_ref) % synodic_month
        phases = age / synodic_month
    

    return pd.DataFrame(
        {
            "time_iso": [t.isot for t in times],
            "separation_deg": separation,
            "moon_phase": phases,
            "above_threshold": separation >= threshold,
        }
    )
