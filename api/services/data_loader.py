from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd
from pandas.errors import ParserError

from api.config import (
    CANDIDATES_CSV,
    LUNAR_CSV,
    META_CSV,
    OPTICAL_DIR,
    TIMELINE_CSV,
)


def _safe_read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except ParserError:
        return pd.read_csv(path, engine="python", on_bad_lines="skip")


@lru_cache(maxsize=10)
def load_candidates() -> pd.DataFrame:
    df = _safe_read_csv(CANDIDATES_CSV)
    if "EP Name" in df.columns:
        df = df.rename(columns={"EP Name": "target"})
    return df


@lru_cache(maxsize=10)
def load_timeline() -> pd.DataFrame:
    return _safe_read_csv(TIMELINE_CSV)


@lru_cache(maxsize=10)
def load_meta() -> pd.DataFrame:
    return _safe_read_csv(META_CSV)


@lru_cache(maxsize=10)
def load_lunar() -> pd.DataFrame:
    df = _safe_read_csv(LUNAR_CSV)
    if "EP Name" in df.columns:
        df = df.rename(columns={"EP Name": "target"})
    return df


def get_target_photometry_csv(target: str) -> Path | None:
    """Find photometry.csv for a given target in the pipeline output."""
    candidates = [OPTICAL_DIR / target / "pipeline" / "photometry.csv"]
    for p in candidates:
        if p.exists():
            return p
    return None


def invalidate_cache():
    """Call after data files are updated to clear cached DataFrames."""
    load_candidates.cache_clear()
    load_timeline.cache_clear()
    load_meta.cache_clear()
    load_lunar.cache_clear()
