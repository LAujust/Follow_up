from __future__ import annotations

from pathlib import Path
import pandas as pd
import streamlit as st
from pandas.errors import ParserError

from config import (
    CANDIDATES_CSV,
    LUNAR_CSV,
    META_CSV,
    OPTICAL_DIR,
    ROOT,
    TIMELINE_CSV,
    FITS_SUFFIXES,
)


def _safe_read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except ParserError:
        # Fallback for occasional malformed lines in generated CSVs.
        return pd.read_csv(path, engine="python", on_bad_lines="skip")


@st.cache_data(show_spinner=False)
def load_candidates() -> pd.DataFrame:
    df = _safe_read_csv(CANDIDATES_CSV)
    if "EP Name" in df.columns:
        df = df.rename(columns={"EP Name": "target"})
    return df


@st.cache_data(show_spinner=False)
def load_timeline() -> pd.DataFrame:
    return _safe_read_csv(TIMELINE_CSV)


@st.cache_data(show_spinner=False)
def load_meta() -> pd.DataFrame:
    return _safe_read_csv(META_CSV)


@st.cache_data(show_spinner=False)
def load_lunar() -> pd.DataFrame:
    df = _safe_read_csv(LUNAR_CSV)
    if "EP Name" in df.columns:
        df = df.rename(columns={"EP Name": "target"})
    return df


@st.cache_data(show_spinner=False)
def build_target_index(root_optical_data: Path = OPTICAL_DIR) -> pd.DataFrame:
    rows = []
    if not root_optical_data.exists():
        return pd.DataFrame(columns=["target", "num_fits", "has_ps_csv"]) 

    for d in sorted(root_optical_data.iterdir()):
        if not d.is_dir():
            continue
        fits_count = 0
        for f in d.rglob("*"):
            name = f.name.lower()
            if f.is_file() and name.endswith(FITS_SUFFIXES):
                fits_count += 1
        rows.append(
            {
                "target": d.name,
                "num_fits": fits_count,
                "has_ps_csv": (d / "ps.csv").exists(),
            }
        )

    return pd.DataFrame(rows)


@st.cache_data(show_spinner=False)
def build_fits_index(root_optical_data: Path = OPTICAL_DIR) -> pd.DataFrame:
    rows = []
    if not root_optical_data.exists():
        return pd.DataFrame(columns=["target", "file", "path"])

    for target_dir in sorted(root_optical_data.iterdir()):
        if not target_dir.is_dir():
            continue
        for p in target_dir.rglob("*"):
            low = p.name.lower()
            if p.is_file() and low.endswith(FITS_SUFFIXES):
                rows.append(
                    {
                        "target": target_dir.name,
                        "file": p.name,
                        "path": str(p),
                    }
                )

    return pd.DataFrame(rows)


@st.cache_data(show_spinner=False)
def build_repo_file_index(root: Path = ROOT) -> pd.DataFrame:
    rows = []
    include_roots = [
        ("code", root / "code"),
        ("lunar_distance", root / "lunar_distance"),
        ("results", root / "results"),
        ("wxtsource", root / "wxtsource"),
        ("optical_data", OPTICAL_DIR),
    ]
    for rel, top in include_roots:
        if not top.exists():
            continue
        for p in top.rglob("*"):
            if not p.is_file():
                continue
            rows.append(
                {
                    "scope": rel,
                    "path": str(p),
                    "suffix": p.suffix.lower(),
                    "size_kb": round(p.stat().st_size / 1024, 2),
                }
            )

    return pd.DataFrame(rows)
