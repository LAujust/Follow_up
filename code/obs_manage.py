#!/home/liangrd/anaconda3/envs/autophot/bin/python
"""
obs_manage.py — EP Follow‑up Observation Manager

Synchronizes the Feishu cloud spreadsheet "Targets" → local Candidates.csv,
filters observable targets by lunar distance and observatory elevation,
and pushes daily reports to a Feishu group chat.

Usage:
    python obs_manage.py --sync          # sync cloud doc → CSV only
    python obs_manage.py --filter        # filter & push report
    python obs_manage.py --all           # sync + filter + push

Intended for cron: 30 07 * * * (UTC) = 15:30 Beijing time.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Optional

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body
from astropy.time import Time

# ── paths ────────────────────────────────────────────────────────────
FOLLOWUP_ROOT = Path("/home/liangrd/Follow_up")
CANDIDATES_CSV = FOLLOWUP_ROOT / "Candidates.csv"
LUNAR_CSV = FOLLOWUP_ROOT / "lunar_distance" / "Candidates_lunar.csv"
LOG_FILE = Path(__file__).with_suffix(".log")

# ── Feishu sheet token ────────────────────────────────────────────────
TARGETS_SHEET_TOKEN = "QX3NshTyVhQ4PsthkQOc7xepnId"
FEISHU_GROUP_CHAT_ID = "oc_34de0318721ab375e02e2afd4ccd495b"

# ── observatory sites ─────────────────────────────────────────────────
OBSERVATORIES = {
    "TNOT": EarthLocation(lat=43.471 * u.deg, lon=87.178 * u.deg, height=2080 * u.m),
    "Sitian": EarthLocation(lat=40.396 * u.deg, lon=117.576 * u.deg, height=900 * u.m),
    "WHUT": EarthLocation(lat=38.330 * u.deg, lon=93.896 * u.deg, height=4200 * u.m),
}

# ── filtering thresholds ─────────────────────────────────────────────
LUNAR_DISTANCE_THRESHOLD = 30.0   # degrees
ELEVATION_THRESHOLD = 50.0        # degrees

# prefer Lunar_Min_Distance over Lunar_Mean_Distance for filtering
LUNAR_COLUMN_PREFERENCE = ["Lunar_Min_Distance", "Lunar_Mean_Distance"]


def get_moon_phase(obs_time: Optional[Time] = None) -> float:
    """Return moon illumination fraction (0=new, 1=full).

    Uses astropy built-in moon position and solar elongation.
    """
    if obs_time is None:
        obs_time = Time.now()
    moon = get_body("moon", obs_time)
    sun = get_body("sun", obs_time)
    elongation = moon.separation(sun)
    # illumination fraction = (1 - cos(elongation)) / 2
    return float((1 - np.cos(elongation.rad)) / 2)

# ── column mapping (cloud sheet → CSV) ────────────────────────────────
COLUMNS = [
    "EP Name", "Priority", "Obs Time", "RA", "Dec", "r_err",
    "o_RA", "o_Dec", "Sx", "Redshift", "Classification",
    "GCNs", "GRB", "Fermi", "Swift", "LCO", "LCO active", "comments",
]

# ── logging ───────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler(sys.stderr),
    ],
)
log = logging.getLogger("obs_manage")


# ═══════════════════════════════════════════════════════════════════════
#  helpers
# ═══════════════════════════════════════════════════════════════════════

def _run_lark_cli(*args: str) -> dict:
    """Run a lark-cli command and return parsed JSON."""
    cmd = ["lark-cli", *args]
    log.debug("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    if result.returncode != 0:
        stderr = result.stderr.strip()
        # skip non-JSON OAuth / TUI messages printed to stderr
        if stderr and stderr.startswith("{"):
            pass
        # try parsing stdout first, then stderr
    for src in (result.stdout, result.stderr):
        try:
            return json.loads(src)
        except json.JSONDecodeError:
            continue
    raise RuntimeError(f"lark-cli failed: {result.stderr.strip()[:200]}")


def _excel_serial_to_datetime(serial: float) -> str:
    """Convert Excel serial date to ISO string.

    Excel epoch: 1899-12-30 for date serials that count from 1900-01-01.
    """
    epoch = datetime(1899, 12, 30)
    # serial may include fractional day for time
    days = int(serial)
    fraction = serial - days
    dt = epoch + timedelta(days=days, seconds=int(fraction * 86400))
    return dt.strftime("%Y-%m-%d %H:%M:%S")


def _extract_gcns_text(value) -> str:
    """Extract text from a GCNs cell (may be a list of URL objects)."""
    if value is None or value == "":
        return ""
    if isinstance(value, list):
        texts = []
        for item in value:
            if isinstance(item, dict):
                texts.append(item.get("text", ""))
            else:
                texts.append(str(item))
        return ", ".join(texts)
    return str(value)


# ═══════════════════════════════════════════════════════════════════════
#  sync: cloud sheet → Candidates.csv
# ═══════════════════════════════════════════════════════════════════════

def sync_cloud_to_csv() -> None:
    """Pull data from the Feishu 'Targets' sheet and write Candidates.csv."""
    log.info("Fetching cloud sheet …")
    resp = _run_lark_cli(
        "sheets", "+read",
        "--spreadsheet-token", TARGETS_SHEET_TOKEN,
    )

    values = resp["data"]["valueRange"]["values"]
    if not values:
        log.warning("Cloud sheet returned no rows")
        return

    header = values[0]
    rows = values[1:]

    records = []
    for row in rows:
        rec = {}
        for i, col in enumerate(COLUMNS):
            if i < len(row):
                val = row[i]
            else:
                val = None

            if col == "Obs Time" and isinstance(val, (int, float)):
                val = _excel_serial_to_datetime(float(val))
            elif col == "GCNs":
                val = _extract_gcns_text(val)
            elif val is None or val == "":
                val = ""
            else:
                val = str(val)

            rec[col] = val
        records.append(rec)

    df = pd.DataFrame(records, columns=COLUMNS)

    # ensure Priority is int
    df["Priority"] = pd.to_numeric(df["Priority"], errors="coerce").fillna(0).astype(int)

    df.to_csv(CANDIDATES_CSV, index=False)
    log.info("Synced %d rows to %s", len(df), CANDIDATES_CSV)


# ═══════════════════════════════════════════════════════════════════════
#  filter: elevation + lunar distance
# ═══════════════════════════════════════════════════════════════════════

def compute_elevations(
    ra_deg: float,
    dec_deg: float,
    obs_time: Time,
) -> dict[str, float]:
    """Compute target elevation angle at each observatory."""
    target = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
    elevations = {}
    for name, loc in OBSERVATORIES.items():
        altaz = AltAz(location=loc, obstime=obs_time)
        alt = target.transform_to(altaz).alt.deg
        elevations[name] = round(float(alt), 1)
    return elevations


def filter_observable_targets(obs_time: Optional[Time] = None) -> pd.DataFrame:
    """Filter candidates and return observable targets.

    Steps:
      1. Load Candidates.csv
      2. Load lunar distance data
      3. Merge
      4. Drop Priority == 0
      5. Drop lunar distance < threshold
      6. Drop targets below elevation threshold at ALL sites
    """
    if obs_time is None:
        # default to 15:30 Beijing time today
        tz_bj = timezone(timedelta(hours=8))
        now_bj = datetime.now(tz_bj)
        obs_dt = now_bj.replace(hour=15, minute=30, second=0, microsecond=0)
        obs_time = Time(obs_dt)

    log.info("Filtering for observation time: %s", obs_time.iso)

    # load candidates
    if not CANDIDATES_CSV.exists():
        log.error("Candidates.csv not found. Run --sync first.")
        return pd.DataFrame()

    cand = pd.read_csv(CANDIDATES_CSV)

    # load lunar
    lunar_cols = ["EP Name"]
    _lunar_dist_col = None  # actual column used for filtering
    if LUNAR_CSV.exists():
        lunar = pd.read_csv(LUNAR_CSV)
        # prefer Lunar_Min_Distance, fall back to Lunar_Mean_Distance
        for col in LUNAR_COLUMN_PREFERENCE:
            if col in lunar.columns:
                lunar_cols.append(col)
                _lunar_dist_col = col
                break
        # also load the other distance column for report display
        for col in ["Lunar_Min_Distance", "Lunar_Mean_Distance"]:
            if col in lunar.columns and col not in lunar_cols:
                lunar_cols.append(col)
        # also load Lunar_Phase if present
        if "Lunar_Phase" in lunar.columns:
            lunar_cols.append("Lunar_Phase")
        lunar = lunar[lunar_cols]
    else:
        log.warning("Lunar CSV not found at %s", LUNAR_CSV)
        lunar = pd.DataFrame(columns=["EP Name"])

    # merge
    merged = cand.merge(lunar, on="EP Name", how="left")

    # rename for clarity
    merged = merged.rename(columns={"EP Name": "EP_Name"})

    # drop Priority <= 1
    merged = merged[merged["Priority"] > 1].copy()
    if merged.empty:
        log.info("No targets with Priority > 0")
        return merged

    # drop by lunar distance (prefer min, fallback to mean)
    if _lunar_dist_col and _lunar_dist_col in merged.columns:
        lunar_mask = (
            merged[_lunar_dist_col].isna()
            | (merged[_lunar_dist_col] > LUNAR_DISTANCE_THRESHOLD)
        )
        n_dropped_lunar = (~lunar_mask).sum()
        if n_dropped_lunar:
            log.info(
                "Dropped %d targets (%s < %s°)",
                int(n_dropped_lunar), _lunar_dist_col, LUNAR_DISTANCE_THRESHOLD,
            )
        merged = merged[lunar_mask].copy()

    if merged.empty:
        log.info("No targets pass lunar distance filter")
        return merged

    # compute elevations
    elevation_rows = []
    for _, row in merged.iterrows():
        try:
            ra = float(row["RA"])
            dec = float(row["Dec"])
        except (ValueError, TypeError):
            continue
        elevs = compute_elevations(ra, dec, obs_time)
        elevation_rows.append(elevs)

    if elevation_rows:
        elev_df = pd.DataFrame(elevation_rows)
        merged = pd.concat([merged.reset_index(drop=True), elev_df], axis=1)

    # elevation filter: must be above threshold at AT LEAST ONE site
    elev_cols = list(OBSERVATORIES.keys())
    present_elev_cols = [c for c in elev_cols if c in merged.columns]
    if present_elev_cols:
        elev_mask = merged[present_elev_cols].max(axis=1) >= ELEVATION_THRESHOLD
        n_dropped_elev = (~elev_mask).sum()
        if n_dropped_elev:
            log.info(
                "Dropped %d targets (max elevation < %s° across all sites)",
                int(n_dropped_elev), ELEVATION_THRESHOLD,
            )
        merged = merged[elev_mask].copy()

    if merged.empty:
        log.info("No targets pass elevation filter")
        return merged

    # sort by priority (desc) then by max elevation (desc)
    if present_elev_cols:
        merged["_max_elev"] = merged[present_elev_cols].max(axis=1)
        merged = merged.sort_values(
            by=["Priority", "_max_elev"], ascending=[False, False]
        )
        merged = merged.drop(columns=["_max_elev"])
    else:
        merged = merged.sort_values("Priority", ascending=False)

    log.info("Final observable targets: %d", len(merged))
    return merged


# ═══════════════════════════════════════════════════════════════════════
#  push: format & send to group chat
# ═══════════════════════════════════════════════════════════════════════

def format_report(df: pd.DataFrame, obs_time: Time) -> str:
    """Build a Markdown report for the group chat."""
    tz_bj = timezone(timedelta(hours=8))
    dt_bj = obs_time.to_datetime(timezone=tz_bj)
    date_str = dt_bj.strftime("%Y-%m-%d")
    time_str = dt_bj.strftime("%H:%M")

    # compute moon phase
    try:
        phase_fraction = get_moon_phase(obs_time)
        phase_pct = round(phase_fraction * 100, 1)
        phase_str = f"🌙 月相: {phase_pct}%"
    except Exception:
        phase_str = ""

    filter_info = f"{phase_str}  |  月距阈值: {LUNAR_DISTANCE_THRESHOLD}°  |  高度角阈值: {ELEVATION_THRESHOLD}°"

    lines = [
        f"🔭 **Observable Targets Report**",
        f"",
        f"📅 {date_str} {time_str} (北京时间)",
        f"",
        f"{filter_info}",
        f"",
    ]

    if df.empty:
        lines.append(f"⚠️ 当前无可观测目标（lunar > {LUNAR_DISTANCE_THRESHOLD}° & elevation > {ELEVATION_THRESHOLD}°）")
        return "\n".join(lines)

    lines.append(f"共 **{len(df)}** 个可观测目标：")
    lines.append("")

    # table header with additional columns
    elev_cols = [c for c in OBSERVATORIES.keys() if c in df.columns]
    header_cols = ["EP Name", "Obs Time", "Sx", "Priority"] + elev_cols
    header = "| " + " | ".join(header_cols) + " |"
    sep = "|" + "|".join(["---"] * len(header_cols)) + "|"
    lines.append(header)
    lines.append(sep)

    for _, row in df.iterrows():
        # build EP Name with optional GCNs indicator
        ep_name = str(row["EP_Name"])
        gcns_val = row.get("GCNs", "")
        if gcns_val and str(gcns_val).strip():
            ep_name = f"{ep_name} 🔗"

        # Obs Time (short format)
        obs_time_val = row.get("Obs Time", "")
        obs_time_short = ""
        if obs_time_val and str(obs_time_val).strip():
            try:
                dt = pd.to_datetime(obs_time_val)
                obs_time_short = dt.strftime("%m-%d %H:%M")
            except Exception:
                obs_time_short = str(obs_time_val)[:10]

        # Sx (significance)
        sx_val = row.get("Sx", "")
        sx_str = str(sx_val) if sx_val and str(sx_val).strip() else "—"

        # Priority
        priority_str = str(int(row["Priority"]))

        cols = [ep_name, obs_time_short, sx_str, priority_str]

        for c in elev_cols:
            val = row.get(c)
            if pd.isna(val):
                cols.append("—")
            else:
                cols.append(f"{float(val):.0f}°")
        lines.append("| " + " | ".join(cols) + " |")

    lines.append("")
    lines.append("---")
    lines.append("🤖 Powered by Turing · obs_manage.py")

    return "\n".join(lines)


def push_report(text: str) -> None:
    """Send a message to the Feishu group chat via lark-cli."""
    # We use --as bot for reliable group-chat sending
    _run_lark_cli(
        "im", "+messages-send",
        "--as", "bot",
        "--chat-id", FEISHU_GROUP_CHAT_ID,
        "--text", text,
    )
    log.info("Report pushed to group chat %s", FEISHU_GROUP_CHAT_ID)


# ═══════════════════════════════════════════════════════════════════════
#  main
# ═══════════════════════════════════════════════════════════════════════

def main() -> None:
    parser = argparse.ArgumentParser(
        description="EP Follow‑up Observation Manager",
    )
    parser.add_argument(
        "--sync", action="store_true",
        help="Sync cloud sheet → Candidates.csv",
    )
    parser.add_argument(
        "--filter", action="store_true",
        help="Filter observable targets and push report",
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Full pipeline: sync → filter → push",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print report without sending",
    )
    parser.add_argument(
        "--no-push", action="store_true",
        help="Print report to stdout without pushing (for agent-relayed delivery)",
    )

    args = parser.parse_args()

    # default to --all
    if not any([args.sync, args.filter, args.all]):
        args.all = True

    try:
        if args.sync or args.all:
            sync_cloud_to_csv()

        if args.filter or args.all:
            obs_time = Time.now()
            # align to 15:30 Beijing time
            tz_bj = timezone(timedelta(hours=8))
            now_bj = datetime.now(tz_bj)
            target_bj = now_bj.replace(hour=15, minute=30, second=0, microsecond=0)
            obs_time = Time(target_bj)

            df = filter_observable_targets(obs_time)
            report = format_report(df, obs_time)

            if args.dry_run or args.no_push:
                print(report)
            else:
                push_report(report)
                print(report)

    except Exception:
        log.exception("Pipeline failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
