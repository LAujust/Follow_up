from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = ROOT / "code"
LUNAR_DIR = ROOT / "lunar_distance"
OPTICAL_DIR = Path(os.environ.get("FOLLOWUP_OPTICAL_DIR", str(Path.home() / "optical_data"))).expanduser()
RESULTS_DIR = ROOT / "results"
WXT_DIR = ROOT / "wxtsource"

for p in (ROOT, CODE_DIR, LUNAR_DIR):
    sp = str(p)
    if sp not in sys.path:
        sys.path.insert(0, sp)

CANDIDATES_CSV = ROOT / "Candidates.csv"
TIMELINE_CSV = RESULTS_DIR / "all_targets_timeline.csv"
TIMELINE_HTML = RESULTS_DIR / "all_targets_timeline.html"
DISTRIBUTION_HTML = RESULTS_DIR / "observation_distribution.html"
META_CSV = RESULTS_DIR / "all_targets_meta.csv"
LUNAR_CSV = LUNAR_DIR / "Candidates_lunar.csv"

PHOTOMETRY_OUTPUT_DIR = RESULTS_DIR / "agent_photometry"
GENERATED_PLAN_DIR = RESULTS_DIR / "generated_plans"

FITS_SUFFIXES = (".fits", ".fits.fz")
