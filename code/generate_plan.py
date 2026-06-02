#!/home/liangrd/anaconda3/envs/autophot/bin/python
"""
generate_plan.py — Generate observation plans for TNOT and/or Sitian telescopes.

Reads target info from Candidates.csv or wxt_candidates CSV, generates
observation plan files in the format required by each telescope, and
optionally sends them to a Feishu group chat.

Usage:
    # Generate plans for specific targets from Candidates.csv
    python generate_plan.py --targets EP260321a EP260110a --tnot --sitian

    # Generate plans from a wxt_candidates CSV file
    python generate_plan.py --from-wxt wxt_candidates_api_example.csv --tnot --sitian

    # TNOT only, with custom exposure parameters
    python generate_plan.py --targets EP260321a --tnot --count 6 --interval 300

    # Sitian only, with custom exposure parameters
    python generate_plan.py --targets EP260321a --sitian --expcount 10 --exptime 180

    # Send generated files to Feishu group chat
    python generate_plan.py --targets EP260321a --all --send

The --send flag will push the plan files to the default Feishu group chat.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

# Add parent to path so we can import utils
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import generate_tnot_plan, generate_tnot_object_json, generate_sitian_plan

# ── Constants ───────────────────────────────────────────────────────
FOLLOWUP_ROOT = Path("/home/liangrd/Follow_up")
CANDIDATES_CSV = FOLLOWUP_ROOT / "Candidates.csv"
WXTSOURCE_DIR = FOLLOWUP_ROOT / "wxtsource"
TNOT_SAVE_DIR = WXTSOURCE_DIR / "tnot"
SITIAN_SAVE_DIR = WXTSOURCE_DIR / "sitian"
FEISHU_GROUP_CHAT_ID = "oc_34de0318721ab375e02e2afd4ccd495b"


def resolve_targets_from_candidates(names: list[str]) -> tuple:
    """Resolve target names from Candidates.csv.

    Returns (target_names, ra_list, dec_list, obs_time_list, significance_list)
    """
    if not CANDIDATES_CSV.exists():
        raise FileNotFoundError(f"Candidates.csv not found: {CANDIDATES_CSV}")

    sources = Table.read(str(CANDIDATES_CSV))
    names_arr = np.array(sources["EP Name"])

    targets, ra_list, dec_list, obs_time_list, sig_list = [], [], [], [], []
    for name in names:
        idx = np.where(names_arr == name)[0]
        if len(idx) == 0:
            print(f"[WARN] Target '{name}' not found in Candidates.csv, skipping.")
            continue
        i = idx[0]
        targets.append(name)
        ra_list.append(float(sources["RA"][i]))
        dec_list.append(float(sources["Dec"][i]))
        # Obs Time — may not exist in all rows
        try:
            obs_time_list.append(str(sources["Obs Time"][i]))
        except (KeyError, IndexError):
            obs_time_list.append("")
        try:
            sig_list.append(float(sources["Sx"][i]))
        except (KeyError, IndexError):
            sig_list.append(0.0)

    return targets, ra_list, dec_list, obs_time_list, sig_list


def resolve_targets_from_wxt(csv_path: str) -> tuple:
    """Resolve targets from a wxt_candidates CSV file.

    Returns (target_names, ra_list, dec_list, obs_time_list, significance_list)
    """
    csv_path = str(WXTSOURCE_DIR / csv_path) if not os.path.isabs(csv_path) else csv_path
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"WXT CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)

    # The wxt CSV has columns: ep_name, ra, dec, obs_time, significance
    targets = df["ep_name"].tolist() if "ep_name" in df.columns else df["simbad_name"].tolist()
    ra_list = df["ra"].tolist()
    dec_list = df["dec"].tolist()
    obs_time_list = df["obs_time"].tolist() if "obs_time" in df.columns else [""] * len(targets)
    sig_list = df["significance"].tolist() if "significance" in df.columns else [0.0] * len(targets)

    return targets, ra_list, dec_list, obs_time_list, sig_list


def send_file_to_group(filepath: str) -> bool:
    """Send a file to the default Feishu group chat.

    Must be run from the directory containing the file (lark-cli requirement).
    """
    filepath = str(filepath)
    if not os.path.exists(filepath):
        print(f"[ERROR] File not found: {filepath}")
        return False

    fdir = os.path.dirname(filepath) or "."
    fname = os.path.basename(filepath)

    try:
        result = subprocess.run(
            ["lark-cli", "im", "+messages-send",
             "--as", "user",
             "--chat-id", FEISHU_GROUP_CHAT_ID,
             "--file", f"./{fname}"],
            cwd=fdir,
            capture_output=True, text=True, timeout=30,
        )
        if result.returncode != 0:
            print(f"[ERROR] Failed to send file: {result.stderr.strip()[:300]}")
            return False
        print(f"[OK] Sent: {fname}")
        return True
    except Exception as e:
        print(f"[ERROR] Send failed: {e}")
        return False


def generate_plans(
    targets: list[str],
    ra_list: list[float],
    dec_list: list[float],
    obs_time_list: list[str] | None = None,
    sig_list: list[float] | None = None,
    tnot: bool = False,
    sitian: bool = False,
    tnot_count: int = 6,
    tnot_interval: int = 300,
    sitian_expcount: int = 10,
    sitian_exptime: int = 180,
    sitian_priority: int = 6,
) -> list[str]:
    """Generate observation plans and return list of generated file paths."""
    generated_files = []
    obs_time_list = obs_time_list or [""] * len(targets)
    sig_list = sig_list or [0.0] * len(targets)

    if tnot:
        TNOT_SAVE_DIR.mkdir(parents=True, exist_ok=True)
        print(f"\n{'='*60}")
        print(f"Generating TNOT plans for {len(targets)} target(s)")
        print(f"{'='*60}")

        generate_tnot_plan(
            target=targets, ra=ra_list, dec=dec_list,
            count=tnot_count, interval=tnot_interval,
            save_path=str(TNOT_SAVE_DIR),
        )

        tnow = datetime.now(timezone(timedelta(hours=8)))
        plan_fname = f"plan_{tnow.year}{tnow.month:02d}{tnow.day:02d}.txt"
        plan_path = TNOT_SAVE_DIR / plan_fname
        if plan_path.exists():
            generated_files.append(str(plan_path))

        # Generate JSON info files
        generate_tnot_object_json(
            target=targets, ra=ra_list, dec=dec_list,
            obs_time=obs_time_list, significance=sig_list,
            save_path=str(TNOT_SAVE_DIR),
        )
        for t in targets:
            json_path = TNOT_SAVE_DIR / f"{t}_info.json"
            if json_path.exists():
                generated_files.append(str(json_path))

    if sitian:
        SITIAN_SAVE_DIR.mkdir(parents=True, exist_ok=True)
        print(f"\n{'='*60}")
        print(f"Generating Sitian plans for {len(targets)} target(s)")
        print(f"{'='*60}")

        generate_sitian_plan(
            target=targets, ra=ra_list, dec=dec_list,
            exptime=sitian_exptime, expcount=sitian_expcount,
            p=sitian_priority, save_path=str(SITIAN_SAVE_DIR),
        )

        tnow = datetime.now(timezone(timedelta(hours=8)))
        plan_fname = f"sitian_plan_{tnow.year}-{tnow.month:02d}-{tnow.day:02d}.txt"
        plan_path = SITIAN_SAVE_DIR / plan_fname
        if plan_path.exists():
            generated_files.append(str(plan_path))

    return generated_files


def format_summary(targets, generated_files, tnot, sitian) -> str:
    """Build a human-readable summary of generated plans."""
    lines = ["🔭 **观测计划已生成**", ""]
    lines.append(f"目标: {', '.join(targets)}")
    lines.append(f"台址: {'TNOT ' if tnot else ''}{'Sitian' if sitian else ''}")
    lines.append("")

    for f in generated_files:
        fpath = Path(f)
        lines.append(f"📄 `{fpath.name}`")

    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate observation plans for TNOT and/or Sitian telescopes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --targets EP260321a --all
  %(prog)s --targets EP260321a EP260110a --tnot --sitian --send
  %(prog)s --from-wxt wxt_candidates_api_example.csv --all --send
        """,
    )

    # Target selection (mutually exclusive group)
    target_group = parser.add_mutually_exclusive_group(required=True)
    target_group.add_argument(
        "--targets", nargs="+",
        help="Target EP names to look up in Candidates.csv"
    )
    target_group.add_argument(
        "--from-wxt",
        help="Read targets from a wxt_candidates CSV file in wxtsource/"
    )

    # Observatory selection
    parser.add_argument("--tnot", action="store_true", help="Generate TNOT plan")
    parser.add_argument("--sitian", action="store_true", help="Generate Sitian plan")
    parser.add_argument("--all", action="store_true", help="Generate both TNOT and Sitian plans")
    parser.add_argument("--send", action="store_true", help="Send generated files to Feishu group chat")

    # TNOT parameters
    parser.add_argument("--count", type=int, default=6, help="TNOT exposure count (default: 6)")
    parser.add_argument("--interval", type=int, default=300, help="TNOT exposure interval in sec (default: 300)")

    # Sitian parameters
    parser.add_argument("--expcount", type=int, default=10, help="Sitian exposure count (default: 10)")
    parser.add_argument("--exptime", type=int, default=180, help="Sitian exposure time in sec (default: 180)")
    parser.add_argument("--priority", type=int, default=6, help="Sitian priority (default: 6)")

    args = parser.parse_args()

    # Resolve observatories
    tnot = args.tnot or args.all
    sitian = args.sitian or args.all
    if not tnot and not sitian:
        print("Please specify --tnot, --sitian, or --all")
        sys.exit(1)

    # Resolve targets
    if args.targets:
        targets, ra_list, dec_list, obs_time_list, sig_list = resolve_targets_from_candidates(args.targets)
    else:
        targets, ra_list, dec_list, obs_time_list, sig_list = resolve_targets_from_wxt(args.from_wxt)

    if not targets:
        print("[ERROR] No targets resolved. Aborting.")
        sys.exit(1)

    print(f"Resolved {len(targets)} target(s): {', '.join(targets)}")
    for i, t in enumerate(targets):
        print(f"  {t}: RA={ra_list[i]:.4f}°, Dec={dec_list[i]:.4f}°")

    # Generate plans
    generated_files = generate_plans(
        targets=targets,
        ra_list=ra_list,
        dec_list=dec_list,
        obs_time_list=obs_time_list,
        sig_list=sig_list,
        tnot=tnot,
        sitian=sitian,
        tnot_count=args.count,
        tnot_interval=args.interval,
        sitian_expcount=args.expcount,
        sitian_exptime=args.exptime,
        sitian_priority=args.priority,
    )

    # Summary
    summary = format_summary(targets, generated_files, tnot, sitian)
    print(f"\n{summary}")

    # Send files
    if args.send and generated_files:
        print("\nSending files to group chat...")
        for f in generated_files:
            send_file_to_group(f)

    # Print JSON list of generated files for programmatic use
    print("\n--- GENERATED_FILES_JSON ---")
    print(json.dumps(generated_files))


if __name__ == "__main__":
    main()
