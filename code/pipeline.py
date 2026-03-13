from __future__ import annotations
import warnings
warnings.filterwarnings("ignore")
import argparse
import hashlib
import json
import logging
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time

from coadd import coadd
from photometry import Photometry


FOLLOWUP_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OPTICAL_ROOT = Path(os.environ.get("FOLLOWUP_OPTICAL_DIR", str(Path.home() / "optical_data"))).expanduser()
DEFAULT_CANDIDATES = FOLLOWUP_ROOT / "Candidates.csv"
DEFAULT_PIPELINE_NAME = "pipeline"
DEFAULT_CONFIG = Path(__file__).with_name("config_file.json")


def _parse_float(v: Any) -> float | None:
    if v is None:
        return None
    s = str(v).strip().replace("\u2013", "-").replace('"', "")
    if s in {"", "nan", "None", "--"}:
        return None
    try:
        return float(s)
    except Exception:
        return None


def _utc_now_str() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _as_bool(v: Any, default: bool = False) -> bool:
    if v is None:
        return default
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s in {"1", "true", "yes", "y", "on"}:
        return True
    if s in {"0", "false", "no", "n", "off"}:
        return False
    return default


def _hash_file_list(files: list[Path]) -> str:
    payload = "|".join(sorted(str(f) for f in files))
    return hashlib.sha1(payload.encode("utf-8")).hexdigest()


def _load_config(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _load_candidates(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "EP Name" in df.columns:
        df = df.rename(columns={"EP Name": "target"})
    req = {"target", "RA", "Dec"}
    if not req.issubset(df.columns):
        missing = sorted(req - set(df.columns))
        raise ValueError(f"Candidates missing required columns: {missing}")
    return df


def _target_coord(row: pd.Series) -> tuple[float, float]:
    ora = _parse_float(row.get("o_RA"))
    odec = _parse_float(row.get("o_Dec"))
    ra = _parse_float(row.get("RA"))
    dec = _parse_float(row.get("Dec"))

    if ora is not None and odec is not None:
        return ora, odec, True
    if ra is None or dec is None:
        raise ValueError(f"Invalid RA/Dec for target={row.get('target')}")
    return ra, dec, False


def _list_telescopes(target_dir: Path, pipeline_name: str) -> list[Path]:
    out: list[Path] = []
    for d in sorted(target_dir.iterdir()):
        if not d.is_dir():
            continue
        if d.name.startswith("."):
            continue
        if d.name == pipeline_name:
            continue
        out.append(d)
    return out


def _fits_meta(path: Path) -> tuple[Time, str]:
    date_obs = None
    band = "UNKNOWN"
    with fits.open(path) as hdul:
        for hdu in hdul:
            hdr = hdu.header
            date_obs = hdr.get("DATE-OBS") or hdr.get("DATE") or date_obs
            filt = hdr.get("FILTER")
            if filt not in (None, "", "NONE"):
                band = str(filt).strip()
            if date_obs:
                break

    if not date_obs:
        t = Time(path.stat().st_mtime / 86400.0 + 40587.0, format="mjd")
    else:
        t = Time(date_obs)

    if band == "UNKNOWN":
        low = path.name.lower()
        for b in ["u", "g", "r", "i", "z", "w", "y", "rp","ip", "zs"]:
            if f"_{b}_" in low:
                band = b
                break
    if band == "rp":
        band = "r"
    elif band == "ip":
        band = "i"
    elif band == 'gp':
        band = 'g'
    elif band == 'zs':
        band = 'z'
        
    for b in ["u", "g", "r", "i", "z", "w", "y"]:
        if f"{b}-" in band:
            band = b
            break
    return t, band


def _scan_raw_fits(telescope_dir: Path, pipeline_name: str) -> list[Path]:
    files = []
    for p in telescope_dir.rglob("*"):
        if not p.is_file():
            continue
        low = p.name.lower()
        if not (low.endswith(".fits") or low.endswith(".fits.fz")):
            continue
        if "ref" in low:
            continue
        if "cutout" in low or "cutout" in str(p):
            continue
        if f"/{pipeline_name}/" in p.as_posix():
            continue
        if low.startswith("._"):
            continue
        files.append(p)
    return sorted(files)


def _group_files_by_band_time(files: list[Path], window_hours: float = 12.0) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for f in files:
        try:
            t, band = _fits_meta(f)
        except Exception:
            continue
        records.append({"file": f, "mjd": float(t.mjd), "band": str(band)})

    if not records:
        return []

    grouped: list[dict[str, Any]] = []
    by_band: dict[str, list[dict[str, Any]]] = {}
    for r in records:
        by_band.setdefault(r["band"], []).append(r)

    for band, rows in by_band.items():
        rows = sorted(rows, key=lambda x: x["mjd"])
        cur = [rows[0]]
        for r in rows[1:]:
            gap_hours = (r["mjd"] - cur[-1]["mjd"]) * 24.0
            if gap_hours <= float(window_hours):
                cur.append(r)
            else:
                mjds = [x["mjd"] for x in cur]
                grouped.append(
                    {
                        "band": band,
                        "files": [x["file"] for x in cur],
                        "mean_mjd": float(np.mean(mjds)),
                    }
                )
                cur = [r]

        mjds = [x["mjd"] for x in cur]
        grouped.append(
            {
                "band": band,
                "files": [x["file"] for x in cur],
                "mean_mjd": float(np.mean(mjds)),
            }
        )

    return sorted(grouped, key=lambda g: g["mean_mjd"])


def _ensure_logger(log_file: Path) -> logging.Logger:
    logger = logging.getLogger(str(log_file))
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fh = logging.FileHandler(log_file, encoding="utf-8")
    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.propagate = False
    return logger


def _read_csv_or_empty(path: Path, cols: list[str]) -> pd.DataFrame:
    if path.exists():
        try:
            return pd.read_csv(path)
        except Exception:
            pass
    return pd.DataFrame(columns=cols)


def _coadd_target(
    target: str,
    target_dir: Path,
    pipeline_dir: Path,
    pipeline_name: str,
    tel_cfg: dict[str, Any],
    logger: logging.Logger,
) -> list[dict[str, Any]]:
    idx_path = pipeline_dir / "coadd_index.csv"
    idx_cols = [
        "target",
        "telescope",
        "group_day_mjd",
        "band",
        "mean_mjd",
        "input_hash",
        "n_input",
        "output_file",
        "status",
        "message",
        "created_utc",
    ]
    idx = _read_csv_or_empty(idx_path, idx_cols)
    out_rows: list[dict[str, Any]] = []

    for tel_dir in _list_telescopes(target_dir, pipeline_name):
        telescope = tel_dir.name
        files = _scan_raw_fits(tel_dir, pipeline_name)
        if not files:
            logger.info(f"[{target}/{telescope}] no raw fits found")
            continue

        coadd_params = tel_cfg.get(telescope, {})
        do_coadd = _as_bool(coadd_params.get("do_coadd", True), default=True)
        if not do_coadd:
            logger.info(f"[{target}/{telescope}] do_coadd=False -> skip coadd stage")
            continue
        window_hours = float(coadd_params.get("coadd_window_hours", 12.0))
        groups = _group_files_by_band_time(files, window_hours=window_hours)
        size = coadd_params.get("size", None)
        if size == 'None':
            size = None

        for g in groups:
            band = str(g["band"])
            flist = list(g["files"])
            mean_mjd = float(g["mean_mjd"])
            day_mjd = int(np.floor(mean_mjd))
            time_tag = Time(mean_mjd, format="mjd").isot.replace(":", "").replace("-", "")
            input_hash = _hash_file_list(flist)
            prefix = f"coadd_{telescope}_{band}_{time_tag}"
            tel_out = pipeline_dir / telescope
            tel_out.mkdir(parents=True, exist_ok=True)
            expected_out = tel_out / f"{prefix}.fits"

            prev = idx[
                (idx.get("telescope") == telescope)
                & (idx.get("band") == band)
                & (idx.get("input_hash") == input_hash)
            ]
            if len(prev) > 0:
                logger.info(f"[{target}/{telescope}] skip coadd already recorded: {prefix}")
                continue

            try:
                out_file = coadd([str(p) for p in flist], size=size, save_dir=str(tel_out), prefix=prefix)
                msg = "ok"
                status = "ok"
                logger.info(f"[{target}/{telescope}] coadd ok: {out_file}")
            except Exception as exc:
                out_file = str(expected_out)
                msg = str(exc)
                status = "failed"
                logger.exception(f"[{target}/{telescope}] coadd failed: {prefix}")

            row = {
                "target": target,
                "telescope": telescope,
                "group_day_mjd": day_mjd,
                "band": band,
                "mean_mjd": mean_mjd,
                "input_hash": input_hash,
                "n_input": len(flist),
                "output_file": out_file,
                "status": status,
                "message": msg,
                "created_utc": _utc_now_str(),
            }
            idx = pd.concat([idx, pd.DataFrame([row])], ignore_index=True)
            out_rows.append(row)

    idx.to_csv(idx_path, index=False)
    return out_rows


def _run_photometry_target(
    target: str,
    ra: float,
    dec: float,
    target_dir: Path,
    pipeline_dir: Path,
    pipeline_name: str,
    tel_cfg: dict[str, Any],
    o_pos: bool,
    logger: logging.Logger,
) -> pd.DataFrame:
    photo_path = pipeline_dir / "photometry.csv"
    photo_cols = [
        "target",
        "telescope",
        "coadd_file",
        "mean_mjd",
        "band",
        "status_psf",
        "status_ap",
        "magpsf",
        "magpsf_err",
        "magap",
        "magap_err",
        "upper_limit",
        "zp",
        "zp_std",
        "message",
        "created_utc",
    ]
    photo_df = _read_csv_or_empty(photo_path, photo_cols)

    cat_dir = None
    for candidate in [target_dir / "ps.csv", target_dir / "ls.csv"]:
        if candidate.exists():
            cat_dir = str(candidate)
            break
    coadd_idx_path = pipeline_dir / "coadd_index.csv"
    mean_mjd_map: dict[str, float] = {}
    if coadd_idx_path.exists():
        try:
            cidx = pd.read_csv(coadd_idx_path)
            if {"output_file", "mean_mjd"}.issubset(cidx.columns):
                for _, r in cidx.iterrows():
                    try:
                        mean_mjd_map[str(r["output_file"])] = float(r["mean_mjd"])
                    except Exception:
                        pass
        except Exception:
            pass

    raw_telescopes = {d.name: d for d in _list_telescopes(target_dir, pipeline_name)}
    pipeline_telescopes = {d.name: d for d in _list_telescopes(pipeline_dir, "__none__")} if pipeline_dir.exists() else {}
    all_telescopes = sorted(set(raw_telescopes.keys()) | set(pipeline_telescopes.keys()))

    for telescope in all_telescopes:
        cfg = dict(tel_cfg.get(telescope, {}))
        methods_cfg = cfg.get("methods", None)
        method_cfg = cfg.get("method", None)
        if isinstance(methods_cfg, list):
            methods = [str(m).lower() for m in methods_cfg if str(m).strip()]
        elif method_cfg is None:
            methods = ["psf", "aperture"]
        else:
            m = str(method_cfg).lower().strip()
            if m in {"both", "all"}:
                methods = ["psf", "aperture"]
            else:
                methods = [m]
        methods = [m for m in methods if m in {"psf", "aperture"}]
        if not methods:
            methods = ["psf", "aperture"]

        pipeline_tel_dir = pipeline_dir / telescope
        # pipeline_tel_cutout_dir = pipeline_tel_dir / "cutouts"
        # if pipeline_tel_cutout_dir.is_dir():
        #     pipeline_tel_dir = pipeline_tel_cutout_dir
        coadd_candidates = sorted(pipeline_tel_dir.glob("coadd_*.fits")) if pipeline_tel_dir.exists() else []

        # For telescopes with pre-coadded products, allow photometry without running coadd.
        if not coadd_candidates and not _as_bool(cfg.get("do_coadd", True), default=True):
            raw_tel_dir = raw_telescopes.get(telescope, None)
            # print(raw_tel_dir)
            patterns = cfg.get(
                "precoadd_patterns",
                ["*stack*.fits", "*stack*.fits.fz", "*coadd*.fits", "*coadd*.fits.fz"],
            )
            if raw_tel_dir is not None:
                seen: set[str] = set()
                for pat in patterns:
                    for p in sorted(raw_tel_dir.rglob(str(pat))):
                        if not p.is_file():
                            continue
                        k = str(p)
                        if k in seen:
                            continue
                        seen.add(k)
                        coadd_candidates.append(p)
                if coadd_candidates:
                    logger.info(
                        f"[{target}/{telescope}] using {len(coadd_candidates)} pre-coadded files from raw telescope dir"
                    )
                else:
                    # If no pre-coadded products are found, do direct photometry on raw non-ref frames.
                    raw_files = _scan_raw_fits(raw_tel_dir, pipeline_name)
                    coadd_candidates.extend(raw_files)
                    if coadd_candidates:
                        logger.info(
                            f"[{target}/{telescope}] do_coadd=False and no stack/coadd files found; "
                            f"using {len(coadd_candidates)} raw files for direct photometry"
                        )

        if not coadd_candidates:
            logger.info(f"[{target}/{telescope}] no photometry candidate files found")
            continue

        for coadd_file in coadd_candidates:
            "skip cutout files"
            if "cutout" in coadd_file.name.lower():
                logger.info(f"[{target}/{telescope}] skip cutout file: {coadd_file.name}")
                continue
            mean_mjd = mean_mjd_map.get(str(coadd_file), None)
            try:
                t, band = _fits_meta(coadd_file)
            except Exception:
                    band = 'r'
                    
            if not band in ['g','r','i','z','y']:
                band = 'r'
                    
            if mean_mjd is None:
                mean_mjd = float(t.mjd)
                
            forced = bool(cfg.get("Forced", cfg.get("forced", False)))
            fwhm = float(cfg.get("fwhm", 3.0))
            sigma = float(cfg.get("sigma", 5.0))
            match_radius = float(cfg.get("match_radius", 1.0))
            fit_shape = cfg.get("fit_shape", (5, 5))
            if isinstance(fit_shape, int):
                fit_shape = (fit_shape, fit_shape)
            elif isinstance(fit_shape, list):
                fit_shape = tuple(fit_shape)

            result = {
                    "target": target,
                    "telescope": telescope,
                    "coadd_file": str(coadd_file),
                    "mean_mjd": mean_mjd,
                    "band": str(band),
                    "status_psf": "failed",
                    "status_ap": "failed",
                    "magpsf": np.nan,
                    "magpsf_err": np.nan,
                    "magap": np.nan,
                    "magap_err": np.nan,
                    "upper_limit": np.nan,
                    "zp": np.nan,
                    "zp_std": np.nan,
                    "message": "",
                    "created_utc": _utc_now_str(),
                }

            for method in methods:
                # prev = photo_df[photo_df.get("coadd_file") == str(coadd_file)]
                # if len(prev) > 0:
                #     logger.info(
                #         f"[{target}/{telescope}] skip photometry already recorded: "
                #         f"{coadd_file.name} method={method}"
                #     )
                #     continue

                out_plot_dir = pipeline_dir / telescope / "photometry" / coadd_file.stem / method
                out_plot_dir.mkdir(parents=True, exist_ok=True)

                # result = {
                #     "target": target,
                #     "telescope": telescope,
                #     "coadd_file": str(coadd_file),
                #     "mean_mjd": mean_mjd,
                #     "band": str(band),
                #     "status": "failed",
                #     "magpsf": np.nan,
                #     "magpsf_err": np.nan,
                #     "magap": np.nan,
                #     "magap_err": np.nan,
                #     "upper_limit": np.nan,
                #     "zp": np.nan,
                #     "zp_std": np.nan,
                #     "message": "",
                #     "created_utc": _utc_now_str(),
                # }

                try:
                    p = Photometry(cat_dir=cat_dir, path=str(out_plot_dir))
                    p.read_image(str(coadd_file))

                    if method == "aperture":
                        out = p.aperture_photometry(
                            ra=ra,
                            dec=dec,
                            r_ap=cfg.get("r_ap", None) if o_pos else 20.0,
                            r_in=cfg.get("r_in", None) if o_pos else 40.0,
                            r_out=cfg.get("r_out", None) if o_pos else 60.0,
                            sigma=sigma,
                            match_radius=match_radius,
                            fwhm=fwhm,
                            fit_shape=fit_shape,
                            mag_col=band,
                            cat_dir=cat_dir,
                            path=str(out_plot_dir),
                            show=False,
                        )
                        if isinstance(out, dict) and "upper_limit" in out:
                            result["upper_limit"] = float(out["upper_limit"])
                        else:
                            odf = out.to_pandas()
                            if len(odf) > 0:
                                if "mag" in odf.columns:
                                    result["magap"] = float(odf.iloc[0]["mag"])
                                if "mag_err" in odf.columns:
                                    result["magap_err"] = float(odf.iloc[0]["mag_err"])
                            result["upper_limit"] = p.uplim
                        result["status_ap"] = "ok"
                    else:
                        out = p.psf_photometry(
                            ra=ra,
                            dec=dec,
                            fwhm=fwhm,
                            sigma=sigma,
                            fit_shape=fit_shape,
                            match_radius=match_radius,
                            mag_col=band,
                            forced=forced,
                            cat_dir=cat_dir,
                            path=str(out_plot_dir),
                            show=False,
                        )
                        if isinstance(out, dict) and "upper_limit" in out:
                            result["upper_limit"] = float(out["upper_limit"])
                        else:
                            odf = out.to_pandas()
                            if len(odf) > 0:
                                if "mag" in odf.columns:
                                    result["magpsf"] = float(odf.iloc[0]["mag"])
                                if "mag_err" in odf.columns:
                                    result["magpsf_err"] = float(odf.iloc[0]["mag_err"])
                            result["upper_limit"] = p.uplim
                        result["status_psf"] = "ok"

                    result["zp"] = float(p.zp) if p.zp is not None else np.nan
                    result["zp_std"] = float(p.zp_std) if p.zp_std is not None else np.nan
                    result["message"] = "ok"
                    logger.info(f"[{target}/{telescope}] photometry ok: {coadd_file.name} method={method}")
                except Exception as exc:
                    result["status"] = "failed"
                    result["message"] = str(exc)
                    logger.exception(f"[{target}/{telescope}] photometry failed: {coadd_file.name} method={method}")
            photo_df = pd.concat([photo_df, pd.DataFrame([result])], ignore_index=True)
    

    # split tables
    # print(photo_df)
    # df_psf = photo_df[photo_df['method'] == 'psf'].copy()
    # df_ap  = photo_df[photo_df['method'] == 'aperture'].copy()

    # # rename columns
    # df_psf = df_psf.rename(columns={
    #     'mag': 'magpsf',
    #     'mag_err': 'magpsf_err'
    # })

    # df_ap = df_ap.rename(columns={
    #     'mag': 'magap',
    #     'mag_err': 'magap_err'
    # })

    # columns used to match rows
    df_merge = photo_df.copy()
    # keys = ['target','telescope','coadd_file','mean_mjd']

    # # merge
    # df_merge = pd.merge(
    #     df_psf,
    #     df_ap[keys + ['magap','magap_err']],
    #     on=keys,
    #     how='outer'
    # )
    df_merge.to_csv(photo_path, index=False)
    return photo_df


def run_pipeline(
    optical_root: Path,
    candidates_csv: Path,
    config_file: Path,
    pipeline_name: str = DEFAULT_PIPELINE_NAME,
    targets: list[str] | None = None,
    redo: bool | None = None,
) -> None:
    cfg = _load_config(config_file)
    config_redo = _as_bool(cfg.get("redo", False), default=False) if isinstance(cfg, dict) else False
    redo_flag = config_redo if redo is None else bool(redo)
    cand = _load_candidates(candidates_csv)
    
    if redo:
        import subprocess
        subprocess.run(['rm','-rf','~/optical_data/**/pipeline'], shell=True)

    if targets:
        cand = cand[cand["target"].isin(targets)]

    for _, row in cand.iterrows():
        target = str(row["target"])
        target_dir = optical_root / target
        if not target_dir.exists():
            continue

        try:
            ra, dec, o_pos = _target_coord(row)
            print(f"Processing target={target} RA={ra} Dec={dec} (o_pos={o_pos})")
        except Exception:
            continue

        pipeline_dir = target_dir / pipeline_name
        if redo_flag and pipeline_dir.exists():
            shutil.rmtree(pipeline_dir, ignore_errors=True)
        log_dir = pipeline_dir / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / f"{_utc_now_str()}.log"
        logger = _ensure_logger(log_file)

        logger.info(f"start target={target} optical_root={optical_root}")
        logger.info(f"using config_file={config_file if config_file.exists() else 'MISSING->defaults'}")
        logger.info(f"redo={redo_flag}")

        _coadd_target(
            target=target,
            target_dir=target_dir,
            pipeline_dir=pipeline_dir,
            pipeline_name=pipeline_name,
            tel_cfg=cfg,
            logger=logger,
        )

        _run_photometry_target(
            target=target,
            ra=ra,
            dec=dec,
            target_dir=target_dir,
            pipeline_dir=pipeline_dir,
            pipeline_name=pipeline_name,
            tel_cfg=cfg,
            o_pos=o_pos,
            logger=logger,
        )
        logger.info(f"done target={target}")


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Coadd + photometry incremental pipeline")
    p.add_argument("--optical-root", default=str(DEFAULT_OPTICAL_ROOT), help="Root optical data path")
    p.add_argument("--candidates", default=str(DEFAULT_CANDIDATES), help="Candidates.csv path")
    p.add_argument("--config", default=str(DEFAULT_CONFIG), help="config_file.json path")
    p.add_argument("--pipeline-name", default=DEFAULT_PIPELINE_NAME, help="Output pipeline folder name")
    p.add_argument("--targets", nargs="*", default=None, help="Optional target names filter")
    p.add_argument("--redo", action="store_true", help="Force remove prior pipeline outputs and rerun all")
    return p


def main() -> None:
    args = build_arg_parser().parse_args()
    run_pipeline(
        optical_root=Path(args.optical_root).expanduser(),
        candidates_csv=Path(args.candidates).expanduser(),
        config_file=Path(args.config).expanduser(),
        pipeline_name=args.pipeline_name,
        targets=args.targets,
        redo=args.redo,
    )


if __name__ == "__main__":
    main()
