from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from reproject.mosaicking import find_optimal_celestial_wcs
from astropy.io import fits
from astropy.table import Table
import sys, os
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
from astropy.io.fits import Header
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.nddata import Cutout2D
import glob
from astropy.time import Time
from astropy.stats import sigma_clipped_stats


__all__ = ['coadd', 'coadd_lco']


def _get_valid_wcs_and_data(fname):
    """
    Iterate over HDUs and return the first (data, header, wcs)
    that has a valid celestial WCS.
    """
    with fits.open(fname) as hdul:
        for i, hdu in enumerate(hdul):
            if hdu.data is None:
                continue
            try:
                wcs = WCS(hdu.header)
                if wcs.has_celestial:
                    data = hdu.data.astype(float)
                    return data, hdu.header, wcs
            except Exception:
                continue
    raise RuntimeError(f"No valid celestial WCS found in {fname}")

"============================================="
"General coadd function"
"============================================="
def coadd(flist, size=None, save_dir='./', prefix='stacked'):
    """
    Coadd images using reproject.

    Parameters
    ----------
    flist : list
        List of FITS filenames to be stacked.
    size : int or tuple, optional
        Cutout size in pixels. If int, use (size, size).
    save_dir : str
        Directory to save stacked image.
    prefix : str
        Output filename prefix (without .fits).
    """

    if len(flist) == 0:
        raise ValueError("flist is empty")

    os.makedirs(save_dir, exist_ok=True)

    images = []
    weights = []

    ref_wcs = None
    ref_shape = None
    ref_header = None

    for i, fname in enumerate(flist):
        try:
            data, header, wcs = _get_valid_wcs_and_data(fname)
        except RuntimeError as e:
            print(f"[SKIP] {e}")
            continue

        # --- Cutout if requested ---
        if size is not None:
            if isinstance(size, int):
                size = (size, size)

            ny, nx = data.shape
            center_pix = (nx / 2, ny / 2)
            center_sky = wcs.pixel_to_world(*center_pix)

            cutout = Cutout2D(
                data,
                position=center_sky,
                size=size,
                wcs=wcs,
                mode='trim'
            )

            data = cutout.data
            wcs = cutout.wcs
            mean, median, std = sigma_clipped_stats(data,sigma=3.0)
            data = data - median

        # --- Use first valid image as reference ---
        if ref_wcs is None:
            ref_wcs = wcs
            ref_shape = data.shape
            ref_header = header.copy()
            ref_header.update(wcs.to_header())

        images.append((data, wcs))
        weights.append(np.ones_like(data))

        print(f"[OK] Loaded {fname}")

    if len(images) == 0:
        raise RuntimeError("No valid images to coadd")

    # --- Reproject & coadd ---
    print("Reprojecting and coadding images...")
    coadd_data, footprint = reproject_and_coadd(
        images,
        ref_wcs,
        ref_shape,
        combine_function='mean',
        reproject_function=reproject_interp,
    )

    # --- Save output ---
    outname = os.path.join(save_dir, f"{prefix}.fits")

    hdu = fits.PrimaryHDU(data=coadd_data, header=ref_header)
    hdu.header['NIMG'] = len(images)
    hdu.header['COMMENT'] = "Stacked with reproject_and_coadd"

    hdu.writeto(outname, overwrite=True)

    print(f"[DONE] Saved stacked image to {outname}")

    return outname


"============================================="
"Pre-defined function to coadd LCO images"
"============================================="
import os
import glob
import numpy as np
from astropy.time import Time
from astropy import units as u
from collections import defaultdict
from reproject import reproject_interp
from astropy.io import fits


def coadd_lco(root, save_dir, dt=1):
    """
    Coadd LCO images grouped by band and observation time.

    Parameters
    ----------
    root : str
        Directory containing *.fits.fz files
    save_dir : str
        Directory to save coadded images
    dt : float
        Time window in hours for grouping images

    Returns
    -------
    out_files : list
        List of output FITS file paths
    """

    os.makedirs(save_dir, exist_ok=True)

    # --------------------------------------------------
    # 1. collect files
    # --------------------------------------------------
    files = sorted(glob.glob(os.path.join(root, "*.fits.fz")))
    if not files:
        raise FileNotFoundError(f"No *.fits.fz files found in {root}")

    records = []

    # --------------------------------------------------
    # 2. read metadata (time, band, prefix)
    # --------------------------------------------------
    for f in files:
        try:
            _, header, _ = _get_valid_wcs_and_data(f)
        except Exception as e:
            print(f"[WARN] Skip {f}: {e}")
            continue

        date_obs = header.get("DATE-OBS")
        band = header.get("FILTER")

        if date_obs is None or band is None:
            print(f"[WARN] Skip {f}: missing DATE-OBS or FILTER")
            continue

        t = Time(date_obs, format="isot", scale="utc")
        prefix = os.path.basename(f).split("-")[0]

        records.append(
            dict(
                file=f,
                time=t,
                band=band,
                prefix=prefix,
            )
        )

    if not records:
        raise RuntimeError("No valid FITS files after header parsing.")

    # --------------------------------------------------
    # 3. group by band
    # --------------------------------------------------
    band_dict = defaultdict(list)
    for r in records:
        band_dict[r["band"]].append(r)

    out_files = []

    # --------------------------------------------------
    # 4. time-based grouping per band
    # --------------------------------------------------
    for band, band_records in band_dict.items():

        # sort by time
        band_records.sort(key=lambda x: x["time"].mjd)

        groups = []
        current_group = [band_records[0]]

        for r in band_records[1:]:
            dt_hour = (r["time"] - current_group[-1]["time"]).to(u.hour).value

            if dt_hour <= dt:
                current_group.append(r)
            else:
                groups.append(current_group)
                current_group = [r]

        groups.append(current_group)

        # --------------------------------------------------
        # 5. coadd each group
        # --------------------------------------------------
        for group in groups:

            if len(group) == 1:
                print(f"[INFO] Single image group, skip: {group[0]['file']}")
                continue

            prefix = group[0]["prefix"]

            times = Time([g["time"].isot for g in group])
            mean_time = Time(np.mean(times.mjd), format="mjd")

            # reference image
            ref_file = group[0]["file"]
            ref_data, ref_header, ref_wcs = _get_valid_wcs_and_data(ref_file)
            ref_data = ref_data.astype(float)

            stack = np.zeros_like(ref_data)
            weight = np.zeros_like(ref_data)

            for g in group:
                data, _, wcs = _get_valid_wcs_and_data(g["file"])
                data = data.astype(float)

                reproj, footprint = reproject_interp(
                    (data, wcs),
                    ref_wcs,
                    shape_out=ref_data.shape
                )

                m = np.isfinite(reproj)
                stack[m] += reproj[m]
                weight[m] += 1

            coadd = np.zeros_like(stack)
            valid = weight > 0
            coadd[valid] = stack[valid] / weight[valid]

            # --------------------------------------------------
            # 6. save
            # --------------------------------------------------
            out_name = (
                f"{prefix}_"
                f"{mean_time.isot.replace(':','').replace('-','')}_"
                f"{band}.fits"
            )
            out_path = os.path.join(save_dir, out_name)

            ref_header["DATE-OBS"] = mean_time.isot
            ref_header["FILTER"] = band
            ref_header["NCOMBINE"] = len(group)

            fits.writeto(
                out_path,
                coadd,
                header=ref_header,
                overwrite=True
            )

            out_files.append(out_path)
            print(f"[OK] {band}: {len(group)} images → {out_path}")

    return out_files