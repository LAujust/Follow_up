import sys, os, subprocess
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_body
import json
from rich import print
import matplotlib.pyplot as plt
from pathlib import Path
from astropy.io import fits
import re
import shutil
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Circle
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval, ImageNormalize, LinearStretch, HistEqStretch
from psimage import *

__all__ = ['read_image','generate_tnot_plan','generate_tnot_object_json','generate_sitian_plan',
           'plot_lunar_distance','get_tnot_data','get_sitian_data',
           'check_source_dirs','fits_plot','show_shift','calculate_observation_stats',
           'show_obs_pies','show_cumulative_observations']


"""
=========================
Constants
=========================
"""

PIXEL_SCALE_ARCSEC = {
    "schmidt": 0.426,
    "LCO": 0.4,
    "GTC": 0.25,
    "TRT": 0.543,
    "Maidanak": 0.2668,
    "skob": 0.29,
    "VST": 0.213,
    "Keck_LRIS": 0.135,
    "REM": 0.57,
    "LT": 0.304,
    "CTIO": 0.526,
    "WFI": 0.237,
    "goodman": 8.0319879 * 0.036,
    "SAAO1m": 0.59,
    "TNOT": 0.781,
    "iTelescope": 0.359
}

TNOT_SERVER = 'tnot@119.78.162.172:/home/tnot/EP/plans'
TNOT_PASSWORD = 'nsqh.800@59726355'
FOLLOWUP_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OPTICAL_ROOT = Path(os.environ.get("FOLLOWUP_OPTICAL_DIR", str(Path.home() / "optical_data"))).expanduser()
DEFAULT_CANDIDATES_FILE = FOLLOWUP_ROOT / "Candidates.csv"
DEFAULT_RESULTS_DIR = FOLLOWUP_ROOT / "results"


"""
=========================
Functions
=========================
"""

def read_image(fdir):
    """Read FITS image and return data and WCS.

    Args:
        fdir (str): path to FITS file

    Returns:
        data (ndarray): image data
        wcs (WCS): WCS object
    """
    hdul = fits.open(fdir)
    wcs = None
    data = None

    for hdu in hdul:
        if hdu.data is None:
            continue
        try:
            w = WCS(hdu.header)
            if w.has_celestial:
                wcs = w
                data = hdu.data
                header = hdu.header
                break
        except Exception as e:
            raise KeyError(f"Error reading WCS: {e}")

    if wcs is None:
        raise RuntimeError("No valid celestial WCS found in FITS.")

    return data, wcs, header


def generate_tnot_plan(target,ra,dec,save_path='./',count=6,interval=300):
    """Generate TNOT observation plan file.

    Args:
        target (str, list): target name(s)
        ra (str,astropy.units, list): Right Ascension(s)
        dec (str,astropy.units, list): Declination(s)
        save_path (str, optional): path to save file. Defaults to './'.
        count (int, optional): exposure counts. Defaults to 6.
        interval (int, optional): exposure time in second. Defaults to 300.
    """
    
    if isinstance(target,str):
        target = [target]
    elif isinstance(target,list):
        pass
    else:
        raise TypeError("target should be str or list of str")
    
    if isinstance(ra,float):
        ra = [ra]*len(target)
        dec = [dec]*len(target)
    elif isinstance(ra,list):
        if len(ra) != len(target) or len(dec) != len(target):
            raise ValueError("length of ra/dec should be equal to target")
    elif isinstance(ra,u.Quantity):
        ra = ra.to(u.deg).value
        dec = dec.to(u.deg).value
        
    if save_path:
        tnow = Time.now()
        yyyy, mm, dd, _, _, _ = tnow.ymdhms
        fname = f"plan_{yyyy}{mm:02d}{dd:02d}.txt"
        
        header_lines = [
            "#dir EP2025",
            "#Filter rp",
            "#Binning 1",
            f"#Count {count}",
            f"#Interval {interval}"
        ]
        
        with open(os.path.join(save_path,fname), "w") as f:
            for line in header_lines:
                f.write(line.rstrip() + "\n")
            f.write("\n")
            
            for i in range(len(target)):
                name = target[i]
                ra_deg = ra[i]
                dec_deg = dec[i]

                # Use SkyCoord to convert to HMS/DMS
                c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='fk5')

                ra_hms = c.ra.to_string(unit=u.hour, sep=':', precision=2, pad=True)
                dec_dms = c.dec.to_string(unit=u.deg, sep=':', precision=2, pad=True, alwayssign=True)

                line = f"{name}\t{ra_hms}\t{dec_dms}"
                f.write(line + "\n")
                print(line)  # For direct viewing
            
            
            
            
def generate_tnot_object_json(target,ra,dec,obs_time=None,significance=None,save_path='./'):
    """Generate TNOT object JSON file(s) and scp to tnot server. This file is used for postprocess in TNOT pipeline.

    Args:
        target (str, list): target name(s)
        ra (str,astropy.units, list): Right Ascension(s)
        dec (str,astropy.units, list): Declination(s)
        obs_time (str, list, optional): EP Obs Time. Defaults to None.
        significance (float, list, optional): EP detection significance. Defaults to None.
        save_path (str, optional): path to save file locally. Defaults to './'.
    """

    if isinstance(target,str):
        target = [target]
    elif isinstance(target,list):
        pass
    else:
        raise TypeError("target should be str or list of str")
    
    if isinstance(ra,float):
        ra = [ra]*len(target)
        dec = [dec]*len(target)
    elif isinstance(ra,list):
        if len(ra) != len(target) or len(dec) != len(target):
            raise ValueError("length of ra/dec should be equal to target")
    elif isinstance(ra,u.Quantity):
        ra = ra.to(u.deg).value
        dec = dec.to(u.deg).value
        
    if obs_time is not None:
        if isinstance(obs_time,str):
            obs_time = [obs_time]*len(target)
        elif isinstance(obs_time,list):
            pass
        else:
            raise TypeError("obs_time should be str or list of str")
    
    if significance is not None:
        if isinstance(significance,(float,int)):
            significance = [significance]*len(target)
        elif isinstance(significance,list):
            pass
        else:
            raise TypeError("significance should be float or list of float")
        
    if save_path:
        for i in range(len(target)):
            # Write individual JSON file for each source
            name = target[i]
            ra_deg = ra[i]
            dec_deg = dec[i]
            t = obs_time[i]
            s = significance[i]
            
            json_outfile = f'{save_path}/{name}_info.json'
            info = {
                "simbad_name": name,
                "ra": ra_deg,
                "dec": dec_deg,
                "obs_time": t,
                "significance": s,
            }
            with open(json_outfile, "w") as json_file:
                json.dump(info, json_file, indent=4)
            print(f"Saved JSON to {json_outfile}")


            #ssh to TNOT server
            subprocess.run([
                "scp", "-P", "5905",
                json_outfile,
                "tnot@119.78.162.172:/home/tnot/EP/plans"
            ], check=True)
      
            
            
            
            
def generate_sitian_plan(target,ra,dec,exptime,expcount,p=6,save_path='./'):
    """Generate Sitian observation plan file.

    Args:
        target (str, list): target name(s)
        ra (str,astropy.units, list): Right Ascension(s)
        dec (str,astropy.units, list): Declination(s)
        exptime (float, int, list): exposure time(s) in second
        expcount (int, list): exposure count(s)
        p (int, optional): Priority. Defaults to 6.
        save_path (str, optional): path to save file. Defaults to './'.
    """
    
    if isinstance(target,str):
        target = [target]
    elif isinstance(target,list):
        pass
    else:
        raise TypeError("target should be str or list of str")
    
    if isinstance(ra,float):
        ra = [ra]*len(target)
        dec = [dec]*len(target)
    elif isinstance(ra,list):
        if len(ra) != len(target) or len(dec) != len(target):
            raise ValueError("length of ra/dec should be equal to target")
    elif isinstance(ra,u.Quantity):
        ra = ra.to(u.deg).value
        dec = dec.to(u.deg).value
        
    if isinstance(exptime,(float,int)):
        exptime = [exptime]*len(target)
    elif isinstance(exptime,list):
        if len(exptime) != len(target):
            raise ValueError("length of exptime should be equal to target")
    else:
        raise TypeError("exptime should be float or list of float/int")
    
    if isinstance(expcount,(float,int)):
        expcount = [expcount]*len(target)
    elif isinstance(expcount,list):
        if len(expcount) != len(target):
            raise ValueError("length of expcount should be equal to target")
    else:
        raise TypeError("expcount should be float or list of float/int")
    
    if save_path:
        tnow = Time.now()
        yyyy, mm, dd, _, _, _ = tnow.ymdhms
        fname = f"sitian_plan_{yyyy}-{mm:02d}-{dd:02d}.txt"
        
        header = "objName\tra\tdec\texpTime\texpCount\ttype\tPriority\n"
        
        with open(os.path.join(save_path,fname), "w") as f:
            f.write(header)
            for i in range(len(target)):
                name = target[i]
                ra_deg = ra[i]
                dec_deg = dec[i]
                exptime_i = exptime[i]
                expcount_i = expcount[i]
                
                # Use SkyCoord to convert to HMS/DMS
                c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='fk5')
                ra_hms = c.ra.to_string(unit=u.hour, sep=':', precision=2, pad=True)
                dec_dms = c.dec.to_string(unit=u.deg, sep=':', precision=2, pad=True, alwayssign=True)
                row = f"{name}\t{ra_hms}\t{dec_dms}\t{exptime_i}\t{expcount_i}\tOBJECT\t{p}\n"

                f.write(row)

            print(f"Observation list saved to {fname}")
        
        


# def get_tnot_data(root_dir="~/optical_data"):
#     """
#     Download and organize TNOT stacked WCS FITS files.
#     """
#     root_dir = Path(root_dir).expanduser().resolve()
#     tmp_dir = root_dir / "tmp_tnot"
#     tmp_dir.mkdir(parents=True, exist_ok=True)

#     REMOTE = "tnot@119.78.162.172:/home/tnot/EP"
#     ssh_port = '5905'

#     print("[1] Querying remote TNOT files ...")

#     # 只列出 yyyyMMdd/stack*wcs.fits
#     cmd_find = [
#         "ssh",
#         '-p',
#         ssh_port,
#         "tnot@119.78.162.172",
#         "find /home/tnot/EP -maxdepth 2 -type f -name 'stack*wcs.fits'"
#     ]

#     result = subprocess.run(
#         cmd_find, capture_output=True, text=True, check=True
#     )

#     remote_files = result.stdout.strip().splitlines()

#     if not remote_files:
#         print("No stack*wcs.fits found on remote server.")
#         return

#     print(f"Found {len(remote_files)} files")

#     # 本地已有目录名（用于匹配）
#     local_dirs = [d for d in root_dir.iterdir() if d.is_dir() and not d.name.startswith("_")]

#     def match_local_dir(name):
#         for d in local_dirs:
#             if d.name.lower() == name.lower():
#                 return d
#         return None

#     unmatched = []

#     for rfile in remote_files:
#         fname = os.path.basename(rfile)
#         local_tmp = tmp_dir / fname

#         print(f"[2] Downloading {fname}")

#         subprocess.run(
#             ["scp",'-P', ssh_port, f"tnot@119.78.162.172:{rfile}", str(local_tmp)],
#             check=True
#         )

#         # 读取 FITS header
#         try:
#             with fits.open(local_tmp) as hdul:
#                 hdr = hdul[0].header
#         except Exception as e:
#             print(f"  [ERROR] Cannot open FITS: {e}")
#             continue

#         name = hdr.get("SIMBAD_NAME") or hdr.get("OBJECT")
#         date_obs = hdr.get("DATE-OBS", "unknown")

#         if name is None:
#             print("  [WARN] No SIMBAD_NAME or OBJECT in header")
#             unmatched.append((fname, "NO_NAME"))
#             continue

#         # 清理名字
#         name_clean = re.sub(r"\s+", "", str(name))
#         date_clean = date_obs.replace(":", "").replace("-", "").replace("T", "_")

#         target_dir = match_local_dir(name_clean)
#         target_dir = target_dir / "TNOT"  if target_dir else None

#         if target_dir is None:
#             print(f"  [UNMATCHED] {name_clean} ({date_obs})")
#             unmatched.append((fname, name_clean))
#             continue
#         elif not os.path.exists(target_dir):
#             os.makedirs(target_dir, exist_ok=True)
        
        
#         new_name = f"TNOT_{name_clean}_{date_clean}.fits"
#         final_path = target_dir / new_name

#         if final_path.exists():
#             print(f"  [SKIP] File already exists: {final_path.name}")
            
#             # 删除临时文件
#             if local_tmp.exists():
#                 local_tmp.unlink()
#                 print(f"         Removed temp file: {local_tmp.name}")
#         else:
#             shutil.move(local_tmp, final_path)

#             # 兜底清理（理论上 move 后不会存在）
#             if local_tmp.exists():
#                 local_tmp.unlink()

#             print(f"  [OK] Saved to {final_path}")

#     # 清理临时目录
#     try:
#         tmp_dir.rmdir()
#     except OSError:
#         subprocess.run(["rm", "-rf", str(tmp_dir)])

#     if unmatched:
#         print("\n===== Unmatched Files =====")
#         for f, n in unmatched:
#             print(f"{f} -> {n}")
    


def get_tnot_data(root_dir=str(DEFAULT_OPTICAL_ROOT), remote_user="tnot", remote_host="119.78.162.172",
                  remote_port=5905, remote_base="/home/tnot/EP"):
    """
    Download TNOT stack FITS files from remote server, organize by target.
    
    Parameters
    ----------
    root_dir : str
        Local root directory to store downloaded data.
    
    Returns
    -------
    local_files : dict
        {target_name: [list of local FITS files]}
    """
    local_files = {}

    # 列出远程服务器所有 stack*_wcs.fits 文件
    remote_cmd = f'ssh -p {remote_port} {remote_user}@{remote_host} "find {remote_base} -type f -name \'*stack*_wcs.fits\'"'
    result = subprocess.run(remote_cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print("Failed to list remote files:", result.stderr)
        return local_files

    remote_files = [f.strip() for f in result.stdout.splitlines()]
    print(f"Found {len(remote_files)} remote files.")

    for rf in remote_files:
        fname = os.path.basename(rf)

        # 尝试从文件名匹配 target
        parts = fname.split("_")
        target = parts[1] if len(parts) >= 3 else None

        out_dir = os.path.join(root_dir, target if target else "UNKNOWN", "TNOT")

        if not os.path.exists(out_dir) or target is None:
            # 如果本地目录不存在，直接跳过
            print(f"[SKIP] Local directory does not exist for target '{rf}': {out_dir}")
            continue
        
        else:
            # 本地目录
            out_dir = os.path.join(root_dir, target, "TNOT")
            os.makedirs(out_dir, exist_ok=True)
            out_file = os.path.join(out_dir, fname)

            # 如果文件已经存在，则跳过下载
            if os.path.exists(out_file):
                print(f"[SKIP] File already exists: {out_file}")
            else:
                # 下载
                scp_cmd = ["scp", "-P", str(remote_port),
                        f"{remote_user}@{remote_host}:{rf}", out_file]
                subprocess.run(scp_cmd, check=True)
                print(f"[Downloading] {fname} → {out_file}")

            # 读取 FITS 文件
            try:
                data, wcs, header = read_image(out_file)
            except Exception as e:
                print(f"[ERROR] Failed to read {out_file}: {e}")
                continue

            local_files.setdefault(target, []).append(out_file)

    return local_files



def get_sitian_data(root=str(DEFAULT_OPTICAL_ROOT)):
    """
    Download Sitian (司天) data via SSH and organize by source name.

    Parameters
    ----------
    root : str or Path
        Local root directory to store data.
        Structure: root/sourcename/sitian/*.wcs.fits
    """

    root = Path(root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)

    ssh_host = "obs@www.xinglong-naoc.cn"
    ssh_port = "20210"

    remote_pattern = "/data/raw/ST001/gal_results/reddir/*/COM_IMAGE/*EP*.fits"

    print("🔍 Querying remote Sitian files...")

    # 1. List remote files
    cmd_ls = [
        "ssh",
        "-p", ssh_port,
        ssh_host,
        f"ls {remote_pattern}"
    ]

    try:
        result = subprocess.run(
            cmd_ls,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
    except subprocess.CalledProcessError as e:
        print("❌ Failed to list remote files")
        print(e.stderr)
        return

    files = [f.strip() for f in result.stdout.splitlines() if f.strip()]

    if not files:
        print("⚠️ No Sitian EP files found.")
        return

    print(f"📦 Found {len(files)} files on server")

    # 2. Process each file
    for remote_file in files:
        fname = os.path.basename(remote_file)

        # Expected format: yyyy-mm-dd_sourcename**.fits
        if "_" not in fname:
            print(f"⚠️ Skip unrecognized filename: {fname}")
            continue

        date_part, rest = fname.split("_", 1)
        sourcename, rest = rest.split("_", 1)
        #sourcename = rest.split(".")[0]

        # Local directory: root/sourcename/sitian
        target_dir = root / sourcename / "sitian"
        target_dir.mkdir(parents=True, exist_ok=True)

        target_file = target_dir / fname.replace(".fits", ".wcs.fits")

        # Skip if file exists
        if target_file.exists():
            print(f"⏭️ Skip existing: {target_file.name}")
            continue

        print(f"⬇️ Downloading {fname} → {target_dir}")

        # scp download
        cmd_scp = [
            "scp",
            "-P", ssh_port,
            f"{ssh_host}:{remote_file}",
            str(target_file)
        ]

        try:
            subprocess.run(cmd_scp, check=True)
        except subprocess.CalledProcessError:
            print(f"❌ Failed to download {fname}")
            continue

    print("✅ Sitian data download finished.")
    
    
def check_source_dirs(root_dir=str(DEFAULT_OPTICAL_ROOT), meta_file=str(DEFAULT_CANDIDATES_FILE)):
    """
    Check if source directories exist for each source in the metadata file.

    Parameters
    ----------
    root_dir : str or Path
        Root directory containing source subdirectories.
    meta_file : str or Path
        Path to the metadata CSV file with a 'name' column.

    Returns
    -------
    missing_sources : list
        List of source names without corresponding directories.
    """
    import pandas as pd

    root_dir = Path(root_dir).expanduser().resolve()
    meta_file = Path(meta_file).expanduser().resolve()

    df = pd.read_csv(meta_file)
    source_names = df['EP Name'].tolist()

    missing_sources = []

    for name in source_names:
        source_dir = root_dir / name
        if not source_dir.exists() or not source_dir.is_dir():
            missing_sources.append(name)
            os.mkdir(source_dir)
            print("Created missing directory:", source_dir)
            
        if not (source_dir / "ps.csv").exists():
            from catalog import PS, LS
            try:
                client = PS()
                ra = df[df['EP Name'] == name]['RA'].values[0]
                dec = df[df['EP Name'] == name]['Dec'].values[0]
                client.get_catalog(ra, dec, save_path=source_dir)
            except:
                client = PS()
                ra = df[df['EP Name'] == name]['RA'].values[0]
                dec = df[df['EP Name'] == name]['Dec'].values[0]
                client.get_catalog(ra, dec, save_path=source_dir)
                
            
            
            


def fits_plot(
    fdir,
    radius: float = 20,   # arcsec
    source_dir='/home/liangrd/Follow_up/Candidates.csv',
    size=100,              # arcsec
    contrast=0.25,
    show_ps=False,
    ps_filter='r',   
    ra=None,
    dec=None,
    **imshow_kwargs,
):
    """
    Visualize FITS image centered on target using Cutout2D.
    Optionally show Pan-STARRS image alongside.

    Parameters
    ----------
    fdir : str
        Path to FITS image
    radius : float
        Circle radius in arcsec
    source_dir : str
        CSV file with columns: EP Name, RA, Dec
    size : int
        Cutout size in arcsecs (size x size)
    show_ps : bool
        Whether to show Pan-STARRS image
    ps_dir : str
        Pan-STARRS FITS path; if None, attempt download
    ra, dec : float
        Optional RA/Dec center (deg)
    """

    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.nddata import Cutout2D
    from astropy.visualization import ZScaleInterval
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    import os
    from reproject import reproject_interp
    import pandas as pd

    if not isinstance(size,u.Quantity):
        size = size * u.arcsec
    
    if not isinstance(radius,u.Quantity):
        radius = radius * u.arcsec

    # ---------- 1. Open FITS & find valid WCS ----------
    hdul = fits.open(fdir)
    wcs = None
    data = None
    header = None

    for hdu in hdul:
        if hdu.data is None:
            continue
        try:
            w = WCS(hdu.header)
            if w.has_celestial:
                wcs = w
                data = hdu.data
                header = hdu.header
                break
        except Exception as e:
            raise KeyError(f"Error reading WCS: {e}")

    if wcs is None:
        raise RuntimeError("No valid celestial WCS found in FITS.")

    # ---------- 2. Get target center ----------
    if ra is not None and dec is not None:
        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    else:
        # fallback to CSV
        cat = pd.read_csv(source_dir)
        if 'OBJECT' in header:
            target = header['OBJECT']
        else:
            target = cat.columns[0]  # just to avoid crash

        row = cat[cat['EP Name'] == target]
        if len(row)==0:
            raise ValueError(f"Target {target} not found in {source_dir}")
        coord = SkyCoord(ra=float(row.iloc[0]['RA'])*u.deg,
                         dec=float(row.iloc[0]['Dec'])*u.deg)

    # ---------- 3. Main cutout ----------
    cutout = Cutout2D(data, coord, size, wcs=wcs)
    zscale = ZScaleInterval(contrast=contrast)
    vmin, vmax = zscale.get_limits(cutout.data)

    # ---------- 4. Pan-STARRS cutout ----------
    ps_cutout = None
    if show_ps is not None:
        ps_dir = str(DEFAULT_OPTICAL_ROOT / str(target))
        ps_fits = os.path.join(ps_dir,f"ps1_{ps_filter}_ref.fits")
        if not os.path.exists(ps_fits):
            # assume ps_dir is directory
            download_ps(coord.ra.deg, coord.dec.deg, ps_filter, size=1024, save_dir=ps_dir)

        with fits.open(ps_fits) as hdu_ps:
            ps_data = hdu_ps[0].data
            ps_wcs = WCS(hdu_ps[0].header)
        ps_cutout0 = Cutout2D(ps_data, coord, size, wcs=ps_wcs)
        # reproject PS onto main cutout WCS for alignment
        ps_reproj, _ = reproject_interp((ps_cutout0.data, ps_cutout0.wcs), cutout.wcs, shape_out=cutout.data.shape)
        ps_cutout = ps_reproj
        
        zscale = ZScaleInterval(contrast=contrast)
        vmin_ps, vmax_ps = zscale.get_limits(ps_cutout.data)

    # ---------- 5. Plotting ----------
    if ps_cutout is not None:
        fig = plt.figure(figsize=(10.5, 5))
        from matplotlib.gridspec import GridSpec
        gs = GridSpec(1, 2, wspace=0.01)

        ax1 = fig.add_subplot(gs[0, 0], projection=cutout.wcs)
        ax2 = fig.add_subplot(gs[0, 1], projection=cutout.wcs, sharex=ax1, sharey=ax1)

        images = [cutout.data, ps_cutout]
        titles = [os.path.basename(fdir), "Pan-STARRS"]

        for i, ax, img, title, v0, v1 in zip([0,1], [ax1, ax2], images, titles, [vmin,vmin_ps],[vmax,vmax_ps]):
            ax.imshow(img, origin='lower', cmap='gray', vmin=v0, vmax=v1, **imshow_kwargs)
            circ = Circle((coord.ra.deg, coord.dec.deg),
                          radius.to(u.deg).value,
                          edgecolor='yellow', facecolor='none', lw=1.5,
                          transform=ax.get_transform('world'))
            ax.add_patch(circ)
            ax.set_title(title)
            ax.coords.grid(color='white', ls='dotted', alpha=0.6)
            ax.set_xlabel('RA')
            ax.tick_params(axis='both', which='both', direction='in')  # 去掉 ticks
            if i == 0:
                ax.set_ylabel('Dec')
            else:
                ax.coords[0].set_ticklabel_visible(False)
                ax.coords[1].set_ticklabel_visible(False)
                ax.set_yticklabels([])        # 右图隐藏 ytick
                ax.set_ylabel('')             # 右图隐藏 ylabel

    else:
        fig = plt.figure(figsize=(6, 5))
        ax1 = fig.add_subplot(111, projection=cutout.wcs)
        ax1.imshow(cutout.data, origin='lower', cmap='gray', vmin=vmin, vmax=vmax, **imshow_kwargs)
        circ = Circle((coord.ra.deg, coord.dec.deg),
                      radius.to(u.deg).value,
                      edgecolor='yellow', facecolor='none', lw=1.5,
                      transform=ax1.get_transform('world'))
        ax1.add_patch(circ)
        ax1.set_title(os.path.basename(fdir))
        ax1.set_xlabel('RA')
        ax1.set_ylabel('Dec')
        ax1.coords.grid(color='white', ls='dotted', alpha=0.6)

    plt.tight_layout()
    plt.show()

    return fig
    
    


def plot_lunar_distance(
    ra,
    dec,
    start_time=None,
    ndays=3,
    step_hours=6,
    name=None,
    threshold=30,
    plot=True
):
    """
    Calculate and plot lunar distance for a given sky position
    over the next month.

    Parameters
    ----------
    ra : float or str
        Right Ascension (deg or sexagesimal string)
    dec : float or str
        Declination (deg or sexagesimal string)
    start_time : str or astropy.time.Time, optional
        Start time (default: now, UTC)
    ndays : int, optional
        Number of days to calculate (default: 30)
    step_hours : float, optional
        Time step in hours (default: 6)

    Returns
    -------
    time : astropy.time.Time
        Time array
    separation : astropy.units.Quantity
        Angular separation (deg)
    """

    # Start time
    if start_time is None:
        t0 = Time.now()
    else:
        t0 = Time(start_time)

    # Time grid
    times = t0 + np.arange(0, ndays * 24, step_hours) * u.hour

    # Target coordinate
    target = SkyCoord(ra=ra, dec=dec, unit=u.deg)

    # Moon position
    moon = get_body("moon", times)

    # Angular separation
    separation = moon.separation(target)

    delta_theta = separation.deg - threshold

    # ---- Plot ----
    if plot:
        plt.figure(figsize=(5,3))
        plt.plot(times.datetime, separation.deg)
        plt.axhline(30,ls='--',lw=1,color='C1')
        plt.axhline(20,ls='--',lw=1,color='r')
        plt.xlabel("Date (UTC)")
        plt.ylabel("Lunar distance (deg)")
        plt.title(f"Lunar Distance of {name} (RA={ra}, Dec={dec})" if name else f"Lunar Distance (RA={ra}, Dec={dec})")
        plt.tight_layout()
        plt.show()

    return times, separation
    
    
    
    
def show_shift(root, save_dir='./', source_dir='/home/liangrd/Follow_up/Candidates.csv',
               ploter='matplotlib'):
    """
    Show observation timelines of ALL targets in ONE figure.

    Directory structure:
        root/target/telescope/**/*.fits or fits.fz
    """

    import os
    import glob
    import numpy as np
    import matplotlib.pyplot as plt

    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.time import Time
    from astropy.table import Table

    os.makedirs(save_dir, exist_ok=True)
    
    sources = Table.read(source_dir)

    targets = sorted([d for d in os.listdir(root)
                      if os.path.isdir(os.path.join(root, d))])

    all_rows = []

    # =========================
    # Collect all observations
    # =========================
    for target in targets:
        target_dir = os.path.join(root, target)
        try:
            t0 = Time(sources[sources['EP Name'] == target]['Obs Time'][0])
        except:
            t0 = None

        fits_files = glob.glob(os.path.join(target_dir, '**', '*.fits'), recursive=True)
        fits_files += glob.glob(os.path.join(target_dir, '**', '*.fits.fz'), recursive=True)
        fits_files = [f for f in fits_files if 'ref' not in f.lower()]  # 排除参考图
        # print(target)
        # print(fits_files)

        for f in fits_files:
            try:
                with fits.open(f) as hdul:
                    valid = False
                    header = None

                    for hdu in hdul:
                        if hdu.data is None:
                            continue
                        try:
                            w = WCS(hdu.header)
                            if w.has_celestial:
                                wcs = w
                                data = hdu.data
                                header = hdu.header
                                break
                        except Exception as e:
                            raise KeyError(f"Error reading WCS: {e}")

                    if wcs is None:
                        raise RuntimeError("No valid celestial WCS found in FITS.")

                    date_obs = header.get('DATE-OBS')
                    if date_obs is None:
                        date_obs = header.get('DATE')
                        if date_obs is None:
                            print(f'[WARN] No DATE-OBS in {f}')
                            continue
                    
                    ts = pd.Timestamp(date_obs)
                    date_obs = ts.strftime("%Y-%m-%dT%H:%M:%S.%f")
                    try:
                        t = Time(date_obs, format='isot', scale='utc')
                    except:
                        t = Time(date_obs)
                        
                    if t0 is None:
                        t0 = t
                    elif t.mjd < t0.mjd:
                        t0 = t
                        
                    dt = t.mjd - t0.mjd
                        
                    

                    band = _get_filter_from_header_or_filename(header, f)
                    if band == 'UNKNOWN':
                        print(f'[WARN] Unknown band in {f}')
                    else:
                        band = band[0] #first alphabet character
                    rel = os.path.relpath(f, target_dir)
                    telescope = rel.split(os.sep)[0]

                    all_rows.append({
                        'target': target,
                        'mjd': t.mjd,
                        'time_iso': t.isot,
                        'telescope': telescope,
                        'band': band,
                        'file': f,
                        'dt':float(dt)
                    })
                    #print(target,all_rows[-1])

            except Exception as e:
                print(f'[ERROR] {f}: {e}')
                continue

    if len(all_rows) == 0:
        print('[WARN] No valid observations found.')
        return None

    # =========================
    # Build table
    # =========================
    tab = Table(rows=all_rows)
    tab.sort('mjd')

    tab_path = os.path.join(save_dir, 'all_targets_timeline.csv')
    tab.write(tab_path, format='csv', overwrite=True)

    # =========================
    # Plot ALL targets together
    # =========================

    targets_unique = list(dict.fromkeys(tab['target']))
    target_y = {t: i for i, t in enumerate(targets_unique)}

    telescopes = sorted(set(tab['telescope']))
    bands = sorted(set(tab['band']))

    if ploter == 'matplotlib':

        import matplotlib.pyplot as plt

        # color for telescope
        tel_color = {tel: plt.cm.tab10(i % 10) for i, tel in enumerate(telescopes)}

        # marker for band
        marker_list = ['o', 's', '^', 'D', 'v', 'P', 'X', '*', '<', '>']
        band_marker = {band: marker_list[i % len(marker_list)]
                       for i, band in enumerate(bands)}

        plt.figure(figsize=(12, max(4, 0.4 * len(targets_unique))))

        for tel in telescopes:
            for band in bands:
                mask = ((tab['telescope'] == tel) &
                        (tab['band'] == band))
                if not np.any(mask):
                    continue

                x = tab['dt'][mask]
                y = [target_y[t] for t in tab['target'][mask]]

                plt.scatter(
                    x, y,
                    color=tel_color[tel],
                    marker=band_marker[band],
                    s=40,
                    alpha=0.8
                )

        plt.yticks(range(len(targets_unique)), targets_unique)
        plt.xlabel(r'$T-T_{0}$ (days)')
        plt.ylabel('Target')
        plt.title('Observation Timeline of All Targets')
        plt.grid(alpha=0.3)

        # telescope legend
        tel_handles = [
            plt.Line2D([], [], marker='o', linestyle='',
                       color=tel_color[tel], label=tel)
            for tel in telescopes
        ]

        # band legend
        band_handles = [
            plt.Line2D([], [], marker=band_marker[band], linestyle='',
                       color='k', label=band)
            for band in bands
        ]

        leg1 = plt.legend(
            handles=tel_handles,
            title='Telescope',
            bbox_to_anchor=(1.02, 1),
            loc='upper left'
        )
        plt.gca().add_artist(leg1)

        plt.legend(
            handles=band_handles,
            title='Band',
            bbox_to_anchor=(1.02, 0.5),
            loc='upper left'
        )

        plt.tight_layout()
        fig_path = os.path.join(save_dir, 'all_targets_timeline.png')
        plt.savefig(fig_path, dpi=200)
        plt.close()

    elif ploter == 'plotly':

        import plotly.graph_objects as go
        from plotly.colors import qualitative

        df = tab.to_pandas()

        targets_unique = list(dict.fromkeys(df['target']))
        target_y = {t: i for i, t in enumerate(targets_unique)}

        df['yval'] = df['target'].map(target_y)

        telescopes = sorted(df['telescope'].unique())
        bands = sorted(df['band'].unique())

        # -------------------------
        # Color map (like matplotlib tab10)
        # -------------------------
        color_list = qualitative.Plotly
        tel_color = {
            tel: color_list[i % len(color_list)]
            for i, tel in enumerate(telescopes)
        }

        # -------------------------
        # Marker map (similar to matplotlib)
        # -------------------------
        marker_list = [
            'circle', 'square', 'triangle-up', 'diamond',
            'triangle-down', 'cross', 'x', 'star',
            'triangle-left', 'triangle-right'
        ]

        band_marker = {
            band: marker_list[i % len(marker_list)]
            for i, band in enumerate(bands)
        }

        fig = go.Figure()

        # ======================================================
        # Draw each telescope + band combination separately
        # ======================================================
        for tel in telescopes:
            for band in bands:

                sub = df[(df['telescope'] == tel) &
                        (df['band'] == band)]

                if len(sub) == 0:
                    continue

                fig.add_trace(
                    go.Scatter(
                        x=sub['dt'],
                        y=sub['yval'],
                        mode='markers',
                        marker=dict(
                            color=tel_color[tel],
                            symbol=band_marker[band],
                            size=8,
                            opacity=0.8,
                            line=dict(width=0)
                        ),
                        name=f"{tel} ({band})",
                        hovertemplate=
                            "<b>%{text}</b><br>" +
                            "T-T0: %{x:.2f} d<br>" +
                            "ISO: %{customdata[0]}<br>" +
                            "File: %{customdata[1]}<extra></extra>",
                        text=sub['target'],
                        customdata=np.stack(
                            (sub['time_iso'], sub['file']), axis=-1
                        ),
                    )
                )

        # -------------------------
        # Axis formatting
        # -------------------------
        fig.update_layout(
            title="Observation Timeline of All Targets",
            xaxis_title="T - T0 (days)",
            yaxis=dict(
                tickmode='array',
                tickvals=list(range(len(targets_unique))),
                ticktext=targets_unique,
            ),
            height=max(400, 30 * len(targets_unique)),
            legend_title="Telescope (Band)",
            template="simple_white"
        )

        fig.update_xaxes(showgrid=True)
        fig.update_yaxes(showgrid=False)

        fig_path = os.path.join(save_dir, 'all_targets_timeline.html')
        fig.write_html(fig_path)

        print(f"[OK] Saved interactive figure to {fig_path}")

    print(f'[OK] Saved:')
    print(f'  - {tab_path}')
    print(f'  - {fig_path}')

    return tab


def _get_filter_from_header_or_filename(header, filepath):
    """
    Get filter name from FITS header or filename.

    Priority:
    1. header['FILTER']
    2. filename pattern _[ugrizyw]_
    """
    # 1. header
    band = header.get('FILTER')
    if band is not None:
        band = str(band).strip()
        if band != '' and band.upper() != 'NONE':
            return band

    # 2. filename fallback
    fname = os.path.basename(filepath).lower()
    m = re.search(r'_([ugrizyw])_', fname)
    if m:
        return m.group(1)

    return 'UNKNOWN'



def calculate_observation_stats(
        obs_file='/home/liangrd/Follow_up/results/all_targets_timeline.csv',
        candidates_file='/home/liangrd/Follow_up/Candidates.csv',
        output_file='/home/liangrd/Follow_up/results/all_targets_meta.csv'):
    """
    Calculate observation statistics for each target (ALL targets in Candidates.csv).
    Targets without observations will have zero values.

    Output columns:
        target, T0, age, latest_watch, nwatch, nwatch5
    """

    import pandas as pd
    import numpy as np
    from astropy.time import Time

    # -------------------------
    # Read files
    # -------------------------
    obs = pd.read_csv(obs_file)
    cand = pd.read_csv(candidates_file)

    cand.rename(columns={'EP Name': 'target'}, inplace=True)

    # Convert T0 to MJD
    cand['T0'] = [Time(t).mjd for t in cand['Obs Time']]

    tnow = Time.now().mjd

    # Ensure required columns
    if 'target' not in obs.columns or 'mjd' not in obs.columns:
        raise ValueError("Observation file must contain columns: target, mjd")

    if 'target' not in cand.columns or 'T0' not in cand.columns:
        raise ValueError("Candidates file must contain: EP Name, Obs Time")

    results = []

    # ======================================================
    # Loop over ALL targets in Candidates.csv
    # ======================================================
    for _, row in cand.iterrows():

        target = row['target']
        T0 = row['T0']

        # Select observations of this target
        group = obs[obs['target'] == target]

        if len(group) == 0:
            # ----------------------------------
            # No observation case
            # ----------------------------------
            nwatch = 0
            nwatch5 = 0
            latest_watch = 0
        else:
            group = group.sort_values('mjd')
            mjds = group['mjd'].values

            nwatch = len(mjds)

            # Latest watch = interval between last two obs
            if nwatch > 1:
                latest_watch = tnow - np.max(mjds)
            else:
                latest_watch = 0

            # Observations after T0+5
            nwatch5 = int(np.sum(mjds >= (T0 + 5)))

        # Age = current time - T0 (always defined)
        age = tnow - T0

        results.append({
            'target': target,
            'T0': T0,
            'age': age,
            'latest_watch': latest_watch,
            'nwatch': nwatch,
            'nwatch5': nwatch5
        })

    # -------------------------
    # Save output
    # -------------------------
    result_df = pd.DataFrame(results)
    result_df.to_csv(output_file, index=False)

    print(f"[OK] Saved results to {output_file}")

    return result_df


def show_obs_pies(
    obs_file='/home/liangrd/Follow_up/results/all_targets_timeline.csv',
    save_path='/home/liangrd/Follow_up/results'
):
    import plotly.express as px

    df = pd.read_csv(obs_file)

    # ===============================
    # Pie 1: Distribution by Telescope
    # ===============================
    telescope_counts = df['telescope'].value_counts().reset_index()
    telescope_counts.columns = ['telescope', 'count']

    fig_tel = px.pie(
        telescope_counts,
        values='count',
        names='telescope',
        title='Observation Distribution by Telescope'
    )

    fig_tel.update_traces(textposition='inside', textinfo='percent+label')
    fig_tel.update_layout(showlegend=False, width=600, height=600)

    fig_tel_path = os.path.join(save_path, 'observation_distribution_telescope.html')
    fig_tel.write_html(fig_tel_path)


    # ===============================
    # Pie 2: Distribution by Target
    # ===============================
    target_counts = df['target'].value_counts().reset_index()
    target_counts.columns = ['target', 'count']

    fig_target = px.pie(
        target_counts,
        values='count',
        names='target',
        title='Observation Distribution by Target'
    )

    fig_target.update_traces(textposition='inside', textinfo='percent+label')
    fig_target.update_layout(showlegend=False, width=600, height=600)

    fig_target_path = os.path.join(save_path, 'observation_distribution_target.html')
    fig_target.write_html(fig_target_path)

    return fig_tel, fig_target
    
    
def show_cumulative_observations(
    candidates_file='/home/liangrd/Follow_up/Candidates.csv',
    obs_file='/home/liangrd/Follow_up/results/all_targets_timeline.csv',
    save_path='/home/liangrd/Follow_up/results'
):
    import plotly.express as px

    # ---- Cumulative number of unique events ----
    candidates = pd.read_csv(candidates_file)
    # Ensure Obs Time is datetime
    candidates['Obs Time'] = pd.to_datetime(candidates['Obs Time'])
    candidates = candidates.sort_values('Obs Time')
    
    # Cumulative count = number of rows so far
    candidates['cum_events'] = range(1, len(candidates)+1)

    fig1 = px.line(candidates, x='Obs Time', y='cum_events',
                   title='Number eFXTs',
                   labels={'Obs Time':'Obs Time (UTC)', 'cum_events':'N'})
    fig1.update_traces(mode='lines')
    fig1.update_layout(width=700, height=500)
    fig1_path = os.path.join(save_path, 'cumulative_events.html')
    fig1.write_html(fig1_path)
    # fig1.show()

    # ---- Cumulative number of observations ----
    obs = pd.read_csv(obs_file)
    obs['time_iso'] = pd.to_datetime(obs['time_iso'])
    obs = obs.sort_values('time_iso')
    
    # Each row is one observation, so cumulative count is 1,2,3,...
    obs['cum_obs'] = range(1, len(obs)+1)

    fig2 = px.line(obs, x='time_iso', y='cum_obs',
                   title='Number of Observations',
                   labels={'time_iso':'Obs Time (UTC)', 'cum_obs':'Ns'})
    fig2.update_traces(mode='lines')
    fig2.update_layout(width=700, height=500)
    fig2_path = os.path.join(save_path, 'cumulative_observations.html')
    fig2.write_html(fig2_path)
    # fig2.show()