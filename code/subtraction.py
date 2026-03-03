import os
import stat
import subprocess
from astropy.io import fits

__all__ = ['run_hotpants_single', 'run_pyzogy']


def run_hotpants_single(science_file, ref_file, hotpants_dir='/Volumes/T7/Shared_Files/EP/Projects/hotpants', ns=7, tu=110000, iu=60000, ssig=5, tl=-2900, il=40000, rss=20, r=7, verbose=True, **kw):
    
    """
    Run HOTPANTS difference imaging on a single science FITS file.

    Parameters
    ----------
    science_file : str
        Path to the science FITS file.
    ref_file : str
        Path to the reference/template FITS file.
    hotpants_dir : str
        Directory containing the HOTPANTS executable.
    ns : int
        Kernel half-size (nsx = nsy = ns).
    tu : int
        Upper threshold.
    iu : int
        Input upper limit.
    ssig : int
        Sigma clipping.
    verbose : bool
        Whether to print the command.

    Returns
    -------
    out_file : str
        Path to the generated difference image.
    """

    # 确保 HOTPANTS 可执行
    os.chmod(hotpants_dir, os.stat(hotpants_dir).st_mode | stat.S_IXUSR)
    os.chdir(os.path.dirname(hotpants_dir))

    out_file = science_file.replace('cutout', 'diff')

    cmd = [
        './hotpants',
        '-inim', science_file,
        '-tmplim', ref_file,
        '-outim', out_file,
        '-n', 't',
        '-c', 't',
        '-ssig', str(ssig),
        '-tu', str(tu),
        '-iu', str(iu),
        '-nsx', str(ns),
        '-nsy', str(ns),
        '-tl', str(tl),
        '-il', str(il),
        '-rss', str(rss),
        '-r', str(r)
    ]
    
    # 添加 **kw 参数
    for k, v in kw.items():
        cmd.append(f'-{k}')
        cmd.append(str(v))

    if verbose:
        print(" ".join(cmd))

    subprocess.run(cmd, check=True)
    return out_file


import os
import subprocess
from pathlib import Path


def run_pyzogy(
    sci_fits,
    ref_fits,
    outdir,
    pyzogy_dir,
    prefix="zogy",
    psf_sci=None,
    psf_ref=None,
    extra_args=None,
):
    """
    Run PyZOGY image subtraction.

    Parameters
    ----------
    sci_fits : str
        Science FITS image path
    ref_fits : str
        Reference FITS image path
    outdir : str
        Output directory
    pyzogy_dir : str
        Path to PyZOGY repository
    prefix : str
        Output file prefix
    psf_sci : str or None
        Science PSF FITS (optional)
    psf_ref : str or None
        Reference PSF FITS (optional)
    extra_args : list or None
        Extra command-line arguments for PyZOGY

    Returns
    -------
    dict
        Paths of main output products
    """

    sci_fits = Path(sci_fits).resolve()
    ref_fits = Path(ref_fits).resolve()
    outdir = Path(outdir).resolve()
    pyzogy_dir = Path(pyzogy_dir).resolve()

    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python",
        str(pyzogy_dir / "pyzogy.py"),
        "--science", str(sci_fits),
        "--reference", str(ref_fits),
        "--outdir", str(outdir),
        "--prefix", prefix,
    ]

    if psf_sci:
        cmd += ["--psf-science", str(psf_sci)]
    if psf_ref:
        cmd += ["--psf-reference", str(psf_ref)]

    if extra_args:
        cmd += extra_args

    print("Running PyZOGY:")
    print(" ".join(cmd))

    subprocess.run(cmd, check=True)

    outputs = {
        "diff": outdir / f"{prefix}_diff.fits",
        "snr": outdir / f"{prefix}_snr.fits",
        "scorr": outdir / f"{prefix}_scorr.fits",
    }

    return outputs