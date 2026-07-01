
ZP_BRIGHT_MAG_LIMIT = 13.5
ZP_FAINT_MAG_LIMIT = 20.5
ZP_THRESHOLD_LIMIT = 3.0
ZP_MAX_FLAGS = 32 # 0 for best, 32
ZP_EXPORT_SOURCE_TABLES = True
bkg_fast_mode = False

def build_config(
    *,
    fits_dir: Path,
    target_name: str,
    ra_deg: float,
    dec_deg: float,
    catalog_name: str,
    custom_catalog_fpath: str | None = None,
    selected_filters: list[str] | None = None,
) -> dict:
    config = AutomatedPhotometry.load()

    # Basic target and I/O settings.
    config["fits_dir"] = str(fits_dir)
    config["target_name"] = target_name
    config["target_ra"] = ra_deg
    config["target_dec"] = dec_deg
    config["wdir"] = wdir
    config["outdir_name"] = f"{fits_dir.name}_{catalog_name}"
    config["nCPU"] = 1
    config["restart"] = True

    if selected_filters:
        config["select_filter"] = True
        config["selected_filters"] = selected_filters
    # Keep executable paths explicit so this script is robust across shells.
    config["wcs"]["sextractor_exe_loc"] = "/opt/homebrew/bin/sex"
    config["wcs"]["scamp_exe_loc"] = "/home/liangrd/anaconda3/envs/autophot/bin/scamp"
    
    # WCS, False if accurate
    config["wcs"]["redo_wcs"] = False
    config["wcs"]["refine_after_solve"] = False
    config["wcs"]["apply_solved_to_fits"] = False
    config["catalog"]["use_catalog"] = catalog_name
    config["catalog"]["catalog_radius"] = 7.0  # arcmin for normal Catalog.download() paths
    if catalog_name == "custom" and custom_catalog_fpath:
        config["catalog"]["catalog_custom_fpath"] = custom_catalog_fpath
    config["cosmic_rays"]["remove_cmrays"] = False
    config["source_detection"]["scale_multipler"] = 3
    config["photometry"]["annulus_gap_fwhm"] = 1.0 # gap between aperture radius and sky annulus in FWHM
    config["photometry"]["apply_colorterms"] = False
    config["photometry"]["psf_snr_min"] = 3.0
    config["photometry"]["psf_threshold_limit"] = 3.0
    config["photometry"]["psf_isolation_radius_fwhm"] = 2.0
    config["photometry"]["psf_oversample"] = 2 #PSF 模型内部用 0.5 pixel 的细网格
    config["photometry"]["crowded_field"] = True
    config["photometry"]["crowded_optimum_radius_fwhm"] = 1.
    config["photometry"]["fitting_xy_bounds"] = 0.5 #pixel
    config["background"]["fast_mode"] = bkg_fast_mode
    config["background"]["local_cutout_annulus_sigma_mult"] = 3 # raise the local cutout floor so the background-like median is >= (this * annulus_std), where annulus_std is sigma-clipped std in an annulus around the target on the locally subtracted cutout.
    config["zeropoint"]["reject_nonlinear_sources"] = True
    config["zeropoint"]["min_source_no"] = 5
    config["zeropoint"]["bright_mag_limit"] = ZP_BRIGHT_MAG_LIMIT
    config["zeropoint"]["faint_mag_limit"] = ZP_FAINT_MAG_LIMIT
    config["zeropoint"]["threshold_limit"] = ZP_THRESHOLD_LIMIT
    config["zeropoint"]["max_flags"] = ZP_MAX_FLAGS
    config["zeropoint"]["export_source_tables"] = ZP_EXPORT_SOURCE_TABLES
    config["zeropoint"]["nonlinear_peak_frac"] = 0.95 # fraction of saturate
    config["zeropoint"]["saturation_peak_frac"] = 0.95 # fraction of saturate
    return config