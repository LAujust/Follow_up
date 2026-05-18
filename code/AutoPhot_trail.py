from autophot import AutomatedPhotometry

# Load default configuration
config = AutomatedPhotometry.load()

# Set basic parameters
config["outdir_name"] = "REDUCED"
config["wdir"] = "/home/liangrd/optical_data/EP260131a/NTT/"
config["fits_dir"] = "/home/liangrd/optical_data/EP260131a/NTT/"
config["target_ra"] = 149.9025
config["target_dec"] = -3.3073
config["template_subtraction"]["do_subtraction"] = False
config["catalog"]["use_catalog"] = {
    "griz": "panstarrs",
    "u": "skymapper", 
    "UBVRI": "apass",
}
config['catalog']

config["cosmic_rays"]["remove_cmrays"] = False
config["wcs"]["redo_wcs"] = False



# Run photometry
output_file = AutomatedPhotometry.run_photometry(default_input=config)
print(f"Results saved to: {output_file}")