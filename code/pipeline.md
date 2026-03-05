# Pipeline Project

## Goal

Perform co-add and photometry for all non-ref data in ~/optical_data using function defined in /code/coadd.py, /code/photometry.py, for which examples can be found in /code/Coadd.ipynb and Unified_Photometry.ipynb. Automatically scan files, skip processed files, only process new files and update results. DO all files if no files are processed to avoid repeated works. Skip failed files, write logs. 

## Input

- `config_file.json`:  pre-defined parameters for different telescopes

## Step

- Read Candidates.csv
- For each target, each telescopes, coadd data within one day, and store in ~/optical_data/[target]/[pipeline]/[telescope], skip process files that have already been coadded; ra, dec can be read from Candidates.csv, use position of optical counterpart if available; cat_dir use ~/optical_data/[target]/ps.csv or ls.csv
- perform photometry for coadded images, and save results to a .csv file stored in ~/optical_data/[target]/[pipeline]/photometry.csv, figures saved to ~/optical_data/[target]/[pipeline]/[telescope]/photometry/; skip failed files and files that have been already done photometry
- Log file should be stored in ~/optical_data/[target]/[pipeline]/logs/[time].log

## config_file.json

format follows:

```[json]
{
    telescope:{
        'fwhm':xx,
        'r_ap':xx,
        'r_in':xx,
        'r_out':xx,
        'fit_shape':xx,
        'Forced':False,
        ...
    }
}
```

that store parameters input to coadd() and Photometry.psf_photometry() and aperture_photometry()