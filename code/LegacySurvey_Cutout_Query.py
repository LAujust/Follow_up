import os
import numpy as np
import pandas as pd
import pyvo

# Example input parameters
RA = 149.9025  # degrees
Dec = -3.3073   # degrees
err = 250          # arcsec search radius
PATH = "/home/liangrd/optical_data/EP260131a"        # directory to save CSV

# Connect to the NOIRLab TAP service
tap = pyvo.dal.TAPService("https://datalab.noirlab.edu/tap")

# Use bounding box instead of cone search
ra_center, dec_center = RA, Dec
radius = err / 3600  # convert arcsec to degrees

query = f"""
SELECT TOP 4000
  ra, dec, ls_id, flux_g, flux_ivar_g,
  flux_i, flux_ivar_i, flux_r, flux_ivar_r,
  flux_z, flux_ivar_z
FROM ls_dr10.tractor
WHERE ra BETWEEN {ra_center - radius} AND {ra_center + radius}
  AND dec BETWEEN {dec_center - radius} AND {dec_center + radius}
"""

# Run query
results = tap.search(query)
df = results.to_table().to_pandas()
print(df.head())

# Filter out non-positive flux or ivar (to avoid log errors)
bands = ['g', 'r', 'i', 'z']
for band in bands:
    flux = f'flux_{band}'
    ivar = f'flux_ivar_{band}'
    mask = (df[flux] > 0) & (df[ivar] > 0)
    df = df[mask].copy()

    # Convert flux to magnitude and calculate errors
    df[f'{band}'] = 22.5 - 2.5 * np.log10(df[flux])
    df[f'{band}_err'] = 2.5 / np.log(10) * (1 / (np.sqrt(df[ivar]) * df[flux]))

# Rename RA/DEC columns
df = df.rename(columns={"ra": "RA", "dec": "DEC"})

# Filter by r-band magnitude < 22
mask = df['r'] < 22
df = df[mask]

# Save to CSV
csv_path = os.path.join(PATH, "ls.csv")
df.to_csv(csv_path, index=False)
print(f"Saved catalog to {csv_path}")