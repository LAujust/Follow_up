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
tap = pyvo.dal.TAPService("https://api.skymapper.nci.org.au/public/tap")

# Use bounding box instead of cone search
ra_center, dec_center = RA, Dec
radius = err / 3600  # convert arcsec to degrees


query = f"""
SELECT TOP 4000
  raj2000, dej2000, object_id, {', '.join([f'{band}_psf, e_{band}_psf' for band in ['u','v','g','r','i','z']])}
FROM dr4.master
WHERE raj2000 BETWEEN {ra_center - radius} AND {ra_center + radius}
  AND dej2000 BETWEEN {dec_center - radius} AND {dec_center + radius}
"""

# Run query
results = tap.search(query)
df = results.to_table().to_pandas()
print(df.head())


# Rename RA/DEC columns
rename_columns = {"raj2000": "RA", "dej2000": "DEC",}
rename_columns.update({f'{band}_psf': f'{band}' for band in ['u','v','g','r','i','z']})
df = df.rename(columns=rename_columns)

# Filter by r-band magnitude < 22
mask = df['r'] < 22
df = df[mask]

# Save to CSV
csv_path = os.path.join(PATH, "skymapper.csv")
df.to_csv(csv_path, index=False)
print(f"Saved catalog to {csv_path}")