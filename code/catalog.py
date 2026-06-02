import requests
import numpy as np
import pandas as pd
from io import StringIO
import sys, os
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import pyvo

__all__ = ["PS","LS","SkyMapper"]


class PS:
    """
    Pan-STARRS DR2 catalog query & region-file generator
    """
    def __init__(self):
        self.BASE_URL = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv"

    def get_catalog(
        self,
        ra,
        dec,
        radius=0.1,
        nDetections_min=4,
        save_path=None,
        filename="ps.csv",
    ):
        """
        Query Pan-STARRS DR2 mean catalog

        Parameters
        ----------
        ra, dec : float
            Sky position in degrees (FK5)
        radius : float
            Search radius in degrees
        nDetections_min : int
            Minimum number of detections
        save_path : str or None
            If provided, save catalog to CSV
        filename : str
            Output CSV filename

        Returns
        -------
        df : pandas.DataFrame
            Processed catalog
        """

        params = {
            "ra": ra,
            "dec": dec,
            "radius": radius,
            "nDetections.gt": nDetections_min,
        }

        response = requests.get(self.BASE_URL, params=params)
        if response.status_code != 200:
            raise RuntimeError(
                f"Pan-STARRS query failed (status {response.status_code})"
            )

        df = pd.read_csv(StringIO(response.text))

        # ---- rename columns ----
        new_columns = {}
        for col in df.columns:
            if col == "raMean":
                new_col = "RA"
            elif col == "decMean":
                new_col = "DEC"
            elif col.endswith("MeanPSFMagErr"):
                new_col = col.replace("MeanPSFMagErr", "_err")
            elif col.endswith("MeanPSFMag"):
                new_col = col.replace("MeanPSFMag", "")
            else:
                new_col = col
            new_columns[col] = new_col

        df = df.rename(columns=new_columns)
        df = df.replace(-999, np.nan)

        if save_path is not None:
            full_path = os.path.join(save_path, filename)
            df.to_csv(f"{full_path}", index=False)
            print(f"Saved PS catalog to {full_path}")

        return df

    def generate_reg(
        self,
        df,
        ra_center,
        dec_center,
        save_path,
        ps_filename="ps.reg",
        wxt_filename="wxt.reg",
        ps_radius_arcsec=5,
        wxt_radius_arcmin=3,
    ):
        """
        Generate DS9 region files

        Parameters
        ----------
        df : pandas.DataFrame
            Catalog dataframe with RA, DEC columns
        ra_center, dec_center : float
            WXT source position (degrees)
        save_path : str
            Output directory
        ps_filename : str
            PS catalog region filename
        wxt_filename : str
            WXT position region filename
        ps_radius_arcsec : float
            Radius of PS source circles (arcsec)
        wxt_radius_arcmin : float
            Radius of WXT circle (arcmin)
        """

        # ---- PS catalog region ----
        with open(f"{save_path}/{ps_filename}", "w", encoding="utf-8") as fh:
            fh.write("fk5\n")
            for ra, dec in zip(df["RA"], df["DEC"]):
                fh.write(f"circle({ra},{dec},{ps_radius_arcsec}\")\n")

        # ---- WXT source region ----
        with open(f"{save_path}/{wxt_filename}", "w", encoding="utf-8") as fh:
            fh.write("fk5\n")
            fh.write(f"circle({ra_center},{dec_center},{wxt_radius_arcmin}')\n")
            



class LS:
    """
    Legacy Survey catalog query via official viewer API
    """

    def __init__(self):
        self.BASE_URL = "https://www.legacysurvey.org/viewer/cat.fits"

    def get_catalog(
        self,
        ra,
        dec,
        radius=0.1,
        save_path=None,
        filename="ls.csv",
    ):
        """
        Query Legacy Survey catalog

        Parameters
        ----------
        ra, dec : float (deg)
        radius : float (deg)
        save_path : str or None
        filename : str

        Returns
        -------
        df : pandas.DataFrame
        """

        tap = pyvo.dal.TAPService("https://datalab.noirlab.edu/tap")
        ra_center, dec_center = ra, dec

        query = f"""
        SELECT TOP 4000
        ra, dec, ls_id, flux_g, flux_ivar_g,
        flux_i, flux_ivar_i, flux_r, flux_ivar_r,
        flux_z, flux_ivar_z
        FROM ls_dr10.tractor
        WHERE ra BETWEEN {ra_center - radius} AND {ra_center + radius}
        AND dec BETWEEN {dec_center - radius} AND {dec_center + radius}
        """

        results = tap.search(query)
        df = results.to_table().to_pandas()
        print(df.head())

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

        if save_path is not None:
            # Save to CSV
            csv_path = os.path.join(save_path, "ls.csv")
            df.to_csv(csv_path, index=False)
            print(f"Saved catalog to {csv_path}")
        else:
            return df



class SkyMapper:
    """
    Legacy Survey catalog query via official viewer API
    """

    def __init__(self):
        self.BASE_URL = "https://www.legacysurvey.org/viewer/cat.fits"

    def get_catalog(
        self,
        ra,
        dec,
        radius=0.1,
        save_path=None,
        filename="ls.csv",
    ):
        """
        Query Legacy Survey catalog

        Parameters
        ----------
        ra, dec : float (deg)
        radius : float (deg)
        save_path : str or None
        filename : str

        Returns
        -------
        df : pandas.DataFrame
        """

        tap = pyvo.dal.TAPService("https://api.skymapper.nci.org.au/public/tap")
        ra_center, dec_center = ra, dec

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

        if save_path is not None:
            # Save to CSV
            csv_path = os.path.join(save_path, "skymapper.csv")
            df.to_csv(csv_path, index=False)
            print(f"Saved catalog to {csv_path}")
        else:
            return df