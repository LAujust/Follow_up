import requests
import numpy as np
import pandas as pd
from io import StringIO
import sys, os
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

__all__ = ["PS","LS"]


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
            
            
    


import requests
import numpy as np
import pandas as pd
import os
from astropy.io import fits
from io import BytesIO

__all__ = ["LS"]


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

        params = {
            "ra": ra,
            "dec": dec,
            "radius": radius
        }

        response = requests.get(self.BASE_URL, params=params)

        if response.status_code != 200:
            raise RuntimeError(
                f"Legacy Survey query failed (status {response.status_code})"
            )

        # 读取内存中的 FITS
        print(response.text[:200])  # 打印前200字符检查响应内容
        hdul = fits.open(BytesIO(response.content))

        if len(hdul) < 2:
            print("No table extension found.")
            return pd.DataFrame()

        data = hdul[1].data

        if data is None or len(data) == 0:
            print("No objects found.")
            return pd.DataFrame()

        # 转成 pandas
        df = pd.DataFrame(np.array(data).byteswap().newbyteorder())

        # ---- 重命名列，保持和 PS 风格一致 ----
        rename_map = {
            "ra": "RA",
            "dec": "DEC",
            "mag_g": "g",
            "mag_r": "r",
            "mag_z": "z",
            "magerr_g": "g_err",
            "magerr_r": "r_err",
            "magerr_z": "z_err",
        }

        df = df.rename(columns=rename_map)

        # 只保留常用列（和 PS 类似）
        keep_cols = ["RA", "DEC",
                     "g", "g_err",
                     "r", "r_err",
                     "z", "z_err"]

        # 有些字段可能不存在（保险写法）
        keep_cols = [c for c in keep_cols if c in df.columns]

        df = df[keep_cols]

        # 替换非法值
        df = df.replace([-999, np.inf, -np.inf], np.nan)

        # 保存
        if save_path is not None:
            full_path = os.path.join(save_path, filename)
            df.to_csv(full_path, index=False)
            print(f"Saved Legacy Survey catalog to {full_path}")

        return df