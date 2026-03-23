import requests
from pathlib import Path
import json
import glob, time
import numpy as np, pandas as pd
from astropy.table import Table

# __all__ = ["get_events", "creat_events", "upload_transients_image", "get_files"]
__all__ = ["TransientHub"]

tel_json = json.load(open("/home/liangrd/nadc_tels.json"))
BASE_URL = "https://nadc.china-vo.org/mwr/transient_hub/api/v1"
CANDIDATES = Table.read("/home/liangrd/Follow_up/Candidates.csv")

import requests
from pathlib import Path
import json
import pandas as pd
from astropy.table import Table


class TransientHub:
    def __init__(
        self,
        base_url="https://nadc.china-vo.org/mwr/transient_hub/api/v1",
        tel_config="/home/liangrd/nadc_tels.json",
    ):
        self.base_url = base_url.rstrip("/")
        self.session = requests.Session()

        # telescope config
        with open(tel_config, "r") as f:
            self.tel_json = json.load(f)

    # -----------------------------
    # Events
    # -----------------------------
    def get_events(self):
        """Fetch all transient events"""
        url = f"{self.base_url}/events"
        r = self.session.get(url)
        r.raise_for_status()
        data = r.json()
        return pd.DataFrame(data.get("events", []))

    def create_event(self, event_id, ra, dec, priority=1):
        """Create a new transient event"""
        url = f"{self.base_url}/events"

        payload = {
            "event_id": event_id,
            "ra": ra,
            "dec": dec,
            "priority": priority
        }

        r = self.session.post(url, json=payload)
        r.raise_for_status()

        try:
            return r.json()
        except Exception:
            return r.text

    # -----------------------------
    # Upload
    # -----------------------------
    def upload_image(
        self,
        tid,
        event_id,
        fits_path,
        filter=None,
        obs_type="reduced",
        photometry_results=None,
        timeout=60
    ):
        """Upload image file

        Args:
            tid (str): Telesciope ID
            event_id (str): Event ID
            fits_path (str): File path
            filter (str, optional): Band name. Defaults to None.
            obs_type (str, optional): Data type (reduced or None). Defaults to "reduced".
            photometry_results (dict, optional): Photometry results, containing (mag, mag_err) or lim_mag keywords. Defaults to None.
            timeout (int, optional): Atempt timout seconds. Defaults to 60.
        """

        if tid not in self.tel_json:
            raise ValueError(f"{tid} not found in telescope config")

        tel_name = self.tel_json[tid]["name"]
        password = self.tel_json[tid]["password"]

        if obs_type == "reduced":
            url = f"{self.base_url}/upload/reduced"
        else:
            url = f"{self.base_url}/upload"

        fits_path = Path(fits_path)
        if not fits_path.exists():
            raise FileNotFoundError(f"{fits_path} not found")

        data = {
            "tid": tid,
            "password": password,
            "event_id": event_id,
            "obs_type": "imaging"
        }

        if filter:
            data["filter_name"] = filter

        if photometry_results:
            data.update(photometry_results)

        with fits_path.open("rb") as f:
            files = {
                "file": (fits_path.name, f, "application/fits")
            }

            r = self.session.post(
                url,
                data=data,
                files=files,
                timeout=timeout
            )

        r.raise_for_status()

        try:
            return r.json()
        except Exception:
            return r.text
        

    # -----------------------------
    # Files
    # -----------------------------
    def get_files(self, event_id):
        """Get uploaded files

        Args:
            event_id (str): Event ID

        Returns:
            list: file names
        """
        url = f"{self.base_url}/events/{event_id}/files"

        r = self.session.get(url)
        r.raise_for_status()

        data = r.json()
        files = data.get("files", [])

        return [item["key"] for item in files]

    # -----------------------------
    # Helper
    # -----------------------------



# def get_events():
#     """
#     Get the list of transient events from transient_hub.

#     Returns
#     -------
#     list of dict
#         List of events with keys: event_id, ra, dec, discovery_time.
#     """
#     url = f"{BASE_URL}/events"
#     r = requests.get(url)
#     r.raise_for_status()
#     data = r.json()
#     events = data["events"]
#     return pd.DataFrame(events)


# def create_event(event_id, ra, dec, priority=1):
#     """
#     Create a new transient event in transient_hub.

#     Parameters
#     ----------
#     event_id : str
#         Transient name.
#     ra : float
#         Right ascension in degrees.
#     dec : float
#         Declination in degrees.
#     priority : int
#         Follow-up priority (1-5).

#     Returns
#     -------
#     dict or str
#         Server response (JSON if possible).
#     """
#     url = f"{BASE_URL}/events"
#     data = {
#         "event_id": event_id,
#         "ra": ra,
#         "dec": dec,
#         "priority": priority
#     }
#     r = requests.post(url, json=data)
#     r.raise_for_status()
#     try:
#         return r.json()
#     except Exception:
#         return r.text



# def upload_transient_image(tid, password, event_id, fits_path, filter:str=None, type='reduced', photometry_results:dict=None, timeout=60):
#     """
#     Upload a FITS image to transient_hub.

#     Parameters
#     ----------
#     tid : str
#         Telescope TID.
#     password : str
#         Telescope password.
#     event_id : str
#         Transient name.
#     fits_path : str or Path
#         Path to FITS file.
#     filter: str
#         band name.
#     photometry_results: dict
#         should contain mag, mag_err or lim_mag.
#     timeout : int
#         Request timeout in seconds.

#     Returns
#     -------
#     dict or str
#         Server response (JSON if possible).
#     """

#     if type == 'reduced':
#         url = f"{BASE_URL}/upload/reduced"
#     else:
#         url = f"{BASE_URL}/upload"
    
#     fits_path = Path(fits_path)
#     if not fits_path.exists():
#         raise FileNotFoundError(f"{fits_path} not found")

#     data = {
#         "tid": tid,
#         "password": password,
#         "event_id": event_id,
#         "obs_type": "imaging"
#     }
    
#     if filter:
#         data.update({'filter_name':filter})
        
#     if photometry_results:
#         data.update(photometry_results)

#     with open(fits_path, "rb") as f:
#         files = {
#             "file": (fits_path.name, f, "application/fits")
#         }

#         r = requests.post(url, data=data, files=files, timeout=timeout)

#     r.raise_for_status()

#     try:
#         return r.json()
#     except Exception:
#         return r.text
    
    
# def get_files(event_id):
#     """get file list for an event
    
#     Args:
#         event_id (str): transient name

#     Returns:
#         list: list of files under event_id
#     """
#     url = f"{BASE_URL}/events/{event_id}/files"
#     r = requests.get(url)
#     files_dictlist =r.json()['files']
#     flist = [item['key'] for item in files_dictlist]
#     return flist
    


# get_files('EP260105a')
    
# events = get_events()

# tel = 'XL100'
# event_id = 'EP251222a'
# if event_id not in events['event_id'].values:
#     cand = CANDIDATES[CANDIDATES['EP Name']==event_id][0]
#     ra, dec = cand['RA'], cand['Dec']
#     print(f"Creating event {event_id} at RA={ra}, Dec={dec}")
#     create_event(event_id, ra, dec, priority=1)
# fpaths = glob.glob(f"/home/liangrd/optical_data/{event_id}/pipeline/sitian/cutouts/*.fits")
# for fpath in fpaths:
#     print('=='*40)
#     print("Uploading", fpath)
#     tstart = time.time()
#     upload_transient_image(tel,tel_json[tel]["password"],event_id,fpath)
#     tend = time.time()
#     print(f"Upload completed in {tend-tstart:.1f} seconds")