import requests
from pathlib import Path
import json
import glob, time
import numpy as np, pandas as pd
from astropy.table import Table

tel_json = json.load(open("/home/liangrd/nadc_tels.json"))
BASE_URL = "https://nadc.china-vo.org/mwr/transient_hub/api/v1"
CANDIDATES = Table.read("/home/liangrd/Follow_up/Candidates.csv")


def get_events():
    """
    Get the list of transient events from transient_hub.

    Returns
    -------
    list of dict
        List of events with keys: event_id, ra, dec, discovery_time.
    """
    url = f"{BASE_URL}/events"
    r = requests.get(url)
    r.raise_for_status()
    data = r.json()
    events = data["events"]
    return pd.DataFrame(events)


def create_event(event_id, ra, dec, priority=1):
    """
    Create a new transient event in transient_hub.

    Parameters
    ----------
    event_id : str
        Transient name.
    ra : float
        Right ascension in degrees.
    dec : float
        Declination in degrees.
    priority : int
        Follow-up priority (1-5).

    Returns
    -------
    dict or str
        Server response (JSON if possible).
    """
    url = f"{BASE_URL}/events"
    data = {
        "event_id": event_id,
        "ra": ra,
        "dec": dec,
        "priority": priority
    }
    r = requests.post(url, json=data)
    r.raise_for_status()
    try:
        return r.json()
    except Exception:
        return r.text



def upload_transient_image(tid, password, event_id, fits_path, type='reduced', timeout=60):
    """
    Upload a FITS image to transient_hub.

    Parameters
    ----------
    tid : str
        Telescope TID.
    password : str
        Telescope password.
    event_id : str
        Transient name.
    fits_path : str or Path
        Path to FITS file.
    timeout : int
        Request timeout in seconds.

    Returns
    -------
    dict or str
        Server response (JSON if possible).
    """

    if type == 'reduced':
        url = f"{BASE_URL}/upload/reduced"
    else:
        url = f"{BASE_URL}/upload"
    
    fits_path = Path(fits_path)
    if not fits_path.exists():
        raise FileNotFoundError(f"{fits_path} not found")

    data = {
        "tid": tid,
        "password": password,
        "event_id": event_id
    }

    with open(fits_path, "rb") as f:
        files = {
            "file": (fits_path.name, f, "application/fits")
        }

        r = requests.post(url, data=data, files=files, timeout=timeout)

    r.raise_for_status()

    try:
        return r.json()
    except Exception:
        return r.text
    


    
events = get_events()

tel = 'XL100'
event_id = 'EP260116a'
if event_id not in events['event_id'].values:
    cand = CANDIDATES[CANDIDATES['EP Name']==event_id][0]
    ra, dec = cand['RA'], cand['Dec']
    print(f"Creating event {event_id} at RA={ra}, Dec={dec}")
    create_event(event_id, ra, dec, priority=1)
fpaths = glob.glob(f"/home/liangrd/optical_data/{event_id}/pipeline/sitian/cutouts/*.fits")
for fpath in fpaths:
    print('=='*40)
    print("Uploading", fpath)
    tstart = time.time()
    upload_transient_image(tel,tel_json[tel]["password"],event_id,fpath)
    tend = time.time()
    print(f"Upload completed in {tend-tstart:.1f} seconds")