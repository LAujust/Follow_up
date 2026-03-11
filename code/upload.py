import requests
from pathlib import Path
import json
import glob, time

tel_json = json.load(open("/home/liangrd/nadc_tels.json"))

URL = "https://nadc.china-vo.org/mwr/transient_hub/api/v1/upload"

def upload_transient_image(tid, password, event_id, fits_path, timeout=60):
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

        r = requests.post(URL, data=data, files=files, timeout=timeout)

    r.raise_for_status()

    try:
        return r.json()
    except Exception:
        return r.text
    
    

tel = 'XL100'
event_id = 'EP260116a'
fpaths = glob.glob(f"/home/liangrd/optical_data/{event_id}/sitian/cutouts/*.fits")
for fpath in fpaths:
    print('=='*40)
    print("Uploading", fpath)
    tstart = time.time()
    upload_transient_image(tel,tel_json[tel]["password"],event_id,fpath)
    tend = time.time()
    print(f"Upload completed in {tend-tstart:.1f} seconds")