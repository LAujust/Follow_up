import numpy as np
from astropy.table import Table
from transient_hub_api import TransientHub
import json, glob
import os, sys, datetime
import logging
from pathlib import Path
import time
sys.path.append("/home/liangrd/Follow_up/code")
from utils import cutout_fits

CANDIDATES = Table.read("/home/liangrd/Follow_up/Candidates.csv")
tel_json = json.load(open("/home/liangrd/nadc_tels.json"))
ROOT = "/home/liangrd/optical_data"

hub = TransientHub()
registered_events = hub.get_events()


def setup_logger():
    log_dir = "./logs"
    os.makedirs(log_dir, exist_ok=True)  # create if not exists
    log_file = os.path.join(log_dir, f"run_{datetime.datetime.now():%Y%m%d_%H%M%S}.log")
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s"
    )

    # File handler
    fh = logging.FileHandler(log_file)
    fh.setFormatter(formatter)

    # Console handler
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger



def upload_event_files(event_id, logger, root=ROOT):
    logger.info(f"Processing {event_id}")
    #Check registraition of event
    if not event_id in registered_events["event_id"].to_list():
        meta = CANDIDATES[CANDIDATES['EP Name']==event_id][0]
        o_RA = meta["o_RA"]
        if isinstance(o_RA,float):
            ra, dec = meta["o_RA"], meta["o_Dec"]
        else:
            ra, dec = meta["RA"], meta["Dec"]
        logger.info(f"Creating Event: {event_id}: ra={ra}, dec={dec}")
        hub.create_event(event_id=event_id, ra=ra, dec=dec, priority=1)
    else:
        logger.info(f"{event_id} registered!")
        
    #Search cutouts and photometry results
    try:
        pipeline_dir = os.path.join(root,event_id,"pipeline")
        photometry = Table.read(os.path.join(pipeline_dir,"photometry.csv"))

        for tid, tel_info in tel_json.items():
            tel_name = tel_info['name']
            logger.info(f"Searching {tel_name} data...")
            cutout_dir = os.path.join(pipeline_dir,tel_name,"cutouts")
            if os.path.exists(cutout_dir):
                logger.info(f"Yes...")
                flist = glob.glob(os.path.join(cutout_dir,'**cutout.**'))
                
                #Search Photometry data
                raw_fnames = list(photometry["coadd_file"])
                cutout_fnames = [Path(rf).name.replace(".fits", "_cutout.fits") for rf in raw_fnames]
                photometry["cutout_name"] = cutout_fnames

                for f in flist:
                    f = Path(f)
                    #Check File Status
                    uploaded_files = hub.get_files(event_id)
                    uploaded_files = [Path(fi).name for fi in uploaded_files]
                    if f.name in uploaded_files:
                        logger.info(f"{f.name} already exists!")
                        continue
                    photo = photometry[photometry['cutout_name']==f.name][0]
                    lim_mag = photo['upper_limit']
                    band = photo['band']
                    photo_res = {'lim_mag':lim_mag}
                    
                    mag, mag_err = None, None
                    if isinstance(photo['magap'], float):
                        mag, mag_err = photo['magap'], photo['magap_err']
                        
                        
                    if isinstance(photo['magpsf'], float):
                        mag, mag_err = photo['magpsf'], photo['magpsf_err']
                        
                    if mag and mag_err:
                        photo_res.update({
                            'mag':mag,
                            'mag_err':mag_err
                        })
                    # else:
                    #     photo_res = {'lim_mag':lim_mag}
                    
                    logger.info(f"{f.name}: {photo_res}")
                    
                    #UPLOAD
                    logger.info(f"Uploadig {f.name}....")
                    hub.upload_image(tid, event_id, f, band, photometry_results=photo_res)
                    logger.info("Done!")
            
            
    except Exception as e:
        logger.error(f"{e}")
            
        
        
        
if __name__ == "__main__":

    events = CANDIDATES["EP Name"]
    logger = setup_logger()
    
    #Cutout Images
    tels = ['sitian','WHUT','LCO','TNOT']
    size = 2000
    redo = False

    logger.info(f"Pre-processing\n")
    logger.info('=='*30)
    for tel in tels:
        logger.info(f"Cuting out {tel} data...")
        cutout_fits(telescope=tel,redo=redo)
    
    for event_id in events:
        upload_event_files(event_id,logger)
    # tel = 'XL100'
    # event_id = 'EP251222a'
    # if event_id not in registered_events['event_id'].values:
    #     cand = CANDIDATES[CANDIDATES['EP Name']==event_id][0]
    #     ra, dec = cand['RA'], cand['Dec']
    #     print(f"Creating event {event_id} at RA={ra}, Dec={dec}")
    #     hub.create_event(event_id, ra, dec, priority=1)
    # fpaths = glob.glob(f"/home/liangrd/optical_data/{event_id}/pipeline/sitian/cutouts/*.fits")
    # for fpath in fpaths:
    #     print('=='*40)
    #     print("Uploading", fpath)
    #     tstart = time.time()
    #     hub.upload_image(tel,event_id,fpath, filter='g', photometry_results={'lim_mag':19.0})
    #     tend = time.time()
    #     print(f"Upload completed in {tend-tstart:.1f} seconds")
