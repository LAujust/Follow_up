import numpy as np
from astropy.table import Table
from transient_hub_api import TransientHub
import json, glob
import os, sys


CANDIDATES = Table.read("/home/liangrd/Follow_up/Candidates.csv")
tel_json = json.load(open("/home/liangrd/nadc_tels.json"))
ROOT = "/home/liangrd/optical_data"


def upload_event_files(event_id,root=ROOT):
    pipeline_dir = os.path.join(ROOT,"pipeline")
    photometry = Table.read(os.path.join(pipeline_dir,"photometry.csv"))

    for tel, tel_info in tel_json.items():
        password = tel_info['password']
        tel_name = tel_info['name']
        
        
        if os.path.exists(os.path.join(pipeline_dir,tel_name)):
            pass
            
        
        
        
upload_event_files('1')
