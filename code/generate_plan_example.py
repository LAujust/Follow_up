import warnings
warnings.filterwarnings("ignore")
from astropy.table import Table, Column, vstack
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import numpy as np
import pandas as pd
import subprocess
import sys, os, json
from ipyaladin import Aladin
from astroquery.ned import Ned
from astroquery.sdss import SDSS
from astropy import coordinates as coords
from utils import *
from catalog import *
from photometry import *
# sys.path.append('/Volumes/T7/Shared_Files/EP/Projects')
# from host_solver import HostGalaxyFinder
# from region import *
check_source_dirs()

sources = Table.read("/home/liangrd/Follow_up/Candidates.csv")

#TNOT
name = ['EP260321a']
idx = np.isin(sources['EP Name'],name)
ra, dec = list(sources['RA'][idx]), list(sources['Dec'][idx])
obs_time = list(sources['Obs Time'][idx])
significance = list(sources['Sx'][idx])

generate_tnot_plan(target=name,ra=ra,dec=dec,save_path='/home/liangrd/Follow_up/wxtsource/tnot')
generate_tnot_object_json(target=name,ra=ra,dec=dec,obs_time=obs_time,significance=significance,save_path='/home/liangrd/Follow_up/wxtsource/tnot')

#Sitian
name = ['EP260321a']
idx = np.isin(sources['EP Name'],name)
ra, dec = list(sources['RA'][idx]), list(sources['Dec'][idx])
obs_time = list(sources['Obs Time'][idx])
significance = list(sources['Sx'][idx])
generate_sitian_plan(target=name,ra=ra,dec=dec,expcount=10,exptime=180,save_path='/home/liangrd/Follow_up/wxtsource/sitian')