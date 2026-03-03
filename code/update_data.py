import warnings
warnings.filterwarnings("ignore")
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.time import Time
from utils import *
from coadd import *
from photometry import *
check_source_dirs()


"============================================="
"                Update data                  "
"============================================="

print("===============================================")
print(f"Updating Data: {Time.now().iso}")
print("===============================================")


get_tnot_data()

get_sitian_data()

show_shift(root='/home/liangrd/Follow_up/optical_data',save_dir='/home/liangrd/Follow_up/results')
show_shift(root='/home/liangrd/Follow_up/optical_data',save_dir='/home/liangrd/Follow_up/results',ploter='plotly')
calculate_observation_stats()