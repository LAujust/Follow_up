import warnings
warnings.filterwarnings("ignore")
import os
from pathlib import Path
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.time import Time
from utils import *
from coadd import *
from photometry import *
import subprocess

FOLLOWUP_ROOT = Path(__file__).resolve().parents[1]
OPTICAL_ROOT = Path(os.environ.get("FOLLOWUP_OPTICAL_DIR", str(Path.home() / "optical_data"))).expanduser()
RESULTS_ROOT = FOLLOWUP_ROOT / "results"
check_source_dirs()


"============================================="
"                Update data                  "
"============================================="

print("===============================================")
print(f"Updating Data: {Time.now().iso}")
print("===============================================")


get_tnot_data()

get_sitian_data()

show_shift(root=str(OPTICAL_ROOT), save_dir=str(RESULTS_ROOT))
show_shift(root=str(OPTICAL_ROOT), save_dir=str(RESULTS_ROOT), ploter='plotly')
calculate_observation_stats()
show_obs_pies()
show_cumulative_observations()


#Push to GitHub
os.chdir('/home/liangrd/Follow_up')
subprocess.run(['git', 'add', '.'])
subprocess.run(['git', 'commit', '-m', f'Update data at {Time.now().iso}'])
subprocess.run(['git', 'push'])