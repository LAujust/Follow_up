from astropy.table import Table
from utils import cutout_fits

tels = ['sitian','WHUT','LCO','TNOT']
size = 2000
redo = False

for tel in tels:
    cutout_fits(telescope=tel,redo=redo)