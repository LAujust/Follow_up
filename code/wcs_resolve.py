from astroquery.astrometry_net import AstrometryNet
AstrometryNet.key = 'rtdbukbqvtiienpu'
from twirl import find_peaks, gaia_radecs, compute_wcs
from twirl.geometry import sparsify
import numpy as np
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from photutils.aperture import CircularAperture
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from rich import print
import sys
import os

# define constants
RA_DEC_DEFINED = True
PEAK_NUM = 35 # number of brighter stars to be used
# 定义像素比例字典 Need to be accurate!
pixel_scale_arcsec = {
    "schmidt": 0.426,
    "LCO": 0.4,
    "GTC": 0.25,
    "TRT": 0.543,
    "Maidanak": 0.2668,
    "skob": 0.29,
    "VST": 0.213,
    "Keck_LRIS": 0.135,
    "REM": 0.57,
    "LT": 0.304,
    "CTIO": 0.526,
    "WFI": 0.237,
    "goodman": 8.0319879*0.036,
    'ultracam':0.076,
    'iTelescope':0.359, #eliot Telescope
    'SALT': 0.14, #1x1 bin
}
PIXEL_SCALE_ARCSEC = pixel_scale_arcsec['GTC']
# PIXEL_SCALE_ARCSEC = 0.222 #Megellan
# PATH = '/Users/lwx-mac/Documents/EP/EP250304a/guozhen/stack/'
# PATH = '/Users/lwx-mac/Documents/EP/EP250304a/0313-0315/'
# PATH = '/Users/lwx-mac/Documents/EP/EP250108a/Keck_yangyi_0226/'#flat_corrected/'
# PATH = '/Users/lwx-mac/Documents/EP/EP250304a/0402-0411/'
# PATH = '/Users/lwx-mac/Documents/wyn/'
# PATH = '/Users/lwx-mac/Documents/EP/EP250506a/0507/stack/'
# PATH = '/Volumes/T7/Shared_Files/EP/Results/SBO/data/EP240506a_AT2024ofs/EP240506a_TRT/'
# PATH = '/Volumes/T7/Shared_Files/EP/Results/SBO/data/EP240506a_AT2024ofs/2024ofs_warps_60418_60472/gpc2/'
# PATH = '/Volumes/T7/Shared_Files/EP/Results/SBO/data/AT2026ndp/LCO/lco_data-20260522-9/'
PATH = '/home/liangrd/optical_data/EP260131a/SOAR/asteris/'

PS_TABLE = None

REDO = False

SAVE_PATH = os.path.join(PATH,'wcs')
if not os.path.exists(SAVE_PATH):
    os.makedirs(SAVE_PATH, exist_ok=True)

# RA = 213.978000; DEC = -16.715000
# RA = 191.506952; DEC = -9.719133 #2024gsa
# RA = 28.84754 ; DEC = 5.938333 # 1021a
# RA, DEC = 59.3567569083, -46.18546 #abfo
# RA, DEC = 55.6182889176, -22.50591 #0108a
# RA, DEC = 208.3945, -42.8046 #0304a
# RA, DEC = 208.3945, -42.8046 #0304a
#RA = 35.016; DEC = +3.329 #1107a
#RA = 191.55901; DEC = -9.61731 #2024gsa
#RA = 239.876 ; DEC = 25.9202 #TCrB
# RA, DEC = 299.1500, 3.4600
# RA, DEC = 213.958859, -16.660795 #0506a
# RA, DEC = 215.2036, -10.100683 #0512a
# RA, DEC = 113.7617,	6.3244 #260116a
# RA, DEC = 46.0605, 39.438
# RA, DEC = 199.4371, 65.8476 #260110a
# RA, DEC = 178.531, 16.486 #260104a
# RA, DEC = 46.0605, 39.438 #260110a
RA, DEC = 149.9024550, -3.3074489 #260131a
# RA, DEC = 299.04067125, -48.63621125 #At2026ndp
# file to be processed and output
#filename = sys.argv[1]


pixel_scale_arcsec = np.array([PIXEL_SCALE_ARCSEC,PIXEL_SCALE_ARCSEC])
pixel_scale_degree = pixel_scale_arcsec/3600
for filename in os.listdir(PATH):
    if filename.find('.fits')>0:
    #if filename=='calibrated-T72-elioth-NGC1493-20241116-035117-R-BIN1-_-060-004.fit':
        if filename[0] == '.':
            continue
        new_filename = os.path.join(SAVE_PATH, filename.replace(".fits","_wcs.fits"))
        if not REDO:
            if os.path.exists(new_filename):
                print(f"WCS file already exists, skip: {new_filename}")
                continue
        # parse object name and acquire ra/dec coordinates
        if RA_DEC_DEFINED==False:
            object_name = filename.split('/')[-1].split('-')[1]
            print(f"Checking Simbad for [{object_name}].")
            result = Simbad.query_object(object_name)
            if result is None:
                print('Simbad query return no result.')
                exit()
            coordinates = SkyCoord(result['RA'][0]+' '+result['DEC'][0], unit='hr, deg')
            center = [coordinates.ra.value, coordinates.dec.value]
        else:
            #center = [float(sys.argv[2]), float(sys.argv[3])]
            center = [RA, DEC]
        print(f"(RA, DEC):  ({center[0]},{center[1]}).", filename)

        # read in data and looking for sources
        with fits.open(PATH + filename) as hdu:
            try:
                data = hdu[0].data
                valid_index = 0
                if data is None or data.shape[0]==0:
                    data = hdu[1].data
                    valid_index = 1
            except:
                raise ValueError(f"Cannot find data in {filename}")
            #data = np.nan_to_num(data0,nan=0.0)

            if type(data) == np.ndarray:
                try:
                    peaks = find_peaks(data)[0:PEAK_NUM]
                    print(f"{len(peaks)} peaks found in the image.")
                    

                    # checking against gaia catalogue and compute WCS
                    fov = data.shape*pixel_scale_degree
                    print(f"Checking Gaia Catalog for {PEAK_NUM} brightest sources within a FOV of {1.5*fov} deg cone.")
                    all_radecs = gaia_radecs(center, 1.5*fov)
                    all_radecs = sparsify(all_radecs, 0.01) # # we only keep stars 0.01 degree apart from each other

                    # all_radecs = PS_TABLE['RA','DEC'].as_array()
                    wcs = compute_wcs(peaks, all_radecs[0:PEAK_NUM],tolerance=3,asterism=3) #asterism: Number of sources to use for matching, 3 or 4
                    # update WCS to the new fits file
                    hdu[valid_index].header.update(wcs.to_header())
                    hdu[valid_index].writeto(new_filename, overwrite=True)
                    print(f"WCS updated, file written to: {new_filename}")
                except Exception as e:
                    print(f"Error processing {filename}: {e}")
'''
# plot out result
plt.figure(figsize=(16,12), dpi=80)
plt.imshow(data, vmin=np.median(data), vmax=3 * np.median(data), cmap="Greys_r")
radec_xy = np.array(wcs.world_to_pixel_values(all_radecs))
CircularAperture(radec_xy, 5).plot(color='r',alpha=0.5)
plt.show()
'''
