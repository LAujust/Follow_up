import sys, os
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from astropy.table import Table
from astropy.time import Time
import pandas as pd
import matplotlib.pyplot as plt
from astropy.cosmology import Planck18
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14

CANDIDATES_DIR = '/home/liangrd/Follow_up/Candidates.csv'
DATA_DIR = '/home/liangrd/optical_data'
BAND_LIST = ['z','i','r','g','w']
SAVE_DIR = '/home/liangrd/Follow_up/results/lcs'



MARKER_MAP = {'sitian':'o','TNOT':'s','LCO':'D','WHUT':'P','SOAR':'X'}
COLOR_LIST = ['#9e0142','#d53e4f','#f46d43','#fdae61','#ff7f00',
          '#b2df8a','#66c2a5','#66c2a5','#3288bd','#5e4fa2',
          '#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e']
idx = np.linspace(0, len(COLOR_LIST)-1, len(BAND_LIST), dtype=int)
COLOR_MAP = {BAND_LIST[i]: COLOR_LIST[b] for i, b in enumerate(idx)}
candidates = Table.read(CANDIDATES_DIR, format='csv')



def plot_photometry(data_dir, target, save_dir='./', **mpl_kwargs):
    
    
    #target info
    meta = candidates[candidates['EP Name'] == target]
    T0 = Time(meta['Obs Time'][0])
    z = meta['Redshift'][0] if not candidates['Redshift'].mask[0] else None
    
    #photometry data
    if not os.path.exists(data_dir):
        print(f"Data file {data_dir} not found.")
        return None
    
    
    data = Table.read(data_dir)
    valid_data = data[data['status']=='ok']
    if len(valid_data) == 0:
        print(f"No valid photometry data for {target}.")
        return None
    
    valid_data['dt'] = valid_data['mean_mjd'] - T0.mjd
    # print(valid_data)
    
    idx_obs = (~valid_data['magpsf'].mask) | (~valid_data['magap'].mask)
    idx_uplim = ~idx_obs
    
    obs = valid_data[idx_obs]
    uplim = valid_data[idx_uplim]
    
    #average masked column
    magpsf = np.ma.array(obs['magpsf'])
    magap  = np.ma.array(obs['magap'])
    
    magpsf_err = np.ma.array(obs['magpsf_err'])
    magap_err = np.ma.array(obs['magap_err'])
    # 计算平均（masked 会自动处理）
    mag = np.ma.mean(np.ma.vstack([magpsf, magap]), axis=0)
    mag_err = np.ma.mean(np.ma.vstack([magpsf_err, magap_err]), axis=0)

    # 只保留至少有一个值的行
    idx = ~mag.mask

    # 写入新列
    obs = obs[idx]
    obs['mag'] = mag[idx]
    obs['mag_err'] = mag_err[idx]
    
    if mpl_kwargs:
        plt.rcParams.update(mpl_kwargs)
    fig, ax = plt.subplots(figsize=(8, 5))
    
    #plot uplim
    marker = 'v'
    for band in np.unique(uplim['band']):
        ax.errorbar(uplim['dt'][uplim['band']==band], uplim['upper_limit'][uplim['band']==band],color=COLOR_MAP.get(band,'gray'), fmt=marker, label=f'{band}', alpha=0.7, markersize=10)
        
    #plot obs
    for band in np.unique(obs['band']):
        for tel in np.unique(obs['telescope']):
            idx = (obs['band']==band) & (obs['telescope']==tel)
            ax.errorbar(obs['dt'][idx], obs['mag'][idx], yerr=obs['mag_err'][idx], color=COLOR_MAP.get(band,'gray'), fmt=MARKER_MAP.get(tel,'o'), mec='k', label=f'{band} {tel}', alpha=0.7, markersize=10)
            
    ax.set_xlabel('Observer-frame Time (days)')
    ax.set_ylabel('Magnitude')
    ax.invert_yaxis()
    ax.tick_params(axis='both', which='both', direction='in')
    ax.set_title(target)
    
    if z:
        distmod = Planck18.distmod(z).value
        ymin, ymax = ax.get_ylim()
        ax2 = ax.twinx()
        ax2.set_ylim(ymin - distmod, ymax - distmod)
            
    # telescope legend
    tel_handles = [
            plt.Line2D([], [], marker=MARKER_MAP.get(tel,'grey'), linestyle='',
                       color='k', label=tel)
            for tel in np.unique(valid_data['telescope'])
        ]

    # band legend
    band_handles = [
            plt.Line2D([], [], marker='v', linestyle='',
                       color=COLOR_MAP[band], label=band)
            for band in np.unique(valid_data['band'])
        ]

    handles = tel_handles + band_handles

    plt.legend(
        handles=handles,
        ncol=2,
        # bbox_to_anchor=(0.8, 1),
        loc='lower right'
    )
    
    fname = f"{target}_lc.png"
    plt.savefig(os.path.join(save_dir,fname), bbox_inches='tight',dpi=300)
    print(f"Saved photometry plot to {os.path.join(save_dir,fname)}")
    plt.close()
    
    

def main():
    targets = candidates['EP Name'].tolist()
    if not os.path.exists(SAVE_DIR):
        os.makedirs(SAVE_DIR, exist_ok=True)
        
    for target in targets:
        data_dir = os.path.join(DATA_DIR, target, 'pipeline', 'photometry.csv')
        print('=='*20)
        print(f"Plotting photometry for {target}...")
        plot_photometry(data_dir, target, SAVE_DIR)
        

if __name__ == "__main__":
    main()
    # plot_photometry('/home/liangrd/optical_data/EP260101a/pipeline/photometry.csv','EP260101a', SAVE_DIR)