from utils import *
from catalog import *
import numpy as np
import pandas as pd
from astropy.table import Table, vstack, QTable
from scipy import optimize,stats
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
from photutils.detection import DAOStarFinder, IRAFStarFinder
from photutils.aperture import (
    SkyCircularAperture,
    SkyCircularAnnulus,
    CircularAperture,
    CircularAnnulus,
    ApertureStats,
    aperture_photometry,
)
from photutils.psf import CircularGaussianPRF, PSFPhotometry, IterativePSFPhotometry, fit_fwhm
from photutils.background import MMMBackground, Background2D, MedianBackground, LocalBackground
import matplotlib.pyplot as plt
from astropy.nddata import Cutout2D

__all__ = ['Photometry']

class Photometry:
    """
    Unified photometry class:
    - source detection
    - catalog matching
    - zeropoint calibration
    - aperture photometry
    - PSF photometry
    """

    def __init__(self,cat_dir=None,cat='PS',path='./',mag_col='r'):
        self.cat_dir = cat_dir
        self.cat = cat
        self.path = path
        self.mag_col = mag_col

        self.data = None
        self.wcs = None
        self.header = None
        self.fhwm = None
        self.zp = None
        self.zp_std = None
        self.mask = None
        
    def read_image(self, image_path, masklow=None):
        """Read FITS image and WCS

        Args:
            image_path (str): Path of FITS image
        """

        with fits.open(image_path) as hdul:
            wcs_found = False
            for i in range(len(hdul)):
                try:
                    wcs = WCS(hdul[i].header)
                    # 检查是否有天球坐标
                    if wcs.has_celestial:
                        self.wcs = wcs
                        self.header = hdul[i].header
                        wcs_found = True
                        break
                except Exception:
                    continue

            if not wcs_found:
                self.wcs = None
                raise ValueError("No valid celestial WCS found in FITS header")
                
                
            try:
                self.data = hdul[0].data.astype(np.float32)
            except:
                try:
                    self.data = hdul[1].data.astype(np.float32)
                except:
                    raise ValueError("No valid image data found in FITS file")

            self.gain = hdul[0].header.get("GAIN", 1.0)
            
        'Check'
        if self.data is not None and self.wcs is not None:
            print('Image and WCS loaded successfully.')
            'Find essential information from header'
            if masklow:
                mask = np.zeros(self.data.shape, dtype=bool)
                mask[self.data < masklow] = True
                self.mask = mask
            #Subtract median background
            mean, median, std = sigma_clipped_stats(self.data)
            self.data -= median
            self.exptime = self.header.get('EXPTIME', None)
            self.band = self.header.get('FILTER', None)
            self.center = self._get_image_center_coord()
            self.error, self.bkg = self.estimate_global_error(self.data)
            self.psf_model = CircularGaussianPRF(flux=1.,fwhm=self.fhwm if self.fhwm is not None else 3.0)

            
            
    def estimate_global_error(self, data, box_size=256, filter_size=3):
        """
        使用 photutils Background2D 估计像素级误差图
        """
        sigma_clip = SigmaClip(sigma=3.0, maxiters=5)
        bkg_stat = MMMBackground(sigma_clip=sigma_clip)
        bkg = Background2D(self.data, self.data.shape, filter_size=(11, 11),
                   bkg_estimator=bkg_stat)

        error = bkg.background_rms
        
        return error, bkg
            

        
    def detect_sources(self, fwhm=3.0, threshold_sigma=5.0, finder=IRAFStarFinder):
        mean, median, std = sigma_clipped_stats(self.data)

        finder = finder(
            fwhm=fwhm,
            threshold=threshold_sigma * std,
        )

        sources = finder(self.data-median, mask=self.mask)
        
        #fit FWHM
        xypos = list(zip(sources['xcentroid'], sources['ycentroid']))
        fwhm = fit_fwhm(self.data, xypos=xypos, error=self.error, fit_shape=(5, 5), fwhm=3)
        self.sources = sources
        return self.sources
    
    def _estimate_fwhm(self):
        if self.sources is None:
            raise RuntimeError("No sources detected yet")
        else:
            mean, median, std = sigma_clipped_stats(self.sources['fwhm'])
            self.fhwm = median
            print(f"Estimated FWHM: {self.fhwm:.1f}+/-{std:.1f} pixels")

    def match_catalog(
        self,
        sources=None,
        match_radius=3.0, #arcsec
        mag_col="r",
    ):
        """
        Detect point sources in image and match to Pan-STARRS catalog

        Parameters
        ----------
        ps_df : pandas.DataFrame
            PS catalog with RA, DEC, and magnitude
        fwhm : float
            FWHM in pixels
        threshold_sigma : float
            Detection threshold in sigma
        match_radius : float
            Matching radius in arcsec
        mag_col : str
            Magnitude column name in PS catalog

        Returns
        -------
        matched_table : astropy.table.Table
        """
        
        if self.cat_dir is None:
            if self.cat == 'PS':
                ps_client = PS()
            elif self.cat == 'LS':
                ps_client = LS()
            ps_df = ps_client.get_catalog(ra=self.center.ra.deg, dec=self.center.dec.deg, radius=0.5,save_path=self.path)
            self.ps_df = ps_df
        
        else:
            self.ps_df = pd.read_csv(self.cat_dir)

        if self.wcs is None:
            raise ValueError("WCS is required for catalog matching")

        # --- 3. pixel → sky
        ra, dec = self.wcs.pixel_to_world_values(
            sources["x_fit"],
            sources["y_fit"]
        )

        det_coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)

        # --- 4. PS catalog → SkyCoord
        ps_coord = SkyCoord(
            ra=self.ps_df["RA"].values * u.deg,
            dec=self.ps_df["DEC"].values * u.deg,
        )

        # --- 5. sky matching
        idx, sep2d, _ = det_coord.match_to_catalog_sky(ps_coord)

        sep_mask = sep2d < match_radius * u.arcsec

        # --- 6. build matched table
        matched = Table()
        matched["x"] = sources["x_fit"][sep_mask]
        matched["y"] = sources["y_fit"][sep_mask]
        matched["ra"] = ra[sep_mask]
        matched["dec"] = dec[sep_mask]
        matched["flux"] = sources["flux_fit"][sep_mask]

        matched["ps_ra"] = self.ps_df.iloc[idx[sep_mask]]["RA"].values
        matched["ps_dec"] = self.ps_df.iloc[idx[sep_mask]]["DEC"].values
        matched["ps_mag"] = self.ps_df.iloc[idx[sep_mask]][mag_col].values
        matched["sep_arcsec"] = sep2d[sep_mask].to(u.arcsec).value
        matched['zp_initial'] = matched["ps_mag"] + 2.5 * np.log10(matched["flux"])

        return matched

        
        
    def aperture_photometry(
            self,
            ra,
            dec,
            r_ap=None,
            r_in=None,
            r_out=None,
            sigma=5.0,
            sigma_clip_val=2.0,
            match_radius=1.0,
            mag_col="r",
            sat_level=50000,
            psf_method=PSFPhotometry,
            finder=IRAFStarFinder,
            psf_model=None,
            fwhm=3.,
            forced=False,
            fit_shape=(5, 5),
            cat_dir=None,
            plot_scale=30,
            path="./",
            show=True,
        ):
        """
        Aperture photometry with internal zeropoint calibration.

        Includes:
        - saturation filtering
        - sigma clipping
        - forced photometry option
        """

        if psf_model is None:
            if fwhm is None:
                fwhm = self.fhwm if self.fhwm is not None else 3.0
            psf_model = CircularGaussianPRF(flux=1.0, fwhm=fwhm)

        ap_radius = 1.25 * fwhm
        finder = finder(sigma, fwhm)
        bkgstat = MMMBackground()
        localbkg_estimator = LocalBackground(3*ap_radius, 6*ap_radius, bkgstat)

        psfphot = psf_method(
            psf_model=psf_model,
            aperture_radius=ap_radius,
            finder=finder,
            fit_shape=fit_shape,
            localbkg_estimator=localbkg_estimator,
        )

        # ==========================================================
        # 1️⃣ FULL FIELD PSF PHOTOMETRY
        # ==========================================================
        phot = psfphot(self.data, error=self.error)
        phot = phot[phot['flux_fit']>0.]
        self.sources = phot

        flux = phot["flux_fit"]
        flux_err = phot["flux_err"]

        # ⭐ 计算 peak 近似值用于饱和判断
        peak_est = flux / (2 * np.pi * (fwhm / 2.355)**2)

        # ==========================================================
        # 2️⃣ CATALOG MATCH
        # ==========================================================
        if cat_dir is not None:
            cat_df = pd.read_csv(cat_dir)
        else:
            ps_client = PS()
            cat_df = ps_client.get_catalog(
                ra=self.center.ra.deg,
                dec=self.center.dec.deg,
                radius=0.5,
                save_path=path,
            )

        det_coord = SkyCoord(
            *self.wcs.pixel_to_world_values(
                phot["x_fit"], phot["y_fit"]
            ),
            unit="deg"
        )

        cat_coord = SkyCoord(
            cat_df["RA"].values * u.deg,
            cat_df["DEC"].values * u.deg,
        )

        idx, sep2d, _ = det_coord.match_to_catalog_sky(cat_coord)

        match_mask = sep2d < match_radius * u.arcsec
        
        self.matched = cat_df.iloc[idx[match_mask]]

        matched_flux = flux[match_mask]
        matched_mag = cat_df.iloc[idx[match_mask]][mag_col].values
        matched_peak = peak_est[match_mask]

        # ==========================================================
        # 3️⃣ 饱和 & 合理性过滤
        # ==========================================================
        good = (
            (matched_flux > 0) &
            np.isfinite(matched_mag) &
            (matched_peak < sat_level)
        )

        print(f"Matched stars before filtering: {len(matched_flux)}")
        print(f"After saturation filtering: {np.sum(good)}")

        matched_mag_good = matched_mag[good]
        zp_values = matched_mag[good] + 2.5 * np.log10(matched_flux[good])

        if len(zp_values) < 5:
            raise RuntimeError("Too few stars after saturation filtering")

        # ==========================================================
        # 4️⃣ Sigma clipping
        # ==========================================================
        zp_mean, zp_med, zp_std = sigma_clipped_stats(
            zp_values,
            sigma=sigma_clip_val
        )

        self.zp = zp_med
        self.zp_std = zp_std
        
        zp_idx = np.where((zp_values > zp_med - sigma_clip_val * zp_std) &
                          (zp_values < zp_med + sigma_clip_val * zp_std))[0]
        zp_values = zp_values[zp_idx]
        matched_mag_good = matched_mag_good[zp_idx]

        print(f"\nZeropoint = {zp_med:.3f} ± {zp_std:.3f}  (N={len(zp_values)})")
        
        # ==========================================================
        # 3️⃣  DIAGNOSTIC PLOT
        # ==========================================================
        fig = plt.figure(figsize=(9,5),dpi=300)
        axes = fig.add_gridspec(1, 2, wspace=0.05, width_ratios=[2, 1])
        ax = axes.subplots()
        # fig, ax = plt.subplots(1, 2, figsize=(8, 4))

        # --- mag vs zp ---
        ax[0].scatter(matched_mag_good, zp_values, s=20, alpha=0.6)
        ax[0].axhline(zp_med, color="r", ls="--")
        ax[0].set_xlabel("Catalog mag")
        ax[0].set_ylabel("ZP per star")
        ax[0].set_title("mag vs ZP")

        # --- histogram ---
        ax[1].hist(zp_values, bins=20, density=True, histtype="step")
        ax[1].axvline(zp_med, color="r", ls="--")
        ax[1].set_xlabel("ZP")
        ax[1].set_title("ZP distribution")

        # --- KDE PDF ---
        kde = stats.gaussian_kde(zp_values)
        xgrid = np.linspace(zp_med - 3*zp_std, zp_med + 3*zp_std, 200)
        ax[1].plot(xgrid, kde(xgrid))
        ax12 = ax[1].twinx()
        ax12.set_ylabel("PDF")
        
        ax[0].tick_params(axis='both', which='both', direction='in')
        ax[1].tick_params(axis='both', which='both', direction='in')
        ax12.tick_params(axis='both', which='both', direction='in')

        plt.tight_layout()
        plt.savefig(f"{path}/zeropoint_diagnostic.png", bbox_inches='tight', dpi=300)
        if show:
            plt.show()
        
        # ==========================================================
        # 3️⃣  PLOT RESIDUAL IMAGE
        # ==========================================================
        pos = SkyCoord(ra*u.deg, dec*u.deg)
        size = int(plot_scale * fwhm)
        cutout = Cutout2D(self.data, pos, size=(size,size), wcs=self.wcs)
        data = cutout.data
        resid = psfphot.make_residual_image(data)
        model = psfphot.make_model_image(shape=data.shape)
        self.model = model
        self.resid = resid
        
        cmap = 'viridis'
        fig, ax = plt.subplots(1,3, figsize=(15,5),sharex=True, sharey=True)
        im0 = ax[0].imshow(data, origin='lower', cmap=cmap, vmin=np.percentile(data,5), vmax=np.percentile(data,99))
        im1 = ax[1].imshow(model, origin='lower', cmap=cmap, vmin=np.percentile(data,5), vmax=np.percentile(data,99))
        im2 = ax[2].imshow(resid, origin='lower', cmap=cmap, vmin=np.percentile(resid,5), vmax=np.percentile(resid,99))
        ax[0].set_title('Original Image')
        ax[1].set_title('PSF Model Image')
        ax[2].set_title('Residual Image')

        # if ra is not None and dec is not None:
        #     #Cutout around target
        #     x, y = self._radec_to_xy(ra, dec)
        #     size = int(plot_scale * fwhm)
        #     ax[0].set_xlim(x - size//2, x + size//2)
        #     ax[0].set_ylim(y - size//2, y + size//2)
        
        ax[0].set_xticklabels([])
        ax[0].set_yticklabels([])
        plt.tight_layout()
        plt.savefig(f"{path}/psf_residual.png", bbox_inches='tight', dpi=300)
        if show:
            plt.show()
        
        #estimage upper-limit
        self.uplim = self.estimate_upperlimit()
        
        # =========================================================
        # 4️⃣  Aperture Setup
        # =========================================================
        if r_ap is None:
            r_ap = 1.25 * fwhm
        if r_in is None:
            r_in = 2.0 * r_ap
        if r_out is None:
            r_out = 3.0 * r_ap
            
        
        position = SkyCoord(ra*u.deg, dec*u.deg)
        aperture = SkyCircularAperture(position,r=r_ap*u.arcsec)
        pix_ap = aperture.to_pixel(self.wcs)
        annulus = SkyCircularAnnulus(position,r_in=r_in*u.arcsec,r_out=r_out*u.arcsec)
        pix_ann = annulus.to_pixel(self.wcs)
        
        aperstats = ApertureStats(self.data, pix_ann)
        bkg_mean = aperstats.mean
        total_bkg = bkg_mean * pix_ap.area
        
        phot_table = aperture_photometry(self.data, pix_ap,error=self.error)
        phot_bkgsub = phot_table['aperture_sum'][0] - total_bkg
        phot_table['total_bkg'] = total_bkg
        phot_table['aperture_sum_bkgsub'] = phot_bkgsub

        
        target_coord = SkyCoord(ra*u.deg, dec*u.deg)
        sep_target = target_coord.separation(det_coord)

        if np.min(sep_target) < match_radius * u.arcsec:
            #convert to magnitude
            phot['mag'] = mag = self.zp - 2.5 * np.log10(phot_table['aperture_sum_bkgsub'][0])
            phot['mag_err'] = mag_err = np.sqrt(
                (1.0857 * phot_table['aperture_sum_err'][0] / phot_bkgsub)**2 +
                self.zp_std**2
            )
            self.phot = phot_table
            
            print(f"Aperture magnitude = {mag:.3f} ± {mag_err:.3f}" )
            
            return phot_table     

        else:
            if forced:
                if phot_bkgsub > 3.0 * total_bkg:
                    #convert to magnitude
                    phot['mag'] = mag = self.zp - 2.5 * np.log10(phot_table['aperture_sum_bkgsub'][0])
                    phot['mag_err'] = mag_err = np.sqrt(
                        (1.0857 * phot_table['aperture_sum_err'][0] / phot_bkgsub)**2 +
                        self.zp_std**2
                    )
                    self.phot = phot_table
                    
                    print(f"Aperture magnitude = {mag:.3f} ± {mag_err:.3f}" )
                    
                    return phot_table  

            print("Target not detected → computing upper limit")
            print(f"3-sigma upper limit = {self.uplim:.3f}")
            return {"upper_limit": self.uplim}




    

    def psf_photometry(
        self,
        ra,
        dec,
        psf_model=None,
        fwhm=None,
        sigma=5.0,
        psf_method=PSFPhotometry,
        finder=IRAFStarFinder,
        forced=False,
        fit_shape=(5, 5),
        match_radius=1.0,
        mag_col="r",
        sigma_clip_val=3.0,
        sat_level=50000,        # ⭐ 新增：饱和阈值（ADU）
        cat_dir=None,
        plot_scale=30, #*fwhm
        path="./",
        show=True,
    ):
        """
        PSF photometry with internal zeropoint calibration.

        Includes:
        - saturation filtering
        - flux filtering
        - sigma clipping
        """

        if psf_model is None:
            if fwhm is None:
                fwhm = self.fhwm if self.fhwm is not None else 3.0
            psf_model = CircularGaussianPRF(flux=1.0, fwhm=fwhm)

        ap_radius = 1.25 * fwhm
        finder = finder(sigma, fwhm)
        bkgstat = MMMBackground()
        localbkg_estimator = LocalBackground(3*ap_radius, 6*ap_radius, bkgstat)

        psfphot = psf_method(
            psf_model=psf_model,
            aperture_radius=ap_radius,
            finder=finder,
            fit_shape=fit_shape,
            localbkg_estimator=localbkg_estimator,
        )

        # ==========================================================
        # 1️⃣ FULL FIELD PSF PHOTOMETRY
        # ==========================================================
        phot = psfphot(self.data, error=self.error)
        phot = phot[phot['flux_fit']>0.]
        self.sources = phot

        flux = phot["flux_fit"]
        flux_err = phot["flux_err"]

        # ⭐ 计算 peak 近似值用于饱和判断
        peak_est = flux / (2 * np.pi * (fwhm / 2.355)**2)

        # ==========================================================
        # 2️⃣ CATALOG MATCH
        # ==========================================================
        if cat_dir is not None:
            cat_df = pd.read_csv(cat_dir)
        else:
            ps_client = PS()
            cat_df = ps_client.get_catalog(
                ra=self.center.ra.deg,
                dec=self.center.dec.deg,
                radius=0.5,
                save_path=path,
            )
            
        self.cat_df = cat_df

        det_coord = SkyCoord(
            *self.wcs.pixel_to_world_values(
                phot["x_fit"], phot["y_fit"]
            ),
            unit="deg"
        )

        cat_coord = SkyCoord(
            cat_df["RA"].values * u.deg,
            cat_df["DEC"].values * u.deg,
        )

        idx, sep2d, _ = det_coord.match_to_catalog_sky(cat_coord)

        match_mask = sep2d < match_radius * u.arcsec
        
        self.matched = cat_df.iloc[idx[match_mask]]

        matched_flux = flux[match_mask]
        matched_mag = cat_df.iloc[idx[match_mask]][mag_col].values
        matched_peak = peak_est[match_mask]

        # ==========================================================
        # 3️⃣ 饱和 & 合理性过滤
        # ==========================================================
        good = (
            (matched_flux > 0) &
            np.isfinite(matched_mag) &
            (matched_peak < sat_level)
        )

        print(f"Matched stars before filtering: {len(matched_flux)}")
        print(f"After saturation filtering: {np.sum(good)}")

        matched_mag_good = matched_mag[good]
        zp_values = matched_mag[good] + 2.5 * np.log10(matched_flux[good])

        if len(zp_values) < 5:
            raise RuntimeError("Too few stars after saturation filtering")

        # ==========================================================
        # 4️⃣ Sigma clipping
        # ==========================================================
        zp_mean, zp_med, zp_std = sigma_clipped_stats(
            zp_values,
            sigma=sigma_clip_val
        )

        self.zp = zp_med
        self.zp_std = zp_std
        
        zp_idx = np.where((zp_values > zp_med - sigma_clip_val * zp_std) &
                          (zp_values < zp_med + sigma_clip_val * zp_std))[0]
        zp_values = zp_values[zp_idx]
        matched_mag_good = matched_mag_good[zp_idx]

        print(f"\nZeropoint = {zp_med:.3f} ± {zp_std:.3f}  (N={len(zp_values)})")
        
        # ==========================================================
        # 3️⃣  DIAGNOSTIC PLOT
        # ==========================================================
        fig = plt.figure(figsize=(9,5),dpi=300)
        axes = fig.add_gridspec(1, 2, wspace=0.05, width_ratios=[2, 1])
        ax = axes.subplots()
        # fig, ax = plt.subplots(1, 2, figsize=(8, 4))

        # --- mag vs zp ---
        ax[0].scatter(matched_mag_good, zp_values, s=20, alpha=0.6)
        ax[0].axhline(zp_med, color="r", ls="--")
        ax[0].set_xlabel("Catalog mag")
        ax[0].set_ylabel("ZP per star")
        ax[0].set_title("mag vs ZP")

        # --- histogram ---
        ax[1].hist(zp_values, bins=10, density=True, histtype="step")
        ax[1].axvline(zp_med, color="r", ls="--")
        ax[1].set_xlabel("ZP")
        ax[1].set_title("ZP distribution")

        # --- KDE PDF ---
        kde = stats.gaussian_kde(zp_values)
        xgrid = np.linspace(zp_med - 3*zp_std, zp_med + 3*zp_std, 200)
        ax[1].plot(xgrid, kde(xgrid))
        ax12 = ax[1].twinx()
        ax12.set_ylabel("PDF")
        
        ax[0].tick_params(axis='both', which='both', direction='in')
        ax[1].tick_params(axis='both', which='both', direction='in')
        ax12.tick_params(axis='both', which='both', direction='in')

        plt.tight_layout()
        plt.savefig(f"{path}/zeropoint_diagnostic.png", bbox_inches='tight', dpi=300)
        if show:
            plt.show()
        
        # ==========================================================
        # 3️⃣  PLOT RESIDUAL IMAGE
        # ==========================================================
        #zoom into target ra, dec
        pos = SkyCoord(ra*u.deg, dec*u.deg)
        size = int(plot_scale * fwhm)
        cutout = Cutout2D(self.data, pos, size=(size,size), wcs=self.wcs)
        data = cutout.data
        resid = psfphot.make_residual_image(data)
        model = psfphot.make_model_image(shape=data.shape)
        self.model = model
        self.resid = resid
        
        cmap = 'viridis'
        fig, ax = plt.subplots(1,3, figsize=(15,5),sharex=True, sharey=True)
        im0 = ax[0].imshow(data, origin='lower', cmap=cmap, vmin=np.percentile(data,5), vmax=np.percentile(data,99))
        im1 = ax[1].imshow(model, origin='lower', cmap=cmap, vmin=np.percentile(data,5), vmax=np.percentile(data,99))
        im2 = ax[2].imshow(resid, origin='lower', cmap=cmap, vmin=np.percentile(resid,5), vmax=np.percentile(resid,99))
        ax[0].set_title('Original Image')
        ax[1].set_title('PSF Model Image')
        ax[2].set_title('Residual Image')
        # if ra is not None and dec is not None:
        #     #Cutout around target
        #     x, y = self._radec_to_xy(ra, dec)
        #     size = int(plot_scale * fwhm)
        #     ax[0].set_xlim(x - size//2, x + size//2)
        #     ax[0].set_ylim(y - size//2, y + size//2)
        
        ax[0].set_xticklabels([])
        ax[0].set_yticklabels([])
        plt.tight_layout()
        plt.savefig(f"{path}/psf_residual.png", bbox_inches='tight', dpi=300)
        if show:
            plt.show()

        #estimage upper-limit
        self.uplim = self.estimate_upperlimit()

        # ==========================================================
        # 5️⃣ TARGET MEASUREMENT
        # ==========================================================
        target_coord = SkyCoord(ra*u.deg, dec*u.deg)
        sep_target = target_coord.separation(det_coord)

        if np.min(sep_target) < match_radius * u.arcsec:

            idx_target = np.argmin(sep_target)

            target_flux = flux[idx_target]
            target_flux_err = flux_err[idx_target]

            if target_flux > 0:
                mag = self.zp - 2.5 * np.log10(target_flux)
                mag_err = np.sqrt(
                    (1.0857 * target_flux_err / target_flux)**2 +
                    self.zp_std**2
                )

                print(f"Target mag = {mag:.3f} ± {mag_err:.3f}")

                phot_target = QTable()
                phot_target["flux"] = [target_flux]
                phot_target["flux_err"] = [target_flux_err]
                phot_target["mag"] = [mag]
                phot_target["mag_err"] = [mag_err]

                return phot_target
            
        else:
            if forced:
                try:
                    x,y = self._radec_to_xy(ra, dec)
                    init_params = QTable()
                    init_params['x'] = [x]
                    init_params['y'] = [y]
                    psfphot = psf_method(
                        psf_model=psf_model,
                        aperture_radius=ap_radius,
                        # finder=finder,
                        fit_shape=fit_shape,
                        localbkg_estimator=localbkg_estimator,
                    )
                    psf_model.x_0.fixed = True
                    psf_model.y_0.fixed = True
                    # psf_model.fwhm.fixed = False

                    phot_target = psfphot(self.data, error=self.error,init_params=init_params)
                    
                    target_flux = phot_target['flux_fit'][0]   
                    target_flux_err = phot_target['flux_fit_err'][0]
                    if target_flux > 0:
                        mag = self.zp - 2.5 * np.log10(target_flux)
                        mag_err = np.sqrt(
                            (1.0857 * target_flux_err / target_flux)**2 +
                            self.zp_std**2
                        )

                        print(f"Target mag = {mag:.3f} ± {mag_err:.3f}")
                        return phot_target
                except:
                    pass

        print("Target not detected → computing upper limit")
        print(f"3-sigma upper limit = {self.uplim:.3f}")
        return {"upper_limit": self.uplim}
    
    def _radec_to_xy(self, ra, dec):
        """
        Convert RA/Dec to pixel coordinates
        """
        if self.wcs is None:
            raise ValueError("WCS is required to use RA/Dec coordinates")

        c = SkyCoord(ra * u.deg, dec * u.deg)
        x, y = self.wcs.world_to_pixel(c)
        return x, y
    
    def _xy_to_radec(self, x, y):
        """
        Convert pixel coordinates to RA/Dec
        """
        if self.wcs is None:
            raise ValueError("WCS is required to use pixel coordinates")

        c = self.wcs.pixel_to_world_values(x, y)
        return c
    
    def _get_image_center_coord(self):
        """
        Get sky coordinates (RA, Dec) of the image center.

        Returns
        -------
        center_coord : astropy.coordinates.SkyCoord
        """

        if self.wcs is None:
            raise ValueError("WCS is required to compute sky coordinates")

        ny, nx = self.data.shape

        # FITS / photutils 使用 0-based pixel convention
        x_center = (nx - 1) / 2.0
        y_center = (ny - 1) / 2.0

        ra, dec = self.wcs.pixel_to_world_values(x_center, y_center)

        return SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
    
    def estimate_upperlimit(self,snr=3,r=5,data=None,bkg=None,zp=None):
        if data is None:
            data = self.data
        if bkg is None:
            bkg = self.bkg
        if zp is None:
            zp = self.zp
            
        Npix = np.pi * r**2 #pix^2
        noise = bkg.background_rms_median * np.sqrt(Npix)
        uplim = zp - 2.5*np.log10(snr*noise)
        return uplim
    
    
    def _plot_zeropoint(self, ):
        """
        Plot zeropoint scatter and distribution.

        Parameters:
        - zp_individual: List of zeropoint values for each source.
        """
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        zp_res = self.zp_bulk
        zp_values = zp_res['zp_values']
        pesu_x = zp_res['mag']
        # Scatter plot
        axes[0].scatter(pesu_x, zp_values, label='ZP values')
        axes[0].fill_between(np.linspace(self.mag_cut,24,100), self.zp - self.zp_std, self.zp + self.zp_std, color='gray', alpha=0.3, label='ZP ± 1σ')
        axes[0].axhline(self.zp, color='k', linestyle='--', label='ZP median')
        axes[0].set_title("Zeropoint Scatter")
        axes[0].set_xlabel("AB mag")
        axes[0].set_ylabel("Zeropoint (mag)")
        axes[0].grid()

        # Histogram
        axes[1].hist(zp_values, bins=10, color='orange', alpha=0.7, density=True)
        zp_x = np.linspace(self.zp - 3*self.zp_std, self.zp + 3*self.zp_std, 100)
        pdf = stats.norm.pdf(zp_x, loc=self.zp, scale=self.zp_std)
        axes[1].plot(zp_x, pdf, color='blue', label='PDF')
        axes[1].set_title("Zeropoint Distribution")
        axes[1].set_xlabel("Zeropoint (mag)")
        axes[1].set_ylabel("Frequency")
        axes[1].grid()

        plt.tight_layout()
        plt.show()