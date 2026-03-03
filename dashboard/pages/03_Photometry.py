from __future__ import annotations

from pathlib import Path
import streamlit as st

from config import PHOTOMETRY_OUTPUT_DIR
from data_access import build_fits_index, load_candidates

try:
    from services.photometry_service import PhotometryRequest, run_photometry
except Exception as exc:  # pragma: no cover
    st.error(
        "Photometry tools are unavailable because required dependencies are missing "
        f"(details: {exc}). Install requirements.txt and redeploy."
    )
    st.stop()

st.title("Photometry")

cand = load_candidates()
fits_idx = build_fits_index()

if cand.empty:
    st.error("Candidates.csv is missing or empty.")
    st.stop()
if fits_idx.empty:
    st.warning("No FITS files found under optical_data/.")
    st.stop()

left, right = st.columns(2)
with left:
    target = st.selectbox("Target", sorted(cand["target"].unique().tolist()))
with right:
    method = st.selectbox("Method", ["psf", "aperture"])

row = cand[cand["target"] == target].iloc[0]
def_ra = float(row["RA"])
def_dec = float(row["Dec"])

ra = st.number_input("RA (deg)", value=def_ra, format="%.6f")
dec = st.number_input("Dec (deg)", value=def_dec, format="%.6f")

subset = fits_idx[fits_idx["target"] == target].sort_values("file")
if subset.empty:
    st.warning("No FITS files found for selected target.")
    st.stop()

selected_file = st.selectbox("FITS file", subset["path"].tolist())

p1, p2, p3, p4 = st.columns(4)
fwhm = p1.number_input("FWHM", min_value=1.0, max_value=20.0, value=3.0)
sigma = p2.number_input("Sigma", min_value=1.0, max_value=20.0, value=5.0)
match_radius = p3.number_input("Match radius (arcsec)", min_value=0.1, max_value=10.0, value=1.0)
forced = p4.checkbox("Forced mode (PSF)", value=True)

if st.button("Run photometry"):
    req = PhotometryRequest(
        target=target,
        fits_file=selected_file,
        ra=float(ra),
        dec=float(dec),
        method=method,
        fwhm=float(fwhm),
        sigma=float(sigma),
        match_radius=float(match_radius),
        forced=bool(forced),
    )

    with st.spinner("Running photometry..."):
        try:
            res = run_photometry(req, Path(PHOTOMETRY_OUTPUT_DIR))
        except Exception as exc:
            st.exception(exc)
            st.stop()

    st.success(f"Completed. Output: {res['output_dir']}")

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Mag", f"{res['mag']:.4f}" if res["mag"] is not None else "-")
    c2.metric("Mag err", f"{res['mag_err']:.4f}" if res["mag_err"] is not None else "-")
    c3.metric("ZP", f"{res['zp']:.4f}" if res["zp"] is not None else "-")
    c4.metric("Upper limit", f"{res['upper_limit']:.4f}" if res["upper_limit"] is not None else "-")

    if res.get("raw_table_preview") is not None:
        st.subheader("Result table preview")
        st.dataframe(res["raw_table_preview"], use_container_width=True, height=280)

    if res.get("diagnostic_paths"):
        st.subheader("Diagnostics")
        for p in res["diagnostic_paths"]:
            st.image(p, caption=p, use_container_width=True)
