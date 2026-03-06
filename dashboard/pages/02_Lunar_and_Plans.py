from __future__ import annotations

import pandas as pd
import streamlit as st

from config import GENERATED_PLAN_DIR
from data_access import load_candidates
from services.plan_service import generate_sitian_plan_text, generate_tnot_plan_text, save_plan_outputs

try:
    from services.lunar_service import compute_lunar_curve
except Exception as exc:  # pragma: no cover
    st.error(
        "Lunar tools are unavailable because required dependencies are missing "
        f"(details: {exc}). Install requirements.txt and redeploy."
    )
    st.stop()

st.title("Lunar Distance & Observation Plans")

cand = load_candidates()
if cand.empty:
    st.error("Candidates.csv is missing or empty.")
    st.stop()

st.header("Lunar Distance")

m1, m2 = st.columns([2, 1])
with m1:
    targets = cand["target"].tolist()
    selected_target = st.selectbox("Target", targets)
with m2:
    manual = st.checkbox("Manual RA/Dec", value=False)

row = cand[cand["target"] == selected_target].iloc[0]
def_ra = float(row["RA"])
def_dec = float(row["Dec"])

ra = st.number_input("RA (deg)", value=def_ra if not manual else def_ra, format="%.6f")
dec = st.number_input("Dec (deg)", value=def_dec if not manual else def_dec, format="%.6f")

c1, c2, c3 = st.columns(3)
ndays = c1.number_input("Days", min_value=1, max_value=30, value=3)
step_hours = c2.number_input("Step hours", min_value=1, max_value=24, value=6)
threshold = c3.number_input("Threshold (deg)", min_value=0.0, max_value=180.0, value=30.0)

if st.button("Compute lunar curve"):
    curve = compute_lunar_curve(ra=ra, dec=dec, ndays=int(ndays), step_hours=int(step_hours), threshold=float(threshold))
    st.line_chart(curve.set_index("time_iso")[["separation_deg", "moon_phase"]])

    min_sep = float(curve["separation_deg"].min())
    mean_sep = float(curve["separation_deg"].mean())
    x1, x2 = st.columns(2)
    x1.metric("Min separation", f"{min_sep:.2f} deg")
    x2.metric("Mean separation", f"{mean_sep:.2f} deg")
    st.dataframe(curve, use_container_width=True, height=260)

st.divider()
st.header("Generate Observation Plans")

sel_targets = st.multiselect("Targets", cand["target"].tolist(), default=[selected_target])
plan_df = cand[cand["target"].isin(sel_targets)][["target", "RA", "Dec"]].copy()

p1, p2, p3, p4, p5 = st.columns(5)
tnot_count = p1.number_input("TNOT count", min_value=1, max_value=30, value=6)
tnot_interval = p2.number_input("TNOT interval (s)", min_value=10, max_value=3000, value=300)
sitian_exptime = p3.number_input("Sitian expTime (s)", min_value=10, max_value=3000, value=300)
sitian_expcount = p4.number_input("Sitian expCount", min_value=1, max_value=30, value=6)
sitian_priority = p5.number_input("Sitian priority", min_value=1, max_value=10, value=6)

if st.button("Generate both plans"):
    if plan_df.empty:
        st.warning("Select at least one target.")
    else:
        tnot_text = generate_tnot_plan_text(plan_df, count=int(tnot_count), interval=int(tnot_interval))
        sitian_text = generate_sitian_plan_text(
            plan_df,
            exptime=int(sitian_exptime),
            expcount=int(sitian_expcount),
            priority=int(sitian_priority),
        )

        st.download_button("Download TNOT plan", tnot_text, file_name="plan.txt", mime="text/plain")
        st.download_button("Download Sitian plan", sitian_text, file_name="sitian_plan.txt", mime="text/plain")

        st.text_area("TNOT preview", tnot_text, height=220)
        st.text_area("Sitian preview", sitian_text, height=220)

        saved = save_plan_outputs(tnot_text, sitian_text, GENERATED_PLAN_DIR)
        st.success(f"Saved: {saved['tnot']} and {saved['sitian']}")
