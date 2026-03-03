from __future__ import annotations

import streamlit as st

from dashboard.data_access import (
    build_target_index,
    load_candidates,
    load_lunar,
    load_meta,
    load_timeline,
)

st.set_page_config(page_title="Follow_up Dashboard", page_icon="🔭", layout="wide")

st.title("Follow_up Dashboard")
st.caption("Targets, observation stats, lunar distance, planning, and photometry")

cand = load_candidates()
meta = load_meta()
timeline = load_timeline()
lunar = load_lunar()
target_index = build_target_index()

c1, c2, c3, c4, c5, c6 = st.columns(6)
c1.metric("Candidates", len(cand))
c2.metric("Observed Targets", meta["target"].nunique() if "target" in meta.columns else 0)
c3.metric("Observations", len(timeline))
c4.metric("Priority > 0", int((cand.get("Priority", 0) > 0).sum()) if len(cand) else 0)
c5.metric("Lunar Entries", len(lunar))
c6.metric("Optical Targets", len(target_index))

st.markdown("Use the left sidebar to open subpages:")
st.markdown("- Main Dashboard")
st.markdown("- Lunar Distance & Plans")
st.markdown("- Photometry")
st.markdown("- Data/Code Explorer")

st.divider()

if not cand.empty:
    preview_cols = [c for c in ["target", "Priority", "Obs Time", "RA", "Dec", "Classification"] if c in cand.columns]
    st.subheader("Candidate Preview")
    st.dataframe(cand[preview_cols], use_container_width=True, height=320)
