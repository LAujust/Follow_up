from __future__ import annotations

from pathlib import Path
import streamlit as st
import streamlit.components.v1 as components

from data_access import (
    build_target_index,
    load_candidates,
    load_lunar,
    load_meta,
    load_timeline,
)
from config import RESULTS_DIR

st.set_page_config(page_title="EP Follow-up Dashboard", page_icon="🔭", layout="centered",)

st.title("Follow up Dashboard")
st.caption("Targets, observation stats, lunar distance, planning, and photometry")

cand = load_candidates()
meta = load_meta()
timeline = load_timeline()
lunar = load_lunar()
target_index = build_target_index()

c1, c2, c3, c4 = st.columns(4)
c1.metric("Candidates", len(cand))
c2.metric("Observed Targets", meta["target"].nunique() if "target" in meta.columns else 0)
c3.metric("Observations", len(timeline))
c4.metric("Alive", int((cand.get("Priority", 0) > 2).sum()) if len(cand) else 0)

st.subheader("Observation Distribution")
dist_candidates = [
    Path(RESULTS_DIR) / "obsrvations_distribution.html",
    Path(RESULTS_DIR) / "observation_distribution.html",
]


dist_html = next((p for p in dist_candidates if p.exists()), None)
if dist_html is None:
    st.info("No observation distribution HTML found in results/.")
else:
    try:
        components.html(dist_html.read_text(encoding="utf-8"), height=600)
    except Exception as exc:
            st.warning(f"Failed to render distribution HTML: {exc}")

st.divider()

if not cand.empty:
    preview_cols = [c for c in ["target", "Priority", "Obs Time", "RA", "Dec", "Classification"] if c in cand.columns]
    st.subheader("Candidate Preview")
    st.dataframe(cand[preview_cols], use_container_width=True, height=360)