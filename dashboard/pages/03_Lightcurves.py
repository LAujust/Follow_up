from __future__ import annotations

from pathlib import Path
import streamlit as st
import streamlit.components.v1 as components

from config import RESULTS_DIR

st.title("Lightcurves")

lc_dir = Path(RESULTS_DIR) / "lcs"
html_files = sorted(lc_dir.glob("*.html")) if lc_dir.exists() else []

if not html_files:
    st.info(f"No lightcurve HTML files found in {lc_dir}.")
    st.stop()

names = [p.name for p in html_files]
selected_name = st.selectbox("Lightcurve", names)
height = st.slider("Viewer height", min_value=400, max_value=1200, value=700, step=50)

selected_path = lc_dir / selected_name
try:
    centered_html = f"""
        <div style=\"display: flex; justify-content: center;\">
            {selected_path.read_text(encoding="utf-8")}
        </div>
        """
    components.html(centered_html, height=height)
except Exception as exc:
    st.warning(f"Failed to render lightcurve HTML: {exc}")
