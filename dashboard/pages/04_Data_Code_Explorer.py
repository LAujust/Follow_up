from __future__ import annotations

from pathlib import Path
import pandas as pd
import streamlit as st

try:
    from astropy.io import fits
except Exception:  # pragma: no cover
    fits = None

from data_access import build_repo_file_index

st.title("Data/Code Explorer")

idx = build_repo_file_index()
if idx.empty:
    st.warning("No indexed files found.")
    st.stop()

c1, c2, c3 = st.columns(3)
scope = c1.selectbox("Scope", ["all"] + sorted(idx["scope"].unique().tolist()))
suffix = c2.selectbox("Suffix", ["all"] + sorted([s for s in idx["suffix"].dropna().unique().tolist() if s]))
q = c3.text_input("Path contains", "")

view = idx.copy()
if scope != "all":
    view = view[view["scope"] == scope]
if suffix != "all":
    view = view[view["suffix"] == suffix]
if q:
    view = view[view["path"].str.contains(q, case=False, na=False)]

st.metric("Matched files", len(view))
st.dataframe(view.sort_values(["scope", "path"]), use_container_width=True, height=360)

if len(view) == 0:
    st.stop()

selected = st.selectbox("Preview file", view["path"].tolist())
path = Path(selected)

st.code(str(path))

low = path.name.lower()
text_like = {".py", ".txt", ".md", ".csv", ".json", ".yaml", ".yml", ".log"}
img_like = {".png", ".jpg", ".jpeg", ".gif"}

if low.endswith((".fits", ".fits.fz")):
    st.subheader("FITS header preview")
    if fits is None:
        st.error("FITS preview requires astropy. Install requirements.txt and redeploy.")
        st.stop()
    with fits.open(path) as hdul:
        ext = st.number_input("HDU index", min_value=0, max_value=len(hdul) - 1, value=0)
        hdr = hdul[int(ext)].header
        rows = [{"key": k, "value": str(v)} for k, v in hdr.items()]
        st.dataframe(pd.DataFrame(rows), use_container_width=True, height=400)
elif path.suffix.lower() in img_like:
    st.image(str(path), caption=str(path), use_container_width=True)
elif path.suffix.lower() == ".csv":
    df = pd.read_csv(path)
    st.dataframe(df.head(200), use_container_width=True, height=400)
elif path.suffix.lower() in text_like:
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
        lines = text.splitlines()
        max_lines = st.slider("Max lines", min_value=20, max_value=500, value=120, step=20)
        st.code("\n".join(lines[:max_lines]))
    except Exception as exc:
        st.exception(exc)
else:
    st.info("Binary or unsupported preview type. Showing metadata only.")
    st.json({"path": str(path), "size_kb": round(path.stat().st_size / 1024, 2)})
