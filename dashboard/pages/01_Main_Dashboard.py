from __future__ import annotations

from pathlib import Path
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
import matplotlib.pyplot as plt

try:
    import plotly.express as px
except Exception:  # pragma: no cover
    px = None

from data_access import build_target_index, load_candidates, load_lunar, load_meta, load_timeline
from config import RESULTS_DIR

st.title("Main Dashboard")

cand = load_candidates()
meta = load_meta()
lunar = load_lunar()
timeline = load_timeline()
index_df = build_target_index()

merged = cand.copy()
if not meta.empty and "target" in meta.columns:
    merged = merged.merge(meta, on="target", how="left")
if not lunar.empty and "target" in lunar.columns:
    keep = [c for c in ["target", "Lunar_Mean_Distance", "Lunar_Min_Distance"] if c in lunar.columns]
    merged = merged.merge(lunar[keep], on="target", how="left")
if not index_df.empty:
    merged = merged.merge(index_df, on="target", how="left")

left, right = st.columns([2, 1])
with left:
    priorities = sorted(merged["Priority"].dropna().unique().tolist()) if "Priority" in merged.columns else []
    prio_selected = st.multiselect("Priority", priorities, default=priorities)
with right:
    search = st.text_input("Target name contains", value="")

filtered = merged.copy()
if prio_selected and "Priority" in filtered.columns:
    filtered = filtered[filtered["Priority"].isin(prio_selected)]
if search:
    filtered = filtered[filtered["target"].str.contains(search, case=False, na=False)]

k1, k2, k3, k4 = st.columns(4)
k1.metric("Shown targets", len(filtered))
k2.metric("Timeline rows", len(timeline))
k3.metric("Telescopes", timeline["telescope"].nunique() if "telescope" in timeline.columns else 0)
k4.metric("Bands", timeline["band"].nunique() if "band" in timeline.columns else 0)

show_cols = [
    c
    for c in [
        "target",
        "Priority",
        "Obs Time",
        "RA",
        "Dec",
        "Classification",
        "nwatch",
        "age",
        "latest_watch",
        "num_fits",
        "Lunar_Mean_Distance",
    ]
    if c in filtered.columns
]
st.dataframe(filtered[show_cols], use_container_width=True, height=700)

st.subheader("Observation Timeline")
if timeline.empty:
    st.warning("Timeline file is empty or missing.")
else:
    if px is not None:
        fig = px.scatter(
            timeline,
            x="dt",
            y="target",
            color="telescope",
            symbol="band",
            hover_data=["time_iso", "file"],
            labels={"dt": "T - T0 (days)", "target": "Target"},
            title="Observation Timeline of All Targets",
        )
        fig.update_layout(height=max(800, 28 * timeline["target"].nunique()))
        st.plotly_chart(fig, use_container_width=True)

        out_html = Path(RESULTS_DIR) / "all_targets_timeline.html"
        fig.write_html(str(out_html), include_plotlyjs="cdn")
        
        try:
            image_path = Path(RESULTS_DIR) / "all_targets_timeline.png"
            st.image(str(image_path))
        except Exception as exc:
            st.warning(f"Could not render saved PNG inline: {exc}")
    else:
        st.warning("plotly is not installed in this environment. Using fallback timeline chart.")
        tel_list = sorted(timeline["telescope"].dropna().unique().tolist()) if "telescope" in timeline.columns else []
        band_list = sorted(timeline["band"].dropna().unique().tolist()) if "band" in timeline.columns else []
        tel_color = {t: c for t, c in zip(tel_list, ["C0", "C1", "C2", "C3", "C4", "C5"])}
        band_marker = {b: m for b, m in zip(band_list, ["o", "s", "^", "D", "v", "P", "X"])}

        targets = list(dict.fromkeys(timeline["target"].tolist()))
        ymap = {t: i for i, t in enumerate(targets)}
        fig, ax = plt.subplots(figsize=(12, max(4, 0.35 * len(targets))))
        for _, row in timeline.iterrows():
            x = row.get("dt", None)
            y = ymap.get(row.get("target"), None)
            if pd.isna(x) or y is None:
                continue
            ax.scatter(
                x,
                y,
                color=tel_color.get(row.get("telescope"), "gray"),
                marker=band_marker.get(row.get("band"), "o"),
                s=22,
                alpha=0.8,
            )
        ax.set_yticks(range(len(targets)))
        ax.set_yticklabels(targets)
        ax.set_xlabel("T - T0 (days)")
        ax.set_ylabel("Target")
        ax.grid(alpha=0.2)
        st.pyplot(fig, use_container_width=True)

st.subheader("Per-target Detail")
if len(filtered):
    selected = st.selectbox("Select target", filtered["target"].tolist(), index=0)
    row = filtered[filtered["target"] == selected].iloc[0]
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Priority", int(row["Priority"]) if "Priority" in row and pd.notna(row["Priority"]) else "-")
    c2.metric("RA", f"{row['RA']:.4f}" if "RA" in row and pd.notna(row["RA"]) else "-")
    c3.metric("Dec", f"{row['Dec']:.4f}" if "Dec" in row and pd.notna(row["Dec"]) else "-")
    c4.metric("FITS files", int(row["num_fits"]) if "num_fits" in row and pd.notna(row["num_fits"]) else 0)

    detail_cols = [c for c in filtered.columns if c not in {"RA", "Dec"}]
    st.dataframe(pd.DataFrame([row[detail_cols].to_dict()]), use_container_width=True)

# st.subheader("Observation Distribution")
# dist_candidates = [
#     Path(RESULTS_DIR) / "obsrvations_distribution.html",
#     Path(RESULTS_DIR) / "observation_distribution.html",
# ]
# dist_html = next((p for p in dist_candidates if p.exists()), None)
# if dist_html is None:
#     st.info("No observation distribution HTML found in results/.")
# else:
#     with st.expander("Show observation distribution", expanded=True):
#         st.code(str(dist_html))
#         try:
#             components.html(dist_html.read_text(encoding="utf-8"), height=700, scrolling=True)
#         except Exception as exc:
#             st.warning(f"Failed to render distribution HTML: {exc}")
