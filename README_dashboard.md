# Follow_up Dashboard

## Run

```bash
pip install -r requirements-dashboard.txt
streamlit run dashboard/app.py
```

## Pages

- Main Dashboard: targets + observation stats + timeline
- Lunar Distance & Plans: lunar curve and TNOT/Sitian plan generation
- Photometry: run PSF or aperture photometry and show diagnostics
- Data/Code Explorer: indexed file browser with safe previews

## Notes

- Uses project data in `Candidates.csv`, `results/`, `lunar_distance/`, and `optical_data/`.
- Optical data root defaults to `~/optical_data` (override with `FOLLOWUP_OPTICAL_DIR`).
- Timeline is rendered as interactive HTML and saved to `results/all_targets_timeline.html`.
- Photometry outputs are saved to `results/agent_photometry/<target>/<timestamp>/`.
- Generated plans are saved to `results/generated_plans/<YYYYMMDD>/`.
