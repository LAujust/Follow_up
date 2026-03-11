# Follow_up Dashboard

## Run

```bash
pip install -r requirements-dashboard.txt
streamlit run dashboard/app.py
```

## Pages

- Main Dashboard: targets + observation stats + timeline
- Lunar Distance & Plans: lunar curve and TNOT/Sitian plan generation
- Lightcurves: view rendered lightcurve HTMLs from results/lcs

## Notes

- Uses project data in `Candidates.csv`, `results/`, `lunar_distance/`, and `optical_data/`.
- Optical data root defaults to `~/optical_data` (override with `FOLLOWUP_OPTICAL_DIR`).
- Timeline is rendered as interactive HTML and saved to `results/all_targets_timeline.html`.
- Lightcurve HTMLs are loaded from `results/lcs/*.html`.
- Generated plans are saved to `results/generated_plans/<YYYYMMDD>/`.
