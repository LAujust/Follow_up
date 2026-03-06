import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import glob
import os, sys
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_body
from astropy.table import Table
sys.path.append("/home/liangrd/Follow_up/code")
from utils import moon_phase, plot_lunar_distance


def calculate_lunar_distance(
    ra,
    dec,
    start_time=None,
    ndays=3,
    step_hours=6,
    name=None,
    threshold=30,
    plot=True,
    save_dir="./lunar.png",
):
    if start_time is None:
        t0 = Time.now()
    else:
        t0 = Time(start_time)

    times = t0 + np.arange(0, ndays * 24, step_hours) * u.hour
    target = SkyCoord(ra=ra, dec=dec, unit=u.deg)
    moon = get_body("moon", times)
    separation = moon.separation(target)
    phases, _ = moon_phase(times)

    if plot:
        plot_lunar_distance(ra,dec,ndays=ndays,step_hours=step_hours,save_dir=save_dir)
        # fig, ax = plt.subplots(figsize=(5, 3))
        # xs = [dt.strftime('%m-%d %H') for dt in times.to_datetime()]
        # ax.plot(xs, separation.deg,color='C0', label='Separation')
        # ax2 = ax.twinx()
        # ax2.plot(xs, phases * 100, color='C1', label="Moon Phase")
        # ax.axhline(threshold, ls="--", lw=1)
        # ax.locator_params(axis='x', nbins=4)
        # ax2.locator_params(axis='x', nbins=4)
        # ax.set_xlabel("Date (UTC)")
        # ax.set_ylabel("Lunar distance (deg)")
        # ax2.set_ylabel("Moon Phase (%)")
        # plt.tight_layout()
        # plt.savefig(save_dir, dpi=300, bbox_inches="tight")
        # plt.close()

    return {
        "min_distance_deg": float(separation.deg.min()),
        "max_distance_deg": float(separation.deg.max()),
        "mean_distance_deg": float(separation.deg.mean()),
        "lunar_phase": float(np.mean(phases) * 100),
    }


def main():

    # 删除旧图
    for file in glob.glob("*.png"):
        os.remove(file)

    objects = Table.read("/home/liangrd/Follow_up/Candidates.csv")
    objects = objects[objects["Priority"] > 2]

    lunar_mean = []
    lunar_min = []
    ep_names = []
    lunar_phases = []

    for obj in objects:
        ep_name = obj["EP Name"]

        print(f"Calculating lunar distance for {ep_name}...")

        result = calculate_lunar_distance(
            ra=obj["RA"],
            dec=obj["Dec"],
            start_time=Time.now(),
            ndays=3,
            step_hours=6,
            name=ep_name,
            threshold=30,
            plot=True,
            save_dir=f"./{ep_name}_lunar.png",  # 修正嵌套字符串问题
        )

        ep_names.append(ep_name)
        lunar_phases.append(result["lunar_phase"])
        lunar_mean.append(result["mean_distance_deg"])
        lunar_min.append(result["min_distance_deg"])

    # === 生成 Markdown 文件 ===
    output_md = "Candidates_lunar.md"

    with open(output_md, "w") as f:
        f.write("| EP Name | Lunar Phase | Lunar_Mean_Distance | Lunar_Min_Distance |\n")
        f.write("|---------|-------------|---------------------|--------------------|\n")

        for name, moon_phase, mean_dist, min_dist in zip(ep_names, lunar_phases, lunar_mean, lunar_min):
            f.write(f"| {name} | {moon_phase:.1f} | {mean_dist:.2f} | {min_dist:.2f}\n")

    print(f"Markdown file saved to {output_md}")
    
    
if __name__ == "__main__":
    main()