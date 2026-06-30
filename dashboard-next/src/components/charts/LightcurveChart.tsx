"use client";

import { useMemo } from "react";
import { PlotlyChart } from "./PlotlyChart";
import type { LightcurveRecord } from "@/lib/types";

interface LightcurveChartProps {
  data: LightcurveRecord[];
  bands: string[];
  telescopes: string[];
  target: string;
}

const BAND_COLORS: Record<string, string> = {
  g: "#19d3f3",
  r: "#ef553b",
  i: "#636efa",
  z: "#00cc96",
  u: "#ab63fa",
  w: "#ffa15a",
  v: "#ff6692",
};

const TELESCOPE_MARKERS: Record<string, string> = {
  sitian: "circle",
  TNOT: "square",
  LCO: "diamond",
  WHUT: "cross",
  SOAR: "x",
};

export function LightcurveChart({ data: lcData, bands, telescopes, target }: LightcurveChartProps) {
  const { data, layout } = useMemo(() => {
    if (lcData.length === 0) return { data: [], layout: {} };

    // Separate detections and upper limits
    const detections = lcData.filter(
      (r) => r.magpsf != null && !String(r.magpsf).startsWith("--")
    );
    const upperLimits = lcData.filter(
      (r) => r.upper_limit != null && String(r.upper_limit).trim() !== ""
    );

    const traces = [];

    // Detection traces: one per band-telescope combo
    const combos = new Set<string>();
    for (const row of detections) {
      combos.add(`${row.band}|${row.telescope}`);
    }

    for (const combo of combos) {
      const [band, telescope] = combo.split("|");
      const rows = detections.filter((r) => r.band === band && r.telescope === telescope);
      if (rows.length === 0) continue;

      traces.push({
        type: "scatter" as const,
        mode: "markers" as const,
        x: rows.map((r) => r.mean_mjd != null ? r.mean_mjd : r.dt),
        y: rows.map((r) => {
          const mag = r.mag || r.magpsf;
          return mag != null ? parseFloat(String(mag)) : null;
        }),
        error_y: {
          type: "data" as const,
          array: rows.map((r) => {
            const err = r.mag_err || r.magpsf_err;
            return err != null ? parseFloat(String(err)) : 0.1;
          }),
          visible: true,
        },
        marker: {
          size: 10,
          color: BAND_COLORS[band] || "#636efa",
          symbol: TELESCOPE_MARKERS[telescope] || "circle",
          line: { width: 1, color: "black" },
        },
        name: `${band} ${telescope}`,
        hovertemplate: "MJD: %{x:.2f}<br>Mag: %{y:.2f}<extra></extra>",
      });
    }

    // Upper limit trace
    if (upperLimits.length > 0) {
      for (const band of [...new Set(upperLimits.map((r) => r.band))]) {
        const rows = upperLimits.filter((r) => r.band === band);
        traces.push({
          type: "scatter" as const,
          mode: "markers" as const,
          x: rows.map((r) => r.mean_mjd != null ? r.mean_mjd : r.dt),
          y: rows.map((r) => parseFloat(String(r.upper_limit))),
          marker: {
            symbol: "triangle-down" as const,
            size: 10,
            color: BAND_COLORS[band] || "gray",
            line: { width: 1, color: "black" },
          },
          name: `${band} upper`,
          hovertemplate: "MJD: %{x:.2f}<br>Limit: %{y:.2f}<extra></extra>",
        });
      }
    }

    const allMags = lcData
      .map((r) => {
        const mag = r.mag || r.magpsf || r.upper_limit;
        return mag != null ? parseFloat(String(mag)) : null;
      })
      .filter((m): m is number => m != null && !isNaN(m));

    const yMin = Math.min(...allMags) - 1;
    const yMax = Math.max(...allMags) + 1;

    const chartLayout = {
      title: target,
      height: 500,
      margin: { t: 40, b: 50, l: 60, r: 60 },
      paper_bgcolor: "transparent",
      plot_bgcolor: "transparent",
      font: { size: 11 },
      xaxis: { title: "Time (MJD)" },
      yaxis: { title: "Magnitude", range: [yMax, yMin], autorange: "reversed" as const },
      legend: { orientation: "h" as const, y: -0.1, font: { size: 9 } },
      hovermode: "closest" as const,
    };

    return { data: traces, layout: chartLayout };
  }, [lcData, bands, telescopes, target]);

  if (lcData.length === 0) {
    return (
      <div className="flex items-center justify-center h-64 text-gray-400 text-sm">
        No lightcurve data for this target
      </div>
    );
  }

  return <PlotlyChart data={data} layout={layout} className="h-[500px]" />;
}
