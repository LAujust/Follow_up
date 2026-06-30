"use client";

import { useMemo } from "react";
import { PlotlyChart } from "./PlotlyChart";
import type { TimelineRow } from "@/lib/types";

interface TimelineChartProps {
  timeline: TimelineRow[];
  selectedTarget: string | null;
  onTargetClick?: (target: string) => void;
}

const TELESCOPE_COLORS: Record<string, string> = {
  LCO: "#636efa",
  sitian: "#ef553b",
  TNOT: "#00cc96",
  WHUT: "#ab63fa",
  SOAR: "#ffa15a",
  Xingming: "#19d3f3",
};

const BAND_SYMBOLS: Record<string, string> = {
  g: "circle",
  r: "square",
  i: "diamond",
  z: "triangle-up",
  w: "cross",
  u: "triangle-down",
  v: "x",
  R: "hexagram",
  I: "star",
};

export function TimelineChart({ timeline, selectedTarget, onTargetClick }: TimelineChartProps) {
  const { data, layout } = useMemo(() => {
    if (timeline.length === 0) return { data: [], layout: {} };

    // Build traces: one per telescope-band combination
    const combos = new Set<string>();
    for (const row of timeline) {
      combos.add(`${row.telescope}|${row.band}`);
    }

    const traces = Array.from(combos).map((combo) => {
      const [telescope, band] = combo.split("|");
      const rows = timeline.filter((r) => r.telescope === telescope && r.band === band);
      return {
        type: "scatter" as const,
        mode: "markers" as const,
        x: rows.map((r) => r.dt),
        y: rows.map((r) => r.target),
        marker: {
          size: 8,
          color: TELESCOPE_COLORS[telescope] || "#636efa",
          symbol: BAND_SYMBOLS[band] || "circle",
          line: { width: 1, color: "rgba(0,0,0,0.3)" },
        },
        name: `${band} ${telescope}`,
        text: rows.map((r) => `${r.target}<br>dt=${r.dt.toFixed(2)}d<br>${r.time_iso}<br>${r.file}`),
        hoverinfo: "text+name",
        legendgroup: telescope,
      };
    });

    const chartLayout = {
      title: "Observation Timeline (T - T₀)",
      height: Math.max(400, timeline.length * 1.5 + 100),
      margin: { t: 40, b: 20, l: 120, r: 20 },
      paper_bgcolor: "transparent",
      plot_bgcolor: "transparent",
      font: { size: 10 },
      xaxis: {
        title: "Days since first observation",
        zeroline: true,
        zerolinecolor: "#ddd",
      },
      yaxis: {
        title: "",
        autorange: "reversed" as const,
        tickfont: { size: 9 },
      },
      legend: {
        orientation: "h" as const,
        y: -0.05,
        font: { size: 9 },
      },
      hovermode: "closest" as const,
    };

    return { data: traces, layout: chartLayout };
  }, [timeline]);

  const clickHandler = {
    event: "plotly_click" as const,
    func: (eventData: { points?: { y?: string }[] }) => {
      const target = eventData.points?.[0]?.y;
      if (target && onTargetClick) {
        onTargetClick(target);
      }
    },
  };

  return (
    <PlotlyChart
      data={data}
      layout={layout}
      config={{
        displayModeBar: true,
        displaylogo: false,
        modeBarButtonsToRemove: ["lasso2d", "select2d"],
        responsive: true,
      }}
    />
  );
}
