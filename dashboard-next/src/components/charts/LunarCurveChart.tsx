"use client";

import { useMemo } from "react";
import { PlotlyChart } from "./PlotlyChart";
import type { LunarCurvePoint } from "@/lib/types";

interface LunarCurveChartProps {
  curve: LunarCurvePoint[];
  target: string;
  threshold?: number;
}

export function LunarCurveChart({ curve, target, threshold = 30 }: LunarCurveChartProps) {
  const { data, layout } = useMemo(() => {
    if (curve.length === 0) return { data: [], layout: {} };

    const chartData = [
      {
        type: "scatter" as const,
        mode: "lines+markers" as const,
        x: curve.map((p) => p.time_iso),
        y: curve.map((p) => p.separation_deg),
        name: "Separation",
        line: { color: "#3b82f6", width: 2 },
        marker: {
          size: 6,
          color: curve.map((p) => (p.above_threshold ? "#22c55e" : "#ef4444")),
        },
        hovertemplate: "%{x|%Y-%m-%d %H:%M}<br>Separation: %{y:.1f}°<br>Moon phase: %{customdata:.0f}%<extra></extra>",
        customdata: curve.map((p) => p.moon_phase),
      },
      {
        type: "scatter" as const,
        mode: "lines" as const,
        x: [curve[0].time_iso, curve[curve.length - 1].time_iso],
        y: [threshold, threshold],
        name: `Threshold ${threshold}°`,
        line: { color: "#ef4444", width: 1, dash: "dash" as const },
      },
    ];

    const chartLayout = {
      title: `Lunar Distance — ${target}`,
      height: 400,
      margin: { t: 40, b: 50, l: 60, r: 60 },
      paper_bgcolor: "transparent",
      plot_bgcolor: "transparent",
      font: { size: 11 },
      xaxis: { title: "Time" },
      yaxis: { title: "Separation (deg)", range: [0, 180] },
      hovermode: "x unified" as const,
    };

    return { data: chartData, layout: chartLayout };
  }, [curve, target, threshold]);

  if (curve.length === 0) {
    return (
      <div className="flex items-center justify-center h-64 text-gray-400 text-sm">
        No lunar curve data — select a target and click compute
      </div>
    );
  }

  return <PlotlyChart data={data} layout={layout} className="h-[400px]" />;
}
