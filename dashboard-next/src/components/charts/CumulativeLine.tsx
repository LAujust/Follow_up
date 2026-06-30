"use client";

import { useMemo } from "react";
import { PlotlyChart } from "./PlotlyChart";
import type { TimelineRow } from "@/lib/types";

interface CumulativeLineProps {
  timeline: TimelineRow[];
  mode: "events" | "observations";
}

export function CumulativeLine({ timeline, mode }: CumulativeLineProps) {
  const { data, layout } = useMemo(() => {
    const sorted = [...timeline].sort(
      (a, b) => new Date(a.time_iso).getTime() - new Date(b.time_iso).getTime()
    );

    // Group by date for events mode (count unique targets per day)
    // or cumulatively add all observations
    const dateMap = new Map<string, Set<string>>();
    for (const row of sorted) {
      const date = row.time_iso.slice(0, 10);
      if (!dateMap.has(date)) dateMap.set(date, new Set());
      dateMap.get(date)!.add(row.target);
    }

    const dates = Array.from(dateMap.keys()).sort();
    let cumSum = 0;
    const cumValues = dates.map((d) => {
      if (mode === "events") {
        cumSum += dateMap.get(d)!.size;
      } else {
        cumSum += sorted.filter((r) => r.time_iso.startsWith(d)).length;
      }
      return cumSum;
    });

    const chartData = [
      {
        type: "scatter" as const,
        mode: "lines+markers" as const,
        x: dates,
        y: cumValues,
        line: { shape: "hv" as const, color: mode === "events" ? "#3b82f6" : "#8b5cf6", width: 2 },
        marker: { size: 4, color: mode === "events" ? "#3b82f6" : "#8b5cf6" },
        fill: "tozeroy",
        fillcolor: mode === "events" ? "rgba(59,130,246,0.08)" : "rgba(139,92,246,0.08)",
        name: mode === "events" ? "Events" : "Observations",
      },
    ];

    const chartLayout = {
      title: mode === "events" ? "Cumulative Events" : "Cumulative Observations",
      height: 350,
      margin: { t: 40, b: 50, l: 60, r: 20 },
      paper_bgcolor: "transparent",
      plot_bgcolor: "transparent",
      font: { size: 11 },
      xaxis: { title: "Date" },
      yaxis: { title: mode === "events" ? "Events" : "Observations" },
    };

    return { data: chartData, layout: chartLayout };
  }, [timeline, mode]);

  if (timeline.length === 0) {
    return (
      <div className="flex items-center justify-center h-64 text-gray-400 text-sm">
        No timeline data
      </div>
    );
  }

  return <PlotlyChart data={data} layout={layout} className="h-[350px]" />;
}
