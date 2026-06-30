"use client";

import { useMemo } from "react";
import { PlotlyChart } from "./PlotlyChart";
import type { TimelineRow } from "@/lib/types";

interface DistributionPieProps {
  timeline: TimelineRow[];
  groupBy: "telescope" | "target";
}

export function DistributionPie({ timeline, groupBy }: DistributionPieProps) {
  const { data, layout } = useMemo(() => {
    const counts: Record<string, number> = {};
    for (const row of timeline) {
      const key = row[groupBy];
      counts[key] = (counts[key] || 0) + 1;
    }

    const sorted = Object.entries(counts)
      .sort((a, b) => b[1] - a[1])
      .slice(0, 15); // Top 15

    const chartData = [
      {
        type: "pie" as const,
        labels: sorted.map(([k]) => k),
        values: sorted.map(([, v]) => v),
        textinfo: "label+percent" as const,
        textposition: "outside" as const,
        hole: 0.4,
        marker: {
          line: { color: "white", width: 1 },
        },
      },
    ];

    const chartLayout = {
      title: `Observations by ${groupBy === "telescope" ? "Telescope" : "Target"}`,
      height: 400,
      margin: { t: 50, b: 20, l: 20, r: 120 },
      paper_bgcolor: "transparent",
      plot_bgcolor: "transparent",
      font: { size: 11 },
      showlegend: true,
      legend: { orientation: "v", x: 1.05 },
    };

    return { data: chartData, layout: chartLayout };
  }, [timeline, groupBy]);

  if (timeline.length === 0) {
    return (
      <div className="flex items-center justify-center h-64 text-gray-400 text-sm">
        No observation data
      </div>
    );
  }

  return <PlotlyChart data={data} layout={layout} className="h-[400px]" />;
}
