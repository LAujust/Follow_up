"use client";

import type { DashboardSummary } from "@/lib/types";

interface SummaryMetricsProps {
  summary: DashboardSummary;
}

export function SummaryMetrics({ summary }: SummaryMetricsProps) {
  const items = [
    { label: "Candidates", value: summary.total_candidates, color: "bg-blue-500" },
    { label: "Observed Targets", value: summary.observed_targets, color: "bg-green-500" },
    { label: "Observations", value: summary.total_observations, color: "bg-purple-500" },
    { label: "Alive", value: summary.alive_targets, color: "bg-amber-500" },
  ];

  return (
    <div className="grid grid-cols-2 lg:grid-cols-4 gap-4">
      {items.map((item) => (
        <div
          key={item.label}
          className="bg-white dark:bg-gray-800 rounded-xl border border-gray-200 dark:border-gray-700 p-5 shadow-sm"
        >
          <div className="flex items-center gap-3">
            <div className={`w-3 h-3 rounded-full ${item.color}`} />
            <span className="text-sm text-gray-500 dark:text-gray-400">{item.label}</span>
          </div>
          <p className="text-3xl font-semibold text-gray-900 dark:text-white mt-2">
            {item.value.toLocaleString()}
          </p>
        </div>
      ))}
    </div>
  );
}
