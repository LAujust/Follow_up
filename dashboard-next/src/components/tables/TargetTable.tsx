"use client";

import { useMemo, useState } from "react";
import type { Candidate, DashboardData } from "@/lib/types";

interface TargetTableProps {
  data: DashboardData;
  selectedTarget: string | null;
  onSelectTarget: (target: string | null) => void;
  searchQuery: string;
  priorityFilter: number[];
}

export function TargetTable({
  data,
  selectedTarget,
  onSelectTarget,
  searchQuery,
  priorityFilter,
}: TargetTableProps) {
  const [sortField, setSortField] = useState<string>("Priority");
  const [sortAsc, setSortAsc] = useState(false);

  const filtered = useMemo(() => {
    let items = data.candidates;

    if (searchQuery) {
      const q = searchQuery.toLowerCase();
      items = items.filter((c) => c.target.toLowerCase().includes(q));
    }

    if (priorityFilter.length > 0) {
      items = items.filter((c) => c.Priority !== undefined && priorityFilter.includes(c.Priority));
    }

    // Add observation count info
    const timelineCounts = new Map<string, number>();
    for (const row of data.timeline) {
      timelineCounts.set(row.target, (timelineCounts.get(row.target) || 0) + 1);
    }

    const metaMap = new Map(data.meta.map((m) => [m.target, m]));

    items = items.map((c) => ({
      ...c,
      _obsCount: timelineCounts.get(c.target) || 0,
      _age: metaMap.get(c.target)?.age,
      _nwatch: metaMap.get(c.target)?.nwatch,
    }));

    items.sort((a, b) => {
      const aVal = (a as unknown as Record<string, unknown>)[sortField];
      const bVal = (b as unknown as Record<string, unknown>)[sortField];
      if (aVal == null && bVal == null) return 0;
      if (aVal == null) return 1;
      if (bVal == null) return -1;
      const cmp = typeof aVal === "number" ? aVal - (bVal as number) : String(aVal).localeCompare(String(bVal));
      return sortAsc ? cmp : -cmp;
    });

    return items;
  }, [data.candidates, data.timeline, data.meta, searchQuery, priorityFilter, sortField, sortAsc]);

  const toggleSort = (field: string) => {
    if (sortField === field) {
      setSortAsc(!sortAsc);
    } else {
      setSortField(field);
      setSortAsc(false);
    }
  };

  const columns = [
    { key: "target", label: "Target", sortable: true },
    { key: "Priority", label: "Pri", sortable: true },
    { key: "_obsCount", label: "Obs", sortable: true },
    { key: "Classification", label: "Class", sortable: true },
    { key: "RA", label: "RA", sortable: true },
    { key: "Dec", label: "Dec", sortable: true },
    { key: "Redshift", label: "z", sortable: true },
    { key: "_age", label: "Age (d)", sortable: true },
    { key: "Obs Time", label: "T₀", sortable: true },
  ] as const;

  const sortArrow = (field: string) => {
    if (sortField !== field) return "";
    return sortAsc ? " ▲" : " ▼";
  };

  return (
    <div className="overflow-x-auto rounded-lg border border-gray-200 dark:border-gray-700">
      <table className="min-w-full text-sm divide-y divide-gray-200 dark:divide-gray-700">
        <thead className="bg-gray-50 dark:bg-gray-800">
          <tr>
            {columns.map((col) => (
              <th
                key={col.key}
                className={`px-3 py-2.5 text-left text-xs font-semibold text-gray-500 dark:text-gray-400 uppercase tracking-wider ${
                  col.sortable ? "cursor-pointer hover:text-gray-700 dark:hover:text-gray-200" : ""
                }`}
                onClick={() => col.sortable && toggleSort(col.key)}
              >
                {col.label}
                {sortArrow(col.key)}
              </th>
            ))}
          </tr>
        </thead>
        <tbody className="bg-white dark:bg-gray-900 divide-y divide-gray-100 dark:divide-gray-800">
          {filtered.map((c) => {
            const isSelected = c.target === selectedTarget;
            const rec = c as unknown as Record<string, unknown>;
            return (
              <tr
                key={c.target}
                className={`cursor-pointer transition-colors ${
                  isSelected
                    ? "bg-blue-50 dark:bg-blue-900/20"
                    : "hover:bg-gray-50 dark:hover:bg-gray-800"
                }`}
                onClick={() => onSelectTarget(isSelected ? null : c.target)}
              >
                <td className="px-3 py-2 font-medium text-gray-900 dark:text-white whitespace-nowrap">
                  {c.target}
                </td>
                <td className="px-3 py-2 whitespace-nowrap">
                  <PriorityBadge priority={c.Priority ?? 0} />
                </td>
                <td className="px-3 py-2 text-gray-600 dark:text-gray-300">
                  {rec._obsCount as number}
                </td>
                <td className="px-3 py-2 text-gray-500 dark:text-gray-400 max-w-[120px] truncate">
                  {c.Classification || "-"}
                </td>
                <td className="px-3 py-2 text-gray-500 dark:text-gray-400 font-mono text-xs">
                  {c.RA?.toFixed(4) ?? "-"}
                </td>
                <td className="px-3 py-2 text-gray-500 dark:text-gray-400 font-mono text-xs">
                  {c.Dec?.toFixed(4) ?? "-"}
                </td>
                <td className="px-3 py-2 text-gray-500 dark:text-gray-400">
                  {c.Redshift != null ? c.Redshift.toFixed(2) : "-"}
                </td>
                <td className="px-3 py-2 text-gray-500 dark:text-gray-400">
                  {rec._age != null ? (rec._age as number).toFixed(1) : "-"}
                </td>
                <td className="px-3 py-2 text-gray-500 dark:text-gray-400 text-xs">
                  {c["Obs Time"] ? c["Obs Time"].slice(0, 10) : "-"}
                </td>
              </tr>
            );
          })}
        </tbody>
      </table>

      {filtered.length === 0 && (
        <div className="text-center py-12 text-gray-400 text-sm">
          No targets match the current filters
        </div>
      )}
    </div>
  );
}

function PriorityBadge({ priority }: { priority: number }) {
  const colorMap: Record<number, string> = {
    1: "bg-red-100 text-red-700 dark:bg-red-900/30 dark:text-red-300",
    2: "bg-orange-100 text-orange-700 dark:bg-orange-900/30 dark:text-orange-300",
    3: "bg-amber-100 text-amber-700 dark:bg-amber-900/30 dark:text-amber-300",
    4: "bg-yellow-100 text-yellow-700 dark:bg-yellow-900/30 dark:text-yellow-300",
    5: "bg-green-100 text-green-700 dark:bg-green-900/30 dark:text-green-300",
  };
  return (
    <span
      className={`inline-flex items-center px-2 py-0.5 rounded-full text-xs font-medium ${
        colorMap[priority] || "bg-gray-100 text-gray-600 dark:bg-gray-800 dark:text-gray-300"
      }`}
    >
      {priority}
    </span>
  );
}
