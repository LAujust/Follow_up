"use client";

import { useEffect, useState } from "react";
import { Sidebar } from "@/components/layout/Sidebar";
import { LightcurveChart } from "@/components/charts/LightcurveChart";
import { getLightcurve, getLightcurveList } from "@/lib/api";
import type { LightcurveResponse } from "@/lib/types";

export default function LightcurvesPage() {
  const [targets, setTargets] = useState<string[]>([]);
  const [selectedTarget, setSelectedTarget] = useState("");
  const [lcData, setLcData] = useState<LightcurveResponse | null>(null);
  const [isLoadingList, setIsLoadingList] = useState(true);
  const [isLoadingLc, setIsLoadingLc] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Load target list on mount
  useEffect(() => {
    getLightcurveList()
      .then((res) => {
        setTargets(res.targets.map((t) => t.target).sort());
        if (res.targets.length > 0) {
          setSelectedTarget(res.targets[0].target);
        }
      })
      .catch((err) => setError(err.message))
      .finally(() => setIsLoadingList(false));
  }, []);

  // Load lightcurve data when target changes
  useEffect(() => {
    if (!selectedTarget) return;
    setIsLoadingLc(true);
    setError(null);
    getLightcurve(selectedTarget)
      .then((res) => setLcData(res))
      .catch((err) => setError(err.message))
      .finally(() => setIsLoadingLc(false));
  }, [selectedTarget]);

  if (isLoadingList) {
    return (
      <WithSidebar>
        <div className="flex items-center justify-center h-64">
          <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-500" />
        </div>
      </WithSidebar>
    );
  }

  return (
    <WithSidebar>
      <div className="space-y-6">
        <div className="flex items-start justify-between">
          <div>
            <h1 className="text-2xl font-bold text-gray-900 dark:text-white">Lightcurves</h1>
            <p className="text-sm text-gray-500 dark:text-gray-400 mt-1">
              {targets.length} targets with lightcurve data
            </p>
          </div>
          <div className="w-72">
            <label className="block text-xs font-medium text-gray-500 dark:text-gray-400 mb-1">
              Select Target
            </label>
            <select
              value={selectedTarget}
              onChange={(e) => setSelectedTarget(e.target.value)}
              className="w-full px-3 py-2 rounded-lg border border-gray-200 dark:border-gray-700 bg-white dark:bg-gray-800 text-sm text-gray-900 dark:text-white focus:outline-none focus:ring-2 focus:ring-blue-500"
            >
              {targets.map((name) => (
                <option key={name} value={name}>
                  {name}
                </option>
              ))}
            </select>
          </div>
        </div>

        {error && (
          <div className="bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-800 rounded-xl p-4 text-sm text-red-700 dark:text-red-300">
            {error}
          </div>
        )}

        {isLoadingLc ? (
          <div className="flex items-center justify-center h-64">
            <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-500" />
          </div>
        ) : lcData ? (
          <div className="bg-white dark:bg-gray-800 rounded-xl border border-gray-200 dark:border-gray-700 p-4 shadow-sm">
            <LightcurveChart
              data={lcData.data}
              bands={lcData.bands}
              telescopes={lcData.telescopes}
              target={lcData.target}
            />
          </div>
        ) : (
          <div className="flex items-center justify-center h-64 text-gray-400 text-sm border-2 border-dashed border-gray-200 dark:border-gray-700 rounded-xl">
            Select a target to view its lightcurve
          </div>
        )}
      </div>
    </WithSidebar>
  );
}

function WithSidebar({ children }: { children: React.ReactNode }) {
  return (
    <div className="flex min-h-screen">
      <Sidebar />
      <main className="flex-1 ml-64 p-6 lg:p-8">{children}</main>
    </div>
  );
}
