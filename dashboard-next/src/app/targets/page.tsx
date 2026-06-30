"use client";

import { useDashboardData } from "@/hooks/useDashboardData";
import { Sidebar } from "@/components/layout/Sidebar";
import { FilterBar } from "@/components/layout/FilterBar";
import { TargetTable } from "@/components/tables/TargetTable";
import { TimelineChart } from "@/components/charts/TimelineChart";
import { useFilterStore } from "@/lib/store";

export default function TargetsPage() {
  const { data, isLoading, error } = useDashboardData();
  const {
    searchQuery,
    priorityFilter,
    selectedTarget,
    setSearchQuery,
    setPriorityFilter,
    setSelectedTarget,
    clearFilters,
  } = useFilterStore();

  if (isLoading) {
    return (
      <WithSidebar>
        <div className="flex items-center justify-center h-64">
          <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-500" />
        </div>
      </WithSidebar>
    );
  }

  if (error || !data) {
    return (
      <WithSidebar>
        <div className="bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-800 rounded-xl p-6 text-red-700 dark:text-red-300">
          <h2 className="font-semibold text-lg">Failed to load data</h2>
          <p className="mt-1 text-sm">{error || "No data available"}</p>
        </div>
      </WithSidebar>
    );
  }

  const selectedMeta = selectedTarget
    ? data.meta.find((m) => m.target === selectedTarget)
    : null;

  return (
    <WithSidebar>
      <div className="space-y-6">
        <div>
          <h1 className="text-2xl font-bold text-gray-900 dark:text-white">Targets</h1>
          <p className="text-sm text-gray-500 dark:text-gray-400 mt-1">
            {data.candidates.length} candidates · {data.summary.observed_targets} observed
            {selectedTarget && ` · Selected: ${selectedTarget}`}
          </p>
        </div>

        <FilterBar
          searchQuery={searchQuery}
          onSearchChange={setSearchQuery}
          priorityFilter={priorityFilter}
          onPriorityChange={setPriorityFilter}
          onClear={clearFilters}
        />

        {selectedMeta && (
          <div className="bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded-xl p-4">
            <div className="flex items-start justify-between">
              <div>
                <h3 className="font-semibold text-blue-900 dark:text-blue-200">
                  {selectedMeta.target}
                </h3>
                <div className="mt-1 flex flex-wrap gap-4 text-sm text-blue-700 dark:text-blue-300">
                  <span>T₀: {selectedMeta.T0?.toFixed(2) ?? "-"} MJD</span>
                  <span>Age: {selectedMeta.age?.toFixed(1) ?? "-"} d</span>
                  <span>Latest: {selectedMeta.latest_watch?.slice(0, 10) ?? "-"}</span>
                  <span>Watches: {selectedMeta.nwatch ?? 0}</span>
                  <span>Watches (5d): {selectedMeta.nwatch5 ?? 0}</span>
                </div>
              </div>
              <button
                onClick={() => setSelectedTarget(null)}
                className="text-blue-500 hover:text-blue-700 text-xs font-medium"
              >
                Dismiss
              </button>
            </div>
          </div>
        )}

        <div className="bg-white dark:bg-gray-800 rounded-xl border border-gray-200 dark:border-gray-700 p-4 shadow-sm">
          <TimelineChart
            timeline={data.timeline}
            selectedTarget={selectedTarget}
            onTargetClick={setSelectedTarget}
          />
        </div>

        <TargetTable
          data={data}
          selectedTarget={selectedTarget}
          onSelectTarget={setSelectedTarget}
          searchQuery={searchQuery}
          priorityFilter={priorityFilter}
        />
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
