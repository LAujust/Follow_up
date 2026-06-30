"use client";

import { useDashboardData } from "@/hooks/useDashboardData";
import { Sidebar } from "@/components/layout/Sidebar";
import { SummaryMetrics } from "@/components/layout/SummaryMetrics";
import { DistributionPie } from "@/components/charts/DistributionPie";
import { CumulativeLine } from "@/components/charts/CumulativeLine";

export default function OverviewPage() {
  const { data, isLoading, error } = useDashboardData();

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

  return (
    <WithSidebar>
      <div className="space-y-6">
        <div>
          <h1 className="text-2xl font-bold text-gray-900 dark:text-white">Overview</h1>
          <p className="text-sm text-gray-500 dark:text-gray-400 mt-1">
            EP Follow-up Observation Dashboard
          </p>
        </div>

        <SummaryMetrics summary={data.summary} />

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <div className="bg-white dark:bg-gray-800 rounded-xl border border-gray-200 dark:border-gray-700 p-4 shadow-sm">
            <DistributionPie timeline={data.timeline} groupBy="telescope" />
          </div>
          <div className="bg-white dark:bg-gray-800 rounded-xl border border-gray-200 dark:border-gray-700 p-4 shadow-sm">
            <DistributionPie timeline={data.timeline} groupBy="target" />
          </div>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <div className="bg-white dark:bg-gray-800 rounded-xl border border-gray-200 dark:border-gray-700 p-4 shadow-sm">
            <CumulativeLine timeline={data.timeline} mode="events" />
          </div>
          <div className="bg-white dark:bg-gray-800 rounded-xl border border-gray-200 dark:border-gray-700 p-4 shadow-sm">
            <CumulativeLine timeline={data.timeline} mode="observations" />
          </div>
        </div>
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
