"use client";

import { useState } from "react";
import { useDashboardData } from "@/hooks/useDashboardData";
import { Sidebar } from "@/components/layout/Sidebar";
import { LunarCurveChart } from "@/components/charts/LunarCurveChart";
import { PlanGeneratorForm } from "@/components/plans/PlanGeneratorForm";
import { PlanPreview } from "@/components/plans/PlanPreview";
import { getLunarCurve, getPlans } from "@/lib/api";
import type { LunarCurvePoint, PlanResponse } from "@/lib/types";

export default function PlansPage() {
  const { data, isLoading } = useDashboardData();

  // Lunar curve state
  const [selectedTarget, setSelectedTarget] = useState("");
  const [lunarCurve, setLunarCurve] = useState<LunarCurvePoint[]>([]);
  const [lunarTarget, setLunarTarget] = useState("");
  const [isComputingLunar, setIsComputingLunar] = useState(false);
  const [lunarError, setLunarError] = useState<string | null>(null);

  // Plan generation state
  const [plans, setPlans] = useState<PlanResponse | null>(null);
  const [isGenerating, setIsGenerating] = useState(false);
  const [planError, setPlanError] = useState<string | null>(null);

  const handleComputeLunar = async () => {
    if (!selectedTarget) return;
    setIsComputingLunar(true);
    setLunarError(null);
    try {
      const result = await getLunarCurve({ target: selectedTarget });
      setLunarCurve(result.curve);
      setLunarTarget(result.target);
    } catch (err) {
      setLunarError(err instanceof Error ? err.message : "Failed to compute lunar curve");
    } finally {
      setIsComputingLunar(false);
    }
  };

  const handleGeneratePlans = async (params: {
    target_names: string[];
    count: number;
    interval: number;
    exptime: number;
    expcount: number;
    priority: number;
  }) => {
    setIsGenerating(true);
    setPlanError(null);
    try {
      const result = await getPlans({ ...params, save_to_disk: true });
      setPlans(result);
    } catch (err) {
      setPlanError(err instanceof Error ? err.message : "Failed to generate plans");
    } finally {
      setIsGenerating(false);
    }
  };

  if (isLoading) {
    return (
      <WithSidebar>
        <div className="flex items-center justify-center h-64">
          <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-500" />
        </div>
      </WithSidebar>
    );
  }

  const targetNames = data?.candidates.map((c) => c.target).sort() || [];

  return (
    <WithSidebar>
      <div className="space-y-8">
        <div>
          <h1 className="text-2xl font-bold text-gray-900 dark:text-white">Lunar &amp; Plans</h1>
          <p className="text-sm text-gray-500 dark:text-gray-400 mt-1">
            Lunar distance computation and observation plan generation
          </p>
        </div>

        {/* Lunar Curve Section */}
        <section className="bg-white dark:bg-gray-800 rounded-xl border border-gray-200 dark:border-gray-700 p-6 shadow-sm">
          <h2 className="text-lg font-semibold text-gray-900 dark:text-white mb-4">
            Lunar Distance Curve
          </h2>
          <div className="flex flex-wrap items-end gap-3 mb-4">
            <div className="flex-1 min-w-[200px]">
              <label className="block text-xs font-medium text-gray-500 dark:text-gray-400 mb-1">
                Target
              </label>
              <select
                value={selectedTarget}
                onChange={(e) => setSelectedTarget(e.target.value)}
                className="w-full px-3 py-2 rounded-lg border border-gray-200 dark:border-gray-700 bg-white dark:bg-gray-800 text-sm text-gray-900 dark:text-white focus:outline-none focus:ring-2 focus:ring-blue-500"
              >
                <option value="">Select a target...</option>
                {targetNames.map((name) => (
                  <option key={name} value={name}>
                    {name}
                  </option>
                ))}
              </select>
            </div>
            <button
              onClick={handleComputeLunar}
              disabled={!selectedTarget || isComputingLunar}
              className="px-4 py-2 bg-blue-600 hover:bg-blue-700 disabled:bg-blue-400 text-white text-sm font-medium rounded-lg transition-colors"
            >
              {isComputingLunar ? "Computing..." : "Compute"}
            </button>
          </div>

          {lunarError && (
            <div className="mb-4 text-sm text-red-600 dark:text-red-400">{lunarError}</div>
          )}

          <LunarCurveChart curve={lunarCurve} target={lunarTarget} threshold={30} />
        </section>

        {/* Plan Generation Section */}
        <section className="bg-white dark:bg-gray-800 rounded-xl border border-gray-200 dark:border-gray-700 p-6 shadow-sm">
          <h2 className="text-lg font-semibold text-gray-900 dark:text-white mb-4">
            Observation Plans
          </h2>
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
            <div>
              <PlanGeneratorForm onGenerate={handleGeneratePlans} isGenerating={isGenerating} />
              {planError && (
                <div className="mt-3 text-sm text-red-600 dark:text-red-400">{planError}</div>
              )}
            </div>
            <div>
              <PlanPreview plans={plans} />
            </div>
          </div>
        </section>
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
