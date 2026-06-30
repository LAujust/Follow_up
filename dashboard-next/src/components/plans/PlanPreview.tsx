"use client";

import type { PlanResponse } from "@/lib/types";

interface PlanPreviewProps {
  plans: PlanResponse | null;
}

export function PlanPreview({ plans }: PlanPreviewProps) {
  if (!plans) {
    return (
      <div className="flex items-center justify-center h-48 text-gray-400 text-sm border-2 border-dashed border-gray-200 dark:border-gray-700 rounded-xl">
        Generate plans to see preview
      </div>
    );
  }

  const copyToClipboard = async (text: string) => {
    try {
      await navigator.clipboard.writeText(text);
    } catch {
      // Fallback for environments without clipboard API
    }
  };

  return (
    <div className="space-y-4">
      <p className="text-xs text-gray-500">
        Generated for {plans.target_count} target{plans.target_count !== 1 ? "s" : ""}
      </p>

      {/* TNOT Plan */}
      <div>
        <div className="flex items-center justify-between mb-1">
          <h4 className="text-sm font-medium text-gray-700 dark:text-gray-300">TNOT Plan</h4>
          <button
            onClick={() => copyToClipboard(plans.tnot)}
            className="text-xs text-blue-500 hover:text-blue-700"
          >
            Copy
          </button>
        </div>
        <pre className="bg-gray-50 dark:bg-gray-900 border border-gray-200 dark:border-gray-700 rounded-lg p-3 text-xs font-mono text-gray-700 dark:text-gray-300 overflow-x-auto max-h-60 overflow-y-auto">
          {plans.tnot}
        </pre>
      </div>

      {/* Sitian Plan */}
      <div>
        <div className="flex items-center justify-between mb-1">
          <h4 className="text-sm font-medium text-gray-700 dark:text-gray-300">Sitian Plan</h4>
          <button
            onClick={() => copyToClipboard(plans.sitian)}
            className="text-xs text-blue-500 hover:text-blue-700"
          >
            Copy
          </button>
        </div>
        <pre className="bg-gray-50 dark:bg-gray-900 border border-gray-200 dark:border-gray-700 rounded-lg p-3 text-xs font-mono text-gray-700 dark:text-gray-300 overflow-x-auto max-h-60 overflow-y-auto">
          {plans.sitian}
        </pre>
      </div>

      {plans.saved_paths && (
        <div className="text-xs text-gray-400 flex gap-2">
          {plans.saved_paths.tnot && (
            <span>Saved: {plans.saved_paths.tnot}</span>
          )}
        </div>
      )}
    </div>
  );
}
