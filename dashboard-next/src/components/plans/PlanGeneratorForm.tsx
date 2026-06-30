"use client";

import { useState } from "react";

interface PlanGeneratorFormProps {
  onGenerate: (params: {
    target_names: string[];
    count: number;
    interval: number;
    exptime: number;
    expcount: number;
    priority: number;
  }) => void;
  isGenerating: boolean;
}

export function PlanGeneratorForm({ onGenerate, isGenerating }: PlanGeneratorFormProps) {
  const [targetInput, setTargetInput] = useState("");
  const [count, setCount] = useState(6);
  const [interval, setInterval] = useState(300);
  const [exptime, setExptime] = useState(300);
  const [expcount, setExpcount] = useState(6);
  const [priority, setPriority] = useState(6);

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    const target_names = targetInput
      .split(/[\n,]+/)
      .map((s) => s.trim())
      .filter(Boolean);

    onGenerate({
      target_names: target_names.length > 0 ? target_names : [],
      count,
      interval,
      exptime,
      expcount,
      priority,
    });
  };

  return (
    <form onSubmit={handleSubmit} className="space-y-4">
      <div>
        <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
          Targets (one per line, comma-separated; empty = all)
        </label>
        <textarea
          value={targetInput}
          onChange={(e) => setTargetInput(e.target.value)}
          rows={3}
          className="w-full px-3 py-2 rounded-lg border border-gray-200 dark:border-gray-700 bg-white dark:bg-gray-800 text-sm text-gray-900 dark:text-white placeholder-gray-400 focus:outline-none focus:ring-2 focus:ring-blue-500 font-mono"
          placeholder="EP260116a&#10;EP260211a"
        />
      </div>

      <div className="grid grid-cols-2 sm:grid-cols-5 gap-3">
        <div>
          <label className="block text-xs font-medium text-gray-500 dark:text-gray-400 mb-1">
            TNOT Count
          </label>
          <input
            type="number"
            value={count}
            onChange={(e) => setCount(Number(e.target.value))}
            min={1}
            max={100}
            className="w-full px-2 py-1.5 rounded border border-gray-200 dark:border-gray-700 bg-white dark:bg-gray-800 text-sm text-gray-900 dark:text-white focus:outline-none focus:ring-2 focus:ring-blue-500"
          />
        </div>
        <div>
          <label className="block text-xs font-medium text-gray-500 dark:text-gray-400 mb-1">
            TNOT Interval (s)
          </label>
          <input
            type="number"
            value={interval}
            onChange={(e) => setInterval(Number(e.target.value))}
            min={1}
            className="w-full px-2 py-1.5 rounded border border-gray-200 dark:border-gray-700 bg-white dark:bg-gray-800 text-sm text-gray-900 dark:text-white focus:outline-none focus:ring-2 focus:ring-blue-500"
          />
        </div>
        <div>
          <label className="block text-xs font-medium text-gray-500 dark:text-gray-400 mb-1">
            Sitian ExpTime (s)
          </label>
          <input
            type="number"
            value={exptime}
            onChange={(e) => setExptime(Number(e.target.value))}
            min={1}
            className="w-full px-2 py-1.5 rounded border border-gray-200 dark:border-gray-700 bg-white dark:bg-gray-800 text-sm text-gray-900 dark:text-white focus:outline-none focus:ring-2 focus:ring-blue-500"
          />
        </div>
        <div>
          <label className="block text-xs font-medium text-gray-500 dark:text-gray-400 mb-1">
            Sitian ExpCount
          </label>
          <input
            type="number"
            value={expcount}
            onChange={(e) => setExpcount(Number(e.target.value))}
            min={1}
            max={100}
            className="w-full px-2 py-1.5 rounded border border-gray-200 dark:border-gray-700 bg-white dark:bg-gray-800 text-sm text-gray-900 dark:text-white focus:outline-none focus:ring-2 focus:ring-blue-500"
          />
        </div>
        <div>
          <label className="block text-xs font-medium text-gray-500 dark:text-gray-400 mb-1">
            Priority
          </label>
          <input
            type="number"
            value={priority}
            onChange={(e) => setPriority(Number(e.target.value))}
            min={1}
            max={10}
            className="w-full px-2 py-1.5 rounded border border-gray-200 dark:border-gray-700 bg-white dark:bg-gray-800 text-sm text-gray-900 dark:text-white focus:outline-none focus:ring-2 focus:ring-blue-500"
          />
        </div>
      </div>

      <button
        type="submit"
        disabled={isGenerating}
        className="px-4 py-2 bg-blue-600 hover:bg-blue-700 disabled:bg-blue-400 text-white text-sm font-medium rounded-lg transition-colors"
      >
        {isGenerating ? "Generating..." : "Generate Plans"}
      </button>
    </form>
  );
}
