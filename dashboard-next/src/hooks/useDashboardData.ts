"use client";

import { useEffect } from "react";
import { getDashboardData } from "@/lib/api";
import { useDataStore } from "@/lib/store";

export function useDashboardData() {
  const { data, isLoading, error, setData, setLoading, setError } = useDataStore();

  useEffect(() => {
    let cancelled = false;
    setLoading(true);

    getDashboardData()
      .then((result) => {
        if (!cancelled) setData(result);
      })
      .catch((err) => {
        if (!cancelled) setError(err instanceof Error ? err.message : "Failed to load data");
      });

    return () => { cancelled = true; };
  }, [setData, setLoading, setError]);

  return { data, isLoading, error };
}
