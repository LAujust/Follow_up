import { create } from "zustand";
import type { DashboardData } from "./types";

interface FilterState {
  searchQuery: string;
  priorityFilter: number[];
  selectedTarget: string | null;
  setSearchQuery: (q: string) => void;
  setPriorityFilter: (p: number[]) => void;
  setSelectedTarget: (t: string | null) => void;
  clearFilters: () => void;
}

export const useFilterStore = create<FilterState>((set) => ({
  searchQuery: "",
  priorityFilter: [],
  selectedTarget: null,
  setSearchQuery: (q) => set({ searchQuery: q }),
  setPriorityFilter: (p) => set({ priorityFilter: p }),
  setSelectedTarget: (t) => set({ selectedTarget: t }),
  clearFilters: () =>
    set({ searchQuery: "", priorityFilter: [], selectedTarget: null }),
}));

interface DataState {
  data: DashboardData | null;
  isLoading: boolean;
  error: string | null;
  setData: (data: DashboardData) => void;
  setLoading: (loading: boolean) => void;
  setError: (error: string | null) => void;
}

export const useDataStore = create<DataState>((set) => ({
  data: null,
  isLoading: true,
  error: null,
  setData: (data) => set({ data, isLoading: false, error: null }),
  setLoading: (isLoading) => set({ isLoading }),
  setError: (error) => set({ error, isLoading: false }),
}));
