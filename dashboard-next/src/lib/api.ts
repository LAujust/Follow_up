import type {
  DashboardData,
  LightcurveListResponse,
  LightcurveResponse,
  LunarCurveResponse,
  PlanResponse,
} from "./types";

// Browser calls the Cloudflare Tunnel URL directly to bypass Vercel's DNS restriction.
// Resolution order:
//  1. NEXT_PUBLIC_API_URL env var (set via Vercel dashboard, inlined at build time)
//  2. Hardcoded tunnel URL (used on Vercel when env var is not set)
//  3. Empty string / relative URL (local dev via next.config.ts rewrites)
const TUNNEL_URL = "https://soa-excel-experimental-cities.trycloudflare.com";
const API_BASE = process.env.NEXT_PUBLIC_API_URL || (
  typeof window !== "undefined"
    && window.location.hostname !== "localhost"
    && window.location.hostname !== "127.0.0.1"
    ? TUNNEL_URL
    : ""
);

async function fetchJSON<T>(url: string, init?: RequestInit): Promise<T> {
  const res = await fetch(`${API_BASE}${url}`, {
    headers: { "Content-Type": "application/json", ...init?.headers },
    ...init,
  });
  if (!res.ok) {
    const text = await res.text().catch(() => "Unknown error");
    throw new Error(`API ${res.status}: ${text}`);
  }
  return res.json();
}

export async function getDashboardData(): Promise<DashboardData> {
  return fetchJSON<DashboardData>("/api/dashboard-data");
}

export async function refreshDashboardData(): Promise<{ status: string }> {
  return fetchJSON<{ status: string }>("/api/dashboard-data/refresh", {
    method: "POST",
  });
}

export async function getLunarCurve(params: {
  target?: string;
  ra?: number;
  dec?: number;
  start_time?: string;
  ndays?: number;
  step_hours?: number;
  threshold?: number;
}): Promise<LunarCurveResponse> {
  return fetchJSON<LunarCurveResponse>("/api/lunar-curve", {
    method: "POST",
    body: JSON.stringify(params),
  });
}

export async function getPlans(params: {
  target_names?: string[];
  count?: number;
  interval?: number;
  exptime?: number;
  expcount?: number;
  priority?: number;
  save_to_disk?: boolean;
}): Promise<PlanResponse> {
  return fetchJSON<PlanResponse>("/api/plans", {
    method: "POST",
    body: JSON.stringify(params),
  });
}

export async function getLightcurveList(): Promise<LightcurveListResponse> {
  return fetchJSON<LightcurveListResponse>("/api/lightcurve");
}

export async function getLightcurve(target: string): Promise<LightcurveResponse> {
  return fetchJSON<LightcurveResponse>(`/api/lightcurve/${encodeURIComponent(target)}`);
}
