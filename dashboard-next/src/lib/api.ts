import type {
  DashboardData,
  LightcurveListResponse,
  LightcurveResponse,
  LunarCurveResponse,
  PlanResponse,
} from "./types";

// In dev: NEXT_PUBLIC_API_URL is empty → relative URLs → next.config.ts rewrites to localhost:8000
// In production (Vercel): set NEXT_PUBLIC_API_URL to the Cloudflare Tunnel URL (with no trailing slash)
//   so the browser calls the tunnel directly, bypassing Vercel's serverless DNS restriction.
const API_BASE = process.env.NEXT_PUBLIC_API_URL || "";

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
