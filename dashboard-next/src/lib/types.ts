// --- Data types matching the FastAPI backend responses ---

export interface Candidate {
  target: string;
  Priority?: number;
  "Obs Time"?: string;
  RA?: number;
  Dec?: number;
  r_err?: number;
  o_RA?: number;
  o_Dec?: number;
  Sx?: number;
  Redshift?: number;
  Classification?: string;
  GCNs?: string;
  GRB?: string;
  Fermi?: string;
  Swift?: string;
  LCO?: string;
  "LCO active"?: string;
  comments?: string;
}

export interface TimelineRow {
  target: string;
  mjd: number;
  time_iso: string;
  telescope: string;
  band: string;
  file: string;
  dt: number;
}

export interface MetaRow {
  target: string;
  T0?: number;
  age?: number;
  latest_watch?: string;
  nwatch?: number;
  nwatch5?: number;
}

export interface LunarRow {
  target: string;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  [key: string]: any;
}

export interface TargetIndexRow {
  target: string;
  num_fits: number;
  has_ps_csv: boolean;
}

export interface DashboardSummary {
  total_candidates: number;
  observed_targets: number;
  total_observations: number;
  alive_targets: number;
}

export interface DashboardData {
  summary: DashboardSummary;
  candidates: Candidate[];
  timeline: TimelineRow[];
  meta: MetaRow[];
  lunar: LunarRow[];
  target_index: TargetIndexRow[];
}

export interface LunarCurvePoint {
  time_iso: string;
  separation_deg: number;
  moon_phase: number;
  above_threshold: boolean;
}

export interface LunarCurveResponse {
  target: string;
  ra: number;
  dec: number;
  curve: LunarCurvePoint[];
}

export interface PlanResponse {
  tnot: string;
  sitian: string;
  target_count: number;
  saved_paths?: {
    tnot: string;
    sitian: string;
  };
}

export interface LightcurveRecord {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  [key: string]: any;
}

export interface LightcurveResponse {
  target: string;
  data: LightcurveRecord[];
  bands: string[];
  telescopes: string[];
  status: string;
}

export interface LightcurveTarget {
  target: string;
  source: string;
}

export interface LightcurveListResponse {
  targets: LightcurveTarget[];
}
