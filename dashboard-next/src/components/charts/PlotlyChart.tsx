"use client";

import { useEffect, useRef } from "react";

interface PlotlyChartProps {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  data: any[];
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  layout?: any;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  config?: any;
  className?: string;
}

export function PlotlyChart({ data, layout, config, className = "" }: PlotlyChartProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const plotlyRef = useRef<typeof import("plotly.js-dist-min") | null>(null);

  useEffect(() => {
    import("plotly.js-dist-min").then((Plotly) => {
      plotlyRef.current = Plotly;
      if (containerRef.current) {
        Plotly.newPlot(containerRef.current, data, layout || {}, {
          responsive: true,
          displayModeBar: true,
          displaylogo: false,
          ...config,
        });
      }
    });

    return () => {
      if (containerRef.current && plotlyRef.current) {
        plotlyRef.current.purge(containerRef.current);
      }
    };
  }, [data, layout, config]);

  return <div ref={containerRef} className={`w-full ${className}`} />;
}
