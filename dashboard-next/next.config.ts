import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  // Proxy API requests to the FastAPI backend.
  // In dev: defaults to localhost:8000
  // In production (Vercel): set API_URL env var to the Render backend URL
  async rewrites() {
    const apiUrl = process.env.API_URL || "http://localhost:8000";
    return [
      {
        source: "/api/:path*",
        destination: `${apiUrl}/api/:path*`,
      },
    ];
  },
};

export default nextConfig;
