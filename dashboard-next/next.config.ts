import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  // Rewrites proxy /api/* to the FastAPI backend.
  // In dev (NEXT_PUBLIC_API_URL is empty): rewrites to localhost:8000
  // In production (NEXT_PUBLIC_API_URL set): rewrites are disabled;
  //   the browser calls the Cloudflare Tunnel URL directly.
  async rewrites() {
    if (process.env.NEXT_PUBLIC_API_URL) return [];
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
