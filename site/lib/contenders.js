// Valid contender keys (must match site/appdata.js) + Elo replay.
// Keeping this authoritative on the server prevents junk keys polluting the store.

export const KEYS = [
  "vantage", "equipoise", "wayfarer", "clarity", "balance",
  "cl_Winkel_Tripel", "cl_Robinson", "cl_Mollweide",
  "cl_Equal_Earth", "cl_Azimuthal_Equidist", "cl_Equirectangular",
];
export const VALID = new Set(KEYS);

// Deterministic Elo from the chronological vote log (start 1000, K=24).
export function elo(votes) {
  const R = {}, W = {}, L = {};
  for (const k of KEYS) { R[k] = 1000; W[k] = 0; L[k] = 0; }
  for (const v of votes) {
    const parts = String(v).split("|");
    const w = parts[0], l = parts[1];
    if (!(w in R) || !(l in R) || w === l) continue;
    const ew = 1 / (1 + Math.pow(10, (R[l] - R[w]) / 400));
    const K = 24;
    R[w] += K * (1 - ew); R[l] -= K * (1 - ew); W[w]++; L[l]++;
  }
  const standings = KEYS.map((k) => ({ key: k, elo: R[k], wins: W[k], losses: L[k], games: W[k] + L[k] }));
  return { standings, total: votes.length };
}
