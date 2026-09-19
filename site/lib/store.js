// Vote persistence. Uses Upstash Redis when its env vars are present (Vercel's
// Upstash/KV Marketplace integration injects them automatically), otherwise a
// per-instance in-memory list so the site still runs before a store is attached.
//
// Design: append-only vote log (one LPUSH per vote, capped by LTRIM). No
// read-modify-write counters, so concurrent votes never clobber each other; the
// Elo standings are replayed from the log on read.

import { Redis } from "@upstash/redis";

const url = process.env.UPSTASH_REDIS_REST_URL || process.env.KV_REST_API_URL;
const token = process.env.UPSTASH_REDIS_REST_TOKEN || process.env.KV_REST_API_TOKEN;

export const hasRedis = Boolean(url && token);
const redis = hasRedis ? new Redis({ url, token }) : null;

const KEY = "pa:votes";
const CAP = 20000; // keep the log bounded (well under Upstash free-tier limits)

// Fallback store survives only within one warm serverless instance.
globalThis.__pa_mem = globalThis.__pa_mem || [];

export async function pushVote(winner, loser) {
  const rec = winner + "|" + loser + "|" + Date.now();
  if (redis) {
    await redis.lpush(KEY, rec);
    await redis.ltrim(KEY, 0, CAP - 1);
  } else {
    globalThis.__pa_mem.unshift(rec);
    if (globalThis.__pa_mem.length > CAP) globalThis.__pa_mem.length = CAP;
  }
}

// Returns the vote log, newest first.
export async function getVotes() {
  if (redis) return await redis.lrange(KEY, 0, -1);
  return globalThis.__pa_mem.slice();
}
