import { getVotes, hasRedis } from "../lib/store.js";
import { elo } from "../lib/contenders.js";

export default async function handler(req, res) {
  try {
    const raw = await getVotes();               // newest first
    const votes = raw.slice().reverse();        // chronological for Elo replay
    const { standings, total } = elo(votes);
    res.setHeader("cache-control", "no-store");
    res.status(200).json({ standings, total, mode: hasRedis ? "live" : "local" });
  } catch (e) {
    res.status(500).json({ error: "store unavailable" });
  }
}
