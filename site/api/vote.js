import { pushVote, hasRedis } from "../lib/store.js";
import { VALID } from "../lib/contenders.js";

export default async function handler(req, res) {
  if (req.method !== "POST") {
    res.status(405).json({ error: "POST only" });
    return;
  }
  let body = req.body;
  if (typeof body === "string") { try { body = JSON.parse(body); } catch (e) { body = {}; } }
  const winner = body && body.winner, loser = body && body.loser;
  if (!VALID.has(winner) || !VALID.has(loser) || winner === loser) {
    res.status(400).json({ error: "invalid winner/loser" });
    return;
  }
  try {
    await pushVote(winner, loser);
    res.status(200).json({ ok: true, persistent: hasRedis });
  } catch (e) {
    res.status(500).json({ error: "store unavailable" });
  }
}
