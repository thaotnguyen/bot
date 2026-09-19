# Pareto Atlas — website

A three-tab web app for the autoresearched map projections:

1. **Explore** — pick an *outline family* and set shape/area/distance priorities
   with a constrained triangle (or three sliders that always sum to 100%); the
   optimizer's best map for that recipe renders live, with Tissot ellipses.
2. **Vote** — pairwise "which map do you prefer?"; every vote is saved.
3. **Leaderboard** — all maps ranked by Elo from everyone's votes.

Static frontend (no build step) + two Vercel serverless functions in `api/`.
Data is precomputed by `../run_site_data.py` into `appdata.js`.

## Architecture

```
site/
  index.html      tabs + styling (no framework, no build)
  app.js          rendering engine, projections, the three tabs, API client
  appdata.js      precomputed per-family objective triangles + contenders + coastlines
  api/vote.js         POST {winner, loser}  -> append to the vote log
  api/leaderboard.js  GET -> Elo standings replayed from the log
  lib/store.js        vote persistence (Upstash Redis, or in-memory fallback)
  lib/contenders.js   valid keys + Elo replay
```

Votes are an **append-only log** (one `LPUSH` per vote, `LTRIM`-capped); Elo is
replayed on read. No read-modify-write counters, so concurrent votes never clobber.

## Deploy to Vercel (≈3 minutes)

**Option A — dashboard (no CLI):**
1. Push this repo to GitHub (already on branch `claude/map-projection-autoresearch-nen5zb`).
2. vercel.com → **Add New… → Project** → import the repo.
3. Set **Root Directory = `site`**. Framework preset: **Other**. Deploy.
4. In the project, **Storage → Create → Upstash for Redis** (free). Vercel injects
   `UPSTASH_REDIS_REST_URL` / `UPSTASH_REDIS_REST_TOKEN` automatically.
5. **Redeploy** (Deployments → ⋯ → Redeploy). Votes now persist and are shared.

**Option B — CLI:**
```bash
npm i -g vercel
cd site
vercel            # link + first deploy (asks to log in)
# add storage:
vercel storage create   # choose Upstash for Redis, link it
vercel --prod
```

Until a store is attached the site still works — votes just live in memory per
serverless instance (non-persistent). The leaderboard shows a "local/live" badge.

## Local dev

```bash
cd site
npm install
npx vercel dev         # serves static + /api locally (needs a Vercel login)
# or open index.html directly: Explore works fully; Vote/Leaderboard fall back
# to browser localStorage when /api isn't reachable.
```

Set env for local persistence (optional): create `site/.env.local` with
`UPSTASH_REDIS_REST_URL=...` and `UPSTASH_REDIS_REST_TOKEN=...`.

## Regenerate the map data

From the repo root (pure-Python, no deps):
```bash
python3 run_site_data.py     # ~5-8 min -> site/appdata.js
```
