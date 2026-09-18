# Pareto Atlas — autoresearching a better world map projection

> You cannot flatten a sphere without distorting it (Gauss's *Theorema Egregium*).
> So "the best map" is not a fact to look up — it is a **multi-objective optimization**.
> This project lets a machine **hill-climb its own projections**, traces the whole
> **Pareto frontier** of shape-vs-area distortion, and finds maps that measurably
> beat the classics — including the Winkel Tripel, National Geographic's chosen
> compromise.

Everything here is **pure-Python standard library** (no numpy/scipy) plus a
**zero-dependency web explorer**, so the whole result is reproducible from source
on a stock Python 3 install.

**▶ Live interactive explorer + poll:** https://claude.ai/artifact/2u277G67CgFVJPmBDXaa9h
*(private to the artifact owner; share from the page's Share menu to let others vote).*

---

## TL;DR results

Scored on a 6,390-point global grid, area-weighted, using Tissot's indicatrix.
Distortion is dimensionless: `shape = ⟨ln(a/b)²⟩`, `area = Var(ln a·b)`
(both **0** for a perfect map; see *Validation*). Lower is better.

| Projection | shape ↓ | area ↓ | **combined** ↓ | mean angular |
|---|---|---|---|---|
| **Equipoise** *(discovered, balanced)* | 0.109 | 0.104 | **0.214** | 18.3° |
| **Vantage** *(discovered)* | **0.201** | **0.040** | 0.241 | 24.4° |
| Winkel Tripel *(Nat Geo, 1998)* | 0.242 | 0.051 | 0.292 | 27.3° |
| Robinson *(Rand McNally)* | 0.300 | 0.045 | 0.345 | 29.1° |
| Equal Earth *(2018)* | 0.516 | 0.000 | 0.516 | 36.4° |
| **Clarity** *(discovered, near-conformal)* | **0.011** | 0.424 | 0.435 | **5.8°** |

* **Equipoise** carries **27 % less total distortion than the Winkel Tripel.**
* **Vantage** *dominates* the Winkel Tripel — **lower shape error *and* lower area
  error simultaneously** (0.201 < 0.242 and 0.040 < 0.051). It is strictly better
  by both measures at once.
* **Clarity** achieves **5.8° mean angular distortion for a whole-world map**
  (Winkel: 27.3°) by spending its budget entirely on shape.
* **Anthropocene** (population-weighted) scores **0.046 combined where people
  actually live**, vs **0.084** for the best classic (Robinson) under the same
  weighting — ~45 % better for humans (while deliberately worse over empty ocean).

> ⚠️ **Honesty note.** "Beats Winkel Tripel" is true *for this shape+area metric on
> this grid*. Other published metrics (e.g. Goldberg–Gott, which also charges for
> flexion, skewness, distances and boundary cuts) could rank differently. The point
> is the *method*: given any differentiable objective, autoresearch finds its
> frontier — you are not limited to the handful of projections a human happened to
> invent.

---

## The idea

1. **A projection is 18 numbers.** We write the map as two symmetric polynomials
   in longitude/latitude (odd in lon, even in lat for `x`; the reverse for `y`),
   so it is automatically symmetric about the equator and central meridian. This
   family is *not* any named projection — equirectangular is a single point in it,
   and Winkel/Mollweide/etc. are transcendental and land on no exact coefficient
   vector.

2. **Measure distortion exactly.** At each point, Tissot's indicatrix gives the two
   principal scale factors `a ≥ b` (singular values of the sphere→plane Jacobian,
   corrected for the `cos(lat)` metric). From them:
   * conformality error `ln(a/b)` — 0 ⟺ angles preserved,
   * equal-area error `ln(a·b)` — constant ⟺ areas preserved (we use its
     scale-invariant variance).

3. **Hill-climb.** Adam gradient descent minimizes `(1−t)·shape + t·area`, with a
   barrier that forbids the map from folding over itself.

4. **Sweep the frontier.** Varying `t` from 0→1 traces a continuum of new
   projections — each optimal for its own shape/area trade — via warm-start
   continuation.

### Creative extensions (beyond rehashing old maps)

* **Population weighting** — weight error by *where people live* instead of by
  acreage, producing a projection optimized to be fair to humans, not to empty
  ocean and ice. (Approximate demographic model in `mapopt/population.py`, trivially
  swappable for a real GPW raster.)
* **The projection election** — a live pairwise poll ("which map do you prefer?")
  feeds an Elo / Bradley–Terry ranking. Human taste is the one objective the math
  can't compute; folding it back in as a third axis is, literally, **RLHF for
  cartography**. See the web explorer.

---

## Reproduce it

```bash
python3 run_research.py     # ~3–4 min: sweeps frontiers, benchmarks classics,
                            #           checks novelty -> results/results.json
python3 export_web.py       # bundles results + Natural Earth coastlines -> web/appdata.js
# then serve web/ (any static server) and open index.html
python3 -m http.server -d web 8000
```

### Validation (why you can trust the numbers)

`run_research.py` scores classic projections with known theory on the same grid:

* every equal-area projection (Sinusoidal, Mollweide, Hammer, Eckert IV, Equal
  Earth, Lambert CEA) scores **`area = 0.00000`**,
* the conformal projection (Mercator) scores **`shape = 0.00000`**,
* the fast analytic path for our family matches an independent finite-difference
  path to 6 decimals.

If the distortion code were wrong, these identities would fail.

### Novelty check

For each discovered map we fit the best possible per-axis rescaling of every
classic and report the residual. Equipoise still differs from its nearest classic
(Robinson) by **6.0 %** of the map's size; the near-conformal Clarity by **11.9 %**.
These coordinates sit on no historical projection.

---

## The web explorer & poll

`web/index.html` is a single self-contained page (only Google Fonts + the local
data bundle) that:

* renders **real Natural Earth coastlines** under any projection,
* lets you **drag along the discovered Pareto frontier** and watch the map morph,
  with live Tissot indicatrices and a shape-vs-area scatter showing our frontier
  sitting *inside* the cloud of classics,
* runs **the projection election** — pairwise voting with a live Elo leaderboard.
  When published as a Claude Artifact it uses the shared `db` capability so votes
  aggregate across everyone who opens it; opened as a plain file it falls back to
  a per-browser tally, so it always works.

---

## Layout

```
mapopt/
  family.py       parametric projection family + grid + precomputed basis
  metrics.py      Tissot distortion (fast analytic path + generic finite-diff)
  classics.py     correctly-formulated baseline projections
  population.py   approximate human-population weighting field
  optimize.py     Adam hill-climber + Pareto sweep (warm-start continuation)
run_research.py   the pipeline -> results/results.json
export_web.py     results + coastlines -> web/appdata.js
web/index.html    interactive explorer + projection-election poll
results/          machine-written outputs (checked in for convenience)
```

Data: coastlines from [Natural Earth](https://www.naturalearthdata.com/) (110m,
public domain). Method inspired by the distortion-metric tradition of Airy,
Kavrayskiy, Tissot, and Goldberg & Gott.
