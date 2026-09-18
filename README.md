# Pareto Atlas — autoresearching a better world map projection

> You cannot flatten a sphere without distorting it (Gauss's *Theorema Egregium*).
> The "best map" is therefore a **three-way trade-off** — you can preserve
> **shapes**, **areas**, or **distances**, but never all three. Instead of
> picking a favorite compromise, this project lets a machine **hill-climb the
> entire 3-objective Pareto surface** and finds maps that measurably beat the
> classics — including a projection that beats the **Winkel Tripel on all three
> axes at once.**

**▶ Live interactive explorer + poll:** https://claude.ai/artifact/2u277G67CgFVJPmBDXaa9h
*(private to the artifact owner; share from the page's Share menu to let others vote).*

Everything here is **pure-Python standard library** (no numpy/scipy) plus a
**zero-dependency web explorer**, so the whole result is reproducible from source
on a stock Python 3 install.

---

## The three objectives

All are dimensionless and **0** for a perfect map (see *Validation*):

| axis | what it means | metric |
|---|---|---|
| **shape** | are local angles/shapes right? (conformality) | `⟨ln(a/b)²⟩` from Tissot's indicatrix (local) |
| **area** | are sizes right? (equal-area) | `Var(ln a·b)` (local) |
| **distance** | do straight-line map distances match true great-circle distances? | `Var(ln( d_map / d_globe ))` over thousands of point pairs (**global**) |

`a, b` are Tissot's two principal scale factors (singular values of the
sphere→plane Jacobian). Shape and area are *local*; distance is *global* — a
genuinely independent third axis (e.g. the Azimuthal Equidistant is great on
distance yet poor on shape and area).

## TL;DR results

Scored on a 6,390-point grid + 11,151 point-pairs. Lower is better; `norm` is the
sum of the three errors each normalized to equirectangular (so equirect = 3.00).

| Projection | shape | area | distance | **norm** |
|---|---|---|---|---|
| **Vantage** *(discovered)* | **0.215** | **0.042** | **0.157** | **1.73** |
| **Equipoise** *(discovered, all-rounder)* | 0.196 | 0.053 | 0.157 | 1.73 |
| Winkel Tripel *(Nat Geo)* | 0.242 | 0.051 | 0.167 | 1.93 |
| Robinson | 0.300 | 0.045 | 0.190 | 2.22 |
| Equal Earth *(2018)* | 0.516 | 0.000 | 0.176 | 2.70 |
| Azimuthal Equidistant | 0.749 | 0.355 | 0.160 | 5.56 |
| **Wayfarer** *(discovered, distance corner)* | 0.746 | 0.073 | **0.131** | – |
| **Clarity** *(discovered, shape corner)* | **0.003** | 0.515 | 0.207 | – |

* **Vantage beats the Winkel Tripel on all three axes simultaneously** — lower
  shape *and* area *and* distance error (−11% / −18% / −6%), for ~10% less total
  distortion. Found by a targeted solve that minimizes the *worst* ratio to Winkel.
* **Wayfarer sets a new distance record: 0.131**, lower than *every* classic —
  including the purpose-built Azimuthal Equidistant (0.160).
* **Clarity** reaches **3.3° mean angular distortion** for a whole-world map.
* **Equipoise** (minimax all-rounder) beats Winkel clearly on shape and distance
  and ties on area — the most *even* map on the surface.

> ⚠️ **Honesty note.** These rankings hold *for these three metrics on this grid*.
> Other published metrics (e.g. Goldberg–Gott, which also charges for flexion,
> skewness and boundary cuts) could rank differently. The durable result is the
> **method**: given any differentiable objectives, autoresearch maps their whole
> Pareto surface — you are not limited to the projections a human happened to invent.

---

## The method (autoresearch, in four moves)

1. **Parametrize.** A projection = two symmetric polynomials in lon/lat (18 free
   numbers). Not any named projection; equirectangular is one point in the space.
2. **Measure ×3.** Shape and area from Tissot's indicatrix; distance from
   great-circle-vs-map comparisons over a fixed pair set.
3. **Hill-climb.** Adam gradient descent minimizes a weighted blend of the three
   (normalized so they're comparable), with a barrier that forbids the map from
   folding over itself.
4. **Sweep the surface.** Walking the weight *simplex* (barycentric grid) with
   warm-start continuation traces the full 3-objective Pareto surface. The
   triangle in the web explorer **is** that surface, made draggable.

### Human-in-the-loop (RLHF for cartography)

Math ranks distortion; humans rank *maps*. The web explorer runs a pairwise
**"projection election"** — vote for the map you prefer, feeding a live Elo /
Bradley–Terry ranking. That human-preference signal is the one objective the
optimizer can't compute, and a future round could fold it back in as a fourth axis.

---

## Reproduce it

```bash
python3 run_research.py     # ~3-4 min: 3-objective simplex sweep + targeted
                            #           Winkel-dominator + benchmarks -> results/results.json
python3 export_web.py       # results + Natural Earth coastlines -> web/appdata.js
python3 -m http.server -d web 8000    # then open http://localhost:8000
```

### Validation (why you can trust the numbers)

`run_research.py` scores classic projections with known theory on the same grid:

* every equal-area projection (Sinusoidal, Mollweide, Hammer, Eckert IV, Equal
  Earth, Lambert CEA) scores **`area = 0.00000`**,
* the conformal projection (Mercator) scores **`shape = 0.00000`**,
* the Azimuthal Equidistant is (correctly) the best *classic* on distance,
* the fast analytic path for our family matches an independent finite-difference
  path to 6 decimals.

### Novelty check

After the best possible per-axis rescaling of every classic, the all-rounder
still differs from its nearest classic (Winkel Tripel) by **3.7%** of the map's
size, and the distance champion from the azimuthal equidistant by **13.4%** —
these coordinates sit on no historical projection.

---

## Layout

```
mapopt/
  family.py       parametric projection family + grid + precomputed basis
  metrics.py      Tissot shape/area + global pairwise distance metric
  classics.py     correctly-formulated baseline projections (incl. az. equidistant)
  optimize.py     Adam hill-climber, 3-objective simplex sweep, targeted solve
run_research.py   the pipeline -> results/results.json
export_web.py     results + coastlines -> web/appdata.js
web/index.html    ternary explorer (shape/area/distance) + projection-election poll
results/          machine-written outputs (checked in for convenience)
```

Coastlines: [Natural Earth](https://www.naturalearthdata.com/) (110m, public
domain). Method in the distortion-metric tradition of Airy, Kavrayskiy, Tissot,
and Goldberg & Gott.
