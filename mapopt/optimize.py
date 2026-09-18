"""Hill-climbing (Adam gradient descent) over the projection family.

The objective is scale-invariant, so between steps we "gauge-fix" the
coefficient vector to keep the mean log-area at 0 (a map of unit average
scale); this keeps the fold barrier meaningful and the numbers well
conditioned.  Gradients use forward finite differences -- cheap and robust
given the objective is smooth and only 18-dimensional.

A Pareto sweep walks the weight t from shape-optimal (t->0) to area-optimal
(t->1) using warm starts (continuation), so each subsequent optimum is a small
move from the last and converges in a few dozen steps.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, List, Sequence, Tuple

from . import metrics
from .family import Grid, init_equirectangular, NPARAMS


def gauge_fix(params: List[float], grid: Grid, weights) -> List[float]:
    """Rescale so the area-weighted mean log-area is 0 (unit average scale)."""
    s = metrics.family_scores(params, grid, weights)
    factor = math.exp(-0.5 * s["mean_ln_area"])
    return [p * factor for p in params]


def _grad(f: Callable[[List[float]], float], x: List[float],
          base: float, eps: float = 1e-5) -> List[float]:
    g = [0.0] * len(x)
    for i in range(len(x)):
        old = x[i]
        x[i] = old + eps
        g[i] = (f(x) - base) / eps
        x[i] = old
    return g


def optimize_single(grid: Grid, weights, t: float, x0: Sequence[float],
                    steps: int = 120, lr: float = 0.02,
                    fold_lambda: float = 50.0) -> Tuple[List[float], Dict[str, float]]:
    """Minimize (1-t)*eps_shape + t*eps_area (+ fold barrier) from x0."""
    x = list(x0)
    m = [0.0] * len(x)
    v = [0.0] * len(x)
    b1, b2, eps_a = 0.9, 0.999, 1e-8

    def f(p):
        return metrics.family_objective(p, grid, weights, t, fold_lambda)

    for step in range(1, steps + 1):
        x = gauge_fix(x, grid, weights)
        base = f(x)
        g = _grad(f, x, base)
        for i in range(len(x)):
            m[i] = b1 * m[i] + (1 - b1) * g[i]
            v[i] = b2 * v[i] + (1 - b2) * g[i] * g[i]
            mhat = m[i] / (1 - b1 ** step)
            vhat = v[i] / (1 - b2 ** step)
            x[i] -= lr * mhat / (math.sqrt(vhat) + eps_a)

    x = gauge_fix(x, grid, weights)
    return x, metrics.family_scores(x, grid, weights)


def pareto_sweep(grid: Grid, weights, t_list: Sequence[float],
                 steps: int = 90, lr: float = 0.02) -> List[Dict]:
    """Trace the shape/area Pareto frontier with warm-start continuation."""
    results: List[Dict] = []
    x = init_equirectangular()
    # Sort ascending and walk upward so each solve warm-starts from the last.
    order = sorted(range(len(t_list)), key=lambda i: t_list[i])
    solved: Dict[int, Dict] = {}
    for oi in order:
        t = t_list[oi]
        x, s = optimize_single(grid, weights, t, x, steps=steps, lr=lr)
        solved[oi] = {
            "t": t,
            "params": list(x),
            "eps_shape": s["eps_shape"],
            "eps_area": s["eps_area"],
            "combined": s["eps_shape"] + s["eps_area"],
            "rms_angular_deg": s["rms_angular_deg"],
            "rms_area_pct": s["rms_area_pct"],
            "n_fold": s["n_fold"],
        }
    for i in range(len(t_list)):
        results.append(solved[i])
    return results


def pareto_filter(points: List[Dict]) -> List[Dict]:
    """Keep only non-dominated points (lower eps_shape AND eps_area is better)."""
    keep = []
    for p in points:
        dominated = False
        for q in points:
            if q is p:
                continue
            if (q["eps_shape"] <= p["eps_shape"] + 1e-9 and
                    q["eps_area"] <= p["eps_area"] + 1e-9 and
                    (q["eps_shape"] < p["eps_shape"] - 1e-9 or
                     q["eps_area"] < p["eps_area"] - 1e-9)):
                dominated = True
                break
        if not dominated:
            keep.append(p)
    return keep


# ---------------------------------------------------------------------------
# Three-objective optimization: shape vs area vs distance
# ---------------------------------------------------------------------------
# Objectives are normalized by a reference map (equirectangular) so the three
# comparable-at-~1 numbers can be blended with barycentric weights that sweep the
# 2-D Pareto *surface* in 3-objective space.

def family_objective3(params, grid, dsamp, refs, w, fold_lambda: float = 50.0) -> float:
    from . import family
    s = metrics.family_scores(params, grid, grid.area_w)
    a, b = family.split(params)
    proj = lambda la, lo: family.forward_ab(a, b, la, lo)
    dist = dsamp.score(proj)
    return (w[0] * s["eps_shape"] / refs[0] +
            w[1] * s["eps_area"] / refs[1] +
            w[2] * dist / refs[2] +
            fold_lambda * s["fold_pen"])


def optimize_single3(grid, dsamp, refs, w, x0, steps: int = 70, lr: float = 0.02,
                     fold_lambda: float = 50.0):
    x = list(x0)
    m = [0.0] * len(x)
    v = [0.0] * len(x)
    b1, b2, eps_a = 0.9, 0.999, 1e-8

    def f(p):
        return family_objective3(p, grid, dsamp, refs, w, fold_lambda)

    for step in range(1, steps + 1):
        x = gauge_fix(x, grid, grid.area_w)
        base = f(x)
        g = _grad(f, x, base)
        for i in range(len(x)):
            m[i] = b1 * m[i] + (1 - b1) * g[i]
            v[i] = b2 * v[i] + (1 - b2) * g[i] * g[i]
            mhat = m[i] / (1 - b1 ** step)
            vhat = v[i] / (1 - b2 ** step)
            x[i] -= lr * mhat / (math.sqrt(vhat) + eps_a)
    x = gauge_fix(x, grid, grid.area_w)
    return x, metrics.family_all(x, grid, dsamp, grid.area_w)


def simplex_weights(n: int):
    """Barycentric weight grid: all (i,j,k)/n with i+j+k=n. n=4 -> 15 points."""
    out = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            k = n - i - j
            out.append((i / n, j / n, k / n))
    return out


def simplex_sweep(grid, dsamp, refs, weights, steps: int = 60, lr: float = 0.02):
    """Trace the shape/area/distance Pareto surface; warm-start each solve from
    the nearest already-solved weight (good continuation across the simplex)."""
    def dist_w(a, b):
        return sum((a[i] - b[i]) ** 2 for i in range(3))
    centroid = (1 / 3, 1 / 3, 1 / 3)
    order = sorted(range(len(weights)), key=lambda i: dist_w(weights[i], centroid))
    solved = {}         # index -> params
    results = [None] * len(weights)
    for idx in order:
        w = weights[idx]
        if not solved:
            x0 = init_equirectangular()
        else:
            nearest = min(solved.keys(), key=lambda j: dist_w(weights[j], w))
            x0 = solved[nearest]
        x, s = optimize_single3(grid, dsamp, refs, w, x0, steps=steps, lr=lr)
        solved[idx] = x
        results[idx] = {
            "w": list(w),
            "params": list(x),
            "eps_shape": s["eps_shape"],
            "eps_area": s["eps_area"],
            "eps_dist": s["eps_dist"],
            "rms_angular_deg": s["rms_angular_deg"],
            "n_fold": int(s["n_fold"]),
        }
    return results


def optimize_obj(f, x0, grid, steps: int = 100, lr: float = 0.02):
    """Adam on an arbitrary objective f(params)->float, gauge-fixed each step."""
    x = list(x0)
    m = [0.0] * len(x)
    v = [0.0] * len(x)
    b1, b2, eps_a = 0.9, 0.999, 1e-8
    best, best_x = float("inf"), list(x)
    for step in range(1, steps + 1):
        x = gauge_fix(x, grid, grid.area_w)
        base = f(x)
        if base < best:
            best, best_x = base, list(x)
        g = _grad(f, x, base)
        for i in range(len(x)):
            m[i] = b1 * m[i] + (1 - b1) * g[i]
            v[i] = b2 * v[i] + (1 - b2) * g[i] * g[i]
            mhat = m[i] / (1 - b1 ** step)
            vhat = v[i] / (1 - b2 ** step)
            x[i] -= lr * mhat / (math.sqrt(vhat) + eps_a)
    x = gauge_fix(x, grid, grid.area_w)
    if f(x) < best:
        best_x = list(x)
    return best_x


# ---------------------------------------------------------------------------
# Generic optimizer over any Family (outline/topology as a search dimension)
# ---------------------------------------------------------------------------

def _clamp_params(x, bounds):
    if bounds:
        for i, bd in enumerate(bounds):
            if bd is not None:
                lo, hi = bd
                if x[i] < lo: x[i] = lo
                if x[i] > hi: x[i] = hi
    return x


def optimize_family(fam, grid, dsamp, refs, w, x0=None, steps: int = 70, lr: float = 0.02,
                    fold_lambda: float = 50.0, anchor: float = 0.003):
    """Minimize w-weighted (shape,area,dist) for an arbitrary Family.

    No gauge-fix: the metrics are scale-invariant and a light anchor keeps the
    map's mean |det| ~ 1 (fixing only the uniform-scale gauge, not aspect or the
    superellipse exponent). Parameters are clamped to the family's bounds.
    """
    x = list(x0 if x0 is not None else fam.p0())
    bounds = fam.bounds
    m = [0.0] * len(x); v = [0.0] * len(x)
    b1, b2, eps_a = 0.9, 0.999, 1e-8

    def f(p):
        s = metrics.famscore(fam, p, grid, dsamp)
        anc = anchor * (math.log(max(1e-9, s["mean_abs_det"]))) ** 2
        return (w[0] * s["eps_shape"] / refs[0] + w[1] * s["eps_area"] / refs[1] +
                w[2] * s["eps_dist"] / refs[2] + fold_lambda * s["fold_pen"] + anc)

    best, best_x = float("inf"), list(x)
    for step in range(1, steps + 1):
        x = _clamp_params(x, bounds)
        base = f(x)
        if base < best:
            best, best_x = base, list(x)
        g = _grad(f, x, base)
        for i in range(len(x)):
            m[i] = b1 * m[i] + (1 - b1) * g[i]
            v[i] = b2 * v[i] + (1 - b2) * g[i] * g[i]
            mhat = m[i] / (1 - b1 ** step)
            vhat = v[i] / (1 - b2 ** step)
            x[i] -= lr * mhat / (math.sqrt(vhat) + eps_a)
    x = _clamp_params(x, bounds)
    if f(x) < best:
        best_x = list(x)
    return best_x, metrics.famscore(fam, best_x, grid, dsamp)
