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
