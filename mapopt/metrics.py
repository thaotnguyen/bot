"""Distortion metrics based on Tissot's indicatrix.

For a projection (lat, lon) -> (x, y) on the unit sphere, the local linear map
from an *orthonormal* frame on the sphere to the plane is

    A = [[ x_lat,  x_lon / cos(lat) ],
         [ y_lat,  y_lon / cos(lat) ]]

(the parallel direction is scaled by 1/cos(lat) because a unit step east spans
cos(lat) of longitude).  The singular values a >= b >= 0 of A are Tissot's two
principal scale factors.  From them:

    area scale  = a * b = |det A|
    conformality error   ln(a / b)   (0  <=>  angles preserved)
    equal-area error     ln(a * b)   (constant  <=>  areas preserved)
    max angular distortion  omega = 2 * asin((a - b) / (a + b))

We aggregate two scale-invariant, area-weighted global errors:

    eps_shape = < ln(a/b)^2 >                    (mean squared conformality err)
    eps_area  = Var( ln(a*b) ) = < (ln(a*b) - <ln(a*b)>)^2 >   (log-area variance)

eps_area is a *variance*, so it ignores the free global scale of the map -- an
equal-area projection scores 0 regardless of how big it is drawn.  These two
numbers are the axes of the Pareto frontier.
"""

from __future__ import annotations

import math
import random
from typing import Callable, Dict, List, Sequence, Tuple

Proj = Callable[[float, float], Tuple[float, float]]

_TINY = 1e-15


def _svd2(xlat: float, xlon: float, ylat: float, ylon: float, cos: float):
    """Return (s1sq, s2sq, det) for the Tissot matrix at one point.

    s1sq >= s2sq are the squared singular values (a^2, b^2); det is the signed
    area scale a*b (negative => the map has folded / flipped orientation).
    """
    p = xlon / cos
    q = ylon / cos
    E = xlat * xlat + ylat * ylat
    G = p * p + q * q
    F = xlat * p + ylat * q
    tmp = E + G
    disc = math.sqrt(max(0.0, (E - G) * (E - G) + 4.0 * F * F))
    s1sq = 0.5 * (tmp + disc)
    s2sq = max(_TINY, 0.5 * (tmp - disc))
    # Right-handed orientation: a non-folded map keeps this positive.
    # (columns are parallel-then-meridian so plate carree gives +1/cos(lat)).
    det = (p * ylat - xlat * q)  # signed a*b
    return s1sq, s2sq, det


def scores_from_partials(xlat: Sequence[float], xlon: Sequence[float],
                         ylat: Sequence[float], ylon: Sequence[float],
                         cos: Sequence[float], weights: Sequence[float],
                         fold_tau: float = 0.05) -> Dict[str, float]:
    """Aggregate global distortion scores from precomputed partials."""
    n = len(xlat)
    ln_ratio2 = 0.0        # sum w * ln(a/b)^2
    ln_area = [0.0] * n
    dets = [0.0] * n
    omega2 = 0.0
    mean_abs_det = 0.0
    for i in range(n):
        s1sq, s2sq, det = _svd2(xlat[i], xlon[i], ylat[i], ylon[i], cos[i])
        w = weights[i]
        lr = 0.5 * math.log(s1sq / s2sq)         # ln(a/b) >= 0
        ln_ratio2 += w * lr * lr
        ln_area[i] = 0.5 * math.log(s1sq * s2sq)  # ln(a*b) == ln|det|
        dets[i] = det
        mean_abs_det += w * abs(det)
        a = math.sqrt(s1sq)
        b = math.sqrt(s2sq)
        om = 2.0 * math.asin(min(1.0, (a - b) / (a + b)))
        omega2 += w * om * om

    mean_la = 0.0
    for i in range(n):
        mean_la += weights[i] * ln_area[i]
    var_la = 0.0
    for i in range(n):
        d = ln_area[i] - mean_la
        var_la += weights[i] * d * d

    # scale-invariant fold barrier: penalize det that dips below a small
    # fraction of the map's typical |det| (0 => orientation flip / overlap)
    m = max(mean_abs_det, _TINY)
    thr = fold_tau * m
    fold_pen = 0.0
    n_fold = 0
    for i in range(n):
        if dets[i] <= thr:
            e = (thr - dets[i]) / m
            fold_pen += weights[i] * e * e
        if dets[i] <= 0.0:
            n_fold += 1

    return {
        "eps_shape": ln_ratio2,
        "eps_area": var_la,
        "mean_ln_area": mean_la,
        "mean_abs_det": mean_abs_det,
        "fold_pen": fold_pen,
        "n_fold": float(n_fold),
        "rms_angular_deg": math.degrees(math.sqrt(max(0.0, omega2))),
        "rms_area_pct": 100.0 * math.sqrt(max(0.0, var_la)),
        "combined": ln_ratio2 + var_la,
    }


# ---- fast path for the parametric family -----------------------------------

def family_scores(params: Sequence[float], grid, weights=None) -> Dict[str, float]:
    if weights is None:
        weights = grid.area_w
    xlat, xlon, ylat, ylon = grid.partials(params)
    return scores_from_partials(xlat, xlon, ylat, ylon, grid.cos, weights)


def family_objective(params: Sequence[float], grid, weights, t: float,
                     fold_lambda: float = 50.0) -> float:
    """Weighted Pareto objective: (1-t)*shape + t*area + fold barrier.

    t in [0,1] slides from a shape-optimal map (t->0) to an equal-area map
    (t->1).  The fold barrier keeps the projection injective (no overlaps).
    """
    s = family_scores(params, grid, weights)
    return (1.0 - t) * s["eps_shape"] + t * s["eps_area"] + fold_lambda * s["fold_pen"]


# ---- generic path for arbitrary projections (finite differences) -----------

def generic_partials(proj: Proj, lat: float, lon: float, h: float = 1e-5):
    x1, y1 = proj(lat + h, lon)
    x2, y2 = proj(lat - h, lon)
    x3, y3 = proj(lat, lon + h)
    x4, y4 = proj(lat, lon - h)
    xlat = (x1 - x2) / (2 * h)
    ylat = (y1 - y2) / (2 * h)
    xlon = (x3 - x4) / (2 * h)
    ylon = (y3 - y4) / (2 * h)
    return xlat, xlon, ylat, ylon


def generic_scores(proj: Proj, grid, weights=None, h: float = 1e-5) -> Dict[str, float]:
    if weights is None:
        weights = grid.area_w
    xlat: List[float] = []
    xlon: List[float] = []
    ylat: List[float] = []
    ylon: List[float] = []
    for i in range(grid.n):
        a, b, c, d = generic_partials(proj, grid.lats[i], grid.lons[i], h)
        xlat.append(a)
        xlon.append(b)
        ylat.append(c)
        ylon.append(d)
    return scores_from_partials(xlat, xlon, ylat, ylon, grid.cos, weights)


# ---------------------------------------------------------------------------
# The THIRD axis: distance error (global, not local)
# ---------------------------------------------------------------------------
# Shape and area are *local* (differential) properties. Distance is *global*:
# does the straight-line distance between two points on the map match their true
# great-circle distance on the globe? No projection can hold all pairwise
# distances (that would be an isometry), so we measure how far the map is from a
# single consistent ruler:
#
#     eps_dist = Var( ln( d_map(P,Q) / d_globe(P,Q) ) )   over sampled pairs P,Q
#
# 0  <=>  every pairwise map distance equals the true distance times one constant
# scale (a perfect ruler). Like eps_area it is a log-variance, so it ignores the
# map's overall size. This is the classic third leg of the projection trilemma
# (cf. Goldberg & Gott's distance error term).

class DistanceSampler:
    """A fixed, reproducible set of point pairs with precomputed geodesic
    distances, for measuring global distance distortion."""

    def __init__(self, n_anchor: int = 120, seed: int = 12345,
                 min_sep_deg: float = 6.0, lat_limit_deg: float = 85.0):
        rng = random.Random(seed)
        smax = math.sin(math.radians(lat_limit_deg))
        self.lat: List[float] = []
        self.lon: List[float] = []
        for _ in range(n_anchor):
            self.lat.append(math.asin(rng.uniform(-smax, smax)))   # equal-area in lat
            self.lon.append(rng.uniform(-math.pi, math.pi))
        self.pairs: List[Tuple[int, int]] = []
        self.ln_dg: List[float] = []
        min_sep = math.radians(min_sep_deg)
        for i in range(n_anchor):
            for j in range(i + 1, n_anchor):
                c = (math.sin(self.lat[i]) * math.sin(self.lat[j]) +
                     math.cos(self.lat[i]) * math.cos(self.lat[j]) *
                     math.cos(self.lon[i] - self.lon[j]))
                d = math.acos(max(-1.0, min(1.0, c)))
                if d > min_sep:
                    self.pairs.append((i, j))
                    self.ln_dg.append(math.log(d))
        self.npairs = len(self.pairs)

    def score(self, proj: Proj) -> float:
        """Var(ln(d_map/d_globe)) over the pair set for a projection."""
        n = len(self.lat)
        xs = [0.0] * n
        ys = [0.0] * n
        for k in range(n):
            xy = proj(self.lat[k], self.lon[k])
            xs[k] = xy[0]
            ys[k] = xy[1]
        pr = self.pairs
        lg = self.ln_dg
        m = len(pr)
        ratios = [0.0] * m
        mean = 0.0
        for t in range(m):
            i, j = pr[t]
            dm = math.hypot(xs[i] - xs[j], ys[i] - ys[j])
            r = math.log(max(_TINY, dm)) - lg[t]
            ratios[t] = r
            mean += r
        mean /= m
        var = 0.0
        for t in range(m):
            d = ratios[t] - mean
            var += d * d
        return var / m


def family_all(params: Sequence[float], grid, dsamp: DistanceSampler,
               weights=None) -> Dict[str, float]:
    """All three raw objectives for a family member: shape, area, dist (+fold)."""
    from . import family
    s = family_scores(params, grid, weights)
    a, b = family.split(params)
    proj = lambda la, lo: family.forward_ab(a, b, la, lo)
    s["eps_dist"] = dsamp.score(proj)
    return s


def generic_all(proj: Proj, grid, dsamp: DistanceSampler, weights=None,
                h: float = 1e-5) -> Dict[str, float]:
    s = generic_scores(proj, grid, weights, h)
    s["eps_dist"] = dsamp.score(proj)
    return s


def famscore(fam, params, grid, dsamp, weights=None) -> Dict[str, float]:
    """Full 3-objective score for any Family object (analytic local + pairwise dist)."""
    if weights is None:
        weights = grid.area_w
    xlat, xlon, ylat, ylon = fam.partials(params, grid)
    s = scores_from_partials(xlat, xlon, ylat, ylon, grid.cos, weights)
    s["eps_dist"] = dsamp.score(lambda la, lo: fam.forward(params, la, lo))
    return s
