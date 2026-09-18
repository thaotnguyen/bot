"""A flexible, symmetry-constrained family of world map projections.

We represent a projection as a pair of *polynomials* in normalized
coordinates::

    u = lon / pi          in [-1, 1]
    v = lat / (pi/2)       in [-1, 1]

with the parities that every sensible whole-world "lenticular" projection
shares (symmetric about the equator and the central meridian):

    x(u, v) = sum_{i,j} a_ij * u^(2i+1) * v^(2j)      # odd in u, even in v
    y(u, v) = sum_{i,j} b_ij * u^(2i)   * v^(2j+1)    # even in u, odd in v

With i, j in {0, 1, 2} this is 9 coefficients each (18 total).  The family is
deliberately *not* any named projection: equirectangular is the single point
(a_00 = pi, b_00 = pi/2, everything else 0), and the classics (Aitoff, Winkel
Tripel, Mollweide, ...) are transcendental and do not sit exactly on any
polynomial coefficient vector.  Hill-climbing these 18 numbers therefore
explores a continuum of genuinely new projections while keeping the map
smooth, symmetric, and (with a fold penalty) non-self-overlapping.

Because x and y are *linear* in the coefficients, all partial derivatives that
the distortion metric needs are also linear, so we precompute a constant basis
per grid point and every objective evaluation is just a handful of dot
products -- fast enough to hill-climb an entire Pareto frontier in pure Python.
"""

from __future__ import annotations

import math
from typing import List, Sequence, Tuple

PI = math.pi
HALF_PI = math.pi / 2.0

# Coefficient index sets.  Degree D_UV controls flexibility; 2 -> {0,1,2}.
D_UV = 2
X_TERMS: List[Tuple[int, int]] = [(i, j) for i in range(D_UV + 1) for j in range(D_UV + 1)]
Y_TERMS: List[Tuple[int, int]] = [(i, j) for i in range(D_UV + 1) for j in range(D_UV + 1)]
NX = len(X_TERMS)
NY = len(Y_TERMS)
NPARAMS = NX + NY


def init_equirectangular() -> List[float]:
    """Parameter vector for the plate carree (equirectangular) map.

    x = lon, y = lat.  This is the hill-climber's starting point: a valid,
    familiar map that is neither conformal nor equal-area.
    """
    a = [0.0] * NX
    b = [0.0] * NY
    # x-term (i=0, j=0) is u^1 v^0 -> a_00 * u = a_00 * lon/pi ; want = lon.
    a[X_TERMS.index((0, 0))] = PI
    # y-term (i=0, j=0) is u^0 v^1 -> b_00 * v = b_00 * lat/(pi/2) ; want = lat.
    b[Y_TERMS.index((0, 0))] = HALF_PI
    return a + b


def split(params: Sequence[float]) -> Tuple[List[float], List[float]]:
    return list(params[:NX]), list(params[NX:])


def forward(params: Sequence[float], lat: float, lon: float) -> Tuple[float, float]:
    """Project a single (lat, lon) in radians to plane coordinates (x, y)."""
    a, b = split(params)
    u = lon / PI
    v = lat / HALF_PI
    x = 0.0
    for k, (i, j) in enumerate(X_TERMS):
        x += a[k] * (u ** (2 * i + 1)) * (v ** (2 * j))
    y = 0.0
    for k, (i, j) in enumerate(Y_TERMS):
        y += b[k] * (u ** (2 * i)) * (v ** (2 * j + 1))
    return x, y


# ---------------------------------------------------------------------------
# Grid + precomputed basis
# ---------------------------------------------------------------------------

class Grid:
    """A lat/lon sample grid with precomputed, parameter-independent basis
    vectors for the four partial derivatives the distortion metric needs.

    For each sample point we store the (constant) vectors so that, given a
    coefficient vector, the partials are simple dot products::

        x_lon = (1/pi)   * dot(a, bxu)
        x_lat = (2/pi)   * dot(a, bxv)
        y_lon = (1/pi)   * dot(b, byu)
        y_lat = (2/pi)   * dot(b, byv)
    """

    def __init__(self, lat_step: float = 3.0, lon_step: float = 5.0,
                 lat_limit: float = 88.5):
        self.lats: List[float] = []
        self.lons: List[float] = []
        self.lat_deg: List[float] = []
        self.lon_deg: List[float] = []
        self.cos: List[float] = []          # cos(lat)
        self.area_w: List[float] = []       # area element weight (normalized)

        # basis matrices, one row (length NX or NY) per sample point
        self.bxu: List[List[float]] = []
        self.bxv: List[List[float]] = []
        self.byu: List[List[float]] = []
        self.byv: List[List[float]] = []

        lat = -lat_limit
        raw_w: List[float] = []
        while lat <= lat_limit + 1e-9:
            lon = -180.0 + lon_step / 2.0
            while lon <= 180.0 - lon_step / 2.0 + 1e-9:
                latr = math.radians(lat)
                lonr = math.radians(lon)
                c = math.cos(latr)
                u = lonr / PI
                v = latr / HALF_PI

                bxu = [0.0] * NX
                bxv = [0.0] * NX
                for k, (i, j) in enumerate(X_TERMS):
                    pu, pv = 2 * i + 1, 2 * j
                    # d/du of u^pu v^pv
                    bxu[k] = pu * (u ** (pu - 1)) * (v ** pv)
                    # d/dv of u^pu v^pv  (zero when pv == 0)
                    bxv[k] = 0.0 if pv == 0 else pv * (u ** pu) * (v ** (pv - 1))
                byu = [0.0] * NY
                byv = [0.0] * NY
                for k, (i, j) in enumerate(Y_TERMS):
                    pu, pv = 2 * i, 2 * j + 1
                    byu[k] = 0.0 if pu == 0 else pu * (u ** (pu - 1)) * (v ** pv)
                    byv[k] = pv * (u ** pu) * (v ** (pv - 1))

                self.lats.append(latr)
                self.lons.append(lonr)
                self.lat_deg.append(lat)
                self.lon_deg.append(lon)
                self.cos.append(c)
                raw_w.append(c)  # area element ~ cos(lat)
                self.bxu.append(bxu)
                self.bxv.append(bxv)
                self.byu.append(byu)
                self.byv.append(byv)
                lon += lon_step
            lat += lat_step

        s = sum(raw_w)
        self.area_w = [w / s for w in raw_w]
        self.n = len(self.lats)

    def partials(self, params: Sequence[float]):
        """Return four lists (x_lat, x_lon, y_lat, y_lon) over the grid."""
        a, b = split(params)
        inv_pi = 1.0 / PI
        two_pi = 2.0 / PI
        xlat: List[float] = []
        xlon: List[float] = []
        ylat: List[float] = []
        ylon: List[float] = []
        bxu, bxv, byu, byv = self.bxu, self.bxv, self.byu, self.byv
        for idx in range(self.n):
            ru = bxu[idx]
            rv = bxv[idx]
            su = byu[idx]
            sv = byv[idx]
            dxu = 0.0
            dxv = 0.0
            for k in range(NX):
                dxu += a[k] * ru[k]
                dxv += a[k] * rv[k]
            dyu = 0.0
            dyv = 0.0
            for k in range(NY):
                dyu += b[k] * su[k]
                dyv += b[k] * sv[k]
            xlon.append(inv_pi * dxu)
            xlat.append(two_pi * dxv)
            ylon.append(inv_pi * dyu)
            ylat.append(two_pi * dyv)
        return xlat, xlon, ylat, ylon
