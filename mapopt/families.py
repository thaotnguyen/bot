"""Projection *families*, each with a characteristic OUTLINE, so the search can
range over the boundary shape/topology itself — not just interior coefficients.

Every family exposes a uniform interface:
    label, outline            : names
    p0()                      : initial parameter list
    bounds                    : per-parameter (lo, hi) clamps (None = free)
    forward(p, lat, lon)      : (x, y)
    partials(p, grid)         : analytic (xlat, xlon, ylat, ylon) over the grid
    describe(p)               : short human summary (e.g. the superellipse exponent)

The metrics only need `partials` (local shape/area) and `forward` (global
distance), so a single generic optimizer drives all families identically.

Outlines covered:
    Cylindrical        -> rectangle
    Pseudocylindrical  -> pointed oval / lens (straight parallels)
    Lenticular         -> rounded lens (curved parallels)   [the polynomial family]
    Azimuthal          -> disc (radial)
    Superellipse       -> boundary morphs rectangle <-> ellipse <-> diamond via p
"""

from __future__ import annotations

import math
from typing import List, Sequence, Tuple

from . import family as _poly

PI = math.pi
HP = PI / 2.0
TWO_OVER_PI = 2.0 / PI


class Cylindrical:
    label = "Cylindrical"
    outline = "rectangle"

    def p0(self):
        return [1.0, HP, 0.0, 0.0]           # A, then y = b1 t + b2 t^3 + b3 t^5
    bounds = [(0.05, 20), None, None, None]

    def forward(self, p, lat, lon):
        t = lat / HP
        A, b1, b2, b3 = p
        return A * lon, b1 * t + b2 * t ** 3 + b3 * t ** 5

    def partials(self, p, g):
        A, b1, b2, b3 = p
        xlat = [0.0] * g.n
        xlon = [A] * g.n
        ylat = []
        ylon = [0.0] * g.n
        for i in range(g.n):
            t = g.lats[i] / HP
            ylat.append((b1 + 3 * b2 * t * t + 5 * b3 * t ** 4) * TWO_OVER_PI)
        return xlat, xlon, ylat, ylon

    def describe(self, p):
        return "rectangular"


class Pseudocylindrical:
    label = "Pseudocylindrical"
    outline = "oval / lens"

    def p0(self):
        return [1.0, 0.0, 0.0, HP, 0.0, 0.0]   # h = a0+a1 t^2+a2 t^4 ; y = b1 t+b2 t^3+b3 t^5
    bounds = [(0.02, 5), None, None, None, None, None]

    def forward(self, p, lat, lon):
        t = lat / HP
        a0, a1, a2, b1, b2, b3 = p
        h = a0 + a1 * t * t + a2 * t ** 4
        return lon * h, b1 * t + b2 * t ** 3 + b3 * t ** 5

    def partials(self, p, g):
        a0, a1, a2, b1, b2, b3 = p
        xlat = []; xlon = []; ylat = []; ylon = [0.0] * g.n
        for i in range(g.n):
            t = g.lats[i] / HP
            lon = g.lons[i]
            h = a0 + a1 * t * t + a2 * t ** 4
            hp = (2 * a1 * t + 4 * a2 * t ** 3) * TWO_OVER_PI
            xlon.append(h)
            xlat.append(lon * hp)
            ylat.append((b1 + 3 * b2 * t * t + 5 * b3 * t ** 4) * TWO_OVER_PI)
        return xlat, xlon, ylat, ylon

    def describe(self, p):
        h_pole = p[0] + p[1] + p[2]
        return "pointed pole" if h_pole < 0.08 else "pole line"


class Azimuthal:
    label = "Azimuthal"
    outline = "disc"

    def p0(self):
        return [1.0, 0.0, 0.0]               # R(c) = d1 c + d2 c^2 + d3 c^3, c = colatitude
    bounds = [(0.02, 5), None, None]

    def forward(self, p, lat, lon):
        d1, d2, d3 = p
        c = HP - lat
        R = d1 * c + d2 * c * c + d3 * c ** 3
        return R * math.sin(lon), -R * math.cos(lon)

    def partials(self, p, g):
        d1, d2, d3 = p
        xlat = []; xlon = []; ylat = []; ylon = []
        for i in range(g.n):
            lon = g.lons[i]
            c = HP - g.lats[i]
            R = d1 * c + d2 * c * c + d3 * c ** 3
            Rp = d1 + 2 * d2 * c + 3 * d3 * c * c        # dR/dc ; dc/dlat = -1
            s, co = math.sin(lon), math.cos(lon)
            xlat.append(-Rp * s)      # d/dlat (R sin) = R'(dc/dlat) sin = -R' sin
            xlon.append(R * co)
            ylat.append(Rp * co)      # y = -R cos ; d/dlat = -R'(dc/dlat) cos = R' cos
            ylon.append(R * s)
        return xlat, xlon, ylat, ylon

    def describe(self, p):
        return "polar disc"


class Lenticular:
    """Adapter around the 18-parameter symmetric polynomial family."""
    label = "Lenticular"
    outline = "rounded lens"

    def p0(self):
        return _poly.init_equirectangular()
    bounds = None

    def forward(self, p, lat, lon):
        return _poly.forward(p, lat, lon)

    def partials(self, p, g):
        return g.partials(p)

    def describe(self, p):
        return "curved-parallel lens"


class Superellipse:
    """Outline is a superellipse |x/A|^p + |y|^p = 1; the optimizer chooses p,
    so the boundary morphs rectangle (p large) <-> ellipse (p=2) <-> diamond (p=1)."""
    label = "Superellipse"
    outline = "tunable (rect<->ellipse<->diamond)"

    def p0(self):
        return [2.0, 2.0, 1.0, HP, 0.0, 0.0]   # p, A(width), B(height), gy = b1 t+b2 t^3+b3 t^5
    bounds = [(1.1, 8.0), (0.05, 20), (0.05, 20), None, None, None]

    def _nmax(self, p):
        nmax = p[3] + p[4] + p[5]
        return nmax if abs(nmax) > 1e-9 else 1e-9

    def forward(self, p, lat, lon):
        pw, A, B, b1, b2, b3 = p
        t = lat / HP
        eta = (b1 * t + b2 * t ** 3 + b3 * t ** 5) / self._nmax(p)
        eta = max(-0.999999, min(0.999999, eta))
        W = max(1e-9, 1.0 - abs(eta) ** pw)
        return (lon / PI) * A * W ** (1.0 / pw), B * eta

    def partials(self, p, g):
        pw, A, B, b1, b2, b3 = p
        nmax = self._nmax(p)
        xlat = []; xlon = []; ylat = []; ylon = [0.0] * g.n
        for i in range(g.n):
            t = g.lats[i] / HP
            lon = g.lons[i]
            eta = (b1 * t + b2 * t ** 3 + b3 * t ** 5) / nmax
            eta = max(-0.999999, min(0.999999, eta))
            detadlat = ((b1 + 3 * b2 * t * t + 5 * b3 * t ** 4) / nmax) * TWO_OVER_PI
            ae = abs(eta)
            W = max(1e-9, 1.0 - ae ** pw)
            xlon.append((1.0 / PI) * A * W ** (1.0 / pw))
            dWinv_deta = -(ae ** (pw - 1)) * (1.0 if eta >= 0 else -1.0) * W ** (1.0 / pw - 1.0)
            xlat.append((lon / PI) * A * dWinv_deta * detadlat)
            ylat.append(B * detadlat)
        return xlat, xlon, ylat, ylon

    def describe(self, p):
        pw = p[0]
        kind = "diamond-ish" if pw < 1.5 else ("ellipse" if pw < 2.6 else ("rounded-rect" if pw < 5 else "rectangle"))
        return "p=%.2f (%s)" % (pw, kind)


ALL_FAMILIES = [Cylindrical(), Pseudocylindrical(), Lenticular(), Azimuthal(), Superellipse()]
BAKE_OFF = [Cylindrical(), Pseudocylindrical(), Lenticular(), Azimuthal()]
