"""An approximate human-population weighting field.

Standard projection error is weighted by *area* -- which spends most of its
"distortion budget" on empty ocean and ice.  A more human-centered question is:
where do the errors fall relative to where people actually live?  We build a
smooth population-density field from a dozen major population centers plus a
latitude profile that reflects the well-known concentration of humanity in the
northern mid-latitudes, then weight each grid cell by (density x cell-area).

This is a deliberately transparent *approximation* (not a census raster): it is
built from public, well-known demographic facts and is easy to swap for a real
gridded population product (e.g. SEDAC GPW) by replacing `density()`.  Its
purpose is to demonstrate that optimizing for people rather than for acreage
yields a measurably different projection.
"""

from __future__ import annotations

import math
from typing import List, Tuple

# (name, lat, lon, weight ~ billions, angular sigma in degrees)
CENTERS: List[Tuple[str, float, float, float, float]] = [
    ("South Asia",        24.0,  80.0, 1.85, 12.0),
    ("East Asia",         32.0, 114.0, 1.45, 12.0),
    ("Southeast Asia",     5.0, 108.0, 0.70, 12.0),
    ("Europe",            48.0,  12.0, 0.60, 12.0),
    ("Eastern N. America",38.0, -83.0, 0.32, 12.0),
    ("West Africa",        9.0,   6.0, 0.45, 10.0),
    ("Middle East",       33.0,  42.0, 0.45, 12.0),
    ("Mexico/Cent. Am.",  20.0, -99.0, 0.24, 9.0),
    ("SE South America",  -23.0,-47.0, 0.32, 11.0),
    ("Japan/Korea",       37.0, 132.0, 0.21, 6.0),
    ("East Africa",        5.0,  38.0, 0.30, 10.0),
    ("Southern Africa",  -26.0,  28.0, 0.10, 9.0),
    ("Australia (SE)",   -34.0, 148.0, 0.03, 6.0),
]


def _ang_dist_deg(lat1, lon1, lat2, lon2):
    """Great-circle distance in degrees between two lat/lon points (degrees)."""
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon1 - lon2)
    c = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    c = max(-1.0, min(1.0, c))
    return math.degrees(math.acos(c))


def density(lat_deg: float, lon_deg: float) -> float:
    """Approximate relative population density at a location (degrees)."""
    d = 0.0
    for _name, clat, clon, wt, sig in CENTERS:
        r = _ang_dist_deg(lat_deg, lon_deg, clat, clon)
        d += wt * math.exp(-0.5 * (r / sig) ** 2)
    # small uniform floor so no region is exactly zero-weighted
    return d + 0.002


def weights_for_grid(grid) -> List[float]:
    """Population weights per grid point = density * cell-area, normalized."""
    w = []
    for i in range(grid.n):
        w.append(density(grid.lat_deg[i], grid.lon_deg[i]) * grid.cos[i])
    s = sum(w)
    return [x / s for x in w]
