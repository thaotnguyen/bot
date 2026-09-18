"""Correctly-formulated classic projections, used as honest baselines.

Each function maps (lat, lon) in *radians* to plane (x, y) on the unit sphere.
Overall scale is irrelevant to our scale-invariant metrics, so no attempt is
made to normalize sizes across projections.

The equal-area members (Lambert CEA, sinusoidal, Mollweide, Hammer, Eckert IV,
Equal Earth) must score eps_area ~ 0 and the conformal member (Mercator) must
score eps_shape ~ 0 -- that is how we validate the distortion code.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, Tuple

PI = math.pi
SQRT2 = math.sqrt(2.0)
SQRT3 = math.sqrt(3.0)


def _clamp(x: float, lo: float = -1.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


def _sinc_factor(alpha: float) -> float:
    """Return alpha / sin(alpha), the reciprocal of sinc, with a safe limit."""
    if abs(alpha) < 1e-8:
        return 1.0
    return alpha / math.sin(alpha)


def equirectangular(lat: float, lon: float) -> Tuple[float, float]:
    return lon, lat


def mercator(lat: float, lon: float) -> Tuple[float, float]:
    return lon, math.asinh(math.tan(lat))  # = ln(tan(pi/4 + lat/2))


def lambert_cea(lat: float, lon: float) -> Tuple[float, float]:
    # Lambert cylindrical equal-area
    return lon, math.sin(lat)


def sinusoidal(lat: float, lon: float) -> Tuple[float, float]:
    return lon * math.cos(lat), lat


def mollweide(lat: float, lon: float) -> Tuple[float, float]:
    if abs(abs(lat) - PI / 2) < 1e-9:
        theta = math.copysign(PI / 2, lat)
    else:
        theta = lat
        for _ in range(30):
            f = 2 * theta + math.sin(2 * theta) - PI * math.sin(lat)
            fp = 2 + 2 * math.cos(2 * theta)
            if abs(fp) < 1e-12:
                break
            step = f / fp
            theta -= step
            if abs(step) < 1e-12:
                break
    x = (2 * SQRT2 / PI) * lon * math.cos(theta)
    y = SQRT2 * math.sin(theta)
    return x, y


def hammer(lat: float, lon: float) -> Tuple[float, float]:
    d = math.sqrt(1 + math.cos(lat) * math.cos(lon / 2))
    x = 2 * SQRT2 * math.cos(lat) * math.sin(lon / 2) / d
    y = SQRT2 * math.sin(lat) / d
    return x, y


def aitoff(lat: float, lon: float) -> Tuple[float, float]:
    alpha = math.acos(_clamp(math.cos(lat) * math.cos(lon / 2)))
    fac = _sinc_factor(alpha)
    x = 2 * math.cos(lat) * math.sin(lon / 2) * fac
    y = math.sin(lat) * fac
    return x, y


_PHI1 = math.acos(2.0 / PI)  # Winkel Tripel standard parallel (~50.46 deg)


def winkel_tripel(lat: float, lon: float) -> Tuple[float, float]:
    alpha = math.acos(_clamp(math.cos(lat) * math.cos(lon / 2)))
    fac = _sinc_factor(alpha)
    x = 0.5 * (lon * math.cos(_PHI1) + 2 * math.cos(lat) * math.sin(lon / 2) * fac)
    y = 0.5 * (lat + math.sin(lat) * fac)
    return x, y


def eckert4(lat: float, lon: float) -> Tuple[float, float]:
    theta = lat / 2.0
    rhs = (2 + PI / 2) * math.sin(lat)
    for _ in range(30):
        f = theta + math.sin(theta) * math.cos(theta) + 2 * math.sin(theta) - rhs
        fp = 1 + math.cos(2 * theta) + 2 * math.cos(theta)
        if abs(fp) < 1e-12:
            break
        step = f / fp
        theta -= step
        if abs(step) < 1e-12:
            break
    x = 2 / math.sqrt(PI * (4 + PI)) * lon * (1 + math.cos(theta))
    y = 2 * math.sqrt(PI / (4 + PI)) * math.sin(theta)
    return x, y


_EE = (1.340264, -0.081106, 0.000893, 0.003796)


def equal_earth(lat: float, lon: float) -> Tuple[float, float]:
    a1, a2, a3, a4 = _EE
    th = math.asin(_clamp(SQRT3 / 2 * math.sin(lat)))
    den = 9 * a4 * th ** 8 + 7 * a3 * th ** 6 + 3 * a2 * th ** 2 + a1
    x = 2 * SQRT3 * lon * math.cos(th) / (3 * den)
    y = a1 * th + a2 * th ** 3 + a3 * th ** 7 + a4 * th ** 9
    return x, y


# Robinson is defined by a lookup table (interpolated), not a closed form.
_ROB_X = [1.0000, 0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600, 0.9427,
          0.9216, 0.8962, 0.8679, 0.8350, 0.7986, 0.7597, 0.7186, 0.6732,
          0.6213, 0.5722, 0.5322]
_ROB_Y = [0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720, 0.4340,
          0.4958, 0.5571, 0.6176, 0.6769, 0.7346, 0.7903, 0.8435, 0.8936,
          0.9394, 0.9761, 1.0000]


def robinson(lat: float, lon: float) -> Tuple[float, float]:
    ad = abs(math.degrees(lat)) / 5.0
    i = int(math.floor(ad))
    if i >= 18:
        px, py = _ROB_X[18], _ROB_Y[18]
    else:
        frac = ad - i
        px = _ROB_X[i] + frac * (_ROB_X[i + 1] - _ROB_X[i])
        py = _ROB_Y[i] + frac * (_ROB_Y[i + 1] - _ROB_Y[i])
    x = 0.8487 * px * lon
    y = 1.3523 * py * (1.0 if lat >= 0 else -1.0)
    return x, y


CLASSICS: Dict[str, Callable[[float, float], Tuple[float, float]]] = {
    "Equirectangular": equirectangular,
    "Mercator": mercator,
    "Lambert CEA": lambert_cea,
    "Sinusoidal": sinusoidal,
    "Mollweide": mollweide,
    "Hammer": hammer,
    "Aitoff": aitoff,
    "Winkel Tripel": winkel_tripel,
    "Robinson": robinson,
    "Eckert IV": eckert4,
    "Equal Earth": equal_earth,
}

# Whether each classic is theoretically equal-area / conformal (for validation).
EQUAL_AREA = {"Lambert CEA", "Sinusoidal", "Mollweide", "Hammer", "Eckert IV", "Equal Earth"}
CONFORMAL = {"Mercator"}
