#!/usr/bin/env python3
"""Build the deployed site's data bundle: a shape/area/distance objective
triangle for EVERY outline family (so both levers are explorable), plus the
champion projections for voting and the coastlines.

Run:  python3 run_site_data.py   ->  site/appdata.js
"""

from __future__ import annotations

import json, math, os

from mapopt import classics, metrics, optimize, families
from mapopt.family import Grid, X_TERMS, Y_TERMS, NX
from export_web import simplify_land, DEFAULT_LAND

HERE = os.path.dirname(os.path.abspath(__file__))
TYPE = {"Cylindrical": "cyl", "Pseudocylindrical": "pseudo", "Lenticular": "lent",
        "Azimuthal": "azim", "Superellipse": "super"}


def family_triangle(fam, OPT, OPT_D, REP, REP_D, refs, n=4, steps=55):
    """Optimize the family across the weight simplex; warm-start each solve from
    the nearest solved weight; re-score on the report grid."""
    weights = optimize.simplex_weights(n)
    def dw(a, b):
        return sum((a[i] - b[i]) ** 2 for i in range(3))
    centroid = (1/3, 1/3, 1/3)
    order = sorted(range(len(weights)), key=lambda i: dw(weights[i], centroid))
    solved = {}
    pts = [None] * len(weights)
    for idx in order:
        w = weights[idx]
        x0 = solved[min(solved, key=lambda j: dw(weights[j], w))] if solved else None
        x, _ = optimize.optimize_family(fam, OPT, OPT_D, refs, w, x0=x0, steps=steps, lr=0.02)
        solved[idx] = x
        s = metrics.famscore(fam, x, REP, REP_D)
        pts[idx] = {"w": [round(v, 4) for v in w], "params": [round(v, 6) for v in x],
                    "shape": round(s["eps_shape"], 4), "area": round(s["eps_area"], 4),
                    "dist": round(s["eps_dist"], 4), "ang": round(s["rms_angular_deg"], 2),
                    "norm": round(s["eps_shape"]/refs[0] + s["eps_area"]/refs[1] + s["eps_dist"]/refs[2], 3)}
    return pts


def main():
    OPT = Grid(lat_step=4.0, lon_step=6.0)
    REP = Grid(lat_step=3.0, lon_step=4.0)
    OPT_D = metrics.DistanceSampler(n_anchor=90, seed=7)
    REP_D = metrics.DistanceSampler(n_anchor=140, seed=101)
    eqr = metrics.generic_all(classics.equirectangular, REP, REP_D)
    refs = (eqr["eps_shape"], eqr["eps_area"], eqr["eps_dist"])
    print(f"grids opt {OPT.n} rep {REP.n}; refs shape={refs[0]:.3f} area={refs[1]:.3f} dist={refs[2]:.3f}")

    fam_tris = []
    for fam in families.ALL_FAMILIES:
        print(f"  triangle for {fam.label} ({fam.outline}) ...")
        pts = family_triangle(fam, OPT, OPT_D, REP, REP_D, refs)
        fam_tris.append({"label": fam.label, "outline": fam.outline, "type": TYPE[fam.label], "points": pts})

    # champions (for voting) from the existing 3-objective results
    r = json.load(open(os.path.join(HERE, "results", "results.json")))
    C = r["champions"]
    def champ(key):
        e = C[key]; p = e["params"]
        return {"a": [round(v, 6) for v in p[:NX]], "b": [round(v, 6) for v in p[NX:]],
                "shape": round(e["eps_shape"], 4), "area": round(e["eps_area"], 4),
                "dist": round(e["eps_dist"], 4)}
    disc = [
        {"key": "vantage", "name": "Vantage", "ours": 1, "type": "lent",
         "why": "Lower shape, area, and distance error than the Winkel Tripel.", **champ("WinkelDominator")},
        {"key": "equipoise", "name": "Equipoise", "ours": 1, "type": "lent",
         "why": "Balanced across shape, area, and distance.", **champ("Triathlon")},
        {"key": "wayfarer", "name": "Wayfarer", "ours": 1, "type": "lent",
         "why": "Lowest distance error of the maps here.", **champ("Equidistant")},
        {"key": "clarity", "name": "Clarity", "ours": 1, "type": "lent",
         "why": "Near-conformal; about 3 degrees mean angular distortion.", **champ("Conformal")},
        {"key": "balance", "name": "Balance", "ours": 1, "type": "lent",
         "why": "Near-equal-area, with less shape distortion than the classics.", **champ("EqualArea")},
    ]
    classic_keys = ["Winkel Tripel", "Robinson", "Mollweide", "Equal Earth", "Azimuthal Equidist.", "Equirectangular"]
    cl_by = {c["name"]: c for c in r["classics"]}
    classic_contenders = [{"key": "cl_" + nm.replace(" ", "_").replace(".", ""), "name": nm, "ours": 0,
                           "type": "classic", "classic": nm,
                           "why": {"Winkel Tripel": "Compromise projection; National Geographic's main world map since 1998.",
                                   "Robinson": "Compromise projection used by Rand McNally.",
                                   "Mollweide": "Equal-area ellipse.",
                                   "Equal Earth": "Equal-area projection, introduced in 2018.",
                                   "Azimuthal Equidist.": "Distances from the center point are correct.",
                                   "Equirectangular": "Longitude and latitude on a plain grid."}[nm],
                           "shape": round(cl_by[nm]["eps_shape"], 4), "area": round(cl_by[nm]["eps_area"], 4),
                           "dist": round(cl_by[nm]["eps_dist"], 4)} for nm in classic_keys]

    data = {
        "xterms": X_TERMS, "yterms": Y_TERMS, "nx": NX,
        "refs": {"shape": round(refs[0], 4), "area": round(refs[1], 4), "dist": round(refs[2], 4)},
        "winkel": {"shape": round(r["winkel_reference"]["eps_shape"], 4),
                   "area": round(r["winkel_reference"]["eps_area"], 4),
                   "dist": round(r["winkel_reference"]["eps_dist"], 4)},
        "families": fam_tris,
        "contenders": disc + classic_contenders,
        "land": simplify_land(DEFAULT_LAND),
    }
    os.makedirs(os.path.join(HERE, "site"), exist_ok=True)
    out = os.path.join(HERE, "site", "appdata.js")
    with open(out, "w") as f:
        f.write("window.APPDATA = ")
        json.dump(data, f, separators=(",", ":"))
        f.write(";\n")
    print(f"wrote {out} ({os.path.getsize(out)//1024} KB, {len(fam_tris)} families x {len(fam_tris[0]['points'])} pts, "
          f"{len(data['contenders'])} contenders, {len(data['land'])} rings)")


if __name__ == "__main__":
    main()
