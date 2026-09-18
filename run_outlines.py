#!/usr/bin/env python3
"""Outline meta-search: treat the projection's OUTLINE / topology as the search
variable. Optimize every outline family (rectangle, oval, lens, disc, and a
free-exponent superellipse) on shape/area/distance and see which outline wins.

Run:  python3 run_outlines.py   ->  results/outlines.json + summary.
"""

from __future__ import annotations

import json, math, os, time

from mapopt import classics, metrics, optimize, families
from mapopt.family import Grid

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")

BLENDS = {"shape": (1, 0, 0), "area": (0, 1, 0), "distance": (0, 0, 1),
          "balanced": (1/3, 1/3, 1/3)}


def rec(s, refs, extra=None):
    d = {"shape": round(s["eps_shape"], 4), "area": round(s["eps_area"], 4),
         "dist": round(s["eps_dist"], 4), "ang": round(s["rms_angular_deg"], 2),
         "folds": int(s["n_fold"]),
         "norm": round(s["eps_shape"]/refs[0] + s["eps_area"]/refs[1] + s["eps_dist"]/refs[2], 3)}
    if extra:
        d.update(extra)
    return d


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    OPT = Grid(lat_step=4.0, lon_step=6.0)
    REP = Grid(lat_step=3.0, lon_step=4.0)
    OPT_D = metrics.DistanceSampler(n_anchor=90, seed=7)
    REP_D = metrics.DistanceSampler(n_anchor=140, seed=101)
    eqr = metrics.generic_all(classics.equirectangular, REP, REP_D)
    refs = (eqr["eps_shape"], eqr["eps_area"], eqr["eps_dist"])
    print(f"opt {OPT.n}pts/{OPT_D.npairs}pairs  report {REP.n}pts/{REP_D.npairs}pairs")

    fam_results = []
    for fam in families.ALL_FAMILIES:
        print(f"optimizing {fam.label} ({fam.outline}) ...")
        blends = {}
        x_bal = None
        for name in ["balanced", "shape", "area", "distance"]:   # balanced first, warm-start corners
            w = BLENDS[name]
            x0 = fam.p0() if name == "balanced" else x_bal
            x, _ = optimize.optimize_family(fam, OPT, OPT_D, refs, w, x0=x0, steps=90, lr=0.02)
            s = metrics.famscore(fam, x, REP, REP_D)          # re-score on report grid
            blends[name] = rec(s, refs, {"note": fam.describe(x), "params": [round(v, 6) for v in x]})
            if name == "balanced":
                x_bal = x
        fam_results.append({"label": fam.label, "outline": fam.outline, "blends": blends})

    # superellipse: fix p, optimize the rest at the balanced blend -> error vs "roundness"
    print("superellipse p-sweep (error vs outline roundness) ...")
    p_sweep = []
    se = families.Superellipse()
    for pval in [1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0]:
        fam = families.Superellipse()
        fam.bounds = [(pval, pval), (0.05, 20), (0.05, 20), None, None, None]
        x0 = [pval, 2.0, 1.0, math.pi/2, 0.0, 0.0]
        x, _ = optimize.optimize_family(fam, OPT, OPT_D, refs, (1/3, 1/3, 1/3), x0=x0, steps=70, lr=0.02)
        s = metrics.famscore(fam, x, REP, REP_D)
        p_sweep.append({"p": pval, **rec(s, refs, {"params": [round(v, 6) for v in x]})})

    results = {
        "meta": {"generated": "2026-09-18", "report_grid_points": REP.n,
                 "report_pairs": REP_D.npairs,
                 "refs_equirect": {"shape": refs[0], "area": refs[1], "dist": refs[2]},
                 "runtime_s": None},
        "families": fam_results,
        "superellipse_p_sweep": p_sweep,
    }
    results["meta"]["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, "outlines.json"), "w") as f:
        json.dump(results, f, indent=2)

    # ---- summary ----
    print("\n=== BEST SCORE PER OUTLINE FAMILY, PER OBJECTIVE (lower better) ===")
    print(f'{"outline":20s} {"shape":>8s} {"area":>8s} {"distance":>9s} {"balanced-norm":>14s}')
    for fr in fam_results:
        b = fr["blends"]
        print(f'{fr["outline"]:20s} {b["shape"]["shape"]:8.3f} {b["area"]["area"]:8.3f} '
              f'{b["distance"]["dist"]:9.4f} {b["balanced"]["norm"]:14.3f}')
    print("\nWinners:")
    for obj, key in [("shape", "shape"), ("area", "area"), ("distance", "dist"), ("balanced (norm)", "norm")]:
        blend = {"shape": "shape", "area": "area", "distance": "distance", "balanced (norm)": "balanced"}[obj]
        metric = key
        best = min(fam_results, key=lambda fr: fr["blends"][blend][metric])
        val = best["blends"][blend][metric]
        print(f'  best {obj:16s}: {best["outline"]:22s} ({val})')
    print("\n=== SUPERELLIPSE: error vs outline roundness (balanced) ===")
    print(f'{"p (exponent)":14s} {"outline":16s} {"balanced-norm":>14s}')
    shapes = {1.2: "diamond-ish", 1.5: "soft diamond", 2.0: "ellipse", 2.5: "rounded", 3.0: "rounded-rect", 4.0: "rect-ish", 6.0: "rectangle"}
    for ps in p_sweep:
        print(f'{ps["p"]:<14.1f} {shapes[ps["p"]]:16s} {ps["norm"]:14.3f}')
    bestp = min(p_sweep, key=lambda x: x["norm"])
    print(f'  -> optimal outline exponent p = {bestp["p"]} ({shapes[bestp["p"]]}), norm {bestp["norm"]}')
    print(f'\nDone in {results["meta"]["runtime_s"]}s -> results/outlines.json')


if __name__ == "__main__":
    main()
