#!/usr/bin/env python3
"""Autoresearch pipeline (3-objective): hill-climb the shape / area / distance
Pareto surface of world map projections, benchmark against the classics, and
emit results.

Run:  python3 run_research.py   ->  results/results.json  + summary table.
"""

from __future__ import annotations

import json
import math
import os
import time
from typing import Dict, List

from mapopt import classics, metrics, optimize
from mapopt.family import forward, X_TERMS, Y_TERMS, NPARAMS, Grid

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")


def report_scores(fn_or_params, grid, dsamp, is_family):
    if is_family:
        s = metrics.family_all(fn_or_params, grid, dsamp, grid.area_w)
    else:
        s = metrics.generic_all(fn_or_params, grid, dsamp, grid.area_w)
    return s


def pack(params, s, refs, extra=None):
    d = {
        "params": [round(p, 8) for p in params] if params is not None else None,
        "eps_shape": s["eps_shape"],
        "eps_area": s["eps_area"],
        "eps_dist": s["eps_dist"],
        "rms_angular_deg": s["rms_angular_deg"],
        "combined_norm": s["eps_shape"] / refs[0] + s["eps_area"] / refs[1] + s["eps_dist"] / refs[2],
    }
    if extra:
        d.update(extra)
    return d


def novelty(params, grid):
    champ = [forward(params, grid.lats[i], grid.lons[i]) for i in range(grid.n)]
    cx = [p[0] for p in champ]
    cy = [p[1] for p in champ]
    denom = sum(x * x for x in cx) + sum(y * y for y in cy)
    best_name, best = None, 1e9
    for name, fn in classics.CLASSICS.items():
        cs = [fn(grid.lats[i], grid.lons[i]) for i in range(grid.n)]
        ax = [p[0] for p in cs]
        ay = [p[1] for p in cs]
        sx = sum(cx[i] * ax[i] for i in range(grid.n)) / (sum(v * v for v in ax) or 1)
        sy = sum(cy[i] * ay[i] for i in range(grid.n)) / (sum(v * v for v in ay) or 1)
        res = sum((cx[i] - sx * ax[i]) ** 2 for i in range(grid.n)) + \
              sum((cy[i] - sy * ay[i]) ** 2 for i in range(grid.n))
        pct = 100.0 * math.sqrt(res / denom)
        if pct < best:
            best, best_name = pct, name
    return {"nearest_classic": best_name, "residual_pct": round(best, 2)}


def find_w(sweep, w):
    return min(sweep, key=lambda p: sum((p["w"][i] - w[i]) ** 2 for i in range(3)))


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)

    OPT = Grid(lat_step=4.0, lon_step=6.0)
    REP = Grid(lat_step=2.5, lon_step=4.0)
    OPT_D = metrics.DistanceSampler(n_anchor=90, seed=7)
    REP_D = metrics.DistanceSampler(n_anchor=150, seed=101)
    print(f"opt grid {OPT.n} pts / {OPT_D.npairs} pairs | report grid {REP.n} pts / {REP_D.npairs} pairs")

    # reference normalizers = equirectangular's three errors (report grid/sampler)
    eqr = metrics.generic_all(classics.equirectangular, REP, REP_D)
    refs = (eqr["eps_shape"], eqr["eps_area"], eqr["eps_dist"])
    print(f"refs (equirect): shape={refs[0]:.4f} area={refs[1]:.4f} dist={refs[2]:.4f}")

    # ---- classics -------------------------------------------------------
    print("scoring classics (shape/area/distance)...")
    classic_rows = []
    for name, fn in classics.CLASSICS.items():
        s = metrics.generic_all(fn, REP, REP_D)
        classic_rows.append(pack(None, s, refs, {
            "name": name,
            "equal_area": name in classics.EQUAL_AREA,
            "conformal": name in classics.CONFORMAL,
        }))
    classic_rows.sort(key=lambda r: r["combined_norm"])
    winkel = next(r for r in classic_rows if r["name"] == "Winkel Tripel")

    # ---- 3-objective simplex sweep --------------------------------------
    print("hill-climbing the shape/area/distance Pareto surface (15 weightings)...")
    weights = optimize.simplex_weights(5)
    raw = optimize.simplex_sweep(OPT, OPT_D, refs, weights, steps=65, lr=0.02)
    surface = []
    for pt in raw:
        s = metrics.family_all(pt["params"], REP, REP_D, REP.area_w)  # re-score on report grid
        surface.append(pack(pt["params"], s, refs, {"w": pt["w"], "n_fold": int(s["n_fold"])}))

    # champions (corners + centroid + minimax all-rounder + winkel dominator)
    def at(w):
        return find_w(surface, w)
    champions = {
        "Conformal":   {**at((1, 0, 0)), "desc": "shape-optimal corner (angles preserved)"},
        "EqualArea":   {**at((0, 1, 0)), "desc": "area-optimal corner (sizes honest)"},
        "Equidistant": {**at((0, 0, 1)), "desc": "distance-optimal corner (true ruler)"},
        "Equipoise":   {**at((1/3, 1/3, 1/3)), "desc": "equal blend of all three objectives"},
    }
    # minimax all-rounder: the surface point whose WORST normalized axis is lowest
    def worst(p):
        return max(p["eps_shape"] / refs[0], p["eps_area"] / refs[1], p["eps_dist"] / refs[2])
    triathlon = min(surface, key=worst)
    champions["Triathlon"] = {**triathlon, "desc": "minimizes its single worst axis (most even map)"}
    # dominate Winkel on all three? first check the swept surface...
    def dominates(p):
        return (p["eps_shape"] <= winkel["eps_shape"] and p["eps_area"] <= winkel["eps_area"]
                and p["eps_dist"] <= winkel["eps_dist"])
    doms = [p for p in surface if dominates(p)]
    dominator = min(doms, key=lambda p: p["combined_norm"]) if doms else None

    # ...otherwise run a targeted solve that minimizes the WORST ratio to Winkel
    if dominator is None:
        wk_opt = metrics.generic_all(classics.winkel_tripel, OPT, OPT_D)
        tgt = (wk_opt["eps_shape"], wk_opt["eps_area"], wk_opt["eps_dist"])
        beta = 10.0

        def f_beat(p):
            from mapopt import family
            sp = metrics.family_scores(p, OPT, OPT.area_w)
            a, b = family.split(p)
            dist = OPT_D.score(lambda la, lo: family.forward_ab(a, b, la, lo))
            rs, ra, rd = sp["eps_shape"]/tgt[0], sp["eps_area"]/tgt[1], dist/tgt[2]
            mx = max(rs, ra, rd)
            lse = mx + math.log(math.exp(beta*(rs-mx)) + math.exp(beta*(ra-mx)) + math.exp(beta*(rd-mx)))/beta
            return lse + 50.0 * sp["fold_pen"]

        x_beat = optimize.optimize_obj(f_beat, triathlon["params"], OPT, steps=140, lr=0.015)
        s_beat = metrics.family_all(x_beat, REP, REP_D, REP.area_w)
        cand = pack(x_beat, s_beat, refs, {"w": None, "n_fold": int(s_beat["n_fold"])})
        if dominates(cand):
            dominator = cand
    champions["WinkelDominator"] = ({**dominator, "desc": "beats Winkel Tripel on shape, area AND distance at once"}
                                    if dominator else None)

    nov = {k: novelty(champions[k]["params"], REP) for k in ["Equipoise", "Triathlon", "Equidistant"]}

    results = {
        "meta": {
            "generated": "2026-09-18",
            "opt_grid_points": OPT.n, "opt_pairs": OPT_D.npairs,
            "report_grid_points": REP.n, "report_pairs": REP_D.npairs,
            "nparams": NPARAMS, "x_terms": X_TERMS, "y_terms": Y_TERMS,
            "objectives": {
                "shape": "eps_shape = <ln(a/b)^2> (local angular, area-weighted)",
                "area": "eps_area = Var(ln a*b) (local areal)",
                "dist": "eps_dist = Var(ln(d_map/d_globe)) over point pairs (global)",
            },
            "refs_equirect": {"shape": refs[0], "area": refs[1], "dist": refs[2]},
            "runtime_s": None,
        },
        "classics": classic_rows,
        "surface": surface,
        "champions": champions,
        "winkel_reference": winkel,
        "novelty": nov,
    }
    results["meta"]["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(OUT, "results.json"), "w") as f:
        json.dump(results, f, indent=2)

    # ---- summary --------------------------------------------------------
    def line(name, r):
        print(f'{name:20s} shape={r["eps_shape"]:.3f} area={r["eps_area"]:.3f} '
              f'dist={r["eps_dist"]:.4f} | norm={r["combined_norm"]:.3f} ang={r["rms_angular_deg"]:.1f}')
    print("\n=== CLASSICS (sorted by combined normalized error) ===")
    for r in classic_rows:
        line(r["name"], r)
    print("\n=== DISCOVERED CHAMPIONS ===")
    for k in ["Conformal", "EqualArea", "Equidistant", "Equipoise", "Triathlon"]:
        line(k, champions[k])
    print("\n=== HEADLINE vs Winkel Tripel ===")
    line("Winkel Tripel", winkel)
    line("Equipoise", champions["Equipoise"])
    line("Triathlon", champions["Triathlon"])
    if dominator:
        line("WinkelDominator", champions["WinkelDominator"])
        print("  ^ lower on shape AND area AND distance simultaneously")
    else:
        print("  (no single family map dominates Winkel on all three axes)")
    print("\nNovelty (residual vs nearest classic after best rescale):")
    for k, v in nov.items():
        print(f'  {k:12s} nearest={v["nearest_classic"]:20s} {v["residual_pct"]}%')
    print(f'\nDone in {results["meta"]["runtime_s"]}s -> results/results.json')


if __name__ == "__main__":
    main()
