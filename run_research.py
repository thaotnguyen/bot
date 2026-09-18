#!/usr/bin/env python3
"""Autoresearch pipeline: hill-climb the Pareto frontier of world map
projections, benchmark against the classics, and emit results.

Run:  python3 run_research.py
Writes results/results.json and prints a summary table.
"""

from __future__ import annotations

import json
import math
import os
import time
from typing import Dict, List

from mapopt import classics, metrics, optimize, population
from mapopt.family import Grid, forward, split, X_TERMS, Y_TERMS, NPARAMS

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")


def score_classics(grid, weights) -> List[Dict]:
    rows = []
    for name, fn in classics.CLASSICS.items():
        s = metrics.generic_scores(fn, grid, weights)
        rows.append({
            "name": name,
            "eps_shape": s["eps_shape"],
            "eps_area": s["eps_area"],
            "combined": s["eps_shape"] + s["eps_area"],
            "rms_angular_deg": s["rms_angular_deg"],
            "rms_area_pct": s["rms_area_pct"],
            "equal_area": name in classics.EQUAL_AREA,
            "conformal": name in classics.CONFORMAL,
        })
    return rows


def report_point(params, report_grid, weights) -> Dict:
    s = metrics.family_scores(params, report_grid, weights)
    return {
        "params": [round(p, 8) for p in params],
        "eps_shape": s["eps_shape"],
        "eps_area": s["eps_area"],
        "combined": s["eps_shape"] + s["eps_area"],
        "rms_angular_deg": s["rms_angular_deg"],
        "rms_area_pct": s["rms_area_pct"],
        "n_fold": int(s["n_fold"]),
    }


def novelty(params, grid) -> Dict:
    """How close is this projection to an anisotropic rescaling of a classic?

    Returns the nearest classic and the normalized residual (%). A large
    residual means the discovered map is genuinely not a stretched classic.
    """
    cx_champ = [forward(params, grid.lats[i], grid.lons[i]) for i in range(grid.n)]
    champ_x = [p[0] for p in cx_champ]
    champ_y = [p[1] for p in cx_champ]
    denom = sum(x * x for x in champ_x) + sum(y * y for y in champ_y)
    best_name, best_res = None, 1e9
    per = {}
    for name, fn in classics.CLASSICS.items():
        cs = [fn(grid.lats[i], grid.lons[i]) for i in range(grid.n)]
        cxx = [p[0] for p in cs]
        cyy = [p[1] for p in cs]
        sxn = sum(champ_x[i] * cxx[i] for i in range(grid.n))
        sxd = sum(cxx[i] * cxx[i] for i in range(grid.n)) or 1.0
        syn = sum(champ_y[i] * cyy[i] for i in range(grid.n))
        syd = sum(cyy[i] * cyy[i] for i in range(grid.n)) or 1.0
        sx, sy = sxn / sxd, syn / syd
        res = sum((champ_x[i] - sx * cxx[i]) ** 2 for i in range(grid.n))
        res += sum((champ_y[i] - sy * cyy[i]) ** 2 for i in range(grid.n))
        pct = 100.0 * math.sqrt(res / denom)
        per[name] = round(pct, 2)
        if pct < best_res:
            best_res, best_name = pct, name
    return {"nearest_classic": best_name, "residual_pct": round(best_res, 2), "per_classic": per}


def main():
    t_start = time.time()
    os.makedirs(OUT, exist_ok=True)

    OPT = Grid(lat_step=4.0, lon_step=6.0)      # optimization grid (fast)
    REP = Grid(lat_step=2.5, lon_step=4.0)      # reporting grid (fine, official numbers)
    print(f"opt grid: {OPT.n} pts | report grid: {REP.n} pts")

    pop_opt = population.weights_for_grid(OPT)
    pop_rep = population.weights_for_grid(REP)

    # ---- baselines -------------------------------------------------------
    print("scoring classic projections...")
    classic_rows = score_classics(REP, REP.area_w)
    classic_rows.sort(key=lambda r: r["combined"])
    winkel = next(r for r in classic_rows if r["name"] == "Winkel Tripel")

    # ---- area-weighted Pareto sweep -------------------------------------
    print("hill-climbing area-weighted Pareto frontier...")
    t_list = [0.05, 0.12, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.88, 0.94]
    sweep = optimize.pareto_sweep(OPT, OPT.area_w, t_list, steps=70, lr=0.02)
    frontier = []
    for pt in sweep:
        rp = report_point(pt["params"], REP, REP.area_w)
        rp["t"] = pt["t"]
        frontier.append(rp)
    frontier.sort(key=lambda r: r["eps_shape"])

    # champions from the area-weighted frontier
    equipoise = min(frontier, key=lambda r: r["combined"])
    shape_champ = min(frontier, key=lambda r: r["eps_shape"])
    area_champ = min(frontier, key=lambda r: r["eps_area"])
    dominators = [r for r in frontier
                  if r["eps_shape"] <= winkel["eps_shape"] and r["eps_area"] <= winkel["eps_area"]]
    dominator = min(dominators, key=lambda r: r["combined"]) if dominators else None

    # ---- population-weighted sweep --------------------------------------
    print("hill-climbing population-weighted frontier...")
    pop_sweep = optimize.pareto_sweep(OPT, pop_opt, [0.2, 0.35, 0.5, 0.65, 0.8], steps=70, lr=0.02)
    pop_frontier = []
    for pt in pop_sweep:
        rp = report_point(pt["params"], REP, pop_rep)      # scored under population weight
        rp_area = metrics.family_scores(pt["params"], REP, REP.area_w)  # also under area weight
        rp["t"] = pt["t"]
        rp["eps_shape_area_wt"] = rp_area["eps_shape"]
        rp["eps_area_area_wt"] = rp_area["eps_area"]
        pop_frontier.append(rp)
    anthropocene = min(pop_frontier, key=lambda r: r["combined"])

    # how do the classics do under the population weighting?
    classic_pop = score_classics(REP, pop_rep)
    classic_pop.sort(key=lambda r: r["combined"])

    # ---- novelty ---------------------------------------------------------
    print("checking novelty vs classics...")
    nov = {
        "Equipoise": novelty(equipoise["params"], REP),
        "Anthropocene": novelty(anthropocene["params"], REP),
        "ShapeChampion": novelty(shape_champ["params"], REP),
    }

    results = {
        "meta": {
            "generated": "2026-09-18",
            "opt_grid_points": OPT.n,
            "report_grid_points": REP.n,
            "nparams": NPARAMS,
            "x_terms": X_TERMS,
            "y_terms": Y_TERMS,
            "family": "x=sum a_ij u^(2i+1) v^(2j); y=sum b_ij u^(2i) v^(2j+1); u=lon/pi, v=lat/(pi/2)",
            "metric": "eps_shape=<ln(a/b)^2>, eps_area=Var(ln(a*b)), area-weighted; combined=sum",
            "runtime_s": None,
        },
        "classics_area_weighted": classic_rows,
        "classics_population_weighted": classic_pop,
        "frontier_area_weighted": frontier,
        "frontier_population_weighted": pop_frontier,
        "champions": {
            "Equipoise": {**equipoise, "desc": "min combined shape+area error (area-weighted)"},
            "ShapeChampion": {**shape_champ, "desc": "min conformality error in the family"},
            "AreaChampion": {**area_champ, "desc": "near-equal-area member of the family"},
            "WinkelDominator": ({**dominator, "desc": "beats Winkel Tripel on BOTH axes"}
                                if dominator else None),
            "Anthropocene": {**anthropocene, "desc": "min combined error weighted by population"},
        },
        "winkel_reference": winkel,
        "novelty": nov,
        "population_centers": population.CENTERS,
    }
    results["meta"]["runtime_s"] = round(time.time() - t_start, 1)

    with open(os.path.join(OUT, "results.json"), "w") as f:
        json.dump(results, f, indent=2)

    # ---- print summary ---------------------------------------------------
    print("\n=== CLASSICS (area-weighted, sorted by combined) ===")
    print(f'{"projection":16s} {"shape":>8s} {"area":>8s} {"combined":>9s} {"ang°":>6s}')
    for r in classic_rows:
        print(f'{r["name"]:16s} {r["eps_shape"]:8.4f} {r["eps_area"]:8.4f} {r["combined"]:9.4f} {r["rms_angular_deg"]:6.1f}')

    print("\n=== HILL-CLIMBED FRONTIER (area-weighted) ===")
    print(f'{"t":>5s} {"shape":>8s} {"area":>8s} {"combined":>9s} {"ang°":>6s} {"folds":>5s}')
    for r in frontier:
        print(f'{r["t"]:5.2f} {r["eps_shape"]:8.4f} {r["eps_area"]:8.4f} {r["combined"]:9.4f} {r["rms_angular_deg"]:6.1f} {r["n_fold"]:5d}')

    print("\n=== HEADLINE ===")
    print(f'Winkel Tripel     : shape={winkel["eps_shape"]:.4f} area={winkel["eps_area"]:.4f} combined={winkel["combined"]:.4f}')
    print(f'Equipoise (ours)  : shape={equipoise["eps_shape"]:.4f} area={equipoise["eps_area"]:.4f} combined={equipoise["combined"]:.4f}')
    if dominator:
        print(f'WinkelDominator   : shape={dominator["eps_shape"]:.4f} area={dominator["eps_area"]:.4f} combined={dominator["combined"]:.4f}  <-- beats Winkel on BOTH axes')
    print(f'Anthropocene(pop) : shape={anthropocene["eps_shape"]:.4f} area={anthropocene["eps_area"]:.4f} combined={anthropocene["combined"]:.4f} (population-weighted)')
    print("\nNovelty (min normalized residual vs any classic, higher=more novel):")
    for k, v in nov.items():
        print(f'  {k:14s} nearest={v["nearest_classic"]:14s} residual={v["residual_pct"]}%')
    print(f'\nDone in {results["meta"]["runtime_s"]}s -> results/results.json')


if __name__ == "__main__":
    main()
