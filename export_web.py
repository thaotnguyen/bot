#!/usr/bin/env python3
"""Bundle results + simplified Natural Earth coastlines into web/appdata.js
for the interactive explorer / poll artifact.

Reads results/results.json and a Natural Earth land GeoJSON (path via argv or
the scratchpad default) and writes web/appdata.js defining window.APPDATA.
"""

from __future__ import annotations

import json
import os
import sys

from mapopt.family import X_TERMS, Y_TERMS, NX
from mapopt.optimize import pareto_filter

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_LAND = "/tmp/claude-0/-home-user-bot/08015404-ca07-5463-b291-70bbc584ea17/scratchpad/ne_110m_land.geojson"


def simplify_land(path, prec=2, min_pts=4):
    d = json.load(open(path))
    rings = []
    for feat in d["features"]:
        geom = feat["geometry"]
        polys = geom["coordinates"] if geom["type"] == "Polygon" else [
            r for mp in geom["coordinates"] for r in mp]
        for ring in polys:
            pts = []
            last = None
            for lon, lat in ring:
                p = [round(lon, prec), round(lat, prec)]
                if p != last:
                    pts.append(p)
                    last = p
            if len(pts) >= min_pts:
                rings.append(pts)
    return rings


def split_params(params):
    return params[:NX], params[NX:]


def champ(entry):
    a, b = split_params(entry["params"])
    return {
        "a": [round(x, 6) for x in a],
        "b": [round(x, 6) for x in b],
        "shape": round(entry["eps_shape"], 4),
        "area": round(entry["eps_area"], 4),
        "combined": round(entry["eps_shape"] + entry["eps_area"], 4),
        "ang": round(entry["rms_angular_deg"], 2),
        "areapct": round(entry["rms_area_pct"], 2),
        "desc": entry.get("desc", ""),
    }


def main():
    land_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_LAND
    r = json.load(open(os.path.join(HERE, "results", "results.json")))

    frontier = []
    for f in sorted(pareto_filter(r["frontier_area_weighted"]), key=lambda x: x["t"]):
        a, b = split_params(f["params"])
        frontier.append({
            "t": f["t"],
            "a": [round(x, 6) for x in a],
            "b": [round(x, 6) for x in b],
            "shape": round(f["eps_shape"], 4),
            "area": round(f["eps_area"], 4),
            "ang": round(f["rms_angular_deg"], 2),
        })

    champions = {}
    for k in ["Equipoise", "WinkelDominator", "ShapeChampion", "AreaChampion", "Anthropocene"]:
        e = r["champions"].get(k)
        if e:
            champions[k] = champ(e)

    classics = [{
        "name": c["name"],
        "shape": round(c["eps_shape"], 4),
        "area": round(c["eps_area"], 4),
        "combined": round(c["combined"], 4),
        "ang": round(c["rms_angular_deg"], 2),
        "equal_area": c["equal_area"],
        "conformal": c["conformal"],
    } for c in r["classics_area_weighted"]]

    data = {
        "xterms": X_TERMS,
        "yterms": Y_TERMS,
        "frontier": frontier,
        "champions": champions,
        "classics": classics,
        "winkel": {
            "shape": round(r["winkel_reference"]["eps_shape"], 4),
            "area": round(r["winkel_reference"]["eps_area"], 4),
            "combined": round(r["winkel_reference"]["combined"], 4),
        },
        "novelty": r["novelty"],
        "population_centers": [
            {"name": c[0], "lat": c[1], "lon": c[2], "w": c[3]} for c in r["population_centers"]
        ],
        "land": simplify_land(land_path),
        "meta": r["meta"],
    }

    os.makedirs(os.path.join(HERE, "web"), exist_ok=True)
    out = os.path.join(HERE, "web", "appdata.js")
    with open(out, "w") as f:
        f.write("window.APPDATA = ")
        json.dump(data, f, separators=(",", ":"))
        f.write(";\n")
    print(f"wrote {out}  ({os.path.getsize(out)//1024} KB, {len(data['land'])} land rings, "
          f"{len(frontier)} frontier pts, {len(champions)} champions)")


if __name__ == "__main__":
    main()
