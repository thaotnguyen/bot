#!/usr/bin/env python3
"""Bundle 3-objective results + simplified Natural Earth coastlines into
web/appdata.js for the interactive ternary explorer / poll artifact."""

from __future__ import annotations

import json
import os
import sys

from mapopt.family import X_TERMS, Y_TERMS, NX

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
            pts, last = [], None
            for lon, lat in ring:
                p = [round(lon, prec), round(lat, prec)]
                if p != last:
                    pts.append(p)
                    last = p
            if len(pts) >= min_pts:
                rings.append(pts)
    return rings


def ab(params):
    return [round(x, 6) for x in params[:NX]], [round(x, 6) for x in params[NX:]]


def emit(entry, refs, extra=None):
    a, b = ab(entry["params"])
    d = {
        "a": a, "b": b,
        "shape": round(entry["eps_shape"], 4),
        "area": round(entry["eps_area"], 4),
        "dist": round(entry["eps_dist"], 4),
        "ang": round(entry["rms_angular_deg"], 2),
        "norm": round(entry["eps_shape"]/refs[0] + entry["eps_area"]/refs[1] + entry["eps_dist"]/refs[2], 3),
    }
    if "desc" in entry:
        d["desc"] = entry["desc"]
    if extra:
        d.update(extra)
    return d


def main():
    land_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_LAND
    r = json.load(open(os.path.join(HERE, "results", "results.json")))
    R = r["meta"]["refs_equirect"]
    refs = (R["shape"], R["area"], R["dist"])

    surface = [emit(p, refs, {"w": [round(x, 4) for x in p["w"]]}) for p in r["surface"]]

    champions = {}
    for k, e in r["champions"].items():
        if e:
            champions[k] = emit(e, refs)

    classics = [{
        "name": c["name"],
        "shape": round(c["eps_shape"], 4), "area": round(c["eps_area"], 4),
        "dist": round(c["eps_dist"], 4), "norm": round(c["combined_norm"], 3),
        "ang": round(c["rms_angular_deg"], 2),
        "equal_area": c["equal_area"], "conformal": c["conformal"],
    } for c in r["classics"]]

    w = r["winkel_reference"]
    data = {
        "xterms": X_TERMS, "yterms": Y_TERMS,
        "refs": {"shape": round(refs[0], 4), "area": round(refs[1], 4), "dist": round(refs[2], 4)},
        "surface": surface,
        "champions": champions,
        "classics": classics,
        "winkel": {"shape": round(w["eps_shape"], 4), "area": round(w["eps_area"], 4),
                   "dist": round(w["eps_dist"], 4), "norm": round(w["combined_norm"], 3)},
        "novelty": r["novelty"],
        "land": simplify_land(land_path),
        "meta": r["meta"],
    }
    os.makedirs(os.path.join(HERE, "web"), exist_ok=True)
    out = os.path.join(HERE, "web", "appdata.js")
    with open(out, "w") as f:
        f.write("window.APPDATA = ")
        json.dump(data, f, separators=(",", ":"))
        f.write(";\n")
    print(f"wrote {out} ({os.path.getsize(out)//1024} KB, {len(data['land'])} rings, "
          f"{len(surface)} surface pts, {len(champions)} champions, {len(classics)} classics)")


if __name__ == "__main__":
    main()
