#!/usr/bin/env python3
"""Coordinate sensitivity experiment (TDS § 10.5): fixed-mode ZTD with perturbed coordinates -> beta.

    python validation/scripts/coordinate_sensitivity.py SITE.o SITE.n --xyz X Y Z [--du 0.01 0.03 0.05 0.10]

The SOI coordinate (ITRF2008 @ 2005.0) is perturbed along the local up / north / east direction before the
standard transformation, the station-day is processed in fixed mode for each perturbation, and the regression
slope of mean dZTD vs dU gives the processor-specific height-to-ZTD coupling beta (expected 0.1-0.6).
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import pwv_ppp  # noqa: E402
from ppp import coords as co  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("obs")
    ap.add_argument("nav")
    ap.add_argument("--xyz", nargs=3, type=float, required=True)
    ap.add_argument("--du", nargs="*", type=float, default=[-0.10, -0.05, -0.03, -0.01, 0.01, 0.03, 0.05, 0.10])
    ap.add_argument("--dh", nargs="*", type=float, default=[0.03, 0.10])
    ap.add_argument("--out", default="validation/runs/sensitivity")
    ap.add_argument("--offline", action="store_true")
    a = ap.parse_args(argv)
    X = np.array(a.xyz)
    lat, lon, _ = co.ecef_to_geodetic(X)
    runs = {}

    def run(tag, d_enu):
        xyz = X + co.enu2xyz(np.array(d_enu), lat, lon)
        args = ["--xyz", *map(str, xyz), "--mode", "fixed", "--out", os.path.join(a.out, tag), "--no-pwv", "--quiet"]
        if a.offline:
            args.append("--offline")
        res = pwv_ppp.run(pwv_ppp.parse_args([a.obs, a.nav] + args))
        fm = res["fm"]
        return (fm.zhd0 + fm.zwd) * 1e3
    ref = run("ref", [0, 0, 0])
    rows = []
    for du in a.du:
        z = run(f"dU{du:+.2f}", [0, 0, du])
        d = z - ref
        rows.append({"dU_m": du, "mean_dZTD_mm": float(np.nanmean(d)), "std_dZTD_mm": float(np.nanstd(d))})
    for dh in a.dh:
        for name, enu in (("dN", [0, dh, 0]), ("dE", [dh, 0, 0])):
            z = run(f"{name}{dh:+.2f}", enu)
            d = z - ref
            rows.append({name + "_m": dh, "mean_dZTD_mm": float(np.nanmean(d)), "std_dZTD_mm": float(np.nanstd(d))})
    du = np.array([r["dU_m"] for r in rows if "dU_m" in r])
    mz = np.array([r["mean_dZTD_mm"] for r in rows if "dU_m" in r])
    beta = float(np.polyfit(du * 1e3, mz, 1)[0])
    out = {"beta": beta, "beta_in_expected_range": 0.1 <= abs(beta) <= 0.6, "rows": rows}
    print(json.dumps(out, indent=1))
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "sensitivity.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
