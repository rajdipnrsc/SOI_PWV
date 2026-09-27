#!/usr/bin/env python3
"""Formal tuning experiment for the ZWD / gradient random-walk noise (TDS § 6.7).

    python validation/scripts/tune_zwd_noise.py FILELIST.txt REFDIR [--values 0.003 0.006 0.010 0.015]

FILELIST.txt: lines "SITE.o SITE.n [X Y Z]"; REFDIR: independent reference troposphere files (IGS/CODE
SINEX_TRO named <STATION>*<YYYYDDD>*.TRO[.gz]). For each sigma_ZWD the station-days are processed (settings
override file), hourly RMSE vs the reference and NIS/dof are collected, and the decision rule of § 6.7 is applied
(smallest sigma within 5 % of the minimum RMSE subject to 0.8 <= NIS/dof <= 1.2).
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import pwv_ppp  # noqa: E402
from ppp import settings  # noqa: E402
from ppp import timesys as ts  # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_ztd as cz  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("filelist")
    ap.add_argument("refdir")
    ap.add_argument("--values", nargs="*", type=float, default=[0.003, 0.006, 0.010, 0.015, 0.024])
    ap.add_argument("--out", default="validation/runs/tuning")
    a = ap.parse_args(argv)
    jobs = [L.split() for L in open(a.filelist) if L.strip() and not L.startswith("#")]
    table = []
    for sv in a.values:
        settings.SIGMA_ZWD_M_SQRT_H = sv
        d_all, nis = [], []
        for job in jobs:
            args = [job[0], job[1], "--out", os.path.join(a.out, f"s{sv * 1e3:.0f}"), "--no-pwv", "--quiet"]
            if len(job) >= 5:
                args += ["--xyz", *job[2:5]]
            try:
                res = pwv_ppp.run(pwv_ppp.parse_args(args))
            except pwv_ppp.Stop:
                continue
            nis.append(res["sol"].nis)
            fm = res["fm"]
            st = res["manifest"]["station"]
            y, doy = ts.year_doy(int(fm.t_gpst[len(fm.t_gpst) // 2]))
            refs = glob.glob(os.path.join(a.refdir, f"*{y:04d}{doy:03d}*TRO*"))
            if not refs:
                continue
            rt, rz, _ = cz.read_tro(refs[0], st)
            t, z = cz.hourly(fm.t_utc, (fm.zhd0 + fm.zwd) * 1e3)
            common, i1, i2 = np.intersect1d(t, rt, return_indices=True)
            d_all.extend(z[i1] - rz[i2])
        d = np.array(d_all)
        table.append({"sigma_zwd_mm_sqrt_h": sv * 1e3, "rmse_mm": float(np.sqrt(np.mean(d ** 2))) if len(d) else None,
                      "n": int(len(d)), "nis_per_dof": float(np.mean(nis)) if nis else None})
    ok = [r for r in table if r["rmse_mm"] is not None and r["nis_per_dof"] and 0.8 <= r["nis_per_dof"] <= 1.2]
    choice = None
    if ok:
        best = min(r["rmse_mm"] for r in ok)
        choice = min((r for r in ok if r["rmse_mm"] <= 1.05 * best), key=lambda r: r["sigma_zwd_mm_sqrt_h"])
    out = {"table": table, "decision": choice}
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "tuning.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
