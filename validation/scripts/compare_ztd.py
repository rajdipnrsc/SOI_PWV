#!/usr/bin/env python3
"""Compare station ZTD with a reference (PRIDE, IGS/CODE troposphere SINEX, ERA5 ...) - TDS § 30.2-30.3.

    python validation/scripts/compare_ztd.py OURS_ZTD.csv REFERENCE [--ref-format tro|csv] [--hourly] [--plot out.png]

REFERENCE: a SINEX_TRO file (IGS final troposphere / CODE COD0OPSFIN_..._TRO.TRO, TROTOT in mm) or a CSV with
columns time (ISO, UTC) and ztd_mm [and sigma_mm]. Reported metrics: bias, STD, RMSE, MAE, Pearson r, P95 |d|,
normalised-error calibration, 24-h/12-h harmonic amplitude of the difference (diurnal composite), and the gate
verdicts for Stage 3 and Stage 5 (harmonised PRIDE benchmark) or Level 2 (independent reference).
PRIDE and ERA5 are references, not truth (TDS § 36.13).
"""
import argparse
import csv
import gzip
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from ppp import timesys as ts  # noqa: E402


def read_ours(path):
    rows = list(csv.DictReader(L for L in open(path) if not L.startswith("#")))
    t, z, s, st = [], [], [], None
    for r in rows:
        if not r["ZTD"]:
            continue
        st = r["station"]
        t.append(_iso(r["timestamp_utc"]))
        z.append(float(r["ZTD"]))
        s.append(float(r["sigma_ZTD"]))
    return st, np.array(t, dtype=np.int64), np.array(z), np.array(s)


def _iso(s):
    s = s.rstrip("Z").replace(" ", "T")
    d, tt = s.split("T")
    y, m, dd = (int(x) for x in d.split("-"))
    hh, mi, sec = tt.split(":")
    return ts.to_ns(y, m, dd, int(hh), int(mi), float(sec))


def read_tro(path, station):
    """SINEX_TRO (1.00/2.00): +TROP/SOLUTION records 'SITE YY:DDD:SSSSS TROTOT STDDEV ...' (GPST epochs)."""
    op = gzip.open if path.endswith(".gz") else open
    t, z, s = [], [], []
    inblk = False
    with op(path, "rt") as fh:
        for L in fh:
            if L.startswith("+TROP/SOLUTION"):
                inblk = True
                continue
            if L.startswith("-TROP/SOLUTION"):
                break
            if not inblk or L.startswith("*") or not L.strip():
                continue
            p = L.split()
            if p[0][:4].upper() != station[:4].upper():
                continue
            y, doy, sod = p[1].split(":")
            y = int(y)
            y = y + (2000 if y < 80 else 1900) if y < 100 else y
            tg = ts.to_ns(y, 1, 1) + (int(doy) - 1) * ts.DAY_NS + int(sod) * ts.NS
            t.append(ts.gpst_to_utc(tg))
            z.append(float(p[2]))
            s.append(float(p[3]) if len(p) > 3 else np.nan)
    return np.array(t, dtype=np.int64), np.array(z), np.array(s)


def read_ref_csv(path):
    rows = list(csv.DictReader(L for L in open(path) if not L.startswith("#")))
    t = np.array([_iso(r["time"]) for r in rows], dtype=np.int64)
    z = np.array([float(r["ztd_mm"]) for r in rows])
    s = np.array([float(r.get("sigma_mm") or "nan") for r in rows])
    return t, z, s


def hourly(t, z):
    h = t // (3600 * ts.NS)
    out_t, out_z = [], []
    for k in np.unique(h):
        m = h == k
        if m.sum() >= 8:
            out_t.append(k * 3600 * ts.NS + 1800 * ts.NS)
            out_z.append(np.mean(z[m]))
    return np.array(out_t, dtype=np.int64), np.array(out_z)


def metrics(d, s_comb=None, t=None):
    out = {"n": int(len(d)), "bias_mm": float(np.mean(d)), "std_mm": float(np.std(d)),
           "rmse_mm": float(np.sqrt(np.mean(d ** 2))), "mae_mm": float(np.mean(np.abs(d))),
           "p95_mm": float(np.percentile(np.abs(d), 95))}
    if s_comb is not None and np.isfinite(s_comb).all():
        z = d / s_comb
        out.update({"z_within_1": float(np.mean(np.abs(z) <= 1)), "z_within_2": float(np.mean(np.abs(z) <= 2)),
                    "z_std": float(np.std(z))})
    if t is not None and len(d) > 24:
        hr = (t % ts.DAY_NS) / ts.NS / 3600.0
        A = np.column_stack([np.ones_like(hr), np.cos(2 * np.pi * hr / 24), np.sin(2 * np.pi * hr / 24),
                             np.cos(2 * np.pi * hr / 12), np.sin(2 * np.pi * hr / 12)])
        c = np.linalg.lstsq(A, d, rcond=None)[0]
        out["amp24_mm"] = float(np.hypot(c[1], c[2]))
        out["amp12_mm"] = float(np.hypot(c[3], c[4]))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("ours")
    ap.add_argument("ref")
    ap.add_argument("--ref-format", choices=["tro", "csv"], default=None)
    ap.add_argument("--hourly", action="store_true", help="compare hourly means (recommended for STO vs PWC)")
    ap.add_argument("--pride", action="store_true", help="reference is a harmonised PRIDE run (Stage 3/5 gates)")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    st, t, z, s = read_ours(a.ours)
    fmt = a.ref_format or ("tro" if a.ref.upper().endswith((".TRO", ".TRO.GZ", ".TRO.Z")) else "csv")
    rt, rz, rs = read_tro(a.ref, st) if fmt == "tro" else read_ref_csv(a.ref)
    if a.hourly:
        t, z = hourly(t, z)
        s = np.full(len(t), np.nan)
    common, i1, i2 = np.intersect1d(t, rt, return_indices=True)
    if len(common) == 0:
        # nearest within 150 s
        k = np.clip(np.searchsorted(rt, t), 1, len(rt) - 1)
        kk = np.where(np.abs(rt[k] - t) < np.abs(rt[k - 1] - t), k, k - 1)
        ok = np.abs(rt[kk] - t) <= 150 * ts.NS
        i1, i2 = np.nonzero(ok)[0], kk[ok]
    d = z[i1] - rz[i2]
    sc = np.sqrt(s[i1] ** 2 + rs[i2] ** 2)
    m = metrics(d, sc, t[i1])
    m["station"] = st
    if a.pride:
        m["stage3_gate"] = bool(abs(m["bias_mm"]) <= 3 and m["std_mm"] <= 5 and m["rmse_mm"] <= 6 and m["p95_mm"] <= 10)
        m["stage5_gate"] = bool(abs(m["bias_mm"]) <= 1.5 and m["std_mm"] <= 3 and m["p95_mm"] <= 6 and
                                m.get("amp24_mm", 0) <= 1 and m.get("amp12_mm", 0) <= 1)
    else:
        m["level2_gate_igs"] = bool(abs(m["bias_mm"]) <= 3 and m["std_mm"] <= 6)
    print(json.dumps(m, indent=1))
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(m, fh, indent=1)
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        hrs = (t[i1] - t[i1][0]) / ts.NS / 3600
        fig, ax = plt.subplots(2, 1, figsize=(10, 5), sharex=True)
        ax[0].plot(hrs, z[i1], color="#2a78d6", lw=2, label="this processor")
        ax[0].plot(hrs, rz[i2], color="#eb6834", lw=2, label="reference")
        ax[0].legend(frameon=False)
        ax[0].set_ylabel("ZTD (mm)")
        ax[1].plot(hrs, d, color="#52514e", lw=1.5)
        ax[1].set_ylabel("difference (mm)")
        ax[1].set_xlabel("hours")
        for x in ax:
            x.grid(True, color="#e4e3df")
        fig.tight_layout()
        fig.savefig(a.plot, dpi=110)
    return 0


if __name__ == "__main__":
    sys.exit(main())
