#!/usr/bin/env python3
"""Spatial validation experiment for System 2 (TDS § 31) and supportable-resolution map (§ 24.4).

    python validation/scripts/spatial_experiment.py RESULTS_DIR --days 2024-001:2024-366 [--era5-dir DIR]
           [--every-h 3] [--realisations 10] [--out validation/stage8_spatial] [--write-support]

For each day and every --every-h hours: LOSO, spatial block CV (settings.BLOCK_CV_KM), network thinning to
settings.THINNING_SPACINGS_KM (random and stratified, settings.THINNING_REALISATIONS each) and the temporal designs
(5/15/30/60-min station windows predicting the withheld stations' 5-min values).  Results are stratified by regime
(JF, MAM, JJAS, OND; dry = PWV < 15 mm) and by distance to the nearest station; E(d) is fitted per regime; the
supportable grid spacing is evaluated per 1-degree tile (criteria 1-4 of § 24.2; criterion 5 is reported as not
evaluated).  ERA5 files are taken from --era5-dir/era5_pl_<YYYYDDD>.nc; without them the HEIGHT_SCALING_ONLY
benchmark background is used and the report says so.

Outputs: report.md, results.json, and (with --write-support) supportable_spacing.json read by pwv_map.py.
The acceptance decisions of § 31 (sigma calibration, skill, national grid spacing) are taken by people from this
report and recorded in docs/decisions.md - the script never claims them.
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from ppp import mapping as mp  # noqa: E402
from ppp import products as prd  # noqa: E402
from ppp import settings  # noqa: E402
from ppp import spatial_validation as sv  # noqa: E402
from ppp import timesys as ts  # noqa: E402

REGIMES = {1: "JF", 2: "JF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJAS", 7: "JJAS", 8: "JJAS", 9: "JJAS",
           10: "OND", 11: "OND", 12: "OND"}


def day_list(spec):
    a, b = spec.split(":") if ":" in spec else (spec, spec)
    y0, d0 = (int(x) for x in a.split("-"))
    y1, d1 = (int(x) for x in b.split("-"))
    t0 = ts.to_ns(y0, 1, 1) + (d0 - 1) * ts.DAY_NS
    t1 = ts.to_ns(y1, 1, 1) + (d1 - 1) * ts.DAY_NS
    return list(range(t0, t1 + 1, ts.DAY_NS))


def run(argv=None):
    ap = argparse.ArgumentParser(description="TDS § 31 spatial validation experiment")
    ap.add_argument("results")
    ap.add_argument("--days", required=True, help="YYYY-DDD or YYYY-DDD:YYYY-DDD")
    ap.add_argument("--era5-dir", default=os.path.join(settings.CACHE_DIR, "ERA5"))
    ap.add_argument("--every-h", type=int, default=3)
    ap.add_argument("--realisations", type=int, default=settings.THINNING_REALISATIONS)
    ap.add_argument("--out", default=os.path.join("validation", "stage8_spatial"))
    ap.add_argument("--write-support", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    rows = {}                      # design -> regime -> list of rows
    temporal = {}                  # window_min -> list of rows
    loso_all, last_sd, cps, bg_sources = [], None, [], set()

    def add(design, regime, rr):
        rows.setdefault(design, {}).setdefault(regime, []).extend(rr)
        rows[design].setdefault("ALL", []).extend(rr)
    for d in day_list(a.days):
        y, doy = ts.year_doy(d)
        tag = f"{y:04d}{doy:03d}"
        stations = mp.read_station_products(a.results, tag)
        if len(stations) < 10:
            continue
        bg = None
        f = os.path.join(a.era5_dir, f"era5_pl_{tag}.nc")
        if os.path.exists(f):
            bg = mp.read_era5(f)
        bg_sources.add(bg.source if bg is not None else "HEIGHT_SCALING_ONLY")
        regime = REGIMES[ts.from_ns(d)[1]]
        pool = []
        slots = [d + h * 3600 * ts.NS for h in range(0, 24, a.every_h)]
        sds = [(sl, sv.slot_data(stations, sl, bg)) for sl in slots]
        for sl, sd in sds:
            if sd is not None and len(sd.names) >= 8:
                pool.append((mp.to_unit(sd.lat, sd.lon), sd.pwv - sd.bg, sd.sig))
        cp = mp.fit_reml(pool) if pool else mp.CovParams()
        cps.append(cp)
        for k, (sl, sd) in enumerate(sds):
            if sd is None or len(sd.names) < 10:
                continue
            last_sd = sd
            rg = [regime] + (["DRY"] if np.median(sd.pwv) < 15 else [])
            n = len(sd.names)
            lo = []
            for i in range(n):
                test = np.arange(n) == i
                lo += sv.predict_withheld(sd, ~test, test, cp)
            for r_ in rg:
                add("LOSO", r_, lo)
            loso_all += lo
            for bk in settings.BLOCK_CV_KM:
                br = sv.block_cv(sd, cp, bk)
                for r_ in rg:
                    add(f"BLOCK_{bk}KM", r_, br)
            for spc in settings.THINNING_SPACINGS_KM:
                for mode in ("random", "stratified"):
                    tr = sv.thinning_cv(sd, cp, spc, a.realisations, mode, seed=int(sl // ts.NS) % 100000 + spc)
                    for r_ in rg:
                        add(f"THIN_{spc}KM_{mode}", r_, tr)
            # temporal designs: train with W-min windows, verify against the withheld station's 5-min value
            if k % 2 == 0:
                truth = sv.slot_data(stations, sl, bg, half_ns=int(2.5 * 60) * ts.NS, min_n=1)
                for w in (5, 15, 30, 60):
                    sw = sv.slot_data(stations, sl, bg, half_ns=int(w * 30) * ts.NS, min_n=max(1, int(w / 5 * 2 / 3)))
                    if sw is None or truth is None:
                        continue
                    common = [nm for nm in truth.names if nm in sw.names]
                    for nm in common[::3]:
                        i = sw.names.index(nm)
                        test = np.arange(len(sw.names)) == i
                        rr = sv.predict_withheld(sw, ~test, test, cp)
                        if rr:
                            j = truth.names.index(nm)
                            r0 = rr[0]
                            temporal.setdefault(w, []).append((nm, float(truth.pwv[j]), r0[2], r0[3], r0[4], r0[5]))
        print(f"{tag}: {len(stations)} stations, regime {regime}, L {cp.L_km:.0f} km")
    if not loso_all:
        print("No day with >= 10 station products: nothing to evaluate")
        return 1
    metrics = {des: {rg: sv.metrics(rr) for rg, rr in by.items()} for des, by in rows.items()}
    all_rows = [r for des, by in rows.items() for r in by.get("ALL", [])]
    fits = {"ALL": sv.fit_error_distance(all_rows)}
    for rg in sorted({rg for by in rows.values() for rg in by} - {"ALL"}):
        fits[rg] = sv.fit_error_distance([r for by in rows.values() for r in by.get(rg, [])])
    tmet = {w: sv.metrics(rr) for w, rr in sorted(temporal.items())}
    L_med = float(np.median([c.L_km for c in cps]))
    cp_med = mp.CovParams(float(np.median([c.sigma_s for c in cps])), L_med, float(np.median([c.nugget for c in cps])))
    tiles = sv.support_by_tile(last_sd, cp_med, fits["ALL"], loso_all) if fits["ALL"].get("fitted") else {}
    res = {"days": a.days, "background": sorted(bg_sources), "covariance_median": vars(cp_med),
           "metrics": metrics, "error_distance_fits": fits, "temporal": tmet,
           "support_tiles": [{"lat0": k[0], "lon0": k[1], **v} for k, v in tiles.items()],
           "software_version": prd.software_version(), "config_hash": prd.config_hash()[0]}
    prd.write_json(os.path.join(a.out, "results.json"), res)
    if a.write_support and tiles:
        prd.write_json(os.path.join(a.out, "supportable_spacing.json"),
                       {"version": f"provisional {a.days} ({prd.software_version()})", "criteria": "TDS § 24.2 1-4",
                        "tiles": [{"lat0": k[0], "lon0": k[1], "spacing_deg": v["spacing_deg"]}
                                  for k, v in tiles.items()]})
    write_report(os.path.join(a.out, "report.md"), res)
    print(f"wrote {a.out}/report.md")
    return 0


def write_report(path, res):
    L = ["# Spatial validation experiment (TDS § 31) — report", "",
         f"Days: {res['days']}; background: {', '.join(res['background'])}; software {res['software_version']}; "
         f"config {res['config_hash'][:12]}.", "",
         "Withheld-station metrics on PWV (mm). Skill = 1 − MSE/MSE_background. Coverage = fraction within ±1σ "
         "(acceptance § 31 (i): 60–76 %).", "",
         "| Design | Regime | n | RMSE | bias | MAE | P95 | skill | ±1σ | ±2σ |", "|---|---|---|---|---|---|---|---|---|---|"]
    for des, by in res["metrics"].items():
        for rg, m in by.items():
            if m.get("n"):
                L.append(f"| {des} | {rg} | {m['n']} | {m['rmse']:.2f} | {m['bias']:+.2f} | {m['mae']:.2f} | "
                         f"{m['p95']:.2f} | {m['skill']:.2f} | {m['within_1sigma']:.2f} | {m['within_2sigma']:.2f} |")
    L += ["", "## Error vs distance to the nearest station, E(d) = √(E0² + (a·d^α)²)", "",
          "| Regime | E0 | a | α | bins (km) | RMSE (mm) |", "|---|---|---|---|---|---|"]
    for rg, f in res["error_distance_fits"].items():
        if f.get("fitted"):
            L.append(f"| {rg} | {f['E0']:.2f} | {f['a']:.3f} | {f['alpha']:.2f} | "
                     f"{', '.join(f'{x:.0f}' for x in f['bins_km'])} | {', '.join(f'{x:.2f}' for x in f['rmse_mm'])} |")
    L += ["", "## Temporal designs (station window → withheld 5-min value)", "", "| Window (min) | n | RMSE | skill |",
          "|---|---|---|---|"]
    for w, m in res["temporal"].items():
        if m.get("n"):
            L.append(f"| {w} | {m['n']} | {m['rmse']:.2f} | {m['skill']:.2f} |")
    L += ["", "## Supportable grid spacing per 1° tile (criteria 1–4 of § 24.2; criterion 5 not evaluated)", ""]
    cnt = {}
    for t in res["support_tiles"]:
        cnt[t["spacing_deg"]] = cnt.get(t["spacing_deg"], 0) + 1
    L += [f"- {k if k is not None else 'background only'}: {v} tiles" for k, v in sorted(cnt.items(), key=str)]
    L += ["", "Decisions (national grid spacing, coverage-class thresholds, temporal resolution) are taken from this "
          "report and recorded in `docs/decisions.md`; grid spacing is not atmospheric resolution."]
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    sys.exit(run())
