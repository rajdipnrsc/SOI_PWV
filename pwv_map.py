#!/usr/bin/env python3
"""System 2 entry point: national PWV/ZWD/ZTD grids from System-1 station products (TDS Part B, Stage 8).

    python pwv_map.py RESULTS_DIR --day 2024-015 [--era5 ERA5_PL.nc] [--dem DEM.nc] [--grid G025] [--out DIR] [--cv]

Inputs: the station PWV CSV files written by pwv_ppp.py (interface contract TDS § 20); optionally an ERA5
pressure-level NetCDF (z, t, q [+ surface geopotential, lsm]) covering the day +- 1 h and the India domain, and a
NetCDF DEM (variables lat, lon, elevation). Without ERA5 the benchmark "height-normalised PWV" method is used and
labelled HEIGHT_SCALING_ONLY. EXPERIMENTAL until the § 31 spatial validation experiment has been completed:
resolution and accuracy claims are not allowed before that (TDS § 36.11, § 36.14).
"""
import argparse
import os
import sys
import time

import numpy as np

from ppp import log as plog
from ppp import mapping as mp
from ppp import products as prd
from ppp import settings
from ppp import timesys as ts

LOG = plog.get()


def main(argv=None):
    ap = argparse.ArgumentParser(description="National PWV maps from station products (System 2).")
    ap.add_argument("results", help="folder with <ST>_<YYYYDDD>_PWV.csv station products")
    ap.add_argument("--day", required=True, help="YYYY-DDD")
    ap.add_argument("--era5", default=None, help="ERA5 pressure-level NetCDF (background)")
    ap.add_argument("--dem", default=None, help="NetCDF DEM with lat, lon, elevation (orthometric m)")
    ap.add_argument("--grid", default="G025", choices=sorted(settings.MAP_GRIDS))
    ap.add_argument("--out", default=os.path.join(settings.RESULTS_DIR, "maps"))
    ap.add_argument("--slot-min", type=int, default=settings.MAP_SLOT_MIN)
    ap.add_argument("--cv", action="store_true", help="leave-one-station-out cross-validation report")
    ap.add_argument("--hourly", action="store_true", help="also write the hourly archive grid (TDS § 25)")
    ap.add_argument("--coarse", action="store_true", help="also write the 0.5 deg variance-aware aggregate (§ 27.1)")
    ap.add_argument("--support", default=settings.SUPPORT_FILE,
                    help="supportable-spacing JSON from validation/scripts/spatial_experiment.py (§ 24.4)")
    ap.add_argument("--as-of", default=None, help="use station products as available at this UTC date/time "
                    "(reproducibility, TDS § 19.3); default best_available")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)
    plog.setup("quiet" if a.quiet else "normal")
    t0 = time.time()
    try:
        y, doy = (int(x) for x in a.day.split("-"))
    except ValueError:
        LOG.error("--day must be YYYY-DDD")
        return 3
    day = ts.to_ns(y, 1, 1) + (doy - 1) * ts.DAY_NS
    tag = f"{y:04d}{doy:03d}"
    os.makedirs(a.out, exist_ok=True)
    plog.add_file(os.path.join(a.out, f"INPWV_{a.grid}_{tag}.log"))
    stations = mp.read_station_products(a.results, tag, as_of=a.as_of)
    LOG.info("Stations with usable PWV products for %s: %d", a.day, len(stations))
    if len(stations) < 3:
        LOG.error("TOO_FEW_STATIONS: at least 3 station products are needed for a map")
        return 3
    tiers = sorted({s.tier for s in stations})
    tier = "FINAL" if tiers == ["FINAL"] else "RAPID"
    bg = None
    if not a.era5 and os.path.exists(os.path.expanduser("~/.cdsapirc")):
        try:
            import importlib.util
            if importlib.util.find_spec("cdsapi") is None:
                raise ImportError("cdsapi not installed (pip install cdsapi)")
            cand = os.path.join(settings.CACHE_DIR, "ERA5", f"era5_pl_{tag}.nc")
            os.makedirs(os.path.dirname(cand), exist_ok=True)
            if not os.path.exists(cand):
                LOG.info("Downloading ERA5 background for %s (Copernicus CDS)", a.day)
                mp.download_era5(day, cand)
            a.era5 = cand
        except Exception as exc:              # noqa: BLE001 - optional source; fall back and say so
            LOG.warning("ERA5_UNAVAILABLE: %s", exc)
    if a.era5:
        bg = mp.read_era5(a.era5)
        LOG.info("Background: %s, %d times, %d levels", a.era5, len(bg.t_utc), len(bg.p_hpa))
    else:
        LOG.warning("NO_BACKGROUND: no ERA5/NWP background supplied; benchmark method HEIGHT_SCALING_ONLY "
                    "(ordinary kriging of height-normalised PWV) is used")
    dem = None
    if a.dem:
        import netCDF4
        with netCDF4.Dataset(a.dem) as nc:
            dem = (np.array(nc["lat"][:]), np.array(nc["lon"][:]), np.array(nc["elevation"][:], dtype=float))
    # station grey/blacklist from the 30-day residual history (TDS § 21 item 2; plain CSV next to the maps)
    hist_path = os.path.join(a.out, "station_residual_history.csv")
    daily = mp.station_daily_residuals(stations, day, bg)
    hist = mp.update_residual_history(hist_path, tag, daily, stations)
    grey, black, lst = mp.station_lists(hist, tag)
    if grey:
        LOG.warning("STATION_GREYLIST (review, still used): %s", ", ".join(grey))
    if black:
        LOG.warning("STATION_BLACKLIST (excluded): %s", ", ".join(f"{n} ({lst[n]['deviation_mm']:+.1f} mm vs "
                                                                f"neighbours)" for n in black))
    support, support_ver = None, "not yet determined (§ 31)"
    if a.support and os.path.exists(a.support):
        import json
        with open(a.support) as fh:
            sj = json.load(fh)
        support = {(t["lat0"], t["lon0"]): t["spacing_deg"] for t in sj.get("tiles", [])}
        support_ver = sj.get("version", os.path.basename(a.support))
        LOG.info("Supportable-spacing map %s (%d tiles)", support_ver, len(support))
    grid = mp.make_grid(a.grid, dem=dem, bg=bg)
    if dem is None:
        LOG.warning("DEM_MISSING: cell heights from %s", grid.dem_source)
    slots = np.arange(day, day + ts.DAY_NS, a.slot_min * 60 * ts.NS, dtype=np.int64)
    # REML hyperparameters pooled over hourly slots of the day (TDS § 22.4; 30-day pooling in operations)
    pool = []
    for sl in slots[::4]:
        vals = [None if s.station in black else mp.slot_values(s, sl) for s in stations]
        sts = [s for s, v in zip(stations, vals) if v is not None]
        if len(sts) < 8:
            continue
        v = np.array([vv for vv in vals if vv is not None])
        lat = np.array([s.lat for s in sts])
        lon = np.array([s.lon for s in sts])
        h = np.array([s.h_orth for s in sts])
        if bg is not None:
            b, _, _ = mp.background_at(bg, lat, lon, h, sl)
        else:
            b = mp.height_scaling_background(v[:, 0], h, h)
        pool.append((mp.to_unit(lat, lon), v[:, 0] - b, v[:, 1]))
    cp = mp.fit_reml(pool) if pool else mp.CovParams()
    LOG.info("Covariance (REML, %d slots): sigma_s %.2f mm, L %.0f km, nugget %.2f mm%s", cp.n_slots, cp.sigma_s,
             cp.L_km, cp.nugget, "" if cp.fitted else " (defaults: too few stations)")
    results = []
    for i, sl in enumerate(slots):
        results.append(mp.map_slot(grid, stations, sl, cp, bg, exclude=black))
        if i % 16 == 0:
            r = results[-1]
            LOG.info("  slot %s: %d stations, classes 1-3 %.0f %% of land, median sigma %.2f mm",
                     ts.iso(int(sl))[11:16], r.n_used, 100 * np.mean((r.cls[grid.land] >= 1) & (r.cls[grid.land] <= 3)),
                     np.nanmedian(r.sigma))
    name = f"INPWV_{a.grid}_{a.slot_min:02d}M_{tag}_{tier}_v{settings.SCHEMA_VERSION}.nc"
    chash, _ = prd.config_hash()
    attrs = {"title": f"GNSS-corrected PWV/ZWD/ZTD grid {a.day}", "summary": "Background + GNSS residual "
             "regression kriging (TDS Part B); EXPERIMENTAL until the spatial validation experiment is complete",
             "product_tier": tier, "background_source": bg.source if bg is not None else "HEIGHT_SCALING_ONLY",
             "background_assimilates_gnss": bg.assimilates_gnss if bg is not None else "n/a",
             "dem_source": grid.dem_source, "covariance": {"sigma_s_mm": cp.sigma_s, "L_km": cp.L_km,
                                                            "nugget_mm": cp.nugget, "fitted": cp.fitted},
             "station_manifests": [s.manifest_id for s in stations], "config_hash": chash,
             "station_selection": f"as_of {a.as_of}" if a.as_of else "best_available",
             "software_version": prd.software_version(), "supportable_resolution_version": support_ver,
             "station_greylist": grey, "station_blacklist": black, "date_created": prd._now_iso()}
    path = os.path.join(a.out, name)
    mp.write_grid_netcdf(path, grid, slots, results, stations, attrs, support)
    LOG.info("wrote %s", path)
    if a.coarse and a.grid != "G050":
        cg, cres = mp.aggregate(grid, results, "G050", cp)
        cpath = os.path.join(a.out, name.replace(f"INPWV_{a.grid}_", "INPWV_G050_"))
        mp.write_grid_netcdf(cpath, cg, slots, cres, stations, dict(attrs, summary=attrs["summary"] +
                             f"; variance-aware aggregation of {a.grid}"), support)
        LOG.info("wrote %s", cpath)
    if a.hourly:
        hslots = np.arange(day, day + ts.DAY_NS, 3600 * ts.NS, dtype=np.int64)
        hres = [mp.map_slot(grid, stations, sl, cp, bg, window="01H", exclude=black) for sl in hslots]
        # consistency: hourly map vs mean of the 15-min maps in [t - 30, t + 30) min (TDS § 25)
        diffs = []
        for sl, hr in zip(hslots, hres):
            sub = [r.pwv for t_, r in zip(slots, results) if sl - 1800 * ts.NS <= t_ < sl + 1800 * ts.NS]
            if sub:
                diffs.append(hr.pwv - np.nanmean(np.array(sub), axis=0))
        rms = float(np.sqrt(np.nanmean(np.array(diffs) ** 2))) if diffs else float("nan")
        LOG.info("Hourly archive vs mean of 15-min maps: RMS difference %.2f mm", rms)
        hpath = os.path.join(a.out, f"INPWV_{a.grid}_01H_{tag}_{tier}_v{settings.SCHEMA_VERSION}.nc")
        mp.write_grid_netcdf(hpath, grid, hslots, hres, stations, dict(attrs, hourly_vs_15min_rms_mm=rms), support)
        LOG.info("wrote %s", hpath)
    if a.cv:
        rows = []
        for sl in slots[::4]:
            rows += mp.loso(stations, sl, cp, bg, exclude=black)
        met = mp.cv_metrics(rows)
        LOG.info("LOSO cross-validation (hourly slots): %s", {k: round(v, 3) for k, v in met.items()})
        prd.write_json(os.path.join(a.out, f"INPWV_{a.grid}_{tag}_cv.json"), {"metrics": met, "rows": rows})
    LOG.info("done in %.0f s", time.time() - t0)
    return 0


if __name__ == "__main__":
    sys.exit(main())
