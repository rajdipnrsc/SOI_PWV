#!/usr/bin/env python3
"""GPS PPP -> 5-min ZTD -> PWV for one station-day (System 1 entry point, TDS § 33.3).

    python pwv_ppp.py SITE.o SITE.n
    python pwv_ppp.py SITE.o SITE.n --xyz X Y Z      (SOI coordinate, ITRF2008 @ 2005.0, metres)

Everything else (CODE products, ANTEX, VMF3, GPT3, leap seconds) is downloaded, verified and cached
automatically.  Unknown station metadata never stop a run (TDS § 0.5, § 9.11): defaults are applied, a
WARNING is printed, and the assumption is recorded in the flags and the manifest.
Exit codes: 0 = products written, 2 = WAIT_FOR_TIER, 3 = FAIL (no PPP possible), 4 = unexpected error.
"""
import argparse
import copy
import csv
import glob
import os
import sys
import time
import traceback

import numpy as np

from ppp import ambiguity as amb
from ppp import antenna as an
from ppp import coords as co
from ppp import estimator as est
from ppp import log as plog
from ppp import model as mdl
from ppp import output as out
from ppp import preprocess as pp
from ppp import products as prd
from ppp import pwv as pw
from ppp import rinex
from ppp import settings
from ppp import timesys as ts
from ppp import troposphere as tr
from ppp.corrections import read_blq

LOG = plog.get()


class Stop(Exception):
    """Clean end of a run (status code + plain-language message), never a traceback on screen."""

    def __init__(self, status, message, code=3):
        super().__init__(message)
        self.status, self.message, self.code = status, message, code


# ============================================================================================ inputs
def parse_args(argv=None):
    p = argparse.ArgumentParser(description="GPS PPP -> 5-min ZTD/PWV (CODE products, float PPP).")
    p.add_argument("obs", help="RINEX 3 observation file (.o/.rnx, .crx/.gz accepted)")
    p.add_argument("nav", help="RINEX navigation file (sanity checks, satellite health)")
    p.add_argument("--xyz", nargs=3, type=float, metavar=("X", "Y", "Z"),
                   help="SOI coordinate (m), ITRF2008 @ 2005.0; default from stations.csv if present")
    p.add_argument("--mode", default="auto", choices=["auto", "fixed", "constrained", "static"])
    p.add_argument("--out", default=settings.RESULTS_DIR)
    p.add_argument("--cutoff", type=float, default=settings.ELEVATION_CUTOFF_DEG)
    p.add_argument("--ar", action="store_true", help="PPP-AR layer (Stage 6; separate FIXED output)")
    p.add_argument("--no-pwv", action="store_true")
    p.add_argument("--offline", action="store_true", help="use the product cache only")
    p.add_argument("--proxy", default=None)
    g = p.add_mutually_exclusive_group()
    g.add_argument("--quiet", action="store_true")
    g.add_argument("--verbose", action="store_true")
    return p.parse_args(argv)


def lookup_station_xyz(station):
    """SOI coordinate from stations.csv / station_coord.csv (flexible column names)."""
    for fn in settings.STATIONS_FILES:
        for base in (os.getcwd(), os.path.dirname(os.path.abspath(__file__))):
            path = os.path.join(base, fn)
            if not os.path.exists(path):
                continue
            with open(path, newline="") as fh:
                rd = csv.DictReader(fh)
                cols = {c.strip().lower(): c for c in (rd.fieldnames or []) if c}
                kid = next((cols[c] for c in ("station", "stationcode", "site", "code", "name", "stationname")
                            if c in cols), None)
                kx = next((cols[c] for c in ("x", "ecef_x", "x_m") if c in cols), None)
                ky = next((cols[c] for c in ("y", "ecef_y", "y_m") if c in cols), None)
                kz = next((cols[c] for c in ("z", "ecef_z", "z_m") if c in cols), None)
                if not (kid and kx and ky and kz):
                    continue
                for row in rd:
                    if (row.get(kid) or "").strip().upper()[:4] == station[:4].upper():
                        try:
                            return np.array([float(row[kx]), float(row[ky]), float(row[kz])]), path
                        except (TypeError, ValueError):
                            continue
    return None, None


def find_adjacent(path, station, day_ns):
    """RINEX files of the same station for D-1 and D+1 in the same folder (for the 30-h window, TDS § 5.10)."""
    d = os.path.dirname(os.path.abspath(path))
    out_ = {}
    for off in (-1, 1):
        dd = day_ns + off * ts.DAY_NS
        y, doy = ts.year_doy(dd)
        pats = [f"{station[:4]}{doy:03d}?.{y % 100:02d}[oOdD]*", f"{station[:4].lower()}{doy:03d}?.{y % 100:02d}[oOdD]*",
                f"{station[:4].upper()}*_{y:04d}{doy:03d}0000_01D_*MO.*"]
        for pat in pats:
            hits = sorted(glob.glob(os.path.join(d, pat)))
            if hits:
                out_[off] = hits[0]
                break
    return out_


def find_met(path, station, day_ns):
    """RINEX meteorological files of the station for D-1, D, D+1 in the folder of SITE.o (TDS § 18.2 priority 1)."""
    d = os.path.dirname(os.path.abspath(path))
    hits = []
    for off in (-1, 0, 1):
        y, doy = ts.year_doy(day_ns + off * ts.DAY_NS)
        for pat in settings.MET_FILE_PATTERNS:
            g = sorted(glob.glob(os.path.join(d, pat.format(st=station[:4].upper(), stl=station[:4].lower(), y=y,
                                                            yy=y % 100, doy=doy))))
            if g:
                hits.append(g[0])
                break
    return hits


def load_barometer(obs_path, station, day, h_ant_orth, undu, t_fallback_k):
    """Station barometer reduced to the ARP from auto-detected RINEX met files; None (with WARNING) if unusable."""
    files = find_met(obs_path, station, day)
    if not files:
        return None, []
    mets = []
    for f in files:
        try:
            mets.append(rinex.read_met(f))
        except (rinex.RinexError, OSError, ValueError) as exc:
            LOG.warning("MET_UNREADABLE: %s: %s", os.path.basename(f), exc)
    if not mets:
        return None, []
    met = mets[0]
    if len(mets) > 1:
        tt = np.concatenate([m.t for m in mets])
        keys = set.intersection(*(set(m.values) for m in mets))
        vals = {k: np.concatenate([m.values[k] for m in mets]) for k in keys}
        o = np.unique(tt, return_index=True)[1]
        met = rinex.MetData(tt[o], {k: v[o] for k, v in vals.items()}, mets[0].sensor_pos, mets[0].sensor_info,
                            mets[0].version, mets[0].path)
    baro, warns = pw.barometer_from_met(met, h_ant_orth, undu, t_fallback_k)
    for w in warns:
        LOG.warning("%s: station met file %s (%s)", w, os.path.basename(met.path), {
            "MET_NO_PRESSURE": "no PR observable", "MET_NO_VALID_PRESSURE": "no plausible pressure values",
            "MET_HEIGHT_ASSUMED": "no SENSOR POS XYZ/H for PR: sensor assumed at antenna height (§ 9.11 style "
                                  "default)", "MET_HEIGHT_NO_GEOID": "no geoid undulation: sensor height used as "
                                                                    "orthometric",
            "MET_SENSOR_TOO_FAR": f"sensor more than {settings.MET_MAX_HEIGHT_DIFF_M:.0f} m from the antenna height; "
                                  "barometer not used",
            "MET_NO_TEMPERATURE": "no TD observable: height reduction with GPT3/standard temperature"}.get(w, ""))
    if baro is not None:
        LOG.info("Station barometer: %s, %d values (%d rejected), sensor H %.1f m (%s) -> ARP H %.1f m",
                 baro.info, len(baro.t), baro.n_rejected, baro.h_sensor, baro.height_source, baro.h_ant)
    return baro, [os.path.basename(f) for f in files] + warns


def era5_for_day(day):
    """ERA5 pressure-level file for the UTC day from the cache (or CDS when enabled); None if unavailable."""
    if not settings.STATION_ERA5:
        return None
    y, doy = ts.year_doy(day)
    cand = os.path.join(settings.CACHE_DIR, "ERA5", f"era5_pl_{y:04d}{doy:03d}.nc")
    if not os.path.exists(cand) and settings.STATION_ERA5_DOWNLOAD and os.path.exists(os.path.expanduser(
            "~/.cdsapirc")):
        try:
            from ppp import mapping as mp
            os.makedirs(os.path.dirname(cand), exist_ok=True)
            LOG.info("Downloading ERA5 for station PWV (Copernicus CDS; may queue)")
            mp.download_era5(day, cand)
        except Exception as exc:              # noqa: BLE001 - optional source; fall back and say so
            LOG.warning("ERA5_UNAVAILABLE: %s", exc)
    if not os.path.exists(cand):
        return None
    try:
        from ppp import mapping as mp
        return mp.read_era5(cand)
    except Exception as exc:                  # noqa: BLE001
        LOG.warning("ERA5_UNREADABLE: %s: %s", os.path.basename(cand), exc)
        return None


def merge_obs(parts, t0, t1):
    """Merge RINEX ObsData objects restricted to [t0, t1] GPST (window with overlap)."""
    sats = sorted({s for p in parts for s in p.sats})
    codes = sorted({c for p in parts for c in p.obs})
    ep = []
    for p in parts:
        ep.extend(int(t) for t in p.epochs if t0 <= t <= t1)
    ep = np.array(sorted(set(ep)), dtype=np.int64)
    ei = {int(t): i for i, t in enumerate(ep)}
    si = {s: j for j, s in enumerate(sats)}
    obs = {c: np.full((len(ep), len(sats)), np.nan) for c in codes}
    lli = {c: np.zeros((len(ep), len(sats)), dtype=np.int8) for c in codes}
    ssi = {c: np.zeros((len(ep), len(sats)), dtype=np.int8) for c in codes}
    clk = np.full(len(ep), np.nan)
    for p in parts:
        rows = np.array([ei.get(int(t), -1) for t in p.epochs])
        keep = rows >= 0
        cols = np.array([si[s] for s in p.sats])
        for c in p.obs:
            obs[c][np.ix_(rows[keep], cols)] = p.obs[c][keep]
            lli[c][np.ix_(rows[keep], cols)] = p.lli[c][keep]
            ssi[c][np.ix_(rows[keep], cols)] = p.snr_flag[c][keep]
        clk[rows[keep]] = p.clock_offset[keep]
    main = parts[0]
    return rinex.ObsData(main.header, ep, sats, obs, lli, ssi, clk, np.zeros(len(ep), dtype=np.int8), main.path)


# ============================================================================================ run
def run(args, products_override=None, tropo_override=None, atx_override=None):
    t_start = time.time()
    plog.setup("quiet" if args.quiet else ("verbose" if args.verbose else "normal"))
    applied = settings.apply_overrides()
    os.makedirs(args.out, exist_ok=True)
    warnings_flags = set()
    manifest_extra = {}

    # ---------------------------------------------------------------- RINEX
    if not os.path.exists(args.obs):
        raise Stop("FAIL", f"RINEX_UNREADABLE: observation file not found: {args.obs}")
    try:
        obs_d = rinex.read_obs(args.obs, decimate_s=settings.PROCESSING_INTERVAL_S)
    except (rinex.RinexError, OSError, ValueError) as exc:
        raise Stop("FAIL", f"RINEX_UNREADABLE: {exc}")
    hdr = obs_d.header
    station = (hdr.marker_name or os.path.basename(args.obs))[:4].upper()
    if len(obs_d.epochs) == 0:
        raise Stop("FAIL", "NO_GPS_DATA: no GPS epochs in the RINEX file")
    first, last = int(obs_d.epochs[0]), int(obs_d.epochs[-1])
    day = ts.day_start(first + (last - first) // 2)
    y, doy = ts.year_doy(day)
    tag = f"{station}_{y:04d}{doy:03d}"
    logpath = os.path.join(args.out, tag + ".log")
    plog.add_file(logpath)
    LOG.info("%s %s  station %s  day %04d-%03d", settings.SOFTWARE_NAME, prd.software_version(), station, y, doy)
    if applied:
        LOG.info("Settings overridden from file: %s", ", ".join(applied))
    gps_sats = [s for s in obs_d.sats if s.startswith("G")]
    LOG.info("RINEX %s: %s, %.2f, %s -> %s (GPST), %d epochs @ %.0f s (file interval %.0f s), %d GPS satellites",
             os.path.basename(args.obs), hdr.receiver, hdr.version, ts.iso(first), ts.iso(last), len(obs_d.epochs),
             settings.PROCESSING_INTERVAL_S, hdr.interval or obs_d.interval_file, len(gps_sats))
    LOG.info("Antenna %r, DELTA H/E/N = %.4f %.4f %.4f m", hdr.antenna_type, *hdr.delta_hen)
    if hdr.rcv_clock_offs_appl:
        LOG.warning("RCV_CLOCK_OFFS_APPL: receiver applied its clock offset (RINEX flag 1); epoch tags are treated "
                    "as corrected reception times")
    nav = None
    try:
        nav = rinex.read_nav(args.nav)
        LOG.info("Navigation file: %d GPS ephemerides for %d satellites", sum(len(v) for v in nav.values()), len(nav))
    except (rinex.RinexError, OSError, ValueError) as exc:
        LOG.warning("NAV_UNREADABLE: %s (continuing; precise products are used for PPP)", exc)

    # ---------------------------------------------------------------- window (30 h if adjacent RINEX exists)
    win0, win1 = first, last
    adj = find_adjacent(args.obs, station, day)
    edge_ok = [False, False]
    parts = [obs_d]
    for off, pth in sorted(adj.items()):
        try:
            o2 = rinex.read_obs(pth, decimate_s=settings.PROCESSING_INTERVAL_S)
            parts.append(o2)
            if off < 0:
                win0 = max(int(o2.epochs[0]), day - int(settings.OVERLAP_H * 3600) * ts.NS)
                edge_ok[0] = True
            else:
                win1 = min(int(o2.epochs[-1]), day + ts.DAY_NS + int(settings.OVERLAP_H * 3600) * ts.NS)
                edge_ok[1] = True
            LOG.info("Adjacent-day RINEX used for the overlap window: %s", os.path.basename(pth))
        except (rinex.RinexError, OSError, ValueError) as exc:
            LOG.warning("ADJACENT_RINEX_UNREADABLE: %s (%s)", pth, exc)
    if len(parts) > 1:
        obs_d = merge_obs(parts, win0, win1)
    else:
        LOG.info("No adjacent-day RINEX: window = data span; outputs near the day edges get the EDGE flag")
    epochs = obs_d.epochs

    # ---------------------------------------------------------------- products
    dl = prd.Downloader(settings.CACHE_DIR, offline=args.offline, proxy=args.proxy)
    leap = prd.get_leap_seconds(dl)
    if leap["source"] == "builtin":
        if ts.builtin_table_valid(last):
            LOG.warning("LEAPSEC_BUILTIN: leap-second file not downloadable; built-in table used (valid through %s)",
                        "-".join(map(str, settings.LEAP_SECONDS_BUILTIN_VALID_UNTIL)))
        else:
            LOG.warning("LEAPSEC_BUILTIN_EXPIRED: built-in leap-second table may be out of date for this date")
    sel_probe = pp.select_observables(obs_d)
    codes_used = sorted({n[0] for n in sel_probe.pair_names.values()} | {n[2] for n in sel_probe.pair_names.values()})
    if not codes_used:
        raise Stop("FAIL", "NO_DUAL_FREQUENCY_GPS: no usable dual-frequency GPS code+phase observations")
    LOG.info("Observables selected (priority lists): %s", "; ".join(
        f"{a}/{b}+{c}/{d}" for (a, b, c, d) in sel_probe.pair_names.values()))
    if products_override is not None:
        ps = products_override
        LOG.warning("PRODUCTS_OVERRIDE: products supplied by the caller (%s) - test/diagnostic run only",
                    getattr(ps, "family", "?"))
    else:
        LOG.info("Product search (CODE only: %s)", " -> ".join(settings.PRODUCT_FAMILIES))
        ps = prd.resolve_products(dl, win0, win1, code_obs_used=codes_used, ar_requested=args.ar)
        if dl.dead_hosts:
            manifest_extra["unreachable_hosts"] = dl.dead_hosts
    LOG.info("Product status: %s %s family=%s tier=%s", ps.status, ("(" + ps.reason + ")") if ps.reason else "",
             ps.family or "-", ps.tier or "-")
    if ps.status == "WAIT_FOR_TIER":
        raise Stop("WAIT_FOR_TIER", "WAIT_FOR_TIER: no PPP-capable CODE product set yet (" + ps.reason +
                   "); rerun later", code=2)
    if ps.status == "FAIL":
        hint = ""
        src_hosts = {u.split("/")[2] for u in settings.PRODUCT_SOURCES}
        if dl.dead_hosts:
            hint = " Unreachable hosts: " + ", ".join(dl.dead_hosts) + " (check network/proxy access)."
        if src_hosts and src_hosts <= set(dl.dead_hosts):
            raise Stop("FAIL", "PRODUCTS_UNREACHABLE: none of the CODE product sources could be contacted." + hint +
                       " Products already in the cache would be used with --offline.")
        raise Stop("FAIL", ps.reason + ". Missing: " + ", ".join(ps.missing[:12]) + hint)
    warnings_flags |= set(ps.flags)
    # frame label
    family = co.label_to_family(ps.frame_label)
    LOG.info("Product frame label %r -> %s; clock header ANTEX %r", ps.frame_label, family or "UNKNOWN",
             ps.antex_name or "-")

    # ---------------------------------------------------------------- ANTEX
    if atx_override is not None:
        atx, atx_meta, atx_status = atx_override, {"sha256": None}, "OVERRIDE"
    else:
        atx_path, atx_meta, atx_status = prd.get_antex(dl, ps.antex_name)
        if atx_path is None:
            raise Stop("FAIL", f"ANTEX_UNAVAILABLE: ANTEX {ps.antex_name or settings.ANTEX_DEFAULT_NAME} not obtainable")
        atx = an.read_antex(atx_path)
        if atx_status == "FAMILY_FALLBACK":
            LOG.warning("ANTEX_VERSION_MISMATCH: clock header names %s, using %s (same family; docs/decisions.md D-004)",
                        ps.antex_name, atx.name)
    rcv_ant, ant_status = an.receiver_entry(atx, hdr.antenna_type, hdr.antenna_serial)
    if ant_status == "RADOME_FALLBACK":
        LOG.warning("RADOME_FALLBACK: %r not in %s; type+NONE calibration used", hdr.antenna_type, atx.name)
        warnings_flags.add("RADOME_FALLBACK")
    elif ant_status == "NO_ANTENNA_CALIBRATION":
        LOG.warning("NO_ANTENNA_CALIBRATION: %r not in %s; PCO/PCV = 0 and mode forced to static",
                    hdr.antenna_type, atx.name)
        warnings_flags.add("NO_ANTENNA_CALIBRATION")
    else:
        LOG.info("Receiver antenna calibration: %s (%s)", rcv_ant.type, ant_status)

    # ---------------------------------------------------------------- station metadata, header registry
    prev = _previous_manifest(args.out, station, day)
    if prev and (prev.get("antenna_type") != hdr.antenna_type or
                 prev.get("delta_hen") != list(hdr.delta_hen)):
        LOG.warning("HEADER_REGISTRY_MISMATCH: antenna/eccentricity differ from the previous run of %s", station)
        warnings_flags.add("HEADER_REGISTRY_MISMATCH")

    # ---------------------------------------------------------------- observables, OSB, IF
    sel = pp.select_observables(obs_d)
    t_mid = int(day + ts.DAY_NS // 2)
    osb_check = pp.verify_osb_sign(obs_d, ps.osb, t_mid)
    if osb_check:
        LOG.info("OSB sign check (C1C/C1W both observed): %s", osb_check)
        manifest_extra["osb_sign_check"] = osb_check
    P1c, P2c, L1c, L2c, osb_info = pp.apply_osb(sel, ps.osb, t_mid)
    if osb_info["missing"] and prd.needs_osb(codes_used):
        LOG.warning("CODE_BIAS_MISSING: OSB missing for %d satellite signals (e.g. %s)", len(osb_info["missing"]),
                    ", ".join(osb_info["missing"][:5]))
    sel_c = copy.copy(sel)
    sel_c.P1, sel_c.P2 = P1c, P2c
    clk_jumps = pp.detect_clock_jumps(sel_c, np.isfinite(sel_c.P1))
    if clk_jumps:
        LOG.warning("CLOCK_JUMP: %d receiver clock jumps detected (%s)", len(clk_jumps),
                    ", ".join(f"epoch {k}: {m:+.3f} ms {'repaired' if r else 'NOT repaired'}" for k, m, r in clk_jumps[:5]))
    Pif, Lif, a1, a2 = pp.if_combination(sel_c)

    # ---------------------------------------------------------------- a priori receiver clock and position (SPP)
    Xapprox = np.array(hdr.approx_xyz)
    X_spp, rclk_spp, spp_rms = mdl.spp(sel_c, epochs, ps, Xapprox)
    if X_spp is None:
        raise Stop("FAIL", "SPP_FAILED: no single-point solution possible (too few satellites with products)")
    LOG.info("Single-point solution: %.3f %.3f %.3f (rms %.2f m); |SPP - RINEX approx| = %.1f m", *X_spp, spp_rms,
             np.linalg.norm(X_spp - Xapprox) if np.linalg.norm(Xapprox) > 0 else -1)

    # ---------------------------------------------------------------- troposphere a priori (needs lat/lon/h)
    lat_s, lon_s, h_s = co.ecef_to_geodetic(X_spp)
    gpt3_grid = None
    if tropo_override is not None:
        vmf_files, oro, gpt3_grid = tropo_override
    else:
        g3path, _m = prd.get_gpt3(dl)
        if g3path:
            try:
                gpt3_grid = tr.read_gpt3(g3path)
            except (OSError, ValueError) as exc:
                LOG.warning("GPT3_UNREADABLE: %s", exc)
        vmf_files, vmf_missing, oro = prd.get_vmf3(dl, win0, win1)
        if vmf_missing:
            LOG.warning("VMF3_MISSING: %d of %d VMF3 grid files not available (%s)", len(vmf_missing),
                        len(vmf_missing) + len(vmf_files), ", ".join(vmf_missing[:3]))
            vmf_files = []
    undu = None
    if gpt3_grid is not None:
        undu = float(tr.gpt3(gpt3_grid, float(ts.mjd(t_mid)), lat_s, lon_s, h_s)["undu"])

    # ---------------------------------------------------------------- coordinates (TDS § 9.10-9.12)
    xyz_soi, src_file = (np.array(args.xyz), "--xyz") if args.xyz else lookup_station_xyz(station)
    mode = args.mode
    coord_soi = None
    if xyz_soi is None:
        LOG.warning("NO_SOI_COORD: no SOI coordinate (--xyz or stations.csv); static mode with a priori from the "
                    "single-point solution")
        if mode in ("fixed", "constrained"):
            LOG.warning("MODE_CHANGED: --mode %s needs an SOI coordinate; using static", mode)
        mode = "static"
    else:
        LOG.info("SOI coordinate (%s): %.3f %.3f %.3f  [ITRF2008 @ 2005.0]", src_file, *xyz_soi)
        if family is None:
            LOG.warning("FRAME_UNKNOWN: product frame label %r not in the table; SOI coordinate cannot be "
                        "transformed", ps.frame_label)
            if mode in ("fixed", "constrained"):
                raise Stop("FAIL", f"FRAME_UNKNOWN: product frame label {ps.frame_label!r} (fixed/constrained mode)")
            mode = "static"
        else:
            lat0, lon0, _ = co.ecef_to_geodetic(xyz_soi)
            zone = co.in_deforming_zone(np.degrees(lat0), np.degrees(lon0))
            LOG.warning("VELOCITY_ASSUMED_PMM: SOI velocity unknown; ITRF2020 plate motion model (India plate + ORB), "
                        "vertical 0 +- %.0f mm/yr", settings.SIGMA_VU_MM_YR)
            if zone:
                LOG.warning("NONLINEAR_MOTION_REGION: station in %s; plate model not valid there", zone)
            LOG.warning("TIDE_SYSTEM_ASSUMED: SOI coordinate assumed conventional tide-free")
            if hdr.delta_hen[0] != 0.0 or hdr.delta_hen[1] != 0.0 or hdr.delta_hen[2] != 0.0:
                LOG.warning("REFPOINT_ASSUMED_MARKER: SOI coordinate assumed to refer to the marker; RINEX "
                            "DELTA H/E/N applied to obtain the ARP")
            coord_soi = co.soi_to_processing(xyz_soi, t_mid, ps.frame_label, hdr.delta_hen, "MARKER", undu)
            coord_soi_arp = co.soi_to_processing(xyz_soi, t_mid, ps.frame_label, hdr.delta_hen, "ARP", undu)
            d = co.discrepancy_enu(coord_soi.xyz, X_spp)
            LOG.info("SOI -> %s (%s) @ %.4f: ARP %.3f %.3f %.3f (moved %.3f m since 2005.0); SOI - SPP (E,N,U) = "
                     "%.2f %.2f %.2f m", family, ps.frame_label, coord_soi.epoch, coord_soi.X, coord_soi.Y,
                     coord_soi.Z, np.linalg.norm(coord_soi.xyz - xyz_soi), *d)
    if ant_status == "NO_ANTENNA_CALIBRATION" and mode != "static":
        mode = "static"
    if not ps.osb and prd.needs_osb(codes_used):
        warnings_flags.add("CODE_BIAS_MISSING")

    # ---------------------------------------------------------------- troposphere + OTL
    X_ap = coord_soi.xyz if coord_soi is not None else X_spp
    lat, lon, h = co.ecef_to_geodetic(X_ap)
    tropo = tr.build_station_tropo(lat, lon, h, win0, win1, vmf_files, oro, gpt3_grid)
    h_orth_ap = h - (undu if undu is not None else 0.0)
    t_fb = tropo.met.get("T") + 273.15 if tropo.met.get("T") is not None else None
    baro, met_info = load_barometer(args.obs, station, day, h_orth_ap, undu, t_fb)
    if met_info:
        manifest_extra["met_files"] = met_info
    if baro is not None:
        grid5 = np.arange(win0 - win0 % (300 * ts.NS), win1 + 300 * ts.NS, 300 * ts.NS, dtype=np.int64)
        p5 = pw.pressure_at(baro, grid5)
        p_ref = tr.saastamoinen_pressure(tropo.at(grid5)[0], lat, h_orth_ap)
        bias = float(np.nanmedian(p5 - p_ref)) if np.isfinite(p5).any() else np.nan
        lim = settings.MET_MAX_BIAS_HPA * (1 if tropo.zhd_source == "VMF3_GRID" else 2)
        if not np.isfinite(bias) or abs(bias) > lim:
            LOG.warning("BAROMETER_REJECTED: median barometer - %s pressure = %.1f hPa (limit %.0f hPa): wrong "
                        "calibration or sensor height suspected; barometer not used", tropo.zhd_source, bias, lim)
            warnings_flags.add("BAROMETER_REJECTED")
            baro = None
        else:
            LOG.info("Barometer vs %s pressure: median difference %.2f hPa; barometer used for ZHD (priority 1)",
                     tropo.zhd_source, bias)
            tropo = tr.apply_barometer(tropo, grid5, p5, h_orth_ap)
            if np.isnan(p5[(grid5 >= first) & (grid5 <= last)]).any():
                warnings_flags.add("BAROMETER_GAPS")
                LOG.warning("BAROMETER_GAPS: barometer gaps > %.0f min filled from %s", settings.MET_MAX_GAP_S / 60,
                            tropo.zhd_source_fallback)
    warnings_flags |= tropo.flags
    LOG.info("Troposphere a priori: ZHD source %s, mapping function %s", tropo.zhd_source, tropo.mapping_function)
    blq = None
    blq_path = os.path.join(settings.BLQ_DIR, f"{station}.blq")
    for cand in (blq_path, os.path.join(settings.BLQ_DIR, f"{station.lower()}.blq")):
        if os.path.exists(cand):
            r = read_blq(cand, station)
            if r is not None:
                blq = (r[0], r[1])
                LOG.info("Ocean loading: %s (%s)", cand, (r[2][1] if len(r[2]) > 1 else "").strip())
                blq_path = cand
                break
    if blq is None:
        LOG.warning("NO_OTL: no ocean-loading file %s; continuing without ocean tide loading. Obtain it once from "
                    "http://holt.oso.chalmers.se/loading/ (model FES2014b, centre-of-mass correction: yes, BLQ "
                    "format) and save it as blq/%s.blq", blq_path, station)
        warnings_flags.add("NO_OTL")

    # ---------------------------------------------------------------- satellite exclusion (§ 13.5)
    excl_mask = np.zeros((len(epochs), len(sel_c.sats)), dtype=bool)
    for j, prn in enumerate(sel_c.sats):
        if prn in settings.SATELLITE_BLACKLIST:
            excl_mask[:, j] = True
            LOG.info("Satellite %s excluded: configured blacklist", prn)
        elif settings.USE_BROADCAST_HEALTH and nav is not None and prn in nav:
            bad = [e for e in nav[prn] if e.health != 0]
            for e in bad:
                m_ = np.abs(epochs - e.toe) <= 2 * 3600 * ts.NS
                if m_.any():
                    excl_mask[m_, j] = True
            if bad:
                LOG.info("Satellite %s: %d epochs excluded (broadcast health flag unhealthy)", prn,
                         int(excl_mask[:, j].sum()))

    # ---------------------------------------------------------------- common preprocessing
    ctx = dict(epochs=epochs, sel=sel_c, Pif=Pif, Lif=Lif, ps=ps, atx=atx, rcv_ant=rcv_ant, tropo=tropo, blq=blq,
               rclk=rclk_spp, cutoff=np.radians(args.cutoff), a1=a1, a2=a2)
    LOG.info("Modelling observations at the a priori ARP (%s)", "SOI transformed" if coord_soi else "SPP")
    obs0, m0 = build_obsset(ctx, X_ap)
    t_s = (epochs - epochs[0]) / ts.NS
    arcs = pp.detect_slips_and_arcs(sel_c, t_s, m0.el, m0.valid & ~excl_mask, ctx["cutoff"])
    if settings.RESET_AMBIGUITIES_AT_DAY_BOUNDARY:
        bounds = [int(np.searchsorted(epochs, b)) for b in range(ts.day_start(int(epochs[0])) + ts.DAY_NS,
                                                                  int(epochs[-1]) + 1, ts.DAY_NS)]
        nb_ = pp.split_at_epochs(arcs, bounds, "DAY_BOUNDARY")
        if nb_:
            LOG.info("Ambiguities reset at %d product-day boundary crossing(s) (TDS § 11.2)", nb_)
    for (kj, _ms, repaired) in clk_jumps:
        if not repaired:                                  # unrepaired receiver clock jump -> all arcs reset (§ 5.7)
            for j in range(len(sel_c.sats)):
                a = arcs.arc[kj, j]
                if a >= 0 and arcs.arc_meta.get(int(a), {}).get("first", kj) < kj:
                    pp.split_arc(arcs, j, kj, "CLOCK_JUMP")
    obs0.usable &= arcs.arc >= 0
    LOG.info("Preprocessing: %d arcs; resets %s; excluded %s; IONO_ACTIVE epochs %d", arcs.n_arcs,
             {k: v for k, v in arcs.reasons.items() if v}, {k: v for k, v in m0.excluded.items() if v},
             int(arcs.iono_active.sum()))
    if arcs.n_arcs == 0:
        raise Stop("FAIL", "NO_USABLE_ARCS: no satellite arc survives preprocessing")

    # ---------------------------------------------------------------- estimation passes (auto mode, § 9.12)
    result = {}
    passA = None
    if mode in ("auto", "static", "fixed") or coord_soi is None:
        LOG.info("Pass A: float PPP, static_estimated coordinates")
        passA = estimate(ctx, X_ap, arcs, "static", obs0, m0)
        XA, PA = est.coordinate_estimate(passA["sol"], passA["obs"])
        sA = np.sqrt(np.diag(PA))
        LOG.info("Pass A coordinate: %.4f %.4f %.4f (sigma %.1f %.1f %.1f mm)", *XA, *(sA * 1e3))
        result["X_A"] = XA
    final_mode, final_source, final = None, None, None
    coord_final = None
    if coord_soi is None:
        final, final_mode, final_source = passA, "static_estimated", "ppp_same_day"
    elif mode == "static":
        final, final_mode, final_source = passA, "static_estimated", "ppp_same_day"
    elif mode == "constrained":
        LOG.info("Constrained mode at the SOI coordinate, sigma N,E,U = %s mm",
                 tuple(round(s * 1e3, 1) for s in settings.CONSTRAINED_SIGMA_NEU_M))
        final = estimate(ctx, coord_soi.xyz, arcs, "constrained", None, None)
        final_mode, final_source = "constrained", "soi_transformed"
    if coord_soi is not None and mode in ("auto", "fixed"):
        dec = co.auto_decision(coord_soi.xyz, coord_soi_arp.xyz, result["X_A"], hdr.delta_hen)
        dE, dN, dU = dec["d_enu"]
        LOG.info("Coordinate check SOI(transformed) - PPP same day: dN %+.1f mm, dE %+.1f mm, dU %+.1f mm "
                 "(tolerance %.0f/%.0f mm)", dN * 1e3, dE * 1e3, dU * 1e3, settings.AUTO_TOL_H_M * 1e3,
                 settings.AUTO_TOL_U_M * 1e3)
        X_fix = coord_soi.xyz
        if dec["refpoint_switched"]:
            LOG.warning("REFPOINT_SWITCHED_TO_ARP: vertical discrepancy matches DELTA H; SOI coordinate "
                        "reinterpreted as ARP")
            coord_soi = coord_soi_arp
            X_fix = coord_soi_arp.xyz
        manifest_extra["coordinate_check"] = {"dN_m": dN, "dE_m": dE, "dU_m": dU, "decision": dec["decision"],
                                              "refpoint_switched": dec["refpoint_switched"]}
        if dec["decision"] == "FIXED_SOI" or mode == "fixed":
            if dec["decision"] != "FIXED_SOI":
                LOG.warning("FORCED_FIXED: --mode fixed although the same-day check failed (see discrepancy above)")
            else:
                LOG.info("Auto mode: SOI coordinate confirmed by the data -> Pass B fixed at the SOI coordinate")
            final = estimate(ctx, X_fix, arcs, "fixed", None, None)
            final_mode, final_source = "fixed", "soi_transformed"
        else:
            dr, dn_ = co.permanent_tide_offset(np.radians(coord_soi.lat))
            LOG.warning("SOI_COORD_NOT_CONFIRMED: keeping the same-day static solution. Possible causes: velocity "
                        "(PMM, unknown vertical rate), tide system (mean-tide would differ by ~%.0f mm up), "
                        "reference point (DELTA H = %.3f m), antenna-model era of the SOI solution",
                        dr * 1e3, hdr.delta_hen[0])
            final, final_mode, final_source = passA, "static_estimated", "ppp_same_day"
    sol, obs_f = final["sol"], final["obs"]
    if final_mode == "static_estimated":
        X_used, P_used = est.coordinate_estimate(sol, obs_f)
        coord_final = co.simple_coordinate(X_used, t_mid, ps.frame_label, "ppp_same_day",
                                           sigma=np.sqrt(np.diag(P_used)), geoid_undulation=undu,
                                           history=[{"step": "static PPP same day", "passes": sol.passes}])
    else:
        coord_final = coord_soi
    coord_pull = False
    if final_mode == "constrained":
        Xc, Pc = est.coordinate_estimate(sol, obs_f)
        d = co.xyz2enu(Xc - coord_soi.xyz, np.radians(coord_soi.lat), np.radians(coord_soi.lon))
        zc = d / np.array([settings.CONSTRAINED_SIGMA_NEU_M[1], settings.CONSTRAINED_SIGMA_NEU_M[0],
                           settings.CONSTRAINED_SIGMA_NEU_M[2]])
        coord_pull = abs(zc[2]) > settings.COORD_PULL_LIMIT
        LOG.info("Constrained coordinate pull z (E,N,U) = %.2f %.2f %.2f", *zc)
        if coord_pull:
            LOG.warning("COORD_PULL: |z_U| = %.1f > %.0f", abs(zc[2]), settings.COORD_PULL_LIMIT)
        manifest_extra["coordinate_pull_enu"] = zc.tolist()
    LOG.info("Final solution: mode %s, coordinate source %s, NIS/dof %.2f, passes %d, edits %s", final_mode,
             final_source, sol.nis, sol.passes, sol.edit_counts)
    conv_t = est.convergence_time(obs_f, sol)
    if conv_t is not None:
        LOG.info("Forward-filter convergence time (|fwd - smoothed| < 10 mm for 30 min): %.0f min", conv_t / 60)

    # ---------------------------------------------------------------- 5-minute extraction (UTC day D)
    day_utc0 = ts.to_ns(*ts.from_ns(day)[:3])            # UTC calendar day with the same date
    epoch_flags = {"ECLIPSE_EXCLUSION": final["model"].eclipse_epochs,
                   "ORBIT_EDGE": final["model"].orbit_edge_epochs if "ORBIT_EDGE" in ps.flags else None,
                   "CLOCK_INTERP": final["model"].clock_interp_epochs}
    gflags = [f for f in warnings_flags if f in settings.QC_BITS and f not in ("ORBIT_EDGE",)]
    if ps.tier != "FINAL":
        gflags.append("RAPID_TIER")
    fm = est.extract_5min(obs_f, sol, day_utc0, day_utc0 + ts.DAY_NS, gflags, arcs.iono_active, epoch_flags,
                          tuple(edge_ok), coord_pull)

    # ---------------------------------------------------------------- PPP-AR layer (separate, never overwrites)
    ar = None
    if args.ar:
        ar = run_ar_layer(ps, sel_c, final, final_mode, fm, arcs, day_utc0, gflags, epoch_flags, edge_ok,
                          coord_pull, t_mid)
        manifest_extra["ppp_ar"] = {"status": ar["status"], "metrics": ar.get("metrics", {}),
                                    "reason": ar.get("reason")}

    # ---------------------------------------------------------------- manifest + outputs
    chash, csnap = prd.config_hash()
    man = prd.new_manifest(station, f"{y:04d}-{doy:03d}")
    files_used = [f["filename"] for f in ps.files]
    meta = {"station": station, "lat": coord_final.lat, "lon": coord_final.lon, "h": coord_final.h_ell,
            "product_AC": "COD", "product_tier": ps.tier, "product_files": _json(files_used),
            "coordinate_mode": final_mode, "coordinate_source": final_source,
            "mapping_function": tropo.mapping_function, "ZWD_process_noise": settings.SIGMA_ZWD_M_SQRT_H * 1e3,
            "software_version": man["software_version"], "solution_type": "FLOAT", "zhd_source": tropo.zhd_source,
            "manifest_id": man["manifest_id"], "config_hash": chash, "elevation_cutoff": args.cutoff,
            "AR_status": "AR_NOT_AVAILABLE" if not args.ar else "FLOAT", "frame": ps.frame_label}
    prev_mid = _archive_previous_outputs(args.out, tag)  # outputs are immutable (TDS § 19.3)
    rows = out.ztd_rows(fm, meta)
    ztd_csv = os.path.join(args.out, tag + "_ZTD.csv")
    out.write_csv(ztd_csv, rows, out.ZTD_COLUMNS, f"{settings.SOFTWARE_NAME} station ZTD product, float PPP, "
                  f"CODE {ps.family}")
    written = [ztd_csv]
    pw_res = None
    if not args.no_pwv:
        pw_res, pw_meta = compute_pwv(fm, tropo, coord_final, final_mode, ant_status, baro, era5_for_day(day))
        prow = []
        for r, i in zip(rows, range(len(rows))):
            rr = dict(r)
            rr.update({"PWV": pw_res.pwv[i], "sigma_PWV_ppp": pw_res.sig_ppp[i], "sigma_PWV_conv": pw_res.sig_conv[i],
                       "sigma_PWV_total": pw_res.sig_total[i], "P_ant": float(pw_res.p_ant[i]),
                       "P_source": pw_res.p_source[i], "T_m": float(pw_res.tm[i]),
                       "Tm_source": pw_res.tm_source[i],
                       "Pi": float(pw_res.pi[i]), "ZHD_met": pw_res.zhd_met[i], "ZWD_met": pw_res.zwd_met[i],
                       "H_ant": coord_final.h_orth, "constants_set": pw_res.constants_set})
            prow.append(rr)
        pwv_csv = os.path.join(args.out, tag + "_PWV.csv")
        out.write_csv(pwv_csv, prow, out.ZTD_COLUMNS + out.PWV_EXTRA,
                      f"{settings.SOFTWARE_NAME} station PWV product ({pw_meta})")
        written.append(pwv_csv)
        rows_nc = prow
        cols_nc = out.ZTD_COLUMNS + out.PWV_EXTRA
    else:
        rows_nc, cols_nc = rows, out.ZTD_COLUMNS
    if ar is not None and ar["status"] in ("FIXED", "PARTIAL"):
        meta_fx = dict(meta, solution_type="FIXED", AR_status=ar["status"])
        rows_fx = out.ztd_rows(ar["fm"], meta_fx)
        for r, nfx in zip(rows_fx, ar["n_fixed"]):
            r["n_amb_fixed"] = int(nfx)
        fx_csv = os.path.join(args.out, tag + "_FIXED_ZTD.csv")
        out.write_csv(fx_csv, rows_fx, out.ZTD_COLUMNS, f"{settings.SOFTWARE_NAME} station ZTD, PPP-AR FIXED layer "
                      f"(experimental until the Stage-6 gate is passed); the FLOAT product is {os.path.basename(ztd_csv)}")
        written.append(fx_csv)
    nc_path = os.path.join(args.out, tag + "_ZTD.nc")
    try:
        out.write_netcdf(nc_path, rows_nc, cols_nc, {
            "title": f"GNSS PPP ZTD/PWV station product {station} {y:04d}-{doy:03d}",
            "summary": "Float PPP (GPS, ionosphere-free) with CODE products; RTS-smoothed 5-min ZTD",
            "institution": "", "source": "GPS observations", "station": station, "manifest_id": man["manifest_id"],
            "config_hash": chash, "software_version": man["software_version"], "product_tier": ps.tier,
            "product_family": ps.family, "solution_type": "FLOAT", "coordinate_mode": final_mode,
            "coordinate_source": final_source, "frame_label": ps.frame_label, "mapping_function": tropo.mapping_function,
            "zhd_source": tropo.zhd_source, "antex": atx.name, "date_created": man["run_timestamp"]},
            coord_final.lat, coord_final.lon, coord_final.h_ell)
        written.append(nc_path)
    except Exception as exc:                  # noqa: BLE001 - NetCDF is a secondary format; never lose the CSV
        LOG.warning("NETCDF_WRITE_FAILED: %s", exc)
    tro = os.path.join(args.out, tag + ".tro")
    out.write_sinex_tro(tro, station, fm, meta, coord_final.xyz)
    written.append(tro)
    ql = os.path.join(args.out, tag + "_quicklook.png")
    try:
        out.quicklook(ql, station, f"{y:04d}-{doy:03d}", fm, pw_res,
                      f"| {final_mode} | CODE {ps.family} | {tropo.mapping_function}")
        written.append(ql)
    except Exception as exc:                  # noqa: BLE001 - plotting problems must not end the run
        LOG.warning("QUICKLOOK_FAILED: %s", exc)
    warns = plog.warnings()
    man.update({
        "rinex_files": [{"name": os.path.basename(p.path), "sha256": prd.sha256(p.path),
                         "header": {"marker": p.header.marker_name, "receiver": p.header.receiver,
                                    "antenna": p.header.antenna_type, "delta_hen": list(p.header.delta_hen)}}
                        for p in parts],
        "antenna_type": hdr.antenna_type, "delta_hen": list(hdr.delta_hen),
        "family": ps.family, "tier": ps.tier, "version_char": settings.FAMILY_VERSION_CHAR.get(ps.family, ""),
        "products": ps.files, "frame_label": ps.frame_label,
        "antex": {"name": atx.name, "sha256": (atx_meta or {}).get("sha256"), "status": atx_status,
                  "receiver_calibration": ant_status, "clock_header_antex": ps.antex_name},
        "blq": {"file": blq_path if blq is not None else None, "model": "user-supplied BLQ (FES2014b+CMC requested)",
                "sha256": prd.sha256(blq_path) if blq is not None else None},
        "vmf3_files": tropo.files if "VMF3_GRID" in (tropo.zhd_source, tropo.zhd_source_fallback) else [],
        "troposphere": {"zhd_source": tropo.zhd_source, "mapping_function": tropo.mapping_function,
                        "gpt3": tropo.files if "GPT3" in (tropo.zhd_source, tropo.zhd_source_fallback) else [],
                        "zhd_source_in_gaps": tropo.zhd_source_fallback},
        "leap_seconds": leap, "processing_mode": "FLOAT", "status": ps.status, "reason": ps.reason,
        "fallbacks": ps.fallbacks, "missing": ps.missing, "flags": sorted(warnings_flags),
        "coordinate_object_id": coord_final.id, "coordinate_object": coord_final.as_dict(),
        "coordinate_soi": coord_soi.as_dict() if coord_soi is not None else None,
        "coordinate_mode": final_mode, "coordinate_source": final_source,
        "config_hash": chash, "window": [ts.iso(win0), ts.iso(win1)], "edge_overlap": edge_ok,
        "preprocessing": {"arcs": arcs.n_arcs, "resets": arcs.reasons, "excluded": m0.excluded,
                          "clock_jumps": clk_jumps, "osb": {"applied": osb_info["applied"],
                                                             "missing": osb_info["missing"][:50]}},
        "estimation": {"nis_per_dof": sol.nis, "passes": sol.passes, "edits": sol.edit_counts,
                       "postfit_rms": sol.postfit_rms, "segments": len(sol.segments),
                       "cholesky_failures": sol.chol_failures, "convergence_time_s": conv_t,
                       "sigma_zwd_m_sqrt_h": settings.SIGMA_ZWD_M_SQRT_H, "attitude": final["model"].attitude_source},
        "warnings": warns, "outputs": [os.path.basename(p) for p in written],
        "statement": out.STATEMENT, "download_attempts": dl.attempts[-500:], **manifest_extra,
    })
    if prev_mid:
        man["supersedes"] = [prev_mid]
        old_m = os.path.join(args.out, "superseded", prev_mid, tag + "_manifest.json")
        try:
            import json as _j
            with open(old_m) as fh:
                om = _j.load(fh)
            om["superseded_by"] = man["manifest_id"]
            prd.write_json(old_m, om)
        except (OSError, ValueError):
            pass
    cfg_dir = os.path.join(args.out, "provenance", "configs")
    os.makedirs(cfg_dir, exist_ok=True)
    prd.write_json(os.path.join(cfg_dir, chash + ".json"), csnap)
    mpath = os.path.join(args.out, tag + "_manifest.json")
    prd.write_json(mpath, man)
    written.append(mpath)
    for p in written:
        LOG.info("  wrote %s", p)
    # ---------------------------------------------------------------- summary
    n_val = int(np.isfinite(fm.zwd).sum())
    conv_frac = np.mean([c == "CONVERGED" for c in fm.conv]) if len(fm.conv) else 0.0
    ztd_mm = (fm.zhd0 + fm.zwd) * 1e3
    LOG.info("SUMMARY %s %04d-%03d: %d/%d 5-min values, %.0f %% converged, mean ZTD %.1f mm%s, flags %s, "
             "run time %.0f s", station, y, doy, n_val, len(fm.zwd), 100 * conv_frac, np.nanmean(ztd_mm),
             (f", mean PWV {np.nanmean(pw_res.pwv):.1f} mm" if pw_res is not None else ""),
             ",".join(sorted(warnings_flags)) or "none", time.time() - t_start)
    return {"fm": fm, "sol": sol, "obs": obs_f, "pwv": pw_res, "manifest": man, "coord": coord_final,
            "mode": final_mode, "files": written, "passA": passA}


def _archive_previous_outputs(out_dir, tag):
    """Move outputs of a previous run of the same station-day to superseded/<manifest_id>/ (never overwrite)."""
    import json
    import shutil
    mpath = os.path.join(out_dir, tag + "_manifest.json")
    if not os.path.exists(mpath):
        return None
    try:
        with open(mpath) as fh:
            mid = json.load(fh).get("manifest_id", "unknown")
    except (OSError, ValueError):
        mid = "unknown"
    dest = os.path.join(out_dir, "superseded", mid)
    os.makedirs(dest, exist_ok=True)
    for f in glob.glob(os.path.join(out_dir, tag + "*")):
        if os.path.isfile(f) and not f.endswith(".log"):
            shutil.move(f, os.path.join(dest, os.path.basename(f)))
    LOG.info("Previous outputs of %s moved to %s (superseded)", tag, dest)
    return mid


def _json(v):
    import json
    return json.dumps(v)


def _previous_manifest(out_dir, station, day):
    import json
    best = None
    for f in sorted(glob.glob(os.path.join(out_dir, f"{station}_*_manifest.json"))):
        try:
            with open(f) as fh:
                m = json.load(fh)
        except (OSError, ValueError):
            continue
        best = m
    return best


def build_obsset(ctx, X0):
    """Compute the observation model at X0 and assemble the estimator input (TDS § 4.2, § 5.4)."""
    epochs = ctx["epochs"]
    ps = ctx["ps"]
    disp, _parts = mdl.station_displacements(X0, epochs, ps.erp, ctx["blq"])
    tro = ctx["tropo"]
    zhd0, zwd0, _, _ = tro.at(epochs)
    for it in range(3):
        m = mdl.compute(epochs, ctx["sel"].sats, X0, ctx["rclk"], ps, ctx["atx"], ctx["rcv_ant"], disp)
        ne, ns = m.el.shape
        mh = np.full((ne, ns), np.nan)
        mw = np.full((ne, ns), np.nan)
        mg = np.full((ne, ns), np.nan)
        ok = m.valid & np.isfinite(m.el) & (m.el > np.radians(1.0))
        if ok.any():
            tt = np.broadcast_to(epochs[:, None], (ne, ns))[ok]
            h_, w_, g_ = tro.mapping(tt, m.el[ok])
            mh[ok], mw[ok], mg[ok] = h_, w_, g_
        # per-epoch code-derived receiver clock for the reception time (TDS § 7.1): refine until < 1 us
        use = ok & np.isfinite(ctx["Pif"]) & (np.nan_to_num(m.el, nan=-1) >= ctx["cutoff"])
        slant = mh * zhd0[:, None] + mw * zwd0[:, None]
        cclk = mdl.receiver_clock_from_code(ctx["Pif"], m.rho, m.dts, slant, use) / settings.C_LIGHT
        good = np.isfinite(cclk)
        new = np.where(good, cclk, ctx["rclk"])
        change = np.nanmax(np.abs(new - ctx["rclk"])) if good.any() else 0.0
        ctx["rclk"] = new
        if change < 1e-6:
            break
        LOG.debug("receiver clock refined (max change %.3g s); re-modelling", change)
    lam_nl = C_NL()
    usable = m.valid & np.isfinite(ctx["Pif"]) & np.isfinite(ctx["Lif"]) & np.isfinite(mh) & \
        (np.nan_to_num(m.el, nan=-1) >= ctx["cutoff"])
    f_if = np.sqrt(ctx["a1"] ** 2 + ctx["a2"] ** 2)
    obs = est.ObsSet(epochs=epochs, t_s=(epochs - epochs[0]) / ts.NS, sats=ctx["sel"].sats, Pif=ctx["Pif"],
                     Lif=ctx["Lif"], rho=m.rho, dts=m.dts, wind=lam_nl * m.windup, el=m.el, az=m.az, los=m.los, mh=mh,
                     mw=mw, mg=mg, zhd0=zhd0, zwd0=zwd0, usable=usable, f_if=f_if, X0=np.asarray(X0, dtype=float))
    return obs, m


def C_NL():
    return settings.C_LIGHT / (settings.FREQ[("G", "1")] + settings.FREQ[("G", "2")])


def estimate(ctx, X0, arcs0, mode, obs=None, m=None):
    """One estimation pass (fixed | static | constrained) incl. static re-linearisation (TDS § 5.6)."""
    X0 = np.asarray(X0, dtype=float)
    if obs is None:
        obs, m = build_obsset(ctx, X0)
    arcs = copy.deepcopy(arcs0)
    obs.usable &= arcs.arc >= 0
    cfg = est.EstConfig(mode=mode)
    sol = est.solve(obs, arcs, cfg)
    if mode == "static":
        for _ in range(settings.STATIC_RELINEARISE_PASSES):
            X1, _P = est.coordinate_estimate(sol, obs)
            if X1 is None:
                break
            LOG.info("  re-linearising at the converged coordinate (shift %.3f m)", np.linalg.norm(X1 - X0))
            X0 = X1
            obs, m = build_obsset(ctx, X0)
            arcs = copy.deepcopy(arcs0)
            obs.usable &= arcs.arc >= 0
            sol = est.solve(obs, arcs, cfg)
    return {"sol": sol, "obs": obs, "model": m, "X0": X0, "arcs": arcs}


def run_ar_layer(ps, sel, final, final_mode, fm_float, arcs, day_utc0, gflags, epoch_flags, edge_ok, coord_pull,
                 t_mid):
    """PPP-AR (TDS § 15): only with product status OK_AR; result is a separate FIXED layer."""
    if ps.status != "OK_AR":
        LOG.warning("AR_NOT_AVAILABLE: product status %s (%s); float solution only", ps.status, ps.reason)
        return {"status": "AR_NOT_AVAILABLE", "reason": ps.reason}
    P1c, P2c, L1c, L2c, info = pp.apply_osb(sel, ps.osb, t_mid, phase=True)
    sel_ar = copy.copy(sel)
    sel_ar.P1, sel_ar.P2, sel_ar.L1, sel_ar.L2 = P1c, P2c, L1c, L2c
    Pif, Lif, _, _ = pp.if_combination(sel_ar)
    obs_ar = copy.copy(final["obs"])
    obs_ar.Pif, obs_ar.Lif = Pif, Lif
    obs_ar.usable = final["obs"].usable & np.isfinite(Lif) & np.isfinite(Pif)
    mode = {"static_estimated": "static"}.get(final_mode, final_mode)
    extract = lambda s_: est.extract_5min(obs_ar, s_, day_utc0, day_utc0 + ts.DAY_NS, gflags, arcs.iono_active,  # noqa: E731
                                          epoch_flags, tuple(edge_ok), coord_pull)
    LOG.info("PPP-AR: wide-lane / narrow-lane resolution (phase OSB applied to %d signals)", info["applied"])
    res = amb.run_ar(obs_ar, arcs, final["sol"], sel_ar, mode, fm_float, extract)
    LOG.info("PPP-AR status %s: %s", res.status, {k: (round(v, 3) if isinstance(v, float) else v)
                                                  for k, v in res.metrics.items()})
    if res.status == "AR_REJECTED_CONSISTENCY":
        LOG.warning("AR_REJECTED_CONSISTENCY: fixed ZTD differs from float beyond the gate; fixed layer not published")
    out_ = {"status": res.status, "metrics": res.metrics}
    if res.status in ("FIXED", "PARTIAL"):
        fixed_arcs = {a for c in res.constraints for a in c[:2]}
        n_fixed = []
        for tg in res.fm.t_gpst:
            kn = int(np.clip(np.searchsorted(obs_ar.epochs, tg), 0, len(obs_ar.epochs) - 1))
            n_fixed.append(len(fixed_arcs & {int(a) for a in res.solution.arcs.arc[kn] if a >= 0}))
        out_.update({"fm": res.fm, "n_fixed": n_fixed})
    return out_


def compute_pwv(fm, tropo, coord, mode, ant_status, baro=None, bg=None):
    """Station PWV (TDS § 18) with the § 18.2 source priority, chosen per 5-min epoch.

    P_ant: station barometer (reduced to the final ARP height) > ERA5 profile > a priori ZHD source (VMF3 grid,
    GPT3 climatology).  Tm: ERA5 profile integration > GPT3 > Bevis Tm-Ts with the station/GPT3 temperature."""
    lat = np.radians(coord.lat)
    h_orth = coord.h_orth
    n = len(fm.zwd)
    grid_src = tropo.zhd_source_fallback if tropo.zhd_source == "BAROMETER" else tropo.zhd_source
    grid_src = {"VMF3_GRID": "VMF3_GRID", "GPT3": "GPT3"}.get(grid_src, grid_src)
    p = tr.saastamoinen_pressure(fm.zhd0, lat, h_orth)
    p_src = np.full(n, grid_src, dtype=object)
    tm = np.full(n, np.nan)
    tm_src = np.full(n, "", dtype=object)
    if bg is not None:
        p_e, tm_e = pw.era5_station_met(bg, coord.lat, coord.lon, h_orth, fm.t_utc)
        okp = np.isfinite(p_e)
        p[okp], p_src[okp] = p_e[okp], "ERA5"
        okt = np.isfinite(tm_e)
        tm[okt], tm_src[okt] = tm_e[okt], "ERA5"
        if not okp.all():
            LOG.warning("ERA5_PARTIAL: ERA5 covers %d of %d epochs", int(okp.sum()), n)
    if baro is not None:
        p_b = pw.pressure_at(baro, fm.t_gpst)
        if abs(h_orth - baro.h_ant) > 0.01:                # a priori -> final antenna height
            tv = np.interp(fm.t_gpst.astype(float), baro.t.astype(float), baro.t_surface_k)
            p_b = pw.reduce_pressure(p_b, baro.h_ant, h_orth, tv)
        okb = np.isfinite(p_b)
        p[okb], p_src[okb] = p_b[okb], "BAROMETER"
    need = ~np.isfinite(tm)
    if need.any():
        if tropo.met.get("Tm"):
            tm[need], tm_src[need] = tropo.met["Tm"], "GPT3"
        else:
            if baro is not None and baro.t_surface_k is not None:
                ts_k = np.interp(fm.t_gpst.astype(float), baro.t.astype(float), baro.t_surface_k)
            elif tropo.met.get("T") is not None:
                ts_k = np.full(n, tropo.met["T"] + 273.15)
            else:
                ts_k = np.full(n, 288.15)
                LOG.warning("TM_DEFAULT: no ERA5/GPT3/station temperature; Bevis Tm with standard temperature")
            tm[need], tm_src[need] = 70.2 + 0.72 * ts_k[need], "BEVIS"
    srcs = {str(k): int(v) for k, v in zip(*np.unique(p_src.astype(str), return_counts=True))}
    tms = {str(k): int(v) for k, v in zip(*np.unique(tm_src.astype(str), return_counts=True))}
    LOG.info("PWV conversion: P %s (at ARP, H = %.1f m), Tm %s (median %.1f K)", srcs, h_orth, tms, np.median(tm))
    if np.any(p_src == "GPT3"):
        LOG.warning("P_CLIMATOLOGY: PWV uses GPT3 climatological pressure for %d epochs (not acceptable for research "
                    "PWV)", int(np.sum(p_src == "GPT3")))
    s_coord = settings.HEIGHT_ZTD_COUPLING_BETA * np.sqrt(coord.covariance[2][2]) if mode != "static_estimated" else 0.0
    if ant_status == "RADOME_FALLBACK":
        s_coord = np.hypot(s_coord, settings.RADOME_FALLBACK_ZTD_SIGMA_M)
    res = pw.convert(fm.zhd0 + fm.zwd, fm.sig, lat, h_orth, p, p_src, tm, tm_src, s_coord)
    return res, f"P {'+'.join(srcs)}, Tm {'+'.join(tms)}, constants {res.constants_set}"


def main(argv=None):
    args = parse_args(argv)
    try:
        run(args)
        return 0
    except Stop as s:
        LOG.error("%s", s.message)
        LOG.error("Run ended with status %s (no ZTD written).", s.status)
        return s.code
    except KeyboardInterrupt:
        LOG.error("Interrupted by user")
        return 4
    except Exception as exc:                  # noqa: BLE001 - never a traceback on screen (TDS § 0.5)
        LOG.error("UNEXPECTED_ERROR: %s: %s. Details are in the log file.", type(exc).__name__, exc)
        for h in LOG.handlers:
            if hasattr(h, "baseFilename"):
                h.stream.write(traceback.format_exc())
        return 4


if __name__ == "__main__":
    sys.exit(main())
