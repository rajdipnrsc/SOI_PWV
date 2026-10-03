#!/usr/bin/env python3
"""Stage 9 reprocessing script (TDS § 19.3): find station-days to re-run and re-run them in parallel processes.

    python network/reprocess.py RESULTS_DIR --rinex-dir DIR [--jobs N] [--dry-run]
    python network/reprocess.py RESULTS_DIR --rinex-dir DIR --campaign 2025A [--settings my.json] [--jobs N]
    python network/reprocess.py RESULTS_DIR --maps MAPS_DIR           (which national grids must be rebuilt)

Triggers (§ 19.3): (i) the best available product of a station-day is not FINAL and the CODE final products for
that day are now verified in the cache/archive; (ii) a declared reprocessing campaign (every station-day whose best
product was not produced in that campaign).  Every re-run is an ordinary `pwv_ppp.py SITE.o SITE.n` call; outputs
are immutable (the previous ones move to superseded/<manifest_id>/), selection is by best_available / as_of.
Grids: a day's grid is rebuilt when the share of FINAL station inputs reaches 100 % (or on a campaign).
"""
import argparse
import concurrent.futures as cf
import glob
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)

from ppp import products as prd  # noqa: E402
from ppp import settings  # noqa: E402
from ppp import timesys as ts  # noqa: E402


def day_ns(day_str):
    y, d = (int(x) for x in day_str.split("-"))
    return ts.to_ns(y, 1, 1) + (d - 1) * ts.DAY_NS


def final_available(dl, d_ns):
    """True if COD0OPSFIN orbit and clock for the day are verified (cache or download)."""
    for ptype in ("SP3", "CLK"):
        p, _ = prd.fetch_product(dl, ptype, "COD0OPSFIN", d_ns)
        if p is None:
            return False
    return True


def locate_inputs(m, rinex_dir):
    """Observation and navigation file of a manifest in rinex_dir (names recorded in the manifest)."""
    obs = (m.get("rinex_files") or [{}])[0].get("name")
    nav = (m.get("nav_file") or {}).get("name")
    hits = {}
    for key, name in (("obs", obs), ("nav", nav)):
        if not name:
            continue
        g = glob.glob(os.path.join(rinex_dir, "**", name), recursive=True)
        if g:
            hits[key] = g[0]
    if "obs" in hits and "nav" not in hits:                 # older manifests: pair by name (SITEdddX.yyo -> .yyn)
        b = hits["obs"]
        for cand in (b[:-1] + "n", b[:-1] + "N", b.replace("_MO.", "_GN."), b.replace("_MO.", "_MN.")):
            if cand != b and os.path.exists(cand):
                hits["nav"] = cand
                break
    return hits.get("obs"), hits.get("nav")


def run_one(cmd, env):
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    return cmd, r.returncode, (r.stderr or "").strip().splitlines()[-1:] or [""]


def maps_to_rebuild(results_dir, maps_dir):
    """Days whose grid was built with < 100 % FINAL station inputs while all inputs are now FINAL."""
    best = prd.select_manifests(prd.scan_manifests(results_dir))
    by_day = {}
    for (_st, day), m in best.items():
        by_day.setdefault(day.replace("-", ""), []).append(m.get("tier"))
    out = []
    for f in sorted(glob.glob(os.path.join(maps_dir, "INPWV_*_*_*_*_v*.nc"))):
        parts = os.path.basename(f).split("_")
        tag, tier = parts[3], parts[4]
        tiers = by_day.get(tag, [])
        if tier != "FINAL" and tiers and all(t == "FINAL" for t in tiers):
            out.append((tag, f))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="Re-run station-days for Rapid->Final or a reprocessing campaign.")
    ap.add_argument("results")
    ap.add_argument("--rinex-dir", help="root folder of the RINEX archive (searched recursively)")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--campaign", default=None, help="declare a reprocessing campaign (label stored in manifests)")
    ap.add_argument("--settings", default=None, help="JSON settings override for the campaign")
    ap.add_argument("--maps", default=None, help="maps folder: list grids that must be rebuilt")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--offline", action="store_true", help="check final products in the cache only")
    a = ap.parse_args(argv)

    if a.maps:
        todo = maps_to_rebuild(a.results, a.maps)
        for tag, f in todo:
            print(f"REBUILD {tag}: {os.path.basename(f)} -> python pwv_map.py {a.results} --day {tag[:4]}-{tag[4:]}")
        print(f"{len(todo)} grid(s) to rebuild")
        if not a.rinex_dir:
            return 0
    if not a.rinex_dir:
        ap.error("--rinex-dir is required for station reprocessing")

    best = prd.select_manifests(prd.scan_manifests(a.results))
    env = dict(os.environ)
    todo, notes = [], []
    if a.campaign:
        over = {}
        if a.settings:
            with open(a.settings) as fh:
                over = json.load(fh)
        over["PROCESSING_CAMPAIGN"] = a.campaign
        cpath = os.path.join(a.results, "provenance", f"campaign_{a.campaign}.json")
        os.makedirs(os.path.dirname(cpath), exist_ok=True)
        prd.write_json(cpath, over)
        env["PWV_PPP_SETTINGS"] = os.path.abspath(cpath)
        todo = [m for m in best.values() if m.get("campaign") != a.campaign]
        notes.append(f"campaign {a.campaign}: {len(todo)} of {len(best)} station-days not yet processed in it")
    else:
        dl = prd.Downloader(settings.CACHE_DIR, offline=a.offline)
        cache = {}
        now = ts.utc_to_gpst(prd._utc_now_ns())
        lat_ns = int(settings.MAX_EXPECTED_LATENCY_H["COD0OPSFIN"] * 3600) * ts.NS
        for m in best.values():
            if m.get("tier") == "FINAL":
                continue
            d = day_ns(m["day"])
            if d not in cache:
                cache[d] = final_available(dl, d)
            if cache[d]:
                todo.append(m)
            elif now - d > lat_ns:
                notes.append(f"{m['station']} {m['day']}: final products overdue and not available")
        notes.append(f"Rapid->Final: {len(todo)} station-day(s) can be upgraded")
    for n in notes:
        print(n)
    cmds = []
    for m in todo:
        obs, nav = locate_inputs(m, a.rinex_dir)
        if not obs or not nav:
            print(f"SKIP {m['station']} {m['day']}: RINEX/nav file not found under {a.rinex_dir}")
            continue
        cmds.append([sys.executable, os.path.join(HERE, "pwv_ppp.py"), obs, nav, "--out", a.results, "--quiet"]
                    + (["--offline"] if a.offline else []))
    if a.dry_run:
        for c in cmds:
            print("WOULD RUN", " ".join(c[1:]))
        return 0
    bad = 0
    with cf.ThreadPoolExecutor(max_workers=max(1, a.jobs)) as ex:
        for cmd, code, last in ex.map(lambda c: run_one(c, env), cmds):
            print(f"{'OK  ' if code == 0 else 'EXIT %d' % code} {os.path.basename(cmd[2])} {last[0] if code else ''}")
            bad += code != 0
    print(f"{len(cmds) - bad} of {len(cmds)} re-runs succeeded")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
