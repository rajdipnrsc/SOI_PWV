#!/usr/bin/env python3
"""Stage 10 monitoring (TDS § 35): summary of a results folder for one day or a range of days.

    python network/monitor.py RESULTS_DIR [--day YYYY-DDD] [--days N] [--expect stations.csv]

Reads the station-day manifests and logs (no database): status per station, tier, converged fraction, flags and
WARNING/ERROR counts; lists expected stations without a product and runs that ended with an ERROR.  Exit code 1 if
any expected station-day is missing or failed, so a cron wrapper can alert.
"""
import argparse
import collections
import csv
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)

from ppp import products as prd  # noqa: E402
from ppp import timesys as ts  # noqa: E402


def expected_stations(path):
    if not path or not os.path.exists(path):
        return set()
    with open(path) as fh:
        rd = csv.DictReader(fh)
        key = next((k for k in rd.fieldnames or [] if k.strip().lower() in ("station", "site", "name", "id")), None)
        return {r[key].strip()[:4].upper() for r in rd if key and r.get(key)}


def summarise(results, day):
    tag = day.replace("-", "")
    best = {k[0]: m for k, m in prd.select_manifests(prd.scan_manifests(results)).items() if k[1] == day}
    rows, flags = [], collections.Counter()
    for st, m in sorted(best.items()):
        f = prd.selected_file(m, "_ZTD.csv")
        conv = ""
        if f:
            with open(f) as fh:
                data = [r for r in csv.DictReader(L for L in fh if not L.startswith("#"))]
            if data:
                conv = f"{100 * sum(r['convergence_flag'] == 'CONVERGED' for r in data) / len(data):.0f}%"
        flags.update(m.get("flags", []))
        rows.append((st, m.get("status"), m.get("tier"), conv, ",".join(m.get("flags", []))))
    errors = {}
    for lg in glob.glob(os.path.join(results, f"*_{tag}.log")):
        st = os.path.basename(lg)[:4]
        with open(lg, errors="replace") as fh:
            txt = fh.read()
        last_run = txt.split(" pwv_ppp ")[-1] if " pwv_ppp " in txt else txt
        n_w = len(re.findall(r" WARNING ", last_run))
        e = re.findall(r" ERROR\s+(.*)", last_run)
        if e:
            errors[st] = (n_w, e[-1])
    return rows, flags, errors


def main(argv=None):
    ap = argparse.ArgumentParser(description="Daily processing monitor (manifests + logs).")
    ap.add_argument("results")
    ap.add_argument("--day", default=None, help="YYYY-DDD (default: yesterday UTC)")
    ap.add_argument("--days", type=int, default=1, help="number of days ending at --day")
    ap.add_argument("--expect", default=None, help="CSV with the expected station list (column station)")
    a = ap.parse_args(argv)
    if a.day:
        y, d = (int(x) for x in a.day.split("-"))
        end = ts.to_ns(y, 1, 1) + (d - 1) * ts.DAY_NS
    else:
        end = ts.day_start(prd._utc_now_ns()) - ts.DAY_NS
    exp = expected_stations(a.expect)
    problem = False
    for k in range(a.days - 1, -1, -1):
        y, d = ts.year_doy(end - k * ts.DAY_NS)
        day = f"{y:04d}-{d:03d}"
        rows, flags, errors = summarise(a.results, day)
        tiers = collections.Counter(r[2] for r in rows)
        print(f"== {day}: {len(rows)} station products ({dict(tiers)}), {len(errors)} run(s) ended with ERROR")
        for r in rows:
            print(f"   {r[0]}  {r[1]:<14} {r[2] or '':<8} converged {r[3]:>4}  {r[4]}")
        for st, (n_w, e) in sorted(errors.items()):
            if st not in {r[0] for r in rows}:
                print(f"   {st}  ERROR: {e} ({n_w} warnings)")
        if flags:
            print("   flags: " + ", ".join(f"{k} x{v}" for k, v in flags.most_common()))
        missing = sorted(exp - {r[0] for r in rows})
        if missing:
            print(f"   MISSING ({len(missing)}): {' '.join(missing)}")
        problem |= bool(missing) or any(st not in {r[0] for r in rows} for st in errors)
    return 1 if problem else 0


if __name__ == "__main__":
    sys.exit(main())
