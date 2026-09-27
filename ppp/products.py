"""CODE product search, download, cache, verification, status codes and provenance manifest.

TDS § 11 (CODE-only policy, status codes), § 12 (date logic, fallback, download mechanics, verification,
cache layout, manifest), § 33.6 (automatic downloads of ANTEX, VMF3, GPT3, leap seconds).

Nothing here hard-codes a directory tree or filename: all templates live in settings.py (TDS § 36.5).
"""
import datetime as _dt
import gzip
import hashlib
import json
import os
import random
import shutil
import time
import uuid
from dataclasses import dataclass, field

import numpy as np

from . import log as plog
from . import orbclk
from . import settings
from . import timesys as ts

LOG = plog.get()


# =============================================================================================== HTTP
class Downloader:
    """HTTP(S) downloads with proxy, CA bundle, .netrc, retries/backoff, atomic writes (TDS § 12.3)."""

    def __init__(self, cache_dir, offline=False, proxy=None):
        self.cache_dir = cache_dir
        self.offline = offline
        self.proxy = proxy
        self.attempts = []                 # log of every attempt (URL, status, bytes, duration)
        self.dead_hosts = {}               # host -> reason (circuit breaker for unreachable hosts)
        self._session = None

    def session(self):
        if self._session is None:
            import requests
            s = requests.Session()
            s.trust_env = True             # HTTPS_PROXY / HTTP_PROXY / ~/.netrc
            if self.proxy:
                s.proxies = {"http": self.proxy, "https": self.proxy}
            if settings.CA_BUNDLE:
                s.verify = settings.CA_BUNDLE
            s.headers["User-Agent"] = f"{settings.SOFTWARE_NAME}/{settings.SOFTWARE_VERSION}"
            self._session = s
        return self._session

    def get(self, url, dest):
        """Download url to dest atomically. Returns True on success, False if not available."""
        if self.offline:
            return False
        import requests
        host = url.split("/")[2]
        if host in self.dead_hosts:
            return False
        staging = os.path.join(self.cache_dir, "staging")
        os.makedirs(staging, exist_ok=True)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        part = os.path.join(staging, os.path.basename(dest) + f".{uuid.uuid4().hex[:8]}.part")
        for attempt in range(settings.HTTP_RETRIES):
            t0 = time.time()
            status, nbytes, err = None, 0, ""
            try:
                with self.session().get(url, stream=True, allow_redirects=True,
                                        timeout=(settings.HTTP_CONNECT_TIMEOUT_S, settings.HTTP_READ_TIMEOUT_S)) as r:
                    status = r.status_code
                    if status == 200:
                        ctype = r.headers.get("Content-Type", "")
                        with open(part, "wb") as fh:
                            for chunk in r.iter_content(1 << 16):
                                fh.write(chunk)
                                nbytes += len(chunk)
                        if "text/html" in ctype and nbytes < 200000 and url.endswith((".gz", ".Z")):
                            # login page / error page instead of the product (e.g. Earthdata redirect)
                            status = "HTML"
                            os.remove(part)
                        else:
                            os.replace(part, dest)
                            self._log(url, status, nbytes, time.time() - t0, "")
                            return True
            except requests.exceptions.RequestException as exc:
                err = type(exc).__name__ + ": " + str(exc)[:160]
            self._log(url, status, nbytes, time.time() - t0, err)
            if os.path.exists(part):
                os.remove(part)
            if status in (404, 410, "HTML", 401):
                return False
            if status == 403 or (err and ("ProxyError" in err or "403" in err)):
                self.dead_hosts[host] = err or "HTTP 403"
                LOG.warning("DOWNLOAD_BLOCKED: host %s refused (%s); it is skipped for the rest of this run",
                            host, err or "HTTP 403")
                return False
            if attempt + 1 < settings.HTTP_RETRIES:
                wait = min(2 ** (attempt + 1) + random.uniform(0, 1), settings.HTTP_BACKOFF_MAX_S)
                LOG.debug("retry %d for %s in %.1f s", attempt + 1, url, wait)
                time.sleep(wait)
        self.dead_hosts[host] = "unreachable after retries"
        LOG.warning("DOWNLOAD_HOST_UNREACHABLE: %s did not respond after %d attempts; skipped for this run",
                    host, settings.HTTP_RETRIES)
        return False

    def _log(self, url, status, nbytes, dur, err):
        rec = {"url": url, "status": status, "bytes": nbytes, "duration_s": round(dur, 2), "error": err}
        self.attempts.append(rec)
        LOG.debug("GET %s -> %s (%d bytes, %.1f s) %s", url, status, nbytes, dur, err)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_meta(path, meta):
    with open(path + ".meta.json", "w") as fh:
        json.dump(meta, fh, indent=1, default=str)


def read_meta(path):
    try:
        with open(path + ".meta.json") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def _now_iso():
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# =========================================================================================== naming
def day_tokens(day_ns):
    y, doy = ts.year_doy(day_ns)
    w, d = ts.gps_week_dow(day_ns)
    yy, mm, dd, *_ = ts.from_ns(day_ns)
    return {"yyyy": f"{y:04d}", "ddd": f"{doy:03d}", "wwww": f"{w:04d}", "d": str(d),
            "yy": f"{y % 100:02d}", "mm": f"{mm:02d}", "dd": f"{dd:02d}"}


def candidate_names(ptype, family, day_ns):
    """Candidate filenames (with compression variants) for product type / family / GPST day."""
    tok = day_tokens(day_ns)
    tok["family"] = family
    names = []
    if (int(tok["yyyy"]), int(tok["mm"]), int(tok["dd"])) >= tuple(settings.LONG_NAME_START):
        tmpl = settings.LONG_NAMES.get(ptype)
        if tmpl:
            base = tmpl.format(**tok)
            names += [base + c for c in settings.COMPRESSION_SUFFIXES]
    else:
        tmpl = settings.LEGACY_NAMES.get(family, {}).get(ptype)
        if tmpl:
            base = tmpl.format(**tok)
            names += [base + c for c in (".Z", ".gz", "")]
    return names


# ======================================================================================= verification
def verify_file(path, ptype, family, day_ns):
    """Verify a product file before use (TDS § 12.4). Returns dict with ok flag, reasons and summary."""
    res = {"ok": False, "reasons": [], "summary": {}}
    try:
        orbclk.open_text(path)                # decompression test
    except Exception as exc:                  # noqa: BLE001 - any decompression failure is a verification failure
        res["reasons"].append(f"decompression failed: {exc}")
        return res
    d0, d1 = day_ns, day_ns + ts.DAY_NS
    try:
        if ptype == "SP3":
            sp = orbclk.read_sp3(path)
            gps = [s for s in sp.sats if s[0] == "G"]
            res["summary"] = {"agency": sp.agency, "frame": sp.frame, "first": ts.iso(sp.epochs[0]),
                              "last": ts.iso(sp.epochs[-1]), "n_gps": len(gps), "interval_s": sp.interval,
                              "time_system": sp.time_system}
            if not sp.agency.upper().startswith("COD"):
                res["reasons"].append(f"agency {sp.agency} is not COD")
            if sp.epochs[0] > d0 or sp.epochs[-1] < d1 - 600 * ts.NS:
                res["reasons"].append("epoch coverage incomplete")
            if len(gps) < settings.MIN_GPS_SATS_IN_PRODUCT:
                res["reasons"].append(f"only {len(gps)} GPS satellites")
            if sp.time_system not in ("GPS",):
                res["reasons"].append(f"time system {sp.time_system}")
        elif ptype in ("CLK", "CLK05S"):
            ck = orbclk.read_clk(path)
            firsts = [v[0][0] for v in ck.sats.values() if len(v[0])]
            lasts = [v[0][-1] for v in ck.sats.values() if len(v[0])]
            res["summary"] = {"agency": ck.agency, "pcvs": ck.pcvs, "trf": ck.trf, "n_gps": len(ck.sats),
                              "interval_s": ck.interval, "version": ck.version}
            if not ck.agency.upper().startswith("COD"):
                res["reasons"].append(f"analysis centre {ck.agency} is not COD")
            if len(ck.sats) < settings.MIN_GPS_SATS_IN_PRODUCT:
                res["reasons"].append(f"only {len(ck.sats)} GPS satellites")
            if not firsts or min(firsts) > d0 or max(lasts) < d1 - 60 * ts.NS:
                res["reasons"].append("epoch coverage incomplete")
        elif ptype == "ERP":
            er = orbclk.read_erp(path)
            res["summary"] = {"mjd": [float(er.mjd[0]), float(er.mjd[-1])]}
        elif ptype == "OSB":
            bi = orbclk.read_bias(path)
            sig = sorted({c for (_p, c) in bi.osb if _p.startswith("G")})
            res["summary"] = {"agency": bi.agency, "gps_signals": sig,
                              "n_gps": len({p for (p, _c) in bi.osb if p.startswith("G")})}
            if bi.agency and not bi.agency.upper().startswith("COD"):
                res["reasons"].append(f"agency {bi.agency} is not COD")
        elif ptype == "ATT":
            at = orbclk.read_obx(path)
            res["summary"] = {"n_sat": len(at.sats)}
            if not at.sats:
                res["reasons"].append("no attitude records")
    except Exception as exc:                  # noqa: BLE001 - unparseable header/body is a verification failure
        res["reasons"].append(f"parse error: {exc}")
    res["ok"] = not res["reasons"]
    return res


# ========================================================================================= fetching
def cache_path(cache_dir, family, day_ns, filename):
    tok = day_tokens(day_ns)
    return os.path.join(cache_dir, "COD", family, tok["yyyy"], tok["ddd"], filename)


def fetch_product(dl, ptype, family, day_ns):
    """Cache -> sources. Returns (path, meta) of a verified file, or (None, reason)."""
    names = candidate_names(ptype, family, day_ns)
    if not names:
        return None, "no filename template"
    # 1. cache
    for nm in names:
        p = cache_path(dl.cache_dir, family, day_ns, nm)
        if os.path.exists(p):
            meta = read_meta(p)
            if meta is None or meta.get("sha256") != sha256(p):
                ver = verify_file(p, ptype, family, day_ns)
                meta = {"filename": nm, "sha256": sha256(p), "source_url": (meta or {}).get("source_url", "cache"),
                        "retrieved": (meta or {}).get("retrieved", _now_iso()), "verification": ver,
                        "family": family, "type": ptype}
                write_meta(p, meta)
            if meta["verification"]["ok"]:
                return p, meta
            LOG.warning("PRODUCT_VERIFY_FAILED: cached %s rejected (%s)", nm, "; ".join(meta["verification"]["reasons"]))
    if dl.offline:
        return None, "not in cache (offline)"
    # 2. sources
    tok = day_tokens(day_ns)
    tok["family"] = family
    for nm in names:
        for tmpl in settings.PRODUCT_SOURCES:
            url = tmpl.format(filename=nm, **tok)
            dest = cache_path(dl.cache_dir, family, day_ns, nm)
            tmp_dest = dest + ".new"
            if not dl.get(url, tmp_dest):
                continue
            ver = verify_file(tmp_dest, ptype, family, day_ns)
            h = sha256(tmp_dest)
            if os.path.exists(dest):
                old = read_meta(dest) or {}
                if old.get("sha256") and old["sha256"] != h:
                    ver_dir = dest + f".v{_dt.datetime.now(_dt.timezone.utc).strftime('%Y%m%d%H%M%S')}"
                    shutil.move(dest, ver_dir)
                    LOG.warning("PRODUCT_CHANGED: remote %s changed (sha256 differs); old copy kept as %s",
                                nm, os.path.basename(ver_dir))
            os.replace(tmp_dest, dest)
            meta = {"filename": nm, "sha256": h, "source_url": url, "retrieved": _now_iso(),
                    "verification": ver, "family": family, "type": ptype}
            write_meta(dest, meta)
            if ver["ok"]:
                LOG.info("  downloaded %s", nm)
                return dest, meta
            LOG.warning("PRODUCT_VERIFY_FAILED: %s from %s rejected (%s)", nm, url, "; ".join(ver["reasons"]))
    return None, "not found at any source"


# ========================================================================================= resolution
@dataclass
class ProductSet:
    status: str = "FAIL"
    reason: str = ""
    family: str = ""
    tier: str = ""
    files: list = field(default_factory=list)          # [{type, filename, sha256, source_url, day, ...}]
    missing: list = field(default_factory=list)
    fallbacks: list = field(default_factory=list)
    flags: set = field(default_factory=set)
    sp3: object = None
    clk: object = None
    erp: object = None
    osb: object = None
    att: object = None
    frame_label: str = ""
    antex_name: str = ""
    orbit_edge_days: list = field(default_factory=list)
    window_end_limit: int = None                       # window shortened to this GPST ns if D+1 missing


def product_days(win_start, win_end):
    """GPST days intersecting [win_start, win_end] (core) and the neighbour days for SP3 interpolation."""
    d0 = ts.day_start(win_start)
    d1 = ts.day_start(win_end - 1)
    core = list(range(d0, d1 + 1, ts.DAY_NS))
    return core, [d0 - ts.DAY_NS] + core + [d1 + ts.DAY_NS]


def needs_osb(code_obs_used):
    """Float mode needs OSBs unless the code observables are exactly the clock-reference signals."""
    ref = set(settings.CLOCK_REFERENCE_CODES["G"])
    return any(c not in ref for c in code_obs_used)


def resolve_products(dl, win_start, win_end, code_obs_used=("C1W", "C2W"), ar_requested=False, now_ns=None):
    """Search families in order; select the first complete one; evaluate status (TDS § 11.4, § 12.2)."""
    core, with_nb = product_days(win_start, win_end)
    want_osb = needs_osb(code_obs_used) or ar_requested
    ps = ProductSet()
    candidates = []
    for fam in settings.PRODUCT_FAMILIES:
        if fam in settings.FORBIDDEN_FAMILIES:
            continue
        LOG.info("Products: trying family %s (%s)", fam, settings.FAMILY_TIER[fam])
        got, miss = {}, []
        for ptype in ["SP3", "CLK", "ERP"] + (["OSB"] if "OSB" in settings.FAMILY_PRODUCTS[fam] else []) + \
                (["ATT"] if "ATT" in settings.FAMILY_PRODUCTS[fam] else []):
            days = with_nb if ptype == "SP3" else core
            for d in days:
                p, meta = fetch_product(dl, ptype, fam, d)
                if p:
                    got[(ptype, d)] = (p, meta)
                else:
                    miss.append((ptype, ts.iso(d)[:10], meta))
        core_ok = all((pt, d) in got for pt in ("SP3", "CLK", "ERP") for d in core)
        osb_ok = all(("OSB", d) in got for d in core)
        nb_ok = all(("SP3", d) in got for d in with_nb)
        n_found = len(got)
        LOG.info("  %s: %d files verified, core %s, neighbour-day orbits %s, OSB %s", fam, n_found,
                 "complete" if core_ok else "INCOMPLETE", "yes" if nb_ok else "no", "yes" if osb_ok else "no")
        if core_ok and (osb_ok or not want_osb):
            candidates.append((fam, got, miss, nb_ok, osb_ok))
            if nb_ok:
                break
        elif core_ok and want_osb and not osb_ok:
            ps.fallbacks.append(f"{fam}: OSB missing for code observables {list(code_obs_used)}")
            candidates.append((fam, got, miss, nb_ok, osb_ok))
        else:
            ps.missing += [f"{fam}:{m[0]}:{m[1]}" for m in miss if m[0] in ("SP3", "CLK", "ERP")]
    if not candidates:
        now_ns = now_ns or ts.utc_to_gpst(_utc_now_ns())
        age_h = (now_ns - win_end) / ts.NS / 3600.0
        max_lat = max(settings.MAX_EXPECTED_LATENCY_H.values())
        if age_h < max_lat:
            ps.status, ps.reason = "WAIT_FOR_TIER", f"no PPP-capable CODE family complete; data only {age_h:.0f} h old"
        else:
            ps.status, ps.reason = "FAIL", "NO_CODE_PRODUCTS: no PPP-capable CODE family complete"
        return ps
    # prefer a candidate with neighbour-day orbits (whole window in the lowest common tier)
    best = next((c for c in candidates if c[3] and (c[4] or not want_osb)), None) or candidates[0]
    fam, got, miss, nb_ok, osb_ok = best
    ps.family, ps.tier = fam, settings.FAMILY_TIER[fam]
    if not nb_ok:
        ps.flags.add("ORBIT_EDGE")
        ps.orbit_edge_days = [ts.iso(d)[:10] for d in with_nb if ("SP3", d) not in got]
        LOG.warning("ORBIT_EDGE: adjacent-day orbits not available in %s (%s); asymmetric SP3 interpolation at "
                    "window edges", fam, ", ".join(ps.orbit_edge_days))
    if want_osb and not osb_ok:
        ps.flags.add("CODE_BIAS_MISSING")
        LOG.warning("CODE_BIAS_MISSING: OSB file missing in %s; code observables %s are used without satellite "
                    "code-bias correction", fam, list(code_obs_used))
    if fam != settings.PRODUCT_FAMILIES[0]:
        ps.flags.add("RAPID_TIER")
    # load
    def _load(ptype, reader, merger, days):
        paths = [got[(ptype, d)][0] for d in days if (ptype, d) in got]
        if not paths:
            return None
        return merger([reader(p) for p in paths])
    ps.sp3 = _load("SP3", orbclk.read_sp3, orbclk.merge_sp3, with_nb)
    ps.clk = _load("CLK", orbclk.read_clk, orbclk.merge_clk, core)
    ps.erp = _load("ERP", orbclk.read_erp, orbclk.merge_erp, core)
    ps.osb = _load("OSB", orbclk.read_bias, orbclk.merge_bias, core)
    ps.att = _load("ATT", orbclk.read_obx, orbclk.merge_att, core)
    for (ptype, d), (p, meta) in sorted(got.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        ps.files.append({"type": ptype, "filename": os.path.basename(p), "sha256": meta["sha256"],
                         "source_url": meta.get("source_url"), "day": ts.iso(d)[:10],
                         "retrieved": meta.get("retrieved"), "header": meta["verification"].get("summary", {})})
    # frame label / ANTEX consistency (TDS § 9.3, § 12.4.6)
    frames = getattr(ps.sp3, "frames", {ps.sp3.frame})
    ps.frame_label = ps.sp3.frame
    if len(frames) > 1:
        ps.status, ps.reason = "FAIL", f"FRAME_INCONSISTENT: SP3 frame labels differ {sorted(frames)}"
        return ps
    pcvs = {p for p in getattr(ps.clk, "pcvs_all", {ps.clk.pcvs}) if p}
    ps.antex_name = _antex_name_from_pcvs(ps.clk.pcvs) if ps.clk.pcvs else ""
    if len({_antex_name_from_pcvs(p) for p in pcvs}) > 1:
        LOG.warning("ANTEX_INCONSISTENT: clock files name different ANTEX versions %s", sorted(pcvs))
    # status
    y, m, d, *_ = ts.from_ns(win_start)
    reasons = []
    if not ar_requested:
        reasons.append("AR_NOT_REQUESTED")
    if (y, m, d) < tuple(settings.INTEGER_CLOCK_START):
        reasons.append("PRE_2018")
    if ps.osb is None or not any(c.startswith("L") for (_p, c) in ps.osb.osb):
        reasons.append("NO_PHASE_OSB")
    if ps.att is None:
        reasons.append("ATT_MISSING_ECLIPSE_EXCLUDED")
    hard = [r for r in reasons if r != "ATT_MISSING_ECLIPSE_EXCLUDED"]
    if not hard:
        ps.status, ps.reason = "OK_AR", ""
    else:
        ps.status, ps.reason = "OK_FLOAT_ONLY", ",".join(reasons)
    return ps


def _antex_name_from_pcvs(pcvs):
    """ANTEX name such as 'igs20_2290' found anywhere in a SYS / PCVS APPLIED field (verify_log V-009)."""
    import re
    m = re.search(r"(igs\d\d(?:_\d{4})?)(?:\.atx)?", pcvs or "", re.IGNORECASE)
    return m.group(1).lower() if m else ""


def _utc_now_ns():
    n = _dt.datetime.now(_dt.timezone.utc)
    return ts.to_ns(n.year, n.month, n.day, n.hour, n.minute, n.second + n.microsecond * 1e-6)


# ================================================================================= non-AC input files
def _fetch_simple(dl, subdir, filename, urls, max_age_days=None, validate=None):
    """Generic cache-or-download for convention/model files (ANTEX, GPT3, VMF3, leap seconds)."""
    dest = os.path.join(dl.cache_dir, subdir, filename)
    if os.path.exists(dest):
        meta = read_meta(dest) or {}
        fresh = True
        if max_age_days and meta.get("retrieved"):
            try:
                age = (_dt.datetime.now(_dt.timezone.utc) -
                       _dt.datetime.strptime(meta["retrieved"], "%Y-%m-%dT%H:%M:%SZ").replace(
                           tzinfo=_dt.timezone.utc)).days
                fresh = age <= max_age_days
            except ValueError:
                pass
        if fresh or dl.offline:
            return dest, meta
    for url in urls:
        tmp = dest + ".new"
        if dl.get(url, tmp):
            if validate is not None:
                try:
                    validate(tmp)
                except Exception as exc:          # noqa: BLE001 - reject any unparseable file
                    LOG.warning("DOWNLOAD_INVALID: %s rejected (%s)", url, exc)
                    os.remove(tmp)
                    continue
            os.replace(tmp, dest)
            meta = {"filename": filename, "sha256": sha256(dest), "source_url": url, "retrieved": _now_iso()}
            write_meta(dest, meta)
            return dest, meta
    if os.path.exists(dest):
        return dest, read_meta(dest) or {}
    return None, None


def get_leap_seconds(dl):
    """Download/refresh the leap-second file monthly; fall back to the built-in table (with warning)."""
    for url in settings.LEAP_SECOND_URLS:
        name = os.path.basename(url)
        p, meta = _fetch_simple(dl, "LEAPSEC", name, [url], settings.LEAP_SECOND_REFRESH_DAYS,
                                validate=lambda f: ts.parse_leap_file(open(f).read()))
        if p:
            ts.set_leap_table(ts.parse_leap_file(open(p).read()), f"{name} sha256={meta.get('sha256', '')[:12]}")
            return {"source": name, "sha256": meta.get("sha256"), "url": meta.get("source_url")}
    ts.set_leap_table(ts._builtin_table(), "builtin (settings.LEAP_SECONDS_BUILTIN)")
    return {"source": "builtin", "sha256": None, "url": None}


def get_antex(dl, wanted_name):
    """ANTEX named in the clock header (exact version, pcv_archive) else current file of the same family.

    Returns (path, meta, status) with status 'EXACT' | 'FAMILY_FALLBACK' | None.
    """
    wanted = (wanted_name or settings.ANTEX_DEFAULT_NAME).lower()
    family = wanted.split("_")[0] if wanted else settings.ANTEX_DEFAULT_NAME
    val = lambda f: _check_antex(f)  # noqa: E731
    if "_" in wanted:
        urls = [u.format(name=wanted, family=family) for u in settings.ANTEX_SOURCES if "{name}" in u]
        p, meta = _fetch_simple(dl, "ANTEX", wanted + ".atx", urls, validate=val)
        if p:
            return p, meta, "EXACT"
    # any cached file of the same family (newest first)
    adir = os.path.join(dl.cache_dir, "ANTEX")
    urls = [u.format(name=wanted, family=family) for u in settings.ANTEX_SOURCES if "{name}" not in u]
    p, meta = _fetch_simple(dl, "ANTEX", family + ".atx", urls, max_age_days=7, validate=val)
    if p:
        return p, meta, ("EXACT" if wanted == family else "FAMILY_FALLBACK")
    if os.path.isdir(adir):
        cands = sorted([f for f in os.listdir(adir) if f.lower().startswith(family) and f.endswith(".atx")],
                       reverse=True)
        if cands:
            p = os.path.join(adir, cands[0])
            return p, read_meta(p) or {"sha256": sha256(p)}, (
                "EXACT" if os.path.splitext(cands[0])[0].lower() == wanted else "FAMILY_FALLBACK")
    return None, None, None


def _check_antex(path):
    with open(path, "rb") as fh:
        head = fh.read(400).decode("ascii", "replace")
    if "ANTEX VERSION / SYST" not in head:
        raise ValueError("not an ANTEX file")


def get_gpt3(dl):
    return _fetch_simple(dl, "GPT3", "gpt3_1.grd", [settings.GPT3_URL], validate=_check_gpt3)


def _check_gpt3(path):
    with open(path) as fh:
        head = fh.readline()
    if "lat" not in head or "lon" not in head:
        raise ValueError("not a GPT3 grid file")


def get_vmf3(dl, win_start, win_end):
    """All 6-hourly VMF3 1x1 grid files covering [start-6h, end+6h] + orography (TDS § 12.1.4)."""
    t0 = win_start - 6 * 3600 * ts.NS
    t1 = win_end + 6 * 3600 * ts.NS
    t = t0 - t0 % (6 * 3600 * ts.NS)
    files, missing = [], []
    while t <= t1:
        y, m, d, hh, *_ = ts.from_ns(t)
        tok = {"yyyy": f"{y:04d}", "mm": f"{m:02d}", "dd": f"{d:02d}", "hh": f"{hh:02d}"}
        url = settings.VMF3_URL.format(**tok)
        name = os.path.basename(url)
        p, meta = _fetch_simple(dl, os.path.join("VMF3", "1x1_OP", tok["yyyy"]), name, [url],
                                validate=_check_vmf3)
        if p:
            files.append((t, p, meta))
        else:
            missing.append(name)
        t += 6 * 3600 * ts.NS
    oro, ometa = _fetch_simple(dl, "VMF3", "orography_ell_1x1", [settings.VMF3_OROGRAPHY_URL])
    return files, missing, oro


def _check_vmf3(path):
    with open(path) as fh:
        txt = fh.read(3000)
    if "!" not in txt:
        raise ValueError("not a VMF3 grid file")


# ============================================================================================ manifest
def config_hash():
    snap = settings.snapshot()
    canon = json.dumps(snap, sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.sha256(canon.encode()).hexdigest(), snap


def software_version():
    v = settings.SOFTWARE_VERSION
    try:
        import subprocess
        here = os.path.dirname(os.path.abspath(__file__))
        commit = subprocess.run(["git", "-C", here, "rev-parse", "--short", "HEAD"], capture_output=True,
                                text=True, timeout=5).stdout.strip()
        dirty = subprocess.run(["git", "-C", here, "status", "--porcelain"], capture_output=True, text=True,
                               timeout=5).stdout.strip()
        if commit:
            v += f"+{commit}" + (".dirty" if dirty else "")
    except (OSError, ValueError):
        pass
    return v


def new_manifest(station, day_str):
    return {"manifest_id": str(uuid.uuid4()), "station": station, "day": day_str,
            "analysis_centre": "COD", "run_timestamp": _now_iso(), "software_version": software_version(),
            "supersedes": None, "superseded_by": None}


def write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, default=_json_default)
    os.replace(tmp, path)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, set):
        return sorted(o)
    return str(o)


def find_superseded(out_dir, station, day_str, new_tier, new_manifest_id):
    """Plain-JSON supersession (TDS § 19.3): mark older manifests of the same station-day."""
    tiers = {"FINAL": 3, "RAPID_M": 2, "RAPID_0": 1}
    old_ids = []
    if not os.path.isdir(out_dir):
        return old_ids
    for f in os.listdir(out_dir):
        if not f.endswith("_manifest.json") or not f.startswith(station):
            continue
        p = os.path.join(out_dir, f)
        try:
            with open(p) as fh:
                m = json.load(fh)
        except (OSError, ValueError):
            continue
        if m.get("day") == day_str and m.get("manifest_id") != new_manifest_id and not m.get("superseded_by"):
            if tiers.get(new_tier, 0) >= tiers.get(m.get("tier"), 0):
                m["superseded_by"] = new_manifest_id
                write_json(p, m)
                old_ids.append(m.get("manifest_id"))
    return old_ids

