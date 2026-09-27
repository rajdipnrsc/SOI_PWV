"""Readers and interpolation for SP3, Clock RINEX, IGS ERP, Bias-SINEX and ORBEX; GPS broadcast orbits.

TDS § 7.2 (SP3 Lagrange, 11 nodes), § 7.4 (clock convention), § 7.14 (ERP), § 11.5 (OSB), § 7.11 (ORBEX).
All times are integer ns GPST (see timesys).  Positions in metres, clocks in seconds.
"""
import gzip
import os
from dataclasses import dataclass, field

import numpy as np

from . import settings
from . import timesys as ts

C = settings.C_LIGHT


def open_text(path):
    """Read a (possibly gzip / Unix-compress) product file as text."""
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    elif raw[:2] == b"\x1f\x9d":
        try:
            import ncompress  # installed with hatanaka
            raw = ncompress.decompress(raw)
        except ImportError:
            import hatanaka
            raw = hatanaka.decompress(raw)
    return raw.decode("ascii", errors="replace")


# ============================================================================================ SP3
@dataclass
class SP3:
    epochs: np.ndarray            # int64 ns GPST
    sats: list
    pos: np.ndarray               # (n_ep, n_sat, 3) m, NaN = bad/missing
    clk: np.ndarray               # (n_ep, n_sat) s, NaN = bad
    maneuver: np.ndarray          # (n_ep, n_sat) bool
    frame: str = ""
    agency: str = ""
    time_system: str = "GPS"
    interval: float = 0.0
    version: str = ""
    files: list = field(default_factory=list)


def read_sp3(path):
    txt = open_text(path)
    lines = txt.splitlines()
    if not lines or lines[0][0] != "#":
        raise ValueError(f"{path}: not an SP3 file")
    version = lines[0][1]
    frame = lines[0][46:51].strip()
    agency = lines[0][56:60].strip()
    interval = float(lines[1][24:38]) if lines[1].startswith("##") else 0.0
    tsys = "GPS"
    for L in lines[:30]:
        if L.startswith("%c"):
            tsys = L[9:12].strip() or "GPS"
            break
    epochs, recs = [], []
    cur = None
    for L in lines:
        if L.startswith("* "):
            t = ts.to_ns(int(L[3:7]), int(L[8:10]), int(L[11:13]), int(L[14:16]), int(L[17:19]), float(L[20:31]))
            epochs.append(t)
            cur = {}
            recs.append(cur)
        elif L.startswith("P") and cur is not None:
            sat = L[1:4].replace(" ", "0")
            try:
                x, y, z = float(L[4:18]), float(L[18:32]), float(L[32:46])
                ck = float(L[46:60]) if L[46:60].strip() else 999999.999999
            except ValueError:
                continue
            man = len(L) > 78 and L[78] == "M"
            cur[sat] = (x, y, z, ck, man)
        elif L.startswith("EOF"):
            break
    sats = sorted({s for r in recs for s in r})
    n, m = len(epochs), len(sats)
    pos = np.full((n, m, 3), np.nan)
    clk = np.full((n, m), np.nan)
    man = np.zeros((n, m), dtype=bool)
    si = {s: j for j, s in enumerate(sats)}
    for i, r in enumerate(recs):
        for s, (x, y, z, ck, mv) in r.items():
            j = si[s]
            if not (x == 0.0 and y == 0.0 and z == 0.0):
                pos[i, j] = (x * 1e3, y * 1e3, z * 1e3)
            if abs(ck) < 999999.0:
                clk[i, j] = ck * 1e-6
            man[i, j] = mv
    return SP3(np.array(epochs, dtype=np.int64), sats, pos, clk, man, frame, agency, tsys, interval,
               version, [os.path.basename(path)])


def merge_sp3(parts):
    """Concatenate consecutive SP3 files (same family) into one time series; duplicate epochs dropped."""
    parts = [p for p in parts if p is not None]
    sats = sorted({s for p in parts for s in p.sats})
    si = {s: j for j, s in enumerate(sats)}
    ep_all = np.unique(np.concatenate([p.epochs for p in parts]))
    ei = {int(t): i for i, t in enumerate(ep_all)}
    pos = np.full((len(ep_all), len(sats), 3), np.nan)
    clk = np.full((len(ep_all), len(sats)), np.nan)
    man = np.zeros((len(ep_all), len(sats)), dtype=bool)
    for p in parts:
        rows = [ei[int(t)] for t in p.epochs]
        cols = [si[s] for s in p.sats]
        for a, r in enumerate(rows):
            good = np.isfinite(p.pos[a, :, 0])
            for b, c in enumerate(cols):
                if good[b] and not np.isfinite(pos[r, c, 0]):
                    pos[r, c] = p.pos[a, b]
                    clk[r, c] = p.clk[a, b]
                man[r, c] |= p.maneuver[a, b]
    frames = {p.frame for p in parts}
    out = SP3(ep_all, sats, pos, clk, man, parts[0].frame, parts[0].agency, parts[0].time_system,
              parts[0].interval, parts[0].version, [f for p in parts for f in p.files])
    out.frames = frames
    return out


def lagrange_weights(x_nodes, x):
    """Lagrange basis weights. x_nodes (n, k) or (k,), x (n,) or scalar -> weights (n, k)."""
    xn = np.atleast_2d(np.asarray(x_nodes, dtype=float))
    x = np.atleast_1d(np.asarray(x, dtype=float))
    d = x[:, None] - xn
    k = xn.shape[1]
    w = np.ones((max(len(x), xn.shape[0]), k))
    for j in range(k):
        for m in range(k):
            if m != j:
                w[:, j] *= d[:, m] / (xn[:, j] - xn[:, m])
    return w


def sp3_interp(sp3, sat, t_ns, nodes=None, want_vel=True):
    """Satellite CoM position (m, ECEF of product frame) at times t_ns (array), TDS § 7.2.

    Lagrange interpolation with `nodes` points (default 11, degree 10) centred on t; near the ends of
    the available data the node set is shifted (asymmetric) and edge=True is returned for those times.
    Node sets spanning a data gap (> 1.5 nominal steps) are rejected (no value).
    Velocity from a central difference of the interpolating polynomial (h = 0.5 s).
    Returns pos (n,3), vel (n,3), ok (n,) bool, edge (n,) bool.
    """
    nodes = nodes or settings.SP3_LAGRANGE_NODES
    t_ns = np.atleast_1d(np.asarray(t_ns, dtype=np.int64))
    n = len(t_ns)
    pos = np.full((n, 3), np.nan)
    vel = np.full((n, 3), np.nan)
    ok = np.zeros(n, dtype=bool)
    edge = np.zeros(n, dtype=bool)
    if sat not in sp3.sats:
        return pos, vel, ok, edge
    j = sp3.sats.index(sat)
    good = np.isfinite(sp3.pos[:, j, 0])
    ep = sp3.epochs[good]
    if len(ep) < nodes:
        return pos, vel, ok, edge
    P = sp3.pos[good, j, :]
    t0 = int(ep[0])
    xs = (ep - t0) / ts.NS
    dt_nom = np.median(np.diff(xs))
    half = nodes // 2
    q = (t_ns - t0) / ts.NS
    inside = (q >= xs[0] - 1e-6) & (q <= xs[-1] + 1e-6)
    k = np.searchsorted(xs, q)
    centre = k - half
    lo = np.clip(centre, 0, len(xs) - nodes)
    idx = lo[:, None] + np.arange(nodes)[None, :]
    xn = xs[idx]
    gap_ok = np.max(np.diff(xn, axis=1), axis=1) <= 1.5 * dt_nom
    use = inside & gap_ok
    if not use.any():
        return pos, vel, ok, edge
    iu = np.nonzero(use)[0]
    rel = xn[iu] - q[iu, None]                 # nodes relative to query time
    w = lagrange_weights(rel, np.zeros(len(iu)))
    Pn = P[idx[iu]]                            # (m, nodes, 3)
    pos[iu] = np.einsum("ij,ijk->ik", w, Pn)
    if want_vel:
        wp = lagrange_weights(rel, np.full(len(iu), 0.5))
        wm = lagrange_weights(rel, np.full(len(iu), -0.5))
        vel[iu] = np.einsum("ij,ijk->ik", wp - wm, Pn)
    ok[iu] = True
    edge[iu] = (centre[iu] < 0) | (centre[iu] > len(xs) - nodes)
    return pos, vel, ok, edge


class OrbitInterpolator:
    """Fast repeated interpolation: caches per-satellite node arrays."""

    def __init__(self, sp3):
        self.sp3 = sp3
        self.nodes = settings.SP3_LAGRANGE_NODES

    def __call__(self, sat, t_ns, want_vel=True):
        return sp3_interp(self.sp3, sat, t_ns, self.nodes, want_vel)


# ============================================================================================ CLK
@dataclass
class ClockProduct:
    sats: dict                    # sat -> (epochs int64 ns, bias s)
    version: float = 0.0
    agency: str = ""
    pcvs: str = ""                # ANTEX named in SYS / PCVS APPLIED
    trf: str = ""                 # # OF SOLN STA / TRF label
    time_system: str = "GPS"
    interval: float = 30.0
    files: list = field(default_factory=list)


def read_clk(path, systems=("G",)):
    txt = open_text(path)
    lines = txt.splitlines()
    version, agency, pcvs, trf, tsys = 0.0, "", "", "", "GPS"
    i = 0
    for i, L in enumerate(lines):
        lab = L[60:80].strip() if len(L) > 60 else ""
        if "RINEX VERSION / TYPE" in L:
            version = float(L[0:9])
        elif lab == "ANALYSIS CENTER":
            agency = L[0:3].strip()
        elif lab == "SYS / PCVS APPLIED":
            if L[0] in systems or not pcvs:
                pcvs = L[1:60].strip()          # program + source (ANTEX name located later by pattern)
        elif lab == "# OF SOLN STA / TRF":
            trf = L[10:60].strip()
        elif lab == "TIME SYSTEM ID":
            tsys = L[3:6].strip() or "GPS"
        elif "END OF HEADER" in L:
            break
    data = {}
    for L in lines[i + 1:]:
        if not L.startswith("AS"):
            continue
        p = L.split()
        sat = p[1]
        if sat[0] not in systems:
            continue
        if len(sat) == 2:
            sat = sat[0] + "0" + sat[1]
        t = ts.to_ns(int(p[2]), int(p[3]), int(p[4]), int(p[5]), int(p[6]), float(p[7]))
        data.setdefault(sat, ([], []))
        data[sat][0].append(t)
        data[sat][1].append(float(p[9]))
    out = {}
    iv = []
    for s, (t, v) in data.items():
        t = np.array(t, dtype=np.int64)
        v = np.array(v)
        o = np.argsort(t)
        out[s] = (t[o], v[o])
        if len(t) > 2:
            iv.append(np.median(np.diff(t[o])) / ts.NS)
    return ClockProduct(out, version, agency, pcvs, trf, tsys, float(np.median(iv)) if iv else 30.0,
                        [os.path.basename(path)])


def merge_clk(parts):
    sats = {}
    for p in parts:
        for s, (t, v) in p.sats.items():
            sats.setdefault(s, ([], []))
            sats[s][0].append(t)
            sats[s][1].append(v)
    out = {}
    for s, (tl, vl) in sats.items():
        t = np.concatenate(tl)
        v = np.concatenate(vl)
        t, idx = np.unique(t, return_index=True)
        out[s] = (t, v[idx])
    p0 = parts[0]
    res = ClockProduct(out, p0.version, p0.agency, p0.pcvs, p0.trf, p0.time_system, p0.interval,
                       [f for p in parts for f in p.files])
    res.pcvs_all = {p.pcvs for p in parts}
    res.trf_all = {p.trf for p in parts}
    return res


def clk_interp(clk, sat, t_ns, max_gap_s=None):
    """Satellite clock bias (s) at t. Exact at nodes; linear between nodes otherwise (flag interp).

    Returns (value array, ok bool array, interp bool array). No extrapolation; gaps > max_gap_s rejected.
    """
    t_ns = np.atleast_1d(np.asarray(t_ns, dtype=np.int64))
    val = np.full(len(t_ns), np.nan)
    ok = np.zeros(len(t_ns), dtype=bool)
    interp = np.zeros(len(t_ns), dtype=bool)
    if sat not in clk.sats:
        return val, ok, interp
    t, v = clk.sats[sat]
    if len(t) == 0:
        return val, ok, interp
    max_gap = int((max_gap_s or 2.5 * clk.interval) * ts.NS)
    k = np.searchsorted(t, t_ns)
    kc = np.clip(k, 0, len(t) - 1)
    exact = t[kc] == t_ns
    val[exact] = v[kc[exact]]
    ok[exact] = True
    a = np.clip(k - 1, 0, len(t) - 1)
    b = kc
    cand = (~exact) & (k > 0) & (k < len(t)) & ((t[b] - t[a]) <= max_gap) & (t[b] > t[a])
    f = (t_ns[cand] - t[a[cand]]) / (t[b[cand]] - t[a[cand]])
    val[cand] = v[a[cand]] + f * (v[b[cand]] - v[a[cand]])
    ok[cand] = True
    interp[cand] = True
    return val, ok, interp


# ============================================================================================ ERP
@dataclass
class ERP:
    mjd: np.ndarray               # UTC-ish MJD of values (IGS ERP uses MJD at 12h/0h)
    xp: np.ndarray                # arcsec
    yp: np.ndarray                # arcsec
    ut1_utc: np.ndarray           # s
    files: list = field(default_factory=list)


def read_erp(path):
    txt = open_text(path)
    mj, xp, yp, du = [], [], [], []
    started = False
    for L in txt.splitlines():
        p = L.split()
        if not p:
            continue
        if p[0].upper() == "MJD":
            started = True
            continue
        if not started:
            continue
        try:
            m = float(p[0])
        except ValueError:
            continue
        if m < 40000 or len(p) < 4:
            continue
        mj.append(m)
        xp.append(float(p[1]) * 1e-6)
        yp.append(float(p[2]) * 1e-6)
        du.append(float(p[3]) * 1e-7)
    if not mj:
        raise ValueError(f"{path}: no ERP records")
    o = np.argsort(mj)
    return ERP(np.array(mj)[o], np.array(xp)[o], np.array(yp)[o], np.array(du)[o], [os.path.basename(path)])


def merge_erp(parts):
    mj = np.concatenate([p.mjd for p in parts])
    mj, idx = np.unique(mj, return_index=True)
    return ERP(mj, np.concatenate([p.xp for p in parts])[idx], np.concatenate([p.yp for p in parts])[idx],
               np.concatenate([p.ut1_utc for p in parts])[idx], [f for p in parts for f in p.files])


def erp_interp(erp, mjd_utc):
    """Linear interpolation of x_p, y_p (arcsec) and UT1-UTC (s) (TDS § 7.14). Clamped at the ends."""
    return (np.interp(mjd_utc, erp.mjd, erp.xp), np.interp(mjd_utc, erp.mjd, erp.yp),
            np.interp(mjd_utc, erp.mjd, erp.ut1_utc))


# ============================================================================================ OSB
@dataclass
class BiasProduct:
    osb: dict                     # (prn, obs) -> list of (start_ns, end_ns, value_m)
    agency: str = ""
    files: list = field(default_factory=list)

    def get(self, prn, code, t_ns):
        lst = self.osb.get((prn, code))
        if not lst:
            return None
        for a, b, v in lst:
            if a <= t_ns < b or (b == a):
                return v
        return lst[0][2] if len(lst) == 1 else None

    def signals(self, prn):
        return {c for (p, c) in self.osb if p == prn}


def _sinex_time(s):
    y, d, sec = s.split(":")
    y = int(y)
    if y < 100:
        y += 2000 if y < 80 else 1900
    return ts.to_ns(y, 1, 1) + (int(d) - 1) * ts.DAY_NS + int(sec) * ts.NS


def read_bias(path):
    """Bias-SINEX 1.00: satellite OSBs (code and phase) converted to metres."""
    txt = open_text(path)
    osb = {}
    agency = ""
    inblk = False
    for L in txt.splitlines():
        if L.startswith("%=BIA"):
            agency = L[15:18]
        if L.startswith("+BIAS/SOLUTION"):
            inblk = True
            continue
        if L.startswith("-BIAS/SOLUTION"):
            break
        if not inblk or L.startswith("*") or not L.strip():
            continue
        if L[1:4] != "OSB":
            continue
        station = L[15:24].strip()
        if station:
            continue
        prn = L[11:14].strip()
        code = L[25:29].strip()
        try:
            start = _sinex_time(L[35:49].strip())
            end = _sinex_time(L[50:64].strip())
            unit = L[65:69].strip()
            val = float(L[70:91])
        except (ValueError, IndexError):
            p = L.split()
            if len(p) < 8:
                continue
            prn, code, start, end, unit, val = p[2], p[3], _sinex_time(p[4]), _sinex_time(p[5]), p[6], float(p[7])
        if unit == "ns":
            v = val * 1e-9 * C
        elif unit == "cyc":
            band = code[1]
            f = settings.FREQ.get((prn[0], band))
            if not f:
                continue
            v = val * C / f
        else:
            continue
        osb.setdefault((prn, code), []).append((start, end, v))
    return BiasProduct(osb, agency, [os.path.basename(path)])


def merge_bias(parts):
    osb = {}
    for p in parts:
        for k, v in p.osb.items():
            osb.setdefault(k, []).extend(v)
    return BiasProduct(osb, parts[0].agency, [f for p in parts for f in p.files])


# ============================================================================================ ORBEX
@dataclass
class Attitude:
    sats: dict                    # sat -> (epochs ns, quaternions (n,4) q0 scalar first)
    convention: str = "T2B"       # q maps terrestrial -> body: (0,B) = q (0,T) q*   [ORBEX 0.09]
    files: list = field(default_factory=list)


def read_obx(path):
    txt = open_text(path)
    cur = None
    data = {}
    inblk = False
    for L in txt.splitlines():
        if L.startswith("+EPHEMERIS/DATA"):
            inblk = True
            continue
        if L.startswith("-EPHEMERIS/DATA"):
            break
        if not inblk:
            continue
        if L.startswith("##"):
            p = L[2:].split()
            cur = ts.to_ns(int(p[0]), int(p[1]), int(p[2]), int(p[3]), int(p[4]), float(p[5]))
        elif L.startswith(" ATT") and cur is not None:
            p = L.split()
            sat = p[1]
            if len(p) < 7:
                continue
            q = [float(x) for x in p[3:7]]
            data.setdefault(sat, ([], []))
            data[sat][0].append(cur)
            data[sat][1].append(q)
    out = {s: (np.array(t, dtype=np.int64), np.array(q)) for s, (t, q) in data.items()}
    return Attitude(out, files=[os.path.basename(path)])


def merge_att(parts):
    sats = {}
    for p in parts:
        for s, (t, q) in p.sats.items():
            sats.setdefault(s, ([], []))
            sats[s][0].append(t)
            sats[s][1].append(q)
    out = {}
    for s, (tl, ql) in sats.items():
        t = np.concatenate(tl)
        q = np.concatenate(ql)
        t, idx = np.unique(t, return_index=True)
        out[s] = (t, q[idx])
    return Attitude(out, parts[0].convention, [f for p in parts for f in p.files])


def quat_to_matrix(q):
    """Rotation matrix R(q) such that v' = R v for the quaternion product (0,v') = q (0,v) q*."""
    w, x, y, z = q / np.linalg.norm(q)
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
                     [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
                     [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])


def att_axes(att, sat, t_ns, max_gap_s=90.0):
    """Body axes (ex, ey, ez) expressed in ECEF from ORBEX quaternions at t (nearest-neighbour slerp-free
    linear quaternion interpolation between 30-s nodes). Returns None if unavailable."""
    if sat not in att.sats:
        return None
    t, q = att.sats[sat]
    k = np.searchsorted(t, t_ns)
    if k < len(t) and t[k] == t_ns:
        qq = q[k]
    else:
        if k == 0 or k >= len(t) or (t[k] - t[k - 1]) > max_gap_s * ts.NS:
            return None
        f = (t_ns - t[k - 1]) / (t[k] - t[k - 1])
        q1, q2 = q[k - 1], q[k]
        if np.dot(q1, q2) < 0:
            q2 = -q2
        qq = (1 - f) * q1 + f * q2
    R = quat_to_matrix(qq)
    if att.convention == "T2B":
        M = R.T                   # body -> terrestrial
    else:
        M = R
    return M[:, 0], M[:, 1], M[:, 2]


# ============================================================================================ broadcast
MU_GPS = 3.986005e14              # IS-GPS-200 [F]
OMEGA_E_GPS = 7.2921151467e-5
F_REL = -4.442807633e-10          # s/sqrt(m) [F] IS-GPS-200


def select_eph(ephs, t_ns, max_age_h=4.0):
    best, bd = None, None
    for e in ephs:
        if e.health != 0:
            continue
        d = abs(t_ns - e.toe)
        if d <= max_age_h * 3600 * ts.NS and (bd is None or d < bd):
            best, bd = e, d
    return best


def broadcast_pos_clk(eph, t_ns, include_rel=True):
    """Satellite ECEF position (m, antenna phase centre per broadcast convention) and clock (s) at GPST
    t_ns from GPS LNAV ephemeris (IS-GPS-200 Table 20-IV). The clock excludes TGD (refers to the L1/L2
    P(Y) IF combination) and includes the relativistic term unless include_rel=False (IGS convention)."""
    tk = (t_ns - eph.toe) / ts.NS
    A = eph.sqrt_a ** 2
    n = np.sqrt(MU_GPS / A ** 3) + eph.delta_n
    M = eph.m0 + n * tk
    E = M
    for _ in range(30):
        dE = (M - E + eph.e * np.sin(E)) / (1 - eph.e * np.cos(E))
        E += dE
        if abs(dE) < 1e-14:
            break
    v = np.arctan2(np.sqrt(1 - eph.e ** 2) * np.sin(E), np.cos(E) - eph.e)
    phi = v + eph.omega
    du = eph.cus * np.sin(2 * phi) + eph.cuc * np.cos(2 * phi)
    dr = eph.crs * np.sin(2 * phi) + eph.crc * np.cos(2 * phi)
    di = eph.cis * np.sin(2 * phi) + eph.cic * np.cos(2 * phi)
    u = phi + du
    r = A * (1 - eph.e * np.cos(E)) + dr
    i = eph.i0 + di + eph.idot * tk
    xp, yp = r * np.cos(u), r * np.sin(u)
    Om = eph.omega0 + (eph.omega_dot - OMEGA_E_GPS) * tk - OMEGA_E_GPS * eph.toe_sow
    x = xp * np.cos(Om) - yp * np.cos(i) * np.sin(Om)
    y = xp * np.sin(Om) + yp * np.cos(i) * np.cos(Om)
    z = yp * np.sin(i)
    tc = (t_ns - eph.toc) / ts.NS
    dts = eph.af0 + eph.af1 * tc + eph.af2 * tc ** 2 + (F_REL * eph.e * eph.sqrt_a * np.sin(E) if include_rel else 0.0)
    return np.array([x, y, z]), dts
