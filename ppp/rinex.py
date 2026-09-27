"""RINEX 3.x/4.x observation and navigation reading (TDS § 13.2).

Only the systems in settings.SYSTEMS are kept (GPS in v1).  Hatanaka (CRX) and gzip/Unix-compress
inputs are decompressed with the ``hatanaka`` package (TDS § 33.1).  Epoch tags are converted to
integer ns since the GPS epoch in GPST (TIME OF FIRST OBS time system must be GPS, TDS § 7.1).
"""
import gzip
import os
from dataclasses import dataclass, field

import numpy as np

from . import settings
from . import timesys as ts


class RinexError(Exception):
    pass


@dataclass
class ObsHeader:
    version: float = 0.0
    marker_name: str = ""
    marker_number: str = ""
    receiver: str = ""
    receiver_serial: str = ""
    antenna_type: str = ""          # full 20-character ANTEX code: type (16) + radome (4)
    antenna_serial: str = ""
    delta_hen: tuple = (0.0, 0.0, 0.0)   # ANTENNA: DELTA H/E/N (m)
    approx_xyz: tuple = (0.0, 0.0, 0.0)
    obs_types: dict = field(default_factory=dict)
    interval: float = 0.0
    time_system: str = "GPS"
    rcv_clock_offs_appl: int = 0
    antenna_zerodir: tuple = None
    phase_shifts: dict = field(default_factory=dict)
    run_by: str = ""
    raw: list = field(default_factory=list)


@dataclass
class ObsData:
    header: ObsHeader
    epochs: np.ndarray               # int64 ns GPST (receiver epoch tags)
    sats: list                       # e.g. ['G01', ...]
    obs: dict                        # code -> float array (n_epoch, n_sat), NaN = missing
    lli: dict                        # code -> int8 array (LLI flag, 0 if blank)
    snr_flag: dict                   # code -> int8 array (signal strength indicator 1-9)
    clock_offset: np.ndarray         # receiver clock offset from epoch record (s), NaN if absent
    epoch_flags: np.ndarray
    path: str = ""
    n_epochs_file: int = 0
    interval_file: float = 0.0


def _read_text(path):
    """Return file content as text; handles .crx/.d Hatanaka and .gz/.Z compression."""
    low = path.lower()
    with open(path, "rb") as fh:
        head = fh.read(4)
    compressed = head[:2] == b"\x1f\x8b" or head[:2] == b"\x1f\x9d"
    if compressed or low.endswith((".crx", ".d", ".crx.gz", ".d.gz", ".d.z", ".crx.z")):
        try:
            import hatanaka
            data = hatanaka.decompress(path)
            return data.decode("ascii", errors="replace") if isinstance(data, bytes) else data
        except ImportError:
            if head[:2] == b"\x1f\x8b":
                with gzip.open(path, "rb") as fh:
                    return fh.read().decode("ascii", errors="replace")
            raise RinexError("hatanaka package required for compressed RINEX: pip install hatanaka")
    with open(path, "rb") as fh:
        txt = fh.read().decode("ascii", errors="replace")
    if "COMPACT RINEX" in txt[:200]:
        import hatanaka
        return hatanaka.rnxdecompress(txt.encode()).decode("ascii", errors="replace")
    return txt


def _f(s, default=0.0):
    s = s.strip().replace("D", "E").replace("d", "e")
    if not s:
        return default
    return float(s)


def parse_obs_header(lines):
    h = ObsHeader()
    cur_sys = None
    n_end = None
    for i, line in enumerate(lines):
        label = line[60:80].strip()
        h.raw.append(line.rstrip())
        if label == "RINEX VERSION / TYPE":
            h.version = _f(line[0:9])
            if line[20:21] not in ("O",):
                raise RinexError("not an observation file (RINEX VERSION / TYPE)")
        elif label == "PGM / RUN BY / DATE":
            h.run_by = line[20:40].strip()
        elif label == "MARKER NAME":
            h.marker_name = line[0:60].strip()
        elif label == "MARKER NUMBER":
            h.marker_number = line[0:20].strip()
        elif label == "REC # / TYPE / VERS":
            h.receiver_serial = line[0:20].strip()
            h.receiver = line[20:40].strip()
        elif label == "ANT # / TYPE":
            h.antenna_serial = line[0:20].strip()
            ant = line[20:40]
            # 20-char ANTEX code: 16-char type + 4-char radome; blank radome -> NONE
            typ, rad = ant[0:16].rstrip(), ant[16:20].strip()
            h.antenna_type = f"{typ:<16s}{(rad or 'NONE'):>4s}"
        elif label == "APPROX POSITION XYZ":
            h.approx_xyz = (_f(line[0:14]), _f(line[14:28]), _f(line[28:42]))
        elif label == "ANTENNA: DELTA H/E/N":
            h.delta_hen = (_f(line[0:14]), _f(line[14:28]), _f(line[28:42]))
        elif label == "ANTENNA: ZERODIR AZI":
            h.antenna_zerodir = (_f(line[0:14]),)
        elif label == "SYS / # / OBS TYPES":
            if line[0] != " ":
                cur_sys = line[0]
                h.obs_types[cur_sys] = []
            codes = line[7:58].split()
            h.obs_types[cur_sys].extend(codes)
        elif label == "INTERVAL":
            h.interval = _f(line[0:10])
        elif label == "TIME OF FIRST OBS":
            h.time_system = line[48:51].strip() or "GPS"
        elif label == "RCV CLOCK OFFS APPL":
            try:
                h.rcv_clock_offs_appl = int(line[0:6])
            except ValueError:
                h.rcv_clock_offs_appl = 0
        elif label == "SYS / PHASE SHIFT":
            if line[0] != " ":
                try:
                    h.phase_shifts[(line[0], line[2:5])] = _f(line[6:14])
                except ValueError:
                    pass
        elif label == "END OF HEADER":
            n_end = i
            break
    if n_end is None:
        raise RinexError("END OF HEADER not found")
    if h.version < 3.0:
        raise RinexError(f"RINEX version {h.version} not supported (RINEX 3.02-4.0x required, TDS § 13.2)")
    return h, n_end + 1


def read_obs(path, systems=None, decimate_s=None):
    """Read a RINEX 3/4 observation file.

    decimate_s: keep only epochs whose GPST second-of-day is a multiple of this (TDS § 13.4).
    Returns ObsData with arrays for the requested systems only.
    """
    systems = tuple(systems or settings.SYSTEMS)
    txt = _read_text(path)
    lines = txt.splitlines()
    hdr, start = parse_obs_header(lines)
    if hdr.time_system not in ("GPS", ""):
        raise RinexError(f"time system {hdr.time_system} not supported (GPS time required)")
    codes = []
    for s in systems:
        for c in hdr.obs_types.get(s, []):
            if c not in codes:
                codes.append(c)
    # column index of each code per system
    col = {s: {c: k for k, c in enumerate(hdr.obs_types.get(s, []))} for s in systems}

    epoch_idx = [i for i in range(start, len(lines)) if lines[i].startswith(">")]
    n_all = len(epoch_idx)
    ep_ns = np.zeros(n_all, dtype=np.int64)
    keep = np.zeros(n_all, dtype=bool)
    flags = np.zeros(n_all, dtype=np.int8)
    clk = np.full(n_all, np.nan)
    for k, i in enumerate(epoch_idx):
        L = lines[i]
        try:
            y, mo, d, hh, mi = int(L[2:6]), int(L[7:9]), int(L[10:12]), int(L[13:15]), int(L[16:18])
            sec = float(L[18:29])
            fl = int(L[31:32] or 0)
        except ValueError:
            continue
        flags[k] = fl
        if fl > 1:
            continue
        ep_ns[k] = ts.to_ns(y, mo, d, hh, mi, sec)
        if len(L) > 41 and L[41:56].strip():
            try:
                clk[k] = float(L[41:56])
            except ValueError:
                pass
        keep[k] = True
    if decimate_s:
        dec = int(round(decimate_s * ts.NS))
        rem = ep_ns % dec
        keep &= (rem < 1_000_000) | (dec - rem < 1_000_000)     # within 1 ms of the grid
    # satellite list
    sats = set()
    for k, i in enumerate(epoch_idx):
        if not keep[k]:
            continue
        n = int(lines[i][32:35])
        for L in lines[i + 1:i + 1 + n]:
            if L[:1] in systems:
                sats.add(L[0:3].replace(" ", "0"))
    sats = sorted(sats)
    sidx = {s: j for j, s in enumerate(sats)}
    kept = np.nonzero(keep)[0]
    n_ep = len(kept)
    obs = {c: np.full((n_ep, len(sats)), np.nan) for c in codes}
    lli = {c: np.zeros((n_ep, len(sats)), dtype=np.int8) for c in codes}
    ssi = {c: np.zeros((n_ep, len(sats)), dtype=np.int8) for c in codes}
    for e, k in enumerate(kept):
        i = epoch_idx[k]
        n = int(lines[i][32:35])
        for L in lines[i + 1:i + 1 + n]:
            s = L[:1]
            if s not in systems:
                continue
            j = sidx[L[0:3].replace(" ", "0")]
            for c, kc in col[s].items():
                a = 3 + 16 * kc
                fld = L[a:a + 14]
                if fld.strip():
                    try:
                        obs[c][e, j] = float(fld)
                    except ValueError:
                        continue
                    fl = L[a + 14:a + 15]
                    if fl.strip():
                        lli[c][e, j] = int(fl)
                    fs = L[a + 15:a + 16]
                    if fs.strip():
                        ssi[c][e, j] = int(fs)
    ep = ep_ns[kept]
    # SNR 'S' observables are kept as obs too (code 'S1C' etc.)
    data = ObsData(header=hdr, epochs=ep, sats=sats, obs=obs, lli=lli, snr_flag=ssi,
                   clock_offset=clk[kept], epoch_flags=flags[kept], path=os.path.abspath(path),
                   n_epochs_file=int(keep.sum()) if not decimate_s else n_all)
    if n_ep > 1:
        data.interval_file = float(np.median(np.diff(ep_ns[keep] if not decimate_s else ep)) / ts.NS)
    if not hdr.interval and n_ep > 1:
        hdr.interval = float(np.median(np.diff(ep)) / ts.NS)
    return data


def first_last_epoch(path):
    """First and last data-record epochs (GPST ns) read from data records, not the file name (TDS § 12.1)."""
    txt = _read_text(path)
    lines = txt.splitlines()
    _, start = parse_obs_header(lines)
    first = last = None
    for L in lines[start:]:
        if L.startswith(">") and int(L[31:32] or 0) <= 1:
            t = ts.to_ns(int(L[2:6]), int(L[7:9]), int(L[10:12]), int(L[13:15]), int(L[16:18]), float(L[18:29]))
            first = t if first is None else first
            last = t
    return first, last


# ------------------------------------------------------------------------------------ navigation
@dataclass
class GpsEph:
    prn: str
    toc: int            # ns GPST
    af0: float
    af1: float
    af2: float
    iode: float
    crs: float
    delta_n: float
    m0: float
    cuc: float
    e: float
    cus: float
    sqrt_a: float
    toe_sow: float
    cic: float
    omega0: float
    cis: float
    i0: float
    crc: float
    omega: float
    omega_dot: float
    idot: float
    week: int
    sv_accuracy: float
    health: float
    tgd: float
    iodc: float
    toe: int = 0        # ns GPST


def read_nav(path):
    """Read GPS LNAV ephemerides from a RINEX 3.x/4.x navigation file. Returns dict prn -> [GpsEph]."""
    txt = _read_text(path)
    lines = txt.splitlines()
    i = 0
    while i < len(lines):
        if "END OF HEADER" in lines[i]:
            break
        i += 1
    i += 1
    out = {}
    while i < len(lines):
        L = lines[i]
        if L.startswith(">"):                          # RINEX 4 record header: > EPH G01 LNAV
            parts = L[1:].split()
            if len(parts) >= 3 and parts[0] == "EPH" and parts[1].startswith("G") and parts[2] == "LNAV":
                i += 1
                continue
            # skip other record types until next '>'
            i += 1
            while i < len(lines) and not lines[i].startswith(">"):
                i += 1
            continue
        if not L or L[0] != "G":
            i += 1
            continue
        blk = lines[i:i + 8]
        if len(blk) < 8:
            break
        try:
            prn = L[0:3].replace(" ", "0")
            y, mo, d, hh, mi, ss = int(L[4:8]), int(L[9:11]), int(L[12:14]), int(L[15:17]), int(L[18:20]), int(L[21:23])
            vals = [_f(L[23 + 19 * k:23 + 19 * (k + 1)]) for k in range(3)]
            for b in blk[1:]:
                vals += [_f(b[4 + 19 * k:4 + 19 * (k + 1)]) for k in range(4)]
            toc = ts.to_ns(y, mo, d, hh, mi, ss)
            eph = GpsEph(prn=prn, toc=toc, af0=vals[0], af1=vals[1], af2=vals[2],
                         iode=vals[3], crs=vals[4], delta_n=vals[5], m0=vals[6],
                         cuc=vals[7], e=vals[8], cus=vals[9], sqrt_a=vals[10],
                         toe_sow=vals[11], cic=vals[12], omega0=vals[13], cis=vals[14],
                         i0=vals[15], crc=vals[16], omega=vals[17], omega_dot=vals[18],
                         idot=vals[19], week=int(vals[21]), sv_accuracy=vals[23], health=vals[24],
                         tgd=vals[25], iodc=vals[26])
            eph.toe = eph.week * 7 * ts.DAY_NS + int(round(eph.toe_sow * ts.NS))
            # week in RINEX nav is continuous GPS week; guard against toe/toc week rollover mismatch
            if abs(eph.toe - toc) > 4 * ts.DAY_NS:
                wk = ts.gps_week_sow(toc)[0]
                cand = [(w * 7 * ts.DAY_NS + int(round(eph.toe_sow * ts.NS))) for w in (wk - 1, wk, wk + 1)]
                eph.toe = min(cand, key=lambda c: abs(c - toc))
            out.setdefault(prn, []).append(eph)
        except (ValueError, IndexError):
            pass
        i += 8
    for prn in out:
        out[prn].sort(key=lambda e: e.toe)
    return out
