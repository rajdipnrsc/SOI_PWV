"""ANTEX 1.4 reading and PCO/PCV evaluation (TDS § 7.6, § 7.7, § 8).

Conventions [F: ANTEX 1.4; checked against RTKLIB antmodel/satantpcv, see verify_log V-006]:
* receiver PCO in the file is (N, E, U) from the ARP, mm;  satellite PCO is (x, y, z) body frame from CoM, mm.
* receiver range correction   d_rho = -PCO_enu . e_enu + PCV(zenith, azimuth)   (e = unit receiver->satellite)
* satellite:  APC = CoM + R_body->ecef . PCO_xyz  (applied to the satellite position), plus PCV(nadir) added
  to the range.
* Ionosphere-free combination of PCO/PCV with alpha1, alpha2 of the selected frequencies.
"""
import os
from dataclasses import dataclass, field

import numpy as np

from . import settings
from . import timesys as ts
from .orbclk import open_text


@dataclass
class AntennaEntry:
    type: str                    # 20-char type+radome (receiver) or block (satellite)
    serial: str                  # receiver serial / satellite PRN code 'G01'
    svn: str = ""                # satellite SVN 'G063'
    cospar: str = ""
    valid_from: int = None       # ns GPST
    valid_until: int = None
    dazi: float = 0.0
    zen: tuple = (0.0, 90.0, 5.0)
    pco: dict = field(default_factory=dict)       # freq code 'G01' -> np.array(3) metres
    noazi: dict = field(default_factory=dict)     # freq -> array over zenith grid (m)
    azi: dict = field(default_factory=dict)       # freq -> (n_azi, n_zen) (m)
    sinex_code: str = ""

    @property
    def is_satellite(self):
        return bool(self.serial) and self.serial[0] in "GRECJS" and len(self.serial.strip()) == 3 \
            and self.serial[1:3].isdigit() and self.type.startswith(("BLOCK", "GLONASS", "GALILEO", "BEIDOU",
                                                                       "QZSS", "IRNSS"))


@dataclass
class Antex:
    name: str
    path: str
    satellites: dict            # prn -> [AntennaEntry]
    receivers: dict             # 20-char type -> AntennaEntry (generic, no serial)
    pcv_type: str = "A"
    refant: str = ""


def _parse_time(s):
    p = s.split()
    if len(p) < 6:
        return None
    y, m, d, hh, mi = (int(x) for x in p[:5])
    return ts.to_ns(y, m, d, hh, mi, float(p[5]))


def read_antex(path):
    txt = open_text(path)
    lines = txt.splitlines()
    sats, recs = {}, {}
    pcv_type, refant = "A", ""
    cur = None
    freq = None
    zgrid = None
    known = {"PCV TYPE / REFANT", "START OF ANTENNA", "TYPE / SERIAL NO", "DAZI", "ZEN1 / ZEN2 / DZEN",
             "VALID FROM", "VALID UNTIL", "SINEX CODE", "START OF FREQUENCY", "NORTH / EAST / UP",
             "END OF FREQUENCY", "END OF ANTENNA", "METH / BY / # / DATE", "# OF FREQUENCIES", "COMMENT",
             "START OF FREQ RMS", "END OF FREQ RMS"}
    in_rms = False
    for L in lines:
        lab = L[60:80].strip()
        if lab not in known:
            lab = ""
        if lab == "START OF FREQ RMS":
            in_rms = True
        elif lab == "END OF FREQ RMS":
            in_rms = False
            continue
        if in_rms:
            continue
        if lab == "PCV TYPE / REFANT":
            pcv_type = L[0:1]
            refant = L[20:40].strip()
        elif lab == "START OF ANTENNA":
            cur = AntennaEntry(type="", serial="")
        elif lab == "TYPE / SERIAL NO" and cur is not None:
            cur.type = L[0:20]
            cur.serial = L[20:40].strip()
            cur.svn = L[40:50].strip()
            cur.cospar = L[50:60].strip()
        elif lab == "DAZI" and cur is not None:
            cur.dazi = float(L[2:8])
        elif lab == "ZEN1 / ZEN2 / DZEN" and cur is not None:
            cur.zen = (float(L[2:8]), float(L[8:14]), float(L[14:20]))
            z1, z2, dz = cur.zen
            zgrid = np.arange(z1, z2 + dz / 2, dz)
        elif lab == "VALID FROM" and cur is not None:
            cur.valid_from = _parse_time(L[0:43])
        elif lab == "VALID UNTIL" and cur is not None:
            cur.valid_until = _parse_time(L[0:43])
        elif lab == "SINEX CODE" and cur is not None:
            cur.sinex_code = L[0:10].strip()
        elif lab == "START OF FREQUENCY" and cur is not None:
            freq = L[3:6].strip()
        elif lab == "NORTH / EAST / UP" and cur is not None and freq:
            cur.pco[freq] = np.array([float(L[0:10]), float(L[10:20]), float(L[20:30])]) * 1e-3
        elif lab == "END OF FREQUENCY":
            freq = None
        elif lab == "END OF ANTENNA" and cur is not None:
            if cur.type.startswith(("BLOCK", "GLONASS", "GALILEO", "BEIDOU", "QZSS", "IRNSS")) or \
                    (len(cur.serial) == 3 and cur.serial[0] in "GRECJI" and cur.serial[1:].isdigit()):
                sats.setdefault(cur.serial, []).append(cur)
            else:
                if not cur.serial:          # generic (type-mean) calibration
                    recs[cur.type] = cur
                else:
                    recs.setdefault(("SERIAL", cur.type, cur.serial), cur)
            cur = None
        elif cur is not None and freq and zgrid is not None and lab == "" and L.strip():
            head = L[0:8].strip()
            vals = L[8:].split()
            try:
                arr = np.array([float(v) for v in vals]) * 1e-3
            except ValueError:
                continue
            if head == "NOAZI":
                cur.noazi[freq] = arr
            else:
                try:
                    az = float(head)
                except ValueError:
                    continue
                cur.azi.setdefault(freq, []).append((az, arr))
    for group in list(sats.values()) + [recs]:
        entries = group if isinstance(group, list) else group.values()
        for e in entries:
            for f, rows in list(e.azi.items()):
                if isinstance(rows, list):
                    rows.sort(key=lambda r: r[0])
                    e.azi[f] = np.array([r[1] for r in rows])
    name = os.path.splitext(os.path.basename(path))[0]
    return Antex(name, path, sats, recs, pcv_type, refant)


# ---------------------------------------------------------------------------------- lookups
def satellite_entry(atx, prn, t_ns):
    """ANTEX satellite entry valid at t for PRN (mapping PRN -> SVN by validity period, TDS § 7.6)."""
    for e in atx.satellites.get(prn, []):
        a = e.valid_from if e.valid_from is not None else -2**62
        b = e.valid_until if e.valid_until is not None else 2**62
        if a <= t_ns < b:
            return e
    return None


def receiver_entry(atx, ant_type, serial=None):
    """Exact 20-character type+radome match; radome fallback to NONE (flag).  TDS § 8, § 9.11.

    Returns (entry or None, status) with status in {'EXACT', 'RADOME_FALLBACK', 'NO_ANTENNA_CALIBRATION'}.
    Individual calibrations (by serial) are used only if present, flagged 'INDIVIDUAL'.
    """
    key = f"{ant_type[:16]:<16s}{ant_type[16:20].strip() or 'NONE':>4s}"
    if serial and ("SERIAL", key, serial) in atx.receivers:
        return atx.receivers[("SERIAL", key, serial)], "INDIVIDUAL"
    if key in atx.receivers:
        return atx.receivers[key], "EXACT"
    none_key = f"{ant_type[:16]:<16s}{'NONE':>4s}"
    if none_key in atx.receivers:
        return atx.receivers[none_key], "RADOME_FALLBACK"
    return None, "NO_ANTENNA_CALIBRATION"


def freq_code(system, band):
    return f"{system}{int(band):02d}"


def if_coefficients(system="G", bands=("1", "2")):
    f1 = settings.FREQ[(system, bands[0])]
    f2 = settings.FREQ[(system, bands[1])]
    a1 = f1 ** 2 / (f1 ** 2 - f2 ** 2)
    a2 = f2 ** 2 / (f1 ** 2 - f2 ** 2)
    return a1, a2


def _pcv_one(entry, freq, zen_deg, azi_deg):
    """Bilinear (azimuth, zenith/nadir) interpolation of one frequency's PCV table (m)."""
    z1, z2, dz = entry.zen
    nz = int(round((z2 - z1) / dz)) + 1
    scalar = np.ndim(zen_deg) == 0
    z = np.clip(np.atleast_1d(np.asarray(zen_deg, dtype=float)), z1, z2)
    fz_all = (z - z1) / dz
    iz = np.clip(np.floor(fz_all).astype(int), 0, nz - 2)
    fz = fz_all - iz
    tab = entry.azi.get(freq)
    if isinstance(tab, np.ndarray) and entry.dazi > 0:
        naz = tab.shape[0]
        az = np.mod(np.broadcast_to(np.atleast_1d(np.asarray(azi_deg, dtype=float)), z.shape), 360.0)
        fa_all = az / entry.dazi
        ia = np.clip(np.floor(fa_all).astype(int), 0, naz - 2)
        fa = fa_all - ia
        v0 = tab[ia, iz] * (1 - fz) + tab[ia, iz + 1] * fz
        v1 = tab[ia + 1, iz] * (1 - fz) + tab[ia + 1, iz + 1] * fz
        out = (1 - fa) * v0 + fa * v1
    elif freq in entry.noazi:
        row = entry.noazi[freq]
        out = row[iz] * (1 - fz) + row[iz + 1] * fz
    else:
        out = np.zeros_like(z)
    return float(out[0]) if scalar else out


def pco_if(entry, system="G", bands=("1", "2")):
    """IF-combined PCO (m): alpha1*PCO1 - alpha2*PCO2 (TDS § 7.6)."""
    a1, a2 = if_coefficients(system, bands)
    f1, f2 = freq_code(system, bands[0]), freq_code(system, bands[1])
    if f1 not in entry.pco or f2 not in entry.pco:
        return None
    return a1 * entry.pco[f1] - a2 * entry.pco[f2]


def pcv_if(entry, zen_deg, azi_deg=0.0, system="G", bands=("1", "2")):
    """IF-combined PCV (m) at zenith (receiver) or nadir (satellite) angle and azimuth (deg)."""
    a1, a2 = if_coefficients(system, bands)
    f1, f2 = freq_code(system, bands[0]), freq_code(system, bands[1])
    return a1 * _pcv_one(entry, f1, zen_deg, azi_deg) - a2 * _pcv_one(entry, f2, zen_deg, azi_deg)


def receiver_correction(entry, e_enu, system="G", bands=("1", "2")):
    """Receiver antenna range correction (m): -PCO_IF(enu).e + PCV_IF(z, a)  (TDS § 7.7).

    e_enu: (n,3) unit vectors receiver->satellite in local E, N, U.  ANTEX gives (N, E, U).
    """
    if entry is None:
        return np.zeros(len(e_enu))
    p = pco_if(entry, system, bands)
    pco_enu = np.array([p[1], p[0], p[2]])
    zen = np.degrees(np.arccos(np.clip(e_enu[:, 2], -1, 1)))
    azi = np.degrees(np.arctan2(e_enu[:, 0], e_enu[:, 1])) % 360.0
    return -(e_enu @ pco_enu) + pcv_if(entry, zen, azi, system, bands)


def block_of(entry):
    return entry.type.strip() if entry is not None else ""
