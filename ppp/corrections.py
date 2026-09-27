"""Deterministic corrections (TDS § 7). Every correction is a small, separately testable function.

Sign convention (TDS § 7): a positive range correction increases the modelled observation.
Station displacements (solid tide, ocean loading, pole tide) are returned as ECEF vectors (m) that are
added to the station position before the geometric range is computed.
"""
import os

import numpy as np

from . import settings
from . import timesys as ts

C = settings.C_LIGHT
D2R = np.pi / 180.0
AS2R = D2R / 3600.0


# ============================================================================== Sun and Moon
def _load_moon_tables():
    lr, b, cur = [], [], None
    with open(os.path.join(settings.DATA_DIR, "moon_meeus47.txt")) as fh:
        for L in fh:
            if L.startswith("#") or not L.strip():
                continue
            if L.startswith("@"):
                cur = L.split()[1]
                continue
            v = [float(x) for x in L.split()]
            (lr if cur == "LR" else b).append(v)
    return np.array(lr), np.array(b)


_MOON_LR, _MOON_B = _load_moon_tables()


def moon_ecliptic(jde):
    """Geocentric ecliptic longitude, latitude (rad, mean equinox of date) and distance (m) of the Moon.

    Meeus (1998) Ch. 47 (truncated ELP-2000/82). Accuracy ~10" in longitude, 4" in latitude [F: Meeus].
    jde: Julian ephemeris date (TT), scalar or array.
    """
    t = (np.asarray(jde, dtype=float) - 2451545.0) / 36525.0
    Lp = 218.3164477 + (481267.88123421 + (-0.0015786 + (1.0 / 538841.0 - t / 65194000.0) * t) * t) * t
    D = 297.8501921 + (445267.1114034 + (-0.0018819 + (1.0 / 545868.0 - t / 113065000.0) * t) * t) * t
    M = 357.5291092 + (35999.0502909 + (-0.0001536 + t / 24490000.0) * t) * t
    Mp = 134.9633964 + (477198.8675055 + (0.0087414 + (1.0 / 69699.9 + t / 14712000.0) * t) * t) * t
    F = 93.2720950 + (483202.0175233 + (-0.0036539 + (-1.0 / 3526000.0 + t / 863310000.0) * t) * t) * t
    A1 = 119.75 + 131.849 * t
    A2 = 53.09 + 479264.290 * t
    A3 = 313.45 + 481266.484 * t
    E = 1.0 + (-0.002516 - 0.0000074 * t) * t
    args = np.stack([np.mod(D, 360), np.mod(M, 360), np.mod(Mp, 360), np.mod(F, 360)]) * D2R   # (4, n)
    Er = E.reshape(1, -1)
    arg_lr = _MOON_LR[:, :4] @ args.reshape(4, -1)
    ecorr_lr = Er ** np.abs(_MOON_LR[:, 1])[:, None]
    sl = np.sum(_MOON_LR[:, 4:5] * ecorr_lr * np.sin(arg_lr), axis=0)
    sr = np.sum(_MOON_LR[:, 5:6] * ecorr_lr * np.cos(arg_lr), axis=0)
    arg_b = _MOON_B[:, :4] @ args.reshape(4, -1)
    ecorr_b = Er ** np.abs(_MOON_B[:, 1])[:, None]
    sb = np.sum(_MOON_B[:, 4:5] * ecorr_b * np.sin(arg_b), axis=0)
    Lpr, Mpr, Fr = (np.mod(Lp, 360) * D2R).ravel(), (np.mod(Mp, 360) * D2R).ravel(), (np.mod(F, 360) * D2R).ravel()
    sl += 3958.0 * np.sin(np.mod(A1, 360).ravel() * D2R) + 1962.0 * np.sin(Lpr - Fr) + \
        318.0 * np.sin(np.mod(A2, 360).ravel() * D2R)
    sb += -2235.0 * np.sin(Lpr) + 382.0 * np.sin(np.mod(A3, 360).ravel() * D2R) + \
        175.0 * np.sin(np.mod(A1, 360).ravel() * D2R - Fr) + 175.0 * np.sin(np.mod(A1, 360).ravel() * D2R + Fr) + \
        127.0 * np.sin(Lpr - Mpr) - 115.0 * np.sin(Lpr + Mpr)
    lam = (np.mod(Lp, 360).ravel() + sl / 1e6) * D2R
    beta = sb / 1e6 * D2R
    dist = (385000.56 + sr / 1000.0) * 1000.0
    shape = np.shape(jde)
    return lam.reshape(shape), beta.reshape(shape), dist.reshape(shape)


def sun_ecliptic(jde):
    """Geometric ecliptic longitude (rad, mean equinox of date) and distance (m) of the Sun.

    Meeus (1998) Ch. 25 (accuracy 0.01 deg) [F: Meeus]. Latitude taken as 0.
    """
    t = (np.asarray(jde, dtype=float) - 2451545.0) / 36525.0
    L0 = 280.46646 + 36000.76983 * t + 0.0003032 * t * t
    M = (357.52911 + 35999.05029 * t - 0.0001537 * t * t) * D2R
    e = 0.016708634 - 0.000042037 * t - 0.0000001267 * t * t
    Cc = (1.914602 - 0.004817 * t - 0.000014 * t * t) * np.sin(M) + (0.019993 - 0.000101 * t) * np.sin(2 * M) + \
        0.000289 * np.sin(3 * M)
    lon = np.mod(L0 + Cc, 360.0) * D2R
    nu = M + Cc * D2R
    R = 1.000001018 * (1 - e * e) / (1 + e * np.cos(nu)) * settings.AU
    return lon, R


def mean_obliquity(jde):
    t = (np.asarray(jde, dtype=float) - 2451545.0) / 36525.0
    return (23.0 + 26.0 / 60 + 21.448 / 3600 + (-46.8150 * t - 0.00059 * t * t + 0.001813 * t ** 3) / 3600) * D2R


def gmst(jd_ut1):
    """Greenwich mean sidereal time (rad), Meeus eq. 12.4 (IAU 1982)."""
    t = (np.asarray(jd_ut1, dtype=float) - 2451545.0) / 36525.0
    th = 280.46061837 + 360.98564736629 * (np.asarray(jd_ut1, dtype=float) - 2451545.0) + \
        0.000387933 * t * t - t ** 3 / 38710000.0
    return np.mod(th, 360.0) * D2R


def _ecl_to_ecef(lam, beta, dist, eps, theta):
    cb = np.cos(beta)
    x = dist * cb * np.cos(lam)
    y = dist * (cb * np.sin(lam) * np.cos(eps) - np.sin(beta) * np.sin(eps))
    z = dist * (cb * np.sin(lam) * np.sin(eps) + np.sin(beta) * np.cos(eps))
    xe = np.cos(theta) * x + np.sin(theta) * y
    ye = -np.sin(theta) * x + np.cos(theta) * y
    return np.stack([xe, ye, z], axis=-1)


def sun_moon_ecef(t_gpst_ns, ut1_minus_utc=0.0):
    """Earth-fixed Sun and Moon positions (m) at GPST instants (array).

    TT from GPST (+51.184 s); UT1 = UTC + (UT1-UTC) from ERP; rotation by GMST about the mean-of-date
    equator; polar motion (< 1e-6 rad) neglected.  Error in solid-tide displacement < 0.05 mm
    (Level-0 test against analytic references, docs/verify_log.md V-015).
    """
    t = np.atleast_1d(np.asarray(t_gpst_ns, dtype=np.int64))
    jde = ts.gpst_to_tt_mjd(t) + 2400000.5
    utc = np.array([ts.gpst_to_utc(int(x)) for x in t])
    jd_ut1 = ts.mjd(utc) + 2400000.5 + np.asarray(ut1_minus_utc) / 86400.0
    eps = mean_obliquity(jde)
    th = gmst(jd_ut1)
    lm, bm, dm = moon_ecliptic(jde)
    ls, ds = sun_ecliptic(jde)
    moon = _ecl_to_ecef(lm, bm, dm, eps, th)
    sun = _ecl_to_ecef(ls, np.zeros_like(ls), ds, eps, th)
    return sun, moon


# ============================================================================== solid Earth tide
def solid_tide(xsta, xsun, xmon, t_utc_ns):
    """Solid Earth tide displacement (m, ECEF), IERS Conventions 2010 § 7.1.1 (DEHANTTIDEINEL port).

    Degree 2 and 3 with latitude-dependent Love/Shida numbers (Step 1), out-of-phase and l^(1) corrections,
    frequency-dependent Step-2 corrections (diurnal + long period). The permanent tide is included, so
    applying it to conventional tide-free coordinates is correct (TDS § 7.10).
    xsta (3,), xsun/xmon (n,3), t_utc_ns (n,) -> (n,3).
    """
    X = np.asarray(xsta, dtype=float)
    xsun = np.atleast_2d(xsun)
    xmon = np.atleast_2d(xmon)
    rsta = np.linalg.norm(X)
    rsun = np.linalg.norm(xsun, axis=1)
    rmon = np.linalg.norm(xmon, axis=1)
    scsun = xsun @ X / rsta / rsun
    scmon = xmon @ X / rsta / rmon
    cosphi = np.hypot(X[0], X[1]) / rsta
    h20, l20, h3, l3 = 0.6078, 0.0847, 0.292, 0.015
    h2 = h20 - 0.0006 * (1 - 1.5 * cosphi ** 2)
    l2 = l20 + 0.0002 * (1 - 1.5 * cosphi ** 2)
    p2sun = 3 * (h2 / 2 - l2) * scsun ** 2 - h2 / 2
    p2mon = 3 * (h2 / 2 - l2) * scmon ** 2 - h2 / 2
    p3sun = 2.5 * (h3 - 3 * l3) * scsun ** 3 + 1.5 * (l3 - h3) * scsun
    p3mon = 2.5 * (h3 - 3 * l3) * scmon ** 3 + 1.5 * (l3 - h3) * scmon
    x2sun, x2mon = 3 * l2 * scsun, 3 * l2 * scmon
    x3sun = 1.5 * l3 * (5 * scsun ** 2 - 1)
    x3mon = 1.5 * l3 * (5 * scmon ** 2 - 1)
    mass_sun, mass_mon, re = 332946.0482, 0.0123000371, 6378136.6
    fac2sun = mass_sun * re * (re / rsun) ** 3
    fac2mon = mass_mon * re * (re / rmon) ** 3
    fac3sun = fac2sun * (re / rsun)
    fac3mon = fac2mon * (re / rmon)
    us = X / rsta
    d = (fac2sun[:, None] * (x2sun[:, None] * xsun / rsun[:, None] + p2sun[:, None] * us)
         + fac2mon[:, None] * (x2mon[:, None] * xmon / rmon[:, None] + p2mon[:, None] * us)
         + fac3sun[:, None] * (x3sun[:, None] * xsun / rsun[:, None] + p3sun[:, None] * us)
         + fac3mon[:, None] * (x3mon[:, None] * xmon / rmon[:, None] + p3mon[:, None] * us))
    d = d + _st1idiu(X, xsun, xmon, fac2sun, fac2mon) + _st1isem(X, xsun, xmon, fac2sun, fac2mon) + \
        _st1l1(X, xsun, xmon, fac2sun, fac2mon)
    t_utc_ns = np.atleast_1d(np.asarray(t_utc_ns, dtype=np.int64))
    mjd_utc = ts.mjd(t_utc_ns)
    fhr = (t_utc_ns % ts.DAY_NS) / ts.NS / 3600.0
    dtt = np.array([ts.tai_minus_utc(m) for m in mjd_utc]) + 32.184
    T = (mjd_utc - 51544.5 + dtt / 86400.0) / 36525.0
    d = d + _step2diu(X, fhr, T) + _step2lon(X, T)
    return d


def _geo(X):
    rsta = np.linalg.norm(X)
    sinphi = X[2] / rsta
    cosphi = np.hypot(X[0], X[1]) / rsta
    sinla = X[1] / cosphi / rsta
    cosla = X[0] / cosphi / rsta
    return rsta, sinphi, cosphi, sinla, cosla


def _to_xyz(dr, dn, de, sinphi, cosphi, sinla, cosla):
    return np.stack([dr * cosla * cosphi - de * sinla - dn * sinphi * cosla,
                     dr * sinla * cosphi + de * cosla - dn * sinphi * sinla,
                     dr * sinphi + dn * cosphi], axis=-1)


def _st1idiu(X, xs, xm, f2s, f2m):
    _, sp, cp, sl, cl = _geo(X)
    rs, rm = np.linalg.norm(xs, axis=1), np.linalg.norm(xm, axis=1)
    dhi, dli = -0.0025, -0.0007
    cos2 = cp * cp - sp * sp
    ks = f2s * xs[:, 2] * (xs[:, 0] * sl - xs[:, 1] * cl) / rs ** 2
    km = f2m * xm[:, 2] * (xm[:, 0] * sl - xm[:, 1] * cl) / rm ** 2
    js = f2s * xs[:, 2] * (xs[:, 0] * cl + xs[:, 1] * sl) / rs ** 2
    jm = f2m * xm[:, 2] * (xm[:, 0] * cl + xm[:, 1] * sl) / rm ** 2
    dr = -3 * dhi * sp * cp * (ks + km)
    dn = -3 * dli * cos2 * (ks + km)
    de = -3 * dli * sp * (js + jm)
    return _to_xyz(dr, dn, de, sp, cp, sl, cl)


def _st1isem(X, xs, xm, f2s, f2m):
    _, sp, cp, sl, cl = _geo(X)
    rs, rm = np.linalg.norm(xs, axis=1), np.linalg.norm(xm, axis=1)
    dhi, dli = -0.0022, -0.0007
    c2l, s2l = cl * cl - sl * sl, 2 * cl * sl
    as_ = f2s * ((xs[:, 0] ** 2 - xs[:, 1] ** 2) * s2l - 2 * xs[:, 0] * xs[:, 1] * c2l) / rs ** 2
    am = f2m * ((xm[:, 0] ** 2 - xm[:, 1] ** 2) * s2l - 2 * xm[:, 0] * xm[:, 1] * c2l) / rm ** 2
    bs = f2s * ((xs[:, 0] ** 2 - xs[:, 1] ** 2) * c2l + 2 * xs[:, 0] * xs[:, 1] * s2l) / rs ** 2
    bm = f2m * ((xm[:, 0] ** 2 - xm[:, 1] ** 2) * c2l + 2 * xm[:, 0] * xm[:, 1] * s2l) / rm ** 2
    dr = -0.75 * dhi * cp ** 2 * (as_ + am)
    dn = 1.5 * dli * sp * cp * (as_ + am)
    de = -1.5 * dli * cp * (bs + bm)
    return _to_xyz(dr, dn, de, sp, cp, sl, cl)


def _st1l1(X, xs, xm, f2s, f2m):
    _, sp, cp, sl, cl = _geo(X)
    rs, rm = np.linalg.norm(xs, axis=1), np.linalg.norm(xm, axis=1)
    l1d, l1sd = 0.0012, 0.0024
    dns = -l1d * sp ** 2 * f2s * xs[:, 2] * (xs[:, 0] * cl + xs[:, 1] * sl) / rs ** 2
    dnm = -l1d * sp ** 2 * f2m * xm[:, 2] * (xm[:, 0] * cl + xm[:, 1] * sl) / rm ** 2
    des = l1d * sp * (cp ** 2 - sp ** 2) * f2s * xs[:, 2] * (xs[:, 0] * sl - xs[:, 1] * cl) / rs ** 2
    dem = l1d * sp * (cp ** 2 - sp ** 2) * f2m * xm[:, 2] * (xm[:, 0] * sl - xm[:, 1] * cl) / rm ** 2
    de, dn = 3 * (des + dem), 3 * (dns + dnm)
    out = _to_xyz(np.zeros_like(de), dn, de, sp, cp, sl, cl)
    c2l, s2l = cl ** 2 - sl ** 2, 2 * cl * sl
    dns = -l1sd / 2 * sp * cp * f2s * ((xs[:, 0] ** 2 - xs[:, 1] ** 2) * c2l + 2 * xs[:, 0] * xs[:, 1] * s2l) / rs ** 2
    dnm = -l1sd / 2 * sp * cp * f2m * ((xm[:, 0] ** 2 - xm[:, 1] ** 2) * c2l + 2 * xm[:, 0] * xm[:, 1] * s2l) / rm ** 2
    des = -l1sd / 2 * sp ** 2 * cp * f2s * ((xs[:, 0] ** 2 - xs[:, 1] ** 2) * s2l - 2 * xs[:, 0] * xs[:, 1] * c2l) / rs ** 2
    dem = -l1sd / 2 * sp ** 2 * cp * f2m * ((xm[:, 0] ** 2 - xm[:, 1] ** 2) * s2l - 2 * xm[:, 0] * xm[:, 1] * c2l) / rm ** 2
    de, dn = 3 * (des + dem), 3 * (dns + dnm)
    return out + _to_xyz(np.zeros_like(de), dn, de, sp, cp, sl, cl)


_DATDI_DIU = np.array([
    [-3, 0, 2, 0, 0, -0.01, 0, 0, 0], [-3, 2, 0, 0, 0, -0.01, 0, 0, 0], [-2, 0, 1, -1, 0, -0.02, 0, 0, 0],
    [-2, 0, 1, 0, 0, -0.08, 0, -0.01, 0.01], [-2, 2, -1, 0, 0, -0.02, 0, 0, 0], [-1, 0, 0, -1, 0, -0.10, 0, 0, 0],
    [-1, 0, 0, 0, 0, -0.51, 0, -0.02, 0.03], [-1, 2, 0, 0, 0, 0.01, 0, 0, 0], [0, -2, 1, 0, 0, 0.01, 0, 0, 0],
    [0, 0, -1, 0, 0, 0.02, 0, 0, 0], [0, 0, 1, 0, 0, 0.06, 0, 0, 0], [0, 0, 1, 1, 0, 0.01, 0, 0, 0],
    [0, 2, -1, 0, 0, 0.01, 0, 0, 0], [1, -3, 0, 0, 1, -0.06, 0, 0, 0], [1, -2, 0, -1, 0, 0.01, 0, 0, 0],
    [1, -2, 0, 0, 0, -1.23, -0.07, 0.06, 0.01], [1, -1, 0, 0, -1, 0.02, 0, 0, 0], [1, -1, 0, 0, 1, 0.04, 0, 0, 0],
    [1, 0, 0, -1, 0, -0.22, 0.01, 0.01, 0], [1, 0, 0, 0, 0, 12.00, -0.80, -0.67, -0.03],
    [1, 0, 0, 1, 0, 1.73, -0.12, -0.10, 0], [1, 0, 0, 2, 0, -0.04, 0, 0, 0], [1, 1, 0, 0, -1, -0.50, -0.01, 0.03, 0],
    [1, 1, 0, 0, 1, 0.01, 0, 0, 0], [0, 1, 0, 1, -1, -0.01, 0, 0, 0], [1, 2, -2, 0, 0, -0.01, 0, 0, 0],
    [1, 2, 0, 0, 0, -0.11, 0.01, 0.01, 0], [2, -2, 1, 0, 0, -0.01, 0, 0, 0], [2, 0, -1, 0, 0, -0.02, 0, 0, 0],
    [3, 0, 0, 0, 0, 0, 0, 0, 0], [3, 0, 0, 1, 0, 0, 0, 0, 0]])
_DATDI_LON = np.array([
    [0, 0, 0, 1, 0, 0.47, 0.23, 0.16, 0.07], [0, 2, 0, 0, 0, -0.20, -0.12, -0.11, -0.05],
    [1, 0, -1, 0, 0, -0.11, -0.08, -0.09, -0.04], [2, 0, 0, 0, 0, -0.13, -0.11, -0.15, -0.07],
    [2, 0, 0, 1, 0, -0.05, -0.05, -0.06, -0.03]])


def _fund_args(T):
    s = 218.31664563 + (481267.88194 + (-0.0014663889 + 0.00000185139 * T) * T) * T
    pr = (1.396971278 + (0.000308889 + (0.000000021 + 0.000000007 * T) * T) * T) * T
    h = 280.46645 + (36000.7697489 + (0.00030322222 + (0.000000020 - 0.00000000654 * T) * T) * T) * T
    p = 83.35324312 + (4069.01363525 + (-0.01032172222 + (-0.0000124991 + 0.00000005263 * T) * T) * T) * T
    zns = 234.95544499 + (1934.13626197 + (-0.00207561111 + (-0.00000213944 + 0.00000001650 * T) * T) * T) * T
    ps = 282.93734098 + (1.71945766667 + (0.00045688889 + (-0.00000001778 - 0.00000000334 * T) * T) * T) * T
    return s, pr, h, p, zns, ps


def _step2diu(X, fhr, T):
    _, sp, cp, sl, cl = _geo(X)
    s, pr, h, p, zns, ps = _fund_args(T)
    tau = fhr * 15 + 280.4606184 + (36000.7700536 + (0.00038793 - 0.0000000258 * T) * T) * T - s
    s = np.mod(s + pr, 360)
    tau, h, p, zns, ps = (np.mod(v, 360) for v in (tau, h, p, zns, ps))
    zla = np.arctan2(X[1], X[0])
    out = np.zeros((len(np.atleast_1d(T)), 3))
    for row in _DATDI_DIU:
        th = (tau + row[0] * s + row[1] * h + row[2] * p + row[3] * zns + row[4] * ps) * D2R
        dr = row[5] * 2 * sp * cp * np.sin(th + zla) + row[6] * 2 * sp * cp * np.cos(th + zla)
        dn = row[7] * (cp ** 2 - sp ** 2) * np.sin(th + zla) + row[8] * (cp ** 2 - sp ** 2) * np.cos(th + zla)
        de = row[7] * sp * np.cos(th + zla) - row[8] * sp * np.sin(th + zla)
        out += _to_xyz(dr, dn, de, sp, cp, sl, cl)
    return out / 1000.0


def _step2lon(X, T):
    _, sp, cp, sl, cl = _geo(X)
    s, pr, h, p, zns, ps = _fund_args(T)
    s = np.mod(s + pr, 360)
    h, p, zns, ps = (np.mod(v, 360) for v in (h, p, zns, ps))
    out = np.zeros((len(np.atleast_1d(T)), 3))
    for row in _DATDI_LON:
        th = (row[0] * s + row[1] * h + row[2] * p + row[3] * zns + row[4] * ps) * D2R
        dr = row[5] * (3 * sp ** 2 - 1) / 2 * np.cos(th) + row[7] * (3 * sp ** 2 - 1) / 2 * np.sin(th)
        dn = row[6] * (cp * sp * 2) * np.cos(th) + row[8] * (cp * sp * 2) * np.sin(th)
        out += _to_xyz(dr, dn, np.zeros_like(dr), sp, cp, sl, cl)
    return out / 1000.0


# ============================================================================== ocean tide loading
def _load_hardisp():
    rows = np.loadtxt(os.path.join(settings.DATA_DIR, "hardisp_constituents.txt"), comments="#")
    return rows[:, :6].astype(int), rows[:, 6]


_HD_IDD, _HD_TAMP = _load_hardisp()
# the 11 BLQ constituents M2 S2 N2 K2 K1 O1 P1 Q1 MF MM SSA (HARDISP IDT) [F]
_BLQ_IDT = np.array([[2, 0, 0, 0, 0, 0], [2, 2, -2, 0, 0, 0], [2, -1, 0, 1, 0, 0], [2, 2, 0, 0, 0, 0],
                     [1, 1, 0, 0, 0, 0], [1, -1, 0, 0, 0, 0], [1, 1, -2, 0, 0, 0], [1, -2, 0, 1, 0, 0],
                     [0, 2, 0, 0, 0, 0], [0, 1, 0, -1, 0, 0], [0, 0, 2, 0, 0, 0]])


def read_blq(path, station):
    """Read the BLQ block of a station (first 4 characters matched, case-insensitive).

    Returns (amp (3,11) m, phase (3,11) deg, header lines) with rows radial, tangential W, tangential S
    [F: BLQ format; TDS § 7.12], or None if not found.
    """
    with open(path) as fh:
        lines = [L.rstrip("\n") for L in fh]
    st = station[:4].upper()
    i = 0
    while i < len(lines):
        L = lines[i]
        if not L.startswith("$$") and L.strip() and L.strip()[:4].upper() == st:
            hdr, nums = [L], []
            j = i + 1
            while j < len(lines) and len(nums) < 6:
                if lines[j].startswith("$$"):
                    hdr.append(lines[j])
                elif lines[j].strip():
                    nums.append([float(v) for v in lines[j].split()])
                j += 1
            if len(nums) == 6 and all(len(n) == 11 for n in nums):
                return np.array(nums[:3]), np.array(nums[3:]), hdr
        i += 1
    return None


def _etutc(year):
    """ET-UTC (s) as in HARDISP ETUTC (only the >= 1972 branch is needed here)."""
    st = [1972.5, 1973.0, 1974.0, 1975.0, 1976.0, 1977.0, 1978.0, 1979.0, 1980.0, 1981.5, 1982.5, 1983.5,
          1985.5, 1988.0, 1990.0, 1991.0, 1992.5, 1993.5, 1994.5, 1996.0, 1997.5, 1999.0, 2006.0, 2009.0,
          2012.5, 2015.5, 2017.0]
    return 42.184 + sum(1.0 for s in st if year >= s)


def _tdfrph(idd, t_utc_ns):
    """Frequency (cycles/day) and phase (deg) of constituents (HARDISP TDFRPH). idd (k,6)."""
    y, m, d, hh, mi, sec = ts.from_ns(int(t_utc_ns))
    doy = ts.mjd_from_date(y, m, d) - ts.mjd_from_date(y, 1, 1) + 1
    dayfr = hh / 24.0 + mi / 1440.0 + int(sec) / 86400.0
    leap = 1 if (y % 4 == 0 and (y % 100 != 0 or y % 400 == 0)) else 0
    year = y + (doy + dayfr) / (365.0 + leap)
    delta = _etutc(year)
    jd = ts.mjd_from_date(y, 1, 1) + 2400000.5 + (doy - 1) + 0.5     # JULDAT (noon JD of the day)
    djd = jd - 0.5 + dayfr
    T = (djd - 2451545.0 + delta / 86400.0) / 36525.0
    F1 = 134.9634025100 + T * (477198.8675605000 + T * (0.0088553333 + T * (0.0000143431 + T * (-0.0000000680))))
    F2 = 357.5291091806 + T * (35999.0502911389 + T * (-0.0001536667 + T * (0.0000000378 + T * (-0.0000000032))))
    F3 = 93.2720906200 + T * (483202.0174577222 + T * (-0.0035420000 + T * (-0.0000002881 + T * 0.0000000012)))
    F4 = 297.8501954694 + T * (445267.1114469445 + T * (-0.0017696111 + T * (0.0000018314 + T * (-0.0000000088))))
    F5 = 125.0445550100 + T * (-1934.1362619722 + T * (0.0020756111 + T * (0.0000021394 + T * (-0.0000000165))))
    Dv = np.array([360.0 * dayfr - F4, F3 + F5, F3 + F5 - F4, F3 + F5 - F1, -F5, F3 + F5 - F4 - F2])
    FD1 = 0.0362916471 + 0.0000000013 * T
    FD2 = 0.0027377786
    FD3 = 0.0367481951 - 0.0000000005 * T
    FD4 = 0.0338631920 - 0.0000000003 * T
    FD5 = -0.0001470938 + 0.0000000003 * T
    DD = np.array([1.0 - FD4, FD3 + FD5, FD3 + FD5 - FD4, FD3 + FD5 - FD1, -FD5, FD3 + FD5 - FD4 - FD2])
    freq = idd @ DD
    phase = np.mod(idd @ Dv, 360.0)
    return freq, phase


def _hd_spline(x, u):
    """HARDISP SPLINE: second derivatives of the interpolating cubic spline (special end conditions)."""
    n = len(x)
    if n <= 3:
        return np.zeros(n)
    q = lambda u1, x1, u2, x2: (u1 / x1 ** 2 - u2 / x2 ** 2) / (1.0 / x1 - 1.0 / x2)  # noqa: E731
    q1 = q(u[1] - u[0], x[1] - x[0], u[2] - u[0], x[2] - x[0])
    qn = q(u[n - 2] - u[n - 1], x[n - 2] - x[n - 1], u[n - 3] - u[n - 1], x[n - 3] - x[n - 1])
    s = np.zeros(n)
    a = np.zeros(n)
    s[0] = 6.0 * ((u[1] - u[0]) / (x[1] - x[0]) - q1)
    for i in range(1, n - 1):
        s[i] = (u[i - 1] / (x[i] - x[i - 1]) - u[i] * (1.0 / (x[i] - x[i - 1]) + 1.0 / (x[i + 1] - x[i]))
                + u[i + 1] / (x[i + 1] - x[i])) * 6.0
    s[n - 1] = 6.0 * (qn + (u[n - 2] - u[n - 1]) / (x[n - 1] - x[n - 2]))
    a[0] = 2.0 * (x[1] - x[0])
    a[1] = 1.5 * (x[1] - x[0]) + 2.0 * (x[2] - x[1])
    s[1] = s[1] - 0.5 * s[0]
    for i in range(2, n - 1):
        c = (x[i] - x[i - 1]) / a[i - 1]
        a[i] = 2.0 * (x[i + 1] - x[i - 1]) - c * (x[i] - x[i - 1])
        s[i] = s[i] - c * s[i - 1]
    c = (x[n - 1] - x[n - 2]) / a[n - 2]
    a[n - 1] = (2.0 - c) * (x[n - 1] - x[n - 2])
    s[n - 1] = s[n - 1] - c * s[n - 2]
    s[n - 1] = s[n - 1] / a[n - 1]
    for j in range(1, n):
        i = n - 1 - j
        s[i] = (s[i] - (x[i + 1] - x[i]) * s[i + 1]) / a[i]
    return s


def _hd_eval(y, x, u, s):
    n = len(x)
    if y <= x[0]:
        return u[0]
    if y >= x[n - 1]:
        return u[n - 1]
    k2 = int(np.nonzero((x[:-1] < y) & (x[1:] >= y))[0][-1]) + 1
    k1 = k2 - 1
    dy, dy1, dk = x[k2] - y, y - x[k1], x[k2] - x[k1]
    deli = 1.0 / (6.0 * dk)
    f1 = (s[k1] * dy ** 3 + s[k2] * dy1 ** 3) * deli
    f2 = dy1 * (u[k2] / dk - s[k2] * dk / 6.0)
    f3 = dy * (u[k1] / dk - s[k1] * dk / 6.0)
    return f1 + f2 + f3


def _admint(amp_in, ph_in, t_utc_ns):
    """HARDISP ADMINT: admittance interpolation from the 11 BLQ constituents to 342 (IERS 2010 § 7.1.2)."""
    rl, aim, rf = [], [], []
    for ll in range(11):
        kk = np.nonzero(np.all(_HD_IDD == _BLQ_IDT[ll], axis=1))[0][0]
        rl.append(amp_in[ll] * np.cos(D2R * ph_in[ll]) / abs(_HD_TAMP[kk]))
        aim.append(amp_in[ll] * np.sin(D2R * ph_in[ll]) / abs(_HD_TAMP[kk]))
        fr, _ = _tdfrph(_BLQ_IDT[ll:ll + 1], t_utc_ns)
        rf.append(fr[0])
    rf, rl, aim = np.array(rf), np.array(rl), np.array(aim)
    o = np.argsort(rf, kind="stable")
    rf, rl, aim = rf[o], rl[o], aim[o]
    lp = rf < 0.5
    di = (rf > 0.5) & (rf < 1.5)
    sd = (rf > 1.5) & (rf < 2.5)
    bands = {}
    for key, m in ((0, lp), (1, di), (2, sd)):
        x = rf[m]
        bands[key] = (x, rl[m], aim[m], _hd_spline(x, rl[m]), _hd_spline(x, aim[m]))
    freq, phase = _tdfrph(_HD_IDD, t_utc_ns)
    amp = np.zeros(len(_HD_TAMP))
    ph = phase.copy()
    for i in range(len(_HD_TAMP)):
        band = _HD_IDD[i, 0]
        if band == 0:
            ph[i] += 180.0
        elif band == 1:
            ph[i] += 90.0
        x, u_re, u_im, s_re, s_im = bands[band]
        re = _hd_eval(freq[i], x, u_re, s_re)
        am = _hd_eval(freq[i], x, u_im, s_im)
        amp[i] = _HD_TAMP[i] * np.hypot(re, am)
        ph[i] += np.degrees(np.arctan2(am, re))
    return amp, freq, ph


def ocean_loading(blq_amp, blq_phase, t_utc_ns):
    """Ocean tide loading displacement (dU, dS, dW) in metres at UTC instants (HARDISP algorithm).

    blq_amp/phase: (3,11) rows radial, west, south; phases are Greenwich lags (deg).
    """
    t = np.atleast_1d(np.asarray(t_utc_ns, dtype=np.int64))
    t0 = int(t[0])
    out = np.zeros((len(t), 3))
    dt_days = (t - t0) / ts.NS / 86400.0
    # HARDISP order: TAMP rows radial(1), west(2), south(3); output dU dS dW
    for comp, row in ((0, 0), (2, 1), (1, 2)):
        amp, freq, ph = _admint(blq_amp[row], -blq_phase[row], t0)
        arg = np.radians(ph)[None, :] + 2 * np.pi * freq[None, :] * dt_days[:, None]
        out[:, comp] = np.cos(arg) @ amp
    return out


def otl_enu_to_xyz(d_usw, lat, lon):
    """(dU, dS, dW) -> ECEF displacement."""
    from .coords import enu2xyz
    enu = np.stack([-d_usw[:, 2], -d_usw[:, 1], d_usw[:, 0]], axis=-1)
    return enu2xyz(enu, lat, lon)


# ============================================================================== pole tide
def secular_pole(decimal_year):
    """IERS conventional secular pole (mas), Conventions 2010 updated (2018) § 7.1.4 [verify_log V-013]."""
    dt = decimal_year - 2000.0
    return 55.0 + 1.677 * dt, 320.5 + 3.460 * dt


def pole_tide(lat, lon, decimal_year, xp_as, yp_as):
    """Solid Earth pole tide displacement (m, ECEF), IERS 2010 eq. 7.26 with the secular pole."""
    from .coords import enu2xyz
    xs, ys = secular_pole(decimal_year)
    m1 = xp_as - xs * 1e-3
    m2 = -(yp_as - ys * 1e-3)
    th = np.pi / 2 - lat
    s_r = -33.0 * np.sin(2 * th) * (m1 * np.cos(lon) + m2 * np.sin(lon)) * 1e-3
    s_th = -9.0 * np.cos(2 * th) * (m1 * np.cos(lon) + m2 * np.sin(lon)) * 1e-3
    s_la = 9.0 * np.cos(th) * (m1 * np.sin(lon) - m2 * np.cos(lon)) * 1e-3
    enu = np.stack([np.atleast_1d(s_la), -np.atleast_1d(s_th), np.atleast_1d(s_r)], axis=-1)
    return enu2xyz(enu, lat, lon)


# ============================================================================== propagation effects
def sagnac_rotate(r_sat, tau):
    """Rotate satellite ECEF position at transmit time into the ECEF frame at reception (TDS § 7.3)."""
    th = settings.OMEGA_E * np.asarray(tau)
    c, s = np.cos(th), np.sin(th)
    x, y, z = r_sat[..., 0], r_sat[..., 1], r_sat[..., 2]
    return np.stack([c * x + s * y, -s * x + c * y, z], axis=-1)


def relativistic_clock(r_sat, v_sat):
    """Periodic relativistic satellite clock correction dt_rel = -2 r.v / c^2 (s), TDS § 7.4."""
    return -2.0 * np.sum(r_sat * v_sat, axis=-1) / C ** 2


def shapiro(r_sat, r_rcv):
    """Gravitational (Shapiro) delay (m), IERS 2010 Ch. 11: (2GM/c^2) ln((rs+rr+rho)/(rs+rr-rho))."""
    rs = np.linalg.norm(r_sat, axis=-1)
    rr = np.linalg.norm(r_rcv, axis=-1)
    rho = np.linalg.norm(r_sat - r_rcv, axis=-1)
    return 2.0 * settings.GM_EARTH / C ** 2 * np.log((rs + rr + rho) / (rs + rr - rho))


# ============================================================================== attitude & wind-up
def nominal_attitude(r_sat, r_sun):
    """Nominal yaw-steering body axes (ex, ey, ez) in ECEF (IGS convention, TDS § 7.11).

    ez toward the Earth centre, ey = ez x e_sun (perpendicular to the Sun-satellite-Earth plane),
    ex = ey x ez (Sun-side). Arrays (n,3).
    """
    ez = -r_sat / np.linalg.norm(r_sat, axis=-1, keepdims=True)
    es = r_sun - r_sat
    es = es / np.linalg.norm(es, axis=-1, keepdims=True)
    ey = np.cross(ez, es)
    ey = ey / np.linalg.norm(ey, axis=-1, keepdims=True)
    ex = np.cross(ey, ez)
    return ex, ey, ez


def beta_and_mu(r_sat, v_sat, r_sun):
    """Sun elevation above the orbital plane beta (rad) and orbit angle from orbit noon mu (rad).

    v_sat must be inertial-like; the ECEF velocity plus w x r is used for the orbit normal.
    """
    w = np.array([0.0, 0.0, settings.OMEGA_E])
    v_in = v_sat + np.cross(w, r_sat)
    n = np.cross(r_sat, v_in)
    n = n / np.linalg.norm(n, axis=-1, keepdims=True)
    es = r_sun / np.linalg.norm(r_sun, axis=-1, keepdims=True)
    beta = np.pi / 2 - np.arccos(np.clip(np.sum(n * es, axis=-1), -1, 1))
    # projection of Sun direction into orbit plane
    sp = es - np.sum(es * n, axis=-1, keepdims=True) * n
    sp = sp / np.linalg.norm(sp, axis=-1, keepdims=True)
    er = r_sat / np.linalg.norm(r_sat, axis=-1, keepdims=True)
    cosmu = np.clip(np.sum(er * sp, axis=-1), -1, 1)
    mu = np.arccos(cosmu)                                   # 0 at orbit noon, pi at orbit midnight
    return beta, mu


def in_shadow(r_sat, r_sun):
    """Cylindrical Earth-shadow test (bool)."""
    es = r_sun / np.linalg.norm(r_sun, axis=-1, keepdims=True)
    proj = np.sum(r_sat * es, axis=-1)
    perp = np.linalg.norm(r_sat - proj[..., None] * es, axis=-1)
    return (proj < 0) & (perp < settings.R_EARTH_SHADOW)


def windup(ex_s, ey_s, r_sat, r_rcv, lat, lon, prev=None):
    """Carrier phase wind-up (cycles), Wu et al. (1993); formulation as RTKLIB windupcorr.

    Arrays over epochs of one arc; unwrapped continuously (initial value arbitrary, absorbed by the
    ambiguity). Returns array of cycles.
    """
    from .coords import rot_xyz2enu
    E = rot_xyz2enu(lat, lon)
    exr, eyr = E[1], -E[0]                                  # receiver x = north, y = west
    ek = r_rcv - r_sat
    ek = ek / np.linalg.norm(ek, axis=-1, keepdims=True)
    eks = np.cross(ek, ey_s)
    ekr = np.cross(ek, eyr[None, :])
    ds = ex_s - ek * np.sum(ek * ex_s, axis=-1, keepdims=True) - eks
    dr = exr[None, :] - ek * (ek @ exr)[:, None] + ekr
    cosp = np.sum(ds * dr, axis=-1) / np.linalg.norm(ds, axis=-1) / np.linalg.norm(dr, axis=-1)
    ph = np.arccos(np.clip(cosp, -1, 1)) / (2 * np.pi)
    drs = np.cross(ds, dr)
    ph = np.where(np.sum(ek * drs, axis=-1) < 0, -ph, ph)
    out = np.empty_like(ph)
    last = prev if prev is not None else ph[0]
    for i, p in enumerate(ph):
        last = p + np.floor(last - p + 0.5)
        out[i] = last
    return out


def eclipse_exclusion(block, beta, mu, shadow, t_ns, post_shadow_s=None):
    """Boolean mask of epochs to exclude when no ORBEX attitude is available (TDS § 7.11).

    Shadow crossing + post-shadow recovery, and noon turns with |beta| < beta0(block) inside the
    configured orbit-angle window around orbit noon.
    """
    post = post_shadow_s if post_shadow_s is not None else settings.ECLIPSE_POST_SHADOW_S
    b0 = settings.ECLIPSE_BETA_NOON_DEG.get(block, settings.ECLIPSE_BETA_NOON_DEG["DEFAULT"]) * D2R
    excl = shadow.copy()
    idx = np.nonzero(shadow)[0]
    if len(idx):
        last = t_ns[idx]
        for k in range(len(t_ns)):
            prev = last[last <= t_ns[k]]
            if len(prev) and (t_ns[k] - prev[-1]) <= post * ts.NS:
                excl[k] = True
    noon = (np.abs(beta) < b0) & (mu < settings.ECLIPSE_NOON_WINDOW_DEG * D2R)
    return excl | noon
