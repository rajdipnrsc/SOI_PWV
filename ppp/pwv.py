"""ZTD -> PWV conversion and uncertainty propagation (TDS § 18). Separate from the PPP estimator.

PWV = Pi * (ZTD - ZHD_met),   Pi = 1e6 / (rho_w R_v (k3/Tm + k2'))          (Bevis et al. 1994)
ZHD_met = 0.0022768 P_ant / (1 - 0.00266 cos 2phi - 0.28e-6 H_ant)          (Saastamoinen / Davis 1985)
"""
import os
from dataclasses import dataclass

import numpy as np

from . import settings
from .troposphere import saastamoinen_zhd


def pi_factor(tm_k, constants=None):
    c = settings.PWV_CONSTANTS[constants or settings.PWV_CONSTANTS_SET]
    return 1e6 / (c["rho_w"] * c["R_v"] * (c["k3"] / np.asarray(tm_k) + c["k2p"]))


def sigma_pi(tm_k, sig_tm, constants=None):
    """sigma_Pi = Pi * (k3/Tm^2)/(k3/Tm + k2') * sigma_Tm (TDS § 18.4)."""
    c = settings.PWV_CONSTANTS[constants or settings.PWV_CONSTANTS_SET]
    tm = np.asarray(tm_k)
    return pi_factor(tm, constants) * (c["k3"] / tm ** 2) / (c["k3"] / tm + c["k2p"]) * sig_tm


def reduce_pressure(p1_hpa, h1_m, h2_m, tv_k):
    """Hypsometric reduction P2 = P1 exp(-g (H2 - H1) / (R_d Tv)) (TDS § 18.3)."""
    return p1_hpa * np.exp(-settings.G0 * (h2_m - h1_m) / (settings.RD_DRY * tv_k))


@dataclass
class PWVResult:
    pwv: np.ndarray              # mm
    sig_ppp: np.ndarray
    sig_conv: np.ndarray
    sig_total: np.ndarray
    p_ant: np.ndarray            # hPa
    p_source: np.ndarray         # label per epoch
    tm: np.ndarray               # K
    tm_source: np.ndarray
    pi: np.ndarray
    zhd_met: np.ndarray          # mm
    zwd_met: np.ndarray          # mm
    constants_set: str
    flags: set


def convert(ztd_m, sig_ztd_m, lat_rad, h_orth_m, p_ant_hpa, p_source, tm_k, tm_source,
            sig_ztd_coord_m=0.0, constants=None):
    """Convert ZTD series (m) to PWV (mm) with the PPP and conversion uncertainties kept separate."""
    flags = set()
    cset = constants or settings.PWV_CONSTANTS_SET
    ztd = np.asarray(ztd_m, dtype=float)
    p = np.broadcast_to(np.asarray(p_ant_hpa, dtype=float), ztd.shape)
    tm = np.broadcast_to(np.asarray(tm_k, dtype=float), ztd.shape)
    zhd = saastamoinen_zhd(p, lat_rad, h_orth_m)
    zwd = ztd - zhd
    Pi = pi_factor(tm, cset)
    pwv = Pi * zwd
    f = 1 - 0.00266 * np.cos(2 * lat_rad) - 0.28e-6 * h_orth_m
    # sources may be one label or one label per epoch (barometer gaps filled from the next source, TDS § 18.2)
    p_src = np.broadcast_to(np.asarray(p_source, dtype=object), ztd.shape)
    tm_src = np.broadcast_to(np.asarray(tm_source, dtype=object), ztd.shape)
    sp = np.array([settings.SIGMA_P_HPA.get(x, 5.0) for x in p_src.ravel()]).reshape(ztd.shape)
    stm = np.array([settings.SIGMA_TM_K.get(x, 5.0) for x in tm_src.ravel()]).reshape(ztd.shape)
    s_zhd = 0.0022768 / f * sp
    s_pi = sigma_pi(tm, stm, cset)
    s_ppp = Pi * np.asarray(sig_ztd_m) * settings.S_ZTD_SCALE
    s_conv = np.sqrt((Pi * s_zhd) ** 2 + (zwd * s_pi) ** 2)
    s_tot = np.sqrt(s_ppp ** 2 + s_conv ** 2 + (Pi * sig_ztd_coord_m) ** 2)
    if np.any(p_src == "GPT3"):
        flags.add("P_CLIMATOLOGY")
    return PWVResult(pwv * 1e3, s_ppp * 1e3, s_conv * 1e3, s_tot * 1e3, np.asarray(p), p_src, np.asarray(tm),
                     tm_src, Pi, zhd * 1e3, zwd * 1e3, cset, flags)


# ================================================================== station meteorology (TDS § 18.2)
def sat_vapour_pressure(t_c):
    """Saturation vapour pressure over water (hPa), Magnus form (WMO 2008, CIMO guide annex 4.B)."""
    t_c = np.asarray(t_c, dtype=float)
    return 6.112 * np.exp(17.62 * t_c / (243.12 + t_c))


def virtual_temperature(t_k, p_hpa, rh_pct=None):
    """T_v = T (1 + 0.608 q) with q from relative humidity (dry air if RH unknown)."""
    t_k = np.asarray(t_k, dtype=float)
    if rh_pct is None:
        return t_k
    e = np.clip(np.asarray(rh_pct, dtype=float), 0, 100) / 100.0 * sat_vapour_pressure(t_k - 273.15)
    q = 0.622 * e / (np.asarray(p_hpa) - 0.378 * e)
    return t_k * (1 + 0.608 * q)


@dataclass
class Barometer:
    """Station barometer pressure reduced to the antenna height (hPa), epochs GPST ns."""
    t: np.ndarray
    p_ant: np.ndarray
    h_sensor: float               # orthometric (m)
    h_ant: float
    height_source: str            # SENSOR_POS | ASSUMED_ARP
    n_rejected: int
    info: str
    t_surface_k: np.ndarray = None


def barometer_from_met(met, h_ant_orth, undulation, t_fallback_k=None):
    """QC + hypsometric reduction of RINEX met pressure to the ARP (TDS § 18.3).

    Returns (Barometer or None, list of warning codes).  Values outside settings.MET_P_RANGE_HPA and spikes larger
    than settings.MET_SPIKE_HPA against a 5-sample running median are rejected."""
    warns = []
    if "PR" not in met.values:
        return None, ["MET_NO_PRESSURE"]
    p = met.values["PR"].astype(float).copy()
    lo, hi = settings.MET_P_RANGE_HPA
    bad = ~np.isfinite(p) | (p < lo) | (p > hi)
    p[bad] = np.nan
    if np.sum(np.isfinite(p)) >= 5:
        pad = np.pad(p, 2, mode="edge")
        med = np.nanmedian(np.lib.stride_tricks.sliding_window_view(pad, 5), axis=1)
        spike = np.isfinite(p) & (np.abs(p - med) > settings.MET_SPIKE_HPA)
        p[spike] = np.nan
        bad |= spike
    ok = np.isfinite(p)
    if ok.sum() < 2:
        return None, ["MET_NO_VALID_PRESSURE"]
    pos = met.sensor_pos.get("PR")
    if pos is not None and (pos[3] != 0.0 or any(pos[:3])):
        h_ell = pos[3]
        if h_ell == 0.0:
            from .coords import ecef_to_geodetic
            h_ell = ecef_to_geodetic(np.array(pos[:3]))[2]
        h_sens = h_ell - (undulation or 0.0)
        hsrc = "SENSOR_POS"
        if undulation is None:
            warns.append("MET_HEIGHT_NO_GEOID")
    else:
        h_sens = h_ant_orth
        hsrc = "ASSUMED_ARP"
        warns.append("MET_HEIGHT_ASSUMED")
    if abs(h_sens - h_ant_orth) > settings.MET_MAX_HEIGHT_DIFF_M:
        return None, warns + ["MET_SENSOR_TOO_FAR"]
    td = met.values.get("TD")
    if td is not None and np.sum(np.isfinite(td)) >= 2:
        tk = np.interp(met.t.astype(float), met.t[np.isfinite(td)].astype(float), td[np.isfinite(td)]) + 273.15
    else:
        tk = np.full(len(p), t_fallback_k if t_fallback_k is not None else 288.15)
        warns.append("MET_NO_TEMPERATURE")
    hr = met.values.get("HR")
    if hr is not None and np.sum(np.isfinite(hr)) >= 2:
        hr = np.interp(met.t.astype(float), met.t[np.isfinite(hr)].astype(float), hr[np.isfinite(hr)])
    else:
        hr = None
    # mean virtual temperature of the layer between sensor and antenna (standard lapse rate 6.5 K/km)
    tv = virtual_temperature(tk - 0.0065 * 0.5 * (h_ant_orth - h_sens), p, hr)
    p_ant = reduce_pressure(p[ok], h_sens, h_ant_orth, tv[ok])
    info = met.sensor_info.get("PR", "")
    return Barometer(met.t[ok], p_ant, float(h_sens), float(h_ant_orth), hsrc, int(bad.sum()),
                     f"{os.path.basename(met.path)} ({info})" if info else os.path.basename(met.path), tk[ok]), warns


def pressure_at(baro, t_query, max_gap_s=None):
    """Barometer pressure at query epochs (GPST ns); NaN inside gaps longer than max_gap_s or outside the record."""
    max_gap = (settings.MET_MAX_GAP_S if max_gap_s is None else max_gap_s) * 1e9
    tq = np.asarray(t_query, dtype=float)
    tb = baro.t.astype(float)
    p = np.interp(tq, tb, baro.p_ant)
    k = np.clip(np.searchsorted(tb, tq), 1, len(tb) - 1)
    gap = tb[k] - tb[k - 1]
    near = np.minimum(np.abs(tq - tb[k - 1]), np.abs(tb[k] - tq))
    bad = (tq < tb[0] - max_gap / 2) | (tq > tb[-1] + max_gap / 2) | ((gap > max_gap) & (near > max_gap / 2))
    p[bad] = np.nan
    return p


def era5_station_met(bg, lat_deg, lon_deg, h_orth, t_utc):
    """ERA5 pressure (hPa) and Tm (K) at the antenna by log-pressure interpolation of the geopotential profile and
    profile integration (TDS § 18.2 priority 2 for P, priority 1 for Tm), linear in time between ERA5 hours."""
    from .mapping import background_at
    tb = bg.t_utc.astype(float)
    ps, tm = [], []
    for t in bg.t_utc:
        _pwv, p_, tm_ = background_at(bg, lat_deg, lon_deg, h_orth, int(t))
        ps.append(float(p_[0]))
        tm.append(float(tm_[0]))
    tq = np.asarray(t_utc, dtype=float)
    inside = (tq >= tb[0]) & (tq <= tb[-1])
    p = np.where(inside, np.interp(tq, tb, ps), np.nan)
    t_m = np.where(inside, np.interp(tq, tb, tm), np.nan)
    return p, t_m
