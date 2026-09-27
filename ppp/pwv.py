"""ZTD -> PWV conversion and uncertainty propagation (TDS § 18). Separate from the PPP estimator.

PWV = Pi * (ZTD - ZHD_met),   Pi = 1e6 / (rho_w R_v (k3/Tm + k2'))          (Bevis et al. 1994)
ZHD_met = 0.0022768 P_ant / (1 - 0.00266 cos 2phi - 0.28e-6 H_ant)          (Saastamoinen / Davis 1985)
"""
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
    p_source: str
    tm: np.ndarray               # K
    tm_source: str
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
    sp = settings.SIGMA_P_HPA.get(p_source, 5.0)
    stm = settings.SIGMA_TM_K.get(tm_source, 5.0)
    s_zhd = 0.0022768 / f * sp
    s_pi = sigma_pi(tm, stm, cset)
    s_ppp = Pi * np.asarray(sig_ztd_m) * settings.S_ZTD_SCALE
    s_conv = np.sqrt((Pi * s_zhd) ** 2 + (zwd * s_pi) ** 2)
    s_tot = np.sqrt(s_ppp ** 2 + s_conv ** 2 + (Pi * sig_ztd_coord_m) ** 2)
    if p_source == "GPT3":
        flags.add("P_CLIMATOLOGY")
    return PWVResult(pwv * 1e3, s_ppp * 1e3, s_conv * 1e3, s_tot * 1e3, np.asarray(p), p_source, np.asarray(tm),
                     tm_source, Pi, zhd * 1e3, zwd * 1e3, cset, flags)
