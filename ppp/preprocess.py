"""Preprocessing and quality control (TDS § 13): observable selection, clock jumps, LLI / Melbourne-Wuebbena
/ geometry-free slip detection, ROTI, gaps, arc construction and short-arc removal, IF formation.

Every detected slip starts a new ambiguity arc (no repair in v1, TDS § 13.11, F11).
"""
import warnings
from dataclasses import dataclass, field

import numpy as np

from . import log as plog
from . import settings

LOG = plog.get()
C = settings.C_LIGHT


@dataclass
class SelectedObs:
    """Per-epoch, per-satellite selected dual-frequency observables (metres), GPS L1/L2."""
    epochs: np.ndarray
    sats: list
    P1: np.ndarray
    P2: np.ndarray
    L1: np.ndarray               # metres
    L2: np.ndarray
    S1: np.ndarray
    S2: np.ndarray
    lli: np.ndarray              # bool: LLI bit 0 on L1 or L2
    pair: np.ndarray             # int code of the (C1, L1, C2, L2) combination (arc break on change)
    pair_names: list             # index -> (C1, L1, C2, L2)
    f1: float
    f2: float

    @property
    def lam1(self):
        return C / self.f1

    @property
    def lam2(self):
        return C / self.f2


def select_observables(obs, system="G"):
    """Per satellite and epoch, choose the highest-priority available code and phase on each band
    (TDS § 13.3). A change of the selected pair within an arc starts a new arc (TDS § 5.7)."""
    b1, b2 = settings.BANDS[system]
    f1, f2 = settings.FREQ[(system, b1)], settings.FREQ[(system, b2)]
    lam = {b1: C / f1, b2: C / f2}
    ne, ns = len(obs.epochs), len(obs.sats)
    sys_mask = np.array([s[0] == system for s in obs.sats])

    def pick(prio):
        """Per satellite: the highest-priority observable whose availability is >= 90 % of the best one
        (stable choice over the window; missing epochs become gaps, not observable switches)."""
        val = np.full((ne, ns), np.nan)
        which = np.full((ne, ns), -1, dtype=int)
        avail = np.zeros((len(prio), ns))
        for k, code in enumerate(prio):
            if code in obs.obs:
                v = obs.obs[code]
                avail[k] = np.sum(np.isfinite(v) & (v != 0.0), axis=0)
        best = avail.max(axis=0)
        for j in range(ns):
            if best[j] == 0:
                continue
            k = next(k for k in range(len(prio)) if avail[k, j] >= 0.9 * best[j])
            v = obs.obs[prio[k]][:, j]
            m = np.isfinite(v) & (v != 0.0)
            val[m, j] = v[m]
            which[m, j] = k
        return val, which

    P1, cp1 = pick(settings.CODE_PRIORITY[system][b1])
    P2, cp2 = pick(settings.CODE_PRIORITY[system][b2])
    L1c, lp1 = pick(settings.PHASE_PRIORITY[system][b1])
    L2c, lp2 = pick(settings.PHASE_PRIORITY[system][b2])
    lli = np.zeros((ne, ns), dtype=bool)
    for band, which, prio in ((b1, lp1, settings.PHASE_PRIORITY[system][b1]),
                              (b2, lp2, settings.PHASE_PRIORITY[system][b2])):
        for k, code in enumerate(prio):
            if code in obs.lli:
                m = which == k
                lli |= m & ((obs.lli[code] & 1) == 1)
    snr = {}
    for band, which in ((b1, lp1), (b2, lp2)):
        prio = settings.PHASE_PRIORITY[system][band]
        s = np.full((ne, ns), np.nan)
        for k, code in enumerate(prio):
            sc = "S" + code[1:]
            if sc in obs.obs:
                m = which == k
                s[m] = obs.obs[sc][m]
        snr[band] = s
    key = cp1 * 1000 + lp1 * 100 + cp2 * 10 + lp2
    valid = np.isfinite(P1) & np.isfinite(P2) & np.isfinite(L1c) & np.isfinite(L2c) & sys_mask[None, :]
    names = {}
    for kk in np.unique(key[valid]):
        a, r = divmod(int(kk), 1000)
        b, r = divmod(r, 100)
        c, d = divmod(r, 10)
        names[int(kk)] = (settings.CODE_PRIORITY[system][b1][a], settings.PHASE_PRIORITY[system][b1][b],
                          settings.CODE_PRIORITY[system][b2][c], settings.PHASE_PRIORITY[system][b2][d])
    for arr in (P1, P2, L1c, L2c):
        arr[~valid] = np.nan
    key = np.where(valid, key, -1)
    return SelectedObs(obs.epochs, obs.sats, P1, P2, L1c * lam[b1], L2c * lam[b2], snr[b1], snr[b2], lli & valid,
                       key, names, f1, f2)


def if_combination(sel):
    a1 = sel.f1 ** 2 / (sel.f1 ** 2 - sel.f2 ** 2)
    a2 = sel.f2 ** 2 / (sel.f1 ** 2 - sel.f2 ** 2)
    return a1 * sel.P1 - a2 * sel.P2, a1 * sel.L1 - a2 * sel.L2, a1, a2


def apply_osb(sel, osb, t_mid, phase=False):
    """Subtract satellite OSBs per exact RINEX code (TDS § 11.5): code always, phase too when phase=True
    (PPP-AR: "apply satellite OSBs to all observations"). Missing phase OSB -> phase set to NaN (no AR).

    Returns (P1c, P2c, L1c, L2c, info) with info listing missing biases.
    """
    P1, P2, L1, L2 = sel.P1.copy(), sel.P2.copy(), sel.L1.copy(), sel.L2.copy()
    missing = set()
    applied = 0
    if osb is None:
        return P1, P2, L1, L2, {"applied": 0, "missing": ["no OSB product"]}
    arrays = [(P1, 0), (P2, 2)] + ([(L1, 1), (L2, 3)] if phase else [])
    for j, prn in enumerate(sel.sats):
        for kk, codes in sel.pair_names.items():
            m = sel.pair[:, j] == kk
            if not m.any():
                continue
            for arr, ci in arrays:
                code = codes[ci]
                v = osb.get(prn, code, t_mid)
                if v is None:
                    missing.add(f"{prn}:{code}")
                    if ci in (1, 3):
                        arr[m, j] = np.nan          # phase without OSB: not usable for AR
                    # code without OSB: kept uncorrected, flagged CODE_BIAS_MISSING by the caller
                else:
                    arr[m, j] -= v
                    applied += 1
    return P1, P2, L1, L2, {"applied": applied, "missing": sorted(missing)}


def verify_osb_sign(obs, osb, t_mid):
    """Data check of the OSB sign convention (TDS § 11.5 VERIFY): where both C1C and C1W are observed,
    (C1C - OSB_C1C) - (C1W - OSB_C1W) must be common to all satellites (receiver bias only)."""
    if osb is None or "C1C" not in obs.obs or "C1W" not in obs.obs:
        return None
    d_minus, d_plus = [], []
    for j, prn in enumerate(obs.sats):
        a, b = osb.get(prn, "C1C", t_mid), osb.get(prn, "C1W", t_mid)
        if a is None or b is None:
            continue
        diff = obs.obs["C1C"][:, j] - obs.obs["C1W"][:, j]
        diff = diff[np.isfinite(diff)]
        if len(diff) < 20:
            continue
        med = np.median(diff)
        d_minus.append(med - (a - b))
        d_plus.append(med + (a - b))
    if len(d_minus) < 5:
        return None
    s_minus, s_plus = np.std(d_minus), np.std(d_plus)
    return {"scatter_subtract_m": float(s_minus), "scatter_add_m": float(s_plus),
            "consistent_with": "subtract" if s_minus <= s_plus else "add", "n_sat": len(d_minus)}


# ================================================================================ clock jumps
def detect_clock_jumps(sel, valid):
    """Detect receiver millisecond clock jumps (TDS § 13.8): common jump in (P_IF - L_IF) of k*c*1ms across
    >= 80 % of satellites. Integer-ms jumps are repaired by adding k*c*1e-3 to the phases after the jump.
    Returns list of (epoch_index, k_ms, repaired)."""
    Pif, Lif, _, _ = if_combination(sel)
    d = Pif - Lif
    out = []
    ne = len(sel.epochs)
    for k in range(1, ne):
        m = valid[k] & valid[k - 1] & np.isfinite(d[k]) & np.isfinite(d[k - 1])
        if m.sum() < 4:
            continue
        jump = d[k, m] - d[k - 1, m]
        ms = jump / (C * 1e-3)
        med = np.median(ms)
        if abs(med) < 0.5:
            continue
        frac = np.mean(np.abs(ms - med) < 0.01)
        if frac < settings.CLOCK_JUMP_FRACTION:
            continue
        kint = int(np.round(med))
        repaired = abs(med - kint) < settings.CLOCK_JUMP_TOL_MS * 100
        if repaired:
            sel.L1[k:] += kint * C * 1e-3
            sel.L2[k:] += kint * C * 1e-3
        out.append((k, float(med), bool(repaired)))
    return out


# ================================================================================ slips and arcs
@dataclass
class ArcInfo:
    arc: np.ndarray                    # (ne, ns) arc id, -1 = no usable data
    reasons: dict = field(default_factory=dict)      # counts of reset causes
    n_arcs: int = 0
    iono_active: np.ndarray = None     # (ne,) bool
    roti: np.ndarray = None            # (ne,) station ROTI (TECU/min)
    arc_meta: dict = field(default_factory=dict)     # arc id -> dict(sat index, first, last, pair)


def _roti(gf, t_s, valid):
    """Rate-of-TEC index per epoch (median over satellites of the 5-min std of ROT), TECU/min."""
    tecu_m = 0.105                     # 1 TECU ~ 0.105 m in L1-L2 [TDS § 13.7, derived]
    ne, ns = gf.shape
    rot = np.full((ne, ns), np.nan)
    dt = np.diff(t_s)
    for j in range(ns):
        g = gf[:, j]
        r = np.diff(g) / tecu_m / (dt / 60.0)
        ok = valid[1:, j] & valid[:-1, j] & (dt <= 60.0)
        r[~ok] = np.nan
        rot[1:, j] = r
    win = max(int(round(settings.ROTI_WINDOW_S / np.median(dt))), 2) if len(dt) else 10
    roti_sat = np.full((ne, ns), np.nan)
    warnings.filterwarnings("ignore", message=".*(Degrees of freedom|All-NaN slice|Mean of empty).*")
    for k in range(ne):
        a = max(0, k - win + 1)
        blk = rot[a:k + 1]
        n = np.sum(np.isfinite(blk), axis=0)
        with np.errstate(invalid="ignore"):
            s = np.nanstd(blk, axis=0)
        s[n < 3] = np.nan
        roti_sat[k] = s
    with np.errstate(all="ignore"):
        station = np.nanmedian(roti_sat, axis=1)
    return station


def detect_slips_and_arcs(sel, t_s, el, valid, cutoff_rad):
    """Build ambiguity arcs (TDS § 5.7, § 13.6-13.7). Returns ArcInfo.

    t_s: epoch times (s, float relative), el: (ne, ns) elevation (rad), valid: (ne, ns) usable data mask
    (products + model available). Elevation cutoff applied here; re-rising after cutoff = new arc.
    """
    ne, ns = sel.P1.shape
    lw = C / (sel.f1 - sel.f2)
    mw = (sel.f1 * sel.L1 - sel.f2 * sel.L2) / (sel.f1 - sel.f2) - (sel.f1 * sel.P1 + sel.f2 * sel.P2) / (sel.f1 + sel.f2)
    mwc = mw / lw
    gf = sel.L1 - sel.L2
    usable = valid & np.isfinite(mwc) & np.isfinite(gf) & (el >= cutoff_rad)
    if settings.SNR_MASK:
        usable &= ~(sel.S1 < settings.SNR_MIN["1"]) & ~(sel.S2 < settings.SNR_MIN["2"])
    roti = _roti(gf, t_s, usable)
    iono = np.nan_to_num(roti, nan=0.0) > settings.ROTI_HIGH_TECU_MIN
    arc = np.full((ne, ns), -1, dtype=int)
    reasons = {"LLI": 0, "MW": 0, "GF": 0, "GAP": 0, "OBS_CHANGE": 0, "RISE": 0}
    nid = 0
    meta = {}
    dt_nom = float(np.median(np.diff(t_s))) if ne > 1 else 30.0
    for j in plog.progress(range(ns), desc="Cycle-slip detection", unit="sat"):
        ks_use = np.nonzero(usable[:, j])[0]
        cur = -1
        last_k = None
        mu = s2 = 0.0
        n = 0
        gf_hist = []                 # (t, gf) of the current arc
        for pos, k in enumerate(ks_use):
            why = None
            new = cur < 0
            if not new:
                gap = t_s[k] - t_s[last_k]
                if gap > settings.GAP_RESET_S:
                    new, why = True, ("RISE" if (el[last_k, j] < cutoff_rad + 0.05 and gap > 600) else "GAP")
                elif sel.pair[k, j] != sel.pair[last_k, j]:
                    new, why = True, "OBS_CHANGE"
                elif sel.lli[k, j]:
                    new, why = True, "LLI"
            if not new and n >= 2:
                # Melbourne-Wuebbena (Blewitt 1990), confirmed by the next epoch (slip vs outlier)
                # running arc std, floored by an elevation-dependent code-noise model [A, D-008]
                sd = max(np.sqrt(s2 / (n - 1)), settings.MW_SIGMA_FLOOR_CYC / max(np.sin(el[k, j]), 0.1))
                thr = max(settings.MW_K * sd, settings.MW_MIN_CYC)
                if abs(mwc[k, j] - mu) > thr:
                    nxt = ks_use[pos + 1] if pos + 1 < len(ks_use) else None
                    confirmed = nxt is not None and (t_s[nxt] - t_s[k]) <= settings.GAP_RESET_S and \
                        abs(mwc[nxt, j] - mwc[k, j]) < thr and abs(mwc[nxt, j] - mu) > thr
                    if confirmed:
                        new, why = True, "MW"
                    else:
                        reasons["MW_OUTLIER"] = reasons.get("MW_OUTLIER", 0) + 1
                        continue                           # isolated outlier: epoch dropped
            if not new and len(gf_hist) >= 3:
                th = np.array([h[0] for h in gf_hist[-settings.GF_POLY_EPOCHS:]])
                gh = np.array([h[1] for h in gf_hist[-settings.GF_POLY_EPOCHS:]])
                deg = min(settings.GF_POLY_DEGREE, len(th) - 2)
                cf = np.polyfit(th - th[-1], gh, deg)
                pred = np.polyval(cf, t_s[k] - th[-1])
                sfit = np.std(gh - np.polyval(cf, th - th[-1])) if len(th) > deg + 1 else 0.0
                # g(dt, e): linear in the sampling gap, sqrt(1/sin e) capped at 2 for elevation [A, D-008]
                scale = max(1.0, (t_s[k] - th[-1]) / 30.0) * min(np.sqrt(1.0 / max(np.sin(el[k, j]), 0.05)), 2.0)
                if iono[k]:
                    scale *= 3.0                           # IONO_ACTIVE: GF thresholds scaled up
                thr_gf = max(settings.GF_K * sfit, settings.GF_MIN_M * scale)
                if abs(gf[k, j] - pred) > thr_gf:
                    # confirm with the next epoch (level shift persists) to separate slips from outliers
                    nxt = ks_use[pos + 1] if pos + 1 < len(ks_use) else None
                    jump = gf[k, j] - pred
                    jump_n = gf[nxt, j] - np.polyval(cf, t_s[nxt] - th[-1]) if nxt is not None else 0.0
                    confirmed = nxt is not None and (t_s[nxt] - t_s[k]) <= settings.GAP_RESET_S and \
                        abs(jump_n) > thr_gf and abs(jump_n - jump) < thr_gf
                    if confirmed or nxt is None:
                        new, why = True, "GF"
                    else:
                        reasons["GF_OUTLIER"] = reasons.get("GF_OUTLIER", 0) + 1
                        continue
            if new:
                nid += 1
                cur = nid
                meta[cur] = {"sat": j, "first": int(k), "pair": int(sel.pair[k, j])}
                mu, s2, n = 0.0, 0.0, 0
                gf_hist = []
                if why:
                    reasons[why] += 1
            arc[k, j] = cur
            n += 1
            d = mwc[k, j] - mu
            mu += d / n
            s2 += d * (mwc[k, j] - mu)
            gf_hist.append((t_s[k], gf[k, j]))
            last_k = k
    # short-arc removal
    n_short = 0
    for aid in list(meta):
        m = arc == aid
        ks = np.nonzero(m.any(axis=1))[0]
        if len(ks) == 0:
            del meta[aid]
            continue
        span = t_s[ks[-1]] - t_s[ks[0]] + dt_nom
        if span < settings.MIN_ARC_S:
            arc[m] = -1
            del meta[aid]
            n_short += 1
        else:
            meta[aid]["first"] = int(ks[0])
            meta[aid]["last"] = int(ks[-1])
    reasons["SHORT_ARCS_REMOVED"] = n_short
    return ArcInfo(arc, reasons, len(meta), iono, roti, meta)


def split_arc(arcinfo, j, k_from, cause):
    """Start a new ambiguity arc for satellite j from epoch index k_from (outlier-induced split)."""
    arc = arcinfo.arc
    aid = arc[k_from, j]
    if aid < 0:
        return None
    m = (arc[:, j] == aid)
    m[:k_from] = False
    new = max(arcinfo.arc_meta) + 1 if arcinfo.arc_meta else 1
    arc[m, j] = new
    ks = np.nonzero(m)[0]
    arcinfo.arc_meta[new] = {"sat": j, "first": int(ks[0]), "last": int(ks[-1]),
                             "pair": arcinfo.arc_meta[aid]["pair"]}
    ks_old = np.nonzero(arc[:, j] == aid)[0]
    if len(ks_old):
        arcinfo.arc_meta[aid]["last"] = int(ks_old[-1])
    else:
        del arcinfo.arc_meta[aid]
    arcinfo.reasons[cause] = arcinfo.reasons.get(cause, 0) + 1
    return new


def split_at_epochs(arcinfo, k_list, cause):
    """Start new arcs for every satellite at each epoch index in k_list (e.g. product-day boundaries,
    TDS § 11.2: ambiguities are reset at product-day boundaries unless continuity is verified)."""
    n = 0
    for kb in k_list:
        for j in range(arcinfo.arc.shape[1]):
            a = int(arcinfo.arc[kb, j]) if kb < arcinfo.arc.shape[0] else -1
            if a < 0:
                # first epoch at/after kb in the arc that spans the boundary
                col = arcinfo.arc[kb:, j]
                nxt = np.nonzero(col >= 0)[0]
                if not len(nxt):
                    continue
                a = int(col[nxt[0]])
                kk = kb + int(nxt[0])
            else:
                kk = kb
            meta = arcinfo.arc_meta.get(a)
            if meta is not None and meta["first"] < kb <= meta["last"]:
                split_arc(arcinfo, j, kk, cause)
                n += 1
    return n
