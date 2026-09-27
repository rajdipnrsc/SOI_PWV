"""Float PPP estimator: extended Kalman filter (Joseph form) + RTS smoother, IGG-III innovation screening,
post-fit residual editing, and 5-minute UTC extraction with QC flags (TDS § 5, § 6, § 13.10, § 16).

State registry (TDS § 5.1): an ordered list of keys ('ZWD', 'GN', 'GE', 'CLK[G]', 'DX', 'DY', 'DZ',
'AMB[<arc>]'); the filter works on keys, so states can appear/disappear between epochs.
Coordinates are estimated as corrections (DX, DY, DZ) to the linearisation point X0.
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from . import coords as co
from . import log as plog
from . import preprocess as pp
from . import settings
from . import timesys as ts

LOG = plog.get()
C = settings.C_LIGHT


@dataclass
class ObsSet:
    """Everything the estimator needs, prepared per epoch/satellite (all metres unless stated)."""
    epochs: np.ndarray           # int64 ns GPST
    t_s: np.ndarray              # float seconds since first epoch
    sats: list
    Pif: np.ndarray              # IF code (OSB corrected)
    Lif: np.ndarray              # IF phase (m)
    rho: np.ndarray              # modelled range at X0 (all deterministic terms except troposphere, clocks)
    dts: np.ndarray              # c * satellite clock (incl. relativistic)
    wind: np.ndarray             # lambda_NL * wind-up (m)
    el: np.ndarray
    az: np.ndarray
    los: np.ndarray              # (ne, ns, 3)
    mh: np.ndarray
    mw: np.ndarray
    mg: np.ndarray
    zhd0: np.ndarray             # (ne,) a priori ZHD
    zwd0: np.ndarray             # (ne,) a priori ZWD (initialisation only)
    usable: np.ndarray           # (ne, ns) bool
    f_if: float = 2.98
    X0: np.ndarray = None


@dataclass
class EstConfig:
    mode: str = "fixed"                      # fixed | static | constrained
    sigma_neu: tuple = settings.CONSTRAINED_SIGMA_NEU_M
    sigma_zwd: float = settings.SIGMA_ZWD_M_SQRT_H
    sigma_grad: float = settings.SIGMA_GRAD_M_SQRT_H
    gradients: bool = settings.ESTIMATE_GRADIENTS
    keep_ambiguities: bool = False           # AR pass: retain ambiguities to the window end (TDS § 5.8)
    sd_constraints: list = field(default_factory=list)    # [(arc_a, arc_b, value_m, sigma_m)]: B_a - B_b = value


@dataclass
class Solution:
    keys: list                   # per epoch key lists
    xs: list                     # smoothed states per epoch
    Ps: list                     # smoothed covariances
    xf: list                     # forward (filtered) states
    Pf: list
    used: np.ndarray             # (ne, ns) 0 none, 1 code only, 2 phase only, 3 both
    rejected: np.ndarray         # (ne, ns) bool, any observation rejected by screening/editing
    n_sat: np.ndarray
    epoch_rejected: np.ndarray
    births: dict                 # arc -> birth epoch index
    arcs: object                 # ArcInfo (after editing)
    postfit_rms: dict = field(default_factory=dict)
    nis: float = np.nan
    passes: int = 0
    edit_counts: dict = field(default_factory=dict)
    segments: list = field(default_factory=list)
    chol_failures: int = 0


def _base_keys(cfg):
    k = ["ZWD", "GN", "GE", "CLK[G]"]
    if cfg.mode in ("static", "constrained"):
        k += ["DX", "DY", "DZ"]
    return k


def _sigma2(el, ab):
    s = np.sin(np.maximum(el, np.radians(1.0)))
    return ab[0] ** 2 + ab[1] ** 2 / s ** 2


def _igg3(vbar):
    a = np.abs(vbar)
    k0, k1 = settings.IGG_K0, settings.IGG_K1
    w = np.ones_like(a)
    mid = (a > k0) & (a <= k1)
    w[mid] = (k0 / a[mid]) * ((k1 - a[mid]) / (k1 - k0)) ** 2
    w[a > k1] = 0.0
    return w


def _segments(t_s, usable):
    """Split into processing windows at gaps > window_break_s (TDS § 5.7)."""
    has = usable.any(axis=1)
    idx = np.nonzero(has)[0]
    if len(idx) == 0:
        return []
    segs = []
    start = idx[0]
    for a, b in zip(idx[:-1], idx[1:]):
        if t_s[b] - t_s[a] > settings.WINDOW_BREAK_S:
            segs.append((start, a))
            start = b
    segs.append((start, idx[-1]))
    return segs


def _forward(obs, arcs, cfg, edit, k0, k1, store):
    """Forward EKF over epochs k0..k1 (inclusive). Appends per-epoch records to store."""
    ne, ns = obs.Pif.shape
    base = _base_keys(cfg)
    nb = len(base)
    sp = (settings.SIGMA_PHASE_AB_M[0] * obs.f_if, settings.SIGMA_PHASE_AB_M[1] * obs.f_if)
    sc = (settings.SIGMA_CODE_AB_M[0] * obs.f_if, settings.SIGMA_CODE_AB_M[1] * obs.f_if)
    qz = cfg.sigma_zwd ** 2 / 3600.0
    qg = cfg.sigma_grad ** 2 / 3600.0
    lat, lon, _ = co.ecef_to_geodetic(obs.X0)
    Renu = co.rot_xyz2enu(lat, lon)
    keys = None
    x = P = None
    consec = {}
    births = store["births"]
    for k in range(k0, k1 + 1):
        m = obs.usable[k] & ~edit["removed"][k]
        arc_k = arcs.arc[k]
        js = np.nonzero(m & (arc_k >= 0))[0]
        # ---------------------------------------------------------------- prediction / initialisation
        slant_h = obs.mh[k] * obs.zhd0[k]
        if keys is None:
            keys = list(base)
            x = np.zeros(nb)
            x[0] = obs.zwd0[k]
            P = np.zeros((nb, nb))
            P[0, 0] = settings.INIT_SIGMA_ZWD_M ** 2
            P[1, 1] = P[2, 2] = settings.INIT_SIGMA_GRAD_M ** 2 if cfg.gradients else 1e-12
            P[3, 3] = settings.SIGMA_CLK0_M ** 2
            if nb > 4:
                if cfg.mode == "static":
                    P[4:7, 4:7] = np.eye(3) * settings.INIT_SIGMA_XYZ_STATIC_M ** 2
                else:
                    sn, se, su = cfg.sigma_neu
                    P[4:7, 4:7] = Renu.T @ np.diag([se ** 2, sn ** 2, su ** 2]) @ Renu
            x[3] = _code_clock(obs, k, js, x, keys, slant_h)
            Phi = None
            x_pred, P_pred = x.copy(), P.copy()
            new_arcs = sorted({int(a) for a in arc_k[js]})
            if new_arcs:
                keys, x_pred, P_pred = _augment(keys, x_pred, P_pred, new_arcs, obs, k, arc_k)
                for a in new_arcs:
                    births.setdefault(a, k)
            dt = 0.0
        else:
            dt = obs.t_s[k] - obs.t_s[kprev]
            alive = []
            for key in keys[nb:]:
                a = int(key[4:-1])
                meta = arcs.arc_meta.get(a)
                if cfg.keep_ambiguities or (meta is not None and meta["last"] >= k):
                    alive.append(key)
            new_arcs = sorted({int(a) for a in arc_k[js] if f"AMB[{int(a)}]" not in keys})
            new_keys = base + alive
            Phi = np.zeros((len(new_keys), len(keys)))
            idx_old = {kk: i for i, kk in enumerate(keys)}
            for i, kk in enumerate(new_keys):
                if kk != "CLK[G]":
                    Phi[i, idx_old[kk]] = 1.0
            Q = np.zeros((len(new_keys), len(new_keys)))
            Q[0, 0] = qz * dt
            Q[1, 1] = Q[2, 2] = qg * dt if cfg.gradients else 0.0
            Q[3, 3] = settings.SIGMA_CLK0_M ** 2
            x_pred = Phi @ x
            P_pred = Phi @ P @ Phi.T + Q
            keys_pred = new_keys
            x_pred[3] = _code_clock(obs, k, js, x_pred, keys_pred, slant_h)
            keys, x_pred, P_pred, Phi = _augment(keys_pred, x_pred, P_pred, new_arcs, obs, k, arc_k, Phi)
            for a in new_arcs:
                births.setdefault(a, k)
        # ---------------------------------------------------------------- measurement update
        n = len(keys)
        kidx = {kk: i for i, kk in enumerate(keys)}
        rows_H, rows_v, rows_R, rows_info = [], [], [], []
        for j in js:
            a = int(arc_k[j])
            h_common = np.zeros(n)
            h_common[0] = obs.mw[k, j]
            if cfg.gradients:
                h_common[1] = obs.mg[k, j] * np.cos(obs.az[k, j])
                h_common[2] = obs.mg[k, j] * np.sin(obs.az[k, j])
            h_common[3] = 1.0
            if n > 4 and nb > 4:
                h_common[4:7] = -obs.los[k, j]
            comp = obs.rho[k, j] - obs.dts[k, j] + obs.mh[k, j] * obs.zhd0[k]
            if not edit["code_removed"][k, j]:
                v = obs.Pif[k, j] - (comp + h_common @ x_pred)
                rows_H.append(h_common)
                rows_v.append(v)
                rows_R.append(_sigma2(obs.el[k, j], sc))
                rows_info.append((j, 0, a))
            hp = h_common.copy()
            ai = kidx.get(f"AMB[{a}]")
            if ai is None:
                continue
            hp[ai] = 1.0
            v = obs.Lif[k, j] - (comp + obs.wind[k, j] + hp @ x_pred)
            rows_H.append(hp)
            rows_v.append(v)
            rows_R.append(_sigma2(obs.el[k, j], sp))
            rows_info.append((j, 1, a))
        # single-difference ambiguity pseudo-observations (fixed solution, TDS § 15.6): applied once, at the
        # first epoch at which both ambiguities are in the state
        new_applied = []
        for ci, (aa, bb, val, sig) in enumerate(cfg.sd_constraints):
            if ci in store["applied"]:
                continue
            ia, ib = kidx.get(f"AMB[{aa}]"), kidx.get(f"AMB[{bb}]")
            if ia is None or ib is None:
                continue
            hc = np.zeros(n)
            hc[ia], hc[ib] = 1.0, -1.0
            rows_H.append(hc)
            rows_v.append(val - (x_pred[ia] - x_pred[ib]))
            rows_R.append(sig ** 2)
            rows_info.append((-1, 2, aa))
            new_applied.append(ci)
        x_f, P_f = x_pred.copy(), P_pred.copy()
        ep_rej = False
        used_rows = []
        if rows_H:
            H = np.array(rows_H)
            v = np.array(rows_v)
            Rd = np.array(rows_R)
            S = H @ P_pred @ H.T + np.diag(Rd)
            vbar = v / np.sqrt(np.diag(S))
            w = np.ones(len(v))
            typ = np.array([r[1] for r in rows_info])
            for t in (0, 1):                  # type 2 (constraints) are never down-weighted
                mt = typ == t
                if mt.any():
                    w[mt] = _igg3(vbar[mt])
            ph = typ == 1
            if ph.sum() >= 3 and np.mean(np.abs(vbar[ph]) > settings.IGG_K1) > settings.EPOCH_REJECT_FRACTION:
                ep_rej = True
            else:
                keep = w > 0
                used_rows = [rows_info[i] for i in np.nonzero(keep)[0]]
                for i in np.nonzero(~keep)[0]:
                    j, t, a = rows_info[i]
                    if t == 1:
                        consec[a] = consec.get(a, 0) + 1
                    if t in (0, 1):
                        store["rejected"][k, j] = True
                for (j, t, a) in used_rows:
                    if t == 1:
                        consec[a] = 0
                H, v, Rd = H[keep], v[keep], Rd[keep] / w[keep]
                if len(v):
                    S = H @ P_pred @ H.T + np.diag(Rd)
                    try:
                        cf = cho_factor(S)
                        K = cho_solve(cf, H @ P_pred).T
                    except np.linalg.LinAlgError:
                        store["chol_failures"] += 1
                        K = P_pred @ H.T @ np.linalg.pinv(S)
                    x_f = x_pred + K @ v
                    IKH = np.eye(n) - K @ H
                    P_f = IKH @ P_pred @ IKH.T + K @ np.diag(Rd) @ K.T
                    P_f = 0.5 * (P_f + P_f.T)
                    store["nis"].append((float(v @ np.linalg.solve(S, v)), len(v)))
        if not ep_rej:
            store["applied"].update(new_applied)
        if ep_rej:
            store["epoch_rejected"][k] = True
            store["rejected"][k, js] = True
        # n consecutive rejected phase observations -> new arc (TDS § 5.7)
        for a, c in list(consec.items()):
            if c >= settings.N_CONSEC_REJECT:
                meta = arcs.arc_meta.get(a)
                if meta is not None and k + 1 <= meta["last"]:
                    pp.split_arc(arcs, meta["sat"], k + 1, "CONSEC_REJECT")
                consec[a] = 0
        for (j, t, a) in used_rows:
            if t in (0, 1):
                store["used"][k, j] |= (1 if t == 0 else 2)
        store["n_sat"][k] = len({r[0] for r in used_rows if r[1] == 1})
        store["recs"].append({"k": k, "keys": list(keys), "x_pred": x_pred, "P_pred": P_pred, "Phi": Phi,
                              "x_f": x_f, "P_f": P_f})
        x, P = x_f, P_f
        kprev = k
    return store


def _augment(keys, x, P, new_arcs, obs, k, arc_k, Phi=None):
    """Ambiguity birth (TDS § 5.8): append new AMB states with prior L_IF - P_IF and sigma_amb0."""
    if not new_arcs:
        return (keys, x, P, Phi) if Phi is not None else (keys, x, P)
    add = []
    vals = []
    for a in new_arcs:
        j = np.nonzero(arc_k == a)[0][0]
        add.append(f"AMB[{a}]")
        vals.append(obs.Lif[k, j] - obs.Pif[k, j] - obs.wind[k, j])
    n0 = len(keys)
    n1 = n0 + len(add)
    x2 = np.r_[x, np.array(vals)]
    P2 = np.zeros((n1, n1))
    P2[:n0, :n0] = P
    P2[n0:, n0:] = np.eye(len(add)) * settings.SIGMA_AMB0_M ** 2
    if Phi is not None:
        Phi2 = np.zeros((n1, Phi.shape[1]))
        Phi2[:n0] = Phi
        return keys + add, x2, P2, Phi2
    return keys + add, x2, P2


def _code_clock(obs, k, js, x, keys, slant_h):
    """Robust (median) code-derived receiver clock (m) given the predicted troposphere (TDS § 5.3)."""
    if len(js) == 0:
        return x[3] if len(x) > 3 else 0.0
    r = obs.Pif[k, js] - obs.rho[k, js] + obs.dts[k, js] - slant_h[js] - obs.mw[k, js] * x[0]
    r = r[np.isfinite(r)]
    return float(np.median(r)) if len(r) else 0.0


def _rts(recs):
    """Rauch-Tung-Striebel smoother over one segment (TDS § 5.10), Cholesky solves, no explicit inverse."""
    N = len(recs)
    xs = [None] * N
    Ps = [None] * N
    xs[-1] = recs[-1]["x_f"]
    Ps[-1] = recs[-1]["P_f"]
    for k in range(N - 2, -1, -1):
        nxt = recs[k + 1]
        Phi = nxt["Phi"]
        Pf = recs[k]["P_f"]
        Pp = nxt["P_pred"]
        A = Phi @ Pf                                  # (n1, n0)
        try:
            cf = cho_factor(Pp)
            Ct = cho_solve(cf, A)                     # = Pp^-1 Phi Pf  -> C^T
        except np.linalg.LinAlgError:
            Ct = np.linalg.lstsq(Pp, A, rcond=None)[0]
        Cm = Ct.T
        xs[k] = recs[k]["x_f"] + Cm @ (xs[k + 1] - nxt["x_pred"])
        Pk = Pf + Cm @ (Ps[k + 1] - Pp) @ Cm.T
        Ps[k] = 0.5 * (Pk + Pk.T)
    return xs, Ps


def solve(obs, arcs, cfg, max_passes=None, edit0=None):
    """Forward EKF + RTS smoother with post-fit residual editing (TDS § 13.10, up to 3 passes).

    edit0: observation edit mask from a previous solution (re-used by the AR / fixed passes)."""
    ne, ns = obs.Pif.shape
    max_passes = settings.MAX_EDIT_PASSES if max_passes is None else max_passes
    if edit0 is not None:
        edit = {k: v.copy() for k, v in edit0.items()}
    else:
        edit = {"removed": np.zeros((ne, ns), dtype=bool), "code_removed": np.zeros((ne, ns), dtype=bool)}
    counts = {"phase_outliers_split": 0, "code_outliers_removed": 0}
    sol = None
    for ps in range(max_passes + 1):
        store = {"recs": [], "births": {}, "rejected": np.zeros((ne, ns), dtype=bool),
                 "used": np.zeros((ne, ns), dtype=np.int8), "n_sat": np.zeros(ne, dtype=np.int16),
                 "epoch_rejected": np.zeros(ne, dtype=bool), "nis": [], "chol_failures": 0, "applied": set()}
        segs = _segments(obs.t_s, obs.usable & (arcs.arc >= 0))
        keys_all = [None] * ne
        xs_all = [None] * ne
        Ps_all = [None] * ne
        xf_all = [None] * ne
        Pf_all = [None] * ne
        for (a, b) in segs:
            store["recs"] = []
            _forward(obs, arcs, cfg, edit, a, b, store)
            xs, Ps = _rts(store["recs"])
            for r, xx, PP in zip(store["recs"], xs, Ps):
                k = r["k"]
                keys_all[k], xs_all[k], Ps_all[k] = r["keys"], xx, PP
                xf_all[k], Pf_all[k] = r["x_f"], r["P_f"]
        nis_num = sum(v for v, _ in store["nis"])
        nis_den = sum(d for _, d in store["nis"])
        sol = Solution(keys_all, xs_all, Ps_all, xf_all, Pf_all, store["used"], store["rejected"], store["n_sat"],
                       store["epoch_rejected"], store["births"], arcs, nis=nis_num / nis_den if nis_den else np.nan,
                       passes=ps + 1, segments=segs, chol_failures=store["chol_failures"])
        if ps == max_passes:
            break
        n_new = _edit_residuals(obs, arcs, cfg, sol, edit, counts)
        LOG.info("  pass %d: NIS/dof %.2f, post-fit phase RMS %.1f mm, code RMS %.2f m, new outliers %d",
                 ps + 1, sol.nis, sol.postfit_rms.get("phase", np.nan) * 1e3, sol.postfit_rms.get("code", np.nan),
                 n_new)
        if n_new == 0:
            break
    sol.edit_counts = counts
    sol.edit = edit
    return sol


def residuals(obs, cfg, sol):
    """Post-fit (smoothed) residuals of code and phase (ne, ns) arrays (NaN where not used)."""
    ne, ns = obs.Pif.shape
    rc = np.full((ne, ns), np.nan)
    rp = np.full((ne, ns), np.nan)
    for k in range(ne):
        keys = sol.keys[k]
        if keys is None:
            continue
        x = sol.xs[k]
        kidx = {kk: i for i, kk in enumerate(keys)}
        js = np.nonzero(sol.used[k] > 0)[0]
        for j in js:
            tro = obs.mh[k, j] * obs.zhd0[k] + obs.mw[k, j] * x[0]
            if cfg.gradients:
                tro += obs.mg[k, j] * (x[1] * np.cos(obs.az[k, j]) + x[2] * np.sin(obs.az[k, j]))
            dpos = -obs.los[k, j] @ x[4:7] if "DX" in kidx else 0.0
            comp = obs.rho[k, j] - obs.dts[k, j] + tro + x[3] + dpos
            rc[k, j] = obs.Pif[k, j] - comp
            ai = kidx.get(f"AMB[{int(sol.arcs.arc[k, j])}]")
            if ai is not None:
                rp[k, j] = obs.Lif[k, j] - comp - obs.wind[k, j] - x[ai]
    return rc, rp


def _edit_residuals(obs, arcs, cfg, sol, edit, counts):
    rc, rp = residuals(obs, cfg, sol)
    sp = np.sqrt(_sigma2(obs.el, (settings.SIGMA_PHASE_AB_M[0] * obs.f_if, settings.SIGMA_PHASE_AB_M[1] * obs.f_if)))
    sc = np.sqrt(_sigma2(obs.el, (settings.SIGMA_CODE_AB_M[0] * obs.f_if, settings.SIGMA_CODE_AB_M[1] * obs.f_if)))
    zp = rp / sp
    zc = rc / sc
    fp = np.isfinite(zp)
    fc = np.isfinite(zc)
    s_p = 1.4826 * np.median(np.abs(zp[fp])) if fp.any() else 1.0
    s_c = 1.4826 * np.median(np.abs(zc[fc])) if fc.any() else 1.0
    sol.postfit_rms = {"phase": float(np.sqrt(np.nanmean(rp ** 2))) if fp.any() else np.nan,
                       "code": float(np.sqrt(np.nanmean(rc ** 2))) if fc.any() else np.nan,
                       "phase_scale": float(s_p), "code_scale": float(s_c)}
    bad_p = fp & (np.abs(rp) > np.maximum(settings.POSTFIT_PHASE_K * s_p * sp, settings.POSTFIT_PHASE_MIN_M))
    bad_c = fc & (np.abs(rc) > np.maximum(settings.POSTFIT_CODE_K * s_c * sc, settings.POSTFIT_CODE_MIN_M))
    n = 0
    for k, j in zip(*np.nonzero(bad_p)):
        if edit["removed"][k, j]:
            continue
        edit["removed"][k, j] = True
        # phase outlier inside an arc -> split the arc there (treated as a slip)
        meta = arcs.arc_meta.get(int(arcs.arc[k, j]))
        if meta is not None and k + 1 <= meta["last"] and k > meta["first"]:
            pp.split_arc(arcs, j, k + 1, "POSTFIT_OUTLIER")
        counts["phase_outliers_split"] += 1
        n += 1
    for k, j in zip(*np.nonzero(bad_c & ~bad_p)):
        if not edit["code_removed"][k, j]:
            edit["code_removed"][k, j] = True
            counts["code_outliers_removed"] += 1
            n += 1
    return n


# ================================================================================ extraction
def state_series(sol, key):
    """Smoothed value and sigma of a base state per epoch (NaN where no solution)."""
    ne = len(sol.keys)
    v = np.full(ne, np.nan)
    s = np.full(ne, np.nan)
    vf = np.full(ne, np.nan)
    for k in range(ne):
        if sol.keys[k] is None or key not in sol.keys[k]:
            continue
        i = sol.keys[k].index(key)
        v[k] = sol.xs[k][i]
        s[k] = np.sqrt(max(sol.Ps[k][i, i], 0.0))
        vf[k] = sol.xf[k][i]
    return v, s, vf


def coordinate_estimate(sol, obs):
    """Smoothed coordinate (ARP) and 3x3 covariance from the first epoch of the first segment."""
    for k in range(len(sol.keys)):
        if sol.keys[k] is not None and "DX" in sol.keys[k]:
            i = sol.keys[k].index("DX")
            return obs.X0 + sol.xs[k][i:i + 3], sol.Ps[k][i:i + 3, i:i + 3]
    return None, None


@dataclass
class FiveMin:
    t_utc: np.ndarray            # int ns (UTC calendar), exact 5-min marks
    t_gpst: np.ndarray
    zwd: np.ndarray              # m (estimated absolute ZWD state)
    sig: np.ndarray
    gn: np.ndarray
    ge: np.ndarray
    zhd0: np.ndarray
    clk_ns: np.ndarray
    n_sat: np.ndarray
    n_rej: np.ndarray
    qc: np.ndarray               # uint32 bitmask
    conv: list
    zwd_fwd: np.ndarray = None


def extract_5min(obs, sol, day_start_utc, day_end_utc, flags_global=(), iono=None, extra_epoch_flags=None,
                 edge_ok=(False, False), coord_pull=False):
    """Linear interpolation of smoothed ZWD/gradients to exact 5-min UTC marks (TDS § 6.5) + flags (§ 6.6, § 16.1).

    edge_ok: (start, end) True if overlapping data exist beyond the day boundary on that side (then no EDGE).
    """
    bits = settings.QC_BITS
    zwd, szwd, zf = state_series(sol, "ZWD")
    gn, sgn, _ = state_series(sol, "GN")
    ge, sge, _ = state_series(sol, "GE")
    clk, _, _ = state_series(sol, "CLK[G]")
    step = settings.OUTPUT_INTERVAL_S * ts.NS
    t_utc = np.arange(day_start_utc, day_end_utc, step, dtype=np.int64)
    t_gps = np.array([ts.utc_to_gpst(int(t)) for t in t_utc], dtype=np.int64)
    ep = obs.epochs
    n = len(t_utc)
    out = {k: np.full(n, np.nan) for k in ("zwd", "sig", "gn", "ge", "zhd0", "clk", "zf")}
    nsat = np.zeros(n, dtype=np.int16)
    nrej = np.zeros(n, dtype=np.int16)
    qc = np.zeros(n, dtype=np.uint32)
    conv = []
    have = np.isfinite(zwd)
    # segment edges (for EDGE flag)
    seg_edges = [(obs.epochs[a], obs.epochs[b]) for a, b in sol.segments]
    dt_nom = float(np.median(np.diff(ep))) if len(ep) > 1 else 30e9
    for i, tg in enumerate(t_gps):
        k1 = int(np.searchsorted(ep, tg))
        k0 = k1 - 1
        if k1 < len(ep) and ep[k1] == tg:
            k0 = k1
        ok = 0 <= k0 < len(ep) and 0 <= k1 < len(ep) and have[k0] and have[k1] and \
            (ep[k1] - ep[k0]) <= 1.5 * dt_nom and sol.keys[k0] is not None
        flag = 0
        for g in flags_global:
            if g in bits:
                flag |= 1 << bits[g]
        if not ok:
            flag |= 1 << bits["NO_DATA"]
            qc[i] = flag
            conv.append("NO_DATA")
            continue
        f = 0.0 if k1 == k0 else (tg - ep[k0]) / (ep[k1] - ep[k0])
        for nm, arr in (("zwd", zwd), ("gn", gn), ("ge", ge), ("zf", zf)):
            out[nm][i] = arr[k0] + f * (arr[k1] - arr[k0])
        out["sig"][i] = max(szwd[k0], szwd[k1])
        out["zhd0"][i] = obs.zhd0[k0] + f * (obs.zhd0[k1] - obs.zhd0[k0])
        kn = k0 if f < 0.5 else k1
        out["clk"][i] = clk[kn] / C * 1e9
        nsat[i] = sol.n_sat[kn]
        w = np.abs(ep - tg) <= 150 * ts.NS
        nrej[i] = int(sol.rejected[w].sum())
        nused = int((sol.used[w] > 0).sum())
        if nsat[i] < settings.FEW_SATS:
            flag |= 1 << bits["FEW_SATS"]
        if nused + nrej[i] > 0 and nrej[i] / (nused + nrej[i]) > settings.HIGH_REJECT_FRACTION:
            flag |= 1 << bits["HIGH_REJECT"]
        if iono is not None and iono[w].any():
            flag |= 1 << bits["IONO_ACTIVE"]
        if extra_epoch_flags:
            for name, arr in extra_epoch_flags.items():
                if arr is not None and arr[w].any():
                    flag |= 1 << bits[name]
        if coord_pull:
            flag |= 1 << bits["COORD_PULL"]
        if out["sig"][i] > settings.SIGMA_HIGH_M:
            flag |= 1 << bits["SIGMA_HIGH"]
        # EDGE: within edge_window of the start/end of a continuous segment with no data beyond it
        edge = False
        for si, (a, b) in enumerate(seg_edges):
            if a <= tg <= b:
                first_seg = si == 0
                last_seg = si == len(seg_edges) - 1
                if (tg - a) < settings.EDGE_WINDOW_S * ts.NS and not (first_seg and edge_ok[0]):
                    edge = True
                if (b - tg) < settings.EDGE_WINDOW_S * ts.NS and not (last_seg and edge_ok[1]):
                    edge = True
        # REINIT: > 50 % of active ambiguities born within the previous 1800 s
        keys = sol.keys[kn]
        amb = [kk for kk in keys if kk.startswith("AMB[")]
        young = 0
        for kk in amb:
            b = sol.births.get(int(kk[4:-1]))
            if b is not None and (ep[kn] - ep[b]) < settings.REINIT_WINDOW_S * ts.NS:
                young += 1
        reinit = len(amb) > 0 and young / len(amb) > 0.5 and not edge
        if edge:
            flag |= 1 << bits["EDGE"]
            conv.append("EDGE")
        elif reinit:
            flag |= 1 << bits["REINIT"]
            conv.append("REINIT")
        elif out["sig"][i] <= settings.SIGMA_CONV_M:
            conv.append("CONVERGED")
        else:
            conv.append("NOT_CONVERGED")
        qc[i] = flag
    return FiveMin(t_utc, t_gps, out["zwd"], out["sig"], out["gn"], out["ge"], out["zhd0"], out["clk"], nsat, nrej,
                   qc, conv, out["zf"])


def convergence_time(obs, sol, tol=0.010, hold_s=1800.0):
    """Forward-only convergence time (s from segment start): first epoch after which |fwd - smoothed| < tol
    for >= hold_s (TDS § 6.6 diagnostic)."""
    z, _, zf = state_series(sol, "ZWD")
    if not sol.segments:
        return None
    a, b = sol.segments[0]
    t = obs.t_s
    good = np.abs(zf - z) < tol
    for k in range(a, b + 1):
        if not good[k]:
            continue
        w = (t >= t[k]) & (t <= t[k] + hold_s) & np.isfinite(z)
        if w.any() and good[w].all() and t[b] - t[k] >= hold_s:
            return float(t[k] - t[a])
    return None
