"""PPP-AR layer (Stage 6, TDS § 15). A separate solution layer: float products are never overwritten.

Steps (TDS § 15.2-15.6):
1. OSB-corrected code and phase (all observations, CODE IAR README); satellite PCO not applied to MW.
2. Wide-lane: per-arc Melbourne-Wuebbena mean (cycles) with sigma = s/sqrt(n_eff); receiver WL bias eliminated
   by between-satellite single differences (SD) forming an independent set (spanning forest of overlapping arcs).
3. Narrow-lane: float SD N1 from the AR forward pass (all ambiguities retained to the window end), eq. (15.1).
4. LAMBDA decorrelation + integer least squares (Teunissen 1995; search as in de Jonge & Tiberius 1996),
   validation by bootstrapped success rate (Teunissen 1998) AND ratio test; partial AR.
5. Accepted SD integers -> IF-ambiguity SD constraints (sigma 1 mm) -> forward + RTS re-run -> fixed ZTD;
   float-consistency gate.
"""
import copy
from dataclasses import dataclass, field
from math import erf, sqrt

import numpy as np

from . import estimator as est
from . import log as plog
from . import settings

LOG = plog.get()
C = settings.C_LIGHT

AR_MIN_ARC_S = 1200.0            # min_arc_ar_s [A, TDS § 15.3]
WL_MAX_FRAC = 0.25               # cycles
WL_MAX_SIGMA = 0.10              # cycles
WL_DECOR_S = 300.0               # decorrelation time for n_eff [A]
PS_MIN = 0.999                   # bootstrapped success rate [A]
RATIO_MIN = 3.0                  # [A]
SIGMA_FIX_M = 0.001              # [A, TDS § 15.6]
GATE_MEAN_M = 0.002              # consistency gate [A, TDS § 15.6]
GATE_MAX_M = 0.010


# ================================================================================ LAMBDA
def _ld(Q):
    """Q = L^T diag(D) L with L unit lower triangular (LAMBDA factorisation)."""
    n = Q.shape[0]
    A = Q.copy()
    L = np.zeros((n, n))
    D = np.zeros(n)
    for i in range(n - 1, -1, -1):
        D[i] = A[i, i]
        if D[i] <= 0.0:
            raise np.linalg.LinAlgError("Q not positive definite")
        a = np.sqrt(D[i])
        L[i, :i + 1] = A[i, :i + 1] / a
        for j in range(i):
            A[j, :j + 1] -= L[i, :j + 1] * L[i, j]
        L[i, :i + 1] /= L[i, i]
    return L, D


def _gauss(L, Z, i, j):
    mu = int(np.round(L[i, j]))
    if mu != 0:
        L[i:, j] -= mu * L[i:, i]
        Z[:, j] -= mu * Z[:, i]


def _perm(L, D, j, dl, Z):
    n = L.shape[0]
    eta = D[j] / dl
    lam = D[j + 1] * L[j + 1, j] / dl
    D[j] = eta * D[j + 1]
    D[j + 1] = dl
    for k in range(j):
        a0, a1 = L[j, k], L[j + 1, k]
        L[j, k] = -L[j + 1, j] * a0 + a1
        L[j + 1, k] = eta * a0 + lam * a1
    L[j + 1, j] = lam
    for k in range(j + 2, n):
        L[k, j], L[k, j + 1] = L[k, j + 1], L[k, j]
    Z[:, [j, j + 1]] = Z[:, [j + 1, j]]


def _reduction(L, D, Z):
    n = L.shape[0]
    j = k = n - 2
    while j >= 0:
        if j <= k:
            for i in range(j + 1, n):
                _gauss(L, Z, i, j)
        dl = D[j] + L[j + 1, j] ** 2 * D[j + 1]
        if dl + 1e-6 < D[j + 1]:
            _perm(L, D, j, dl, Z)
            k = j
            j = n - 2
        else:
            j -= 1


def _sgn(x):
    return -1.0 if x <= 0.0 else 1.0


def _search(L, D, zs, m=2, loopmax=100000):
    n = len(zs)
    S = np.zeros((n, n))
    dist = np.zeros(n)
    zb = np.zeros(n)
    z = np.zeros(n)
    step = np.zeros(n)
    zn = np.zeros((n, m))
    s = np.full(m, np.inf)
    nn = imax = 0
    maxdist = 1e99
    k = n - 1
    zb[k] = zs[k]
    z[k] = np.round(zb[k])
    y = zb[k] - z[k]
    step[k] = _sgn(y)
    for _ in range(loopmax):
        newdist = dist[k] + y * y / D[k]
        if newdist < maxdist:
            if k != 0:
                k -= 1
                dist[k] = newdist
                S[k, :k + 1] = S[k + 1, :k + 1] + (z[k + 1] - zb[k + 1]) * L[k + 1, :k + 1]
                zb[k] = zs[k] + S[k, k]
                z[k] = np.round(zb[k])
                y = zb[k] - z[k]
                step[k] = _sgn(y)
            else:
                if nn < m:
                    if nn == 0 or newdist > s[imax]:
                        imax = nn
                    zn[:, nn] = z
                    s[nn] = newdist
                    nn += 1
                else:
                    if newdist < s[imax]:
                        zn[:, imax] = z
                        s[imax] = newdist
                        imax = int(np.argmax(s))
                    maxdist = s[imax]
                z[0] += step[0]
                y = zb[0] - z[0]
                step[0] = -step[0] - _sgn(step[0])
        else:
            if k == n - 1:
                break
            k += 1
            z[k] += step[k]
            y = zb[k] - z[k]
            step[k] = -step[k] - _sgn(step[k])
    o = np.argsort(s)
    return zn[:, o], s[o]


def lambda_ils(a, Q, m=2):
    """Integer least squares by LAMBDA. Returns (candidates (n,m), squared norms (m,), D of decorrelated Q)."""
    a = np.asarray(a, dtype=float)
    L, D = _ld(np.asarray(Q, dtype=float))
    Z = np.eye(len(a))
    _reduction(L, D, Z)
    zs = Z.T @ a
    E, s = _search(L, D, zs, m)
    F = np.linalg.solve(Z.T, E)
    return np.round(F), s, D


def bootstrap_success_rate(D):
    """P_s = prod(2 Phi(1/(2 sqrt(d_i))) - 1) (Teunissen 1998), d_i conditional variances (cycles^2)."""
    phi = lambda x: 0.5 * (1 + erf(x / sqrt(2)))  # noqa: E731
    p = 1.0
    for d in D:
        p *= 2 * phi(1.0 / (2 * sqrt(d))) - 1
    return p


# ================================================================================ wide-lane
@dataclass
class ARResult:
    status: str = "AR_NOT_AVAILABLE"
    metrics: dict = field(default_factory=dict)
    constraints: list = field(default_factory=list)
    solution: object = None
    obs: object = None


def wl_per_arc(sel, arcs, t_s):
    """Per-arc MW mean (cycles), sigma of the mean with n_eff = duration / decorrelation time."""
    lw = C / (sel.f1 - sel.f2)
    mw = ((sel.f1 * sel.L1 - sel.f2 * sel.L2) / (sel.f1 - sel.f2) -
          (sel.f1 * sel.P1 + sel.f2 * sel.P2) / (sel.f1 + sel.f2)) / lw
    out = {}
    for a, meta in arcs.arc_meta.items():
        m = arcs.arc == a
        v = mw[m]
        v = v[np.isfinite(v)]
        if len(v) < 10:
            continue
        ks = np.nonzero(m.any(axis=1))[0]
        dur = t_s[ks[-1]] - t_s[ks[0]]
        neff = max(dur / WL_DECOR_S, 1.0)
        out[a] = {"mu": float(np.mean(v)), "sig": float(np.std(v) / np.sqrt(neff)), "sat": meta["sat"],
                  "first": int(ks[0]), "last": int(ks[-1]), "dur": float(dur)}
    return out


def independent_sd(wl, t_s, min_overlap_s=AR_MIN_ARC_S):
    """Independent between-satellite SD set: maximum-overlap spanning forest over overlapping arcs
    (Blewitt 1989 / Ge et al. 2008 approach). Returns list of (a, b, overlap_s)."""
    ids = sorted(wl)
    edges = []
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            if wl[a]["sat"] == wl[b]["sat"]:
                continue
            k0 = max(wl[a]["first"], wl[b]["first"])
            k1 = min(wl[a]["last"], wl[b]["last"])
            if k1 <= k0:
                continue
            ov = t_s[k1] - t_s[k0]
            if ov >= min_overlap_s:
                edges.append((ov, a, b))
    edges.sort(reverse=True)
    parent = {a: a for a in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    out = []
    for ov, a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
            out.append((a, b, ov))
    return out


def fix_wl(sd, wl):
    """Round SD WL when |frac| <= 0.25 cycle and sigma <= 0.1 cycle (TDS § 15.3)."""
    fixed = {}
    fr = []
    for a, b, ov in sd:
        v = wl[a]["mu"] - wl[b]["mu"]
        s = np.hypot(wl[a]["sig"], wl[b]["sig"])
        n = int(np.round(v))
        fr.append(v - n)
        if abs(v - n) <= WL_MAX_FRAC and s <= WL_MAX_SIGMA:
            fixed[(a, b)] = n
    return fixed, fr


# ================================================================================ narrow-lane
def nl_float(x, P, keys, pairs_wl, f1, f2):
    """Float SD N1 (cycles) and covariance from IF ambiguities (eq. 15.1)."""
    lam_nl = C / (f1 + f2)
    k_wl = C / (f1 - f2) * f2 / (f1 + f2)
    idx = {k: i for i, k in enumerate(keys)}
    rows, vals, used = [], [], []
    for (a, b), nwl in pairs_wl.items():
        ia, ib = idx.get(f"AMB[{a}]"), idx.get(f"AMB[{b}]")
        if ia is None or ib is None:
            continue
        j = np.zeros(len(keys))
        j[ia], j[ib] = 1.0 / lam_nl, -1.0 / lam_nl
        rows.append(j)
        vals.append((x[ia] - x[ib] - k_wl * nwl) / lam_nl)
        used.append((a, b))
    if not rows:
        return None, None, []
    J = np.array(rows)
    return np.array(vals), J @ P @ J.T, used


def partial_ar(N1, Q, order_key=None):
    """LAMBDA with validation; iteratively drop the ambiguity with the largest variance (TDS § 15.5)."""
    keep = list(range(len(N1)))
    while len(keep) >= 4:
        Qs = Q[np.ix_(keep, keep)]
        try:
            F, s, D = lambda_ils(N1[keep], Qs, 2)
        except np.linalg.LinAlgError:
            return None
        ps = bootstrap_success_rate(D)
        ratio = s[1] / s[0] if s[0] > 0 else np.inf
        if ps >= PS_MIN and ratio >= RATIO_MIN:
            return {"idx": keep, "ints": F[:, 0], "ps": ps, "ratio": ratio}
        drop = int(np.argmax(np.diag(Qs)))
        keep.pop(drop)
    return None


# ================================================================================ driver
def run_ar(obs, arcs_float, sol_float, sel_osb, cfg_mode, fm_float, extract):
    """Complete AR layer. Returns ARResult (status, metrics, fixed solution) - float products untouched.

    obs / arcs_float / sol_float: float solution inputs (after residual editing); sel_osb: SelectedObs with
    code AND phase OSBs applied; extract: callable(sol) -> FiveMin for the consistency gate.
    """
    f1, f2 = sel_osb.f1, sel_osb.f2
    res = ARResult()
    t_s = obs.t_s
    arcs = copy.deepcopy(sol_float.arcs)
    wl = wl_per_arc(sel_osb, arcs, t_s)
    sd = independent_sd(wl, t_s)
    wl_fixed, frac = fix_wl(sd, wl)
    m = {"n_amb_total": len(arcs.arc_meta), "n_WL_candidates": len(sd), "n_WL_fixed": len(wl_fixed),
         "wl_frac_rms": float(np.sqrt(np.mean(np.square(frac)))) if frac else None}
    if len(wl_fixed) < 4:
        res.status, res.metrics = "AR_FAILED", {**m, "reason": "fewer than 4 WL-fixed SD ambiguities"}
        return res
    # AR forward pass: all ambiguities retained to the window end (TDS § 5.8)
    cfg = est.EstConfig(mode=cfg_mode, keep_ambiguities=True)
    sol_ar = est.solve(obs, copy.deepcopy(arcs), cfg, max_passes=0, edit0=sol_float.edit)
    last = max(k for k in range(len(sol_ar.keys)) if sol_ar.keys[k] is not None)
    N1, Q, used = nl_float(sol_ar.xf[last], sol_ar.Pf[last], sol_ar.keys[last], wl_fixed, f1, f2)
    m["n_NL_candidates"] = len(used)
    if N1 is None or len(used) < 4:
        res.status, res.metrics = "AR_FAILED", {**m, "reason": "no NL candidates"}
        return res
    sel = partial_ar(N1, Q)
    if sel is None:
        res.status, res.metrics = "AR_FAILED", {**m, "n_NL_fixed": 0, "reason": "validation failed"}
        return res
    lam_nl = C / (f1 + f2)
    k_wl = C / (f1 - f2) * f2 / (f1 + f2)
    cons = []
    for i, n1 in zip(sel["idx"], sel["ints"]):
        a, b = used[i]
        cons.append((a, b, lam_nl * n1 + k_wl * wl_fixed[(a, b)], SIGMA_FIX_M))
    m.update({"n_NL_fixed": len(cons), "fixing_ratio": len(cons) / len(used), "ratio_test_value": sel["ratio"],
              "success_rate": sel["ps"], "n_rejected_fixes": len(used) - len(cons)})
    # fixed re-solution
    cfg_fix = est.EstConfig(mode=cfg_mode, sd_constraints=cons)
    sol_fix = est.solve(obs, copy.deepcopy(arcs), cfg_fix, max_passes=0, edit0=sol_float.edit)
    fm_fix = extract(sol_fix)
    d = (fm_fix.zwd - fm_float.zwd)
    d = d[np.isfinite(d)]
    m["ztd_fixed_minus_float_mean_mm"] = float(np.mean(d) * 1e3) if len(d) else None
    m["ztd_fixed_minus_float_max_mm"] = float(np.max(np.abs(d)) * 1e3) if len(d) else None
    if len(d) == 0 or abs(np.mean(d)) > GATE_MEAN_M or np.max(np.abs(d)) > GATE_MAX_M:
        res.status = "AR_REJECTED_CONSISTENCY"
    else:
        res.status = "FIXED" if len(cons) == len(used) else "PARTIAL"
    res.metrics, res.constraints, res.solution = m, cons, sol_fix
    res.fm = fm_fix
    return res
