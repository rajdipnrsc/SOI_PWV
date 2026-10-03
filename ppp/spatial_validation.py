"""System 2 spatial validation experiment and supportable resolution (TDS § 24, § 31).

Cross-validation designs on the station residual field r = PWV_GNSS - PWV_bg (the quantity that is mapped):
leave-one-station-out (mapping.loso), spatial block CV (blocks of 100 / 200 km) and network thinning to target
spacings (random and spatially stratified realisations).  From the withheld-station errors the error-distance
relation E(d) = sqrt(E0^2 + (a d^alpha)^2) is fitted, skill over the background is computed and the supportable
grid spacing of § 24.2 is evaluated per region (1-degree tiles).

This module only produces evidence; thresholds and the national grid decision are taken after the experiment and
recorded in docs/decisions.md (TDS § 31 acceptance).  Grid spacing is never called atmospheric resolution.
"""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from . import settings
from .mapping import background_at, gc_km, height_scaling_background, krige, slot_values, to_unit


@dataclass
class SlotData:
    """Station values of one slot: positions, heights, PWV, sigma, background at the station."""
    names: list
    lat: np.ndarray
    lon: np.ndarray
    h: np.ndarray
    pwv: np.ndarray
    sig: np.ndarray
    bg: np.ndarray


def slot_data(stations, slot_ns, bg=None, exclude=(), half_ns=None, min_n=2):
    """Station values for one slot; half_ns/min_n select the temporal window (5/15/30/60-min designs, § 31)."""
    half_ns = int(7.5 * 60) * 10 ** 9 if half_ns is None else half_ns
    vals = [None if s.station in exclude else slot_values(s, slot_ns, half_ns, min_n) for s in stations]
    sts = [s for s, v in zip(stations, vals) if v is not None]
    if not sts:
        return None
    v = np.array([vv for vv in vals if vv is not None])
    lat = np.array([s.lat for s in sts])
    lon = np.array([s.lon for s in sts])
    h = np.array([s.h_orth for s in sts])
    if bg is not None:
        b, _, _ = background_at(bg, lat, lon, h, slot_ns)
    else:
        b = height_scaling_background(v[:, 0], h, h)
    return SlotData([s.station for s in sts], lat, lon, h, v[:, 0], v[:, 1], np.asarray(b, dtype=float))


def predict_withheld(sd, train, test, cp, sigma_repr=None):
    """Predict withheld stations from the training stations (residual kriging + background).

    Returns rows (station, obs, pred, sigma, d_nearest_km, bg_only_pred)."""
    sigma_repr = settings.MAP_SIGMA_BG_REPR_MM if sigma_repr is None else sigma_repr
    if train.sum() < 3 or test.sum() == 0:
        return []
    u = to_unit(sd.lat, sd.lon)
    r = sd.pwv - sd.bg
    s2 = np.sqrt(sd.sig ** 2 + sigma_repr ** 2)
    ti = np.nonzero(train)[0]
    qi = np.nonzero(test)[0]
    rh, var = krige(u[ti], sd.h[ti], r[ti], s2[ti], u[qi], sd.h[qi], cp)
    dn = gc_km(u[qi], u[ti]).min(1)
    return [(sd.names[q], float(sd.pwv[q]), float(sd.bg[q] + rh[k]), float(np.sqrt(var[k] + s2[q] ** 2)),
             float(dn[k]), float(sd.bg[q])) for k, q in enumerate(qi)]


def block_ids(lat, lon, block_km):
    """Square blocks of ~block_km on a latitude-scaled grid (spatial block CV)."""
    dlat = block_km / 111.2
    ilat = np.floor(np.asarray(lat) / dlat).astype(int)
    lat_c = (ilat + 0.5) * dlat
    dlon = block_km / (111.2 * np.cos(np.radians(lat_c)))
    ilon = np.floor(np.asarray(lon) / dlon).astype(int)
    return ilat * 100000 + ilon


def block_cv(sd, cp, block_km):
    """Spatial block CV: each block is withheld in turn and predicted from all other blocks."""
    ids = block_ids(sd.lat, sd.lon, block_km)
    rows = []
    for b in np.unique(ids):
        test = ids == b
        rows += predict_withheld(sd, ~test, test, cp)
    return rows


def thin(lat, lon, spacing_km, rng, mode="random"):
    """Subset of stations with nominal spacing: 'random' = random-order greedy selection with minimum distance
    spacing_km; 'stratified' = one random station per square cell of size spacing_km."""
    n = len(lat)
    if mode == "stratified":
        ids = block_ids(lat, lon, spacing_km)
        keep = np.zeros(n, dtype=bool)
        for b in np.unique(ids):
            keep[rng.choice(np.nonzero(ids == b)[0])] = True
        return keep
    u = to_unit(lat, lon)
    keep = np.zeros(n, dtype=bool)
    for i in rng.permutation(n):
        if not keep.any() or gc_km(u[i:i + 1], u[keep]).min() >= spacing_km:
            keep[i] = True
    return keep


def thinning_cv(sd, cp, spacing_km, n_real=10, mode="random", seed=0):
    """Network thinning (TDS § 31): keep a thinned network, predict all withheld stations; n_real realisations."""
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n_real):
        keep = thin(sd.lat, sd.lon, spacing_km, rng, mode)
        rows += predict_withheld(sd, keep, ~keep, cp)
    return rows


def metrics(rows):
    """RMSE, MAE, bias, r, P95, skill over background, +-1/2 sigma coverage (TDS § 31)."""
    if not rows:
        return {"n": 0}
    o = np.array([r[1] for r in rows])
    p = np.array([r[2] for r in rows])
    s = np.array([r[3] for r in rows])
    b = np.array([r[5] for r in rows])
    e = p - o
    mse_bg = float(np.mean((b - o) ** 2))
    return {"n": len(rows), "bias": float(e.mean()), "rmse": float(np.sqrt(np.mean(e ** 2))),
            "mae": float(np.mean(np.abs(e))), "p95": float(np.percentile(np.abs(e), 95)),
            "r": float(np.corrcoef(o, p)[0, 1]) if len(rows) > 2 and np.std(p) > 0 else float("nan"),
            "rmse_background": float(np.sqrt(mse_bg)),
            "skill": float(1 - np.mean(e ** 2) / mse_bg) if mse_bg > 0 else float("nan"),
            "within_1sigma": float(np.mean(np.abs(e / s) <= 1)), "within_2sigma": float(np.mean(np.abs(e / s) <= 2))}


def fit_error_distance(rows, bins_km=(0, 10, 20, 30, 50, 75, 100, 150, 200, 300)):
    """E(d) = sqrt(E0^2 + (a d^alpha)^2) fitted to the binned RMSE vs distance to the nearest training station
    (TDS § 24.1).  Returns dict with E0, a, alpha, the bins and a callable-free description."""
    d = np.array([r[4] for r in rows])
    e = np.array([r[2] - r[1] for r in rows])
    centres, rmse, cnt = [], [], []
    for lo, hi in zip(bins_km[:-1], bins_km[1:]):
        m = (d >= lo) & (d < hi)
        if m.sum() >= 5:
            centres.append(float(np.median(d[m])))
            rmse.append(float(np.sqrt(np.mean(e[m] ** 2))))
            cnt.append(int(m.sum()))
    out = {"bins_km": centres, "rmse_mm": rmse, "n": cnt, "fitted": False}
    if len(centres) < 3:
        return out
    c, y, w = np.array(centres), np.array(rmse), np.sqrt(np.array(cnt))

    def res(p):
        e0, a, al = p
        return w * (np.sqrt(e0 ** 2 + (a * c ** al) ** 2) - y)
    sol = least_squares(res, [max(y.min(), 0.1), 0.05, 0.7], bounds=([0, 0, 0.1], [50, 10, 2.0]))
    out.update({"E0": float(sol.x[0]), "a": float(sol.x[1]), "alpha": float(sol.x[2]), "fitted": bool(sol.success)})
    return out


def E_of(fit, d_km):
    return float(np.sqrt(fit["E0"] ** 2 + (fit["a"] * d_km ** fit["alpha"]) ** 2))


def nearest_neighbour_km(lat, lon):
    u = to_unit(np.asarray(lat), np.asarray(lon))
    if len(u) < 2:
        return float("inf")
    d = gc_km(u, u)
    np.fill_diagonal(d, np.inf)
    return float(np.median(d.min(1)))


def supportable_spacing(d_nn_km, L_km, fit, skill_at_dnn, lat_c=20.0, candidates=None, k_d=None, s_min=None,
                        e_target=None):
    """Finest candidate spacing satisfying criteria 1-4 of TDS § 24.2 (criterion 5, the incremental benefit with
    block bootstrap, needs RMSE per spacing from point reconstruction and is reported as NOT_EVALUATED here).

    Returns (spacing_deg or None for background-only, per-candidate criteria table)."""
    candidates = candidates or settings.SUPPORT_CANDIDATES_DEG
    k_d = settings.SUPPORT_K_D if k_d is None else k_d
    s_min = settings.SUPPORT_S_MIN if s_min is None else s_min
    e_target = settings.SUPPORT_E_TARGET_MM if e_target is None else e_target
    table = {}
    best = None
    for dg in sorted(candidates, reverse=True):              # coarse -> fine
        dkm = dg * 111.2 * np.sqrt(np.cos(np.radians(lat_c)))  # geometric-mean cell size
        c1 = d_nn_km <= k_d * dkm
        c2 = d_nn_km <= L_km
        c3 = skill_at_dnn >= s_min
        c4 = fit.get("fitted", False) and E_of(fit, dkm / 2) <= e_target
        table[dg] = {"density": bool(c1), "correlation": bool(c2), "skill": bool(c3), "accuracy": bool(c4),
                     "incremental_benefit": "NOT_EVALUATED"}
        if c1 and c2 and c3 and c4:
            best = dg
    return best, table


def support_by_tile(sd, cp, fit, rows_loso, tile_deg=1.0, **kw):
    """Supportable spacing per tile (TDS § 24.4): d_nn of the stations within the tile and its 1-tile margin, skill
    from the LOSO rows of those stations."""
    by_name = {r[0]: r for r in rows_loso}
    out = {}
    ti = np.floor(sd.lat / tile_deg).astype(int)
    tj = np.floor(sd.lon / tile_deg).astype(int)
    for a, b in sorted(set(zip(ti, tj))):
        m = (np.abs(ti - a) <= 1) & (np.abs(tj - b) <= 1)
        if m.sum() < 3:
            continue
        rows = [by_name[n] for n in np.array(sd.names)[m] if n in by_name]
        sk = metrics(rows).get("skill", float("nan")) if rows else float("nan")
        dnn = nearest_neighbour_km(sd.lat[m], sd.lon[m])
        lat_c = (a + 0.5) * tile_deg
        best, table = supportable_spacing(dnn, cp.L_km, fit, sk if np.isfinite(sk) else -1.0, lat_c, **kw)
        out[(float(a * tile_deg), float(b * tile_deg))] = {"d_nn_km": dnn, "skill": sk, "spacing_deg": best,
                                                           "criteria": table}
    return out
