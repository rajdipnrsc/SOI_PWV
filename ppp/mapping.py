"""System 2: national PWV/ZWD/ZTD mapping from station products (TDS Part B, § 21-27, § 31).

Baseline method (TDS § 22, F14): NWP/ERA5 background evaluated at station and cell *heights*, GNSS residuals
r = PWV_GNSS - PWV_bg regression-kriged (exponential covariance, REML hyperparameters, local universal kriging),
PWV_grid = PWV_bg(x0, H_cell) + r_hat with sigma^2 = sigma_UK^2 + sigma_bg,repr^2 + sigma_height^2.
Raw ZTD/PWV are never interpolated (TDS § 22.1, § 36.10).

If no ERA5/NWP background is supplied, the "OK on height-normalised PWV" catalogue method (TDS § 23) is used with a
2-km exponential height scaling; such grids are labelled background_source = HEIGHT_SCALING_ONLY (benchmark method,
not the baseline product).

System 2 consumes only the published System-1 interface (station CSV + manifests, TDS § 20).
"""
import csv
import glob
import json
import os
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree

from . import log as plog
from . import pwv as pw
from . import settings
from . import timesys as ts
from .troposphere import saastamoinen_zhd

LOG = plog.get()
R_EARTH_KM = 6371.0
G0 = settings.G0
SUPPORTED_SCHEMA_MAJOR = "1"


# ================================================================================ station ingestion
@dataclass
class StationSeries:
    station: str
    lat: float
    lon: float
    h_ell: float
    h_orth: float
    t_utc: np.ndarray            # int ns UTC calendar
    pwv: np.ndarray              # mm
    sig: np.ndarray              # mm (sigma_PWV_total, after QC inflation)
    ztd: np.ndarray
    tier: str
    manifest_id: str
    solution_type: str


def product_paths(folder, day_str=None, as_of=None, campaign=None):
    """PWV CSV files to use: best_available (or as_of) selection over the manifests (TDS § 19.3); folders without
    manifests fall back to the current *_PWV.csv files."""
    from .products import scan_manifests, select_manifests, selected_file
    sel = select_manifests(scan_manifests(folder), as_of=as_of, campaign=campaign)
    paths = []
    for (st, day), m in sorted(sel.items()):
        if day_str and day.replace("-", "") != day_str.replace("-", ""):
            continue
        p = selected_file(m, "_PWV.csv")
        if p:
            paths.append(p)
    if not sel and as_of is None:
        paths = sorted(glob.glob(os.path.join(folder, f"*_{day_str}_PWV.csv" if day_str else "*_PWV.csv")))
    return paths


def read_station_products(folder, day_str=None, as_of=None, campaign=None):
    """Read System-1 PWV CSV products (TDS § 20). Applies ingestion QC (TDS § 21):
    drop NO_DATA / FEW_SATS / SIGMA_HIGH rows and NOT_CONVERGED; EDGE/REINIT rows get sigma x 2."""
    bits = settings.QC_BITS
    drop_mask = (1 << bits["NO_DATA"]) | (1 << bits["FEW_SATS"]) | (1 << bits["SIGMA_HIGH"])
    out = []
    for path in product_paths(folder, day_str, as_of, campaign):
        with open(path) as fh:
            lines = [L for L in fh if not L.startswith("#")]
        rd = csv.DictReader(lines)
        rows = list(rd)
        if not rows:
            continue
        ver = _schema_of(path)
        if ver and ver.split(".")[0] != SUPPORTED_SCHEMA_MAJOR:
            LOG.warning("SCHEMA_REJECTED: %s has schema %s (supported major %s)", path, ver, SUPPORTED_SCHEMA_MAJOR)
            continue
        t, p, s, z = [], [], [], []
        for r in rows:
            try:
                qc = int(r["QC_flag"])
                if qc & drop_mask or r["convergence_flag"] == "NOT_CONVERGED" or not r["PWV"]:
                    continue
                sig = float(r["sigma_PWV_total"])
                if r["convergence_flag"] in ("EDGE", "REINIT"):
                    sig *= 2.0
                t.append(_iso_ns(r["timestamp_utc"]))
                p.append(float(r["PWV"]))
                s.append(sig)
                z.append(float(r["ZTD"]))
            except (KeyError, ValueError):
                continue
        if not t:
            continue
        r0 = rows[0]
        out.append(StationSeries(r0["station"], float(r0["latitude"]), float(r0["longitude"]), float(r0["height"]),
                                 float(r0.get("H_ant") or r0["height"]), np.array(t, dtype=np.int64), np.array(p),
                                 np.array(s), np.array(z), r0.get("product_tier", ""), r0.get("manifest_id", ""),
                                 r0.get("solution_type", "FLOAT")))
    return out


def _schema_of(path):
    with open(path) as fh:
        for _ in range(3):
            L = fh.readline()
            if "schema_version=" in L:
                return L.split("schema_version=")[1].split(";")[0].strip()
    return None


def _iso_ns(s):
    s = s.rstrip("Z")
    d, t = s.split("T")
    y, m, dd = (int(x) for x in d.split("-"))
    hh, mi, sec = t.split(":")
    return ts.to_ns(y, m, dd, int(hh), int(mi), float(sec))


def slot_values(st, slot_ns, half_ns=int(7.5 * 60) * ts.NS, min_n=2):
    """Station value for a 15-min slot: mean of the 5-min values in [tau - 7.5, tau + 7.5) min, >= 2 of 3;
    sigma_window = mean of member sigmas (full correlation assumed, TDS § 25)."""
    m = (st.t_utc >= slot_ns - half_ns) & (st.t_utc < slot_ns + half_ns)
    if m.sum() < min_n:
        return None
    return float(np.mean(st.pwv[m])), float(np.mean(st.sig[m]))


# ================================================================================ background
@dataclass
class Background:
    """Pressure-level NWP/ERA5 fields on a regular lat/lon grid (TDS § 21-22)."""
    t_utc: np.ndarray            # (nt,) int ns
    lat: np.ndarray              # (ny,) degrees, any order
    lon: np.ndarray              # (nx,)
    p_hpa: np.ndarray            # (nl,) pressure levels
    z: np.ndarray                # (nt, nl, ny, nx) geopotential height (m)
    T: np.ndarray                # (nt, nl, ny, nx) K
    q: np.ndarray                # (nt, nl, ny, nx) kg/kg
    zs: np.ndarray               # (ny, nx) surface height (m) (orography)
    lsm: np.ndarray = None       # (ny, nx) land-sea mask 0..1
    source: str = "ERA5"
    assimilates_gnss: str = "UNKNOWN"


def read_era5(path):
    """Read an ERA5 pressure-level NetCDF (variables z, t, q; optional z_surface/orography, lsm).

    VERIFY V-031: geopotential (m^2/s^2) is converted to geopotential height z/g0 and used as orthometric
    height (difference < 0.1 % of height)."""
    import netCDF4
    with netCDF4.Dataset(path) as nc:
        v = nc.variables
        tname = "valid_time" if "valid_time" in v else "time"
        tv = v[tname]
        times = netCDF4.num2date(tv[:], tv.units, only_use_cftime_datetimes=False)
        t_ns = np.array([ts.to_ns(d.year, d.month, d.day, d.hour, d.minute, d.second) for d in times], dtype=np.int64)
        lev = "pressure_level" if "pressure_level" in v else "level"
        lat = np.array(v["latitude"][:])
        lon = np.array(v["longitude"][:])
        z = np.array(v["z"][:]) / G0
        T = np.array(v["t"][:])
        q = np.array(v["q"][:])
        zs = None
        for nm in ("z_surface", "orography", "zs"):
            if nm in v:
                a = np.array(v[nm][:])
                zs = (a[0] if a.ndim == 3 else a) / (G0 if a.max() > 10000 else 1.0)
        lsm = np.array(v["lsm"][:])[0] if "lsm" in v and v["lsm"].ndim == 3 else (
            np.array(v["lsm"][:]) if "lsm" in v else None)
        p = np.array(v[lev][:], dtype=float)
    o = np.argsort(-p)                                      # surface first (high pressure)
    if zs is None:
        zs = z[0, o[0]]                                     # crude: height of the highest-pressure level
    return Background(t_ns, lat, lon, p[o], z[:, o], T[:, o], q[:, o], zs, lsm, "ERA5")


ERA5_LEVELS = [1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300,
               250, 225, 200, 175, 150, 125, 100, 70, 50]


def download_era5(day_ns, path, area=None):
    """Download ERA5 pressure-level z/t/q (+ surface geopotential and land-sea mask) for one UTC day (+1 h)
    with the Copernicus CDS API. Needs `cdsapi` and a key in ~/.cdsapirc (TDS § 33.6, optional).
    VERIFY V-031: dataset/variable names follow the CDS catalogue at implementation time."""
    import cdsapi
    import netCDF4
    la0, la1, lo0, lo1 = area or settings.MAP_DOMAIN
    y, m, d, *_ = ts.from_ns(day_ns)
    c = cdsapi.Client(quiet=True)
    base = {"product_type": ["reanalysis"], "year": [f"{y:04d}"], "month": [f"{m:02d}"], "day": [f"{d:02d}"],
            "time": [f"{h:02d}:00" for h in range(24)], "area": [la1 + 0.5, lo0 - 0.5, la0 - 0.5, lo1 + 0.5],
            "data_format": "netcdf"}
    pl = path + ".pl.nc"
    sl = path + ".sl.nc"
    c.retrieve("reanalysis-era5-pressure-levels", dict(base, variable=["geopotential", "temperature",
                                                                      "specific_humidity"],
                                                      pressure_level=[str(v) for v in ERA5_LEVELS]), pl)
    c.retrieve("reanalysis-era5-single-levels", dict(base, variable=["geopotential", "land_sea_mask"],
                                                     time=["00:00"]), sl)
    # merge surface geopotential and lsm into the pressure-level file
    with netCDF4.Dataset(sl) as s_, netCDF4.Dataset(pl, "a") as p_:
        zs = np.array(s_["z"][:])
        zs = zs[0] if zs.ndim == 3 else zs
        v = p_.createVariable("z_surface", "f4", ("latitude", "longitude"))
        v[:] = zs
        w = p_.createVariable("lsm", "f4", ("latitude", "longitude"))
        ls = np.array(s_["lsm"][:])
        w[:] = ls[0] if ls.ndim == 3 else ls
    os.replace(pl, path)
    os.remove(sl)
    return path


def column_quantities(z, T, q, p_hpa, H):
    """Vectorised column integrals at height H for profiles (..., nl) ordered surface -> top.

    Returns p_s (hPa), PWV (mm), Tm (K) for each column (TDS eq. 22.1-22.3)."""
    z = np.asarray(z, dtype=float)
    lnp = np.log(np.broadcast_to(p_hpa, z.shape))
    H = np.asarray(H, dtype=float)[..., None] if np.ndim(H) else np.full(z.shape[:-1] + (1,), float(H))
    # pressure at H by log-p interpolation in height; below the lowest level: hypsometric extrapolation
    k = np.clip(np.sum(z < H, axis=-1, keepdims=True) - 1, 0, z.shape[-1] - 2)
    z0 = np.take_along_axis(z, k, -1)
    z1 = np.take_along_axis(z, k + 1, -1)
    l0 = np.take_along_axis(lnp, k, -1)
    l1 = np.take_along_axis(lnp, k + 1, -1)
    lnps = l0 + (H - z0) * (l1 - l0) / (z1 - z0)
    below = H < z[..., :1]
    if below.any():
        Tv = T[..., :1] * (1 + 0.608 * q[..., :1])
        lnps = np.where(below, lnp[..., :1] + G0 * (z[..., :1] - H) / (settings.RD_DRY * Tv), lnps)
    ps = np.exp(lnps)[..., 0]
    # PWV = (1/g) int_0^{p_s} q dp (Pa) -> kg/m2 = mm; integrate levels above H plus partial segment
    p_pa = np.broadcast_to(p_hpa, z.shape) * 100.0
    above = z >= H
    qs = np.where(below[..., 0], q[..., 0], _interp_at(z, q, H[..., 0]))
    Ts = np.where(below[..., 0], T[..., 0] + 0.0065 * (z[..., 0] - H[..., 0]), _interp_at(z, T, H[..., 0]))
    # stack (p_s, q_s) as the bottom node then the levels above H
    pp_ = np.where(above, p_pa, np.nan)
    qq_ = np.where(above, q, np.nan)
    P = np.concatenate([ps[..., None] * 100.0, pp_], -1)
    Q = np.concatenate([qs[..., None], qq_], -1)
    order = np.argsort(-np.nan_to_num(P, nan=-1.0), axis=-1)
    P = np.take_along_axis(P, order, -1)
    Q = np.take_along_axis(Q, order, -1)
    dp = -np.diff(P, axis=-1)
    qm = 0.5 * (Q[..., 1:] + Q[..., :-1])
    pwv = np.nansum(np.where(np.isfinite(dp) & np.isfinite(qm), qm * dp, 0.0), axis=-1) / G0
    # Tm = int e/T dz / int e/T^2 dz from H upward
    e = q * np.broadcast_to(p_hpa, z.shape) / (0.622 + 0.378 * q)
    es = qs * ps / (0.622 + 0.378 * qs)
    Z = np.concatenate([H, np.where(above, z, np.nan)], -1)
    E = np.concatenate([es[..., None], np.where(above, e, np.nan)], -1)
    TT = np.concatenate([Ts[..., None], np.where(above, T, np.nan)], -1)
    o2 = np.argsort(np.nan_to_num(Z, nan=1e9), axis=-1)
    Z, E, TT = (np.take_along_axis(a, o2, -1) for a in (Z, E, TT))
    dz = np.diff(Z, axis=-1)
    f1 = 0.5 * (E[..., 1:] / TT[..., 1:] + E[..., :-1] / TT[..., :-1])
    f2 = 0.5 * (E[..., 1:] / TT[..., 1:] ** 2 + E[..., :-1] / TT[..., :-1] ** 2)
    ok = np.isfinite(dz) & np.isfinite(f1)
    tm = np.nansum(np.where(ok, f1 * dz, 0), -1) / np.nansum(np.where(ok, f2 * dz, 0), -1)
    return ps, pwv, tm


def _interp_at(z, v, H):
    k = np.clip(np.sum(z < H[..., None], axis=-1, keepdims=True) - 1, 0, z.shape[-1] - 2)
    z0, z1 = np.take_along_axis(z, k, -1)[..., 0], np.take_along_axis(z, k + 1, -1)[..., 0]
    v0, v1 = np.take_along_axis(v, k, -1)[..., 0], np.take_along_axis(v, k + 1, -1)[..., 0]
    return v0 + (H - z0) * (v1 - v0) / (z1 - z0)


def _bilinear_weights(bg, lat, lon):
    la, lo = bg.lat, bg.lon
    ia = np.interp(lat, la if la[0] < la[-1] else la[::-1], np.arange(len(la)) if la[0] < la[-1]
                   else np.arange(len(la))[::-1])
    io = np.interp(lon, lo, np.arange(len(lo)))
    i0 = np.clip(np.floor(ia).astype(int), 0, len(la) - 2)
    j0 = np.clip(np.floor(io).astype(int), 0, len(lo) - 2)
    return i0, j0, ia - i0, io - j0


def background_at(bg, lat, lon, H, t_ns):
    """PWV_bg (mm), p_s (hPa), Tm (K) at arbitrary points and heights (bilinear profiles, linear in time)."""
    lat, lon, H = np.atleast_1d(lat), np.atleast_1d(lon), np.atleast_1d(H)
    i0, j0, fa, fo = _bilinear_weights(bg, lat, lon)
    tt = bg.t_utc.astype(float)
    it = int(np.clip(np.searchsorted(tt, t_ns) - 1, 0, max(len(tt) - 2, 0)))
    wt = 0.0 if len(tt) == 1 else float(np.clip((t_ns - tt[it]) / (tt[it + 1] - tt[it]), 0, 1))
    res = []
    for ti, w in ((it, 1 - wt), (min(it + 1, len(tt) - 1), wt)):
        prof = {}
        for nm in ("z", "T", "q"):
            a = getattr(bg, nm)[ti]                       # (nl, ny, nx)
            prof[nm] = ((1 - fa) * (1 - fo))[:, None] * a[:, i0, j0].T + (fa * (1 - fo))[:, None] * a[:, i0 + 1, j0].T + \
                ((1 - fa) * fo)[:, None] * a[:, i0, j0 + 1].T + (fa * fo)[:, None] * a[:, i0 + 1, j0 + 1].T
        ps, pwv, tm = column_quantities(prof["z"], prof["T"], prof["q"], bg.p_hpa, H)
        res.append((w, ps, pwv, tm))
    ps = sum(w * a for w, a, _, _ in res)
    pwv = sum(w * b for w, _, b, _ in res)
    tm = sum(w * c for w, _, _, c in res)
    return pwv, ps, tm


def height_scaling_background(stations_pwv, H_st, H_q, scale_h=2000.0):
    """Benchmark 'background' (no NWP): PWV0 * exp(-H / 2 km) with PWV0 from the median normalised value."""
    p0 = np.nanmedian(np.asarray(stations_pwv) * np.exp(np.asarray(H_st) / scale_h))
    return p0 * np.exp(-np.asarray(H_q) / scale_h)


# ================================================================================ covariance / kriging
def to_unit(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def gc_km(u1, u2):
    """Great-circle distance matrix (km) between unit-vector sets."""
    d = np.clip(u1 @ u2.T, -1, 1)
    return R_EARTH_KM * np.arccos(d)


def cov_exp(d_km, dh_m, sigma_s, L_km, L_h_m=np.inf):
    c = sigma_s ** 2 * np.exp(-d_km / L_km)
    if np.isfinite(L_h_m):
        c = c * np.exp(-np.abs(dh_m) / L_h_m)
    return c


@dataclass
class CovParams:
    sigma_s: float = settings.MAP_DEFAULT_COV["sigma_s_mm"]
    L_km: float = settings.MAP_DEFAULT_COV["L_km"]
    nugget: float = settings.MAP_DEFAULT_COV["nugget_mm"]
    L_h_m: float = np.inf
    fitted: bool = False
    n_slots: int = 0


def fit_reml(slots, max_n=300, seed=0):
    """REML fit (constant mean) of sigma_s, L, nugget pooled over slots (TDS § 22.4).

    slots: list of (unit_vectors (n,3), residuals (n,), sigma_i (n,)). L bounded to [20, 500] km."""
    rng = np.random.default_rng(seed)
    use = []
    for u, r, s in slots:
        if len(r) < 8:
            continue
        if len(r) > max_n:
            idx = rng.choice(len(r), max_n, replace=False)
            u, r, s = u[idx], r[idx], s[idx]
        use.append((gc_km(u, u), r, s))
    if not use:
        return CovParams()
    lo, hi = settings.MAP_L_BOUNDS_KM

    def nll(th):
        ss, L, tau = np.exp(th[0]), lo + (hi - lo) / (1 + np.exp(-th[1])), np.exp(th[2])
        tot = 0.0
        for D, r, s in use:
            Cm = cov_exp(D, 0, ss, L) + np.diag(tau ** 2 + s ** 2)
            try:
                Lc = np.linalg.cholesky(Cm)
            except np.linalg.LinAlgError:
                return 1e12
            F = np.ones(len(r))
            Ci_r = np.linalg.solve(Lc.T, np.linalg.solve(Lc, r))
            Ci_F = np.linalg.solve(Lc.T, np.linalg.solve(Lc, F))
            FtCiF = float(F @ Ci_F)
            beta = float(F @ Ci_r) / FtCiF
            res = r - beta
            Ci_res = Ci_r - beta * Ci_F
            tot += 2 * np.sum(np.log(np.diag(Lc))) + res @ Ci_res + np.log(FtCiF)
        return 0.5 * tot
    x0 = np.array([np.log(2.0), 0.0, np.log(0.5)])
    r = minimize(nll, x0, method="Nelder-Mead", options={"maxiter": 400, "xatol": 1e-3, "fatol": 1e-3})
    ss, L, tau = np.exp(r.x[0]), lo + (hi - lo) / (1 + np.exp(-r.x[1])), np.exp(r.x[2])
    return CovParams(float(ss), float(L), float(tau), fitted=True, n_slots=len(use))


def krige(u_st, h_st, r, s_st, u_q, h_q, cp, k_max=None, r_max_km=None, drift=None, drift_q=None, batch=2000):
    """Local universal kriging of residuals (TDS eq. 22.8), batched over query points.

    Neighbourhood: nearest k_max stations within R_max = max(3L, 300 km); stations beyond R_max get zero
    covariance with the target and are decoupled (weight 0). Constant drift always; optional extra drift
    covariates (n, p) / (m, p). Returns r_hat, sigma2_UK (prediction variance of the smooth field)."""
    k_max = k_max or settings.MAP_K_MAX
    r_max = r_max_km or max(3 * cp.L_km, 300.0)
    n, m = len(r), len(u_q)
    rhat = np.zeros(m)
    c0 = cp.sigma_s ** 2
    var = np.full(m, c0 + cp.nugget ** 2)
    if n == 0 or m == 0:
        return rhat, var
    k = min(k_max, n)
    tree = cKDTree(u_st)
    dist, idx = tree.query(u_q, k=k)
    if k == 1:
        dist, idx = dist[:, None], idx[:, None]
    dkm = 2 * np.arcsin(np.clip(dist / 2, 0, 1)) * R_EARTH_KM
    inside = dkm <= r_max
    p = 1 + (drift.shape[1] if drift is not None else 0)
    for a0 in range(0, m, batch):
        sl = slice(a0, min(a0 + batch, m))
        I = idx[sl]
        ins = inside[sl]
        B = I.shape[0]
        us = u_st[I]                                            # (B, k, 3)
        D = R_EARTH_KM * np.arccos(np.clip(np.einsum("bik,bjk->bij", us, us), -1, 1))
        dh = h_st[I][:, :, None] - h_st[I][:, None, :]
        Cm = cov_exp(D, dh, cp.sigma_s, cp.L_km, cp.L_h_m)
        Cm = Cm * (ins[:, :, None] & ins[:, None, :])
        diag = np.where(ins, cp.nugget ** 2 + s_st[I] ** 2, 1e12)
        Cm[:, np.arange(k), np.arange(k)] = np.where(ins, Cm[:, np.arange(k), np.arange(k)] + diag, diag)
        c = cov_exp(dkm[sl], h_q[sl][:, None] - h_st[I], cp.sigma_s, cp.L_km, cp.L_h_m) * ins
        F = np.ones((B, k, p))
        f0 = np.ones((B, p))
        if drift is not None:
            F[:, :, 1:] = drift[I]
            f0[:, 1:] = drift_q[sl]
        F = F * ins[:, :, None]
        A = np.zeros((B, k + p, k + p))
        A[:, :k, :k] = Cm
        A[:, :k, k:] = F
        A[:, k:, :k] = np.transpose(F, (0, 2, 1))
        rhs = np.concatenate([c, f0], axis=1)
        enough = ins.sum(1) >= p + 1
        sol = np.zeros((B, k + p))
        if enough.any():
            sol[enough] = np.linalg.solve(A[enough], rhs[enough][..., None])[..., 0]
        w, mu = sol[:, :k], sol[:, k:]
        rh = np.einsum("bk,bk->b", w, r[I])
        vv = c0 - np.einsum("bk,bk->b", w, c) - np.einsum("bp,bp->b", mu, f0)
        rhat[sl] = np.where(enough, rh, 0.0)
        var[sl] = np.where(enough, np.maximum(vv, 0.0), c0 + cp.nugget ** 2)
    return rhat, var


# ================================================================================ grid, classes
@dataclass
class Grid:
    name: str
    step: float
    lat: np.ndarray
    lon: np.ndarray
    h: np.ndarray                # (ny, nx) cell mean orthometric height (m)
    h_std: np.ndarray
    land: np.ndarray             # (ny, nx) bool
    dem_source: str = ""


def make_grid(name="G025", domain=None, dem=None, bg=None):
    """National grid with ERA5-coincident cell centres (multiples of the step, decisions D-010)."""
    step = settings.MAP_GRIDS[name]
    la0, la1, lo0, lo1 = domain or settings.MAP_DOMAIN
    lat = np.arange(np.ceil(la0 / step) * step, la1 + 1e-9, step)
    lon = np.arange(np.ceil(lo0 / step) * step, lo1 + 1e-9, step)
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    h = np.zeros(LA.shape)
    hs = np.zeros(LA.shape)
    land = np.ones(LA.shape, dtype=bool)
    src = "none (flat 0 m; supply a DEM)"
    if dem is not None:
        dl, do, dh = dem
        h, hs = _aggregate_dem(dl, do, dh, lat, lon, step)
        src = "user DEM"
    elif bg is not None:
        i0, j0, fa, fo = _bilinear_weights(bg, LA.ravel(), LO.ravel())
        zs = bg.zs
        h = ((1 - fa) * (1 - fo) * zs[i0, j0] + fa * (1 - fo) * zs[i0 + 1, j0] + (1 - fa) * fo * zs[i0, j0 + 1] +
             fa * fo * zs[i0 + 1, j0 + 1]).reshape(LA.shape)
        src = f"{bg.source} orography (coarse; DEM_FROM_BACKGROUND)"
        if bg.lsm is not None:
            lsm = bg.lsm
            land = ((1 - fa) * (1 - fo) * lsm[i0, j0] + fa * (1 - fo) * lsm[i0 + 1, j0] + (1 - fa) * fo * lsm[i0, j0 + 1]
                    + fa * fo * lsm[i0 + 1, j0 + 1]).reshape(LA.shape) >= 0.5
    return Grid(name, step, lat, lon, h, hs, land, src)


def _aggregate_dem(dl, do, dh, lat, lon, step):
    h = np.full((len(lat), len(lon)), np.nan)
    hs = np.zeros_like(h)
    for i, la in enumerate(lat):
        mi = np.abs(dl - la) <= step / 2
        if not mi.any():
            continue
        for j, lo in enumerate(lon):
            mj = np.abs(do - lo) <= step / 2
            blk = dh[np.ix_(mi, mj)]
            if blk.size:
                h[i, j] = np.nanmean(blk)
                hs[i, j] = np.nanstd(blk)
    return np.nan_to_num(h), hs


def coverage_class(d1_km, n50, n100, dh_m, info, land):
    """Coverage classes 0-4 (TDS § 26, provisional thresholds from settings)."""
    c1d, c1n, c1h, c1i = settings.MAP_CLASS_RULES["c1"]
    c2d, c2n, c2i = settings.MAP_CLASS_RULES["c2"]
    c3d, c3i = settings.MAP_CLASS_RULES["c3"]
    cls = np.full(np.shape(d1_km), 4, dtype=np.int8)
    cls[(d1_km <= c3d) & (info >= c3i)] = 3
    cls[(d1_km <= c2d) & (n100 >= c2n) & (info >= c2i)] = 2
    cls[(d1_km <= c1d) & (n50 >= c1n) & (np.abs(dh_m) <= c1h) & (info >= c1i)] = 1
    # height-difference demotions
    cls[(cls == 1) & (np.abs(dh_m) > 300)] = 2
    cls[(cls == 2) & (np.abs(dh_m) > 800)] = 3
    cls[~land] = 0
    return cls


# ================================================================================ one slot
@dataclass
class SlotResult:
    pwv: np.ndarray
    sigma: np.ndarray
    bg: np.ndarray
    resid: np.ndarray
    ztd: np.ndarray
    zwd: np.ndarray
    zhd: np.ndarray
    sig_ztd: np.ndarray
    info: np.ndarray
    cls: np.ndarray
    d1: np.ndarray
    n50: np.ndarray
    n100: np.ndarray
    res_eff: np.ndarray
    ps: np.ndarray
    tm: np.ndarray
    avail: np.ndarray            # (n_station,) 0/1
    n_used: int = 0
    outliers: list = field(default_factory=list)


SLOT_WINDOWS = {"15M": (int(7.5 * 60) * ts.NS, 2),        # [tau - 7.5, tau + 7.5) min, >= 2 of 3 (TDS § 25)
                "01H": (30 * 60 * ts.NS, 8)}                # hourly archive: [t - 30, t + 30) min, >= 8 of 12


def map_slot(grid, stations, slot_ns, cp, bg=None, sigma_repr=settings.MAP_SIGMA_BG_REPR_MM, res_bg_km=100.0,
             window="15M", exclude=()):
    """Build all layers for one 15-min (or hourly) slot (TDS § 22.3-22.6, § 24.3, § 26).
    exclude: blacklisted station names (§ 21 item 2)."""
    half, min_n = SLOT_WINDOWS[window]
    vals = [None if s.station in exclude else slot_values(s, slot_ns, half, min_n) for s in stations]
    avail = np.array([v is not None for v in vals], dtype=np.int8)
    sts = [s for s, v in zip(stations, vals) if v is not None]
    v = np.array([vv for vv in vals if vv is not None]) if any(avail) else np.zeros((0, 2))
    LA, LO = np.meshgrid(grid.lat, grid.lon, indexing="ij")
    shape = LA.shape
    uq = to_unit(LA.ravel(), LO.ravel())
    hq = grid.h.ravel()
    if bg is not None:
        bg_q, ps_q, tm_q = background_at(bg, LA.ravel(), LO.ravel(), hq, slot_ns)
    else:
        bg_q = ps_q = tm_q = None
    if len(sts):
        lat_s = np.array([s.lat for s in sts])
        lon_s = np.array([s.lon for s in sts])
        h_s = np.array([s.h_orth for s in sts])
        us = to_unit(lat_s, lon_s)
        if bg is not None:
            bg_s, _, _ = background_at(bg, lat_s, lon_s, h_s, slot_ns)
        else:
            bg_s = height_scaling_background(v[:, 0], h_s, h_s)
            bg_q = height_scaling_background(v[:, 0], h_s, hq)
        r = v[:, 0] - bg_s
        s2 = np.sqrt(v[:, 1] ** 2 + sigma_repr ** 2)
        # per-slot spatial consistency (leave-one-out standardised residual > 5 -> outlier, iterate once)
        keep = np.ones(len(r), dtype=bool)
        outl = []
        for _ in range(2):
            idx = np.nonzero(keep)[0]
            if len(idx) < 5:
                break
            new_out = []
            for ii in idx:
                oth = idx[idx != ii]
                rh, vv = krige(us[oth], h_s[oth], r[oth], s2[oth], us[ii:ii + 1], h_s[ii:ii + 1], cp)
                z = abs(r[ii] - rh[0]) / np.sqrt(vv[0] + cp.nugget ** 2 + s2[ii] ** 2)
                if z > 5:
                    new_out.append(ii)
            if not new_out:
                break
            keep[new_out] = False
            outl += [sts[i].station for i in new_out]
        rhat, var = krige(us[keep], h_s[keep], r[keep], s2[keep], uq, hq, cp)
        tree = cKDTree(us[keep])
        d, nn = tree.query(uq, k=1)
        d1 = 2 * np.arcsin(np.clip(d / 2, 0, 1)) * R_EARTH_KM
        dh = hq - h_s[keep][nn]
        dmat = gc_km(uq, us[keep])
        n50 = (dmat <= 50).sum(1)
        n100 = (dmat <= 100).sum(1)
        n_used = int(keep.sum())
    else:
        rhat = np.zeros(len(uq))
        var = np.full(len(uq), cp.sigma_s ** 2 + cp.nugget ** 2)
        d1 = np.full(len(uq), np.inf)
        dh = np.zeros(len(uq))
        n50 = n100 = np.zeros(len(uq), dtype=int)
        n_used, outl = 0, []
        if bg_q is None:
            bg_q = np.full(len(uq), np.nan)
    prior = cp.sigma_s ** 2 + cp.nugget ** 2
    info = np.clip(1 - var / prior, 0, 1)
    s_height = 0.0005 * np.abs(dh) + 0.5 * 0.0005 * grid.h_std.ravel()        # [E] coefficients, refined by CV
    pwv = bg_q + rhat
    sig = np.sqrt(var + sigma_repr ** 2 + s_height ** 2)
    cls = coverage_class(d1, n50, n100, dh, info, grid.land.ravel())
    if tm_q is None:
        tm_q = np.full(len(uq), 280.0)
    if ps_q is None:
        ps_q = 1013.25 * (1 - 2.2557e-5 * hq) ** 5.2568
    Pi = pw.pi_factor(tm_q)
    zwd = pwv / Pi
    zhd = saastamoinen_zhd(ps_q, np.radians(LA.ravel()), hq) * 1e3
    ztd = zhd + zwd
    s_ztd = np.sqrt((sig / Pi) ** 2 + (0.0022768 * 1e3 * 1.0) ** 2)            # 1 hPa background pressure [A]
    res_eff = np.where(info >= 0.5, 2 * d1, res_bg_km)
    fill = cls == 4
    for a in (pwv, sig, zwd, ztd, s_ztd):
        a[fill | (cls == 0)] = np.nan
    rs = lambda a: np.asarray(a, dtype=float).reshape(shape)  # noqa: E731
    return SlotResult(rs(pwv), rs(sig), rs(bg_q), rs(rhat), rs(ztd), rs(zwd), rs(zhd), rs(s_ztd), rs(info),
                      cls.reshape(shape), rs(d1), n50.reshape(shape), n100.reshape(shape), rs(res_eff), rs(ps_q),
                      rs(tm_q), avail, n_used, outl)


# ================================================================================ cross-validation
def loso(stations, slot_ns, cp, bg=None, exclude=()):
    """Leave-one-station-out residual CV for one slot (TDS § 31): returns list of (station, obs, pred, sigma,
    d_nearest_km)."""
    vals = [None if s.station in exclude else slot_values(s, slot_ns) for s in stations]
    sts = [s for s, v in zip(stations, vals) if v is not None]
    v = np.array([vv for vv in vals if vv is not None])
    if len(sts) < 5:
        return []
    lat = np.array([s.lat for s in sts])
    lon = np.array([s.lon for s in sts])
    h = np.array([s.h_orth for s in sts])
    u = to_unit(lat, lon)
    if bg is not None:
        b, _, _ = background_at(bg, lat, lon, h, slot_ns)
    else:
        b = height_scaling_background(v[:, 0], h, h)
    r = v[:, 0] - b
    s2 = np.sqrt(v[:, 1] ** 2 + settings.MAP_SIGMA_BG_REPR_MM ** 2)
    out = []
    for i in range(len(sts)):
        o = np.arange(len(sts)) != i
        rh, var = krige(u[o], h[o], r[o], s2[o], u[i:i + 1], h[i:i + 1], cp)
        dn = gc_km(u[i:i + 1], u[o]).min()
        out.append((sts[i].station, float(v[i, 0]), float(b[i] + rh[0]), float(np.sqrt(var[0] + s2[i] ** 2)),
                    float(dn)))
    return out


def cv_metrics(rows):
    if not rows:
        return {}
    o = np.array([r[1] for r in rows])
    p = np.array([r[2] for r in rows])
    s = np.array([r[3] for r in rows])
    e = p - o
    z = e / s
    return {"n": len(rows), "bias": float(np.mean(e)), "rmse": float(np.sqrt(np.mean(e ** 2))),
            "mae": float(np.mean(np.abs(e))), "p95": float(np.percentile(np.abs(e), 95)),
            "r": float(np.corrcoef(o, p)[0, 1]) if len(rows) > 2 else np.nan,
            "within_1sigma": float(np.mean(np.abs(z) <= 1)), "within_2sigma": float(np.mean(np.abs(z) <= 2))}


# ================================================================================ station grey/blacklist
def station_daily_residuals(stations, day_ns, bg=None, step_h=1):
    """Median (PWV_GNSS - PWV_bg) per station for one day from hourly 15-min slot values (input to § 21 item 2)."""
    res = {s.station: [] for s in stations}
    for k in range(0, 24, step_h):
        sl = day_ns + k * 3600 * ts.NS
        vals = [slot_values(s, sl) for s in stations]
        sts = [s for s, v in zip(stations, vals) if v is not None]
        if not sts:
            continue
        v = np.array([vv[0] for vv in vals if vv is not None])
        h = np.array([s.h_orth for s in sts])
        if bg is not None:
            b, _, _ = background_at(bg, np.array([s.lat for s in sts]), np.array([s.lon for s in sts]), h, sl)
        else:
            b = height_scaling_background(v, h, h)
        for s, r in zip(sts, v - b):
            res[s.station].append(float(r))
    return {k: (float(np.median(v)), len(v)) for k, v in res.items() if v}


def update_residual_history(path, day_tag, daily, stations):
    """Plain CSV history station,day,lat,lon,median_resid_mm,n (one row per station-day; the day is replaced)."""
    rows = []
    if os.path.exists(path):
        with open(path) as fh:
            rows = [r for r in csv.DictReader(fh) if r["day"] != day_tag]
    pos = {s.station: (s.lat, s.lon) for s in stations}
    for st, (med, n) in sorted(daily.items()):
        rows.append({"station": st, "day": day_tag, "lat": f"{pos[st][0]:.5f}", "lon": f"{pos[st][1]:.5f}",
                     "median_resid_mm": f"{med:.3f}", "n": n})
    tmp = path + f".{os.getpid()}.tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, ["station", "day", "lat", "lon", "median_resid_mm", "n"])
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)
    return rows


def station_lists(history_rows, day_tag, days=None, radius_km=None, grey_mm=None, black_mm=None, min_days=10):
    """Greylist/blacklist (TDS § 21 item 2): the station's 30-day median residual minus the median of its
    neighbours' (within 150 km) medians; > 2 mm greylist (review, still used), > 4 mm blacklist (excluded) [A].
    Uses the days before day_tag only.  Returns (grey, black, table)."""
    days = days or settings.MAP_LIST_DAYS
    radius_km = radius_km or settings.MAP_LIST_RADIUS_KM
    grey_mm = grey_mm or settings.MAP_GREY_MM
    black_mm = black_mm or settings.MAP_BLACK_MM
    y, d = int(day_tag[:4]), int(day_tag[4:])
    d_ns = ts.to_ns(y, 1, 1) + (d - 1) * ts.DAY_NS
    win = set()
    for k in range(1, days + 1):
        yy, dd = ts.year_doy(d_ns - k * ts.DAY_NS)
        win.add(f"{yy:04d}{dd:03d}")
    per, pos = {}, {}
    for r in history_rows:
        if r["day"] in win:
            per.setdefault(r["station"], []).append(float(r["median_resid_mm"]))
            pos[r["station"]] = (float(r["lat"]), float(r["lon"]))
    med = {k: float(np.median(v)) for k, v in per.items() if len(v) >= min_days}
    names = sorted(med)
    grey, black, table = [], [], {}
    if len(names) < 2:
        return grey, black, table
    u = to_unit(np.array([pos[n][0] for n in names]), np.array([pos[n][1] for n in names]))
    dm = gc_km(u, u)
    for i, n in enumerate(names):
        nb = [names[j] for j in np.nonzero((dm[i] <= radius_km) & (np.arange(len(names)) != i))[0]]
        if not nb:
            continue
        dev = med[n] - float(np.median([med[m] for m in nb]))
        table[n] = {"median_mm": med[n], "deviation_mm": dev, "neighbours": len(nb)}
        if abs(dev) > black_mm:
            black.append(n)
        elif abs(dev) > grey_mm:
            grey.append(n)
    return grey, black, table


# ================================================================================ aggregation (TDS § 27.1)
def _agg_weights(fine, coarse_step):
    """Area weights of fine cells inside each coarse cell (1-D): overlap length / fine step."""
    w = np.zeros((len(coarse_step[0]), len(fine)))
    for i, c in enumerate(coarse_step[0]):
        lo, hi = c - coarse_step[1] / 2, c + coarse_step[1] / 2
        flo, fhi = fine - coarse_step[2] / 2, fine + coarse_step[2] / 2
        w[i] = np.clip(np.minimum(hi, fhi) - np.maximum(lo, flo), 0, None) / coarse_step[2]
    return w


def aggregate(grid, results, name="G050", cp=None):
    """Variance-aware aggregation of a fine grid to a coarse one (0.25 deg -> 0.5 deg, TDS § 27.1).

    Values: area-weighted means of the fine cells overlapping the coarse cell (ERA5-coincident centres on both grids,
    so edge cells count with weight 1/2).  sigma: sqrt(w^T C w)/sum(w) with C_ij = rho_ij sigma_i sigma_j and
    rho_ij = exp(-d_ij/L) from the fitted residual covariance (errors of neighbouring cells are strongly correlated;
    averaging does not reduce sigma by sqrt(n)).  A coarse cell is filled only if >= 50 % of its area has values;
    its class is the worst class of the contributing fine cells (conservative)."""
    cp = cp or CovParams()
    cg = make_grid(name)
    cg.lat = cg.lat[(cg.lat >= grid.lat.min()) & (cg.lat <= grid.lat.max())]
    cg.lon = cg.lon[(cg.lon >= grid.lon.min()) & (cg.lon <= grid.lon.max())]
    wa = _agg_weights(grid.lat, (cg.lat, cg.step, grid.step))
    wo = _agg_weights(grid.lon, (cg.lon, cg.step, grid.step))
    ny, nx = len(cg.lat), len(cg.lon)
    W = lambda a: np.einsum("ia,jb,ab->ij", wa, wo, a)  # noqa: E731
    area = W(np.ones_like(grid.h))
    cg.h = W(grid.h) / area
    cg.h_std = np.sqrt(np.maximum(W(grid.h_std ** 2 + grid.h ** 2) / area - cg.h ** 2, 0))
    cg.land = W(grid.land.astype(float)) / area >= 0.5
    cg.dem_source = grid.dem_source + f" (aggregated from {grid.name})"
    # children of each coarse row/column (<= K per axis, padded with weight 0)
    def kids(w):
        K = int((w > 0).sum(1).max())
        idx = np.zeros((w.shape[0], K), dtype=int)
        ww = np.zeros((w.shape[0], K))
        for i_ in range(w.shape[0]):
            nz = np.nonzero(w[i_])[0]
            idx[i_, :len(nz)] = nz
            idx[i_, len(nz):] = nz[-1] if len(nz) else 0
            ww[i_, :len(nz)] = w[i_, nz]
        return idx, ww
    ia, wia = kids(wa)
    jo, wjo = kids(wo)
    Ka, Ko = ia.shape[1], jo.shape[1]
    # correlation of fine-cell errors within a coarse cell: rho = exp(-d/L) between child cell centres
    LAk = grid.lat[ia][:, None, :, None] + 0 * grid.lon[jo][None, :, None, :]
    LOk = grid.lon[jo][None, :, None, :] + 0 * grid.lat[ia][:, None, :, None]
    U = to_unit(LAk.ravel(), LOk.ravel()).reshape(LAk.shape + (3,))
    rho = np.exp(-R_EARTH_KM * np.arccos(np.clip(np.einsum("ijabx,ijcdx->ijabcd", U, U), -1, 1)) / cp.L_km)
    Wt = wia[:, None, :, None] * wjo[None, :, None, :]                   # (ny, nx, Ka, Ko)
    out = []
    for r in results:
        layers = {}
        for nm in ("pwv", "bg", "resid", "ztd", "zwd", "zhd", "info", "d1", "ps", "tm", "res_eff"):
            a = getattr(r, nm)
            ok = np.isfinite(a)
            num = W(np.where(ok, a, 0.0))
            den = W(ok.astype(float))
            layers[nm] = np.where(den >= 0.5 * area, num / np.maximum(den, 1e-12), np.nan)
        sg = {}
        for key, src in (("sig", r.sigma), ("sig_ztd", r.sig_ztd)):
            sv = src[ia[:, None, :, None], jo[None, :, None, :]]           # (ny, nx, Ka, Ko)
            ww = Wt * np.isfinite(sv)
            ws = ww * np.nan_to_num(sv)
            q = np.einsum("ijab,ijcd,ijabcd->ij", ws, ws, rho)
            sg[key] = np.sqrt(q) / np.maximum(ww.sum((2, 3)), 1e-12)
        cc = r.cls[ia[:, None, :, None], jo[None, :, None, :]].astype(int)
        cc = np.where((Wt > 0) & (cc > 0), cc, 0)
        cls = cc.max((2, 3)).astype(np.int8)
        sig, sig_ztd = sg["sig"], sg["sig_ztd"]
        sig[~np.isfinite(layers["pwv"])] = np.nan
        sig_ztd[~np.isfinite(layers["ztd"])] = np.nan
        n50 = np.rint(W(r.n50.astype(float)) / area).astype(int)
        n100 = np.rint(W(r.n100.astype(float)) / area).astype(int)
        cls[~cg.land] = 0
        for a in (layers["pwv"], layers["ztd"], layers["zwd"], sig, sig_ztd):
            a[(cls == 4) | (cls == 0)] = np.nan          # class-4/0 cells carry fill values (TDS § 26, § 36)
        out.append(SlotResult(layers["pwv"], sig, layers["bg"], layers["resid"], layers["ztd"], layers["zwd"],
                              layers["zhd"], sig_ztd, layers["info"], cls, layers["d1"], n50, n100,
                              layers["res_eff"], layers["ps"], layers["tm"], r.avail, r.n_used, r.outliers))
    return cg, out


# ================================================================================ output
SUPPORT_CODES = {0.1: 1, 0.125: 2, 0.25: 3, 0.5: 4, None: 0}     # supportable_grid_spacing flag values


def support_layer(grid, support):
    """Static supportable_grid_spacing layer from a {(lat0, lon0) 1-deg tile: spacing_deg|None} table; -1 where
    not determined (before the § 31 experiment)."""
    lay = np.full((len(grid.lat), len(grid.lon)), -1, dtype=np.int8)
    if not support:
        return lay
    for (la0, lo0), sp in support.items():
        mi = (grid.lat >= la0) & (grid.lat < la0 + 1.0)
        mj = (grid.lon >= lo0) & (grid.lon < lo0 + 1.0)
        lay[np.ix_(mi, mj)] = SUPPORT_CODES.get(sp, 0)
    return lay


def write_grid_netcdf(path, grid, slots_t, results, stations, attrs, support=None):
    """NetCDF-4 CF-1.10 + ACDD-1.3 grid file (TDS § 27.3)."""
    import netCDF4
    ny, nx = len(grid.lat), len(grid.lon)
    with netCDF4.Dataset(path, "w", format="NETCDF4") as nc:
        nc.createDimension("time", len(slots_t))
        nc.createDimension("lat", ny)
        nc.createDimension("lon", nx)
        nc.createDimension("bnds", 2)
        nc.createDimension("station", len(stations))
        t = nc.createVariable("time", "i4", ("time",))
        t.units = "minutes since 2000-01-01 00:00:00"
        t.calendar = "standard"
        t.standard_name = "time"
        t0 = ts.to_ns(2000, 1, 1)
        t[:] = [(s - t0) // (60 * ts.NS) for s in slots_t]
        la = nc.createVariable("lat", "f8", ("lat",))
        la.units, la.standard_name, la.bounds = "degrees_north", "latitude", "lat_bnds"
        la[:] = grid.lat
        lo = nc.createVariable("lon", "f8", ("lon",))
        lo.units, lo.standard_name, lo.bounds = "degrees_east", "longitude", "lon_bnds"
        lo[:] = grid.lon
        lb = nc.createVariable("lat_bnds", "f8", ("lat", "bnds"))
        lb[:] = np.stack([grid.lat - grid.step / 2, grid.lat + grid.step / 2], -1)
        ob = nc.createVariable("lon_bnds", "f8", ("lon", "bnds"))
        ob[:] = np.stack([grid.lon - grid.step / 2, grid.lon + grid.step / 2], -1)
        crs = nc.createVariable("crs", "i4")
        crs.grid_mapping_name = "latitude_longitude"
        crs.epsg_code = "EPSG:4326"
        kw = dict(zlib=True, complevel=4, shuffle=True, chunksizes=(min(4, len(slots_t)), ny, nx))
        layers = [("PWV", "pwv", "mm", "lwe_thickness_of_atmosphere_mass_content_of_water_vapor"),
                  ("ZWD", "zwd", "mm", None), ("ZTD", "ztd", "mm", None), ("ZHD", "zhd", "mm", None),
                  ("sigma_PWV", "sigma", "mm", None), ("sigma_ZTD", "sig_ztd", "mm", None),
                  ("background_PWV", "bg", "mm", None), ("gnss_residual", "resid", "mm", None),
                  ("distance_to_nearest_station", "d1", "km", None),
                  ("gnss_information_fraction", "info", "1", None), ("surface_pressure", "ps", "hPa", None),
                  ("Tm", "tm", "K", None), ("effective_resolution_km", "res_eff", "km", None)]
        for name, attr, units, sn in layers:
            v = nc.createVariable(name, "f4", ("time", "lat", "lon"), fill_value=np.float32(-9999.0), **kw)
            v.units = units
            v.grid_mapping = "crs"
            if sn:
                v.standard_name = sn
            arr = np.array([getattr(r, attr) for r in results], dtype=np.float32)
            v[:] = np.where(np.isfinite(arr), arr, -9999.0)
        c = nc.createVariable("confidence_class", "i1", ("time", "lat", "lon"), **kw)
        c.flag_values = np.array([0, 1, 2, 3, 4], dtype=np.int8)
        c.flag_meanings = "outside_domain_or_sea well_covered moderate sparse_background_dominated unsupported"
        c[:] = np.array([r.cls for r in results], dtype=np.int8)
        for name, attr in (("n_stations_50km", "n50"), ("n_stations_100km", "n100")):
            v = nc.createVariable(name, "i2", ("time", "lat", "lon"), **kw)
            v[:] = np.array([getattr(r, attr) for r in results], dtype=np.int16)
        for name, arr, un in (("dem_height", grid.h, "m"), ("dem_height_std", grid.h_std, "m")):
            v = nc.createVariable(name, "f4", ("lat", "lon"))
            v.units = un
            v[:] = arr
        sg = nc.createVariable("supportable_grid_spacing", "i1", ("lat", "lon"))
        sg.flag_values = np.array([-1, 0, 1, 2, 3, 4], dtype=np.int8)
        sg.flag_meanings = "not_determined background_only 0.1deg 0.125deg 0.25deg 0.5deg"
        sg.comment = "TDS § 24.4; from the § 31 spatial validation experiment (version in global attributes)"
        sg[:] = support_layer(grid, support)
        lm = nc.createVariable("land_sea_mask", "i1", ("lat", "lon"))
        lm[:] = grid.land.astype(np.int8)
        sn = nc.createVariable("station_name", str, ("station",))
        for i, s in enumerate(stations):
            sn[i] = s.station
        sa = nc.createVariable("station_availability", "i1", ("time", "station"))
        sa[:] = np.array([r.avail for r in results], dtype=np.int8)
        nc.Conventions = "CF-1.10, ACDD-1.3"
        nc.comment = ("Grid spacing is not atmospheric resolution: see effective_resolution_km and "
                      "confidence_class (TDS § 0.3). Class-4 cells carry fill values in GNSS-corrected layers.")
        for k, v in attrs.items():
            setattr(nc, k, v if isinstance(v, (str, int, float)) else json.dumps(v, default=str))
