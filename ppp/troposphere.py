"""Troposphere models: VMF3 grids and mapping function, GPT3, GMF (benchmark), Saastamoinen, gradients.

TDS § 6.1-6.3.  References:
* VMF3 / GPT3: Landskron & Boehm (2018) J. Geod. 92:349-360; TU Wien reference code vmf3_ht.m, gpt3_1_fast.m.
* Height reduction of gridded ZHD/ZWD: Kouba (2008) J. Geod. 82:193-205.
* Saastamoinen (1972) / Davis et al. (1985): ZHD = 0.0022768 P / (1 - 0.00266 cos 2phi - 0.28e-6 H).
* Askne & Nordius (1987) wet delay from e, Tm, lambda.
* Chen & Herring (1997) gradient mapping function, C = 0.0032 (IERS 2010 eq. 9.12).
* GMF: Boehm et al. (2006), IERS 2010 GMF.F (benchmark/emergency only, TDS § 6.1).
"""
import os
from dataclasses import dataclass, field

import numpy as np

from . import log as plog
from . import settings
from . import timesys as ts

LOG = plog.get()


# ======================================================================== basic formulas
def saastamoinen_zhd(p_hpa, lat, h_m):
    """Zenith hydrostatic delay (m); p in hPa, lat rad, h (orthometric) m [TDS eq. 6.2 / 18.1]."""
    return 0.0022768 * np.asarray(p_hpa) / (1.0 - 0.00266 * np.cos(2 * lat) - 0.28e-6 * np.asarray(h_m))


def saastamoinen_pressure(zhd, lat, h_m):
    return np.asarray(zhd) * (1.0 - 0.00266 * np.cos(2 * lat) - 0.28e-6 * np.asarray(h_m)) / 0.0022768


def askne_nordius_zwd(e_hpa, tm_k, la):
    """Askne & Nordius (1987) zenith wet delay (m) as in TU Wien asknewet.m."""
    k1, k2, k3 = 77.604, 64.79, 377600.0
    k2p = k2 - k1 * 18.0152 / 28.9644
    rd = 8314.0 / 28.9644
    gm = 9.80665
    return 1e-6 * (k2p + k3 / np.asarray(tm_k)) * rd / (np.asarray(la) + 1.0) / gm * np.asarray(e_hpa)


def chen_herring_mg(el):
    """Gradient mapping function m_g(e) = 1 / (sin e tan e + C) (TDS § 6.1)."""
    return 1.0 / (np.sin(el) * np.tan(el) + settings.GRADIENT_MF_C)


def _cf(a, b, c, s):
    return (1 + a / (1 + b / (1 + c))) / (s + a / (s + b / (s + c)))


def niell_height_corr(el, h_m):
    """Niell (1996) hydrostatic height correction used by VMF3/GMF."""
    s = np.sin(el)
    return (1.0 / s - _cf(2.53e-5, 5.49e-3, 1.14e-3, s)) * np.asarray(h_m) / 1000.0


def std_pressure(h_m):
    """Standard-atmosphere pressure (hPa) (Berg 1948) - emergency fallback only."""
    return 1013.25 * (1 - 2.2557e-5 * np.asarray(h_m)) ** 5.2568


# ======================================================================== spherical harmonics helper
def _legendre_vw(lat, lon, nmax):
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    V = np.zeros((nmax + 2, nmax + 2))
    W = np.zeros((nmax + 2, nmax + 2))
    V[0, 0] = 1.0
    V[1, 0] = z
    for n in range(2, nmax + 1):
        V[n, 0] = ((2 * n - 1) * z * V[n - 1, 0] - (n - 1) * V[n - 2, 0]) / n
    for m in range(1, nmax + 1):
        V[m, m] = (2 * m - 1) * (x * V[m - 1, m - 1] - y * W[m - 1, m - 1])
        W[m, m] = (2 * m - 1) * (x * W[m - 1, m - 1] + y * V[m - 1, m - 1])
        if m < nmax:
            V[m + 1, m] = (2 * m + 1) * z * V[m, m]
            W[m + 1, m] = (2 * m + 1) * z * W[m, m]
        for n in range(m + 2, nmax + 1):
            V[n, m] = ((2 * n - 1) * z * V[n - 1, m] - (n + m - 1) * V[n - 2, m]) / (n - m)
            W[n, m] = ((2 * n - 1) * z * W[n - 1, m] - (n + m - 1) * W[n - 2, m]) / (n - m)
    vv, ww = [], []
    for n in range(nmax + 1):
        for m in range(n + 1):
            vv.append(V[n, m])
            ww.append(W[n, m])
    return np.array(vv), np.array(ww)


def _load_blocks(fname):
    blocks, cur = {}, None
    with open(os.path.join(settings.DATA_DIR, fname)) as fh:
        for L in fh:
            if L.startswith("#") or not L.strip():
                continue
            if L.startswith("@"):
                cur = L.split()[1]
                blocks[cur] = []
                continue
            blocks[cur].extend(float(v) for v in L.split())
    return {k: np.array(v) for k, v in blocks.items()}


_VMF3 = None
_GMF = None


def _doy_frac(mjd):
    """Day of year + fraction (1-based) as in the TU Wien codes."""
    mjd = float(mjd)
    y, m, d = ts.date_from_mjd(int(np.floor(mjd)))
    doy = ts.mjd_from_date(y, m, d) - ts.mjd_from_date(y, 1, 1) + 1
    return doy + mjd - np.floor(mjd)


def vmf3_bc(mjd, lat, lon):
    """Empirical VMF3 b_h, b_w, c_h, c_w (Landskron & Boehm 2018, vmf3_ht.m)."""
    global _VMF3
    if _VMF3 is None:
        raw = _load_blocks("vmf3_bc_coefficients.txt")
        _VMF3 = {k: v.reshape(91, 5) for k, v in raw.items()}
    V, W = _legendre_vw(lat, lon, 12)
    doy = _doy_frac(mjd)
    tv = np.array([1.0, np.cos(doy / 365.25 * 2 * np.pi), np.sin(doy / 365.25 * 2 * np.pi),
                   np.cos(doy / 365.25 * 4 * np.pi), np.sin(doy / 365.25 * 4 * np.pi)])
    out = {}
    for q in ("bh", "bw", "ch", "cw"):
        coef = V @ _VMF3["anm_" + q] + W @ _VMF3["bnm_" + q]      # (5,)
        out[q] = float(coef @ tv)
    return out["bh"], out["bw"], out["ch"], out["cw"]


def vmf3_ht(ah, aw, mjd, lat, lon, h_ell, el, bc=None):
    """VMF3 hydrostatic and wet mapping functions (with Niell height correction for the hydrostatic part).

    ah, aw from the VMF3 grid (or GPT3); el elevation (rad), array allowed.
    """
    bh, bw, ch, cw = bc if bc is not None else vmf3_bc(mjd, lat, lon)
    s = np.sin(el)
    mfh = _cf(ah, bh, ch, s) + niell_height_corr(el, h_ell)
    mfw = _cf(aw, bw, cw, s)
    return mfh, mfw


def gmf(mjd, lat, lon, h_ell, el):
    """Global Mapping Function (IERS 2010 GMF.F). Benchmark / emergency fallback only (TDS § 6.1)."""
    global _GMF
    if _GMF is None:
        _GMF = _load_blocks("gmf_coefficients.txt")
    V, W = _legendre_vw(lat, lon, 9)
    doy = mjd - 44239.0 + 1 - 28
    bh, c0h = 0.0029, 0.062
    if lat < 0:
        phh, c11h, c10h = np.pi, 0.007, 0.002
    else:
        phh, c11h, c10h = 0.0, 0.005, 0.001
    ch = c0h + ((np.cos(doy / 365.25 * 2 * np.pi + phh) + 1) * c11h / 2 + c10h) * (1 - np.cos(lat))
    ahm = V @ _GMF["ah_mean"] + W @ _GMF["bh_mean"]
    aha = V @ _GMF["ah_amp"] + W @ _GMF["bh_amp"]
    ah = (ahm + aha * np.cos(doy / 365.25 * 2 * np.pi)) * 1e-5
    awm = V @ _GMF["aw_mean"] + W @ _GMF["bw_mean"]
    awa = V @ _GMF["aw_amp"] + W @ _GMF["bw_amp"]
    aw = (awm + awa * np.cos(doy / 365.25 * 2 * np.pi)) * 1e-5
    s = np.sin(el)
    mfh = _cf(ah, bh, ch, s) + niell_height_corr(el, h_ell)
    mfw = _cf(aw, 0.00146, 0.04391, s)
    return mfh, mfw


# ======================================================================== GPT3 grid
@dataclass
class GPT3Grid:
    cols: dict                    # name -> (n,5) or (n,) arrays
    path: str = ""
    has_mf: bool = False


def read_gpt3(path):
    """Read a GPT3 1-deg grid; columns located from the header tokens (robust to reduced variants)."""
    with open(path) as fh:
        header = fh.readline()
    data = np.loadtxt(path, comments="%")
    toks = header.lstrip("%").split()
    cols, k = {}, 0
    i = 0
    while i < len(toks):
        t = toks[i]
        if ":" in t:
            name = t.split(":")[0].lower().replace("_", "")
            cols[name] = data[:, k:k + 5].copy()
            k += 5
            i += 5
        else:
            cols[t.lower()] = data[:, k].copy()
            k += 1
            i += 1
    if "q" in cols:
        cols["q"] = cols["q"] / 1000.0
    if "dt" in cols:
        cols["dt"] = cols["dt"] / 1000.0
    for nm in ("ah", "aw"):
        if nm in cols:
            cols[nm] = cols[nm] / 1000.0
    g = GPT3Grid(cols, path, has_mf=("ah" in cols and "aw" in cols))
    return g


def gpt3(grid, mjd, lat, lon, h_ell):
    """GPT3 (1 deg) meteorology and mapping-function coefficients at a site (gpt3_1_fast.m logic).

    Returns dict p (hPa), T (degC), dT (deg/km), Tm (K), e (hPa), ah, aw (None if absent), la, undu (m).
    """
    doy = _doy_frac(mjd)
    cs = np.array([1.0, np.cos(doy / 365.25 * 2 * np.pi), np.sin(doy / 365.25 * 2 * np.pi),
                   np.cos(doy / 365.25 * 4 * np.pi), np.sin(doy / 365.25 * 4 * np.pi)])
    plon = np.degrees(lon) % 360.0
    ppod = 90.0 - np.degrees(lat)
    ipod = int(np.floor(ppod + 1))
    ilon = int(np.floor(plon + 1))
    diffpod = ppod - (ipod - 0.5)
    difflon = plon - (ilon - 0.5)
    ipod = min(ipod, 180)
    ilon = 1 if ilon == 361 else (360 if ilon == 0 else ilon)
    ipod1 = ipod + int(np.sign(diffpod) or 1)
    ilon1 = ilon + int(np.sign(difflon) or 1)
    ilon1 = 1 if ilon1 == 361 else (360 if ilon1 == 0 else ilon1)
    bilinear = 0.5 < ppod < 179.5
    idx = [(ipod - 1) * 360 + ilon - 1]
    if bilinear:
        idx += [(ipod1 - 1) * 360 + ilon - 1, (ipod - 1) * 360 + ilon1 - 1, (ipod1 - 1) * 360 + ilon1 - 1]
    c = grid.cols
    gm, dmtr, rg = 9.80665, 28.965e-3, 8.3143
    res = {k: [] for k in ("p", "T", "dT", "Tm", "e", "ah", "aw", "la", "undu")}
    for ix in idx:
        undu = c["undu"][ix]
        hgt = h_ell - undu
        T0 = c["t"][ix] @ cs
        p0 = c["p"][ix] @ cs
        Q = c["q"][ix] @ cs
        dT = c["dt"][ix] @ cs
        redh = hgt - c["hs"][ix]
        T = T0 + dT * redh - 273.15
        Tv = T0 * (1 + 0.6077 * Q)
        p = p0 * np.exp(-gm * dmtr / (rg * Tv) * redh) / 100.0
        la = c["lambda"][ix] @ cs
        Tm = c["tm"][ix] @ cs
        e0 = Q * p0 / (0.622 + 0.378 * Q) / 100.0
        e = e0 * (100.0 * p / p0) ** (la + 1)
        res["p"].append(p)
        res["T"].append(T)
        res["dT"].append(dT * 1000.0)
        res["Tm"].append(Tm)
        res["e"].append(e)
        res["la"].append(la)
        res["undu"].append(undu)
        res["ah"].append(c["ah"][ix] @ cs if grid.has_mf else np.nan)
        res["aw"].append(c["aw"][ix] @ cs if grid.has_mf else np.nan)
    if bilinear:
        dp1, dl1 = abs(diffpod), abs(difflon)
        dp2, dl2 = 1 - dp1, 1 - dl1
        out = {}
        for k, v in res.items():
            r1 = dp2 * v[0] + dp1 * v[1]
            r2 = dp2 * v[2] + dp1 * v[3]
            out[k] = dl2 * r1 + dl1 * r2
    else:
        out = {k: v[0] for k, v in res.items()}
    if not grid.has_mf:
        out["ah"] = out["aw"] = None
    return out


# ======================================================================== VMF3 grids
@dataclass
class VMF3Epoch:
    t_ns: int
    ah: np.ndarray                # (180, 360) lat 89.5..-89.5, lon 0.5..359.5
    aw: np.ndarray
    zhd: np.ndarray
    zwd: np.ndarray
    path: str = ""


def read_vmf3_grid(path, t_ns):
    arr = np.loadtxt(path, comments="!")
    lat, lon = arr[:, 0], arr[:, 1]
    i = np.round(89.5 - lat).astype(int)
    j = np.round(lon - 0.5).astype(int) % 360
    out = []
    for c in range(2, 6):
        g = np.full((180, 360), np.nan)
        g[i, j] = arr[:, c]
        out.append(g)
    return VMF3Epoch(t_ns, *out, path=path)


def read_orography(path):
    vals = np.array(open(path).read().split(), dtype=float)
    if vals.size == 180 * 360:
        return vals.reshape(180, 360)
    if vals.size == 181 * 361:
        # node-registered 90..-90, 0..360 -> average to cell centres
        g = vals.reshape(181, 361)
        return 0.25 * (g[:-1, :-1] + g[1:, :-1] + g[:-1, 1:] + g[1:, 1:])
    raise ValueError(f"unexpected orography size {vals.size}")


def _grid_corners(lat_deg, lon_deg):
    lon_deg = lon_deg % 360.0
    fi = 89.5 - lat_deg
    i0 = int(np.clip(np.floor(fi), 0, 178))
    fj = lon_deg - 0.5
    j0 = int(np.floor(fj)) % 360
    wi = fi - i0
    wj = (fj - np.floor(fj))
    return i0, j0, wi, wj


def vmf3_station(epochs, oro, lat, lon, h_ell):
    """Station time series (per 6-h epoch) of ah, aw, zhd, zwd from VMF3 grids with Kouba (2008) height
    reduction of zhd/zwd at each grid point before bilinear interpolation."""
    latd, lond = np.degrees(lat), np.degrees(lon)
    i0, j0, wi, wj = _grid_corners(latd, lond)
    pts = [(i0, j0, (1 - wi) * (1 - wj)), (i0 + 1, j0, wi * (1 - wj)), (i0, (j0 + 1) % 360, (1 - wi) * wj),
           (i0 + 1, (j0 + 1) % 360, wi * wj)]
    rows = []
    for ep in epochs:
        ah = aw = zhd = zwd = 0.0
        for (ii, jj, w) in pts:
            hg = oro[ii, jj]
            latg = np.radians(89.5 - ii)
            # hydrostatic: grid ZHD -> pressure at grid height -> station height -> ZHD (Kouba 2008)
            pg = saastamoinen_pressure(ep.zhd[ii, jj], latg, hg)
            ps = pg * (1.0 - 0.0000226 * (h_ell - hg)) ** 5.225
            zhd += w * saastamoinen_zhd(ps, lat, h_ell)
            zwd += w * ep.zwd[ii, jj] * np.exp(-(h_ell - hg) / settings.ZWD_SCALE_HEIGHT_M)
            ah += w * ep.ah[ii, jj]
            aw += w * ep.aw[ii, jj]
        rows.append((ep.t_ns, ah, aw, zhd, zwd))
    return np.array(rows, dtype=float)


# ======================================================================== station troposphere
@dataclass
class StationTropo:
    """A priori troposphere at a station: ZHD0/ZWD0 and mapping-function coefficients vs time."""
    t_ns: np.ndarray
    zhd: np.ndarray
    zwd: np.ndarray
    ah: np.ndarray
    aw: np.ndarray
    lat: float
    lon: float
    h_ell: float
    zhd_source: str               # VMF3_GRID | GPT3 | STANDARD_ATMOSPHERE | BAROMETER
    mapping_function: str         # VMF3 | GPT3 (VMF3 with GPT3 a) | GMF
    flags: set = field(default_factory=set)
    files: list = field(default_factory=list)
    met: dict = field(default_factory=dict)   # GPT3 p, T, e, Tm, la, undu at window midpoint

    def at(self, t_ns):
        tq = np.asarray(t_ns, dtype=float)
        tt = self.t_ns.astype(float)
        f = lambda v: np.interp(tq, tt, v)  # noqa: E731
        return f(self.zhd), f(self.zwd), f(self.ah), f(self.aw)

    def mapping(self, t_ns, el):
        """m_h, m_w, m_g for elevations el (rad) at times t_ns (same shape)."""
        mjd = float(ts.mjd(int(np.median(t_ns))))
        mg = chen_herring_mg(el)
        if self.mapping_function == "GMF":
            mh, mw = gmf(mjd, self.lat, self.lon, self.h_ell, el)
            return mh, mw, mg
        _, _, ah, aw = self.at(t_ns)
        bc = vmf3_bc(mjd, self.lat, self.lon)
        mh, mw = vmf3_ht(ah, aw, mjd, self.lat, self.lon, self.h_ell, el, bc)
        return mh, mw, mg


def build_station_tropo(lat, lon, h_ell, win_start, win_end, vmf3_files=None, oro_path=None, gpt3_grid=None,
                        barometer=None):
    """A priori ZHD/ZWD + mapping function following the TDS priority (§ 6.2-6.3, § 9.11)."""
    flags = set()
    files = []
    mid = (win_start + win_end) // 2
    mjd_mid = float(ts.mjd(mid))
    met = {}
    g = None
    if gpt3_grid is not None:
        g = gpt3(gpt3_grid, mjd_mid, lat, lon, h_ell)
        met = {k: (float(v) if v is not None else None) for k, v in g.items()}
    # 1) VMF3 grids
    if vmf3_files and oro_path:
        try:
            oro = read_orography(oro_path)
            eps = [read_vmf3_grid(p, t) for (t, p, _m) in vmf3_files]
            rows = vmf3_station(eps, oro, lat, lon, h_ell)
            files = [os.path.basename(p) for (_t, p, _m) in vmf3_files] + [os.path.basename(oro_path)]
            st = StationTropo(rows[:, 0].astype(np.int64), rows[:, 3], rows[:, 4], rows[:, 1], rows[:, 2], lat, lon,
                              h_ell, "VMF3_GRID", "VMF3", flags, files, met)
            if barometer is not None:
                st = _apply_barometer(st, barometer)
            return st
        except (OSError, ValueError) as exc:
            LOG.warning("VMF3_UNREADABLE: %s; falling back to GPT3", exc)
    t = np.array([win_start - 6 * 3600 * ts.NS, win_end + 6 * 3600 * ts.NS], dtype=np.int64)
    if g is not None:
        # 2) GPT3: pressure -> Saastamoinen (ZHD_CLIMATOLOGY), Askne-Nordius ZWD, GPT3 ah/aw for VMF3
        h_orth = h_ell - g["undu"]
        vals = []
        for tq in t:
            gg = gpt3(gpt3_grid, float(ts.mjd(int(tq))), lat, lon, h_ell)
            vals.append((saastamoinen_zhd(gg["p"], lat, h_orth), askne_nordius_zwd(gg["e"], gg["Tm"], gg["la"]),
                         gg["ah"] if gg["ah"] is not None else np.nan, gg["aw"] if gg["aw"] is not None else np.nan))
        vals = np.array(vals, dtype=float)
        flags.add("ZHD_CLIMATOLOGY")
        mf = "GPT3" if gpt3_grid.has_mf else "GMF"
        if not gpt3_grid.has_mf:
            flags.add("MF_FALLBACK")
            LOG.warning("MF_FALLBACK: GPT3 grid without a_h/a_w coefficients and no VMF3 grid -> GMF mapping "
                        "function used (benchmark model; see docs/decisions.md D-005)")
        LOG.warning("ZHD_CLIMATOLOGY: VMF3 grids not available; a priori ZHD from GPT3 climatology "
                    "(acceptable for ZTD, not for PWV)")
        st = StationTropo(t, vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3], lat, lon, h_ell, "GPT3", mf, flags,
                          [os.path.basename(gpt3_grid.path)], met)
    else:
        # 3) emergency: standard atmosphere + GMF (never silent)
        flags |= {"ZHD_CLIMATOLOGY", "MF_FALLBACK"}
        LOG.warning("MF_FALLBACK: neither VMF3 nor GPT3 available; standard-atmosphere ZHD and GMF mapping "
                    "(emergency fallback, docs/decisions.md D-005)")
        p = std_pressure(h_ell)
        z = saastamoinen_zhd(p, lat, h_ell)
        st = StationTropo(t, np.array([z, z]), np.array([0.1, 0.1]), np.array([np.nan] * 2), np.array([np.nan] * 2),
                          lat, lon, h_ell, "STANDARD_ATMOSPHERE", "GMF", flags, [], met)
    if barometer is not None:
        st = _apply_barometer(st, barometer)
    return st


def _apply_barometer(st, baro):
    """Replace ZHD0 by Saastamoinen with barometer pressure reduced to the ARP (priority 1, TDS § 6.2)."""
    t_ns, p_ant, h_orth = baro
    z = saastamoinen_zhd(p_ant, st.lat, h_orth)
    st.zhd = np.interp(st.t_ns.astype(float), np.asarray(t_ns, dtype=float), z)
    st.zhd_source = "BAROMETER"
    st.flags.discard("ZHD_CLIMATOLOGY")
    return st
