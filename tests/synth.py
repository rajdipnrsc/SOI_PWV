"""Synthetic GNSS observation generator for Level-0 estimator / preprocessing tests (known truth)."""
import numpy as np

from ppp import coords as co
from ppp import estimator as est
from ppp import preprocess as pp
from ppp import settings
from ppp import troposphere as tr

C = settings.C_LIGHT
F1, F2 = settings.FREQ[("G", "1")], settings.FREQ[("G", "2")]
A1 = F1 ** 2 / (F1 ** 2 - F2 ** 2)
A2 = F2 ** 2 / (F1 ** 2 - F2 ** 2)
FIF = np.sqrt(A1 ** 2 + A2 ** 2)


def geometry(ne=2880, dt=30.0, ns=12, seed=1):
    """Elevation/azimuth tracks: each satellite has 6-h passes repeating every ~12 h."""
    rng = np.random.default_rng(seed)
    t = np.arange(ne) * dt
    el = np.full((ne, ns), -0.5)
    az = np.zeros((ne, ns))
    T = 6 * 3600.0
    for j in range(ns):
        emax = np.radians(rng.uniform(25, 88))
        a0 = rng.uniform(0, 2 * np.pi)
        for start in (j * 3600.0 - 6 * 3600, j * 3600.0 + 6 * 3600 - 6 * 3600 + 12 * 3600 * 0.99, j * 3600.0 +
                      6 * 3600):
            ph = (t - start) / T
            m = (ph > 0) & (ph < 1)
            el[m, j] = np.arcsin(np.sin(emax) * np.sin(np.pi * ph[m]))
            az[m, j] = (a0 + np.pi * ph[m]) % (2 * np.pi)
    return t, el, az


def make_obs(mode="fixed", seed=3, dX=(0.0, 0.0, 0.0), sig_zwd=0.006, sig_g=0.0005, ne=2880):
    rng = np.random.default_rng(seed)
    t, el, az = geometry(ne=ne)
    ne_, ns = el.shape
    X0 = co.geodetic_to_ecef(np.radians(17.4), np.radians(78.55), 420.0)
    lat, lon, _ = co.ecef_to_geodetic(X0)
    R = co.rot_xyz2enu(lat, lon)
    enu = np.stack([np.cos(el) * np.sin(az), np.cos(el) * np.cos(az), np.sin(el)], axis=-1)
    los = enu @ R                                          # ENU -> ECEF unit vectors
    vis = el > np.radians(7.0)
    mh = 1.0 / np.sin(np.maximum(el, 0.05))
    mw = mh * 1.001
    mg = tr.chen_herring_mg(np.maximum(el, 0.05))
    dt = np.diff(t, prepend=t[0])
    zwd = 0.15 + np.cumsum(rng.normal(0, sig_zwd * np.sqrt(dt / 3600.0)))
    gn = np.cumsum(rng.normal(0, sig_g * np.sqrt(dt / 3600.0)))
    ge = np.cumsum(rng.normal(0, sig_g * np.sqrt(dt / 3600.0)))
    clk = rng.normal(0, 3e4, ne_)
    zhd0 = np.full(ne_, 2.2)
    rho = 2.2e7 + 1e6 * np.sin(el)
    tropo = mh * zhd0[:, None] + mw * zwd[:, None] + mg * (gn[:, None] * np.cos(az) + ge[:, None] * np.sin(az))
    geo = rho - los @ np.asarray(dX)
    sp = np.sqrt((0.003 ** 2 + 0.003 ** 2 / np.sin(np.maximum(el, 0.05)) ** 2)) * FIF
    sc = np.sqrt((0.3 ** 2 + 0.3 ** 2 / np.sin(np.maximum(el, 0.05)) ** 2)) * FIF
    # arcs: one per continuous visibility period
    arc = np.full((ne_, ns), -1, dtype=int)
    meta = {}
    nid = 0
    B = np.zeros((ne_, ns))
    for j in range(ns):
        prev = False
        for k in range(ne_):
            if vis[k, j] and not prev:
                nid += 1
                meta[nid] = {"sat": j, "first": k, "pair": 0}
                amb = rng.uniform(-50, 50)
            if vis[k, j]:
                arc[k, j] = nid
                B[k, j] = amb
                meta[nid]["last"] = k
            prev = vis[k, j]
    Pif = geo + tropo + clk[:, None] + rng.normal(0, 1, (ne_, ns)) * sc
    Lif = geo + tropo + clk[:, None] + B + rng.normal(0, 1, (ne_, ns)) * sp
    Pif[~vis] = np.nan
    Lif[~vis] = np.nan
    obs = est.ObsSet(epochs=(t * 1e9).astype(np.int64), t_s=t, sats=[f"G{j + 1:02d}" for j in range(ns)], Pif=Pif,
                     Lif=Lif, rho=rho, dts=np.zeros((ne_, ns)), wind=np.zeros((ne_, ns)), el=el, az=az, los=los,
                     mh=mh, mw=mw, mg=mg, zhd0=zhd0, zwd0=np.full(ne_, 0.12), usable=vis.copy(), f_if=FIF, X0=X0)
    arcs = pp.ArcInfo(arc, {}, len(meta), np.zeros(ne_, dtype=bool), np.zeros(ne_), meta)
    truth = {"zwd": zwd, "gn": gn, "ge": ge, "clk": clk, "dX": np.asarray(dX)}
    return obs, arcs, truth


def make_raw(ne=2880, seed=21, sig_p=0.15, sig_l=0.0015, slips=()):
    """Raw dual-frequency observables (SelectedObs) with smooth ionosphere and optional injected slips
    slips: list of (epoch, sat, dN1, dN2)."""
    rng = np.random.default_rng(seed)
    t, el, az = geometry(ne=ne)
    ne_, ns = el.shape
    vis = el > np.radians(7.0)
    lam1, lam2 = C / F1, C / F2
    g = (F1 / F2) ** 2
    rho = 2.2e7 + 2e6 * (1 - np.sin(np.maximum(el, 0)))
    tec = 20 + 10 * np.sin(2 * np.pi * t / 86400.0)[:, None]          # TECU, smooth
    I1 = 0.162 * tec / np.sin(np.maximum(el, 0.1))                    # m on L1 (40.3e16/f1^2 per TECU)
    N1 = rng.integers(-1000, 1000, ns).astype(float)[None, :] * np.ones((ne_, 1))
    N2 = rng.integers(-1000, 1000, ns).astype(float)[None, :] * np.ones((ne_, 1))
    for (k, j, d1, d2) in slips:
        N1[k:, j] += d1
        N2[k:, j] += d2
    wfac = 1 + 1 / np.sin(np.maximum(el, 0.1))                        # sigma = s (1 + 1/sin e)
    P1 = rho + I1 + rng.normal(0, 1, (ne_, ns)) * sig_p * wfac
    P2 = rho + g * I1 + rng.normal(0, 1, (ne_, ns)) * sig_p * wfac
    L1 = rho - I1 + lam1 * N1 + rng.normal(0, 1, (ne_, ns)) * sig_l * wfac
    L2 = rho - g * I1 + lam2 * N2 + rng.normal(0, 1, (ne_, ns)) * sig_l * wfac
    for a in (P1, P2, L1, L2):
        a[~vis] = np.nan
    sel = pp.SelectedObs((t * 1e9).astype(np.int64), [f"G{j + 1:02d}" for j in range(ns)], P1, P2, L1, L2,
                         np.full((ne_, ns), 45.0), np.full((ne_, ns), 40.0), np.zeros((ne_, ns), dtype=bool),
                         np.where(vis, 0, -1), {0: ("C1W", "L1C", "C2W", "L2W")}, F1, F2)
    return sel, t, el, vis


def make_ar(seed=13, sig_p=0.15, sig_l=0.0015):
    """Dual-frequency synthetic data with integer N1/N2 per arc and fractional receiver phase biases,
    returning (ObsSet, ArcInfo, SelectedObs, truth) for PPP-AR tests."""
    obs, arcs, truth = make_obs(seed=seed)
    rng = np.random.default_rng(seed + 100)
    ne, ns = obs.el.shape
    vis = obs.usable.copy()
    lam1, lam2 = C / F1, C / F2
    g = (F1 / F2) ** 2
    geo = obs.rho + obs.mh * obs.zhd0[:, None] + obs.mw * truth["zwd"][:, None] + \
        obs.mg * (truth["gn"][:, None] * np.cos(obs.az) + truth["ge"][:, None] * np.sin(obs.az)) + truth["clk"][:, None]
    I1 = 0.162 * 25.0 / np.sin(np.maximum(obs.el, 0.1))
    N1 = np.zeros((ne, ns))
    N2 = np.zeros((ne, ns))
    ints = {}
    for a, meta in arcs.arc_meta.items():
        n1, n2 = int(rng.integers(-5000, 5000)), int(rng.integers(-5000, 5000))
        m = arcs.arc == a
        N1[m], N2[m] = n1, n2
        ints[a] = (n1, n2)
    phi1, phi2 = 0.31, -0.17                      # receiver phase biases (cycles), common to all satellites
    w = 1 + 1 / np.sin(np.maximum(obs.el, 0.1))
    P1 = geo + I1 + 1.3 + rng.normal(0, 1, (ne, ns)) * sig_p * w
    P2 = geo + g * I1 + 2.1 + rng.normal(0, 1, (ne, ns)) * sig_p * w
    L1 = geo - I1 + lam1 * (N1 + phi1) + rng.normal(0, 1, (ne, ns)) * sig_l * w
    L2 = geo - g * I1 + lam2 * (N2 + phi2) + rng.normal(0, 1, (ne, ns)) * sig_l * w
    for a_ in (P1, P2, L1, L2):
        a_[~vis] = np.nan
    obs.Pif = A1 * P1 - A2 * P2
    obs.Lif = A1 * L1 - A2 * L2
    sel = pp.SelectedObs(obs.epochs, obs.sats, P1, P2, L1, L2, np.full((ne, ns), 45.0), np.full((ne, ns), 40.0),
                         np.zeros((ne, ns), dtype=bool), np.where(vis, 0, -1), {0: ("C1W", "L1W", "C2W", "L2W")},
                         F1, F2)
    truth["ints"] = ints
    return obs, arcs, sel, truth
