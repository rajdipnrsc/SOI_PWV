"""Modelled observations: assembles geometry and all deterministic corrections of TDS § 7 term by term.

(Added to the § 33.2 layout as one small module so that corrections.py stays a library of pure functions;
docs/decisions.md D-002.)  Also provides the code-only single-point solution (SPP) used for the a priori
receiver clock (reception time, TDS § 7.1) and for the a priori position when no SOI coordinate exists.
"""
from dataclasses import dataclass, field

import numpy as np

from . import antenna as an
from . import coords as co
from . import corrections as cr
from . import log as plog
from . import orbclk
from . import settings
from . import timesys as ts

LOG = plog.get()
C = settings.C_LIGHT


@dataclass
class Model:
    valid: np.ndarray            # (ne, ns) model available
    rho: np.ndarray              # modelled range incl. all non-tropospheric, non-clock corrections (m)
    dts: np.ndarray              # c * satellite clock incl. relativistic term (m)
    el: np.ndarray               # rad
    az: np.ndarray               # rad
    los: np.ndarray              # (ne, ns, 3) unit vector receiver -> satellite (ECEF)
    windup: np.ndarray           # cycles
    terms: dict = field(default_factory=dict)
    excluded: dict = field(default_factory=dict)
    eclipse_epochs: np.ndarray = None
    orbit_edge_epochs: np.ndarray = None
    clock_interp_epochs: np.ndarray = None
    attitude_source: str = "nominal"
    disp: np.ndarray = None      # (ne, 3) station displacement applied (m)


def station_displacements(X, epochs, erp=None, blq=None, want=("solid", "otl", "pole")):
    """Solid tide + ocean loading + pole tide displacement (ECEF, m) at each epoch (TDS § 7.10-7.13)."""
    lat, lon, _ = co.ecef_to_geodetic(X)
    utc = np.array([ts.gpst_to_utc(int(t)) for t in epochs], dtype=np.int64)
    out = {}
    d = np.zeros((len(epochs), 3))
    if erp is not None:
        xp, yp, dut = orbclk.erp_interp(erp, ts.mjd(utc))
    else:
        xp = yp = dut = np.zeros(len(epochs))
    if "solid" in want:
        sun, moon = cr.sun_moon_ecef(epochs, dut)
        out["solid"] = cr.solid_tide(X, sun, moon, utc)
        d += out["solid"]
    if "otl" in want and blq is not None:
        dusw = cr.ocean_loading(blq[0], blq[1], utc)
        out["otl"] = cr.otl_enu_to_xyz(dusw, lat, lon)
        d += out["otl"]
    if "pole" in want and erp is not None:
        dy = np.array([ts.decimal_year(int(t)) for t in epochs])
        out["pole"] = np.array([cr.pole_tide(lat, lon, dy[i], xp[i], yp[i])[0] for i in range(len(epochs))])
        d += out["pole"]
    return d, out


def compute(epochs, sats, X0, rclk_s, prod, atx, rcv_ant, disp, att_convention_check=True, windup_on=True):
    """Model all satellites at all epochs for receiver ARP X0 (ECEF, m).

    rclk_s: (ne,) a priori receiver clock (s) used for the reception time; prod: ProductSet-like object with
    sp3, clk, att; atx: Antex; rcv_ant: receiver AntennaEntry or None; disp: (ne,3) station displacement.
    """
    ne, ns = len(epochs), len(sats)
    lat, lon, _ = co.ecef_to_geodetic(X0)
    R = co.rot_xyz2enu(lat, lon)
    valid = np.zeros((ne, ns), dtype=bool)
    rho = np.full((ne, ns), np.nan)
    dts = np.full((ne, ns), np.nan)
    el = np.full((ne, ns), np.nan)
    az = np.full((ne, ns), np.nan)
    los = np.full((ne, ns, 3), np.nan)
    wu = np.zeros((ne, ns))
    terms = {k: np.full((ne, ns), np.nan) for k in ("geometric", "sagnac", "relativity", "shapiro", "rcv_antenna",
                                                     "sat_pcv", "sat_pco", "tides", "windup_m")}
    excluded = {"no_orbit": 0, "no_clock": 0, "no_antex": 0, "maneuver": 0, "eclipse": 0}
    ecl_ep = np.zeros(ne, dtype=bool)
    edge_ep = np.zeros(ne, dtype=bool)
    cint_ep = np.zeros(ne, dtype=bool)
    Xr = X0[None, :] + disp                                         # (ne,3)
    t_rx = epochs - np.round(np.asarray(rclk_s) * ts.NS).astype(np.int64)
    # Sun at every epoch: interpolating a 10-min subsample caused up to 1.7 mm range error through the nominal yaw
    # near noon/midnight turns (Level-0 decomposition test); the full computation costs ~40 ms per day
    sun, _ = cr.sun_moon_ecef(epochs)
    a1, a2 = an.if_coefficients("G", settings.BANDS["G"])
    lam_nl = C / (settings.FREQ[("G", "1")] + settings.FREQ[("G", "2")])
    # receiver epochs not aligned to clock epochs -> CLOCK_INTERP flag (TDS § 13.4)
    civ = int(round(prod.clk.interval * ts.NS)) if prod.clk is not None else 30 * ts.NS
    cint_ep[:] = (epochs % civ) != 0
    att_src = "nominal"
    use_att = prod.att is not None
    conv_checked = False
    for j, sat in plog.progress(list(enumerate(sats)), desc="Modelling satellites", unit="sat"):
        if sat not in prod.sp3.sats:
            excluded["no_orbit"] += 1
            continue
        ent = an.satellite_entry(atx, sat, int(epochs[ne // 2]))
        if ent is None or an.pco_if(ent) is None:
            excluded["no_antex"] += 1
            continue
        tau = np.full(ne, 0.075)
        for _ in range(3):
            t_tx = t_rx - np.round(tau * ts.NS).astype(np.int64)
            pos, vel, ok, edge = orbclk.sp3_interp(prod.sp3, sat, t_tx, want_vel=True)
            rs = cr.sagnac_rotate(pos, tau)
            tau_new = np.linalg.norm(rs - Xr, axis=1) / C
            if np.nanmax(np.abs(tau_new - tau)) < settings.LIGHT_TIME_TOL_S * 1e3:
                tau = tau_new
                break
            tau = np.where(np.isfinite(tau_new), tau_new, tau)
        t_tx = t_rx - np.round(tau * ts.NS).astype(np.int64)
        pos, vel, ok, edge = orbclk.sp3_interp(prod.sp3, sat, t_tx, want_vel=True)
        ck, cok, cip = orbclk.clk_interp(prod.clk, sat, t_tx)
        good = ok & cok
        excluded["no_clock"] += int(np.sum(ok & ~cok))
        if not good.any():
            continue
        # maneuver flags (SP3 event flags) -> exclude whole day for that satellite
        jj = prod.sp3.sats.index(sat)
        if prod.sp3.maneuver[:, jj].any():
            excluded["maneuver"] += 1
            continue
        # attitude: ORBEX (same family) else nominal + eclipse exclusion
        ex, ey, ez = cr.nominal_attitude(pos, sun)
        if use_att and sat in prod.att.sats:
            if att_convention_check and not conv_checked:
                _check_att_convention(prod.att, sat, t_tx, ex, ey, ez, good)
                conv_checked = True
            exn = ex.copy()
            got = np.zeros(ne, dtype=bool)
            for k in np.nonzero(good)[0]:
                axes = orbclk.att_axes(prod.att, sat, int(t_tx[k]))
                if axes is not None:
                    ex[k], ey[k], ez[k] = axes
                    got[k] = True
            # body-frame convention check (IGS vs manufacturer frame): outside eclipse ORBEX x must agree with
            # the nominal (IGS) x-axis; a systematic 180-deg yaw offset is removed (logged)
            beta, _mu = cr.beta_and_mu(pos, vel, sun)
            chk = got & (np.abs(beta) > np.radians(15.0))
            if chk.sum() > 10 and np.median(np.sum(ex[chk] * exn[chk], axis=1)) < -0.5:
                ex[got], ey[got] = -ex[got], -ey[got]
                LOG.info("ORBEX %s: x/y axes rotated by 180 deg to the IGS body-frame convention", sat)
            att_src = "ORBEX"
        else:
            beta, mu = cr.beta_and_mu(pos, vel, sun)
            shadow = cr.in_shadow(pos, sun)
            excl = cr.eclipse_exclusion(an.block_of(ent), beta, mu, shadow, t_tx)
            excl &= good
            if excl.any():
                excluded["eclipse"] += int(excl.sum())
                ecl_ep |= excl
                good &= ~excl
        # satellite PCO (IF) -> antenna phase centre
        p = an.pco_if(ent)
        apc = pos + ex * p[0] + ey * p[1] + ez * p[2]
        rs = cr.sagnac_rotate(apc, tau)
        vec = rs - Xr
        r = np.linalg.norm(vec, axis=1)
        e = vec / r[:, None]
        enu = e @ R.T
        elev = np.arcsin(np.clip(enu[:, 2], -1, 1))
        azim = np.arctan2(enu[:, 0], enu[:, 1]) % (2 * np.pi)
        # satellite PCV (nadir angle), receiver antenna, Shapiro, relativity
        nad = np.degrees(np.arccos(np.clip(np.sum(-e * ez, axis=1), -1, 1)))
        spcv = an.pcv_if(ent, nad)
        rant = an.receiver_correction(rcv_ant, enu) if rcv_ant is not None else np.zeros(ne)
        shp = cr.shapiro(rs, Xr)
        rel = cr.relativistic_clock(pos, vel)
        rho_j = r + rant + spcv + shp
        wind = np.zeros(ne)
        if windup_on and good.any():
            idx = np.nonzero(good)[0]
            wind[idx] = cr.windup(ex[idx], ey[idx], rs[idx], Xr[idx], lat, lon)
        valid[:, j] = good & np.isfinite(rho_j)
        rho[:, j] = rho_j
        dts[:, j] = C * (ck + rel)
        el[:, j] = elev
        az[:, j] = azim
        los[:, j] = e
        wu[:, j] = wind
        edge_ep |= edge & good
        terms["geometric"][:, j] = np.linalg.norm(pos - X0[None, :], axis=1)
        terms["sagnac"][:, j] = np.linalg.norm(cr.sagnac_rotate(pos, tau) - X0[None, :], axis=1) - \
            terms["geometric"][:, j]
        terms["sat_pco"][:, j] = np.linalg.norm(rs - X0[None, :], axis=1) - \
            np.linalg.norm(cr.sagnac_rotate(pos, tau) - X0[None, :], axis=1)
        terms["tides"][:, j] = r - np.linalg.norm(rs - X0[None, :], axis=1)
        terms["relativity"][:, j] = -C * rel
        terms["shapiro"][:, j] = shp
        terms["rcv_antenna"][:, j] = rant
        terms["sat_pcv"][:, j] = spcv
        terms["windup_m"][:, j] = lam_nl * wind
    m = Model(valid, rho, dts, el, az, los, wu, terms, excluded, ecl_ep, edge_ep, cint_ep, att_src, disp)
    return m


def _check_att_convention(att, sat, t_tx, ex, ey, ez, good):
    """Resolve the ORBEX quaternion direction with data (TDS § 7.11 VERIFY): outside eclipse the ORBEX
    body z-axis must agree with the nominal z-axis (toward Earth centre). Sets att.convention."""
    idx = np.nonzero(good)[0][::40][:20]
    best = None
    for conv in ("T2B", "B2T"):
        att.convention = conv
        ang = []
        for k in idx:
            axes = orbclk.att_axes(att, sat, int(t_tx[k]))
            if axes is None:
                continue
            ang.append(np.degrees(np.arccos(np.clip(np.dot(axes[2], ez[k]), -1, 1))))
        if ang and (best is None or np.median(ang) < best[1]):
            best = (conv, float(np.median(ang)))
    if best is None or best[1] > 5.0:
        LOG.warning("ATT_CONVENTION_UNRESOLVED: ORBEX z-axis disagrees with nominal by %.1f deg; ORBEX ignored",
                    best[1] if best else -1)
        att.sats.clear()
    else:
        att.convention = best[0]
        LOG.info("ORBEX quaternion convention resolved from data: %s (z-axis vs nominal %.3f deg)", *best)


# ============================================================================ single point positioning
def spp(sel, epochs, prod, X_apriori=None, max_epochs=240):
    """Code-only ionosphere-free single point solution with precise orbits/clocks.

    Returns (X (3,), per-epoch receiver clock (s) array (NaN where unavailable), rms_m).
    A simple a priori troposphere (2.3 m / sin e) is removed; the position is the median of epoch solutions.
    """
    ne, ns = sel.P1.shape
    a1 = sel.f1 ** 2 / (sel.f1 ** 2 - sel.f2 ** 2)
    a2 = sel.f2 ** 2 / (sel.f1 ** 2 - sel.f2 ** 2)
    Pif = a1 * sel.P1 - a2 * sel.P2
    X = np.array(X_apriori if X_apriori is not None and np.linalg.norm(X_apriori) > 6e6 else [0.0, 0.0, 0.0])
    step = max(1, ne // max_epochs)
    idx = np.arange(0, ne, step)

    def sat_state(t_rx_ns, j, Xr):
        tau = 0.075
        for _ in range(3):
            pos, vel, ok, _e = orbclk.sp3_interp(prod.sp3, sel.sats[j], np.array([t_rx_ns - int(tau * 1e9)]),
                                                 want_vel=True)
            if not ok[0]:
                return None, None
            rs = cr.sagnac_rotate(pos[0], tau)
            tau = np.linalg.norm(rs - Xr) / C
        ck, cok, _ = orbclk.clk_interp(prod.clk, sel.sats[j], np.array([t_rx_ns - int(tau * 1e9)]))
        if not cok[0]:
            return None, None
        return rs, ck[0] + cr.relativistic_clock(pos[0], vel[0])

    sols = []
    for it in range(2):
        sols = []
        for k in plog.progress(idx, desc=f"Single-point solution {it + 1}/2", unit="epoch"):
            m = np.isfinite(Pif[k]) & np.array([s in prod.sp3.sats for s in sel.sats])
            js = np.nonzero(m)[0]
            if len(js) < 5:
                continue
            Xk = X.copy()
            dt = 0.0
            for _ in range(8):
                H, y = [], []
                for j in js:
                    rs, ck = sat_state(int(epochs[k] - dt * 1e9), j, Xk)
                    if rs is None:
                        continue
                    vec = rs - Xk
                    r = np.linalg.norm(vec)
                    lat, lon, h = co.ecef_to_geodetic(Xk) if np.linalg.norm(Xk) > 6e6 else (0, 0, 0)
                    elv = np.arcsin(np.dot(vec / r, Xk / np.linalg.norm(Xk))) if np.linalg.norm(Xk) > 6e6 else np.pi / 2
                    if np.linalg.norm(Xk) > 6e6 and elv < np.radians(10):
                        continue
                    trop = 2.3 / max(np.sin(elv), 0.1) if np.linalg.norm(Xk) > 6e6 else 0.0
                    y.append(Pif[k, j] - (r + C * dt - C * ck + trop))
                    H.append(np.r_[-vec / r, 1.0])
                if len(y) < 5:
                    break
                H, y = np.array(H), np.array(y)
                dx, *_ = np.linalg.lstsq(H, y, rcond=None)
                Xk = Xk + dx[:3]
                dt += dx[3] / C
                if np.linalg.norm(dx[:3]) < 1e-3:
                    break
            if len(y) >= 5:
                res = y - H @ dx
                sols.append((k, Xk, dt, np.sqrt(np.mean(res ** 2))))
        if not sols:
            break
        X = np.median(np.array([s[1] for s in sols]), axis=0)
    if not sols:
        return None, np.full(ne, np.nan), np.nan
    # per-epoch clock with fixed position (all epochs)
    clk = np.full(ne, np.nan)
    ks = np.array([s[0] for s in sols])
    cs = np.array([s[2] for s in sols])
    clk[:] = np.interp(np.arange(ne), ks, cs)
    rms = float(np.median([s[3] for s in sols]))
    return X, clk, rms


def receiver_clock_from_code(Pif, rho, dts, tropo_slant, valid):
    """Robust (median over satellites) code-derived receiver clock per epoch (m)."""
    r = Pif - rho + dts - tropo_slant
    r[~valid] = np.nan
    with np.errstate(all="ignore"):
        return np.nanmedian(r, axis=1)
