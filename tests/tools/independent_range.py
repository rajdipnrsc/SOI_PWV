"""TEST-ONLY independent implementation of the GPS IF modelled range (TDS § 30.1 'Full modelled range').

Written separately from ppp/model.py (different algorithms where possible): satellite positions/clocks directly
from the broadcast ephemeris polynomials (no SP3/CLK interpolation), velocity by central differences, light time in a
quasi-inertial frame aligned with ECEF at reception, attitude from an explicit yaw-steering construction, antenna
models from analytic definitions of a synthetic ANTEX file.  Not imported by the program.
"""
import numpy as np

from ppp import orbclk

C = 299792458.0
OMEGA = 7.2921151467e-5
GM = 3.986004418e14


def _rot_z(v, ang):
    """Passive rotation R3(ang): coordinates of a fixed vector in axes rotated by +ang about z."""
    c, s = np.cos(ang), np.sin(ang)
    return np.array([c * v[0] + s * v[1], -s * v[0] + c * v[1], v[2]])


def ephemeris_source(eph):
    """Continuous truth from the broadcast polynomials: t_ns -> (ECEF position m, clock s)."""
    return lambda t_ns: orbclk.broadcast_pos_clk(eph, int(t_ns), include_rel=False)


def node_source(node_t, node_pos, clk_t, clk_v, n=10):
    """Interpolation of product nodes with scipy's barycentric Lagrange form (n+1 nearest nodes) and linear clock
    interpolation - an implementation independent of orbclk.sp3_interp / clk_interp."""
    from scipy.interpolate import BarycentricInterpolator
    tn = (np.asarray(node_t) - node_t[0]) / 1e9

    def f(t_ns):
        tq = (t_ns - node_t[0]) / 1e9
        k = int(np.clip(np.searchsorted(tn, tq) - (n + 1) // 2, 0, len(tn) - n - 1))
        sl = slice(k, k + n + 1)
        pos = np.array([BarycentricInterpolator(tn[sl], node_pos[sl, i])(tq) for i in range(3)], dtype=float)
        return pos, float(np.interp(t_ns, clk_t, clk_v))
    return f


def model_one(src, t_rx_ns, Xr, sun, sat_pco_body, sat_pcv_fn, rcv_fn, enu_rot):
    """Return dict of terms (m) for one satellite and reception time (GPST ns) at receiver position Xr (ECEF).
    src: t_ns -> (CoM position, clock s), e.g. ephemeris_source or node_source."""
    _sat = lambda _e, t: src(t)  # noqa: E731
    eph = None
    # light time in the inertial frame I that coincides with ECEF at t_rx: ECEF(t) = R3(theta(t)) I with
    # theta(t_rx) = 0, theta(t_tx) = -omega*tau, hence X_I(t_tx) = R3(omega*tau) x_ECEF(t_tx)
    tau = 0.07
    for _ in range(10):
        r_e, _ = _sat(eph, t_rx_ns - int(round(tau * 1e9)))
        r_i = _rot_z(r_e, OMEGA * tau)
        tau_new = np.linalg.norm(r_i - Xr) / C
        if abs(tau_new - tau) < 1e-13:
            tau = tau_new
            break
        tau = tau_new
    t_tx = t_rx_ns - int(round(tau * 1e9))
    r_e, clk = _sat(eph, t_tx)
    h = 0.5
    rp, _ = _sat(eph, t_tx + int(h * 1e9))
    rm, _ = _sat(eph, t_tx - int(h * 1e9))
    v_e = (rp - rm) / (2 * h)
    # yaw-steering axes (IGS): z to Earth centre, y = z x (unit vector to Sun), x completes the right-handed set
    z = -r_e / np.linalg.norm(r_e)
    s = (sun - r_e) / np.linalg.norm(sun - r_e)
    y = np.cross(z, s)
    y /= np.linalg.norm(y)
    x = np.cross(y, z)
    apc_e = r_e + sat_pco_body[0] * x + sat_pco_body[1] * y + sat_pco_body[2] * z
    apc_i = _rot_z(apc_e, OMEGA * tau)
    geo_i = np.linalg.norm(apc_i - Xr)
    los = (apc_i - Xr) / geo_i
    nadir = np.degrees(np.arccos(np.clip(np.dot(-los, z), -1, 1)))
    enu = enu_rot @ los
    rs, rr = np.linalg.norm(apc_i), np.linalg.norm(Xr)
    return {"range_geometry": geo_i, "sat_pcv": sat_pcv_fn(nadir), "rcv_antenna": rcv_fn(enu),
            "shapiro": 2 * GM / C ** 2 * np.log((rs + rr + geo_i) / (rs + rr - geo_i)),
            "sat_clock": C * clk, "relativity": -2.0 * np.dot(r_e, v_e) / C,      # c * dt_rel
            "elevation": np.arcsin(enu[2])}
