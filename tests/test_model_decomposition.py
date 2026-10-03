"""Level-0 'full modelled range' (TDS § 30.1): ppp/model.py vs an independent implementation, term by term.

Products are built from broadcast ephemerides (5-min SP3 nodes, 30-s clocks) so that the independent code can use
the continuous ephemeris itself as truth; a synthetic ANTEX gives analytic satellite and receiver antenna models.
Thresholds: < 1 mm per term, < 2 mm total."""
import os

import numpy as np

from ppp import antenna as an
from ppp import coords as co
from ppp import corrections as cr
from ppp import model as mdl
from ppp import rinex
from ppp import timesys as ts
from tests.conftest import SAMPLE
from tests.test_antenna import _fix, _freq_block
from tests.tools import broadcast_products as bp
from tests.tools import independent_range as ir

SATS = [f"G{k:02d}" for k in range(2, 33)]


def _atx(path):
    zen = np.arange(0, 91, 5.0)
    L = ["     1.4            M                                       ANTEX VERSION / SYST",
         "A                                                           PCV TYPE / REFANT",
         "                                                            END OF HEADER",
         "                                                            START OF ANTENNA",
         "TSTANT1         NONE                                        TYPE / SERIAL NO",
         "     5.0                                                    DAZI",
         "     0.0  90.0   5.0                                        ZEN1 / ZEN2 / DZEN"]
    # 5-deg azimuth rows so that the ANTEX grid represents the analytic cos(az) pattern to < 0.01 mm
    L += _freq_block("G01", (1.0, 2.0, 100.0), -0.1 * zen,
                     [(az, -0.1 * zen + 0.05 * zen * np.cos(np.radians(az))) for az in np.arange(0, 361, 5.0)])
    L += _freq_block("G02", (0.5, 1.0, 120.0), -0.2 * zen,
                     [(az, -0.2 * zen + 0.05 * zen * np.cos(np.radians(az))) for az in np.arange(0, 361, 5.0)])
    L += ["                                                            END OF ANTENNA"]
    nad = np.arange(0, 15, 1.0)
    for prn in SATS:
        L += ["                                                            START OF ANTENNA",
              f"BLOCK IIF           {prn}                 G0{prn[1:]}      2010-022A  TYPE / SERIAL NO",
              "     0.0                                                    DAZI",
              "     0.0  14.0   1.0                                        ZEN1 / ZEN2 / DZEN",
              "  2014     5    17     0     0    0.0000000                 VALID FROM"]
        for code in ("G01", "G02"):
            L += [f"   {code}                                                    START OF FREQUENCY",
                  f"{394.0:10.2f}{0.0:10.2f}{1500.0:10.2f}                              NORTH / EAST / UP",
                  "   NOAZI" + "".join(f"{v:8.2f}" for v in 0.5 * nad),
                  f"   {code}                                                    END OF FREQUENCY"]
        L += ["                                                            END OF ANTENNA"]
    path.write_text("\n".join(_fix(x) for x in L) + "\n")


def run_decomposition(tmp_path):
    """Max |model - independent| per term (m) over 2 h of HYDE 2024-015 geometry (every 10th epoch)."""
    import pathlib
    tmp_path = pathlib.Path(tmp_path)
    nav = rinex.read_nav(os.path.join(SAMPLE, "HYDE015.24n"))
    hdr = rinex.read_obs(os.path.join(SAMPLE, "HYDE015.24o"), decimate_s=3600).header
    X0 = np.array(hdr.approx_xyz)
    t0 = ts.to_ns(2024, 1, 15, 6)
    epochs = t0 + np.arange(0, 2 * 3600, 30, dtype=np.int64) * ts.NS
    ps = bp.build(nav, int(epochs[0]), int(epochs[-1]))
    p = tmp_path / "t.atx"
    _atx(p)
    atx = an.read_antex(str(p))
    rcv, _ = an.receiver_entry(atx, "TSTANT1         NONE")
    disp = np.tile(np.array([0.012, -0.021, 0.034]), (len(epochs), 1))         # tides etc. (m)
    rclk = np.full(len(epochs), 2.5e-4)
    sats = [s for s in SATS if s in nav]
    m = mdl.compute(epochs, sats, X0, rclk, ps, atx, rcv, disp, windup_on=False)
    # internal consistency: the stored terms add up to the modelled range
    tsum = sum(m.terms[k] for k in ("geometric", "sagnac", "sat_pco", "tides", "rcv_antenna", "sat_pcv", "shapiro"))
    ok = m.valid & (m.el > np.radians(5))
    assert ok.sum() > 1000
    assert np.nanmax(np.abs(tsum - m.rho)[ok]) < 1e-6
    # independent implementation
    lat, lon, _ = co.ecef_to_geodetic(X0)
    R = co.rot_xyz2enu(lat, lon)
    a1, a2 = an.if_coefficients()

    def rcv_fn(enu):
        zen = np.degrees(np.arccos(enu[2]))
        az = np.degrees(np.arctan2(enu[0], enu[1])) % 360
        pco = a1 * np.array([2.0, 1.0, 100.0]) - a2 * np.array([1.0, 0.5, 120.0])           # E, N, U (mm)
        pcv = a1 * (-0.1 * zen + 0.05 * zen * np.cos(np.radians(az))) - \
            a2 * (-0.2 * zen + 0.05 * zen * np.cos(np.radians(az)))
        return (-(enu @ pco) + pcv) * 1e-3
    sun, _ = cr.sun_moon_ecef(epochs)
    err = {k: [] for k in ("range_geometry", "sat_pcv", "rcv_antenna", "shapiro", "sat_clock", "relativity",
                           "total")}
    for j, sat in enumerate(sats):
        for k in np.nonzero(ok[:, j])[0][::10]:
            t_rx = int(epochs[k] - round(rclk[k] * 1e9))
            jj = ps.sp3.sats.index(sat)
            ct, cv = ps.clk.sats[sat]
            src = ir.node_source(ps.sp3.epochs, ps.sp3.pos[:, jj], ct, cv)
            ind = ir.model_one(src, t_rx, X0 + disp[k], sun[k], (0.394, 0.0, 1.5), lambda n: 0.5e-3 * n, rcv_fn, R)
            geo_model = m.terms["geometric"][k, j] + m.terms["sagnac"][k, j] + m.terms["sat_pco"][k, j] + \
                m.terms["tides"][k, j]
            err["range_geometry"].append(geo_model - ind["range_geometry"])
            err["sat_pcv"].append(m.terms["sat_pcv"][k, j] - ind["sat_pcv"])
            err["rcv_antenna"].append(m.terms["rcv_antenna"][k, j] - ind["rcv_antenna"])
            err["shapiro"].append(m.terms["shapiro"][k, j] - ind["shapiro"])
            err["relativity"].append(m.terms["relativity"][k, j] + ind["relativity"])   # range term = -c dt_rel
            err["sat_clock"].append(m.dts[k, j] + m.terms["relativity"][k, j] - ind["sat_clock"])
            tot_model = m.rho[k, j] - m.dts[k, j]
            tot_ind = ind["range_geometry"] + ind["sat_pcv"] + ind["rcv_antenna"] + ind["shapiro"] - \
                (ind["sat_clock"] + ind["relativity"])
            err["total"].append(tot_model - tot_ind)
    return {k: float(np.max(np.abs(v))) for k, v in err.items()}, int(sum(len(v) for v in err.values()) / len(err))


def test_full_modelled_range_decomposition(tmp_path):
    mx, n = run_decomposition(tmp_path)
    assert n > 100
    for k, v in mx.items():
        assert v < (2e-3 if k == "total" else 1e-3), (k, v)
