"""Stage-7 unit tests: PWV conversion constants and uncertainty propagation (TDS § 18)."""
import numpy as np
import pytest

from ppp import pwv as pw
from ppp import troposphere as tr


def test_pi_range():
    """Pi ~ 0.155-0.167 for Tm 270-295 K (Bevis 1994 constants) [TDS § 18.1]."""
    assert 0.150 < pw.pi_factor(270.0) < 0.160
    assert 0.160 < pw.pi_factor(295.0) < 0.169            # 0.1680 (TDS quotes "~0.167")
    assert pw.pi_factor(280.0) == pytest.approx(1e6 / (1000 * 461.5 * (3739.0 / 280.0 + 0.221)))


def test_conversion_and_budget():
    lat, h = np.radians(17.4), 497.0
    p = 960.0
    zhd = tr.saastamoinen_zhd(p, lat, h)
    ztd = np.array([zhd + 0.25])
    r = pw.convert(ztd, np.array([0.004]), lat, h, p, "VMF3_GRID", np.array([285.0]), "GPT3")
    assert r.zwd_met[0] == pytest.approx(250.0, abs=1e-6)
    assert r.pwv[0] == pytest.approx(250.0 * pw.pi_factor(285.0), rel=1e-12)
    assert r.sig_ppp[0] == pytest.approx(4.0 * pw.pi_factor(285.0), rel=1e-9)
    # conversion term: pressure (1 hPa -> 2.3 mm ZHD) and Tm (4 K) combined
    s_zhd = 0.0022768 / (1 - 0.00266 * np.cos(2 * lat) - 0.28e-6 * h) * 1.0 * 1e3
    s_pi = pw.sigma_pi(285.0, 4.0)
    assert r.sig_conv[0] == pytest.approx(np.hypot(pw.pi_factor(285.0) * s_zhd, 250.0 * s_pi), rel=1e-9)
    assert r.sig_total[0] == pytest.approx(np.hypot(r.sig_ppp[0], r.sig_conv[0]), rel=1e-9)


def test_pressure_reduction():
    assert pw.reduce_pressure(1000.0, 0.0, 100.0, 288.0) == pytest.approx(1000 * np.exp(-9.80665 * 100 / (287.05 * 288)))


def _met(p, td=None, hr=None, pos=None, dt_s=300):
    from ppp import rinex
    from ppp import timesys as ts
    t = ts.to_ns(2024, 1, 15) + np.arange(len(p), dtype=np.int64) * dt_s * 10 ** 9
    vals = {"PR": np.asarray(p, dtype=float)}
    if td is not None:
        vals["TD"] = np.asarray(td, dtype=float)
    if hr is not None:
        vals["HR"] = np.asarray(hr, dtype=float)
    return rinex.MetData(t, vals, {"PR": pos} if pos else {}, {}, 2.11, "TEST015A.24m")


def test_barometer_reduction_and_qc():
    p = np.full(20, 955.0)
    p[7] = 962.0                                     # spike -> rejected
    p[9] = 300.0                                     # implausible -> rejected
    m = _met(p, td=np.full(20, 20.0), hr=np.full(20, 50.0), pos=(0.0, 0.0, 0.0, 450.0))
    b, w = pw.barometer_from_met(m, h_ant_orth=500.0, undulation=-70.0)
    assert w == [] and b.height_source == "SENSOR_POS"
    assert b.h_sensor == pytest.approx(520.0)        # H_ell 450 - N (-70)
    assert len(b.t) == 18 and b.n_rejected == 2
    tv = pw.virtual_temperature(293.15 - 0.0065 * 0.5 * (500 - 520), 955.0, 50.0)
    assert b.p_ant[0] == pytest.approx(pw.reduce_pressure(955.0, 520.0, 500.0, tv))
    assert b.p_ant[0] > 955.0                       # antenna 20 m below the sensor
    # no sensor position: assumed at antenna height (warn, continue)
    b2, w2 = pw.barometer_from_met(_met(np.full(5, 955.0)), 500.0, None, 290.0)
    assert "MET_HEIGHT_ASSUMED" in w2 and "MET_NO_TEMPERATURE" in w2 and b2.p_ant == pytest.approx(955.0)
    # sensor far from the antenna -> not used
    b3, w3 = pw.barometer_from_met(_met(np.full(5, 955.0), pos=(0, 0, 0, 900.0)), 500.0, 0.0)
    assert b3 is None and "MET_SENSOR_TOO_FAR" in w3


def test_pressure_at_gaps():
    p = np.r_[np.full(6, 950.0), np.full(6, 951.0)]
    m = _met(p)
    m.t[6:] += 3600 * 10 ** 9                        # 1-h gap after sample 5
    b, _ = pw.barometer_from_met(m, 500.0, 0.0, 290.0)
    tq = m.t[0] + np.arange(0, 3 * 3600, 300, dtype=np.int64) * 10 ** 9
    q = pw.pressure_at(b, tq)
    assert q[0] == pytest.approx(950.0)
    assert np.isnan(q[12])                           # inside the gap (> MET_MAX_GAP_S)
    assert np.isfinite(q[5]) and np.isnan(q[-1])     # last sample + 30 min ... then outside the record


def test_convert_per_epoch_sources():
    lat, h = np.radians(17.4), 497.0
    z = tr.saastamoinen_zhd(960.0, lat, h)
    r = pw.convert(np.full(2, z + 0.2), np.full(2, 0.004), lat, h, 960.0, np.array(["BAROMETER", "GPT3"], object),
                   np.full(2, 285.0), np.array(["ERA5", "BEVIS"], object))
    assert list(r.p_source) == ["BAROMETER", "GPT3"] and "P_CLIMATOLOGY" in r.flags
    assert r.sig_conv[0] < r.sig_conv[1]


def test_apply_barometer_keeps_resolution_and_gaps():
    from ppp import timesys as ts
    t0 = ts.to_ns(2024, 1, 15)
    st = tr.StationTropo(np.array([t0, t0 + 6 * 3600 * 10 ** 9]), np.array([2.2, 2.2]), np.array([0.1, 0.1]),
                         np.array([1e-3] * 2), np.array([5e-4] * 2), np.radians(17.4), np.radians(78.5), 500.0,
                         "GPT3", "GPT3", {"ZHD_CLIMATOLOGY"})
    g = t0 + np.arange(0, 6 * 3600 + 1, 300, dtype=np.int64) * 10 ** 9
    p = 955.0 + 0.001 * np.arange(len(g))
    p[30:40] = np.nan
    st = tr.apply_barometer(st, g, p, 480.0)
    z, *_ = st.at(g)
    assert st.zhd_source == "BAROMETER" and st.zhd_source_fallback == "GPT3"
    assert z[10] == pytest.approx(tr.saastamoinen_zhd(p[10], st.lat, 480.0))
    assert z[35] == pytest.approx(2.2)              # gap: previous source kept
    assert "ZHD_CLIMATOLOGY" in st.flags            # not fully replaced


def test_era5_station_met():
    from ppp import mapping as mp
    from ppp import timesys as ts
    lev = np.array([1000.0, 850.0, 700.0, 500.0, 300.0])
    zl = np.array([110.0, 1500.0, 3100.0, 5800.0, 9600.0])
    Tl = np.array([300.0, 290.0, 280.0, 265.0, 235.0])
    ql = np.array([0.015, 0.010, 0.006, 0.002, 0.0002])
    t = ts.to_ns(2024, 1, 15) + np.array([0, 3600], dtype=np.int64) * 10 ** 9
    shp = (2, 5, 3, 3)
    bg = mp.Background(t, np.array([16.0, 17.0, 18.0]), np.array([77.0, 78.0, 79.0]), lev,
                       np.broadcast_to(zl[None, :, None, None], shp).copy(),
                       np.broadcast_to(Tl[None, :, None, None], shp).copy(),
                       np.broadcast_to(ql[None, :, None, None], shp).copy(), np.zeros((3, 3)))
    p, tm = pw.era5_station_met(bg, 17.4, 78.5, 500.0, t[0] + np.array([0, 1800, 7200], dtype=np.int64) * 10 ** 9)
    p_ref = np.exp(np.log(1000) + (500 - 110) * (np.log(850) - np.log(1000)) / (1500 - 110))
    assert p[0] == pytest.approx(p_ref, rel=1e-9) and p[1] == pytest.approx(p_ref, rel=1e-9)
    assert np.isnan(p[2])                            # outside the ERA5 hours
    assert 255 < tm[0] < 295
