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
