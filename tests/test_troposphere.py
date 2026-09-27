"""Level-0: mapping functions and a priori delays (TDS § 6)."""
import numpy as np
import pytest

from ppp import troposphere as tr


def test_gmf_iers_test_vector():
    mh, mw = tr.gmf(55055.0, 0.6708665767, -1.393397187, 844.715, np.pi / 2 - 1.278564131)
    assert mh == pytest.approx(3.425245519339138678, abs=1e-12)
    assert mw == pytest.approx(3.449589116182419257, abs=1e-12)


def test_vmf3_properties():
    lat, lon = np.radians(17.4), np.radians(78.55)
    mh, mw = tr.vmf3_ht(0.00123, 0.00056, 60324.5, lat, lon, 420.0, np.radians([90.0, 30.0, 7.0]))
    assert mh[0] == pytest.approx(1.0, abs=1e-12) and mw[0] == pytest.approx(1.0, abs=1e-12)
    assert 1.99 < mh[1] < 2.0 and 7.0 < mh[2] < 8.5 and mw[2] > mh[2]
    bh, bw, ch, cw = tr.vmf3_bc(60324.5, lat, lon)
    assert 0.002 < bh < 0.0035 and 0.001 < bw < 0.002 and 0.04 < ch < 0.08 and 0.03 < cw < 0.06


def test_saastamoinen_and_inverse():
    z = tr.saastamoinen_zhd(1013.25, np.radians(45.0), 0.0)
    assert z == pytest.approx(0.0022768 * 1013.25, rel=1e-12)
    assert tr.saastamoinen_pressure(z, np.radians(45.0), 0.0) == pytest.approx(1013.25)


def test_gradient_mf():
    el = np.radians(10.0)
    assert tr.chen_herring_mg(el) == pytest.approx(1.0 / (np.sin(el) * np.tan(el) + 0.0032))


def test_kouba_height_reduction_vs_linearised():
    """Kouba (2008) power-law pressure reduction vs PRIDE's linearised standard-atmosphere form: the two
    conventions differ by ~1 % of the pressure change, i.e. < 1.5 mm ZHD for |dh| <= 300 m (verify_log V-019)."""
    lat = np.radians(17.4)
    zg, hg = 2.20, 500.0
    for hs in (250.0, 500.0, 800.0):
        pg = tr.saastamoinen_pressure(zg, lat, hg)
        ps = pg * (1 - 0.0000226 * (hs - hg)) ** 5.225
        lin = pg - 1013.25 * 5.225 * (1 - 0.226e-4 * hg) ** 4.225 * 0.226e-4 * (hs - hg)
        assert abs(tr.saastamoinen_zhd(ps, lat, hs) - tr.saastamoinen_zhd(lin, lat, hs)) < 1.5e-3


def test_gpt3_reader_reduced_and_full(tmp_path):
    header = "%  lat    lon   p:a0    A1   B1   A2   B2  T:a0   A1   B1   A2   B2  Q:a0   A1   B1   A2   B2  " \
             "dT:a0   A1   B1   A2   B2  undu      Hs   a_h:a0   A1   B1   A2   B2  a_w:a0   A1   B1   A2   B2  " \
             "lambda:a0  A1   B1   A2   B2  Tm:a0   A1   B1   A2   B2\n"
    rows = []
    for la in np.arange(89.5, -90, -1.0):
        for lo in np.arange(0.5, 360, 1.0):
            rows.append(f"{la} {lo} 101300 0 0 0 0 290 0 0 0 0 10 0 0 0 0 -6.5 0 0 0 0 -70 0 1.2 0 0 0 0 "
                        f"0.55 0 0 0 0 3.0 0 0 0 0 285 0 0 0 0")
    p = tmp_path / "gpt3_1.grd"
    p.write_text(header + "\n".join(rows) + "\n")
    g = tr.read_gpt3(str(p))
    assert g.has_mf
    out = tr.gpt3(g, 60324.5, np.radians(17.4), np.radians(78.55), 0.0)
    assert out["ah"] == pytest.approx(0.0012) and out["aw"] == pytest.approx(0.00055)
    assert out["undu"] == pytest.approx(-70.0) and out["Tm"] == pytest.approx(285.0)
    assert 1003 < out["p"] < 1006                                    # reduced from Hs=0 to h_orth=70 m
