"""Level-0: solid tide, ocean loading, pole tide, Sun/Moon, relativity, Shapiro, Sagnac, wind-up
(TDS § 7, § 30.1), against official IERS routine outputs and RTKLIB where required."""
import os

import numpy as np
import pytest

from ppp import corrections as cr
from ppp import coords as co
from ppp import settings
from ppp import timesys as ts
from tests.conftest import DATA


def test_solid_tide_iers_published_cases():
    cases = [((4075578.385, 931852.890, 4801570.154), (137859926952.015, 54228127881.4350, 23509422341.6960),
              (-179996231.920342, -312468450.131567, -169288918.592160), (2009, 4, 13),
              (0.07700420357108125891, 0.06304056321824967613, 0.05516568152597246810)),
             ((1112189.660, -4842955.026, 3985352.284), (-54537460436.2357, 130244288385.279, 56463429031.5996),
              (300396716.912, 243238281.451, 120548075.939), (2012, 7, 13),
              (-0.02036831479592075833, 0.05658254776225972449, -0.07597679676871742227)),
             ((1112200.5696, -4842957.8511, 3985345.9122), (100210282451.6279, 103055630398.3160, 56855096480.4475),
              (369817604.4348, 1897917.5258, 120804980.8284), (2015, 7, 15),
              (.00509570869172363845, .0828663025983528700, -.0636634925404189617))]
    for xs, su, mo, (y, m, d), exp in cases:
        r = cr.solid_tide(np.array(xs), np.array([su]), np.array([mo]), np.array([ts.to_ns(y, m, d)]))
        assert np.max(np.abs(r[0] - exp)) < 1e-7


def test_solid_tide_vs_iers_routine_1000_cases():
    """10 stations x 100 epochs vs official DEHANTTIDEINEL output: max < 0.1 mm (TDS § 30.1)."""
    rows = np.loadtxt(os.path.join(DATA, "iers_dehanttideinel_ref.txt"), comments="#")
    errs = []
    for r in rows:
        d = cr.solid_tide(r[0:3], r[None, 3:6], r[None, 6:9], np.array([int(r[9])], dtype=np.int64))
        errs.append(np.max(np.abs(d[0] - r[10:13])))
    assert max(errs) < 1e-4


def _hardisp_cases():
    cases, cur = [], None
    for L in open(os.path.join(DATA, "iers_hardisp_ref.txt")):
        if L.startswith("#"):
            continue
        if L.startswith("@"):
            p = L.split()
            cur = {"t": tuple(int(x) for x in p[2:8]), "n": int(p[8]), "step": int(p[9]), "blq": [], "out": []}
            cases.append(cur)
        elif len(cur["blq"]) < 6:
            cur["blq"].append([float(x) for x in L.split()])
        elif L.strip():
            cur["out"].append([float(x) for x in L.split()])
    return cases


def test_ocean_loading_vs_hardisp():
    """HARDISP algorithm port vs official program output: max < 0.1 mm (TDS § 30.1)."""
    for c in _hardisp_cases():
        blq = np.array(c["blq"])
        t0 = ts.to_ns(*c["t"][:3], c["t"][3], c["t"][4], c["t"][5])
        t = t0 + np.arange(c["n"]) * c["step"] * ts.NS
        d = cr.ocean_loading(blq[:3], blq[3:], t)
        assert np.max(np.abs(d - np.array(c["out"]))) < 1e-4


def test_blq_reader(tmp_path):
    p = tmp_path / "x.blq"
    p.write_text("$$ header\n  HYDE\n$$ Computed by OLMPP\n" +
                 "\n".join(" " + " ".join(f"{0.001 * (i + 1):.5f}" for _ in range(11)) for i in range(3)) + "\n" +
                 "\n".join(" " + " ".join(f"{10.0 * (i + 1):6.1f}" for _ in range(11)) for i in range(3)) + "\n")
    amp, ph, hdr = cr.read_blq(str(p), "HYDE")
    assert amp.shape == (3, 11) and amp[2, 0] == 0.003 and ph[1, 5] == 20.0


def test_sun_moon_meeus_examples():
    lam, beta, dist = cr.moon_ecliptic(np.array([2448724.5]))      # Meeus example 47.a
    assert np.degrees(lam[0]) == pytest.approx(133.162655, abs=2e-6)
    assert np.degrees(beta[0]) == pytest.approx(-3.229126, abs=2e-6)
    assert dist[0] / 1e3 == pytest.approx(368409.7, abs=0.1)
    lon, R = cr.sun_ecliptic(np.array([2448908.5]))                # Meeus example 25.a
    assert np.degrees(lon[0]) == pytest.approx(199.90988, abs=2e-4)
    assert R[0] / settings.AU == pytest.approx(0.99766, abs=2e-5)


def test_sun_moon_ecef_geometry():
    t = np.array([ts.to_ns(2024, 3, 20, 12, 0, 0)])                   # near equinox, local noon at lon 0
    sun, moon = cr.sun_moon_ecef(t)
    s = sun[0] / np.linalg.norm(sun[0])
    assert abs(np.degrees(np.arcsin(s[2]))) < 0.5                     # declination ~ 0
    assert abs(np.degrees(np.arctan2(s[1], s[0]))) < 2.5              # sub-solar longitude ~ 0 (eq. of time)
    assert 3.5e8 < np.linalg.norm(moon[0]) < 4.1e8


def test_pole_tide_formula():
    lat, lon = np.radians(17.4), np.radians(78.55)
    d = cr.pole_tide(lat, lon, 2024.0, 0.2, 0.4)
    xs, ys = cr.secular_pole(2024.0)
    assert (xs, ys) == pytest.approx((55.0 + 1.677 * 24, 320.5 + 3.46 * 24))
    m1, m2 = 0.2 - xs * 1e-3, -(0.4 - ys * 1e-3)
    th = np.pi / 2 - lat
    sr = -33 * np.sin(2 * th) * (m1 * np.cos(lon) + m2 * np.sin(lon)) * 1e-3
    enu = co.xyz2enu(d[0], lat, lon)
    assert enu[2] == pytest.approx(sr, abs=1e-12) and abs(enu[2]) < 0.025


def test_relativity_shapiro_sagnac():
    rs = np.array([[20579808.3, 11885481.2, 12607999.5]])
    vs = np.array([[-1500.0, 2000.0, 1200.0]])
    assert cr.relativistic_clock(rs, vs)[0] == pytest.approx(-2 * np.sum(rs * vs) / settings.C_LIGHT ** 2)
    rr = np.array([[1208298.2, 5967329.1, 1895460.8]])
    shp = cr.shapiro(rs, rr)[0]
    assert 0.010 < shp < 0.020                                        # ~18-19 mm at low elevation [L]
    # Sagnac: exact rotation vs RTKLIB geodist linear term: agree < 0.1 mm
    tau = np.linalg.norm(rs - rr) / settings.C_LIGHT
    r_rot = np.linalg.norm(cr.sagnac_rotate(rs, tau) - rr)
    r_rtk = np.linalg.norm(rs - rr) + settings.OMEGA_E * (rs[0, 0] * rr[0, 1] - rs[0, 1] * rr[0, 0]) / settings.C_LIGHT
    assert abs(r_rot - r_rtk) < 1e-4


def test_windup_vs_rtklib():
    """Wind-up vs RTKLIB windupcorr (nominal attitude): < 0.001 cycle after removing the arc constant
    (the initial value is absorbed by the ambiguity; TDS § 7.9)."""
    rows = np.loadtxt(os.path.join(DATA, "rtklib_windup_ref.txt"), comments="#")
    starts = list(np.nonzero(rows[:, 8] == 1)[0]) + [len(rows)]
    worst = 0.0
    for a, b in zip(starts[:-1], starts[1:]):
        r = rows[a:b]
        t = (r[:, 0] * 7 * 86400 + r[:, 1]).astype(np.int64) * ts.NS
        rs, rr = r[:, 2:5], r[:, 5:8]
        sun = r[:, 10:13]                           # RTKLIB's own Sun -> isolates the wind-up formula
        ex, ey, ez = cr.nominal_attitude(rs, sun)
        lat, lon, _ = co.ecef_to_geodetic(rr[0])
        w = cr.windup(ex, ey, rs, rr, lat, lon)
        d = w - r[:, 9]
        worst = max(worst, np.max(np.abs(d - np.round(np.median(d)))))
    assert worst < 1e-3
    # NOTE: RTKLIB's own Sun differs from ours by ~0.33 deg in 2024 (= precession since J2000 applied to a
    # mean-of-date position inside RTKLIB 2.4.3 sunmoonpos); our Sun is checked against Meeus instead.


def test_sun_radec_and_gmst_meeus():
    jde = np.array([2448908.5])                                    # Meeus example 25.a (apparent values)
    lon, R = cr.sun_ecliptic(jde)
    eps = cr.mean_obliquity(jde)
    ra = np.degrees(np.arctan2(np.cos(eps) * np.sin(lon), np.cos(lon))) % 360
    dec = np.degrees(np.arcsin(np.sin(eps) * np.sin(lon)))
    assert ra[0] == pytest.approx(198.38083, abs=0.02)           # 13h13m31.4s (apparent; aberr.+nutation ~0.01 deg)
    assert dec[0] == pytest.approx(-7.78361, abs=0.02)
    # Meeus example 12.a: 1987 April 10, 0h UT -> GMST 13h10m46.3668s
    assert np.degrees(cr.gmst(2446895.5)) == pytest.approx(197.693195, abs=1e-6)


def test_eclipse_shadow():
    sun = np.array([[1.5e11, 0, 0]])
    assert cr.in_shadow(np.array([[-2.6e7, 0, 0]]), sun)[0]
    assert not cr.in_shadow(np.array([[2.6e7, 0, 0]]), sun)[0]
    assert not cr.in_shadow(np.array([[-2.6e7, 1.0e7, 0]]), sun)[0]
