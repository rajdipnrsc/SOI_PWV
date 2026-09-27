"""Level-0: frame transformation, propagation, PMM, reference point, auto decision (TDS § 9, § 30.1)."""
import numpy as np
import pytest

from ppp import coords as co
from ppp import timesys as ts

X_HYDE = np.array([1208298.991, 5967329.098, 1895460.127])     # SOI HYDE, ITRF2008 @ 2005.0 (sample data)


def test_geodetic_roundtrip():
    lat, lon, h = co.ecef_to_geodetic(X_HYDE)
    assert np.allclose(co.geodetic_to_ecef(lat, lon, h), X_HYDE, atol=1e-6)
    assert np.degrees(lat) == pytest.approx(17.40201222, abs=1e-7)       # station_coord.csv values
    assert np.degrees(lon) == pytest.approx(78.55318912, abs=1e-7)
    assert h == pytest.approx(419.573, abs=1e-3)


def test_helmert_roundtrip_and_magnitude():
    key = ("ITRF2020", "ITRF2008")
    for t in (2005.0, 2015.0, 2026.0):
        X20 = co.helmert_inverse(X_HYDE, key, t)
        assert np.allclose(co.helmert_forward(X20, key, t), X_HYDE, atol=1e-7)
    # magnitude check from TDS § 9.5: at 2026.0 T = (0.2, -0.1, 4.4) mm
    T, D, R = co._helmert_params(key, 2026.0)
    assert np.allclose(T * 1e3, [0.2, -0.1, 4.4], atol=1e-9)
    assert D * 1e9 == pytest.approx(0.04, abs=1e-9)


def test_equivalence_of_propagation_orders():
    """'propagate in ITRF2008 then transform at t' == 'transform at t0 then propagate in ITRF2020' < 0.1 mm."""
    t0, t = 2005.0, 2024.5
    V08 = co.velocity_itrf2008_from_pmm(X_HYDE)
    a = co.helmert_inverse(X_HYDE + V08 * (t - t0), ("ITRF2020", "ITRF2008"), t)
    V20 = V08 - co.velocity_rate_terms(X_HYDE, ("ITRF2020", "ITRF2008"))
    b = co.helmert_inverse(X_HYDE, ("ITRF2020", "ITRF2008"), t0) + V20 * (t - t0)
    assert np.max(np.abs(a - b)) < 1e-4


def test_pmm_india_velocity():
    lat, lon, _ = co.ecef_to_geodetic(X_HYDE)
    venu = co.xyz2enu(co.pmm_velocity_itrf2020(X_HYDE), lat, lon) * 1e3
    # India plate at Hyderabad: ~40 mm/yr east, ~35 mm/yr north in ITRF2020 [L]
    assert 35 < venu[0] < 45 and 30 < venu[1] < 40 and abs(venu[2]) < 1.5


def test_soi_to_processing_and_refpoint():
    t = ts.to_ns(2024, 1, 15, 12)
    c = co.soi_to_processing(X_HYDE, t, "IGS20", (0.1, 0.0, 0.0), "MARKER")
    c_arp = co.soi_to_processing(X_HYDE, t, "IGS20", (0.1, 0.0, 0.0), "ARP")
    d = co.xyz2enu(c.xyz - c_arp.xyz, np.radians(c.lat), np.radians(c.lon))
    assert np.allclose(d, [0, 0, 0.1], atol=1e-6)            # DELTA H applied once, marker only
    assert c.frame == "ITRF2020" and c.tide_system == "TIDE_FREE" and c.reference_point == "ARP"
    moved = co.xyz2enu(c_arp.xyz - X_HYDE, np.radians(c.lat), np.radians(c.lon))
    assert 0.7 < moved[0] < 0.85 and 0.6 < moved[1] < 0.75 and abs(moved[2]) < 0.01
    with pytest.raises(KeyError):
        co.soi_to_processing(X_HYDE, t, "XYZ99", (0, 0, 0))


def test_igs14_era_chain():
    c14 = co.soi_to_processing(X_HYDE, ts.to_ns(2019, 1, 1), "IGS14", (0, 0, 0))
    assert c14.frame == "ITRF2014"
    assert co.frame_chain("ITRF2008", "ITRF2014") == [("inverse", ("ITRF2014", "ITRF2008"))]
    assert co.frame_chain("ITRF2008", "ITRF2008") == []


def test_auto_decision_rules():
    lat, lon, _ = co.ecef_to_geodetic(X_HYDE)
    up = co.enu2xyz(np.array([0, 0, 1.0]), lat, lon)
    XA = X_HYDE
    ok = co.auto_decision(X_HYDE + 0.01 * up, X_HYDE + 0.01 * up, XA, (0, 0, 0))
    assert ok["decision"] == "FIXED_SOI"
    bad = co.auto_decision(X_HYDE + 0.05 * up, X_HYDE + 0.05 * up, XA, (0, 0, 0))
    assert bad["decision"] == "KEEP_STATIC"
    # SOI coordinate was actually the ARP: marker-based interpretation is 0.5 m too high
    sw = co.auto_decision(X_HYDE + 0.5 * up, X_HYDE + 0.004 * up, XA, (0.5, 0, 0))
    assert sw["refpoint_switched"] and sw["decision"] == "FIXED_SOI"


def test_deforming_zones_and_permanent_tide():
    assert co.in_deforming_zone(17.4, 78.5) is None
    assert co.in_deforming_zone(11.6, 92.7) == "ANDAMAN_NICOBAR"
    dr, dn = co.permanent_tide_offset(np.radians(17.4))
    assert 0.03 < dr < 0.06                                  # several cm at Indian latitudes [L]
