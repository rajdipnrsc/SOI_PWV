"""Stage-8 tests: kriging maths, background column integrals, synthetic reconstruction, CLI (TDS Part B)."""
import os

import numpy as np
import pytest

from ppp import mapping as mp
from ppp import output as out
from ppp import settings
from ppp import timesys as ts

DAY = ts.to_ns(2024, 7, 15)


def test_uk_variance_matches_dense_formula():
    rng = np.random.default_rng(0)
    n = 25
    lat = rng.uniform(15, 20, n)
    lon = rng.uniform(75, 80, n)
    u = mp.to_unit(lat, lon)
    h = rng.uniform(0, 800, n)
    s = rng.uniform(0.5, 1.5, n)
    r = rng.normal(0, 2, n)
    cp = mp.CovParams(2.0, 150.0, 0.5)
    uq = mp.to_unit(np.array([17.3]), np.array([77.1]))
    rh, var = mp.krige(u, h, r, s, uq, np.array([300.0]), cp, k_max=100, r_max_km=5000)
    C = mp.cov_exp(mp.gc_km(u, u), 0, cp.sigma_s, cp.L_km) + np.diag(cp.nugget ** 2 + s ** 2)
    c = mp.cov_exp(mp.gc_km(uq, u)[0], 0, cp.sigma_s, cp.L_km)
    Ci = np.linalg.inv(C)
    F = np.ones(n)
    beta = F @ Ci @ r / (F @ Ci @ F)
    pred = beta + c @ Ci @ (r - beta)
    v = cp.sigma_s ** 2 - c @ Ci @ c + (1 - F @ Ci @ c) ** 2 / (F @ Ci @ F)
    assert rh[0] == pytest.approx(pred, rel=1e-6) and var[0] == pytest.approx(v, rel=1e-6)   # chord vs arccos distance rounding


def test_column_integral_analytic():
    p = np.array([1000.0, 850.0, 700.0, 500.0, 300.0, 100.0])
    z = np.array([100.0, 1450.0, 3000.0, 5600.0, 9200.0, 16000.0])
    T = np.array([300.0, 290.0, 282.0, 268.0, 240.0, 200.0])
    q = np.where(p >= 500, 0.01, 0.0)
    ps, pwv, tm = mp.column_quantities(z[None], T[None], q[None], p, 100.0)
    assert ps[0] == pytest.approx(1000.0)
    # trapezoid between 500 and 300 hPa adds half a layer of q
    assert pwv[0] == pytest.approx((0.01 * 50000 + 0.5 * 0.01 * 20000) / settings.G0, rel=1e-9)
    assert 250 < tm[0] < 300


def _synthetic_background():
    lat = np.arange(38.0, 5.75, -0.25)
    lon = np.arange(68.0, 98.25, 0.25)
    p = np.array([1000.0, 925, 850, 700, 500, 400, 300, 200, 100])
    zref = np.array([110.0, 760, 1460, 3010, 5580, 7180, 9160, 11800, 16200])
    nt = 25
    z = np.broadcast_to(zref[None, :, None, None], (nt, len(p), len(lat), len(lon))).copy()
    T = np.broadcast_to((300 - 0.0065 * zref)[None, :, None, None], z.shape).copy()
    qref = 0.018 * np.exp(-zref / 2000.0)
    q = np.broadcast_to(qref[None, :, None, None], z.shape).copy()
    t = DAY - 3600 * ts.NS + np.arange(nt) * 3600 * ts.NS
    return mp.Background(t, lat, lon, p, z, T, q, np.zeros((len(lat), len(lon))), None, "SYNTHETIC")


def _stations(n=60, seed=1, bg=None, bias_fn=None):
    rng = np.random.default_rng(seed)
    lat = rng.uniform(12, 30, n)
    lon = rng.uniform(72, 88, n)
    h = rng.uniform(0, 1500, n)
    st = []
    t = DAY + np.arange(288) * 300 * ts.NS
    for i in range(n):
        base, _, _ = mp.background_at(bg, lat[i:i + 1], lon[i:i + 1], h[i:i + 1], DAY) if bg else (
            [50 * np.exp(-h[i] / 2000)], None, None)
        truth = base[0] + (bias_fn(lat[i], lon[i]) if bias_fn else 0.0)
        st.append(mp.StationSeries(f"S{i:03d}", lat[i], lon[i], h[i], h[i], t, np.full(288, truth) +
                                   rng.normal(0, 0.3, 288), np.full(288, 0.8), np.full(288, 2400.0), "FINAL", "m",
                                   "FLOAT"))
    return st


def test_background_residual_reconstruction():
    bg = _synthetic_background()
    grid = mp.make_grid("G050", bg=bg)
    bias = lambda la, lo: 3.0 * np.sin(np.radians(la * 10)) + 0.1 * (lo - 80)  # noqa: E731
    st = _stations(bg=bg, bias_fn=bias)
    cp = mp.CovParams(3.0, 300.0, 0.3)
    r = mp.map_slot(grid, st, DAY + 6 * 3600 * ts.NS, cp, bg)
    i = np.argmin(np.abs(grid.lat - 20.0))
    j = np.argmin(np.abs(grid.lon - 80.0))
    assert r.cls[i, j] in (1, 2, 3)
    truth = r.bg[i, j] + bias(20.0, 80.0)
    assert abs(r.pwv[i, j] - truth) < 3 * r.sigma[i, j] + 0.5
    far = (np.argmin(np.abs(grid.lat - 37.5)), np.argmin(np.abs(grid.lon - 68.5)))
    assert r.cls[far] == 4 and np.isnan(r.pwv[far])                  # unsupported -> fill value
    assert np.isfinite(r.bg[far])                                      # background layer still available
    rows = mp.loso(st, DAY + 6 * 3600 * ts.NS, cp, bg)
    m = mp.cv_metrics(rows)
    assert m["rmse"] < 2.0 and 0.5 < m["within_1sigma"] <= 1.0


def test_reml_recovers_range():
    rng = np.random.default_rng(3)
    slots = []
    for k in range(6):
        n = 120
        lat, lon = rng.uniform(8, 32, n), rng.uniform(70, 95, n)
        u = mp.to_unit(lat, lon)
        C = mp.cov_exp(mp.gc_km(u, u), 0, 2.0, 200.0) + np.eye(n) * 0.25
        r = rng.multivariate_normal(np.zeros(n), C) + 1.0
        slots.append((u, r, np.full(n, 0.3)))
    cp = mp.fit_reml(slots)
    assert cp.fitted and 100 < cp.L_km < 400 and 1.3 < cp.sigma_s < 2.8


def test_pwv_map_cli_height_scaling(tmp_path):
    import pwv_map
    st = _stations(n=30)
    res = tmp_path / "res"
    res.mkdir()
    for s in st:
        rows = []
        for i in range(0, 288, 1):
            rows.append({"timestamp_utc": ts.iso(int(s.t_utc[i])) + "Z", "station": s.station, "latitude": s.lat,
                         "longitude": s.lon, "height": s.h_ell, "H_ant": s.h_orth, "PWV": s.pwv[i],
                         "sigma_PWV_total": 0.8, "ZTD": 2400.0, "QC_flag": 0, "convergence_flag": "CONVERGED",
                         "product_tier": "FINAL", "manifest_id": "m", "solution_type": "FLOAT"})
        out.write_csv(str(res / f"{s.station}_2024197_PWV.csv"), rows, out.ZTD_COLUMNS + out.PWV_EXTRA, "synthetic")
    rc = pwv_map.main([str(res), "--day", "2024-197", "--grid", "G050", "--out", str(tmp_path / "maps"),
                       "--slot-min", "180", "--cv", "--hourly", "--quiet"])
    assert rc == 0
    import netCDF4
    files = sorted(x for x in os.listdir(tmp_path / "maps") if x.endswith(".nc"))
    assert any("_01H_2024197_" in x for x in files)
    assert os.path.exists(tmp_path / "maps" / "station_residual_history.csv")
    f = [x for x in files if "_180M_" in x][0]
    assert f.startswith("INPWV_G050_180M_2024197_FINAL")
    with netCDF4.Dataset(str(tmp_path / "maps" / f)) as nc:
        assert nc.background_source == "HEIGHT_SCALING_ONLY"
        assert nc["PWV"].shape == (8, len(nc["lat"]), len(nc["lon"]))
        assert set(np.unique(nc["confidence_class"][:])) <= {0, 1, 2, 3, 4}
        assert (nc["supportable_grid_spacing"][:] == -1).all()            # not determined before § 31


def test_aggregation_variance_aware():
    """0.25 -> 0.5 deg (TDS § 27.1): area-weighted means; sigma of a fully correlated field is not reduced."""
    g = mp.make_grid("G025", domain=(10.0, 14.0, 76.0, 80.0))
    ny, nx = len(g.lat), len(g.lon)
    LA, LO = np.meshgrid(g.lat, g.lon, indexing="ij")
    f = lambda a: np.array(a, dtype=float)  # noqa: E731
    pwv = 20.0 + LA + 0.5 * LO
    cls = np.full((ny, nx), 1, dtype=np.int8)
    cls[0, :] = 4                                                # unsupported first row
    r = mp.SlotResult(f(pwv), np.full((ny, nx), 2.0), f(pwv), np.zeros((ny, nx)), f(pwv) * 6.4, f(pwv) * 6.4,
                      np.full((ny, nx), 2000.0), np.full((ny, nx), 12.0), np.ones((ny, nx)), cls,
                      np.full((ny, nx), 10.0), np.full((ny, nx), 3), np.full((ny, nx), 5), np.full((ny, nx), 20.0),
                      np.full((ny, nx), 950.0), np.full((ny, nx), 280.0), np.array([1]), 1)
    cg, out = mp.aggregate(g, [r], "G050", mp.CovParams(2.0, 1e6, 0.5))   # L huge -> full correlation
    o = out[0]
    i = int(np.argmin(np.abs(cg.lat - 12.0)))
    j = int(np.argmin(np.abs(cg.lon - 78.0)))
    assert o.pwv[i, j] == np.float64(20.0 + 12.0 + 39.0) or abs(o.pwv[i, j] - 71.0) < 1e-9   # linear field
    assert abs(o.sigma[i, j] - 2.0) < 1e-3                       # fully correlated: no sqrt(n) reduction
    cg2, out2 = mp.aggregate(g, [r], "G050", mp.CovParams(2.0, 1.0, 0.5))  # tiny L -> independent cells
    assert out2[0].sigma[i, j] < 1.3
    assert o.cls[0, 2] == 4 and np.isnan(o.pwv[0, 2])            # worst class kept, fill value


def test_station_grey_and_blacklist(tmp_path):
    rng = np.random.default_rng(0)
    st = [mp.StationSeries(f"S{k:03d}", 17 + 0.3 * (k % 5), 78 + 0.3 * (k // 5), 500, 500, np.zeros(1, np.int64),
                           np.zeros(1), np.zeros(1), np.zeros(1), "FINAL", "m", "FLOAT") for k in range(20)]
    path = str(tmp_path / "hist.csv")
    for d in range(1, 31):
        daily = {s.station: (rng.normal(0, 0.3) + (5.0 if s.station == "S007" else 0)
                             + (3.0 if s.station == "S012" else 0), 24) for s in st}
        rows = mp.update_residual_history(path, f"2024{d:03d}", daily, st)
    rows = mp.update_residual_history(path, "2024030", {"S000": (0.0, 24)}, st[:1])   # day replaced, not added
    assert sum(r["day"] == "2024030" for r in rows) == 1
    grey, black, table = mp.station_lists(rows, "2024031")
    assert black == ["S007"] and grey == ["S012"]
    assert abs(table["S007"]["deviation_mm"] - 5.0) < 0.5


def test_hourly_window():
    s = mp.StationSeries("A", 17, 78, 500, 500, DAY + np.arange(0, 3600, 300, dtype=np.int64) * ts.NS,
                         np.arange(12.0), np.ones(12), np.zeros(12), "FINAL", "m", "FLOAT")
    half, n = mp.SLOT_WINDOWS["01H"]
    assert mp.slot_values(s, DAY + 1800 * ts.NS, half, n) == (5.5, 1.0)
    assert mp.slot_values(s, DAY, half, n) is None                       # only 6 of 12 values
