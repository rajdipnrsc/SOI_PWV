"""§ 24 / § 31 experiment tools on a synthetic residual field with known covariance."""
import numpy as np

from ppp import mapping as mp
from ppp import spatial_validation as sv


def _field(n=150, L=150.0, seed=2):
    rng = np.random.default_rng(seed)
    lat, lon = rng.uniform(12, 24, n), rng.uniform(72, 86, n)
    u = mp.to_unit(lat, lon)
    C = mp.cov_exp(mp.gc_km(u, u), 0, 3.0, L) + np.eye(n) * 0.3 ** 2
    r = rng.multivariate_normal(np.zeros(n), C)
    bg = np.full(n, 30.0)
    return sv.SlotData([f"S{i:03d}" for i in range(n)], lat, lon, np.zeros(n), bg + r, np.full(n, 0.5), bg)


def test_cv_designs_and_error_distance():
    sd = _field()
    cp = mp.CovParams(3.0, 150.0, 0.3)
    loso = sv.predict_withheld(sd, np.ones(len(sd.names), bool), np.zeros(len(sd.names), bool), cp)
    assert loso == []
    rows = []
    for i in range(len(sd.names)):
        test = np.arange(len(sd.names)) == i
        rows += sv.predict_withheld(sd, ~test, test, cp)
    m = sv.metrics(rows)
    assert m["skill"] > 0.3 and 0.55 < m["within_1sigma"] < 0.85          # calibrated sigma, skill over bg
    b100 = sv.metrics(sv.block_cv(sd, cp, 100))
    b200 = sv.metrics(sv.block_cv(sd, cp, 200))
    assert m["rmse"] < b100["rmse"] < b200["rmse"] * 1.1                   # larger blocks -> larger errors
    thin_rows = []
    for spc in (30, 75, 150):
        thin_rows += sv.thinning_cv(sd, cp, spc, n_real=3, seed=spc)
    fit = sv.fit_error_distance(rows + thin_rows)
    assert fit["fitted"] and sv.E_of(fit, 150) > sv.E_of(fit, 20)          # error grows with distance


def test_thinning_modes():
    sd = _field(n=200)
    rng = np.random.default_rng(1)
    k = sv.thin(sd.lat, sd.lon, 100, rng)
    u = mp.to_unit(sd.lat[k], sd.lon[k])
    d = mp.gc_km(u, u)
    np.fill_diagonal(d, np.inf)
    assert d.min() >= 100
    ks = sv.thin(sd.lat, sd.lon, 100, rng, "stratified")
    assert len(np.unique(sv.block_ids(sd.lat[ks], sd.lon[ks], 100))) == ks.sum()


def test_supportable_spacing_rules():
    fit = {"E0": 0.5, "a": 0.05, "alpha": 1.0, "fitted": True}           # E(d) = sqrt(0.25 + (0.05 d)^2)
    best, table = sv.supportable_spacing(25.0, 150.0, fit, 0.5, lat_c=20.0)
    # 0.125 deg ~ 13.5 km: density 25 <= 27, accuracy E(6.7) = 0.6 mm -> supported; 0.1 deg: 25 <= 21.6 fails
    assert best == 0.125 and table[0.1]["density"] is False
    best2, _ = sv.supportable_spacing(25.0, 150.0, fit, 0.1)
    assert best2 is None                                                   # no skill -> background only


def test_experiment_script_end_to_end(tmp_path, monkeypatch):
    import importlib.util
    import os
    from ppp import output as out
    from ppp import settings
    from ppp import timesys as ts
    from tests.test_mapping import _stations
    monkeypatch.setattr(settings, "THINNING_SPACINGS_KM", [50, 150])
    res = tmp_path / "res"
    res.mkdir()
    for s in _stations(n=30):
        rows = [{"timestamp_utc": ts.iso(int(s.t_utc[i])) + "Z", "station": s.station, "latitude": s.lat,
                 "longitude": s.lon, "height": s.h_ell, "H_ant": s.h_orth, "PWV": s.pwv[i], "sigma_PWV_total": 0.8,
                 "ZTD": 2400.0, "QC_flag": 0, "convergence_flag": "CONVERGED", "product_tier": "FINAL",
                 "manifest_id": "m", "solution_type": "FLOAT"} for i in range(288)]
        out.write_csv(str(res / f"{s.station}_2024197_PWV.csv"), rows, out.ZTD_COLUMNS + out.PWV_EXTRA, "synthetic")
    spec = importlib.util.spec_from_file_location("spx", os.path.join(os.path.dirname(__file__), "..", "validation",
                                                                      "scripts", "spatial_experiment.py"))
    spx = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(spx)
    rc = spx.run([str(res), "--days", "2024-197", "--every-h", "12", "--realisations", "1", "--out",
                  str(tmp_path / "exp"), "--write-support", "--era5-dir", str(tmp_path / "none")])
    assert rc == 0
    txt = (tmp_path / "exp" / "report.md").read_text()
    assert "LOSO" in txt and "BLOCK_100KM" in txt and "HEIGHT_SCALING_ONLY" in txt
