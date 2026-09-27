"""Station product writers: CSV / NetCDF round trip (TDS § 16-17)."""
import numpy as np

from ppp import estimator as est
from ppp import output as out
from ppp import timesys as ts


def _fm(n=12):
    t = ts.to_ns(2024, 1, 15) + np.arange(n) * 300 * ts.NS
    return est.FiveMin(t, t + 18 * ts.NS, np.full(n, 0.15), np.full(n, 0.004), np.zeros(n), np.zeros(n),
                       np.full(n, 2.2), np.zeros(n), np.full(n, 8, dtype=np.int16), np.zeros(n, dtype=np.int16),
                       np.zeros(n, dtype=np.uint32), ["CONVERGED"] * n)


META = {"station": "TEST", "lat": 17.4, "lon": 78.5, "h": 420.0, "product_AC": "COD", "product_tier": "FINAL",
        "product_files": "[]", "coordinate_mode": "fixed", "coordinate_source": "soi_transformed",
        "mapping_function": "VMF3", "ZWD_process_noise": 6.0, "software_version": "x", "solution_type": "FLOAT",
        "zhd_source": "VMF3_GRID", "manifest_id": "id", "config_hash": "h", "elevation_cutoff": 7.0, "frame": "IGS20"}


def test_csv_and_netcdf(tmp_path):
    import netCDF4
    fm = _fm()
    rows = out.ztd_rows(fm, META)
    assert rows[0]["ZTD"] == 2350.0 and rows[0]["timestamp_utc"] == "2024-01-15T00:00:00Z"
    p = tmp_path / "a.csv"
    out.write_csv(str(p), rows, out.ZTD_COLUMNS, "test")
    lines = p.read_text().splitlines()
    assert lines[0].startswith("#") and lines[2].split(",") == out.ZTD_COLUMNS and len(lines) == 15
    nc = tmp_path / "a.nc"
    out.write_netcdf(str(nc), rows, out.ZTD_COLUMNS, {"station": "TEST"}, 17.4, 78.5, 420.0)
    with netCDF4.Dataset(str(nc)) as d:
        assert d.Conventions.startswith("CF-1.10") and d["ZTD"].units == "mm"
        assert np.allclose(d["ZTD"][:], 2350.0) and d["time"][0] == (ts.to_ns(2024, 1, 15) - ts.to_ns(1970, 1, 1)) // ts.NS
        assert "NO_DATA" in d["QC_flag"].flag_meanings
    tro = tmp_path / "a.tro"
    out.write_sinex_tro(str(tro), "TEST", fm, META, (1.0, 2.0, 3.0))
    assert "+TROP/SOLUTION" in tro.read_text()
