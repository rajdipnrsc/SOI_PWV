"""Stage 9/10 helpers: reprocessing triggers and grid-rebuild rule (TDS § 19.3)."""
import importlib.util
import os

from ppp import products as prd

_spec = importlib.util.spec_from_file_location("reprocess", os.path.join(os.path.dirname(__file__), "..", "network",
                                                                         "reprocess.py"))
rp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rp)


def _man(folder, st, tier, nav=True):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f"{st}_2024015_PWV.csv").write_text("x\n")
    m = {"manifest_id": st + tier, "station": st, "day": "2024-015", "tier": tier, "status": "OK_FLOAT_ONLY",
         "software_version": "0.9.0", "run_timestamp": "2024-02-01T00:00:00Z", "outputs": [f"{st}_2024015_PWV.csv"],
         "rinex_files": [{"name": f"{st}015A.24o"}]}
    if nav:
        m["nav_file"] = {"name": f"{st}015A.24n"}
    prd.write_json(str(folder / f"{st}_2024015_manifest.json"), m)
    return m


def test_maps_to_rebuild(tmp_path):
    res = tmp_path / "res"
    _man(res, "AAAA", "FINAL")
    _man(res, "BBBB", "FINAL")
    maps = tmp_path / "maps"
    maps.mkdir()
    (maps / "INPWV_G025_15M_2024015_RAPID_v1.0.nc").write_bytes(b"")
    (maps / "INPWV_G025_15M_2024016_RAPID_v1.0.nc").write_bytes(b"")
    todo = rp.maps_to_rebuild(str(res), str(maps))
    assert [t for t, _ in todo] == ["2024015"]
    _man(res / "superseded" / "x", "CCCC", "RAPID_M")
    _man(res, "CCCC", "RAPID_M")
    assert rp.maps_to_rebuild(str(res), str(maps)) == []          # not all inputs FINAL yet


def test_locate_inputs(tmp_path):
    rin = tmp_path / "rinex" / "2024" / "015"
    rin.mkdir(parents=True)
    for f in ("AAAA015A.24o", "AAAA015A.24n", "BBBB015A.24o", "BBBB015A.24n"):
        (rin / f).write_text("")
    m1 = _man(tmp_path / "r", "AAAA", "RAPID_M")
    m2 = _man(tmp_path / "r2", "BBBB", "RAPID_M", nav=False)        # older manifest without nav_file
    for m in (m1, m2):
        o, n = rp.locate_inputs(m, str(tmp_path / "rinex"))
        assert o.endswith("015A.24o") and n.endswith("015A.24n")
