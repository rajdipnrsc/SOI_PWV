"""Stage-1 product manager scenarios (TDS § 11-12, § 35 Stage 1): naming, cache, download via HTTP,
verification (corrupt file), family fallback without mixing, WAIT_FOR_TIER / FAIL, offline mode."""
import gzip
import http.server
import os
import threading

import numpy as np
import pytest

from ppp import products as prd
from ppp import settings
from ppp import timesys as ts
from tests.test_orbclk import EPH, truth

DAY = ts.to_ns(2024, 1, 15)


def test_candidate_names():
    n = prd.candidate_names("SP3", "COD0OPSFIN", DAY)
    assert n[0] == "COD0OPSFIN_20240150000_01D_05M_ORB.SP3.gz"
    assert prd.candidate_names("CLK", "CODMOPSRAP", DAY)[2] == "CODMOPSRAP_20240150000_01D_30S_CLK.CLK"
    legacy = prd.candidate_names("SP3", "COD0OPSFIN", ts.to_ns(2020, 1, 1))
    w, d = ts.gps_week_dow(ts.to_ns(2020, 1, 1))
    assert legacy[0] == f"COD{w}{d}.EPH.Z"
    assert prd.needs_osb(["C1C", "C2W"]) and not prd.needs_osb(["C1W", "C2W"])


def _sp3_text(day, sats=24):
    L = [f"#dP2024  1 {ts.from_ns(day)[2]:2d}  0  0  0.00000000     289 ORBIT IGS20 FIT  COD",
         "## 2297      0.00000000   300.00000000 60324 0.0000000000000", "%c G  cc GPS"]
    for i in range(289):
        t = day + i * 300 * ts.NS
        y, m, d, hh, mi, s = ts.from_ns(t)
        L.append(f"*  {y:4d} {m:2d} {d:2d} {hh:2d} {mi:2d} {s:11.8f}")
        r = truth(t) / 1e3
        for k in range(sats):
            L.append(f"PG{k + 1:02d}{r[0] + k:14.6f}{r[1]:14.6f}{r[2]:14.6f}{-1.0:14.6f}")
    L.append("EOF")
    return "\n".join(L) + "\n"


def _clk_text(day, sats=24):
    h = lambda a, b: f"{a:<60s}{b}"  # noqa: E731
    L = [h("     3.04           C                   G", "RINEX VERSION / TYPE"),
         h("COD  CODE", "ANALYSIS CENTER"),
         h("G BERNESE GNSS SOFTWARE 5.4 IGS20_2290", "SYS / PCVS APPLIED"),
         h("   300    IGS20", "# OF SOLN STA / TRF"),
         h("", "END OF HEADER")]
    for i in range(0, 2880, 60):
        y, m, d, hh, mi, s = ts.from_ns(day + i * 30 * ts.NS)
        for k in range(sats):
            L.append(f"AS G{k + 1:02d}       {y} {m:02d} {d:02d} {hh:02d} {mi:02d} {s:9.6f}  1    1.0E-04")
    y, m, d, hh, mi, s = ts.from_ns(day + 2879 * 30 * ts.NS)
    for k in range(sats):
        L.append(f"AS G{k + 1:02d}       {y} {m:02d} {d:02d} {hh:02d} {mi:02d} {s:9.6f}  1    1.0E-04")
    return "\n".join(L) + "\n"


def _erp_text():
    return "version 2\n  MJD Xpole\n 60323.50 1 2 3 0\n 60324.50 1 2 3 0\n 60325.50 1 2 3 0\n 60326.50 1 2 3 0\n"


@pytest.fixture()
def server(tmp_path):
    root = tmp_path / "srv"
    root.mkdir()

    class H(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *a, **k):
            super().__init__(*a, directory=str(root), **k)

        def log_message(self, *a):
            pass
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), H)
    th = threading.Thread(target=srv.serve_forever, daemon=True)
    th.start()
    yield root, f"http://127.0.0.1:{srv.server_address[1]}"
    srv.shutdown()


def _publish(root, family, days, products=("SP3", "CLK", "ERP"), corrupt=()):
    d = root / "2024"
    d.mkdir(exist_ok=True)
    for day in days:
        for p in products:
            name = prd.candidate_names(p, family, day)[0]
            txt = {"SP3": lambda: _sp3_text(day), "CLK": lambda: _clk_text(day), "ERP": _erp_text}[p]()
            data = gzip.compress(txt.encode())
            if (p, day) in corrupt:
                data = data[:100]
            (d / name).write_bytes(data)


@pytest.fixture()
def local_sources(server, monkeypatch):
    root, url = server
    monkeypatch.setattr(settings, "PRODUCT_SOURCES", [url + "/{yyyy}/{filename}"])
    monkeypatch.setattr(settings, "HTTP_RETRIES", 2)
    monkeypatch.setenv("NO_PROXY", "127.0.0.1,localhost")
    monkeypatch.setenv("no_proxy", "127.0.0.1,localhost")
    return root


def test_download_verify_cache_and_rapid_fallback(tmp_path, local_sources):
    root = local_sources
    days = [DAY - ts.DAY_NS, DAY, DAY + ts.DAY_NS]
    _publish(root, "COD0OPSFIN", [DAY - ts.DAY_NS, DAY])                  # FIN incomplete (no D+1 orbit)
    _publish(root, "CODMOPSRAP", days)                                     # final rapid complete
    dl = prd.Downloader(str(tmp_path / "cache"))
    ps = prd.resolve_products(dl, DAY, DAY + ts.DAY_NS - 30 * ts.NS, code_obs_used=("C1W", "C2W"))
    assert ps.status == "OK_FLOAT_ONLY" and "AR_NOT_REQUESTED" in ps.reason
    assert ps.family == "CODMOPSRAP" and ps.tier == "RAPID_M" and "RAPID_TIER" in ps.flags
    assert {f["filename"][:10] for f in ps.files} == {"CODMOPSRAP"}      # never mixed
    assert ps.frame_label == "IGS20" and ps.antex_name == "igs20_2290"
    cached = [os.path.join(dp, f) for dp, _, fs in os.walk(tmp_path / "cache") for f in fs]
    assert any(f.endswith(".meta.json") for f in cached)
    assert not any(f.endswith(".part") for f in cached)                   # atomic writes
    # offline mode uses the cache only and gives the same result
    dl2 = prd.Downloader(str(tmp_path / "cache"), offline=True)
    ps2 = prd.resolve_products(dl2, DAY, DAY + ts.DAY_NS - 30 * ts.NS, code_obs_used=("C1W", "C2W"))
    assert ps2.family == "CODMOPSRAP" and len(dl2.attempts) == 0


def test_orbit_edge_when_neighbour_missing_everywhere(tmp_path, local_sources):
    _publish(local_sources, "COD0OPSFIN", [DAY - ts.DAY_NS, DAY])
    dl = prd.Downloader(str(tmp_path / "cache"))
    ps = prd.resolve_products(dl, DAY, DAY + ts.DAY_NS - 30 * ts.NS, code_obs_used=("C1W", "C2W"))
    assert ps.family == "COD0OPSFIN" and "ORBIT_EDGE" in ps.flags


def test_corrupt_file_rejected(tmp_path, local_sources):
    _publish(local_sources, "COD0OPSFIN", [DAY - ts.DAY_NS, DAY, DAY + ts.DAY_NS], corrupt=[("CLK", DAY)])
    dl = prd.Downloader(str(tmp_path / "cache"))
    ps = prd.resolve_products(dl, DAY, DAY + ts.DAY_NS - 30 * ts.NS, code_obs_used=("C1W", "C2W"),
                              now_ns=DAY + 400 * ts.DAY_NS)
    assert ps.status == "FAIL" and any("CLK" in m for m in ps.missing)


def test_wait_for_tier_and_fail_without_products(tmp_path):
    dl = prd.Downloader(str(tmp_path / "cache"), offline=True)
    young = prd.resolve_products(dl, DAY, DAY + ts.DAY_NS - 30 * ts.NS, now_ns=DAY + 2 * ts.DAY_NS)
    assert young.status == "WAIT_FOR_TIER"
    old = prd.resolve_products(dl, DAY, DAY + ts.DAY_NS - 30 * ts.NS, now_ns=DAY + 60 * ts.DAY_NS)
    assert old.status == "FAIL" and "NO_CODE_PRODUCTS" in old.reason


def test_ultra_rapid_never_used():
    assert "COD0OPSULT" not in settings.PRODUCT_FAMILIES
    assert "COD0OPSULT" in settings.FORBIDDEN_FAMILIES


def test_manifest_supersession(tmp_path):
    m_old = prd.new_manifest("HYDE", "2024-015")
    m_old["tier"] = "RAPID_M"
    prd.write_json(str(tmp_path / "HYDE_2024015_manifest.json"), m_old)
    old = prd.find_superseded(str(tmp_path), "HYDE", "2024-015", "FINAL", "new-id")
    assert old == [m_old["manifest_id"]]
    h, snap = prd.config_hash()
    assert len(h) == 64 and snap["PRODUCT_FAMILIES"] == settings.PRODUCT_FAMILIES
