"""Stage-1 product manager scenarios (TDS § 11-12, § 35 Stage 1): naming, cache, download via HTTP,
verification (corrupt file), family fallback without mixing, WAIT_FOR_TIER / FAIL, offline mode."""
import gzip
import http.server
import os
import threading

import pytest

from ppp import products as prd
from ppp import settings
from ppp import timesys as ts
from tests.test_orbclk import truth

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
    L = [f"#dP2024  1 {ts.from_ns(day)[2]:2d}  0  0  0.00000000     289 ORBIT IGS20 FIT AIUB",   # CODE writes AIUB here
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


def test_outputs_archived_not_overwritten(tmp_path):
    import json

    import pwv_ppp
    (tmp_path / "HYDE_2024015_ZTD.csv").write_text("old")
    (tmp_path / "HYDE_2024015_manifest.json").write_text(json.dumps({"manifest_id": "abc"}))
    mid = pwv_ppp._archive_previous_outputs(str(tmp_path), "HYDE_2024015")
    assert mid == "abc"
    assert (tmp_path / "superseded" / "abc" / "HYDE_2024015_ZTD.csv").read_text() == "old"
    assert not (tmp_path / "HYDE_2024015_ZTD.csv").exists()


def test_published_checksum_manifest(tmp_path, server, monkeypatch):
    """TDS § 12.4 item 2: files are compared with the source's SHA512SUMS; a mismatch rejects the file."""
    import hashlib
    root, url = server
    tmpl = url + "/{yyyy}/{filename}"
    monkeypatch.setattr(settings, "PRODUCT_SOURCES", [tmpl])
    monkeypatch.setattr(settings, "CHECKSUM_MANIFESTS", {tmpl: [("sha512", url + "/{yyyy}/SHA512SUMS")]})
    monkeypatch.setattr(settings, "HTTP_RETRIES", 1)
    monkeypatch.setenv("NO_PROXY", "127.0.0.1,localhost")
    monkeypatch.setenv("no_proxy", "127.0.0.1,localhost")
    _publish(root, "COD0OPSFIN", [DAY])
    d = root / "2024"
    lines = []
    for f in sorted(d.iterdir()):
        h = hashlib.sha512(f.read_bytes()).hexdigest()
        if "CLK" in f.name:
            h = "0" * 128                                   # wrong digest for the clock file
        lines.append(f"{h}  {f.name}")
    (d / "SHA512SUMS").write_text("\n".join(lines) + "\n")
    assert prd.parse_checksum_manifest("ab12  *x.gz\nnot a line\n") == {"x.gz": "ab12"}
    dl = prd.Downloader(str(tmp_path / "cache"))
    p, meta = prd.fetch_product(dl, "SP3", "COD0OPSFIN", DAY)
    assert p is not None and meta["verification"]["published_checksum"]["result"] == "OK"
    p2, why = prd.fetch_product(dl, "CLK", "COD0OPSFIN", DAY)
    assert p2 is None
    assert sum(a["url"].endswith("SHA512SUMS") for a in dl.attempts) == 1      # manifest fetched once per run


def test_best_available_and_as_of(tmp_path):
    """TDS § 19.3: FINAL > RAPID_M > RAPID_0 > ..., then version, then run time; as_of reproduces an earlier choice."""
    def man(mid, tier, ver, when, d="", status="OK_FLOAT_ONLY"):
        folder = tmp_path / d if d else tmp_path
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "HYDE_2024015_PWV.csv").write_text("x\n")
        prd.write_json(str(folder / "HYDE_2024015_manifest.json"),
                       {"manifest_id": mid, "station": "HYDE", "day": "2024-015", "tier": tier, "status": status,
                        "software_version": ver, "run_timestamp": when, "campaign": "default",
                        "outputs": ["HYDE_2024015_PWV.csv"]})
    man("a", "RAPID_0", "0.9.0+abc", "2024-01-16T10:00:00Z", "superseded/a")
    man("b", "FINAL", "0.9.0+abc", "2024-02-05T10:00:00Z", "superseded/b")
    man("c", "RAPID_M", "1.0.0", "2024-03-01T10:00:00Z")          # later, lower tier: must not win
    ms = prd.scan_manifests(str(tmp_path))
    assert len(ms) == 3
    assert prd.select_manifests(ms)[("HYDE", "2024-015")]["manifest_id"] == "b"
    assert prd.select_manifests(ms, as_of="2024-01-20")[("HYDE", "2024-015")]["manifest_id"] == "a"
    assert prd.select_manifests(ms, as_of="2024-01-01") == {}
    p = prd.selected_file(prd.select_manifests(ms)[("HYDE", "2024-015")], "_PWV.csv")
    assert p.endswith(os.path.join("superseded", "b", "HYDE_2024015_PWV.csv"))


def test_foreign_agency_rejected(tmp_path):
    p = tmp_path / "x.sp3"
    p.write_text(_sp3_text(DAY).replace("FIT AIUB", "FIT  GFZ"))
    v = prd.verify_file(str(p), "SP3", "COD0OPSFIN", DAY)
    assert not v["ok"] and any("not CODE" in r for r in v["reasons"])


def test_earthdata_login_and_isolation(tmp_path, monkeypatch):
    """Credentials are sent to Earthdata hosts only; a login page triggers EARTHDATA_LOGIN_REQUIRED, not a crash."""
    import base64
    import http.server
    import threading
    want = "Basic " + base64.b64encode(b"alice:s3cret").decode()
    seen = []

    class H(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            seen.append(self.headers.get("Authorization"))
            if self.headers.get("Authorization") == want:
                self.send_response(200)
                self.end_headers()
                self.wfile.write(b"0123abcd  file.gz\n")
            else:
                self.send_response(200)
                self.send_header("Content-Type", "text/html")
                self.end_headers()
                self.wfile.write(b"<!DOCTYPE html><html>Earthdata Login</html>")

        def log_message(self, *a):
            pass
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    url = f"http://127.0.0.1:{srv.server_address[1]}/SHA512SUMS"
    monkeypatch.setenv("NO_PROXY", "127.0.0.1,localhost")
    monkeypatch.setenv("no_proxy", "127.0.0.1,localhost")
    monkeypatch.setattr(settings, "HTTP_RETRIES", 1)
    monkeypatch.delenv("EARTHDATA_USERNAME", raising=False)
    monkeypatch.delenv("EARTHDATA_PASSWORD", raising=False)
    monkeypatch.chdir(tmp_path)
    from ppp import credentials
    monkeypatch.setattr(credentials, "HERE", str(tmp_path))
    monkeypatch.setenv("HOME", str(tmp_path))                          # no ~/.netrc
    # host not an Earthdata host -> never receives credentials
    monkeypatch.setenv("EARTHDATA_USERNAME", "alice")
    monkeypatch.setenv("EARTHDATA_PASSWORD", "s3cret")
    dl = prd.Downloader(str(tmp_path / "c"))
    assert dl.get(url, str(tmp_path / "a")) is False and seen[-1] is None
    # Earthdata host -> credentials used, file downloaded
    monkeypatch.setattr(settings, "EARTHDATA_HOSTS", ["127.0.0.1:%d" % srv.server_address[1]])
    dl = prd.Downloader(str(tmp_path / "c"))
    assert dl.get(url, str(tmp_path / "b")) is True and seen[-1] == want
    # wrong password -> login page -> clear hint, host skipped
    monkeypatch.setenv("EARTHDATA_PASSWORD", "wrong")
    dl = prd.Downloader(str(tmp_path / "c"))
    assert dl.get(url, str(tmp_path / "d")) is False
    assert "Earthdata login" in list(dl.dead_hosts.values())[0]
    srv.shutdown()
    # credentials never enter the configuration snapshot / hash
    monkeypatch.setattr(settings, "EARTHDATA_PASSWORD", "s3cret")
    assert "EARTHDATA_PASSWORD" not in settings.snapshot() and "s3cret" not in str(settings.snapshot())
