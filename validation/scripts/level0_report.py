#!/usr/bin/env python3
"""Regenerate validation/level0_report.md: Level-0 component tests (TDS § 30.1) with measured numbers."""
import copy
import datetime
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from ppp import coords as co  # noqa: E402
from ppp import corrections as cr  # noqa: E402
from ppp import estimator as est  # noqa: E402
from ppp import orbclk  # noqa: E402
from ppp import preprocess as pp  # noqa: E402
from ppp import settings  # noqa: E402
from ppp import timesys as ts  # noqa: E402
from ppp import troposphere as tr  # noqa: E402
from tests import synth  # noqa: E402
from tests.test_corrections import _hardisp_cases  # noqa: E402
from tests.test_orbclk import EPH, truth, write_sp3  # noqa: E402
from tests.test_preprocess import SLIPS, run_injection  # noqa: E402

DATA = os.path.join(ROOT, "tests", "data")
rows = []


def add(test, dataset, threshold, measured, ok):
    rows.append((test, dataset, threshold, measured, "PASS" if ok else "FAIL"))


# SP3 interpolation
import tempfile  # noqa: E402
import pathlib  # noqa: E402
with tempfile.TemporaryDirectory() as d:
    p = pathlib.Path(d) / "t.sp3"
    t0 = EPH.toe - 6 * 3600 * ts.NS
    write_sp3(p, t0, 145)
    sp = orbclk.read_sp3(str(p))
    tq = t0 + np.arange(60, 12 * 3600 - 60, 30, dtype=np.int64)[5:-5] * ts.NS
    pos, _, ok, _ = orbclk.sp3_interp(sp, "G02", tq)
    err = np.linalg.norm(pos - np.array([truth(t) for t in tq]), axis=1) * 1e3
add("SP3 Lagrange interpolation (11 nodes)", "analytic GPS orbit, 5-min nodes, 1420 epochs", "RMS < 1 mm, max < 5 mm",
    f"RMS {np.sqrt(np.mean(err ** 2)):.4f} mm, max {err.max():.4f} mm", np.sqrt(np.mean(err ** 2)) < 1 and err.max() < 5)
# solid tide
r = np.loadtxt(os.path.join(DATA, "iers_dehanttideinel_ref.txt"), comments="#")
e = max(np.max(np.abs(cr.solid_tide(x[0:3], x[None, 3:6], x[None, 6:9], np.array([int(x[9])]))[0] - x[10:13])) for x in r)
add("Solid Earth tide", "official IERS DEHANTTIDEINEL, 10 stations x 100 epochs", "< 0.1 mm", f"max {e * 1e3:.2e} mm",
    e < 1e-4)
# ocean loading
worst = 0.0
for c in _hardisp_cases():
    blq = np.array(c["blq"])
    tt = ts.to_ns(*c["t"][:3], c["t"][3], c["t"][4], c["t"][5]) + np.arange(c["n"]) * c["step"] * ts.NS
    worst = max(worst, np.max(np.abs(cr.ocean_loading(blq[:3], blq[3:], tt) - np.array(c["out"]))))
add("Ocean tide loading (HARDISP)", "official IERS HARDISP, 6 cases x 48 epochs", "< 0.1 mm",
    f"max {worst * 1e3:.4f} mm (reference printed to 1 um)", worst < 1e-4)
# wind-up
w = np.loadtxt(os.path.join(DATA, "rtklib_windup_ref.txt"), comments="#")
starts = list(np.nonzero(w[:, 8] == 1)[0]) + [len(w)]
wr = 0.0
for a, b in zip(starts[:-1], starts[1:]):
    x = w[a:b]
    ex, ey, ez = cr.nominal_attitude(x[:, 2:5], x[:, 10:13])
    la, lo, _ = co.ecef_to_geodetic(x[0, 5:8])
    dd = cr.windup(ex, ey, x[:, 2:5], x[:, 5:8], la, lo) - x[:, 9]
    wr = max(wr, np.max(np.abs(dd - np.round(np.median(dd)))))
add("Phase wind-up", "RTKLIB 2.4.3 windupcorr, 7 arcs / 1060 epochs", "< 0.001 cycle", f"max {wr:.1e} cycle", wr < 1e-3)
# Sagnac: model formulation (ECEF position at t_tx rotated by w*tau) vs an independent inertial-frame light-time
# solution (satellite propagated in inertial space: r_I(t) = R3(-w (t - t_rx)) r_ECEF(t)); plus RTKLIB info
rr = np.array([1208298.2, 5967329.1, 1895460.8])
dmax = drtk = 0.0
W = settings.OMEGA_E
for k in range(0, 43200, 600):
    t_rx = EPH.toe - 6 * 3600 * ts.NS + k * ts.NS
    rs_e = truth(t_rx)
    if np.dot(rs_e - rr, rr) < 0:
        continue
    # model: iterate tau with ECEF positions and rotation
    tau = 0.07
    for _ in range(6):
        rs = truth(t_rx - int(tau * 1e9))
        rho = np.linalg.norm(cr.sagnac_rotate(rs[None], tau)[0] - rr)
        tau = rho / settings.C_LIGHT
    # independent: inertial frame = ECEF at t_rx; solve |r_I(t_tx) - rr| = c (t_rx - t_tx)
    tt = 0.07
    for _ in range(8):
        th = -W * tt
        R = np.array([[np.cos(th), np.sin(th), 0], [-np.sin(th), np.cos(th), 0], [0, 0, 1.0]])
        rI = R.T @ truth(t_rx - int(tt * 1e9))
        tt = np.linalg.norm(rI - rr) / settings.C_LIGHT
    rho_i = np.linalg.norm(rI - rr)
    dmax = max(dmax, abs(rho - rho_i))
    rs_tx = truth(t_rx - int(tau * 1e9))
    rho_rtk = np.linalg.norm(rs_tx - rr) + W * (rs_tx[0] * rr[1] - rs_tx[1] * rr[0]) / settings.C_LIGHT
    drtk = max(drtk, abs(rho - rho_rtk))
add("Sagnac / light time (rotation of r(t_tx) by w*tau)", "independent inertial-frame light-time solution, 12 h pass",
    "< 0.1 mm", f"max {dmax * 1e3:.2e} mm (RTKLIB first-order geodist differs by {drtk * 1e3:.3f} mm: its truncation)",
    dmax < 1e-4)
# GMF / VMF3 / Sun-Moon
mh, mw = tr.gmf(55055.0, 0.6708665767, -1.393397187, 844.715, np.pi / 2 - 1.278564131)
g = max(abs(mh - 3.425245519339138678), abs(mw - 3.449589116182419257))
add("GMF mapping function", "IERS GMF.F test vector", "exact (1e-12)", f"{g:.1e}", g < 1e-12)
lam, beta, dist = cr.moon_ecliptic(np.array([2448724.5]))
add("Moon position (Meeus 47)", "Meeus example 47.a", "< 1e-5 deg", f"dlon {abs(np.degrees(lam[0]) - 133.162655):.1e} deg",
    abs(np.degrees(lam[0]) - 133.162655) < 1e-5)
# frame transformation equivalence
X = np.array([1208298.991, 5967329.098, 1895460.127])
V08 = co.velocity_itrf2008_from_pmm(X)
a = co.helmert_inverse(X + V08 * 19.5, ("ITRF2020", "ITRF2008"), 2024.5)
b = co.helmert_inverse(X, ("ITRF2020", "ITRF2008"), 2005.0) + (V08 - co.velocity_rate_terms(X, ("ITRF2020", "ITRF2008"))) * 19.5
add("Frame transformation (propagation-order equivalence)", "HYDE SOI coordinate, 2005.0 -> 2024.5", "< 0.1 mm",
    f"{np.max(np.abs(a - b)) * 1e3:.1e} mm", np.max(np.abs(a - b)) < 1e-4)
# estimator synthetic truth
obs, arcs, truth_ = synth.make_obs()
sol = est.solve(obs, copy.deepcopy(arcs), est.EstConfig(mode="fixed"), max_passes=0)
z, s, _ = est.state_series(sol, "ZWD")
okk = np.isfinite(z)
zs = (z[okk] - truth_["zwd"][okk]) / s[okk]
add("Kalman/RTS numerics (synthetic truth)", "24 h, 12 satellites, ZWD RW 6 mm/sqrt(h)", "z-STD 0.9-1.1",
    f"z-STD {np.std(zs):.3f}, RMSE {np.sqrt(np.mean((z[okk] - truth_['zwd'][okk]) ** 2)) * 1e3:.2f} mm, "
    f"NIS/dof {sol.nis:.3f}", 0.9 <= np.std(zs) <= 1.1)
# slip injection
tab = {k: [0, 0] for k in SLIPS}
for sd in (3, 4, 6):
    for k, (h, n) in run_injection(seed=sd).items():
        tab[k][0] += h
        tab[k][1] += n
fa = 0
days = 0.0
for sd in (5, 8, 9):
    sel, t, el, vis = synth.make_raw(seed=sd)
    ai = pp.detect_slips_and_arcs(sel, t, el, vis, np.radians(7))
    fa += ai.reasons["MW"] + ai.reasons["GF"]
    days += vis.sum() * 30 / 86400
slip_txt = ", ".join(f"({a},{b}): {h}/{n}" for (a, b), (h, n) in tab.items())

# full modelled range vs independent implementation
from tests.test_model_decomposition import run_decomposition  # noqa: E402
with tempfile.TemporaryDirectory() as d:
    mx, nobs = run_decomposition(d)
add("Full modelled range, term by term", f"HYDE 2024-015 geometry, broadcast-node products, {nobs} sat-epochs; "
    "independent test-only implementation (`tests/tools/independent_range.py`)", "< 1 mm per term, total < 2 mm",
    ", ".join(f"{k} {v * 1e3:.3f} mm" for k, v in mx.items()),
    all(v < (2e-3 if k == "total" else 1e-3) for k, v in mx.items()))
git = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
L = ["# Level-0 component tests (TDS § 30.1) — measured results", "",
     f"Generated {datetime.date.today().isoformat()} by `validation/scripts/level0_report.py` at commit `{git}`. "
     "All numbers are measured by this script; references are official IERS routine outputs (gfortran builds), "
     "RTKLIB 2.4.3 (compiled test-only), published test vectors, or analytic truth.", "",
     "| Test | Dataset | Threshold | Measured | Result |", "|---|---|---|---|---|"]
L += [f"| {a} | {b} | {c} | {d} | {e} |" for a, b, c, d, e in rows]
L += ["", "## Injected cycle slips (synthetic, § 30.4 method; real-data run still required)", "",
      f"Detected/injected per (dN1, dN2) over 3 datasets: {slip_txt}.", "",
      f"False alarms on clean data: {fa} in {days:.1f} satellite-days = {fa / days:.2f} per satellite-day "
      "(gate ≤ 1).", "",
      "(1,1)-type slips change neither the Melbourne–Wübbena combination nor the geometry-free combination by more "
      "than ~5 cm and are not reliably detectable by any MW/GF detector at 30 s; the § 30.4 gate refers to "
      "MW-detectable slips.", "",
      "## Not yet covered at Level 0", "",
      "* The modelled-range decomposition above uses broadcast-node products (CODE products were not reachable); "
      "repeat it with one CODE station-day (same script) when products are available.",
      "* Satellite PCO/PCV against RTKLIB on real igs20 entries (the convention is verified with a synthetic "
      "asymmetric antenna, `tests/test_antenna.py`).",
      "* Pole tide: checked against the IERS formula re-implemented independently in the test (no official "
      "routine output available)."]
with open(os.path.join(ROOT, "validation", "level0_report.md"), "w") as fh:
    fh.write("\n".join(L) + "\n")
print("\n".join(L))
