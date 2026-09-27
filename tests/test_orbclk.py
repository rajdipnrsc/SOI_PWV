"""Level-0: SP3 interpolation, clock handling, ERP, Bias-SINEX, ORBEX (TDS § 7.2, § 30.1)."""
import numpy as np
import pytest

from ppp import orbclk
from ppp import timesys as ts
from ppp.rinex import GpsEph

EPH = GpsEph(prn="G02", toc=0, af0=-5.0e-4, af1=5.6e-12, af2=0.0, iode=91, crs=116.125, delta_n=4.2652e-09,
             m0=-2.98868, cuc=6.0331e-06, e=0.0163, cus=8.6147e-06, sqrt_a=5153.9024, toe_sow=86400.0,
             cic=1.7695e-07, omega0=-1.9639, cis=2.4773e-07, i0=0.96748, crc=212.72, omega=-1.26324,
             omega_dot=-7.7675e-09, idot=2.1537e-10, week=2297, sv_accuracy=2, health=0, tgd=-1.77e-8, iodc=91)
EPH.toc = EPH.toe = 2297 * 7 * ts.DAY_NS + 86400 * ts.NS


def truth(t):
    return orbclk.broadcast_pos_clk(EPH, int(t))[0]


def write_sp3(path, t0, n, step_s=300, sats=("G02",)):
    y0, m0, d0, h0, mi0, s0 = ts.from_ns(t0)
    L = [f"#dP{y0:4d} {m0:2d} {d0:2d} {h0:2d} {mi0:2d} {s0:11.8f} {n:7d} ORBIT IGS20 FIT  COD",
         "## 2297      0.00000000   300.00000000 60324 0.0000000000000",
         "+    1   G02", "%c G  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc"]
    for i in range(n):
        t = t0 + i * step_s * ts.NS
        y, m, d, hh, mi, s = ts.from_ns(t)
        L.append(f"*  {y:4d} {m:2d} {d:2d} {hh:2d} {mi:2d} {s:11.8f}")
        for sat in sats:
            r = truth(t) / 1e3
            L.append(f"P{sat}{r[0]:14.6f}{r[1]:14.6f}{r[2]:14.6f}{-500.123456:14.6f}")
    L.append("EOF")
    path.write_text("\n".join(L) + "\n")


def test_sp3_interpolation_accuracy(tmp_path):
    """Interpolated 5-min SP3 vs analytic orbit at 30-s epochs: RMS < 1 mm, max < 5 mm (TDS § 7.2)."""
    t0 = EPH.toe - 6 * 3600 * ts.NS
    p = tmp_path / "t.sp3"
    write_sp3(p, t0, 12 * 12 + 1)
    sp = orbclk.read_sp3(str(p))
    assert sp.frame == "IGS20" and sp.agency == "COD" and len(sp.epochs) == 145
    tq = t0 + np.arange(60, 12 * 3600 - 60, 30, dtype=np.int64)[5:-5] * ts.NS
    pos, vel, ok, edge = orbclk.sp3_interp(sp, "G02", tq)
    err = np.linalg.norm(pos - np.array([truth(t) for t in tq]), axis=1)
    assert ok.all()
    assert np.sqrt(np.mean(err ** 2)) < 1e-3 and err.max() < 5e-3
    # velocity vs numerical derivative of the truth
    v_true = (truth(tq[100] + ts.NS // 2) - truth(tq[100] - ts.NS // 2))
    assert np.linalg.norm(vel[100] - v_true) < 1e-3
    # edge flag for asymmetric node sets, no extrapolation
    _, _, ok2, edge2 = orbclk.sp3_interp(sp, "G02", np.array([t0 + 60 * ts.NS, t0 - 60 * ts.NS]))
    assert ok2[0] and edge2[0] and not ok2[1]


def test_sp3_gap_rejected(tmp_path):
    t0 = EPH.toe - 3 * 3600 * ts.NS
    p = tmp_path / "t.sp3"
    write_sp3(p, t0, 73)
    sp = orbclk.read_sp3(str(p))
    sp.pos[30] = np.nan                                  # bad record (0.000000 in file)
    _, _, ok, _ = orbclk.sp3_interp(sp, "G02", np.array([sp.epochs[30] + 10 * ts.NS]))
    assert not ok[0]


def test_clock_exact_and_interp():
    t = np.arange(0, 3600, 30, dtype=np.int64) * ts.NS
    v = 1e-4 + 1e-9 * np.arange(len(t))
    ck = orbclk.ClockProduct({"G01": (t, v)}, interval=30.0)
    val, ok, ip = orbclk.clk_interp(ck, "G01", np.array([t[5], t[5] + 15 * ts.NS, t[-1] + 30 * ts.NS]))
    assert val[0] == v[5] and not ip[0]
    assert val[1] == pytest.approx((v[5] + v[6]) / 2) and ip[1]
    assert not ok[2]


def test_read_clk_erp_bias_obx(tmp_path):
    clk = tmp_path / "c.clk"
    clk.write_text("""     3.04           C                   G                   RINEX VERSION / TYPE
COD  CODE: CENTER FOR ORBIT DETERMINATION IN EUROPE          ANALYSIS CENTER
G    BERNESE GNSS SOFTWARE 5.4 IGS20_2290                    SYS / PCVS APPLIED
   300    IGS20                                             # OF SOLN STA / TRF
                                                            END OF HEADER
AS G01       2024 01 15 00 00  0.000000  1    1.234567890123E-04
AS G01       2024 01 15 00 00 30.000000  1    1.234567990123E-04
AS R01       2024 01 15 00 00  0.000000  1    1.0E-04
""")
    c = orbclk.read_clk(str(clk))
    assert c.agency == "COD" and c.trf == "IGS20" and "IGS20_2290" in c.pcvs
    assert list(c.sats) == ["G01"] and c.sats["G01"][1][0] == pytest.approx(1.234567890123e-4)
    erp = tmp_path / "e.erp"
    erp.write_text("""version 2
  MJD      Xpole   Ypole  UT1-UTC    LOD  Xsig  Ysig   UTsig LODsig  Nr Nf Nt     Xrt    Yrt  Xrtsig Yrtsig
             10**-6"        .1us    .1us/d    10**-6"     .1us  .1us/d                10**-6"/d    10**-6"/d
60324.50   40000  300000  -120000    0     10    10     5     5    100  0  0       0      0     0     0
60325.50   42000  302000  -121000    0     10    10     5     5    100  0  0       0      0     0     0
""")
    e = orbclk.read_erp(str(erp))
    xp, yp, du = orbclk.erp_interp(e, 60325.0)
    assert xp == pytest.approx(0.041) and yp == pytest.approx(0.301) and du == pytest.approx(-0.01205)
    erp.write_text(erp.read_text().replace("-120000", "-370120000").replace("-121000", "-370121000"))
    assert orbclk.erp_interp(orbclk.read_erp(str(erp)), 60325.0)[2] == pytest.approx(-0.01205)   # UT1-TAI input
    bia = tmp_path / "b.bia"
    bia.write_text("""%=BIA 1.00 COD 2024:016:00000 COD 2024:015:00000 2024:016:00000 A 00000
+BIAS/SOLUTION
*BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___
 OSB  G063 G01           C1C       2024:015:00000 2024:016:00000 ns                  -1.0000      0.0100
 OSB  G063 G01           L1C       2024:015:00000 2024:016:00000 ns                   0.5000      0.0100
-BIAS/SOLUTION
""")
    b = orbclk.read_bias(str(bia))
    t = ts.to_ns(2024, 1, 15, 12)
    assert b.get("G01", "C1C", t) == pytest.approx(-1e-9 * 299792458.0)
    assert b.signals("G01") == {"C1C", "L1C"}
    obx = tmp_path / "a.obx"
    obx.write_text("""%=ORBEX  0.09
+EPHEMERIS/DATA
## 2024 01 15 00 00  0.000000000   1
 ATT G01                  4  1.000000000000  0.000000000000  0.000000000000  0.000000000000
-EPHEMERIS/DATA
""")
    a = orbclk.read_obx(str(obx))
    ex, ey, ez = orbclk.att_axes(a, "G01", ts.to_ns(2024, 1, 15))
    assert np.allclose(ez, [0, 0, 1])


def test_quaternion_rotation():
    # 90 deg about z: v' = R v maps x -> y
    q = np.array([np.cos(np.pi / 4), 0, 0, np.sin(np.pi / 4)])
    assert np.allclose(orbclk.quat_to_matrix(q) @ [1, 0, 0], [0, 1, 0])


def test_broadcast_orbit_radius():
    r, dt = orbclk.broadcast_pos_clk(EPH, EPH.toe)
    assert 26.0e6 < np.linalg.norm(r) < 27.0e6
    r2, dt2 = orbclk.broadcast_pos_clk(EPH, EPH.toe, include_rel=False)
    assert abs(dt - dt2) < 5e-8
