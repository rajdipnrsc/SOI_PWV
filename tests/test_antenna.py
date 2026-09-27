"""Level-0: ANTEX reading, PCO/PCV, IF combination, radome fallback (TDS § 7.6-7.7, § 8)."""
import numpy as np
import pytest

from ppp import antenna as an
from ppp import timesys as ts

LABELS = ("START OF FREQUENCY", "END OF FREQUENCY", "NORTH / EAST / UP", "ANTEX VERSION / SYST",
          "PCV TYPE / REFANT", "END OF HEADER", "START OF ANTENNA", "TYPE / SERIAL NO", "DAZI",
          "ZEN1 / ZEN2 / DZEN", "VALID FROM", "END OF ANTENNA")


def _fix(line):
    for lab in LABELS:
        if line.rstrip().endswith(lab) and not line.startswith("   NOAZI"):
            body = line.rstrip()[: -len(lab)].rstrip()
            return f"{body:<60s}{lab}"
    return line


def _freq_block(code, neu, noazi, azi_rows):
    L = [f"   {code}                                                    START OF FREQUENCY",
         f"{neu[0]:10.2f}{neu[1]:10.2f}{neu[2]:10.2f}                              NORTH / EAST / UP",
         "   NOAZI" + "".join(f"{v:8.2f}" for v in noazi)]
    for az, row in azi_rows:
        L.append(f"{az:8.1f}" + "".join(f"{v:8.2f}" for v in row))
    L.append(f"   {code}                                                    END OF FREQUENCY")
    return L


def make_atx(path):
    zen = np.arange(0, 91, 5.0)
    noazi1 = -0.1 * zen                      # mm, linear in zenith
    noazi2 = -0.2 * zen
    # asymmetric: PCV depends on azimuth: + 0.05*zen*cos(az)
    rows1 = [(az, noazi1 + 0.05 * zen * np.cos(np.radians(az))) for az in np.arange(0, 361, 90.0)]
    rows2 = [(az, noazi2 + 0.05 * zen * np.cos(np.radians(az))) for az in np.arange(0, 361, 90.0)]
    L = ["     1.4            M                                       ANTEX VERSION / SYST",
         "A                                                           PCV TYPE / REFANT",
         "                                                            END OF HEADER",
         "                                                            START OF ANTENNA",
         "TSTANT1         NONE                                        TYPE / SERIAL NO",
         "    90.0                                                    DAZI",
         "     0.0  90.0   5.0                                        ZEN1 / ZEN2 / DZEN"]
    L += _freq_block("G01", (1.0, 2.0, 100.0), noazi1, rows1)
    L += _freq_block("G02", (0.5, 1.0, 120.0), noazi2, rows2)
    L += ["                                                            END OF ANTENNA",
          "                                                            START OF ANTENNA",
          "BLOCK IIF           G01                 G063      2010-022A  TYPE / SERIAL NO",
          "     0.0                                                    DAZI",
          "     0.0  14.0   1.0                                        ZEN1 / ZEN2 / DZEN",
          "  2014     5    17     0     0    0.0000000                 VALID FROM"]
    nad = np.arange(0, 15, 1.0)
    for code in ("G01", "G02"):
        L += [f"   {code}                                                    START OF FREQUENCY",
              f"{394.0:10.2f}{0.0:10.2f}{1500.0:10.2f}                              NORTH / EAST / UP",
              "   NOAZI" + "".join(f"{v:8.2f}" for v in 0.5 * nad),
              f"   {code}                                                    END OF FREQUENCY"]
    L += ["                                                            END OF ANTENNA"]
    path.write_text("\n".join(_fix(x) for x in L) + "\n")


def test_antex_parse_and_pcv(tmp_path):
    p = tmp_path / "test.atx"
    make_atx(p)
    atx = an.read_antex(str(p))
    e, st = an.receiver_entry(atx, "TSTANT1         NONE")
    assert st == "EXACT"
    assert np.allclose(e.pco["G01"], [0.001, 0.002, 0.100])
    # node values and bilinear interpolation (asymmetric PCV)
    assert an._pcv_one(e, "G01", 30.0, 0.0) == pytest.approx((-3.0 + 1.5) * 1e-3)
    assert an._pcv_one(e, "G01", 30.0, 180.0) == pytest.approx((-3.0 - 1.5) * 1e-3)
    assert an._pcv_one(e, "G01", 32.5, 45.0) == pytest.approx((-3.25 + 0.5 * 1.625) * 1e-3)
    # radome fallback and missing type
    e2, st2 = an.receiver_entry(atx, "TSTANT1         SCIS")
    assert st2 == "RADOME_FALLBACK" and e2 is e
    assert an.receiver_entry(atx, "UNKNOWN         NONE") == (None, "NO_ANTENNA_CALIBRATION")


def test_receiver_correction_convention(tmp_path):
    """d_rho = -PCO_IF.e + PCV_IF (ANTEX 1.4 / RTKLIB antmodel convention)."""
    p = tmp_path / "test.atx"
    make_atx(p)
    atx = an.read_antex(str(p))
    e, _ = an.receiver_entry(atx, "TSTANT1         NONE")
    a1, a2 = an.if_coefficients()
    el, az = np.radians(40.0), np.radians(90.0)
    u = np.array([[np.cos(el) * np.sin(az), np.cos(el) * np.cos(az), np.sin(el)]])
    pco = a1 * np.array([0.002, 0.001, 0.100]) - a2 * np.array([0.001, 0.0005, 0.120])   # E, N, U
    pcv = a1 * an._pcv_one(e, "G01", 50.0, 90.0) - a2 * an._pcv_one(e, "G02", 50.0, 90.0)
    assert an.receiver_correction(e, u)[0] == pytest.approx(-(u[0] @ pco) + pcv, abs=1e-12)
    # zenith: correction = -PCO_U(IF)
    assert an.receiver_correction(e, np.array([[0, 0, 1.0]]))[0] == pytest.approx(-pco[2], abs=1e-12)


def test_satellite_entry_and_if(tmp_path):
    p = tmp_path / "test.atx"
    make_atx(p)
    atx = an.read_antex(str(p))
    assert an.satellite_entry(atx, "G01", ts.to_ns(2013, 1, 1)) is None     # before VALID FROM
    s = an.satellite_entry(atx, "G01", ts.to_ns(2024, 1, 1))
    assert s.svn == "G063" and an.block_of(s) == "BLOCK IIF"
    a1, a2 = an.if_coefficients()
    assert a1 - a2 == pytest.approx(1.0)
    assert an.pco_if(s)[2] == pytest.approx(1.5)                  # equal PCO on L1/L2 -> IF = same
    assert an.pcv_if(s, 10.0) == pytest.approx(5e-3)
