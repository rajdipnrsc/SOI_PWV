"""Level-0: RINEX 3 observation/navigation reading (TDS § 13.2)."""
import os

import numpy as np
import pytest

from tests.conftest import SAMPLE
from ppp import rinex
from ppp import timesys as ts

OBS = """     3.04           OBSERVATION DATA    M                   RINEX VERSION / TYPE
TEST                                                        MARKER NAME
                    TRM159900.00    SCIS                    ANT # / TYPE
  1208297.7644  5967332.9995  1895462.8434                  APPROX POSITION XYZ
        0.1000        0.0000        0.0000                  ANTENNA: DELTA H/E/N
G    6 C1C L1C C2W L2W S1C S2W                              SYS / # / OBS TYPES
E    2 C1X L1X                                              SYS / # / OBS TYPES
    30.000                                                  INTERVAL
  2024     1    15     0     0    0.0000000     GPS         TIME OF FIRST OBS
                                                            END OF HEADER
{DATA}"""


def _f(v, lli=" ", ssi=" "):
    return (f"{v:14.3f}" if v is not None else " " * 14) + lli + ssi


def _data():
    L = ["> 2024 01 15 00 00  0.0000000  0  3",
         "G02" + _f(24501904.680, " ", "7") + _f(128756420.123) + _f(24501916.375) + _f(100331253.372, "1")
         + _f(45.0) + _f(40.0),
         "E05" + _f(22000000.0) + _f(115610000.0),
         "G05" + _f(22501904.680) + _f(118246420.123) + _f(22501916.375) + _f(92141253.372),
         "> 2024 01 15 00 00 15.0000000  0  1",
         "G02" + _f(24501904.680) + _f(128756420.123) + _f(24501916.375) + _f(100331253.372),
         "> 2024 01 15 00 00 30.0000000  0  1",
         "G02" + _f(24501904.680) + _f(128756420.123) + _f(24501916.375) + _f(100331253.372)]
    return "\n".join(L) + "\n"


def test_read_obs_synthetic(tmp_path):
    p = tmp_path / "TEST0150.24o"
    p.write_text(OBS.replace("{DATA}", _data()))
    d = rinex.read_obs(str(p), decimate_s=30)
    assert d.sats == ["G02", "G05"]                      # Galileo dropped (GPS only v1)
    assert len(d.epochs) == 2                            # 15-s epoch decimated
    assert d.epochs[0] == ts.to_ns(2024, 1, 15)
    assert d.header.antenna_type == "TRM159900.00    SCIS"
    assert d.header.delta_hen == (0.1, 0.0, 0.0)
    assert d.obs["L1C"][0, 0] == pytest.approx(128756420.123)
    assert d.lli["L2W"][0, 0] == 1 and d.snr_flag["C1C"][0, 0] == 7
    assert np.isnan(d.obs["S1C"][0, 1])


def test_rinex2_rejected(tmp_path):
    p = tmp_path / "x.o"
    p.write_text("     2.11           OBSERVATION DATA    G                   RINEX VERSION / TYPE\n"
                 "                                                            END OF HEADER\n")
    with pytest.raises(rinex.RinexError):
        rinex.read_obs(str(p))


@pytest.mark.skipif(not os.path.exists(os.path.join(SAMPLE, "HYDE015.24o")), reason="sample data absent")
def test_read_sample_files():
    d = rinex.read_obs(os.path.join(SAMPLE, "HYDE015.24o"), decimate_s=30)
    assert len(d.epochs) == 2880 and d.header.marker_name == "HYDE"
    nav = rinex.read_nav(os.path.join(SAMPLE, "HYDE015.24n"))
    assert len(nav) >= 28
    first, last = rinex.first_last_epoch(os.path.join(SAMPLE, "HYDE015.24o"))
    assert first == ts.to_ns(2024, 1, 15) and last == ts.to_ns(2024, 1, 15, 23, 59, 30)
