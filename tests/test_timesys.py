"""Level-0: time conversions incl. leap seconds and GPS week rollovers (TDS § 30.1, exact)."""
from ppp import timesys as ts


def test_gps_epoch_and_rollovers():
    assert ts.to_ns(1980, 1, 6) == 0
    assert ts.gps_week_dow(ts.to_ns(1999, 8, 22)) == (1024, 0)     # first rollover [F]
    assert ts.gps_week_dow(ts.to_ns(2019, 4, 7)) == (2048, 0)      # second rollover [F]
    assert ts.gps_week_dow(ts.to_ns(2022, 11, 27)) == (2238, 0)    # IGS20 / long names [F]
    assert ts.gps_week_dow(ts.to_ns(2018, 6, 3)) == (2004, 0)      # CODE integer clocks [F]


def test_calendar_roundtrip_and_doy():
    t = ts.to_ns(2024, 2, 29, 23, 59, 30.5)
    assert ts.from_ns(t) == (2024, 2, 29, 23, 59, 30.5)
    assert ts.year_doy(ts.to_ns(2024, 12, 31)) == (2024, 366)
    assert ts.year_doy(ts.to_ns(2023, 12, 31)) == (2023, 365)
    assert abs(ts.decimal_year(ts.to_ns(2005, 1, 1)) - 2005.0) < 1e-12
    assert ts.mjd(ts.to_ns(2000, 1, 1, 12)) == 51544.5


def test_leap_seconds():
    ts.set_leap_table(ts._builtin_table(), "builtin")
    # GPST - UTC: 13 s in 2005, 17 s in 2016, 18 s since 2017-01-01 [F]
    assert ts.gpst_minus_utc(ts.to_ns(2005, 6, 1)) == 13.0
    assert ts.gpst_minus_utc(ts.to_ns(2016, 12, 31, 23, 59, 59)) == 17.0
    assert ts.gpst_minus_utc(ts.to_ns(2017, 1, 1)) == 18.0
    t_utc = ts.to_ns(2024, 1, 15, 0, 0, 0)
    t_gps = ts.utc_to_gpst(t_utc)
    assert t_gps - t_utc == 18 * ts.NS
    assert ts.gpst_to_utc(t_gps) == t_utc
    # across the 2016/2017 leap second
    t = ts.utc_to_gpst(ts.to_ns(2017, 1, 1, 0, 0, 5))
    assert ts.gpst_to_utc(t) == ts.to_ns(2017, 1, 1, 0, 0, 5)
    # TT = GPST + 51.184 s
    assert abs((ts.gpst_to_tt_mjd(0) - 44244.0) * 86400 - 51.184) < 1e-6


def test_parse_leap_files():
    iers = "\n".join(["#  File expires on 28 December 2026"] +
                     [f"   {m}.0    1  1 {y}      {v}" for (y, m, v) in
                      [(1972, 41317, 10), (1972, 41499, 11), (1973, 41683, 12), (1974, 42048, 13), (1975, 42413, 14),
                       (1976, 42778, 15), (1977, 43144, 16), (1978, 43509, 17), (1979, 43874, 18), (1980, 44239, 19),
                       (1981, 44786, 20), (1982, 45151, 21), (1983, 45516, 22), (1985, 46247, 23), (1988, 47161, 24),
                       (1990, 47892, 25), (1991, 48257, 26), (1992, 48804, 27), (1993, 49169, 28), (1994, 49534, 29),
                       (1996, 50083, 30), (1997, 50630, 31), (1999, 51179, 32), (2006, 53736, 33), (2009, 54832, 34),
                       (2012, 56109, 35), (2015, 57204, 36), (2017, 57754, 37)]])
    tab = ts.parse_leap_file(iers)
    assert tab[-1] == (57754, 37.0)
    ietf = "\n".join(["#$ 3676924800"] + [f"{(m - 15020) * 86400}\t{int(v)}\t# x" for m, v in tab])
    assert ts.parse_leap_file(ietf)[-1] == (57754, 37.0)
