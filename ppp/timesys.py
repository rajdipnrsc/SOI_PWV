"""Time systems: GPST / UTC / TAI / TT, leap seconds, GPS week / day-of-year (TDS § 7.1).

Internal time representation [TDS § 7.1]: integer nanoseconds since the GPS epoch
1980-01-06 00:00:00 GPST (numpy int64 or Python int).  UTC instants are represented the same way on the
UTC *calendar* (a uniform count that skips leap seconds), so that GPST_ns = UTC_ns + (TAI-UTC - 19 s).

Relations [F]: GPST = TAI - 19 s;  TT = TAI + 32.184 s  ->  TT = GPST + 51.184 s;  UTC = TAI - (TAI-UTC).
"""
import numpy as np

from . import settings

NS = 1_000_000_000
DAY_NS = 86400 * NS
GPS_EPOCH_MJD = 44244                  # 1980-01-06 [F]
TT_MINUS_GPST_S = 51.184               # [F] 19 s + 32.184 s
GPST_MINUS_TAI_S = -19.0

_leap = {"table": None, "source": "builtin"}


# ---------------------------------------------------------------------------- calendar helpers
def _days_from_civil(y, m, d):
    """Days since 1970-01-01 for a proleptic Gregorian date (H. Hinnant's algorithm), integers."""
    y = y - (m <= 2)
    era = (y if y >= 0 else y - 399) // 400
    yoe = y - era * 400
    doy = (153 * (m + (-3 if m > 2 else 9)) + 2) // 5 + d - 1
    doe = yoe * 365 + yoe // 4 - yoe // 100 + doy
    return era * 146097 + doe - 719468


def _civil_from_days(z):
    z += 719468
    era = (z if z >= 0 else z - 146096) // 146097
    doe = z - era * 146097
    yoe = (doe - doe // 1460 + doe // 36524 - doe // 146096) // 365
    y = yoe + era * 400
    doy = doe - (365 * yoe + yoe // 4 - yoe // 100)
    mp = (5 * doy + 2) // 153
    d = doy - (153 * mp + 2) // 5 + 1
    m = mp + (3 if mp < 10 else -9)
    return y + (m <= 2), m, d


_MJD_1970 = 40587


def mjd_from_date(y, m, d):
    return _days_from_civil(y, m, d) + _MJD_1970


def date_from_mjd(mjd_int):
    return _civil_from_days(int(mjd_int) - _MJD_1970)


def to_ns(y, m, d, hh=0, mi=0, sec=0.0):
    """Calendar date/time (in whatever scale) -> integer ns since 1980-01-06 00:00 of that scale."""
    days = mjd_from_date(int(y), int(m), int(d)) - GPS_EPOCH_MJD
    whole = int(np.floor(sec))
    frac_ns = int(round((sec - whole) * NS))
    return ((days * 86400 + int(hh) * 3600 + int(mi) * 60 + whole) * NS) + frac_ns


def from_ns(ns):
    """Integer ns -> (y, m, d, hh, mi, sec_float) on the same scale's calendar."""
    ns = int(ns)
    days, rem = divmod(ns, DAY_NS)
    y, m, d = date_from_mjd(days + GPS_EPOCH_MJD)
    hh, rem = divmod(rem, 3600 * NS)
    mi, rem = divmod(rem, 60 * NS)
    return y, m, d, int(hh), int(mi), rem / NS


def iso(ns, sep="T", frac=0):
    y, m, d, hh, mi, s = from_ns(ns)
    if frac:
        return f"{y:04d}-{m:02d}-{d:02d}{sep}{hh:02d}:{mi:02d}:{s:0{3+frac}.{frac}f}"
    return f"{y:04d}-{m:02d}-{d:02d}{sep}{hh:02d}:{mi:02d}:{int(round(s)) if s < 59.5 else int(s):02d}"


def mjd(ns):
    """MJD (float) on the calendar of the scale of ns. Works on arrays."""
    return GPS_EPOCH_MJD + np.asarray(ns, dtype=np.float64) / DAY_NS


def gps_week_sow(ns):
    ns = int(ns)
    week, rem = divmod(ns, 7 * DAY_NS)
    return int(week), rem / NS


def gps_week_dow(ns):
    week, sow = gps_week_sow(ns)
    return week, int(sow // 86400)


def year_doy(ns):
    y, m, d, *_ = from_ns(ns)
    return y, mjd_from_date(y, m, d) - mjd_from_date(y, 1, 1) + 1


def day_start(ns):
    ns = int(ns)
    return ns - ns % DAY_NS


def decimal_year(ns):
    """Decimal year (e.g. 2024.0396) of an instant, from its calendar date (scale difference negligible)."""
    y, m, d, hh, mi, s = from_ns(ns)
    start = mjd_from_date(y, 1, 1)
    length = mjd_from_date(y + 1, 1, 1) - start
    frac = (mjd_from_date(y, m, d) - start + (hh * 3600 + mi * 60 + s) / 86400.0) / length
    return y + frac


def from_decimal_year(dy):
    y = int(np.floor(dy))
    start = mjd_from_date(y, 1, 1)
    length = mjd_from_date(y + 1, 1, 1) - start
    days = (dy - y) * length
    return (start - GPS_EPOCH_MJD) * DAY_NS + int(round(days * DAY_NS))


# ---------------------------------------------------------------------------- leap seconds
def _builtin_table():
    return [(mjd_from_date(*d), float(v)) for d, v in settings.LEAP_SECONDS_BUILTIN]


def set_leap_table(table, source):
    """table: list of (MJD_UTC from which TAI-UTC applies, TAI-UTC seconds), ascending."""
    _leap["table"] = sorted(table)
    _leap["source"] = source


def leap_source():
    return _leap["source"]


def leap_table():
    if _leap["table"] is None:
        _leap["table"] = _builtin_table()
    return _leap["table"]


def parse_leap_file(text):
    """Parse IERS Leap_Second.dat (MJD d m y TAI-UTC) or IETF/NTP leap-seconds.list. Returns table."""
    table = []
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        try:
            if len(parts) >= 5 and "." in parts[0]:            # IERS: 41317.0  1  1 1972  10
                table.append((int(float(parts[0])), float(parts[4])))
            elif len(parts) >= 2 and parts[0].isdigit() and len(parts[0]) >= 9:  # IETF: NTP seconds
                ntp = int(parts[0])
                mjd_utc = 15020 + ntp // 86400                # NTP epoch 1900-01-01 = MJD 15020
                table.append((mjd_utc, float(parts[1])))
        except ValueError:
            continue
    table = [t for t in table if t[1] >= 10]
    if len(table) < 20:
        raise ValueError("leap-second file not recognised")
    return sorted(table)


def tai_minus_utc(mjd_utc):
    """TAI-UTC (s) valid at a UTC MJD (float or int)."""
    val = 10.0
    for m0, v in leap_table():
        if mjd_utc >= m0:
            val = v
        else:
            break
    return val


def gpst_minus_utc(ns_utc):
    """GPST-UTC (s) at a UTC instant (ns on the UTC calendar)."""
    return tai_minus_utc(mjd(ns_utc)) - 19.0


def utc_to_gpst(ns_utc):
    ns_utc = int(ns_utc)
    return ns_utc + int(round(gpst_minus_utc(ns_utc) * NS))


def gpst_to_utc(ns_gps):
    """GPST ns -> UTC-calendar ns (iterated once to resolve the leap boundary)."""
    ns_gps = int(ns_gps)
    guess = ns_gps - int(round(gpst_minus_utc(ns_gps) * NS))
    return ns_gps - int(round(gpst_minus_utc(guess) * NS))


def gpst_to_tt_mjd(ns_gps):
    """TT expressed as MJD (float) for a GPST instant (array ok)."""
    return mjd(ns_gps) + TT_MINUS_GPST_S / 86400.0


def utc_hours_of_day(ns_utc):
    ns_utc = int(ns_utc)
    return (ns_utc % DAY_NS) / NS / 3600.0


def builtin_table_valid(ns_gps):
    y, m, d, *_ = from_ns(ns_gps)
    return (y, m, d) <= tuple(settings.LEAP_SECONDS_BUILTIN_VALID_UNTIL)
