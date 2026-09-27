"""Station products: CSV, NetCDF-4 (CF-1.10 + ACDD-1.3), SINEX_TRO (optional), quick-look plot (TDS § 16-17, § 33.5).

Units: delays, gradients and PWV in mm; clock in ns; heights m; angles degrees (TDS § 17).
"""
import json

import numpy as np

from . import settings
from . import timesys as ts

STATEMENT = ("Adjacent 5-minute values are strongly correlated through the random-walk model and smoothing. "
             "They are samples of a smoothed estimate, not independent 5-minute atmospheric observations. The "
             "effective independent time resolution is estimated at 15-30 min [E] and must be quantified in "
             "validation.")                              # mandatory statement (TDS § 6.5)

ZTD_COLUMNS = ["timestamp_utc", "timestamp_gpst", "station", "latitude", "longitude", "height", "ZTD", "sigma_ZTD",
               "ZHD", "ZWD", "sigma_ZWD", "G_N", "G_E", "receiver_clock", "n_sat", "n_rejected", "n_amb_fixed",
               "convergence_flag", "AR_status", "QC_flag", "product_AC", "product_tier", "product_files",
               "coordinate_mode", "coordinate_source", "mapping_function", "ZWD_process_noise", "software_version",
               "solution_type", "zhd_source", "manifest_id", "config_hash", "elevation_cutoff"]
PWV_EXTRA = ["PWV", "sigma_PWV_ppp", "sigma_PWV_conv", "sigma_PWV_total", "P_ant", "P_source", "T_m", "Tm_source",
             "Pi", "ZHD_met", "ZWD_met", "H_ant", "constants_set"]


def base_name(station, day_ns, solution, tier, kind, ext):
    y, doy = ts.year_doy(day_ns)
    return f"{station}_{y:04d}{doy:03d}_{kind}.{ext}" if solution is None else \
        f"{station}_{y:04d}{doy:03d}_{solution}_{tier}_{kind}_v{settings.SCHEMA_VERSION}.{ext}"


def ztd_rows(fm, meta):
    """Build the list of dict rows of the § 16 schema from a FiveMin extraction and constant metadata."""
    rows = []
    for i in range(len(fm.t_utc)):
        ztd = fm.zhd0[i] + fm.zwd[i]
        r = {
            "timestamp_utc": ts.iso(int(fm.t_utc[i])) + "Z",
            "timestamp_gpst": ts.iso(int(fm.t_gpst[i])),
            "station": meta["station"], "latitude": meta["lat"], "longitude": meta["lon"], "height": meta["h"],
            "ZTD": ztd * 1e3, "sigma_ZTD": fm.sig[i] * 1e3, "ZHD": fm.zhd0[i] * 1e3, "ZWD": fm.zwd[i] * 1e3,
            "sigma_ZWD": fm.sig[i] * 1e3, "G_N": fm.gn[i] * 1e3, "G_E": fm.ge[i] * 1e3,
            "receiver_clock": fm.clk_ns[i], "n_sat": int(fm.n_sat[i]), "n_rejected": int(fm.n_rej[i]),
            "n_amb_fixed": 0, "convergence_flag": fm.conv[i], "AR_status": meta.get("AR_status", "AR_NOT_AVAILABLE"),
            "QC_flag": int(fm.qc[i]),
        }
        for k in ("product_AC", "product_tier", "product_files", "coordinate_mode", "coordinate_source",
                  "mapping_function", "ZWD_process_noise", "software_version", "solution_type", "zhd_source",
                  "manifest_id", "config_hash", "elevation_cutoff"):
            r[k] = meta[k]
        rows.append(r)
    return rows


def _fmt(v):
    if isinstance(v, (float, np.floating)):
        return "" if not np.isfinite(v) else f"{v:.4f}"
    return str(v)


def write_csv(path, rows, columns, comment):
    with open(path, "w") as fh:
        fh.write(f"# {comment}\n")
        fh.write(f"# schema_version={settings.SCHEMA_VERSION}; units: delays/gradients/PWV mm, clock ns, heights m, "
                 f"angles deg; {STATEMENT}\n")
        fh.write(",".join(columns) + "\n")
        for r in rows:
            vals = []
            for c in columns:
                s = _fmt(r.get(c, ""))
                if "," in s or '"' in s:
                    s = '"' + s.replace('"', '""') + '"'
                vals.append(s)
            fh.write(",".join(vals) + "\n")


_NC_UNITS = {"ZTD": "mm", "sigma_ZTD": "mm", "ZHD": "mm", "ZWD": "mm", "sigma_ZWD": "mm", "G_N": "mm", "G_E": "mm",
             "receiver_clock": "ns", "PWV": "mm", "sigma_PWV_ppp": "mm", "sigma_PWV_conv": "mm",
             "sigma_PWV_total": "mm", "P_ant": "hPa", "T_m": "K", "Pi": "1", "ZHD_met": "mm", "ZWD_met": "mm",
             "n_sat": "1", "n_rejected": "1", "n_amb_fixed": "1", "QC_flag": "1"}
_NC_LONG = {"ZTD": "zenith total delay at ARP height (smoothed)", "sigma_ZTD": "formal 1-sigma of ZTD (unscaled)",
            "ZHD": "a priori zenith hydrostatic delay used in PPP", "ZWD": "ZTD minus a priori ZHD",
            "sigma_ZWD": "formal 1-sigma of ZWD", "G_N": "north horizontal gradient (Chen and Herring)",
            "G_E": "east horizontal gradient (Chen and Herring)", "receiver_clock": "receiver clock estimate",
            "PWV": "precipitable water vapour", "QC_flag": "quality control bitmask (see flag_meanings)"}


def write_netcdf(path, rows, columns, attrs, lat, lon, alt):
    """NetCDF-4 CF-1.10 + ACDD-1.3 timeSeries file (TDS § 17)."""
    import netCDF4
    with netCDF4.Dataset(path, "w", format="NETCDF4") as nc:
        nc.createDimension("time", len(rows))
        tv = nc.createVariable("time", "i8", ("time",))
        tv.units = "seconds since 1970-01-01T00:00:00Z"
        tv.calendar = "standard"
        tv.standard_name = "time"
        tv.comment = "UTC; leap seconds not represented"
        epoch1970 = ts.to_ns(1970, 1, 1)
        tv[:] = [(ts.to_ns(*_parse_iso(r["timestamp_utc"])) - epoch1970) // ts.NS for r in rows]
        for nm, val, un, sn in (("lat", lat, "degree_north", "latitude"), ("lon", lon, "degree_east", "longitude"),
                                ("alt", alt, "m", "height")):
            v = nc.createVariable(nm, "f8")
            v.units, v.standard_name = un, sn
            v[...] = val
        sid = nc.createVariable("station_name", str)
        sid.cf_role = "timeseries_id"
        sid[0] = str(attrs.get("station", ""))
        for c in columns:
            if c in ("timestamp_utc", "timestamp_gpst", "latitude", "longitude", "height"):
                continue
            vals = [r.get(c) for r in rows]
            if all(isinstance(x, (int, np.integer)) for x in vals) and c in ("n_sat", "n_rejected", "n_amb_fixed",
                                                                               "QC_flag"):
                dt = "u4" if c == "QC_flag" else "i2"
                v = nc.createVariable(c, dt, ("time",))
                v[:] = np.array(vals, dtype=np.uint32 if dt == "u4" else np.int16)
            elif all(isinstance(x, (float, int, np.floating, np.integer)) for x in vals):
                v = nc.createVariable(c, "f8", ("time",), fill_value=-9999.0)
                arr = np.array(vals, dtype=float)
                v[:] = np.where(np.isfinite(arr), arr, -9999.0)
            else:
                v = nc.createVariable(c, str, ("time",))
                for i, x in enumerate(vals):
                    v[i] = str(x)
            if c in _NC_UNITS:
                v.units = _NC_UNITS[c]
            if c in _NC_LONG:
                v.long_name = _NC_LONG[c]
            v.coordinates = "time lat lon alt station_name"
            if c == "QC_flag":
                names = sorted(settings.QC_BITS, key=settings.QC_BITS.get)
                v.flag_masks = np.array([1 << settings.QC_BITS[n] for n in names], dtype=np.uint32)
                v.flag_meanings = " ".join(names)
        nc.Conventions = "CF-1.10, ACDD-1.3"
        nc.featureType = "timeSeries"
        nc.schema_version = settings.SCHEMA_VERSION
        nc.comment = STATEMENT
        for k, v in attrs.items():
            setattr(nc, k, v if isinstance(v, (str, int, float)) else json.dumps(v, default=str))


def _parse_iso(s):
    s = s.rstrip("Z")
    d, t = s.split("T")
    y, m, dd = (int(x) for x in d.split("-"))
    hh, mi, sec = t.split(":")
    return y, m, dd, int(hh), int(mi), float(sec)


def write_sinex_tro(path, station, fm, meta, xyz):
    """SINEX_TRO 2.00-style file (optional product). Field names TROTOT/STDDEV/TGNTOT/TGETOT in mm.

    VERIFY item V-030: block layout follows the SINEX_TRO 2.00 draft as understood at implementation; the
    file is an optional interoperability product, not the primary archive.
    """
    def ep(t):
        y, doy = ts.year_doy(int(t))
        sod = int(round((int(t) % ts.DAY_NS) / ts.NS))
        return f"{y:04d}:{doy:03d}:{sod:05d}"
    t0, t1 = fm.t_gpst[0], fm.t_gpst[-1]
    lines = [f"%=TRO 2.00 {settings.SOFTWARE_NAME[:3].upper()} {ep(t1)} COD {ep(t0)} {ep(t1)} P "
             f"{len(fm.t_utc):05d} 0 T"]
    lines += ["+FILE/REFERENCE", f" DESCRIPTION        {settings.SOFTWARE_NAME} float PPP troposphere",
              " OUTPUT             Station ZTD, 5-min", f" SOFTWARE           {meta['software_version']}",
              "-FILE/REFERENCE", "+TROP/DESCRIPTION",
              f" ELEVATION CUTOFF ANGLE     {meta['elevation_cutoff']:5.1f}",
              f" SAMPLING INTERVAL          {int(settings.PROCESSING_INTERVAL_S):5d}",
              f" SAMPLING TROP              {settings.OUTPUT_INTERVAL_S:5d}",
              f" TROP MAPPING FUNCTION      {meta['mapping_function']}",
              " SOLUTION_FIELDS_1          TROTOT STDDEV TGNTOT STDDEV TGETOT STDDEV",
              "-TROP/DESCRIPTION", "+TROP/STA_COORDINATES",
              "*SITE PT SOLN T __STA_X_____ __STA_Y_____ __STA_Z_____ SYSTEM REMRK",
              f" {station[:4]:4s}  A    1 P {xyz[0]:12.3f} {xyz[1]:12.3f} {xyz[2]:12.3f} {meta.get('frame', ''):6s} "
              f"{meta.get('coordinate_mode', '')[:5]}", "-TROP/STA_COORDINATES", "+TROP/SOLUTION",
              "*SITE ____EPOCH_____ TROTOT STDDEV TGNTOT STDDEV TGETOT STDDEV"]
    for i in range(len(fm.t_utc)):
        if not np.isfinite(fm.zwd[i]):
            continue
        lines.append(f" {station[:4]:4s} {ep(fm.t_gpst[i])} {(fm.zhd0[i] + fm.zwd[i]) * 1e3:6.1f} {fm.sig[i] * 1e3:6.1f} "
                     f"{fm.gn[i] * 1e3:6.2f} {fm.sig[i] * 1e3:6.2f} {fm.ge[i] * 1e3:6.2f} {fm.sig[i] * 1e3:6.2f}")
    lines += ["-TROP/SOLUTION", "%=ENDTROP"]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def quicklook(path, station, day_label, fm, pwv=None, title_extra=""):
    """Quick-look PNG: ZTD, PWV, satellites and flags as single-series small multiples on a shared time axis."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, muted, grid, surface = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    series = "#2a78d6"                              # categorical slot 1 (reference palette)
    band = "#2a78d6"
    hrs = ((fm.t_utc - fm.t_utc[0]) / ts.NS / 3600.0)
    panels = [("ZTD (mm)", (fm.zhd0 + fm.zwd) * 1e3, fm.sig * 1e3)]
    if pwv is not None:
        panels.append(("PWV (mm)", pwv.pwv, pwv.sig_total))
    panels.append(("Satellites used", fm.n_sat.astype(float), None))
    n = len(panels) + 1
    fig, axes = plt.subplots(n, 1, figsize=(10, 2.1 * n), sharex=True, facecolor=surface,
                             gridspec_kw={"height_ratios": [3] * (n - 1) + [1.2]})
    for ax, (lab, y, s) in zip(axes, panels):
        ax.set_facecolor(surface)
        if s is not None:
            ax.fill_between(hrs, y - s, y + s, color=band, alpha=0.18, linewidth=0)
        ax.plot(hrs, y, color=series, linewidth=2 if s is not None else 1.5,
                drawstyle="steps-mid" if s is None else "default")
        ax.set_title(lab, loc="left", fontsize=10, color=ink)
        ax.grid(True, color=grid, linewidth=0.8)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(grid)
        ax.tick_params(colors=muted, labelsize=8)
    ax = axes[-1]
    ax.set_facecolor(surface)
    shade = {"CONVERGED": "#c3c2b7", "NOT_CONVERGED": "#eda100", "EDGE": "#52514e", "REINIT": "#eb6834",
             "NO_DATA": "#ffffff"}
    for i, c in enumerate(fm.conv):
        ax.axvspan(hrs[i] - 2.5 / 60, hrs[i] + 2.5 / 60, color=shade.get(c, "#ffffff"), linewidth=0)
    ax.set_yticks([])
    ax.set_title("Convergence flag (grey converged, dark edge, amber not converged, orange re-init, blank no data)",
                 loc="left", fontsize=9, color=ink)
    for sp in ax.spines.values():
        sp.set_color(grid)
    ax.set_xlabel("hours since " + ts.iso(int(fm.t_utc[0]))[:16] + " UTC", color=muted, fontsize=9)
    ax.tick_params(colors=muted, labelsize=8)
    fig.suptitle(f"{station} {day_label} {title_extra}", x=0.01, ha="left", fontsize=11, color=ink)
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=surface)
    plt.close(fig)
