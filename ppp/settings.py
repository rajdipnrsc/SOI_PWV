"""ALL defaults, thresholds, URL/filename templates and constants in one place (TDS § 33.2).

Units are stated for every value.  References: TDS section numbers in square brackets, external
sources by name.  Labels follow the TDS: [F] fact, [A] engineering assumption, [E] estimate.
VERIFY items that were resolved during implementation are recorded in docs/verify_log.md.

A small YAML (or JSON) file may override any UPPER_CASE value here (TDS § 33.2 simplicity rule):
place ``pwv_ppp_settings.yaml`` / ``.json`` in the working directory; see ``apply_overrides``.
"""
import json
import os

SOFTWARE_NAME = "pwv_ppp"
SOFTWARE_VERSION = "0.9.0"          # semantic version; git commit is appended at run time
SCHEMA_VERSION = "1.0"              # station product schema (TDS § 16, § 20)

# ----------------------------------------------------------------------------------------------
# Physical constants [TDS § 3]
# ----------------------------------------------------------------------------------------------
C_LIGHT = 299792458.0               # m/s [F] SI
OMEGA_E = 7.2921151467e-5           # rad/s [F] IS-GPS-200 / WGS84
GM_EARTH = 3.986004418e14           # m^3/s^2 [F] IERS 2010 / WGS84
GM_SUN = 1.32712442099e20           # m^3/s^2 [F] IERS 2010 (TDB-compatible, used only for eclipse geometry)
RE_WGS84 = 6378137.0                # m, semi-major axis (GRS80 = WGS84 to < 0.1 mm)
F_GRS80 = 1.0 / 298.257222101       # GRS80 flattening [F]
AU = 149597870700.0                 # m [F] IAU 2012
R_EARTH_SHADOW = 6378137.0          # m, Earth radius used for cylindrical shadow test [A]

# Signal definition table [TDS § 3]: frequencies are never hard-coded in the estimator.
# key = (system, band) -> Hz
FREQ = {
    ("G", "1"): 1575.42e6,          # GPS L1 [F] IS-GPS-200
    ("G", "2"): 1227.60e6,          # GPS L2 [F]
    ("G", "5"): 1176.45e6,          # GPS L5 [F] (not used in v1)
    ("E", "1"): 1575.42e6,          # Galileo E1 (future)
    ("E", "5"): 1176.45e6,          # Galileo E5a (future)
}
SYSTEMS = ("G",)                    # constellation v1 = GPS only [TDS F1]
BANDS = {"G": ("1", "2")}           # IF combination bands per system

# Observable selection priority (TDS § 13.3), configurable, no vendor hard-coding
CODE_PRIORITY = {"G": {"1": ["C1W", "C1C", "C1X"], "2": ["C2W", "C2L", "C2S", "C2X"]}}
PHASE_PRIORITY = {"G": {"1": ["L1C", "L1W", "L1X"], "2": ["L2W", "L2L", "L2S", "L2X"]}}
# clock-reference signals of CODE GPS clocks: their OSBs cancel in the IF combination (TDS § 11.2 note 1)
CLOCK_REFERENCE_CODES = {"G": ("C1W", "C2W")}

# ----------------------------------------------------------------------------------------------
# Time [TDS § 7.1]
# ----------------------------------------------------------------------------------------------
# Built-in leap-second table (UTC date from which TAI-UTC applies, TAI-UTC in s) [F: IERS Bulletin C].
# Used only when no downloaded/cached leap-second file is available; the source used is recorded.
LEAP_SECONDS_BUILTIN = [
    ((1980, 1, 1), 19), ((1981, 7, 1), 20), ((1982, 7, 1), 21), ((1983, 7, 1), 22),
    ((1985, 7, 1), 23), ((1988, 1, 1), 24), ((1990, 1, 1), 25), ((1991, 1, 1), 26),
    ((1992, 7, 1), 27), ((1993, 7, 1), 28), ((1994, 7, 1), 29), ((1996, 1, 1), 30),
    ((1997, 7, 1), 31), ((1999, 1, 1), 32), ((2006, 1, 1), 33), ((2009, 1, 1), 34),
    ((2012, 7, 1), 35), ((2015, 7, 1), 36), ((2017, 1, 1), 37),
]
# [A] no leap second after 2017-01-01 was known at implementation; the downloaded IERS/IETF file is
# authoritative. After this date the built-in table is used only with a WARNING.
LEAP_SECONDS_BUILTIN_VALID_UNTIL = (2026, 6, 30)
LEAP_SECOND_URLS = [
    "https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat",
    "https://data.iana.org/time-zones/tzdb/leap-seconds.list",
]
LEAP_SECOND_REFRESH_DAYS = 30        # re-download monthly [TDS § 33.6]

# ----------------------------------------------------------------------------------------------
# Products: CODE-only policy [TDS § 11]
# ----------------------------------------------------------------------------------------------
PRODUCT_FAMILIES = ["COD0OPSFIN", "CODMOPSRAP", "COD0OPSRAP"]     # fallback order [TDS § 11.1]
FAMILY_TIER = {"COD0OPSFIN": "FINAL", "CODMOPSRAP": "RAPID_M", "COD0OPSRAP": "RAPID_0"}
FAMILY_VERSION_CHAR = {"COD0OPSFIN": "0", "CODMOPSRAP": "M", "COD0OPSRAP": "0"}
FORBIDDEN_FAMILIES = ["COD0OPSULT"]  # never used for PPP [TDS § 11.1]

# Long filenames (since GPS week 2238, 27 Nov 2022). {family} {yyyy} {ddd}. Compression variants tried
# in order (rapid files are listed uncompressed in CODE's table) [TDS § 11.2 VERIFY -> docs/verify_log.md].
LONG_NAMES = {
    "SP3": "{family}_{yyyy}{ddd}0000_01D_05M_ORB.SP3",
    "CLK": "{family}_{yyyy}{ddd}0000_01D_30S_CLK.CLK",
    "CLK05S": "{family}_{yyyy}{ddd}0000_01D_05S_CLK.CLK",
    "ERP": "{family}_{yyyy}{ddd}0000_01D_01D_ERP.ERP",
    "OSB": "{family}_{yyyy}{ddd}0000_01D_01D_OSB.BIA",
    "ATT": "{family}_{yyyy}{ddd}0000_01D_30S_ATT.OBX",
    "TRO": "{family}_{yyyy}{ddd}0000_01D_01H_TRO.TRO",
}
# ATT is not provided by CODMOPSRAP [TDS § 11.2]
FAMILY_PRODUCTS = {
    "COD0OPSFIN": ["SP3", "CLK", "ERP", "OSB", "ATT"],
    "CODMOPSRAP": ["SP3", "CLK", "ERP", "OSB"],
    "COD0OPSRAP": ["SP3", "CLK", "ERP", "OSB", "ATT"],
}
# Legacy short names (before 27 Nov 2022) [TDS § 11.2; VERIFY -> verify_log]. {wwww}{d}
LEGACY_NAMES = {
    "COD0OPSFIN": {"SP3": "COD{wwww}{d}.EPH", "CLK": "COD{wwww}{d}.CLK", "ERP": "COD{wwww}{d}.ERP",
                   "OSB": "COD{wwww}{d}.BIA", "ATT": "COD{wwww}{d}.OBX"},
    "CODMOPSRAP": {"SP3": "COD{wwww}{d}.EPH_M", "CLK": "COD{wwww}{d}.CLK_M", "ERP": "COD{wwww}{d}.ERP_M",
                   "OSB": "COD{wwww}{d}.BIA_M"},
    "COD0OPSRAP": {"SP3": "COD{wwww}{d}.EPH_R", "CLK": "COD{wwww}{d}.CLK_R", "ERP": "COD{wwww}{d}.ERP_R",
                   "OSB": "COD{wwww}{d}.BIA_R"},
}
LONG_NAME_START = (2022, 11, 27)     # GPS week 2238 [F]
COMPRESSION_SUFFIXES = [".gz", ".Z", ""]
# Ordered source list: URL templates with tokens {yyyy} {ddd} {wwww} {d} {family} {filename}
# [TDS § 12.2, § 33.6]. CDDIS (Earthdata login via ~/.netrc) is tried last.
PRODUCT_SOURCES = [
    "http://ftp.aiub.unibe.ch/CODE/{yyyy}/{filename}",
    "https://zhw-b.s3.cloud.switch.ch/aiub/CODE/{yyyy}/{filename}",
    "https://cddis.nasa.gov/archive/gnss/products/{wwww}/{filename}",
]
# Published checksum manifests per source (TDS § 12.4 item 2): source template -> [(algorithm, manifest URL template)].
# CDDIS publishes SHA512SUMS (earlier MD5SUMS) per weekly product directory (VERIFY V-035).
CHECKSUM_MANIFESTS = {
    "https://cddis.nasa.gov/archive/gnss/products/{wwww}/{filename}": [
        ("sha512", "https://cddis.nasa.gov/archive/gnss/products/{wwww}/SHA512SUMS"),
        ("md5", "https://cddis.nasa.gov/archive/gnss/products/{wwww}/MD5SUMS")],
}
REQUIRED_FLOAT = ["SP3", "CLK", "ERP"]          # + OSB when observables are not clock-reference [TDS § 11.3]
REQUIRED_AR = ["SP3", "CLK", "ERP", "OSB"]
MAX_EXPECTED_LATENCY_H = {"COD0OPSFIN": 21 * 24, "CODMOPSRAP": 48, "COD0OPSRAP": 48}  # [A, TDS § 11.4]
INTEGER_CLOCK_START = (2018, 6, 3)               # GPS week 2004 [F]
MIN_GPS_SATS_IN_PRODUCT = 24                     # [A, TDS § 12.4]
ATT_FROM_OTHER_SOLUTION = False                  # [TDS § 11.2 note 2]

# ANTEX [TDS § 33.6]: the exact version named in the clock header is tried first (pcv_archive), then
# the current file of the same family (recorded as ANTEX_VERSION_MISMATCH, see docs/decisions.md D-004).
ANTEX_SOURCES = [
    "https://files.igs.org/pub/station/general/pcv_archive/{name}.atx",
    "https://files.igs.org/pub/station/general/{family}.atx",
]
ANTEX_DEFAULT_NAME = "igs20"

# Troposphere grids [TDS § 6, § 33.6]
VMF3_URL = "https://vmf.geo.tuwien.ac.at/trop_products/GRID/1x1/VMF3/VMF3_OP/{yyyy}/VMF3_{yyyy}{mm}{dd}.H{hh}"
VMF3_OROGRAPHY_URL = "https://vmf.geo.tuwien.ac.at/station_coord_files/orography_ell_1x1"
GPT3_URL = "https://vmf.geo.tuwien.ac.at/codes/gpt3_1.grd"
VMF3_HOURS = (0, 6, 12, 18)

# Download mechanics [TDS § 12.3]
HTTP_CONNECT_TIMEOUT_S = 20
HTTP_READ_TIMEOUT_S = 120
HTTP_RETRIES = 5
HTTP_BACKOFF_MAX_S = 300
CA_BUNDLE = os.environ.get("REQUESTS_CA_BUNDLE") or os.environ.get("SSL_CERT_FILE") or None

# ----------------------------------------------------------------------------------------------
# Processing window, epochs, cutoff [TDS § 5.10, § 13]
# ----------------------------------------------------------------------------------------------
PROCESSING_INTERVAL_S = 30.0         # epochs coincide with CODE 30-s clocks [TDS § 13.4]
OVERLAP_H = 3.0                      # window [D-3h, D+1+3h] when adjacent RINEX exists [TDS § 5.10]
ELEVATION_CUTOFF_DEG = 7.0           # [TDS § 6.1, F6]
SP3_LAGRANGE_NODES = 11              # degree 10 [TDS § 7.2]
LIGHT_TIME_TOL_S = 1e-12             # [TDS § 7.1]

# ----------------------------------------------------------------------------------------------
# Estimator [TDS § 5, § 6]
# ----------------------------------------------------------------------------------------------
SIGMA_ZWD_M_SQRT_H = 0.006           # random walk, m/sqrt(h) [A, TDS § 6.4]
SIGMA_GRAD_M_SQRT_H = 0.0005         # m/sqrt(h) [A]
ESTIMATE_GRADIENTS = True
SIGMA_CLK0_M = 100.0                 # receiver clock white noise [A]
INIT_SIGMA_ZWD_M = 0.10
INIT_SIGMA_GRAD_M = 0.001
SIGMA_AMB0_M = 20.0                  # [A, TDS § 5.5]
INIT_SIGMA_XYZ_STATIC_M = 1.0
CONSTRAINED_SIGMA_NEU_M = (0.003, 0.003, 0.005)   # [A, TDS § 5.2]
SIGMA_PHASE_AB_M = (0.003, 0.003)    # raw phase sigma_a, sigma_b [A, TDS § 5.4]
SIGMA_CODE_AB_M = (0.3, 0.3)         # raw code sigma_a, sigma_b [A]
SNR_WEIGHTING = False                # optional SNR variance model, off by default [TDS § 5.4]
SNR_WEIGHT_REF_DBHZ = None           # None: per-day normalisation to the elevation model (D-011) [A]
SNR_WEIGHT_MIN_EL_DEG = 5.0          # SNR variance factor clipped to [1, 1/sin^2(5 deg)] [A]
GAP_RESET_S = 120.0                  # [A, TDS § 5.7]
GAP_RESET_ALL_S = 600.0
WINDOW_BREAK_S = 3600.0
N_CONSEC_REJECT = 3
MIN_ARC_S = 600.0                    # [A, TDS § 5.8]
IGG_K0 = 2.0                         # [A, TDS § 13.10]
IGG_K1 = 5.0
EPOCH_REJECT_FRACTION = 0.30
POSTFIT_PHASE_K = 4.0
POSTFIT_PHASE_MIN_M = 0.04
POSTFIT_CODE_K = 4.0
POSTFIT_CODE_MIN_M = 4.0
MAX_EDIT_PASSES = 3
STATIC_RELINEARISE_PASSES = 1        # re-run once with converged coordinate [TDS § 5.6]

# Output extraction and flags [TDS § 6.5-6.6, § 16.1]
OUTPUT_INTERVAL_S = 300              # 5-min UTC marks
EDGE_WINDOW_S = 3600.0
REINIT_WINDOW_S = 1800.0
SIGMA_CONV_M = 0.010
SIGMA_HIGH_M = 0.015
FEW_SATS = 5
HIGH_REJECT_FRACTION = 0.20
COORD_PULL_LIMIT = 3.0

QC_BITS = {
    "NO_DATA": 0, "FEW_SATS": 1, "EDGE": 2, "REINIT": 3, "IONO_ACTIVE": 4, "HIGH_REJECT": 5,
    "COORD_PULL": 6, "RADOME_FALLBACK": 7, "ORBIT_EDGE": 8, "CLOCK_INTERP": 9, "ZHD_CLIMATOLOGY": 10,
    "ECLIPSE_EXCLUSION": 11, "RAPID_TIER": 12, "HEADER_REGISTRY_MISMATCH": 13, "SIGMA_HIGH": 14,
    # reserved bits used by this implementation (documented in docs/decisions.md D-006)
    "NO_OTL": 15, "NO_ANTENNA_CALIBRATION": 16, "MF_FALLBACK": 17,
}

# ----------------------------------------------------------------------------------------------
# Preprocessing / QC [TDS § 13]
# ----------------------------------------------------------------------------------------------
MW_K = 4.0                           # k_MW
MW_MIN_CYC = 1.0                     # n_MW,min (cycles)
MW_SIGMA_FLOOR_CYC = 0.25
MW_CONFIRM_EPOCHS = 1
GF_POLY_DEGREE = 2
GF_POLY_EPOCHS = 10
GF_K = 4.0
GF_MIN_M = 0.05                      # for 30-s data [A]
ROTI_HIGH_TECU_MIN = 0.5             # [A]
ROTI_WINDOW_S = 300.0
CLOCK_JUMP_FRACTION = 0.8            # of satellites [A, TDS § 13.8]
CLOCK_JUMP_TOL_MS = 1e-4             # |k - round(k)| tolerance (ms) for integer-ms jump
SNR_MASK = False
RESET_AMBIGUITIES_AT_DAY_BOUNDARY = True   # CODE product-day continuity not verified -> reset [TDS § 11.2]
SATELLITE_BLACKLIST = []             # e.g. ["G04"]; excluded with a logged reason [TDS § 13.5]
USE_BROADCAST_HEALTH = True          # exclude epochs where the broadcast ephemeris flags the satellite unhealthy
SNR_MIN = {"1": 30.0, "2": 25.0}     # dB-Hz [A]

# Attitude / eclipse [TDS § 7.11]. Block-specific noon-turn limits derived from maximum hardware yaw
# rates (Kouba 2009; IIF: Dilssner 2010): beta0 = atan(mu/R), mu = 0.00836 deg/s [A, verify_log V-021].
ECLIPSE_BETA_NOON_DEG = {"BLOCK IIA": 4.0, "BLOCK IIR-A": 2.4, "BLOCK IIR-B": 2.4, "BLOCK IIR-M": 2.4,
                         "BLOCK IIF": 4.4, "BLOCK IIIA": 4.4, "DEFAULT": 4.5}
ECLIPSE_NOON_WINDOW_DEG = 10.0       # orbit angle from noon inside which the noon turn happens [A]
ECLIPSE_POST_SHADOW_S = 1800.0       # post-shadow recovery (IIA) [A, Kouba 2009]
ECLIPSE_SHADOW_BETA_DEG = 14.0       # shadow crossing possible only for |beta| < ~13.9 deg (GPS) [F]

# ----------------------------------------------------------------------------------------------
# Troposphere [TDS § 6]
# ----------------------------------------------------------------------------------------------
GRADIENT_MF_C = 0.0032               # Chen & Herring (1997) [F: IERS 2010 Ch. 9 eq. 9.12]
ZHD_SOURCE_PRIORITY = ["BAROMETER", "VMF3_GRID", "GPT3"]
MAPPING_FUNCTION = "VMF3"            # 'VMF3' (GPT3 fallback) or 'GMF' (benchmark only) [TDS F5]
ZWD_SCALE_HEIGHT_M = 2000.0          # Kouba (2008) wet-delay height reduction [F/L]

# ----------------------------------------------------------------------------------------------
# Coordinates [TDS § 9]
# ----------------------------------------------------------------------------------------------
SOI_FRAME = "ITRF2008"
SOI_EPOCH = 2005.0
# Transformation parameters FROM the key frame TO the target frame (IGN official files) [F].
# T mm, D ppb, R mas, rates per year, epoch decimal year. Convention X_S = X + T + D X + R X.
# Verified against PROJ data/ITRF2020 and data/ITRF2014 (copies of IGN files) - verify_log V-011.
ITRF_TRANSFORMS = {
    ("ITRF2020", "ITRF2008"): {"T": (0.2, 1.0, 3.3), "D": -0.29, "R": (0.0, 0.0, 0.0),
                               "dT": (0.0, -0.1, 0.1), "dD": 0.03, "dR": (0.0, 0.0, 0.0), "epoch": 2015.0,
                               "source": "IGN Transfo-ITRF2020_TRFs.txt"},
    ("ITRF2014", "ITRF2008"): {"T": (1.6, 1.9, 2.4), "D": -0.02, "R": (0.0, 0.0, 0.0),
                               "dT": (0.0, 0.0, -0.1), "dD": 0.03, "dR": (0.0, 0.0, 0.0), "epoch": 2010.0,
                               "source": "IGN Transfo-ITRF2014_ITRFs.txt"},
}
# Product frame label -> ITRF family [TDS § 9.2-9.3] (identity transformation within a family)
FRAME_LABELS = {
    "IGS08": "ITRF2008", "IGB08": "ITRF2008", "ITR08": "ITRF2008", "ITRF08": "ITRF2008",
    "IGS14": "ITRF2014", "IGB14": "ITRF2014", "ITR14": "ITRF2014", "ITRF14": "ITRF2014",
    "IGS20": "ITRF2020", "IGB20": "ITRF2020", "IGC20": "ITRF2020", "ITR20": "ITRF2020",
    "ITRF20": "ITRF2020",
}
# ITRF2020 Plate Motion Model, India plate [F: Altamimi et al. 2023, GRL, doi:10.1029/2023GL106373;
# ITRF2020-PMM.dat]. Rotation rates mas/yr; ORB translation rates mm/yr added (V = w x X + ORB).
PMM_INDIA = {"wx": 1.137, "wy": 0.013, "wz": 1.444, "orb": (0.37, 0.35, 0.74),
             "source": "ITRF2020-PMM (Altamimi et al. 2023) via PROJ data/ITRF2020 <INDI_T>"}
SIGMA_VU_MM_YR = 3.0                 # assumed vertical rate uncertainty [A, TDS § 9.4]
SIGMA_VH_MM_YR = 1.5                 # PMM horizontal residual uncertainty [L/E]
SIGMA_SOI_COORD_M = (0.005, 0.005, 0.010)   # assumed N, E, U sigma of SOI coordinates at 2005.0 [A]
# Deforming zones where the rigid-plate model is not valid [TDS § 9.11]; simple lat/lon boxes [A]
DEFORMING_ZONES = {
    "HIMALAYAN_FRONT": (27.0, 37.5, 72.0, 97.5),     # lat_min, lat_max, lon_min, lon_max
    "NORTH_EAST": (22.0, 29.5, 89.5, 97.5),
    "ANDAMAN_NICOBAR": (6.0, 14.5, 91.5, 94.5),
    "KUTCH": (22.5, 24.8, 68.0, 71.5),
}
# Auto-mode thresholds [TDS § 9.12]
AUTO_TOL_H_M = 0.02
AUTO_TOL_U_M = 0.02
AUTO_REFPOINT_TOL_M = 0.02
RADOME_FALLBACK_ZTD_SIGMA_M = 0.005  # [A, TDS § 10.6]
HEIGHT_ZTD_COUPLING_BETA = 0.3       # provisional beta until measured (§ 10.5) [L]

# Geoid for orthometric height: GPT3 grid undulation (EGM2008-based) [A, docs/decisions.md D-007]

# ----------------------------------------------------------------------------------------------
# PWV [TDS § 18]
# ----------------------------------------------------------------------------------------------
PWV_CONSTANTS = {  # Bevis et al. (1994) [F]
    "BEVIS1994": {"rho_w": 1000.0, "R_v": 461.5, "k2p": 0.221, "k3": 3739.0},   # SI: K/Pa, K^2/Pa
    "RUEGER2002": {"rho_w": 1000.0, "R_v": 461.51, "k2p": 0.2295, "k3": 3754.63},
}
PWV_CONSTANTS_SET = "BEVIS1994"
SIGMA_P_HPA = {"BAROMETER": 0.3, "ERA5": 1.0, "NWP": 1.5, "VMF3_GRID": 1.0, "GPT3": 5.0,
               "STANDARD_ATMOSPHERE": 10.0}
SIGMA_TM_K = {"ERA5": 2.0, "GPT3": 4.0, "BEVIS": 5.0}
# Station meteorology (RINEX met file next to SITE.o, found automatically; no extra CLI input, TDS § 33.3)
MET_FILE_PATTERNS = ["{st}{doy:03d}?.{yy:02d}[mM]", "{st}{doy:03d}?.{yy:02d}[mM].*", "{stl}{doy:03d}?.{yy:02d}[mM]*",
                     "{st}*_{y:04d}{doy:03d}0000_01D_*MM.rnx*"]
MET_P_RANGE_HPA = (500.0, 1100.0)    # plausible station pressure in India (sites below ~5.5 km) [A]
MET_SPIKE_HPA = 2.0                  # |P - 5-sample running median| rejection limit [A]
MET_MAX_GAP_S = 1800.0               # longer barometer gaps are filled from the next P source (flagged) [A]
MET_MAX_HEIGHT_DIFF_M = 200.0        # sensor-ARP height difference above which the barometer is not used [A]
MET_MAX_BIAS_HPA = 6.0               # |median(P_baro - P_reference)| above which the barometer is rejected [A]
# ERA5 for station PWV (TDS § 18.2: P priority 2, Tm priority 1). Used when cache/ERA5/era5_pl_<YYYYDDD>.nc exists
# (e.g. written by pwv_map.py); downloaded by the station run only if enabled (CDS queues can take hours, D-014).
STATION_ERA5 = True
STATION_ERA5_DOWNLOAD = False
S_ZTD_SCALE = 1.0                    # formal-error scale factor (to be derived, O14)
RD_DRY = 287.05                      # J/kg/K
G0 = 9.80665                         # m/s^2

# ----------------------------------------------------------------------------------------------
# System 2 mapping [TDS Part B]
# ----------------------------------------------------------------------------------------------
MAP_DOMAIN = (6.0, 38.0, 68.0, 98.0)
MAP_GRIDS = {"G050": 0.5, "G025": 0.25, "G0125": 0.125}
MAP_SLOT_MIN = 15
MAP_K_MAX = 40
MAP_L_BOUNDS_KM = (20.0, 500.0)
MAP_SIGMA_BG_REPR_MM = 0.5
MAP_DEFAULT_COV = {"sigma_s_mm": 2.0, "L_km": 150.0, "nugget_mm": 0.5}
MAP_CLASS_RULES = {"c1": (25.0, 3, 300.0, 0.7), "c2": (50.0, 2, 0.4), "c3": (150.0, 0.1)}
MAP_COASTAL_BUFFER_KM = 10.0

# ----------------------------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------------------------
CACHE_DIR = os.environ.get("PWV_PPP_CACHE", "cache")
RESULTS_DIR = "results"
BLQ_DIR = "blq"
STATIONS_FILES = ["stations.csv", "station_coord.csv"]
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def apply_overrides(path=None):
    """Override UPPER_CASE settings from a small YAML/JSON file (TDS § 33.2). Returns dict applied."""
    candidates = [path] if path else ["pwv_ppp_settings.yaml", "pwv_ppp_settings.yml", "pwv_ppp_settings.json"]
    g = globals()
    for p in candidates:
        if not p or not os.path.exists(p):
            continue
        with open(p) as fh:
            txt = fh.read()
        if p.endswith(".json"):
            data = json.loads(txt)
        else:
            try:
                import yaml  # optional
            except ImportError as exc:
                raise RuntimeError(f"{p}: PyYAML not installed; use a .json override file") from exc
            data = yaml.safe_load(txt) or {}
        applied = {}
        for k, v in data.items():
            if k.isupper() and k in g:
                g[k] = v
                applied[k] = v
        return applied
    return {}


def snapshot():
    """Canonical dict of all UPPER_CASE settings (for configuration hashing, TDS § 19.2)."""
    g = globals()
    out = {}
    for k in sorted(g):
        if k.isupper():
            v = g[k]
            if isinstance(v, dict):
                v = {str(kk): vv for kk, vv in v.items()}
            out[k] = v
    return out
