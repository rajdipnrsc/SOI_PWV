"""Station coordinates: SOI ITRF2008 @ 2005.0 -> processing frame at the observation epoch (TDS § 9).

Pipeline (TDS § 9.10):  official XYZ -> velocity model (SOI | PMM default) -> propagate in ITRF2008 (9.1)
-> table-driven frame transformation to the product's ITRF family (9.3) -> tide system (assumed
tide-free, § 9.11) -> reference point to ARP (RINEX DELTA H/E/N) -> immutable coordinate object.

Also: geodetic conversions (GRS80), ENU rotations, § 9.11 default assumptions, § 9.12 auto decision.
"""
import uuid
from dataclasses import dataclass, field

import numpy as np

from . import log as plog
from . import settings
from . import timesys as ts

LOG = plog.get()
MAS = np.pi / 180.0 / 3600.0 / 1000.0          # rad per mas (4.8481368e-9) [TDS § 9.5]


# ======================================================================================== geodesy
def ecef_to_geodetic(xyz):
    """GRS80 geodetic latitude, longitude (rad) and ellipsoidal height (m). Accepts (3,) or (n,3)."""
    xyz = np.asarray(xyz, dtype=float)
    x, y, z = xyz[..., 0], xyz[..., 1], xyz[..., 2]
    a = settings.RE_WGS84
    f = settings.F_GRS80
    e2 = f * (2 - f)
    lon = np.arctan2(y, x)
    p = np.hypot(x, y)
    lat = np.arctan2(z, p * (1 - e2))
    for _ in range(8):
        N = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
        h = p / np.cos(lat) - N
        lat = np.arctan2(z, p * (1 - e2 * N / (N + h)))
    N = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    h = p / np.cos(lat) - N
    return lat, lon, h


def geodetic_to_ecef(lat, lon, h):
    a = settings.RE_WGS84
    f = settings.F_GRS80
    e2 = f * (2 - f)
    N = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    return np.array([(N + h) * np.cos(lat) * np.cos(lon), (N + h) * np.cos(lat) * np.sin(lon),
                     (N * (1 - e2) + h) * np.sin(lat)])


def rot_xyz2enu(lat, lon):
    """Rotation matrix with rows = unit vectors E, N, U in ECEF: enu = R @ dxyz."""
    sl, cl, sp, cp = np.sin(lon), np.cos(lon), np.sin(lat), np.cos(lat)
    return np.array([[-sl, cl, 0.0], [-sp * cl, -sp * sl, cp], [cp * cl, cp * sl, sp]])


def xyz2enu(dxyz, lat, lon):
    return np.asarray(dxyz) @ rot_xyz2enu(lat, lon).T


def enu2xyz(denu, lat, lon):
    return np.asarray(denu) @ rot_xyz2enu(lat, lon)


# ============================================================================ frame transformation
def _helmert_params(key, t):
    p = settings.ITRF_TRANSFORMS[key]
    dt = t - p["epoch"]
    T = (np.array(p["T"]) + np.array(p["dT"]) * dt) * 1e-3
    D = (p["D"] + p["dD"] * dt) * 1e-9
    R = (np.array(p["R"]) + np.array(p["dR"]) * dt) * MAS
    return T, D, R


def _rmat(R):
    r1, r2, r3 = R
    return np.array([[0.0, -r3, r2], [r3, 0.0, -r1], [-r2, r1, 0.0]])


def helmert_forward(X, key, t):
    """X_S = X + T + D X + R X  (IGN convention, from key[0] to key[1]) at epoch t."""
    T, D, R = _helmert_params(key, t)
    return X + T + D * X + _rmat(R) @ X


def helmert_inverse(X_S, key, t):
    """First-order inverse (TDS eq. 9.3): X = X_S - T - D X_S - R X_S."""
    T, D, R = _helmert_params(key, t)
    return X_S - T - D * X_S - _rmat(R) @ X_S


def velocity_rate_terms(X, key):
    """d(X_S)/dt contribution from parameter rates: Tdot + Ddot X + Rdot X (m/yr)."""
    p = settings.ITRF_TRANSFORMS[key]
    return np.array(p["dT"]) * 1e-3 + p["dD"] * 1e-9 * X + _rmat(np.array(p["dR"]) * MAS) @ X


def frame_chain(source, target):
    """Table-driven chain source -> target (TDS § 9.5 'never composed ad hoc'). Returns list of steps."""
    if source == target:
        return []
    if (target, source) in settings.ITRF_TRANSFORMS:
        return [("inverse", (target, source))]
    if (source, target) in settings.ITRF_TRANSFORMS:
        return [("forward", (source, target))]
    raise KeyError(f"no transformation {source} -> {target} in settings.ITRF_TRANSFORMS")


def apply_chain(X, chain, t):
    for kind, key in chain:
        X = helmert_forward(X, key, t) if kind == "forward" else helmert_inverse(X, key, t)
    return X


def pmm_velocity_itrf2020(X, pmm=None):
    """ITRF2020 plate-motion-model velocity (m/yr): V = omega x X + ORB (TDS § 9.4, § 9.11)."""
    pmm = pmm or settings.PMM_INDIA
    w = np.array([pmm["wx"], pmm["wy"], pmm["wz"]]) * MAS
    return np.cross(w, X) + np.array(pmm["orb"]) * 1e-3


def velocity_itrf2008_from_pmm(X08):
    """PMM velocity expressed in ITRF2008: V08 = V20 + rate terms of the ITRF2020->ITRF2008 transform."""
    v20 = pmm_velocity_itrf2020(X08)
    return v20 + velocity_rate_terms(X08, ("ITRF2020", "ITRF2008"))


def label_to_family(label):
    lab = (label or "").strip().upper()
    return settings.FRAME_LABELS.get(lab)


# =========================================================================== coordinate object
@dataclass(frozen=True)
class Coordinate:
    X: float
    Y: float
    Z: float
    lat: float                   # deg
    lon: float                   # deg
    h_ell: float                 # m
    h_orth: float
    geoid_model: str
    frame: str
    realisation_label: str
    epoch: float                 # decimal year
    velocity: tuple
    velocity_frame: str
    velocity_source: str
    reference_epoch_of_source: float
    tide_system: str
    reference_point: str
    source_reference_point: str
    covariance: tuple            # 3x3 nested tuple (m^2)
    source: str                  # soi_transformed | ppp_same_day | ppp_validated | official_raw | rinex_approx
    transformation_history: tuple
    validation_status: str
    original_official_record: tuple
    id: str = field(default_factory=lambda: str(uuid.uuid4()))

    @property
    def xyz(self):
        return np.array([self.X, self.Y, self.Z])

    def as_dict(self):
        d = dict(self.__dict__)
        d["transformation_history"] = list(self.transformation_history)
        return d


def in_deforming_zone(lat_deg, lon_deg):
    for name, (a, b, c, d) in settings.DEFORMING_ZONES.items():
        if a <= lat_deg <= b and c <= lon_deg <= d:
            return name
    return None


def permanent_tide_offset(lat):
    """Mean-tide minus conventional tide-free position (radial, north) in m, IERS 2010 eq. 7.14 (hypothesis
    check only, TDS § 9.11)."""
    P2 = (3 * np.sin(lat) ** 2 - 1) / 2
    dr = (-0.1206 + 0.0001 * P2) * P2
    dn = (-0.0252 - 0.0001 * P2) * np.sin(2 * lat)
    return dr, dn


def soi_to_processing(xyz08, t_proc_ns, product_label, delta_hen=(0.0, 0.0, 0.0), refpoint="MARKER",
                      geoid_undulation=None, velocity08=None, velocity_source=None, sigma0_neu=None):
    """SOI ITRF2008@2005.0 coordinate -> ARP in the product frame at t (TDS § 9.5-9.10).

    Returns a Coordinate. Raises KeyError('FRAME_UNKNOWN') if the product label is not in the table.
    """
    X0 = np.asarray(xyz08, dtype=float)
    t0 = settings.SOI_EPOCH
    t = ts.decimal_year(t_proc_ns)
    family = label_to_family(product_label)
    if family is None:
        raise KeyError(f"FRAME_UNKNOWN: product frame label '{product_label}'")
    lat0, lon0, _ = ecef_to_geodetic(X0)
    history = [{"step": "official", "frame": settings.SOI_FRAME, "epoch": t0, "xyz": X0.tolist()}]
    if velocity08 is None:
        V08 = velocity_itrf2008_from_pmm(X0)
        vsrc = "ITRF2020-PMM India plate (+ORB), vertical 0"
        venu = xyz2enu(V08, lat0, lon0)
        venu[2] = 0.0                                   # vertical: no model -> 0 [TDS § 9.4]
        V08 = enu2xyz(venu, lat0, lon0)
    else:
        V08 = np.asarray(velocity08, dtype=float)
        vsrc = velocity_source or "user"
    X08t = X0 + V08 * (t - t0)                          # (9.1)
    history.append({"step": "propagate", "frame": settings.SOI_FRAME, "from_epoch": t0, "to_epoch": t,
                    "velocity_m_per_yr": V08.tolist(), "velocity_source": vsrc})
    chain = frame_chain(settings.SOI_FRAME, family)
    Xt = apply_chain(X08t, chain, t)                    # (9.3)
    for kind, key in chain:
        history.append({"step": "helmert_" + kind, "parameters": key[0] + "->" + key[1], "epoch": t,
                        "source": settings.ITRF_TRANSFORMS[key]["source"]})
    history.append({"step": "realisation", "label": product_label, "family": family,
                    "transformation": "identity (IGS realisation aligned to ITRF family)"})
    history.append({"step": "tide_system", "assumed": "TIDE_FREE", "conversion": "none"})
    lat, lon, h = ecef_to_geodetic(Xt)
    dH, dE, dN = delta_hen
    if refpoint == "MARKER":
        Xarp = Xt + enu2xyz(np.array([dE, dN, dH]), lat, lon)
        history.append({"step": "marker_to_arp", "delta_hen": list(delta_hen)})
    else:
        Xarp = Xt
        history.append({"step": "reference_point", "value": refpoint, "delta_hen_applied": False})
    lat, lon, h = ecef_to_geodetic(Xarp)
    # covariance (9.2)
    s0 = np.array(sigma0_neu or settings.SIGMA_SOI_COORD_M)
    sv = np.array([settings.SIGMA_VH_MM_YR, settings.SIGMA_VH_MM_YR, settings.SIGMA_VU_MM_YR]) * 1e-3
    s_enu = np.array([s0[1], s0[0], s0[2]]) ** 2 + ((t - t0) * np.array([sv[1], sv[0], sv[2]])) ** 2
    R = rot_xyz2enu(lat, lon)
    cov = R.T @ np.diag(s_enu) @ R
    N = geoid_undulation if geoid_undulation is not None else 0.0
    V20 = V08 - velocity_rate_terms(X08t, ("ITRF2020", "ITRF2008")) if family == "ITRF2020" else V08
    return Coordinate(
        X=float(Xarp[0]), Y=float(Xarp[1]), Z=float(Xarp[2]), lat=float(np.degrees(lat)), lon=float(np.degrees(lon)),
        h_ell=float(h), h_orth=float(h - N), geoid_model="GPT3 1-deg grid undulation" if geoid_undulation is not None
        else "none", frame=family, realisation_label=product_label, epoch=float(t),
        velocity=tuple(float(v) for v in V20), velocity_frame=family, velocity_source=vsrc,
        reference_epoch_of_source=t0, tide_system="TIDE_FREE", reference_point="ARP",
        source_reference_point=refpoint, covariance=tuple(tuple(float(c) for c in row) for row in cov),
        source="soi_transformed", transformation_history=tuple(history), validation_status="UNCHECKED",
        original_official_record=(settings.SOI_FRAME, t0, float(X0[0]), float(X0[1]), float(X0[2])))


def simple_coordinate(xyz, t_proc_ns, frame_label, source, sigma=1.0, geoid_undulation=None, history=None):
    """Coordinate object for an ARP estimated/approximated in the processing frame itself."""
    X = np.asarray(xyz, dtype=float)
    lat, lon, h = ecef_to_geodetic(X)
    N = geoid_undulation if geoid_undulation is not None else 0.0
    s = np.atleast_1d(sigma)
    cov = np.diag(np.broadcast_to(s ** 2, (3,)))
    fam = label_to_family(frame_label) or frame_label
    return Coordinate(
        X=float(X[0]), Y=float(X[1]), Z=float(X[2]), lat=float(np.degrees(lat)), lon=float(np.degrees(lon)),
        h_ell=float(h), h_orth=float(h - N), geoid_model="GPT3 1-deg grid undulation" if geoid_undulation is not None
        else "none", frame=fam, realisation_label=frame_label, epoch=float(ts.decimal_year(t_proc_ns)),
        velocity=(0.0, 0.0, 0.0), velocity_frame=fam, velocity_source="none (same-epoch estimate)",
        reference_epoch_of_source=float(ts.decimal_year(t_proc_ns)), tide_system="TIDE_FREE",
        reference_point="ARP", source_reference_point="ARP",
        covariance=tuple(tuple(float(c) for c in row) for row in cov), source=source,
        transformation_history=tuple(history or []), validation_status="UNCHECKED",
        original_official_record=())


# ============================================================================ § 9.12 auto decision
def discrepancy_enu(X_soi, X_A):
    lat, lon, _ = ecef_to_geodetic(X_A)
    return xyz2enu(np.asarray(X_soi) - np.asarray(X_A), lat, lon)


def auto_decision(X_soi_marker_based, X_soi_arp_based, X_A, delta_hen):
    """Apply the § 9.12 rules. Returns dict(decision, use, d_enu, refpoint_switched, message)."""
    dH = delta_hen[0]
    d = discrepancy_enu(X_soi_marker_based, X_A)
    out = {"refpoint_switched": False, "d_enu": d}
    if dH != 0.0:
        # data say the SOI coordinate is already the ARP if the vertical discrepancy equals dH
        d_arp = discrepancy_enu(X_soi_arp_based, X_A)
        if abs(abs(d[2]) - abs(dH)) <= settings.AUTO_REFPOINT_TOL_M and abs(d_arp[2]) < abs(d[2]):
            out["refpoint_switched"] = True
            d = d_arp
            out["d_enu"] = d
    ok = abs(d[0]) <= settings.AUTO_TOL_H_M and abs(d[1]) <= settings.AUTO_TOL_H_M and \
        abs(d[2]) <= settings.AUTO_TOL_U_M
    out["decision"] = "FIXED_SOI" if ok else "KEEP_STATIC"
    return out
