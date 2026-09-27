# GPS PPP → ZTD → PWV: Technical Design Specification (TDS)

| Field | Value |
|---|---|
| Document | `GPS_PPP_ZTD_PWV_Technical_Design_Specification.md` |
| Version | 1.1 (design baseline + Revision 1.1: Python-only, never-stop defaults, minimal CLI — see § 0.5) |
| Date | 27 Sep 2026 (v1.0: 25 Sep 2026) |
| Status | **Authoritative design contract** between the research/design phase and implementation |
| Baseline documents | `PPP_ZTD_Feasibility_Assessment.md` (main assessment; its §G product policy is superseded) and `PPP_Addendum_CODE_and_National_Mapping.md` (active product policy + national mapping) |
| Scope | System 1: GPS-only PPP → 5-min ZTD/ZWD/PWV at SOI CORS. System 2: network station products → India-wide PWV/ZWD/ZTD fields with uncertainty and confidence |
| Code | **None.** This document specifies; it does not implement |

---

## 0. How to read this document

### 0.1 Evidence and requirement labels

| Label | Meaning |
|---|---|
| **[F]** | Established fact from an authoritative specification or documentation |
| **[L]** | Literature finding; may not transfer to India without testing |
| **[A]** | Engineering assumption / design decision |
| **[E]** | Estimate or hypothesis that must be validated experimentally |
| **VERIFY DURING IMPLEMENTATION** | A detail (format label, constant, file path, external behaviour) that must be confirmed against the primary source before code depends on it. Implementation must not guess these |
| **SHALL / SHALL NOT** | Mandatory requirement |
| **SHOULD** | Strong recommendation; deviation requires a documented reason |
| **MAY** | Optional |

### 0.2 Four accuracy concepts that must never be merged

| Term | Definition in this project |
|---|---|
| **Theoretical capability** | What the method can achieve with ideal inputs, per literature |
| **Engineering target** | What we aim for during development |
| **Acceptance criterion** | The numerical gate a subsystem must pass to be declared complete |
| **Measured result** | What the validation experiment actually produced — the only number that may be published as performance |

### 0.3 Four spatial-resolution concepts that must never be merged

| Term | Definition |
|---|---|
| **Grid spacing** | The distance between output pixel centres (e.g., 0.25°). A property of the file |
| **Station density / spacing** | Local nearest-neighbour distance between valid stations (km). A property of the network at a time slot |
| **Effective atmospheric resolution** | The smallest horizontal scale of real atmospheric structure the product resolves at a location. Roughly 2 × local station spacing where GNSS dominates; the background model's effective resolution where it does not |
| **Supportable resolution** | The finest candidate grid spacing at which cross-validation demonstrates skill in a given region and weather regime (§ 24). Determined empirically, never assumed |

A 0.125° raster is **not** a 14-km atmospheric product unless the supportable-resolution test says so for that region and regime.

### 0.4 Two temporal-resolution concepts that must never be merged

| Term | Definition |
|---|---|
| **Sampling interval** | The spacing of output time stamps (station: 5 min; map: 15 min) |
| **Independent temporal information** | The time scale over which estimates are statistically independent. For random-walk-smoothed ZTD this is of order 15–30 min **[E]**; adjacent 5-min station values are strongly correlated and are **not** independent atmospheric observations |

### 0.5 Revision 1.1 — binding changes (27 Sep 2026)

These three decisions by the project owner **override any conflicting text elsewhere in this document**:

1. **Python only.** All of System 1 and System 2 are implemented in Python (NumPy/SciPy). There is no C++ core and no compiled extension. Speed-ups, if ever needed, use vectorised NumPy (optionally Numba), never a second language. (§ 33)
2. **Never stop on unknown station metadata.** The only known station fact is: **coordinates in ITRF2008 at epoch 2005.0**. Tide system, reference point, velocities, antenna-model era, ocean-loading coefficients and meteorological data are unknown. The software SHALL apply the documented default assumptions of § 9.11, print a WARNING for each, record it in the output flags and manifest, run the automatic consistency checks of § 9.12, and **continue**. Acceptance tests (§ 10.3) become automatic checks and optional diagnostics, **not blockers**.
3. **Minimal, simple program.** One command runs everything: `python pwv_ppp.py SITE.o SITE.n` (optional flags in § 33.3). The program downloads every product it needs by itself, caches it, and shows progress logs on screen (and in a log file). The code structure is deliberately flat and simple (§ 33.2): scientific rigour lives in the equations, not in software complexity.

The only situations in which a run ends without ZTD output are those in which PPP is physically impossible: no usable dual-frequency GPS data in the RINEX file, or no CODE precise orbit/clock products available for the date (status `WAIT_FOR_TIER` or `FAIL`, § 11.4). Even then the program exits cleanly with a clear message, never with a traceback.

---

## 1. System overview

```
                         ┌──────────────────────────── SYSTEM 1 ───────────────────────────┐
 SOI RINEX 3 ──► Station │ Coordinate  ─► Product  ─► Preproc/QC ─► Float PPP  ─► 5-min    │──► Station ZTD product
 station logs    registry│ reference      manager     (slips,       KF + RTS      extraction│    (CSV/NetCDF/SINEX_TRO) 
 SOI XYZ (ITRF2008@2005) │ module         (CODE-only)  outliers)    [PPP-AR later]          │
                         │                                                                  │──► PWV module ──► Station PWV product
                         └──────────────────────────────────────────────────────────────────┘
                                                            │  Interface contract (§ 20)
                         ┌──────────────────────────── SYSTEM 2 ───────────────────────────┐
 Station products ─────► │ Ingest + QC ─► Temporal windows ─► NWP background at station     │
 NWP / ERA5 profiles     │ height ─► residuals ─► regression kriging ─► grid residual        │──► National grid products
 DEM, geoid, masks       │ ─► + background at cell height ─► PWV/ZWD/ZTD + σ + confidence    │    (NetCDF archive, Zarr, COG)
                         └──────────────────────────────────────────────────────────────────┘
```

Design principles **[A]**:
1. **Troposphere first.** Every choice is judged by its effect on ZTD/PWV, not on coordinates.
2. **One Analysis Centre (CODE), one solution family per processing window.**
3. **No silent substitution, no silent tier change, no silent coordinate approximation.**
4. **Float PPP is the validated baseline.** PPP-AR is a separate, later layer that can never overwrite float products.
5. **Every number is traceable** to inputs, configuration and software version.
6. **System 2 consumes only the published System 1 interface**, never System 1 internals.

---

## 2. Project facts and inputs

| Item | Status |
|---|---|
| Stations | SOI CORS, static; SOI reports > 1,000 (announced > 1,105) **[F]** |
| Data | RINEX 3.x observation files, GPS dual-frequency code + phase (sampling to be confirmed per station) |
| **Official SOI coordinates** | **ITRF2008 at reference epoch 2005.0** **[F — confirmed by project owner]** |
| SOI velocities | **Unknown** → default: ITRF2020 plate-motion model horizontal, zero vertical (§ 9.11) |
| SOI tide system, reference point (marker/ARP/APC), ANTEX used to derive coordinates | **Unknown** → documented default assumptions + automatic checks; processing never stops (§ 9.11–9.12) |
| Implementation | **Python only**; one command `python pwv_ppp.py SITE.o SITE.n` (§ 33) |
| Product policy | CODE/AIUB only for Analysis-Centre products (`COD0OPSFIN` → `CODMOPSRAP` → `COD0OPSRAP`) **[A — project decision]** |
| Constellation (v1) | GPS only |
| Mapping target | 0.25° national, 0.125° regional (if supported), 0.5° coarse; 15-min maps; hourly archive — **initial design targets, not proven resolutions** |

# PART A — SYSTEM 1: STATION PPP → ZTD

## 3. Constants, signals and notation

| Symbol | Value | Source |
|---|---|---|
| c | 299 792 458 m s⁻¹ | [F] SI |
| f₁ (GPS L1) | 1575.42 MHz | [F] IS-GPS-200 |
| f₂ (GPS L2) | 1227.60 MHz | [F] IS-GPS-200 |
| λᵢ = c/fᵢ | λ₁ ≈ 0.1903 m, λ₂ ≈ 0.2442 m | derived |
| α₁ = f₁²/(f₁²−f₂²) | ≈ 2.5457 | derived |
| α₂ = f₂²/(f₁²−f₂²) | ≈ 1.5457 (α₁ − α₂ = 1) | derived |
| λ_WL = c/(f₁−f₂) | ≈ 0.8619 m | derived |
| λ_NL = c/(f₁+f₂) | ≈ 0.1070 m | derived |
| ω_E | 7.2921151467 × 10⁻⁵ rad s⁻¹ | [F] IS-GPS-200 / WGS84 |
| GM | 3.986004418 × 10¹⁴ m³ s⁻² | [F] IERS 2010 / WGS84 |

Notation: superscript *s* = satellite, subscript *r* = receiver, *k* = epoch index, *e* = elevation, *a* = azimuth (from north, clockwise), *z* = zenith angle = 90° − e. All internal quantities in **SI units (m, s, rad)**; output units are specified in § 17.

Frequency-dependent values are always taken from the signal definition table (system, band, RINEX code), never hard-coded in the estimator, so that later systems (Galileo E1/E5a, etc.) plug in without estimator changes **[A]**.

---

## 4. Observation equations

### 4.1 Raw observables (per satellite s, epoch k)

```
Pᵢ = ρ + c(dt_r − dt^s) + T + Iᵢ + b_{r,Pᵢ} − b^s_{Pᵢ} + ε_{Pᵢ}                       [m]
Lᵢ = λᵢφᵢ = ρ + c(dt_r − dt^s) + T − Iᵢ + λᵢ(Nᵢ + β_{r,Lᵢ} − β^s_{Lᵢ}) + λᵢ w + ε_{Lᵢ}   [m]
```

| Term | Definition | Unit |
|---|---|---|
| Pᵢ | Pseudorange on band i, from RINEX code observable (e.g., C1W) | m |
| Lᵢ | Carrier phase converted to metres: λᵢ × (RINEX phase in cycles) | m |
| ρ | Geometric distance between satellite antenna phase centre at transmit time and receiver antenna phase centre at receive time, in the processing frame (inertial-effects handled via Sagnac rotation, § 7.3) | m |
| dt_r | Receiver clock error w.r.t. GPS time | s |
| dt^s | Satellite clock error w.r.t. GPS time (from CODE CLK, plus relativistic periodic correction § 7.4) | s |
| T | Slant tropospheric delay (§ 6) | m |
| Iᵢ | First-order ionospheric delay on band i (code delay, phase advance) | m |
| b, β | Code / phase hardware biases (receiver r, satellite s) | m / cycles |
| Nᵢ | Integer carrier-phase ambiguity | cycles |
| w | Phase wind-up (§ 7.9) | cycles |
| ε | Noise + multipath + unmodelled error | m |

### 4.2 Ionosphere-free (IF) combination (v1 formulation) [F/A]

```
P_IF = α₁P₁ − α₂P₂
L_IF = α₁L₁ − α₂L₂
```

After applying all deterministic corrections (§ 7), including satellite OSBs (§ 11.5):

```
P_IF = ρ̃ + c·dt_r − c·dt^s + T + ε_P,IF                                            (4.1)
L_IF = ρ̃ + c·dt_r − c·dt^s + T + B_IF + λ_NL·w + ε_L,IF                            (4.2)
```

- ρ̃ = ρ plus all deterministic range corrections (antenna PCO/PCV, tides, loading, Shapiro, etc.), § 7.
- **Receiver code bias** of the IF code is absorbed into c·dt_r (the receiver clock is defined by the IF code) **[F]**.
- **B_IF** (m) = IF ambiguity including receiver and uncalibrated phase biases. Float: real-valued, constant per arc.
- Phase wind-up in the IF combination is λ_NL·w (since w is equal in cycles on both bands) **[F, derivation: α₁λ₁ − α₂λ₂ = c/(f₁+f₂)]**.
- **Higher-order ionosphere** is not modelled in v1 (§ 7.17).
- Noise: σ²_IF = α₁²σ₁² + α₂²σ₂² ≈ (2.98 σ_raw)² for equal raw noise **[F, derivation]**.

Relation of B_IF to integer ambiguities (used in PPP-AR, § 15):

```
B_IF = λ_NL·N₁ + (λ_WL·f₂/(f₁+f₂))·N_WL + (bias terms),    N_WL = N₁ − N₂        (4.3)
```

### 4.3 Classification of quantities

| Category | Quantities |
|---|---|
| **Directly observed** | P₁, P₂, L₁, L₂ (RINEX), SNR, LLI flags, receiver epoch tags |
| **Deterministically computed** | Satellite position/velocity (SP3), satellite clock (CLK), relativistic clock term, Shapiro delay, Sagnac rotation, satellite and receiver PCO/PCV, phase wind-up, solid Earth tide, ocean tide loading, pole tide, ARP eccentricity, station coordinates (fixed mode), a priori ZHD, a priori ZWD, mapping functions, OSBs, Earth orientation (for Sun/Moon/pole tide) |
| **Stochastic parameters (states)** | ZWD (random walk), G_N, G_E (random walk), receiver clock (white noise) |
| **Constants (states)** | IF ambiguities (per arc), station coordinates (static_estimated / constrained modes) |
| **Nuisance parameters** | Receiver clock, ambiguities, (coordinates in constrained/static mode) — estimated but not products of interest |
| **Primary products** | ZTD = ZHD₀ + ZWD̂ (smoothed), σ_ZTD, G_N, G_E |

---

## 5. Estimator specification

### 5.1 Parameter registry (multi-GNSS-ready) [A]

Every state has a unique typed key (e.g., `ZWD`, `GRAD_N`, `GRAD_E`, `CLK_REC[GPS]`, `AMB[G05,arc=17]`, `POS_X`…). The filter operates on the registry, not on fixed indices. Adding Galileo later adds `CLK_REC[GAL]` (or an inter-system bias) and `AMB[E..]` keys without estimator changes.

### 5.2 State vectors

**Fixed mode (`fixed`)**

```
x_k = [ ZWD, G_N, G_E, c·dt_r, B_IF¹, …, B_IF^{n_k} ]ᵀ          dimension 4 + n_k
```

Coordinates are not states; ρ̃ is computed from the processing coordinate (§ 9).

**Static estimated mode (`static_estimated`)**

```
x_k = [ X, Y, Z, ZWD, G_N, G_E, c·dt_r, B_IF¹, …, B_IF^{n_k} ]ᵀ  dimension 7 + n_k
```

X, Y, Z = ARP position in the processing frame at the processing epoch (§ 9.5), constant over the processing window.

**Constrained mode (`constrained`)**

Same state as static_estimated. The coordinate state is initialised at the processing coordinate X₀ with covariance

```
P₀,xyz = R_enu→xyz · diag(σ_N², σ_E², σ_U²) · R_enu→xyzᵀ
```

and zero process noise. Because the coordinate is constant, this is mathematically identical to one pseudo-observation X = X₀ with that covariance **[F]**. Defaults: σ_N = σ_E = 3 mm, σ_U = 5 mm **[A]** (configurable).

Diagnostic (all coordinate modes with states): normalised coordinate pull `z_c = (X̂ − X₀)_enu / σ_prior,enu` from the smoothed solution. |z_c,U| > 3 → QC flag `COORD_PULL`.

### 5.3 Transition model (epoch k → k+1, Δt = t_{k+1} − t_k in s)

| State | Φ | Q (process noise) | Notes |
|---|---|---|---|
| ZWD | 1 | q_ZWD · Δt, q_ZWD = (σ_ZWD)² / 3600 | σ_ZWD in m/√h; default 0.006 (tuning § 6.7) |
| G_N, G_E | 1 | q_G · Δt, q_G = (σ_G)² / 3600 | σ_G default 0.0005 m/√h [A], tuning § 6.7 |
| c·dt_r | **0** | σ²_clk,0 | White noise: fully re-estimated each epoch; σ_clk,0 = 100 m [A]; a priori value = robust (median) code-derived clock at k+1 |
| B_IF (existing arc) | 1 | 0 | Constant |
| B_IF (new arc) | augmentation | prior σ²_amb,0 | § 5.8 |
| B_IF (ended arc) | removal (selection) | — | § 5.8 |
| X, Y, Z | I | 0 | Static/constrained |

Φ is a generalised transition matrix (identity for kept states, zero row for clock, selection for removed states, augmentation for new states). The RTS smoother (§ 5.10) is defined for such dimension-changing transitions provided P_{k+1|k} is non-singular, which the priors above guarantee **[A — verify numerically in Level-0 tests]**.

### 5.4 Measurement model and covariance

For each valid satellite s at epoch k, two rows (IF code, IF phase):

```
Partial derivatives (fixed mode):
  ∂/∂ZWD   = m_w(e)
  ∂/∂G_N   = m_g(e)·cos(a)
  ∂/∂G_E   = m_g(e)·sin(a)
  ∂/∂(c·dt_r) = 1
  ∂/∂B_IF^s  = 1 (phase row only)
  ∂/∂X     = −(r^s − r_r)ᵀ/ρ (unit line-of-sight, static/constrained modes)
Computed observation:
  h = ρ̃ − c·dt^s + m_h(e)·ZHD₀ + m_w(e)·ZWD + m_g(e)(G_N cos a + G_E sin a) + c·dt_r [+ B_IF + λ_NL·w]
```

Measurement covariance R_k: diagonal (code and phase uncorrelated; satellites uncorrelated) **[A]**:

```
σ²_obs(e) = (σ_a² + σ_b²/sin²e) · F_IF² / w_robust
```

| Parameter | Default [A] | Range to test |
|---|---|---|
| σ_a, σ_b (raw phase) | 0.003 m, 0.003 m | 0.002–0.005 |
| σ_a, σ_b (raw code) | 0.3 m, 0.3 m | 0.2–1.0 |
| F_IF | √(α₁² + α₂²) ≈ 2.98 | derived |
| w_robust | 1 (IGG-III factor, § 13.10) | — |
| Optional SNR model | σ² ∝ 10^(−SNR/10) (Hartinger & Brunner 1999) | option, off by default |

### 5.5 Initialisation (first epoch of a processing window)

| State | Initial value | Initial σ |
|---|---|---|
| ZWD | a priori ZWD₀ from the same source as ZHD₀ (VMF3 grid ZWD reduced to station height, or GPT3) | 0.10 m [A] |
| G_N, G_E | 0 | 0.001 m [A] |
| c·dt_r | robust code-derived value | 100 m |
| B_IF | L_IF − P_IF at arc start (mean over first valid epoch) | 20 m [A] |
| X, Y, Z (static) | processing coordinate X₀ | 1 m per axis [A] |
| X, Y, Z (constrained) | processing coordinate X₀ | σ_N, σ_E, σ_U (§ 5.2) |

### 5.6 Update algorithm

- Linearised EKF about the current state (the model is nearly linear once geometry is computed; one relinearisation per epoch is sufficient for fixed mode; for static mode, re-run the window once with the converged coordinate as new linearisation point) **[A]**.
- Covariance update in **Joseph form**; symmetrise after each update; check positive definiteness (Cholesky) and log failures **[A]**. A square-root (UD / SRIF) implementation MAY replace this if Level-0 numerical tests show loss of precision.
- Observations processed as a vector update per epoch (not sequential scalar updates) to allow innovation-based outlier testing on the full set.

### 5.7 State reset conditions

| Condition | Action |
|---|---|
| Cycle slip detected (§ 13.6–13.7) | New ambiguity arc for that satellite |
| LLI bit 0 set on L1 or L2 | New arc |
| Data gap for satellite > `gap_reset_s` (default 120 s) [A] | New arc |
| Observable code change within an arc (e.g., L2W → L2L) | New arc (arcs SHALL use one observable pair throughout) |
| Satellite below cutoff then re-rises | New arc |
| ≥ `n_consec_reject` (default 3) consecutive rejected phase observations | New arc |
| Unrepaired receiver clock jump (§ 13.8) | All arcs reset |
| Station-wide gap > `gap_reset_all_s` (default 600 s) [A] | All arcs reset; ZWD and gradients kept with process noise propagated over the gap (random walk grows with √Δt) |
| Gap > `window_break_s` (default 3600 s) [A] | Processing window split; smoother run separately per segment; 5-min outputs inside the gap = missing |

### 5.8 Ambiguity birth / death

- **Birth** at the first epoch of an arc: augment state with B_IF, prior value § 5.5, prior σ = σ_amb,0; zero cross-covariance with existing states.
- **Death** at the last epoch of an arc: in the **float pass**, the ambiguity is removed from the state after its last epoch (selection transition), keeping the smoother memory small.
- **AR pass** (§ 15): a separate forward pass retains all ambiguities to the end of the window (Φ = 1, Q = 0 after death), so that the final covariance contains the full joint ambiguity covariance. For constant parameters, the final forward estimate equals the full-window estimate **[F]**.
- **Short arcs:** arcs shorter than `min_arc_s` (default 600 s; PRIDE's strict mode excludes arcs < 10 min [F]) are excluded before estimation.
- **Satellite reappearance** is always a new arc (no bridging in v1).

### 5.9 Forward pass outputs to store (per epoch)

x̂_{k|k}, P_{k|k}, x̂_{k+1|k}, P_{k+1|k}, Φ_k (sparse), the state-key list, innovations, their covariance, rejected-observation list.

### 5.10 RTS backward smoother

```
C_k   = P_{k|k} Φ_kᵀ P_{k+1|k}⁻¹
x̂_{k|N} = x̂_{k|k} + C_k (x̂_{k+1|N} − x̂_{k+1|k})
P_{k|N} = P_{k|k} + C_k (P_{k+1|N} − P_{k+1|k}) C_kᵀ
```

- Computed with Cholesky solves (no explicit inverse).
- The smoothed solution is the **primary product**. Forward-only values are retained as a diagnostic (convergence analysis) and for later near-real-time use.
- Processing windows SHOULD overlap day boundaries: default window = [D 00:00 − 3 h, D+1 00:00 + 3 h] (30 h) when adjacent-day RINEX and products of the same family/tier exist; only outputs for day D are published **[A]**. Without adjacent data, publish with edge flags (§ 6.6).


## 6. Troposphere estimation specification

### 6.1 Model

```
T(e, a) = m_h(e)·ZHD₀ + m_w(e)·ZWD + m_g(e)·(G_N·cos a + G_E·sin a)            (6.1)
ZTD = ZHD₀ + ZWD                                                              (6.2)
```

| Item | Frozen choice |
|---|---|
| Estimated tropospheric state | **ZWD** (absolute wet zenith delay, initialised at ZWD₀) — not ZTD. ZHD₀ is held fixed |
| Primary product | **ZTD** = ZHD₀ + ZWD̂ (smoothed). ZTD is insensitive to the ZHD₀ choice to first order because m_h ≈ m_w; the estimated ZWD absorbs ZHD₀ error **[F]** |
| ZHD/ZWD in System 1 output | The *a priori* split: ZHD = ZHD₀ (with its source recorded), ZWD = ZTD − ZHD₀. **The PWV module re-splits ZTD with its own best pressure** (§ 18); System 1's ZWD is not the scientific ZWD unless ZHD₀ came from an accurate barometer |
| Hydrostatic MF m_h | **VMF3** (Landskron & Böhm 2018) with gridded a_h from the TU Wien VMF data server; b_h, c_h from the VMF3 empirical formulation; station-height correction per the VMF3 reference implementation — VERIFY DURING IMPLEMENTATION (exact height-correction formula and grid interpolation as in the official VMF3 code) |
| Wet MF m_w | VMF3 (a_w from grid, b_w, c_w per VMF3) |
| Fallback MF | GPT3 (empirical VMF3 coefficients from the GPT3 grid). Fallback use is logged; VMF1/GMF selectable only for benchmarking against PRIDE |
| Gradient MF m_g | Chen & Herring (1997): m_g(e) = 1 / (sin e · tan e + C), C = 0.0032 — VERIFY DURING IMPLEMENTATION (constant as given in IERS Conventions 2010 Ch. 9) |
| Gradient units | G_N, G_E in metres (reported in mm); delay at elevation e, azimuth a as in (6.1) |
| Elevation cutoff | 7° default [A; PRIDE default is 7° F]; tested at 7/10/15° |

### 6.2 A priori ZHD₀ (priority order, configurable)

1. **Station barometer** at known height: pressure reduced to ARP orthometric height (§ 18.3), then Saastamoinen (Davis et al. 1985):
   ```
   ZHD₀ = 0.0022768 · P / (1 − 0.00266·cos 2φ − 0.28×10⁻⁶·H)       [m; P in hPa; φ latitude; H orthometric height in m]
   ```
2. **VMF3 grid ZHD** reduced from grid height to station height (Kouba 2008, J. Geod. 82 — VERIFY DURING IMPLEMENTATION the height-reduction equations for gridded ZHD/ZWD) **(default)**.
3. **GPT3** pressure → Saastamoinen (flagged `ZHD_CLIMATOLOGY`; acceptable for ZTD, not for PWV).

### 6.3 A priori ZWD₀

Same source as ZHD₀ where available (VMF3 grid ZWD with height reduction), otherwise GPT3 + Askne & Nordius (1987). ZWD₀ affects only initialisation.

### 6.4 Stochastic model

- ZWD: **random walk at the observation epoch rate**, σ_ZWD default **6 mm/√h** (in the recommended 5–10 mm/√h band) [A]; this value is **not assumed universally correct** and is set by the tuning experiment (§ 6.7), possibly per season/region.
- Gradients: random walk, σ_G default 0.5 mm/√h [A], tuned (§ 6.7). Gradient estimation ON by default in all modes.
- Initial variances: § 5.5.

### 6.5 Smoothing and 5-minute extraction

1. Forward EKF over the processing window (§ 5.6).
2. RTS backward pass (§ 5.10) → x̂_{k|N}, P_{k|N} at every observation epoch.
3. **Output epochs**: exact UTC multiples of 5 min (hh:00, hh:05, …) **[A]**. Observation epochs are GPST-aligned (GPST − UTC = leap-second offset, 18 s since 2017 [F]), so each UTC mark lies between two observation epochs.
4. Value at output epoch t*: **linear interpolation** of smoothed ZWD, G_N, G_E between the bracketing epochs (Δ ≤ sampling interval, e.g., 30 s); σ = the larger of the two bracketing smoothed σ values (conservative) **[A]**. Receiver clock output = value at the nearest epoch.
5. If either bracketing epoch is missing or lies in a gap: no value; row written with `QC_flag` bit `NO_DATA`.
6. ZTD = ZHD₀(t*) + ZWD̂(t*); σ_ZTD = σ_ZWD (ZHD₀ treated as fixed) — the ZHD₀ error is carried in the PWV module, not here.

**Statement for all product documentation (mandatory):** *Adjacent 5-minute values are strongly correlated through the random-walk model and smoothing. They are samples of a smoothed estimate, not independent 5-minute atmospheric observations. The effective independent time resolution is estimated at 15–30 min [E] and must be quantified in validation.*

### 6.6 Convergence and edge flags

| Flag | Condition |
|---|---|
| `EDGE` | Output epoch within `edge_window_s` (default 3600 s) of the start or end of a continuous data segment where no overlapping data exist beyond it |
| `REINIT` | More than 50 % of active ambiguities born within the previous 1800 s |
| `CONVERGED` | Otherwise, and smoothed σ_ZTD ≤ `sigma_conv` (default 10 mm) |
| `NOT_CONVERGED` | Smoothed σ_ZTD > `sigma_conv` |

Forward-only convergence time (diagnostic, § 30): first epoch after which |ZTD_forward − ZTD_smoothed| < 10 mm for ≥ 30 min.

### 6.7 Formal tuning experiment (ZWD and gradient process noise)

| Element | Specification |
|---|---|
| Dataset | Stage 0/2 pilot stations; later the 10–15 station year set; stratified by season (winter, pre-monsoon, SW monsoon, post-monsoon) |
| Values | σ_ZWD ∈ {3, 6, 10, 15} mm/√h (plus PRIDE default 24 mm/√h for reference); σ_G ∈ {0.1, 0.3, 0.5, 1.0} mm/√h |
| Fixed | All else per frozen configuration |
| Metrics | (1) Hourly and 15-min RMSE/bias vs independent references (IGS troposphere product at IGS stations, radiosonde, MWR where available; ERA5 as secondary); (2) normalised innovation squared (NIS) of phase residuals — mean NIS/dof within [0.8, 1.2]; (3) high-frequency power of ZTD (structure function at 5–30 min lags) vs radiosonde/MWR where available; (4) phase post-fit residual RMS by elevation |
| Decision rule | Choose the smallest σ_ZWD whose RMSE vs independent references is within 5 % of the minimum, subject to NIS criterion; separate values per season allowed only if the difference is significant (block bootstrap over days, 95 %) |
| Output | Frozen value(s) recorded in configuration; report with plots of RMSE vs σ_ZWD per season |


## 7. Deterministic correction specification

Every correction is implemented as a separately testable function returning its own term, so that the modelled range can be decomposed term by term (Level-0 tests, § 30.1). Sign convention for this section: **a positive range correction increases the modelled observation** (it is added to ρ̃ in (4.1)/(4.2)).

Authorities: IERS Conventions 2010 (IERS TN 36) and its official updates **[F]**; IS-GPS-200; ANTEX 1.4; RINEX 3.05/4.0x; SP3-d; Clock RINEX 3.04; Bias-SINEX 1.00; ORBEX 0.09. Where the IERS Conventions have been updated after 2010 (e.g., the secular pole model for pole tide), use the current official version — VERIFY DURING IMPLEMENTATION which chapter updates apply.

### 7.1 Time systems and epochs

| Aspect | Specification |
|---|---|
| Formulation | GPST = TAI − 19 s; UTC = GPST − ΔLS(t) (ΔLS = 18 s since 2017-01-01 [F]); TT = TAI + 32.184 s; UT1 = UTC + (UT1−UTC) from ERP |
| Inputs | Leap-second table (IERS Bulletin C; maintained as data file with version), CODE ERP |
| Units | Internal time: integer GPS week + seconds of week (float64), or integer nanoseconds since GPS epoch **[A]** — never float64 seconds since 1980 for sub-ns work |
| Reception time | t_rx = epoch_tag − dt_r (if RINEX `RCV CLOCK OFFS APPL` = 0); if = 1, the receiver has applied its own offset — handle per RINEX 3.05 definition, VERIFY |
| Transmit time | t_tx = t_rx − τ, τ = ρ/c, iterated until |Δτ| < 10⁻¹² s |
| Mandatory v1 | Yes |
| Validation | Known conversions (tabulated test vectors across leap-second boundaries and GPS week rollovers) |
| Common errors | Forgetting leap seconds in UTC output; using epoch tag directly as reception time; float64 precision loss; mixing TT/UTC when computing Sun/Moon positions |

### 7.2 Precise satellite orbit (SP3)

| Aspect | Specification |
|---|---|
| Formulation | Lagrange interpolation, 11 nodes (degree 10) centred on t_tx, on CODE 5-min SP3 positions (Earth-fixed, satellite centre of mass); velocity from analytic derivative of the interpolating polynomial |
| Inputs | CODE SP3 for D−1, D, D+1 (same family/tier) |
| Frame | SP3 coordinate system from header (§ 9.3) — ECEF of the product's reference-frame realisation |
| Units | SP3 km → m |
| Application | r^s(t_tx) centre of mass, then satellite PCO (§ 7.6) |
| Magnitude | — |
| Mandatory | Yes |
| Validation | Leave-one-node-out: interpolate a withheld SP3 node → error < 1 mm RMS; continuity across day boundary < 5 mm (report) |
| Common errors | Edge extrapolation without adjacent days; using satellites with bad values (999999.999999 / 0.000000 per SP3-d); ignoring SP3 clock/position event flags; mixing tiers across days |

### 7.3 Earth rotation during signal travel (Sagnac)

| Aspect | Specification |
|---|---|
| Formulation | r^s_rx-frame = R₃(ω_E·τ)·r^s(t_tx), with R₃(θ) the rotation about Z by angle θ such that the satellite position is expressed in the ECEF frame at reception time |
| Inputs | τ, ω_E |
| Magnitude | Up to ~30 m in range [L] |
| Mandatory | Yes |
| Validation | Compare with RTKLIB geometric range function (known implementation) < 0.1 mm |
| Common errors | Wrong rotation sign; applying both the rotation and a separate Sagnac range term (double counting); not iterating τ |

### 7.4 Relativistic satellite clock correction

| Aspect | Specification |
|---|---|
| Formulation | dt^s = dt^s_CLK − 2·(r^s · v^s)/c²  (s) |
| Reference | IERS 2010 Ch. 10; IGS clock convention (products exclude the periodic term) — VERIFY in CODE clock header/README |
| Inputs | Interpolated r^s, v^s (ECEF; the ECEF-velocity approximation is standard practice) |
| Magnitude | Up to ~±45 ns (~±14 m) for GPS eccentricities [L] |
| Mandatory | Yes |
| Validation | Against RTKLIB/independent computation < 1 mm; sign check: residuals must not show a once-per-revolution sinusoid |
| Common errors | Sign error; applying the broadcast-style F·e·√A·sinE term *and* the r·v term |

### 7.5 Gravitational (Shapiro) delay

| Aspect | Specification |
|---|---|
| Formulation | Δρ_rel = (2GM/c²) · ln[(r^s + r_r + ρ)/(r^s + r_r − ρ)] (m), r^s, r_r geocentric distances |
| Reference | IERS 2010 Ch. 11 |
| Magnitude | ~18–19 mm, elevation-dependent [L] |
| Mandatory | Yes (elevation-dependent → otherwise aliases into ZTD) |
| Validation | Independent computation < 0.01 mm |
| Common errors | Omission; using ρ at wrong epoch (negligible) |

### 7.6 Satellite antenna PCO and PCV

| Aspect | Specification |
|---|---|
| Formulation | r^s_APC = r^s_CoM + R_body→ECEF · PCO_IF;  PCO_IF = α₁PCO_L1 − α₂PCO_L2 (body frame x,y,z). PCV: nadir-angle (and azimuth if given) dependent; PCV_IF = α₁PCV_L1 − α₂PCV_L2, applied as range correction |
| Reference | ANTEX 1.4; IGS igs20.atx (frame-consistent ANTEX, § 9.3) |
| Inputs | ANTEX (satellite by SVN/PRN validity period), satellite attitude (§ 7.11) |
| Frame | Satellite body frame (IGS convention: z toward Earth, y along solar panel axis, x completing, per IGS attitude convention) — VERIFY against ANTEX/ORBEX definitions |
| Magnitude | PCO z ≈ 0.5–2.7 m depending on block [L]; PCV mm-level |
| Mandatory | Yes |
| Validation | Compare APC positions with an independent implementation (RTKLIB) < 0.5 mm; check PRN→SVN mapping on dates with PRN reassignment |
| Common errors | Using PRN instead of SVN for antenna lookup; wrong ANTEX (e.g., igs14 with IGS20 products); applying PCO to Melbourne–Wübbena in AR (CODE says do not, § 15); sign of PCV; IF combination with wrong frequencies |

### 7.7 Receiver antenna PCO and PCV

| Aspect | Specification |
|---|---|
| Formulation | Range correction Δρ_ant = −(PCO_IF,enu · ê_enu) + PCV_IF(z, a), where ê is the unit vector receiver→satellite in local ENU, PCO relative to ARP (ANTEX gives N, E, U) |
| Reference | ANTEX 1.4 — VERIFY DURING IMPLEMENTATION the sign convention of PCV application against the ANTEX 1.4 text and a reference implementation (RTKLIB `antmodel`) using a test antenna with known asymmetric PCV |
| Inputs | Antenna type + radome (20-character ANTEX code) from RINEX header `ANT # / TYPE`, cross-checked with station log (§ 8) |
| Frame | Local ENU at station; antenna assumed north-oriented (RINEX `ANTENNA: ZERODIR` if present overrides) |
| Magnitude | PCO up ~5–15 cm; PCV up to ~1–2 cm at low elevation [L] |
| Mandatory | Yes — **one of the largest ZTD bias sources** |
| Validation | Round-trip test; comparison with RTKLIB; ZTD sensitivity test with deliberately wrong radome |
| Common errors | Radome ignored (using `NONE` calibration silently); PCO applied in XYZ instead of ENU; N/E/U order swapped; azimuth convention; double application of ARP height |

### 7.8 ARP eccentricity and antenna height

| Aspect | Specification |
|---|---|
| Formulation | X_ARP = X_marker + R_enu→xyz(φ, λ) · [ΔE, ΔN, ΔH]ᵀ, from RINEX `ANTENNA: DELTA H/E/N` (H = ARP height above marker) |
| Rule | Applied **only** if the station coordinate refers to the marker (§ 8). If the coordinate refers to the ARP, ΔH/E/N SHALL NOT be applied again |
| Magnitude | cm to m |
| Mandatory | Yes |
| Validation | Station-log cross-check; consistency across all RINEX files of a station |
| Common errors | Double application; RINEX ΔH containing slant height; order H/E/N confusion |

### 7.9 Phase wind-up

| Aspect | Specification |
|---|---|
| Formulation | Wu et al. (1993): effective dipoles D' (satellite) and D (receiver) from body-frame and local-frame unit vectors; Δφ = sign(k̂·(D'×D)) · arccos(D'·D / (|D'||D|)) (cycles), unwrapped continuously along each arc; IF term = λ_NL · w |
| Inputs | Satellite attitude (§ 7.11), receiver local frame, LOS unit vector k̂ |
| Magnitude | Up to one cycle; large jumps during eclipse yaw manoeuvres [L] |
| Mandatory | Yes |
| Validation | Against RTKLIB windupcorr < 0.001 cycles; continuity (no ±1 cycle jumps except genuine) |
| Common errors | Not unwrapping (jumps of 1 cycle); nominal attitude during eclipse turns; forgetting the initial value is absorbed by the ambiguity (so only continuity matters in float, but wind-up errors matter in AR) |

### 7.10 Solid Earth tide

| Aspect | Specification |
|---|---|
| Formulation | IERS 2010 Ch. 7.1.1, degree-2/3 Love/Shida numbers with latitude dependence (Step 1) + frequency-dependent corrections (Step 2), as in official routine DEHANTTIDEINEL |
| Inputs | Station ECEF, Sun and Moon ECEF positions (JPL DE ephemeris or an analytic model whose error is shown < 0.1 mm in displacement — VERIFY), time (TT/UT1) |
| Tide system | Applied fully so that the resulting coordinate is **conventional tide-free**, matching ITRF/IGS **[F]**. Coordinates in any other tide system must be converted first (§ 9.6) |
| Magnitude | Radial up to ~30–40 cm; horizontal up to ~5 cm [L] |
| Mandatory | Yes (essential in fixed mode) |
| Validation | Against official IERS routine outputs for test stations/epochs < 0.1 mm |
| Common errors | Mixing mean-tide and tide-free coordinates; wrong time scale for Sun/Moon; applying only Step 1 |

### 7.11 Satellite attitude and eclipse handling

| Aspect | Specification |
|---|---|
| Formulation | Priority: (1) CODE ORBEX quaternions from the **same family** (FIN, COD0-RAP) → body-frame axes; (2) nominal yaw attitude (z to Earth centre, y ⟂ Sun–satellite–Earth plane) outside eclipse seasons; (3) during eclipse/noon turns without ORBEX: **exclude the satellite** (v1 default) |
| Exclusion window | Satellites whose |β| < β_crit (block-specific, from Kouba 2009 and later block IIF/III literature — VERIFY DURING IMPLEMENTATION values) during Earth-shadow passage plus post-shadow recovery, and during noon turns. Defaults held in configuration, not code |
| Mandatory | Yes (attitude); eclipse modelling optional (exclusion acceptable) |
| Validation | Wind-up and PCO continuity during eclipse seasons; residual inspection for eclipsing PRNs; comparison ORBEX vs nominal outside eclipse (< 0.1° yaw) |
| Common errors | Using nominal attitude through turns; quaternion convention (scalar-first vs scalar-last, rotation direction) — VERIFY with ORBEX spec |

### 7.12 Ocean tide loading

| Aspect | Specification |
|---|---|
| Formulation | IERS 2010 Ch. 7.1.2; 11 BLQ constituents expanded to 342 via admittance (official HARDISP algorithm) |
| Inputs | Per-station BLQ file from the Onsala Space Observatory loading service, model **FES2014b**, **with centre-of-mass correction (CMC)** consistent with CM-based orbits — VERIFY DURING IMPLEMENTATION the IGS/CODE convention |
| Order | BLQ rows: radial, tangential west, tangential south; amplitudes then phases (Greenwich phase lags) **[F: BLQ format]** |
| Magnitude | mm to several cm (larger at Indian coasts with high tidal range, e.g., Gulf of Khambhat) [L] |
| Mandatory | Yes |
| Validation | HARDISP official output for test cases < 0.1 mm |
| Common errors | Interpreting W/S as E/N (sign errors); phase sign; wrong model/CMC choice; station coordinate mismatch in BLQ request (> a few km) |

### 7.13 Pole tide (solid Earth)

| Aspect | Specification |
|---|---|
| Formulation | IERS Conventions Ch. 7.1.4 with the current conventional secular pole model — VERIFY which update applies |
| Inputs | x_p, y_p from CODE ERP |
| Magnitude | Few mm (radial up to ~2 cm) [L] |
| Mandatory | Yes (cheap) |
| Validation | Against IERS reference formulas < 0.05 mm |
| Common errors | Using x_p, y_p without subtracting the secular pole; arcsec vs mas units |
| Ocean pole tide loading | v2 (sub-mm to mm) |

### 7.14 Earth orientation

SP3 orbits are Earth-fixed in the product frame, so geometry needs no EOP **[F]**. EOP (CODE ERP: x_p, y_p, UT1−UTC) is needed for: Sun/Moon Earth-fixed positions (tides, attitude, wind-up) and the pole tide. Linear interpolation of daily ERP values is sufficient **[A]**; the ERP SHALL come from the same CODE family as the orbits.

### 7.15 Frame transformations (station coordinates)

Specified in § 9. All geometry is computed in the **product frame at the observation epoch**.

### 7.16 Atmospheric loading (tidal/non-tidal), hydrological loading

Not applied in v1 **[A]**. Their neglect is part of the fixed-coordinate error budget (§ 10.6) and must be assessed if coordinate pulls show seasonal signals.

### 7.17 Higher-order ionosphere

Not applied in v1 **[A]**. Expected sub-mm to ~mm effect on ZTD [L]; revisit for low-latitude high-TEC periods using CODE GIM (`COD0OPSFIN_…_GIM.INX`, same family), which does not violate the CODE-only policy.

### 7.18 Correction summary

| Correction | v1 | Typical size | Primary risk if wrong |
|---|---|---|---|
| SP3 orbit | M | — | edge errors, tier mixing |
| Satellite clock | M | — | interpolation, family mixing |
| Sagnac | M | ≤ 30 m | sign |
| Relativistic clock | M | ≤ 14 m | sign |
| Shapiro | M | ~19 mm | omission → ZTD bias |
| Satellite PCO/PCV | M | m / mm | SVN mapping, ANTEX mismatch |
| Receiver PCO/PCV + radome | M | cm | **ZTD bias** |
| ARP eccentricity | M | cm–m | double application |
| Phase wind-up | M | ≤ 1 cycle | eclipse errors |
| Solid tide | M | ≤ 40 cm | tide system mismatch |
| Ocean loading | M | mm–cm | W/S sign |
| Pole tide | M | mm | secular pole |
| Attitude/eclipse | M (exclusion allowed) | — | wind-up/PCO |
| Atmos./hydro loading | v2 | mm | seasonal coord. error |
| Higher-order iono | v2 | sub-mm–mm | — |


## 8. Station geometry definitions

| Term | Definition (this project) |
|---|---|
| **Marker** | The physical geodetic monument point (RINEX `MARKER NAME`/`MARKER NUMBER`) |
| **ARP (Antenna Reference Point)** | The antenna's mechanical reference point as defined in the IGS antenna definition (usually bottom of the antenna mount); ANTEX PCOs are relative to the ARP **[F]** |
| **RINEX antenna height/eccentricity** | `ANTENNA: DELTA H/E/N` — vertical height of ARP above marker (H) and horizontal eccentricities (E, N), metres **[F]** |
| **PCO** | Frequency-specific mean phase-centre offset from ARP (receiver: N, E, U; satellite: body x, y, z from CoM) |
| **PCV** | Elevation (and azimuth) dependent phase-centre variation relative to the mean phase centre |
| **APC** | Antenna phase centre = ARP + PCO (+ PCV per direction). The geometric range ρ is between satellite APC and receiver APC |
| **Radome** | 4-character code in ANTEX columns 17–20 of the antenna type field; calibrations are type+radome specific **[F]** |
| **Antenna type matching** | Exact 20-character match of RINEX `ANT # / TYPE` (type + radome) with ANTEX `TYPE / SERIAL NO`. Individual calibrations (by serial) MAY be used if present and flagged |
| **Radome fallback** | If type+radome is absent but type+`NONE` exists: allowed **only** with flag `RADOME_FALLBACK`, σ inflation (§ 10.6), and the station may not enter `production_fixed` unless the acceptance tests (§ 10.3) pass. If type is absent: `FAIL` (no guessing) |
| **Station coordinate reference point** | One of `MARKER`, `ARP`, `APC_L1`, `APC_IF`. Unknown → default `MARKER` with automatic check (§ 9.11–9.12); never blocking in fixed/constrained mode |
| **Reference frame / epoch** | Frame label (e.g., ITRF2008) + epoch (decimal year) of the coordinate |
| **Tide system** | `TIDE_FREE` (conventional tide-free, ITRF convention), `MEAN_TIDE`, or `ZERO_TIDE` |

---

## 9. SOI CORS ITRF2008@2005.0 → Current CODE/IGS Processing Frame

### 9.1 The problem

| Side | Frame | Epoch |
|---|---|---|
| SOI official coordinates | ITRF2008 **[F, project fact]** | 2005.0 **[F]** |
| CODE products | The IGS realisation used by CODE at the observation date, read from the product headers (§ 9.3) | Observation epoch (orbits are instantaneous) |

These are **not** interchangeable. Differences arise from (i) station motion between 2005.0 and the observation epoch (dominant: India-plate horizontal motion of the order of 5 cm/yr in ITRF [L — to be computed per station]), (ii) the frame transformation ITRF2008 → ITRF2014/ITRF2020 family (mm level), (iii) possible tide-system and reference-point differences, and (iv) the ANTEX calibration era used when SOI computed its coordinates (mm–cm vertical).

### 9.2 IGS realisations and their ITRF families **[F unless marked]**

| IGS realisation | ITRF family | Used in IGS products from | Transformation to ITRF of same family |
|---|---|---|---|
| IGS08 / IGb08 | ITRF2008 | 2011 / 2012 — VERIFY exact GPS weeks | Zero (by construction) — VERIFY |
| IGS14 / IGb14 | ITRF2014 | 2017 / 2020 — VERIFY exact GPS weeks | Zero — VERIFY |
| **IGS20** | ITRF2020 | GPS week 2238 (27 Nov 2022) | Zero |
| **IGb20** | ITRF2020-u2023 (aligned to ITRF2020) | GPS week 2352 (2 Feb 2025) | **All transformation parameters zero** (IGSMAIL-8543) |
| **IGc20** | ITRF2020-u2024 (aligned to ITRF2020) | GPS week 2401 (11 Jan 2026) | **All transformation parameters zero** (IGSMAIL-8634) |

Consequence: for a non-IGS station such as an SOI CORS, coordinates expressed in ITRF2020 at the observation epoch are directly consistent with IGS20, IGb20 and IGc20 products. The realisation label must still be read and recorded, because future realisations may not be zero-transformation.

### 9.3 Determining the product frame (no hard-coding)

The processor SHALL determine the frame from the products for each processing window:

1. SP3-d header "coordinate system" field (line 1) — e.g., `IGS20`, `IGb20`, `IGc20` — VERIFY DURING IMPLEMENTATION the exact labels CODE writes for each period.
2. Clock RINEX header (`# OF SOLN STA / TRF` or equivalent) — VERIFY field name in Clock RINEX 3.04.
3. Clock RINEX `SYS / PCVS APPLIED` → ANTEX name (e.g., igs20_WWWW) — must match the local ANTEX file.
4. A configuration table `frame_labels.yaml` maps label → ITRF family → transformation chain. **Unknown label → coordinate status `FRAME_UNKNOWN` → processing `FAIL` for fixed/constrained modes.** Orbit and clock labels must agree.

### 9.4 Velocities

| Case | Procedure | Status |
|---|---|---|
| **SOI official velocities available** (ITRF2008-consistent, per station) | Use them (V₀₈) | Preferred; request from SOI (**unresolved input**) |
| SOI velocities unavailable, **station has ≥ 2.5 years of RINEX** | Estimate velocity from the processor's own daily static PPP time series (linear trend + annual/semi-annual terms, with discontinuities at equipment changes/earthquakes). ≥ 2.5 yr is recommended to limit seasonal bias in velocity (Blewitt & Lavallée 2002) [L] | Scientifically preferred fallback |
| Unavailable and short history | **Horizontal:** ITRF2020 Plate Motion Model (Altamimi et al. 2023), v = Ω_India × X (+ ORB as specified by the PMM documentation), parameters read from the official ITRF2020-PMM file (itrf-ftp.ign.fr/pub/itrf/itrf2020 — VERIFY file and ORB usage). **Vertical:** no model → v_U = 0 with σ_vU = 3 mm/yr [A] | Degraded; usable only for **a priori** coordinates and for checks, not for fixed mode unless validated (§ 10.3) |
| Stations not on the rigid India plate (Himalayan front, NE India/Shillong region, Andaman–Nicobar with post-2004 post-seismic deformation, Kutch) | PMM **not valid** [L]; use own time series or nearby ITRF2020 stations with post-seismic models | Flag `NONLINEAR_MOTION_REGION` |

Plate-model accuracy limits **[L/E]**: rigid-plate residuals are typically ~1–2 mm/yr in stable interiors; over 20 years (2005 → 2026) a 1.5 mm/yr error gives ~3 cm horizontal; an unknown vertical rate of ±3 mm/yr gives ±6 cm vertical — **too large for fixed-mode ZTD** (≈ 12–24 mm ZTD). Hence the PMM fallback cannot by itself qualify a station for fixed mode.

### 9.5 Equations

**Step 1 — propagate in the source frame (ITRF2008) to the processing epoch t** (decimal year, midpoint of the processing window):

```
X₀₈(t) = X₀₈(t₀) + V₀₈ · (t − t₀),     t₀ = 2005.0                                  (9.1)
Σ_X(t) = Σ_X(t₀) + (t − t₀)² Σ_V   (plus cross terms if provided)                    (9.2)
```

**Step 2 — transform ITRF2008 → ITRF2020 at epoch t.** The official ITRF2020 → ITRF2008 parameters (IGN, `Transfo-ITRF2020_TRFs.txt`) **[F]**:

| Parameter | T1 (mm) | T2 (mm) | T3 (mm) | D (ppb) | R1, R2, R3 (mas) | Epoch |
|---|---|---|---|---|---|---|
| ITRF2020 → ITRF2008 | 0.2 | 1.0 | 3.3 | −0.29 | 0.00, 0.00, 0.00 | 2015.0 |
| rates (per yr) | 0.0 | −0.1 | 0.1 | 0.03 | 0.00, 0.00, 0.00 | |

Official model (from ITRF2020, X, to ITRF2008, X_S) **[F]**:

```
X_S = X + T + D·X + R·X,   R = [ 0  −R3  R2 ; R3  0  −R1 ; −R2  R1  0 ]
P(t) = P(2015.0) + Ṗ·(t − 2015.0)   for each parameter P
```

Inverse used here (ITRF2008 → ITRF2020), first order (second-order terms < 10⁻¹⁷ relative, negligible):

```
X₂₀(t) = X₀₈(t) − T(t) − D(t)·X₀₈(t) − R(t)·X₀₈(t)                                   (9.3)
V₂₀     = V₀₈ − Ṫ − Ḋ·X₀₈ − Ṙ·X₀₈                                                     (9.4)
```

Units: T in mm → m (×10⁻³); D in ppb → ×10⁻⁹; R in mas → rad (×4.8481368×10⁻⁹).

Magnitude check **[derived from the official parameters]**: at t = 2026.0, T ≈ (0.2, −0.1, 4.4) mm and D ≈ +0.04 ppb (~0.3 mm); at t = 2005.0, T ≈ (0.2, 2.0, 2.3) mm and D ≈ −0.59 ppb (~−3.8 mm radial). **The frame transformation is millimetre-level; the propagation (Step 1) is metre-level.** Both are mandatory.

**Step 3 — ITRF2020 → IGS20/IGb20/IGc20:** identity (§ 9.2).

**Historical data in the IGS14/IGb14 era:** transform ITRF2008 → ITRF2014 using the official ITRF2014 → ITRF2008 parameters (IGN `Transfo-ITRF2014_ITRFs.txt`) **[F]**: T = (1.6, 1.9, 2.4) mm, D = −0.02 ppb, R = 0 at epoch 2010.0; rates (0.0, 0.0, −0.1) mm/yr, Ḋ = 0.03 ppb/yr, with the same sign convention (from ITRF2014 to ITRF2008). **IGS08/IGb08 era:** no frame transformation (propagation only). Chains SHALL be table-driven, never composed ad hoc.

Equivalence test (Level 0): "propagate in ITRF2008 then transform at t" vs "transform position and velocity at t₀ then propagate in ITRF2020" must agree to < 0.1 mm.

**Step 4 — tide system and reference point** (§ 9.6), then ARP/APC geometry (§ 7.7–7.8).

### 9.6 Tide system and reference-point verification

- ITRF/IGS coordinates are **conventional tide-free** [F]. SOI's tide convention for the ITRF2008@2005.0 coordinates is **unknown**; default assumption: conventional tide-free (§ 9.11).
- If SOI coordinates are mean-tide, convert to tide-free using the IERS 2010 Ch. 7 permanent-tide expressions (VERIFY DURING IMPLEMENTATION the exact equations); the radial difference scales with P₂(sin φ) and is of order several cm at Indian latitudes [L] — far above ZTD tolerance.
- The coordinate reference point (marker/ARP/APC) must be declared (§ 8). If unknown, the default assumption is the marker, checked automatically (§ 9.11–9.12); processing never stops.
- **ANTEX era**: SOI's 2005.0 solution was necessarily computed with an older IGS antenna model (igs05/igs08 era — VERIFY with SOI). Changes in receiver/satellite antenna models between eras can shift estimated heights by mm to cm [L]. This cannot be corrected by transformation parameters; it is detected by the PPP comparison (§ 10.3, test C4).

### 9.7 What the transformation affects

| Component | Affected? | Reason |
|---|---|---|
| Receiver position vector | **Yes (metre-level)** | Propagation + frame |
| Satellite geometry (elevation/azimuth) | Negligible | 1 m / 20,000 km ≈ 10⁻⁷ rad |
| PCO/PCV application | Negligible | Local ENU orientation change ~10⁻⁷ rad |
| Solid tide, ocean loading, pole tide | Negligible | Smooth functions of position; m-level shifts change them by ≪ 0.01 mm |
| ZTD | **Yes** | Through range errors (§ 9.8) |

### 9.8 Consequences of using raw SOI ITRF2008@2005.0 coordinates for 2020–2026 data

| Error source | Size **[L/E]** | Effect |
|---|---|---|
| Unpropagated horizontal motion (15–21 yr × ~4–5.5 cm/yr) | **~0.6–1.2 m** | Azimuth- and elevation-dependent range errors up to ~1 m; far beyond what gradients (mm) absorb; ambiguities absorb only arc-constant parts; residuals, rejections, biased or failed ZTD |
| Unpropagated vertical motion | mm/yr × 20 yr (cm-level possible where subsidence occurs [L]) | ZTD bias ≈ β·δU |
| Frame ITRF2008 vs ITRF2020 | ~2–5 mm | ≈ 1 mm ZTD |
| Tide-system mismatch (if any) | several cm vertical | cm-level ZTD bias |
| ANTEX-era/reference-point mismatch | mm–cm (or dm–m for reference point) | mm–cm ZTD bias |

**Rule:** raw ITRF2008@2005.0 coordinates SHALL NOT be used as processing coordinates in any mode. They SHALL be retained as the `official` coordinate record.

### 9.9 Vertical coordinate error → ZTD (hypothesis to be measured)

Literature rule of thumb: ZTD bias ≈ β · δU with β ≈ 0.2–0.4 (depends on elevation cutoff and weighting; a fixed height that is too high produces a positive ZTD bias) [L; Santerre 1991 and subsequent work]:

| δU | ZTD bias (β = 0.2–0.4) | PWV bias (Π ≈ 0.16) |
|---|---|---|
| 1 cm | 2–4 mm | 0.3–0.6 mm |
| 3 cm | 6–12 mm | 1.0–1.9 mm |
| 5 cm | 10–20 mm | 1.6–3.2 mm |
| 10 cm | 20–40 mm | 3.2–6.4 mm |

The processor-specific β SHALL be measured (§ 10.5) and stored per station; it is used to convert coordinate uncertainty into a ZTD uncertainty-budget term.

### 9.10 Coordinate module pipeline and coordinate object

```
SOI official XYZ [ITRF2008 @ 2005.0]
  → metadata validation (frame, epoch, velocity, tide system, reference point, antenna)
  → velocity/trajectory model (SOI velocity | own time series | PMM fallback)
  → propagate to processing epoch t (9.1–9.2)
  → frame transformation to product frame family (9.3–9.4, table-driven)
  → tide-system conversion (if needed) → reference-point conversion to ARP
  → processing coordinate → estimator
```

The **coordinate object** (immutable, one per station per processing epoch) contains:

| Field | Content |
|---|---|
| `X, Y, Z` | Processing coordinate (m), reference point = ARP |
| `lat, lon, h_ell` | Geodetic latitude/longitude (deg) and ellipsoidal height (m), GRS80 |
| `h_orth`, `geoid_model` | Orthometric height (m) and geoid used (e.g., EGM2008) — needed for pressure reduction |
| `frame` | e.g., ITRF2020 family; `realisation_label` from product (e.g., IGc20) |
| `epoch` | Decimal year of the processing coordinate |
| `velocity`, `velocity_frame`, `velocity_source` | (m/yr) and its origin (SOI / own time series / PMM) |
| `reference_epoch_of_source` | 2005.0 for SOI |
| `tide_system` | `TIDE_FREE` after conversion |
| `reference_point` | `ARP` after conversion; source reference point recorded |
| `covariance` | 3 × 3 (m²), propagated |
| `source` | `soi_transformed`, `ppp_validated`, `official_raw` (never used for processing) |
| `transformation_history` | Ordered list of steps with parameters, parameter-file names and hashes |
| `validation_status` | From § 10.3 |
| `original_official_record` | Unmodified SOI values and metadata |


### 9.11 Default assumptions for unknown SOI metadata (never-stop policy, Revision 1.1)

| Unknown item | Default assumption | Basis | Printed / recorded as | Automatic check (§ 9.12) |
|---|---|---|---|---|
| Frame / epoch | ITRF2008 @ 2005.0 | **Known project fact** | — | — |
| Station coordinates themselves | Taken from `--xyz` or `stations.csv`; if absent: no SOI coordinate → `static` mode with a priori from RINEX `APPROX POSITION XYZ` (or single-point solution from the `.n` file) | — | WARNING `NO_SOI_COORD` | — |
| Velocity | Horizontal: ITRF2020 Plate Motion Model, India plate rotation pole (+ ORB as specified by the PMM documentation), parameters stored in `settings.py` with citation — values to be copied from the official PMM table (VERIFY DURING IMPLEMENTATION; never invented). Vertical: 0 mm/yr with σ = 3 mm/yr | Rigid-plate model is the best available generic model for stations without velocities [L] | WARNING `VELOCITY_ASSUMED_PMM` | Same-day static comparison |
| Stations in deforming zones (Himalayan front, north-east, Andaman–Nicobar, Kutch) | Same PMM default, plus warning, based on simple lat/lon boxes in `settings.py` | PMM not valid there [L] | WARNING `NONLINEAR_MOTION_REGION` | Same-day static comparison |
| Tide system | **Conventional tide-free** | ITRF and IGS coordinates are conventional tide-free [F]; SOI coordinates were published in an ITRF frame, so the same convention is the most probable [A] | WARNING `TIDE_SYSTEM_ASSUMED` | Discrepancy compared with the predicted mean-tide/tide-free offset at the station latitude (logged as a hypothesis only) |
| Coordinate reference point | **Marker**; RINEX `ANTENNA: DELTA H/E/N` applied to obtain the ARP | Geodetic station coordinates conventionally refer to the marker, and RINEX eccentricities are defined from marker to ARP [F: RINEX 3.05]. When ΔH/E/N = 0 the question is moot [F] | WARNING `REFPOINT_ASSUMED_MARKER` (only if ΔH ≠ 0) | Auto-switch to ARP interpretation if the vertical discrepancy matches ΔH (§ 9.12) |
| Antenna type + radome | From RINEX header | — | — | ANTEX lookup |
| Radome not in ANTEX | Use type + `NONE` calibration | Common practice; residual error mm–cm [L] | WARNING `RADOME_FALLBACK` | — |
| Antenna type not in ANTEX | PCO/PCV = 0; mode forced to `static` (height absorbs part of the mean offset) | Only option without calibration | WARNING `NO_ANTENNA_CALIBRATION`; ZTD flagged | — |
| Antenna model era of SOI coordinates | Assumed consistent with igs20.atx | Cannot be recovered | Recorded | Same-day static comparison |
| Ocean-loading BLQ | If `blq/<STATION>.blq` exists: use it. Otherwise **continue without ocean loading** | The Onsala service delivers BLQ by e-mail; no documented scripted interface was found [F: provider FAQ] | WARNING `NO_OTL` with one-line instructions for obtaining the file once per station | — |
| Surface pressure (for PWV) | VMF3 grid ZHD (ECMWF-based) reduced to antenna height | Automatic, anonymous download | INFO `ZHD_SOURCE=VMF3` | — |
| T_m (for PWV) | GPT3 T_m | Automatic, no account needed | INFO `TM_SOURCE=GPT3` | — |
| ERA5 (optional upgrade) | Used only if the user has a Copernicus CDS key configured | ERA5 needs a CDS account | INFO | — |

Magnitude of what the assumptions can cost (for the record) [L/E]:
- Tide-system error (if SOI were mean-tide): several cm vertical → several mm to ~1 cm ZTD in fixed mode — caught by § 9.12.
- Reference-point error: equal to ΔH (cm–m) → caught by § 9.12.
- Missing ocean loading: mm inland (e.g., Deccan), up to cm at macro-tidal coasts; mostly semi-diurnal/diurnal and partly absorbed, but it can leave mm-level sub-daily ZTD artefacts.
- PMM velocity error: ~1–2 mm/yr horizontally in the stable interior over 21 years ≈ 2–4 cm horizontal (minor for ZTD); unknown vertical rate is the main risk → caught by § 9.12.

### 9.12 Automatic coordinate consistency check and `auto` mode (default)

Run inside every station-day run, without any user action:

```
X_soi  = SOI coordinate → propagate (9.1, PMM default) → transform (9.3) → tide-free → ARP (§ 9.11)
Pass A = float PPP in static_estimated mode on the same data → X_A (σ_A)
Δ      = X_soi − X_A  in local N, E, U
```

Decision rules (thresholds in `settings.py`) [A]:

| Condition | Action | `coordinate_source` |
|---|---|---|
| No SOI coordinate given | Use Pass A as the final solution | `ppp_same_day` |
| ΔH ≠ 0 in RINEX and |Δ_U − ΔH| ≤ 2 cm | Reinterpret SOI coordinate as ARP, recompute X_soi and Δ, log `REFPOINT_SWITCHED_TO_ARP` | — |
| |Δ_N|, |Δ_E| ≤ 2 cm **and** |Δ_U| ≤ 2 cm | Pass B: **fixed** mode at X_soi (project preference) | `soi_transformed` |
| Otherwise | Keep Pass A (static) as the final solution; print the discrepancy table and likely causes (velocity, tide system, reference point, antenna era) | `ppp_same_day` |

User override: `--mode fixed` forces fixed at X_soi (warning printed if the check failed); `--mode static` skips Pass B; `--mode constrained` uses X_soi with σ = (3, 3, 5) mm.

Why this is scientifically safe: for a 24-hour static solution, estimating the coordinate changes ZTD by only ~1–3 mm STD [E], whereas fixing a coordinate that is wrong by a few cm biases ZTD by several mm (§ 9.9). So the program fixes the SOI coordinate only when the data confirm it, and otherwise uses the safer static solution. Every run prints and records Δ, so a multi-day history of Δ per station builds up automatically (useful later to estimate station velocities and to settle the unknown metadata).


---

## 10. Coordinate modes and coordinate validation

### 10.1 Modes

| Mode | Coordinates | Use |
|---|---|---|
| `static_estimated` | One constant XYZ per processing window, a priori σ = 1 m | Validation, coordinate derivation, diagnostics |
| `constrained` | Constant XYZ, prior σ_N, σ_E, σ_U (default 3, 3, 5 mm) | Recommended production mode once validated |
| `fixed` | Not estimated; processing coordinate from § 9 | Production mode (project preference) once validated |

### 10.2 Production coordinate source hierarchy (configurable)

1. `soi_transformed` — SOI official coordinate propagated and transformed (§ 9), **only if** acceptance tests C1–C8 pass. This honours the project preference for official SOI coordinates.
2. `ppp_validated` — the processor's own multi-week static PPP solution (CODE Final), propagated with the best available velocity, used when (1) fails or velocities are unavailable. Refreshed monthly and after any discontinuity.
3. `ppp_same_day` — the same-day static solution chosen by `auto` mode (§ 9.12) when the SOI coordinate is not confirmed by the data.

**Revision 1.1:** in the program this hierarchy is applied automatically by `auto` mode (§ 9.12); nothing blocks a run.

### 10.3 Multi-day coordinate diagnostics (optional; never blocking in Revision 1.1)

These tests are run as an offline validation study (Stage 4) on accumulated runs; the per-run automatic check is § 9.12.

| Test | Procedure | Pass | Warn | Fail |
|---|---|---|---|---|
| **C1 Metadata completeness** | Frame, epoch, velocity (or declared absence), tide system, reference point, antenna+radome, ΔH/E/N, receiver | All present | Velocity absent (PMM fallback) | Reference point or antenna unknown |
| **C2 Transformation** | Apply § 9; propagate covariance | Completed | — | Frame label unknown |
| **C3 Static PPP** | ≥ 14 (preferably 28) consecutive days, CODE Final, `static_estimated`, frozen config; robust daily mean X_PPP | ≥ 14 valid days | 7–13 days | < 7 days |
| **C4 SOI vs PPP** | ΔN, ΔE, ΔU = X_SOI,proc − X_PPP [A thresholds] | |ΔH| ≤ 10 mm and |ΔU| ≤ 15 mm | |ΔU| ≤ 30 mm | larger |
| **C5 Repeatability** | Daily scatter of X_PPP | σ_H ≤ 4 mm, σ_U ≤ 8 mm [A] | ≤ 2× | larger (station unstable) |
| **C6 ZTD consistency** | ZTD_fixed(X_selected) − ZTD_static over the C3 period | |mean| ≤ 2 mm, STD ≤ 3 mm [A] | |mean| ≤ 4 mm | larger |
| **C7 Antenna/ARP** | Exact ANTEX match (type+radome), ΔH/E/N consistent over all files and with station log | Exact | `RADOME_FALLBACK` | Type missing or inconsistent heights |
| **C8 Frame/epoch/tide** | Product frame label known; tide system tide-free after conversion; epoch within 1 day of window midpoint | All true | — | Any false |

Outcome:
- `PRODUCTION_FIXED_OK` — all Pass (C1 Warn allowed if C4 passes with `ppp_validated`).
- `CONSTRAINED_ONLY` — any Warn; constrained mode with σ inflated by the observed |ΔU| (σ_U = max(5 mm, |ΔU|)).
- `REJECTED` — any Fail; the station is flagged in reports and its runs use `auto`/`static` behaviour.

The processor SHALL NOT silently accept an inconsistent coordinate: inconsistencies are printed, flagged and handled by § 9.12 (fixed only when the data confirm the coordinate).

### 10.4 Diagnostic comparison workflow

For each station, report and plot:
- SOI official (raw) vs SOI transformed vs X_PPP (N, E, U differences, with uncertainties).
- ZTD from `fixed(SOI transformed)`, `fixed(ppp_validated)`, `constrained`, `static_estimated`: differences time series, mean/STD per day.
- Phase residuals by elevation and azimuth (skyplot) for each mode.
- Coordinate pull `z_c` time series in constrained mode.

### 10.5 Coordinate sensitivity experiment (mandatory in Stage 4)

| Element | Specification |
|---|---|
| Stations | All pilot stations, then ≥ 5 representative SOI stations (flat inland, coastal, mountain, NE, IGP) |
| Period | ≥ 14 days per station, including one monsoon and one dry period |
| Variants | **A** raw SOI ITRF2008@2005.0; **B** SOI propagated to epoch (ITRF2008); **C** propagated + transformed to product frame; **D** processor-estimated (`static_estimated`, then `ppp_validated` fixed) |
| Vertical perturbations | Applied to C (or D): δU = ±1, ±3, ±5, ±10 cm |
| Horizontal perturbations | δN, δE = ±1, ±3, ±10 cm (separately) |
| Outputs | ΔZTD, ΔZWD, ΔPWV (mean, STD, diurnal composite), convergence time, phase residual RMS by elevation, rejection counts, fixing rate (when AR exists) |
| Derived | Height-to-ZTD coupling coefficient β = regression slope of mean ΔZTD vs δU (per station, cutoff and weighting); horizontal coupling to gradients |
| Acceptance | Experiment completed, β estimated with uncertainty; β outside [0.1, 0.6] triggers investigation [A] |

### 10.6 Coordinate error budget for ZTD

σ_ZTD,coord = β · σ_U,coord (σ_U,coord from the coordinate object), plus radome-fallback inflation (default 5 mm ZTD [A], to be refined). Reported in metadata, and included in PWV uncertainty as a systematic term.


## 11. CODE-only product management

### 11.1 Policy **[A — project decision]**

- Analysis-Centre-estimated products (orbit, clock, ERP, satellite code/phase biases, attitude) SHALL come from **one CODE solution family** for the whole processing window.
- Families (in fallback order): **`COD0OPSFIN`** (Final) → **`CODMOPSRAP`** (final rapid) → **`COD0OPSRAP`** (early rapid). Families are never mixed, including COD0 vs CODM rapid and OPS vs MGX.
- Non-AC inputs (ANTEX, BLQ, VMF3/GPT3, leap seconds, meteorological data, DEM, geoid) are not AC-mixing; each is recorded in provenance.
- **CODE Ultra-Rapid (`COD0OPSULT`) SHALL NOT be used for PPP**: CODE's product table lists ultra-rapid orbits and ERP but no clock RINEX and no OSB **[F: CODE 2024 IGS Technical Report]**. If a future CODE release adds them, adoption requires a design change and Level-1/2 re-validation, not a configuration edit.
- No silent substitution: a missing product produces an explicit status (§ 11.4).

### 11.2 Required products per family **[F: CODE 2024 IGS TR table; VERIFY DURING IMPLEMENTATION filenames/compression against live archive]**

| Product | `COD0OPSFIN` | `CODMOPSRAP` | `COD0OPSRAP` | Content / sampling | Float | AR | Days needed |
|---|---|---|---|---|---|---|---|
| Orbit SP3 | `COD0OPSFIN_YYYYDDD0000_01D_05M_ORB.SP3.gz` | `CODMOPSRAP_YYYYDDD0000_01D_05M_ORB.SP3` | `COD0OPSRAP_YYYYDDD0000_01D_05M_ORB.SP3` | CoM positions (+clocks), 5 min, GPS+GLO+GAL | **Req.** | **Req.** | D−1, D, D+1 (window) |
| Clock | `COD0OPSFIN_…_01D_30S_CLK.CLK.gz` (also `_05S_`) | `CODMOPSRAP_…_01D_30S_CLK.CLK` | `COD0OPSRAP_…_01D_30S_CLK.CLK` | Satellite clocks 30 s (5 s FIN option); Clock RINEX 3.04; integer-cycle (phase) clocks since 3 Jun 2018 | **Req.** | **Req.** | All window days |
| ERP | `COD0OPSFIN_…_01D_01D_ERP.ERP.gz` | `CODMOPSRAP_…_01D_01D_ERP.ERP` | `COD0OPSRAP_…_01D_01D_ERP.ERP` | x_p, y_p, UT1−UTC, daily | **Req.** | **Req.** | All window days |
| OSB (Bias-SINEX 1.00) | `COD0OPSFIN_…_01D_01D_OSB.BIA.gz` | `CODMOPSRAP_…_01D_01D_OSB.BIA` | `COD0OPSRAP_…_01D_01D_OSB.BIA` | Satellite code + phase OSBs consistent with the family's clocks | Cond.¹ | **Req.** | All window days |
| Attitude (ORBEX) | `COD0OPSFIN_…_01D_30S_ATT.OBX.gz` | **Not provided** | `COD0OPSRAP_…_01D_30S_ATT.OBX` | Quaternions, 30 s | Opt.² | Rec.² | All window days |
| GIM (IONEX) | `COD0OPSFIN_…_01D_01H_GIM.INX.gz` | — | `COD0OPSRAP_…_01D_01H_GIM.INX.gz` | 1 h | Not v1 | Not v1 | — |
| Troposphere SINEX | `COD0OPSFIN_…_01D_01H_TRO.TRO.gz` | — | `COD0OPSRAP_…_01D_01H_TRO.TRO` | CODE ZTD at CODE network stations | **Validation only** | — | — |

¹ Float mode requires OSBs unless the code observables used are exactly the clock-reference signals (C1W, C2W), whose OSBs cancel in the IF combination by CODE's datum convention — VERIFY DURING IMPLEMENTATION in the Bias-SINEX header/CODE README. Otherwise missing OSB → `OK_FLOAT_ONLY` is not allowed with C1C etc.; the processor must switch to C1W or flag `CODE_BIAS_MISSING`.
² Missing ORBEX → nominal attitude + eclipse exclusion (§ 7.11); borrowing COD0OPSRAP attitude for a CODMOPSRAP solution is disallowed by default (explicit flag `ATT_FROM_OTHER_SOLUTION` if enabled).

Other facts **[F]**: final products available after ~2 weeks; rapid series after ~18 h (AIUB); ultra-rapid 4×/day, ~3 h. Integer-cycle clocks and phase biases from GPS week 2004 (3 Jun 2018) (Schaer et al. 2021). Legacy short names before 27 Nov 2022 (`CODwwwwd.EPH.Z`, `CODwwwwd_v3.CLK.Z`, `CODwwwwd.CLK_05S.Z`, `CODwwwwd.BIA.Z`, `CODwwwwd.OBX.Z`, `CODwwwwd.ERP.Z`; rapid `.EPH_R`, `.CLK_R`, `_M` for final rapid) — held in configuration, VERIFY. CODE's daily clock/phase-bias continuity across midnight for GPS in the Final series is documented (Schaer et al. 2021) — VERIFY per tier before bridging ambiguities across day boundaries; default: ambiguities are **reset at product-day boundaries** unless continuity is verified by a Level-0 test.

### 11.3 Float vs AR product requirements

| Requirement | FLOAT | AR |
|---|---|---|
| SP3 (D−1..D+1), CLK 30 s, ERP — single family | Required | Required |
| OSB code (non-clock-reference signals) | Required if used | Required |
| OSB phase for every phase observable used | Not used | **Required ("no phase bias, no AR" — CODE IAR README)** |
| Integer-cycle clocks (date ≥ 2018-06-03) | Not required | Required |
| ORBEX | Optional | Recommended (wind-up during eclipse); without it, eclipsing satellites excluded from AR |
| ANTEX matching clock header | Required | Required |

### 11.4 Processing status codes (exact conditions)

Evaluated per station per processing window, after product search and verification:

| Status | Condition (all must hold) |
|---|---|
| **`OK_AR`** | (a) A single family F ∈ {COD0OPSFIN, CODMOPSRAP, COD0OPSRAP} provides SP3 + CLK + ERP + OSB for every window day, all verified (§ 12.4); (b) window date ≥ 2018-06-03; (c) phase OSBs exist for all selected phase observables of ≥ 4 satellites at ≥ 90 % of epochs [A]; (d) ANTEX and frame label verified; (e) AR requested in configuration |
| **`OK_FLOAT_ONLY`** | (a) As above for SP3 + CLK + ERP (+ code OSB if required by observables), all verified; and (b) at least one AR condition fails. A **reason code** is mandatory: `AR_NOT_REQUESTED`, `NO_PHASE_OSB`, `PRE_2018`, `OSB_SIGNAL_MISMATCH`, `ATT_MISSING_ECLIPSE_EXCLUDED` (informational) |
| **`WAIT_FOR_TIER`** | No PPP-capable family is complete **and** the observation window end is younger than `max_expected_latency` for the next tier (default: rapid 48 h, final 21 days [A]); e.g., only ultra-rapid exists. No processing is run |
| **`FAIL`** | Any of: no PPP-capable family complete after `max_expected_latency`; verification failure that cannot be recovered by re-download (corrupt, wrong AC/family in header, epoch coverage incomplete); ANTEX named in clock header not obtainable; frame label unknown (fixed/constrained modes); RINEX unreadable or no usable GPS dual-frequency data. A **failure code** is mandatory |

Tier upgrades: a window processed with `CODMOPSRAP` or `COD0OPSRAP` receives status `OK_*` with `product_tier = RAPID` and is **queued for reprocessing** when `COD0OPSFIN` becomes complete (§ 19.3).

Adjacent-day rule: the whole window uses the **lowest common tier** that is complete for all window days. If D+1 is not yet available in any PPP-capable tier, the window is shortened to end at D 24:00 (no mixing), outputs near the end get `EDGE`, and the SP3 interpolation uses asymmetric nodes flagged `ORBIT_EDGE` **[A]**.

### 11.5 OSB application rule **[F: CODE IAR README]**

Satellite OSBs are applied to **all** observations used (code and phase), per exact RINEX observation code: O_corrected = O_measured − OSB^s(code) with the sign convention defined in Bias-SINEX 1.00 and CODE's README (`O_measured = O_computed + B_satellite`) — VERIFY DURING IMPLEMENTATION the sign by reproducing CODE's IF clock consistency. Satellite PCO is **not** applied to the Melbourne–Wübbena combination (CODE IAR README).

---

## 12. Product manager, downloader and cache

### 12.1 Date logic

1. Read first and last observation epochs from **data records** (not filename, not header only).
2. Processing days: all GPS days intersecting [first, last]; window days = those ± 1 (orbit/clock/ERP/OSB/ATT) for overlap windows (§ 5.10).
3. GPS week/day, year/DOY computed from GPST; product days are GPST days (as in CODE naming) — VERIFY that CODE daily files span 00:00–23:59:30/24:00 GPST.
4. VMF3 grids: all 6-hourly epochs covering [window start − 6 h, window end + 6 h] — VERIFY file naming/epochs on the VMF data server.
5. ERP coverage: daily files for window days (FIN weekly collection file optional).
6. Output: a **ProductRequest** listing each product type, family candidates, days, and required time span.

### 12.2 Search and fallback algorithm

```
for family in configured_order [COD0OPSFIN, CODMOPSRAP, COD0OPSRAP]:
    for each required product and window day:
        look in cache → else try each configured source (URL template list) → verify
    if all required (for requested mode) present and verified: select family; break
evaluate status (§ 11.4); write manifest; return
```

Sources are **configuration**: an ordered list of URL templates per family and product with tokens `{yyyy} {ddd} {wwww} {d} {family} {filename}`. Candidate sources (VERIFY DURING IMPLEMENTATION): AIUB HTTP (`http://ftp.aiub.unibe.ch/CODE/{yyyy}/`), AIUB S3 mirror (`https://zhw-b.s3.cloud.switch.ch/aiub/CODE/{yyyy}/`), CDDIS (`https://cddis.nasa.gov/archive/gnss/products/{wwww}/`, Earthdata login). **No directory tree or filename is hard-coded.**

### 12.3 Download mechanics

| Requirement | Specification |
|---|---|
| Proxy | HTTP(S)_PROXY environment + explicit config; custom CA bundle path; no TLS verification bypass |
| Auth | `.netrc` for Earthdata (CDDIS), with redirect/cookie handling |
| Timeouts / retry | connect 20 s, read 120 s; 5 retries, exponential backoff 2ⁿ s with jitter, max 300 s [A] |
| Atomicity | download to `*.part` in cache staging, verify, then atomic rename |
| Rate limiting | per-host concurrency limit (default 2) |
| Offline mode | cache only; missing → status per § 11.4 without network |
| Logging | every attempt (URL, HTTP status, bytes, duration) in the product log |

### 12.4 Verification (every file, before use)

1. Decompression test (gzip/Unix-compress/Hatanaka where applicable).
2. Checksum where published (e.g., CDDIS SHA512/MD5 manifest); always compute SHA-256 for provenance.
3. Header parse: AC = COD; family/campaign/tier consistent with filename; format version.
4. Epoch coverage: first/last epoch cover the product day; SP3 epoch interval 5 min; CLK satellite epoch interval 30 s (or 5 s).
5. GPS satellites present (≥ 24 [A]); satellite list consistent across SP3/CLK/OSB.
6. Frame label (SP3) and ANTEX name (CLK `SYS / PCVS APPLIED`) extracted and compared across files.
7. OSB file: signals present for the station's selected observables; time validity covers the day.
8. Result stored as a sidecar JSON next to the cached file.

### 12.5 Cache layout (configurable root)

```
products/
  cache/
    COD/{family}/{yyyy}/{ddd}/{original_filename}
    COD/{family}/{yyyy}/{ddd}/{original_filename}.meta.json      # sha256, source URL, retrieval time, header summary
    ANTEX/{antex_name}.atx (+ .meta.json)
    VMF3/{grid_type}/{yyyy}/{filename}
    GPT3/{filename}
    BLQ/{station}.blq
    LEAPSEC/{filename}
  staging/                                                        # *.part files
  (no database: each file has a .meta.json sidecar)
```

Files are stored under their original names; the SHA-256 in each sidecar prevents silent replacement. A file whose remote version changes (same name, different hash) is stored as a new version, and the change is logged.

### 12.6 Provenance manifest per station-day (JSON)

`manifest_id` (UUID), `station`, `day`, `rinex_files` (+ SHA-256, header summary), `analysis_centre = "COD"`, `family`, `tier`, `version_char` (`0`/`M`), `products[]` (type, filename, SHA-256, sampling, time span, header creation date, source URL), `frame_label`, `antex` (name + SHA-256), `blq` (model, CMC flag, SHA-256), `vmf3_files[]`, `processing_mode` (FLOAT/AR), `status` + reason/failure code, `fallbacks[]` / `missing[]`, `coordinate_object_id`, `coordinate_mode`, `config_hash` (SHA-256 of canonicalised configuration), `software_version` (semantic version + git commit), `run_timestamp`, `supersedes` / `superseded_by`.

---

## 13. Preprocessing and quality control

### 13.1 Pipeline order

```
RINEX read → header validation → observable selection → epoch alignment/decimation
→ satellite exclusion (products, health, ANTEX) → receiver clock jump detection
→ elevation computation & cutoff → cycle-slip detection (LLI, MW, GF, gaps)
→ arc construction & short-arc removal → [OSB application] → IF formation
→ estimation with innovation screening (IGG-III) → smoothing
→ post-fit residual editing (iterate) → final smoothing → outputs
```

### 13.2 RINEX reading

- RINEX 3.02–3.05 (4.0x SHOULD be supported for observation files; VERIFY differences).
- Accept Hatanaka-compressed files (decompress via a verified CRX2RNX implementation or library).
- Header fields required: `MARKER NAME`, `ANT # / TYPE`, `ANTENNA: DELTA H/E/N`, `SYS / # / OBS TYPES`, `TIME OF FIRST OBS` (time system), `INTERVAL` (derived from data if absent), `REC # / TYPE / VERS`.
- Header values (antenna, ΔH/E/N) are compared with the previous run of the same station (from its last manifest); a change → WARNING and flag `HEADER_REGISTRY_MISMATCH` (non-blocking).

### 13.3 Observable selection (configurable priority lists, no vendor hard-coding)

| Band | Code priority (default) | Phase priority (default) |
|---|---|---|
| L1 | C1W, C1C, C1X | L1C, L1W, L1X |
| L2 | C2W, C2L, C2S, C2X | L2W, L2L, L2S, L2X |

Rules:
- Selection is made **per satellite per arc**; an arc SHALL use one code and one phase observable per band throughout. A change starts a new arc.
- In AR mode, only phase observables with a phase OSB in the family's OSB file are eligible.
- Selected observables are recorded per arc in the diagnostic output.
- Phase alignment (quarter-cycle between L2 tracking modes) is assumed per RINEX 3 conventions and OSBs; VERIFY with CODE OSB signal list.

### 13.4 Epoch alignment and decimation

Processing epochs SHALL coincide with CODE 30-s clock epochs by default. Higher-rate RINEX is decimated to 30 s. RINEX epochs not aligned to clock epochs → clock interpolation (linear) with flag `CLOCK_INTERP` (5-s FIN clocks preferred when available).

### 13.5 Satellite exclusion

Excluded (per epoch or arc), with a logged reason: missing/flagged SP3 position or clock; missing CLK; SP3/CLK event flags (manoeuvre); no valid ANTEX entry (SVN/PRN validity); eclipse window without attitude (§ 7.11); OSB missing for selected signal (AR mode: excluded from AR only); configured blacklist; broadcast health flag unhealthy (optional, if navigation file available).

### 13.6 Melbourne–Wübbena slip detection

```
MW = (f₁L₁ − f₂L₂)/(f₁ − f₂) − (f₁P₁ + f₂P₂)/(f₁ + f₂)            [m];   in cycles: MW / λ_WL
```

Running arc mean μ and standard deviation s (recursive, Blewitt 1990). Slip if |MW_k/λ_WL − μ| > max(k_MW · s, n_MW,min), then confirmation by the next epoch (to separate outliers from slips).

| Parameter | Default [A] | Validation |
|---|---|---|
| k_MW | 4 | Injected-slip test (§ 30.4) |
| n_MW,min | 1.0 cycle | idem |
| s floor | 0.25 cycle | idem |
| confirmation epochs | 1 | idem |

### 13.7 Geometry-free slip detection

```
GF = L₁ − L₂      [m]  (= ionospheric term + λ₁N₁ − λ₂N₂ + const)
```

Predict GF_k by a low-order polynomial (default degree 2) fitted to the previous n (default 10) epochs of the arc; slip if |GF_k − GF_pred| > max(k_GF · σ_fit, GF_min · g(Δt, e)) with GF_min = 0.05 m for 30-s data [A] scaled for sampling interval and elevation.

Low-latitude ionosphere: India lies near the equatorial ionisation anomaly crest; post-sunset TEC gradients and scintillation (especially equinox seasons and high solar activity) raise false-alarm rates [L]. A rate-of-TEC index (ROTI) is computed per station-epoch; when ROTI exceeds `roti_high` (default 0.5 TECU/min [A]) GF thresholds are scaled up and MW becomes the primary detector; periods are flagged `IONO_ACTIVE`. (1 TECU corresponds to ≈ 0.105 m in L₁ − L₂ [derived].)

### 13.8 Receiver clock jumps

Detect a common jump in (P_IF − L_IF) of ≈ k·(c × 1 ms) across ≥ 80 % of satellites [A]. If integer-ms: repair by adding k·c·10⁻³ m to the phase (or code) observations after the jump to restore consistency; log `CLOCK_JUMP_REPAIRED`. Otherwise reset all arcs (§ 5.7).

### 13.9 Elevation screening and SNR

Cutoff 7° (configurable). Optional SNR mask (default off): reject if S1 < 30 dB-Hz or S2 < 25 dB-Hz [A]; optional SNR-based weighting (§ 5.4).

### 13.10 Outlier detection and robust weighting

**Forward pass (innovations):** standardised innovation v̄ᵢ = vᵢ/√S_ii. IGG-III weight factor:
```
w = 1                                   if |v̄| ≤ k₀
w = (k₀/|v̄|) · ((k₁ − |v̄|)/(k₁ − k₀))² if k₀ < |v̄| ≤ k₁
w = 0 (reject)                          if |v̄| > k₁
```
Defaults k₀ = 2.0, k₁ = 5.0 [A]; code and phase screened separately; if > 30 % of phase observations at an epoch exceed k₁, suspect undetected clock jump/common problem → epoch rejected and logged.

**Post-fit residual editing (after smoothing):** smoothed residuals rᵢ; flag phase |rᵢ| > max(4 σ̂_phase(e), 0.04 m) and code |rᵢ| > max(4 σ̂_code(e), 4 m) [A]; phase outliers inside an arc → split arc there (treat as slip); isolated code outliers → remove; iterate estimation up to 3 passes or until no new flags.

### 13.11 Post-processing of a cycle slip

Slip repair is **not** required in v1: every detected slip creates a new ambiguity (scientifically sufficient for static 24-h ZTD with smoothing) **[A]**. Slip repair MAY be added later only if validation shows ZTD degradation from frequent resets.

### 13.12 Configurable QC parameters summary

All thresholds in §§ 5, 6, 13 live in `config/qc.yaml` with the defaults above, units stated, and a reference to the validation experiment that tunes them.


## 14. Processing flow per station-day (float)

1. Build the coordinate object (§ 9.10) from `--xyz`/`stations.csv` with the § 9.11 defaults; the mode is chosen by § 9.12 (never blocking).
2. Product request → product status (§ 11.4). Stop on `WAIT_FOR_TIER` / `FAIL`.
3. Build the processing window (30 h default; shorter at data/product edges).
4. Preprocess (§ 13).
5. Forward EKF (§ 5), innovation screening.
6. RTS smoother.
7. Residual editing loop (≤ 3 passes).
8. Extract 5-min outputs (§ 6.5) for day D only; compute flags.
9. Write station products (§ 16–17), manifest (§ 12.6), log and quick-look plot (§ 33.5).

---

## 15. PPP-AR subsystem (separate, later)

### 15.1 Preconditions (hard)

- Float PPP has passed its Stage-3 and Stage-5 gates (§ 35).
- Product status `OK_AR` (single CODE family; integer-cycle clocks; phase OSBs) — **absolute family consistency is a hard requirement**; any inconsistency → AR disabled for the window.
- AR output is written as a **separate solution layer** (`solution_type = FIXED`); float products are never overwritten.

### 15.2 Observation preparation

OSB-corrected code and phase (§ 11.5); satellite PCO not applied to MW (CODE IAR README).

### 15.3 Wide-lane (WL)

Per arc: MW mean in cycles, μ_MW, with standard error σ_μ = s/√n_eff (n_eff accounting for time correlation; default decorrelation time 300 s [A]).

Receiver WL bias is common to all satellites; it is eliminated by **between-satellite single differences (SD)** relative to a reference arc, choosing an **independent set** of SD ambiguities over the window (reference arc = longest overlapping, high-elevation arc; e.g., Blewitt 1989 / Ge et al. 2008 approach).

WL fixing criterion per SD: |frac(μ_SD)| ≤ 0.25 cycle and σ_μ,SD ≤ 0.1 cycle and arc overlap ≥ `min_arc_ar_s` (default 1200 s) [A]; optionally the Dong & Bock (1989) probabilistic test. Record fixed/unfixed and fractional parts.

### 15.4 Narrow-lane (NL)

From the **AR forward pass** (all ambiguities retained to the window end, § 5.8), obtain float B_IF and full covariance. For WL-fixed SDs:

```
N₁,SD(float) = ( B_IF,SD − (λ_WL f₂ /(f₁ + f₂)) · N_WL,SD ) / λ_NL        (15.1)
```

with covariance propagated from B_IF.

### 15.5 Integer estimation and validation

- **LAMBDA** (Teunissen 1995) decorrelation + integer least-squares search on the NL SD vector.
- Validation: (i) bootstrapped success-rate P_s ≥ 0.999 (Teunissen 1998); (ii) ratio test with fixed failure-rate threshold (Verhagen & Teunissen 2013) or, as v1 fallback, ratio ≥ 3.0 [A]; both required.
- **Partial AR:** if full-set validation fails, iteratively remove the ambiguity with the largest conditional variance (or lowest elevation / shortest arc) and retry until validation passes or fewer than 4 remain → then `AR_FAILED` (float retained).

### 15.6 Fixed solution

Accepted SD integers are converted back to IF constraints (4.3) and applied as pseudo-observations with σ_fix = 1 mm [A] at the corresponding ambiguity birth epochs; forward + RTS re-run → fixed ZTD series.

**Consistency gate per window:** ZTD_fixed − ZTD_float: |mean| ≤ 2 mm and max |Δ| ≤ 10 mm [A]; else `AR_REJECTED_CONSISTENCY` and the fixed layer is not published for that window.

### 15.7 Reinitialisation

A slip ends the arc; its ambiguity is independent. Integer fixes never propagate across arcs; ambiguities are reset at product-day boundaries unless continuity is verified (§ 11.2).

### 15.8 Why a wrong fix is worse than float for ZTD

A wrong NL integer introduces a range bias of k·λ_NL ≈ k × 10.7 cm on that satellite for the whole arc, **imposed as a hard constraint with a very small variance**. The estimator must then explain it through the remaining parameters; because the ZTD partial is m(e), the bias leaks into ZTD (and gradients) with an elevation-dependent weighting — several mm to cm over the arc — while the formal σ becomes *smaller*. The float solution simply absorbs such offsets into a free ambiguity without bias. Therefore AR is only allowed with strict validation, partial fixing and the float-consistency gate.

### 15.9 AR metrics (per window and per 5-min row)

`n_amb_total`, `n_WL_candidates`, `n_WL_fixed`, `n_NL_candidates`, `n_NL_fixed`, `fixing_ratio` (NL fixed / NL candidates), `ratio_test_value`, `success_rate`, `n_rejected_fixes`, `ambiguity_age_s` (per arc), `reset_causes` (counts: slip MW, slip GF, LLI, gap, observable change, clock jump, outlier, day boundary), `AR_status` ∈ {`FLOAT`, `FIXED`, `PARTIAL`, `AR_FAILED`, `AR_REJECTED_CONSISTENCY`, `AR_NOT_AVAILABLE`}.

---

## 16. Station-level output schema (System 1)

One row per station per 5-min UTC epoch. Constant per-window fields are stored once in NetCDF global/variable attributes and as repeated columns in CSV.

| Column | Type | Unit | Description |
|---|---|---|---|
| `timestamp_utc` | string ISO 8601 / int64 | UTC | Output epoch (exact 5-min UTC mark) |
| `timestamp_gpst` | string ISO 8601 | GPST | Same instant in GPS time (UTC + leap seconds) |
| `station` | string | — | 4- or 9-character station ID |
| `latitude` | float64 | degree_north | Geodetic latitude of ARP (processing coordinate) |
| `longitude` | float64 | degree_east | Geodetic longitude of ARP |
| `height` | float64 | m | Ellipsoidal height of ARP (GRS80) |
| `ZTD` | float64 | **mm** | Smoothed zenith total delay at ARP height |
| `sigma_ZTD` | float64 | mm | Formal 1σ (smoothed), unscaled |
| `ZHD` | float64 | mm | A priori ZHD₀ used in PPP (source in `zhd_source`) |
| `ZWD` | float64 | mm | ZTD − ZHD₀ (a priori split) |
| `sigma_ZWD` | float64 | mm | = sigma_ZTD (ZHD₀ held fixed) |
| `G_N`, `G_E` | float64 | mm | Horizontal gradient parameters (Chen & Herring) |
| `receiver_clock` | float64 | ns | Receiver clock estimate (nearest epoch) |
| `n_sat` | int16 | — | Satellites used at the nearest epoch |
| `n_rejected` | int16 | — | Observations rejected in the ±2.5-min window |
| `n_amb_fixed` | int16 | — | Fixed ambiguities active (0 in float layer) |
| `convergence_flag` | enum string | — | `CONVERGED`, `NOT_CONVERGED`, `EDGE`, `REINIT` |
| `AR_status` | enum string | — | § 15.9 (`AR_NOT_AVAILABLE` in float-only runs) |
| `QC_flag` | uint32 bitmask | — | § 16.1 |
| `product_AC` | string | — | `COD` |
| `product_tier` | enum | — | `FINAL`, `RAPID_M` (CODMOPSRAP), `RAPID_0` (COD0OPSRAP) |
| `product_files` | string (JSON list) | — | Exact filenames used (also in manifest) |
| `coordinate_mode` | enum | — | `fixed`, `constrained`, `static_estimated` |
| `coordinate_source` | enum | — | `soi_transformed`, `ppp_validated` |
| `mapping_function` | string | — | e.g., `VMF3` (`GPT3` if fallback) |
| `ZWD_process_noise` | float64 | mm/√h | σ_ZWD used |
| `software_version` | string | — | Semantic version + git commit |
| `solution_type` | enum | — | `FLOAT` or `FIXED` |
| `zhd_source` | enum | — | `BAROMETER`, `VMF3_GRID`, `GPT3` |
| `manifest_id` | string (UUID) | — | Link to provenance manifest |
| `config_hash` | string | — | SHA-256 of canonical configuration |
| `elevation_cutoff` | float32 | degree | Cutoff used |

### 16.1 QC_flag bits

| Bit | Name | Meaning |
|---|---|---|
| 0 | `NO_DATA` | No bracketing observation epochs |
| 1 | `FEW_SATS` | n_sat < 5 |
| 2 | `EDGE` | Segment edge |
| 3 | `REINIT` | Majority of ambiguities recently reset |
| 4 | `IONO_ACTIVE` | High ROTI |
| 5 | `HIGH_REJECT` | Rejection fraction > 20 % in window |
| 6 | `COORD_PULL` | Constrained-mode coordinate pull |z_c,U| > 3 |
| 7 | `RADOME_FALLBACK` | Antenna calibration fallback |
| 8 | `ORBIT_EDGE` | Asymmetric SP3 interpolation |
| 9 | `CLOCK_INTERP` | Clock interpolation beyond 30-s nodes |
| 10 | `ZHD_CLIMATOLOGY` | ZHD₀ from GPT3 |
| 11 | `ECLIPSE_EXCLUSION` | Satellites excluded for attitude |
| 12 | `RAPID_TIER` | Product tier not Final |
| 13 | `HEADER_REGISTRY_MISMATCH` | Station metadata inconsistency |
| 14 | `SIGMA_HIGH` | sigma_ZTD > 15 mm |
| 15–31 | reserved | |

---

## 17. Station product formats and units

- **Units:** all delays and gradients in **mm**; PWV in **mm** (numerically equal to kg m⁻² for ρ_w = 1000 kg m⁻³); clock in ns; heights in m; angles in degrees. No mixed units within a file type.
- **CSV** (primary tabular output, ISO timestamps, header comment line with schema version). Parquet MAY be added later for large archives.
- **NetCDF-4** (CF-1.10 + ACDD-1.3, `featureType = "timeSeries"`), one file per station-day: dimension `time`; scalar station coordinates (`lat`, `lon`, `alt`), variables as in § 16 with `units`, `long_name`, `_FillValue`; `time` encoded as `seconds since 1970-01-01T00:00:00Z` (int64), calendar `standard` (UTC; leap seconds not represented — documented in `comment`); global attributes: provenance fields (§ 12.6), `schema_version`, `product_tier`, `solution_type`.
- **SINEX_TRO 2.00** per station-day (5-min records; TROTOT, STDDEV, TGNTOT, TGETOT in mm) — VERIFY DURING IMPLEMENTATION field names and header blocks against the current SINEX_TRO 2.00 specification.
- File naming: `{STATION}_{YYYY}{DDD}_{SOLUTION}_{TIER}_ZTD_v{schema}.{nc|parquet|csv|tro}`.

---

## 18. PWV conversion module (separate)

### 18.1 Equations

```
ZHD_met = 0.0022768 · P_ant / (1 − 0.00266 cos 2φ − 0.28×10⁻⁶ H_ant)        [m; P hPa; H m]   (18.1)
ZWD_met = ZTD − ZHD_met                                                                        (18.2)
Π       = 10⁶ / [ ρ_w · R_v · (k₃/T_m + k₂′) ]                                                  (18.3)
PWV     = Π · ZWD_met                                                                           (18.4)
```

Constants (frozen default, Bevis et al. 1994) **[F]**: ρ_w = 1000 kg m⁻³; R_v = 461.5 J kg⁻¹ K⁻¹; k₂′ = 22.1 K hPa⁻¹ = 0.221 K Pa⁻¹; k₃ = 3.739 × 10⁵ K² hPa⁻¹ = 3739 K² Pa⁻¹ (SI units in (18.3)). Π ≈ 0.155–0.167 for T_m 270–295 K. Alternative refractivity constant sets (e.g., Rüeger 2002) MAY be configured for sensitivity studies; the set used is recorded.

### 18.2 Required inputs and sources (priority)

| Input | Priority 1 | Priority 2 | Priority 3 | Not acceptable for research PWV |
|---|---|---|---|---|
| P_ant | Station barometer (calibrated), reduced to ARP height | ERA5 (hourly) interpolated in the vertical profile to H_ant | Operational NWP (same method) | GPT3 (climatology) — flagged only |
| T_m | ERA5 profile integration: T_m = ∫(e/T)dz / ∫(e/T²)dz | GPT3 T_m | Regional India T_m–T_s relation (to be fitted with radiosondes) / Bevis T_m = 70.2 + 0.72 T_s | — |
| T_s | Station sensor | ERA5 2-m temperature (height-adjusted) | — | — |
| H_ant | h_ell − N_geoid (EGM2008) | — | — | Using h_ell as orthometric (India geoid undulation ~ −100 to −30 m → up to ~12 hPa error) |

### 18.3 Pressure height reduction

Hypsometric equation with virtual temperature: P₂ = P₁ · exp(−g·(H₂ − H₁) / (R_d · T̄_v)), R_d = 287.05 J kg⁻¹ K⁻¹, g = 9.80665 m s⁻² (or latitude-dependent) [F]. For ERA5, the preferred method is interpolation in log-pressure along the model/pressure-level geopotential profile to H_ant; geopotential height to orthometric height conversion per ECMWF definitions — VERIFY DURING IMPLEMENTATION.

### 18.4 Uncertainty propagation

```
σ²_PWV,ppp  = Π² · σ²_ZTD                                   (from System 1, scaled by validated factor s_ZTD)
σ²_ZHD      = (0.0022768 / f)² · σ²_P
σ²_PWV,conv = Π² · σ²_ZHD + ZWD² · σ²_Π,   σ_Π/Π ≈ (k₃/T_m²)/(k₃/T_m + k₂′) · σ_Tm
σ²_PWV,tot  = σ²_PWV,ppp + σ²_PWV,conv (+ Π² σ²_ZTD,coord, § 10.6)
```

PPP uncertainty and conversion uncertainty are reported **separately** (`sigma_PWV_ppp`, `sigma_PWV_conv`, `sigma_PWV_total`). Default σ_P: barometer 0.3 hPa, ERA5 1.0 hPa, NWP 1.5 hPa; σ_Tm: ERA5 2 K, GPT3 4 K, Bevis 5 K [A — to be validated with radiosondes].

### 18.5 Station PWV product schema (additional columns)

`PWV`, `sigma_PWV_ppp`, `sigma_PWV_conv`, `sigma_PWV_total` (mm); `P_ant` (hPa), `P_source`; `T_m` (K), `Tm_source`; `Pi` (dimensionless); `ZHD_met`, `ZWD_met` (mm); `constants_set`; plus all System 1 keys (`timestamp_utc`, `station`, `manifest_id`).

---

## 19. Provenance and reprocessing

### 19.1 Provenance objects

| Object | Keyed by | Contains |
|---|---|---|
| Input file record | SHA-256 | name, source, retrieval time, verification result |
| Station-day manifest | UUID | § 12.6 |
| Coordinate object | UUID | § 9.10 |
| PWV run record | UUID | met sources, constants, files, System 1 manifest ID |
| Grid time-step record | UUID | station product IDs used (with tiers), background files, DEM/geoid/mask versions, covariance parameters, config hash, software version |

### 19.2 Configuration hashing

Configuration is canonicalised (sorted keys, normalised numbers, explicit units) and hashed (SHA-256). Every output carries the hash; the full configuration is archived under `provenance/configs/{hash}.yaml`.

### 19.3 Reprocessing (Rapid → Final)

- Outputs are **immutable**. A reprocessing run writes new files with a new `manifest_id`; the new manifest records `supersedes` and the old manifest is updated with `superseded_by` (plain JSON files, no database).
- Selection rule (applied by scanning manifests): `best_available` (Final > Rapid-M > Rapid-0, highest software version within a declared processing campaign) and `as_of(date)` for reproducibility.
- Triggers: (i) a higher tier becomes complete; (ii) a declared reprocessing campaign (new software/config version).
- System 2 grids built from station products are rebuilt when the share of Final station inputs changes from < 100 % to 100 % for that day, or on campaign reprocessing; old grids retained per retention policy (default: keep all Final; keep Rapid for 90 days after Final exists [A]).

---

## 20. Interface contract System 1 → System 2

| Item | Contract |
|---|---|
| Transport | Station ZTD/PWV CSV + NetCDF files and manifests in the results folder |
| Required fields | `timestamp_utc`, `station`, `latitude`, `longitude`, `height`, `H_ant` (orthometric), `PWV`, `sigma_PWV_total`, `sigma_PWV_ppp`, `ZTD`, `sigma_ZTD`, `ZHD_met`, `ZWD_met`, `T_m`, `P_ant`, `QC_flag`, `convergence_flag`, `product_tier`, `solution_type`, `coordinate_mode`, `manifest_id` |
| Schema versioning | Semantic `schema_version`; System 2 rejects unknown major versions |
| Independence | System 2 never calls System 1 code; it can be tested with synthetic station products |
| Station information | Read from the manifests (coordinates, heights, coordinate source, antenna/radome, flags) |


# PART B — SYSTEM 2: NATIONAL PWV/ZWD/ZTD MAPPING

## 21. Inputs and station ingestion

| Input | Specification |
|---|---|
| Station products | System 1/PWV products via § 20; `solution_type = FLOAT` by default (FIXED allowed only after AR validation) |
| Background (research/reprocessing) | **ERA5** hourly, 0.25°, pressure-level or model-level profiles (T, q, geopotential, surface pressure, orography) — VERIFY exact variables/levels at implementation |
| Background (operational/rapid) | An operational NWP analysis/forecast for India (e.g., NCMRWF NCUM or ECMWF IFS open data) — VERIFY availability, resolution, licensing and latency; source recorded per grid |
| Double counting check | Record whether the background assimilated ground-based GNSS ZTD from these stations; if yes, residuals are not independent and CV results must be reported separately |
| DEM | Copernicus GLO-30 or GLO-90 (orthometric, EGM2008) aggregated to each grid: cell mean height and cell height standard deviation — VERIFY vertical datum |
| Geoid | EGM2008 |
| Masks | Land/sea mask (consistent with DEM), India analysis domain 6–38°N, 68–98°E |

Station ingestion QC:
1. Drop rows with `QC_flag` bits {NO_DATA, FEW_SATS, SIGMA_HIGH} or `convergence_flag ∈ {NOT_CONVERGED}`; `EDGE`/`REINIT` rows down-weighted (σ × 2) [A].
2. Station blacklist/greylist: long-term (30-day) median of (PWV_GNSS − PWV_bg) for the station deviating from the median of its neighbours (within 150 km) by > 2 mm → greylist (review); > 4 mm → blacklist [A].
3. Per-slot spatial consistency: leave-one-out standardised residual |r_i − r̂_i(−i)| / √(σ²_pred + σ²_i) > 5 → slot outlier (excluded, logged); iterate once [A].

---

## 22. Mapping formulation (baseline)

### 22.1 Why not interpolate raw ZTD

ZHD decreases by ≈ 0.27 mm per metre of height (≈ −0.12 hPa m⁻¹) [F-derived] and PWV decreases roughly exponentially with a scale height of ~2 km [L]. Interpolating raw ZTD or raw PWV between stations at different heights produces errors of tens to hundreds of mm in ZTD over the Himalaya, Western Ghats and north-east hills. **Raw ZTD SHALL NOT be interpolated.**

### 22.2 Background at arbitrary height

For horizontal position x, height H (orthometric), time t (background time-interpolated linearly between analysis times):

```
p_s(x, H, t)   = pressure at height H from the background profile (log-p interpolation in geopotential height; below model orography: hypsometric extrapolation)
PWV_bg(x,H,t)  = (1 / (ρ_w · g)) ∫₀^{p_s(x,H,t)} q dp        [kg m⁻² = mm]                  (22.1)
ZHD_bg(x,H,t)  = Saastamoinen(p_s(x,H,t), φ, H)                                                 (22.2)
T_m,bg(x,H,t)  = ∫(e/T)dz / ∫(e/T²)dz  from H upward                                             (22.3)
```

Horizontal interpolation of the background to station positions: bilinear on the native grid **[A]**.

### 22.3 Residuals at stations (15-min slot τ)

```
r_i(τ) = PWV_GNSS,i(τ) − PWV_bg(x_i, H_i, τ)                                                    (22.4)
σ²_i(τ) = σ²_PWV,total,i(τ) + σ²_bg,repr,i                                                      (22.5)
```

σ_bg,repr,i covers background interpolation to the station (default 0.5 mm [A], estimated in CV).

### 22.4 Regression kriging / kriging with external drift (baseline)

Model for the residual field:
```
r(x) = β₀ + Σ_k β_k z_k(x) + ε(x),    ε ~ GP(0, C)                                              (22.6)
```
Default drift covariates z_k: none beyond β₀; optional: H/1000 (height-dependent background bias), PWV_bg (proportional bias). Covariates are admitted only if they reduce CV error (§ 31).

Covariance (default exponential, Matérn ν = 1/2; ν = 3/2 as alternative):
```
C(x_i, x_j) = σ_s² · exp(−d_ij / L) · exp(−|ΔH_ij| / L_H)  +  δ_ij (τ² + σ²_i)                 (22.7)
```
- d_ij: great-circle distance (km).
- L: horizontal e-folding range; L_H: vertical range (optional; default off → L_H = ∞).
- τ²: nugget (representativeness/micro-scale); σ²_i: heteroscedastic station noise (22.5).
- **Anisotropy** (option): geometric anisotropy in local LAEA coordinates (angle θ, ratio a) — research option for Western Ghats/Himalayan fronts.
- **Regional covariance:** parameters estimated per covariance region (initially climate/terrain zones, to be refined) and season.

Hyperparameter estimation: restricted maximum likelihood (REML) over a rolling 30-day window per region and season (pooled over slots), with per-slot re-scaling of σ_s² by the robust variance of that slot's residuals; bounds L ∈ [L_min, L_max] = [20, 500] km [A]. Per-slot fits are diagnostics only.

Prediction at grid cell centre x₀ (local neighbourhood: nearest `k_max` = 40 stations within `R_max` = max(3L, 300 km) [A]):
```
[ C  F ] [ w ]   [ c₀ ]
[ Fᵀ 0 ] [ μ ] = [ f₀ ]          (universal kriging system; F = [1, z_k(x_i)], f₀ = [1, z_k(x₀)])
r̂(x₀) = wᵀ r                                                                                     (22.8)
σ²_UK(x₀) = C(0) − wᵀ c₀ − μᵀ f₀   (sign convention per standard UK derivation; VERIFY in Level-0 test vs a reference kriging library)
```

### 22.5 Reconstruction on the grid

```
H_cell            = mean DEM height of the cell (orthometric)
PWV_grid(x₀,τ)    = PWV_bg(x₀, H_cell, τ) + r̂(x₀, τ)                                           (22.9)
σ²_PWV,grid       = σ²_UK + σ²_bg,repr,cell + σ²_height,cell                                  (22.10)
Π_grid            = Π(T_m,bg(x₀, H_cell, τ))
ZWD_grid          = PWV_grid / Π_grid
ZHD_grid          = ZHD_bg(x₀, H_cell, τ)
ZTD_grid          = ZHD_grid + ZWD_grid
σ²_ZTD,grid       ≈ (σ_PWV,grid / Π_grid)² + σ²_ZHD,bg
```

σ_height,cell accounts for applying the residual at H_cell rather than at station height (default: proportional to |H_cell − H_nearest| and to cell height std; coefficients from CV [E]).

Assumption (to be tested in CV stratified by ΔH): the residual r is approximately height-independent after the background has supplied the vertical structure.

### 22.6 GNSS information fraction

```
I(x₀) = 1 − σ²_UK(x₀) / σ²_prior(x₀),   σ²_prior = σ_s² + τ² (+ trend uncertainty)            (22.11)
```

I → 0 where no station constrains the cell (product ≈ background); I → 1 near dense stations.

### 22.7 Point service

PWV at arbitrary (x, H): PWV_bg(x, H, τ) + r̂(x, τ), with σ per (22.10); used for validation at withheld stations and for users needing site-specific values in complex terrain.

---

## 23. Method catalogue

| Method | Role | Inputs | Assumptions | Covariance | Complexity (n stations, m cells) | Uncertainty | Sparse areas | Validation |
|---|---|---|---|---|---|---|---|---|
| **IDW** | Quick-look only | Residuals | None explicit | None | O(n·m) | None | Bull's-eyes, no flag basis | Not validated as product |
| Ordinary kriging | Component/benchmark | Height-normalised residuals | Stationary mean | Fitted (22.7) | O(k³·m) local | Kriging variance | Reverts to mean; flagged by σ | LOSO, block CV |
| **Regression kriging / UK with NWP external drift** | **Baseline production** | Residuals vs background, covariates | Residual stationarity within region/season | (22.7), REML | O(k³·m) local | σ_UK + terms (22.10) | Reverts to background; I → 0 | Full § 31 |
| Gaussian process regression | Research (equivalent maths; learned anisotropy/non-stationarity) | As baseline | Kernel form | Learned kernels | O(n³) exact; sparse/local approx. | Posterior σ | As baseline | § 31, compare to baseline |
| Spatio-temporal OI / Kalman filter | **Advanced (first upgrade)** | Residual time series + background | Separable or advective space-time covariance | C(d, Δt) | O(k³·m·T) local | Posterior σ | Temporal persistence helps short gaps | § 31 + temporal CV |
| NWP data assimilation of ZTD | Advanced, external (NWP centres) | Station ZTD + obs operator | NWP system | Background error covariance | NWP-scale | Analysis error | Model dynamics | Forecast skill |
| ML residual correction (RF/GBM/NN) | Research: learn background bias vs terrain/season | Residuals + covariates | Stationarity of learned relations | Implicit | Training-heavy | Only via ensembles/quantiles | Poor extrapolation | Strict spatial CV; never default |
| RBF / thin-plate splines | Not recommended | — | Smoothness | None | O(n³) | None | Overshoot | — |
| Natural neighbour | Not recommended | — | — | None | O(n log n) | None | No extrapolation | — |
| Tomography | Regional research in dense clusters only | Slant delays | Voxel model, constraints | Regularisation | High | Posterior (model-dependent) | Not applicable | Radiosonde profiles |

---

## 24. Supportable resolution and effective resolution

### 24.1 Quantities (per region R, weather regime s)

| Symbol | Definition |
|---|---|
| d_nn | Median nearest-neighbour distance among valid stations in R (km) |
| L_s | Fitted covariance range (e-folding, km) of residual field in R, regime s |
| E(d) | CV RMSE at withheld stations as a function of distance d to the nearest remaining station, from LOSO and thinning; fitted as E(d) = √(E₀² + (a·d^α)²) |
| S(d) | Skill over background: S = 1 − MSE_method / MSE_background |
| RMSE_Δ | CV RMSE (at withheld stations, via point reconstruction using grid spacing Δ for the height/representativeness term) |

### 24.2 Supportable grid spacing (candidates Δ ∈ {0.5°, 0.25°, 0.125°, 0.1°})

A region R supports spacing Δ in regime s if **all** hold [A — thresholds configurable, fixed after Stage 8 experiment]:
1. **Density:** d_nn ≤ k_d · Δ_km, k_d = 2 (at least ~one station per 2 × 2 cells along a line).
2. **Correlation:** d_nn ≤ L_s (stations lie within the correlation range).
3. **Skill:** S(d_nn) ≥ S_min = 0.2.
4. **Accuracy:** E(Δ_km / 2) ≤ E_target (initial hypothesis 3 mm PWV [E]).
5. **Incremental benefit:** RMSE_Δ < RMSE_{next coarser} − δ_min, significant at 95 % by spatial block bootstrap (δ_min = 0.1 mm [A]).

Δ_sup(R, s) = the finest Δ satisfying all; if none, R is supported only at the background level (`Class 4` for GNSS-corrected product).

### 24.3 Effective atmospheric resolution (reported layer)

```
Res_eff(x) = 2 · d_nn,local(x)           if I(x) ≥ 0.5
Res_eff(x) = Res_bg                       if I(x) < 0.5
```

Res_bg = effective resolution of the background (for ERA5 several grid lengths, order 100 km [L] — to be quantified via spectral analysis in Stage 8 [E]). Grid spacing never enters Res_eff.

### 24.4 Supportable resolution map

Produced per season from the § 31 experiment on 1° tiles (and aggregated to covariance regions): layer `supportable_grid_spacing` (categorical: 0.1°, 0.125°, 0.25°, 0.5°, background-only) and `effective_resolution_km`. National products are always generated on the fixed national grid; the map tells users where the grid is supported.

---

## 25. Temporal design

| Item | Specification |
|---|---|
| Station sampling | 5 min (System 1) |
| Map slot times | :00, :15, :30, :45 UTC |
| Station window | Slot τ uses station 5-min values at τ − 5, τ, τ + 5 min (window [τ − 7.5, τ + 7.5) min) |
| Averaging | Mean of available values; requires ≥ 2 of 3; else station absent for that slot |
| Window σ | σ_window = mean of member σ (conservative: assumes full correlation, since adjacent values are strongly correlated) [A] |
| Hourly archive | Separate kriging run using hourly station means (window [t − 30, t + 30) min, ≥ 8 of 12 values); consistency check against the mean of the four 15-min maps |
| Temporal covariance (advanced method) | Estimated from station residual series (e.g., exponential in Δt with time scale T_c, advective coupling option) — Stage 8 research |
| 5-min national maps | **Experimental only**, in regions where Δ_sup ≤ 0.125° and temporal CV (§ 31) shows 5-min maps have skill over 15-min maps; labelled `EXPERIMENTAL` |
| Missing stations | Excluded for that slot; `station_availability` recorded; coverage classes recomputed every slot |
| Latency | RAPID maps: after CODE rapid products (~18 h) and the PPP run for the day; background = operational NWP. FINAL maps: after CODE Final (~2 weeks) and ERA5 availability. True near-real-time (< 1–2 h) is not achievable under the CODE-only policy |

---

## 26. Coverage classes (per cell, per slot)

| Class | Provisional criteria (all configurable; final values from § 31) | Main PWV/ZWD/ZTD layers |
|---|---|---|
| **0** | Outside domain or sea beyond coastal buffer (default 10 km [A]) | Fill value |
| **1 — Well covered** | nearest valid station d₁ ≤ 25 km; n₅₀ ≥ 3; |H_cell − H_nearest| ≤ 300 m; I ≥ 0.7 | Filled |
| **2 — Moderate** | d₁ ≤ 50 km; n₁₀₀ ≥ 2; I ≥ 0.4 | Filled |
| **3 — Sparse** | d₁ ≤ d* (empirical support limit; initial 150 km [E]); I ≥ 0.1 | Filled, flagged "background-dominated" |
| **4 — Unsupported** | Otherwise | **Fill value** in GNSS-corrected layers; `background_PWV` layer still available |

Cells in class 1 with |ΔH| > 300 m are demoted to class 2; class 2 with |ΔH| > 800 m are demoted to class 3 [A].

---

## 27. Grid products

### 27.1 Grids

| Grid | Spacing | Extent | Alignment |
|---|---|---|---|
| National | 0.25° | 6–38°N, 68–98°E (128 × 120 cells) | Cell centres at k·0.25° + 0.125° — VERIFY alignment choice against ERA5 grid (ERA5 points lie at multiples of 0.25°); choose ERA5-coincident centres if background sampling is to be exact |
| Regional nests | 0.125° | Tiles where Δ_sup ≤ 0.125° | Nest exactly inside 0.25° cells |
| Coarse | 0.5° | Full domain | Aggregation of 0.25° (variance-aware) |

CRS: geographic WGS84 (EPSG:4326), CF `grid_mapping` variable `crs` with `grid_mapping_name = "latitude_longitude"`. Distances computed on the sphere/ellipsoid; LAEA projection only internally for anisotropy.

### 27.2 Layers

| Variable | Units | Type | Class |
|---|---|---|---|
| `PWV` | mm | float32 | **Primary** (CF standard name `lwe_thickness_of_atmosphere_mass_content_of_water_vapor` — VERIFY exact CF name) |
| `ZWD` | mm | float32 | Primary |
| `ZTD` | mm | float32 | Primary |
| `ZHD` | mm | float32 | Primary (background-derived) |
| `sigma_PWV` | mm | float32 | **Primary** |
| `sigma_ZTD` | mm | float32 | Primary |
| `confidence_class` | 0–4 | int8 | **Primary** (flag_values/flag_meanings) |
| `background_PWV` | mm | float32 | Diagnostic |
| `gnss_residual` | mm | float32 | Diagnostic |
| `distance_to_nearest_station` | km | float32 | QA |
| `n_stations_50km`, `n_stations_100km` | count | int16 | QA |
| `gnss_information_fraction` | 1 | float32 | QA |
| `dem_height`, `dem_height_std` | m | float32 | Ancillary (static) |
| `surface_pressure` (at H_cell) | hPa | float32 | Ancillary |
| `Tm` | K | float32 | Ancillary |
| `effective_resolution_km` | km | float32 | QA |
| `station_availability` | per-station 0/1 on `station` dimension | int8 | QA |
| `land_sea_mask` | 0/1 | int8 | Ancillary (static) |
| `supportable_grid_spacing` | categorical | int8 | QA (static per season/version) |

### 27.3 File formats

| Role | Format | Structure |
|---|---|---|
| **Primary scientific archive** | **NetCDF-4/HDF5, CF-1.10 + ACDD-1.3** | One file per grid per day per tier: dims `time` (96 for 15-min; 24 for hourly), `lat`, `lon`, `station` (for availability and contributing-station list), `bnds`; `time` = `minutes since 2000-01-01 00:00:00` UTC, int32, calendar `standard`; `lat_bnds`, `lon_bnds`; zlib level 4 + shuffle; float32; chunks (time = 4, lat = full, lon = full) for 0.25°, (4, 128, 120) tiles for 0.125°; global attributes: ACDD + provenance (tier, CODE families of contributing stations, background source, config hash, software version, supportable-resolution version) |
| Analysis store | Zarr (v2 or v3 — VERIFY library support) | Same variables and CF attributes; rechunked for time-series access (time = 2880 (30 days of 15-min), lat = 32, lon = 32) |
| GIS distribution | Cloud-Optimised GeoTIFF | Selected derived layers (hourly PWV, sigma_PWV, confidence_class), one band per file, tagged with tier/time; **not** an archive |
| File naming | `INPWV_{GRID}_{STEP}_{YYYYDDD}_{TIER}_v{X.Y}.nc`, e.g., `INPWV_G025_15M_2026154_FINAL_v1.0.nc` | |

---

## 28. Computational requirements (1,000 stations) [E]

| Item | Estimate |
|---|---|
| Station PPP solutions | 1,000 station-days/day; 288 outputs each → 288,000 ZTD values/day |
| PPP compute | Python/NumPy ~0.5–3 min per station-day [E] → ~10–50 CPU-h/day → ~0.5–2 h wall-clock on a 32-core server with parallel processes |
| Grid cells | 0.5°: 3,840 (~1,150 land); 0.25°: 15,360 (~4,600 land); 0.125°: 61,440 (~18,300 land); 0.1°: 96,000 (~28,700 land) |
| Kriging | local k = 40: ~0.1–0.2 s per 0.25° map; 96 maps/day → < 1 min/day; REML fitting: minutes/day |
| Background preparation | ERA5 profile integration for all cells and stations: dominated by I/O; minutes per day |
| Storage (12 float32 layers, uncompressed) | 0.25°/15-min ≈ 70 MB/day (~26 GB/yr); 0.125° nationwide equivalent ≈ 280 MB/day; compression 2–4× |
| RINEX | ~1–2 TB/yr (30-s, Hatanaka) |
| Hardware | One 32–64-core server, 128 GB RAM, ≥ 10 TB storage |
| Bottleneck | Data access, product latency and QC — not computation |


# PART C — DATA, VALIDATION, ARCHITECTURE, ROADMAP

## 29. Station metadata and development datasets

### 29.1 Metadata required before processing a station

| Item | Source | If unknown (Revision 1.1: never blocking — § 9.11 default applies) |
|---|---|---|
| Station ID, name, DOMES (if any) | SOI | Yes |
| Receiver type, firmware, change history | Station log / RINEX headers | Yes (history) |
| Antenna type + radome (ANTEX code), serial, change history | Station log / RINEX | **Yes** |
| ARP height and eccentricities (ΔH/E/N), history | Station log / RINEX | **Yes** |
| Official coordinates: X, Y, Z | SOI | Yes |
| Frame and epoch | SOI: **ITRF2008 @ 2005.0 [F]** | Yes |
| Velocities (and their frame) | SOI — **unresolved** | No (fallback § 9.4) but limits fixed mode |
| Tide system of coordinates | SOI — **unknown** | **Yes** for fixed/constrained |
| Coordinate reference point (marker/ARP/APC) | SOI — **unknown** | **Yes** for fixed/constrained |
| ANTEX/antenna model used by SOI to derive coordinates | SOI | No (diagnostic) |
| Monument type, mounting (rooftop/pillar), obstruction/multipath notes | SOI / site photos | No |
| RINEX sampling rate, file span, data access route and terms | SOI | Yes |
| Meteorological sensors (P, T, RH), heights, calibration | SOI | No (fallback ERA5) |
| BLQ file (FES2014b, CMC) | Onsala service | **Yes** |
| Equipment change / discontinuity list | SOI | Yes (for velocity and coordinate validation) |

### 29.2 Datasets

| Phase | Stations | Period | Purpose |
|---|---|---|---|
| **Pilot (Stage 0–5)** | 3: one IGS station in India with IGS troposphere products (e.g., HYDE or IISC — VERIFY data and product availability) + two SOI CORS (one inland plateau, one coastal or IGP) | 30 consecutive days, recent year with CODE Final + ERA5 available; SHOULD include monsoon or post-monsoon conditions | Algorithm development, PRIDE harmonised benchmark |
| **Validation (Stage 5–7)** | 10–15: ≥ 2 IGS stations; ≥ 3 near IMD radiosonde sites; coastal, IGP, Deccan, Western Ghats, Himalaya, north-east | ≥ 1 full year | Level-1/2/3 gates; tuning experiments |
| **Network (Stage 8–9)** | All SOI stations passing § 10.3 | ≥ 1 full year | Mapping experiment, supportable-resolution map |

---

## 30. Validation and acceptance gates (System 1)

General rules: PRIDE-PPPAR is a **reference implementation, not truth**. Independent references (IGS troposphere products, radiosondes, microwave radiometers, ERA5 as secondary) decide scientific correctness. Metrics (all levels as applicable): bias, STD, RMSE, MAE, Pearson r, P95 |error|, convergence time, AR fixing rate, normalised error z = Δ/σ_comb, uncertainty calibration (share of |z| ≤ 1 and ≤ 2; STD of z).

### 30.1 Level 0 — component tests

| Test | Dataset | Threshold | Pass/Fail |
|---|---|---|---|
| Time conversions | Test vectors incl. leap seconds, week rollover | Exact | All pass |
| SP3 interpolation | Leave-one-node-out, all GPS sats, 7 days | RMS < 1 mm, max < 5 mm | Both |
| CLK handling | Epoch matching, interpolation at node | Exact at nodes | — |
| Solid tide | IERS routine outputs, 10 stations × 100 epochs | < 0.1 mm | Max |
| Ocean loading | HARDISP outputs, test BLQ | < 0.1 mm | Max |
| Pole tide | Reference formula | < 0.05 mm | Max |
| Relativity, Shapiro, Sagnac | Independent implementation (RTKLIB) | < 0.1 mm | Max |
| Satellite/receiver PCO/PCV | RTKLIB + synthetic antenna with asymmetric PCV | < 0.5 mm | Max |
| Wind-up | RTKLIB `windupcorr` | < 0.001 cycle | Max |
| Frame transformation | ITRF examples; equivalence of propagation orders | < 0.1 mm | Max |
| Full modelled range | Term-by-term decomposition vs independent implementation for 1 station-day | < 1 mm per term, total < 2 mm | Max |
| Kalman/RTS numerics | Synthetic data with known truth (simulated ZWD random walk, clock, ambiguities) | Recovered ZWD RMSE consistent with formal σ (z-STD 0.9–1.1) | Both |
| Required logs | Per-test report with inputs, hashes, max/RMS error | | |

### 30.2 Level 1 — harmonised PRIDE benchmark

| Element | Specification |
|---|---|
| Dataset | Pilot set (30 days × 3 stations); later validation set |
| Harmonisation | Same RINEX, same processing coordinates (PRIDE F mode with the same XYZ file), same CODE product files (VERIFY PRIDE accepts them), same ANTEX, BLQ, cutoff, mapping function (both VMF3 or both GMF), ZTD stochastic model (PRIDE STO with the same noise), gradients, sampling |
| Compare | ZTD, ZWD, ZHD, σ_ZTD, receiver clock (after removing mean/trend), WL/NL integers (AR stage), convergence time, rejection and slip counts, continuity, PWV (same met inputs) |
| Metrics | Bias, STD, RMSE, MAE, r, P95; diurnal composite amplitude of Δ; day-boundary jump statistics |
| **Stage 3 gate** | |bias| ≤ 3 mm, STD ≤ 5 mm, RMSE ≤ 6 mm, P95 ≤ 10 mm [A] |
| **Stage 5 gate** | |bias| ≤ 1.5 mm, STD ≤ 3 mm, P95 ≤ 6 mm; 24-h and 12-h harmonic amplitude of Δ ≤ 1 mm; median day-boundary jump ≤ 2 mm [A] |
| Concerning (investigate regardless) | station- or season-dependent bias > 3 mm; Δ correlated with PWV, elevation cutoff, or satellite count |
| Plots | Overlaid time series + Δ; scatter with 1:1 and regression; histogram; monthly bias/STD; diurnal composite; σ ratio; convergence curves; slip/rejection counts |
| Logs | Both configurations, product hashes, PRIDE version, run logs |

### 30.3 Level 2 — independent atmospheric references

| Reference | Use | Gate [A/E] |
|---|---|---|
| IGS final troposphere product (IGS stations in India) | ZTD at hourly matches | |bias| ≤ 3 mm, STD ≤ 6 mm; and custom RMSE ≤ 1.1 × PRIDE RMSE vs same reference |
| CODE troposphere SINEX (validation only) | ZTD at CODE network stations | Reported |
| Radiosonde (IMD; nearest-in-time launches, height-corrected) | PWV | Custom RMSE within ±0.3 mm of PRIDE RMSE; bias reported per season |
| Microwave radiometer (where available) | PWV, high temporal resolution | Reported; used for temporal-resolution assessment |
| ERA5 | ZTD/PWV all stations, all conditions | Reported; not a pass/fail truth |
| Uncertainty calibration | z = Δ/√(σ²_custom + σ²_ref) vs IGS tropo/radiosonde | Derive scale factor s_ZTD; after scaling, 60–76 % of |z| ≤ 1 |

### 30.4 Level 3 — stress and sensitivity tests

| Test | Metric | Gate |
|---|---|---|
| Coordinate perturbations (§ 10.5) | β, horizontal coupling | Completed; β ∈ [0.1, 0.6] or investigated |
| Cutoff 7/10/15° | ΔZTD mean/STD | Reported; chosen cutoff justified |
| Mapping functions VMF3/VMF1/GMF/GPT3 | ΔZTD, RMSE vs references | VMF3 not worse than alternatives |
| ZWD/gradient noise (§ 6.7) | RMSE, NIS | Decision recorded |
| Gradients on/off | RMSE vs references | Reported |
| Injected cycle slips (1, 2, 5, 10 cycles on L1/L2 combinations; random epochs) | Detection rate, false-alarm rate | ≥ 99 % detection for MW-detectable slips; false alarms ≤ 1 per satellite-day in quiet ionosphere [A] |
| Data gaps (5 min – 2 h) | ZTD recovery | No bias > 2 mm after recovery [A] |
| Final vs Rapid products | ΔZTD | STD ≤ 2 mm [E]; reported |
| Eclipse seasons, storm days, cyclone/heavy-rain days | RMSE vs references, rejection rates | Reported; no failures unexplained |
| Radome-fallback simulation | ΔZTD | Quantifies σ inflation |

### 30.5 PPP-AR gates (Stage 6)

WL fixing rate ≥ 90 %, NL ≥ 80 % of candidates [A — compare to PRIDE on same data]; ≥ 95 % agreement of fixed SD integers with PRIDE where comparable; ZTD_fixed − ZTD_float |mean| ≤ 1 mm; no degradation vs Level-2 references; zero published windows failing § 15.6.

### 30.6 PWV gates (Stage 7)

Station PWV vs radiosonde: bias and RMSE reported per season; RMSE not worse than PRIDE-based PWV (same met inputs) by > 0.3 mm; σ_PWV,total calibration as in § 30.3.

---

## 31. Spatial validation experiment (System 2, formal module)

| Element | Specification |
|---|---|
| Data | All stations passing § 10.3 with ≥ 80 % completeness for ≥ 1 year; CODE Final; ERA5 background |
| CV designs | (1) Leave-one-station-out (LOSO); (2) spatial block CV (blocks 100 km and 200 km); (3) network thinning to target spacings **20, 30, 50, 75, 100, 150 km** (random and spatially stratified, ≥ 10 realisations each) |
| Methods | Background only; IDW (reference only); OK on height-normalised PWV; baseline regression kriging (with/without covariates); GPR; space-time OI |
| Grids | **0.5°, 0.25°, 0.125°, 0.1°** (via point reconstruction and cell-representativeness term) |
| Temporal | **5, 15, 30, 60 min** maps vs withheld stations' 5-min series; test whether consecutive 5-min maps differ by more than their σ |
| Regimes | **Winter (JF), pre-monsoon (MAM), SW monsoon (JJAS), post-monsoon (OND), dry (PWV < 15 mm), heavy rainfall (IMD heavy-rain days), cyclones, convective events** (≥ 10 events for extreme categories) |
| Stratification | Region (terrain zones: Himalaya, IGP, Deccan, Western Ghats, coast, north-east, Thar), station density, nearest-station distance, |ΔH|, PWV magnitude, regime |
| Metrics | RMSE, MAE, bias, r, P95; skill vs background; error vs density, distance, |ΔH|, PWV, regime; σ calibration (coverage of ±1σ, ±2σ) |
| Deliverable relationships | E(d, regime, terrain) fits; d* (support limit); Δ_sup map per season (§ 24); calibrated coverage-class thresholds; temporal-resolution decision |
| Acceptance | (i) σ calibration: 60–76 % within ±1σ in classes 1–3; (ii) skill > 0 in classes 1–3 in all regimes; (iii) 0.25° national retained only if classes 1–2 cover a reported fraction of land and meet E_target; otherwise the national primary grid becomes 0.5° — all decisions recorded |
| Plots | Error vs distance per regime; skill maps; supportable-resolution map; reliability diagrams; example maps with confidence overlays; temporal-resolution curves |

---

## 32. Performance: capability, target, acceptance, measured

| Quantity | Theoretical capability [L] | Engineering target [A/E] | Acceptance criterion | Measured |
|---|---|---|---|---|
| Station ZTD vs independent refs (hourly) | ~3–5 mm (IGS tropo ~4 mm) | 4–6 mm RMS | § 30.3 | TBD |
| Station ZTD 5-min | — | 5–10 mm RMS | Reported | TBD |
| Custom vs PRIDE (harmonised) | ~0–2 mm STD for identical models | |bias| ≤ 1 mm, STD ≤ 2–3 mm | § 30.2 gates | TBD |
| Station PWV | ~1–2 mm (normal), 2–3 mm (convective) | 1–2 mm | § 30.6 | TBD |
| National PWV, class 1 | — | ~1.5–2.5 mm **(hypothesis)** | § 31 | TBD |
| National PWV, class 2 | — | ~2–3.5 mm (hypothesis) | § 31 | TBD |
| National PWV, class 3 | background-limited | ~3–5 mm (hypothesis) | § 31 | TBD |
| Convective/mountain cells | — | larger (× 1.5–2) (hypothesis) | § 31 | TBD |

Only the "Measured" column may be quoted as performance.

---

## 33. Software architecture — Python only, deliberately simple (Revision 1.1)

### 33.1 Language and libraries [frozen]

- **Python ≥ 3.10 only.** No C++, no compiled extensions, no pybind11.
- Required: `numpy`, `scipy`, `pandas`, `requests`, `netCDF4` (or `xarray`), `matplotlib`, `hatanaka` (Hatanaka/CRX decompression).
- Optional: `numba` (speed, only if profiling shows need), `cdsapi` (ERA5, only if the user has a CDS key), `pytest` (tests).
- Nothing else without an ADR. No web frameworks, databases, message queues, async code, GUI or plugin systems.

Performance expectation [E]: the state has ~15–40 parameters, so a vectorised NumPy implementation needs of the order of 0.5–3 min per 24-h, 30-s station-day on one core. That is adequate for single stations; the 1,000-station network is handled by running independent station-days in parallel processes (§ 28).

### 33.2 Code structure (flat)

```
pwv_ppp/
├── pwv_ppp.py            # THE entry point for System 1: RINEX in → ZTD/PWV out
├── pwv_map.py            # Entry point for System 2 (Stage 8 only)
├── ppp/                  # One flat package, one module per topic, plain functions + dataclasses
│   ├── settings.py       # ALL defaults, thresholds, URLs/filename templates, constants — one place, commented, with units and references
│   ├── log.py            # Console + file logging setup
│   ├── timesys.py        # GPST/UTC/TT, leap seconds, GPS week/DOY
│   ├── rinex.py          # RINEX 3 obs + nav reading, observable selection
│   ├── products.py       # CODE product search, download, cache, verification, status codes, manifest
│   ├── orbclk.py         # SP3, CLK, ERP, Bias-SINEX, ORBEX readers + interpolation
│   ├── antenna.py        # ANTEX reading, PCO/PCV
│   ├── coords.py         # ITRF2008@2005.0 → processing frame, defaults (§ 9.11), auto check (§ 9.12)
│   ├── corrections.py    # Sagnac, relativity, Shapiro, wind-up, attitude/eclipse, solid tide, ocean loading, pole tide
│   ├── troposphere.py    # VMF3/GPT3, mapping functions, Saastamoinen, gradients
│   ├── preprocess.py     # Arcs, MW/GF slips, clock jumps, ROTI, screening
│   ├── estimator.py      # EKF (Joseph) + RTS smoother + IGG-III + residual editing
│   ├── ambiguity.py      # PPP-AR (Stage 6 only)
│   ├── pwv.py            # ZWD → PWV and uncertainty
│   ├── output.py         # CSV, NetCDF, SINEX_TRO (optional), quick-look plot, summary
│   └── mapping.py        # System 2 (Stage 8 only)
├── tests/                # pytest: small, one test file per module
├── blq/                  # Optional user-supplied ocean-loading files (STATION.blq)
├── stations.csv          # Optional: station, X, Y, Z (ITRF2008 @ 2005.0)
├── cache/                # Auto-created product cache
└── results/              # Auto-created outputs
```

Simplicity rules for implementation [binding]:
- Plain functions and `@dataclass` containers; no inheritance hierarchies, abstract base classes, factories, registries of plugins or dependency injection.
- One settings module instead of a configuration system; an optional small YAML file MAY override `settings.py` values, nothing more.
- A module should stay readable (guideline: ≤ ~800 lines); split only when it clearly helps.
- The parameter registry of § 5.1 is a simple ordered list of state names with index lookup — enough to add Galileo later without redesign.
- Provenance is a JSON file written next to the outputs (no database). The product cache is plain files with a `.meta.json` sidecar (hash, source URL, time).
- Every scientific formula is a small, pure, separately testable function that cites its reference in its docstring.

### 33.3 Command-line interface (minimal)

```
python pwv_ppp.py SITE.o SITE.n
python pwv_ppp.py SITE.o SITE.n --xyz X Y Z
```

| Argument | Required | Default | Meaning |
|---|---|---|---|
| `SITE.o` (RINEX 3 obs; `.crx`/`.gz` accepted) | **Yes** | — | Observations |
| `SITE.n` (RINEX nav) | Yes (as requested by the user) | — | Used for data sanity checks, satellite health and the initial single-point solution; precise CODE products are always downloaded for the PPP itself |
| `--xyz X Y Z` | No | from `stations.csv` if present | SOI coordinate in metres, ITRF2008 @ 2005.0 |
| `--mode` | No | `auto` | `auto`, `fixed`, `constrained`, `static` (§ 9.12) |
| `--out DIR` | No | `./results` | Output folder |
| `--cutoff DEG` | No | 7 | Elevation cutoff |
| `--ar` | No | off | PPP-AR (only after Stage 6) |
| `--no-pwv` | No | off | Skip PWV conversion |
| `--offline` | No | off | Use cache only |
| `--proxy URL` | No | environment `HTTPS_PROXY` | HTTP(S) proxy |
| `--quiet` / `--verbose` | No | normal | Log detail on screen |

Everything else (products, ANTEX, VMF3, GPT3, leap seconds) is found, downloaded, verified and cached automatically. Multiple days or stations are handled by calling the command repeatedly (a trivial shell loop or `xargs -P` for parallel runs) — no scheduler is built into the program.

### 33.4 Logging (visible progress)

- Python `logging`, same messages to the screen and to `results/<STATION>_<YYYYDDD>.log`.
- One INFO line per stage: RINEX summary (station, time span, interval, epochs, GPS satellites, selected observables); products (family/tier tried, files found/downloaded/verified, final status); metadata assumptions (each WARNING of § 9.11); coordinate check (Δ N/E/U, chosen mode); preprocessing (arcs, slips by cause, rejections); filter progress (every 10 % of epochs: epoch, satellites, current ZTD ± σ); smoother and residual-editing passes; output files written; final summary (number of 5-min values, converged fraction, mean ZTD/PWV, flags raised, run time).
- WARNING for every assumption or fallback; ERROR only for conditions that end the run, always with a plain-language explanation and the status code.

### 33.5 Outputs per run (in `--out`)

`<STATION>_<YYYYDDD>_ZTD.csv`, `<STATION>_<YYYYDDD>_ZTD.nc`, `<STATION>_<YYYYDDD>_PWV.csv` (unless `--no-pwv`), `<STATION>_<YYYYDDD>_manifest.json`, `<STATION>_<YYYYDDD>.log`, `<STATION>_<YYYYDDD>_quicklook.png` (ZTD, PWV, number of satellites, flags), and optionally `<STATION>_<YYYYDDD>.tro` (SINEX_TRO).

### 33.6 Automatic downloads (anonymous sources preferred; all URLs/templates in `settings.py`, VERIFY DURING IMPLEMENTATION)

| Item | Primary source | Notes |
|---|---|---|
| CODE SP3/CLK/ERP/OSB/ORBEX | AIUB `http://ftp.aiub.unibe.ch/CODE/{yyyy}/` (anonymous) | Mirror templates may be added; CDDIS only if the user sets up Earthdata `.netrc` |
| ANTEX (igs20) | IGS files server (`files.igs.org`, station/general) | Version must match CODE clock header |
| VMF3 grids | TU Wien VMF data server, `trop_products/GRID/1x1/VMF3/VMF3_OP/{yyyy}/` | 6-hourly files |
| GPT3 grid | TU Wien VMF data server (`codes/`) | One-time download |
| Leap seconds | IERS / IETF leap-second file | One-time, refreshed monthly |
| BLQ | **Not downloadable automatically** | User places `blq/<STATION>.blq` once; otherwise `NO_OTL` warning |

---

## 34. Repository structure

As in § 33.2, plus:

```
docs/
├── TDS.md                # This document
├── decisions.md          # One short entry per changed decision (date, reason, evidence)
└── verify_log.md         # Outcome of every VERIFY DURING IMPLEMENTATION item, with source
CLAUDE.md                 # Short rules for Claude Code (§ 41)
requirements.txt
```

Validation material (PRIDE runs, reference data, reports) lives in `validation/` with one subfolder per gate; it is not part of the program.

---

## 35. Staged roadmap with gates

| Stage | Modules | Tests | Gate | MUST NOT yet |
|---|---|---|---|---|
| **0 Reference data & PRIDE baseline** | none (data collection) | — | Pilot RINEX (.o/.n) for 3 stations × 30 days available; SOI coordinates collected; PRIDE runs with CODE products archived | Write PPP code |
| **1 I/O, products, time** | `settings`, `log`, `timesys`, `rinex`, `orbclk`, `antenna`, `products`, CLI skeleton that runs end-to-end up to "products ready" with logs | Level-0 time/reader/SP3 tests; downloader scenarios (offline, proxy, corrupt file, missing tier → WAIT/FAIL) | 30 pilot days: products found, verified, cached, statuses correct, logs readable | Any estimation |
| **2 Deterministic corrections** | `coords`, `corrections`, `troposphere` | All § 30.1 tests; frame transformation tests; § 9.11 defaults exercised | Level-0 thresholds met | Estimation |
| **3 Float PPP-ZTD** | `preprocess` (basic), `estimator`, `output` | Synthetic-truth test; Level-1 PRIDE benchmark | Stage-3 gate (§ 30.2) | AR, mapping |
| **4 Coordinate modes & auto check** | `auto`/`fixed`/`constrained`/`static`; § 9.12 | Perturbation experiment (§ 10.5); auto-check logic tests | β measured; auto mode behaves as specified on pilot stations | Claims about SOI metadata |
| **5 QC/robustness** | full `preprocess`, eclipse handling, overlap windows, tuning (§ 6.7) | Level-3; Level-1 Stage-5 gate; Level-2 | Gates passed on the 10–15-station year set | AR |
| **6 PPP-AR** | `ambiguity` | § 30.5 | § 30.5 | Replacing float products |
| **7 PWV** | `pwv` | § 30.6 | § 30.6 | Mapping with unvalidated PWV |
| **8 National mapping** | `mapping`, `pwv_map.py` | § 31 | § 31 | Resolution claims before § 31 |
| **9 Network processing** | parallel shell/`xargs` recipe, reprocessing script | Reproducibility, throughput | Full network ≥ 1 yr reproducible | Operational promises |
| **10 Operationalisation** | daily cron recipe, monitoring of logs | 30-day stable run | Stable daily RAPID/FINAL | NRT promises (< ~18 h) |

---

## 36. Hard engineering and scientific prohibitions

1. Do not mix CODE with WUM, GFZ, CNES/CLS, IGS combined or any other AC product in a CODE-only run.
2. Do not mix CODE families (COD0OPSFIN / CODMOPSRAP / COD0OPSRAP; OPS / MGX) within a processing window.
3. Do not use CODE Ultra-Rapid as a PPP product.
4. Do not silently substitute missing products or silently change product tiers.
5. Do not hard-code product directory structures, filenames, URLs or frame labels.
6. Do not hard-code receiver observation codes or vendor conventions.
7. Do not use raw SOI ITRF2008@2005.0 coordinates as processing coordinates; do not mix coordinate frames, epochs, tide systems or reference points.
8. Do not fix a coordinate that the same-day automatic check (§ 9.12) contradicts, unless the user forces `--mode fixed` (then warn).
9. Do not apply the ARP eccentricity twice or ignore the radome.
10. Do not interpolate raw ZTD (or raw PWV) across India without height treatment.
11. Do not call grid spacing "atmospheric resolution"; do not present a fine raster as a fine atmospheric product.
12. Do not fill unsupported regions without a confidence class; class-4 cells carry fill values in GNSS-corrected layers.
13. Do not treat PRIDE (or ERA5) as absolute truth.
14. Do not claim any accuracy before the corresponding gate has produced a measured result.
15. Do not allow PPP-AR results to overwrite or contaminate the validated float products.
16. Do not present 5-min station values as independent 5-min observations or 5-min national maps as an operational product without § 31 evidence.
17. Do not hide product, configuration or software provenance; every output carries its manifest ID and configuration hash.
18. Do not invent constants, transformation parameters or undocumented PRIDE behaviour; mark and resolve VERIFY items.
19. Do not promise near-real-time (< ~18 h) products under the CODE-only policy.
20. Do not stop a run because station metadata are unknown; apply § 9.11 defaults, warn and flag.
21. Do not introduce C++ or any other language, frameworks, databases or complex class hierarchies (§ 33).
22. Do not require inputs beyond `SITE.o SITE.n` for a standard run; everything else is optional or automatic.

---

## 37. Frozen scientific decisions

| # | Decision | Value |
|---|---|---|
| F1 | Constellation v1 | GPS only; architecture multi-GNSS-ready |
| F2 | Observation model | Undifferenced IF code + phase, L1/L2 |
| F3 | Estimator | EKF (Joseph form) + RTS smoother; smoothed solution is primary |
| F4 | Troposphere state | ZWD random walk at observation rate; ZHD₀ fixed; ZTD = ZHD₀ + ZWD |
| F5 | Mapping function | VMF3 (GPT3 fallback, flagged); Chen & Herring gradients, estimated by default |
| F6 | Cutoff | 7° default (tested 7/10/15°) |
| F7 | Output sampling | 5-min UTC marks, interpolated from smoothed states |
| F8 | Coordinate modes | `auto` (default) / fixed / constrained / static; `auto` fixes the SOI coordinate only when the same-day check confirms it (§ 9.12) |
| F9 | SOI coordinates | ITRF2008@2005.0 → propagate (9.1) → transform (9.3, official IGN parameters) → product frame; never raw |
| F10 | Product policy | CODE only; COD0OPSFIN → CODMOPSRAP → COD0OPSRAP; no ULT for PPP; status codes § 11.4 |
| F11 | Slips | Detect (LLI, MW, GF, gaps) and reset; no repair in v1 |
| F12 | PPP-AR | Separate later layer; hard family consistency; never overwrites float |
| F13 | PWV | Separate module; Bevis 1994 constants default; P and T_m sources ranked (§ 18.2) |
| F14 | Mapping formulation | NWP background at height + GNSS residual regression kriging; no raw ZTD interpolation |
| F15 | National product targets | 0.25° national, 0.125° regional where supported, 0.5° coarse; 15-min maps; hourly archive; 5-min station retained — **initial targets subject to § 31** |
| F16 | Formats | NetCDF-4 CF/ACDD archive; Zarr analysis store; COG distribution; station CSV/NetCDF (+ optional SINEX_TRO) |
| F17 | Language | **Python only** (NumPy/SciPy), flat simple structure (§ 33) |
| F18 | Units | mm for delays/PWV in all products |
| F19 | Unknown metadata | Documented defaults (§ 9.11), warnings, flags; never stop |
| F20 | User interface | `python pwv_ppp.py SITE.o SITE.n` + optional flags; automatic downloads; on-screen logs |

Changing a frozen decision requires an ADR in `docs/decisions/` with evidence.

---

## 38. Open scientific decisions (require experiments)

| # | Question | Resolved by |
|---|---|---|
| O1 | σ_ZWD and σ_G values (per season/region?) | § 6.7 |
| O2 | Processor-specific height-to-ZTD coupling β | § 10.5 |
| O3 | SOI velocities, tide system, reference point, antenna model era (defaults in use, § 9.11) | Accumulated § 9.12 discrepancies over many days; SOI information if it becomes available |
| O4 | Stations off the rigid India plate: velocity treatment | Own time series; § 9.4 |
| O5 | Code/phase noise and IGG-III thresholds | Level-1/3 tests |
| O6 | Slip thresholds under high ionospheric activity | § 30.4 injected-slip + storm-day tests |
| O7 | Whether AR measurably improves ZTD for this network | § 30.5 vs Level-2 |
| O8 | Best T_m source for India | Radiosonde comparison (§ 30.6) |
| O9 | Covariance model, regions, anisotropy | § 31 |
| O10 | Coverage-class thresholds and d* | § 31 |
| O11 | Supportable resolution by region/season | § 31 → § 24.4 |
| O12 | Whether 5-min or 30-min national maps are justified anywhere | § 31 temporal tests |
| O13 | Operational NWP background source and its GNSS assimilation status | Stage 8 |
| O14 | ZTD formal-error scaling factor s_ZTD | § 30.3 |
| O15 | Ambiguity continuity across CODE product days | Level-0 test (§ 11.2) |

---

## 39. Required external references

**Standards and formats:** IERS Conventions 2010 (IERS TN 36) and official updates (Ch. 5, 7, 9, 10, 11); IS-GPS-200 (current); RINEX 3.05 and 4.0x; SP3-d; Clock RINEX 3.04; ANTEX 1.4; Bias-SINEX 1.00; ORBEX; SINEX_TRO 2.00; CF Conventions 1.10; ACDD 1.3.
**Reference frames:** ITRF2008, ITRF2014, ITRF2020 documentation and official transformation files (IGN: `Transfo-ITRF2020_TRFs.txt`, `Transfo-ITRF2014_ITRFs.txt`); ITRF2020-u2023/u2024; IGSMAIL-8238 (IGS20), IGSMAIL-8543 (IGb20), IGSMAIL-8634 (IGc20); Altamimi et al. (2023) ITRF2020 Plate Motion Model; Blewitt & Lavallée (2002) on velocity estimation.
**CODE/IGS products:** CODE IGS Technical Report 2024 (product table); CODE `IAR_README.TXT`; Schaer et al. (2021) J. Geod. 95; IGS long-filename guideline v2.0; Kouba (2015) Guide to Using IGS Products.
**PPP/AR:** Zumberge et al. (1997); Kouba & Héroux (2001); Wu et al. (1993); Blewitt (1989, 1990); Ge et al. (2008); Laurichesse et al. (2009); Collins et al. (2010); Geng et al. (2012, 2019); Teunissen (1995, 1998); Verhagen & Teunissen (2013); Dong & Bock (1989); Kouba (2009) eclipse yaw; Yang et al. (IGG-III robust estimation).
**Troposphere/meteorology:** Saastamoinen (1972); Davis et al. (1985); Askne & Nordius (1987); Chen & Herring (1997); Boehm et al. (2006); Kouba (2008) gridded VMF; Landskron & Böhm (2018) VMF3/GPT3; Bevis et al. (1992, 1994); Ning et al. (2016); Stępniak et al. (2022); Santerre (1991).
**Spatial:** Emardson et al. (2003); Xu et al. (2011); Cressie (1993) kriging; Rasmussen & Williams (2006) GPR.
**Reference software (consult, do not copy):** IERS routines (DEHANTTIDEINEL, HARDISP), SOFA, TU Wien VMF3/GPT3 code, RTKLIB (BSD-2), PRIDE-PPPAR (GPLv3 — algorithm reference and benchmark only).

---

## 40. Acceptance checklist

**Stage 0**
- [ ] Pilot stations selected; RINEX .o/.n and SOI ITRF2008@2005.0 coordinates collected (other metadata: § 9.11 defaults)
- [ ] BLQ files obtained where possible (otherwise `NO_OTL` runs accepted); ANTEX entries checked
- [ ] PRIDE harmonised runs with CODE products archived

**Stage 1**
- [ ] All readers pass Level-0 tests
- [ ] Product manager produces OK_AR / OK_FLOAT_ONLY / WAIT_FOR_TIER / FAIL on test scenarios
- [ ] Offline mode, proxy, retries, atomic writes, verification, manifests working
- [ ] No hard-coded filenames/paths/URLs (config-driven, reviewed)

**Stage 2**
- [ ] Every correction in § 7 passes its Level-0 threshold
- [ ] Coordinate module: equations (9.1)–(9.4) implemented table-driven; equivalence test < 0.1 mm; frame labels read from products
- [ ] Modelled-range decomposition agrees with independent implementation

**Stage 3**
- [ ] Synthetic-truth filter/smoother test passes
- [ ] 5-min UTC extraction and flags implemented
- [ ] Stage-3 PRIDE gate passed

**Stage 4**
- [ ] `auto` mode check (§ 9.12) verified on pilot stations; optional C1–C8 diagnostics reported
- [ ] Perturbation experiment done; β per station stored

**Stage 5**
- [ ] Full QC implemented; tuning experiments decided and recorded
- [ ] Stage-5 PRIDE gate and Level-2 gates passed; Level-3 report complete

**Stage 6**
- [ ] AR gates (§ 30.5) passed; float layer untouched

**Stage 7**
- [ ] PWV gates (§ 30.6) passed; σ components reported separately

**Stage 8**
- [ ] § 31 experiment complete; supportable-resolution map and class thresholds published
- [ ] σ calibration and skill criteria met

**Stage 9–10**
- [ ] Reproducibility (identical hashes) and throughput demonstrated
- [ ] 30-day stable RAPID/FINAL operation with reprocessing/supersession working

---

## 41. HANDOFF TO CLAUDE CODE

**Read first:** § 0 (incl. § 0.5 Revision 1.1), § 9.11–9.12, § 33, § 35, § 36, § 37 of this document. The two earlier feasibility documents are background only.

**Binding instructions:**
1. **Python only.** No C++ or other languages (§ 33.1).
2. **Keep it simple.** Flat package exactly as § 33.2; plain functions and dataclasses; one `settings.py`; no frameworks, no databases, no class hierarchies, no plugin systems. Prefer clear, commented scientific code over clever abstractions.
3. **One command:** `python pwv_ppp.py SITE.o SITE.n` must work with no other input. Optional flags only as in § 33.3.
4. **Download everything automatically** (CODE products, ANTEX, VMF3, GPT3, leap seconds), cache it, verify it, and never mix CODE families or use other Analysis Centres.
5. **Never stop on unknown metadata.** Apply § 9.11 defaults, print a WARNING, flag outputs and continue. End a run only when PPP is impossible (no dual-frequency GPS data, no CODE orbit/clock) — cleanly, with a clear message and status code.
6. **Show logs** on screen and in a log file as in § 33.4.

**Create first:** `CLAUDE.md` (these six rules + § 36 prohibitions + current stage), `requirements.txt`, the § 33.2 skeleton, `settings.py` with the defaults of this TDS (units and references in comments), `docs/verify_log.md`.

**Start with Stage 1** (Stage 0 is data collection by the project owner): `timesys`, `rinex`, `orbclk`, `antenna`, `products`, `log`, and a `pwv_ppp.py` that runs end-to-end up to "all products ready", with Level-0 tests alongside.

**Do NOT, before the Stage 1 and 2 gates are passed:** write the estimator, smoother, PPP-AR, PWV or mapping code; use raw ITRF2008@2005.0 coordinates in any geometry; tune thresholds on validation data; report any accuracy number; guess any VERIFY item (resolve from the primary source and record it in `docs/verify_log.md`).

**Working rules:** each module ships with small pytest tests; each run writes a manifest; each gate produces a short report in `validation/`; deviations from this TDS are recorded in `docs/decisions.md`.

---

## Appendix A — Sources consulted for this TDS (current documentation)

- CODE Analysis Center, IGS Technical Report 2024 (product table): http://ftp.aiub.unibe.ch/CODE/CODE_IGSTR.pdf
- CODE IAR README: http://ftp.aiub.unibe.ch/CODE/IAR_README.TXT
- AIUB CODE Analysis Center page: https://www.aiub.unibe.ch/research/code___analysis_center/index_eng.html
- Schaer et al. (2021): https://link.springer.com/article/10.1007/s00190-021-01521-9
- IGN transformation parameters ITRF2020 → past ITRFs: https://itrf.ign.fr/docs/solutions/itrf2020/Transfo-ITRF2020_TRFs.txt
- IGN transformation parameters ITRF2014 → past ITRFs: https://itrf.ign.fr/docs/solutions/itrf2014/Transfo-ITRF2014_ITRFs.txt
- ITRF2020 solution page (PMM, post-seismic models): https://itrf.ign.fr/en/solutions/ITRF2020
- IGSMAIL-8543 (IGb20): https://lists.igs.org/pipermail/igsmail/2024/008539.html
- IGSMAIL-8634 (IGc20): https://lists.igs.org/pipermail/igsmail/2025/008630.html
- IGS long-filename guideline v2.0: https://files.igs.org/pub/resource/guidelines/Guideline_for_the_transition_of_the_IGS_products_to_IGS20_and_long_filenames_v2.0.pdf
- PRIDE PPP-AR II manual and repository (reference behaviour cited only where documented)
