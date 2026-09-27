# GPS PPP → ZTD → PWV Processor: Feasibility Assessment

*Prepared 21 Sep 2026 · Scope: SOI CORS, GPS-only, fixed/constrained coordinates, 5-min ZTD · No code written yet*

---

## 0. Short answers to your eight questions

| # | Question | Short answer |
|---|---|---|
| 1 | Can it be built from literature + IGS standards + PRIDE as reference? | **Yes.** Every model needed is published (IERS 2010 Conventions, IGS formats, ANTEX, VMF3/GPT3, LAMBDA, PPP-AR papers). No proprietary knowledge is required. |
| 2 | Easy / hard / accuracy-deciding parts | Easy: file I/O, orbit/clock interpolation, most deterministic corrections. Hard: preprocessing/QC (especially at Indian low-latitude ionosphere), estimator bookkeeping, PPP-AR. **Accuracy-deciding:** antenna model + site metadata, mapping function/elevation weighting, product consistency, outlier handling, ZWD stochastic model, correctness of fixed coordinates. AR matters least. |
| 3 | Expected ZTD accuracy vs PRIDE | After validation: custom − PRIDE ≈ **bias ≤ 1–3 mm, STD 2–5 mm** (harmonised settings). Against independent truth, both should land around **4–8 mm RMS** for hourly-scale ZTD. Not demonstrated until tested. |
| 4 | Implied PWV | ≈ 0.16 × ZTD error → **~0.6–1.3 mm PWV** from ZTD alone; realistic total **1–2 mm** (up to ~3 mm in monsoon/convective conditions or without a barometer). |
| 5 | Do known coordinates help? | **Yes, modestly for a 24-h static solution, strongly for short windows, convergence and re-convergence** — *if* the coordinates are correct in the same frame, epoch, tide system and antenna model as the processing. |
| 6 | Reasons not to fix / bias risk | **Yes.** A height error δh leaks into ZTD as roughly 0.2–0.4·δh (cutoff-dependent; fixing too high a height gives a positive ZTD bias). Frame/epoch mismatch (India moves ~5 cm/yr), ANTEX mismatch, ARP/marker confusion, tide-system mismatch and seasonal loading are all realistic sources of cm-level error. |
| 7 | Absolutely required inputs | Consistent precise orbit + clock (30-s clocks ideally), ERP, igs20.atx (+ correct receiver antenna/radome), ocean-loading BLQ per station, a good mapping function (VMF3 grids or GPT3), a priori ZHD source (barometer or NWP), code/phase OSB for AR, correct station metadata. |
| 8 | Validation | Model-by-model unit tests → harmonised PRIDE benchmark → independent truth (IGS troposphere product at HYDE/IISC-type stations, ERA5, radiosonde) → coordinate-perturbation and stress tests → ≥10 stations × 1 year. |

---

## Part 1 — How PPP-for-troposphere works (plain-language architecture)

### 1.1 The observation equations

For each satellite *s* and epoch *t*, the receiver gives code (pseudorange, P) and carrier phase (L) on two frequencies. Combining L1 and L2 in the **ionosphere-free (IF) combination** removes ~99.9 % of the ionospheric delay (the first-order term):

```
P_IF = ρ + c·dt_r − c·dt^s + m_h(e)·ZHD + m_w(e)·ZWD + G-terms + b_P + ε_P
L_IF = ρ + c·dt_r − c·dt^s + m_h(e)·ZHD + m_w(e)·ZWD + G-terms + λ·N_IF + wind-up + ε_L
```

- ρ — geometric range between satellite (from the precise orbit, SP3) and receiver antenna.
- dt_r — receiver clock error (unknown), dt^s — satellite clock (from the precise clock file).
- m_h, m_w — hydrostatic and wet mapping functions (≈ 1/sin(elevation)).
- G-terms — optional horizontal gradients (north/east asymmetry of the atmosphere).
- N_IF — carrier-phase ambiguity, constant along an unbroken tracking arc.
- ε — noise, multipath, residual model error.

Everything that can be computed is computed and subtracted: orbit, satellite clock, antenna phase centres, tides, loading, relativity, Sagnac, phase wind-up, a priori ZHD. What remains is a small linear problem.

### 1.2 What is estimated, epoch by epoch (fixed-coordinate mode)

| Parameter | Count | Behaviour in the estimator |
|---|---|---|
| Receiver clock | 1 per epoch | White noise (re-estimated freely every epoch) |
| ZWD (residual zenith delay) | 1 | Random walk (slowly varying) |
| Gradients G_N, G_E | 2 (optional) | Very slow random walk |
| IF ambiguities | 1 per tracked satellite arc (~8–12) | Constant until a cycle slip; reset on slip |
| Coordinates | 0 (fixed) or 3 (static / constrained) | Constant |

So the state is only ~12–16 parameters. The *physics* is not in the estimator; it is in the deterministic model that computes ρ and all the corrections. **That is where most accuracy is won or lost.**

### 1.3 How ZTD is separated from the other unknowns

The estimator separates parameters by how their **partial derivatives vary across satellites and time**:

- **Receiver clock**: same effect (1.0) on every satellite at a given epoch.
- **ZWD**: effect scales with m(e): ≈1 at zenith, ≈3.8 at 15°, ≈7–8 at 7°. The contrast between high and low satellites is what isolates ZTD from the clock. This is why low-elevation data (7°–10° cutoff) matter so much and why mapping-function quality matters so much.
- **Ambiguities**: constant per satellite arc, while geometry changes as the satellite moves; carrier phase alone cannot separate them from other constant offsets until enough geometry change has accumulated. This is the **convergence period**.
- **Station height**: effect ∝ sin(e) — the *opposite* elevation signature to troposphere. Height, clock and ZTD are strongly correlated (typical correlation magnitudes 0.8–0.9 in short windows). **This is the core reason fixing the coordinates helps**: it removes the parameter that competes most directly with ZTD.

### 1.4 How orbit and clock errors reach ZTD

The satellite's radial orbit error and clock error both project directly into the range. A component common to all satellites is absorbed by the receiver clock. A satellite-specific constant part is absorbed by that satellite's ambiguity. What is left — errors that vary along the pass and correlate with elevation — leaks into ZTD. With IGS-quality final products (orbits ~2–3 cm, clocks ~tens of ps, internally consistent) the leak is at the mm level; published comparisons find final-product quality is no longer the dominant ZTD error. Two rules follow: **never mix orbits and clocks from different analysis centres**, and **avoid interpolating 5- or 15-min clocks to 30-s epochs** (process at the clock sampling instead).

### 1.5 How ambiguities are handled

1. Each new or slipped arc gets a new float ambiguity with a large initial variance.
2. The filter/smoother estimates it as real-valued (**float PPP**). For 24-h static tropospheric work, float PPP already gives good ZTD.
3. **PPP-AR** (optional, stage 2): the IF ambiguity is split into wide-lane (λ ≈ 86 cm, fixed from Melbourne–Wübbena averages + satellite WL biases) and narrow-lane (λ ≈ 10.7 cm, fixed with LAMBDA/bootstrapping using satellite phase biases). The phase biases (OSB) must come from the **same analysis centre** as the clocks. Fixed ambiguities are then imposed as constraints and the solution recomputed.
4. After a slip: a new ambiguity is created (or the slip is repaired if the geometry-free and MW tests can determine the integer jump reliably).

### 1.6 Why the fixed coordinates help — and when they don't

- They remove the height–ZTD–clock correlation, so ZTD becomes determinable from fewer epochs → faster convergence, better ZTD after data gaps, more stable 5-min values.
- In a 24-h static solution the daily coordinate is itself determined to a few mm, so the gain for day-long ZTD is **modest** (expect ≲1–3 mm STD change). The big gains are in short windows, near-real-time use and re-convergence.
- The cost: **any coordinate error becomes a ZTD bias**, with no way for the estimator to flag it.

### 1.7 What ultimately limits ZTD accuracy

In rough order of importance for a good static CORS: (1) receiver antenna calibration/radome/ARP metadata and site multipath; (2) mapping function, elevation cutoff and elevation weighting; (3) product quality and consistency; (4) preprocessing and outlier handling; (5) ZWD stochastic model; (6) coordinate error (fixed mode); (7) ambiguity resolution (smallest effect for 24-h static).

---

## A. Feasibility

**Can this processor be built? Yes.** A research-grade GPS float PPP engine for ZTD, with fixed/constrained/estimated coordinates, full IERS-2010-level corrections, Kalman filter + backward smoother, automated IGS product management and 5-minute ZTD output, is a well-defined engineering task. Every component has an authoritative specification:

- Observation model and PPP: Zumberge et al. (1997); Kouba & Héroux (2001); Kouba (2015, *A Guide to Using IGS Products*).
- Corrections: IERS Conventions 2010 (Petit & Luzum), Ch. 7 (solid tide, ocean loading, pole tide), with reference Fortran code; Wu et al. (1993) phase wind-up; ANTEX 1.4 format; Kouba (2009) GPS yaw during eclipse.
- Troposphere: Saastamoinen (1972) / Davis et al. (1985) ZHD; Boehm et al. (2006) VMF1/GMF; Landskron & Böhm (2018) VMF3/GPT3; Chen & Herring (1997), Bar-Sever et al. (1998) gradients.
- GNSS meteorology: Bevis et al. (1992, 1994); Ning et al. (2016) PWV uncertainty budget.
- PPP-AR: Ge et al. (2008); Laurichesse et al. (2009); Collins et al. (2010); Geng et al. (2012, 2019); Teunissen (1995) LAMBDA; Schaer et al. (2021) CODE OSB/integer clocks; Geng et al. (2019) PRIDE PPP-AR.
- Formats: RINEX 3.05/4.0x, SP3-d, Clock RINEX 3.04, Bias-SINEX 1.00, SINEX_TRO 2.00, IGS long-filename convention.

**PRIDE-PPPAR** (GPLv3, Fortran, v3.2.x as of Aug 2026) is an excellent reference and benchmark. Key facts relevant to the design (from its manual and wiki):

- Estimator: batch least squares with parameter elimination/recovery ("process" parameters: clock, ZTD; "state": ambiguities; "constant": coordinates in S/F modes). Receiver clock white noise per epoch.
- Positioning modes: **S** (static), **P** (piecewise), **K** (kinematic), **F** (fixed — coordinates from a SINEX or a user-supplied file), L (LEO).
- ZTD: **STO** (random walk, default noise 0.0004 m/√s ≈ 24 mm/√h) or **PWC** (piecewise constant, default 0.02 m/√h, minimum interval 1 min). Gradients: default PWC:720 (12 h) in S/F modes.
- Mapping functions: GMF (default), NMF, VMF1, VMF3. Default cutoff 7°.
- Corrections: solid, ocean and pole tides; BLQ from the Onsala loading service; ANTEX PCO/PCV; wind-up; optional 2nd-order ionosphere with CODE GIM.
- AR: WL/NL, LAMBDA with ratio/AI validation for short sessions, rounding for >6 h sessions.
- Default products: Wuhan **WUM0MGXRAP** (orbit, clock, ERP, OSB, attitude quaternions).

**Licensing note (important):** PRIDE is GPLv3. If your processor must be distributable under a different licence (or is to be owned by your organisation without GPL obligations), it must be a **clean implementation from the literature**, using PRIDE only for algorithm understanding, configuration comparison and numerical benchmarking — not code translation. That also matches your requirement for independence.

### Two claims, kept separate

| Claim | Assessment |
|---|---|
| **"Claude can implement this."** | **Yes, with high confidence** for float PPP-ZTD, the product manager, fixed/constrained/static coordinate modes and QC. **Yes, with moderate confidence** for PPP-AR (the conventions for biases and integer properties are unforgiving). |
| **"Claude can implement this with accuracy comparable to an established geodetic package."** | **Plausible, not guaranteed, and only after an iterative validation loop on real data.** Writing the code is the smaller part; finding the sign errors, convention mismatches, and QC weaknesses by comparing against PRIDE and independent truth is the larger part. I would expect the first working version to disagree with PRIDE by more than the target; reaching few-mm agreement is a matter of systematic debugging, which I can do but you must run (or give me) the data. |

My specific limitations: I cannot observe your SOI data or run PRIDE unless the data are made available in the working environment; subtle convention errors (sign of a tide term, frame of a PCO, units of a bias) do not announce themselves and are found only by differencing against references; and I have no memory across sessions beyond what is saved in the project, so design decisions must be documented as we go.

---

## B. Scientific scope — what research-grade ZTD requires

Status key: **M** = mandatory in v1 · **R** = strongly recommended in v1 · **L** = later version.

### Observation & preprocessing
| Item | Status | Notes |
|---|---|---|
| Dual-frequency code + phase, IF combination | M | L1/L2. IGS clocks refer to the C1W/C2W (P1/P2) IF combination. Undifferenced-uncombined is a valid alternative but brings nothing for ZTD in GPS-only v1; keep the design open for it. |
| Observable selection by priority (not by vendor) | M | e.g. L1: C1W > C1C; L2: C2W > C2L/C2S/C2X. Record which signals were used. Apply C1C→C1W code bias (from OSB/DCB) when C1W absent. |
| Elevation cutoff | M | Default 7°; configurable (sensitivity tests at 7°, 10°, 15°). |
| Elevation-dependent weighting | M | σ² = a² + b²/sin²e; phase ~3 mm, code ~0.3 m (code/phase ≈ 100). |
| Cycle-slip detection (MW + geometry-free, TurboEdit-type) | M | Blewitt (1990). Thresholds must be robust to strong ionospheric gradients (see risks). |
| Loss-of-lock indicator (LLI) handling | M | From RINEX flags. |
| Post-fit residual screening / robust weighting | M | Iterative residual editing (like PRIDE's) plus IGG-III-type down-weighting. |
| Receiver clock jump (ms) handling | M | Common with some receivers. |
| SNR-based weighting/screening | R | Useful against low-elevation multipath; keep as option. |
| Cycle-slip *repair* | L | Detection + new ambiguity is sufficient for static 24-h ZTD. |

### Orbits, clocks and satellite side
| Item | Status | Notes |
|---|---|---|
| Precise SP3 orbit interpolation | M | Lagrange (order ~10) or Neville; needs adjacent-day files at day boundaries. |
| Precise clock (30 s preferred) | M | Linear interpolation only between dense clock epochs; otherwise process at clock sampling. |
| Satellite health / exclusion | M | Missing/flagged SP3 or clock values, manoeuvres, bad residual arcs. |
| Satellite PCO + PCV (igs20.atx) | M | Frequency-specific, combined into IF. Must match the product's ANTEX. |
| Satellite attitude (nominal yaw) | M | Needed for PCO direction and wind-up. |
| Eclipse-season yaw | R | Either use ORBEX attitude quaternions from the product AC, implement Kouba (2009)-type models, or exclude eclipsing satellites during noon/midnight turns (simplest acceptable v1). |
| Relativity: periodic clock term (−2 r·v/c²) | M | IGS clocks exclude it by convention. |
| Shapiro (gravitational) delay | M | ~cm level, elevation-dependent. |
| Earth rotation during signal travel (Sagnac) | M | Rotate satellite position by ω_E·τ. |
| Phase wind-up | M | Wu et al. (1993). |

### Station side
| Item | Status | Notes |
|---|---|---|
| Receiver antenna PCO + PCV incl. radome | M | **One of the largest ZTD bias sources if wrong.** Confirm SOI antenna/radome codes exist in igs20.atx. |
| ARP eccentricity (RINEX ΔH/E/N) | M | Make sure the SOI coordinates refer to the same point (marker vs ARP). |
| Solid Earth tide (IERS 2010) | M | Up to ~30 cm vertical; essential in fixed mode. |
| Ocean tide loading (BLQ, FES2014b) | M | mm to cm; larger at coastal SOI stations. Per-station BLQ from the Onsala service. |
| Pole tide (solid) | M | Few mm; cheap. Ocean pole tide loading: L. |
| Tide system consistency | M | ITRF/IGS coordinates are conventional tide-free; SOI coordinate convention must be confirmed. |
| Atmospheric tidal / non-tidal loading | L | mm-level; relevant to fixed-coordinate accuracy but can wait. |

### Troposphere
| Item | Status | Notes |
|---|---|---|
| A priori ZHD | M | Saastamoinen with pressure from barometer, or NWP (VMF3 grid ZHD / ERA5). GPT3 as last fallback. |
| Mapping function | M | **VMF3** (gridded, from TU Wien VMF data server) preferred; **VMF1** acceptable; GPT3/empirical VMF3 fallback. Stępniak et al. (2022) found VMF1 notably better than GMF for ZTD. |
| ZWD estimation | M | Random walk (see §B.1). |
| Horizontal gradients | R (on by default) | Chen & Herring gradient mapping; slow random walk or 6–12 h piecewise-linear. With 7° cutoff, unmodelled gradients alias into ZTD. |

### Estimation
| Item | Status | Notes |
|---|---|---|
| Sequential estimator (Kalman filter) | M | Numerically stable form (square-root/UD or Joseph form). |
| Backward smoother (RTS) | M (for post-processing) | Removes the convergence problem for archival data; equivalent in information to PRIDE's batch LSQ. |
| Ambiguity initialisation / reset | M | Large a priori variance; reset on slip/gap. |
| Formal uncertainties | M | From smoothed covariance; plus an empirical scaling factor determined in validation (formal σ is typically optimistic). |
| Convergence/quality flags | M | Per-epoch QC flags in output. |
| PPP-AR (WL/NL, LAMBDA, validation, partial fixing) | R → stage 6 | Implement after float is validated. |

### B.1 ZWD parameterisation — recommendation

| Option | Pros | Cons | Verdict |
|---|---|---|---|
| Piecewise constant, 5 min, no constraints | Simple | Too few observations per window → noisy, poorly determined | ✗ |
| Piecewise constant, 1 h (PRIDE PWC:60 style) | Stable | Misses convective variability; jumps at boundaries; 5-min output would be artificial | ✗ for your goal |
| **Random walk at observation rate + RTS smoother** | Physically reasonable (troposphere as a continuous stochastic process); natural in a Kalman filter; smoother gives each 5-min value information from before and after | Result depends on the process-noise choice | **✓ Recommended** |
| Piecewise linear, 5-min nodes, with relative constraints (batch) | Continuous; standard in Bernese-style processing; clean "value per 5-min node" | Equivalent to RW+smoother when constraints match; batch form is less convenient for later NRT | ✓ equivalent alternative |

**Choice:** ZWD as a random walk updated every observation epoch (e.g., 30 s), Kalman forward pass + RTS backward smoother, output the smoothed ZTD at exact 5-min epochs with its smoothed σ. Start with process noise ≈ **5–10 mm/√h**, and determine the value for Indian conditions (monsoon convection may need 10–15 mm/√h) through a sensitivity experiment against ERA5/radiosonde and PRIDE. Note clearly in the product: **adjacent 5-min values are strongly correlated**; the effective independent time resolution is ~15–30 min, and 5-min variations below a few mm are not resolved signal.

Output record: `timestamp (GPST and UTC) | ZTD | σ_ZTD | ZHD_apriori (+ source) | ZWD = ZTD − ZHD | G_N | G_E | n_sat | n_rejected | ambiguity status | QC flag`, plus a **SINEX_TRO 2.00** file for interoperability with PRIDE/IGS tools.

---

## C. What can safely be simplified or omitted in v1

| Can omit/simplify in v1 | Why it is safe for ZTD | When it matters |
|---|---|---|
| PPP-AR | For 24-h static GPS with final products, literature consistently reports only modest ZTD gains from AR (Stępniak et al. 2022 found little influence). | Short sessions, NRT, positioning. |
| Cycle-slip repair | New ambiguity after slip is fine for 24-h smoothed solution. | Short arcs, frequent slips. |
| 2nd/3rd-order ionosphere | Sub-mm to ~mm effect on ZTD. | High-TEC periods at low latitude — revisit in v2 given India's position near the equatorial anomaly crest. |
| Atmospheric/non-tidal/hydrological loading | mm-level coordinate effects. | Fixed-coordinate accuracy at the 1-cm level (Indian monsoon loading can have ~cm annual vertical amplitude — measure it first). |
| Ocean pole tide loading, J2 relativistic clock term | Sub-mm. | Never critical for ZTD. |
| Full eclipse yaw models | Exclude eclipsing satellites during turns instead. | If data loss is significant. |
| SNR weighting | Elevation weighting is standard and sufficient. | Multipath-heavy sites. |
| RINEX 2.11, multi-GNSS | Out of scope, but keep the architecture ready. | Later. |
| Real-time streams (SSR) | Out of scope. | NRT meteorology. |

**Not to be simplified:** antenna PCO/PCV (receiver + satellite), solid tide, ocean loading, consistent products, a good mapping function, residual screening, a correct ZHD source.

---

## D. Accuracy expectation

### D.1 ZTD

Distinguishing what the method can achieve, what I expect in practice, and what must be demonstrated:

| Condition | Theoretical/model capability (mature software, good site) | Expected practical — custom processor after validation | Status |
|---|---|---|---|
| **Convergence** (forward filter only, start of data or after long gap) | Fixed coordinates: ZTD < 1 cm after ~20–40 min; estimated coordinates: ~1–2 h. Initial errors 2–5 cm. | Similar or somewhat longer in early versions. **With the RTS smoother, convergence largely disappears for post-processing** except at data edges; use overlapping windows (e.g., 30 h) to protect day boundaries. | To be demonstrated |
| **After convergence, 5-min ZTD** | Formal σ ~2–4 mm; true error vs. independent truth ~5–8 mm | **5–10 mm RMS** | To be demonstrated |
| **Hourly-scale / 24-h static solution** | ~3–5 mm (IGS final troposphere product quotes ~4 mm) | **4–8 mm RMS**, mean bias < 2–3 mm vs IGS troposphere product at an IGS station | To be demonstrated |
| **Worst reasonable conditions** (monsoon convection, ionospheric scintillation, heavy multipath, frequent gaps, Rapid products) | ~8–15 mm | **10–20 mm RMS**, with occasional outliers > 30 mm that QC should flag | To be demonstrated |
| **Using Ultra-Rapid (observed half)** | ~1–2 cm | 1–2 cm+ | Not recommended for research product |

Justification: IGS final troposphere product accuracy (~4 mm) represents what mature PPP (GIPSY-based) achieves; published PPP vs. reference comparisons typically show STD of a few mm and biases of 1–3 mm; Stępniak et al. (2022) showed that model choices (especially mapping function) dominate over product and AR differences, and that sub-daily spurious signals at GPS-repeat harmonics remain even in established processing. These are literature-level magnitudes, not measurements of our code.

### D.2 PWV implied by ZTD error

Π = 10⁶ / [ρ_w R_v (k₃/T_m + k₂′)] ≈ **0.155–0.167** (T_m ≈ 270–295 K; Indian humid conditions ~0.16–0.165). So PWV ≈ ZWD / 6.2.

| ZTD error | PWV error (ZHD assumed perfect) |
|---|---|
| 1 mm | ~0.16 mm |
| 5 mm | ~0.8 mm |
| 10 mm | ~1.6 mm |
| 20 mm | ~3.2 mm |

### D.3 Other PWV error sources

| Source | Mechanism | Magnitude |
|---|---|---|
| **Surface pressure** | ZHD ≈ 2.2768 mm/hPa × P / f(φ,h) | 1 hPa → 2.3 mm ZHD → **~0.37 mm PWV**. Barometer at site (with height reduction ~0.12 hPa/m to antenna): < 0.2 mm. ERA5 interpolated + reduced: ~0.2–0.5 mm. GPT3 climatology: several hPa → **1–3 mm PWV, larger in cyclones/monsoon depressions** — not acceptable for research PWV. |
| **Mean temperature T_m** | Π ∝ ~1/T_m; 1 % Π per ~2.8 K | Bevis global model error 3–5 K in tropics → ~1–2 % → ~0.5–1 mm at 50 mm PWV. Use ERA5-derived T_m or GPT3 T_m, or a regional Indian T_m–T_s relation. |
| ZHD model itself | Saastamoinen is excellent given correct P | < 0.1 mm PWV |
| Antenna PCV/radome errors | Constant, elevation-dependent → ZTD bias | mm to ~1 cm ZTD → **0.2–1.5 mm PWV bias** |
| Mapping function | Low-elevation errors | ~1–5 mm ZTD (GMF vs VMF-class) → 0.2–0.8 mm PWV |
| Multipath | Site-specific, sidereal-repeating | 1–5 mm ZTD → 0.2–0.8 mm PWV, with sub-daily artefacts |
| Coordinate error (fixed mode) | ~0.2–0.4·δh into ZTD | δh = 1 cm → ~0.5 mm PWV; δh = 10 cm → **~5 mm PWV** (catastrophic, silent) |

**Realistic total PWV accuracy** vs radiosonde/radiometer: **1–2 mm RMS** in normal conditions, **2–3 mm** in monsoon/convective regimes or without on-site pressure. Ning et al. (2016) found ZTD uncertainty contributes >75 % of the GNSS IWV uncertainty budget, with the conversion-factor term growing with moisture — consistent with this budget. (Radiosonde errors themselves contribute to the observed RMS.)

---

## E. PRIDE comparison — realistic agreement

Two comparisons are needed, because they answer different questions:

1. **Harmonised** (same RINEX, coordinates, products [WUM], ANTEX, BLQ, cutoff, mapping function, ZHD a priori, ZTD stochastic model, gradients, sampling) → tests implementation correctness.
2. **Each processor at its best configuration** → tests scientific choices.

Expected ZTD_custom − ZTD_PRIDE (harmonised, after convergence/smoothing, 24-h sessions):

| Grade | Mean bias | STD | RMSE | P95 |abs| |
|---|---|---|---|---|
| **Excellent** | ≤ 1 mm | ≤ 2–3 mm | ≤ 3 mm | ≤ 6 mm |
| **Acceptable** | ≤ 2–3 mm | ≤ 4–5 mm | ≤ 5 mm | ≤ 10 mm |
| **Concerning** | > 3–5 mm, or station-/season-dependent | > 6–8 mm | > 8 mm | > 15 mm |

Also concerning even if statistics look good: a **diurnal or 12-h pattern** in the difference (tide/loading/wind-up convention error), **day-boundary jumps**, dependence of the difference on **PWV level or elevation cutoff** (mapping/antenna issue), a **correlation of the difference with satellite count/geometry**, and a consistently smaller or larger formal σ ratio.

Why not ~0 mm: batch-LSQ vs Kalman+smoother, different QC decisions (which arcs are rejected), different ambiguity handling and slightly different numerical implementations of the same models produce a few-mm STD even between correct programs. Harmonising the ZTD stochastic model is essential: PRIDE's default STO noise (~24 mm/√h) is looser than typical, and comparing a RW solution with a PWC:60 solution at 5-min resolution measures the parameterisation, not the code. Also compare hourly means, where parameterisation differences largely average out.

**Caution on circularity:** agreeing with PRIDE when both use the same products and models proves the implementation is consistent with PRIDE, not that the ZTD is correct. Common-mode errors (products, ANTEX, BLQ, mapping function, site multipath) are invisible in that comparison. Independent truth is required (§ Validation).

---

## F. Fixed-coordinate assessment

### Advantages
- Removes the parameter (height) most correlated with ZTD and receiver clock → better conditioned problem.
- Faster convergence and re-convergence after gaps/slips; more stable short-window and 5-min ZTD.
- Removes day-to-day height scatter as a ZTD noise source; improves consistency over long time series.
- Standard practice in operational GNSS meteorology networks (coordinates fixed or tightly constrained to a validated long-term solution).

### Risks — coordinate errors become ZTD biases
Rule of thumb (Santerre 1991 and subsequent geodetic literature; exact factor depends on cutoff and weighting): **ZTD bias ≈ 0.2–0.4 × height error**. Fixing a height that is too high gives a positive ZTD bias. Horizontal errors mostly alias into gradients/residuals and matter less. Realistic sources of height/position error for SOI coordinates:

| Error source | Possible size |
|---|---|
| **Epoch not propagated** (India plate ~5 cm/yr, mostly horizontal; vertical rates smaller but non-zero) | cm to dm horizontal after years |
| **Frame mismatch** (e.g., ITRF2008/2014 vs IGS20 orbits) | mm–cm |
| **Different ANTEX used to derive SOI coordinates** (e.g., igs08/igs14 vs igs20, or radome calibration differences) | mm to > 1 cm vertical |
| **ARP vs marker vs phase centre confusion** | Up to the antenna height (dm–m!) |
| **Tide system** (tide-free vs mean-tide) | cm-level, latitude dependent |
| Seasonal hydrological loading (monsoon), groundwater subsidence | mm to cm annual/secular vertical |
| Equipment changes, earthquakes | Arbitrary steps |

### Scientific reasons not to fix blindly
- Unmodelled real vertical motion (loading, subsidence) at cm level turns into a ZTD bias that looks like a humidity signal with seasonal structure — exactly the timescale of interest for monsoon PWV climatology.
- A fixed solution loses the diagnostic: in estimated mode, an antenna/metadata problem shows up as a coordinate offset; in fixed mode it hides in ZTD.

### Recommended approach: estimate, verify, then fix (or tightly constrain)
1. Process each station in **static mode** (one coordinate set per day) with the same processor, products and ANTEX for several weeks.
2. Form a robust mean (or a linear-trend fit for long series) in **IGS20 at the processing epoch**.
3. Compare with SOI coordinates propagated with SOI velocities into IGS20/ITRF2020 at the same epoch. Agreement within ~5 mm horizontal / ~1 cm vertical → consistent; otherwise investigate antenna/ARP/frame/tide conventions before trusting either.
4. Use the validated coordinate in **fixed mode** (your preference — scientifically sound once validated) or **tightly constrained** (σ ≈ 2–3 mm horizontal, 3–5 mm vertical). With such small σ the two are practically equivalent for ZTD; the constrained form additionally reports whether the data are pulling against the coordinate — a free quality check. I would implement both and default to constrained with a "fixed" option.
5. Refresh the coordinate periodically (e.g., monthly rolling solution) and after any equipment change.

Implementation must support: `fixed` · `constrained(σ_N, σ_E, σ_U)` · `static-estimated` · (later) `kinematic`, with the coordinate source recorded in every output file.

A **coordinate-perturbation experiment** (±1, ±3, ±10 cm vertical) should be part of validation to measure the actual ZTD sensitivity for your cutoff/weighting.

---

## G. IGS product strategy

### G.1 Which products are actually required

| Product | Float PPP | PPP-AR | Purpose / note |
|---|---|---|---|
| Orbit (SP3, 15 min or 5 min) | ✔ | ✔ | Day D plus D−1 and D+1 for interpolation at boundaries |
| Satellite clock (Clock RINEX, **30 s** preferred) | ✔ | ✔ (integer/phase clocks) | Must be same AC as orbit |
| ERP | ✔ | ✔ | Polar motion for pole tide; UT1 for Sun/Moon positions in ECEF |
| ANTEX (igs20.atx; igs14.atx for data before GPS week 2238) | ✔ | ✔ | Satellite + receiver PCO/PCV; must match products |
| Code OSB/DCB (C1C–C1W etc.) | ✔ if receiver lacks C1W/C2W | ✔ | Code consistency with clocks |
| Phase OSB (Bias-SINEX) | — | ✔ | From **same AC** as clocks |
| Attitude (ORBEX .OBX quaternions) | Optional | Recommended | Eclipse yaw; PRIDE uses them |
| VMF3 grids (ah, aw, ZHD, ZWD; 6-hourly) | ✔ (recommended) | ✔ | Mapping function + a priori ZHD; from TU Wien VMF data server; GPT3 grid file as fallback |
| BLQ (ocean loading coefficients) | ✔ | ✔ | Per station, once, from Onsala loading service (FES2014b) |
| Met data / ERA5 (pressure, T_m) | For PWV | For PWV | Separate PWV module |
| Broadcast navigation | Optional | Optional | Health flags, sanity checks only |
| GIM (IONEX) | — (v1) | — | Only if higher-order ionosphere added |

Not needed: IGS troposphere products (except as validation), station SINEX (except for IGS stations in validation), MGEX files for GPS-only.

### G.2 Naming (IGS long filenames, since GPS week 2238 / 27 Nov 2022)

| Tier | Orbit | Clock | ERP |
|---|---|---|---|
| Final | `IGS0OPSFIN_YYYYDDD0000_01D_15M_ORB.SP3.gz` | `IGS0OPSFIN_YYYYDDD0000_01D_30S_CLK.CLK.gz` (also 05M) | `IGS0OPSFIN_YYYYDDD0000_07D_01D_ERP.ERP.gz` (weekly; DDD = first day of GPS week) |
| Rapid | `IGS0OPSRAP_YYYYDDD0000_01D_15M_ORB.SP3.gz` | `IGS0OPSRAP_YYYYDDD0000_01D_05M_CLK.CLK.gz` | `IGS0OPSRAP_YYYYDDD0000_01D_01D_ERP.ERP.gz` |
| Ultra-rapid | `IGS0OPSULT_YYYYDDDHH00_02D_15M_ORB.SP3.gz` (HH = 00/06/12/18) | clocks only inside SP3 | `IGS0OPSULT_YYYYDDDHH00_02D_01D_ERP.ERP.gz` |

AC products follow the same pattern: `AAAVPPPTTT_YYYYDDDHHMM_LEN_SMP_CNT.FMT.gz` (e.g., `COD0OPSFIN_…_OSB.BIA.gz`, `WUM0MGXFIN_…_30S_CLK.CLK.gz`, `…_ATT.OBX.gz`). Pre-2238 data need the legacy short names (`igsWWWWD.sp3.Z`, `igsWWWWD.clk_30s.Z`, `igsWWWW7.erp.Z`). The resolver must support both and store the naming table as configuration, not code.

### G.3 A scientific correction to the proposed fallback logic

Your Final → Rapid → Ultra-Rapid logic is right for **float PPP with IGS combined products**, but there are two issues:

1. **IGS combined products cannot be used for PPP-AR.** As of 2026 the IGS PPP-AR Working Group combines AC clock/bias products experimentally, but there is **no official IGS combined integer-clock/OSB product**. AR therefore requires a single AC's consistent set (e.g., CODE `COD0OPSFIN/RAP` with OSB, Wuhan `WUM0MGXFIN/RAP` with OSB, or CNES/CLS `GRG`). Mixing IGS combined clocks with any AC's phase biases breaks AR.
2. **Rapid IGS clocks are 5-min; Ultra-rapid clocks are 15-min (in SP3), with a predicted half that is unusable** (ns-level clock errors). Interpolating these to 30 s adds error that goes straight into ZTD.

**Recommended product policy (configurable "product profiles"):**

- **Profile "AR-consistent" (default, and identical to PRIDE for benchmarking):** WUM (or COD) Final → WUM (or COD) Rapid → [Ultra only if an AC ultra product with consistent biases is verified to exist; otherwise float-only on IGS ULT]. Orbit, clock, OSB, ERP, attitude all from the same AC and tier.
- **Profile "IGS-combined float":** IGS0OPSFIN (30 s clocks) → IGS0OPSRAP (process at 5-min epochs) → IGS0OPSULT observed half only (process at 15-min epochs, flag "preliminary").
- Every output file records the tier, AC and filenames used. **Tier is a quality flag**: when a Final product becomes available, the day is automatically reprocessed and the Rapid/ULT result superseded.
- Never mix tiers or ACs within one processing window; if D−1/D/D+1 are not all available in one tier, drop to the lower tier for the whole window.

Latency (approximate): Final ~12–19 days; Rapid ~17–41 h; Ultra observed half ~3–9 h.

### G.4 Access mechanism

- **CDDIS** (`https://cddis.nasa.gov/archive/gnss/products/WWWW/`): HTTPS requires a NASA **Earthdata Login** (`.netrc` with `urs.earthdata.nasa.gov`, cookie handling, redirect following). Behind a proxy, both `cddis.nasa.gov` and `urs.earthdata.nasa.gov` must be allowed. CDDIS publishes checksum files per directory.
- **Mirrors without login** (useful under proxy restrictions): IGN (`igs.ign.fr`), BKG (`igs.bkg.bund.de`), Wuhan IGS data centre (`igs.gnsswhu.cn`, also the natural source for WUM products), CODE/AIUB (`ftp.aiub.unibe.ch/CODE`, source for COD OSB/DCB). Exact directory paths to be verified at implementation time and kept in configuration.
- Downloader design: ordered mirror list per product; HTTP(S)/FTP-TLS; configurable proxy (HTTP_PROXY/HTTPS_PROXY + explicit config, CA bundle); timeouts, retries with backoff; atomic writes (download to `.part`, verify, rename); content verification (decompress test, header parse, epoch coverage, satellite count, checksum where published); local cache keyed by canonical filename; an "offline" mode that only uses the cache; a manifest of what was used for each processed day.

### G.5 Date logic
Read first/last epoch from RINEX header and data (not only filename) → convert to GPS week/DOW and year/DOY → required product days = [first day − 1, last day + 1] for orbits; clocks for the processed days (+ neighbours if using an overlap window); weekly ERP file for each GPS week touched; VMF3 6-hourly grids spanning the window ± 6 h.

---

## H. Development difficulty by subsystem

Scale 1 (routine) – 5 (research-hard). "Accuracy impact" = how much this subsystem can damage ZTD if subtly wrong.

| Subsystem | Difficulty | Accuracy impact | Comment |
|---|---|---|---|
| RINEX 3 parser + observable selection | 2 | Medium | Many vendor variations; must be robust, not clever |
| Product downloader, cache, proxy | 2 | Low (but operationally critical) | Engineering; auth/proxy details are fiddly |
| SP3/CLK/ATX/ERP/BIA/OBX parsers, interpolation | 2 | Medium | Well-specified formats |
| Time systems & frames (GPST/UTC/TT/UT1, ECEF) | 2–3 | High if wrong | Leap seconds, week rollover, sub-ms tags |
| Deterministic corrections (tides, BLQ via IERS HARDISP, relativity, Sagnac, wind-up, PCO/PCV) | 3 | **High** | Documented but convention/sign errors are easy and silent |
| Satellite attitude & eclipse handling | 3–4 | Medium | Can be simplified by exclusion |
| Troposphere models (VMF3 grids, GPT3, gradients) | 2–3 | **High** | Grid interpolation, height corrections |
| Estimator (KF + RTS smoother, parameter bookkeeping, numerical stability) | 3–4 | High | Ambiguities appearing/disappearing; covariance handling |
| **Preprocessing / QC (slips, outliers, robust weights)** | **4** | **High** | Hardest practical part; Indian stations sit near the equatorial ionization anomaly crest, with steep TEC gradients and scintillation after sunset, especially around equinoxes and solar maximum |
| **PPP-AR (WL/NL, OSB conventions, validation, partial AR)** | **4–5** | Low–medium for 24-h ZTD; high if wrong fixes accepted | Conventions and integer-property consistency are unforgiving |
| Stochastic tuning & formal-error realism | 4 | **High** | Determined empirically |
| Validation framework | 3 | — | Where the credibility comes from |
| Multi-GNSS-ready architecture | 3 | — | Design discipline, not algorithms |

**Hardest overall:** robust preprocessing under low-latitude ionospheric conditions, PPP-AR bias handling, and getting all conventions exactly consistent (the "last 3 mm").

---

## I. Development roadmap

Each stage ends with a measurable gate.

| Stage | Content | Gate |
|---|---|---|
| **0. Reference setup** | Pick 3 pilot stations (1 IGS station in India, e.g. HYDE or IISC, + 2 SOI CORS), 30 days. Collect antenna/radome/ARP metadata, SOI coordinates with frame/epoch/velocity/tide convention, BLQ files. Run PRIDE in a harmonised configuration (F mode, WUM, VMF3/GMF, fixed cutoff, STO with chosen noise). | Reference dataset archived |
| **1. I/O foundation** | RINEX 3 parser; SP3/CLK/ATX/ERP/BIA parsers; time systems; product manager with cache, mirrors, proxy, verification, fallback profiles | Unit tests; automatic retrieval for 30 days works offline-after-cache |
| **2. Deterministic model engine** | Computes "modelled range" per satellite-epoch with every correction broken out | Each correction agrees with IERS reference software / independent implementation at < 0.1 mm; total modelled range vs PRIDE/RTKLIB intermediate output at mm level |
| **3. Float PPP-ZTD v0.1** | IF combination, fixed coordinates, KF + RTS smoother, clock WN, ZWD RW, gradients, ambiguities, basic MW/GF slip detection, residual screening, 5-min output + SINEX_TRO | vs PRIDE (harmonised): bias ≤ 3 mm, STD ≤ 5 mm on pilot set |
| **4. Coordinate modes** | Static-estimated and constrained modes; coordinate validation workflow vs SOI; perturbation experiments | Daily static coordinates repeat at few-mm; SOI–own offsets explained |
| **5. Robustness** | Robust weighting, SNR option, eclipse handling or ATT quaternions, overlap windows for day boundaries, receiver clock jumps, low-latitude ionosphere tuning, QC flags; ZWD noise sensitivity study | No day-boundary jumps; difference with PRIDE shows no diurnal structure |
| **6. PPP-AR** | WL (MW + OSB), NL (LAMBDA/bootstrapping), validation, partial AR, fixed re-solution | Fixing rate and fixed ambiguities agree with PRIDE; ZTD change float→fixed consistent with PRIDE's |
| **7. Validation campaign** | 10–15 stations × 1 year vs PRIDE, IGS troposphere product, ERA5, radiosonde | Report with all metrics; formal-error scaling factor |
| **8. PWV module** | Pressure (station barometer / ERA5 with height reduction), T_m (ERA5/GPT3/regional), Π, uncertainty propagation | PWV vs radiosonde/ERA5 within literature range |
| **9. Multi-GNSS (later)** | Galileo first (clean E1/E5a, good satellite metadata), then BDS-3 | Galileo-only and GPS+Gal ZTD consistent |

**Architecture for later multi-GNSS:** signals identified by (system, RINEX 3-char code); frequency-agnostic IF (or uncombined) builder; per-system receiver-clock/ISB parameters; biases always handled as OSB per signal; satellite metadata (block, PCO, frequencies) from ANTEX/SINEX metadata files; a named-parameter registry in the estimator so parameters can be added/removed without touching filter code.

### Language recommendation: **C++17/20 core + Python layer**

- **C++ core** (Eigen for linear algebra; no heavy dependencies): RINEX parsing, product parsing/interpolation, deterministic models, estimator/smoother, AR. Reasons: processing SOI's network (hundreds of stations × years = 10⁵–10⁶ station-days) and reprocessing campaigns; deterministic numerical behaviour; square-root/UD filter implementations; memory control for high-rate data.
- **Python layer** (via pybind11): configuration, product manager/downloader, batch orchestration, validation statistics, plotting, PWV module, ERA5 handling.
- Honest note: for a single station-day, pure Python/NumPy would be fast enough (the state vector is only ~15 parameters); the case for C++ is network-scale throughput and long-term maintainability of a numerical core, not single-run speed. RINEX parsing in pure Python is the practical bottleneck at scale.
- Reference code to consult (not copy): IERS Conventions software (DEHANTTIDEINEL, HARDISP), SOFA, TU Wien GPT3/VMF3 code, RTKLIB (BSD-2, useful for cross-checking intermediate values), PRIDE (GPLv3, algorithm reference only).

---

## J. Scientific risks — how a custom implementation could produce biased or unstable ZTD

1. **Antenna modelling errors**: wrong antenna/radome code, missing calibration, ARP offset mis-applied, IF-combined PCO/PCV computed wrongly → constant elevation-dependent bias (mm to > 1 cm ZTD).
2. **Wrong fixed coordinates**: frame, epoch, velocity, tide system, ANTEX-inconsistency, marker/ARP confusion → silent ZTD bias ≈ 0.2–0.4 × δh.
3. **Product inconsistency**: mixing ACs, tiers or ANTEX versions; interpolating sparse clocks.
4. **Convention/sign errors** in solid tide, ocean loading (phase-lag convention), pole tide, wind-up, relativity → spurious 12-h/24-h/sidereal signals in ZTD.
5. **Mapping function/a priori errors** amplified by low cutoff; empirical (GPT) mapping instead of NWP-based.
6. **Stochastic mis-tuning**: RW too tight → lagged, smoothed ZTD that underestimates convective peaks; too loose → noisy ZTD absorbing multipath.
7. **Undetected slips/outliers** → ZTD jumps; **over-aggressive editing** → many ambiguity resets → weak ZTD. Low-latitude ionospheric scintillation is a specific hazard for Indian stations.
8. **Wrong integer fixes** in PPP-AR → mm–cm ZTD offsets, worse than float.
9. **Day-boundary and edge effects** without overlap windows.
10. **Time-tag errors**: receiver ms clock jumps, GPST/UTC confusion, leap seconds in output timestamps.
11. **Eclipse-season yaw mismodelling** → wind-up/PCO errors on specific satellites.
12. **Receiver signal quirks**: L2C vs L2P tracking, quarter-cycle alignment, half-cycle flags, missing C1W.
13. **Unrealistic formal σ** (errors are time-correlated) → overconfident PWV uncertainties.
14. **Validation circularity**: agreement with PRIDE using shared products/models hides common-mode errors.
15. **Real station motion** (monsoon loading, groundwater subsidence) aliasing into ZTD when coordinates are fixed long-term.

---

## Validation strategy (detailed)

### Level 0 — component tests
Each correction tested against reference software or published test values (IERS test cases, SOFA); modelled-range decomposition compared term-by-term with an independent implementation. Target < 0.1 mm per term.

### Level 1 — PRIDE benchmark (harmonised, then best-config)
Same RINEX, coordinates, products, ANTEX, BLQ, cutoff, mapping function, ZHD a priori, ZTD stochastic model, gradients, sampling.

| Quantity | How to compare | Note |
|---|---|---|
| Receiver clock | Epoch-wise difference after removing a constant/linear term | Reference-dependent; look for jumps/noise, not level |
| Ambiguities | WL and NL integers after AR, per arc | Float IF values are only comparable with identical arc definitions |
| ZHD | Direct | Should be identical if same a priori; sanity check |
| ZWD, ZTD | 5-min and hourly-mean | Primary |
| ZTD formal σ | Ratio σ_custom/σ_PRIDE | Should be similar in shape and scale |
| Convergence | Forward-only run: time to |ΔZTD| < 1 cm and < 5 mm | Also from cold starts inserted mid-day |
| Data rejection, slips | Counts per satellite/arc/day | Large disagreements localise QC bugs |
| Continuity | Gaps, day-boundary jumps | |
| PWV | Using identical P and T_m | Checks propagation only |

Metrics: bias, STD, RMSE, MAE, Pearson r, P95 |Δ|, plus normalised differences Δ/σ, monthly bias/STD, diurnal composites of Δ.

Plots: overlaid ZTD time series with difference panel; scatter with 1:1 line and regression; histogram of Δ; monthly bias/STD bars; diurnal composite of Δ; Δ vs PWV level; Δ vs elevation cutoff; convergence curves; slip/rejection counts per day; residual stack vs elevation/azimuth (skyplot) for multipath/antenna diagnosis.

### Level 2 — independent truth
- **IGS final troposphere product** at IGS stations in India (e.g., HYDE, IISC; check availability) — process them exactly like SOI stations.
- **ERA5** ZTD/PWV at every station (full year, all conditions; known limits in convective cases).
- **Radiosondes** (IMD sites near SOI stations; account for radiosonde humidity biases and spatial separation).
- **Collocated/nearby stations** (< 10 km): ZTD difference after height correction should be small and stable.

### Level 3 — sensitivity and stress tests
Coordinate perturbation (±1/3/10 cm vertical, ±5 cm horizontal); fixed vs constrained vs static-estimated; cutoff 7/10/15°; mapping function VMF3/VMF1/GMF; ZWD noise 3/6/10/15 mm/√h; gradients on/off; injected slips and gaps; Final vs Rapid products; eclipse seasons; geomagnetic storm days; cyclone/heavy-rain days.

### How much data is convincing
- **Development pilot:** 3 stations × 30 days.
- **Publication-grade validation:** **≥ 10–15 stations × ≥ 1 full year** covering coastal humid, inland semi-arid, north-east monsoon, and high-altitude sites, including ≥ 2 IGS stations and ≥ 3 stations near radiosonde sites. One year is the minimum to cover pre-monsoon, monsoon, post-monsoon and winter regimes and the equinoctial ionospheric seasons. 3–5 years at a subset if climate-trend use is intended.

---

## K. Final recommendation

**Is it achievable?** Yes. A rigorous, independent GPS PPP processor optimised for ZTD, with fixed/constrained coordinates, 5-min smoothed output, automated product management and a PRIDE-benchmarked validation framework, is realistic. PPP-AR is achievable but is the highest-risk component and should come after float PPP is validated — it is not needed to reach research-grade 24-h ZTD.

**What I expect before any code is written (conservative):**

| Stage | Expected ZTD agreement with PRIDE (harmonised) | Expected ZTD vs independent truth |
|---|---|---|
| First working float version | 5–15 mm STD, possibly biased by several mm (bugs/conventions not yet found) | 1–2 cm |
| After component validation and QC hardening | bias ≤ 2–3 mm, STD 3–5 mm | 5–8 mm (hourly) |
| Mature, fully validated | bias ≤ 1 mm, STD ≤ 2–3 mm | 4–6 mm (hourly), comparable to PRIDE |
| PWV (with barometer or ERA5 pressure, NWP T_m) | — | 1–2 mm RMS normal; 2–3 mm monsoon/convective |

**Unrealistic claims to avoid:** sub-2-mm absolute ZTD accuracy at 5-min resolution; 5-min ZTD values as independent measurements; ultra-rapid-based "research-grade" PWV; PWV better than ~1 mm without an accurate pressure source.

**Decisions I would take:** float PPP first; KF + RTS smoother with ZWD random walk; VMF3 mapping and NWP-based ZHD; gradients on; 7° cutoff with elevation weighting (tested against 10°); coordinates estimated → validated against SOI → fixed or tightly constrained; one AC's consistent product set (same as PRIDE for benchmarking) with IGS combined products as a float-only alternative; C++ core with Python orchestration; clean-room implementation respecting PRIDE's GPL licence.

**Information needed from you before implementation:**
1. SOI coordinate metadata: frame, reference epoch, velocities, tide system, and whether they refer to marker, ARP or phase centre.
2. Antenna and radome models per station, and whether they are in igs20.atx.
3. Whether SOI stations have meteorological sensors (pressure, temperature).
4. Sampling rate of the RINEX (30 s / 15 s / 1 s) and receiver types.
5. Access to CDDIS (Earthdata account) or which mirrors your proxy allows.
6. The PRIDE version and configuration you currently use.

---

## Key references

- Bevis M. et al. (1992) GPS meteorology. *JGR* 97(D14); Bevis et al. (1994) *J. Appl. Meteorol.* 33.
- Zumberge J.F. et al. (1997) Precise point positioning. *JGR* 102(B3).
- Kouba J., Héroux P. (2001) PPP using IGS orbit and clock products. *GPS Solutions* 5(2); Kouba J. (2009) GPS yaw attitude during eclipse. *GPS Solutions* 13; Kouba J. (2015) *A Guide to Using IGS Products*.
- Wu J.T. et al. (1993) Phase wind-up. *Manuscripta Geodaetica* 18.
- Blewitt G. (1990) TurboEdit. *GRL* 17(3).
- Saastamoinen J. (1972); Davis J.L. et al. (1985) *Radio Science* 20.
- Chen G., Herring T.A. (1997) *JGR* 102(B9); Bar-Sever Y. et al. (1998) *JGR* 103(B3).
- Boehm J. et al. (2006) VMF1 / GMF. *JGR* 111; *GRL* 33. Landskron D., Böhm J. (2018) VMF3/GPT3. *J. Geod.* 92.
- Petit G., Luzum B. (2010) IERS Conventions 2010, IERS TN 36.
- Teunissen P.J.G. (1995) LAMBDA. *J. Geod.* 70.
- Ge M. et al. (2008) *J. Geod.* 82; Laurichesse D. et al. (2009) *Navigation* 56; Collins P. et al. (2010) *Navigation* 57; Geng J. et al. (2012, 2019) PPP-AR and PRIDE PPP-AR, *J. Geod.*/*GPS Solutions*; Schaer S. et al. (2021) CODE OSB/integer clocks, *J. Geod.* 95.
- Santerre R. (1991) Impact of GPS satellite sky distribution. *Manuscripta Geodaetica* 16.
- Ning T. et al. (2016) Uncertainty of GNSS IWV. *AMT* 9, 79–92.
- Stępniak K., Bock O., Bosser P., Wielgosz P. (2022) Outliers and uncertainties in GNSS ZTD from DD and PPP. *GPS Solutions* 26.
- IGS (2022) Guideline for the transition to IGS20 and long filenames v2.0; IGS PPP-AR Working Group pages.
- PRIDE PPP-AR II manual; PRIDE-PPPAR GitHub repository and wiki (v3.2.x).
