# Addendum 1 — CODE-only Products & Nationwide ZTD/PWV Mapping from SOI CORS

*Prepared 22 Sep 2026 · Extends* `PPP_ZTD_Feasibility_Assessment.md` *· Supersedes its Section G product policy (WUM default) · No code written*

**Evidence labels used throughout**

- **[F] Established fact** — from authoritative documentation (CODE/AIUB, IGS, formats).
- **[L] Literature finding** — published result; transferability to India is not guaranteed.
- **[A] Engineering assumption** — my design judgement.
- **[E] Estimate to be validated** — a number that must be demonstrated experimentally.

---

## Part I — CODE-only product policy

### 1. What the "CODE-only" policy should cover

The policy applies to **Analysis-Centre-estimated products**: satellite orbits, satellite clocks, Earth rotation parameters, satellite code/phase biases and satellite attitude. These are jointly estimated and must come from **one CODE solution** to remain consistent.

A few required inputs are **not Analysis-Centre products** and cannot come from CODE by nature. Using them does **not** mix Analysis Centres, but each is recorded in provenance:

| Input | Source | Why it is not an AC-mixing issue |
|---|---|---|
| ANTEX (igs20.atx) | IGS convention file | CODE products are generated with the IGS ANTEX; the clock file header (`SYS / PCVS APPLIED`) names the ANTEX used, and the software must check it matches the local file **[F/A]** |
| Ocean-loading BLQ | Onsala loading service (FES2014b) | Station-side geophysical model, not a satellite product |
| Mapping function / a priori ZHD grids | TU Wien VMF3 (or GPT3) | Atmospheric model |
| Pressure, T_m (for PWV) | Barometer / ERA5 / NWP | Meteorological input |
| IERS conventions code (tides) | IERS 2010 | Standards |

### 2. CODE product inventory (verified against CODE's 2024 IGS Technical Report product table and CODE's `IAR_README.TXT`)

**Product family = AC (COD) + version character + campaign (OPS) + tier (FIN/RAP/ULT).** CODE runs *two* rapid solutions: **early rapid** (`COD0OPSRAP`) and **final rapid** (`CODMOPSRAP`, version character "M"). They are separate solutions and must not be mixed with each other. **[F]**

#### Final tier — `COD0OPSFIN` (GPS+GLONASS+Galileo; weekly submission, available after ~2 weeks) [F]

| Product | Filename | Sampling | Needed for |
|---|---|---|---|
| Orbit | `COD0OPSFIN_YYYYDDD0000_01D_05M_ORB.SP3.gz` | 5 min | Float + AR |
| Clock (RINEX clock 3.04) | `COD0OPSFIN_YYYYDDD0000_01D_30S_CLK.CLK.gz` | 30 s sat/ref clocks | Float + AR |
| Clock, high rate | `COD0OPSFIN_YYYYDDD0000_01D_05S_CLK.CLK.gz` | 5 s | Only for > 30-s data (not needed for 30-s RINEX) |
| ERP (daily) | `COD0OPSFIN_YYYYDDD0000_01D_01D_ERP.ERP.gz` | 1 d | Float + AR |
| ERP (weekly collection) | `COD0OPSFIN_YYYYDDD0000_07D_01D_ERP.ERP.gz` (Sunday only) | 1 d | Alternative to daily |
| Code + phase OSB (Bias-SINEX 1.00) | `COD0OPSFIN_YYYYDDD0000_01D_01D_OSB.BIA.gz` | daily | Float (code OSB for C1C etc.) + **mandatory for AR** |
| Attitude (ORBEX) | `COD0OPSFIN_YYYYDDD0000_01D_30S_ATT.OBX.gz` | 30 s | Recommended (eclipse yaw, wind-up, PCO) |
| Troposphere SINEX | `COD0OPSFIN_YYYYDDD0000_01D_01H_TRO.TRO.gz` | 1 h | **Validation only** (CODE's own ZTD at its network stations) |
| GIM (IONEX) | `COD0OPSFIN_YYYYDDD0000_01D_01H_GIM.INX.gz` | 1 h | Only if higher-order ionosphere is added |

#### Final-rapid tier — `CODMOPSRAP` [F]

| Product | Filename | Sampling |
|---|---|---|
| Orbit | `CODMOPSRAP_YYYYDDD0000_01D_05M_ORB.SP3` | 5 min |
| Clock | `CODMOPSRAP_YYYYDDD0000_01D_30S_CLK.CLK` | 30 s |
| ERP | `CODMOPSRAP_YYYYDDD0000_01D_01D_ERP.ERP` | 1 d |
| OSB | `CODMOPSRAP_YYYYDDD0000_01D_01D_OSB.BIA` | daily |
| Attitude | **Not listed in CODE's product table** | — |

#### Early-rapid tier — `COD0OPSRAP` (computed daily; AIUB states rapid series available after ~18 h) [F]

| Product | Filename | Sampling |
|---|---|---|
| Orbit | `COD0OPSRAP_YYYYDDD0000_01D_05M_ORB.SP3` | 5 min |
| Clock | `COD0OPSRAP_YYYYDDD0000_01D_30S_CLK.CLK` | 30 s |
| ERP | `COD0OPSRAP_YYYYDDD0000_01D_01D_ERP.ERP` | 1 d |
| OSB | `COD0OPSRAP_YYYYDDD0000_01D_01D_OSB.BIA` | daily |
| Attitude | `COD0OPSRAP_YYYYDDD0000_01D_30S_ATT.OBX` | 30 s |

#### Ultra-rapid tier — `COD0OPSULT` (4 updates/day, ~3 h latency) [F]

| Product | Filename |
|---|---|
| Orbit (latest, static name) | `COD0OPSULT.SP3` (5 min, includes prediction) |
| Orbit (24-UT solution, until early rapid is available) | `COD0OPSULT_YYYYDDD0000_01D_05M_ORB.SP3` |
| ERP | `COD0OPSULT.ERP`, `COD0OPSULT_YYYYDDD0000_01D_01D_ERP.ERP` |
| Clock RINEX | **Not listed** |
| OSB | **Not listed** |
| Attitude | **Not listed** |

#### Other CODE bias files [F]
Sliding 30-day and monthly DCB files (`P1C1.DCB`, `P1P2.DCB`, `CODE.BIA`, `P1C1yymm.DCB.Z`, …). **Not needed** when the daily OSB file of the same tier exists, because the daily OSB is the one consistent with that day's clocks. Use them only for pre-2018 float processing (see below).

**Historical data [F]:** CODE ambiguity-fixed ("integer-cycle") clocks and phase biases exist from **3 June 2018 (GPS week 2004)** for Rapid and Final (Schaer et al. 2021). Before that, CODE products support float PPP only. Files before the long-filename switch (27 Nov 2022) use legacy names (`CODwwwwd.EPH.Z`, `CODwwwwd_v3.CLK.Z`, `CODwwwwd.BIA.Z`, `CODwwwwd.OBX.Z`, `CODwwwwd.ERP.Z`); the resolver must hold both naming tables as configuration.

**Receiver code/phase biases [F/A]:** no Analysis Centre provides them, and none are needed: receiver code biases are absorbed by the receiver clock, receiver phase biases by the ambiguities (and cancel in between-satellite AR).

### 3. Answers to the seven checks

| Check | Answer |
|---|---|
| 1. Does CODE provide each product? | Orbit, clock, ERP, OSB (code+phase) — **yes at FIN, CODM-RAP and COD0-RAP**. Attitude — **FIN and COD0-RAP only** (not CODM-RAP). Ultra-rapid — **orbit + ERP only**. |
| 2. Naming | IGS long filenames as above; legacy short names before Nov 2022. Rapid files are listed without `.gz` in CODE's table — the resolver must try compression variants. |
| 3. Tier availability/latency | Final ≈ 2 weeks; rapid ≈ 18 h (early rapid; final-rapid later — exact timing must be measured by the downloader's availability log); ultra ≈ 3 h, 4×/day. |
| 4. Sampling | Orbit 5 min; clock 30 s (5 s Final); ERP 1 d; OSB daily; ATT 30 s. |
| 5. Family relationship | All products with the same `COD0OPSFIN` (or `CODMOPSRAP`, or `COD0OPSRAP`) prefix and day belong to one solution. OSB files are explicitly "related to" the orbit and clock of the same series. |
| 6. Internally consistent for PPP-AR? | **Yes, within one family.** CODE's IAR README: bias files are consistent with the associated ambiguity-fixed clock product; apply satellite OSBs to *all* observations; **"no phase bias, no AR"**; do not apply satellite PCO to the Melbourne–Wübbena combination. Schaer et al. (2021) report between-satellite clock consistency ~5 ps (~1.5 mm) and GPS clock/bias continuity across day boundaries in the Final series. |
| 7. If something is missing | See § 4 — explicit reporting, never silent substitution. |

### 4. Missing-product handling (no silent substitution)

| Situation | Can processing continue? | Scientifically acceptable fallback | Violates CODE-only? |
|---|---|---|---|
| **FIN not yet available** | Yes | Use `CODMOPSRAP` (complete set except ATT) → else `COD0OPSRAP` | No |
| **ATT missing** (CODM-RAP, or pre-ORBEX era) | Yes | (a) internal nominal-yaw + published eclipse yaw models, or (b) exclude GPS satellites during eclipse noon/midnight turns. Both are models, not AC products. | No. (Borrowing COD0-RAP attitude for a CODM-RAP solution is *same AC but different solution* — disallowed by default; allowed only with an explicit flag, recorded as `ATT_FROM_OTHER_SOLUTION`.) |
| **OSB missing** in a tier | Float: yes (code OSB only matters for C1C-type observables; use C1W/C2W where available, else flag `CODE_BIAS_MISSING`). AR: **no** | Float PPP only | No |
| **Clock 30 s missing but 5-min clocks available** | Not applicable for CODE (all tiers provide 30 s) | — | — |
| **Only ULT available** (< ~18 h after observation) | **No, for research-grade PPP.** CODE ULT has no clock RINEX and no OSB. | None within CODE. Wait for early rapid. | Any other clock source would violate the policy |
| **Data before June 2018** | Float only | CODE monthly DCB for code biases | No |
| **Mixed tiers across D−1/D/D+1** (e.g., D+1 only in rapid) | Yes | Drop the **whole window** to the lowest common tier, or process without the D+1 orbit (small edge-interpolation cost, flagged) | No |
| **ANTEX in clock header ≠ local ANTEX** | Stop | Fetch the matching igs20 version | No (convention file) |

The product manager returns one of: `OK_AR`, `OK_FLOAT_ONLY (reason)`, `WAIT_FOR_TIER`, `FAIL (missing list)`.

### 5. Revised downloader logic

```
RINEX → first/last epoch (from data records, not filename)
      → processing window [D_first − 1, D_last + 1]  (orbits), clocks/OSB/ATT/ERP for processed days (+ neighbours if overlap windows are used)
      → for tier in [COD0OPSFIN, CODMOPSRAP, COD0OPSRAP]:
            find every required file for the whole window in that one family
            (cache → AIUB HTTP → AIUB S3 mirror → CDDIS [Earthdata login])
            if complete: break
      → verify: decompress; parse headers; AC = COD; family/tier/version identical;
               epoch coverage; GPS satellites present; clock ANTEX == local ANTEX;
               OSB present for the selected signals (AR mode); checksum where published
      → decide mode: AR if OSB+integer clocks OK, else FLOAT (with reason)
      → write provenance manifest → process
      → schedule reprocessing when a higher tier appears (RAP → FIN)
```

ULT is never used for PPP under this policy. It could be used only to *display* preliminary orbit-based information, which has no value for ZTD.

**Provenance manifest per station-day [A]:** `analysis_centre=COD`, `family` (e.g., COD0OPSFIN), `tier`, `version_char` (0/M), exact filenames + SHA-256, sampling of each product, time span of each product, product header creation date, ANTEX name/version, processing mode (`FLOAT`/`AR`), AR fixing rate, fallbacks/missing flags, processor version and configuration hash.

**Access points [F/A]:** AIUB `http://ftp.aiub.unibe.ch/CODE/YYYY/` (anonymous; latest ultra/rapid files in `/CODE/`), an AIUB S3 mirror (`zhw-b.s3.cloud.switch.ch/aiub/CODE/`), and CDDIS (`/archive/gnss/products/WWWW/`, Earthdata login) for FIN/RAP. Which CODE series (especially CODM) are mirrored at CDDIS must be checked at implementation. Note for proxy planning: this cloud workspace could not reach `ftp.aiub.unibe.ch` directly (host not allow-listed), so allow-listing must be confirmed in your environment too.

**PRIDE benchmark impact [A]:** The harmonised PRIDE comparison must now run PRIDE with the **same CODE files** (PRIDE accepts user-specified orbit/clock/bias products; confirm with your installed version), not PRIDE's default WUM products.

**GPS-only in a multi-GNSS file [F/A]:** CODE OPS products contain GPS+GLONASS+Galileo; the processor simply reads GPS records. Future multi-GNSS inside the same family is natural for Galileo; BeiDou would need the separate `COD0MGXFIN` family (a family change, never mixed with OPS).

---

## Part II — Nationwide spatial ZTD/PWV mapping

### Key facts that drive everything

- **Network size [F]:** SOI reports **> 1,000 CORS** (launched Oct 2023; SOI states more than 1,105). India's area is ~3.29 million km².
- **Mean spacing if uniform [F, arithmetic]:** 1,000 stations → ~57 km; 1,105 → ~55 km; 1,500 → ~47 km. The real network is **not uniform**: expect ~20–40 km in dense states/plains and 100–200+ km gaps in the Himalaya, Thar, central forests and north-east **[A — must be mapped from the actual station list]**.
- **Water vapour varies strongly in space [L]:** Emardson et al. (2003, GPS network, S. California) found the RMS difference of ZWD between two sites follows σ ≈ c·L^0.5 + k·H with c ≈ 2.5 mm, L in km, k ≈ 4.8 mm/km (height difference), valid 10–800 km. The exponent was largely site-independent; c scales with regional humidity.
- **Height dependence is huge [F/L]:** ZHD falls ~0.27 mm per metre of height (pressure ~ −0.12 hPa/m). PWV falls roughly exponentially with a scale height of ~2 km. Between two stations 500 m apart, ZTD differs by ~135 mm from ZHD alone.

Illustrative raw ZWD variability from that structure function, converted to PWV (÷6.3). **The India row is my extrapolation [E]** (c raised to ~4.5 mm for monsoon humidity):

| Separation | ZWD RMS difference (c=2.5, California) | PWV equiv. | India monsoon (c≈4.5) PWV equiv. [E] |
|---|---|---|---|
| 10 km | 7.9 mm | 1.3 mm | 2.3 mm |
| 25 km | 12.5 mm | 2.0 mm | 3.6 mm |
| 50 km | 17.7 mm | 2.8 mm | 5.1 mm |
| 100 km | 25 mm | 4.0 mm | 7.1 mm |
| 200 km | 35 mm | 5.6 mm | 10 mm |

These are *station-to-station differences in the raw field*, not interpolation errors. Interpolating from several neighbours typically reduces the error to roughly 0.5–0.8 of this at the nearest-station distance, and removing an NWP background reduces it further **[A/E]**. Beyond ~100–200 km, GNSS adds little and the error approaches the background model's own error.

### 1. Is nationwide mapping feasible?

Three different products, with different answers:

| Product | Definition | Feasible? |
|---|---|---|
| **Station map** | Points coloured by station PWV | **Yes, fully.** Station PWV accuracy ~1–2 mm (from the main assessment). |
| **Interpolated surface** | Values between stations by any interpolator | Always *computable*; **meaningful only where stations are closer than the atmospheric correlation scale** and after height normalisation. |
| **Scientifically meaningful continuous field** | An estimate at every cell with a defensible, validated uncertainty | **Yes for PWV/ZWD in well- and moderately-covered regions, if built as NWP background + GNSS residual, with explicit height handling and an uncertainty/confidence layer.** Not for sparse/unsupported regions, where it can only be a background model with weak GNSS influence. |

Per variable:
- **ZTD**: directly interpolating absolute ZTD is **not** appropriate — it is dominated by the height-driven ZHD. A ZTD field is feasible only as ZHD(grid, from pressure + DEM) + ZWD(grid).
- **ZWD / PWV**: feasible via residual interpolation; this is the scientifically valuable field.

### 2. Recommended spatial resolution

Distinguish the **output grid spacing** (pixel size) from the **effective atmospheric resolution** (smallest real feature the data resolve). For an interpolated field, effective resolution is roughly **1.5–2× the local station spacing** for the GNSS-constrained part **[A]**; the NWP background contributes its own effective resolution (for ERA5, several grid lengths, i.e. ~100 km) **[L]**. The DEM can legitimately add finer detail only to the **height-dependent** part of PWV.

| Grid | Station spacing needed for the grid to carry GNSS information | Expected interpolation uncertainty (PWV) [E] | Monsoon gradients (e.g., Western Ghats windward/lee, ~50–100 km) | Coverage gaps | Defensible? |
|---|---|---|---|---|---|
| **0.5° (~55 km)** | ≤ ~75–100 km | ~2–4 mm (cells near stations); background-limited far away | Smeared | Few | **Yes nationally**, but coarse; wastes information in dense regions |
| **0.25° (~27 km)** | ≤ ~40–55 km | ~1.5–3.5 mm near stations; 3–5+ mm where spacing > 50 km | Partly resolved in dense regions | Moderate; explicit flags needed | **Yes — matches mean station spacing and ERA5 grid** |
| **0.125° / 0.1° (~11–14 km)** | ≤ ~20–30 km | ~1.5–3 mm only in dense clusters | Resolved where dense | Large fraction of India under-sampled | **Regional only**, or nationally as a *rendering grid* for DEM-driven height structure (must be stated as such) |
| **0.05° (~5.5 km)** | ≤ ~10 km | Not better than 0.1° except height term | — | Most of India unsupported | **No** for a national product; possible for local case studies |
| **0.025° (~2.8 km)** | ≤ ~5 km | Same | — | — | **No** — pixels would carry DEM information only |
| ~50 km | = 0.5° row | | | | |
| ~25 km | = 0.25° row | | | | |
| ~10 km | = 0.1° row | | | | |
| ~5 km | = 0.05° row | | | | |

**Recommendation [A, to be confirmed by the § 11 experiment]:**
1. **Primary national grid: 0.25°** (regular lat/lon, cell centres aligned with ERA5), effective resolution documented per cell (from local station spacing).
2. **Regional high-resolution option: 0.125°** (nests exactly in 0.25°; 0.1° is equivalent if nesting is not required), only in regions where station spacing ≤ ~25–30 km *and* cross-validation shows skill over the 0.25° product.
3. **Coarse/fallback: 0.5°** aggregates for sparse regions and for climatological products; cells beyond the supported distance are flagged or masked (§ 8).

For mountainous regions, report PWV at the **cell's mean DEM height** (stored as a layer), and offer a point-query tool that reconstructs PWV at any user height — this is more honest than a finer pixel grid.

### 3. Recommended temporal resolution

| Scale | Considerations | Assessment |
|---|---|---|
| **Station level: 5 min** | Random-walk ZTD with smoothing; effective independent resolution ~15–30 min; formal noise a few mm ZTD | Keep 5 min as the station product (as planned) |
| **Map: 5 min** | Adjacent maps differ by less than their own uncertainty except during convection; background NWP is hourly (ERA5) so the background part is time-interpolated; ~290 maps/day | **Not recommended as the primary national product** |
| **Map: 15 min** | Built from 15-min window means of station ZTD (averaging three 5-min values reduces noise); resolves the evolution of 50-km features (advection ~10 m/s ≈ 36 km/h → a 50-km feature changes over ~1 h); matches satellite imager cadence (INSAT-class 15–30 min) for joint studies | **Recommended primary map cadence** |
| **Map: 30 min** | Safe; loses part of convective onset timing | Good NRT/operational default if 15-min proves noisy in CV |
| **Map: 1 h** | Aligns with ERA5; ideal for climatology/validation archive | **Recommended archive/validation cadence** (derived from 15-min maps) |
| **Map: 3 h** | Only synoptic/climate use | Derived aggregate only |

### 4. Expected spatial and temporal accuracy [E — must be demonstrated]

| Coverage class | Typical conditions | PWV error (grid cell, 15-min/hourly) |
|---|---|---|
| Station itself | Any | 1–2 mm (2–3 mm monsoon/convective) |
| Well covered (nearest ≤ 25 km, small height difference) | Normal | ~1.5–2.5 mm |
| Moderately covered (25–50 km) | Normal | ~2–3.5 mm |
| Sparse (50–150 km) | Normal | ~3–5 mm, approaching the background error |
| Any class, monsoon convection | Active convection | × ~1.5–2 of the above |
| Mountain cells with large height difference to stations | Himalaya, NE hills, Ghats | Larger unless DEM-aware reconstruction proves effective; could exceed 5 mm |
| Unsupported (> ~150 km, or beyond empirical limit) | — | Background-model error only (NWP PWV error over India is typically a few mm and larger in convection **[L, magnitude to be measured]**) |

In relative terms, for PWV of 40–60 mm (monsoon), a 3-mm error is 5–8 %.

Temporal: maps at 15 min should have an additional ~0.5–1 mm error in convective conditions relative to hourly averages because the fast part of the signal is poorly constrained between stations **[E]**.

### 5. Interpolation methodology

Assessment specific to GNSS PWV over India:

| Method | Strengths for GNSS PWV | Weaknesses | Role |
|---|---|---|---|
| IDW | Simple, fast | No uncertainty, no height term, bull's-eyes around stations, ignores clustering | Quick-look only |
| Ordinary kriging | Uncertainty from covariance model | Stationary mean; must be applied to height-normalised residuals | Component of baseline |
| **Universal kriging / kriging with external drift (regression-kriging)** | Covariates (height, NWP PWV, latitude) in the trend; residual covariance gives uncertainty | Needs good variogram fitting per time/regime | **Baseline production method** (applied to GNSS − NWP residuals) |
| Gaussian process regression | Same maths as kriging, principled hyperparameter learning, easy anisotropy/non-stationarity | Scaling for space-time needs sparse/local approximations | Equivalent to baseline; research extension |
| RBF / thin-plate splines | Smooth fields; TPS with height covariate works for climatologies | No native uncertainty; overshoot in sparse areas | Not recommended for the operational field |
| Natural neighbour | No overshoot | No uncertainty, no extrapolation, no covariates | Not recommended |
| **Spatio-temporal Kalman filter / optimal interpolation (OI) against background** | Uses time continuity; OI with NWP background covariances is what meteorological analysis does | More complex; needs background-error statistics | **Advanced method** (natural upgrade of baseline) |
| Tomography (3-D wet refractivity) | Vertical structure | Needs dense networks (~10–30 km) and slant delays; ill-posed; not nationwide | Regional research in dense clusters only |
| Machine learning (RF, GBM, NN) | Can learn NWP bias patterns with terrain/season predictors | Uncalibrated uncertainty unless ensembles/quantiles; overfits clustered stations; must use spatial CV | Research: NWP bias-correction feeding the baseline, not a replacement |
| **NWP data assimilation of ZTD** | The scientifically best way to put ZTD into a 4-D weather field | Needs an NWP system (e.g., NCMRWF/IMD) | Recommend offering station ZTD to NWP centres; out of scope to build |

**Baseline pipeline (per 15-min slot) [A]:**

```
Station ZTD (CODE-based PPP, 15-min window mean, QC'd)
 → station ZHD from NWP/ERA5 surface pressure reduced to antenna height (or barometer)
 → station ZWD → station PWV (Π from T_m, ERA5/NWP)
 → background PWV_NWP at station coordinates AND station height (integrate NWP humidity profile from the antenna height)
 → residual r = PWV_GNSS − PWV_NWP
 → regression-kriging of r (covariance fitted per slot, constrained by seasonal/regional priors; local neighbourhood ~30 stations)
 → PWV_grid = PWV_NWP(at cell DEM height) + r̂ ;  σ_grid from kriging variance ⊕ station σ ⊕ background σ where applicable
 → ZWD_grid = PWV_grid / Π_grid ; ZHD_grid from NWP pressure at DEM height ; ZTD_grid = ZHD_grid + ZWD_grid
```

**Covariates/background use:**

| Item | Use |
|---|---|
| **ERA5** | Research/reprocessing background (hourly, 0.25°, ~5-day latency) and T_m/pressure |
| **NWP (IFS open data; NCMRWF NCUM for India)** | Operational background when ERA5 is not yet available; check whether the model already assimilates these GNSS stations (residuals would then not be independent) |
| **VMF3 grids** | Adequate for ZHD fallback; too coarse as a ZWD background |
| **GPT3** | Climatology — only for height-scaling defaults, never as the weather background |
| **Surface pressure** | Drives ZHD everywhere (grid and stations) |
| **DEM** (e.g., Copernicus GLO-30 aggregated) | Cell mean height; reconstruction of PWV/ZHD at cell height |
| **Latitude/longitude** | Weak trend terms; mostly captured by NWP background |
| **Station height** | Essential for height normalisation of each station value |

### 6. ZTD vs PWV mapping strategy

| Option | Pros | Cons (India) | Verdict |
|---|---|---|---|
| 1. Interpolate ZTD → convert grid to PWV | One field | ZTD dominated by ZHD's height dependence (~0.27 mm/m); Himalaya/Ghats height contrasts of km mean 100s of mm; grid conversion needs ZHD grid anyway | ✗ (unless ZTD height-normalised, which amounts to Option 2/3) |
| 2. Station PWV → interpolate PWV | Target variable; pressure/T_m handled per station | PWV scale height ~2 km → still needs height normalisation; no information in sparse regions | ✓ fallback when no NWP background available |
| **3. GNSS residual + NWP background (in PWV/ZWD space)** | NWP supplies terrain/vertical structure and large-scale field; GNSS corrects it near stations; uncertainty transitions naturally to background in sparse areas | Depends on NWP availability and biases; risk of double counting if assimilated | **✓ Recommended** |

**Topography:**
- **Himalaya**: km-scale height changes within tens of km; stations in valleys are unrepresentative of ridges. Interpolate only height-normalised residuals; reconstruct at DEM height; expect low confidence.
- **Western Ghats**: sharp windward/leeward moisture contrast during the SW monsoon across ~50 km; isotropic kriging smears across the crest. The NWP background captures part of this; anisotropic/non-stationary covariances are a research item.
- **Deccan plateau**: moderate relief (300–900 m), homogeneous — good for interpolation.
- **Indo-Gangetic plain**: flat, humid; likely densest network — best place for the 0.125° regional product.
- **Coastal regions**: land–sea contrast and sea breeze; no offshore stations → do not extend the field seaward beyond a short coastal buffer.
- **North-east hills**: complex terrain + extreme monsoon + likely sparse network → low confidence.

### 7./8. Sparse regions and uncertainty/confidence layers

**Coverage classes [A; thresholds set by the § 11 experiment]:**

| Class | Provisional criteria | Main PWV layer |
|---|---|---|
| 1 — Well covered | nearest station ≤ 25 km, ≥ 3 stations within 50 km, height difference to nearest ≤ 300 m | Filled |
| 2 — Moderate | nearest ≤ 50 km, ≥ 2 stations within 100 km | Filled |
| 3 — Sparse | nearest ≤ ~150 km (or empirical limit d*) | Filled, flagged "mostly background" |
| 4 — Unsupported | beyond d* | **Masked in the main layer**; background-only value available in a separate layer |

**Output layers:** PWV; ZWD; ZTD; ZHD; σ_PWV (total); σ_ZTD; distance to nearest station; number of stations within 50 km and 100 km; GNSS information fraction (1 − σ²_posterior/σ²_background); background PWV; confidence class; cell DEM height and height difference to nearest station; T_m and surface pressure used; land/sea mask; station-availability mask for that time slot.

### 9./10. Computational feasibility (1,000 stations)

| Item | Estimate |
|---|---|
| Station PPP solutions | 1,000 station-days/day; 288 ZTD epochs each → 288,000 ZTD values/day |
| PPP compute | C++ core ~1–10 s per station-day → ≤ 3 CPU-h/day (≈ 5–10 min on a 32-core server); Python-only ~30–60 s each → ~10–17 CPU-h/day |
| Grid cells (India bbox 6–38°N, 68–98°E) | 0.5°: 3,840 (~1,150 land) · 0.25°: 15,360 (~4,600 land) · 0.125°: 61,440 (~18,300 land) · 0.1°: 96,000 (~28,700 land) · 0.05°: 384,000 |
| Kriging per map | Local kriging, ~30 neighbours per cell: 0.25° ≈ 0.1–0.2 s; 0.125° ≈ 0.5–1 s per map. 96 maps/day → seconds to ~2 min per day |
| Storage (12 float32 layers, uncompressed) | 0.25°/15-min: ~70 MB/day (~26 GB/yr); 0.125°/15-min: ~280 MB/day (~100 GB/yr); 0.1°/5-min: ~1.3 GB/day (~480 GB/yr). Compression + land mask typically reduce 2–4× |
| RINEX input | ~1–2 TB/yr at 30 s (multi-GNSS Hatanaka-compressed; more at 1 s) |
| Hardware | One 32–64-core server, 128 GB RAM, ~10 TB storage covers years of operation |

**Near-real-time:** compute is not the constraint. **Product latency is**: under the CODE-only policy, the earliest PPP-capable CODE tier is early rapid (~18 h). So the operational product is **next-day** (rapid), superseded by a **final** reprocessing after ~2 weeks. True NRT (< 1–2 h) is **not achievable under a CODE-only policy** because CODE ULT has no clocks/biases — it would require real-time correction streams or another Analysis Centre, which the policy forbids. SOI data delivery latency and bulk-access terms also need confirming.

### 11. Cross-validation experiment to set the empirical resolution

**Design:**
1. **Data:** all SOI stations with ≥ 80 % completeness for one full year (plus IGS stations in India), CODE Final products, QC'd 5-min PWV.
2. **Station QC before mapping:** long-term GNSS − NWP difference per station to detect station biases (antenna/multipath), which would otherwise create spatial artefacts.
3. **Leave-one-station-out (LOSO)** for every station and time slot.
4. **Spatial block CV:** withhold 100–200 km blocks (prevents optimistic results from clustered stations).
5. **Thinning experiments:** thin the network to target spacings of 20, 30, 50, 75, 100, 150 km (random and spatially stratified), then predict the withheld stations.
6. **Methods compared:** background-only (ERA5/NWP), IDW, OK on height-normalised PWV, regression-kriging of residuals (baseline), GPR, space-time OI/KF.
7. **Grids compared:** 0.5°, 0.25°, 0.125°, 0.1° — evaluated at withheld stations after reconstruction at station height (also using the cell value at cell height to quantify representativeness error).
8. **Temporal:** maps at 5, 15, 30, 60 min, compared with the withheld stations' 5-min series; also test whether consecutive 5-min maps differ by more than their uncertainty.
9. **Regimes:** winter (JF), pre-monsoon (MAM), SW monsoon (JJAS), post-monsoon (OND), extreme rainfall days (IMD heavy-rain days, cyclones), dry conditions (PWV < 15 mm).
10. **Independent checks:** radiosondes, microwave radiometers where available, and spatial pattern checks against satellite PWV (e.g., clear-sky NIR products).

**Metrics:** RMSE, MAE, bias, correlation, P95 |error|, skill score vs background (1 − MSE_method/MSE_background); stratified by region, station density, nearest-station distance, elevation difference, PWV magnitude and regime; calibration of σ (fraction of errors within ±1σ, ±2σ).

**Defining the supportable resolution empirically:** fit error vs nearest-station distance per regime and region. Define **d*** as the distance where (a) GNSS skill over the background falls below ~0.2, or (b) RMSE exceeds a target (e.g., 3 mm PWV). A grid resolution is supported in a region where local station spacing ≤ d* and the finer grid improves CV skill over the coarser one by more than its sampling noise. Output a **"supportable resolution map"** of India per season.

**Convincing scale:** full network × ≥ 1 year, all four seasons, plus ≥ 10 extreme events.

### 12. National product specification (proposal)

| Element | Specification |
|---|---|
| Grid | Regular lat/lon (WGS84/EPSG:4326), **0.25°** national; **0.125°** regional nests; 0.5° aggregate |
| Domain | 6–38°N, 68–98°E; land mask + coastal buffer; islands flagged |
| Time | **15-min** maps (from 15-min station window means); **hourly** archive/validation product; 5-min station time series retained |
| Variables | PWV, ZWD, ZTD, ZHD, σ layers, confidence class, distance to nearest station, station counts, GNSS information fraction, background PWV, DEM height, T_m, pressure, masks |
| Tiers | `RAPID` (next-day, CODE rapid) → superseded by `FINAL` (CODE final, ~2 weeks) |
| **Primary archive format** | **NetCDF-4/HDF5, CF-1.x + ACDD metadata**, daily files, chunked (time × lat × lon), compressed; full provenance (CODE family/tier, processor version) |
| Analysis store | **Zarr** (same CF metadata), rechunked for long time-series access |
| GIS/visual delivery | **Cloud-Optimised GeoTIFF** for selected layers (hourly PWV, confidence) — derived product, not archive |
| Station product | NetCDF (CF discrete sampling geometry) + CSV/Parquet + SINEX_TRO |

### 13. Risks and limitations

1. Station biases (antenna, radome, rooftop multipath) turn into spatial bull's-eyes.
2. Heterogeneous and non-uniform network; clustered stations make CV look better than reality unless block CV is used.
3. Height handling errors in mountains dominate the error budget there.
4. NWP background biases in convection; possible double counting if the background assimilated the same GNSS stations.
5. Representativeness: a cell value ≠ a point value, especially in complex terrain.
6. PWV conversion at grid level depends on T_m and pressure grids.
7. Latency: CODE-only policy precludes true NRT.
8. Data access: bulk archival access to 1,000 SOI stations and terms of use must be confirmed.
9. Non-stationary covariance (monsoon vs winter, windward vs leeward) — a single variogram is insufficient.
10. Uncertainty calibration: kriging variance assumes the covariance model is right; must be checked with CV.

### 14. Final recommendation

- **CODE-only:** feasible and scientifically clean for **float PPP and PPP-AR** using one of three complete families: `COD0OPSFIN` → `CODMOPSRAP` → `COD0OPSRAP`. Attitude is missing for `CODMOPSRAP` (use an internal yaw model or eclipse exclusion; no policy violation). **CODE Ultra-Rapid cannot support PPP** (no clocks/biases), so the fastest CODE-only product is **next-day**. PPP-AR with CODE is possible from **June 2018** onward.
- **National mapping:** scientifically feasible for **PWV/ZWD** as an **NWP background + GNSS residual** field with explicit height handling and an honest confidence layer. Not by directly interpolating ZTD.
- **Resolution:** **0.25° national** (effective resolution ~50–100 km, documented per cell), **0.125° regional** where spacing ≤ ~25–30 km and CV proves skill, **0.5°** aggregates; **no national 0.05°/0.025° product** — at those spacings the pixels carry DEM information, not new atmospheric information.
- **Timing:** **5-min station ZTD + 15-min national PWV maps**, with an hourly archive. 5-min national maps are not supported by the network's spatial information content except possibly in dense clusters, and only if the experiment shows skill.
- **Accuracy [E]:** ~1.5–2.5 mm PWV in well-covered cells, ~2–3.5 mm moderate, ~3–5 mm sparse, larger in convection and mountains — all to be demonstrated by the LOSO/block-CV/thinning experiment, which also determines the final resolution map.

---

## References (Addendum)

- CODE Analysis Center, IGS Technical Report 2024 (product table), `ftp.aiub.unibe.ch/CODE/CODE_IGSTR.pdf`.
- CODE, `IAR_README.TXT` (integer ambiguity resolution with CODE products), `ftp.aiub.unibe.ch/CODE/IAR_README.TXT`.
- AIUB CODE Analysis Center web page (latencies, station counts).
- Schaer S. et al. (2021) The CODE ambiguity-fixed clock and phase bias analysis products. *J. Geod.* 95.
- IGS (2022) Guideline for the transition to IGS20 and long filenames v2.0.
- Emardson T.R., Simons M., Webb F.H. (2003) Neutral atmospheric delay in InSAR applications. *JGR* 108(B5).
- Xu W.B. et al. (2011) Interpolating atmospheric water vapor delay by incorporating terrain elevation information. *J. Geod.* 85.
- Survey of India / PIB (Oct 2023) National CORS network launch (> 1,000 stations).
