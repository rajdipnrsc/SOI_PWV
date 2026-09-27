# Sandbox end-to-end runs on the HYDE sample (diagnostic only — NOT products, NOT accuracy)

Date 2026-09-27. Purpose: prove that the complete System-1 chain (RINEX → preprocessing → model → EKF/RTS with
residual editing → auto coordinate mode → 5-min extraction → PWV → CSV/NetCDF/SINEX_TRO/quick-look/manifest)
runs robustly on real SOI RINEX 3 data. CODE products, VMF3 grids and the full GPT3 grid could not be downloaded in the
development sandbox, so these runs used:

* orbits/clocks: **broadcast ephemerides from the RINEX navigation file** (`tests/tools/broadcast_products.py`,
  test-only; ~1 m orbits, ~1-2 ns clocks, control-segment antenna offsets that are inconsistent with igs20 PCOs);
* troposphere: GPT3 meteorology (reduced grid without a_h/a_w) → **GMF emergency mapping** (`MF_FALLBACK`);
* no ocean-loading file (`NO_OTL`), no OSBs (`CODE_BIAS_MISSING`, receiver tracks C1C).

Because of the broadcast orbits/clocks the absolute ZTD/PWV values below are **not** a measure of the processor's
accuracy (a systematic ZTD offset of several cm, correlated with a ~0.1–0.3 m height offset in the same-day static
solution, is expected from control-segment vs igs20 satellite antenna offsets). They are recorded only to document
run-time behaviour.

| Day | Arcs | Slip/reset causes | Eclipse exclusions (sat-epochs) | IONO_ACTIVE epochs | NIS/dof | Phase outliers split | 5-min values | CONVERGED |
|---|---|---|---|---|---|---|---|---|
| 2024-015 (Jan) | 75 | LLI 58, MW 43, GF 6, GAP 116 | 1097 | 0 | 20.9 | 130 | 288/288 | 91 % |
| 2024-197 (Jul, monsoon) | 66 | LLI 59, MW 37, GF 9, GAP 64 | 1442 | 0 | 10.2 | 232 | 288/288 | 91 % |
| 2024-229 (Aug) | 91 | LLI 54, MW 74, GF 193, GAP 74 | 453 | 82 | 8.7 | 254 | 288/288 | 80 % |

Observations useful for the real-data stages:
* The HYDE Trimble Alloy tracks L2W intermittently (many LLI/GAP resets and short arcs); the stable per-satellite
  observable choice (docs/decisions.md D-003) removed 61 artificial L2W↔L2X switches per day.
* 2024-229 shows post-sunset equatorial-anomaly ionospheric activity (ROTI > 0.5 TECU/min, 82 epochs) with many GF
  slips/outliers: the § 30.4 storm-day tests should include such days.
* In mid-January and mid-July one GPS orbital plane is in eclipse season (β ≈ 0°): without CODE ORBEX attitude
  (COD0OPSFIN provides it) these satellites are excluded around shadow and noon turns (§ 7.11 default).
* Run time ≈ 40 s per 24-h, 30-s station-day on one core (auto mode = static pass with re-linearisation + check).
