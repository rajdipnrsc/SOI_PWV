# Decisions and deviations from the TDS

One short entry per changed or added decision (date, decision, reason, evidence). Frozen decisions (TDS § 37) are
not changed by any entry below.

**D-001 (2026-09-27) — Implementation scope beyond the Stage-1 handoff.**
The project owner asked for a full Python implementation. All of System 1 (Stages 1–5, 7) and the separate layers
for PPP-AR (Stage 6) and national mapping (Stage 8) were implemented, each with Level-0 tests. The stage *gates*
that need real CODE products, PRIDE runs, IGS troposphere products, radiosondes or ERA5 (Stage-3 PRIDE gate onwards)
have **not** been run (the sandbox could not download CODE products); `validation/README.md` lists every gate and its
status. PPP-AR and mapping are marked experimental and must not be used for products before their gates pass.

**D-002 (2026-09-27) — Added module `ppp/model.py`.**
The observation model (geometry + all § 7 terms, term-by-term decomposition, SPP) is assembled in one small module
so that `corrections.py` remains a library of pure, separately testable functions (§ 33.2 simplicity rules allow a
split "when it clearly helps").

**D-003 (2026-09-27) — Observable selection is stable per satellite over the window.**
Per satellite and band, the highest-priority observable whose availability is ≥ 90 % of the best one is used for the
whole window; epochs where it is missing become gaps. Per-epoch priority selection made receivers with intermittent
L2W (semi-codeless) tracking flip between L2W and L2X, creating many artificial arcs (HYDE sample: 61 resets/day).
The TDS rule "one code and one phase observable per arc; a change starts a new arc" is kept.

**D-004 (2026-09-27) — ANTEX version fallback within the same family.**
If the exact ANTEX named in the CODE clock header (e.g. `igs20_2290`) cannot be obtained, the current file of the
same family (`igs20.atx`) is used with WARNING `ANTEX_VERSION_MISMATCH` instead of `FAIL` (TDS § 11.4), following the
Revision-1.1 never-stop principle; satellite calibrations of existing satellites do not change between igs20
versions. A different family (e.g. igs14 vs IGS20 products) still ends the run.

**D-005 (2026-09-27) — Mapping-function emergency fallback.**
Priority remains VMF3 grid → GPT3 (VMF3 with GPT3 a-coefficients). If neither VMF3 grids nor a GPT3 grid with a_h/a_w
are available, the IERS GMF is used with flag `MF_FALLBACK` (QC bit 17) and a WARNING, and the ZHD comes from GPT3
(or a standard atmosphere, `ZHD_CLIMATOLOGY`). The TDS allows GMF only for PRIDE benchmarking; it is used here only so
that a run never stops for a missing climatology file.

**D-006 (2026-09-27) — Reserved QC bits used.**
Bits 15 `NO_OTL`, 16 `NO_ANTENNA_CALIBRATION`, 17 `MF_FALLBACK` (TDS § 16.1 bits 15–31 reserved) so that § 9.11
assumptions are visible per row, as § 0.5 requires ("record it in the output flags").

**D-007 (2026-09-27) — Geoid for orthometric heights.**
The GPT3 1° grid undulation is used for H_ant (needed only for pressure reduction; 1 m error ≈ 0.1 hPa ≈ 0.04 mm PWV).
EGM2008 at full resolution can replace it without interface changes.

**D-008 (2026-09-27) — Slip-detection details [A].**
(a) The MW running standard deviation is floored by an elevation-dependent code-noise model
0.25 cycle / sin(e); (b) the GF threshold scaling g(Δt, e) = max(1, Δt/30 s) · min(√(1/sin e), 2) (×3 when
`IONO_ACTIVE`); (c) GF exceedances are confirmed by the next epoch (the residual jump must persist) to separate slips
from outliers. Synthetic tests: 100 % detection of slips changing MW by ≥ 2 cycles or GF by ≥ 0.3 m, ≥ 80 % for
single-cycle L1- or L2-only slips, ≤ 1 false alarm per satellite-day on clean data. Real-data tuning belongs to the
§ 30.4 injected-slip experiment.

**D-009 (2026-09-27) — Test-only broadcast product harness.**
`tests/tools/broadcast_products.py` builds SP3/CLK-like objects from the RINEX navigation file so the full chain can
be exercised where CODE products cannot be downloaded. It is not imported by the program; runs using it are
labelled `PRODUCTS_OVERRIDE`/`TEST_BROADCAST` and are never products.

**D-010 (2026-09-27) — National grid alignment.**
The 0.25° national grid uses cell centres at multiples of 0.25° (ERA5-coincident, TDS § 27.1 option) so that the
ERA5 background is sampled exactly at cell centres.
