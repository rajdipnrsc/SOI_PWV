# Validation and acceptance gates — status

The program is complete for System 1 (float PPP → 5-min ZTD/PWV), with separate experimental layers for PPP-AR
(Stage 6) and national mapping (Stage 8). **No accuracy number may be quoted** until the gate that measures it has
been run (TDS § 0.2, § 32, § 36.14). The table shows which gates have a measured result.

| Stage | Gate (TDS) | Status | Evidence |
|---|---|---|---|
| 0 | Pilot data, SOI coordinates, PRIDE runs | **Partly available** | 3 HYDE days (2024-015/197/229) + SOI coordinate supplied; PRIDE runs not yet available |
| 1 | Readers Level-0; product manager statuses; offline/proxy/retries/atomic writes/verification/manifests | **PASS (unit level)** | `tests/test_timesys.py`, `test_rinex.py`, `test_orbclk.py`, `test_antenna.py`, `test_products.py` (local HTTP server scenarios: fallback without mixing, ORBIT_EDGE, corrupt file → FAIL, WAIT_FOR_TIER, offline). Live CODE archive not reachable from the development sandbox → re-run the 30-pilot-day product gate in an environment with access. |
| 2 | Every § 7 correction at its Level-0 threshold; frame transformation; modelled-range decomposition | **PASS except decomposition vs independent implementation on CODE data** | `level0_report.md` (IERS DEHANTTIDEINEL/HARDISP official outputs, RTKLIB wind-up, analytic SP3, frame equivalence, GMF vector, Meeus) |
| 3 | Synthetic-truth filter/smoother; 5-min UTC extraction + flags; **PRIDE Stage-3 gate** | **Synthetic PASS; PRIDE gate NOT RUN** | `level0_report.md` (z-STD 0.999, NIS 0.92); PRIDE benchmark: `scripts/compare_ztd.py OURS.csv PRIDE.csv --pride` |
| 4 | Auto-mode check on pilot stations; perturbation experiment, β per station | **Logic tested; experiment NOT RUN** | `tests/test_coords.py::test_auto_decision_rules`; run `scripts/coordinate_sensitivity.py` with CODE products |
| 5 | Full QC, tuning (§ 6.7), Stage-5 PRIDE gate, Level-2 gates, Level-3 report | **NOT RUN** (needs CODE products, IGS/CODE troposphere, radiosondes) | tools: `scripts/tune_zwd_noise.py`, `scripts/compare_ztd.py` (IGS/CODE SINEX_TRO reader, harmonic amplitudes, gate verdicts); synthetic slip-injection results in `level0_report.md` |
| 6 | PPP-AR gates (§ 30.5) | **Synthetic PASS; real-data gate NOT RUN** | `tests/test_ambiguity.py` (all NL SD integers correct, ratio 357, P_s ≈ 1, fixed−float 0.2 mm) |
| 7 | PWV gates vs radiosonde (§ 30.6) | **NOT RUN** | conversion and uncertainty budget unit-tested (`tests/test_pwv.py`) |
| 8 | Spatial validation experiment (§ 31) | **NOT RUN** (needs ≥ 1 year of network products + ERA5) | kriging maths (UK variance vs dense formula), REML recovery, LOSO/block-CV/thinning/E(d)/supportable-spacing machinery tested (`tests/test_mapping.py`, `tests/test_spatial_validation.py`); run `scripts/spatial_experiment.py RESULTS --days 2024-001:2024-366 --write-support` → `stage8_spatial/` |
| 9–10 | Network reproducibility, operations | **Tools ready; NOT RUN** | `network/run_network.sh`, `network/reprocess.py` (Rapid→Final, campaigns, grid rebuild list), `network/daily.sh`, `network/monitor.py`; selection rules tested (`tests/test_products.py`, `tests/test_network.py`) |

## What the development sandbox could and could not do

* **Could:** run every Level-0 component test against official reference software outputs, synthetic-truth
  estimator tests, and the complete processing chain on the real HYDE RINEX files using **test-only broadcast
  orbits/clocks** (`tests/tools/broadcast_products.py`). Those runs prove that the chain executes robustly on real
  data (see `sandbox_runs.md`); because broadcast orbits/clocks are ~100× worse than CODE products, **their ZTD values
  are not a measure of accuracy** and are not products.
* **Could not:** download CODE products, IGS ANTEX, VMF3/GPT3 grids, IERS leap seconds or ITRF files (egress policy
  denied `ftp.aiub.unibe.ch`, `files.igs.org`, `vmf.geo.tuwien.ac.at`, `hpiers.obspm.fr`, `data.iana.org`,
  `itrf.ign.fr`, `cddis.nasa.gov`).

## First steps in an environment with network access

```bash
pip install -r requirements.txt
python -m pytest -q                                   # Level-0 suite
python pwv_ppp.py Data/HYDE015.24o Data/HYDE015.24n   # downloads CODE FIN products, igs20.atx, VMF3, GPT3
# HYDE is an IGS station: compare with the CODE (validation-only) or IGS final troposphere product
python validation/scripts/compare_ztd.py results/HYDE_2024015_ZTD.csv COD0OPSFIN_20240150000_01D_01H_TRO.TRO.gz --hourly
```
Then run the Stage-3 PRIDE benchmark (harmonised: same CODE files, ANTEX, BLQ, cutoff, VMF3, ZTD random walk
6 mm/√h, gradients) with `compare_ztd.py --pride`, the coordinate sensitivity experiment and the tuning experiment,
and record each outcome in a new report file in this folder.
