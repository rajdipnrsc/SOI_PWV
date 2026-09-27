# Level-0 component tests (TDS § 30.1) — measured results

Generated 2026-09-27 by `validation/scripts/level0_report.py` at commit `7a6caa7`. All numbers are measured by this script; references are official IERS routine outputs (gfortran builds), RTKLIB 2.4.3 (compiled test-only), published test vectors, or analytic truth.

| Test | Dataset | Threshold | Measured | Result |
|---|---|---|---|---|
| SP3 Lagrange interpolation (11 nodes) | analytic GPS orbit, 5-min nodes, 1420 epochs | RMS < 1 mm, max < 5 mm | RMS 0.4578 mm, max 1.9618 mm | PASS |
| Solid Earth tide | official IERS DEHANTTIDEINEL, 10 stations x 100 epochs | < 0.1 mm | max 5.73e-10 mm | PASS |
| Ocean tide loading (HARDISP) | official IERS HARDISP, 6 cases x 48 epochs | < 0.1 mm | max 0.0005 mm (reference printed to 1 um) | PASS |
| Phase wind-up | RTKLIB 2.4.3 windupcorr, 7 arcs / 1060 epochs | < 0.001 cycle | max 5.0e-10 cycle | PASS |
| Sagnac / light time (rotation of r(t_tx) by w*tau) | independent inertial-frame light-time solution, 12 h pass | < 0.1 mm | max 3.73e-06 mm (RTKLIB first-order geodist differs by 0.095 mm: its truncation) | PASS |
| GMF mapping function | IERS GMF.F test vector | exact (1e-12) | 0.0e+00 | PASS |
| Moon position (Meeus 47) | Meeus example 47.a | < 1e-5 deg | dlon 3.1e-07 deg | PASS |
| Frame transformation (propagation-order equivalence) | HYDE SOI coordinate, 2005.0 -> 2024.5 | < 0.1 mm | 0.0e+00 mm | PASS |
| Kalman/RTS numerics (synthetic truth) | 24 h, 12 satellites, ZWD RW 6 mm/sqrt(h) | z-STD 0.9-1.1 | z-STD 1.032, RMSE 2.11 mm, NIS/dof 0.830 | PASS |

## Injected cycle slips (synthetic, § 30.4 method; real-data run still required)

Detected/injected per (dN1, dN2) over 3 datasets: (1,0): 17/18, (0,1): 17/18, (1,-1): 18/18, (2,1): 18/18, (2,2): 16/18, (5,5): 18/18, (10,9): 18/18, (3,0): 18/18, (1,1): 4/18.

False alarms on clean data: 8 in 20.2 satellite-days = 0.40 per satellite-day (gate ≤ 1).

(1,1)-type slips change neither the Melbourne–Wübbena combination nor the geometry-free combination by more than ~5 cm and are not reliably detectable by any MW/GF detector at 30 s; the § 30.4 gate refers to MW-detectable slips.

## Not yet covered at Level 0

* Full modelled-range decomposition against an independent implementation for one real station-day with CODE products (requires CODE products; the term decomposition is stored by `ppp/model.py`).
* Satellite PCO/PCV against RTKLIB on real igs20 entries (the convention is verified with a synthetic asymmetric antenna, `tests/test_antenna.py`).
* Pole tide: checked against the IERS formula re-implemented independently in the test (no official routine output available).
