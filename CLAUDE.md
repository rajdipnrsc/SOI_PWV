# CLAUDE.md — rules for Claude Code in this repository

Authoritative specification: `docs/TDS.md` (GPS PPP → ZTD → PWV Technical Design Specification v1.1).
Background only: `docs/background/*.md`.

## Six binding rules (TDS § 41)
1. **Python only.** No C++ or any other language, no compiled extensions (TDS § 33.1).
2. **Keep it simple.** Flat package as in § 33.2; plain functions + dataclasses; one `ppp/settings.py`; no frameworks,
   databases, class hierarchies or plugin systems. Clear, commented scientific code over clever abstractions.
3. **One command:** `python pwv_ppp.py SITE.o SITE.n` must work with no other input; optional flags only as in § 33.3.
4. **Download everything automatically** (CODE products, ANTEX, VMF3, GPT3, leap seconds), cache and verify it;
   never mix CODE families; never use other Analysis Centres.
5. **Never stop on unknown metadata.** Apply § 9.11 defaults, print a WARNING, flag outputs, continue. End a run only
   when PPP is impossible (no dual-frequency GPS data, no CODE orbit/clock) — cleanly, with message + status code.
6. **Show logs** on screen and in `results/<STATION>_<YYYYDDD>.log` (§ 33.4).

## Hard prohibitions (TDS § 36, abbreviated — read the full list)
No mixing CODE with other ACs; no mixing CODE families/tiers in a window; no CODE Ultra-Rapid for PPP; no silent
substitution or tier change; no hard-coded product paths/filenames/URLs/frame labels/observation codes (all in
`settings.py`); never use raw SOI ITRF2008@2005.0 coordinates as processing coordinates; never fix a coordinate the
same-day check contradicts unless `--mode fixed` (then warn); never apply ARP eccentricity twice or ignore the radome;
never interpolate raw ZTD/PWV across India without height treatment; never call grid spacing "atmospheric
resolution"; class-4 cells carry fill values; PRIDE/ERA5 are not truth; no accuracy claims without a measured gate
result; PPP-AR never overwrites float products; 5-min values are not independent observations; every output carries
manifest ID + config hash; never invent constants — resolve VERIFY items and record them in `docs/verify_log.md`;
no near-real-time promises; no C++/frameworks/databases; no required inputs beyond `SITE.o SITE.n`.

## Current stage
Stages 1–3 (+ auto coordinate mode of Stage 4, PWV module of Stage 7) are implemented and pass their Level-0
component tests (`pytest`). The Stage-3 **PRIDE gate and all Level-1/2 gates are NOT yet run** (they need CODE
products and PRIDE runs, which the development sandbox could not download). PPP-AR (`ppp/ambiguity.py`, `--ar`) and
national mapping (`ppp/mapping.py`, `pwv_map.py`) exist as separate, experimental layers that must not be used for
products before their gates (§ 30.5, § 31) are passed. See `validation/README.md` for gate status.

## Working rules
- Every module ships with small pytest tests (`tests/test_<module>.py`); run `python -m pytest -q` before committing.
- Every run writes a manifest; each gate produces a short report in `validation/`.
- Deviations from the TDS go into `docs/decisions.md`; VERIFY outcomes into `docs/verify_log.md`.
