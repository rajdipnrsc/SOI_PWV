# SOI_PWV — GPS PPP → ZTD → PWV (Python)

Implementation of `docs/TDS.md` (Technical Design Specification v1.1): a Python-only GPS precise-point-positioning
processor for Survey of India CORS that produces **5-minute zenith total delay (ZTD) and precipitable water vapour
(PWV)** per station-day using **CODE products only**, plus a national mapping layer.

> Status: all of System 1 and the separate PPP-AR and mapping layers are implemented and pass their Level-0
> component tests (`validation/level0_report.md`). The real-data gates (PRIDE benchmark, IGS/radiosonde comparisons)
> have **not** been run yet, so no accuracy may be claimed. See `validation/README.md`.

## Install

```bash
pip install -r requirements.txt          # numpy scipy pandas requests netCDF4 matplotlib hatanaka (+pytest)
```
Python ≥ 3.10. No compiled code, no other language.

## Run one station-day (System 1)

```bash
python pwv_ppp.py SITE.o SITE.n                       # e.g. python pwv_ppp.py Data/HYDE015.24o Data/HYDE015.24n
python pwv_ppp.py SITE.o SITE.n --xyz X Y Z           # SOI coordinate, ITRF2008 @ 2005.0 (metres)
```

| Option | Default | Meaning |
|---|---|---|
| `--xyz X Y Z` | from `stations.csv` / `station_coord.csv` | SOI coordinate (ITRF2008 @ 2005.0) |
| `--mode` | `auto` | `auto` (fix the SOI coordinate only if the same-day static check confirms it, § 9.12), `fixed`, `constrained`, `static` |
| `--out DIR` | `./results` | output folder |
| `--cutoff DEG` | 7 | elevation cutoff |
| `--ar` | off | PPP-AR layer (separate `*_FIXED_ZTD.csv`; experimental until the Stage-6 gate) |
| `--no-pwv` | off | skip PWV |
| `--offline` | off | product cache only |
| `--proxy URL` | `HTTPS_PROXY` | HTTP(S) proxy |
| `--quiet` / `--verbose` | normal | screen log detail |

Everything else is automatic: CODE orbits/clocks/ERP/OSB/attitude (`COD0OPSFIN` → `CODMOPSRAP` → `COD0OPSRAP`, never
mixed, never Ultra-Rapid), igs20 ANTEX (version named in the clock header), VMF3 grids, GPT3 grid and leap seconds are
downloaded, verified (decompression, header, coverage, satellites, SHA-256) and cached in `cache/` with `.meta.json`
sidecars. Unknown station metadata never stop a run: the defaults of TDS § 9.11 are applied, printed as WARNINGs and
recorded in the QC flags and the manifest. The run ends without ZTD only when PPP is impossible (no dual-frequency GPS
data, or no CODE orbit/clock: status `WAIT_FOR_TIER`, exit 2, or `FAIL`, exit 3).

**Network hosts** the program needs: `ftp.aiub.unibe.ch` (CODE, anonymous), `files.igs.org` (ANTEX),
`vmf.geo.tuwien.ac.at` (VMF3/GPT3), `hpiers.obspm.fr` or `data.iana.org` (leap seconds); optional mirrors
`zhw-b.s3.cloud.switch.ch`, `cddis.nasa.gov` (Earthdata login via `~/.netrc`). All URL templates are in
`ppp/settings.py`.

**Ocean loading:** the Onsala loading service delivers BLQ files by e-mail only. Request model FES2014b with
centre-of-mass correction once per station at http://holt.oso.chalmers.se/loading/ and save it as
`blq/<STATION>.blq`; without it the run continues with the `NO_OTL` flag.

**Station meteorology (optional, automatic):** a RINEX meteorological file of the same station and day placed next
to `SITE.o` (`HYDE015A.24m` or `*_MM.rnx`) is detected and used as priority-1 pressure (reduced to the antenna
height, QC'd against VMF3/GPT3, docs/decisions.md D-014). An ERA5 file in `cache/ERA5/` (written by `pwv_map.py`)
supplies pressure in barometer gaps and T_m. The source of each 5-min value is in the `P_source`/`Tm_source` columns.

### Outputs (in `--out`)

`<ST>_<YYYYDDD>_ZTD.csv`, `_ZTD.nc` (CF-1.10/ACDD-1.3), `_PWV.csv`, `.tro` (SINEX_TRO, optional format),
`_quicklook.png`, `_manifest.json` (provenance: every input file with SHA-256, family/tier, frame label, ANTEX,
coordinate object and transformation history, warnings, configuration hash, software version), `.log`, and
`provenance/configs/<hash>.json`. Delays, gradients and PWV are in mm; clock in ns (TDS § 17).
*Adjacent 5-minute values are strongly correlated (random walk + smoothing); they are not independent 5-minute
observations.*

## National maps (System 2, experimental)

```bash
python pwv_map.py results --day 2024-197 --era5 era5_pl_2024197.nc --dem dem_india.nc --grid G025 --cv
```
Options: `--hourly` (hourly archive grid), `--coarse` (0.5° variance-aware aggregate of the 0.25° grid),
`--as-of DATE` (product selection for reproducibility), `--support FILE` (supportable-spacing map from the § 31
experiment, `validation/scripts/spatial_experiment.py`). Station grey/blacklists (30-day residual history vs
neighbours) are maintained automatically in the maps folder.
ERA5/NWP background at station and cell height + GNSS residual regression kriging, σ and confidence classes, CF NetCDF
`INPWV_{GRID}_{STEP}_{YYYYDDD}_{TIER}_v1.0.nc`. Without `--era5`, ERA5 is downloaded automatically if `cdsapi` is
installed and a Copernicus key is in `~/.cdsapirc`; otherwise the benchmark height-normalised method is used and
labelled `HEIGHT_SCALING_ONLY`. Resolution/accuracy claims require the § 31 experiment.

## Network processing, reprocessing, daily operation (Stages 9-10)

```bash
network/run_network.sh /data/rinex/2024/197 32 --out results     # all station-days of a folder, 32 processes
python network/reprocess.py results --rinex-dir /data/rinex       # Rapid -> Final where CODE finals now exist
python network/reprocess.py results --rinex-dir /data/rinex --campaign 2025A --settings new.json   # campaign
python network/reprocess.py results --maps results/maps           # grids to rebuild (100 % FINAL inputs)
python network/monitor.py results --day 2024-197 --expect stations.csv
network/daily.sh /data/rinex results 32                           # cron recipe (see header of the script)
```
Each run is an ordinary `pwv_ppp.py` call; the product cache is safe for parallel processes. Outputs are
immutable: re-runs move earlier files to `results/superseded/<manifest_id>/`, and consumers (e.g. `pwv_map.py`)
select products by `best_available` (FINAL > RAPID_M > RAPID_0, then software version, then run time) or
`--as-of DATE` for reproducibility (TDS § 19.3). A campaign label (`PROCESSING_CAMPAIGN`) is stored in each manifest;
`PWV_PPP_SETTINGS=<file.json>` points a run to a settings override file.

## Layout

```
pwv_ppp.py          System 1 entry point          pwv_map.py   System 2 entry point
ppp/settings.py     all defaults, thresholds, URLs, constants (with references)
ppp/timesys.py rinex.py orbclk.py antenna.py products.py log.py         Stage 1
ppp/coords.py corrections.py troposphere.py model.py                     Stage 2
ppp/preprocess.py estimator.py output.py pwv.py ambiguity.py mapping.py  Stages 3-8
ppp/spatial_validation.py  § 24/§ 31 experiment tools (LOSO, block CV, thinning, E(d), supportable spacing)
ppp/data/           coefficient tables (VMF3 b/c, GMF, HARDISP 342 constituents, Meeus Moon)
tests/              pytest suite (+ reference fixtures from official IERS routines and RTKLIB)
network/            parallel run, reprocessing, daily cron and monitoring recipes (Stages 9-10)
validation/         gate status, Level-0 report, benchmark/sensitivity/tuning scripts
docs/TDS.md         the specification; docs/decisions.md, docs/verify_log.md
Data/               sample RINEX (HYDE, 2024 DOY 015/197/229)
```

## Tests

```bash
python -m pytest -q                        # ~70 tests, ~30 s
PWV_PPP_TEST_REF=<dir with igs20*.atx and gpt3_1.grd> python -m pytest -q --run-slow   # + end-to-end on HYDE
python validation/scripts/level0_report.py # regenerate the Level-0 report
```
