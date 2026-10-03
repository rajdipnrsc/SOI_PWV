#!/bin/sh
# Stage 10 daily recipe (TDS § 35): process yesterday's RINEX, upgrade Rapid -> Final, rebuild maps, monitor.
# Example crontab (06:17 UTC; CODE rapid products are normally available ~12-18 h after the day ends, so the
# previous day is processed; no near-real-time promise, TDS § 36):
#   17 6 * * *  /opt/SOI_PWV/network/daily.sh /data/rinex /data/results 32 >> /data/results/daily.log 2>&1
#
# RINEX layout assumed: RINEX_ROOT/YYYY/DDD/ (adapt the DIR line below if different).
set -u
ROOT=${1:?usage: daily.sh RINEX_ROOT RESULTS_DIR [JOBS]}
RES=${2:?usage: daily.sh RINEX_ROOT RESULTS_DIR [JOBS]}
JOBS=${3:-8}
HERE=$(cd "$(dirname "$0")/.." && pwd)
PY=${PYTHON:-python3}
cd "$HERE" || exit 1
Y=$(date -u -d yesterday +%Y); D=$(date -u -d yesterday +%j)
DIR="$ROOT/$Y/$D"
echo "=== $(date -u +%FT%TZ) daily run for $Y-$D"
[ -d "$DIR" ] && sh network/run_network.sh "$DIR" "$JOBS" --out "$RES" || echo "no RINEX folder $DIR"
"$PY" network/reprocess.py "$RES" --rinex-dir "$ROOT" --jobs "$JOBS"
"$PY" network/reprocess.py "$RES" --maps "$RES/maps"
"$PY" network/monitor.py "$RES" --day "$Y-$D" ${EXPECT:+--expect "$EXPECT"} || echo "MONITOR: problems found"
