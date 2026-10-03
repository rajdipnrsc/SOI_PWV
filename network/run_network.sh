#!/bin/sh
# Stage 9 recipe (TDS § 33.3, § 28): process every station-day RINEX file of a folder in parallel processes.
# Each station-day is an independent run of `python pwv_ppp.py SITE.o SITE.n`; the shared product cache is safe
# for concurrent runs (atomic writes). No scheduler is built into the program.
#
#   network/run_network.sh RINEX_DIR [JOBS] [extra pwv_ppp.py options ...]
#   e.g. network/run_network.sh /data/rinex/2024/197 32 --out /data/results
#
# Navigation file pairing: SITEdddX.yyo -> SITEdddX.yyn (or .yyN/.yyp); RINEX 3 long names
# XXXX00CCC_R_YYYYDDDHHMM_01D_30S_MO.crx.gz -> the same prefix with *_GN.rnx* / *_MN.rnx*.
set -u
DIR=${1:?usage: run_network.sh RINEX_DIR [JOBS] [options]}
JOBS=${2:-4}
shift; [ $# -gt 0 ] && shift
HERE=$(cd "$(dirname "$0")/.." && pwd)
PY=${PYTHON:-python3}
find "$DIR" -maxdepth 1 -type f \( -name '*.[0-9][0-9][oOdD]' -o -name '*.[0-9][0-9][oOdD].gz' -o -name '*.[0-9][0-9][oOdD].Z' \
     -o -name '*_MO.rnx*' -o -name '*_MO.crx*' \) | sort | while read -r OBS; do
  B=$(basename "$OBS")
  case "$B" in
    *_MO.*) P=$(echo "$B" | cut -c1-23); NAV=$(ls "$DIR"/"$P"*_GN.rnx* "$DIR"/"$P"*_MN.rnx* 2>/dev/null | head -1) ;;
    *)      S=$(echo "$B" | sed 's/\(\.[0-9][0-9]\)[oOdD]\(.*\)$/\1/');
            NAV=$(ls "$DIR"/"$S"[nNpP] "$DIR"/"$S"[nNpP].gz "$DIR"/"$S"[nNpP].Z 2>/dev/null | head -1) ;;
  esac
  if [ -z "$NAV" ]; then echo "SKIP $B: no navigation file" >&2; continue; fi
  printf '%s\n%s\n' "$OBS" "$NAV"
done | xargs -n 2 -P "$JOBS" sh -c '"$0" "'"$HERE"'/pwv_ppp.py" "$1" "$2" --quiet '"$*"' || echo "EXIT $? $1" >&2' "$PY"
