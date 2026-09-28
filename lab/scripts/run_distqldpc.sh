#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
# Runs DistQLDPC (built in WSL at ~/src/DistQLDPC) on exported codes in parallel, one solver process each.
# Usage (from Windows): wsl.exe -d Ubuntu -- bash /mnt/c/.../lab/scripts/run_distqldpc.sh CPU_LIM PREFIX...
# A PREFIX is a path stem from export_distqldpc.py (results/distqldpc/<stem>_X) or a bundled name (BB_144_12_12),
# optionally followed by :FLAG for the solver, e.g. results/distqldpc/180-18-14_X_r0:-one-z=0.
# At most $JOBS solvers run at once (default nproc). Logs: lab/results/distqldpc/logs/<basename>.log, ending with a "wall N s" line.
set -u
BIN="$HOME/src/DistQLDPC/bin/distqldpc"
LAB="$(cd "$(dirname "$0")/.." && pwd)" || exit 1
LOGS="$LAB/results/distqldpc/logs"
mkdir -p "$LOGS" || exit 1
lim="$1"; shift
cd "$HOME/src/DistQLDPC" || exit 1              # bundled names resolve against data/matrices here
jobs_max="${JOBS:-$(nproc)}"
for arg in "$@"; do
    while [ "$(jobs -rp | wc -l)" -ge "$jobs_max" ]; do wait -n; done
    p="${arg%%:*}"; flag=""; [ "$p" != "$arg" ] && flag="${arg#*:}"
    case "$p" in /*|results/*) p="$(cd "$LAB" && realpath -m "$p")" ;; esac
    log="$LOGS/$(basename "$p").log"
    ( echo "c args: -cpu-lim=$lim $flag $p"; /usr/bin/time -f "wall %e s" nice -n 10 "$BIN" -cpu-lim="$lim" $flag "$p" ) > "$log" 2>&1 &
done
wait
for arg in "$@"; do p="${arg%%:*}"; echo "== $(basename "$p")"; tail -n 3 "$LOGS/$(basename "$p").log"; done
