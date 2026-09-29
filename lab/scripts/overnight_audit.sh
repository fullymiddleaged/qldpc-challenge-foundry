#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
# Second night (29-30 Sep 2026), from Git Bash in lab/:
#   1. exact distances for uncertified 2D-local board codes (scripts/audit_board_distances.py): any overclaim is a
#      distance revision, which the board accepts as a submission (TRACKS.md);
#   2. a weight-6 hunt round on l x 6 tori only (the only shape that laid out on 28-29 Sep, notes/hunt-2026-09-29.md),
#      new seed, skipping every genome screened before (results/hunt/*/screened.txt);
#   3. alongside, one core: the [[216,4,18]] MILP with a 14 h cap (4 h reached 13 <= d <= 18).
# Resumable: each stage skips what its output already holds.
cd "$(dirname "$0")/.." || exit 1
PY=.venv/Scripts/python.exe
H="$PY scripts/hunt_bilayer.py"
export MSYS_NO_PATHCONV=1
echo "$(date -Iseconds) start"
[ -f results/logs/milp_216_4_18_long.log ] || nohup $PY -m qec_search.exact_distance \
    --genome '{"l":18,"m":6,"A":[[0,2],[8,0],[13,0]],"B":[[0,1],[7,0],[11,0]]}' --claimed 18 --hours 14 \
    --out results/proofs/raw_216_4_18_long.json > results/logs/milp_216_4_18_long.log 2>&1 &
$PY scripts/audit_board_distances.py --local-only --retry-open --n-max 700 --d-max 40 --secs 900 --jobs 7 --hours 6
$H screen --out results/hunt/w6b --weight 6 --sizes 12x6,14x6,16x6,18x6,20x6,22x6,24x6,26x6,28x6,30x6 \
    --n-max 360 --count 200000 --hours 2 --seed 2 --workers 7
$H layout --out results/hunt/w6b --max 40 --workers 7
$H prove --out results/hunt/w6b --max 12 --cube-secs 1800 --jobs 7
echo "$(date -Iseconds) done"
grep -c '"verdict": "LOWER"' results/logs/board_distance_audit.jsonl
grep -h '"new_entry": true' results/hunt/w6b/proofs.jsonl 2>/dev/null | cut -c1-200
