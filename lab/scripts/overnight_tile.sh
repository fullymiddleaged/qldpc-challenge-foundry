#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
# Single-layer tile hunt (scripts/hunt_tile.py, radius <= 4): 2 x 2 tiles at weights 6 and 8, then 3 x 3 tiles
# exhaustively at weights 6 and 8, then exact distances for the best of each. The 2 x 2 weight-4 tiles were all
# dominated on 30 Sep (results/hunt/tile_s_B2). Resumable: screened tiles are skipped via results/hunt/*/screened.txt.
cd "$(dirname "$0")/.." || exit 1
PY=.venv/Scripts/python.exe
T="$PY scripts/hunt_tile.py"
export MSYS_NO_PATHCONV=1
echo "$(date -Iseconds) start"
for w in 6 8; do $T screen --out results/hunt/tile_s_B2 --B 2 --w $w --hours 0.5; done
$T screen --out results/hunt/tile_s_B3w6 --B 3 --w 6 --hours 3.5
$T screen --out results/hunt/tile_s_B3w8 --B 3 --w 8 --hours 3.5
for run in tile_s_B2 tile_s_B3w6 tile_s_B3w8; do $T prove --out results/hunt/$run --max 12 --cube-secs 1200; done
echo "$(date -Iseconds) done"
grep -h '"new_entry": true' results/hunt/tile_s_*/proofs.jsonl 2>/dev/null | cut -c1-240
