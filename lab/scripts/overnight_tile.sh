#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
# Single-layer tile hunt (scripts/hunt_tile.py, radius <= 4): every 2 x 2 and 3 x 3 tile of weight 4..8, and
# 100,000 random 4 x 4 tiles per weight, screened with a probe size first (about 0.13 s per tile), then exact
# distances for the best 15 of each box size. Resumable: screened tiles are skipped via results/hunt/*/screened.txt.
cd "$(dirname "$0")/.." || exit 1
PY=.venv/Scripts/python.exe
T="$PY scripts/hunt_tile.py"
export MSYS_NO_PATHCONV=1
echo "$(date -Iseconds) start"
for B in 2 3; do for w in 4 5 6 7 8; do $T screen --out results/hunt/tile_s_B$B --B $B --w $w --hours 1; done; done
for w in 4 5 6 7 8; do $T screen --out results/hunt/tile_s_B4 --B 4 --w $w --samples 100000 --seed $w --hours 1; done
for run in tile_s_B2 tile_s_B3 tile_s_B4; do $T prove --out results/hunt/$run --max 15 --cube-secs 1200; done
echo "$(date -Iseconds) done"
grep -h '"new_entry": true' results/hunt/tile_s_B*/proofs.jsonl 2>/dev/null | cut -c1-240
