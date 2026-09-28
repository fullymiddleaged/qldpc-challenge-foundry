#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
# Overnight hunt for new 2D-local bilayer entries (scripts/hunt_bilayer.py), from Git Bash in lab/:
#   covers  2- and 3-fold covers of our own proven bilayer codes (arXiv:2511.13560), weight 6
#   w8      random weight-8 BB codes on thin tori, n <= 336
#   w6      random weight-6 BB codes, n <= 336, where others submit too
# then lays out and proves the best survivors of each. Needs results/hunt/frontier_bilayer_w{6,8}.txt
# (scripts/dump_board_cells.py). Resumable: every stage skips what its output file already holds.
cd "$(dirname "$0")/.." || exit 1
PY=.venv/Scripts/python.exe
H="$PY scripts/hunt_bilayer.py"
export MSYS_NO_PATHCONV=1
for w in 6 8; do [ -f "results/hunt/frontier_bilayer_w$w.txt" ] || { echo "missing frontier w$w"; exit 1; }; done
echo "$(date -Iseconds) start"
$H screen --out results/hunt/covers --covers bb_144_8_12_12x6,bb_216_8_16_18x6,bb_288_12_16_24x6,bb_288_8_20_24x6 --h 2,3 --hours 1.5
$H screen --out results/hunt/w8 --weight 8 --sizes 8x6,10x6,12x6,14x6,16x6,18x6,20x6,22x6,24x6,26x6,28x6,8x8,10x8,12x8,14x8,16x8,18x8,20x8 --n-max 336 --count 150000 --hours 2.5
$H screen --out results/hunt/w6 --weight 6 --sizes 14x6,16x6,18x6,20x6,22x6,24x6,26x6,28x6,12x8,14x8,16x8,18x8,20x8 --n-max 336 --count 150000 --hours 1.5 --seed 1
$H layout --out results/hunt/covers --max 8
for run in w8 w6; do $H layout --out results/hunt/$run --max 40; done
for run in covers w8 w6; do $H prove --out results/hunt/$run --max 10 --cube-secs 2400; done
echo "$(date -Iseconds) done"
for run in covers w8 w6; do echo "== $run"; grep -h '"new_entry": true' results/hunt/$run/proofs.jsonl 2>/dev/null | cut -c1-200; done
