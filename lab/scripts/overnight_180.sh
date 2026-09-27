#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
# Waits for the depth benchmark (results/logs/certify_sym_bench.log, 5 arms) to finish, picks the orbital-branching
# depth for 180-20-14 from the 240-12-12 arms, then certifies upstream's two regression codes. Resumable: rerun it.
cd "$(dirname "$0")/.." || exit 1
PY=.venv/Scripts/python.exe
BENCH=results/logs/certify_sym_bench.log
until [ "$(grep -c 'wall_s=' "$BENCH" 2>/dev/null)" -ge 5 ]; do sleep 60; done
depth=$($PY - <<'PY'
import json
def cpu(t):
    s = json.load(open(f"results/logs/240-12-12.sym_bench{t}_d12.status.json"))["summary"]
    return s["cpu_hours"] if s["verdict"] == "PROVEN" else float("inf")
print(2 if cpu(2) < cpu(1) else 1)
PY
)
echo "$(date -Iseconds) bench done; 180-20-14 at depth $depth, 180-18-14 at depth 1"
$PY -m qec_search.certify_sym ../codes/180-20-14.json --depth "$depth" --tlim 600 --workers 4 &
$PY -m qec_search.certify_sym ../codes/180-18-14.json --depth 1 --tlim 600 --workers 4 &
wait
echo "$(date -Iseconds) done"
$PY -m qec_search.certify_sym ../codes/180-20-14.json ../codes/180-18-14.json --status
