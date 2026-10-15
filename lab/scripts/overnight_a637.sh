#!/bin/bash
# a637 tile scale-up queue (WSL, from lab/): grafts at the ladder's d, then fresh ladders on every grafted code, then
# a 20M rung on the sizes whose d was still falling. Every stage is resumable: ris_ladder skips logged rungs and a
# graft whose output exists is skipped. Waits for a running ris_ladder.py first so the two don't share cores.
# Grafts are written as _g3 (hashed seeds, 20k + 2 x 200k filter, final 2 x 1M check with rollback). _graft (2k + 20k)
# and _g2 (small consecutive seeds) each let a distance-losing graft through, at n = 392 and n = 512.
#   bash scripts/overnight_a637.sh > results/hunt/a637_scale/overnight.out 2>&1 &
set -u
H=results/hunt/a637_scale
GF2=~/build/gf2_fast
while pgrep -f "[s]cripts/ris_ladder.py" > /dev/null; do sleep 60; done

graft() {  # bulk side, distance to keep
    [ -f $H/cand/a637_$1x$1_g3.npz ] ||
        python3 scripts/graft_reduce.py $H/cand/a637_$1x$1.npz --d $2 --gf2 $GF2 --out $H/cand/a637_$1x$1_g3.npz
}
ladder() {  # trial ladder, npz files
    python3 scripts/ris_ladder.py --gf2 $GF2 --log $H/ladder.jsonl --ladder "$1" "${@:2}"
}
FULL=20000,200000,1000000,5000000

# the staged [[360,8,17]] and [[404,8,18]] (_g2, grafted before hashed seeds) get a 20M rung first
ladder 20000000 $H/cand/a637_12x12_g2.npz $H/cand/a637_13x13_g2.npz
graft 14 20
ladder $FULL $H/cand/a637_14x14_g3.npz
# above n ~ 500 grafts mostly lose distance (512: rolled back from 500 to 511), so the plain codes' 20M rungs go first
ladder 20000000 $H/cand/a637_15x15.npz $H/cand/a637_16x16.npz $H/cand/a637_17x17.npz
graft 15 21
ladder $FULL $H/cand/a637_15x15_g3.npz
echo OVERNIGHT_DONE
