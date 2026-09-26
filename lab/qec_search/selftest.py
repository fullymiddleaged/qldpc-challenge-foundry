"""Run this first on any new machine:  python -m qec_search.selftest

Checks, in order:
  1. the five IBM reference codes come out with the published [[n,k,d]]
  2. the decoder corrects errors (code-capacity sanity check)
  3. if stim is installed: the circuit builder produces deterministic
     detectors and observables, and decoding the circuit beats not decoding
"""
from __future__ import annotations

import sys
import time

import numpy as np

from .bbcode import build_bb
from .decoders import backend_name
from .evaluate import code_capacity_ler
from .known_codes import KNOWN


def main():
    ok = True
    print(f"decoder backend: {backend_name()}")
    print("1) reference codes")
    for name, g, exp in KNOWN:
        c = build_bb(g)
        d = c.distance_upper_bound(trials=400, stop_at=exp[2])
        got = (c.n, c.k, d)
        good = c.check_css() and got == exp
        ok &= good
        print(f"   {name:22s} got [[{got[0]},{got[1]},{got[2]}]]  {'OK' if good else 'FAIL expected ' + str(exp)}")

    print("2) decoder sanity ([[72,12,6]], code capacity)")
    c = build_bb(KNOWN[0][1])
    r = code_capacity_ler(c, 0.01, shots=3000, max_fails=1000)
    raw = 1 - (1 - 0.01) ** c.n  # chance of any error at all
    good = r["ler"] < 0.2 * raw
    ok &= good
    print(f"   p=0.01: LER={r['ler']:.2e} vs P(any error)={raw:.2e}  {'OK' if good else 'FAIL'}")

    print("3) circuit level")
    try:
        import stim  # noqa: F401
    except ImportError:
        print("   stim not installed - skipped (pip install stim)")
    else:
        from .circuits import memory_z_circuit, dem_to_matrices
        from .evaluate import circuit_level_ler

        code = build_bb(KNOWN[0][1])
        circ = memory_z_circuit(code, rounds=3, p=0.0)
        det, obs = circ.compile_detector_sampler().sample(200, separate_observables=True)
        good = not det.any() and not obs.any()
        ok &= good
        print(f"   noiseless circuit: {circ.num_detectors} detectors, {circ.num_observables} observables, "
              f"no spurious events: {'OK' if good else 'FAIL'}")
        noisy = memory_z_circuit(code, rounds=3, p=0.002)
        try:
            dem = noisy.detector_error_model(decompose_errors=False)
            H, Lm, pri = dem_to_matrices(dem)
            print(f"   detector error model: {H.shape[0]} detectors x {H.shape[1]} mechanisms  OK")
        except Exception as e:  # non-deterministic detectors raise here
            ok = False
            print(f"   detector error model FAILED: {e}")
        else:
            _, raw_obs = noisy.compile_detector_sampler(seed=3).sample(2000, separate_observables=True)
            raw_rate = raw_obs.any(axis=1).mean()
            t0 = time.time()
            code.meta["d_ub"] = 6
            res = circuit_level_ler(code, 0.002, rounds=3, shots=2000, max_fails=2000)
            good = res["ler"] < 0.5 * raw_rate
            ok &= good
            print(f"   p=0.002, 3 rounds: decoded LER={res['ler']:.2e} vs undecoded {raw_rate:.2e}  "
                  f"{'OK' if good else 'FAIL'}  ({time.time() - t0:.0f}s)")

    print("\nALL CHECKS PASSED" if ok else "\nSOME CHECKS FAILED - send this output to Claude")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
