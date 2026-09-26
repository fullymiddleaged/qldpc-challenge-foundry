"""Scoring a code.

Two levels, mirroring the fast/slow loop used in recent code-discovery work:

* code_capacity_ler  (fast): i.i.d. X errors on data qubits, perfect
  syndromes, BP-OSD on H_Z. Cheap, good for screening thousands of codes.
* circuit_level_ler  (slow): full noisy syndrome-extraction circuit built in
  Stim, decoded with BP-OSD on the detector error model. This is the number
  that actually matters, so only run it on the survivors of the fast loop.

BB codes are symmetric under X <-> Z (H_Z is H_X with A, B transposed and
swapped), so the X-error rate is representative of both.
"""
from __future__ import annotations

import math
import time

import numpy as np

from .bbcode import CSSCode
from .decoders import make_decoder, backend_name


def wilson(fails: int, shots: int, z: float = 1.96):
    if shots == 0:
        return (0.0, 1.0)
    p = fails / shots
    d = 1 + z * z / shots
    c = p + z * z / (2 * shots)
    s = z * math.sqrt(p * (1 - p) / shots + z * z / (4 * shots * shots))
    return (max(0.0, (c - s) / d), min(1.0, (c + s) / d))


def code_capacity_ler(code: CSSCode, p: float, shots: int = 2000, max_fails: int = 200,
                      batch: int = 250, seed: int = 0, force_numpy: bool = False) -> dict:
    """Block logical error rate under i.i.d. X noise with probability p."""
    rng = np.random.default_rng(seed)
    _, lz = code.logicals
    hz = code.hz
    dec = make_decoder(hz, p, force_numpy=force_numpy)
    fails = done = 0
    t0 = time.time()
    while done < shots and fails < max_fails:
        b = min(batch, shots - done)
        err = (rng.random((b, code.n)) < p).astype(np.uint8)
        synd = (err.astype(np.int64) @ hz.T.astype(np.int64) & 1).astype(np.uint8)
        corr = dec.decode_batch(synd)
        resid = err ^ corr
        # residual must have zero syndrome (decoder always returns a valid
        # correction when OSD succeeds); failure = it flips a logical
        bad = ((resid.astype(np.int64) @ lz.T.astype(np.int64)) & 1).any(axis=1)
        fails += int(bad.sum())
        done += b
    lo, hi = wilson(fails, done)
    ler = fails / done if done else float("nan")
    return {
        "level": "code_capacity", "p": p, "shots": done, "fails": fails,
        "ler": ler, "ler_lo": lo, "ler_hi": hi,
        # per-logical-qubit rate, so codes with different k compare fairly
        "ler_per_lq": 1 - (1 - ler) ** (1 / max(code.k, 1)) if ler < 1 else 1.0,
        "decoder": backend_name(force_numpy), "seconds": round(time.time() - t0, 3),
    }


def circuit_level_ler(code: CSSCode, p: float, rounds: int | None = None, shots: int = 5000,
                      max_fails: int = 100, seed: int = 0) -> dict:
    """Memory-Z experiment at circuit level. Requires stim + ldpc."""
    from .circuits import memory_z_circuit, dem_to_matrices  # needs stim
    import stim  # noqa: F401

    d_guess = code.meta.get("d_ub", 6)
    rounds = rounds or max(2, int(d_guess))
    circ = memory_z_circuit(code, rounds=rounds, p=p)
    dem = circ.detector_error_model(decompose_errors=False)
    h, obs, priors = dem_to_matrices(dem)
    dec = make_decoder(h, priors)
    sampler = circ.compile_detector_sampler(seed=seed)
    fails = done = 0
    t0 = time.time()
    while done < shots and fails < max_fails:
        b = min(500, shots - done)
        det, actual = sampler.sample(b, separate_observables=True)
        corr = dec.decode_batch(det.astype(np.uint8))
        pred = (corr.astype(np.int64) @ obs.T.astype(np.int64)) & 1
        fails += int((pred != actual).any(axis=1).sum())
        done += b
    lo, hi = wilson(fails, done)
    ler = fails / done
    per_round = 1 - (1 - ler) ** (1 / rounds) if ler < 1 else 1.0
    return {
        "level": "circuit", "p": p, "rounds": rounds, "shots": done, "fails": fails,
        "ler": ler, "ler_lo": lo, "ler_hi": hi, "ler_per_round": per_round,
        "num_detectors": int(h.shape[0]), "num_error_mechanisms": int(h.shape[1]),
        "decoder": backend_name(), "seconds": round(time.time() - t0, 3),
    }
