"""Stim circuits for circuit-level evaluation of any CSS code.

The schedule is deliberately generic rather than hand-optimised: all X checks
are measured, then all Z checks, with CNOTs packed into layers by greedy edge
colouring. IBM's gross-code paper uses a tuned depth-8 interleaved schedule,
so absolute numbers here will be somewhat worse than theirs. What matters for
the search is that every candidate is scored with the SAME schedule and noise
model, so comparisons between candidates are fair. Schedule optimisation for
a promising code is a separate, later step.

Noise model (uniform circuit-level depolarising, parameter p):
  * resets: flip with prob p
  * CX: DEPOLARIZE2(p)
  * measurements: flip with prob p
  * idle qubits during CX layers: DEPOLARIZE1(p_idle), default 0
"""
from __future__ import annotations

import numpy as np
import stim

from .bbcode import CSSCode


def _layers(pairs):
    """Greedy edge colouring: split (control, target) pairs into layers where
    no qubit is used twice."""
    layers = []
    used = []
    for a, b in pairs:
        for li, u in enumerate(used):
            if a not in u and b not in u:
                layers[li].append((a, b))
                u.update((a, b))
                break
        else:
            layers.append([(a, b)])
            used.append({a, b})
    return layers


def memory_z_circuit(code: CSSCode, rounds: int, p: float, p_idle: float = 0.0) -> stim.Circuit:
    hx, hz = code.hx, code.hz
    n, mx, mz = code.n, hx.shape[0], hz.shape[0]
    _, lz = code.logicals
    data = list(range(n))
    xanc = list(range(n, n + mx))
    zanc = list(range(n + mx, n + mx + mz))
    allq = data + xanc + zanc

    x_pairs = [(xanc[i], q) for i in range(mx) for q in np.nonzero(hx[i])[0].tolist()]
    z_pairs = [(q, zanc[i]) for i in range(mz) for q in np.nonzero(hz[i])[0].tolist()]
    x_layers = _layers(x_pairs)
    z_layers = _layers(z_pairs)

    c = stim.Circuit()
    total = 0  # measurements so far
    xm, zm = [], []

    def rec(abs_indices):
        return [stim.target_rec(a - total) for a in abs_indices]

    c.append("R", allq)
    c.append("X_ERROR", data, p)

    for r in range(rounds):
        # ---- X checks
        c.append("RX", xanc)
        c.append("Z_ERROR", xanc, p)
        for layer in x_layers:
            flat = [q for pair in layer for q in pair]
            c.append("CX", flat)
            c.append("DEPOLARIZE2", flat, p)
            if p_idle > 0:
                busy = set(flat)
                idle = [q for q in data + xanc if q not in busy]
                if idle:
                    c.append("DEPOLARIZE1", idle, p_idle)
        c.append("Z_ERROR", xanc, p)
        c.append("MX", xanc)
        xm.append(list(range(total, total + mx)))
        total += mx
        # ---- Z checks
        c.append("R", zanc)
        c.append("X_ERROR", zanc, p)
        for layer in z_layers:
            flat = [q for pair in layer for q in pair]
            c.append("CX", flat)
            c.append("DEPOLARIZE2", flat, p)
            if p_idle > 0:
                busy = set(flat)
                idle = [q for q in data + zanc if q not in busy]
                if idle:
                    c.append("DEPOLARIZE1", idle, p_idle)
        c.append("X_ERROR", zanc, p)
        c.append("M", zanc)
        zm.append(list(range(total, total + mz)))
        total += mz
        # ---- detectors
        if r == 0:
            for i in range(mz):
                c.append("DETECTOR", rec([zm[0][i]]))
        else:
            for i in range(mx):
                c.append("DETECTOR", rec([xm[r][i], xm[r - 1][i]]))
            for i in range(mz):
                c.append("DETECTOR", rec([zm[r][i], zm[r - 1][i]]))
        c.append("TICK")

    # ---- final data readout
    c.append("X_ERROR", data, p)
    c.append("M", data)
    dm = list(range(total, total + n))
    total += n
    for i in range(mz):
        supp = np.nonzero(hz[i])[0].tolist()
        c.append("DETECTOR", rec([dm[q] for q in supp] + [zm[-1][i]]))
    for j, row in enumerate(lz):
        supp = np.nonzero(row)[0].tolist()
        c.append("OBSERVABLE_INCLUDE", rec([dm[q] for q in supp]), j)
    return c


def dem_to_matrices(dem: stim.DetectorErrorModel):
    """Turn a detector error model into (H, L, priors) for BP-OSD.

    H[d, e] = 1 if error mechanism e flips detector d; L[o, e] = 1 if it flips
    observable o. Mechanisms with identical footprints are merged.
    """
    merged: dict = {}
    for inst in dem.flattened():
        if inst.type != "error":
            continue
        p = inst.args_copy()[0]
        dets, obs = [], []
        for t in inst.targets_copy():
            if t.is_relative_detector_id():
                dets.append(t.val)
            elif t.is_logical_observable_id():
                obs.append(t.val)
        key = (tuple(sorted(set(dets))), tuple(sorted(set(obs))))
        if not key[0]:
            continue  # undetectable; cannot be decoded either way
        q = merged.get(key, 0.0)
        merged[key] = q * (1 - p) + p * (1 - q)
    keys = list(merged)
    H = np.zeros((dem.num_detectors, len(keys)), dtype=np.uint8)
    L = np.zeros((dem.num_observables, len(keys)), dtype=np.uint8)
    pri = np.zeros(len(keys))
    for e, (d, o) in enumerate(keys):
        H[list(d), e] = 1
        if o:
            L[list(o), e] = 1
        pri[e] = merged[(d, o)]
    return H, L, pri
