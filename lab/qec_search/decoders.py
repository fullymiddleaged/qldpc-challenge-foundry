# SPDX-FileCopyrightText: 2026 Pete Salmond
# SPDX-License-Identifier: Apache-2.0
"""Decoders.

Preferred backend: the `ldpc` package (Roffe et al.), whose BP-OSD is the
standard baseline in qLDPC papers. If it is not installed, a pure-numpy
min-sum BP + OSD-0 fallback is used. The fallback is slower and weaker
(OSD-0 rather than OSD-CS order 7), so numbers from it are only for smoke
tests; `backend_name()` is recorded in every result so the two never mix.

All decoders expose: decode_batch(syndromes: (shots, m) uint8) -> (shots, n) uint8
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from . import gf2

try:  # ldpc >= 2
    from ldpc import BpOsdDecoder as _LdpcBpOsd  # type: ignore
    _LDPC = "ldpc2"
except Exception:  # pragma: no cover
    try:
        from ldpc import bposd_decoder as _LdpcBpOsdV1  # type: ignore
        _LDPC = "ldpc1"
    except Exception:
        _LDPC = None


def backend_name(force_numpy: bool = False) -> str:
    if _warned:
        return "numpy-bp-osd0 (ldpc failed)"
    return "numpy-bp-osd0" if (force_numpy or _LDPC is None) else _LDPC


class NumpyBpOsd:
    """Batched normalised min-sum BP with OSD-0 post-processing."""

    def __init__(self, h: np.ndarray, priors: np.ndarray, max_iter: int = 50, alpha: float = 0.75):
        self.h = gf2.as_gf2(h)
        self.m, self.n = self.h.shape
        self.priors = np.clip(np.broadcast_to(np.asarray(priors, float), (self.n,)), 1e-12, 0.5 - 1e-12)
        self.llr0 = np.log((1 - self.priors) / self.priors)
        self.max_iter = max_iter
        self.alpha = alpha
        rows, cols = np.nonzero(self.h)
        self.edge_check = rows
        self.edge_var = cols
        self.E = rows.size
        wmax = np.bincount(rows, minlength=self.m).max()
        # padded (m, wmax) table of edge ids; -1 = padding
        self.check_edges = -np.ones((self.m, wmax), dtype=np.int64)
        fill = np.zeros(self.m, dtype=np.int64)
        for e, c in enumerate(rows):
            self.check_edges[c, fill[c]] = e
            fill[c] += 1
        self.pad = self.check_edges < 0
        self.e2v = sp.csr_matrix((np.ones(self.E), (np.arange(self.E), cols)), shape=(self.E, self.n))
        self.h_sp = sp.csr_matrix(self.h.astype(np.int64))

    def _bp(self, synd: np.ndarray):
        S = synd.shape[0]
        sgn_s = 1.0 - 2.0 * synd[:, self.edge_check]  # (S, E)
        r = np.zeros((S, self.E))
        post = np.broadcast_to(self.llr0, (S, self.n)).copy()
        done = np.zeros(S, dtype=bool)
        hard = np.zeros((S, self.n), dtype=np.uint8)
        idx = np.where(self.pad, 0, self.check_edges)
        for _ in range(self.max_iter):
            vsum = self.llr0[None, :] + np.asarray(self.e2v.T.dot(r.T)).T  # (S, n)
            q = vsum[:, self.edge_var] - r  # var->check, (S, E)
            qt = q[:, idx]  # (S, m, w)
            mag = np.where(self.pad[None], np.inf, np.abs(qt))
            sg = np.where(self.pad[None], 1.0, np.where(qt < 0, -1.0, 1.0))
            tot = np.prod(sg, axis=2, keepdims=True)
            order = np.argsort(mag, axis=2)
            min1 = np.take_along_axis(mag, order[:, :, :1], axis=2)
            min2 = np.take_along_axis(mag, order[:, :, 1:2], axis=2)
            is_min = np.zeros_like(mag, dtype=bool)
            np.put_along_axis(is_min, order[:, :, :1], True, axis=2)
            other_min = np.where(is_min, min2, min1)
            rt = self.alpha * tot * sg * other_min  # sign excluding self: tot*sg (sg=+-1)
            new_r = np.zeros_like(r)
            valid = ~self.pad
            new_r[:, self.check_edges[valid]] = rt[:, valid]
            r = np.where(done[:, None], r, new_r * sgn_s)
            post = self.llr0[None, :] + np.asarray(self.e2v.T.dot(r.T)).T
            hard_new = (post < 0).astype(np.uint8)
            hard = np.where(done[:, None], hard, hard_new)
            ok = ~((self.h_sp.dot(hard.T).T.astype(np.int64) & 1) ^ synd).any(axis=1)
            done |= ok
            if done.all():
                break
        return hard, post, done

    def decode_batch(self, synd: np.ndarray) -> np.ndarray:
        synd = gf2.as_gf2(synd)
        hard, post, done = self._bp(synd)
        for s in np.nonzero(~done)[0]:
            order = np.argsort(post[s])  # least reliable (most likely flipped) first
            x = gf2.solve(self.h, synd[s], col_order=order)
            if x is not None:
                hard[s] = x
        return hard


class LdpcBpOsd:
    def __init__(self, h, priors, max_iter=100, osd_order=7):
        hd = gf2.as_gf2(h)
        h = sp.csr_matrix(hd)
        pri = np.broadcast_to(np.asarray(priors, float), (h.shape[1],)).copy()
        if _LDPC == "ldpc2":
            self.dec = _LdpcBpOsd(h, error_channel=list(pri), max_iter=max_iter, bp_method="minimum_sum",
                                  ms_scaling_factor=0.625, osd_method="osd_cs", osd_order=osd_order,
                                  schedule="parallel")
        else:
            self.dec = _LdpcBpOsdV1(hd, channel_probs=list(pri), max_iter=max_iter, bp_method="ms",
                                    ms_scaling_factor=0.625, osd_method="osd_cs", osd_order=osd_order)
        self.n = h.shape[1]

    def decode_batch(self, synd):
        out = np.zeros((synd.shape[0], self.n), dtype=np.uint8)
        for i, s in enumerate(gf2.as_gf2(synd)):
            out[i] = self.dec.decode(s)
        return out


_warned = False


def make_decoder(h, priors, force_numpy=False, **kw):
    global _warned
    if not (force_numpy or _LDPC is None):
        try:
            return LdpcBpOsd(h, priors, max_iter=kw.get("max_iter", 100), osd_order=kw.get("osd_order", 7))
        except Exception as e:  # API drift between ldpc versions
            if not _warned:
                print(f"WARNING: ldpc decoder failed to initialise ({e}); using numpy fallback. "
                      f"Send this message to Claude.")
                _warned = True
    return NumpyBpOsd(h, priors, max_iter=kw.get("max_iter", 50))
