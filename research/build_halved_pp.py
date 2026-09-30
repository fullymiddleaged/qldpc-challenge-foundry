"""Build symplectic-halved CPM pair-partition codes (arXiv:2609.30069 recipe).

Given (J, L, P), a fixed-point-free involution sigma of the L block columns and
eta = +/-1, Prop 4 says a CSS parent with exponent arrays E, D

    D[j][l] = eta * E[j][sigma(l)]   (rho = id, alpha = beta = 0 gauge)

admits symplectic halving under pi(l, t) = (sigma(l), eta * t) whenever each
off-diagonal difference array T[i][j](l) = E[i][l] - eta * E[j][sigma(l)] takes
each value an even number of times.  We impose that as pair-partition equations
(choose a perfect matching M_ij of the L columns per cell i<j, require
T[i][j] constant on its pairs); diagonal cells are automatic for eta = -1.

Outputs (under research/scratch_out/):
  halved_<tag>_parent.npz   hx, hz   the CSS parent (2-fold of the stabilizer)
  halved_<tag>_fold.npz     s        the folded stabilizer S = (A | B)

Both are then checkable with ./qldpc submit <file> --dry-run.
"""

import argparse
import itertools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for p in (os.path.join(HERE, "kit"), HERE, os.path.join(os.path.dirname(HERE), "verify")):
    if p not in sys.path:
        sys.path.insert(0, p)

import heuristic_distance as hd  # noqa: E402
from css import compute_k, rank, verify_css  # noqa: E402

OUT_DIR = os.path.join(HERE, "scratch_out")


# ---------------------------------------------------------------- F_p linear algebra
def nullspace_mod(rows, P):
    """Basis of {x : rows @ x = 0 mod P} over the field F_P (P prime)."""
    A = [list(int(v) % P for v in r) for r in np.atleast_2d(np.asarray(rows))]
    m, n = len(A), len(A[0])
    piv_cols, r = [], 0
    for c in range(n):
        piv = next((i for i in range(r, m) if A[i][c] % P), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], P - 2, P)
        A[r] = [(v * inv) % P for v in A[r]]
        for i in range(m):
            if i != r and A[i][c] % P:
                f = A[i][c]
                A[i] = [(a - f * b) % P for a, b in zip(A[i], A[r])]
        piv_cols.append(c)
        r += 1
        if r == m:
            break
    free = [c for c in range(n) if c not in piv_cols]
    basis = []
    for fc in free:
        x = [0] * n
        x[fc] = 1
        for i, pc in enumerate(piv_cols):
            x[pc] = (-A[i][fc]) % P
        basis.append(np.array(x, dtype=np.int64))
    return basis, len(free)


def equal_col_pairs(H):
    """Pairs of identical columns (each gives a weight-2 kernel vector)."""
    packed = np.packbits(np.asarray(H, dtype=np.uint8).T, axis=1)
    seen, out = {}, []
    for idx, row in enumerate(packed):
        key = row.tobytes()
        if key in seen:
            out.append((seen[key], idx))
        else:
            seen[key] = idx
    return out


def cycles4(X, P):
    """Row-pair differences X[i][j]-X[i'][j] must be distinct over j."""
    J, L = X.shape
    for i in range(J):
        for ip in range(i + 1, J):
            diffs = [(int(X[i, j]) - int(X[ip, j])) % P for j in range(L)]
            if len(set(diffs)) != L:
                return True
    return False


def cycles6(X, P):
    """Alternating exponent sums over row triples must never vanish."""
    J, L = X.shape
    for a, b, c in itertools.permutations(range(J)):
        for j0 in range(L):
            for j1 in range(L):
                if j1 == j0:
                    continue
                for j2 in range(L):
                    if j2 in (j0, j1):
                        continue
                    s = (
                        (int(X[a, j0]) - int(X[b, j0]))
                        + (int(X[b, j1]) - int(X[c, j1]))
                        + (int(X[c, j2]) - int(X[a, j2]))
                    ) % P
                    if s == 0:
                        return True
    return False


# ---------------------------------------------------------------- construction
def circulant_build(E, D, P):
    """H_X block (i, l) = C(E[i][l]), H_Z block (j, l) = C(D[j][l]) over Z_P.

    C(s): row r has its 1 at column (r - s) mod P.
    """
    J, L = E.shape
    HX = np.zeros((J * P, L * P), dtype=np.int8)
    HZ = np.zeros((J * P, L * P), dtype=np.int8)
    for i in range(J):
        for ell in range(L):
            for r in range(P):
                HX[i * P + r, ell * P + (r - int(E[i, ell])) % P] = 1
                HZ[i * P + r, ell * P + (r - int(D[i, ell])) % P] = 1
    return HX, HZ


def random_matching(rng, L, sig):
    """Perfect matching of the L block columns with NO pair {l, sigma(l)}.

    A sigma-pair {u, sigma(u)} in M_ij forces E[i][u]-E[j][u] to repeat at
    sigma(u) (pairing equation with v = sigma(u)), a guaranteed 4-cycle, so
    such matchings can never pass the cycles4 screen.
    """
    while True:
        perm = rng.permutation(L)
        pairs = [(int(perm[t * 2]), int(perm[t * 2 + 1])) for t in range(L // 2)]
        if all(v != sig(u) and u != sig(v) for u, v in pairs):
            return pairs


def make_sig(kind, L):
    if kind == "shift":
        return lambda ell: (ell + L // 2) % L  # noqa: E731
    if kind == "reflect":
        return lambda ell: (-ell) % L  # noqa: E731
    raise ValueError(kind)


def draw_instance(J, L, P, seed, eta=-1, keep=1, min_k=1, verbose=True, max_attempts=4000, sig=None, cycles=False):
    """Sample CSS parents satisfying the halving criterion; yield (E, D, HX, HZ, k).

    rho = id, alpha = beta = 0.  Each attempt draws fresh pair partitions; an
    instance is kept only if neither H_X nor H_Z has two equal columns (that
    structural weight-2 kernel vector kills distance).  ``cycles=True`` adds
    the (3,8) fieldnote's no-4-/6-cycle screen on both exponent arrays.
    """
    rng = np.random.default_rng(seed)
    if sig is None:
        sig = make_sig("shift", L)
    cells = [(i, j) for i in range(J) for j in range(i + 1, J)]

    D = np.empty((J, L), dtype=np.int64)
    found = []
    attempts = 0
    match_cache = None
    while len(found) < keep and attempts < max_attempts:
        attempts += 1
        if match_cache is None:
            matchings = {}
            for i, j in cells:
                matchings[(i, j)] = random_matching(rng, L, sig)
            rows = []
            for (i, j), pairs in matchings.items():
                for u, v in pairs:
                    row = np.zeros(J * L, dtype=np.int64)
                    row[i * L + u] += 1
                    row[i * L + v] -= 1
                    row[j * L + sig(u)] += -eta
                    row[j * L + sig(v)] += eta
                    rows.append(row)
            basis, nullity = nullspace_mod(np.array(rows) % P, P)
            match_cache = (basis, nullity)
            if verbose:
                print(f"  system: {len(rows)} eqs, {J * L} vars, nullity {nullity}")
        basis, _nullity = match_cache

        coeff = rng.integers(0, P, size=len(basis))
        vec = sum((int(c) * b for c, b in zip(coeff, basis)), np.zeros(J * L, dtype=np.int64)) % P
        E = vec.reshape(J, L)
        if E.sum() == 0:
            continue
        for j in range(J):
            for ell in range(L):
                D[j, ell] = (eta * int(E[j, sig(ell)])) % P
        HX, HZ = circulant_build(E, D, P)
        if not verify_css(HX, HZ):
            raise RuntimeError("pairing equations did not give CSS orthogonality")
        if (
            equal_col_pairs(HX)
            or equal_col_pairs(HZ)
            or (cycles and (cycles4(E, P) or cycles4(D, P) or cycles6(E, P) or cycles6(D, P)))
        ):
            # degenerate solution space or (optionally) short cycles: fresh
            # matching
            if attempts % 25 == 0:
                match_cache = None
            continue
        k = compute_k(HX, HZ)
        if k < min_k:
            continue
        found.append((E.copy(), D.copy(), HX, HZ, k))
        match_cache = None  # diversity: next kept draw gets a fresh matching
    if not found:
        raise RuntimeError("no instance with k >= min_k found")
    if verbose:
        print(f"  drew {len(found)} instance(s) after {attempts} attempts")
    return found


def halve(HX, HZ, J, L, P, eta=-1, sig=None):
    """Fold the parent under pi(l, t) = (sigma(l), eta t); return S = (A | B)."""
    if sig is None:
        sig = make_sig("shift", L)
    n_parent = L * P
    partner = np.empty(n_parent, dtype=np.int64)
    for ell in range(L):
        for t in range(P):
            partner[ell * P + t] = sig(ell) * P + (eta * t) % P
    assert np.array_equal(partner[partner], np.arange(n_parent)), "pi not involutive"
    firsts = [q for q in range(n_parent) if q < partner[q]]
    assert 2 * len(firsts) == n_parent, "pi has fixed points"
    order = firsts + [int(partner[q]) for q in firsts]
    n = len(firsts)

    HXo = HX[:, order]
    A, B = HXo[:, :n], HXo[:, n:]

    # row (j, b) of H_Z <-> row (j, a) of H_X with b = eta*a (alpha = 0)
    rp = np.empty(J * P, dtype=np.int64)
    for j in range(J):
        for b in range(P):
            rp[j * P + b] = j * P + (eta * b) % P
    HZ2 = HZ[rp][:, order]
    if not (np.array_equal(HZ2[:, :n], B) and np.array_equal(HZ2[:, n:], A)):
        raise RuntimeError("halving identity H_Z' = (B | A) failed")

    S = np.concatenate([A, B], axis=1)
    isotropy = bool(((A @ B.T + B @ A.T) % 2).any())
    if isotropy:
        raise RuntimeError("A B^T + B A^T != 0: not a valid stabilizer code")
    k_fold = n - rank(S)
    return S, n, k_fold


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--J", type=int, default=4)
    ap.add_argument("--L", type=int, default=10)
    ap.add_argument("--P", type=int, default=19)
    ap.add_argument("--eta", type=int, default=-1, choices=[1, -1])
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--draws", type=int, default=4, help="CSS parents to sample")
    ap.add_argument("--trials", type=int, default=4000, help="RIS trials per screen call")
    ap.add_argument("--ris-seconds", type=float, default=90.0, help="wall-clock cap per RIS screen call")
    ap.add_argument("--out", default=OUT_DIR)
    ap.add_argument("--sigma", default="shift", choices=["shift", "reflect"])
    ap.add_argument("--cycles", action="store_true", help="also enforce no-4-/6-cycle screens on E and D")
    args = ap.parse_args()

    J, L, P = args.J, args.L, args.P
    assert P > 2 and L % 2 == 0
    tag = f"{J}_{L}_{P}"
    os.makedirs(args.out, exist_ok=True)
    sig = make_sig(args.sigma, L)

    print(f"(J,L,P)=({J},{L},{P}) eta={args.eta} sigma={args.sigma}: parent n={L * P}, fold n={L * P // 2}")
    found = draw_instance(J, L, P, args.seed, eta=args.eta, keep=args.draws, min_k=1, sig=sig, cycles=args.cycles)

    for idx, (E, D, HX, HZ, k) in enumerate(found):
        S, n_f, k_f = halve(HX, HZ, J, L, P, eta=args.eta, sig=sig)
        w = int(np.logical_or(S[:, :n_f], S[:, n_f:]).sum(axis=1).max())
        dX, _ = hd.ris_min_logical(HX, HZ, trials=args.trials, seed=args.seed, max_seconds=args.ris_seconds)
        dZ, _ = hd.ris_min_logical(HZ, HX, trials=args.trials, seed=args.seed, max_seconds=args.ris_seconds)
        dP, _ = hd.ris_min_pauli_logical(
            S[:, :n_f], S[:, n_f:], trials=args.trials, seed=args.seed, max_seconds=args.ris_seconds
        )
        d_par = min(dX, dZ)
        w_par = int(max(HX.sum(axis=1).max(), HZ.sum(axis=1).max()))
        p_path = os.path.join(args.out, f"halved_{tag}_d{idx}_parent.npz")
        f_path = os.path.join(args.out, f"halved_{tag}_d{idx}_fold.npz")
        np.savez(p_path, hx=HX, hz=HZ)
        np.savez(f_path, s=S)
        print(f"  draw {idx}: parent [[{L * P},{k},<={d_par}]] w<={w_par} -> {p_path}")
        print(f"  draw {idx}: fold   [[{n_f},{k_f},<={dP}]] w={w} -> {f_path}")


if __name__ == "__main__":
    main()
