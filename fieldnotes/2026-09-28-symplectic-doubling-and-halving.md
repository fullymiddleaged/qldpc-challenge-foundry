---
title: "Symplectic doubling and its inverse: turning a non-CSS code into its CSS double"
date: 2026-09-28
author: "@MathysRennela"
model: "Mimo-V2.6-Flash (opencode agent harness)"
topics: [symplectic-doubling, pair-partition, cpm-codes, calibration, negative-results]
related:
  - 2026-09-01-multiband-dense-packing-method.md
---

## Summary

One screening pass over the CPM pair-partition family of arXiv:2609.30069 produced
two board entries from a single construction: a non-CSS (stabilizer) code
[[145,32,6]] (k = 32, d <= 6, kd^2/n = 7.945) and its CSS double [[290,64,6]]
(k = 64, d <= 6, kd^2/n = 7.945). The pair sits in the equality case: the
parent's weight-6 distance witness (X on {0, 16, 88, 91}, Z on {58, 68}) is
Y-free — X and Z supports disjoint — so its doubled image is a weight-6 CSS
logical of the double, and the gate's d_X = d_Z = 6 witnesses on [[290,64,6]]
confirm d' = d. Both passed the trusted gate
(`verify/validate_candidate.py`: passed, board_advancing, empty dominator list; an
8000-trial refutation found nothing lighter). This note records the operation, its
inverse, and — explicitly — what is known versus what we only characterised.

## What this is, and is not

**Not a discovery of the operation.** Symplectic doubling is standard, and the
CSS -> non-CSS (halving) direction is the content of Prop. 4 of arXiv:2609.30069.
The repo's own verifier already uses the doubling: `verify/heuristic_distance.py`
searches the doubled CSS code H'_X = (A | B), H'_Z = (B | A). We did not invent it.

**What we found out (calibration-grade).** The operational cost of *using* it in
this family: which draws survive, that girth optimisation is incompatible with the
recipe's pairing, that the fold is genuinely non-CSS yet still halves out of a CSS
parent, and that the doubling is efficiency non-decreasing (kd^2/n >= the
parent's, with equality iff a minimum-weight logical is Y-free). Read the
numbers below as the boundary of this route, not as a new theorem.

## The operation: symplectic doubling (always available)

Take any stabilizer code S = (A | B), isotropic (A B^T + B A^T = 0). Define

    H'_X = (A | B)      H'_Z = (B | A)

Then H'_X H'_Z^T = A B^T + B A^T = 0, so the result is CSS. Both sides have rank
rank(S), so k' = 2n - 2 rank(S) = 2k. The distance is *not* carried across
unchanged in general: a Y occupies one qubit but two bits in the doubled vector,
so Hamming(doubled logical) = PauliWeight + (#Y qubits), giving d' >= d, with
equality iff some minimum-Pauli-weight logical is Y-free. Hence n -> 2n,
k -> 2k, d' >= d, and kd^2/n is non-decreasing — invariant exactly in the
Y-free case. This direction has no preconditions — every stabilizer code doubles
to a CSS code at twice the blocklength — which is why the non-CSS side of a
family always has a CSS twin twice the size, at *at least* the same efficiency
(and strictly better when no minimum-weight logical is Y-free).

## The inverse: halving (the constrained direction)

In the CPM pair-partition family a CSS code is fixed by (J, L, P) and exponent
arrays E, D over Z_P, with block (i, l) = circ(E[i][l]) in H_X and block
(j, l) = circ(D[j][l]) in H_Z. It halves only when

- D[j][l] = eta * E[j][sigma(l)] (rho = id, alpha = beta = 0 gauge; eta = +/-1),
- sigma is a fixed-point-free involution of the L columns,
- each off-diagonal difference array E[i][l] - eta * E[j][sigma(l)] takes every
  value an even number of times (the pair-partition equations).

That last condition is a *linear* system over F_P, so the family is sampled by
choosing a perfect matching per off-diagonal cell and exponentiating over its
nullspace basis. The halving map pi(l, t) = (sigma(l), eta * t) reorders the L*P
columns so H_X -> (A | B) and H_Z -> (B | A); S = (A | B) is the fold, with
k_fold = k_parent / 2 and isotropy A B^T + B A^T = 0 checked during the build. A
The builder is `research/build_halved_pp.py`, included in this PR: it samples the family from the pair-partition equations, screens each draw, halves the parent, and writes both matrices.

## The sweep

(J, L, P) = (4, 10, 29), sigma(l) = (l + 5) mod 10 (shift by L/2), eta = -1,
rho = id, alpha = beta = 0. Draws are sigma-pair-free matchings (a pair
{u, sigma(u)} provably forces a 4-cycle), screened for equal columns, then RIS at
thousands of trials per draw under a wall-clock cap. Draw 6 gave the best pair:
parent [[290,64,<=6]] and fold [[145,32,<=6]], with k = 64 and 32. The submitted
matrices reconstruct bit-for-bit from E and D (rebuilt H_X/H_Z equal the submission
exactly), and k matches the paper's scaling k_fold = k_parent / 2.

## Dead ends and boundaries

- Equal columns: **87% of all draws** have repeated columns in H_X/H_Z, giving a
  weight-2 kernel vector and d <= 2. Screening for this is mandatory; naive draws
  almost never survive.
- Short cycles: any matching pairing {l, sigma(l)} forces a repeated row
  difference. With the no-4-/6-cycle filters on, **0 of 164** matchings survived,
  and cycles4(E)-freedom failed in **61 of 61** equal-column-clean draws. Girth is
  not reachable under the sigma-pairing, so the filters stayed off (girth is not
  part of the gate).
- Scaling P to chase the record family fails: the P = 47 lift gave k = 100, not the
  ~109 required, so the cell threshold (d > 14) was out of reach while the screen
  only reached 4.
- Size (4, 12), P = 31 (n = 372) is dominated: the board already holds
  [[372,130,16]] from the related pair-partition CPM family of arXiv:2607.14091, so
  at check weight 12 the bar is d > 17.

## Non-CSS is not the same as "half of a CSS code"

The fold [[145,32,6]] is **not** CSS up to local Hadamards — the gate's
equivalence check found no Hadamard subset making every generator pure, so it is
not "secretly CSS". Yet it doubles to a CSS code. The doubling route is therefore
strictly different from relabelling qubits to expose a CSS structure: a code can be
irreducibly non-CSS under local Cliffords and still be the half of a CSS code. That
distinction is the main conceptual takeaway and the reason the two boards are
related at all.

## CSS in, CSS out: the map degenerates on a CSS input

Doubling a CSS code to get another CSS code does not help — the output is CSS
(the map always is), but it is degenerate. A CSS code in its canonical form has
pure generators, `A = [H_X; 0]` and `B = [0; H_Z]`, so the double is a block
direct sum with no check touching both halves: two disconnected copies of the
original. Checked directly, the [[7,1,3]] Steane code doubles to a CSS code whose
combined X/Z Tanner graph has **2 components**, whereas the genuinely non-CSS
[[5,1,3]] code doubles to a **connected** CSS code. The disconnected result is
inadmissible (the verifier requires a single connected component) and
efficiency non-decreasing (`n -> 2n`, `k -> 2k`, `d' >= d` with equality iff a
minimum-weight logical is Y-free, so `kd^2/n` never drops), and
the CSS -> halve -> double round trip returns the same code. So the map is one-way:
it only does anything on a *non-CSS* input, and connectivity of the double is the
signature of genuine non-CSSness.

## Reopen conditions

This route is bounded to pair-partition CPM codes with the eta/sigma structure and
a *solvable* even-multiplicity system. Reopen with a genuinely new mechanism: other
involutions sigma (we only used the L/2 shift), gauges other than alpha = beta = 0,
or non-circulant pair partitions (e.g. non-abelian groups), where the matching
system is no longer linear over F_P. Stop a variant if, as here, the cycle screen
kills essentially every matching and scaling P cannot reach the cell's bar.

## Provenance

Findings and the search are from the [[145,32,6]] / [[290,64,6]] submission
session (model Mimo-V2.6-Flash, opencode harness); the write-up was assembled
afterwards. The family was reconstructed independently because the paper's data
repository (github.com/ultra-high-rate-qec/design-principles-data) was empty when
searched, so every code here is our own draw of the published recipe. All distances
are witness-backed upper bounds, not exact certificates.
