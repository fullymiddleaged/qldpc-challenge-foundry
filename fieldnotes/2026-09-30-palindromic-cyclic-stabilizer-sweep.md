---
title: "One-block palindromic cyclic stabilizer codes: the whole sweep, and why its seven finds were all CSS behind a Hadamard"
date: 2026-09-30
author: "@MathysRennela"
model: "Claude Fable 5.1"
topics: [stabilizer, cyclic, palindromic, group-algebra, clifford-relabel, board-hygiene]
---

**This sweep is @FarLab's work, not the reclassifier's.** All seven entries carry
`provenance.authors: ["@FarLab"]`, and @FarLab is the git author of all seven
note commits: #2470 (`48-10-6-s`), #2480 (`154-34-6`), #2496 (`110-30-5`),
#2501 (`130-34-5`), #2549 (`52-4-7`), #2552 (`28-4-5-s`) and #2577
(`182-38-6`). Everything below — the construction, the sweep, the numbers, the
dead ends, the tooling and the witnesses — is @FarLab's, carried over from those
seven notes. What this reclassification adds is only the observation that binds
them together: the reason none of the seven is on the board any more is a
property they all share.

Reclassified from the submission notes of those seven board entries, removed in
`91d689a4` and `836566c9`: `110-30-5`, `130-34-5`, `154-34-6`, `28-4-5-s`,
`48-10-6-s`, `52-4-7` and `182-38-6`. The seven notes were 96-98 percent
identical, differing only in the title and the reproduction line, so they
collapse to one account. The construction and the sweep are worth keeping; the
outcome is a caution about which board a code belongs on.

## Direction & hypothesis

The general-stabilizer board's low-weight cells (w <= 8) held only a handful of
small codes and the two Chamon baselines (`kd^2/n = 4`), while the higher-k
entries sat at weight 32. The cheapest genuinely non-CSS LDPC family is the
one-block cyclic code `S = (circ(a) | circ(b))`: one generator `X^a Z^b` and
its n cyclic shifts on n qubits, a Camara-Ollivier-Tillich-type group-algebra
code (arXiv:quant-ph/0502086) over `Z_n`.

The shifts commute iff `a b* + b a* = 0` (where `b*(x) = b(1/x)`); that holds
whenever `a = g p` and `b = g q` with `p, q` palindromic (`p = p*`) and `g`
arbitrary, and for even n also after `b -> x^(n/2) b`. Then
`k = deg gcd(a, b, x^n - 1)` and the check weight is `|supp a u supp b|`.
Hypothesis: sparse palindromic pairs would fill the weight-4/6/8 cells at
`n = 18-200` with codes at `kd^2/n` well above the existing entries.

## What was searched

- All n from 17 to 200; `p, q` palindromic of weight <= 4 (`p` taken up to the
  multiplier group `Z_n^*`), `g in {1} u {1 + x^j}`, both shift variants for
  even n, check weight <= 8. `k >= 1` by polynomial gcd — most pairs have
  `k = 0`; per n the 6000 highest-k, lowest-w candidates were screened.
- Screening cascade with the verifier's own Pauli-weight RIS
  (`verify/heuristic_distance.py`): 8 trials (drop `d < 5`, and `kd^2/n < 3`
  for `n > 40`), then 300 trials; kept only if not dominated on
  `(n, k, d, w)` by any stabilizer board entry or earlier find.
- Deep pass on every Pareto point: 20,000 RIS trials on each of two fresh
  seeds.

## Evidence trail

The ladder ran identically for all seven; the distance it reached differs per
code, and the `d` column of the table below is the one each entry was filed at.

- Screen: `d <= 5` at 300 trials for the four entries filed at `d = 5`,
  `d <= 6` for `154-34-6`, `48-10-6-s` and `182-38-6`, `d <= 7` for `52-4-7`.
- Deep: 20,000 trials on each of two fresh seeds found nothing lighter in any of
  the seven (held).
- `qldpc submit` re-searched at 20,000 trials; the filed witness has the Pauli
  weight given in the `d` column, and its explicit support is the `witness`
  column.
- Claim: `d` as in the table, witness-backed upper bound; the board's certifier
  does not minimise Pauli weight, so no stabilizer entry is exact on the board.

The seven survivors, with their reproduction parameters and filed witnesses
(generator `i` is X on `{i + e mod n : e in supp a}` and Z on
`{i + e mod n : e in supp b}`, both to Y where they meet;
`k = n - rank_2(S)`; witness supports are as recorded in each filed note):

| slug | n | k | d | a(x) | b(x) | witness (X on / Z on) |
|---|---|---|---|---|---|---|
| `28-4-5-s` | 28 | 4 | 5 | `x^12 + x^16` | `x^7 + x^13 + x^15 + x^21` | `[5, 13, 25]` / `[16, 22]` |
| `48-10-6-s` | 48 | 10 | 6 | `1 + x^20 + x^24 + x^28` | `x^19 + x^21 + x^27 + x^29` | `[14, 18, 26, 34]` / `[13, 15]` |
| `52-4-7` | 52 | 4 | 7 | `x^24 + x^28` | `x^13 + x^19 + x^33 + x^39` | `[9, 25, 41]` / `[14, 20, 30, 36]` |
| `110-30-5` | 110 | 30 | 5 | `x^39 + x^49 + x^61 + x^71` | `x^18 + x^48 + x^62 + x^92` | `[7, 29, 51, 73, 95]` / `[]` |
| `130-34-5` | 130 | 34 | 5 | `x^47 + x^57 + x^73 + x^83` | `x^24 + x^54 + x^76 + x^106` | `[29, 107]` / `[58, 68, 78]` |
| `154-34-6` | 154 | 34 | 6 | `x^59 + x^73 + x^81 + x^95` | `x^34 + x^76 + x^78 + x^120` | `[20, 34]` / `[17, 37, 59, 149]` |
| `182-38-6` | 182 | 38 | 6 | `x^71 + x^85 + x^97 + x^111` | `x^44 + x^86 + x^96 + x^138` | `[47, 61]` / `[20, 46, 62, 88]` |

Note that the witness is a mixed Pauli operator in every case except `110-30-5`,
whose witness is X-only — the stabilizer code has one logical of each type, and
which side a search lands on depends on the seed, not on the code.

## Why none of them is on the board any more

**Every one of the seven is a CSS code with a Hadamard on half the qubits.**
Filed as general-stabilizer entries, they were removed when the rule landed
that such a code must be typed CSS (as that image) — the verifier rejects it
otherwise. Their CSS images were then checked and found **dominated on the CSS
board, or tying `28-4-5`**, so nothing was re-filed.

This is the transferable part: the sweep did not fail at finding codes, it
failed at placing them. A candidate recovered as `S = (circ(a) | circ(b))` with
both `a` and `b` real is one Clifford-twist away from a CSS code, and the CSS
board already had the equivalent points. **Screen for CSS-behind-local-Clifford
before filing on the stabilizer board**, because a code that survives the
distance search and loses the placement rule costs a full CI refutation pass
and leaves no entry behind.

## Dead ends

- Two-generator and non-palindromic pairs were not searched; the family is
  exactly the one above.
- Composite `n` (63, 68, ...) produce 50-80k `k >= 1` candidates and dominate
  the run time; prime `n` produce few and mostly `k = 1` codes (the
  `X_{i,i+1} Z_{i-j,i+j}` weight-4 family, `d` up to 11 at `n = 61`).
- Candidates that screened at 300 trials and dropped under the deep pass were
  discarded; the deep pass moved a minority of finds by one or two.

## Tools

@FarLab's, as recorded in all seven original notes: Claude Fable 5.1 in Claude
Code; numpy; `verify/heuristic_distance.py` (Pauli-weight RIS),
`cli/qldpc.py`, `verify/`. @FarLab ran the whole sweep on a laptop CPU over a
few hours. The reclassification added no compute.

## Reproduction

Build `S = (circ(a) | circ(b))` for any row of the table above, take
`k = n - rank_2(S)`, and check the commutation condition
`a b* + b a* = 0`. Witness search with Pauli-weight RIS at 20,000 trials on
each of two fresh seeds, as described above.
