# Spec (Job B): growth of the two-site bracket counts with the boxes — [2,1,0], [3,1,0], [3,2,0], [4,1,0], [4,2,0] (2026-10-11)

Delegated to a wolfbook-builder subagent (sonnet); planned and reviewed by the main session. **New notebook:**
`Paul/Mathematica/new-conventions/L2-bracket-count-growth.wb` (rep-generic with `Λ`; build/run as the previous stages:
JSON, `code` CLI, sources in K1, printed results as markdown cells; say so). Paul's acceptance criterion: a Q-number's cost
must scale with the boxes, not with dim V. The two-site block formula (Stage B) has a count that is polynomial in
(λ₁,λ₂) by construction — (λ₂+1)² blocks, each a matrix on C(λ₁+2,2) ladder minors — but the constants matter (54 for
[3,1,0] vs dim V² = 225). This job measures the actual counts over five reps so the growth law is empirical, not argued.

Read first: HANDOFF top section; the Stage B notebook `L2-two-site-rotation-brackets.wb` (cells 2–12 setup, 20–25 K tensor
and tools, 37–41 blocks and counts, 47–49 the (∞,∞) family) and Stage C cells 26–30 (pivot scan) and 68 ([3,1,0] scan +
greedy). Dumps: `ref-two-site-all.txt`, `ref-general-all.txt` (cells 3–6 rep-generic setup, 9–10 labels/monomials),
`ref-310-two-rotations-all.txt` (dim-15 setup); for [4,1,0]/[4,2,0] the exact transfer matrices follow from the same
rep-generic code (`JJ[k,l,rep]`, `tL2`); if the 729-dim eigensystem at prec 60 is slow, use prec 40 and say so (the
block ranks are exact rationals after rationalisation — check the rationalisation error).

## 0. Reps and sizes

| Λ | dim V | dim V² | |Λ| | ladder points λ₁+2 | nMin = C(λ₁+2,2) | blocks (λ₂+1)² | one-site counts (known) |
|---|---|---|---|---|---|---|---|
| [2,1,0] | 8 | 64 | 3 | 4 | 6 | 4 | 5 / 4 |
| [3,1,0] | 15 | 225 | 4 | 5 | 10 | 4 | 7 / 5 |
| [3,2,0] | 15 | 225 | 5 | 5 | 10 | 9 | 11 / 9 |
| [4,1,0] | 24 | 576 | 5 | 6 | 15 | 4 | ? |
| [4,2,0] | 27 | 729 | 6 | 6 | 15 | 9 | ? |

Conventions as in Stage B (companion frame, h = 1/2, θ = (2/7,−3/5), z = (3/2,5/4,2/3)); w = (5/2,7/3,3/5),
w′ = (7/4,−2/5,9/7); (a,a′) = (5/7,−4/9); generic pair g₁ = (2/3,−3/7,5/4), g₂ = (−1/2,4/9,3/5).

## 1. Tasks (per rep, in the order of the table; each rep one block of cells; a `Grid` per rep; save after each rep)

**Task 0 — setup per rep.** Generators, exact sparse transfer matrices, eigenstates (time the eigensystem: for dim 576/729
at prec 60 expect minutes; if > 10 min, prec 40), dual Q's (degrees μ_j; Baxter residual), SoV basis of B (left
eigenvectors; label rule (a|bc) check: count = dim V per site, vacuum label, all pairs distinct), wavefunction monomial
check (rows {2,3}, k = 3; max deviation), N_m. Report the sizes and times.

**Task 1 — K tensor and blocks.** The two-site K tensor (dim V² × dim V² linear solve; rationalise; error; denominators;
z-independence only for the two new reps at one z′ if the solve is < 2 min). Blocks A^{(c,c′)} in the Stage B convention.
For each Φ in {generic pair; random covector; vacuum Ωc[w]⊗Ωc[w]; Ωc[w]⊗Ωc[w′]; (∞,∞) at (a,a′); c₁(g₁)⊗e₁ᵀ and e₁ᵀ⊗c₂(g₂)
(one-site regressions: 5/7/11 for the first three reps, and the new one-site numbers for [4,1,0], [4,2,0])}: per block
rank r, the product count P for the default pivots in both orientations (min), and a bounded pivot scan (all pivot sets
if C(nMin, r) ≤ 3000 per block, else 300 random pivot sets per orientation); Σr and ΣP (best found); the overlap identity
with explicit bracket products on all states (10⁻³⁸ or better, one Φ per rep is enough, plus the vacuum).
Also the automatic upper bound per block (Σ over a rank decomposition of ⌈rank/2⌉·⌈rank/2⌉ worst case) and the trivial
bound (λ₂+1)²·nMin·⌈(λ₁+2)/2⌉², for the table.

**Task 2 — reach of the (∞,∞) family for the new reps (exact, finite twist).** The undone covector of
Ωc[(w₁,∞)]⊗Ωc[(w₁′,∞)] in the diagonal frame (Stage D cell 21 machinery: e₁ᵀρ(n(z,w)n(z,w′)⁻¹) with w₂ = w₂′ = ∞,
w₁ ≠ w₁′): weights reached and ranks on multi-dimensional weight spaces, vs the slice weights of Λ⊗Λ (list the HW states
of Λ⊗Λ for [3,2,0], [4,1,0], [4,2,0] — compute them, don't type them). Is the family CG-complete for every rep?

**Task 3 — the growth table and fits.** One table: rep, |Λ|, dim V, dim V², one-site count (generic/best), two-site Σr and
ΣP for: generic pair, vacuum (det₆), Ωc[w]⊗Ωc[w′], (∞,∞); the trivial bound. Then: (i) is ΣP(generic) a polynomial in
(λ₁,λ₂) of low degree? Fit the five values with the monomials {1, λ₁, λ₂, λ₁², λ₁λ₂, λ₂²} (6 unknowns, 5 data: report
the 1-parameter family or fix one monomial; also try the form (λ₂+1)²·p(λ₁) with p of degree ≤ 2); (ii) the same for
Σr(generic), ΣP(vacuum), ΣP((∞,∞)); (iii) the ratios ΣP/dim V² and ΣP/|Λ|². Say plainly whether the counts grow like a
polynomial in the boxes of degree ≤ 4 with the constants observed, and how they compare with dim V².

**Task 4 — markdown summary** `### bracket-count growth with the boxes — result (2026-10-11)`: the table, the fits, the
reach results, the one-site numbers for the new reps, timings per rep, anything that did not finish (e.g. if [4,2,0] is
too slow, stop after [4,1,0] and say so). Save; no HANDOFF.

## 2. Rules and gotchas

As before (K1 shared; first cell wipes; no `*)` in comments; no `D`,`E`,`N`,`I`,`C`,`K` symbols; time one state / one
block before loops; `TimeConstrained`; save after each rep; no other notebooks, no `paper/*.tex`, no git, no HANDOFF).
Memory: a 729×729 dense matrix at prec 60 is ~40 MB — fine; keep the K tensor as a rational `SparseArray`. The B operator
eigensystem for dim 729: use `Eigensystem` on the numeric matrix (dense); if the spectrum is degenerate at one u, use
B(u₁) + πB(u₂). The label rule may need checking for λ₁ = 4 (the general notebook verified it only up to λ₁ = 3): if a
triple is not of the form {−1−a,−b,−c}, record the actual label set and adapt the monomial rule (the wfTest will tell).
