# gl(4) L=2 untwisted transfer-matrix eigenvalues via Casimirs — design

**Date:** 2026-07-19
**Subproject:** XXX (gl(4))
**Notebook:** `Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb` — new section appended after the
existing twisted acceptance run (cells 1–34, commit `fcd1349`)
**Status:** design approved 2026-07-19

## Goal

In a new section of the existing gl(4) L=2 notebook, find the **closed-form eigenvalues `Λ_a(u)`**
(`a=1,2,3,4`, one per fundamental transfer matrix `τ1,τ2,τ3,qdetT`) of the **untwisted** (`z[k]=1`
for all `k`) two-site transfer matrices, for the two-row rectangular rep family
`λ1=[S,S,0,0]`, `λ2=[P,P,0,0]`, expressed in terms of **su(4) Casimir invariants** of the irreps `μ`
appearing in the tensor decomposition `λ1⊗λ2 = ⊕_μ N_μ · μ`.

This is the deliberate `z→1` follow-on flagged in `HANDOFF.md` after the twisted-case build
(`docs/superpowers/plans/2026-07-19-su4-l2-twisted-baxter.md`), now scoped to the tractable
two-row-rectangular rep family rather than fully general reps.

## Why this is a different problem from the twisted case

At `z[k]=1` the diagonal twist that broke su(4) symmetry is gone, so the untwisted transfer matrices
commute with the **full** su(4) (equivalently gl(4)) symmetry algebra acting on the two-site space.
By Schur's lemma:

- On a **multiplicity-free** irrep block (`N_μ=1`), the transfer matrix acts as a pure scalar — this
  is exactly the regime where a clean `Λ_a(u) = f(Casimir(μ), Casimir(λ1), Casimir(λ2); u)` closed
  form can hold, generalizing the familiar su(2) Heisenberg-chain relation
  `H = P − 1 ∝ C2(total) − C2(1) − C2(2)`.
- On a block where `μ` occurs with multiplicity `N_μ>1`, the transfer matrix acts as a genuine
  `N_μ×N_μ` matrix on the multiplicity space — states can share identical weight *and* identical
  Casimir invariants yet carry different transfer-matrix eigenvalues. No pure-Casimir formula can
  distinguish them; that would require finer data (Racah/6j-type coefficients), which is out of scope
  here.

Consequently `(H1,H2,H3)` weight labels alone (sufficient for the twisted case, since distinct twist
eigenvalues fully lift degeneracy) are **not** sufficient at `z=1`: different `μ`'s can overlap in
weight space, and states within one `μ` already have interior weight multiplicity. Resolving the
spectrum requires genuine su(4) Casimir operators.

**Honest scope acknowledgement:** whether `λ1⊗λ2=[S,S,0,0]⊗[P,P,0,0]` is multiplicity-free (for the
swept `(S,P)` range) is an empirical question this work answers, not an assumption. If some blocks
have `N_μ>1`, that is a reportable finding (in the same spirit as the existing Section M report on
the θ1-tower local-rotations ansatz) — the deliverable is a **precise, verified** closed form on
whatever fraction of the spectrum it actually covers, not a forced-fit formula.

## Existing state reused (no changes to committed cells 1–34)

- GT engine `ps/dim/J0/Jp/Jm/JJ`, `Ee/id/EE/idd` (N-general, already gl(4)).
- Non-memoized, z-read-live building blocks: `L`, `T`, `Tg`, `ν`, `ν1`, `qdetT`, `H1op/H2op/H3op`.
  These need **no** changes — they recompute correctly the instant `z[k]` is reassigned.
- Memoized, z-dependent building blocks that **do** need clearing at the section boundary: `qmin`,
  `qming`, `tg`, `tgCoeffs` (all memoized on `(λ1,λ2)`, not on `z` — a later call with different `z`
  would silently return the stale twisted-case value without an explicit `ClearAll`).

## Section boundary

A single, clearly commented cell marks the transition from the twisted-verification section to the
untwisted-exploration section:

```wolfram
(*Section boundary: switch to the untwisted (z=1) regime for the Casimir-eigenvalue exploration
  below. qmin/qming/tg/tgCoeffs are memoized on (λ1,λ2) but read z live inside their definitions —
  without this ClearAll they would silently return stale twisted-z values. T/Tg/ν/ν1/qdetT/H1op/
  H2op/H3op are NOT memoized and need no clearing; they recompute correctly automatically.*)
ClearAll[qmin, qming, tg, tgCoeffs];
z[1] = 1; z[2] = 1; z[3] = 1; z[4] = 1;
{χ1, χ2, χ3, χ4}
```

(`χ1..χ4` recompute automatically to `{4,6,4,1}` — the elementary symmetric functions of four 1's —
confirming the reassignment took hold.)

A fresh-kernel top-to-bottom run of the **whole** notebook after this work therefore first reproduces
the twisted acceptance report (cells 1–34, `z` at 70-digit twisted values), then crosses the boundary
cell, then runs the untwisted section — both halves correct in one linear pass, cell order carrying
the meaning (no out-of-order-evaluation footgun).

## New machinery

### 1. Casimir operators

**Quadratic Casimir** (needed first; cheap — reuses `EE`):

```wolfram
ClearAll[C2op];
C2op[λ1_, λ2_] := C2op[λ1, λ2] = Sum[EE[λ1, λ2][i, j] . EE[λ1, λ2][j, i], {i, 4}, {j, 4}];
```

**Cubic/quartic Casimirs** (`C3op`, `C4op` — analogous triple/quadruple traces of `EE`) are built
**only if** the multiplicity check (below) shows `C2` alone fails to resolve all degeneracies within
some weight sector for the swept `(S,P)` range. Not built speculatively.

### 2. Casimir-resolved eigensystem

A new routine (own name — no collision with `TauEigensystem4`, so both remain independently callable
in the same kernel session) generalizing `TauEigensystem4`'s proven "generic linear combination of
mutually-commuting operators, diagonalize once" technique. The operator family being mixed changes
from `{τ-coefficient matrices, H1,H2,H3}` (right for the twisted case) to
`{C2op (and C3op/C4op if needed), H1op,H2op,H3op}` (right for the untwisted case, since at `z=1` the
`τ`-coefficient matrices themselves are among the commuting family too and can optionally be folded
in, but the Casimir+Cartan family is what actually separates different `μ`'s sharing weight space).

From each resolved simultaneous eigenvector: Rayleigh-quotient against `τ1,τ2,τ3,qdetT`'s coefficient
matrices to get `Λ_a(u)` (same technique as `tau1/tau2/tau3` accessors), and against `C2op` (`C3op`,
`C4op` if built) to get the Casimir label(s).

### 3. Multiplicity check

For each `(S,P)` in the sweep: after resolving states by `(H1,H2,H3,C2,...)`, check whether any
`(H1,H2,H3,C2,...)` sector contains more than one linearly-independent eigenvector with **different**
`Λ_a(u)` (a direct symptom of `N_μ>1` on that block). Record per-`(S,P)` whether the decomposition is
fully multiplicity-free.

### 4. Fit

On multiplicity-free blocks: tabulate `Λ_a(u)` vs. `(S,P,C2(μ))` (and `C3`,`C4` if used) across the
sweep; propose and hard-assert a closed-form candidate generalizing the Heisenberg-type relation to
the full `u`-dependent, gl(4), two-site case. If any block fails to be multiplicity-free, report that
finding precisely (which `(S,P)`, which sector, the actual dimension of the multiplicity space) rather
than omitting it.

## Sweep range

Small `(S,P)` pairs, escalating only as needed to distinguish candidate closed forms and to probe for
multiplicity: start with `(1,1), (2,1), (1,2), (2,2)`. `dim[[S,S,0,0]]` grows with `S` (Weyl dimension
formula for gl(4)); the actual per-pair tensor-product dimension will be read from the existing `dim`
function during execution and used to judge whether to extend the sweep further, rather than fixed in
advance.

## Verification convention

Same as the rest of the notebook: hard-assert with `If[worst > tol, Print[...]; Abort[]]` where a
closed form is claimed; explicit, precise reporting (not silent omission) where a block does not
reduce to a pure-Casimir formula.

## Out of scope

- General (non-two-row-rectangular) reps at `z=1` — this is the two-row-rectangular special case only.
- Racah/6j-type multiplicity-resolving data for any `N_μ>1` blocks found — reported, not derived.
- Any change to the already-committed twisted-case cells (1–34).
- Promotion to `Clean/`.
