---
# Left eigenvectors of the twisted transfer matrix and the <Psi|Psi>/det[Q,Q] ratio — design

## Context

`Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb` Section B's `TauEigensystem[λ1,λ2]`
diagonalizes the companion-twisted transfer matrix `t2[λ1,λ2][1,1][u] = C0 + C1 u + C2 u^2`
(Section A) via a generic combination `C0 + r C1`, producing, per genuine state `{M,n}`, the
**right** eigenvector `Psi[λ1,λ2][M,n]` (normalized so its last component is 1) and the τ
eigenvalue `τ[λ1,λ2][M,n][u]`. Section C then ported the FSoV measure/bracket/determinant
machinery from `FSoV_L2.wb` and built `det[λ1,λ2][Q_A(u)Q_B(u)]`, verifying numerically (via
`GramMatrix`/`G21`) that it vanishes for distinct states and is nonzero on the diagonal for
`(λ1,λ2)=(2,1)`.

This design covers the natural next step: also compute the **left** eigenvectors of the same
twisted transfer matrix (same normalization convention — last component 1), pair each with its
same-state right eigenvector to form `<Psi|Psi>`, and compute the ratio `<Psi|Psi>/det[Q,Q]`
per state. `paper/xxxCG.tex` only covers the single-site (`L=1`) SoV construction so far — there
is no two-site scalar-product formula written up yet to check this ratio against. This
calculation is exploratory: report the values, don't assert against a predicted closed form.

## Scope

- Add a new Section D to `XXX_CG_L2_V2.wb`, after Section C. Section B (`TauEigensystem`,
  `Psi`, `τ`, `QSolve`, `Q1`) and Section C (`mu`/`Nn`/`nn`/`br`/`det`/`GramMatrix`/`states21`/
  `G21`) are **not modified** — Section D only reads from them.
- `TauEigensystemLeft[λ1,λ2]`: an independent routine (not an edit of `TauEigensystem`) that
  diagonalizes `Transpose[C0 + r C1]` for a generic `r` (own `rCandidates`/`SelectFirst`
  selection, duplicated from `TauEigensystem`'s logic rather than importing it, since `r` only
  needs to be *some* value that resolves the eigenspace with a simple spectrum — it need not be
  numerically identical to the `r` `TauEigensystem` happened to pick).
- Matching left eigenvectors to existing `{M,n}` states: for each state already known from
  `TauEigensystem[λ1,λ2]`, compute the target combo-eigenvalue
  `τ[λ1,λ2][M,n][r] - (κ1+κ2) r^2` (valid because `t2(r).v = τ(r) v` for the right eigenvector
  and `C2 = (κ1+κ2) Id` exactly, already verified in Section B), and match it to the closest
  eigenvalue of `Transpose[C0 + r C1]`. Hard-assert if the closest match is not unique or not
  within tolerance.
- Normalize each matched left eigenvector so its last component is 1, same convention and
  same-style hard assert (`lastzero`) as `Psi`.
- `PsiLeft[λ1,λ2][M,n]`: on-demand accessor, same bounds/spurious-guard pattern as `Psi`/`τ`.
- `braket[λ1,λ2][M,n] := PsiLeft[λ1,λ2][M,n] . Psi[λ1,λ2][M,n]` — plain bilinear pairing, no
  complex conjugation (matches this notebook's existing non-Hermitian convention, e.g.
  `TauEigensystem`'s `t0 = (v.C0.v)/(v.v)`).
- `normRatio[λ1,λ2][M,n] := braket[λ1,λ2][M,n] / G21[[i,i]]`, where `i` is `states21`'s index of
  `{M,n}` — reuses the diagonal of the already-computed Gram matrix, no new determinant
  machinery.

Out of scope: sweeping over other `(λ1,λ2)` representations; comparing `normRatio` to a
predicted closed-form formula (nothing to compare against yet — this calculation is how one
would be discovered); touching Section B or Section C.

## Test case

Same as Section C: `(λ1,λ2) = (2,1)`, reusing the already-computed `states21`
(`{{0,0},{1,0},{1,1},{2,0},{2,1},{3,0}}`) and `G21` (6x6 Gram matrix, verified numerically
diagonal).

## Architecture

New cells appended after Section C in `XXX_CG_L2_V2.wb`:

1. **Section D header comment** — states purpose, references this design doc.
2. **`TauEigensystemLeft[λ1,λ2]`** — as described above: own generic-`r` selection, diagonalize
   `Transpose[C0+r C1]`, match each left eigenvector to a `{M,n}` state via the
   `τ[λ1,λ2][M,n][r] - (κ1+κ2) r^2` target eigenvalue, hard-assert uniqueness/tolerance of the
   match, normalize to last-component-1. Returns an `Association` keyed by `{M,n} -> "PsiLeft"`,
   mirroring `TauEigensystem`'s shape. Memoized per `(λ1,λ2)`.
3. **`PsiLeft[λ1,λ2][M,n]`** — accessor, same bounds/spurious pattern as `Psi`.
4. **Biorthogonality smoke test** — `PsiLeft[λ1,λ2][M,n] . Psi[λ1,λ2][M',n']` for all
   `(M,n)≠(M',n')` pairs among `states21`, hard-assert max|value| at the numeric noise floor
   (matching Section C's `10^-15` convention). This is the correctness check on the left-
   eigenvector construction itself (a genuine left/right eigenvector pair of a diagonalizable
   operator with simple spectrum must be biorthogonal across distinct states) — the natural
   companion to Section C's off-diagonal Gram check.
5. **`braket[λ1,λ2][M,n]` and `normRatio[λ1,λ2][M,n]`** — as described above.
6. **Final report cell** — for each of the 6 states of `(2,1)`: `{M, n, braket, G21 diagonal
   entry, normRatio}`, printed/returned as an `Association`/table. Hard-assert `braket` nonzero
   for every state (non-degenerate pairing) before reporting. No comparison to a predicted
   formula.
7. Save the notebook.

## Risks

- The generic-`r` selection for `TauEigensystemLeft` could, in principle, pick a value making
  `Transpose[C0+r C1]` fail to resolve the full `d`-dimensional eigenspace (same failure mode
  `TauEigensystem` already guards against) — same `rCandidates`/`SelectFirst`/hard-assert
  pattern handles this identically.
- The eigenvalue-matching step (`τ[λ1,λ2][M,n][r] - (κ1+κ2) r^2` vs. `Transpose[C0+r C1]`'s
  eigenvalues) relies on the combo eigenvalues being pairwise distinct across the `d` states —
  already guaranteed by `TauEigensystem`'s own construction (`Length[Union[Chop[ee,10^-12]]] ==
  d` for its chosen `r`; `TauEigensystemLeft`'s independent `r` selection re-derives this same
  guarantee for its own `r`, so a match should always be unique to numeric precision, but the
  plan hard-asserts this rather than assuming it).
