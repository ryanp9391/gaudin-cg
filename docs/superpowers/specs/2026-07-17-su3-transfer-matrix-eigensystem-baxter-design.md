# su(3) transfer-matrix eigensystem + Baxter Q-functions — design

**Date:** 2026-07-17
**Subproject:** XXX (gl(3))
**Notebook:** `Paul/Mathematica/XXX/Experiments/su3_V1.wb` (extend in place)
**Status:** design approved 2026-07-17

## Goal

Build the gl(3) analog of the su(2) Clean notebook's **Section B** (Baxter/TQ solver): an
eigenvector builder and eigenvalues for the diagonal-twist transfer matrices of the two-site
gl(3) XXX chain (each site a symmetric rep `[S,0,0]`), plus the three Baxter Q-functions.

Eigenstates are labeled by a triple `(H1,H2,n)`:
- `H1 = E11 − E22`, `H2 = E22 − E33` (Cartan weights of the eigenspace),
- `n` indexes solutions within a fixed `(H1,H2)` sector.

The 3rd-order Baxter equation has three solutions `Q_k(u) = z_k^{u/h} q_k(u)`, `k=1,2,3`, with
`z_k` the companion-twist eigenvalues and `q_k` a polynomial of degree `M_k`. Determining how the
`M_k` relate to `(H1,H2,S1,S2)` is an explicit deliverable.

## Fixed decisions (context)

- **Twist frame = diagonal (`tg`).** The `(H1,H2,n)` labeling requires a frame in which the
  transfer matrices commute with the Cartan generators. The companion twist `TG`/`tG` mixes weight
  sectors, so the eigenvector builder uses diagonal-twist `tg[S1,S2]`. (Matches su(2), which built
  `Psi` from `t2diag`.) Companion-twist SoV covectors are a later step (the CG-overlap program),
  out of scope here.
- **Eigenvector method = direct simultaneous diagonalization** (user-chosen) of the commuting
  transfer-matrix + Cartan family — the generalization of su(2)'s `TauEigensystem`. Not a
  constructive Bethe/SoV builder; Q-functions come separately from the Baxter solver.
- **Numeric spectrum, symbolic-in-`u` eigenvalues.** Within a multiplicity-`>1` sector the
  eigenvalues are roots of a characteristic polynomial not generally solvable in radicals, so
  diagonalize numerically (numeric `z[i]`, `θ1,θ2`, `h`) while keeping `τ1(u)`, `τ2(u)`, `q_k(u)`
  symbolic in `u` with numeric coefficients. (Matches su(2).)
- **General reps, validated on a sweep** `{(1,1),(2,1),(1,2),(2,2),(3,1),(1,3)}`.
- **Extend `su3_V1.wb` in place** (Experiments/), mirroring Section B naming. No Clean promotion in
  this task.

## Existing state (`su3_V1.wb`, the Section A analog — already present)

gl(3) generators via GT patterns (`ps`,`dim`,`J0`,`Jp`,`Jm`,`JJ`); single-/two-site generator
embeddings (`Ee`,`id`,`idd`,`EE`); central charge `C1` (checked `= S1+S2`); Lax `L[S]` for
symmetric reps `[S,0,0]`; companion twist `G` with eigenvalues `z[1],z[2],z[3]` (and characters
`χ1,χ2,χ3`); monodromy in three frames `T` (untwisted), `Tg` (diagonal), `TG` (companion);
Yangian weights `ν[S1,S2][i][u]` and HWS (checked); quantum minors `qmin`/`qming`/`qminG`;
transfer matrices `tg[S1,S2][a,1][u]` (diagonal) and `tG` (companion); quantum determinant value
`qdetT` (checked `= tg[·][3,1]`); asymptotic `tasympt[1,1]`; a **sketched** 3rd-order Baxter cell
`Bax` (references placeholder names `τg`/`ν1`/`qdet` that must be reconciled to `tg`/`ν`/`qdetT`).

## Component 1 — Eigenvector builder & spectrum

**Commuting family** (all pairwise-commuting, verified once):
- coefficient matrices of `τ1(u) = tg[S1,S2][1,1](u)` (degree 2 in `u`),
- coefficient matrices of `τ2(u) = tg[S1,S2][2,1](u)` (degree 4 in `u`),
- Cartan operators `H1op[S1,S2] = EE[S1,S2][1,1] − EE[·][2,2]`,
  `H2op[S1,S2] = EE[·][2,2] − EE[·][3,3]`.

**`tgCoeffs[S1,S2]`** — memoized; entrywise `CoefficientList` (in a dummy variable, not point
sampling) of `τ1` and `τ2`, returning the coefficient-matrix lists `{C0,C1,C2}` and
`{D0,…,D4}`. Sanity check: leading coefficients are the expected `χ_a · 𝟙` multiples (`C2 = χ1·𝟙`,
`D4 = χ2·𝟙`).

**`wtMult[S1,S2][H1,H2]`** — weight-multiplicity predictor: with `(n1,n2,n3)` recovered from
`(H1,H2)` and `n1+n2+n3 = S1+S2`, count splittings `(a1,a2,a3)+(b1,b2,b3)=(n1,n2,n3)` with
`a_i,b_i ≥ 0`, `Σa=S1`, `Σb=S2` (each symmetric-rep weight has multiplicity 1). Pure combinatorics.
This is the exact per-sector count the diagonalizer must reproduce.

**`TauEigensystem3[S1,S2]`** — memoized. Diagonalizes one **generic real combination** of the
non-trivial commuting members (su(2)'s `C0 + r·C1` trick, extended; retry an `r`-vector list until
the combination's eigenvectors are full rank `d = dim[{S1,0,0}]·dim[{S2,0,0}]` with `d` distinct
eigenvalues — guards accidental degeneracy). For each eigenvector `v`:
- extract every conserved quantity (all `τ1`/`τ2` coefficient matrices, `H1op`, `H2op`) as a
  Rayleigh quotient `vᵀ M v / vᵀ v`;
- hard-`Abort` if any residual `‖M v − (quotient) v‖ > 1e-15` (this residual check **is** the
  genuineness test — no post-hoc spurious-state filter);
- read `H1,H2` directly from the `H1op,H2op` quotients, hard-asserted within `1e-10` of integers.

States within each `(H1,H2)` sector are sorted ascending by `(Re,Im)` of a fixed conserved
quantity to assign `n`; each sector's count is hard-asserted `== wtMult[S1,S2][H1,H2]`.

**Returns** an `Association` keyed `{H1,H2,n} → <|"tau1coef"->…, "tau2coef"->…, "Psi"->v|>`
(`Psi` normalized to a fixed component). Accessors:
- `Psi[S1,S2][H1,H2,n]` — the eigenvector,
- `tau1[S1,S2][H1,H2,n][u]`, `tau2[S1,S2][H1,H2,n][u]` — eigenvalue polynomials (symbolic in `u`),
with bounds/spurious guards analogous to su(2)'s `Psi`/`τ`.

## Component 2 — Baxter solver & degree relation

**Reconcile `Bax`** with Section A's real names, yielding the 3rd-order TQ operator
```
Q(u+h) − A1(u) Q(u) + A2(u) Q(u−h) − A3(u) Q(u−2h) = 0
A1 = τ1/ν1[u],  A2 = τ2/(ν1[u] ν1[u−h]),  A3 = qdetT /(ν1[u] ν1[u−h] ν1[u−2h])
```
with `ν1[S1,S2][u] = ν[S1,S2][1][u] = (u−θ1−h S1)(u−θ2−h S2)` and `qdetT` the twisted quantum
determinant (which already carries the `χ3` factor — the sketched `Bax` cell's explicit `χ3` must
not be double-counted). Reconcile the sketch carefully: confirm the overall sign/normalization of
the `A3` term (and whether the intended `qdet` is the twisted `qdetT` or the untwisted
`ν1[u−2h] ν2[u−h] ν3[u]`) so the operator annihilates the known Q's. Once `τ1,τ2` are known from Component 1, `A1,A2,A3` are explicit functions of `u`.

**Three exponents** `z_1,z_2,z_3` = roots of `z³ − χ1 z² + χ2 z − χ3 = 0` (the companion-twist /
`G` eigenvalues), assigned numerically.

**`QSolve3[S1,S2][H1,H2,n]`** — for each `k`, substitute `Q_k = z_k^{u/h} q_k(u)`, obtain a linear
recurrence for the coefficients of the monic degree-`M_k` polynomial `q_k`, and solve as a
`NullSpace` problem at the known eigenvalue (generalizing su(2)'s `QSolve`; `Chop` the
precision-tracked matrix before `NullSpace`). Returns the three coefficient lists `{q_1,q_2,q_3}`.

**Degree prediction.** Master formula from asymptotic matching (writing
`Q_k = z_k^{u/h} u^{M_k}(1 + …)` and balancing the `1/u` term of the TQ equation):
```
M_k = ( α1 − α2/z_k + α3/z_k² ) / ( h ( z_k − χ2/z_k + 2 χ3/z_k² ) )
```
where `α1,α2,α3` are the `u¹`-coefficients of `A1,A2,A3` (equivalently of `τ1,τ2,qdet`) on the
state — numeric functions of `(n1,n2,n3)`. Steps:
1. compute the `α_i` from the eigenvalue asymptotics (generalize `tasympt` to `a=2` and the qdet);
2. predict `M_k`, confirm the solved `q_k` actually has that degree;
3. tabulate `M_k` vs `(H1,H2,S1,S2)` across the sweep and extract the **closed-form** relation
   (the "figure it out" deliverable).

**Accessor** `Q[S1,S2][H1,H2,n][k][u] = z_k^{u/h} qpoly[q_k, u]`.

## Component 3 — Verification suite

Across the sweep `{(1,1),(2,1),(1,2),(2,2),(3,1),(1,3)}`, per state:
1. **TQ residual** `≈ 0` (numeric, symbolic-in-`u`) for each `Q_k`.
2. **Wronskian / QQ check** — the 3×3 Casoratian (discrete Wronskian) of `{Q_1,Q_2,Q_3}`
   reproduces the quantum determinant `qdetT` (up to the standard normalization).
3. **Eigenvalue consistency** — `τ1,τ2` reconstructed from the `{Q_k}` (via the symmetric-function
   / fusion relations) match the `TauEigensystem3` eigenvalues.
4. **Sector counts** — states per `(H1,H2)` sector `== wtMult`.
5. **Degree formula** — solved `deg q_k == M_k(H1,H2,S1,S2)` from the extracted closed form.

Each check follows the notebook's hard-assert convention (a real failure aborts).

## Naming summary (mirrors su(2) Section B, gl(3)-generalized)

| su(2) Section B | gl(3) here |
|---|---|
| `t2Coeffs` | `tgCoeffs` (both `τ1`,`τ2`) |
| `dPred[M,λ1,λ2]` | `wtMult[S1,S2][H1,H2]` |
| `TauEigensystem` | `TauEigensystem3` |
| `Psi`, `τ` | `Psi`, `tau1`, `tau2` |
| `BaxterMatrix`/`QSolve` | `QSolve3` (per exponent `z_k`) |
| `Q1` | `Q[S1,S2][H1,H2,n][k]` |
| (deg fixed = `M`) | `M_k` degree formula (derived + verified) |

## Out of scope

- Companion-twist SoV covectors / FSoV scalar product / CG-overlap dictionary (the later program).
- Clean promotion and companion `.md`/`.tex` reference docs.
- Lengths other than `L=2`; non-symmetric reps.
- `paper/xxxCG.tex` write-up.

## Open questions

None. Two minor choices confirmed with the user: `n` ordered by `(Re,Im)` of the eigenvalue within
each sector; sweep set `{(1,1),(2,1),(1,2),(2,2),(3,1),(1,3)}`.
