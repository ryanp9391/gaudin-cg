# gl(4) L=2 XXX chain — twisted Baxter solver — design

**Date:** 2026-07-19
**Subproject:** XXX (gl(4))
**Notebook:** `Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb` (currently empty; build the full
pipeline here by adapting each `su3_V2_general.wb` cell to gl(4))
**Status:** design approved 2026-07-19

## Goal

Port the two-site gl(3) XXX construction of `su3_V2_general.wb` (Section A — reps/Lax/twist/transfer
matrices — and Section B — eigensystem + Baxter Q-functions) to **gl(4)**: two independent generic
Young-diagram reps `λ = {λ1,λ2,λ3,λ4}` (`λ1 ≥ λ2 ≥ λ3 ≥ λ4 ≥ 0`, integers), one per site,
inhomogeneities `θ1, θ2`, and **four distinct twist eigenvalues**. Build the eigensystem and four
Baxter Q-functions and pass the same hard-assert TQ / Casoratian verification convention used in the
gl(3) work.

This is the **twisted** version — the deliberate, verifiable stepping stone toward the untwisted
(`z=1`) physical chain. Distinct twist eigenvalues fully break the su(4) symmetry, so the diagonal-
twist transfer matrix is non-degenerate, every eigenstate carries a clean `(H1,H2,H3,n)` label, and
the four Baxter functions `Q_k = z_k^{u/h} q_k(u)` are linearly independent because their exponential
prefactors differ. The `z→1` untwisted limit (where su(4) symmetry is restored, the spectrum
degenerates into irrep multiplets, and the `Q_k` ansatz collapses) is an explicit **follow-on**, out
of scope for this spec.

## Fixed decisions (confirmed with the user)

- **Home = `Baxter_L2_XXX_SU4.wb`.** The full pipeline (Section A reps/Lax/twist/transfer *and*
  Section B eigensystem/Q/verification) lives in this one file — not split into a bare solver. It is
  populated by adapting the `su3_V2_general.wb` cells to gl(4), not re-derived from scratch.
- **Interface = λ-vectors everywhere**, exactly as in gl(3): every function takes full 4-vectors
  `λ1, λ2` (one per site). No scalar-`S` wrappers.
- **Twist = four distinct unit-modulus eigenvalues**, matching the gl(3) convention:
  `z[k] = N[Exp[I Zeta[2k+1]], 30]` for `k = 1,2,3,4` (i.e. `Zeta[3], Zeta[5], Zeta[7], Zeta[9]`).
  `θ1 = 1/3`, `θ2 = 1/7`, `h = 1` carry over unchanged. `χ1..χ4` are the elementary symmetric
  functions of `z[1..4]`.
- **Degree formula = verify the gl(3) closed form generalizes.** Assume `M_k = Λ_1 − n_k` (with
  `Λ = λ1 + λ2` the pseudovacuum highest weight and `Λ_1 = (λ1+λ2)[[1]]`) carries to gl(4); hard-
  assert `deg q_k == M_k` over the whole sweep. Only fall back to a full candidate-fit re-derivation
  if the assumed form fails somewhere.
- **Verification sweep = small set + adjoint⊗fund** (see below), including reps with genuine interior
  weight-multiplicity `> 1` to exercise the generalized `wtMult` and the degenerate
  `TauEigensystem4`.

## Existing state

`su3_V2_general.wb` — a validated 40-cell gl(3) pipeline (full-sweep verified: 102 states,
worstTQ = worstCasoratian = 0). `Baxter_L2_XXX_SU4.wb` — empty. The GT engine `ps/dim/J0/Jp/Jm/JJ`
is **N-general** and already produces gl(4) generators for any Young diagram (as used in
`Gaudin/Clean/Hamiltonians_SU4.wb`); only the layers above it are gl(3)-specific.

## Change classification (gl(3) → gl(4))

### (a) Pure copy / re-index (structurally unchanged, indices now run 1..4)

- GT engine `ps/dim/J0/Jp/Jm/JJ` — copied verbatim.
- `Ee[λ][i,j] := JJ[i,j,λ]`, `id[λ]`, `EE[λ1,λ2][i,j]`, `idd`, `C1[λ1,λ2]`.
- Lax `L[λ][i,j][u] := IdentityMatrix[dim[λ]] u KroneckerDelta[i,j] − h JJ[j,i,λ]`, `i,j` in 1..4.
- Monodromies `T/Tg/TG[λ1,λ2]`, HWS, `qmin/qming/qminG` — same structure, 4×4.
- Section B: `tgCoeffs`, `TauEigensystem4` (renamed from `…3`), `Psi/tau1/tau2/tau3`, `BaxOp`,
  `alphaCoeffs`, `Mdeg`, `QSolve4`, `Qfun` — signatures and internals as in gl(3), extended by one
  index/term where noted below.

### (b) Genuine gl(4) math changes

1. **Twist `G` and `χ`.** 4×4 companion matrix
   `G = {{χ1,-χ2,χ3,-χ4},{1,0,0,0},{0,1,0,0},{0,0,1,0}}`; add `z[4]`; `χ1..χ4` = elementary symmetric
   functions of `z[1..4]`.

2. **Yangian weights** `ν[λ1,λ2][i][u] = (u − θ1 − h λ1[[i]]) (u − θ2 − h λ2[[i]])` for **i = 1..4**
   (from `L[λ][i,i][u]|hw⟩ = (u − h λ_i)|hw⟩` on the two-site HWS). `ν1 := ν[·][1]`.

3. **Transfer matrices.** `tg[λ1,λ2][a,1][u] = Σ` over size-`a` principal quantum minors, for
   **a = 1,2,3**: `τ1 = tg[1,1]` (deg 2 in `u`), `τ2 = tg[2,1]` (deg 4), `τ3 = tg[3,1]` (deg 6).
   Leading coefficients `τ_a → χ_a · Id`. Keep `tG` (companion-twist) for convention/future use.

4. **Quantum determinant.**
   `qdetT[λ1,λ2][u] = χ4 · ν[·][1][u−3h] · ν[·][2][u−2h] · ν[·][3][u−h] · ν[·][4][u]`.
   Guard: `qdetT` central value `== tg[λ1,λ2][4,1]` (the top minor) — the ν shift-ordering check on
   non-symmetric reps.

5. **4th-order Baxter/TQ equation.**
   ```
   Bax = f[u+h] − A1 f[u] + A2 f[u−h] − A3 f[u−2h] + A4 f[u−3h]
   A1 = τ1 / ν1[u]
   A2 = τ2 / (ν1[u] ν1[u−h])
   A3 = τ3 / (ν1[u] ν1[u−h] ν1[u−2h])
   A4 = qdetT / (ν1[u] ν1[u−h] ν1[u−2h] ν1[u−3h])
   ```
   `qdetT` already carries `χ4`, so no explicit `χ4` in `A4` (would double-count) — same convention
   as gl(3)'s `A3`.

6. **`tgCoeffs`.** Extract coefficient matrices of `τ1, τ2, τ3` (`PadRight` to 3, 5, 7 respectively).

7. **Cartan operators.** Three gl(4) weight generators `H1 = E11−E22`, `H2 = E22−E33`,
   `H3 = E33−E44` on the two-site space; state sectors labelled `(H1,H2,H3)`.

8. **`wtMult` — honest gl(4) tensor-product weight tally.**
   - `siteWeightMult[λ]` — tally the gl(4) weight `(w1,w2,w3,w4)` of every GT pattern in `ps[{λ}]`,
     returning `weight → multiplicity`.
   - `wtMult[λ1,λ2][H1,H2,H3]` — recover `(n1,n2,n3,n4)` via `H1=n1−n2`, `H2=n2−n3`, `H3=n3−n4`,
     `n1+n2+n3+n4 = N` (total boxes `= Σ λ1_i + Σ λ2_i`); reject non-integer/negative; return
     `Σ_{a+b=(n1,n2,n3,n4)} siteWeightMult[λ1][a] · siteWeightMult[λ2][b]`.

9. **`TauEigensystem4`.** Simultaneously diagonalize the three commuting matrices `τ1, τ2, τ3` per
   `(H1,H2,H3)` sector (distinct twist ⇒ generically 1-dim simultaneous eigenspaces). Extract
   symbolic-in-`u` eigenvalues `Ψ, tau1, tau2, tau3`.

10. **Asymptotics + degree formula.** `alphaCoeffs → {α1,α2,α3,α4}`; asymptotic master formula runs
    over the four twist eigenvalues `z[k]`. **Verify** `M_k = Λ_1 − n_k` (`k=1..4`,
    `Σ_k M_k = 4Λ_1 − N`); hard-assert `deg q_k == M_k`.

11. **`QSolve4` / `Qfun`.** Solve the 4th-order Baxter equation for the four monic polynomials
    `q_k` (degrees `M_k`); `Qfun` returns `Q_k(u) = z[k]^{u/h} q_k(u)`.

## Verification (acceptance)

Full sweep over these rep pairs `(λ1, λ2)`:

| `λ1` | `λ2` | tensor dim | note |
|---|---|---|---|
| `{1,0,0,0}` | `{1,0,0,0}` | 16 | baseline (4⊗4 = 10 ⊕ 6), all mult 1 |
| `{1,0,0,0}` | `{1,1,1,0}` | 16 | 4 ⊗ 4̄ = 15 ⊕ 1 |
| `{1,1,0,0}` | `{1,0,0,0}` | 24 | 6 ⊗ 4 |
| `{1,1,0,0}` | `{1,1,0,0}` | 36 | interior multiplicity |
| `{2,1,1,0}` | `{1,0,0,0}` | 60 | adjoint (zero-weight mult 3) ⊗ fund — **required**, stresses `TauEigensystem4` |

Required hard-asserts (all machine-zero / `< 10^-9`):

- `wtMult` sum-rule: `Σ_{(H1,H2,H3)} wtMult == dim[λ1] · dim[λ2]` for every rep.
- `qdetT == tg[4,1]` shift-ordering guard.
- `τ_a` leading coeff `== χ_a · Id` (a=1,2,3); `{τ1,τ2,τ3}` pairwise-commute.
- `deg q_k == M_k` with `M_k = Λ_1 − n_k`.
- **worstTQ = 0**: every `Q_k` solves the 4th-order TQ equation.
- **worstCasoratian = 0**: the Casoratian reconstruction of `A1..A4` from the four `Q_k` matches
  `τ1, τ2, τ3, qdetT`.

Final acceptance = a fresh-kernel top-to-bottom run (restart → run all cells in order) with every
hard-assert passing in-kernel and no `Abort`.

## Out of scope

- The `z→1` untwisted / physical limit (su(4) irrep multiplets, degenerate spectrum, CG overlaps) —
  separate follow-on spec.
- Companion-twist SoV covectors / the CG-overlap program.
- Promotion to `Clean/`.

## Gotchas (from the gl(3) work)

- Multi-statement Wolfram cells need `;` between statements.
- Memoized definitions (`f[x_]:=f[x]=…`) cache stale results — after fixing a definition, re-run the
  defining cell and everything downstream.
- wolfbook can leave stale/blank `Out[]` boxes after `insertCells`/edits — re-run affected cells, do
  not trust blank output boxes.
- Cell edits/evaluations don't persist until the notebook is explicitly saved.
