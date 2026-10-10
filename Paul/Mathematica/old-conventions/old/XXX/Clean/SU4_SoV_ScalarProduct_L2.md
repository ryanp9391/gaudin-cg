# `SU4_SoV_ScalarProduct_L2.wb` — companion reference

The FSoV scalar product for the twisted gl(4), `L=2` chain with rectangular site representations
`[S,S,0,0] ⊗ [P,P,0,0]`, and the computational proof that distinct eigenstates are orthogonal.

- **Derivation:** `SU4_SoV_ScalarProduct_L2.tex` (same directory).
- **Exploratory origin:** `../Experiments/SU4_twisted_Baxter_L2.wb`, Sections 11–12 (including the
  dead ends — worth reading before trying a variant).
- **Status:** every result below was verified in a fresh kernel by running the notebook top to
  bottom. The orthogonality table takes ~48 s.

## 0. Read this first

Two conditions are not cosmetic; the construction fails without them.

1. **The twist must be unimodular**, `z1 z2 z3 z4 = 1`. Otherwise the pairing `q̂^j q̂_i` retains a
   factor `Ztot^{u/h}` from the Hodge-dual normalisation and nothing vanishes.
2. **The contour encircles all poles** of `E_1`, both ladders, `S+P+2` of them. Taking only the
   `S+P` poles of the bare measure `μ` leaves a nonzero remainder.

Also: the `z_j` must be distinct (Vandermonde in the duals), and `θ1 ≠ θ2` generically.

## 1. Conventions

| symbol | meaning |
|---|---|
| `S`, `P` | site representations `[S,S,0,0]`, `[P,P,0,0]`; `S_1 = S`, `S_2 = P` |
| `θ1`, `θ2` | inhomogeneities |
| `h` | shift unit; `f^{[±2n]}(u) = f(u ± n h)` |
| `z[j]`, `j=1..4` | diagonal twist eigenvalues, unimodular and distinct |
| `μ = {μ_1..μ_4}` | gl(4) weight of the state; `Σ μ_k = 2(S+P)` |
| `M_k = (S+P) - μ_k` | degree of the Baxter polynomial `q_k` |
| `numU` | test substitution: `θ1→1/3, θ2→1/7, h→1, z→{2,3,1/2,1/3}` (unimodular, exact rationals) |

Everything is exact and non-perturbative in the twist — no `eps` expansion anywhere.

## 2. Symbol reference

### Cell 2 — twist, weights, transfer matrices

| symbol | signature | purpose |
|---|---|---|
| `zvars` | list | `{z[1],…,z[4]}` |
| `χ[a]` | `χ[a]` | elementary symmetric polynomial of the twists; memoised |
| `ν[S,P][i][u]` | function | Yangian weight `(u-θ1-h λ1_i)(u-θ2-h λ2_i)`; rectangular reps collapse these to `ν1=ν2` and `ν3=ν4` |
| `ν1[S,P][u]` | function | `ν[S,P][1][u] = (u-θ1-Sh)(u-θ2-Ph)` |
| `qdet[S,P][u]` | function | quantum determinant `χ4 ∏_i ν_i(u-(4-i)h)`; central, state-independent |
| `Ξ[a,n]`, `subl[a,n]` | functions | `Ξ_a = Σ_i n_i z_i ∂χ_a/∂z_i`; `subl[a,n] = [u^{2a-1}]τ_a = -a(θ1+θ2+(a-1)h)χ_a - h Ξ_a` |
| `tc[1..5]`, `tcvars` | symbols | the only free transfer-matrix coefficients |
| `τ[a][S,P,μ][u]` | function | transfer matrix, `a=1..4`; degree `2a`, leading `χ_a`, subleading `subl`, free parts `tc[1]` in `τ1`, `tc[2..4]` in `τ2`, `tc[5]` in `τ3` |
| `τ3quad[S,P,μ]` | value | the `u`-coefficient of `τ3`'s quadratic, fixed by matching `subl[3,μ]`; memoised |

`τ3` is written in factorised form `ν1(u-2h) ν4(u) × (quadratic)` — the trivial factor used in the
derivation.

### Cell 3 — Baxter equation and solver

| symbol | signature | purpose |
|---|---|---|
| `Mdeg[S,P,μ]` | list | `M_k = (S+P) - μ_k` |
| `qmonic[k,M][x]` | function | monic ansatz `x^M + Σ qc[k,j] x^j` |
| `qpoly[clist,u]` | value | polynomial from an ascending coefficient list |
| `BaxRes[S,P,μ][k][qf][x]` | expression | branch-`k` Baxter residual, exponential stripped (hence the explicit `z[k]` powers) |
| `TwistedSolve[S,P,μ,sub]` | list of associations | exact joint solve for `tc[1..5]` and the four monic `q_k`; one association per state of weight `μ`, keys `"M"`, `"tc"`, `"q"` |
| `numU` | rules | the unimodular test substitution |

The number of solutions equals the weight multiplicity in `[S,S,0,0] ⊗ [P,P,0,0]`.

### Cell 5 — Casoratians, vacuum factor, Q-functions, operators

| symbol | signature | purpose |
|---|---|---|
| `tcOf[r]` | rules | `tc[i] → r["tc"][[i]]` for a solution record |
| `Wronsk[rec,II][x]` | value | `|II|×|II|` Casoratian `det[z_{i_a}^{k-b} q_{i_a}(x-(b-1)h)]` |
| `Vand[II]` | value | Vandermonde `∏_{a<b}(z_{i_a} - z_{i_b})`; the leading coefficient of `Wronsk` |
| `qWronsk[rec,II][x]` | value | `Wronsk/Vand`, monic of degree `Σ_{i∈II} M_i` |
| `AA[S,P][x]` | value | vacuum factor `∏_{k=1..S}(x-θ1-kh) ∏_{k=1..P}(x-θ2-kh)`, degree `S+P` |
| `QQfun[r,sub][i][y]` | value | `q̂_i(y) = z_i^{y/h} q_i(y)` |
| `qDual[S,P,r,sub][j][v]` | value | monic dual polynomial `(W_ȷ̄/V_ȷ̄/A)(v+h)`, degree `μ_j` |
| `QDual[S,P,r,sub][j][v]` | value | `(-1)^j ζ_ȷ̄^{v/h} qDual[j][v]` |
| `qhat[r,sub][i][y]` | value | `= QQfun`; the family annihilated by `O` (after the `[2]` shift) |
| `qdualhat[S,P,r,sub][j][y]` | value | `QDual[j][y-h]`; the family annihilated by `O†` |
| `OpO[S,P,μ,r][G][x]` | expression | the fourth-order operator `O`, with `T_{a,1} = τ_a` |
| `OpOd[S,P,μ,r][G][x]` | expression | the dual operator `O†`, `ν3 → ν4` |

**The `h` shift in `qdualhat` matters.** `QDual` is centred (it is the object appearing at argument
`u-h` in the QQ-relation `F(u) Q^4(u-h) = det[Q_1(u),Q_2(u-h),Q_3(u-2h)]`); the function solving the
dual Baxter equation is the *uncentred* one, `QDual[j](u-h)`. Using `QDual` directly satisfies
`O†q = 0` only for `j=4`, and then only accidentally (its polynomial part is constant).

### Cell 6 — check of the two defining relations

`relCheck[S,P,μ,sub]` evaluates `O q̂_k^{[2]}` and `O† q̂^k` at a generic point, and compares
`deg q̂^k` with `μ_k`. All four branches vanish, at `(1,1)` and `(2,2)`.

### Cell 8 — measure and scalar product

| symbol | signature | purpose |
|---|---|---|
| `mPoles[S_α,θ_α,sub]` | list | the `S_α` poles of the bare measure, `θ_α + n h`, `n = 0..S_α-1` |
| `mMeas[S,P,sub][x]` | value | `μ(u) = Γ[ν1^{[2]}]/Γ[ν4^{[2]}]`, telescoped to `1/∏_α ∏_{n=0}^{S_α-1}(u-θ_α-nh)` |
| `E1meas[S,P,sub][y]` | value | the common effective measure `E_1 = μ/ν1` |
| `ladE1[S,P,sub]` | list | poles of `E_1`: `θ_α + n h`, `n = 0..S_α`, i.e. `S+P+2` of them |
| `perFac[α,sub][x]` | value | `1` for `α=1`, the `h`-periodic `Exp[2πi x/h]` for `α=2` |
| `brkE1[S,P,sub][f,α,x]` | value | the contour integral: sum of residues of `E_1 · perFac · f` over the whole ladder |
| `spMatFSoV[S,P,rA,rB,sub]` | 8×8 matrix | `M^{AB}`; rows `(α,i,j)`, columns `(r,β)` |

Row and column index sets:

```
rows: (α, i, j)   α ∈ {1,2},  i ∈ {1,2},  j ∈ {3,4}
cols: (r, β)      r ∈ {1, 21, 22, 3},  β ∈ {0,1}

r = 1  :  q̂^j_A(u)    q̂^B_i(u)
r = 21 :  q̂^j_A(u)    q̂^B_i(u-h)
r = 22 :  q̂^j_A(u+h)  q̂^B_i(u)
r = 3  :  q̂^j_A(u+h)  q̂^B_i(u-h)
```

The four `r` are the four shift combinations `(0,0), (0,-h), (+h,0), (+h,-h)` of the left dual
against the right Q; `β` is the insertion `u^β`.

### Cell 9 — orthogonality

`orthTest[S,P,μ_A,k_A,μ_B,k_B,sub]` solves both states (taking solution `k` of each weight sector)
and returns whether `det M^{AB}` is exactly zero.

| case | result |
|---|---|
| `(1,1)` `{2,2,0,0}` vs `{2,1,1,0}`, both orderings | `det = 0` |
| `(1,1)` `{2,1,1,0}` #1 vs #2, both orderings (**same weight**) | `det = 0` |
| `(1,1)` each state against itself | `det ≠ 0` |
| `(2,2)` `{4,4,0,0}` vs `{4,3,1,0}`, both orderings | `det = 0` |
| `(2,2)` `{4,3,1,0}` #1 vs #2, both orderings (**same weight**) | `det = 0` |
| `(2,2)` each state against itself (3 states) | `det ≠ 0` |

8/8 cross determinants vanish, 5/5 norms do not. The same-weight rows are the sharp test: those
states share every `M_k` and differ only in their integrals of motion, so nothing but the SoV
structure can force the determinant to vanish.

The test is `Simplify[...] === 0` on exact algebraic numbers, **not** numerical. `N[]` cannot
certify these determinants: they are identically zero, so a numerical evaluation just exhausts
`$MaxExtraPrecision` and warns.

## 3. Things that do not work (checked, in the Experiments notebook)

- Measure poles at `θ_α ± (2l-1)h/2` with the paper's 8×8 determinant: vanishes at `(1,1)` by
  coincidence, fails at `(2,1)` and `(2,2)`.
- One functional per pole with `β = 1..S+P`: the matrix is degenerate, norm included.
- Fitting `(a,b)` in a pole set `θ-a h, θ, θ+b h` by demanding the determinant vanish: the 2D
  root-find escapes to `a ≈ 10³`, where the determinant underflows trivially.
- Using `QDual` without the `-h` shift as the dual: `O†q ≠ 0` for `j = 1,2,3`.

## 4. Caveat on empty checks

If `M_i = 0` the polynomial `q_i` is constant and the corresponding integrand can lose all its
poles, so that `(i,j)` pair vanishes trivially rather than as a test. Prefer states with `M_1 > 0`
when adding new checks.
