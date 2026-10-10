# `SU3_310_L2_Untwisting.wb` — companion reference

gl(3) XXX chain, **L = 2**, both sites in **[3,1,0]** (the 15). It repeats `SU3_L2_Untwisting.wb`
(adjoint [2,1,0]): first it re-establishes the twisted identity
`⟨L_A|R_B⟩/(⟨L_A|0⟩⟨0|R_B⟩) = C_vac det N^{AB}/(c^ref_A Ψ^ref_B)`, then it takes `z_j = e^{εφ_j}`, `Σφ_j = 0`, to zero twist.

- **Adjoint parent:** `SU3_L2_Untwisting.wb` / `.md` / `.tex`.
- **Status.** Sections 1–11 were run in order in one live kernel (K2), at 200 digits. There has been **no cold run**.
  Section 12 is incomplete (§5). Stopped at the user's request on 2026-09-16.
- Symbols carry the suffix `F`, so nothing collides with the adjoint notebooks (suffix `U`).

## 0. What changes relative to the adjoint

| item | adjoint [2,1,0] | [3,1,0] |
|---|---|---|
| states / hw states | 64 / 6 | 225 / 10 |
| `ν_i(u)` | `∏_s(u−θ_s−λ_i)` | same, with λ=(3,1,0); `T_3 = ν1(u−2)ν2(u−1)ν3(u)` checked exactly |
| twisted degrees | `4−μ_i`, `μ_j` | `6−μ_i`, `μ_j` |
| untwisted degrees | `(4−λ_i)+(i−1)`, `λ_j+3−j` | `(6−λ_i)+(i−1)`, `λ_j+3−j` |
| **measure poles** | `θ_s+{0,1,2}` (= poles of 1/(ν1ν2ν3)) | **`θ_s+{0,1,2,3}`** (full ladder; 1/(ν1ν2ν3) fails) |
| H column basis | `(c1+c3, c5−c6+5c3, c2−c3)` | `(c1+c3, c5−c6+7c3, c2−c3)`; the constant is `2(λ1−λ3+1)−1` of the site rep |
| O column basis | `(c4, c6+c3, c3)` | same |
| `E_λ` | `−(λ1−λ2+1)(λ2−λ3+1)/[(λ1−λ3+1)!(λ1−λ2)!(λ2−λ3)!]` | **same function of λ** |

Everything else carries over unchanged: the row triple `{(1,3),(1,2),(2,3)}`, the columns
`u^β(−Q^jQ_i, Q^jQ_i^{[−2]}, Q^{j[2]}Q_i)`, the branch-I product reference states (`a13=3/5`, `a23=4/7`), and the
L=1 forms of `Ψ^ref` and `c^ref`.

**States** (`m = 1..10`, `nameF`): 60 (6,2,0) = vacuum, 42 (5,3,0), 21 (6,1,1), 24₁ and 24₂ (5,2,1),
15' (4,4,0), 15̄₁ and 15̄₂ (4,3,1), 6 (4,2,2), 3̄ (3,3,2). Doublets are ordered by `N[[u⁰]T₂]`. At P1 their
spectra lie in ℚ(√100921) (24) and ℚ(√151) (15̄), so the untwisted q's are handled exactly with radicals.

## 1. Cell map

| cell | content |
|---|---|
| 2 | GT generators; `buildChainF[θ]` gets the quantum-minor coefficients by exact interpolation (nodes 1..2|I|+1); `TcoF`, `nu*F`, `qdetF`. Checks: interpolation, commutativity, `T3 = qdet`, gl(3) invariance at z=1 |
| 3 | `BaxDirF`, `BaxDualF`, `polyKernelF` (exact, degree ≤ 14, RootReduce), `hwBasisF`, `buildUntwF`, `untwF`, `nameF`, `psiRefLocF`, `cRefLocF` |
| 4 | `twF[ε,φ]` (memoized): all 225 twisted states, per-sector eigenvectors, left eigenvectors = rows of the blockwise inverse, Q's, `braR = ⟨0|R⟩`, `Lket = ⟨L|0⟩`. About 2 s per call |
| 5 | Measure test (`setMeasureF["nu"/"ladder"]`), `rowF`, `detNF`, `psiRefF`, `cRefF`, orthogonality and the 225-state match. **Leaves the measure set to "ladder"** |
| 6 | Divergence powers, `kPredF`, eigenvector limits |
| 7 | `qSideF[ε,φ]`, ε-powers of every factor |
| 8 | `limDataF[φ]`, `untwFactorsF`, `nF`, `ndF` |
| 9 | `DhatF`, `D0F`, `kappaF`, `psi0NormF`, `matUntwF`, rank table |
| 10 | markdown: findings; H/O derivation |
| 11 | `splitT2F`, `wVecF`, `wCheckF`, `zeroTwistF`, `colChangeF`, `untwNF` (250-digit q's), `regQF`, `matRegF`, `entryPowF` |
| 12 | `splitResF`, `untwHF`, `rowPowPredF`, `EformulaF`, `EnumF`, `relationTabF` |
| 13 | P2 check (**rebuilds the chain at θ=(1/5,46/55)**), `checkP2F` |
| 14 | markdown: the precise relation |
| 15 | `SdiagF` (**fails**, §5), `buildAtF[t]`, `bracketExactF`, `BexactF`, `BratioF`, `compIdxF`, `dataAtF` (**rebuilds the chain at θ=(t,0)**) |
| 16 | markdown: status |

## 2. Twisted identity (cells 4–5; ε = 1e-3, φ = (1, 1/3, −4/3))

- Eigen residual 2e-193; Baxter residual 2e-177 with degrees `6−μ_i`, `μ_j`.
- Multiplet counts are 60, 42, 21, 24, 24, 15, 15, 15, 6, 3. Each multiplet has exactly one hw-limit state.
- **Orthogonality** over the 14 states of weight (4,2,2), max|cross det|/min|diag det|:
  - 1/(ν1ν2ν3), poles θ+{0,1,3}: **1e+34 (fails)**;
  - ladder θ+{0,1,2,3}: **1e-126**.
- **Match:** `det N/(c^ref Ψ^ref)·⟨L|0⟩⟨0|R⟩` is the same for all 225 states, relative spread 1.5e-151.

## 3. Untwisting (cells 6–12, P1 θ=(1/3,−8/57); limits from ε = 1e-6, 1e-7, 1e-8; three φ)

- **Divergence.** `⟨0|R_A⟩ ~ ⟨L_A|0⟩ ~ ε^k`, with `k = ht(λ_vac−λ_A) = 2d1+d2 = 0,1,1,2,2,2,3,3,4,5`. The eigenvectors
  converge to the untwisted hw vectors at O(ε).
- **ε-powers.** `C_vac ~ ε^{−12}`, `c^ref ~ Ψ^ref ~ ε^2`, `det N ~ ε^{16−2k}`, `S_vec ~ ε^{−2k}`. The identity holds at every ε to ≲ 1e-116.
- **Twist factors** (residuals ≤ 1e-16):
  - `Ψ^ref → (z1−z2)²Ψ^ref_untw/n²` and `c^ref → (z2−z3)²c^ref_untw/n'²`, with `n = λ1−λ2+1`, `n' = λ2−λ3+1`;
  - `det N_reg → V_λA(z) D̂⁰_A`, with `D̂⁰` φ-independent to 1e-18;
  - `C_vac → κ/V_λvac(z)`, with `κ = (cΨ)^untw_vac/D̂⁰_vac` to 6e-19.
- **Untwisted 6×6.** Rank 3,3,4,4,4,4,4,4,4,5, against vanishing orders `2(λ1−λ3)` = 12,10,10,8,8,8,6,6,4,2.
  It is never full rank, so the block split is needed for every state.
- **Column identity.**
  - `Σ_c w_c⟨c⟩ = 0` with `w = (Δc1, b1, a2+b2, Δs1, a1, a2)`, where
    `ΔT2 = ∏_s(u−θ_s−λ1−1)(a1u+b1) + ∏_s(u−θ_s−λ3)(a2u+b2)`. Checked on 348 twisted cross rows to 0 at 162 digits.
  - Zero-twist relations `Δs1 = Δs2 = 0`, `a1+a2 = 0`, `b1+b2 = Δc1 + 8a1` hold exactly for all 90 untwisted pairs.
    This gives H.
- **Block powers.** The bare-row profile in (H1..3, O1..3) is `(n+2,n+2,n+2,n,n,n+1)`, as for the adjoint. The O1 power of
  bare row (i,j) is `max(0, 6−d_ij)`, where `d_ij = deg q_i + deg q^j` (untwisted); this matches all 10 states.
  The adjoint had `4−d_ij`.
- **Laplace split.** `det N_reg = −det[A1]_O det[A2]_H (1+O(ε))`. The residual is 1e-6.1 to 1e-6.7 at ε=1e-7 and one decade smaller at 1e-8.
- **Periodic block.** `det[A2]_H → det[A2]^untw_H` at exactly O(ε): residual 1e-4.7 to 1e-7.4 at ε=1e-7, one decade smaller at 1e-8.
- **`E_λ`**, formula vs numeric, relative agreement ~1e-18 for all 10:

  | state | 60 | 42 | 21 | 24 ×2 | 15' | 15̄ ×2 | 6 | 3̄ |
  |---|---|---|---|---|---|---|---|---|
  | `E_λ` | −1/16128 | −1/720 | −1/14400 | −1/90 | −1/576 | −1/8 | −1/4 | −1 |

- **Boxed result.** It holds in the same form as for the adjoint, matching the physical limit to ~1e-19 for all 10 states:
  ```
  lim V_λvac(z)/V_λA(z) · S(A,A) = B_λA/B_λvac,   B_λ = E_λ det𝒩/(c^ref Ψ^ref)
  ```
  𝒩 is the untwisted periodic 3×3 matrix over the ladder poles. Its rows are (1,3), (1,2), (2,3), and its columns are
  `⟨q^j(u+1)q_i − q^jq_i⟩`, `⟨u q^j q_i(u−1) − (u−7) q^j(u+1) q_i⟩`, `⟨q^j q_i(u−1) − q^j(u+1) q_i⟩`.
- **Physical limit at P1** (`psi0NormF`): 60 → 1, 42 → 5776/1363, 21 → 361/35, 24₂ → −95.02…
  The doublets are conjugate pairs, e.g. 24: `(−4153200896027459443 ± 33203096237782613√100921)/154711843570771200`.

## 4. Independent check (cell 13; P2 θ=(1/5,46/55), φ=(3/2,1/2,−2), ε = 1e-7, 1e-8)

- Split residual 1e-7.0 to 1e-7.5, and `det[A2]_H` against the untwisted minor 1e-6.6 to 1e-8, for all 10 states.
- `E_λ` formula against the two-point extrapolation: 1e-12.4 to 1e-13.4 for all 10 states.
- The exact identity `S_Q/S_vec − 1` is below 1e-80.
- **Evidence level.** `E_λ` was taken from the adjoint notebook without refitting. It is **independently validated** on 8 new weights at two chain points,
  but it is **not derived**.

## 5. Incomplete / open

1. **Cell 15 (bare eigenvectors).**
   - The exact explicit formula at θ=(t,0) is ω=e^{2πit}-independent and even in t, and at t=9/19 it matches the physical limit to 2e-19.
   - **`SdiagF` has no solution.** No diagonal S with `E_ji = S⁻¹E_ijᵀS` exists in this GT basis, even for `E21` alone. So the
     adjoint identification `l_0 = (S⊗S)r_0`, `⟨l_0|r_0⟩ = ‖r_0‖²_S`, cannot be carried over as is. The reason was not
     diagnosed.
2. **Not done:** fitting the t-dependence of `B_A/B_vac` (for the doublets, fit `B₁+B₂` and `B₁B₂`), normalizing the bare
   eigenvectors (a canonical component exists only for 60, 42, 21, 15', 3̄), and the separate right/left limits (adjoint cell 19).
3. **Derive `E_λ`.** It is now known to be the same function of λ for two site representations.
4. **Cold top-to-bottom run.**

## 6. Pitfalls hit

- **Wolfbook output truncation.** Replies stop at the first `{`; read results through `showF` (braces replaced).
- **`0^0` in the Vandermonde matrix** when interpolating at node 0; use nodes 1..d+1.
- **`N[...]` of radical expressions** hits `$MaxExtraPrecision`. Evaluate the q coefficients numerically first (`untwNF`).
- **Module initializers use global symbols.** `Module[{x, d = f /. u -> x}, …]` substitutes the global `x`.
- **Cells 13 and 15 overwrite global chain state** (`thetaF`, `polesF`, `untwF`). The kernel is currently at θ=(2/7,0).
