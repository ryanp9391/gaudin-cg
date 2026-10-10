# `SU3_Adjoint_Twisted_SoV_L2.wb` — companion reference

Twisted functional SoV scalar product for the gl(3) chain, `L = 2`, sites `[2,1,0] ⊗ [2,1,0]`, and
its zero-twist limit. The **untwisted 3×3 norm** of `SU3_Adjoint_SoV_ScalarProduct_L2.wb` is
identified as the **upper-left minor of the twisted 6×6 norm**, once the six integrals of motion
are split into three that survive at zero twist and three that are `O(φ)`.

- **Write-up:** `SU3_Adjoint_Twisted_SoV_L2.tex`.
- **Untwisted companion:** `SU3_Adjoint_SoV_ScalarProduct_L2.wb` / `.md` / `.tex`.
- **gl(4) precedent:** `untwisted-norm.md`, Part 3.
- **Status:** self-contained; the cold run (fresh kernel) takes about 10 s after the cell-4 build,
  with no kernel messages. Numerics at 150 digits. The untwisted reference is exact.

## 0. Read this first

1. **Unimodular twist.** `z_j = e^{εψ_j}`, `ψ = (1, 1/3, -4/3)`, `Σψ = 0`. Hence
   `T_3 = z1z2z3 · qdet = qdet`, and the Baxter operators below use `T_3 = qdet`.
2. **Same parameter point as the companion.** `θ1 = 1/3, θ2 = -8/57, h = 1`. `h = 1` is hard-coded
   in `ν1, ν2, ν3, qdet`, in the shifts `x-1`, and in `colFuns`. Changing `h` means editing those.
3. **State matching.** A twisted state is picked in the weight sector `μ_U` as the eigenvector whose
   `T_1, T_2` coefficients are closest to the untwisted `U`. All six untwisted states have distinct
   `(T_1, T_2)`, so this is unambiguous for small `ε`.
4. **Raw twisted Q's cannot converge individually.** Their degrees `4-μ_k`, `μ_j` are smaller than
   the untwisted staggered ones. Only the regularised combinations `R_k`, `R^j` tend to the
   untwisted `q`.
5. Wolfbook replies truncate at `{`; read results from stored variables.

## 1. Conventions

| symbol | meaning |
|---|---|
| `T_a(u)` | `Σ_{|I|=a} z_I qmin_{I,I}(u)`, `z_I = ∏_{i∈I} z_i` |
| `T_1` | `χ_1 u² + s_1 u + c_1` |
| `T_2` | `χ_2 u⁴ + s_2 u³ + c_22 u² + c_21 u + c_20` |
| `Q_k` | `z_k^u q_k`, `deg q_k = 4-μ_k`, monic (raw) |
| `Q^j` | `z_j^{-u} q^j`, `deg q^j = μ_j`, monic (raw) |
| `w` | `(Δc_1, Δc_21, Δc_20, Δc_22-Δc_1, Δs_1, Δs_2)` |
| `F_1..F_6` | `Q^j(u²Q_i^{[-2]} - ν1^{[-2]}Q_i)`, `Q^j u Q_i^{[-2]}`, `Q^j Q_i^{[-2]}`, `Q^j u² Q_i^{[-2]}`, `-Q^j u ν1^{[-2]} Q_i`, `Q^j u³ Q_i^{[-2]}` |
| `⟨F⟩` | full contour: `Σ_p F(p)/∏_{p'≠p}(p-p')` over the 8 poles `θ_α + n`, `n = 0..3` |
| row `(i,j)` | `(⟨F_1⟩, …, ⟨F_6⟩)` for dual `j` of `A`, direct `i` of `B` |

Twisted identity: `Σ_c w_c ⟨F_c⟩ = 0` for `A ≠ B`. Columns 1–3 are the untwisted matrix, because
at zero twist `Δc_22 = Δc_1` and `Δs_a = 0`.

**Equivalent form over `M = 1/(ν1ν2ν3)` (cell 13).** `M = ∏_α 1/[(u-θ_α)(u-θ_α-h)(u-θ_α-2h)]`,
6 simple poles.

| symbol | meaning |
|---|---|
| split | `ΔT_2 = ν1^{[-2]}(a_1u+b_1) + ν3(a_2u+b_2)` (unique: `ν3`, `ν1^{[-2]}` share no root) |
| `w'` | `(Δs_1, Δc_1, a_1, b_1, a_2, b_2)` |
| `G_1..G_6` | `-u Q^jQ_i`, `-Q^jQ_i`, `u Q^jQ_i^{[-2]}`, `Q^jQ_i^{[-2]}`, `(u+h) Q^{j[2]}Q_i`, `Q^{j[2]}Q_i` |
| `⟨G⟩_M` | full contour over the 6 poles `θ_α + n`, `n = 0,1,2` |

Identity: `Σ_c w'_c ⟨G_c⟩_M = 0`. The `ΔT_1` and `ν1^{[-2]}δ_1` terms sit on `M` directly
(`E_2 ν1^{[-2]} = μ/ν1 = M`); the `ν3 δ_2` term has poles `θ_α+{1,2,3}h` and is moved onto `M` by
the full-contour shift `u → u+h`. At zero twist `a_2 = -a_1`, `b_2 = Δc_1 + 6h a_1 - b_1`, so the
untwisted columns are `H = (G_2+G_6, G_3-G_5+6G_6, G_4-G_6)` for `(Δc_1, a_1, b_1)`. The `O(φ)`
directions are `Δs_1`, `a_1+a_2`, `b_1+b_2-6a_1-Δc_1`.

## 2. Symbol reference

### Cell 2 — generators
Gelfand–Tsetlin `ps, dim, J0, Jp, Jm, JJ`, copied unchanged from `SU4_L1.wb`.

### Cell 4 — twisted transfer matrices
| symbol | purpose |
|---|---|
| `E1g`, `E2g`, `Etot` | site and total generators (64×64) |
| `prm` | `{θ1 -> 1/3, θ2 -> -8/57, h -> 1}` |
| `Tm[i,j][x]` | monodromy at `prm` |
| `qminM[I]` | principal quantum minor for index set `I` |
| `qco[I,k]` | `[u^k] qmin_{I,I}` (exact rational 64×64) |
| `Tco[a,k,z]` | `Σ_{|I|=a} z_I qco[I,k]` for any `z` |
| `ν1, ν2, ν3, qdet` | `ν_i(x) = (x-θ1-λ_i)(x-θ2-λ_i)`, `qdet = ν1(x-2)ν2(x-1)ν3(x)` |

Checks at exact `z = (2,3,1/6)`: commutativity, `T_3 = qdet`, Cartan invariance (and `E_12`
non-invariance).

### Cell 6 — Baxter operators, untwisted reference, rows
| symbol | purpose |
|---|---|
| `BaxOp[T1,T2,z][q][x]` | `z ν1ν1^{[-2]}ν1^{[-4]} q^{[2]} - T1 ν1^{[-2]}ν1^{[-4]} q + z⁻¹ T2 ν1^{[-4]} q^{[-2]} - z⁻² qdet q^{[-4]}` |
| `BaxOpD[T1,T2,y][g][x]` | `ν3^{[2]}ν3^{[4]}ν3^{[6]} g - y T1^{[2]}ν3^{[4]}ν3^{[6]} g^{[2]} + y² T2^{[4]}ν3^{[6]} g^{[4]} - y³ qdet^{[6]} g^{[6]}`, `y = 1/z_j` |
| `polyKernel[op]` | exact kernel, staggering gauge |
| `wtOf`, `hwMus`, `hwBasis` | weights, the 5 highest weights, highest-weight spaces |
| `untw` | 6 untwisted records `<|"mu","T1","T2","q","qd"|>`; `qd[[j]] = q^j` (deg `μ_j+3-j`); octets ordered by `[u⁰]T2` |
| `poles`, `pw`, `rowInt[F]` | the 8 poles, their weights, the full-contour functional |
| `colFuns[Qa,Qb]` | the six integrands `F_c` |
| `rowU[A,B,{i,j}]` | untwisted 6-vector |
| `TaRows` | `{(1,1),(2,1),(2,2)}` |
| `NhatU` | vacuum-normalised untwisted 3×3 norms, `= {1, 1805/11748, 2800/2937, 4085/5874, -4085/5874, 1805/11748}` (matches companion) |

### Cell 8 — twisted states and Q-functions
| symbol | purpose |
|---|---|
| `prec`, `psi` | 150, `(1, 1/3, -4/3)` |
| `rr` | fixed rational weights for the generic combination diagonalised in a sector |
| `solveQ[op, d]` | monic degree-`d` least-squares solve; returns `{coeffs, residual/operator scale}` |
| `twState[U, eps]` | `<|"mu","eps","z","T1","T2","eigRes","q","qd","bRes","U"|>` |
| `regFun[S,"direct",k]` | `R_k ∈ span(Q_1..Q_k)` matching untwisted `q_k` at `{1/10, 7/20, 3/5}[[1;;k]]` |
| `regFun[S,"dual",j]` | `R^j ∈ span(Q^3..Q^j)` matching untwisted `q^j` at `4-j` points |
| `rawFun[S, …]` | raw `Q_k`, `Q^j` |
| `rowT[A,B,{i,j},mode]` | twisted 6-vector, `mode = "raw"` or `"reg"` |
| `wT[A,B]` | the twisted `w` vector |
| `relRes[row, w]` | `|row·w|/Σ|row_c w_c|` (0 if the row vanishes) |
| `tw[eps]` | the six twisted states (memoised) |

The regularisation is well defined for `ε → 0`: the regular combinations span the same space as
`Q_1..Q_k`, and the untwisted values at the matching points form an invertible matrix, so
`R_k = q_k + O(ε)` uniformly on compacts.

### Cell 10 — block structure
| symbol | purpose |
|---|---|
| `epsList` | `{10⁻³, 10⁻⁴, 10⁻⁵}` |
| `rowsAt[eps]` | all diagonal (`A = B`) rows, raw and reg, all `(i,j)`, all states (memoised) |
| `slope[n,mode,ij]` | `ε`-power per column (from `10⁻⁴ → 10⁻⁵`) |
| `detU` | untwisted `det` of `T_a` rows, columns 1–3 |
| `goodA[eps,n]` | reg `T_a` rows |
| `blockRatio[eps,n,tr]` | `det M_6 / (det A · det D)` for raw rows `tr` |
| `surv`, `domTriples` | survey of 84 raw triples; the 76 with dominance |
| `detD[eps,n,tr]` | lower-right `det D` |
| `powTab` | powers of `det D` for `i<j` and `i>j` triples, with `Σ_{i<j}(μ_i-μ_j)` |

### Cell 13 — rows over `M = 1/(ν1ν2ν3)`
| symbol | purpose |
|---|---|
| `polesM`, `pwM`, `rowIntM[F]` | the 6 poles of `M`, their weights, the full-contour functional |
| `colFunsM[Qa,Qb]` | the six integrands `G_c` |
| `splitT2[dT2]` | `{a_1, b_1, a_2, b_2}` from the `u⁰..u³` coefficients of `ΔT_2` |
| `wM[A,B]` | the vector `w'` |
| `rowTM[A,B,{i,j},mode]` | twisted 6-vector over `M`, `mode = "raw"` or `"reg"` |
| `rowUM[A,B,{i,j}]` | untwisted 6-vector over `M` |
| `colsH[g]` | untwisted columns `(g_2+g_6, g_3-g_5+6g_6, g_4-g_6)` |
| `arbCheck[S,{i,j}]` | `E_2` form vs `M` form for fixed arbitrary `ΔT_1` (deg 1), `ΔT_2` (deg 3) |
| `detRatio[eps,n,tr]` | `det M_6` over `M` / `det M_6` over `E_2`, rows `T_a` (reg) + raw `tr` |
| `NhatM` | untwisted `N̂` from `H`-columns over `M` |

## 3. Results

### Twisted system (validated)
- Baxter residuals ~`10⁻¹³¹` with `deg q_k = 4-μ_k`, `deg q^j = μ_j`. Eigenvector residuals ~`10⁻¹⁴³`.
- Twisted `T → ` untwisted `T` linearly in `ε` (0.036 at `ε = 10⁻³`).
- **Identity:** `|M^{AB}·w|/Σ|…| ≲ 10⁻¹¹⁷` for all 30 cross pairs, all 9 raw rows and the 3 reg rows.
  So `F_1..F_6` is the complete twisted system in 6 differences.
- `max|w_{4,5,6}| / max|w_{1,2,3}|` = `9.5·10⁻⁴`, `9.5·10⁻⁵` at `ε = 10⁻³, 10⁻⁴`.

### ε-profile of raw rows (vacuum)
A raw row `(i≠j)` scales across columns 1..6 as `ε^{(a, a, a+1, a-1, a-2, a-2)}`. For example
`(1,2)`: `4,4,5,3,2,2`; `(1,3)`: `6,6,7,5,4,4`; `(2,1)`: `0,0,1,-1,-2,-2`; `(3,1)`: `-2,-2,-1,-3,-4,-4`.
Diagonal raw rows `(i,i)` vanish identically in columns 1–4. Columns 4–6 are therefore always the
least suppressed.

### The minor
`det A / det N_untw - 1` (max over states) = `0.093, 0.0093, 0.00093` at `ε = 10⁻³, 10⁻⁴, 10⁻⁵`.

### Block dominance
`det M_6^{AA} = det N^{AA}_untw · det D^{AA} · (1 + O(ε))`:

- **Holds** for 76/84 triples of raw rows, for all 6 states. Selection criterion: residual
  `< 10⁻³` at `ε = 10⁻⁵` and at least 5× smaller than at `10⁻⁴`. The largest residual among the 76 is
  `9.1·10⁻⁴`; for the `i<j` and `i>j` triples it is ≤ `2.3·10⁻⁴`.
- **Fails** for 8 triples:
  - six share a direct index `i` or a dual index `j`: `{(1,1),(1,2),(1,3)}`, `{(2,1),(2,2),(2,3)}`,
    `{(3,1),(3,2),(3,3)}`, `{(1,1),(2,1),(3,1)}`, `{(1,2),(2,2),(3,2)}`, `{(1,3),(2,3),(3,3)}`;
  - `{(1,1),(1,2),(2,3)}`;
  - `{(1,1),(2,2),(3,3)}`, which has `det D = 0` (diagonal raw rows vanish in columns 1–4).

### Order of the lower-right block

| state | `Σ_{i<j}(μ_i-μ_j)` | power, rows `i<j` | power, rows `i>j` |
|---|---|---|---|
| `{4,2,0}` | 8 | 10 | -6 |
| `{3,3,0}` | 6 | 8 | -4 |
| `{4,1,1}` | 6 | 8 | -4 |
| `{3,2,1}` ×2 | 4 | 6 | -2 |
| `{2,2,2}` | 0 | 2 | 2 |

`det D ∝ ε^{2 ± Σ_{i<j}(μ_i-μ_j)}`. The `±Σ` part is the monic normalisation of raw `Q`'s (the
`∏(φ_i-φ_j)^{μ_i-μ_j}` factor); the universal `ε²` is the suppression carried by the three `O(φ)`
integrals of motion. The rows `i<j`, `{(1,2),(1,3),(2,3)}`, vanish at zero twist.

### Measure `1/(ν1ν2ν3)` (validated)
- Rewriting `E_2 → M` holds term by term for arbitrary `ΔT_1`, `ΔT_2`, all 9 raw `(i,j)`, all 6
  states: relative difference ~`10⁻¹²⁵` (`ε = 10⁻³`).
- Identity `⟨G⟩_M·w' = 0`: ~`10⁻¹¹⁸` for all 30 cross pairs, 9 raw + 3 reg rows.
- `O(φ)` directions / `max|w'|` = `0.0070`, `0.00070` at `ε = 10⁻³, 10⁻⁴`.
- `det M_6[M] / det M_6[E_2] = -28512/361 ≈ -78.98` for all 6 states, both `i<j` and `i>j` triples
  (`ε = 10⁻⁴`). `w ↦ w'` is a fixed linear map, so all block-structure results above carry over
  unchanged.
- Untwisted: `NhatM === NhatU` exactly.

## 4. Open
1. Leading coefficient of `det D` on several rays `ψ`: product of `(φ_i-φ_j)` powers?
2. Normalised twist rows making `det D` state-independent, so that
   `det M_6/det M_6^{vac} → N̂` directly.
3. Relation between `det D` and the Hilbert norms `H` of the companion notebook.
