# `SU3_Adjoint_SoV_ScalarProduct_L2.wb` — companion reference

Functional SoV scalar product at **zero twist** for the gl(3) XXX chain with `L = 2` and **adjoint**
sites `[2,1,0] ⊗ [2,1,0]` (non-rectangular). Orthogonality of transfer-matrix eigenstates is proved
by a 3×3 determinant. The notebook also compares the SoV norm with the Hilbert-space norm of the
explicit eigenvectors.

- **Derivation / write-up:** `SU3_Adjoint_SoV_ScalarProduct_L2.tex` (same directory).
- **Twisted continuation:** `SU3_Adjoint_Twisted_SoV_L2.wb` / `.md` / `.tex`. There the 3×3 below
  is identified as a minor of the twisted 6×6.
- **Template:** `SU4_SoV_Orthogonality_Untwisted_L2.wb` (gl(4), rectangular sites) and
  `untwisted-norm.md` at the repo root.
- **Status:** every result below comes from a cold top-to-bottom run (fresh kernel, no kernel
  messages). Sections 1–6 are exact. Section 7 is exact except the Q-value search, which is a
  30-digit numerical screen.
- **Runtime:** a few minutes. The slow cells are 12 (816-triple survey) and 18/22 (three parameter
  points, symbolic residues).

## 0. Read this first

1. **Name clash, fixed.** The measure is `E2meas`. An earlier version called it `E2`, which is
   also the site-2 generator `E2[i,j]`. `ClearAll[E2]` then silently destroyed the generators, and
   `hwBasis` returned wrong spaces. Do not reuse `E1`, `E2`, `Etot`, `tt`, `tco` as names.
2. **Parameter point.** `θ1 = 1/3, θ2 = -8/57, h = 1`. The two octets are split by
   `3h³√(4h² + 5(θ1-θ2)²)`, so a generic rational point gives irrational Q-functions. `θ1-θ2 = 9/19`
   makes the root rational (`43/19`) and keeps `θ1-θ2 ∉ hℤ`, so the pole ladders stay separate.
   Further rational points: `t = (θ1-θ2)/h ∈ {9/19, -7/11, -33/19}`. General family:
   `t = (m-1)(m-5)/(m²-5)`.
3. **Memoisation.** `tco[a,k]` and `honest[s]` are memoised. Re-run from cell 4 if the chain changes.
4. **Degenerate states.** At zero twist the chain is gl(3)-invariant. Only highest-weight vectors
   are separated; descendants share eigenvalues.

## 1. Conventions

| symbol | meaning |
|---|---|
| `f^{[±2n]}(u)` | `f(u ± n h)` |
| `L_ij(u)` | `u δ_ij - h E_ji` (per site) |
| `T(u)` | `L^{(1)}(u-θ1) L^{(2)}(u-θ2)` in auxiliary space |
| `ν_i(u)` | `(u-θ1-hλ_i)(u-θ2-hλ_i)`, `λ = {2,1,0}` |
| `T_{a,1}` | `Σ_{|I|=a} qmin[I,I](u)`; `T_3 = ν1^{[-4]} ν2^{[-2]} ν3` (quantum determinant) |
| `O` | `1 - T1/ν1 D^{-1} + T2/(ν1 ν1^{[-2]}) D^{-2} - T3/(ν1 ν1^{[-2]} ν1^{[-4]}) D^{-3}` |
| `O†` | `1 - T1^{[2]}/ν3^{[2]} D + T2^{[4]}/(ν3^{[2]}ν3^{[4]}) D² - T3^{[6]}/(ν3^{[2]}ν3^{[4]}ν3^{[6]}) D³` |
| `μ(u)` | `∏_α 1/((u-θ_α)(u-θ_α-h))`, proportional to `Γ[ν1^{[2]}]/Γ[ν3^{[2]}]` |
| `E_2` | `μ/(ν1 ν1^{[-2]}) = 1/∏_α∏_{n=0}^{3}(u-θ_α-nh)` |
| `q_i` | direct solutions, `O q_i^{[2]} = 0`, `i = 1,2,3` by increasing degree |
| `q^j` | dual solutions, `O† q^j = 0`, `deg q^j = μ_j + 3 - j` (`q^3` lowest) |
| `L_s[·]` | sum of residues of `E_2·(·)` over the ladder `θ_s + nh`, `n = 0..3` |
| `Δt1, Δt21, Δt20` | `ΔT1`, `[u¹]ΔT2`, `[u⁰]ΔT2` (the three independent differences) |

## 2. Symbol reference

### Cell 2 — generators
`ps`, `dim`, `J0`, `Jp`, `Jm`, `JJ[k,l,rep]`: Gelfand–Tsetlin construction, ported unchanged from
`SU4_L1.wb` (identical to `Gaudin/Clean/Hamiltonians_SU4.wb`). It is rank-independent; the cell checks
`dim[{2,1,0}] = 8` and all gl(3) commutation relations.

### Cell 4 — chain
| symbol | purpose |
|---|---|
| `rep`, `d1`, `id1` | `{2,1,0}`, 8, `IdentityMatrix[8]` |
| `E1[i,j]`, `E2[i,j]`, `Etot[i,j]` | site-1, site-2, total generators on the 64-dim space (memoised) |
| `Tm[i,j][x]` | monodromy entries, symbolic in `x, θ1, θ2, h` |
| `qmin[I,K][x]` | quantum minor `Σ_σ sgn σ ∏_k T_{I_σk, K_k}(x-(k-1)h)` |
| `tt[1]`, `tt[2]`, `tt[3]` | transfer matrices as 64×64 matrices polynomial in `u` |

### Cell 5 — checks
`ν[i]`, `ν1`, `ν2`, `ν3`, `qdet`; `tco[a,k]` = `[u^k] tt[a]` (memoised). Three checks, all `True`:
`[T_a(u),T_b(v)] = 0`, `T_3 = ν1^{[-4]}ν2^{[-2]}ν3`, and `[T_a, E^{tot}_{ij}] = 0`.

### Cell 7 — highest-weight states
| symbol | purpose |
|---|---|
| `prm` | `{θ1 -> 1/3, θ2 -> -8/57, h -> 1}` |
| `strChk` | operator identities: `[u⁴]T2 = 3`, `[u³]T2 = -6(3h+θ1+θ2)`, `[u²]T2 = [u⁰]T1 + 3(10h²+θ1²+3θ1θ2+θ2²+8h(θ1+θ2))`; `[u¹]T2`, `[u⁰]T2` are not in `span{[u⁰]T1, 1}` |
| `wtOf` | 64×3 list of total weights |
| `hwBasis[w]` | basis of `ker E^{tot}_{12} ∩ ker E^{tot}_{23}` in weight `w` (64×k) |
| `hwStates` | list of records `<|"mu","copy","vec","T1","T2"|>` |
| `eigChk` | every `"vec"` is an eigenvector of the full `tt[1]`, `tt[2]` at `prm` |

Highest weights: `{4,2,0}` (27), `{3,3,0}` (10), `{4,1,1}` (10̄), `{3,2,1}` ×2 (8, 8), `{2,2,2}` (1).
`T1` is the same on 10 and 10̄ and on both octets; `T2` separates them.

### Cell 9 — Q-functions
| symbol | purpose |
|---|---|
| `Dmax` | 12, degree bound for the kernel search |
| `BaxRes[s][q][x]` | `ν1ν1^{[-2]}ν1^{[-4]} q^{[2]} - T1 ν1^{[-2]}ν1^{[-4]} q + T2 ν1^{[-4]} q^{[-2]} - T3 q^{[-4]}` |
| `BaxResD[s][g][x]` | `ν3^{[2]}ν3^{[4]}ν3^{[6]} g - T1^{[2]}ν3^{[4]}ν3^{[6]} g^{[2]} + T2^{[4]}ν3^{[6]} g^{[4]} - T3^{[6]} g^{[6]}` |
| `polyKernel[res]` | exact polynomial kernel, row-reduced from the top (staggering gauge), ascending coefficient lists sorted by degree |
| `qpoly[cl, y]` | polynomial from coefficient list |
| `states` | `hwStates` records plus `"q"` (direct) and `"qd"` (duals, **sorted by increasing degree**, so `qd[[4-j]] = q^j`) |
| `qChk` | all residuals vanish identically |

Degrees (`deg q_1..q_3`; `deg q^1, q^2, q^3`):

| state | direct | dual |
|---|---|---|
| `{4,2,0}` | 0, 3, 6 | 6, 3, 0 |
| `{3,3,0}` | 1, 2, 6 | 5, 4, 0 |
| `{4,1,1}` | 0, 4, 5 | 6, 2, 1 |
| `{3,2,1}` (both) | 1, 3, 5 | 5, 3, 1 |
| `{2,2,2}` | 2, 3, 4 | 4, 3, 2 |

So `deg q_i = (4-μ_i) + (i-1)` and `deg q^j = μ_j + (3-j)`.

### Cell 11 — measure and rows
| symbol | purpose |
|---|---|
| `mMeas[x]`, `E2meas[x]` | `μ` and `E_2` (symbolic) |
| `ladder[s]`, `Lsum[f,s,x]` | the 4 poles of ladder `s` at `prm`; residue sum |
| `qB[B,i][y]`, `qA[A,j][y]` | direct `q_i` of `B`; dual `q^j` of `A` |
| `rowR[A,B,s,i,j]` | 3-vector `(L_s[E_2 q^j_A(-ν1^{[-2]}q^B_i + u² q_i^{B[-2]})], L_s[E_2 q^j_A u q_i^{B[-2]}], L_s[E_2 q^j_A q_i^{B[-2]}])` |
| `dInv[A,B]` | `(Δt1, Δt21, Δt20)` |
| `d2Chk` | `ΔT1` is constant and `[u²]ΔT2 = ΔT1`, all pairs |
| `fullId[A,B,s,i,j]` | direct check `L_s[μ q^j_A O^A q_i^{B[2]}] = 0` (poles `n = 0..4`) |
| `idTab`, `rowTab` | 540/540 and 540/540 (30 ordered cross pairs × 2 ladders × 9 `(i,j)`) |

### Cell 12 — per-ladder triples
`rowKeys` = 18 triples `(s,i,j)`; `normRow`; `zeroRows` (empty); `triples` (816); `normSurvey`;
`goodTriples` (732: nonzero norm for all 6 states; the other 84 are zero for all 6); `Nmat[A,B,t]`;
`orthTab[t]`. At the first good triple `(1,1,1),(1,1,2),(1,2,1)`: norms 6/6 nonzero, cross
determinants 30/30 zero, `N·Δ = 0` for 30/30, no identically zero cross matrix.

### Cell 15 — explicit matrix (Section 6)
`rowsSimple`, `bracket[f,x]` (independent implementation of the displayed contour formula, ladder 1
only), `Nexplicit[A,B]`, `oct1`, `oct2` (copies 1, 2 in cell-7 order), `sameAsRowR` (`True`, 36 pairs).

With `⟨f⟩ = ∮_{C_1} du/(2πi) f/∏_α∏_{n=0}^{3}(u-θ_α-nh)`:

```
N^{AB} = [ <q^1_A (u² q_1^{B[-2]} - ν1^{[-2]} q_1^B)>   <q^1_A u q_1^{B[-2]}>   <q^1_A q_1^{B[-2]}> ]
         [ <q^2_A (u² q_1^{B[-2]} - ν1^{[-2]} q_1^B)>   <q^2_A u q_1^{B[-2]}>   <q^2_A q_1^{B[-2]}> ]
         [ <q^1_A (u² q_2^{B[-2]} - ν1^{[-2]} q_2^B)>   <q^1_A u q_2^{B[-2]}>   <q^1_A q_2^{B[-2]}> ]
```

Octet example: `T(8_1) - T(8_2) = (0, 0, -129/19)`,
`det N^{8_1 8_1} = 53291113880694676031/410708765061000000`, `det N^{8_1 8_2} = 0`,
`rank N^{8_1 8_2} = 2`. For `8_1`: `q_1 = u - 434/285`, `q_2 = u³ - 41u²/19 + 49084/48735`,
`q^1 = u⁵ - 305u⁴/57 + 325295u²/9747 - 151058891/8333685`, `q^2 = u³ - 84u²/19 + 104791/48735`.

### Cell 18 — full contour, three points (Section 7)
| symbol | purpose |
|---|---|
| `pts` | `P1 = (1/3, -8/57, 1)`, `P2 = (1/5, 46/55, 1)`, `P3 = (1/4, 107/76, 2/3)` for `(θ1, θ2, h)` |
| `mkStates[pp]` | cells 7 + 9 pipeline at an arbitrary point (`Block[{prm = pp}, …]`) |
| `st[P]` | states at point `P`, octets ordered by `[u⁰]T2`, key `"n"` = position |
| `vacOf[P]` | the `{4,2,0}` record |
| `rowFull[A,B,i,j]` | `rowR` summed over both ladders (full contour) |
| `fullRows`, `fullZero` | diagonal full rows; `(i,j)` with `i<j` vanish identically (18/18) |
| `fullTriples`, `fullSurvey`, `fullGood` | 84 triples; 19 nonzero for all 18 state×point; 65 zero for all |
| `crossRows`, `fullOrth` | cross rows at P1; all 19 good triples give 30/30 cross `det = 0` |

### Cell 20 — gauge-invariant triples, Hilbert norm
| symbol | purpose |
|---|---|
| `Ta`, `Tb` | `{(1,1),(2,1),(2,2)}`, `{(2,2),(3,2),(3,3)}` |
| `NhatF[tr,P,n]` | `det tr^{AA}/det tr^{vac}` |
| `pivotN[v]` | first nonzero component set to 1 |
| `honest[s]` | `<|"dims","R","L","H","Rhw"|>`: right/transpose-left joint null vectors of all `T_a` coefficients in the weight sector, pivot gauge, `H = Ψ_L·Ψ_R` |
| `normTab` | per point and state: dims, `Ψ_R` = hw vector, `H`, `N̂(T_a)`, `N̂(T_b)`, `H/N̂` |

Residual gauge: `q_2 → q_2 + a q_1`, `q_3 → q_3 + b q_1 + c q_2`, `q^1 → q^1 + d q^2 + e q^3`,
`q^2 → q^2 + f q^3`. Using `R_{i<j} ≡ 0`, only `{R11,R22,R33}` (degenerate), `T_a` and `T_b` are
invariant triples. `N̂(T_b)` is `N̂(T_a)` with 10 ↔ 10̄ exchanged and the octet signs flipped.

### Cell 22 — symbolic `t`, Q-values, search
| symbol | purpose |
|---|---|
| `mkStatesW[pp, ws]` | non-octet pipeline, 1-dim hw spaces only |
| `stSym` | vac, 10, 10̄, 1 at `θ1 = t, θ2 = 0, h = 1` |
| `symRows`, `symNhat` | symbolic `T_a` rows and `N̂(t)` |
| `casD`, `casU`, `monicP` | monic Casoratians `f(u)g(u-1)-f(u-1)g(u)`, `f(u)g(u+1)-f(u+1)g(u)` |
| `qfuns[s]`, `qfSym` | `q1, q2, q12, q13, d1, d2, d3, d12, d23, d13` (`d_j = q^j`, `d3` lowest) |
| `symChk` | `θ1↔θ2` symmetry at ladder points: true exactly for `q1, q12, d3, d23` |
| `xs`, `tvs`, `logv`, `tgt`, `dualPairs`, `dirPairs`, `searchHits` | the bounded search (Section 3 below) |

## 3. Results

### Orthogonality (exact)
- Per ladder: 732/816 triples work for all states; at the simplest one, norms 6/6 nonzero and cross
  determinants 30/30 zero.
- Full contour: `R_{i<j} ≡ 0`; 19/84 triples work; all 19 orthogonal.

### Vacuum-normalised norm vs Hilbert norm

| state | `N̂(T_a)` (symbolic, `t = (θ1-θ2)/h`) | `H = Ψ_L·Ψ_R` |
|---|---|---|
| `{4,2,0}` | 1 | 1 |
| `{3,3,0}` | `20/(3(49-25t²))` | 2 |
| `{4,1,1}` | `160(t²-1)/(3(25t²-49))` | 2 |
| `{2,2,2}` | `20/(3(49-25t²))` | 8 |

Numerical values at P1: `1, 1805/11748, 2800/2937, ±4085/5874, 1805/11748`.

- Octets (numerical, not symbolic): `H = 645/128, 86/3` (P1), `405/32, 162/25` (P2), `830/27, 249/50` (P3).
  They depend on θ because pivot gauge is not adapted to the 2-dim octet space.
- **No triple gives a θ-independent `N̂`** (0/19). This is unlike gl(4).

### Search (numerically rejected within the stated ansatz)
Ansatz: `N̂_A Π_vac/Π_A = c_A H_A`, with `c_A` t-independent, where
`Π = q^a(x1) q^{bc}(x2) q_1(x3) q_12(x4)`:
- `a ∈ {1,2,3}`, `bc ∈ {12, 23, 13}`;
- `x_i ∈ {θ_α + nh : n = -2..4}`;
- both `Π_vac/Π_A` and its inverse;
- 1254 × 168 × 2 candidates, 30 digits, 4 values of `t`.

**No hit.** Structural reason: the gauge-invariant `q_1, q_12, q^3, q^{23}` have only trivial roots
(`θ_α`, `θ_α+2h`) and midpoint roots `(θ1+θ2+kh)/2`. At ladder points they give products of `(t±k)/2`
and cannot produce `7 ± 5t`.

## 4. Pitfalls
- Wolfbook replies truncate at `{`: store results in variables and read small pieces.
- `Position[v, _?(# =!= 0 &), {1}, 1]` returned 0 here. Use `Position[Unitize[v], 1]`.
- `"T1" -> expr &[p]` parses as a pure function of the rule. Use a `Module` local instead.
- `Table[..., {tv, P} ∈ {...}]` is not a valid iterator. Use `{pr, {{…},{…}}}`.
- Residue sums over a zero vacuum value: guard the division before taking logs.

## 5. Open
1. Why the per-ladder 84 and full-contour 65 degenerate triples fail (cf. gl(4) `i ≥ j` rule).
2. A θ-independent row combination, or `u^β` insertion rows, that could match `H`.
3. The octets symbolically (rational parametrisation `m`) with an octet-adapted Hilbert gauge.
4. See the twisted notebook for the minor identification.
