# `SU3_L2_Baxter_Twisted_Final.wb` — companion reference

gl(3) XXX chain at **L = 2**, both sites in the adjoint `λ = {2,1,0}`, diagonal **unimodular** twist.
Shows that the functional SoV determinant, normalized by the product of the L=1 reference
coefficients, equals the physical scalar product of eigenvectors normalized against the product
reference states, for **all 64x64 pairs** and **both** a/b branches. The equality holds for every
admissible choice of Q-functions in the determinant rows (36 triples), with a triple-independent
vacuum constant.

- **L = 1 parent:** `SU3_L1_Baxter_Twisted_Final.wb` (exact, symbolic, 8x8).
- **Untwisted / older L = 2 relatives:** `SU3_Adjoint_SoV_ScalarProduct_L2.wb` (zero twist, 3x3 norm),
  `SU3_Adjoint_Twisted_SoV_L2.wb` (6x6 twisted norm, 3x3 minor = untwisted norm, eps-expansion).
- **gl(4) template for the determinant layout:** `../Clean/SU4_SoV_ScalarProduct_L2.wb` (`spMatFSoV`,
  rows `(α,i,j)`, `perFac`).
- **Status:** cells 1–8 verified from a cold kernel with message policy `stop` (no kernel messages).
  Cells 9–10 were run in the same warm kernel immediately afterwards (no messages). Numerical at
  90 digits; exact rational chain point and twist.
- **Runtime:** a few minutes for cells 1–8; cell 10 about 60 s.

## 0. Read this first

1. **The twist must be unimodular**, `z1 z2 z3 = 1`. The test twist is `{5/2, 7/3, 6/35}`. (The
   pre-existing `{5/2, 7/3, 1/6}` has product 35/36 and was replaced.) The Baxter equations still
   carry `(Times @@ zTwist) qdetL2`, which is now just `qdetL2`.
2. **Every bracket is one contour around both ladders**, all six poles `θ_s + nh`, `n = 0,1,2`, of
   `M = 1/(ν1ν2ν3)`. Rows are **not** per-site residue sums. The second copy of the row triple
   carries the `h`-periodic factor `e^{2πiu/h}` instead.
3. **The reference coefficient is a product of L=1 coefficients**, one per inhomogeneity. It is *not*
   extracted from the determinant, and nothing sums over a basis of wave functions.
4. **`bra02`/`ket02` must be flattened.** `KroneckerProduct` of two vectors is their 8x8 outer
   product, not a 64-vector. The earlier version of cell 3 had this bug, so its physical side
   could never evaluate (`Dot::dotsh`). `Flatten` (row-major) matches the `Esite1 = Kron[E, id]`
   ordering.
5. **The notebook shares kernel K2** with other notebooks in this directory; restart before a run.
6. Wolfbook replies truncate at `{`. Read results through stored variables, as strings with braces
   replaced.

## 1. Conventions

| symbol | meaning |
|---|---|
| `θ = (1/3, -8/57)`, `h = 1` | chain point (`thetaValues`, `hValue`) |
| `z = (5/2, 7/3, 6/35)` | unimodular twist (`zTwist`) |
| `T_ij(u)` | `L^{(1)}(u-θ1) L^{(2)}(u-θ2)`, `L_ij(u) = u δ_ij - h E_ji` |
| `T_a(u)` | `Σ_{|I|=a} z_I qmin_{I,I}(u)`; `T_3 = qdet` |
| `ν_i(u)` | `∏_s (u - θ_s - hλ_i)`, `λ = {2,1,0}` |
| `Q_i(u)` | `z_i^{u/h} q_i(u)`, `q_i` monic of degree `4 - μ_i` (direct) |
| `Q^j(u)` | `z_j^{-u/h} q^j(u)`, `q^j` monic of degree `μ_j` (dual) |
| `f^{[±2]}` | `f(u ± h)` |
| `<0_1|` | `<1,0,…,0| e^{a12 E12} e^{a13 E13} e^{a23 E23}` |
| `|0_1>` | `e^{b21 E21} e^{b31 E31} e^{b32 E32} |1,0,…,0>` |
| `<0_2|`, `|0_2>` | `<0_1| ⊗ <0_1|`, `|0_1> ⊗ |0_1>` |

### Physical scalar product

Right eigenvectors `R_B` and left eigenvectors `L_A = (R^{-1})_A` (so `L_A R_B = δ_AB`):

```
S_vec(A,B) = <L_A|R_B> / ( <L_A|0_2> <0_2|R_B> )
```

This is independent of the eigenvector gauge.

### FSoV determinant

```
<f>_α = Σ_{p ∈ poles of M} Res_{u=p} M(u) e^{2πi(α-1)u/h} f(u)

rows    (α, i, j),  α ∈ {1,2},  (i,j) ∈ triple
columns (r, β),     β ∈ {0,1}:  u^β · ( -Q^j_A Q_{B,i},  Q^j_A Q_{B,i}^{[-2]},  Q^{j[2]}_A Q_{B,i} )
```

The column order `(β=0: r=1,2,3; β=1: r=1,2,3)` fixes an overall sign, which `C_vac` absorbs.

### Normalization for the triple `{(a,c),(a,d),(b,c)}`

Direct content `{a,a,b}`, dual content `{c,c,d}`:

```
Ψ_B^ref = ∏_α Q_{B,a}(θ_α+h) [ Q_{B,a}(θ_α) Q_{B,b}(θ_α-h) - Q_{B,b}(θ_α) Q_{B,a}(θ_α-h) ]

c_A^ref = ∏_α (1/4h^6) [ Q^c_A(θ_α+h) ( Q^c_A(θ_α) Q^d_A(θ_α+h) - Q^d_A(θ_α) Q^c_A(θ_α+h) )
                        - Q^c_A(θ_α)   ( Q^c_A(θ_α) Q^d_A(θ_α+2h) - Q^d_A(θ_α) Q^c_A(θ_α+2h) ) ]

S_Q(A,B) = C_vac det N^{AB} / ( c_A^ref Ψ_B^ref ),   S_Q(vac,vac) = 1
```

The L=1 triple `{(1,1),(1,2),(2,1)}` is `a=c=1`, `b=d=2`.

### a/b branches (the L=1 relations, both sites)

```
b21 = -z1(z1-z3) / (a12 z2 (z2-z3)),     b32 = -z2^2 (z1-z3) / (a23 (z1-z2) z3^2)

Branch I:   b31 = -z1^2/(a13 z3^2),   a12 a23/a13 = -z2(z1-z3)/(z1(z2-z3))
Branch II:  b31 = -z1/(a13 z3),       a12 a23/a13 = -(z1-z3)/(z2-z3)
```

For the numerical check, `a13 = 3/5`, `a23 = 4/7`, and `a12` is fixed by the branch constraint.

## 2. Symbol reference

### Cell 2 — generators and twisted chain
| symbol | purpose |
|---|---|
| `ps, dim, J0, Jp, Jm, JJ` | Gelfand–Tsetlin construction (same as the other gl(3)/gl(4) notebooks) |
| `rep2`, `dsite`, `dchain`, `idsite` | `{2,1,0}`, 8, 64, `IdentityMatrix[8]` |
| `EElocal[i,j]` | 8x8 site generator |
| `Esite1`, `Esite2`, `Etot2` | 64x64 site-1, site-2 and total generators |
| `thetaValues`, `hValue`, `zTwist` | the chain point and twist |
| `Tmon[i,j][x]` | monodromy entries (64x64, numeric point) |
| `qminor2[I]`, `qcoeff2[I,k]` | principal quantum minor, its `u^k` coefficient |
| `Tcoeff2[a,k]` | `[u^k] T_a` including twist factors |
| `nu1L2`, `nu2L2`, `nu3L2`, `qdetL2` | `ν_i` and the quantum determinant |

Checks: twist unimodular; all `T_1, T_2` coefficients commute; `T_3 = qdet`.

### Cell 3 — eigenstates and physical scalar product
| symbol | purpose |
|---|---|
| `precisionL2` | 90 |
| `weightList2`, `sectors2` | total weights, the distinct weight sectors |
| `rho2[a,k]` | fixed rational weights of the generic combination diagonalized per sector |
| `stateRecords2` | 64 records `<|"mu","vec","T1","T2","eigResidual"|>`, sorted by weight (descending lexicographic) |
| `rightMatrix2`, `leftMatrix2` | eigenvectors as columns; its inverse |
| `vacuumState2` | index of the `{4,2,0}` state |
| `alphaConstraints2`, `betaRules2`, `alphaValues2[q]`, `betaValues2[q]` | branch data |
| `bra01[q]`, `ket01[q]` | one-site reference covector/vector (numeric) |
| `bra02[q]`, `ket02[q]` | flattened product references (64-vectors) |
| `rightReferencePhysical2[q]`, `leftReferencePhysical2[q]` | `<0_2|R_B>`, `<L_A|0_2>` |
| `physicalScalar2[q]` | `S_vec` as a 64x64 matrix (diagonal, because `L_A R_B = δ_AB`) |

### Cell 4 — Baxter polynomials
| symbol | purpose |
|---|---|
| `BaxterDirect2[T1,T2,z][q][x]` | `z ν1ν1^{[-2]}ν1^{[-4]} q^{[2]} - T1 ν1^{[-2]}ν1^{[-4]} q + z^{-1} T2 ν1^{[-4]} q^{[-2]} - z^{-2} (z1z2z3) qdet q^{[-4]}` |
| `BaxterDual2[T1,T2,y][g][x]` | `ν3^{[2]}ν3^{[4]}ν3^{[6]} g - y T1^{[2]} ν3^{[4]}ν3^{[6]} g^{[2]} + y^2 T2^{[4]} ν3^{[6]} g^{[4]} - y^3 (z1z2z3) qdet^{[6]} g^{[6]}`, `y = 1/z_j` |
| `solveMonic2[op, d]` | monic degree-`d` least-squares solve; returns `{coefficients, relative residual}` |
| `statesQ2` | `stateRecords2` plus `"q"` (direct, `deg 4-μ_i`), `"qd"` (dual, `deg μ_j`), `"baxterResidual"` |
| `qRight2[S,i,x]` | `z_i^x q_i(x)` |
| `qLeft2[S,j,x]` | `z_j^{-x} q^j(x)` |

### Cell 5 — product right reference wave function
| symbol | purpose |
|---|---|
| `gtPatterns2` | the 8 one-site labels |
| `referenceConfiguration2` | `{{1,0,0},{1,0,0}}` |
| `psiRightLocal2[S,th,{aa,bb,cc}]` | one-site `Q_1(θ+aa)[Q_1(θ+bb)Q_2(θ+cc-1) - Q_2(θ+bb)Q_1(θ+cc-1)]` |
| `psiRightReference2` | `Ψ_B^ref` for all 64 states (L=1 triple) |

### Cell 7 — determinant and the 64x64 match (L=1 triple)
| symbol | purpose |
|---|---|
| `measurePoles2`, `measureWeights2` | the 6 poles of `M` and their residue weights |
| `rowTriple2`, `rowKeys2` | `{(1,1),(1,2),(2,1)}`; the six keys `{α,{i,j}}` |
| `periodicFactor2[α,x]` | `e^{2πi(α-1)x/h}` |
| `determinantRow2[A,B,{α,{i,j}}]` | the 6-vector row |
| `determinantMatrix2` | `det N^{AB}`, 64x64 |
| `crossRelative2` | `max|off-diagonal det| / min|diagonal det|` |
| `cReferenceLocal2[S,th]`, `cReference2` | L=1 coefficient at `th`; its product over sites |
| `rawQScalar2`, `vacuumConstant2`, `qScalar2` | unnormalized, `C_vac`, normalized `S_Q` |
| `matchResidual2`, `diagonalRelative2` | `max|S_Q - S_vec|`, max relative diagonal residual (per branch) |

### Cell 10 — all row triples
| symbol | purpose |
|---|---|
| `rowCache2[[A,B]]` | association of all 18 rows `{α,{i,j}}` for the pair `(A,B)` |
| `triples2` | the 36 triples `{(a,c),(a,d),(b,c)}`, `a≠b`, `c≠d` |
| `psiRefLocalGen2[S,th,a,b]`, `cRefLocalGen2[S,th,c,d]` | normalization factors for general content |
| `logAbs2[x]` | `log10|x|`; an arbitrary-precision zero `0``p` gives `-accuracy` |
| `tripleData2[tr]` | `<|"triple","ok","cross","diagSpread","Cvac","res"|>` |
| `tripleTable2` | the 36 records |

## 3. Results

### L=1 triple (cell 7, cold run)

| check | value |
|---|---|
| simultaneous-eigenvector residual | ~1e-88 |
| Baxter residual (direct and dual, all 64) | 0 to ~79 digits |
| max off-diagonal det / min diagonal det | 3e-66 |
| min \|c_A^ref\| | 2.5e-13 (all 64 nonzero) |
| `S_Q(vac,vac)` | 1 |
| max \|S_Q - S_vec\|, branch I / II | 4e-76 / 4e-76 |
| max relative diagonal residual, branch I / II | 0 to ~72 digits |
| range of \|S_Q(A,A)\| | 1.5e-8 … 1.9e3 |

### All 36 triples (cell 10)

| check | value |
|---|---|
| triples with no zero norm, `Ψ^ref` or `c^ref` | 36/36 |
| worst relative 64x64 residual (both branches) | 1e-78 |
| worst max\|cross det\|/min\|diag det\| | 1e-34 |
| worst \|C_vac / C_vac(L=1 triple) - 1\| | 1e-75 |
| widest spread min/max \|diag det\| | 1e-51 |

**`C_vac` is the same for every triple.** All triple dependence sits in `c^ref Ψ^ref`, exactly as
the Q-content of the rows dictates.

The widest spread belongs to `{(1,2),(1,3),(3,2)}` and `{(2,1),(2,3),(3,1)}`. For 46 states (the
vacuum included) the diagonal determinant is ~1e-34 relative, and `Ψ^ref`, `c^ref` are each ~1e-17
relative, so the ratio is well-conditioned: ~77 digits remain and the match is exact to working
precision. These are not failures. A relative cutoff (`< 1e-40`) wrongly flagged them in a first
version; the cell now rejects only exact zeros.

## 4. Pitfalls hit

- **Distinct global rows fail.** Six distinct `(i,j)` rows cannot carry the direct content
  `{1,1,1,1,2,2}` of the product reference; the extracted reference coefficient vanished for all 84
  distinct sextuples. The fix is repeating the triple with the periodic label α.
- **Do not sum over wave functions.** An earlier draft normalized via `det N . Inverse[W]`; rejected
  and removed.
- **Do not divide by the left product reference wave function.** It vanishes for some states.
- **Non-unimodular twist.** The earlier test twist had `z1z2z3 = 35/36`.
- **`KroneckerProduct` of vectors** (see §0.4).
- **Relative-magnitude cutoffs** on diagonal determinants (see §3).

## 5. Open

1. ~~**Untwisting.**~~ **Done (2026-09-16)** in `SU3_L2_Untwisting.wb`; see its `.md` and `.tex`. The
   twist factor is `V_{λvac}(z)/V_{λA}(z)`. The limit is `B_{λA}/B_{λvac}`, with
   `B_λ = E_λ det𝒩/(c^ref Ψ^ref)` built from untwisted q's (Appendix N mechanism). The right and left
   eigenvectors have separate untwisted limits.
2. A symbolic (exact) version of the 64x64 match, as at L=1.
