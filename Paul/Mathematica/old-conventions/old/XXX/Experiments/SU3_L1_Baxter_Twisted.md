# `SU3_L1_Baxter_Twisted.wb` — companion reference

gl(3) XXX chain at **L = 1**, single site in the adjoint `lambda = {2,1,0}`, inhomogeneity `theta[1]`,
**diagonal** twist `K = diag(z1,z2,z3)`. Builds the twisted fused transfer matrices, the Baxter
operators `O`, `O†` as written in (12.38)–(12.39), the vacuum `Q`-functions, the SoV covector basis,
and expands the abstract `3×3` SoV scalar-product determinant into wave functions.

- **Companion `.nb`:** `su3_210_L1_companion.nb` (same directory) — source of the rescaled-Lax /
  fused-transfer-matrix / SoV-basis conventions, ported in cells 10–12.
- **L = 2 relatives:** `SU3_Adjoint_SoV_ScalarProduct_L2.wb`, `SU3_Adjoint_Twisted_SoV_L2.wb`.
- **Status:** cold top-to-bottom run from a fresh kernel completes with no errors. The only kernel
  messages are two benign `Solve::svars` (cells 13, 14) — `Solve` reporting the expected
  one-parameter normalisation freedom in the Baxter solution.
- **Everything is exact and symbolic** in `u, th, h, z1, z2, z3`. The single exception is cell 15,
  which is a deliberately independent check at random rational parameters.
- 30 cells: 21 code, 9 markdown.

## 0. Read this first

1. **Two parameter namings coexist, on purpose.** Cells 3–6 use `th, z1, z2, z3`; the ported
   companion cells 10–12 keep the companion's own `theta[1], z[i]`. The rule `bridge` maps the
   latter onto the former. Nothing is assigned across the two — `bridge` is only ever applied as a
   replacement, so neither set is clobbered.
2. **`tG` is NOT a Baxter input.** `tG[a] = T_{a,1}/(nu_1 nu_1^{[-2]} ...)` is already divided by the
   `nu_1` ladder, and `O` divides by it again. Feed `O` the polynomial `TT[a]` / `TvacDiag[a]`.
   Cell 11 checks both that the relation holds and that substituting `tG` gives a nonzero residual.
3. **The SoV pseudovacuum is generic**, `xv[1],...,xv[8]`. The `{0,...,0,1}` in the companion `.nb`
   came from an earlier notebook and is not the right object here.
4. **`x` carries DownValues** (`x[{{1},{0,0}}] = ...` etc. in cell 12). Do not `ClearAll[x]` after
   cell 12, and do not reuse `x` as a scratch symbol.
5. **`ns` runs over `{0,1,2}^3`, not over the 8 GT labels.** The determinant expansion of cells
   20–21 needs the full lattice; the 8 Gelfand–Tsetlin labels span only an 8-dimensional subspace
   whereas the determinant's direct side is 16-dimensional.
6. Wolfbook replies truncate at `{`; read results from the stored variables rather than the cell
   output when querying from an agent.
7. The notebook shares kernel `K2` with other notebooks in this directory. It defines `nu1, nu2,
   nu3, qdet, T, t, x, G` among others; run it in a clean kernel.

## 1. Conventions

| symbol | meaning |
|---|---|
| `f^{[2n]}(u)` | `f(u + n h)`; `D f(u) = f(u+h)`, `D^{-1} f(u) = f(u-h)` |
| `rep`, `Lam` | `{2,1,0}`, the single-site highest weight (both names, cells 3 and 10) |
| `T_ij(u)` | `(u - theta) delta_ij - h E_ji`, one site |
| `qmin_{I,K}(u)` | `sum_sigma sgn(sigma) prod_k T_{I_sigma(k), K_k}(u - (k-1)h)` |
| `T_{a,1}(u)` | `sum_{|I|=a} z_I qmin_{I,I}(u)`, `z_I = prod_{i in I} z_i`; degree `a` |
| `nu_i(u)` | `u - theta - h lambda_i`; `nu_1 = u-th-2h`, `nu_2 = u-th-h`, `nu_3 = u-th` |
| `qdet` | `T_{3,1} = z1 z2 z3 nu_1^{[-4]} nu_2^{[-2]} nu_3` |
| `tG[a](u)` | rescaled: `T_{a,1}(u) / (nu_1 nu_1^{[-2]} ... nu_1^{[-2(a-1)]})`, **not polynomial** |
| `M(u)` | `mu/nu_1 = 1/(nu_1 nu_2 nu_3)`, three simple poles `th, th+h, th+2h` |
| `x_ab` | separated variable at `theta[1] + n_ab h`; labels `{{n21},{n11,n12}}` |
| `au[i,n]`, `ad[j,n]` | abstract `Q_i(th + n h)` and `Q^j(th + n h)` |
| `ns`, `ms` | right / left SoV labels, each a triple in `{0,1,2}^3` |

**Baxter operators.**

```
O      = 1 - T1/nu1 D^-1 + T2/(nu1 nu1^[-2]) D^-2 - T3/(nu1 nu1^[-2] nu1^[-4]) D^-3
O†     = 1 - T1^[2]/nu3^[2] D + T2^[4]/(nu3^[2] nu3^[4]) D^2 - T3^[6]/(nu3^[2] nu3^[4] nu3^[6]) D^3
```

acting as `O qhat_i^{[2]}(u) = 0` and `O† qhat^j(u) = 0`.

**Wave functions.**

```
Psi_B({n21},{n11,n12}) = Q_1(x21) [ Q_1(x11) Q_2(x12 - h) - Q_2(x11) Q_1(x12 - h) ]
Psi_A({n21},{n11,n12}) = Q^1(x21) [ Q^1(x11) Q^2(x12 + h) - Q^2(x11) Q^1(x12 + h) ]
```

`Psi_A` is the mirror of `Psi_B`: `Q_i -> Q^i` and `-h -> +h`. With these, `Psi_B` samples the direct
`Q` at `n = -1,0,1,2` and `Psi_A` samples the duals at `n = 0,1,2,3` — exactly the two sets the
residues of `M` produce.

## 2. Symbol reference

### Cell 2 — gl(N) generators (`---- 0 ----`)
Ported verbatim from `SU4_L1.wb`; rank-independent.

| symbol | purpose |
|---|---|
| `ps[rep]` | Gelfand–Tsetlin patterns of `rep`, as `{{m11},{m12,m22},{lambda}}`, in basis order (memoised) |
| `dim[rep]` | `Length[ps[{rep}]]` (memoised); `dim[{2,1,0}] = 8` |
| `J0[k][pat]` | Cartan weight component `k` of a pattern |
| `Jp[k][pat]`, `Jm[k][pat]` | raising / lowering action, returned as a combination of `v[pattern]` |
| `JJ[k,l,rep]` | the matrix of `E_kl` in the GT basis (memoised); off-diagonal cases built by nested commutators |
| `v[pat]` | formal basis-vector head used by `Jp`, `Jm` |

Checks: `dim = 8`; the printed GT patterns and their weights; all gl(3) commutation relations.
**Basis order matters:** pattern 1 is `{{2},{2,1},{2,1,0}}`, of weight `(2,1,0)` — the highest
weight vector is basis vector 1.

### Cell 3 — chain (`---- 1 ----`)

| symbol | purpose |
|---|---|
| `rep`, `dd`, `idd`, `zmat` | `{2,1,0}`, `8`, `IdentityMatrix[8]`, the `8×8` zero matrix |
| `EE[i,j]` | `JJ[i,j,rep]`, the site generators (memoised) |
| `zz` | `{z1, z2, z3}`, the diagonal twist |
| `Tm[i,j][x]` | monodromy `(x - th) delta_ij idd - h EE[j,i]` |
| `qmin[II,KK][x]` | quantum minor, `k`-th factor at `x - (k-1)h` |
| `TT[a][x]` | twisted fused transfer matrix `sum_{|I|=a} z_I qmin[I,I][x]` (memoised) |
| `nu[i][x]`, `nu1`, `nu2`, `nu3` | `x - th - h rep[[i]]` and its three instances |
| `qdet[x]` | `z1 z2 z3 nu1[x-2h] nu2[x-h] nu3[x]` |

Checks: degrees `{1,2,3}`; `T_3 = qdet`; all `T_a` coefficients mutually commute; all commute with
the Cartan; `T_1` does **not** commute with `E_12` (confirming the twist is acting).

### Cell 4 — Baxter operators (`---- 2 ----`)

| symbol | purpose |
|---|---|
| `shift[f,n][x]` | `f[x + n h]`, i.e. `D^n` |
| `OpO[T1,T2,T3][f][x]` | `O` as in (12.38), applied to `f`; call it on `u -> qhat(u+h)` |
| `OpOdag[T1,T2,T3][g][x]` | `O†` as in (12.39) |
| `OpOpoly[...]`, `OpOdagPoly[...]` | the same two, cleared of denominators (multiplied by `nu_1 nu_1^{[-2]} nu_1^{[-4]}` and `nu_3^{[2]} nu_3^{[4]} nu_3^{[6]}`) |

Checks: the shifted `nu`'s are displayed explicitly, and each cleared form is verified equal to the
operator times its ladder, for symbolic `T1, T2, T3` and an arbitrary function.

### Cell 5 — the vacuum by direct diagonalization (`---- 3 ----`)

| symbol | purpose |
|---|---|
| `wts` | `8×3` table of Cartan weights of the basis vectors |
| `hwPos` | position of weight `{2,1,0}`; evaluates to `1` |
| `vac` | `UnitVector[8, hwPos]` = `(1,0,0,0,0,0,0,0)` |
| `eigQ[a]` | is `vac` an eigenvector of every `u`-coefficient of `TT[a]`? |
| `TvacDiag[a]` | the eigenvalue, read off component `hwPos` of `TT[a][u] . vac` |
| `TvacPred[a]` | the closed form `sum_{i1<...<ia} z_I prod_k nu_{I_k}(u-(a-k)h)` |

Checks: `eigQ[a]` true for `a = 1,2,3`; `TvacDiag[a] === TvacPred[a]` for all three.

### Cell 6 — Baxter check on the vacuum (`---- 4 ----`)

| symbol | purpose |
|---|---|
| `T1v, T2v, T3v` | `TvacDiag[1,2,3]`; the vacuum eigenvalues used everywhere downstream |
| `Q1vac[x]` | `z1^(x/h)` |
| `Q1sh[x]` | `Q1vac[x+h]`, the argument `O` actually takes |
| `resO`, `resOpoly` | `O Q_1^{[2]}(u)` and its cleared form, divided by `Q1vac[u]` |
| `Q3vac[x]` | `z3^(-x/h)`, the dual counterpart |
| `resOdag` | `O† Q^3(u) / Q^3(u)` |

Checks: `resO === resOpoly === 0` and `resOdag === 0`, identically in `u, th, h, z1, z2, z3`. The
dual check is what validates the `O†` implementation.

### Cell 8 — GT / SoV dictionary (`---- 5 ----`)

| symbol | purpose |
|---|---|
| `codePats` | `Most /@ ps[{{2,1,0}}]`, i.e. patterns as `{{m11},{m12,m22}}` |
| `toSoV[{{a},{b,c}}]` | `{{b},{a,c}}` — the relabelling |
| `mapped` | `toSoV /@ codePats` |
| `userList` | the eight SoV labels as written by hand |

Check: `Sort[mapped] === Sort[userList]`. The dictionary is `x21 = m12`, `x11 = m11`, `x12 = m22` —
a **permutation of roles, not a shift**. Note the pair `{x11, x12}` deliberately mixes GT levels.

### Cell 10 — companion conventions (`---- 6 ----`)
Ported from `su3_210_L1_companion.nb`, keeping its symbol names.

| symbol | purpose |
|---|---|
| `Lam` | `{2,1,0}` (the companion's `Lambda`) |
| `Lop[i,j][w]` | the companion's script-`L`: `w idd delta_ij - h JJ[j,i,Lam]` |
| `G` | `DiagonalMatrix[{z[1],z[2],z[3]}]` (kept for fidelity; `TG` uses `z[i]` directly) |
| `T[{i,j}][w]` | `Lop[i,j][w - theta[1]]`, the unrescaled monodromy |
| `TG[{i,j}][w]` | **rescaled**: `z[i] T[{i,j}][w] / (w - theta[1] - 2h)` — not polynomial |
| `QM[Ii,Jj][w]`, `QMG[Ii,Jj][w]` | quantum minors of `T` and of `TG`; the `Position` trick selects the argument `w - (k-1)h` for the `k`-th factor |
| `tG[a][w]` | `sum_{|I|=a} QMG[I,I][w]` |
| `theta[1]`, `z[i]` | the companion's parameter names |
| `bridge` | `{theta[1] -> th, z[1] -> z1, z[2] -> z2, z[3] -> z3}` |

**Deviation from the source:** the companion writes `T[{i_,j_}][u_] := T[i,j][u] = ...`, memoising
under a head that never matches again. Fixed here to `T[{i,j}][w]` so the cache works. Semantics
unchanged.

### Cell 11 — bridge, and why `tG` must stay out of Baxter (`---- 7 ----`)

| symbol | purpose |
|---|---|
| `nuLadder[a][w]` | `prod_{k=0}^{a-1} nu1[w - k h]` |
| `tGasT` | per `a`: does `tG[a] * nuLadder[a] == TT[a]`? |
| `baxWithTG` | the residual obtained by wrongly feeding `tG` into `O` |

Checks: `QM[I,I] == qmin[I,I]` for all `I`; `tGasT` true for `a = 1,2,3`; the correct Baxter residual
is `0` while `baxWithTG =!= 0`; and `tG[a]` equals `TvacDiag[a]/nuLadder[a]` on the vacuum.

### Cell 12 — fused transfer matrices and the SoV basis (`---- 8 ----`)

| symbol | purpose |
|---|---|
| `t[{1,0,0}]`, `t[{1,1,0}]` | `tG[1]`, `tG[2]` |
| `t[{2,0,0}][w]` | `tG[1][w] . tG[1][w+h] - tG[2][w+h]` |
| `t[{2,1,0}][w]` | `tG[2][w] . tG[1][w+h] - tG[3][w+h]` |
| `gtpats` | the 8 SoV labels `{{n21},{n11,n12}}` |
| `xv[i]` | components of the **generic** pseudovacuum |
| `x[label]` | SoV covectors, `1×8` row matrices, generated from `xvac` by acting with `t[...]` at `u = theta[1]` |
| `xvac` | `x[{{1},{0,0}}]` = `{Array[xv, 8]}` |
| `xmat` | the `8×8` matrix of all covectors |
| `normPsi[psi]` | `psi/First[xvac . psi]` — the convention `<x_vac|Psi> = 1` |
| `PsiVac` | `normPsi[vac]` = `(1/xv[1], 0, ..., 0)` |

Note `x[{{2},{0,0}}]` is built with an `Inverse[t[{1,0,0}][theta[1]]]`; the other three `{{2},...}`
covectors descend from it. Checks: `Dimensions[xmat] = {8,8}`, `MatrixRank[xmat] = 8` (a generic
pseudovacuum still gives a basis), and `<x_vac|Psi_vac> = 1`.

### Cell 13 — the linear equation for the polynomial part (`---- 9 ----`)

| symbol | purpose |
|---|---|
| `redBax[zk, Pf][xx]` | the Baxter equation after `Q = zk^{u/h} P`, divided by `zk^{u/h}`; linear in `P` |
| `degMax`, `cc[m]`, `Ppoly` | search degree `4` and the generic polynomial `sum_m cc[m] xx^m` |
| `res`, `eqs`, `sol`, `solSpace` | residual, the 7 linear equations, the solution, and `Ppoly` on it |

`redBax` reads

```
zk nu1 nu1^[-2] nu1^[-4] P^[2] - T1 nu1^[-2] nu1^[-4] P + (T2/zk) nu1^[-4] P^[-2] - (T3/zk^2) P^[-4]
```

Checks: `redBax[z1, 1] === 0` (validating it against the known `Q_1`); the residual has degree 6 in
`u`, giving 7 equations; for `zk = z2` the solution is one-dimensional, `deg = 1`.

### Cell 14 — the three `Q`-functions (`---- 10 ----`)

| symbol | purpose |
|---|---|
| `polyPart[zk, dmax]` | solves `redBax[zk, P] = 0` for `deg P <= dmax`, returns it **monic**; `"NO SOLUTION"` if none |
| `PP[k]` | monic polynomial part of `Q_k`, `k = 1,2,3` |
| `Qfun[k][xx]` | the full `Q_k(xx) = z_k^{xx/h} PP[k]` |
| `degTab`, `rootTab`, `baxTab` | degrees, Bethe roots, and the symbolic Baxter check |

Results:

```
P_1(u) = 1
P_2(u) = u - (z1(th+h) - z2(th+2h))/(z1-z2)
P_3(u) = (th-u)^2 + h(th-u)(z1+2z2-3z3)z3/((z1-z3)(z3-z2))
                  + h^2(z1^2 z2 - 3 z1 z2 z3 + 2 z3^3)/((z1-z3)^2(z3-z2))
```

`degTab = {0,1,2}`, matching `deg P_k = lambda_1 - mu_k = 2 - mu_k` for `mu = (2,1,0)`. The root of
`P_2` tends to `th+h` as `z2 -> 0` and to `th+2h` as `z1 -> 0`, i.e. onto the zeros of `nu_2`, `nu_1`.

### Cell 15 — independent Baxter check (`---- 11 ----`)

| symbol | purpose |
|---|---|
| `num` | `{th -> 1/3, h -> 1, z1 -> 5/2, z2 -> 7/3, z3 -> 1/6}` |
| `nuN[w]` | `nu_1` at those parameters |
| `baxDirect[zk, Pk, u0]` | the Baxter residual assembled **directly** from the transfer matrices, not via `redBax` |
| `indepChk` | `3 × 3` table of residuals, `k = 1,2,3` at `u = 2/7, -3/5, 11/4` |

All nine residuals are exactly `0`. This is the non-circular confirmation of cell 14, which would
otherwise be checking the solution against the equation it came from.

### Cell 17 — `<x|Psi>` versus the wave function (`---- 12 ----`)

| symbol | purpose |
|---|---|
| `Qred[k][n]` | `z_k^n PP[k]` at `u = th + n h` — the reduced `Q` with the common `z_k^{th/h}` dropped |
| `PsiWF[{{a},{b,c}}]` | `Qred[1][a] (Qred[1][b] Qred[2][c-1] - Qred[2][b] Qred[1][c-1])` |
| `lhsX` | `First[(x[g] /. bridge) . PsiVac]` for each label |
| `rhsWF` | `PsiWF[g]/PsiWF[{{1},{0,0}}]` |
| `xvFree` | any `xv[i]` surviving in `lhsX` — must be empty |

Dropping the common exponential is legitimate because every term of `Psi` carries exactly two `Q_1`
and one `Q_2`, so `z_1^{2 th/h} z_2^{th/h}` cancels in the ratio.

### Cell 19 — the abstract determinant and the wave functions (`---- 13 ----`)

| symbol | purpose |
|---|---|
| `Mmeas[w]` | `1/(nu1 nu2 nu3)` |
| `mPoles`, `mW` | the three poles `th + {0,1,2}h` and their residues `{1/(2h^2), -1/h^2, 1/(2h^2)}` |
| `au[i,n]`, `ad[j,n]` | abstract direct / dual `Q` values at `th + n h` |
| `gvecA[{i,j}, n]` | the three column integrands at pole `n`: `{-ad[j,n] au[i,n], ad[j,n] au[i,n-1], ad[j,n+1] au[i,n]}` |
| `rowA[{i,j}]` | `sum_n mW[[n+1]] gvecA[{i,j}, n]` — one row of the determinant |
| `detA[tr]` | `Det[rowA /@ tr]` for a triple `tr` of `(i,j)` pairs |
| `allIJ`, `trG` | all 9 `(i,j)` pairs; the triples whose `Q`-content matches `Psi` on both sides |
| `PsiBa[{a,b,c}]`, `PsiAa[{a,b,c}]` | right and left wave functions in the abstract `au`/`ad` |
| `cube` | `Tuples[{0,1,2},3]`, the 27 `ns` labels |

**`trG` is forced, not chosen.** Each row supplies one `Q_i` and one `Q^j`; matching `Psi`'s content
(two `Q_1`, one `Q_2`) requires both index multisets to be `{1,1,2}`, and of the 84 triples exactly
one is non-degenerate: `{(1,1),(1,2),(2,1)}`.

The three columns come from the `L = 1` form of the functional orthogonality identity. Since
`chi_1, chi_2` are state-independent, `Delta T_1` is constant and `Delta T_2` linear, so the split
`Delta T_2 = nu_1^{[-2]} delta_1 + nu_3 delta_2` has **constant** `delta_1 = a_1`, `delta_2 = a_2`,
and the six columns of the `L = 2` notebooks collapse to

```
G_1 = -Q^j Q_i,    G_2 = Q^j Q_i^[-2],    G_3 = Q^{j[2]} Q_i        against   w' = (Dc_1, a_1, a_2)
```

### Cell 20 — the expansion (`---- 14 ----`)

| symbol | purpose |
|---|---|
| `monsU` | the direct monomials appearing in `PsiBa` and in `detA` |
| `Amat`, `bvec` | columns = `PsiBa[ns]` for every `ns`; `bvec` = the dual-valued coefficients of `detA` |
| `cSol`, `cAssoc` | `LinearSolve[Amat, bvec]`, and the same as an association `ns -> c[ns]` |
| `suppC` | the `ns` with `c[ns] != 0` |
| `resid` | `detA - sum_ns c[ns] PsiBa[ns]`, must be `0` |

Gives `det N = sum_ns c[ns] Psi_B[ns]` exactly, with the `c[ns]` polynomials in the **dual** `Q^j`
alone. Dividing by `Psi_B({1},{0,0})` converts `Psi_B[ns]` into the normalised `<x_ns|Psi_B>` and
leaves the `c[ns]` unchanged, since the factor cancels from both sides.

### Cell 21 — the non-diagonal measure (`---- 15 ----`)

| symbol | purpose |
|---|---|
| `A27` | `PsiAa /@ cube`, the left wave functions |
| `monsD`, `Adual` | dual monomials and the corresponding coefficient matrix |
| `Mrow` | per `ns`, the solution of `Adual . m = c[ns]` |
| `Mtab` | association `{ns, ms} -> M[ns,ms]` over the nonzero entries |
| `suppM` | `Keys[Mtab]` |
| `residM` | `detA - sum M[ns,ms] PsiAa[ms] PsiBa[ns]`, must be `0` |

### Cell 23 — charge conservation (`---- 16 ----`)

| symbol | purpose |
|---|---|
| `grade` | `{au[i,n] -> t^n au[i,n], ad[j,n] -> t^-n ad[j,n]}`, i.e. the twist rescaling `z -> t z` |
| `claimOK`, `viol` | does every entry satisfy `Total[ns] == Total[ms]`, and which do not |
| `colHom` | are the three columns homogeneous of degree `0, -1, -1` for every `n`? |
| `dHom`, `bHom`, `aHom` | homogeneity of `detA`, `PsiBa[ns]`, `PsiAa[ms]` |
| `blockTab` | nonzero entries per charge sector `q = Total[ns]` |

### Cells 26–29 — all states, two-sided normalisation, and matching (`---- 17–20 ----`)

The final section implements

```
<0| = <1,0,...,0| Exp[a12 E12] Exp[a13 E13] Exp[a23 E23]
|0> = Exp[b21 E21] Exp[b31 E31] Exp[b32 E32] |1,0,...,0>
```

and diagonalises all coefficients of the three transfer matrices. The constant coefficient of
`TT[2]` resolves the two-dimensional weight-`{1,1,1}` block. For every eigenstate, cell 27 solves
all three direct and all three dual Baxter equations; all 48 residuals vanish symbolically.

For the right reference configuration `nB = {1,0,0}`, treating the sampled Q-values as formally
independent gives the off-shell coefficient

```
cOff[{1,0,0}] = (Psi_A[{1,0,0}] - Psi_A[{0,0,1}])/(4 h^6).
```

Every one of the eight on-shell dual Baxter solutions satisfies

```
Q^1(theta) Q^2(theta+h) - Q^2(theta) Q^1(theta+h) = 0,
```

so `Psi_A[{1,0,0}] = 0` for every state, not only for the vacuum. The exact on-shell coefficient
therefore has the unique dual configuration `nA = {0,0,1}`:

```
c[{1,0,0}] = -Psi_A[{0,0,1}]/(4 h^6),
Psi_A[{0,0,1}] = Q^1(theta) (
  Q^1(theta) Q^2(theta+2h) - Q^2(theta) Q^1(theta+2h)).
```

Hence

```
Nhat[s] = Cvac detN[s]/(Psi_A[s,{0,0,1}] Psi_B[s,{1,0,0}]),
Cvac = -4 h^6 z2 z3^3/((z1-z2) (z1-z3)^2 (z2-z3)).
```

This is one on the vacuum solution `Q_1 = z1^(u/h)`. All `b` parameters can be written in terms of
the `a` parameters as

```
b21 = -z1 (z1-z3)/(a12 z2 (z2-z3))
b32 = -z2^2 (z1-z3)/(a23 (z1-z2) z3^2)
```

There are two branches for `b31`, each with one constraint among the three `a` parameters:

```
I:
 b31 = -z1^2/(a13 z3^2)
 a12 a23/a13 = -z2 (z1-z3)/(z1 (z2-z3))

II:
 b31 = -z1/(a13 z3)
 a12 a23/a13 = -(z1-z3)/(z2-z3)
```

The constraint is necessary because the matching solution has two continuous Cartan-gauge
freedoms, while three `a` parameters were introduced. Both branches match all eight norms exactly.

## 3. Results

### Chain and Baxter (exact, symbolic)
- `deg T_{a,1} = a`; `T_{3,1} = z1 z2 z3 nu_1^{[-4]} nu_2^{[-2]} nu_3`; all `T_a` coefficients
  commute, commute with the Cartan, and `T_1` does not commute with `E_12`.
- The vacuum is basis vector `1`, `(1,0,0,0,0,0,0,0)`, a joint eigenvector of every `T_a`
  coefficient, with `T_{a,1}|vac = sum_{i1<...<ia} z_I nu_{i1}(u-(a-1)h) ... nu_{ia}(u)`.
- `O Q_1^{[2]} = 0` for `Q_1 = z_1^{u/h}`, and `O† Q^3 = 0` for `Q^3 = z_3^{-u/h}`, both identically.

### `Q`-functions (exact)
`deg P_k = 2 - mu_k = 0, 1, 2`; the solutions are unique up to normalisation; verified both
symbolically and by the independent numeric check of cell 15.

### `<x|Psi>` (exact, all 8 labels)

```
<x[{{n21},{n11,n12}}] | Psi> = Psi({n21},{n11,n12}) / Psi({1},{0,0})
```

| `{{n21},{n11,n12}}` | `<x\|Psi>` |
|---|---|
| `{{1},{0,0}}` | `1` |
| `{{1},{1,0}}` | `z1 + z2/2` |
| `{{1},{1,1}}` | `z1 z2/2` |
| `{{2},{0,0}}` | `z1` |
| `{{2},{1,0}}` | `z1(2 z1 + z2)/2` |
| `{{2},{1,1}}` | `z1^2 z2/2` |
| `{{2},{2,0}}` | `z1^2(2 z1 + z2)/2` |
| `{{2},{2,1}}` | `z1^3 z2/2` |

`theta[1]` and `h` cancel completely; the `n21` dependence factorises as `z1^{n21}`, as the
`Q_1(x21)` prefactor demands. The `xv[2..8]` drop out entirely — `Psi_vac` is the first basis vector,
so only `xv[1]` enters and the normalisation removes it. **This therefore does not yet test the
generic pseudovacuum**; that needs the other seven states.

### The determinant expansion (exact, all `Q` abstract)

```
det N / Psi_B({1},{0,0}) = sum_ns c[ns] <x_ns|Psi_B>        (17 nonzero c[ns])
det N                    = sum_{ns,ms} M[ns,ms] Psi_A[ms] Psi_B[ns]     (29 nonzero entries)
```

`M` is **strongly non-diagonal** — only 4 of the 29 entries have `ns = ms`. Every entry is a dyadic
rational over `h^6`, with no dependence on `theta[1]` or on the twist; the values occurring are
`{±1/8, ±1/4, ±3/8, ±1/2, 3/4} / h^6`.

| `ns` | `ms` | `h^6 M` | | `ns` | `ms` | `h^6 M` |
|---|---|---|---|---|---|---|
| `{0,0,0}` | `{0,0,0}` | `1/8` | | `{1,2,0}` | `{0,1,2}` | `-1/4` |
| `{0,0,2}` | `{0,1,1}` | `3/8` | | `{1,2,0}` | `{1,0,2}` | `1/4` |
| `{0,0,2}` | `{1,0,1}` | `-3/8` | | `{1,2,0}` | `{2,0,1}` | `-1/4` |
| `{0,1,0}` | `{1,0,0}` | `-1/4` | | `{1,2,1}` | `{0,2,2}` | `-1/8` |
| `{0,2,0}` | `{0,1,1}` | `-1/8` | | `{1,2,1}` | `{1,1,2}` | `1/2` |
| `{0,2,0}` | `{1,0,1}` | `1/8` | | `{1,2,2}` | `{1,2,2}` | `1/4` |
| `{0,2,1}` | `{0,1,2}` | `1/4` | | `{1,2,2}` | `{2,1,2}` | `-1/4` |
| `{0,2,1}` | `{1,0,2}` | `-1/4` | | `{2,0,0}` | `{0,0,2}` | `1/8` |
| `{0,2,2}` | `{1,1,2}` | `-1/2` | | `{2,0,0}` | `{1,0,1}` | `-1/8` |
| `{0,2,2}` | `{2,0,2}` | `1/8` | | `{2,1,0}` | `{1,0,2}` | `-1/4` |
| `{0,2,2}` | `{2,1,1}` | `3/8` | | `{2,1,0}` | `{2,0,1}` | `1/4` |
| `{1,0,0}` | `{0,0,1}` | `-1/4` | | `{2,2,0}` | `{0,2,2}` | `1/8` |
| `{1,0,0}` | `{1,0,0}` | `1/4` | | `{2,2,1}` | `{1,2,2}` | `-1/4` |
| `{1,0,2}` | `{1,1,1}` | `3/4` | | `{2,2,2}` | `{2,2,2}` | `1/8` |
| `{1,1,0}` | `{0,1,1}` | `1/2` | | | | |

### Charge conservation (proved, not just observed)
`M[ns,ms] != 0` only if `Total[ns] == Total[ms]`, with zero violations. Under the grading
`deg Q_i(th+nh) = n`, `deg Q^j(th+nh) = -n` (the twist rescaling `z -> t z`):

| object | degree |
|---|---|
| `G_1, G_2, G_3` | `0, -1, -1`, independent of `n` |
| `det N` | `-2` |
| `Psi_B[ns]` | `Total[ns] - 1` |
| `Psi_A[ms]` | `-(Total[ms] + 1)` |

so `Psi_A[ms] Psi_B[ns]` has degree `Total[ns] - Total[ms] - 2` and matching `-2` forces the claim.
Because `det N` is purely of grade `-2`, any expansion may be projected onto that grade and stays
valid, so a charge-conserving representative always exists.

Block sizes by `q = Total[ns]`:

| `q` | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| nonzero entries | 1 | 3 | 7 | 8 | 6 | 3 | 1 |

`q = 0` and `q = 6` are one-dimensional and diagonal, both with `h^6 M = 1/8`.

### Non-uniqueness
`Psi_A` and `Psi_B` each span only 17 of the 27 labels, so `c[ns]` and `M[ns,ms]` are **not unique**;
`LinearSolve` returns one representative and the table above is that one. What is
representative-independent: the direct side of `det N` is 16-dimensional and sits inside the
17-dimensional `Psi_B` span, and every admissible representative can be taken charge-conserving.

### All-state norm matching

The extension in cells 26–29 is exact and symbolic. The direct and dual Q-functions are derived
from each state's transfer-matrix eigenvalues. After the corrected double reference division and
vacuum calibration, the physical two-sided norms and functional norms agree on both
reference-parameter branches for all eight states, including both weight-`{1,1,1}` states.

## 4. Pitfalls

- **Never feed `tG` to `O`.** See §0.2; cell 11 guards it.
- **`Solve::svars` in cells 13 and 14 is expected**, not an error: the Baxter solution is fixed only
  up to normalisation, so one coefficient stays free. `polyPart` then divides by the leading
  coefficient to make it monic.
- **Do not `ClearAll[x]` after cell 12** — `x` holds the SoV covectors as DownValues.
- **`Together`, not `Simplify`, in cell 12.** With 8 generic `xv[i]` the covectors are linear in
  them and `Simplify` is wasted work. The companion `.nb` uses `Simplify`; this is a deliberate
  deviation.
- **Cell 21 depends on `cSol` from cell 20**, and cell 23 on `suppM` from cell 21. Run 19–23 in order.
- **`ns` is not the GT lattice.** Restricting the expansion to the 8 GT labels fails: their span is
  8-dimensional against the determinant's 16.

## 5. Open

1. Re-test the individual coordinate identity `<x|Psi> = Psi_B(x)/Psi_B({1},{0,0})` on the
   remaining seven states after substituting the exponential reference covector into `xvac`.
2. Whether charge conservation plus the block sizes pins `M` down uniquely inside each sector. The
   `q = 0` and `q = 6` blocks already are, being one-dimensional.
3. Whether the `h^{-6}` dyadic structure of `M` survives at `L = 2`, and what it becomes for a
   general single-site rep `lambda`.
