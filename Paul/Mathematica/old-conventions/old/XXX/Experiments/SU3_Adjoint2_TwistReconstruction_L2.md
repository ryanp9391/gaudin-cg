# `SU3_Adjoint2_TwistReconstruction_L2.wb` — reference

gl(3) XXX chain at **L = 2**, **companion twist**, twists fully general (`chi3 = det G` *not* set
to 1). One question, answered for five pairs of site reps, with the mechanism derived and a cheap
fix found.

> Given `<Om(w)| tau_a(y,u) |Psi(z)>` for **all** vacuum twists `w`, **all** transfer-matrix twists
> `y`, **all** spectral parameters `u`, and `a = 0,1,2` (`tau_0 = Id`), how much of `|Psi(z)>` is
> determined?

**Short answer.** Not all of it. The data spans exactly the Cartan isotypic component
`V_{lamA+lamB}` plus everything one root step from it. The reason is that `tau_1` and `tau_2` are
`adj (+) 1` tensor operators of the diagonal gl(3) and `<Om(w)|` sits in a single isotypic
component. **Rotating the vacuum on one site only, with one extra continuous parameter, completes
the space.**

---

## 0. Read this first

1. **The question is `Psi`-free.** The equations `<Om(w)|tau_a(y,u)|v> = <Om(w)|tau_a(y,u)|Psi_A(z)>`
   are linear in `v`, so the solution set is `|Psi_A(z)> + W^perp` with
   `W = span{<Om(w)|tau_a(y,u)}`. Only the right-hand side depends on `A`. The whole question is the
   **rank of a covector family**; no eigenvector at twist `z` is ever built, and the answer is one
   number, the same for every state.
2. **Nothing is carried symbolically.** `chi` and `u` enter through a finite spanning list (section 6),
   and `w` is sampled. Both reductions are *verified*, not assumed.
3. **`<Om(w)|` needs no diagonalisation.** The vacuum eigenvalue is known in closed form, so
   `<Om(w)|` is an exact rational null vector, one linear solve per twist (~30 ms at dim 64).
4. **Exact rational arithmetic throughout.** There are no residuals; the orthogonality relations are
   literal zeros. Where a rank is computed mod a large prime it is a *certificate*, since
   `rank mod p <= rank over Q`.
5. **Kernel traps.** Cell 17 leaves the kernel at the second chain point `theta = (1/5, 46/55)` and
   clears the memoised `Tmat`, `qminM`, `kappa`, `Op`; re-run cells 5, 7, 9, 11, 13, 15 to get back.
   The rep-generic section is keyed on `(lamA, lamB, thA, thB)` and is unaffected. After any restart,
   cells 2-9, 20-21 and 24 restore the definitions.
6. Wolfbook replies truncate at `{`; read results through stored variables with braces replaced.

---

## 1. Results at a glance

| `lamA (x) lamB` | dim | `dim W` | `dim W^perp` | invisible irreps |
|---|---|---|---|---|
| `[2,1,0] (x) [1,0,0]` | 24 | 24 | **0** | none |
| `[2,1,0] (x) [2,1,0]` | 64 | 63 | **1** | `(2,2,2)` — the singlet |
| `[3,0,0] (x) [3,0,0]` | 100 | 63 | **37** | `(4,2,0):27`, `(3,3,0):10` |
| `[3,1,0] (x) [2,1,0]` | 120 | 111 | **9** | `(3,3,1):6`, `(3,2,2):3` |
| `[3,1,0] (x) [3,1,0]` | 225 | 171 | **54** | `(4,4,0):15`, `(4,3,1):15` x2, `(4,2,2):6`, `(3,3,2):3` |

`W^perp` is `gl(3)`-invariant in every case — a sum of *whole* isotypic components.

**With one per-site rotation parameter added, `adj (x) adj` goes 63 -> 64: complete.**

---

## 2. The rule

### 2.1 Sharp form

`<Om(w)|` has support on **exactly one** isotypic component, `V_{lamA+lamB}`. Therefore

```
W  =  V_{lamA+lamB}  (+)  (everything ONE ROOT STEP away),   with multiplicity
```

Verified against all five cases:

| `V1 (x) V2` | Cartan component + one root step | total | measured |
|---|---|---|---|
| `[2,1,0](x)[1,0,0]` | `(3,1,0)15 + (2,2,0)6 + (2,1,1)3` | 24 | 24 |
| `[2,1,0]^2` | `(4,2,0)27 + (4,1,1)10 + (3,3,0)10 + (3,2,1)16` | 63 | 63 |
| `[3,0,0]^2` | `(6,0,0)28 + (5,1,0)35` | 63 | 63 |
| `[3,1,0](x)[2,1,0]` | `(5,2,0)42 + (4,3,0)24 + (4,2,1)x2 30 + (5,1,1)15` | 111 | 111 |
| `[3,1,0]^2` | `(6,2,0)60 + (6,1,1)21 + (5,3,0)42 + (5,2,1)x2 48` | 171 | 171 |

### 2.2 Coarse form (the `lambda_1` shadow)

Because a root step changes `lambda_1` by at most one and `lamA_1 + lamB_1` is maximal,

```
lam is VISIBLE   <=>   lam_1  >=  lamA_1 + lamB_1 - 1
```

i.e. only the top two levels of the `lambda_1` grading survive. This is how the rule was first
found, and it agrees with the sharp form on all five cases. **Use the root-step form** — it is the
derived one; `lambda_1` could in principle differ for a rep pair not tested here.

### 2.3 Two things the rule is *not*

- **Not about singlets.** `[3,1,0]` has Dynkin labels `(2,1)`, is not self-dual, and `[3,1,0]^2`
  contains **no singlet** — yet loses 54 directions. The coincidence `dim W^perp = #singlets` in
  `adj (x) adj` is an accident of that case: the singlet is simply the only component more than one
  root step out.
- **Not the rectangular-vs-generic distinction.** `[3,0,0]` is rectangular and loses 37 of 100. This
  obstruction is independent of the one `HANDOFF.md` records for localization in generic reps.

---

## 3. Why: `adj (+) 1` tensor operators on a one-component seed

### 3.1 The operators

The quantum minors `qmin_{I,K}(u)` with `|I| = |K| = a` transform under the **diagonal** `gl(3)` in
`Lam^a(C^3) (x) Lam^a(C^3)*`:

| `a` | `Lam^a(C^3)` | transforms as |
|---|---|---|
| 1 | `3` | `3 (x) 3bar = adj (+) 1` |
| 2 | `3bar (x) det` | `3bar (x) 3 = adj (+) 1` — the `det` factors cancel |
| 3 | `det` | `1 (x) 1` — pure singlet, hence `tau_3` is a multiple of `Id` |

Measured: the `gl(3)` module generated under `ad` by the nine `T_ij(u)` at one `u` closes at
**dimension 9 with exactly 1 invariant**, and likewise for the nine `2x2` quantum minors. The
singlet of the `a = 1` module is `Sum_i T_ii(u)`, the **untwisted** transfer matrix.

A twist only chooses *which element* of that fixed 9-dimensional module you take; it cannot leave
the module. That is why a general non-companion twist buys nothing (5.2).

### 3.2 The seed

`<Om(w)|` lies in the single component `V_{lamA+lamB}` — verified as a Casimir statement. All six
polarized vacua do too, which is why none of them helps (5.3).

### 3.3 Conclusion

By Wigner-Eckart an `adj (+) 1` operator connects `V_mu` to `V_nu` only when `nu - mu` is a root or
zero. The roots of `gl(3)` are `e_i - e_j`, so **one insertion moves at most one root step**. From a
one-component seed, `k` insertions reach `k` steps. The `adj (x) adj` singlet `(2,2,2)` dies because
it is **two** steps from `(4,2,0)` — one insertion cannot reach it, two can, which is exactly what
5.1 measured.

### 3.4 Still unexplained

The per-twist rank is **9** in every rep tested. `dim(adj (+) 1) = 9` is suggestive but the
identification is *not* established: one 9-dimensional module applied to `<Om(w)|` gives only **4**
— the same as `gl(3)` itself gives — and 9 is reached only by the union over all `a`, `y`, `u`.

---

## 4. Where free covectors sit

| covector | isotypic support (`adj^2`) |
|---|---|
| `<hw (x) hw|`, `<lw (x) lw|`, `<Om(w)|` | Cartan component `(4,2,0)` **only** |
| `<hw (x) lw|`, `<lw (x) hw|` | **every** component, singlet included |
| `<Om(w)| (exp(c E^(1)_{31}) (x) 1)` | **every** component |

The vacuum is the one covector concentrated in a single component. That is the entire problem, and
`<hw (x) lw|` — a free unit vector — already has support everywhere.

---

## 5. What fixes it

### 5.1 More insertions — WORKS

```
<Om(w)| tau(y1,u1) ... tau(yk,uk) |Psi(z)>     reaches k root steps
```

Measured on `[3,0,0]^2`:

| insertions | `k <= 1` | `k <= 2` | `k <= 3` |
|---|---|---|---|
| measured rank | 63 | 90 | **100** |
| predicted | 63 | 90 | 100 |

`T_{a,s}(u)` with a larger auxiliary rep is a Cherednik-Bazhanov-Reshetikhin determinant of the
`T_{a,1}`, i.e. exactly a product, so "the rest of the fusion hierarchy" is the same fix said
invariantly, and it is what FSoV does. **Cost: `k` grows linearly with the rep size.**

### 5.2 A general, non-companion twist `y` — DOES NOT HELP

All nine `T_ij(u)` and all nine `2x2` quantum minors — 54 operators instead of 32 — still annihilate
the whole 37-dimensional `W^perp`. Explained by 3.1: the module is twist-independent.

### 5.3 The other polarized vacua — DOES NOT HELP

All six Weyl images give a genuine one-dimensional eigen-covector, with eigenvalue
`Lam^sigma_a(w,u) = Sum_{|I|=a} (Prod_{i in I} w_{sigma(i)}) kappa_I(u)` — still closed-form, no
Bethe equations. But all six are Cartan-component covectors, so all six share the same `W^perp`.

### 5.4 One per-site rotation parameter — WORKS, AND IS CHEAP

Rotate the vacuum on **one site only**: `<Om(w)| (g1 (x) 1)` with `g1 = exp(c E^(1)_{31})`. This is
an exponential of an adjoint tensor operator, so it takes arbitrarily many root steps and has full
isotypic support. Measured on `adj (x) adj`, rank vs number of sampled `c` values:

| # rotation parameters | 1 (`c = 0`) | 2 | 3 | 4 |
|---|---|---|---|---|
| rank | 63 | **64** | 64 | 64 |

**One nonzero `c` completes the space.** `O(1)` extra data, not `O(dim V)`. And `rotations/` section
19 established that per-site rotations **are** Q-determinants, so these overlaps are computable —
this is the bridge to that thread.

---

## 6. How the infinite data is made finite

`tau_a(y,u)` is polynomial in `u` of degree `2a` and depends on `y` only through `chi`. Both the
`1x1` and `2x2` minors of a companion matrix are **linear** in `chi`, so each `tau_a` decomposes
over exactly four monomials `{1, c1, c2, c3}`:

```
opList = { Id }  u  { Op[a,m][y_p] },   Op[a,m][y] = Sum_{K,I} [chi^m](det G_{K,I}) qminM[I,K][y]
         1  +  4 monomials x 3 u-points (a=1)  +  4 x 5 (a=2)   =   33 operators
```

Both reduction steps are checked at random `chi` and random `u`:

1. `tau_a(chi,y) == Sum_m chi^m Op[a,m][y]` exactly;
2. `tau_a(chi,u) == Sum_p Lagrange_p(u) tau_a(chi,y_p)` exactly.

---

## 7. Conventions

| symbol | meaning |
|---|---|
| `th1 = 1/3`, `th2 = -8/57`, `hh = 1` | chain point (cell 17 repeats at `(1/5, 46/55)`) |
| `L_ij(u) = u d_ij - h E_ji` | Lax; `T_ij(u) = Sum_k L^(1)_ik(u-th1) L^(2)_kj(u-th2)` |
| `Gsym = {{c1,-c2,c3},{1,0,0},{0,1,0}}` | companion twist, `det(t-G) = t^3 - c1 t^2 + c2 t - c3` |
| `chiOf[w]` | `(e1(w), e2(w), e3(w))`; **no unimodularity imposed** |
| `qminM` / `qmG` | quantum minors, exact rational sparse, memoised on `y` |
| `tauM` / `tauG` | `Sum_{|I|=|K|=a} qmin_{I,K}(u) det G_{K,I}` (Cauchy-Binet) |
| `nu[i][y]` | `(y-th1-h lamA_i)(y-th2-h lamB_i)` |
| `kappa` / `kapG` | `kappa_i = nu_i`, `kappa_{ij}(u) = nu_i(u-h) nu_j(u)` |
| `vacEig` / `vacEG` | `Sum_{|I|=a} (Prod_{i in I} w_i) kappa_I(y)` |
| `omegaVec` / `omegaG` | `<Om(w)|`, the vacuum **left eigenvector** at twist `w` |
| `opList` / `tagsG` | the 33-element spanning list of `{Id, tau_1(y,u), tau_2(y,u)}` |
| `uPQ` | the rep-generic `u`-points, 3 for `tau_1`, 5 for `tau_2` (defined in cell 21) |

---

## 8. Notebook map

| cell | contents |
|---|---|
| 2 | `ps`, `dim`, `J0`, `Jp`, `Jm`, `JJ` — GT construction, rank-independent |
| 3 | `lam1`, `lam2`, `d1`, `d2`, `dtot`, `E1g`, `E2g`, `Etot`, `glRelQ` |
| 5 | `th1`, `th2`, `hh`, `subsetsA`, `Tmat`, `qminM`, `nu`, `hwIndex`, `ehw`, `kappa` |
| 7 | `Gsym`, `c1..c3`, `minorGsym`, `tauM`, `tauDiagM`, `chiOf`, `vacEig` |
| 9 | `uPts`, `minRules`, `monoms`, `coeffKI`, `Op`, `lagr`, `opList`, `opTag` |
| 11 | `tauFromOps`, `omegaVec`, `wSamples`, `omegaSamples`, `nW` (79) |
| 13 | `rowsFor`, `idIdx`, `a1Idx`, `a2Idx`, `rkOf`, `satCurve`, `rankW` |
| 15 | `nullBasis`, `nullVec`, `singletNS`, fresh-point validation, `sOrbit`, `U0` |
| 17 | `runAtTheta`, `thetaCheck` — second chain point |
| 20 | `prim`, `repDim`, `EAg`, `EBg`, `EtG`, `TmG`, `qmG`, `hwOf`, `kapG` |
| 21 | `uPQ`, `opG`, `tauG`, `vacEG`, `omegaG`, `tagsG` |
| 22 | `analyse` — rank of `W` for an arbitrary rep pair |
| 24 | `dimOfLam`, `c2OfLam`, `lamTable`, `lamOfC2`, `irrepSplit`, `splitByRule` |
| 25 | the five-rep table |
| 27 | the three candidate fixes on `[3,0,0]^2` |
| 29 | `adClose`, `invCount` — the mechanism (`adj (+) 1`, one-component seed) |
| 31 | `proj`, `support`, `g1`, `seedRot` — free covectors and the rotation fix |

Cell 7 also confirms: all `tau_a` coefficients commute; `tau_3 = chi3 qdet`; companion spectrum =
diagonal spectrum (so the vacuum eigenvalue may be read off in the diagonal frame); **simple
spectrum**; and the companion twist genuinely breaks the Cartan, so no residual weight grading is
hiding the result.

---

## 9. Evidence

- Exact rationals everywhere — no floating point, no tolerances.
- `adj (x) adj` validated on **400** `(w, y, u)` points never used to build the span, plus 12 fresh
  `w` for `<Om(w)|s> = 0`.
- `<Om(w)|` confirmed to be a left eigenvector at `u`-points *not* used to construct it, with a
  one-dimensional null space at every sampled twist, pairwise non-proportional across `w`.
- Rebuilt from scratch at `theta = (1/5, 46/55)` with fresh `u`-points and fresh twists (cell 17):
  same `dim W = 63`, same null direction.
- Saturation curves plateau well before the last sample in every case (`adj^2`: 63 from 8 twists,
  flat to 79; `[3,1,0]^2`: 171 from 25 twists, flat to 36).
- Where computed mod a large prime, the exact rank was recomputed and agreed in every case.

---

## 10. Traps and mistakes made

- **`Join @@ Table[..., {a,..}, {b,..}]` does not stack rows.** With two iterators `Table` returns a
  nested list and `Join @@` gives a list of *matrices*. Use `Join @@ Flatten[Table[...], 1]` (level 2
  for three iterators). This bit **four** separate times, twice silently leaving a `MatrixRank`
  unevaluated.
- **`Position[tags, {0,_,_}]` searches all levels** and matched the monomial `{0,0,0}` *inside* a
  tag, scrambling the filtration. Select on `tags[[All,1]]` instead.
- **`Coefficient[expr, 1]`** is not how you get a constant term; use `CoefficientRules`.
- **Module-local result variables are lost.** A name declared in the `Module` variable list cannot be
  read back afterwards.
- **Hoist `om.O` out of inner loops**, and reduce rows to primitive integers with `prim` before
  `MatrixRank`. Two computations had to be killed for want of this — exact-rational rank on
  unreduced rows does not converge in reasonable time at dim 100+.
- **`C_2` is NOT injective on partitions of fixed size.** At `|lam| = 6` both `(3,3,0)` and `(4,1,1)`
  give `C_2 = 24`. That is the `adj (x) adj` case, where it happens not to matter (both visible, and
  `W^perp` there is the unambiguous `C_2 = 12`); `|lam| = 4, 7, 8` have no collisions. `lamOfC2` in
  cell 24 takes the first match and would **silently mislabel** a collision — check before reusing it
  on a new rep pair.
- **Predicting `dim W^perp` needs multiplicities.** Summing `dimOfLam` over partitions with small
  `lam_1` gives 39, not 54, for `[3,1,0]^2`, because `(4,3,1)` occurs twice.

---

## 11. Open

1. **Push the rotation fix to the hard cases.** `[3,1,0]^2` and `[3,0,0]^2` must lose 54 and 37
   directions respectively; check that one per-site rotation parameter still suffices there, or find
   how many are needed as a function of the rep. **This is the first thing to do.**
2. **Which generator, and which site?** The fix was tested with `exp(c E^(1)_{31})`. `rotations/`
   sections 17 and 19 found `E_32` easy to realise as a Q-determinant and `E_21`/`E_31` hard. The
   practically decisive question is whether the **cheap** `E_32` rotations already suffice.
3. **The ladder, made precise.** Seeds `<Om(w)| (E^(1)_{31})^k` sit `k` root steps down; since one
   insertion covers one step, seeds at `k = 0, 2, 4, ...` suffice, i.e. `ceil(steps/2)` families.
   Measure how `steps` grows for `[S,S,0]^2` to get the real scaling against SoV's `dim V`.
4. **`<hw (x) lw|` as a seed family.** One unit vector gives one dimension; find what parametrises it
   cheaply, and whether its overlaps are computable.
5. **What is the 9?** See 3.4.
6. **Multiplicity.** The reach hits whole isotypic components *including multiplicity* in all five
   cases. Check this survives at larger reps — if only part of a multiplicity space were reached, the
   counting in 2.1 breaks.
7. **`L > 2` and `gl(N > 3)`.** Everything here is `L = 2`, `gl(3)`. The mechanism
   (`Lam^a (x) Lam^a*` is `adj (+) 1` for all `0 < a < N`) has an obvious `gl(N)` statement; nothing
   here tests it.
8. **Does the CG application need the low part at all?** If the coefficients of interest live in the
   Cartan component and its neighbours, the obstruction costs nothing and no fix is needed.
