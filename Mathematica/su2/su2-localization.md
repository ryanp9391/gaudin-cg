# `su2-localization.wb` — reference

Companion reference for `su2-localization.wb`. The notebook is a **teaching notebook**: Paul wrote
it to convey the localization structure of the SoV basis, not as production machinery. This file
records what it establishes, the exact conventions, and the two places where the obvious guess is
wrong.

**Scope constraint: this project only ever deals with TWO sites.** Do not build or propose `L>2`
generalizations. Paul has confirmed the results below generalize, but that generalization is not
wanted here.

---

## 1. The result in one paragraph

For the two-site gl(2) XXX chain in the companion twist frame, the SoV `B` operator is just
`T[{1,1}][u]`, and the SoV basis is built by acting on a vacuum with the symmetric-rep transfer
matrix tower. `B` has 16 distinct eigenvalues, so its left eigenvectors are **unique up to
normalization** — there is only one SoV basis. What varies is which end of it a given construction
reaches cheaply, and correspondingly which states factorize across the two sites:

- Build **up** from the lowest-weight vacuum: the states with the **second** label unexcited
  factorize, `x[{mu,0}] = x1[{mu}] (x) vac2`, with site 2 in its minimal vacuum.
- Build **down** from the maximal `B`-eigenvector: the states with the **first** label unexcited
  factorize, `xb[{0,nu}] = vacb1 (x) x2b[{nu}]`, with site 1 in its maximal vacuum.

Both factorizations are **exact — no leftover scalar**. Nothing factorizes in between.

Localization is stronger than a statement about the vacuum: it holds on the entire orbit of the
vacuum under the **global** su(2) rotation, where the two-site transfer matrix at `w_n` acts as the
one-site transfer matrix on site `n`. See the mechanism note in section 4.

The asymmetry is not from the SoV construction (transfer matrices commute, so the order of the two
factors is irrelevant). It comes from the monodromy ordering `T = L1(u-w1).L2(u-w2)` together with
the choice of vacuum: site 1 is the site adjacent to the vacuum in that product.

---

## 2. Conventions that matter

| Convention | Value / consequence |
|---|---|
| Reps | `Lam[1] = Lam[2] = {3,0}`, so `S = 3`, `dim = 4` per site, 16-dimensional chain |
| GT ordering | `ps` returns patterns `Reverse@Sort`ed, so basis runs **highest weight first**. `Ee[1][{1,1}] = diag(3,2,1,0)` confirms it. Last basis vector = lowest weight |
| Lax index flip | `LL[n][{i,j}][u] = u KroneckerDelta[i,j] - h Ee[n][{j,i}]` — note the **transposed** generator index. Hence `T[{2,1}]` contains `E12`, and `T[{1,2}]` contains `E21` |
| Which lowers a bra | `T[{2,1}]` (contains `E12`) lowers bras; `T[{1,2}]` raises them. `T[{2,1}]` annihilates `<0|`, `T[{1,2}]` annihilates `<0bar|` |
| Twist | `G = {{z1+z2, -z1 z2},{1,0}}`, companion matrix of `(z-z1)(z-z2)`. Consequence: `t[{1,1}] = tr(G Tt)` and the SoV `B` operator is plainly `T[{1,1}][u]`, no dressing |
| `nu[1][u]` | first Yangian weight, `prod_n (u - w_n - h Lam[n]_1)` = the `B`-eigenvalue on the **maximal** vacuum |
| `nu[2][u]` | second Yangian weight, `prod_n (u - w_n - h Lam[n]_2)` = the `B`-eigenvalue on the **minimal** vacuum. For `{S,0}` this is just `prod_n (u - w_n)` |

Note the crossed pairing: construction A uses the *minimal* vacuum but normalizes by `nu[1]`;
construction B uses the *maximal* vacuum and normalizes by `nu[2]`.

---

## 3. Symbols

### Shared machinery

| Symbol | Signature | Purpose |
|---|---|---|
| `ps` | `ps[rep]` | Gelfand-Tsetlin patterns for a gl(n) diagram, highest weight first |
| `dim` | `dim[rep]` | `Length[ps[{rep}]]` |
| `J0`, `Jp`, `Jm` | `J0[k][pat]` etc. | GT matrix elements of the Cartan / raising / lowering generators |
| `JJ` | `JJ[k,l,rep]` | Full generator matrix `E_kl`; off-neighbour cases by nested commutators |
| `Ee` | `Ee[n][{i,j}]` | `JJ[i,j,Lam[n]]` — local generator on site `n` |
| `LL` | `LL[n][{i,j}][u]` | Local Lax matrix, **transposed generator index** (see above) |
| `G` | — | Companion twist |
| `T` | `T[{i,j}][u]` | Polynomial monodromy, `Sum_k LL[1][{i,k}][u-w1] (x) LL[2][{k,j}][u-w2]` |
| `BB` | `BB[u]` | The SoV `B` operator, `= T[{1,1}][u]` |
| `EE` | `EE[{i,j}]` | **Global** (coproduct) generator, `Ee[1][{i,j}] (x) 1 + 1 (x) Ee[2][{i,j}]`. Generates the rotation orbit on which localization holds |

### Construction A — build up from the minimal vacuum

| Symbol | Signature | Purpose |
|---|---|---|
| `Tt` | `Tt[{i,j}][u]` | `T[{i,j}][u] / nu[1][u]` |
| `t` | `t[{a,s}][u]` | Transfer matrix, `a` rows / `s` columns. `t[{1,1}] = tr(G Tt)`; `t[{2,1}]` is the twisted quantum determinant; `t[{1,s}]` by fusion with **upward** shifts |
| `x` | `x[{mu1,mu2}]` | SoV basis. `x[{0,0}]` = lowest weight (x) lowest weight (last basis vector); `x[{mu1,mu2}] = x[{0,0}].t[{1,mu1}][w1].t[{1,mu2}][w2]` |
| `Tt1`, `t1`, `x1` | as above, one site | Length-1 chain on **site 1**, inhomogeneity `w[1]`, rep `Lam[1]`, normalized by its own first Yangian weight |
| `vac1`, `vac2` | — | Minimal (lowest weight) vacua on each site |

### Construction B — build down from the maximal vacuum

| Symbol | Signature | Purpose |
|---|---|---|
| `nu[2]` | `nu[2][u]` | Second Yangian weight (defined in this section, not with `nu[1]`) |
| `Ttb` | `Ttb[{i,j}][u]` | `T[{i,j}][u] / nu[2][u]` — mirror of `Tt` |
| `tb` | `tb[{a,s}][u]` | Mirror transfer matrix tower; **downward** fusion shifts |
| `vacb1`, `vacb2` | — | Maximal (highest weight) vacua, i.e. the first basis vector on each site |
| `xb` | `xb[{nu1,nu2}]` | Mirror SoV basis. `xb[{0,0}] = vacb1 (x) vacb2`; `xb[{nu1,nu2}] = xb[{0,0}].tb[{1,nu1}][w1 + h Lam[1]_1].tb[{1,nu2}][w2 + h Lam[2]_1]`. Labels count **down** from `Lam[n]_1`. Memoized |
| `Tt2b`, `tb2`, `x2b` | as above, one site | Length-1 chain on **site 2**, inhomogeneity `w[2]`, rep `Lam[2]`, normalized by its own second Yangian weight |

**The mirror tower is Paul's `t` tower with exactly two changes:**

1. normalize by `nu[2]` instead of `nu[1]`;
2. run the fusion shifts downward:
   `tb[{1,s}][u] = tb[{1,s-1}][u].tb[{1,1}][u-(s-1)h] - tb[{2,1}][u-(s-2)h].tb[{1,s-2}][u]`

and evaluate at `w_n + h Lam[n]_1` (roots of the maximal eigenvalue) rather than at `w_n` (roots of
the minimal one).

---

## 4. Verified results

All checks below are symbolic in `h, w[1], w[2], z[1], z[2]`, and the whole notebook re-runs clean
top to bottom (26 code cells, ~4.5 s, no kernel messages) as of 2026-08-27.

| Claim | Check | Output |
|---|---|---|
| A-basis diagonalizes `B`, eigenvalue `prod_n (u - w_n - h mu_n)` | 4x4 grid | `{0}` |
| `x[{mu,0}] = x1[{mu}] (x) vac2`, exactly | 4 states | `{0}` |
| `x[{0,1}]` admits no rank-one form | `Solve` | `{}` |
| `vacb1 (x) vacb2` is the maximal `B`-eigenvector, eigenvalue `prod_n (u - w_n - h Lam[n]_1)` | — | `{0}` |
| B-basis diagonalizes `B`, eigenvalue `prod_n (u - w_n - h (Lam[n]_1 - nu_n))` | 4x4 grid | `{0}` |
| `xb[{0,nu}] = vacb1 (x) x2b[{nu}]`, exactly | 4 states | `{0}` |
| `xb[{1,0}]` admits no rank-one form | `Solve` | `{}` |
| `xb[{nu1,nu2}]` proportional to `x[{Lam1-nu1, Lam2-nu2}]` | `MatrixRank`, 4x4 grid | `{1}` |
| Localization on the whole global-rotation orbit, `<0| Exp[r EE[{2,1}]]` | 4 states, symbolic `r` | `{0}` |
| Same, mirrored: `xb[{0,0}] Exp[r EE[{1,2}]]` (scratch) | 4 states, symbolic `r` | `{0}` |
| Independent local rotations `r1 != r2` break it; restored at `r2 = r1` (scratch) | 4 states | fails / `{0}` |

**One-site normalizations.** For `L=1`, `BB[u] = (u-w1) I - h Ee[1][{1,1}]` is already diagonal in
the GT basis, so `x1[{mu}]` is just a rescaled basis vector. The scalars come out `1, 1/3, 1/6, 1/6`,
which is `(S-mu)!/S!` at `S=3`, with no `h`, `w` or `z` dependence. **The closed form is a
conjecture checked at `S=3` only** — changing `Lam[1] = {S,0}` for `S = 2,4` would settle it.

**Cross-relation.** In construction A's labels the mirror statement reads
`x[{3,mu2}] = (1/6) vacb1 (x) x2[{mu2}]`, the `1/6` being `(S-mu1)!/S!` at `mu1 = S`.

**Mechanism: the protected set is the global-rotation orbit of the vacuum.** Localization is not
an operator identity, but it holds on far more than the vacuum alone. The sharp statement (Paul,
notebook cells 39-41) is

```
<0| Exp[r EE[{2,1}]] . t[{1,mu}][w1]  ==  <0| Exp[r EE[{2,1}]] . (t1[{1,mu}][w1] (x) 1)
```

for all `r` and all `mu`, where `EE[{i,j}] = Ee[1][{i,j}] (x) 1 + 1 (x) Ee[2][{i,j}]` is the
**global** (coproduct) generator. So on the whole orbit, the two-site transfer matrix at `w1` acts
as the **one-site** transfer matrix on site 1. Since `EE` is a sum of commuting terms,
`Exp[r EE[{2,1}]] = Exp[r Ee[1][{2,1}]] (x) Exp[r Ee[2][{2,1}]]`, the orbit consists of product
states and the factorized form survives — site 2 now sits in the *rotated* vacuum rather than the
bare one. The orbit is genuinely one-parameter: `<0| Ee[{1,2}] = 0`, so the raising direction acts
trivially and the Cartan only rescales; `Exp[r EE[{2,1}]]` already sweeps the full group orbit.

**Globality is essential.** Replacing the rotation by independent parameters,
`Exp[r1 Ee[1][{2,1}]] (x) Exp[r2 Ee[2][{2,1}]]`, breaks the identity; it is restored exactly at
`r2 = r1`. So this is the diagonal su(2) action — the one under which the tensor product
decomposes into irreps — and not a product of independent local rotations.

The mirror holds too (verified in scratch, not in the notebook): with the maximal vacuum and the
raising direction,
`xb[{0,0}].Exp[r EE[{1,2}]].tb[{1,nu}][w2 + h Lam[2]_1] == xb[{0,0}].Exp[r EE[{1,2}]].(1 (x) tb2[{1,nu}][...])`.

**Earlier misreading, recorded so it is not repeated.** Probing with a generic
`{{a1,a2,a3,a4}} (x) vac2` — site 1 varied, site 2 held at the bare vacuum — leaks, with leakage
proportional to `a1, a2, a3`, vanishing only on `<0|`. That is true but it is *not* a probe of the
orbit, because the orbit rotates **both** sites by the same parameter. It shows only that the
subspace `(anything) (x) vac2` is not preserved; it says nothing about the orbit, which is.

The leakage coefficients carry denominators `3h - w1 + w2 = h Lam[1]_1 - (w1 - w2)`. That pole sits
exactly where the `B` spectrum collides (`w1 + h mu1 = w2 + h mu2` with `mu` in range forces
`w1 - w2 = h (mu2 - mu1)`), so the SoV basis degenerates there anyway — a consistency check rather
than a separate restriction.

### Checked in scratch, not recorded in the notebook

- `T[{1,2}][w_n]` annihilates `<0bar|` (mirroring `T[{2,1}]` annihilating `<0|`).
- `T[{2,1}][w[1]]` on `<0bar|` does lower `mu1` in its leading component but is **not** a `B`
  eigenvector. Bare monodromy entries are not the lowering operators.
- The unchanged `t` tower evaluated on `<0bar|` at any of the obvious points is not an eigenvector
  either.

---

## 5. Traps

**The quantum-determinant trap (cost most of the search).** In the mirror tower, `tb[{2,1}]` keeps
its own `-h` internal shift — it stays the genuine quantum determinant
`T11(v)T22(v-h) - T21(v)T12(v-h)`. Mirroring `h -> -h` there as well gives
`T11(v)T22(v+h) - T21(v)T12(v+h)`, which is **not central**, and the fusion then fails. It fails
only at `s >= 2`: `s = 0` and `s = 1` still look perfect, so a shallow test passes and the bug hides.
Whenever a mirrored/dual fusion tower is being built, check `s >= 2` explicitly.

**Poles at the evaluation point.** `nu[1]` vanishes exactly at `w_n + h Lam[n]_1`, which is where
the mirror tower must be evaluated. This is why the mirror normalizes by `nu[2]`; an attempt to
reuse `Tt` there hits `0/0`.

**Notebook hygiene** (from the repo-level `CLAUDE.md`, and live here):
- Cell 3's `ClearAll[Ee, LL, nu]` also wipes `nu[2]`, which is defined much further down. A full
  top-to-bottom run is fine; re-running cell 3 alone is not.
- `Ttb`, `Tt2b`, `xb` are memoized (`f[x_] := f[x] = ...`). After editing a definition, re-run the
  defining cell **and** everything downstream that used the stale cache.
- `x` and `xb` each have an explicit `[{0,0}]` rule plus a general pattern rule; the `ClearAll` sits
  in the earlier cell, so re-running that cell alone drops the pattern rule.

---

## 6. What generalizes (Paul, confirmed — not re-derived here)

- Unequal reps `Lam[1] != Lam[2]` behave identically. The equal-rep choice is convenience, not a
  special case.
- There is a further mirror swapping left and right `B`-eigenvectors, with the opposite site
  localized in each case. Four variants in total; two are implemented in the notebook.
- The two variants in the notebook between them factorize **7 of the 16 states**, not all of
  them: the `mu2 = 0` edge (4 states, construction A) and the `mu1 = Lam[1]_1` edge (4 states,
  construction B), overlapping at `x[{3,0}]`. The other two edges do **not** factorize — cell 16
  checks `x[{0,1}]`, on the `mu1 = 0` edge, and finds no rank-one form. The 9 remaining states are
  genuinely entangled across the two sites. The left/right mirror gives the corresponding
  statements for the *right* eigenbasis, which is a different set of objects, not more factorizing
  left states.
- So the usable statement is: one *edge* of the SoV grid can always be split off as a bare
  single-site object, in either of two ways. See the `project_cg_overlap_dictionary` and
  `reference_fsov_determinant_normalization` session memories.

## 7. Open

### What localization is for

Sketched by Paul on 2026-08-27, **to be set up explicitly in a later session**. Recorded here so it
is not re-derived or guessed at; none of it is implemented or verified yet.

1. There are globally rotated states `<Omega(r)|` for which the overlap `<Omega(r) | Psi>` is known
   **exactly in terms of Q-functions**.
2. Localization (section 4) says that acting on such a state with `t[{1,mu}][w1]` produces something
   that acts on **site 1 only**.
3. There are **combinations of the `t`'s that sum up to global rotations**. Combined with (2), such a
   combination acting on `<Omega(r)|` rotates only site 1 — producing `<Omega(r'), Omega(r)|`, a
   product state with **independent** rotation parameters on the two sites.
4. `Psi` is a transfer-matrix eigenvector, `t|Psi> = tau|Psi>`, so the same overlap also equals
   `tau` times the already-known `<Omega(r) | Psi>`.
5. Therefore `<Omega(r'), Omega(r) | Psi>` is known in terms of Q-functions — a two-parameter family
   of overlaps against independently rotated per-site states, which is the object CG coefficients
   are read off from.

**Why localization is the crux.** The input family is **one**-parameter (a global rotation, the
diagonal su(2)); the target `<Omega(r'), Omega(r)|` is **two**-parameter. Localization is precisely
what converts a global object into a one-site one, breaking the diagonal constraint and buying the
extra parameter. This is why the "globality is essential" check in section 4 matters in both
directions: the orbit must be global for localization to hold, and the payoff is that acting on it
with a transfer matrix takes you off that orbit.

Open question for that session: only an edge of the SoV grid factorizes (7 of 16 states) — is that
edge all the construction ever needs, or is the entanglement of the interior 9 the obstruction the
local-rotations thread is working around?
- The `(S-mu)!/S!` one-site normalization, checked only at `S=3`.
- The left/right mirror (variant C) is not implemented in the notebook.
