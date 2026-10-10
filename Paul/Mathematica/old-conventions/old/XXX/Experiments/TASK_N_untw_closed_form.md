# TASK: closed form for `N_untw` — status and briefing for the next agent

**Read this file completely before touching anything.** Updated 2026-09-11 (session 4).

Work lives in `Experiments/SU4_SoV_PerturbativeOverlap_L2.wb`, now **68 cells**, explicitly saved.
Session 4 is in cells 44-68. The main result and its proof are in
**[N_untw_balanced_closed_form.md](N_untw_balanced_closed_form.md)**; read that and Section 9 first.
Session 3 is in cells 24-43. Cell 19's incorrect SOLVED heading was corrected, and cell 22 now
includes the previously missing `Wuntw` definition.

---

## 0. Bottom line — what is and is not done

**Asked for:** `N_untw` in closed form in terms of untwisted (regular at `phi = 0`) Q-functions.

**Session-4 result:** for the balanced family `mu={L-m,L-m,m,m}`, `L=S+P`, `d=2m`, there is now
a **derived 8x8 determinant entirely in untwisted Q-functions**, implemented by `balancedQMat`
and `balancedQNorm` (cell 60). The product simplifications at **d=2 and d=4 are proved for
general spins**, not fitted. A common all-m product is independently validated at d=6 and d=8,
but its all-m factorisation proof remains open. **Unbalanced irreps are not covered.**

**Answered in full:** the previously tested equal-weight cases (Section 2.1), and now the vacuum
for **arbitrary positive integer S,P**, with a proof and explicit product (Section 8.1).

**New exact Q-function result:** for `mu={L,L-1,1,0}`, `L=S+P`, the identity
`N/Nvac = L [Qbar12(theta1+S h) Qbar12(theta2+P h)]^3` is established by exact analytic
determinants for six spin pairs with symbolic separation, including `(2,2)` and `(3,2)`.
It remains a conjecture for arbitrary spins. The exact twisted Q-functions themselves are proved
for general S,P. See Sections 8.2-8.3.

**Historical empirical table:** Section 1 describes six previously measured states, not all mu.
The balanced d=2 entries are now superseded by a general-spin derivation in Section 9. The old
degree-only inhomogeneity pattern does **not** extend unchanged to d=4; the exact factor is
`A[-1]^3 A[0]^4 A[1]`, not `A[0]^3 A[1] A[2] A[3]`, with `A[k]=(S+k+t)(P+k-t)`.

**Newly ruled out (Section 4):** the whole class of first-order reduced-minor / Appendix-N (N4)
formulas, because `corank < vanishing order` for every degenerate state.

**Original lead (Section 5; advanced in Section 8):** for all six states studied the monic branch polynomials have degree
`<= 1`, so the twisted Baxter system is *exactly* solvable in the twists. That gives `det M` as an
explicit function of `z_1..z_4` and turns the confluent limit into an honest Taylor extraction.
The d=1 family is now solved exactly. Eliminate its quadratic root by an on-shell rational
parametrisation **before** taking the determinant; keeping three symbolic variables was slow.

---

## 1. The empirical closed form

Write `t = (theta1 - theta2)/h` and `d = M_1 + M_2 = 2(S+P) - mu_1 - mu_2 = deg Q_12`
(degree of the untwisted Casoratian). Normalising by the vacuum `mu = {S+P, S+P, 0, 0}`:

```
N_untw(mu)                                3    d-1
----------------  =  c(S,P,d) [(S+t)(P-t)]   PROD (S+k+t)(P+k-t)
N_untw(vacuum)                                k=1
```

The first cross-ladder difference pair enters **cubed**, every later one linearly. Here
`S+k+t = (theta1 + (S+k)h - theta2)/h` and `P+k-t = (theta2 + (P+k)h - theta1)/h`.

### 1.1 The verified data

Every ratio is **exactly real** and **exactly rational**. Method: Richardson-extrapolate
`Det[spMat]/eps^ord` to `eps = 0` over `eps = 1e-6 .. 1e-12` at 500-600 digits (cell 23), recognise
the 30-digit limit as a rational, then **reverify at a second inhomogeneity** `theta2 = 1/5`
(the main runs use `1/7`). Both `theta2` agree to 30-32 digits, so the `t`-dependence is established,
not fitted.

| (S,P) | mu | d | ratio to vacuum (exact, theta2 = 1/7) | closed form |
|---|---|---|---|---|
| (1,1) | {2,2,0,0} | 0 | 1 | 1 |
| (1,1) | {2,1,1,0} | 1 | 76765625/2744515872 | `(1/32)(1+t)^3(1-t)^3` |
| (1,1) | {1,1,1,1} | 2 | 1073490500000/1021217202747 | `(8/27)(1+t)^3(1-t)^3(2+t)(2-t)` |
| (2,1) | {3,3,0,0} | 0 | 1 | 1 |
| (2,1) | {3,2,1,0} | 1 | 3825694144/20841167403 | `(8/243)(2+t)^3(1-t)^3` |
| (2,1) | {2,2,1,1} | 2 | 2435054322656/63038098935 | `(6/5)(2+t)^3(1-t)^3(3+t)(2-t)` |

At `theta2 = 1/5` (`t = 2/15`) the same formulas hold; e.g. `r({2,2,1,1}) = 189481091072/4271484375`.

`c(1,1,0)=1`, `c(1,1,1)=1/32`, `c(1,1,2)=8/27`, `c(2,1,0)=1`, `c(2,1,1)=8/243`, `c(2,1,2)=6/5`.

### 1.2 The twist factor it multiplies (session 1, unchanged, verified on two rays)

```
det M^AA = (phi_1-phi_2)^4 (phi_3-phi_4)^4 PROD_{i<j}(phi_i-phi_j)^(mu_i-mu_j) x N_untw + ...
```

Vanishing order `= Sum_{i<j}(mu_i-mu_j) + 8` in the monic basis, `= Sum_{i<j}(mu_i-mu_j)` in the
regular (divided-difference) basis. **Correction to session 1:** the `(phi_1-phi_2)^4(phi_3-phi_4)^4`
is *not* a monic-normalisation artifact — it is the direct/dual confluence, four rows each.

### 1.3 The partial Q-function rewriting — read the caveat

The monic untwisted Casoratian at the **top of each pole ladder** factorises exactly (verified, cell 22):

```
Qbar_12(theta1 + S h) = kappa_1 PROD_{k=0}^{d-1} (S+k+t)
Qbar_12(theta2 + P h) = kappa_2 PROD_{k=0}^{d-1} (P+k-t)
```

| (S,P) | mu | kappa_1, kappa_2 |
|---|---|---|
| (1,1) | {2,2,0,0} | 1, 1 |
| (1,1) | {2,1,1,0} | 1/2, 1/2 |
| (1,1) | {1,1,1,1} | 1/6, 1/6 |
| (2,1) | {3,3,0,0} | 1, 1 |
| (2,1) | {3,2,1,0} | 2/3, 1/3 |
| (2,1) | {2,2,1,1} | 3/10, 1/10 |

so the ratio can be written

```
N_untw(mu)/N_untw(vac) = const(S,P,d) x [(S+t)(P-t)]^2 Qbar_12(theta1 + S h) Qbar_12(theta2 + P h)
```

**But this is a repackaging, not a formula.** Two things go wrong:

1. `[(S+t)(P-t)]^2` is **not** a Q-value except accidentally at `d = 1` (where `Qbar_12` happens to
   be linear).
2. `const = c/(kappa_1 kappa_2)` is **not constant**: it is `1/8, 32/3, 4/27, 40` for
   `(1,1,d=1), (1,1,d=2), (2,1,d=1), (2,1,d=2)`. No natural renormalisation of `Q_12` fixes this —
   using the un-normalised Casoratian `W = -(mu_1-mu_2+1) Qbar_12` gives `1/32, 32/3, 1/27, 40`,
   which is worse.

Useful fact regardless: **`Q^34(u) = (13/5) Q_12(u - h)` exactly**, on every state checked. The dual
Casoratian is the same function shifted, so the "four Q's" of the original guess collapse to two.

---

## 2. What IS `N_untw` in untwisted Q-functions (the parts that are done)

The *untwisted* 8x8 — `spMatU` in cell 20 — is built from the
untwisted Baxter polynomials with the same measure/residue functional as `spMat`. Rows `(alpha,i,a)`
use the direct `P_i` from `UntwSolve` and the dual `Pd_a` from the new `UntwSolveD`.
Its determinant equals `N_untw` in the tested equal-weight cases, but vanishes for the other
states. It is not a formula for their nonzero `N_untw`.

### 2.1 `mu = {m,m,m,m}`: exact, verified

There `Sum_{i<j}(mu_i-mu_j) = 0`, `detU` is nonzero, and **`detU = N_untw` exactly**. For
`{1,1,1,1}` at `(1,1)`:

```
detU = 1425.34685148326739971849627658529 - 439.66127673479367585438512664067 I
```

and `detU/N_untw_measured - 1` scales **exactly linearly in eps** (`0.0078729`, `0.00078730`,
`0.000078730` at `eps = 1e-3, 1e-4, 1e-5`), i.e. the residual is entirely the *measurement's* own
`O(eps)` error. `detU[2,2,{2,2,2,2}] = 1223.86142488821645769646867390055 - 377.51125352602203800640858352704 I`
is computed but its `N_untw` has not been measured for comparison (see Section 6).

Note the row-gauge freedom is harmless: adding a lower-degree `P` to a higher one is a simultaneous
row operation on `(alpha, i, a)` for all `alpha, a`, so it cannot change the determinant. Only the
monic normalisation matters.

### 2.2 The vacua `mu = {S+P,S+P,0,0}`: exact, closed form

There `M_1 = M_2 = 0` and `mu_3 = mu_4 = 0`, so the monic branch polynomials are `q = 1` **exactly**
on both sides and every Q is a **pure exponential**. The whole 8x8 is then an explicit analytic
function of the twist (`spMatE`, cell 21) — no perturbative solve, no truncation:

```
N_untw(vacuum, (1,1)) = 1355.94001494748764131246536678 - 418.252103007227979590132389910 I
```

with vanishing order 16. Reading it off at `eps = 1e-30` with 400 digits agrees with the Richardson
extrapolation of the perturbative determinant to 29 digits — a strong cross-check of both.

Explicitly, with `x_{ia} = phi_i - phi_{j(a)}` (`j(1)=4, j(2)=3`), `w_a = e^{-phi_{j(a)} h}`,
`v_i = e^{-phi_i h}`, `sigma = s_a/h in {0,1}`, `tau = -s_b/h in {0,1}`:

```
M[(alpha,i,a),(sigma,tau,beta)] = w_a^sigma v_i^tau G_alpha(x_{ia}, beta),
G_alpha(x,beta) = Sum_p perV(alpha,p) w(p) p^beta e^{x p}
```

This is fully explicit and is the natural starting point for an analytic confluent limit.

---

## 3. Machinery (all in the notebook)

Cells 1-18 (session 1): `SetState`, `Bax`/`BaxD`, `UntwSolve`, `pushQ12`/`pushQ34`, `dirQ`/`dulQ`,
`spMat`, `ordNumD`, `leadC`, `L1fac`, `cand`. Conventions (twist ray, monic `q_k`, dual indexing
`a=1,2 <-> j=4,3`, the measure `E1`, the 8x8 row/column ordering) are in the notebook markdown and
in `SU4_PerturbativeExpansion_L2.md`. **Kernel restart: run cells 1-15, then 19-22.**

New this session:

| symbol | cell | what |
|---|---|---|
| `UntwSolveD[S,P,mu]` | 20 | untwisted **dual** solver, mirror of `UntwSolve`; staggered degrees `MhatD = mu[[{4,3,2,1}]] + {0,1,2,3}`; returns `"Pd"` |
| `spMatU[S,P,mu]`, `detU` | 20 | the untwisted 8x8 and its determinant |
| `spMatE[S,P]` | 21 | the exact pure-exponential vacuum matrix |
| `closed`, `Cst`, `q12mon`, `qfac`, `meas` | 22 | the empirical closed form + verification table |
| `detAt`, `richard`, `NuntwExact` | 23 | high-precision extraction (**not evaluated**, ~20 min) |

`Wuntw[S,P,mu]` (Casoratian of `P_1,P_2`) was also added in session 2 exploration; it is defined in
cell 22's dependencies via `q12mon`.

---

## 4. Dead ends — do NOT repeat

**Session-4 qualification:** the searches below excluded the scanned product ansatz, not
quotients of Q-values. The new d=2 and d=4 formulas use
`[Qbar12(top1) Qbar12(top2)]^3/[P1(top1) P1(top2)]^2`. Also, the successful new determinant
first divides out analytic cross-pair zeros; it is not a first-order reduced minor of spMatU.

1. **Pure products of Q-values.** Uniform shifts (session 1), *independent* shifts per inhomogeneity
   (all four scanned over `{0,+-1/2,+-1,+-3/2,+-2}`), both sign conventions of the `I^nu`
   bookkeeping, multisets of up to 8 factors, against the exact rational targets: **none** works. The
   scanned pure-product ansatz does not work across those states. This does not exclude Q-function
   quotients, which now give the balanced d=2 and d=4 answers.
2. **The literal (N4) reduced minor** — drop the four `alpha=1` rows and the `beta=1` columns, keep
   the 4x4 `alpha=2` minor of the untwisted matrix — gives *complex* ratios. Wrong.
3. **`spMatX`**: staggered untwisted `P_k` dressed with one exponential `e^{phi_i u}` per branch. It
   gets the vanishing order right (`8, 6, 0` = `Sum_{i<j}(mu_i-mu_j)`) but the coefficient is
   **ray-dependent** (ray B / ray A = `-1.35` and `7.20` where it must be `1`). Reason: the physical
   regular basis for branch 2 is a *divided difference* mixing both branches, which is not of the
   form `e^{phi u} x polynomial`.
4. **`spMatB`**: the opposite fix — same leading polynomial `P_1` on both rows, divided-difference
   exponential. Correct for the vacua (trivially, `P_1 = 1`) but off by `1/0.00244` for `{1,1,1,1}`
   and `1/0.0782` for `{2,2,1,1}`.
5. **Any first-order reduced-minor / degenerate-perturbation formula.** Measured coranks of
   `spMatU` versus the vanishing order:

   | mu | (S,P) | vanishing order | corank of `spMatU` |
   |---|---|---|---|
   | {2,2,0,0} | (1,1) | 8 | 3 |
   | {2,1,1,0} | (1,1) | 6 | 2 |
   | {1,1,1,1} | (1,1) | 0 | 0 |
   | {3,3,0,0} | (2,1) | 12 | 4 |
   | {3,2,1,0} | (2,1) | 10 | 3 |
   | {2,2,1,1} | (2,1) | 4 | 2 |

   `corank < order` always (and the order is not even a function of the corank: corank 3 gives both 8
   and 10; corank 2 gives both 6 and 4). So the degeneracy is **not** lifted at first order — the
   kernel directions have Jordan-type structure and the leading term needs the matrix expanded to
   order `~ order - corank + 1` in those directions. A single reduced minor cannot express it.
6. **Reading `N_untw` at `eps = 1e-4`** gives only ~8 digits of the ratio — not enough to recognise a
   rational with denominator `~1e12`. Always Richardson-extrapolate (cell 23).

---

## 5. The best open lead: exact solution in the twists

**Session-3 update:** steps 1-3 below are complete for d=1 at the tested spin pairs, and step 4
is proved for general spins. The exact solution is algebraic, not rational in generic twists.
Section 8 records the formulas and the much faster on-shell parametrisation.

**The observation.** `M_k = (S+P) - mu_k`, and the monic branch polynomial `q_k` has degree `M_k`.
For all six states studied the two *used* branches have degree `<= 1`:

| mu | (S,P) | M | direct `q_1, q_2` | dual `q^4, q^3` (degrees `mu_4, mu_3`) |
|---|---|---|---|---|
| {2,2,0,0} | (1,1) | {0,0,2,2} | `1`, `1` (exact) | `1`, `1` (exact) |
| {2,1,1,0} | (1,1) | {0,1,1,2} | `1`, `u+b` | `1`, `u+c` |
| {1,1,1,1} | (1,1) | {1,1,1,1} | `u+b1`, `u+b2` | `u+c1`, `u+c2` |
| {3,3,0,0} | (2,1) | {0,0,3,3} | `1`, `1` (exact) | `1`, `1` (exact) |
| {3,2,1,0} | (2,1) | {0,1,2,3} | `1`, `u+b` | `1`, `u+c` |
| {2,2,1,1} | (2,1) | {1,1,2,2} | `u+b1`, `u+b2` | `u+c1`, `u+c2` |

So the entire twisted problem has **at most four scalar unknowns** plus the five `tc`. That is small
enough to solve **exactly in `z_1..z_4`**, which turns the confluent limit from a numerical
extrapolation into an analytic Taylor extraction — exactly what is needed to *derive* Section 1 and
to see what untwisted object it is.

**Already verified to work.** For `mu = {2,1,1,0}` at `(1,1)`, branch 1 is `q_1 = 1` exactly, and

```
Solve[(# == 0 &) /@ DeleteCases[CoefficientList[Expand[(BaxRes[1,1,{2,1,1,0}][1][Function[y,1]][x] /. prm)], x], 0], tcvars]
```

returns a **unique** solution: 7 coefficient equations fix `tc[3], tc[4], tc[5]` as explicit rational
functions of `tc[1], tc[2]` and `z_1..z_4`, leaving `tc[1], tc[2]` free. Branch 2 (`q_2 = u + b`,
one new unknown) should then close the system. `Solve::svars` is expected and harmless here.

**Suggested route.**
1. Finish the exact solve: branch 1 (`q=1`) plus branch 2 (`q = u+b`) for `tc[1], tc[2], b` in terms
   of the `z_j`. Do the same on the dual side with `BaxD` for `c`. Cross-check against
   `pushQ12`/`pushQ34` by expanding on the ray — they must agree order by order.
2. Build `spMat` with these exact `q`'s. Now `det M` is an explicit function of `z_1..z_4`, exact.
3. Take the confluent limit analytically: divide by
   `(phi_1-phi_2)^4 (phi_3-phi_4)^4 PROD_{i<j}(phi_i-phi_j)^{mu_i-mu_j}` and set `phi -> 0`. The
   result must be the Section 1 rational; the point is to see *which untwisted determinant* it is.
4. Start with the vacuum, where step 1 is empty (`q = 1` on both sides) and the matrix is already
   explicit (Section 2.2). Getting `N_untw(vacuum)` in analytic closed form for general `(S,P)` would
   already settle much of the `c(S,P,d)` question.

**Caution.** Whether the *exact-in-`z`* solution is even necessary is unclear: it may be enough to
expand `q_k(u; phi)` to the order the leading term needs. But do not assume the Taylor coefficients
of the twisted family are the staggered `P_k` — they are not. For `mu = {2,1,1,0}` the `eps^1`
coefficient of `Q_1` has degree `M_1 + 1 = 1`, while the untwisted solution degrees are
`Mhat = {0,2,3,5}`; there is no untwisted solution of degree 1. The reason is that `tau` is corrected
at the same order, so the `eps`-corrections to `Q` are **not** untwisted Baxter solutions.

---

## 6. Other open items

**Session-3 update:** `c(d=1)=(SP)^3/(S+P)^5` is now established exactly for six spin pairs
with symbolic t; its arbitrary-spin proof is open. The old P=1 fit below is superseded. Also,
degree d alone need not uniquely label an irrep for larger spins: keep mu in the labels.

1. **`c(S,P,d)` and `kappa` are tabulated, not derived** — fitted from `(1,1)` and `(2,1)` only, both
   with `P = 1`. `c(d=1)` is consistent with `S^3/(S+P)^5`, but `S = 1` hides the numerator, and
   `c(d=2)` (`8/27`, `6/5`) has no fit. `kappa` at `(1,1)` is `1, 1/2, 1/6` which looks like
   `1/(d+1)!`; at `(2,1)` it is `1,1`, `2/3,1/3`, `3/10,1/10`.
   **`(1,2)` gives nothing new** — the problem is symmetric under `theta1 <-> theta2, S <-> P,
   t -> -t`. You need `(2,2)` or `(3,1)`.
2. **Only `d = 0, 1, 2` tested.** `{2,2,2,2}` at `(2,2)` has `d = 4` and is the cheapest test of the
   `PROD_{k=1}^{d-1}` tail; its `detU` is already computed (Section 2.1), so it only needs its
   `N_untw` measured (order 8) to compare. `{4,3,1,0}` at `(2,2)` has `d = 1`.
3. **The exponent pattern `3, 1, 1, ...`.** The `3` is plausibly `N-1` for `gl(N)` but that is a
   guess; `gl(3)` would test it.

---

## 7. Traps

Full session-1 list in `SU4_PerturbativeExpansion_L2.md` (back-substitution of accumulated solution
rules; `Solve` genericity; never name a pattern variable `N`; read determinant orders numerically;
distinguish real vanishing from truncation by varying `Ntr`; `pushQ12/pushQ34` read `S,P,mu` through
`Bax`'s globals; `psi` is global; one expression per line in the front end; `abort` often fails on a
long `Solve`). Added this session:

- **`prm` is baked into memoised results.** `UntwSolve`, `UntwSolveD` and `tau3quad` memoise on their
  arguments only, so changing `theta1/theta2/h` silently returns stale answers. Strip the cached
  `DownValues` first:
  `DownValues[f] = Select[DownValues[f], ! FreeQ[#[[1]], Pattern] &]`
  for each of the three, and re-issue the `chi[a_] := chi[a] = ...` definition (it memoises too).
  Restore `prm` afterwards — the notebook's documented state is `{theta1 -> 1/3, theta2 -> 1/7, h -> 1}`.
- **Working precision.** With a leading order of 12-20 the determinant is `~1e-120` at `eps = 1e-10`;
  use 500-600 digits in `detAt`, and `Block[{$MaxExtraPrecision = 2000}, Quiet[...]]` for `spMatE`.
- **`I^nu` bookkeeping.** `leadC` divides by `cand[psi,mu]` which uses `psi`, not `phi = I eps psi`.
  So the tabulated `N_untw` equals `I^nu` times the `phi`-form, `nu = Sum_{i<j}(mu_i-mu_j) + 8`.
  Ratios to the vacuum must be taken as `(N_a/I^{nu_a})/(N_vac/I^{nu_vac})`; otherwise `{2,1,1,0}`
  comes out with a spurious minus sign. All ratios in Section 1 are in the corrected convention and
  are exactly real (imaginary parts vanish to 60 digits).
- **Ratios converge far better than absolute values.** The universal phase (`Arg = -17.0977 deg` at
  `(1,1)`, from the measure and the `Exp[2 Pi I theta_alpha]` factors of the `alpha=2` rows) and most
  of the truncation error cancel.

---

## 8. Session 3: exact vacuum product and d=1 Q-function identity

### 8.1 Vacuum: proved for arbitrary positive integer S,P

Let `L=S+P`, `t=(theta1-theta2)/h`, and `e_alpha=Exp[2 Pi I theta_alpha/h]`. In the phi convention,

```
Nvac = (e2-e1)^4 h^(4-4L) / [S! P! Product[t+k,{k,-P,S}]]^4.
```

This is a proof, not an extrapolation. The notebook's cell 38 contains the full argument:

1. Separate the two ladder functionals. The change of basis contributes `(e2-e1)^4`.
2. Remove the row exponentials `Exp[(phi_i-phi_j) theta_alpha]`, whose product tends to 1.
   Multiply row `(alpha,i,j)` by `z_i z_j^(S_alpha+1)`, where `z_i=Exp[h phi_i]`.
   The resulting determinant is a homogeneous polynomial of degree `4L+8` in the four z's.
3. Coincident direct/dual branches give `(z1-z2)^4(z3-z4)^4`. For each cross pair i=1,2 and
   j=3,4, recombine the ladder rows. The full residue functional kills polynomials of degree
   <=L, so every entry has a zero of order >=L at `z_i=z_j`. Thus the determinant contains
   `Delta(z)=(z1-z2)^4(z3-z4)^4 Product[(z_i-z_j)^L,{i,1,2},{j,3,4}]`.
   These factors exhaust its degree.
4. The coefficient of `z1^(2L+4) z2^(2L)` uses only the two top poles and factorises into two
   4x4 endpoint determinants. It fixes the remaining constant to
   `K=[(p2top-p1top) w(p1top) w(p2top)]^4`
   `=h^(-8L-4)/[S! P! Product[t+k,{k,-P,S}]]^4`.
5. Divide out the row multipliers, which tend to 1, and use `z_i-z_j ~ h(phi_i-phi_j)`.

In fact the finite-twist ladder determinant is exactly

```
K Delta(z) / [(z1 z2)^4 (z3 z4)^(2L+4)].
```

Checks: symbolic t at `(1,1),(2,1),(2,2),(3,1)`; arbitrary rational nonconfluent twists;
and `h=2/3`. All exact identity tests pass. The `(1,1)` limit matches the old 400-digit
finite-eps evaluation to 29 digits, as expected from its O(eps) error.

### 8.2 Exact twisted d=1 family: proved for general S,P

For `mu={L,L-1,1,0}`, introduce the root r of

```
z2 (r-theta1-S h)(r-theta2-P h) = z3 (r-theta1)(r-theta2)
```

that is regular at confluence, with

```
r0 = (SP h + P theta1 + S theta2)/L.
q1(u) = 1
q2(u) = u-r + h z1/(z1-z2)
q^4(u) = 1
q^3(u) = u-r + h z4/(z3-z4).
```

The five transfer coefficients are in `exTc` (cell 32). All four Baxter residuals vanish
identically modulo the quadratic equation, with symbolic S,P,h,theta1,theta2,z1..z4.
The other root diverges as z2 approaches z3; it does not continue the selected untwisted irrep.
Both branch Laurent series match the independently saved `rA12`/`rAd12` through eps^8 exactly.

Regularising the direct pair gives `P1=1`, `P2=u^2+(h-2r0)u` up to a constant, hence the monic
Casoratian is **`Qbar12(u)=u-r0`** for this family. In particular,

```
Qbar12(theta1+S h) = h S(S+t)/L
Qbar12(theta2+P h) = h P(P-t)/L.
```

### 8.3 Norm identity: exact for six spin pairs; general-spin proof still open

Exact determinant extraction, with symbolic t rather than rational recognition, gives

```
N_d=1/Nvac = (SP)^3 h^6/L^5 [(S+t)(P-t)]^3
          = L [Qbar12(theta1+S h) Qbar12(theta2+P h)]^3.
```

The spin pairs proved by these exact symbolic computations are:

| (S,P) | coefficient at h=1 |
|---|---|
| (1,1) | 1/32 |
| (2,1) | 8/243 |
| (2,2) | 1/16 |
| (3,1) | 27/1024 |
| (3,2) | 216/3125 |
| (4,1) | 64/3125 |

The h dependence follows by rescaling u and the monic Q polynomials; calculations use h=1.
An independent rational ray at `(2,2)` gives the same coefficient exactly. This d=1 rewriting
uses the linearity of Qbar12 and does not solve d>=2 or the general-irrep problem.

**Do not describe the all-spin d=1 identity as proved.** Its twisted Q solution is proved for
general spins; its norm has been established for the six listed pairs, at symbolic separation.

### 8.4 Efficient machinery and next step

| cells | symbols / purpose |
|---|---|
| 26,31-33 | `exGeneralSol`, `exRootEq`, `exRed`, `exTc`, `exDualSol`: exact d=1 Baxter solution |
| 28-30 | `exVacMat`, `exVacLead`, symbolic vacuum determinants |
| 35-37 | `exD1Mat`, `exD1Det`, `exRootJet`, `exD1Lead`: exact rational-twist determinant and root jet |
| 38 | arbitrary-spin vacuum proof |
| 39-40 | `exD1OnShell[S,P,True]`: exact norm coefficient without the universal phase, with symbolic t; six symbolic identities |
| 41-42 | independent perturbative-series, ray, and finite-twist identity checks; `exVacFormula[S,P]` |
| 43 | Q-function result and scope |

The crucial speedup: eliminate the quadratic root **before** `Det`. At h=1 set

```
r=r0+eps;
delta=-L eps/[(r-theta1)(r-theta2)];
z=1+delta {3,0,1,-3}.
```

Then the Bethe equation is an identity and all matrix entries are rational in eps and t.
The six symbolic-separation checks take seconds. A determinant retaining eps, r and t as
three independent variables exceeded 40 seconds for each pair; do not repeat that route.
`exD1Mat`, `exD1Lead`, `exD1OnShell` implement h=1; do not use them at other h without adapting
the polynomial shifts. `exVacMat` and the general vacuum proof support other nonzero h.

For a cold kernel, load the existing dependencies as before, then run the session-3 cells in
order. `Wuntw` is now defined in cell 22, so that missing dependency is fixed. Cell 23 is not
needed. The session-3 code does not require any of the old unsaved exploratory definitions.

**Next research:** prove the d=1 norm formula for arbitrary spins from the exact on-shell
determinant, and extend the exact method to the degree-one-on-both-sides d=2 family
`mu={L-1,L-1,1,1}`. Beyond P=1, do not assume d alone specifies an irrep or its normalisation.

---

## 9. Session 4: balanced d=2,4 proved; all even degrees have an untwisted determinant

Full derivation, assumptions and formulas:
**[N_untw_balanced_closed_form.md](N_untw_balanced_closed_form.md)**.

Let `mu={L-m,L-m,m,m}`, `L=S+P`, `d=2m`, `n=L-2m`, and use the Q-functions of that state:

```
A = P1(theta1+S h) P1(theta2+P h)
B = Qbar12(theta1+S h) Qbar12(theta2+P h).
```

The two **proved general-spin** product identities are

```
d=2: N/Nvac = 8 (L+1)^3 (L+2) B^3/A^2
d=4: N/Nvac = 432 (L-1)^3 L^4 (L+1) B^3/A^2.
```

The proposed common product factorisation is

```
N_2m/Nvac = m! (m+1)!^3 Pochhammer[n+3,m]^3 Pochhammer[n+4,m] B^3/A^2.
```

It is proved at m=0,1,2, and independently validated at m=3,4 by exact arithmetic:
`(S,P,m)=(3,3,3),(4,3,3),(4,4,4)`, including a second inhomogeneity for `(4,3,3)`.
**The all-m product-factorisation proof is still open.** This distinction is recorded in the
candidate ledger (cell 67).

### 9.1 The all-balanced-degree determinant is derived

`balancedQMat[S,P,m]` (cell 60) uses only the first two direct and first two dual untwisted
polynomials, finite residue sums, and factorials. With `e_alpha=Exp[2 Pi I theta_alpha/h]`,

```
N_untw = (e2-e1)^4 Det[balancedQMat[S,P,m]].
```

The crucial observation is `deg(q_i q^j u^beta)<=2m+1`. The total residue functional kills
all degrees <=L, so the four total-functional rows have analytic factors
`(phi_i-phi_j)^(L-2m)`. Divide these out first. The remaining direct/dual divided differences
only need the first branch differences, which are encoded by the untwisted basis:

```
R = P2 - (u-a) P1,                  a = Coefficient[P1,u,m-1]
Rd = -Pd2 + (u-ad) Pd1,             ad = Coefficient[Pd1,u,m-1].
```

These have degree <=m-1. Taking the four coefficients `1,w,v,vw` of the resulting analytic
row-generating polynomials gives the 8x8 in the linked proof. No unknown twisted Taylor
coefficients remain. For m=0 set a=ad=0.

Independent checks: d=2 agrees at jet orders 1,2,3,4,6,8 and on a second ray; d=4 `(3,2)` agrees
exactly with an independently solved order-4 twisted Baxter jet and determinant order 12;
d=4 `(2,2)` agrees exactly with the equal-weight `detU`.

### 9.2 Inhomogeneity products and correction to the old guess

With `A[k]=(S+k+t)(P+k-t)`,

```
N_2/Nvac = 8 (SP)^3 (S+1)(P+1) h^8 /[(L+1)^3(L+2)] A[0]^3 A[1]

N_4/Nvac = 432 (S-1)^3 S^4 (S+1) (P-1)^3 P^4 (P+1) h^16
          /[(L-1)^3 L^4 (L+1)] A[-1]^3 A[0]^4 A[1].
```

The d=4 formula is not the old `3,1,1,1` tail. Both expressions were obtained by factorising
the determinant with symbolic S,P,t. The remaining denominator A in the Q formula is exactly
what the earlier searches over products of Q-values missed.

### 9.3 Where the code is, and safe replay

| cells | content |
|---|---|
| 44-47 | initial exact d=2 elimination; the full generic Solve in 47 is time-bounded and was too slow |
| 48-51,54 | `d2DualRes`, `d2FastJet`, `jetSPMat`, `jetDet`, exact d=2 benchmarks |
| 52-58 | symbolic d=2 untwisted seed, first derivatives, confluent matrix and factorisation |
| 59-60 | general balanced-family construction, `balancedQMat`, `balancedQNorm` |
| 61 | `balancedFastJet`; independent d=4 twisted-jet verification |
| 62-63 | `balUntwSeed`, general-spin d=4 direct polynomials and dual residual checks |
| 64 | `balNormalizedMat`, normalised ladder moments, symbolic d=4 determinant factorisation |
| 65-66 | product conjecture, exact d=6/d=8 and second-inhomogeneity checks |
| 67-68 | candidate ledger, rechecked d=2/d=4 identities, final scope |

**Minimal replay from a cold kernel:** run cells **2-5, 7, 12-13, 20, 60, 62, 64-67**.
Check execution evidence after each stage. In particular, cell 20 needs `pairF` from cell 13.
A startup timeout does not mean definitions loaded; do not proceed with dependent cells.
The current notebook is explicitly saved with 68 cells. A VS Code restart between sessions
preserved the earlier 63 cells; the required definitions and final checks were restored/replayed
in the live kernel on 2026-09-11. `prm` remains the documented `{1/3,1/7,1}` configuration.

`balancedQMat` includes general h. The numeric helpers `balEndpoint`, `balVacK`,
`balQCandidate` and the symbolic seed/normalised-moment helpers implement h=1; use the written
scaling formulas for other h. Cell 66 demonstrates scoped memoisation caches when changing prm.

### 9.4 Next work: unbalanced degrees and the all-m product proof

1. Prove the common product factorisation for arbitrary m, starting from the already derived
   finite untwisted determinant. The endpoint factors and factorial-moment representation are
   useful inputs.
2. Treat `mu={L,L-2,2,0}` (another d=2 family) and `mu={L-1,L-2,2,1}` (d=3). The balanced proof
   does not cover their Laurent branches or unequal cross-pair exponents.
3. Do not state that all d>=2 or all irreps are solved: the current general determinant theorem
   concerns `mu={L-m,L-m,m,m}` only.
