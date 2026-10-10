# Perturbative twist expansion — gl(4), L = 2

Reference for `SU4_PerturbativeExpansion_L2.wb` (clean, self-contained, runs cold in a few seconds)
and its scratch predecessor `SU4_PerturbativeBaxter_L2.wb` (same results plus the dead ends).

Sites `[S,S,0,0] x [P,P,0,0]`, inhomogeneities `theta1, theta2`, shift `h`. The twist is switched on
along a **unimodular** ray `z_j = Exp[I eps psi_j]` with `Sum_j psi_j = 0` (`psi = {1, 1/3, 1/5, -23/15}`),
so `z_1 z_2 z_3 z_4 = 1` exactly and `phi_j = I eps psi_j`. Test parameters `theta1 = 1/3`,
`theta2 = 1/7`, `h = 1`; everything is exact and rational (times `I`), never numerical.

---

## Result

**Two branches alone determine the whole transfer matrix, with a one-order lag.** Pushing to `eps^N`
using only `Q_1, Q_2` fixes all five `tc` completely at orders `1 .. N-1`; only the top order keeps
free directions, and one more order removes them. Verified to `N = 8` against the exact `tau`.

The **dual** side behaves identically: `Q^4, Q^3` alone fix the same `tau`, same counts, same lag.

| order m | #tc fixed (of 5), direct | dual |
|---|---|---|
| 1–7 | 5 | 5 |
| 8 (= N, top) | 0, free `tcn[1,8], tcn[2,8]` | same |

Both agree with the exact `tau` at orders 1–7 exactly.

## The structure that works

`Q_k = z_k^{u/h} q_k`, with `q_k` **monic of fixed degree** `M_k = (S+P) - mu_k`. The exponential is
kept **exact** and expanded by the operator; it is never re-expanded into the polynomial. Only the
coefficients are expanded, and they are **Laurent**: branch `k` diverges as `eps^{-(M_k - M_1)}`.

Dual side: on the unimodular ray `Prod_{i != j} z_i = 1/z_j`, so `Q^j = Exp[-phi_j u] q^j` with
`deg q^j = mu_j` — the same structure with `mu_j` for `M_k` and the shifts reversed. `BaxD`
annihilates `Exp[-phi_j u]` for `j = 3,4` (and not `j = 1,2`) to order 8, mirroring `Bax` on `j = 1,2`.

The staggered untwisted degrees `Mhat_k = M_k + (k-1)` are the *shadow* of this: expanding
`Exp[phi u]` into polynomials is what raises the degree, so the untwisted `P_k` are limits of
**combinations** of branches, not of individual `Q_k`. Confirmed directly for `k = 2`:
`Q2reg = a Exp[phi_1 u] q_1 + b Exp[phi_2 u] q_2` is regular iff `b = -a` with `a ~ 1/eps` and
nothing worse; the limit is `I(psi_1 - psi_2) av[-1] u + const`, degree `M_1 + 1 = Mhat_2`. It is the
divided difference. `Bax` annihilates `Q1`, `Q2` and `Q2reg` (its `1/eps` piece included).

For `k >= 3` regularity alone does **not** fix the pole depth: with `q_3`'s own coefficients Laurent
the reachable limit degree just follows however deep the poles are allowed to go (`-1 -> 2`,
`-3 -> 3`, `-4 -> 4`), so `Mhat_3 = 4` needs `eps^-4` and nothing internal selects it. A monic `q_3`
also cannot be regular: assuming it is, the limit only reaches degree `2 = M_3`. The `k >= 3` pole
order has to come from the Baxter equation, not from regularity of the combination.

## The independent yardstick: exact `tau`

`tau_a(u) = Sum_{|I|=a} z_I Prod_m nu_{i_m}(u - (a-m) h)` — the pseudovacuum eigenvalue, exact in
the `z`. **The first index carries the largest shift**; the other convention does not match the
ansatz (that cost time once). Matching it to the five-`tc` ansatz succeeds, giving `tc` exactly in
the `z`, and expanding on the ray gives `tau` to all orders. Checks:

- `eps^0` reproduces `UntwSolve`'s `tc0` exactly;
- `eps^1` gives `{124 I/63, 416 I/21, -32384 I/1323, 13312 I/1323, 292 I/63}`, the same values the
  four-branch solve produced;
- `Bax` annihilates `Exp[phi_k u]`, `k = 1,2`, to order 8; `BaxD` annihilates `Exp[-phi_j u]`,
  `j = 3,4`, to order 8.

## API

```
SetState[S, P, mu]              (* sets SS, PP, MU and TC0; Bax/BaxD/push* all read this *)
Bax[n][f]                       (* Baxter operator on f, to order eps^n; expression in the global u *)
BaxD[n][f]                      (* the dual operator, same convention *)
UntwSolve[S, P, mu]             (* untwisted problem: "Mhat", "tc0", "P" (four monic coeff lists) *)
pushQ12[S, P, mu, N]            (* direct branches 1,2 to order N: "M", "report", "sol" *)
pushQ34[S, P, mu, N]            (* dual branches 4,3 to order N: "deg", "report", "sol" *)
qs[k][S, P, mu][n, m][u]        (* seed P_k plus eps^1..eps^n, each an arbitrary degree-m polynomial *)
```

`Bax[n][f]` uses the exactly-known parts of `tau_a` rather than carrying them as unknowns:
leading `chi_a`; subleading `-a (th1 + th2 + (a-1) h) chi_a - h Xi_a` with
`Xi_a = Sum_{|I|=a} z_I Sum_{i in I} n_i = Sum_i n_i z_i d(chi_a)/d(z_i)` and `n = mu`; and
`tau_4 = qdet = chi_4 Prod_i nu_i(u-(4-i)h)` with `chi_4 = 1` on the unimodular ray. Only the five
free `tc[i]` carry unknowns, `tc[i] -> TC0[i] + Sum_m eps^m tcn[i,m]`.

---

## Bugs and traps — the expensive ones

**1. `acc` must be back-substituted. This is the one that cost the most.**
Rules accumulated at earlier orders keep the then-free symbols on their right-hand sides, and `/.`
does **not** re-substitute into an RHS when a later order solves for those symbols. The loop then
silently freezes a wrong branch and a later order reports `INCONSISTENT` when it is not. Fix:

```
acc = Join[acc /. First[sol], First[sol]]        (* NOT Join[acc, First[sol]] *)
```

Symptom that identified it: the *exact* `tau` satisfied the order-3 equations while `Solve` still
returned `{}`. Whenever a per-order solve reports an inconsistency, check the known solution against
those very equations before believing it.

**2. `Solve` genericity.** `Solve[eqns, vars]` treats every symbol outside `vars` as a *parameter*
and returns `{}` whenever the system fails for **generic** parameter values — even if particular
values solve it. Coefficients left undetermined at lower orders must therefore be put **back into
the unknown list**, not carried as loose symbols. Reading the unknowns off the equations
(`Cases[Variables[eqs], _tcn | _qd]`) makes this automatic.

**3. Do not zero leftover `tau` directions.** Leftover *Q-coefficients* are gauge (add `eps^n` times
a lower untwisted solution) and may be set to zero. Leftover `tcn` are **not** — `tau` is shared and
physical, and zeroing it kills the solution three or four orders later.

**4. Do not impose a growing-degree polynomial ansatz.** A Taylor series of polynomials of degree
`Mhat_k + n` is just `Exp[phi u]` written out; imposing it as an *independent* ansatz on top of the
stripped equation over-constrains and dies at `eps^3`. `BaxRes` carries `z[k]^m` prefactors and so
presupposes `Q = z_k^{u/h} q_k`; only the unstripped `Bax` accepts a plain polynomial series.
Concretely the obstruction is `h z_k Prod_{j != k}(z_k - z_j) (mu_k - S - P + deg)` with
`Prod_{j != k}(z_k - z_j) = O(eps^3)`.

**5. Solve order by order, never all orders at once.** The graded system is quadratic across orders
(`tcn[i,1] qc[k,j,1]` enters at order 2), so a single `Solve` is nonlinear, slow, and returns
spurious branches. One order at a time it is linear: minutes become seconds.

**6. `pushQ12`/`pushQ34` read `S, P, mu` through `Bax`'s globals.** Call `SetState` for the *same*
state first. A mismatched call returns an **empty** solution rather than erroring.

**7. Front-end line rule (this notebook's editor).** Keep every expression on **one line, comments
included**. A line break starts a new statement, so a definition, a list, or a multi-line comment
continued on the next line is either silently truncated or reported as `Syntax::sntxi`. Check
bracket balance per line when a cell mysteriously produces nothing.

**8. `Exponent`/`FirstPosition`/`PossibleZeroQ` when reading off eps-orders.** `PossibleZeroQ`
silently misses exact zeros in these coefficients (`Simplify` each first); `FirstPosition` searches
heads and returns `{0}` for a list that itself matches the pattern (use `LengthWhile`); and an entry
order equal to the truncation order is truncation noise — solve two orders deeper than you trust.

**9. Kernel `abort` often does not take on a long `Solve`.** Budget for a restart; this notebook
reloads cold in seconds.

---

## Plugged into the SoV scalar product — orthogonality holds order by order

`SU4_SoV_PerturbativeOverlap_L2.wb` = this notebook plus four cells. Both sides of
`M^{AB}_{(alpha,i,j),(r,beta)} = Contour du E1(u) Exp[2 Pi I u (alpha-1)/h] u^beta [q^j_A q^B_i]_r(u)`
come straight from the perturbative solvers: `pushQ12` gives `q_i` (i = 1,2), `pushQ34` gives `q^j`
(j = 3,4). **The Wronskian construction of the duals is not needed.** `E1` and the periodic factor
are twist-independent, all poles are simple, so the contour integral is a residue sum
`Sum_p perFac(p) g(p) / Prod_{q != p} (p - q)`.

States at `(S,P) = (1,1)`: `A = {2,1,1,0}`, `V = {2,2,0,0}`, both pushed to `eps^12`.

| truncation Ntr | det A\|A (norm) | det V\|V (norm) | det A\|V (cross) | det V\|A (cross) |
|---|---|---|---|---|
| 5 | 14 | 16 | 17 | 17 |
| 6 | 14 | 16 | 18 | 18 |
| 7 | 14 | 16 | 19 | 19 |
| 8 | 14 | 16 | 20 | 20 |
| 10 | 14 | 16 | 22 | 22 |

**Norms sit at a fixed order** (14 for A, 16 for V) whatever the truncation. **Cross terms track
`Ntr + 12` exactly**, i.e. they vanish to every order the expansion resolves and what remains is
pure truncation artifact. That is orthogonality holding order by order in the twist.

Notes:

- The entry orders are additive, `order = row + column`, with `beta = 1` one order below `beta = 0`.
  For `A|A` the row orders are `{3,2,2,1,0,-1,-1,-2}`, so the determinant cannot start below
  `eps^-4` — yet it starts at `eps^14`, i.e. ~18 orders of cancellation. That is why the leading
  blocks looked rank-deficient in the earlier attempt: the naive leading term cancels.
- Read the determinant order **numerically** (log-log slope over `eps = 10^-2, 10^-3, 10^-4`). A
  symbolic `Series` scan would need a truncation deeper than the entries support, and it silently
  reports `Infinity` (everything vanishing) when the true leading order is past the scan window --
  which is exactly what happened at first and looked like a structural singularity. The matrix has
  full rank 8 at finite `eps`; check that before concluding a determinant vanishes.
- Normalisation here (monic `q`, no Wronskian prefactors) differs from the Clean notebook, so the
  absolute orders 14 and 16 are normalisation-dependent; the *comparison* norm vs cross is not.

## Leading behaviour of the norm — the Appendix N statement for this setup

Measured on **nine states across three `(S,P)`** — `(1,1)`, `(2,1)`, `(2,2)` — the norm determinant
vanishes as

```
det M^AA ~ eps^nu ,     nu = Sum_{i<j}(mu_i - mu_j) + 8 = 3 mu_1 + mu_2 - mu_3 - 3 mu_4 + 8
```

| (S,P) | mu | M | order | Sum(mu_i-mu_j) | difference |
|---|---|---|---|---|---|
| (1,1) | {2,2,0,0} | {0,0,2,2} | 16 | 8 | 8 |
| (1,1) | {2,1,1,0} | {0,1,1,2} | 14 | 6 | 8 |
| (1,1) | {1,1,1,1} | {1,1,1,1} | 8 | 0 | 8 |
| (2,1) | {3,3,0,0} | {0,0,3,3} | 20 | 12 | 8 |
| (2,1) | {3,2,1,0} | {0,1,2,3} | 18 | 10 | 8 |
| (2,1) | {2,2,1,1} | {1,1,2,2} | 12 | 4 | 8 |
| (2,2) | {4,4,0,0} | {0,0,4,4} | 24 | 16 | 8 |
| (2,2) | {4,3,1,0} | {0,1,3,4} | 22 | 14 | 8 |
| (2,2) | {2,2,2,2} | {2,2,2,2} | 8 | 0 | 8 |

`Sum_{i<j}(mu_i - mu_j)` is exactly the Appendix L/N exponent of `2603.26617`
(`nu = -Sum_{j<k}(lambda_j - lambda_k)`); the constant `+8` is the price of the monic normalisation
used here, and it disappears in the regular basis (below). The formula holds whether or not `q_2`
has poles (`M_2 > M_1` for `{2,1,1,0}`, `{3,2,1,0}`, `{4,3,1,0}`; not for the others).

### The twist prefactor factorises

Verified against **two independent unimodular rays** (`psiA = {1,1/3,1/5,-23/15}` and
`psiB = {2,-1/2,1/3,-11/6}`) on four states including one at `(2,1)`:

```
det M^AA = (phi_1-phi_2)^4 (phi_3-phi_4)^4 Prod_{i<j}(phi_i-phi_j)^(mu_i-mu_j) x N_untw + ...
```

with `phi_i - phi_j = I eps (psi_i - psi_j)`. `Prod (phi_i-phi_j)^(mu_i-mu_j)` is the paper's `(L1)`
prefactor; the extra `(phi_1-phi_2)^4 (phi_3-phi_4)^4` are the two divided differences, one per
side, each entering the four rows that carry the second branch. Ratio of leading coefficients
between the two rays matches the prediction to `1e-4` (the truncation noise level) on every state.

**`N_untw` is twist-independent** — computed from the two rays it agrees to the same noise level:

| mu | N_untw (ray A) | N_untw (ray B) |
|---|---|---|
| {2,2,0,0} | 1356.2689 - 417.18443 I | 1356.3099 - 417.05095 I |
| {2,1,1,0} | -37.935589 + 11.668879 I | -37.936737 + 11.665145 I |
| {1,1,1,1} | 1425.6926 - 438.53896 I | 1425.7357 - 438.39865 I |

All three share one complex phase (ratio `-3.251`), a universal constant coming from the measure and
the `Exp[2 Pi I theta_alpha]` of the `alpha = 2` rows, so `N_untw = K x (real, state-dependent)`.

### Regular basis: the Appendix N pattern, explicitly

Rebuilding the same 8x8 from the **divided-difference** (regular) basis instead of the monic
branches removes the `+8` exactly — the order becomes `Sum_{i<j}(mu_i-mu_j)`, i.e. `8, 6, 0` for the
three `(1,1)` states — because the basis change contributes `eps^(4(p_dir - k_dir) + 4(p_dual - k_dual)) = eps^-8`.
The regularised limits are the untwisted Q-functions: for the vacuum the second direct one tends to
`I(psi_2-psi_1) P_2` with `P_2 = u`, and for `{2,1,1,0}` two iterations give a degree-2 polynomial,
i.e. `P_2` up to normalisation and an admixture of `P_1` (the expected freedom).

Entry orders in the regular basis (vacuum), rows `(alpha,i,a)`, columns `(r,beta)`:

```
alpha = 1 :  3 2 3 2 3 2 3 2        <- suppressed, beta-dependently
             2 1 2 1 2 1 2 1
             2 1 2 1 2 1 2 1
             1 0 1 0 1 0 1 0
alpha = 2 :  0 0 0 0 0 0 0 0        <- all O(1)
             ...
```

This is Appendix N's structure directly: the `alpha = 1` functionals vanish with a `beta`-dependent
power while the `alpha = 2` rows stay finite, so the leading term is carried by a reduced minor.

**Still open:** identifying `N_untw` in closed form as an explicit determinant of the untwisted
`P_k` (the analogue of `(N4)`'s reduced determinant). The entry orders above say which minor should
carry it; the numbers in the table are the targets to reproduce.
