# `SU4_twisted_Baxter_L2.wb` — companion reference

Twisted gl(4) `L=2` Baxter (TQ) solver for site reps `[S,S,0,0] x [P,P,0,0]`, plus a
perturbation theory in the twist. **Self-contained and functional only**: no GT patterns, no
monodromy matrices, no eigensystems — everything is a polynomial identity in `u`.

Written for future agents picking this up cold. Sections 1–7 predate 2026-09-09; sections 8–10
were added on 2026-09-09.

---

## 0. Read this first

**If you only want to solve something, use `JointSolve` (Section 10.2).** It is the most
capable and best-validated entry point. Everything before it is either the exact
(non-perturbative) machinery it is built on, or the record of how it was arrived at.

**Loading into a fresh kernel.** Cells are *not* all needed and some are slow. The minimal
definition set is:

| cells | what they define |
|---|---|
| 3–8 | `z`, `χ`, `ν`, `qdet`, `Ξ`, `subl`, `τ`, `τ3quad`, `tc`, `tcvars` |
| 9 | `τvac` (also runs a structural check; ~0.5 s) |
| 11 | `Mdeg`, `qmonic`, `BaxRes` |
| 14 | `schur`, `tensorDecomp`, `dimGL4`, `irreps` — **easy to miss; skipping it makes `irreps[S,P]` stay unevaluated and downstream maps silently return `irreps[False]`** |
| 15–17 | `TwistedSolve`, `qpoly`, `verifyTwisted`, `num` |
| 39–40 | `psi`, `th`, `u0`, `zser`, `qnorm`, `PertSetup` |
| 44–47 | the ported untwisted machinery (Section 8.1) |
| 49 | `Order1` |
| 68 | `JointSolve` |

Cells 20–22, 29–34, 36–37 are validation sweeps and can be skipped.

**`Solve::svars` warnings are benign** and appear all over this notebook. They mean the system
is underdetermined, which is expected almost everywhere here (see *the lag*, below). Many cells
therefore show as `[err]` in the notebook listing while being perfectly fine. Do not "fix" them.

---

## 1. Conventions

- `z[j]`, `j=1..4` are the primary twist symbols; `φ[j]` is the exponential parametrisation
  `z[j] = Exp[h φ[j]]`, used only via the substitution rules `zφ` / `φz`.
- `Q_j(u) = z_j^(u/h) q_j(u) = Exp[φ_j u] q_j(u)`; `q_j` is a polynomial.
- `χ[a]` = elementary symmetric polynomials in the `z[j]`, so the `z[j]` are the roots of
  `z^4 - χ1 z^3 + χ2 z^2 - χ3 z + χ4`. This identity holds automatically (verified, cell 4) —
  never impose it by hand.
- Sites carry inhomogeneities `θ1, θ2`; shift `h`. `μ` labels a gl(4) **weight** (the irrep, for
  a highest-weight state).
- `ν[S,P][i][u] = (u-θ1-h λ1_i)(u-θ2-h λ2_i)`. For rectangular reps these collapse to two
  distinct functions: `ν1 = ν2` and `ν3 = ν4`. That degeneracy is what makes `τ3` factorise.
- **Degrees:** `M_k = Mdeg[S,P,μ] = (S+P) - μ_k`, and `Mhat_k = M_k + (k-1)` are the *untwisted*
  degrees.
- **Perturbation theory** uses `z[j] = Exp[I eps psi[j]]` with `psi = {1, 1/3, 1/5, 1/7}`,
  `th = {θ1 -> 1/3, θ2 -> 1/7, h -> 1}`, and normalisation point `u0 = 0`.

---

## 2. Sections 1–7: the exact (non-perturbative) solver

### Transfer matrices — only five unknowns

`τ_a(u)` has degree `2a` with leading coefficient `χ_a`. Two structural facts reduce 12 free
coefficients to 5:

- `τ3 = ν1(u-2h) ν4(u) x (quadratic with leading χ3)` — because `ν1=ν2`, `ν3=ν4`.
- **Subleading coefficient, universal:**
  `[u^(2a-1)] τ_a = -a (θ1+θ2+(a-1)h) χ_a - h Ξ_a`, with `Ξ_a = Σ_i n_i z_i ∂χ_a/∂z_i = D_n χ_a`
  (the Cartan-weight derivation acting on `χ_a`). Verified exactly for `a = 0..4`.

Leftover unknowns `tc[1..5]`: `tc[1]` = constant of `τ1`; `tc[2],tc[3],tc[4]` = the `u^2,u^1,u^0`
of `τ2`; `tc[5]` = the constant in `τ3`'s quadratic. `τ4` is the quantum determinant, fully known.

**Why exactly 5** (Paul, 2026-09-09): the system has 8 integrals of motion; 3 are Cartans (the
weight `μ`, which is input via `subl`), leaving 5. The `u^1` coefficient of `τ2` is *not*
determined in the twisted case, unlike the untwisted one — do not try to eliminate it.

### Key symbols

| symbol | meaning |
|---|---|
| `τ[a][S,P,μ][u]` | transfer matrix, `a=1..4`, in terms of `tc[1..5]` |
| `τvac[S,P][a][u]` | pseudovacuum eigenvalue, closed form, no rep data needed |
| `BaxRes[S,P,μ][k][qf][x]` | branch-`k` Baxter residual for `Q_k = z_k^(u/h) qf`. **This is the undivided operator with the exponential stripped** — see Section 9 |
| `Mdeg[S,P,μ]` | `M_k = (S+P) - μ_k` |
| `TwistedSolve[S,P,μ,sub]` | exact joint bilinear solve for `tc` + monic `q_k` |
| `verifyTwisted[...]` | counts nonzero residual coefficients; `0` = pass |
| `qpoly[clist, u]` | polynomial from an ascending coefficient list |
| `Wronsk`, `Vand`, `qWronsk`, `AA` | Casoratian machinery, Vandermonde, vacuum factor `A(u)` |

### Established (verified, not conjecture)

- **Degree formula derived, not assumed.** For arbitrary `q` of degree `M`, the top residual
  coefficient is `charPoly[z_k] = 0` automatically; the next is
  `h z_k Π_{j≠k}(z_k - z_j)(μ_k - S - P + M)`, zero iff `M_k = S+P-μ_k`.
- **Number of solutions = weight multiplicity** (the twist breaks su(4)).
- **Validation on the six permuted vacua** at `(1,1),(2,1),(2,2),(3,3),(5,4)`: 30 states, all
  unique, Baxter residual **exactly 0** on all four branches. Repeated at `(1,1)` with *no*
  substitution at all (everything symbolic).
  - Shift ordering for the permuted vacua is **not** index order: `ord = Join[A, complement A]`
    with `A = {i : n_i = S+P}`. Found empirically.
- **Wronskians:** `N_I = V_I / (Π_{i∈I} z_i)^(|I|-1)`, `q_I = W_I/V_I` monic; `q_I` regular as
  `φ → 0` exactly when `I` is a lower set for the degrees `M`.
- **`q_1234` closed form**, state- and twist-independent:
  `q_1234(u) = A(u) A(u-h)`, `A(u) = Π_{j=1..S}(u-θ1-jh) Π_{j=1..P}(u-θ2-jh)`.
- **Divergence count** under *monic* normalisation: `q_k = O(eps^(-(M_k - M_1)))`, finer pattern
  `[u^j] q_k = O(eps^(-(M_k - max(j,M_1))))`. (Input verified explicitly for `k=3`; asserted for
  `k=4`.)

### Test states (reused throughout)

With `z_j = 1 + eps psi_j`, `psi = {1,1/3,1/5,1/7}`, `θ1=1/3, θ2=1/7, h=1`:

| | `(S,P)` | `μ` | `M` | why |
|---|---|---|---|---|
| A | (1,1) | `{2,1,1,0}` | `{0,1,1,2}` | generic |
| B | (1,1) | `{0,0,2,2}` | `{2,2,0,0}` | lowest weight; flag reversed |
| C | (1,1) | `{1,1,1,1}` | `{1,1,1,1}` | control — all `M_k` equal, nothing diverges |
| D | (3,1) | `{3,3,1,1}` | `{1,1,3,3}` | `M_1>0` **and** real divergence |

---

## 3. Section 8: perturbation theory in the twist

### The normalisation change

`q_k(u0) = 1` instead of monic. Under monic the coefficients diverge like `eps^(-(M_k-M_1))`;
normalising at a point makes every coefficient **regular** (plain Taylor, no Laurent), at the
price of letting the *degree drop* in the limit — the leading coefficients switch off smoothly.

**Measured pattern** (`μ={4,3,1,0}` at `(2,2)`, `M={0,1,3,4}`, `M_1=0`): `[u^j] p_k = O(eps^j)`
exactly, for every branch. In general `[u^j] p_k = O(eps^(max(j,M_1) - M_1))`.

### 8.1 The leading order is already known exactly — do not re-solve it

`Experiments/baxter_SU4_L2_XXX.wb` solves the **untwisted** problem in closed form for arbitrary
`(S,P)` and `μ`, verified there up to `(S,P)=(8,8)`. Its `degQ[mu,S,P]` **is** this notebook's
`Mhat`. Ported here with a `U` suffix (that notebook keeps `t`, `S`, `P`, `ν_i`, `irreps`,
`dimGL4` as globals and shares kernel K1 with this one, so it cannot simply be run alongside).

| symbol | meaning |
|---|---|
| `Cas2U[μ]`, `Cas4U[μ]` | Casimirs from the shifted weights |
| `degU[μ,S,P]` | `= Mhat`; asserted equal to `Mdeg + Range[0,3]`, checked over all irreps to `(3,2)` |
| `tau1U`, `tau2cU`, `tau3U`, `tau4U` | closed-form untwisted transfer matrices |
| `QsolveU[μ,S,P]` | four gauge-fixed monic untwisted `P_1..P_4`, degrees `degU` |
| `verifyU[μ,S,P]` | checks `P_k` against this notebook's own `BaxRes` at `z→1` |
| `tc0U[μ,S,P]` | **the bridge**: this notebook's five `tc` at order `eps^0`, no solving |

Verified: `verifyU` True on **15/15** irreps across `(S,P) = (1..2, 1..2)`; `tc0U` reproduces
`τvac` at the pseudovacuum and `τ4` agrees identically.

> **Why you must use `tc0U` rather than solving order 0.** The `Mhat` degree bound also admits
> the lower-degree solutions of *other* irreps in the same tensor product, so an order-0 solve
> returns spurious extra branches. At `(1,1)` with `μ={2,2,0,0}` it returns four, of which one is
> the pseudovacuum and one is actually `μ={2,1,1,0}`.

### 8.2–8.3 First and second order (superseded by Section 10.2)

`Order1[S,P,μ]` — order-0 data fixed from `tc0U` and `P_1`; order 1 then linear. Order 2 is
genuinely **quadratic** (degree 2 in the order-1 free parameters, 1 in the order-2 unknowns):
that is the secular equation of degenerate perturbation theory, not a setup artifact.

Kept for the record and as a cross-check on `JointSolve`. Prefer `JointSolve`.

### 8.5 The exact twisted `q_k` is unique

`ExactQ[S,P,μ,tcv,dbound]`: given `τ`, each branch's Baxter equation is **linear** in `q_k`, so
the exact twisted solution is cheap. Result, symbolic in all four `z[j]`:

- unique, **zero free parameters**, at both degree bounds `M` and `Mhat`;
- coefficients beyond `M_k` vanish **identically in `z`**.

So `M_k` is the right degree bound; `Mhat` is harmless over-parametrisation.

---

## 4. Section 9: the multi-exponential ansatz — a dead end, and why

`BaxOp` is the **undivided** Baxter operator, with `Z[j]` standing for `z[j]^(u/h)`.

**Result: the ansatz `Q = Σ_j Exp[φ_j u] p_j(u)` decouples exactly.** The coefficient of `Z[j]`
involves only `p_j`, and equals `z_j^(-3) BaxRes_j[p_j]` identically. Verified, including that no
cross terms appear.

Consequence: at generic twist each `p_j` is forced to be a multiple of `q_j`, so the general
solution is `Q = Σ_j c_j(eps) Exp[φ_j u] q_j(u)` with **scalar** `c_j` — no new polynomial
freedom, and the normalisation of `q_j` is irrelevant (absorbed into `c_j`).

**What was tried and did not work.** With eps-independent `c_j`, imposing vanishing to order
`eps^(k-1)` reaches only `span{P_1, P_2}`: at the pseudovacuum the limits have degrees
`{0,1,1,1}` against `Mhat = {0,1,4,5}`. Allowing Laurent `c_j` with poles up to `eps^-3` did not
reach `P_3`, `P_4` either (cell 60). **This is unresolved** — the arithmetic says degree
`M_k + (k-1)` should be reachable, so either more pole orders are needed or one of those two
cells has a bug. Treat cell 60's negative result as *not established*.

---

## 5. Section 10: the working scheme

### 10.1 Deforming the untwisted solution

`BaxU2` = undivided operator on a plain polynomial (nothing stripped). Push `P_k` out in `eps`
with `deg Q^(n) = deg P_k + n` — the degree **grows one per order**, because expanding
`Exp[φ u]` manufactures powers of `u`. Each order is linear.

Validated to order 3 against the exact `Q = Exp[I eps psi_1 u]` (the `k=1` pseudovacuum branch,
where `q_1 = 1`): `Q^(n)` and all 15 `tc` corrections match exactly.

**But a single branch can never fix `τ`.** Feeding in the exact `Q = Exp[φ_1 u]·1` and imposing
the branch-1 equation at orders 1, 2 *and* 3 gives 21 conditions on 15 unknowns and still leaves
`tc[1]`, `tc[2]` free. This is the fact already recorded in `HANDOFF.md` — "branch 1 alone leaves
a 2-parameter family; all four together are unique" — and `tc[1]`, `tc[2]` are exactly those two
parameters.

### 10.2 `JointSolve` — use this

```
JointSolve[S, P, μ, N, u]
```

Four branches, **one shared τ**, exponential ansatz

```
Q_k(u) = Exp[φ_k u] p_k(u),   deg p_k = M_k,   p_k(u0) = 1,   p_k^(0) = P_1 normalised,
```

all orders `1..N` solved together. Returns an association with keys
`"M"`, `"q0"`, `"#eqs"`, `"#unk"`, `"#sol"`, `"sol"`, `"unk"`, `"free"`.

**The form of the answer.** Each solution is exactly `Exp[φ_k u] x` a polynomial of degree `M_k`,
**fixed non-perturbatively** — the degree never grows. Only the `M_k` coefficients are expanded,
each as a plain Taylor series. In the limit the *leading coefficient vanishes* at rate
`eps^(M_k - M_1)`, so the polynomial degenerates onto `P_1` without the degree ever having grown.

**The lag.** Solving to order `N` determines everything through roughly order `N-1`; the top
orders keep a few free directions, and the exact answer always lies inside that family. The lag
**scales with the degree spread**: the pseudovacuum (`M={0,0,2,2}`) lags 1 order, while
`μ={4,3,1,0}` (`M={0,1,3,4}`) lags 2, and the `p_k` free coefficients reach back further for the
branches with the larger `M_k - M_1`. Budget `N ≈ target + (M_max - M_1)/2`, and always solve at
least one order beyond what you need.

`residualOrders[st,S,P,μ,N,u]` validates a solution by substituting back into all four branch
residuals; use it when no closed form is available.

### Validation status

| state | `(S,P)` | `M` | result |
|---|---|---|---|
| pseudovacuum `{2,2,0,0}` | (1,1) | `{0,0,2,2}` | `τ` **exact through order 5** vs closed form (`N=6`); `p_k^(1)` exact |
| state D `{3,3,1,1}` | (3,1) | `{1,1,3,3}` | `tc^(1)` definite, satisfies Section 8's order-1 relations exactly |
| `{4,3,1,0}` | (2,2) | `{0,1,3,4}` | residual vanishes identically at **every order 0–6 on all four branches**; `τ` determined through order 4 |

`{4,3,1,0}` at `(2,2)` is the most demanding state available at that size: the only irrep in the
whole decomposition with all four Q-degrees distinct, dim 175, and genuinely non-rectangular.

**Timing** (warm cache): `N=6` takes ~3.3 s for the pseudovacuum and ~9.6 s for `{4,3,1,0}`.
Runtime is not the constraint. Cost grows faster than linearly in `N` (`Solve` on a polynomial
system).

---

## 6. Traps — all of these cost real time

**Wolfram / notebook**

- **Line continuation.** Breaking a definition where the expression is *already complete* makes
  the front end start a **new statement**. This silently dropped `tau3U`'s quadratic factor and
  `tau2U`'s `Cas4` term. Keep continued lines ending with their operator. Symptom: a `τ` of the
  wrong degree, plus stray `Part::partd` messages from the orphaned statement evaluating with an
  unbound pattern variable.
- **`CoefficientList[0, x]` returns `{}`, not `{0}`.** A test like `Union[...] === {0}` therefore
  rejects exactly the orders that *did* vanish. Test the coefficient against `0` directly.
- **Module-localised spectral variable** leaks a scoped `x$nnn` into returned expressions; the
  caller's `CoefficientList` then reads every coefficient as degree 0. Pass the variable in.
- **`Exponent[list, list]` threads componentwise.** For total degree in a *set* of variables,
  scale them all by one marker: `Exponent[expr /. Thread[vars -> s vars], s]`.
- **`&` binds looser than `/@`** — but before blaming precedence, check the function is actually
  defined. `irreps[1,1]` staying unevaluated (cell 14 not run) produces the identical symptom.
- **`Limit[]` on a whole rational expression is unreliable** on multiplicity>1 states — it
  returned different branches on repeated runs of identical code. Take limits coefficient by
  coefficient.
- **Detecting admissible untwisted degrees by nullspace dimension** silently misses a jump when
  two admissible degrees are adjacent. Test **monic solvability** instead.
- **`RootReduce`** is exact but far too slow to sweep with; bare `PossibleZeroQ` falls back to an
  undecided numeric test. Do rigorous checks on the **rational** (vacuum) states, where
  `Expand[Together[...]]` collapses to a literal 0.
- **Do not confirm anything by `N[]`** where a Vandermonde has been divided out — dividing by
  `V ~ eps^6` is catastrophic numerically.
- **Kernel `abort` has repeatedly failed** on long `Solve`/`NSolve` (2026-09-02 and 2026-09-09);
  a `restart` was needed, losing state.

**Notebook-specific**

- `baxter_SU4_L2_XXX.wb` shares kernel K1 with this notebook and collides on `t`, `S`, `P`,
  `ν1..ν4`, `irreps`, `dimGL4`. **Do not run it in the same kernel** — the port in Section 8.1
  exists precisely to avoid this.
- Substitute `th` into *everything* derived from `QsolveU`, which is symbolic in `θ1,θ2,h`, not
  just into `tc0U`. Missing this leaves the order-0 residual silently nonzero — and it is
  **invisible at the pseudovacuum**, where `P_1 = 1` carries no `θ` dependence. Test new code on
  state D, not the pseudovacuum.
- When perturbing, never solve for order-`n` unknowns *together with* leftover order-`(n-1)`
  free parameters: `BaxRes` is bilinear in `tc x qc`, so that makes the system quadratic and
  `Solve` returns spurious branches.

---

## 7. Open threads

- **Section 9 is unresolved** (see above): whether `P_3`, `P_4` are reachable as limits of
  Laurent combinations of the twisted branches. Cell 60 says no up to `eps^-3`; the degree
  arithmetic says it should be yes. Recheck before building on either answer.
- Nothing here is wired to `SU4_norm_untwisting_L2_general.wb` or the FSoV norm.
- No scan over reps.
- Nothing promoted to `Clean/`; `paper/xxxCG.tex` is still gl(2)/single-site, so writing any of
  this up is a scope decision for Paul.
