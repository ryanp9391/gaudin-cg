# Small-twist expansion of the FSoV scalar product — the two possible routes

Companion note to `SU4_SoV_SmallTwist_L2.wb` (gl(4), L=2, sites `[S,S,0,0] x [P,P,0,0]`).

**Goal.** Feed the twist expansion, on a unimodular ray `phi_1+phi_2+phi_3+phi_4 = 0`, into the SoV
scalar product of `Clean/SU4_SoV_ScalarProduct_L2.wb`, and find the leading behaviour as the twist
is switched off: untwisted Q-functions times overall twist factors whose powers are fixed by the
quantum numbers of the state. First checkpoint: **orthogonality of distinct states must survive the
expansion**, `det M^{AB} = 0` for `A != B` and `det M^{AA} != 0`, order by order in the twist.

**The ray.** `z_j = Exp[I eps psi_j]`, `psi = {1, 1/3, 1/5, -23/15}`, `Sum_j psi_j = 0`, so
`z_1 z_2 z_3 z_4 = 1` exactly (the scalar product requires a unimodular twist; the Experiments
notebook's `psi = {1,1/3,1/5,1/7}` is **not** unimodular and must not be used here). Everything is a
rational power series in `eps`: the leading, untwisted problem is unique and rational, so no
algebraic branches ever appear — unlike the exact finite-twist solver `TwistedSolve`, whose
multiplicity-`>1` weight sectors produce nested radicals.

---

## The fork

The four twisted Q-functions can be organised in two inequivalent ways, and the small-twist limit of
`M^{AB}` means something different in each. Both solve the same twisted Baxter equation with one
shared `tau`; they differ in which basis of the four-dimensional solution space is followed as
`phi -> 0`.

### Route 1 — fixed-degree exponential basis (`JointSolve`)

`Q_k(u) = Exp[phi_k u] q_k(u)` with `deg q_k = M_k = (S+P) - mu_k` **fixed at every order in eps**;
normalisation `q_k(0) = 1`; all four branches expanded around the *same* untwisted solution `P_1`.

- This is the basis the Clean scalar product currently consumes, and it is a genuine solution:
  validated against the exact twisted solution in cell 70 of
  `Experiments/SU4_twisted_Baxter_L2.wb`.
- As `phi -> 0` the four branches **collide**: every `q_k -> P_1`, and likewise the two duals. There
  is no contradiction with `P_1..P_4` being distinct — a monic `q_k` of degree `M_k` cannot tend to
  `P_k` of degree `Mhat_k = M_k + (k-1) > M_k`; the only untwisted solutions of degree `<= M_k`
  are the low ones.
- Measured consequence, `(S,P) = (1,1)`, `V = {2,2,0,0}` vs `A = {2,1,1,0}`, solved to `eps^7` and
  truncated at `eps^5`:

  | rows | `beta = 0` cols | `beta = 1` cols |
  |---|---|---|
  | `alpha = 1` | `O(eps^5)` (at truncation, so `>= 5`) | `O(eps^4)` |
  | `alpha = 2` | `O(eps^2)` | `O(eps^2)` |

  i.e. `M^{AB} = [[Z, X],[Y, ...]]` with the `alpha=1` rows suppressed and the suppression
  `beta`-dependent — the gl(4) L=2 analogue of Appendix N of `2603.26617v1` (`References/`), where
  the `alpha=1` row vanishes as `phi^{2L-M-1} = phi^{2(L-lambda_1-lambda_2)-beta}` and the leading
  determinant is carried by the reduced minor, eq. (N4):

      B -> i (-1)^L (2 phi)^{L-M_A-M_B} / (L-M_A-M_B)! x det_{1<=alpha,beta<=L-1} <Q_A u^{beta-1} Q_B>_{L,alpha+1}

  — one row (`alpha=1`) and one column (`beta=L`) dropped, the twist prefactor restoring finiteness.
- **Obstruction as it stands:** both candidate leading blocks (the `4x4` `(alpha=1, beta=1)` block
  and the `4x4` `(alpha=2, beta=0)` block) have **rank 1**, for the norm exactly as for the cross.
  All four rows `(i,j)` degenerate because `q_1, q_2` collide and `q^3, q^4` collide. So the leading
  determinant vanishes for every pair and orthogonality is not yet decidable at leading order; the
  information sits in the orders that lift the degeneracy.

### Route 2 — rotated / regular basis, deforming `P_1..P_4` (`BranchSolve3`) — **being tried first**

Follow the four **distinct** untwisted solutions `P_k` (staggered degrees
`Mhat_k = degU = M_k + (k-1)`, e.g. `{0,1,4,5}` for `V`) and switch the twist on.

- With `deg p_k = Mhat_k` held fixed the system is **inconsistent**: branch 4 already fails at order
  `eps^1` (12 equations, 11 unknowns); branches 1–3 are fine.
- It becomes consistent as soon as the degree is allowed to grow by one per order,
  `deg p_k^(n) = Mhat_k + n`. All four branches then solve simultaneously with one shared `tau`,
  starting from the four distinct `P_k`.
- That degree growth **is** the rotation: the solution tending to `P_k` is not
  `Exp[phi_k u] x (fixed-degree polynomial)` but a combination of branches with `1/eps`-type
  coefficients, and expanding `Exp[(phi_j - phi_k) u]` raises the degree by one per order. Same
  mechanism as the `lambdahat_j = lambda_j + (j-1)` degree jump and the `Q_2 -> Q_2 + r Q_1`,
  `r = 1/phi_2` rotation of the 2026-09-02 thread (`SU4_norm_untwisting_L2*.wb`).
- **Where it stops (2026-09-10).** Solving order by order (linear at each order, seconds — a single
  `Solve` over all orders is quadratic across orders, slow, and returns spurious branches), the
  system is consistent at orders `eps^1` and `eps^2` and **inconsistent at `eps^3`**
  (42 equations, 31 unknowns, no solution). Neither dropping the gauge conditions nor faster degree
  growth helps: `grow = 2` already fails at order `eps^2`.

  Diagnosis: the mixing coefficients are **Laurent**, not Taylor. A rotated branch is
  `Q_2 + (a/eps) Q_1 = Exp[phi_2 u](q_2 + (a/eps) Exp[(phi_1-phi_2)u] q_1)`, so after normalising
  the leading order to `eps^0` the degree grows as `deg q_1 + n + 1`, with an offset — not the
  `Mhat_k + n` assumed. Solving the Baxter equation directly in this ansatz is therefore the wrong
  way round.

  **Proposed fix (not yet done):** do not re-solve. Build the rotated basis explicitly from the
  Route-1 solutions, which are already known to `eps^7`: the Baxter equation is linear, so
  `Qtilde_2 = (Q_2 - Q_1)/eps` is a solution too, and it is finite with a leading term of degree 1,
  i.e. `P_2` — and similarly up the tower. Rotating the `q_i` by a triangular, `eps`-dependent
  matrix multiplies `det M^{AB}` by a known factor, so it cannot manufacture orthogonality, but it
  does redistribute the `eps` orders and should lift the rank-1 degeneracy of the leading blocks.
  Note the **dual** side already uses wedge objects (`q^j = W_{jbar}/(V_{jbar} A)`); it is only the
  `i in {1,2}` rows that sit in the colliding basis.

- **Open point:** leftover mixing freedom (`pcn[1,1,n]`, `pcn[3,5,n]`, ...) — how much of a lower
  `P_j` sits inside branch `k`. It is not gauge: it selects which twisted solution is being
  followed. A convention must be fixed before this feeds the scalar product; the natural one is the
  Wronskian/wedge combinations `q_1`, `W_12/V_12`, `W_123/V_123` of Sections 6–7 of
  `SU4_twisted_Baxter_L2.wb`, which are exactly the objects that stay regular as `phi -> 0`
  (`q_I` regular iff `I` is a lower set for the degrees `M`).

---

## Plan

1. **Route 2 first.** Fix the mixing convention (wedge basis), rebuild `q_i` and the duals `q^j`
   from it, feed into `spMatFSoV`, expand in `eps`, and check `det M^{AB} = 0` for `A != B` with
   `det M^{AA} != 0` at leading order — with the leading blocks now full rank, not rank 1.
2. **If Route 2 does not work, fall back to Route 1**: keep the colliding basis and extract the
   answer the Appendix-N way — drop the suppressed row and column, keep the reduced minor, and read
   off the overall twist power from the quantum numbers (`M_k`).

## Traps already paid for

- The twist must be unimodular; `TwistedSolve` and the Clean scalar product both assume it.
- The exact finite-twist solver produces radicals for multiplicity `> 1` weights. The perturbative
  route is rational at every order — use it.
- `PossibleZeroQ` silently misses the exact zeros of these series coefficients: `Simplify` every
  coefficient before testing.
- `FirstPosition` searches heads and returns `{0}` for a list that itself matches the pattern; use
  `LengthWhile`.
- An entry order equal to the truncation order is truncation noise. Solve two orders deeper than
  the order you intend to trust (solve to `eps^7`, trust `eps^5`).
- `ToExpression` on a raw `.wb` cell string splits at every newline that ends a complete
  expression, silently truncating a definition whose RHS continues on the next line (`BaxRes`).
  `joinCont` in cell 0 glues continuation lines back on; verify afterwards that `BaxRes` depends on
  the `tc`.
- Kernel `abort` does not always take on a long `Solve`; budget for a restart.
- `QsolveU` gives only the untwisted **highest-weight** solution per weight, so the second state of
  a multiplicity-2 weight (the descendant) is not reachable this way — the same-weight sharp test
  needs separate work.
