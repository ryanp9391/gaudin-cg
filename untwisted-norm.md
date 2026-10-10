# Untwisted norm & orthogonality — gl(4), L=2, rectangular sites (session 2026-09-11)

Handoff for the untwisted (`z_j = 1`, `phi = 0`) SoV scalar product and the attempt to untwist the
twisted 8x8 norm. Written for a future session picking this up cold. Everything here is at
`(S,P) = (1,1), (2,1), (2,2)` with `theta1 = 1/3, theta2 = 1/7, h = 1` unless stated otherwise.
Notation: `f^{[±2n]}(u) = f(u ± n h)`; sites `[S,S,0,0] x [P,P,0,0]`;
`nu1 = nu2 = (u-theta1-Sh)(u-theta2-Ph)`, `nu3 = nu4 = (u-theta1)(u-theta2)`.

**Files**

| file | what |
|---|---|
| `Paul/Mathematica/XXX/Experiments/SU4_SoV_Orthogonality_Untwisted_L2.wb` | **the clean, final zero-twist construction** (runs cold, top to bottom) |
| `Paul/Mathematica/XXX/Experiments/SU4_SoV_Orthogonality_Untwisted_L2.tex` / `.pdf` | write-up of the above (5 pp) |
| `Paul/Mathematica/XXX/Experiments/SU4_IoM_Count_L3.wb` | counting independent integrals of motion at L = 2, 3, 4 (explicit transfer matrices) |
| `Paul/Mathematica/XXX/Experiments/SU4_Untwisting_Norm_L2.wb` | exploratory: untwisting the twisted norm (Appendix-N analysis, exact route). Messy, but every result and dead end is recorded in markdown cells |
| `References/2603.26617v1.pdf`, Appendix N | the su(2) untwisting mechanism we tried to repeat |

---

## Part 1 — THE RESULT: the 2x2 scalar product at zero twist

### 1.1 The idea

The twisted orthogonality identity (Clean `SU4_SoV_ScalarProduct_L2`) is

    0 = ∮ mu q^j_A (O^A - O^B) q_i^{B[2]}

and is linear in **eight** differences of integrals of motion: the five `Delta tc` and three
weight combinations sitting in the subleading coefficients of `tau_1, tau_2, tau_3`. At zero twist
six of these become dependent, leaving **two**. So the 8x8 system collapses to a **2x2** system.

### 1.2 The dependent integrals of motion (proved)

At `z = 1` the subleading coefficient
`[u^{2a-1}] tau_a = -a(theta1+theta2+(a-1)h) C(4,a) - h C(3,a-1) Sum(mu)` no longer depends on the
individual weights. The **three weight directions drop out**. The transfer-matrix eigenvalues depend
on the irrep only through the Casimirs `C2, C4`. The closed forms come from `baxter_SU4_L2_XXX.wb`,
cells 5, 7, 9, 10, 17. With `tau1vac = sum_i nu_i = 2 nu1 + 2 nu4`, for any two states:

    Delta T1 = Delta tc1
    Delta T2 = (1/4) Delta tc1 (tau1vac + tau1vac^{[-2]}) + Delta c
    Delta T3 = nu4 nu1^{[-4]} Delta tc1

with the two independent invariants

    Delta tc1 = (h^2/2) Delta C2
    Delta c   = h^4 ( (1/8) Delta(C2^2) - (1/4) Delta C4 + (1+2S+2P)/4 Delta C2 )

Component-wise, `Delta tc5 = Delta tc1`, `Delta tc2 = 2 Delta tc1`,
`Delta tc3 = -(2(theta1+theta2) + (S+P+2)h) Delta tc1`, and `Delta tc4 = Delta c + ...`.

These are proved **symbolically** in `theta1, theta2, h, S, P` and two generic weights with
`Sum mu = 2(S+P)`. They are also checked against the Baxter-solved `tc` of all 12 irreps, and
directly as `Delta tc = Delta tc1 (1,2,kappa3,0,1) + Delta tc4 (0,0,0,1,0)` for 54/54 pairs at two
unrelated points (`theta = (1/3,1/7), h=1` and `(2/5,1/11), h=2/3`).

Casimirs, with shifted weights `x_i = mu_i + 5/2 - i` and `p_k = sum x_i^k`:
`C2 = p2 - 5` and `C4 = p4 + 4 p3 - p1 p2 + p2/2 - 2 p1^2 - 10 p1 - 51/4`.

### 1.3 The 2x2 matrix — how to build it

**Measure and functionals.** `E2 = mu/(nu1 nu1^{[-2]})`, with
`mu = prod_alpha prod_{n=0}^{S_alpha - 1} 1/(u - theta_alpha - n h)`. `L_s[.]` is the **sum of residues
over the pole ladder of site s**: `theta_s + n h`, `n = 0..S_s+1`. The two ladders are independent
functionals, because an h-periodic insertion `e^{2 pi i u/h}` is constant on each ladder and may be
added freely. This needs `theta1 - theta2 not in hZ`.

**Q-functions at zero twist** (no Hodge duals: the Vandermonde vanishes at `z = 1`):

- Direct: four monic polynomials `q_k` solving `O q_k^{[2]} = 0`, of staggered degrees
  `(S+P) - mu_k + (k-1)`.
- Dual: four monic `q^j` solving `O^dagger q^j = 0`, of degrees `mu_4, mu_3+1, mu_2+2, mu_1+3` for
  `j = 4,3,2,1`.
- Gauge: in each family, set to zero the coefficient of every lower family degree. Then one joint
  exact solve fixes `tc` and the `q_k`; the duals then follow linearly and uniquely.
- Irreps: `{S+P-j, S+P-k, k, j}`, `0 <= j <= k <= min(S,P)`, each with multiplicity one.

**The matrix, nu-form** (no T2 split needed). Fix a direct `q_i` of B and a dual `q^j` of A:

    N^{AB}_{s,1} = L_s[ E2 q^j_A ( -nu1^{[-2]} q^B_i + (1/2)(nu1 + nu1^{[-2]} + nu4 + nu4^{[-2]}) q^{B[-2]}_i - nu4 q^{B[-4]}_i ) ]
    N^{AB}_{s,2} = L_s[ E2 q^j_A q^{B[-2]}_i ]

    N^{AB} . (Delta tc1, Delta c)^T = 0

So `det N^{AB} = 0` for distinct irreps, and `det N^{AA}` is the norm.

**Difference form** (a constant column operation, determinant unchanged). Because
`(1/2)(nu1 + nu1^{[-2]} + nu4 + nu4^{[-2]}) = nu1^{[-2]} + nu4 - (h^2/2)(S+P)`, replace
column 1 by column 1 + `(h^2/2)(S+P)` times column 2. With `nabla f = f - f^{[-2]}`:

    N'^{AB}_{s,1} = L_s[ E2 q^j_A ( -nu1^{[-2]} nabla q^B_i + nu4 nabla q^{B[-2]}_i ) ]
    N'^{AB}_{s,2} = L_s[ E2 q^j_A q^{B[-2]}_i ]
    N' . (Delta tc1, Delta c')^T = 0,   Delta c' = Delta c - (h^2/2)(S+P) Delta tc1
                                         = h^4((1/8)Delta(C2^2) - (1/4)Delta C4 + (1+S+P)/4 Delta C2)

Equivalently, `Delta T2 = Delta tc1 (nu1^{[-2]} + nu4) + Delta c'`. Both the identity and the
`Delta c'` formula are proved symbolically.

### 1.4 Which Q-pairings work

A **per-ladder** row pair (same `(i,j)`, `s = 1,2`) gives `det N^{AA} != 0` for every irrep exactly
for

    (i,j) in {(2,1), (2,2), (3,1), (3,2), (3,3), (4,1), (4,2), (4,3)}

and `det = 0` for every irrep otherwise. **The four twisted-Clean pairings `i in {1,2}, j in {3,4}`
all give vanishing norms at zero twist.** Default choice: `(i,j) = (2,1)`.

**Full-contour variant** (each row over both ladders, no periodic insertion; a different `(i,j)` per
row):
- a row vanishes identically for every irrep exactly when `i < j`;
- `(1,1) ∝ (4,4)` and `(2,2) ∝ (3,2) ∝ (3,3)` are dependent;
- 41 of the 120 pairs of rows work for all 12 irreps; the other 79 fail for all.

### 1.5 Verification

At `(i,j) = (2,1)`, all irreps of `(1,1), (2,1), (2,2)`:

| check | result |
|---|---|
| norms `det N^{AA} != 0` | 12 / 12 |
| cross `det N^{AB} = 0` | 42 / 42 |
| `N^{AB} . (Delta tc1, Delta c) = 0` | 42 / 42 |
| `N^{AB} ≡ 0` (orthogonality holds but untested) | 6 |

Example at `(1,1)`:
- `q_2^{2110} = u^2 - (10/21)u` and `q^1_{2110} = u^5 - (130/21)u^4 + (3097/63)u`.
- `N^{AA}` for `A = {2,1,1,0}`: `{{-150144211/334305, 3235736/334305}, {692896376/1820105, -41930251/1820105}}`,
  with `det = 1039313041/156009`.
- Cross `{2,2,0,0}` vs `{2,1,1,0}`: `det = 0`, and `(Delta tc1, Delta c) = (2, 2)` spans the kernel.
- Norms (the `(2,1)` 2x2, difference form) used below: `N'({2,2,0,0}) = -145570417/111435 ≈ -1306.33`,
  `N'({2,1,1,0}) = 1039313041/156009 ≈ 6661.97`, `N'({1,1,1,1}) = -427/3`.

**Open caveat — vacuous pairs.** At `(2,2)` the pairs `{4,4,0,0}↔{2,2,2,2}`,
`{4,3,1,0}↔{3,2,2,1}` and `{4,2,2,0}↔{3,3,1,1}` have **every** row identically zero: all 16
full-contour and all 32 per-ladder rows, for every `(i,j)`. They are not tested by any choice. This
looks structural, like a selection rule. In the labelling `{4-j,4-k,k,j}` each pair satisfies
`(j+k)_A + (j+k)_B = 4`. It is not a Casimir degeneracy. Unexplained.

**Degenerate states.** Descendants share `tc` with their highest weight, so only one state per irrep
is separated at zero twist.

---

## Part 2 — Counting integrals of motion (`SU4_IoM_Count_L3.wb`)

Method: build the transfer matrices explicitly (GT Lax operators, monodromy, quantum minors, as in
`SU4_untwisted.wb`), sample `tau_a(u_k)` at `aL+1` points, and take the rank of `{1, tau_a(u_k) v}` for
one generic vector `v`. This is valid because the `tau_a` commute and are diagonalisable. Checks:
the `tau_a` commute, `tau4 ∝ 1`, and two random vectors agree.

| L | twisted | zero twist | lost |
|---|---|---|---|
| 2 | 8 | 2 | 6 |
| 3 | 12 | 6 | 6 |
| 4 | 16 | 10 | 6 |

Twisted `= 4L` (3 weights + `4L-3` free coefficients). Zero twist `= 4L - 6`. The number lost is 6 at
every L tested, equal to the number of positive roots of gl(4). This is observed, not derived. The
L=4 zero-twist count is unchanged between dimensions 1296 and 4320, so it is not capped. The explicit
L=3 relations (analogue of 1.2) have **not** been written down.

---

## Part 3 — Untwisting the twisted 8x8 norm (`SU4_Untwisting_Norm_L2.wb`, exploratory)

Goal (Paul): recombine the twisted 8x8 system so that exactly the two surviving IoM directions
remain, and show the untwisted 2x2 is the leading contribution as `phi -> 0`.

### 3.1 Transformed twisted system (validated)

- **Columns:** the unaligned twisted identity in coordinates
  `w = (Delta tc1, Delta c', Delta tc2 - 2 Delta tc1, Delta tc3 - kappa3 Delta tc1, Delta tc5 - Delta tc1, Delta s1, Delta s2, Delta tau3quad)`.
  Here `s_a` are the weight-dependent subleading coefficients, `kappa3 = -(2(theta1+theta2)+(S+P+2)h)`,
  and `k0 = nu1(-h) + nu4(0)`. The change is triangular with unit diagonal. Columns 1-2 are exactly the
  untwisted difference form. The other six columns are the insertions
  `u^2 q^j q_i^{[-2]}`, `u q^j q_i^{[-2]}`, `-nu4 q^j q_i^{[-4]}`, `-u nu1^{[-2]} q^j q_i`,
  `u^3 q^j q_i^{[-2]}`, `-u nu4 q^j q_i^{[-4]}`.
- **Rows:** any pairing of twisted solutions. "Regularised" `R_k = sum_l c_l Q_l`, `R^k = sum_l d_l Q^l`
  are fixed by matching the untwisted `q_k`, `q^k` at `u = 0..3`. These are exact twisted solutions
  tending to the untwisted ones.
- **Twist:** unimodular rays, e.g. `z = (1+eps, 1+eps/3, 1-eps/2, 1/prod)`.
- **Validated at eps = 1/5:** cross `|M'.w| = 0` to 69 digits, and the norm is non-zero.

### 3.2 Zero-twist 8x8 is singular; forced kernel = degenerate multiplet

- The per-ladder 8-row `M'(0)`, built with untwisted Q's, has `det M'(0) = 0` for all three `(1,1)`
  states. The kernel dimensions are 3, 2, 2 for `{2,2,0,0}`, `{2,1,1,0}`, `{1,1,1,1}`, and the kernel
  mixes `(w1,w2)` with the vanishing coordinates.
- **Forced kernel.** With all 32 rows, `dim ker` is 3, 2, **0**. So the singlet's kernel was an
  artefact of the row choice; the others are forced.
- **Mechanism, verified for `{2,1,1,0}`.** Its 14 twisted partners are the rest of the 15-dim
  multiplet (11 permutations of `{2,1,1,0}` plus 3 states at `{1,1,1,1}`). Their `w(A,A')/eps` lie in
  `ker M'(0)`: the residual falls like `eps` (0.061 → 0.0061), and they span exactly the 2-dim kernel.

### 3.3 Appendix-N analysis: block dominance (the positive result)

Rows: two **good** rows (regularised `(2,1)` pairing, per ladder, which tend to `N'`) plus six
**full-contour** rows built from **raw** twisted `Q_i, Q^j` (the Appendix-N `alpha = 1` analogue;
each product carries `e^{(phi_i-phi_j)u}`).

- **Power matrix.** In every full row the eps-powers across the eight `w` columns follow the fixed
  offset profile `(0, 0, -2, -1, -2, -3, -3, -3)`. The base power depends on the pairing:
  `(1,2)`:4, `(1,3)`:4, `(1,4)`:5, `(2,3)`:3, `(2,4)`:4, `(3,4)`:4, `(2,1)`:2, `(3,1)`:2, `(4,1)`:1,
  `(3,2)`:3, `(4,2)`:2, `(4,3)`:2. Diagonal pairings are zero in columns 1-5 and O(1) in columns 6-8.
  So columns 1-2 are the most suppressed, and the minimal Laplace assignment is the block one.
- **Block dominance, verified** with full rows `(2,1),(3,1),(4,1),(3,2),(4,2),(4,3)` ("set 1"):
  `det M'(eps) = det[6 full rows x cols 3..8] . det[2 good rows x cols 1..2] . (1 + O(eps))`, and the
  2x2 tends to `N'`. Ratios at eps = 1e-3, 1e-4:

  | state | ratio |
  |---|---|
  | `{2,2,0,0}` | 1.0013, 1.00013 |
  | `{2,1,1,0}` | 1.0029, 1.00029 |
  | `{1,1,1,1}` | 1.0060, 1.00060 |

  Two other row sets for `{2,1,1,0}` agree. With these raw rows, the order of `det M'` is 0, 2, 8,
  which equals `8 - sum_{i<j}(mu_i - mu_j)`: observed only.

### 3.4 Does the 6x6 block simplify? (open)

- **Singlet `{1,1,1,1}`: yes.** Over ten generic rays the leading coefficient is exactly
  `m6 = -(1/2) (phi1-phi2)^2 (phi1-phi3) (phi2-phi3)^2 (phi2-phi4) (phi3-phi4)^2`. Exponents are
  integers to 1e-7, and `K = -1/2` on every ray.
- **`{2,1,1,0}`, `{2,2,0,0}` with set 1: no simple form.** The leading coefficient is:
  - not a monomial in `(phi_i - phi_j)`;
  - still not a monomial after normalising the raw Q's by their most divergent coefficient ("option
    1"; singlet control unchanged);
  - not a polynomial;
  - not a twist function common to all states (state ratios are not monomials).

  It behaves like a rational function with zeros and poles where twists coincide.
- **Row choices with product form exist but lose dominance.** Scanning all C(16,6) = 8008 full-row
  sets finds pure products (321, 282, 5176 sets for the three states; 82 common to all three). But for
  these, `det M'/(6x6 x 2x2)` is ≈ 87 or ≈ 0.005, not 1: they reach product form by cancelling the
  leading 6x6 terms. **Product form and block dominance did not coexist for the non-singlets.**
  Rational constants I first read off these products (`1/90720`, `1/20160`) were **retracted**: they
  are not distinguishable from arbitrary fractions at 1e-7.

### 3.5 Exact route for `(1,1) {2,1,1,0}` (in progress when stopped)

- Exact twisted solution (TASK file §8.2), rational after eliminating the root. Take `(z1, z2, r)` as
  parameters, with `z3 = z2 (r-theta1-1)(r-theta2-1)/((r-theta1)(r-theta2))` and `z4 = 1/(z1 z2 z3)`.
  Then `q1 = 1`, `q2 = u - r + z1/(z1-z2)`; `q3`, `q4` and all four duals follow by linear solves or
  exact Wronskians.
  - All four Baxter residuals vanish identically.
  - Our `q^3` equals the TASK formula with `u -> u-1` (centring convention).
  - `q3 = u - r' + z1/(z1-z3)`, with `r'` the root of the relation with `z2 <-> z3`.
- Exact leading coefficient of the set-1 6x6 on the ray `z1 = 1+eps`, `z2 = 1+eps/3`,
  `r = 31/42 + eps/5` (`r0 = 31/42`): a rational, `C6 = 140.28972856788119555703...`, power 2. It is
  identical at working orders 16 and 20 and matches the numerics to 1.75e-8. Fast method: log-series
  of `z_j` once, `z^w = exp(w log z)`, truncated polynomial products, row scaling by eps-powers before
  `Det` (`exRayC6fast`).
- Useful fact: the first-order ray has `a3 = a2 - L rho/((r0-theta1)(r0-theta2))`, and
  `(r0-theta1)(r0-theta2) = Qbar12(theta1+1) Qbar12(theta2+1) = (25/42)(17/42)`, with
  `Qbar12(u) = u - r0`. **`Q12(theta)`-type values enter naturally here** (Paul's hint).
- On the slice `a1 = 1, a3 = -1/2, a2 = x`: zeros at `x = 1` (order 1), `-1/2` (order 2) and `0`
  (order 1); a pole of order 2 at `x = -1/4`; the eps-power changes near `x = -3/2` (`a1 = a4`).
  - 16 exact points: **no** rational function with ≤14 unknowns fits.
  - A ~70-point run with fixed denominator `prod_p (x-p)^10` was **still computing when the session
    ended** (cell 52, variables `exPts2`, `polyRecon`). Its result is unknown. The kernel may still hold
    it if it wasn't restarted.
- Symbolic-parameter exact runs (one twist parameter left symbolic) are too slow: over 800 s. Don't
  retry naively.

### 3.6 Suggested next steps

1. Finish the slice reconstruction (3.5): exact points are cheap (~4 s each), and reconstruction can
   be done by fixed-denominator interpolation or Thiele. Then factor and look for `Qbar12(theta+nh)`
   factors. Then a second slice to get the two-variable form.
2. Try the full-contour rows with **finite** Q's (Wronskian-type `Q1, Q12, Q123` and duals) instead
   of raw ones. Raw-Q divergence is the leading suspect for the non-product structure, though "option
   1" did not fix it.
3. Relate `M'` to the Clean 8x8 (same identity, different rows and columns) before claiming anything
   about "the" twisted norm.
4. Explain the vacuous pairs (1.5) and the `4L-6` count (Part 2); write the L=3 IoM relations.

---

## Part 4 — Pitfalls hit this session (read before touching the notebooks)

- **The laptop crashed once.** 20 exact twisted solves at tiny eps plus exact `Cancel/Together` on
  algebraic numbers in the Hodge duals ran out of memory. Wrap heavy cells in
  `MemoryConstrained[TimeConstrained[...]]`. Solve twisted states **numerically** (`NSolve` at
  60-100 digits) or exactly **only after eliminating the root** (rational parametrisation).
- **Never use `ν` (or any global function name) as a loop or pattern variable.** `Do`/`Table`
  temporarily reassign their iterator, which clobbered `ν[S,P][i][u]` and poisoned the memoised
  `τ3quad`. This made non-dominant weights look unsolvable for two turns. Clear with
  `DownValues[τ3quad] = Select[DownValues[τ3quad], !FreeQ[#[[1]], Pattern] &]`.
- `Table[..., {k, f}, {{1,q1},{2,q2}}]` is **not** a loop over pairs; use `{kf, {{1,q1},{2,q2}}}`.
- **Zero tests on tiny numbers:** use `== 0 || Precision[#] < 10`, not an absolute threshold. The
  singlet's block is ~1e-40 at eps = 1e-5.
- **Wolfbook quirks.**
  - Tool replies truncate at `{`: store results in variables and read small pieces.
  - Don't export to scratch `.txt` (Paul wants everything in notebook cells).
  - Insert-and-evaluate can report "completed" while the cell is still running; check
    `list_clients` for `active_operation`.
  - `newNotebook(overwrite=true)` on a notebook **open in VS Code** does not overwrite (the editor
    buffer wins). Delete and insert cells through the editor instead.
- **Twisted scalar product requires a unimodular twist**, and the contour must enclose all poles
  (see Clean `.tex`).
- `Rationalize` at a 1e-7 tolerance always "succeeds" with large denominators. Don't read physics
  into it.

---

## Part 5 — Translation to the conventions of arXiv:2603.26617

The construction has been repeated directly in the paper's conventions in
`Paul/Mathematica/XXX/Experiments/SU4_SoV_Orthogonality_Untwisted_L2_2603.wb` for the paper's
two-site chain with local representation `[1,1,0,0]`. The notebook runs cold from top to bottom.

The dictionary is `h = i`, `f^[a](u) = f(u + a i/2)`, and a recentering of the contour variable
`u_old = u + i`. Thus, with `Qtheta(u) = (u-theta1)(u-theta2)`,

    E2_old(u+i) = 1/(Qtheta^[-2](u) Qtheta(u) Qtheta^[2](u)),

which is exactly the measure in equation (E1). The transfer matrices are identified as

    tau+(u) = T1(u),    tau0(u) = T2(u+i/2),    tau-(u) = t3(u+i),

where `T3 = Qtheta(u) Qtheta(u-3i) t3(u)`. In these variables the two surviving differences are

    Delta tau+ = Delta tau- = Delta t,
    Delta tau0^+ = Delta t (Qtheta^[-2] + Qtheta^[2]) + Delta c',

with `Delta t = -(1/2) Delta C2` and
`Delta c' = (1/8) Delta(C2^2) - (1/4) Delta C4 + (3/4) Delta C2`.

Substitution in (E3) gives the two columns

    F1 = Q_A^j [Qtheta^[-2](Q_i^{B[2]}-Q_i^B)
                 + Qtheta^[2](Q_i^{B[-2]}-Q_i^B)],
    F2 = -Q_A^j Q_i^B.

For `(i,j)=(2,1)`, the notebook verifies exactly all direct and dual Baxter equations, all Casimir
relations, three nonzero norms, all six cross determinants, and the kernel identity for every
ordered pair of the three irreps. The paper's two rows with insertions `1,e^(2 pi u)` are a
state-independent row transformation of the two pole-ladder rows, with determinant
`exp(2 pi theta2)-exp(2 pi theta1)`.

Easy-to-miss paper detail: (C1) contains `tau-^[2]`, while (C15) contains `tau+^[2]`. Dropping these
shifts makes the polynomial Baxter system inconsistent.

---

## Part 6 — Relation to transfer-matrix eigenvector norms (2026-09-12)

Implemented in `SU4_SoV_Orthogonality_Untwisted_L2.wb`, cells 18-21. The notebook now builds the
GT-basis transfer matrices explicitly and constructs, for every benchmark irrep, a right eigenvector
`PsiR` and a transpose-left eigenvector `PsiL`. Since zero-twist transfer matrices label a whole
`gl(4)` multiplet, the calculation is restricted to the total-weight sector equal to the irrep's
highest weight. The target eigenline is one-dimensional there. Both vectors are put in **pivot
gauge**: their first nonzero component in the ordered tensor-product GT basis is one.

The exact bilinear norms `PsiL.PsiR` are:

| `(S,P)` | irrep | `PsiL.PsiR` |
|---|---|---:|
| `(1,1)` | `{2,2,0,0}` | `1` |
| `(1,1)` | `{2,1,1,0}` | `2` |
| `(1,1)` | `{1,1,1,1}` | `6` |
| `(2,1)` | `{3,3,0,0}` | `1` |
| `(2,1)` | `{3,2,1,0}` | `3/2` |
| `(2,1)` | `{2,2,1,1}` | `10/3` |
| `(2,2)` | `{4,4,0,0}` | `1` |
| `(2,2)` | `{4,3,1,0}` | `2` |
| `(2,2)` | `{4,2,2,0}` | `3` |
| `(2,2)` | `{3,3,1,1}` | `5` |
| `(2,2)` | `{3,2,2,1}` | `8` |
| `(2,2)` | `{2,2,2,2}` | `20` |

Every vector satisfies the exact right/transpose-left eigenvalue equations for both independent
zero-twist transfer matrices `tau1` and `tau2`. The relation to the SoV determinant is

    det N^{AA} = C_A (PsiL.PsiR),
    C_A = det N^{AA}/(PsiL.PsiR).

The notebook tabulates `C_A` exactly for all 12 irreps. It is state-dependent, so the SoV
determinant is not the honest norm in the pivot gauge up to one universal chain constant. This is a
normalisation statement, not a contradiction: `PsiR -> a PsiR` and `PsiL -> b PsiL` multiply the
honest bilinear norm by `a b`, while the monic/gauge-fixed Q-functions impose a different state
normalisation. A closed formula for `C_A`, or an explicit map between the Q-gauge and Hilbert-space
state gauge, remains open.

### 6.1 Vacuum normalisation of the SoV determinant

The same notebook now defines, for the default good pairing `(i,j) = (2,1)`,

    Nvac[S,P] = det N^{vac,vac},
    mu_vac = {S+P,S+P,0,0},
    Nnormalized[S,P,A] = det N^{AA}/Nvac[S,P].

The exact vacuum denominators at the benchmark point are

| `(S,P)` | `Nvac[S,P]` |
|---|---:|
| `(1,1)` | `-145570417/111435` |
| `(2,1)` | `-5960610348/497743` |
| `(2,2)` | `-1183291756455123/14683418500` |

The notebook tabulates all 12 normalized norms in `normalizedNormData`; the field is
`"normalizedSoVNorm"`. All three vacuum entries are verified exactly to equal one. Optional `i,j`
arguments on `Nvac` and `Nnormalized` allow the same convention to be applied to other good
Q-pairings.

With the Hilbert-space comparison included, the normalized relation is

    Nnormalized_A = (C_A/C_vac)
                    (PsiL_A.PsiR_A)/(PsiL_vac.PsiR_vac).

The pivot-gauge vacuum honest norm is one for all three benchmark chains.

### 6.2 All allowed Q-pairings

The notebook now defines

    allowedPairings = {{2,1},{2,2},{3,1},{3,2},{3,3},{4,1},{4,2},{4,3}}

and stores the complete `8 x 12 = 96` raw/vacuum-normalized ledger in `pairingNormData`. Each record
contains `"pairing"`, `"SP"`, `"mu"`, `"NAA"`, `"Nvac"`, and `"normalizedNorm"`. Every one of the
96 diagonal determinants is nonzero, and all 24 pairing-specific vacuum ratios equal one.

Important result: vacuum normalization does **not** remove pairing dependence. The three vacuum
states agree across pairings trivially; every one of the nine non-vacuum states has eight distinct
normalized values. Hence a Q-pairing, or a further pairing-dependent conversion to Hilbert-space
normalization, remains part of the norm convention.

### 6.3 Symbolic first residue-pair term

Cell 29 of the orthogonality notebook leaves the Q-functions completely symbolic. Write
`Qd(u) = q_A^j(u)`, `Q(u) = q_i^A(u)`, `Delta Q(u) = Q(u)-Q(u-h)`, and
`delta = theta1-theta2`. The determinant has the residue-pair expansion

    det N = Sum_{k=0}^{S+1} Sum_{l=0}^{P+1} det({r_{1,k},r_{2,l}}).

The first term is `D00 = det({r_{1,0},r_{2,0}})`, with poles at `theta1,theta2`. Define

    rho1 = (-1)^(S+1)/(h^(S+1) (S+1)! Product_{m=0}^{P+1}(delta-m h)),
    rho2 = (-1)^(P+1)/(h^(P+1) (P+1)! Product_{m=0}^{S+1}(theta2-theta1-m h)).

Then

    r10 = rho1 Qd(theta1) {
            (S+1) h (delta-(P+1)h) Delta Q(theta1),
            Q(theta1-h)},

    r20 = rho2 Qd(theta2) {
           -(P+1) h (delta+(S+1)h) Delta Q(theta2),
            Q(theta2-h)}.

Consequently

    D00 = rho1 rho2 Qd(theta1) Qd(theta2) h [
            (S+1)(delta-(P+1)h) Delta Q(theta1) Q(theta2-h)
          + (P+1)(delta+(S+1)h) Q(theta1-h) Delta Q(theta2)].

Thus the dual Q is evaluated unshifted at the first pole of each ladder. The direct Q occurs at a
backward shift and through a backward finite difference. This formula is universal for all eight
allowed pairings; only the substitution `Qd=q^j`, `Q=q_i` changes.

### 6.4 Normalisation by the first residue-pair

For every allowed pairing the notebook now defines

    ND00normalized[S,P,state,i,j] = NAA[S,P,state,i,j]/D00[S,P,state,i,j],
    D00remainderRatio[S,P,state,i,j] = (NAA-D00)/D00.

Hence, whenever `D00 != 0`,

    ND00normalized = 1 + D00remainderRatio
                     = 1 + Sum_{(k,l)!=(0,0)} Dkl/D00.

The complete ledger is stored in `d00NormData`. At the benchmark point `D00` is nonzero for all
96 allowed pairing/state records, and the `1 + remainderRatio` identity is verified exactly in all
96 cases. This remains a conditional normalization at generic parameters: special loci can make
`D00` vanish even if the full determinant is nonzero.

### 6.5 Vacuum-normalised first-residue ratio versus the honest norm

The notebook now defines the double ratio

    ND00vacuumNormalized[S,P,A,i,j]
      = (NAA_A/D00_A)/(Nvac/D00_vac).

This is the consistent way to retain the state-specific `D00` normalisation while fixing the
vacuum value to one. Dividing every raw determinant by the fixed number `D00_vac` and then
vacuum-normalising would instead cancel `D00_vac` and simply reproduce `NAA_A/Nvac`.

Let

    H_A = PsiL_A.PsiR_A,
    C_{A;i,j} = NAA_{A;i,j}/H_A.

Then the exact relation is

    ND00vacuumNormalized_A
      = Gamma_{A;i,j} (H_A/H_vac),

    Gamma_{A;i,j}
      = (C_{A;i,j}/C_{vac;i,j}) (D00_{vac;i,j}/D00_{A;i,j}).

Equivalently,

    H_A/H_vac = ND00vacuumNormalized_A/Gamma_{A;i,j}.

The full 96-record ledger is stored in `d00VacuumNormData`. All 24 vacuum records equal one, and
all 96 conversion identities above are verified exactly using the independently computed raw
SoV-to-honest factors. But `Gamma=1` only for the vacuum records: for each non-vacuum state the
eight allowed pairings give eight distinct double ratios and eight distinct conversion factors.
Thus the `D00` double ratio is not yet the honest norm. It fixes the vacuum and isolates a residual
state- and pairing-dependent Q-gauge conversion factor. In the present transfer-matrix pivot gauge
`H_vac=1`, so the honest norm is `H_A = ND00vacuumNormalized_A/Gamma_{A;i,j}`.

### 6.6 Superseded per-ladder analytic examples (no D00)

**Not the target norm construction.** The formulas in this subsection came from assigning the two
determinant rows to different inhomogeneity ladders. They are retained only as a diagnostic. The
correct construction, in which both rows use the same full contour and differ only by their
Q-functions, is in Section 6.8 below.

The notebook now contains a symbolic Baxter solve with `theta1`, `theta2`, and `h` left unfixed for
the full `(S,P)=(1,1)` multiplet. To keep the formulas compact it uses the allowed pairing
`(i,j)=(3,3)`, i.e. direct `q_3` and dual `q^3`, and defines

    sigma = theta1 + theta2,
    delta = theta1 - theta2,
    D = 3 delta^4 - 6 delta^2 sigma^2 - 12 h^2 delta^2
        - 4 h^2 sigma^2 + 3 sigma^4.

For the simple vacuum normalization `Nhat_A=NAA_A/Nvac`, without any `D00` division, the exact
results are

    Nhat_[2,2,0,0] = 1,

    Nhat_[2,1,1,0]
      = -sigma (delta^2-sigma^2)
          (delta^2+8 h^2+12 h sigma+3 sigma^2)
        /((2 h+sigma) D),

    Nhat_[1,1,1,1]
      = (4 h^2-delta^2) (4 h^2-3 delta^2+3 sigma^2)/(3 D).

At the benchmark `theta1=1/3`, `theta2=1/7`, `h=1`, these specialize to
`1`, `-7955/18564`, and `-1748/357`. All three agree exactly with the pre-existing numerical
calculation for the same `(3,3)` pairing. These are still Q-gauge determinant norms rather than
the pivot-gauge transfer-matrix norms.

### 6.7 Analytic audit of the transpose-left/right norm

The earlier transfer-matrix calculation had `theta1=1/3`, `theta2=1/7`, and `h=1`. The notebook
now reconstructs the GT transfer matrices for the `(S,P)=(1,1)` chain with all three parameters
symbolic, using

    L_lambda(u-theta_s)_{ij}
      = (u-theta_s) delta_ij 1 - h J_{ji}^{(s)}

and quantum-minor shifts by `h`. It uses the fully analytic Baxter eigenvalues, restricts to the
highest-weight sector of each irrep, solves the right and transpose-left null spaces, and puts both
vectors independently in first-nonzero-component-equals-one pivot gauge.

The exact results are

| `mu` | `PsiR` in its highest-weight sector | `PsiL` in the same sector | `PsiL.PsiR` |
|---|---|---|---:|
| `[2,2,0,0]` | `(1)` | `(1)` | `1` |
| `[2,1,1,0]` | `(1,-1)` | `(1,-1)` | `2` |
| `[1,1,1,1]` | `(1,-2,2,2,-2,1)` | `(1,-1/2,1/2,1/2,-1/2,1)` | `6` |

The vectors and norms contain no `theta1`, `theta2`, or `h`. For all three states, both the right
and transpose-left equations are verified symbolically for both `tau1` and `tau2` (12 exact
checks), using the same analytic Baxter eigenvalues. Specializing the analytic vectors reproduces
the earlier benchmark vectors and norms exactly.

Therefore `1,2,6` are correct for the explicitly stated independent pivot gauges. They are not
gauge-independent absolute Hilbert norms: `PsiR -> a PsiR`, `PsiL -> b PsiL` changes the bilinear
norm by `a b`. The analytic result verifies the convention actually defined here—transpose-left,
right, independent pivot normalization, and the bilinear pairing `PsiL.PsiR`.

### 6.8 Transfer-first Q solve and corrected full-contour norm

The exact repository solver was located in
`Paul/Mathematica/XXX/Experiments/baxter_SU4_L2_XXX.wb`. Its `QsolveG` fixes the exact untwisted
transfer eigenvalues first and then solves a linear, monic, staggered-gauge Baxter problem. The norm
notebook now implements the same fixed-eigenvalue architecture locally, with an additional linear
solver for the dual Q-functions.

For `(S,P)=(1,1)` the independent pipeline is:

1. Find the global highest-weight kernel at each `mu`, without using Baxter data. Every kernel has
   dimension one.
2. Read the exact symbolic `tau1(u),tau2(u),tau3(u)` eigenvalues from the transfer-matrix action.
3. Solve the four direct and four dual Baxter equations linearly. Each of the eight solves has
   exactly one monic/staggered-gauge solution for every irrep.
4. Check that the resulting Q-polynomials agree exactly with the previous formulas.

The corrected determinant has two full-contour rows

    row 1: (q_3, q^3), summed over every pole on both ladders,
    row 2: (q_4, q^4), summed over every pole on both ladders.

The rows are distinguished only by their Q-functions. For the three `(1,1)` irreps the exact
vacuum-normalised results are

| `mu` | full-contour `NAA/Nvac` | `PsiL.PsiR` | conversion |
|---|---:|---:|---:|
| `[2,2,0,0]` | `1` | `1` | `1` |
| `[2,1,1,0]` | `3/2` | `2` | `4/3` |
| `[1,1,1,1]` | `1/2` | `6` | `12` |

All three full-contour ratios are exactly independent of `theta1`, `theta2`, and `h`. Therefore the
inhomogeneity dependence in Section 6.6 was an artifact of the wrong per-ladder row construction;
no special-point Q-product is required to cancel it. A Q-normalisation formula may still explain
the remaining state-dependent constants `1,4/3,12`, but that is now a parameter-independent
normalisation problem.

### 6.9 Extended representations and the universal full-contour factor

The corrected calculation was extended exactly to

    (S,P) = (1,1), (2,1), (2,2), (3,1), (3,2), (4,1), (4,2), (3,3), (4,3),

covering 50 irreps.  With `L=S+P` and

    mu = (L-j, L-k, k, j),   0 <= j <= k <= min(S,P),

every state satisfies

    NAA_mu/NAA_vac = (k-j+1) (L+2-j-k)/(L+2).

All 50 checks are exact.  Thus the vacuum-normalised full-contour SoV norm depends only on `L,j,k`,
not on the split `(S,P)`, and contains no inhomogeneity or `hbar` dependence.  This sharply reduces
the comparison problem: special-point Q-values are not needed to remove theta dependence.

The transpose-transfer Hilbert norm still depends on the site representations.  In the current
independent pivot gauges the exact relation is

    PsiL.PsiR = F_{S,P;j,k} (k-j+1) (L+2-j-k)/(L+2),

where the notebook tabulates `F=(PsiL.PsiR)/(NAA/Nvac)` for all 50 states.  Three nontrivial
Hilbert-norm subfamilies already have simple combinatorial forms:

    j=0:       PsiL.PsiR = Binomial[L-k+1,k]/Binomial[S,k],
    j=k=1:     PsiL.PsiR = Binomial[L+2,2]/Binomial[S+1,2],
    j=k=P:     PsiL.PsiR = dim[S,S]/dim[S-P,S-P],

with

    dim[n,n] = (n+1)(n+2)^2(n+3)/12.

Every applicable row in the extended data set passes these formulas exactly.  The remaining
interior families, beginning with `j=1,k=2`, are the right place to test a general hook/factorial
formula or, only if necessary, a compact monic-Q endpoint expression for the residual gauge factor.

### 6.10 Comparing Q choices and the first special-point conversion formula

In the corrected construction a determinant is labelled by two full-contour rows.  Write

    Nhat_{r1,r2} = det{R_r1,R_r2}_state / det{R_r1,R_r2}_vac,
    r=(i,j) meaning the row built from (q_i,q^j).

The four universally nonzero choices made only from diagonal rows are `(1,2)`, `(1,3)`, `(2,4)`,
and `(3,4)`, where `(a,b)` abbreviates the row pair `{(a,a),(b,b)}`.  On all 50 exact states,

    Nhat_12 = Nhat_13 = Nhat_24 = Nhat_34.

Thus every ratio between these choices is exactly one.  Before vacuum normalization the
corresponding determinant ratios are state independent and equal their vacuum ratios.  Generic
mixed row pairs do not share this property: among the 41 universally nonzero pairs, the original
non-vacuum benchmark states give 24 or 27 distinct normalized values.

A natural nested-Q object gives the Hilbert normalization on the complete 50-state data set.  Define
the monic direct Casoratian

    Qbar12(u) = monic[q1(u) q2(u-h) - q1(u-h) q2(u)],
    d = deg Qbar12 = j+k,

and

    Q12hat_top1 = Qbar12(theta1+S h) /
                  (h^(j+k) (S-k+1+t)_k (S-j+2+t)_j),
    t = (theta1-theta2)/h.

Then all 50 tested states satisfy exactly

    PsiL.PsiR = 1/Q12hat_top1.

Consequently their conversion factor is

    (PsiL.PsiR)/(NAA/Nvac)
      = (L+2)/[(k-j+1)(L+2-j-k) Q12hat_top1].

This is independently validated, not merely pair-fitted.  Two simpler normalizations were retained
as useful rejected candidates in the notebook ledger: `(S+t)_(j+k)` fails 23 states, while
`(S-k+1+t)_(j+k)` fails 9 states with `j>=2`.  Splitting the kinematic denominator into the two
content strings `(S-k+1+t)_k (S-j+2+t)_j` removes every residual.  An all-spin analytic proof remains
open, but no additional Q-minor is required by the current exact data.

### 6.11 First explicit `Q12hat_top1` table

Section 11 of the notebook now tabulates all 21 states in
`(S,P)=(1,1),(2,1),(2,2),(3,1),(3,2)`, including `mu`, `(j,k)`,
`Q12hat_top1`, `1/(PsiL.PsiR)`, and an exact equality check.  Every check is true.
