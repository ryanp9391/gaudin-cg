# sl(3) Gaudin, round 5 plan: from "the adjoint works" to "the formula" — master kernel, general-Λ limit, the mechanism of the ε-series, the norm, and the write-up (2026-10-16)

Status: DONE 2026-10-16 ~23:00 — all eight SONNET agents returned (Paul: "read it and implement with sonnet subagents"). Agents and their reports: M1 `Paul/Python/gaudin_sl3/RESULTS-master-kernel.md` (A1+A3: no four-index closed form; orders t⁰–t² closed for all λ₁ in symbol form, [5,0,0]-confirmed; Casimir at order 1; translation covariance), M2 `RESULTS-jet-basis-sym.md` (A2: the jet basis DERIVED as an operator identity ⟨i,j| = ⟨x₀|ν_i(Ĥ,b̂)ν_j(Ĥ,b̂); master kernels λ₁ = 3, 4 derived, [5,0,0] available), M3 `RESULTS-master-kernel-2.md` (follow-up: [5,0,0] acceptance passed (6e-56 on 441 states), exact [5,0,0]/[6,0,0] kernels, one universal matrix N for all λ, the kernel = normal-ordered exponential of a universal bidifferential operator, closed through t⁵, Hamilton–Jacobi tree recursion + one-loop terms; the single-second-operator span test fails at λ = 3), N1 `RESULTS-nu-general.md` (B1+B2: ν-reduction for general Λ, constant 1 on all 20 branches, rules k_B = 2(λ₁−λ₃) − (μ₁−μ₃), N_B = 4λ₁ − 2λ₂ − k_B, all 20 CG tables), W5 `Experiments/Gaudin-Laurent-series.wb` + `ref/laurent-series-{210,310,320}.txt` (B3+C3+C2: exact series to ε¹⁰, H₁ exactly linear in ζ, χ_B = char poly of a linear pencil, d_B = weight-space dimension except the adjoint singlet, the 2×2 closed form, the duality map), S1 `RESULTS-bethe-series.md` + `cg_from_Q.py` (C1+F3: the Bethe equations alone determine every series, 26/26 branches = W5 exactly; the end-to-end script, 3 s, ~1e-91), Q1 `RESULTS-gaudin-norm.md` and Q2 `RESULTS-gaudin-norm-general.md` (D1+D2 and general Λ: N = (−1)^{μ₂}C(Λ)det J₆Π(ζ_a−ζ_b)^{μ_c−2|Λ|}, C(Λ) = [λ₂!(λ₁−λ₂)!(λ₁+1)!/((λ₂+1)(λ₁−λ₂+1))]², sign rule, 15 reps; the SoV measure Γ and completeness). Decisions taken by the session: A and C/D in parallel; [4,2,0] backlog; paper (E2) after this round — see HANDOFF.md for the per-agent result blocks and the decision on `paper/gaudinCG.tex` still open for Paul.
Project target: Clebsch–Gordan coefficients of gl(3) from separation of variables in the two-site Gaudin model. Rounds 1–4
(`Gaudin-sl3-CG-plan.md`, `-round2-`, `-round3-`, `-round4-plan.md`) are summarised in HANDOFF.md and, once written, in the notes.

## 0. Where we are (one paragraph per pillar; everything here is verified by ≥ 2 agents unless marked)

**P1 — the overlap is a function of the eigenvalue data.** For the CG-complete two-vacua family Φ∞(ω₁,ω₁′) and any eigenstate Ψ with
Laurent data (H_s, a_s, b_s) of S₂, S₃ at θ_s: c·⟨Φ∞|Ψ⟩/⟨x₀^G|Ψ⟩ = Σ_ij G_ij(ω₁,ω₁′,θ₁₂) ν_i(θ₁)ν_j(θ₂), with ν rational in the Laurent data
(adjoint: ν = (1, κ, π, κπ, ε_s), κ = a_s/3, π = 2H_s − a_s − 2Σζ − 6/δ, ε_s explicit; b_s = a_s(H_s − a_s)/3). The Q polynomials never enter.
Proved by exact algebra for the adjoint (U2, U1); the mechanism (Frobenius at θ_s, exponents {0, λ₂+1, λ₁+2}, all Plücker coordinates ∝ one
factor) is rep-generic but the ν-formulas for λ₁ ≥ 3 are not yet written.

**P2 — the kernel.** G = (ω₁ω₁′)^{2λ₁−λ₂}Θ, Θ polynomial in (a, a′, t) = (1/ω₁, 1/ω₁′, 1/θ₁₂); Θ(t = 0) = θ(a)⊗θ(a′) with the one-site kernel in
closed form for every Λ (θ_kl(a) = c_kl a^{k+l−λ₂−1}, c_kl = (−1)^{k+l−λ₂−1}C(λ₁−λ₂, l−λ₂−1)C(λ₁+1,k)(l−k)/((λ₂+1)C(l,λ₂+1)); a dressed Wronskian);
BLOCK THEOREM Θ^{[λ₁,λ₂,0]} = [λ₂!(λ₁−λ₂)!/λ₁!]²(aa′)^{−λ₂}Θ^{[λ₁,0,0]}|_{non-forced pairs}, i.e. Θ = θ_i(a)θ_j(a′)R^{λ₁}_ij(at, a′t) with R universal per
λ₁; R closed to first order in t for all λ₁ (G1's four-index formula) and to ALL orders for the single-index reps [l,l,0] (Γ-form; l = 5 predicted
and confirmed); top-state sum rule Σ Θ m^{top}m^{top} = MM′R(a)R(a′), R(a) = (1−ζ₂a)^{λ₁−λ₂}(1−ζ₃a)^{λ₁}; λ₁ = 1: Θ = (1+z)θ⊗θ′ − zθ′⊗θ, z = t/(ω−ω′).
Structural origin (D3): Θ is the Gaudin K tensor — Φ_c∞ expanded in the canonical, ζ-independent, rational jet basis Y defined by the SoV
wavefunction rule ⟨i,j|Ψ_c⟩/⟨x₀|Ψ_c⟩ = m_i(θ₁)m_j(θ₂)/(M₁M₂). OPEN: the general four-index master kernel R^{λ₁} beyond first order.

**P3 — the untwisting limit.** Σ G_ij(εω̂)ν_iν_j = ε^{6−k_B}·ĉ·⟨φ_u⊗e₁|Ψ₀⟩/X_B + O(ε^{7−k_B}) as a polynomial identity in (ω̂₁, ω̂₁′), k_B = 0,1,1,2,2,4 for the
six adjoint branches, φ_u = e₁ᵀρ(n′) with ε-free Gauss parameters, ĉ = c/ε⁸, X_B = ⟨x(φ)⊗x(φ)|Ψ_{k_B}⟩ in closed form (W4). Every CG coefficient
of adj⊗adj follows, φ-independent, to 77 digits. Needs the ε-series of (H_s, a_s) to order ≤ 5; these are exact and polynomial in φ (U1, W4):
top linear, (4,1,1)/(3,3,0) algebraic (√(1 + ε²(φ_i−φ_j)²/4)), singlet to n = 6, octets to n = 5 over Q(√5). Untwisted facts: the ζ = 0 oper is
Fuchsian at ∞ with three polynomial solutions of degrees μ + (2,1,0); ₃F₂ iff a consecutive exponent pair; a₁ = −a₂ = (μ₁+2)(μ₂+1)μ₃/(2θ₁₂);
b_s ≠ 0 only in multiplicity-2 sectors (one exceptional sector, [3,2,0] (5,3,2)). OPEN: general Λ; a derivation (not fit) of the ε-series; the
normalisation of Ψ₀ from the Q side.

**P4 — what else is in hand.** Localisers (ℓ₁ = λ₁+1; every D₁ Laurent coefficient localises; X₄ for λ₁ = 3), det N for equal twists with
Lagrange columns (all reps), no-log conditions (Gaudin Bethe equations), the Gaudin B hierarchy (b₁ Cartan, b₂ with cross-site terms, solvable
per level block, no single polynomial B_G), no cleaning anywhere, universal multiplicity-2 mixing matrices, duality [3,2,0] ≅ [3,1,0]*⊗det³ acting
literally on the data.

## 1. Round-5 work packages (ordered by value; each is self-contained and sonnet-sized)

### WP-A. The master kernel R^{λ₁} in closed form (the general formula)

Goal: a four-index closed form R^{λ₁}_{(kl),(k′l′)}(x, y), x = at, y = a′t, valid for all λ₁ (then every Λ by the block theorem and θ ⊗ θ′).
Data in hand: exact R for λ₁ ≤ 4 as polynomials (G1's `general_kernel_*.pkl`; [4,0,0] 225 entries), the first-order formula, the single-index
Γ-form (which is R restricted to pairs (k, λ₁+1)), the extreme-hop coefficients (1/m!)(n−k)_m(n−l)_m(k′)_m(l′)_m/((n)_m)², the exchange rule
R_ij(x,y) = R_ji(y,x) with t → −t, homogeneity, the rank pattern 1,2,4,8 of the t-expansion, and the top-state sum rule.
- A1 (fit-guided ansatz, Python exact): generalise the Γ-form: R_ij = Σ_{m,m′} F_{mm′}(k,l;k′,l′) x^m y^{m′} with F a ratio of products of
  Pochhammer symbols in the six integers (n, k, l, k′, l′, m, m′); constrain by the known m = 0, m′ = 0 rows (one-site), the first-order row, the
  extreme hops, the single-index block; fit the remaining structure on [3,0,0], [4,0,0]; validate on a NEW rep ([5,0,0]: 1225 states — G1 found
  the build too slow at 80 digits; do it at 40 digits with the ν-form, which needs only the Laurent data of the 1225 states and no Q solving —
  this is where P1 pays off computationally).
- A2 (structural derivation, Python exact): the jet basis Y for one-row reps. V_{(λ₁,0,0)} = Sym^{λ₁}(C³) has a polynomial realisation (E_ij =
  x_i∂_j on degree-λ₁ polynomials); the two-site space is bihomogeneous polynomials; the SoV wavefunction rule and D3's construction of Y should
  become explicit differential/integral operators in this realisation. Deliverable: Y^{[λ₁,0,0]} in closed form (or at least its 1/θ₁₂ structure),
  hence R = (Φ_c∞ in the Y basis) — a derivation rather than a fit. Start with λ₁ = 1, 2 where Y is small (9×9, 36×36) and the R-matrix form is known.
- A3 (the auxiliary sl(2)): G1 observed that the jet index space carries an sl(2) of spin (λ₁+1)/2 under which θ^{(a+s)} is the Möbius translate
  of θ^{(a)}; the two-site kernel is NOT covariant under a common shift. Ask: is R covariant under a TWISTED action (shift accompanied by t → t′)?
  If so, R is fixed by its value at a = 0 or a′ = 0, i.e. by the one-site data plus a finite number of coefficients — the simplest conceivable
  closed form. Also test whether R^{λ₁} is the spin-(λ₁+1)/2 fusion of the λ₁ = 1 R-matrix form (G1 refuted one "box-exchange" version; the twisted
  action may be what was missing).
- Acceptance: R^{λ₁} reproduces all exact kernels for λ₁ ≤ 4 and predicts [5,0,0] (and [5,1,0] via the block theorem) on ≥ 200 states to ≥ 30 digits.

### WP-B. The ν-formulas for general Λ, and the untwisting limit for [3,1,0] and [3,2,0]

- B1 (algebra, sympy exact): Frobenius at θ_s with exponents {0, λ₂+1, λ₁+2}: the order-(λ₁+1) jet of every solution lies in a 2-plane
  spanned by the exponent-0 and exponent-(λ₂+1) Frobenius vectors, whose coefficients are rational in the Laurent data (H_s, a_s, b_s, the other
  site's data through the regular parts, e_k(ζ)); hence every non-forced Plücker coordinate p_kl = D_s·ν_kl with ν_kl rational and D_s = the
  2×2 determinant of the (P², P³) components. Derive ν_kl for general (λ₁, λ₂) (symbolic recursion; verify on P3/P4's 225 states for [3,1,0],
  [3,2,0] to 1e-50), and prove that P³·D_s cancels between the left and right sides (M = P³p_{0,λ₂+1} = P³D_s·ν_{0,λ₂+1}): the general statement
  "c·Qnum = Σ G νν depends only on the Laurent data". Also the general no-log conditions in closed form (b_s = a_s(H_s − a_s)/c₂ generalises? P3/P4
  have them as polynomials).
- B2 (numerics, Python 200 digits): U2's pipeline for [3,1,0] (G1's new polynomial laws give the ω-degree split) and [3,2,0] (P4's grid laws):
  k_B per branch (10 each), the leading polynomial identity with the exact CG tables (W4's data files), the depth per branch. Expect k_B to follow
  a weight rule (adjoint: k_B = 0,1,1,2,2,4 — is it the "depth" of μ below the top weight in some grading? test).
- B3 (exact series): the ε-series of the Laurent data for the 10 + 10 branches of [3,1,0] and [3,2,0] to the needed depth (W4's exact
  Rayleigh–Schrödinger extended; U1's no-log solver with the degenerate case fixed — see WP-C). Duality check: [3,2,0]'s series should be the
  dual images of [3,1,0]'s.
- Acceptance: all 20 CG tables of [3,1,0]² and [3,2,0]² from the Q side to ≥ 30 digits, φ-independent.

### WP-C. The mechanism and closed form of the ε-series of the Laurent data (the "perturbative Bethe equations")

- C1 (the degenerate solver): U1's order-by-order solution of the six no-log conditions + two global rules failed where the Jacobian has rank 5
  (top state and octets) — implement the correct degenerate-solvability step (Fredholm alternative at each order; the kernel direction is fixed
  by the next order's solvability), symbolic in φ. Deliverable: exact series for ALL adjoint branches from the Bethe side alone (no chain), to
  n = 8; then [3,1,0], [3,2,0].
- C2 (the algebraic forms): H₁(ε) for (4,1,1) and (3,3,0) is 6 + linear + (√(1 + x²/4) − 1), the eigenvalue of a 2×2 problem. Identify the two
  states that mix (the HW-limit branch and which neighbour?), derive the 2×2 matrix from the twisted Casimir structure (the exact rule H₁ + H₂ =
  Σ_{i≠j}ζ_iμ_j says the twisted Hamiltonians are Casimir-like), and generalise: which branches are algebraic of which degree? (the singlet and
  octets are not algebraic of degree 2 — are they of higher degree, or transcendental?) A conjecture for the structure of H_s(ζ) on a HW-limit
  branch for general Λ (e.g. eigenvalues of a finite matrix built from the ζ-twisted gl(3) Casimirs restricted to the μ-isotypic component).
- C3 (closed form via Casimirs, chain side, Wolfram or Python exact): at finite ζ the Laurent data of an eigenstate are eigenvalues of the
  Laurent-coefficient operators; on the μ-isotypic component of V⊗V these operators are expressible through the Casimirs of the total gl(3), the
  twist element Σζ_kE_kk and their products — write them in that algebra (exact), diagonalise symbolically on each isotypic component (dimension =
  multiplicity × dim μ, but the twist only mixes within weight spaces: small matrices), and read off H_s(ζ), a_s(ζ) as algebraic functions. This
  would be the closed form of the ε-series and the Gaudin analogue of "solving the Bethe equations in the twisted Casimir basis".
- Acceptance: a formula (or algorithm with a proof) producing the ε-series of every branch of [2,1,0], [3,1,0], [3,2,0] exactly, matching U1/W4.

### WP-D. The normalisation: ⟨Ψ₀|Ψ₀⟩ from the Q side (the Gaudin norm)

The CG formula gives the generating function up to one constant per branch (G-norm). Classical result (Reshetikhin–Varchenko / Mukhin–Varchenko):
the norm of a Gaudin Bethe vector is the Hessian of the Yang–Yang function at the Bethe roots. Here the "Bethe equations" are the no-log
conditions on the Laurent data (six per state). Hypothesis: ⟨Ψ|Ψ⟩_G/⟨x₀|Ψ⟩² ∝ det(Jacobian of the no-log conditions w.r.t. the Laurent data)
times explicit factors, or ∝ the equal-twist det N (P1) suitably normalised (D3 found Qnum_eq·M₁M₂/det N_B state-independent — that IS a norm-type
statement: ⟨Ω_ω⊗Ω_ω|Ψ⟩²-type quantities vs det N).
- D1: test numerically on the 64 adjoint states at finite twist: ⟨Ψ|Ψ⟩_G·(normalisation-free combination) vs (i) det N_B(ω) with the Lagrange
  columns, (ii) the Jacobian determinant of the six no-log conditions, (iii) the Wronskian-type determinants of the quantum-Weyl formulas; find the
  one that is state-independent after dividing.
- D2: its untwisting limit (the G-norm of Ψ₀ from the Q side) — then the CG coefficients are fully determined from Q-data, no chain input at all.
- Acceptance: a Q-side norm formula verified on all states of [2,1,0] and on the HW branches in the regime.

### WP-E. The write-up

- E1: the NOTES (dispatched now, documentation only): `Paul/Mathematica/Gaudin/Gaudin-sl3-CG-notes.{tex,md,pdf}`, hep-th conversational style,
  complete data, readability first — see the notes agent's brief.
- E2 (after WP-A/B at least partially): the paper. Decision for Paul: `paper/gaudinCG.tex` is a gl(4)-rectangular stub — replace its body with the
  sl(3) Gaudin CG paper (recommended), or start a new file `paper/gl3GaudinCG.tex`. Proposed structure: 1 Introduction (CG from SoV; why Gaudin);
  2 The two-site Gaudin model, cdet Hamiltonians, Laurent data, opers and Q-functions (gauge, exponents, Wronskians); 3 Vacua, frames, the
  two-vacua family and its CG-completeness (chain side, no cleaning); 4 The jet formula and its structure (jet basis = Gaudin SoV, K tensor, block
  theorem, one-site kernel, the master kernel); 5 The overlap as a function of the Laurent data; 6 The untwisting limit and the CG generating
  function; 7 Worked example: adj⊗adj; 8 Outlook (general Λ, norm, higher rank). Appendices: conventions, kernels, series.

### WP-F. Housekeeping and smaller items

- F1: Paul reruns the four Gaudin notebooks in VS Code (no native outputs stored by the agents); gitignore `*.pkl` (≈ 230 MB in `Paul/Python/
  gaudin_sl3/`) and commit the Gaudin tree (plans, specs, notebooks, Python, RESULTS, ref data); update CLAUDE.md's goal statement (Paul).
- F2: promote `Gaudin-L2-chain.wb` (and -2, -320, untwisted-exact) to `Paul/Mathematica/Gaudin/Clean/` once Paul has rerun them; docs agent.
- F3: the "worked example" script `cg_from_Q.py`: input = G (closed-form pieces where available) + the Laurent-data series + Gauss parametrisation;
  output = the adj⊗adj CG table — a one-file demonstration of the whole chain of results (good for the paper and for sanity).
- F4 (backlog): [4,2,0] (dim 27, multiplicity 3) — the last cleaning check; the (5,3,2) exception of the multiplicity-2 law and β(Λ, μ); the
  auxiliary sl(2) on the jet space; a remark on the Gaudin B hierarchy for the paper; the λ₂ = 0 reps as master objects (also for the paper).

## 2. Suggested agent allocation (NOT dispatched; all SONNET)

| WP | agent | type | kernel |
|---|---|---|---|
| A1 + A3 | M1 | Python exact (sympy/Fractions), 40–80 digits | none |
| A2 | M2 | Python exact (polynomial realisation of Sym^{λ₁}) | none |
| B1 + B2 | N1 | Python (sympy + 200-digit regime) | none |
| B3 + C1 + C3 | W5 | wolfbook-builder (exact perturbation theory, Casimir algebra on isotypic components) | K1 |
| C2 | S1 | Python exact (2×2 mixings, algebraic degrees) | none |
| D1 + D2 | Q1 | Python (norm tests on 64 states + regime) | none |
| E2 | after A/B | writer | none |

Dependencies: B2 needs G1's [3,1,0] laws (exist) and P4's [3,2,0] laws (exist); B3/C1 feed B2's exact series; A2 is independent; D is
independent; E2 waits for A (at least the status) and B. Run ≤ 4 Python agents at once; W5 alone on K1.

## 3. Decision points for Paul

1. Paper file: rewrite `paper/gaudinCG.tex` or a new `gl3GaudinCG.tex`? (E2)
2. Priority between WP-A (the closed-form kernel — the "general formula") and WP-C/D (derivation of the series, the norm — the "no chain input"
   completion)? Both are needed for a paper; A is the more striking result, C+D the more complete one.
3. Whether to spend a round on [4,2,0] (multiplicity 3) before writing, or accept the three reps as the evidence base.
4. CLAUDE.md: update the goal statement to gl(3) Gaudin CG (and mention that `old/Gaudin/` is the old gl(4) material).

## 4. Risks

- WP-A may have no closed form in Pochhammer symbols (the rank pattern 1,2,4,8 then capped hints at a "fusion" structure rather than a product
  formula); A2's derivation route is the fallback and also the proof route.
- WP-C3's Casimir diagonalisation may give algebraic functions of high degree for multiplicity-2 sectors (Galois conjugate copies already appear);
  that is still a closed form.
- WP-D's norm hypothesis is a guess; if none of (i)–(iii) is state-independent, the Gaudin norm needs the Hessian of the actual Yang–Yang function
  in terms of Bethe roots (Q-function zeros), which we have avoided so far — a different formalism, flag and stop.
