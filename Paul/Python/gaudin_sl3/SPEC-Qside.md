# Spec P1: the Q side of the sl(3) Gaudin chain in Python — opers, Q-functions, Wronskians, residue brackets, det N (2026-10-15)

Delegated to a general-purpose subagent (sonnet), Python only (the Wolfram kernel is taken by another agent). Read first:
`Paul/Mathematica/Gaudin/Gaudin-sl3-CG-plan.md` §0, §1 (binding conventions), §3, §7. Reference: `Paul/Mathematica/Gaudin/ref/
ref-Companion-B-C-diagonalize-clean-all.txt` (the XXX L = 1 Baxter/dual Baxter/QQ/Wronskian/orthogonality/det N cells — the
objects whose Gaudin limits you are building) and `Paul/Mathematica/new-conventions/Companion-B-C-diagonalize-clean.md`,
`L2-diagonal-frame-overlap.md` (the L = 2 det N: rows (s,(i,j)), columns, Ψ_Q, the x₀-free ratio form in Section 13).
Existing Python: `Paul/Python/gl3_reps.py` (GT generators E_ij of gl(3) in sympy; `test_gl3_reps.py` shows usage) — reuse it.

**Files you own** (create under `Paul/Python/gaudin_sl3/`): `gaudin_chain.py` (model + eigenstates), `qside.py` (Tasks Q1–Q3),
`RESULTS-Qside.md` (all results, numbers with precision, conventions fixed). Do not touch `xxx_*.py` or `RESULTS-xxx-limit.md`
(another agent's). Run with `python -X utf8` (Python 3.13; sympy 1.14, mpmath 1.3, numpy available; no `uv`). For exact
64×64 operator algebra use `fractions.Fraction` in numpy object arrays or sympy `DomainMatrix` over QQ (dense sympy `Matrix`
is too slow); for eigenproblems use mpmath at 60–100 digits per weight sector (sectors are ≤ 10-dim in the HW weights, ≤ 10 in
general for the adjoint). Time one 64×64 exact product before loops.

## 0. The model (identical to the W1 spec; a separate implementation is the point)

Λ = (2,1,0) on both sites (then (3,1,0) if everything passes), θ = (2/7, −3/5), E⁽¹⁾ = E⊗1, E⁽²⁾ = 1⊗E (Kronecker order, same as
Mathematica `KroneckerProduct`), covectors act from the right (row vectors), HW ket e₁ = first GT pattern as ordered by
`gl3_reps` — check that e₁ is annihilated by E_12, E_13, E_23 (if the ordering differs from the notebook's `ps` order, say so and
use the actual HW index). Twist g = diag(ζ), finite-twist point ζ = (3/2, 5/4, 2/3), bra twists ω = (5/2, 7/3, 3/5), ω′ = (7/4, −2/5, 9/7).
L_g(u)_ij = g_ij + Σ_s E⁽ˢ⁾_ij/(u−θ_s). cdet(∂_u − L_g) = Σ_σ sgn σ M_{σ(1)1}M_{σ(2)2}M_{σ(3)3}, M = ∂ − L_g, ∂ moved right (∂f = f∂ + f′)
=: ∂³ − S₁∂² + S₂∂ − S₃. Do the normal ordering symbolically once (27 noncommuting symbols l_ij, l′_ij, l″_ij; e.g. represent
each factor as a pair (coefficient of ∂, operator part) and multiply with the Leibniz rule), print S₂, S₃ in symbols into the
results file, then substitute matrices at rational u. Checks: [S₂(u₀), S₃(u₁)] = 0 exactly (u₀ = 3/7, u₁ = −2/9), S₁ scalar,
weight conservation. If commutativity fails try the transposed placement E_ji and the row-ordered product; record which works.
Eigenstates: per weight sector, eigenvectors of S₂(u₀) + π S₃(u₁) (mpmath); eigenvalue FUNCTIONS S₂^B(u), S₃^B(u) as exact
partial fractions (evaluate at enough rational u, fit the pole ansatz: S₂ poles ≤ 2, S₃ poles ≤ 3 at each θ_s, plus constants;
or get them as rational functions from the eigenvector by Rayleigh quotients at symbolic u — the partial-fraction fit is
simpler). Bra vacuum: Ω_{ζω} = e₁ᵀρ(n), n = exp(r₁E₁₂)exp(r₂E₁₃)exp(r₃E₂₃) with the r_k of plan §1 (ρ of a nilpotent via the
finite exponential series); g′ = P_aux(ζ)g_c(ω)P_aux(ζ)⁻¹ (plan §1); verify Ω⊗Ω is a left eigencovector of S₂^{(g′)}(u₀) and
S₃^{(g′)}(u₀) (try the variants listed in the plan if not) and record its eigenvalue functions S_k^{vac,ω}(u).

## 1. Tasks

**Q1 — the opers and the Q-functions.** For each of the 64 eigenstates B (weight μ): test the two candidate ODEs
  (D) (∂³ − S₁^B(u)∂² + S₂^B(u)∂ − S₃^B(u)) Q = 0 and (A) its formal adjoint (−∂)³Q − (−∂)²(S₁^B Q) + (−∂)(S₂^B Q) − S₃^B Q = 0
  (i.e. −Q‴ − (S₁Q)″ − (S₂Q)′ − S₃Q = 0), for solutions of the form e^{κu}q(u) with q polynomial: for κ ∈ {±ζ₁, ±ζ₂, ±ζ₃} find
  the smallest polynomial degree admitting a solution (linear algebra on the coefficients; exact rationals; allow a tolerance
  only if eigenvalues are numeric — then use 60-digit mpmath and report residuals). Expectation (plan §1): one of (D),(A) has
  solutions e^{+ζ_j u}q^j with deg q^j = μ_j (dual Q^j), the other e^{−ζ_i u}q_i with deg q_i = 4 − μ_i (direct Q_i; for Λ = (2,1,0),
  L = 2: 2λ₁ − μ_i). Fix the assignment (which ODE is "dual", which "direct") and the sign of the exponent — record it as THE
  convention in RESULTS. Gauge fixing: as in `old-conventions/old/Gaudin/Clean/GaudinSolver.wb` (strictly increasing degrees
  → kill the coefficient of u^{deg q_j} in q_k for j < k — check whether degrees are strictly ordered here; for the dual Q's
  with degrees μ₁ ≥ μ₂ ≥ μ₃ they are not — then no gauge mixing is possible between different exponentials anyway since
  e^{ζ_j u} differ; state the solution-space dimension per κ). Verify also for the vacuum eigenvalue functions of twist ω:
  direct Q_i^{vac,ω} = e^{−ω_i u}q_i with degrees (0, 2, 4) expected.

**Q2 — Wronskian (QQ) relation and quantum-Weyl formulas.** (a) For each state, W_ab := Q_aQ_b′ − Q_bQ_a′ (direct Q's) vs Q^c
  (dual, c ∉ {a,b}): find the ratio W_ab/Q^c — expected e^{−(ζ_a+ζ_b)u}·(…)·F(u) hmm: check whether W_ab/(e^{−(ζ_a+ζ_b)u}Q^c·
  e^{−ζ_c u}) is a state-independent polynomial F(u) times a constant c_ab(ζ); find F (expected ∏_s(u−θ_s)^{λ₁−λ₂} = (u−θ₁)(u−θ₂)
  for the adjoint) and c_ab (expected ∝ (ζ_b − ζ_a); compare with the XXX (1/z_b − 1/z_a)χ₃^{−u/h}F — note χ₃^{−u/h} → e^{−Σζ u}).
  Same with direct ↔ dual swapped (dual Wronskians vs direct Q's). (b) Quantum-Weyl: with the 3×3 Wronskian-type determinants
  D[k] := det[Q^j, Q^{j′}, Q^{j″}] (dual Q's and derivatives) and the shifted variants, express S₂^B(u) and S₃^B(u): the XXX
  rows 10–11 (t_a = ν-factors × D_{…}/D_{0,−1,−2}) limit to ratios of determinants of derivatives; find the Gaudin form by
  fitting a small ansatz (S₃ = c·det[Q, Q′, Q″]′/det[Q,Q′,Q″]-type combinations — try S₁ = (ln W)′ with W = det[Q^j,Q^j′,Q^j″],
  and for S₂, S₃ the standard formula for the monic 3rd-order operator with solution basis Q^j: ∂³ − S₁∂² + S₂∂ − S₃ =
  W⁻¹·det[[Q,Q′,Q″,Q‴]-type bordered determinant] — i.e. S_k are the coefficients of the unique monic operator annihilating the
  three Q^j, computable directly from the Q's: compute them and compare with the eigenvalue functions). (c) Report degrees,
  residuals, and the explicit F, c_ab.

**Q3 — residue brackets, functional orthogonality, det N (x₀-free).** Brackets per site: ⟨f⟩_s := Res_{u=θ_s} f(u)/(u−θ_s)^{λ₁+1}
  = f^{(λ₁)}(θ_s)/λ₁! (λ₃ = 0). (a) Orthogonality search: for pairs of states A ≠ B (same ζ), build for each (i,j) ∈ {1,2,3}² the
  candidate functionals F₁ = ⟨P^jQ_i⟩, F₂ = ⟨P^jQ_i′⟩, F₃ = ⟨P^{j′}Q_i⟩, F₄ = ⟨P^jQ_i″⟩, F₅ = ⟨P^{j″}Q_i⟩, F₆ = ⟨P^{j′}Q_i′⟩, … with
  P^j = dual Q of the ket B and Q_i = direct Q of the bra A (the XXX orientation: ket dual, bra direct), summed over s = 1,2 with
  per-site weights allowed (c_s). Find the linear relations Σ_k c_k F_k = 0 that hold for ALL pairs A ≠ B with coefficients c_k
  depending only on (i, j) and the twist — rank/nullspace over the stacked pairs (exact if all Q's are exact; else 60-digit SVD).
  The XXX relation −δI₁⟨P^jQ_i⟩ + a⟨P^jQ_i^{[2]}⟩ + b⟨P^{j[−2]}Q_i⟩ = 0 limits to a relation among F₁, F₂, F₃ (and possibly F₄–F₆
  from the second order of the expansion). Report what you find; derive the coefficients in closed form if they are simple
  (a, b in XXX were functions of i, j, z). (b) det N: rows (s, (i,j)) for s = 1,2 and (i,j) ∈ {(1,3),(1,2),(2,3)}, columns = the
  functionals that make det N vanish for A ≠ B and not for A = B (try the per-site residue versions of the XXX columns
  {−⟨PQ⟩, ⟨PQ′⟩, ⟨P′Q⟩} and the moment versions {⟨PQ⟩ with pole orders λ₁+1, λ₁+2, λ₁+3}); report the rank of the 18×6 stacked
  matrix for A ≠ B (expect 5) and which column choice works. (c) The overlap identity, x₀-free: for the vacuum direct Q's of
  twists ω and ω′ (Q1), D_B(ω) := det N[Q^j_B, Q_i^{vac,ω}], and the chain overlaps R_B(ω) := ⟨Ω_ω⊗Ω_ω|Ψ_B⟩ (any normalisation of
  Ψ_B — it cancels): test R_B(ω)/R_B(ω′) = [D_B(ω)/D_B(ω′)]·[D_{B₀}(ω′)/D_{B₀}(ω)] with B₀ = e₁⊗e₁ (the top state) on all 64
  states; report max relative deviation and its precision. If it fails for the first column choice, try the others of (b); if
  the per-site brackets fail, try the XXX "all-poles" form in the Gaudin guise Res at both θ_s of P·Q·μ(u) with μ = ∏_s(u−θ_s)^{−λ₁−1}
  and with a second insertion u (the x^{β−1} columns) — i.e. a 6×6 with rows (i,j) × {insertion 1, insertion u}... the structure
  is 6 columns × 6 rows; report what closes. (d) With x₀^G if W1's notebook has defined it by then (not required): the
  Ψ_Q-normalised form; otherwise skip.

**Q4 (only if Q1–Q3 are done).** Regime ζ = εφ, ω = εω̂, ω′ = εω̂′ (plan §1 directions, one grid point, ε = 10⁻⁶): the HW-limit
  branches (smallest raising residual per HW sector), their Q's, D_B(ω), D_B(ω′), and the ratio test of Q3(c) on the branches;
  ε-exponents of the det N entries (the cancellation depth). Record for the later two-vacua Q-side stage.

## 2. Reporting

`RESULTS-Qside.md`: §0 conventions fixed (which ODE is direct/dual, exponent signs, the explicit symbolic S₂, S₃, any variant
used), §1–§4 per task with tables (degrees, residuals with precision, the relations found, the det N column choice, the
64-state ratio test), §5 surprises and open points. Keep scripts runnable top-to-bottom (`python -X utf8 qside.py` prints the
tables); no notebooks, no git, no edits outside `Paul/Python/gaudin_sl3/`. Numbers, not adjectives. If a task's expectation is
refuted and two attempts do not fix it, stop at that point, keep what works, report what you saw.
