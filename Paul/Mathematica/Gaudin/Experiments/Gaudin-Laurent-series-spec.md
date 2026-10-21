# Spec W5: exact ε-series of the Laurent data for [3,1,0] and [3,2,0], and their algebraic structure (the cyclic-subspace / Casimir closed form) (round 5, WP-B3 + C3 + C2; 2026-10-16)

Delegated to a wolfbook-builder subagent (SONNET); planned by the main session. **Notebook:** `Paul/Mathematica/Gaudin/Experiments/
Gaudin-Laurent-series.wb` — created by the main session (two cells: title, `Quiet[ClearAll["Global`*"]]`), opened in VS Code and bound to the
default kernel K1 — so `wolfbook_runCell`/`wolfbook_insertCells` with `evaluate:true` should work and store native outputs. If a run answers
"NOT executed"/"controller not associated", fall back to the known workaround: evaluate the cell's exact saved source in K1 with
`wolfbook_evaluateExpression` (large `timeoutSeconds`) and store the printed result as a markdown cell under the code cell. Save often
(`wolfbook_getNotebookContext` action "save"). K1 is shared with nobody else this round; the Python agents (M1, M2, N1, Q1) must not be touched —
you must NOT run Python computations of your own beyond reading/dumping notebook JSON.
Read first: `../Gaudin-sl3-round5-plan.md` §0 (P3) and §1 WP-B3, WP-C (C1 belongs to Python agent S1 later; your C-items are C3 and C2),
`../Gaudin-sl3-round4-plan.md` §0 (notation — binding) and its W4 result block, the notes `../Gaudin-sl3-CG-notes.md` §1 (conventions), §10 (the
ε-series), §11 (untwisted model), Appendix D, and `Gaudin-untwisted-exact-spec.md` (W4's spec). **Templates:** `Gaudin-untwisted-exact.wb` (W4:
setup — generators/metric/Lax/cdet/`S[k]`/Laurent operators `Hres`, `S3c`, exact HW kets `hwKets`/`PHWc`, universal copies, E4 = exact
Rayleigh–Schrödinger in ε on the weight-μ space with the isolating operator H₁ + a₁/7, E5 = x(ζ), k_B) and `Gaudin-L2-chain-320.wb` (W3: [3,2,0]
specifics). Dump their cell sources with `python -X utf8 -c "import json,sys; ..."` to a UTF-8 text file in your scratchpad and COPY the machinery you
need into the new notebook (self-contained; same symbol names where you reuse them; first code cell `Quiet[ClearAll["Global`*"]]`); do not
evaluate the old notebooks.

## 0. Conventions

Two-site gl(3) Gaudin, Λ on both sites, θ = (θ₁, θ₂) symbolic where feasible (scaling: H_s, a_s at order εⁿ ∝ θ₁₂^{n−1}, b_s ∝ θ₁₂^{n−2} — verify once
per rep, then use θ₁₂ = 1 or numeric θ = (2/7, −3/5)); twist ζ = εφ with φ SYMBOLIC where feasible and otherwise at the Python agents' numeric points
φ = (1/3, −2/7, 5/11) and φ′ = (−2/5, 3/8, 7/11) (both are needed by N1). Laurent data: S₂ = e₂(ζ) + Σ_s[H_s/(u−θ_s) + c₂/(u−θ_s)²], S₃ = e₃(ζ) +
Σ_s[b_s/(u−θ_s) + a_s/(u−θ_s)²], c₂ = λ₂(λ₁+1), H_s, a_s, b_s = eigenvalues of W4's operators `Hres[s]`, `S3c[s,2]`, `S3c[s,1]` (check W4's naming).
HW-limit branch of sector μ = the eigenstate of the twisted family tending to the HW vector Ψ₀ (G-norm 1, first nonzero GT component positive;
multiplicity-2 copies = eigenvectors of the θ-independent A₁ = θ₁₂a₁ on the HW space, in Q(√β): β = 105, 192 for [3,1,0]'s (5,2,1), (4,3,1) and for
[3,2,0]'s (5,4,1), (5,3,2) — W4's E2; copy 1 = smaller α₁ as in W4's files). All perturbation theory is on the weight-μ subspace of V⊗V (the twisted
Hamiltonians commute with the total Cartan) — exact rational linear algebra, never numerics.

## 1. Tasks (one code cell + one markdown result cell per task and rep; a `checks` association of named booleans; save after each)

**B3 — exact ε-series for every HW-limit branch of [3,1,0] (10 branches) and [3,2,0] (10).** Port W4's E4 (degenerate Rayleigh–Schrödinger with the
isolating operator — check that H₁ + a₁/7 still separates every branch at ε → 0 in these reps; if not, use H₁ + a₁/7 + b₁θ₁₂²/11 or diagonalise the
commuting pair exactly on the degenerate eigenspace) to both reps: (H_s, a_s, b_s) for s = 1, 2 to order ε⁵ at φ and at φ′ (exact rationals, and
in Q(√β) for the multiplicity-2 copies), symbolic φ to order 2 (order 3 if it runs in < 20 min per rep). Checks per branch: H₁ + H₂ =
Σ_{i≠j}ζ_iμ_j (exact, all orders — it is linear in ζ), b₁ + b₂ = ε²(φ₁φ₂μ₃ + φ₁φ₃μ₂ + φ₂φ₃μ₁)? (adjoint rule — test), the no-log relations
(b_s = a_s(H_s − a_s)/4 for [3,1,0]; P4's inhomogeneous quadratic for [3,2,0] — find it in `Gaudin-L2-chain-320.wb`/the notes §2.3 and verify order by
order), ε⁰ values = W4's untwisted data, the θ-scaling. Duality: [3,2,0] ≅ [3,1,0]*⊗det³ — state the map on sectors and on ζ (the dual of the twist
is ζ → −ζ composed with…; derive it from the chain) and check that the [3,2,0] series are the images of the [3,1,0] ones. **Write the data files
`Paul/Mathematica/Gaudin/ref/laurent-series-310.txt` and `-320.txt` as soon as each rep's numeric-φ series exist** (N1 and S1 read them): same
TAB-separated InputForm format and header style as W4's `ref/untwisted-exact-*.txt` E4_SERIES/E4_SYMBOLIC records (one record per branch, per φ-point,
per order; sqrt's as (vA, vB) pairs with the radicand; document the branch labelling and copy convention in the header). Then also (B3b, cheap
once the series exist): the eigenstate series Ψ_n to order 4, x(ζ) for both reps (W4's XVEC, degree-λ₂ polynomial covector — reuse), k_B := the
first n with ⟨x(φ)⊗x(φ)|Ψ_n⟩ ≠ 0 and the exact coefficient X_B, written as E5_KB-type records; report the k_B pattern (adjoint: 0,1,1,2,2,4) and
test the rule "k_B = the number of lowering steps from the top weight to μ in some grading" against both reps.

**C3 — the closed form: the cyclic subspace of the HW vector.** On a branch the Laurent data are eigenvalues of commuting operators H₁(ε), a₁(ε),
b₁(ε) (and site 2) on the weight-μ space. The eigenvalue on the branch is an algebraic function of ε whose minimal polynomial is the
characteristic polynomial of H₁(ε) restricted to the cyclic subspace C_B = span{Ψ₀, H₁Ψ₀, H₁²Ψ₀, …} over Q(ε, φ) (use the isolating operator or
the commutative algebra generated by H₁, a₁ — the cyclic subspace of the algebra is what matters; compute both and compare). For every branch of
the ADJOINT first (six branches; φ symbolic if feasible, else numeric φ and φ′ with ε symbolic): d_B = dim C_B; the characteristic polynomial
χ_B(λ; ε, φ) of H₁ on C_B (exact; factor it over Q(φ) or Q(√5)(φ)); the series of its branch root vs U1/W4's series (must agree to the available
order); same for a₁ and b₁ (or give the whole Laurent data as polynomials in the H₁-root over Q(ε,φ) — the "primitive element" form). Expectations:
top: d = 1; (4,1,1), (3,3,0): d = 2 with the √(1 + ε²(φ_i − φ_j)²/4) forms — identify the partner state in C_B explicitly (which weight-μ state of
which irrep; it should be the ε-independent second basis vector of C_B after Gram–Schmidt) and derive the 2×2 matrix from the structure of H₁ =
−[C_ζ⁽¹⁾ + Ω₁₂/θ₁₂] + |Λ|σ₁ + |Λ|²/θ₁₂ (W4/round-1 form; Ω₁₂ = Σ_ik E⁽¹⁾_ikE⁽²⁾_ki = ½(C₂^{tot} − C₂⁽¹⁾ − C₂⁽²⁾) — a Casimir, constant on irreps of the
total gl(3) — and C_ζ⁽¹⁾ the twist term, linear in ζ, mapping an irrep component to neighbouring irreps); octets and singlet: d = ? (the plan's
open question: algebraic of which degree, or is d the full weight-space dimension 10?). Then the general statement: d_B and χ_B for every branch of
[3,1,0] and [3,2,0] at numeric φ (ε symbolic), and a conjecture for d_B in terms of (Λ, μ) (e.g. the number of irreps of V⊗V reachable from μ by the
twist term within the weight-μ space, i.e. the length of the chain μ → μ ± (e_i − e_j) → … inside the HW decomposition). This is the "closed form
of the ε-series": the branch eigenvalue is a root of an explicit polynomial.

**C2 (if time remains) — the matrix form.** For the d = 2 branches of all three reps, write the 2×2 matrix of H₁ on C_B in closed form in (ε, φ, θ₁₂,
Λ, μ) and check it reproduces the sqrt forms; for d > 2, the d×d matrix in a canonical basis (the Krylov/Lanczos basis with the G-metric) — the
entries should be polynomials in (ε φ) of degree ≤ 2 with coefficients rational in the Casimirs: test and report.

## 2. Reporting

Final markdown cell in the notebook AND the data files. In your report: the series (headline orders), the duality verdict, the k_B pattern, d_B
per branch for all three reps, the characteristic polynomials (adjoint fully), surprises, not done. No native outputs are expected to survive
if the fallback was used — say which mode you used. Never edit `paper/*.tex`, never commit.
