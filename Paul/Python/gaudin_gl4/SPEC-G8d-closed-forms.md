# Spec G8d: closed forms of the nesting scalars a(κ), b(κ), the V₄ᵇ within-path parameters, their origin in the word structure, and the gl(4) recursion assembled from closed forms alone (2026-10-24)

Delegated to a general-purpose subagent (SONNET), **Python only** (no Wolfram, no wolfbook tools, no git, no notebooks). Work package G8d of
`Paul/Mathematica/Gaudin/gl4-jacobi-plan.md` — read its §4 (the verdicts) first, then `RESULTS-G8b-nested-gl3.md` (§0, §2 incl. §2.3 "What a, b are", §4, §8 usage),
`RESULTS-G8a-lattice-recursion.md` (§0, §1 the recursion formula, §9 usage), `RESULTS-G4-structure-gl4.md` §3 (closed forms: the universal product formula of the 12 hopping
classes, the cubic Q of the six Δm₁₂ = 0 classes of V₃, the within-path classes = linear factors × Q₆(κ,λ), the diagonal), and `Paul/Python/gaudin_sl3/RESULTS-jacobi-cg.md` §3
(the gl(3) chain closed forms for λ = λ′: u, d, α with the pole term K/ξ). Environment: Python 3.13, numpy, scipy, sympy, mpmath; run from `Paul/Python/gaudin_gl4/` as
`python -X utf8 ...`; ≤ 2 heavy processes; `perf_counter`. **Never edit existing files; import G1–G5, G8a, G8b read-only.** Files you own: `closed_forms_gl4.py`,
`test_closed_forms_gl4.py`, `RESULTS-G8d-closed-forms.md`, `closed_forms_gl4_output.txt`, `closed_forms_gl4_tables.txt`. Time-box ≈ 3 hours; deliver what is finished and say
precisely what is not. Tooling trap: Python source written through Bash heredocs loses `\n` escapes — write code with the file tools and check with `ast.parse`.

## 0. What is established (exact, G8a/G8b) and what you must make explicit

W = V_λ[ν] = ⊔_κ P_κ (ν = μ − λ′; κ = row 3; P_κ = the path in m₁₂). Notation: y_ik = m_ik − i (l-coordinates), Y_s = λ_s − s, q = λ′, ξ = x(x+2) with x = m₁₂ − m₂₂.
* **G8b B1:** the κ-block of V₃ is J^{(3)}_κ = a(κ;λ)·V^{(3)}_κ + b₃·1, of V₄ᵃ is a₄ₐ(κ;λ)·V^{(3)}_κ + b₄ₐ·1, where V^{(3)}_κ = `slice_ops_gl4.reduced_op('V', lamp[:3])` on
  `glN_reps.GLNRep.of(kappa)` at weight ν₁₂₃ (the gl(3) reduced vertex operator = `slice_engine.vred_apply` + η₁/2). a, a₄ₐ depend only on (κ, λ) (not on λ′, ν); a(λ₁₂₃) = 1,
  a₄ₐ(λ₁₂₃) = 0 for λ₄ = 0. b₃ is quadratic in ν and affine in λ′ (example for λ = (6,4,2,0), κ = (5,3,0) in RESULTS-G8b §2.3); b₄ₐ = 80 for all κ with κ₃ = 0 of that λ,
  (ν, λ′)-dependent otherwise. The table of (a, a₄ₐ) for all 21 κ of λ = (6,4,2,0) is in §2.3 (e.g. κ = (5,3,1): 32/35, 128/35; (5,2,1): 464/315, 2792/315 — not products of
  linear forms). G8b's API: `nested_gl3_gl4.within_data(lam, lamp, mu)`, `fit_ansatz(J, V, D)`, `within_V4b_fit(d)`.
* **G8b V₄ᵇ:** J^{(4b)}_κ = a′V^{(3)} + e[ξ, V^{(3)}] + b′ + cξ + gξ², g = −e/4 (305/305, 369/369 held-out); equivalently the off-diagonal entry (k′ ← k) is (p + r·y_a)·V^{(3)}[k′,k]
  with y_a the row-2 entry raised by the move — the same shape as G4's hopping-class factor −(y_a − y_i + 2y_k − q₃ + 5).
* **G4 §3.3:** the 12 hopping classes (Δm₁₂ = ±1, row-3 root move) have the universal product formula `structure_gl4.universal_F(op, pt, d)` (pt = `pt_of_pattern(pattern, lam, lamp)`,
  d ∈ `f_classes()`); the within-path classes (±1,0,0) are "linear factors × Q₆(κ, λ)" (Q₆ of degree 6 in (m₁₃, m₂₃, m₃₃, λ) only); the six Δm₁₂ = 0 root classes are "linear
  factors × cubic Q" with the six cubics of V₃ explicit in `structure_gl4_tables.txt` (V₄ᵃ has the same Q on them; V₄ᵇ a quartic, explicit for one class); the diagonal has degree 11/12.
  G4's machinery: `closed_form_one(name, d)`, `reconstruct_q(rec)`, `eval_Qexp`, `verify_closed`.
* **G8a:** `lattice_recursion_gl4.slice_from_triple(lam, lamp, mu, trip, grading='best', space='W1')`, `sector`, `sym_recursion`; the recursion F_d ψ_{d+1} = [v_k ψ_d − Σ_known M_k ψ]_k
  with F_d = [M_k[L_d, L_{d+1}]]_k; grading rule per active-string class A: +m₁₂ for A ∈ {∅, {1}, {2}, {2,3}, {1,2,3}}, −m₁₂ for {1,2}, {1,3}, +m₃₃ or −m₁₃ for {3}.

Goal: make every coefficient of the gl(4) recursion a closed-form function of the GT labels, λ, λ′, μ — so that, as for gl(3) (`jacobi_cg`: "the core imports no
representation matrix"), the slices of all copies of a sector are computed from the labels alone — and explain where a(κ) and the V₄ᵃ/V₄ᵇ dressings come from.

## 1. Questions

**D1. Closed form of a(κ; λ) and a₄ₐ(κ; λ).** Compute a(κ) exactly for every κ ≺ λ (interlacing) of λ ∈ {[6,4,2,0], [9,6,3,0], [5,3,2,1], [7,4,1,0], [8,5,3,1], [6,6,3,0]} with
G8b's `within_data` (any ν with |P_κ| ≥ 2 and any λ′, since a is independent of both — verify once more on 3 random λ′ per κ). Then identify a as a rational function of
(κ, λ): (i) the Wigner–Eckart-shaped ansatz suggested by the index-4 words, a(κ) = 1 + Σ_{i=1}^{3} R_i(κ, λ)·(affine in y_i, with unknown coefficients), where
R_i = Π_{s=1}^{4}(y_i − Y_s) / Π_{i′≠i}(y_i − y_{i′})(y_i − y_{i′} + 1) is the GT factor of "raise row-3 entry i by E₃₄ then lower it back" (also try the lowering-first version
Π_s(y_i − Y_s − 1)/Π(y_i − y_{i′})(y_i − y_{i′} − 1) and the symmetric combination); fit the few unknowns on [6,4,2,0] by exact linear algebra and VERIFY on the other five λ
(≥ 100 κ values, 0 mismatches required; report the number of free parameters and the held-out count); (ii) if (i) fails, exact rational interpolation: a·D(κ) with D = the
product of the denominators Π_{i<i′}(y_i − y_{i′})(y_i − y_{i′} + 1)…, is a polynomial in (κ, λ) — find its degree, interpolate it exactly in the 7 variables (κ₁,κ₂,κ₃, λ₁..λ₄)
with sympy on an exact grid (G4's `LineProbe`/line-reconstruction method is the template), factor it over ℚ, and present it in l-coordinates. Same for a₄ₐ (expected: a₄ₐ = Σ_i
(y_k-type factor)·…, cf. M^{V₄ᵃ} = (y_k + 3)M^{V₃} on the hopping classes — test a₄ₐ(κ) = Σ_i R_i·(y_i + c) type ansätze and the relation between a₄ₐ and a: is
a₄ₐ − c₀·a a simple function?). Give the final formulas in l-coordinates and in the labels, and the verification counts.

**D2. Closed forms of b₃, b₄ₐ and of the V₄ᵇ parameters (a′, e, b′, c).** With exact samples over (ν, λ′) for each κ (G8b did 72 samples for one κ; do all κ of [6,4,2,0] and
[5,3,2,1]): fit b₃(κ; ν, λ′) = quadratic in ν + bilinear (ν, λ′) + affine λ′ + constant with coefficients depending on (κ, λ), then identify the κ-dependence of each
coefficient (it should be expressible through gl(3)/gl(4) Casimir-type combinations: |κ|, e₂(κ) or Σκ_i(κ_i + 4 − 2i), |λ|, …, and the depth d = ν₄ − λ₄). Verify on the
held-out λ. Same for b₄ₐ, a′, e, c, b′ (V₄ᵇ; recall g = −e/4). State which of these are "Casimir scalars" (functions of (|κ|, C₂(κ), |λ|, C₂(λ), ν, λ′)) and which need more.
Cross-check against G5's scalars: on depth-0 sectors with λ₄ > 0 G8b found b₄ₐ − b₃ = 18 for [5,3,2,1]² — explain it from the formulas.

**D3. Where it comes from (a derivation sketch, not a proof).** Split the 39 universal reduced words of V₃ (`slice_ops_gl4.reduce_terms(so.op_terms('V3'), lamp)` or
`universal_reduction('V3')`) into three groups: (α) words with all indices in {1,2,3} ∪ Cartan E₄₄ (they act within κ-blocks); (β) words containing exactly one E_{i4} and one
E_{4j} (they pass through the neighbour blocks κ + e_i and come back; their κ-block is the Wigner–Eckart piece); (γ) the rest. Compute the κ-block of each group's operator
separately (exact, with `weight_space_mats` on a group-restricted ReducedOp — construct a ReducedOp from a sub-dict of `.red`) on all paths of [4,2,1,0]² and [6,4,2,0]²:
confirm group (α) = V^{(3)} + scalar, group (β) = (a − 1)V^{(3)} + scalar, group (γ) = scalar or 0 (report whatever it is), and that the SAME split for V₄ᵃ explains
a₄ₐ (which group carries the (y_k + 3)-type dressing — the words with a Cartan E_kk inserted?). Write down the mechanism in one paragraph with the formula for (β)'s κ-block
if you can read it off (a sum over i of R_i × the gl(3) operator), and say whether it proves D1's formula modulo the gl(3)-covariance argument.

**D4. The recursion from closed forms alone.** Implement `closed_form_mats(lam, lamp, mu)` returning the three exact w×w matrices on W built ONLY from: the gl(3) chain on each
path (derive or take the general-λ′ gl(3) closed forms: G4's structure gives the twist factor (y₁₁ − y_a + q₁ − q₂); the λ = λ′ forms of `jacobi_cg` (u, d, α, K/ξ) must be
generalised — obtain the general-λ′ N = 3 chain entries symbolically with G4's method on `so.reduced_op('V', lamp3)` for N = 3, or by the same line-reconstruction, and
verify on ≥ 500 gl(3) weight spaces with generic λ′) times a(κ) plus b (D1, D2); `universal_F` on the 12 hopping classes; the explicit cubic Q's on the six Δm₁₂ = 0 classes for
V₃ and V₄ᵃ (from `structure_gl4_tables.txt` / `reconstruct_q`) and, for V₄ᵇ, the quartic Q's (reconstruct the remaining five with G4's `closed_form_one` + `reconstruct_q`, or
state that V₄ᵇ is left out and give the recursion with (V₃, V₄ᵃ) only where those two suffice — G8a's determinacy check with two operators: run it and report the fraction).
Acceptance: `closed_form_mats` == G3's `weight_space_mats` EXACTLY on every weight space of P_small ∪ P_big ([2,1,0,0]², [2,1,1,0]², [3,2,1,0]², [4,2,0,0]², [4,2,1,0]²,
[3,2,1,0]⊗[4,2,1,0] both orders, [4,2,2,0]², [4,3,1,0]², [6,4,2,0]², [5,3,1,0]², [5,3,2,1]²) and on 20 random weight spaces of [9,6,3,0]² with generic λ′. Then
`slices_closed(lam, lamp, mu)` = all copies of a sector from the closed-form matrices + G8a's recursion (eigenvalues from the m×m PRV matrices or from the closed-form
W₁ matrices' joint spectrum; say which), validated against `so.copies` on [6,4,2,0]² (all 297 copies, 50 digits) and timed per sector and for the whole product against G3's
`copies` (which builds the operators by applying the 88/152/324-term words to GT vectors) and against G8a's `slice_from_triple` (which still uses G3's matrices); report the
breakdown: closed-form matrix assembly / eigenvalues / recursion. Also [9,6,3,0]² (1592 copies) if it runs in < 10 min; else the largest subset that does.

## 2. Deliverables and acceptance

* `closed_forms_gl4.py` with a CLI `python -X utf8 closed_forms_gl4.py a | b | words | gl3chain | mats | slices | all` (and `run <stage,...>`), tables to
  `closed_forms_gl4_tables.txt` (all formulas printed in both l-coordinates and labels; the (κ → a, a₄ₐ) tables for all six λ), printed output of the final run to
  `closed_forms_gl4_output.txt`; public `a_of(kappa, lam)`, `a4a_of(kappa, lam)`, `b_of(name, kappa, lam, nu, lamp)`, `gl3_chain(kappa, lamp3, nu123)` (u, d, α for general λ′),
  `closed_form_mats`, `slices_closed`.
* `test_closed_forms_gl4.py`: ≥ 6 exact regression tests (< 60 s): a(κ) formula on 5 held-out κ; a₄ₐ on 3; b₃ on 2 samples; the general-λ′ gl(3) chain on 2 weight spaces;
  `closed_form_mats` == G3 on 2 weight spaces of [6,4,2,0]²; `slices_closed` == `copies` on 1 sector.
* `RESULTS-G8d-closed-forms.md`: §0 headline (the formulas for a, a₄ₐ, b₃, b₄ₐ, the V₄ᵇ parameters with their held-out counts; the word-group mechanism in two sentences; the
  closed-form recursion's validation count and timings), §1 conventions, §2–§5 per question, §6 surprises and corrections (to G8a/G8b/G4 or this brief), §7 not done / limits,
  §8 usage. Every claim exact or "50 digits". If a formula resists identification, give the exact table and the degree/denominator information you did establish.
