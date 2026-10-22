# Spec G2: the two-site gl(4) Gaudin missing-label operators — cdet(∂ − L(u)) for gl(N), its Laurent coefficients, commutativity, the three mixed invariants, simple joint spectrum (2026-10-22)

Delegated to a general-purpose subagent (SONNET), **Python only** (no Wolfram kernel, no wolfbook tools, no git, no notebooks). This is work package
WP-G2 of `Paul/Mathematica/Gaudin/gl4-slice-engine-plan.md` (read its §0–§2 and §4 "Risks" first). Environment: Python 3.13.13 (`python -X utf8`),
numpy 2.5.2, scipy 1.18.1, sympy 1.14, mpmath 1.3 (no gmpy2, no `uv`); i7-12700KF, 32 GB, Windows 11. `time.perf_counter()`. ≤ 2 heavy processes.
**Never edit existing files; import them read-only.**

**Directory:** `Paul/Python/gaudin_gl4/` (new; this file lives there). Run from there. Files you own: `gaudin_gl4_ops.py`, `test_gaudin_gl4_ops.py`,
`RESULTS-G2-gaudin-gl4-ops.md`, `gaudin_gl4_ops_output.txt`, and `../../Mathematica/Gaudin/ref/gl4-gaudin-ops.txt` (the term lists, §3 G2.5).

Another Sonnet agent (WP-G1, `SPEC-G1-slice-engine-gl4.md`, same directory) is **concurrently** writing `glN_reps.py` (class `GLNRep(lam)`: exact sparse gl(N)
GT generators `E[(i,j)]` as lists over columns of dicts {row: Fraction}, `pats`, `dim`, `top` = HW index, `wsp`, metric `g`, `float_ops()`, `dense(ij)`) and the
gl(4) slice engine. **Do not wait for it and do not touch its files.** Where you need gl(4) representation matrices, first check whether
`Paul/Python/gaudin_gl4/glN_reps.py` exists and `python -X utf8 test_glN_reps.py` passes; if not (or not yet), build the matrices yourself with a private
scaffold `_scaffold_rep(lam)` inside `gaudin_gl4_ops.py` from `Paul/Python/gl3_reps.py`'s `gt_patterns`, `_raise`, `_lower`, `_weight` — **these GT formulas are
general-N** (verified today: with E_{k,k+1} from `_raise(k,·)`, E_{k+1,k} from `_lower(k,·)` (targets outside the pattern set dropped, as `GL3Rep._adjacent_up/_down`
do), E_kk = diag(`_weight(k,·)`) and E_ij = [E_{i,j−1}, E_{j−1,j}], E_ji = [E_{j,j−1}, E_{j−1,i}] for |i−j| ≥ 2, all 256 relations [E_ij,E_kl] = δ_jk E_il − δ_li E_kj
hold exactly on [1,0,0,0], [1,1,0,0], [2,1,0,0], [2,1,1,0]; the HW pattern is the last index). Everything you compute (commutators, spectra, term lists) is
basis-independent, so the scaffold is harmless; say in RESULTS which one you used.

## 0. Why this exists

For gl(3), the multiplicity of V_μ in V_λ⊗V_λ′ is resolved by ONE operator, the two-site Gaudin "vertex operator" V = ½ Σ_{σ,τ∈S₃} sgn σ sgn τ e¹_{σ₁τ₁} e¹_{σ₂τ₂} e²_{σ₃τ₃}
(e^s_ij = E_ij on site s), which is — up to a μ-dependent constant — the double-pole Laurent coefficient a₁ of the cubic cdet coefficient S₃ at θ₁:
**θ₁₂·a₁ = V + η₁/2, η₁ = |Λ|² + cas₂(Λ) − cas₂(μ)/2** (exact on all 45 branches of the gl(3) reference data; `RESULTS-slice-engine.md`). V has simple spectrum on
every multiplicity space and its eigenvectors are the canonical copies. For gl(4) the missing-label problem needs (N−1)(N−2)/2 = 3 commuting operators; the
two-site gl(4) Gaudin algebra (the Laurent coefficients of the quantum characteristic polynomial) supplies exactly one mixed cubic and two mixed quartic
invariants — classically: S_k(u) = Σ_{a+b=k} P_{a,b} (u−θ₁)^{−a}(u−θ₂)^{−b} with P_{a,b} the polarised invariant with a site-1 and b site-2 generators; at
θ₁ the x^{−k}, x^{−(k−1)}, …, x^{−1} coefficients are triangular combinations of P_{k,0} (site-1 Casimir), P_{k−1,1}, …, P_{1,k−1}, and
Σ_a P_{a,k−a} is a combination of Casimirs, leaving k − 2 mixed invariants per degree k (0 for k = 2, 1 for k = 3, 2 for k = 4). The quantum operators
carry corrections (derivative terms in the cdet); **only the exact cdet combinations commute** — the naive trace invariants do not. This package builds
them exactly, verifies commutativity, fixes the canonical triple, and measures its joint spectrum on multiplicity spaces. WP-G3 then reduces the triple
to the PRV slice (normal-ordering the site-2 factors) using your term lists.

## 1. Conventions (binding — the gl(3) notes §1.1, `Paul/Mathematica/Gaudin/Gaudin-sl3-CG-notes.md`, generalised verbatim)

* Two sites, V_λ (site 1) ⊗ V_λ′ (site 2), generators E⁽¹⁾ = E⊗1, E⁽²⁾ = 1⊗E, [E_ij,E_kl] = δ_jk E_il − δ_li E_kj; inhomogeneities θ₁, θ₂, **θ₁₂ = θ₁ − θ₂ =: t**;
  **untwisted** (ζ = 0). Lax matrix **L(u)_ij = E⁽¹⁾_ij/(u−θ₁) + E⁽²⁾_ij/(u−θ₂)**.
* **Column-ordered quantum determinant with the derivative moved to the right** (∂f = f∂ + f′):
  **cdet(∂_u − L) = Σ_{σ∈S_N} sgn σ · M_{σ(1)1} M_{σ(2)2} ⋯ M_{σ(N)N}, M = ∂_u·1 − L(u), = ∂^N − S₁∂^{N−1} + S₂∂^{N−2} − ⋯ + (−1)^N S_N.**
  Operator products are kept in the order written. For N = 3, with l_ij = L(u)_ij, l′ = ∂_u l, l″ = ∂_u² l, the notes' eq. (4) is
  S₁ = l₁₁ + l₂₂ + l₃₃,
  S₂ = Σ_{i<j} (l_ii l_jj − l_ji l_ij) − l′₂₂ − 2 l′₃₃,
  S₃ = l″₃₃ − l₁₁l′₃₃ − l₂₂l′₃₃ − l′₂₂l₃₃ + l₃₁l′₁₃ + l₃₂l′₂₃ + l′₃₂l₂₃ + l₁₁l₂₂l₃₃ − l₁₁l₃₂l₂₃ − l₂₁l₁₂l₃₃ + l₂₁l₃₂l₁₃ + l₃₁l₁₂l₂₃ − l₃₁l₂₂l₁₃.
  **Your symbolic expansion must reproduce exactly this for N = 3** (as noncommutative polynomials; this is the convention test). The gl(3) notes verified
  [S_a(u), S_b(v)] = 0 for this convention; if (and only if) your N = 4 commutativity test (G2.2) fails with it, try the row-ordered variant
  Σ sgn σ M_{1σ(1)} ⋯ M_{Nσ(N)} and the "derivative moved to the left" variant, report which one commutes, and keep the one that does.
* Laurent data at θ₁: with **x = u − θ₁** (so u − θ₂ = x + t), S_k(u) = Σ_{j≥1} S_{k,j} x^{−j} + regular. We call S_{k,j} the Laurent coefficients at θ₁ (operators
  on V_λ⊗V_λ′ whose coefficients are Laurent polynomials in t). The θ₂ data follow by the site swap 1 ↔ 2, t → −t.
  gl(3) facts to reproduce (N = 3, [2,1,0]² and [3,1,0]²): S₂'s double pole S_{2,2} = the scalar c₂ = λ₂(λ₁+1) (3 and 4); the simple pole
  **H₁ = S_{2,1} = −(1/t) Σ_{i,k} E⁽¹⁾_ik E⁽²⁾_ki + |Λ|²/t** (eq. (7) at ζ = 0), i.e. on HW vectors of weight μ: H₁ = (|Λ|² − ½[cas₂(μ) − 2cas₂(Λ)])/t with
  cas₂(μ) = Σ_i μ_i(μ_i + 4 − 2i); S₃ has no triple pole on states (check whether S_{3,3} vanishes as an operator or only on HW vectors — report);
  **t·S_{3,2} = V + η₁/2 on every HW vector**, V the 36-term vertex operator (`slice_engine.VTERMS`, `apply_V_full`; its ordering: e¹_{σ₁τ₁} e¹_{σ₂τ₂} e²_{σ₃τ₃}
  as written, e¹ acting last), η₁ = |Λ|² + cas₂(Λ) − cas₂(μ)/2.
* GT basis and matrices as in G1 (`gt_patterns` order, HW last, columns = kets); two-site flat index ν·dim(V_λ′) + q, site 1 slow. For N = 3 the exact
  sparse matrices are `slice_engine.Rep.of(lam).E` (read-only import, `sys.path.insert(0, '../gaudin_sl3')`).
* Quadratic Casimir of gl(N): C₂(λ) = Σ_i λ_i(λ_i + N + 1 − 2i). Multiplicities: Racah–Speiser (G1's `racah_speiser`, or your own copy of
  `slice_engine.racah_speiser` generalised with ρ = (N−1,…,0)); HW vectors of weight μ: the null space of the total simple raising operators on the weight-μ
  subspace (Munich-style; small products only — that is all you need), orthonormalised.

## 2. The symbolic engine (`gaudin_gl4_ops.py`, general N)

Implement a small noncommutative algebra: an element is a dict {word: coefficient}, word = tuple of factor symbols, coefficient a sympy expression in t (or a
Fraction-valued Laurent polynomial in t — your choice, exact either way). Stage A: factor symbols ('l', i, j, d) = the d-th u-derivative of l_ij; elements of the
form Σ_k C_k ∂^k with C_k such polynomials; multiplication (A ∂^a)(B ∂^b) = Σ_j C(a,j) A·B^{(j)} ∂^{a+b−j} (Leibniz; B^{(j)} by the product rule, raising d
on one factor at a time). Expand the column-ordered cdet, read off S₁..S_N as polynomials in ('l',i,j,d). Stage B: substitute l_ij^{(d)} = (−1)^d d! [e¹_ij x^{−d−1} +
e²_ij (x+t)^{−d−1}], expand each word's scalar prefactor x^{−A}(x+t)^{−B} as a Laurent series in x with coefficients in ℚ[t^{−1}] ((x+t)^{−B} = Σ_n C(−B,n) t^{−B−n} x^n),
and collect the coefficients of x^{−j}, j = 1..N (and x⁰, for completeness): the operator words are now tuples of (site, i, j). Stage C: builders
`op_matrix(terms, rep1, rep2, t_value)` → sparse matrix on V_λ⊗V_λ′ (exact over Fractions via dict-of-columns products for the small products; scipy float64 CSR for
the larger ones), and `restrict(op, basis)` → m×m matrix on an orthonormal HW basis. Keep every term list exact; never simplify numerically.

## 3. Tasks and acceptance (every number goes into `RESULTS-G2-gaudin-gl4-ops.md`)

**G2.1 Convention test.** N = 3 symbolic S₁, S₂, S₃ = eq. (4) exactly (print both; assert equality as polynomials). Then the gl(3) numerical facts of §1 on [2,1,0]²
(dim 64) and [3,1,0]² (dim 225), exact arithmetic: all S_{2,j}, S_{3,j} (θ₁ and θ₂) commute pairwise and with every E^{tot}_ij; S_{2,2} = λ₂(λ₁+1)·1; H₁ as stated;
t·S_{3,2} = V + η₁/2 on every HW vector of every irrep (use `slice_engine.hw_from_slice_exact` / `reference_irrep`, or your own null space); report the
triple-pole question. **Acceptance: exact equalities, 0 violations.**

**G2.2 Commutativity for N = 4.** All Laurent coefficients S_{k,j}, k = 2,3,4, j = 1..k, at θ₁ AND at θ₂, commute pairwise and with every E^{tot}_ij, exactly, on
[1,0,0,0]⊗[1,1,0,0] (dim 24) and [1,1,0,0]⊗[2,1,0,0] (dim 120); on [2,1,0,0]² (dim 400) and [2,1,1,0]² (dim 225) exactly if affordable, else float64 at 1e-12
relative, at two rational values of t (31/35 and −3/7). Report the scalar double pole S_{2,2} for gl(4) as a function of λ (guess from the data over several λ:
candidates λ₂(λ₁+1) + λ₃(λ₁+λ₂+2)+…, i.e. the quantum-corrected "e₂" — state the formula you find and verify on ≥ 6 reps) and whether S₃, S₄ have higher poles
as operators. **Acceptance: 0 violations** (if the column-ordered cdet fails, see §1 for the fallback; then report which convention is right).

**G2.3 The canonical triple.** Define **V₃ := t·S_{3,2}, V₄ᵃ := t·S_{4,3}, V₄ᵇ := t²·S_{4,2}** (classically P_{2,1}, P_{3,1}, P_{2,2} − P_{3,1}; see §0). On the HW
space of each irrep μ (m ≥ 2) of [2,1,0,0]², [2,1,1,0]², [3,2,1,0]² (dim 64² = 4096), [2,2,1,0]⊗[2,1,1,0] (60·15 = 900), and [3,2,0,0]², [2,2,1,0]² if m ≥ 3
sectors are scarce (compute Racah–Speiser first and print the multiplicity table of each product): restrict V₃, V₄ᵃ, V₄ᵇ (exact m×m in an orthonormal HW basis;
eigenvalues with mpmath ≥ 30 digits); check (i) they commute on the HW space (exact), (ii) the joint spectrum is simple on every sector (min distance between
distinct joint eigenvalue triples), (iii) whether V₃ ALONE has simple spectrum on every sector (min gap), (iv) whether the spectra are real and t-independent
(two t values; if the eigenvalues depend on t, the t-dependent part must be a scalar on the sector or a Casimir — report what you see), (v) that every other
Laurent coefficient (θ₁ and θ₂, S₂..S₄) is diagonal in the joint eigenbasis (so the triple generates everything). **Acceptance:** a table per product: μ, m,
V₃-spectrum simple?, joint simple?, min gaps; **if the joint spectrum is degenerate anywhere, report the sector and the degeneracy and STOP that line of
work** (the fallback — twisted/shift-of-argument operators — is a decision for Paul; do not invent one).

**G2.4 Structure of V₃.** Compare V₃ with the natural cubic mixed invariants: T^{(1,1,2)} = Σ_{i,j,k} e¹_{ji} e¹_{kj} e²_{ik} and T^{(1,2,2)} = Σ e¹_{ji} e²_{kj} e²_{ik}
(CPV's polarised traces), and with the gl(3)-style ½ Σ_{σ,τ∈S₄} sgn σ sgn τ e¹_{σ₁τ₁} e¹_{σ₂τ₂} e²_{σ₃τ₃} δ_{σ₄τ₄}-type antisymmetrised form if that is what the
cdet produces: give V₃ − (a·T^{(1,1,2)} + b·T^{(1,2,2)}) as an explicit polynomial and state whether the remainder is a combination of Casimirs/scalars on HW
spaces (numerically on the test sectors: is the remainder restricted to each M_μ a multiple of the identity?). Same for the quartic pair against the degree-4
polarised traces T^{(1,1,1,2)}, T^{(1,1,2,2)}, T^{(1,2,1,2)}, T^{(1,2,2,2)}. For N = 3 the answer must be V₃ = V + η₁/2.

**G2.5 Term lists** `Paul/Mathematica/Gaudin/ref/gl4-gaudin-ops.txt`: for N = 4, the exact Laurent coefficients S_{k,j} (k = 2,3,4; j = 1..k; θ₁) and the triple
V₃, V₄ᵃ, V₄ᵇ as term lists: one line per term, TAB-separated: `OP name`, then `TERM coefficient word` with coefficient an exact Laurent polynomial in t written
as `{{power,p/q},...}` and word as `{{site,i,j},...}` in product order (leftmost factor acts last). `#` header: conventions (§1 verbatim in short), the definition
of x and t, the site swap rule, the N = 3 check, counts. Also a loader `load_terms(path)` in your module and a re-verification that the loaded V₃ reproduces the
matrices (commutativity + spectrum on one sector). Say in the header that e¹ and e² commute, so G3 may move all site-2 factors to the right of the site-1
factors keeping each site's internal order.

**G2.6 `RESULTS-G2-gaudin-gl4-ops.md`**: files, conventions, the convention test, the commutativity results, the double-pole scalar, the multiplicity tables of
the test products, the simple-spectrum table, the structure of V₃ and the quartic pair, term counts, timings, surprises, what is not done. Output of all runs in
`gaudin_gl4_ops_output.txt`. CLI: `python -X utf8 gaudin_gl4_ops.py <stage>` with stages `sym | gl3 | comm | triple | struct | terms | all`. Report failures with
numbers; do not paper over a failing test.

## 4. Pitfalls

* The ordering inside each term matters (operators do not commute within a site); keep the cdet order; e¹ and e² factors commute across sites.
* Keep t symbolic in the term lists; substitute only in the numerical builders. 1/t terms from the derivative corrections are expected.
* Fractions through 4-fold sparse products at dim 400 are fine with dict-of-columns arithmetic (`slice_engine._compose` style); avoid dense sympy matrices above dim ~100.
* HW spaces must be computed in a basis where the restriction is faithful: orthonormalise the null-space basis in the ORTHONORMAL GT basis (metric √g per site) or
  use the plain basis with the invariant form G⊗G — either is fine, state which; the spectra do not depend on it.
* The gl(3) vertex operator's ordering convention (e¹ e¹ e² with the first factor acting last) is the one `apply_V_full` implements — reuse it for the N = 3 check
  instead of re-deriving.
