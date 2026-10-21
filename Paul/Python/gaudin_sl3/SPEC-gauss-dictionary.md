# Spec D: the Gauss dictionary for general Λ — ĉ, φ_u = e₁ᵀρ_Λ(n′), the PRV slice space M_μ, and the (m+1)-point dictionary (round 6, WP-D; 2026-10-21)

Delegated to a general-purpose subagent (SONNET), Python only (no Wolfram, no wolfbook tools, no kernel). This brief is derived from
`Paul/Mathematica/Gaudin/Gaudin-sl3-round6-plan.md` (§0–§2: your ingredient is row (e) of the §1 table; the formula block at the top of §1 defines
φ_u, n′, ĉ and the identity; §2.1 defines the PRV slice space M_μ, §2.2 step (e) the m+1-point dictionary; §3 row D; §5 risk D; §6 conventions)
and from the round-4/5 sources listed in §4. It must be readable without them. **Files you own:** `gauss_dictionary.py`, `RESULTS-gauss-dictionary.md`
(this file is also yours; keep it accurate), optional cache `gauss_dictionary_*.pkl` (gitignored, rebuildable). Existing scripts are READ-ONLY (import them, never
edit); no notebooks, no git commits/pushes; ≤ 2 heavy processes; no n_W-sized object anywhere (never build the two-site space, never `Chain(...)` for
a large Λ — the covector lives in the one-site space of dimension dim V_Λ = 512 at [14,7,0]). Environment: Python 3.13 (`python -X utf8`), sympy 1.14,
mpmath 1.3, numpy 2.x, scipy; run from `Paul/Python/gaudin_sl3/`; import pattern of `gaudin_chain.py` (`sys.path` gets the script dir and `..` for `gl3_reps`).

## 0. Where this sits, and what the dictionary is

The Q-route (plan §1) ends in the identity (constant exactly 1 on all 26 branches of [2,1,0], [3,1,0], [3,2,0]; N1, W5, S1)

    L_B(ω̂₁, ω̂₁′) := Σ_{d+n+n′=N_B} G^{[d]}_ij(ω̂₁,ω̂₁′;θ₁₂) ν_i^{(n)}(θ₁) ν_j^{(n′)}(θ₂)  =  ĉ(φ,ω̂) · ⟨φ_u(ω̂)⊗e₁|Ψ₀⟩ / X_B ,

one polynomial in (ω̂₁, ω̂₁′) per branch B = (μ, copy) (degree ≤ 2λ₁−λ₂ in each variable), X_B one unknown constant per branch (drops out of the
direction). The right-hand side is a **generating function of the slice**: with C_ν := ⟨e_ν⊗e₁|Ψ₀⟩ (the coefficient of the HW vector Ψ₀ ∈ V_Λ⊗V_Λ at (site-1
GT index ν, site-2 index = the HW state e₁)),

    ⟨φ_u(ω̂)⊗e₁|Ψ₀⟩ = Σ_ν [φ_u(ω̂)]_ν C_ν ,    φ_u(ω̂) := e₁ᵀ ρ_Λ(n′(ω̂)) ∈ V_Λ^* ,

and only ν in the weight space W_μ := V_Λ[μ−Λ] contribute (Ψ₀ has weight μ, site 2 has weight Λ). W_μ has dimension w (≤ 8 at [14,7,0]² near
a = 7, ≤ 15 at a = 14); the slice s_B = (C_ν)_{ν∈W_μ}/X_B is a vector of w numbers, and it lives in the m-dimensional PRV space (m = multiplicity of μ in V_Λ⊗V_Λ).
**D's job: turn the polynomial L_B into the slice, cheaply and without any n_W-sized object**: m+1 evaluations of L_B at generic points ω̂_k, an m×m
(least-squares (m+1)×m) linear solve, no overdetermined polynomial fit (N1/U2 solved the identity coefficientwise on 64 grid points — that route needs 64 numbers and
the HW space H_μ of dimension n_W-ish; the dictionary needs m+1 numbers and M_μ ⊂ W_μ).

### 0.1 Definitions (binding)

* Λ = [λ₁,λ₂,0] (the code may accept λ₃ ≠ 0 for the covector, but ĉ is only defined for λ₃ = 0). One-site GT basis = **`gl3_reps.GL3Rep(Λ)` order** ("my" order;
  patterns = (row₁, row₂, Λ), rows increasing lexicographically in (m₁,m₂,k)); the HW state is `hw := index of ((λ₁),(λ₁,λ₂),Λ)` (the LAST index in this order). The
  basis is **non-orthonormal** (invariant diagonal metric G; do not use the G-norm in D). W4's `ref/untwisted-exact-<Λ>.txt` use the notebook order
  (`Reverse[Sort[patterns]]`; HW first); `nu_general.w4_to_mine(lam, w4, ch)` maps W4 one-site index (0-based) → my index (pattern-wise).
* Generators: E₁₂ (k=1) and E₂₃ (k=2) are the raising matrices of `gl3_reps` (`_raise(k, pattern)`); **E₁₃ := [E₁₂, E₂₃]** (as in `gaudin_chain.Chain.__init__`; assert
  [E₁₃,E₃₁] = E₁₁−E₃₃ with E₃₁ := [E₃₂,E₂₁] on the small reps). Matrices act on kets as columns, `M[row,col]`; covectors act from the right: (e₁ᵀM)_c = M[hw,c]
  (`gaudin_chain.SM.rapply`). Weights = (E₁₁,E₂₂,E₃₃) eigenvalues: E₁₂ raises the weight by α₁ = (1,−1,0), E₂₃ by α₂ = (0,1,−1), E₁₃ by α₁+α₂; a component ν of φ_u with Λ−wt(ν) = n₁α₁ + n₂α₂
  has exactly n₁ = p+q, n₂ = q+r in every term a′ᵖg′^q c′^r of its expansion.
* n′ = exp(a′E₁₂) exp(g′E₁₃) exp(c′E₂₃), with (round-4 plan §0; `nu_general.Rep.phi_u`, `untwist_limit.gauss_closed`)

      a′ = (φ₁−φ₂)(ω̂₁−ω̂₁′) / ((φ₂−ω̂₁)(φ₂−ω̂₁′)),
      c′ = (φ₁−φ₃)(φ₂−φ₃)(ω̂₁−ω̂₁′) / ((φ₁−φ₂)(φ₃−ω̂₁)(φ₃−ω̂₁′)),
      g′ = (φ₃−φ₁)(φ₃−φ₂)(ω̂₁−ω̂₁′) / ((φ₃−ω̂₁)(φ₂−ω̂₁′)(φ₃−ω̂₁′)).

  These are **Λ-independent group data** (ω̂₁ ≡ `w1`, ω̂₁′ ≡ `w1p`); only ρ_Λ changes. φ_u is a function of (φ; ω̂₁, ω̂₁′) and of Λ through the nilpotent exponentials
  only (E₁₂ⁿ = 0 for n > λ₁−λ₃, ...). Nice side remark to verify: the 3×3 unipotent matrix n′ has entries x₁₂ = a′, x₂₃ = c′ and x₁₃ = g′ + a′c′ =
  (φ₃−φ₁)(φ₃−φ₂)(ω̂₁−ω̂₁′)/((φ₃−ω̂₁)(φ₃−ω̂₁′)(φ₂−ω̂₁)).
* ĉ(Λ; φ; ω̂) = (φ₁−φ₂)^{2λ₂} · Π_{ω∈{ω̂₁,ω̂₁′}} (φ₂−ω)^{λ₁−λ₂}(φ₃−ω)^{λ₁} (RESULTS-nu-general §0; adjoint: `untwist_limit.chat_expr`).
* PRV slice space (plan §2.1; Zhelobenko/PRV): **M_μ = { v ∈ V_Λ[μ−Λ] : E₁₂^{λ₁−λ₂+1} v = 0, E₂₃^{λ₂−λ₃+1} v = 0 }**, dim M_μ = m = mult of μ in V_Λ⊗V_Λ (Racah–Speiser:
  m = Σ_{σ∈S₃} sgn σ · dim V_Λ[σ·μ − Λ], σ·μ = σ(μ+ρ)−ρ, ρ = (2,1,0), dim V_Λ[wt] = #GT patterns of weight wt). v ↦ Ψ₀(v) is the inverse of the slice map: the exact slices
  of every HW vector of weight μ (in particular W4's `cg` records and the slices of `HWBASIS` kets) lie in M_μ and span it.
* The dictionary equation (per branch): for m+1 or more generic points ω̂_k = (w1_k, w1p_k): **L_B(ω̂_k)/ĉ(ω̂_k) = Σ_j x_j ⟨φ_u(ω̂_k)|b_j⟩**, b_j a basis of M_μ (unknowns x_j, j = 1..m; their overall
  scale absorbs 1/X_B), solved by least squares; the residual of the redundant equation(s) is the consistency check; slice := Σ_j x_j b_j, normalised (first nonzero entry in my order = 1; the
  full-vector normalisation/sign is a separate convention, plan §0 row (g)).

## 1. Tasks (all Python; sympy for exact, mpmath for numeric)

**D1 — `gauss_covector` and `chat`.** (a) `rep1(lam)`: one-site sparse exact generators from `gl3_reps` (build E₁₂, E₂₃ from `gl3_reps._raise`, E₁₃ by the sparse commutator — NOT the dense sympy `GL3Rep.generator`,
which is hopeless at dim 512; you may reuse `gaudin_chain.SM`/`comm`). (b) `gauss_covector(lam, phi, w1, w1p)` = e₁ᵀρ_Λ(n′) as a list over my GT index, **exact `Fraction` for rational inputs, `mpf` for numeric inputs**
(nilpotent exponentials applied to the covector; also a Λ-only tensor T[ν;(p,q,r)] = ⟨hw|E₁₂^p E₁₃^q E₂₃^r|e_ν⟩/(p!q!r!) so that φ_ν = Σ T a′^p g′^q c′^r — your choice which route is the default, but implement both
and assert they agree exactly). `chat(lam, phi, w1, w1p)`. (c) `gauss_covector_poly(lam, phi)` → (list over ν of the sympy `Poly` in (W1,W2) of ĉ·φ_ν, or `None` where ĉ does not clear the
denominator; ĉ as a Poly) and **reproduce `untwist_limit.gf_polys(Chain((2,1,0)), phi)` exactly** (same polynomials, same `None` pattern; both φ = (1/3,−2/7,5/11), φ′ = (−2/5,3/8,7/11); the adjoint gf_polys uses gl3_reps order, so no map is needed).
(d) U2's "Surprise 1" generalised: for [2,1,0], [3,1,0], [3,2,0] list the components ν where ĉφ_ν is NOT polynomial (GT pattern, weight, which factor (φ₂−ω)/(φ₃−ω) is left over), explain the pattern (it should be a pole
order condition in terms of the weight difference Λ−ν = n₁α₁ + n₂α₂), and **assert that every such component vanishes on W4's HW kets** (all `HWBASIS` kets and all `BRANCH` vectors of all sectors, slice
at site-2 = HW index; load with `nu_general.load_w4` / `hw_kets` or parse yourself); also assert that no component with ν ∈ W_μ, M_μ ≠ 0, is left uncleared. (e) cross-check `gauss_covector` exactly against `nu_general.Rep.phi_u`
(N1; Chain-based) at a few points for the three reps. (f) **timing at [14,7,0]** (dim 512): `rep1` setup, `chat`, `gauss_covector` per point (exact and mp 50 digits), and the table for Λ = [2a,a,0], a = 1..14. Target: < 1 s per point; report.

**D2 — the PRV space.** `prv_basis(lam, mu)` = exact basis of M_μ (list of sparse vectors `{my GT index: Fraction}`, RREF nullspace form, free columns in increasing my-index order; document it — agent R builds the same
space, you are independent), `multiplicity(lam, mu)` (Racah–Speiser from the GT weight counts), `irreps(lam)` (all dominant μ with m ≥ 1). Checks: dim M_μ = Racah–Speiser m **for every irrep of [2,1,0], [3,1,0], [3,2,0], [4,2,0],
[6,3,0], [8,4,0]** (and [14,7,0]: all 113 irreps if the time is reasonable; report the time per irrep), = W4's `SECTOR … mult` for the three W4 reps, Σ_μ m·dim V_μ = dim V_Λ² (dim V_μ by the Weyl formula), and **W4's exact slices (`cg` of every BRANCH,
slices of every `HWBASIS` ket) satisfy the PRV conditions and span M_μ** (rank = m, exact).

**D3 — the (m+1)-point dictionary.** `slice_from_points(lam, mu, phi, points, values)`: build A[k][j] = Σ_ν [φ_u(ω̂_k)]_ν b_j[ν] (only ν ∈ W_μ, i.e. w components per point), y_k = values[k]/ĉ(ω̂_k), solve A x = y by exact-structure least squares in mpmath
(QR; ≥ m+1 points; report the relative residual ‖Ax−y‖/‖y‖ and cond(A)), slice = Σ x_j b_j normalised. **Validate on all 26 known branches** (adjoint 6, [3,1,0] 10, [3,2,0] 10; both φ and φ′ — N1's records
`nu_general_res_<lam>.pkl` (φ) and `nu_general_res_<lam>2.pkl` (φ′): a dict with `br`, a list of records with `mu`, `copy`, `LB` = {(a,b): mpf} the leading polynomial Σ_{(a,b)} c_ab ω̂₁^a ω̂₁′^b, `slex` = the exact slice in my index, `phi`; load with `pickle` after `import mpmath`):
use L_B from the records (state it) and, as an independent self-consistency test, L from the exact slices through the identity itself (say which is which); m+1 points, then m+2, 2m+2; recover all 26 slices
**to ≥ 30 digits** (compare with W4's exact slice after fixing the scale on the first nonzero entry in my order; `nu_general.w4_cg(R, mu, copy)` gives it), **multiplicity-2 copies included** (the two copies of a sector must come out as two
different vectors of the 2-dimensional M_μ, each proportional to its own exact slice). Then: (i) conditioning study — cond(A) vs m for the real multiplicities, for naive vs well-chosen points; digits lost; (ii) noise study — perturb the values by
relative noise 1e-10, 1e-20, 1e-40 and report the slice error (amplification factor); (iii) scale test with SYNTHETIC data (clearly labelled: random exact element of M_μ → L(ω̂_k) = ĉ·⟨φ_u|s⟩) at Λ = [4,2,0], [6,3,0], [8,4,0], [14,7,0], every irrep: recovery error
(exact data → exact recovery; mp data at 50 digits → error), time per irrep, worst cond(A), worst m, worst w. Also (iv) a held-out check (the recovered slice must predict L_B at N1's 64-point grid through the identity), (v) a negative control
(L_B of one irrep fed to another irrep's dictionary: the residual must separate right from wrong), (vi) which points to use at large m (naive first m+1 points vs a per-irrep greedy choice vs ONE shared nested list — worst cond(A) over all irreps), and
(vii) once agent R's `ref/cg-slices-<Λ>.txt` exist: compare the PRV bases exactly (both use the RREF nullspace form, free columns in increasing GT index) and run the dictionary on R's V_red eigen-copies (the real SoV copies) at [4,2,0] … [14,7,0].

**D4 — the interface round 7 will call** (documented in the module docstring and in RESULTS):

    gauss_dictionary.slice_from_points(lam, mu, phi, points, values, dps=None) -> (x, residual, slice_vector)

lam = (λ₁,λ₂,0), mu = (μ₁,μ₂,μ₃) a dominant weight with M_μ ≠ 0, phi = (φ₁,φ₂,φ₃) (`Fraction` or `mpf`; the same φ as in the ε-regime of the L_B data), points = list of ≥ m+1 pairs (w1_k, w1p_k) (`Fraction`/int
exact, or `mpf`; w1 ≠ w1p, never equal to a φ_i), values = list of L_B(w1_k, w1p_k) (mpf/Fraction/float, one per point). Returns x (list of m mp numbers: coordinates in `prv_basis(lam, mu)` of RAW = Σ x_j b_j = C/X_B, i.e. the overall scale 1/X_B is
kept), residual (mpf, relative, row-scaled: ‖D(Ax−y)‖/‖Dy‖, D = 1/‖A_k‖), slice_vector (list of length dim V_Λ in my GT order, mp numbers, zero outside W_μ, first nonzero entry (|·| > 10^(−dps/2) max|·|) = 1). Also provide:
`slices_from_points(...)` (all copies of one irrep at once), `DictionaryProblem(lam, mu, phi, points).solve(values)` (one matrix for many right-hand sides; `.cond()`), `slice_from_points_ex(...) -> dict` (A, c-hat, singular values, cond, raw slice, pivot, exact path,
timing), `default_points(phi, n)`, `shared_points(lam, phi)` (one nested ordered list for the whole Λ; use `[:m+1]`) and `choose_points(lam, mu, phi, n)`, `chat`, `gauss_params`, `gauss_covector`, `gauss_covector_poly`, `uncleared_components`, `prv_basis`,
`prv_free_indices`/`coords_of` (coordinates x_j = v[f_j]), `multiplicity`, `irreps`, `synthetic_values` (the identity itself). Per-branch cost must be O(m·w·(#terms)) after a once-per-Λ setup; nothing may scale with n_W.

## 2. Acceptance (numbers to report in RESULTS)

1. `gauss_covector_poly((2,1,0), φ)` ≡ `untwist_limit.gf_polys` exactly (8 components, both φ).
2. Non-cleared components listed and asserted to vanish on W4's HW kets for [2,1,0], [3,1,0], [3,2,0] (max |C_ν| = 0 exactly: rational data).
3. dim M_μ = multiplicity for all irreps of the six small reps (exact); W4 slices ∈ M_μ and spanning.
4. 26/26 branches recovered to ≥ 30 digits (target: ≈ 90 digits, the accuracy of N1's L_B), residual of the redundant equation reported; both φ.
5. Timing at [14,7,0]: ĉ and φ_u per point < 1 s (exact and mp); scaling table a = 1..14; synthetic recovery over all 113 irreps of [14,7,0].

## 3. Conventions, rules, traps

* Run from `Paul/Python/gaudin_sl3`; no edits to existing files; heavy imports (`nu_general` 0.5 s, `untwist_limit` 0.5 s) only inside the validation functions, so that `import gauss_dictionary` stays light (needs only `gl3_reps`, `gaudin_chain.SM` optional, sympy, mpmath).
* The GT basis is non-orthonormal; compare vectors up to scale, slices componentwise after fixing the scale on the first nonzero entry (in my index order — state the order whenever a normalised vector is printed). Sign and G-norm are not available without the full vector and are not part of D.
* Memoised helpers: if you cache by (lam, ...) keys make the key include everything that matters (φ, point, number type, mp.dps for mp results).
* mpmath precision is global: set `mp.mp.dps` explicitly in every entry point that returns mp numbers (or use `mp.workdps`), and never mix a 200-digit pickle value with a lower-precision context without saying so.
* Failures are reported with numbers, not adjectives: if a conditioning problem appears at large m say what m, what cond, how many digits are lost.

## 4. Sources (all read-only)

`Paul/Mathematica/Gaudin/Gaudin-sl3-round6-plan.md` §0–2, 3, 5, 6; `…/Gaudin-sl3-round4-plan.md` §0 (n′, Gauss parameters, jet formula, the canonical kernel); `RESULTS-nu-general.md` (N1: §0 ĉ, §2 the 20 exact slices of [3,1,0]/[3,2,0] with GT indices, §2 (iv) the dictionary
by overdetermined fit — the thing D replaces by m+1 points), `RESULTS-untwist-limit.md` (U2: the adjoint, Surprise 1), `untwist_limit.py` (`gf_polys(ch, phi)`, `gauss_closed`, `chat_expr`, `expm_nil`, `sym_E`), `nu_general.py` (`Rep.phi_u`, `Rep.chat`, `load_w4`, `w4_to_mine`,
`w4_cg`, `hw_kets`, `analyse`), `gaudin_chain.py` (`SM`, `comm`, `exp_nil`, `Chain.__init__` for E₁₃ = [E₁₂,E₂₃] and the import pattern), `Paul/Python/gl3_reps.py` (`GL3Rep`, `_raise`, `_lower`, `_weight`), W4's `Paul/Mathematica/Gaudin/ref/untwisted-exact-{210,310,320}.txt`
(header documents every record), N1's regime caches `nu_general_res_<lam>[2].pkl`.

## 5. Reporting

`RESULTS-gauss-dictionary.md` in the style of the round-5 RESULTS files: files and how to run, conventions, headline numbers, the tables (non-cleared components, PRV dimensions, the 26-branch recovery table with m, #points, residual, digits, conditioning, noise amplification, scaling a = 1..14, synthetic [14,7,0]), surprises, not done. End with a concise summary.

## 6. As built (2026-10-21)

`gauss_dictionary.py` (CLI `python -X utf8 gauss_dictionary.py [--stages polys,prv,timing,branches,cond,points,scale,misc,crossR]`, ~2.5 min; printout `gauss_dictionary_output.txt`) implements all of the above; numbers and discussion in `RESULTS-gauss-dictionary.md`.
Deviations from this brief: none in the interface; the optional pkl caches were not needed (everything is cheap); the default point rule is `shared_points(...)[:m+1]` rather than the naive first m+1 `default_points` (the latter loses up to 9 digits at m = 15).
