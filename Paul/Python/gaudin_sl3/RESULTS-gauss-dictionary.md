# RESULTS D: the Gauss dictionary for general Λ — ĉ, φ_u = e₁ᵀρ_Λ(n′), the PRV slice space M_μ and the (m+1)-point dictionary (agent D, round 6, 2026-10-21)

Files: `gauss_dictionary.py` (module + CLI: `python -X utf8 gauss_dictionary.py [--stages polys,prv,timing,branches,cond,points,scale,misc,crossR] [--amax 14]`, ~2.5 min in total; **its full printout is `gauss_dictionary_output.txt`**, every number below is from it),
`SPEC-gauss-dictionary.md` (the brief, updated to the as-built interface). Python only (sympy exact, mpmath numeric, numpy/scipy only for the point choice); no caches written (everything is cheap), no two-site object anywhere (the largest object is a dim V_Λ vector: 512 at [14,7,0], 3375 at [28,14,0]).
Read-only inputs: `gl3_reps.py` (`GL3Rep`, `_raise`, `_lower`, `_weight`), `gaudin_chain.py` (`SM`, `comm`), and — validation stages only, lazily imported — `untwist_limit.py` (`gf_polys`), `nu_general.py` (`Rep.phi_u`, `load_w4`, `w4_cg`, `hw_kets`, `GRID8_*`),
N1's records `nu_general_res_<lam>[2].pkl` (L_B, exact slices, X_B), W4's `ref/untwisted-exact-{210,310,320}.txt`, and WP-R's `ref/cg-slices-{420,630,840,1470}.txt` (cross-check).

## 0. Conventions and what the module does

* One-site GT basis = `gl3_reps.GL3Rep` order ("my" order; HW = last index; non-orthonormal). E₁₂, E₂₃ = the `gl3_reps` raising matrices, **E₁₃ := [E₁₂,E₂₃]** (checked: [E₁₃,E₃₁] = E₁₁−E₃₃ and [E₁₂,E₂₁] = E₁₁−E₂₂ on [2,1,0], [3,1,0], [3,2,0], [4,2,0], [8,4,0] and on the λ₃ > 0 smoke tests); sparse exact (`Fraction`) matrices built straight from `_raise`, never the dense sympy generators.
* `gauss_params(phi, w1, w1p)` = (a′,g′,c′) of the round-4 plan; `chat` = ĉ (λ₃ = 0); `gauss_covector` = e₁ᵀρ_Λ(n′): exact `Fraction`s for rational input, `mpf` for numeric input. Default route: the Λ-only tensor T[ν;(p,q,r)] = ⟨hw|E₁₂ᵖE₁₃^qE₂₃^r|e_ν⟩/(p!q!r!)
  (φ_ν = Σ T a′ᵖg′^q c′^r; once per Λ, 0.010 s at [14,7,0]); second route: direct nilpotent exponentials of the covector. Third (check): the other Gauss ordering e₁ᵀexp(x₂₃E₂₃)exp(x₁₃E₁₃)exp(x₁₂E₁₂) with the matrix entries x₁₂ = a′, x₂₃ = c′, **x₁₃ = g′+a′c′ = (φ₃−φ₁)(φ₃−φ₂)(ω̂₁−ω̂₁′)/((φ₃−ω̂₁)(φ₃−ω̂₁′)(φ₂−ω̂₁))**
  (no cross terms: E₂₃E₁₂ = E₂₃E₁₃ = E₁₃E₁₂ = 0 in the 3×3 matrices). All three agree **exactly** (Fractions) on [2,1,0], [3,1,0], [3,2,0], [4,2,0], [8,4,0], [14,7,0] at 3 points each; `gauss_covector` equals N1's Chain-based `Rep.phi_u` exactly (0.0 difference, 12 (φ, point) combinations × 3 reps).
* PRV: `prv_basis(lam, mu)` = exact RREF-nullspace basis of M_μ (free columns in increasing my-index order, b_j[f_j] = 1, zero at the other free indices) — **so the coordinates of any v ∈ M_μ are simply x_j = v[f_j]** (`prv_free_indices`, `coords_of`). `multiplicity` = Racah–Speiser from the GT weight counts, `irreps(lam)` = all dominant μ with m ≥ 1.
* Dictionary: L_B(ω̂_k)/ĉ(ω̂_k) = Σ_j x_j ⟨φ_u(ω̂_k)|b_j⟩ (sum over W_μ only: w components per point), least squares in mpmath (QR after scaling each equation by 1/‖A_k‖; residual = ‖D(Ax−y)‖/‖Dy‖), exact rational least squares when phi, points and values are all rational. **RAW := Σ x_j b_j = C/X_B** (scale kept); the returned `slice_vector` is RAW scaled so that its first nonzero entry (my order) is 1 —
  the same normalisation as WP-R's COPY records. Interface: `slice_from_points(lam, mu, phi, points, values, dps=None) -> (x, residual, slice_vector)`; also `slices_from_points` (all copies of an irrep, one matrix), `DictionaryProblem(...).solve/.cond`, `slice_from_points_ex`, `default_points`, `shared_points`, `choose_points`, `synthetic_values` (the identity itself; module docstring has the details).
  `dps` defaults to max(mp.mp.dps, 30): **pass dps ≥ the digits of your data** (the data precision is not inferred).

## 1. Headline results

1. **`untwist_limit.gf_polys` reproduced exactly** for the adjoint (8 components, sympy `Poly` equality component by component, ĉ equal; both φ = (1/3,−2/7,5/11) and φ′ = (−2/5,3/8,7/11)): 7 polynomials and `None` on ν = 3 (weight (0,2,1)) in both. ĉφ_ν has bidegree ≤ 2λ₁−λ₂ in (ω̂₁,ω̂₁′) on every cleared component ([2,1,0]: 3; [3,1,0]: 5; [3,2,0]: 4; [4,2,0]: 6), as N1 found for L_B.
2. **U2's Surprise 1, generalised and explained.** ĉ fails to clear φ_ν only through a pole **at ω̂₁′ = φ₂ and nowhere else** (sympy: the leftover denominator is a pure power of (ω̂₁′−φ₂) on every uncleared component of 7 reps × 2 φ; a fast exact pole-order routine from the tensor T agrees with sympy in every case). Reason: with Λ−wt(ν) = n₁α₁+n₂α₂, every term of φ_ν carries (φ₂−ω̂₁′)^{−(p+q)} = (φ₂−ω̂₁′)^{−n₁} against the (φ₂−ω̂₁′)^{λ₁−λ₂} of ĉ,
   all other exponents are ≥ 0 automatically (p ≤ λ₁−λ₂ because ⟨hw|E₁₂ᵖ = 0 beyond the string; n₂ ≤ λ₁−λ₃), so the leftover pole has order n₁−(λ₁−λ₂) minus the vanishing order of the numerator at ω̂₁′ = φ₂. **Why the slices never see it:** in the other Gauss ordering x₁₂ = a′ is the ONLY entry with a pole at ω̂₁′ = φ₂ and ⟨hw|ρ(n′)|v⟩ has x₁₂-degree ≤ λ₁−λ₂ whenever E₁₂^{λ₁−λ₂+1}v = 0 (E₁₂ acts first on v in that ordering), which is the PRV condition — so on M_μ the pole is cleared by ĉ itself.
   Checked: the uncleared sets (ν, pole order) are φ-independent, and **their intersection with the union of the PRV supports of all irreps is empty for every Λ tested, up to [28,14,0]** (560 uncleared components of 3375, none in the 2507-component PRV support):

   | Λ | dim | #uncleared (max order) | #ν with n₁>λ₁−λ₂ | PRV support | ∩ |
   |---|---|---|---|---|---|
   | [2,1,0] | 8 | 1 (1) | 2 | 6 | 0 |
   | [3,1,0] | 15 | 1 (1) | 3 | 11 | 0 |
   | [3,2,0] | 15 | 4 (2) | 6 | 10 | 0 |
   | [4,2,0] | 27 | 4 (2) | 9 | 20 | 0 |
   | [6,3,0] | 64 | 10 (3) | 24 | 47 | 0 |
   | [8,4,0] | 125 | 20 (4) | 50 | 92 | 0 |
   | [14,7,0] | 512 | 84 (7) | 224 | 378 | 0 |
   | [28,14,0] | 3375 | 560 (14) | 1575 | 2507 | 0 |

   The uncleared ν themselves (my index, GT pattern k\|m₁,m₂, weight, pole order): [2,1,0]: 3 (0\|2,0; (0,2,1); 1). [3,1,0]: 8 (0\|3,0; (0,3,1); 1). [3,2,0]: 0 (0\|2,0; (0,2,3); 1), 6 (0\|3,0; (0,3,2); 2), 7 (1\|3,0; (1,2,2); 1), 10 (1\|3,1; (1,3,1); 1). **Asserted: C_ν = 0 on every W4 HW ket** (all HWBASIS kets and all BRANCH slices of the three reps, both φ; max |C_ν| = 0.0, rational data) —
   note ν = 7 of [3,2,0] sits in the weight space (1,2,2) of μ = (4,4,2): that space has two GT states, the PRV vector lives on the cleared one.
3. **PRV space.** dim M_μ = Racah–Speiser m for every irrep of [2,1,0], [3,1,0], [3,2,0], [4,2,0], [6,3,0], [8,4,0], **[14,7,0] (113 irreps, 344 copies, max m = max w = 8) and [28,14,0] (421 irreps, 2255 copies, max m = max w = 15)**; Σ m dim V_μ = dim V_Λ² in all eight; PRV bases for ALL irreps in 0.03 s ([14,7,0]) / 0.4–0.7 s ([28,14,0], depending on machine load).
   W4's exact slices (every HWBASIS ket and every BRANCH `cg`) satisfy E₁₂^{λ₁−λ₂+1}v = E₂₃^{λ₂+1}v = 0 exactly, lie in W_μ, and span M_μ (exact ranks = m in all 21 sectors; surd branches: residual of the expansion in the PRV basis ≤ 2e-63 at 60 digits). The λ₃ > 0 PRV/RS code was smoke-tested ([2,1,1], [3,2,1], [4,2,2], [5,3,1]: dims = RS; ĉ raises for λ₃ ≠ 0 as documented).
   Observed pattern (a = 1,2,3,4,7,14): Λ = [2a,a,0] has 2a(a+1)+1 irreps, of which 4(a+1−m) have multiplicity m for 1 ≤ m ≤ a and exactly one has m = a+1; copies = (2/3)a(a+1)(a+2)+(a+1). **For m = 1 the slice IS the PRV vector — the Q-route data only fix the scale 1/X_B and test consistency; the direction in M_μ is genuinely determined by Q-data only for m ≥ 2:** 316 of 344 copies at [14,7,0] (28 irreps have m = 1), 2199 of 2255 at [28,14,0] (56 irreps).
4. **The (m+1)-point dictionary on the 26 known branches.** L_B taken from N1's records (200-digit leading polynomials, both φ), m+1, m+2 and 2m+2 points, dps = 120, compared with W4's exact slice (both normalised on its first nonzero entry): all 26 branches × 2 φ = 52 records recovered; the multiplicity-2 copies (20 records) each match their own exact slice and are different vectors of the 2-dimensional M_μ.

   | rep | records | worst residual (m+1 pts) | worst slice error, m = 2 sectors (m+1 / m+2 / 2m+2 pts) | worst abs. scale max\|X_B·RAW − exact slice\|/max | worst held-out prediction of L_B at N1's 64-point grid | cond(A) (m = 2) |
   |---|---|---|---|---|---|---|
   | [2,1,0] | 12 | 1.1e-53 | 4.3e-65 / 4.0e-65 / 1.6e-63 (≥ 62 digits) | 8.4e-54 | 8.4e-54 | 1.5–3.8 |
   | [3,1,0] | 20 | 1.3e-91 | 5.7e-100 / 5.8e-100 / 7.8e-98 (≥ 97 digits) | 9.3e-92 | 9.3e-92 | 1.4–4.0 |
   | [3,2,0] | 20 | 2.1e-113 | 1.5e-119 / 1.5e-119 / 1.9e-117 (≥ 116 digits) | 2.3e-113 | 2.3e-113 | 1.5–4.5 |

   The accuracy is the accuracy of N1's L_B (adjoint: 12-point U2 fits, ≈ 1e-53…1e-66; [3,1,0]: ≈ 1e-92…1e-114; [3,2,0]: 1e-113…1e-121), not of the dictionary: the **identity-only test** (L generated from the exact slice through the identity, then recovered; the self-consistency test the plan allowed) returns the slice to the working precision (3e-120 at 120 digits) on all 52 records.
   `max|X_B·RAW − exact slice|` (X_B from N1's records) is the independent check that the **constant in the identity is exactly 1** through the m+1-point route. Negative control: feeding L_B of one irrep to the dictionary of ANOTHER irrep with the same (m,w) gives residual ≥ 0.35 (smallest of 88 wrong cases) against ≤ 2.1e-53 for the right irrep — the residual is a sharp consistency test. Rounding the values to float64 (what a machine-precision kernel would deliver) gives slice errors 1e-17…2e-16 at m = 2.
5. **Independent agreement with WP-R (`slice_engine.py`, `ref/cg-slices-*.txt`).** Same irrep sets; dim M_μ, w and WSPACE identical; **PRVBASIS identical to mine, exactly, in every sector** (13 + 25 + 41 + 113 sectors, 19 + 44 + 85 + 344 basis vectors; both agents use the RREF form, so the coordinates x_j are directly comparable); R's COPY slices (the V_red eigenvectors = the SoV copies) lie in M_μ to the 30-digit data accuracy (0 / 2e-30 / 6e-30 / 7e-29); the **dictionary applied to the identity-generated L of R's 344 copies at [14,7,0]** (shared points, m up to 8) returns R's slices to 1.5e-26 (worst residual 9e-31; the floor is R's 30-digit decimals × cond), at [4,2,0], [6,3,0], [8,4,0]: 7e-40, 6e-29, 1.4e-28.
6. **Timing (< 1 s per point everywhere).** Λ = [2a,a,0], dim (a+1)³, medians over 5 points, seconds (machine shared with the other agents, so ±2×):

   | a | dim | rep1 setup | tensor T | covector exact | covector mp-50 | ĉ | #irreps | PRV for all irreps |
   |---|---|---|---|---|---|---|---|---|
   | 1 | 8 | 0.000 | 0.000 | 0.0000 | 0.0001 | 1e-5 | 5 | 0.00 |
   | 3 | 64 | 0.001 | 0.001 | 0.0005 | 0.0005 | 1e-5 | 25 | 0.00 |
   | 5 | 216 | 0.004 | 0.003 | 0.0027 | 0.0024 | 2e-5 | 61 | 0.00 |
   | **7 ([14,7,0])** | **512** | **0.009** | **0.010** | **0.010** | **0.008** | **2e-5** | **113** | **0.02** |
   | 9 | 1000 | 0.020 | 0.044 | 0.045 | 0.035 | 4e-5 | 181 | 0.07 |
   | 11 | 1728 | 0.136 | 0.051 | 0.067 | 0.052 | 4e-5 | 265 | 0.18 |
   | 14 ([28,14,0]) | 3375 | 0.108 | 0.255 | 0.256 | 0.182 | 4e-5 | 421 | 0.52 |

   (The covector column is the FULL φ_u over all dim components; the dictionary only needs the w ≤ 15 components in W_μ, which is what `DictionaryProblem` evaluates.) A dictionary solve for one irrep (matrix at m+1 points + QR, all copies sharing it): numeric 2.3 ms mean / 5.7 ms max at [14,7,0], 11 ms / 41 ms at [28,14,0]; exact-rational path 3 ms / 18 ms and 48 ms / 0.8 s. Setup per Λ ≈ 0.02 s ([14,7,0]) / 0.4 s ([28,14,0]).
7. **Which points to use at large m (the one real pitfall).** cond(DA) over ALL irreps of multiplicity m (m ≥ 2), four point rules — `default` = first m+1 of `default_points`, `chosen` = per-irrep greedy (column-pivoted QR over a pool of 40), `shared` = prefix m+1 of ONE nested list `shared_points(lam, phi)` for the whole Λ, `default2` = first 2m+2 default points:

   | Λ | worst cond default | chosen | **shared prefix** | default 2m+2 |
   |---|---|---|---|---|
   | [8,4,0] (m ≤ 5) | 5.7e2 | 3.8 | 4.4 | 17 |
   | [14,7,0] (m ≤ 8) | 1.8e5 | 7.9 | **43** | 85 |
   | [28,14,0] (m ≤ 15) | 6.4e9 | 92 | **3.0e3** | 5.2e3 |

   (per m at [14,7,0], shared: 6.7, 20, 43, 40, 12, 14, 10 for m = 2…8; at [28,14,0] up to 3.0e3 at m = 11.) Caveat: the shared list comes from a float64 pivoted QR (heuristic): its worst cond at [28,14,0] varies between 1.2e3 and 3.6e5 with the pool size (40/48/60/100) — use `choose_points` per irrep if every digit matters.
   `shared_points` has the practical advantage that the kernel pieces G^{[d]}(ω̂_k) are needed at only the first m_max+1 = 9 points ([14,7,0]) / 16 ([28,14,0]).
8. **Noise.** Relative noise δ on every L_B value, direction error = ‖s−ce‖/‖s‖ (scale free), amplification = worst error/δ: real m = 2 branches (12 records, 4 trials, δ = 1e-10, 1e-20, 1e-40): ≤ 1.55 (m+1 default), ≤ 0.68 (m+1 chosen), ≤ 1.24 (2m+2) — no amplification. Synthetic at scale, δ = 1e-20, one irrep per m:
   m+1 default points amplify up to 3.9e7 (m = 15 at [28,14,0]; 2.9e3 at m = 8, [14,7,0]), chosen points ≤ 8.5 (m ≤ 15), 2m+2 default ≤ 33. (The componentwise error relative to the *first nonzero entry* can be up to ~100× larger than the direction error when that entry is small; compare vectors with `dir_err` or normalise on the largest entry.)
9. **Synthetic scale test** (L = ĉ⟨φ_u|s⟩ for a random exact s ∈ M_μ, m+1 default points, EVERY irrep; labelled synthetic): exact path recovers the exact PRV coordinates in **113/113 irreps of [14,7,0] and 421/421 of [28,14,0]** (also 13/13, 25/25, 41/41 at [4,2,0], [6,3,0], [8,4,0]; residual exactly 0, so the (m+1)×m matrix has full column rank for all of them);
   numeric path at 50 digits: worst direction error 2e-50 ([4,2,0]), 8e-50, 2e-49, 1.1e-47 ([14,7,0]), 6e-43 ([28,14,0]; cond 6.4e9), with the chosen points 6e-51 and 3e-50.

## 2. Usage for round 7

    import gauss_dictionary as gd
    pts = gd.shared_points(lam, phi)[:m+1]                       # m = gd.multiplicity(lam, mu); phi = the eps-regime direction of the L_B data
    vals = [L_B(w1, w1p) for (w1, w1p) in pts]                   # leading coefficient eps^N_B of sum G nu nu at the points (exact Fractions or mpf)
    x, residual, s = gd.slice_from_points(lam, mu, phi, pts, vals, dps=60)    # or gd.slices_from_points(lam, mu, phi, pts, [vals_copy1, ...])
    # s: first nonzero entry (my GT order) = 1, directly comparable with WP-R's COPY k (same normalisation); x: coordinates of C/X_B in gd.prv_basis(lam, mu) (x_j = RAW[free index j]);
    # residual ~ data precision x cond(A) when the identity holds; O(1) if mu or the branch is wrong.

Pair a Q-route branch with an R copy by the V_red eigenvalue = a₁ − η₁/2 (R's header). Points must avoid w1 = w1p and the φ_i (checked; ValueError).

## 3. Surprises / remarks

1. For m = 1 (28 of 113 irreps at [14,7,0]) the dictionary has nothing to decide: the slice is the PRV vector. All the content of the Q-route sits in the m ≥ 2 sectors; the real Q-data for m ≥ 2 so far exist only for m = 2 (10 branches) — m ≥ 3 is validated only on synthetic L and on R's V_red copies.
2. ĉ's only blind spot is the simple/multiple pole at ω̂₁′ = φ₂ (never at ω̂₁ = φ₂ or φ₃): x₁₃ = g′+a′c′ has (φ₂−ω̂₁) but not (φ₂−ω̂₁′) in its denominator, which is the whole asymmetry between the two points.
3. The "a = 14" of the plan is Λ = [28,14,0] (dim 3375), not [14,14,0]; [14,7,0] is a = 7 (dim 512). Both are covered; the Q-route needs only the PRV bases (≤ 0.7 s for all 421 irreps), ĉ and φ_u (0.26 s) there.
4. The naive first m+1 generic points lose up to 9 digits at m = 15 — the conditioning is a property of the point set, not of the problem (shared/chosen points: ≤ 3e3 / ≤ 92).
5. N1's adjoint L_B has only ~1e-53…1e-66 accuracy (12-point fits); the adjoint rows of the table measure the data, not the dictionary.
6. The Λ-only tensor T has 2708 entries at [14,7,0] and 33530 at [28,14,0]; φ_u for one point is then a sum of ≤ 15 monomials per component.

## 4. Not done

* No real (Q-route) data for m ≥ 3: S/Q0 produce [4,2,0]/[6,3,0] branches in this round; the round-7 test is `slices_from_points` on them against R's copies. No proof that the (m+1)×m matrix has full column rank at generic points (it holds for all 421 irreps of [28,14,0] with the default points: the exact solve never hit a singular system).
* ĉ and the pole analysis are for λ₃ = 0 only (covector, PRV, Racah–Speiser are general in λ₃; smoke-tested). Absolute normalisation/sign of the slice (the G-norm of the full vector, Q2's X_B) is not part of D: RAW = C/X_B keeps the scale.
* `shared_points` is a heuristic ordering (see the caveat in item 7); no optimal (Fekete-type) point sets were attempted.
* Timings were taken on a machine shared with the other round-6 agents (±2×).
