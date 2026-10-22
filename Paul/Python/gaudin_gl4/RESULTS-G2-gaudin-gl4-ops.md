# RESULTS-G2 (agent G2, 2026-10-22): the two-site gl(4) Gaudin missing-label operators — cdet(∂ − L(u)) for gl(N), its Laurent coefficients, commutativity, the canonical triple, simple joint spectrum

Files (all new; nothing existing was edited; all in `Paul/Python/gaudin_gl4/` unless stated):
`gaudin_gl4_ops.py` (all code; CLI `python -X utf8 gaudin_gl4_ops.py sym | scaffold | gl3 | comm | triple | struct | terms | big | all`),
`test_gaudin_gl4_ops.py` (11 exact regression tests, ≈ 5 s, `python -X utf8 test_gaudin_gl4_ops.py`: **11 tests, 0 failed**),
`gaudin_gl4_ops_output.txt` (the printed output of one run of `all`, 301 s, 1437 lines), this file, and
`../../Mathematica/Gaudin/ref/gl4-gaudin-ops.txt` (the exact term lists, §7; loader `load_terms` in the module).
Python 3.13.13, numpy 2.5.2, scipy 1.18.1, sympy 1.14 (only for exact discriminants/Sturm counts and one representation-free check), mpmath 1.3 (40 digits); i7-12700KF, Windows 11; ≤ 2 heavy processes at any time.
**Everything quantitative below is exact over ℚ** (Fractions; exact discriminants and Sturm counts for the simplicity claims); mpmath eigenvalues are only used to print numbers and measure gaps.

Representations: a **private scaffold** (`Rep`, from `gl3_reps.gt_patterns/_raise/_lower/_weight`, E_ij = [E_{i,j−1},E_{j−1,j}] etc., exactly the brief's recipe) was used for everything. G1's `glN_reps.GLNRep` appeared meanwhile and its own tests pass (10 tests, 0 failed, run by me), so I **cross-checked** the scaffold against it read-only (`scaffold` stage, test 11): patterns, weights, HW index, all 16 generators and the metric are identical entry by entry on [1,0,0,0], [1,1,0,0], [2,1,0,0], [2,1,1,0], [2,2,1,0], [3,2,0,0], [3,2,1,0], [4,2,1,0]; Racah–Speiser identical on 5 products. For N = 3 the scaffold equals `slice_engine.Rep.of(lam)` entry by entry (E, metric) on [2,1,0], [3,1,0], [3,2,0], [4,2,0]; `slice_engine.VTERMS/eta1/racah_speiser` are imported read-only.

## 0. Headline results

1. **Convention test passed.** The symbolic engine (column-ordered cdet(∂−L), derivative moved to the right) gives for N = 3 exactly S₁, S₂, S₃ of eq. (4) of the notes as noncommutative polynomials (3, 8, 13 words). For N = 4: S₁..S₄ have 4, 15, 42, 73 words (printed in the output). The **row-ordered variant fails** (not even gl(N)-invariant with L_ij = E_ij/(u−θ)): 37/45 (N=3) and 136/153 resp. 137/153 (N=4) pairs do not commute; the "derivative-to-the-left" re-normal-ordering is the same commutative family. So the binding convention is the right one and no fallback was needed.
2. **Every Laurent coefficient is homogeneous**: S_{k,j} = t^{j−k}·Ŝ_{k,j} with Ŝ_{k,j} t-independent (exact). Hence the canonical triple **V₃ = t S_{3,2}, V₄ᵃ = t S_{4,3}, V₄ᵇ = t² S_{4,2} is exactly t-independent**; term counts 88 / 152 / 324 (`S_{4,1}`: 412). All residue relations (S_k = O(u^{−k}) at ∞) hold identically in U(gl_N)⊗U(gl_N), which checks stage B and the θ₂ rule (site swap, t → −t) without any representation.
3. **Commutativity, N = 4, exact (Fractions, whole space):** all 18 operators S_{k,j}, k = 2,3,4, j = 1..k, at θ₁ and θ₂ commute pairwise and with all 16 E^tot_ij: **0 violations** on [1,0,0,0]⊗[1,1,0,0] (24), [1,1,0,0]⊗[2,1,0,0] (120), [2,1,0,0]² (400), [2,1,1,0]² (225), at t = 31/35 and −3/7; float64 matrix-free cross-checks (also [3,2,0,0]², [3,2,1,0]², dim 3600, 4096) at 1e-15.
4. **Double pole:** S_{2,2} = λ₂(λ₁+1) + λ₃(λ₁+λ₂+2) + λ₄(λ₁+λ₂+λ₃+3) = e₂(λ) + Σ_k (k−1)λ_k (13 reps verified; central scalar); in general S_{m,m} = (−1)^m a_{N−m} with Π_c(θ − y_c) = Σ_k a_k θ(θ−1)…(θ−k+1), y_c = λ_c + N − c (closed forms for S_{3,3}, S_{4,4} verified).
5. **The triple on every multiplicity space** (14 products, up to dim 531 441; 180 sectors with m ≥ 2, Σm = 678 copies, m up to 12): HW spaces exact; V₃, V₄ᵃ, V₄ᵇ commute exactly, are self-adjoint for the invariant form (exact), spectra real; **the joint spectrum is simple on every sector** (exact certificate: a generic rational combination has nonzero discriminant and Sturm count m). No degenerate joint spectrum anywhere ⇒ nothing to stop.
6. **V₃ alone is NOT enough for gl(4)**: it is degenerate (exact discriminant = 0) on some sectors of all the larger λ = λ′ products ([3,2,1,0]², [4,2,1,0]², [4,3,1,0]², [5,3,1,0]², [5,3,2,0]², [6,4,2,0]²: 2, 2, 6, 6, 8, 19 sectors), e.g. [3,2,1,0]²: (5,4,2,1) (m = 4, 3 distinct values) and (4,4,2,2) (m = 2, V₃ is a scalar there); it is simple on every sector of the λ ≠ λ′ products and of the small λ = λ′ ones. The pair **{V₃, V₄ᵇ} separates every sector of every product tested**; {V₃, V₄ᵃ} fails exactly once, on the m = 12 sector (9,7,5,3) of [6,4,2,0]² (V₄ᵇ resolves it); no single quartic and not {V₄ᵃ, V₄ᵇ} suffices. The degeneracy of V₃ is explained in part by the site-exchange symmetry (§6e).
7. **Structure (exact identities in U(gl₄)⊗U(gl₄)):** V₃ = Σ_I D(2,1)_I + (n₁n₂ − Ω) = T(1,1,2) + ½n₁²n₂ + ½n₁n₂ − n₁Ω − ½n₂C₂¹ + Ω = −T(1,2,2) + [Casimirs]; V₄ᵃ = D(3,1) + V₃ − η/3, η = n₁n₂ − Ω = t S_{2,1}; V₄ᵇ = D(2,2) − D(3,1) + ½(ΣD(1,2) − ΣD(2,1)) − η/6; trace form in §6. For N = 3 the same engine gives V₃ = V + η/2 exactly.
8. **The triple generates the whole family:** every one of the 18 Laurent coefficients is an exact integer combination of V₃, V₄ᵃ, V₄ᵇ plus polynomials in Casimir-type scalars (table in §6c); in particular all other Laurent coefficients are diagonal in the joint eigenbasis (checked also sector by sector).
9. **Naive traces:** T(1,1,2,2) and T(1,2,1,2) do **not** commute with V₃ on HW spaces (20 of 22 sectors) but T(1,1,2,2) + ½T(1,2,1,2) does; T(1,1,2), T(1,2,2), T(1,1,1,2), T(1,2,2,2) commute with the whole triple.
10. **Term lists written** (`ref/gl4-gaudin-ops.txt`, 12 operators, 1861 lines, loader + re-verification on two sectors).

## 1. Conventions and the convention test (G2.1a)

Binding conventions of the brief §1 (generalised verbatim from the gl(3) notes §1.1): [E_ij,E_kl] = δ_jk E_il − δ_li E_kj, two sites V_λ (site 1) ⊗ V_λ′ (site 2), e^s_ij = E_ij on site s, e¹ and e² commute, **L(u)_ij = e¹_ij/(u−θ₁) + e²_ij/(u−θ₂)**, t = θ₁₂ = θ₁ − θ₂, x = u − θ₁ (so u − θ₂ = x + t),
cdet(∂_u − L) = Σ_σ sgn σ M_{σ(1)1}⋯M_{σ(N)N} = ∂^N − S₁∂^{N−1} + S₂∂^{N−2} − … ; Laurent data S_k(u) = Σ_j S_{k,j} x^{−j} + regular; θ₂ data = site swap + t → −t. GT basis as `gl3_reps.gt_patterns` (HW last), plain rational GT matrices on kets as columns, two-site index ν·dim₂ + q.

Stage A (`cdet_S`) is a noncommutative differential-operator algebra (words of factors (i,j,d) = d-th u-derivative of l_ij; (A∂^a)(B∂^b) = Σ_j C(a,j) A B^{(j)} ∂^{a+b−j}). Stage B (`laurent`) substitutes l^{(d)}_ij = (−1)^d d![e¹_ij x^{−d−1} + e²_ij (x+t)^{−d−1}] and expands (x+t)^{−B} = Σ_n C(−B,n) t^{−B−n} xⁿ.

**N = 3 convention test:** engine = eq. (4) exactly for S₁ (3 words), S₂ (8), S₃ (13) — both printed in the output, asserted equal. (S₃ = l″₃₃ − l₁₁l′₃₃ − l₂₂l′₃₃ − l′₂₂l₃₃ + l₃₁l′₁₃ + l₃₂l′₂₃ + l′₃₂l₂₃ + l₁₁l₂₂l₃₃ − …, order of factors as written.) For N = 4 the engine's S₂ = Σ_{i<j}(l_ii l_jj − l_ji l_ij) − l′₂₂ − 2l′₃₃ − 3l′₄₄ (the pattern −Σ(k−1)l′_kk), S₃, S₄ (73 words, up to l⁽³⁾₄₄) are printed in the output.

Variants (information, exact, `comm` stage): column-ordered (binding) 0 non-commuting pairs and 0 [S,E^tot] violations (N=3 [2,1,0]²: 28 pairs; N=4 [1,1,0,0]⊗[2,1,0,0]: 66; [2,1,1,0]²: 120); row-ordered Σ sgn M_{1σ(1)}⋯M_{Nσ(N)}: **37/45, 136/153, 137/153 pairs non-commuting, 48/190/192 [S,E^tot] violations**; column-ordered with coefficients read off in the derivative-left normal form: commutes (a re-basing of the same commutative algebra).

## 2. Laurent data at θ₁: term counts, homogeneity, representation-free checks

Merged terms (w₁,w₂) (site-1 factors first, each site's internal order = cdet order; e¹, e² commute) per S_{k,j}, N = 4 (N = 3 in brackets); j = 0 is the regular part x⁰ at θ₁:

| | j = 4 | j = 3 | j = 2 | j = 1 | j = 0 |
|---|---|---|---|---|---|
| S₂ | | | 15 (8) | 24 (12) | 39 (20) |
| S₃ | | 36 (11) | 88 (24) | 156 (40) | 198 (55) |
| S₄ | 49 | 152 | 324 | 412 | 493 |

(i.e. S_{2,2}: 15, S_{2,1}: 24; S_{3,3}: 36, S_{3,2}: 88, S_{3,1}: 156; S_{4,4}: 49, S_{4,3}: 152, S_{4,2}: 324, S_{4,1}: 412; the raw (word, site-assignment) counts before merging are 54, 260, 722 for S₂, S₃, S₄.) The merged lists are minimal: PBW normal-ordering does not shorten them (same counts).
* **Homogeneity (asserted for all k, j, N = 3, 4):** every term of S_{k,j} carries exactly t^{j−k}. (Reason: S_k(u) has total weight k, since Σ(d+1) over the factors of every word is k.) Thus S_{3,3}, S_{4,4} are the top poles with t⁰, V₃ = t S_{3,2} etc. have t⁰.
* **Classical check of stage B** (commuting random integer matrices instead of e^s_ij, sympy series of the rational functions): N = 3: 9 coefficients × 2 seeds, N = 4: 14 × 2 seeds, 0 mismatches.
* **Residue relations** Σ_s Σ_{i≤j} C(j−1,i−1) θ_s^{j−i} S^{(s)}_{k,i} = 0 for j = 1..k−1 (S_k = O(u^{−k}) at ∞), with the θ₂ operators defined by site swap and t → −t: **hold identically** in U(gl_N)⊗U(gl_N) (PBW normal form, 0 residual terms) for N = 3 (3 relations) and N = 4 (6 relations). E.g. S^{(2)}_{k,1} = −S^{(1)}_{k,1} for all k ≥ 2.
* N = 3 operator identities (PBW, exact, whole algebra): t·S_{2,1} = n₁n₂ − Ω (this is H₁ = −(1/t)Σ E¹_ik E²_ki + |Λ|²/t with |Λ|² = n₁n₂), and **t·S_{3,2} = V + ½(n₁n₂ − Ω)** with V = ½Σ sgn σ sgn τ e¹e¹e² (36 terms), i.e. the notes' θ₁₂a₁ = V + η₁/2 as an identity in U(gl₃)⊗U(gl₃) (n_s = Σ_i e^s_ii, Ω = Σ e¹_ij e²_ji, η₁ = |Λ|² − Ω).

## 3. G2.1: gl(3) numerical facts (exact arithmetic)

[2,1,0]² (dim 64) and [3,1,0]² (dim 225), t = 31/35 and −3/7: all S_{2,j}, S_{3,j} at θ₁ and θ₂ (10 operators; 8 nonzero) commute pairwise (45 pairs) and with the 9 E^tot_ij (90 tests): **0 violations**.

| | [2,1,0]² | [3,1,0]² |
|---|---|---|
| S_{2,2} (θ₁, θ₂) | 3·1 = λ₂(λ₁+1) | 4·1 |
| H₁ = S_{2,1} vs −(1/t)ΣE¹_ik E²_ki + |Λ|²/t | 0 differing entries (whole space) | 0 |
| S_{3,3} (the triple pole) | **0 as an operator** on V_λ⊗V_λ (central scalar, vanishes identically, not only on HW vectors) | 0 |
| t·S_{3,2} = V + η₁/2 on HW vectors | 5 irreps (Racah–Speiser = `slice_engine.racah_speiser`), 6 exact HW vectors, 0 violations | 8 irreps, 10 vectors, 0 |
| same, whole-space operator identity t S_{3,2} = V + ½(n₁n₂ − Ω) | 0 differing entries | 0 |

The triple pole is the central scalar S_{3,3} = (λ₁+2)(λ₂+1)λ₃ (zero iff λ₃ = 0; checked for [2,1,1], [3,2,1], [4,2,2], [2,2,2]: 8, 15, 36, 24), S_{2,2} = λ₂(λ₁+1) + λ₃(λ₁+λ₂+2).

## 4. G2.2: commutativity for N = 4

Exact (Fractions, whole space), t = 31/35 and t = −3/7, all pairs of the nonzero operators among the 18 and all 16 E^tot_ij:

| product | dim | nonzero operators | pairs | violations | E^tot tests | violations | float64 matrix-free (2 random vectors) |
|---|---|---|---|---|---|---|---|
| [1,0,0,0]⊗[1,1,0,0] | 24 | 6 | 15 | **0** | 96 | **0** | 6e-16 / 4e-16 |
| [1,1,0,0]⊗[2,1,0,0] | 120 | 12 | 66 | **0** | 192 | **0** | 2.6e-15 / 8e-16 |
| [2,1,0,0]² | 400 | 12 | 66 | **0** | 192 | **0** | 2.3e-15 / 7e-16 |
| [2,1,1,0]² | 225 | 16 | 120 | **0** | 256 | **0** | 8e-16 / 6e-16 |
| [3,2,0,0]² | 3600 | 12 | float only | | | | 1.2e-15 / 9e-16 |
| [3,2,1,0]² | 4096 | 16 | float only | | | | 1.1e-15 / 9e-16 |

(Operators that vanish identically on a product are listed in the output: e.g. on [2,1,0,0]² and [1,1,0,0]⊗[2,1,0,0] the six S_{3,3}, S_{4,3}, S_{4,4} (θ₁ and θ₂) — S_{3,3}, S_{4,4} are the central scalars, zero for λ₃ = 0 resp. λ₄ = 0, and V₄ᵃ itself vanishes there —, on [2,1,1,0]² only S_{4,4}. The float residuals are relative to ‖S_aX‖‖S_bX‖/‖X‖.) The whole-space exact arithmetic is cheap because the operators conserve the total weight (≈ 6–10 nonzeros per column).

**Double-pole scalar for gl(4)** (S_{2,2} = central scalar of site 1, measured on V_λ as an operator, vs formula vs Newton closed form; all equal on the 13 reps [1,0,0,0], [1,1,0,0], [2,1,0,0], [2,1,1,0], [2,2,0,0], [3,1,0,0], [3,2,1,0], [2,2,2,0], [3,3,0,0], [4,2,1,0], [2,2,1,1], [3,1,1,1], [4,3,2,1]):

  **S_{2,2} = λ₂(λ₁+1) + λ₃(λ₁+λ₂+2) + λ₄(λ₁+λ₂+λ₃+3) = e₂(λ) + Σ_k (k−1)λ_k**  (0, 2, 3, 8, 6, 4, 15, 18, 12, 18, 20, 18, 45 on those reps).

Derivation (and verification of S_{3,3}, S_{4,4}): cdet(∂ − E/x) = x^{−N}Π_c(θ − y_c) with θ = x∂_x, y_c = λ_c + N − c on the HW vector, hence S_{m,m} = (−1)^m a_{N−m}, a_k = Δ^k f(0)/k!, f(θ) = Π(θ − y_c). Examples: S_{3,3} = (λ₁+2)(λ₂+1)λ₃ for λ₄ = 0 (e.g. 8, 15, 18 on [2,1,1,0], [3,2,1,0], [4,2,1,0]; 0 if λ₃ = 0), S_{4,4} = (λ₁+3)(λ₂+2)(λ₃+1)λ₄ (0 for λ₄ = 0; 40, 36, 105 on [2,2,1,1], [3,1,1,1], [4,3,2,1]).
**Higher poles:** S₃ has the triple pole S_{3,3} (central scalar); S₄ has S_{4,4} (central scalar) and **S_{4,3} = V₄ᵃ/t, a genuine two-site operator** (every term has a site-2 factor); no pole of order > k occurs.

## 5. G2.3: the canonical triple on the multiplicity spaces

Method (all exact): HW space of weight μ in V_λ⊗V_λ′ = null space over ℚ (sparse RREF, canonical basis: 1 at a free pair of the weight-μ subspace, 0 at the other free pairs) of the total simple raising operators in the plain GT basis; **dim = Racah–Speiser multiplicity on every sector** (asserted); the operators are applied through the Kronecker-sum form O = Σ_g P_g ⊗ Q_g (grouping on the site with fewer distinct words), restricted by O b_j = Σ_i R_ij b_i verified entry by entry (faithful: 0 violations everywhere); self-adjointness for the invariant form G⊗G (diagonal GT metric): Rᵀ Gram = Gram R exactly; exact characteristic polynomials (Faddeev–LeVerrier), exact discriminants and Sturm counts (sympy) as simplicity certificates; mpmath (40 digits) eigenvalues/eigenvectors of the generic combination M = V₃ + (7/19)V₄ᵃ + (11/13)V₄ᵇ for the joint eigenvalue triples. **t-independence:** the term list evaluated at t = 31/35 and −3/7 and multiplied by t reproduces the t-free operator exactly (checked for V₃ on every product; it is exact by homogeneity).

**Brief typo:** dim[2,2,1,0] = 20 (dual of [2,1,0,0]), not 60; [2,2,1,0]⊗[2,1,1,0] has dim 300. To honour the "60·15 = 900" intent I also ran [3,2,0,0]⊗[2,1,1,0] (60·15 = 900).

Multiplicity tables (μ : m; dim V_μ; n_W; Σ m dim V_μ = dim·dim′ true for all):
* [2,1,0,0]² (400; 7 irreps, Σm = 8): 4200:1, 4110:1, 3300:1, **3210:2**, 3111:1, 2220:1, 2211:1.
* [2,1,1,0]² (225; 6, 7): 4220, 4211, 3320, 3311 (all 1), **3221:2**, 2222:1.
* **[3,2,1,0]² (4096; 16 irreps, Σm = 30, max m = 4)**: 6420:1, 6411:1, 6330:1, **6321:2**, 6222:1, 5520:1, 5511:1, **5430:2**, **5421:4**, **5331:3**, **5322:3**, 4440:1, **4431:3**, **4422:2**, **4332:3**, 3333:1.
* [2,2,1,0]⊗[2,1,1,0] (300; 6, 7): 4320:1, 4311:1, 4221:1, 3330:1, **3321:2**, 3222:1.   [3,2,0,0]⊗[2,1,1,0] (900; 7, 8): 5310:1, 5211:1, 4410:1, 4320:1, **4311:2**, 4221:1, 3321:1.
* [3,2,0,0]² (3600; 15, 18): m = 2 at 5410, 5320, 4321, else 1.   [2,2,1,0]² (400; 7, 8): m = 2 at 4321 only (dual to [2,1,0,0]²).
* Larger products (extras, because gl(4) multiplicities with m ≥ 3 are scarce in the required list): [4,2,1,0]² (19600; 27 irreps, Σm = 56, max m = 5), [3,2,1,0]⊗[4,2,1,0] (8960; 20, 38, max 4), [4,2,2,0]² (7056; 19, 26, max 3), [4,3,1,0]² (30625; 31, 69, max 5), **[5,3,1,0]² (129 600; 49 irreps, max m = 7)**, **[5,3,2,0]² (90 000; 45, max m = 8)**, **[6,4,2,0]² (531 441; 76 irreps, Σm over the 56 sectors with m ≥ 2 = 277, max m = 12, n_W up to 4050)**. Full tables are in the output file.

**Per-product summary** (every sector with m ≥ 2; "comm" = triple pairwise commuting on the HW space (exact), "self-adj" = exact; V₃ simple = exact discriminant ≠ 0 and Sturm count m; joint simple = same for the generic combination; gaps = smallest distance between distinct eigenvalues, resp. between joint eigenvalue triples in max-norm):

| product | dim | # sectors m ≥ 2 | Σ m | max m | comm | self-adj | V₃ simple | **joint simple** | min V₃ gap (simple sectors) | min joint distance | V₃-degenerate sectors (μ, m, #distinct V₃ values) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| [2,1,0,0]² | 400 | 1 | 2 | 2 | yes | yes | yes | **yes** | 6.7082 | 6.7082 | none |
| [2,1,1,0]² | 225 | 1 | 2 | 2 | yes | yes | yes | **yes** | 13.856 | 27.713 | none |
| **[3,2,1,0]²** | 4096 | 8 | 22 | 4 | yes | yes | **NO** | **yes** | 6.7082 | 6.7082 | (5,4,2,1) m=4 (3); (4,4,2,2) m=2 (1) |
| [2,2,1,0]⊗[2,1,1,0] | 300 | 1 | 2 | 2 | yes | yes | yes | **yes** | 8.544 | 51.264 | none |
| [3,2,0,0]⊗[2,1,1,0] | 900 | 1 | 2 | 2 | yes | yes | yes | **yes** | 10.247 | 20.494 | none |
| [3,2,0,0]² | 3600 | 3 | 6 | 2 | yes | yes | yes | **yes** | 6.7082 | 6.7082 | none |
| [2,2,1,0]² | 400 | 1 | 2 | 2 | yes | yes | yes | **yes** | 6.7082 | 26.833 | none |
| [4,2,1,0]² | 19600 | 14 | 43 | 5 | yes | yes | NO | **yes** | 3.1749 | 10.247 | (6,4,2,2) m=4 (3); (5,5,2,2) m=2 (1) |
| [3,2,1,0]⊗[4,2,1,0] | 8960 | 10 | 28 | 4 | yes | yes | yes | **yes** | 1.7130 | 8.544 | none |
| [4,2,2,0]² | 7056 | 6 | 13 | 3 | yes | yes | yes | **yes** | 19.596 | 58.788 | none |
| [4,3,1,0]² | 30625 | 18 | 56 | 5 | yes | yes | NO | **yes** | 5.8977 | 5.8977 | 6 of 18: (7,6,2,1) m=4 (3); (7,5,3,1) m=5 (3); (6,6,2,2) m=2 (1); (6,5,3,2) m=5 (3); (6,4,4,2) m=4 (3); (5,5,3,3) m=2 (1) |
| **[5,3,1,0]²** | 129600 | 32 | 119 | 7 | yes | yes | NO | **yes** | 1.7571 | 9.2036 | 6 of 32: (8,6,2,2) m=4; (8,4,4,2) m=4; (7,7,2,2) m=2; (7,5,3,3) m=4; (6,6,4,2) m=4; (6,6,3,3) m=2 |
| **[5,3,2,0]²** | 90000 | 28 | 104 | 8 | yes | yes | NO | **yes** | 4.3301 | 10.247 | 8 of 28: (9,6,4,1) m=4; (8,7,3,2) m=4; (8,6,4,2) m=8; (8,6,3,3) m=4; (7,7,4,2) m=4; (7,7,3,3) m=2; (7,6,4,3) m=6; (6,6,4,4) m=2 |
| **[6,4,2,0]²** | 531441 | 56 | 277 | **12** | yes | yes | NO | **yes** | 2.6867 | 10.247 | 19 of 56 (e.g. (10,8,4,2) m=10; (10,7,5,2) m=11; (9,7,5,3) m=12; (9,8,4,3) m=8; (8,7,5,4) m=7; full list in the output) |

**Which subsets of the triple separate every sector** (exact: a generic rational combination of the subset has nonzero discriminant on every sector; "NO" = discriminant 0 for two independent generic combinations on some sector):

| product | {V₃} | {V₄ᵃ} | {V₄ᵇ} | {V₃,V₄ᵃ} | {V₃,V₄ᵇ} | {V₄ᵃ,V₄ᵇ} | {V₃,V₄ᵃ,V₄ᵇ} |
|---|---|---|---|---|---|---|---|
| [2,1,0,0]² | yes | NO | NO | yes | yes | NO | yes |
| [2,1,1,0]² | yes | yes | yes | yes | yes | yes | yes |
| [3,2,1,0]² | **NO** | NO | NO | yes | yes | NO | yes |
| [2,2,1,0]⊗[2,1,1,0] | yes | yes | yes | yes | yes | yes | yes |
| [3,2,0,0]⊗[2,1,1,0] | yes | NO | yes | yes | yes | yes | yes |
| [3,2,0,0]² | yes | NO | NO | yes | yes | NO | yes |
| [2,2,1,0]² | yes | yes | yes | yes | yes | yes | yes |
| [4,2,1,0]² | NO | NO | NO | yes | yes | NO | yes |
| [3,2,1,0]⊗[4,2,1,0] | yes | NO | NO | yes | yes | NO | yes |
| [4,2,2,0]² | yes | yes | yes | yes | yes | yes | yes |
| [4,3,1,0]² | NO | NO | NO | yes | yes | NO | yes |
| [5,3,1,0]² | NO | NO | NO | yes | yes | NO (5 sectors, m = 2,2,2,3,2, all with μ₄ = 0) | yes |
| [5,3,2,0]² | NO | NO | NO | yes | yes | NO ((9,6,5,0), (8,7,5,0)) | yes |
| [6,4,2,0]² | NO (19 sectors) | NO (8) | NO (6) | **NO: fails only on (9,7,5,3), m = 12** | yes | NO (6 sectors, μ₄ = 0) | yes |

**The headline table [3,2,1,0]²** (V_ρ ⊗ V_ρ, dim 4096; all 8 sectors with m ≥ 2; joint eigenvalue triples (V₃, V₄ᵃ, V₄ᵇ), 18-digit values are in the output):

| μ | m | n_W | V₃ simple? | joint simple? | min V₃ gap | min joint distance | eigenvalue triples (rounded) |
|---|---|---|---|---|---|---|---|
| (6,3,2,1) | 2 | 6 | yes | yes | 6.708 | 33.54 | (55.854, 54.271, 5.729), (49.146, 20.729, 39.271) |
| (5,4,3,0) | 2 | 6 | yes | yes | 6.708 | 6.708 | (34.146, 0, 0), (40.854, 0, 0) |
| (5,4,2,1) | 4 | 18 | **NO** (57 twice) | yes | 0 | 15.716 | (67.392, 70.981, −16.981), (46.608, 19.019, 34.981), **(57, 50.697, 21.303), (57, 21.303, 50.697)** |
| (5,3,3,1) | 3 | 28 | yes | yes | 12.649 | 28.377 | (77.649, 81.623, −21.623), (52.351, 18.377, 41.623), (65, 30, 70) |
| (5,3,2,2) | 3 | 44 | yes | yes | 8.718 | 26.153 | (68.282, 55.923, 46.077), (85.718, 82.077, 19.923), (77, 50, 90) |
| (4,4,3,1) | 3 | 44 | yes | yes | 8.718 | 26.488 | (77.718, 79.512, −9.512), (60.282, 18.488, 51.512), (69, 30, 78) |
| (4,4,2,2) | 2 | 68 | **NO** (V₃ = 81·1) | yes | 0 | 82.486 | (81, 76.748, 21.757), (81, 49.252, 104.243) |
| (4,3,3,2) | 3 | 102 | yes | yes | 8.718 | 35.206 | (89, 84, 28), (97.718, 86.794, 63.206), (80.282, 43.206, 106.794) |

(μ with m = 1 are also tabulated in the output: e.g. the top irrep (6,4,2,0): (V₃, V₄ᵃ, V₄ᵇ) = (25, 0, 0); the singlet (3,3,3,3): (105, 90, 90).)

**Reality, t-independence, diagonal:** (iv) all spectra are real (self-adjointness is exact for every operator on every sector, and the exact Sturm count equals m for the joint certificate); t-independent exactly (homogeneity; also checked at two t values); (v) the 15 other Laurent coefficients (θ₁ and θ₂, S₂..S₄, 18 operators in all, built at t = 31/35) commute **exactly** with the generic combination on every sector (all products except the very largest [6,4,2,0]², where only the triple was analysed) and their worst off-diagonal entry in the joint eigenbasis is ≤ 1e-36 (numerical zero at 40 digits). Stronger: they are exact combinations of the triple and scalars (§6c).

## 6. G2.4: structure of the triple

Notation: n_s = Σ_i e^s_ii; C₂ˢ = Σ e^s_ij e^s_ji; C₃ˢ = Σ e^s_ji e^s_kj e^s_ik (central on site s); Ω = Σ e¹_ij e²_ji; η := n₁n₂ − Ω (= t·S_{2,1}, = |Λ|² − Ω for equal reps); polarised traces (CPV) T(s₁…s_r) = Σ e^{s₁}_{i₂i₁} e^{s₂}_{i₃i₂}⋯e^{s_r}_{i₁i_r}, so T(1,1,2) = Σ e¹_ji e¹_kj e²_ik, T(1,2,2) = Σ e¹_ji e²_kj e²_ik; antisymmetrised forms D(a,b)_I = (1/a!b!) Σ_{σ,τ∈S(I)} sgn σ sgn τ e¹_{σ₁τ₁}…e¹_{σ_aτ_a} e²_{σ_{a+1}τ_{a+1}}…e²_{σ_{a+b}τ_{a+b}} (first factor acts last), D(a,b) = D(a,b)_{1234}, ΣD(2,1) := Σ_{|I|=3} D(2,1)_I. All identities below are **exact in the PBW normal form of U(gl₄)⊗U(gl₄)** (rational coefficients, solved by exact linear algebra over libraries of invariants whose PBW images are linearly independent; consistency = no inconsistent row) and are verified on 81 sectors (all irreps with n_W ≤ 260 of 7 products, m = 1 and m ≥ 2) with exact matrices.

**(a) cubic.**
* **V₃ = Σ_I D(2,1)_I + (n₁n₂ − Ω)** (general N: V₃ = Σ_I D(2,1)_I + ((N−2)/2)·η, exactly solved for N = 3, 4, 5, 6 with c_N = ½, 1, 3/2, 2; terms of V₃: 24, 88, 210, 408), where Σ_I D(2,1)_I = ½ Σ_{σ,τ∈S₄} sgn σ sgn τ e¹_{σ₁τ₁} e¹_{σ₂τ₂} e²_{σ₃τ₃} δ_{σ₄τ₄} (the "δ_{σ₄τ₄}-type" form of the brief; equality with the sum over 3-subsets of the gl(3) vertex operators checked as polynomials, 144 terms). So the cdet **does** produce the gl(3)-type antisymmetrised vertex operator summed over the four 3-subsets, plus η = t·S_{2,1}; the correction coefficient is (N−2)/2 (½ for N = 3, 1 for N = 4).
* **V₃ = T(1,1,2) + ½n₁²n₂ + ½n₁n₂ − n₁Ω − ½n₂C₂¹ + Ω**  (a = 1, b = 0). The remainder is a polynomial in the scalars n₁, n₂, C₂¹, Ω: **a multiple of the identity on every M_μ**: checked on all 81 sectors (29 with m ≥ 2) against the group-theory value Ω = (C₂(μ) − C₂(λ) − C₂(λ′))/2 — e.g. the value is 9 on (6,3,2,1) and (5,4,3,0), 24 on (5,4,2,1), 54 on (4,4,2,2) of [3,2,1,0]² (0 violations; value = prediction on all).
* T(1,1,2) + T(1,2,2) = ⅓(C₃^{tot} − C₃¹ − C₃²) + ⅔n₁n₂ − (8/3)Ω exactly (C₃^{tot} := the same trace form of the total generators e¹+e², central; 75 sectors checked with C₃(μ) evaluated on V_μ). Hence **V₃ = −T(1,2,2) + ⅓C₃^{tot} − ⅓C₃¹ − ⅓C₃² + ½n₁²n₂ + (7/6)n₁n₂ − n₁Ω − ½n₂C₂¹ − (5/3)Ω** (a = 0, b = −1): on a sector V₃ ≡ T(1,1,2) ≡ −T(1,2,2) modulo scalars; only a − b = 1 matters. For N = 3: V₃ = V + ½η.

**(b) quartic pair.**
* **V₄ᵃ = D(3,1) + ΣD(2,1) + ⅔η = D(3,1) + V₃ − η/3.**
* **V₄ᵇ = D(2,2) − D(3,1) + ½ΣD(1,2) − ½ΣD(2,1) − η/6** (classically P₂,₂ − P₃,₁; unique in the span of D(2,2), D(3,1), ΣD(1,2), ΣD(2,1), n₁n₂, Ω: rank 6, consistent; ΣD(1,2) is the site swap of ΣD(2,1)).
* Trace basis: **V₄ᵃ = −T(1,1,1,2) + (n₁ − 3)T(1,1,2) + [scalar polynomial in n₁, n₂, C₂¹, C₃¹, Ω]** (scalar part: ⅙n₁³n₂ − ⅙n₁²n₂ − ½n₁²Ω + n₁n₂ − ½n₁n₂C₂¹ + (3/2)n₁Ω + ⅔n₂C₂¹ + ⅓n₂C₃¹ + ½C₂¹Ω − 2Ω);
  **V₄ᵇ = T(1,1,1,2) − T(1,1,2,2) − ½T(1,2,1,2) + (n₂ − n₁ + 3/2)T(1,1,2) + (n₁ − 3/2)T(1,2,2) + [scalar polynomial]** (printed in full in the output). Both verified on all 81 sectors: remainder = predicted scalar, 0 violations.
* Hence modulo scalars the commutative family is spanned by T(1,1,2) [≅ V₃], T(1,1,1,2) and **only the combination T(1,1,2,2) + ½T(1,2,1,2)** of the two quartic traces of bidegree (2,2).

**(c) the triple generates everything.** In the PBW algebra, each Laurent coefficient lies in span_ℤ{V₃, V₄ᵃ, V₄ᵇ} + (polynomials of weight ≤ 4 in the scalar invariants n_s, C₂ˢ, C₃ˢ, C₄ˢ, Ω, C₃^{tot}, C₄^{tot}; 51 monomials, rank 54 with the triple). Coefficients (V₃, V₄ᵃ, V₄ᵇ) of the t-free hat (S^{(1)}_{k,j} = t^{j−k}·hat; S^{(2)}_{k,j} = (−t)^{j−k}·swap(hat)):

| operator | at θ₁ | at θ₂ (hat swapped) |
|---|---|---|
| S_{2,1} | (0,0,0) [= η/t, a pure scalar] | (0,0,0) |
| S_{3,1} | (−2, 0, 0) | (2, 0, 0) |
| S_{3,2} | (1, 0, 0) | (−1, 0, 0) |
| S_{4,1} | (0, −3, −3) | (0, −3, −3) |
| S_{4,2} | (0, 0, 1) | (0, 3, 2) |
| S_{4,3} | (0, 1, 0) | (0, −2, −1) |

(θ₂ coefficients are therefore linear combinations of the θ₁ ones; the independent non-scalar content of the whole family is exactly V₃, V₄ᵃ, V₄ᵇ, i.e. one cubic and two quartic invariants as predicted classically.) Polarised traces in the module spanned by V₃, V₄ᵃ, V₄ᵇ, n₁V₃, n₂V₃ (mod scalars): T(1,2,2) ≡ −V₃; T(1,1,1,2) ≡ −V₄ᵃ + (n₁−3)V₃; T(1,2,2,2) ≡ 3V₃ + 2V₄ᵃ + V₄ᵇ − n₂V₃; T(1,1,2,2) + ½T(1,2,1,2) ≡ −V₄ᵃ − V₄ᵇ + (n₂ − n₁)V₃; **T(1,1,2,2) and T(1,2,1,2) are not in the span individually.**

**(d) naive traces do not commute** (exact, 22 sectors with m ≥ 2 of [3,2,1,0]² and [4,2,1,0]²; entries = number of sectors where [T, V] ≠ 0): T(1,1,2), T(1,2,2), T(1,1,1,2), T(1,2,2,2) and T(1,1,2,2) + ½T(1,2,1,2): 0 for V₃, V₄ᵃ, V₄ᵇ; **T(1,1,2,2): 20, 17, 17; T(1,2,1,2): 20, 17, 17.** The quantum correction is the relative weight ½.

**(e) site exchange (λ = λ′).** P: (ν,q) → (q,ν) preserves every HW space (exact) and conjugates the family: **P V₃ P = −V₃ + c, P V₄ᵃ P = −2V₄ᵃ − V₄ᵇ + c′, P V₄ᵇ P = 3V₄ᵃ + 2V₄ᵇ + c″** with scalars c, c′, c″ (57/57 sectors, exact). So V₃ and 3V₄ᵃ + V₄ᵇ are P-odd and V₄ᵃ + V₄ᵇ is P-even (mod scalars). On M_μ with p₊ symmetric and p₋ antisymmetric copies the odd operator V₃ − c/2 has a symmetric spectrum and **at least |p₊ − p₋| zero modes: it is forced to be degenerate when |p₊ − p₋| ≥ 2** — observed for all 6 such sectors among the 57 ((4,4,2,2), (6,4,2,2), (5,5,2,2), (6,6,2,2), (6,4,4,2), (5,5,3,3)); the 4 other degenerate sectors ((5,4,2,1) of [3,2,1,0]², (7,6,2,1), (7,5,3,1), (6,5,3,2) of [4,3,1,0]²) have an additional, non-forced degeneracy. For λ ≠ λ′ there is no such symmetry and V₃ was simple on every sector tested ([3,2,1,0]⊗[4,2,1,0], [2,2,1,0]⊗[2,1,1,0], [3,2,0,0]⊗[2,1,1,0]).

## 7. G2.5: term lists

`Paul/Mathematica/Gaudin/ref/gl4-gaudin-ops.txt` (87 kB, 1861 lines): `OP name` followed by `TERM {{power,p/q},…} {{site,i,j},…}` lines (TAB-separated; coefficient = exact Laurent polynomial in t; word = product in order, leftmost factor acts last; site-1 factors first, each site's internal order = cdet order; e¹ and e² commute so WP-G3 may normal-order the site-2 factors alone). Operators and numbers of terms: S_2_1: 24, S_2_2: 15, S_3_1: 156, S_3_2: 88, S_3_3: 36, S_4_1: 412, S_4_2: 324, S_4_3: 152, S_4_4: 49 (the θ₁ Laurent coefficients S_{k,j} including their powers of t), and **V3: 88, V4a: 152, V4b: 324** (the t-scaled triple, power t⁰). Header: conventions, definition of x and t, the θ₂ site-swap rule, homogeneity, the N = 3 check, counts, maximal word lengths per operator (site 1/site 2): V₃ 2/1, V₄ᵃ 3/1, V₄ᵇ 3/2; V₃ and V₄ᵃ have exactly one site-2 factor per term (16 distinct site-2 words: the 16 generators e²_ik), V₄ᵇ one or two (88 distinct site-2 words: 16 single + 72 pairs).
`load_terms(path)` returns {name: {(w₁,w₂): {power: Fraction}}}. Re-verification (stage `terms`): the file reloads to exactly the generated lists; the loaded V₃, V₄ᵃ, V₄ᵇ reproduce the restricted matrices on [3,2,1,0]², μ = (5,4,2,1) (m = 4, V₃ degenerate) and [4,2,1,0]², μ = (6,4,3,1) (m = 5), commute exactly, and have a simple joint spectrum (exact discriminant, all roots real; min joint distance 15.72 and 14.92); the nine S_{k,j}(θ₁) rebuilt from the file at t = 31/35 equal the directly built operators on [2,1,1,0]².

## 8. Timings (this machine; exact Fractions throughout)

Stage `sym` ≈ 3 s, `scaffold` 2 s, `gl3` 2.4 s, `comm` ≈ 25 s (incl. the convention variants; the dim-400 exact products 3–6 s per t), `triple` ≈ 25 s for the 11 products up to [4,3,1,0]² (dim 30 625: 11.8 s, [4,2,1,0]²: 4.7 s, all 7 required products together 2 s), `struct` ≈ 30 s, `terms` 1 s; `big`: [5,3,1,0]² 62 s (32 sectors, n_W ≤ 990, all 18 operators restricted), [5,3,2,0]² 44 s, [6,4,2,0]² 111 s (56 sectors, n_W ≤ 4050, triple only). The large cases are cheap because only the HW null spaces (sparse RREF over ℚ, n_W ≤ 4050) and m × m matrices are needed, never the full two-site matrices. A full `all` run takes 301 s.

## 9. Surprises, caveats, what is not done

* **Surprise 1 (the main one):** for gl(4) the cubic operator alone does not resolve the multiplicity spaces (it did for gl(3)); the degeneracy is partly forced by the site-exchange symmetry (§6e) and shows up as soon as a sector has |p₊ − p₋| ≥ 2 (first at [3,2,1,0]², (4,4,2,2)). The joint spectrum of the triple is simple everywhere we looked; {V₃, V₄ᵇ} suffices in all 14 products, {V₃, V₄ᵃ} not (sector (9,7,5,3), m = 12, of [6,4,2,0]²).
* **Surprise 2:** everything is a very short exact calculation: the family is spanned by V₃, V₄ᵃ, V₄ᵇ and Casimir scalars with small integer coefficients (§6c), and V₃, V₄ᵃ, V₄ᵇ equal the polarised antisymmetrised forms ΣD(2,1), D(3,1), D(2,2) − D(3,1) up to the lower antisymmetrised forms and multiples of η = t S_{2,1} (§6a,b); the cubic correction coefficient is (N−2)/2 for gl(N) (exactly solved for N = 3, 4, 5, 6: ½, 1, 3/2, 2).
* **Surprise 3:** the row-ordered determinant is not invariant (it would be for the transposed Lax matrix); the derivative-left normal ordering is harmless.
* Brief typo: dim[2,2,1,0] = 20 (not 60).
* Evidence is exhaustive on the listed products (all sectors with m ≥ 2, m up to 12), **not a proof for all λ, λ′**; simple joint spectrum for generic evaluation points is the MTV-type theorem, here only θ₁ − θ₂ = t symbolic (the triple is t-independent). V₃-degeneracy is certified exactly (discriminant 0); the "NO" entries for subsets use two independent generic combinations.
* Not done: (i) no comparison with the Baxter/Bethe solver of July 2026 (WP-G5) — the exact triples (output file) are the data for it; (ii) no slice reduction (WP-G3) — the term lists have ≤ 1 (V₃, V₄ᵃ) resp. ≤ 2 (V₄ᵇ) site-2 factors per term, so the reduction needs the normal ordering of at most two site-2 generators; (iii) the `others` check (all 18 coefficients on every sector) was skipped for [6,4,2,0]² (time); (iv) no float64/mod-p variant of the HW spaces was needed (exact arithmetic suffices up to n_W = 4050); (v) the θ₂ coefficients are tested through the swap rule and the residue relations, not by an independent re-derivation.

## 10. Hand-off notes for WP-G3 and WP-G5

* Use `Paul/Mathematica/Gaudin/ref/gl4-gaudin-ops.txt` (`V3`, `V4a`, `V4b`; t-free). Reduction to the slice: site-2 factors e²_{ab}: a > b (lowering) → 0, a = b → λ′_a, a < b (raising) → −e¹_{ab} (rightmost first after normal ordering); only V₄ᵇ has two-factor site-2 words.
* Everything needed to test a slice-reduced triple is in the module: `HWExact(R1, R2, mu)` (exact HW basis), `triple_ops(R1, R2)` (Kronecker-sum operators), `restrict_exact`, `analyse_sector` (exact certificates, joint eigenvalue triples).
* For WP-G5: eigenvalue triples for every sector with m ≥ 2 of the products above (18 digits) and the m = 1 values of [3,2,1,0]² are in `gaudin_gl4_ops_output.txt`; the dictionary of the gl(3) case was V = a₁ − η/2, η = |Λ|² − Ω; for gl(4) the analogue is V₃ = t S_{3,2}, V₄ᵃ = t S_{4,3}, V₄ᵇ = t² S_{4,2} with S_{k,j} the Laurent coefficients of the oper's characteristic polynomial; the exact relations of §6 fix the scalars.
