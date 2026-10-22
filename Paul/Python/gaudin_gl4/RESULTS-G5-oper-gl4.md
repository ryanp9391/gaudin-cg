# RESULTS-G5 (agent G5, 2026-10-22): the oper side of the gl(4) missing labels — the exact dictionary (V₃, V₄ᵃ, V₄ᵇ) ↔ the 4th-order Fuchsian oper, its polynomial solutions, the independent exact Q-system solve, and the July-2026 gl(4) Gaudin ODE

Files (all new; nothing existing was edited; all in `Paul/Python/gaudin_gl4/` unless stated):
`oper_gl4.py` (all code; CLI `python -X utf8 oper_gl4.py dict | exps | solve | search | wb | big | all`, plus `wbsym` = the 150 s symbolic derivation of §5.3, not in `all`),
`test_oper_gl4.py` (7 regression tests, ≈ 5 s, `python -X utf8 test_oper_gl4.py`: **7 tests, 0 failed**), `oper_gl4_output.txt` (the printed output of one run of `all`), this file, and
`../../Mathematica/Gaudin/ref/gl4-oper-dictionary.txt` (the dictionary, 135 lines, §2.6).
The printed output of `all` (`oper_gl4_output.txt`) took 516 s in total on this machine (the stages: dict 11 s, exps ≈ 45 s, solve ≈ 175 s, search ≈ 60 s, wb ≈ 170 s, big ≈ 100 s).
Python 3.13.13, numpy 2.5.2, scipy 1.18.1, sympy 1.14 (polynomial rings, Gröbner bases, `linsolve`), mpmath 1.3 (60 digits); i7-12700KF, Windows 11; ≤ 2 heavy processes at any time (one extra worker process with a time-box for Gröbner bases of the single conditions, §4.3).
Read-only imports: `gaudin_gl4_ops` (G2: `Rep`, `HWExact`, `build_ops`, `triple_ops`, `restrict_exact`, `joint_spectrum`, `racah_speiser`, `GENERIC_COMBO`), `glN_reps` (G1: `GLNRep` for the Casimir eigenvalues, `racah_speiser`, `weyl_dim`).
The July-2026 notebook `GaudinSolver.wb` is read as JSON text. **Everything quantitative below is exact (Fractions / sympy rationals / exact Gröbner bases) unless a number is quoted with digits (mpmath, 60 digits).**

## 0. Headline results

1. **The dictionary is forced, not fitted, and every one of the 18 scalars has a closed form.** With y_c(ν) = ν_c + 4 − c and F_ν(θ) = Π_c(θ − y_c) = Σ_k a_k(ν) θ(θ−1)…(θ−k+1), put **T_k(ν) = (−1)^k a_{4−k}(ν)** (T₁ = |ν|, T₂ = e₂(y) − 3e₁(y) + 7, T₃ = e₃(y) − e₂(y) + e₁(y) − 1, T₄ = e₄(y) = Π y_c) and **τ_k = T_k(μ) − T_k(λ) − T_k(λ′)**. Then (t = θ₁ − θ₂; hats are t-free)
   S⁽¹⁾_{k,j} = t^{j−k}(a·V₃ + b·V₄ᵃ + c·V₄ᵇ + scalar), S⁽²⁾_{k,j} = (−t)^{j−k}(a′·V₃ + b′·V₄ᵃ + c′·V₄ᵇ + scalar′), with the integer coefficients of G2's table (§6c) and the scalars **T_k(λ), T_k(λ′), ±τ_k** (table in §2.1). The scalars are the top poles (central) and the three differences τ₂ = n₁n₂ − Ω, τ₃, τ₄ of the central characters of V_μ, V_λ, V_λ′ — in trace-form Casimirs T₂ = −C₂/2 + n²/2 + 3n/2, T₃ = C₃/3 − nC₂/2 − C₂/3 + n³/6 + 5n²/6 + 2n, T₄ = (§2.3). The table is a theorem: it is the solution of the linear relations "S_k = O(u⁻ᵏ) at ∞ and the coefficient of u⁻ᵏ is the total Casimir T_k(μ)", whose free parameters are exactly V₃, V₄ᵃ, V₄ᵇ (k − 2 per k: 0, 1, 2). **Verified on 95 sectors (Σm = 146 copies) of nine products (the six required, [3,2,0,0]⊗[2,1,1,0], [2,2,1,0]², [3,2,0,0]²): 1710 (sector, operator) tests, remainder = the scalar matrix of the closed form, 0 violations**; 0 violations also at t = 31/35 and t = −3/7 (46 further sector tests).
2. **Exponents, degrees, polynomial solutions (G5.2): 210 copies of 134 sectors (11 products; 88 copies with rational triple treated in exact arithmetic), all pass.** The indicial polynomials at θ₁, θ₂, ∞ are v-independent (symbolically, as affine forms) and equal Π(e − (λ_c+4−c)), Π(e − (λ′_c+4−c)), Π(n − (μ_c+4−c)) in 134/134 sectors; the residue relations hold identically in v; the null space of D(v) on polynomials of degree ≤ d_j has dimension #{k : d_k ≤ d_j} = 5 − j for all j in **210/210** copies (exact rank, or 60-digit SVD with a gap ≥ 3·10⁵⁹ between kept and dropped singular values); the four polynomial solutions have Wronskian ∝ (u−θ₁)^{|λ|}(u−θ₂)^{|λ′|} in 210/210 (deviation ≤ 2·10⁻⁵⁹), and their **measured** vanishing orders at θ₁, θ₂ are {λ_c+4−c}, {λ′_c+4−c} in 210/210. t-independence: 56/56 further exact tests at t = 31/35, −3/7. Negative tests: all 630 perturbed and 210 random triples rejected; 352/352 exact (threshold-free) rejections on the rational triples; **the brief's literal criterion (σ_min/σ_max ≥ 10⁻⁶) is NOT met by 615 of the 630 perturbed triples** (§3.5: the criterion is not meaningful at that scale; the degree-exactly-d₁ residual separates true and perturbed triples by > 47 orders of magnitude).
3. **The independent solve (G5.3): Σ = {v : all four polynomial-solution conditions hold} equals G2's joint spectrum EXACTLY in all 165 sectors solved (Σm = 279 copies, 64 sectors with m ≥ 2, m up to 5), and in the 15 largest-multiplicity sectors of [5,3,1,0]², [5,3,2,0]², [6,4,2,0]² (m up to 12).** Method: Frobenius recursion at ∞ → polynomial conditions in v (resonances = conditions, free parameters gauge-fixed) → Gröbner basis (sympy); every condition vanishes exactly on the commuting G2 matrices (so J ⊆ K), and dim_ℚ ℚ[v]/J = m = the number of distinct points (trace-form rank) ⇒ J = K: Σ is exactly the m triples, each a simple solution (Jacobian rank 3 in every case). The points computed from J alone agree with G2's 60-digit triples to ≤ 2.8·10⁻⁵⁸ (exactly for m = 1). The interesting sector (5,4,2,1) of [3,2,1,0]² (m = 4, V₃ = 57 twice) is solved like all others: 4 distinct points. An independent numerical search (240 random real starts per sector, Levenberg–Marquardt) finds **no solution other than G2's triples** in 12 large sectors (§4.4).
4. **gl(4) pencil-bug analysis (§4.3).** On the 68 sectors of the eight small products no single condition cuts out the spectrum in general (j = 4 alone: a curve in all 68 sectors; {3}: 0/68 exact; {2}: 6/68; {1}: 7/68), the pairs reach 17–32/68 ({4,3} best: 32/68, 6 of the 16 sectors with m ≥ 2), the triples 38–60/68 ({4,3,2}: 60/68), and **only all four conditions cut out exactly the spectrum in every sector (68/68)**; in the 8 sectors where {4,3,2} fails (all have μ₁ = μ₂) the condition for P₁ is indispensable ((4,4,3,1) of [3,2,1,0]²: P₁ alone gives exactly the 3 triples) — the exact analogue of "gcd(P₁,P₂) over-counts and P₃ is needed" of gl(3).
5. **Structural by-product:** on every sector with μ₄ = 0 (and λ₄ = λ′₄ = 0) the oper forces **V₄ᵃ = V₄ᵇ = 0 exactly** (the constant polynomial must solve D·1 = s₄ = 0), so there the 3-parameter solve is the 1-parameter gl(3)-type solve in V₃ (checked on all 47 such sectors of the 11 products of §3, 54 copies, 7 with m ≥ 2: G2's triples are (v₃, 0, 0) to < 10⁻⁴⁰ — a one-off check; `test_mu4_zero` in the test file repeats it for the five small products).
6. **July-2026 gl(4) Gaudin ODE (G5.4):** transcribed mechanically from the notebook JSON; polynomial solutions of degrees degQ exist for all 28 irreps of the five (S,P) = (1,1), (2,1), (2,2), (3,2), (3,3) and only at those degrees (56/56 tests). **gaudinODE = 8 (z−w₁)²(z−w₂)² · Π^p ∘ D† ∘ Π^{−p} with p = (S, P)** (S at w₁ = θ₁, P at w₂ = θ₂, D† = formal adjoint of OUR oper): the matching system is linear in our unknowns (v₃, v₄ᵃ, v₄ᵇ, T_k(μ)), has a **unique** consistent gauge among the 16 (resp. 50) scanned (p₁, p₂, site assignment), and its unique solution equals G2's exact triple in 28/28 irreps; the two operators coincide coefficient by coefficient (exact rational functions of z, 28/28). **Dictionary (derived with S, P, w₁, w₂, C2, C4 all symbolic):** with N = S+P, **v₃ = (S+1)(2N(N+2) − Cas2)/2, v₄ᵃ = 0, v₄ᵇ = T₄(μ) = (Cas2² + 2(N+1)Cas2 − 2Cas4 + 4N(N+2)(N+3))/8**, T₂(μ) = (4N²+6N−Cas2)/2, T₃(μ) = (N+2)(2N(N+2)−Cas2)/2. The degrees map as degQ_k = (S+P+3) − (μ_k+4−k) = d_{5−k}; the adjoint solutions W₃[P_k, k≠j]/((z−w₁)^S(z−w₂)^P) built from OUR polynomial solutions are polynomials of degree degQ_j and solve the notebook ODE exactly (112/112: 4 per irrep, 28 irreps).

## 1. Conventions and what is taken from G2

Binding conventions of G2 (= notes §1.1 generalised): [E_ij,E_kl] = δ_jk E_il − δ_li E_kj; two sites V_λ (site 1, pole θ₁) ⊗ V_λ′ (site 2, pole θ₂); L(u) = e¹/(u−θ₁) + e²/(u−θ₂); t = θ₁ − θ₂; cdet(∂ − L) = ∂⁴ − S₁∂³ + S₂∂² − S₃∂ + S₄ (column ordered, derivative to the right); on a joint eigenvector the S_k become the scalar functions
**s_k(u) = Σ_{s=1,2} Σ_{j=1}^{k} S⁽ˢ⁾_{k,j}/(u−θ_s)^j** (exact; s₁ = |λ|/(u−θ₁) + |λ′|/(u−θ₂)), and **D = ∂⁴ − s₁∂³ + s₂∂² − s₃∂ + s₄** is the oper. G2 facts used: the top poles are central scalars S_{m,m} = (−1)ᵐ a_{4−m}; S_k = O(u⁻ᵏ) (six residue relations, identities in U(gl₄)⊗U(gl₄)); θ₂ = site swap and t → −t; V₃ = tS_{3,2}, V₄ᵃ = tS_{4,3}, V₄ᵇ = t²S_{4,2} are t-free, commute, have simple joint spectrum on every sector (the data of G2's output). In the code **θ₂ = 0, θ₁ = t** (t = 1 by default; the scaling u → u/t is exact, tests at t = 31/35, −3/7 included). All Laurent data, expansion coefficients at ∞ and pencil entries are **affine forms in v = (v₃, v₄ᵃ, v₄ᵇ)** (tuples of Fractions), so every statement "for all v" is exact.

## 2. G5.1: the dictionary

### 2.1 Derivation (exact) and the table

Unknowns: S⁽¹⁾_{k,j} =: A_{k,j}, S⁽²⁾_{k,j} =: B_{k,j}, j = 1..k−1 (2(k−1) numbers per k); known: A_{k,k} = T_k(λ), B_{k,k} = T_k(λ′) (central). The expansion of s_k at ∞ has the coefficient of u⁻ʳ equal to Σ_s Σ_{j≤r} C(r−1,j−1) θ_s^{r−j} S⁽ˢ⁾_{k,j}: it vanishes for r < k (G2) and equals **T_k(μ) for r = k** (at ∞ L → E^tot/u + O(u⁻²), so the leading coefficient of S_k is the single-site word polynomial with E → E^tot, i.e. S_{k,k} with λ → μ; on V_μ it is T_k(μ)). These are k relations, triangular in the highest Laurent index, so 2(k−1) − k = **k − 2 free parameters: 0, 1, 2 for k = 2, 3, 4** — exactly one cubic and two quartic missing labels, which we choose as A_{3,2} = V₃/t, A_{4,3} = V₄ᵃ/t, A_{4,2} = V₄ᵇ/t² (their definition). Solving (with θ₁ − θ₂ = t; the θ's drop out) gives:

| k | j | at θ₁: t^{k−j} S⁽¹⁾_{k,j} = | at θ₂: (−t)^{k−j} S⁽²⁾_{k,j} = |
|---|---|---|---|
| 2 | 2 | T₂(λ) | T₂(λ′) |
| 2 | 1 | τ₂ | τ₂ |
| 3 | 3 | T₃(λ) | T₃(λ′) |
| 3 | 2 | V₃ | −V₃ + τ₃ |
| 3 | 1 | −2V₃ + τ₃ | 2V₃ − τ₃ |
| 4 | 4 | T₄(λ) | T₄(λ′) |
| 4 | 3 | V₄ᵃ | −2V₄ᵃ − V₄ᵇ + τ₄ |
| 4 | 2 | V₄ᵇ | 3V₄ᵃ + 2V₄ᵇ − τ₄ |
| 4 | 1 | −3V₄ᵃ − 3V₄ᵇ + τ₄ | −3V₄ᵃ − 3V₄ᵇ + τ₄ |

(the integer parts are exactly G2's table §6c, which is thereby explained: the θ₂ column is the θ₁ column rewritten with the relations; S⁽²⁾_{2,1} = −S⁽¹⁾_{2,1}, S⁽²⁾_{3,1} = −S⁽¹⁾_{3,1}, S⁽²⁾_{4,1} = −S⁽¹⁾_{4,1}). All 18 scalars are therefore closed forms in T_k(λ), T_k(λ′), τ_k; none remains tabulated. (`DICT_TABLE`, `dictionary()` in the code.)

### 2.2 Verification (every sector, every operator)

For every irrep μ (m ≥ 1) of the products below the 18 Laurent coefficients (G2's `build_ops`, θ₂ by swap/(−t)) are restricted to M_μ (G2's `HWExact`, `restrict_exact`, faithful: 0 violations), divided by t^{j−k}, the integer combination of the restricted triple is subtracted, and the remainder is required to be a scalar matrix equal to the closed form. Result (`oper_gl4.py dict`, 10 s):

| product | dim | sectors | Σm | violations |
|---|---|---|---|---|
| [2,1,0,0]² | 400 | 7 | 8 | 0 |
| [2,1,1,0]² | 225 | 6 | 7 | 0 |
| [3,2,1,0]² | 4096 | 16 | 30 | 0 |
| [2,2,1,0]⊗[2,1,1,0] | 300 | 6 | 7 | 0 |
| [1,1,0,0]⊗[2,1,0,0] | 120 | 4 | 4 | 0 |
| [4,2,1,0]² | 19600 | 27 | 56 | 0 |
| [3,2,0,0]⊗[2,1,1,0] | 900 | 7 | 8 | 0 |
| [2,2,1,0]² | 400 | 7 | 8 | 0 |
| [3,2,0,0]² | 3600 | 15 | 18 | 0 |
| **total** | | **95** | **146** | **0 of 1710 (sector, operator) tests** |

and at t = 31/35 and t = −3/7 on the 23 sectors of the four small products: 0 violations each (the t-power bookkeeping). For m ≥ 2 the remainder is a scalar matrix (not only a scalar diagonal); e.g. [3,2,1,0]², μ = (5,4,2,1), m = 4, all 18 tests.

### 2.3 Closed forms in the Casimirs

The central characters are explicit in the shifted weights: **T₂ = e₂(y) − 3e₁(y) + 7, T₃ = e₃(y) − e₂(y) + e₁(y) − 1, T₄ = e₄(y) = Π_c y_c** (e_k = elementary symmetric functions of y_c = ν_c + 4 − c; conversion of Π(θ − y_c) to falling factorials; asserted on all 70 partitions with ν₁ ≤ 4 in `test_oper_gl4.py`). In the trace-form Casimirs of G2 (n = |ν|, C_k = Σ e_{i₂i₁}e_{i₃i₂}⋯e_{i₁i_k} evaluated on the HW vector of `GLNRep`; exact fit on 63 partitions with ν₁ ≤ 4, |ν| ≤ 12 — ranks 4, 7, 12 of the monomial libraries of weight ≤ k in (n, C₂, C₃, C₄) — and verified on 20 held-out partitions with ν₁ ∈ {5, 6}):

* T₂ = −C₂/2 + n²/2 + 3n/2,
* T₃ = C₃/3 − n C₂/2 − C₂/3 + n³/6 + 5n²/6 + 2n,
* T₄ = −C₄/4 + C₂²/8 + n C₃/3 − n²C₂/4 − C₃/2 + n C₂/6 − 3C₂/4 + n⁴/24 + n³/12 + 7n²/8 + 3n/2.

Hence τ₂ = n₁n₂ − Ω with **Ω = (C₂(μ) − C₂(λ) − C₂(λ′))/2** (G2's tS_{2,1} = n₁n₂ − Ω, now a consequence), and τ₃, τ₄ are the same differences of T₃, T₄ (n_μ = n₁ + n₂); all verified on all 95 sectors (0 violations). For gl(3) the analogue was V = a₁ − η/2; here V₃ = tS_{3,2} exactly and the "η" are τ₂, τ₃.

### 2.4 Worked example

[3,2,1,0]², μ = (6,3,2,1) (m = 2; T(λ) = (15, 15, 0), T(μ) = (57, 135, 135), τ = (27, 105, 135)); t = 1; the two copies (v₃, v₄ᵃ, v₄ᵇ) = (55.854…, 54.271…, 5.729…) and (49.146…, 20.729…, 39.271…). Copy 1: S⁽¹⁾_{2,·} = (15; 27), S⁽²⁾_{2,·} = (15; −27); S⁽¹⁾_{3,3,2,1} = (15, 55.854, −6.708), S⁽²⁾_{3,3,2,1} = (15, −49.146, +6.708); S⁽¹⁾_{4,4..1} = (0, 54.271, 5.729, −45), S⁽²⁾_{4,4..1} = (0, −20.729, 39.271, +45); copy 2 has the θ₁ and θ₂ data exchanged up to the signs (the two copies are the symmetric/antisymmetric pair under site exchange, G2 §6e).

### 2.5 Which scalars are closed forms, which are tabulated

**All 18 are closed forms** (tops: T_k(λ), T_k(λ′); the other 8 non-triple coefficients: ±τ_k, plus the integer combinations of the triple); none is only tabulated. Per-sector tables of T_k(λ), T_k(λ′), T_k(μ), τ_k are in the ref file. The only non-closed objects are, as in G2, the joint eigenvalues of the triple themselves (the point of §4).

### 2.6 The file

`ref/gl4-oper-dictionary.txt`: header (conventions, derivation, format), 18 `DICT k j site tpow a b c scalar` lines (S⁽ˢⁱᵗᵉ⁾_{k,j} = (±t)^{tpow}(a v₃ + b v₄ᵃ + c v₄ᵇ + scalar)), 3 `CASFORM` lines (T_k in n, C₂, C₃, C₄), the closed-form dictionary of the multiplicity-free family of §5 in the header, and 95 `SECTOR` lines (λ, λ′, μ, m | T_k(λ) | T_k(λ′) | T_k(μ) | τ_k | the m joint triples of G2 to 25 digits) — the data of the oper for the 146 copies.

## 3. G5.2: exponents, polynomial solutions, Wronskian, negative tests

D(v) is built from the dictionary for a sector and a triple v (θ₂ = 0, θ₁ = 1; every number exact or 60-digit). Products: the six required plus [3,2,0,0]⊗[2,1,1,0], [2,2,1,0]², [3,2,0,0]², [3,2,1,0]⊗[4,2,1,0], [4,2,2,0]² (11 products, **134 sectors, 210 copies**, 88 of them with rational triple). `oper_gl4.py exps`, 42 s.

**(i) local exponents.** The indicial polynomial at θ_s from the most singular terms is e^{(4)} − S_{1,1}e^{(3)} + S_{2,2}e^{(2)} − S_{3,3}e + S_{4,4} (falling factorials) = Π_c(e − y_c(λ or λ′)): the coefficients are the central top poles, hence **independent of v** (checked as an identity of affine forms in v): 134/134 sectors at θ₁ and at θ₂. Measured: the vanishing orders of the span of the four polynomial solutions at θ₁ and θ₂ (row reduction of the Taylor matrices at u = θ_s, exact or with a 10⁻²⁵ rank tolerance at 60 digits) are exactly {λ_c+4−c} and {λ′_c+4−c} in **210/210** copies (e.g. {0,2,4,6} for [3,2,1,0]²).
**(ii) infinity.** D uⁿ = I(n)u^{n−4} + …, I(n) = n^{(4)} − σ₁n^{(3)} + σ₂n^{(2)} − σ₃n + σ₄ (falling factorials n^{(k)}) with σ_k the coefficient of u⁻ᵏ in s_k, which equals T_k(μ) (§2.1) whatever v: **I(n) = Π_c(n − (μ_c+4−c)) identically in v** (134/134), so polynomial solutions have degrees d_j = μ_j + 4 − j (d₁ > d₂ > d₃ > d₄ strictly). The Fuchs relation Σ exponents = (|λ|+6) + (|λ′|+6) − (|μ|+6) = 6 holds.
**(iii) polynomial solutions of every copy.** The pencil of D̃ = (u−θ₁)⁴(u−θ₂)⁴D on polynomials of degree ≤ d is a (d+5) × (d+1) matrix, affine in v: D̃p is a polynomial of degree ≤ d+4 (each (u−θ₁)⁴(u−θ₂)⁴ s_k ∂^{4−k} has degree ≤ (8−k) + d − (4−k) = d+4), whose top coefficient is I(d)p_d, vanishing identically for d = d_j; the null space must have dimension #{k : d_k ≤ d_j} = 5 − j. **210/210 copies pass**: exact rank on the 88 rational triples (dimensions (4,3,2,1) for j = 1..4, gauge-fixed solutions of exact degree d_j), and for the 122 irrational copies the singular values at 60 digits: worst gap between the smallest kept and the largest dropped singular value **3.1·10⁵⁹** (criterion ≥ 10²⁰), worst dropped singular value / entrywise scale 3.5·10⁻⁶², worst residual of the gauge-fixed least-squares solutions 1.8·10⁻⁵⁶.
**(iv) Wronskian.** The Wronskian of the four gauge-fixed solutions (p_{d_j} = 1, p_{d_k} = 0 for k > j) has degree |μ| = Σd_j − 6 and equals c·(u−θ₁)^{|λ|}(u−θ₂)^{|λ′|}: exactly (coefficient by coefficient) on the 88 exact copies, deviation ≤ 1.9·10⁻⁵⁹ on the others: **210/210** (Abel's identity holds trivially since s₁ = |λ|/(u−θ₁) + |λ′|/(u−θ₂); the content is that the four solutions are polynomial and independent).
**t-independence:** the same exact triples at t = 31/35 and t = −3/7 give the null dimensions (4,3,2,1): 56/56 tests.
**Residue relations:** σ_{k,r}(v) = 0 for r < k and S⁽²⁾_{2,1} = −S⁽¹⁾_{2,1} as identities in v: 134/134.

### 3.5 Negative tests (and a criterion that fails as stated)

The brief's measure — smallest singular value of M₁(v′) divided by the largest ≥ 10⁻⁶ for v′ = known triple + 10⁻³ in one component or a random triple — is **not met by 615 of the 630 perturbed triples** (random triples: 2.9·10⁻⁶ at worst, so it holds for those). Two reasons, neither a failure of the rejection: (a) when the perturbed component is v₃ (or v₄ᵃ, v₄ᵇ) and μ₄ = 0, the lower-degree polynomial solutions that do not depend on that component (the constant 1 solves D·1 = s₄ = 0 whatever v₃) survive, so σ_min = 0 exactly although no solution of degree d₁ exists; (b) the first-order response of the four lifted null singular values to a shift δ of one component of v is δ·‖M_i P_null‖ ≈ 10⁻³ × O(1–10), while σ_max of the pencil is 10²–10⁴ (entries 5 to 7·10³ at d₁ = 6–11), i.e. σ_min/σ_max ≈ 10⁻⁶–10⁻¹⁰ for the larger sectors — below 10⁻⁶ by construction, although the null singular values of the true triples are 10⁻²⁸–10⁻⁶² (kept ones ≥ 60). The correct test is that **no solution of degree exactly d₁ exists** (the column uᵈ¹ is not in the span of the other columns): least-squares residual relative to ‖M₁‖_F: **true triples ≤ 1.3·10⁻⁶⁰; the 630 perturbed triples ≥ 1.0·10⁻¹² (a separation of > 47 orders of magnitude); the 210 random triples ≥ 2.1·10⁻⁵**; and, threshold-free, **exact rank tests on the 88 rational triples: the true triples admit the d₁-solution in 88/88, all 352 tests (each triple + 1/1000 in each component, and one random rational triple per copy) are rejected in 352/352**.

## 4. G5.3: the independent solve

### 4.1 Setting and method

D(v) is affine in v; the condition "a polynomial solution of degree exactly d_j exists" is a rank condition on the (d_j+5) × (d_j+1) pencil M_j(v) (sizes explained in §3 (iii)); Σ = {v : all four conditions}. By Frobenius at ∞: D uⁿ = Σ_R I_R(n) u^{n−4−R}, I_R(n) affine in v; for a polynomial p = Σ p_i uⁱ of degree d the d+5 equations c_m := Σ_{q=0}^{m} p_{d−q} I_{m−q}(d−q) = 0 (m = 0..d+4) are equivalent to D̃p = 0. m = 0 is I₀(d)p_d = 0 (trivial); for 1 ≤ m ≤ d the equation determines p_{d−m} unless d − m = d_k (k > j) is another exponent — a **resonance**: I₀(d_k) = 0, the equation is one **consistency condition** on v (a polynomial of degree ≤ 9 in the 165 sectors of §4.2, ≤ 10 in the big ones) and p_{d_k} is a free parameter, **gauge-fixed to 0** (legitimate when the solution of degree d_k exists, i.e. in the intersection; this is the "echelon gauge" of the gl(3) Q-system reports); for m > d the 4 **bottom equations** are conditions. Per j: 3, 2, 1, 0 resonance conditions for j = 1..4 plus 4 bottom equations: 7, 6, 5, 4 polynomial conditions in v = (v₃, v₄ᵃ, v₄ᵇ) (22 in all; degrees ≤ 9 in the 165 sectors, ≤ 10 in the big ones, while d₁ ≤ 11 resp. 13), computed in sympy rings over ℚ. Σ = V(J), J = the ideal of all 22 conditions.
**Exactness argument.** Let K = ker(ℚ[v] → Mat_m(ℚ), v ↦ (V₃, V₄ᵃ, V₄ᵇ restricted to M_μ)): since the joint spectrum is simple (G2, exact certificates), K is a radical ideal of m points and dim ℚ[v]/K = m. Each generator of J is evaluated **exactly on the commuting matrices** (V₃, V₄ᵃ, V₄ᵇ): it vanishes (so J ⊆ K, and the G2 triples lie in Σ) in every sector. Then **dim_ℚ ℚ[v]/J = #standard monomials of the Gröbner basis (grevlex) = m ⟹ J = K**, i.e. Σ is exactly the m triples. The number of distinct points is the rank of the trace form on the quotient algebra (= m), the points are recomputed from J alone (joint eigenvalues of the multiplication matrices of ℚ[v]/J, 60 digits) and compared with G2, and the Jacobian of the 22 conditions has rank 3 at every G2 triple (60-digit SVD): each solution is simple.

### 4.2 Results (`oper_gl4.py solve`, 173 s; `big`, ≈ 2 min)

"Agreement" = the number of the 60 digits in which the points of Σ computed from J alone (joint eigenvalues of the multiplication matrices of ℚ[v]/J) agree with G2's joint triples, worst sector of the product ("identical": all digits, exact rational points for m = 1); "Jacobian rank" is 3 at every triple of every sector.

| product | dim | sectors | Σm | max m | Σ = G2 exactly | agreement (worst sector) | Σ over sectors of Gröbner time |
|---|---|---|---|---|---|---|---|
| [2,1,0,0]² | 400 | 7 | 8 | 2 | 7/7 | identical (60 digits) | < 0.1 s |
| [2,1,1,0]² | 225 | 6 | 7 | 2 | 6/6 | ≥ 60 digits | 0.1 s |
| [3,2,1,0]² | 4096 | 16 | 30 | 4 | 16/16 | ≥ 59 digits | 2.1 s |
| [2,2,1,0]⊗[2,1,1,0] | 300 | 6 | 7 | 2 | 6/6 | identical | 0.3 s |
| [1,1,0,0]⊗[2,1,0,0] | 120 | 4 | 4 | 1 | 4/4 | identical (rational) | < 0.1 s |
| [4,2,1,0]² | 19600 | 27 | 56 | 5 | 27/27 | ≥ 59 digits | 6.2 s |
| [3,2,0,0]⊗[2,1,1,0] | 900 | 7 | 8 | 2 | 7/7 | ≥ 60 digits | 0.4 s |
| [2,2,1,0]² | 400 | 7 | 8 | 2 | 7/7 | ≥ 60 digits | 0.5 s |
| [3,2,0,0]² | 3600 | 15 | 18 | 2 | 15/15 | ≥ 59 digits | 1.0 s |
| [3,2,1,0]⊗[4,2,1,0] | 8960 | 20 | 38 | 4 | 20/20 | ≥ 59 digits | 3.5 s |
| [4,2,2,0]² | 7056 | 19 | 26 | 3 | 19/19 | ≥ 59 digits | 5.7 s |
| [4,3,1,0]² | 30625 | 31 | 69 | 5 | 31/31 | ≥ 58 digits | 11.2 s |
| **total** | | **165** | **279** | 5 | **165/165** | **2.8·10⁻⁵⁸ (worst over all)** | |
| [5,3,1,0]², 5 largest-m sectors (m = 7,7,7,6,6) | 129600 | 5 | 33 | 7 | 5/5 | ≥ 57 digits | 2–4 s each |
| [5,3,2,0]², 5 largest-m sectors (m = 8,6,6,6,5) | 90000 | 5 | 31 | 8 | 5/5 | ≥ 56 digits | 2–5 s each |
| [6,4,2,0]², 5 largest-m sectors (m = 12,11,10,10,10) | 531441 | 5 | 53 | **12** | 5/5 | ≥ 55 digits | 12–17 s each |

Every required sector is covered: all sectors of [2,1,0,0]², [2,1,1,0]², [2,2,1,0]⊗[2,1,1,0], [1,1,0,0]⊗[2,1,0,0], **all 16 sectors of [3,2,1,0]²** (including (5,4,2,1), m = 4, where V₃ is degenerate: Σ has 4 distinct points, (v₃, v₄ᵃ, v₄ᵇ) = (67.392, 70.981, −16.981), (46.608, 19.019, 34.981), (57, 50.697, 21.303), (57, 21.303, 50.697), Jacobian rank 3 at each), and **all 27 sectors of [4,2,1,0]²** including the m = 4, 5 sectors (7,4,2,1), (6,5,2,1), (6,4,2,2), (5,4,3,2) (m = 4) and (6,4,3,1) (m = 5): exactly solved, no certification needed. The largest case, **(9,7,5,3) of [6,4,2,0]² (m = 12, n_W = 916, d = (12,9,6,3); G2: the only sector where {V₃, V₄ᵃ} fails to separate)**: dim ℚ[v]/J = 12 = distinct points, deviation from G2's 12 triples ≤ 10⁻⁵⁵ (7·10⁻⁵⁵), 17 s including the HW space and the three restricted operators; the four next-largest sectors of [6,4,2,0]² (m = 11, 10, 10, 10; conditions of degree ≤ 10) are equally exact (12–15 s each). Cost: the Gröbner computation of the full system takes 0.1–0.8 s per sector up to m = 5, a few seconds for m = 6–8, ≈ 10 s for m = 10–12 (HW space and restricted triple included in the times quoted for the big sectors).

### 4.3 The gl(4) version of the pencil-bug analysis

The 15 non-empty subsets S of {1,2,3,4}: J_S = ideal generated by the gauge-fixed conditions of the j ∈ S (so {4} = "P₄ alone", {3} = "P₃ alone with its u^{d₄} coefficient fixed to 0", … — the echelon systems of the gl(3) Q0/J reports). Computed on the 68 sectors of the eight small products (the four required small ones, [3,2,1,0]², [3,2,0,0]⊗[2,1,1,0], [2,2,1,0]², [3,2,0,0]²) with a time-box of 5 s per Gröbner basis (worker process; "T/O" = not finished). Entries: number of sectors in which J_S = K, i.e. S cuts out exactly the spectrum (exact), the same for the 16 sectors with m ≥ 2, number of sectors with a positive-dimensional solution set, with a finite but too large set:

| subset | exact (of 68) | exact, m ≥ 2 (of 16) | dim > 0 | finite, too many | T/O |
|---|---|---|---|---|---|
| {4} | 0 | 0 | 68 | 0 | 0 |
| {3} | 0 | 0 | 54 | 14 | 0 |
| {2} | 6 | 2 | 55 | 7 | 0 |
| {1} | 7 | 2 (of 15) | 36 | 20 | 5 |
| {4,3} | 32 | 6 | 0 | 36 | 0 |
| {4,2} | 24 | 6 | 0 | 44 | 0 |
| {4,1} | 21 | 5 | 0 | 46 | 1 |
| {3,2} | 20 | 4 | 0 | 48 | 0 |
| {3,1} | 27 | 5 | 0 | 39 | 2 |
| {2,1} | 17 | 4 | 0 | 49 | 2 |
| {4,3,2} | 60 | 13 | 0 | 8 | 0 |
| {4,3,1} | 49 | 10 | 0 | 19 | 0 |
| {4,2,1} | 38 | 9 | 0 | 30 | 0 |
| {3,2,1} | 38 | 8 | 0 | 30 | 0 |
| **{4,3,2,1}** | **68** | **16** | 0 | 0 | 0 |

Findings: (a) in gl(4) a single condition essentially never gives a finite answer: {4} is a **curve in all 68 sectors** (on μ₄ = 0 sectors the line v₄ᵃ = v₄ᵇ = 0), unlike gl(3) where one unknown makes each system a polynomial; (b) minimal exact subsets (counted over the 68 sectors, a sector can have several): {4,3} in 32, {3,1} in 20, {4,2} in 18, {4,1} and {3,2} in 14 each, {1} alone in 7, {2} alone in 6, {4,3,2} in 5, {2,1} in 4, {4,3,1} in 2; (c) **no proper subset works in all sectors**; the triple {4,3,2} fails in 8 sectors — (4,4,3,1) of [3,2,1,0]², (4,4,2,2) of [3,2,1,0]², (3,3,2,1) of [2,2,1,0]⊗[2,1,1,0], (2,2,1,0) of [1,1,0,0]⊗[2,1,0,0], (3,3,2,1) of [3,2,0,0]⊗[2,1,1,0], (3,3,3,1) and (3,3,2,2) of [2,2,1,0]², (4,4,2,0) of [3,2,0,0]² — **all of them have μ₁ = μ₂** (adjacent top exponents d₁ − d₂ = 1; but only 8 of the 26 sectors with μ₁ = μ₂ fail), and there the condition for P₁ is needed: P₁ alone gives exactly the spectrum in (4,4,3,1) (3 solutions of 3), (3,3,2,1) (2 of 2 and 1 of 1); {1,3} or {1,4} in the others. This is the gl(4) analogue of Q0/J: "gcd(P₁,P₂) over-counts, P₃ is needed" (E₁₂-cut sectors); (d) the cheapest conditions are the low-degree ones: {4,3} suffices in 32/68 and 6/16 of the m ≥ 2 sectors. The honest "existence of a solution of degree d_j" conditions (with the free parameters p_{d_k} kept) contain the gauge-fixed loci, so every statement "subset S already cuts out exactly the spectrum" above is a fortiori true for the honest versions; their counts of too-many solutions are larger (not computed).

### 4.4 Numerical, independent of the Gröbner computation

`oper_gl4.py search`: Levenberg–Marquardt (float64, residual = the 22 conditions scaled by Σ|coefficient|·(3·max|spectrum|)^{degree}) from 240 (m ≥ 8: 2400) random real starts in |v_i| ≤ 3·max|spectrum|, on the sectors with the largest multiplicity: m = 4 and 3 of [3,2,1,0]² (5,4,2,1), (5,3,3,1), (5,3,2,2), (4,4,3,1), (4,3,3,2) and the V₃-degenerate (4,4,2,2); (7,4,2,1), (6,5,2,1), (6,4,3,1) m = 5, (6,4,2,2), (5,4,3,2) of [4,2,1,0]²; and (9,7,5,3) m = 12 of [6,4,2,0]². Result: in every sector all starts that converge (residual < 10⁻¹⁰: 240/240 except 169/240 at (7,4,2,1); 2400/2400 at the m = 12 sector) land on G2's triples and **no solution other than G2's is found (0 extra solutions in all 12 sectors)**; all m triples are reproduced in the 11 sectors with m ≤ 5, and 10 of the 12 triples of the m = 12 sector (the other two are never reached from 2400 random starts; the exact solve finds all 12). Together with the Jacobian rank 3 at every triple this supports (and the exact elimination proves) that the triples are isolated and simple.

## 5. G5.4: the July-2026 gl(4) Gaudin ODE

### 5.1 Transcription

`gaudinODE` is read from the cell of `GaudinSolver.wb` (JSON text) and translated mechanically (Derivative[n][Q][z] → Q_n, ^ → **, sympy `parse_expr`; the one-line normalisation of the wrapped source is the only edit): ODE = Σ_n c_n(z) Q⁽ⁿ⁾ with
c₄ = 8(z−w₁)²(z−w₂)², c₃ = −16(z−w₁)(z−w₂)(P(z−w₁) + S(z−w₂)), c₁ = 4(Cas2 − 2N(N+2))(…) etc. (all five in the output; c₀ = Cas2² + 2Cas2(N+1) − 2Cas4 + 4N(N+2)(N+3)). Cas2, Cas4, muhat, xw, Pw, degQ are transcribed by hand (the Mathematica source lines are quoted in the docstrings): xw_i = μ_i + 5/2 − i, Pw_k = Σ xw_i^k, Cas2 = Pw₂ − 5, Cas4 = Pw₄ + 4Pw₃ − Pw₁Pw₂ + Pw₂/2 − 2Pw₁² − 10Pw₁ − 51/4, muhat_k = μ_k + 4 − k, degQ_k = (S+P+3) − muhat_k. (Cas2 = Σ_i μ_i(μ_i + 5 − 2i) = C₂^{gl}(μ) exactly, since Σ_i (5/2 − i)² = 5 for N = 4 — [corrected by the session 2026-10-22 after G6's check; the original text said "C₂^{gl}(μ) − |μ|²/4", which contradicts the examples below, e.g. the singlet (1,1,1,1) has Cas2 = 4 = C₂^{gl}].) The S ↔ w₁ assignment is read off the residues: the Q‴ coefficient gives s₁-like residues 2S at w₁ and 2P at w₂.

### 5.2 (a) The notebook's claim

For (S,P) ∈ {(1,1), (2,1), (2,2), (3,2), (3,3)} (3, 3, 6, 6, 10 irreps, all multiplicity 1, G1's Racah–Speiser) and w = (0,1) (the notebook's `wsub`) as well as w = (2/7, −3/5): the dimension of the polynomial solutions of degree ≤ n (exact, n ≤ 2(S+P)+8) jumps exactly at degQ[S,P,μ] and the total is 4: **56 of 56 tests**.

### 5.3 (b), (c) Our oper = formal adjoint + gauge

Matching 8(z−w₁)²(z−w₂)² Π^p D† Π^{−p} (D† = ∂⁴ + ∂³s₁ + ∂²s₂ + ∂s₃ + s₄; the conjugation by Π^p replaces ∂ by ∂ − p₁/(z−w₁) − p₂/(z−w₂)) with gaudinODE coefficient by coefficient — polynomial identities in z, **linear** in the unknowns (v₃, v₄ᵃ, v₄ᵇ, T₂(μ), T₃(μ), T₄(μ)) — scanned over (p₁, p₂) ∈ {0..S+P+1}² and both assignments of [S,S,0,0] to the sites (w = (2/7, −3/5), top irrep): **exactly one candidate is consistent: p = (S, P) with S at w₁**: 16 candidates for (1,1) (only (1,1)), 50 for (2,1) (only (2,1), S at w₁). For all 28 irreps of the five (S,P) the system with p = (S,P) has a unique solution and **it equals G2's exact triple (v₃, v₄ᵃ, v₄ᵇ) of the sector and the central characters T_k(μ)** (28/28); with these substituted the two operators coincide as differential operators, all five coefficient functions being identical rational functions of z (exact, 28/28). Symbolic solution (`wbsym`, 150 s, S, P, w₁, w₂, C2, C4 all symbolic, p = (S,P)): with N = S + P

* **v₃ = (S+1)(2N(N+2) − Cas2)/2,   v₄ᵃ = 0,   v₄ᵇ = (Cas2² + 2(N+1)Cas2 − 2Cas4 + 4N(N+2)(N+3))/8 = T₄(μ)** (= the constant term c₀ of the ODE divided by 8),
* T₂(μ) = (4N² + 6N − Cas2)/2, T₃(μ) = (N+2)(2N(N+2) − Cas2)/2 = ((N+2)/(S+1))·v₃, T₄(μ) = v₄ᵇ,

verified against G2 on all 28 irreps (examples: (1,1), μ = (1,1,1,1) singlet: Cas2 = 4, Cas4 = 4, (v₃, v₄ᵃ, v₄ᵇ) = (12, 0, 24); (3,3), μ = (3,3,3,3): (120, 0, 360); (2,2), μ = (4,4,0,0): (0,0,0)). The family has μ₁ + μ₄ = μ₂ + μ₃ = N, so the exponents at ∞, μ̂ = (μ₁+3, μ₂+2, μ₃+1, μ₄), are symmetric about (N+3)/2: the odd characters are slaved to the even ones and only (Cas2, Cas4) appear; one of our three labels (v₄ᵃ) vanishes identically in this multiplicity-free family. The gauge power p = (S, P) is the gl(4) version of gl(3)'s p = λ₁ − λ₃ (here λ₁ − λ₄ at each site).

### 5.4 (d) Degrees

Our degrees d_j = μ_j + 4 − j ↔ their degQ_k = (S+P+3) − (μ_k + 4 − k) = d_{5−k} (symmetry of the exponent set). Moreover the adjoint solutions are built from OUR polynomial solutions: ψ_j = W₃[P_k, k ≠ j]/W₄, and Q_j = Π^pψ_j = W₃[P_k, k≠j]/((z−w₁)^S(z−w₂)^P) (W₄ ∝ (z−w₁)^{2S}(z−w₂)^{2P}) is a polynomial of degree S+P+3−d_j = degQ_j (= d_{5−j}) and solves gaudinODE exactly: **112 of 112** (4 per irrep, all 28 irreps; 12, 12, 24, 24, 40 for the five (S,P)).

## 6. Surprises

1. **The dictionary needs no fit at all**: every scalar is a central character difference (§2.1); G2's "scalars of Casimir type" are T_k(λ), T_k(λ′), τ_k, and the table of G2 §6c is exactly the solution of the residue relations.
2. **μ₄ = 0 ⟹ V₄ᵃ = V₄ᵇ = 0** (exactly), a consequence of the constant polynomial being the d₄ = 0 solution; and in the multiplicity-free family [S,S,0,0]⊗[P,P,0,0] V₄ᵃ ≡ 0 and V₄ᵇ = T₄(μ).
3. **No single condition cuts out the spectrum** in gl(4): P₄ alone is a curve in every sector; and unlike the expectation of "a few conditions suffice", all four are needed in general (8 of 68 small sectors need P₁ even given {4,3,2}); the sectors with μ₁ = μ₂ (adjacent top exponents) are the analogue of gl(3)'s E₁₂-cut.
4. **The exact solve is cheap**: the conditions have degree ≤ 9 (≤ 10 for the m ≥ 10 sectors) although d₁ ≤ 13, and the full Gröbner basis takes ≤ 1 s per sector up to m = 5 and ≈ 10 s at m = 10–12 (HW space included); the exact certification of G2's triples therefore needs no numerics at all (the numerical search is only a cross-check).
5. **The brief's negative-test criterion is mis-scaled** (§3.5): a 10⁻³ perturbation of triples of size ~50 in a pencil of dynamic range > 10¹⁰ cannot give σ_min/σ_max ≥ 10⁻⁶, and σ_min is exactly 0 whenever a lower-degree solution is independent of the perturbed component. Use the degree-exactly-d₁ residual / exact rank test.
6. The gl(3) relation V = a₁ − η/2 becomes V₃ = tS_{3,2} (no shift) with the dictionary of §2.1; the quantity that plays the role of η is τ₂ = n₁n₂ − Ω, and τ₃, τ₄ enter at the next orders.

## 7. Not done / caveats

* The "honest" single-condition loci (existence of a solution of degree d_j with the free resonance parameters kept) were not computed; the gauge-fixed (echelon) conditions of §4.3 are the exact gl(4) analogue of the gl(3) reports, and they are subsets of the honest loci.
* The subset analysis (§4.3) is restricted to the 68 sectors of eight small products and time-boxed (5 s per Gröbner basis): 5 + 1 + 2 + 2 of the 15 × 68 entries (sets containing j = 1) did not finish and are listed as T/O. The full-system solve (Σ = G2) is complete for all 165 + 15 sectors listed.
* Evidence, not a proof for all λ, λ′, μ: the exact statements hold on the sectors tested; the dictionary itself is a proof (it only uses G2's identities) and holds for all λ, λ′, μ for which the G2 identities hold (all N = 4 irreps).
* The θ₂ coefficients are not independently derived beyond G2's swap rule; they are tested (restricted operators, residue relations, Wronskian/exponent checks at θ₂).
* In the numerical search of §4.4 only 10 of the 12 solutions of the m = 12 sector are reached from 2400 random starts (no extra solution appears); the exact solve finds all 12.
* The Casimir dictionary of §5 is verified on the five (S,P) listed (28 irreps) and derived symbolically for all (S,P) by the `wbsym` stage; the notebook's `irreps[S,P]` (Schur-polynomial peeling) was replaced by G1's Racah–Speiser (identical lists, multiplicity 1 asserted).

## 8. Hand-off notes (WP-G3, WP-G4, WP-G6)

* **For G3/G4:** the joint eigenvalue triples of the slice-reduced operators must equal G2's triples (they do by construction); the oper side gives an independent characterisation of the spectrum: the points of V(J), with J generated by the gauge-fixed Frobenius conditions (`frob_conditions`, `solve_sector`); on μ₄ = 0 sectors only V₃ is non-trivial, and the characteristic polynomial of V₃ there is the gl(3)-type Q-system polynomial.
* **For the notes:** the dictionary table of §2.1 and the Casimir formulas T₂, T₃, T₄ are the gl(4) version of "a₁ = V + η/2"; the polynomial solutions have degrees μ + (3,2,1,0); Wronskian ∝ Π_s (u−θ_s)^{|λ_s|}.
* **For the July-2026 family:** (v₃, v₄ᵃ, v₄ᵇ) = ((S+1)(2N(N+2) − Cas2)/2, 0, T₄(μ)); the gauge Π^{(S,P)} and the factor 8Π² turn the adjoint of our oper into `gaudinODE`.
* Usage: `from oper_gl4 import Product, Oper, solve_sector, copy_check, dictionary`; `Oper(l1, l2, mu, t)` gives affine Laurent data, `.pencil(d)`, `.I(R, n)`; `solve_sector(l1, l2, mu, 60, subsets=False)` returns the exact solve of one sector (dict with `full`, `distinct`, `dev`, `jac_ranks`).
