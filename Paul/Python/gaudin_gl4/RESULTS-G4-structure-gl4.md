# RESULTS-G4-structure-gl4 (agent G4, 2026-10-22): structure of the slice-reduced gl(4) missing-label operators V_red,3, V_red,4a, V_red,4b on the GT lattice — stencils, closed forms, the Weyl/BGG characterisation of the physical spectrum, ordering, chain structure, polynomiality

Files (all new, in `Paul/Python/gaudin_gl4/`; nothing existing was edited, no git, no Wolfram, no notebooks; G6's files untouched):
`structure_gl4.py` (all code; CLI `python -X utf8 structure_gl4.py stencil | symbolic | symq | univ | weyl | rule | chain | poly | cycW | all` and `run <stage,...>` / `merge` which write/assemble per-stage part files),
`test_structure_gl4.py` (9 exact regression tests, ≈ 7 s: **9 tests, 0 failed**: one stencil, the universal formula and one symbolic class on held-out reps, a non-factorising class, the Weyl identity on two sectors incl. the double-string one, its sensitivity, gl(3), the Krylov polynomials, graph statistics), `structure_gl4_output.txt` (the printed output of the final runs of all stages, in the order of the tables),
`structure_gl4_tables.txt` (every table: stencils, all closed forms with their factor lists, the explicit polynomials, the Weyl tables per product, rank tables, chain and polynomiality tables) and this file.
Imports, read-only: `glN_reps`, `slice_engine_gl4` (G1), `gaudin_gl4_ops` (G2), `slice_ops_gl4` (G3: `reduced_op`, `universal_reduction`, `weight_space_mats`, `copies`, `triple_mats`); for one cross-check `gaudin_sl3/jacobi_cg.py`.
Python 3.13, numpy, scipy (LP, sparse graphs), sympy (exact factorisation, rref over Q), mpmath (50 digits); `time.perf_counter()`; at most two heavy processes at any time. Every statement below is exact (Fractions / exact rref) unless "50 digits" is quoted.

## 0. Headline results

1. **The Weyl/BGG characterisation holds exactly on every sector tested.** For every joint triple τ, Σ_{w∈S₄} sgn(w)·mult_τ(V_red(w·λ′) on V_λ[μ − w·λ′]) = [τ ∈ Spec(M_μ)], w·λ′ = w(λ′+ρ) − ρ: **832 sectors of 22 products** (Σ m = 1944 physical copies; 1586 non-empty translate spaces; 3408 distinct joint triples matched at 1e-30 from 50-digit spectra), **0 failures**. Products: the six of the brief ([2,1,0,0]², [2,1,1,0]², [3,2,1,0]², [4,2,0,0]², [4,2,1,0]², **[6,4,2,0]²** — all 76 sectors, m ≤ 12), the held-out [3,2,1,0]⊗[4,2,1,0] (both orders), [4,2,2,0]², [4,3,1,0]², and twelve more (λ ≠ λ′, λ₄ > 0, up to [8,4,0,0]² with 131 sectors and [7,4,1,0]² with 109). "All sectors" includes every dominant μ with some non-empty translate space (the irreps and formal multiplicity-0 candidates; for these products they coincide with the irreps). The test can fail: with all signs +1 it holds on only 11/16, 14/27, 32/76 sectors of [3,2,1,0]², [4,2,1,0]², [6,4,2,0]², with the plain Weyl action w(λ′) (no ρ-shift) on 1/16, 2/27, 8/76. The joint spectrum is simple on **every** translate space (also at the non-dominant weights w·λ′; closest distinct triples ≥ 4.0; every eigen-decomposition accepted at 1e-40). Not proved (§2.5 gives the argument and what is missing).
2. **One-step version and depth.** Let S(τ) = strings E_i^{k_i} violated by the eigenvector of τ ∈ Spec(W), T(τ) = {i : τ ∈ Spec of V_red(s_i·λ′) on V_λ[μ−s_i·λ′]}. Observed on all 3408 triples: **S(τ) = T(τ)** (the BGG intertwining gives ⊆; the converse is an observation) and **τ physical ⇔ S(τ) = ∅**. The union of the one-step inherited sets is **not disjoint**: of the 1464 spurious triples (w > m) 1383 violate one string, **81 violate two** — each of those lies in exactly the four translate spaces {e, s_i, s_j, w_ij} (w_ij of length 2: s₁s₃ in 38 of the 66 cases classified, one of the two products s₂s₃ / s₃s₂ in 26, one of s₁s₂ / s₂s₁ in 2) and is cancelled by inclusion–exclusion (1 − 2 + 1). **No triple needs a translate of length ≥ 3** (the longest non-empty translate space has length 2 in all 22 products). Minimal truncation depth: 329 sectors need only the identity, 438 need translates of length ≤ 1, 65 need length 2.
3. **gl(3) sanity check.** The same statement (S₃, ρ = (2,1,0), one operator V = tS₃,₂) holds exactly on all sectors of [4,2,0]², [6,3,0]², [8,4,0]², [6,4,2]², [14,7,0]², [4,2,0]⊗[3,1,0], [5,2,1]⊗[4,3,0] (0 failures; 113 sectors / 344 copies at [14,7,0]²); one-step S = T everywhere, no double-string triple, depth ≤ 1. The N = 3 family equals jacobi_cg's V_red plus the sector constant η₁/2 on all sectors of [4,2,0], [6,3,0], [6,4,2], [8,4,0] (50 digits). **Ordering in gl(3): the physical eigenvalues are the m largest of V on W when only E₂₃ is active and the m SMALLEST when only E₁₂ is active — exactly, 100 % of the sectors with w > m (including λ ≠ λ′, λ₃ > 0)**; the "m largest" of RESULTS-jacobi-cg T-statements is the E₂₃ half.
4. **Stencils.** All three operators have the same stencil of **21 displacements**: {Δm₁₂ ∈ {0, ±1}} × {row-3 move 0 or one of the 6 roots of gl(3)} = 3 × 7 (gl(3): {0, ±(1,−1)}). Structural bound from the words: 21 (V₃), 35 (V₄ᵃ, V₄ᵇ: the |Δm₁₂| = 2 classes), whose coefficients sum to **exactly zero** (a Pauli-type exclusion).
5. **Closed forms.** Symbolic GT formulas through the universal words, exact. **12 of the 20 off-diagonal classes (all with Δm₁₂ = ±1 and a non-zero row-3 move) are given by ONE product formula for all three operators.** In l-coordinates y_ik = m_ik − i, Y_s = λ_s − s, q = λ′, with a, b = the raised / lowered entry of row 2 and i, j, k = the raised / lowered / spectator entry of row 3 of the displacement:
   **M^{V₃} = (y₁₁−y_b)(y₁₁−y_a+q₁−q₂)(y_a−y_j)(y_a−y_j+1)(y_a−y_k)·Π_{s=1..4}(y_i−Y_s) / [(y_a−y_b)(y_a−y_b+1)(y_i−y_j)(y_i−y_j+1)(y_j−y_k)(y_i−y_k)],  M^{V₄ᵃ} = (y_k+3)·M^{V₃},  M^{V₄ᵇ} = −(y_a−y_i+2y_k−q₃+5)·M^{V₃}**
   (17 730 exact tests on 7 held-out reps, 16 253 non-zero entries, 0 mismatches; the only λ′-dependence is the factor (y₁₁−y_a+q₁−q₂), and q₃ in V₄ᵇ). The other 8 off-diagonal classes (6 row-3 moves with Δm₁₂ = 0; the two pure row-2 moves (±1,0,0)) and the diagonal are NOT products of linear forms: linear factors × an irreducible polynomial Q (degree 3 / 6 / 11 for V₃; the six cubics of V₃ were reconstructed explicitly, 46–50 terms, and V₄ᵃ has the SAME Q as V₃ on them). Verified against G3 on held-out reps (§3). The off-diagonal elements of V₄ᵃ are an affine multiple of those of V₃ in 18 of the 20 classes (all except (±1,0,0)), V₄ᵇ in the 12 factorising ones (§3.4).
6. **Intrinsic ordering.** No single ordering (V₃, V₄ᵃ, V₄ᵇ, G2's combination) picks Spec(M_μ): the physical set is the top-m / bottom-m of a given ordering in only 20 + 9 of 78 training sectors by V₃ (§4). Sorted by the **active strings** A = {i : E_i^{k_i} non-void on W} (intrinsic): for A = {3} and A = {1} a **linear functional ℓ_A = a V₃ + b V₄ᵃ + c V₄ᵇ exists with physical = top m of ℓ_A**, found by LP on the training products and exact on the **held-out** sectors (7/7 and 5/5); for A = {2} (the middle node, the gl(3) E₂₃ analogue) and for all mixed classes no linear functional works (LP infeasible; 21 training sectors).
7. **Chain structure.** V_red,3 is **not** Jacobi: its sparsity graph on R_ν is connected in 150/150 weight spaces (dim ≥ 2) of 9 products, is a path (tridiagonal in a suitable order) in **29 of the 44 spaces with m = w** and 19 of 106 with m < w, otherwise contains cycles (every tree that occurs is a path; max degree 8 resp. 13; reverse-Cuthill–McKee bandwidth ≤ 6 resp. 8, mean 1.84 / 3.34); in gl(3) 100 % are paths. The zero pattern is symmetric for m = w and asymmetric (cut) in 43 of 106 spaces with m < w.
8. **Polynomiality.** Some coordinate of M_μ is cyclic for the triple in **181 of 181** sectors with m ≥ 2 (cyclic coordinates 668 of 681), hence every component of every copy is a polynomial in (v₃, v₄ᵃ, v₄ᵇ) with rational coefficients (exact; reproduces G3's copies to 7e-46 at 50 digits); on the whole weight space W a cyclic coordinate vector (lattice point) exists in 150 of 150 spaces (the first GT point in 146), 87 % of all lattice points are cyclic. The minimal total degree equals the dimension-counting minimum d_gen(m) = min{d : C(d+3,3) ≥ m} in 162 sectors, d_gen+1 in 19, never smaller: components are **affine in 122 of the 138 sectors with m ≤ 4, quadratic for 5 ≤ m ≤ 10 (three sectors need 3), cubic for m = 11, 12**, versus degree m − 1 for any single operator (m = 12: degree 3 vs 11).
9. **Signed lattice count** m = Σ_w sgn(w)·dim V_λ[μ − w·λ′] holds on all 832 sectors (and all 227 gl(3) sectors); it is the Kostant/Racah–Speiser formula, so it also checks the implementation of the dot action.

## 1. Setting and conventions

Site 1 = V_λ, GT patterns ((m₁₁), (m₁₂,m₂₂), (m₁₃,m₂₃,m₃₃), λ), plain rational GT basis (G1). A weight space W = V_λ[ν], ν = μ − λ′, fixes the row sums; free lattice coordinates (m₁₂, m₁₃, m₂₃), w = dim W. V_red,k(λ′) are G3's universal polynomials in the site-1 generators with coefficients in ℚ[λ′]; for non-dominant or rational λ′ they are the same polynomials evaluated there (G3 §7: they commute on every weight space). M_μ = ∩ ker E_i^{k_i} (PRV), physical triples = joint eigenvalues on M_μ (G3's `copies`, exact PRV matrices + 50-digit eigen-decomposition). Matrices act on coordinate columns, the eigenvectors are in the plain GT basis. G2's generic combination G = V₃ + 7/19 V₄ᵃ + 11/13 V₄ᵇ is used for every joint eigen-decomposition; each operator is then conjugated to the eigenbasis and the largest off-diagonal entry checked (< 1e-40).

## 2. Q3 + Q6: the Weyl/BGG characterisation (`structure_gl4.py weyl`)

### 2.1 What was computed
For every sector μ and each of the 24 permutations p of {1,..,4}: w·λ′ = (λ′+ρ)_{p(·)} − ρ, ν_w = μ − w·λ′; if ν_w is a weight of V_λ the exact matrices of V₃, V₄ᵃ, V₄ᵇ at the weight w·λ′ on V_λ[ν_w] (G3's `reduced_op` + `weight_space_mats`), joint spectrum at 50 digits; all triples of all translate spaces clustered (tolerance 1e-30; the closest pair of distinct clusters is ≥ 4.0 apart, worst intra-cluster deviation 2.7e-48); signed count Σ sgn(w)·(multiplicity) per cluster compared with the indicator of the m physical triples. Also recorded per sector: truncated alternating sums (only translates of length ≤ L), the strings violated by each eigenvector of the identity space (exact images E_i^{k_i} of the 50-digit eigenvector), the signed lattice count.

### 2.2 Main table (all sectors of each product; failures = sectors where the signed count differs from the indicator)

| product | sectors | Σ m | translate spaces | distinct triples | signed = indicator | S = T (of triples) | triples with 2 violated strings | longest translate |
|---|---|---|---|---|---|---|---|---|
| [2,1,0,0]² | 7 | 8 | 9 | 10 | 7/7 | 10/10 | 0 | 1 |
| [2,1,1,0]² | 6 | 7 | 7 | 8 | 6/6 | 8/8 | 0 | 1 |
| [3,2,1,0]² | 16 | 30 | 22 | 36 | 16/16 | 36/36 | 0 | 1 |
| [4,2,0,0]² | 26 | 32 | 45 | 59 | 26/26 | 59/59 | 1 | 2 |
| [4,2,1,0]² | 27 | 56 | 44 | 77 | 27/27 | 77/77 | 0 | 1 |
| **[6,4,2,0]²** | 76 | 297 | 143 | 424 | 76/76 | 424/424 | 3 | 2 |
| [3,2,1,0]⊗[4,2,1,0] (held-out) | 20 | 38 | 24 | 42 | 20/20 | 42/42 | 0 | 1 |
| [4,2,1,0]⊗[3,2,1,0] | 20 | 38 | 38 | 62 | 20/20 | 62/62 | 0 | 1 |
| [4,2,2,0]² (held-out) | 23 | 26 | 37 | 42 | 23/23 | 42/42 | 0 | 1 |
| [4,3,1,0]² (held-out) | 31 | 69 | 52 | 94 | 31/31 | 94/94 | 2 | 2 |
| [3,2,0,0]⊗[2,1,1,0] | 10 | 8 | 19 | 19 | 10/10 | 19/19 | 0 | 1 |
| [2,1,1,0]⊗[2,2,1,0] | 6 | 7 | 7 | 8 | 6/6 | 8/8 | 0 | 1 |
| [5,3,1,0]² | 49 | 136 | 88 | 198 | 49/49 | 198/198 | 2 | 2 |
| [5,3,2,0]² | 45 | 121 | 79 | 171 | 45/45 | 171/171 | 0 | 1 |
| [6,3,0,0]² | 65 | 90 | 132 | 217 | 65/65 | 217/217 | 9 | 2 |
| [5,2,1,1]⊗[3,2,2,1] | 12 | 8 | 27 | 23 | 12/12 | 23/23 | 3 | 2 |
| [5,4,2,1]² | 31 | 69 | 52 | 94 | 31/31 | 94/94 | 2 | 2 |
| [6,3,2,1]² | 41 | 82 | 73 | 131 | 41/41 | 131/131 | 1 | 2 |
| [8,4,0,0]² | 131 | 205 | 292 | 614 | 131/131 | 614/614 | 36 | 2 |
| [4,4,2,0]⊗[4,3,2,1] | 21 | 27 | 41 | 51 | 21/21 | 51/51 | 1 | 2 |
| [7,4,1,0]² | 109 | 396 | 222 | 681 | 109/109 | 681/681 | 15 | 2 |
| [6,4,2,0]⊗[5,3,1,0] | 60 | 194 | 133 | 347 | 60/60 | 347/347 | 6 | 2 |
| **TOTAL** | **832** | **1944** | **1586** | **3408** | **832/832** | **3408/3408** | **81** | |

Every physical triple occurs in the identity space only (signed count +1); every spurious triple of W occurs in the identity space and in the translates s_i·λ′ of the violated strings (pattern of lengths (0,1): 1383 triples) or additionally in w_ij (pattern (0,1,1,2): 81 triples); the signed count of each spurious triple is 0 in all 1464 cases. The per-sector listings (m, w, translate dimensions by length, minimal depth) for the first five products are in `structure_gl4_tables.txt`. Sensitivity (sectors where the identity survives a deliberately wrong rule): [3,2,1,0]²: 16/16 (conjecture), 11/16 (all signs +1), 1/16 (plain w(λ′)); [4,2,1,0]²: 27/27, 14/27, 2/27; [6,4,2,0]²: 76/76, 32/76, 8/76.

### 2.3 (ii) the one-step version; (iii) gl(3)
See headline 2 and 3. Quantitatively (all 22 gl(4) products): S ⊆ T and S = T on 3408/3408 triples; "physical ⇔ no string violated" on 3408/3408 (this also confirms that the eigenvectors in Spec(M_μ) are exactly those satisfying the PRV conditions); of 1464 spurious triples 1383 have |S| = 1 and 81 have |S| = 2 (never 3). In gl(3) (7 products, 227 sectors): signed count = indicator everywhere, S = T everywhere, |S| ≤ 1 always, depth ≤ 1.

### 2.4 Q6: the signed lattice count
m = Σ_w sgn(w) |R_{μ−w·λ′}| (|R_ν| = dim V_λ[ν]) equals the Racah–Speiser multiplicity on 832/832 sectors (and on all gl(3) sectors); for the formal sectors it is 0. It is the Kostant–Racah–Speiser formula with GT weight multiplicities, so it is a consistency check of the dot action and of the weight-space dimensions, not new information.

### 2.5 What is verified and what is argued
*Verified (exact/50 digits):* the identity above on 832 gl(4) + 227 gl(3) sectors; the BGG intertwining E_i^{k_i} V(λ′) = V(s_i·λ′) E_i^{k_i} (G3, tested 5814 times) is what makes S ⊆ T. *Argued, not proved:* V_red(w·λ′) on V_λ[μ−w·λ′] is the Gaudin family on the singular vectors of V_λ ⊗ M(w·λ′) (G3's interpretation). The Gaudin algebra commutes with gl(4), so the decomposition by joint eigenvalue τ is a decomposition of g-modules and an exact functor; applied to V_λ ⊗ (BGG resolution of L(λ′)) it gives an exact resolution of (V_λ ⊗ L(λ′))_τ by (V_λ ⊗ M(w·λ′))_τ with signs (−1)^{ℓ(w)}; the identity is then the statement that taking singular vectors of weight μ (Hom(M(μ), ·), only left exact) preserves the Euler characteristic τ by τ, as it does for the dimensions (Kostant's formula, §2.4). That last step is what the data confirm and what would have to be proved (e.g. via the Verma flag and the vanishing of Ext¹ between the relevant Verma modules for dominant μ). The table also shows how little depth is needed: only translates of length ≤ 2 ever contribute for λ, λ′ ≤ 8, and the double-string spurious triples occur in exactly the Bruhat interval [e, w_ij] = {e, s_i, s_j, w_ij}.

## 3. Q1 + Q2: stencils and closed forms (`stencil`, `symbolic`, `symq`)

### 3.1 Stencils (Q1)
Generator moves as GLNRep builds them: E_{k,k+1} raises ONE entry of row k, E_{k+1,k} lowers one, E_ii does not move, E_ij (|i−j| > 1) = commutator of two factors moves one entry in each row min..max−1. The weight-zero part of the Minkowski sum over each universal word gives the structural stencil; data = nonzero entries of G3's exact matrices on **every** weight space of V_λ for λ ∈ {[2,1,0,0], [3,2,1,0], [4,2,1,0], [4,3,1,0], [4,2,2,0], [5,3,2,1], [6,4,2,0]}, with 3 generic rational non-dominant λ′ (e.g. (1/3, 5/7, −2, 4/5)) and with λ′ = λ.

| operator | reduced words | structural stencil | data stencil, generic λ′ | data, dominant λ′ = λ (union) | symmetric under Δ → −Δ |
|---|---|---|---|---|---|
| V₃ | 39 | 21 | **21** | 21 | yes |
| V₄ᵃ | 78 | 35 | **21** | 21 | yes |
| V₄ᵇ | 85 | 35 | **21** | 21 | yes |

The 21 displacements are (Δm₁₂, Δm₁₃, Δm₂₃) ∈ {−1,0,1} × ({(0,0)} ∪ {(±1,0), (0,±1), ±(1,−1)}), i.e. "row 2 moves by 0 or one root of gl(2), row 3 by 0 or one root of gl(3)" (|Δ|₁ ≤ 3, 20 off-diagonal classes); for dominant λ′ individual reps can have fewer (V₄ᵃ ≡ 0 on λ₃ = 0 reps; [4,2,2,0]: 9). The structural bound of V₄ᵃ, V₄ᵇ contains the 14 extra classes with Δm₁₂ = ±2; **their coefficient is exactly 0 after summation** (symbolic: all 14 classes sum to the zero rational function for both operators) — a Pauli-type exclusion (an entry of row 2 is never raised twice). Every nonzero entry of G3's matrices on the tested spaces lies in the stencil (0 uncovered).

### 3.2 How the closed forms were derived (Q2)
The generic GT formulas (l-coordinates y_ik = m_ik − i, y_j4 = λ_j − j; raising: −Π_j (y_ik − y_{j,k+1}) / Π_{j≠i}(y_ik − y_jk), lowering with the row below; E_ij by commutators) are applied to a generic pattern with symbolic labels through G3's universal reduced words, **all intermediate states kept** (no boundary treatment), λ′-monomials as single-variable forms. A matrix element is then a sum of 148 – 3384 terms c·Π(affine forms)^e per displacement class (5052 / 25342 / 26118 terms in all for V₃ / V₄ᵃ / V₄ᵇ); naive expansion is hopeless (sympy fields: > 2 min for one word of V₃), so the sums are simplified by **exact line reconstruction**: N = M·D (D = maximal denominator) is a polynomial in a line parameter, interpolated exactly and factored over ℚ; the roots on the 14 coordinate axes are clustered into affine forms (a form with coefficient vector c and value F at the base point has the roots −F/c_i), kept only if they are genuine factors on generic lines; what is left is the non-linear part Q. A class is "fully factorising" if nothing is left. The identification is checked on independent random lines (checks True for every class; a failed random draw is retried with fresh points — one class of V₄ᵇ needed it). Forms are printed in the GT labels (m_ik, l_j = λ_j, q_j = λ′_j).

### 3.3 Results
| operator | classes with terms | identically 0 after summation | fully factorising (product of linear forms) | with a non-linear part Q |
|---|---|---|---|---|
| V₃ | 21 | 0 | **12** | 9 (6 row-3 moves with Δm₁₂ = 0: Q cubic; 2 pure row-2 moves (±1,0,0): degree 6; diagonal: degree 11) |
| V₄ᵃ | 35 | 14 (|Δm₁₂| = 2) | **12** | 9 (cubic ×6; degree 7 ×2; diagonal 12) |
| V₄ᵇ | 35 | 14 | **12** | 9 (quartic ×6; degree 7 ×2; diagonal 12) |

*Fully factorising classes (Δm₁₂ = ±1, row-3 move ≠ 0).* Their 36 factor lists (`structure_gl4_tables.txt`, printed in the labels m_ik, l_j) become uniform in the l-coordinates; the identification gave 9 (V₃) / 10 (V₄ᵃ) / 10–15 (V₄ᵇ) linear numerator factors over 6–9 denominator factors, with κ = ±1, and in the orientation chosen below the sign is +1 for all 12 classes of V₃. **The universal formula** (stage `univ`; a, b = raised / lowered row-2 entry, i, j, k = raised / lowered / spectator row-3 entry, y_ik = m_ik − i, Y_s = λ_s − s):
M^{V₃} = (y₁₁−y_b)(y₁₁−y_a+q₁−q₂)·(y_a−y_j)(y_a−y_j+1)(y_a−y_k)·Π_{s=1..4}(y_i−Y_s) / [(y_a−y_b)(y_a−y_b+1)·(y_i−y_j)(y_i−y_j+1)(y_j−y_k)(y_i−y_k)],
M^{V₄ᵃ} = (y_k+3)·M^{V₃},   M^{V₄ᵇ} = −(y_a−y_i+2y_k−q₃+5)·M^{V₃}.
E.g. Δ = (+1,+1,0) (a = m₁₂, b = m₂₂, i = m₁₃, j = m₃₃, k = m₂₃) in the labels: M^{V₃} = −(m₁₃−l₁)(m₁₃−l₂+1)(m₁₃−l₃+2)(m₁₃−l₄+3)(m₁₂−m₂₃+1)(m₁₂−m₃₃+2)(m₁₂−m₃₃+3)(m₁₁−m₂₂+1)(m₁₁−m₁₂+q₁−q₂) / [(m₁₂−m₂₂+1)(m₁₂−m₂₂+2)(m₁₃−m₂₃+1)(m₁₃−m₃₃+2)(m₂₃−m₃₃+1)(m₁₃−m₃₃+3)] (the sign −1 in the labels is the orientation of the factors, e.g. (y_j − y_k) = −(m₂₃ − m₃₃ + 1)…).
Structure: the raised row-3 entry carries the four GT factors Π_s(y_i − Y_s); the lowered one the pair (y_a − y_j)(y_a − y_j + 1) with the raised row-2 entry; the row-1/row-2 factors (y₁₁−y_b)(y₁₁−y_a+q₁−q₂)/[(y_a−y_b)(y_a−y_b+1)] are those of the gl(3) chain u, d of jacobi_cg with the twist q₁ − q₂ (with m₁₁ = μ₁ − λ′₁ these read (μ₁ − λ′₂ − y_a …)). The three operators therefore have the same off-diagonal structure on these classes and differ by the factors 1, (y_k + 3), −(y_a − y_i + 2y_k − q₃ + 5).

*Not products of linear forms.* The 6 row-3 moves with Δm₁₂ = 0 have 6 – 7 linear numerator factors and a **cubic** Q (V₃, V₄ᵃ; support m₁₁, m₁₂, m₂₂ and the two moved row-3 entries and q₁, q₂, q₃) resp. quartic Q (V₄ᵇ); the pure row-2 moves (±1,0,0) a Q of degree 6 (V₃) / 7 (V₄) in (m₁₃, m₂₃, m₃₃, l₁..l₄) (+ q₃, q₄ for V₄); the diagonals degree 11 / 12 in all labels. The cubic Q of V₃ were **reconstructed explicitly** (exact interpolation in the 8 support variables, 119 monomials, exact rref over ℚ): 46 – 50 terms each, irreducible, e.g. for Δ = (0,−1,0) (printed in full in the tables) Q = m₁₁²m₁₂ + m₁₁²m₂₂ − 2m₁₁²m₂₃ − m₁₁m₁₂² − 2m₁₁m₁₂m₂₂ + … + 4m₁₂ − … − 4m₂₂ (47 terms). For V₄ᵃ the (0,−1,0) polynomial is the same 47-term Q as for V₃, for V₄ᵇ it has 163 terms (degree 4). The row-2 moves and the diagonal were not reconstructed explicitly (too many monomials for exact interpolation within the time box); they are verified as "linear factors × degree-qdeg polynomial" (below).

### 3.4 Held-out verification against G3 (exact)
Protocol: for every weight space W of V_λ, every λ′ in the list, every source pattern x ∈ W and every class Δ of the stencil with a **valid target pattern** (the formal coefficient need not vanish when the target is no GT pattern — the finite-dimensional module is a quotient of the formal one — so only valid targets are tested), closed form(x) = G3's matrix entry; additionally every nonzero entry of G3's matrices must belong to a class (completeness). Removable poles at boundary patterns are handled by the exact Laurent series along a line (no zero-denominator skips remained).
* factorising classes, V₃/V₄ᵃ/V₄ᵇ, λ ∈ {[3,2,1,0], [4,3,1,0], [5,3,2,1], [6,4,2,0]} (none used in the derivation: the derivation is symbolic), λ′ = λ and two generic rational non-dominant λ′: **15 546 tests, 12 471 nonzero entries, 0 mismatches, 0 pole skips.**
* all classes incl. the non-factorising ones and the diagonal (Q evaluated as a degree-qdeg polynomial on a line through the point, degree checked by an extra point), [3,2,1,0], λ′ = λ and generic: V₃ 152 + 152 tests, V₄ᵃ 152 + 152, V₄ᵇ 152 + 152, **0 mismatches, 0 uncovered nonzero entries**.
* explicit Q of the six V₃ row-3 classes (and the example class of V₄ᵃ, V₄ᵇ), evaluated from the printed polynomial (independent of the term sums): [3,2,1,0], [4,3,1,0], [5,3,2,1], λ′ = λ and generic: 40 + 40, 164 + 164, 120 + 120 tests for V₃ (V₄ᵃ, V₄ᵇ: 6–20 tests per run), **0 mismatches**.

*Relations between the operators.* On the 12 factorising classes M^{V₄ᵃ} = (y_k + 3) M^{V₃} and M^{V₄ᵇ} = −(y_a − y_i + 2y_k − q₃ + 5) M^{V₃} (above, 17 730 tests). Beyond them (observed exactly along three random lines per class): the ratio M^{V₄ᵃ}/M^{V₃} is an **affine function of the labels in 18 of the 20 off-diagonal classes** (all except the pure row-2 moves (±1,0,0)); for the six Δm₁₂ = 0 classes V₄ᵃ and V₃ even have the same polynomial Q (checked for (0,−1,0): identical 47-term polynomials). M^{V₄ᵇ}/M^{V₃} is not affine for the six Δm₁₂ = 0 classes and the row-2 moves (V₄ᵇ's Q is a different quartic, 163 terms for (0,−1,0)).


## 4. Q4: which eigenlines are physical? (`rule`)

Sectors with w > m: 78 in the training products ([2,1,0,0]², [2,1,1,0]², [3,2,1,0]², [4,2,0,0]², [4,2,1,0]², [6,4,2,0]²), 28 in the held-out ones ([3,2,1,0]⊗[4,2,1,0], [4,2,2,0]², [4,3,1,0]²). The rank (1 = largest) of every physical triple in the four orderings is tabulated per sector in `structure_gl4_tables.txt`; example [3,2,1,0]² (ranks of the physical triples in V₃ | V₄ᵃ | V₄ᵇ | G, w = 4, m = 3 or 2):

| μ | m | active strings | ranks V₃ | V₄ᵃ | V₄ᵇ | G | directions |
|---|---|---|---|---|---|---|---|
| (5,3,3,1) | 3 | {2} | 1,2,4 | 1,3,4 | 1,2,4 | 1,3,4 | none/none/none/none |
| (5,3,2,2) | 3 | {3} | 1,2,3 | 1,2,3 | 1,3,4 | 1,2,3 | top/top/none/top |
| (4,4,3,1) | 3 | {1} | 2,3,4 | 2,3,4 | 1,2,3 | 1,2,3 | bottom/bottom/top/top |
| (4,4,2,2) | 2 | {1,3} | 2,2 | 2,3 | 1,3 | 1,4 | none |
| (4,3,3,2) | 3 | {2} | 1,2,4 | 1,2,4 | 2,3,4 | 2,3,4 | none/none/bottom/bottom |

**No single ordering works.** Training, per ordering (top = physical set is exactly the top m; bottom likewise; none = interleaved): V₃: top 20, bottom 9, none 49; V₄ᵃ: top 11, bottom 6, tie 20, none 41; V₄ᵇ: top 18, bottom 4, tie 6, none 50; G2's combination: top 24, bottom 7, none 47 (of 78). An LP over all real functionals a V₃ + b V₄ᵃ + c V₄ᵇ finds none that works on all 78 training sectors (margin 0).

**Classes by active strings A (intrinsic: A = {i : ν + k_iα_i is a weight of V_λ}); LP for a linear functional with physical = top m on every training sector of the class, then tested on held-out sectors:**

| class A | training sectors | linear functional separates? | functional l = (a,b,c) for (V₃, V₄ᵃ, V₄ᵇ), max-margin | held-out sectors picked correctly |
|---|---|---|---|---|
| {1} | 10 | **yes**, margin 0.023 | (−1.000, +0.269, +0.103) | **5/5** |
| {3} | 22 | **yes**, margin 0.144 | (−0.487, +1.000, +0.329) | **7/7** |
| {2} | 21 | **no** (LP: no functional with positive margin) | – | – (a quadratic functional with 9 monomials separates the training sectors, margin 0.003, and picks 10 of the 11 held-out sectors) |
| {1,2} | 5 | no | – | – |
| {1,3} | 9 | no (quadratic: no) | – | – |
| {2,3} | 9 | no (quadratic: no) | – | – |
| {1,2,3} | 2 | no | – | – |
| all sectors | 78 | no (quadratic: no) | – | – |

So for the two **end-node** strings the physical eigenvalues are separated from the spurious ones by a **linear** functional of the triple (in gl(3) the E₁₂ condition gives 'bottom of V', the E₂₃ condition 'top of V'; in gl(4) the two end nodes are separable with different functionals, ℓ₁ ≈ −V₃ + 0.27 V₄ᵃ + 0.10 V₄ᵇ and ℓ₃ ≈ V₄ᵃ − 0.49 V₃ + 0.33 V₄ᵇ, i.e. 'bottom' by V₃ for node 1 and 'top' by a V₄ᵃ-dominated combination for node 3); for the **middle** string and every mixed class there is no linear rule (the quadratic fit for {2} is a 9-parameter fit to 149 pairs with margin 0.003 and should not be trusted). The functionals found are maximum-margin LP solutions, not canonical: the cone of valid functionals was not computed. A pure "sign rule" (sign structure of the matrix, Perron–Frobenius) is not available: the off-diagonal elements of V₃ in the plain GT basis have both signs within the same displacement class (11 of 20 classes mixed for [6,4,2,0]).
**gl(3), for comparison:** the same tabulation gives a two-sided but exact rule: A = {2} (E₂₃ active): physical = top m of V; A = {1} (E₁₂ active): physical = bottom m (7 products, all 70 sectors with w > m, 0 exceptions; e.g. [14,7,0]²: 21 top + 21 bottom).

## 5. Q5: chain structure and polynomiality (`chain`, `poly`)

### 5.1 Sparsity graph (a)
Over the weight spaces (dim ≥ 2) of 9 products ([2,1,0,0]², [2,1,1,0]², [3,2,1,0]², [4,2,0,0]², [4,2,1,0]², [6,4,2,0]², [3,2,1,0]⊗[4,2,1,0], [4,2,2,0]², [4,3,1,0]²): 150 spaces, of which 44 with m = w:

| operator | kind | spaces | connected | paths (= tridiagonal in some order) | trees | mean RCM bandwidth | max bandwidth | max degree | max w |
|---|---|---|---|---|---|---|---|---|---|
| V₃ | m = w | 44 | 44 | **29** | 29 | 1.84 | 6 | 8 | 10 |
| V₃ | m < w | 106 | 106 | 19 | 19 | 3.34 | 8 | 13 | 15 |
| V₄ᵇ alone | m = w | 44 | 27 | 29 | 29 | 1.34 | 6 | 8 | 10 |
| V₄ᵇ alone | m < w | 106 | 100 | 19 | 19 | 3.24 | 8 | 13 | 15 |

The graphs of V₃, V₄ᵃ and V₄ᵇ together have exactly the graph of V₃ (V₃'s nonzero pattern contains the others). Every tree that occurs is a path; the other graphs contain cycles. So V_red,3 is tridiagonal in an appropriate ordering in 29 of the 44 spaces with m = w (in the remaining 15 the graph contains cycles; the smallest example is (5,4,2,1) of [3,2,1,0]², m = w = 4, where the graph is the complete graph K₄: 6 edges, bandwidth 3), but **it is not a chain in general**: the lattice is genuinely 3-dimensional (21-point stencil). The zero pattern is symmetric for all m = w spaces and asymmetric in 43 of 106 spaces with m < w (the one-sided cut: the analogue of the E₁₂ cut of gl(3), where u = 0 at the cut). gl(3) (the same code, [8,4,0]², [6,4,2]²): all 38 spaces are paths (bandwidth 1). V₄ᵇ alone often disconnects the graph (27/44 connected for m = w): its zero pattern is sparser.

### 5.2 Polynomiality and cyclic coordinate vectors (b, c)
On M_μ (exact m×m matrices in the PRV coordinates, m ≥ 2). A coordinate x₀ is cyclic iff every copy has a nonzero component there; then u_j[y]/u_j[x₀] = P_y(v₃, v₄ᵃ, v₄ᵇ) with P_y ∈ ℚ[v] **exactly** (left Krylov space of e_x₀ᵀ under the commutative algebra of the triple; the graded minimal total degree deg(y) and the coefficients are exact; the polynomials are then evaluated at the 50-digit triples of G3's copies and compared with their coordinates: worst relative error 7.4e-46). Products: the six training products, the held-out ones, [5,3,1,0]², [5,3,2,0]².

* **181 sectors with m ≥ 2** (Σ m = 681): a cyclic coordinate exists in **181**, all m coordinates are cyclic in 171; coordinates cyclic for the triple ⇔ for the generic single operator G in every sector.
* Degree of the best cyclic coordinate (max over y of deg(y)) by multiplicity m, compared with d_gen(m) = min{d : C(d+3,3) ≥ m} (polynomials of degree ≤ d in three variables can interpolate any function on m generic points): m = 2: d = 1 (57 sectors); m = 3: 1 (44), 2 (6); m = 4: 1 (21), 2 (10); m = 5: 2 (9); m = 6: 2 (16); m = 7: 2 (6), 3 (1); m = 8: 2 (4); m = 9: 2 (2); m = 10: 2 (1), 3 (2); m = 11: 3; m = 12: 3. **Equal to d_gen in 162 sectors, d_gen + 1 in 19, never smaller** (in all 19 — m = 3, 4, 7, 10, μ₄ > 0 in 16 of them — the degree-d_gen evaluation matrix on the m joint triples has rank < m, checked numerically: three collinear triples (m = 3), four coplanar (m = 4), 7 or 10 triples on a quadric). For a single generic operator the Krylov degree is m − 1 (mean over [6,4,2,0]²: 3.95 vs 1.62 for the triple; m = 12: 11 vs 3).
* **Interpretation:** the answer to "is the k-th slice component a polynomial in the triple?" is **yes, with rational coefficients, and the degree is the minimum that dimension counting allows** — grows like (6m)^{1/3}, not like m − 1 as in the gl(3) chain (where slice_k = p_{k−1}(V)). The three-parameter lattice gives no further degree reduction (the components are "generic" polynomials of the triple of that degree): the multivariate orthogonal-polynomial structure is that of m generic points. In particular for m ≤ 4 the components are affine functions of (v₃, v₄ᵃ, v₄ᵇ) in 122 of the 138 sectors (a consequence of counting for non-degenerate point sets, not a structural coincidence; the 16 others have collinear / coplanar triples; for m ≥ 5 affine is impossible).
* **(c) Coordinate vectors of the lattice R_ν on the whole weight space W** (`cycW`; exact graded Krylov spaces of e_xᵀ under the triple acting on W, spurious eigenvectors included): **a cyclic lattice point exists in 150 of 150 weight spaces** (dim ≥ 2) of the 9 products, in particular in all 106 spaces with m < w; **620 of 709 lattice points (87.4 %) are cyclic** (all points cyclic in 101 spaces); **the first lattice point in GT order is cyclic in 146 of the 150 spaces** (the last one in 109); the best graded degree equals d_gen(w) in 141 spaces and exceeds it in 9 (7 in [4,2,0,0]², where V₄ᵃ ≡ 0 and V₄ᵇ is a function of V₃). So the copies of M_μ ⊂ W are, on the lattice, Lanczos-type vectors: every eigenvector of W (physical or not) is a polynomial in its triple applied to one coordinate component, with the degree growth of §5.2.

## 6. Surprises and corrections

0. **One product formula covers 12 of the 20 off-diagonal classes of all three operators** (§3.3): V₄ᵃ and V₄ᵇ are V₃ times a linear factor ((y_k+3), resp. −(y_a−y_i+2y_k−q₃+5)); V₃'s own factor structure is the gl(3) chain factor × the four GT factors of the raised row-3 entry × simple row-3 couplings — it was found by converting the machine-identified factor lists to l-coordinates and then proved exact only by 17 730 held-out tests (not derived).
1. The Weyl identity is *exact*, needs no tolerance beyond the 50-digit linear algebra, and the translates carry exactly the spurious part: the physical spectrum is the Euler characteristic of the BGG complex of Gaudin spectra. Depth is tiny (≤ 2).
2. The double-string triples are the Bruhat-interval patterns {e, s_i, s_j, w_ij}: the one-step union is not disjoint, but a single further translate repairs it.
3. V₄ᵃ and (partly) V₄ᵇ off-diagonal elements are *affine multiples* of V₃'s; the stencil of all three is the same 21 although the words of V₄ reach |Δm₁₂| = 2 (exact cancellation).
4. The λ′-dependence of the 12 factorising classes sits in ONE linear factor (q₁ − q₂ in V₃, V₄ᵃ; V₄ᵇ has in addition one factor with q₃); q₃ also enters the non-linear parts Q of the V₃ row-3 classes.
5. For gl(3) the physical eigenvalues are the m largest of V when only E₂₃ is active **and the m smallest when only E₁₂ is active** (RESULTS-jacobi-cg §4 states the first half only).
6. In gl(4) the end-node classes are linearly separable, the middle node is not.
7. Brief inaccuracies/corrections: "the off-diagonal elements should factor into products of linear factors" holds for 12 of the 20 classes only; "is M_μ cut out by coordinate conditions" — only the E₁₂-type string is a coordinate cut in the GT basis (the asymmetric zero patterns for m < w), the others are spectral conditions; the explicit minimum degree for polynomiality is not m − 1 but the counting bound; physical-triple tables of G3 §8 (V₃-degeneracy) are consistent with this report (simple joint spectra on every translate space).
8. Tooling trap (as in G1–G3): Python source written through Bash heredocs loses `\n` escapes; all code was written with the file tools or checked by `ast.parse` (twice repaired).
9. A bug of my own caught by the tests: comparing a 60-digit mpf with a rounded sum at default precision (`v > x + 1e-20` at 15 digits) counted a value as larger than itself; fixed (`v − x > tol`), all rank tables are from the fixed code.

## 7. Not done / limits

* No proof of the Weyl identity (§2.5), of S = T, of the exclusion of |Δm₁₂| = 2 (exact cancellation verified symbolically for V₄ᵃ, V₄ᵇ only as the vanishing of the summed rational function), or of the closed forms (they are derived symbolically and verified exactly, not proved by induction; the reconstruction of linear factors is exact but its "genuine factor" test uses generic lines).
* Explicit closed forms: the 12 factorising classes of V₃, V₄ᵃ, V₄ᵇ and the six cubic Q of V₃ (+ one class each for V₄ᵃ, V₄ᵇ). **Not explicit:** the cubic/quartic Q of the other V₄ classes, the Q of the pure row-2 moves (degree 6/7) and the diagonals (degree 11/12, a polynomial plus poles); these are given as "linear factors × degree-d polynomial" and verified through the exact term-sum / line extrapolation. The diagonal has removable poles at m₁₂ = m₂₂ etc. (the gl(3) convention K/ξ := 0 generalises).
* The weight-space/GT lattice structure for the largest products ([9,6,3,0]², m ≤ 27, w ≤ 32) was not run in the Weyl test (cost: 50-digit eigen-decomposition of w×w matrices on up to 24 translate spaces per sector); [8,4,0,0]², [7,4,1,0]², [6,4,2,0]² ⊗ [5,3,1,0] and [6,4,2,0]² are the largest.
* Ordering rules: the cone of separating functionals (not only one LP point) and a rule for the middle string and the mixed classes were not found; a rule that is independent of the (derived) class S was found only through the intrinsic class A for the end nodes.
* Only N = 3, 4. λ₄ > 0 only in the listed products.

## 8. Usage

```python
import sys; sys.path.insert(0, 'Paul/Python/gaudin_gl4')
import structure_gl4 as sg
r = sg.weyl_sector((4,2,0,0), (4,2,0,0), (5,5,1,1))     # all 24 translates of one sector: r['groups'] (triples, lens, signed, phys), r['fails'], r['depth_ok'], r['onestep'] ((S, T) per eigenvector), r['mult_signed']
res, S = sg.weyl_product((6,4,2,0), (6,4,2,0), kind='all')     # all sectors of a product (fam=sg.FAM[3] for gl(3)); sg.sector_list(lam, lamp, 'rs'|'all')
sg.dot_action(p, lamp)                                  # w.lam' for a permutation tuple p
sg.structural_stencil('V3'); sg.data_stencil(lam, [lamp])
rec = sg.closed_form_one('V3', d)                       # d = 6-vector of entry changes (sg.project(d) = (dm12, dm13, dm23)); rec['net'] = factor list, rec['kappa'], rec['qdeg'], sg.reconstruct_q(rec) -> rec['Qstr']
sg.verify_closed('V3', lam, [lamp], random.Random(1), which='F', CF={d: rec})     # exact comparison with G3's matrices
sg.universal_F('V4a', pt, d)                            # the product formula of the 12 factorising classes (pt = sg.pt_of_pattern(pattern, lam, lamp), d in sg.f_classes()); sg.verify_universal(lam, [lamp])
sec = sg.ident_sector(lam, lamp, mu)                    # whole weight space: matrices, 50-digit joint spectrum, physical flags, violated strings
sg.left_krylov(mats, x0)                                # exact graded Krylov degrees / polynomials of a coordinate vector
sg.graph_stats(M)                                       # sparsity graph, path test, RCM bandwidth
sg.stage_cycW(); sg.stage_rule(); sg.stage_poly()       # the tables of Q5(c), Q4, Q5(b)
```
CLI: `python -X utf8 structure_gl4.py <stage>`; each stage prints its tables (also written to `structure_gl4_tables.txt` by `run`/`merge`). Run times (this machine, from the output): stencil 7 s, weyl 45 s (22 gl(4) products + 7 gl(3) products + sensitivity), rule 5 s, chain 1.4 s, poly 4.9 s, cycW 2.9 s, univ 7 s, symq 294 s, symbolic 1061 s (V₃ 34 s, V₄ᵃ 230 s, V₄ᵇ 241 s for the line reconstruction; the exact verification against G3 takes the rest).
