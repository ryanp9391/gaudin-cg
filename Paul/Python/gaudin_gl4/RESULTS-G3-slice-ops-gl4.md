# RESULTS-G3-slice-ops-gl4 (agent G3, 2026-10-22): the two-site gl(4) Gaudin missing-label operators V₃, V₄ᵃ, V₄ᵇ reduced to the PRV slice space, and the canonical copies

Files (all new, in `Paul/Python/gaudin_gl4/` unless stated; nothing existing was edited, no git, no Wolfram, no notebooks; G5's files untouched):
`slice_ops_gl4.py` (all code; CLI `python -X utf8 slice_ops_gl4.py order | reduce | gl3 | direct | copies | bench [ext] | ref | verify | all [ext]`),
`test_slice_ops_gl4.py` (9 exact regression tests, ≈ 1.5 s: **9 tests, 0 failed**),
`slice_ops_gl4_output.txt` (the printed output of ONE run of `all ext`, total 476.6s; every table below is parsed from it),
`../../Mathematica/Gaudin/ref/gl4-cg-copies-{2100,2110,3210,2210-2110,4200,6420}.txt` (the six reference files of G3.6, §10) and this file.
Python 3.13.13, numpy 2.5.2, scipy 1.18.1, mpmath 1.3, no gmpy2; i7-12700KF, 32 GB, Windows 11; single process (≤ 2 heavy processes of mine at any time; the machine was shared with agent G5, so timings carry ≈ 5–10 % noise; they are min over cold runs).
Imports, read-only: `glN_reps` (`GLNRep`, `racah_speiser`, …) and `slice_engine_gl4` (`prv_space`, `hw_from_slice`, `hw_from_slice_exact`, `hw_residual`, `Recursion`, `parse_ref`) of G1; `gaudin_gl4_ops` of G2 (`load_terms`, `Op2`, `Rep`, `HWExact`, `restrict_exact`, `joint_spectrum`, `laurent`, `antisym_form`, `pbw_word`, …) and its term-list file `ref/gl4-gaudin-ops.txt`; `gaudin_sl3/slice_engine.py` (the gl(3) engine) for the N = 3 convention test.

## 0. Headline results

1. **The reduction rule works, with the sign/order rule as stated.** N = 3: the term list of t·S_{3,2} (G2's symbolic engine, N = 3) reduced by the new engine equals `slice_engine.vred_apply` + (η₁/2)·1 **exactly on all 35 PRV basis vectors** of every irrep of [2,1,0]², [3,1,0]², [4,2,0]²; S_{3,1} (whose words have two site-2 factors) satisfies t²·S_{3,1,red} = −2·(reduced tS_{3,2}) + c·1 with an integer c on all 26 sectors, and the full two-site S_{3,1}, S_{3,2} on the exact HW vectors confirm it (70/70, slice and full vector). gl(4): for V₃, V₄ᵃ, V₄ᵇ **and** S_{3,1}, S_{4,1} (up to three site-2 factors) on all 46 irreps (64 PRV vectors) of the six products of G3.3 the slice read off the FULL two-site operator applied exactly to the exact HW vectors equals the slice-reduced operator **exactly, and the whole two-site vector O Ψ_j equals Σ_i A_ij Ψ_i exactly** (320 + 320 checks, 0 mismatches); at scale 598 further PRV vectors × 3 operators (11 more products incl. λ ≠ λ′ in both orders, up to [6,4,2,0]² with 297 vectors) exact, 0 mismatches; float route ≤ 1.7e-14; the θ₂ Laurent coefficients (sites exchanged; 2275 slice checks) and G2's antisymmetrised representative of V₄ᵇ (449 vectors, a different 1464-term polynomial) agree, §5.
2. **The reduced operators are universal polynomials in the site-1 generators with coefficients in ℚ[λ′]**: V₃, V₄ᵃ, V₄ᵇ, S_{3,1}, S_{4,1} reduce to **39, 78, 85, 40, 86 words** (88 / 152 / 324 / 156 / 412 original terms; at λ′ = [3,2,1,0]: 39 / 67 / 65 / 40 / 69 words; gl(3): 24). V₃ and V₄ᵃ have closed forms (the gl(3) reduced vertex operator summed over the four 3-subsets plus η-terms; §3). The reduction itself costs 0.3–1.2 ms (universal) and ≈ 4 ms per λ′ (all three operators incl. building the application tries): negligible.
3. **They preserve M_μ (0 failures on every sector tested), commute exactly, and reproduce G2's spectrum**: exact characteristic polynomials of V₃, V₄ᵃ, V₄ᵇ and of G2's generic combination are identical to those of G2's full-space restriction on all 46 sectors; the joint eigenvalue triples agree to ≤ 1.3e-39 (G2's 40 digits). The invariance of M_μ is explained by an exact **BGG intertwining relation** E_i^{k_i} V_red(λ′) = V_red(s_i·λ′) E_i^{k_i} (5814 tests, 1130 non-trivial, 0 failures); the operators even commute on the **whole** weight space V_λ[μ−λ′], for generic non-dominant λ′ too (§7): they are the Gaudin operators on Sing(V_λ⊗M(λ′)) and M_μ is a spectral subspace.
4. **The scalars** of G2's relations: t²S_{3,1} = −2V₃ + c, t³S_{4,1} = −3(V₄ᵃ+V₄ᵇ) + c′ (tables §6). I found **exact closed forms** for c and c′ as polynomials in |λ|, |λ′| and the Perelomov–Popov Casimirs C₂, C₃, C₄ of λ, λ′, μ (fitted on random small products, tested on 252 irreps incl. products not used for the fit, 0 mismatches; not proved); c is also the constant in the site-exchange relation P V₃ P = −V₃ + c, and c′ the one in P V₄ᵃ P = −2V₄ᵃ − V₄ᵇ + c′ (all λ = λ′ sectors tested).
5. **Canonical copies** (joint eigenvectors, 50 digits, ordering by (ev₃, ev₄ᵃ, ev₄ᵇ), first nonzero slice component = 1): every V_red,k is diagonal in the eigenbasis to ≤ 8.5e-43 on all sectors of all products (worst case [12,8,4,0]², m = 51; [9,6,3,0]²: 5.0e-43; criterion 10^{−(dps−10)} = 1e-40), all eigenvalues real, each V_red,k G-self-adjoint (exact, 0 violations), the copies G⊗G-orthogonal (≤ 5.8e-47 with the exact Gram, ≤ 1.0e-15 with the float full vectors), HW residuals ≤ 9.9e-15, and the eigenvalue triples (≤ 1.3e-39) and slices (≤ 6.7e-39) agree with G2's exact HW eigenvectors on all sectors of the six G3.3 products. The site exchange acts on the copies as predicted by G2: the affine map T permutes the joint triples on all 156 sectors tested (0 failures), the full-vector overlap matrices are signed permutation matrices (≤ 9.6e-16); V₃-degeneracy is forced by |p₊−p₋| ≥ 2 exactly as G2 explained (13 forced + 6 accidental of the 19 V₃-degenerate sectors of [6,4,2,0]²).
6. **Separation of the copies at scale** (50-digit eigenvalues): the replication of G2's table on [6,4,2,0]² is exact (V₃ degenerate on 19 sectors, {V₃,V₄ᵃ} on 1, {V₃,V₄ᵇ} and the triple on none) and extends it to [9,6,3,0]² (179 sectors with m ≥ 2, m ≤ 27: V₃ degenerate on 69, **{V₃,V₄ᵃ} and {V₃,V₄ᵇ} on none**) and [12,8,4,0]² (409 sectors, m ≤ 51: {V₃,V₄ᵃ} fails on 1, {V₃,V₄ᵇ} on none); on the scaling family [2a,a,0,0]² (a ≤ 8, up to 484 sectors, m ≤ 9) V₃ alone separates every sector (V₄ᵃ ≡ 0 there).
7. **Timing** (§9): t_a = 0.009 / 0.041 / 0.17 / 0.52 / 1.38 / 3.16 s for [2a,a,0,0]², a = 1..6 (dim V up to 6860, 728 copies, n_W up to 29 028; t_b = 0.010 ... 3.89 s) and 0.049 / 1.47 / 23.4 s for [3a,2a,a,0]², a = 1, 2, 3 (1592 copies, m up to 27; 6.4 s with a float64 eigen-decomposition instead of the 50-digit one); a = 4 ([12,8,4,0]², 5945 copies, m up to 51): t_a = 258 s (37.2 s with float64 eigen-decomposition). The reduction is ≈ 4 ms per λ′; the reduced operators (exact) and, for m ≥ 10, the 50-digit eigen-decomposition dominate, not the recursion.
8. **Reference files**: six `gl4-cg-copies-*.txt` (§10) re-verified from the text with code independent of the reduction (the full two-site operators of G2 on the exact HW vectors, G1's exact and float recursions): 0 mismatches; a mutation test detects 17 of 17 single-record corruptions while the pristine file passes.
9. **Not done / limits**: see §12.

## 0b. Acceptance against brief §2

| item | requirement | result |
|---|---|---|
| G3.1 | `normal_order`; ≥ 200 random words of length ≤ 4 in gl(4), exact matrix on [2,1,0,0], [1,1,0,0]; the length-3 site-2 words of S_{4,1} | 656 words (+ PBW cross-check) and the 96 words: 0 mismatches (§2) |
| G3.2 | `reduce_terms`, `vred_apply`, `V_red`; N = 3: reduced tS_{3,2} = `slice_engine.vred_apply` + η₁/2 exactly on [2,1,0]², [3,1,0]², [4,2,0]²; S_{3,1} | 35 PRV vectors, 0 mismatches; S_{3,1}: 26 sectors, integer c (§4) |
| G3.3 | direct route on the six products (exact; V₃, V₄ᵃ, V₄ᵇ, S_{3,1}, S_{4,1}); M_μ preserved, commute, G2's spectrum, c, c′ | 0 mismatches, exact (slice and full vector); 0 non-commuting pairs; identical char. polys, joint eigenvalues ≤ 1.3e-39; c, c′ tables and closed forms (§5, §6) |
| G3.4 | `copies(dps=50)`, diagonality to 10^{−(dps−10)}, ordering, normalisation; real eigenvalues; self-adjointness; G⊗G-orthogonality; agreement with G2's eigenbasis (all sectors of [3,2,1,0]² and [2,1,0,0]²); exchange structure for λ = λ′ | all met (worst diagonality 8.5e-43 at m = 51); exchange table with pairings and parities (§8) |
| G3.5 | t_a / t_b for [2a,a,0,0]² (a ≤ 6), [3a,2a,a,0]² (a = 1, 2, 3), [S,S,0,0]² (S = 4, 8); t_red counted separately; honest comparison with the 24-term gl(3) V_red | §9 (+ extension a = 7, 8 and stress a = 4) |
| G3.6 | six `ref/gl4-cg-copies-*.txt` and `verify_ref` with 0 mismatches at 1e-12 | §10: 0 mismatches in all six; mutation test 17/17 detected |
| G3.7 | this file, `slice_ops_gl4_output.txt`, `test_slice_ops_gl4.py` (few-second exact regression), CLI stages `order | reduce | gl3 | direct | copies | bench | ref | verify | all` | done |

## 1. The reduction theorem as implemented (in my own words)

**Setting** (G1/G2 conventions): site 1 = V_λ (first factor), site 2 = V_λ′, e¹_ij / e²_ij = E_ij acting on site 1 / 2 (they commute), Ψ ∈ V_λ⊗V_λ′ a highest-weight (HW) vector of weight μ, its slice ψ = (id⊗⟨top|)Ψ ∈ M_μ ⊂ V_λ[μ−λ′] (⟨top| = coefficient functional of the HW state e_top of site 2). For a word W² of site-2 generators write its PBW normal form in U(gl_N), with the class order lowering < Cartan < raising,

  W² = Σ_w c_w · L_w · H_w · R_w ,  L_w = product of lowering generators E_ab (a > b), H_w = product of Cartan generators E_aa, R_w = E_{r₁}E_{r₂}⋯E_{r_k} product of raising generators (a < b).

**Theorem (reduction rule).** For every word W¹ of site-1 generators

  (id⊗⟨top|) (W¹ W²) Ψ = Σ_{w : L_w = 1}  c_w · Π_{E_aa ∈ H_w} λ′_a · (−1)^{k_w} · W¹ E_{r_k} ⋯ E_{r₂} E_{r₁} ψ        (site-1 generators acting on ψ ∈ V_λ[μ−λ′], right to left).

So the rightmost raising factor of W² is transferred first and the transferred factors stand in **reversed order**, to the right of W¹, with the sign (−1)^k; every Cartan factor of site 2 is evaluated to λ′ (site 2's HW), every term with a lowering factor in the normal form dies — but only after normal ordering (E_12 E_21 = E_21 E_12 + E_11 − E_22 contributes λ′₁ − λ′₂).

**Proof.** (i) ⟨top|E_ab = 0 for a > b (a lowering operator never produces the top state: λ′ + (positive root) is not a weight of V_λ′); (ii) ⟨top|E_aa = λ′_a⟨top|; (iii) E^tot_ab Ψ = 0 for every a < b (HW), i.e. e²_ab Ψ = −e¹_ab Ψ, and since e² commutes with e¹, e²_ab (X¹Ψ) = X¹ e²_abΨ = −X¹ e¹_abΨ for any site-1 word X¹. Hence R Ψ = e²_{r₁}⋯e²_{r_k}Ψ = (−1)^k e¹_{r_k}⋯e¹_{r₁} Ψ (the rightmost factor acts first on Ψ and each subsequent transfer passes through the site-1 factors already produced), and ⟨top| L H R Ψ = [L=1] Πλ′_a ⟨top|RΨ. Site-1 operators commute with ⟨top|. ∎

**Consequences.** (a) Red_λ′ : U(gl_N)⊗U(gl_N) → U(gl_N) is a well-defined linear map (PBW normal forms are unique), polynomial in λ′ (the universal polynomials are what `reduce_terms_poly` returns), independent of λ and μ except through the action on V_λ[μ−λ′]. (b) For O in the commutant of the diagonal gl_N, OΨ is again HW, so slice(OΨ) = Red(O)ψ ∈ M_μ and O ↦ Red(O)|_{M_μ} is the matrix of O on the multiplicity space in slice coordinates (algebra homomorphism: slice(O₁O₂Ψ) = Red(O₁) slice(O₂Ψ)). (c) For one site-2 factor the rule is the gl(3) one-liner (lowering → 0, Cartan → λ′_c, raising → −E¹); the new ingredient for gl(4) is the PBW normal ordering for words with ≥ 2 site-2 factors (V₄ᵇ, S_{3,1}, S_{4,1}, the antisymmetrised forms).


## 2. G3.1 `normal_order(word)` — the PBW engine

`normal_order(word)` reorders a word of gl(N) generators (tuple of `(i, j)`, product order, leftmost factor acts last) into the class order lowering < Cartan < raising (lexicographic inside a class) by adjacent transpositions E_a E_b = E_b E_a + [E_a, E_b], [E_ij, E_kl] = δ_jk E_il − δ_li E_kj; memoised (`lru_cache`), returns `{word: integer coefficient}`.

| test | result |
|---|---|
| 656 words of gl(4): all 256 words of length 2 + 400 random words of length 1..4 (seed 11): Σ c·(matrix of the output words) = matrix of the input word, **exactly**, on [2,1,0,0] (dim 20) and [1,1,0,0] (dim 6) | 0 mismatches |
| the same identity in the **PBW normal form of G2** (`gaudin_gl4_ops.pbw_word`: lexicographic order, a different total order, representation-free) | 0 mismatches |
| every output word is (lowering)(Cartan)(raising), each class sorted; `normal_order` is idempotent on its own output | 0 violations |
| the 96 distinct length-3 site-2 words occurring in S_{4,1} (exact matrices + PBW) | 0 mismatches; all 96 are single words (their normal form needs no commutator: see §3) |
| e12 e21 = e21 e12 + e11 − e22 (the pitfall of the brief); e21 e12 already normal | both reproduced |

(max number of terms of a normal-ordered random word of length ≤ 4: 11.)


## 3. G3.2 The reduction: engine, universal polynomials, closed forms, sparsity

**Engine** (`reduce_terms_poly`, `reduce_terms`, `ReducedOp`, `vred_apply`, `V_red`). `reduce_terms_poly(terms, N)` maps a term list {(W¹, W²): c} to the universal polynomial {site-1 word: {λ′-monomial: coefficient}} (normal-order W², drop words with a lowering factor, Cartan E_aa → λ′_a as a monomial, raising factors appended reversed with the sign); `reduce_terms(terms, lamp)` evaluates it at λ′ = lamp (exact Fractions; any weight, also non-dominant or rational). A `ReducedOp` stores the word list as a prefix trie and evaluates W_p v = c_p v + Σ_g E_g (W_{pg} v) (each trie node costs ONE generator application on a sum of vectors, dead branches are cut, Cartan generators act as the weight scalar, the non-adjacent generators are G1's lazily built exact matrices); `V_red(lam, lamp, mu, name)` applies it to the PRV basis (`PRV.coords` raises if an image leaves M_μ) and returns the exact m × m matrix, column j = V_red b_j. Operator names: `V3`, `V4a`, `V4b` (the t-free triple of the file), `S31`, `S41`, `S_k_j`, and `S_k_j_sw` (the θ₂ data: sites exchanged).

**Term counts** (from the output; "killed outright" = every normal-ordered site-2 word of the term has a lowering factor):

| operator | terms in the file | distinct site-2 words (max length) | killed outright | universal reduced words (max length; time) | words at λ′ = [2,1,1,0] / [3,2,1,0] / [6,4,2,0] / [12,6,0,0] / [8,8,0,0] | trie nodes at [3,2,1,0] |
|---|---|---|---|---|---|---|
| V₃ | 88 | 16 (1) | 29 | **39** (3; 0.3 ms) | 39 / 39 / 39 / 36 / 36 | 70 |
| V₄ᵃ | 152 | 16 (1) | 57 | **78** (4; 0.5 ms) | 59 / 67 / 67 / 61 / 61 | 146 |
| V₄ᵇ | 324 | 88 (2) | 173 | **85** (4; 0.9 ms) | 63 / 65 / 69 / 67 / 67 | 144 |
| S_{3,1} | 156 | 76 (2) | 71 | **40** (3; 0.6 ms) | 38 / 40 / 40 / 36 / 36 | 70 |
| S_{4,1} | 412 | 184 (3) | 246 | **86** (4; 1.2 ms) | 66 / 69 / 68 / 65 / 65 | 150 |

(gl(3): the 36-term V reduces to 24 words, Φ_ab = 0 for a > b removes 12.) V₃ reduces to the polynomial printed in the output (39 words; coefficients linear in λ′ for the quadratic words E_ij E_ji, constants ±1 for the 24 cubic words; V₄ᵃ, V₄ᵇ are printed there too). All of them have total weight 0 in site 1 (every word preserves the weight space V_λ[μ−λ′]).

**Closed forms (independent cross-check for the operators with one site-2 factor).** G2's exact identities V₃ = Σ_{|I|=3} D(2,1)_I + η, V₄ᵃ = D(3,1) + V₃ − η/3, η = n₁n₂ − Ω (D(a,b)_I = (1/a!b!) Σ_{σ,τ∈S(I)} sgn σ sgn τ e¹…e¹ e²…e²) have one site-2 factor per term, so the gl(3) one-line rule (`hand_reduce_single`, no normal ordering) gives

  V_red,3 = Σ_{|I|=3} ½ Σ_{σ,τ∈S(I)} sgn σ sgn τ E_{σ₁τ₁} E_{σ₂τ₂} Φ_{σ₃τ₃} + |λ||λ′| − Σ_i λ′_i E_ii + Σ_{i>j} E_ij E_ji ,
  V_red,4a = (1/6) Σ_{σ,τ∈S₄} sgn σ sgn τ E_{σ₁τ₁} E_{σ₂τ₂} E_{σ₃τ₃} Φ_{σ₄τ₄} + V_red,3 − (1/3)(|λ||λ′| − Σ_i λ′_i E_ii + Σ_{i>j} E_ij E_ji) ,

with Φ_cc = λ′_c, Φ_ab = −E_ab (a < b), Φ_ab = 0 (a > b), i.e. the **gl(3) reduced vertex operator on each 3-subset of {1,2,3,4}** plus the η-terms. These (71/337 hand-reduced words) agree with the engine's reduction of the file's lists (88 / 152 terms, a different polynomial) on every PRV vector of nine products (898 (operator, vector) pairs, 0 mismatches).

**Where the normal ordering matters — and where it does not.** In a column-ordered cdet word a raising factor E_{ac} (a < c) never stands left of a lowering factor E_{rc′} (c′ > c) it fails to commute with (r > c′ > c > a ⇒ [E_{ac}, E_{rc′}] = 0), the same holds for the θ₂ words (same word order), so for **all operators of the file the naive rule (lowering factor anywhere in the word ⇒ 0, no reordering) gives the same operator** as the full normal ordering (checked: deliberately wrong "naive_order" rule passes every test on the cdet lists, table below). The normal-ordering step is nevertheless essential for any other representative: **G2's antisymmetrised V₄ᵇ = D(2,2) − D(3,1) + ½ΣD(1,2) − ½ΣD(2,1) − η/6 (1464 terms, site-2 words e²_{s₃t₃}e²_{s₄t₄} in arbitrary index order)** reduced with the engine reproduces the file's V₄ᵇ on **449 of 449** PRV vectors (9 products), whereas the naive rule differs on 441 of 449. The brief's example is the word e²_{12}e²_{21} (raising left of lowering, contributes λ′₁−λ′₂), not e²_{21}e²_{12} (already normal ordered, dies): the roles are exchanged in the brief (§11).

**Sensitivity of the tests to the rule** (deliberately wrong rules on [2,1,0,0]², [3,2,1,0]², [2,2,1,0]⊗[2,1,1,0]; failing PRV vectors of all, same two tests as for the real rule — "image outside M_μ" via `PRV.coords`, "wrong slice" against the full two-site operator): omitting the sign (−1)^{|R|}: V₃ 29/30 wrong slices on [3,2,1,0]² (13/30 images outside M_μ), V₄ᵃ 24/30, V₄ᵇ 24/30, S_{3,1} 25/30, S_{4,1} 24/30; transferring in the original instead of the reversed order: V₃, V₄ᵃ unaffected (one raising factor), V₄ᵇ 22/30 wrong (13/30 outside M_μ), S_{3,1} 25/30, S_{4,1} 22/30; no normal ordering: undetected on the cdet lists (above), detected on the antisymmetrised V₄ᵇ. Both the "image leaves M_μ" test and the "slice differs from the direct route" test see the sign and order errors.

**Sparsity** (how often the reduced words act trivially; `ReducedOp.apply`): trie nodes evaluated / dead (zero vector or zero image): [3,2,1,0]²: V₃ 2100 / 668 (32 %), V₄ᵃ 4380 / 1609 (37 %), V₄ᵇ 4320 / 1593 (37 %); [4,2,0,0]²: 50 % / 62 % / 61 %; [6,4,2,0]²: 13 % / 15 % / 16 %; (word, weight-space basis vector) pairs with w e_c = 0: 31–40 % for [3,2,1,0]², 36–46 % for [4,2,0,0]², 13–17 % for [6,4,2,0]² (the larger the weight spaces, the denser the action).

## 4. G3.2 The N = 3 convention test

| product | irreps | PRV vectors | reduced t S_{3,2} = `slice_engine.vred_apply` + η₁/2 (mismatches) | t²S_{3,1,red} + 2·(t S_{3,2,red}) = c·1 on every M_μ |
|---|---|---|---|---|
| [2,1,0]² | 5 | 6 | 0 | all 5; c = 0, 12, 0, 15, 24 for μ = (4,2,0), (4,1,1), (3,3,0), (3,2,1), (2,2,2) |
| [3,1,0]² | 8 | 10 | 0 | all 8; c = 0, 16, 0, 21, 0, 24, 36, 40 |
| [4,2,0]² | 13 | 19 | 0 | all 13; c = 0, 40, 60, 0, 45, 72, 0, 48, 80, 96, 84, 105, 120 |

(η₁ = |Λ|² + cas₂(Λ) − cas₂(μ)/2; c − η₁ = −4, 6, −6, 6, 12 for [2,1,0]²; the integers c are those of the relation with the *reduced t·S_{3,2}* = V_red + η₁/2.) Total: 35 PRV vectors, 0 mismatches; 26 sectors, 0 failures of the scalar relation; direct route (the full two-site tS_{3,2} and S_{3,1} of G2's N = 3 engine applied exactly to the exact HW vectors, slice AND full vector): 0 mismatches of 70. So t·S_{3,2} = V + η/2 (G2) and the convention of the reduction rule are confirmed for N = 3 by the independent gl(3) engine; the N = 3 integer −2 of G2's pattern is confirmed by the direct route.

## 5. G3.3 The direct route on gl(4)

The direct route: the exact HW vectors Ψ_k of the PRV basis (G1's exact recursion; HW verified exactly), the full two-site operator (G2's Kronecker-sum builder `Op2`, term lists read from `ref/gl4-gaudin-ops.txt` at t = 1, exact Fractions, G2's private GT scaffold, entry-by-entry identical to `GLNRep`) applied to Ψ_j, (i) the slice (components (ν, top)) compared with V_red b_j = Σ_i A_ij b_i exactly, (ii) the whole vector O Ψ_j compared with Σ_i A_ij Ψ_i exactly (this also proves that O maps the HW space into itself). Nothing of the reduction code is used by the direct route.

### 5.1 The six products of G3.3 (exact; V₃, V₄ᵃ, V₄ᵇ, S_{3,1}, S_{4,1})

| product | irreps | PRV vectors | slice mismatches V3 / V4a / V4b / S31 / S41 | full-vector mismatches (same order) | image left M_μ | non-commuting pairs | scalar-relation failures | char-poly mismatches vs G2 | max joint-eigenvalue diff vs G2 |
|---|---|---|---|---|---|---|---|---|---|
| [2,1,0,0]^2 | 7 | 8 | 0 / 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 / 0 | 0 | 0 of 70 | 0 | 0 | 1.75e-41 |
| [2,1,1,0]^2 | 6 | 7 | 0 / 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 / 0 | 0 | 0 of 60 | 0 | 0 | 5.37e-40 |
| [2,2,1,0]x[2,1,1,0] | 6 | 7 | 0 / 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 / 0 | 0 | 0 of 60 | 0 | 0 | 5.01e-40 |
| [1,1,0,0]x[2,1,0,0] | 4 | 4 | 0 / 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 / 0 | 0 | 0 of 40 | 0 | 0 | 0.0 |
| [3,2,0,0]x[2,1,1,0] | 7 | 8 | 0 / 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 / 0 | 0 | 0 of 70 | 0 | 0 | 3.47e-40 |
| [3,2,1,0]^2 | 16 | 30 | 0 / 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 / 0 | 0 | 0 of 160 | 0 | 0 | 1.33e-39 |

Per sector: the five operators commute exactly (all 10 pairs); the relations t²S_{3,1} = −2V₃ + c·1, t³S_{4,1} = −3(V₄ᵃ+V₄ᵇ) + c′·1 hold with scalar matrices on every sector (0 failures, 46 sectors); the exact characteristic polynomials of V₃, V₄ᵃ, V₄ᵇ and of G2's generic combination V₃ + 7/19 V₄ᵃ + 11/13 V₄ᵇ equal those of G2's route (HWExact + restrict_exact, 0 restriction violations), and the 50-digit joint eigenvalue triples equal G2's 40-digit ones to the quoted maximum (= G2's noise level).

### 5.2 The direct route at scale (exact; slice read-off of the full operators V₃, V₄ᵃ, V₄ᵇ)

| product | irreps | PRV vectors × 3 operators | mismatches | time [s] |
|---|---|---|---|---|
| [4,2,0,0]^2 | 24 | 32 | 0 | 0.2 |
| [4,2,1,0]^2 | 27 | 56 | 0 | 0.3 |
| [3,2,1,0]x[4,2,1,0] | 20 | 38 | 0 | 0.1 |
| [4,2,1,0]x[3,2,1,0] | 20 | 38 | 0 | 0.2 |
| [6,3,0,0]^2 | 58 | 90 | 0 | 1.2 |
| [6,4,2,0]^2 | 76 | 297 | 0 | 8.5 |
| [2,1,0,0]x[1,1,0,0] | 4 | 4 | 0 | 0.0 |
| [2,1,1,0]x[2,2,1,0] | 6 | 7 | 0 | 0.0 |
| [2,1,1,0]x[3,2,0,0] | 7 | 8 | 0 | 0.0 |
| [2,1,0,0]x[3,2,1,0] | 10 | 14 | 0 | 0.0 |
| [3,2,1,0]x[2,1,0,0] | 10 | 14 | 0 | 0.1 |

Float variant (G1 float recursion + G2 float operators on the whole two-site vector; max over irreps of ‖slice(OΨ_j) − V_red b_j‖/max(‖V_red b_j‖, ‖b_j‖) in the orthonormal basis of site 1; V₄ᵃ ≡ 0 on λ₃ = 0 products, so the denominator is the vector norm):

* [4,2,0,0]^2 24 irreps: worst deviation V3 3.5e-16, V4a 1.7e-14, V4b 3.4e-16 (0.5s)
* [4,2,1,0]^2 27 irreps: worst deviation V3 3.4e-16, V4a 1.7e-15, V4b 1.1e-15 (0.9s)

### 5.3 The θ₂ Laurent coefficients and G2's relations (sites exchanged)

The θ₂ data are the same term lists with the sites exchanged (t-free hat swapped, S^{(2)}_{k,j}(t) = (−t)^{j−k} swap(hat)): their site-2 words are the former site-1 words (up to three generators in cdet column order). For (k,j) = (3,1), (3,2), (4,1), (4,2), (4,3): the reduced swapped operators preserve M_μ and agree exactly with the slice read-off of the full swapped two-site operators; G2's relations swap(hat S_{k,j}) = aV₃ + bV₄ᵃ + cV₄ᵇ + scalar with (a,b,c) = (2,0,0), (−1,0,0), (0,−3,−3), (0,3,2), (0,−2,−1) (G2 §6c, θ₂ column) hold with scalar remainders on every M_μ:

| product | irreps | PRV vectors | image left M_μ | slice vs FULL swapped operator (mismatches / checks) | G2 relations (failures / checks) | naive rule: wrong slices |
|---|---|---|---|---|---|---|
| [2,1,0,0]^2 | 7 | 8 | 0 | 0 / 40 | 0 / 35 | 0 of 40 |
| [2,1,1,0]^2 | 6 | 7 | 0 | 0 / 35 | 0 / 30 | 0 of 35 |
| [2,2,1,0]x[2,1,1,0] | 6 | 7 | 0 | 0 / 35 | 0 / 30 | 0 of 35 |
| [1,1,0,0]x[2,1,0,0] | 4 | 4 | 0 | 0 / 20 | 0 / 20 | 0 of 20 |
| [3,2,0,0]x[2,1,1,0] | 7 | 8 | 0 | 0 / 40 | 0 / 35 | 0 of 40 |
| [3,2,1,0]^2 | 16 | 30 | 0 | 0 / 150 | 0 / 80 | 0 of 150 |
| [4,2,1,0]^2 | 27 | 56 | 0 | 0 / 280 | 0 / 135 | 0 of 280 |
| [3,2,1,0]x[4,2,1,0] | 20 | 38 | 0 | 0 / 190 | 0 / 100 | 0 of 190 |
| [6,4,2,0]^2 | 76 | 297 | 0 | 0 / 1485 | 0 / 380 | 0 of 1485 |

(the naive rule gives the same operators here, see §3).

## 6. Commutativity, spectrum and the scalars c, c′

**On M_μ**: the three (indeed the five) reduced operators commute exactly on all 46 sectors of the six products (§5.1), their joint spectrum is that of G2, and it is **simple on every sector** (smallest max-norm distance between two joint triples: 6.708 for [2,1,0,0]² and [3,2,1,0]², 10.25 for [4,2,0,0]²/[4,2,1,0]²/[6,4,2,0]², 12.52 for [9,6,3,0]², 17.2 for [12,8,4,0]², 24.2 / 31.1 for [12,6,0,0]² / [16,8,0,0]²; see §9 for all families).

**The scalars** c (t²S_{3,1} = −2V₃ + c) and c′ (t³S_{4,1} = −3(V₄ᵃ+V₄ᵇ) + c′), per irrep (m):

* [2,1,0,0]^2: {4,2,0,0} (m=1): c=0, c'=0;  {4,1,1,0} (m=1): c=12, c'=0;  {3,3,0,0} (m=1): c=0, c'=0;  {3,2,1,0} (m=2): c=15, c'=0;  {3,1,1,1} (m=1): c=36, c'=36;  {2,2,2,0} (m=1): c=24, c'=0;  {2,2,1,1} (m=1): c=40, c'=40
* [2,1,1,0]^2: {4,2,2,0} (m=1): c=20, c'=0;  {4,2,1,1} (m=1): c=40, c'=56;  {3,3,2,0} (m=1): c=24, c'=0;  {3,3,1,1} (m=1): c=44, c'=60;  {3,2,2,1} (m=2): c=56, c'=72;  {2,2,2,2} (m=1): c=80, c'=120
* [2,2,1,0]x[2,1,1,0]: {4,3,2,0} (m=1): c=28, c'=0;  {4,3,1,1} (m=1): c=50, c'=70;  {4,2,2,1} (m=1): c=64, c'=84;  {3,3,3,0} (m=1): c=40, c'=0;  {3,3,2,1} (m=2): c=70, c'=90;  {3,2,2,2} (m=1): c=94, c'=144
* [1,1,0,0]x[2,1,0,0]: {3,2,0,0} (m=1): c=0, c'=0;  {3,1,1,0} (m=1): c=10, c'=0;  {2,2,1,0} (m=1): c=12, c'=0;  {2,1,1,1} (m=1): c=30, c'=30
* [3,2,0,0]x[2,1,1,0]: {5,3,1,0} (m=1): c=20, c'=0;  {5,2,1,1} (m=1): c=56, c'=64;  {4,4,1,0} (m=1): c=22, c'=0;  {4,3,2,0} (m=1): c=40, c'=0;  {4,3,1,1} (m=2): c=62, c'=70;  {4,2,2,1} (m=1): c=76, c'=84;  {3,3,2,1} (m=1): c=82, c'=90
* [3,2,1,0]^2: {6,4,2,0} (m=1): c=50, c'=0;  {6,4,1,1} (m=1): c=78, c'=108;  {6,3,3,0} (m=1): c=66, c'=0;  {6,3,2,1} (m=2): c=105, c'=135;  {6,2,2,2} (m=1): c=138, c'=216;  {5,5,2,0} (m=1): c=54, c'=0;  {5,5,1,1} (m=1): c=82, c'=112;  {5,4,3,0} (m=2): c=75, c'=0;  {5,4,2,1} (m=4): c=114, c'=144;  {5,3,3,1} (m=3): c=130, c'=160;  {5,3,2,2} (m=3): c=154, c'=240;  {4,4,4,0} (m=1): c=90, c'=0;  {4,4,3,1} (m=3): c=138, c'=168;  {4,4,2,2} (m=2): c=162, c'=252;  {4,3,3,2} (m=3): c=178, c'=280;  {3,3,3,3} (m=1): c=210, c'=360

**Closed forms** (found here; exact polynomial fits, **not proved**). With n₁ = |λ|, n₂ = |λ′| and C_r(ν) = Σ_i l_i^r Π_{j≠i}(1 − 1/(l_i − l_j)), l_i = ν_i + 4 − i (the Perelomov–Popov eigenvalue of Σ E_{i₁i₂}E_{i₂i₃}⋯E_{i_r i₁} on V_ν; C₂ is the quadratic Casimir):

  c = ½n₁²n₂ + (7/3)n₁n₂ + ½n₁n₂² + ½n₁C₂(λ) + ½n₂C₂(λ′) + (5/3)(C₂(λ)+C₂(λ′)) − ⅓(C₃(λ)+C₃(λ′)) − (5/3 + (n₁+n₂)/2) C₂(μ) + ⅓C₃(μ),

  c′ = (31 monomials of weight ≤ 4 in n₁, n₂, C_r(λ), C_r(λ′), C_r(μ), listed in `_CP_TERMS` and printed in the output; the terms without a factor n are (11/4)(C₂(λ)+C₂(λ′)−C₂(μ)) − ⅛(C₂(λ)²+C₂(λ′)²−C₂(μ)²) − (3/2)(C₃(λ)+C₃(λ′)−C₃(μ)) + ¼(C₄(λ)+C₄(λ′)−C₄(μ)); the rest has the factors n₁, n₂ and the leading part (1/6)n₁n₂(n₁²+n₂²) + ¼n₁²n₂²).

They were obtained by exact linear algebra over ℚ on the scalars of 90 (c) / 140 (c′) random small products (dim·dim′ ≤ 1500 / 2500, 552 sectors; the fit for c′ has rank 51 = the number of monomials and no inconsistent row) and then checked against the scalars computed from the reduced operators on **252 irreps of 12 products including products not used in the fit** ([3,2,1,0]², [4,2,1,0]², [3,2,1,0]⊗[4,2,1,0], [4,3,1,0]², [5,3,2,0]², [6,4,2,0]², [5,2,1,1]⊗[3,2,2,1]): 0 mismatches. For λ = λ′, c/2 = tr V₃/m and these c, c′ are the constants of the site-exchange relations of G2 §6e (verified on all λ = λ′ sectors, §8).

## 7. Structure of the reduced operators beyond the brief (exact arithmetic; observations, not proofs)

Three exact facts found while testing, all on the weight space W = V_λ[μ−λ′] ⊃ M_μ:

**(a) BGG intertwining.** Let k_i = λ′_i − λ′_{i+1} + 1 and s_i·λ′ = λ′ − k_i α_i (a non-dominant weight; the reduction is polynomial in λ′, so V_red(s_i·λ′) makes sense). For every operator k ∈ {V₃, V₄ᵃ, V₄ᵇ}, i = 1, 2, 3 and every GT basis vector of every weight space V_λ[μ−λ′] of a candidate μ

  E_i^{k_i} ∘ V_red,k(λ′) = V_red,k(s_i·λ′) ∘ E_i^{k_i}     (as maps V_λ[μ−λ′] → V_λ[μ−s_i·λ′]).

Tests: 90 / 72 / 81 / 36 / 171 / 324 / 531 / 693 / 3816 vectors-with-operators on the nine products of the table below (of which 8 / 6 / 8 / 0 / 28 / 36 / 85 / 134 / 825 are non-trivial, i.e. one side is nonzero): **0 failures**. Since M_μ = ∩_i ker E_i^{k_i} (the PRV theorem), this is the reason why the reduced operators preserve M_μ (until now only observed through `PRV.coords`). The maps E_i^{k_i}: V_λ[μ−λ′] → V_λ[μ−s_i·λ′] are the maps of the BGG/Jantzen analysis of V_λ ⊗ M(λ′) (the conditions E_i^{k_i}ψ = 0 are the PRV conditions); the relation says that the reduced Gaudin operators are compatible with the Verma-module embeddings M(s_i·λ′) ↪ M(λ′).

**(b) Commutation on the whole weight space, for every λ′.** The three reduced operators commute on the **whole** weight space V_λ[μ−λ′] of every candidate μ (not only on M_μ): 0 failures in 3 × (number of weight spaces) pairs for the nine products below (e.g. [6,4,2,0]²: 76 weight spaces, w up to 15 > m = 12). They also commute on every weight space of dimension ≥ 2 of V_λ for **generic, non-dominant, rational λ′** (four values such as (1/3, 5/7, −2, 4/5), including λ′ = 0): 168 / 456 / 12 tests for V_λ = [3,2,1,0], [4,2,0,0], [2,1,1,0], 0 failures. So V_red,k(λ′) are, for every λ′, commuting operators on V_λ[ν]. Natural interpretation (not proved here): for generic λ′ they are the Gaudin missing-label operators acting on the singular vectors of V_λ ⊗ M(λ′) (M(λ′) = Verma module; the slice identifies Sing_μ with V_λ[μ−λ′]), and the finite-dimensional case is the subspace M_μ cut out by the conditions E_i^{k_i}ψ = 0 of (a).

**(c) M_μ is a spectral subspace.** The joint spectrum of the triple on each whole weight space V_λ[μ−λ′] is **simple** (w distinct joint triples, 50 digits) in every one of the weight spaces tested, and the m joint triples of M_μ are among them (the other w − m eigenvectors violate the string conditions: "spurious"): [3,2,1,0]²: 16 weight spaces, Σ w = 36, Σ m = 30; [4,2,0,0]²: 26, 59, 32; [4,2,1,0]²: 27, 77, 56; [6,4,2,0]²: 76, 424, 297. So M_μ is spanned by joint eigenlines of the commuting family on V_λ[μ−λ′] — answering in the affirmative the WP-G4 question "is M_μ a spectral subspace?" for these cases.

| product | weight spaces | commutator failures (pairs) | intertwining tests (non-zero) | failures | joint spectrum simple on W | M_μ triples ⊂ spectrum | Σ w | Σ (w − m) |
|---|---|---|---|---|---|---|---|---|
| [2,1,0,0]^2 | 7 | 0 of 21 | 90 (8) | 0 | 7 of 7 | 7 of 7 | 10 | 2 |
| [2,1,1,0]^2 | 6 | 0 of 18 | 72 (6) | 0 | 6 of 6 | 6 of 6 | 8 | 1 |
| [2,2,1,0]x[2,1,1,0] | 6 | 0 of 18 | 81 (8) | 0 | 6 of 6 | 6 of 6 | 9 | 2 |
| [1,1,0,0]x[2,1,0,0] | 4 | 0 of 12 | 36 (0) | 0 | 4 of 4 | 4 of 4 | 4 | 0 |
| [3,2,0,0]x[2,1,1,0] | 10 | 0 of 30 | 171 (28) | 0 | 10 of 10 | 10 of 10 | 19 | 11 |
| [3,2,1,0]^2 | 16 | 0 of 48 | 324 (36) | 0 | 16 of 16 | 16 of 16 | 36 | 6 |
| [4,2,0,0]^2 | 26 | 0 of 78 | 531 (85) | 0 | 26 of 26 | 26 of 26 | 59 | 27 |
| [4,2,1,0]^2 | 27 | 0 of 81 | 693 (134) | 0 | 27 of 27 | 27 of 27 | 77 | 21 |
| [6,4,2,0]^2 | 76 | 0 of 228 | 3816 (825) | 0 | 76 of 76 | 76 of 76 | 424 | 127 |


## 8. G3.4 The canonical copies

`copies(lam, lamp, mu, dps=50)`: joint eigenvectors of (V_red,3, V_red,4a, V_red,4b) on M_μ from the mpmath eigen-decomposition of the generic rational combination M = V₃ + (7/19)V₄ᵃ + (11/13)V₄ᵇ (G2's), after which **every** V_red,k is conjugated to the eigenbasis and its largest off-diagonal entry measured (acceptance 10^{-(dps−10)} relative to max(1,|A|), precision raised by 30 digits if not met — never needed); copies ordered lexicographically ascending by (ev₃, ev₄ᵃ, ev₄ᵇ) (ties compared with tolerance 10^{-25}: V₃-degenerate sectors are ordered by ev₄ᵃ), slice normalised to the first nonzero component (ascending GT index) = 1, full vector by G1's float level recursion (and, for the files, by linearity from the exact vectors of the PRV basis). Checks per sector: real eigenvalues; each V_red,k G-self-adjoint for the restricted invariant form (exact: R_kᵀ·Gram = Gram·R_k with the exact G⊗G Gram matrix of the PRV basis' full HW vectors); G⊗G-orthogonality of the copies' full vectors (from the exact Gram with the 50-digit coordinates, and from the float recursion); HW residuals (all E_ij, float); agreement with G2's exact joint eigenbasis (slices of Σ u_i·(G2's exact HW basis vectors) with u from G2's restricted matrices, normalised identically).

| product | sectors | copies | joint-diag failures | worst off-diagonal of V_red,k in the eigenbasis | smallest joint-triple distance | self-adjointness violations (exact) | G-orth. (exact Gram) | G-orth. (float vectors) | worst HW residual | vs G2: max triple diff | vs G2: max slice diff |
|---|---|---|---|---|---|---|---|---|---|---|---|
| [2,1,0,0]^2 | 7 | 8 | 0 | 2.71e-51 | 6.708 | 0 | 4.45e-52 | 7.7e-17 | 3.5e-16 | 1.75e-41 | 1.17e-41 |
| [2,1,1,0]^2 | 6 | 7 | 0 | 7.88e-49 | 27.71 | 0 | 5.82e-50 | 6.8e-17 | 2.5e-15 | 5.37e-40 | 4.16e-39 |
| [2,2,1,0]x[2,1,1,0] | 6 | 7 | 0 | 5.56e-50 | 51.26 | 0 | 2.83e-51 | 7.0e-16 | 2.6e-15 | 5.01e-40 | 8.55e-41 |
| [1,1,0,0]x[2,1,0,0] | 4 | 4 | 0 | 0.0 | +inf | 0 | 0.0 | 0.0e+00 | 2.0e-16 | 0.0 | 0.0 |
| [3,2,0,0]x[2,1,1,0] | 7 | 8 | 0 | 4.55e-50 | 20.49 | 0 | 1.43e-51 | 1.0e-15 | 2.9e-15 | 3.47e-40 | 1.08e-40 |
| [3,2,1,0]^2 | 16 | 30 | 0 | 1.77e-48 | 6.708 | 0 | 1.24e-49 | 3.2e-16 | 2.4e-15 | 1.33e-39 | 6.66e-39 |
| [4,2,0,0]^2 | 24 | 32 | 0 | 1.4e-49 | 10.25 | 0 | 1.66e-50 | 2.2e-16 | 1.7e-15 | – | – |
| [4,2,1,0]^2 | 27 | 56 | 0 | 7.52e-48 | 10.25 | 0 | 7.0e-49 | 4.5e-16 | 3.7e-15 | – | – |
| [6,4,2,0]^2 | 76 | 297 | 0 | 1.16e-44 | 10.25 | 0 | 5.76e-47 | 8.9e-16 | 9.9e-15 | – | – |

The G2 comparison is made on all six products of G3.3 (the table's last two columns: max |ev − ev_G2|, max |slice − slice_G2| over all copies; G2 works at 40 digits).

### 8.1 Exchange structure (λ = λ′)

The site exchange P: (ν, q) → (q, ν) preserves every HW space (G2) and conjugates the triple: P V₃ P = −V₃ + c, P V₄ᵃ P = −2V₄ᵃ − V₄ᵇ + c′, P V₄ᵇ P = 3V₄ᵃ + 2V₄ᵇ − c′ (c, c′ fixed by the traces: c = 2 tr V₃/m, c′ = (3 tr V₄ᵃ + tr V₄ᵇ)/m). On the copies: (i) the spectrum of V₃ is symmetric about c/2 (to ≤ 1.4e-48); (ii) the affine map T(e₃, e₄ᵃ, e₄ᵇ) = (c − e₃, −2e₄ᵃ − e₄ᵇ + c′, 3e₄ᵃ + 2e₄ᵇ − c′) **permutes the joint eigenvalue triples** (partner σ(k): P Ψ_k ∝ Ψ_σ(k); fixed copies σ(k) = k have a P-parity ±, the others come in 2-cycles); (iii) the overlap matrix ⟨Ψ_j|PΨ_k⟩/norms of the float full vectors is the signed permutation matrix of σ (max deviation in the table); (iv) the exchange scalars (c, c′) are exactly the closed-form scalars of §6 on every sector tested. Per product:

* [2,1,0,0]^2 — exchange: 7 sectors, T-map failures 0, fixed copies 6, 2-cycles 1, max |overlap - signed permutation| 4.4e-16, exchange scalars (c, c') = the closed-form scalars of the S31 / S41 relations on 7 of 7 sectors; of the 1 sectors with m >= 2: |p+ - p-| >= 2 (V3 degeneracy forced): 0, V3 degenerate: 0 (not forced: 0), forced but V3 simple: 0 
* [2,1,1,0]^2 — exchange: 6 sectors, T-map failures 0, fixed copies 5, 2-cycles 1, max |overlap - signed permutation| 1.3e-15, exchange scalars (c, c') = the closed-form scalars of the S31 / S41 relations on 6 of 6 sectors; of the 1 sectors with m >= 2: |p+ - p-| >= 2 (V3 degeneracy forced): 0, V3 degenerate: 0 (not forced: 0), forced but V3 simple: 0 
* [3,2,1,0]^2 — exchange: 16 sectors, T-map failures 0, fixed copies 14, 2-cycles 8, max |overlap - signed permutation| 8.9e-16, exchange scalars (c, c') = the closed-form scalars of the S31 / S41 relations on 16 of 16 sectors; of the 8 sectors with m >= 2: |p+ - p-| >= 2 (V3 degeneracy forced): 1, V3 degenerate: 2 (not forced: 1), forced but V3 simple: 0 
* [4,2,0,0]^2 — exchange: 24 sectors, T-map failures 0, fixed copies 18, 2-cycles 7, max |overlap - signed permutation| 4.4e-16, exchange scalars (c, c') = the closed-form scalars of the S31 / S41 relations on 24 of 24 sectors; of the 7 sectors with m >= 2: |p+ - p-| >= 2 (V3 degeneracy forced): 0, V3 degenerate: 0 (not forced: 0), forced but V3 simple: 0 
* [4,2,1,0]^2 — exchange: 27 sectors, T-map failures 0, fixed copies 22, 2-cycles 17, max |overlap - signed permutation| 4.4e-16, exchange scalars (c, c') = the closed-form scalars of the S31 / S41 relations on 27 of 27 sectors; of the 14 sectors with m >= 2: |p+ - p-| >= 2 (V3 degeneracy forced): 2, V3 degenerate: 2 (not forced: 0), forced but V3 simple: 0 
* [6,4,2,0]^2 — exchange: 76 sectors, T-map failures 0, fixed copies 69, 2-cycles 114, max |overlap - signed permutation| 9.6e-16, exchange scalars (c, c') = the closed-form scalars of the S31 / S41 relations on 76 of 76 sectors; of the 56 sectors with m >= 2: |p+ - p-| >= 2 (V3 degeneracy forced): 13, V3 degenerate: 19 (not forced: 6), forced but V3 simple: 0 

The explicit exchange table of [3,2,1,0]² (all sectors with m ≥ 2; copies numbered 1..m in the lexicographic order of (ev₃, ev₄ᵃ, ev₄ᵇ); c, c′ are the scalars of §6; "fixed" = copies with P Ψ = ±Ψ, parity in brackets; "2-cycles" = pairs (k, σ(k)) with P Ψ_k ∝ Ψ_σ(k)):

| mu | m | c | c' | fixed copies (P-parity) | 2-cycles (partner copies) | max abs(overlap - signed permutation) |
|---|---|---|---|---|---|---|
| (6,3,2,1) | 2 | 105.0 | 135.0 | []  | [(1,2)] | 2.3e-16 |
| (5,4,3,0) | 2 | 75.0 | 0.0 | []  | [(1,2)] | 2.3e-16 |
| (5,4,2,1) | 4 | 114.0 | 144.0 | []  | [(1,4),(2,3)] | 2.2e-16 |
| (5,3,3,1) | 3 | 130.0 | 160.0 | [2] (+) | [(1,3)] | 4.6e-16 |
| (5,3,2,2) | 3 | 154.0 | 240.0 | [2] (-) | [(1,3)] | 5.3e-16 |
| (4,4,3,1) | 3 | 138.0 | 168.0 | [2] (-) | [(1,3)] | 7.4e-16 |
| (4,4,2,2) | 2 | 162.0 | 252.0 | [1, 2] (++) | [] | 2.2e-16 |
| (4,3,3,2) | 3 | 178.0 | 280.0 | [2] (-) | [(1,3)] | 8.9e-16 |

and of the sectors of [6,4,2,0]² with m ≥ 6 (the other sectors, including all m < 6, are in the output):

| mu | m | c | c' | fixed copies (P-parity) | 2-cycles (partner copies) | max abs(overlap - signed permutation) |
|---|---|---|---|---|---|---|
| (11,7,4,2) | 7 | 684.0 | 1260.0 | [4] (-) | [(1,7),(2,6),(3,5)] | 7.1e-16 |
| (11,6,5,2) | 6 | 729.0 | 1344.0 | []  | [(1,6),(2,5),(3,4)] | 3.5e-16 |
| (11,6,4,3) | 6 | 792.0 | 1680.0 | []  | [(1,6),(2,5),(3,4)] | 6.4e-16 |
| (10,8,5,1) | 7 | 620.0 | 780.0 | [4] (-) | [(1,7),(2,6),(3,5)] | 7.7e-16 |
| (10,8,4,2) | 10 | 708.0 | 1300.0 | [5, 6] (++) | [(1,10),(2,9),(3,8),(4,7)] | 5.2e-16 |
| (10,7,6,1) | 6 | 659.0 | 819.0 | []  | [(1,6),(2,5),(3,4)] | 6.3e-16 |
| (10,7,5,2) | 11 | 764.0 | 1404.0 | [7] (-) | [(1,11),(2,10),(3,9),(4,8),(5,6)] | 6.9e-16 |
| (10,7,4,3) | 10 | 827.0 | 1755.0 | []  | [(1,10),(2,9),(3,8),(4,7),(5,6)] | 6.0e-16 |
| (10,6,6,2) | 6 | 792.0 | 1456.0 | [3, 4] (++) | [(1,6),(2,5)] | 6.6e-16 |
| (10,6,5,3) | 9 | 872.0 | 1872.0 | [5] (-) | [(1,9),(2,8),(3,7),(4,6)] | 6.2e-16 |
| (10,6,4,4) | 6 | 912.0 | 2080.0 | [3, 4] (++) | [(1,6),(2,5)] | 6.8e-16 |
| (9,8,6,1) | 6 | 680.0 | 840.0 | []  | [(1,6),(2,5),(3,4)] | 5.9e-16 |
| (9,8,5,2) | 10 | 785.0 | 1440.0 | []  | [(1,10),(2,9),(3,8),(4,7),(5,6)] | 6.1e-16 |
| (9,8,4,3) | 8 | 848.0 | 1800.0 | []  | [(1,8),(2,7),(3,4),(5,6)] | 6.1e-16 |
| (9,7,6,2) | 9 | 824.0 | 1512.0 | [5] (-) | [(1,9),(2,8),(3,7),(4,6)] | 6.4e-16 |
| (9,7,5,3) | 12 | 904.0 | 1944.0 | [5, 7] (++) | [(1,12),(2,10),(3,9),(4,11),(6,8)] | 6.0e-16 |
| (9,7,4,4) | 6 | 944.0 | 2160.0 | [3, 4] (--) | [(1,6),(2,5)] | 6.3e-16 |
| (9,6,6,3) | 7 | 932.0 | 2016.0 | [4] (-) | [(1,7),(2,6),(3,5)] | 6.5e-16 |
| (9,6,5,4) | 8 | 989.0 | 2304.0 | []  | [(1,8),(2,7),(3,6),(4,5)] | 9.1e-16 |
| (8,8,6,2) | 6 | 840.0 | 1540.0 | [3, 4] (++) | [(1,6),(2,5)] | 3.5e-16 |
| (8,8,5,3) | 6 | 920.0 | 1980.0 | [3, 4] (--) | [(1,6),(2,5)] | 5.5e-16 |
| (8,7,6,3) | 8 | 959.0 | 2079.0 | []  | [(1,8),(2,7),(3,6),(4,5)] | 8.2e-16 |
| (8,7,5,4) | 7 | 1016.0 | 2376.0 | [5] (-) | [(1,7),(2,6),(3,4)] | 6.4e-16 |
| (8,6,6,4) | 6 | 1044.0 | 2464.0 | [3, 4] (++) | [(1,6),(2,5)] | 3.8e-16 |

**V₃-degeneracy**: with p₊/p₋ the numbers of P-symmetric/antisymmetric copies (p₊ − p₋ = Σ of the P-parities of the fixed copies), V₃ − c/2 is odd under P, has a symmetric spectrum and at least |p₊ − p₋| zero modes: **|p₊ − p₋| ≥ 2 forces V₃ to be degenerate** (13 sectors of [6,4,2,0]², 1 of [3,2,1,0]², 2 of [4,2,1,0]²; in every one V₃ is indeed degenerate: "forced but simple: 0"). The remaining V₃-degenerate sectors are accidental: (5,4,2,1) of [3,2,1,0]² (m = 4, p₊ = p₋ = 2, no P-fixed copy: the 2-cycle (2,3) pairs the copies (57, 21.303, 50.697) ↔ (57, 50.697, 21.303), both with V₃ = c/2 = 57), and 6 sectors of [6,4,2,0]². In the tables above, e.g. (4,4,2,2): fixed copies 1, 2, parities ++ (p₊ − p₋ = 2), V₃ = 81 = c/2 on both; (9,7,5,3) (m = 12): fixed copies 5, 7, parities ++.

### 8.2 Separation of the copies by subsets of the triple, at scale

Number of sectors (with m ≥ 2) on which the spectrum of the given subset is degenerate (two copies with equal eigenvalues of the subset to 1e-20; 50-digit eigenvalues):

| family | a | λ | sectors with m ≥ 2 | degenerate on {V₃} / {V₄ᵃ} / {V₄ᵇ} / {V₃,V₄ᵃ} / {V₃,V₄ᵇ} / {V₄ᵃ,V₄ᵇ} / {V₃,V₄ᵃ,V₄ᵇ} |
|---|---|---|---|---|
| [2a,a,0,0]² | 1 | [2,1,0,0] | 1 | 0 / 1 / 1 / 0 / 0 / 1 / 0 |
| [2a,a,0,0]² | 2 | [4,2,0,0] | 7 | 0 / 7 / 7 / 0 / 0 / 7 / 0 |
| [2a,a,0,0]² | 3 | [6,3,0,0] | 24 | 0 / 24 / 24 / 0 / 0 / 24 / 0 |
| [2a,a,0,0]² | 4 | [8,4,0,0] | 58 | 0 / 58 / 58 / 0 / 0 / 58 / 0 |
| [2a,a,0,0]² | 5 | [10,5,0,0] | 115 | 0 / 115 / 115 / 0 / 0 / 115 / 0 |
| [2a,a,0,0]² | 6 | [12,6,0,0] | 201 | 0 / 201 / 201 / 0 / 0 / 201 / 0 |
| [3a,2a,a,0]² | 1 | [3,2,1,0] | 8 | 2 / 1 / 1 / 0 / 0 / 1 / 0 |
| [3a,2a,a,0]² | 2 | [6,4,2,0] | 56 | 19 / 8 / 6 / 1 / 0 / 6 / 0 |
| [3a,2a,a,0]² | 3 | [9,6,3,0] | 179 | 69 / 13 / 13 / 0 / 0 / 13 / 0 |
| [2a,a,0,0]² (extension, a = 7, 8) | 7 | [14,7,0,0] | 322 | 0 / 322 / 322 / 0 / 0 / 322 / 0 |
| [2a,a,0,0]² (extension, a = 7, 8) | 8 | [16,8,0,0] | 484 | 0 / 484 / 484 / 0 / 0 / 484 / 0 |
| [3a,2a,a,0]² (extension, a = 4) | 4 | [12,8,4,0] | 409 | 169 / 27 / 25 / 1 / 0 / 25 / 0 |

For [6,4,2,0]² this reproduces G2's exact-discriminant table exactly (V₃: 19, V₄ᵃ: 8, V₄ᵇ: 6, {V₃,V₄ᵃ}: 1 [the m = 12 sector (9,7,5,3)], {V₃,V₄ᵇ}: 0, {V₄ᵃ,V₄ᵇ}: 6, triple: 0). New: [9,6,3,0]² (179 sectors, m ≤ 27) and [12,8,4,0]² (409 sectors, m ≤ 51) are separated by {V₃,V₄ᵇ} on **every** sector; {V₃,V₄ᵃ} fails on none at a = 3 and on one sector at a = 4 (the sector (18,14,10,6), m = 50, where four copies with V₃ = c/2 = 3084 share (V₃, V₄ᵃ) = (3084, 5292) and are separated by V₄ᵇ alone (5308.34, 5950.55, 6300 [P-fixed], 6873.11); in the m = 51 sector (19,14,10,5) V₃ is degenerate but {V₃,V₄ᵃ} separates); on the scaling family [2a,a,0,0]² the cubic V₃ alone separates every sector (V₄ᵃ vanishes identically for λ₃ = 0 and V₄ᵇ is a function of V₃ on these, so {V₄ᵃ}, {V₄ᵇ}, {V₄ᵃ,V₄ᵇ} fail everywhere). No degenerate joint spectrum occurs anywhere: on all 1864 sectors with m ≥ 2 of the twelve family members of the table (1212 for [2a,a,0,0]², a ≤ 8; 652 for [3a,2a,a,0]², a ≤ 4) the triple (V₃, V₄ᵃ, V₄ᵇ) separates all copies (smallest joint-triple distance 6.7), i.e. it is a complete set of commuting missing-label operators there.

## 9. G3.5 Timing: the n_W-free part t_a against the full reconstruction t_b

**Protocol.** `time.perf_counter()`, single process, **min over cold runs** (before each run: G1's caches, the normal-ordering memo, the universal and evaluated reductions and G2's scaffold are cleared; 5 / 3 / 2 runs for dim V ≤ 500 / ≤ 5000 / larger, one run for the extension rows). Components:

* **t_a** (the n_W-free part; no n_W-sized object is touched) = t_gen (GT generators of both sites: adjacent + diagonal) + t_gen_na (the six non-adjacent generators of site 1, forced: the reduced words use them) + t_rs (Racah–Speiser incl. weight multiplicities) + t_prv (exact PRV spaces of all irreps) [G1's t_a = t_gen + t_rs + t_prv] + **t_red** (the reduction of the three operators for this λ′: universal polynomial + evaluation at λ′ + trie construction, **counted once per λ′ and cold**, ≈ 4–5 ms for every size) + **t_ops** (the three reduced operators applied exactly to every PRV basis vector of every irrep: the exact m × m matrices) + **t_eig** (the joint eigen-decomposition in mpmath at 50 digits, with the slices of all copies).
* **t_b** = t_a + t_setup (G1's float setup of the level recursion: metric, orthonormal float matrices, left inverses) + t_recon (the float reconstruction of EVERY full copy, batched per irrep).
* "float64 eig" columns: the 50-digit mpmath eigen-decomposition replaced by numpy float64 (`joint_diag_float`): t_a,f = t_a − t_eig + t_eig,f. Accuracy against the 50-digit result: max |ev − ev_mp|/max(1,|ev|) ≤ 2.1e-15 and worst off-diagonal ≤ 7.0e-12 (relative to max(1,|A|)) on all sectors of all rows below. This is the like-for-like number against the gl(3) engine (`slice_engine.bench_one`, whose t_a also uses a float eigen-decomposition).
* t_b for the extension rows is not measured (slices only).

**[2a,a,0,0]^2**

| a | λ | dim V | irreps | copies | max m | max w | max n_W | t_a [s] | t_b [s] | t_a, float64 eig [s] | t_b, float64 eig [s] | vector entries | worst HW residual |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | [2,1,0,0] | 20 | 7 | 8 | 2 | 2 | 20 | 0.009 | 0.010 | 0.008 | 0.009 | 59 | 3.5e-16 |
| 2 | [4,2,0,0] | 126 | 24 | 32 | 3 | 4 | 174 | 0.041 | 0.046 | 0.035 | 0.041 | 1287 | 1.7e-15 |
| 3 | [6,3,0,0] | 480 | 58 | 90 | 4 | 7 | 920 | 0.167 | 0.193 | 0.142 | 0.168 | 14208 | 6.3e-15 |
| 4 | [8,4,0,0] | 1375 | 115 | 205 | 5 | 11 | 3541 | 0.518 | 0.611 | 0.432 | 0.525 | 101092 | 9.1e-15 |
| 5 | [10,5,0,0] | 3276 | 201 | 406 | 6 | 15 | 10968 | 1.375 | 1.659 | 1.147 | 1.431 | 531186 | 1.5e-14 |
| 6 | [12,6,0,0] | 6860 | 322 | 728 | 7 | 20 | 29028 | 3.160 | 3.885 | 2.614 | 3.338 | 2236946 | 9.8e-14 |

| a | t_gen | t_gen_na | t_rs | t_prv | t_red | t_ops | t_eig (mpmath, 50 digits) | t_eig (float64) | t_setup | t_recon | words V₃/V₄ᵃ/V₄ᵇ | worst joint-diag off-diag. | min joint-triple distance | cold runs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.0005 | 0.0003 | 0.0001 | 0.0001 | 0.0041 | 0.0029 | 0.0007 | 0.0002 | 0.0006 | 0.0003 | 36/53/66 | 2.7e-51 | 6.708 | 5 |
| 2 | 0.0034 | 0.0030 | 0.0003 | 0.0005 | 0.0039 | 0.0236 | 0.0064 | 0.0007 | 0.0021 | 0.0032 | 36/61/67 | 1.4e-49 | 10.25 | 5 |
| 3 | 0.0128 | 0.0163 | 0.0009 | 0.0025 | 0.0043 | 0.1032 | 0.0268 | 0.0017 | 0.0074 | 0.0187 | 36/61/67 | 3.0e-48 | 13.75 | 5 |
| 4 | 0.0406 | 0.0518 | 0.0024 | 0.0100 | 0.0041 | 0.3192 | 0.0897 | 0.0039 | 0.0184 | 0.0751 | 36/61/66 | 8.9e-48 | 17.23 | 3 |
| 5 | 0.1309 | 0.1416 | 0.0050 | 0.0318 | 0.0043 | 0.8253 | 0.2362 | 0.0078 | 0.0434 | 0.2410 | 36/61/67 | 4.1e-47 | 20.71 | 3 |
| 6 | 0.2214 | 0.4296 | 0.0108 | 0.0913 | 0.0042 | 1.8423 | 0.5606 | 0.0140 | 0.0882 | 0.6367 | 36/61/67 | 9.7e-47 | 24.19 | 2 |

**[3a,2a,a,0]^2**

| a | λ | dim V | irreps | copies | max m | max w | max n_W | t_a [s] | t_b [s] | t_a, float64 eig [s] | t_b, float64 eig [s] | vector entries | worst HW residual |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | [3,2,1,0] | 64 | 16 | 30 | 4 | 4 | 152 | 0.049 | 0.052 | 0.035 | 0.038 | 1069 | 2.4e-15 |
| 2 | [6,4,2,0] | 729 | 76 | 297 | 12 | 15 | 4743 | 1.467 | 1.515 | 0.709 | 0.757 | 185606 | 6.5e-15 |
| 3 | [9,6,3,0] | 4096 | 213 | 1592 | 27 | 32 | 59008 | 23.418 | 23.904 | 6.422 | 6.909 | 8889611 | 3.9e-14 |

| a | t_gen | t_gen_na | t_rs | t_prv | t_red | t_ops | t_eig (mpmath, 50 digits) | t_eig (float64) | t_setup | t_recon | words V₃/V₄ᵃ/V₄ᵇ | worst joint-diag off-diag. | min joint-triple distance | cold runs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.0016 | 0.0016 | 0.0002 | 0.0003 | 0.0039 | 0.0273 | 0.0141 | 0.0007 | 0.0013 | 0.0016 | 39/67/65 | 1.8e-48 | 6.708 | 5 |
| 2 | 0.0213 | 0.0300 | 0.0012 | 0.0065 | 0.0040 | 0.6412 | 0.7630 | 0.0048 | 0.0101 | 0.0380 | 39/67/69 | 1.2e-44 | 10.25 | 3 |
| 3 | 0.1343 | 0.2525 | 0.0067 | 0.0793 | 0.0044 | 5.9174 | 17.0233 | 0.0274 | 0.0553 | 0.4313 | 39/67/69 | 5.0e-43 | 12.52 | 3 |

**[2a,a,0,0]^2 (extension, a = 7, 8)**

| a | λ | dim V | irreps | copies | max m | max w | max n_W | t_a [s] | t_b [s] | t_a, float64 eig [s] | t_b, float64 eig [s] | vector entries | worst HW residual |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 7 | [14,7,0,0] | 13056 | 484 | 1212 | 8 | 26 | 68208 | 6.355 | – | 5.251 | – | – | – |
| 8 | [16,8,0,0] | 23085 | 693 | 1905 | 9 | 33 | 146025 | 11.986 | – | 9.887 | – | – | – |

| a | t_gen | t_gen_na | t_rs | t_prv | t_red | t_ops | t_eig (mpmath, 50 digits) | t_eig (float64) | t_setup | t_recon | words V₃/V₄ᵃ/V₄ᵇ | worst joint-diag off-diag. | min joint-triple distance | cold runs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 7 | 0.5066 | 0.7687 | 0.0740 | 0.2299 | 0.0050 | 3.6442 | 1.1265 | 0.0224 | – | – | 36/61/67 | 2.0e-46 | 27.66 | 1 |
| 8 | 0.9188 | 1.6161 | 0.0317 | 0.5523 | 0.0048 | 6.7277 | 2.1341 | 0.0353 | – | – | 36/61/67 | 5.6e-46 | 31.13 | 1 |

**[3a,2a,a,0]^2 (extension, a = 4)**

| a | λ | dim V | irreps | copies | max m | max w | max n_W | t_a [s] | t_b [s] | t_a, float64 eig [s] | t_b, float64 eig [s] | vector entries | worst HW residual |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 4 | [12,8,4,0] | 15625 | 459 | 5945 | 51 | 65 | 426425 | 257.550 | – | 37.152 | – | – | – |

| a | t_gen | t_gen_na | t_rs | t_prv | t_red | t_ops | t_eig (mpmath, 50 digits) | t_eig (float64) | t_setup | t_recon | words V₃/V₄ᵃ/V₄ᵇ | worst joint-diag off-diag. | min joint-triple distance | cold runs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 0.7057 | 1.2213 | 0.0208 | 0.6978 | 0.0041 | 34.3659 | 220.5339 | 0.1364 | – | – | 39/67/69 | 8.5e-43 | 17.23 | 1 |

(the last two blocks are the extension: a = 7, 8 of the scaling family and a = 4 of the stress family, slices only: t_b and the HW residuals are not measured.)

**Breakdown at the top of each family** (best run): scaling family a = 6 ([12,6,0,0]², dim 6860, 322 irreps, 728 copies, n_W ≤ 29 028): t_a = 3.16 s = generators 0.22 + non-adjacent generators 0.43 + Racah–Speiser 0.01 + PRV spaces 0.09 + reduction 0.004 + **reduced operators 1.84** + eigen 0.56; t_b = 3.89 s with 0.09 s recursion setup and 0.64 s reconstruction of all 728 vectors (2.2·10⁶ coefficients). Stress family a = 3 ([9,6,3,0]², dim 4096, 213 irreps, 1592 copies, m ≤ 27, w ≤ 32, n_W ≤ 59 008): t_a = 23.4 s = **operators 5.9 + 50-digit eigen 17.0** + generators 0.39 + PRV 0.08; with float64 eigen: t_a,f = 6.4 s, t_b,f = 6.9 s (reconstruction 0.43 s). a = 4 ([12,8,4,0]², dim 15 625, 459 irreps, **5945 copies, m ≤ 51**, w ≤ 65, n_W ≤ 426 425; one run): t_a = 258 s = operators 34.4 s + 50-digit eigen 221 s; with float64 eigen t_a,f = 37.2 s. a = 7, 8 of the scaling family (dim 13 056 / 23 085, 484 / 693 irreps, 1212 / 1905 copies): t_a = 6.4 / 12.0 s. [S,S,0,0]² (m = 1 everywhere: the operators are not needed): S = 4: t_a = 0.017 s, t_b = 0.022 s; S = 8 (dim 825, 45 irreps): t_a = 0.083 s, t_b = 0.148 s.

**Where the reduction cost sits.**
* **Not in the reduction itself**: normal ordering + evaluation at λ′ + trie construction is ≈ 4–5 ms per λ′ for all three operators and does not grow with the representation; the universal polynomials cost 0.3–1.2 ms.
* **In the application of the reduced polynomials to the PRV basis**: t_ops is 56–62 % of t_a on the scaling family and ≈ 25 % of the stress family at a = 3 (the rest being the mpmath eigen-decomposition); it grows like copies × w (a = 8: 1905 copies, w ≤ 33: 6.7 s). With exact Fractions the three reduced operators now cost **more than the whole level recursion** (t_ops 1.84 s vs t_recon 0.64 s at [12,6,0,0]²: resolving the multiplicities costs ≈ 3× the reconstruction of all HW vectors at this size), i.e. the n_W-free part is "free" only relative to Munich's O(n_W³) (where the same family needs n_W up to 29 028), not relative to the recursion.
* **Against the gl(3) 24-term V_red** (same machine and protocol; ONE operator, exact application; µs per copy): gl(3) [28,14,0]² (3375-dim, max w = 15): 265 µs per copy (1.42 ms per irrep); gl(4) [10,5,0,0]² (max w = 15): 2.03 ms per copy for the triple = **678 µs per operator and copy** = 2.6× the gl(3) cost at the same w; at smaller w the factor is smaller (1.5× at w = 4–5, 2.2× at w = 11), it grows with the size of the weight space (more words per operator, 36/61/67 reduced words at λ′ = [2a,a,0,0] against 24, and longer words, up to 4 generators against 3):

* gl(3) [ 8, 4,0]^2: dim 125, irreps 41, copies 85, max w 5 | t_Vred 0.0139 s = 163 us per copy ( 338 us per irrep) | t_a 0.0194 s, t_b 0.0291 s
* gl(3) [14, 7,0]^2: dim 512, irreps 113, copies 344, max w 8 | t_Vred 0.0739 s = 215 us per copy ( 654 us per irrep) | t_a 0.1002 s, t_b 0.1621 s
* gl(3) [20,10,0]^2: dim 1331, irreps 221, copies 891, max w 11 | t_Vred 0.2071 s = 232 us per copy ( 937 us per irrep) | t_a 0.2872 s, t_b 0.5081 s
* gl(3) [28,14,0]^2: dim 3375, irreps 421, copies 2255, max w 15 | t_Vred 0.5977 s = 265 us per copy (1420 us per irrep) | t_a 0.8806 s, t_b 1.7238 s
* gl(4) [ 2, 1,0,0]^2: dim 20, irreps 7, copies 8, max w 2 | t_ops 0.0029 s = 365 us per copy for the triple = 122 us per operator and copy ( 417 us per irrep) | t_a 0.0087 s (t_a_f 0.0082 s with float64 eig), t_b 0.0097 s
* gl(4) [ 4, 2,0,0]^2: dim 126, irreps 24, copies 32, max w 4 | t_ops 0.0236 s = 737 us per copy for the triple = 246 us per operator and copy ( 982 us per irrep) | t_a 0.0411 s (t_a_f 0.0353 s with float64 eig), t_b 0.0463 s
* gl(4) [ 6, 3,0,0]^2: dim 480, irreps 58, copies 90, max w 7 | t_ops 0.1032 s = 1146 us per copy for the triple = 382 us per operator and copy (1779 us per irrep) | t_a 0.1669 s (t_a_f 0.1418 s with float64 eig), t_b 0.1930 s
* gl(4) [ 8, 4,0,0]^2: dim 1375, irreps 115, copies 205, max w 11 | t_ops 0.3192 s = 1557 us per copy for the triple = 519 us per operator and copy (2775 us per irrep) | t_a 0.5176 s (t_a_f 0.4319 s with float64 eig), t_b 0.6112 s
* gl(4) [10, 5,0,0]^2: dim 3276, irreps 201, copies 406, max w 15 | t_ops 0.8253 s = 2033 us per copy for the triple = 678 us per operator and copy (4106 us per irrep) | t_a 1.3751 s (t_a_f 1.1468 s with float64 eig), t_b 1.6595 s
* gl(4) [12, 6,0,0]^2: dim 6860, irreps 322, copies 728, max w 20 | t_ops 1.8423 s = 2531 us per copy for the triple = 844 us per operator and copy (5721 us per irrep) | t_a 3.1602 s (t_a_f 2.6136 s with float64 eig), t_b 3.8851 s

* **The eigen-decomposition** (50 digits) dominates only for large multiplicities: m³-like growth, 17.0 of 23.4 s at [9,6,3,0]² (213 sectors, m ≤ 27) and 221 of 258 s at [12,8,4,0]² (m ≤ 51, 459 sectors); the float64 version costs 0.027 s and 0.14 s (|ev − ev_mp| ≤ 7.2e-16 / 2.1e-15 and off-diagonal ≤ 6.8e-12 / 7.0e-12 (relative) at m = 27 / 51).
* The float residuals of the reconstructed copies (G1's recursion): 9.8e-14 at [12,6,0,0]² (the float64 recursion's residual grows with the representation, G1 §4.4), 3.9e-14 at [9,6,3,0]², ≤ 1e-14 elsewhere.
* The total n_W-free cost of **all canonical copies** (operators + 50-digit eigen-decomposition + slices) for [12,6,0,0]² (728 copies, vectors with up to 29 028 entries, dim 6860) is 3.16 s, and with the full reconstruction of all copies t_b = 3.89 s; Munich's null-space method would need an n_W × n_W null space per irrep here (n_W up to 29 028, n_W³ = 2.4·10¹³ flops for the largest) and was not run (no cross-machine comparison is made here; G1 §6 has the Munich table for the multiplicity-free family).

**Multiplicity-free family [S,S,0,0]²** (S = 4, 8; m = 1 on every irrep, so the operators are not needed, but they are scalars): the scalars (V_red,3, V_red,4a, V_red,4b) on M_μ, all 15 irreps of [4,4,0,0]² (the 45 irreps of S = 8 are printed in the output, e.g. (0,0,0) on the top irrep (16,16,0,0), (144,0,0) on (16,15,1,0), (1440,0,7920) on the singlet (8,8,8,8)):

  {8,8,0,0}: (0, 0, 0); {8,7,1,0}: (40, 0, 0); {8,6,2,0}: (70, 0, 0); {8,5,3,0}: (90, 0, 0); {8,4,4,0}: (100, 0, 0); {7,7,1,1}: (90, 0, 180); {7,6,2,1}: (120, 0, 240); {7,5,3,1}: (140, 0, 280); {7,4,4,1}: (150, 0, 300); {6,6,2,2}: (160, 0, 432); {6,5,3,2}: (180, 0, 504); {6,4,4,2}: (190, 0, 540); {5,5,3,3}: (210, 0, 672); {5,4,4,3}: (220, 0, 720); {4,4,4,4}: (240, 0, 840)

**V_red,4a = 0 identically on [S,S,0,0]²**, and with C_r(ν) the Perelomov–Popov Casimirs of §6 and top = (2S,2S,0,0):

  V_red,3(μ) = ((S+1)/2)·(C₂(top) − C₂(μ)) = 4S(S+1)² − ((S+1)/2)·C₂(μ),    V_red,4b(μ) = F(μ) − F(top),  F = ((2S+1)/4)·C₂ − ¼C₄ + ⅛C₂²,

checked on **S = 1..10 (285 irreps, 0 mismatches)**; for λ = λ′ and m = 1 these are the exchange relations of §8.1 (V₃ = c/2 and V₄ᵇ = c′, the closed forms of §6 specialised to the rectangle). Timing of the two multiplicity-free runs and the closed-form check (output lines):

* S = 4: dim V 105, irreps 15, max w 3, max nW 147 | t_a 0.0171 s (gen 0.0027, gen_na 0.0022, rs 0.0003, prv 0.0002, red 0.0040, ops 0.0073, eig 0.0005), t_b 0.0225 s | worst HW resid 8.2e-16
* S = 8: dim V 825, irreps 45, max w 5, max nW 1749 | t_a 0.0834 s (gen 0.0223, gen_na 0.0247, rs 0.0018, prv 0.0011, red 0.0040, ops 0.0278, eig 0.0016), t_b 0.1480 s | worst HW resid 5.6e-15
* closed forms on [S,S,0,0]^2 (m = 1): V_red,3 = (S+1)/2 (C2(top) - C2(mu)), V_red,4a = 0, V_red,4b = F(mu) - F(top), F = (2S+1)/4 C2 - C4/4 + C2^2/8 (Perelomov-Popov chain Casimirs, top = (2S,2S,0,0)): checked on S = 1..10: 285 irreps, 0 mismatches


## 10. G3.6 Reference files `ref/gl4-cg-copies-*.txt` and their independent re-verification

Six files for [2,1,0,0]², [2,1,1,0]², [3,2,1,0]², [2,2,1,0]⊗[2,1,1,0], [4,2,0,0]², [6,4,2,0]² (slices only for the last, and for [3,2,1,0]², [4,2,0,0]², whose dim·dim′ exceeds 2000). TAB-separated records with a `#` header (conventions, the reduction rule, the operator definitions with G2's normalisation, the ordering and normalisation of the copies, script, date): `REP REP2 DIM DIM2 GTBASIS GTBASIS2 WEIGHTS1 WEIGHTS2 GMET1 GMET2 NSECTORS NCOPIES GENERIC`, then per sector (descending μ) `SECTOR μ m w nW`, `WSPACE`, `PRVBASIS`×m (exact, as in the G1 files), `PRVFREE` (the free positions = PRV coordinates), `PRVGRAM` (exact G⊗G Gram matrix of the full HW vectors of the PRV basis; GNORM2 of a copy = uᵀ·Gram·u), `VRED3`, `VRED4A`, `VRED4B` (exact m × m in the PRV basis, column j = V_red b_j), `CHARPOLY3` (exact characteristic polynomial of VRED3, ascending powers) and `CHARPOLYJ` (of the generic combination V₃ + 7/19 V₄ᵃ + 11/13 V₄ᵇ used for the eigen-decomposition, stated in the header and in `GENERIC`), and per copy k = 1..m `COPY μ k ev3 ev4a ev4b {{ν,val},…}` (30 digits; slice normalised to first nonzero component 1; copies ordered lexicographically ascending by (ev₃, ev₄ᵃ, ev₄ᵇ)), `HWRESID` (float level recursion, all E_ij), `GNORM2` (30 digits), and `VEC` (30 digits; where dim·dim′ ≤ 2000: three files).

| file | sectors | copies | VEC records | max HWRESID | worst joint-eigenbasis off-diagonal | bytes | time [s] |
|---|---|---|---|---|---|---|---|
| `gl4-cg-copies-2100.txt` | 7 | 8 | 8 | 3.47e-16 | 2.71e-51 | 12895 | 0.0 |
| `gl4-cg-copies-2110.txt` | 6 | 7 | 7 | 2.50e-15 | 7.88e-49 | 11992 | 0.0 |
| `gl4-cg-copies-3210.txt` | 16 | 30 | 0 | 2.39e-15 | 1.77e-48 | 26028 | 0.1 |
| `gl4-cg-copies-2210-2110.txt` | 6 | 7 | 7 | 2.63e-15 | 5.56e-50 | 12167 | 0.0 |
| `gl4-cg-copies-4200.txt` | 24 | 32 | 0 | 1.72e-15 | 1.4e-49 | 31805 | 0.1 |
| `gl4-cg-copies-6420.txt` | 76 | 297 | 0 | 9.89e-15 | 1.16e-44 | 309001 | 9.2 |

**Independent re-verification** (`verify_copies_ref`: parses the text, nothing of the reduction code is used; the G1/G2 builders are):
(1) GTBASIS/WEIGHTS/GMET/DIM equal `glN_reps`; sector set, m = Racah–Speiser, w = |WSPACE| = the weight space computed from WEIGHTS1, nW;
(2) PRVBASIS: exact string conditions E_i^{k_i}v = 0, rank m, canonical structure at PRVFREE, and the dimension of the exact nullspace of the string-condition matrix = m;
(3) VRED3/4A/4B against the **full two-site operators of G2** (term lists of `gl4-gaudin-ops.txt`) applied exactly to the exact HW vectors of the PRV basis (G1's exact recursion; HW verified exactly): the slice read-off Σ_i A_ij b_i (exact), PRVGRAM against the exact Gram of those vectors, CHARPOLY3/CHARPOLYJ recomputed exactly, the three matrices commute exactly and are G-self-adjoint (exact);
(4) COPY: slice in span(PRVBASIS), the eigen-equation A_k u = ev_k u (u = slice at PRVFREE) in 60-digit arithmetic, first nonzero component exactly 1, lexicographic order, G-orthogonality of the copies and GNORM2 from PRVGRAM;
(5) the FLOAT full vector of every copy from the level recursion: HW residual of all E_ij ≤ 1e-12, GNORM2 against the float norm, VEC against the recursion, and the **full** operators V₃, V₄ᵃ, V₄ᵇ of G2 (sparse float) applied to it: ‖OΨ − evΨ‖/(‖Ψ‖(1+|ev|)) ≤ 1e-12. Acceptance for every check: exact checks exact, 30-digit checks 1e-24…1e-27, float checks 1e-12.

| file | sectors | copies | TOTAL MISMATCHES | copy-in-span | eigen-equation (30 digits) | G-orthogonality (30 digits) | HW residual (float) | full-operator eigen-equation | GNORM2 vs float | VEC vs recursion | time [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `gl4-cg-copies-2100.txt` | 7 | 8 | **0** | 0.0e+00 | 6.8e-31 | 2.7e-31 | 3.5e-16 | 1.4e-15 | 5.9e-16 | 5.9e-16 | 0.1 |
| `gl4-cg-copies-2110.txt` | 6 | 7 | **0** | 6.7e-30 | 4.1e-30 | 3.0e-31 | 2.5e-15 | 1.7e-14 | 9.9e-16 | 7.8e-16 | 0.1 |
| `gl4-cg-copies-3210.txt` | 16 | 30 | **0** | 1.5e-29 | 4.3e-30 | 9.8e-31 | 2.4e-15 | 4.4e-15 | 1.1e-15 | 0.0e+00 | 0.4 |
| `gl4-cg-copies-2210-2110.txt` | 6 | 7 | **0** | 0.0e+00 | 1.1e-30 | 4.2e-32 | 2.6e-15 | 2.4e-14 | 1.4e-15 | 1.5e-15 | 0.1 |
| `gl4-cg-copies-4200.txt` | 24 | 32 | **0** | 1.0e-29 | 1.5e-30 | 9.5e-31 | 1.7e-15 | 1.1e-14 | 9.6e-16 | 0.0e+00 | 0.5 |
| `gl4-cg-copies-6420.txt` | 76 | 297 | **0** | 4.7e-28 | 9.9e-29 | 3.4e-30 | 9.9e-15 | 1.2e-13 | 2.7e-15 | 0.0e+00 | 13.6 |

(copy-in-span / eigen-equation / G-orthogonality are measured in 60-digit arithmetic on the 30-digit decimals of the file, hence ≈ 1e-28…1e-31; the HW residual, full-operator eigen-equation, GNORM2 and VEC columns are float64 quantities.)

**Consistency with G1's files** (`ref/gl4-cg-slices-*.txt`): SECTOR fields, WSPACE and PRVBASIS identical (exact), and on every m = 1 sector the COPY slice (30 digits vs G1's exact rationals), the GNORM2 and, where both have it, the full vector VEC agree:

* [2,1,0,0]^2 vs gl4-cg-slices-2100.txt: 7 sectors (SECTOR fields, WSPACE, PRVBASIS exact), 6 m = 1 copies (slice, GNORM2), 6 full vectors: 0 mismatches
* [2,1,1,0]^2 vs gl4-cg-slices-2110.txt: 6 sectors (SECTOR fields, WSPACE, PRVBASIS exact), 5 m = 1 copies (slice, GNORM2), 5 full vectors: 0 mismatches
* [3,2,1,0]^2 vs gl4-cg-slices-3210.txt: 16 sectors (SECTOR fields, WSPACE, PRVBASIS exact), 8 m = 1 copies (slice, GNORM2), 0 full vectors: 0 mismatches
* [2,2,1,0]x[2,1,1,0] vs gl4-cg-slices-2210-2110.txt: 6 sectors (SECTOR fields, WSPACE, PRVBASIS exact), 5 m = 1 copies (slice, GNORM2), 5 full vectors: 0 mismatches
* [4,2,0,0]^2 no G1 file (gl4-cg-slices-4200.txt) to compare with
* [6,4,2,0]^2 no G1 file (gl4-cg-slices-6420.txt) to compare with

**Mutation test of the verifier** (a verifier that cannot fail proves nothing): ONE record at a time of a copy of `gl4-cg-copies-2100.txt` is corrupted (in the m = 2 sector (3,2,1,0), or in the header) and `verify_copies_ref` is run; every corruption must be detected, the pristine copy must pass:

| corruption of ONE record | failed checks | verdict |
|---|---|---|
| pristine copy | 0 | OK (passes) |
| VRED3 entry changed | 6 | DETECTED |
| VRED4A entry changed | 6 | DETECTED |
| VRED4B entry changed | 6 | DETECTED |
| CHARPOLY3 coefficient changed | 1 | DETECTED |
| CHARPOLYJ coefficient changed | 1 | DETECTED |
| PRVGRAM entry changed | 5 | DETECTED |
| PRVBASIS value changed | 6 | DETECTED |
| PRVFREE index changed | 8 | DETECTED |
| COPY eigenvalue ev3 changed by 1e-6 | 2 | DETECTED |
| COPY eigenvalue ev4b changed by 1e-6 | 2 | DETECTED |
| COPY slice value changed by 1e-6 | 6 | DETECTED |
| GNORM2 changed by 1e-6 (relative) | 2 | DETECTED |
| VEC value changed by 1e-6 | 1 | DETECTED |
| SECTOR multiplicity changed | 1 | DETECTED (verifier raised IndexError) |
| WSPACE index changed | 1 | DETECTED (verifier raised KeyError) |
| GMET1 entry changed | 1 | DETECTED |
| GENERIC changed | 1 | DETECTED |


## 11. Surprises, pitfalls, inaccuracies in the brief

* **Everything agreed on the first try**: no failing test, no discrepancy with G2 or G1 anywhere (46 + 598 + 2275 + 449 exact direct checks, 6 reference files). The checks that could have failed were made sensitive on purpose (§3: deliberately wrong sign / order / ordering rules; §10: mutation test of the file verifier).
* **The normal ordering is not needed for the cdet term lists** (§3): the column order of the cdet puts raising factors only to the left of lowering factors they commute with, for θ₁ and θ₂ data alike. It is needed (and tested) for any other representative of the same operators, e.g. G2's antisymmetrised V₄ᵇ. The pitfall as written in the brief has the roles of e²_{12}e²_{21} / e²_{21}e²_{12} exchanged (the word with the raising factor on the left is the one that contributes).
* **The reduction is much smaller than its input**: 324 → 85 universal words for V₄ᵇ, 88 → 39 for V₃; and the reduced polynomials are *universal in λ′*.
* **The reduced operators commute on the whole weight space and obey the BGG intertwining relation** (§7), which was not asked for and explains the invariance of M_μ; they are the Gaudin operators on the singular vectors of V_λ ⊗ M(λ′).
* **The scalars c, c′ of G2's relations have closed forms in the Casimirs** (§6), and equal the site-exchange scalars of G2 §6e.
* **The cost has moved**: the operators are n_W-free, but with exact Fractions the three reduced operators now cost more than the whole level recursion (§9: t_ops = 1.84 s vs t_recon = 0.64 s at [12,6,0,0]²), and for large multiplicities the 50-digit mpmath eigen-decomposition dominates (17.0 s of 23.4 s at [9,6,3,0]², 221 of 258 s at [12,8,4,0]²); the float64 eigen-decomposition is accurate to 2.1e-15 in the eigenvalues and 40× (m ≤ 7) to 1600× (m ≤ 51) cheaper.
* **G2's explanation of the V₃ degeneracy is exact but incomplete**: |p₊−p₋| ≥ 2 forces it, a coincident 2-cycle (V₃ = c/2 on both partners) produces additional accidental ones (§8.1).
* **Tooling trap** (G1, G2): `\n`, `\t` in Python source written through a Bash heredoc are un-escaped by the tool; I hit it twice in patches of the module (each caught at once by an `ast.parse` and repaired with the file tools) and twice in scratch scripts; every file of the package was written with the file tools.
* **Inaccuracies in the brief (harmless)**: the pitfall example above; "V₃ … forced degenerate … e.g. μ = (5,4,2,1)": the (5,4,2,1) degeneracy of [3,2,1,0]² is accidental (§8.1), the forced one is (4,4,2,2); otherwise no inaccuracy affected the work.

## 12. Not done / limits

* **No proof** of the closed forms of c, c′ (§6), of the rectangular-family scalars (§9), nor of (a)–(c) of §7: they are exact verifications on all tested cases (the BGG intertwining and the Verma-module interpretation are natural but not proved here).
* The reduced operators are applied in **exact Fractions** (Python); no float or modular fast path for the operators themselves (only for the joint eigen-decomposition, §9). The benchmark therefore reports the honest cost of the exact n_W-free pipeline; a float64 version of the operator application would be ≈ 10–100× faster but was not written.
* No gl(N), N ≠ 3, 4 beyond the normal-ordering engine (general N) and the N = 3 test; λ₄ > 0 only in the listed products; products with m up to 51 only for the stress family (slices only).
* The comparison with Munich's HW vectors remains deferred (G1 §9); the Gaudin eigenvalue ↔ oper dictionary (G5) and the structure of the copies' entries (G4) are not part of this package.
* The reference files carry VEC only for dim·dim′ ≤ 2000 (three small products), as asked; the others reconstruct by linearity from the PRVBASIS exact vectors (G1) and the copies' coordinates.

## 13. Usage

```python
import sys; sys.path.insert(0, 'Paul/Python/gaudin_gl4')
import glN_reps as gr, slice_engine_gl4 as se, slice_ops_gl4 as so
lam = lamp = (3, 2, 1, 0); mu = (5, 4, 2, 1)
prv, mats = so.triple_mats(lam, lamp, mu)          # prv = G1's PRV; mats = {'V3','V4a','V4b'}: exact m x m (rows of Fractions), column j = V_red b_j in PRV coordinates
A  = so.V_red(lam, lamp, mu, 'V3')                 # one operator; 'V4a','V4b','S31','S41','S_k_j','S_k_j_sw' (theta_2); raises ValueError if an image leaves M_mu
cps = so.copies(lam, lamp, mu, dps=50)             # canonical copies: cp['ev'] = (ev3, ev4a, ev4b) mpf; cp['slice'] = {global GT index: mpf} (first nonzero component 1); cp['coords'] = PRV coordinates; cps.info = diagnostics
F  = so.copy_float_slices(prv, cps)                # (w, m) float slices -> se.hw_from_slice(lam, lamp, mu, F) = the full vectors (float recursion), hw.vals (plain GT basis)
red = so.reduce_terms(so.op_terms('V4b'), lamp)    # {site-1 word: Fraction} at lam'; so.reduce_terms_poly(so.op_terms('V4b'), 4) = the universal polynomial in lam'
so.normal_order(((1, 2), (2, 1)))                  # {((2,1),(1,2)): 1, ((1,1),): 1, ((2,2),): -1}
so.write_copies_ref(lam, lamp); so.verify_copies_ref('Paul/Mathematica/Gaudin/ref/gl4-cg-copies-3210.txt')
```
Reading a copies file: `so.parse_copies_ref(path)` → sectors[mu]['copies'][k]['ev' | 'slice' | 'gnorm2' | 'vec'] as exact Fractions of the 30-digit decimals.
