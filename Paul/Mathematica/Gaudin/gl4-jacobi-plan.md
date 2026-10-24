# Plan G8: the Jacobi structure of gl(4) — nesting, lattice gradings, relabelled chains (written 2026-10-24)

Paul (2026-10-24): "The beautiful Jacobi structure we had for gl(3) fails for gl(4). Surely life is not so cruel. Maybe we need to do some
nesting, then each operator becomes Jacobi on a certain subspace and we build up eigenvectors recursively? Or maybe we need linear combinations
of generators? Set up a plan and implement with sonnet subagents; I really want something beautiful here too."

This plan turns that into three falsifiable structures, each with a decisive computation, and dispatches them as three parallel Sonnet work
packages G8a/G8b/G8c (briefs `Paul/Python/gaudin_gl4/SPEC-G8{a,b,c}-*.md`). Everything below builds on G1–G7 of `gl4-slice-engine-plan.md`
and on the reports `RESULTS-G3-slice-ops-gl4.md` §7, `RESULTS-G4-structure-gl4.md` (stencils, closed forms, chain structure, polynomiality),
and the gl(3) template `Paul/Python/gaudin_sl3/RESULTS-jacobi-cg.md`.

## 0. What "the gl(3) Jacobi structure" was, and exactly how it fails for gl(4)

gl(3): the weight space W = V_λ[μ−λ′] is a **path** in the single free GT coordinate m₁₂; the reduced vertex operator V_red is tridiagonal
on it with closed-form entries u, d, α (products of linear factors; the diagonal has one pole term K/ξ); the E₁₂ string condition is a
**coordinate cut** m₁₂ ≤ μ₁ − λ′₂ (so M_μ = the initial segment W₁, invariant because u vanishes at the cut) and the E₂₃ string is the same
cut in the w₀-reflected chain; the slice of a copy is the **forward three-term recurrence** p₀ = 1, p_{k+1} = ((V−α_k)p_k − u_{k−1}p_{k−1})/d_k
at its eigenvalue, and the eigenvalues are the roots of the continuant (= the Q-system polynomial). One operator, one coordinate, one
recurrence, one polynomial.

gl(4) (G3, G4, exact): W has **three** free coordinates (m₁₂, m₁₃, m₂₃) — a 3-parameter region of the GT lattice, w ≤ 15 on [6,4,2,0]², ≤ 32 on
[9,6,3,0]². The three commuting operators V₃, V₄ᵃ, V₄ᵇ share the same **21-point stencil** {Δm₁₂ ∈ {0,±1}} × {row-3 move 0 or a root of gl(3)}.
So: (i) within a fixed row 3 (κ = (m₁₃,m₂₃,m₃₃)) each operator IS tridiagonal in m₁₂ — the weight space is a disjoint union of gl(3)-type
paths P_κ, one per κ, and the operators are "tridiagonal within paths + hopping between paths that differ by a gl(3) root"; (ii) the sparsity
graph of V₃ is a path in only 29 of 44 spaces with m = w (smallest counter-example: (5,4,2,1) of [3,2,1,0]², m = w = 4, graph K₄); (iii) 12 of
the 20 off-diagonal classes (all with Δm₁₂ = ±1 AND a row-3 move, i.e. the hopping classes) have ONE product formula for all three operators,
M^{V₄ᵃ} = (y_k+3)·M^{V₃}, M^{V₄ᵇ} = −(y_a−y_i+2y_k−q₃+5)·M^{V₃} (k = the spectator row-3 entry), while the within-path classes (±1,0,0) and
the Δm₁₂ = 0 hopping classes carry an irreducible polynomial Q (degree 6 resp. 3); (iv) only E₁₂ is a coordinate cut, E₂₃ and E₃₄ are
spectral; (v) every slice component is a polynomial in (v₃,v₄ᵃ,v₄ᵇ) of the generic counting degree d_gen(m) ~ (6m)^{1/3} (Krylov from a cyclic
coordinate), but no transparent recurrence produces those polynomials; (vi) M_μ is a spectral subspace of W; BGG intertwining
E_i^{k_i} V_red(λ′) = V_red(s_i·λ′) E_i^{k_i}; the triple commutes on all of W for every λ′.

Facts (iii)–(iv) are the clues. The within-path parts are the complicated ones and the hopping parts are the simple ones, and the three
operators differ on the hopping classes by affine factors of the spectator entry. The 21-point stencil is the Minkowski product of the gl(2)
stencil {0,±1} (row 2) and the gl(3) root stencil (row 3) — the lattice is "a gl(3) weight lattice of paths".

## 1. Three candidate structures and their decisive tests

### H-A. Block-Jacobi with respect to a lattice grading → a forward recursion (G8a)
The multivariate analogue of a Jacobi matrix (Xu's three-term relation for orthogonal polynomials in several variables): a grading g of the
lattice such that each operator is block-tridiagonal with respect to the levels L_d = {x : g(x) = d} and the stacked inter-level blocks
[M_k(L_d, L_{d+1})]_{k=3,4a,4b} (3n_d × n_{d+1}) have full column rank at every level. Then the eigenvector equations at level d determine
ψ_{d+1} from ψ_d, ψ_{d−1}: **ψ is built level by level from its bottom-level component**, every component is a polynomial in the triple
(explicitly, by the recursion), and the unused equations are the "characteristic equations" whose common solutions are the joint spectrum
(the gl(4) continuant). For a stencil within {0,±1}³ the band-1 gradings are exactly g ∈ {±m₁₂, ±m₁₃, ±m₂₃, ±m₃₃}; band-2 gradings
(e.g. m₁₂+m₁₃) are allowed if the recursion is run greedily (equations at level d determine everything not yet known up to level d+2).
Side fact (derived from G4's universal formula, to be checked on data): for a fixed source x and the three Δm₁₂ = +1 root moves, the
3×3 matrix of (V₃, V₄ᵃ, V₄ᵇ)-entries has determinant Π_α c_α · (y₂−y₁)(y₂−y₃) ≠ 0 — the three operators are "Vandermonde-independent" on
the hopping classes, which is why full column rank is plausible.
Decisive data: for every weight space (W, the E₁₂-cut W₁, and M̃ of H-B) and every grading: band, level sizes, bottom-level size n₀, rank
of every stacked block; the recursion implemented and validated against G3's copies; the characteristic equations for n₀ = 1 sectors.
"Beautiful" outcome: a grading (or a short list depending on the active strings) with full rank everywhere and n₀ = 1 in a natural class of
sectors, so that slice = P_x(v₃,v₄ᵃ,v₄ᵇ)·1 by an explicit lattice recurrence with closed-form coefficients.

### H-B. Nesting through gl(3) ⊂ gl(4): the paths carry the gl(3) chains (G8b)
W = ⊕_κ P_κ with P_κ = the gl(3) weight space V_κ[ν₁₂₃] of the gl(3)-irrep κ ≺ λ (interlacing, |κ| = |λ| − ν₄). The GT basis of P_κ IS the gl(3)
GT basis of V_κ (the GT formulas for E₁₂, E₂₁, E₂₃, E₃₂ do not see row 4). Two exact questions:
(B1) Is the κ-diagonal block of V₃ (the within-path tridiagonal part) equal to a(κ)·V^{(3)}_red(κ; λ′₁₂₃) + b(κ)·1, with V^{(3)}_red the gl(3)
reduced vertex operator of the sector (κ, λ′₁₂₃, μ₁₂₃)? G4's closed forms make this possible: the within-path classes (±1,0,0) are
"linear factors × Q₆" with Q₆ a polynomial in (κ, λ) only, so the m₁₂-dependence could be exactly the gl(3) chain's. Same question for
V₄ᵃ, V₄ᵇ (with the 1- and 2-parameter ansätze a·V^{(3)} + b, and a·V^{(3)} + b + c·(gl(2) Casimir)).
(B2) The nested subspace M̃ = ker E₁₂^{k₁} ∩ ker E₂₃^{k₂} ⊂ W equals ⊕_κ M^{(3)}_κ, the gl(3) PRV slices of the κ-sectors (dim = Σ_κ m₃(κ)), is
invariant under the triple (BGG), contains M_μ as a spectral subspace, and carries a canonical **nested basis** |κ, s⟩ = the gl(3) Gaudin copies
of each κ (the gl(3) Jacobi-chain vectors). In this basis the gl(4) triple is "diagonal within κ (if B1) + hopping between κ and κ ± root".
What do the hopping blocks look like — rank, sparsity in (s, s′), closed form in the gl(3) data? The physical copies expanded in the nested
basis: c(κ, s) — any product structure?
(B3) Depth: the number of paths grows with n = ν₄ = μ₄ − λ′₄ (n = 0: one path = the gl(3) problem; n = 1: ≤ 3 paths). Solve the depth-1
sectors completely (three gl(3) chains + explicit hopping from the universal formula) — the first genuinely gl(4) case in closed form.
"Beautiful" outcome: B1 true with closed-form a(κ), b(κ), and the depth-1 sectors solved by "three gl(3) chains coupled by one explicit
product-formula hopping"; the general case as a recursion over κ with gl(3) chains as the units.

### H-C. A different Gelfand–Tsetlin chain, and the smallest counter-examples by hand (G8c)
The GT basis is adapted to gl(4) ⊃ gl(3)_{123} ⊃ gl(2)_{12}; relabelling the indices by σ ∈ S₄ gives 24 GT bases of V_λ (the σ-chain
gl(4) ⊃ gl(σ₁σ₂σ₃) ⊃ gl(σ₁σ₂)). The matrix of V_red in the σ-basis is the matrix of the relabelled word (E_ab ↦ E_{σ(a)σ(b)}) in the
standard basis on the permuted weight space — the gl(3) "chain 2" (σ = w₀) was exactly this. In the σ-basis the string E_{σ(i)σ(i+1)}^{k_i}
is a coordinate cut iff {σ(i), σ(i+1)} = {1, 2}. Questions: in how many weight spaces does SOME σ make V₃ (or the whole triple) a path /
reduce the bandwidth; is there a σ-basis in which the K₄ example (5,4,2,1) of [3,2,1,0]² is tridiagonal; which strings are cuts in which σ;
and the complete anatomy of the smallest non-path sectors (patterns, exact 4×4 matrices, joint eigenvectors, exchange parity).
"Beautiful" outcome: a rule "sector class (active strings) → σ-chain in which the triple is Jacobi / has the coordinate cuts", generalising
the two chains of gl(3).

### Why constant linear combinations of V₃, V₄ᵃ, V₄ᵇ cannot help, and what would
On every hopping class the three operators are c·(1, y_k+3, affine) with the spectator y_k varying along the class, so as functions on the
lattice the three restrictions are linearly independent: no constant combination kills a class or shrinks the stencil (polynomials of
degree ≥ 2 in the triple only enlarge it). Position-dependent combinations (coefficients depending on the GT labels) are not in the
commutative algebra. The useful "linear combinations" are therefore of a different kind: (a) combinations whose coefficients depend on the
SECTOR (λ, λ′, μ) only — G4's LP functionals ℓ_A for the end-node classes are of this kind and only order the spectrum; (b) the "spectator
weighting" V₄ᵃ ≈ (y_k+3)V₃ suggests writing V₄ᵃ = ½{V₃, Y} + … with Y a diagonal operator — a structural identity to test in G8b (B1/B2),
not a new generator; (c) changing the BASIS rather than the operators (H-B, H-C).

## 2. Work packages (parallel, Sonnet, Python only, ≈ 2.5 h each; never edit existing files; no git; no Wolfram)

| package | brief | owns | decisive outputs |
|---|---|---|---|
| G8a lattice recursion | `SPEC-G8a-lattice-recursion.md` | `lattice_recursion_gl4.py`, `test_lattice_recursion_gl4.py`, `RESULTS-G8a-lattice-recursion.md`, `lattice_recursion_gl4_{output,tables}.txt` | grading tables (band, n₀, rank at every level) on W / W₁ / M̃; `slice_from_triple` validated against `copies`; characteristic equations for n₀ = 1; separability rank of Ψ[κ][m₁₂] |
| G8b nested gl(3) | `SPEC-G8b-nested-gl3.md` | `nested_gl3_gl4.py`, `test_nested_gl3_gl4.py`, `RESULTS-G8b-nested-gl3.md`, `nested_gl3_gl4_{output,tables}.txt` | B1 verdict with a(κ), b(κ); M̃ = ⊕ M^{(3)}_κ; the triple in the nested basis (hopping-block structure); depth-1 sectors in closed form; c(κ,s) tables |
| G8c relabelled chains | `SPEC-G8c-relabelled-bases.md` | `relabelled_gl4.py`, `test_relabelled_gl4.py`, `RESULTS-G8c-relabelled-bases.md`, `relabelled_gl4_{output,tables}.txt` | per-σ bandwidth/path tables; cut table; the K₄ example and the next smallest cycle sectors in full; a σ-rule if one exists |

Shared test products: P_small = [2,1,0,0]², [2,1,1,0]², [3,2,1,0]², [4,2,0,0]², [4,2,1,0]², [3,2,1,0]⊗[4,2,1,0] (both orders), [4,2,2,0]²,
[4,3,1,0]²; P_big = [6,4,2,0]² (w ≤ 15, m ≤ 12), [5,3,1,0]², [5,3,2,1]² (λ₄ > 0). Every statement exact (Fractions / exact rref) unless "50 digits" is said.

## 3. Synthesis (after the three return)
A verdict table H-A / H-B / H-C (true / true on a class / false, with the sector classes), the single best "beautiful" formulation that
survives, and the follow-up G8d: closed forms and proofs for whatever survived (e.g. the κ-recursion coefficients, the characteristic
equations vs G5's oper conditions, the Q-side meaning), plus the gl(N) statement. Nothing goes into the notes or the paper until Paul has seen
the verdict table.

## 4. VERDICTS (2026-10-24, all three agents returned; 18 regression tests pass; reports `RESULTS-G8{a,b,c}-*.md`)

| structure | verdict | the numbers |
|---|---|---|
| H-A lattice recursion (G8a) | **TRUE** (per-sector grading; bottom point unique in 85 %) | forward-determined on 222/222 weight spaces of W₁ and of W with the best grading per space; n₀ = 1 on 188/222 (W₁), n₀ = 2 on 33, n₀ = 3 on 1; no single grading works everywhere (+m₁₂: 218/222); `slice_from_triple` = `copies` on 358/358 (≤ 2.9e-43), 0.3 ms/copy (×8 vs eigen); deg P_x = level of x = G4's Krylov degree (14/14); the unused equations generate the vanishing ideal of Spec(W₁) exactly (dim ℚ[v]/J = dim W₁, 14/14) — **the gl(4) continuant**; n₀ > 1: ψ = T(v)ψ₀, C(v)ψ₀ = 0 verified on 203/203 copies of all 34 such sectors, the 2×2 minors generate the vanishing ideal (5/5 symbolic) |
| H-B nesting gl(3) ⊂ gl(4) (G8b) | **TRUE for V₃, V₄ᵃ; V₄ᵇ needs one commutator** | J^{(3)}_κ = a(κ;λ)·V^{(3)}_red(κ; λ′₁₂₃) + b₃(ν,λ′)·1 and J^{(4a)}_κ = a₄ₐ(κ;λ)·V^{(3)} + b₄ₐ EXACTLY on 305/305 (weight space, κ) pairs + 369/369 held-out on [9,6,3,0] with generic non-dominant λ′; a depends only on (κ, λ), a(λ₁₂₃) = 1 (depth 0 = the gl(3) problem, 90/90); J^{(4b)}_κ = aV + e[ξ,V] + b + cξ − (e/4)ξ² (ξ = gl(2)₁₂ Casimir) 305/305 + 369/369; M̃ = ker E₁₂^{k₁} ∩ ker E₂₃^{k₂} = ⊕_κ (gl(3) PRV slices), dim = Σ m₃(κ), triple-invariant, M_μ ⊂ M̃, dim M̃ = m iff the E₃₄ string is void (322/322); nested basis: V₃, V₄ᵃ diagonal within κ, hopping blocks dense and of full rank — no factorisation, no end-point mechanism; depth-1 sectors solved by `depth1_copies` (100/100 sectors) but not smaller than G3's problem; the Δm₁₂ = ±1 hopping = `universal_F` (5214 tests) |
| H-C relabelled chains (G8c) | **FALSE** (the path property is a property of the sector) | V₃ is a path in 41/63 (m = w) and 30/158 (m < w) spaces in 20 of the 24 chains, 31/158 in the other 4 (one space); 12 inequivalent chains (σ ≃ σ∘(12)); same 20-displacement stencil in every chain; cut rule exact: E_{i,i+1}^{k_i} is a coordinate cut iff {σ₁,σ₂} = {i,i+1}, kernel m₁₂ ≤ μ_i − λ′_{i+1}; two active strings never both cuts (49/49); the exchange parity is no lattice involution; K₄ example: copy polynomial splits over ℚ as 2+2 by parity (ℚ(√6), ℚ(√3)), exchange = Galois conjugation inside each block (λ = λ′ only) |

**The statement that survives (verified, not proved):** on W = V_λ[μ−λ′] = ⊔_κ P_κ,
(1) *nesting*: on every path P_κ the cubic and the first quartic gl(4) label are the gl(3) missing-label operator of the sub-sector (κ, λ′₁₂₃, μ₁₂₃) up to a scale a(κ; λ) and a shift — the gl(3) Jacobi chains sit inside the gl(4) operators, with a = 1 on the top path;
(2) *coupling*: paths differing by a gl(3) root are coupled by G4's one product formula (Δm₁₂ = ±1) and a cubic (Δm₁₂ = 0);
(3) *recursion*: with the grading by m₁₂ (−m₁₂, or m₃₃/−m₁₃ when only E₃₄ is active) the triple is block-Jacobi with injective forward blocks; the slice is built level by level from one lattice point (85 %) or from a ≤ 3-dimensional bottom level, its components are polynomials in (v₃,v₄ᵃ,v₄ᵇ) of degree = level, and the unused equations are the gl(4) continuant — the Gröbner basis of the vanishing ideal of the joint spectrum of W₁.
What is NOT there: no basis in which the triple is a chain; no factorised slices; no small nested secular problem (the hopping is dense).

**Follow-up G8d (dispatched 2026-10-24, `SPEC-G8d-closed-forms.md`):** closed forms of a(κ;λ), a₄ₐ, b₃, b₄ₐ and the V₄ᵇ within-path parameters, their derivation from the word structure, and the recursion assembled entirely from closed forms (gl(3) chain × a + b on the paths, `universal_F` and the explicit Q's between them) with no representation matrices — the gl(4) analogue of `jacobi_cg`'s matrix-free core. Open after that: proofs (B1 by Wigner–Eckart on the index-4 words; determinacy of the m₁₂ grading per class), the relation of the characteristic ideal to G5's four oper conditions, and the notes.

### 4.1 G8d returned (2026-10-24 evening; `RESULTS-G8d-closed-forms.md`, `closed_forms_gl4.py`, 8 tests pass)

| item | result |
|---|---|
| a(κ; λ) | **a = 1 − R₁ − R₂ − R₃**, R_i = Π_{s=1..4}(y_i − Y_s) / Π_{j≠i}(y_i − y_j)(y_i − y_j + 1) — the GT norm of the neighbour κ + e_i (y_i = κ_i − i, Y_s = λ_s − s); one fitted parameter (c = −1) on 21 κ, 129/129 held-out κ on five other λ, 1503 further samples on 10 unseen λ: 0 mismatches; a = 1 at κ = λ₁₂₃ and at κ = (λ₂,λ₃,λ₄) |
| a₄ₐ(κ; λ) | **a₄ₐ = d + Σ_i R_i (2y_i − |κ| + 4)**, d = |λ| − |κ| (4 fitted parameters; same held-out counts); a₄ₐ = λ₄ at the top path |
| b₃, b₄ₐ | exact TABLES (degree-2 polynomials in (ν₁,ν₂,λ′) with 28 coefficient functions each of the form P(κ) + Σ_i R_i S(y_i, …); 183 / 487 terms; 2004 held-out samples, 0 mismatches), not a derivation; the ν-quadratic part of b₃ is −(d + ΣR_i(y_i+3))·[ν₁² + ν₁ν₂ + ν₂² + λ′₁ν₁ + λ′₂ν₂ − λ′₃(ν₁+ν₂)]; not Casimir scalars |
| V₄ᵇ | e = (1 − a)/8 = ΣR_i/8 exactly (378/378); a′, c, b′ are affine / quadratic / quartic in (ν, λ′) with unidentified κ-coefficients — **no closed form for V₄ᵇ yet** |
| mechanism (D3) | the words split α (indices in {1,2,3} ∪ E₄₄; 21/30/36) and β (one E_{i4}, one E_{4j}; 18/48/49), γ = ∅; on every path block α = (1, d)·V^{(3)} + scalar and β = (a − 1, a₄ₐ − d)·V^{(3)} + scalar: β passes through the neighbours κ + e_k and returns, which is where the R_k come from (read off, not proved) |
| D4 | `closed_form_mats` (V₃, V₄ᵃ from the labels alone) == G3's matrices on all 328 weight spaces + 20 random spaces of [9,6,3,0]² with generic λ′; [V₃, V₄ᵃ] = 0 on the whole W; matrices at λ = (30,20,10,0) (dim 1.77·10⁶) in 0.05 s; `slices_closed` = `copies` on 295/297 of [6,4,2,0]² (the 2 others share (v₃,v₄ᵃ): V₄ᵇ needed) and 1592/1592 of [9,6,3,0]²; the general-λ′ gl(3) chain in closed form (614 spaces); **no speed gain in exact arithmetic** (2.97 s vs G3's 1.65 s on [6,4,2,0]²: the eigenvalues still come from G3's PRV matrices and the recursion runs at 50 digits) |

**Verdict after G8d:** the nesting is now explicit — V₃ on a path is the gl(3) chain scaled by 1 minus the sum of the three neighbour norms; V₄ᵃ the same with a weight-dressed sum. Open: the V₄ᵇ closed form (a′, c, b′ and five quartic Q's), a derivation of b₃/b₄ₐ, proofs. The fast-CG algorithm built on this is `SPEC-G9-fast-cg-ubuntu.md` (float64 engine; physical slices from the w₁ × w₁ eigenproblem on the E₁₂-cut space with the string-kernel test; G1's level recursion for the full vectors; Munich and G3 as the same-machine benchmarks).
