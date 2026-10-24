# RESULTS G8c: the 24 relabelled Gelfand–Tsetlin chains of gl(4) — the triple and the strings in every σ-basis (2026-10-24)

Code `relabelled_gl4.py` (CLI `check | bandwidth | cuts | small | parity | all`, `run <stage,...>`), tests `test_relabelled_gl4.py` (6 exact tests, 5 s), tables `relabelled_gl4_tables.txt`,
full printed output `relabelled_gl4_output.txt` (the run is ≈ 60 s cold; the per-product record caches `relabelled_gl4_rec_*.pkl` are rebuildable and gitignored).
Everything is exact (Fractions, exact charpolys via Faddeev–LeVerrier, sympy for factorisations/radicals) except where "50 digits" is said (only the nested-basis graphs of C3).

## 0. Headline

**Verdict: the 24 relabelled chains do not rescue the Jacobi structure of gl(4).** The path property of the triple is a property of the *sector*, not of the chain.

1. **Per-σ path table** (V₃ sparsity graph a path; the union graph of the whole triple is a path in exactly the same spaces; 221 weight spaces with w ≥ 2 of P_small ∪ P_big, 63 with m = w, 158 with m < w):

   | σ-class | paths, m = w | paths, m < w | mean RCM bandwidth |
   |---|---|---|---|
   | id (1234) | 41 / 63 | 30 / 158 | 2.85 |
   | w₀ (4321) | 41 / 63 | 30 / 158 | 2.84 |
   | 2431, 4231 | 41 / 63 | **31** / 158 | 2.81 |
   | the other 20 σ | 41 / 63 | 30 / 158 | 2.81 – 2.85 |

   The path property is the same in all 24 bases in **220 of 221** spaces (one exception: (5,3,3,1) of [4,2,0,0]², w = 3, m = 1, a triangle at id that is a path at σ = 2431, 4231). "Some σ makes V₃ a path": 41/63 (m = w) and 31/158 (m < w) — i.e. **0 new spaces with m = w, 1 new with m < w**. Per class of active strings A: A = ∅ 41/63; {1} 6/22; {2} 10/45; {3} 14/42; |A| ≥ 2: 0/49 except one ({2,3}). The K₄ example is K₄ (6 edges, bw 3) in all 24 bases. The best σ only lowers the bandwidth mildly: min over σ of mean bw 2.56 vs 2.85 at id (m = w: 1.73 vs 1.84; m < w: 2.89 vs 3.25); strictly smaller at some σ in 60 of 221 spaces (7 of 63, 53 of 158); max bandwidth unchanged (8). (σ = id reproduces G4's 29/44 and 19/106, mean bw 1.84/3.34, exactly, on G4's nine products.)
2. **No rule "active strings → σ" exists for the path property** (it is σ-independent). The only exact rule is the **cut rule**: σ and σ∘(12) give identical graphs (0 differences in 221 × 24), so there are 12 inequivalent chains, and the string E_{i,i+1}^{k_i} is a **coordinate cut in the σ-basis iff {σ₁, σ₂} = {i, i+1}** (51/51, 76/76, 83/83 non-void spaces for i = 1, 2, 3; 0 failures; the converse holds except 120 accidental cases, all with w ≤ 3 and a 1-dimensional kernel). The cut is the same inequality for every such σ: **ker E_{i,i+1}^{k_i} = {σ-patterns with m₁₂ ≤ μ_i − λ′_{i+1}}** (m₁₂ = highest weight of the gl({i,i+1}) irrep; equivalently m₂₂ ≥ μ_{i+1} − λ′_i), 840 (space, string, σ) cases, 0 failures; for σ = w₀ the E₃₄ string is E₂₁ with kernel {m₂₂ > ν₄ − k₃}.
3. **Two active strings are never both cuts in one σ** (49 spaces with |A| ≥ 2; 0 counter-examples, also allowing accidental kernels): as in gl(3), the second string stays spectral. For one active string, M_μ = {m₁₂ ≤ thr} is a coordinate sub-lattice of W^σ in the cut σ (dim = m, triple-invariant, 22/22, 45/45, 42/42), and the triple restricted to it is a path in 7/22 (A = {1}), 19/45 ({2}), 19/42 ({3}) spaces — versus 6, 10, 14 for the whole W at id (M path ⊇ W path, never the converse: 0 cases W path, M not). So the cut helps (it removes spurious points) but M_μ is still not a chain in ≈ 60 % of the one-string sectors.
4. **Stencil**: in every σ the support of V₃ is the same 20 off-diagonal displacements ({0,±1} × ({0} ∪ gl(3)-roots) minus 0), nothing new.
5. **K₄ anatomy** ((5,4,2,1) of [3,2,1,0]², m = w = 4): lattice (m₁₂,m₁₃,m₂₃) = (2,2,2), (3,3,1), (2,3,2), (3,3,2); V₃ charpoly (x−57)²(x²−114x+3141); the generic copy polynomial factors over ℚ as **2 + 2**: a rational block A (V₃ = 57·𝟙, V₄ᵃ+V₄ᵇ = 72·𝟙, joint eigenvalues (57, 36 ∓ 6√6, 36 ± 6√6), components in ℚ(√6)) and a rational block B (joint eigenvalues (57 ∓ 6√3, 45 ∓ 15√3, 9 ± 15√3), components in ℚ(√3)); the exchange P swaps the two eigenlines inside each block. In the within-path (κ-eigen) basis the graph is still K₄-like (5 edges for V₃, 6 for the triple). **No σ, and not the nested basis, makes it a path.**
6. **Parity**: no monomial lattice involution realises the exchange in any σ-basis for m ≥ 3 (0 of the 12 sectors with m ≥ 3 of [3,2,1,0]², [4,2,1,0]² that have a σ with coordinate M^σ; only m = 2 sectors admit one, trivially) — §6.
7. Side finding (§5): for λ = λ′ the copy polynomial splits over ℚ in about half the sectors with m ≥ 2, for λ ≠ λ′ in none (0/20).

## 1. Conventions (C0)

σ = (σ₁,σ₂,σ₃,σ₄) ∈ S₄ is the chain gl(4) ⊃ gl({σ₁σ₂σ₃}) ⊃ gl({σ₁σ₂}) ⊃ gl({σ₁}) of the **original index labels**; the σ-GT basis b_p of V_λ is the GT basis adapted to it, i.e. E_{σ_aσ_b} acts on b_p by the standard GT matrix of E_ab. Hence with τ = σ⁻¹:

* the matrix of π(E_cd) in the σ-basis = the matrix of E_{τ(c)τ(d)} in the standard basis; a word X relabels letter by letter (the `red` dict keys are tuples of (i,j) pairs, so this is a dict comprehension; ReducedOp is rebuilt);
* the σ-pattern has E_aa-eigenvalue ν_{σ(a)}, i.e. it lives in `R.wsp[(ν_{σ1},…,ν_{σ4})]` (ν = μ − λ′);
* the string E_{i,i+1}^{k_i} becomes E_{τ(i)τ(i+1)}^{k_i} (a lowering operator when τ(i) > τ(i+1); `power_images` works for any (a,b));
* the brief's "E_ab → E_{σ(a)σ(b)}" is the same thing for the representation π∘φ_σ; for involutions σ = τ.

**Check (exact, all 16 sectors of [3,2,1,0]² × all 24 σ = 384 cases)**: dim M^σ = m, M^σ invariant under the triple, the exact charpolys of V₃, V₄ᵃ, V₄ᵇ and of the generic combination V₃ + 7/19 V₄ᵃ + 11/13 V₄ᵇ on M^σ equal those on M_μ (identity chain), and the restricted triple commutes: **384/384 ok**. (The brief asks for 50 digits; charpolys of a combination with simple joint spectrum are a stronger exact statement — the 50-digit jspec was not run.) The **wrong** pairing (generators relabelled by τ = σ⁻¹ but weight permuted by σ⁻¹ instead of σ) is accepted for the 10 involutions (160/160) and fails for 208 of the 224 non-involution cases (16 pass trivially: M^σ = W^σ of dimension 1 or equal spaces); the pairing (σ, σ⁻¹) with the roles exchanged is not wrong, it is the chain of σ⁻¹. σ = id reproduces G3's weight-space matrices, G1's PRV basis and G3's V_red on M_μ exactly (test), and G4's sparsity table (§2).
σ = w₀: E₃₄^{k₃} → E₂₁^{k₃}, kernel = {m₂₂ > ν₄ − k₃} = {m₁₂ ≤ μ₃ − λ′₄} exactly on all 29 non-void sectors of [3,2,1,0]², [4,2,1,0]², [6,4,2,0]².

Products: P_small (9 incl. both orders of [3,2,1,0]⊗[4,2,1,0]: 102 spaces with w ≥ 2) and P_big = [6,4,2,0]² (62 spaces), [5,3,1,0]² (39), [5,3,2,1]² (18): 221 spaces in all.

## 2. C1: bandwidth and chain structure per σ

Full per-σ table (all 24 rows; `relabelled_gl4_output.txt`): paths m = w always 41/63; m < w 30/158 except the pair {2431, 4231} with 31/158 (the exception above); mean bandwidth of V₃ between 2.81 and 2.85. σ and σ∘(12) have identical statistics: the gl(2)_{σ₁σ₂} Weyl reflection permutes (up to scale) the weight vectors of every gl(2)-irrep, so the two GT bases are monomially equivalent — 12 inequivalent chains (ordered pair (σ₄, σ₃) of the omitted index and the leftover index). The number of edges of the V₃ graph is σ-independent in 125/221 spaces; where it varies the path/non-path verdict still does not (apart from the one exception).

By class of active strings A (non-void strings) — spaces / V₃ path at id / path at some σ / mean bw at id / mean bw at the best σ per space:

| A | spaces | id | some σ | bw id | bw optimum |
|---|---|---|---|---|---|
| ∅ (m = w) | 63 | 41 | 41 | 1.84 | 1.73 |
| {1} | 22 | 6 | 6 | 3.14 | 2.64 |
| {2} | 45 | 10 | 10 | 2.89 | 2.69 |
| {3} | 42 | 14 | 14 | 2.60 | 2.43 |
| {1,2} | 8 | 0 | 0 | 4.00 | 3.50 |
| {1,3} | 18 | 0 | 0 | 4.17 | 3.83 |
| {2,3} | 20 | 0 | 1 | 4.00 | 3.10 |
| {1,2,3} | 3 | 0 | 0 | 6.33 | 5.67 |

Candidate rules tested (hits = spaces with a path, mean bw): "A = {i} → the cut σ of i": A = {1} → σ = 1234 (6 paths, bw 3.14 = id), {2} → 2314 (10, 2.89 = id), {3} → 3412 (14, 2.69 vs id 2.60, w₀ 2.71); "A = {3} → w₀" (14, 2.71), "A = {1} → w₀" (6, 2.91): none beats id on paths, none is better than the unconditional optimum on bandwidth; no σ is a path for all spaces of any class (the lists "all-σ-path sigmas" are empty). So there is **no rule A → σ** for path structure; the non-trivial content of the 24 chains is only the cut rule of §3.

The sub-lattice cut out by the coordinate strings of a σ (all non-void strings that are coordinate there) is always triple-invariant (0 failures over all σ) and has mean relative size 0.48–0.70 of W; its triple graph is a path in 13/55 (σ = 1234), 24/77 (2314), 32/89 (3421), … — see tables C2b.

## 3. C2: coordinate cuts

Per (σ, i): the coordinate classification over non-void spaces (counts coordinate/non-void): for the four σ with {σ₁,σ₂} = {i,i+1}: i = 1: 51/51, i = 2: 76/76, i = 3: 83/83 (the table of 24 × 3 cells is in the output); outside the prediction only 120 accidental coordinate kernels, all with w ∈ {2,3} and nullity 1. The cut inequality is stated in §0.2; E_{12}-type (σ₁ = i) and E_{21}-type (σ₁ = i+1) chains give the same threshold in m₁₂ because m₁₂ is intrinsic to the pair {i,i+1}. In terms of the BGG picture: E_i^{k_i}V(λ′) = V(s_i·λ′)E_i^{k_i} makes ker E_i^{k_i} triple-invariant, and in the cut σ it is the initial segment m₁₂ ≤ thr of every κ-fibre — the exact analogue of gl(3)'s "u = 0 at the cut".
**Two active strings**: |A| = 2: A = {1,2}: 8 spaces, {1,3}: 18, {2,3}: 20; |A| = 3: 3 spaces; in none is there a σ in which all active kernels are coordinate (with or without allowing 0-/full-dimensional kernels): the expectation (never) is confirmed, no counter-example.
**One active string** (table in §0.3): shape of M^σ = {m₁₂ ≤ thr}; triple restricted to it connected in 100 % and a path in 32 % (A={1}), 42 % ({2}), 45 % ({3}).

## 4. C3: the smallest non-chain sectors by hand

(Complete listings of patterns, exact 4×4 matrices, factorisations and the per-σ graph table: `relabelled_gl4_output.txt`, section C3.)

**(a) K₄: (5,4,2,1) of [3,2,1,0]²** (m = w = 4, A = ∅). Lattice points (coordinates (m₁₂,m₁₃,m₂₃); full patterns in the output): 1:(2,2,2), 2:(3,3,1), 3:(2,3,2), 4:(3,3,2). Exact V₃ (column j = image of point j):

```
V3 = [[381/8, -9/8, -11/8, 3/8], [-15/4, 231/4, 1/4, -9/4], [-45/8, 9/8, 459/8, -27/8], [15/4, -15/4, -9/4, 261/4]]
```
V₄ᵃ, V₄ᵇ in the output; all three have the full K₄ graph (6 edges, bw 3) at every σ (all 24: the graph of a 4-point space with every off-diagonal pair connected cannot change). Characteristic polynomials: V₃: (x−57)²(x²−114x+3141); V₄ᵃ: (x²−90x+1350)(x²−72x+1080); V₄ᵇ: (x²−72x+1080)(x²−18x−594). The generic combination has two irreducible quadratic factors over ℚ, hence two **rational invariant planes** (in lattice coordinates):

* block A = span{(−1/15, −2/3, 1, 0), (−1/5, 2, 0, 1)}: V₃ = 57·𝟙, V₄ᵃ = [[48,9],[8,24]], V₄ᵇ = [[24,−9],[−8,48]] (V₄ᵃ+V₄ᵇ = 72·𝟙); joint eigenvectors u₁ + t u₂ with t = −4/3 ∓ (2/3)√6, joint eigenvalues (57, 36 ∓ 6√6, 36 ± 6√6);
* block B = span{(7/3, 2/3, 1, 0), (1, 0, 0, 1)}: V₃ = [[45,−9],[4,69]], V₄ᵃ = [[15,−45/2],[10,75]], V₄ᵇ = [[39,45/2],[−10,−21]] (V₄ᵃ+V₄ᵇ = 54·𝟙); joint eigenvectors with t = −4/3 + (2/3)√3 (eigenvalues (57−6√3, 45−15√3, 9+15√3)) and t = −4/3 − (2/3)√3 (eigenvalues (57+6√3, 45+15√3, 9−15√3)).

So the four copies are two conjugate pairs, one over ℚ(√6) (the V₃-degenerate pair, V₃ = c/2 = 57) and one over ℚ(√3). The exchange P (λ = λ′; c = 114, c′ = 144 reproduced from the traces) acts through T(e₃, e₄ᵃ, e₄ᵇ) = (114−e₃, −2e₄ᵃ−e₄ᵇ+144, 3e₄ᵃ+2e₄ᵇ−144): exactly the Galois conjugation inside each block (it swaps 36∓6√6 in block A and 57∓6√3 in block B); this reproduces G3's 2-cycles (1,4),(2,3), and the V₃ double eigenvalue 57. Change of basis: the κ-groups (κ = (m₁₃,m₂₃)) are {1}, {2}, {3,4} (the single within-path pair 3,4 with κ = (3,2)); diagonalising V₃'s κ-block (50 digits, entries < 10⁻³⁰ counted as zero) leaves V₃ with 5 edges and the whole triple with 6 — the nested basis does **not** make K₄ a path; the rational 2+2 basis does block-diagonalise it, but that is the spectral decomposition, not a lattice statement.

**(b) the next smallest non-path m = w spaces** (all w = 4, K₄ graph, all three have λ ≠ λ′): (5,4,3,1) and (6,4,2,1) of [3,2,1,0]⊗[4,2,1,0], (6,4,2,1) of [4,2,1,0]⊗[3,2,1,0]. For all three: V₃ is K₄ in the 12 inequivalent chains except (6,4,2,1) of [4,2,1,0]⊗[3,2,1,0], where σ = 1324 gives bandwidth 2 (still not a path); for the other two the best bandwidth is 3 and no σ is a path; the generic copy polynomial is an **irreducible quartic over ℚ** (no rational blocks), listed in the output. Further K₄'s with λ = λ′: [4,2,1,0]² (6,5,2,1) (irreducible), (7,4,2,1) (irreducible), [4,3,1,0]² (7,6,2,1) (2 + 2, like the example). The smallest m < w cycle: (5,3,3,1) of [4,2,0,0]² (w = 3, m = 1, A = {2,3}): a triangle at id and at all σ except 2431/4231 where it is a path (bw 1); the single copy is the rational line spanned by the third lattice point (coordinates (2,4,1), pattern ((1),(2,0),(4,1,0))) with joint eigenvalues (V₃, V₄ᵃ, V₄ᵇ) = (80, 0, 160). Other m < w cycles of w = 3 (6 in the list): (5,5,1,1), (7,3,1,1) of [4,2,0,0]², etc.

**(c) description of the copies.** For λ = λ′ the exact structure is the Galois splitting of §5 (copies in blocks of degree 1, 2, …, with components in quadratic fields in the smallest cases) and the T-involution acting as conjugation inside blocks; for λ ≠ λ′ the copy polynomial is irreducible of degree m.

## 5. Galois block structure (extra; `small` stage)

Degrees of the irreducible factors over ℚ of the charpoly of the generic combination on M_μ (sectors with m ≥ 2): [3,2,1,0]²: 3 irreducible, 5 split; [4,2,0,0]²: 5 irreducible, 1 split, 1 totally rational; [4,2,1,0]²: 8, 6, 0; [4,2,2,0]²: 5, 1, 0; [4,3,1,0]²: 8, 9, 1; [5,3,1,0]²: 14, 17, 1; [5,3,2,1]²: 8, 6, 0; **[3,2,1,0]⊗[4,2,1,0] and the reverse: 10 irreducible, 0 split each**; [2,1,0,0]², [2,1,1,0]²: their one sector irreducible. Examples: degrees [2,2] at (5,4,2,1), (6,4,2,2) of [4,2,1,0]²; [1,2] for m = 3; [1,4], [1,6] (e.g. (8,6,3,1) of [5,3,1,0]²: m = 7 = 1 + 6). Splitting happens only for λ = λ′, presumably an effect of the exchange symmetry (rational structure of P and its fixed lines); this was not proved and not analysed further.

## 6. C4: exchange parity as a lattice symmetry

On the sectors of [3,2,1,0]² (8 with m ≥ 2) and [4,2,1,0]² (14) the question was: in a σ-basis in which M^σ is a coordinate sub-lattice, is there a monomial S (S e_a = d_a e_{π(a)}) with S V₃ S⁻¹ = −V₃ + c, S V₄ᵃ S⁻¹ = −2V₄ᵃ − V₄ᵇ + c′, S V₄ᵇ S⁻¹ = 3V₄ᵃ + 2V₄ᵇ − c′ (exhaustive backtracking over π, exact scales)? Result: **only for m = 2** (all the m = 2 sectors that have a σ with coordinate M^σ: 2 of 8 and 4 of 14 sectors counted over all m ≥ 2; every solution involutive — a 2 × 2 accident); for every one of the 12 sectors with m ≥ 3 that have a σ with coordinate M^σ (5 of 5 in [3,2,1,0]², 7 of 9 in [4,2,1,0]²; the K₄ sectors (5,4,2,1), (7,4,2,1), (6,5,2,1) have all 12 chains available) **no** monomial implementation exists. Sectors with ≥ 2 active strings have no σ with coordinate M^σ (§3), so are untestable here. So P is not a lattice reflection in any of the 12 chains (necessary condition; even a monomial S would equal P only up to a function of the triple). The true P (from the full HW vectors) was not computed.

## 7. Surprises and corrections

* The path property (and almost the whole statistics) is σ-independent: 220/221. I expected at least the A = {3}, A = {1} classes to profit from w₀ / the cut chain; they do not (A = {3}: 14 paths at every σ; bandwidth at the cut σ is worse than at id).
* σ ≃ σ∘(12) (monomial equivalence): only 12 distinct bases, not 24 (not in the brief).
* "Wrong convention" in the brief: the pair (E → E_{σ(a)σ(b)}, weight ν_{σ⁻¹}) is not wrong, it is the chain of σ⁻¹; only mixing the two (σ-relabelling with σ⁻¹-weights) fails, and only for non-involutions.
* G4's 29/44 and 19/106 are reproduced; P_small ∪ P_big gives 41/63 and 30/158 at id (different product set).
* The stencil is 20 off-diagonal points in all σ (the "21" counts the diagonal).
* K₄ (5,4,2,1) has an exact 2 + 2 rational block structure with fields ℚ(√6), ℚ(√3); the other three smallest K₄'s (λ ≠ λ′) are irreducible quartics.

## 8. Not done / limits

* The 50-digit joint-spectrum check of C0 was replaced by exact charpolys (+ exact commutation).
* C3(b): for the three further K₄'s only graphs, matrices and the irreducible quartic are given, no hand eigenvectors (degree 4 over ℚ); the m < w cycle sector is m = 1. No eigenvector tables for non-K₄ cycles beyond w = 3.
* Accidental coordinate kernels (120) are only classified (w ≤ 3, nullity 1), not derived; no proof that the path property is σ-independent (it is a measurement on 221 spaces); P_big were not used in the by-hand stage; the λ ≠ λ′ and λ₄ > 0 cases are covered only through [3,2,1,0]⊗[4,2,1,0] (both orders) and [5,3,2,1]².
* C4: the true exchange P was not computed (needs the full HW vectors of G3's `hw`); only the monomial-necessity test was done, and only in σ with coordinate M^σ.

## 9. Usage

```
cd Paul/Python/gaudin_gl4
python -X utf8 test_relabelled_gl4.py                 # 6 exact tests, ~5 s
python -X utf8 relabelled_gl4.py check                # C0
python -X utf8 relabelled_gl4.py run bandwidth,cuts   # C1, C2 (records cached in relabelled_gl4_rec_*.pkl)
python -X utf8 relabelled_gl4.py small                # C3 + Galois blocks;   parity  -> C4;   all -> everything
```
API: `sigma_mats(lam, lamp, mu, sigma)` returns W^σ (global GT indices), σ-patterns, coordinates (m₁₂,m₁₃,m₂₃), exact matrices of V₃, V₄ᵃ, V₄ᵇ, and per string i (k_i, generator, void?, exact kernel rows, coordinate?, zero columns); `kernel_basis`, `restrict_to` give M^σ and the restricted triple.
