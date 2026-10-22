# RESULTS-G1-slice-engine-gl4 (agent G1, 2026-10-22): `glN_reps.py` and the PRV slice engine for gl(N) — GT representations, slice spaces with N−1 string conditions, the level recursion, the SU(4) benchmark against Munich, reference files

Files (all new, in `Paul/Python/gaudin_gl4/`; nothing existing was edited, no git, no notebooks, no Wolfram):
`glN_reps.py` (exact sparse Gelfand–Tsetlin irreps of gl(N), any N; Racah–Speiser, Weyl, Casimir; the GT-basis side of the Munich bridge),
`test_glN_reps.py` (10 tests, pytest-style functions and `python -X utf8 test_glN_reps.py`),
`slice_engine_gl4.py` (all engine code; CLI `python -X utf8 slice_engine_gl4.py tests|prv|hw|exact|neg|brute|bench [Smax [ext]]|families|ref [all|l1 l2 l3 l4 [l1p l2p l3p l4p]]|verify|munich|munichcmp S m1 m2 m3 m4 file|n5|extras|all`),
`slice_engine_gl4_output.txt` (the printed output of one run of `all`, 7.4 min on this machine, followed by two appended re-runs of the timing stages, see §6; every table below is parsed from it),
and ten reference files `../../Mathematica/Gaudin/ref/gl4-cg-slices-{2100,2110,3210,2100-1100,2210-2110,2200,3300,4400,6600,8800}.txt` (250 kB in total, table in §8).
Python 3.13.13, numpy 2.5.2, scipy 1.18.1, sympy 1.14, no gmpy2; i7-12700KF (12 cores / 20 threads), 32 GB, Windows 11; everything single-process (numpy BLAS threads untouched; they matter only for the dense `eigh` of the brute-force check);
imports `gl3_reps` (`gt_patterns`, `_raise`, `_lower`, `_weight`) and, for the N = 3 regression only, `gaudin_sl3/slice_engine.py`, read-only. The Munich binary cannot be built here: the literal comparison of HW vectors is deferred (§9), but §9 shows what could already be settled without it.

## 0. Conventions actually used (and where I refined the brief)

* **One site** (as `GL3Rep`/`slice_engine.Rep`, generalised): GT pattern = tuple of N rows, `rows[0]` the singleton, `rows[N−1] = λ`; basis order = `gl3_reps.gt_patterns(λ)` (itertools.product order = ascending lexicographic in (row N−1, row N−2, …, row 1)); **HW state = LAST index (asserted in the constructor)**; Python indices 0-based, files 1-based. `E[(i,j)]` = list over columns of dicts `{row: Fraction}` (the `slice_engine.Rep` layout) for ALL N² pairs: the diagonal and the 2(N−1) adjacent matrices are built with the representation (targets outside the pattern set are dropped, never clipped), every other one **lazily** on first access by `E_ij = [E_{i,j−1}, E_{j−1,j}]`, `E_ji = [E_{j,j−1}, E_{j−1,i}]` (j > i+1). All 256 commutation relations hold exactly (§2), and along every commutator path.
* **Metric**: diagonal G with `E_ijᵀ G = G E_ji`, `G[HW] = 1`, exact (walk from the HW state through the simple raising operators; `check_metric()` verifies every matrix element of every off-diagonal pair). Orthonormal basis `ẽ_q = e_q/√g_q`: `Eh[r,q] = E[r,q]√(g_r/g_q)`, `Eh_ji = Eh_ijᵀ`; `float_ops()` is a lazy dict of CSR matrices (incl. the diagonal). The recursion runs in the orthonormal basis and `HW.vals` converts back to the plain GT basis.
* **Two sites**: Ψ = ΣΨ[ν,q] e_ν⊗e_q, **site 1 = V_λ (ν, slow), site 2 = V_λ′ (q, fast)**, flat index ν·dim₂ + q; **slice** = Ψ[·, q_top] ∈ V_λ[μ−λ′] (weight space of SITE 1; the string exponents λ′_i − λ′_{i+1} + 1 come from SITE 2). `prv_space(lam, lamp, mu)` takes the two representations in this order; a string condition whose image weight is not a weight of V_λ is void.
* **Level recursion** (any N): site-2 weights labelled by `a ∈ ℤ≥0^{N−1}` (`a_i = Σ_{k≤i}(λ′_k − w_k)`), processed in lexicographic order; one SVD left inverse of the stacked site-2 raising blocks `A_a` per `(λ′, a)`, independent of μ and of the slice. **Refinement of the brief**: the exact twin verifies TWO kinds of conditions — (T) the consistency rows `T·RHS = 0` of the overdetermined stacked systems, and (D) the *dangling* rows `E_i ψ_r = 0` for `r ∈ Q(b)` with `b + e_i` not a weight of site 2 (they are the PRV string conditions themselves and are not rows of any `A_a`). Together T + D are exactly the rows of `E_i^{tot}Ψ = 0`, i = 1..N−1. §4.3 shows that for gl(4) the T-rows alone do NOT imply D (the gl(3) exact twin has no D check).
* **Racah–Speiser**: mult(μ) = Σ over the weights w of the SECOND factor V_λ′, `m_λ′(w)·sgn(s)`, over the w for which λ + w + ρ is regular, μ = s(λ+w+ρ) − ρ, ρ = (N−1,…,0) (`slice_engine.racah_speiser` sums over the weights of the first factor; the two agree on every N = 3 product tested, incl. λ ≠ λ′, and the result is symmetric in the factors for all products). **Candidates**: the brief's `candidates(lam, lamp)` = dominant λ + wt(V_λ′); the slice space is non-empty exactly for μ ∈ λ′ + wt(V_λ), so the PRV test runs over the UNION of both candidate sets (`candidates(..., both=True)`), including all those of multiplicity 0 (M_μ must be {0}).
* **Munich**: SU(4) labels = our μ minus μ₄·(1,1,1,1) (`munich_label_to_ours`); Munich's pattern index Q(M) = our index + 1; its orthonormal GT basis = ours rescaled by 1/√g (the sign gauge is trivial, §9).
* **Slice normalisation in the files**: first nonzero component (ascending GT index) = 1; the full vector carries the same scale; `GNORM2` = exact `Σ g_ν g_q Ψ[ν,q]²`.
* **Residuals**: `hw_residual` = max over the simple raising operators (or ALL E_ij, i < j, with `extra=True`; all residual tables below use `extra=True` unless said otherwise) of ‖E^{tot}Ψ‖_F/‖Ψ‖_F in the orthonormal basis (Frobenius norm = invariant norm), plus the max deviation of the total weight from μ over all entries (always 0).
* **Timing**: `time.perf_counter()`, **min of 5 cold runs** (before each run `GLNRep._cache`, the weight-multiplicity cache, `Recursion._cache`, `ExactRecursion._cache` are cleared), single process. **t_a** = `GLNRep` (patterns, weights, diagonal + adjacent exact generators; the 6 non-adjacent generators are lazy and NOT included, their cost is reported in §6) + Racah–Speiser (incl. the GT weight multiplicities) + exact PRV spaces of all irreps; no n_W-sized object is touched. **t_b** = t_a + metric, orthonormal float matrices, dense raising blocks and left inverses of all `A_a` (`t_setup`) + the float reconstruction of EVERY full HW vector (for m ≥ 2 all m PRV basis vectors, batched per irrep; `t_recon`). `entries` = Σ n_W·m = number of coefficients written. `worst HW resid` = untimed pass over all irreps (all E_ij for dim ≤ 400, the simple raising operators otherwise). **The machine was shared with agent G2** (a Python process of its own kept at least one core busy during my runs): the first block of timings is up to 10–15% slower than the quietest re-run; §6 quotes the quietest of three blocks and shows all three for the key S.

## 1. Headline results

1. **`glN_reps.py` works for any N, and its GT basis is Munich's.** Exact sparse GT irreps of gl(N) (`GLNRep.of(lam)`); N = 3 identical to `slice_engine.Rep` / `GL3Rep` (all nine E_ij, 4 irreps); N = 4: all 256 commutation relations exact on 9 irreps, N = 5: all 625 on 2; dimensions = Weyl, HW state last, weight multiset S_N-symmetric, metric exact (0 violations), C₂ scalar on every basis vector; Racah–Speiser Σ m dim V_μ = dim·dim′ on 26 products and Munich's SU(4) table (S = 1..8: irreps, dimensions, multiplicities) reproduced; **10/10 tests pass**. The pattern order is Munich's Q order, 6024 orthonormal lowering matrix elements equal Munich's `lowering_coeff²` exactly, and the sign gauge between the two bases is trivial (§9).
2. **dim M_μ = Racah–Speiser multiplicity for every dominant candidate μ — including those of multiplicity 0, where M_μ = {0} — of all 26 products**: the 16 products of G1.2 (λ ≠ λ′ in both orders), [S,S,0,0]² for S = 1..10, and the multiplicity families [4,2,0,0]², [6,4,2,0]², [9,6,3,0]² (max m = 3, 12, 27): **0 mismatches**. N = 3 regression: the PRV bases (W, basis vectors, free columns) are identical to `slice_engine.prv_space` on [4,2,0]², [6,3,0]², [14,7,0]² (151 candidates), and the float recursion reproduces `slice_engine.hw_from_slice` to 5e-15.
3. **The level recursion reconstructs every HW vector.** Float64: the HW residual of ALL E_ij (i < j) is ≤ 9.1e-15 on every irrep of [S,S,0,0]² for S ≤ 10 (incl. the singlet (8,8,8,8): 8.8e-15, the top (16,16,0,0): 0), ≤ 9.8e-15 on [6,4,2,0]² (m ≤ 12), ≤ 4.4e-14 on [9,6,3,0]² (1592 vectors, m ≤ 27), ≤ 3e-15 on all other products; the total weight is exactly μ and nW agrees with Racah–Speiser everywhere. Exact twin: **0 exact HW violations (E_ij^{tot}Ψ = 0 for all i < j) on EVERY irrep of EVERY product listed**, with all consistency (T) and dangling (D) conditions verified exactly (81 528 conditions for [10,10,0,0]², 176 256 for [6,4,2,0]²). Negative tests: **1655 of 1655 random slices outside M_μ rejected, 708 of 708 inside accepted**, and "some exact condition violated ⇔ the vector is not exactly HW" never fails.
4. **Munich-style brute force** (null space of the stacked total raising operators on the n_W-dimensional weight subspace, no PRV, no recursion): null dimension = multiplicity on every irrep of all products with dim·dim′ ≤ 10⁴, of [S,S,0,0]² for S ≤ 10 (S ≤ 8 is Munich's range; S = 10: n_W up to 4290) and of [4,2,0,0]², [6,4,2,0]² (71 irreps with n_W ≤ 2600, m up to 12): **0 mismatches**; our vectors lie in the null space to ≤ 1.6e-13 (largest principal angle) and span it.
5. **Benchmark against Munich's SU(4) table** (§6; different machine: Munich i5-1345U, ours i7-12700KF): all HW vectors of [8,8,0,0]² (dim 825, 45 irreps, n_W ≤ 1749): **t_a = 26 ms** (generators + Racah–Speiser + exact PRV spaces), **t_b = 93 ms** (+ metric, float matrices, left inverses, all 45 float vectors, 15 625 coefficients) against Munich's **80.19 s (×860)**; S = 10 (not runnable by Munich): 0.23 s; S = 12: 0.44 s; S = 20 (dim 19 481, 231 irreps, 3.06·10⁶ coefficients): 4.6 s. Munich grows ×5–7 per unit of S, ours ×1.4–2; per irrep the singlet (Munich's slowest, 18.7 s) takes 4.5 ms warm / 45 ms cold. (For S ≤ 2 Munich's C++ is faster in absolute terms: Python overhead.)
6. **Multiplicity families for G3** (§7): [3a,2a,a,0]² has max m = 4, 12, 27, 51 for a = 1..4 (a single irrep each time), [2a,a,0,0]² and [2a,a,a,0]² have max m = a + 1; the two Pieri-type families are multiplicity-free. Even the stress case [9,6,3,0]² (dim 4096, 213 irreps, 1592 vectors, m up to 27, w up to 32, n_W up to 59 008, 8.9·10⁶ coefficients) takes t_a = 0.23 s, t_b = 0.73 s.
7. **Ten reference files** `gl4-cg-slices-*.txt` (250 kB, exact rationals throughout, exact GNORM2, exact full vectors for the small products), re-verified from the TEXT with code independent of the engine: **0 mismatches at 1e-12** (exact checks exact, float ones ≤ 9.1e-15); a mutation test of the verifier detects 11 of 11 single-record corruptions while the pristine copy passes.
8. **Findings**: (a) for gl(4) the dangling PRV conditions are necessary in the exact twin — 253 of 1655 random out-of-M_μ slices (15%) pass every consistency (T) condition and are rejected only by the dangling (D) conditions (gl(3): none of ~500); (b) the smallest non-null eigenvalue of AᵀA (Munich's null-space problem) is not 1 as in gl(3): it is 2 for [1,0,0,0]² and 0.586, 0.438, 0.354, 0.298, … ≈ 1/S for [S,S,0,0]², S = 1, 2, 3, 4, … (0.156 at S = 10); (c) the float64 recursion's residual grows ×2.4–3.7 per unit of S beyond S = 14 (2.6e-13 at S = 16, 3.9e-11 at S = 20), `hw_from_slice_rounded` keeps 1e-15; (d) for λ ≠ λ′ the choice of the factor at site 1 changes w (e.g. 4 vs 2 for [3,2,1,0]⊗[2,1,0,0]).
9. **Not done**: the literal comparison with Munich's HW vectors (binary not buildable here); §9 gives the exact `hwcgc.sh` commands, the conversion (`CGC = ±Ψ·√(g_ν g_q/GNORM2)`) and a ready comparer (`munichcmp`), after showing that the two GT bases coincide.

## 2. G1.1 `glN_reps.py` and `test_glN_reps.py`

| test | result | what was checked |
|---|---|---|
| test_n3_regression | PASS | 4 irreps: all nine E_ij identical to slice_engine.Rep (336 off-diagonal entries) and to GL3Rep (adjacent+diagonal), metric identical, Racah-Speiser identical to slice_engine on 6 products incl. lam != lam' |
| test_commutation_n4 | PASS | N = 4: all 256 relations exact on 9 irreps (2304 relations x all basis vectors, 0 violations): [[1, 0, 0, 0], [1, 1, 0, 0], [2, 0, 0, 0], [2, 1, 0, 0], [2, 1, 1, 0], [2, 2, 0, 0], [3, 2, 1, 0], [1, 1, 1, 0], [3, 3, 0, 0]] |
| test_commutation_n5 | PASS | N = 5: all 625 relations exact on [[1, 0, 0, 0, 0], [2, 1, 0, 0, 0]] (1250 relations, 0 violations) |
| test_path_independence | PASS | E_14, E_41, E_13, E_31, E_24, E_42 (N = 4, three irreps) and E_15 (N = 5) are the same matrix along every commutator path |
| test_dims_weights_hw | PASS | 15 irreps (N = 3, 4, 5): dim = Weyl, HW pattern last with weight lam, weight multiset S_N-symmetric, sum_w m(w) = dim; dims e.g. {'[2, 1, 0, 0]': 20, '[3, 2, 1, 0]': 64, '[4, 4, 0, 0]': 105, '[10, 10, 0, 0]': 1716} |
| test_metric | PASS | 11 irreps (N = 3, 4, 5): E_ij^T G = G E_ji exact on every element of all off-diagonal pairs (0 violations), g[HW] = 1, g > 0, orthonormal Eh_ji = Eh_ij^T |
| test_casimir | PASS | C_2 = Sum_ij E_ij E_ji acts as the scalar Sum_i lam_i(lam_i + N + 1 - 2i) on every basis vector of 11 irreps: {'[1, 0, 0, 0]': 4, '[1, 1, 0, 0]': 6, '[2, 1, 0, 0]': 12, '[2, 1, 1, 0]': 12, '[2, 2, 0, 0]': 16, '[3, 2, 1, 0]': 24, '[3, 3, 0, 0]': 30, '[2, 1, 0, 0, 0]': 15, '[2, 1, 0]': 9, '[4, 2, 0]': 28, '[3, 2, 1]': 18} |
| test_racah_speiser_dimsums | PASS | sum_mu m dim V_mu = dim V_lam dim V_lam' on 26 products (G1.2 list in both orders + [S,S,0,0]^2, S = 1..10), RS symmetric in the two factors, all m > 0 |
| test_munich_table | PASS | Munich hw_timings.txt reproduced for S = 1..8: same irreps (our mu minus mu_4), same dimensions, all multiplicities 1; #irreps = [3, 6, 10, 15, 21, 28, 36, 45] = (S+1)(S+2)/2 |
| test_munich_gt_basis | PASS | 17 irreps (N = 3, 4): pattern order = Munich Q order, 6024 orthonormal lowering matrix elements equal Munich lowering_coeff^2 exactly, sign gauge trivial |

Notes: (i) the nine E_ij of `GLNRep` equal `slice_engine.Rep.E` entry for entry (the six off-diagonal ones; diagonal: `Rep.wt`) and `GL3Rep.generator(i,j)` (adjacent + diagonal) on [2,1,0], [3,1,0], [3,2,0], [4,2,0]; the exact metric, `pats`, `wt`, `wsp`, `top` are identical; (ii) the relations are checked on every basis vector with exact Fractions; E_14 = [E_13,E_34] = [E_12,E_24] (and the lower-triangular analogues, and E_15 along three paths for N = 5) are the same matrix, so the lazily built generators do not depend on the path; (iii)–(v) as listed; (vi) the Munich table is parsed from `hw_timings.txt` (SU(4) labels with μ₄ = 0): all irreps, dimensions and multiplicities of S = 1..8 are reproduced and all multiplicities are 1; (vii) is the Munich-basis test of §9. `weight_mults` counts GT patterns (O(dim V)): Racah–Speiser incl. the weight multiplicities costs 0.034 s at dim 19 481 (S = 20).

## 3. G1.2 `prv_space`: dim M_μ = Racah–Speiser multiplicity

`prv_space(lam, lamp, mu)` = exact nullspace over ℚ (RREF, Fractions; canonical basis: one vector per free column, 1 there, 0 at the other free columns, coordinates = free entries; `PRV` has `.W .basis .free .m .w .vec(k) .coords(v) .float_basis()`) of the N−1 stacked images E_{i,i+1}^{λ′_i−λ′_{i+1}+1} on V_λ[μ−λ′]. **Every dominant candidate μ of every product below, including those of multiplicity 0 (where M_μ = {0}), has dim M_μ = Racah–Speiser multiplicity: 0 mismatches in all 26 products.** `max w` = max over the irreps of dim V_λ[μ−λ′] (site 1 = first factor), `max n_W` = max over the irreps of the dimension of the weight-μ subspace of V_λ⊗V_λ′. Both orders are listed for every λ ≠ λ′. (The last three rows, [4,2,0,0]², [6,4,2,0]², [9,6,3,0]², are the multiplicity families of §7 with max m = 3, 12, 27.)

| product | dim x dim' | #candidates | #irreps | sum m | max m | max w | max n_W | sum m dim V = dim dim' | mismatches |
|---|---|---|---|---|---|---|---|---|---|
| [1,0,0,0]^2 | 4 x 4 | 2 | 2 | 2 | 1 | 1 | 2 | True | 0 |
| [1,1,0,0]^2 | 6 x 6 | 3 | 3 | 3 | 1 | 1 | 6 | True | 0 |
| [2,1,0,0]^2 | 20 x 20 | 7 | 7 | 8 | 2 | 2 | 20 | True | 0 |
| [2,1,1,0]^2 | 15 x 15 | 6 | 6 | 7 | 2 | 3 | 21 | True | 0 |
| [2,2,0,0]^2 | 20 x 20 | 8 | 6 | 6 | 1 | 2 | 22 | True | 0 |
| [3,2,1,0]^2 | 64 x 64 | 16 | 16 | 30 | 4 | 4 | 152 | True | 0 |
| [2,2,1,0]^2 | 20 x 20 | 7 | 7 | 8 | 2 | 2 | 20 | True | 0 |
| [3,3,0,0]^2 | 50 x 50 | 16 | 10 | 10 | 1 | 2 | 62 | True | 0 |
| [2,1,0,0]x[1,1,0,0] | 20 x 6 | 4 | 4 | 4 | 1 | 2 | 9 | True | 0 |
| [1,1,0,0]x[2,1,0,0] | 6 x 20 | 4 | 4 | 4 | 1 | 1 | 9 | True | 0 |
| [3,1,0,0]x[2,2,0,0] | 45 x 20 | 11 | 7 | 7 | 1 | 2 | 24 | True | 0 |
| [2,2,0,0]x[3,1,0,0] | 20 x 45 | 11 | 7 | 7 | 1 | 2 | 24 | True | 0 |
| [2,2,1,0]x[2,1,1,0] | 20 x 15 | 6 | 6 | 7 | 2 | 2 | 18 | True | 0 |
| [2,1,1,0]x[2,2,1,0] | 15 x 20 | 6 | 6 | 7 | 2 | 3 | 18 | True | 0 |
| [3,2,1,0]x[2,1,0,0] | 64 x 20 | 10 | 10 | 14 | 2 | 4 | 52 | True | 0 |
| [2,1,0,0]x[3,2,1,0] | 20 x 64 | 10 | 10 | 14 | 2 | 2 | 52 | True | 0 |
| [4,4,0,0]^2 | 105 x 105 | 29 | 15 | 15 | 1 | 3 | 147 | True | 0 |
| [5,5,0,0]^2 | 196 x 196 | 47 | 21 | 21 | 1 | 3 | 308 | True | 0 |
| [6,6,0,0]^2 | 336 x 336 | 72 | 28 | 28 | 1 | 4 | 588 | True | 0 |
| [7,7,0,0]^2 | 540 x 540 | 104 | 36 | 36 | 1 | 4 | 1044 | True | 0 |
| [8,8,0,0]^2 | 825 x 825 | 145 | 45 | 45 | 1 | 5 | 1749 | True | 0 |
| [9,9,0,0]^2 | 1210 x 1210 | 195 | 55 | 55 | 1 | 5 | 2794 | True | 0 |
| [10,10,0,0]^2 | 1716 x 1716 | 256 | 66 | 66 | 1 | 6 | 4290 | True | 0 |
| [4,2,0,0]^2 | 126 x 126 | 26 | 24 | 32 | 3 | 4 | 174 | True | 0 |
| [6,4,2,0]^2 | 729 x 729 | 76 | 76 | 297 | 12 | 15 | 4743 | True | 0 |
| [9,6,3,0]^2 | 4096 x 4096 | 213 | 213 | 1592 | 27 | 32 | 59008 | True | 0 |

**N = 3 regression** against `slice_engine.prv_space` / `hw_from_slice` / `racah_speiser`:

* [4,2,0]^2      N = 3 regression: 13 dominant candidates, PRV (W, basis, free) different from slice_engine.prv_space: 0; float recursion vs slice_engine.hw_from_slice on 19 vectors: max relative deviation 2.66e-16
* [6,3,0]^2      N = 3 regression: 25 dominant candidates, PRV (W, basis, free) different from slice_engine.prv_space: 0; float recursion vs slice_engine.hw_from_slice on 44 vectors: max relative deviation 3.81e-16
* [14,7,0]^2     N = 3 regression: 113 dominant candidates, PRV (W, basis, free) different from slice_engine.prv_space: 0; float recursion vs slice_engine.hw_from_slice on 344 vectors: max relative deviation 5.19e-15
* [3,2,0]x[3,1,0] N = 3 lam != lam': 8 irreps (Racah-Speiser identical to slice_engine), dim M_mu mismatches in both orders: 0

(The PRV bases coincide with the gl(3) engine exactly — same W, same basis vectors, same free columns — for all 13 + 25 + 113 dominant candidates; the float recursion reproduces `slice_engine.hw_from_slice` entry for entry to 5e-15 on all 19 + 44 + 344 vectors.)

## 4. G1.3 The level recursion

### 4.1 Float recursion (all m PRV basis vectors of every irrep; HW residual of ALL E_ij with i < j; total weight)

| product | irreps | vectors | max HW residual (all E_ij, i<j) | at mu | weight dev. | n_W mismatches | sectors > 1e-12 | min sv ratio of A_a |
|---|---|---|---|---|---|---|---|---|
| [1,0,0,0]^2 | 2 | 2 | 0.00e+00 | (2, 0, 0, 0) | 0 | 0 | 0 | 1.000 |
| [1,1,0,0]^2 | 3 | 3 | 1.92e-16 | (1, 1, 1, 1) | 0 | 0 | 0 | 1.000 |
| [2,1,0,0]^2 | 7 | 8 | 3.47e-16 | (2, 2, 2, 0) | 0 | 0 | 0 | 0.577 |
| [2,1,1,0]^2 | 6 | 7 | 2.94e-15 | (3, 2, 2, 1) | 0 | 0 | 0 | 0.414 |
| [2,2,0,0]^2 | 6 | 6 | 2.91e-16 | (2, 2, 2, 2) | 0 | 0 | 0 | 0.518 |
| [3,2,1,0]^2 | 16 | 30 | 2.49e-15 | (4, 3, 3, 2) | 0 | 0 | 0 | 0.342 |
| [2,2,1,0]^2 | 7 | 8 | 3.47e-16 | (4, 2, 2, 2) | 0 | 0 | 0 | 0.577 |
| [3,3,0,0]^2 | 10 | 10 | 5.89e-16 | (3, 3, 3, 3) | 0 | 0 | 0 | 0.509 |
| [2,1,0,0]x[1,1,0,0] | 4 | 4 | 1.48e-16 | (2, 1, 1, 1) | 0 | 0 | 0 | 1.000 |
| [1,1,0,0]x[2,1,0,0] | 4 | 4 | 2.05e-16 | (2, 1, 1, 1) | 0 | 0 | 0 | 0.577 |
| [3,1,0,0]x[2,2,0,0] | 7 | 7 | 3.77e-16 | (3, 2, 2, 1) | 0 | 0 | 0 | 0.518 |
| [2,2,0,0]x[3,1,0,0] | 7 | 7 | 5.48e-16 | (4, 2, 2, 0) | 0 | 0 | 0 | 0.500 |
| [2,2,1,0]x[2,1,1,0] | 6 | 7 | 2.60e-15 | (3, 3, 2, 1) | 0 | 0 | 0 | 0.414 |
| [2,1,1,0]x[2,2,1,0] | 6 | 7 | 4.21e-16 | (3, 2, 2, 2) | 0 | 0 | 0 | 0.577 |
| [3,2,1,0]x[2,1,0,0] | 10 | 14 | 5.68e-16 | (3, 3, 2, 1) | 0 | 0 | 0 | 0.577 |
| [2,1,0,0]x[3,2,1,0] | 10 | 14 | 1.71e-15 | (3, 2, 2, 2) | 0 | 0 | 0 | 0.342 |
| [4,4,0,0]^2 | 15 | 15 | 8.24e-16 | (5, 4, 4, 3) | 0 | 0 | 0 | 0.451 |
| [5,5,0,0]^2 | 21 | 21 | 1.80e-15 | (8, 6, 4, 2) | 0 | 0 | 0 | 0.436 |
| [6,6,0,0]^2 | 28 | 28 | 1.77e-15 | (9, 7, 5, 3) | 0 | 0 | 0 | 0.400 |
| [7,7,0,0]^2 | 36 | 36 | 2.61e-15 | (13, 7, 7, 1) | 0 | 0 | 0 | 0.385 |
| [8,8,0,0]^2 | 45 | 45 | 9.14e-15 | (9, 8, 8, 7) | 0 | 0 | 0 | 0.362 |
| [9,9,0,0]^2 | 55 | 55 | 7.73e-15 | (11, 9, 9, 7) | 0 | 0 | 0 | 0.349 |
| [10,10,0,0]^2 | 66 | 66 | 4.72e-15 | (11, 10, 10, 9) | 0 | 0 | 0 | 0.333 |
| [4,2,0,0]^2 | 24 | 32 | 1.64e-15 | (5, 3, 2, 2) | 0 | 0 | 0 | 0.434 |
| [6,4,2,0]^2 | 76 | 297 | 9.79e-15 | (7, 6, 6, 5) | 0 | 0 | 0 | 0.234 |
| [9,6,3,0]^2 | 213 | 1592 | 4.35e-14 | (11, 9, 9, 7) | 0 | 0 | 0 | 0.179 |

Weights are exactly μ and `hw.nW` equals the Racah–Speiser n_W for every sector (0 mismatches). The recursion is well conditioned: the smallest singular-value ratio of any stacked site-2 block A_a (`min sv ratio`) is 0.33 at S = 10, 0.25 at S = 20, 0.23 for [6,4,2,0]² and 0.18 for [9,6,3,0]² (1.0 where the A_a are single rows). **All residuals are ≤ 1e-12 for every product of this table (worst 4.4e-14, [9,6,3,0]²) and up to S = 14 (4.6e-14); beyond S ≈ 16 the float64 recursion degrades (§4.4) and the exact-rounded path is the cure.**

**[8,8,0,0]², all 45 irreps** (the whole range of Munich's table; Munich's SU(4) labels; the singlet is Munich's slowest row, the top irrep (16,16,0,0) its fastest; weights exactly μ, residuals ≤ 9.1e-15 including the singlet (8,8,8,8) and the top (16,16,0,0)):

| mu (ours) | mu (Munich) | dim V_mu | n_W | w | HW residual (all E_ij) | weight dev. |
|---|---|---|---|---|---|---|
| (16,16,0,0) | (16 16 0 0) | 8721 | 1 | 1 | 0.00e+00 | 0 |
| (16,15,1,0) | (16 15 1 0) | 27455 | 2 | 1 | 0.00e+00 | 0 |
| (16,14,2,0) | (16 14 2 0) | 47424 | 3 | 1 | 0.00e+00 | 0 |
| (16,13,3,0) | (16 13 3 0) | 62700 | 4 | 1 | 2.48e-16 | 0 |
| (16,12,4,0) | (16 12 4 0) | 69825 | 5 | 1 | 1.77e-16 | 0 |
| (16,11,5,0) | (16 11 5 0) | 67431 | 6 | 1 | 3.01e-16 | 0 |
| (16,10,6,0) | (16 10 6 0) | 55860 | 7 | 1 | 2.63e-16 | 0 |
| (16,9,7,0) | (16 9 7 0) | 36784 | 8 | 1 | 0.00e+00 | 0 |
| (16,8,8,0) | (16 8 8 0) | 12825 | 9 | 1 | 0.00e+00 | 0 |
| (15,15,1,1) | (14 14 0 0) | 5440 | 8 | 2 | 5.71e-16 | 0 |
| (15,14,2,1) | (14 13 1 0) | 16575 | 14 | 2 | 5.74e-16 | 0 |
| (15,13,3,1) | (14 12 2 0) | 27489 | 20 | 2 | 6.09e-16 | 0 |
| (15,12,4,1) | (14 11 3 0) | 34476 | 26 | 2 | 6.43e-16 | 0 |
| (15,11,5,1) | (14 10 4 0) | 35700 | 32 | 2 | 7.58e-16 | 0 |
| (15,10,6,1) | (14 9 5 0) | 30855 | 38 | 2 | 2.13e-15 | 0 |
| (15,9,7,1) | (14 8 6 0) | 20825 | 44 | 2 | 1.15e-15 | 0 |
| (15,8,8,1) | (14 7 7 0) | 7344 | 48 | 1 | 1.26e-15 | 0 |
| (14,14,2,2) | (12 12 0 0) | 3185 | 34 | 3 | 7.72e-16 | 0 |
| (14,13,3,2) | (12 11 1 0) | 9295 | 54 | 3 | 7.96e-16 | 0 |
| (14,12,4,2) | (12 10 2 0) | 14580 | 74 | 3 | 9.17e-16 | 0 |
| (14,11,5,2) | (12 9 3 0) | 16940 | 94 | 3 | 1.57e-15 | 0 |
| (14,10,6,2) | (12 8 4 0) | 15625 | 114 | 3 | 1.81e-15 | 0 |
| (14,9,7,2) | (12 7 5 0) | 10935 | 132 | 2 | 1.75e-15 | 0 |
| (14,8,8,2) | (12 6 6 0) | 3920 | 142 | 1 | 3.18e-15 | 0 |
| (13,13,3,3) | (10 10 0 0) | 1716 | 104 | 4 | 9.20e-16 | 0 |
| (13,12,4,3) | (10 9 1 0) | 4719 | 154 | 4 | 1.59e-15 | 0 |
| (13,11,5,3) | (10 8 2 0) | 6825 | 204 | 4 | 2.01e-15 | 0 |
| (13,10,6,3) | (10 7 3 0) | 7020 | 252 | 3 | 3.37e-15 | 0 |
| (13,9,7,3) | (10 6 4 0) | 5200 | 292 | 2 | 1.92e-15 | 0 |
| (13,8,8,3) | (10 5 5 0) | 1911 | 312 | 1 | 3.26e-15 | 0 |
| (12,12,4,4) | (8 8 0 0) | 825 | 259 | 5 | 1.80e-15 | 0 |
| (12,11,5,4) | (8 7 1 0) | 2079 | 362 | 4 | 2.02e-15 | 0 |
| (12,10,6,4) | (8 6 2 0) | 2640 | 457 | 3 | 2.36e-15 | 0 |
| (12,9,7,4) | (8 5 3 0) | 2156 | 532 | 2 | 2.52e-15 | 0 |
| (12,8,8,4) | (8 4 4 0) | 825 | 567 | 1 | 2.79e-15 | 0 |
| (11,11,5,5) | (6 6 0 0) | 336 | 548 | 4 | 2.12e-15 | 0 |
| (11,10,6,5) | (6 5 1 0) | 735 | 714 | 3 | 2.64e-15 | 0 |
| (11,9,7,5) | (6 4 2 0) | 729 | 840 | 2 | 3.30e-15 | 0 |
| (11,8,8,5) | (6 3 3 0) | 300 | 896 | 1 | 4.26e-15 | 0 |
| (10,10,6,6) | (4 4 0 0) | 105 | 980 | 3 | 2.71e-15 | 0 |
| (10,9,7,6) | (4 3 1 0) | 175 | 1176 | 2 | 6.14e-15 | 0 |
| (10,8,8,6) | (4 2 2 0) | 84 | 1260 | 1 | 6.69e-15 | 0 |
| (9,9,7,7) | (2 2 0 0) | 20 | 1464 | 2 | 8.32e-15 | 0 |
| (9,8,8,7) | (2 1 1 0) | 15 | 1584 | 1 | 9.14e-15 | 0 |
| (8,8,8,8) | (0 0 0 0) | 1 | 1749 | 1 | 8.82e-15 | 0 |

### 4.2 Exact twin (Fractions; T- and D-conditions verified; exact check of E_ij^{tot}Ψ = 0 for ALL i < j)

The brief asked for the exact twin on the products with dim·dim′ ≤ 2000 and on ≥ 5 irreps of [3,2,1,0]² and [4,4,0,0]²; it is cheap enough to run on **every irrep and every PRV basis vector of every product listed, including [S,S,0,0]² up to S = 10 (2 s) and [6,4,2,0]² (297 vectors, 185 606 coefficients, 8 s)**; [9,6,3,0]² is only run through the float recursion (one exact basis vector of its m = 27 sector, 4158 entries, takes 1.2 s).

| product | irreps | PRV basis vectors | entries | exact conditions verified (T + dangling) | exact HW violations (all E_ij, i<j) |
|---|---|---|---|---|---|
| [1,0,0,0]^2 | 2 of 2 | 2 | 3 | 0 | 0 |
| [1,1,0,0]^2 | 3 of 3 | 3 | 9 | 1 | 0 |
| [2,1,0,0]^2 | 7 of 7 | 8 | 59 | 17 | 0 |
| [2,1,1,0]^2 | 6 of 6 | 7 | 50 | 13 | 0 |
| [2,2,0,0]^2 | 6 of 6 | 6 | 48 | 24 | 0 |
| [3,2,1,0]^2 | 16 of 16 | 30 | 1069 | 542 | 0 |
| [2,2,1,0]^2 | 7 of 7 | 8 | 59 | 17 | 0 |
| [3,3,0,0]^2 | 10 of 10 | 10 | 186 | 149 | 0 |
| [2,1,0,0]x[1,1,0,0] | 4 of 4 | 4 | 16 | 4 | 0 |
| [1,1,0,0]x[2,1,0,0] | 4 of 4 | 4 | 16 | 2 | 0 |
| [3,1,0,0]x[2,2,0,0] | 7 of 7 | 7 | 57 | 22 | 0 |
| [2,2,0,0]x[3,1,0,0] | 7 of 7 | 7 | 57 | 19 | 0 |
| [2,2,1,0]x[2,1,1,0] | 6 of 6 | 7 | 48 | 11 | 0 |
| [2,1,1,0]x[2,2,1,0] | 6 of 6 | 7 | 48 | 11 | 0 |
| [3,2,1,0]x[2,1,0,0] | 10 of 10 | 14 | 211 | 88 | 0 |
| [2,1,0,0]x[3,2,1,0] | 10 of 10 | 14 | 211 | 72 | 0 |
| [4,4,0,0]^2 | 15 of 15 | 15 | 580 | 592 | 0 |
| [5,5,0,0]^2 | 21 of 21 | 21 | 1544 | 1822 | 0 |
| [6,6,0,0]^2 | 28 of 28 | 28 | 3646 | 4736 | 0 |
| [7,7,0,0]^2 | 36 of 36 | 36 | 7838 | 10894 | 0 |
| [8,8,0,0]^2 | 45 of 45 | 45 | 15625 | 22832 | 0 |
| [9,9,0,0]^2 | 55 of 55 | 55 | 29277 | 44455 | 0 |
| [10,10,0,0]^2 | 66 of 66 | 66 | 52090 | 81528 | 0 |
| [4,2,0,0]^2 | 24 of 24 | 32 | 1287 | 845 | 0 |
| [6,4,2,0]^2 | 76 of 76 | 297 | 185606 | 176256 | 0 |

### 4.3 Negative tests

Random exact vectors of V_λ[μ−λ′] outside M_μ (m < w sectors; 20 per sector, 5 for [6,4,2,0]²; entries p/q with |p| ≤ 9, q ≤ 5) must be rejected by the exact conditions; random exact vectors inside M_μ must be accepted; and the verdict of the conditions must agree with the exact HW check of the reconstructed vector (cross-validation of the condition logic: *some condition violated ⇔ the vector is not exactly HW*). The last column is the float residual of the float recursion run on the out-of-M_μ slices (the float recursion itself cannot see an inconsistent slice; the residual check can).

| product | sectors w > m | slices outside M_mu | rejected | by a T-condition | ONLY by a dangling condition | slices inside M_mu | accepted | cond. vs exact-HW inconsistencies | min float HW residual of a rejected slice |
|---|---|---|---|---|---|---|---|---|---|
| [2,1,0,0]^2 | 2 | 40 | 40 | 20 | 20 | 21 | 21 | 0 | 6.3e-02 |
| [2,1,1,0]^2 | 1 | 20 | 20 | 0 | 20 | 18 | 18 | 0 | 3.7e-03 |
| [2,2,0,0]^2 | 1 | 18 | 18 | 0 | 18 | 18 | 18 | 0 | 5.8e-02 |
| [3,2,1,0]^2 | 5 | 96 | 96 | 96 | 0 | 48 | 48 | 0 | 2.9e-02 |
| [2,1,0,0]x[1,1,0,0] | 2 | 40 | 40 | 0 | 40 | 12 | 12 | 0 | 6.0e-02 |
| [3,1,0,0]x[2,2,0,0] | 4 | 79 | 79 | 60 | 19 | 21 | 21 | 0 | 1.4e-02 |
| [2,2,1,0]x[2,1,1,0] | 2 | 38 | 38 | 19 | 19 | 18 | 18 | 0 | 7.5e-02 |
| [3,3,0,0]^2 | 3 | 59 | 59 | 40 | 19 | 30 | 30 | 0 | 3.4e-02 |
| [3,2,1,0]x[2,1,0,0] | 5 | 100 | 100 | 60 | 40 | 30 | 30 | 0 | 7.8e-02 |
| [4,4,0,0]^2 | 6 | 115 | 115 | 96 | 19 | 45 | 45 | 0 | 3.6e-02 |
| [6,6,0,0]^2 | 15 | 287 | 287 | 268 | 19 | 84 | 84 | 0 | 2.2e-02 |
| [8,8,0,0]^2 | 28 | 546 | 546 | 526 | 20 | 135 | 135 | 0 | 4.5e-02 |
| [6,4,2,0]^2 | 44 | 217 | 217 | 217 | 0 | 228 | 228 | 0 | 2.4e-02 |

**All 1655 random slices outside M_μ are rejected and all 708 inside are accepted; the condition logic and the exact HW check never disagree.** **Finding**: 253 of the 1655 out-of-M_μ slices (15%) pass ALL T-conditions and are rejected only by the dangling D-conditions: in every product with a sector of w − m = 1 (e.g. [2,1,1,0]²: 20 of 20, [2,2,0,0]²: 18 of 18, [2,1,0,0]⊗[1,1,0,0]: 40 of 40, [8,8,0,0]²: 20 of 546). For gl(3) the T-conditions alone rejected every one of ~500 random out-of-M_μ slices that I tried on [3,1,0], [3,2,0], [4,2,0], [6,3,0], [4,1,0], [5,2,0], [6,2,0] (0 accepted), so the gl(3) exact twin is adequate there, but **for gl(4) (and for N = 2, where there are no T-rows at all) the D-check is necessary**. The float residual of an out-of-M_μ slice is ≥ 3.7e-3 against ≤ 1e-13 inside: the residual check separates the two cases by ten orders of magnitude.

### 4.4 The float recursion's residual vs the exact-rounded path

The float64 recursion's HW residual grows with the size of the representation (the rounding of the double slice is amplified by the chain of left inverses): 5.6e-15 at S = 8, 1.3e-14 at S = 12, 4.6e-14 at S = 14, 2.6e-13 at S = 16, 3.9e-11 at S = 20, 7.5e-9 at S = 24 (×2.4 per unit of S between 14 and 16, ×3.5–3.7 per unit between 16 and 24); **the threshold 1e-12 holds up to S = 16 and fails at S = 20**. `hw_from_slice_rounded` (exact recursion rounded to double) keeps ~1e-15 at every size (the loss is in the float arithmetic of the recursion, not in the rounding of the slice):

| S | mu | n_W | HW residual float64 recursion | time [s] | HW residual exact recursion rounded to double | time [s] | max rel. difference of the vectors |
|---|---|---|---|---|---|---|---|
| 16 | (16, 16, 16, 16) | 31977 | 2.6e-13 | 0.040 | 9.5e-16 | 1.93 | 1.1e-12 |
| 16 | (17, 16, 16, 15) | 31008 | 1.7e-13 | 0.041 | 1.0e-15 | 1.24 | 4.1e-13 |
| 16 | (18, 16, 16, 14) | 28696 | 1.1e-13 | 0.032 | 1.1e-15 | 1.10 | 1.7e-13 |
| 20 | (20, 20, 20, 20) | 86779 | 3.9e-11 | 0.094 | 1.2e-15 | 6.16 | 2.4e-10 |
| 20 | (21, 20, 20, 19) | 85008 | 3.5e-11 | 0.071 | 1.3e-15 | 4.13 | 1.2e-10 |
| 20 | (22, 20, 20, 18) | 80598 | 2.6e-11 | 0.070 | 1.4e-15 | 3.96 | 5.8e-11 |
| 24 | (25, 24, 24, 23) | 196560 | 7.5e-09 | 0.140 | 1.5e-15 | 15.58 | 3.1e-08 |
| 24 | (24, 24, 24, 24) | 199485 | 7.1e-09 | 0.150 | 1.4e-15 | 10.91 | 4.8e-08 |
| 24 | (26, 24, 24, 22) | 189060 | 6.1e-09 | 0.125 | 1.6e-15 | 10.73 | 1.5e-08 |

(Cost of the exact path: 1.1–15.6 s per vector at n_W = 2.9·10⁴ … 2·10⁵ in Python Fractions; a batched/vectorised exact path was not written. At every size of the tables of this package except the extension row S = 20 (3.9e-11) the float path is within 1e-12; the S = 14 and S = 16 rows of §6 (4.6e-14, 2.6e-13) are inside.)

## 5. G1.4 Munich-style brute-force null space (no PRV theorem, no recursion)

For every irrep μ: the stacked total simple raising operators E_{i,i+1}^{tot} = E⊗1 + 1⊗E (scipy `kron`, orthonormal basis) restricted to the n_W columns of the weight-μ subspace, Gram matrix AᵀA, dense `eigh`; null dimension = multiplicity; our reconstructed vectors (orthonormalised) must lie in the null space and span it (largest principal angle = ‖(1 − NNᵀ)Q‖₂). **All products of G1.2 with dim·dim′ ≤ 10⁴, [S,S,0,0]² for S ≤ 5 (as asked) and, because it is cheap, S = 6..10 (S ≤ 8 is Munich's whole range; S = 10: 66 irreps, n_W up to 4290), plus [4,2,0,0]² and the 71 irreps of [6,4,2,0]² with n_W ≤ 2600 (multiplicities up to m = 12): 0 mismatches.** (The brute-force timings depend strongly on the machine load: dense multithreaded `eigh`; S = 10 took 68 s in one run and 125 s in another; they are not a benchmark.)

| product | irreps | max n_W | dim null != m | smallest non-null eigenvalue of A^T A | max distance of our vectors from the null space | largest principal angle | time [s] (load-dependent) |
|---|---|---|---|---|---|---|---|
| [1,0,0,0]^2 | 2 | 2 | 0 | 2.000e+00 | 2.2e-16 | 2.2e-16 | 0.0 |
| [1,1,0,0]^2 | 3 | 6 | 0 | 5.858e-01 | 2.3e-16 | 2.3e-16 | 0.0 |
| [2,1,0,0]^2 | 7 | 20 | 0 | 4.859e-01 | 1.8e-15 | 1.8e-15 | 0.0 |
| [2,1,1,0]^2 | 6 | 21 | 0 | 4.384e-01 | 2.9e-15 | 2.9e-15 | 0.0 |
| [2,2,0,0]^2 | 6 | 22 | 0 | 4.384e-01 | 1.4e-15 | 1.4e-15 | 0.0 |
| [3,2,1,0]^2 | 16 | 152 | 0 | 3.542e-01 | 1.7e-14 | 1.7e-14 | 0.0 |
| [2,2,1,0]^2 | 7 | 20 | 0 | 4.859e-01 | 2.1e-15 | 2.1e-15 | 0.0 |
| [3,3,0,0]^2 | 10 | 62 | 0 | 3.542e-01 | 4.6e-15 | 4.6e-15 | 0.0 |
| [2,1,0,0]x[1,1,0,0] | 4 | 9 | 0 | 5.188e-01 | 9.8e-16 | 9.8e-16 | 0.0 |
| [1,1,0,0]x[2,1,0,0] | 4 | 9 | 0 | 5.188e-01 | 1.3e-15 | 1.3e-15 | 0.0 |
| [3,1,0,0]x[2,2,0,0] | 7 | 24 | 0 | 4.384e-01 | 1.3e-15 | 1.3e-15 | 0.0 |
| [2,2,0,0]x[3,1,0,0] | 7 | 24 | 0 | 4.384e-01 | 1.8e-15 | 1.8e-15 | 0.0 |
| [2,2,1,0]x[2,1,1,0] | 6 | 18 | 0 | 5.188e-01 | 2.8e-15 | 2.8e-15 | 0.0 |
| [2,1,1,0]x[2,2,1,0] | 6 | 18 | 0 | 5.188e-01 | 1.1e-15 | 1.2e-15 | 0.0 |
| [3,2,1,0]x[2,1,0,0] | 10 | 52 | 0 | 3.983e-01 | 3.1e-15 | 3.1e-15 | 0.0 |
| [2,1,0,0]x[3,2,1,0] | 10 | 52 | 0 | 3.983e-01 | 2.7e-15 | 2.9e-15 | 0.0 |
| [1,1,0,0]^2 | 3 | 6 | 0 | 5.858e-01 | 2.3e-16 | 2.3e-16 | 0.0 |
| [2,2,0,0]^2 | 6 | 22 | 0 | 4.384e-01 | 1.4e-15 | 1.4e-15 | 0.0 |
| [3,3,0,0]^2 | 10 | 62 | 0 | 3.542e-01 | 4.6e-15 | 4.6e-15 | 0.0 |
| [4,4,0,0]^2 | 15 | 147 | 0 | 2.984e-01 | 8.0e-15 | 8.0e-15 | 0.0 |
| [5,5,0,0]^2 | 21 | 308 | 0 | 2.583e-01 | 1.3e-14 | 1.3e-14 | 0.1 |
| [6,6,0,0]^2 | 28 | 588 | 0 | 2.280e-01 | 3.7e-14 | 3.7e-14 | 0.6 |
| [7,7,0,0]^2 | 36 | 1044 | 0 | 2.042e-01 | 6.0e-14 | 6.0e-14 | 3.7 |
| [8,8,0,0]^2 | 45 | 1749 | 0 | 1.849e-01 | 8.7e-14 | 8.7e-14 | 14.6 |
| [9,9,0,0]^2 | 55 | 2794 | 0 | 1.690e-01 | 4.2e-14 | 4.2e-14 | 39.1 |
| [10,10,0,0]^2 | 66 | 4290 | 0 | 1.557e-01 | 7.5e-14 | 7.5e-14 | 124.6 |
| [4,2,0,0]^2 | 24 | 174 | 0 | 3.451e-01 | 2.2e-14 | 2.2e-14 | 0.0 |
| [6,4,2,0]^2 | 71 | 2520 | 0 | 2.280e-01 | 9.6e-14 | 1.6e-13 | 29.3 |

The smallest non-null eigenvalue is **not 1** as in gl(3) (there it was exactly 1 in every irrep: the sl(2) minimal singular value): it is 2 for [1,0,0,0]², 2 − √2 = 0.586 for [1,1,0,0]² and then decreases like ≈ 1/S (0.354, 0.298, 0.258, 0.228, 0.204, 0.185, 0.169, 0.156 at S = 3, 4, …, 10), so the gap that any null-space solver (Munich's included) has to resolve shrinks, but stays far above double precision: the Munich slowness is not a conditioning problem, it is the O(n_W³) dense linear algebra.

## 6. G1.5 Benchmark: [S,S,0,0]² against Munich's SU(4) table

Munich's column is `Paul/Cpp/AlexCGC/results/hw_timings.txt`: total `compute_highest_weight_coeffs()` time over all irreps (reps = 1), **13th-gen i5-1345U, g++ 11.4 + LAPACK 3.10 — a DIFFERENT machine from ours (i7-12700KF)**; the ratio is therefore only indicative (single-thread speed differences between the two CPUs are of order 1–2×, which does not change the orders of magnitude). Munich's code was not run for S ≥ 9 (and cannot run S = 10 on 16 GB according to the brief). Both sides produce the same object: the HW vector (all n_W coefficients) of every irrep of [S,S,0,0]²; ours is the float64 plain recursion result (HW residual column), the family is multiplicity-free so no missing-label operator is involved. **Noise**: three independent blocks of the whole table are in the output file (the first one inside `all`, two re-runs of `bench 12 ext` and `families`); the quietest block (smallest Σ t_b over S = 4..12) is quoted. Run-to-run differences are ≤ 3% between the two re-runs and ≤ 15% for the first block (more load from G2):

| S | block 0: t_a / t_b [s] | block 1: t_a / t_b [s] | block 2: t_a / t_b [s] |
|---|---|---|---|
| 4 | 0.0044 / 0.0105 | 0.0034 / 0.0088 | 0.0034 / 0.0090 |
| 6 | 0.0113 / 0.0341 | 0.0110 / 0.0326 | 0.0107 / 0.0320 |
| 8 | 0.0291 / 0.1012 | 0.0258 / 0.0932 | 0.0257 / 0.0918 |
| 10 | 0.0748 / 0.2579 | 0.0579 / 0.2276 | 0.0578 / 0.2213 |
| 12 | 0.1217 / 0.4964 | 0.1027 / 0.4427 | 0.1061 / 0.4653 |
| 14 | 0.2579 / 1.0165 | 0.1780 / 0.8831 | 0.2233 / 0.9088 |
| 16 | 0.3780 / 1.7573 | 0.3218 / 1.6202 | 0.3325 / 1.6544 |
| 20 | 0.9270 / 4.8989 | 0.8242 / 4.5883 | 0.7432 / 4.6196 |

| S | dim V | irreps (= copies) | max w | max n_W | t_gen | t_rs | t_prv | t_a | t_setup | t_recon | t_b | entries written | worst HW resid | Munich total hw-solve (x = Munich / t_b) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 6 | 3 | 1 | 6 | 0.0002 | 0.0000 | 0.0000 | 0.0002 | 0.0004 | 0.0001 | 0.0007 | 9 | 1.9e-16 | 2.785e-05 s (x0.04) |
| 2 | 20 | 6 | 2 | 22 | 0.0005 | 0.0001 | 0.0001 | 0.0007 | 0.0007 | 0.0004 | 0.0018 | 48 | 2.9e-16 | 0.0008577 s (x0.48) |
| 3 | 50 | 10 | 2 | 62 | 0.0013 | 0.0002 | 0.0002 | 0.0016 | 0.0012 | 0.0012 | 0.0040 | 186 | 5.9e-16 | 0.01072 s (x2.7) |
| 4 | 105 | 15 | 3 | 147 | 0.0027 | 0.0003 | 0.0003 | 0.0034 | 0.0023 | 0.0031 | 0.0088 | 580 | 8.2e-16 | 0.0574 s (x6.5) |
| 5 | 196 | 21 | 3 | 308 | 0.0054 | 0.0006 | 0.0004 | 0.0064 | 0.0041 | 0.0074 | 0.0179 | 1544 | 1.8e-15 | 0.3551 s (x20) |
| 6 | 336 | 28 | 4 | 588 | 0.0094 | 0.0010 | 0.0007 | 0.0110 | 0.0065 | 0.0150 | 0.0326 | 3646 | 1.8e-15 | 2.486 s (x76) |
| 7 | 540 | 36 | 4 | 1044 | 0.0150 | 0.0013 | 0.0009 | 0.0171 | 0.0105 | 0.0296 | 0.0572 | 7838 | 1.3e-15 | 15.31 s (x268) |
| 8 | 825 | 45 | 5 | 1749 | 0.0227 | 0.0020 | 0.0012 | 0.0258 | 0.0160 | 0.0514 | 0.0932 | 15625 | 5.6e-15 | 80.19 s (x860) |
| 9 | 1210 | 55 | 5 | 2794 | 0.0340 | 0.0026 | 0.0015 | 0.0381 | 0.0198 | 0.0838 | 0.1417 | 29277 | 4.6e-15 | n/a (not in hw_timings.txt) |
| 10 | 1716 | 66 | 6 | 4290 | 0.0520 | 0.0037 | 0.0022 | 0.0579 | 0.0301 | 0.1396 | 0.2276 | 52090 | 2.8e-15 | n/a (not in hw_timings.txt) |
| 11 | 2366 | 78 | 6 | 6370 | 0.0762 | 0.0059 | 0.0027 | 0.0847 | 0.0377 | 0.1988 | 0.3212 | 88700 | 1.0e-14 | n/a (not in hw_timings.txt) |
| 12 | 3185 | 91 | 7 | 9191 | 0.0923 | 0.0070 | 0.0034 | 0.1027 | 0.0489 | 0.2911 | 0.4427 | 145456 | 1.3e-14 | n/a (not in hw_timings.txt) |
| 14 | 5440 | 120 | 8 | 17816 | 0.1628 | 0.0100 | 0.0052 | 0.1780 | 0.0821 | 0.6230 | 0.8831 | 356052 | 4.6e-14 | n/a (not in hw_timings.txt) |
| 16 | 8721 | 153 | 9 | 31977 | 0.2977 | 0.0164 | 0.0077 | 0.3218 | 0.1282 | 1.1703 | 1.6202 | 787257 | 2.6e-13 | n/a (not in hw_timings.txt) |
| 20 | 19481 | 231 | 11 | 86779 | 0.7744 | 0.0343 | 0.0155 | 0.8242 | 0.3728 | 3.3913 | 4.5883 | 3064644 | 3.9e-11 | n/a (not in hw_timings.txt) |

(block 1 of 3 quoted; sums of t_b over S = 4..12 per block: block 0: 1.522 s, block 1: 1.343 s, block 2: 1.362 s)

(S = 1..12 as asked; S = 14, 16, 20 are extensions without Munich data: S = 20 is dim V = 19 481, 231 irreps, vectors of up to 86 779 entries, 3.06·10⁶ coefficients, 4.6 s. For S ≤ 2 Munich's C++ is faster in absolute terms (28 µs vs 0.7 ms, 0.86 ms vs 1.8 ms: Python overhead); the crossover is at S = 3.)

**Breakdown at S = 8 and S = 10** (quoted block): at S = 8 (dim 825, 45 irreps, n_W ≤ 1749) the exact generators take 22.7 ms, Racah–Speiser 2.0 ms and the 45 exact PRV spaces 1.2 ms: **t_a = 25.8 ms**; the float setup (metric, orthonormal matrices, dense blocks, 489 SVD left inverses) 16.0 ms and the reconstruction of all 45 vectors (15 625 coefficients) 51.4 ms: **t_b = 93.2 ms against Munich's 80.19 s (×860)**. At S = 10 (dim 1716, 66 irreps, n_W ≤ 4290): t_a = 57.9 ms (generators 52.0, Racah–Speiser 3.7, PRV 2.2), setup 30.1 ms, recon 139.6 ms, **t_b = 0.228 s** (52 090 coefficients; Munich cannot run it). The n_W-free part t_a is dominated by the O(dim V) construction of the GT generators, never by the PRV nullspaces (w ≤ ⌊S/2⌋ + 1: the slices are at most 6-dimensional at S = 10).

**Scaling**: Munich's total grows by ×5.2–7.0 per unit of S between S = 4 and 8 (0.0574 → 0.355 → 2.49 → 15.3 → 80.2 s; its cost is O(n_W³) per irrep with n_W,max ≈ S^4.3); ours by ×2.0 (S = 4 → 5) decreasing to ×1.4 (S = 11 → 12) (t_b: 8.8 ms at S = 4, 93 ms at S = 8, 0.23 s at S = 10, 0.44 s at S = 12, 0.88 s at S = 14, 4.6 s at S = 20); t_b ∝ entries ∝ S^{≈5.9} (580 → 52 090 → 3.06·10⁶). The time per coefficient is 15 µs at S = 4, 4.4 µs at S = 10 and 1.5 µs at S = 20 (batching improves with block size).

**Per irrep, Munich's slowest rows** (S = 8; ours WARM = after the per-λ setup of 16 ms, i.e. `prv_space` + `hw_from_slice` of that irrep only; Munich's per-μ column of `hw_timings.txt`, i5-1345U):

| mu (Munich label) | dim V_mu | n_W | Munich hw-solve [s] | ours: prv_space + hw_from_slice, warm [s] | speed-up |
|---|---|---|---|---|---|
| (16 16 0 0) | 8721 | 1 | 0.024730 | 0.00024 | 101 |
| (16 15 1 0) | 27455 | 2 | 0.029378 | 0.00024 | 122 |
| (16 14 2 0) | 47424 | 3 | 0.031028 | 0.00025 | 126 |
| (16 13 3 0) | 62700 | 4 | 0.032595 | 0.00025 | 131 |
| (16 12 4 0) | 69825 | 5 | 0.034100 | 0.00027 | 128 |
| (16 11 5 0) | 67431 | 6 | 0.035693 | 0.00026 | 137 |
| (16 10 6 0) | 55860 | 7 | 0.038969 | 0.00026 | 149 |
| (16 9 7 0) | 36784 | 8 | 0.040596 | 0.00027 | 152 |
| (16 8 8 0) | 12825 | 9 | 0.042272 | 0.00027 | 157 |
| (14 14 0 0) | 5440 | 8 | 0.040521 | 0.00030 | 135 |
| (14 13 1 0) | 16575 | 14 | 0.050295 | 0.00032 | 158 |
| (14 12 2 0) | 27489 | 20 | 0.060747 | 0.00036 | 167 |
| (14 11 3 0) | 34476 | 26 | 0.070540 | 0.00038 | 185 |
| (14 10 4 0) | 35700 | 32 | 0.082441 | 0.00041 | 203 |
| (14 9 5 0) | 30855 | 38 | 0.095320 | 0.00044 | 218 |
| (14 8 6 0) | 20825 | 44 | 0.103035 | 0.00046 | 222 |
| (14 7 7 0) | 7344 | 48 | 0.110981 | 0.00046 | 239 |
| (12 12 0 0) | 3185 | 34 | 0.085641 | 0.00040 | 215 |
| (12 11 1 0) | 9295 | 54 | 0.128170 | 0.00047 | 274 |
| (12 10 2 0) | 14580 | 74 | 0.158822 | 0.00054 | 292 |
| (12 9 3 0) | 16940 | 94 | 0.195808 | 0.00065 | 303 |
| (12 8 4 0) | 15625 | 114 | 0.237354 | 0.00070 | 340 |
| (12 7 5 0) | 10935 | 132 | 0.270702 | 0.00077 | 352 |
| (12 6 6 0) | 3920 | 142 | 0.293158 | 0.00079 | 370 |
| (10 10 0 0) | 1716 | 104 | 0.217479 | 0.00062 | 352 |
| (10 9 1 0) | 4719 | 154 | 0.313077 | 0.00076 | 414 |
| (10 8 2 0) | 6825 | 204 | 0.414750 | 0.00090 | 458 |
| (10 7 3 0) | 7020 | 252 | 0.516294 | 0.00104 | 494 |
| (10 6 4 0) | 5200 | 292 | 0.612051 | 0.00120 | 513 |
| (10 5 5 0) | 1911 | 312 | 0.658413 | 0.00128 | 516 |
| (8 8 0 0) | 825 | 259 | 0.535824 | 0.00103 | 519 |
| (8 7 1 0) | 2079 | 362 | 0.788352 | 0.00129 | 614 |
| (8 6 2 0) | 2640 | 457 | 1.056958 | 0.00155 | 684 |
| (8 5 3 0) | 2156 | 532 | 1.307627 | 0.00168 | 778 |
| (8 4 4 0) | 825 | 567 | 1.437635 | 0.00186 | 771 |
| (6 6 0 0) | 336 | 548 | 1.389312 | 0.00158 | 881 |
| (6 5 1 0) | 735 | 714 | 2.105674 | 0.00191 | 1101 |
| (6 4 2 0) | 729 | 840 | 2.840146 | 0.00233 | 1220 |
| (6 3 3 0) | 300 | 896 | 3.235421 | 0.00247 | 1313 |
| (4 4 0 0) | 105 | 980 | 3.986252 | 0.00235 | 1696 |
| (4 3 1 0) | 175 | 1176 | 6.065269 | 0.00273 | 2224 |
| (4 2 2 0) | 84 | 1260 | 7.219066 | 0.00315 | 2291 |
| (2 2 0 0) | 20 | 1464 | 10.971580 | 0.00363 | 3022 |
| (2 1 1 0) | 15 | 1584 | 13.520279 | 0.00434 | 3119 |
| (0 0 0 0) | 1 | 1749 | 18.702406 | 0.00454 | 4116 |

* sum Munich 80.187 s, sum ours (warm, excl. setup) 0.0520 s
* singlet (S,S,S,S) of [6,6,0,0]^2: nW = 588, cold (everything incl. generators, setup) 0.0176 s, warm (prv + recon only) 0.0020 s
* singlet (S,S,S,S) of [8,8,0,0]^2: nW = 1749, cold (everything incl. generators, setup) 0.0447 s, warm (prv + recon only) 0.0045 s
* singlet (S,S,S,S) of [10,10,0,0]^2: nW = 4290, cold (everything incl. generators, setup) 0.0954 s, warm (prv + recon only) 0.0094 s

(The singlet — Munich's slowest row, 18.7 s because n_W = 1749 is largest there — takes 4.5 ms warm (×4100) and 45 ms for a single COLD irrep including the generators and the setup (×420); the speed-up of the individual rows is between ×100 (the top irrep (16 16 0 0), where Munich's own solve is already trivial) and ×4100.)

**Cost of what t_a does not include** (`extras`): the six non-adjacent exact generators and the metric + all 12 float CSR matrices, needed by G2/G3 and by `hw_residual(extra=True)`:

* [4,4,0,0]: exact E_13,E_24,E_14,E_31,E_42,E_41 0.0023 s;  metric + all 12 float CSR matrices 0.0026 s
* [8,8,0,0]: exact E_13,E_24,E_14,E_31,E_42,E_41 0.0257 s;  metric + all 12 float CSR matrices 0.0171 s
* [10,10,0,0]: exact E_13,E_24,E_14,E_31,E_42,E_41 0.0583 s;  metric + all 12 float CSR matrices 0.0367 s

## 7. Multiplicity families (for G3's choice)

Racah–Speiser only (no vectors); `max w (site 1 = λ)` = max dim V_λ[μ−λ′] with site 1 = the first factor, `max w (site 1 = λ′)` with the factors exchanged (the exchange can matter: e.g. [3,2,1,0]⊗[2,1,0,0] has max w = 4 but [2,1,0,0]⊗[3,2,1,0] has 2 — put the factor with the smaller weight spaces at site 1 when λ ≠ λ′).

**[2a,a,0,0]^2**

| a | lam (x) lam' | dim | dim' | irreps | sum m | max m | max w (site 1 = lam) | max w (site 1 = lam') | max n_W | #irreps with max m | first irrep with max m |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | [2,1,0,0]^2 | 20 | 20 | 7 | 8 | 2 | 2 | 2 | 20 | 1 | [3,2,1,0] |
| 2 | [4,2,0,0]^2 | 126 | 126 | 24 | 32 | 3 | 4 | 4 | 174 | 1 | [6,4,2,0] |
| 3 | [6,3,0,0]^2 | 480 | 480 | 58 | 90 | 4 | 7 | 7 | 920 | 1 | [9,6,3,0] |
| 4 | [8,4,0,0]^2 | 1375 | 1375 | 115 | 205 | 5 | 11 | 11 | 3541 | 1 | [12,8,4,0] |
| 5 | [10,5,0,0]^2 | 3276 | 3276 | 201 | 406 | 6 | 15 | 15 | 10968 | 1 | [15,10,5,0] |
| 6 | [12,6,0,0]^2 | 6860 | 6860 | 322 | 728 | 7 | 20 | 20 | 29028 | 1 | [18,12,6,0] |
| 7 | [14,7,0,0]^2 | 13056 | 13056 | 484 | 1212 | 8 | 26 | 26 | 68208 | 1 | [21,14,7,0] |
| 8 | [16,8,0,0]^2 | 23085 | 23085 | 693 | 1905 | 9 | 33 | 33 | 146025 | 1 | [24,16,8,0] |

**[3a,2a,a,0]^2**

| a | lam (x) lam' | dim | dim' | irreps | sum m | max m | max w (site 1 = lam) | max w (site 1 = lam') | max n_W | #irreps with max m | first irrep with max m |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | [3,2,1,0]^2 | 64 | 64 | 16 | 30 | 4 | 4 | 4 | 152 | 1 | [5,4,2,1] |
| 2 | [6,4,2,0]^2 | 729 | 729 | 76 | 297 | 12 | 15 | 15 | 4743 | 1 | [9,7,5,3] |
| 3 | [9,6,3,0]^2 | 4096 | 4096 | 213 | 1592 | 27 | 32 | 32 | 59008 | 1 | [14,11,7,4] |
| 4 | [12,8,4,0]^2 | 15625 | 15625 | 459 | 5945 | 51 | 65 | 65 | 426425 | 1 | [19,14,10,5] |

**[2a,a,a,0]^2**

| a | lam (x) lam' | dim | dim' | irreps | sum m | max m | max w (site 1 = lam) | max w (site 1 = lam') | max n_W | #irreps with max m | first irrep with max m |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | [2,1,1,0]^2 | 15 | 15 | 6 | 7 | 2 | 3 | 3 | 21 | 1 | [3,2,2,1] |
| 2 | [4,2,2,0]^2 | 84 | 84 | 19 | 26 | 3 | 6 | 6 | 186 | 1 | [6,4,4,2] |
| 3 | [6,3,3,0]^2 | 300 | 300 | 44 | 70 | 4 | 10 | 10 | 1002 | 1 | [9,6,6,3] |
| 4 | [8,4,4,0]^2 | 825 | 825 | 85 | 155 | 5 | 15 | 15 | 3927 | 1 | [12,8,8,4] |
| 5 | [10,5,5,0]^2 | 1911 | 1911 | 146 | 301 | 6 | 21 | 21 | 12363 | 1 | [15,10,10,5] |
| 6 | [12,6,6,0]^2 | 3920 | 3920 | 231 | 532 | 7 | 28 | 28 | 33188 | 1 | [18,12,12,6] |
| 7 | [14,7,7,0]^2 | 7344 | 7344 | 344 | 876 | 8 | 36 | 36 | 78948 | 1 | [21,14,14,7] |
| 8 | [16,8,8,0]^2 | 12825 | 12825 | 489 | 1365 | 9 | 45 | 45 | 170829 | 1 | [24,16,16,8] |

**[a,a,a,0]x[a,a,0,0]**

| a | lam (x) lam' | dim | dim' | irreps | sum m | max m | max w (site 1 = lam) | max w (site 1 = lam') | max n_W | #irreps with max m | first irrep with max m |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | [1,1,1,0]x[1,1,0,0] | 4 | 6 | 2 | 2 | 1 | 1 | 1 | 3 | 2 | [2,2,1,0] |
| 2 | [2,2,2,0]x[2,2,0,0] | 10 | 20 | 3 | 3 | 1 | 1 | 1 | 6 | 3 | [4,4,2,0] |
| 3 | [3,3,3,0]x[3,3,0,0] | 20 | 50 | 4 | 4 | 1 | 1 | 1 | 10 | 4 | [6,6,3,0] |
| 4 | [4,4,4,0]x[4,4,0,0] | 35 | 105 | 5 | 5 | 1 | 1 | 1 | 15 | 5 | [8,8,4,0] |

**[a,a,0,0]x[a,0,0,0]**

| a | lam (x) lam' | dim | dim' | irreps | sum m | max m | max w (site 1 = lam) | max w (site 1 = lam') | max n_W | #irreps with max m | first irrep with max m |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | [1,1,0,0]x[1,0,0,0] | 6 | 4 | 2 | 2 | 1 | 1 | 1 | 3 | 2 | [2,1,0,0] |
| 2 | [2,2,0,0]x[2,0,0,0] | 20 | 10 | 3 | 3 | 1 | 1 | 1 | 6 | 3 | [4,2,0,0] |
| 3 | [3,3,0,0]x[3,0,0,0] | 50 | 20 | 4 | 4 | 1 | 1 | 1 | 10 | 4 | [6,3,0,0] |
| 4 | [4,4,0,0]x[4,0,0,0] | 105 | 35 | 5 | 5 | 1 | 1 | 1 | 15 | 5 | [8,4,0,0] |
| 5 | [5,5,0,0]x[5,0,0,0] | 196 | 56 | 6 | 6 | 1 | 1 | 1 | 21 | 6 | [10,5,0,0] |
| 6 | [6,6,0,0]x[6,0,0,0] | 336 | 84 | 7 | 7 | 1 | 1 | 1 | 28 | 7 | [12,6,0,0] |

**Reading for G3.** [2a,a,0,0]² and [2a,a,a,0]² have max m = a + 1 (linear growth, attained in a single irrep, (3a,2a,a,0) resp. (3a,2a,2a,a)), max w growing like ≈ a²/2 resp. (a+1)(a+2)/2; **[3a,2a,a,0]² is the fast family: max m = 4, 12, 27, 51 for a = 1..4** (a single irrep each time: (5,4,2,1), (9,7,5,3), (14,11,7,4), (19,14,10,5)) with 297 copies and m = 12 already at a = 2 (dim 729, w ≤ 15) and 1592 copies at a = 3 (dim 4096, w ≤ 32): that is the stress family for the missing-label operators; [2a,a,0,0]² is the scaling family (m = a + 1 up to 9 at a = 8, dim 23 085, the gl(4) analogue of [2a,a,0]). The two Pieri-type families ([a,a,a,0]⊗[a,a,0,0], [a,a,0,0]⊗[a,0,0,0]) are multiplicity-free (all m = 1, a + 1 irreps) and not useful for G3.

**t_a / t_b along the families** (all irreps, all m PRV basis vectors per irrep; for m ≥ 2 the PRV basis vectors are not yet the canonical copies — the recursion is linear, so the canonical copies cost the same once their slices are known; min of 5 cold runs):

**[2a,a,0,0]^2**

| a | lam | dim V | irreps | copies (= PRV basis vectors) | max m | max w | max n_W | t_a [s] | t_b [s] | entries | worst HW resid | breakdown t_gen / t_rs / t_prv / t_setup / t_recon [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | [2,1,0,0] | 20 | 7 | 8 | 2 | 2 | 20 | 0.0007 | 0.0019 | 59 | 3.5e-16 | 0.0005 / 0.0001 / 0.0001 / 0.0008 / 0.0004 |
| 2 | [4,2,0,0] | 126 | 24 | 32 | 3 | 4 | 174 | 0.0044 | 0.0103 | 1287 | 1.6e-15 | 0.0034 / 0.0004 / 0.0006 / 0.0026 / 0.0033 |
| 3 | [6,3,0,0] | 480 | 58 | 90 | 4 | 7 | 920 | 0.0163 | 0.0428 | 14208 | 6.3e-15 | 0.0129 / 0.0009 / 0.0024 / 0.0073 / 0.0192 |
| 4 | [8,4,0,0] | 1375 | 115 | 205 | 5 | 11 | 3541 | 0.0541 | 0.1515 | 101092 | 8.7e-15 | 0.0422 / 0.0023 / 0.0096 / 0.0195 / 0.0779 |
| 5 | [10,5,0,0] | 3276 | 201 | 406 | 6 | 15 | 10968 | 0.1459 | 0.4783 | 531186 | 1.6e-14 | 0.1049 / 0.0063 / 0.0347 / 0.0462 / 0.2862 |
| 6 | [12,6,0,0] | 6860 | 322 | 728 | 7 | 20 | 29028 | 0.3901 | 1.1084 | 2236946 | 1.3e-13 | 0.2915 / 0.0106 / 0.0879 / 0.0902 / 0.6281 |

**[3a,2a,a,0]^2**

| a | lam | dim V | irreps | copies (= PRV basis vectors) | max m | max w | max n_W | t_a [s] | t_b [s] | entries | worst HW resid | breakdown t_gen / t_rs / t_prv / t_setup / t_recon [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | [3,2,1,0] | 64 | 16 | 30 | 4 | 4 | 152 | 0.0020 | 0.0048 | 1069 | 2.5e-15 | 0.0016 / 0.0002 / 0.0003 / 0.0013 / 0.0015 |
| 2 | [6,4,2,0] | 729 | 76 | 297 | 12 | 15 | 4743 | 0.0295 | 0.0761 | 185606 | 6.5e-15 | 0.0208 / 0.0013 / 0.0074 / 0.0095 / 0.0371 |
| 3 | [9,6,3,0] | 4096 | 213 | 1592 | 27 | 32 | 59008 | 0.2347 | 0.7275 | 8889611 | 3.5e-14 | 0.1406 / 0.0070 / 0.0871 / 0.0588 / 0.4340 |

**[2a,a,a,0]^2**

| a | lam | dim V | irreps | copies (= PRV basis vectors) | max m | max w | max n_W | t_a [s] | t_b [s] | entries | worst HW resid | breakdown t_gen / t_rs / t_prv / t_setup / t_recon [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | [2,1,1,0] | 15 | 6 | 7 | 2 | 3 | 21 | 0.0005 | 0.0014 | 50 | 2.9e-15 | 0.0004 / 0.0001 / 0.0001 / 0.0006 / 0.0003 |
| 2 | [4,2,2,0] | 84 | 19 | 26 | 3 | 6 | 186 | 0.0029 | 0.0071 | 1021 | 1.0e-15 | 0.0023 / 0.0002 / 0.0004 / 0.0017 / 0.0025 |
| 3 | [6,3,3,0] | 300 | 44 | 70 | 4 | 10 | 1002 | 0.0108 | 0.0300 | 10978 | 4.5e-15 | 0.0087 / 0.0007 / 0.0015 / 0.0049 / 0.0143 |
| 4 | [8,4,4,0] | 825 | 85 | 155 | 5 | 15 | 3927 | 0.0304 | 0.0983 | 77366 | 4.7e-15 | 0.0242 / 0.0015 / 0.0047 / 0.0121 / 0.0557 |
| 5 | [10,5,5,0] | 1911 | 146 | 301 | 6 | 21 | 12363 | 0.0758 | 0.2768 | 405522 | 3.1e-14 | 0.0585 / 0.0033 / 0.0139 / 0.0263 / 0.1747 |

(block 0 of 3 quoted)

Even the stress case [9,6,3,0]² — dim 4096, 213 irreps, 1592 vectors, 8.9·10⁶ coefficients, n_W up to 59 008, m up to 27, w up to 32 — takes t_a = 0.23 s and t_b = 0.73 s (float residual 3.5e-14); exact PRV spaces with w = 32 and one left inverse per site-2 weight are not a bottleneck. **Exact vs float for G3**: the exact recursion + exact HW check of all 297 PRV basis vectors of [6,4,2,0]² takes 7.9 s (185 606 coefficients), one exact basis vector of the m = 27 sector (14,11,7,4) of [9,6,3,0]² 1.2 s (4158 entries), the float recursion for all 27 of them 0.066 s: use the exact path where exact vectors are wanted (it is not limited to dim·dim′ ≤ 2000), the float one everywhere else.

## 8. G1.6 Reference files and their independent re-verification

Records exactly as in the brief (`REP, REP2, DIM, DIM2, GTBASIS, GTBASIS2, WEIGHTS1, WEIGHTS2, GMET1, GMET2, NSECTORS, NCOPIES`; per sector in descending lexicographic μ `SECTOR mu m w nW dimVmu cas2mu`, `WSPACE`, `PRVBASIS`×m; for m = 1 `COPY mu 1 0 {{nu,val},..}`, `HWRESID`, `GNORM2`, and `VEC` when DIM·DIM2 ≤ 2000), TAB-separated, `#` header documenting every record and convention; every number in the files is an **exact rational** (the `COPY` slice too; the placeholder 0 is the eigenvalue slot for G3). Two deliberate additions: (a) `GNORM2` is the EXACT G-norm² for every file (also for the files without `VEC`: the exact recursion is cheap), (b) for m ≥ 2 sectors of the small products an extra record `PRVVEC mu k {{idx,val},..}` carries the exact full HW vector of PRV basis vector k (so that G3 can test its slice-reduced operators against full vectors without re-running the recursion). Files without `VEC` (DIM·DIM2 > 2000): [3,2,1,0]² (4096), [3,3,0,0]² (2500), [4,4,0,0]², [6,6,0,0]², [8,8,0,0]².

| file | sectors | sum m (copies) | m = 1 irreps (COPY) | exact full vectors written (VEC/PRVVEC) | exact conditions verified | exact HW violations | max HWRESID | GNORM2 float vs exact mismatches | bytes |
|---|---|---|---|---|---|---|---|---|---|
| gl4-cg-slices-2100.txt | 7 | 8 | 6 | 8 | 17 | 0 | 3.47e-16 | 0 | 9442 |
| gl4-cg-slices-2110.txt | 6 | 7 | 5 | 7 | 13 | 0 | 2.25e-15 | 0 | 8705 |
| gl4-cg-slices-3210.txt | 16 | 30 | 8 | 0 | 162 | 0 | 2.39e-15 | 0 | 13798 |
| gl4-cg-slices-2100-1100.txt | 4 | 4 | 4 | 4 | 4 | 0 | 1.48e-16 | 0 | 8009 |
| gl4-cg-slices-2210-2110.txt | 6 | 7 | 5 | 7 | 11 | 0 | 2.41e-15 | 0 | 8944 |
| gl4-cg-slices-2200.txt | 6 | 6 | 6 | 6 | 24 | 0 | 2.91e-16 | 0 | 9245 |
| gl4-cg-slices-3300.txt | 10 | 10 | 10 | 0 | 149 | 0 | 5.89e-16 | 0 | 11947 |
| gl4-cg-slices-4400.txt | 15 | 15 | 15 | 0 | 592 | 0 | 8.24e-16 | 0 | 18327 |
| gl4-cg-slices-6600.txt | 28 | 28 | 28 | 0 | 4736 | 0 | 1.77e-15 | 0 | 47003 |
| gl4-cg-slices-8800.txt | 45 | 45 | 45 | 0 | 22832 | 0 | 9.14e-15 | 0 | 114844 |

**`verify_ref(path)`** re-parses each file from the TEXT and checks it with code independent of the engine's PRV/recursion code (it uses `GLNRep.dense` dense float matrices — for all N² E_ij — plus float SVDs and the exact `GLNRep.apply`): GTBASIS/WEIGHTS/GMET/dims/sector set; sector fields (m = Racah–Speiser, dimVmu = Weyl, cas2mu, w, nW, WSPACE recomputed from the parsed weights); PRVBASIS: string conditions E_i^{k_i}v = 0 with dense matrix powers AND in exact arithmetic, rank m, and an independent float SVD of the string-condition matrix on the weight space whose null dimension must equal m; COPY = PRVBASIS/first nonzero (exact); VEC/PRVVEC: slice at q_top = COPY/PRVBASIS (exact), total weight, exact HW conditions for all E_ij (i < j), dense float HW residual, exact GNORM2; for the files without VEC the vector is reconstructed from the file's slice (float recursion) and checked with the dense independent matrices for all E_ij, plus GNORM2 vs the float value. **0 mismatches at 1e-12 in all ten files** (the exact checks are exact; the float ones are ≤ 9.1e-15):

| file | sectors | PRVBASIS vectors | COPY | full vectors checked (exact) | metric E^T G = G E (dense) dev | string conditions rel. dev (float) | float null dim != m | other field mismatches | exact HW violations | slice/weight mismatches | float HW residual (dense indep.) | GNORM2 vs float | TOTAL MISMATCHES |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| gl4-cg-slices-2100.txt | 7 | 8 | 6 | 8 | 0.0e+00 | 0.0e+00 | 0 | {} | 0 | 0 | 0.0e+00 | 0.0e+00 | 0 |
| gl4-cg-slices-2110.txt | 6 | 7 | 5 | 7 | 0.0e+00 | 0.0e+00 | 0 | {} | 0 | 0 | 6.6e-17 | 0.0e+00 | 0 |
| gl4-cg-slices-3210.txt | 16 | 30 | 8 | 0 | 0.0e+00 | 0.0e+00 | 0 | {} | 0 | 0 | 2.4e-15 | 4.0e-16 | 0 |
| gl4-cg-slices-2100-1100.txt | 4 | 4 | 4 | 4 | 0.0e+00 | 0.0e+00 | 0 | {} | 0 | 0 | 0.0e+00 | 0.0e+00 | 0 |
| gl4-cg-slices-2210-2110.txt | 6 | 7 | 5 | 7 | 0.0e+00 | 0.0e+00 | 0 | {} | 0 | 0 | 0.0e+00 | 0.0e+00 | 0 |
| gl4-cg-slices-2200.txt | 6 | 6 | 6 | 6 | 0.0e+00 | 0.0e+00 | 0 | {} | 0 | 0 | 0.0e+00 | 0.0e+00 | 0 |
| gl4-cg-slices-3300.txt | 10 | 10 | 10 | 0 | 0.0e+00 | 0.0e+00 | 0 | {} | 0 | 0 | 6.7e-16 | 4.3e-16 | 0 |
| gl4-cg-slices-4400.txt | 15 | 15 | 15 | 0 | 4.2e-25 | 0.0e+00 | 0 | {} | 0 | 0 | 8.5e-16 | 6.9e-16 | 0 |
| gl4-cg-slices-6600.txt | 28 | 28 | 28 | 0 | 1.5e-16 | 0.0e+00 | 0 | {} | 0 | 0 | 1.7e-15 | 1.1e-15 | 0 |
| gl4-cg-slices-8800.txt | 45 | 45 | 45 | 0 | 5.2e-20 | 0.0e+00 | 0 | {} | 0 | 0 | 9.1e-15 | 1.2e-15 | 0 |

**Mutation test of the verifier** (a verifier that cannot fail proves nothing): corrupting one record at a time in a copy of `gl4-cg-slices-2100.txt` — every corruption is detected, the pristine copy passes:

| corruption of a copy of gl4-cg-slices-2100.txt | verify_ref total mismatches | verdict |
|---|---|---|
| pristine copy | 0 | OK (0 expected) |
| PRVBASIS value changed (vector with >= 2 entries) | 3 | DETECTED |
| PRVBASIS scaled (m = 2 sector, pinned by PRVVEC) | 1 | DETECTED |
| COPY value changed | 2 | DETECTED |
| VEC value changed | 4 | DETECTED |
| GNORM2 changed | 1 | DETECTED |
| GMET1 entry changed | 1 | DETECTED |
| SECTOR multiplicity changed | 3 | DETECTED |
| WSPACE index changed | 2 | DETECTED |
| WEIGHTS1 entry changed | 1 | DETECTED |
| VEC index moved | 12 | DETECTED |
| PRVBASIS record dropped | 2 | DETECTED |

(An early version of the mutation test "failed" on rescaling a PRVBASIS vector with a single entry in an m = 1 sector: that is not a hole, the scale of a PRVBASIS vector is irrelevant there because COPY is its normalisation; in m ≥ 2 sectors the PRVVEC record pins the scale and the corruption is detected, second row.)

## 9. The Munich bridge: what could be settled without the binary, and what the Ubuntu laptop must run

The Munich code (`ClebschGordan.cpp`) works in an orthonormal GT basis; I re-implemented, from its source, (i) the pattern index (`clebsch::pattern::operator++`, i.e. Q(M)) and (ii) the lowering coefficients (`clebsch::pattern::lowering_coeff(k,l)`, Eq. 28) and compared them with the exact `GLNRep` (`glN_reps.munich_check`, also test (vii)):

| irrep | pattern order = Munich Q order | # lowering matrix elements compared | Munich lowering_coeff^2 != ours | extra nonzero entries of ours | gauge inconsistent | gauge trivial (all +1) |
|---|---|---|---|---|---|---|
| [1,0,0,0] | True | 3 | 0 | 0 | 0 | True |
| [1,1,0,0] | True | 6 | 0 | 0 | 0 | True |
| [2,1,0,0] | True | 33 | 0 | 0 | 0 | True |
| [2,1,1,0] | True | 22 | 0 | 0 | 0 | True |
| [3,2,1,0] | True | 146 | 0 | 0 | 0 | True |
| [2,0,0,0] | True | 12 | 0 | 0 | 0 | True |
| [1,1,0,0] | True | 6 | 0 | 0 | 0 | True |
| [2,2,0,0] | True | 32 | 0 | 0 | 0 | True |
| [3,3,0,0] | True | 100 | 0 | 0 | 0 | True |
| [4,4,0,0] | True | 240 | 0 | 0 | 0 | True |
| [5,5,0,0] | True | 490 | 0 | 0 | 0 | True |
| [6,6,0,0] | True | 896 | 0 | 0 | 0 | True |
| [7,7,0,0] | True | 1512 | 0 | 0 | 0 | True |
| [8,8,0,0] | True | 2400 | 0 | 0 | 0 | True |
| [10,10,0,0] | True | 5280 | 0 | 0 | 0 | True |

* comparer self-test: 19 synthetic hwcgc.sh outputs accepted with 0 failures; corrupted outputs rejected: 39 of 39

Consequences: **Q(M) = our GT index + 1 in every irrep; Munich's orthonormal basis vectors are ours rescaled by 1/√g (the diagonal sign gauge between the two bases is trivial: all plain-basis simple raising/lowering coefficients are ≥ 0, as Munich's √(…))**. The HW coefficients printed by `hwcgc.sh` (rows `Q(M)  Q(M')  Q(M'')=dim(μ)  CGC`) must therefore equal, up to ONE global sign σ (the LAPACK null vector's sign; for λ = λ′ the exchange of the two factors can only produce such a sign),

    CGC[Q(M), Q(M')] = σ · Ψ[ν,q] · sqrt(g_ν g_q / GNORM2),   ν = Q(M), q = Q(M'),

with Ψ the vector of the `VEC` record (or the recursion output) and `GNORM2` its exact G-norm² (so Σ CGC² = 1; the files give `GMET1/GMET2` and `GNORM2`). `munich_compare`/`munichcmp` does exactly this comparison from the text printed by `hwcgc.sh` (same support, equal values up to one global sign, tolerance 1e-9 as Munich prints 15 decimals); I tested it on synthetic outputs generated from our vectors with the sign flipped (19/19 accepted) and on corrupted outputs (value changed, index changed, row dropped: 39/39 rejected). **To run on the Ubuntu laptop** (after `bash build.sh` in `Paul/Cpp/AlexCGC`; the comparison script needs `slice_engine_gl4.py`, `glN_reps.py` and `Paul/Python/gl3_reps.py` in the repo layout; the Munich label of the irrep is the argument of `hwcgc.sh`, and `munichcmp` converts it):

    cd Paul/Cpp/AlexCGC
    bash hwcgc.sh 4 "S S 0 0" "S S 0 0" "m1 m2 m3 m4" > hw_S_m1_m2_m3_m4.txt
    cd ../../Python/gaudin_gl4          # with the output files copied here
    python -X utf8 slice_engine_gl4.py munichcmp S m1 m2 m3 m4 hw_S_m1_m2_m3_m4.txt

for all irreps of S = 1, 2, 3 (the 19 command pairs below are also printed by `python -X utf8 slice_engine_gl4.py munich`); the acceptance is "MATCH" on all 19, then S = 4..8 for a few irreps (e.g. the singlet `0 0 0 0` of S = 8, whose Munich solve takes 18.7 s; all irreps of [S,S,0,0]² are multiplicity-free, so `alpha` = 1 is the only copy):

    bash hwcgc.sh 4 "1 1 0 0" "1 1 0 0" "2 2 0 0" > hw_1_2_2_0_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 1 2 2 0 0 hw_1_2_2_0_0.txt
    bash hwcgc.sh 4 "1 1 0 0" "1 1 0 0" "2 1 1 0" > hw_1_2_1_1_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 1 2 1 1 0 hw_1_2_1_1_0.txt
    bash hwcgc.sh 4 "1 1 0 0" "1 1 0 0" "0 0 0 0" > hw_1_0_0_0_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 1 0 0 0 0 hw_1_0_0_0_0.txt
    bash hwcgc.sh 4 "2 2 0 0" "2 2 0 0" "4 4 0 0" > hw_2_4_4_0_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 2 4 4 0 0 hw_2_4_4_0_0.txt
    bash hwcgc.sh 4 "2 2 0 0" "2 2 0 0" "4 3 1 0" > hw_2_4_3_1_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 2 4 3 1 0 hw_2_4_3_1_0.txt
    bash hwcgc.sh 4 "2 2 0 0" "2 2 0 0" "4 2 2 0" > hw_2_4_2_2_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 2 4 2 2 0 hw_2_4_2_2_0.txt
    bash hwcgc.sh 4 "2 2 0 0" "2 2 0 0" "2 2 0 0" > hw_2_2_2_0_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 2 2 2 0 0 hw_2_2_2_0_0.txt
    bash hwcgc.sh 4 "2 2 0 0" "2 2 0 0" "2 1 1 0" > hw_2_2_1_1_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 2 2 1 1 0 hw_2_2_1_1_0.txt
    bash hwcgc.sh 4 "2 2 0 0" "2 2 0 0" "0 0 0 0" > hw_2_0_0_0_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 2 0 0 0 0 hw_2_0_0_0_0.txt
    bash hwcgc.sh 4 "3 3 0 0" "3 3 0 0" "6 6 0 0" > hw_3_6_6_0_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 3 6 6 0 0 hw_3_6_6_0_0.txt
    bash hwcgc.sh 4 "3 3 0 0" "3 3 0 0" "6 5 1 0" > hw_3_6_5_1_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 3 6 5 1 0 hw_3_6_5_1_0.txt
    bash hwcgc.sh 4 "3 3 0 0" "3 3 0 0" "6 4 2 0" > hw_3_6_4_2_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 3 6 4 2 0 hw_3_6_4_2_0.txt
    bash hwcgc.sh 4 "3 3 0 0" "3 3 0 0" "6 3 3 0" > hw_3_6_3_3_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 3 6 3 3 0 hw_3_6_3_3_0.txt
    bash hwcgc.sh 4 "3 3 0 0" "3 3 0 0" "4 4 0 0" > hw_3_4_4_0_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 3 4 4 0 0 hw_3_4_4_0_0.txt
    bash hwcgc.sh 4 "3 3 0 0" "3 3 0 0" "4 3 1 0" > hw_3_4_3_1_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 3 4 3 1 0 hw_3_4_3_1_0.txt
    bash hwcgc.sh 4 "3 3 0 0" "3 3 0 0" "4 2 2 0" > hw_3_4_2_2_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 3 4 2 2 0 hw_3_4_2_2_0.txt
    bash hwcgc.sh 4 "3 3 0 0" "3 3 0 0" "2 2 0 0" > hw_3_2_2_0_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 3 2 2 0 0 hw_3_2_2_0_0.txt
    bash hwcgc.sh 4 "3 3 0 0" "3 3 0 0" "2 1 1 0" > hw_3_2_1_1_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 3 2 1 1 0 hw_3_2_1_1_0.txt
    bash hwcgc.sh 4 "3 3 0 0" "3 3 0 0" "0 0 0 0" > hw_3_0_0_0_0.txt    && python -X utf8 slice_engine_gl4.py munichcmp 3 0 0 0 0 hw_3_0_0_0_0.txt

## 10. Surprises, pitfalls, inaccuracies in the brief

* **Everything is far faster than the plan's guess** ("seconds to minutes for all vectors at dim ~10⁴"): all vectors of [10,10,0,0]² in 0.23 s, of [20,20,0,0]² (dim 19 481) in 4.6 s, of [9,6,3,0]² (1592 vectors, m up to 27) in 0.73 s; the rectangular family has tiny slices (w_max = ⌊S/2⌋ + 1 = 6 at S = 10) while n_W grows like S^4.3 (4290 at S = 10): the exact PRV null space is trivial and the whole cost is Python overhead + output size.
* **The dangling conditions are necessary for N = 4** (§4.3): the exact twin of the gl(3) engine (T-conditions only) would accept slices outside M_μ here (15% of my random out-of-M_μ slices; the float residual then flags them, ≥ 3.7e-3).
* **The gap of AᵀA is not 1** for gl(4) (§5): 2 for [1,0,0,0]², then 0.586, 0.438, 0.354, 0.298, … for [S,S,0,0]², S = 1, 2, 3, 4, …, decreasing like 1/S (0.156 at S = 10).
* **Munich's GT basis = ours** (§9): pattern order, orthonormal lowering coefficients (6024 matrix elements in test (vii), 17 irreps, exact) and trivial sign gauge — so the deferred literal comparison is a formality up to ONE global sign.
* **Float64 recursion residual** grows ×2.4–3.7 per unit of S beyond S = 14 (§4.4); the exact-rounded path costs 1–16 s per vector.
* **Choice of site 1 for λ ≠ λ′** changes w (§7): [3,2,1,0]⊗[2,1,0,0] w = 4, reversed w = 2.
* Inaccuracies in the brief (harmless): "n_W ≤ dim V_λ = 336 at S = 5" — dim V_{[5,5,0,0]} = 196 and max n_W = 308 (336 is the dimension of the SU(4) irrep (6 6 0 0)); `hw_residual(hw)` for gl(4) has N−1 = 3 simple raising operators, I report all six E_ij (i < j) in the tables; the brief's candidate set μ = λ + wt(V_λ′) is only half of the μ for which the slice space is non-empty (I test the union).
* **Tooling trap** (as in the gl(3) report): inside the Bash tool's heredocs `\t`, `\n` in Python source are un-escaped; the files with such literals were written with the file tools (I fell into it twice and repaired both).

## 11. Not done / limits

* **The literal comparison with Munich's HW vectors** (binary not buildable here): §9 gives the exact commands, the conversion and the comparer; the GT basis identity is already verified.
* Float64 residual > 1e-12 beyond S ≈ 17 (cure: `hw_from_slice_rounded`, 1–16 s per vector; not batched).
* t_a does not include the six non-adjacent generators (lazy; cost in §6) — they are not needed by the slice engine, only by tests, `hw_residual(extra=True)` and G2/G3.
* Only N = 3 and N = 4 are tested at scale; N = 2, 5, 6 and λ₄ > 0 only as smoke tests (PRV = RS on all candidates, float + exact recursion, exact HW check for all E_ij; N = 5: 9 products, section "generality" of the output file, with negative tests and brute force).
* No vertex operator / missing-label operator (WP-G2/G3); the PRV basis for m ≥ 2 is the canonical RREF basis, not the canonical copies (`COPY` has the eigenvalue placeholder 0 and exists for m = 1 only).
* No parallelism; Windows only; Python Fractions (exact twin) are the only exact arithmetic.
* CLAUDE.md and the memory still say "gl(4) dropped on 2026-10-16" (the plan asks to update both when this starts); I did not touch them.

## 12. Usage

```python
import sys; sys.path.insert(0, 'Paul/Python/gaudin_gl4')            # gl3_reps is picked up from '..'
import glN_reps as gr, slice_engine_gl4 as se
R  = gr.GLNRep.of((3, 2, 1, 0))                    # exact sparse GT irrep: R.E[(i,j)] (columns of {row: Fraction}), R.g, R.wsp, R.top, R.float_ops(), R.dense((i,j)), R.apply((i,j), {idx: Fraction})
lam = lamp = (3, 2, 1, 0); mu = (5, 4, 2, 1)       # site 1 = V_lam (slow), site 2 = V_lamp (fast), slice = Psi[., top] in V_lam[mu - lamp]
p  = se.prv_space(lam, lamp, mu)                   # p.W (global GT indices of V_lam[mu-lamp]), p.basis (exact, aligned with p.W), p.free, p.m = 4, p.w, p.vec(k), p.coords(v), p.float_basis()
hw = se.hw_from_slice(lam, lamp, mu, p.float_basis())    # float recursion, all 4 vectors at once: hw.rows (site-1 index), hw.cols (site-2 index), hw.vals (nW, 4) in the plain GT basis
res, wdev = se.hw_residual(hw, extra=True)         # max ||E_ij^tot Psi||/||Psi|| over all i < j, and the weight deviation (0)
psi, nchk, viol = se.hw_from_slice_exact(lam, lamp, mu, p.basis[0])   # exact twin: {(nu, q): Fraction}, #exact conditions verified, [T, dangling] violation counts (0, 0)
hr = se.hw_from_slice_rounded(lam, lamp, mu, p.basis[0])               # exact recursion rounded to double: the accurate path for large S
D  = se.parse_ref('Paul/Mathematica/Gaudin/ref/gl4-cg-slices-3210.txt')  # D['sectors'][mu] -> m, w, nW, wspace, prv, copy, hwres, gnorm2, vec, prvvec (1-based indices, exact Fractions)
se.bruteforce_check(lam, lamp)                     # Munich-style null space check (no PRV, no recursion)
```
For G3: `se.Recursion.of(lam, lamp)` holds the left inverses (`.plan`), the dense raising blocks (`.B1`, `.B2`) and the site-2 labelling (`.ab2`, `.order`); `se.ExactRecursion.of(lam, lamp)` the exact twin; both accept λ ≠ λ′ (site 1 = λ).
