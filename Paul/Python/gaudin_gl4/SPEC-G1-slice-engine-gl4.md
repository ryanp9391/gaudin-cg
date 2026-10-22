# Spec G1: `glN_reps.py` and the gl(4) slice engine — GT representations of gl(N), PRV slice spaces with three string conditions, the level recursion, the SU(4) benchmark against Munich (2026-10-22)

Delegated to a general-purpose subagent (SONNET), **Python only** (no Wolfram kernel, no wolfbook tools, no git, no notebooks). This is work package
WP-G1 of `Paul/Mathematica/Gaudin/gl4-slice-engine-plan.md` (read its §0–§2 first: it says what carries over from gl(3) and what changes).
Environment: Python 3.13.13 (`python -X utf8`), numpy 2.5.2, scipy 1.18.1, sympy 1.14, mpmath 1.3 (no gmpy2, no `uv`); i7-12700KF (20 threads), 32 GB, Windows 11.
Use `time.perf_counter()`. ≤ 2 heavy processes at a time. **Never edit existing files; import them read-only.**

**Directory:** `Paul/Python/gaudin_gl4/` (new; this file lives there). Run everything from that directory. Files you own:
`glN_reps.py`, `test_glN_reps.py`, `slice_engine_gl4.py`, `RESULTS-G1-slice-engine-gl4.md`, `slice_engine_gl4_output.txt`,
and the reference files `../../Mathematica/Gaudin/ref/gl4-cg-slices-*.txt`. Caches (`*.pkl`, `*.npz`) are gitignored — fine to write, never required.

Another Sonnet agent (WP-G2, `SPEC-G2-gaudin-gl4-ops.md`, same directory) is working **concurrently** on the two-site gl(4) Gaudin operators. It may
import your `glN_reps.py` once it exists; it does not depend on you otherwise. Do not touch its files (`gaudin_gl4_ops.py`, `RESULTS-G2-*`).

## 0. Why this exists

For gl(3) the highest-weight (HW) Clebsch–Gordan problem of V_λ⊗V_λ′ collapses to a w-dimensional slice (w = a weight multiplicity of V_λ, ≤ 8 at
[14,7,0]²) and the gl(3) engine `Paul/Python/gaudin_sl3/slice_engine.py` does all 113 irreps of [14,7,0]² in 0.1 s against 206 s for the general-N
null-space code of Munich (`Paul/Cpp/AlexCGC`, arXiv:1009.0437). The theorem behind it (Kostant 1959 / PRV 1967) is general: for gl(N) the slice
space has N−1 string conditions. For N = 4 the only general code is Munich's, and its own SU(4) benchmark exists: `Paul/Cpp/AlexCGC/results/hw_timings.txt`
([S,S,0,0]², S = 1..8; S = 8: 45 irreps, 80.19 s total hw-solve, singlet 18.7 s; the code cannot run S = 10 on 16 GB). The family [S,S,0,0]² is
**multiplicity-free**, so for it the slice alone (a 1-dimensional null space per irrep, no missing-label operator) gives every HW vector. **This work
package is the SU(4) headline:** the gl(N) GT machinery, the PRV slice engine for gl(4) with its tests, the level recursion for the full vectors, the
benchmark table against Munich, and reference files for the later packages (G2/G3 add the Gaudin missing-label operators on the multiplicity spaces).

## 1. Assets (read before writing)

* `Paul/Python/gl3_reps.py` — GT patterns and the GT formulas `_raise(k, pattern)`, `_lower(k, pattern)`, `_weight(k, pattern)`, `gt_patterns(lam)`.
  **They are already general-N** (they only use the row lengths). Verified today by the main session: building E_{k,k+1}, E_{k+1,k} from them (dropping
  targets outside the pattern set, exactly as `GL3Rep._adjacent_up/_down` do), E_kk from `_weight`, and the non-adjacent generators by commutators,
  ALL 256 relations [E_ij, E_kl] = δ_jk E_il − δ_li E_kj hold exactly for the gl(4) irreps [1,0,0,0], [1,1,0,0], [2,1,0,0], [2,1,1,0], and the HW
  pattern (rows λ[:k]) is the LAST index. Reuse these four functions read-only.
* `Paul/Python/gaudin_sl3/slice_engine.py` — the gl(3) engine to generalise (1955 lines; read at least the module docstring, §1 `Rep` (exact sparse
  representation, metric `g`, `float_ops`, `block_exact`), §2 Racah–Speiser, §3 `prv_space`/`_rref_nullspace`/`PRV`, §5 `Recursion`/`hw_from_slice`/
  `ExactRecursion`/`hw_from_slice_exact`/`hw_residual`, §8 `bench_one`, §9 `write_ref`, §10 `verify_ref`/`bruteforce_check`). Its brief
  `SPEC-slice-engine.md` and report `RESULTS-slice-engine.md` (§0 conventions, the pitfalls: GT-order traps, float64 recursion residual growing like
  2.2^a and the `hw_from_slice_rounded` cure, timing methodology) are the model for yours. Import it read-only for the N = 3 regression tests
  (`sys.path.insert(0, '../gaudin_sl3')`; it inserts its own paths for `gl3_reps`).
* The reference-file format: `Paul/Mathematica/Gaudin/ref/cg-slices-420.txt` (read its `#` header: every record and convention is documented there).
* Munich: `Paul/Cpp/AlexCGC/README.md`, `results/hw_timings.txt` (per-μ tables with dim, mult, hw-solve time), `results/summary.csv`. The binary is
  **not** available on this Windows machine (no compiler, no WSL): the literal comparison of HW vectors with Munich is deferred to the Ubuntu laptop;
  you replace it by the Munich-style brute-force null space in Python (§3, G1.5) and you write in RESULTS the exact `hwcgc.sh` commands to run there.

## 2. Conventions (binding; the same as gl(3) wherever gl(3) had a convention)

* **GT patterns for gl(N):** tuple of N rows, `rows[0]` the singleton, `rows[N−1] = λ` (the `gl3_reps` convention). Order of the basis = the order of
  `gl3_reps.gt_patterns(lam)` (itertools.product order; for N = 3 identical to `GL3Rep(lam).patterns`); **the HW state is the LAST index** (assert it).
  All indices 0-based in Python, 1-based in the reference files.
* **Generators:** E_{k,k+1} from `_raise(k, ·)`, E_{k+1,k} from `_lower(k, ·)` (targets outside the pattern set dropped), E_kk = diag(`_weight(k, ·)`),
  E_ij for j > i+1 by E_ij = [E_{i,j−1}, E_{j−1,j}], E_ji = [E_{j,j−1}, E_{j−1,i}] (for N = 3 this is E₁₃ = [E₁₂,E₂₃], E₃₁ = [E₃₂,E₂₁], the gl(3)
  convention; since the GT matrices form a representation the commutators ARE the generators, so any path gives the same matrix — test it).
  Plain rational non-unitary GT normalisation (Fractions). Matrices act on kets as columns, E e_col = Σ_row M[row,col] e_row.
  Weights = (E₁₁,…,E_NN) eigenvalues; simple roots α_i = e_i − e_{i+1}.
* **Invariant metric:** the diagonal G with E_ijᵀ G = G E_ji, G[HW] = 1 (BFS from the top through the simple raising operators, as `Rep.g`); check it
  on every off-diagonal E_ij. Orthonormal basis ẽ_q = e_q/√g_q for the float recursion (Ê[r,q] = E[r,q]√(g_r/g_q), Ê_ijᵀ = Ê_ji).
* **Two sites:** Ψ = Σ Ψ[ν,q] e_ν⊗e_q, **site 1 = V_λ (ν, slow), site 2 = V_λ′ (q, fast)**, flat index ν·dim(V_λ′) + q. The **slice** of a vector of
  weight μ is ψ_top = Ψ[·, q_top], q_top = HW state of V_λ′ (last index), a vector in the weight space V_λ[μ − λ′] of site 1.
* **Theorem (Kostant / PRV):** Ψ ↦ ψ_top is an isomorphism of the HW vectors of weight μ in V_λ⊗V_λ′ onto
  **M_μ = { v ∈ V_λ[μ−λ′] : E_{i,i+1}^{λ′_i − λ′_{i+1} + 1} v = 0 for i = 1,…,N−1 }**, so dim M_μ = multiplicity of V_μ in V_λ⊗V_λ′.
  Note the roles: the string exponents come from λ′ (site 2, the site whose HW state defines the slice), the weight space lives in V_λ (site 1).
  A constraint whose image weight μ − λ′ + k α_i is not a weight of V_λ is void (as in gl(3)).
* **Level recursion (general N):** label the weights of site 2 by a = (a₁,…,a_{N−1}) ∈ ℤ≥0^{N−1}, wt(q) = λ′ − Σ a_i α_i; Q(a) = the GT states of V_λ′ of that
  weight; X_a = (ψ_q)_{q∈Q(a)} with ψ_q ∈ V_λ[μ − wt(q)]. The HW conditions (E_i⊗1 + 1⊗E_i)Ψ = 0, E_i = E_{i,i+1}, read at a site-2 component
  r ∈ Q(a − e_i):  E_i ψ_r + Σ_{q∈Q(a)} (E_i)_{rq} ψ_q = 0.  Hence for a ≠ 0:
  **A_a X_a = −[ E_i^{(site 1)} X_{a−e_i} ]_{i=1..N−1} (stacked over the i with a − e_i ∈ weights), A_a = [ (E_i^{(site 2)})_{Q(a−e_i), Q(a)} ]_i stacked.**
  A_a has full column rank (a vector of V_λ′ killed by all simple raising operators is the HW vector, and a ≠ 0) and is independent of μ and of the slice:
  one left inverse per (λ′, a) (SVD pseudo-inverse in the orthonormal basis; assert the smallest singular value is not tiny). Process a in lexicographic
  order (every a − e_i precedes a). The exact twin (Fractions) computes the left inverse by elimination and verifies the leftover consistency conditions
  T·RHS = 0 — they are the PRV conditions; a slice not in M_μ must fail them.
* **Racah–Speiser for gl(N):** ρ = (N−1, N−2, …, 0); mult(μ) = Σ_{w ∈ wt(V_λ′)} m_λ′(w)·sgn(s) over the w for which λ + w + ρ has distinct entries, s the
  permutation sorting it decreasingly, μ = s(λ+w+ρ) − ρ (as `slice_engine.racah_speiser`, which sums over the weights of the first argument and shifts the
  second — either way is correct; keep one and say which). Checks: Σ_μ m_μ dim V_μ = dim V_λ · dim V_λ′ for every product; no negative multiplicities.
  Weyl dimension formula for gl(N): dim V_λ = Π_{i<j} (λ_i − λ_j + j − i)/(j − i). Quadratic Casimir C₂ = Σ_ij E_ij E_ji = Σ_i λ_i(λ_i + N + 1 − 2i) (scalar
  on V_λ; for N = 3 this is W4's cas2 = Σ μ_i(μ_i + 4 − 2i)).
* **Munich's SU(4) labels** are our μ minus μ₄·(1,1,1,1) (e.g. the singlet of [1,1,0,0]² is our (1,1,1,1), Munich's (0 0 0 0)); their dims and multiplicities
  must be reproduced exactly; # irreps of [S,S,0,0]² = (S+1)(S+2)/2 = 3, 6, 10, 15, 21, 28, 36, 45 for S = 1..8.
* **Normalisation of slices** in reference files: first nonzero component (ascending GT index) = 1; the full vector carries the same scale. Compare vectors
  only up to scale (the GT basis is not orthonormal).

## 3. Tasks and acceptance (every number goes into `RESULTS-G1-slice-engine-gl4.md`)

**G1.1 `glN_reps.py`** — `class GLNRep(lam)` for any N (cache `GLNRep.of(lam)`): exact sparse generators `E[(i,j)]` (list over columns of dicts {row: Fraction},
the `slice_engine.Rep` layout, so that the recursion code ports directly), `pats`, `dim`, `index`, `wt`, `wsp` (weight → array of indices), `top`, metric `g`,
`check_metric()`, `apply(ij, vec)`, `float_ops()` (CSR in the orthonormal basis incl. diagonal), `block_exact(ij, target_w)`, `dense(ij)` (numpy float) for tests;
module functions `weyl_dim(lam)`, `weight_mults(lam)`, `racah_speiser(lam, lamp)`, `candidates(lam, lamp)` (all dominant μ = λ + w, w ∈ wt(V_λ′)), `casimir2(lam)`.
`test_glN_reps.py` (pytest-style functions AND runnable as `python -X utf8 test_glN_reps.py`): (i) N = 3: all nine E_ij equal `slice_engine.Rep.of(lam).E` entry for
entry and `GL3Rep(lam).generator(i,j)` for the adjacent ones, for [2,1,0], [3,1,0], [3,2,0], [4,2,0]; (ii) N = 4: all 256 commutation relations exact for
[1,0,0,0], [1,1,0,0], [2,0,0,0], [2,1,0,0], [2,1,1,0], [2,2,0,0], [3,2,1,0], [1,1,1,0], [3,3,0,0]; N = 5: [1,0,0,0,0], [2,1,0,0,0] (the code is general);
(iii) dims = Weyl formula, HW last, weight multiset permutation-symmetric, Σ_w m(w) = dim; (iv) metric: 0 violations; (v) C₂ scalar = formula; (vi) Racah–Speiser:
Σ m dim = dim·dim′ on every product of §G1.3 and the Munich table reproduced for S = 1..8. **Acceptance: all tests pass; print a one-line summary per test.**

**G1.2 `prv_space(lam, lamp, mu)`** in `slice_engine_gl4.py` — exact basis of M_μ: RREF nullspace over ℚ of the stacked power images E_{i,i+1}^{λ′_i−λ′_{i+1}+1}
on V_λ[μ−λ′] (canonical basis: one vector per free column, 1 there, 0 at the other free columns, so coordinates = free entries; `PRV` object with `.W`, `.basis`,
`.free`, `.m`, `.w`, `.vec(k)`, `.coords(v)`). **Acceptance:** dim M_μ = Racah–Speiser multiplicity for EVERY dominant candidate μ (including those with
multiplicity 0, where M_μ must be {0}) of: [1,0,0,0]², [1,1,0,0]², [2,1,0,0]², [2,1,1,0]², [2,2,0,0]², [3,2,1,0]², [2,2,1,0]², [3,3,0,0]²,
and λ ≠ λ′: [2,1,0,0]⊗[1,1,0,0], [1,1,0,0]⊗[2,1,0,0], [3,1,0,0]⊗[2,2,0,0], [2,2,1,0]⊗[2,1,1,0], [3,2,1,0]⊗[2,1,0,0] — in BOTH orders where λ ≠ λ′ (the
slice lives on different sites; both must give the multiplicity); then [S,S,0,0]², S = 1..10 (all m = 1), and the N = 3 regression: on [4,2,0]², [6,3,0]², [3,2,0]⊗[3,1,0]
the PRV bases coincide with `slice_engine.prv_space` exactly (same W, same basis; for λ ≠ λ′ `slice_engine` has no routine — just the multiplicity).
Table per product: #candidates, #irreps, Σm, max m, max w, max n_W, mismatches (= 0).

**G1.3 `hw_from_slice(lam, lamp, mu, v)`** — float level recursion in the orthonormal basis (left inverses by SVD, batched over k slices: v of shape (w,) or (w,k)),
returning COO data in the plain GT basis (rows = site-1 index, cols = site-2 index, vals (nW, k)); `hw_from_slice_exact` (Fractions, exact consistency conditions,
returns the number of conditions verified); `hw_from_slice_rounded` (exact → float, the accurate path); `hw_residual(hw)` = max over i of ‖E_{i,i+1}^{tot}Ψ‖/‖Ψ‖
(also check E_{i,j}^{tot} for all i < j, and the weight). **Acceptance:** relative HW residual ≤ 1e-12 and weight exactly μ for ALL irreps of all products in G1.2
(float path), the exact twin gives exactly 0 (E_ij^{tot}Ψ = 0 for every i < j, exact) on every irrep of the products with dim·dim′ ≤ 2000 and on ≥ 5 irreps each of
[3,2,1,0]² and [4,4,0,0]²; ≥ 20 random slices NOT in M_μ (random vectors of the weight space, m < w sectors) are rejected by the exact consistency conditions
(negative tests); [8,8,0,0]²: all 45 irreps, residuals ≤ 1e-12 incl. the singlet (8,8,8,8), the top (16,16,0,0), (10,5,5,0)+8·… — whatever Munich's slowest rows are.

**G1.4 Munich-style brute force `bruteforce_check(lam, lamp)`** — in the orthonormal two-site basis, the null space of the stacked total simple raising
operators E_{i,i+1}^{tot} restricted to the weight-μ subspace (n_W columns; dense SVD or sparse eigen-decomposition of AᵀA as `slice_engine.bruteforce_check`):
null dimension = multiplicity; our reconstructed HW vectors (orthonormalised) lie in it and span it (principal angles ≤ 1e-10). Do it for every irrep of every
product of G1.2 with dim·dim′ ≤ 10⁴ and for [S,S,0,0]², S ≤ 5 (n_W ≤ dim V_λ = 336 at S = 5). **Acceptance: 0 mismatches.** This replaces the Munich binary
(same mathematics: Munich solves exactly this null space).

**G1.5 Benchmark** (`bench`; min of 5 cold runs, caches cleared between runs, single process, `perf_counter`): [S,S,0,0]², S = 1..10 (go to 12 if S = 10 takes
< 2 min): (a) t_a = the n_W-free part (GT generators + RS + all PRV spaces, exact), (b) t_b = additionally metric/float matrices/left inverses + every full HW vector
by the float recursion (batched). Table: S, dim V, #irreps, max w, max n_W, t_a, t_b, entries written, worst HW residual, and Munich's total hw-solve for S ≤ 8
(0.0000279, 0.000858, 0.0107, 0.0574, 0.355, 2.49, 15.3, 80.2 s; machine i5-1345U — different from ours; say so). Also a timing breakdown at S = 8 and S = 10.
Then **multiplicity families** (for G3's choice; cheap: Racah–Speiser only, no vectors): tables of (#irreps, Σm, max m, max w) vs a for [2a,a,0,0]² (a = 1..5),
[3a,2a,a,0]² (a = 1..3), [2a,a,a,0]² (a = 1..4), [a,a,a,0]⊗[a,a,0,0] (a = 1..4), [a,a,0,0]⊗[a,0,0,0] (a = 1..6); plus t_a/t_b for [2a,a,0,0]², a = 1..4
(all irreps, all copies of M_μ — for m ≥ 2 the PRV basis vectors are not yet the canonical copies; reconstruct the basis vectors anyway: the recursion is linear).

**G1.6 Reference files** `Paul/Mathematica/Gaudin/ref/gl4-cg-slices-<l1l2l3l4>[-<l1'l2'l3'l4'>].txt` (e.g. `gl4-cg-slices-2100.txt` for [2,1,0,0]²,
`gl4-cg-slices-2100-1100.txt` for [2,1,0,0]⊗[1,1,0,0]) for: [2,1,0,0]², [2,1,1,0]², [3,2,1,0]², [2,1,0,0]⊗[1,1,0,0], [2,2,1,0]⊗[2,1,1,0], [2,2,0,0]², [3,3,0,0]²,
[4,4,0,0]², [6,6,0,0]², [8,8,0,0]² (the last two: slices only, no VEC). TAB-separated, `#` header documenting every record and convention (copy the style of
`cg-slices-420.txt`; state the GT order, the slice definition, the theorem, the site roles, the metric, the normalisation, N = 4, script name, date). Records:
`REP`, `REP2`, `DIM`, `DIM2`, `GTBASIS` (site-1 patterns in order, as nested braces {{k},{m1,m2},{..},{l1,l2,l3,l4}}), `GTBASIS2`, `WEIGHTS1`, `WEIGHTS2`, `GMET1`, `GMET2`,
`NSECTORS`, `NCOPIES`; per sector (descending lexicographic μ): `SECTOR mu m w nW dimVmu cas2mu`, `WSPACE mu {nu,...}` (1-based site-1 indices of V_λ[μ−λ′]),
`PRVBASIS mu {{nu,val},...}` × m (exact rationals), and for m = 1 `COPY mu 1 0 {{nu,val},...}` (the normalised slice; the 0 is a placeholder for the eigenvalue
that G3 will fill), `HWRESID mu 1 residual`, `GNORM2 mu 1 value` (exact G-norm² of the full vector, G = GMET1⊗GMET2), and `VEC mu 1 {{idx,val},...}` (exact
rationals, idx = (ν−1)·DIM2 + q, 1-based) when DIM·DIM2 ≤ 2000. `verify_ref(path)`: re-parse from the TEXT and check with code paths independent of the engine
(dense float matrices from `GLNRep.dense`, the string conditions, the HW conditions on VEC, GNORM2): **0 mismatches at 1e-12**.

**G1.7 `RESULTS-G1-slice-engine-gl4.md`** (style of `RESULTS-slice-engine.md`): files, conventions actually used, the tables of G1.1–G1.6 with numbers, the
benchmark table against Munich, what surprised you, what is not done (the literal Munich comparison: give the `bash hwcgc.sh 4 "S S 0 0" "S S 0 0" "mu"` commands
and the conversion needed — Munich uses the orthonormal GT basis, so our COPY/VEC must be multiplied by √(GMET) per site before comparison up to scale — so the
next Ubuntu session can run it). Save the printed output of all runs to `slice_engine_gl4_output.txt`. CLI: `python -X utf8 slice_engine_gl4.py <stage>` with
stages `tests | prv | hw | neg | brute | bench [Smax] | families | ref [all | l1 l2 l3 l4 [l1' l2' l3' l4']] | verify | all`. Report failures with numbers; never
paper over a failing test. Keep the code general in N (it must also run N = 3 and N = 5) but only N = 3 and 4 are tested at scale.

## 4. Pitfalls known from gl(3) (RESULTS-slice-engine.md)

* The raw GT formulas propose targets outside the interlacing range with nonzero coefficient — drop them (as `GL3Rep` does), do not clip.
* The float64 recursion's HW residual grows with the rep size (2.2^a for [2a,a,0]); at [8,8,0,0]² expect ~1e-13..1e-11; if a residual exceeds 1e-12, use
  `hw_from_slice_rounded` for that sector and say so. Normalise residuals relative to ‖Ψ‖ in the orthonormal basis.
* Orthonormal vs plain basis: all theorem statements are in the plain GT basis; the recursion runs in the orthonormal one; `HW.vals` converts back.
* Weight-space dictionaries keyed by tuples: build `wsp` once per rep; `nW` of a sector = Σ_{q} dim V_λ[μ − wt(q)] over the site-2 states q.
* Exact arithmetic with Fractions is fine up to w ~ 100; the PRV nullspace is w×w-sized; the exact recursion is only for small products.
* Timing on Windows: `perf_counter`, several repetitions, caches cleared (`GLNRep._cache.clear()` etc.) between cold runs.
