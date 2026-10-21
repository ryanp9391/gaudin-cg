# RESULTS-slice-engine (agent R, round 6, 2026-10-21): the PRV/Zhelobenko slice engine — exact slice spaces, V_red, full HW vectors, benchmark, reference files

Files (all new; nothing existing was edited): `SPEC-slice-engine.md` (the brief), `slice_engine.py` (all code; CLI `python -X utf8 slice_engine.py table1|hw|neg|w4|bench [amax [amin]]|ref [l1 l2|all]|verify|brute|tridiag|closed|rounded|exact <a>|eigcheck|all`),
`slice_engine_output.txt` (the printed output of every run quoted below: part 1 = one run of `all`, part 2 = the supplementary single-command runs incl. the quiet benchmark batch, `exact 14`, `closed`, `rounded`), `../../Mathematica/Gaudin/ref/cg-slices-{420,630,840,1470}.txt` (the reference files, 19 / 69 / 243 / 154 kB).
Python 3.13.13, numpy 2.5.2, scipy 1.18.1, sympy 1.14, mpmath 1.3 (no gmpy2); i7-12700KF, 32 GB, Windows 11; imports `gl3_reps` (`gt_patterns`, `_raise`, `_lower`, `_weight`), `nu_general.load_w4/w4_to_mine`, `gaudin_chain.Chain` (independent two-site check) read-only.

## 0. Conventions and what was assumed

* **One site**: GT states in the `GL3Rep(Λ).patterns` order (ascending lexicographic in (m₁,m₂,k); **HW state = LAST index**); plain rational non-unitary GT matrices; E₁₃ = [E₁₂,E₂₃], E₃₁ = [E₃₂,E₂₁]. Exact sparse matrices (Fractions) are built from `gl3_reps._raise/_lower` (out-of-interlacing targets dropped like `GL3Rep._adjacent_up/_down`) and agree with `GL3Rep.generator` entry for entry on [2,1,0], [3,1,0], [3,2,0], [4,2,0].
* **Metric**: diagonal G, E_ijᵀG = GE_ji, G[HW] = 1, computed exactly by a walk from the HW state; **equal to W4's `GMET1` (after the order map) for [1,1,0], [2,1,0], [2,2,0], [3,1,0], [3,2,0], [4,1,0]**, and `check_metric()` verifies g_r E_ij[r,q] = g_q E_ji[q,r] for ALL matrix elements of E₁₂, E₂₃, E₁₃ (0 violations on [2,1,0], [3,1,0], [3,2,0], [4,2,0]; for [4,2,0] … [14,7,0] the metric written into the reference files is re-checked as E_ijᵀG = GE_ji with `GL3Rep` float matrices, 7e-17 … 3.6e-16). The floating-point recursion runs in the orthonormal basis ẽ_q = e_q/√g_q (Ê[r,q] = E[r,q]√(g_r/g_q)); output converted back to the plain GT basis. g reaches 1.7e37 at [14,7,0], 1e103 at [28,14,0], 1e170 at [40,20,0] (still inside double range up to a ≈ 30).
* **Two sites**: Ψ = ΣΨ[ν,q] e_ν⊗e_q, site 1 = ν (slow), site 2 = q (fast), 0-based flat index ν·dim + q; slice = Ψ[·,q_top] (= W4's `cg`). **W4 order** = descending lexicographic in (k,m₁,m₂): different from the GL3Rep order, map by pattern (`w4_index_map`, asserted equal to `nu_general.w4_to_mine`).
* **V** = ½Σ_{σ,τ} sgnσ sgnτ e¹_{σ₁τ₁}e¹_{σ₂τ₂}e²_{σ₃τ₃}, the 36 terms exactly as the notebook's `vTriples`. Copies: slice scaled so that its first nonzero component (ascending GT index) is exactly 1; ordered by increasing eigenvalue.
* **Timing**: `time.perf_counter()` (Windows `time.time()` has 15 ms resolution — the first quick timings I took with it were meaningless); benchmark = min of 5 cold runs (caches cleared); a first pass with 3 repetitions showed ±30% noise on single rows (a = 11,12), hence 5.

## 1. Headline results

1. **dim M_μ (PRV nullspace) = Racah–Speiser multiplicity for EVERY irrep, and dim M_μ = 0 for every dominant candidate μ = Λ+w that is not an irrep**: 15 reps incl. [2,1,0] … [14,7,0], mismatches 0 (table §2). [14,7,0]: 113 irreps, 344 copies, max m = 8, max w = 8, max n_W = 2080, Σ m dim V_μ = 512² (all as expected).
2. **The slice determines the vector, exactly.** The level recursion in exact arithmetic reproduces W4's **45/45 exact HW kets entry for entry** (rational equality, 0 mismatches), and for every irrep of [2a,a,0] with a = 7 (344 vectors, 155 176 entries, 104 240 exact consistency conditions) and **a = 14 (2255 vectors, 10 739 904 entries, 8 765 049 exact consistency conditions, dim V⊗V = 11.4·10⁶) the exact HW check (E₁₂,E₂₃,E₁₃) leaves 0 violations** (561 s for a = 14 in Fractions). A slice outside M_μ is rejected by the exact consistency conditions (78/78 random cases).
3. **V_red needs only the slice** (the plan's "level-1 ψ_q" are not needed): V_red = ½Σ sgnσ sgnτ E_{σ₁τ₁}E_{σ₂τ₂}Φ_{σ₃τ₃}, Φ_cc = λ_c, Φ_ab = −E_ab (a<b), Φ_ab = 0 (a>b) (the raising generators annihilate Ψ, projected on q_top). The full 36-term V applied to the reconstructed full vectors equals it to ≤ 1.4e-14 (all entries; ≤ 6.8e-15 at [14,7,0]).
4. **W4 comparison, all 45 branches of the six reps (37 sectors): 0 mismatches.** The exact relation: **α₁(W4) = eigenvalue(V_red) + η₁/2, η₁ = |Λ|² + cas2(Λ) − cas2(μ)/2 = |Λ|² − Ω, Ω = Σ_ij E¹_ijE²_ji = (cas2(μ) − 2cas2(Λ))/2, cas2(μ) = Σμ_i(μ_i+4−2i)** (= W4's `cas2`, |Λ| = λ₁+λ₂+λ₃), i.e. W4's LAURMAT matrix A1 = V_red + η₁/2·1 exactly (rational matrices, in W4's HWBASIS slice basis, 37/37 sectors); eigenvalues agree to 6.4e-50, slices to 8.0e-51, full vectors to 1.0e-15 (double) / 1.0e-50 (through the exact vectors), the G-norm² of the full vector to 2.9e-49; copy 1 = smaller eigenvalue = smaller α₁ (45/45).
5. **Benchmark**: Λ = [14,7,0] (a = 7, dim V = 512, 113 irreps, 344 copies, n_W up to 2080): **n_W-free part t_a = 0.10 s, everything incl. all 344 full vectors t_b = 0.17 s** (context: Munich original 206 s, Munich fixed 10.9 s, `fastEig` 12.7 s — other machines). a = 8 (Munich out of memory): 0.24 s. **a = 14 ([28,14,0], dim V = 3375, 421 irreps, 2255 copies, n_W up to 25 425, 10.7·10⁶ vector entries): t_a = 1.04 s, t_b = 2.02 s.** Extension to a = 20 (9261, 841 irreps, 6181 copies, 1.07·10⁸ entries): t_b = 9.1 s (§5).
6. **Reference files** `cg-slices-{420,630,840,1470}.txt` written, then re-parsed from the TEXT and re-verified with code independent of the engine (§6–7), plus a Munich-like brute-force null space (n_W-sized, no PRV, no recursion) for all irreps of [4,2,0], [6,3,0], [8,4,0], [14,7,0] (§7): all agree.
7. **Structure found on the way**: V_red is **tridiagonal** in the GT basis of the weight space for every irrep tested (421/421 sectors of [28,14,0]); for w = m sectors the off-diagonal products are strictly positive (1106/1106 at [28,14,0]) ⇒ V_red is a Jacobi matrix with real **simple** spectrum there; singlet eigenvalue V(Λ⊗Λ → (2a,2a,2a)) = 2a²(2a+1) = λ₁λ₂(λ₁+1) (a = 1..5, 7: 6, 40, 126, 288, 550, 1470), top irrep V = −λ₁λ₂ (a = 7: −98) (§8).

## 2. R1: PRV dimension vs Racah–Speiser

`prv_space(lam, mu)` = exact nullspace over ℚ (RREF, Fractions) of [E₁₂^{λ₁−λ₂+1}; E₂₃^{λ₂−λ₃+1}] on V_Λ[μ−Λ] (constraint skipped when the image weight does not exist); canonical basis = one vector per free column (1 there, 0 at the other free columns), so the coordinates of any v ∈ M_μ are its free entries (`PRV.coords` checks membership exactly). `racah_speiser(lam)` = Klimyk/RS from the GT weight multiplicities (ρ = (2,1,0)).

| Λ | dim V | irreps | Σ m | max m | max w | max n_W | Σ m·dim V_μ = dim² | dominant candidates tested | mismatches |
|---|---|---|---|---|---|---|---|---|---|
| [2,1,0] | 8 | 5 | 6 | 2 | 2 | 10 | yes | 5 | 0 |
| [3,1,0] | 15 | 8 | 10 | 2 | 2 | 18 | yes | 8 | 0 |
| [3,2,0] | 15 | 8 | 10 | 2 | 2 | 18 | yes | 8 | 0 |
| [4,2,0] | 27 | 13 | 19 | 3 | 3 | 45 | yes | 13 | 0 |
| [6,3,0] | 64 | 25 | 44 | 4 | 4 | 136 | yes | 25 | 0 |
| [8,4,0] | 125 | 41 | 85 | 5 | 5 | 325 | yes | 41 | 0 |
| **[14,7,0]** | 512 | **113** | **344** | **8** | **8** | **2080** | yes | 113 | **0** |
| extra: [3,0,0] | 10 | 4 | 4 | 1 | 1 | 4 | yes | 6 (2 non-irreps → dim M = 0 ✓) | 0 |
| extra: [4,1,0] | 24 | 11 | 14 | 2 | 2 | 26 | yes | 12 (1 non-irrep) | 0 |
| extra: [5,2,0] | 42 | 18 | 28 | 3 | 3 | 72 | yes | 18 | 0 |
| extra: [5,5,0] | 21 | 6 | 6 | 1 | 1 | 6 | yes | 12 (6 non-irreps) | 0 |
| extra: [6,2,0] | 60 | 23 | 37 | 3 | 3 | 99 | yes | 24 (1) | 0 |
| extra: [7,3,0] | 90 | 32 | 60 | 4 | 4 | 200 | yes | 32 | 0 |
| extra: [9,1,0] | 99 | 26 | 34 | 2 | 2 | 66 | yes | 38 (12) | 0 |
| extra: [10,5,0] | 216 | 61 | 146 | 6 | 6 | 666 | yes | 61 | 0 |

(For [2a,a,0], a ≤ 7, every dominant Λ+w turned out to be an irrep, so the "zero multiplicity" half of the test is vacuous there — hence the extra reps. λ₃ > 0 also works: [3,2,1], [5,3,2], [4,4,2], [6,4,1] PRV = RS; [3,2,1], [5,3,2], [4,3,1] full recursion + V checks pass.) Counting laws seen in the benchmark: #irreps([2a,a,0]) = 2a²+2a+1, max w = a+1, #copies(a) − #copies(a−1) = #irreps(a), max n_W = 10, 45, 136, 325, 666, 1225, 2080, 3321, … (quartic).

## 3. R2/R3: reconstruction, HW residuals, V

`hw_from_slice(lam, mu, v)` (v of shape (w,) or (w,k) in the plain GT basis on the weight space, ascending global GT index) → `HW` object (COO vectors in the plain basis; `.dense(k)`, `.flat()`; `hw_residual`, `check_V`). **The recursion, as implemented**: site-2 weights labelled (a,b) (weight labels here, not the benchmark parameter a), wt(q) = Λ − aα₁ − bα₂, unknown block X_ab = (ψ_q)_{q∈Q(a,b)} (≤ w_max² numbers, w_max = a+1 for Λ = [2a,a,0]); the HW conditions give **A_ab X_ab = −[E₁₂^{site1}X_{a−1,b}; E₂₃^{site1}X_{a,b−1}]** with A_ab the stacked site-2 raising blocks. **A_ab has full column rank and does not depend on μ or on the slice**, so its left inverse (SVD pseudo-inverse in the orthonormal basis, smallest A_ab singular-value ratio 0.26 at a = 7, 0.19 at a = 14, 0.16 at a = 20 — the recursion is well conditioned) is computed ONCE per Λ; each irrep then costs one pass over the ≤ 3a²+3a+1 weight blocks with tiny dense matmuls, batched over the m copies. The overdetermined leftover rows are exactly the PRV constraints; `ExactRecursion` (Fractions, same recursion) computes the exact left inverse AND the exact consistency conditions T·RHS = 0 and raises if one fails.

| Λ | irreps | copies | max HW residual (E₁₂,E₂₃,E₁₃) | weight deviation | 36-term V vs ev·copy, all entries (slice only) | exact twin: vectors / exact consistency conditions / HW violations |
|---|---|---|---|---|---|---|
| [2,1,0] | 5 | 6 | 3.6e-16 | 0 | 1.3e-15 (9.4e-16) | 6 / 3 / 0 |
| [3,1,0] | 8 | 10 | 5.3e-16 | 0 | 9.4e-16 (4.3e-16) | 10 / 12 / 0 |
| [3,2,0] | 8 | 10 | 4.9e-16 | 0 | 5.2e-16 (3.2e-16) | 10 / 12 / 0 |
| [4,2,0] | 13 | 19 | 6.7e-16 | 0 | 1.4e-14 (7.7e-15) | 19 / 83 / 0 |
| [6,3,0] | 25 | 44 | 8.7e-15 | 0 | 3.7e-15 (2.1e-15) | 44 / 716 / 0 |
| [8,4,0] | 41 | 85 | 1.5e-14 | 0 | 4.0e-15 (1.7e-15) | 85 / 3632 / 0 |
| **[14,7,0]** | 113 | 344 | **2.6e-14** | 0 | 6.8e-15 (1.3e-15) | (a = 7 separately: 344 / 104 240 / 0) |

Named irreps of [14,7,0]² (copies = eigenvectors of V_red; residual of the reconstructed copies): (14,14,14) m = 1, w = 1, n_W = 2080, HW 2.3e-14, V = 1470; (15,14,13) m = 2, n_W = 2016, 2.5e-14, V = 1435.1707, 1462.8293; (28,14,0) m = 1, n_W = 1, 0, V = −98; (21,14,7) m = 8, w = 8, n_W = 540, 1.8e-14; also (17,13,12), (16,14,12), (17,14,11), (20,12,10), (18,12,12), (24,12,6), (27,14,1), (22,13,7): all ≤ 2.5e-14, V vs full V ≤ 1.8e-15 (full table in the output file). Float64 eigenvalues of V_red vs mpmath (50 digits), 316 copies in the sectors with m ≥ 2 at [14,7,0]: 2.1e-15.

**HW residual vs size (float64)**: 1.6e-14 (a = 7) → 6.5e-14 (8) → 2.2e-13 (10) → 8.4e-13 (12) → 1.9e-12 (13) → **5.3e-12 (a = 14)** → 4.8e-10 (a = 20): grows by ≈ ×2.2 per unit of a (the rounding of the float slice is amplified by the norm of the chain operators E^{a+1}, a property of double precision, not of the algorithm). **So the acceptance threshold 1e-12 holds up to a = 12 (8.4e-13) and is exceeded at a = 13, 14 (1.9e-12, 5.3e-12)**; the exact twin has residual exactly 0 at a = 14 (all 2255 PRV basis vectors, 561 s of Fractions). **Cure: `hw_from_slice_rounded(lam, mu, v)` = the exact recursion rounded to double**, HW residual 1.5e-15 … 1.6e-15 at a = 14 (the three worst irreps: (37,24,23): float recursion 1.3e-12 vs rounded-exact 1.6e-15; (39,23,22): 2.1e-12 vs 1.6e-15; (40,22,22): 2.9e-12 vs 1.5e-15; [8,4,0] singlet 1.4e-14 vs 4.2e-16), at 0.7–2 s per vector on this loaded machine (singlet (28,28,28): n_W = 25 425, 1.1–2.1 s; `rounded`). So the loss is in the float arithmetic of the 600-step recursion, not in the rounding of the input slice.

Negative tests (`neg`): random exact slices outside M_μ: 78 generated ([3,1,0] 20, [3,2,0] 19, [4,2,0] 39), 78 rejected by the exact consistency conditions; random combinations inside M_μ: 102/102 accepted.

## 4. R4: validation against W4 (`ref/untwisted-exact-{110,210,220,310,320,410}.txt`), all 45 branches

Per branch: (i) multiplicity; (ii) W4's exact HWBASIS kets: slices ∈ M_μ and the exact reconstruction from the slice reproduces the whole ket exactly; (iii) V_red eigenvalue vs α₁ − η₁/2 and the exact LAURMAT identity; (iv) slice vs `cg` (up to scale, 50 digits); (v) full vector vs `vec` up to scale (double precision recursion, and 50 digits through the exact HW vectors); (vi) W4's exact G-norm² (`gnorm2`) vs ours. Copies matched by eigenvalue (and found to be in W4's order).

| Λ | branches (sectors) | mult. mismatches | exact HW-ket mismatches | LAURMAT A1 ≠ V_red+η₁/2 | cas2 / η₁ formula mismatches | max \|ev − (α₁−η₁/2)\| | max slice dev | max vec dev (double) | max vec dev (50 digits) | max G-norm² dev | copy-order mismatches |
|---|---|---|---|---|---|---|---|---|---|---|---|
| [1,1,0] | 2 (2) | 0 | 0 | 0 | 0 / 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| [2,1,0] | 6 (5) | 0 | 0 | 0 | 0 / 0 | 2.1e-50 | 2.1e-51 | 4.4e-16 | 2.1e-51 | 2.4e-50 | 0 |
| [2,2,0] | 3 (3) | 0 | 0 | 0 | 0 / 0 | 0 | 0 | 2.2e-16 | 0 | 0 | 0 |
| [3,1,0] | 10 (8) | 0 | 0 | 0 | 0 / 0 | 2.1e-50 | 5.4e-51 | 4.4e-16 | 3.3e-51 | 6.4e-50 | 0 |
| [3,2,0] | 10 (8) | 0 | 0 | 0 | 0 / 0 | 4.3e-50 | 2.4e-51 | 1.0e-15 | 3.0e-51 | 6.4e-50 | 0 |
| [4,1,0] | 14 (11) | 0 | 0 | 0 | 0 / 0 | 6.4e-50 | 8.0e-51 | 9.5e-16 | 1.0e-50 | 2.9e-49 | 0 |
| **total** | **45 (37)** | **0** | **0 / 45 kets** | **0 / 37** | **0 / 0** | **6.4e-50** | **8.0e-51** | **1.0e-15** | **1.0e-50** | **2.9e-49** | **0** |

Pitfall found by the counters: my first guess of the shift, η₁ = 2cas2(Λ) − cas2(μ)/2, is right for λ₂ = 1 ([1,1,0], [2,1,0], [3,1,0], [4,1,0]) but wrong for [2,2,0] (3 mismatches) and [3,2,0] (10); the matrix identity and eigenvalues were right all along with W4's own η₁, and η₁ = |Λ|² − Ω reproduces W4's `eta` on all 45 branches (2λ₁λ₂ for the top irrep, etc.). The identity A1 = V_red + η₁/2 means: **the V eigenvalue of the notebook (= Q0's pencil root, = `fastEig` V-spectrum) is α₁ − η₁/2, and the Laurent datum a₁ of S3 at θ₁ is a₁ = (V + η₁/2)/θ₁₂** (needed by S/D). The slice-vs-`cg` and vector-vs-`vec` scale factors are just the sign/normalisation conventions (W4: first nonzero GT component positive in its own order; here: first nonzero = 1 in the GL3Rep order).

## 5. R5: benchmark Λ = [2a,a,0] (i7-12700KF, `perf_counter`, min of 5 cold runs per row, all times in s)

**Caveat on the noise**: other round-6 agents (K: dense 14400² linear algebra, Q0, ...) were running on the same machine, so repeated batches of this benchmark differ by up to 15–40% per row (e.g. a = 14: t_b = 2.02 s / 2.36 s in two 5-repetition batches, a = 12: 1.05 s / 1.54 s). The table is the quietest 5-repetition batch (`bench 14`, block "quiet batch" in `slice_engine_output.txt`; rows 15–20 from `bench 20 15`); the batch of the final `all` run is printed in the same file (a = 7: t_a = 0.101, t_b = 0.162; a = 14: t_a = 1.144, t_b = 2.364). Timings are therefore upper bounds on the intrinsic cost, accurate to ~±30%.

(a) = exact GT generators (`gen`, O(dim V)) + Racah–Speiser (`rs`) + PRV (`prv`) + exact V_red (`Vred`) + float64 eigendecomposition for ALL irreps — no n_W-sized object is touched. (b) = (a) + metric, orthonormal float matrices and left inverses of all A_ab (`setup`) + reconstruction of all copies' full vectors (`recon`, batched over the copies of an irrep, eigenvector slices from the float eigendecomposition). `resid` = worst HW residual over all irreps (E₁₂,E₂₃; not timed). `entries` = Σ n_W·m = total number of coefficients written.

| a | dim V | irreps | copies | max w | max n_W | t_gen | t_rs | t_prv | t_Vred | t_eig | **t_a** | t_setup | t_recon | **t_b** | entries | resid |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 8 | 5 | 6 | 2 | 10 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.001 | 0.001 | 0.000 | 0.002 | 27 | 2.9e-16 |
| 2 | 27 | 13 | 19 | 3 | 45 | 0.001 | 0.000 | 0.000 | 0.002 | 0.000 | 0.003 | 0.001 | 0.001 | 0.006 | 288 | 5.2e-16 |
| 3 | 64 | 25 | 44 | 4 | 136 | 0.002 | 0.000 | 0.000 | 0.006 | 0.000 | 0.009 | 0.002 | 0.003 | 0.013 | 1 712 | 6.3e-15 |
| 4 | 125 | 41 | 85 | 5 | 325 | 0.003 | 0.000 | 0.001 | 0.014 | 0.000 | 0.020 | 0.003 | 0.007 | 0.030 | 7 113 | 7.6e-15 |
| 5 | 216 | 61 | 146 | 6 | 666 | 0.006 | 0.001 | 0.002 | 0.027 | 0.001 | 0.036 | 0.005 | 0.015 | 0.055 | 23 247 | 8.1e-15 |
| 6 | 343 | 85 | 231 | 7 | 1 225 | 0.010 | 0.001 | 0.004 | 0.046 | 0.001 | 0.061 | 0.007 | 0.028 | 0.096 | 64 032 | 1.1e-14 |
| **7** | **512** | **113** | **344** | **8** | **2 080** | 0.020 | 0.001 | 0.008 | 0.072 | 0.002 | **0.103** | 0.017 | 0.050 | **0.170** | 155 176 | 1.6e-14 |
| 8 | 729 | 145 | 489 | 9 | 3 321 | 0.021 | 0.002 | 0.013 | 0.110 | 0.002 | 0.147 | 0.014 | 0.082 | 0.243 | 340 497 | 6.5e-14 |
| 9 | 1 000 | 181 | 670 | 10 | 5 050 | 0.031 | 0.002 | 0.022 | 0.154 | 0.002 | 0.212 | 0.018 | 0.128 | 0.358 | 690 171 | 6.5e-14 |
| 10 | 1 331 | 221 | 891 | 11 | 7 381 | 0.041 | 0.003 | 0.035 | 0.219 | 0.003 | 0.300 | 0.024 | 0.205 | 0.529 | 1 311 200 | 2.2e-13 |
| 11 | 1 728 | 265 | 1 156 | 12 | 10 440 | 0.067 | 0.004 | 0.052 | 0.301 | 0.004 | 0.427 | 0.041 | 0.369 | 0.837 | 2 360 328 | 2.9e-13 |
| 12 | 2 197 | 313 | 1 469 | 13 | 14 365 | 0.078 | 0.005 | 0.076 | 0.418 | 0.005 | 0.581 | 0.041 | 0.426 | 1.049 | 4 059 705 | 8.4e-13 |
| 13 | 2 744 | 365 | 1 834 | 14 | 19 306 | 0.111 | 0.006 | 0.110 | 0.463 | 0.006 | 0.696 | 0.048 | 0.594 | 1.339 | 6 715 519 | 1.9e-12 |
| **14** | **3 375** | **421** | **2 255** | **15** | **25 425** | 0.140 | 0.007 | 0.237 | 0.638 | 0.014 | **1.037** | 0.115 | 0.872 | **2.024** | 10 739 904 | 5.3e-12 |
| *15* | 4 096 | 481 | 2 736 | 16 | 32 896 | 0.166 | 0.008 | 0.220 | 0.775 | 0.010 | 1.179 | 0.078 | 1.183 | 2.439 | 16 676 336 | 7.5e-12 |
| *16* | 4 913 | 545 | 3 281 | 17 | 41 905 | 0.190 | 0.011 | 0.333 | 0.983 | 0.022 | 1.538 | 0.110 | 1.788 | 3.436 | 25 228 833 | 2.1e-11 |
| *17* | 5 832 | 613 | 3 894 | 18 | 52 650 | 0.273 | 0.024 | 0.441 | 1.338 | 0.014 | 2.090 | 0.112 | 2.156 | 4.357 | 37 295 163 | 3.4e-11 |
| *18* | 6 859 | 685 | 4 579 | 19 | 65 341 | 0.322 | 0.016 | 0.600 | 1.518 | 0.018 | 2.473 | 0.137 | 3.028 | 5.638 | 54 004 384 | 1.3e-10 |
| *19* | 8 000 | 761 | 5 340 | 20 | 80 200 | 0.376 | 0.028 | 0.726 | 1.808 | 0.027 | 2.965 | 0.254 | 4.242 | 7.461 | 76 758 912 | 2.0e-10 |
| *20* | 9 261 | 841 | 6 181 | 21 | 97 461 | 0.527 | 0.018 | 1.036 | 2.233 | 0.032 | 3.847 | 0.237 | 4.989 | 9.073 | 107 281 449 | 4.8e-10 |

(rows 15–20 in italics = extension beyond the requested range.) Context at a = 7 (HANDOFF.md, other machines, solving the n_W-sized null space): Munich original 206 s, Munich with the linear algebra fixed (`cgc_hw_prof` eig + OpenBLAS) 10.9 s, our sparse null space + V (`fastEig`) 12.7 s — i.e. the slice engine is ≈ 1200×, 64×, 75× faster for the same output (all 344 HW vectors in the SoV eigenbasis), with the caveat of different machines. At a = 8 the Munich code runs out of memory (> 12 GB), here 0.24 s; at a = 14 its dense n_coeff × dim² array would have 25 425 × 3375² doubles = 2.3 TB, here 2.0 s.

**Scaling** (a = 7 → 14 → 20): t_a ×10.1 then ×3.7 (≈ a^{3.3}, a^{3.7}; dominated by the exact Python-Fraction V_red, 58–70% of t_a; `t_gen` is the O(dim V) GT-matrix construction, 14% at a = 14, avoidable by generating only the needed weight blocks); t_b ×11.9, ×4.5; entries ×69, ×10 (≈ a^6); t_recon ×17, ×5.7 (≈ linear in the output: 0.32 µs/entry at a = 7, 0.08 µs/entry at a = 14, 0.05 µs/entry at a = 20 — batching efficiency improves with block size). The n_W-free part scales like the number of copies × (w ≈ a) × the number of operator applications (poly(a), tiny) and **never sees n_W**; only the optional reconstruction is output-size bound.

Float64 eigen-decomposition quality (not timed): all eigenvalues real (max imaginary part 0.0), minimum relative gap between consecutive eigenvalues of a sector 0.355 (a = 2), 0.019 (a = 7), 0.0049 (a = 14), 0.0017 (a = 20) — the spectra crowd (the gap shrinks like a^{−2…−3}), so beyond a ≈ 20 use the mpmath path (`copies`) or the exact char poly; max eigen-residual of numpy's eigenvectors 2e-12 (a = 7), 3e-11 (a = 14).

## 6. R6: reference files (`ref/cg-slices-<l1><l2>0.txt`; header documents every record; regeneration is byte-identical)

| file | Λ | sectors | copies | max HWRESID | max eigen-residual (60 digits) | min relative eigenvalue gap | VEC | bytes |
|---|---|---|---|---|---|---|---|---|
| cg-slices-420.txt | [4,2,0] | 13 | 19 | 6.7e-16 | 2.4e-59 | 0.355 | yes (19, 30 digits) | 18 938 |
| cg-slices-630.txt | [6,3,0] | 25 | 44 | 8.7e-15 | 8.4e-59 | 0.125 | yes (44) | 68 579 |
| cg-slices-840.txt | [8,4,0] | 41 | 85 | 1.5e-14 | 2.2e-58 | 0.064 | yes (85) | 243 112 |
| cg-slices-1470.txt | [14,7,0] | 113 | 344 | 2.6e-14 | 2.0e-57 | 0.019 | no | 153 582 |

Records: `REP`, `DIM`, `GTBASIS` (patterns in the GL3Rep order, W4's format; the map to the notebook order is by pattern), `WEIGHTS1`, `GMET1` (exact), `NSECTORS`, `NCOPIES`; per sector (descending lexicographic μ) `SECTOR mu m w nW dimVmu cas2mu eta1` (**w = slice length, not W4's third field, which is n_W**; the plan's `SECTOR mu mult w` plus four trailing fields), `WSPACE mu {nu..}`, `PRVBASIS mu {{nu,val},..}` × m (exact rationals; k = line order), `VRED mu {{..}}` (exact, rows, column j = V_red b_j), `CHARPOLY mu {c0..cm}` (exact characteristic polynomial, ascending; its roots = V-spectrum = pencil roots), per copy `COPY mu k ev {{nu,val},..}` (30 digits, first nonzero component 1, increasing eigenvalue), `HWRESID mu k r`, `GNORM2 mu k n2` (G-norm² of the full vector with that slice, double), and for the small reps `VEC mu k {{idx,val},..}` (full vector, 30 digits, 1-based flat index (ν−1)·DIM + q). Indices are 1-based in the GL3Rep order. Numbers are plain decimals (no exponent) except the diagnostics (Python float syntax). Example uses: `se.parse_ref(path)` returns all records parsed (no eval); `se.reference_irrep(lam, mu)` recomputes one irrep (PRV basis, V_red, copies, full float vectors).

## 7. Independent re-verification (nothing in these checks uses the engine's recursion or its exact generators)

* **Written files re-parsed from the text** (`verify`): GTBASIS = `GL3Rep` order and WEIGHTS1 correct in all four; GMET1: E_ijᵀG = GE_ji with `GL3Rep` float matrices to 7e-17 … 3.6e-16; PRVBASIS satisfies both chain conditions (`GL3Rep` matrices) to ≤ 2.8e-17 relative; WSPACE = the weight space (0 wrong); VRED = the generic V_red formula on `GL3Rep` floats to ≤ 9.3e-16 (relative to the entries); CHARPOLY = exact charpoly(VRED) in all sectors and CHARPOLY(ev)/scale ≤ 1.3e-30 at every COPY eigenvalue (limited by the 30 printed digits); for every COPY the full vector reconstructed from the FILE's slice satisfies the HW conditions (plain basis, G-norm) to ≤ 2.6e-14 and the full 36-term V Ψ = ev Ψ to 8.3e-15 / 1.9e-13 / 4.9e-13 / **2.5e-12** ([14,7,0]; the plain-basis `GL3Rep` float matrices span ≈ 40 orders of magnitude, the same check in the orthonormal basis gives 6.8e-15), slice entries = COPY to ≤ 2.8e-16.
* **VEC records vs `gaudin_chain.Chain`** (exact kron two-site operators): all 19 / 44 / 85 vectors satisfy E₁₂,E₂₃,E₁₃ totals to ≤ 7.9e-15 and V Ψ = ev Ψ (36 terms from `Chain.Es`) to ≤ 1.3e-13; their site-2 = top entries equal the COPY slice exactly (dev 0.0).
* **Brute-force null space** (`brute`; Munich-like: eigen-decomposition of AᵀA of the total E₁₂,E₂₃ on the weight-μ subspace, n_W columns, orthonormal basis): for ALL irreps of [4,2,0], [6,3,0], [8,4,0] **and [14,7,0] (113 irreps, n_W up to 2080)**: null dimension = RS multiplicity (0 exceptions), the first non-null eigenvalue is exactly 1.0 (no near-degeneracy), our copies lie in the null space (distance ≤ 8.2e-13, limited by the Gram eigen-solver), the coordinate matrix of our m copies has smallest singular value 1.000000 (span = the whole null space), and the eigenvalues of the full 36-term V restricted to the null space equal ours to 3.3e-15. (Time for scale, not a benchmark: 76–93 s for the 113 irreps of [14,7,0] with dense `eigh`, on a loaded machine.)
* **Exact twin on all PRV basis vectors**: a = 7: 344 vectors / 155 176 entries / 104 240 exact consistency conditions / 0 HW violations (4.5 s); **a = 14: 2255 vectors / 10 739 904 entries / 8 765 049 exact consistency conditions / 0 HW violations (561 s, `exact 14`)**.

## 8. Structure and closed forms found on the way

* **V_red is tridiagonal in the GT basis of the weight space V[μ−Λ]** — as the operator defined by the formula on the whole weight space (not only on M_μ): tridiagonal in 13/13, 41/41, 113/113, 421/421 sectors of [4,2,0], [8,4,0], [14,7,0], [28,14,0] and in 53/53, 17/17 of [10,3,0], [6,1,0] (diagonal for [7,7,0]). The weight space is labelled by one integer (m₁ at fixed weight), so this is a three-term recurrence. In the **w = m sectors (71/113 at [14,7,0], 239/421 at [28,14,0]) the off-diagonal products M[i,i+1]M[i+1,i] are strictly positive in all cases (161/161, 1106/1106)** ⇒ V_red is a Jacobi matrix: real SIMPLE spectrum and a continuant characteristic polynomial; in the w > m sectors the full-space matrix has zero/negative products (it is not self-adjoint off M_μ) and M_μ is an invariant subspace, there the simplicity is only numerical (gaps above). An easy source of the three-term structure for the Q-side: the V-spectrum of a w = m sector is the spectrum of a w×w Jacobi matrix with explicit rational entries (e.g. (21,14,7) of [14,7,0]: diagonal 777, 801, 849, 921, 1017, 1137, …, off-diagonals −189/2, −441/2, −120, …).
* **Singlet**: V(μ = (2a,2a,2a)) = 2a²(2a+1) = λ₁λ₂(λ₁+1) for a = 1,2,3,4,5 (exact V_red, 6, 40, 126, 288, 550) and a = 7 (1470); **top irrep** μ = (2λ₁,2λ₂,0): V = −λ₁λ₂ (α₁ = 0, η₁ = 2λ₁λ₂; a = 7: −98). Many sectors have rational V-spectra (e.g. (7,4,1) of [4,2,0]: 5/2, 35/2), others quadratic surds as in W4.
* The singular values of the stacked site-2 raising blocks A_ab, hence the conditioning of the recursion, decay only slowly (0.58 at a = 1, 0.26 at a = 7, 0.16 at a = 20 for the ratio smallest/largest).

## 9. Surprises, pitfalls

* The plan's statement "only the slice and the level-1 ψ_q are needed" for V_red can be sharpened: **only the slice** (Φ_ab = −E_ab ψ_top from the raising conditions); the full 36-term V check confirms it to machine precision on whole vectors.
* The per-μ "overdetermined sparse system at each level" collapses to one **μ-independent** left inverse per site-2 weight, computed once per Λ (≤ 0.12 s at a = 14, 0.24 s at a = 20, incl. the float matrices and the metric). The recursion needs no sparse solver at all; everything is tiny dense numpy; most of the benchmark is Python overhead and exact Fractions.
* Everything is far faster than the plan's guess (≪ 1 s at a = 7 — got 0.17 s incl. all vectors; "seconds" at a = 8 and 14 — got 0.24 s and 2.0 s). The exact rational part (PRV, V_red) is never the bottleneck: 0.1 s for 113 irreps at a = 7.
* η₁ pitfall (§4): `2cas2(Λ)` coincides with the right constant only for λ₂ = 1; correct: |Λ|² + cas2(Λ). The eigenvalue the Q-system's pencil reproduces is the 36-term V eigenvalue (`ev`), the W4/S3 Laurent datum is α₁ = ev + η₁/2.
* The first non-null eigenvalue of AᵀA is exactly 1 in every irrep tested (orthonormal GT basis; sl(2) minimal singular value) — a large, robust spectral gap for any n_W-sized null-space method; the Munich/fastEig slowness is purely the linear algebra, as HANDOFF found.
* Float64 limit: HW residual ≈ 2.2^a·1e-16·const (§3); the spec's 1e-12 holds up to a = 12. For [14,7,0] (a = 7) 2.6e-14 over all irreps.
* Shell trap (tooling): inside the Bash tool heredocs `\\n`/`\\t`/`\\"` in Python source were un-escaped (a `'\n'` literal became a real newline), which silently broke two patch scripts; edits to code containing backslashes were made with the file tools.

## 10. Not done / limits

* Absolute normalisation and sign of the copies: only the scale "first nonzero slice component = 1" and the double-precision `GNORM2` (G-norm² of the full vector with that slice; checked against W4's exact `gnorm2` to 2.9e-49 on 45 branches) are provided; the Q-side norm formula (round 7, plan §1 (g)) is not touched.
* The float reconstruction is double precision; HW residual > 1e-12 for a ≥ 13 (cure: `hw_from_slice_rounded`, exact recursion rounded to double, ~1 s per vector at a = 14, residual ~1.5e-15; a batched/vectorised version of that path was not written). Eigenvectors of V_red at 30 digits are computed with mpmath (`copies`), benchmarked only for m ≤ 8 ([14,7,0] files); the benchmark's float64 `numpy.linalg.eig` has eigen-residual up to 3e-11 at a = 14.
* The GT generators are built for the whole representation (O(dim V), 0.14 s at a = 14); a variant that builds only the weight blocks along the PRV chains would make t_a independent of dim V as well (not needed at these sizes).
* λ₃ > 0 only smoke-tested ([3,2,1], [5,3,2], [4,4,2], [6,4,1], [4,3,1]); no run with a > 20; no parallelism; Windows only (no Munich runs here — the Munich/`fastEig` numbers are Paul's from HANDOFF).
* `[8,4,0]` and `[14,7,0]` have no independent exact reference beyond the checks of §7 (W4 covers a ≤ 4 reps of other shapes); Q0's pencil roots at a = 7 (round 7: reconcile with `CHARPOLY`/`COPY` eigenvalues) are the first genuinely independent check of the V-spectrum at [14,7,0] — except that the brute-force null space (§7) already reproduces all of it.
* No git operations; no notebooks; no existing file modified (the engine reads `gl3_reps.py`, `nu_general.py`, `gaudin_chain.py` and W4's files read-only).

## 11. Usage (for D, K, S, round 7)

```python
import sys; sys.path.insert(0, 'Paul/Python/gaudin_sl3')       # gl3_reps is picked up from '..'
import slice_engine as se
r  = se.reference_irrep((14,7,0), (21,14,7))      # r['prv'] (exact PRV basis), r['Vm'] (exact 8x8 V_red), r['copies'] (mpmath: ev, coords, slice dict {global GT index: value}),
                                                  # r['hw'] (float full vectors, plain GT basis), r['eta1'] (alpha_1 = ev + eta1/2)
D  = se.parse_ref('Paul/Mathematica/Gaudin/ref/cg-slices-1470.txt')   # D['sectors'][mu] -> prv, vred, charpoly, copies [(ev, {nu(1-based): val})], hwres, gnorm2
hw = se.hw_from_slice((14,7,0), (21,14,7), [0.0]*7 + [1.0])           # slice on V[mu-Lam] (w = 8 numbers, ascending GT index)
# for a >= 13 or whenever 1e-15 is wanted: se.hw_from_slice_rounded(lam, mu, exact_slice_as_list_of_Fractions)   (exact recursion, rounded to double)
```
Everything the Q-route needs per branch is the pair (ev, slice): compare slices componentwise after fixing the scale on the first nonzero entry; indices in the files are 1-based, GL3Rep order.
