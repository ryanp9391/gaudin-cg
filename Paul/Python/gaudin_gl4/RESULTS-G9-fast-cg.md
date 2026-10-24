# RESULTS G9: the fast float64 gl(4) Clebsch–Gordan engine on the Ubuntu laptop (2026-10-24)

Brief: `SPEC-G9-fast-cg-ubuntu.md` (+ the answers in `questions-for-fable.md`). Machine: i5-1345U (12 threads, cpufreq governor `powersave`), 15 GB, Ubuntu 22.04, Python 3.10.12, numpy 2.2.6
(bundled OpenBLAS), scipy 1.15.3, sympy 1.14.0, mpmath 1.3.0 (`g9/environment_g9.txt`). HEAD = 63f7892 (round G8 + the answers). Nothing committed; no Wolfram; no Windows-verified
module, test, RESULTS or SPEC file was edited (no Python-3.10 shim was needed). Files: `cg_fast_gl4.py` (module + CLI + all acceptance stages), `test_cg_fast_gl4.py` (11 tests, 3.2 s),
`cg_fast_gl4_output.txt` (T1, T2, T4, T5, T6 printed output of `cg_fast_gl4.py validate; t2; t4s; t6`), `cg_fast_gl4_timings.txt` (T3), raw files in
`Paul/Cpp/AlexCGC/results/su4-same-machine/g9/` (`bench_main_stdout.txt`, `bench_lean_stdout.txt`, `environment_g9.txt`, Munich's full menu-5 tables `full_*.txt` of T6, and two
`superseded_*` logs of aborted benchmark runs, see §6). Generated, gitignored caches: `cg_fast_gl4_refcache_*.pkl` (the 50-digit reference slices of G3's `copies`, so that T1 re-runs in
2 min instead of 10) and `closed_forms_gl4_qcache.pkl` (G8d's symbolic Q records; its one-off build took ≈ 5 min on this machine).

## 0. Headline

**What the engine is.** `cg_fast_gl4.hw_vectors(lam, lamp)` returns, for every irrep μ and every copy, the canonical triple (v₃, v₄ᵃ, v₄ᵇ), the slice (G3's normalisation) and the HW vector
(plain GT basis, or orthonormal, or in Munich's rows via `se.to_munich`), in float64, with no Fractions in the production path: the three operators are applied as float reduced words
to the w₁ unit vectors of the E₁₂-cut space W₁ (Route B) — or evaluated from the G8d closed forms (Route A, V₃ and V₄ᵃ only) —, the physical subspace is the null space of the E₂₃^{k₂}, E₃₄^{k₃} powers
on W₁ and the triple is diagonalised on it (Rayleigh–Ritz), then G1's float level recursion gives the full vector. **Fallbacks to the exact path: 0 in 1796 sectors (10 504 copies) of all runs of T1/T5, 0 in
the benchmark runs.** Tests: 11 pass, 3.2 s.

**Timings, one machine, cold, min of 5 (2 for [15,10,5,0]²; G3 columns: min of 3 for small products, one run otherwise; Munich = G7's same-machine numbers, not re-run).** t_slices = t_ops + t_eig
(the n_W-free time); G3 "n_W-free" = t_prv + t_ops (exact Fractions) + t_diag; Munich = hw-solve sum over irreps (`orig` = `cgc_hw_bench`, `eig` = repaired profiler). Full table in §3.

| product | ours t_slices | G3 exact n_W-free (speed-up) | G3 float n_W-free (speed-up) | ours t_total | G3 exact t_total (speed-up) | G3 float t_total (speed-up) | Munich orig / eig |
|---|---|---|---|---|---|---|---|
| [8,8,0,0]² (S = 8, m = 1) | 0.0542 s | 0.0742 s (x1.4) | 0.0725 s (x1.3) | 0.1443 s | 0.1623 s (x1.1) | 0.1606 s (x1.1) | 84.36 s (**x585** of ours) / 3.518 s (x24) |
| [12,12,0,0]² (S = 12, m = 1) | 0.1326 s | 0.2546 s (x1.9) | 0.2511 s (x1.9) | 0.5881 s | 0.7138 s (x1.2) | 0.7103 s (x1.2) | not run (S = 11: eig 201.8 s = x466 of our S = 11 t_total 0.4329 s) |
| [6,4,2,0]² (297 copies) | 0.0972 s | 1.727 s (**x17.8**) | 0.9269 s (x9.5) | 0.1706 s | 1.798 s (**x10.5**) | 0.9982 s (x5.9) | orig not run / eig 39.39 s (**x231**) |
| [9,6,3,0]² (1592 copies) | 0.3574 s | 26.33 s (**x73.7**) | 8.584 s (x24.0) | 0.9684 s | 26.94 s (**x27.8**) | 9.194 s (x9.5) | eig ≥ 385.1 s on 166 of 213 irreps (≥ **x398**) |
| [12,8,4,0]² (5945 copies, m ≤ 51) | 1.185 s | 276.1 s (**x233**) | 50.81 s (x42.9) | 6.284 s | 282.5 s (**x44.9**) | 57.21 s (x9.1) | not run |
| [15,10,5,0]² (17 550 copies, 2.49·10⁹ coefficients) | 3.986 s | did not finish in 600 s | 223.7 s (x56.1) | 58.18 s | did not finish | 285.3 s (x4.9) | not run |
| [16,8,0,0]² (1905 copies) | 1.167 s | 14.12 s (x12.1) | 11.90 s (x10.2) | 5.699 s | 19.59 s (x3.4) | 17.38 s (x3.0) | not run |

Reading: the gain is in the n_W-free stage and it is large exactly where there are many copies of large multiplicity (x18 at [6,4,2,0]², x74 at [9,6,3,0]², x233 at [12,8,4,0]² against G3's exact
50-digit pipeline; x10–56 against G3's float64 path, which still applies the operators in exact Fractions). On the multiplicity-free rectangular family G3 has nothing expensive to do
(w = 1 spaces), so the new stage gains only x1.0–2.9 in t_slices and x1.1–1.3 in t_total: there t_total is setup (exact `GLNRep`) plus the output-bound level recursion
(t_hw / t_total = 0.12 at S = 4 … 0.62 at S = 20; 0.64 at [12,8,4,0]², 0.87 at [15,10,5,0]²). **Lean mode** (multiplicity-one sectors skip the operators) takes t_total of S = 8 to 0.0903 s, of S = 12 to 0.4512 s (§3.2).

**T1 agreement** (float64 vs G3's 50-digit copies, 18 products in §2): slices max relative deviation **≤ 1.2e-14 on [2,1,0,0]²…[4,2,2,0]² (except 4.1e-13 for [4,2,1,0]x[3,2,1,0]), 3.9e-14 [4,3,1,0]², 5.8e-14 [6,4,2,0]², 6.0e-14 [5,3,2,1]², 5.1e-12 [9,6,3,0]²**; triples ≤ 5.8e-14 (relative to max(1,|v|));
HW residual ≤ 3.1e-14 (≤ 3.0e-11 at [16,8,0,0]²); G⊗G orthogonality of copies ≤ 4.2e-13; #physical = m in every sector (1796 of 1796); Route B vs G3's exact `weight_space_mats` (conjugated by √g) ≤ 6.9e-14 relative to max(|M̃|, 1)
on the 13 products it was run on (the exact matrices leave W₁ invariant exactly, leak 0); Route A vs Route B ≤ 6.9e-14 on the same 13 products. **The brief's 1e-12 on the slices is met for all products except the largest sectors: 5.1e-12 ([9,6,3,0]², w₁ = 26, m = 18),
2.0e-11 and 1.9e-11 (the two largest-m sectors of [12,8,4,0]²), 4.8e-11 ([16,8,0,0]², λ₃ = λ₄ = 0)**; the floor set by rounding the *exact* matrices to double is 0.6–1.3e-12 on the three worst [9,6,3,0]² sectors (§4).
**T2** (Munich literal comparison, new engine as the source, G7's cached raw outputs; nothing re-run): 42 of 42 multiplicity-free irreps MATCH, max deviation after one global sign **1.4e-15** (tol 1e-9); 94 of 94 multiplicity
sectors of G7's 8 products OK, worst max|1 − σ| **2.4e-15** (tol 1e-9; sum m = 148, max m = 5); the 6 big [6,4,2,0]² sectors run through the repaired profiler (m = 9…12) OK with max|1 − σ| ≤ **1.6e-12** (Munich prints 12 decimals). **String-test gap**
(row-normalised statistic, §1): physical ≤ 4.5e-12, spurious ≥ 0.20 at [9,6,3,0]² (gap factor 4.5e10), ≥ 0.07 at [12,8,4,0]²/[16,8,0,0]² (gap ≥ 1.4e11), 5.6e11–4e15 on the small products; the brief's criterion alone
(count of eigenvectors with s < 1e-6 equals m) holds in 1796 of 1796 sectors. T6 (descendants, optional) done: full CG tables vs Munich's menu-5 tables ≤ 1.95e-13 ([3,2,1,0]², 4096 vectors), unitarity ≤ 3.5e-13.

## 1. The algorithm as implemented (deviations from the brief §1, with reasons)

S0 (`Setup`): `GLNRep.of` ×2, `se.Recursion.of` (G1's left inverses; counted in t_setup), the float orthonormal generators of site 1: adjacent and diagonal ones from `float_ops()`, **non-adjacent ones by float commutators
(Eh_ij = [Eh_{i,j-1}, Eh_{j-1,j}], Eh_ji = Eh_ijᵀ)** instead of `float_ops()`'s lazy exact-Fraction commutators (deviation; validated: every operator matrix agrees with G3's exact `weight_space_mats` to ≤ 6.9e-14 relative, T1). The reduced words
(`so.reduced_op(name, lamp)` tries, G3's) are evaluated at λ′ on first use, so the reduction is inside t_ops exactly as in G3's accounting (G3: inside its t_ops).

S1 Route B (`operators_B`): W₁ = {m₁₂ ≤ μ₁ − λ′₂} ⊂ W = V_λ[μ − λ′] (checked that it is exactly ker E₁₂^{k₁} with `check_cut` in `test_03` on every third sector of [4,2,1,0]² and of both orders of [3,2,1,0]×[4,2,1,0]; for a void cut W₁ = W; moreover the exact matrices of G3 have zero entries from W₁ to W∖W₁ on every sector of the 13 products with an exact twin). The trie of each reduced word is evaluated on the w₁ unit vectors as
**dense arrays per weight space** ({weight: (states × w₁)}); every generator application is a cached dense weight block (cut out of the generator's COO form once per generator) — deviation from "sparse @ dense on the whole dim × w₁", chosen because
building blocks by CSR fancy indexing cost more than the words themselves at [9,6,3,0]². Shared leaf images between the three operators. Invariance of W₁ measured per sector (`operators_B(check=True)`): the norm of the image outside W₁ over max(‖M‖, 1) is ≤ 5.9e-15 on all T1 products and 1.0e-11 at [16,8,0,0]² (entries of order 10⁵); for the exact matrices it is exactly 0.
Route A (`operators_A`): float re-implementation of `closed_form_mats` for V₃, V₄ᵃ restricted to W₁ and conjugated by √g_r/√g_c: `a`, `a₄ₐ`, the `gl3_entries` chain, the b-tables and the 12 factorising hopping classes (`universal_F`) in float; **the 6 Δm₁₂ = 0 (cubic Q) classes still go
through G4's exact evaluator (Fractions) entry by entry** (counted: 5924 entries at [9,6,3,0]²) — so Route A is *not* the vectorised n_W-free engine the brief hoped for, and it is slower than Route B at every product I timed (§3). It has no V₄ᵇ (no closed form, G8d); `physical_slices(route='A')` takes V₄ᵇ from Route B.

S2 (`_slices_float`) — **the main deviation**: the brief's selection (eigenvectors of G with s(u) = max_i‖E_{i,i+1}^{k_i}u‖/‖u‖ < τ) is implemented and reported (statistics, count check) but is not the production selector and its statistic is row-normalised:
 (a) the raw statistic is useless at λ = (16,8,0,0): the entries of E^{k} span ≥ 7 orders of magnitude (k₂ = 9). First symptom: with a null-space acceptance test on the *raw* stacked powers (first non-null singular value / σ_max = 5e-7 < my threshold 1e-6) 127 of 693 sectors fell back to the exact path (`superseded_bench_run2_killed_before_nullspace_fix.txt`); after row scaling the null space was accepted, but 26 sectors were still miscounted by the eigenvector test s < τ (spurious eigenvectors with raw s ≈ 1e-6 against 6e-7 of eigenvector error on the physical ones). The fix: scale every row of every string-power matrix to unit norm (does not change the null space) — then all 693 sectors pass with a gap ≥ 1e9;
 (b) the physical *subspace* Q is taken as the null space of the stacked row-normalised powers (SVD of a (rows × w₁) matrix; required: exactly m singular values < 1e-10 σ_max); eigenvectors of G are classified physical iff their distance to Q is < 1e-3 (physical ones: ~1e-10, spurious: O(1); 0 eigenvectors in the intermediate range 1e-8…0.1 in all runs);
 (c) the triple is then diagonalised **inside Q** (Rayleigh–Ritz, T_k = QᵀA_kQ, m × m): among 11 fixed combinations of (V₃, V₄ᵃ, V₄ᵇ) the one with the best (smallest eigen-gap)/‖T_c‖ is diagonalised, the joint diagonality of the three T_k in its eigenbasis is the acceptance test (≤ 1e-7 relative). Reason: eigenvectors of G alone have error ~ε‖G‖/gap_G(w₁ × w₁ spectrum including the spurious part) = 2e-10 at (8,6,6,4) of [6,4,2,0]² (G-gap 0.037) and 1.7e-10 at [9,6,3,0]²;
 after (b)+(c): 1.3e-12 / 5.1e-12. Two numerical traps found on the way: V₄ᵃ ≡ 0 on some sectors (a "gap" of rounding noise must not select the combination: floor 1e-6·max‖T‖) and T_k that come out numerically triangular when Q = qr(Up) (LAPACK's balancing then loses 1e-8: fixed by a fixed generic rotation of Q).
 The canonical order (v₃, v₄ᵃ, v₄ᵇ) ascending with tolerance 1e-8·max(1,|v|); the slice = plain components u_c = ũ_c/√g_c, first nonzero component (> 1e-9 of the max, ascending GT index, structurally-zero components excluded) = 1. G3's `joint_diag_float` is not used.
 τ = 1e-6 (not 1e-8) on the row-normalised statistic for the cross-check count only. Fallback: wrong null-space dimension, wrong eigenvector count, complex physical eigenvector, or failed joint diagonality → G3's `copies` (50 digits) for that sector, logged in `sector['fallback']`.
S2′ (high-precision refinement by G8a's recursion at 30 digits): **not implemented** (`refine=...` raises `NotImplementedError`).
S3: unchanged (`se.hw_from_slice`, all copies of a sector at once). `hw_vectors(..., keep_hw=False)` streams (the first benchmark run was OOM-killed at [15,10,5,0]² with all 2.49·10⁹ coefficients held, 11 GB).
S4 (done): `descendants` / `cg_table` (§5). Extra: **lean mode** (`lean=True`): a multiplicity-one sector skips S1 entirely (slice = null vector of the string powers, or the unit vector if w₁ = 1; the triple of such a sector is only a label and is returned as nan).
CLI: `python -X utf8 cg_fast_gl4.py run l1 l2 l3 l4 [l1p l2p l3p l4p] [--mus a,b,c,d;e,f,g,h] [--out file] [--orthonormal] [--descendants]`, stages `validate`, `t1`, `t2`, `t4s`, `t5`, `t6`, `bench`, `benchlean`.

## 2. T1 and T2: agreement with the existing engine and with Munich

T1 per product (all sectors; reference = G3 `copies` at 50 digits reduced to float; "dev" = max over copies of max|slice − ref| / max|ref|; triples relative to max(1,|v|); HW residual = `se.hw_residual` (simple raising operators, relative); orth = max off-diagonal of the Gram matrix of the unit-norm copies of a sector (orthonormal basis);
s_phys / s_spur = max over physical / min over spurious eigenvectors of the row-normalised string statistic):

| product | irreps / copies | max w₁ | fallbacks | max dev slices | max dev triples | HW resid | G⊗G orth | s_phys max | s_spur min | route B vs exact twin | route A vs B |
|---|---|---|---|---|---|---|---|---|---|---|---|
| [2,1,0,0]² | 7 / 8 | 2 | 0 | 2.1e-16 | 7.9e-16 | 4.1e-16 | 2.2e-16 | 2.2e-16 | 9.45e-01 | 4.0e-15 | 4.0e-15 |
| [2,1,1,0]² | 6 / 7 | 3 | 0 | 6.1e-16 | 1.8e-15 | 1.8e-15 | 2.2e-16 | 8.3e-16 | 8.94e-01 | 1.8e-15 | 3.7e-17 |
| [3,2,1,0]² | 16 / 30 | 4 | 0 | 2.1e-15 | 1.1e-15 | 2.4e-15 | 9.1e-16 | 6.1e-15 | 8.61e-01 | 2.2e-16 | 3.5e-16 |
| [4,2,0,0]² | 24 / 32 | 4 | 0 | 1.2e-14 | 5.8e-14 | 1.3e-15 | 1.6e-15 | 2.3e-15 | 6.83e-01 | 6.9e-14 | 6.9e-14 |
| [4,2,1,0]² | 27 / 56 | 6 | 0 | 6.2e-15 | 1.4e-15 | 2.4e-15 | 2.2e-15 | 8.9e-15 | 7.13e-01 | 6.8e-16 | 4.2e-16 |
| [3,2,1,0]×[4,2,1,0] | 20 / 38 | 4 | 0 | 7.3e-15 | 9.6e-16 | 1.6e-15 | 8.5e-16 | 5.9e-15 | 8.95e-01 | 3.9e-16 | 3.2e-16 |
| [4,2,1,0]×[3,2,1,0] | 20 / 38 | 6 | 0 | 4.1e-13 | 2.0e-15 | 2.7e-15 | 3.4e-15 | 2.2e-14 | 6.81e-01 | 4.6e-16 | 4.6e-16 |
| [4,2,2,0]² | 19 / 26 | 6 | 0 | 1.9e-15 | 1.1e-15 | 8.3e-16 | 9.2e-16 | 5.5e-15 | 6.74e-01 | 5.6e-16 | 2.2e-16 |
| [4,3,1,0]² | 31 / 69 | 6 | 0 | 3.9e-14 | 1.2e-15 | 1.2e-15 | 4.1e-15 | 1.4e-13 | 8.45e-01 | 7.3e-16 | 5.5e-16 |
| [6,4,2,0]² | 76 / 297 | 14 | 0 | 5.8e-14 | 1.8e-14 | 6.0e-15 | 1.6e-14 | 1.1e-13 | 5.01e-01 | 6.8e-16 | 1.4e-15 |
| [5,3,1,0]² | 49 / 136 | 9 | 0 | 3.7e-14 | 2.4e-15 | 3.1e-15 | 6.5e-15 | 7.7e-13 | 5.87e-01 | 9.3e-16 | 1.1e-15 |
| [5,3,2,1]² | 27 / 56 | 6 | 0 | 6.0e-14 | 6.7e-16 | 2.4e-15 | 9.3e-15 | 1.3e-12 | 7.13e-01 | 4.4e-16 | 5.8e-16 |
| [9,6,3,0]² | 213 / 1592 | 31 | 0 | **5.1e-12** | 1.8e-14 | 3.1e-14 | 4.2e-13 | 4.5e-12 | 2.01e-01 | 1.3e-15 | 1.6e-15 |
| [6,4,2,0]×[5,3,1,0] | 60 / 194 | 14 | 0 | 7.5e-14 | 7.0e-15 | 4.1e-15 | 1.4e-14 | 4.0e-13 | 4.71e-01 | – | – |
| [7,4,1,0]×[4,2,2,0] | 39 / 61 | 13 | 0 | 6.6e-15 | 1.3e-15 | 1.8e-15 | 1.9e-15 | 6.9e-13 | 4.92e-01 | – | – |
| [4,3,2,1]×[3,2,1,1] (λ₄ > 0) | 10 / 14 | 4 | 0 | 6.3e-16 | 1.4e-15 | 3.9e-16 | 3.4e-16 | 1.8e-14 | 7.50e-01 | – | – |
| [16,8,0,0]² | 693 / 1905 | 33 | 0 | **4.8e-11** | 1.2e-11 | 3.0e-11 | 2.2e-13 | 4.6e-13 | 6.45e-02 | – | – |
| [12,8,4,0]² (no full exact reference) | 459 / 5945 | 60 | 0 | (2 sectors: 2.0e-11, 1.9e-11) | (1.3e-15, 1.5e-15) | 4.3e-13 (every 15th sector + the 5 largest n_W) | 5.7e-13 (same) | 4.6e-11 | 7.00e-02 | – | – |

"route B vs exact twin" = max over V₃, V₄ᵃ, V₄ᵇ of max|M̃ − M̃_ref| / max(max|M̃_ref|, 1) with M̃_ref[r,c] = float(M[r,c])√(g_r/g_c) from G3's exact `weight_space_mats` restricted to W₁; "route A vs B" the same for V₃, V₄ᵃ. The weight deviation of every HW vector from μ is 0 (all products). Relative to the brief:
slices ≤ 1e-12: **not met** for [9,6,3,0]² (5.1e-12), [12,8,4,0]² (2.0e-11), [16,8,0,0]² (4.8e-11) (met everywhere else; largest of the others 4.1e-13) — see §4 for why this is conditioning; triples ≤ 1e-11: met except [16,8,0,0]² (1.2e-11); HW residual ≤ 1e-11: met except [16,8,0,0]² (3.0e-11);
orthogonality ≤ 1e-11: met; #physical = m: met (0 mismatches); fallbacks 0.

**T2** (`cg_fast_gl4.py t2`, G7's cached raw outputs `hw_*.txt`, `bench_eig_*.txt`; nothing of Munich was re-run, no cached file was missing): 2a vector by vector (same support |C| > 1e-9, equal after one global sign to 1e-9): **42 of 42 MATCH**, max deviation 1.39e-15 (41 hwcgc.sh outputs; the 42nd, the top irrep of [8,8,0,0]², from G7's repaired-profiler output because `hwcgc.sh` timed out there: single vector, trivially OK).
2b subspace per sector (G7's `compare_sector` logic with our sector as the source; m Munich α-vectors vs our m canonical copies): **94 of 94 OK**, worst max|1 − σ_i(Q)| = 2.44e-15, worst Munich Gram deviation 5.6e-15, worst residual 9.7e-15 (tolerances 1e-9, 1e-9, 1e-7). The 6 largest [6,4,2,0]² sectors from the profiler (m = 12, 11, 10, 10, 10, 9 at μ = (9,7,5,3), (10,7,5,2), (10,7,4,3), (9,8,5,2), (10,8,4,2), (10,6,5,3)): max|1 − σ| = 1.6e-12, 1.4e-12, 1.4e-12, 1.3e-12, 1.4e-12, 1.1e-12 (profiler prints 12 decimals), all OK.
Total 42 + 100 sectors as in G7.

## 3. T3: same-machine benchmark

Discipline (G7): `LC_ALL=C OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`, one process, caches cleared before every timed run (`Setup.clear()` = G1's `_clear_caches` + normal-order memo + universal/evaluated reductions + `GLNRep` cache; Python `gc.collect()`), `time.perf_counter`. Ours: min of 5 cold runs (2 at [15,10,5,0]², as soon as one run exceeds 20 s). G3/G1 (same venv, same process, run after ours): min of 3 for small products, 1 run otherwise (limit 900 s, 600 s for the two largest);
the exact and float variants share the `prv_space` and the exact operator application of one run (only the diagonalisation differs: mpmath 50 digits vs `joint_diag_float`). Munich numbers are copied from G7 (`rect_table.md`, `products_summary.txt`; same laptop, 1 thread; orig = min of 3 at S ≤ 8 as in G7's table): not re-run.
`/proc/loadavg` (1/5/15 min) before each product is in `cg_fast_gl4_timings.txt`: 1-minute 0.50–1.04 throughout; **the 5- and 15-minute averages at the start (7.0 / 13.7) are the residue of my own earlier aborted benchmark run and an out-of-memory episode ten minutes before (load 88, `[15,10,5,0]²` first attempt), not of a concurrent job; `top` showed 88 % idle CPU**. Firefox, VS Code and the Claude client were running (as in G7); governor `powersave`. Column definitions: t_setup = `GLNRep` ×2 + float generators of the adjacent operators + Racah–Speiser + G1's left inverses; t_ops = S1 (Route B, incl. the reduction of the operators, the float non-adjacent generators and the weight blocks, all cold);
t_eig = S2 (G eig, string powers, SVD, classification, Ritz, canonical order); t_hw = S3; G3: t_prv + t_ops(exact) + t_diag = "n_W-free". Route A in brackets = warm time of `operators_A` over all sectors (45 s budget).

### 3.1 Full table (seconds; `n/f` = not finished in the limit)

| product | irreps / copies / max w₁ | t_setup | t_ops [route A] | t_eig | **t_slices** | t_hw | **t_total** | G3 exact: n_W-free / total | G3 float: n_W-free / total | speed-up of t_slices (vs exact / float) | speed-up of t_total (vs exact / float) | t_hw / t_total | Munich orig / eig [s] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| [4,4,0,0]² | 15 / 15 / 1 | 0.0062 | 0.0155 [0.02] | 0.0024 | **0.0179** | 0.0032 | **0.0273** | 0.0183 / 0.0275 | 0.0178 / 0.027 | x1.0 / x1.0 | x1.0 / x1.0 | 0.12 | 0.03126 / 0.02583 |
| [5,5,0,0]² | 21 / 21 / 1 | 0.011 | 0.0207 [0.022] | 0.0033 | **0.024** | 0.007 | **0.042** | 0.0265 / 0.0442 | 0.0257 / 0.0434 | x1.1 / x1.1 | x1.1 / x1.0 | 0.17 | 0.3533 / 0.09284 |
| [6,6,0,0]² | 28 / 28 / 1 | 0.0186 | 0.0274 [0.026] | 0.0045 | **0.0319** | 0.014 | **0.0646** | 0.0383 / 0.0703 | 0.0372 / 0.0692 | x1.2 / x1.2 | x1.1 / x1.1 | 0.22 | 2.447 / 0.2271 |
| [7,7,0,0]² | 36 / 36 / 1 | 0.0293 | 0.036 [0.032] | 0.0058 | **0.0418** | 0.0259 | **0.0972** | 0.0524 / 0.1072 | 0.0511 / 0.1058 | x1.3 / x1.2 | x1.1 / x1.1 | 0.27 | 15 / 1.018 |
| [8,8,0,0]² | 45 / 45 / 1 | 0.0442 | 0.0468 [0.035] | 0.0073 | **0.0542** | 0.0458 | **0.1443** | 0.0742 / 0.1623 | 0.0725 / 0.1606 | x1.4 / x1.3 | x1.1 / x1.1 | 0.32 | 84.36 / 3.518 |
| [9,9,0,0]² | 55 / 55 / 1 | 0.0642 | 0.0601 [0.041] | 0.009 | **0.069** | 0.074 | **0.2075** | 0.1182 / 0.2588 | 0.116 / 0.2566 | x1.7 / x1.7 | x1.2 / x1.2 | 0.36 | – / 15.47 |
| [10,10,0,0]² | 66 / 66 / 1 | 0.093 | 0.0755 [0.044] | 0.0107 | **0.0862** | 0.1169 | **0.2964** | 0.1527 / 0.3603 | 0.1501 / 0.3576 | x1.8 / x1.7 | x1.2 / x1.2 | 0.39 | – / 56.83 |
| [11,11,0,0]² | 78 / 78 / 1 | 0.1285 | 0.1135 [0.049] | 0.0127 | **0.1262** | 0.1779 | **0.4329** | 0.217 / 0.5234 | 0.2139 / 0.5204 | x1.7 / x1.7 | x1.2 / x1.2 | 0.41 | – / 201.8 |
| [12,12,0,0]² | 91 / 91 / 1 | 0.1931 | 0.1175 [0.054] | 0.0152 | **0.1326** | 0.2621 | **0.5881** | 0.2546 / 0.7138 | 0.2511 / 0.7103 | x1.9 / x1.9 | x1.2 / x1.2 | 0.45 | – / – |
| [14,14,0,0]² | 120 / 120 / 1 | 0.3383 | 0.1762 [0.071] | 0.0207 | **0.1969** | 0.5312 | **1.067** | 0.4331 / 1.302 | 0.4283 / 1.297 | x2.2 / x2.2 | x1.2 / x1.2 | 0.50 | – / – |
| [16,16,0,0]² | 153 / 153 / 1 | 0.5491 | 0.2468 [0.076] | 0.0276 | **0.2744** | 0.9939 | **1.818** | 0.7338 / 2.277 | 0.7275 / 2.271 | x2.7 / x2.7 | x1.3 / x1.2 | 0.55 | – / – |
| [20,20,0,0]² | 231 / 231 / 1 | 1.208 | 0.5193 [0.101] | 0.0449 | **0.5642** | 2.919 | **4.693** | 1.601 / 5.688 | 1.59 / 5.677 | x2.8 / x2.8 | x1.2 / x1.2 | 0.62 | – / – |
| [4,2,1,0]² | 27 / 56 / 6 | 0.0076 | 0.0212 [0.116] | 0.0105 | **0.0318** | 0.0044 | **0.0439** | 0.142 / 0.1538 | 0.1017 / 0.1136 | x4.5 / x3.2 | x3.5 / x2.6 | 0.10 | 0.291 / 0.0553 |
| [6,3,0,0]² | 58 / 90 / 7 | 0.0242 | 0.0403 [0.143] | 0.0211 | **0.0614** | 0.0185 | **0.1043** | 0.194 / 0.235 | 0.1652 / 0.2063 | x3.2 / x2.7 | x2.3 / x2.0 | 0.18 | 12.68 / 1.028 |
| [6,4,2,0]² | 76 / 297 / 14 | 0.0381 | 0.0577 [0.645] | 0.0394 | **0.0972** | 0.035 | **0.1706** | 1.727 / 1.798 | 0.9269 / 0.9982 | x17.8 / x9.5 | x10.5 / x5.9 | 0.21 | – / 39.39 |
| [8,4,0,0]² | 115 / 205 / 11 | 0.0706 | 0.0832 [0.399] | 0.0473 | **0.1306** | 0.0701 | **0.2719** | 0.6347 / 0.772 | 0.5388 / 0.6761 | x4.9 / x4.1 | x2.8 / x2.5 | 0.26 | 73.82* / 40.34 |
| [9,6,3,0]² | 213 / 1592 / 31 | 0.2487 | 0.1959 [4.09] | 0.1615 | **0.3574** | 0.3614 | **0.9684** | 26.33 / 26.94 | 8.584 / 9.194 | x73.7 / x24.0 | x27.8 / x9.5 | 0.37 | – / 385.1** |
| [12,8,4,0]² | 459 / 5945 / 60 | 1.083 | 0.6095 [17.8] | 0.575 | **1.185** | 4.012 | **6.284** | 276.1 / 282.5 | 50.81 / 57.21 | x233.1 / x42.9 | x44.9 / x9.1 | 0.64 | – / – |
| [15,10,5,0]² | 846 / 17550 / 103 | 3.476 | 1.959 [45+, 697 of 846 sectors] | 2.027 | **3.986** | 50.71 | **58.18** | n/f / n/f | 223.7 / 285.3 | n/f / x56.1 | n/f / x4.9 | 0.87 | – / – |
| [16,8,0,0]² | 693 / 1905 / 33 | 1.51 | 0.752 [7.56] | 0.4147 | **1.167** | 3.019 | **5.699** | 14.12 / 19.59 | 11.9 / 17.38 | x12.1 / x10.2 | x3.4 / x3.0 | 0.53 | – / – |

\* [8,4,0,0]² orig: 80 of 115 irreps only (memory guard; on those rows eig = 7.57 s). \*\* [9,6,3,0]² eig: 166 of 213 irreps (n_W ≤ 6000), a lower bound; [6,4,2,0]² orig not run. For the two largest products the Munich code was not run at all.
Quoted `/proc/loadavg` lines before every product: `cg_fast_gl4_timings.txt`. "[12,12,0,0]²" is the S = 12 row of the rectangular family (listed once). Munich's `eig` is marginally faster at S = 4 (0.0258 s vs our 0.0273 s) and for the smaller products of G7's tables ([2,1,1,0]²: eig 0.0003 s, [3,2,1,0]²: 0.011 s; not re-run here, ours ≥ 0.01 s), ours is x1.26 faster at [4,2,1,0]² (0.0439 vs 0.0553 s) and x2.2 at S = 5: the crossover with Munich's `eig` is at about S = 4 / [4,2,1,0]², as in G7.

Observations. (i) The structure of the cost: t_setup (exact `GLNRep` in Fractions + G1's left inverses, identical in G3) is 6–26 % of t_total; t_slices is 7–66 % (66 % at S = 4 rectangular, 57 % at [6,4,2,0]², 19 % at [12,8,4,0]², 7 % at [15,10,5,0]²); t_hw is 10–87 %. (ii) t_ops of the new engine is dominated at
small sectors by the Python overhead of the trie (≈ 0.8 ms per sector at [4,2,1,0]²: 0.0212 s / 27 sectors) and by the cold construction of the non-adjacent float generators (0.14 s of the 0.20 s at [9,6,3,0]² in a cProfile run of the same code), not by the arithmetic; the rectangular family was not profiled. (iii) G3's cost explodes with m (t_diag in mpmath: 17.8 s at [9,6,3,0]², 225 s at [12,8,4,0]²; the exact operators 8.4 s / 49.6 s); the float path removes the diagonalisation but keeps the exact operators. (iv) Route A (warm, entry by entry in Python, 71 044 Q-class entries through exact Fractions at [15,10,5,0]²)
is 3.5–29x **slower** than Route B at every multiplicity product (4.09 s vs 0.196 s at [9,6,3,0]²); it is only comparable/faster (0.02–0.10 s vs 0.02–0.52 s) on the rectangular family where it has nothing to compute.

### 3.2 Lean mode (multiplicity-one sectors skip the operators; triples of those sectors not computed), cold, min of 5

| product | lean: t_setup | t_slices | t_hw | **t_total** | full pipeline t_total | full / lean |
|---|---|---|---|---|---|---|
| [4,4,0,0]² | 0.0062 | 0.0005 | 0.0030 | **0.0098** | 0.0273 | 2.79 |
| [6,6,0,0]² | 0.0186 | 0.0009 | 0.0135 | **0.0332** | 0.0646 | 1.95 |
| [8,8,0,0]² | 0.0444 | 0.0015 | 0.0442 | **0.0903** | 0.1443 | 1.60 |
| [10,10,0,0]² | 0.0927 | 0.0024 | 0.1144 | **0.2098** | 0.2964 | 1.41 |
| [12,12,0,0]² | 0.1915 | 0.0035 | 0.2558 | **0.4512** | 0.5881 | 1.30 |
| [16,16,0,0]² | 0.5421 | 0.0075 | 0.9748 | **1.5252** | 1.818 | 1.19 |
| [20,20,0,0]² | 1.2418 | 0.0149 | 2.8645 | **4.1228** | 4.693 | 1.14 |
| [4,2,1,0]² | 0.0077 | 0.0251 | 0.0046 | **0.0375** | 0.0439 | 1.17 |
| [6,4,2,0]² | 0.0381 | 0.0873 | 0.0344 | **0.1601** | 0.1706 | 1.07 |
| [9,6,3,0]² | 0.2473 | 0.3377 | 0.3578 | **0.9440** | 0.9684 | 1.03 |
| [12,8,4,0]² | 1.0667 | 1.1542 | 3.9940 | **6.2187** | 6.284 | 1.01 |
| [16,8,0,0]² | 1.5255 | 0.9749 | 2.9828 | **5.4871** | 5.699 | 1.04 |

(all rows in `cg_fast_gl4_timings.txt`, section T3b.) In the rectangular family the n_W-free stage is then 0.5–15 ms (x36–38 below the full pipeline's t_slices); the product is setup + output bound. (G1's `t_b` of S = 8 in G7 is 0.0929 s: the lean t_total 0.0903 s is the same object; the non-lean 0.1443 s additionally computes the triples.)

## 4. T4 and T5: precision and robustness

**Deviation of the float64 slices from the 50-digit copies** (`cg_fast_gl4.py t1`; the 13 T1 products of §2, 644 sectors, all with an exact reference): by w₁: ≤ 4: 7.3e-14, ≤ 8: 4.1e-13, ≤ 16: 5.3e-13, ≤ 24: 4.7e-12, ≤ 40: 5.1e-12. By the min joint gap of the sector (G3's `copies.info['jgap']`, max-norm distance of the joint eigenvalue triples): gap in (1, 10]: 3.9e-14 (6 sectors), (10, 100]: 5.1e-12 (233), > 100: 1.8e-12 (117); no sector has a gap below 1.
The worst eight sectors are all in [9,6,3,0]² with w₁ = 20–30, m = 15–24 and "smallest relative slice component" 4e-6…8e-5 (the slice is normalised to its first nonzero component; the tiny components amplify the relative error, measured here relative to the max component).
**The floor:** feeding the *exact* PRV matrices (G3's `triple_mats`) rounded to double into the same Ritz routine gives 1.3e-12, 6.2e-13, 1.1e-12 on the three worst [9,6,3,0]² sectors ((13,10,9,4), (12,11,8,5), (13,10,8,5)), against 5.1e-12, 4.7e-12 and ≤ 2.1e-12 in the pipeline — i.e. 1–4x above the floor that double-precision input matrices already impose; the operator entries themselves are accurate to ≤ 6.9e-14 relative to max|M̃| (exact twin), ≈ 1e-16·max at [9,6,3,0]² (checked on (13,10,9,4): B vs exact 1.5e-16 / 3.5e-16 / 5.9e-16). The larger deviations at [12,8,4,0]² (2.0e-11, 1.9e-11, m = 50 / 51) and [16,8,0,0]² (4.8e-11; leak out of W₁ of the float operators 1.0e-11 relative; operator entries of order 10⁵) are the same effect with larger entries and more copies; no refinement stage (S2′) exists to go below it.
**String-test statistic** (row-normalised, orthonormal basis, 2-norm; min over spurious / max over physical per product, tolerance for "physical": distance of the G-eigenvector to the null space < 1e-3): [2,1,0,0]² 0.945 / 2.2e-16 (4.3e15); [3,2,1,0]² 0.861 / 6.1e-15 (1.4e14); [4,2,1,0]² 0.713 / 8.9e-15 (8.0e13); [4,3,1,0]² 0.845 / 1.4e-13 (5.9e12); [6,4,2,0]² 0.501 / 1.1e-13 (4.4e12); [5,3,1,0]² 0.587 / 7.7e-13 (7.6e11); [5,3,2,1]² 0.713 / 1.3e-12 (5.6e11); [9,6,3,0]² 0.201 / 4.5e-12 (4.5e10); [12,8,4,0]² 0.070 / 4.6e-11 (1.5e9); [16,8,0,0]² 0.0645 / 4.6e-13 (1.4e11) (min spurious / max physical (gap factor); full list in `cg_fast_gl4_output.txt`). **The gap is ≥ 1e6 everywhere (smallest 1.5e9 at [12,8,4,0]²), τ = 1e-6 lies inside it.** The brief's expectation "spurious ≥ 16, physical 0" was for the plain basis, max norm; in this basis/norm the spurious values are 0.06–0.95 and drift down with the size of λ (rows scaled to 1); the raw statistic without row scaling gave 1e-6 for spurious eigenvectors at [16,8,0,0]² (§1 S2).
**Large S of the level recursion (S3)** (`cg_fast_gl4.py t4s`; worst HW residual over all irreps of [S,S,0,0]², all multiplicity one; float slices → `hw_from_slice`): S = 4: 6.9e-16, S = 8: 5.4e-15, S = 12: 1.2e-14, S = 14: 5.7e-14, S = 16: 2.9e-13, S = 20: 2.9e-11 (G1 measured 4.6e-14 / 2.6e-13 / 3.9e-11 at S = 14 / 16 / 20: same). Against G1's exact recursion rounded to double (`hw_from_slice_rounded`) on the two worst irreps: residual 8.5e-16 / 9.5e-16 / 1.3e-15 at S = 14 / 16 / 20 (1.0–10 s per vector against 0.02–0.07 s), max |difference| of the two unit vectors 2.4e-15 (S = 14), 1.3e-14 (S = 16), 1.3e-12 (S = 20). Complete table in the output file.
**S2′** (recursion at 30 digits against the float eigenvector): not done.

**T5 robustness.** (a) The two (V₃,V₄ᵃ)-degenerate sectors: (9,7,5,3) of [6,4,2,0]² (m = 12, w = 15, w₁ = 14): min (V₃,V₄ᵃ) gap of the copies 1.7e-13 (V₃ alone: 0), min joint gap 38.4, min G gap 6.5, the Ritz step chose combination (0.3, 1.0, −0.7) (V₄ᵇ enters): slices vs G3 4.8e-14, triples 6.4e-16, no fallback. (18,14,10,6) of [12,8,4,0]² (this is the m = 50 sector the brief names; the largest multiplicity of the product is m = 51 at (19,14,10,5), not degenerate): w = 65, w₁ = 60, min (V₃,V₄ᵃ) gap 4.6e-13, min joint gap 77.8, min G gap 0.416, combination (0.3, 1.0, −0.7), no fallback, slices 2.0e-11 vs G3, triples 1.3e-15; (19,14,10,5): slices 1.9e-11, triples 1.5e-15, min G gap 2.4 (G3's 50-digit copies of that sector took 12.2 s in the first run of this stage; the (18,14,10,6) reference was taken from the cache and its time was not printed).
(b) λ ≠ λ′: [3,2,1,0]×[4,2,1,0] 7.3e-15, [4,2,1,0]×[3,2,1,0] 4.1e-13, [6,4,2,0]×[5,3,1,0] 7.5e-14, [7,4,1,0]×[4,2,2,0] 6.6e-15 (slices vs G3, §2); λ₄, λ′₄ > 0: [5,3,2,1]² 6.0e-14, [4,3,2,1]×[3,2,1,1] 6.3e-16; no fallback anywhere.
(c) Sector classes over the 1796 sectors of T1 + T5 (dev = max slice deviation vs G3 where a reference exists): void E₁₂ cut (W₁ = W): 1348 sectors (6212 copies), dev ≤ 4.8e-11 ([16,8,0,0]²); active cut (W₁ < W): 448 (4292), ≤ 5.1e-12; 0 non-void strings i = 2, 3: 614 sectors, 1 string: 878, 2 strings: 304 (4.7e-12); **every string active** (cut active and both i = 2, 3 non-void, |A| = 3): 64 sectors (1428 copies), ≤ 4.7e-12, triples ≤ 2.1e-15; fallbacks 0 in all classes.
(d) Complex-pair guard (eig of G on W₁ returned complex eigenvalues; physical eigenvectors must be real or the sector falls back): **0 of 1796 sectors** had a complex eigenvalue at all; so the guard was never exercised (a forced fallback is tested in `test_09`). The brief's criterion alone (count of eigenvectors with row-normalised s < 1e-6 equals m) held in 1796 of 1796 sectors; eigenvectors of G at an intermediate distance (1e-8…0.1) from the physical subspace: 0.
(e) `strings_of` void strings: at [S,S,0,0]² every sector has w₁ = m = 1 (no string active, s_phys = 0 and s_spur = ∞ in the table of the output file).

## 5. T6: descendants (S4, done)

`descendants(lam, lamp, mu, sector, k)` builds |M⟩ for every GT state M of V_μ from the HW vector Ψ_μ (orthonormal basis, unit norm) by E^{tot}_{l+1,l} = E⁽¹⁾_{l+1,l} ⊗ 1 + 1 ⊗ E⁽²⁾_{l+1,l} (orthonormal CSR on both sites) with the orthonormal lowering matrix of V_μ (`GLNRep.float_ops`): from every known state q and level l with exactly one unknown target t,
Ψ_t = (E^{tot}Ψ_q − Σ_{t′ known} L_l[t′,q] Ψ_{t′}) / L_l[t,q] (worklist; all states of V_μ reached in every case). `cg_table(lam, lamp)` stacks (irrep, copy, state) → sparse CG matrix. Result (`cg_fast_gl4.py t6`):

| product | CG matrix | time (all HW + descendants) | unitarity max(\|CCᵀ − 1\|, \|CᵀC − 1\|) | vs Munich menu 5 (every row of `clebsch`, raw outputs `g9/full_*.txt`) |
|---|---|---|---|---|
| [2,1,0,0]² | 400 × 400 | 0.16 s | 2.1e-15 | 7 irreps, 400 (copy, state) vectors, max deviation 1.4e-15 |
| [2,1,1,0]² | 225 × 225 | 0.08 s | 4.0e-15 | 6 irreps, 225 vectors, 3.6e-15 |
| [3,2,1,0]² | 4096 × 4096 | 1.33 s | 3.5e-13 | 16 irreps, 4096 vectors, 1.95e-13 |

The comparison with Munich is subspace-wise for m ≥ 2: the rotation Q between our canonical copies and Munich's α-copies is fixed on the HW rows (SVD → nearest orthogonal matrix), then every state's table must equal ours·Q: max deviation as above (a trivial sign gauge, `munich_check`, is assumed and confirmed by the agreement). Cost: about one sparse matvec pair per descendant; the whole table of [3,2,1,0]² (4096 vectors of 4096 components) in 1.3 s; it is output-bound (dim V_λ · dim V_λ′ entries over all irreps) and not part of the headline benchmark.

## 6. Surprises and corrections (to the brief and to G8)

1. **The brief's string statistic (§1 S2) fails without row scaling** at λ with λ₃ = λ₄ = 0 and large k (k₂ = 9 at [16,8,0,0]²): spurious G-eigenvectors with s ≈ 1e-6, 127 of 693 sectors routed to the exact path by the raw null-space acceptance test, then 26 more miscounted by the eigenvector test (the logged fallback count was the alarm that exposed it; that run's partial log is `g9/superseded_bench_run2_killed_before_nullspace_fix.txt`). Row normalisation of the string-power matrices (a null-space-invariant scaling) restores a gap ≥ 1e9.
2. Eigenvectors of G alone are less accurate than the brief assumes (1.7e-10 at [9,6,3,0]², 2e-10 at (8,6,6,4) of [6,4,2,0]²: the G-eigenvalue gap on the whole of W₁ is small, 0.037, even though the triple gap of the *physical* copies is ≥ 12): the subspace null space + Ritz step brings it to 5e-12 / 1e-13 (§1).
3. The brief's "[12,8,4,0]², m = 50" degenerate sector is μ = (18,14,10,6); the largest multiplicity of that product is m = 51 at (19,14,10,5) (not degenerate).
4. G8d's Route A is not n_W-free-fast in this implementation: the 6 Q classes need the exact evaluator, 5924 entries at [9,6,3,0]² and 71 044 at [15,10,5,0]²; the symbolic Q records (`closed_forms_gl4_qcache.pkl`) take ≈ 5 min to build the first time on this machine (the test `test_02` skips if the cache is absent). A vectorised float Q (polynomial of degree 3 in the 14 labels, `eval_Qexp` with a removable-pole branch) was not written.
5. In the rectangular (multiplicity-free) family the engine is **not** substantially faster than G3 (x1.1–1.3 in t_total, x1.0–2.9 in t_slices): the gain over Munich on that family (x24 over the repaired Munich at S = 8, x466 at S = 11, x585 over the original at S = 8, with our full pipeline including the triples) was already present in G1 (G7: x38 / x109 / x908 with G1's `t_b`); what is new is for multiplicity products (§0 table). Lean mode (x1.1–2.8 in t_total) is the right interface for that family.
6. The first benchmark run died with the kernel's OOM-killer at [15,10,5,0]² (all 2.49·10⁹ HW coefficients, 11 GB, were kept in a dict of sectors); `keep_hw=False` streams. The load-average residue of that episode is visible in the 5/15-minute columns of the quoted loadavg lines; the timings of the final runs are from the later quiet run (`bench_main_stdout.txt`, 1-minute load 0.5–1.0). A first, discarded benchmark run (`superseded_bench_run1_partial_OOM_at_15_10_5_0.txt`) had the same t_total to ≤ 1 % for the products it covered ([6,4,2,0]², [8,4,0,0]², [9,6,3,0]², [12,8,4,0]²).
7. `float_ops()` of `GLNRep` builds the non-adjacent generators through the lazy exact commutators; the float commutator replacement of this module agrees with the exact matrices to ≤ 1e-12 relative (T1 twin) and is cheaper (0.14 s for 6 generators at [9,6,3,0]² instead of the exact build). No change to the G1–G8 modules was necessary.
8. The exact references (`copies`, 50 digits) are slow: 27 s for [9,6,3,0]², 12 s for one m = 50 sector of [12,8,4,0]², a few minutes for the 2 000-copy [16,8,0,0]²; hence the gitignored `cg_fast_gl4_refcache_*.pkl`. Tests use 30 digits.

## 7. Not done / limits

* **S2′** (30-digit lattice-recursion refinement of the slices) — not implemented; the float slices are limited to the 1e-12…5e-11 of §4 on the largest sectors. T4's "S2′ accuracy and cost on the 20 longest-grading sectors" therefore has no data.
* T1 slices ≤ 1e-12 is **not** met at [9,6,3,0]² (5.1e-12), at the two largest sectors of [12,8,4,0]² (2.0e-11), at [16,8,0,0]² (4.8e-11); triples and HW residual of [16,8,0,0]² are 1.2e-11 and 3.0e-11. [12,8,4,0]² has no complete 50-digit reference (G3 exact takes 282 s plus the diagonalisation); only its two largest-m sectors were compared with the exact copies, and its HW checks cover every 15th sector and the 5 largest n_W (full vectors cost 193·10⁶ coefficients).
* Route A: V₃, V₄ᵃ only (no V₄ᵇ closed form), entry-by-entry Python, Q classes through G4's exact evaluator; not vectorised, slower than Route B; not timed at products where building `GLNRep` is the bottleneck (the [30,20,10,0]² case of the brief: dim 1.77·10⁶, not run; the exact `GLNRep` itself would need hours and GB).
* Munich: no new Munich run (none of the cached files was missing); no Munich numbers exist for [12,8,4,0]², [15,10,5,0]², [16,8,0,0]², S ≥ 12; [9,6,3,0]² eig is partial (166/213 irreps) and [8,4,0,0]² orig partial (80/115); the Munich columns are G7's (min of 3 at S ≤ 8, one run elsewhere), taken at G7's time, not re-measured in the same run as ours. G3 numbers were measured in the same process directly after ours; min of 3 only for products with t_total < 0.5 s.
* [15,10,5,0]²: G3's exact variant did not finish within the 600 s limit (its float variant, which shares the exact operator pass, finished: 285 s); Route A's timing was budget-limited to 45 s (697 of 846 sectors). Ours: min of 2 runs.
* No multithreaded run, no C++; the cpufreq governor was `powersave` (as in G7); Firefox/VS Code/Claude were running.
* S4 validated on [2,1,0,0]², [2,1,1,0]², [3,2,1,0]² only; no timing of descendants at larger products.
* `to_munich` uses `munich_check`'s sign gauge (trivial) — not re-verified for N ≠ 4 or beyond the products named.
* Every test of T1 ran on the *cached 50-digit references* produced by the same session (`copies(dps=50)`); the cache is keyed by the sector only and is not invalidated by changes in G3 (gitignored; delete it to regenerate).

## 8. Usage

```python
import cg_fast_gl4 as cg
secs = cg.hw_vectors((9,6,3,0), (9,6,3,0))                    # dict mu -> Sector (float64; ~1 s)
sec = secs[(13,10,9,4)]
sec['trip']        # (m x 3) canonical triples (v3, v4a, v4b), ascending
sec['slices']      # (w x m) plain-GT slices on V_lam[mu - lamp], first nonzero component 1
sec['hw']          # se.HW: .rows (site-1 state), .cols (site-2 state), .vals (nW x m, plain GT), .vals_hat (orthonormal)
cg.hw_vectors(lam, lamp, mus=[mu1, mu2], orthonormal=True)    # only some irreps; unit-norm orthonormal-basis vectors
cg.hw_vectors(lam, lamp, keep_hw=False, lean=True)            # streaming; multiplicity-one sectors skip the operators (triple = nan)
cg.physical_slices(lam, lamp, mu)                             # S1+S2 only (the n_W-free part)
cg.operators_W1(lam, lamp, mu, route='A')                     # V3, V4a from the closed forms (orthonormal basis, on W1)
C = cg.cg_table((3,2,1,0), (3,2,1,0))                         # full CG matrix, sparse, orthonormal bases
```
CLI: `Paul/Python/.venv/bin/python -X utf8 cg_fast_gl4.py run 4 2 1 0 3 2 1 0 --mus "5,4,2,1" --out out.txt [--descendants]`; acceptance stages `validate` (T1 + T4 + T5, ≈ 2 min with the reference cache, ≈ 10 min without), `t2`, `t4s`, `t6`, `bench`, `benchlean`;
tests: `python -X utf8 test_cg_fast_gl4.py` (3 s). For timings use `LC_ALL=C OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`.
Each Sector has `fallback` (None or the reason), `s_phys`, `s_spur`, `info['ritz_*']`, `t_ops`, `t_eig`, `t_hw`.

## 9. Is a C++ port worth it?

Measured fractions of t_total (§3.1): at [6,4,2,0]² 22 % setup / 57 % slices / 21 % hw; [9,6,3,0]² 26 / 37 / 37; [12,8,4,0]² 17 / 19 / 64; [16,8,0,0]² 26 / 20 / 53; [15,10,5,0]² 6 / 7 / 87; rectangular S = 20: 26 / 12 / 62. The n_W-free stage — the one the algorithm of G8 is about — is already 18–233 x
cheaper than G3 and ≥ x398 below Munich's `eig` where Munich ran, and it is now 7–57 % of the total at the multiplicity products; a port of it (even to zero cost) gains at most 1.1–2.3 x at the products above and ≤ 1.1 x beyond [12,8,4,0]². The remaining stages are (1) the exact `GLNRep` construction (Fractions; 17–26 % of the total at the mid products, 1.1–3.5 s at [12,8,4,0]² / [15,10,5,0]²):
a direct float construction of the orthonormal GT generators (from the closed GT formulas, no exact generators and no exact metric) in C++ — or even vectorised numpy — would cut it by a large factor but the exact generators are also what G1 and the reductions use; (2) the level recursion `hw_from_slice` (87 % at [15,10,5,0]², 20 ns per output coefficient there; 1–6 µs per coefficient on the rectangular family, where Python overhead dominates): output-bound, the dense block products are already BLAS;
a C++ port would help mostly for small products (per-call overhead) and cannot beat the cost of writing 2.5·10⁹ coefficients (20 GB for [15,10,5,0]², more than the laptop's RAM; streaming is the only way). **Verdict: a port of the slice stage is not worth it (≤ 1.1–2x overall); if latency per product matters, port the float generator construction and the recursion's per-level loop (together 43–93 % of the total at the products above) — expected overall gain 2–4x at the mid products, ≈ 1x at the output-bound giants — and keep the slice stage in Python.** The estimates for the ports are extrapolations from the measured fractions, not measurements.
