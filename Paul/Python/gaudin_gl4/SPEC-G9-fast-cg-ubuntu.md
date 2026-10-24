# Spec G9: the fast gl(4) Clebsch–Gordan engine on the Ubuntu laptop — algorithm, implementation, tests, same-machine benchmark (written 2026-10-24 on Windows, after round G8; to be run on the Ubuntu laptop)

For a general-purpose subagent (SONNET) in a session on the **Ubuntu laptop** (i5-1345U, 12 threads, 15 GB, g++ 11.4, reference LAPACK/BLAS 3.10; venv `Paul/Python/.venv`,
Python 3.10.12, numpy 2.2.6 with bundled OpenBLAS, scipy 1.15.3, sympy 1.14.0, mpmath 1.3.0; the Munich code `Paul/Cpp/AlexCGC` builds there (`build.sh`, the repaired
profiler `cgc_hw_prof.cpp`, `HW_MODE` = orig | compact | eig) and its raw outputs of G7 are cached in `Paul/Cpp/AlexCGC/results/su4-same-machine/`). **No Wolfram, no
notebooks, no git commits unless Paul asks.** Never edit the Windows-verified modules (`glN_reps.py`, `slice_engine_gl4.py`, `gaudin_gl4_ops.py`, `slice_ops_gl4.py`,
`oper_gl4.py`, `structure_gl4.py`, `lattice_recursion_gl4.py`, `nested_gl3_gl4.py`, `relabelled_gl4.py`, `closed_forms_gl4.py`, `munich_bridge.py`, the tests, any
RESULTS/SPEC) except for a Python-3.10 shim if unavoidable (say exactly what in RESULTS). Files you own: `Paul/Python/gaudin_gl4/cg_fast_gl4.py`, `test_cg_fast_gl4.py`,
`RESULTS-G9-fast-cg.md`, `cg_fast_gl4_output.txt`, `cg_fast_gl4_timings.txt`, new raw files under `Paul/Cpp/AlexCGC/results/su4-same-machine/g9/`, and one new block at the top
of `HANDOFF.md` when done. Prerequisite: round G8 (`SPEC/RESULTS-G8{a,b,c,d}-*`, the four `*_gl4.py` modules and their tests) must be committed on Windows and pulled here —
`git log -1` must show it; if not, stop and say so. Time-box: one long session (≈ 4–5 h of compute and writing); deliver what is finished and say precisely what is not.

Read first: this file; `Paul/Mathematica/Gaudin/gl4-jacobi-plan.md` §4 (verdicts of G8); `RESULTS-G8d-closed-forms.md` §0, §5, §8 (`closed_form_mats`, `slices_closed`,
the D4 timings and why the exact closed-form pipeline is slower than G3's `copies`); `RESULTS-G8a-lattice-recursion.md` §0, §1, §9 (`slice_from_triple`, grading rule);
`RESULTS-G3-slice-ops-gl4.md` §8–§9, §13 (canonical copies, `copies`, `joint_diag_float`, `copy_float_slices`, timings; "Not done: float fast path for the operators");
`RESULTS-G1-slice-engine-gl4.md` §4, §6, §9, §12 (`hw_from_slice`, `hw_from_slice_rounded`, `hw_residual`, the Munich bridge, timings breakdown t_gen/t_rs/t_prv/t_setup/t_recon);
`RESULTS-G4-structure-gl4.md` §2.1–2.2 ("τ physical ⇔ no string violated"); `RESULTS-G7-munich-comparison.md` and `SPEC-G7-munich-comparison-ubuntu.md` §1 (setup, timing
discipline), §2 (`munich_bridge.py`: `compare_mf`, `compare_sector`, `run_hwcgc` with `reuse=True`); `benchmark_gl4.md`; `Paul/Python/gaudin_sl3/RESULTS-jacobi-cg.md` §7
(the gl(3) matrix-free pipeline and the float64 "twisted" two-ended recurrence that was needed for accuracy).

## 0. What the G8 technology changes in the algorithm, and what it does not

The engine of G1–G3 computes, for every irrep μ of V_λ ⊗ V_λ′ and every copy, the highest-weight (HW) CG vector in the canonical Gaudin basis: (1) the PRV slice M_μ ⊂ W =
V_λ[μ−λ′] (exact nullspace of the three string powers); (2) the three missing-label operators V₃, V₄ᵃ, V₄ᵇ on M_μ by applying their 39/78/85 reduced words EXACTLY
(Fractions) to the PRV basis — **56–62 % of the n_W-free time, 23 s for [9,6,3,0]²**; (3) joint diagonalisation (mpmath 50 digits, or float64); (4) the full HW vector by
the level recursion from the slice (float64, output-bound: 0.73 s for the 8.9·10⁶ coefficients of [9,6,3,0]²). Munich computes the same HW vectors by an O(n_W³) null space.

Round G8 established (exact, 222 weight spaces, 12 products; details in the plan §4): (a) **closed-form operators** — on every path P_κ (fixed row 3 = κ) of W the κ-block of
V₃ is a(κ;λ)·V^{(3)} + b₃ with V^{(3)} the gl(3) reduced vertex operator of (κ, λ′₁₂₃, μ₁₂₃) and **a = 1 − R₁ − R₂ − R₃**, R_i = Π_{s=1..4}(y_i − Y_s)/Π_{j≠i}(y_i − y_j)(y_i − y_j + 1)
(y_i = κ_i − i, Y_s = λ_s − s; the GT norm of the neighbour κ + e_i), likewise V₄ᵃ with a₄ₐ = d + Σ_i R_i(2y_i − |κ| + 4), d = |λ| − |κ|; between paths the 12 hopping classes have
G4's product formula and the 6 Δm₁₂ = 0 classes a cubic; `closed_forms_gl4.closed_form_mats(lam, lamp, mu)` builds V₃ and V₄ᵃ on W from the labels alone (exact on all 328
weight spaces; **V₄ᵇ has no closed form yet**); (b) **M_μ is a spectral subspace of the E₁₂-cut space W₁ = {m₁₂ ≤ μ₁ − λ′₂}** (a coordinate sub-lattice of W, triple-invariant):
the joint eigenvectors of the triple on W₁ are the m physical slices plus w₁ − m spurious ones, and **an eigenvector is physical iff it is killed by E₂₃^{k₂} and E₃₄^{k₃}**
(k_i = λ′_i − λ′_{i+1} + 1; G4/G8b: the spurious ones violate a string by a factor ≥ 16 in norm, the physical ones by 0) — so the slices can be read off a w₁ × w₁ eigenproblem
without computing the PRV nullspace; (c) the **forward lattice recursion** `lattice_recursion_gl4.slice_from_triple` builds a slice from its eigenvalue triple level by level
(0.3 ms per copy at 50 digits), and the unused equations are the characteristic equations (the gl(4) continuant); (d) [V₃, V₄ᵃ] alone determine the recursion everywhere but
separate the copies in all but two known sectors ((9,7,5,3) of [6,4,2,0]², m = 12; [12,8,4,0]², m = 50): **V₄ᵇ is needed for completeness**.

What it does NOT change: the full HW vector is output-bound (n_W coefficients per copy) and stays G1's level recursion; the eigenproblem is tiny (w₁ ≤ 32 at [9,6,3,0]², ≤ 46
at [30,20,10,0]); so the whole gain is in step (2), and the real engineering question is double precision: G8d's exact closed-form pipeline is SLOWER than G3's `copies`
(2.97 vs 1.65 s on [6,4,2,0]²) because it stays in Fractions/mpmath. **The job here is the float64 engine, with certified residuals, benchmarked on one machine against
Munich and against the G1/G3 engine.**

## 1. The algorithm (implement exactly this; every stage has a validation twin)

Inputs: λ, λ′ (dominant gl(4) weights, λ₄, λ′₄ ≥ 0), optionally a list of μ. Output: for every irrep μ (Racah–Speiser multiplicity m ≥ 1) and every copy, the HW vector in the
plain GT ⊗ GT basis (rows = site-1 index, cols = site-2 index, values) in G3's canonical normalisation (copies ordered by (v₃, v₄ᵃ, v₄ᵇ) ascending, first nonzero slice
component = 1) AND, on request, G⊗G-orthonormal (CG) normalisation and Munich's row format (`slice_engine_gl4.to_munich`).

**S0 — setup per (λ, λ′)** (reuse G1): `gr.GLNRep.of(lam)`, `.of(lamp)` (patterns, weights, exact adjacent generators, metric √g, `float_ops()` = orthonormal CSR generators);
`gr.candidates`/`racah_speiser` (all μ with multiplicities); `se.Recursion.of(lam, lamp)` (left inverses of the site-2 raising blocks, once per λ′); the universal reduced
words `so.reduced_op(name, lamp).red` for the three operators (polynomial in λ′, evaluated once). Record t_setup.

**S1 — operators on W₁ in float64, two routes, both implemented:**
  * Route B (complete, simplest): for each of the w₁ basis vectors e_c of W₁ (orthonormal basis), apply every reduced word (a product of ≤ 4 generators; the dict keys are
    tuples of (i, j) pairs, coefficient Fraction → float) with the CSR generators of `float_ops()`; accumulate the image restricted to W₁ (the operators preserve W₁; assert that
    the image has no component outside W₁ beyond 1e-12 relative). Vectorise: apply each word to the (dim × w₁) matrix of basis vectors at once (sparse @ dense), or better,
    group the words by their leftmost generator. Cost ≈ Σ_words (length × nnz) — expect well under a millisecond per sector up to [9,6,3,0]². This gives V₃, V₄ᵃ, V₄ᵇ.
  * Route A (closed form, V₃ and V₄ᵃ only): a float64 re-implementation of `closed_forms_gl4.closed_form_mats` (vectorised numpy: a(κ), a₄ₐ(κ), b₃/b₄ₐ from the tables,
    the general-λ′ gl(3) chain entries, `universal_F` on the hopping classes, the cubic Q's — read the exact rational records from `closed_forms_gl4` once and evaluate in
    float); restrict to W₁. Validation twin: Route A == Route B to 1e-12 relative on every sector of the test products; both == G3's exact `weight_space_mats` (converted) on P_small.
  Use Route B for the production engine (it has V₄ᵇ); report Route A's timing separately (it is the n_W-free, representation-free variant — the one that will matter for a
  C++ port or for products where even building `GLNRep` is the bottleneck: [30,20,10,0] has dim 1.77·10⁶).

**S2 — physical spectrum and slices from a w₁ × w₁ eigenproblem:** G = V₃ + (7/19)V₄ᵃ + (11/13)V₄ᵇ (G2's generic combination) on W₁; `numpy.linalg.eig` (non-symmetric: the
reduced operators are not symmetric on the slice); for each eigenvector u: the string test s(u) = max_i ||E_{i,i+1}^{k_i} u|| / ||u|| over the non-void strings i = 2, 3
(E₁₂ is the coordinate cut — verify once per sector that W₁ is exactly its kernel); physical iff s(u) < τ with τ = 1e-8 (report the distribution of s over physical and
spurious eigenvectors: the gap must be ≥ 10⁶; adapt τ if it is not and say so). Assert #physical = m (Racah–Speiser) in every sector; on failure (degenerate G, complex pair,
gap too small) fall back to the exact path (G3 `copies`) for that sector and log it. Then the eigenvalue triples (v₃, v₄ᵃ, v₄ᵇ) = Rayleigh quotients of the physical u (or the
diagonal of U⁻¹V_kU), the canonical ordering and normalisation of G3, and the slice = u (plain basis = u_c/√g_c). Optional refinement (S2′): one step of inverse iteration on
G − v, or G8a's recursion from the triple (at 30 digits: the certified high-precision route; measure its cost and accuracy against the float eigenvector — the gl(3) lesson is
that a plain forward recurrence loses digits for long chains, the two-ended "twisted" version does not; test the same on the longest gradings here).

**S3 — full HW vectors:** `se.hw_from_slice(lam, lamp, mu, slices)` (all copies of a sector at once, float64), `hw_residual` ≤ 1e-11 relative (G1's levels: 1e-15 at S ≤ 8,
growing to 1e-11 at S = 20); for S ≥ 14 of the rectangular family compare with `hw_from_slice_rounded`. This stage is unchanged; its time is reported separately (output-bound).

**S4 — (optional, only if S0–S3 and the tests are done) descendants = full CG tables:** for a copy with HW vector Ψ_μ, the GT basis vector |M⟩ of V_μ is reached from the
pattern M⁺ (one entry of row k raised) by the total lowering operator E^{tot}_{k+1,k} = E^{(1)}_{k+1,k} ⊗ 1 + 1 ⊗ E^{(2)}_{k+1,k} (orthonormal CSR on both sites):
E^{tot}_{k+1,k}|M⁺⟩ = Σ_{M′} c_{M⁺→M′}|M′⟩ with Munich's/our orthonormal lowering coefficients c (`gr.munich_lowering_sq`, exact), so |M⟩ = (E^{tot}|M⁺⟩ − Σ_{M′≠M} c_{M⁺→M′}|M′⟩)/c_{M⁺→M}
in a triangular sweep over the patterns of V_μ in decreasing weight (choose, for each M, the raised entry with the largest coefficient; every |M′⟩ in the sum has already been
computed if the sweep is by weight). Cost = Σ_M nnz(Ψ_M) sparse matvecs — output-bound (the whole table has dim V_λ · dim V_λ′ entries over all μ and copies). Validate on
[2,1,0,0]², [2,1,1,0]², [3,2,1,0]² against Munich's full tables (menu 5 of `clebsch` prints every row; `hwcgc.sh` keeps only the HW rows — write a small parser for the full
output; conversion CGC = σ Ψ √(g_ν g_q/GNORM2) as in G1 §9; for m ≥ 2 the comparison is subspace-wise per weight as in G7's `compare_sector`) and by the unitarity of the whole
CG matrix of a product (Σ_μ,copies,M |⟨…⟩|² = 1 per basis vector of V_λ ⊗ V_λ′, ≤ 1e-12). Report its cost separately; it is not part of the headline benchmark.

Output writer: `write_ref` format of G1 (`ref/gl4-cg-*.txt`) for the slices and `to_munich` rows for the vectors on request; a `main` CLI
`python -X utf8 cg_fast_gl4.py run l1 l2 l3 l4 [l1' l2' l3' l4'] [--mus ...] [--out file] [--orthonormal] [--descendants]` and `bench`, `validate`, `munich` stages.

## 2. Tests (acceptance; every number into RESULTS)

T1 **Agreement with the existing engine** (the Windows-verified reference): slices vs `so.copies(lam, lamp, mu, dps=50)` and HW vectors vs `se.hw_from_slice` of the PRV basis
  rotated to the copies, on P_small ([2,1,0,0]², [2,1,1,0]², [3,2,1,0]², [4,2,0,0]², [4,2,1,0]², [3,2,1,0]⊗[4,2,1,0] both orders, [4,2,2,0]², [4,3,1,0]²), [6,4,2,0]² (297 copies),
  [5,3,1,0]², [5,3,2,1]², [9,6,3,0]² (1592 copies): max relative deviation of every slice component ≤ 1e-12 (report the actual max per product), the eigenvalue triples to 1e-11,
  #physical = m in every sector (0 fallbacks expected — list any), HW residual ≤ 1e-11, G⊗G-orthogonality of the copies of each sector ≤ 1e-11, and the exact twin: Route B's
  matrices vs G3's exact `weight_space_mats` on W₁ (≤ 1e-12 relative) on P_small; Route A vs Route B on everything.
T2 **Munich literal comparison** (G7's pipeline with the new engine as the source): all 42 multiplicity-free irreps of G7 (`munich_bridge.compare_mf`-style, reusing the cached
  `hw_*.txt`: MATCH ≤ 1e-9 after one global sign) and the 100 multiplicity sectors (`compare_sector`, subspace agreement max |1 − sv| ≤ 1e-9); additionally the sectors of
  [6,4,2,0]² that G7 ran through the profiler. Nothing of Munich needs to be re-run unless a cached file is missing (then `run_hwcgc` with `timeout 1800`, memory rule of G7).
T3 **Same-machine benchmark** (the G7 discipline: `LC_ALL=C`, `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`, nothing else running, `/proc/loadavg` quoted, min of 5 cold runs,
  caches cleared as G1's `bench` does, `/usr/bin/time -v` for Munich): per product the columns t_setup / t_ops (S1, Route B; Route A in brackets) / t_eig+select (S2) /
  t_slices = sum of the three (the n_W-free time) / t_hw (S3) / t_total, against (i) G3 `copies` exact 50-digit and float64 (`joint_diag_float`) + G1 `hw_from_slice` — the
  current engine, same venv, same run — and (ii) Munich: the G7 same-machine numbers (orig, eig) where they exist, re-run only what is cheap and missing. Products:
  rectangular [S,S,0,0]² S = 4..12 (+ 14, 16, 20 for ours; Munich orig S ≤ 8, eig S ≤ 11 from G7), multiplicity products [4,2,1,0]², [6,3,0,0]², [6,4,2,0]², [8,4,0,0]²,
  [9,6,3,0]², [12,8,4,0]² (m ≤ 50; G3: 258 s exact / 37 s float), and the largest you can run in < 10 min each ([15,10,5,0]², [16,8,0,0]², [12,12,0,0]²). One table in
  `benchmark_gl4.md` style (do not edit that file; write yours into RESULTS §3 and `cg_fast_gl4_timings.txt`). State the speed-up of t_slices over G3's n_W-free part and of
  t_total over Munich, per product; and the fraction t_hw / t_total (how output-bound the engine is).
T4 **Precision:** float64 vs the 50-digit copies (T1's numbers by sector size: max deviation as a function of w₁ and of the min joint gap; `so.min_joint_gap`); the string-test
  statistic s(u) distributions (physical vs spurious, per product: min spurious / max physical); the large-S behaviour of S3 (residual vs S for [S,S,0,0]², S ≤ 20, float vs
  rounded); S2′ if implemented: recursion-at-30-digits vs float eigenvector accuracy and cost on the 20 longest-grading sectors.
T5 **Robustness:** the two (V₃,V₄ᵃ)-degenerate sectors ((9,7,5,3) of [6,4,2,0]², and [12,8,4,0]² m = 50): separation by G, min joint gap, deviation; λ ≠ λ′ products incl.
  [6,4,2,0]⊗[5,3,1,0] and [7,4,1,0]⊗[4,2,2,0]; λ₄, λ′₄ > 0 ([5,3,2,1]², [4,3,2,1]⊗[3,2,1,1]); a sector with a void E₁₂ cut (W₁ = W, 171 of 222 spaces) and one with every
  string active (|A| = 3); the complex-pair guard of `eig` (count how often G has a complex pair numerically and how it was resolved).
T6 **Descendants** (only if S4 was done): the three small products vs Munich full tables, unitarity, cost.
`test_cg_fast_gl4.py`: ≥ 8 tests (< 90 s): one sector each of T1 (slices, HW residual, orthogonality, multiplicity), Route A == Route B, the string-test gap, a Munich
multiplicity-free MATCH from a cached file, a λ ≠ λ′ and a λ₄ > 0 sector, the degenerate sector of [6,4,2,0]².

## 3. Deliverables

`cg_fast_gl4.py` (module + CLI; public `hw_vectors(lam, lamp, mus=None, orthonormal=False, route='B', refine=None)` returning per sector the copies (triples, slices, HW
vectors as (rows, cols, vals) + diagnostics: s-statistics, residuals, fallbacks), `operators_W1(lam, lamp, mu, route)`, `physical_slices(...)`, `descendants(...)` if done);
`test_cg_fast_gl4.py`; `RESULTS-G9-fast-cg.md` with §0 headline (the benchmark table in short: t_slices and t_total vs G3 and Munich on the rectangular family at S = 8, 12
and on [6,4,2,0]², [9,6,3,0]², [12,8,4,0]²; the T1/T2 agreement numbers; the string-test gap; the fallback count), §1 algorithm as implemented (deviations from §1 above, with
reasons), §2 T1–T2, §3 T3 table + breakdown, §4 T4–T5, §5 T6, §6 surprises and corrections (to this brief or to G8), §7 not done / limits, §8 usage, §9 one paragraph on
whether a C++ port is worth it (which stage would gain, by how much, given the measured fractions); raw timing logs in `cg_fast_gl4_timings.txt`; a HANDOFF block at the top.
Every statement with its tolerance; "exact" only for Fractions.
