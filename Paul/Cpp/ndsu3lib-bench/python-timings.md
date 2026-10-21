# Same-machine timings of the Python engines and the Munich code (2026-10-21, Ubuntu laptop)

Machine: Intel i5-1345U (12 threads), 15 GB RAM, Ubuntu (Linux 6.8). Python 3.10.12, numpy 2.2.6 (bundled scipy-openblas 0.3.29, multithreaded),
scipy 1.15.3, sympy 1.14.0, mpmath 1.3.0, from the existing venv `Paul/Python/.venv` (nothing missing; pyproject.toml/uv.lock NOT changed).
Load: every timing below was taken while the other agent was building/benchmarking ndsu3lib; `/proc/loadavg` (1 min) is quoted per batch.
The load figure includes the benchmark itself (numpy BLAS uses several threads, user time ~3.5x wall time at a = 7), so a reading of ~1.7 means
roughly one other busy core. Quoted time = minimum over 5 repeats; spread = (max - min) over the same 5.

## Commands (run from `Paul/Python/gaudin_sl3`; driver scripts are in the session scratchpad, the repo scripts are untouched)

The repo CLIs (`slice_engine.py bench [amax [amin]]`, `jacobi_cg.py bench`) loop over many a (jacobi_cg bench also runs the 30-digit mpmath
pipelines and cannot be restricted to one a), so the same functions were called directly, one a at a time:

    ../.venv/bin/python -X utf8 se.py <a> 5     # loops: slice_engine.bench_one(a, 1)  (cold caches each time; 5 separate calls)
    ../.venv/bin/python -X utf8 jc.py <a> 5     # loops: jacobi_cg.stage_bench([a], verbose=False, with_R=False)  (each call itself keeps the best of 5 inner runs)

(`se.py`/`jc.py` are 5-line wrappers that print the returned timing dict entries.) Equivalent to the CLI: `python -X utf8 slice_engine.py bench 14 7`
(min of 5 cold runs per row; our min of 5 outer repeats of reps=1 is the same statistic).

## What each timing covers

* slice_engine t_a ("n_W-free part"): build gl3 rep generators, Racah-Speiser decomposition, PRV spaces for every irrep, exact (Fraction) V_red for every
  irrep, float64 eigen-decomposition of V_red -> the m eigenvalues and the multiplicity-space (slice) eigenvectors of every copy. No n_W-sized object is touched.
* slice_engine t_b = t_a + metric/orthonormal float matrices/left inverses (Recursion.of) + reconstruction of the FULL highest-weight vector (all n_W
  components in the GT basis) of every copy (the "entries" column: 155 176 numbers at a = 7, 10 739 904 at a = 14).
* jacobi_cg t_float: closed-form Jacobi chains of V_red (no PRV space, no exact linear algebra, no dim V), float64 eigenvalues of every chain, and the forward
  three-term recurrence giving the slice (weight-space vector) of every copy. This is the slice only (the HW vector on the weight space W), not the n_W-sized
  vector. t_float = irreps + chain + eig_f + rec_f. (3rd column: "twisted" recurrence, more accurate, shown separately.)

## Table (Lambda = [2a,a,0], Lambda (x) Lambda)

| code | a | what is computed | #irreps / #copies | min time | spread (max-min) | load avg during batch |
|---|---|---|---|---|---|---|
| slice_engine t_a | 7 | HW slices of all copies, n_W-free | 113 / 344 | 0.128 s | 0.034 s (0.128-0.162) | 1.65 |
| slice_engine t_b | 7 | + every full HW vector (155 176 numbers) | 113 / 344 | 0.185 s | 0.036 s (0.185-0.221) | 1.65 |
| slice_engine t_a | 14 | HW slices, n_W-free | 421 / 2255 | 1.306 s | 0.080 s (1.306-1.386) | 1.66-1.77 |
| slice_engine t_b | 14 | + every full HW vector (10 739 904 numbers) | 421 / 2255 | 2.105 s | 0.123 s (2.105-2.228) | 1.66-1.77 |
| jacobi_cg t_float | 7 | slices of all copies, three-term recurrence, float64 | 113 / 344 | 7.7 ms | 4.7 ms (7.7-12.4; the max is the first call, warm-up) | 2.87 |
| jacobi_cg t_float | 14 | same | 421 / 2255 | 49.0 ms | 5.0 ms (49.0-54.0) | 2.66-2.87 |
| Munich `cgc_hw_bench` (original) | 7 | HW CG coeffs, sum of hw-solve over irreps, reps=1 | 113 / 344 | 205.9 s | single run, no spread | not recorded (earlier session) |
| Munich `cgc_hw_prof_ob` (HW_MODE=eig + OpenBLAS) | 7 | same | 113 / 344 | 10.93 s | single run, no spread | not recorded |
| Munich (any variant) | 8 | - | - | out of memory (original allocates n_coeff x dim^2 doubles, 4.4 GB for the a = 7 singlet) | | |

Detail of the jacobi_cg t_float minimum (ms): a = 7: irreps 0.7 + chains 4.2 + eigenvalues 2.5 + recurrences 0.3; a = 14: 2.4 + 27.0 + 15.5 + 3.6.
Twisted (more accurate) recurrence variant: 1.3 ms (a = 7), 14.7 ms (a = 14) for the recurrences alone instead of 0.3 / 3.6 ms.
The a = 14 jacobi_cg slice-only time is NOT comparable to slice_engine t_b (different object); the like-for-like pair is jacobi_cg t_float vs slice_engine t_a.

## Munich sources (not rerun)

* Original: `Paul/Cpp/AlexCGC/results/full_decomp_2aa0_a7.txt`, per-irrep `hw_solve_s` (reps=1) summed over the 113 irreps with `LC_ALL=C awk` = 205.878 s;
  also stated in HANDOFF ("Done (2026-10-21, Ubuntu)": orig 20.5 s singlet / 206 s sum). Same i5-1345U laptop, g++ 11.4, reference LAPACK/BLAS.
* Repaired: `results/full_decomp_2aa0_a7_cgc_hw_prof_ob_eig.txt` (binary `cgc_hw_prof_ob`, HW_MODE=eig, OpenBLAS) sum = 10.932 s (singlet 0.52 s; HANDOFF: 0.45 s / 10.9 s).
  Intermediate variants from the same directory (sums): eig with reference BLAS 17.67 s; compact + OpenBLAS 68.5 s.
* The load during those runs was not recorded; they were single runs (reps=1) with no repeats, so they carry unknown noise (HANDOFF/README report these as the best numbers available).
  Note the Munich totals are the HW-solve only (coefficients of the HW state, all n_W components), i.e. comparable to slice_engine t_b, not to t_a.

## Cross-machine ratios (this laptop / Windows i7-12700KF numbers quoted in RESULTS-slice-engine.md and RESULTS-jacobi-cg.md)

| quantity | Windows (quoted) | this laptop (min) | ratio |
|---|---|---|---|
| slice_engine t_a, a = 7 | 0.090-0.101 s | 0.128 s | 1.3-1.4 |
| slice_engine t_b, a = 7 | 0.146-0.17 s | 0.185 s | 1.1-1.3 |
| slice_engine t_a, a = 14 | 1.0-1.144 s (quiet batch ~1.0) | 1.306 s | 1.1-1.3 |
| slice_engine t_b, a = 14 | 2.0-2.36 s | 2.105 s | 0.9-1.05 |
| jacobi_cg t_float, a = 7 | 5.9 ms ("this machine", R run in same session) | 7.7 ms | 1.3 |
| jacobi_cg t_float, a = 14 | 39 ms | 49 ms | 1.26 |
| Munich repaired, a = 7 | 12.7 s was the Mathematica pipeline on Windows (not Munich) | 10.9 s Munich on this laptop | - |

So this laptop is about 1.0-1.4x slower than the Windows desktop for the Python engines (the quoted Windows figures themselves carry +-30% noise from other agents,
as RESULTS-slice-engine.md says, and the present batch was also run under load ~1.7-2.9), i.e. the machines are within noise of each other for the larger runs.
Conclusion for the same-machine comparison at a = 7 (all 113 irreps / 344 copies, HW CG data): Munich original 206 s, Munich repaired 10.9 s,
slice_engine 0.19 s (full vectors) / 0.13 s (slices only), jacobi_cg 8 ms (slices only).
