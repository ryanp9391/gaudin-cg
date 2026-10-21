# ndsu3lib on (7,7) x (7,7) -> 113 irreps (GT [14,7,0]^2), same-machine comparison (2026-10-21, Ubuntu laptop)

Machine: Intel i5-1345U (12 threads), 15 GB, Linux 6.8. Load at ndsu3lib timing: `uptime` load average 1.00 (1 min) and a single-threaded benchmark,
i.e. about one other busy core (the Python agent was working). The Python numbers (`python-timings.md`) were taken at load 1.65-2.9 with multithreaded
BLAS, so they are slightly pessimistic relative to the ndsu3lib numbers. Ratios of order 1.5x are therefore within noise; the factors below are larger.

## 1. Build (reproducible: `./build.sh`)

* Library: ndsu3lib (github.com/nd-nuclear-theory/ndsu3lib, J. Herko et al.), commit `dbd38b30c3342d3deeb1e9c2f1885590926ecc0e` (2026-01-05). Built OK; SU3lib fallback not needed.
* Toolchain (user space, no sudo): micromamba 2.9.0 -> conda-forge env in `env/`: gfortran/gcc/gxx 16.2.0, cmake 4.4.4, make 4.4.1, GSL 2.8,
  libopenblas 0.3.34, liblapack 3.11.0. SU(2) coefficients from GSL (`-DSU2COEF_LIBRARY=gsl`), double precision (`-DPRECISION=double`), no OpenMP, no SU(3)-SO(3)/MPFUN.
* Flags: `-O3` for Fortran/C/C++ (CMake Release), driver `g++ -O3 -std=c++17`, linked with the static `libndsu3lib.a`, `-lgsl -llapack -lblas`. Single thread
  (`OPENBLAS_NUM_THREADS=1`; ndsu3lib itself is serial unless built with OPENMP=ON). The library's own example programs build and run.
* `env/`, `tools/`, `src_dl/`, `build/`, `bench_driver` are gitignored (rebuild with `build.sh`).

## 2. Driver (`driver.cpp`, C++ through the library's `ndsu3lib.h`, `irreps.txt`)

`irreps.txt` = the 113 GT labels mu and multiplicities m from the `SECTOR` records of `Paul/Mathematica/Gaudin/ref/cg-slices-1470.txt` (sum m = 344).
Label map: GT mu=[m1,m2,m3] -> ndsu3lib irrep (lambda,mu)_nd = (m2-m3, m1-m2). This is the transpose of the HANDOFF's (lambda3,mu3)=(m1-m2,m2-m3): ndsu3lib's highest-weight
state has eps = 2*lambda+mu and Lambda = mu/2 (header: Lambda = (mu+p-q)/2), which equals the GT highest-weight state (Lambda = (m12-m22)/2, eps = E11+E22-2E33 = m1+m2-2*m3) only with
this assignment. (The HANDOFF's "Lambda3 = lambda3/2" is in GT/Dynkin terms. The set of irreps and m are symmetric under the swap because (7,7) is self-conjugate, so the multiplicity check is
insensitive; the span check below is what confirms the map.) Modes: `mult` (rho_max check), `time [reps]`, `dump out.txt`.

* **rho_max check: ndsu3lib `outer_multiplicity` equals m for all 113 irreps, total 344 (0 mismatches).** (`timings_raw.txt`)

What each timing produces (reduced SU(3) > SU(2)xU(1) Wigner coefficients, i.e. WITHOUT the SU(2) Clebsch-Gordan factor; one call = fixed (eps3, Lambda3), returns all
(eps1,Lambda1,eps2,Lambda2) and all rho):

* (a) "full set": for each of the 113 irreps, every state (eps3,Lambda3) of the coupled irrep = (lambda_nd+1)(mu_nd+1) calls (GT pairs (m12,m22)); 8 408 calls in total,
  **4 042 816 reduced coefficients** (rho included).
* (b) "HW row": one call per irrep at eps3 = eps3max, Lambda3 = mu_nd/2 (the coupled highest-weight state), 113 calls, **26 086 coefficients** (rho included). This contains, as a subset,
  the 378 (eps1,Lambda1) rows (site 2 at its highest weight, p2=q2=0; each carries rho_max values) that are our "slice" data; the library cannot restrict to that subset, so the whole row is computed.
* Initialisation (`initialize_ndsu3lib`, once) is not included; it is not significant here.

## 3. Timings (5 repeats, ndsu3lib; two batches of 5, `timings_raw.txt`)

| object | calls / coefficients | min | spread (5 reps, batch 1 / batch 2) |
|---|---|---|---|
| (a) full set, all 113 irreps | 8 408 / 4 042 816 | **0.1157 s** | 0.1175-0.1264 / 0.1157-0.1206 |
| (b) HW row only, all 113 irreps | 113 / 26 086 | **0.78 ms** | 0.78-0.81 ms / 0.79-0.86 ms |

(Per irrep: (b) averages 7 us; the full set averages 1.0 ms. Largest irreps dominate (a).)

## 4. Correctness cross-check (`span_check.py`, output `span_check_output.txt`, uses `dump_hwrow.txt` from `bench_driver irreps.txt dump dump_hwrow.txt`)

Procedure: for each irrep mu and each rho, take ndsu3lib's reduced coefficients with site 2 at its highest weight (p2=q2=0: eps2=21, Lambda2=7/2, M2=7/2) and every site-1 (eps1,Lambda1),
eps1 = eps3max-21; multiply by the SU(2) CG <Lambda1 M1; 7/2 7/2 | Lambda3 Lambda3> (M1 = Lambda3-7/2, sympy `clebsch_gordan`); map (eps1,2*Lambda1,M1) to the site-1 GT pattern
(2*Lambda1 = m12-m22, eps1 = 3(m12+m22)-42, M1 = k-(m12+m22)/2) in the weight space V[mu-Lambda]; divide by sqrt(G_nu) (GT metric G of `slice_engine.Rep.g`, i.e. GMET1 up to ordering) to go
from ndsu3lib's orthonormal basis to the plain GT basis. No sign or phase adjustments of any kind were fitted: all signs +1. The m vectors per irrep are then compared with the exact PRV space
M_mu = { v in V[mu-Lambda] : E12^8 v = 0, E23^8 v = 0 } of `slice_engine.prv_space` (dim M_mu = m, asserted for all 113).

* Max relative distance of any of the 344 ndsu3lib vectors from M_mu (float64): **6.7e-16**; 0 of 113 irreps above 1e-8.
* Independence: min singular value of the normalised m vectors per irrep = 0.049 (full rank m); min cos of the principal angles between span(ndsu3lib) and M_mu = 0.999999999999999. So the span
  of the m ndsu3lib vectors equals M_mu in all 113 irreps (the rho-resolution and our Gaudin eigenbasis differ by an m x m invertible rotation, as expected).
* Strength of the test: 42 of 113 irreps have dim V[mu-Lambda] > m (up to m+3), so membership in M_mu is a nontrivial condition there (for the rest it is automatic and only independence is tested).
  Negative control: flipping the sign of one random slice component destroys membership in 21 of those 42 irreps (median residual 2.7e-3; in the other 21 the flip was not detected, not investigated)
  34 states of the slice weight spaces have no ndsu3lib entry; all are forbidden by the SU(2) triangle rule, i.e. CG = 0.
* Not done: the full HW vectors (all site-2 states) were not reconstructed from ndsu3lib; the check covers the slice (site 2 at highest weight), which by the PRV/Zhelobenko theorem
  (slice -> HW vector is an isomorphism) determines the HW vectors.

## 5. Same-machine table, a = 7 ([14,7,0]^2, 113 irreps, 344 copies)

Numbers other than ndsu3lib are from `python-timings.md` (same laptop, load 1.65-2.9, numpy multithreaded) and `Paul/Cpp/AlexCGC/results` (Munich, single runs, load not recorded).

| code | object | time |
|---|---|---|
| Munich original `cgc_hw_bench` | HW CG coefficients of all 344 copies, full vectors (all n_W entries) | 205.9 s |
| Munich repaired `cgc_hw_prof_ob` (HW_MODE=eig, OpenBLAS) | same | 10.93 s |
| ndsu3lib (a) | reduced coefficients for ALL states (eps3,Lambda3) of all 113 irreps, 4.04 M numbers, all rho | 0.116 s |
| ndsu3lib (b) | reduced coefficients, HW row only, all 113 irreps, 26 086 numbers, all rho (rho-resolution of Draayer-Akiyama, not the Gaudin eigenbasis) | 0.78 ms |
| slice_engine t_b (a=7) | full HW vector (GT basis, 155 176 numbers) of every copy, Gaudin eigen-copies | 0.185 s |
| slice_engine t_a (a=7) | HW slices of every copy in the Gaudin eigenbasis, no n_W-sized object | 0.128 s |
| jacobi_cg t_float (a=7) | slices of every copy via closed-form Jacobi chains, float64 | 7.7 ms |
| slice_engine t_a / t_b (a=14, [28,14,0]^2, 421 irreps / 2255 copies) | as above | 1.306 s / 2.105 s |

(ndsu3lib was not run beyond a = 7.)

## 6. Verdict

For the highest-weight Clebsch-Gordan data of (7,7)x(7,7) the SU(3)-specific Notre Dame library is not slower than our engine; it is faster or equal in every like-for-like comparison I can make. The HW row (the
coefficients of the coupled highest-weight state, which contain the slice data) takes 0.78 ms for all 113 irreps and verifiably spans exactly the same multiplicity spaces M_mu (deviation 6.7e-16), against 0.128 s for
our slice_engine t_a (about 160x slower than ndsu3lib at this size) and 7.7 ms for jacobi_cg (about 10x slower, and it does not produce the ndsu3lib-style rho-resolved row either). For the "full set"
the comparison is closer: ndsu3lib's 0.116 s produces every reduced coefficient of every state of every irrep (4.0 M numbers, from which all CG coefficients follow via SU(2) CG), while our 0.185 s produces only the
highest-weight vector of each copy (155 k numbers); equal order of magnitude per run, but per number or per object ndsu3lib is far ahead, and I did not time converting its output into explicit GT-basis HW vectors, which
would add Python-side cost. What our engine has that ndsu3lib does not: the Gaudin eigenbasis of the multiplicity space (the SoV copies; ndsu3lib's rho-labels are an arbitrary m x m rotation of it), exact rational arithmetic,
and a structure (V_red tridiagonal, Jacobi chains, Q-system link) that is the research point. The claim "orders of magnitude faster than existing codes" is true against the Munich/general-N code (206 s / 10.9 s) but
not against ndsu3lib: against that, the engine is about 1.5x slower for the full set of per-copy HW vectors and about two orders of magnitude slower for the HW row, with the caveats that the objects differ and the
Python numbers are from a loaded, multithreaded run while ndsu3lib is compiled Fortran.


# Scaling to a = 20 (Lambda = [2a,a,0]^2, SU(3) (a,a)x(a,a); 2026-10-21 night, same laptop)

Everything here is single-threaded (`OPENBLAS_NUM_THREADS=1`), run sequentially; `uptime` load average was 0.9-1.5 (about one other busy core, no heavy job of ours in parallel;
the `load` value printed per repetition is in `ours_raw.txt`). Raw outputs: `timings_scaling_raw.txt` (ndsu3lib, min of 5), `timings_a20_subset.txt` (a = 20 full-set subsets),
`ours_raw.txt` (slice_engine / jacobi_cg), `span_check_a{7,10,14,20,25,30}.txt`, `span_per_irrep_a*.txt` (per-irrep residuals). Driver: `driver.cpp` generalised
(env `SITE=a`, default 7 = old behaviour, reproduced: 0.1159 s / 0.79 ms at a = 7; env `HWONLY=1` skips the full set, `STRIDE`/`OFFSET` time an irrep subset, mode `one` times the largest irrep, prints maxrss);
irrep lists `irreps_a{7,10,14,20,25,30}.txt` from `gen_irreps.py`; wrapper `bench_ours.py`; `span_check.py` takes `a` as 2nd argument (adds direct E^k residual, orthonormal-basis residual, per-irrep file).

## 1. Irrep lists (cross-check)

Racah-Speiser (`slice_engine.racah_speiser`) == closed form (`jacobi_cg.irreps`) as dictionaries for a = 7, 10, 14, 20 (assert in `gen_irreps.py`; a = 25, 30 only the closed form),
and ndsu3lib's `outer_multiplicity` equals m for every irrep at every a (0 mismatches, a = 7...30, `mult` output in `timings_scaling_raw.txt`):

| a | 7 | 10 | 14 | 20 | 25 | 30 |
|---|---|---|---|---|---|---|
| #irreps / #copies (sum m) | 113 / 344 | 221 / 891 | 421 / 2255 | 841 / 6181 | 1301 / 11726 | 1861 / 19871 |

## 2. Timings (min of 5 for ndsu3lib and jacobi_cg [its own inner best-of-5, 3 outer calls]; slice_engine = min of 3 cold `bench_one(a,1)`)

ndsu3lib objects as above: (a) full set = every (eps3,Lambda3) state of every irrep, all rho; (b) HW row = coupled highest-weight state only, all rho. Our objects as in `python-timings.md`
(slice_engine t_a: HW slices in the Gaudin eigenbasis, no n_W-sized object; t_b: + all full HW vectors in the GT basis; jacobi_cg t_float: slices only, float64, forward recurrence).

| a | code / object | time | #calls | #numbers produced |
|---|---|---|---|---|
| 7 | ndsu3lib (a) full set | 0.1151 s | 8 408 | 4.04e6 |
| 7 | ndsu3lib (b) HW row | 0.79 ms | 113 | 2.6e4 |
| 7 | slice_engine t_a / t_b | 0.125 / 0.180 s | - | t_b: 1.55e5 |
| 7 | jacobi_cg t_float | 7.8 ms | - | slices |
| 10 | ndsu3lib (a) | 1.503 s | 31 141 | 5.12e7 |
| 10 | ndsu3lib (b) | 5.52 ms | 221 | 1.63e5 |
| 10 | slice_engine t_a / t_b | 0.404 / 0.612 s | - | t_b: 1.31e6 |
| 10 | jacobi_cg t_float | 19.3 ms | - | slices |
| 14 | ndsu3lib (a) | 19.6 s (single earlier run 21.5 s; 5-rep spread 1.9 s) | 110 405 | 6.09e8 |
| 14 | ndsu3lib (b) | 33.8 ms | 421 | 9.9e5 |
| 14 | slice_engine t_a / t_b | 1.250 / 2.047 s | - | t_b: 1.07e7 |
| 14 | jacobi_cg t_float | 48.5 ms | - | slices |
| 20 | ndsu3lib (a) | **about 3.5-4 s x 100 = 350-400 s, extrapolated, not run in full** (see below) | about 7.9e5 (est.) | about 8.95e9 (est.) |
| 20 | ndsu3lib (b) | 314.7 ms (spread 60 ms) | 841 | 7.17e6 |
| 20 | slice_engine t_a / t_b | 4.34 / 8.18 s | - | t_b: 1.07e8 |
| 20 | jacobi_cg t_float | 136.7 ms | - | slices |
| 25 | ndsu3lib (b) | 1.372 s | 1 301 | 2.53e7 |
| 25 | slice_engine t_a / t_b (1 run) | 9.81 / 20.8 s | - | t_b: 4.66e8 |
| 25 | jacobi_cg t_float (1 run) | 277 ms | - | slices |
| 30 | ndsu3lib (b) | 4.747 s | 1 861 | 7.16e7 |
| 30 | jacobi_cg t_float (1 run) | 490 ms | - | slices |

How the a = 20 full-set number was obtained: one full call set would be about 6 min, so I timed every 8th irrep (offset 0: 106 irreps, 50.2 s, 1.1325e9 coefficients, 56 290 calls;
offset 3: 105 irreps, 48.3 s, 1.1142e9 coefficients, 51 212 calls) and scaled by 841/106 and 841/105 -> 398 s and 387 s (coefficients 8.98e9 and 8.93e9). The same estimator applied at a = 14
(53 irreps each, 2.79 s and 2.72 s) gives 22.3 and 21.7 s against the measured full 21.5 s (single run) / 19.6 s (min of 5), i.e. it overestimates by 1-3% (single-run noise is about 10%).
So the a = 20 estimate is about 350-400 s. A single largest irrep ([67,40,13], rho = 14, 784 calls, 3.30e7 coefficients) takes 1.51 s measured. ndsu3lib's throughput is a steady 31-35 M coefficients/s at every a.
Memory: ndsu3lib max RSS 16.6 MB (a = 7), 17.6 (14), 21.3 (20), 27.4 (25), 38.2 (30) (buffers only, the 8.95e9 numbers are not stored); slice_engine at a = 20 (t_b run): 149 MB; jacobi_cg at a = 20 (including its mp30 checks): 130 MB (Python+numpy baseline about 90 MB).

Ratios (time of our engine / time of ndsu3lib; "-" = ndsu3lib faster; objects differ as stated, so read the pairs, not the single number):

| a | ndsu3lib (b) HW row vs jacobi_cg t_float | ndsu3lib (b) vs slice_engine t_a | ndsu3lib (a) vs slice_engine t_b |
|---|---|---|---|
| 7 | ndsu3lib 10x faster | ndsu3lib 160x faster | ndsu3lib 1.6x faster |
| 10 | ndsu3lib 3.5x faster | ndsu3lib 73x faster | ndsu3lib 2.5x SLOWER |
| 14 | ndsu3lib 1.4x faster | ndsu3lib 37x faster | ndsu3lib 9.6x slower |
| 20 | **jacobi_cg 2.3x faster** | ndsu3lib 14x faster | ndsu3lib about 45x slower (est.) |
| 25 | jacobi_cg 4.9x faster | ndsu3lib 7x faster (t_a 9.81 s vs 1.37 s) | - |
| 30 | jacobi_cg 9.7x faster | - | - |

## 3. Fitted scaling exponents (time ~ a^p; least squares on log-log over a = 7,10,14,20, and local slopes)

| object | fit 7..20 | local 14 -> 20 | local at larger a |
|---|---|---|---|
| ndsu3lib (a) full set | 7.7 | 8.3 (8.1 with the 350 s estimate) | - |
| ndsu3lib (b) HW row | 5.7 | 6.3 | 6.6 (20->25), 6.8 (25->30) |
| slice_engine t_a | 3.4 | 3.5 | 3.7 (20->25) |
| slice_engine t_b | 3.6 | 3.9 | 4.2 (20->25) |
| jacobi_cg t_float | 2.7 | 2.9 | 2.7 (20->25: 0.137->0.277 s), 2.7 (25->30) |

Reasons (plausible, not separately verified): ndsu3lib computes the entire row of (p1,p2,q2) for a given (eps3,Lambda3), so the HW row alone already has about a^6 numbers (26 086 -> 7.2e6 from a = 7 to 20: exponent 5.6)
and the full set about a^8; our engines only touch the m-dimensional multiplicity space (slice) and only the weight space W (dim w about a), hence the much smaller exponents.

## 4. Accuracy of ndsu3lib at large labels (double precision; `span_check.py`, same procedure as section 4)

For every irrep: ndsu3lib HW row (site 2 at its highest weight) x SU(2) CG, / sqrt(GT metric), no phases fitted -> m slice vectors in V[mu - Lambda]; residual = distance of each (relative to its norm)
from the exact PRV space M_mu = {E12^(a+1) v = E23^(a+1) v = 0} of `slice_engine.prv_space` (exact rational basis), which is equivalent to the requested E12^k v = E23^k v = 0 test.
The residual of E^k v itself, measured relative to |v| ||E^k|_W||_2, is also given.

| a | max membership resid. (plain GT basis) | max resid. orthonormal basis | max direct E^k resid. (rel.) | irreps > 1e-8 | worst 1 - cos (principal angle) | worst irreps |
|---|---|---|---|---|---|---|
| 7 | 6.7e-16 | 2.0e-15 | 1.6e-16 | 0 / 113 | 1e-15 | (22,10,10), (21,12,9) |
| 10 | 5.5e-15 | 1.9e-14 | 5.4e-15 | 0 / 221 | 1e-15 | (30,19,11), (30,18,12) |
| 14 | 1.2e-12 | 5.3e-12 | 2.1e-13 | 0 / 421 | 1e-15 | (42,26,16), (42,25,17) |
| 20 | 4.9e-09 | 2.1e-08 | 2.0e-09 | 0 / 841 | 1e-15 | (60,38,22), (60,37,23), (60,36,24) |
| 25 | 2.1e-06 | 1.1e-05 | 1.5e-06 | 9 / 1301 | 2.1e-11 | (75,48,27), (75,49,26) |
| 30 | 7.8e-04 | 4.3e-03 | 7.8e-04 | 31 / 1861 | 2.4e-04 | (90,59,31), (90,58,32), (90,57,33) |

(The "worst 1 - cos" column is 1e-15 up to a = 20 because the principal angle test only resolves 1e-15; its first real deviation shows at a = 25. Per-irrep numbers: `span_per_irrep_a*.txt`.)

Finding: **ndsu3lib's double-precision coefficients lose accuracy rapidly with a: worst-case residual 7e-16 (a = 7), 5e-15 (10), 1.2e-12 (14), 4.9e-9 (20), 2.1e-6 (25), 7.8e-4 (30)**, i.e. a factor about 3.3 per unit of a
(about 12 decades between a = 7 and 30, roughly exponential: 10^(0.5 a)); a = 14 -> 20 loses 4 digits. All a <= 20 irreps still pass 1e-8 (a = 20 is on the edge: the worst is 4.9e-9 in the plain GT basis,
2.1e-8 in the orthonormal basis, so measured in that basis it is already above 1e-8). The worst irreps are always the same family: m1 = 3a (the top row, Lambda3 maximal), second label near 2a - 2: the irreps
with the highest multiplicity m about a-1 and w = m + 1 or so, i.e. the ones with the largest rho_max (the Draayer-Akiyama rho-resolution is the ill-conditioned step there; the singular values of the m ndsu3lib
vectors, normalised, shrink to 0.049 (a = 7), 0.0098 (10), 0.0011 (14), about 1e-5 or below (20, printed as 0.0000), so the rho-labelled basis itself becomes nearly degenerate, which amplifies rounding).
Caveat: the residual is a lower bound on the coefficient error in the worst irrep (membership in M_mu does not see errors inside M_mu); the negative-control sensitivity also degrades
(at a = 20 a sign flip is detected in 188 / 380 non-vacuous irreps only, median residual of a flip 3.6e-11) because the vectors are nearly parallel.
I did not check the coefficients of the full set (all states), only the HW row.

Higher precision: **ndsu3lib has CMake options `-DPRECISION=quad|multi|multiquad` (quad via gfortran REAL128 + libquadmath, multi via MPFUN20-Fort), but they apply only to the SU(3)-SO(3) coefficients**
(INSTALL.md: "for SU(3)-SO(3) reduced coupling coefficients"; the only precision `#if`s are in `ndsu3lib_tools.F90`, precomputed binomial/I/S tables, and none in `ndsu3lib_coupling_canonical.F90`).
I built `-DPRECISION=quad` (`build_quad/`, `bench_driver_quad`, 2 min, gitignored): for the canonical SU(3) > SU(2)xU(1) coefficients used here the output has **the same accuracy** (max residuals 1.246e-12 at a = 14 and 4.928e-9 at a = 20, identical to the double build)
and the time is unchanged (HW row 34.5 ms at a = 14, 326 ms at a = 20, vs 33.8 / 315 ms). So there is no higher-precision mode for this branch; one would have to port the extremal/non-extremal recursion in canonical.F90 to REAL128 or use the SU3lib (Dytrych) alternative.
Not tried: multi/multiquad (needs MPFUN20 download, and would not help for the same reason).

Our float paths, against the same kind of reference (slice_engine: exact rational PRV space; jacobi_cg: 30-digit mpmath pipeline, itself validated against slice_engine's 50-digit eigenvalues up to a = 14, `RESULTS-jacobi-cg.md`):

| a | slice_engine float64: max HW residual (E12, E23, full vectors) | slice_engine float eig residual | jacobi_cg float64 forward recurrence, max rel. dev. | jacobi_cg float64 twisted recurrence |
|---|---|---|---|---|
| 7 | 1.7e-14 | 2.0e-12 | 3.9e-13 | 8.1e-15 |
| 10 | 2.2e-13 | 6.9e-12 | 7.2e-11 | 2.2e-14 |
| 14 | 4.2e-12 | 3.0e-11 | 2.9e-8 | 2.7e-14 |
| 20 | 5.0e-10 | 2.2e-10 | 5.0e-5 | 3.7e-14 |
| 25 | 5.1e-8 | 2.8e-9 | 8.2e-2 | 5.5e-14 |
| 30 | - | - | 72 (breaks down) | 7.4e-14 |

So slice_engine's float64 path also degrades (1.7e-14 -> 5e-10 -> 5e-8 for a = 7, 20, 25: a factor about 2 per unit of a, against about 3.3 for ndsu3lib; cause not investigated, presumably the float64 eigenvectors of V_red and the triangular reconstruction); at a = 20 it is
about 10x more accurate than ndsu3lib in this measure and at a = 25 about 40x (5.1e-8 vs 2.1e-6). **jacobi_cg's plain forward float64 recurrence fails completely by a = 25-30; its twisted variant is accurate to <= 7e-14 at every a tested (a = 30 incl.)**,
the only path here whose accuracy does not degrade. Our own numbers: the exact paths (Fractions, the slice_engine reference files, 30-60-digit mpmath) carry no float error.
(jacobi_cg: the "twisted" accuracy is the max relative deviation of the float64 slices from the 30-digit mp slices, over all 6 181 copies at a = 20.)

## 5. Verdict at large a

* **ndsu3lib's lead does not persist; it shrinks and then reverses against jacobi_cg, and it is overtaken on the full-set comparison against slice_engine already at a = 10.** HW row vs jacobi_cg's slice
  (the like-for-like pair for the slice data): ndsu3lib 10x faster at a = 7, 1.4x at a = 14, jacobi_cg 2.3x faster at a = 20, 4.9x at a = 25, 9.7x at a = 30. Against slice_engine t_a, ndsu3lib's lead falls from 160x (a = 7) to 37x (14), 14x (20), 7x (25).
  The reason is the exponent: ndsu3lib HW row about a^6.5, our slice engines a^3.4 - a^3.7, jacobi_cg a^2.7; extrapolating the local slopes, slice_engine t_a would overtake the ndsu3lib HW row near a = 50 (14x lead at a = 20, exponent gap about 2.9; extrapolation only).
* **Full set**: ndsu3lib's time grows as a^8 (0.12 s -> 19.6 s -> about 6 min for a = 20, 9e9 coefficients); slice_engine t_b a^3.6-3.9 (0.18 -> 2.0 -> 8.2 s). ndsu3lib is 1.6x faster at a = 7, slower from a = 10 on (2.5x), 9.6x slower at a = 14 and about 45x slower at a = 20 (objects differ: 9e9 reduced coefficients
  of all states versus 1.07e8 numbers of the HW vectors; ndsu3lib's per-number rate, 3e7/s, is still excellent, it is just producing a^8 numbers).
* **Accuracy**: ndsu3lib in double precision is trustworthy (<= 1e-12) up to a = 14, borderline at a = 20 (5e-9 to 2e-8) and unusable beyond a = 25-30 (1e-6, 1e-3), with no quad/multi option for this coefficient branch (verified by building it). Our twisted jacobi_cg path
  stays at 1e-14 up to a = 30; slice_engine float64 degrades more slowly than ndsu3lib (5e-10 at a = 20, 5e-8 at a = 25). So for rep sizes beyond (14,14) our engines are faster (jacobi_cg) or equally cheap per useful object AND more accurate; for a <= 10 ndsu3lib remains the faster code.
* What ndsu3lib still does that ours does not: all states / all ordinary CG, rho-resolved reduced coefficients, general (non-self-conjugate, unequal) irrep pairs, a mature library; our engines do the HW/slice data of Lambda x Lambda in the Gaudin eigenbasis only.
* Caveats: timings single-threaded on a laptop at load about 1; the a = 20 full set is extrapolated from two subsets (+-10%); the accuracy statement is the residual of the HW row against the exact multiplicity space (needs the PRV space; error inside M_mu not seen);
  ndsu3lib could be run in quad only by rewriting its canonical routine; I did not test OpenMP or other ndsu3lib build options for speed.
