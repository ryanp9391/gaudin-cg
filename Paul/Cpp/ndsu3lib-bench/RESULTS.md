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
