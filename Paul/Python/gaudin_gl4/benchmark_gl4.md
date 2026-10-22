# Benchmark: the gl(4) slice engine against the Munich code (compiled 2026-10-22 from RESULTS-G1, -G3, -G7; no new runs)

*Standalone page for `Paul/Mathematica/Gaudin/gl4-slice-engine-notes.{tex,md,pdf}` (Section 8). Every number is copied from a section of a RESULTS file listed in the last line; the few derived numbers are marked "computed here". Engine: Python 3.13, numpy 2.5.2, scipy 1.18.1, mpmath 1.3, single process, i7-12700KF (12 cores / 20 threads, 32 GB, Windows 11), minimum of cold runs. Munich: `Paul/Cpp/AlexCGC/results/hw_timings.txt`, 13th-gen i5-1345U, g++ 11.4, LAPACK 3.10, one run per irrep, generated 2026-07-22, the code as distributed. **Parts A and B compare the i7 desktop with the July Munich timings (two machines, two code bases: read those ratios as orders of magnitude); parts E and F are the same-machine comparison on the i5 laptop (RESULTS-G7), original and repaired Munich against the engine, which supersedes the cross-machine ratios.***

## A. SU(4), multiplicity-free family $[S,S,0,0]^{\otimes 2}$: all highest-weight vectors of all irreps

$(S+1)(S+2)/2$ irreps, all of multiplicity 1, so no missing-label operator enters. Both sides produce the same object: the HW vector (all $n_W$ coefficients) of every irrep. **Munich** = sum over the irreps of `compute_highest_weight_coeffs()` (the raising-operator null space; no Littlewood–Richardson step, no descendants). **$t_a$** = the part that never touches $n_W$: exact GT generators (simple + Cartan) + Racah–Speiser + exact PRV spaces. **$t_b$** = $t_a$ + metric, orthonormal float matrices, dense blocks and left inverses of the level recursion + the float64 reconstruction of every full vector. "worst resid." = largest HW residual over all irreps (all $E_{ij}$ for $\dim\le400$, the simple raising operators otherwise; untimed). Quietest of three blocks of the whole table (block 1 of G1 §6).

| S | dim V | irreps | max w | max $n_W$ | Munich [s] | $t_a$ [s] | $t_b$ [s] | Munich / $t_b$ | entries written | worst resid. |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 6 | 3 | 1 | 6 | 0.0000279 | 0.0002 | 0.0007 | 0.04 | 9 | 1.9e-16 |
| 2 | 20 | 6 | 2 | 22 | 0.000858 | 0.0007 | 0.0018 | 0.48 | 48 | 2.9e-16 |
| 3 | 50 | 10 | 2 | 62 | 0.0107 | 0.0016 | 0.0040 | 2.7 | 186 | 5.9e-16 |
| 4 | 105 | 15 | 3 | 147 | 0.0574 | 0.0034 | 0.0088 | 6.5 | 580 | 8.2e-16 |
| 5 | 196 | 21 | 3 | 308 | 0.355 | 0.0064 | 0.0179 | 20 | 1 544 | 1.8e-15 |
| 6 | 336 | 28 | 4 | 588 | 2.49 | 0.0110 | 0.0326 | 76 | 3 646 | 1.8e-15 |
| 7 | 540 | 36 | 4 | 1 044 | 15.3 | 0.0171 | 0.0572 | 268 | 7 838 | 1.3e-15 |
| 8 | 825 | 45 | 5 | 1 749 | 80.2 | 0.0258 | 0.0932 | 860 | 15 625 | 5.6e-15 |
| 9 | 1 210 | 55 | 5 | 2 794 | n/a | 0.0381 | 0.1417 | – | 29 277 | 4.6e-15 |
| 10 | 1 716 | 66 | 6 | 4 290 | n/a | 0.0579 | 0.2276 | – | 52 090 | 2.8e-15 |
| 11 | 2 366 | 78 | 6 | 6 370 | n/a | 0.0847 | 0.3212 | – | 88 700 | 1.0e-14 |
| 12 | 3 185 | 91 | 7 | 9 191 | n/a | 0.1027 | 0.4427 | – | 145 456 | 1.3e-14 |
| 14 | 5 440 | 120 | 8 | 17 816 | n/a | 0.1780 | 0.8831 | – | 356 052 | 4.6e-14 |
| 16 | 8 721 | 153 | 9 | 31 977 | n/a | 0.3218 | 1.6202 | – | 787 257 | 2.6e-13 |
| 20 | 19 481 | 231 | 11 | 86 779 | n/a | 0.8242 | 4.5883 | – | 3 064 644 | 3.9e-11 |

`n/a`: not in `hw_timings.txt` (Munich's code was not run for $S\ge9$ in July; the same-machine run of part E has the original at $S\le8$, 33 GB needed at $S=9$, and the repaired variant to $S=11$). Munich per-S totals: $2.785\times10^{-5}$, $8.577\times10^{-4}$, 0.010718, 0.057397, 0.355094, 2.486182, 15.306753, 80.186760 s for $S=1,\dots,8$ (`hw_timings.txt`).

**Run-to-run noise** (three blocks of the whole table, $t_a$ / $t_b$ [s]; block 0 shared the machine with another job, blocks 1–2 are re-runs):

| S | block 0 | block 1 (quoted) | block 2 |
| ---: | ---: | ---: | ---: |
| 4 | 0.0044 / 0.0105 | 0.0034 / 0.0088 | 0.0034 / 0.0090 |
| 8 | 0.0291 / 0.1012 | 0.0258 / 0.0932 | 0.0257 / 0.0918 |
| 12 | 0.1217 / 0.4964 | 0.1027 / 0.4427 | 0.1061 / 0.4653 |
| 16 | 0.3780 / 1.7573 | 0.3218 / 1.6202 | 0.3325 / 1.6544 |
| 20 | 0.9270 / 4.8989 | 0.8242 / 4.5883 | 0.7432 / 4.6196 |

Differences are at most 3% between the two re-runs and at most 15% for block 0. The sum of $t_b$ over $S=4..12$ is 1.522, 1.343, 1.362 s in the three blocks.

**Munich's per-irrep rows at $S=8$** (Munich label $\mu_4=0$; engine "warm" = after the per-$\lambda'$ setup of 16 ms, i.e. `prv_space` + `hw_from_slice` of that irrep only; selected from the 45 rows of G1 §6):

| $\mu$ (Munich) | dim $V_\mu$ | $n_W$ | Munich hw-solve [s] | engine, warm [s] | speed-up |
| :--- | ---: | ---: | ---: | ---: | ---: |
| (16 16 0 0), top | 8 721 | 1 | 0.024730 | 0.00024 | 101 |
| (16 8 8 0) | 12 825 | 9 | 0.042272 | 0.00027 | 157 |
| (14 7 7 0) | 7 344 | 48 | 0.110981 | 0.00046 | 239 |
| (12 6 6 0) | 3 920 | 142 | 0.293158 | 0.00079 | 370 |
| (8 8 0 0) | 825 | 259 | 0.535824 | 0.00103 | 519 |
| (6 3 3 0) | 300 | 896 | 3.235421 | 0.00247 | 1 313 |
| (4 2 2 0) | 84 | 1 260 | 7.219066 | 0.00315 | 2 291 |
| (2 1 1 0) | 15 | 1 584 | 13.520279 | 0.00434 | 3 119 |
| **(0 0 0 0), singlet** | 1 | 1 749 | **18.702406** | **0.00454** | **4 116** |
| sum over the 45 irreps | | | 80.187 | 0.0520 | |

The singlet, Munich's slowest row, takes 4.5 ms warm and 44.7 ms cold (a single irrep including generators and setup); cold/warm for the singlet of $[6,6,0,0]^{\otimes2}$: 17.6 / 2.0 ms, of $[10,10,0,0]^{\otimes2}$: 95.4 / 9.4 ms. At $S=8$ the breakdown of the quoted block is: generators 22.7 ms, Racah–Speiser 2.0 ms, 45 exact PRV spaces 1.2 ms ($t_a=25.8$ ms); float setup 16.0 ms (489 SVD left inverses), reconstruction of the 45 vectors (15 625 coefficients) 51.4 ms ($t_b=93.2$ ms).

## B. Multiplicity families: canonical copies (the three missing-label operators applied, jointly diagonalised)

Here $t_a$ also contains the reduction of $V_3,V_4^a,V_4^b$ (once per $\lambda'$, counted cold), the **exact** (Python fractions) application of the three reduced operators to the PRV basis of every irrep, and the joint eigen-decomposition (mpmath, 50 digits; or numpy float64, columns marked $f$); $t_b$ adds the recursion setup and the float64 reconstruction of every full copy. Slices only for the rows marked $*$. `copies` = Σ over irreps of $m_\mu$. Minimum over 5 / 3 / 2 cold runs for $\dim V\le500$ / $\le5000$ / larger; one run for rows $*$ (G3 §9).

| family | a | $\lambda$ | dim V | irreps | copies | max m | max w | max $n_W$ | $t_a$ [s] | $t_b$ [s] | $t_a^{f}$ [s] | $t_b^{f}$ [s] | worst HW resid. |
| :--- | ---: | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| $[2a,a,0,0]^{\otimes2}$ | 1 | [2,1,0,0] | 20 | 7 | 8 | 2 | 2 | 20 | 0.009 | 0.010 | 0.008 | 0.009 | 3.5e-16 |
| | 2 | [4,2,0,0] | 126 | 24 | 32 | 3 | 4 | 174 | 0.041 | 0.046 | 0.035 | 0.041 | 1.7e-15 |
| | 3 | [6,3,0,0] | 480 | 58 | 90 | 4 | 7 | 920 | 0.167 | 0.193 | 0.142 | 0.168 | 6.3e-15 |
| | 4 | [8,4,0,0] | 1 375 | 115 | 205 | 5 | 11 | 3 541 | 0.518 | 0.611 | 0.432 | 0.525 | 9.1e-15 |
| | 5 | [10,5,0,0] | 3 276 | 201 | 406 | 6 | 15 | 10 968 | 1.375 | 1.659 | 1.147 | 1.431 | 1.5e-14 |
| | 6 | [12,6,0,0] | 6 860 | 322 | 728 | 7 | 20 | 29 028 | 3.160 | 3.885 | 2.614 | 3.338 | 9.8e-14 |
| | 7 $*$ | [14,7,0,0] | 13 056 | 484 | 1 212 | 8 | 26 | 68 208 | 6.355 | – | 5.251 | – | – |
| | 8 $*$ | [16,8,0,0] | 23 085 | 693 | 1 905 | 9 | 33 | 146 025 | 11.986 | – | 9.887 | – | – |
| $[3a,2a,a,0]^{\otimes2}$ | 1 | [3,2,1,0] | 64 | 16 | 30 | 4 | 4 | 152 | 0.049 | 0.052 | 0.035 | 0.038 | 2.4e-15 |
| | 2 | [6,4,2,0] | 729 | 76 | 297 | 12 | 15 | 4 743 | 1.467 | 1.515 | 0.709 | 0.757 | 6.5e-15 |
| | 3 | [9,6,3,0] | 4 096 | 213 | 1 592 | 27 | 32 | 59 008 | 23.418 | 23.904 | 6.422 | 6.909 | 3.9e-14 |
| | 4 $*$ | [12,8,4,0] | 15 625 | 459 | 5 945 | 51 | 65 | 426 425 | 257.550 | – | 37.152 | – | – |

**Breakdown of $t_a$ and $t_b$** [s] (G3 §9; $t_{\rm red}$ = reduction of all three operators for this $\lambda'$, universal polynomials 0.3–1.2 ms included; $t_{\rm ops}$ = exact application; $t_{\rm eig}$ = joint eigen-decomposition, 50 digits (mp) and float64 (f); $t_{\rm gen}$ = simple + Cartan GT generators, $t_{\rm gen,na}$ = the six non-adjacent generators, which the reduced words use):

| family, a | $t_{\rm gen}$ | $t_{\rm gen,na}$ | $t_{\rm rs}$ | $t_{\rm prv}$ | $t_{\rm red}$ | $t_{\rm ops}$ | $t_{\rm eig}$ (mp) | $t_{\rm eig}$ (f) | $t_{\rm setup}$ | $t_{\rm recon}$ |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| $[2a,a,0,0]$, 1 | 0.0005 | 0.0003 | 0.0001 | 0.0001 | 0.0041 | 0.0029 | 0.0007 | 0.0002 | 0.0006 | 0.0003 |
| 2 | 0.0034 | 0.0030 | 0.0003 | 0.0005 | 0.0039 | 0.0236 | 0.0064 | 0.0007 | 0.0021 | 0.0032 |
| 3 | 0.0128 | 0.0163 | 0.0009 | 0.0025 | 0.0043 | 0.1032 | 0.0268 | 0.0017 | 0.0074 | 0.0187 |
| 4 | 0.0406 | 0.0518 | 0.0024 | 0.0100 | 0.0041 | 0.3192 | 0.0897 | 0.0039 | 0.0184 | 0.0751 |
| 5 | 0.1309 | 0.1416 | 0.0050 | 0.0318 | 0.0043 | 0.8253 | 0.2362 | 0.0078 | 0.0434 | 0.2410 |
| 6 | 0.2214 | 0.4296 | 0.0108 | 0.0913 | 0.0042 | 1.8423 | 0.5606 | 0.0140 | 0.0882 | 0.6367 |
| 7 $*$ | 0.5066 | 0.7687 | 0.0740 | 0.2299 | 0.0050 | 3.6442 | 1.1265 | 0.0224 | – | – |
| 8 $*$ | 0.9188 | 1.6161 | 0.0317 | 0.5523 | 0.0048 | 6.7277 | 2.1341 | 0.0353 | – | – |
| $[3a,2a,a,0]$, 1 | 0.0016 | 0.0016 | 0.0002 | 0.0003 | 0.0039 | 0.0273 | 0.0141 | 0.0007 | 0.0013 | 0.0016 |
| 2 | 0.0213 | 0.0300 | 0.0012 | 0.0065 | 0.0040 | 0.6412 | 0.7630 | 0.0048 | 0.0101 | 0.0380 |
| 3 | 0.1343 | 0.2525 | 0.0067 | 0.0793 | 0.0044 | 5.9174 | 17.0233 | 0.0274 | 0.0553 | 0.4313 |
| 4 $*$ | 0.7057 | 1.2213 | 0.0208 | 0.6978 | 0.0041 | 34.3659 | 220.5339 | 0.1364 | – | – |

Reading: the reduction costs about 4 ms per $\lambda'$ at every size; the exact operator application ($t_{\rm ops}$) is 56–62% of $t_a$ on the scaling family (1.84 s against 0.64 s for the whole level recursion at $a=6$, i.e. about three times the recursion), and for $m\ge10$ the 50-digit eigen-decomposition dominates (17.0 of 23.4 s at $[9,6,3,0]^{\otimes2}$, 221 of 258 s at $[12,8,4,0]^{\otimes2}$); the float64 eigen-decomposition is 40 ($m\le7$) to 1600 ($m\le51$) times cheaper, accurate to $2.1\times10^{-15}$ in the eigenvalues and $7.0\times10^{-12}$ (relative) in the off-diagonal entries of the reduced operators in the eigenbasis (at $m=27$ and $51$). Per operator and copy the gl(4) reduced operators cost 678 µs at $[10,5,0,0]^{\otimes2}$ ($w\le15$), 2.6 times the 24-word gl(3) operator at the same $w$ (G3 §9). The Munich SU(4) data are multiplicity-free, so there is **no Munich comparison for this table**; its null-space method would need an $n_W\times n_W$ null space per irrep ($n_W$ up to 29 028 at $a=6$, $n_W^3=2.4\times10^{13}$) and was not run.

**Without operators** (G1 §7: the PRV basis vectors in place of the canonical copies, `[2a,a,0,0]^2` and `[3a,2a,a,0]^2`; $t_a$ / $t_b$ [s]): $a=1..6$: 0.0007/0.0019, 0.0044/0.0103, 0.0163/0.0428, 0.0541/0.1515, 0.1459/0.4783, 0.3901/1.1084; $[3a,2a,a,0]$, $a=1,2,3$: 0.0020/0.0048, 0.0295/0.0761, 0.2347/0.7275 (8.9·10⁶ coefficients at $a=3$).

## C. The gl(3) context (three lines)

* $[14,7,0]^{\otimes2}$ (113 irreps, 344 copies): slice engine $t_a=0.10$ s, $t_b=0.17$ s on the i7 (Jacobi chain, slices only: 0.006 s) against the Munich code as distributed 206 s and repaired (compact columns, $AA^T$ eigenvectors, OpenBLAS) 10.9 s on the i5 (`RESULTS-slice-engine.md` §5; `Slice-engine-notes` Table 2).
* On the same i5 laptop (HANDOFF, ndsu3lib block): ndsu3lib (SU(3)-specific, Herko et al., EPJA 62 (2026) 95) needs 0.116 s for all reduced coefficients of all states and 0.80 ms for the HW row only; our slice engine 0.128 / 0.185 s ($t_a$ / $t_b$), the Jacobi chain 7.7 ms there; the Jacobi chain overtakes ndsu3lib's HW row at about $a=15$.
* So the advantage holds against **general-$N$ null-space codes only** (Munich ×60–1200 for gl(3); for SU(4) at $S=8$, on one machine, ×908 against the original and ×38 against the repaired code, part E); there is no ndsu3lib analogue for $N=4$ (one missing label, a prescribed zero pattern and closed-form seeds are SU(3)-specific); the SU(4)-specific CGSU4 (Kuhn–Walliser 2008) was not obtained and not benchmarked.

## D. How to read this page

* **What each column measures.** Munich = the highest-weight solve only (one run per irrep, i5, the code as distributed; `cgc_hw_bench`). $t_a$ = everything that does not scale with $n_W$ (generators of both sites, Racah–Speiser incl. GT weight multiplicities, exact PRV spaces; in table B also the reduction, the exact operator application and the joint eigen-decomposition). $t_b$ = $t_a$ + float64 setup of the recursion + writing out every full vector (`entries written` = Σ $n_W\cdot m$ coefficients).
* **$n_W$-free steps.** Generators, Racah–Speiser, PRV null spaces, reduction, operator application on $M_\mu$, eigen-decomposition (steps 1–4 of the engine): no object of the size of $n_W$ is touched. The recursion (step 5) is bound by the output size: $t_b\propto$ entries $\propto S^{5.9}$ (580 at $S=4$, 52 090 at $S=10$, $3.06\times10^6$ at $S=20$); time per coefficient 15 µs at $S=4$, 4.4 µs at $S=10$, 1.5 µs at $S=20$ (G1 §6).
* **Exact and float.** Generators, Racah–Speiser, PRV spaces and the reduced operators are exact (Python fractions); the recursion and the eigen-decomposition of the multiplicity-free table are float64. The float64 recursion's HW residual stays at or below $1.3\times10^{-14}$ up to $S=12$ and grows by ×2.4–3.7 per unit of $S$ beyond $S=14$ ($2.6\times10^{-13}$ at $S=16$, $3.9\times10^{-11}$ at $S=20$, $7.5\times10^{-9}$ at $S=24$); the exact recursion rounded at the end keeps $\sim10^{-15}$ at 1.1–15.6 s per vector ($n_W=2.9\times10^4$ to $2\times10^5$, not vectorised) (G1 §4.4).
* **Machine difference (parts A and B).** i5-1345U (Ubuntu, g++ 11.4, LAPACK 3.10, one run) against i7-12700KF (Windows 11, Python, minimum of five cold runs); single-thread speeds differ by a factor of order 1–2 (G1 §6), so the ratios of parts A and B are orders of magnitude, not a measured speed-up. Part E removes the machine difference (one laptop for all codes) and adds the repaired Munich variant: for SU(4) the repair buys ×1.9 ($S=4$) to ×22.6 ($S=8$), against ×19 for gl(3) at $[14,7,0]^{\otimes2}$.
* **Python-overhead regime.** For $S\le2$ the Munich C++ is faster in absolute terms (28 µs against 0.7 ms at $S=1$; 0.86 ms against 1.8 ms at $S=2$); the crossover is at $S=3$ (10.7 ms against 4.0 ms).
* **Scaling.** Munich's total grows by ×5.2–7.0 per unit of $S$ between $S=4$ and $S=8$ (0.0574 → 0.355 → 2.49 → 15.3 → 80.2 s; cost $O(n_W^3)$ per irrep); the engine's $t_b$ grows by ×2.0 ($S=4\to5$) decreasing to ×1.4 ($S=11\to12$). Computed here from the table: against $t_a$ the ratio at $S=8$ is $80.19/0.0258\approx3100$; against $t_b$ 860 (cross-machine; ×908 on one machine, part E).
* **Munich HW vectors compared (RESULTS-G7 §3, §5).** On 22 October the HW coefficients of the original Munich code were compared with ours on the i5 laptop. Multiplicity-free: 42 of 42 irreps MATCH up to one global sign (all 34 irreps of $[S,S,0,0]^{\otimes2}$ for $S\le4$ and eight at $S=6,8$; $(16,16,0,0)$ at $S=8$ via the repaired profiler's HW output), maximal deviation $1.4\times10^{-15}$. Multiplicity sectors: 100 sectors (sum of multiplicities 210, $m\le12$: 94 sectors of eight products run through the original code, six sectors of $[6,4,2,0]^{\otimes2}$ from the profiler's 12-decimal output), the spans of Munich's Gram–Schmidt copies and of our canonical copies agree to $2.4\times10^{-15}$ ($1.6\times10^{-12}$ for the six profiler sectors); the rotation between the two orthonormal bases is never a signed permutation (38 sectors with $m\ge2$). The GT bases coincide (pattern order = Munich's $Q(M)$, trivial sign gauge).

## E. Same machine: the multiplicity-free family $[S,S,0,0]^{\otimes2}$ against the original and the repaired Munich code (RESULTS-G7 §2)

All numbers of this part and of part F were measured on one laptop: 13th-gen i5-1345U (12 threads), 15 GB, Ubuntu 22.04, g++ 11.4; `LC_ALL=C`, `OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1`, `/usr/bin/time -v`, `timeout 1800`, `ulimit -v 12000000`. **orig** = `cgc_hw_bench` (Munich as distributed, reference LAPACK/BLAS 3.10); **eig** / **compact** = the repaired profiler `cgc_hw_prof_ob` with `HW_MODE=eig|compact` (OpenBLAS 0.3.34, one thread; eig takes the HW vectors as eigenvectors of $AA^{T}$, compact stores only the used columns of the same `dgesvd`); **ours** = the slice engine, Python 3.10.12, single process ($t_a$, $t_b$ as in part A; "multithread" = the same code with threads allowed). One run per irrep (reps 1) unless stated; "min of 3" for orig and eig at $S=4,6,8$. The load average sat at 1.0–1.3 before and after every batch (our one process plus about one background load: Firefox, VS Code). Ratios use the best (min) Munich number and our single-thread $t_b$. [s]:

| S | Munich published (July, same laptop) | orig now (reps 1) | orig min of 3 | eig (reps 1) | eig min of 3 | compact (reps 1) | ours $t_a$ | ours $t_b$ | ours $t_b$, multithread | orig / $t_b$ | eig / $t_b$ | compact / $t_b$ | max RSS MB orig / eig / compact |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2.785e-05 | 5.388e-05 | – | 0.0001385 | – | 5.493e-05 | 0.0003 | 0.0008 | 0.0009 | 0.1 | 0.2 | 0.1 | 5 / 7 / 7 |
| 2 | 0.0008577 | 0.0006639 | – | 0.0006106 | – | 0.0006228 | 0.0008 | 0.0019 | 0.0023 | 0.3 | 0.3 | 0.3 | 5 / 7 / 7 |
| 3 | 0.01072 | 0.007665 | – | 0.004528 | – | 0.006763 | 0.0018 | 0.0043 | 0.0051 | 1.8 | 1.1 | 1.6 | 6 / 8 / 8 |
| 4 | 0.0574 | 0.05368 | 0.03126 | 0.02891 | 0.02583 | 0.03526 | 0.0039 | 0.0096 | 0.0107 | 3.3 | 2.7 | 3.7 | 18 / 8 / 8 |
| 5 | 0.3551 | 0.3533 | – | 0.09284 | – | 0.1507 | 0.0076 | 0.0189 | 0.0196 | 19 | 4.9 | 8.0 | 96 / 11 / 11 |
| 6 | 2.486 | 2.481 | 2.447 | 0.2605 | 0.2271 | 0.7541 | 0.0131 | 0.0349 | 0.0347 | 70 | 6.5 | 22 | 515 / 19 / 19 |
| 7 | 15.31 | 15.00 | – | 1.018 | – | 4.139 | 0.0200 | 0.0588 | 0.0604 | 255 | 17 | 70 | 2338 / 42 / 42 |
| 8 | 80.19 | 84.62 | 84.36 | 3.738 | 3.518 | 22.50 | 0.0299 | 0.0929 | 0.0947 | 908 | 38 | 242 | 9116 / 98 / 102 |
| 9 | not run | not run (33 GB) | – | 15.47 | – | not run | 0.0442 | 0.1421 | 0.1431 | – | 109 | – | – / 230 / – |
| 10 | – | – | – | 56.83 | – | – | 0.0677 | 0.2203 | 0.2254 | – | 258 | – | – / 519 / – |
| 11 | – | – | – | 201.8 | – | – | 0.0894 | 0.3246 | 0.3715 | – | 622 | – | – / 1119 / – |
| 12 | – | – | – | not run (skipped on a change of priorities) | – | – | 0.1562 | 0.4842 | 0.5112 | – | – | – | – |

The original reproduces the July table to $1.0\pm0.06$ for $S\ge4$ ($S=8$: 84.6 s now, 80.2 s then); the small-$S$ rows are microsecond noise. Our code does not use more than one core meaningfully (multithread $t_b$ equals single-thread $t_b$); our $t_b$ on this laptop at $S=8$, 0.0929 s, is identical to the i7 table of part A. Consistency check of the profiler's `HW_MODE=orig` (`proforig_S4`, `proforig_S6`): 0.0449 s and 1.759 s against 0.0537 s and 2.481 s of `cgc_hw_bench` (OpenBLAS against reference BLAS, so not identical code speed). **The repair factor orig / eig in SU(4)** is 1.9 ($S=4$), 3.8, 9.5, 14.7 ($S=7$) and 22.6 ($S=8$; 24 with the best-of-3 numbers), against ×19 for gl(3) at $[14,7,0]^{\otimes2}$; the repaired code does not hit a memory wall ($S=9,10,11$ need 230, 519, 1119 MB, the original 33 GB at $S=9$) and $S=11$ takes 202 s (eig grows by about ×3.6 per step of $S$). The compact variant gains ×3.6–3.9 over the original at $S=8$ and uses 102 MB instead of 9.1 GB. **The "×860" of part A (two machines) becomes ×908 on one machine** (orig 84.4–84.6 s against $t_b=0.0929$ s at $S=8$), while the **repaired Munich (eig) is ×38** slower than ours at $S=8$ and ×109, ×258, ×622 at $S=9,10,11$: the ratio keeps growing because eig is roughly cubic in $n_W$ while $t_b$ is almost $n_W$-free.

## F. Same machine: multiplicity products (RESULTS-G7 §1)

Whole products $P\otimes P$, $P=[\ell_1,\ell_2,\ell_3,\ell_4]$: Munich HW-solve times summed over all irreps; one run per irrep for $n_W$ above the repetition thresholds, minimum of 5 / 3 / 1 repetitions for $n_W<300$ / $<1000$ / larger (the number is in the per-irrep files `products_*.txt`). Guards: orig only if $n_W\dim^2\cdot8$ B $\le10$ GB (and not at all at $[6,4,2,0]^{\otimes2}$ and $[9,6,3,0]^{\otimes2}$, as briefed); eig only if $n_W\le6000$; compact only if $n_W\le1500$. **ours** = `slice_ops_gl4.bench_one` (minimum of 2–3 cold runs; the whole product: all irreps, all copies, canonical copies with the three missing-label operators and the joint diagonalisation, **full vectors** in $t_b$), Python 3.10, single thread. [s]:

| product $P$ | irreps | $\sum m$ | max $m$ | max $n_W$ | orig | eig | compact | ours $t_a$ | ours $t_b$ | orig / $t_b$ | eig / $t_b$ | max RSS MB orig / eig / compact |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| [2,1,1,0] | 6 | 7 | 2 | 21 | 0.00024 | 0.00030 | 0.00031 | 0.0107 | 0.0119 | 0.02 | 0.02 | 5 / 7 / 7 |
| [3,2,1,0] | 16 | 30 | 4 | 152 | 0.0168 | 0.0110 | 0.0133 | 0.0609 | 0.0643 | 0.26 | 0.17 | 10 / 8 / 8 |
| [4,2,0,0] | 24 | 32 | 3 | 174 | 0.0558 | 0.0339 | 0.0426 | 0.0510 | 0.0575 | 0.97 | 0.59 | 27 / 9 / 9 |
| [4,2,1,0] | 27 | 56 | 5 | 368 | 0.291 | 0.0553 | 0.131 | 0.147 | 0.154 | 1.9 | 0.36 | 61 / 12 / 12 |
| [4,3,1,0] | 31 | 69 | 5 | 577 | 0.944 | 0.107 | 0.346 | 0.187 | 0.197 | 4.8 | 0.54 | 143 / 18 / 17 |
| [6,3,0,0] | 58 | 90 | 4 | 920 | 12.7 | 1.03 | 3.06 | 0.212 | 0.238 | 53 | 4.3 | 1631 / 33 / 32 |
| [6,4,2,0] | 76 | 297 | 12 | 4743 | not run | 39.4 | 8.02 (65 irreps, $n_W\le1500$) | 1.74 | 1.79 | n/a | 22 | – / 565 / 54 |
| [8,4,0,0] | 115 | 205 | 5 | 3541 | 73.8 (80 irreps with ≤ 10 GB) | 40.3 | 35.0 (101 irreps) | 0.693 | 0.785 | ≥ 94 on the 80 irreps | 51 | 9442 / 325 / 79 |
| [9,6,3,0] | 213 | 1592 | 27 | 59008 | not run | 385 (166 irreps, $n_W\le6000$) | 120 (128 irreps) | 26.2 | 26.8 | n/a | 14 (partial) | – / 873 / 179 |

Partial-row comparisons (same irreps): $[8,4,0,0]^{\otimes2}$ on the 80 irreps all three ran: orig 73.8, eig 7.57, compact 8.89 s, ours $t_b$ for the full product 0.785 s, so on those rows orig / $t_b\ge94$. $[6,4,2,0]^{\otimes2}$: eig 2.75 s against compact 8.02 s on the 65 common irreps (compact is slower than eig, and eig grows steeply with the largest $n_W$: its biggest sectors $n_W=4743, 4050, 3456$ take 13.3, 8.0, 4.9 s; compact does not run there). $[9,6,3,0]^{\otimes2}$: on the 128 common irreps eig 108 s against compact 120 s; the 47 irreps with $n_W>6000$ (up to 59 008) were not run by any Munich variant (the dense $n_W\times n_W$ matrix of eig would need up to 28 GB and cubic time), while our engine does all 213 irreps with all 1592 copies in 26.8 s ($t_a$ 26.2 s: the exact missing-label reduction dominates, not $n_W$). So **on multiplicity products the original Munich code is slower than ours from $[4,2,1,0]^{\otimes2}$ on (×1.9), ×5 at $[4,3,1,0]^{\otimes2}$, ×53 at $[6,3,0,0]^{\otimes2}$, and exceeds the memory limit** ($n_W\dim^2\cdot8$ B up to 20 GB at $[6,4,2,0]^{\otimes2}$, 54 GB at $[8,4,0,0]^{\otimes2}$, 8 TB at $[9,6,3,0]^{\otimes2}$); **the repaired eig variant is faster than ours up to $[4,3,1,0]^{\otimes2}$ and slower from $[6,3,0,0]^{\otimes2}$ (×4.3) on, ×22 at $[6,4,2,0]^{\otimes2}$.** Caveats that go with these numbers: (i) Munich's hw-solve produces only the $m$ orthonormal HW vectors per irrep (not canonical: for $m\ge2$ an arbitrary orthonormal basis of the multiplicity space), ours produces canonical copies (joint eigenvectors of the Gaudin missing-label operators) and the full vectors; (ii) ours is Python on a different code path from a very different algorithm, Munich is C++; (iii) for tiny products our fixed overhead (reduction of the operators) dominates ($t_a$ 0.011 s at $[2,1,1,0]^{\otimes2}$).

**Sources.** Parts E and F: `RESULTS-G7-munich-comparison.md` §1–2 (the correctness comparison: §3, §5), `Paul/Cpp/AlexCGC/results/su4-same-machine/` (`rect_table.md`, `products_summary.txt`, raw Munich outputs, `load_log.txt`, `environment.txt`); no number of parts E and F was recomputed here. Parts A–D: Table A and the noise table: `RESULTS-G1-slice-engine-gl4.md` §6 (block 1 of 3; per-irrep table and singlet cold/warm; breakdown at $S=8$); Munich column: `Paul/Cpp/AlexCGC/results/hw_timings.txt` (totals per $S$) and `README.md`; residual growth, exact-rounded path: G1 §4.4; Munich-basis identity and the 19 commands: G1 §9; PRV-basis numbers for the families: G1 §7. Table B, the breakdown and the cost remarks: `RESULTS-G3-slice-ops-gl4.md` §9 (all rows; the $\ast$ rows are the extension rows). Table C: `Paul/Python/gaudin_sl3/RESULTS-slice-engine.md` §5, `Slice-engine-notes.md` Table 2, `HANDOFF.md` ("DONE (2026-10-21 night, Ubuntu): ndsu3lib benchmark"). Only derived number: the ratio 3100 (computed here).
