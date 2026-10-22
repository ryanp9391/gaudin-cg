# Benchmark: the gl(4) slice engine against the Munich code (compiled 2026-10-22 from RESULTS-G1, -G3; no new runs)

*Standalone page for `Paul/Mathematica/Gaudin/gl4-slice-engine-notes.{tex,md,pdf}` (Section 8). Every number is copied from a section of a RESULTS file listed in the last line; the few derived numbers are marked "computed here". Engine: Python 3.13, numpy 2.5.2, scipy 1.18.1, mpmath 1.3, single process, i7-12700KF (12 cores / 20 threads, 32 GB, Windows 11), minimum of cold runs. Munich: `Paul/Cpp/AlexCGC/results/hw_timings.txt`, 13th-gen i5-1345U, g++ 11.4, LAPACK 3.10, one run per irrep, generated 2026-07-22, the code as distributed. **Two different machines, two different code bases; read the ratios as orders of magnitude.***

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

`n/a`: not in `hw_timings.txt` (Munich's code was not run for $S\ge9$). Munich per-S totals: $2.785\times10^{-5}$, $8.577\times10^{-4}$, 0.010718, 0.057397, 0.355094, 2.486182, 15.306753, 80.186760 s for $S=1,\dots,8$ (`hw_timings.txt`).

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
* So "orders of magnitude" holds against **general-$N$ null-space codes only** (Munich ×60–1200 for gl(3), ×860 against $t_b$ at SU(4) $S=8$); there is no ndsu3lib analogue for $N=4$ (one missing label, a prescribed zero pattern and closed-form seeds are SU(3)-specific), and no repaired Munich run exists for SU(4).

## D. How to read this page

* **What each column measures.** Munich = the highest-weight solve only (one run per irrep, i5, the code as distributed; `cgc_hw_bench`). $t_a$ = everything that does not scale with $n_W$ (generators of both sites, Racah–Speiser incl. GT weight multiplicities, exact PRV spaces; in table B also the reduction, the exact operator application and the joint eigen-decomposition). $t_b$ = $t_a$ + float64 setup of the recursion + writing out every full vector (`entries written` = Σ $n_W\cdot m$ coefficients).
* **$n_W$-free steps.** Generators, Racah–Speiser, PRV null spaces, reduction, operator application on $M_\mu$, eigen-decomposition (steps 1–4 of the engine): no object of the size of $n_W$ is touched. The recursion (step 5) is bound by the output size: $t_b\propto$ entries $\propto S^{5.9}$ (580 at $S=4$, 52 090 at $S=10$, $3.06\times10^6$ at $S=20$); time per coefficient 15 µs at $S=4$, 4.4 µs at $S=10$, 1.5 µs at $S=20$ (G1 §6).
* **Exact and float.** Generators, Racah–Speiser, PRV spaces and the reduced operators are exact (Python fractions); the recursion and the eigen-decomposition of the multiplicity-free table are float64. The float64 recursion's HW residual stays at or below $1.3\times10^{-14}$ up to $S=12$ and grows by ×2.4–3.7 per unit of $S$ beyond $S=14$ ($2.6\times10^{-13}$ at $S=16$, $3.9\times10^{-11}$ at $S=20$, $7.5\times10^{-9}$ at $S=24$); the exact recursion rounded at the end keeps $\sim10^{-15}$ at 1.1–15.6 s per vector ($n_W=2.9\times10^4$ to $2\times10^5$, not vectorised) (G1 §4.4).
* **Machine difference.** i5-1345U (Ubuntu, g++ 11.4, LAPACK 3.10, one run) against i7-12700KF (Windows 11, Python, minimum of five cold runs); single-thread speeds differ by a factor of order 1–2 (G1 §6), so the ratios are orders of magnitude, not a measured speed-up. For gl(3) the same routine became 19 times faster when its linear algebra was repaired; no such run exists for SU(4), so every Munich number here is for the code as distributed.
* **Python-overhead regime.** For $S\le2$ the Munich C++ is faster in absolute terms (28 µs against 0.7 ms at $S=1$; 0.86 ms against 1.8 ms at $S=2$); the crossover is at $S=3$ (10.7 ms against 4.0 ms).
* **Scaling.** Munich's total grows by ×5.2–7.0 per unit of $S$ between $S=4$ and $S=8$ (0.0574 → 0.355 → 2.49 → 15.3 → 80.2 s; cost $O(n_W^3)$ per irrep); the engine's $t_b$ grows by ×2.0 ($S=4\to5$) decreasing to ×1.4 ($S=11\to12$). Computed here from the table: against $t_a$ the ratio at $S=8$ is $80.19/0.0258\approx3100$; against $t_b$ 860.
* **Not run.** The literal comparison with the Munich HW vectors (the binary cannot be built on the Windows machine): the GT bases coincide (pattern order = Munich's $Q(M)$, 6024 orthonormal lowering coefficients equal exactly, trivial sign gauge), so it should be a formality up to one global sign; G1 §9 has the 19 `hwcgc.sh` commands for $S\le3$ and a tested comparer.

**Sources.** Table A and the noise table: `RESULTS-G1-slice-engine-gl4.md` §6 (block 1 of 3; per-irrep table and singlet cold/warm; breakdown at $S=8$); Munich column: `Paul/Cpp/AlexCGC/results/hw_timings.txt` (totals per $S$) and `README.md`; residual growth, exact-rounded path: G1 §4.4; Munich-basis identity and the 19 commands: G1 §9; PRV-basis numbers for the families: G1 §7. Table B, the breakdown and the cost remarks: `RESULTS-G3-slice-ops-gl4.md` §9 (all rows; the $\ast$ rows are the extension rows). Table C: `Paul/Python/gaudin_sl3/RESULTS-slice-engine.md` §5, `Slice-engine-notes.md` Table 2, `HANDOFF.md` ("DONE (2026-10-21 night, Ubuntu): ndsu3lib benchmark"). Only derived number: the ratio 3100 (computed here).
