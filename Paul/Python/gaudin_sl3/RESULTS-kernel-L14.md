# RESULTS-kernel-L14 (agent K, round 6 WP-K, Python only) -- COMPLETE: the exact one-row master kernel through lambda_1 = 14; Theta^{[8,4,0]} and Theta^{[14,7,0]} delivered, exact and as float64

Files: `kernel_L14.py` (library + CLI, ~1500 lines), `SPEC-kernel-L14.md` (brief, revised at the end to the route actually used), data (all `.npz`, rebuildable, 0.48 GB in all):
`kernel_N_L8.npz`, `kernel_N_L14.npz` (N as float64 CSR + the exact integer M), `kernel_Theta_840.npz`, `kernel_Theta_1470.npz` (float64 coefficients + exact digits), `kernel_Theta_420.npz` ([4,2,0] test file),
`kernel_M_L{1..7}.npz` (small caches of the exact M; M^(8) and M^(14) live inside the N files, L = 9..13 are rebuilt on demand by `get_M`). Environment: `python -X utf8`, numpy 2.5, scipy 1.18 (OpenBLAS, threads capped at 8 inside the module),
sympy/mpmath only in the one-site exact part and in checks. No Wolfram, no git, existing scripts untouched.
**Warning: `.gitignore` ignores `Paul/Python/gaudin_sl3/*.pkl` only; `*.npz` is NOT ignored (kernel_N_L14.npz is 0.42 GB, kernel_Theta_1470.npz 48 MB) -- add a pattern before any `git add -A`.**
Regeneration (all from `Paul/Python/gaudin_sl3`): `python -X utf8 kernel_L14.py M <L> [batch]` | `N <L> [nprimes|- [nofloat] [nosave]]` | `Theta <L> <l2> [nprimes|- [nosave]]` | `univ <L...>` | `accept [Lbig]`.

## 0 Conventions
* The M2/M3 one-row model: V = C[x1,x2,x3]_L, E_ij = x_i d_j, jet pairs (k,l), 0 <= k < l <= L+1 (`pairs_L(L)` == `general_kernel.pairs_of((L,0,0))`, lexicographic, d = C(L+2,2) = 120 at L = 14),
  two-site row index r = i*d + j (i = site 1 = the `a` site, j = site 2 = the `a'` site), jet level rho(k,l) = k+l-1, two-site rho_r = rho(r1)+rho(r2); t = 1/theta12, a = 1/w1, a' = 1/w1'.
* M := M(1) = Y(1) (y x y)^-1 (integer, unit block-lower-triangular in rho), N := M^-1. `Theta^{[L,0,0]}_{pq}(a,a',t) = sum_{r,s} c_r c_s a^rho_r a'^rho_s t^{rho_r+rho_s-rho_p-rho_q} N_{rs;pq}`, c_kl = `ckl_closed`;
  `Theta^{[L,l2,0]}` = [l2!(L-l2)!/L!]^2 (a a')^-l2 x the entries of `Theta^{[L,0,0]}` on the non-forced pairs (l >= l2+1), which by label-triangularity needs only N[S x S, S x S] = (M[S x S, S x S])^-1 (8464 x 8464 for [14,7,0]).
* Computation point: zeta = 0 (g = E21 + E32), th1 - th2 = 1 (delta = +1 at site 1, -1 at site 2), monomial basis. By the grading nothing is lost: M(t)_{rp} = t^{rho_r-rho_p} M(1)_{rp} was verified in float64 at th12 = 2, -3, 1/2 for L = 2, 3, 4
  (absolute deviations 2e-16 .. 3.4e-9 on entries up to 4e3, i.e. rounding).
* File formats. `kernel_M_L<L>.npz`: L, pairs (n_pairs,2) int16, rho (n,) int16, M_{indptr,indices,data,shape} (CSR, int64, exact). `kernel_N_L<L>.npz`: the same arrays for M plus N_{...} (CSR float64 = the exact integers
  correctly lifted and rounded to double, original row order; N[r,c] != 0 => rho_r >= rho_c) and scalars depth, nprimes, maxabsN, build_seconds, float_solve_max_rel_err, float_solve_seconds, restriction6_max_rel_dev; loader `load_N(path)` -> dict with csr `N`, `M`.
  `kernel_Theta_<L><l2>0.npz`: lam, pairs (nS,2) int16 (== `general_kernel.pairs_of(lam)`, asserted), entry_ptr (nS^2+1), n/p/q int16, coef float64 with
  Theta_ij = sum_{m in entry_ptr[e]:entry_ptr[e+1]} coef[m] t^n[m] a^p[m] a'^q[m], e = i*nS + j (T-format of `general_kernel`; block-theorem prefactor and the (a a')^-l2 shift already included), prefactor, and the
  EXACT coefficients: `exact_digits` (depth, nmon) int32 = symmetric Garner digits of the integer C~_m = sum_k digits[k,m] prod_{j<k} primes[j], exact coefficient = C~_m * exact_num / exact_den (decimal strings in the file).
  Loaders: `load_Theta(path)` (adds the dict `T[(i,j)] = {(n,p,q): float}`), `theta_coefficients_exact(arrs, slice)` (Fractions), `theta_T_from_file(path)` (exact T-format), `eval_Theta_entry_mp(arrs, e, a, a', t, dps)` (mpmath with the exact coefficients), `G_eps_coefficients(arrs, e, wh, wh', t, dps)` (the exact eps-polynomial of G_ij(eps wh, eps wh'; th12), agrees with the direct evaluation to 1e-60),
  `N_exact_ints(M, rho, cols)` (exact Python-int columns of N on demand, any L).

## 1 Headline results
1. **The operator route scaled to L = 14, exactly.** The covector Frobenius recursion (never an expanded polynomial) in exact modular arithmetic reproduces the integer matrix M^(L) for every L = 2..14. For L = 2..6 it equals the exact M3 data
   (`master_kernel2_N6.pkl` inverted) entry by entry, each L computed independently and restricted to labels <= L+1; L = 7..14 are new (L = 14: 14400 x 14400, 13 414 280 nonzeros = 6.5 %, max |M| = 2 529 843 369 384 038 400 = 2.5e18 < 2^62).

   | L | dim | nnz(M) | fill | max abs M | build (4 primes) | peak WS |
   |---|---|---|---|---|---|---|
   | 2 | 36 | 157 | 12.1 % | 4 | | |
   | 3 | 100 | 949 | 9.5 % | 40 | | |
   | 4 | 225 | 4295 | 8.5 % | 540 | | |
   | 5 | 441 | 15098 | 7.8 % | 12864 | | |
   | 6 | 784 | 46025 | 7.5 % | 241200 | 0.3 s | |
   | 7 | 1296 | 121199 | 7.2 % | 6.16e6 | 1.2 s | |
   | 8 | 2025 | 288316 | 7.0 % | 2.16e8 | 4.0 s | 0.46 GB |
   | 9 | 3025 | 629178 | 6.9 % | 8.39e9 | 9.9 s | |
   | 10 | 4356 | 1286626 | 6.8 % | 3.87e11 | 28 s | 1.8 GB (one batch) |
   | 11 | 6084 | 2470294 | 6.7 % | 1.51e13 | 62 s | 3.4 GB (one batch) |
   | 12 | 8281 | 4529290 | 6.6 % | 7.49e14 | 118 s | 6.2 GB (one batch) |
   | 13 | 11025 | 7939935 | 6.5 % | 4.31e16 | 222 s | 2.4 GB (batches of 20) |
   | 14 | 14400 | 13414280 | 6.5 % | 2.53e18 | 392 s (first build 458 s under load) | 3.6 GB (batches of 20) |

   (The structure is asserted, not assumed: `check_unit_triangular` finds M_{rc} != 0 only for rho_r > rho_c or r = c with M_rr = 1, i.e. the level-preserving part of Y(1) equals y (x) y exactly, and the maximal level difference is exactly L; asserted for L = 1..12 and 14, L = 13 only through the restriction property.)
   ("build" = 4 full passes: 3 primes for the CRT lift and a 4th prime that verifies EVERY entry; wall clock on the shared 20-thread machine. M(t) is a polynomial of degree L in t: nnz of the level-difference-n parts M_n at L = 14, n = 0..14:
   14400, 518056, 874602, 1177848, 1401279, 1531414, 1564762, 1507862, 1369938, 1169328, 927648, 671478, 421263, 207286, 57116; max|M_n| grows from 468 (n = 1) to 2.5e18 (n = 11, 12) and falls to 6.6e16 (n = 14).)
2. **Universality of M holds exactly through L = 14** (M3 had L <= 6): M^(L)|_{labels <= L'+1} == M^(L') as integer matrices for all 78 pairs 2 <= L' < L <= 14 (`kernel_L14.py univ 2 3 ... 14`). The inverse of a block-triangular matrix restricts the same way, so
   acceptance (ii) holds: N^(L)|_{labels <= 7} == `master_kernel2_N6.pkl` with relative deviation 0.0 and no spurious nonzero, for L = 8, 9, 10, 11, 12, 14 (integer-for-integer for L = 8 with exact Python ints: 140275 nonzeros, 0 mismatches).
3. **N exactly; N^(14) in 5 minutes.** Exact block forward substitution modulo primes p < 790000 with float64 BLAS and a Garner lift:

   | L | primes needed (depth) / used | nnz(N) | fill | max abs N | time | peak WS | plain float64 forward substitution from the exact M: max rel error of N |
   |---|---|---|---|---|---|---|---|
   | 6 | 4 | 140275 | 22.8 % | 2.56e17 | 0.2 s | | 0 |
   | 7 | 4 | 389164 | 23.2 % | 5.79e22 | 0.6 s | | 7.1e-15 |
   | 8 | 5 / 10 | 956691 | 23.3 % | 2.58e28 | 1.8 s | 0.38 GB | 5.9e-12 |
   | 9 | 6 / 11 | 2140433 | 23.4 % | 2.05e34 | 5.3 s | 0.64 GB | 4.7e-10 |
   | 10 | 7 / 12 | 4438905 | 23.4 % | 2.71e40 | 13.8 s | 1.07 GB | 6.8e-9 |
   | 11 | 8 / 13 | 8647251 | 23.4 % | 5.58e46 | 33.7 s | 1.84 GB | 3.2e-7 |
   | 12 | 10 / 14 | 15985769 | 23.3 % | 1.72e53 | 74.1 s | 3.13 GB | 2.0e-5 |
   | 14 | 12 / 16 | 48074739 | 23.2 % | 4.7503e66 | 302 s | 8.80 GB | 1.1e-5 |

   (L = 13 not run.) log10 max|N| = 17.4, 22.8, 28.4, 34.3, 40.4, 46.7, 53.2, 66.7 for L = 6..12, 14: the increments grow from 5.4 to 6.7 per unit of L.
   Exactness certificate independent of the modular arithmetic: the exact integer columns of N (Garner digits -> Python ints) satisfy M N[:,c] = e_c exactly in Python-int arithmetic: 40 random columns at L = 8 and 24 random columns at L = 14 (max |N| in those columns 3.99e63), 0 failures.
4. **Acceptance (i): the exact kernels are reproduced EXACTLY** (every monomial coefficient as an exact Fraction, identical monomial support, 0 mismatching entries): Theta^{[L,0,0]}, L = 1..6 (rep_300, rep_400 from `general_kernel`; rep_500, rep_600 from `master_kernel2`), and the block-theorem kernels
   [1,1,0] [2,1,0] [2,2,0] [3,1,0] [3,2,0] [3,3,0] [4,4,0]; built both from the independently computed M^(L') and from M^(14) restricted to labels <= L'+1 (`kernel_L14.py accept 14`). The float64 coefficients of the same kernels deviate from the exact ones by <= 2.7e-16 relative (required: 1e-10).
5. **The new kernels** (`Theta^{[8,4,0]}`, `Theta^{[14,7,0]}` delivered; [14,0,0] and [4,2,0] also computed):

   | kernel | pairs | entries | monomials | max abs coef | min abs coef | a,a'-degree | t-degree | build / peak WS |
   |---|---|---|---|---|---|---|---|---|
   | [4,2,0] | 12 | 144 | 2721 | 2.02e6 | 2.8e-2 | 6 | 11 | 0.0 s |
   | [8,4,0] | 35 | 1225 | 84273 | 1.92e22 | 2.0e-4 | 12 | 23 | 1.2 s / 0.21 GB |
   | [14,7,0] | 92 | 8464 | 1691086 | 1.2966e51 | 8.49e-8 | 21 | 41 | 122-134 s / 3.25 GB |
   | [14,0,0] (not saved) | 120 | 14400 | 4025449 | 4.7503e66 | 1.0 | 28 | 52 | 586 s / 5.8 GB |

   The degree laws of RESULTS-general-kernel (a-degree 2 l1 - l2, t-degree 4(l1-1) - max(0, 2 l2 - 3)) hold at l1 = 8 and 14 (predictions 23 and 41). Monomials per entry of [14,7,0]: min 1, median 194, max 482; per-entry max/min coefficient ratio: median 3.9e28, max 9.9e49.
6. **Acceptance (iii): exchange rule and top-state sum rule, EXACT.**
   * Exchange rule Theta_ij(a,a',t) = Theta_ji(a',a,-t): 0 violations in all residues of all primes (5.7e6 residues for [8,4,0], 1.57e8 for [14,7,0] with 22 primes, 2.9e8 for [14,0,0] with 24 primes); the stored float64 coefficients satisfy it bitwise (max relative deviation 0.0, identical supports).
   * Top-state sum rule sum_ij Theta_ij m_i m_j = M M R(a) R(a'): evaluated in F_p at random residues for every prime used in the CRT (an exact identity of rational functions; 22 primes x 2 trials for [14,7,0], 16 x 2 for [8,4,0], 24 x 2 for [14,0,0]): 0 failures. G1's own exact checker
     (`general_kernel.top_sum_rule_check`: sympy rationals at random rational points, applied to the exact coefficients reconstructed from the files): worst deviation exactly 0 for [4,2,0] (3 points), [8,4,0] (3 points) and [14,7,0] (2 points). The [14,0,0] check exercises all of N^(14), including the rows/columns with l < 8 that [14,7,0] never uses.
   * Negative control ([8,4,0], 8 primes): shifting ONE coefficient C~ by +1 at one monomial is caught by the exchange rule (16 violations) and by the sum rule (16 of 16 (prime, trial) checks fail), so the residue checks are not vacuous.
   * Independent closed forms of G1 (not used in the construction): the first-order term `theta1_closed` (conjecture C4, fitted on l1 <= 4) agrees with the t^1 part of the new kernels to 7.4e-15 ([14,7,0], 600-800 random entries at random rational (a,a')) and 1.7e-14 ([8,4,0]); the corner block (pairs (k,l1+1)) equals the single-index closed form
     (C3: (-1)^{k+k'} C(l,k) C(l,k') ... F_{mm'}) times [l2!(l1-l2)!/l1!]^2 (a a')^{l1-l2} MONOMIAL BY MONOMIAL: 225 entries / 13994 monomials at [14,7,0] (identical support, max relative deviation 4.4e-16) and 81 / 1889 at [8,4,0] (2.9e-16); G1 had verified C3 only for l <= 5.
   * Kernels with NO stored exact counterpart that were generated and pass the exact exchange rule and the exact (residue) sum rule, 12 primes x 2 trials, 0 failures: [4,2,0], [5,1,0], [6,3,0] (block theorem on the new N), [6,6,0]; the stored ones ([2,1,0], [3,1,0], [3,2,0], [4,0,0], [6,0,0], ...) pass the same checks.
7. **Acceptance (iv): conditioning, accuracy, and what double precision can and cannot do.**
   * *The pure float64 route for M fails from L = 8.* Same code with float64 instead of residues:

     | L | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
     |---|---|---|---|---|---|---|---|---|---|
     | max abs error of the integer M | 8.9e-16 | 2.4e-13 | 5.3e-11 | 1.5e-8 | 4.7e-6 | 1.5e-2 | 4.8 | 1.3e4 | 1.7e7 |
     | relative to the row maximum | 2e-16 | 1e-14 | 2e-13 | 7e-12 | 4e-10 | 4e-8 | 4e-6 | 4e-4 | 6e-2 |
     | max abs Y(1) (= L!) / max abs M | 2 / 4 | 6 / 40 | 24 / 540 | 120 / 1.3e4 | 720 / 2.4e5 | 5040 / 6.2e6 | 4.0e4 / 2.2e8 | 3.6e5 / 8.4e9 | 3.6e6 / 3.9e11 |

     about 2.5 digits are lost per unit of L; the integers cannot be recovered by rounding from L = 8. Y(1) in the monomial basis is small (max = L!), but M = Y (y x y)^-1 is up to 1e5 larger and the recursion's error grows by 300-3000 per unit of L.
   * *Plain float64 forward substitution for N, even from the exact integer M*: max relative error 5.9e-12 (L = 8), 6.8e-9 (10), 2.0e-5 (12), 1.1e-5 (14, worst of 48 million entries). The stored float64 pair (M, N) has componentwise relative residual max|(M N - I)_ij| / (|M||N|)_ij = 2.5e-16 (6.9e-17 at L = 8) but absolute residual
     max|M N - I| = 9e52 (2.4e13 at L = 8) because (|M||N|)_ij reaches ~4e68: N has to be obtained exactly (our lift is correct to <= 1.2e-15 relative per entry), and M N = I cannot be reproduced in floating point.
   * *The coefficients of Theta are accurate to double precision.* Stored float64 vs the correctly rounded exact rational: max relative difference 6.8e-16 over all 1 691 086 coefficients of [14,7,0] (4.3e-16 over 84273 of [8,4,0]); 51 % / 67 % are bitwise equal to the correctly rounded value. A plain float64 contraction of the stored float64 N reproduces
     the exact [14,7,0] coefficients to median 1.4e-16, 99 % <= 7e-15, max 4.1e-11 (no spurious monomials): the contraction is benign, the difficulty is N.
   * *Evaluating one entry in plain double precision* (numpy power sum of the stored coefficients vs mpmath with the exact coefficients; 30-50 random entries per point; WLOG t = 1 by the homogeneity Theta_ij(la,la',t/la) = la^{e_i+e_j} Theta_ij(a,a',t); grid a = s, a' = +-0.7 s):
     relative error <= 2e-15 (kappa = sum|terms|/|value| <= 24) for a a' < 0 at every s = 0.01 .. 1000 and for a a' > 0 at s = 0.01 and s >= 1; for a a' > 0 and s = 0.03 .. 0.3 the entries pass near their zeros: kappa up to 4.7e8 (median 5.5e5 at s = 0.1), relative error up to 1.5e-8.
     So the polynomial COEFFICIENTS can be trusted to <= 8e-16 everywhere; a double-precision EVALUATION of an entry is trustworthy (<= 1e-15) except near the zeros (a a' > 0, |a| t ~ 0.03-0.3), where up to 8 digits are lost.
   * *Contractions over the jet indices are different.* The top-state sum rule sum_ij Theta_ij m_i m_j cancels by 1e11..3e26 for [8,4,0] (6 random points with |a|,|a'| <= 1, |z|,|theta| ~ 1) and by 2.5e28..2.9e29 for [14,7,0] (2 points): in mpmath with the stored double coefficients the identity holds with relative deviation
     (0.007..0.03) * 1.1e-16 * kappa, i.e. exactly the rounding of the coefficients, so a numerical check of such identities in double precision is meaningless (hence the exact residue checks above). **The eps-series of the untwisting limit,
     G_ij(eps w^, eps w^'; th12) = sum_m c_m t^n w^^(2l1-l2-p) w^'^(2l1-l2-q) eps^(2(2l1-l2)-p-q) (a polynomial in eps with exponents >= 0, 42 - p - q >= 0 at [14,7,0]), must use the exact coefficients (the stored digits; `eval_Theta_entry_mp`, `theta_coefficients_exact`)
     in high-precision arithmetic wherever order-by-order cancellations occur (the sum over (i,j) at order eps^{N_B}); the double-precision coefficients suffice only for orders without cancellation.**

## 2 Time and memory (wall clock on the shared machine; peak = process working set)

| step | L = 8 | L = 14 |
|---|---|---|
| exact M (4 primes, batches) | 4.0 s, 0.46 GB | 392 s, 3.6 GB (batches of 20 site-1 rows) |
| exact N, all columns (10 / 16 primes) | 1.8 s, 0.38 GB | 302 s, 8.8 GB (chunks of 1000 columns; int32 level blocks for 16 primes 2.7 GB + residues/digits of one chunk 3.5 GB + the 1.7 GB result) |
| float64 N (plain forward substitution) | 0.0 s | 6.4 s (1.7 GB dense) |
| Theta [8,4,0] / [14,7,0] with all checks | 1.2 s, 0.21 GB | 122-134 s, 3.25 GB (22 primes x 8464 columns of N[S x S]) |
| Theta [14,0,0] (not saved) | | 586 s, 5.8 GB |
| exact certificate (Python-int check of N columns) | 1.6 s (40 columns) | 270 s (24 columns) |

Scaling: the M build grows by a factor 1.8-2.8 per unit of L (4.0 s at L = 8, 392 s at L = 14: ~L^7..L^8), the exact N by 2.0-2.9 per unit of L (1.8 s -> 302 s); N dominates the storage (23 % fill of a lower block-triangle, 48 million exactly-lifted entries at L = 14).
The first L = 14 build of M (458 s, 3.61 GB) ran while other tests of mine were active; the quiet rebuild (392 s, 3.61 GB) is bitwise identical (same CSR arrays). The 20-thread BLAS default had to be avoided (see Surprises).

## 3 Usage (round 7)
```python
import sys; sys.path.insert(0, 'Paul/Python/gaudin_sl3')
import kernel_L14 as K
th = K.load_Theta('Paul/Python/gaudin_sl3/kernel_Theta_1470.npz')   # th['pairs'] = non-forced (k,l), th['T'][(i,j)] = {(n,p,q): float}, Theta_ij = sum c t^n a^p a'^q
nS = len(th['pairs']); e = i * nS + j
val = K.eval_Theta_entry_mp(th, e, a, ap, t, dps=80)                  # exact coefficients, mpmath
# G_ij(eps*wh, eps*wh'; theta12): monomial (n,p,q) -> c * t^n * wh^(21-p) * wh'^(21-q) * eps^(42-p-q)   (2*l1 - l2 = 21 at [14,7,0])
g = K.G_eps_coefficients(th, e, wh, whp, t, dps=80)                  # the exact eps-polynomial of G_ij: list g[k], k = 0 .. 42 (exact coefficients, mpmath)
fr = K.theta_coefficients_exact(th, slice(th['entry_ptr'][e], th['entry_ptr'][e + 1]))   # Fractions of entry e
d = K.load_N('Paul/Python/gaudin_sl3/kernel_N_L14.npz')              # d['N'] (csr float64), d['M'] (csr int64 exact), d['pairs'], d['rho']
Nint = K.N_exact_ints(d['M'], K.rho_vector(14), [col0, col1])        # exact Python-int columns of N (70 s for 2 columns at L = 14, dominated by the one-off preparation of the 16 prime blocks; agrees with the stored float64 N to 2.9e-16)
```

## 4 Surprises
* M is an INTEGER matrix, sparse (6.5 % at L = 14), whose entries stay below 2^62 through L = 14, although N = M^-1 reaches 4.75e66 (the inverse of a unit block-triangular integer matrix of depth L).
* The pure floating-point operator route is hopeless beyond L = 7, while the exact modular route costs only 4x one float pass (392 s for the whole M^(14) including the verification prime). Even the final triangular solve loses all accuracy in floating point beyond L ~ 10; the Theta coefficients computed from the exact N are benign.
* The entire kernel is determined by 14400 x 14400 *integers*; the exact coefficients of Theta^{[14,7,0]} are rationals with denominator dividing 3432^2 * 360360^2 (3432 = C(14,7)); the largest denominator in a sample of 2000 coefficients is 226512.
* Universality, the degree laws, the first-order closed form C4 and the corner closed form C3 -- all conjectured from data with l1 <= 5 -- hold exactly (universality, degree laws) or to 1e-14 (C3, C4) at l1 = 8 and l1 = 14, so the structure found in round 5 is real.
* OpenBLAS with 20 threads is ~100x slower than with 8 for small-m dgemm (80x1000 @ 1000x1000: 1.1 GFlop/s vs 250-300 for the transposed product); fixed by the transposed layout of the solver and a thread cap of 8 (set in the module before numpy is imported).
* The top-state sum rule is an extremely ill-conditioned identity (condition 1e29 at [14,7,0]); the 'numerical precision' check requested in the brief had to be done in exact residue arithmetic (done); a double-precision check only confirms agreement to ~0.01-0.03 u kappa.
* The brief says `.npz`/`.pkl` are gitignored; only `*.pkl` is.

## 5 Not done / notes for round 7
* The kernel for other lambda_2 at l1 = 14 is one CLI call away (`python -X utf8 kernel_L14.py Theta 14 <l2>`, 2-10 minutes each; [14,0,0] was computed for the checks but not saved: 4e6 monomials), as is exact N on demand (`N_exact_ints`).
* No fit of the loop part of the universal operator exp(sum t^n D_n) (Research 2); the exact N^(14) and the exact M of every L <= 14 are the data for it (M is the better object: sparse, integer, polynomial of degree L in t; Theta^{[14,0,0]} contains the top entry of N as its largest coefficient 4.7503e66); the D_n beyond n = 5 were not touched.
* Not tested: another twist/theta point in exact arithmetic (only the float64 grading check at th12 = 2, -3, 1/2), the L = 14 float64 Y-route (known to fail), L = 13 for N, Theta for l1 > 14 (S0 = lcm(1..15) clears the denominators l <= 15 of c_kl and would have to be extended to lcm(1..l1+1)).
