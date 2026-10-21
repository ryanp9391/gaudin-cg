# RESULTS-bethe-series-scale (agent S, round 6, WP-S, 2026-10-21): the eps-series of the Laurent data at scale — order 0 from the Q-system pencil, Frobenius recursion on series, numeric degenerate solver

Files (all mine): `bethe_series_scale.py` (module + CLI; `python -X utf8 bethe_series_scale.py --stage order0|validate|family|depth|timing --lam l1 l2 --N N --dps D --phis A|B|AB`),
`SPEC-bethe-series-scale.md` (the brief, written first; an "as built" addendum at its end), caches `bethe_series_scale_{420,630,840}.pkl` (every branch: x_0, the series as 100-digit strings, rank, singular values, lags,
times; gitignored, rebuilt in 0.3 / 4 / 10 minutes) and `bethe_series_scale_1470_samples.pkl` (summary of the [14,7,0] samples), logs `bethe_series_scale_*.log`. Python only (mpmath, sympy for exact polynomial gcd/factor in the pencil, no Wolfram).
Imports read-only: `nu_general` (`_parse`, `W4DIR`, `load_w4`), `bethe_series` is NOT needed any more (its Q(sqrt D) solver is replaced), `gl3_reps` (GT weights), `general_kernel.load_rep`, W5's `ref/laurent-series-{210,310,320}.txt`,
W4's `ref/untwisted-exact-*.txt`, agent R's `ref/cg-slices-{420,630,840,1470}.txt`, agent Q0's `ref/qsystem-spectrum-{630,840,1470}.txt` (cross-checks of order 0 only), `master_kernel2_rep_600.pkl`, agent K's `kernel_Theta_{840,1470}.npz` (exact kernels, read through `kernel_L14.theta_T_from_file`; used only for the depth analysis). No chain data and no external order-0 data enter any series; no git.

## 0. Headline numbers

* **Order-0 source used for EVERY series in this report: my own pencil (no files read, no Gröbner)** — the untwisted Q-system is a tridiagonal pencil in b (§1), derived independently of the shipped `pencilSolve` and **not a port of it**: it takes the exact GCD of **all three** echelon consistency polynomials, and a system without unknowns
  (d_1 = mu3 = 0, the analogue of M1 = 0) enters through its equation itself (a 1 x 1 augmented matrix), so it contains both fixes of the coordinator's correction and does not inherit the shipped bug (explicit test: [14,7,0]^2, mu = (16,13,13): my polynomial is `b`, one root b = 0, V = a_1 - eta/2 = 1431 — not {1383, 1431, 1479}).
  One unknown per branch, b = b_1^0; counts of roots = Racah–Speiser multiplicities for **all 113 irreps / 344 copies of [14,7,0]** (and 210, 310, 320, 410, 420, 630, 840);
  agrees with W4's exact eta/alpha/beta on all 40 branches of 210, 310, 320, 410 (difference 0 / 2.5e-60 at 60 digits), with agent R's V_red eigenvalues on **all 492 branches** of 420, 630, 840, 1470 (<= 4.6e-27, the 30-digit data of R) and with agent Q0's CORRECTED spectrum `ref/qsystem-spectrum-{630,840,1470}.txt`
  (V, H1, a1, b1, H2, a2, b2 of every copy: all 44 + 85 + 344 branches, <= 4.6e-27 = their 30 digits); the eight equations vanish on x_0 (<= 1e-55, 60 digits, up to [8,4,0]; 1e-50 at [14,7,0]). Cost: 0.5 s / 3 s / 12 s / 238 s for ALL irreps of 420 / 630 / 840 / 1470 incl. the checks (exact rational determinants, 1–4.5 s per irrep at a = 7).
* **Acceptance 1 (exact agreement with W5)**: [2,1,0], [3,1,0], [3,2,0], all **26 branches x 2 directions**, orders 0..10, all six data of both sites: **max abs deviation 4.4e-99 / 8.1e-97 / 3.4e-97, max rel. deviation <= 2.5e-90** (100 digits), no failures; 0.01–0.33 s per branch.
* **[4,2,0]**: all 19 branches (multiplicity up to 3) x 2 directions to order 14: 0 failures, unique surviving candidate in all 38 runs, residual of the eight equations (recomputed from scratch) **<= 7e-94**, 17 s for the whole family.
* **[6,3,0]**: all 44 branches x 2 directions to order 20: 0 failures, unique survivors, residual **<= 1.2e-88**, 127 s of branch time (136 s wall) for the family (mean 1.7 s per rank-5 run, max 3.5 s, 0.08 s per rank-6 run). The two chain-free tests of §4 (leading polynomial with the exact kernel) also pass for all of them.
* **Stress: [8,4,0]** all 85 branches (phi_A, order 26): 0 failures, residual <= 2.0e-85, 346 s of branch time (365 s wall; mean 4.5 s per rank-5 branch, max 10.2 s). **[14,7,0] samples** (the near-singlet, the multiplicity-8 irrep (21,14,7) all 8 copies, the multiplicity-3 irrep (26,14,2), (28,13,1), the top state), orders up to N_B + 2 = 44:
  residuals <= 4e-97, 67–83 s for the top state to order 44 (dps 120 / dps 80, different machine load). **Tolerance lesson (§3): at (28,13,1) of [14,7,0] a near-solution with relative inconsistency ~1e-41 survived the old 1e-40 pruning tolerance to order 80 (two "survivors", one rational); with the tolerance tied to the precision (default 10^-(0.6 dps), 1e-72 at 120 digits) the physical rational solution is the unique survivor.**
* **Depth (measured, N1's method: nu-series x exact kernel -> leading polynomial)**: for **all 174 branches** of 210, 310, 320, 420, 630, 840 (kernels: exact; 420 and 630 from the block theorem on the exact one-row kernels, 840 from agent K's exact `kernel_Theta_840.npz`) **and 14 branches of [14,7,0]** (K's exact `kernel_Theta_1470.npz`: the near-singlet, all 8 copies of the multiplicity-8 irrep (21,14,7), the 3 copies of (26,14,2), (28,13,1), the top state) the first non-vanishing order of sum G nu nu equals **N_B = 4 lam1 - 2 lam2 - k_B** (lower orders cancel to <= 9e-94 relative at 100 digits, <= 2e-91 at 120 digits for [14,7,0]) and the depth (largest eps-order of nu^(n) entering the leading polynomial) equals **N_B for every branch with k_B >= 1 and N_B/2 = 2 lam1 - lam2 for the top state (k_B = 0)**: top states 3 / 5 / 4 / 6 / 9 / 12 / **21** for 210 / 310 / 320 / 420 / 630 / 840 / **1470**; deepest measured branch of [14,7,0]: (28,13,1), depth 41 (N_B = 41).
* **Time**: per rank-5 branch at order 10: 0.15 / 0.28 / 0.28 / 0.50 / 1.3 / 2.5 s for 210 / 310 / 320 / 420 / 630 / 840 (rank-6 branches 0.01–0.07 s); top state to order 40: 1.2 / 1.7 / 1.6 / 3.8 / 9.1 / 17.8 s and **67 s at [14,7,0]** (dps 80); about 0.1–0.6 s per further order on the top state.

## 1. Order 0 (eps = 0): closed forms + the echelon Baxter pencil (WP item 1)

**Closed forms** (units theta12 = 1; U1/W4, checked here on every irrep up to [14,7,0] through the agreement with R's V_red): H_1 = eta = |Lambda|^2 - (C2(mu) - 2 C2(Lambda))/2, H_2 = -H_1, b_2 = -b_1,
a_1 + a_2 + b_1 = 0, (a_1 - a_2)/2 = q = (mu1+2)(mu2+1) mu3/2. Hence x_0(b) = (eta, q - b/2, b, -eta, -q - b/2, -b): one unknown b per branch. (Derivation: at infinity P ~ u^k with k^3 - (3+2|L|)k^2 + (2+2|L|+2c2+eta)k - 2q = 0 and
the three exponents are (mu3, mu2+1, mu1+2) [monodromy at 0 and 1 trivial => at infinity trivial]; both closed forms are Vieta's relations, and the Casimir form of eta is identically the e_2 form for lam3 = 0.)

**The oper as a three-term recursion.** Multiplying the untwisted dual oper by u(u-1) (site 1 at u = 1, site 2 at u = 0):
u^2(u-1)^2 P''' - |L|(2u-1)u(u-1) P'' + [eta u(u-1) + c2(2u^2-2u+1)] P' + [b/2 - q(2u-1)] P = 0, i.e. for P = sum p_j u^j
**alpha_{N-1} p_{N-1} + beta_N p_N + rho(N+1) p_{N+1} = 0**, alpha_j = j(j-1)(j-2) - 2|L|j(j-1) + (eta+2c2)j - 2q (zero exactly at the three exponents), beta_j = -2j(j-1)(j-2) + 3|L|j(j-1) - (eta+2c2)j + q + b/2,
rho_j = j(j-lam2-1)(j-lam1-2): **the spectral parameter b/2 enters linearly on the diagonal** (checked: adjoint octet gives b^2 = 45, [3,1,0] (5,2,1): 105, (4,3,1): 192, [3,2,0]: 105, 192).
**Echelon pencil** (own implementation of the idea of the .wl `pencilSolve` — the Baxter equations are affine in the eigenvalue and the spectrum is a GCD of consistency polynomials — but with the three systems treated uniformly): the three polynomial solutions have degrees d_1 < d_2 < d_3 = (mu3, mu2+1, mu1+2); P_j monic of degree d_j without u^{d_i} terms (i < j);
the d_j + 1 equations (N = 0..d_j) are affine in b in d_j - (j-1) unknown coefficients; b is admissible iff [A|v] has a null vector with last entry 1; polynomial in b = gcd of exact determinants of three random integer projections of [A|v] (a system with NO unknowns, j = 1 with mu3 = 0, is the 1 x 1 matrix of its equation, i.e. the equation itself is the consistency polynomial);
the branch values = the squarefree part of the gcd over **all three** j (exact Fractions + sympy gcd, `nroots`). **Pitfalls found on the way:** (a) the gcd of the three *plain* characteristic polynomials of the tridiagonal matrices T_{d_j}(b) is useless — since alpha_{d_1} = 0 the continuant of T_{d_1+1} contains that of T_{d_1} as a factor, so every root of det T_{d_1} survives
(e.g. b in {0, +-12} for the adjoint singlet instead of {0}); the echelon normalisation removes exactly these; (b) **the gcd over only two of the echelon systems is still not enough**, the same failure type as the coordinator's report on the shipped code: gcd(1,2) has more roots than the multiplicity on 3 of 25 irreps of [6,3,0] ((8,8,2), (8,7,3), (7,7,4)) and on 6 of 41 irreps of [8,4,0]
and on **21 of the 113 irreps of [14,7,0]** (e.g. (20,20,2), (20,19,3), ..., (20,15,7)); gcd(1,3) fails on 7 / 14 / 45 irreps, gcd(2,3) on 9 / 16 / 49, each single system on 13–15 / 24–26 / 73–77 (of 25 / 41 / 113); the three-way gcd has exactly the multiplicity on every irrep ([6,3,0], [8,4,0], [14,7,0]; `bethe_series_scale_pencil_systems_1470.log`). Roots are real; copies are ordered by a_1 = q - b/2 ascending (copy 1 = smaller a_1, W4/W5/Q0 convention).

| Lambda | irreps / branches | root count = Racah–Speiser mult. | vs W4 (exact) | vs agent R V_red eigenvalue (30 digits) | vs agent Q0 (corrected spectrum; V, H1, a1, b1, H2, a2, b2; 30 digits) | max eight-equation residual at eps = 0 (60 digits) | time, all irreps (incl. checks) |
|---|---|---|---|---|---|---|---|
| [2,1,0] | 5 / 6 | all | 0 | – | – | 5e-60 | 0.1 s |
| [3,1,0] | 8 / 10 | all | 0 | – | – | 2e-59 | 0.2 s |
| [3,2,0] | 8 / 10 | all | 0 | – | – | 5e-59 | 0.2 s |
| [4,1,0] | 11 / 14 | all | 2.5e-60 | – | – | 1e-58 | 0.4 s |
| [4,2,0] | 13 / 19 | all | – | 4.0e-29 | (no file) | 3e-58 | 0.5 s |
| [6,3,0] | 25 / 44 | all | – | 4.1e-28 | 4.6e-28 | 1e-56 | 2.9 s |
| [8,4,0] | 41 / 85 | all | – | 4.8e-28 | 4.8e-28 | 6e-55 | 12 s |
| [14,7,0] | 113 / 344 | all (max mult 8) | – | 4.6e-27 | 4.6e-27 | 1e-50 | 238 s |

(a_1^0 - eta/2 equals R's V_red eigenvalue, the relation documented in R's header. Racah–Speiser reproduces the plan's 113 irreps / 344 copies at a = 7.) Sector (mu3 = 0) is always multiplicity 1 with b = 0.
Not needed here, but the same pencil gives the Q-polynomials (null vectors) if round 7 wants them.

## 2. The series engine (WP items 2–3)

`Engine(lam, mu, phi, th12)`: the Frobenius recursion of N1 at theta_s (exponents {0, lam2+1, lam1+2}, solutions u, v, regular parts of S1, S2, S3 from the other site and sigma_1, e_2, e_3 of zeta = eps*phi) is run on **truncated eps-series, one eps-order per pass**. Per recursion step m (m <= lam1+2, two solutions, two sites):
five series products (H_s(m-1) - a_s) c_{m-1}, b_s c_{m-2}, H_t W_m, b_t X_m, a_t Y_m (W, X, Y: fixed scalar-weighted sums of the lower c's, the (-1)^k delta^(-k-1) weights) plus shifted twist terms; cost per pass O(lam1 * n) multiplications, i.e. **O(lam1^2) series coefficients per order** (no sympy polynomial anywhere).
Arrays are stored per order, so `truncate(n)` + re-pass and `clone()` (list slicing) are cheap. Outputs per order: the eight equations F_n (A, B1, B2 of site 1, A, B1, B2 of site 2, the two global rules), and `nu_series()` (nu_kl = k! l! (u_k v_l - u_l v_k)/(lam2+1)! as series, two series products per pair).
The Jacobian J = dF/dx at (x_0, 0) is J[:, j] = F_1(x_1 = e_j) - F_1(0) (six order-1 passes).
Validation: (i) **the eight equations vanish on W5's exact series for all orders <= 10, all 52 (branch, direction) runs of 210/310/320** (max |F_n| 1.1e-78 / 2.2e-78 / 8.6e-78 at 80 digits, theta12 = 31/35); (ii) nu-series vs N1's polynomial formulas (`nu_numeric`) on the same series at eps = 1e-3 / 1e-2:
relative difference 3.8e-34 / 3.8e-25 ([2,1,0]), 1.3e-28 / 1.3e-20 ([3,1,0]), 3.2e-30 / 3.2e-22 ([3,2,0]) — pure truncation (order 11); the sharper test is §4 (all lower orders of sum G nu nu cancel); (iii) the singular values of J reproduce S1's rank table (rank 5 on the multiple roots, rank 6 elsewhere), see §3.
Pass cost: 0.5 ms (lam1 = 3, n <= 10) ... 10 ms (lam1 = 14, n ~ 30).

## 3. The numeric degenerate solver (WP item 2)

`BranchSolver`: rank of J from a pivoted rref with threshold 1e-35 and the singular values (reported); rank 6 -> x_n unique, the two extra equations are checked consistent at every order. **Rank 5** (kernel dimension 1; rank <= 4 raises `NotImplementedError` — it never occurred in any of the 277 runs):
x_n = x_n^p + c_n k. **Stage m** (x_0..x_m fixed): the pending coordinate c_{m+1} is fixed by the first later order n > m+1 at which the solvability residual g_n = N r_n (the three left-null-vector combinations) is not identically zero ("lag" = n - m - 1). *Scan*: three trial runs
(c = 0.37, -1.21, 2.43) on cloned engines advance in lockstep until g_n is non-zero (this costs 2 passes per run per order); g_n independent of the later pending coordinates is verified (two extra runs); the polynomial g_n(c) is reconstructed from the three scan points plus extra trial points with adaptive degree (verified on one more point);
the **real common roots of the non-zero components** are the candidates (roots of the lowest-degree component by `polyroots`, clusters of nearby roots -> centroid, Newton polish on the component of largest slope *until converged*, test against all components), then refined by Newton steps on g_n itself (the fit is spread over c ~ O(1..5), where g is huge, and loses up to 40 digits of the root at deep lag / high degree;
the direct residual at the root does not). Candidates are followed depth first; a wrong one dies at the next stage (no common root of the next-stage components, or the no-pending-unknown consistency test g_{m+1} = 0 fails). Complex roots are discarded (the physical branches are real; S1 used K-rational roots only).
Working precision 100 (family runs) / 120 ([14,7,0]) digits, zero tests at 10^-tol_exp of the natural scale of the residual, **tol_exp = max(40, 0.6 dps)** (60 at 100 digits, 72 at 120; see the [14,7,0] trap below). Defaults look = maxdeg = 2 lam1 + 4.

**Structure found (all 178 rank-5 runs of 420/630/840): the lag is the same at every stage; the degree of g_n in c_{m+1} is lag + 1 at the first stage and 1 at every later stage.**
* **Top state** (mu = 2 Lambda): lag = 2 min(lam2, lam1 - lam2): 2 ([2,1,0], [3,1,0], [3,2,0], [4,1,0], [4,3,0], [5,1,0]), 4 ([4,2,0], [5,2,0], [5,3,0], [6,2,0]), 6 ([6,3,0]), 8 ([8,4,0]), **14 ([14,7,0], degree 15 polynomial, 4 real candidates of which 3 die at the next stage)**. S1's "lag 2" was the small-rep value.
* On the line mu2 = 2 lam2 of the reps with lam1 = 2 lam2 the lag is lam1 - mu3 and there is exactly one real common root at c_1 ((12,6,0), (11,6,1), ..., (7,6,5) of [6,3,0]: lags 6, 5, 4, 3, 2, 1; (8,4,0), (7,4,1), (6,4,2), (5,4,3) of [4,2,0]: 4, 3, 2, 1).
* All other rank-5 branches have lag 1..4 ([6,3,0]) and up to lam1 + 1 real candidates at the first stage; **in 178 runs of 420/630/840 the consistency tests alone left exactly one survivor** (candidate deaths: 48 in 12 of 38 runs (420), 294 in 44 of 88 (630), 516 in 52 of 85 (840)); no multiple roots among the candidates.
* **Pruning tolerance at [14,7,0] (a trap found at (28,13,1), k_B = 1, multiplicity 1)**: 13 real candidates at c_1 (lag 12, degrees [13, 13, -1]; after refinement all three components of g_13 vanish at every one of the 13 roots to 1e-121, the two non-zero components are proportional, stage 2 does not discriminate, the candidates die at later stages); with the zero-test tolerance 1e-40 **two candidates survived every consistency test up to order 80** (N = 24, 43, 80; 33 deaths): c_1 = -83.39393939... = -2752/33 (rational at phi_A; eight-equation residual 2.7e-109) and c_1 = -74.10044338964941... (irrational; residual **1.6e-41**, i.e. a near-solution whose inconsistency at some stage lies between 1e-60 and 1e-40 of the natural scale). With the tolerance tied to the working precision (1e-60 or 1e-80 at 120 digits: 36 deaths) exactly one survivor remains, the rational one, and it is the physical one: its series satisfies sum G nu nu = eps^{N_B}(...) with agent K's exact [14,7,0] kernel to 4.7e-96 (§5). The default is now `tol_exp = max(40, 0.6 dps)`; with it nothing changed for 210..840 (the same 178 unique survivors, the same candidate deaths). `BranchSolver.series` keeps a safety net for genuine ties — PSLQ membership of c_1 in K_b = Q(1, b, ..., b^(d-1)) (d = degree of b over Q, `field_degrees`; the physical solution is Galois stable because the chain is defined over Q) — which never fires with the tight tolerance. Consistency alone therefore determined the series in every branch tested, with a tolerance that must be much tighter than the 1e-40 that sufficed up to [8,4,0].
* Jacobian ranks and gaps: rank-5 branches have sigma_6/sigma_1 <= 1.4e-98 and sigma_5/sigma_1 >= 1.3e-3 (420), 5.0e-5 (630), 3.4e-6 (840); rank-6 branches sigma_6/sigma_1 >= 5e-4 / 1.3e-5 / 4.9e-7. Ranks reproduce S1's table for 210/310/320; rank 5 on 14 of 19 branches of 420 (5 rank 6), 37 of 44 of 630 (7 rank 6), 76 of 85 of 840 (9 rank 6) — all multiplicity >= 2 irreps are rank 5.

### Acceptance 1 (W5's exact series)

| Lambda | (branch, direction) runs | failures | max abs deviation | max rel. deviation (|ref| > 1e-12) | time (all runs, 100 digits) |
|---|---|---|---|---|---|
| [2,1,0] | 12 | 0 | 4.4e-99 | 4.8e-92 | 1 s |
| [3,1,0] | 20 | 0 | 8.1e-97 | 5.5e-91 | 3 s |
| [3,2,0] | 20 | 0 | 3.4e-97 | 2.5e-90 | 3 s |

(units: my series are in theta12 = 1 and rescaled by theta12^(n-1) (H, a), theta12^(n-2) (b) to W5's theta12 = 31/35.) phi-homogeneity x_n(2 phi) = 2^n x_n(phi): <= 1.9e-74 (80 digits) on the 4 tested branches of 420 and the adjoint.

## 4. [4,2,0] and [6,3,0] (and [8,4,0]): all branches

Per-irrep data (phi_A and phi_B, N = 14 resp. 20, 100 digits; "real common roots" = candidates at the first stage; times are the maximum over the copies; residual = max over both directions of the eight equations recomputed from scratch on the final series).

**[4,2,0]** (13 irreps, 19 branches, 38 runs, 0 failures, 17 s):

| irrep mu | mult | k_B | N_B | rank | first-stage lag / degree | real common roots at c_1 | candidates that die later | max time per copy (s) | max residual (both phi) |
|---|---|---|---|---|---|---|---|---|---|
| (8, 4, 0) | 1 | 0 | 12 | 5 | 4 / 5 | 1 | 0 | 0.98 | 7.1e-94 |
| (8, 3, 1) | 1 | 1 | 11 | 5 | 2 / 3 | 3 | 12 | 0.70 | 1.1e-95 |
| (8, 2, 2) | 1 | 2 | 10 | 6 | - | - | 0 | 0.04 | 1.3e-98 |
| (7, 5, 0) | 1 | 1 | 11 | 5 | 2 / 3 | 3 | 12 | 0.69 | 3.7e-95 |
| (7, 4, 1) | 2 | 2 | 10 | 5 | 3 / 4 | 1 | 0 | 0.71 | 3.6e-95 |
| (7, 3, 2) | 2 | 3 | 9 | 5 | 1 / 2 | 2 | 12 | 0.31 | 1.2e-97 |
| (6, 6, 0) | 1 | 2 | 10 | 6 | - | - | 0 | 0.03 | 9.8e-98 |
| (6, 5, 1) | 2 | 3 | 9 | 5 | 1 / 2 | 2 | 12 | 0.29 | 3.8e-97 |
| (6, 4, 2) | 3 | 4 | 8 | 5 | 2 / 3 | 1 | 0 | 0.51 | 1.1e-96 |
| (6, 3, 3) | 1 | 5 | 7 | 6 | - | - | 0 | 0.03 | 2.6e-97 |
| (5, 5, 2) | 1 | 5 | 7 | 6 | - | - | 0 | 0.03 | 2.0e-97 |
| (5, 4, 3) | 2 | 6 | 6 | 5 | 1 / 2 | 1 | 0 | 0.25 | 5.0e-97 |
| (4, 4, 4) | 1 | 8 | 4 | 6 | - | - | 0 | 0.03 | 3.2e-97 |

**[6,3,0]** (25 irreps, 44 branches, 88 runs, 0 failures, 136 s wall):

| irrep mu | mult | k_B | N_B | rank | lag / degree | real common roots at c_1 | dying later | max time (s) | max residual (both phi) |
|---|---|---|---|---|---|---|---|---|---|
| (12, 6, 0) | 1 | 0 | 18 | 5 | 6 / 7 | 1 | 0 | 3.54 | 2.7e-90 |
| (12, 5, 1) | 1 | 1 | 17 | 5 | 4 / 5 | 5 | 24 | 3.28 | 2.9e-91 |
| (12, 4, 2) | 1 | 2 | 16 | 5 | 2 / 3 | 3 | 24 | 1.90 | 7.9e-94 |
| (12, 3, 3) | 1 | 3 | 15 | 6 | - | - | 0 | 0.09 | 4.1e-96 |
| (11, 7, 0) | 1 | 1 | 17 | 5 | 4 / 5 | 3, 5 | 18 | 3.35 | 3.1e-91 |
| (11, 6, 1) | 2 | 2 | 16 | 5 | 5 / 6 | 1 | 0 | 2.82 | 2.0e-90 |
| (11, 5, 2) | 2 | 3 | 15 | 5 | 3 / 4 | 4 | 36 | 2.36 | 7.5e-93 |
| (11, 4, 3) | 2 | 4 | 14 | 5 | 1 / 2 | 2 | 24 | 0.81 | 2.5e-96 |
| (10, 8, 0) | 1 | 2 | 16 | 5 | 2 / 3 | 3 | 24 | 2.00 | 9.5e-93 |
| (10, 7, 1) | 2 | 3 | 15 | 5 | 3 / 4 | 2, 4 | 24 | 2.46 | 6.3e-93 |
| (10, 6, 2) | 3 | 4 | 14 | 5 | 4 / 5 | 1 | 0 | 2.41 | 3.5e-91 |
| (10, 5, 3) | 3 | 5 | 13 | 5 | 2 / 3 | 3 | 36 | 1.58 | 9.8e-92 |
| (10, 4, 4) | 1 | 6 | 12 | 6 | - | - | 0 | 0.08 | 2.5e-96 |
| (9, 9, 0) | 1 | 3 | 15 | 6 | - | - | 0 | 0.08 | 2.2e-96 |
| (9, 8, 1) | 2 | 4 | 14 | 5 | 1 / 2 | 2 | 24 | 0.82 | 1.8e-95 |
| (9, 7, 2) | 3 | 5 | 13 | 5 | 2 / 3 | 3 | 36 | 1.61 | 5.5e-93 |
| (9, 6, 3) | 4 | 6 | 12 | 5 | 3 / 4 | 1 | 0 | 1.78 | 1.2e-88 |
| (9, 5, 4) | 2 | 7 | 11 | 5 | 1 / 2 | 2 | 12 | 0.72 | 2.5e-95 |
| (8, 8, 2) | 1 | 6 | 12 | 6 | - | - | 0 | 0.08 | 2.5e-96 |
| (8, 7, 3) | 2 | 7 | 11 | 5 | 1 / 2 | 2 | 12 | 0.72 | 1.2e-95 |
| (8, 6, 4) | 3 | 8 | 10 | 5 | 2 / 3 | 1 | 0 | 1.45 | 4.8e-92 |
| (8, 5, 5) | 1 | 9 | 9 | 6 | - | - | 0 | 0.09 | 4.1e-95 |
| (7, 7, 4) | 1 | 9 | 9 | 6 | - | - | 0 | 0.08 | 2.3e-95 |
| (7, 6, 5) | 2 | 10 | 8 | 5 | 1 / 2 | 1 | 0 | 0.68 | 6.7e-95 |
| (6, 6, 6) | 1 | 12 | 6 | 6 | - | - | 0 | 0.09 | 7.3e-95 |

([8,4,0]: 41 irreps, 85 branches (phi_A, N = 26, 100 digits): 76 rank 5 / 9 rank 6, lags up to 8, 516 candidate deaths in 52 runs, residual <= 2.0e-85, 4.5 s mean per rank-5 branch, 10.2 s max, 365 s wall.)

**Independent check of the new series (no chain, no fit): the identity sum G nu nu = eps^{N_B}(...) with the EXACT kernel.** The kernel of [4,2,0] is the block theorem Theta^{[l1,l2,0]} = [l2!(l1-l2)!/l1!]^2 (a a')^(-l2) Theta^{[l1,0,0]} applied to `general_kernel_rep_400.pkl`, that of [6,3,0] the same
applied to `master_kernel2_rep_600.pkl` (exact [6,0,0] kernel of M3); my construction reproduces N1's `G_terms` exactly for 210/310/320 and the block theorem reproduces the exact [2,1,0], [3,1,0], [3,2,0] kernels exactly. With the nu-series of the new Laurent series at theta12 = 31/35,
for **all 19 + 44 branches the orders t < N_B of sum_{d+n+m=t} G^[d]_ij nu_i^(n) nu_j^(m) cancel to <= 8e-98 (420) / <= 9e-94 (630) of the size of the individual terms and the first non-vanishing order is exactly N_B = 4 lam1 - 2 lam2 - k_B** (k_B = 2 lam1 - (mu1 - mu3)): a cancellation over up to 17 orders that any error
in x_n (n <= N_B) would destroy. Together with the eight-equation residuals, the unique survivors, the homogeneity test and the order-0 check against R, this validates the [4,2,0] and [6,3,0] series without chain data.

## 5. The depth actually needed (measured)

Method (N1's): leading polynomial R_{N_B}(w, w') = sum_{d+n+m=N_B} G^[d]_ij nu_i^(n)(theta_1) nu_j^(m)(theta_2); **depth = the largest n (or m) of a contribution with |term| > 1e-25 of the largest** at t = N_B. phi_A, 100 digits.

| Lambda | branches | N_B measured = rule | depth = N_B | depth of the top state | lower-order cancellation (relative) |
|---|---|---|---|---|---|
| [2,1,0] | 6 | 6/6 | 5 of 5 | 3 (N_B = 6) | <= 1.2e-101 |
| [3,1,0] | 10 | 10/10 | 9 of 9 | 5 (N_B = 10) | <= 9e-100 |
| [3,2,0] | 10 | 10/10 | 9 of 9 | 4 (N_B = 8) | <= 2e-99 |
| [4,2,0] | 19 | 19/19 | 18 of 18 | 6 (N_B = 12) | <= 8e-98 |
| [6,3,0] | 44 | 44/44 | 43 of 43 | 9 (N_B = 18) | <= 9e-94 |
| [8,4,0] | 85 | 85/85 | 84 of 84 | 12 (N_B = 24) | <= 4e-86 |
| [14,7,0] (14 sample branches) | 14 | 14/14 | 13 of 13 | 21 (N_B = 42) | <= 2e-91 |

(N1's table for 310/320 and U2's depths for the adjoint (3, 5, 5, 4, 4, 2) are reproduced exactly; the deepest branches: [3,1,0] (6,1,1), (5,3,0): 9; [4,2,0] (8,3,1), (7,5,0): 11; [6,3,0] (12,5,1) ... (11,7,0): 17.)
**Law (measured on 188 branches, 6 reps up to [14,7,0]): depth = N_B = 4 lam1 - 2 lam2 - k_B for every branch except the top state, whose depth is N_B/2 = 2 lam1 - lam2.** For the benchmark [14,7,0]: N_B = 42 - k_B (41 at k_B = 1 is the deepest, e.g. (28,13,1) measured), the top state needs only 21; the near-singlet (14,14,14) (k_B = 28) needs 14, the multiplicity-8 irrep (21,14,7) (k_B = 14) needs 28 for all 8 copies (so the plan's "~40" is right for the near-top irreps only).
([8,4,0]: all 85 branches with `kernel_Theta_840.npz`, `bethe_series_scale_depth_840.log`; [14,7,0]: 14 branches, `bethe_series_scale_depth_1470.log`; per-branch lines (N_B rule/measured, lower-order residue, depth) are in the logs `bethe_series_scale_depth_*.log`.)
The solver reaches these depths (§6: orders 40–44 at [14,7,0] in 45–83 s per branch).

## 6. Time per branch vs order and vs lambda

Time per branch at fixed N = 10 (80 digits, phi_A, all branches, includes the Jacobian; the pencil is not included — it is 0.01–0.2 s per irrep up to [6,3,0], 1–4.5 s at [14,7,0]):

| Lambda | branches (rank 5 / rank 6) | mean | max | rank-5 mean | rank-6 mean | top state | family total |
|---|---|---|---|---|---|---|---|
| [2,1,0] | 6 (3 / 3) | 0.08 s | 0.21 s | 0.15 s | 0.012 s | 0.21 s | 0.5 s |
| [3,1,0] | 10 (6 / 4) | 0.18 s | 0.49 s | 0.28 s | 0.025 s | 0.49 s | 1.8 s |
| [3,2,0] | 10 (6 / 4) | 0.18 s | 0.48 s | 0.28 s | 0.023 s | 0.48 s | 1.8 s |
| [4,2,0] | 19 (14 / 5) | 0.38 s | 0.95 s | 0.50 s | 0.035 s | 0.95 s | 7.1 s |
| [6,3,0] | 44 (37 / 7) | 1.08 s | 3.1 s | 1.27 s | 0.060 s | 2.3 s | 47 s |
| [8,4,0] | 85 (76 / 9) | 2.27 s | 8.8 s | 2.53 s | 0.069 s | 4.8 s | 193 s |

Cumulative time to fix order n (80 digits; "top" = mu = 2 Lambda, "oct" = the irrep of largest multiplicity (copy 1), "rank 6" = the singlet-like irrep):

| Lambda | irrep | n = 5 | 10 | 20 | 30 | 40 | (44) |
|---|---|---|---|---|---|---|---|
| [2,1,0] | top (4,2,0) | 0.08 | 0.18 | 0.44 | 0.77 | 1.18 | |
| [2,1,0] | oct (3,2,1) | 0.04 | 0.09 | 0.21 | 0.37 | 0.57 | |
| [3,1,0] | top (6,2,0) | 0.12 | 0.26 | 0.64 | 1.14 | 1.74 | |
| [3,2,0] | top (6,4,0) | 0.12 | 0.26 | 0.61 | 1.05 | 1.59 | |
| [4,2,0] | top (8,4,0) | 0.29 | 0.61 | 1.42 | 2.47 | 3.84 | |
| [4,2,0] | oct (6,4,2) mult 3 | 0.15 | 0.36 | 0.84 | 1.51 | 2.27 | |
| [6,3,0] | top (12,6,0) | 0.71 | 1.44 | 3.38 | 5.86 | 9.09 | |
| [6,3,0] | oct (9,6,3) mult 4 | 0.36 | 0.77 | 1.94 | 3.37 | 5.00 | |
| [8,4,0] | top (16,8,0) | 1.47 | 3.12 | 7.06 | 11.9 | 17.8 | |
| [8,4,0] | oct (12,8,4) | 0.95 | 1.85 | 4.17 | 6.73 | 10.0 | |
| [14,7,0] | top (28,14,0) | 8.79 | 14.4 | 28.2 | 44.9 | 66.8 | 82.9 |
| [14,7,0] | oct (21,14,7) mult 8 | 4.49 | 8.82 | 18.5 | 30.3 | 40.5 | 45.1 |
| [14,7,0] | rank 6 (14,14,14) | 0.05 | 0.11 | 0.24 | 0.45 | 0.66 | 0.8 |

Rank-6 irreps cost 0.06–0.26 s to order 40 for all lam (80 passes). The cost per further order grows mildly with n (series products are O(n)), the first-stage cost is dominated by the polynomial reconstruction (top state at [14,7,0]: 8.8 s to order 5). Time scales roughly like lam1^2 (top state to order 40: 1.2 -> 17.8 s from lam1 = 2 to 8).
[14,7,0] samples at their N_B + 2 (120 digits, phi_A, residuals of the eight equations recomputed from scratch; the unique survivor was found by consistency alone in all of them, no tie-break): near-singlet (14,14,14) N = 16: 0.2 s, 6.7e-109; (21,14,7) mult 8, N = 30: 19–26 s per copy (178 s for all 8), <= 3.8e-97; (26,14,2) mult 3, N = 40: 49–50 s per copy (152 s), <= 3.1e-109; (28,13,1), N = 43: 85 s, 2.7e-109; top (28,14,0), N = 44: 67 s, 1.3e-110.
A full [14,7,0] family would be ~344 branches x 20–190 s (embarrassingly parallel over branches and directions: a few hours on one core, ~15–20 minutes on the 20 cores here); not run.

## Surprises

1. The eps = 0 problem is a tridiagonal eigenproblem in the single number b/2 (alpha, beta, rho above); Gröbner bases are not needed at all, and the pencil also gives the Q-polynomials. The naive "gcd of the three characteristic polynomials" does NOT work (inherited factors); the echelon systems are essential, and all THREE of them are needed (gcd of two over-counts on 3/25, 6/41 irreps of [6,3,0], [8,4,0]) — the same trap as the shipped `pencilSolve` (coordinator's correction), which the corrected data of Q0 and agent R confirm my pencil avoids on all 344 copies of [14,7,0].
2. The degenerate lag of the top-type branches grows with lambda: lag = 2 min(lam2, lam1 - lam2) (14 at [14,7,0]) with a polynomial of degree lag + 1 for c_1 — S1's "lag 2, cubic" is the small-rep case. Everything else (lag constant over the stages, degree 1 afterwards) is regular enough to automate.
3. With many real candidates (13 at (28,13,1) of [14,7,0]) the first-stage polynomial system is exactly degenerate (proportional components), so ALL roots survive stage 1 and the wrong ones die only through tiny inconsistencies at later stages (here ~1e-41 relative for the worst ghost): a pruning tolerance of 1e-40 produced two "formal solutions" to order 80, 10^-(0.6 dps) produces one. In 178 runs up to [8,4,0] and in the [14,7,0] samples the consistency tests alone (with the tight tolerance) isolated the physical solution, which is rational on multiplicity-1 irreps (Galois stable); no extra equation or selection rule was needed.
4. Numerical pitfalls hit and fixed: (a) importing `nu_general`/`general_kernel` resets `mp.mp.dps` — my first "100-digit" validation silently ran at 80 digits (deviations 1e-81); every entry point now re-applies the precision; (b) root polishing: Newton from a polyroots root that is not a common root of all components can converge slowly onto the common root and block
   the accurate duplicate (c_1 correct to 1e-50 only, residual 4e-52 at the (8,4,0) top state, phi_B) — now iterate to convergence and keep the more accurate of coinciding roots; (c) the fitted polynomial loses up to 40 digits of the root at lag 7+ (degree 8 at (21,14,7)): Newton steps on g_n itself restore 1e-106.
5. Every run (52 validation runs, 38 + 88 + 85 family runs, 14 samples at [14,7,0]) had rank 5 or 6 of the 8 x 6 Jacobian: a kernel of dimension 2 never occurred. Which multiplicity-1 irreps are rank 5 has no simple rule in k_B ((5,5,0) of [3,2,0] with k_B = 1 is rank 6, (8,2,2) of [4,2,0] with k_B = 2 is rank 6, (12,4,2) of [6,3,0] with k_B = 2 is rank 5).

## Not done / open

* The depth of [14,7,0] was measured on 14 of the 344 branches only (5 irreps), because each branch needs its own series (20–140 s) and a 45–80 s analysis; the other 330 are covered by the law, not by measurement. The solver reached order 26 on all of [8,4,0] and 44 on 5 irreps of [14,7,0].
* The full [14,7,0] family (344 branches) was not run; the pencil for all of [14,7,0] is the slowest part of the order-0 route (267 s, exact rational determinants) and could be replaced by modular/numeric root finding.
* Uniqueness of the formal series through x_0 is established numerically (tolerance-dependent: the ghost of (28,13,1) had an inconsistency of 1e-41), not proved; the PSLQ field-membership safety net (heuristic: height limit, only c_1) was never needed with the tight tolerance. A Bethe-side proof that the consistency conditions single out the physical solution is open.
* Rank <= 4 (kernel dimension >= 2) is not handled (raises); not met. Complex candidate roots are discarded by assumption (physical data are real).
* The nu-series and the series were never compared with a chain at [4,2,0] or larger directly (only through the exact identity of §4 and the independent checks listed there); no X_B / dictionary / CG slices (round 7).

## Reproduce (all from `Paul/Python/gaudin_sl3`, `python -X utf8`)

* Acceptance 1: `bethe_series_scale.py --stage validate --N 10 --dps 100 --phis AB` (3 s per rep).
* Order 0: `--stage order0 --lam 6 3 --dps 60` (also 4 2 / 8 4 / 14 7; compares with W4, R, Q0 whichever files exist; `14 7` takes 4 min).
* Families: `--stage family --lam 4 2 --N 14 --phis AB --dps 100` (17 s); `--lam 6 3 --N 20` (2.3 min); `--lam 8 4 --N 26 --phis A` (6 min); cache `bethe_series_scale_<l1><l2>0.pkl`, per-branch log lines.
* Depth: `--stage depth --lam 6 3 --dps 100 --phis A` (kernel found automatically: `general_kernel_rep`, block theorem, or agent K's `kernel_Theta_<l1><l2>0.npz`; 8 4 takes ~3 min, 14 7 is meant for single irreps via `kernel_summary` + `leading_analysis_fast`).
* Time vs order: `--stage timing --lam 8 4 --N 40 --dps 80`.
* From Python: `BranchSolver(lam, mu, x0, phi).series(N)` with `x0 = order0_branches(lam, mu)[copy-1]`; `solve_sector_branches(lam, mu, phi, N)` does all copies of an irrep (order 0, rank, series, field degrees); `Engine.nu_series()` gives the nu_kl(theta_s) series;
  units theta12 = 1, convert with `to_theta_mp` (H, a ~ theta12^(n-1), b ~ theta12^(n-2)).
