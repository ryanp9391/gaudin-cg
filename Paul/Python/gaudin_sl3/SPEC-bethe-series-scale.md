# Spec S (round 6, WP-S): the eps-series of the Laurent data at scale — numeric (mpmath) degenerate Fredholm solver, series-arithmetic Frobenius recursion, order 0 from the untwisted Q-system (2026-10-21)

Written by agent S (SONNET, Python only, no Wolfram kernel) as its own brief, from `Paul/Mathematica/Gaudin/Gaudin-sl3-round6-plan.md` §0–§3 (row S), §5, §6 and the sources
listed in §1. Deliverables: this file, `bethe_series_scale.py`, `RESULTS-bethe-series-scale.md` (style of the existing `RESULTS-*.md`). Run from
`Paul/Python/gaudin_sl3` with `python -X utf8`. mpmath 1.3 (pure-python backend, ~0.9 us per 80-digit multiply), sympy 1.14, numpy; no gmpy2, no Wolfram. Existing
scripts are read-only (import them). No git commits. <= 2 heavy processes.

## 0. Goal

Round 5 (agent S1, `bethe_series.py`) produced the eps-series of the six Laurent data x = (H_1, a_1, b_1, H_2, a_2, b_2) of every highest-weight-limit branch of
two sites with Lambda = [lam1, lam2, 0] **from the Bethe side alone** (six no-log conditions + two global rules), but (i) in exact arithmetic over Q(sqrt D) (impossible
for the copies of a multiplicity-m sector: degree-m algebraic numbers), (ii) with the order-0 points from Groebner bases (did not finish for [3,2,0]), (iii) with the
no-log conditions as expanded sympy polynomials in 10 variables (`nu_general.frobenius`), which does not scale with lam1. The Q-route of round 7 needs the series
(and the nu-series built from them) for [4,2,0], [6,3,0], [8,4,0], [14,7,0]. Task: make the solver work for general Lambda at scale:

1. **Order 0 without Groebner**: H^0 = eta/theta12 (closed form), a^0 and b^0 per copy from the untwisted polynomial Q-system (own port of the linear-pencil idea of
   `Minimal_complete_working-fast.wl::pencilSolve`; Q0's `ref/qsystem-spectrum-*.txt` are used only as a cross-check if they exist).
2. **Numbers instead of Q(sqrt D)**: mpmath at 60-100 digits; keep the degenerate Fredholm logic (rank-5 Jacobians, lag 1 and lag 2) but with numerical zero tests,
   real candidate roots of the solvability polynomials, depth-first following of all candidates and pruning by consistency at the next order.
3. **Series arithmetic**: the no-log conditions and the nu_kl are evaluated by running the Frobenius recursion of N1 (`nu_general.frobenius`) directly on truncated
   power series in eps (cost per order O(lam1^2) series coefficients), not through expanded polynomials.
4. **Measure**: time per branch vs order and vs lam; the depth actually needed.

## 1. Read first (all read-only)

`RESULTS-bethe-series.md` + `bethe_series.py` (S1: the equations, the ranks, the degenerate solver `Solver`, `solve_series`, the mechanism); `RESULTS-nu-general.md` +
`nu_general.py` (N1: `frobenius(lam)`, `nu_symbolic`, the depth table §2 (vi), `G_terms`/`leading_poly`); `Paul/Mathematica/Gaudin/ref/laurent-series-{210,310,320}.txt`
(W5: exact series, both directions phiA = (1/3,-2/7,5/11), phiB = (-2/5,3/8,7/11), orders 0..10, theta12 = 31/35, header documents format and branch labelling:
copy 1 = smaller a_1 at eps = 0); `RESULTS-untwisted-Q.md` (U1: closed forms of the untwisted data, polynomial solutions, exponents at infinity); W4's
`ref/untwisted-exact-{210,310,320,410}.txt` (BRANCH records: exact eta, alpha, beta per branch, theta12 = 1 after rescaling H = eta/th12, a = alpha/th12, b = beta/th12^2);
`Paul/Mathematica/Gaudin/Experiments/twist-pert-Baxter.md` (a Q-side perturbative solver, for orientation only); plan §3 row S.

## 2. The equations (restated; units theta1 - theta2 = 1, zeta = eps*phi)

Dual oper (d^3 - S1 d^2 + S2 d - S3) P = 0, x = u - theta_s, S1 = |L|/x + sum s_k x^k, S2 = c2/x^2 + H/x + sum h_k x^k, S3 = a/x^2 + b/x + sum t_k x^k,
|L| = lam1 + lam2, c2 = lam2 (lam1 + 1); regular parts from the other site t (delta = theta_s - theta_t = +-1) and the twist constants:
s_k = |L| (-1)^k delta^(-k-1) + sigma1 [k=0], h_k = H_t (-1)^k delta^(-k-1) + c2 (-1)^k (k+1) delta^(-k-2) + e2 [k=0], t_k = b_t (-1)^k delta^(-k-1) + a_t (-1)^k (k+1) delta^(-k-2) + e3 [k=0]
with sigma1 = eps*s1(phi), e2 = eps^2 e2(phi), e3 = eps^3 e3(phi). Frobenius P = sum c_m x^m, exponents {0, lam2+1, lam1+2}, rho(m) = m (m-lam2-1)(m-lam1-2),
rho(m) c_m = sum_k s_k (m-1-k)(m-2-k) c_{m-1-k} - (H (m-1) - a) c_{m-1} + b c_{m-2} - sum_k h_k (m-2-k) c_{m-2-k} + sum_k t_k c_{m-3-k}.
Solutions u (c_0 = 1, c_{lam2+1} = 0) and v (c_0 = 0, c_{lam2+1} = 1) for m <= lam1+1; nu_kl = k! l! (u_k v_l - u_l v_k)/(lam2+1)! on the non-forced pairs (k < l <= lam1+1, not both <= lam2).
No-log conditions per site: A = the equation at m = lam2+1 (solution u), B1, B2 = the equation at m = lam1+2 for u and v. Eight equations in the six unknowns:
(A,B1,B2) of site 1, (A,B1,B2) of site 2, H_1 + H_2 = eps sum_{i != j} phi_i mu_j, b_1 + b_2 = eps^2 (phi_1 phi_2 mu_3 + phi_1 phi_3 mu_2 + phi_2 phi_3 mu_1).
Order by order: J x_n = -r_n(x_0..x_{n-1}), J = dF/dx at (x_0, eps = 0), 8 x 6. Rank 6: unique. Rank 5 (multiple roots of the eps = 0 system: top-type and octet-type sectors):
x_n = x_n^p + c_n k, the 3 left-null-vector conditions N r_n = 0 are polynomial in the pending kernel coordinates; the pending c is fixed by the next order (lag 1: octet-like,
a quadratic at c_1, then linear) or the next-but-one order (lag 2: top-type, a cubic at c_1, then linear). S1 found a unique surviving K-rational candidate in all 15 cases.

## 3. Task S-0: the order-0 points (eps = 0)

Closed forms (U1/W4, verified in all sectors of five reps; derived below from the exponents at infinity): H_1 = eta/theta12, H_2 = -H_1, eta = |Lambda|^2 - (C2(mu) - 2 C2(Lambda))/2,
C2(l) = sum l_i (l_i + 4 - 2i); b_2 = -b_1; a_1 + a_2 + theta12 b_1 = 0; (a_1 - a_2)/2 = q/theta12, q = (mu1+2)(mu2+1) mu3/2. So in units theta12 = 1 each branch is
x_0(b) = (eta, q - b/2, b, -eta, -q - b/2, -b) with ONE unknown b = b_1 (the multiplicity-space eigenvalue; mult-1 sectors: b = 0; adjoint octet b = +-sqrt45).
(Derivation: at infinity P ~ u^k with k^3 - (3+2|L|) k^2 + (2 + 2|L| + 2 c2 + eta) k - 2q = 0 and the three exponents are (mu3, mu2+1, mu1+2); this reproduces both closed forms.)
**Pencil (the Q-system).** For P = sum p_j u^j the untwisted oper multiplied by u(u-1) reads u^2(u-1)^2 P''' - |L|(2u-1)u(u-1) P'' + [eta u(u-1) + c2(2u^2-2u+1)] P' + [b/2 - q(2u-1)] P = 0
(site 1 at u = 1, site 2 at u = 0), i.e. the three-term recursion alpha_{N-1} p_{N-1} + beta_N p_N + rho(N+1) p_{N+1} = 0 with alpha_j = j(j-1)(j-2) - 2|L| j(j-1) + (eta+2c2) j - 2q,
beta_j = -2j(j-1)(j-2) + 3|L| j(j-1) - (eta+2c2) j + b/2 + q: **linear in b/2**. The three polynomial solutions have degrees d = (mu3, mu2+1, mu1+2) (alpha_d = 0 for each). Echelon ansatz as in the
notebook (P_j monic of degree d_j, no u^{d_i} term for i < j): the equations N = 0..d_j are linear in the coefficients and affine in b; b is admissible iff the augmented matrix [A|v] has a
null vector with last entry 1; polynomial in b = gcd of two exact determinants of random integer projections (as `charPoly` in the .wl); the branch values are the common roots for j = 1, 2, 3 (exact
rational arithmetic, sympy Poly gcd), real. Number of roots must equal the multiplicity of mu in Lambda x Lambda (Racah-Speiser from GT weights). Check: reproduces W4's eta/alpha/beta for every branch of 210, 310, 320, 410.
Then the eight equations at eps = 0 are evaluated on x_0(b) (residual <= 1e-60) and J is analysed (rank, singular values).

## 4. Task S-1: the series engine

A truncated series is a list of mpf (eps^0..eps^N). One pass per ORDER n (incremental): given x_0..x_n (x_n provisional) compute the order-n coefficient of every c_m^{u,v} (m <= lam1+2 recursion
for the two solutions, both sites) and of the eight equations F_n; arrays are stored per order so that `truncate(n)` (drop orders >= n) followed by a re-pass is cheap. Per m the recursion
needs 5 series products (H_s(m-1) - a_s) c_{m-1}, b_s c_{m-2}, H_t W_m, b_t X_m, a_t Y_m with W, X, Y fixed scalar-weighted sums of lower c's (the (-1)^k delta^(-k-1) weights);
the twist constants add shifted terms. The Jacobian at eps = 0: J[:, j] = F_1(x_1 = e_j) - F_1(x_1 = 0). nu-series: nu_kl^{(n)} from the stored u, v series (two series products per pair).
Validation of the engine: (i) the eight residuals vanish on W5's exact series (all orders <= 10, both phi, all 26 branches); (ii) the nu-series equal N1's polynomial formulas (`nu_numeric`) at a numeric eps;
(iii) the Jacobian ranks of S1 (6 or 5 as tabulated) are reproduced.

## 5. Task S-2: the numeric degenerate solver

Rank decision by singular values / pivot tolerance (gap reported). Rank 6: x_n unique, the two extra equations consistent (checked at every order). Rank 5 (kernel dim 1; if rank < 5 is met the
code must say so and report, not silently continue): x_n = x_n^p + c_n k. Pending coordinate c_a: evaluate the solvability residual g_n(c_a) (the 3 left-null-vector combinations) at several trial values
of c_a (each = re-pass of the orders a..n), fit polynomials in c_a (degree detected adaptively, <= 5, verified at an extra point), take the real roots of the lowest-degree non-zero component that
are also roots (to tolerance) of the other components. g independent of later pending c's required; identically zero -> go one order further (lag 2). Follow every candidate (depth first), prune by the
consistency (no pending unknown) test at the next orders and by the final residual check; if several candidates survive to the end, break the tie with the order-0 data and the exact conditions
(and report that this happened); if none survives report the candidate roots and residuals — never fall back to chain data. Tolerances: zero tests at 10^-tol_exp of the natural scale, tol_exp = max(40, 0.6 dps) (the brief's fixed 1e-40 was found too loose at [14,7,0], see the addendum).

## 6. Acceptance (plan §3 row S)

1. Exact agreement (to the working precision, <= 1e-50 relative) with `ref/laurent-series-{210,310,320}.txt`: all 26 branches, both directions, orders 0..10, all six data of both sites (units rescaled by theta12 = 31/35).
2. All branches of [4,2,0] (multiplicity up to 3, new) and of [6,3,0]: order-0 data, ranks, series to the depth N_B (and beyond if cheap), unique surviving candidate everywhere (or an explicit failure list),
   residuals, eps^0 data cross-checked against the pencil solver (and against Q0's `ref/qsystem-spectrum-630.txt` if present).
3. Tables: time per branch vs order and vs lam (210, 310, 320, 420, 630); the depth needed per branch: measured (nu series + exact kernel `general_kernel.load_rep`, block theorem for [4,2,0] from rep_400)
   next to N1's rule N_B = 4 lam1 - 2 lam2 - k_B, k_B = 2 lam1 - (mu1 - mu3); report what was measured, not the rule.

## 7. Files

`bethe_series_scale.py` (module + CLI `--stage order0|validate|family|depth|timing`), `RESULTS-bethe-series-scale.md` (files, conventions, headline numbers, tables, surprises, not done),
caches `bethe_series_scale_*.pkl` (gitignored, rebuildable). Sources are never edited.

## 8. As built (addendum, written at the end; where reality differed from this brief)

* **Order 0** is my own echelon Baxter pencil in b (a tridiagonal three-term recursion, linear in b/2; the closed forms eta, q reduce each branch to ONE unknown) in exact rational arithmetic — no Groebner, no external order-0 data used in any series. It is NOT a port of the shipped `pencilSolve`: it takes the gcd over ALL THREE echelon systems and treats a system
  without unknowns through its equation itself, i.e. it contains the two fixes of the coordinator's correction (explicit test: [14,7,0]^2 mu = (16,13,13) gives the single root V = 1431). Checked against W4 (exact), agent R's V_red eigenvalues (`ref/cg-slices-*.txt`, all 492 branches of 420/630/840/1470), agent Q0's corrected `ref/qsystem-spectrum-{630,840,1470}.txt` (all branches, all seven numbers V, H1, a1, b1, H2, a2, b2) and the eight
  equations at eps = 0. (Plain gcd of the characteristic polynomials of the three tridiagonal matrices does not work: inherited factors; gcd of only two of the echelon systems over-counts on some irreps.)
* **Degenerate solver**: lags are much deeper than S1's 1–2 for large Lambda (top state: lag 2 min(lam2, lam1 - lam2), polynomial of degree lag + 1 at the first stage, degree 1 afterwards), so the candidate search is a lockstep scan over cloned engines with adaptive polynomial reconstruction and direct Newton refinement (§3 of the RESULTS).
  Consistency alone left a unique survivor in every branch tested (178 runs up to [8,4,0], the [14,7,0] samples) PROVIDED the zero-test tolerance is tight: at (28,13,1) of [14,7,0] (13 real candidates at the first stage) a near-solution with relative inconsistency ~1e-41 survived the brief's 1e-40 to order 80; the default is now 10^-(0.6 dps). A PSLQ field-membership test of c_1 in Q(b) (numeric analogue of S1's K-rationality) is kept as a safety net for ties and never fires. The brief's "use the order-0 data to pick" could not decide anyway: the candidates share x_0.
* **Depth**: measured with exact kernels for 210, 310, 320, 420, 630 (89 branches; 420: block theorem on `general_kernel_rep_400.pkl`; 630: block theorem on `master_kernel2_rep_600.pkl`), for all 85 branches of [8,4,0] and for 14 branches of [14,7,0] (agent K's exact `kernel_Theta_{840,1470}.npz`, available late in the round).
* Beyond the brief: [8,4,0] complete (85 branches), [14,7,0] samples to order 44, a phi-homogeneity test, and the nu-series (`Engine.nu_series`) validated through the leading-polynomial identity.
