# Spec K: the one-row master kernel at lambda_1 = 14 (and 8) by the operator route, scaled up (round 6, WP-K; 2026-10-21)

Delegated to a SONNET agent (agent K), Python only (numpy/scipy; mpmath/sympy for checks; no Wolfram, no wolfbook tools), run from
`Paul/Python/gaudin_sl3/` with `python -X utf8`. Plan: `Paul/Mathematica/Gaudin/Gaudin-sl3-round6-plan.md` §0-§1 (ingredient (d)), §3 (row K), §5
(K fallbacks), §6 (conventions). The construction being scaled up is round 5's M2/M3: `RESULTS-master-kernel-2.md` §1 ("K5" and "K5 extra:
universality"), `master_kernel2.py`, `RESULTS-jet-basis-sym.md` + `jet_basis_sym.py`, `RESULTS-general-kernel.md` (T-format, block theorem, exchange
rule, top-state sum rule) + `general_kernel.py`. Other agents (R, S, D, Q0) run concurrently: <= 2 heavy processes of your own, total memory
< ~12 GB (the machine has 32 GB shared), nobody commits to git, existing scripts are READ-ONLY (import them; never edit). **Files you own:**
`kernel_L14.py`, `RESULTS-kernel-L14.md`, `kernel_N_L8.npz`, `kernel_N_L14.npz`, `kernel_Theta_840.npz`, `kernel_Theta_1470.npz` (+ `kernel_M_L*.npz` and other
`kernel_*.npz` caches).  NOTE: the repository's `.gitignore` ignores `Paul/Python/gaudin_sl3/*.pkl` only -- `*.npz` is NOT ignored (the files below are up to 0.4 GB).

(Revision note, end of the task: §2 below was written first as the plan and has been rewritten to the route that was actually implemented and verified;
every statement in §1-§2 is now backed by a check listed in `RESULTS-kernel-L14.md`.)

## 0. Goal and deliverables

Needed by the Q-route at [14,7,0]^2 (round 6 plan §1 row (d)): the kernel Theta^{[14,7,0]} on the 92 non-forced jet pairs per site (l >= 8), i.e. 92^2 = 8464
entries, each a polynomial in (a, a', t) (a-degree <= 21, t-degree <= 4*13 - 11 = 41), about 1.7*10^6 coefficients; likewise Theta^{[8,4,0]} (35 pairs, 1225
entries). Everything follows from the ONE universal integer matrix N (below) on jet-pair labels <= L+1, which for L = 14 is 14400 x 14400 (exact
rational/sympy inversion is impossible; M3 did L <= 6 only).

Deliverables: `SPEC-kernel-L14.md` (this file, written first), `kernel_L14.py` (library + CLI: build M^(L), N^(L), Theta^{[L,l2,0]} via the block theorem, checks),
`RESULTS-kernel-L14.md` (style of the round-5 RESULTS), `kernel_N_L8.npz`, `kernel_N_L14.npz` (N on labels <= L+1 and the exact M; formats in RESULTS),
`kernel_Theta_840.npz`, `kernel_Theta_1470.npz` (Theta on the non-forced pairs: pair list, per entry the monomial exponents (n,p,q) and coefficients, for
Theta_ij = sum c t^n a^p a'^q, plus the EXACT coefficients as Garner digits).

**Acceptance (all four, with numbers in RESULTS):**
(i) the pipeline reproduces the exact kernels `general_kernel_rep_{300,400}.pkl` (via `general_kernel.load_rep`) and `master_kernel2_rep_{500,600}.pkl` (via
`master_kernel2.get_T`) to >= 1e-10 relative on EVERY monomial coefficient; also the block-theorem kernels [1,1,0], [2,1,0], [2,2,0], [3,1,0], [3,2,0], [3,3,0], [4,4,0] (stored `general_kernel_rep_*.pkl`);
(ii) universality: N^(L)|_{labels <= 7} equals `master_kernel2_N6.pkl` (exact integers) for every L built (8, 10, 12, 14), and M^(L)|_{<= L'+1} = M^(L') exactly;
(iii) the new kernels satisfy the exchange rule Theta_ij(a,a',t) = Theta_ji(a',a,-t) and the top-state sum rule
sum_ij Theta_ij m_i^top(th1) m_j^top(th2) = M^top(th1) M^top(th2) R(a) R(a'), R(a) = (1 - z2 a)^{l1-l2}(1 - z3 a)^{l1}
(`general_kernel.top_monos_exact`, `top_sum_rule_check` are the exact reference implementations);
(iv) memory/time for L = 8 and L = 14, the conditioning (largest coefficient of N and of Theta, the accuracy of the double-precision objects against exact references,
the residual of M N - 1), and a precise statement of the (a, a', t) ranges in which the double-precision kernel can be trusted (the eps-scaling a = a^/eps of the untwisting
limit needs the polynomial COEFFICIENTS, not values: say how accurate each coefficient is).

## 1. Facts you can use (all exact, from rounds 4-5; do not re-derive)

* One-row Lambda = [L,0,0]: V = C[x1,x2,x3]_L (E_ij = x_i d_j, column convention E x^alpha = alpha_j x^{alpha+e_i-e_j}); two sites: E^(1) = E x 1, E^(2) = 1 x E.
  Jet pairs (k,l), 0 <= k < l <= L+1 (count C(L+2,2) = dim V = 120 at L = 14), jet level rho(k,l) = k+l-1, two sites: rho_r = rho(r1) + rho(r2).
* Y(t) (rows = jet-pair covectors <i,j|, columns = kets) is graded: Y(t)_{rc} = t^{rho_r - kappa_c} Y(1)_{rc}, kappa_c = alpha_2 + 2 alpha_3 of the two kets summed.
  Y(0) = Y_0 = the level-preserving part = y (x) y (one-site jet matrices). M(1) := Y(1) Y_0^{-1} is UNIT block-lower-triangular in the level rho, has
  INTEGER entries (L = 6: 46025 nonzeros of 784^2, max 241200), is label-triangular (M_{r,p} != 0 only if l(r1) >= l(p1) and l(r2) >= l(p2)),
  and is UNIVERSAL: the restriction of M^(L)(1) to labels <= L'+1 equals M^(L')(1), same for N := M^{-1} (an integer matrix too; max entry 2.6e17 at L = 6).
  N^(6) exact: `master_kernel2_N6.pkl` (dict {(row, col): sympy Integer}, index r = i*28+j over `general_kernel.pairs_of((6,0,0))`).
* Kernel from N (exact, checked L = 2..6): Theta^{[L,0,0]}_{pq}(a,a',t) = sum_{r,s} c_r c_s a^{rho_r} a'^{rho_s} t^{rho_r+rho_s-rho_p-rho_q} N_{rs;pq},
  c_kl(L) = (-1)^{k+l-1}(l-k) C(L,l-1) C(L+1,k)/l (= `general_kernel.ckl_closed((L,0,0),k,l)`), r = (r1,r2) etc., (p,q) over ALL pairs of [L,0,0].
* Block theorem: Theta^{[l1,l2,0]}_{ij} = [l2!(l1-l2)!/l1!]^2 (a a')^{-l2} Theta^{[l1,0,0]}_{ij} for the non-forced pairs (0 <= k < l <= l1+1, l >= l2+1) i,j.
  By label-triangularity of N this needs only the block N[S,S] with S = {pairs with l >= l2+1} (= (M[S,S])^{-1}).
* T-format: T[(i,j)] = {(n,p,q): coeff}, Theta_ij = sum coeff t^n a^p a'^q, (i,j) indices into `pairs_of(lam)`^2; homogeneity p+q-n = e_i+e_j (e = k+l-l2-1),
  0 <= p,q <= 2 l1 - l2, t-degree 4(l1-1) - max(0,2 l2 - 3); exchange rule exact; Theta(t=0) = theta(a) x theta(a') with theta_kl = c_kl a^{k+l-l2-1}.
* The operator route of M2 (`jet_basis_sym.Y_from_operators`): row (i,j) of Y = <x0| nu_i(site-1 data) nu_j(site-2 data), x0 = coefficient of x1^L y1^L,
  nu_kl = k! l! (gamma_k delta_l - gamma_l delta_k), gamma/delta = Frobenius coefficients of the two local solutions f0 = 1 + 0 x + ..., f1 = x + ... of the dual
  equation psi''' - S1 psi'' + S2 psi' - S3 psi = 0 at the pole theta_s (exponents {0,1,L+2}; `frobenius_planes(L)`), with Laurent data S1 = e1 + L/(u-th_s) + L/(u-th_t),
  S2 = H_s/x + (regular, H_t/(x+delta) + e2), S3 = b_s/x + b_t/(x+delta) + e3 (a_s = a_t = 0 and the S2 double pole, S3 triple pole vanish for one-row reps),
  x = u - th_s, delta = th_s - th_t; H_s, b_s are the Laurent coefficients of the quantum characteristic determinant (`gaudin_chain.symbolic_cdet`,
  `Operators.fit`, `jet_basis_sym.companion_ops`); Y is independent of the twist zeta (polynomial in e_k) and polynomial of degree L in t = 1/th12.

## 2. The route that was implemented (all steps verified, see RESULTS)

Work in the MONOMIAL basis (integer generators, no GT normalisation scalars), at the point **zeta = 0 (nilpotent twist g = E21 + E32, e_k = 0), th1 - th2 = 1 (t = 1,
delta = +1 at site 1, -1 at site 2)**; by the grading M(t)_{rp} = t^{rho_r-rho_p} M(1)_{rp} no information is lost (checked at th12 = 2, -3, 1/2 in float64, L = 2..4).

1. **Operators** (`build_ops`). Single-site E_ij (d x d, scipy.sparse, one nonzero per column), two-site E^(s) by Kronecker products. Laurent expansion (class `LS`, matrix coefficients,
   ordered products) of l_ij(u) = g_ij + E^(s)_ij/x + E^(t)_ij/(x+delta) inserted into the 8 monomials of S2 and 13 of S3 (`gaudin_chain.symbolic_cdet()`):
   H_s = [x^-1] S2, b_s = [x^-1] S3; a_s = [x^-2] S3, [x^-2] S2, [x^-3] S3 are asserted to vanish. Verified against `jet_basis_sym.companion_ops` (exact, GT basis) at L = 2, 3:
   0 mismatching entries (H1, H2, b1, b2), identical nonzero patterns. One-site analogues h, b for y.
2. **nu as a vector recursion** (`frob_batch`; never an expanded polynomial). With covectors as columns and the right action by transposed operators, for f0 (rho = 0, c0 = v, c1 = 0) and
   f1 (rho = 1, c0 = v), n_k = rho + k, I(n) = n(n-1)(n-2-L), inv1(p) = (-1)^p delta^(-p-1):
   c_m I(n_m) = - [ n_{m-1} c_{m-1}H_s - c_{m-2} b_s - sum_{k<m} L inv1(m-k-1) n_k(n_k-1) c_k + sum_{k<=m-2} inv1(m-k-2) n_k c_k H_t - sum_{k<=m-3} inv1(m-k-3) c_k b_t ];
   gamma_k = c_k (f0, gamma_1 = 0), delta_l = c_{l-1} (f1, delta_0 = 0). Checked against the symbolic `nu_formulas`/`frobenius_planes` with rational scalar data (L = 2, 3, 4; delta = 1, -1, 2).
   All rows: (site 1) one f1-run from x0 and one f0-run started from the resulting L+2 vectors (batch), (site 2) the same on the batch of the C(L+2,2) site-1 rows: 17 batched runs.
3. **y and M** (`y_exact`, `y_inverse_exact`, `build_Y_rows`, `M_from_Y`). y_{kl} = x0 nu^{(0)}_{kl}(h, b) from the one-site t = 0 recursion in exact Fractions (block diagonal in the principal weight,
   blocks <= 8 x 8); M(1)_r = (y^-1)^T Y_r y^-1 (Y_r reshaped d x d) = Y K with the sparse K = kron(y^-1, y^-1).
4. **Exactness is essential.** The plain float64 route (same code, `mod=None`) gives M with absolute error 9e-16, 2e-13, 5e-11, 1.5e-8, 4.7e-6, 1.5e-2, 4.8, 1.3e4, 1.7e7 for L = 2..10 (digits lost ~ 2.5 per L):
   integers are no longer recoverable at L >= 8. The recursion is therefore done in exact modular arithmetic (int64, primes just below 2^28, so that sums of <= 64 products stay < 2^63;
   kron(Z,Z) must be reduced mod p BEFORE the product), 3-prime CRT lift (|M| < 2^62; max|M^(14)| = 2.5e18) + a 4th prime that verifies every entry (`build_M_exact`, batches of site-1 rows).
5. **N** (`ModSolver`, `N_exact_columns`). N = M^{-1} exactly by block forward substitution over the levels (window of the last L levels) with float64 BLAS modulo primes p < 790000
   (all GEMM sums are integers < 2^53), Garner mixed-radix lift with symmetric digits on the union support (`garner_signed`), the number of primes actually needed is reported as `depth` (must be <= #primes - 2).
   Practical: OpenBLAS with 20 threads is ~100x slower for small-m dgemm (80 x 1000 @ 1000 x 1000: 1 GFlop/s vs 250-300 for the transposed product) -- the solver keeps the transposed layout
   (columns of N as rows) and caps OPENBLAS_NUM_THREADS at 8. The plain float64 forward substitution from the exact M is kept only for the error statistics.
6. **Theta** (`theta_exact`, `theta_to_arrays`). Only N[S x S, S x S], S = {l >= l2+1}; weights c~_r = S0 c_r (S0 = lcm(1..15) = 360360, integers); C~[x,y,(p,q)] = sum_{rho_r=x, rho_s=y} c~_r c~_s N_{rs;pq}
   exactly mod each prime (two small GEMMs), Garner; coefficient of t^{x+y-rho_p-rho_q} a^{x-l2} a'^{y-l2} is [l2!(l1-l2)!/l1!]^2 C~/S0^2 (exact Fraction); stored as float64 AND as the exact Garner digits.

## 3. Fallbacks (plan §5) -- not needed

L = 14 turned out feasible (M: 6.5 min, 3.6 GB; N: 5 min, 8.8 GB; Theta^{[14,7,0]}: 2 min, 3.2 GB). The stopgap named in the plan (state-independent K-tensor solve on the chain at the m+1
needed points) is NOT part of this spec and was not used.

## 4. Reporting

`RESULTS-kernel-L14.md`: files, conventions, headline results (numbers), tables (time/memory/conditioning vs L; acceptance residuals), surprises, not done. End with a concise
summary of the key numbers.
