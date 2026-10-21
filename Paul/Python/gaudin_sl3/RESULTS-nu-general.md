# RESULTS N1: the nu-formulas for general Lambda (B1) and the untwisting limit in nu-form for [3,1,0] and [3,2,0] (B2) — agent N1, round 5, 2026-10-16

Files: `nu_general.py` (all code; `python -X utf8 nu_general.py B1` = B1 verification on cached states, `python -X utf8 nu_general.py RUN 310 24` = full B2 analysis of a rep with 24 eps-points, both phi),
caches `nu_general_reg_<lam>[<tag>]_<k>.pkl` (regime points: 24 eps x 2 phi x {210, 310, 320}, 200 digits, ~2.3 s each), `nu_general_res_<lam>[2].pkl` (analysis records, phi / phi'). Python only; imports `general_kernel`, `gaudin_chain`, `twovac_fit.cov_inf` read-only; W4's `ref/untwisted-exact-{210,310,320}.txt` parsed by `load_w4`.
`ref/laurent-series-{310,320}.txt` (W5) did not exist when I finished; the only exact cross-check of my numerical eps-series is therefore the eps^0 term against W4 (below).

## 0. Conventions and what was assumed

* Dual oper (d^3 - S1 d^2 + S2 d - S3) P = 0, x = u - theta_s, delta = theta_s - theta_t, sigma_1 = zeta_1+zeta_2+zeta_3, e_2, e_3 elementary symmetric functions of zeta,
  S1 = |L|/x + sum s_k x^k, S2 = c2/x^2 + H/x + sum h_k x^k, S3 = a/x^2 + b/x + sum t_k x^k, c2 = lam2 (lam1+1), |L| = lam1+lam2.
  Regular parts (exact, from the other site's poles and the constants): s_k = |L| (-1)^k delta^{-k-1} + sigma_1 d_k0, h_k = H_t (-1)^k delta^{-k-1} + c2 (-1)^k (k+1) delta^{-k-2} + e_2 d_k0, t_k = b_t (-1)^k delta^{-k-1} + a_t (-1)^k (k+1) delta^{-k-2} + e_3 d_k0.
* **Partial-fraction form verified exactly on all 225 states** of [3,1,0], [3,2,0] (and 64 of [2,1,0]): S2 constant = e_2(zeta), S3 constant = e_3(zeta), both double-pole coefficients of S2 equal c2 = 4, 8, and S3 has NO triple pole (max deviation 4e-81 at 80 digits). Indicial polynomial rho(m) = m(m-1)(m-2) - |L| m(m-1) + c2 m = m (m - lam2 - 1)(m - lam1 - 2).
* Laurent data of a state in the code: `lam2 = [e2, H_1, c2, H_2, c2]`, `lam3 = [e3, b_1, a_1, 0, b_2, a_2, 0]` (`gaudin_chain` basis); `laurent_of(lam2, lam3, s)` -> (H_s, a_s, b_s).
* Pair order (`pairs_of`): (k,l), 0<=k<l<=lam1+1, not both <= lam2. [2,1,0]: (02)(03)(12)(13)(23); [3,1,0]: (02)(03)(04)(12)(13)(14)(23)(24)(34); [3,2,0]: (03)(04)(13)(14)(23)(24)(34).
* B2 regime: theta = (2/7,-3/5); zeta = eps phi, phi = (1/3,-2/7,5/11), phi' = (-2/5,3/8,7/11); omega = eps w^1, omega' = eps w^1'; eps_k = 1e-6 2^-k; 200 digits; **k = 0..23 (24 points, not 12)**: with 12 points the leading coefficients of the high orders (n up to 9) are only good to 1e-26 ... 1e-66 (numbers below), 24 points give 1e-96 ... 1e-163; points are cheap (2.3 s).
  Branches: per HW weight (HW sectors from W4's files, 10 branches per rep) the weight space is diagonalised separately (`mp.eig`, largest block 18x18), branch = state(s) of smallest raising residual (selected residual <= 0.8 eps, the next state's >= 1.41 at all 24 eps, both phi, both reps). **Copy convention (multiplicity 2): copy 1 = smaller a_1** (a_1 = coefficient of 1/(u-theta_1)^2 in S3). Verified against W4's `BRANCH` records: the eps^0 coefficients of (H_s, a_s, b_s) of all 10 branches x 2 sites equal W4's eta/theta12, alpha/theta12, beta/theta12^2 to 9e-198 / 2e-197 ([3,1,0] / [3,2,0]); the vacua Phi_inf are exactly eps-free (exact Fraction check) in both reps.
  Psi_B(eps): G-normalised (G = Gram diag of the GT matrices), sign made continuous in eps (max |Psi(eps_k)-Psi(eps_{k+1})|/eps_k <= 0.06 for all 20 branches, no flips) and global sign = first nonzero GT component (W4's order) > 0 at eps^0.
* x(zeta): W4's `XVEC` strings parsed and compared with `general_kernel.x_cov`: **identical as exact rationals** for 3 zeta points, all three reps; homogeneous of degree lam2 in zeta (checked symbolically). x0^G = x (x) x.
* Kernel: `G = (w1 w1')^{2 lam1 - lam2} Theta`, Theta = `load_rep(lam)` (exact), t = 1/theta12 = 35/31; in the regime G_ij = sum c t^n eps^{4 lam1-2 lam2-p-q} w^^{2lam1-lam2-p} w'^^{2lam1-lam2-q}. c-hat = c/eps^{4 lam1} = (phi1-phi2)^{2 lam2} prod_{w in {w^1,w^1'}} (phi2-w)^{lam1-lam2}(phi3-w)^{lam1}.

## 1. B1 — nu for general Lambda (exact, sympy)

**Setting.** P = sum c_m x^m. For m not in {lam2+1, lam1+2} the recursion fixes c_m:

  rho(m) c_m = sum_k s_k (m-1-k)(m-2-k) c_{m-1-k} - (H(m-1) - a) c_{m-1} + b c_{m-2} - sum_k h_k (m-2-k) c_{m-2-k} + sum_k t_k c_{m-3-k}.

c_0 and c_{lam2+1} are free; every solution is c_m = c_0 u_m + c_{lam2+1} v_m + c_{lam1+2} w_m, and for the jets needed (m <= lam1+1) only u, v matter (u_0 = 1, u_{lam2+1} = 0; v_m = 0 for m <= lam2, v_{lam2+1} = 1). Hence
  **nu_kl = k! l! (u_k v_l - u_l v_k)/(lam2+1)!  (non-forced pairs), nu_{0,lam2+1} = 1**, and for k <= lam2 < l the minor factorises: **nu_kl = k! u_k nu_{0l}**; only pairs lam2 < k < l are genuine 2x2 minors.
  All nu are POLYNOMIALS in (H_s, a_s, b_s, H_t, a_t, b_t, sigma_1, e_2, e_3, 1/delta) with rational coefficients (no division by Laurent data). Code: `nu_symbolic(lam1, lam2)` (returns pairs, nu list, three no-log polynomials), `nu_numeric(lam, laurent_s, laurent_t, zeta, theta_s, theta_t)` (mpmath), `frobenius(lam)` (u, v, no-log).
  The factor D_s = c_0^(2) c_{lam2+1}^(3) - c_0^(3) c_{lam2+1}^(2) common to every Plucker coordinate p_kl = D_s k! l!(u_k v_l - u_l v_k) cancels against M = P^3 D_s (lam2+1)!, so c Qnum = sum Theta_ij (w w')^{2lam1-lam2} nu_i(theta_1) nu_j(theta_2) with Laurent data only (this is exact algebra; verified numerically below, and the full limit statement in section 2).

**Adjoint [2,1,0]** (c2 = 3): nu = (nu02, nu03, nu12, nu13, nu23) = (1, pi, kappa, kappa pi, eps_s), kappa = a/3, pi = 2H - a - 2 sigma_1 - 6/delta; eps_s identical to U2's (symbolic difference 0).

**[3,1,0]** (c2 = 4, |L| = 4), b_s eliminated with b_s = a_s(H_s - a_s)/4:
* nu02 = 1, nu12 = kappa = a/4, nu03 = pi = H - a/2 - sigma_1 - 4/delta, nu13 = kappa pi,
* nu04 = [12 H^2 d^2 - 13 H a d^2 - 36 H d^2 s1 - 144 H d + 24 H_t d + 5 a^2 d^2 + 16 a d^2 s1 + 64 a d + 24 d^2 e2 + 24 d^2 s1^2 + 192 d s1 + 576]/(8 d^2)  (d = delta, s1 = sigma_1), nu14 = kappa nu04,
* nu23 = -[ -H a^2 d^2 + 4 H_t a d + a^3 d^2 + 4 a d^2 e2 + 16 a - 16 a_t - 16 b_t d - 16 d^2 e3 ]/(16 d^2),
* nu24, nu34: longer (nu24 involves a_t, b_t, H_t, e_2, e_3, s1; nu34 is of degree 5 in the data) — `nu_symbolic(3,1)`.
Un-eliminated: nu03 = (2 H d - a d - 2 d s1 - 8)/(2 d), nu04 contains -6 b_s d^2 in the numerator (/(4 d^2)): `nu_symbolic(3,1)[1][2]`.

**[3,2,0]** (c2 = 8, |L| = 5; b_s NOT eliminable by a site-s-only relation):
* nu03 = 1, nu13 = kappa = a/8, nu04 = pi = 3 H - a - 6 sigma_1 - 30/delta, nu14 = kappa pi,
* nu23 = 2 u_2 = -(H a - a^2 - 8 b)/24,  nu24 = nu23 nu04,
* nu34 = [2 H H_t a d^2 - H a b d^3 + 2 H a d^3 e2 + 26 H a d - 2 H_t a^2 d^2 + 6 H_t a d - 16 H_t b d^2 + a^2 b d^3 - 2 a^2 d^3 e2 - 26 a^2 d + 6 a a_t d + 6 a b_t d^2 + 6 a d^3 e3 + 96 a - 96 a_t + 8 b^2 d^3 - 16 b d^3 e2 - 208 b d - 48 b_t d]/(8 d^3).
(Note that for lam2 = 2 the order-2 recursion is not resonant, so nu23 != 0 and b enters through u_2.)

**No-log conditions** (3 per site: one at order lam2+1 [call A], two at order lam1+2 [B1: coefficient of c_0, B2: coefficient of c_{lam2+1}]; exactly the "3 independent conditions" of P3/P4; they hold on all 225 states, residual <= 6e-77 relative, 80 digits):
* **A, lam2 = 1, any lam1:  (lam1+1) b_s = a_s (H_s - a_s)**  (adjoint /3, [3,1,0] /4, [4,1,0] /5, [5,1,0] /6, checked symbolically). Only (H_s, a_s, b_s) enter.
* **A, lam2 = 2, any lam1** (checked against lam1 = 3,4,5 by the exact recursion; d = delta, s1 = sigma_1, l1 = lam1):
  4 d^2 l1 (l1+1) A = -2 H^2 a d^2 + 3 H a^2 d^2 + 2 H a d^2 s1 + 2 H a d l1 + 4 H a d + 4 H b d^2 (l1+1) + 2 H_t a d l1 - a^3 d^2 - 2 a^2 d^2 s1 - 2 a^2 d l1 - 4 a^2 d - 2 a b d^2 (2 l1 + 1) + 2 a d^2 e2 l1 + 4 a l1 (l1+1) - 4 a_t l1 (l1+1) - 4 b d^2 s1 (l1+1) - 4 b d (l1^2 + 3 l1 + 2) - 4 b_t d l1 (l1+1) - 4 d^2 e3 l1 (l1+1)  (= 0).
  It is cubic in (H, a, b), linear in b_s, and involves the OTHER site (H_t, a_t, b_t) and e_2, e_3, sigma_1. For [3,2,0] (l1 = 3) it is P4's order-x^3 condition (identical coefficients up to the factor 192).
  **Which variables enter A: lam2 <= 1: only (H_s,a_s,b_s); lam2 >= 2: also (H_t, a_t, b_t, delta, e_2, e_3, sigma_1); degree in b_s: 1 for lam2 <= 2, 2 for lam2 = 3** (lam2 = 3 expressions exist via `frobenius((l1,3))`, not written out).
* B1, B2 of [3,2,0]: B2 = [12 H^2 d^2 - 7 H a d^2 - 60 H d^2 s1 - 300 H d + 12 H_t d + a^2 d^2 + 18 a d^2 s1 + 90 a d - 4 b d^2 + 12 d^2 e2 + 72 d^2 s1^2 + 720 d s1 + 2016]/(4 d^2) = 0: **quadratic, inhomogeneous, only (H_s, a_s, b_s, H_t)** (this is the quadratic P4 found); B1 is cubic with both sites (1602 characters, `frobenius((3,2))['nolog'][1]`).
  [3,1,0]: B2 = cubic in (H,a,b) with H_t (604 characters), B1 = degree 4 with both sites (1901 characters). Adjoint: B2 = 3 d^2 B2 = 6 H^2 d^2 - 5 H a d^2 - 18 H d^2 s1 - 54 H d + 6 H_t d + a^2 d^2 + 8 a d^2 s1 + 24 a d - 3 b d^2 + 6 d^2 e2 + 12 d^2 s1^2 + 72 d s1 + 144.

**Verification on actual states** (`verify_B1`; m_i = P3 p_kl from `site_monos` at G1's cached base points, M = P3 p_{0,lam2+1}; ALL 225 states x 2 sites x all pairs; 80 digits):

| rep | data set (zeta tag, theta) | max |m_i/M - nu_i| | no-log A, B1, B2 residual |
|---|---|---|---|
| [3,1,0] | P, (2/7,-3/5) | 6.9e-74 | 7e-79, 4e-77, 6e-77 |
| [3,1,0] | Pp, (2/7,-3/5) | 1.5e-73 | 3e-78, 2e-77, 2e-77 |
| [3,1,0] | P, (2/7,-113/315) | 1.3e-74 | 8e-79, 2e-77, 7e-78 |
| [3,2,0] | P, (2/7,-3/5) | 1.4e-73 | 7e-79, 2e-78, 2e-78 |
| [3,2,0] | Pp, (2/7,-3/5) | 9.5e-75 | 6e-79, 3e-78, 8e-79 |
| [3,2,0] | P, (2/7,-113/315) | 2.9e-74 | 1e-79, 1e-77, 3e-78 |
| [2,1,0] | P, (2/7,-3/5) (64 states) | 1.6e-76 | 2e-79, 4e-78, 5e-78 |

(limited by the 80-digit data, i.e. agreement to the working precision.) The adjoint nu symbolically equals U2's `nu_funcs` (difference 0 for kappa, pi, eps).

## 2. B2 — the untwisting limit in nu-form

**Headline: the nu-reduction holds for both reps, with the SAME constant 1 for every branch.** For all 20 branches (10 per rep), both phi:
  **sum_ij G_ij(eps w^1, eps w^1'; theta12) nu_i(theta_1) nu_j(theta_2) = eps^{N_B} c-hat <phi_u (x) e1|Psi_0> / X_B + O(eps^{N_B+1})**,
  with **constant C = 1 (|C - 1| <= 2.3e-96 over the 10 branches of [3,1,0], <= 4e-117 over those of [3,2,0], both phi, 24 points)**, i.e. NO rep- or branch-dependent constant appears in this normalisation (G = (ww')^{2lam1-lam2} Theta with Theta from G1's `load_rep`; for the adjoint the same code returns C = 1 too, k_B and N_B as U2).
  * **k_B = height of 2Lambda - mu in the root lattice = 2(lam1 - lam3) - (mu1 - mu3)**; **N_B = 4 lam1 - 2 lam2 - k_B** (all 20 + the 6 adjoint branches; N_B is the first eps-order where sum G nu nu is nonzero: lower orders vanish to 1e-103 ... 1e-170).
  * X_B = <x(phi)(x)x(phi)|Psi_{k_B}>, k_B = first n with it nonzero (|X_n| > 1e-25); on all exact HW kets it vanishes, as in the adjoint.

### [3,1,0] (c2 = 4, exponents {0,2,5}); 10 branches; phi = (1/3,-2/7,5/11), phi' in brackets

| branch (mu) | k_B | N_B (pred 10-k_B) | X_B(phi) | X_B(phi') | |C-1| | identity residual |L_B - C c-hat <phi_u e1|Psi_0>/X_B| / |L_B| (24 pts) | same with 12 pts | max n of nu^(n) needed | nu_i lead exponents (i = 02,03,04,12,13,14,23,24,34) |
|---|---|---|---|---|---|---|---|---|---|
| (6,2,0) | 0 | 10 | 0.38321995 | 0.600625 | 6e-98 | 7e-97 | 2.5e-129 | 5 | 0,0,1,1,1,2,2,3,4 |
| (6,1,1) | 1 | 9 | -0.017394137 | 0.065736668 | 2e-96 | 9e-96 | 2.7e-26 | 9 | 0,1,2,0,1,2,2,3,4 |
| (5,3,0) | 1 | 9 | -0.05252981 | 0.10307154 | 6e-108 | 3e-107 | 6.4e-31 | 9 | 0,0,0,1,1,1,2,2,3 |
| (5,2,1) copy 1 | 2 | 8 | 0.0010985815 | 0.0051977418 | 1e-106 | 5e-106 | 4.1e-38 | 8 | 0,0,1,0,0,1,1,2,3 |
| (5,2,1) copy 2 | 2 | 8 | (same) | (same) | 9e-107 | 3e-106 | 6.0e-38 | 8 | (same) |
| (4,4,0) | 2 | 8 | 0.011085918 | 0.027232169 | 9e-118 | 3e-117 | 1.3e-32 | 8 | 0,0,0,1,1,1,2,2,3 |
| (4,3,1) copy 1 | 3 | 7 | -0.00013694922 | 0.00081118448 | 4e-118 | 2e-117 | 8e-46 | 7 | 0,0,0,0,0,0,1,1,2 |
| (4,3,1) copy 2 | 3 | 7 | (same) | (same) | 2e-117 | 1e-116 | 1.4e-45 | 7 | (same) |
| (4,2,2) | 4 | 6 | 3.3547722e-6 | 4.7915124e-5 | 3e-130 | 8e-130 | 9e-53 | 6 | 0,0,1,0,0,1,0,1,2 |
| (3,3,2) | 5 | 5 | -6.3360545e-7 | 1.132938e-5 | 4e-132 | 4e-151 | 1.5e-66 | 5 | 0,0,0,0,0,0,0,0,1 |

The two copies of each multiplicity-2 sector have identical k_B, N_B, X_B and exponents (as in the adjoint). The L_B are polynomials of degree <= 2 lam1 - lam2 = 5 in each of w^1, w^1' (4 for the last two branches), total degree between 0..3 (lowest) and N_B (highest), 12 ... 36 monomials. Branch exponents of the Laurent data are the same at both sites and nu_i^(n)(theta_s) at the SAME order for both sites.

### [3,2,0] (c2 = 8, exponents {0,3,5}); 10 branches

| branch (mu) | k_B | N_B (pred 8-k_B) | X_B(phi) | X_B(phi') | |C-1| | identity residual (24 pts) | 12 pts | max n of nu^(n) | nu_i lead exponents (i = 03,04,13,14,23,24,34) |
|---|---|---|---|---|---|---|---|---|---|
| (6,4,0) | 0 | 8 | 0.14685753 | 0.36075039 | 1e-117 | 2e-117 | 9e-143 | 4 | 0,0,1,1,2,2,3 |
| (6,3,1) | 1 | 7 | -0.0047134186 | 0.027918758 | 2e-118 | 7e-118 | 2.6e-44 | 7 | 0,0,0,0,1,1,2 |
| (6,2,2) | 2 | 6 | 0.00023290773 | 0.0033265456 | 9e-128 | 2e-127 | 1.7e-47 | 6 | 0,1,0,1,0,1,2 |
| (5,5,0) | 1 | 7 | -0.028468786 | 0.087550205 | 4e-128 | 1e-127 | 9.9e-41 | 7 | 0,0,1,1,2,2,3 |
| (5,4,1) copy 1 | 2 | 6 | 0.00042099833 | 0.0031218937 | 1e-128 | 4e-128 | 8e-53 | 6 | 0,0,0,0,1,1,2 |
| (5,4,1) copy 2 | 2 | 6 | (same) | (same) | 6e-128 | 2e-127 | 1.3e-52 | 6 | (same) |
| (5,3,2) copy 1 | 3 | 5 | -1.2288241e-5 | 0.00021972373 | 8e-141 | 3e-140 | 8e-59 | 5 | 0,0,0,0,0,0,1 |
| (5,3,2) copy 2 | 3 | 5 | (same) | (same) | 2e-140 | 1e-139 | 1.3e-58 | 5 | (same) |
| (4,4,2) | 4 | 4 | 1.2856156e-6 | 2.8779022e-5 | 3e-144 | 2e-151 | 1.9e-69 | 4 | 0,0,0,0,0,0,1 |
| (4,3,3) | 5 | 3 | -5.6852437e-8 | 3.0687638e-6 | 4e-132 | 3e-163 | 7.4e-74 | 3 | 0,0,0,0,0,0,0 |

L_B has degree <= 4 in each variable (3 for the last two).

### (ii)-(vi) items

* **(i) leading polynomial**: L_B(w^1,w^1') = sum over the contributions with d + n + m = N_B exactly (exact kernel polynomial x fitted nu-series), lower orders cancel: max |R_t| for t < N_B is 1e-103 ... 1e-170 ([3,1,0]: 3e-103, 6e-102, 2e-113, 1e-111/4e-112, 1e-123, 5e-123/2e-122, 5e-136, 2e-157; [3,2,0]: 4e-125, 9e-126, 6e-135, 2e-135, 3e-135/8e-135, 2e-147/5e-147, 2e-158, 1e-170), against O(1) at t = N_B.
* **(iii) the identity**: checked at the 64 grid points (8 x 8: w^1 in {3/4,-5/3,1/2,-7/3,9/7,-4/9,11/5,-13/6} x w^1' in {-12/5,7/2,1/8,5/4,-8/3,13/9,-3/7,17/6}, which contains U2's 16-point grid) as a polynomial/rational-function identity (bidegree <= 5 per variable, 64 points suffice); **the ratio L_B / (c-hat <phi_u e1|Psi_0>/X_B) is constant over the 64 points to <= 8e-87 ([3,1,0]) / <= 3e-108 ([3,2,0]) relative (`spread`) and equals 1 to <= 2.3e-96 / <= 4e-117**.
* **(iv) CG dictionary**: solving L_B = c-hat sum_nu C_nu [e1^T rho(n')]_nu in the exact HW space (W4's HWBASIS), G-normalising and fixing W4's sign returns W4's exact CG slice: max |C_nu^Q - C_nu^exact| over all nu, all 20 branches, **0 (for the 12 branches with a rational exact slice: <= 3e-203) and <= 4e-106 ([3,1,0]) / <= 3e-128 ([3,2,0]) for the 8 branches in multiplicity-2 sectors (surd entries)**; the whole 225-component vector agrees with the chain's Psi_0 to the same level (<= 6e-106 / 5e-128; and the chain's Psi_0 equals W4's exact vector to 3e-200, so the copy convention matches).
  The 20 tables (C_nu = <e_nu (x) e_1|Psi_0>, W4's GT index nu = 1..15, only nonzero entries, 12 digits):
  [3,1,0]: (6,2,0): nu1=1; (6,1,1): nu2=-0.707106781187; (5,3,0): nu3=-0.5; (5,2,1)#1: nu4=0.0564580373849, nu5=0.0983443475906; (5,2,1)#2: nu4=0.299890285165, nu5=-0.0538008072719; (4,4,0): nu7=0.288675134595; (4,3,1)#1: nu8=-0.00744303312309, nu9=-0.0758903002102; (4,3,1)#2: nu8=-0.103668077988, nu9=0.0203347446547; (4,2,2): nu10=0.037267799625, nu11=-0.012422599875; (3,3,2): nu14=-0.0228217732294.
  [3,2,0]: (6,4,0): nu1=1; (6,3,1): nu2=-0.5; (6,2,2): nu3=0.288675134595; (5,5,0): nu4=-0.707106781187; (5,4,1)#1: nu5=0.0423435280387, nu7=0.203745949854; (5,4,1)#2: nu5=0.224917713874, nu7=-0.0701153288982; (5,3,2)#1: nu6=-0.00744303312309, nu8=-0.0758903002102; (5,3,2)#2: nu6=-0.103668077988, nu8=0.0203347446547; (4,4,2): nu12=0.07453559925; (4,3,3): nu13=-0.0228217732294.
  (GT index/weight lists: [3,1,0]: nu1 (3,1,0), nu2 (3,0,1), nu3 (2,2,0), nu4,5 (2,1,1), nu7 (1,3,0), nu8,9 (1,2,1), nu10,11 (1,1,2), nu14 (0,2,2); [3,2,0]: nu1 (3,2,0), nu2 (3,1,1), nu3 (3,0,2), nu4 (2,3,0), nu5,7 (2,2,1), nu6,8 (2,1,2), nu12 (1,2,2), nu13 (1,1,3).) Observation: under the duality [3,2,0] = [3,1,0]* x det^3 (mu -> (6-mu3, 6-mu2, 6-mu1)) the numbers of 5 pairs coincide: (6,1,1)/(5,5,0), (5,3,0)/(6,3,1), (4,4,0)/(6,2,2), (4,3,1)#a,b/(5,3,2)#a,b, (3,3,2)/(4,3,3) (and (4,2,2) -> (4,4,2): 2 x 0.0373 = 0.0745 for the dominant component), the (5,2,1)/(5,4,1) ones do not (different second-site covector).
* **(v) phi-independence**: phi vs phi': C_nu differ by <= 1e-106 ([3,1,0]; 0 or <= 3e-203 for the rational ones) and <= 4e-128 ([3,2,0]); X_B, c-hat, the Gauss parameters change, k_B, N_B, nu exponents, depths do not.
* **Phi_inf limit**: <Phi_inf(w^1,w^1')|Psi_B(eps)> ~ eps^0 on all 16 grid points on every branch, with lead / <phi_u e1|Psi_0> - 1 <= 1e-189 ... 1e-198 (24 points).
* **(vi) depth** (the largest eps-order n of a nu-coefficient nu_i^(n)(theta_s) entering the leading polynomial; equals the eps-order needed of the Laurent data (H_s, a_s [, b_s]); same for both sites):
  [3,1,0]: (6,2,0) 5; (6,1,1) 9; (5,3,0) 9; (5,2,1) 8; (4,4,0) 8; (4,3,1) 7; (4,2,2) 6; (3,3,2) 5. [3,2,0]: (6,4,0) 4; (6,3,1) 7; (6,2,2) 6; (5,5,0) 7; (5,4,1) 6; (5,3,2) 5; (4,4,2) 4; (4,3,3) 3.
  Details: the n's per nu_i (all > 1e-25 of the largest term at order N_B): the three branches of height 0 (top) are SPARSE ((6,2,0): nu_i^(n) only for n = (0,1,2,1,2,3,3,4,5) for i = 02..34, one n each, i.e. the top state is polynomial); the height-1 branches use n = 0..9 with parity patterns (e.g. (6,1,1) nu_i: i=03: n in {1,2,4,6,8}); heights >= 2 use essentially every n <= N_B. The exact series therefore must go to eps^9 (not 5 as for the adjoint) for the two height-1 branches of [3,1,0], eps^7 for [3,2,0].
  Since for lam2 = 1 b_s = a_s(H_s - a_s)/(lam1+1) only (H_s, a_s) series are needed ([3,1,0]); for [3,2,0] b_s enters through nu23, nu34 and is fixed by the no-log system (A is linear in b_s but couples the two sites).
* **Cross-check of the eps-series against W5's `ref/laurent-series-*.txt`**: not possible (files absent). Against W4's exact eps^0 data: 9e-198 / 2e-197 (above). The numerical series are in `nu_general_reg_*` (use `laurent_series(pts, bi, 24)`).

## 3. Rules

* **k_B = ht(2Lambda - mu) = n_1 + n_2 with 2Lambda - mu = n_1 alpha_1 + n_2 alpha_2**: n_1 = 2 lam1 - mu1, n_2 = mu3 - 2 lam3; for lam3 = 0: k_B = 2 lam1 - (mu1 - mu3). Holds for all 6 + 10 + 10 branches (adjoint 0,1,1,2,2,4; [3,1,0] 0,1,1,2,2,2,3,3,4,5; [3,2,0] 0,1,2,1,2,2,3,3,4,5).
* **N_B = 4 lam1 - 2 lam2 - k_B** (total eps-power of sum G nu nu), forced by <x0^G|Psi_B> ~ eps^{2 lam2 + k_B} X_B and c ~ eps^{4 lam1}.
* **The constant is 1** (no branch, rep or phi dependence) with G = (w w')^{2 lam1 - lam2} Theta, Theta = G1's kernel, nu = m/M with M = P3 p_{0,lam2+1}.
* Degree of L_B per variable <= 2 lam1 - lam2; total degree <= N_B.

## 4. Surprises / remarks

1. P4's "inhomogeneous quadratic" for [3,2,0] is the order-(lam1+2) condition B2 (quadratic in (H_s,a_s,b_s,H_t) only); the order-(lam2+1) condition A is cubic and couples both sites (the spec's wording "polynomial in (H_s,a_s,b_s) only when lam2 <= 1" is confirmed; for lam2 >= 2 it contains e_2, e_3, sigma_1, delta, and the other site's data).
2. The factorisation nu_kl = k! u_k nu_{0l} (k <= lam2 < l) means only the pairs lam2 < k < l are non-trivial functions; e.g. nu_{1l} = kappa nu_{0l} in every rep with lam2 >= 1, kappa = a/c2 = u_1 (c2 = lam2 (lam1+1)).
3. No rep-dependent constant: C = 1 for all branches (the adjoint's G^s = G_D2/4 and G1's Theta are the same normalisation).
4. 12 points are not enough for the deep branches (identity residual 1e-26 ... 1e-31 for the height-1 branches of [3,1,0]); the ratio is already visible (rounding of the spread), but I used 24 points for the final numbers. The needed depth is N_B-ish (up to 9), much deeper than the adjoint's 5.
5. Branch selection by smallest raising residual is rock-stable (0.3-0.8 eps vs >= 1.41).
6. Every pairing of a copy with W4's copy is by (smaller a_1 first): the copy-1 chain vector equals W4's copy 1 to 1e-200 in all 8 multiplicity-2 cases.

## 5. Not done

* B3 (the exact eps-series) is W5's; I produced numerical series only. No closed form of X_B(phi, theta12) (the values at phi, phi' are listed). No proof that the ν-reduction + C = 1 holds for general Lambda (only [2,1,0], [3,1,0], [3,2,0] tested; the exponent pattern k_B, N_B is a pattern from 26 branches). Multiplicity-3 ([4,2,0]) not attempted; lam3 > 0 not attempted. The general-lam2 closed form of A is given for lam2 = 1, 2 only.
