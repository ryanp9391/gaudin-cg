# RESULTS U2: the analytic untwisting limit of the Gaudin jet formula, adjoint [2,1,0]^2 (agent U2, 2026-10-16)

Files: `untwist_limit.py` (`python -X utf8 untwist_limit.py [--stages R,N,NU]`, about 1 min with caches; its full printout is `untwist_limit_output.txt`), `untwist_limit_frob.py`
(symbolic Frobenius at theta_s), caches `untwist_limit_reg_<k>.pkl` and `untwist_limit_reg2_<k>.pkl` (12 regime points for phi and for phi', 200 digits, ~1 s each). Python only; reuses `twovac_fit.py`, `qside.py`, `gaudin_chain.py` read-only.

Setting: theta = (2/7,-3/5), phi = (1/3,-2/7,5/11) (second direction phi' = (-2/5,3/8,7/11)), eps_k = 1e-6 2^-k (k = 0..11), 200 digits, 16-point grid w^1 in {3/4,-5/3,1/2,-7/3} x w^1' in {-12/5,7/2,1/8,5/4}.
Branches: smallest raising residual per HW sector; the two (3,2,1) copies ordered by the Hamiltonian eigenvalue (a: smaller). Psi_B G-normalised (Gram diagonal), sign = first nonzero GT component > 0, continuous in eps (|Psi(eps_k) - Psi(eps_{k+1})| = O(eps), no flips).
**Kernel normalisation:** the spec's G is D2's printed G divided by 4 (D2: Qnum = (1/4c) sum G N N/(M1M2)); with G^s := G_D2/4, c Qnum M1 M2 = sum G^s_ij N_i N_j exactly (used throughout).
eps-series: exact interpolation through the 12 points (x = 2^-k scaling). The 11- and 12-point fits agree to 1e-80 at the leading order of each N_i^{(n)} and to 1e-20 at lead+7, so every order used below (n <= 6) is good to better than 1e-45.
Chain-side facts used and checked: (i) the vacua are exactly eps-free (cov_inf(eps phi, eps w^) == cov_inf(phi, w^) as exact Fractions); (ii) <Phi_inf|Psi_B> = <phi_u (x) e1|Psi_B>(1+O(eps)) (closed-form Gauss parameters of plan section 0); chat = (phi2-w^1)(phi3-w^1)^2(phi2-w^1')(phi3-w^1')^2(phi1-phi2)^2.

## Headline results

1. **Order table (L1):** k_B = 0, 1, 1, 2, 2, 4 for (4,2,0), (4,1,1), (3,3,0), (3,2,1)a,b, (2,2,2); <Phi|Psi_B> ~ eps^0 at all 16 grid points on every branch; sum G N N ~ eps^{N_B}, N_B = 8,7,7,6,6,4 = 8 - k_B.
2. **Leading-polynomial identity (L2):** all six branches, both phi: L_B(w^1,w^1') = R_B * chat * <phi_u (x) e1|Psi_0> as polynomials, relative coefficient residual 1e-39 ... 1e-75, and **kappa_B = 1** with R_B = m_B/X_B (m_B = eps^-2 M1M2 at leading order, X_B = eps^-(2+k_B) <x0^G|Psi_B> at leading order), agreement 1e-41 ... 1e-63.
3. **Cleaner form (new):** dividing by M1M2 the monomials become nu = (1, kappa, pi, kappa pi, eps_s) = N/M, and **nu is an exact function of the Laurent data (H_s, a_s, b_s) of the oper at theta_s alone** (Frobenius, verified to 1e-175). The limit reads
   `sum_ij G^s_ij(eps w^) nu_i(theta1) nu_j(theta2) = eps^{6-k_B} chat(w^) <phi_u(w^) (x) e1|Psi_0> / <x(phi)(x)x(phi)|Psi_{k_B}> + O(eps^{7-k_B})`,
   verified coefficientwise to 1e-55 ... 1e-157, lower orders vanish to 1e-62 ... 1e-165. No M, no kappa_B, no Q-normalisation enter.
4. **CG dictionary:** solving the polynomial identity for the HW ket (least squares over the polynomial coefficients, unknown in H_mu, then G-normalise + round-1 sign) returns the exact CG slices; phi and phi' give identical C_nu (difference <= 2e-77).

## L1 - eps-order bookkeeping (phi; phi' has identical exponents)

| branch | k_B | <x0^G|Psi_B> ~ eps^(2+k_B) | <Phi|Psi_B> | M1M2 | N_i exponents (i=1..5, both sites) | N_B (sum G N N) | first n with <x(phi)x(phi)|Psi_n> != 0 |
|---|---|---|---|---|---|---|---|
| (4,2,0) | 0 | 2 | 0 | 2 | 1,2,1,2,3 | 8 | 0 |
| (4,1,1) | 1 | 3 | 0 | 2 | 1,1,2,2,3 | 7 | 1 (odd n only: even n vanish) |
| (3,3,0) | 1 | 3 | 0 | 2 | 1,2,1,2,3 | 7 | 1 (odd n only) |
| (3,2,1)a | 2 | 4 | 0 | 2 | 1,1,1,1,2 | 6 | 2 |
| (3,2,1)b | 2 | 4 | 0 | 2 | 1,1,1,1,2 | 6 | 2 |
| (2,2,2) | 4 | 6 | 0 | 2 | 1,1,1,1,1 | 4 | 4 (n = 0..3 vanish to 1e-74) |

* Consistency N_B = 8 + 0 - (2+k_B) + 2 = 8 - k_B: 8,7,7,6,6,4, all six OK.
* **Interpretation of k_B.** On the exact HW kets of every non-top mu, <x(phi)(x)x(phi)|v> = 0 identically (<= 4e-202 at 200 digits, both phi, all of H_mu, both octet copies); for (4,2,0) it equals (phi2-phi1)^2 = 0.38322.
  So <x0^G|Psi_B> gets its first contribution from the k_B-th order vector: X_B = <x(phi)x(phi)|Psi_{k_B}>. Values (phi): 0.383220, -0.017394, -0.074288, 0.0013930 (both copies), 6.9931e-6; (phi'): 0.600625, 0.065737, 0.145765, 0.0065907, 9.9880e-5.
  X_B is expected homogeneous of degree 2+k_B in phi (not independently fitted). The singlet needs Psi to O(eps^4), the octets O(eps^2): this is the required depth of the degenerate perturbation theory.
* The two copies in (3,2,1) have identical exponents, identical X_B and identical R_B.

## L2 - leading polynomial and the identity

(a) N_i^{(n)}(theta_s) extracted for 6 branches x 2 sites x 5 monomials (leading exponents in the table above); example (4,2,0), site 1: N_1 = 0.655659 eps + 0.116778 eps^2 + 0.0103995 eps^3 + ...
(b) L_B(w^) = sum_{d+n+m = N_B} G^s_ij^{[d]} N_i^{(n)}(theta1) N_j^{(m)}(theta2). All lower orders cancel: max coefficient of R_t at t = N_B - 1 is 4.5e-51 (4,2,0), 1.2e-55 (4,1,1), 1.1e-55 (3,3,0), 6.6e-66 (3,2,1)a/b, 2.9e-84 (2,2,2) (and 1e-103 ... 1e-60 below), versus 0.2-0.4 at t = N_B. (The fall-off at lower t is the fit's truncation, not physics.)
**Depth (contributing (d,n,m) at order N_B, terms > 1e-30):** (4,2,0): d in 2..6, n,m in 1..4 (largest (4,2,2), (3,2,3), (5,1,2), (4,1,3)); (4,1,1) and (3,3,0): d in 1..5, n,m <= 5 (largest (3,2,2), (4,1,2), (2,2,3)); (3,2,1)a,b: d in 0..4, n,m <= 4 ((2,2,2), (3,1,2), (2,1,3), (0,3,3), (0,2,4), ...);
(2,2,2): only six triples (1,1,2), (2,1,1), (1,2,1), (0,2,2), (0,1,3), (0,3,1). Full per-branch lists are in `untwist_limit_output.txt`.
(c) Identity coefficient by coefficient, as polynomials in (w^1,w^1'). chat*[e1^T rho(n')]_nu are polynomials for every component that occurs in a slice (see Surprises 1). Residual |L_B - R_B chat GF|/|L_B| (GF built from the chain's own Psi_0, each octet copy separately), R_B = best-fit constant, compared with m_B/X_B from the series:

| branch | R_B (phi) | residual (phi) | rel. diff to m_B/X_B (phi) | R_B (phi') | residual (phi') | rel. diff (phi') |
|---|---|---|---|---|---|---|
| (4,2,0) | -1.1217793595 | 2.4e-43 | 6e-46 | -0.0892224658 | 7.7e-39 | 4e-41 |
| (4,1,1) | -4.8470775935 | 2.7e-47 | 3e-57 | 0.1598811435 | 2.8e-46 | 6e-56 |
| (3,3,0) | -1.1349116102 | 1.8e-47 | 3e-57 | 0.0721026407 | 1.6e-45 | 8e-55 |
| (3,2,1)a | -7.5969660751 | 2.9e-57 | 1e-63 | -0.2001611630 | 2.1e-54 | 3e-63 |
| (3,2,1)b | -7.5969660751 | 7.3e-58 | 1e-63 | -0.2001611630 | 3.4e-55 | 6e-63 |
| (2,2,2) | -90.949002471 | 7.0e-75 | 4e-45 | -0.7937987684 | 3.8e-72 | 3e-44 |

So **kappa_B = 1**, with M-factor := m_B = lim eps^-2 M1M2 and X_B as above: L_B = chat * <phi_u (x) e1|Psi_0> * m_B/X_B. The two octet copies pair with their own HW ket (different generating functions, same constant R_B).
For the top state m_B = -(phi3-phi2)^2 theta12^2 (Q_1 = 1) to 13 digits (-0.42988823527), so R_B(top) = -(phi3-phi2)^2 theta12^2/(phi2-phi1)^2. In general M_s = P3 A_s Q_1 with A_s = (z3-z2) e^{sum z theta_s}(theta_s-theta_t), so m_B = -(phi3-phi2)^2 theta12^2 Q_1^0(theta1) Q_1^0(theta2) (not separately checked beyond the top).
(d) Normalisation-free check: rho_B(w^;w^0) = [chat(w^0)/chat(w^)] L(w^)/L(w^0), w^0 = (3/4,-12/5), at all 16 grid points vs the exact ratio from <phi_u (x) e1|Psi_0>: max relative deviation (nu-form L'), branches in the usual order:
phi: 2.5e-153, 3.1e-52, 4.4e-51, 3.0e-64, 1.2e-64, 3.4e-82; phi': 1.2e-152, 2.5e-57, 3.8e-51, 8.3e-64, 3.5e-62, 1.2e-80.

**CG dictionary** (solve L/chat in the span of chat*[e1^T rho(n')]_nu over the HW space H_mu, G-normalise, round-1 sign; slice C_nu = <e_nu (x) e1|Psi_0>, gl3_reps index nu):

| branch | recovered C_nu | exact (round 1 / W1) | max |v_Q - Psi_0| over all 64 components |
|---|---|---|---|
| (4,2,0) | C_7 = 1 | 1 | 0 |
| (4,1,1) | C_5 = -0.70710678119 | -1/sqrt2 | 6e-100 |
| (3,3,0) | C_6 = -0.70710678119 | -1/sqrt2 | 7e-101 |
| (3,2,1)a | C_2 = 0.208887289734, C_4 = 0.049311600013 | (0.0493116000, 0.2088872897) for (a=4, a=2) | 3e-77 (nu-form), 1e-66 (N-form) |
| (3,2,1)b | C_2 = -0.0797878448606, C_4 = 0.337986734608 | (0.3379867346, -0.0797878449) for (a=4, a=2) | 3e-77, 2e-68 |
| (2,2,2) | C_0 = 0.117851130198 | 0.1178511302 | 3e-102 |

phi and phi' give the same C_nu: max differences 0, 0, 0, 2e-77, 1e-77, 1e-202 (branch order as above). The exact HW slices in H_mu for (3,2,1) are (C_2,C_4) proportional to (-1/2,3/2) and (1/2,1/2); the two Q-side vectors are the two chain copies.

## L3 - the formula in terms of Q-functions and what it needs

**Reduction of the jet monomials.** At theta_s the dual solutions are P^j = A_j f_0 + B_j f_2 + C_j f_4 (Frobenius exponents 0, 2, 4); to derivative order 3 only f_0 = 1 + kappa x + c2 x^2 + c3 x^3 and f_2 = x^2 + d3 x^3 matter, so
p_kl(theta_s) = [A,B]_23 w_kl with w_02 = 2, w_03 = 6 d3, w_12 = 2 kappa, w_23 = 12 (c2 d3 - c3), w_01 = 0. Hence N = M (1, kappa, pi, kappa pi, eps_s) with
kappa_s = P'/P(theta_s) (any dual Q), pi_s = p03/p02 = 3 d3, eps_s = p23/p02 = 6(c2 d3 - c3), and **M cancels** from c Qnum M1 M2 = sum G N N:

  **c Qnum = sum_ij G^s_ij(omega1, omega1', theta12) nu_i(theta1) nu_j(theta2),**  Qnum = <Phi_inf|Psi>/<x0^G|Psi>.

The Frobenius recursion with S1 = tr zeta + 3/x1 + 3/x2, S2 = e2(zeta) + sum_s [H_s/x_s + 3/x_s^2], S3 = e3(zeta) + sum_s [b_s/x_s + a_s/x_s^2] (qside conventions, delta = theta_s - theta_t, sigma_1 = zeta1+zeta2+zeta3) gives exactly (sympy, `untwist_limit_frob.py`):

* kappa_s = a_s/3,
* pi_s = 2 H_s - a_s - 2 sigma_1 - 6/delta,
* eps_s = -(2/3)[ H_t a_s/delta - a_s b_s + a_s e2 + 3 a_s/delta^2 - 3 a_t/delta^2 - 3 b_t/delta - 3 e3 ],
* the m = 2 no-log (Bethe) condition H_s a_s - a_s^2 - 3 b_s = 0, i.e. b_s = a_s (H_s - a_s)/3 (verified to 1e-199 on all branches).

Verified: nu(H,a,b) equals N_i/M built from the Q's to 1e-174 ... 1e-193 on all six branches. **So the Q-number is a function of (H_1, a_1, H_2, a_2), zeta, theta12, omega only**: the Q polynomials, their monic normalisation, m_B and kappa_B never enter.
Exact side data: H_1 + H_2 = sum_{i != j} zeta_i mu_j (checked on the eps-coefficients, e.g. (3,2,1): 2.12987); H_1^{(0)} = -H_2^{(0)} = [9 - (C_2(mu)-18)/2]/theta12 = 4, 6, 6, 9, 12 (over theta12 = 31/35) for (4,2,0), (4,1,1)/(3,3,0), (3,2,1), (2,2,2) (checked: 4.516129, 6.774194, 10.16129, 13.548387);
a^{(0)} = 0 ((4,2,0), (3,3,0)), = H^{(0)} ((4,1,1), (2,2,2)), and for the octet copies (a_1,a_2)^{(0)} = (4.680852619, -12.25463125) (a) and (12.25463125, -4.680852619) (b). The (4,2,0) series are exactly polynomial: H_s = +-140/31 + 1.12554 eps, a_s = 1.36364 eps, b_{1,2} = +-2.05279 eps - 0.108225 eps^2.

**Depth needed (nu-form, order 6 - k_B).** (i, n) pairs entering the leading polynomial, i = 0..4 labelling (1, kappa, pi, kappa pi, eps_s), n = eps-order of nu_i, over all contributing G-pieces (terms > 1e-30):

| branch | nu exponents | n needed per site | max n |
|---|---|---|---|
| (4,2,0) | 0,1,0,1,2 | only the leading coefficient of each: (0,0),(1,1),(2,1),(3,2),(4,3) | 3 |
| (4,1,1) | 0,0,1,1,2 | i=1: 0,1,2,4; i=2: 1,2,4; i=3: 1,2,3,5; i=4: 2,3,4 | 5 |
| (3,3,0) | 0,1,0,1,2 | i=1: 1; i=2: 0,1,2,4; i=3: 1,2,3,5; i=4: 2,3,4 | 5 |
| (3,2,1)a,b | 0,0,0,0,1 | i=1..3: 0..4, i=4: 1..4 | 4 |
| (2,2,2) | 0,0,0,0,0 | i=1..4: 0..2 | 2 |

In terms of the Laurent data this means eps-orders up to eps^5 ((4,1,1),(3,3,0)), eps^4 (octets), eps^3 (top), eps^2 (singlet) of (H_s, a_s), much shallower than N_B = 8 for the unreduced monomials. The non-top series are genuinely infinite (e.g. (3,3,0): H_1 = 6.774194 + 1.435065 eps + 0.042428 eps^2 + 0 eps^3 - 0.000797 eps^4 + ...), so a closed form of the limit needs the exact perturbation theory of the Hamiltonians to those orders; I did not attempt it.

**Explicit generating function** (Psi_0 G-norm 1):

  <phi_u(w^1,w^1') (x) e1|Psi_0> = X_B lim_{eps->0} eps^-(6-k_B) [sum_ij G^s_ij(eps w^) nu_i(theta1;eps) nu_j(theta2;eps)] / chat(phi,w^1,w^1'),  X_B = <x(phi)x(phi)|Psi_{k_B}>.

X_B is an overall constant (drops out of every ratio and of the dictionary, where the vector is fixed by the polynomial identity up to scale and then G-normalised). In "CG from Q" form: C_nu are the coefficients of
[sum_ij sum_{d+n+m = 6-k_B} G^s_ij^{[d]}(w^1,w^1',theta12) nu_i^{(n)}(theta1) nu_j^{(m)}(theta2)] / chat in the basis chat*[e1^T rho(n'(w^))]_nu (the polynomials `gf_polys`), up to one overall constant (fixed by the G-norm / sign convention).
**The final formula needs exactly:** (i) the kernel G (the 25 polynomials of `RESULTS-twovac-limit.md`, G^s = G_D2/4) split by (w1,w1')-degree; (ii) the one-site numbers nu(theta_s; eps) = (1, a/3, pi, a pi/3, eps_s), i.e. (H_s, a_s) per site and branch as eps-series up to the orders in the table
(b_s from the no-log relation); (iii) the Gauss parametrisation (a',g',c') of phi_u (plan section 0) and chat. The Q polynomials are not needed beyond these numbers (kappa = P'/P, pi, eps_s are Plucker ratios at theta_s).
U1's `untwisted_Q_series_210.pkl` / `RESULTS-untwisted-Q.md` and W4's `ref/untwisted-exact-*.txt` were not present when I finished, so the eps-series of the Laurent data are my own numerical ones (200 digits); the formula is stated in the Laurent-data form precisely so that U1's/W4's exact series can be dropped in.

**phi-independence.** Repeated L1-L2 and the dictionary with phi' = (-2/5, 3/8, 7/11): N^{(n)}, chat, the Gauss parameters, R_B, X_B and m_B change (tables above); the recovered C_nu are identical (max difference <= 2e-77; 0 for the 1-dim sectors), and the exponents, k_B, N_B and depths are unchanged.

## L4 - not done

[3,1,0]/[3,2,0] need P3/P4's regime machinery and polynomial laws of the kernel (not available for [3,1,0]); skipped. The nu-form should generalise directly (non-forced Plucker pairs, Frobenius exponents {0, lambda2+1, lambda1+2} at theta_s, p_kl as minors of the local coefficients), giving the kernel contraction in terms of Laurent data; not attempted.

## Surprises / remarks

1. chat does not clear all components of e1^T rho(n'): the weight-(0,2,1) component -2 a' g' (nu = 3) carries (phi2-w^1')^-2. Harmless: C_3 = <e_3 (x) e1|Psi_0> = 0 on every HW ket of the sectors involved (asserted < 1e-40 in the code).
2. kappa_B = 1 exactly, with the M-factor defined as m_B = leading coefficient of M1M2/eps^2; the whole difficulty of the monic normalisation sits in m_B and disappears in the nu-form.
3. The nu-form improves the identity by 14+ digits (1e-55..1e-157 vs 1e-43..1e-75) and needs eps-orders <= 5 instead of 8.
4. <x(phi)x(phi)|Psi_n> has definite parity in n for (4,1,1), (3,3,0) (odd n only); for (2,2,2) the first four orders vanish identically.
5. a_s^{(0)} separates the two octet copies at zeroth order; their H_s^{(0)} coincide.
6. Branch selection by smallest raising residual is stable for all 12 eps (1e-7..1e-10 vs gap 1.4-2.1); the (3,2,1) copies keep their order (a below b).

## Not done

L4; exact Casimir/closed-form eps-series of the Laurent data beyond the top branch and zeroth orders (U1/W4); use of U1's pkl and W4's exact tables (absent at finish; I used the exact HW kets of `twovac_fit.hw_space`); closed form of X_B in (phi, theta12); independent check of m_B = -(phi3-phi2)^2 theta12^2 Q_1^0 Q_1^0 beyond the top state.
