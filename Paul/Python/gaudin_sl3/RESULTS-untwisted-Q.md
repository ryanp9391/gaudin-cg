# RESULTS-untwisted-Q (agent U1, 2026-10-16): the untwisted Q-system and the eps-series of the HW-limit branches (two-site sl(3) Gaudin, Lambda = (2,1,0))

Code `untwisted_Q.py` (`python -X utf8 untwisted_Q.py --stages T0,T1,EXPORT`; reuses `gaudin_chain.py` read-only; no Wolfram). Units below: theta1 - theta2 = 1 ("units"); at general theta
H_s, a_s scale as theta12^(n-1) at order eps^n (zeta = eps*phi), b_s as theta12^(n-2) (verified at two theta points; W4's files `ref/untwisted-exact-*.txt` agree exactly with all T0 data,
copy 1 = smaller alpha_1 = my copy 0). Caches: `untwisted_Q_reg*.pkl` (branch data/Q's per eps), `untwisted_Q_laurent_*.pkl`, `untwisted_Q_octet_fit.pkl`, `untwisted_Q_pert_sym.pkl`.

## Summary of what changed in the plan (per the main session's redirect)
Priority was moved to the eps-series of the Laurent data (H_s, a_s; b_s = a_s(H_s - a_s)/3). Result: **all six branches have EXACT Laurent-data series, polynomial in phi** (homogeneous of degree n at order eps^n,
coefficients rational, or in Q(sqrt5) for the octets); top state, (4,1,1), (3,3,0) in closed algebraic form. U2's nu-formulas were re-verified independently (below).

## T0: untwisted Laurent data (exact, Fractions on exact HW kets), five reps all sectors
Operator identities (global, on the whole V x V, exact): H1 + H2 = 0, b1 + b2 = 0, a1 + a2 + th1 b1 + th2 b2 = 0 (reps 110, 210, 220, 310, 320; two theta points each). So S2 = O(u^-2), S3 = O(u^-3): Fuchsian at infinity.
Closed forms (verified in every sector of all five reps, eta = H1 th12 etc.):
* H1 = (|L|^2 - (C2(mu) - 2 C2(L))/2)/th12, C2(l) = sum l_i(l_i+4-2i), the E_ik E_ki normalisation (|L| = sum L_i). Adjoint: 4, 6, 6, 9, 12 for (4,2,0),(4,1,1),(3,3,0),(3,2,1),(2,2,2); H2 = -H1.
* (a1 - a2)/2 = q/th12 with **q = (mu1+2)(mu2+1) mu3 / 2** (found by an exact fit over 23 sectors; in terms of C3 = sum E_ij E_jk E_ki, eigenvalue p3 - 3p2/2 - p1^2/2 + p1 with l = (mu1+2, mu2+1, mu3),
  p_k = sum l_i^k; equivalently q = C3/6 - (|L|+1) C2(mu)/2 + (2/3)|L|(|L|+1)^2.) q = 0, 6, 0, 15/2, 12 for the adjoint sectors; q = 0 whenever mu3 = 0.
* (a1 + a2) = -b1 th12 (sum rule), b2 = -b1. b1 = 0 in multiplicity-1 sectors. Multiplicity-2 sectors (the only ones occurring): b1 th12^2 = +-sqrt(beta), beta = 45 [210: (3,2,1)], 105, 192 [310: (5,2,1),(4,3,1)], 105, 192 [320: (5,4,1),(5,3,2)];
  a1 th12 = q -+ sqrt(beta)/2: adjoint (15 -+ 3 sqrt5)/2 (= universal matrix (3/(2 th12))[[3,-1],[-1,7]] eigenvalues). No closed form for beta in (Lambda, mu) found (not a function of C2, C3 alone: needs the multiplicity-space operator).
* Sum rules: b1+b2 = 0 and a1+a2+th1 b1+th2 b2 = 0 hold identically => untwisted oper Fuchsian at infinity.
* C3 eigenvalues (adjoint): 144, 132, 96, 69, 24. Both Casimirs computed exactly on each HW ket.

## T1: the untwisted oper
* Exponents: {0,2,4} at theta_s (indicial rho(rho-2)(rho-4), S3 has no triple pole); at infinity solutions grow like u^k with k^3 - 9k^2 + (14 + h)k - 2q = 0 (h = H1 th12): **k = (mu3, mu2+1, mu1+2)** in all 6 branches
  (e.g. top (0,3,6), (4,1,1): (1,2,6), singlet (2,3,4)). Fuchs relation: 12 - sum k = 12 - (|mu|+3) = 3 OK.
* **Refutation of the expected "fewer than 3 polynomial solutions": all three solutions are polynomials** (no-log at theta_s => trivial local monodromy there, integer exponents at infinity => trivial monodromy), of degrees (mu3, mu2+1, mu1+2),
  Wronskian c (u-th1)^3 (u-th2)^3. Exact dimension of polynomial solutions of degree <= d (checked d <= 8): top 1,1,1,2,2,2,3; (4,1,1) 0,1,2,2,2,2,3; etc., always 3.
* Exact echelon bases V0 (monic, x = (u-th2)/th12; site 2 at x=0, site 1 at x=1), all verified to solve D0 (sympy):
  top (4,2,0): {1, x^2(2x-3)/2, x^4(2x-3)^2/4} = {1, F, F^2} with F = int Pi (a Sym^2); (4,1,1): {(2x-1)/2, x^2, x^4(2x^2-6x+5)/2}; (3,3,0): {1, x^2(x-1)^2, x^2(4x^3-10x+5)/4};
  (2,2,2): {(2x-1)^2/4, (2x-1)(4x^2+2x-1)/8, x^4}; (3,2,1) copy 0/1: {(10x-5+sqrt5)/10, x^2(2x-3-sqrt5)/2, x^4(2x-5+sqrt5)/2} (copy 1: sqrt5 -> -sqrt5) (see `v0_bases()`; copies are Galois conjugate).
* Hypergeometric verdict: **not a 3F2 in general.** Only when the exponents at some singular point have the pattern {0,1,s} up to a common shift: infinity of (4,1,1) (1,2,6), (3,3,0) (0,4,5), (2,2,2) (2,3,4) -- for these three the oper is exactly
  (Moebius X = (u-th2)/(u-th1) and gauge (u-th1)^c, c = 2, 5, 3 (or 4)) the 3F2 equation with parameters a = (rho - c), b = (-1,-3) (degenerate/reducible integer parameters; checked by applying the operator to the V0 basis). For top (0,3,6) and the octets (1,3,5) no singular point has that pattern. The ode is a trivial-monodromy Fuchsian equation, so "hypergeometric" is not the natural frame; the polynomial bases above are the closed forms.
* **eps -> 0 limit of the branch Q's (the key structural finding).** The monic p_j(u;eps) do NOT have finite coefficient-wise limits: coefficient c_m of p_j ~ eps^(-min(mu_j - mu3, mu_j - m)) (top: P1 coefficients eps^-4..0, P2 eps^-2..0, P3 finite; (4,1,1): P1 eps^-3, P2,P3 finite; (3,3,0): P1,P2 eps^-3; (3,2,1): P1 eps^-2, P2 eps^-1; singlet finite).
  Expand P^j in the exact V0 basis (v_3,v_2,v_1 of degrees (mu3, mu2+1, mu1+2)), P^j = sum_a M_ja(eps) v_a(eps) (v_a(eps): solutions of D_eps with the same jets at u = 0): P^3 -> v_3 (finite), P^2, P^1 have pole order mu_j - mu3 on v_3 and
  O(eps^-(N_j-1)) on v_2, v_1; det M = c_W Pi(0)^3 ~ eps^3 exactly. Normalised rows (1, rho_2, rho_1)_j have rho = eps r^(1) + ... with **the first order direction common to all j** (and its scale tracking phi3 in the top state), the rows
  differing only at order eps^2 (octets, singlet) .. eps^4 (top: R1-R3 = O(eps^4), R2-R3 = O(eps^3)): the limiting 3-plane is V0 and the zeta-jet basis is the Smith/saturation basis of this matrix (elementary-divisor valuations were not computed -- script prepared but not run). No closed form for r^(1) obtained, only the observed structure.
* Limits of differences/derivatives: leading Laurent coefficients of P^j solve D0 (in V0); lower ones solve D0 p^(n+1) = -D1 p^(n) - ... (inhomogeneous), consistent with the numerics.

## T2: eps-series
Series fit of f(eps_k), eps_k = 1e-6 2^-k: (a) K = 12, 200 digits, n <= 10 [`untwisted_Q_series_210.pkl`]; (b) extended K = 24, 400 digits, n <= 20 [`untwisted_Q_series_210_ext.pkl`]; (c) second direction phi' = (-2/5,3/8,7/11), K = 12, 200 digits [`..._phi2.pkl`].
Stability = relative difference of fits dropping the first/last point: order n reliable to 1e-90 (n=1) ... 1e-16 (n=9), ~1e-8 (n=10) for K=12; for the extended set rel.stab 1e-239 (n=1), 1e-62 (n=16), 1e-25 (n=20): **reliable depth n <= 18 (1e-34), <= 20 (1e-14)**.
Low orders vanish to the noise floor (~1e-100 at 200 digits) and are set to exactly 0 below `n0`.
First nonzero orders n0 of N_i = P3 p_kl [pairs (02),(12),(03),(13),(23)], identical on both sites: top 1,2,1,2,3; (4,1,1) 1,1,2,2,3; (3,3,0) 1,2,1,2,3; (3,2,1) (both copies) 1,1,1,1,2; (2,2,2): 1,1,1,1,1 (site 1) [see pkl]; M = P3 p02 ~ eps^1; p01 = 0 (1e-195, Frobenius).
**pkl format** (all three files the same): dict(meta=dict(phi, theta, dps, K, eps, nmax, branches=[(mu,copy)...]), branches={(mu,copy): dict(N, M, P3jets, plucker, data, p)}) where
 N[(s,(k,l))], M[s], P3jets[(s,k)] (k = 0..4, jets of P^3 at theta_s), plucker[(s,(k,l))] (0<=k<l<=3), data[nm] (nm in H1,a1,b1,H2,a2,b2), p[(j,m)] (coefficient of u^m of p_j, j = 0,1,2 for P^1,P^2,P^3; Laurent: coefficient series of eps^pole * c_m)
 are each dict(n0 = first nonzero order, coef = list of series coefficients n = 0..nmax as 200-digit strings, relstab = per-order relative stability, n_reliable = largest n with relstab < 1e-10 (pole shifts for p)); s = 0,1 sites.
 Branch labels: (4,2,0),(4,1,1),(3,3,0),(3,2,1) copy 0 (smaller a1^(0): (15-3sqrt5)/2 th12^-1 ... a1 = 4.6808526 at th12 = 31/35),(3,2,1) copy 1,(2,2,2).
* **U2's nu-formulas verified independently**: N_i = M nu_i with nu = (1, kappa, pi, kappa pi, eps_s), kappa = a_s/3, pi = 2H_s - a_s - 2 sigma1 - 6/delta, eps_s = -(2/3)[H_t a_s/delta - a_s b_s + a_s e2 + 3 a_s/delta^2 - 3a_t/delta^2 - 3 b_t/delta - 3 e3] holds to 1e-185 for all 6 branches, both sites, at eps = 1e-6/8. (Not needing the Q polynomials at all.)
* The no-log conditions (derived and verified to 1e-148): per site, (1) 3b = a(H - a); (2),(3) two more involving the other site's Taylor data (`frobenius_conditions()`), 6 equations on (H_s,a_s,b_s) plus the two exact global rules
  H1 + H2 = eps sum_{i != j} phi_i mu_j and b1 + b2 = eps^2 (phi1 phi2 mu3 + phi1 phi3 mu2 + phi2 phi3 mu1) (checked at all orders; first reduces to U2's rule).

## T2(b)/T3: EXACT Laurent-data series (units theta12 = 1; zeta = eps phi; b_s = a_s(H_s-a_s)/3; H2 = eps sum_{i!=j} phi_i mu_j - H1 [all branches]; also a2 listed in pkl)
Method: order-by-order solution of the no-log system + the two global rules (Jacobian of 8 equations on 6 unknowns of rank 6 for (4,1,1),(3,3,0),(2,2,2): exact symbolic phi to n = 6, `perturb_data`); for the top state and octets the Jacobian has rank 5
(kernel dimension 1; my degenerate-solvability implementation failed -- cubic with no rational root, not resolved), so top = exact closed form from the Wronskian structure, octets = polynomial fit in phi on 29 generic phis (200-300 digits) followed by exact Q(sqrt5) reconstruction (PSLQ), n <= 5.
All series validated against the 200/400-digit numerics at phi = (1/3,-2/7,5/11): agreement 1e-60 (n <= 5, all branches, H1, a1, b1). They are polynomial in phi, so rational in phi in the weak sense (no denominators); non-rational only through sqrt5 at the octets.
* (4,2,0) top (exact, all orders): H1 = 4 + eps(phi1 + 2 phi2 + 3 phi3); a1 = a2 = 3 eps phi3 (= 3 zeta3: kappa = zeta3, P^3 = e^{zeta3 u}); b_s = a_s(H_s-a_s)/3. At phi: H1 = 4 + (260/231) eps in units.
* (4,1,1), with x = eps(phi2 - phi3): H1 = 6 + eps(2phi1 + 5phi2 + 5phi3)/2 + (sqrt(1 + x^2/4) - 1), a1 = 6 + (3/2) eps(phi2 + phi3) + 3(sqrt(1 + x^2/4) - 1) (n <= 6: x^2/8 - x^4/128 + x^6/1024 resp. 3x^2/8, ...; odd orders n >= 3 vanish).
* (3,3,0), y = eps(phi1 - phi2): H1 = 6 + (3/2) eps(phi1 + phi2 + 2phi3) + (sqrt(1 + y^2/4) - 1), a1 = 3 eps phi3 (exact), a2 = 3 eps phi3 + ... (a1 = a2).
* (2,2,2): with Q = phi1^2 + phi2^2 + phi3^2 - phi1phi2 - phi1phi3 - phi2phi3, D = (phi1-2phi2+phi3)(phi1+phi2-2phi3)(2phi1-phi2-phi3), s1 = phi1+phi2+phi3 (all in units, eps powers implied):
  H1 = 12 + 2 s1 e + Q e^2/6 + 0 e^3 - Q^2 e^4/288 + 0 e^5 + (degree-6 poly) e^6 ; a1 = 12 + s1 e + Q e^2/4 + D e^3/72 - Q^2 e^4/192 - (11/10368) D Q e^5 + (degree-6 poly) e^6 (poly printed in `untwisted_Q_pert_sym.pkl`; same sextic for H1 up to factor 2/3 relation: H1 e^6 coefficient = (2/3) a1's).
* (3,2,1) copy 0 (sqrt5 -> -sqrt5 for copy 1), orders in phi, units theta12 = 1:
  H1^(0) = 9, H1^(1) = [(15+sqrt5) p1 + (20 - 2 sqrt5) p2 + (25 + sqrt5) p3]/10, H1^(2) = 2(2p1^2 - 3p1p2 - p1p3 + 3p2^2 - 3p2p3 + 2p3^2)/75 (rational), H1^(3) = -sqrt5 (p1-p3)^2 (p1 - 2p2 + p3)/375,
  H1^(4) = -2 (p1^2 - 4p1p2 + 2p1p3 + 9p2^2 - 14p2p3 + 6p3^2)(6p1^2 - 14p1p2 + 2p1p3 + 9p2^2 - 4p2p3 + p3^2)/84375, H1^(5) = sqrt5 (p1 - 2p2 + p3) (quartic)/6328125 (printed in the pkl / `t2o` log);
  a1^(0) = (15 - 3sqrt5)/2, a1^(1) = [(3 - sqrt5) p1 + (9 - 3sqrt5) p2 + (18 + 4 sqrt5) p3]/10, a1^(2), ..., a1^(5) and b1^(1..5) exact in `untwisted_Q_octet_fit.pkl` (b1^(0) = 3 sqrt5, b1^(1) = (12p1+11p2+22p3)/5, H1^(1)+... sum rules hold).
  Even orders H^(2n) are rational (Galois-invariant), odd orders carry sqrt5; the two copies are Galois conjugates.
* Exact series (n <= 5) for all six branches: `untwisted_Q_laurent_exact.pkl` (dict[(branch, 'H1'|'a1'|'b1')] -> list of sympy strings in p1,p2,p3, orders 0..5).
Orders reached: exact symbolic in phi: top all orders, (4,1,1),(3,3,0),(2,2,2) n <= 6 (algebraic for the first two), octets n <= 5 (limited only by the number of phis in the fit grid; n = 6 needs >= 30 phis).

## T3
* The nu-vector is a function of the Laurent data only (verified), so nu^(n) follow from the series above: nu_2 = kappa = a_s/3, pi = 2H_s - a_s - 2 sigma1 - 6/delta (sigma1 = eps s1), eps_s via H_t, a_s, b_s, a_t, b_t, e2, e3 = eps^2, eps^3 times symmetric polynomials. Closed forms for the three algebraic branches follow by substitution of the sqrt expressions.
  Example (top): kappa = eps phi3; pi = 2(4 + eps(phi1+2phi2+3phi3)) - 3 eps phi3 - 2 eps s1 - 6/delta with delta = +-1.
* Which combinations are phi-independent: the zeta = 0 data (T0) and V0; the first-order direction of the jet rows (up to scale) per branch; the leading orders of N_i are rational functions of phi with poles from (phi_j - phi_k)^-1 (top: P^2's coefficients ~ (zeta3-zeta2)^-k).
* Top state in closed form: P^3 = e^{zeta3 u}, P^2 = e^{zeta2 u}(Pi + Pi'/d + Pi''/d^2), d = zeta3 - zeta2 (Wronskian relation; checked consistent with the 200-digit numerics).

## Surprises
1. All three untwisted solutions are polynomials (degrees mu3, mu2+1, mu1+2); the "fewer than 3 polynomial solutions" expectation is refuted. 2. Exponents at infinity are exactly (mu_i + 3 - i) reordered: the Weyl-shifted weight. 3. q = (mu1+2)(mu2+1)mu3/2 (a_s is a cubic-Casimir quantity = the "determinant" of l). 4. The monic p_j have poles in eps (coefficient-wise limit does not exist); finite objects are the Plucker/N_i monomials and the Laurent data. 5. (4,1,1),(3,3,0) Laurent data are algebraic (sqrt(1+(zeta_j - zeta_k)^2 theta12^2/4)): the Bethe equations solve in closed form there.

## Not done / limits
* Degenerate perturbation theory (rank-5 Jacobian) for the top and octets symbolically (replaced by closed form / fit; octet series beyond n = 5 not produced; singlet beyond n = 6; no closed algebraic form identified for (2,2,2) or octets).
* ZETA-jet basis not given in closed form (only structure, valuations and directions); beta(Lambda,mu) closed form; reps other than (2,1,0) for T1/T2 (T0 is done for 110, 220, 310, 320).
* The exact p_j^(n) polynomial series (T2 exact orders of the Q polynomials) were replaced by the Laurent-data route (per the redirect); the numerical Laurent series of p_j are in the pkl.
