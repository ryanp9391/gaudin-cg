# RESULTS-gaudin-norm (agent Q1, round 5 WP-D, 2026-10-16): the Gaudin norm from the Q side

Files: `gaudin_norm.py` (`python -X utf8 gaudin_norm.py --stages D0,D1,G,D2`, ~1-2 min with caches), `gaudin_norm_output.txt` (full printout), caches `gaudin_norm_D0_<tag>.pkl` (per-state table: sector, Laurent data, m_i(theta_s), M_s, det N_B, Rv, x0, norm, N, N_Omega, det J6), `gaudin_norm_c0.pkl`, `gaudin_norm_G1sym.pkl`.
Python only; reuses twovac_fit, gaudin_B, untwisted_Q, untwist_limit read-only. Adjoint [2,1,0]^2, theta = (2/7,-3/5) unless stated.

## 0. Ground truth and conventions
N(Psi) := Psi^T (G(x)G) Psi / <x0^G|Psi>^2 (G = gram_diag; x0^G = x(zeta)(x)x(zeta) = x0_D2; x0_onesite gives the same). N > 0 on all 64 states at P, Pp, X1 (ranges e.g. P: 0.62 ... 1e9).
With the twisted vacuum <Omega_0 (x) Omega_0| (Rv of detdata) instead of x0: N_Omega = N (x0/Rv)^2, and **(Rv/x0) M1 M2 / det N_B is state independent** (P: 12.5381319888976, Pp: 0.76580951791315, X1: 0.0056553714130498; spread 1e-74): D3's identity, so N_Omega is N times a det-N_B factor, nothing new. Per-state table of Laurent data, monomials, M_s, det N_B in the D0 pickles.

## 1. Candidates
| candidate | result |
|---|---|
| (i) det N_B (Lagrange columns), alone or with M1, M2 | FAIL (log-regression residual 22 of 37 spread) |
| (iii) Q-values at theta_s (P^j, Q_j, M_s) times det N_B, 15-17 free exponents | FAIL (residual 2.8) |
| **(ii) Jacobian of the six no-log conditions** | **PASS** |
| (iv) SoV measure Gamma | PASS by construction (section 2) |

**(ii) result (exact, found by regression then verified to 73 digits):** with J6 = det d(F_1..F_6)/d(H_1,a_1,b_1,H_2,a_2,b_2) (the 3 no-log conditions per site of `untwisted_Q.frobenius_conditions`, the other site entering through the regular parts; implemented as `jac_funcs`, F = 0 on all states to 1e-79)

    N(Psi) = (9/4) |det J6| (zeta_1-zeta_2)^(mu_3-6) (zeta_1-zeta_3)^(mu_2-6) (zeta_2-zeta_3)^(mu_1-6),   mu = dual-Q degrees (sector), |.| on the zeta differences

Equivalently N/(|J6| d23^mu1 d13^mu2 d12^mu3) = c0 = (9/4)/(d12 d13 d23)^6. State independent to 1e-73 over all 64 states at P, Pp, X1; **independent of theta** (5 thetas at P give the identical c0 = 698440.12045259454819...); c0 vs zeta checked at 34 points: |c0 (d12 d13 d23)^6/(9/4) - 1| <= 1.7e-73. The same holds for the reduced 4x4 (b eliminated; dF_1/db = -1). No M_s, P^3, Wronskian or det N_B factor is needed. This is the Gaudin "Yang-Yang Hessian in Laurent variables". Only open point: det J6 has no definite sign (32 positive, 32 negative at every zeta; not explained by the sign of any Q-value, Laurent datum, nu, difference thereof; GF(2) searches failed), so the formula uses |det J6| (N > 0 forces it).

## 2. Gamma, the SoV measure
Setup: Psi_c = <x0|Psi> sum_ij nu_i(theta_1) nu_j(theta_2) |ij>, |ij> = columns of Y^-1 (D3's zeta-independent exact Y, `comp_jet_matrix`), nu = m/M. All 8 one-site nu follow from Frobenius: **nu = (1, k, p, k p, k^2, e, k^2 p, k e)** (k = a_s/3, p = pi_s, e = eps_s of U2; since P' = k P at theta_s, p12 = k p02, p13 = k p03; verified 3e-76), so the norm is a polynomial in the six numbers (k_s, p_s, e_s) per state.
Gamma = Y^-T (G1(x)G1)_c Y^-1, (G1)_c = S^T D S, S = P_phys^{G}(zeta) (exact in zeta). **The norm identity N = nu^T Gamma nu holds: 5e-68 (P), 2e-71 (Pp), 1e-75 (X1), on all 64 states**, and at 5 other thetas (1e-67..1e-69). Completeness: sum_A nu_A nu_A^T / N_A = Gamma^-1 to 5e-74 (1/N_A is the SoV measure on the spectrum, with moments Gamma^-1).
Structure:
* **t = 1/theta12 dependence is exact scaling:** Y(t) = t^L Y(1) t^-L, L = (level of row jet = w_i+w_j-4 with one-site levels (0,1,1,2,2,3,3,4); depth of column). All 833 nonzero entries of Y are pure monomials c t^k, k = 0..4, verified at 3 thetas; it also reproduces the norm at 2 further thetas. Hence Gamma_ab(t) = sum_{c,d} (Y^-1)_{ca}(G1(x)G1)_{cd}(Y^-1)_{db} t^{(L_c-L_a)+(L_d-L_b)} is a polynomial in t of degree <= 10 (lowest 0); Gamma is dense (4096 nonzero).
* **Gamma(t=0) = Gamma^1 (x) Gamma^1** exactly (reshaped rank 1), Gamma^1 = Y0^-T G1 Y0^-1 with Y(0) = Y^1 (x) Y^1, Y^1 an explicit 8x8 (block diag, entries -6..6). For t != 0 Gamma does NOT factorise (reshaped singular values 1, .41, .41, .17, ...); coupling is Gamma = B^-T [C(t)^-T (G1(x)G1) C(t)^-1] B^-1, C(t) = Y(t)Y(0)^-1 graded unipotent (log C has shifts 1..7, so not exp(tK) of a single level-1 generator).
* **zeta dependence:** (G1)_c entries are rational, **S3-symmetric in zeta, denominator Delta^4 (Delta = Vandermonde)**, homogeneous of degree -2-(L_i+L_j); G1 positive definite. So Gamma_ab^{(k)} (coefficient of t^k) is homogeneous of degree -4-L_a-L_b-k in zeta, i.e. Gamma depends on t/zeta only (up to the weight prefactor).
* **5x5 block does not suffice:** fitting N by the 25 products of (1,k,p,kp,e) leaves residual 1.5e7 of 9.9e8 (adding k^2: 7e5); all 8 per site are needed (k^2, k^2 p, k e are genuinely independent monomials).
* Relation to (ii): the ratio Gamma-norm / [(9/4)|J6| prod ...] = 1 state independently (same N), so the SoV measure and the no-log Jacobian are the same functional; (ii) is the compact closed form.

## 3. Untwisting limit and the normalised CG table (D2)
Regime caches (12 eps, 200 digits, phi and phi2). N formula vs direct G-norm: ratio-1 = 1e-165..1e-180 at every branch.
**lim eps^(4+2k_B) N = 1/X_B^2** (k_B = 0,1,1,2,2,4): checked against W4's exact E5_KB (lead^2/gnorm0, `ref/untwisted-exact-210.txt`) to 60 digits (sympy-limited; 1e-150 for the polynomial (4,2,0)), phi and phi2. E.g. phi = (1/3,-2/7,5/11): 1/X^2 = 6.80932040194671, 3305.17309616675, 181.200062269537, 515347.383730086 (both (3,2,1) copies), 20448728043.9319; N ~ eps^-(4+2k_B) confirmed.
Depth: J-route needs det J6 to relative accuracy eps^(8-2k_B) (cancellation of det terms: ~1e54, 1e42, 1e41, 1e26, 0 digits at eps=1e-6, i.e. Laurent data of the branch to ~eps^9, 7, 7, 5, 0 order; singlet none); the Gamma-route cancels ~eps^-14, -12, -12, -14 ... (max term / result 1e84, 1e72, ..., 1e58 at eps = 1e-6) so it needs nu to ~14 orders; the J-formula is the economical one. (The nu-form leading-polynomial identity itself still needs U2's depths 3,5,5,4,2.)
**Fully normalised CG from Q-data only:** Psi_0^Q = X_Q v_q, X_Q = 1/sqrt(lim eps^(4+2k_B) N), v_q = U2's dictionary vector (shape from the nu-series leading polynomial); the closure ||X_Q v_q||_G = 1 is a genuine check (G-norm taken on the chain side only for that check). Result (phi; phi2 identical):

| branch | CG slice (Q-side = chain/exact) | max over all 64 components of |Psi_0^Q - Psi_0| |
|---|---|---|
| (4,2,0) | C_7 = 1 | 3e-148 |
| (4,1,1) | C_5 = -0.70710678118655 | 1e-65 |
| (3,3,0) | C_6 = -0.70710678118655 | 1e-65 |
| (3,2,1)a | C_2 = 0.208887289734197, C_4 = 0.0493116000129644 | 3e-77 |
| (3,2,1)b | C_2 = -0.0797878448606162, C_4 = 0.337986734607777 | 3e-77 |
| (2,2,2) | C_0 = 0.117851130197758 | 4e-102 |

||X_Q v_q||_G = 1.0 to all printed digits on every branch. The overall sign cannot come from the Q side (aligned to the chain).

## Surprises / not done
* The norm is NOT a monomial in Q-values (det N_B, M_s, P^j(theta_s) fail); it is a Jacobian times a pure power of zeta differences, theta-independent.
* Y(t) = t^L Y(1) t^-L: the whole theta12 coupling is a grading.
* Not done: sign of det J6 (open), a closed form of Gamma^1/C(t) (only exact numerics), symbolic factorisation of det J6 (sympy too slow, aborted), general Lambda, the theta-extrapolated Y(t) is verified only through the norm identity. No Yang-Yang Hessian in Bethe roots attempted.
