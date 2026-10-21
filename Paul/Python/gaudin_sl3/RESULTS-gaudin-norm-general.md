# RESULTS-gaudin-norm-general (agent Q2, round 5, 2026-10-16): the Gaudin norm formula for general Lambda

Files: `gaudin_norm_gen.py` (library: `table(lam, tag, th)` per-state table cached in `gaudin_norm_gen_tab_*.pkl`; `jac_funcs`; `Cpred`), drivers `q2_fit.py`, `q2_final.py`, `q2_theta.py`, `q2_t3.py`, `q2_sign*.py`, `q2_reg.py`,
printouts `q2_final_output.txt`, `q2_sign_output.txt`, `q2_reg_output.txt`. Python only, 80 digits (regime check 200 digits), reuses `nu_general`, `general_kernel` read-only.
Data: [3,1,0], [3,2,0] (225 states each) at 13 zeta points (P, Pp, X1, X2, X3 and 8 new ones `q2a..q2h` in `gaudin_norm_gen.ZQ`, theta = (2/7,-3/5)), 13 thetas each at zeta = P; extra reps in 1.3.

## 0. Conventions

* N(Psi) := <Psi|Psi>_G / <x0^G|Psi>^2, G = invariant diagonal one-site metric (`nu_general.gram_diag`, G[hw] = 1), two-site metric = Kronecker square, x0^G = x(zeta) (x) x(zeta), x = `general_kernel.x_cov`. **Metric check:** the Gram diagonal, in W4's GT order, equals W4's `GMET1` record exactly for [2,1,0], [3,1,0], [3,2,0]. N > 0 on all states.
* Condition normalisation (the constant depends on it): the six no-log conditions are N1's raw polynomials `frobenius(lam)['nolog'] = [A, B1, B2]` per site (coefficient of x^(m-3) of the dual ODE at the resonant order m = lam2+1 for the u-solution (c_0 = 1), and at m = lam1+2 for u and for v (c_{lam2+1} = 1)); no rescaling. Order: (A_1,B1_1,B2_1,A_2,B1_2,B2_2) against the variables (H_1,a_1,b_1,H_2,a_2,b_2). Site 2 = site 1 with (H,a,b)_s <-> (H,a,b)_t and delta -> -delta, delta = theta_1 - theta_2. The other site enters through the regular parts, so all cross-derivatives (condition A for lam2 = 2, and B1) are in the full 6x6 J6. **This normalisation reproduces Q1's adjoint J6 exactly (same 9/4, same c0 = 698440.12045... at P).** Rescaling a condition by f rescales C by f^2.

## 1. The general formula

    N(Psi) = C(Lambda) |det J6| (zeta_2-zeta_3)^(mu_1-c) (zeta_1-zeta_3)^(mu_2-c) (zeta_1-zeta_2)^(mu_3-c),    mu = dual-Q degrees (sector)

* **Exponents: p_23 = mu_1 - c, p_13 = mu_2 - c, p_12 = mu_3 - c, with c = 2(lam1+lam2) = 2|Lambda| (lam3 = 0)**: adjoint c = 6, [3,1,0] c = 8, [3,2,0] c = 10.
  How found: (i) at each single zeta, log(N/|det J6|) regressed on (1, mu_1, mu_2, mu_3) over all 225 states: the coefficient differences (mu_1-mu_3, mu_2-mu_3) equal (log d23 - log d12, log d13 - log d12) to 1e-14 (double-precision lstsq) at all 13 zeta points, both reps (exponent of d_jk is mu_i + const, slope exactly 1);
  (ii) c0 = N/(|J| d23^mu1 d13^mu2 d12^mu3) is state independent (spread <= 7e-68 over 225 states, 13 zeta points, 2 reps), and over the 13 zeta points c0 = C (d12 d13 d23)^(-c) with fitted slopes -8.00000000 ([3,1,0]) and -10.00000000 ([3,2,0]).
* **Constant (in the above normalisation): C = [ lam2! (lam1-lam2)! (lam1+1)! / ((lam2+1)(lam1-lam2+1)) ]^2**, i.e. sqrt(C) = (lam2! (lam1-lam2)!)^2 (lam1+1)! / ((lam2+1)! (lam1-lam2+1)!), manifestly invariant under the duality lam2 -> lam1-lam2.
  Values: [2,1,0] 9/4; [3,1,0] 64; [3,2,0] 64. Found from the 15 reps of 1.3 (lam2 = 0: sqrt C = (lam1!)^2; lam2 = 1: ((lam1-1)!)^2 (lam1+1)/2; then the symmetric guess) and verified on all of them, not only the two target reps.
  C is NOT (c_2/2)^2 with c_2 = lam2(lam1+1) (that would give 9/4, 4, 16): the adjoint 9/4 = (3/2)^2 is a coincidence of the factorial formula.

### 1.1 Residuals: max over all 225 states of |N / (C |J6| prod d^(mu-c)) - 1|, 80 digits
| zeta tag | [3,1,0] | [3,2,0] |
|---|---|---|
| P (3/2,5/4,2/3) | 1.6e-69 | 3.3e-70 |
| Pp (2,5/3,-1/2) | 8.6e-71 | 4.9e-71 |
| X1 (7/5,-3/4,1/6) | 3.4e-73 | 1.6e-73 |
| X2 (-4/3,2/9,5/7) | 1.2e-72 | 3.8e-72 |
| X3 (1/3,-2/7,5/11) | 7.3e-68 | 3.1e-68 |
| q2a..q2h (8 more) | <= 5.9e-72 | <= 2.6e-70 |

Worst over the 13 zeta points: 7.4e-68 and 3.1e-68 (limited by the 80-digit data).

### 1.2 theta independence (zeta = P)
13 thetas each: [3,1,0] residuals 9.5e-72 ... 2.4e-68 (theta2 in {-107/91, -113/315, -124/105, -13/28, -23/42, -27/35, -3/5, -31/63, -65/42, -75/56, 11/7, 139/70, 41/21}, theta1 = 2/7); [3,2,0] 1.2e-71 ... 3.5e-69 (also (-1/3,-128/105), (5/4,51/140)). Same C at every theta: N is theta independent.

### 1.3 Other representations (zeta = P, 80 digits; 'spread' = state independence of c0; 'resid' = formula with c = 2(lam1+lam2))
| Lambda | states | C | spread | resid |
|---|---|---|---|---|
| [1,0,0], [1,1,0] | 9 | 1 | 7e-78, 6e-78 | 7e-78, 4e-78 |
| [2,0,0], [2,2,0] | 36 | 16 | 4e-74, 5e-74 | 4e-74, 4e-74 |
| [3,0,0], [3,3,0] | 100 | 1296 | 7e-71, 1e-69 | 7e-71, 1e-69 |
| [4,0,0], [4,4,0] | 225 | 331776 | 1e-66, 6e-65 | 1e-66, 6e-65 |
| [2,1,0] | 64 | 9/4 | 1e-73 | 1e-73 |
| [4,1,0], [4,3,0] | 576 | 8100 | 1.9e-64, 4e-65 | same |
| [4,2,0] | 729 | 25600/9 | 7e-66 | 7e-66 |
| [5,1,0] | 1225 | 2985984 | 2.6e-60 | 2.6e-60 |

Small reps also at X1: residuals <= 1.5e-69, which confirms c independently of the constant. Dual pairs (lam1,lam2) <-> (lam1,lam1-lam2) have the same C but different c (e.g. [4,1,0] c = 10, [4,3,0] c = 14). The formula also holds for lam2 = 0 and lam2 = lam1, so it is not tied to the lam2 = 1, 2 structure of the conditions.

## 2. Sign of det J6 (solved; hypothesis (a) false)
* (a) sign(det J6) = sign(det J3^(1) det J3^(2)) (diagonal blocks): **false**: agreement 108/225, 114/225 ([3,1,0], [3,2,0] at P), i.e. coin flip; det J6/(B1 B2) ranges over [-0.4, 0.56] ... [-758, 61]. J6 is not block triangular.
* **Rule: sign(det J6) = (-1)^(mu_2) * sgn(zeta_2-zeta_3)^(mu_1) sgn(zeta_1-zeta_3)^(mu_2) sgn(zeta_1-zeta_2)^(mu_3).** It depends only on the sector mu, not on the copy or the state within a sector (all 36 sectors of [3,1,0], [3,2,0] and 19 of the adjoint have constant sign at every zeta); the zeta-ordering factor is the sign of the signed power prod (zeta_a-zeta_b)^(mu_c). Tested on [2,1,0], [3,1,0], [3,2,0], [2,2,0], [3,0,0] at 9 zeta points covering all sign patterns of the three differences: 0 violations (`q2_sign_output.txt`). For the larger reps ([4,*,0], [5,1,0]; 3556 states at the ordered P) plain (-1)^(mu_2) has 0 violations. At P for [3,1,0]: (1,2,5)+, (1,3,4)-, (1,4,3)+, (2,3,3)-, (3,2,3)+.
* Hence |.| can be dropped: **N = (-1)^(mu_2) C det J6 (zeta_2-zeta_3)^(mu_1-c) (zeta_1-zeta_3)^(mu_2-c) (zeta_1-zeta_2)^(mu_3-c)** with signed differences (c even). The (-1)^(mu_2) refers to the condition/variable ordering of section 0 (a transposition of two conditions flips it). The "(b) Bethe-ordering" test was unnecessary.

## 3. Regime (untwisting limit) check: done
All 20 branches of [3,1,0], [3,2,0], both phi and phi' (24 eps points, 200 digits, N1's `nu_general_reg_*` caches, theta = (2/7,-3/5)):
* formula (C, c as above, |det J6| from the branch's Laurent data) vs direct G-norm/overlap^2: max relative deviation 4e-100 ... 1e-147 over the 24 eps points.
* **lim eps^p N = 1/X_B^2 with p = 2(2 lam2 + k_B) = 4 lam2 + 2 k_B** ([3,1,0]: 4 + 2k_B; [3,2,0]: 8 + 2k_B; adjoint 4 + 2k_B), from <x0^G|Psi_B> = eps^(2 lam2 + k_B) X_B + ..., ||Psi||_G = 1. The limit (24-point series fit of eps^p N) agrees with N1's X_B: relative difference 1e-99 ... 1e-146 on all 40 branch x phi combinations (`q2_reg_output.txt`). Example [3,1,0], phi: 1/X_B^2 = 6.80932040195, 3305.17309617, 362.400124539, 828581.965404 (both copies of (5,2,1)), 8136.85725012, 53318866.0717 (both (4,3,1)), 88853375578.1, 2.49093385945e12.

## Surprises / not done
* c = 2|Lambda| exactly, and C is a pure factorial product; the adjoint 9/4 is its lam = (2,1) value.
* The sign is a function of mu_2 (times a zeta-ordering sign), so the apparent "no definite sign" came from mixing sectors.
* Not done: lam3 > 0 (presumably a det-twist of the same formula; untested); no proof, the formula is verified numerically on 15 reps with lam3 = 0; no closed form of det J6; the "natural" condition normalisation that makes C = 1 was not identified (the factorials suggest k!/l!-type factors from the Plucker coordinates); Gamma/SoV-measure route not repeated for general Lambda; theta variation tested only at zeta = P.
