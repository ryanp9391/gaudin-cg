# RESULTS-Qside-320: the Q side of the sl(3) Gaudin chain for Lambda = (3,2,0) (L = 2) and the [3,2,0] jet formula (agent P4, 2026-10-16)

Files: `qside_320.py` (copy of P3's `qside_310.py`, LAM = (3,2,0), everything parametrised; `python -X utf8 qside_320.py [0 1 2 3 4]`, about 4 min with caches, ~8 min cold; printed output of the last run in `qside_320_output.txt`; caches `qside_320_{states,pscan,recs,q4_*}.pkl`), `twovac_fit_320.py` (J; `python -X utf8 twovac_fit_320.py [--stages T0,...,T7]`, caches `twovac_fit_320_*.pkl`, printed output in `twovac_fit_320_output.txt`, stage logs `twovac_fit_320_T*.txt`), helper builders `twovac_fit_320_build*.py` (parallel point builders, no results of their own). Reuses `qside.py`, `gaudin_chain.py`, `qside_310.py`, `twovac_fit.py` read-only. No Wolfram. 60 digits (Q), 80 digits (J), 100 digits (Q4). Python only.

Headline: every Q-side rule of (3,1,0) holds for its dual with the expected duality swaps, with **three corrections to the plan's table** (marked CORRECTED): vacuum a_s, dual Wronskian F', direct local exponents. The jet formula **converges on the first hypothesis level**: with c = prod_w (z2-w)(z3-w)^3 (z1-z2)^4 and M = P3 p03 the target is exactly (residual 0.0 at 80 digits) a bilinear form in **P3 times the seven non-forced Pluecker coordinates {p03,p04,p13,p14,p23,p24,p34}** per site; unique, rational, polynomial kernel; **7 monomials per site** (adjoint 5).

## 0. Conventions, Laurent structure, vacuum

* Rep dim 15, 225 two-site states, 36 weight sectors (largest 18). **HW ket index (0-based, gl3_reps order) = 14**, GT pattern ((3),(3,2),(3,2,0)) (last index, as for (2,1,0), (3,1,0)); annihilated by E12, E13, E23. Nilpotency index 4 for E12, E13, E23 (as (3,1,0)). Top two-site state e1 x e1 = index 224, weight (6,4,0).
* cdet exact checks: [S2(3/7),S3(-2/9)] = 0, all [S_k(u),S_l(v)] = 0 on a 3x3 grid, weight conservation; S1 = e1(zeta) + 5/(u-th1) + 5/(u-th2) (scalar); S2 constant = e2(zeta) = 89/24, S3 constant = e3(zeta) = 5/4.
* **c2 = 8 = lam2(lam1+1) exactly** (scalar double pole of S2 at both sites; Casimir/2 would be 7). **S3 has no triple pole** (operators identically zero). Six non-scalar data (H1,H2,a1,a2,b1,b2). Fits exact (0 mismatches). H_s has the same explicit formula as for (2,1,0).
* Indicial polynomial e(e-1)(e-2) - 5e(e-1) + 8e = e(e-3)(e-5): **local exponents {0,3,5} = {0, lam2+1, lam1+2}** (prediction confirmed).
* Vacuum (r = (4/5,-7/11,65/44), variant 0 of g'): exact eigenvalue functions
  * omega: S1 = 163/30 + 5/x1 + 5/x2; S2 = 131/15 + (885/31)/x1 + 8/x1^2 + (45/31)/x2 + 8/x2^2; S3 = 7/2 + (2376/155)/x1 + (24/5)/x1^2 + (-144/155)/x2 + (24/5)/x2^2
  * omega': S1 = 369/140 + 5/x1 + 5/x2; S2 = 29/28 + (48341/2170)/x1 + 8/x1^2 + (-10459/2170)/x2 + 8/x2^2; S3 = -9/10 + (44217/2170)/x1 + (72/7)/x1^2 + (-31383/2170)/x2 + (72/7)/x2^2
  * **H_s^vac = sum_i(|L|-lam_i) omega_i +- 12/theta12: confirmed** (H_s - 2w1 - 3w2 - 5w3 = +420/31, -420/31 = +-12/theta12, theta12 = 31/35; both omega and omega').
  * **b_s^vac = omega3 (H_s - 5 omega3): confirmed exactly (difference 0), both omega, omega'.**
  * **a_s^vac = 8 omega3 = c2 omega3 (CORRECTED; the plan's 5 omega3 = |L| omega3 is refuted).** (a = 24/5 at omega3 = 3/5, 72/7 at 9/7.) For (2,1,0) and (3,1,0) c2 = |L| so the earlier closed form could not distinguish them. Consistency: kappa_vac = a/c2 = omega3 = P^3 exponent.
  * Vacuum Q degrees: dual (6,4,0) (= weight (6,4,0)), direct (0,2,6) = 6 - dual; kernel dimension 1.

## Q1. Opers, degrees, gauge, exponents (225 states, 60 digits, eigen-residual 9.4e-57)

* Dual oper P^j = e^{zeta_j u} p_j, **deg p_j = mu_j** in all 36 sectors (kernel dim 1, six exponentials, all states); direct Q_i = e^{-zeta_i u} q_i, **deg q_i = 6 - mu_i**; max ODE residual 7.9e-54, largest smallest-singular-value at the degree found 1.05e-55. Per-sector degree table in `qside_320_output.txt` (all 36 rows OK, e.g. (6,4,0): dual (6,4,0), direct (0,2,6); (3,3,4): dual (3,3,4), direct (3,3,2)).
* **Gauge power scan** (one state per sector, p = 0..5, smallest-degree polynomial solution): p = 0,1,2: no solution up to degree 10 in any sector; **p = 3: degrees exactly 6 - mu_i, kernel dim 1, 36/36 sectors**; p = 4: degrees +2 in all 36; p = 5: degrees +4 in all 36. **p = lam1 - lam3 = 3 confirmed. As (3,1,0) has the same lam1 - lam3 (and lam1), this does NOT separate "lam1 - lam3" from lam1 (or lam1 - lam2 + ... forms): lam3 = 0 in all runs so far.**
* **Local exponents at theta_s** (vanishing orders of the span; 450 of 450 (state,site) cases): dual **(0,3,5)** as predicted; direct Q (CORRECTED: not (0,3,5) as for (3,1,0)) **(0,2,5)**: duality swaps lam2+1 <-> lam1-lam2+1 for the direct side ((3,1,0): dual (0,2,5), direct (0,3,5)).
* **Consequences, derived from the indicial recursion** (P = c0 + c1 x + c2 x^2 + ..., x = u - theta_s; orders x^-2 and x^-1 of (d^3 - S1 d^2 + S2 d - S3)P with S1 = 5/x + .., S2 = 8/x^2 + H/x + .., S3 = a/x^2 + b/x + ..; only the leading singular parts and a, H, b enter):
  * **kappa_s := P'(theta_s)/P(theta_s) = a_s / c2 = a_s/8**
  * **kappa2_s := P''(theta_s)/P(theta_s) = (a_s^2 - H_s a_s + 8 b_s)/24** (= 2 c2/c0 with c2 = (a^2 - H a + 8b) c0/48)
  * (adjoint: kappa = a_s/3 = a_s/c2; in general kappa = a/c2 and kappa2 involves H, a, b: the regular parts of S1,S2,S3 first enter at x^3.)
  Verified for all 3 P^j on all 225 states x 2 sites: max relative deviation 3.1e-57 (kappa), 3.6e-57 (kappa2); in the 80-digit data of J: 5.8e-75, 7.4e-75.
* **Jet matrices** [P^j(k)(theta_s)/k!]: 3x3 (k <= 2): rank 1 (450/450); **3x5 (k <= 4): rank 2 (450/450)**, min sigma_2/sigma_1 = 3.0e-6; **3x6 (k <= 5): rank 3 (450/450)**, min sigma_3/sigma_1 = 3.1e-11. Hence all Pluecker coordinates p_01 = p_02 = p_12 of (P^2,P^3) vanish (max relative 1.2e-72) and every p_kl, l <= 4 is proportional to Delta = c_2 d_3 - c_3 d_2 times a universal function of the exponent-0 and exponent-3 Frobenius solutions.
* **No-log (Frobenius) conditions**: exponent-0 series needs one condition at order x^3 (resonance with exponent 3) and two at order x^5 (coefficients of the free c0 and c3; resonance with exponent 5): **3 independent conditions per site, 6 in total on the 6 data** (as for (3,1,0)), polynomials printed in `qside_320_output.txt` (Section 1; `frobenius_polys()` generalised to c2 = 8, |L| = 5, exponents {0,3,5}), order x^3: -8 H^2 a + 12 H a^2 + 8 H a s0 + 64 H b - 4 a^3 - 8 a^2 s0 - 56 a b + 24 a h0 - 64 b s0 - 192 t0 = 0 (x192; here a, b, H, s_k, h_k, t_k as in `frobenius_polys`). Verified on all 225 states x 2 sites: max relative violation 3.0e-59 (x^3), 2.5e-60 (x^5, c0), 2.8e-59 (x^5, c3). (Number of solutions of the system vs 225 states not counted.)
* Negative controls (36 states, degree <= 9): dual ODE with e^{-zeta u}: none; plain adjoint (p = 0) with e^{-zeta u}: none; Pi^-3 gauge with e^{+zeta u}: none.

## Q2. Wronskians and quantum Weyl

* Direct: Q_a Q_b' - Q_b Q_a' = c_ab e^{-(z_a+z_b)u} F(u) p_c: **F = Pi^1 = Pi^(lam1-lam2) confirmed** (deg 2; relative remainder 7.6e-51, |F - Pi| 8.1e-51, 225 states x 3 pairs); **c_ab = zeta_a - zeta_b exactly (deviation 0.0)**.
* Dual: P^a P^b' - P^b P^a' = c'_ab e^{+(z_a+z_b)u} F'(u) q_c: **F' = Pi^2 (CORRECTED: the plan's table says Pi)**: deg F' = 4 in all cases, remainder 5.7e-49, |F' - Pi^2| 4.8e-49, c'_ab = zeta_b - zeta_a exactly. So **F' = Pi^(lam2) = Pi^(|L| - p)** ((3,1,0): Pi; adjoint: Pi; [3,2,0]: Pi^2: the duality swaps F <-> F'); consistent with P3's F' = Pi^(|L|-p).
* **W = det[P,P',P''] = c_W e^{(z1+z2+z3)u} Pi^5** (Pi^|L| as predicted), c_W = prod_{i<j}(zeta_j - zeta_i) = -35/288 (state independent, deviation 1.7e-45).
* Quantum Weyl: S1 = det[P,P',P''']/W, S2, S3 from the 225 states x 3 u: max relative deviations 1.8e-45, 1.7e-45, 1.7e-45. Direct side psi_i = Pi^-3 Q_i: kernel of the monic operator (-S1, 2S1'+S2, -(S1''+S2'+S3)): deviations 1.9e-43, 3.1e-41, 3.3e-42 (numerical derivatives).

## Q3. Lagrange columns, rank tallies, det N, ratio identity (gauge p = 3)

* Same six residue columns as (2,1,0), (3,1,0): N[(s,i,j),H_t] = -Res[psi P'/(u-th_t)], N[.,a_t] = Res[psi P/(u-th_t)^2], N[.,b_t] = Res[psi P/(u-th_t)]; Pi^-3 weights; sum_m Delta_m N_m = 0 (c2 scalar, S3 without triple pole).
* Identity residual max |sum Delta N|/sum |Delta N| = **8.4e-57** over 600 random ordered pairs x 18 rows.
* **Rank tally** (stacked 18x6, rank = #{sigma_k/sigma_1 > 1e-40}), 600 pairs A != B: same sector rank 5 (31/31); different sectors: rank 5 (439), **rank 4 (130)**. **P3's Delta-mu rule holds: the 130 rank-4 pairs are exactly those with (1 positive, 2 negative) components of mu_B - mu_A** (tally by (#pos, #neg, rank): (0,0,5): 31; (1,1,5): 318; (1,2,4): 130; (2,1,5): 121). Gaps: max sigma_6/sigma_1 = 1.1e-57, min sigma_5/sigma_1 among rank-5 pairs 1.4e-11, max sigma_5/sigma_1 among rank-4 pairs 2.2e-58. **A = B: rank 6 for all 225 states** (min sigma_6/sigma_1 = 1.65e-12).
* det N (6x6, rows (s,(i,j)) in {(1,3),(1,2),(2,3)}): |det N|/prod(row norms): A != B (150 pairs) max **2.2e-62**; A = B (225 states) min 2.1e-37, max 2.9e-6 (all non-zero).
* **x0-free ratio identity** R_B(w)/R_B(w') = [D_B(w)/D_B(w')][D_B0(w')/D_B0(w)], B0 = e1 x e1 (index 224, weight (6,4,0)): **max relative deviation 1.6e-51 over ALL 225 states** (all non-degenerate; worst states (4,2,4), (4,1,5), (5,3,2)); e.g. B = 0: 1.08259118072039 on both sides.
* 84-triple scan (same triple both sites): **72 pass (< 1e-30), 6 fail (deviation 1.1 - 2.0), 6 degenerate**; the failing are exactly the 6 permutation triples; the degenerate ones the 6 triples with all i or all j equal. The XXX triple passes.
* Controls (must fail, do): a2 column replaced by Res[psi P/(u-th2)^3]: deviation 1.89; orientation swapped: 1.64; **gauge weight Pi^-2: 1.56; Pi^-4: 1.98**.

## Q4. Regime (zeta = eps phi, omega = eps om^, omega' = eps om^'; eps = 1e-6, 2e-6; 100 digits; eigen-residual 1.8e-97)

* 10 HW branches (sectors = dual images of (3,1,0)^2, mu* = (6-mu3,6-mu2,6-mu1)) picked by the smallest raising residual; at eps = 1e-6 residual vs the next state: (6,4,0) 0 | (5,5,0) 1.9e-7 vs 1.41 | (6,3,1) 2.3e-7 vs 2.8 | (5,4,1) x2 2.4e-7 vs 2.1 | (6,2,2) 5.4e-7 vs 2.0 | (5,3,2) x2 4.1e-7 vs 2.1 | (4,4,2) 3.6e-7 vs 1.7 | (4,3,3) 7.2e-7 vs 1.5: clean gaps (multiplicity 2 in (5,4,1) and (5,3,2) as predicted); residuals double at 2e-6.
* Degrees at the regime: dual = mu, direct = 6 - mu, nullity 1 (kernel threshold 1e-50 needed at 100 digits). Vacuum eigencovectors ok, direct Q degrees (0,2,6).
* Ratio identity on the 10 branches: max deviation 2.2e-52 (eps = 1e-6), 3.8e-54 (eps = 2e-6).
* eps-exponents of |det N(w)|: **16 for (6,4,0), (5,5,0), (6,3,1), (6,2,2); 15 for both (5,4,1) and both (5,3,2); 14 for (4,4,2), (4,3,3)** — exactly the exponents of (3,1,0) transported by duality (16/15/14 pattern); |R(w)| ~ eps^{0,1,1,1,1,2,2,2,2,3}. For both multiplicity-2 sectors the two branches have equal det N (e.g. (5,4,1): 2.01e-98 both).

## J. The [3,2,0] jet formula (80 digits; points zeta in {P,P',X1,X2,X3,(X4,X5)}, each 225 states, built in ~90 s per point in parallel)

### Objects, normalisation (T0)
* Phi_inf = Omega^inf[zeta,w1] x Omega^inf[zeta,w1'] with the rep-independent r_k^inf (`twovac_fit.onesite_inf`); x0^G = (e1^T P_phys^{-1})^{x2} from `twovac_fit.x0_onesite`, which is rep-generic (Gauss-factor construction with the rep's own E's). <x0^G|Psi_B> and <Phi_inf|Psi_B> are non-zero for all 225 states (P: 9e-6 <= |<x0|Psi>| <= 104, 0.024 <= |<Phi|Psi>| <= 273).
* **Frame test of x0** (the analogue of D1's F0): <x0^G|Psi_B>/n_B/[(Q1P3)(th1)(Q1P3)(th2)], n_B = R_B/D_B(omega0): state independent to **1.5e-72** (without the M1 M2 factor: spread 1.02). So x0^G is the right covector up to a state-independent, zeta-dependent scalar, as for the adjoint.
* **Frame-tie scalar** (derived): kappa^G = Delta1^(lam1-lam2) Delta2^(lam2-lam3) Delta3^lam3 = Delta1 Delta2^2 with D1 = w^2/((z2-w)(z3-w)), D2 = w/((z1-z2)(z3-w)), so kappa^G = w^4/((z2-w)(z3-w)^3(z1-z2)^2) and
  **c = w^4 w'^4/(kappa^G kappa^G') = (z2-w)(z3-w)^3 (z2-w')(z3-w')^3 (z1-z2)^4** (general: (z2-w)^(lam1-lam2) (z3-w)^(lam1-lam3) per w, and (z1-z2)^(2(lam2-lam3)) in total; adjoint: (z2-w)(z3-w)^2 per w, (z1-z2)^2; [3,1,0]: (z2-w)^2(z3-w)^3, (z1-z2)^2). Verified: with this c the **per-zeta-point fits (225 rows, 49 unknowns, rank 49) give the same kernel as the 5-point stacked fit to 6.8e-72, 1.5e-74, 1.8e-73, 2.7e-73, 7.6e-71 (P, P', X1, X2, X3)**: no free scalar per zeta. Wrong exponents fail: (1,2,2): resid 0.18; (1,3,1): 0.51; (0,3,2): 0.69; (2,3,2): 0.28; (1,3,3): 0.39; (1,4,2): 0.23.
* Target: Y = c Qnum M1 M2, Qnum = <Phi_inf|Psi>/<x0^G|Psi> at (w1,w1') = (5/2,7/4), theta = (2/7,-3/5); jets of P^j at theta_s to order 4; Pluecker p_kl of (P^2,P^3), 0 <= k < l <= 4 (10); p01 = p02 = p12 zeroed explicitly (they are 1e-72 noise: forced by P' = kappa P, P'' = kappa2 P); the lone factor P3^(c), c <= 2, kept as separate monomials.

### Hypothesis ladder (5 zeta points, 1125 rows; residual = relative ||Ak-y||/||y||, 0.0 = below the 1e-90 floor)
| level | basis per site (same on both sites) | unknowns | rank | null | resid |
|---|---|---|---|---|---|
| a | B1 = {P3^(c) p_kl : c <= 2, (kl) in {03,04,13,14,23,24,34}} (21 -> per-site rank 15) | 225 | 225 | 0 | **0.0** |
| a | B1 with lone c <= 1 (14 -> rank 12) | 144 | 144 | 0 | 0.0 |
| a | **lone c = 0 only: {P3 p_kl}, 7 pairs** | **49** | **49** | **0** | **0.0** |
| b | p_kl with l <= 3 only {03,13,23} x c <= 2 (9 -> rank 6) | 36 | 36 | 0 | 0.50 |
| b | B1 without p34 (18 -> rank 12) | 144 | 144 | 0 | 0.0084 |
| M | M = P3 p03 | 225 | 225 | 0 | 0.0 |
| M | M = P3' p03 (= kappa P3 p03) | 225 | 225 | 0 | 0.0 |
| M | M = P3 p04 | 225 | 225 | 0 | 0.144 |
| M | no M (target c Qnum) | 225 | 225 | 0 | 0.91 |
| ctl | lone index P2 / P1 instead of P3 | 225 | 225 | 0 | 0.988 / 0.99 |
| ctl | both factors from the site-1 jets / site-2 jets | 225 | 71 | 154 | 0.30 / 0.44 |

* **The first level that fits is the first one tried**, and the minimal basis is the **pure Pluecker basis P3 p_kl with all seven non-forced pairs (k<l<=4)**: dropping any one element of the 7 (same on both sites) leaves a non-zero residual (without p03: 0.41; p04: 0.41; p13: 0.16; p14: 0.25; p23: 0.078; p24: 0.14; p34: 0.061; rank 36/36). The lone P3^(c) factors are redundant: P3' p_0l = P3 p_1l and P3'' p_0l = P3 p_2l (l = 3,4) exactly, because p_kl = (c2 d3 - c3 d2) 6 f0_k (l=3) etc. In the nu-variables nu = (kappa, kappa2, pi4 := p04/p03, pi34 := p34/p03): the seven monomials are **(1, kappa, kappa2) x (1, pi4) and pi34**, i.e. **N7 = {P3^(c) p03, P3^(c) p04 (c = 0,1,2), P3 p34}**. nu-monomials of weight <= 3 (13) fail (resid 0.0046), weight <= 4 fit (the support is in weights {0,1,1,2,2,3,4}).
* **Sparse monomial count: 7 per site (the adjoint has 5)**: unique kernel (49 unknowns, rank 49 on 1125 rows; no null space). Per-site rank of the full space of candidates is 15 (21 features), but the kernel lives on 7 of the 15 (the other 8 rows/columns of the 15x15 kernel in the nu-monomial basis vanish identically).
* **Pattern across reps** (this run + the [3,1,0] comparison): the formula is P3 x (all non-vanishing Pluecker coordinates p_kl, 0 <= k < l <= lam1+1) with the vanishing ones those inside the lam2+1 forced jets {0..lam2}: **count = C(lam1+2,2) - C(lam2+1,2)**: adjoint (2,1,0): 6 - 1 = 5 (N5 = P3 {p02,p12,p03,p13,p23}); [3,1,0]: 10 - 1 = 9 (below); [3,2,0]: 10 - 3 = 7 (above).

### Kernel at the base point (w1,w1') = (5/2,7/4), theta12 = 31/35 (all 49 entries rational; basis P3{p03,p04,p13,p14,p23,p24,p34} = P3{p03,p04} x {1,kappa,kappa2} order (P3 p03, P3 p04, P3' p03, P3' p04, P3'' p03, P3'' p04, P3 p34); rows i, columns j; the same matrix in either basis)
```
 -2477531875/3936256       -15649375/95232      -102534604765838875/42259375274496   201835815824125/113600471168   -1066229095827875/227200942336   841817515875/1832265664   1241317375/2770563
  8703625/190464            42875/4608           213468283144375/1022404240512        -1720011432125/10993593984      29908228000/85887453             -4528239625/88658016      -66781625/2144952
  92531938123379125/21129687637248  -1336404322859375/1022404240512  2973626853294125/1022404240512  -7269823449875/5496796992  12484966188275/5496796992  -1902394025/5541126  -179699825/1072476
 -319234542990475/908803769344      2843820794675/21987187968  -651463216075/2748398496  7355513725/59105344  -5263492675/29552672  16750475/476656  147925/11532
 -347834314964675/454401884672      1716815095975/5496796992  -5119851097655/10993593984  13720239855/59105344  -17517770585/59105344  14281435/238328  227285/11532
  1123221350545/14658125312         -14416082485/354632064  4168591945/88658016  -27137425/953312  6967205/238328  -57785/7688  -365/186
 -3124391459/709264128              60638627/17159616  -2906729/1072476  100009/46128  -37121/23064  203/372  1/9
```
Last entry K_77 = 1/9 (also 1/9 for [3,1,0] at its last entry). Kernel is **not symmetric** at fixed theta12 even for w1 = w1' (max |K - K^T|/max|K| = 1.4 at w1 = w1' = 5/2); the **site-exchange symmetry holds in the form K_ij(w,w,thi) = K_ji(w,w,-thi)** (thi = 1/theta12): deviation 5.8e-63 (degree-7 polynomial in thi evaluated at -thi; the relation K + K^T is not the symmetry: 2.0).

### Kernel laws (T4–T6; zeta independent)
* **(w1,w1') grid**: 9 values (5/2, -2/5, 9/7, 3, -5/3, 7/4, 11/3, -7/2, 13/5), 72 ordered pairs, 5 zeta points each: every fit rank 49, resid 0.0 (72/72). Each of the 49 entries is a **polynomial in (w1,w1') with all rational coefficients**, minimal total degree table d_ij (row, col in the basis order above): row 1: 8 7 7 6 6 5 4; rows 2,3: 7 6 6 5 5 4 3; rows 4,5: 6 5 5 4 4 3 2; row 6: 5 4 4 3 3 2 1; row 7: 4 3 3 2 2 1 0, i.e. **d_ij = 8 - (e_i + e_j) with e = (0,1,1,2,2,3,4)** (e = derivative weights of the basis monomials relative to P3 p03: p04, P3', ... ); max degree in each of w1, w1' separately is 4. No poles.
* **theta dependence**: 14 theta variants (single zeta point P; 12 distinct theta12 values), all rank 49 resid 0.0. **Only theta12 matters** (two pairs of thetas with the same theta12 but different theta: kernels agree to 4.8e-72 and 8e-70). Each entry is a **polynomial in thi = 1/theta12 of degree <= 7** (8 coefficients on 12 points: resid 0.0 at degree 7, 4.8e-6 at degree 6, so degree exactly 7 somewhere) (adjoint: degree <= 4). Poles only at theta12 = 0.

### [3,1,0] comparison (same pipeline; points P,P',X1,X2,X3 built for qside_310 at 80 digits; T7)
* [3,1,0] (exponents {0,2,5}: only P' = kappa P forced, p01 = 0, first nonzero p02; lone c <= 1; c = (z2-w)^2(z3-w)^3 (z1-z2)^2 per w, M = P3 p02): full basis {P3^(c) p_kl, c <= 1, 9 pairs (k<l<=4, (01) removed)} (per-site rank 15): resid 0.0 (225 unknowns, rank 225). The adjoint-like N5 = {P3^(c)p02, P3^(c)p03, P3 p23} fails (0.316) and N5 + {P3 p04, P3' p04} fails (0.051); wrong c exponents (1,3,2): 0.32.
* **Greedy backward elimination** on the 18 features ends with the **pure Pluecker basis {P3 p_kl: 9 pairs}, lone c = 0 only**: unknowns 81, rank 81, resid 0.0, unique rational kernel (0 non-rational entries of 81; K_99 = 1/9); the l <= 3 subset (5 coordinates, the adjoint pattern) fails (0.316); dropping any one of the 9 fails (0.013 - 0.27).
* **Monomial count: [3,1,0] has 9 per site, [3,2,0] has 7, the adjoint 5**; all equal C(lam1+2,2) - C(lam2+1,2) (the number of Pluecker coordinates p_kl, k<l<=lam1+1, not forced to zero by the lam2+1 leading jets being proportional). [3,1,0] and [3,2,0] share lam1 and hence the jet order 4 and the ladder length 5, but differ in the forced zeros.

## Surprises / corrections to the plan's table

1. **a_s^vac = c2 omega3 = 8 omega3**, not |L| omega3 = 5 omega3 (b_s^vac = omega3 (H_s - 5 omega3) IS |L|-type, so the closed form has both c2 and |L| in it).
2. **Dual Wronskian F' = Pi^2** (not Pi): F' = Pi^lam2, F = Pi^(lam1-lam2); with W = Pi^|L| and F' = Pi^(|L|-p).
3. **Direct local exponents are (0,2,5)**, not (0,3,5): the duality swaps (0,2,5) <-> (0,3,5) between the dual and direct opers when lam2 <-> lam1-lam2.
4. kappa2_s = (a^2 - H a + 8 b)/24 involves H and b; kappa_s = a/c2 alone.
5. The jet formula needs no lone P3^(c) at all in the Pluecker form: P3 p_kl with the 7 (resp. 9) non-forced coordinates; the "lone c <= 2" monomials of the plan are P3 p_1l, P3 p_2l.
6. The (w1,w1') law has total degree <= 8 - e_i - e_j (adjoint: 10 - w_i - w_j with weights (2,3,3,4,5)); the theta12 law has degree 7 (adjoint 4).
7. p = 3 does not separate p = lam1 - lam3 from p = lam1 (lam3 = 0 throughout).

## Not done

* Frobenius-system solution count vs 225 states; a closed form of the kernel entries beyond the polynomial reconstruction (degree table only; coefficients are rational and stored in `twovac_fit_320_grid.pkl`, not printed); scaling u -> cu check of the homogeneity degree (not run; the table d_ij = 8 - e_i - e_j is the observed w-degree); the regime (untwisting) check of the Q-side formula against the chain (D1's F4 analogue) and the CG extraction for [3,2,0] (W3's job on the chain side); [3,1,0] kernel laws (grid/theta).
