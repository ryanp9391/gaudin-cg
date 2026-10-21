# RESULTS-bethe-series (agent S1, round 5, WP-C1 + WP-F3, 2026-10-16): the eps-series of the Laurent data from the Bethe side alone, and `cg_from_Q.py`

Files (all mine): `bethe_series.py` (solver; `python -X utf8 bethe_series.py --stage sym210|sym310|sym320` writes `bethe_series_<lam>_sym.pkl`), `cg_from_Q.py` (F3; ~3 s), `bethe_series_{210,310,320}_sym.pkl`
(exact series polynomial in phi: `{(mu, copy): {'H1': [str per order], 'a1', 'b1', 'H2', 'a2', 'b2'}}`, sympy strings in p1,p2,p3, units theta1-theta2 = 1; H,a at order n scale as theta12^(n-1), b as theta12^(n-2)).
Python only (sympy for the symbolic/eps=0 parts, own exact arithmetic over Q(sqrt D) for everything else). Reuses `nu_general.frobenius/nu_symbolic/load_w4/rep_of/w4_cg`, `untwist_limit.gf_polys`, `general_kernel.load_rep`, `untwisted_Q.frobenius_conditions` read-only.

## 0. The equations and their ranks

Unknowns x = (H_1,a_1,b_1,H_2,a_2,b_2); units theta_1 - theta_2 = 1, zeta = eps*phi (phi a direction; x_n is homogeneous of degree n in phi). Eight equations F(x; eps) = 0:
* 6 no-log conditions: per site A (order lam2+1), B1, B2 (order lam1+2) from `nu_general.frobenius(lam)['nolog']` (S1,S2,S3 regular parts through the other site's data and sigma_1, e_2, e_3 of zeta; for lam2 = 1, A is b_s (lam1+1) = a_s(H_s - a_s); for [3,2,0] A is cubic and couples the sites);
* 2 global rules, **general in (Lambda, mu)**: H_1 + H_2 = sum_{i != j} zeta_i mu_j, b_1 + b_2 = zeta_1 zeta_2 mu_3 + zeta_1 zeta_3 mu_2 + zeta_2 zeta_3 mu_1.
  (Derivation: P^j = e^{zeta_j u} p_j with deg p_j = mu_j; order u^(m_j-1) of the oper gives q(zeta_j) = Pi'(zeta_j) mu_j for q(z) = 2|L| z^2 - (H_1+H_2) z + (b_1+b_2), Pi = prod (z - zeta_k); Lagrange interpolation gives the two rules and 2|L| = |mu|.
  A third rule a_1 + a_2 + theta_1 b_1 + theta_2 b_2 = ... does NOT exist at finite zeta: at the next order of 1/u the subleading coefficients c_j of the p_j enter with the non-vanishing coefficient -Pi'(zeta_j), so that order only determines them. It is not needed.)
  All six no-log conditions and the two rules are independent of nothing else: they are the complete input. The rules were confirmed on all 26 branches (exact agreement with W5's chain-derived series below), for all three reps.

**Order 0** (zeta = 0). The eps = 0 system alone has a zero-dimensional solution set that is **exactly the set of branches**: 6 points for [2,1,0], 10 points for [3,1,0] (`order0_points`: Groebner/lex + exact solve; the eps^0 values of W4/U1 (H_1 = 4,6,6,9,9,12; octet a_1 = (15 -+ 3 sqrt5)/2, b_1 = +-3 sqrt5; [3,1,0]: H_1 = 6,8,10,12,13,13,16,16,18,20 ...) all follow from the Bethe equations; the point is also checked (`check_order0`) for [3,2,0] with W4's values; the lex Groebner basis of [3,2,0] with all six variables did not finish in the time I gave it, see "not done"). H_2 = -H_1, b_2 = -b_1 at eps = 0.
Total multiplicity of the zero-dimensional ideal of [2,1,0] is 10 for 6 distinct points: the top state is a triple root, each octet copy a double root, the other three simple. In the lex basis the H_1-polynomial of [3,1,0] is (X-6)^2 (X-10)^2 (X-13)^2 (X-16)^2 (X-8)(X-12)(X-18)(X-20).

**Rank of the Jacobian J = dF/dx (8 x 6) at (x_0, eps = 0):** rank 6 on the simple roots, **rank 5 exactly on the multiple roots**:

| rep | rank 5 (kernel dim 1, 3 left null vectors) | rank 6 |
|---|---|---|
| [2,1,0] | (4,2,0) top, (3,2,1) x2 | (4,1,1), (3,3,0), (2,2,2) |
| [3,1,0] | (6,2,0), (5,3,0), (5,2,1) x2, (4,3,1) x2 | (6,1,1), (4,4,0), (4,2,2), (3,3,2) |
| [3,2,0] | (6,4,0), (6,3,1), (5,4,1) x2, (5,3,2) x2 | (6,2,2), (5,5,0), (4,4,2), (4,3,3) |

**Which equations are independent (rank of J restricted to the six no-log conditions / all eight, at eps = 0, units theta12 = 1):** each site's three conditions are independent (rank 3 per site, always); the six together have rank
3 (top-type: (4,2,0), (6,2,0), (5,3,0), (6,4,0), (6,3,1)), 4 (all octets, (4,1,1), (3,3,0), (6,1,1), (4,4,0), (4,2,2), (6,2,2), (5,5,0), (4,4,2)), 5 ((3,3,2), (4,3,3)), 6 (singlet (2,2,2)); the two global rules add
2 (top-type: 3+2 = 5; (4,1,1)-type: 4+2 = 6; (3,3,2): 5+1), 1 for the octets (4+1 = 5) and 0 for the singlet. So at eps = 0 the six no-log conditions are NOT independent and the two sum rules are essential; at eps > 0 the 6x6 Jacobian J6 of the no-log conditions
alone is non-singular with det J6 ~ eps^(8-2k_B) (Q1's norm Jacobian, section 3), i.e. the rank deficiency at eps = 0 is exactly the zero of the norm 1/<Psi|Psi>-type determinant.

## 1. The degenerate solver and the adjoint series

**Algorithm (`bethe_series.py`).** Order n: J x_n = -r_n(x_0..x_{n-1}) with r_n = [F(x_{<n}; eps)]_n (exact monomial-series arithmetic in Q(sqrt D), phi a rational direction).
* rank 6: x_n unique; the 2 extra equations are automatically consistent (checked at every order).
* rank 5: row-reduce J: x_n = x_n^p + c_n k (k = kernel vector, one free coordinate c_n per order), and the 3 consistency conditions g_n = N r_n = 0 (N = the 3 left null vectors) are polynomial in the *pending* c's. The c_{n-1} is fixed by the common root of the three polynomials g_n(c_{n-1}) (gcd over Q(sqrt D)); if the order-n conditions are empty (identically zero) the next order is used (second-order degenerate perturbation theory). Roots are taken in K = Q(sqrt D) (x_n are polynomial in phi, so rational at rational phi); if several K-roots exist, each is followed (depth-first) and the ones inconsistent at the next order are discarded. **In all 15 rank-5 branches exactly one candidate survives** (up to 2 K-rational candidates for c_1; the loser dies at the very next solvability), so the Bethe equations alone determine the series; no chain data, no fit, no input but the equations and mu.

**Mechanism (what the rank-5 degeneracy is).** At a multiple root of the eps = 0 system the 8 linearised conditions are dependent: there are 3 linear relations (left null vectors), of the form  (1/3) dB1_(site1) + dA_(site2) = 0, 3 dA_(site1) + dB1_(site2) = 0, dA_1 + (1/3) dB1_1 - dB2_1 + dB2_2 = 0  for the adjoint top
(eq. order A_1,B1_1,B2_1,A_2,B1_2,B2_2,R_1,R_2; for the octets the same pairing of B1/B2 of one site with A of the other plus the b-rule enters), i.e. the cross-site B conditions are first-order consequences of the A conditions of the other site. The kernel direction is
k = (0, -3/4, -1, 0, -3/4, +1) for the adjoint top: (a_s, b_s) moving along the one-parameter family b_s = a_s(H_s - a_s)/3 at fixed H_s with b_1 + b_2 = 0, i.e. exactly the direction in which the twist acts to first order (a_1 = a_2 = 3 eps phi_3, b_1 = 4 eps phi_3, b_2 = -4 eps phi_3, H first-order shift orthogonal to k): J cannot see the leading physical response, whose coefficient c_1 is fixed one or two orders later. Two classes:
* **multiplicity-2 (octet-like) sectors, lag 1:** c_n is fixed by the solvability of order n+1; c_1 comes from a QUADRATIC (adjoint: the three conditions have exactly one common root c_1 = -(12 p_1 + 11 p_2 + 22 p_3)/5 = -b_1^(1) (the sum-rule value); [3,1,0], [3,2,0]: only one of the three conditions is nonzero, two K-roots, one dies at order 3), for n >= 2 linear.
* **single sectors with a multiple root (top-type: (4,2,0), (6,2,0), (5,3,0), (6,4,0), (6,3,1)), lag 2:** the solvability at order n+1 is *identically empty* in c_n (second-order degenerate); c_n is fixed at order n+2: c_1 by a CUBIC (adjoint top: the three conditions have degrees 3,2,3 in c_1 with the single common root c_1 = -4 phi_3 (= -20/11 at phi); [3,1,0] top: cubic = (11c + 30)(26411 c^2 + 4508 c - 15664)/..., the quadratic factor has no K-root; c_1 = -6 phi_3), later c_n linear.
  U1's attempt used a single left null vector (a cubic with no rational root); using all three conditions (their gcd) the root is unique and rational.
Free parameters per order: rank 6: 0; rank 5: 1 per order (c_n), all fixed by the next (resp. next-but-one) order; no free parameter survives.

**Results (adjoint).** All six branches, every Laurent datum of both sites, exact to **n = 10** (`compare_w5`) and to n = 8 as polynomials in phi (`bethe_series_210_sym.pkl`, 45-point triangular-grid interpolation, verified at 2 extra directions, exact):
* vs **U1** `untwisted_Q_laurent_exact.pkl`: exact equality of H_1, a_1, b_1, n <= 5, all six branches at phi = (1/3,-2/7,5/11); vs U1's 400-digit numerics `untwisted_Q_series_210_ext.pkl`: agreement through n = 8 to the working 60 digits;
* vs **W4/W5** (`ref/laurent-series-210.txt`, E4_SERIES, n <= 10, both directions phiA, phiB, all six data of both sites): **0 mismatching entries** of 12 runs; E4_SYMBOLIC (orders 0..3, symbolic phi): 144 series, 0 mismatches.
* closed forms recovered as resummed series (checked to n = 8 exactly): top H_1 = 4 + eps(p_1+2p_2+3p_3), a_1 = a_2 = 3 eps p_3, b_1 = eps p_3 (4 + eps(p_1 + 2 p_2)) (all orders: exact polynomial, series terminates at n = 2);
  (4,1,1): H_1 = 6 + eps(2p_1+5p_2+5p_3)/2 + (sqrt(1+x^2/4) - 1), a_1 = 6 + (3/2) eps(p_2+p_3) + 3(sqrt(1+x^2/4) - 1), x = eps(p_2-p_3);
  (3,3,0): H_1 = 6 + (3/2) eps(p_1+p_2+2p_3) + (sqrt(1+y^2/4) - 1), a_1 = 3 eps p_3, y = eps(p_1-p_2): difference of the n <= 8 series to the sqrt expansions: 0;
  singlet: eq. (71) of the notes to eps^5 exactly (0 difference), and new: H_1^(6) = (degree-6, 62208 denominator; printed in the pkl), **H_1^(7) = 0**, H_1^(8) = Q(phi) x (sextic)/17915904; octets: H_1^(1..3), a_1^(0,1), b_1^(0,1) identical to (72).
* runtime: 0.6 s per (branch, direction) to n = 10 for the adjoint, 1-2.4 s for [3,1,0], [3,2,0].

## 2. The other representations [3,1,0], [3,2,0]

The same code (no change but `lam`) with the no-log conditions of `nu_general.frobenius`: all 20 branches solved from the Bethe side alone, unique survivor in every case. **Exact equality with W5's `laurent-series-310.txt`/`-320.txt` (E4_SERIES, all six data of both sites, both directions, orders 0..10): 0 mismatches in 20 + 20 runs.**
The depth needed by N1/Q1 (n = 9 for the height-1 branches (6,1,1), (5,3,0) of [3,1,0], n = 7 for [3,2,0]) is covered with 3 spare orders, at 0.1-2.4 s per branch/direction.
(See `bethe_series_310_sym.pkl`, `bethe_series_320_sym.pkl` for the symbolic-in-phi series: 10 + 10 branches, n = 7 for all, n = 9 for (6,1,1), (5,3,0) of [3,1,0]; symbolic comparison with W5's E4_SYMBOLIC (orders 0..3, symbolic phi): [3,1,0] 240 series, [3,2,0] 240 series, 0 mismatches.)
In [3,2,0] the condition A involves the other site's data (and b_s only linearly) and does not change the structure: the same two classes (lag 1 for the four multiplicity-2 copies pairs, lag 2 for (6,4,0), (6,3,1)).

## 3. `cg_from_Q.py` (F3): adj (x) adj CG table from Q-side data only

`python -X utf8 cg_from_Q.py` (about 3 s; `... phi2` for the second direction): steps in the docstring. It uses only (i) the series of section 1 (computed inside, order 8; ε⁰ values from W4's file but verified to solve the ε=0 equations by `check_order0` and equal to the Groebner solutions), (ii) `nu_symbolic(2,1)` (the 5 nu_i polynomials of the Laurent data), (iii) the committed kernel `general_kernel_rep_210.pkl` (at t = 1/theta12 = 1), (iv) the exact Gauss parametrisation / chat (`untwist_limit.gf_polys`), (v) Q1's norm N(Psi) = (9/4)|det J6| |z12|^(mu3-6)|z13|^(mu2-6)|z23|^(mu1-6) with J6 built from `untwisted_Q.frobenius_conditions` and evaluated on the exact series (det of a 6x6 matrix of series), so |X_B| = lim^(-1/2). No chain diagonalisation, no fit. (W4's exact table is used only as the comparison.)
Output (phi = (1/3,-2/7,5/11), theta12 = 1; identical at phi' = (-2/5,3/8,7/11); 90-digit arithmetic, agreement is at the working precision):

| branch | k_B | N_B | det J6 ~ eps^(8-2k_B) | \|X_B\| (J6 route = W4 (73)) | slice C_nu (gl3_reps index) | max \|Q - exact\| |
|---|---|---|---|---|---|---|
| (4,2,0) | 0 | 6 | eps^8 | 0.38321995464853 | C_7 = 1 | 4e-91 |
| (4,1,1) | 1 | 5 | eps^6 | 0.019638542188461 | C_5 = -0.70710678118655 | 0 |
| (3,3,0) | 1 | 5 | eps^6 | 0.083873966000801 | C_6 = -0.70710678118655 | 2e-91 |
| (3,2,1) a | 2 | 4 | eps^4 | 0.0017756716595655 | C_2 = 0.2088872897342, C_4 = 0.049311600012964 | 1e-91 |
| (3,2,1) b | 2 | 4 | eps^4 | 0.0017756716595655 | C_2 = -0.079787844860616, C_4 = 0.33798673460778 | 6e-92 |
| (2,2,2) | 4 | 2 | eps^0 | 1.1362980267775e-5 | C_0 = 0.11785113019776 | 2e-92 |

* "N_B = 6 - k_B": all lower total orders of sum G nu nu vanish **exactly** (exact arithmetic), the first non-zero order gives k_B = 0,1,1,2,2,4. The coefficient-by-coefficient dictionary (16 equations, 7 unknowns) is exactly consistent in Q(sqrt5).
* The no-log conditions (U1's normalisation) vanish identically along the Bethe-side series (consistency of N1's and U1's equations); det J6 starts at exactly eps^(8-2k_B).
* The sign of a slice is a phase convention; it is aligned to the exact table (not fixed by the Q side). The magnitude (X_B, hence the absolute normalisation, G-norm 1) comes entirely from det J6.

## Surprises
1. The eps = 0 Bethe system alone has exactly the branch points as its solutions (6 of 6, 10 of 10), so even the ε⁰ data are Bethe-side outputs.
2. Rank 5 <=> multiple root of the eps = 0 system; the two sorts of degeneracy (lag 1 octets, lag 2 top-type) and the use of ALL three left null vectors (U1 used one, hence the "cubic without rational root").
3. The first-order physical response lies entirely in the Jacobian's kernel (a_s = 3 eps phi_3 on the top state).
4. Nothing beyond the 6 no-log conditions + 2 global rules is needed; in particular no extra equation chooses the branch: K-rationality plus consistency at the next order selects a unique solution in all 15 degenerate branches.
5. The whole CG table (including absolute normalisation) is a ~3 s computation from the Bethe series.

## Not done
* Order-0 solution set of [3,2,0] by Groebner basis (6-variable lex did not finish); W4's eps^0 values were verified to satisfy the equations (ranks above) but not proved to be all solutions.
* No proof that the selection (unique K-rational survivor) holds for all Lambda; no closed algebraic form of the octet/singlet series; [3,1,0], [3,2,0] not packaged in `cg_from_Q.py` (needs a Q-side norm for general Lambda: Q1's J6 formula is verified for the adjoint only) — their slices from the same ingredients were done numerically by N1.
* No independent theoretical derivation of the (9/4)|det J6| norm (Q1's regression).
