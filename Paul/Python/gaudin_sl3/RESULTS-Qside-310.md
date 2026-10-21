# RESULTS-Qside-310: the Q side of the sl(3) Gaudin chain for Lambda = (3,1,0) (L = 2), and the rank-4 puzzle of (2,1,0) — agent P3, 2026-10-16

Files: `qside_310.py` (everything; `python -X utf8 qside_310.py`, ~4 min with caches, ~15 min cold; optional section arguments `0 1 2 3 4 R`), `qside_310_output.txt` (full printed output of the last run; all numbers below are copied from it or from the logged runs), caches `qside_310_{states,pscan,recs,R,q4_*}.pkl`. Imports P1's `qside.py`/`gaudin_chain.py`/`gl3_reps.py` read-only (hard-wired parts re-written in `qside_310.py`). No Wolfram, no notebooks. 60 digits (eigen-residual 2.3e-56), 100 digits in Q4.

## 0. Conventions for (3,1,0)

* Rep dim 15, 225 two-site states, 36 weight sectors (largest 18). **HW ket index (0-based) = 14** (GT pattern ((3),(3,1),(3,1,0)), the LAST index; for (2,1,0) it was 7 = dim-1 as well); annihilated by E12, E13, E23 (checked). Top two-site state e1 x e1 = index 224, sector (6,2,0).
* cdet: same symbolic S1,S2,S3 as round 1. Exact checks: [S2(3/7),S3(-2/9)] = 0; all [S_k(u),S_l(v)] = 0 on a 3x3 rational grid; weight conservation; **S1 = e1(zeta) + 4/(u-th1) + 4/(u-th2)** (scalar); S2 constant e2(zeta) = 89/24, S3 constant e3(zeta) = 5/4 (scalars; S_k(inf) = e_k).
* **Laurent structure**: S2 = e2 + sum_s[H_s/x_s + c2/x_s^2], S3 = e3 + sum_s[b_s/x_s + a_s/x_s^2]; **S3 has no triple pole** (operators zero); fits exact (0 mismatches at extra points); six non-scalar data (H1,H2,a1,a2,b1,b2) (all six operators non-scalar). H_s has the same explicit form as for (2,1,0) (checked exactly).
* **c2(3,1,0) = 4, NOT the Casimir/2 (= 6)**. Scan over reps (exact fits, two sites): c2 = lam2 (lam1 + 1) for lam3 = 0, lam2 > 0: (1,1,0): 2, (2,1,0): 3, (3,1,0): 4, (2,2,0): 6, (3,2,0): 8; for lam2 = 0 (the symmetric reps (1,0,0), (2,0,0), (3,0,0)) the double pole of S2 is NOT scalar; (2,1,1): c2 = 8 but there S3 has a triple pole. So c2 is not the quadratic Casimir (Casimir/2 would give 3, 3, 6 for (2,1,0),(3,1,0),...: only coincides for (2,1,0)).
* **Indicial polynomial** of the dual oper at theta_s: e(e-1)(e-2) - |Lambda| e(e-1) + c2 e = e(e-2)(e-5): **local exponents {0, 2, 5}** = {0, lam2+1, lam1+2} (round 1, (2,1,0): {0,2,4}). Verified on all 225 states x 2 sites from the vanishing orders of the span of the three P^j: always (0,2,5).
* Vacuum: Omega_{zeta omega} = e1^T rho(n), r = (4/5, -7/11, 65/44) (rep-independent formulas), exp series terminates; **nilpotency index (smallest k with E^k = 0 on the 15-dim rep) is 4 for E12, E23 and E13 (not 3 as the spec guessed)**. g' = P gc P^-1 (variant 0) works: Omega x Omega is an exact left eigencovector of all S2,S3 coefficient operators for omega and omega' (variants P^-1 gc P and (P gc P^-1)^T fail). Exact vacuum eigenvalue functions:
  * omega: S1 = 163/30 + 4/x1 + 4/x2; S2 = 131/15 + (5789/310)/x1 + 4/x1^2 + (1589/310)/x2 + 4/x2^2; S3 = 7/2 + (3027/310)/x1 + (12/5)/x1^2 + (507/310)/x2 + (12/5)/x2^2 (no triple pole)
  * omega': S1 = 369/140 + 4/x1 + 4/x2; S2 = 29/28 + (54107/4340)/x1 + 4/x1^2 + (-4693/4340)/x2 + 4/x2^2; S3 = -9/10 + (40869/4340)/x1 + (36/7)/x1^2 + (-34731/4340)/x2 + (36/7)/x2^2
  * vacuum direct Q degrees (0,4,6), vacuum dual degrees (6,2,0), kernel dimension 1, both twists (so the dual vacuum degrees are the weight (6,2,0), direct = 6 - weight).

## Q1. Opers, degrees, the gauge power

* Eigenstates: per weight sector of S2(3/7) + (pi/3) S3(-2/9) (the problem is real; `mp.eig` returns 1e-59 imaginary noise in some sectors, removed). 225 states, residual 2.3e-56.
* Dual oper (d^3 - S1 d^2 + S2 d - S3) P = 0: P^j = e^{+zeta_j u} p_j with **deg p_j = mu_j** in every one of the 36 sectors (monic, kernel dim exactly 1 for each of the 6 exponentials, all 225 states). Direct: Q_i = e^{-zeta_i u} q_i with **deg q_i = 6 - mu_i = 2 lam1 - mu_i**; max ODE residual (u = 0.37) 2.6e-53; largest smallest-singular-value at the degree found 7.2e-56.
* **Gauge-power scan** (one state per sector, 36 sectors, smallest-degree polynomial solution of Pi^p o A o Pi^-p, A = formal adjoint, p = 0..4):
  | p | sectors with degrees == 6 - mu_i | kernel dim 1 | notes |
  |---|---|---|---|
  | 0,1,2 | 0 / 36 | — | no solution up to degree 9 in any sector |
  | **3** | **36 / 36** | 36 / 36 | degrees exactly 6 - mu_i |
  | 4 | 0 / 36 | 36 / 36 | degrees 8 - mu_i (= Pi x the p = 3 solutions: degree +2 in every sector) |
  So **p = 3** is the smallest (and the only exactly-matching) gauge. With (2,1,0) (p = 2) this gives p = 2, 3 for lam1 = 2, 3. **The rule is p = lam1 - lam3** (= lam1 here), and it is explained by the local exponents: the adjoint solutions psi have exponents 2 - e_i = {2, 0, -3} at theta_s, so Pi^-p is needed to absorb the pole of order lam1 - lam3 = e_max - 2; the direct Q = Pi^p psi then has exponents (0,3,5) (observed, all states). **Round-1 (2,1,0) is consistent: p = 2 = lam1; the candidates "p = lam1 - lam3" and "p = lam1" cannot be separated with lam3 = 0 data** (nor from |Lambda| - 1 = 2, 3, or |Lambda| - lam2); the exponent argument selects lam1 - lam3 (e_max - 2 with e_max = lam1 - lam3 + 2, assuming exponents {0, lam2 - lam3 + 1, lam1 - lam3 + 2}, which is the sl(3) shifted-weight pattern; the lam3 dependence itself is untested here).
* Local exponents at theta_s (dual): (0,2,5), 450 of 450 (state, site) cases; direct Q: (0,3,5), 450 / 450.
* **Frobenius (no-logarithm) constraints on jets (for D1)**: write x = u - theta_s, S1 = 4/x + sum s_k x^k, S2 = 4/x^2 + H/x + sum h_k x^k, S3 = a/x^2 + b/x + sum t_k x^k (s_k, h_k, t_k are the Taylor coefficients of the regular part = constants e_k(zeta) + the other site's poles). Exponents {0,2,5}: the exponent-0 series needs two conditions, the exponent-2 series one more, but the order-5 coefficient of the free constant c_2 coincides with the exponent-2 condition, so **3 independent conditions per site (6 in total on the 6 unknowns)**. Derived symbolically and verified on all 225 states x 2 sites (max relative violation 1.3e-58):
  * C1 (order x^2): **a^2 - H a + 4 b = 0**
  * C3 (order x^5, c_2 coefficient and the exponent-2 series): 96 H^3 - 104 H^2 a - 576 H^2 s0 + 36 H a^2 + 440 H a s0 - 160 H b + 384 H h0 + 1056 H s0^2 - 576 H s1 - 4 a^3 - 80 a^2 s0 + 56 a b - 144 a h0 - 432 a s0^2 + 240 a s1 + 352 b s0 - 768 h0 s0 + 384 h1 - 576 s0^3 + 960 s0 s1 - 384 s2 - 192 t0 = 0
  * C2 (order x^5, c_0 coefficient): printed in `qside_310_output.txt` (Section 1, "order x^5, coefficient of free constant c_0", needs t1, t2, h2, s1) = 0.
  (`frobenius_polys()` in `qside_310.py` returns the three polynomials with the variables (H,a,b,s0..s2,h0..h2,t0..t2).)
* Negative controls (36 states, one per sector, all to degree 9): dual ODE with e^{-zeta u}: none; plain adjoint (p = 0) with e^{-zeta u}: none; Pi^-3 gauge with e^{+zeta u}: none.
* Per-sector degree table (all 36 sectors): dual deg = mu, direct deg = 6 - mu; max residual per sector between 1.1e-56 and 2.6e-53 (see `qside_310_output.txt`).

## Q2. Wronskians and quantum Weyl

* Direct: Q_a Q_b' - Q_b Q_a' = c_ab e^{-(z_a+z_b)u} F(u) p_c (c complementary, all monic): **F = Pi^2 = (u-th1)^2 (u-th2)^2 = prod_s (u-theta_s)^{lam1-lam2}** (deg F = 4 in all cases); 225 states x 3 pairs: relative remainder of the polynomial division <= 9.2e-49, |F - Pi^2| <= 8.0e-49; **c_ab = zeta_a - zeta_b exactly** (monic normalisation, deviation 0).
* Dual: P^aP^b' - P^bP^a' = c'_ab e^{+(z_a+z_b)u} F'(u) q_c with **F' = Pi^1** (deg 2), remainder <= 7.0e-51, |F' - Pi| <= 6.6e-51, c'_ab = zeta_b - zeta_a exactly. (So F' = Pi^{|Lambda| - p}: consistent with the (2,1,0) result Pi and with psi = Pi^-p Q = W[P,P]/W.)
* Quantum Weyl: S1 = det[P,P',P''']/W, S2 = -det[P,P''',P'']/W, S3 = det[P''',P',P'']/W with W = det[P,P',P'']: max relative deviation from the eigenvalue functions S1 3.8e-47, S2 4.1e-47, S3 4.1e-47 (225 states x 3 u values). **W = c_W e^{(z1+z2+z3)u} Pi^4** (the power of Pi is |Lambda| = 4 = residue of S1; for (2,1,0) it was 3), c_W = prod_{i<j}(zeta_j - zeta_i) = -35/288 = -0.1215277778 (state independent, deviation 4.1e-47). Direct side: psi_i = Pi^-3 Q_i span the kernel of the monic operator with coefficients (-S1, 2S1'+S2, -(S1''+S2'+S3)): deviations 1.4e-40, 5.1e-40, 7.3e-39 (numerical derivatives).

## Q3. Brackets, orthogonality, det N, ratio identity (gauge p = 3)

* Lagrange identity: with psi = Pi^-3 Q_i^A, P = P^j_B: P (dS2 psi)' + dS3 P psi = C'; the double pole of S2 is the scalar c2 = 4 (cancels in dS2), S3 has no triple pole, so Res_{theta_s} gives **the same six columns** as for (2,1,0): N[(s,i,j),H_t] = -Res[psi P'/(u-th_t)], N[.,a_t] = Res[psi P/(u-th_t)^2], N[.,b_t] = Res[psi P/(u-th_t)] (pole order p + m <= 5, Taylor to order 4), and sum_m Delta_m N_m = 0, Delta = (dH1,dH2,da1,da2,db1,db2).
* Identity residual: max |sum Delta N| / sum |Delta N| = **5.5e-56** over 600 random ordered pairs A != B x 18 rows.
* Rank tally of the stacked 18x6 (columns normalised, rank = #{sigma_k/sigma_1 > 1e-40}), 600 random pairs A != B: **same sector: rank 5 (28/28); different sectors: rank 5 (439), rank 4 (133)**. max sigma_6/sigma_1 = 1.6e-57, min sigma_5/sigma_1 among rank-5 pairs 9.8e-12, max sigma_5/sigma_1 among rank-4 pairs 5.5e-57 (clean gaps). **A = B: rank 6 for all 225 states** (min sigma_6/sigma_1 = 7.1e-12). Rank 4 occurs exactly as for (2,1,0): see R (tally by sign pattern of mu_B - mu_A: (1+,2-): 133 pairs all rank 4; (1+,1-): 303 rank 5; (2+,1-): 136 rank 5; (0,0): 28 rank 5).
* det N (6x6, rows (s,(i,j)), (i,j) in {(1,3),(1,2),(2,3)}, columns H1,H2,a1,a2,b1,b2), |det N|/prod(row norms): A != B (150 pairs) max **5.9e-65**; A = B (all 225 states) min 2.0e-37, max 5.4e-5.
* x0-free ratio identity R_B(w)/R_B(w') = [D_B(w)/D_B(w')][D_B0(w')/D_B0(w)], B0 = e1 x e1 (index 224, weight (6,2,0)): **max relative deviation 6.8e-52 over ALL 225 states** (all non-degenerate; worst states (3,3,2), (2,3,3), (3,2,3), (3,1,4), (1,3,4)); e.g. B = 0: 13.7396044771896 both sides.
* 84-triple scan (same triple on both sites): **72 pass (< 1e-30), 6 fail (deviation 1.3 - 1.4), 6 degenerate**; the failing are exactly the 6 permutation triples ((1,1),(2,2),(3,3)), ((1,1),(2,3),(3,2)), ((1,2),(2,1),(3,3)), ((1,2),(2,3),(3,1)), ((1,3),(2,1),(3,2)), ((1,3),(2,2),(3,1)); the degenerate ones are the 6 triples with all i equal or all j equal (det identically zero). Same as (2,1,0); the XXX triple (1,3),(1,2),(2,3) passes.
* Controls (must fail, do): a2 column replaced by Res[psi P/(u-th2)^3]: deviation 1.96; orientation swapped (vacuum dual in ket slot, state direct in bra slot): 1.64; **gauge weight Pi^-2 instead of Pi^-3 in the functionals: 1.70; Pi^-4: 1.91**.

## Q4. Regime (zeta = eps phi, omega = eps om^, omega' = eps om^'; eps = 1e-6, 2e-6; 100 digits; eigen-residual 1.7e-97)

* 10 HW branches picked by the smallest raising residual (same sectors as the spec): eps = 1e-6 residual vs the next state: (6,2,0) 0 | (6,1,1) 2.3e-7 vs 1.41 | (5,3,0) 1.9e-7 vs 2.83 | (5,2,1) 2.25e-7, 2.25e-7 vs 2.41 | (4,4,0) 4.5e-7 vs 2.0 | (4,3,1) 3.5e-7, 3.5e-7 vs 2.41 | (4,2,2) 4.0e-7 vs 2.28 | (3,3,2) 6.6e-7 vs 1.75. Clean gaps of 6-7 orders; the residuals double at 2e-6.
* Degrees at the regime (kernel threshold 1e-50): dual = mu and direct = 6 - mu (e.g. (6,1,1): dual (6,1,1), direct (0,5,5)), kernel dim 1. Vacua: eigencovector ok, direct Q degrees (0,4,6). (Pitfall noted: with the default 1e-35 kernel threshold the 100-digit regime gives spurious lower degrees, e.g. (4,2,0) for the (6,2,0) dual; the threshold must be set to 1e-50, done in `q4_point`.)
* Ratio identity on the 10 branches: max deviation **6.6e-52** (eps = 1e-6) and 8.9e-54 (eps = 2e-6) (0.0 for the top branch by construction).
* eps-exponents of det N (f ~ eps^p from 1e-6 / 2e-6): **|det N(w)| ~ eps^16 for (6,2,0), (6,1,1), (5,3,0), (4,4,0); eps^15 for both (5,2,1) and both (4,3,1); eps^14 for (4,2,2), (3,3,2)**, identical for omega'; (round 1 (2,1,0): 12, 11, 10); det/prod(col norms) ~ eps^{38,40,38,39,39,38,39,39,38,38}; sigma_6/sigma_1 of the column-normalised N ~ eps^{11,11,11,11,11,11,11,11,10,10}; max |N(w)| ~ eps^-4; |R(w)| = |<Omega x Omega|Psi_B>| ~ eps^{0,1,1,1,1,2,2,2,2,3}. For both multiplicity-2 sectors the two branches have det N equal up to sign to 6 digits (e.g. (5,2,1): 5.04230709e-98 vs -5.04230858e-98) and equal ratio-test quality.

## R. The rank-4 puzzle of (2,1,0)

Data: all 64*63 = 4032 ordered pairs (A bra, B ket) of P1's (2,1,0) states recomputed (gauge p = 2, 60 digits).
* **Tally: rank 5: 3318; rank 4: 714; same-sector pairs always rank 5 (282/282).** The rank is a function of Dmu = mu_B - mu_A ONLY (61 distinct Dmu, each with a single rank; 87 ordered sector pairs have all their pairs at rank 4, no sector pair mixes).
* **Characterisation: rank 4 <=> Dmu has exactly one positive and two negative components** (sign pattern (+,-,-) in any order, e.g. (2,-1,-1), (3,-2,-1), (4,-2,-2), (-1,-1,2), (-2,3,-1), ...). Rank 5 for the patterns (0,0), (+,-,0), (+,+,-). In terms of "which column vanishes identically" or "dH1 = dH2": not the case (no column is zero; the rank is not about the weights differing by a root: (1,-1,0)-type roots give rank 5, (2,-1,-1) gives 4).
* **Mechanism (verified)**: for every k with Dmu_k <= -1 the product psi_k^A P^k_B is a RATIONAL function (the exponentials e^{-zeta_k u} e^{+zeta_k u} cancel) of degree Dmu_k <= -1, so the sum of its residues (with any of the six weights 1/(u-th_t)^m) vanishes: **row(1,k,k) + row(2,k,k) = 0 in all six columns**. Max relative violation 5.4e-58 over the 4464 (pair, k) cases with Dmu_k <= -1; for Dmu_k >= 0 the relation fails (smallest violation 1.2e-6). Each site's 9 rows span exactly a 3-dim space S_s (rank (3,3) for all 400 pairs checked; locally only the exponent pairs (-2,0), (-2,2), (0,0) of psi x P give poles), so rank = 6 - dim(S_1 n S_2); each negative component of Dmu puts one independent row vector into S_1 n S_2: **one negative component (or Dmu = 0, where the corresponding relation is the residue-at-infinity relation) gives rank 5 (this null vector is the Lagrange vector Delta); two negative components give two independent elements, rank 4**. Three negatives are impossible (sum zero). Same rule for (3,1,0): of 600 random pairs the (1+,2-) pairs are exactly the 133 rank-4 pairs (see Q3).
* **Second null vector n2**: for rank-4 pairs the null space is 2-dim and contains Delta. n2 is NOT a fixed linear form in the data tested: not a unit vector (best fit residual 0.32), not dataA = (H1,H2,a1,a2,b1,b2)^A or dataB (residuals 0.2 - 0.8), not Delta restricted to one site (0.1 - 0.9; the plane differs from state pair to state pair within the same sector pair). I did not find a closed form (two attempts, stopped as instructed); what is established is that n2 spans, with Delta, the orthogonal complement of the rowspace S_1 + S_2 (4-dim).

## Surprises / open points

1. c2 = 4 for (3,1,0), not 6: c2 = lam2 (lam1 + 1) (lam3 = 0), not the Casimir; local exponents {0,2,5} = {0, lam2+1, lam1+2}.
2. Nilpotency index of the vacuum exponential series is 4 for E12, E23, E13 (the spec expected 3).
3. The gauge power is p = 3 = lam1 - lam3 (exponent argument); (2,1,0) is consistent. lam3 = 0 in both, so lam1 vs lam1 - lam3 is not decided by data; also dual Wronskian F' = Pi^(|Lambda|-p) = Pi^1.
4. W = Pi^|Lambda| (power 4) but the Q-Wronskian F = Pi^(lam1-lam2) = Pi^2 and F' = Pi.
5. The Frobenius conditions are 3 per site (6 equations in the 6 unknowns): the Laurent data are cut out by a polynomial system; their number of solutions vs the 225 states was not counted.
6. Mixed-weight pairs of the same sign-pattern class as in (2,1,0) have rank 4; n2 unidentified.
7. 100-digit regime solver needs the kernel threshold 1e-50 (1e-35 gives wrong degrees).
8. The mp eigen solver returns complex numbers with 1e-59 imaginary parts in some sectors (e.g. (4,2,2)); they are cleaned in `build_states`.

## Not done

* Closed form of n2 (see R); dependence of the gauge rule on lam3 (needs a rep with lam3 > 0, e.g. (3,2,1) or (2,1,1); (2,1,1) has a triple pole in S3, so the Laurent structure changes).
* Number of solutions of the Frobenius system vs states; the (4,3,1) CG-mixing functionals were not part of this task.
