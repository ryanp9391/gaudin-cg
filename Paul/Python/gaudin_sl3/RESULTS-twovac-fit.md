# RESULTS-twovac-fit: direct fit of the Q-side formula for the two-vacua overlap of the sl(3) Gaudin chain (agent D1, 2026-10-16)

Files: `twovac_fit.py` (`python -X utf8 twovac_fit.py [--stages F0,F1,F2,F3,F4,F5]`, about 2.5 min with warm caches, 80 digits for F0-F3 and 120 for F4; independent of Wolfram),
`twovac_fit_output.txt` (full printed output of the last run; every number below is copied from it), caches `twovac_fit_*.pkl` (delete to recompute).
Python only. Imports `gaudin_chain.py`, `qside.py` (Q solvers, `Nrows`/`build_N_det`, `raising_residual`) read-only. The 25-parameter linear least squares uses an own fixed-point
(320-bit integer) column-pivoted Gram-Schmidt QR (rank threshold 1e-40 relative to the largest pivot); a 512 x 400 problem takes 30 s.

## 0. Objects, orders, normalisations

* Basis order: **gl3_reps order, HW ket index 7** (as `gaudin_chain`/`qside`). Notebook (descending-GT) index a maps to gl3 index NB2G = [7,5,6,4,2,1,3,0]; this is **not** a reversal.
  The x0^G direction of RESULTS-xxx-limit X1c, mapped with NB2G, is reproduced to 1e-81.
* Phi_inf = Omega^inf[zeta,w1] (x) Omega^inf[zeta,w1'], r_k^inf as in the plan; against `vacuum_covector_1site` at omega = (w1, 1e6, 3/5): projective difference 2.2e-6, 9.0e-7 (O(1/omega2)).
* States: all 64 per zeta point (eigenvectors of S2(3/7) + (pi/3) S3(-2/9), 80 digits, eigen-residual < 4e-77; Q's from P1's solvers, ODE residual 1.2e-74). Points: P, P' and
  three extras of the plan (X1, X2, X3), further X4..X6 and 26 real-spectrum random rationals R1..R26 (for the large fits). For some zeta triples the spectrum of S2 + pi S3 is complex: such triples are skipped.
* x0^G: both my Gauss-factor construction and D2's closed form x^G(zeta) = (z2-z1, z3-z1, z2-z1, 2z3-z1-z2, 3(z2-z1), 3(z3-z1), 2(z3-z2), 3(z3-z2)) (notebook order): the overlaps agree
  to 1e-77 at every zeta point and for every state (ratio 1.0 at P, P', X1, X2, X3). So **x0^G is normalised exactly as D2 and W2** (the frame-tie scalar kappa^G does not hide in x0).
* **Normalisation used for the fits: Y0 := 4 c(zeta,w1,w1') Qnum(B), c = (z2-w1)(z3-w1)^2(z2-w1')(z3-w1')^2(z1-z2)^2 (= w1^3 w1'^3/(kappa^G kappa^G') with W2's kappa^G = Delta1 Delta2 = w^3/((z2-w)(z3-w)^2(z1-z2))).**
  Qnum is divided by kappa^G kappa^G' (up to the zeta-independent w^3 w'^3), i.e. the companion-frame number, as W2 prescribed. Dividing only by kappa kappa' and using the Q1 P3 factor (the plan's M)
  does NOT make the kernel zeta-independent: A_s := (z3-z2) e^{sigma theta_s} Pi'(theta_s) (p02 = A_s Q1) carries a zeta-dependent factor.
  **Per-zeta scalar check:** per-zeta-point fits (64 rows, 25 unknowns, rank 25) give the same kernel as the 5-point stacked fit to 2e-73, 2e-73, 2e-73, 4e-73, 3e-72 (P, P', X1, X2, X3): no free scalar per zeta point is needed once c and M = P3 p02 are used.

## F0. Sanity

* <x0^G|Psi_B> != 0 for all 64 states at P (max-component normalised kets: 1.7e-4 <= |.| <= 8.3); <Phi_inf|Psi_B> at (5/2,7/4): 0.0139 <= |.| <= 8.5, no vanishing state.
* The 64 Qnum(B) at P, (5/2,7/4), with (S2(3/7), S3(-2/9)) labels to 30 digits: in `twovac_fit_output.txt` (F0 table; e.g. top state (4,2,0): 16 exactly, (2,4,0): 6.4, (4,0,2): 5.034965034965...). **They agree with W2's A5-1 table (`Gaudin-L2-chain-2.wb`, [2,1,0] run) for both (5/2,7/4) and (-2/5,3): max relative difference 3.3e-31 and 2.8e-31 (= the 30 printed digits), all 64 states matched by the eigenvalue labels.**
* T(B) = <Phi_inf|Psi> D_B(w0)/<Omega_D(w0)(x)Omega_D(w0)|Psi> (w0 = (5/2,7/3,3/5), D_B = P1's det N): **T/Qnum is NOT state independent** (state spread 1.01 at P, 2.6e4 at P'). The reason (log-regression over 64 states, residual 2e-77): <x0|Psi_B>/n_B (n_B = R_B/D_B) = const x Q1(th1)P3(th1)Q1(th2)P3(th2) exactly (exponent 1 on each, 0 on the other 10 numbers).
  Hence **T(B) = c_T(zeta) Qnum(B) M1 M2 with state-independent, zeta-dependent c_T** (spread 1e-75; c_T/(Qnum P3p02 P3p02) at P = 0.0797567, at P' = 1.3058...). T and "Qnum x M1 M2" are the same target; the "T without M1 M2" and "Qnum with M1 M2" columns of the plan coincide. Qnum alone is used below.

## F1. Jet facts (point P, 64 states, 2 sites; 80 digits)

(a) p01 = W(P2,P3)(theta_s) = 0: max |p01|/(|J2_0 J3_1| + |J3_0 J2_1|) = 1.9e-78. p02 = W'(theta_s) = (z3-z2) e^{sigma theta_s} Pi'(theta_s) Q1(theta_s): max relative deviation 1.4e-76.
(b) Frobenius: P^j'(theta_s) = kappa_s P^j(theta_s) with the SAME kappa_s for j = 1,2,3: **kappa_s = a_s/3, a_s = coefficient of (u-theta_s)^-2 in S3** (max relative deviation over states, sites, j: 2.7e-78).
Derivation: indicial polynomial rho(rho-1)(rho-2) - 3 rho(rho-1) + 3 rho = rho(rho-2)(rho-4). For P = c0 + c1 x + ..., x = u - theta_s, the x^-2 coefficient of (d^3 - S1 d^2 + S2 d - S3)P is I(1) c1 - a_s c0 with I(1) = 1(-1)(-3) = 3:
only the leading singular parts (3/x in S1, 3/x^2 in S2) and the double pole a_s/x^2 of S3 enter; H_s, b_s and the regular parts of S1, S2 first appear in c2. (The spec's expectation "in terms of regular parts of S1, S2, H_s, a_s" is refuted: kappa depends on a_s alone.) E.g. kappa = 3/2 for the (0,2,4), (0,3,3)_1, (0,4,2) states; the rest in the output.
(c) 3x4 jet matrix [J^j_k]_{k=0..3}: sigma3/sigma1 <= 4.5e-79 (rank 2); 3x5 (k=0..4): sigma3/sigma1 >= 5.0e-10 (rank 3). Structure: the exponent-4 solution has J_0..J_3 = 0, the exponent-2 solution has J_0 = J_1 = 0, all three have J_1 = kappa_s J_0.
Dependent monomials across states: per site, the 12 {p_kl P3^(c)} (k<l<=3, c<=1) have rank 8 (the two p01 P3^(c) vanish; the other two relations of the 4 lost features were not identified individually); 24 -> rank 16 (c <= 3); the Q1-basis {Q1^(k) P3^(c)} (k<=3, c<=3) has full rank 16.

## F2. The hypothesis ladder (base point (w1,w1') = (5/2,7/4), theta = (2/7,-3/5); rows = 64 x (5 or 8 zeta points); target Y0 = 4 c Qnum, "M" = times (P3 p02)(th1)(P3 p02)(th2))

Residual = relative ||A k - y||/||y|| (row-normalised; 0.0 means below the 1e-90 fixed-point floor). A_s Q1 = p02, A_s := (z3-z2) e^{sigma theta_s} Pi'(theta_s) (needed so that the direct basis is on the same footing).

| level | basis per site (same on both sites) | unknowns | rank | null | rows | resid with M | resid without M |
|---|---|---|---|---|---|---|---|
| a | **N5 = {P3 p02, P3' p02, P3 p03, P3' p03, P3 p23}** (D2) | 25 | 25 | 0 | 320 | **0.0** | 0.966 |
| a | A_s Q1^(k) P3^(c), k<=2, c<=1 (direct, 6) | 36 | 36 | 0 | 320 | 0.294 | 0.944 |
| a | A_s Q1^(k) P3^(c), k<=1, c<=3 (8) | 64 | 64 | 0 | 512 | 0.562 | 0.921 |
| a | D2's nine direct monomials with kappa = P3'/P3 | 81 | 81 | 0 | 320 | 0.363 | 0.871 |
| b | Plucker {p_kl P3^(c)}, k<l<=3, c<=1 (12) | 144 | 64 | 80 | 320 | **0.0** | 0.916 |
| c | Plucker k<l<=3, c<=3 (24) | 576 | 256 | 320 | 512 | **0.0** | 0.692 |
| c | A_s Q1^(k) P3^(c), k<=3, c<=3 (16) | 256 | 256 | 0 | 512 | 0.0643 | 0.747 |
| d | N5 with K = K0 + sum_m Delta_m K_m (six Laurent data) | 175 | 119 | 56 | 512 | 0.0 | - |

* **The first level that fits is the 5-element basis of D2**: residual 0 (< 1e-90) over 5 zeta points (320 rows), rank 25 (all 25 products independent), null space 0, so the kernel is **unique** in this basis. (Earlier, with the zeta-independence assumed and the plan's M = Q1 P3, nothing fitted at any level up to 400 unknowns on 512 rows and even the full degree-(1,2) jet space {J2_a J3_b J3_c, a,b,c<=3} (rank 36/site, 1295 unknowns on 1792 training rows) gave train residual 4e-2, hold-out error 8e2: that was the c(zeta) / A_s(zeta) normalisation problem, not a missing basis.)
* The twelve-element Plucker basis fits because it contains N5 (rank 8 per site, 64 independent products, null space 80 = the 144 - 64 relations among the 12 x 12 products). Minimal-norm representative: |K|_F = 471, support 49 of 144, all entries rational; sparse representative = the N5 kernel (25 nonzero entries, 5 of the 12 Plucker elements per site) - minimal because the 25 products are independent.
* The direct bases do NOT fit with a zeta-independent kernel (0.29, 0.56, 0.36, 0.064). Reason (checked): with W = (z3-z2) e^{sigma u} Pi Q1, p03 = A_s[(Lambda_s - kappa) Q1 + 2 Q1'] with Lambda_s = 2 sigma + 2/Pi'(theta_s): the coefficient contains sigma = z1+z2+z3. At a single zeta point the nine monomials express all five N_i (residual 0.0 for all 5 N_i, both sites) but the coefficients depend on sigma. The Plucker/N5 form is the zeta-free one.
* Level d: the kernel is not itself a Q-number: the K0-only fit has residual 0 at rank 25. In the extended design the minimal-norm K_m are non-zero (|K_m| = 12-19 vs |K0| = 652, the extended design has a 56-dim null space) only because the Laurent data are tied to the jets (kappa_s = a_s/3, H1 + H2 = const - W).
* **Kernel at the base point (basis order N = (P3 p02, P3' p02, P3 p03, P3' p03, P3 p23), thi = 1/theta12 = 35/31), rows i, columns j; all 25 entries are rational (reconstructed from 80-digit values with denominators <= 1e12; max |fit - D2's printed G| = 2.8e-74)**:

```
 i\j      1                      2                      3                     4                      5
 1   -2786875/123008      698893125/1906624       -361375/1984        453618375/1847042      -23065925/119164
 2   -609539175/3813248   -20215755/61504         5269159875/59105344  -5815005/238328         218655/3844
 3    77175/3968          -5328599325/59105344     1225/64             -4894925/238328          96825/3844
 4    981927415/29552672   34151355/476656         -7935795/476656      -2275/1922               -295/31
 5   -6878473/953312      -157773/15376            43169/15376          63/62                      1
```
* **Controls (must fail; all do)**, N5-shaped basis, same target: lone index P2 instead of P3: resid 0.988; lone index P1: 0.989; both factors built from the site-1 jets: 0.365 (rank 14); both from the site-2 jets: 0.58 (rank 14); jets at the wrong points u = (0, 1/3): 0.482; right basis but target without M1 M2: 0.966.
* T(B) as target: T = c_T(zeta) Qnum M1 M2 (F0), so with one free scalar per zeta point it is the same fit; not repeated.

## F3. Kernel laws (basis N5, 5 zeta points per fit, rank 25, residual 0.0 everywhere)

* zeta independence: per-zeta fits equal the stacked kernel to <= 3.3e-72 (see section 0).
* 3x3 grid w1 in {5/2,-2/5,9/7} x w1' in {7/4,3,-5/3}: kernel vs D2's analytic G: max |fit - G| between 7e-75 and 3.2e-72 at all nine points; all 25 entries rational at each point.
* Polynomial dependence on (w1,w1'): each G_ij fitted as a polynomial of total degree <= 6 on the 72 ordered pairs of 9 values (w1 != w1'; 28 coefficients), residual 0.0, all coefficients rational, **identical (exact rational comparison) to D2's analytic G_ij(w1,w1',1/theta12) in 25/25 entries**. No poles: every entry is a polynomial in (w1,w1') (none in w1-w1', none in w1-zeta_j); the maximal total degree is 6 (G_11 contains 4 w1^3 w1'^3).
* theta dependence: five variants (1/5,-4/9), (3/7,-1/4), (1/2,-1/3), (-1/7,-6/5), (4/5,-2/3) (theta12 = 29/45, 19/28, 5/6, 37/35, 22/15) at (5/2,7/4): each fit matches D2's G to <= 2.3e-73; with the base theta (6 values) each entry is a polynomial in thi = 1/theta12 of degree <= 4 (5 coefficients, one redundant check, residual 0.0), **identical to D2's polynomial in 25/25 entries**. Poles only at theta12 = 0 (polynomial in 1/theta12).
* Scaling u -> c u, g -> g/c, c = 2 (data at 2 theta, zeta/2, omega/2): kernel fit = D2's G at (w/2, 2 theta12) to 5e-72; entrywise G'_ij/G_ij = 2^-(10-w_i-w_j) with w = (2,3,3,4,5) exactly (table in the output; matrix entries 1/64, 1/32, 1/32, 1/16, 1/8 / ... / 1 at (5,5)); D2's homogeneity G(w/2,w'/2,thi/2) = 2^-(10-w_i-w_j) G holds exactly (True). So **the kernel depends on (w1,w1',thi) homogeneously of degree 10 - w_i - w_j, no jet-rescaling factor**.
* epsilon-order of each entry under (w1,w1') -> eps (w1,w1') (lowest total degree of the polynomial in w1, w1'; thi fixed):

```
  4 2 4 0 0
  2 2 0 0 0
  4 0 4 0 0
  0 0 0 0 0
  0 0 0 0 0
```
  (the entries are not homogeneous in (w1,w1') alone at fixed thi; this is the lowest power).
* XXX cell-62-style law: Qnum_G(B) = (1/(4c)) sum_ij G_ij nu_i(th1) nu_j(th2), nu = (1, kappa, pi, kappa pi, eps) with pi = p03/p02, eps = p23/p02. The fitted G is exactly D2's (printed in `RESULTS-twovac-limit.md`); I did not independently re-derive a closed form beyond the polynomial reconstruction above, which reproduces all 25 printed entries.

## F4. Untwisting regime (zeta = eps phi, phi = (1/3,-2/7,5/11); w = eps w^, w^1 in {3/4,-5/3}, w^1' in {-12/5,7/2}; eps = 1,2,4 e-6; 120 digits; eigen-residual 2e-118)

Branches: smallest raising residual per HW sector (4,2,0): 0, (4,1,1): 2.3e-7 [next 1.41], (3,3,0): 1.9e-7 [1.41], (3,2,1) x2: 2.2e-7, 2.2e-7 [2.05], (2,2,2): 3.5e-7 [1.37] at eps = 1e-6 (residuals double with eps). Q degrees as at finite twist, kernel dimension 1 for all six Q's of all branches.

| eps | max over the 4 grid points of \|Qform/Qnum_chain - 1\| for branches (4,2,0), (4,1,1), (3,3,0), (3,2,1), (3,2,1), (2,2,2) | cancellation in sum G N N (digits) |
|---|---|---|
| 1e-6 | 1.2e-90, 5.2e-98, 4.1e-91, 6.3e-97, 1.4e-96, 4.5e-91 | 30.8, 23.6, 25.3, 24.4, 25.1, 13.6 |
| 2e-6 | 1.6e-92, 5.7e-99, 3.1e-91, 2.9e-98, 9.3e-98, 2.9e-91 | 29.3, 22.4, 24.0, 23.2, 23.9, 13.0 |
| 4e-6 | 3.5e-94, 2.4e-100, 5.1e-93, 1.7e-99, 1.3e-98, 3.5e-93 | 27.8, 21.2, 22.8, 22.0, 22.7, 12.4 |

Qform = (1/(4c)) sum G_ij N_i(th1) N_j(th2)/(M1 M2) uses only the Q's (jets of the six branches' dual Q's); the chain value uses <Phi_inf|Psi>/<x0^G|Psi>. The deviation is at the 1e-90..1e-100 level: the formula holds in the regime with a loss of 13-31 digits to cancellation (the largest |G_ij N_i N_j| is eps^2..eps^3, the sum eps^(4..8)), well inside 120 digits.
eps-exponents p (from eps = 1e-6 -> 2e-6 | 2e-6 -> 4e-6, equal to the digit in both columns), first grid point, branches (4,2,0), (4,1,1), (3,3,0), (3,2,1), (3,2,1), (2,2,2):
sum_ij G N N: 8, 7, 7, 6, 6, 4; largest single term |G N N|: 3, 3, 3, 2, 2, 2; c: 8 (all); M1 M2: 2 (all); N_max (site 1): 1 (all); p02(1): 1 (all); p03(1): 1, 2, 1, 1, 1, 1; p23(1): 3, 3, 3, 2, 2, 1; P3(th1): 0 (exponential, no eps power); Qnum_chain: -2, -3, -3, -4, -4, -6; <x0|Psi>: 2, 3, 3, 4, 4, 6 (so Qnum <x0|Psi> = eps^0: the overlap <Phi|Psi> is O(1), as W2's S = O(1)). Kernel entries: see F3 table.

**CG from the Q side alone.** Procedure: for each branch the Q-side numbers Qform at the four grid points (eps = 1,2,4 e-6) give rho_a = Qform_a/Qform_0; Richardson (eliminating eps, eps^2) gives rho_a(eps -> 0). With Psi_B -> v_B in the exact HW space H_mu of the untwisted diagonal action and Phi_a eps-free, <Phi_a|v_B> = Qnum_a <x0|v_B>, so (f_a - rho_a f_0)(v_B) = 0 (f_a(v) = <Phi_a|v>): 3 equations on dim H_mu = 1 or 2 unknowns. Solve, G-normalise (Gram diagonal of the plain GT generators, d_a E_ij[a,b] = E_ji[b,a] d_b; sign: first nonzero component in descending GT order positive), read the slice <a,hw|v> (second factor = HW ket, gl3 index 7).

| branch | dim H | consistency (smallest singular value of the 3 x dim system, relative) | rho_a (Q side) | slice of the Q-side vector (first-factor gl3 index a: component) | exact value (W1) |
|---|---|---|---|---|---|
| (4,2,0) | 1 | 5.8e-90 | 1, 1, 1 | a=7: 1.0 | 1 |
| (4,1,1) | 1 | 6.0e-20 | 0.818289504857, -0.0324263038549, -0.214136798998 | a=5: -0.707106781187 | -1/sqrt2 |
| (3,3,0) | 1 | 2.8e-19 | 0.487571129081, -0.174603174603, -0.687032045523 | a=6: -0.707106781187 | -1/sqrt2 |
| (3,2,1)_1 | 2 | 3.6e-20 | 0.733223361576, 0.02053459544, 0.157867160066 | a=2: 0.208887289734, a=4: 0.049311600013 | (0.0493116000, 0.2088872897) |
| (3,2,1)_2 | 2 | 5.7e-19 | -0.16434558767, -0.0194039421019, 0.129004381313 | a=2: -0.0797878448606, a=4: 0.337986734608 | (0.3379867346, -0.0797878449) |
| (2,2,2) | 1 | 1.5e-18 | -0.373964006175, -0.000788598886267, 0.0192070080153 | a=0: 0.117851130198 | 0.1178511302 |

The Q-side rho_a equal the chain's rho_a to all printed digits (12). The consistency singular values 1e-18..1e-20 are the O(eps^3) Richardson error (the (4,2,0) one is exact: rho = 1). The Q-side vectors (G-norm 1) reproduce W1's exact values: -1/sqrt2, -1/sqrt2, 0.1178511302 and both octet vectors to the printed digits: the slice entries (a=4, a=2) are (0.049311600013, 0.208887289734) and (0.337986734608, -0.0797878448606), W1's pairs (0.0493116000, 0.2088872897), (0.3379867346, -0.0797878449). Against the chain-side eigenvectors (eps = 1e-6) projected on H_mu with the same normalisation: max difference 0 (1-dim top and singlet), 1e-121 (10, 10bar), 1.1e-9 (octets, O(eps^2..eps) from the chain vector's own HW leakage |v - Pv|/|v| = 8e-8). Octet mixing: the Q-side directions are the two eigenvectors of S3c[1,2] on the HW space (the mixing roots 5.854101966, -0.854101966 were not recomputed here; the two Q-side vectors coincide with the chain-side copies and with W1's numbers, which is the same statement in components).
The x0-based normalisation drops out of rho_a; absolute normalisation of v is fixed by the G-norm of the algebraically determined HW vector (no chain eigen-solve enters the Q-side vector).

## F5. Cross-checks

* **D2 (`RESULTS-twovac-limit.md`)**: D2's five monomials, constant 4c and the vacuum factor M = P3 p02 are confirmed by an independent fit; the fitted kernel equals D2's printed G at the base point (2.8e-74) and at the 8 other points of the 3x3 grid, the 72-point polynomial reconstruction (28 coefficients) equals D2's symbolic G_ij in 25/25 entries (exact rationals), 5 theta variants (25/25 identical polynomials in 1/theta12) and the c = 2 scaling (homogeneity degree 10 - w_i - w_j exact). G rank 5 is consistent with the 25 independent products. D2's statement that only 5 of the 12 Plucker elements occur is the null-space statement of level b.
* **W2**: the 64 Qnum(B) at P for (5/2,7/4) and (-2/5,3) agree to 3e-31 (the printed digits); W2's undone covector closed forms a', c', g' were not needed (I used Phi_a directly, eps-free); W2's CG-completeness and K/S = O(eps) statements are consistent with section F4 (the Q-side CG reproduced W1's values).

## Surprises

1. Kappa_s = a_s/3 (only the S3 double pole), not a function of H_s, regular parts.
2. The plan's M = Q1 P3 and "kernel zeta-independent" are incompatible: p02 = A_s Q1 with A_s zeta-dependent, and the direct-Q1 bases carry a hidden sigma = z1+z2+z3 in Lambda_s; with M = P3 p02 and the scalar c(zeta,w,w') the kernel is exactly zeta-free (all zeta points, both normalisations of x0 equal).
3. T/Qnum is not constant: <x0^G|Psi> = n_B x const x M1 M2 (monomial), so T = c_T(zeta) Qnum M1 M2.
4. The (inf,inf) vacua span only a 4-dimensional subspace of the one-site 8-dimensional rep (rank of {Omega^inf[zeta,w]}_w = 4), yet the family is CG-complete on HW slices (F4).
5. Fixed-point QR is ~6x faster than mpmath mpf for 500 x 400 problems. A numerically-vanishing feature (p01 = W(theta_s) = 0 at the 1e-78 level) must be zeroed explicitly or it is counted as an independent feature.
6. With the 5 monomials the fit is rank-25-exact; larger independent-feature spaces (V3 x V3 of rank 1296) show no additional structure beyond N5 once c and M are right (not re-run after the normalisation fix, not needed).

## Not done

* The plan's T(B)-target fit with one free scalar per zeta point and the L1-type sparse selection beyond the N5 identification (N5 is already minimal).
* An independent derivation of closed forms of G beyond polynomial reconstruction (D2 supplies them; they were reproduced exactly).
* Octet mixing roots R^2 - 5R - 5 = 0 were not recomputed explicitly (equivalent information: the two Q-side octet vectors equal the chain copies and W1's numbers).
* Plan-level variants (P1 jets, jets at a third point) were not needed since the base level fit; "P1 jets" appear only as the lone-index control.
