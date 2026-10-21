# RESULTS-jet-basis-sym (agent M2, round 5 WP-A2, 2026-10-16)

Files: `jet_basis_sym.py` (library; J0 dictionary, companion frame, Y reconstruction), `jet_basis_sym_*.pkl` (Y, A, Theta caches), scratch drivers `t_*.py`.
Python only, exact (Fractions/sympy); numerics only for the 80-digit cached eigenstates (`general_kernel_pt_*`), rationally reconstructed (0 failures).

## §0 Dictionary and H0 (J0) — H0 CONFIRMED exactly
* Realisation: V = C[x1,x2,x3]_L, E_ij = x_i d_j (column convention, E x^alpha = alpha_j x^{alpha+e_i-e_j}). GT basis vector v_p = s_p x^{alpha(p)} with alpha(p) = weight of the GT pattern
  (weights are multiplicity-free for [L,0,0]); s_p fixed by matching ALL gl3 matrix elements of `GL3Rep` (checked for L=1..4): HW vector <-> x1^L with s=1; e.g. L=1: s=(2,1,1) for (x3,x2,x1); L=2: s=(12,6,6,2,2,1) for
  (x3^2,x2x3,x1x3,x2^2,x1x2,x1^2); L=3: (144,72,72,24,24,12,6,6,3,1). (`dictionary(L)` returns alpha, s, hw index.)
  Coefficient dictionary: a ket psi = sum psi_p v_p is the polynomial with coefficients c_alpha = s_p psi_p.
* rho acts as (rho(g)f)(x) = f(x g); the covector e_hw^T rho(g) is "f -> f(e_1 g)" (coefficient of x1^L = value at (1,0,0)).
* **H0:** Omega_c^inf[omega] = e_hw^T rho(e^{aE12} e^{aE23}) (a = 1/omega) is EXACTLY evaluation at the conic point (1, a, a^2): component p equals s_p a^{alpha2+2 alpha3}
  (verified exact Fractions L=1..4, two values of a). No transpose, no signs. x0^c = e_hw^T = "coefficient of x1^L y1^L" = Psi_c(e1;e1).
* Frame tie checked: Omega_c[a=1/w] P_phys^{-1}(zeta) = Omega_D^inf[zeta,w] / kappa, kappa = Delta_1^{L} = (w^2/((z2-w)(z3-w)))^L exactly (L=1,2), with
  P_phys^{-1} the one-site matrix built as in `general_kernel.x_cov` (`pinv_rows`), x0^G = e_hw P_phys^{-1} = x_cov.

## §1 Y exactly for [1,0,0], [2,0,0] (and [3,0,0], [4,0,0]); t-structure (J1)
Method: Y = W Psi_c^{-1} at 80 digits from the cached states (tags P, Pp, X1; theta12 variants), rationally reconstructed (denominators <= 1e12; **0 failures** for L=1,2,3 at every point).
Companion kets Psi_c = (P_phys^{-1})^{(x)2} Psi built from `pinv_rows` (same recipe as `general_kernel.x_cov`); x0^c = e_hw(x)e_hw.
* **zeta-independence** exact (identical rationals at P, Pp, X1; L=1,2).
* **Y depends on theta only through t = 1/theta12 and is a POLYNOMIAL in t of degree exactly L** (L=1: degree 1, L=2: degree 2, L=3: degree 3; fits on 7-11 distinct t with 3-4 held-out checks, exact).
  Y(t) = Y_0 (1 + t A_1 + t^2 A_2 + ... + t^L A_L) with det Y = const, T(t) = 1 + ... unipotent (T^{-1} polynomial of degree L as well).
* **t = 0**: Y_0 = y (x) y (product of one-site jets): for L=1: <01|=c11 (coefficient of x1 y1 ...); one-site jet covectors in the monomial coefficient basis c_alpha
  (c11 = coefficient of x1^2 etc.)
  - L=1:  j01 = c1, j02 = -c2, j12 = c3.   [c_i = coefficient of x_i]
  - L=2:  j01 = c11, j02 = -c12/2, j03 = c22 - c13/2, j12 = c13/2, j13 = -c23/2, j23 = c33.   [c_ij = coefficient of x_i x_j]
  i.e. jets of f/x1^L at e_1 in the coordinates u=x2/x1, v=x3/x1; the weight w = k+l-1 (principal grading, H_pr = E22+2E33) of j_kl is the power of a in theta_kl(a) = c_kl a^{k+l-1}
  (verified: sum_kl c_kl a^{k+l-1} j_kl = f(1,a,a^2), e.g. L=2: j03 + 3 j12 = c22 + c13 = [a^2] f(1,a,a^2)). So (a) of J1 is verified including the constants c_kl.
* **t != 0, L=1** (all 9 covectors, monomial coefficients c_ij = coeff of x_i y_j):
    <(01)(01)| = c11, <01|02| = -c12, <01|12| = c13, <02|01| = -c21, <12|01| = c31,
    <02|02| = c22 + t (c21 - c12),   <02|12| = -c23 + t (c13 - c31),   <12|02| = -c32 + t (c13 - c31),   <12|12| = c33 + t (c32 - c23).
  (rows (i,j) with jets i=(01),(02),(12) at site 1 and j at site 2.)  Only the 4 rows with both... couple the sites.
* **The operator A_1 is universal in L** (L = 1,2,3 checked by exact fit with L-independent coefficients, in the E_ij = x_i d_j basis with E_11,E_22,E_33 independent):
      A_1 = (E22+E33)^(1) E21^(2) - E21^(1) (E22+E33)^(2)  +  E23^(1) E31^(2) - E31^(1) E23^(2)  +  E33^(1) E32^(2) - E32^(1) E33^(2)
  (= (E22+E33) ^ E21 + E23 ^ E31 + E33 ^ E32 with X^Y := X^(1)Y^(2) - Y^(1)X^(2)); it changes the principal weight H_pr by +1 and the gl3 weight by -alpha_1 or -alpha_2 (simple roots),
  it commutes neither with the total gl3 nor with the total Cartan, NOT proportional to Omega_12 (and not in span{C_ij, C'_ij, single-site E}), but cross-site only.
  At the level of jets (M_1 := Y_1 Y_0^{-1}, rows/cols = jet pairs [kl|k'l']) the matrix has **small integer entries and is universal in L** (restriction of M_1(L=3) to labels <= 3 equals M_1(L=2), and
  of M_1(L=2) to labels <= 2 equals M_1(L=1); same for M_2 (L=2 vs 3)). So Y(t) Y_0^{-1} =: M(t) is the truncation of ONE infinite matrix on Lambda^2 U (x) Lambda^2 U, U = span{e_0,..}.
  Operator-Schmidt rank (site1 vs site2 factor) of M_1, M_2, M_3 = 6, 18, 40 (L-independent).
* **A_2 != A_1^2/2** (L=2: ratios A_2/A_1^2 take values 0,1/2,1/3,1/4,1/5 entrywise); M(t) is not exp(t M_1), and T^{-1}T' does not truncate.
  A_2 is a universal order-(2,2) differential operator (fit exists with L-independent coefficients, stacked L=2,3, but not unique with only two L's) - no pattern found.

(Correction of the L=1 list above: the nine rows are  <01|01>=c11, <01|02>=-c12, <01|12>=c13, <02|01>=-c21, <02|02>=c22+t(c21-c12), <02|12>=-c23+t(c13-c31),
<12|01>=c31, <12|02>=-c32+t(c13-c31), <12|12>=c33+t(c32-c23); c_ij = coefficient of x_i y_j.)
**Grading (confirmed after Q1's hint):** every entry of Y(t) is a pure monomial: t-power = (jet level of the row, sum of k+l-1 over both sites) - (principal weight alpha2+2alpha3 summed over the two kets of the column), L=1,2,3
exactly (pairs (t^n, level-weight) = (n,n) only). So Y(t) = t^rho Y(1) t^-kappa; all content is the fixed rational Y(1); Y(0) = y(x)y.

## §2 DERIVATION of Y (not a fit): Y rows = <x0| nu_i(site 1) nu_j(site 2) evaluated on the Gaudin Hamiltonians
Key facts (all exact):
* In the companion frame the operators H_s,a_s,b_s (Laurent coefficients of S2,S3 with twist g_c(zeta), `Operators(ch, gc(zeta))`) have the Psi_c as eigenvectors (residual 1e-77) with the Q-side Laurent data.
* For one-row Lambda=[L,0,0] (c2=0) the dual oper psi'''-S1 psi''+S2 psi'-S3 psi has exponents {0,1,L+2}; every solution has jet (order <= L+1) in span(f0,f1) (Frobenius series f0=1+0x+g2 x^2+..., f1=x+d2 x^2+...),
  so p_kl/p_01 = nu_kl := k! l! (g_k d_l - g_l d_k), a polynomial in (H_s,a_s,b_s) and the partner data (H_t,a_t,b_t), e_i(zeta), delta = theta_s - theta_t. The recursion is `frobenius_planes(L)`/`nu_formulas(L)`.
  The no-log condition (resonance 0 -> 1) is  **a_s = 0**  (for all one-row states; a_1,a_2 are the zero operators on shell).  Examples: L=1: nu = (1, H_s, b_s) for (01),(02),(12).
  L=2: nu_02=H_s/2, nu_12=b_s/2, nu_03 = -(-2H^2 d+3.. )/(4d) etc. (full list printed by t_nu1.py; stored in jet_basis_sym_nu*.pkl).
* Since <x0|f(H,b)|Psi_c> = f(eigenvalues)<x0|Psi_c> on a basis, **<i,j|^c = <x0| nu_i(H_1,b_1; H_2,b_2) nu_j(H_2,b_2; H_1,b_1)** as an operator identity. Verified EXACTLY (all rows, all entries, Fractions):
  L=1,2,3 at zeta = P, Pp, X1 AND at zeta = 0 (nilpotent companion twist g = E21-type, e_k = 0); the operator formula is polynomial in e_k, so Y is e-independent ( = zeta-independence explained).
  Y(t) then follows with ONLY operator matrices (no eigenstates, no Q's): `Y_from_operators`.
* t-polynomiality explained: H_1 = h_1 - t(Omega_12 - L^2) (+ const), delta^{-1} = t enters nu polynomially.

## §3 Theta derived without state data (J2)
Theta = ev(a,a') Y(t)^{-1}, ev = Omega_c[a] (x) Omega_c[a'] = evaluation at ((1,a,a^2),(1,a',a'^2)) in the monomial dictionary (H0).
* L=1: Y = Y0(1+tA_1), A_1^2 = 0, and exactly  E(a,a') A_1 = (a a'/(a-a')) E(a,a')(1-P)  (P = site exchange = Omega_12 for the fundamental), hence
  **Theta = theta (x) theta' + z (theta(x)theta' - theta'(x)theta), z = t a a'/(a'-a)**, i.e. (1+z) theta(x)theta' - z theta'(x)theta, theta_kl = c_kl a^{k+l-1}, c = (1,-1,1): derived symbolically and checked entrywise (9/9) against general_kernel_rep_100.
* L=2: Theta(Y from states) = rep_200 exactly (36/36 symbolic polynomials). L=3: exact at two random rational points (100/100), from state-derived Y and again from the OPERATOR-derived Y(t) (7 theta12 values, degree 3 fit + checks).
* **L=4 (225 entries): Theta from the operator formula only = general_kernel_rep_400 exactly** (two random points). So the master kernels for lambda_1 = 3,4 are DERIVED (the cached fits are confirmed), and the same code gives
  Y(t) for any L: [5,0,0] in progress (see below). `t_op4.py L` -> `jet_basis_sym_Yops<L>.pkl` (exact rows at 8 values of t), `poly_fit_in_t`, `theta_check`.
* **[5,0,0] (dim V = 21, 441x441)**: Y(t) from operators only, exact, degree 5 in t (6 fit points + 1 check; `jet_basis_sym_Ycoef_ops5.pkl`, list of sympy coefficient matrices Y_0..Y_5 in the GT companion basis).
  M_k = Y_k Y_0^-1 restricted to labels <= 5 equals the L=4 matrices for k=1..4 (universality confirmed again). Theta^{[5,0,0]} = ev Y^-1 is therefore available (not compared: no cache); M1 can test its fit against it.

## §4 [3,0,0] (J3)
Done by two routes (state data, 100x100, rational reconstruction with 0 failures; and operator formula): identical Y, t-degree 3, same universal M_k, Theta = rep_300 exactly. [4,0,0] state-derived Y points cached in
`jet_basis_sym_Ypts4.pkl` (8 t values, 0 reconstruction failures); operator-derived Y agrees with the kernel rep_400.

## Surprises / not done
* H0 holds with no corrections (even the GT scalars only a diagonal dictionary).  A_1 universal & explicit, but A_n (n>=2) found only as non-unique order-(n,n) fits; M(t) is not an exponential; no closed form of Y(1) or of the master kernel beyond n=1 was found.
* a_s = 0 as an operator identity for one-row reps (no-log), so the whole jet basis is generated by only (H_s, b_s): Y rows = <x0| nu_i nu_j, nu_kl = Frobenius Plucker ratios. Not done: b_1-eigencovector test (c), symbolic-in-L nu_kl, the closed form of Y(1).
