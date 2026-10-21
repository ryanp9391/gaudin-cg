# RESULTS: the Gaudin two-vacua Q-side formula as the formal h -> 0 limit of the XXX ten-product formula (agent D2, 2026-10-16)

Files: `twovac_limit.py` (`python -X utf8 twovac_limit.py`, about 1 minute, exact Fractions / sympy until the mpmath checks and L3; prints everything below), this file.
Reads (never edits) `xxx_chain.py`, `gaudin_chain.py`, `qside.py`, and the extract `Paul/Mathematica/Gaudin/ref/xxx-210-infinf-formula-cells.txt`.
Verdict: **case (i)**. The full RHS of the XXX identity has leading order **N = 10** (relative to the naive h^4), its coefficient involves **only the n = 0 jets** of the
Gaudin dual Q's (J_{0,0}, J_{0,1} = kappa J_{0,0}, J_{0,2}, J_{0,3}); the q_1 (n = 1) jets enter first at order 11. The formula is verified on all 64 Gaudin states to 1e-55 with
state-independent constant exactly 1 (the constant predicted by L1). Two corrections to the plan's heuristics are recorded in "Surprises".

Conventions used throughout: Qnum_G(B) := <Phi_inf|Psi_B>/<x0^G|Psi_B>, Phi_inf = Omega_D^inf[zeta,omega1] (x) Omega_D^inf[zeta,omega1'], r_k^inf as in the round-2 plan, x0^G = x^G(zeta) (x) x^G(zeta) with x^G
of L1(b); P^j = e^{zeta_j u} p_j (monic p_j of degree mu_j), jets J^j_{n,k} = (d/du)^k P^j_n(theta_s) of the h-expansion Q^j(u;h) = sum_n h^n P^j_n(u), P^j_n = e^{zeta_j u} q^j_n (so the factor e^{zeta_j m h} of
z_j^{m} is inside the jets and the Frobenius/Wronskian rules below are those of P^j). theta12 = theta1 - theta2.

## L0. Parse and verification

Parsed from cell 62 (`sympy.parsing.mathematica.parse_mathematica`, sS = s, tt = t): the 40 kernel laws kappa^{(cc')}_{(mm'),(nn')}(s,t) and the ten products {c, c', f, g, f', g'}
(blocks (0,0),(0,1),(1,0),(1,1) with 3,3,2,2 products).
* Conventions read off the setup cells: ladder points theta_s + m h, m = -3,-2,-1,0; minor_s(m,m') = Q2(m')Q3(m) - Q3(m')Q2(m) for m < m'; Bk2_s[F,G] = <Q3 F><Q2 G> - <Q2 F><Q3 G> = sum_{m<m'} (F^G)_{mm'} minor_s(m,m'),
  (F^G)_{mm'} = F_m G_m' - F_m' G_m; F_m = a^m f_m, G_m = a^{m+e} g_m, so (F^G)_{mm'} = a^{m+m'+e}(f^g)_{mm'}; M_vac = minor_s(-3,-2) Q3(theta_s - h) (cell 49: `mVac`); (e,e') = (7,5),(7,4),(6,5),(6,4).
* **sum_k (f_k^g_k)_{mm'} (f'_k^g'_k)_{nn'} reproduces the printed kernel law identically in (s,t) for all 144 block entries (4 blocks x 6x6 slots): 0 mismatches** (`sp.cancel` of the difference). The slot (-1,0) never occurs.
* Kernel laws are polynomials in s (degree <= 4) with t-denominators (t-1) t (t+1) only. Individual wedge entries are Laurent in h (orders -2..+2: 1/(s-1) = O(1/h) in products 2 and 7, t = theta12/h in product 2), the wedge sums are O(1).
* **L0b (extra, strong): the whole identity checked exactly (rational arithmetic, finite h = 1/2) on the top state e1 (x) e1** (weights (4,2,0), t_a = vacuum eigenvalues, q^3 = 1, q^2 = x^2 - 29x/35 + 144/245) at three generic (a,a'):
  <Phi_c|Psi_c>/<x0_c|Psi_c> M1 M2 == sum_{cc'} Q3 Q3 sum_k Bk2 Bk2 as exact rationals, with Psi_c = (P^-1 (x) P^-1) psi_d, Phi_c = gaussCov{a,0,a} (x) gaussCov{a',0,a'}. This fixes all conventions above (including that M_vac uses Q3(theta_s - h)).
  Also exact: kappa(z,w=(1/a,inf)) (L1a) * <Omega_c(a)|P^-1 e1> = 1, and minor_s(-1,0) = 0.

## L1. h-dependent ingredients

**(a)** kappaTie(z,w) = Delta_1(nP3)^(lam1-lam2) Delta_2(nP3)^(lam2-lam3) Delta_3^lam3 with Lambda = (2,1,0), n = nOf at w2 = infinity, P3 = Pphys3 (cell 252). Exact (sympy):
Delta_1 = w1^2/((z2-w1)(z3-w1)), Delta_2 = w1/((z3-w1)(z1-z2)), so

    kappa(z,w) = w1^3 / ((z2-w1)(z3-w1)^2 (z1-z2)),   z = e^{h zeta}, w1 = e^{h omega}:   kappa = h^-4 / ((zeta2-omega)(zeta3-omega)^2 (zeta1-zeta2)) (1 + O(h)).

So kappa kappa' = h^-8 / c, **c(zeta,omega,omega') = (zeta2-omega)(zeta3-omega)^2 (zeta2-omega')(zeta3-omega')^2 (zeta1-zeta2)^2** (a pure product of ladder-type differences; no omega-omega' coupling).

**(b)** x0(z) = e1 P_phys(z)^-1 (one site, notebook order = descending GT patterns, index 0 = highest weight) has components exactly linear in z (all orders in h):
(z2-z1, z3-z1, z2-z1, 2z3-z1-z2, 3(z2-z1), 3(z3-z1), 2(z3-z2), 3(z3-z2)). Hence x0(e^{h zeta}) = h x^G(zeta) + O(h^2) with **c0 = 1 exactly** and x^G(zeta) = the same linear forms in zeta
(P2's "parallel" statement made exact: x^G is what the Gauss-factor construction gives at z = zeta), and **<x0^d|psi_d> = h^2 <x0^G|psi> (1 + O(h))**. (In gl3_reps order the components are permuted by perm = [7,5,6,4,2,1,3,0]: gl3 index perm[a] <- notebook index a.)
Also Phi_d -> Phi_inf: r_k(z,(e^{h omega},inf)) = r_k^inf(zeta,omega) + O(h) (ratios of O(h) quantities), so <Phi_d|psi_d> = <Phi_inf|psi> (1 + O(h)).

**h-order table** (exponent e: factor ~ h^e):

| factor | order | remark |
|---|---|---|
| kappa(z,w), kappa' | h^-4 each, product h^-8 | closed form above |
| <x0^d|psi_d> | h^2 | x ~ h per site |
| <Phi_d|psi_d> | h^0 | |
| Qnum_XXX = <Phi_d|psi_d>/(kappa kappa' <x0^d|psi_d>) | **h^6** | h^8/h^2; limit h^6 c Qnum_G |
| a^m, a'^n | 1 + O(h) | a = e^{-h omega} |
| s - 1, 1/t | O(h), O(h) | s = e^{-h(omega'-omega)}, t = theta12/h |
| kernel entries kappa^{cc'}(s,t) | h^0 | polynomial in s, t-degree 0; singular values of the 10x10 K(h): 1, h^2, h^2, h^4, h^4, h^4 (rank 6); K(0) = U (x) U, U = (u_0|u_1), u_c = (c, 1-3c, -2+5c) on slots [(-3,-2),(-3,-1),(-2,-1)] |
| wedge entries f,g,f',g' | h^-2 .. h^2 | cancel inside the kernel laws |
| Q^j(theta_s + m h), Q3(theta_s - c h) | h^0 | jets J_{n,k} h^{n+k} m^k/k! |
| minor_s(m,m') | **h^2** | (the h^1 term is -(m'-m) W(theta_s) = 0) |
| M_vac = minor_s(-3,-2) Q3(theta_s - h) | h^2 per site, **= 2 h^2 W'(theta_s) P3(theta_s)** | not (5/2) h^2 W' (see Surprises) |
| naive order of RHS | h^4 | kernel O(1) x minor h^2 x minor h^2 |
| exact leading order of RHS = Qnum_XXX M1 M2 | **h^10** | cancellations through h^9 |

Independent check of (a),(b) with no Gaudin input (L1c): top state, exact t-eigenvalues at h = 1e-3..1e-6 (60 digits), F(h) = <Phi_c|Psi_c>/<x0_c|Psi_c> ~ h^6 x 2.4653983408...(Richardson) vs the prediction c Qnum_G(top) = (zeta2-w)(zeta3-w)^2 (zeta2-w')(zeta3-w')^2 = 102245/41472 = 2.4653983410 (zeta = (3/2,5/4,2/3), omega = 5/2, omega' = 7/4; Qnum_G(top) = 1/(zeta2-zeta1)^2).

## L2. The formal limit

**(a) The exact ladder constraint (the generalisation of the Wronskian zero W(theta_s) = 0).** In the dual Baxter equation nu3(x)nu3(x-h)nu3(x-2h)Q(x+h) - t1 nu3(x-h)nu3(x-2h)Q(x) + t2 nu3(x-2h)Q(x-h) - t3 Q(x-2h) = 0, at x = theta_s one has nu3(theta_s) = 0 AND
t3(theta_s) = chi3 nu1 nu2(x-h) nu3(x-2h)|_{x=theta_s} = 0 (nu2(x-h) = prod_s (x - theta_s)). So Q^j(theta_s - h) = rho_s Q^j(theta_s), rho_s = t1(theta_s) nu3(theta_s-h)/t2(theta_s), **the same rho_s for every j**, hence
**minor_s(-1,0) = 0 exactly, to all orders in h, for every XXX eigenstate** (checked exactly on the top state). Consequences: the 4 ladder values of (Q2,Q3) lie in a 2-plane; with Q(-3) = alpha Q(0) + beta Q(-2), all minors are D times
(-3,-2): -alpha, (-3,-1): rho beta, (-2,-1): rho, (-3,0): beta, (-2,0): 1, (-1,0): 0 (D = minor(-2,0)); alpha, beta, rho are explicit rational functions of t1,t2 at theta_s, theta_s - h (verified against the actual Q's). Leading values (alpha,beta,rho) = (-2,3,1).
In jets: with Q(-1) = rho Q(0) imposed, the generic weight-2 minor (m-m')[(m+m')p/2 + Q] (p = W'(theta_s), Q = J2_00 J3_11 - J3_00 J2_11 + J2_10 J3_01 - J3_10 J2_01, an n=1 combination) becomes (m-m')(m+m'+1) p/2.

**(b) Mechanism of the cancellations.** K(h) = U (x) U + O(h^2) with the site contraction Z_s = sum_c Q3(-c) u_c . (minor(-3,-2), minor(-3,-1), minor(-2,-1)) = 0 identically, in all jets and in rho_i, **through weight 4**
(symbolic, all J^j_{n,k} with n+k <= 6 and rho_1..rho_6 free; only minor(-1,0) = 0 imposed); first nonzero at weight 5.

**(c) Leading order, exact (Fractions, 16-term series) with generic jets.** Random rational jets J^j_{n,k} (all n,k), random rho_r, with only the ladder constraint imposed: the RHS series has **valuation exactly 10** (coefficients of h^4..h^9 vanish
identically; three random trials, minor(-3,-2) has valuation 2). Perturbing one jet at a time: the coefficient of h^10 depends only on J^j_{0,0}, J^j_{0,2}, J^j_{0,3} and rho_1 = -kappa (J^j_{0,1} = kappa J^j_{0,0}). The first order at which each other variable enters:
J_{n,k} (k != 1, n >= 1): 10 + n (k <= 3), 10 + n + (k-3) for k >= 4; rho_r: 9 + r (r >= 1); J_{0,k}, k >= 4: 7 + k. (J_{n,1}, n >= 1, are fixed by the ladder relation: J_{1,1} = J_{0,2}/2 - rho_2 J_{0,0} - rho_1 J_{1,0}, ...) So **no q_1 (n=1) or q_2 jet enters at order 10; they enter at order 11 and 12.**
This is invariant: the h^10 coefficient is a function of the q_0 jets alone for arbitrary values of all other jets.

**(d) The formula.** Effective local ladder function (all other jets can be set to zero at order 10): Q^j(theta_s + m h) = J^j_0 (1 + kappa m h) + (J^j_2/2) h^2 m(m+1) + (J^j_3/6) h^3 m(m^2-1). With kappa_s := P3'/P3 = P2'/P2 (Frobenius; W(theta_s) = 0),
p_kl := J2_k J3_l - J3_k J2_l (p_01 = 0, p_12 = kappa p_02, p_13 = kappa p_03), P3^(k) := J3_k, W = P2 P3' - P3 P2' (W' = p02, W'' = p12 + p03), the coefficient of h^10 of the full RHS (symbolic in omega1, omega1', 1/theta12) is

    RHS_10 = sum_{i,j=1}^{5} G_ij(omega1,omega1',theta12) N_i(theta1) N_j(theta2),    N = (P3 p02, P3' p02, P3 p03, P3' p03, P3 p23)   [Pluecker basis {P3^(c) p_kl}: 5 of the 12 elements; P3' p02 = kappa P3 p02]

and, with LHS = Qnum_XXX M1 M2, M_vac = 2 h^2 (P3 p02)(1 + O(h)) per site, kappa kappa' = h^-8/c, <x0^d|psi_d> = h^2 <x0^G|psi>:

    **4 c Qnum_G(B) (P3 p02)(theta1) (P3 p02)(theta2) = sum_ij G_ij N_i(theta1) N_j(theta2)**,     equivalently     **Qnum_G(B) = (1/(4c)) sum_ij G_ij nu_i(theta1) nu_j(theta2)**,  nu = (1, kappa, pi, kappa pi, eps),  pi = p03/p02,  eps = p23/p02 = (P3''/P3) pi - P3'''/P3

(so Qnum_G depends on each site only through the four local numbers kappa_s, pi_s, P3''/P3, P3'''/P3; the b2 p03 and b3 p02 columns of the raw 6x6 form are exactly opposite, and P3 p23 = P3'' p03 - P3''' p02 is the Pluecker identity).
Properties (exact): G has rank 5; site-exchange symmetry G_ji(omega1,omega1',theta12) = G_ij(omega1',omega1,-theta12); homogeneity deg 10 - w_i - w_j in (omega, 1/theta12) with weights w = (2,3,3,4,5).
**zeta does not appear in the kernel**; it enters only through the state-independent normalisation c(zeta,omega,omega') (the frame-tie scalar kappa kappa').

Kernel, untwisted jets (w1 = omega1, w1p = omega1', thi = 1/theta12; row/column order N = (P3 p02, P3' p02, P3 p03, P3' p03, P3 p23)):

```
G[1,1] = -4 w1^2 w1p^2 (3 thi^2 + thi w1 - thi w1p - w1 w1p)
G[1,2] = 2 (6 thi^3 w1^2 + 24 thi^3 w1 w1p + 6 thi^3 w1p^2 + 2 thi^2 w1^3 + 12 thi^2 w1^2 w1p - 2 thi^2 w1p^3 + 2 thi w1^3 w1p - 6 thi w1^2 w1p^2 - 2 thi w1 w1p^3 - 3 w1^3 w1p^2)
G[1,3] = -2 w1^2 w1p^2 (2 thi + w1)
G[1,4] = 4 (-3 thi^4 - 4 thi^3 w1 - 2 thi^3 w1p - thi^2 w1^2 + 3 thi^2 w1 w1p + thi^2 w1p^2 + 3 thi w1^2 w1p + thi w1 w1p^2 + w1^3 w1p)
G[1,5] = -2 (6 thi^3 + 8 thi^2 w1 + 4 thi^2 w1p + 4 thi w1^2 + 2 thi w1 w1p + w1^3)
G[2,1] = -2 (6 thi^3 w1^2 + 24 thi^3 w1 w1p + 6 thi^3 w1p^2 + 2 thi^2 w1^3 - 12 thi^2 w1 w1p^2 - 2 thi^2 w1p^3 - 2 thi w1^3 w1p - 6 thi w1^2 w1p^2 + 2 thi w1 w1p^3 + 3 w1^2 w1p^3)
G[2,2] = -12 thi^2 w1^2 - 48 thi^2 w1 w1p - 12 thi^2 w1p^2 - 4 thi w1^3 - 12 thi w1^2 w1p + 12 thi w1 w1p^2 + 4 thi w1p^3 + 9 w1^2 w1p^2
G[2,3] = 12 thi^4 + 16 thi^3 w1 + 8 thi^3 w1p + 4 thi^2 w1^2 - 12 thi^2 w1 w1p - 4 thi^2 w1p^2 - 4 thi w1^2 w1p + 4 thi w1 w1p^2 + 3 w1^2 w1p^2
G[2,4] = 2 (6 thi^3 + 8 thi^2 w1 + 4 thi^2 w1p + 2 thi w1^2 - 6 thi w1 w1p - 2 thi w1p^2 - 3 w1^2 w1p)
G[2,5] = 6 thi^2 + 8 thi w1 + 4 thi w1p + 3 w1^2
G[3,1] = 2 w1^2 w1p^2 (2 thi - w1p)
G[3,2] = 12 thi^4 - 8 thi^3 w1 - 16 thi^3 w1p - 4 thi^2 w1^2 - 12 thi^2 w1 w1p + 4 thi^2 w1p^2 - 4 thi w1^2 w1p + 4 thi w1 w1p^2 + 3 w1^2 w1p^2
G[3,3] = w1^2 w1p^2
G[3,4] = -2 (-6 thi^3 - 2 thi^2 w1 + 2 thi^2 w1p + 2 thi w1 w1p + w1^2 w1p)
G[3,5] = 6 thi^2 + 4 thi w1 + w1^2
G[4,1] = 4 (-3 thi^4 + 2 thi^3 w1 + 4 thi^3 w1p + thi^2 w1^2 + 3 thi^2 w1 w1p - thi^2 w1p^2 - thi w1^2 w1p - 3 thi w1 w1p^2 + w1 w1p^3)
G[4,2] = 2 (-6 thi^3 + 4 thi^2 w1 + 8 thi^2 w1p + 2 thi w1^2 + 6 thi w1 w1p - 2 thi w1p^2 - 3 w1 w1p^2)
G[4,3] = -2 (6 thi^3 + 2 thi^2 w1 - 2 thi^2 w1p - 2 thi w1 w1p + w1 w1p^2)
G[4,4] = -4 (3 thi^2 + thi w1 - thi w1p - w1 w1p)
G[4,5] = -2 (2 thi + w1)
G[5,1] = -2 (-6 thi^3 + 4 thi^2 w1 + 8 thi^2 w1p - 2 thi w1 w1p - 4 thi w1p^2 + w1p^3)
G[5,2] = 6 thi^2 - 4 thi w1 - 8 thi w1p + 3 w1p^2
G[5,3] = 6 thi^2 - 4 thi w1p + w1p^2
G[5,4] = 2 (2 thi - w1p)
G[5,5] = 1
```

**Twisted-jet form** (compact; "the two vacua couple only through delta = omega' - omega", the Gaudin image of s = a'/a). Use jets of e^{-omega_s u} P^j at theta_s (omega_s = omega1 for site 1, omega1' for site 2): J~_k = (d/du - omega_s)^k P at theta_s,
kappa~ = kappa - omega_s, p~_kl, N~ = (P~3 p~02, P~3' p~02, P~3 p~03, P~3' p~03, P~3 p~23) (p~02 = p02, p~03 = p03 - 3 omega p02). The same RHS_10 = sum G~_ij N~_i N~_j with G~ a function of delta and thi only (verified symbolically and on the 64 states):

```
G~[1,1] = -4 d^2 thi^2 (d + 3 thi)(2 d + 3 thi)      G~[1,2] = -4 d thi^2 (2 d + 3 thi)^2             G~[1,3] = -4 d thi^2 (d^2 + 5 d thi + 3 thi^2)
G~[1,4] = -4 thi^2 (3 d^2 + 8 d thi + 3 thi^2)        G~[1,5] = -4 thi^2 (2 d + 3 thi)
G~[2,1] =  4 d thi^2 (2 d + 3 thi)^2                  G~[2,2] =  2 d thi (d + 6 thi)(2 d + 3 thi)      G~[2,3] =  2 thi^2 (5 d^2 + 10 d thi + 6 thi^2)
G~[2,4] =  4 thi (d^2 + 5 d thi + 3 thi^2)            G~[2,5] =  2 thi (2 d + 3 thi)
G~[3,1] =  4 d thi^2 (d^2 + 5 d thi + 3 thi^2)        G~[3,2] =  2 thi^2 (5 d^2 + 10 d thi + 6 thi^2)  G~[3,3] =  2 d thi^2 (d + 6 thi)
G~[3,4] =  4 thi^2 (2 d + 3 thi)                      G~[3,5] =  6 thi^2
G~[4,1] = -4 thi^2 (3 d^2 + 8 d thi + 3 thi^2)        G~[4,2] = -4 thi (d^2 + 5 d thi + 3 thi^2)       G~[4,3] = -4 thi^2 (2 d + 3 thi)
G~[4,4] = -4 thi (d + 3 thi)                          G~[4,5] = -4 thi
G~[5,1] =  4 thi^2 (2 d + 3 thi)                      G~[5,2] =  2 thi (2 d + 3 thi)                   G~[5,3] =  6 thi^2
G~[5,4] =  4 thi                                      G~[5,5] =  1
```
(d = omega1' - omega1.)

**Direct Q_1 basis.** Dual Wronskian (RESULTS-Qside section 2): W = P2 P3' - P3 P2' = (zeta3-zeta2) e^{sigma u} Pi(u) Q1(u), Q1 = e^{-zeta1 u} q_1 (monic direct Q), sigma = zeta1+zeta2+zeta3, Pi = (u-theta1)(u-theta2).
With W(theta_s) = 0: p02 = W'(theta_s) = A_s Q1(theta_s), p03 = W'' - kappa W' = A_s [(Lambda_s - kappa) Q1 + 2 Q1'](theta_s), A_s = (zeta3-zeta2) e^{sigma theta_s} Pi'(theta_s), Lambda_s = 2 sigma + 2/Pi'(theta_s) (Pi'(theta_1) = theta12, Pi'(theta_2) = -theta12), and p23 via
P3 p23 = P3'' p03 - P3''' p02. Hence nu = (1, kappa, pi, kappa pi, eps) with pi = (Lambda_s - kappa) + 2 Q1'/Q1, eps = (P3''/P3) pi - P3'''/P3, i.e. the formula in the basis {P3^(c)(theta_s) Q1^(k)(theta_s)} is

    Qnum_G(B) = (1/(4c)) sum_ij G_ij nu_i(theta1) nu_j(theta2),   nu_s = (1, kappa_s, pi_s, kappa_s pi_s, eps_s),   pi_s = Lambda_s - kappa_s + 2 Q1'(theta_s)/Q1(theta_s).

Expanded, each site uses the 9 monomials {P3 Q1, kappa P3 Q1, kappa^2 P3 Q1, P3 Q1', kappa P3 Q1', P3'' Q1, kappa P3'' Q1, P3'' Q1', P3''' Q1} (divided by P3 Q1): k = 0,1 in Q1^(k), c <= 3 in P3^(c). kappa_s (= P3'/P3 at theta_s, equivalently the Frobenius constant) is left as a symbol, as in the spec.
Numerically verified (below): W'(theta_s) = A_s Q1(theta_s) to 1e-55 and the Q1-form reproduces Qnum on all 64 states.

**(e) The q_1 equation is not needed** (case (i)); no q_1 jet enters at the leading order. (For the next layer, order 11: it receives J_{1,0}, J_{1,2}, J_{1,3} and rho_2, which is where the order-h^4 Baxter equation for q_1 would be needed; not done.)

## L3. Numerical test (60 digits, 64 Gaudin states, Lambda = (2,1,0), L = 2)

Qnum_G(B) = <Phi_inf|Psi_B>/<x0^G|Psi_B> from the P1 eigenstates (all 64 nonzero, |Qnum_G| in [0.25, 337] at the base point), P^j jets at theta_s from P1's dual polynomials, direct Q1 from P1's direct polynomials.
Tested quantity: R10/(c Qnum_G 4 (P3 p02)_1 (P3 p02)_2), R10 = sum G_ij N_i N_j, over all 64 states; also Qformula/Qnum with the Q1-basis formula, and the twisted-jet kernel against the untwisted.

| point | (omega1, omega1') | ratio (all 64 states) | max rel. deviation | Q1-basis form | twisted kernel |
|---|---|---|---|---|---|
| zeta = (3/2,5/4,2/3), theta = (2/7,-3/5) | (5/2, 7/4) | **1.0** | 2.8e-56 | 1.0 (1.3e-55) | 9e-56 |
| same | (-2/5, 3) | 1.0 | 2.2e-56 | 1.0 (2.2e-56) | 1.5e-56 |
| same | (9/7, -5/3) | 1.0 | 9.2e-56 | 1.0 (4.6e-56) | 1.5e-54 |
| zeta same, theta = (1/5,-4/9) | (5/2, 7/4) | 1.0 | 4.2e-56 | 1.0 (1.5e-55) | 1.1e-55 |
| zeta = (7/5,-3/4,1/6), theta = (2/7,-3/5) | (5/2, 7/4) | 1.0 | 3.4e-55 | 1.0 (3.1e-56) | 2.2e-54 |

* The state-independent constant between LHS and formula is **exactly 1 = the constant predicted by L1** (c from kappa kappa', x^G with c0 = 1): no leftover normalisation.
* Gaudin facts used, measured on all 64 states: Frobenius |P2'/P2 - P3'/P3|(theta_s) <= 2.3e-57, Wronskian |W(theta_s)| <= 1.6e-55.
* Kernel uniqueness: in the basis N (5 monomials per site, 25 products) the 64 x 25 matrix of product monomials has full numerical rank 25 (normalised singular values from 1 down to 1e-8 at 60 digits; per site the 64 x 5 matrices have rank 5), so on the 64 states G is **the unique** kernel in this basis (no null space).
  Embedded in the 12 x 12 basis {P3^(c) p_kl} (c = 0,1, kl in {02,03,12,13,23,...}) with p_12 = kappa p_02, p_13 = kappa p_03 the representatives differ by the null space of that larger basis.
* Comparison with D1's fitted kernel: D1's `RESULTS-twovac-fit.md` was not available when this file was written (see the note at the end).

## L4. Naive companion-frame limit (requested by W2's update)

lim_{t->inf} kappa^{(cc')}(s,t): 25 of the 40 laws survive (the rest vanish), all equal to **s * u_c[mm'] u_c'[nn']** with u_0 = (0,1,-2), u_1 = (1,-2,3) on slots (-3,-2),(-3,-1),(-2,-1) (the same rank-one U (x) U as K(0); the laws with slots (-3,0),(-2,0) vanish).
With a = 1/omega, a' = 1/omega' fixed and minors -> h^2 (m-m')(m+m'+1) p/2, the naive RHS/(M1_vac M2_vac) = (1-a)^3 (1-a')^3: **state-independent**. W2's kappa^G (omega2 -> inf) = omega^3/((zeta2-omega)(zeta3-omega)^2(zeta1-zeta2)) is the XXX kappa with z -> zeta, w -> omega,
so the Gaudin companion Q-number is c Qnum_G/(omega^3 omega'^3), which varies with the state (L3). **The two do not match**: the fixed-a limit is not the Gaudin scaling limit (there a = e^{-h omega} -> 1); it reproduces only the power h^6 ((1-a)^3(1-a')^3 ~ h^6 omega^3 omega'^3),
while all state dependence comes from the cancellation hierarchy of L2 (K(h) singular values 1, h^2, h^2, h^4, h^4, h^4 and the ladder constraint), invisible to a rank-one naive kernel. W2's reference tables (notebook cells 51-52) were not used; L3 uses P1's states and the same x^G.

## Surprises

1. The cancellations of the RHS from the naive h^4 to h^10 are not between independent pieces: they are all consequences of one exact local relation, minor_s(-1,0) = 0 (Baxter at x = theta_s with t3(theta_s) = 0). Imposing only this relation, with ALL other jets (q_1, q_2, ... jets, rho_r) random, gives valuation exactly 10.
2. **M_vac = 2 h^2 W'(theta_s) P3(theta_s), not (5/2) h^2 W'** (plan section 1): the weight-2 minor contains the n=1 combination Q (q_1 jets), which the exact relation fixes to p/2. The formal "Taylor of the n=0 jets only" expansion is wrong at the leading order of the minors, although the final formula only involves n=0 jets.
3. The kernel in twisted jets depends on omega, omega' only through delta = omega' - omega (the Gaudin image of "the vacua couple only through s = a'/a"); in untwisted jets the omega-dependence is the binomial twist of this kernel.
4. Only 5 of the 12 Pluecker-basis elements {P3^(c) p_kl} occur (c = 0 with p02, p03, p23; c = 1 with p02, p03); Qnum_G depends on a site through 4 local numbers (kappa, p03/p02, P3''/P3, P3'''/P3).
5. x0 is exactly linear in z (x^G exact, c0 = 1); kappa kappa' = h^-8/c exactly with no omega-omega' coupling; the Gaudin limit of LHS/RHS fixes N = 10 independently of the jet analysis (L1c: top state to 1e-10).
6. The plan's expectation "(i) N = 0" in the spec refers to the O(h^0) normalisation ratio <Phi_d|psi_d>; the full RHS (including M1 M2 ~ h^4 and kappa kappa'/<x0|psi> ~ h^-6) has leading order N = 10.

## Not done / limits

* The order-11 layer (q_1 corrections) and the q_1 ODE from the order-h^4 Baxter expansion were not derived: not needed since q_1 does not enter the leading coefficient.
* Only Lambda = (2,1,0), L = 2 (the XXX closed forms exist only for [2,1,0]); the effective-model symbolic computation is exact in (omega1, omega1', theta12); the 64-state tests are at 5 parameter points.
* No Wolfram, no notebook, no edits outside the two files. The reduction "n >= 1 jets do not enter at order 10" was established by exact random-jet series (all jets perturbed one at a time) and symbolically in the effective model, not by a closed-form proof for general states.
