# RESULTS-master-kernel-2 (agent M3, round 5, Python only) -- COMPLETE for K5/K6/K7 as scoped (see Not done)

Files: `master_kernel2.py` (library: `theta_T(L)` exact polynomial kernel from M2's Y, `symbol_orders_exact(L)` symbol Q_n in E1,E2,H1,H2,a,ap),
`master_kernel2_Y1inv{3,4,5,6}.pkl` (exact Y(1)^-1, kappa, rho), `master_kernel2_rep_500.pkl` (exact [5,0,0] kernel, general_kernel T-format),
`mk2_*.py` drivers. Environment: python -X utf8, sympy 1.14.

## 0 Conventions
* T-format as in `general_kernel`: T[(i,j)] = {(n,p,q): c} meaning Theta_{ij} = sum c t^n a^p a'^q, (i,j) indices into `pairs_of(lam)^2`.
* Symbol: Q(xi;eta) with sum Theta_pq (x^k..)=(x2-x1)(y2-y1)Q, i.e. Q = sum Theta_{(k,l),(k',l')} e2^k h_{l-k-1}(xi) h2^k' h_{l'-k'-1}(eta)
  (h_m = complete homogeneous in x1,x2); variables E1,E2 (=e1,e2 of xi), H1,H2 (of eta). Pi=(1-aE1+a^2E2)(1-a'H1+a'^2H2), K=aa'/(a'-a)(Pi - F_a'G_a) (M1).

## 1 K5: exact [5,0,0] kernel from the operator route
**Key structural fact (new, exact, checked L=3,4,5):** M2's Y(t) is graded: Y(t)_{rc} = t^{rho_r - kappa_c} Y(1)_{rc}, with rho_r = (k+l-1)+(k'+l'-1) the jet level of the
row r=((k,l),(k',l')) and kappa_c = alpha_2+2alpha_3 summed over the two GT kets (principal weight). (0 grading conflicts in all Y_0..Y_5 entries.) Hence
 Y(t)^-1 = t^kappa Y(1)^-1 t^-rho  and  **Theta_r(a,a',t) = sum_c ev_c(a,a') t^{kappa_c-rho_r} (Y(1)^-1)_{cr}** with ev_c = s_p s_q a^{w_p} a'^{w_q} (monomial in a,a').
So ONE exact rational inverse of the 441x441 integer matrix Y(1)=sum_k Y_k (15414 nonzeros, entries <= 3e10) gives Theta as exact polynomials in (a,a',t) -- no per-point solves,
no grid, no reconstruction. Sparse DomainMatrix(QQ).inv(): 32 s (L=5), 4 s (L=4).
* Check: the same code reproduces the cached state-fit kernels rep_300 and rep_400 EXACTLY (100/100 and 225/225 entries, every monomial coefficient).
* [5,0,0]: 441 entries, 15875 monomials total, t-degree max 16 = 4(l1-1), a-degree max 10 = 2 l1: both agree with M1's support rule. Stored `master_kernel2_rep_500.pkl`.
* **Acceptance test (K5a)**, exact Theta^{[5,0,0]} vs the 441-state jet LHS `c*Qnum*M1*M2` (M1's `state_data`/`residuals`; w1=5/2, w1'=7/4, a=1/w1, a'=1/w1'; NO fit, NO truncation):

  | point | digits | max rel. resid | median | note |
  |---|---|---|---|---|
  | P, theta12=100 (t=.01) | **80** | **6.0e-56** | **5.4e-70** | 0 states > 1e-40 (the clean test: ~1e-70 = working precision) |
  | P, theta12=100 | 40 | 1.3e-15 | 7.1e-30 | worst states are numerical, see below |
  | P, theta12=1000 | 40 | 8.0e-7 | 1.7e-24 | |
  | P, theta12=31/35 | 40 | 1.9e-20 | 1.8e-25 | |
  | Pp, theta12=31/35 | 40 | 1.8e-15 | 3.2e-24 | |

  The 40-digit maxima are conditioning noise of the eigenstates (error grows like theta12^9 and the same states are accurate at 80 digits: 1.3e-15 -> 6e-56). At 40 digits M1's order<=2 truncation had 1.3e-2 / 2e-5 at the first two points; the exact kernel gains 13 / 18 orders of magnitude. So WP-A's "predict [5,0,0] on >= 200 states" is closed by the operator route: all 441 states agree to working precision. (Scripts `mk2_accept.py`, `mk2_accept80.py`; the 80-digit point was built with `master_kernel_bld.py 500 P 80 100`, 72 s.)
* [6,0,0] (dim 28, 784x784): Y(1) from `Y_from_operators` at TWO values of t (t=31/35, 6/5; 135 s each after a 153 s Frobenius step) via the grading (only Y(1) is needed: `mk2_Y6.py`); the two Y(1) agree entry by entry (46543 nonzeros), which checks grading + operator route at L=6. 784x784 inverse 203 s. Stored `master_kernel2_rep_600.pkl` (784 entries, 40931 monomials, t-degree max 20 = 4(l1-1)). No state data exists at L=6, so L=6 is validated through the universality below and the cumulant/operator checks of section 2 (all pass at lambda=6).

### K5 extra: universality (new, exact; this is the real structure)
Write Y(t)=Y_0 T(t), M(t):=Y(t)Y_0^{-1}, N(t):=M(t)^{-1} on the jet-pair space (rows/cols = pairs of one-row jets (k<l), two sites). Y_0 = entries of Y(1) with rho_r = kappa_c.
1. M(t)_{rp} = t^{rho_r-rho_p} M(1)_{rp}, rho = jet level k+l-1 summed over the two sites (a consequence of the grading).
2. **M(1) and N(1)=M(1)^{-1} are universal in l1**: M^{(L)}(1) restricted to labels <= L'+1 equals M^{(L')}(1), AND the same for the inverse (L'=3,4 vs L=4,5; and `N^{(6)}|_5 == N^{(5)}` exactly, 784x784 -> 441x441). So there is ONE infinite matrix N on (Lambda^2 U)^{(x)2}.
3. **All master kernels at once.** With theta_kl(a)=c_kl(l1) a^{k+l-1}, c_kl = (-1)^{k+l-1}(l-k) C(l1,l-1) C(l1+1,k)/l (`ckl_closed`):
   Theta^{[l1,0,0]}_{pq}(a,a',t) = sum_{r,s} c_r c_s a^{rho_r} a'^{rho_s} t^{rho_r+rho_s-rho_p-rho_q} N_{rs;pq},
   equivalently R_{pq}(x,y) = sum_{rs}(c_r c_s/(c_p c_q)) x^{rho_r-rho_p} y^{rho_s-rho_q} N_{rs;pq} (x=at, y=a't).
   **Checked: N^{(5)} plus this formula reproduces the exact kernels for l1=2,3,4,5 (36/100/225/441 entries, every monomial coefficient exact)** (`mk2_univTheta.py`). The whole l1-dependence is the explicit binomial c_kl(l1), which also truncates labels automatically (C(l1,l-1)=0 for l>l1+1). The t-grade of N is n=rho_r+rho_s-rho_p-rho_q >= 0 (order in t); N_n lowers the level by n.

## 2 K6: orders 3, 4, 5 closed; the master kernel as a normal-ordered exponential (all lambda)

### 2.1 What was fitted/verified (exact data lambda=1..6, assuming only polynomiality in lambda)
Data: exact kernels lambda=1..5 (state fit + operator), lambda=6 (operator). Symbol Q_n in (E1,E2,H1,H2,a,a'): `mk2_sym.py` (symbolic a,a'), `mk2_cumfit.py` (numeric a=1/3, a'=1/5, lambda<=6, n<=5).
* **Divisibility** Q_n = Pi^{l-n} P_n holds EXACTLY for l>=n: (l,n)=(3,3),(4,3),(4,4),(5,3),(5,4),(5,5) symbolic in a,a'; all n<=5 for l<=6 at a=1/3,a'=1/5. P_n is polynomial in lambda.
* **Order 3 in closed form (all l):** P_3 = l^2(l-1)(c_0+c_1 l+c_2 l^2+c_3 l^3) with **c_3 = K^3/6 exactly** (degree 6 in lambda as M1 expected; leading coefficient K^3/6). Fitted from l=2,3,4,5 with symbolic a,a' (c_0..c_3 have 240/246/194/108 terms in E,H,a,a'; stored `master_kernel2_P3coef.pkl`); verified at l=6 (cumulant test below) and by the operator formula of 2.2 at l=1..6.
* **Cumulant structure (key to the pattern).** With W_n := Pi^n [t^n] log( sum_m Q_m t^m / Pi^l ): W_n is a polynomial in lambda of degree exactly n+1 with no lambda^0, lambda^1 term. Verified by fitting lambda^2..lambda^{n+1} from the lowest lambdas and checking on the rest: n=1 (5 checks), n=2 (4), n=3 (3), n=4 (2), n=5 (1 check, at lambda=6) -- all exact. In the Newton basis
  W_n = l^2 * sum_{j=0}^{n-1} (l-1)_j Gamma_{n,j},  (l-1)_j = (l-1)(l-2)...(l-j),
  with Gamma_{n,0} = (-1)^{n-1} K^n/n (the lambda=1 kernel Pi+tK), Gamma_{1,0}=K, and W_2 = l^2(l-1) C - l^2 K^2/2 where, in the coordinates x1=F_a(xi)=1-aE1+a^2E2, x2=F_a'(xi), y1=G_a(eta), y2=G_a'(eta), Pi=x1y2, Pi'=x2y1 (x0=y0=1; Q is bihomogeneous of degree (l,l) in the three values (1,x1,x2),(1,y1,y2) at the vacuum point 0 and at a, a'):
  C = A_2 - K^2/2 = [ a^2a'^2(Pi^2-Pi'^2) - (Pi-Pi')(a^4 x2y2 + a'^4 x1y1) + (a-a')^2( a^2 x2y2 (x1-y1) + a'^2 x1y1 (y2-x2) ) ] / (a-a')^2,   K = aa'(Pi-Pi')/(a'-a).
  **W_3 (order 3)** = l^2(l-1)(l-2) A_3 + l^2(l-1) B_3 + l^2 K^3/3, with (terms `(coef) x1^i x2^j y1^k y2^m`, each polynomial = (a'-a)^3 * A_3 resp. B_3):
  A_3 (a'-a)^3 = (-a**3*(a - ap)**4/ap)*x1^0 x2^1 y1^1 y2^2 + (2*a**5*(a - ap)**2/ap)*x1^0 x2^2 y1^1 y2^2 + (-3*a**3*ap*(a - ap)**2)*x1^0 x2^2 y1^2 y2^1 + (-a**7/ap)*x1^0 x2^3 y1^1 y2^2 + (3*a**5*ap)*x1^0 x2^3 y1^2 y2^1 + (-5*a**3*ap**3/3)*x1^0 x2^3 y1^3 y2^0 + (a*ap*(a - ap)**4)*x1^1 x2^0 y1^1 y2^2 + (a**3*(a - ap)**4/ap)*x1^1 x2^1 y1^0 y2^2 + (-a**5*(a - ap)**2/ap)*x1^1 x2^1 y1^0 y2^3 + (-a*ap*(a - ap)**4)*x1^1 x2^1 y1^1 y2^1 + (-2*a**3*ap*(a - ap)**2)*x1^1 x2^1 y1^1 y2^2 + (4*a*ap**3*(a - ap)**2)*x1^1 x2^1 y1^2 y2^1 + (-a**5*(a - ap)**2/ap)*x1^1 x2^2 y1^0 y2^2 + (a**7/ap)*x1^1 x2^2 y1^0 y2^3 + (4*a**3*ap*(a - ap)**2)*x1^1 x2^2 y1^1 y2^1 + (-a**5*ap)*x1^1 x2^2 y1^1 y2^2 + (-3*a*ap**3*(a - ap)**2)*x1^1 x2^2 y1^2 y2^0 + (-4*a**3*ap**3)*x1^1 x2^2 y1^2 y2^1 + (3*a*ap**5)*x1^1 x2^2 y1^3 y2^0 + (-a*ap*(a - ap)**4)*x1^2 x2^0 y1^0 y2^2 + (a**3*ap*(a - ap)**2)*x1^2 x2^0 y1^0 y2^3 + (ap**3*(a - ap)**4/a)*x1^2 x2^0 y1^1 y2^1 + (-ap**5*(a - ap)**2/a)*x1^2 x2^0 y1^2 y2^1 + (a*ap*(a - ap)**4)*x1^2 x2^1 y1^0 y2^1 + (-2*a**5*ap)*x1^2 x2^1 y1^0 y2^3 + (-ap**3*(a - ap)**4/a)*x1^2 x2^1 y1^1 y2^0 + (-2*a*ap**3*(a - ap)**2)*x1^2 x2^1 y1^1 y2^1 + (5*a**3*ap**3)*x1^2 x2^1 y1^1 y2^2 + (2*ap**5*(a - ap)**2/a)*x1^2 x2^1 y1^2 y2^0 + (-a*ap**5)*x1^2 x2^1 y1^2 y2^1 + (-ap**7/a)*x1^2 x2^1 y1^3 y2^0 + (a*ap**3*(a - ap)**2)*x1^3 x2^0 y1^0 y2^2 + (2*a**3*ap**3/3)*x1^3 x2^0 y1^0 y2^3 + (-ap**5*(a - ap)**2/a)*x1^3 x2^0 y1^1 y2^1 + (-2*a*ap**5)*x1^3 x2^0 y1^1 y2^2 + (ap**7/a)*x1^3 x2^0 y1^2 y2^1
  B_3 (a'-a)^3 = (-4*a**3*ap*(a - ap)**2)*x1^0 x2^2 y1^2 y2^1 + (4*a**5*ap)*x1^0 x2^3 y1^2 y2^1 + (-3*a**3*ap**3)*x1^0 x2^3 y1^3 y2^0 + (a*ap*(a - ap)**4)*x1^1 x2^0 y1^1 y2^2 + (-a*ap*(a - ap)**4)*x1^1 x2^1 y1^1 y2^1 + (5*a*ap**3*(a - ap)**2)*x1^1 x2^1 y1^2 y2^1 + (5*a**3*ap*(a - ap)**2)*x1^1 x2^2 y1^1 y2^1 + (-5*a**5*ap)*x1^1 x2^2 y1^1 y2^2 + (-4*a*ap**3*(a - ap)**2)*x1^1 x2^2 y1^2 y2^0 + (-2*a**3*ap**3)*x1^1 x2^2 y1^2 y2^1 + (4*a*ap**5)*x1^1 x2^2 y1^3 y2^0 + (-a*ap*(a - ap)**4)*x1^2 x2^0 y1^0 y2^2 + (a**3*ap*(a - ap)**2)*x1^2 x2^0 y1^0 y2^3 + (-2*a*ap**3*(a - ap)**2)*x1^2 x2^0 y1^1 y2^2 + (a*ap*(a - ap)**4)*x1^2 x2^1 y1^0 y2^1 + (-2*a**3*ap*(a - ap)**2)*x1^2 x2^1 y1^0 y2^2 + (a**5*ap)*x1^2 x2^1 y1^0 y2^3 + (7*a**3*ap**3)*x1^2 x2^1 y1^1 y2^2 + (-5*a*ap**5)*x1^2 x2^1 y1^2 y2^1 + (a*ap**3*(a - ap)**2)*x1^3 x2^0 y1^0 y2^2 + (-2*a**3*ap**3)*x1^3 x2^0 y1^0 y2^3 + (a*ap**5)*x1^3 x2^0 y1^1 y2^2
  The l^4 coefficient A_3 is the order-3 term of the semiclassical action: Q ~ exp( l [log Pi + tau K/Pi + tau^2 C/Pi^2 + tau^3 A_3/Pi^3 + ...] ), tau = l t.
* **M1's A_l(z)** (coefficient of Pi^l = (UV)^l in Q, z = t aa'/(a'-a)) extended: l=1: 1+z; l=2: 1+4z+10z^2+12z^3; l=3: 1+9z+54z^2+222z^3+576z^4+720z^5; **l=4: 1+16z+168z^2+1296z^3+7416z^4+30240z^5+79200z^6+100800z^7** (z^2 coefficient l^2(l-1)(l+3)/2 holds; top coefficient = (2l-1)! C(2l-2,l-1) = 1, 12, 720, 100800; no small-factor factorisation for l=3,4; l=2: (1+2z)(1+2z+6z^2)).

### 2.2 The master kernel is a NORMAL-ORDERED EXPONENTIAL of a bidifferential operator (new, exact)
From the universality (section 1): with iota(k,l) = e2^k h_{l-k-1}(xi) (so level = weight with wt E1 = 1, wt E2 = 2; unitriangular change to monomials),
  **Q^{(l)}(t) = N(t)[ F_a(xi)^l G_{a'}(eta)^l ],  N(t) = sum_n t^n N_n universal (independent of l, a, a'),  N_n lowers the weight by n.**
N_n is a differential operator of bi-order (n,n) with polynomial coefficients. Its normal symbol (coefficients left, p_{E_i}=d/dE_i, p_{H_j}=d/dH_j) was found by exact linear fits (unique solutions, 0 free parameters; `mk2_opfit.py`, `mk2_fast.py`, `mk2_expfit.py`):
  **Phi(t) = exp(S),  S = t D_1 + t^2 D_2 + t^3 D_3 + t^4 D_4 + t^5 D_5 + ...,**  D_n polynomial in (E,H,p) of weight -n.
* D_1 = x := c_11 pE1 pH1 + c_12 (pE1 pH2 + pE2 pH1) + c_22 pE2 pH2 with c_11 = E1-H1, c_12 = E2-H2, c_22 = E2 H1 - E1 H2 (8 monomials, no lower-order terms). Applied to Pi^l: x[F^l G^l] = l^2 Pi^{l-1} K (= M1's order 1 = Casimir on the coherent state), for ALL l.
* D_2 = <x,x>_1 where <f,g>_1 := sum_{k in {E1,E2,H1,H2}} (df/dp_k)(dg/dq_k) (28 monomials, pure p-degree 3);
  D_3 = (1/2)( <D_2,D_1>_1 + <D_1,D_2>_1 ) (76 monomials, pure p-degree 4).
* **Tree-level recursion (verified n=2..5):** the top p-degree part (degree n+1) of D_n obeys (n-1) D_n^{tree} = sum_{i+j=n} < D_i^{tree}, D_j^{tree} >_1, i.e. generating function  t dS/dt - S = sum_k (dS/dp_k)(dS/dq_k) (a Hamilton-Jacobi-type PDE; with S = sigma(tp;q)/t: u.sigma_u - 2 sigma = sum_k sigma_{u_k} sigma_{q_k}). This is the l -> infinity action, matching the degree n+1 of W_n.
* **One-loop corrections start at n=4.** They are why the plain exponential e^{tX} of one operator X is wrong (exp(tX) with symbol x would give D_2 = (1/2) sum x_p x_q, the true D_2 is twice that).
  D_4 = D_4^{tree} + 3 (c_12^2 - c_11 c_22) pE2^2 pH2^2 (7 extra monomials of p-degree 4; equals (1/2)<D_1,D_3>_2 with <f,g>_2 := (1/2) sum_{k,m} f_{p_k p_m} g_{q_k q_m});
  D_5 = D_5^{tree} + 42 monomials of p-degree 5 (448 in all; expressible through <.,.>_2 brackets of D_1..D_4^{tree}, non-uniquely).
  Stored: `master_kernel2_Dsym.pkl` (strings of D_1..D_5), `master_kernel2_D3.pkl`, `master_kernel2_D{4,5}_gen.pkl` ({(coef exponents (E1,E2,H1,H2), p exponents): rational}).
* **Closed form of every order for every l:** the t^n coefficient of exp(S) is a sum of monomials coef(E,H) pE1^al1 pE2^al2 pH1^be1 pH2^be2, each of which acts on F_a^l G_{a'}^l as
  (l)_{|al|} (l)_{|be|} (-a)^{al1} a^{2 al2} (-a')^{be1} a'^{2 be2} F_a^{l-|al|} G_{a'}^{l-|be|}  ((l)_k falling factorial),
  hence **P_n(l) = sum_{i,j<=n} (l)_i (l)_j [weight-n polynomial]**, which explains "P_n polynomial of degree 2n in l" and the factors l^2(l-1)... (`mk2_opQ.py`).
  **Verified exactly against the independent exact kernels (symbolic E,H; a=1/3, a'=1/5): n=1..4 for l=1,2,3,4,5; n=5 for l=5 and l=6; n=1..5 at l=6 all True (l=6 kernel from the operator route, a=1/3,a'=1/5; `mk2_opQ56.log`).**
* D_4 fitted from N data of L=4 and L=5 (identical), D_5 from L=5 and L=6 (identical): 1544 unknowns, 39962 equations at L=6.
* Not found: a closed form of all D_n (the loop part); the one-loop source at n>=4 and a solution of the HJ equation (would give the classical action) are the next steps.

## 3 K7: the second operator and the span test
* Space V_l (x) V_l = Sym^l(X) (x) Sym^l(Y) in symbol coordinates (X0,X1,X2)=(1,E1,E2); Omega = sum X_i d_j (x) Y_j d_i (eigenvalue on [2l-k,k,0] = (c(mu)-2c(l))/2, c(mu)=sum mu_i(mu_i+4-2i); top irrep l^2). Checked: Q_1 = aa'/(a'-a)(l^2-Omega) v EXACT for l=2 (36-dim) and l=3 (100-dim).
* Component spans of {P_mu Q_n}_n reproduce M1: l=2: 1,2,2; l=3: 1,3,5,3 (W_Q = span Omega^k Q_n has dim 5 and 12).
* **Stabiliser test:** no operator of the bilinear family {1, E1_ij, E2_ij, E1_ij E2_kl} (rank 81 as operators) preserves W_Q beyond span{1, Omega}: l=3: 21 solutions in the family but only 2 independent operators (1 and Omega); l=2 degenerate (W_Q only 5-dim). So W_Q is not invariant under any simple extra operator.
* **The second operator is A := N_1 = D_1** (the universal first-order coupling = M2's A_1 transported through iota). It is not a gl3 operator and does not commute with Omega, but A v_{a,a'} = Q_1 = z(l^2-Omega) v for ALL (a,a'): on the coherent states it coincides with the Casimir. Cyclic subspace generated from v by words in {Omega, A}, graded by the number of A's (per-irrep dims):
  l=2: total dims 3,4,5,5,..; per irrep (1,1,1) -> (1,2,1) -> **(1,2,2) = W_Q exactly, closes at 2 A's**;
  l=3: total dims 4,6,9,12,13,13,..; per irrep (1,1,1,1) -> (1,2,2,1) -> (1,3,3,2) -> (1,3,5,3) -> (1,3,6,3): **W_Q (1,3,5,3) is contained but the closure has one extra vector in [4,2] (13 vs 12)**; the minimal number of A's to contain W_Q is 2(l-1).
  Verdict: Theta(z)v lies in the span of words in {Omega, A} applied to v for l=2,3 (A = D_1), but at l=3 with one extra vector, so a single A (with Omega) does not reproduce the component spans 1,3,5,3 exactly. The exact statement is the operator exponential of section 2.2 (A is only its first term; D_2 etc. are not polynomials in A and Omega). Not done: the z-dependence per irrep (R-matrix form), separate tests of N,F and difference operators E^(1)-E^(2) (the stabiliser search over the whole bilinear family covers them as first-order pieces and finds nothing new).

## Surprises
* Y(t) is graded, so ONE exact inverse gives all polynomial kernels (no grids), and the grading collapses everything to a single universal infinite matrix N with the l-dependence an explicit binomial.
* The kernel is a normal-ordered exp(S) of a bidifferential symbol with a tree-level HJ equation; the first loop correction appears only at t^4 (invisible at orders 0..3, which is why M1's low orders look "classical").
* "Plain exponential of one operator" is false exactly by the factor 2 in D_2.

## Not done
* Closed form of D_n for all n (only D_1..D_5 and the tree recursion); the one-loop source; the HJ solution.
* [5,1,0] and non-one-row; per-irrep z-dependence in K7; order 6+.
