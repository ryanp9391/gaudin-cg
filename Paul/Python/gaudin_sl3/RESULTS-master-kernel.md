# RESULTS-master-kernel (agent M1, round 5, Python only). PARTIAL: orders 0,1,2 closed for all lambda; full 4-index closed form NOT found.

Files: master_kernel.py (tab, symbol_orders, theta_from_symbol_orders, state_data, residuals), master_kernel_bld.py (point builder with
precision-scaled SVD tolerance), master_kernel_tab.pkl (R tables l1=1..4), mk_scratch*.py (exploration), mk_symbol_{1..4}.pkl, mk_Pn_*.pkl.

## 0 Conventions / support
R_{pq} = Theta^{[l1,0,0]}_{pq} / (theta_p(a) theta_q(a')) = sum F_{mm'} x^m y^m', x=at, y=a't, N=m+m'. Tabulated exactly for l1=1..4
(17, 193, 1231, 5062 entries). Verified: support = {N>=0, -N<=m<=2N (equivalently -N<=m'<=2N), m,m'<=2 l1, N<=4(l1-1) (l1>=2; N<=1 for l1=1)} (0 violations);
exchange R_ij(x,y)=R_ji(-y,-x) (0 violations). Single-index block = Gamma-form (G1).

## 1 Closed form through order 2 (all lambda, all pairs)
Symbol: sum_pq Theta_pq (x1^k x2^l - x1^l x2^k)(y1^k' y2^l' - ..) = (x2-x1)(y2-y1) Q(xi;eta). With F_c(xi)=prod(1-c xi_i), G_c(eta)=prod(1-c eta_j),
Pi=F_a G_a', K = aa'/(a'-a) (Pi - F_a' G_a), e=sym fns of xi, h of eta:
 Q_0=Pi^l;  Q_1 = l^2 Pi^(l-1) K;  Q_2 = l^2 (l-1) Pi^(l-2) (A + l K^2/2),
 A = -aa'(a-a')(e1-h1) Pi + K( K/2 + (a-a')(aa'(e1h1-e2-h2) - a e1 - a' h1 + 2) ).
Found as: Q_n^{(l)} is exactly divisible by Pi^(l-n) (n=1,2,3; l<=4); P_n:=Q_n/Pi^(l-n) is polynomial in l (n=1: l^2 K; n=2: l^2(l-1)(A+lB), B=K^2/2 exactly,
fit on l=2,3 and confirmed on l=4). Exact match with the pkl kernels l=2,3,4 (0 mismatches, orders<=1 and <=2). Order 1 equals the old alpha,beta,gamma,eps formula.
Order 3+: P_3 would need l=5 data (degree-6 polynomial in l, factor l^2(l-1)); not done.
Structure discovered (useful for the next step): (i) Q is a polynomial in the values U=P(w),U'=P(w'),V=Q(w'),V'=Q(w) of the monic quadratics P=prod(x-xi), Q=prod(x-eta)
at w=1/a,w'=1/a' plus the leading coeff 1 (a third point, w=inf = the vacuum a0=0); lambda=1: (w w')^2 Q = (1+z)UV - zU'V', z=t/(w-w'). Coefficients depend only on z (nu=w-w' fixed by weights).
(ii) Theta lives on V_[l,0,0] x V_[l,0,0]: the pair space is Sym^l(C^3) (dim C(l+2,2) = number of pairs) and the symbol variables (1,e1,e2) are its coordinates; theta(a)=l_a^l, l_a=(1,-a,a^2) coherent state.
Order 1 is exactly the gl(3) Casimir: Q_1 = aa'/(a'-a) (l^2 - Omega) v, Omega=sum x_i d_j (x) y_j d_i, v=l_a^l l_a'^l. But Theta(z) is NOT a function of Omega acting on v (l=2: components of Theta in the 3 irreps span 1,2,2 dims; l=3: 1,3,5,3): a second vector per irrep is needed.
(iii) top order: Q_{4(l-1)} = c_l (aa')^(2l-2) prod_{ij}(xi_i-eta_j) (c = 12, 8640 for l=2,3); every t^n, n>=1, vanishes on {xi=eta as sets}.
(iv) A_l(z) = coefficient of (UV)^l: l=1: 1+z; 2: 1+4z+10z^2+12z^3; 3: 1+9z+54z^2+222z^3+576z^4+720z^5 (z^2 coeff l^2(l-1)(l+3)/2).

## 2 [5,0,0] test (no fit; 441 states, 40 digits, build ~67 s per point after fixing SVD tolerance to 1e-25; G1's 1e-50 fails at 40 digits)
Predicted Theta truncated at t^2 vs exact jet LHS, state-wise max relative residual (median):
 theta12=100 (t=.01): order<=0: 0.40; <=1: 9.8e-2 (4.3e-3); <=2: 1.3e-2 (2.9e-4).
 theta12=1000 (t=.001): <=0: 4.5e-2; <=1: 1.3e-3 (4.7e-5); <=2: 2.0e-5 (3.2e-7).
Residual scaling t^1, t^2, t^3 respectively (ratios ~10, 78, 650): orders 0,1,2 of the closed form are CORRECT for [5,0,0], and the first missing order is t^3. Full [5,0,0] prediction NOT possible (order>=3 unknown). Not tested: [5,1,0]. 
(Coordinator note: M2's operator route gives exact [5,0,0] kernel at rational points; use it as target for orders>=3.)

## 3 A3
(a) theta(a+s)=exp(s N^{wedge2}) theta(a) exactly (l<=5) with N e_k = -(n-k) e_{k+1}, n=l1+1 (N_{k+1,k}=-(n-k)); F e_k = k e_{k-1}, H=n-2k complete an sl(2) (spin n/2).
(b) Diagonal E-shift a->a+s, a'->a'+s: fails for t'=t for all l; for l=1 holds with t'=t aa'/((a+s)(a'+s)); for l>=2 no scalar t' works (ratio test). Reason found: Theta depends on a third point (vacuum a0=0). What IS exact (l=1,2,3, checked): translation covariance of w=1/a with t FIXED:
 [(1-as)(1-a's)]^(2l) Theta(a/(1-as), a'/(1-a's), t) = (e^{sF}x e^{sF}) Theta(a,a',t) (both signs of s). Together with scaling this reduces Theta to a one-variable matrix family C^(l)(z), z=t/(w-w').
(c) Fusion: the naive "product of l=1 exchange factors (1+z_su(1-S_su)) with coincident limit taken afterwards" has no derivative terms and cannot work: the l=2 symbol is not in the span of f(a)^2g(a')^2, f f' g g', f'^2 g^2 (checked). Weighted-subset-sum ansatz (weights (1+z), -z by crossing) also fails numerically (finite limit, wrong value). Casimir-function ansatz Theta=f(Omega)v fails from order 2. Hirota/Toda bilinear tests in (a,a'),(t) failed. Not resolved.

## Surprises / not done
Not done: general 4-index closed form, orders>=3, [5,1,0], l=5 full reconstruction. Suggested next: get Theta^{[5,0,0]} exactly (M2 operator route), fix P_3(l), look for the second-vector/operator completing the Casimir picture (state space V_l x V_l, three-point/Moebius structure above).
