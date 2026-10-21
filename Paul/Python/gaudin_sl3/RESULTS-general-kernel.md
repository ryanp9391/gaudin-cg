# RESULTS-general-kernel: the general-Lambda structure of the kernel of the two-vacua jet formula, sl(3) Gaudin (agent G1, 2026-10-16)

Files: `general_kernel.py` (`python -X utf8 general_kernel.py --stages K1|REP|ZC|BLK|SGL|ANALYZE ...`; Python only, 80 digits, exact Fractions for the one-site part),
caches `general_kernel_pt_*.pkl` (states + dual Q's per rep/zeta/theta), `general_kernel_rep_<lam>.pkl` (reconstructed exact polynomial kernels), `general_kernel_K1.pkl`,
printed output of the analysis stage `general_kernel_output_analyze.txt`. Helper `general_kernel_bld.py` builds one point. Imports (read-only): `qside.py` (ODE solver/PF),
`gaudin_chain.py` (two-site chain, eigenstates), `twovac_fit.py` (`PivQR`, `rat_recon`, `mpf`), `gl3_reps.py`. No Wolfram, no notebooks, no git.

## 0. Conventions (binding for everything below)

* Jet formula: `c Qnum M1 M2 = sum_{(kl),(k'l')} K_{(kl),(k'l')} m_(kl)(th1) m_(k'l')(th2)`, `m_(kl) = P3 p_kl`, `p_kl = P2^(k) P3^(l) - P3^(k) P2^(l)` at theta_s,
  pairs `0 <= k < l <= l1+1` NOT forced (`l >= l2+1`; count `C(l1+2,2) - C(l2+1,2)`), `M = P3 p_{0,l2+1}` EXACTLY (no factor 4), `c` as in the task,
  `Qnum = <Phi_inf|Psi>/<x0^G|Psi>`. **Canonical kernel `Theta := K/(w1 w1')^(2 l1 - l2)`**, `(a,a') = (1/w1,1/w1')`, `t = 1/theta12`. The common factor `e^{(z2+2 z3) theta_s}` of all monomials
  of a site is dropped (formula homogeneous of degree 1 in the monomials), so everything is rational. Jet weight of a pair `e_kl := k + l - l2 - 1`.
* **The "4" of D2/D3**: D2's `G = 4K` (its `M_vac = 2 h^2 P3 p02` per site); with `M = P3 p_{0,l2+1}` exactly one has **Theta = G/(4 w^3 w'^3)** (adjoint, checked symbolically in all 25 entries, below). I use no 4 anywhere.
* Pair order of the pairs (kl): lexicographic. Adjoint: `(02),(03),(12),(13),(23)` = D2's `N5 = (P3p02, P3p03, P3'p02, P3'p03, P3p23)` re-ordered (`P3'p0l = kappa P3 p0l = P3 p1l` since `P3' = kappa P3`, `p_1l = kappa p_0l`).
  [3,2,0]: `(03),(04),(13),(14),(23),(24),(34)` = P4's order. `Plucker pair` and `N5`-type labels are therefore the same thing.
* Pipeline: states from `Chain(lam, theta)` (80 digits, eigen-residual <= 4e-75), dual Q's `P^j = e^{z_j u} p_j` (monic, deg p_j = mu_j checked for every state, kernel dim 1: 0 bad states in every point),
  kernel by the column-pivoted fixed-point QR of D1 (row-scaled), QR once per (zeta-points, theta), then `kernel(w1,w1')` for many omegas by re-using the factorisation. **Polynomial reconstruction**: tensor-product
  least squares on a grid (w1 = 5/2 fixed) x (w1' in 7..11 values) x (theta12 in 7..19 values, theta = (2/7, 2/7 - theta12)), bidegree `(2 l1 - l2, Dy)` in `(1/w1', 1/theta12)`, rational reconstruction of the coefficients;
  the dependence on the third variable `a` is then fixed by the homogeneity law (below) and **independently checked at a = 1/w with w = -2/5 and 17/3 (not used in the fit)**: worst relative deviation 5e-78 ([1,1,0]) ... 1.3e-64 ([4,4,0]), see section K0.

## Headline results

1. **One-site kernel in closed form, exact, 17 representations** (K1): `theta_kl(a) = c_kl a^{k+l-l2-1}`,
   `c_kl = (-1)^{k+l-l2-1} C(l1-l2, l-l2-1) C(l1+1, k) (l-k) / ((l2+1) C(l, l2+1))`; equivalently the generating symbol
   `sum_{k<l} theta^{[l1,0,0]}_kl (x1^k x2^l - x1^l x2^k) = (x2 - x1) [(1-a x1)(1-a x2)]^{l1}`, i.e. **`sum_kl theta_kl(a) m_kl = eps a^{-l2} P3 W[(1-a d)^{l1} P2, (1-a d)^{l1} P3](theta)`**, `eps = (-1)^{l2} l2!(l1-l2)!/l1!`, `W` the Wronskian: the one-site kernel is a *decomposable bivector* (a 2-plane), a dressed Wronskian of the two dual Q's.
2. **Block theorem (new, exact, verified)**: `Theta^{[l1,l2,0]}_{ij}(a,a',t) = [l2!(l1-l2)!/l1!]^2 (a a')^{-l2} Theta^{[l1,0,0]}_{ij}(a,a',t)` for all non-forced pairs i,j. All kernels of a given l1 are sub-blocks of the symmetric-representation kernel. Verified symbolically (exact polynomial identity) for l1 <= 3 (all l2) and
   for l1 = 4: l2 = 4 symbolically (polynomial identity with the reconstructed [4,0,0]), and by direct independent fits at two theta12 x two (w,w') points for l2 = 1,2,3 (dev <= 3e-58). Also the *corner block* (pairs `(k,l1+1)`) of every `[l1,l2,0]` equals `(l2!(l1-l2)!/l1!)^2 (a a')^{l1-l2}` times the closed form of the single-index kernel `Theta^{[l1,l1,0]}` (C3): direct fits [3,1,0], [3,2,0], [4,1,0] (dev <= 1.3e-62).
3. **`Theta(t=0) = theta(a) x theta(a')` exactly** in all computed two-site reps (symbolic polynomial identity); `Theta` is a polynomial in `(a,a',t)` with rational coefficients, homogeneous: every monomial `a^p a'^q t^n` of the (i,j) entry has `p + q - n = e_i + e_j`, `0 <= p,q <= 2 l1 - l2` (max attained), exchange rule `Theta_ij(a,a',t) = Theta_ji(a',a,-t)` exact.
4. **Closed form for the whole family `[l,l,0]`** (single Plucker index, verified exactly for l = 1..4, all entries, i.e. 4+9+16+25 entries, and **predicted correctly for l = 5**: direct fit of [5,5,0] agrees to <= 9.2e-58 at 6 points):
   `Theta_{kk'} = (-1)^{k+k'} C(l,k) C(l,k') a^k a'^{k'} sum_{m=0}^{l-k} sum_{m'=0}^{l-k'} C(l-k,m) C(l-k',m') F_{mm'}(k,k') (a t)^m (-a' t)^{m'}`,
   `F_{00} = 1, F_{m0} = (k')_m, F_{0m'} = (k)_{m'}, F_{mm'} = (k+m+1)_{m'-1} (k'+m'+1)_{m-1} (k k' + m' k + m k')`
   (uniformly `F = Gamma(k+m+m') Gamma(k'+m+m') [(k+m)(k'+m') - m m'] / (Gamma(k+m+1) Gamma(k'+m'+1))`).
5. **Top-state sum rule, exact for all computed reps and with the correct right-hand side**: `sum_ij Theta_ij m_i^top(th1) m_j^top(th2) = M^top(th1) M^top(th2) R(a) R(a')`, `R(a) = (1 - z2 a)^{l1-l2} (1 - z3 a)^{l1}`. **The spec's literal right-hand side (without R) is refuted** (the top state is not `e1 x e1` in the companion frame; `c a^{2l1-l2}/x_hw = R(a)` exactly, `x_hw = (z2-z1)^{l2}`).
6. The l2 = 0 probes work: structure unchanged, exponents `{0,1,l1+2}`, `c2 = 0` (the S2 double-pole eigenvalue is the scalar 0, *not* non-scalar), kernel unique, same laws, and they are the "master" kernels of the block theorem.
7. **Exact first-order term for all pairs and all (l1,l2)** (K3b, `theta1_closed`): `Theta_1 = t theta_p theta_q [alpha a + beta a' + gamma a^2/a' + eps a'^2/a]` with explicit `alpha, beta` (rank 2), `gamma = k'l'(n-k)(n-l)/n^2`, `eps = -kl(n-k')(n-l')/n^2` (`n = l1+1`); 0 mismatches over 1062 + 675 tests, held-out reps included. Together with C1 and C3 this fixes the kernel to first order everywhere and to all orders on the corner block.
8. Degree laws: a-degree `2 l1 - l2` (attained), t-degree `4(l1-1) - max(0, 2 l2 - 3)` (l1 >= 2; `2l-1` for `[l,l,0]`); `Theta_n` ranks `1,2,4,8,...`; `Theta(xi = eta)` (plane waves) is t-independent in the single-index symbol.

## K0. The three known kernels in the common labelling

* **Adjoint [2,1,0]** (pairs `(02),(03),(12),(13),(23)`): `Theta = G/(4 w^3 w'^3)` as an exact rational-function identity in 25/25 entries (D2's G parsed from `RESULTS-twovac-limit.md`, index map D2-> mine: 1->(02), 2->(12), 3->(03), 4->(13), 5->(23)). Fit at the base point: 25 unknowns, rank 25, 192 rows (3 zeta), residual 0.0 (floor 1e-90). t-degree 4, a-degree 3. Full table (`Theta[(kl),(k'l')]`, a = 1/w1, ap = a', t = 1/theta12):

```
Theta[(02),(02)] = 1 + (a - ap) t - 3 a ap t^2                                   Theta[(23),(23)] = a^3 ap^3 / 4
Theta[(02),(03)] = -ap (2 a t + 1)/2                                              Theta[(03),(03)] = a ap/4
Theta[(03),(02)] =  a (2 ap t - 1)/2                                              Theta[(13),(13)] = -a^2 ap^2 (3 a ap t^2 - a t + ap t - 1)
Theta[(02),(12)] = (6 a^3 ap t^3 - 2 a^3 t^2 + 24 a^2 ap^2 t^3 - 2 a^2 t + 6 a ap^3 t^3 + 12 a ap^2 t^2 - 6 a ap t + 2 ap^3 t^2 + 2 ap^2 t - 3 ap)/2
Theta[(02),(13)] = -ap (3 a^3 ap^2 t^4 + 2 a^3 ap t^3 - a^3 t^2 + 4 a^2 ap^2 t^3 - 3 a^2 ap t^2 - a^2 t + a ap^2 t^2 - 3 a ap t - ap)
Theta[(02),(23)] = -ap^2 (6 a^3 ap t^3 + 4 a^3 t^2 + 8 a^2 ap t^2 + 2 a^2 t + 4 a ap t + ap)/2
Theta[(03),(12)] = a ap (12 a^2 ap^2 t^4 - 16 a^2 ap t^3 + 4 a^2 t^2 - 8 a ap^2 t^3 - 12 a ap t^2 + 4 a t - 4 ap^2 t^2 - 4 ap t + 3)/4
Theta[(03),(13)] = a ap^2 (6 a^2 ap t^3 - 2 a^2 t^2 + 2 a ap t^2 - 2 a t - 1)/2 ;  Theta[(03),(23)] = a ap^3 (6 a^2 t^2 + 4 a t + 1)/4
Theta[(12),(12)] = -(12 a^3 ap t^2 - 4 a^3 t + 48 a^2 ap^2 t^2 - 12 a^2 ap t + 12 a ap^3 t^2 + 12 a ap^2 t - 9 a ap + 4 ap^3 t)/4
Theta[(12),(13)] = a ap (6 a^2 ap^2 t^3 + 4 a^2 ap t^2 - 2 a^2 t + 8 a ap^2 t^2 - 6 a ap t + 2 ap^2 t - 3 ap)/2 ;  Theta[(12),(23)] = a ap^2 (6 a^2 ap t^2 + 4 a^2 t + 8 a ap t + 3 ap)/4
... the remaining 11 entries follow from the exchange rule Theta_ij(a,ap,t) = Theta_ji(ap,a,-t) (all 25 stored exactly in general_kernel_rep_210.pkl).
```
  Normalisation check: the t=0 limit `Theta_(02)(02) = 1`, `Theta_(23)(23) = a^3 a'^3/4 = (1/4)(a a')^3` is `theta_(23)(a) theta_(23)(a')` with `theta_(23) = -a^3/2`.
* **[3,2,0]** (pairs `(03),(04),(13),(14),(23),(24),(34)`, 49 unknowns, 225 rows x 3 zeta, rank 49, residual 0.0): my reconstructed polynomial law evaluated at `(w1,w1',theta12) = (5/2,7/4,31/35)` reproduces P4's printed base-point kernel exactly
  (`K11 = -2477531875/3936256, K12 = -15649375/95232, K17 = 1241317375/2770563, K71 = -3124391459/709264128, K77 = 1/9`, `K = Theta (w w')^4`). Laws: a-degree 4 (= 2 l1 - l2), `t`-degree 7, entry (i,j) homogeneous of degree `e_i + e_j` with `e = (0,1,1,2,2,3,4)`, `#monomials 458`. By the block theorem `Theta^{[3,2,0]} = (1/9)(a a')^{-2} Theta^{[3,0,0]}|`.
* **[3,1,0]** (new polynomial laws; pairs `(02),(03),(04),(12),(13),(14),(23),(24),(34)`, 81 unknowns, rank 81, residual 0.0, 8 x 13 grid, bidegree (5,10) fit residual 4.9e-68, all coefficients rational): a-degree 5, t-degree 8, `e = (0,1,2,1,2,3,3,4,5)`, 998 monomials, `Theta^{[3,1,0]} = (1/9)(a a')^{-1} Theta^{[3,0,0]}|`.
* Polynomiality / homogeneity / exchange verdicts for all computed reps (output of `--stages ANALYZE`):

| Lambda | pairs | a-deg = a'-deg (= 2l1-l2) | t-deg | exchange rule | Theta(t=0) = theta x theta' | rank Theta generic | ranks of Theta_n (n = 0,1,..) | top-state sum rule |
|---|---|---|---|---|---|---|---|---|
| [1,0,0] | 3 | 2 | 1 | exact | exact | 2 | 1,2 | exact |
| [1,1,0] | 2 | 1 | 1 | exact | exact | 2 | 1,2 | exact |
| [2,0,0] | 6 | 4 | 4 | exact | exact | 6 | 1,2,4,6,6 | exact |
| [2,1,0] | 5 | 3 | 4 | exact | exact | 5 | 1,2,4,5,4 | exact |
| [2,2,0] | 3 | 2 | 3 | exact | exact | 3 | 1,2,3,2 | exact |
| [3,0,0] | 10 | 6 | 8 | exact | exact | 8 | 1,2,4,8,10,10,10,8,6 | exact |
| [3,1,0] | 9 | 5 | 8 | exact | exact | 8 | 1,2,4,8,9,9,8,6,4 | exact |
| [3,2,0] | 7 | 4 | 7 | exact | exact | 7 | 1,2,4,7,7,6,4,2 | exact |
| [3,3,0] | 4 | 3 | 5 | exact | exact | 4 | 1,2,3,4,3,2 | exact |
| [4,0,0] | 15 | 8 | 12 | exact | exact | 12 | 1,2,4,8,13,15,15,15,15,12,10,8,6 | exact |
| [4,4,0] | 5 | 4 | 7 | exact | exact | 5 | 1,2,3,4,5,4,3,2 | exact |

  ("exact" = polynomial identity in the stored rational coefficients; top-state sum rule at 3-4 random rational points of `(z2,z3,th1,th2,a,a')`, deviation exactly 0.)
* **Homogeneity (XXX torus law) per entry**: `Theta_ij = sum_n t^n Theta_{ij,n}(a,a')`, `Theta_{ij,n}` homogeneous of degree `e_i + e_j + n` in `(a,a')`, i.e. `Theta_ij(la,la',lam theta12) = lam^{e_i+e_j} Theta_ij` (w1, w1' scale like zeta, theta12 like 1/zeta; `e_kl = k+l-l2-1`, weight of the monomial `P3 p_kl` relative to M).
  Independent check at `w = -2/5` and `17/3` (not used in the reconstruction): worst rel. deviation 5.0e-78 [1,1,0], 5.4e-77 [1,0,0], 3.4e-72 [2,1,0], 6.4e-74 [2,2,0], 1.7e-72 [2,0,0], 7.6e-70 [3,3,0], 4.8e-67 [3,0,0], 5.8e-67 [3,1,0], 7.4e-68 [3,2,0], 1.3e-64 [4,4,0], 3.0e-62 [4,0,0] (relative differences between the fixed-point-QR kernels and the polynomial law at held-out (w1',theta12) points and two values of w1; the law for [3,2,0] was in addition confirmed against P4's printed numbers).
* **Dependence on theta** only through `theta12`: three pairs `(theta1,theta2) = (2/7,-3/5), (-1/3,-128/105), (5/4,51/140)` with the same `theta12 = 31/35` give the same Theta: max abs deviation 1.9e-74 ([2,1,0]), 8.1e-77 ([2,0,0], 2 zeta stacked), 2.4e-68 ([3,2,0]).

## K1. The one-site kernel (exact)

Method: one site, twist `diag(z)`, `L(u) = g + E/(u-theta)`, `S_k` from the quantum characteristic determinant (same symbolic `cdet` as the chain), partial-fraction coefficient operators fitted exactly (Fractions, residual 0),
eigenstates = GT weight vectors (an exact eigenvector of all coefficient operators in the GT basis: 6/8, 9/15, 12/24, 12/27, ... states; **the other states (weights with multiplicity > 1) are not GT eigenvectors and were not used**), dual Q's from the exact ODE (monic, nullity 1, `deg p_j = nu_j` for every used state),
jets in Fractions, `Omega^inf[z,w]` and `x(z)` in Fractions. Fit of `c1 Qnum M = sum_kl K_kl m_kl` over 8 zeta points and all diagonal states: **consistent (residual exactly 0) with rank = number of pairs for every Lambda, K_kl = w^{2 l1 - l2 - e_kl} c_kl with c_kl independent of w (4 values) and of zeta (8 values)**. One-site local exponents/Frobenius: `{0,l2+1,l1+2}` at theta; `c2 = l2(l1+1)`.

Verified exactly for 17 representations: [1,0,0] [1,1,0] [2,0,0] [2,1,0] [2,2,0] [3,0,0] [3,1,0] [3,2,0] [3,3,0] [4,0,0] [4,1,0] [4,2,0] [4,3,0] [4,4,0] [5,1,0] [5,2,0] [5,3,0] (the nine asked-for ones [1,0,0] [1,1,0] [2,0,0] [2,1,0] [2,2,0] [3,1,0] [3,2,0] [4,1,0] [4,2,0] are all included); the closed form was guessed from the Lambda <= [4,2,0] data
and **then correctly predicted [3,0,0], [3,3,0], [4,0,0], [4,3,0], [4,4,0], [5,1,0], [5,2,0], [5,3,0]**:

`c_kl = (-1)^{k+l-l2-1} * C(l1-l2, l-l2-1) * C(l1+1, k) * (l-k) / ((l2+1) * C(l, l2+1))`, `theta_kl(a) = c_kl a^{k+l-l2-1}`  (`c_{0,l2+1} = 1`).

Taylor-normalised form: with `p~_kl = p_kl/(k! l!)`, `c~_kl = (-1)^{l2+1} (l1-l2)! l2! (l1+1)! (-1)^{k+l} (l-k) / ((n-k)! (n-l)!)`, `n = l1+1`:
`(l-k)/((n-k)!(n-l)!) = 1/((n-k-1)!(n-l)!) - 1/((n-k)!(n-l-1)!)`, so the kernel is `g^ ^ h^` with `g_k = (-a)^k/(n-1-k)!`, `h_k = (-a)^k/(n-k)!`: a decomposable bivector. In jet terms
`sum_kl theta_kl m_kl = eps a^{-l2} P3 W[(1-a d)^{l1} P2, (1-a d)^{l1} P3](theta)`, `eps = (-1)^{l2} l2!(l1-l2)!/l1!` (symbol check `sum theta^{[l1,0,0]}_kl (x1^k x2^l - x1^l x2^k) = (x2-x1)((1-a x1)(1-a x2))^{l1}`, exact for l1 = 1..5, and `c^{(l2)}_kl/c^{(0)}_kl = (-1)^{l2} l2!(l1-l2)!/l1!` for all pairs, l1 <= 5).
Examples (pairs in lexicographic order):
```
[2,1,0]: (02) 1, (03) -1/2, (12) -3/2, (13) 1, (23) -1/2          [1,1,0]: (02) 1, (12) -1          [2,2,0]: (03) 1, (13) -2, (23) 1       [1,0,0]: (01) 1, (02) -1, (12) 1
[3,1,0]: (02) 1,(03) -1,(04) 1/3,(12) -2,(13) 8/3,(14) -1,(23) -2,(24) 1,(34) -1/3                        [3,2,0]: (03) 1,(04) -1/3,(13) -8/3,(14) 1,(23) 2,(24) -1,(34) 1/3
[4,2,0]: (03) 1,(04) -2/3,(05) 1/6,(13) -10/3,(14) 5/2,(15) -2/3,(23) 10/3,(24) -10/3,(25) 1,(34) 5/3,(35) -2/3,(45) 1/6
[2,0,0]: (01) 1,(02) -2,(03) 1,(12) 3,(13) -2,(23) 1
```
Exactly diagonalised GT states used per zeta (dim in brackets): [1,0,0] 3 (3), [1,1,0] 3 (3), [2,0,0] 6 (6), [2,1,0] 6 (8), [2,2,0] 6 (6), [3,0,0] 10 (10), [3,1,0] 9 (15), [3,2,0] 9 (15), [3,3,0] 10 (10), [4,0,0] 15 (15), [4,1,0] 12 (24), [4,2,0] 12 (27), [4,3,0] 12 (24), [4,4,0] 15 (15), [5,1,0] 15 (35), [5,2,0] 15 (42), [5,3,0] 15 (42); 8 zeta points each, K unknowns = number of pairs, always rank-complete and consistent.
(full list for all 17 reps and the equality with the closed form: `general_kernel_output_analyze.txt` and `python general_kernel.py --stages K1`.)

**One-site top-state sum rule** (exact, symbolic in `z2, z3, theta`, 10 reps): `sum_kl theta_kl(a) m^top_kl = M^top (1 - z2 a)^{l1-l2} (1 - z3 a)^{l1}` (with `P3 = e^{z3 u}`, `P2 = e^{z2 u} p2`, `p2 = sum_k (z3-z2)^{-k} (Pi^{l2})^{(k)}`, `Pi = u - theta`).
Sign/normalisation of `p2` as proposed in the spec is right (monic, `p2' - (z3-z2) p2 = -(z3-z2) Pi^{l2}`); the right-hand side is not `M^top` because `<x0^c| ` is not the top state: `x_hw = (z2-z1)^{l2}` and `c1 a^{2l1-l2}/x_hw = R(a)` (checked exactly for 10 reps).

## K2. More two-site representations

Fits (all: residual **0.0** = below the 1e-90 fixed-point floor, 80-digit data, rank = all `pairs^2` unknowns so the kernel is unique, all kernels rational, 0 bad states, eigen-residuals <= 1e-74):

| Lambda | dim | states | pairs | unknowns | rows/zeta | rank | grid (w1' x theta12) | bidegree fit, worst rel. residual | non-rational coefficients |
|---|---|---|---|---|---|---|---|---|---|
| [1,0,0] | 3 | 9 | 3 | 9 | 9 | 9 | 5 x 7 | (2,4): 4.8e-77 | 0 |
| [1,1,0] | 3 | 9 | 2 | 4 | 9 | 4 | 6 x 6 | (1,3): 8e-78 | 0 |
| [2,0,0] | 6 | 36 | 6 | 36 | 36 | 36 | 7 x 11 | (4,8): 1.1e-73 | 0 |
| [2,2,0] | 6 | 36 | 3 | 9 | 36 | 9 | 5 x 7 | (2,4): 1.3e-74 | 0 |
| [2,1,0] | 8 | 64 | 5 | 25 | 64 | 25 | 6 x 9 | (3,6): 5.7e-73 | 0 |
| [3,0,0] | 10 | 100 | 10 | 100 | 100 | 100 | 9 x 13 | (6,10): 7.8e-68 | 0 |
| [3,3,0] | 10 | 100 | 4 | 16 | 100 | 16 | 6 x 11 | (3,8): 3.1e-69 | 0 |
| [4,0,0] | 15 | 225 | 15 | 225 | 225 (x2 zeta = 450) | 225 | 11 x 19 | (8,15): 1.4e-66 | 0 |
| [3,1,0] | 15 | 225 | 9 | 81 | 225 | 81 | 8 x 13 | (5,10): 4.9e-68 | 0 |
| [3,2,0] | 15 | 225 | 7 | 49 | 225 | 49 | 7 x 11 | (4,8): 1.2e-68 | 0 |
| [4,4,0] | 15 | 225 | 5 | 25 | 225 | 25 | 7 x 11 | (4,8): 7.7e-66 | 0 |
| [4,1,0] | 24 | 576 | 14 | 196 | 576 | 196 | base theta12 + block-theorem points | -- | -- |
| [4,2,0] | 27 | 729 | 12 | 144 | 729 | 144 | idem | -- | -- |
| [4,3,0] | 24 | 576 | 9 | 81 | 576 | 81 | idem | -- | -- |

* Build times per point: [3,1,0]/[3,2,0] ~ 35 s; [4,1,0] 3 s (ops) + 68 s (eigen) + 58 s (Q) ~ 2 min, QR 8 s; [4,2,0] ~ 4 min, QR 6 s. The expected cost problem of [4,1,0] (576 states, 196 unknowns) did not materialise: rank 196, residual 0.0 at the base point in 8 s.
  The full (a',t) polynomial grid for [4,1,0] was **not** done (replaced by the block-theorem check against [4,0,0], which is the stronger statement; see K4), nor was a direct [4,x,0] reconstruction except [4,4,0] and [4,0,0] (below).
* **zeta-independence and the formula for l2 = 0** (stage ZC): kernels from 3 zeta points stacked (`P, P', X1`; overdetermined; for [1,0,0] 27 rows vs 9 unknowns, [2,0,0] 108 vs 36, [3,0,0] 300 vs 100) have residual 0.0, rank full, and the single-zeta kernels agree with the stacked one to <= 1.7e-68 ([1,0,0]: 4.8e-78, [2,0,0]: 1.5e-74, [3,0,0]: 1.7e-68, [3,1,0]: 8e-69, [3,2,0]: 4.5e-69, [3,3,0]: 2.6e-69, [1,1,0], [2,2,0], [2,1,0] <= 1e-74). So the jet formula holds for the one-row (l2 = 0) representations with the same monomials.
* **l2 = 0 probes** ([1,0,0], [2,0,0], and [3,0,0]): the oper has `c2 = l2(l1+1) = 0`: the S2 double-pole eigenvalue is the scalar 0 on all states (not non-scalar in the sense of the spec; checked on 9/36/100 states) and S3 has no triple pole; local exponents (rank jumps of the 3 x m dual jet matrix on 18+ (state,site) cases) `{0,1,3}` for [1,0,0], `{0,1,4}` for [2,0,0] (= `{0,l2+1,l1+2}`); no Plucker coordinate forced (`p_01 != 0`, `M = P3 p_01`); count `C(l1+2,2)` (3, 6, 10). Nothing breaks. (Compare [2,1,0]: `{0,2,4}`, c2 = 3.)
* **theta-dependence**: Theta depends on (theta1,theta2) only through theta12 (checks in K0 above).


## K3. Structure across representations

**(a) Factorisation at 1/theta12 = 0.** `Theta(t=0) = theta(a) x theta(a')` with `theta` the exact one-site kernel of K1: symbolic identity entry by entry for all 10 reps in the table (including `Theta_{(23),(23)}^{adj}(0) = (a a')^3/4`, and the l2 = 0 ones).

**(b) Expansion `Theta = sum_n t^n Theta_n`.** Ranks of `Theta_n` (generic (a,a')), listed in the table above: `1, 2, 4, 8, ...` at the first orders (rank 2^n until the number of pairs caps it), then a plateau and a symmetric decrease
(for `[l,l,0]`: `rank Theta_n = min(n+1, 2l-n)`, n = 0..2l-1, i.e. `1,2,3,..` -- l = 1..4 checked). Supports: nonzero entries of `Theta_n` per order are listed by `python ranks.py` style in the notes (e.g. [2,1,0]: 25,23,19,12,4 of 25 entries; [3,2,0]: 49,47,43,37,29,19,9,2 of 49): the top order `t^{Dt}` has only 2-4 nonzero entries.
t-degree `Dt` (all computed reps + the l2-blocks of [4,0,0] via the block theorem): `Dt = 4(l1-1) - max(0, 2 l2 - 3)` for `l1 >= 2` (`l1 = 1: Dt = 1`): `[2,l2,0]: 4,4,3`, `[3,l2,0]: 8,8,7,5`, `[4,l2,0]: 12,12,11,9,7`; for `[l,l,0]`: `2l-1`; the one-row sequence `1,4,8,12` is `4(l-1)` from l = 2 on (my earlier guess `C(l+2,2)-2 = 13` for [4,0,0] is refuted: 12). Prediction: `[5,0,0], [5,1,0]: 16`, `[5,2,0]: 15`, `[5,3,0]: 13`.
**Exact first-order structure (new closed form, verified: 0 mismatches in all computed reps incl. the held-out [4,4,0] and [4,0,0] (675 entries), 1062 random-point entry tests for l1 <= 4).** With `n = l1+1`, pairs `p = (k,l)`, `q = (k',l')`, `rho(p) = (n-k)(n-l)`, `sigma(p) = k+l-1`, `f2(p) = 2(n-1) - sigma(p)`, `N0 = n(n-1)`:
```
Theta_1,pq = t theta_p(a) theta_q(a') [ alpha_pq a + beta_pq a' + gamma_pq a^2/a' + eps_pq a'^2/a ]
alpha_pq = -(3/n^2) rho(p)(rho(q)-N0) - (2/n) rho(p) sigma(q) + (1/n) f2(p)(rho(q)-N0) + f2(p) sigma(q)
beta_pq  = -alpha_qp
gamma_pq =  k' l' rho(p)/n^2           eps_pq = - k l rho(q)/n^2
```
(`general_kernel.theta1_closed`; the formula was fitted on l1 = 1,2,3 and then reproduced the l1 = 4 single-index kernel and all l2-blocks exactly, which is the l2-independence C1 at work). As matrices on the pair space `alpha, beta` have **rank 2** and `gamma, eps` rank 1 (so `Theta_1` has rank 2: the four terms recombine); the factors `(n-k)(n-l)` and `k'l'` are the products of the raising amplitudes `n-p` / lowering amplitudes `p` of the spin-`n/2` representation of an auxiliary sl(2) acting on the jet index (two-body, double hops: `gamma` raises both indices of site 1 and lowers both of site 2). Single-index specialisation (`l = l' = n`): `alpha = (l1-k) k'`, `beta = -(l1-k') k`, `gamma = eps = 0`.
Order 2 (l1 = 3 data): the exponent pairs `x^a y^b` (a+b = 2) are `(a,b) = (-2,4),(-1,3),(0,2),(1,1),(2,0),(3,-1),(4,-2)` with coefficient matrices of rank `1,2,4,5,4,2,1`; the extreme ones are `C^{(4,-2)}_{pq} = (1/2) (n-k)_2 (n-l)_2 k'_2 l'_2/(n(n-1))^2`, `C^{(-2,4)}_{pq} = (1/2) k_2 l_2 (n-k')_2 (n-l')_2/(n(n-1))^2` (falling factorials `(z)_2 = z(z-1)`; exact on all entries for l1 = 2,3). **General m (verified exactly on all 225 entries of [4,0,0] for m = 1,2,3,4): the coefficient of `x^{2m} y^{-m}` in `R_pq` is `(1/m!) (n-k)_m (n-l)_m (k')_m (l')_m / ((n)_m)^2 = (m!)^3 C(n-k,m) C(n-l,m) C(k',m) C(l',m)/C(n,m)^2`** (falling factorials; mirror image for `y^{2m} x^{-m}` with the roles of p and q exchanged and a sign `(-1)^m`).
Further `Theta_1` facts: rank exactly 2, `theta(a)` in its column space and `theta(a')` in its row space (l1 = 2, 3; for l1 = 1 both in span{theta(a),theta(a')}).
* **l1 = 1 (exact R-matrix form).** `Theta^{[1,l2,0]} = (1+z) theta(a) x theta(a') - z theta(a') x theta(a)`, `z = t/(omega-omega') = a a' t/(a'-a)`, i.e. `(omega - omega' + t) theta(a) theta(a')^T/(omega-omega') - t theta(a') theta(a)^T/(omega-omega')`. (a pure permutation-exchange of the two coherent states with the rational R-matrix `R(u) = u + t P`.)
* **Symbol of the single-index family.** Let `A = 1 - a xi`, `B = 1 - a' eta`, `Thetahat(xi,eta) = sum_{kk'} Theta_{kk'} xi^k eta^{k'}` (the pairing with the jets of the "exponent-0 local solution"):
  `Thetahat^{[l,l,0]} = A^l B^l + t l^2 a a' (xi-eta) A^{l-1}B^{l-1} + ...`; every `t^n`, `n >= 1`, coefficient is divisible by `(xi - eta)` (so on `xi = eta` -- plane waves, e.g. the Q `P3 = e^{z3 u}` -- there is **no coupling**), the top coefficient is `t^{2l-1}: (-1)^{l+1} 4^{l-1} (2l-3)!! (2l-1)!! (a a')^l (xi - eta) = (-1)^{l+1} (2l-2)! (2l)! / (2 (l-1)! l!) (a a')^l (xi - eta)` (values 1, 12, 720, 100800, 25401600, 10059033600 for l = 1..6; l <= 4 from the fitted kernels, l = 5, 6 from the closed form C3). In the variables `Xi = omega - xi`, `H = omega' - eta`, `nu = omega - omega'`, `delta = xi - eta = nu - Xi + H`, `Khat = (omega omega')^l Thetahat` is omega-free:
  `K1 = Xi H + t delta`, `K2 = Xi^2H^2 + 4 t delta Xi H + 2 t^2 delta (3 delta - 2 nu) - 12 t^3 delta`, `K3 = Xi^3H^3 + 9 t delta Xi^2H^2 + 18 t^2 delta Xi H (2 delta - nu) + ...`.
* **Singular locus (single-index family).** With `nu = omega - omega' = (a' - a)/(a a')`, `T = l(l+1)/2`, `d_l = ceil(l/2)`: `det Theta^{[l,l,0]} = c_l (a a')^{s_l} t^{T} nu^{T - d_l} p_l(nu,t)`, `p_l` homogeneous of degree `d_l` in `(nu,t)`:
  `p_1 = nu + t`, `p_2 = nu + 3t`, `p_3 = nu^2 + 6 t nu + 3 t^2`, `p_4 = nu^2 + 10 t nu + 15 t^2`, `p_5 = nu^3 + 15 t nu^2 + 45 t^2 nu + 15 t^3` (l = 5 from the closed form C3; l <= 4 from the fitted kernels: e.g. `l=2: 8 (aa')^3 t^3 (a-a')^2 (3 a a' t - a + a')`, `l=3: 972 (aa')^6 t^6 (a-a')^4 (3 a^2 a'^2 t^2 - 6 a^2 a' t + a^2 + 6 a a'^2 t - 2 a a' + a'^2)`, `l=4: 2654208 (aa')^10 t^10 (a-a')^8 (15 a^2a'^2 t^2 - 10 a^2 a' t + a^2 + 10 a a'^2 t - 2 a a' + a'^2)`). `Theta` is rank deficient at a generic point for `[1,0,0]`, `[3,0,0]`, `[3,1,0]`, `[4,0,0]` (rank 2/3, 8/10, 8/9, 12/15: a pattern of lost directions, not understood); all other computed reps have full rank. The single-index kernel loses rank at `nu = 0` (omega = omega', where CG-completeness fails) and on the roots of `p_l`.

**(c) Top-state sum rule.** Exact identity: **fully symbolic in `(z2,z3,theta1,theta2,a,a')` for [1,0,0], [1,1,0], [2,2,0] and the adjoint [2,1,0]** (sympy), exact rational arithmetic at 3-4 random rational points for all other reps (incl. [4,0,0], [4,4,0]); one-site version symbolic for 10 reps: see headline 5. The two-site sums use the exact top-state Q's: `P3 = e^{z3 u}`, `P2 = e^{z2 u} p2`, `p2 = sum_{k>=0} (z3-z2)^{-k} (Pi^{l2})^{(k)}`, `Pi = (u-th1)(u-th2)`.
It is a strong constraint: at `a = 0` it gives `R(0) = 1`, and every power of `a` gives one linear relation between the top-state jets and `Theta`; every candidate kernel must satisfy it.

**(d) Determinantal tests.**
* `theta(a)` is exactly a decomposable bivector (Plucker vector of a 2-plane) in the Taylor-normalised index space `V = C^{l1+2}` (the space of jets of order `<= l1+1`): `theta^ = a (xi1 - xi2) [(1-a xi1)(1-a xi2)]^{l1}`, i.e. the image of the lowest-weight vector `(xi1 - xi2)` of the spin-`l1` irrep of sl(2) in `Lambda^2 V` (V = spin `(l1+1)/2`) under `exp(-a S_+)`.
* Plucker-relation tests on `Theta_n` viewed as a matrix on `Lambda^2 V x Lambda^2 V` ([2,0,0], [3,0,0]): all columns and rows of `Theta_0` and `Theta_1` are decomposable (Plucker relations exactly satisfied); `Theta_2, Theta_3,...` and the full `Theta(t)` are NOT (so `Theta` is **not** a compound matrix `Lambda^2 M`; for [2,0,0] top order t^4 is again decomposable). The `[l,l,0]` kernel is not of compound type either (rank 1 at t = 0 would need `rank M = 2`).
* "Box-exchange" hypothesis tested and refuted: for `l1 = 1` the kernel is a 2-term sum of permuted-parameter products (above, exact); for `l1 = 2` (single-index kernel against the 6 products `prod_{r in S}(u_r - xi) prod_{s notin S}(u_s - eta)` with parameter sets `u = (omega + c t, omega' + c' t)`, shifts in `t/2` units up to +-3) the kernel is NOT in the span: the fused-R picture is not literally a sum over parameter permutations.
* Rank count of `Theta_n`: `1, 2, 4, 8` for `n <= 3` (l1 >= 3; capped by the number of pairs for smaller reps), then `13` at n = 4 for [4,0,0] (not 16): the "two channels per order" doubling stops at n = 3. `Theta_1 = sum_{r=1,2} u_r (x) v_r` has two products of decomposable bivectors (all columns/rows of `Theta_1` are decomposable).

## K4. The conjecture (general Lambda = [l1,l2,0]) and its status

Let `n = l1+1`, pairs `(kl)` with `0 <= k < l <= n`, `l >= l2+1`; `a = 1/omega`, `a' = 1/omega'`, `t = 1/theta12`, `x = a t`, `y = a' t`.

**C0 (one-site, theorem-level, verified 17 reps).** `theta_kl(a) = c_kl a^{k+l-l2-1}` with the closed form above; equivalently `sum theta_kl m_kl = eps a^{-l2} P3 W[(1-a d)^{l1}P2,(1-a d)^{l1}P3]`.

**C1 (block theorem, verified).** `Theta^{[l1,l2,0]}_{ij} = (l2!(l1-l2)!/l1!)^2 (a a')^{-l2} Theta^{[l1,0,0]}_{ij}` for non-forced `i,j`. Equivalently `Theta_ij = theta_i(a) theta_j(a') R^{l1}_{ij}(a t, a' t)`, with a **universal, l2-independent, Laurent polynomial `R^{l1}_{(kl),(k'l')}(x,y)`**, homogeneous of degree 0 in `(a,a',theta12)` (a function of `x, y` only). All of the content of l2 is in `theta`.

**C2 (structure).** `Theta` is a polynomial in `(a,a',t)` with rational coefficients, `0 <= p,q <= 2 l1 - l2`, `p + q - n = e_i + e_j` for each monomial `a^p a'^q t^n` of entry (i,j); `R^{l1}(x,y)` has `R(0,0) = 1`, `R_ij(x,y) = Laurent`, exchange `R_ij(x,y) = R_ji(-y,-x)`; `Theta(t=0) = theta x theta'`; `Theta_1` has rank 2; **every `t^n` (n >= 1) term of the symbol vanishes at `xi = eta`**; top-state sum rule C5; t-degree `Dt` as in the table.

**C3 (closed form for `[l,l,0]`, verified l <= 4 exactly, l = 5 predicted and confirmed).** (Headline 4.) In terms of the sl(2) language: `Thetahat = sum_{K,K'} C(l,K) C(l,K') { E[(u+vs)^K (u'+v's')^{K'}] - K K' v v' E[(u+vs2)^{K-1}(u'+v's2')^{K'-1}] }`, `u = -a xi`, `v = a t`, `u' = -a' eta`, `v' = -a' t`,
with the expectations over independent Gamma variables `s ~ Gamma(K')`, `s' ~ Gamma(K)`, `s2 ~ Gamma(K'+1)`, `s2' ~ Gamma(K+1)` (a Laplace/Borel representation: `(k')_m x^m` are the moments); `Thetahat^{[l,l,0]}|_{t=0} = A^l B^l`.

**C4 (first order, closed form, verified in all reps).** The `theta1_closed` formula of K3(b) (alpha, beta of rank 2; gamma, eps of rank 1).

**C5 (sum rule, verified).** `sum_ij Theta_ij m^top_i(th1) m^top_j(th2) = M1 M2 R(a) R(a')`, `R(a) = (1-z2 a)^{l1-l2}(1-z3 a)^{l1}`.

**Verification status.** C0: 17 reps exact (K1). C1: exact for l1 <= 3 (11 kernels, all entries), direct fits for l1 = 4 (l2 = 1,2,3 vs l2 = 0, deviation <= 3e-58 at 2 theta12 x 2 (w,w')). C2: all reps. C3: l = 1..4 exact polynomial identity with the fitted kernels; **l = 5 was a genuine prediction and is confirmed by direct fits of [5,5,0] (441 states, 36 unknowns, rank 36, residual 0.0): max relative deviation 1.5e-59 .. 9.2e-58 at theta12 = 31/35, 29/45 and (w1,w1') = (5/2,7/4), (-2/5,3), (13/5,-7/2)**. C4: all computed reps (0 mismatches, incl. held-out [4,4,0]). C5: 10 reps (two-site) + 10 (one-site).

**Predictions** (next reps):
* `[4,2,0]` (12 pairs: `(03),(04),(05),(13),(14),(15),(23),(24),(25),(34),(35),(45)`): `theta` from C0 (listed above), `Theta^{[4,2,0]} = (1/6)^2 (a a')^{-2} Theta^{[4,0,0]}|` (N = (2!2!/4!)^2 = 1/36), a-degree 6, t-degree 11, full rank 12 at a generic point, ranks of `Theta_n` `1,2,4,8,12,12,12,11,9,6,4,2`; **verified directly** at two theta12 and two (w,w'): deviation <= 8e-62 (BLK output); its corner block `(k,5)` is `(1/36)(a a')^2 Theta^{[4,4,0]}` (C3). Same block-theorem data for the others: [4,1,0] (14 pairs, t-degree 12, rank 12 (deficient), ranks `1,2,4,8,13,14,14,14,13,10,8,6,4`), [4,3,0] (9 pairs, t-degree 9, rank 9, ranks `1,2,4,7,9,9,8,6,4,2`).
* `[5,1,0]` (20 pairs, e.g. `theta` from C0, N = (1!4!/5!)^2 = 1/25, a-degree 9): `Theta^{[5,1,0]} = (1/25)(a a')^{-1} Theta^{[5,0,0]}|`; t-degree 16; the first-order structure C4 with n = 6; the sum rule C5 with `R(a) = (1-z2 a)^4 (1-z3 a)^5`; `Theta(t=0)` from C0. `[5,5,0]` (36 entries from C3: **confirmed**, see above); `[6,6,0]` from C3 (t-degree 11, top coefficient 10059033600 (a a')^6 (xi - eta)).

## Surprises / corrections to the plan

1. The top-state sum rule of the plan (`= M^top M^top`) is **wrong by the factor `R(a)R(a')`**; companion-frame top state is not `e1 x e1` (`P_phys^{-1} e_hw` is not parallel to `e_hw`); corrected statement exact.
2. The S2 double pole is the scalar `c2 = l2(l1+1)` also for the l2 = 0 reps (value 0): no loss of structure for one-row reps.
3. Not 4 from D2 (convention only), and `Theta` is rational in `omega` ONLY through the single `omega^{2l1-l2}` normalisation; the kernel needs no `omega - omega'` poles (the XXX R-matrix poles `(omega - omega' + m t)` appear only in determinants of Theta, not in entries).
4. The jet index space (jets of order `<= l1+1`, dimension `l1+2`) carries an sl(2) of spin `(l1+1)/2` acting on symbols by the Moebius action `(U_s f)(xi1,xi2) = prod_i (1 - s xi_i)^{l1+1} f(xi_i/(1 - s xi_i))`; the one-site kernel is `theta^(a+s) = U_s theta^(a)` (exact, l1 <= 4 symbolic): a coherent state of the spin-`l1` irrep of `Lambda^2 V` generated from `(xi2 - xi1)`. The two-site kernel is NOT covariant under the diagonal shift `a -> a+s, a' -> a'+s` (tested, l1 = 1, 2).
5. The number of t-orders is only linear in l1 (not quadratic as the XXX fusion count `l1^2` of R-matrix factors would suggest): `Dt = 4(l1-1) - max(0, 2 l2 - 3)` (l1 >= 2).

## Not done / limits

* A closed form for the **general** (four-index) kernel `R^{l1}_{(kl),(k'l')}` beyond first order: established are the single-index family (C3, all orders), the first-order term (C4, all pairs), the extreme second-order terms (m-fold double hops) and the universality (C1). A box-permutation R-matrix form and a compound-matrix form were tested and refuted (K3d).
* The full (a',t) polynomial grid for [4,1,0], [4,2,0], [4,3,0] (replaced by the block theorem against [4,0,0]); l1 = 5 direct fits except the single-index [5,5,0] (a corner-block check of [5,1,0] was started but not finished: the 1225-state build was too slow and I stopped it).
* One-site kernel for weights with multiplicity > 1 (their eigenvectors are not GT vectors; they would need irrational eigen-decompositions); not needed since the exactly diagonal states give rank-complete data.
