# RESULTS: the XXX -> Gaudin limit (agent P2, 2026-10-16)

Files: `xxx_chain.py` (XXX chain port: exact sparse 64x64 operators, truncated h-series `Ser`), `xxx_limit.py` (X1-X4 + checks; `python -X utf8 xxx_limit.py`, ~2 min, 0 failures; X4 is the slow part, ~95 s of sympy).
Numbers: zeta = (3/2,5/4,2/3) (s = sum zeta = 41/12), omega = (5/2,7/3,3/5), theta = (2/7,-3/5) (th12 = 31/35), Lam = (2,1,0), h in {1e-2,1e-3,1e-4}.

**Exactness policy.** X2/X3 coefficients are exact rationals: the whole XXX chain is evaluated with h the formal series variable (series to h^15, Fraction coefficients) and z = exp(h zeta), w = exp(h omega) as exact series, so there is no Richardson error. The finite-h cross-checks use rational h with the exponential series truncated at h^15. X1 uses mpmath, 80 digits (250 at h = 1e-30), exact exponentials. X4 is sympy-symbolic.

**Basis note.** `gl3_reps.py` lists GT patterns in the opposite order from the notebooks (its index 0 is the LOWEST-weight vector, weight (0,1,2)). `xxx_chain.py` re-orders to the notebook order (patterns sorted descending), so basis vector 0 = e1 = HW ket, annihilated by E12, E13, E23 (checked), weights of the 8 states = (2,1,0),(2,0,1),(1,2,0),(1,1,1),(1,1,1),(1,0,2),(0,2,1),(0,1,2). E13, E31 come from commutators exactly as in the notebook `JJ`. All vectors below are in this order.

## Port checks (h = 1/2 notebook point)
[E_ij,E_kl] relations: 0 violations of 81. qdet t3 = chi3 nu1(u)nu2(u-h)nu3(u-2h): exact, True (one exact 64x64 quantum minor: 0.014 s - sparse Fraction operators are fast, no caching needed). [t1,t2] = 0 True. Omega_zw(x)Omega_zw is an exact left eigencovector of t_a(u=2/9), a = 1,2,3, with bra twist G' = P_aux G_c(w) P_aux^-1: [True,True,True]; control diag(w): [False,False,True].

## X1 frames and vacuum
**(a)** r_k(zeta,omega) = (4/5, -7/11, 65/44). |r_k(z,w) - r_k(zeta,omega)| at h = 1e-2/1e-3/1e-4: (9.99e-4, 2.65e-3, 4.30e-3) / (9.999e-5, 2.651e-4, 4.308e-4) / (9.99999e-6, 2.65146e-5, 4.3086e-5); successive ratios 9.99,9.98,9.98 -> exactly O(h). Exact series (Fractions): r1 = 4/5 + h/10 - h^2/80 - h^3/480 + ..., r2 = -7/11 - (35/132)h + (175/3168)h^2 + ..., r3 = 65/44 + (455/1056)h - (455/4224)h^2 + .... So r_k(z,w) = r_k(zeta,omega) + O(h), with O(h) coefficients (1/10, -35/132, 455/1056).

**(b)** The notebook's P_aux (exponential product) equals the rows (1, z_i - chi1, z_i^2 - chi1 z_i + chi2) to 2e-81. (G'-1)/h with G' = P_aux(z) G_c(w) P_aux(z)^-1 CONVERGES, error O(h), to

g' = P_aux(zeta) g_c(omega) P_aux(zeta)^-1 = [[-21/10, 169/28, -44/105], [-18/5, 51/7, -44/105], [-18/5, 169/28, 26/105]]

(no transpose, no inverse). max-entry distance: h=1e-2: 0.1458, 1e-3: 0.01443, 1e-4: 0.001441, 1e-30: 1.4e-29 (rate 10x/decade). Wrong variants stay O(1): (g')^T 9.6, (g')^-1 8.6, P^-1 g_c P 17.9, -g' 14.6. The Jordan block of G_c(1,1,1) is harmless: the conjugation by P_aux(z) (which degenerates like diag(1,h,h^2)) restores G' = 1 + h g' + O(h^2).
Gaudin side, exact (u = 2/9, 64x64): Omega_{zeta omega}(x)Omega_{zeta omega} is a left eigencovector of S2(u), S3(u) for the plan's cdet (Lax E_ij, column-ordered) with g = g' (and also g'^-1, which has the same flag), NOT with (g')^T, P^-1 g_c P, or the control diag(omega). (The mirror convention Lax E_ji + row-ordered cdet works with (g')^T; the other two mixed conventions fail.)

**(c)** x0(z) = e1^T P_phys(z)^-1 (one site; the two-site x0 is x(x)x). ||x0(e^{h zeta})|| ~ h^(1.0000) (two-site h^2). Leading direction (projective, component 5 = 1): (1/10, 1/3, 1/10, 17/30, 3/10, 1, 7/15, 7/10) (components in the basis order above). The Gaudin x0^G(zeta) built by the same Gauss-factor construction has the identical direction: projective distance 8.8e-32. **PARALLEL** (the plan expected "not in general"; here they are parallel). Correction is O(h) (direction stable to 8.8e-22 between h=1e-20 and 1e-30).

## X2 Hamiltonian dictionary (exact, u in {3/7,-2/9,5/3,-7/4}, plus zeta-sets (1/3,-2/5,7/4) and (2,-1/2,3/7))
With t1/nu = 3 + hA1 + h^2 A2 + h^3 A3 + ..., t2/(nu(u)nu(u-h)) = 3 + hB1 + h^2 B2 + h^3 B3 + ..., nu(u) = (u-th1)(u-th2), G = diag(e^{h zeta}) (A0 = B0 = 3 exactly, scalars):

- A1 = S1 and B1 = 2 S1, scalars, S1(u) = s + 3 sum_s 1/(u-th_s) = tr L_g (s = sum zeta).
- **S2 = -A2 + (1/2)S1^2 + (1/2)S1'**
- **S3 = A3 - B3 + (1/6)S1^3 + (1/2)S1 S1' + (1/2)S1''**
- (B2 = A2 + (1/2)S1^2 + (1/2)S1'; no relation fixes A2, A3-B3 beyond the above.)
The c-numbers c2 = S1^2/2 + S1'/2 and c3 = S1^3/6 + S1S1'/2 + S1''/2 are exactly these rational functions of u and depend on zeta only through s (verified exactly at 4 u's for zeta = (3/2,5/4,2/3) and at 2 u's for two more zeta; the free exact fit with {A2,1,A1^2,B2} resp. {A3,B3,1,A2,B2} gives coefficients (-1, c2, 0, 0) resp. (1,-1,c3,0,0)). With sigma_k = sum_s (u-th_s)^-k: c2 = (s+3 sigma_1)^2/2 - (3/2) sigma_2, c3 = (s+3 sigma_1)^3/6 - (3/2)(s+3 sigma_1) sigma_2 + 3 sigma_3. For our numbers:
S1 = 41/12 + 21/(7u-2) + 15/(5u+3);
c2 = 1681/288 + 17717/(124(7u-2)) + 147/(7u-2)^2 + 55/(124(5u+3)) + 75/(5u+3)^2;
c3 = 68921/10368 + 1088017/(2976(7u-2)) + 124019/(124(7u-2)^2) + 1029/(7u-2)^3 - 256045/(2976(5u+3)) + 275/(124(5u+3)^2) + 375/(5u+3)^3.
(The plan's guessed structure S2 = S1^2/2 - A2 + c(u) holds with c = S1'/2.)
Finite-h check at u = 3/7: max-entry residual of (t1/nu-3)/h vs A1+hA2+h^2A3: 4.80e-5 / 4.78e-8 / 4.78e-11 at h = 1e-2/1e-3/1e-4 (ratio 1e3 = O(h^3)); same for the B-combination (1.7e-2/1.6e-5/1.6e-8). u -> infinity (u = 1e9): S2 -> e2(zeta) = 89/24, S3 -> e3(zeta) = 5/4 (1.5e-8, 9.5e-9).
**Which cdet convention the limit selects:** the plan's (Lax L_g = g + sum E_ij/(u-th), column-ordered cdet, d to the right): S1 scalar, [S_k(u),S_l(v)] = 0 exactly, and the dictionary closes. Row-ordered with E_ij and column-ordered with E_ji (the mixed variants) do NOT commute and do not fit; row-ordered + E_ji (the transposed mirror) commutes and fits identically. Nothing in the plan needs changing.

## X3 localised covectors (one site, exact series; <v_lambda| = h^|lambda| (v0 + h v1 + h^2 v2 + ...))
Omega_{zeta omega} = [1, 65/44, 4/5, 20/11, 18/11, 585/242, 56/55, 126/121]. Notation: C = sum_k zeta_k E_kk; S3[-k] = coefficient of (u-th1)^-k in the one-site Gaudin S3(u) = cdet(d - zeta - E/(u-th1)) (constant S3[0] = 5/4 = e3, S3[-3] = 0 as an operator; S2[-2] = 3).
**T(1,1)(th1-2h) = t2:** valuation h^2.
v0 = 3 Omega_{zeta omega} (exact).
v1 = [23/4, 3445/352, 11/2, 599/44, 663/44, 48555/1936, 91/11, 2709/242]
v2 = [529/96, 14495/1408, 1499/240, 17111/1056, 7941/352, 161915/3872, 13993/1320, 2324/121]
v1 = Omega . ( 215/8 - (9/2) C - S3[-2] ) EXACTLY (LS residual 0; independent generators needed: 1, C, S3[-2]; S3[-1], S3[-3] and Y = sum (g-g')_ij E_ji are not needed / dependent; equivalently v1 = Omega.(-49/8 - (9/2) Y - S3[-2])). v2 is NOT in Omega.span{1,C,S3[-k],Y} (residual^2 = 1.38) but IS in Omega.{degree <= 2 products} (rank 8 of 8 - uninformative).
**T(2,1)(th1-2h) = t2(v)t1(v+h) - t3(v+h):** valuation h^3.
v0 = -3 Omega_{zeta omega}; v1 = [-31/4, -5395/352, -71/10, -230/11, -207/11, -66105/1936, -707/55, -3843/242]; v2 = [-961/96, -36595/1408, -2507/240, -20725/528, -6263/176, -38155/484, -34909/1320, -4781/121].
v1 = Omega . ( -215/8 + (9/2) C ) EXACTLY (no S3 term; residual 0). (v2 again needs degree-2 products: residual^2 = 9.2 with linear generators.)
Finite-h check (rational h): |<Omega_zw|T/h^|lambda| - (v0+hv1+h^2v2+h^3v3)| = 3.4e-7, 3.4e-11, 3.4e-15 (T11) and 1.3e-6, 1.3e-10, 1.3e-14 (T21) at h = 1e-2,1e-3,1e-4: O(h^4) as expected.
**Spectator factors:** f(1,1) = (th12-h)(th12-3h) = 961/1225 - (124/35)h + 3h^2 -> th12^2 (relative first-order coefficient -4/th12 = -140/31); f(2,1) = (th12-h)^2(th12-3h) = 29791/42875 - (961/245)h + (31/5)h^2 - 3h^3 -> th12^3 (relative -5/th12 = -175/31).
**Site 2** (one-site chain at th2): <Omega|t2(th2+h): valuation h^2, w0 = 3 Omega, w1 = [33/4, 3835/352, 69/10, 599/44, 663/44, 40365/1936, 413/55, 189/22] = Omega.(133/8 - (3/2)C - S3^{th2}[-2]) exactly (S3 coefficient at th2). <Omega|t2(th2+h)t1(th2-h): valuation h^4, w0 = [5/2, 195/176, 7/5, 0, 0, -4095/968, -42/55, -315/121] = 3 Omega.(C - s) = Omega.(-41/4 + 3C), NOT parallel to Omega (components 3,4 vanish); w1 = Omega.(-943/24 + (23/2)C). The effective covector of the second site-2 identity, -f(2,1)<Omega|t2t1 + 6h^4 chi3 f(1,1)<Omega|, also has valuation h^4 with w0 = Omega.(2028671/171500 - (89373/42875) C), not parallel to Omega either.
**Two-site sanity (64x64 exact, h = 1/1000, truncated series):** residuals of all four identities <Omega(x)Omega|T(1,1)(th1-2h) = f11 <v11|(x)<Omega|, same for T(2,1) with f21, site-2 t2(th2+h) = f11 Omega(x)(Omega t2^{(1)}), and the t2t1 identity with -f21 w + 6h^4 chi3 f11 Omega: all exactly 0 (sizes of the left sides 1.4e-5, 1.2e-8, 1.4e-5, 3.5e-11). The same identity also holds as a formal h-series on all 64 components to order h^11.
**X1 operator** X1 := lim (T(1,1)(th1-2h)/(h^2 th12^2) - c)/h, c = 3 (the h^0 and h^1 coefficients vanish; the h^2 coefficient is exactly the scalar 3 th12^2): 64x64, 216 nonzero entries. **X1 = -(4631/124) 1 + 3 H1 - S3[th1,-2]** EXACTLY, in the linear span of {1, W, H1, S3 Laurent coefficients at th1, th2} (basis rank 7 of 9; H1 = residue of S2(u) at th1, S3[th1,-2] = coefficient of (u-th1)^-2 of S3, W = sum zeta_k (E_tot)_kk, H1+H2 = (41/2) 1 - W). The quadratic extension (73 operators, rank 24) finds nothing new. Moreover Omega0(x)Omega0 . X1 (Omega0 = Omega_{zeta omega}) has 8x8 rank 1: = u (x) Omega0 with u = Omega0 . ( 431/62 - 3C - S3^{1-site}[-2] ) exactly (u = [-967/124, -62855/5456, -874/155, -18131/1364, -12477/1364, -384345/30008, -12236/1705, -41391/7502]).
Prediction for W1: the site-1 localiser in the Gaudin D2 is 3 H1 - S3[th1,-2] (mod 1); the localised one-site covector is in Omega0.span{1, C_zeta, S3^{1-site}[-2]} (no (g-g')-type element needed). The constant -4631/124 depends on th12, zeta.

## X4 Baxter -> oper (sympy, symbolic h; L = 2 adjoint; direct and dual identical constraints)
Unknowns t1 = nu[3 + hA1 + h^2A2 + h^3A3], t2 = nu nu^[-1][3 + hB1 + ...], chi3 = e^{hs}, Q(x+kh) and shifted A_k, B_k Taylor-expanded, nu(x) = (x-a)(x-b).
- h^0: vanishes identically (no condition).
- h^1: only Q appears, coefficient -nu^2[(s+A1-B1) nu + 3 nu']  => **B1 = A1 + s + 3 nu'/nu**.
- h^2: Q' coefficient ~ nu^2[(A1 - s) nu - 3 nu'] => **A1 = s + 3 nu'/nu = S1**, hence **B1 = 2 A1 = 2 S1** (no free A1); Q'', Q''' vanish; Q coefficient vanishes iff **B2 = A2 + S1^2/2 + S1'/2**. No condition on A2, A3, B3.
- h^3: ODE  [-nu^3 (direct), +nu^3 (dual)] ( Q''' + p2 Q'' + p1 Q' + p0 Q ) = 0 with
 direct: p2 = s - 3 nu'/nu ; dual: p2 = -(s + 3 nu'/nu) = -S1 ;
 p1(dual) = -A2 + S1^2/2 + S1'/2 ; p0(dual) = B3 - A3 - S1^3/6 - S1S1'/2 - S1''/2 (verified symbolically; the direct ones are given in the script).
 With the X2 dictionary: **dual**: (d^3 - S1 d^2 + S2 d - S3) Q = 0 EXACTLY (coefficients differ by 0).
 **direct**: with Q = nu^2 R the equation for R is EXACTLY the formal adjoint -(d^3 + d^2 S1 + d S2 + S3) = (-d)^3 - (-d)^2 S1 + (-d) S2 - S3 acting on R (differences 0 for Q''', Q'', Q', Q); the gauge powers nu^0, nu^1 fail.
**One sentence:** At order h^3 the dual Baxter equation on Q^j = z_j^{+x/h} q^j becomes exactly the Talalaev oper (d^3 - S1 d^2 + S2 d - S3) Q^j = 0 with solutions e^{+zeta_j x} q^j, while the direct Baxter equation on Q_i = z_i^{-x/h} q_i becomes exactly the formal adjoint oper (-d)^3 - (-d)^2 S1 + (-d) S2 - S3 acting not on Q_i but on R_i = Q_i/nu(x)^2 with nu = prod_s (x-th_s) (so Q_i = nu^2 R_i, solutions e^{-zeta_i x}, the 1/nu^2 poles being absorbed in the polynomial q_i). The exponent signs follow from S_k(u -> infinity) = e_k(zeta) (checked, 1e-8 at u = 1e9): oper kills e^{+zeta_j x}, adjoint kills e^{-zeta_i x}.
(Consistency, inferred not separately tested: local exponents of the oper at th_s are {0,2,4}; those of the adjoint {2,0,-2}, i.e. R has double poles, which nu^2 clears.)

## Surprises
1. x0 one-site direction IS parallel to the Gaudin x0^G(zeta) (scale h^1 per site); plan expected non-parallel.
2. A1 and B1 are scalars with A1 = S1 exactly, B1 = 2S1 - forced by the h^1, h^2 orders of Baxter itself; the (c-number) dictionary has only S1-derivatives in it.
3. The direct Baxter equation is the adjoint oper up to the gauge Q = nu^2 R (not the adjoint oper on Q_i). P1 should test the adjoint oper on Q_i/nu^2, or equivalently nu^2 (adjoint) nu^-2 on Q_i.
4. X1 is a 3-term combination (3 H1 - S3[th1,-2] + const) and Omega0(x)Omega0.X1 is exactly localised on site 1; v1 of T(1,1) needs S3[-2] but v1 of T(2,1) does not; the site-2 (2,1) covector has leading term 3 Omega.(C - s), not parallel to Omega.
5. G_c Jordan block: no obstruction to (G'-1)/h -> g'.

## Not done / limits
- Richardson over the three h's was replaced by the exact series (stronger); finite-h residuals are reported instead.
- v2 was only tested for linear-span membership (fails) and degree-2 products (uninformative: spans all 8 dims); no closed form for v2 or for the constants (215/8, 431/62, -4631/124 ...) as functions of (th12, zeta, omega) - only at the single numerical point.
- Only Lambda = (2,1,0), L = 2, one theta, one (zeta, omega); X4's gauge factor nu^2 = prod(x-th_s)^(lam1-lam3) is established for this rep only. [3,1,0] not attempted. No Wolfram used.
