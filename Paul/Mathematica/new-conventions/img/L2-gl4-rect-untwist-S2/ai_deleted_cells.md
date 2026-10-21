## Deleted Cell 26 — 2026-10-10 11:06:30
**Kind:** markdown
```markdown
### gl(4) [S,S,0,0] untwist result (2026-10-10), SS = 1

- Vacuum: Omega_zw = e1 exp(X) with X in the abelian nilradical (exact); closed form x_ij = (-1)^i (z_i' - z_j)/(z1-z2) prod_{k=1,2}(w_k-z_i)/(w_k-z_j); w3, w4 drop out, w1, w2 torus directions; torus invariant = z cross-ratio (so the old kappa-torus family has det X' = 0 and never sees the singlet: (0,0,1,1) slice exactly 0 at finite twist, cell 8).
- Paper's family W_om = Omega(1 - om rho(z))/(1 - om z1z2) = e1 exp(R), R_ij = x_ij(1 - om z_i'z_j)/(1 - om z1z2) (exact); v2 = W_om, v1 = W_om rho(wk - z); undo exact: phi_u = delta^Lambda e1 exp(R'), R'_ij = R_ij (delta_j/delta_i - 1).
- Untwisting regime (z, w, om, wk all 1 + eps(.)): all six weights of the undone covector O(1), R'(0) generic 2x2 block (det = 5.6e8/2.5e10 != 0 at the primary point, nonzero at all 12 grid points). Control with om fixed at 3/7: singlet slice O(eps^2), kappa -> 0.9838 (K/S = 0.01617 constant in eps).
- Chain side (eps = 1,2,4e-6, 12 grid points): S ~ eps^2 (delta^Lambda), K ~ eps^3 for [1,1,1,1] and [2,1,1,0], K = 0 for [2,2,0,0]; max |K/S| = 1.9e-7 / 6.1e-7 at eps = 1e-6, doubling with eps. CG slices read off from R_C alone with relative error 1e-6..5e-6 (linear in eps), Richardson agreement 2e-12 / 6e-13.
- Localisation: only t1, t2 localise (rank 1), at v = theta1 - h (not theta1), site-2 factor = Omega exactly. The old kappa family is in span{Om2.t1,t2 at six points}; W(x)W in the span of all shapes; the W-family product v1(x)v2 is in none of the spans tried (273 generators, rank 30, +1). All 36 states have nonzero overlaps and distinct eigenvalues, so Q-number-ness is only meaningful for a bounded ansatz.
- gl(3) contrast confirmed symbolically (cell 14).
```

## Deleted Cell 25 — 2026-10-10 11:06:30
**Kind:** markdown
```markdown
**Results, cells 15-17 (Task 1, SS = 1, finite twist z = (3/2,5/4,2/3,4/7), w = (5/2,7/3,3/5,9/8); executed in K1 via evaluateExpression, text below copied from the kernel output).**

Rank scan (rank of <Om(x)Om|T_lambda(v) as a 6x6 matrix; v = theta1 + m h):

| lambda \ m | -2 | -1 | 0 | 1 | 2 |
|---|---|---|---|---|---|
| (1,0) | 2 | **1** | 2 | 2 | 2 |
| (1,1) | 2 | **1** | 2 | 2 | 2 |
| (2,0) | 3 | 0 | 3 | 3 | 3 |
| (2,1) | 4 | 0 | 4 | 4 | 4 |
| (2,2) | 2 | 0 | 0 | 6 | 6 |

Only t1 and t2 localise, at v = theta1 - h (m = -1), not at theta1 as expected from the paper's App. M (convention shift by h: in our Lax L(x) = x + hE the quantum-minor shifts run downward). At m = -1 the site-2 factor is Omega_zw exactly (checked: the column vector of the rank-1 matrix equals Omega), the site-1 factor c1, c2 differs from Omega (and rank{Om,c1,c2} = 3).

Span tests (rank of generators {Om2, Om2.T_lambda(theta1+m h)}, and rank after appending each target; targets: 1 = W(x)W om=3/7, 2 = W(x)W om=-5/3, 3 = v1(x)v2 of the W family (om = 3/7, wk = 7/5), 4 = old kappa-torus family Om rho(delta)(x)Om):

| generators | #gens | rank | rank+target 1,2,3,4 |
|---|---|---|---|
| t1,t2 at m=-1 only | 3 | 3 | 4,4,4,4 |
| t1,t2, m=-3..2 | 13 | 7 | 8,8,8,**7** |
| t1..t4, m=-3..2 | 25 | 9 | 10,10,10,**9** |
| all 34 shapes (<=4 rows, parts <=3), m=-4..3 | 273 | 30 | **30,30**,31,**30** |

So: the old kappa family (Omega rho(delta) (x) Omega) is in the span of Om2.{t1,t2 at six points} (needs more than the single localising point); W(x)W is in the span of all shapes; the W-family product v1(x)v2 (W rho(delta) (x) W) is NOT in the span of any generator set tried (rank +1 each time) - the paper's Q-determinant property for it is therefore not of the linear 'T_lambda acting on <Om(x)Om|' type in my ansatz. No exact c(w2') polynomial was extracted: the old-family solution in the 13-generator set is non-unique (rank 7 of 13).

Task 1(c) (all 36 states at the finite twist, pr = 200): 36 distinct eigenvalues of t1(3/7)+Pi t2(-2/9); #states with zero overlap with Om(x)Om, W(x)W, v1(x)v2: 0, 0, 0 (min |<Om2|Psi>|/||Psi|| = 7.3e-6, max 8.8). Hence with simple spectrum and nonvanishing overlaps every covector's overlap IS a function of the eigenvalues; a 'Q-number' statement is only non-trivial for a bounded ansatz (the span table). Not done: direct eigenvalue-combination check on eigenstates (c) beyond this.
```

## Deleted Cell 21 — 2026-10-10 11:06:30
**Kind:** markdown
```markdown
**Results, cell 14 (Task 5).** (a) gl(3): n'13 - n'12 n'23 = (a c (d1-d2) d3 + b d2 (d3-d1))/(d1 d2); for d_k = 1/(w2-z_k) this is (b(z3-z1) + a c(z1-z2))/(w2-z3) = (d3/d1 - 1)(b - a c z21/z31) (confirmed, spec formula correct). With the gl(3) vacuum n(w1,y): b' - a'c' - a'c'(z2-w2)/(y-z2) = 0 exactly (symbolic), and b' - a'c' = (w1-z1)(z1-z3)(z3-z2)/((w1-z3)(z3-w2)(z3-y)) (carries the factor from y; O(eps) in the gl(3) regime because y cannot untwist). (b) gl(4) abelian: n' = 1 + X', X'_ij = x_ij (d_j/d_i - 1) exactly and (X')^2 = 0 (both True, symbolic): no commutator direction, an entry vanishes only if a twist difference does.

Comparison: the paper's limit (4) untwists kappa (the torus point), omega and z together; in gl(4)-rect there are two independent localising deformations (torus Lambda^2(1-kappa G) and the Lambda^2-polynomial 1 - omega rho(z), which preserves the rotation form because of the special companion-vacuum cross-ratio), the undone block at eps = 0 is a generic 2x2 matrix and everything is entrywise. In gl(3) the E13 coefficient is the commutator b' - a'c' proportional to (z2 - w2)/(y - z2) with y pinned at y0 = (t-1)/(t-2) != 1, hence O(eps). The pure kappa-torus family has the finite-twist blind spot det X' = 0 (cell 8) - the gl(4) shadow of 'b' = 0 iff omega = y'.
```

## Deleted Cell 19 — 2026-10-10 11:06:30
**Kind:** markdown
```markdown
**Results, cells 11-13 (Tasks 3, 4; SS = 1; eps in {1,2,4}e-6; phi primary, omw = (3/4,-1/5,2/9,-7/4), (o,k) grid 4 x 3; pr = 200).** Timings: 3 eigen-decompositions + 12 points x 3 eps loop 0.58 s, control 0.13 s. Branches at eps = 1e-6: 9 states, 3 HW branches ([1,1,1,1], [2,1,1,0], [2,2,0,0]); max rr of branches 9.9e-8 / 2.0e-7 / 0, min rr of the other states 2.97 / 1.41 (clean separation); hwFrac = 1 - 3e-16, 1 - 1e-14, 1. Top branch = e1(x)e1 to 10^-198.

JOINT regime (o untwists with z):

| branch | S eps-exponent | K eps-exponent | max_grid |K/S| at eps = 1,2,4e-6 | undo rel dev |
|---|---|---|---|---|
| [1,1,1,1] | 2.0000 (1.99999..2.00001) | 3.0000 | 1.89e-7, 3.78e-7, 7.55e-7 | 0 (acc 194) |
| [2,1,1,0] | 2.0000 | 3.0000 | 6.10e-7, 1.22e-6, 2.44e-6 | 0 (acc 198) |
| [2,2,0,0] | 2.0000 | K = 0 exactly | 0 | 0 |

(S carries delta^Lambda = O(eps^2); the CG slice itself is O(1).) So K/S = O(eps) for every branch including the singlet, kappa -> 1 linearly. CONTROL (om = 3/7 fixed, not untwisting): singlet S ~ eps^4.0000, K ~ eps^4.0000, K/S = 0.016171, 0.016171, 0.016172 at the three eps (kappa -> 0.984, NOT 1); the other branches unchanged (K/S = 6.1e-7...2.4e-6).

Task 4 (CG from R_C alone): estimator rhs_j = R_C/delta^Lambda against phi0_slice (exact eps^0 coefficient): grid spread of C over the 12 (o,k) points / rel err vs chain-side slice <m(x)e1|P Psi>, at eps = 1,2,4e-6 (both linear in eps): singlet: spread 3.6e-6, 7.2e-6, 1.4e-5; err 1.18e-6, 2.37e-6, 4.73e-6. [2,1,1,0]: spread 5.9e-7, 1.2e-6, 2.3e-6; err 6.7e-7, 1.3e-6, 2.7e-6. Top: 0. Richardson (2 C(eps) - C(2 eps)) vs chain-side Richardson: rel diff 2.1e-12 (singlet), 5.9e-13 ([2,1,1,0]), 0 (top) - O(eps^2). All slices are 1-dim at SS = 1 (no 2x2 system needed).
```

## Deleted Cell 15 — 2026-10-10 11:06:30
**Kind:** markdown
```markdown
**Results, cells 9-10 (Task 2; SS = 1).** Single point (phi primary, o = -12/5, k = -1/5): symbolic-eps series of phi_u/delta^Lambda: every one of the six weights is O(1), leading coefficients (state weight: coefficient) (0,0,1,1): 2258375000/25411092993; (0,1,0,1): 3451000/4556409; (0,1,1,0): -21475/72501; (1,0,0,1): 2293900/1518803; (1,0,1,0): -38570/72501; (1,1,0,0): 1. R'(0) = {21475/72501, -1725500/4556409, -38570/72501, 1146950/1518803} (x13,x14,x23,x24 order) equals the predicted limit R'_ij(0) = x0_ij (o+phi_i'+phi_j)/(o+phi1+phi2) (phi_i-phi_j)/(k-phi_i) exactly (also at all 12 grid points x {phi, phi'}); phi_u(0)/delta^Lambda = e1 exp(R'(0)) exactly; det R'(0) = 564593750/25411092993 != 0 (nonzero on all 12 grid points). Per-weight ranks of the eps^0 coefficient over the grid: all six weights rank 1 of dim 1, nonzero at 12/12 points, identical for phi' and for a second omega-direction (-2/3,5/7,-3/8,4/5). CONTROL om = 3/7 fixed: leading orders (0,0,1,1): O(eps^2); all other five weights O(1) - the gl(4) analogue of the gl(3) problem. (Our delta_j/delta_i - 1 -> (phi_i - phi_j)/(k - phi_i).)
```

## Deleted Cell 12 — 2026-10-10 11:06:30
**Kind:** markdown
```markdown
**Results, cell 8 (family and undo, exact, SS = 1, point zP, wP, om = 3/7, wk = 7/5), all True:** undo exact on all 3 HW kets (<v1(x)v2|h> = <phi_u(x)e1|h>); v2 rho(n)^-1 = e1 and (rho(n)(x)rho(n)).h = h; phi_u = delta^Lambda e1 exp(R'); v1 = delta^Lambda e1 exp(R_C) with R_C = R delta_j/delta_i (paper G14 shape); CR(R) = CR(z)(1-om z2z3)(1-om z1z4)/((1-om z2z4)(1-om z1z3)) at 3 fresh points; det R' != 0 at 3 points. **Blind-spot record:** for the old kappa-torus family (v2 = Omega_zw) det X' = 0 and the (0,0,SS,SS) slice of the undone covector is exactly 0 at 3 points x 3 values of wk.
```

## Deleted Cell 10 — 2026-10-10 11:06:30
**Kind:** markdown
```markdown
**Results, cell 7 (SS = 1):** HW sectors of V(x)V (weight, multiplicity, sector dim): (1,1,1,1): 1, 6; (2,1,1,0): 1, 2; (2,2,0,0): 1, 1 - the expected three irreps [1,1,1,1], [2,1,1,0], [2,2,0,0], each multiplicity 1 (dims 1 + 15 + 20 = 36). HW kets annihilated by all six total raising operators; G(x)G projector idempotent and fixes HW kets (exact).
```

## Deleted Cell 8 — 2026-10-10 11:06:30
**Kind:** markdown
```markdown
**Results, cell 6 (SS = 1, exact):** closed form x_ij = (-1)^i (z_i' - z_j)/(z1 - z2) prod_{k=1,2} (w_k - z_i)/(w_k - z_j), i' = other of {1,2} (found from a symbolic-w/numeric-z solve plus an integer-exponent fit over 16 random z; verified exactly at 3 fresh points). At P1: x = (-7/11, -665/999, 65/44, 5915/3996) = read-off. w3, w4 drop out; w1, w2 are pure torus directions (rescale x_ij by d_j/d_i); torus invariant x13 x24/(x14 x23) = (z1-z4)(z2-z3)/((z1-z3)(z2-z4)) (= 91/95 at P1; checked at 3 fresh points). Generic e1 exp(X) reaches all 6 weights. W_om = Omega (1 - om rho(z))/(1 - om z1 z2) = e1 exp(R) with R_ij = x_ij (1 - om z_i' z_j)/(1 - om z1 z2) exactly (3 points x om in {3/7, -5/3}), and its (0,0,1,1) component equals det R.
```

## Deleted Cell 6 — 2026-10-10 11:06:30
**Kind:** markdown
```markdown
**Results, cells 2-5 (SS = 1; executed in K1, VS Code does not execute cells of this file):** d1 = 6, d2 = 36, 6 distinct weights, all multiplicity 1; gl(4) relations of JJ True; e1 HW ket of weight (1,1,0,0); e1^T Ee[i,j] != 0 exactly for i <= 2 < j (i != j); nilradical E13,E14,E23,E24 commutes; G = diag(1,1,4,1,4,16) unique, positive, invariant. Cell 4 (exact rational, 0.57 s): quantum determinant t4 = (Prod z) nu1(u)nu2(u-h)nu3(u-2h)nu4(u-3h) deviation 0 (2 points); [t1,t2], [t1,t3], [t2,t3] deviations 0. Cell 5: Paux Gc(z) Paux^-1 = diag(z) and G'(z,z) = diag(z) symbolic True; Omega_zw: null space of t1^(1)(u0; G') - ev dim 1, u0-independent (u0 = 3/7, -5/11), eigen-bra of t1..t4^(1) with deviation 0 exactly, Omega_zw = e1^T exp(sum x_ij E_ij) exactly with the four exponentials commuting. Omega_zw(P1) = (1, 65/44, 5915/1998, 7/11, 1330/999, 1820/10989).
```

