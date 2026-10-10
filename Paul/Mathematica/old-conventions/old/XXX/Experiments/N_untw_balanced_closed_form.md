# Untwisted norm: the balanced even-degree family

Derived in `SU4_SoV_PerturbativeOverlap_L2.wb`, cells 44-66, 2026-09-10/11.
All computations use the live Wolfbook kernel and exact arithmetic.

Let

\[
L=S+P,\qquad \mu=(L-m,L-m,m,m),\qquad d=2m,\qquad n=L-2m,
\]

where (S,P) are positive integers and (0\le m\le\min(S,P)). The pole ladders
are assumed disjoint and simple, with nonzero (h). Put

\[
t=\frac{\theta_1-\theta_2}{h},\qquad a_1=\theta_1+Sh,\qquad a_2=\theta_2+Ph.
\]

These results concern this family. They do not cover, for example,
\(\mu=(L,L-2,2,0)\), which also has (d=2), or the odd-degree families.

## Results and their scope

The vacuum normalisation, proved previously, is

\[
N_{\rm vac}=\frac{(e_2-e_1)^4h^{4-4L}}
 {\left[S!P!\prod_{k=-P}^{S}(t+k)\right]^4},
\qquad e_\alpha=e^{2\pi i\theta_\alpha/h}.
\]

Let (U=P_1) and (V=P_2) be the first two monic **untwisted** direct Baxter
polynomials, of degrees (m,m+1). Choose the coefficient of (u^m) in (V)
to be zero. The monic Casoratian is

\[
W(u)=\bar Q_{12}(u)
 =-\frac{U(u)V(u-h)-U(u-h)V(u)}h.
\]

Write (A=U(a_1)U(a_2)) and (B=W(a_1)W(a_2)). Exact symbolic determinant
factorisation gives the following **general-spin identities**:

\[
\boxed{\frac{N_{d=2}}{N_{\rm vac}}
 =8(L+1)^3(L+2)\frac{B^3}{A^2}}
\]

and

\[
\boxed{\frac{N_{d=4}}{N_{\rm vac}}
 =432(L-1)^3L^4(L+1)\frac{B^3}{A^2}}.
\]

Here the Q-functions are those of the corresponding state; (A,B) differ
between (m=1) and (m=2). Quotients at exceptional zeros are interpreted by
the simplified rational expression, where the limit exists.

For arbitrary (m), the untwisted (8\times8) determinant below is derived.
Its simpler product factorisation is currently the conjecture

\[
\boxed{\frac{N_{2m}}{N_{\rm vac}}
 \stackrel{?}=m!(m+1)!^3(n+3)_m^3(n+4)_m\frac{B^3}{A^2}},
\]

where ((x)_m) is the rising factorial. It is proved for (m=0,1,2), and
independently validated against the derived determinant at (m=3,4) in the
cases listed below. An all-(m) product-factorisation proof remains open.

## Explicit inhomogeneity formulas

Set (A_k=(S+k+t)(P+k-t)). The two proven ratios are

\[
\frac{N_2}{N_{\rm vac}}
 =\frac{8(SP)^3(S+1)(P+1)}{(L+1)^3(L+2)}h^8 A_0^3A_1,
\]

\[
\frac{N_4}{N_{\rm vac}}
 =\frac{432(S-1)^3S^4(S+1)(P-1)^3P^4(P+1)}
 {(L-1)^3L^4(L+1)}h^{16} A_{-1}^3A_0^4A_1.
\]

Thus the old extrapolated (A_0^3A_1A_2A_3\) pattern does **not** extend to
this (d=4) family.

For (m=1), the endpoint polynomials are

\[
U(a_1)=\frac{(S+1)h(S+1+t)}{L+2},\qquad
W(a_1)=\frac{S(S+1)h^2(S+t)(S+1+t)}{(L+1)(L+2)},
\]

with the second endpoint obtained by (S\leftrightarrow P,t\mapsto-t).
For (m=2), they are

\[
U(a_1)=\frac{S(S+1)h^2(S+t)(S+1+t)}{L(L+1)},
\]

\[
W(a_1)=\frac{(S-1)S^2(S+1)h^4(S-1+t)(S+t)^2(S+1+t)}
 {(L-1)L^2(L+1)}.
\]

These follow from the symbolic untwisted Baxter solve, not from fitted values.

## Closed determinant for every balanced even degree

Let (widetilde U=P^4) and (widetilde V=P^3) denote the first two monic
untwisted dual polynomials in the notebook's ordering, with degrees (m,m+1)
and the same coefficient gauge. Define

\[
a=[u^{m-1}]U(u),\quad \tilde a=[u^{m-1}]\widetilde U(u),
\]

\[
R(u)=V(u)-(u-a)U(u),\qquad
\widetilde R(u)=-\widetilde V(u)+(u-\tilde a)\widetilde U(u).
\]

For (m=0), set (a=\tilde a=0); both (R,\widetilde R) vanish.
Otherwise their degrees are at most (m-1).

Let (\mathcal F) be the residue sum over both ladders and (\mathcal H)
the residue sum over the second ladder alone, with the original weights

\[
w(p)=\prod_{q\ne p}(p-q)^{-1}.
\]

For each of the eight columns ((\sigma,\tau,\beta)\in\{0,1\}^3), define
the polynomial in two **auxiliary coefficient-extraction variables** (v,w):

\[
G(u;v,w)=u^\beta
 [U(u-\tau h)+vR(u-\tau h)]
 [\widetilde U(u+\sigma h)+w\widetilde R(u+\sigma h)],
\]

and

\[
E(v,w)=1-h(\tau v+\sigma w)+\tfrac12h^2(\tau v+\sigma w)^2.
\]

The two row-generating polynomials, needed only through total degree two, are

\[
F(v,w)=E(v,w)\sum_{\ell=0}^{2}
 \frac{(v-w)^\ell}{(n+\ell)!}
 \mathcal F[u^{n+\ell}G(u;v,w)],
\]

\[
H(v,w)=E(v,w)\sum_{\ell=0}^{2}
 \frac{(v-w)^\ell}{\ell!}
 \mathcal H[u^\ell G(u;v,w)].
\]

For each column, take the four coefficients ([1],[w],[v],[vw]) of (F),
then the same four coefficients of (H). This defines an explicit matrix

\[
D_{(\alpha,r),(\sigma,\tau,\beta)},\qquad
\boxed{N_{\rm untw}=(e_2-e_1)^4\det D}.
\]

Every entry uses only the untwisted Q-polynomials, finite residue sums, and
factorials. There are no fitted constants or twisted-family derivatives to
solve for in this formula. The variables (v,w) only implement the four
specified coefficient extractions. The implementation is `balancedQMat` /
`balancedQNorm` in cell 60.

### Why this determinant is the leading coefficient

The total residue functional annihilates polynomials of degree at most (L).
Every original integrand has polynomial degree at most (2m+1). Consequently
each of the four total-functional rows is divisible by

\[
(\phi_i-\phi_j)^{L-2m}=(\phi_i-\phi_j)^n.
\]

Divide by these four factors first. The remaining matrix entries are analytic
at confluence. Within each functional block, replace the four branch rows by
their value, dual divided difference, direct divided difference, and mixed
divided difference. Their determinants extract precisely the remaining

\[
(\phi_1-\phi_2)^4(\phi_3-\phi_4)^4.
\]

The first branch differences are already encoded in the untwisted basis:

\[
\frac{q_1-q_2}{\phi_1-\phi_2}\longrightarrow R,
\qquad
\frac{q^4-q^3}{\phi_4-\phi_3}\longrightarrow\widetilde R.
\]

This follows by taking the corresponding differences of the full exponential
Q-functions, which solve the same Baxter equation, and fixing their monic
degree-(m+1) limits and coefficient gauges. Higher branch corrections
disappear from these divided-difference limits. Taylor expansion of the
analytic residue quotient now gives (F,H) above. The change from the periodic
functional to the second ladder alone contributes ((e_2-e_1)^4).

This is not a first-order reduced minor of the singular untwisted matrix.
The higher cross-pair zeros are removed analytically before taking the
confluent row differences.

## Symbolic factorisation and checks

At (h=1,\theta_2=0,\theta_1=t), normalised second-ladder factorial moments are

\[
\frac{\mathcal H[u^{\underline k}]}{\mathcal H[1]}
 =\frac{P^{\underline k}(S+t)^{\underline k}}{L^{\underline k}}.
\]

Here the underline denotes falling factorials. Ordinary moments follow by
Stirling numbers. For the total functional, only the moments of degrees

\[
L+1,L+2,L+3
\]

are needed after the row division: they equal (1,s_1,s_2), where (s_1,s_2)
are the first two complete homogeneous symmetric polynomials in all pole
locations. These identities produce the rational symbolic matrices in cells
57 and 64. Their determinants factor to the formulas above.

The generic moment calculation avoids small-(L) zero denominators by taking
(L\ge2m+3). For each fixed integer (P\ge m), both the finite-ladder moments
and the resulting formula are rational functions of (S,t); the identity on
all sufficiently large integer (S) extends to all admissible (S) where
the original sums are regular. Terms (k>P) must be set to zero before
specialising small (L). Direct finite-sum checks additionally cover small
spins.

| Claim/check | Evidence | Exact residual |
|---|---|---|
| Balanced (d=2), arbitrary spins | Symbolic determinant factorisation, cell 58 | 0 |
| Balanced (d=4), arbitrary spins | Symbolic determinant factorisation, cells 64-65 | 0 |
| (d=2), (2,2), truncations 1,2,3,4,6,8 | Independent reduced Baxter jets; a second ray at order 6 | 0 |
| (d=4), (3,2) | Independent order-4 Baxter jets, determinant order 12, cell 61 | 0 |
| (d=4), (2,2) | Agreement with the equal-weight `detU`, cell 60 | 0 |
| Product conjecture (d=6), (3,3),(4,3) | Derived untwisted determinant, cell 65 | 0 |
| Product conjecture (d=6), (4,3), theta2=1/5 | Independent inhomogeneity with isolated caches, cell 66 | 0 |
| Product conjecture (d=8), (4,4) | Derived untwisted determinant, cell 66 | 0 |

All residuals are exact rational/symbolic zeros; no working-precision or
rational-recognition tolerance is involved.

## Remaining work

1. Prove the proposed product factorisation for arbitrary (m). The finite
   untwisted determinant is already derived; its product simplification is
   the remaining question within this family.
2. Treat unbalanced irreps, starting with (\mu=(L,L-2,2,0)) at (d=2) and
   (\mu=(L-1,L-2,2,1)) at (d=3). Their direct/dual degree gaps introduce
   Laurent branches and unequal cross-pair powers, so the balanced proof
   must not be applied unchanged.
