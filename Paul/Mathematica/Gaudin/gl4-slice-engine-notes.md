# The slice engine for $\mathrm{gl}(4)$: highest-weight Clebsch–Gordan vectors and the two-site Gaudin missing labels

*Paul Ryan (project notes; draft compiled from the work packages G1–G5 of 22 October 2026) — 22 October 2026*

*(Markdown twin of `gl4-slice-engine-notes.tex`, generated from it by a small script: same sections, same numbering of equations, tables and footnotes. Math is written in `$…$` / `$$…$$`, tables are pipe tables, footnotes are collected at the end.)*

## Abstract

For two irreducible representations $V_\lambda$, $V_{\lambda'}$ of $\mathrm{gl}(4)$ the highest-weight (HW) Clebsch–Gordan vectors of $V_\lambda\otimes V_{\lambda'}$ are normally found as a null space on a weight space of dimension $n_W$ (up to 1749 in the SU(4) benchmark $[8,8,0,0]^{\otimes2}$ of the Munich code).
The criterion of Kostant and of Parthasarathy, Ranga Rao and Varadarajan says that a HW vector is determined by a *slice* of dimension $w\le5$ there, subject to three $\mathrm{sl}(2)$ string conditions.
Used as an algorithm it gives all 45 HW vectors of $[8,8,0,0]^{\otimes2}$ in 0.09 s on our machine, against 80 s for the Munich code on a different machine, and it reaches sizes the Munich code does not (dimension 19 481 in 4.6 s).
For multiplicities $m\ge2$ gl(4) has three missing labels.
We show that they are the $t$-independent Laurent coefficients $V_3,V_4^{a},V_4^{b}$ of the two-site Gaudin characteristic determinant: they commute and have a simple joint spectrum on every multiplicity space tested (180 exactly, 1864 at 50 digits, $m\le51$), whereas the cubic operator alone, which suffices for gl(3), does not.
They reduce to universal polynomials on the slice space, which gives canonical copies without any object of the size of $n_W$.
The oper side is forced: the Laurent data of the fourth-order Fuchsian oper are explicit in the triple and the central characters, and an exact Q-system solve reproduces the joint spectrum.
On the GT lattice the physical part of the spectrum on a weight space is cut out by an exact Weyl-type alternating sum over translated weights (verified on 832 sectors).
We say which statements are theorems, which are exact verifications on tested cases, and which are measurements.

## 1. Punchline

Take two copies of the SU(4) irrep $[8,8,0,0]$ (dimension 825, the largest case of the SU(4) sweep of the Munich code of Alex, Kalus, Huckleberry and von Delft [1]) and ask for the highest-weight (HW) vector of each of the 45 irreps $\mu$ in $V\otimes V$, written in the product of Gelfand–Tsetlin (GT) bases.
All 45 are multiplicity-free.
The largest, the singlet, has $n_W=1749$ components, and the standard algorithm finds it as the null space of the raising operators on that weight space.
The whole problem costs 80.2 s in the Munich code on an i5 laptop (18.7 s for the singlet), growing by a factor 5–7 per unit of $S$ in $[S,S,0,0]$.[^1]

*The problem is much smaller than that.*
A HW vector $\Psi$ of $V_\lambda\otimes V_{\lambda'}$ is determined by its *slice* $\psi_{\rm top}$: the components in which the second factor sits on its highest-weight state.
The slice lies in a single weight space of the first factor, of dimension $w\le5$ here (for the singlet, $w=1$).
The HW conditions fix every other component from the slice by a triangular recursion, and the conditions on the slice itself are $N-1=3$ $\mathrm{sl}(2)$ string conditions; they cut the slice space down to an $m$-dimensional subspace, $m$ being the multiplicity.
This is the theorem of Kostant and of Parthasarathy, Ranga Rao and Varadarajan [2, 3], and what we add is its use as an algorithm for gl(4).
One detail differs from gl(3): besides the consistency conditions of the recursion, the exact verification of a reconstructed vector needs the “dangling” conditions at the ends of the site-2 strings (Section 3).

*Missing labels.*
For $m\ge2$ the multiplicity space needs a distinguished basis, and gl(4) has $(N-1)(N-2)/2=3$ missing labels where gl(3) has one.
The two-site Gaudin model supplies them: the Laurent coefficients of the column determinant $\mathrm{cdet}(\partial_u-L(u))$ contain, besides central scalars, exactly one cubic and two quartic operators $V_3$, $V_4^a$, $V_4^b$, exactly independent of the pole separation $t=\theta_{12}$.
They commute and have a simple joint spectrum on the 180 multiplicity spaces with $m\ge2$ of 14 products tested exactly ($m\le12$) and on 1864 sectors at 50 digits ($m\le51$) (Section 4).
The cubic operator alone, which resolves the multiplicities for gl(3), is degenerate on 43 of the 180 sectors, in part for a reason of site-exchange symmetry; the pair $\{V_3,V_4^b\}$ separates every sector tested.

*Reduction, copies, opers.*
Each operator reduces to a universal polynomial in the site-1 generators (39, 78 and 85 words) that acts on the slice space and preserves the multiplicity space; joint diagonalisation gives the canonical copies (Section 5).
On the oper side the Laurent data are forced, with no fit, the polynomial solutions have degrees $\mu_j+4-j$, and an exact Q-system solve returns the joint spectrum on 180 sectors (Section 6).
On the GT lattice the three operators share one 21-point stencil, and the physical part of the spectrum on a weight space satisfies an exact Weyl-type alternating-sum identity over the translates $w\cdot\lambda'$ (832 sectors, 0 failures; verified, not proved) (Section 7).

*Numbers.*
All 45 HW vectors of $[8,8,0,0]^{\otimes2}$ take 0.026 s without writing the vectors out and 0.093 s with all of them written in double precision, on an i7 desktop, against 80.2 s for the Munich code on an i5 laptop; $S=10$ takes 0.23 s and $S=20$ (dimension 19 481, $3.1\cdot10^6$ coefficients) 4.6 s (Section 8).
With multiplicities, the slices of all 1592 canonical copies of $[9,6,3,0]^{\otimes2}$ ($m\le27$) take 6.4 s (23 s with 50-digit eigenvalues); the Munich SU(4) data are multiplicity-free, so there is no comparison there.

*Three points up front.*
First, the Munich code is the only general-$N$ code we know; the comparison is across two machines and against the code as distributed, and the SU(3)-specific codes that match or beat our engine for $N=3$ have no counterpart for $N=4$.
Second, the slice theorem and the reduction rule are theorems; commutation and simple spectrum on tested sectors, the intertwining relations, the Weyl identity and the closed forms are exact verifications, not proofs for all $\lambda,\lambda'$ (Section 9 tabulates the status of each statement).
Third, the Gaudin algebra as a missing-label resolver is classical [7, 8, 9, 10, 11] and for $N=3$ the operator is known [16, 17]; what is new is in Section 9, with the caveat that our literature check concerned gl(3).

## 2. The problem and the standard algorithm

**Punchline.** The standard algorithm solves an $n_W$-sized null space with dense linear algebra; on $[S,S,0,0]^{\otimes2}$ its cost grows by a factor 5–7 per unit of $S$ and is dominated by the singlet.

### 2.1. The problem

Let $V_\lambda$ be the irrep of $\mathrm{gl}(4)$ with highest weight $\lambda=[\lambda_1,\dots,\lambda_4]$, in the GT basis with the plain rational generators $E_{ij}$, $[E_{ij},E_{kl}]=\delta_{jk}E_{il}-\delta_{li}E_{kj}$ (Appendix A), and write $V_\lambda\otimes V_{\lambda'}=\bigoplus_\mu V_\mu\otimes\mathbb{C}^{m_\mu}$.
For every $\mu$ and each of its $m_\mu$ copies we want the HW vector: the vector $\Psi$ of weight $\mu$ annihilated by $E^{\rm tot}_{i,i+1}=E_{i,i+1}\otimes1+1\otimes E_{i,i+1}$, $i=1,2,3$ (the other raising operators follow by commutators).
The rest of $V_\mu$ follows by lowering.
The first factor is site 1 and the second site 2, as in the Gaudin chain; for $m_\mu\ge2$ one also needs a distinguished basis of the multiplicity space (Section 4).

The Munich benchmark family is $\lambda=\lambda'=[S,S,0,0]$, with $\dim V_\lambda=(S+1)(S+2)^2(S+3)/12$: 6, 20, 50, 105, 196, 336, 540, 825 for $S=1,\dots,8$.
The product has $(S+1)(S+2)/2$ irreps, all of multiplicity one.
The cost driver is $n_W$, the dimension of the weight-$\mu$ subspace of $V_\lambda\otimes V_{\lambda'}$ in which $\Psi$ lives; it is largest for the singlet $\mu=(S,S,S,S)$: 6, 22, 62, 147, 308, 588, 1044, 1749 for $S=1,\dots,8$, and 4290 at $S=10$, 86 779 at $S=20$.[^2]

### 2.2. What the standard code costs

The Munich code works with GT patterns and finds the HW coefficients as the null space of the raising operators on that weight space; everything else follows by lowering.[^3]
Its SU(4) sweep gives, as the sum of the highest-weight solves over all irreps, $2.8\times10^{-5}$ s, 0.86 ms, 10.7 ms, 57 ms, 0.36 s, 2.5 s, 15.3 s and 80.2 s for $S=1,\dots,8$ (13th-gen i5-1345U, g++ 11.4, LAPACK 3.10, one run per irrep).[^4]
The singlet alone takes 18.7 s at $S=8$, 23% of the total, because $n_W=1749$ is largest there.
The cost per irrep tracks the size of the linear system, of order $n_W^3$; the profile of the same routine for gl(3) traced it to a dense singular value decomposition (gl(3) notes, Section 2.2), and we did not repeat the profile for SU(4).
For gl(3) the routine became 19 times faster, 206 s to 10.9 s at $[14,7,0]^{\otimes2}$, after the SVD was replaced and an optimised BLAS linked; no such repair was made for SU(4), so every Munich number below is for the code as distributed.
The code was not run for $S\ge9$.

## 3. The slice theorem for gl(N) and the level recursion

**Punchline.** A HW vector is determined by its components in one weight space of the first factor, the slice; the slice must satisfy $N-1$ string conditions, and those are sufficient.
The multiplicity space is the joint null space of $N-1$ small matrices.
For $N=4$ the exact verification of a reconstructed vector needs conditions that are absent for $N=3$.

### 3.1. The slice and the level recursion

Write a HW vector as

$$
\Psi=\sum_q\psi_q\otimes e_q,\qquad\psi_q\in V_\lambda, \tag{1}
$$

with $q$ running over the GT basis of site 2, and call $\psi_{\rm top}$ the *slice*: the coefficient of the highest-weight state $e_{\rm top}$ of site 2.
If $\Psi$ has weight $\mu$ then $\psi_q$ has weight $\mu-\mathrm{wt}(q)$, so $\psi_{\rm top}$ lies in the weight space $V_\lambda[\mu-\lambda']$ of site 1, whose dimension $w$ is a weight multiplicity of $V_\lambda$.
(For $\lambda\ne\lambda'$ the choice of the factor placed at site 1 matters for $w$: for $[3,2,1,0]\otimes[2,1,0,0]$ the largest $w$ is 4, with the factors exchanged 2.)

Label the weights of site 2 by $\lambda'-\sum_ia_i\alpha_i$, $a\in\mathbb{Z}_{\ge0}^{N-1}$, let $Q(a)$ be the set of site-2 states of that weight and collect the $\psi_q$, $q\in Q(a)$, in a block $X_a$.
The HW condition for $E_i=E_{i,i+1}$ at the site-2 component $r\in Q(a)$ reads

$$
E_i\,\psi_r+\sum_{q\in Q(a+e_i)}(E_i)_{rq}\,\psi_q=0,\qquad i=1,\dots,N-1, \tag{2}
$$

where $E_i$ acts on site 1 and $(E_i)_{rq}$ is the site-2 matrix element.
Equation (2) determines the block one level down from those above: $A_bX_b=-\big(E_iX_{b-e_i}\big)_i$, with $A_b$ the stacked site-2 raising blocks from the weight $b$ to the weights $b-e_i$.
The matrix $A_b$ is injective for $b\ne0$, because a vector of $V_{\lambda'}$ killed by all simple raising operators below the top weight would be a second highest-weight vector; it depends neither on $\mu$ nor on the slice, so its left inverse is computed once per $(\lambda',b)$.
Processing the levels in lexicographic order determines $\Psi$ from the slice, and writing $\Psi$ out is one pass over the weight blocks of site 2, at a cost linear in $n_W$.[^5]

### 3.2. The theorem

Which slices are allowed?
Fix $i$ and follow the $\alpha_i$-string of site 2 through the top state, $e_{\rm top},E_{i+1,i}e_{\rm top},\dots,E_{i+1,i}^{k_i-1}e_{\rm top}$ with $k_i=\lambda'_i-\lambda'_{i+1}+1$.
Along it (2) makes the $k$-th component proportional to $E_i^k\psi_{\rm top}$, and at the end of the string, where there is no state further down, the sum in (2) is empty and the equation reads $E_i\psi_{\rm end}=0$, that is $E_i^{k_i}\psi_{\rm top}=0$.
These $N-1$ necessary conditions are also sufficient.

**Theorem** (Kostant; Parthasarathy–Ranga Rao–Varadarajan).
The map $\Psi\mapsto\psi_{\rm top}$ is an isomorphism from the HW vectors of weight $\mu$ in $V_\lambda\otimes V_{\lambda'}$ onto

$$
M_\mu=\{\,v\in V_\lambda[\mu-\lambda']\;:\;E_{i,i+1}^{\lambda'_i-\lambda'_{i+1}+1}v=0,\ \ i=1,\dots,N-1\,\}. \tag{3}
$$

In particular $\dim M_\mu=m_\mu$.

The statement is classical [2, 3, 4] (valid for any two dominant weights of any semisimple Lie algebra; Zhelobenko [5] in the physics literature).
We use it for $N=4$ in the form of three string conditions and a level recursion, and we check the implementation rather than the theorem (Section 3.3).
The conditions (2) at the ends of strings that do not pass through the top state, $E_i\psi_r=0$ for $r$ with $Q(a+e_i)=\emptyset$, are further conditions on components already fixed by the recursion; by the theorem they hold automatically when $\psi_{\rm top}\in M_\mu$.
We call them the *dangling* conditions (D), and the rows of the overdetermined stacked systems that the left inverse cannot absorb the *consistency* conditions (T); together T and D are exactly the rows of $E_i^{\rm tot}\Psi=0$.

### 3.3. Checks of the implementation

*Counting.*
$\dim M_\mu$ equals the Racah–Speiser multiplicity for every dominant candidate $\mu$, including those of multiplicity $0$ where $M_\mu=\{0\}$, on 26 products ($\lambda\ne\lambda'$ in both orders, $[S,S,0,0]^{\otimes2}$ for $S\le10$, and the multiplicity families $[4,2,0,0]^{\otimes2}$, $[6,4,2,0]^{\otimes2}$, $[9,6,3,0]^{\otimes2}$ with maximal $m=3,12,27$): 0 mismatches.
For $N=3$ the PRV bases coincide exactly with those of the gl(3) engine on $[4,2,0]^{\otimes2}$, $[6,3,0]^{\otimes2}$, $[14,7,0]^{\otimes2}$ (151 candidates).[^6]

*Recursion and exact twin.*
In double precision the HW residual of all $E_{ij}$, $i<j$, is at most $9.1\times10^{-15}$ on every irrep of $[S,S,0,0]^{\otimes2}$ for $S\le10$, $9.8\times10^{-15}$ on $[6,4,2,0]^{\otimes2}$ (297 vectors) and $4.4\times10^{-14}$ on $[9,6,3,0]^{\otimes2}$ (1592 vectors); weights are exact.
The same recursion in rational arithmetic (the *exact twin*), with all T and D conditions verified, gives 0 exact HW violations on every irrep of every product (81 528 conditions for $[10,10,0,0]^{\otimes2}$, 176 256 for $[6,4,2,0]^{\otimes2}$).[^7]

*Negative tests, and why D is needed.*
Of 1655 random exact vectors of $V_\lambda[\mu-\lambda']$ outside $M_\mu$ all are rejected by the exact conditions, and all 708 random elements of $M_\mu$ are accepted.
But 253 of the 1655 (15%) satisfy all consistency conditions T and are rejected only by D, in every product with a sector of $w-m=1$ (e.g. 20 of 20 for $[2,1,1,0]^{\otimes2}$).
For gl(3) none of about 500 such slices passed T, so the exact twin of the gl(3) engine does not check D; for gl(4) it must.[^8]

*Brute force.*
A Munich-style null space (stacked total raising operators on the weight space, no slice, no recursion) has dimension $m$ on every irrep of all products with $\dim V_\lambda\dim V_{\lambda'}\le10^4$, of $[S,S,0,0]^{\otimes2}$ for $S\le10$ and of $[4,2,0,0]^{\otimes2}$, $[6,4,2,0]^{\otimes2}$ (71 irreps with $n_W\le2600$, $m\le12$): 0 mismatches; our vectors lie in it to $1.6\times10^{-13}$ and span it.
Its smallest non-null eigenvalue of $A^{T}A$ is not 1 as for gl(3) but about $1/S$ for $[S,S,0,0]$ (0.586, 0.298 and 0.156 at $S=1,4,10$): far above double precision, so the slowness of the standard code is dense linear algebra and not conditioning.[^9]

### 3.4. How big is the problem now?

Table 1 shows what the theorem buys.[^10]
The singlet of $[8,8,0,0]^{\otimes2}$ has $n_W=1749$ but $w=1$: its slice is the lowest-weight vector of $V_\lambda$, a single number, and the recursion produces the rest.
For the rectangular family $w\le\lfloor S/2\rfloor+1$ while $n_W$ grows by a factor about 2.5 from $S=8$ to $S=10$.
With multiplicities the string conditions bite: at $(9,7,5,3)$ of $[6,4,2,0]^{\otimes2}$ the space $M_\mu$ is 12-dimensional inside a 15-dimensional weight space, while $n_W=916$.

**Table 1.** Problem sizes for seven irreps of three products: multiplicity $m$, slice length $w=\dim V_\lambda[\mu-\lambda']$ and number $n_W$ of components of the full vector.

| product | irrep $\mu$ | $m$ | $w$ | $n_W$ | the slice space $M_\mu$ |
| :--- | :--- | ---: | ---: | ---: | :--- |
| $[8,8,0,0]^{\otimes2}$ | $(8,8,8,8)$ (singlet) | 1 | 1 | 1749 | the lowest-weight vector of $V_\lambda$ |
|  | $(12,12,4,4)$ | 1 | 5 | 259 | a line in a 5-dimensional weight space |
| $[6,4,2,0]^{\otimes2}$ | $(6,6,6,6)$ (singlet) | 1 | 1 | 4743 | the lowest-weight vector of $V_\lambda$ |
|  | $(9,7,5,3)$ | 12 | 15 | 916 | 12 of 15 dimensions |
| $[3,2,1,0]^{\otimes2}$ | $(5,4,2,1)$ | 4 | 4 | 18 | all of $V_\lambda[\mu-\lambda']$ |
|  | $(4,4,2,2)$ | 2 | 4 | 68 | 2 of 4 dimensions |
|  | $(3,3,3,3)$ (singlet) | 1 | 1 | 152 | the lowest-weight vector of $V_\lambda$ |

## 4. The missing-label operators: the two-site Gaudin triple

**Punchline.** The Laurent coefficients of the two-site Gaudin $\mathrm{cdet}$ contain exactly one cubic and two quartic non-scalar operators; they are $t$-independent, commute and have a simple joint spectrum on every multiplicity space tested.
The cubic operator alone, sufficient for gl(3), is degenerate on 43 of 180 sectors for gl(4); the triple is the canonical definition.

### 4.1. The quantum determinant and its Laurent data

Let $e^s_{ij}$ be $E_{ij}$ acting on site $s$ ($e^1$ and $e^2$ commute), $L_{ij}(u)=e^1_{ij}/(u-\theta_1)+e^2_{ij}/(u-\theta_2)$, $M=\partial_u-L$ ($\partial_u$ times the unit matrix), and

$$
\mathrm{cdet}\big(\partial_u-L(u)\big)=\sum_{\sigma\in S_4}\mathrm{sgn}\,\sigma\;M_{\sigma(1)1}M_{\sigma(2)2}M_{\sigma(3)3}M_{\sigma(4)4}=\partial_u^4-S_1\partial_u^3+S_2\partial_u^2-S_3\partial_u+S_4, \tag{4}
$$

column-ordered, with every derivative moved to the right, as in the gl(3) Gaudin notes.[^11]
The $S_k(u)$ are rational with poles at $\theta_1,\theta_2$ of order at most $k$ and $S_k=O(u^{-k})$ at infinity; write $S_k(u)=\sum_{j=1}^kS_{k,j}\,(u-\theta_1)^{-j}+\text{regular}$ at $\theta_1$, the data at $\theta_2$ being obtained by exchanging the sites and $t\to-t$, with $t=\theta_1-\theta_2$.
This gives the 18 operators $S_{k,j}$, $k=2,3,4$, $j\le k$, at the two poles.
The convention is binding: for $N=3$ the engine reproduces $S_1,S_2,S_3$ of the gl(3) notes word for word (3, 8, 13 words), and the row-ordered variant is not even gl(N)-invariant (for $N=4$, 136 and 137 of 153 pairs of coefficients fail to commute on the two products tested).[^12]

*Homogeneity.*
Each word of $S_k$ has total weight $k$ ($\sum(d+1)=k$ over its factors $l^{(d)}_{ij}$), so every term of $S_{k,j}$ carries exactly $t^{\,j-k}$: $S_{k,j}=t^{\,j-k}\hat S_{k,j}$ with $\hat S_{k,j}$ independent of $t$.
Hence

$$
V_3=t\,S_{3,2},\qquad V_4^{a}=t\,S_{4,3},\qquad V_4^{b}=t^{2}S_{4,2} \tag{5}
$$

are exactly $t$-independent (proved; also checked at two values of $t$), with 88, 152 and 324 merged terms.[^13]
The top poles are central scalars, $S_{k,k}=T_k(\lambda)$ (Section 6.1); for instance $S_{2,2}=\lambda_2(\lambda_1+1)+\lambda_3(\lambda_1+\lambda_2+2)+\lambda_4(\lambda_1+\lambda_2+\lambda_3+3)$, while for $\lambda_4=0$ one has $S_{3,3}=(\lambda_1+2)(\lambda_2+1)\lambda_3$ and $S_{4,4}=0$, so for two-row irreps only the double pole is non-zero.
The independent non-scalar content of the 18 operators is one cubic and two quartic invariants, in line with the count $(N-1)(N-2)/2=3$ of missing labels: each $S_{k,j}$ at either pole is an integer combination of $V_3,V_4^a,V_4^b$ plus scalars (Table 3).

### 4.2. Commutativity and structure

Commutativity of all the Laurent coefficients is the theorem of Talalaev [14]; what we check is that the explicit term lists of the column-ordered determinant commute.
In exact arithmetic on the whole space ($t=31/35$ and $-3/7$) all 18 operators commute pairwise and with the 16 total generators on $[1,0,0,0]\otimes[1,1,0,0]$, $[1,1,0,0]\otimes[2,1,0,0]$, $[2,1,0,0]^{\otimes2}$ and $[2,1,1,0]^{\otimes2}$ (0 violations); in float64 also on $[3,2,0,0]^{\otimes2}$ and $[3,2,1,0]^{\otimes2}$ (dimensions 3600 and 4096, residual about $10^{-15}$).[^14]

The structure is exact in the PBW normal form of $U(\mathrm{gl}_4)\otimes U(\mathrm{gl}_4)$.
With $n_s=\sum_ie^s_{ii}$, $\Omega=\sum e^1_{ij}e^2_{ji}$, $\eta=n_1n_2-\Omega=t\,S_{2,1}$, the antisymmetrised forms

$$
D(a,b)_I=\frac1{a!\,b!}\sum_{\sigma,\tau\in S(I)}\mathrm{sgn}\,\sigma\,\mathrm{sgn}\,\tau\;e^1_{\sigma_1\tau_1}\cdots e^1_{\sigma_a\tau_a}\,e^2_{\sigma_{a+1}\tau_{a+1}}\cdots e^2_{\sigma_{a+b}\tau_{a+b}},\qquad|I|=a+b,
$$

$D(a,b)=D(a,b)_{1234}$ and $\Sigma D(a,b)=\sum_{|I|=a+b}D(a,b)_I$, and the polarised traces $T(s_1\dots s_r)=\sum e^{s_1}_{i_2i_1}e^{s_2}_{i_3i_2}\cdots e^{s_r}_{i_1i_r}$, one finds

$$
V_3=\Sigma D(2,1)+\eta=T(1,1,2)+\tfrac12n_1^2n_2+\tfrac12n_1n_2-n_1\Omega-\tfrac12n_2C_2^{1}+\Omega, \tag{6}
$$

$$
V_4^{a}=D(3,1)+V_3-\tfrac13\eta, \tag{7}
$$

$$
V_4^{b}=D(2,2)-D(3,1)+\tfrac12\big(\Sigma D(1,2)-\Sigma D(2,1)\big)-\tfrac16\eta. \tag{8}
$$

Here $C_2^1=\sum e^1_{ij}e^1_{ji}$; the second form of $V_3$ differs from $T(1,1,2)$ by a scalar on every multiplicity space, so on a sector $V_3\equiv T(1,1,2)\equiv-T(1,2,2)$ modulo scalars.
For general $N$ the coefficient of $\eta$ in (6) is $(N-2)/2$ (solved exactly for $N=3,\dots,6$), so for $N=3$ the same engine gives $V_3=V+\tfrac12\eta$, with $V$ the 36-term vertex operator of the gl(3) notes (an identity in $U(\mathrm{gl}_3)\otimes U(\mathrm{gl}_3)$), and the two sets of notes agree.
Equation (6) says that the cdet produces the gl(3)-type antisymmetrised vertex operator summed over the four 3-subsets of $\{1,2,3,4\}$ plus $\eta$.
In the trace basis the quartic operators contain the two traces of bidegree $(2,2)$ only in the combination $T(1,1,2,2)+\tfrac12T(1,2,1,2)$ (explicit forms in G2 §6): taken separately they do *not* commute with $V_3$ on 20 of the 22 sectors with $m\ge2$ of $[3,2,1,0]^{\otimes2}$ and $[4,2,1,0]^{\otimes2}$, nor with $V_4^a,V_4^b$ on 17.
The naive trace invariants are therefore not the missing-label operators; the quantum correction is the relative weight $\tfrac12$.[^15]

### 4.3. Spectrum: the cubic operator alone is not enough

On every sector tested the HW space is built exactly (null space over $\mathbb{Q}$), and the triple restricted to it commutes exactly and is self-adjoint for the invariant form $G\otimes G$ (exactly), hence has a real spectrum.
The joint spectrum is simple: an exact certificate (a generic rational combination has non-zero discriminant and Sturm count $m$) on 180 sectors with $m\ge2$ of 14 products ($\sum m=678$ copies, $m\le12$, up to $[6,4,2,0]^{\otimes2}$), and at 50 digits on 1864 sectors of twelve family members ($m\le51$).
Simplicity of the spectrum of the *whole* Gaudin algebra on the multiplicity spaces of a tensor product of finite-dimensional modules is a known theorem [8, 9, 10, 11]; since the triple generates that algebra on a multiplicity space modulo scalars, it holds for the triple as well, and our certificates check the identification of the triple and the conventions, not the theorem.[^16]

The new fact for $N=4$ is that the cubic operator alone does not separate the copies.
The exact discriminant of $V_3$ vanishes on some sectors of every larger product with $\lambda=\lambda'$: 2, 2, 6, 6, 8 and 19 sectors of $[3,2,1,0]^{\otimes2}$, $[4,2,1,0]^{\otimes2}$, $[4,3,1,0]^{\otimes2}$, $[5,3,1,0]^{\otimes2}$, $[5,3,2,0]^{\otimes2}$ and $[6,4,2,0]^{\otimes2}$; it is simple on every sector of the products with $\lambda\ne\lambda'$ and of the small $\lambda=\lambda'$ ones.
Two examples are in Table 2: at $(4,4,2,2)$ of $[3,2,1,0]^{\otimes2}$, $m=2$ and $V_3=81$ is a scalar; at $(5,4,2,1)$, $m=4$ and $V_3=57$ occurs twice; in both the triple has distinct joint values.

**Table 2.** Two sectors of $[3,2,1,0]^{\otimes2}$ where $V_3$ is degenerate: the joint eigenvalue triples $(V_3,V_4^a,V_4^b)$ (rounded).

| $\mu$ | $m$ | triples | site exchange $P$ |
| :--- | ---: | :--- | :--- |
| $(4,4,2,2)$ | 2 | $(81,\,76.748,\,21.757)$, $(81,\,49.252,\,104.243)$ | both copies fixed, parities $++$ |
| $(5,4,2,1)$ | 4 | $(67.392,\,70.981,\,-16.981)$, $(46.608,\,19.019,\,34.981)$, | 2-cycles $(1,4)$, $(2,3)$; |
|  |  | $(57,\,50.697,\,21.303)$, $(57,\,21.303,\,50.697)$ | $V_3=57=c/2$ on the pair $(2,3)$ |

*Part of the degeneracy is forced.*
For $\lambda=\lambda'$ the exchange $P:(\nu,q)\mapsto(q,\nu)$ preserves every HW space and conjugates the family: $PV_3P=-V_3+c$, $PV_4^aP=-2V_4^a-V_4^b+c'$, $PV_4^bP=3V_4^a+2V_4^b-c'$ with scalars $c,c'$ (exact, 57 of 57 sectors; the relations follow from the dictionary of Section 6.1, where $c=\tau_3$, $c'=\tau_4$).
So $V_3-c/2$ is odd under $P$.
If $m_+$ and $m_-$ are the multiplicities of $V_\mu$ in the symmetric and antisymmetric parts of $V_\lambda\otimes V_\lambda$, an odd operator maps the $P$-even subspace of the multiplicity space into the $P$-odd one and back, so it has a symmetric spectrum and a kernel of dimension at least $|m_+-m_-|$.
Hence $V_3$ is degenerate whenever $|m_+-m_-|\ge2$: this is the case at $(4,4,2,2)$ ($m_+=2$, $m_-=0$) and, among the 19 degenerate sectors of $[6,4,2,0]^{\otimes2}$, for 13; the other 6 degeneracies, like $(5,4,2,1)$ above, are accidental (a 2-cycle whose two partners both have $V_3=c/2$).
For $\lambda\ne\lambda'$ there is no such symmetry.[^17]

*Which subsets separate.*
Counting the sectors with $m\ge2$ whose joint spectrum is degenerate (50-digit eigenvalues, two counted equal below $10^{-20}$): $V_3$ alone fails on 2 of 8, 19 of 56, 69 of 179 and 169 of 409 sectors at $[3,2,1,0]^{\otimes2}$, $[6,4,2,0]^{\otimes2}$, $[9,6,3,0]^{\otimes2}$, $[12,8,4,0]^{\otimes2}$; each single quartic and $\{V_4^a,V_4^b\}$ fail on 1 of 8, 6 to 8 of 56, 13 of 179 and 25 to 27 of 409.
The pair $\{V_3,V_4^b\}$ and the triple separate every sector of every product tested (1864 sectors at 50 digits; all 14 products exactly), whereas $\{V_3,V_4^a\}$ fails on the $m=12$ sector $(9,7,5,3)$ of $[6,4,2,0]^{\otimes2}$ and on one sector at $[12,8,4,0]^{\otimes2}$.
On the two-row families $[2a,a,0,0]^{\otimes2}$ the operator $V_4^a$ vanishes and $V_4^b$ is a function of $V_3$, so the problem is gl(3)-like and $V_3$ alone separates (1212 sectors, $a\le8$).
We therefore take the triple, not $V_3$, as the definition of the canonical multiplicity basis: the eigenbasis of the Gaudin algebra on the multiplicity space, which is also the basis separation of variables selects.

## 5. Reduction to the slice and the canonical copies

**Punchline.** Applied to a HW vector and read off at the slice, each missing-label operator becomes a universal polynomial in the site-1 generators acting on the slice space: after normal ordering of the site-2 factors, lowering factors kill, Cartan factors take their value at $\lambda'$ and raising factors are transferred to site 1 with a sign.
The reduced triple preserves $M_\mu$, commutes, has the spectrum of the full operators, and is compatible with the PRV conditions through an intertwining relation.

### 5.1. The rule

For a word $W^2$ of site-2 generators write its normal form in $U(\mathrm{gl}_4)$ in the class order lowering, Cartan, raising, $W^2=\sum_wc_w\,L_wH_wR_w$, with $L_w$ a product of lowering generators $E_{ab}$ ($a>b$), $H_w$ of Cartan generators and $R_w=E_{r_1}\cdots E_{r_k}$ of raising generators ($a<b$).

**Proposition** (reduction rule).
For every word $W^1$ of site-1 generators and every HW vector $\Psi$ with slice $\psi$,

$$
(\mathrm{id}\otimes\langle\mathrm{top}|)\,W^1W^2\,\Psi=\sum_{w:\;L_w=1}c_w\prod_{E_{aa}\in H_w}\lambda'_a\;(-1)^{k_w}\;W^1E_{r_k}\cdots E_{r_2}E_{r_1}\,\psi . \tag{9}
$$

*Proof.* (i) $\langle\mathrm{top}|E_{ab}=0$ for $a>b$, because a lowering operator never produces the top state; (ii) $\langle\mathrm{top}|E_{aa}=\lambda'_a\langle\mathrm{top}|$; (iii) $E^{\rm tot}_{ab}\Psi=0$ for $a<b$ means $e^2_{ab}\Psi=-e^1_{ab}\Psi$, and since $e^2$ commutes with $e^1$ this holds on any site-1 word applied to $\Psi$, so $R\Psi=(-1)^ke^1_{r_k}\cdots e^1_{r_1}\Psi$: the rightmost raising factor is transferred first and the transferred factors stand in reversed order. Site-1 operators commute with $\langle\mathrm{top}|$. $\square$

Hence the reduction $\mathrm{Red}_{\lambda'}:U(\mathrm{gl}_4)\otimes U(\mathrm{gl}_4)\to U(\mathrm{gl}_4)$ is a well-defined linear map, polynomial in $\lambda'$ and independent of $\lambda$ and $\mu$.
If $O$ commutes with the diagonal $\mathrm{gl}_4$, then $O\Psi$ is again a HW vector, so $\mathrm{Red}(O)$ maps $M_\mu$ into itself and $O\mapsto\mathrm{Red}(O)|_{M_\mu}$ is the matrix of $O$ on the multiplicity space in slice coordinates (an algebra homomorphism).
For one site-2 factor the rule is the gl(3) one-liner of the gl(3) notes (lowering $\to0$, Cartan $\to\lambda'_c$, raising $E_{cd}\to-E_{cd}$); the new ingredient is the normal ordering for words with two or three site-2 factors.[^18]
The ordering must not be skipped: $E_{12}E_{21}=E_{21}E_{12}+E_{11}-E_{22}$ contributes $\lambda'_1-\lambda'_2$ to the reduction of a site-2 word $e^2_{12}e^2_{21}$, while $e^2_{21}e^2_{12}$ is already in normal order and dies.
On the column-ordered cdet term lists it is a no-op, because the column order puts a raising factor only to the left of lowering factors it commutes with; but for the antisymmetrised representative of $V_4^b$ in (8) (1464 terms) the rule with normal ordering reproduces the cdet operator on 449 of 449 PRV vectors, and the naive rule differs on 441.[^19]

### 5.2. Universal polynomials and closed forms

The reduction gives universal polynomials in the site-1 generators with coefficients in $\mathbb{Q}[\lambda']$: $V_3$ (88 terms) reduces to 39 words, $V_4^a$ (152) to 78 and $V_4^b$ (324) to 85, of length at most 3, 4 and 4 (the 36-term vertex operator of gl(3) reduces to 24 words).
The reduction costs 0.3–1.2 ms for the universal polynomials and about 4 ms per $\lambda'$ for all three operators including the preparation for their application, and does not grow with the representation.[^20]
The 39 words of $V_3$ are in Appendix B.
For $V_3$ and $V_4^a$ every term has one site-2 factor, so (6)–(7) and the one-line rule give closed forms.
With $\Phi_{cc}=\lambda'_c$, $\Phi_{ab}=-E_{ab}$ for $a<b$, $\Phi_{ab}=0$ for $a>b$ and $|\lambda|=\sum_i\lambda_i$,

$$
V_{\mathrm{red},3}=\sum_{|I|=3}\tfrac12\sum_{\sigma,\tau\in S(I)}\mathrm{sgn}\,\sigma\,\mathrm{sgn}\,\tau\;E_{\sigma_1\tau_1}E_{\sigma_2\tau_2}\Phi_{\sigma_3\tau_3}+|\lambda||\lambda'|-\sum_i\lambda'_iE_{ii}+\sum_{i>j}E_{ij}E_{ji}, \tag{10}
$$

$$
V_{\mathrm{red},4a}=\tfrac16\sum_{\sigma,\tau\in S_4}\mathrm{sgn}\,\sigma\,\mathrm{sgn}\,\tau\;E_{\sigma_1\tau_1}E_{\sigma_2\tau_2}E_{\sigma_3\tau_3}\Phi_{\sigma_4\tau_4}+V_{\mathrm{red},3}-\tfrac13\Big(|\lambda||\lambda'|-\sum_i\lambda'_iE_{ii}+\sum_{i>j}E_{ij}E_{ji}\Big). \tag{11}
$$

The first term of (10) is the gl(3) reduced vertex operator, eq. (6) of the gl(3) notes, on each 3-subset of $\{1,2,3,4\}$.
Expanded, (10) has 71 words and differs as a word list from the 39 words of the engine, but the two agree on every PRV vector of nine products (898 operator–vector pairs) and, in a check made for these notes, on every basis vector of $V_\lambda$ for three $\lambda$ and five $\lambda'$ (990 comparisons): they are the same operator.[^21]
The scalars of the relations $t^2S_{3,1}=-2V_3+c$ and $t^3S_{4,1}=-3(V_4^a+V_4^b)+c'$ were first found as exact polynomial fits in the Casimirs (252 irreps); they are $c=\tau_3$ and $c'=\tau_4$ of Section 6.1 and the exchange constants of Section 4.3.[^22]

### 5.3. Verification against the full operators

For $N=3$ the reduced $t\,S_{3,2}$ equals the gl(3) engine's $V_{\rm red}$ plus $\tfrac12\eta$ ($\eta$ of eq. (5) of the gl(3) notes) exactly on all 35 PRV basis vectors of $[2,1,0]^{\otimes2}$, $[3,1,0]^{\otimes2}$, $[4,2,0]^{\otimes2}$.
For $N=4$, take the exact HW vectors $\Psi_j$ of the PRV basis, apply the full two-site operator (the term lists of the file) exactly and read off the slice: it equals the reduced operator applied to the slice, and the whole vector $O\Psi_j$ equals $\sum_iA_{ij}\Psi_i$.
This holds exactly for $V_3,V_4^a,V_4^b,S_{3,1},S_{4,1}$ on 46 irreps (64 PRV vectors) of six products (320 + 320 checks), on 598 further vectors of 11 products including $\lambda\ne\lambda'$ in both orders and $[6,4,2,0]^{\otimes2}$ (297 vectors), on 2275 slice checks of the $\theta_2$ Laurent coefficients (sites exchanged), and in float64 to $1.7\times10^{-14}$: 0 mismatches.
The tests are sensitive: omitting the sign $(-1)^k$ gives wrong slices on 29 of 30 vectors of $[3,2,1,0]^{\otimes2}$ for $V_3$, and transferring in the original instead of the reversed order gives wrong slices for $V_4^b$ on 22 of 30.[^23]

### 5.4. Three exact structural facts on the weight space

Write $W=V_\lambda[\mu-\lambda']\supset M_\mu$, $k_i=\lambda'_i-\lambda'_{i+1}+1$ and $s_i\cdot\lambda'=\lambda'-k_i\alpha_i$.
Because the reduction is polynomial in $\lambda'$, $V_{\mathrm{red}}(s_i\cdot\lambda')$ makes sense for the non-dominant weight $s_i\cdot\lambda'$.

1. *BGG intertwining.* For each $k\in\{3,4a,4b\}$, $i=1,2,3$ and every GT basis vector of every weight space,

   $$
   E_i^{k_i}\circ V_{\mathrm{red},k}(\lambda')=V_{\mathrm{red},k}(s_i\cdot\lambda')\circ E_i^{k_i}\qquad(W\to V_\lambda[\mu-s_i\cdot\lambda']). \tag{12}
   $$

   5814 tests on nine products, 1130 of them non-trivial (one side non-zero), 0 failures.
   Since $M_\mu=\bigcap_i\ker E_i^{k_i}$, this is a second, representation-theoretic reason why the reduced operators preserve $M_\mu$ (the first is Section 5.1); it involves the non-dominant weight $s_i\cdot\lambda'$, for which there is no finite-dimensional module.
2. *Commutation on the whole weight space, for every $\lambda'$.* The three reduced operators commute on all of $W$, not only on $M_\mu$ (nine products, 0 failures), and on every weight space of dimension $\ge2$ of $V_\lambda=[3,2,1,0],[4,2,0,0],[2,1,1,0]$ for generic non-dominant rational $\lambda'$ such as $(\tfrac13,\tfrac57,-2,\tfrac45)$ and for $\lambda'=0$ (168, 456 and 12 tests, 0 failures).
3. *$M_\mu$ is a spectral subspace.* The joint spectrum on each whole weight space $W$ is simple ($w$ distinct triples at 50 digits) and the $m$ triples of $M_\mu$ are among them; the other $w-m$ eigenvectors violate a string condition.
   For $[6,4,2,0]^{\otimes2}$ the 76 weight spaces have $\sum w=424$ and $\sum m=297$, so 127 of the 424 eigenvectors are spurious.

A natural interpretation, *not proved here*: for generic $\lambda'$ the operators $V_{\mathrm{red}}(\lambda')$ are the Gaudin missing-label operators on the singular vectors of $V_\lambda\otimes M(\lambda')$, $M(\lambda')$ a Verma module, the slice identifies the singular vectors of weight $\mu$ with $W$, and (12) says the operators are compatible with the embeddings $M(s_i\cdot\lambda')\hookrightarrow M(\lambda')$.[^24]

### 5.5. Canonical copies and the engine

The canonical copies are the joint eigenvectors of the reduced triple on $M_\mu$, from the 50-digit eigen-decomposition of the generic combination $V_3+\tfrac7{19}V_4^a+\tfrac{11}{13}V_4^b$, after which every reduced operator is conjugated to the eigenbasis and its largest off-diagonal entry measured: at most $8.5\times10^{-43}$ on all sectors of all products (worst case $[12,8,4,0]^{\otimes2}$, $m=51$).
Copies are ordered lexicographically by $(V_3,V_4^a,V_4^b)$ and normalised so that the first non-zero slice component is 1; the full vectors follow from the level recursion.
All eigenvalues are real, each reduced operator is self-adjoint for the invariant form (exact), the copies are orthogonal for $G\otimes G$ ($5.8\times10^{-47}$ with the exact Gram matrix, $10^{-15}$ with float vectors), HW residuals are at most $9.9\times10^{-15}$, and triples and slices agree with the exact HW eigenbasis of the full operators on six products (to $1.3\times10^{-39}$ and $6.7\times10^{-39}$).
For $\lambda=\lambda'$ the affine map $T(e_3,e_4^a,e_4^b)=(c-e_3,\,-2e_4^a-e_4^b+c',\,3e_4^a+2e_4^b-c')$ permutes the joint triples on all 156 sectors tested, and the overlap matrices of the full vectors are signed permutation matrices.[^25]

On the multiplicity-free family the triple is a set of scalars with closed forms: for $[S,S,0,0]^{\otimes2}$, $V_{\rm red,4a}=0$ and, with $C_2(\nu)=\sum_i\nu_i(\nu_i+5-2i)$,

$$
V_{3}(\mu)=\tfrac{S+1}2\big(C_2(2S,2S,0,0)-C_2(\mu)\big),\qquad V_4^{b}(\mu)=T_4(\mu)=\prod_{c=1}^4(\mu_c+4-c), \tag{13}
$$

checked for $S\le10$ (285 irreps), not proved.
G3 states $V_4^b$ through Perelomov–Popov Casimirs and G5 as the central character of Section 6.1; their equality was checked for these notes on the same 285 irreps.[^26]

*The engine in five steps.*
(1) Build exact sparse GT matrices for $V_\lambda$, $V_{\lambda'}$ (simple generators and Cartan matrices; the six non-adjacent ones by commutators) and list irreps and multiplicities by Racah–Speiser.
(2) For each $\mu$ compute $M_\mu\subset W$, the exact joint null space of the three string operators.
(3) Reduce $V_3,V_4^a,V_4^b$ once per $\lambda'$ and apply them to the basis of $M_\mu$: exact $m\times m$ matrices.
(4) Diagonalise jointly (float64 or 50 digits): the eigenvectors are the slices of the canonical copies.
(5) Optionally write out the full vectors by the level recursion.
Steps 1–4 never see $n_W$; step 5 is bound by the size of the output.

## 6. The oper side

**Punchline.** On a joint eigenvector the cdet becomes a fourth-order Fuchsian operator with three singular points.
Its Laurent data are forced: every coefficient is an integer combination of the triple plus central characters, with no fit.
Its polynomial solutions have the degrees $\mu_j+4-j$, and an exact elimination in the three unknowns $(V_3,V_4^a,V_4^b)$ returns exactly the joint spectrum of the triple; all four degree conditions are needed.

### 6.1. The dictionary

On a joint eigenvector the operators $S_k(u)$ of (4) become scalar functions $s_k(u)=\sum_{s=1,2}\sum_{j=1}^kS^{(s)}_{k,j}(u-\theta_s)^{-j}$ (no regular part, since $S_k\to0$ at infinity), and the oper is $D=\partial^4-s_1\partial^3+s_2\partial^2-s_3\partial+s_4$.
Put $y_c(\nu)=\nu_c+4-c$, expand $\prod_c(\theta-y_c)=\sum_ka_k(\nu)\,\theta(\theta-1)\cdots(\theta-k+1)$ in falling factorials and set $T_k(\nu)=(-1)^ka_{4-k}(\nu)$:

$$
\begin{aligned}
T_1&=|\nu|,\qquad T_2=e_2(y)-3e_1(y)+7,\\
T_3&=e_3(y)-e_2(y)+e_1(y)-1,\qquad T_4=e_4(y)=\textstyle\prod_cy_c,
\end{aligned} \tag{14}
$$

$e_k$ the elementary symmetric functions of $y$; $T_k(\nu)$ is the central character of the single-site top pole, $S_{k,k}=T_k(\nu)$ on $V_\nu$.
Let $\tau_k=T_k(\mu)-T_k(\lambda)-T_k(\lambda')$; then $\tau_2=n_1n_2-\Omega$ with $\Omega=\tfrac12\big(C_2(\mu)-C_2(\lambda)-C_2(\lambda')\big)$.

**Proposition** (dictionary).
On $M_\mu$ the Laurent data are those of Table 3.
*Derivation.* At fixed $k$ the unknowns are the $2(k-1)$ non-top Laurent coefficients at the two poles; the top poles are central.
Since $s_k=O(u^{-k})$ and, because $L\to E^{\rm tot}/u$, the coefficient of $u^{-k}$ is the single-site word polynomial with $E\to E^{\rm tot}$, i.e. $T_k(\mu)$ on $V_\mu$, the expansion at infinity gives $k$ linear relations (the coefficient of $u^{-r}$ vanishes for $r<k$ and equals $T_k(\mu)$ for $r=k$).
That leaves $2(k-1)-k=k-2$ free parameters, namely $0,1,2$ for $k=2,3,4$: one cubic and two quartic labels, which we call $V_3,V_4^a,V_4^b$ by (5).
Solving gives Table 3; its integer part is the table found independently in G2 §6 by exact linear algebra in the PBW algebra.

**Table 3.** The dictionary. Left: $t^{k-j}S^{(1)}_{k,j}$ at $\theta_1$; right: $(-t)^{k-j}S^{(2)}_{k,j}$ at $\theta_2$; $\tau_k=T_k(\mu)-T_k(\lambda)-T_k(\lambda')$.

| $k$ | $j$ | at $\theta_1$ | at $\theta_2$ |
| :---: | :---: | :--- | :--- |
| 2 | 2 | $T_2(\lambda)$ | $T_2(\lambda')$ |
| 2 | 1 | $\tau_2$ | $\tau_2$ |
| 3 | 3 | $T_3(\lambda)$ | $T_3(\lambda')$ |
| 3 | 2 | $V_3$ | $-V_3+\tau_3$ |
| 3 | 1 | $-2V_3+\tau_3$ | $2V_3-\tau_3$ |
| 4 | 4 | $T_4(\lambda)$ | $T_4(\lambda')$ |
| 4 | 3 | $V_4^a$ | $-2V_4^a-V_4^b+\tau_4$ |
| 4 | 2 | $V_4^b$ | $3V_4^a+2V_4^b-\tau_4$ |
| 4 | 1 | $-3V_4^a-3V_4^b+\tau_4$ | $-3V_4^a-3V_4^b+\tau_4$ |

The proposition is a theorem given the central characters (14), which come from the Capelli-type evaluation $\mathrm{cdet}(\partial-E/x)=x^{-4}\prod_c(\theta-y_c)$ on the HW vector ($\theta=x\partial_x$) and were checked on all 70 partitions with $\nu_1\le4$.
It was also verified sector by sector: the 18 Laurent coefficients, restricted exactly to $M_\mu$ and corrected by the integer combination of the table, leave the scalar matrix of the closed form on all 95 sectors ($\sum m=146$ copies) of nine products (1710 sector–operator tests, 0 violations; also at $t=31/35$ and $-3/7$).
In trace-form Casimirs, $T_2=-\tfrac12C_2+\tfrac12n^2+\tfrac32n$ and $T_3=\tfrac13C_3-\tfrac12nC_2-\tfrac13C_2+\tfrac16n^3+\tfrac56n^2+2n$ (fitted on 63 partitions, verified on 20 held-out; $T_4$ in G5 §2.3).[^27]
In gl(3) the corresponding relation was $V=a_1-\tfrac12\eta$; here $V_3=tS_{3,2}$ with no shift, and $\tau_2,\tau_3,\tau_4$ play the role of $\eta$.

### 6.2. Exponents, polynomial solutions, Wronskian

The indicial polynomials of $D$ at $\theta_1$, $\theta_2$ and $\infty$ do not depend on the triple: they are $\prod_c(e-(\lambda_c+4-c))$, $\prod_c(e-(\lambda'_c+4-c))$ and $\prod_c(n-(\mu_c+4-c))$ (134 of 134 sectors, as identities of affine forms in the triple), the Fuchs relation (sum of exponents $=6$) holds, and polynomial solutions can only have the degrees $d_j=\mu_j+4-j$, $d_1>d_2>d_3>d_4$.
For all 210 copies of the 134 sectors of 11 products, the null space of the pencil of $D$ on polynomials of degree $\le d_j$ has dimension $5-j$ (exact rank on the 88 rational triples; a 60-digit singular value decomposition with a gap of at least $3\times10^{59}$ on the 122 irrational copies), the four polynomial solutions have Wronskian proportional to $(u-\theta_1)^{|\lambda|}(u-\theta_2)^{|\lambda'|}$ (deviation $\le2\times10^{-59}$), and their vanishing orders at $\theta_1,\theta_2$ are $\{\lambda_c+4-c\}$, $\{\lambda'_c+4-c\}$; all 630 perturbed and 210 random triples are rejected.[^28]
That the oper has these polynomial solutions is what the Gaudin/oper correspondence leads one to expect [12, 13, 8]; the statement above is its explicit verification, copy by copy.

### 6.3. The exact Q-system solve

Let $v=(v_3,v_4^a,v_4^b)$ be unknown.
That a polynomial solution of exact degree $d_j$ exists is a rank condition on a $(d_j+5)\times(d_j+1)$ pencil affine in $v$.
The Frobenius recursion at infinity, with the free parameters at resonances fixed to zero, turns it into 7, 6, 5, 4 polynomial conditions in $v$ for $j=1,\dots,4$ (22 in all, of degree $\le9$, and $\le10$ in the largest sectors, while $d_1\le13$); let $J$ be the ideal they generate.
Each generator vanishes exactly on the commuting restrictions of $(V_3,V_4^a,V_4^b)$ to $M_\mu$, so $J\subseteq K$, the ideal of the $m$ joint triples; a Gröbner basis gives $\dim_{\mathbb{Q}}\mathbb{Q}[v]/J=m$, the number of distinct points, hence $J=K$: the solution set $\Sigma$ is exactly the joint spectrum, each point simple (Jacobian of rank 3).
This holds on all 165 sectors of 12 products ($\sum m=279$, $m\le5$, 64 sectors with $m\ge2$) and on the 15 largest-multiplicity sectors of $[5,3,1,0]^{\otimes2}$, $[5,3,2,0]^{\otimes2}$, $[6,4,2,0]^{\otimes2}$ ($m\le12$): 180 sectors in all.
Points computed from $J$ alone agree with the triples of the full operators to $2.8\times10^{-58}$ at worst on the 165 sectors and to $10^{-55}$ on the $m=12$ sector $(9,7,5,3)$ of $[6,4,2,0]^{\otimes2}$ ($n_W=916$, $d=(12,9,6,3)$), the one where $\{V_3,V_4^a\}$ fails, which takes 17 s including the HW space; sectors with $m\le5$ take at most 1 s.[^29]
For $(5,4,2,1)$ of $[3,2,1,0]^{\otimes2}$, where $V_3=57$ occurs twice, $\Sigma$ consists of the four distinct triples of Table 2.

*All four degree conditions are needed.*
Take a subset $S$ of the four conditions and ask whether its solution set is exactly the joint spectrum, on the 68 sectors of eight small products (Gröbner basis, 5 s limit).
It is for $S=\{4\},\{3\},\{2\},\{1\}$ in 0, 0, 6 and 7 of the 68 sectors ($\{4\}$ alone is a curve in all 68), for the best pair $\{4,3\}$ in 32, for $\{4,3,2\}$ in 60, and for all four in all 68 (all 16 with $m\ge2$): no proper subset works in every sector.
For gl(3) the analogous lesson was that two of the three Baxter systems over-count roots in the $E_{12}$-cut sectors; here $\{4,3,2\}$ fails in 8 sectors, all with $\mu_1=\mu_2$ (adjacent top exponents), where the condition for $j=1$ is indispensable ($j=1$ alone gives exactly the three triples at $(4,4,3,1)$ of $[3,2,1,0]^{\otimes2}$); only 8 of the 26 sectors with $\mu_1=\mu_2$ fail.[^30]

### 6.4. $\mu_4=0$ and the rectangular family

On every sector with $\mu_4=0$ (and $\lambda_4=\lambda'_4=0$) the oper forces $V_4^a=V_4^b=0$: the constant is the polynomial solution of degree $d_4=\mu_4=0$, so $s_4\equiv0$, and in Table 3 $\tau_4=0$ because every $y_4$ vanishes.
This was checked on all 47 such sectors of the 11 products of Section 6.2 (54 copies, 7 with $m\ge2$: the triples of the full operators are $(v_3,0,0)$ to $10^{-40}$), and it follows from the polynomial-solution statement above; the three-parameter solve is then the one-parameter, gl(3)-type solve in $V_3$.
In the rectangular family $[S,S,0,0]\otimes[P,P,0,0]$ one finds $V_4^a\equiv0$ and $V_4^b=T_4(\mu)$ (Section 6.5).[^31]

### 6.5. The July-2026 gl(4) Gaudin ODE

The notebook `GaudinSolver.wb` of July 2026 solves the two-site gl(4) Gaudin model for $[S,S,0,0]\otimes[P,P,0,0]$ by a fourth-order ODE, `gaudinODE`, with singular points $w_1,w_2$ and polynomial solutions $Q_k$ of degrees $\mathrm{degQ}_k=(S+P+3)-(\mu_k+4-k)$.
Its claim (polynomial solutions of these degrees and only these, for all 28 irreps of $(S,P)=(1,1),(2,1),(2,2),(3,2),(3,3)$) was reproduced in 56 of 56 tests.
We find that

$$
\mathrm{gaudinODE}=8\,(z-w_1)^2(z-w_2)^2\;\Pi^{p}\circ D^{\dagger}\circ\Pi^{-p},\qquad\Pi^{p}=(z-w_1)^{S}(z-w_2)^{P}, \tag{15}
$$

with $D^\dagger=\partial^4+\partial^3s_1+\partial^2s_2+\partial s_3+s_4$ the formal adjoint of our oper: an adjoint plus a gauge.
Matching coefficient by coefficient is linear in $(v_3,v_4^a,v_4^b,T_2(\mu),T_3(\mu),T_4(\mu))$; among the 16 (for $(1,1)$) and 50 (for $(2,1)$) gauges $p$ and site assignments scanned exactly one is consistent, $p=(S,P)$ with $S$ at $w_1=\theta_1$, and its unique solution equals the exact triple of the full operators and the central characters on all 28 irreps; with these substituted the two operators coincide as exact rational functions of $z$ (28 of 28).
With $N=S+P$ and $\mathrm{Cas2}=C_2(\mu)$, the symbolic solution is

$$
\begin{aligned}
v_3&=\tfrac12(S+1)\big(2N(N+2)-\mathrm{Cas2}\big),\qquad v_4^a=0,\\
v_4^b&=T_4(\mu)=\tfrac18\big(\mathrm{Cas2}^2+2(N+1)\mathrm{Cas2}-2\,\mathrm{Cas4}+4N(N+2)(N+3)\big),
\end{aligned} \tag{16}
$$

and $T_2(\mu)=\tfrac12(4N^2+6N-\mathrm{Cas2})$, $T_3(\mu)=\tfrac12(N+2)(2N(N+2)-\mathrm{Cas2})$, with $\mathrm{Cas4}$ the notebook's quartic Casimir.
The notebook's degrees are ours reversed, $\mathrm{degQ}_k=d_{5-k}$, and its solutions are built from ours: $Q_j=W_3[P_k,\,k\ne j]\big/\big((z-w_1)^S(z-w_2)^P\big)$, the $3\times3$ Wronskians of three of our four polynomial solutions, are polynomials of degree $\mathrm{degQ}_j$ and solve `gaudinODE` exactly (112 of 112).
The gauge power $p=(S,P)$ is the gl(4) version of $p=\lambda_1-\lambda_3$ of gl(3).
This is verified on the five $(S,P)$ listed, not proved for all.[^32]

## 7. Structure on the GT lattice

**Punchline.** For gl(3) the weight space is a chain and the reduced vertex operator a Jacobi matrix.
For gl(4) it is a three-parameter region of the GT lattice on which the three reduced operators share one 21-point stencil and, on 12 of the 20 off-diagonal classes, one product formula, but are not Jacobi in general.
The physical part of the spectrum on a weight space is cut out by an exact Weyl-type alternating sum over the translates $w\cdot\lambda'$ (832 sectors, 0 failures), and the components of every copy are polynomials in the triple of the degree that counting allows; all of it is verification on tested cases, not proof.[^33]

### 7.1. Stencil and closed forms

A GT pattern of gl(4) has rows $(m_{11})$, $(m_{12},m_{22})$, $(m_{13},m_{23},m_{33})$ and $\lambda$; a weight $\nu$ fixes the three row sums and leaves three free coordinates, so $W=V_\lambda[\nu]$ is a region $R_\nu$ of the GT lattice with $w=|R_\nu|$ points (for gl(3) one coordinate, a chain).
*Stencil.* All three reduced operators have the same stencil of 21 displacements $(\Delta m_{12},\Delta m_{13},\Delta m_{23})\in\{-1,0,1\}\times\big(\{(0,0)\}\cup\{(\pm1,0),(0,\pm1),\pm(1,-1)\}\big)$: row 2 moves by 0 or one root of gl(2), row 3 by 0 or one root of gl(3).
The words of $V_4^a,V_4^b$ alone reach 35 displacements; the 14 extra classes with $|\Delta m_{12}|=2$ cancel exactly after summation, a Pauli-type exclusion.
Every non-zero entry of the exact matrices on all weight spaces of seven representations, at $\lambda'=\lambda$ and at generic rational non-dominant $\lambda'$, lies in the stencil.

*Closed forms.* The matrix elements were derived symbolically from the GT formulas through the universal words (all intermediate states kept, exact line reconstruction).
Twelve of the 20 off-diagonal classes, those with $\Delta m_{12}=\pm1$ and a non-zero row-3 move, are given by one product formula for all three operators.
With $y_{ik}=m_{ik}-i$, $Y_s=\lambda_s-s$, $q=\lambda'$, let $y_a,y_b$ be the raised and lowered entry of row 2 and $y_i,y_j,y_k$ the raised, lowered and spectator entry of row 3:

$$
\begin{aligned}
M^{V_3}&=\frac{(y_{11}-y_b)(y_{11}-y_a+q_1-q_2)(y_a-y_j)(y_a-y_j+1)(y_a-y_k)\prod_{s=1}^4(y_i-Y_s)}{(y_a-y_b)(y_a-y_b+1)(y_i-y_j)(y_i-y_j+1)(y_j-y_k)(y_i-y_k)},\\
M^{V_4^a}&=(y_k+3)\,M^{V_3},\qquad M^{V_4^b}=-(y_a-y_i+2y_k-q_3+5)\,M^{V_3}.
\end{aligned} \tag{17}
$$

The row-1 and row-2 factors are those of the gl(3) chain with the twist $q_1-q_2$, the only dependence on $\lambda'$ besides $q_3$ in $V_4^b$; the raised row-3 entry carries the four GT factors $\prod_s(y_i-Y_s)$.
The formula was read off the machine-identified factor lists of the 12 classes and is supported by 17 730 exact held-out tests (seven representations, 16 253 non-zero entries, 0 mismatches); it is verified, not proved.
The other 8 off-diagonal classes (six row-3 moves with $\Delta m_{12}=0$, the pure row-2 moves $(\pm1,0,0)$) and the diagonal are linear factors times an irreducible polynomial $Q$ (degrees 3, 6, 11 for $V_3$); the six cubics of $V_3$ were reconstructed explicitly (46–50 terms), and $V_4^a$ has the same $Q$ on them (compared term by term on one class).
The off-diagonal elements of $V_4^a$ are an affine multiple of those of $V_3$ in 18 of the 20 classes (all but $(\pm1,0,0)$); all forms agree with G3's exact matrices on held-out representations (0 mismatches, 0 uncovered entries).[^34]

### 7.2. The Weyl/BGG identity

By (12) a joint eigenvector on $W$ that violates a string condition is carried to a joint eigenvector at $s_i\cdot\lambda'$ with the same triple; the data say more, that the physical spectrum is the Euler characteristic over the whole Weyl orbit.
With $w\cdot\lambda'=w(\lambda'+\rho)-\rho$, $\rho=(3,2,1,0)$, and $V_\lambda[\mu-w\cdot\lambda']=0$ when $\mu-w\cdot\lambda'$ is not a weight of $V_\lambda$, for every joint triple $\tau$

$$
\sum_{w\in S_4}\mathrm{sgn}(w)\;\mathrm{mult}_\tau\big(V_{\mathrm{red}}(w\cdot\lambda')\ \text{on}\ V_\lambda[\mu-w\cdot\lambda']\big)=\begin{cases}1,&\tau\in\mathrm{Spec}(M_\mu),\\0,&\text{otherwise.}\end{cases} \tag{18}
$$

This holds exactly on 832 sectors of 22 products (all dominant $\mu$ with a non-empty translate space; $\sum m=1944$ physical copies, 3408 distinct joint triples matched at $10^{-30}$ from 50-digit spectra), with 0 failures: the six products of the programme (all 76 sectors of $[6,4,2,0]^{\otimes2}$, $m\le12$), the held-out $[3,2,1,0]\otimes[4,2,1,0]$ in both orders, $[4,2,2,0]^{\otimes2}$, $[4,3,1,0]^{\otimes2}$, and twelve more (including $\lambda\ne\lambda'$ and $\lambda_4>0$) up to $[8,4,0,0]^{\otimes2}$ (131 sectors).
The test can fail: with all signs $+1$ it holds on only 11 of 16, 14 of 27 and 32 of 76 sectors of $[3,2,1,0]^{\otimes2}$, $[4,2,1,0]^{\otimes2}$, $[6,4,2,0]^{\otimes2}$, and with the plain Weyl action $w(\lambda')$ (no $\rho$-shift) on 1 of 16, 2 of 27 and 8 of 76.
The joint spectrum is simple on every translate space, also at non-dominant $w\cdot\lambda'$ (closest distinct triples at least 4.0); taking dimensions gives the Kostant–Racah–Speiser formula $m=\sum_w\mathrm{sgn}(w)\dim V_\lambda[\mu-w\cdot\lambda']$, true on all 832 sectors (a check of the dot action).

*One step and depth.* Let $S(\tau)$ be the set of strings violated by the eigenvector of $\tau$ on $W$ and $T(\tau)$ the set of $i$ for which $\tau$ lies in the spectrum at $s_i\cdot\lambda'$.
On all 3408 triples $S(\tau)=T(\tau)$ ((12) gives $\subseteq$, the converse is an observation), and $\tau$ is physical iff $S(\tau)=\emptyset$.
Of the 1464 spurious triples ($w>m$), 1383 violate one string and 81 violate two; each of the latter lies in exactly the four translate spaces $\{e,s_i,s_j,w_{ij}\}$, $w_{ij}$ of length 2, and is cancelled by inclusion–exclusion ($1-2+1$).
So the one-step union is not disjoint, but no triple needs a translate of length $\ge3$: no such translate space is non-empty on any of the 832 sectors (the longest non-empty one has length 1 in ten products and 2 in twelve); the minimal depth is 0, 1, 2 on 329, 438, 65 sectors.
*gl(3).* The same identity ($S_3$, $\rho=(2,1,0)$, $V=tS_{3,2}$) holds on all 227 sectors of seven products (among them $[14,7,0]^{\otimes2}$ with 113 sectors and 344 copies), with $S=T$, $|S|\le1$ and depth $\le1$.[^35]

*Status.* Identity (18) is *argued*, not proved.
If $V_{\mathrm{red}}(w\cdot\lambda')$ on $V_\lambda[\mu-w\cdot\lambda']$ are the Gaudin operators on the singular vectors of $V_\lambda\otimes M(w\cdot\lambda')$ (the interpretation of Section 5.4), then decomposition by joint eigenvalue is an exact functor, because the Gaudin algebra commutes with gl(4), and $V_\lambda\otimes(\text{BGG resolution of }L(\lambda'))$ resolves the $\tau$-part of $V_\lambda\otimes L(\lambda')$ by the $\tau$-parts of $V_\lambda\otimes M(w\cdot\lambda')$ with signs $(-1)^{\ell(w)}$.
(18) says that taking singular vectors of weight $\mu$ (only left exact) preserves this Euler characteristic $\tau$ by $\tau$, as it does for dimensions: that last step is what the data confirm and what is not proved.

### 7.3. Which eigenlines are physical

For gl(3) the physical eigenvalues of $V$ on $W$ are the $m$ largest when only $E_{23}$ is active and the $m$ smallest when only $E_{12}$ is active (all 70 sectors with $w>m$ of seven products, no exception); the “$m$ largest” of `PY3/RESULTS-jacobi-cg.md` §4 is the $E_{23}$ half.
For gl(4) no single ordering works: by $V_3$ the physical set is the top $m$ in 20 and the bottom $m$ in 9 of the 78 training sectors with $w>m$ (six products), similarly for $V_4^a$, $V_4^b$ and G2's generic combination, and no linear functional $aV_3+bV_4^a+cV_4^b$ works on all 78.
Sorted by the active strings $A=\{i:E_i^{k_i}\ \text{non-void on}\ W\}$ (intrinsic), a linear functional with physical $=$ top $m$ exists for $A=\{1\}$ and for $A=\{3\}$ (linear programme on the training sectors; exact on the held-out sectors, 5 of 5 and 7 of 7): $\ell_1\approx-V_3+0.27V_4^a+0.10V_4^b$ and $\ell_3\approx V_4^a-0.49V_3+0.33V_4^b$, maximum-margin points, not canonical.
For the middle string $A=\{2\}$ (21 training sectors) and every mixed class none exists.[^36]

### 7.4. Chain structure and polynomiality

*The graph.* $V_{\mathrm{red},3}$ is not a Jacobi matrix in general: on the 150 weight spaces of dimension $\ge2$ of nine products its sparsity graph is connected in all 150, but a path (tridiagonal in a suitable order) in only 29 of the 44 spaces with $m=w$ and 19 of the 106 with $m<w$; otherwise it has cycles (the smallest example, $(5,4,2,1)$ of $[3,2,1,0]^{\otimes2}$ with $m=w=4$, is the complete graph $K_4$).
For $m=w$ the zero pattern is symmetric; for 43 of the 106 spaces with $m<w$ it is asymmetric, the one-sided cut that is the analogue of the $E_{12}$ cut of gl(3), where all 38 spaces tested are paths.

*Polynomiality.* Some coordinate of $M_\mu$ is cyclic for the triple in all 181 sectors with $m\ge2$, so every component of every copy is a polynomial in $(v_3,v_4^a,v_4^b)$ with rational coefficients, exactly (reproducing G3's 50-digit copies to $7\times10^{-46}$).
The minimal total degree equals the dimension-counting minimum $d_{\rm gen}(m)=\min\{d:\binom{d+3}{3}\ge m\}$ in 162 sectors and $d_{\rm gen}+1$ in 19, never less: affine in 122 of the 138 sectors with $m\le4$, quadratic for $5\le m\le10$ (three sectors need cubic), cubic for $m=11,12$, against degree $m-1$ for a single operator (for $m=12$: 3 against 11).
So $\mathrm{slice}_k=p_k(V)$ of the gl(3) chain becomes a polynomial in the triple of the degree of $m$ generic points; on the whole weight space a cyclic lattice point exists in 150 of 150 spaces.[^37]

## 8. The benchmark and its caveats

**Punchline.** On the Munich family $[S,S,0,0]^{\otimes2}$ the engine is about 860 times faster than the Munich code at $S=8$, where the latter takes 80 s, and it continues to $S=20$; on multiplicity families it resolves up to 51 copies per sector in seconds.
Every number has a caveat: two machines, the Munich code as distributed, Python overhead for $S\le2$, exact against float arithmetic.

### 8.1. The multiplicity-free family

Table 4 puts the Munich sweep next to the engine on $[S,S,0,0]^{\otimes2}$.
Both produce the same object, the HW vector with all $n_W$ coefficients of every irrep, and the family is multiplicity-free, so no operator enters.
$t_a$ is the part that never touches $n_W$ (exact GT generators, Racah–Speiser, exact PRV spaces); $t_b$ adds the metric, the float matrices and left inverses, and the float64 reconstruction of every full vector.
Engine: i7-12700KF, minimum of five cold runs, the quietest of three blocks of the whole table (run-to-run differences at most 3% between two re-runs and 15% for the first block, during which another job was running); Munich: i5-1345U, one run per irrep.[^38]

**Table 4.** $[S,S,0,0]^{\otimes2}$: highest-weight solve of all irreps. Munich (i5-1345U, as distributed, hw-solve only) against the engine (i7-12700KF); ratio = Munich${}/t_b$. `n/a`: not run.

| $S$ | $\dim V$ | irreps | max $w$ | max $n_W$ | Munich [s] | $t_a$ [s] | $t_b$ [s] | ratio | worst resid. |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 6 | 3 | 1 | 6 | 0.0000279 | 0.0002 | 0.0007 | 0.04 | $1.9\times10^{-16}$ |
| 2 | 20 | 6 | 2 | 22 | 0.000858 | 0.0007 | 0.0018 | 0.48 | $2.9\times10^{-16}$ |
| 3 | 50 | 10 | 2 | 62 | 0.0107 | 0.0016 | 0.0040 | 2.7 | $5.9\times10^{-16}$ |
| 4 | 105 | 15 | 3 | 147 | 0.0574 | 0.0034 | 0.0088 | 6.5 | $8.2\times10^{-16}$ |
| 5 | 196 | 21 | 3 | 308 | 0.355 | 0.0064 | 0.0179 | 20 | $1.8\times10^{-15}$ |
| 6 | 336 | 28 | 4 | 588 | 2.49 | 0.0110 | 0.0326 | 76 | $1.8\times10^{-15}$ |
| 7 | 540 | 36 | 4 | 1044 | 15.3 | 0.0171 | 0.0572 | 268 | $1.3\times10^{-15}$ |
| 8 | 825 | 45 | 5 | 1749 | 80.2 | 0.0258 | 0.0932 | 860 | $5.6\times10^{-15}$ |
| 10 | 1716 | 66 | 6 | 4290 | n/a | 0.0579 | 0.228 | – | $2.8\times10^{-15}$ |
| 12 | 3185 | 91 | 7 | 9191 | n/a | 0.103 | 0.443 | – | $1.3\times10^{-14}$ |
| 16 | 8721 | 153 | 9 | 31 977 | n/a | 0.322 | 1.62 | – | $2.6\times10^{-13}$ |
| 20 | 19 481 | 231 | 11 | 86 779 | n/a | 0.824 | 4.59 | – | $3.9\times10^{-11}$ |

Munich's total grows by a factor 5.2–7.0 per unit of $S$ between $S=4$ and $S=8$; the engine's $t_b$ by 2.0 from $S=4$ to 5, decreasing to 1.4 from $S=11$ to 12, in proportion to the number of coefficients written ($\propto S^{5.9}$: 580 at $S=4$, 52 090 at $S=10$, $3.06\times10^6$ at $S=20$).
At $S=8$ the singlet, Munich's slowest row (18.7 s), takes 4.5 ms once the per-$\lambda'$ setup of 16 ms has been paid and 45 ms cold; per-irrep speed-ups run from 101 (the top irrep) to 4100 (the singlet), and the sum over the 45 irreps is 80.19 s against 0.052 s (warm, setup excluded).
For $S\le2$ the Munich C++ is faster in absolute terms (28 $\mu$s against 0.7 ms at $S=1$, 0.86 ms against 1.8 ms at $S=2$): the Python overhead regime, with the crossover at $S=3$.

### 8.2. Multiplicity families

For $m\ge2$ the engine also resolves the multiplicities, so $t_a$ now contains the reduction (4 ms per $\lambda'$), the exact application of the three reduced operators to the PRV basis, and the joint eigen-decomposition; $t_b$ adds the recursion setup and the reconstruction of every full copy.
Table 5 gives the scaling family $[2a,a,0,0]^{\otimes2}$, the gl(4) analogue of $[2a,a,0]$ ($\max m=a+1$), and the stress family $[3a,2a,a,0]^{\otimes2}$ ($\max m=4,12,27,51$, one irrep each time).
The Munich SU(4) data are multiplicity-free, so no comparison exists; its null-space method would need an $n_W\times n_W$ null space per irrep with $n_W$ up to 29 028 at $a=6$, and it was not run.

**Table 5.** Canonical copies on multiplicity families (i7-12700KF; minimum of 5/3/2 cold runs, one run for the rows marked $\ast$; the rows $a=1,3,5,7$ of the first family are in `GL/benchmark_gl4.md`). $t_a$ with a 50-digit eigen-decomposition and $t_a^f$ with float64; $t_b^f$: float64 eigen-decomposition and all full vectors written (n/a: slices only).

| family, $a$ | $\dim V$ | irreps | copies | max $m$ | max $w$ | max $n_W$ | $t_a$ [s] | $t_a^{f}$ [s] | $t_b^{f}$ [s] |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| $[2a,a,0,0]^{\otimes2}$, 2 | 126 | 24 | 32 | 3 | 4 | 174 | 0.041 | 0.035 | 0.041 |
| 4 | 1375 | 115 | 205 | 5 | 11 | 3541 | 0.518 | 0.432 | 0.525 |
| 6 | 6860 | 322 | 728 | 7 | 20 | 29 028 | 3.16 | 2.61 | 3.34 |
| 8$^\ast$ | 23 085 | 693 | 1905 | 9 | 33 | 146 025 | 12.0 | 9.89 | n/a |
| $[3a,2a,a,0]^{\otimes2}$, 1 | 64 | 16 | 30 | 4 | 4 | 152 | 0.049 | 0.035 | 0.038 |
| 2 | 729 | 76 | 297 | 12 | 15 | 4743 | 1.47 | 0.709 | 0.757 |
| 3 | 4096 | 213 | 1592 | 27 | 32 | 59 008 | 23.4 | 6.42 | 6.91 |
| 4$^\ast$ | 15 625 | 459 | 5945 | 51 | 65 | 426 425 | 258 | 37.2 | n/a |

*Where the time goes.*
The reduction is negligible.
At $[12,6,0,0]^{\otimes2}$ ($a=6$, 728 copies) $t_a=3.16$ s is generators 0.65 s, Racah–Speiser and PRV spaces 0.10 s, *exact application of the reduced operators 1.84 s* and the 50-digit eigen-decomposition 0.56 s; $t_b=3.89$ s adds 0.09 s of setup and 0.64 s for all 728 vectors ($2.2\times10^6$ coefficients).
With exact Python fractions the three reduced operators thus cost about three times the level recursion: the $n_W$-free part is cheap relative to Munich's $n_W^3$ but not relative to the recursion, and a float or modular application was not written.
For large multiplicity the 50-digit eigen-decomposition dominates (17.0 of 23.4 s at $[9,6,3,0]^{\otimes2}$, 221 of 258 s at $[12,8,4,0]^{\otimes2}$); the float64 version is 40 to 1600 times cheaper and accurate to $2.1\times10^{-15}$ in the eigenvalues.
Per operator and copy the gl(4) reduced operators cost 678 $\mu$s at $[10,5,0,0]^{\otimes2}$ ($w\le15$), 2.6 times the 24-word gl(3) operator at the same $w$.[^39]

### 8.3. Caveats

*Two machines, one code state.*
Munich: 13th-gen i5-1345U (Ubuntu, g++ 11.4, LAPACK 3.10), one run per irrep, the code as distributed; engine: i7-12700KF (Windows 11, Python 3.13, numpy 2.5.2, scipy 1.18.1), single process, minimum of five cold runs.
The single-thread speeds of the two CPUs differ by a factor of order 1–2, so the ratios of Table 4 are orders of magnitude, not a measured speed-up.
For gl(3) the repair of the same routine bought a factor 19 (gl(3) notes, Table 2); for SU(4) no such run exists, so the 860 at $S=8$ is against the code as distributed, and a repair would reduce it by an amount that was not measured.

*What is compared, and precision.*
Munich's routine returns the null vector of each irrep in an orthonormal GT basis, up to sign; ours is the plain GT basis, converted by $\mathrm{CGC}=\pm\Psi\sqrt{g_\nu g_q/\mathrm{GNORM2}}$.
Our vectors are float64 with HW residual at most $1.3\times10^{-14}$ up to $S=12$, $2.6\times10^{-13}$ at $S=16$ and $3.9\times10^{-11}$ at $S=20$: beyond $S\approx16$ the double-precision recursion degrades by a factor 2.4–3.7 per unit of $S$, while the exact recursion rounded at the end keeps $10^{-15}$ at 1–16 s per vector (Python fractions, not vectorised).[^40]

*What was not run.*
The literal comparison of our HW vectors with the Munich ones: the binary cannot be built on the Windows machine.
The two GT bases do coincide (pattern order equal to Munich's index $Q(M)$ on 17 irreps, 6024 lowering matrix elements equal to Munich's `lowering_coeff`$^2$ exactly), so the comparison should be a formality up to one global sign; G1 §9 gives the 19 `hwcgc.sh` commands for $S\le3$ and a tested comparer.

*The gl(3) context.*
At $[14,7,0]^{\otimes2}$ (113 irreps, 344 copies) the gl(3) slice engine takes 0.10 s ($t_a$) and 0.17 s ($t_b$) on the i7, against 206 s for the Munich code as distributed and 10.9 s repaired on the i5 (gl(3) notes, Table 2).
The SU(3)-specific ndsu3lib [20], on the same i5 laptop, needs 0.80 ms for the HW row only (0.116 s for all reduced coefficients of all states): it is faster than the slice engine at $a=7$ for the HW data, and the Jacobi-chain route overtakes it at about $a=15$.[^41]
So “orders of magnitude” holds against general-$N$ null-space codes only; ndsu3lib is specific to SU(3) and has no counterpart for $N=4$.

## 9. What is classical, what is new, what is next

### 9.1. Status of the statements

Table 6 lists each statement with the way it is established.
“Verified” means exact or high-precision agreement on the tested cases, not a proof for all $\lambda,\lambda',\mu$.

**Table 6.** Status of the statements. *Theorem*: classical or proved here. *Exact identity*: a finite computation in $U(\mathrm{gl}_4)^{\otimes2}$. *Verified*: tested cases only.

| statement | status | Section |
| :--- | :--- | ---: |
| slice theorem with $N-1$ string conditions, $\dim M_\mu=m_\mu$ | classical theorem [2, 3] | 3.2 |
| level recursion fixes $\Psi$ from the slice; $A_b$ injective | proved; implementation verified (26 products, exact twin) | 3.1, 3.3 |
| the dangling conditions are needed in the exact twin for $N=4$ | verified (253 of 1655 slices) | 3.3 |
| Laurent coefficients of the cdet commute; $S_{k,j}\propto t^{\,j-k}$, so $V_3,V_4^a,V_4^b$ are $t$-independent | commutation: theorem [14], term lists verified; homogeneity proved | 4.1, 4.2 |
| identities (6)–(8); naive traces fail | exact; scalar remainders verified (81 sectors) | 4.2 |
| simple joint spectrum of the triple | theorem for the Gaudin algebra [8, 9, 10, 11]; verified: 180 sectors exact, 1864 at 50 digits | 4.3 |
| $V_3$ alone degenerate (43 of 180); $\{V_3,V_4^b\}$ separates; $\{V_3,V_4^a\}$ fails on 2 | exact discriminants; tested products | 4.3 |
| $PV_3P=-V_3+c$; forced degeneracy for $\lvert m_+-m_-\rvert\ge2$ | from the dictionary; forced/accidental split observed | 4.3 |
| reduction rule; universality in $\lambda'$; $M_\mu$ preserved | proved | 5.1 |
| direct route agrees with the reduced operators; closed forms (10), (11) | exact (46 + 598 irreps, 2275 $\theta_2$ checks; 898 pairs); closed forms follow from identities and rule | 5.2, 5.3 |
| BGG intertwining; commutation on $W$ for all $\lambda'$; $M_\mu$ spectral | verified (5814 tests, nine products); not proved | 5.4 |
| Eq. (13) on $[S,S,0,0]^{\otimes2}$ | verified, $S\le10$ (285 irreps) | 5.5 |
| the triple generates the 18 Laurent coefficients mod scalars; dictionary of Table 3, $c=\tau_3$, $c'=\tau_4$ | exact identity; dictionary proved given $T_k$ (70 partitions), verified on 95 sectors | 6.1 |
| polynomial solutions of degrees $\mu_j+4-j$, Wronskian, exponents; $\mu_4=0\Rightarrow V_4^a=V_4^b=0$ | verified on 210 copies; $\mu_4=0$ on 47 sectors, follows from the former | 6.2, 6.4 |
| $\Sigma$ = joint spectrum; all four degree conditions needed | exact certificate per sector (180 sectors); pencil analysis on 68 | 6.3 |
| `gaudinODE` is adjoint plus gauge of the oper | verified, 28 irreps of five $(S,P)$ | 6.5 |
| speed against the Munich code | measurement, two machines, code as distributed | 8 |
| Weyl/BGG identity (18); $S=T$; translates of length $\le2$ suffice; signed lattice count $=m$ | verified: 832 sectors, 22 products, 0 failures (gl(3): 227); argued, one step unproved | 7.2 |
| 21-point stencil; product formula (17) on 12 of 20 off-diagonal classes; the rest: linear factors $\times$ irreducible $Q$ | derived symbolically; 17 730 exact held-out tests; not proved | 7.1 |
| no single ordering of the triple picks the physical eigenlines; linear functionals for $A=\{1\},\{3\}$, none for $\{2\}$; gl(3) top/bottom rule; $V_{\mathrm{red},3}$ not Jacobi in general | verified (78 + 28 sectors; 70 gl(3) sectors; 150 weight spaces) | 7.3, 7.4 |
| components of every copy are polynomials in the triple, of the counting degree | verified, exact (181 sectors with $m\ge2$) | 7.4 |

### 9.2. Classical and new

*Classical.*
The slice criterion is the lemma of Kostant and the theorem of Parthasarathy, Ranga Rao and Varadarajan [2, 3] (statement in Kumar [4]; Zhelobenko [5] in the physics literature); related reductions have precedents as algorithms for SU(3): the extremal-projector method of Asherova, Smirnov and Tolstoy [6] and, in the nuclear-physics literature, recursions from extremal coefficients with a canonical zero pattern (Draayer–Akiyama and its descendants [19, 20]).
That the Gaudin (argument-shift, Bethe) algebra resolves multiplicities, with simple joint spectrum on every multiplicity space, is known [7, 8, 9, 10, 11]; its eigenvalues are encoded by opers [12, 13]; the quantum characteristic polynomial is Talalaev's [14, 15].
For $N=3$ the operator is the missing-label operator of Pluhař, Smirnov and Tolstoy and of Crampé, Poulain d'Andecy and Vinet [16, 17], who also identify it with the two-site Gaudin model (with Campoamor-Stursberg and Musso [18]) and propose computing Clebsch–Gordan coefficients from Gaudin correlation functions.

*New, as far as we know.*
The explicit gl(4) engine with $N-1$ string conditions and a level recursion whose exact twin needs the dangling conditions, with the SU(4) benchmark;
the exact triple $V_3,V_4^a,V_4^b$ with its $t$-independence and structure identities, the failure of the naive trace invariants (the quantum correction has weight $\tfrac12$) and of the cubic operator alone, with the exchange-parity argument;
the reduction rule with normal ordering, the universal polynomials, the canonical copies and the intertwining facts;
the Weyl-type identity on a weight space, the common stencil and product formula of the three operators, and the polynomiality of the copies in the triple;
the forced dictionary with $T_k$, $\tau_k$, the exact Q-system solve (all four degree conditions needed), $\mu_4=0\Rightarrow V_4^a=V_4^b=0$, and the July-2026 gl(4) ODE as the adjoint of the oper up to a gauge.
The literature check of 21 October concerned gl(3) and the SU(4) literature was not searched separately: “new” is relative to that check and the references cited here.

### 9.3. Next

1. A proof of the Weyl identity (18): the missing step is that taking singular vectors of weight $\mu$ preserves the Euler characteristic joint eigenvalue by joint eigenvalue (for instance through the Verma flag and the vanishing of $\mathrm{Ext}^1$ between the relevant Verma modules for dominant $\mu$); likewise for $S=T$, the exclusion of $|\Delta m_{12}|=2$ and the closed forms of Section 7.1.
2. What Section 7 leaves open: explicit $Q$ for the remaining classes and the diagonal, a rule for the middle string and the mixed classes (and the full cone of separating functionals), and the Weyl test beyond $[8,4,0,0]^{\otimes2}$ ($[9,6,3,0]^{\otimes2}$ was not run).
   Further proofs: commutation of the reduced operators on the whole weight space for every $\lambda'$; the intertwining relation; simple spectrum of the triple beyond the cited theorem; the closed forms (13).
3. The literal comparison with the Munich HW vectors on the machine where the binary builds ($S\le3$ first), a same-machine Munich SU(4) run including a repaired variant, and a search of the SU(4) Clebsch–Gordan literature.
4. A float or modular application of the reduced operators (10–100 times faster than exact fractions by estimate; not written), and a compiled version.
5. gl(N), $N\ge5$: the slice engine runs for $N=5$ on nine products as a smoke test and the coefficient $(N-2)/2$ in (6) is solved for $N\le6$, but the $(N-1)(N-2)/2$ missing-label operators of higher $N$ were not constructed.
6. Where this material goes: into the gl(3) slice-engine notes, or into the paper.

## A. Conventions

- *One-site basis.* GT patterns of $\mathrm{gl}(N)$ with rows of length $1,\dots,N$, the last row equal to $\lambda$, in the order of `gl3_reps.gt_patterns($\lambda$)`: ascending lexicographic in (row $N-1$, …, row 1), highest-weight state *last*; Python indices are 0-based, the reference files 1-based (the Wolfram notebooks use another order; the map is by pattern).
- *Generators and metric.* Plain rational, non-unitary GT matrices acting on kets as columns for $E_{i,i+1}$, $E_{i+1,i}$ and the diagonal ones; the others by commutators, $E_{ij}=[E_{i,j-1},E_{j-1,j}]$, $E_{ji}=[E_{j,j-1},E_{j-1,i}]$ for $j>i+1$ (path-independent; all 256 relations hold exactly on nine irreps).
  Weights are the eigenvalues of $(E_{11},\dots,E_{44})$, $\alpha_i=\epsilon_i-\epsilon_{i+1}$.
  The invariant metric is the diagonal $G$ with $E_{ij}^TG=GE_{ji}$, $G[\mathrm{HW}]=1$; the float recursion runs in the orthonormal basis $\tilde e_q=e_q/\sqrt{g_q}$, which is Munich's (pattern index $Q(M)$ = our index $+1$; its SU(4) label is $\mu-\mu_4(1,1,1,1)$).
- *Two sites, slices, copies.* $\Psi=\sum\Psi[\nu,q]\,e_\nu\otimes e_q$, site 1 $=V_\lambda$ (index $\nu$, slow), site 2 $=V_{\lambda'}$ ($q$, fast); the slice is $\Psi[\,\cdot\,,q_{\rm top}]\in V_\lambda[\mu-\lambda']$.
  The basis of $M_\mu$ is the reduced row echelon form of the exact null space; copies are ordered lexicographically by $(V_3,V_4^a,V_4^b)$ and normalised so that the first non-zero slice component is 1.
  Multiplicities by Racah–Speiser with $\rho=(3,2,1,0)$ over the weights of the second factor; the PRV test runs over the union of the candidate sets of both factors, including those of multiplicity 0.
- *Term lists* (`GD/ref/gl4-gaudin-ops.txt`): a line `OP` with the operator name, then lines `TERM`: an exact Laurent polynomial in $t$ (list of `{power,p/q}`) and a word (list of `{site,i,j}`, leftmost factor acting last, site-1 factors first, each site's internal order that of the cdet).

## B. The reduced universal polynomial of $V_3$

Table 7 and the list below give the 39 words of $V_{\mathrm{red},3}$, the reduction of the 88 terms of $V_3$ (word $=$ product of site-1 generators $E_{ij}$, leftmost acting last; $\lambda'_i$ written $l_i$): three words of one generator, twelve of two and twenty-four of three.
For $\lambda'_3=\lambda'_4=0$ the three words with coefficient $l_3+l_4$ drop out (36 words).
The closed form (10) is a different presentation of the same operator (Section 5.2).

**Table 7.** Universal $V_{\mathrm{red},3}$: the 15 words of degree $\le2$ with their coefficients (polynomials in $l_i=\lambda'_i$); the 24 cubic words are listed below.

| word | coefficient | word | coefficient | word | coefficient |
| :--- | :--- | :--- | :--- | :--- | :--- |
| $E_{22}$ | $l_3+l_4$ | $E_{11}E_{44}$ | $l_2+l_3$ | $E_{32}E_{23}$ | $1-l_1-l_4$ |
| $E_{33}$ | $l_1+l_2+2l_4$ | $E_{21}E_{12}$ | $-(l_3+l_4)$ | $E_{33}E_{44}$ | $l_1+l_2$ |
| $E_{44}$ | $2(l_1+l_2+l_3)$ | $E_{22}E_{33}$ | $l_1+l_4$ | $E_{41}E_{14}$ | $-(l_2+l_3)$ |
| $E_{11}E_{22}$ | $l_3+l_4$ | $E_{22}E_{44}$ | $l_1+l_3$ | $E_{42}E_{24}$ | $1-l_1-l_3$ |
| $E_{11}E_{33}$ | $l_2+l_4$ | $E_{31}E_{13}$ | $-(l_2+l_4)$ | $E_{43}E_{34}$ | $2-l_1-l_2$ |


The 24 cubic words all have coefficient $\pm1$.
$+1$: $E_{11}E_{32}E_{23}$, $E_{11}E_{42}E_{24}$, $E_{11}E_{43}E_{34}$, $E_{21}E_{33}E_{12}$, $E_{21}E_{44}E_{12}$, $E_{22}E_{43}E_{34}$, $E_{31}E_{22}E_{13}$, $E_{31}E_{44}E_{13}$, $E_{32}E_{44}E_{23}$, $E_{41}E_{22}E_{14}$, $E_{41}E_{33}E_{14}$, $E_{42}E_{33}E_{24}$.
$-1$: $E_{21}E_{32}E_{13}$, $E_{21}E_{42}E_{14}$, $E_{31}E_{12}E_{23}$, $E_{31}E_{23}E_{12}$, $E_{31}E_{43}E_{14}$, $E_{32}E_{43}E_{24}$, $E_{41}E_{12}E_{24}$, $E_{41}E_{13}E_{34}$, $E_{41}E_{24}E_{12}$, $E_{41}E_{34}E_{13}$, $E_{42}E_{23}E_{34}$, $E_{42}E_{34}E_{23}$.

## C. Reproduction

All code is Python 3.13 (numpy, scipy, sympy, mpmath) in `GL/`, run from that directory with `python -X utf8`; each module has a CLI and a few-second regression test.
Times are for the i7-12700KF; the output of one run of `all` is stored beside each module as `*_output.txt`.

- G1, `slice_engine_gl4.py` with `glN_reps.py`: `tests`, `prv`, `hw`, `exact`, `neg`, `brute`, `bench`, `families`, `ref`, `verify`, `munich`, `munichcmp`, `n5`, `extras`, `all` (7.4 min); `test_glN_reps.py`: 10 tests.
- G2, `gaudin_gl4_ops.py`: `sym`, `scaffold`, `gl3`, `comm`, `triple`, `struct`, `terms`, `big`, `all` (301 s; `big` is 62, 44 and 111 s for $[5,3,1,0]^{\otimes2}$, $[5,3,2,0]^{\otimes2}$, $[6,4,2,0]^{\otimes2}$); 11 tests, 5 s.
- G3, `slice_ops_gl4.py`: `order`, `reduce`, `gl3`, `direct`, `copies`, `bench`, `ref`, `verify`, `all` (477 s with the option `ext`); 9 tests, 1.5 s.
- G5, `oper_gl4.py`: `dict`, `exps`, `solve`, `search`, `wb`, `big`, `all` (516 s: 11, 45, 175, 60, 170, 100 s), and `wbsym` (150 s, not in `all`); 7 tests, 5 s.
- G4, `structure_gl4.py`: `stencil`, `symbolic`, `symq`, `univ`, `weyl`, `rule`, `chain`, `poly`, `cycW`, `all` (`run`, `merge` assemble per-stage parts; `stencil` 7 s, `weyl` 45 s, `rule` 5 s, `chain` 1.4 s, `poly` 4.9 s, `cycW` 2.9 s, `univ` 7 s, `symq` 294 s, `symbolic` 1061 s); 9 tests, 7 s; all tables in `structure_gl4_tables.txt`.

The reference data, in `GD/ref/`, are exact machine-readable files with a documented header:

- `gl4-cg-slices-{2100,2110,3210,2100-1100,2210-2110,2200,3300,4400,6600,8800}.txt` (G1, 250 kB): GT bases, weights, metric; per sector $(\mu,m,w,n_W,\dim V_\mu,C_2(\mu))$, weight space, PRV basis, $m=1$ slices, residuals, exact $G$-norms and full vectors (for $\dim V_\lambda\dim V_{\lambda'}\le2000$); re-verified from the text with independent code (0 mismatches at $10^{-12}$).
- `gl4-gaudin-ops.txt` (G2, 87 kB): 12 operators $S_{k,j}$ at $\theta_1$, `V3`, `V4a`, `V4b` as term lists (loader `load_terms`).
- `gl4-cg-copies-{2100,2110,3210,2210-2110,4200,6420}.txt` (G3): the slice records plus the exact matrices of the reduced triple and the canonical copies to 30 digits (0 mismatches).
- `gl4-oper-dictionary.txt` (G5, 135 lines): 18 dictionary lines (Table 3), the Casimir formulas, and the joint triples of 146 copies to 25 digits.

## References

[1] A. Alex, M. Kalus, A. Huckleberry and J. von Delft, *A numerical algorithm for the explicit calculation of SU(N) and SL(N,C) Clebsch–Gordan coefficients*, J. Math. Phys. 52 (2011) 023507, arXiv:1009.0437.

[2] B. Kostant, *A formula for the multiplicity of a weight*, Trans. Amer. Math. Soc. 93 (1959) 53–73 (Lemma 4.1).

[3] K. R. Parthasarathy, R. Ranga Rao and V. S. Varadarajan, *Representations of complex semi-simple Lie groups and Lie algebras*, Ann. of Math. 85 (1967) 383–429.

[4] S. Kumar, *Tensor product decomposition*, Proc. ICM Hyderabad 2010, vol. IV (statement of Kostant's lemma).

[5] D. P. Zhelobenko, *Compact Lie Groups and Their Representations*, AMS Transl. Math. Monographs 40, 1973.

[6] R. M. Asherova, Yu. F. Smirnov and V. N. Tolstoy, Theor. Math. Phys. 8 (1971) 813–825 (extremal projectors).

[7] E. B. Vinberg, *On certain commutative subalgebras of a universal enveloping algebra*, Izv. Akad. Nauk SSSR Ser. Mat. 54 (1990) 3–25 [Math. USSR-Izv. 36 (1991) 1–22].

[8] L. G. Rybnikov, *Argument shift method and the Gaudin model*, Funct. Anal. Appl. 40 (2006) 188–199, arXiv:math/0606380; *A proof of the Gaudin Bethe Ansatz conjecture*, arXiv:1608.04625.

[9] B. Feigin, E. Frenkel and L. Rybnikov, Duke Math. J. 155 (2010) 337–363, arXiv:0712.1183.

[10] E. Mukhin, V. Tarasov and A. Varchenko, J. Amer. Math. Soc. 22 (2009) 909–940, arXiv:0711.4079.

[11] I. Halacheva, J. Kamnitzer, L. Rybnikov and A. Weekes, Duke Math. J. 169 (2020) 2337–2419, arXiv:1708.05105.

[12] B. Feigin, E. Frenkel and N. Reshetikhin, Commun. Math. Phys. 166 (1994) 27–62, hep-th/9402022.

[13] E. Frenkel, *Opers on the projective line, flag manifolds and Bethe Ansatz*, Mosc. Math. J. 4 (2004) 655–705, arXiv:math/0308269; *Gaudin model and opers*, Progr. Math. 237 (2005) 1–58, arXiv:math/0407524.

[14] D. Talalaev, *Quantization of the Gaudin system*, Funct. Anal. Appl. 40 (2006) 73–77, arXiv:hep-th/0404153.

[15] A. Chervov and D. Talalaev, *Quantum spectral curves, quantum integrable systems and the geometric Langlands correspondence*, arXiv:hep-th/0604128.

[16] Z. Pluhař, Yu. F. Smirnov and V. N. Tolstoy, *Clebsch–Gordan coefficients of SU(3) with simple symmetry properties*, J. Phys. A 19 (1986) 21–28.

[17] N. Crampé, L. Poulain d'Andecy and L. Vinet, *The missing label of $\mathrm{su}_3$ and its symmetry*, Commun. Math. Phys. 400 (2023) 179–213, arXiv:2110.03521; *A Calabi–Yau algebra with $E_6$ symmetry and the Clebsch–Gordan series of $\mathfrak{sl}_3$*, J. Lie Theory 31 (2021) 1085–1112, arXiv:2005.13444.

[18] R. Campoamor-Stursberg and F. Musso, *Two-body homogeneous rational Gaudin models and the missing label problem*, J. Phys. A 46 (2013) 335201.

[19] J. P. Draayer and Y. Akiyama, J. Math. Phys. 14 (1973) 1904–1912.

[20] Herko, Caprio, McCoy and Fasano (ndsu3lib), Eur. Phys. J. A 62 (2026) 95, arXiv:2505.08993.

Talalaev and Chervov–Talalaev are quoted from memory, and Rybnikov 2020 and Kumar were seen only by abstract in the literature check, `GD/Literature-check-missing-label.md`.

---

[^1]: `CPP/results/hw_timings.txt`; G1 §6. Convention for all footnotes: paths are relative to the repository root, with GL = `Paul/Python/gaudin_gl4`, GD = `Paul/Mathematica/Gaudin`, CPP = `Paul/Cpp/AlexCGC`, PY3 = `Paul/Python/gaudin_sl3`; G1 to G5 refer to `GL/RESULTS-G1-slice-engine-gl4.md`, `GL/RESULTS-G2-gaudin-gl4-ops.md`, `GL/RESULTS-G3-slice-ops-gl4.md`, `GL/RESULTS-G4-structure-gl4.md`, `GL/RESULTS-G5-oper-gl4.md`; § refers to a section of the file named. “The gl(3) notes” are `GD/Slice-engine-notes.tex` and its `.md` and `.pdf` twins.

[^2]: G1 §3 (column max $n_W$), G1 §6; per-irrep $n_W$ at $S=8$: G1 §4.1.

[^3]: `CPP/README.md`, sections “Highest-weight-state CGCs only” and “Benchmarking just the highest-weight state”; the algorithm is described in Ref. [1]. The routine timed in the SU(4) sweep is `compute_highest_weight_coeffs()`, without the Littlewood–Richardson step and without the descendants.

[^4]: `CPP/results/hw_timings.txt` (generated 2026-07-22); the same numbers in G1 §6 and in Table 4.

[^5]: G1 §0, G1 §4. The smallest singular-value ratio of the blocks $A_b$, which controls the conditioning, is 0.33 at $S=10$, 0.25 at $S=20$, 0.23 for $[6,4,2,0]^{\otimes2}$ and 0.18 for $[9,6,3,0]^{\otimes2}$.

[^6]: G1 §3.

[^7]: G1 §4.1–G1 §4.2.

[^8]: G1 §4.3.

[^9]: G1 §5.

[^10]: Rows 1–2: `GD/ref/gl4-cg-slices-8800.txt`, records SECTOR (fields $\mu$, $m$, $w$, $n_W$, $\dim V_\mu$, $C_2(\mu)$); rows 3–4: `GD/ref/gl4-cg-copies-6420.txt`; rows 5–7: `GD/ref/gl4-cg-copies-3210.txt`; G1 §4.1, G3 §10.

[^11]: `GD/Gaudin-sl3-CG-notes.md` §1.1, eqs. (3)–(4); the quantum characteristic polynomial of the Gaudin Lax matrix is Talalaev's [14, 15].

[^12]: G2 §1.

[^13]: G2 §2.

[^14]: G2 §4.

[^15]: G2 §6: the identities are exact in the PBW normal form (rational coefficients, exact linear algebra over libraries of invariants; no inconsistent row); the scalar remainders are verified on 81 sectors.

[^16]: G2 §5, G3 §8.2; generation modulo scalars: G2 §6.

[^17]: G2 §6e and G3 §8.1. G2 and G3 state the argument in terms of the numbers $p_\pm$ of symmetric and antisymmetric copies, $p_\pm=m_\pm$.

[^18]: G3 §1, G3 §2 (the normal-ordering engine: 656 words of gl(4), exact on $[2,1,0,0]$ and $[1,1,0,0]$ and in a second PBW order; 0 mismatches).

[^19]: G3 §3.

[^20]: G3 §3.

[^21]: G3 §3 (“Closed forms”).

[^22]: G3 §6 lists $c,c'$ for 46 sectors; the identification with $\tau_3,\tau_4$ of G5 §2.1 was checked against all 46 for these notes (0 mismatches).

[^23]: G3 §3–G3 §5.

[^24]: G3 §7, with the table of the nine products.

[^25]: G3 §8, G3 §8.1.

[^26]: G3 §9 (last block); G5 §5.3.

[^27]: G5 §2.1–G5 §2.3.

[^28]: G5 §3, G5 §3.5 (a singular-value-ratio criterion first proposed for the rejection test is not meaningful at this scale; the test used is that no solution of exact degree $d_1$ exists).

[^29]: G5 §4.1, G5 §4.2; a numerical search from random starts on 12 of the largest sectors finds no other solution, G5 §4.4.

[^30]: G5 §4.3, where the full table (also the numbers of positive-dimensional and of too-large finite solution sets) is. Subsets containing $j=1$ had time-outs on a few sectors ($\{1\}$: 5, $\{4,1\}$: 1, $\{3,1\}$: 2, $\{2,1\}$: 2 of 68); they are not counted as exact.

[^31]: G5 §6 (item 2), G5 §0 (item 5).

[^32]: G5 §5.1–G5 §5.4. $\mathrm{Cas2}=\sum_i\mu_i(\mu_i+5-2i)=C_2(\mu)$, which equals $\sum_i(\mu_i+\tfrac52-i)^2-5$ in the notebook's variables; e.g. 4 for the singlet $(1,1,1,1)$.

[^33]: This section: G4 §0; every table in `GL/structure_gl4_tables.txt`.

[^34]: G4 §3.1–G4 §3.4.

[^35]: G4 §2.2–G4 §2.5.

[^36]: G4 §4.

[^37]: G4 §5.1, G4 §5.2.

[^38]: G1 §6; the Munich column is `CPP/results/hw_timings.txt`. Full table (also $S=14$) with the breakdown of $t_a$, $t_b$ and the entries written: `GL/benchmark_gl4.md`.

[^39]: G3 §9; with PRV basis vectors in place of copies (no operators) the numbers are in G1 §7.

[^40]: G1 §4.4.

[^41]: `HANDOFF.md`, “DONE (2026-10-21 night, Ubuntu): ndsu3lib benchmark”; `Paul/Cpp/ndsu3lib-bench/RESULTS.md`.

