# The slice engine: why the highest-weight Clebsch–Gordan problem of $\mathrm{gl}(3)$ is tiny

*Paul Ryan (project notes, compiled from the round-6 and round-7 computations of 2026-10-21) — 21 October 2026*

*(Markdown twin of `Slice-engine-notes.tex`: same sections, same numbering of equations, tables and footnotes. Math is written in `$…$` / `$$…$$`, tables are pipe tables, footnotes are collected at the end.)*

## Abstract

For two copies of an irreducible representation $V_\Lambda$ of $\mathrm{gl}(3)$ the highest-weight (HW) Clebsch–Gordan vectors of $V_\Lambda\otimes V_\Lambda$ are normally found as a null space on a weight space of dimension $n_W$ (about 2000 for $\Lambda=[14,7,0]$).
The classical Zhelobenko (Parthasarathy–Ranga Rao–Varadarajan) criterion says that a HW vector is determined by a *slice* of dimension $w\le8$, and that the allowed slices form an explicit $m$-dimensional null space, $m$ being the multiplicity.
Used as an algorithm, together with a reduction to the slice of the Gaudin vertex operator (whose eigenbasis is the one separation of variables singles out), it gives all 344 HW vectors of the 113 irreps of $[14,7,0]\otimes[14,7,0]$ in 0.10 s (0.17 s with every full vector written out), against 206 s for the standard dense algorithm, which runs out of memory one size later.
A closed-form refinement shows that the reduced vertex operator is a Jacobi matrix on an explicit chain in the GT basis (of length at most 8 for $[14,7,0]$): the multiplicity space is a sub-chain, the slice of each copy is the three-term recurrence of that chain evaluated at the copy's eigenvalue, and the slices of all 344 copies take 6 ms.
The Q-functions supply only that eigenvalue, itself the spectrum of the Jacobi matrix: the speed-up is representation theory, not integrability.
We describe the construction, its checks and the consequences for the Gaudin/SoV programme behind it.

## 1. Punchline

Take two copies of the $\mathrm{gl}(3)$ irrep with highest weight $\Lambda=[14,7,0]$ (dimension 512) and ask for the highest-weight (HW) vectors of $V_\Lambda\otimes V_\Lambda$: for each of the 113 irreps $\mu$ in the decomposition, and for each of its $m_\mu\le8$ copies, the vector killed by the raising operators, written in the product of Gelfand–Tsetlin (GT) bases.
There are 344 of them.
The largest, the singlet, has $n_W=2080$ components, and the standard algorithm, the Munich code of Alex, Kalus, Huckleberry and von Delft [1], finds it as the null space of the raising operators on the whole 2080-dimensional weight space.
The full problem costs 206 s, and one size up the code runs out of memory.[^1]

*The problem is much smaller than that.*
A HW vector $\Psi$ is determined by its *slice* $\psi_{\mathrm{top}}$: the components in which the second tensor factor sits on its highest-weight state.
The slice is a vector in a single weight space of $V_\Lambda$, of dimension $w\le8$ (for the singlet, $w=1$).
The HW conditions fix every other component from the slice by a triangular recursion, and the only conditions on the slice itself are two $\mathrm{sl}(2)$ string conditions.
They cut the $w$-dimensional space of slices down to an $m$-dimensional subspace $M_\mu$, and that subspace *is* the multiplicity space, written down explicitly.
This is the classical criterion of Zhelobenko and of Parthasarathy, Ranga Rao and Varadarajan [2, 3]; what we add is only its use as an algorithm.[^2]
The Gaudin vertex operator, whose eigenvectors are the multiplicity copies that separation of variables (SoV) singles out, is linear in the generators of the second factor and therefore reduces to an $m\times m$ matrix on the slice, again without any reference to $n_W$.

The result: all 113 irreps and 344 copies of $[14,7,0]^{\otimes2}$, in the SoV eigenbasis, take 0.10 s, or 0.17 s with every full vector written out; the standard code needs 206 s, or 10.9 s once its linear algebra is repaired.
At $[16,8,0]$, where the standard code runs out of memory, the engine needs 0.24 s, and at $[28,14,0]$ (421 irreps, 2255 copies, $n_W$ up to 25 425) 2.0 s.[^3]
A closed-form refinement of the reduced vertex operator (Section 6) goes further: the slices of all 344 copies of $[14,7,0]^{\otimes2}$ take 5.9 ms in double precision, without any Q-function and without any object of the size of $V_\Lambda$ or of the weight space.[^4]
The timings come from different machines and are to be read as orders of magnitude (Section 5).

*Context and two points up front.*
The project behind these notes asks for the Clebsch–Gordan (CG) coefficients of $\mathrm{gl}(3)$ as explicit functions of Q-functions, obtained by SoV in the two-site rational Gaudin model, and in October 2026 it set itself a benchmark: *beat the standard Munich code using Q-function technology, for general representations, with $[14,7,0]^{\otimes2}$ as the test case*.[^5]
The engine of these notes was built as a reference and a lower bound, and it turned out to answer the benchmark by itself.
First, the speed-up is representation theory, not integrability: the question that started this was whether the Gaudin/SoV machinery makes our own pipeline faster than the standard code, and the answer was no (Section 2).[^6]
Second, the result still matters for the programme: the SoV overlap formula, which writes the overlap of rotated vacua with a Gaudin eigenstate as a bilinear form in jets of Q-functions and, as the twist goes to zero, turns it into the generating function of the CG coefficients, produces exactly the slice components,[^7] so the slice is what a Q-function formula should deliver, and the engine is both the baseline against which that formula must be timed and the reference that says what a correct answer looks like (Section 7).
A follow-up computation then showed that the slice of a copy is the three-term recurrence of an explicit chain evaluated at its eigenvalue, so that the Q-functions are needed only for that eigenvalue, itself the spectrum of a small Jacobi matrix (Sections 6 and 7).

## 2. The problem and the standard algorithm

**Punchline.** The standard algorithm solves an $n_W$-sized null space with dense linear algebra: time of order $n_W^3$ and, in the implementation we timed, memory $n_W\dim^2$.
Repairing the linear algebra buys a factor 19 and shows that the 16-fold gap between our own pipeline and the standard code was never about integrability.

### 2.1. The problem

Let $V_\Lambda$ be the irrep of $\mathrm{gl}(3)$ with highest weight $\Lambda=[\lambda_1,\lambda_2,\lambda_3]$ ($\lambda_3=0$ throughout), in the GT basis $\{e_\nu\}$ with the plain rational generators $E_{ij}$ (Appendix A), and write $V_\Lambda\otimes V_\Lambda=\bigoplus_\mu V_\mu\otimes\mathbb{C}^{m_\mu}$.
For every $\mu$ and each of its $m_\mu$ copies we want the HW vector: the vector $\Psi$ of weight $\mu$ annihilated by $E_{12}^{\mathrm{tot}}=E_{12}\otimes1+1\otimes E_{12}$ and by $E_{23}^{\mathrm{tot}}$ ($E_{13}=[E_{12},E_{23}]$ follows).
The rest of $V_\mu$ follows by lowering, so the HW vectors are the hard part and what the standard code is benchmarked on.
We call the first tensor factor site 1 and the second site 2, as in the Gaudin chain; for $m_\mu\ge2$ one also needs a distinguished basis of the multiplicity space (Section 4).

The benchmark family is $\Lambda=[2a,a,0]$, with $\dim V_\Lambda=(a+1)^3$.
At $a=7$, $\Lambda=[14,7,0]$, there are 113 irreps and 344 copies, and one irrep, $\mu=(21,14,7)$, has multiplicity 8.
The cost driver is $n_W$, the dimension of the weight-$\mu$ subspace of $V_\Lambda\otimes V_\Lambda$ in which $\Psi$ lives; it is largest for the singlet, $\binom{(a+1)^2+1}{2}$: 10, 45, 136, 325, 666, 1225, 2080, 3321 for $a=1,\dots,8$, and 25 425 at $a=14$.[^8]

### 2.2. What the standard code costs, and why

The Munich code works with GT patterns and finds the HW coefficients as the null space of the raising operators on that weight space; everything else follows by lowering.[^9]
Timing just that step on $[2a,a,0]^{\otimes2}$, the sum over all irreps is 0.1 ms, 2.8 ms, 24 ms, 0.34 s, 3.9 s, 31 s and 206 s for $a=1,\dots,7$, consistent with $n_W^3\sim a^{12}$; at $a=8$ the fifteen lowest irreps, the singlet among them, are killed by the out-of-memory killer on the 15 GB laptop.
For the singlet at $a=7$ (20.5 s in all), allocating and zeroing the system takes 1.25 s and filling it 0.013 s, and 19.1 s go into a full dense singular value decomposition (LAPACK `dgesvd`, all left singular vectors) from which a handful of null vectors are read off.
The array is stored as $n_{\mathrm{coeff}}\times\dim^2=2080\times512^2$ doubles, 4.4 GB, although the system itself is $2080\times4032$ (67 MB); that layout is what kills $a=8$, and at $a=14$ it would need 2.3 TB.[^10]

Our own pipeline (sparse null space on the weight space, then diagonalisation of the vertex operator; Mathematica) solves the same problem in 12.7 s on a different machine, 16 times faster, and the first thing we wanted to know was whether integrability was the reason.
It is not.[^11]
Keeping only the columns that are used, replacing the SVD by the lowest $m+1$ eigenvectors of $AA^{T}$ (`dsyrk` plus `dsyevr`) and linking OpenBLAS brings the singlet from 20.5 s to 0.45 s and the whole $a=7$ problem from 206 s to 10.9 s (Table 2), about as fast as our pipeline: the gap was linear algebra.
But all these codes solve the same $n_W$-sized problem, with $n_W\sim a^4$: dense methods die around $a\approx8$, sparse ones reach further, and all of them pay for the weight space.

## 3. The slice map and the theorem

**Punchline.** A HW vector is determined by its components in one weight space of $V_\Lambda$, the slice; the slice must satisfy two $\mathrm{sl}(2)$ string conditions, and those are sufficient.
The multiplicity space is the null space of two small matrices.

### 3.1. The slice and the level recursion

Write a HW vector as

$$
\Psi=\sum_q\psi_q\otimes e_q,\qquad\psi_q\in V_\Lambda, \tag{1}
$$

with $q$ running over the GT basis of site 2, and call $\psi_{\mathrm{top}}$ the *slice*: the coefficient of the highest-weight state $e_{\mathrm{top}}$ of site 2.
If $\Psi$ has weight $\mu$ then $\psi_q$ has weight $\mu-\mathrm{wt}(q)$, so $\psi_{\mathrm{top}}$ lies in the weight space $V_\Lambda[\mu-\Lambda]$, whose dimension $w$ is a weight multiplicity of $V_\Lambda$.
(For $\mathrm{gl}(2)$ this is the textbook recursion that builds the Clebsch–Gordan coefficients of $|J,J\rangle$ from $J_+|J,J\rangle=0$, with $w=1$.)

Now read the HW condition for $E_i\in\{E_{12},E_{23}\}$ at the site-2 component $r$:

$$
E_i\,\psi_r+\sum_q(E_i)_{rq}\,\psi_q=0, \tag{2}
$$

where $E_i$ acts on site 1 and $(E_i)_{rq}$ is the site-2 matrix element, the sum running over the states $q$ one step below $r$.
Equation (2) determines the components one level down from those above.
Label the weights of site 2 by $\Lambda-s\alpha_1-t\alpha_2$ and collect the $\psi_q$ of that weight in a block $X_{st}$; then (2) reads
$A_{st}X_{st}=-\big(E_{12}X_{s-1,t};\,E_{23}X_{s,t-1}\big)$, with $A_{st}$ the stacked site-2 raising blocks.
The matrix $A_{st}$ is injective for $(s,t)\neq(0,0)$, because a vector killed by both raising operators below the top weight would be a second highest-weight vector, and it depends neither on $\mu$ nor on the slice: its left inverse is computed once per $\Lambda$.
So the slice determines $\Psi$, and writing $\Psi$ out is one pass over the at most $3a^2+3a+1$ weight blocks of site 2, at a cost linear in $n_W$.[^12]

### 3.2. The theorem

Which slices are allowed?
Follow the $\alpha_1$-string of site 2: $e_{\mathrm{top}}$, $E_{21}e_{\mathrm{top}}$, …, $E_{21}^{p}e_{\mathrm{top}}$ with $p=\lambda_1-\lambda_2$.
Along it (2) makes the $k$-th component proportional to $E_{12}^k\psi_{\mathrm{top}}$, and at the end of the string, where there is no state further down, the equation reads $E_{12}\psi_{\mathrm{end}}=0$, that is $E_{12}^{p+1}\psi_{\mathrm{top}}=0$.
The $\alpha_2$-string gives $E_{23}^{\lambda_2-\lambda_3+1}\psi_{\mathrm{top}}=0$.
The classical theorem says that these two necessary conditions are also sufficient: every other equation in the overdetermined system (2) then holds automatically.

**Theorem** (Zhelobenko; Parthasarathy–Ranga Rao–Varadarajan).
The map $\Psi\mapsto\psi_{\mathrm{top}}$ is an isomorphism from the HW vectors of weight $\mu$ in $V_\Lambda\otimes V_\Lambda$ onto

$$
M_\mu=\{\,v\in V_\Lambda[\mu-\Lambda]\;:\;E_{12}^{\lambda_1-\lambda_2+1}v=0,\ \ E_{23}^{\lambda_2-\lambda_3+1}v=0\,\}. \tag{3}
$$

In particular $\dim M_\mu=m_\mu$.

The statement is classical [2, 3] (valid for any two dominant weights of any semisimple Lie algebra); what is new is only its use as an algorithm: the multiplicity space is the null space of two matrices with $w$ columns and at most $a+1$ rows, computed exactly over $\mathbb{Q}$, once per irrep.
We checked the statement on our data: $\dim M_\mu$ equals the Racah–Speiser multiplicity for every irrep in 15 representations (Section 5.1).

### 3.3. How big is the problem now?

Table 1 shows what the theorem buys at $[14,7,0]$.[^13]
The singlet has $n_W=2080$ but $w=1$: its slice is the lowest-weight vector of $V_\Lambda$, a single number, and the recursion produces the rest.
The irrep $(21,14,7)$, with $m=8$, has $w=m=8$: the string conditions are empty and $M_\mu$ is the whole weight space.
At $(20,12,10)$ the conditions bite: $M_\mu$ is 3-dimensional inside a 5-dimensional weight space.
In general $w\le a+1$ while $n_W\sim a^4$; Appendix B gives the distribution of $(w,m)$ at $[14,7,0]$.

**Table 1.** Problem sizes for five irreps of $[14,7,0]^{\otimes2}$: multiplicity $m$, slice length $w=\dim V_\Lambda[\mu-\Lambda]$ and number $n_W$ of components of the full vector.

| irrep $\mu$ | $m$ | $w$ | $n_W$ | the slice space $M_\mu$ |
| :--- | ---: | ---: | ---: | :--- |
| $(14,14,14)$ (singlet) | 1 | 1 | 2080 | the lowest-weight vector of $V_\Lambda$ |
| $(28,14,0)$ (top) | 1 | 1 | 1 | the highest-weight vector of $V_\Lambda$ |
| $(21,14,7)$ | 8 | 8 | 540 | all of $V_\Lambda[\mu-\Lambda]$ |
| $(18,12,12)$ | 1 | 3 | 1498 | a line in a 3-dimensional weight space |
| $(20,12,10)$ | 3 | 5 | 979 | a 3-dimensional subspace of a 5-dimensional one |

## 4. The multiplicity copies

**Punchline.** The vertex operator that selects the SoV copies has exactly one site-2 generator per term, so it reduces to an $m\times m$ matrix acting on slices; its eigenvalues are the roots of the Q-system's Baxter pencil, which we compute independently, in 3 s for all irreps at $a=7$.

### 4.1. The vertex operator

For $m_\mu\ge2$ (316 of the 344 copies at $[14,7,0]$ sit in such sectors) the choice of basis of the multiplicity space is real.
SoV chooses the eigenbasis of the leading-order Gaudin Hamiltonians, which for two sites in the untwisted limit are all fixed by the Casimirs of $\mu$ except one number per copy, the eigenvalue of the vertex operator

$$
V=\tfrac12\sum_{\sigma,\tau\in S_3}\mathrm{sgn}\,\sigma\;\mathrm{sgn}\,\tau\;E^{(1)}_{\sigma_1\tau_1}E^{(1)}_{\sigma_2\tau_2}E^{(2)}_{\sigma_3\tau_3}, \tag{4}
$$

36 terms, $E^{(s)}_{ij}$ being $E_{ij}$ acting on site $s$.[^14]
It commutes with the diagonal $\mathrm{gl}(3)$, so it acts on every multiplicity space; its eigenvectors are the SoV copies and its eigenvalues label them.
In terms of the Laurent data of the Gaudin Hamiltonians at site 1 (units $\theta_{12}=1$),

$$
H_1=\eta,\qquad A_1=V+\tfrac12\eta,\qquad\eta=|\Lambda|^2+c_2(\Lambda)-\tfrac12c_2(\mu),\qquad c_2(\mu)=\sum_i\mu_i(\mu_i+4-2i), \tag{5}
$$

with $|\Lambda|=\lambda_1+\lambda_2+\lambda_3$.[^15]
So $\eta$ depends only on $\Lambda$ and $\mu$, and $V$ carries all the dependence on the copy.

### 4.2. Reduction to the slice

The point is that every term of (4) contains exactly one site-2 generator.
Applying $V$ to $\Psi(v)$ and reading off the slice, for the site-2 factor $E^{(2)}_{cd}$ we need only the component of $E^{(2)}_{cd}\Psi$ at $e_{\mathrm{top}}$, and the HW conditions turn that into a site-1 operator: for $c=d$ it is $\lambda_c\psi_{\mathrm{top}}$, since $e_{\mathrm{top}}$ has weight $\Lambda$; for $c<d$ the HW condition $(E^{(1)}_{cd}+E^{(2)}_{cd})\Psi=0$ read at $e_{\mathrm{top}}$ gives $-E_{cd}\psi_{\mathrm{top}}$; for $c>d$ it vanishes, because a lowering operator never produces $e_{\mathrm{top}}$.
Hence

$$
V_{\mathrm{red}}\,v=\tfrac12\sum_{\sigma,\tau}\mathrm{sgn}\,\sigma\,\mathrm{sgn}\,\tau\;E_{\sigma_1\tau_1}E_{\sigma_2\tau_2}\,\Phi_{\sigma_3\tau_3}\,v,\quad
\Phi_{cc}=\lambda_c,\ \ \Phi_{cd}=-E_{cd}\ (c<d),\ \ \Phi_{cd}=0\ (c>d), \tag{6}
$$

an operator on the weight space $V_\Lambda[\mu-\Lambda]$ built from site-1 generators alone; 24 of the 36 terms survive.[^16]
It preserves $M_\mu$, because $V$ maps HW vectors to HW vectors, so on $M_\mu$ it is an exact rational $m\times m$ matrix in the canonical basis of $M_\mu$ (the reduced-row-echelon basis, in which the coordinates of $v$ are its free entries); its eigenvectors are the copies.
(A Hamiltonian quadratic in the site-2 generators would need components two levels below the slice; linearity in site 2 is what makes the reduction work.)

An example.
For $\mu=(21,14,7)$ of $[14,7,0]$, where $m=w=8$, the matrix of $V_{\mathrm{red}}$ in the GT basis of the weight space is tridiagonal, with

$$
\begin{array}{ll}
\text{diagonal:}&777,\ 801,\ 849,\ 921,\ 1017,\ 1137,\ 1281,\ 1449,\\[2pt]
\text{above:}&-\tfrac{189}{2},\ -120,\ -\tfrac{1815}{14},\ -128,\ -\tfrac{2535}{22},\ -\tfrac{1176}{13},\ -\tfrac{105}{2},\\[2pt]
\text{below:}&-\tfrac{441}{2},\ -120,\ -\tfrac{165}{2},\ -\tfrac{384}{7},\ -\tfrac{65}{2},\ -\tfrac{168}{11},\ -\tfrac{105}{26}.
\end{array}
$$

Its eight eigenvalues are $607.713$, $768.308$, $897.791$, $993.417$, $1064.583$, $1160.209$, $1289.692$, $1450.287$ (to three decimals), summing to the trace, 8232.[^17]

### 4.3. The same spectrum from the Q-system

The same eigenvalues come out of the Q-system, without any reference to the representation.
In the untwisted limit the Baxter equation of the two-site model is a third-order linear ODE with regular singular points at $\theta_1=1$, $\theta_2=0$ and three solutions $Q_k=l(u)\,q_k(u)$, $q_k$ a polynomial, whose coefficients depend on the state through $\mu$ and one free number, which can be taken to be $V$.
The equations for the coefficients of $q_k$ are linear in them and affine in $V$, a linear pencil, and $V$ is admissible exactly when all three systems are consistent, that is, when it is a root of the greatest common divisor of three exact consistency polynomials.
The $q_k$ have degrees $(4a-\mu_1,\ 4a-\mu_2+1,\ 4a-\mu_3+2)$, in total $6a+3$ (45 at $a=7$), against $n_W$ up to 2080.[^18]

For every irrep of $[2a,a,0]^{\otimes2}$, $a=1,\dots,14$, the number of roots equals $m_\mu$, and the roots are real and simple; on all 492 branches of $[4,2,0]$, $[6,3,0]$, $[8,4,0]$ and $[14,7,0]$ they coincide with the eigenvalues of $V_{\mathrm{red}}$ to the 30 digits printed.
Two independent computations, Baxter polynomials in Mathematica and GT matrix elements in Python, agree on 492 spectra.[^19]
The Q-system costs 3.0 s for all 113 irreps at $a=7$ and 38 s for all 421 irreps at $a=14$: 0.05–0.28 s per irrep while $n_W$ ranges over $1,\dots,25\,425$.[^20]

For multiplicity one the eigenvalue can be read off.
With $\kappa_\mu=\tfrac12(\mu_1+2)(\mu_2+1)\mu_3$,

$$
V=\kappa_\mu-\tfrac12\eta\quad(m_\mu=1),\qquad V_{\mathrm{singlet}}=\lambda_1\lambda_2(\lambda_1+1),\qquad V_{\mathrm{top}}=-\lambda_1\lambda_2, \tag{7}
$$

the last two for $\Lambda=[2a,a,0]$ and its top irrep $\mu=(2\lambda_1,2\lambda_2,0)$ ($1470$ and $-98$ at $a=7$).
The first holds to $10^{-55}$ on all 28, 40 and 56 multiplicity-one sectors at $a=7,10,14$, the other two exactly for $a=1,\dots,10$ (singlet: 6, 40, 126, 288, 550, 936, 1470, 2176, 3078, 4200); none of the three has a proof for general $\Lambda$ yet (the chain of Section 6 reproduces the first on all multiplicity-one sectors tested).[^21]

### 4.4. The engine in five steps

1. Build exact sparse GT matrices for $V_\Lambda$ ($E_{12},E_{23},E_{21},E_{32}$, and $E_{13}=[E_{12},E_{23}]$; $O(\dim V_\Lambda)$) and list the irreps and multiplicities (Racah–Speiser).
2. For each $\mu$, take the weight space $V_\Lambda[\mu-\Lambda]$ and compute $M_\mu$ as the exact null space (row reduction over $\mathbb{Q}$) of the two string operators of (3).
3. Build the exact $m\times m$ matrix of $V_{\mathrm{red}}$, eq. (6), on $M_\mu$.
4. Diagonalise it: the eigenvectors are the slices of the copies (float64, or 30 digits with mpmath).
5. Optionally, write out the full vectors by the triangular recursion of Section 3.1, using the left inverses of the $A_{st}$ computed once per $\Lambda$.

Steps 1–4 never see $n_W$; step 5 is bound by the size of the output.[^22]

## 5. What was verified, and the benchmark

**Punchline.** Every claim above was checked against an independent route, mostly in exact arithmetic; the engine is about three orders of magnitude faster than the standard code at $a=7$ and keeps going where the standard code cannot.

### 5.1. Checks

- *Counting.* $\dim M_\mu$ equals the Racah–Speiser multiplicity for every irrep in 15 representations from $[2,1,0]$ to $[14,7,0]$ (0 mismatches), and $\dim M_\mu=0$ for every dominant $\Lambda+w$ that is not an irrep, where such candidates exist. At $[14,7,0]$: 113 irreps, 344 copies, $\max m=\max w=8$, $\sum_\mu m_\mu\dim V_\mu=512^2$; at $[28,14,0]$: 421 irreps, 2255 copies, $\max m=\max w=15$.[^23]
- *Exact HW vectors from slices.* For the exact data of six small representations, $[1,1,0]$ to $[4,1,0]$ (45 branches in 37 sectors), computed earlier by an entirely different, symbolic route, the exact recursion reproduces all 45 HW vectors entry for entry; the eigenvalues of $V_{\mathrm{red}}$ satisfy (5) to $6\times10^{-50}$, the matrix identity “Laurent-datum matrix $=V_{\mathrm{red}}+\tfrac12\eta$” holds exactly in all 37 sectors, slices agree to $8\times10^{-51}$ and invariant norms to $2.9\times10^{-49}$.[^24]
- *Exact twin and negative control.* The recursion run in rational arithmetic on all PRV basis vectors: at $a=7$, 344 vectors, 155 176 entries, 104 240 exact consistency conditions, 0 HW violations; at $a=14$, 2255 vectors, 10 739 904 entries, 8 765 049 conditions, 0 violations (561 s in Python fractions). Of 78 random slices outside $M_\mu$ all 78 are rejected by the exact conditions; 102 of 102 random elements of $M_\mu$ are accepted.[^25]
- *Brute force.* A Munich-style null space written by us in numpy (dense eigendecomposition of $A^TA$ on the weight space; no slice, no recursion) for every irrep of $[4,2,0]$, $[6,3,0]$, $[8,4,0]$, $[14,7,0]$: the null dimension equals the multiplicity (0 exceptions), the first non-null eigenvalue is exactly 1, the engine's copies lie in the null space (distance $\le8.2\times10^{-13}$) and span all of it, and the eigenvalues of the full 36-term $V$ restricted to it agree with ours to $3.3\times10^{-15}$.[^26]

### 5.2. The benchmark

Table 2 puts the standard algorithm, its repaired versions, the engine and the Jacobi chain of Section 6 side by side at $a=7$, and Table 3 shows the scaling of the engine in $a$.[^27]
Between $a=7$ and $a=14$ the $n_W$-free part grows by a factor 10.1, and between 14 and 20 by 3.7 (about $a^{3.3}$ and $a^{3.7}$); it is dominated by exact rational arithmetic in Python.
The full output grows like the number of coefficients written, about $a^6$ ($1.6\times10^{5}$ at $a=7$, $1.1\times10^{7}$ at $a=14$, $1.1\times10^{8}$ at $a=20$).

**Table 2.** Highest-weight solve of $[14,7,0]^{\otimes2}$ ($a=7$): 113 irreps, 344 copies. “i5” = Ubuntu laptop (i5-1345U, 15 GB), “i7” = i7-12700KF desktop (Windows 11, 32 GB). Not like for like: rows 1–4 use different linear algebra for the same $n_W$-sized null space, row 5 adds the diagonalisation of $V$, the engine and the chain produce the SoV copies (the chain only their slices), the Q-system only their labels.

| method | singlet [s] | all irreps [s] | machine |
| :--- | ---: | ---: | :--- |
| Munich code as distributed (dense `dgesvd`, reference BLAS) | 20.5 | 206 | i5 |
| &emsp; compact columns, `dgesvd`, OpenBLAS | 7.0 | 69 | i5 |
| &emsp; lowest $m+1$ eigenvectors of $AA^T$, reference BLAS | 1.7 | 18 | i5 |
| &emsp; the same with OpenBLAS (“repaired”) | 0.45 | 10.9 | i5 |
| our sparse null space plus diagonalisation of $V$ (Mathematica) | – | 12.7 | i7 |
| slice engine, without writing the vectors | – | 0.10 | i7 |
| slice engine, all 344 full vectors written | – | 0.17 | i7 |
| Jacobi chain of Section 6, double precision, Q-free (slices only) | – | 0.006 | i7 |
| Q-system spectrum only (Baxter pencil, Mathematica) | – | 3.0 | i7 |
| Q-route end to end (jet formula, $\epsilon$-series) | – | hours (not run) | i7 |

**Table 3.** Scaling for $\Lambda=[2a,a,0]$: the Munich code (sum of the highest-weight solves over all irreps; OOM = killed by the out-of-memory killer) against the engine ($t_a$: the part that never touches $n_W$; $t_b$: everything, including all full vectors). Munich: i5, single runs for $a\ge6$; engine: i7, minimum of five cold runs, accurate to about $\pm30\%$ because other jobs shared the machine. Rows $a=9,11,13,15,17,19$ are omitted.

| $a$ | $\dim V_\Lambda$ | irreps | copies | max $n_W$ | Munich [s] | $t_a$ [s] | $t_b$ [s] |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 8 | 5 | 6 | 10 | 0.0001 | 0.001 | 0.002 |
| 2 | 27 | 13 | 19 | 45 | 0.0028 | 0.003 | 0.006 |
| 3 | 64 | 25 | 44 | 136 | 0.024 | 0.009 | 0.013 |
| 4 | 125 | 41 | 85 | 325 | 0.34 | 0.020 | 0.030 |
| 5 | 216 | 61 | 146 | 666 | 3.9 | 0.036 | 0.055 |
| 6 | 343 | 85 | 231 | 1225 | 31 | 0.061 | 0.096 |
| 7 | 512 | 113 | 344 | 2080 | 206 | 0.103 | 0.170 |
| 8 | 729 | 145 | 489 | 3321 | OOM | 0.147 | 0.243 |
| 10 | 1331 | 221 | 891 | 7381 | – | 0.300 | 0.529 |
| 12 | 2197 | 313 | 1469 | 14365 | – | 0.581 | 1.049 |
| 14 | 3375 | 421 | 2255 | 25425 | – | 1.037 | 2.024 |
| 16 | 4913 | 545 | 3281 | 41905 | – | 1.538 | 3.436 |
| 18 | 6859 | 685 | 4579 | 65341 | – | 2.473 | 5.638 |
| 20 | 9261 | 841 | 6181 | 97461 | – | 3.847 | 9.073 |

### 5.3. Caveats

*Different machines.*
The standard code ran on the i5 laptop and the engine on the i7 desktop, so the ratios (1200 against 206 s, 64 against 10.9 s) are orders of magnitude, not measurements; the only same-machine, same-language comparison is the numpy brute force above, 76–93 s for the 113 irreps of $[14,7,0]$ on a loaded machine (for scale, not a benchmark), against 0.10 s.[^28]

*Output.*
The engine returns the SoV eigenbasis (the Munich code's SVD returns an arbitrary orthonormal basis of each multiplicity space), normalised by “first non-zero slice component $=1$”; the absolute normalisation and sign from the Q side are not part of it.

*Precision.*
In float64 the largest HW residual over all irreps of $[14,7,0]$ is $2.6\times10^{-14}$ and grows by a factor of about 2.2 per unit of $a$ ($5.3\times10^{-12}$ at $a=14$, $4.8\times10^{-10}$ at $a=20$), a property of double precision in the 600-step recursion; the exact recursion rounded at the end gives $1.5\times10^{-15}$ at $a=14$.
The eigenvalues of a sector crowd (smallest relative gap between consecutive eigenvalues: 0.019 at $a=7$, 0.0049 at $a=14$, 0.0017 at $a=20$), so beyond $a\approx20$ one needs the arbitrary-precision path or the exact characteristic polynomial.

*Scope.*
$\lambda_3>0$ was only smoke-tested, there is no run beyond $a=20$, and for $[8,4,0]$ and $[14,7,0]$ there is no independent exact reference beyond the brute-force and exact-twin checks.
The reference files (`GD/ref/cg-slices-*.txt`) were re-parsed from the text and re-checked with code paths that use neither the engine's recursion nor its exact generators.[^29]

## 6. The copies are a Jacobi chain

**Punchline.** In the GT basis the reduced vertex operator is a tridiagonal (Jacobi) matrix on an explicit chain, with closed-form entries; the multiplicity space is a sub-chain cut out by the two string conditions; and the slice of a copy is the three-term recurrence of that chain evaluated at the copy's eigenvalue $V$.
The Q-functions supply only $V$, and $V$ is itself the spectrum of a Jacobi matrix of size at most $a+1$.

### 6.1. The chain and its closed forms

The weight $\mu-\Lambda$ fixes two GT entries, the bottom one $t=\mu_1-\lambda_1$ and the sum of the middle row $s=m_1+m_2=(\mu_1-\lambda_1)+(\mu_2-\lambda_2)$, so the weight space $V_\Lambda[\mu-\Lambda]$ is a chain in the single coordinate $m_1$, running from $\mathrm{lo}=\max(\lambda_2,\,s-\lambda_2,\,t,\,s-t)$ to $\mathrm{hi}=\min(\lambda_1,\,s-\lambda_3)$, of length $w=\mathrm{hi}-\mathrm{lo}+1$.[^30]
Of the 36 terms of (6), 24 survive, and in each of them exactly one factor has first index 3 and exactly one has second index 3 (either $E_{33}$ times two $\mathrm{gl}(2)$ generators, or $E_{a3}E_{3b}$ times one); since $E_{23}$ and $E_{32}$ move exactly one of $m_1,m_2$ by $\pm1$, the net shift of $(m_1,m_2)$ is $0$ or $\pm(1,-1)$, and $V_{\mathrm{red}}$ is tridiagonal for a structural reason.
Its entries are explicit.
With $m_2=s-m_1$, $x=m_1-m_2$ and $\xi=x(x+2)$, one has $V_{\mathrm{red}}\,e(m_1)=d(m_1)\,e(m_1-1)+\alpha(m_1)\,e(m_1)+u(m_1)\,e(m_1+1)$ with

$$
\begin{aligned}
u&=-\frac{(\lambda_1-m_1)(\lambda_2-m_1-1)(\lambda_3-m_1-2)(t+1-m_2)(\mu_1-\lambda_2-m_1)}{(x+1)(x+2)},\\
d&=-\frac{(\lambda_2-m_2)(t-m_1)(\lambda_1-m_2+1)(\lambda_3-m_2-1)(\mu_1-\lambda_2-m_2+1)}{x(x+1)},\\
\alpha&=A_0+\tfrac{3}{16}(s-2\lambda_3+2)\,\xi-\frac{K}{16\,\xi},
\end{aligned} \tag{8}
$$

where $K=(s-2\lambda_2)(s-2t)(s-2\lambda_1-2)(s-2\lambda_3+2)(s-2\mu_1+2\lambda_2-2)$ ($K/\xi:=0$ at $x=0$) and $A_0(\lambda;s,t)$ is an explicit cubic polynomial.
These formulas were derived symbolically from the GT matrix elements (36 terms, generic labels), not fitted.
They reproduce the engine's exact $V_{\mathrm{red}}$ entry for entry on 39 608 dominant weights of 968 representations (every $[\lambda_1,\lambda_2,\lambda_3]$ with $\lambda_1\le16$, $\lambda_3>0$ included), and the exact chains of all 192 irreps of $[4,2,0]$, $[6,3,0]$, $[8,4,0]$, $[14,7,0]$ are tabulated.[^31]
(In the example of Section 4.2, $K=0$ and the diagonal is $777+12k(k+1)$.)
Below the cut of the next subsection $u,d<0$, and $ud$ is a product of ten linear factors over $(x+1)(x+2)^2(x+3)$, each factor a positive integer (a distance to a wall of the GT range): $V_{\mathrm{red}}$ is a Jacobi matrix with simple real spectrum.

### 6.2. The multiplicity space is a sub-chain

$E_{12}$ is diagonal in the chain basis, $E_{12}\,e(t,m_1,m_2)=(m_1-t)(t+1-m_2)\,e(t+1,m_1,m_2)$, so the condition $E_{12}^{\lambda_1-\lambda_2+1}v=0$ is a coordinate cut, $m_1\le\mu_1-\lambda_2$.
The entry $u$ contains the factor $\mu_1-\lambda_2-m_1$, which vanishes exactly at the cut: the truncated chain $W_1$ is invariant under $V_{\mathrm{red}}$, the recurrence below stops by itself, and beyond the cut the products $ud$ are negative.
The $E_{23}$ condition is the same coordinate cut, $m_1\le\mu_2-\lambda_3$, in the Weyl-reflected chain $V'=UV_{\mathrm{red}}U^{-1}$, where $U$ is the longest Weyl element ($E_{ij}\to E_{4-i,4-j}$).
The two conditions are never both active, because every $E_{ii}$-eigenvalue of a weight of $V_\Lambda$ lies in $[\lambda_3,\lambda_1]$ (the image weights would need $\mu_2\ge\lambda_1+\lambda_3+1$ and $\mu_2\le\lambda_1+\lambda_3-1$ respectively).
Hence for $\mu_2\ge\lambda_1+\lambda_3$ the multiplicity space is the sub-chain $W_1$ and $m=w_1$ (its length), and for $\mu_2\le\lambda_1+\lambda_3$ it corresponds to a sub-chain $W_2$ of the second chain and $m=w_2$.
This closed-form multiplicity equals the Racah–Speiser one on all 39 608 dominant weights tested.[^32]
At $[14,7,0]$ the 113 irreps split into 71 sectors in which the whole chain is the multiplicity space ($w=m$), 21 in which the $E_{12}$ cut truncates it ($w>m$), and 21 in which the $E_{23}$ condition is active ($w>m$ again, but now the first chain is Jacobi throughout and the physical eigenvalues are those of the $m\times m$ Jacobi matrix of the second chain); for $[2a,a,0]$ the three counts are $a^2+3a+1$, $\binom a2$ and $\binom a2$ ($a=1,\dots,20$), the equality of the last two being the $E_{12}\leftrightarrow E_{23}$ duality made visible.

### 6.3. For one irrep in four the answer is symmetry alone

If $m_\mu=1$ then $M_\mu$ is a line, and the slice map and the theorem fix the HW vector up to scale: no operator needs diagonalising and no Q-function needs solving.
At $[14,7,0]$ this is the case for 28 of the 113 irreps; the other 316 of the 344 copies live in sectors with $m\ge2$ (at $[28,14,0]$: 56 irreps, and 2199 of 2255 copies).[^33]
For $m=1$ even the eigenvalue is explicit, eq. (7), and the chain reproduces it ($\lambda_3=0$, $\lambda_1\le15$): $V=\kappa_\mu-\tfrac12\eta$ on all 2510 multiplicity-one sectors, and the mean of $V$ over the copies equals $\kappa_\mu-\tfrac12\eta$ in all 6255 sectors, being the trace of the explicit matrix.[^34]
Everything that distinguishes a Q-function approach from the PRV theorem therefore sits in the sectors with $m\ge2$, and there only in the choice of basis of $M_\mu$.

### 6.4. The slice is a three-term recurrence

For every eigenvalue $V$ of the sector, the slice of the corresponding copy is

$$
p_0=1,\qquad p_1=\frac{V-\alpha_0}{d_0},\qquad p_{k+1}=\frac{(V-\alpha_k)\,p_k-u_{k-1}\,p_{k-1}}{d_k}\quad(k\ge1), \tag{9}
$$

on the chain below the cut and zero beyond it, where $\alpha_k=\alpha(\mathrm{lo}+k)$, $u_k=u(\mathrm{lo}+k)$ and $d_k=d(\mathrm{lo}+k+1)$ are the diagonal, sub-diagonal and super-diagonal entries of (8).
It is unique up to scale because the spectrum is simple, and it holds in every sector, with $V$ taken from either chain or from the Q-system.
The component at chain position $k$ (counted from 0) is thus a polynomial of degree $k$ in $V$, equivalently in the accessory parameter $A_1=V+\tfrac12\eta$ of the copy's oper, and $V=\alpha_0+d_0\,p_1$: the accessory parameter is a linear function of the first slice ratio.

*Checks.*
Against the slice engine at 50 digits, on all 3 638 copies of $[4,2,0]$, $[6,3,0]$, $[8,4,0]$, $[14,7,0]$, $[20,10,0]$ and $[28,14,0]$, the maximal relative deviation lies between $4\times10^{-50}$ and $4\times10^{-42}$ (the growth with $a$ is the engine's eigensolver: $7\times10^{-78}$ and $1\times10^{-73}$ at 80 digits for $[14,7,0]$ and $[28,14,0]$); a further 848 copies of 12 other representations, nine of them with $\lambda_3>0$, agree with the engine to $6\times10^{-48}$.
Without the engine, the recurrence slice satisfies both string conditions and $V_{\mathrm{red}}v=Vv$ exactly, up to the working precision, for 6 796 copies of 285 representations ($\lambda_1\le10$, $\lambda_3\ge0$; worst residual $10^{-39}$, 0 failures), and the 587 non-physical eigenvalues of $V_{\mathrm{red}}$ on the chain all violate the $E_{23}$ condition.[^35]

*Cost.*
The chain code builds no representation matrix and no object of the size of $V_\Lambda$ or of the weight space ($O(a^2)$ set-up per $\Lambda$, $O(w)$ per sector).
In double precision all 344 copies of $[14,7,0]^{\otimes2}$ take 5.9 ms (closed-form chains 3.0 ms, eigenvalues 2.5 ms, forward recurrences 0.33 ms) and all 2255 copies of $[28,14,0]^{\otimes2}$ take 39 ms, which is 13–15 and 17–23 times faster than the $n_W$-free part of the slice engine in the same session (0.090 s and 0.89 s); given $V$, one copy costs 14 $\mu$s at 30 digits and 1 $\mu$s in double precision.
The plain forward recurrence loses accuracy in double precision at large $a$ ($3\times10^{-8}$ at $a=14$, $5\times10^{-5}$ at $a=20$); the twisted variant (forward from the bottom, backward from the top) stays below $8\times10^{-14}$ up to $a=30$, at 6.8 ms and 51 ms for the two cases above.[^36]
The full HW vectors still need the triangular recursion of Section 3.1 and are bound by the size of the output.

### 6.5. What the Q-system adds, and what it does not

The Q-system and the chain share the spectrum, exactly.
The consistency polynomial of the Baxter pencil (the gcd of the three echelon systems, with the spectral parameter $b$ mapped by $V=\kappa_\mu-\tfrac12b-\tfrac12\eta$) equals the exact characteristic polynomial of $V_{\mathrm{red}}$ in all 192 sectors of $[4,2,0]$, $[6,3,0]$, $[8,4,0]$, $[14,7,0]$.
The continuant of the first chain below the cut equals it whenever the multiplicity space is the whole truncated chain (the $71+21=92$ sectors with $w_1=m$ at $[14,7,0]$) and is divisible by it in the 21 $E_{23}$-active sectors, where the continuant of the second chain's sub-chain equals it instead.
The three single Baxter systems also explain an earlier surprise: in the $E_{23}$-active sectors the second system alone already has degree $m$, in the $E_{12}$-cut sectors the third alone, so the gcd of the first two systems, the one our first solver used, over-counts roots exactly in the 21 $E_{12}$-cut sectors, which together with the $a+1$ sectors with $\mu_1=2\lambda_1$ (where a system has no unknowns) are the $\binom a2+a+1=29$ of 113 irreps of the footnote to Section 4.3.
But the chain recurrence is not the Baxter recurrence.
The chains have length at most $a+1$ where the Q-polynomials have degree up to $4a+2$; the off-diagonal products are a product of ten linear factors divided by $(x+1)(x+2)^2(x+3)$ against a polynomial sextic in the index; the diagonals differ.
A regression of the slices on the coefficients, on the Taylor coefficients at both sites and on the values of the Q-polynomials finds no copy-independent map to the interior slice components; only the last component is proportional to a ratio of two Q-data, in 39 of the 85 multi-copy sectors at $[14,7,0]$, and these ratios look like normalisation data rather than direction data.[^37]
The Q-function information that fixes the direction of a slice is $V$ alone.
In one sentence: the eigenvector is representation theory (the chain), the Q-functions supply $V$, and $V$ is itself the spectrum of a Jacobi matrix of size at most 15 at $[28,14,0]$, so the Q-system step is replaceable, at the price of using no Q-function at all.

## 7. What this means for the integrability programme

**Punchline.** The benchmark was beaten by representation theory, and the Jacobi chain of Section 6 sharpens this: the slice of a copy is a three-term recurrence evaluated at its eigenvalue $V$, the Q-functions supply only $V$, and $V$ is itself the spectrum of a Jacobi matrix of size at most $a+1$, so even the Q-system step is replaceable.
The Q-function formula is not a faster way to compute Clebsch–Gordan coefficients, and we do not present it as one; what it offers is a structural statement and a precise point of contact with the slice.

### 7.1. The balance sheet

*Munich's cost was never intrinsic* (an $n_W$-sized null space in dense linear algebra, 10.9 s once repaired), *and the Q-route as formulated cannot compete with the engine on wall-clock.*
To extract the slice of one copy from the SoV formula one needs the Laurent data of the Q-system, which is cheap (Section 4.3); their $\epsilon$-series to a depth equal to $N_B=4\lambda_1-2\lambda_2-k_B$ with $k_B=2\lambda_1-(\mu_1-\mu_3)$ (at $[14,7,0]$: order 21 for the top state, 28 for the eight copies of $(21,14,7)$, 41 for $(28,13,1)$), at 20–85 s per copy for the series alone on the deeper of the 14 sample branches we ran; the kernel, which at $\lambda_1=14$ is 8464 polynomial entries with 1.69 million monomials in all and coefficients up to $10^{51}$ (2 minutes and 3 GB to build exactly); and a contraction of these that cancels by $10^{26}$ to $10^{29}$, so that it has to be done with some 60 digits.[^38]
Extrapolated to all 344 copies this is hours on one core (15–20 minutes on 20 cores), against 10.9 s for the repaired Munich code and 0.17 s for the engine; the family was not run (`PY/RESULTS-bethe-series-scale.md` §6, last paragraph before “Surprises”).
The Q-route would only “win” where the standard code cannot run at all, and there the engine is faster still, by three to four orders of magnitude even if the Q-route is spread over 20 cores.

*The Q-system spectrum can compete, and does:* 3.0 s for all irreps at $a=7$ (against 206 s and 10.9 s), 38 s at $a=14$, independently of $n_W$; but it delivers the labels of the copies, not their vectors, and the same labels are the spectra of the explicit Jacobi matrices of Section 6, a few milliseconds for all irreps at $a=7$ in double precision and without any Q-function.

*What the formula is for.*
The relation that the programme established reads

$$
\sum_{i,j}G_{ij}(\epsilon\hat\omega_1,\epsilon\hat\omega_1';\theta_{12})\;\nu_i(\theta_1;\epsilon)\,\nu_j(\theta_2;\epsilon)
=\epsilon^{N_B}\;\hat c(\phi,\hat\omega)\;\frac{\langle\varphi_u(\hat\omega)\otimes e_1|\Psi_0\rangle}{X_B}+O(\epsilon^{N_B+1}). \tag{10}
$$

The left side is built from the Q-functions of the twisted two-site Gaudin eigenstate, with twist $\zeta=\epsilon\phi$: the $\nu_i(\theta_s)$ are Plücker-type ratios of their jets at the two sites, $G$ is a universal kernel, and $\hat\omega_1,\hat\omega_1'$ are two rotation parameters; $N_B$ is the order at which the left side first becomes non-zero, all lower orders cancelling.
The right side is the generating function of the slice, $\langle\varphi_u(\hat\omega)\otimes e_1|\Psi_0\rangle=\sum_\nu[e_1^{T}\rho_\Lambda(n')]_\nu\,C_\nu$ with $C_\nu=\langle e_\nu\otimes e_1|\Psi_0\rangle$ the components of $\psi_{\mathrm{top}}$, times an explicit prefactor $\hat c$ and one constant $X_B$ per branch.
Relation (10) holds with no further free constant on all 26 branches of $[2,1,0]$, $[3,1,0]$, $[3,2,0]$ that were checked; it is verified, not proved.[^39]
Its value is structural: the CG coefficients of $\mathrm{gl}(3)$ are Q-numbers of the Gaudin model, depending on the state only through a handful of Laurent data, and the multiplicity basis it selects is the eigenbasis of $V$, whose spectrum the Q-system reproduces exactly.

*Where the two routes meet.*
Both objects are $w$-dimensional.
In (10) the unknown is the slice, that is the $m$ coordinates of $\psi_{\mathrm{top}}$ in the explicit basis of $M_\mu$, and $m+1$ evaluations of the Q-side expression at different $\hat\omega$ determine them: the “dictionary”.
On the 26 known branches it recovers every slice to 53–116 digits and confirms independently that the constant is 1; at $[14,7,0]$ it has been tested only on synthetic data generated from known slices (random elements of $M_\mu$ in all 113 irreps, and the engine's own 344 copies), not yet on Q-side data.[^40]
The engine is thus both the baseline against which the Q-route is timed and the reference against which any Q-route output is compared; and it sharpens what the Q-route can add: for $m=1$ nothing, since the slice is fixed by symmetry (Section 6.3), for $m\ge2$ the choice of the eigenbasis of $V$, whose only Q-side datum is $V$ itself (Section 6.5).

### 7.2. Next

1. Explain, rather than check, why the consistency polynomial of the Baxter pencil equals the characteristic polynomial of the chain (Section 6.5; verified in all 192 sectors up to $[14,7,0]$), and ask whether the interior slice components have any Q-side expression at all: it would have to involve Wronskian-type minors of the jets, as in eq. (10), and has not been explored.
2. Assemble the jet-formula route end to end at $[4,2,0]$, $[6,3,0]$, $[8,4,0]$ and on a sample of $[14,7,0]$ copies, validated against the engine, with its cost stated as in Section 7.1.
3. One benchmark table with all columns (Munich as distributed and repaired, sparse null space, engine, Jacobi chain, Q-spectrum, Q-route) from $a=7$ to $a=14$.

How to restate the project's goal in the light of this is open.[^41]

## A. Conventions

- *One-site basis.* GT patterns of $\mathrm{gl}(3)$ with rows $(\lambda_1,\lambda_2,\lambda_3)$, $(r_1,r_2)$, $(k)$, in the order of `gl3_reps.GL3Rep(Lambda).patterns`: ascending lexicographic in $(r_1,r_2,k)$, highest-weight state *last*. (The Wolfram notebooks use the opposite order, descending in $(k,r_1,r_2)$, highest-weight state first; the map is by pattern.)
- *Generators and metric.* Plain rational, non-unitary GT matrices acting on kets as columns; $E_{13}=[E_{12},E_{23}]$, $E_{31}=[E_{32},E_{21}]$; weights are the eigenvalues of $(E_{11},E_{22},E_{33})$, $\alpha_1=(1,-1,0)$, $\alpha_2=(0,1,-1)$. The invariant metric is the diagonal matrix $G$ with $E_{ij}^TG=GE_{ji}$, $G[\mathrm{HW}]=1$; the floating-point recursion runs in the orthonormal basis $\tilde e_q=e_q/\sqrt{g_q}$ for conditioning ($g$ reaches $1.7\times10^{37}$ at $[14,7,0]$, $10^{103}$ at $[28,14,0]$) and the output is converted back.
- *Two sites, slices, copies.* $\Psi=\sum\Psi[\nu,q]\,e_\nu\otimes e_q$, site 1 the first factor (slow index), flat index $(\nu-1)\dim V_\Lambda+q$, 1-based; the slice is $\Psi[\,\cdot\,,q_{\mathrm{top}}]$. The basis of $M_\mu$ is the reduced row echelon form of the null space (one vector per free column, so the coordinates of $v\in M_\mu$ are its free entries); copies are normalised so that the first non-zero slice component is 1 and ordered by increasing eigenvalue of $V_{\mathrm{red}}$. Multiplicities by Racah–Speiser (Klimyk), $\rho=(2,1,0)$.

Sources: `PY/RESULTS-slice-engine.md` §0; `GD/ref/cg-slices-1470.txt` (header).

## B. The $(w,m)$ table for $[14,7,0]$

Table 4 counts the 113 irreps of $[14,7,0]^{\otimes2}$ by slice length $w$ and multiplicity $m$.[^42]
The diagonal $w=m$ holds 71 sectors, where the slice space is the whole weight space; the remaining 42 have $w>m$.
The 28 irreps with $m=1$ are those of Section 6.3.
In general $\Lambda=[2a,a,0]$ has $2a(a+1)+1$ irreps, $4(a+1-m)$ of multiplicity $m\le a$ and one with $m=a+1$, and $\tfrac23a(a+1)(a+2)+(a+1)$ copies (checked at $a=1,2,3,4,7,14$).

**Table 4.** Irreps of $[14,7,0]^{\otimes2}$ by slice length $w$ (rows) and multiplicity $m$ (columns); blank means none.

| $w\backslash m$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | irreps | copies |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | ---: | ---: |
| 1 | 16 |  |  |  |  |  |  |  | 16 | 16 |
| 2 | 4 | 14 |  |  |  |  |  |  | 18 | 32 |
| 3 | 4 | 4 | 12 |  |  |  |  |  | 20 | 48 |
| 4 | 4 | 4 | 4 | 10 |  |  |  |  | 22 | 64 |
| 5 |  | 2 | 4 | 4 | 8 |  |  |  | 18 | 72 |
| 6 |  |  |  | 2 | 4 | 6 |  |  | 12 | 64 |
| 7 |  |  |  |  |  | 2 | 4 |  | 6 | 40 |
| 8 |  |  |  |  |  |  |  | 1 | 1 | 8 |
| total | 28 | 24 | 20 | 16 | 12 | 8 | 4 | 1 | 113 | 344 |

## C. Reproduction

The engine (Python 3.13 with numpy, scipy, sympy, mpmath), from `Paul/Python/gaudin_sl3`:

```
python -X utf8 slice_engine.py all       # every check and table quoted above
python -X utf8 slice_engine.py bench 20  # Table 3 (a = 1..20)
# further stages: w4, brute, tridiag, closed, exact 14, ref all
```

In Python, `import slice_engine as se; se.reference_irrep((14,7,0),(21,14,7))` returns the exact basis of $M_\mu$, the exact $V_{\mathrm{red}}$, the copies and the full vectors (`PY/RESULTS-slice-engine.md` §11); reference data: `GD/ref/cg-slices-*.txt`.

The Q-system (Mathematica): `Get` `GD/Experiments/gl3-GT-tools.wl`, `Minimal_complete_working-fast.wl`, `GD/Experiments/Gaudin-Qsystem-scale.wl`, then `runAll[7]`. The Munich variants: `CPP/cgc_hw_prof.cpp` with `HW_MODE=orig|compact|eig` (`CPP/README.md`).

## References

[1] A. Alex, M. Kalus, A. Huckleberry and J. von Delft, *A numerical algorithm for the explicit calculation of SU(N) and SL(N,C) Clebsch–Gordan coefficients*, arXiv:1009.0437.

[2] D. P. Zhelobenko, *Compact Lie Groups and Their Representations*, AMS Transl. Math. Monographs 40, 1973.

[3] K. R. Parthasarathy, R. Ranga Rao and V. S. Varadarajan, *Representations of complex semi-simple Lie groups and Lie algebras*, Ann. of Math. 85 (1967) 383–429.

---

[^1]: `CPP/results/sweep_2aa0.csv`, `CPP/results/full_decomp_2aa0_a7.txt`, `CPP/results/full_decomp_2aa0_a8.txt`; summarised in `HANDOFF.md`, section “Done (2026-10-21): SU(3) [2a,a,0] timing sweep”. Times are for the highest-weight solve only, on an Ubuntu laptop (i5-1345U, 15 GB). Convention for all footnotes: paths are relative to the repository root, with PY = `Paul/Python/gaudin_sl3`, GD = `Paul/Mathematica/Gaudin`, CPP = `Paul/Cpp/AlexCGC`; § refers to a section of the file named.

[^2]: Statement and use: `PY/SPEC-slice-engine.md` §2; `GD/Gaudin-sl3-round6-plan.md` §2.1.

[^3]: `PY/RESULTS-slice-engine.md` §1 (item 5), §5. The 0.10 s never touches $n_W$ (exact generators, Racah–Speiser, null spaces, reduced vertex operator, eigendecomposition); the 0.17 s adds the metric, the left inverses of Section 3 and all 344 reconstructions (155 176 coefficients). Repaired Munich code: `HANDOFF.md`, “Done (2026-10-21, Ubuntu): experiments 2–3”.

[^4]: `PY/RESULTS-jacobi-cg.md` §1 (items 2, 5), §7 (complexity table).

[^5]: Paul Ryan, 2026-10-21, quoted in `GD/Gaudin-sl3-round6-plan.md`, header, and `HANDOFF.md`, “ROUND 6” block, “Original dispatch note”; the notes: `GD/Gaudin-sl3-CG-notes.md` §0, §5, §9.

[^6]: `HANDOFF.md`, “NEWEST (2026-10-21): why is our sl(3) CG pipeline $\sim$16$\times$ faster than the Munich C++ — and is integrability the reason?” and “Done (2026-10-21, Ubuntu): experiments 2–3”.

[^7]: `GD/Gaudin-sl3-round6-plan.md` §1, the paragraph after the displayed formula; `PY/RESULTS-gauss-dictionary.md` §1 (item 4).

[^8]: `PY/RESULTS-slice-engine.md` §2 (last paragraph, “counting laws”); `GD/Experiments/Gaudin-Qsystem-scale.md` §6.2.

[^9]: `CPP/README.md`, sections “Highest-weight-state CGCs only” and “Benchmarking just the highest-weight state”; the algorithm is described in Ref. [1].

[^10]: `HANDOFF.md`, “Done (2026-10-21): SU(3) [2a,a,0] timing sweep” (sum at $a=8$: $>284$ s, 15 irreps missing) and “Done (2026-10-21, Ubuntu): experiments 2–3” (breakdown); the 2.3 TB: `PY/RESULTS-slice-engine.md` §5; per-irrep data: `CPP/results/full_decomp_2aa0_a7.txt`, `CPP/results/full_decomp_2aa0_a8.txt`.

[^11]: `HANDOFF.md`, “NEWEST (2026-10-21): why is our sl(3) CG pipeline …” (“What our approach is”, “Current answer”); repaired variants: “Done (2026-10-21, Ubuntu): experiments 2–3”, data in `CPP/results/full_decomp_2aa0_a7_cgc_hw_prof_*.txt`, code `CPP/cgc_hw_prof.cpp`.

[^12]: `PY/SPEC-slice-engine.md` §2 (“Level recursion”); `PY/RESULTS-slice-engine.md` §3. The smallest singular-value ratio of the blocks $A_{st}$, which controls the conditioning of the recursion, is 0.26 at $a=7$ and 0.19 at $a=14$.

[^13]: `PY/slice_engine_output.txt`, block “named irreps of [14,7,0] x [14,7,0]”; `GD/ref/cg-slices-1470.txt`, records SECTOR (fields $\mu$, $m$, $w$, $n_W$, …).

[^14]: `PY/SPEC-slice-engine.md` §1 (“The vertex operator”); `Minimal_complete_working-fast.wl` (`vTriples`); `GD/Experiments/gl3-GT-tools.wl`.

[^15]: `PY/RESULTS-slice-engine.md` §4 (relation found and verified on 45 exact branches); `GD/Experiments/Gaudin-Qsystem-scale.md` §0 (item 4) and §3.1 ($H_1=\eta$ on every sector). $A_1$ is called $a_1$ elsewhere in the project; here $a$ is the family parameter.

[^16]: `PY/SPEC-slice-engine.md` §2 (“$V_{\mathrm{red}}$”); `PY/RESULTS-slice-engine.md` §1 (item 3) and §3 (the full 36-term $V$ on the reconstructed vectors agrees with (6) to $1.4\times10^{-14}$ on all entries).

[^17]: `GD/ref/cg-slices-1470.txt`, records VRED, CHARPOLY and COPY for $\mu=(21,14,7)$; the same eigenvalues appear in `GD/Experiments/Gaudin-Qsystem-scale.md` §10, row $(21,14,7)$.

[^18]: `GD/Experiments/Gaudin-Qsystem-scale.md` §3 (“The corrected solver”), §5 (degrees); the pencil is also written as a three-term recursion in `PY/RESULTS-bethe-series-scale.md` §1. Using only two of the three systems returns too many roots on $\binom a2+a+1$ irreps of $[2a,a,0]^{\otimes2}$ (29 of 113 at $a=7$); see §2 of the first file.

[^19]: `GD/Experiments/Gaudin-Qsystem-scale.md` §0 (items 2–3), §4 (table). A separate implementation of the pencil agrees on the same 492 branches to $4.6\times10^{-27}$: `PY/RESULTS-bethe-series-scale.md` §1 (table).

[^20]: `GD/Experiments/Gaudin-Qsystem-scale.md` §0 (item 6), §6.1 (CPU times; wall times 3.4 s and 48 s).

[^21]: `PY/RESULTS-slice-engine.md` §8; `PY/RESULTS-jacobi-cg.md` §1 (item 7); `PY/slice_engine_output.txt`, block “closed forms”; `GD/Experiments/Gaudin-Qsystem-scale.md` §0 (item 4), where the mean of $A_1$ over the copies of a sector is also found to equal $\kappa_\mu$ (numerical identity). $\kappa_\mu$ is called $\alpha$ in the project files.

[^22]: `PY/RESULTS-slice-engine.md` §5 (breakdown of $t_a$ and $t_b$), §11 (usage). Pure Python with exact rationals; the exact $V_{\mathrm{red}}$ is 58–70% of the $n_W$-free time.

[^23]: `PY/RESULTS-slice-engine.md` §2 (table; for $[2a,a,0]$ with $a\le7$ every dominant $\Lambda+w$ is an irrep, so the zero-multiplicity half of the test needs extra representations such as $[3,0,0]$, $[4,1,0]$, $[5,5,0]$, $[9,1,0]$); `PY/RESULTS-gauss-dictionary.md` §1 (item 3).

[^24]: `PY/RESULTS-slice-engine.md` §4 (table); data `GD/ref/untwisted-exact-110.txt` and the files for 210, 220, 310, 320, 410, produced by `GD/Experiments/Gaudin-untwisted-exact.wb`.

[^25]: `PY/RESULTS-slice-engine.md` §3 (table, last paragraph) and §7 (last item).

[^26]: `PY/RESULTS-slice-engine.md` §7 (“Brute-force null space”).

[^27]: Sources: `HANDOFF.md` (the two “Done (2026-10-21)” sections and the round-6 blocks); `PY/RESULTS-slice-engine.md` §5 (all rows $a=1,\dots,20$); `GD/Experiments/Gaudin-Qsystem-scale.md` §6.1; `PY/RESULTS-bethe-series-scale.md` §5–§6 for the Q-route.

[^28]: `PY/RESULTS-slice-engine.md` §5 (“Caveat on the noise” and the context paragraph), §7.

[^29]: `PY/RESULTS-slice-engine.md` §3 (“HW residual vs size”), §5 (last paragraph), §6–§7, §10; `HANDOFF.md`, “NEWEST (2026-10-21): why …” (the SVD returns an arbitrary basis).

[^30]: `PY/RESULTS-jacobi-cg.md` §0 (conventions), §2 (chain coordinate).

[^31]: `PY/RESULTS-jacobi-cg.md` §1 (item 1), §2 (bandwidth table), §3 (closed forms; $A_0$ is the cubic printed there); exact chains: `PY/jacobi_cg_tables.txt`.

[^32]: `PY/RESULTS-jacobi-cg.md` §1 (item 3), §4, §7 (T2–T4).

[^33]: `PY/RESULTS-gauss-dictionary.md` §1 (item 3), §3 (item 1).

[^34]: `PY/RESULTS-jacobi-cg.md` §1 (item 7).

[^35]: `PY/RESULTS-jacobi-cg.md` §1 (item 2), §5 (deviation table, 80-digit rerun), §7 (T6).

[^36]: `PY/RESULTS-jacobi-cg.md` §1 (item 5), §5 (timing and float64 accuracy), §7 (complexity table).

[^37]: `PY/RESULTS-jacobi-cg.md` §1 (items 4, 6, 8), §3 (characteristic polynomial), §6 (J5).

[^38]: `PY/RESULTS-bethe-series-scale.md` §0, §5, §6 (depth law $N_B$ measured on 188 branches; times per branch); `PY/RESULTS-kernel-L14.md` §1–§2 (kernel size, build time and memory; cancellation in the contraction); `HANDOFF.md`, ROUND-6 block, items K and S.

[^39]: `GD/Gaudin-sl3-round6-plan.md` §1 (the displayed formula and the paragraph after it); `GD/Gaudin-sl3-CG-notes.md` §5, §7, §9; `HANDOFF.md`, round-5 blocks N1, W5, S1.

[^40]: `PY/RESULTS-gauss-dictionary.md` §1 (items 4, 5, 9), §3 (item 1), §4.

[^41]: `HANDOFF.md`, ROUND-6 block, “ASSESSMENT GIVEN TO PAUL” (recommended round-7 priorities, open decisions); `GD/Gaudin-sl3-round6-plan.md` §4, §7.

[^42]: Tabulated from the SECTOR records of `GD/ref/cg-slices-1470.txt`; the pattern is quoted in `PY/RESULTS-gauss-dictionary.md` §1 (item 3).
