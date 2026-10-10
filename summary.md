> **SUPERSEDED (2026-09-16).** The L=2 determinant described below was built and verified, and the untwisting limit was done. Current state: `HANDOFF.md` (ACTIVE section), `Paul/Mathematica/XXX/Experiments/SU3_L2_Untwisting.tex`, and the `.md` files of `SU3_L2_Baxter_Twisted_Final.wb` and `SU3_L2_Untwisting.wb`. Note: the "per-site/per-ladder rows" suggestion below turned out to be WRONG. The correct rows are `(α,i,j)`, with α labelling bare vs `e^{2πiu/h}` measure, each a full contour.

# Handoff: SU(3) twisted Baxter/SoV normalization, L=1 and L=2

## User's actual goal

Build clean Wolfbook notebooks for the adjoint SU(3) XXX chain with reference states

\[
\langle0_1|=\langle1,0,\ldots,0|
e^{a_{12}E_{12}}e^{a_{13}E_{13}}e^{a_{23}E_{23}},
\qquad
|0_1\rangle=e^{b_{21}E_{21}}e^{b_{31}E_{31}}e^{b_{32}E_{32}}|1,0,\ldots,0\rangle,
\]

normalize every right eigenstate by \(\langle0|\Psi\rangle=1\), every left eigenstate by
\(\langle\Psi|0\rangle=1\), and match the resulting Hilbert-space scalar products to a normalized
functional SoV determinant made from direct and dual Q-functions.

For L=2, the reference state is the tensor product

\[
\langle0_2|=\langle0_1|\otimes\langle0_1|,
\qquad
|0_2\rangle=|0_1\rangle\otimes|0_1\rangle.
\]

The user explicitly permits numerical diagonalization for L=2 and feeding the transfer-matrix
eigenvalues into Baxter.

Most important latest instruction: **do not sum over a basis of wave functions to obtain the
normalization. Compute the Q-determinant itself and normalize it directly.**

## Mandatory workflow

- This is Wolfram/Wolfbook work. Read `/home/paul/.agents/skills/wolfbook/SKILL.md` completely.
- Use Wolfbook cell tools and the live VS Code kernel for `.wb` files. Do not edit notebook JSON.
- Read notebook context before editing, run relevant cells visibly, inspect messages/results, and
  explicitly save.
- Display user-facing results as `Grid`, not bare `Association`.

## Completed L=1 notebook

Final file:

`Paul/Mathematica/XXX/Experiments/SU3_L1_Baxter_Twisted_Final.wb`

This notebook was cold-run successfully after a fresh kernel. Cells 2--7 all passed with no kernel
messages. It was explicitly saved. It contains only the final construction:

- the 8-dimensional adjoint representation `[2,1,0]`;
- twisted fused transfer matrices;
- all eight simultaneous left/right eigenstates;
- direct and dual Baxter Q-functions for every state;
- the exactly normalized 3x3 functional determinant;
- both branches relating `aij` and `bij`;
- exact equality of the complete 8x8 physical and Q-determinant scalar-product matrices.

The final cold run reported `True` for both complete 8x8 matches.

### Exact L=1 determinant normalization

The determinant rows are

```wl
rowTriple = {{1, 1}, {1, 2}, {2, 1}};
```

and

\[
\mathcal S_Q(A,B)=C_{\rm vac}\,
\frac{\det N^{AB}}{c_A^{\rm ref}\,\Psi_B^{\rm ref}},
\qquad
\mathcal S_Q({\rm vac},{\rm vac})=1.
\]

The right reference wave function is

\[
\Psi_B^{\rm ref}=
Q_{B,1}(\theta+h)
\left[
Q_{B,1}(\theta)Q_{B,2}(\theta-h)
-Q_{B,2}(\theta)Q_{B,1}(\theta-h)
\right].
\]

The exact coefficient of that right reference wave function, with no on-shell deletion of either
term, is

\[
\begin{aligned}
c_A^{\rm ref}=\frac{1}{4h^6}\Big[&
Q^1_A(\theta+h)
\big(Q^1_A(\theta)Q^2_A(\theta+h)-Q^2_A(\theta)Q^1_A(\theta+h)\big)\\
&-Q^1_A(\theta)
\big(Q^1_A(\theta)Q^2_A(\theta+2h)-Q^2_A(\theta)Q^1_A(\theta+2h)\big)
\Big].
\end{aligned}
\]

This is exactly the user's previously requested off-shell coefficient of
`Psi_B[{1},{0,0}]`. Do not simplify it by imposing an on-shell vanishing relation.

The vacuum calibration found in the notebook is

\[
C_{\rm vac}=\frac{z_2z_3^3}
{(z_1-z_2)(z_1-z_3)^2(z_2-z_3)}.
\]

### L=1 relations between alphas and betas

The common relations are

\[
b_{21}=-\frac{z_1(z_1-z_3)}{a_{12}z_2(z_2-z_3)},
\qquad
b_{32}=-\frac{z_2^2(z_1-z_3)}
{a_{23}(z_1-z_2)z_3^2}.
\]

Branch I:

\[
b_{31}=-\frac{z_1^2}{a_{13}z_3^2},
\qquad
\frac{a_{12}a_{23}}{a_{13}}
=-\frac{z_2(z_1-z_3)}{z_1(z_2-z_3)}.
\]

Branch II:

\[
b_{31}=-\frac{z_1}{a_{13}z_3},
\qquad
\frac{a_{12}a_{23}}{a_{13}}
=-\frac{z_1-z_3}{z_2-z_3}.
\]

Both branches give exact full 8x8 equality.

## Existing L=2 reference material

Read these before changing conventions:

- `Paul/Mathematica/XXX/Experiments/SU3_Adjoint_Twisted_SoV_L2.wb`
- `Paul/Mathematica/XXX/Experiments/SU3_Adjoint_Twisted_SoV_L2.md`
- `Paul/Mathematica/XXX/Experiments/SU3_Adjoint_Twisted_SoV_L2.tex`
- `Paul/Mathematica/XXX/Experiments/SU3_Adjoint_SoV_ScalarProduct_L2.wb`
- `Paul/Mathematica/XXX/Experiments/SU3_Adjoint_SoV_ScalarProduct_L2.md`

Useful established L=2 facts:

- Hilbert-space dimension is 64.
- Site representation is `[2,1,0]` at both sites.
- The tested generic chain point is
  `thetaValues = {1/3, -8/57}`, `hValue = 1`.
- Current twisted test point is `zTwist = {5/2, 7/3, 1/6}`.
- This twist is **not unimodular**: `Times @@ zTwist = 35/36`.
- Therefore the last term in both reduced Baxter equations must contain the twisted quantum
  determinant `(Times @@ zTwist) qdetL2`, not bare `qdetL2`.
- Direct polynomial degrees for a state of weight `mu` are `4-mu_i`.
- Dual polynomial degrees are `mu_j`.
- With the corrected twisted quantum determinant, all 64 direct and dual Baxter solves gave a
  maximum relative residual `0.*10^-79` at 90-digit working precision.

The established full-contour six columns over `M=1/(nu1 nu2 nu3)` are

\[
(-uQ^jQ_i,-Q^jQ_i,uQ^jQ_i^{[-2]},Q^jQ_i^{[-2]},
(u+h)Q^{j[2]}Q_i,Q^{j[2]}Q_i).
\]

## Current L=2 notebook and what is trustworthy

Current work-in-progress file:

`Paul/Mathematica/XXX/Experiments/SU3_L2_Baxter_Twisted_Final.wb`

Cells 2--5 are useful and were successfully evaluated in a clean kernel:

1. self-contained adjoint generators and two-site twisted transfer matrices;
2. numerical simultaneous diagonalization of all 64 states, weight sector by weight sector;
3. corrected direct/dual Baxter solves for all 64 states;
4. the product right reference-point wave function.

Cell 7 is **not correct yet** and must be replaced. It currently tries a full-contour 6x6
determinant with six distinct Q-pair rows and then extracts a reference monomial coefficient. That
coefficient vanishes identically for every possible distinct row sextuple. The cell therefore emits
division-by-zero/indeterminate messages. Do not treat its displayed result as valid.

Also edit cell 1 and the final markdown cell after the correct construction is known; they currently
describe the unfinished full-contour attempt.

## Failed approach that must not be repeated

An earlier draft built the full 64x64 product wave-function matrix `W`, formed

```wl
coefficientMatrix2 = determinantMatrix2 . Inverse[W];
```

and used its reference column. The user explicitly rejected this: **do not sum over wave functions**.
That code has been removed from the notebook and must not be reintroduced.

A second naive attempt divided by the product left reference wave function. It failed because the
left product reference value vanishes for some of the 64 eigenstates. This is the L=2 analogue of
the L=1 warning that the correct divisor is the coefficient coupled to the right minimum
configuration, not an arbitrarily guessed left wave-function value.

## Key structural discovery for the correct L=2 determinant

The failed full-contour determinant used six distinct Q-pair rows, initially

```wl
{{1,1},{1,2},{1,3},{2,1},{2,2},{3,1}}
```

and then all 84 subsets of six distinct pairs from `Tuples[Range[3],2]` were screened.

For the product right reference wave function, the leading monomial contains four `Q_1` factors and
two `Q_2` factors: direct Q-index content `{1,1,1,1,2,2}`. A determinant with six **distinct** global
Q-pair rows cannot have four rows with direct index 1, because only `(1,1),(1,2),(1,3)` exist. Direct
mixed-coefficient extraction confirmed that the desired reference monomial coefficient is zero for
all 84 distinct row sextuples.

This is strong evidence that the product-reference L=2 determinant must retain the **site/layer
label**, allowing the L=1 row triple to repeat once per site:

```wl
siteRows = Flatten[
  Table[{site, pair}, {site, 2},
    {pair, {{1,1},{1,2},{2,1}}}], 1];
```

Equivalently, the six rows should be `(site, i, j)`, not six distinct global `(i,j)` pairs. This
gives exactly the required direct and dual Q-index content `{1,1,1,1,2,2}`.

The next agent should therefore derive/use **per-ladder residue functionals**, one for each site,
rather than first summing the two ladders into a single full-contour row. This is also consistent
with the older untwisted L=2 notebook, whose rows are labeled `(s,i,j)` and whose simplest working
triple on one ladder is `(s,1,1),(s,1,2),(s,2,1)`.

Likely structure to test first:

- rows: for each site `s=1,2`, the repeated L=1 triple
  `{{1,1},{1,2},{2,1}}`;
- each row is the residue sum restricted to that site's pole ladder;
- six columns are the two `u` moments (`beta=0,1`) of the three L=1 column structures, or the
  equivalent six-column basis obtained from the L=2 transfer-matrix difference coefficients;
- then extract the coefficient of the four monomials in the product right reference wave function
  directly from this determinant by multilinearity, with no wave-function-basis inversion;
- divide by that coefficient and by `psiRightReference2`, fix one `Cvac`, and compare the complete
  64x64 matrix against `physicalScalar2[q]` for both alpha/beta branches.

The exact choice and order of the six columns must be derived from the per-ladder Baxter identity,
not guessed from the already-summed full-contour determinant. The existing L2 `.tex` and the
`rowR[A,B,s,i,j]` construction in `SU3_Adjoint_SoV_ScalarProduct_L2.wb` are the best starting points.

## Direct coefficient extraction without wave-function sums

Once the correct site-labeled determinant is built, the coefficient of a degree-six right-Q
monomial can be extracted efficiently by multilinearity.

For each of the six target right-Q factors, form its 6x6 coefficient matrix `M_k` by differentiating
the determinant matrix entries with respect to that Q sample (or by direct Kronecker-delta selection
inside the residue sum). Then

\[
[x_1\cdots x_6]\det\Big(\sum_{k=1}^6x_kM_k\Big)
=\sum_{\pi\in S_6}
\det\big(M_{\pi(1)}[:,1],\ldots,M_{\pi(6)}[:,6]\big).
\]

This is determinant-only normalization and obeys the user's instruction. For the product reference
wave function there are four monomials with signs `+,-,-,+`. Their signed determinant coefficients
should agree if exactly one reference configuration couples. Use the first coefficient for all 64
states and check the other three at least on the vacuum (preferably on several states) as a direct
consistency test.

## Physical L=2 normalization already implemented

The current notebook constructs `rightMatrix2` from simultaneous eigenvectors and
`leftMatrix2 = Inverse[rightMatrix2)`. Check/fix the closing bracket if copying this text; the actual
notebook code uses the correct syntax `Inverse[rightMatrix2]`.

For branch `q=1,2` it defines

```wl
bra02[q_] := KroneckerProduct[bra01[q], bra01[q]];
ket02[q_] := KroneckerProduct[ket01[q], ket01[q]];
rightReferencePhysical2[q_] := bra02[q] . rightMatrix2;
leftReferencePhysical2[q_] := leftMatrix2 . ket02[q];
physicalScalar2[q_] := IdentityMatrix[64]/
  Outer[Times, leftReferencePhysical2[q], rightReferencePhysical2[q]];
```

This correctly represents the fully normalized physical scalar-product matrix, because
`leftMatrix2.rightMatrix2 = IdentityMatrix[64]` numerically.

The notebook uses exact rational sample alpha values satisfying each branch only to perform the
numerical check:

```wl
a13 -> 3/5;
a23 -> 4/7;
a12 -> alphaConstraint * (3/5)/(4/7);
```

The displayed final relations must remain symbolic as written in the L=1 section above.

## Required final verification

Before handing back:

1. Replace the incorrect L2 determinant cell with the site-labeled/per-ladder construction.
2. Cold-restart the Wolfbook kernel.
3. Run the L2 notebook top to bottom with message policy `stop`.
4. Require no kernel messages.
5. Report:
   - maximum simultaneous-eigenvector residual;
   - maximum Baxter residual;
   - maximum relative off-diagonal Q-determinant entry;
   - agreement of the four signed reference-monomial coefficients;
   - nonvanishing of all 64 `c_A^ref`;
   - vacuum entry exactly/numerically 1;
   - maximum absolute residual of the complete 64x64 match for branch I and branch II.
6. Explicitly save both final notebooks.

Do not claim L=2 is complete until those two full 64x64 residuals are small at high precision.
