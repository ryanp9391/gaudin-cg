# Conversation transcript: closed-form rectangular CG overlaps

This transcript begins at the requested sentence. It contains the visible user and assistant messages only; internal tool calls and hidden instructions are omitted.

---

## Assistant

The exact small cases expose a much simpler scaling law than the matrix implementation: every allowed output is

\[
\mu=(s+p-a,\ s+p-b,\ b,\ a),\qquad 0\le a\le b\le \min(s,p),
\]

and its complete overlap sector has only \(a+1\) entries. Their projective values are

\[
\langle r,a-r,b-r,r|\Psi_\mu\rangle/\langle0,a,b,0|\Psi_\mu\rangle
=(-1)^r/\binom br.
\]

I’ve checked this exactly against the existing GT-matrix reconstruction through `s=3,p=2`. The corresponding exponential overlap collapses to one Plücker power, so extraction is linear in `min(s,p)`, not in representation dimension.

---

## User

you crashed my pc, whats the conclusion?

---

## Assistant

Sorry — the live Wolfbook editing bridge accumulated stuck processes while I was updating the notebook. I stopped them, but I should have backed off sooner. The workbook was not successfully updated.

The mathematical conclusion is clean:

For

\[
[s,s,0,0]\otimes[p,p,0,0],
\]

every output irrep is uniquely parameterized by

\[
\mu=(s+p-a,\ s+p-b,\ b,\ a),
\qquad 0\le a\le b\le\min(s,p).
\]

Only these first-site bras can overlap:

\[
(n_{13},n_{14},n_{23},n_{24})
=(r,a-r,b-r,r),\qquad r=0,\ldots,a.
\]

After fixing

\[
\langle0,a,b,0|\Psi_\mu\rangle=1,
\]

all overlaps are

\[
\boxed{
\langle r,a-r,b-r,r|\Psi_\mu\rangle
=\frac{(-1)^r}{\binom br}
}.
\]

Equivalently, the complete exponential generating function is

\[
\boxed{
\frac{u_{23}^{\,b-a}}{a!\,b!}
\left(u_{14}u_{23}-u_{13}u_{24}\right)^a
}.
\]

All other \(\langle n_{13},n_{14},n_{23},n_{24},0,0,0,0|\Psi_\mu\rangle\) vanish by Cartan weight.

I verified this exactly against the existing explicit GT/tensor reconstruction through \(s=3,p=2\), including:

- \(a=b=1:\ \{1,-1\}\)
- \(a=1,b=2:\ \{1,-\tfrac12\}\)
- \(a=b=2:\ \{1,-\tfrac12,1\}\)

Thus the overlap extraction itself is representation-free and costs only \(O(a)\). The Q-functions and bilinear-residue determinant can still be computed alongside it as Q-side verification, without constructing `dim[rep]`, generator matrices, Kronecker products, or `Psi`.

What remains is only safely placing this API into the workbook and running the requested clean-kernel \(s=p=5\) timing. I did not complete that benchmark, so I cannot honestly claim the measured sub-second result yet.

---

## User

Wait so youre saying there is an exact closed form expression for the CG coefficients?

---

## Assistant

Yes—but specifically for the restricted first-site overlaps you asked for, not automatically for every CG coefficient in a conventional GT tensor basis.

For an output irrep

\[
\mu=(s+p-a,\ s+p-b,\ b,\ a),
\qquad 0\le a\le b\le\min(s,p),
\]

the only Cartan-compatible overlaps are

\[
\langle r,a-r,b-r,r,0,0,0,0|\Psi_\mu\rangle,
\qquad r=0,\ldots,a.
\]

Up to one overall normalization, they have the exact closed form

\[
\frac{\langle r,a-r,b-r,r|\Psi_\mu\rangle}
{\langle0,a,b,0|\Psi_\mu\rangle}
=
\frac{(-1)^r}{\binom br}.
\]

So their generating function is simply

\[
\frac{u_{23}^{b-a}}{a!\,b!}
(u_{14}u_{23}-u_{13}u_{24})^a
\]

when \(\langle0,a,b,0|\Psi_\mu\rangle=1\).

This completely determines the HWS projectively because the output irrep occurs with multiplicity one. However:

- it does not fix the absolute normalization;
- it is not yet a closed formula for every component of the HWS in Matteo’s/GT tensor basis;
- the Q-determinant should supply/check the remaining normalization, but that identification still needs to be established carefully.

So the surprising result is: the entire family of overlaps you requested appears to have a universal exact closed form, while the bilinear-residue machinery is mainly needed for the Q-normalization and comparison with Matteo’s convention.

---

## User

Is this formula known anywhere in the literature?

---

## Assistant

I’ll check whether this precise binomial overlap formula is published explicitly, and distinguish that from the well-known ingredients behind it—rectangular Littlewood–Richardson decomposition, Plücker relations, and standard monomial theory.

---

## Assistant

The underlying formula is classical, but I have not found this exact four-index overlap formula stated verbatim as an SU(4) Clebsch–Gordan formula.

It is a direct consequence of the Plücker/standard-monomial realization:

- \([s,s,0,0]\) is the degree-\(s\) component of the homogeneous coordinate ring of \(\mathrm{Gr}(2,4)\).
- \(\mathrm{Gr}(2,4)\) has one quadratic Plücker relation.
- Products of rectangular Schur representations are multiplicity-free, so each permitted highest-weight vector is unique up to scale. This rectangular multiplicity-free result is classical; see Stembridge’s [Multiplicity-Free Products of Schur Functions](https://deepblue.lib.umich.edu/bitstream/handle/2027.42/41839/26-5-2-113_10050113.pdf%3Bsequence%3D1).
- Standard monomial theory constructs these representations and their bases using Plücker coordinates; see Seshadri’s [Introduction to the Theory of Standard Monomials](https://link.springer.com/book/10.1007/978-981-10-1813-8).

In our variables, the relevant standard monomial is simply

\[
u_{23}^{\,b-a}
\det\begin{pmatrix}
u_{14}&u_{13}\\
u_{24}&u_{23}
\end{pmatrix}^{a}
=
u_{23}^{\,b-a}(u_{14}u_{23}-u_{13}u_{24})^a.
\]

Expanding it and accounting for the factorials in the exponential basis gives

\[
\frac{\langle r,a-r,b-r,r|\Psi_\mu\rangle}
{\langle0,a,b,0|\Psi_\mu\rangle}
=
(-1)^r\frac{r!(b-r)!}{b!}
=
\frac{(-1)^r}{\binom br}.
\]

So I would describe it as:

> a simple specialization of the classical Plücker/standard-monomial highest-weight construction, expressed in our unusual raising-operator overlap convention.

General SU(\(N\)) CG literature usually presents algorithms or GT-basis coefficients rather than this restricted coherent-bra slice—for example the [Alex–Kalus–Huckleberry–von Delft algorithm](https://arxiv.org/abs/1009.0437). I did not locate a paper writing precisely

\[
\langle r,a-r,b-r,r|\Psi\rangle=(-1)^r/\binom br
\]

with our normalization.

Thus the closed form itself is probably not a fundamentally new representation-theory theorem. What may be new and relevant to this project is showing that Matteo’s bilinear Q-residue determinant selects and normalizes exactly this standard-monomial state directly from the Gaudin Q-functions.

---

## User

Now i want a new wb which takes these CG coefficients, builds the state with irrep mu and imposes that it is highest weight, then show me that the overlap is precisely that you wrote. so do direct diagonalization. maybe do sparse matrices etc. goal is to verify your formula explicitly

---

## Assistant

I’ll make this as a separate, standalone verification workbook. It will construct sparse GT generators, restrict immediately to the requested Cartan sector, form the total raising-operator constraint, obtain the HWS both from its exact nullspace and from diagonalizing \(A^\dagger A\), then normalize by the first requested overlap and compare every overlap against \((-1)^r/\binom br\). No Q-function assumptions will enter this check.

---

## Assistant

The independent construction is now working exactly. In a fresh live kernel, the sparse-sector nullspace and the lowest eigenvector of \(A^\dagger A\) agree projectively, and for the nontrivial \(s=p=2,\mu=(2,2,2,2)\) channel the direct overlaps are exactly \(\{1,-\tfrac12,1\}\), with zero raising residual. I also swept every irrep for all \(1\le s\le3,\ 1\le p\le2\), plus all ten channels at \(s=p=3\); every check passed exactly.

---

## Assistant

Created [CGClosedFormDirectVerification.wb](C:/Users/pryan/git-projects/gaudin-cg/Paul/Mathematica/Gaudin/Experiments/CGClosedFormDirectVerification.wb) with companion implementation [CGClosedFormDirectVerification.wl](C:/Users/pryan/git-projects/gaudin-cg/Paul/Mathematica/Gaudin/Experiments/CGClosedFormDirectVerification.wl).

It independently:

- builds sparse GT generators;
- restricts to the requested Cartan sector;
- imposes all three total raising constraints;
- finds the exact one-dimensional nullspace;
- independently diagonalizes \(A^\dagger A\);
- constructs the sparse tensor-product HWS;
- evaluates the monomial bras directly;
- compares against
  \[
  \frac{(-1)^r}{\binom br}.
  \]

For \(s=p=2,\ \mu=(2,2,2,2)\), it obtains

\[
\{1,-\tfrac12,1\}
\]

exactly. The numerical diagonalized vector has projective error zero, while the smallest Gram eigenvalue is \(2.17\times10^{-15}\).

Exact sweeps passed for:

- every channel with \(1\le s\le3,\ 1\le p\le2\);
- all ten channels at \(s=p=3\), in about \(0.31\) seconds in the live kernel.

The construction uses no Q-functions or determinant identities, so this is an independent explicit verification of the formula.
