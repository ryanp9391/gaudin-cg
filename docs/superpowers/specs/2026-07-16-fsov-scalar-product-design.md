---
# FSoV scalar product on numerical Q-functions — design

## Context

`Paul/Mathematica/XXX/Experiments/FSoV_L2.wb` builds the Functional Separation of Variables
(FSoV) measures $\mu_i(u)$ for the two-site XXX chain — one measure per site $i=1,2$, each
supported on that site's inhomogeneity tower $\theta_i+hk$ — together with a bracket pairing
$\langle f\rangle_i$ (residue sum of $f\cdot\mu_i$ over site $i$'s poles) and a $2\times2$
determinant `det[λ1,λ2][f]` built from $\langle f\rangle_i,\langle uf\rangle_i$. For
$f=Q_A(u)Q_B(u)$, this determinant is the SoV scalar product $\langle Q_A,Q_B\rangle$: it must
vanish whenever $A\neq B$ label distinct eigenstates, and expanding it turns the product into a
sum over SoV basis states. `FSoV_L2.wb` only tests this machinery against a generic
$f(u)z^{u/h}$ placeholder, not against real Baxter Q-functions.

`Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb` is a copy of
`Clean/XXX_CG_L2_Clean.wb` (Section A: SoV construction; Section B: Baxter/TQ solver ending in
`Q1[λ1,λ2][M,n][u] = κ1^(u/h) q_n(u)`, the actual numerical Baxter Q-function for weight-$M$
state $n$). This design covers porting the FSoV machinery into a new Section C of that notebook
and using it to numerically verify orthogonality/non-degeneracy of the `Q1` Q-functions already
computed there.

## Scope

- Port `mu`, `Nn`/`nn`, `br`, `det` from `FSoV_L2.wb` into a new Section C of
  `XXX_CG_L2_V2.wb`, verbatim (they close over the same global `θ1,θ2,h` Section A/B already
  set in cell 0; `λ1,λ2` remain local pattern variables, no adaptation needed).
- Build the full Gram matrix of `det[λ1,λ2][Q1(...) Q1(...)]` over **all** states of one
  representation, covering multiple $M$ sectors (including cross-$M$ pairs) in a single matrix.
- Verify numerically: off-diagonal entries vanish to (report) high precision; diagonal entries
  are nonzero. No comparison to a predicted closed-form norm — existence/vanishing only.

Out of scope: symbolic/general-$f$ proof of orthogonality (that's what `FSoV_L2.wb` already
explored); testing multiple different $(\lambda_1,\lambda_2)$ representations; comparing norms
to a predicted formula.

## Test case

$(\lambda_1,\lambda_2) = (2,1)$ — already partially exercised in `XXX_CG_L2_V2.wb`'s existing
smoke tests (`Psi[2,1][2,2]`, `QSolve[2,1][2,0]`). Dimension $d=(\lambda_1+1)(\lambda_2+1)=6$,
splitting into weight sectors $M=0,1,2,3$ with `dPred[M,2,1] = 1,2,2,1` — two genuinely
degenerate $M$-sectors ($M=1,2$, each multiplicity 2) in one representation, which covers
same-$M$ degeneracy and (via the full Gram matrix) cross-$M$ pairs at the same time.

## Architecture

New cells appended after Section B (cell 26+) in `XXX_CG_L2_V2.wb`:

1. **Section C header comment** — states purpose, references `FSoV_L2.wb` as the source.
2. **`mu[λ1,λ2][i][u]`** — ported verbatim from `FSoV_L2.wb` cell 1.
3. **`Nn`/`nn`** — ported verbatim from `FSoV_L2.wb` cell 2.
4. **`br[λ1,λ2][i][f]`** — ported verbatim from `FSoV_L2.wb` cell 3.
5. **`det[λ1,λ2][f]`** — ported verbatim from `FSoV_L2.wb` cell 4.
6. **Port sanity-check cell** — re-evaluate `det[1,1][f[u] z^(u/h)]` (the same test as
   `FSoV_L2.wb`'s cell 5) and confirm it reproduces the same closed form, catching any
   transcription error before building on top of it.
7. **`GramMatrix[λ1,λ2]`** — collects `TauEigensystem[λ1,λ2]`'s keys `{M,n}` (sorted by
   `{M,n}`), and builds the $d\times d$ table
   $G_{(M,n),(M',n')} = \det[\lambda_1,\lambda_2]\big[Q_1[\lambda_1,\lambda_2][M,n][u]\cdot
   Q_1[\lambda_1,\lambda_2][M',n'][u]\big]$.
8. **Verification cell** — computes `G = GramMatrix[2,1]`, splits into diagonal vs
   off-diagonal entries, and (matching the notebook's existing `Message`+`Abort` hard-assert
   style used by `BaxterMatrix`/`TauEigensystem`/`QSolve`) asserts
   `Max[Abs[off-diagonal]]` is at the numeric noise floor of the notebook's existing 30-digit
   working precision, reporting the number of digits it vanishes to (same reporting style as
   the L3 Baxter-residual check: "exact zero at N tracked digits"). Separately reports all 6
   diagonal values, confirming they are well above that floor (nonzero norms) — no comparison
   to a predicted formula.

## Risks

- `det[λ1,λ2][f]` divides by `f` evaluated at $\theta_1,\theta_2$ — a genuine risk only if some
  `Q1` state happens to vanish at $\theta_1$ or $\theta_2$, which is not expected generically for
  the notebook's fixed numeric $\theta_1=1/3,\theta_2=1/7,\kappa_1=N[\exp(i\zeta(3)),30]$, but the
  verification cell should surface a clear error (not a silent wrong answer) if it happens.
- `Residue` inside `br` needs simple, distinct poles — true generically for these fixed
  numeric $\theta_1,\theta_2,h$, but not re-derived from scratch; relies on the same reasoning
  already implicit in `FSoV_L2.wb`.
