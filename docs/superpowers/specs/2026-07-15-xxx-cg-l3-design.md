# XXX CG L=3 notebook — design

## Context

`Paul/Mathematica/XXX/Clean/XXX_CG_L2_Clean.wb` builds the two-site rational gl(2) (XXX) spin
chain — independent reps $[\lambda_1,0],[\lambda_2,0]$, independent inhomogeneities
$\theta_1,\theta_2$, a shared companion twist $G$ — through the SoV construction and a Baxter
TQ-equation solver. This design covers porting that construction to three sites, into an already
existing but empty file, `Paul/Mathematica/XXX/Clean/XXX_CG_L3.wb`.

Unlike the L2 notebook (which was assembled by merging two validated `Experiments/` notebooks),
this is a from-scratch build directly in the target `Clean/` file — a deliberate deviation from
the project's usual Experiments-then-Clean promotion convention, per user decision.

## Scope

Full port of both L2 sections:
- **Section A** — SoV construction (GT basis through covector basis).
- **Section B** — Baxter/TQ solver (`TauEigensystem`, `QSolve`, `Q1`).

Chain setup: straight 3-site generalization — independent reps $[\lambda_1,0],[\lambda_2,0],[\lambda_3,0]$,
independent inhomogeneities $\theta_1,\theta_2,\theta_3$, one shared companion twist $G$ (same
$2\times2$ form as L2, since $G$ lives in the 2-dimensional auxiliary space regardless of chain
length).

## Architecture

Single notebook, two sections mirroring `XXX_CG_L2_Clean.wb`'s structure:

- **Section A**: GT basis/generators (reused verbatim from `XXX_su2.wb`, unchanged — L-independent),
  three-site Lax operators `L1`,`L2`,`L3`, companion twist, monodromy $T(u)=L_1(u)L_2(u)L_3(u)$,
  transfer matrix $t_{1,1}(u)$ (now cubic in $u$), fusion/Hirota hierarchy, SoV $B$-operator,
  pseudovacuum/covector basis indexed by $(s_1,s_2,s_3)$.
- **Section B**: Baxter-equation ingredients (triple-product $a(u)$, $Q_\theta(u)$), `BaxterMatrix`,
  four-coefficient extraction ($C_0,\ldots,C_3$), `dPred` as a triple convolution,
  `TauEigensystem` re-derived for the cubic case, `QSolve`, `Q1`.

## Components & key formulas

### Section A (low risk — mechanical 3-site generalization)

- `L1`, `L2`, `L3`: each site's Lax matrix shifted by its own $\theta_i$, tensored with identity
  on the other two factors.
- `T(u) = L1(u) L2(u) L3(u)`: monodromy, same ordered-product pattern as L2.
- `t11(u) = tr_a[G T(u)]`: cubic in $u$. Commutativity check at different spectral parameters
  generalizes directly (same pattern as L2 Cell 12).
- `qdetT(u)`, fusion hierarchy (`t21`, Hirota recursion): these formulas are **L-independent** in
  algebraic form (reference $T(u)$ only, not $L$ explicitly) — port verbatim; only boundary/domain
  checks need $\lambda_1,\lambda_2,\lambda_3$ ranges extended.
- SoV pseudovacuum/covector basis `x3(s1,s2,s3)`: pseudovacuum is a triple tensor product;
  covector formula generalizes the L2 pattern by adding a third fused-transfer-matrix/normalization
  factor for $\theta_3$. Validation: $B$-operator eigenvalue must be the triple product
  $(u-\theta_1-hs_1)(u-\theta_2-hs_2)(u-\theta_3-hs_3)$.

### Section B (re-derivation risk concentrated here)

- `a(u) = ∏_i (u-\theta_i-h\lambda_i)`, `Qθ(u) = ∏_i (u-\theta_i)` over $i=1,2,3$: mechanical
  extension.
- `t3Coeffs[λ1,λ2,λ3]`: extracts $C_0,C_1,C_2,C_3$ via `CoefficientList[t11(u), u]`; checks
  $C_3 = (\kappa_1+\kappa_2)\mathbb{1}$ (trivial leading term — analog of L2's $C_2$ check).
- `dPred[M,λ1,λ2,λ3]`: triple convolution count of $(m_1,m_2,m_3)$ with $\sum m_i = M$,
  $0\le m_i\le\lambda_i$ (brute sum over the small ranges; no closed form required).
- `TauEigensystem[λ1,λ2,λ3]` — **the genuinely new derivation**. $t_{1,1}(u)$ is cubic, so there
  are three genuinely operator-valued coefficients ($C_0,C_1,C_2$; only leading $C_3\propto\mathbb{1}$
  is trivial) versus L2's two ($C_0,C_1$). Jointly diagonalize a generic combination
  $C_0 + r_1 C_1 + r_2 C_2$ (two free parameters, not one), same random-parameter/full-rank/
  distinct-eigenvalue retry strategy as L2. Extract $(t_0,t_1,t_2)$ per eigenvector, then
  re-derive the $M$-extraction formula algebraically from $\tau_0(u)$'s coefficients (generalizing
  L2's eq. for $M$) and confirm via the same hard-assert pattern as L2: residuals $<10^{-15}$,
  $M$ integer within $10^{-10}$, per-sector counts matching `dPred` exactly. The extraction may
  need both $t_1,t_2$, not just one coefficient — solve for it however the algebra falls out, but
  it must still hard-assert against `dPred`.
- `QSolve`, `Q1`: same linear nullspace / Chop-before-NullSpace pattern as L2, just larger
  matrices.

**No check gets dropped**: every validation present in L2 gets an L3 counterpart (commutativity,
central charge identity, $B$-eigenvalue, `BaxterMatrix` degree-vanishing asserts, `TauEigensystem`
residuals + per-sector counts, end-to-end TQ residual).

## Parameters & testing

- **Numeric setup**: reuse L2's $\kappa_1=e^{i\zeta(3)}$ (30-digit), $\kappa_2=1/\kappa_1$, $h=1$,
  and $\theta_1,\theta_2$ as-is (residual thresholds carry over unchanged). Pick a new rational
  $\theta_3$, distinct from $\theta_1$ and $\theta_2$ modulo integers (checked explicitly), e.g.
  $\theta_3=7/3$ pending whatever's convenient once the builder sees the actual $\theta_1,\theta_2$
  values in Cell 1.
- **Representative test cases**: $(\lambda_1,\lambda_2,\lambda_3)=(1,1,1)$ (smallest nontrivial
  case, state space dimension 8) as the primary regression case for every check, plus one
  asymmetric case (e.g. $(2,1,1)$) mirroring L2's use of $(2,1)$ to catch bugs that only show up
  with unequal reps.
- **Error handling**: same as L2 — hard `Assert`/abort on failed structural checks rather than
  silent fallback; memoized routines (`ps`, `dim`, `T`, `x3`, `t3Coeffs`, `TauEigensystem`,
  `QSolve`) need the usual re-run-defining-cell-after-upstream-fix discipline.
- **Provenance comments**: since this is a from-scratch build (not a merge of validated
  Experiments notebooks like L2 was), cells should carry a short `(* why *)` comment wherever the
  L=3 generalization required a real judgment call (most importantly the $M$-extraction
  re-derivation), per the project's comment convention.

## Out of scope

- No `Experiments/` counterpart notebook — this design deliberately builds directly in `Clean/`.
- No changes to `paper/xxxCG.tex` (still states only the L=1 case) or to the L2 notebook/docs.
- No `wolfbook-docs-agent` run yet — documentation of `XXX_CG_L3.wb` is a separate follow-up once
  the notebook is validated.
