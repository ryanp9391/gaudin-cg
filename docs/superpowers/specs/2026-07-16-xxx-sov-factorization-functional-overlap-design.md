---
# L=2 SoV-basis factorization and the functional CG overlap — design

## Context

`Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb` already builds, for the two-site XXX chain
(physical space $V_{\lambda_1}\otimes V_{\lambda_2}$, independent inhomogeneities $\theta_1,\theta_2$,
shared companion twist $G$):

- **Section A** — single-site building blocks (`ps`/`dim`/`J0`/`Jp`/`Jm`/`JJ`, `Ee[λ][i,j]`,
  single-site Lax `L[λ][i,j][u]`), the two-site monodromy `T`, the fundamental and fused transfer
  matrices `t2[λ1,λ2][a,s][u]` (via `qdetT` + Hirota), the SoV $B$-operator `B2 = T[1,1]`, the
  single-site lowest-weight covectors `xSingle[λ][{0}]`, and the **L=2 SoV covector basis**
  `x2[λ1,λ2][{s1,s2}] = (xSingle[λ1][{0}]⊗xSingle[λ2][{0}]) . t2[1,s1][θ1] . t2[1,s2][θ2] / norm`,
  already verified to diagonalize `B2` with eigenvalue $(u-\theta_1-h s_1)(u-\theta_2-h s_2)$ for
  $\lambda_1,\lambda_2\in\{1,2,3\}$. The `x2` normalization divides by the **two-site**
  $Q_\theta(u)=(u-\theta_1)(u-\theta_2)$.
- **Section B** — the Baxter/TQ solver: `TauEigensystem` → right eigenvectors `Psi[λ1,λ2][M,n]`
  (last component 1) and eigenvalues `τ`, `QSolve`, and the numerical Baxter Q-function
  `Q1[λ1,λ2][M,n][u] = κ1^(u/h) q_n(u)`.
- **Section C** — the FSoV measures/brackets/determinant `mu`/`Nn`/`nn`/`br`/`det` (ported from
  `FSoV_L2.wb`), the Gram matrix `GramMatrix`/`G21` over states `states21`, verified numerically
  diagonal for $(\lambda_1,\lambda_2)=(2,1)$.
- **Section D** — left eigenvectors `PsiLeft`, biorthogonality, `braket`/`normRatio`.

`paper/xxxCG.tex` currently writes up only the **single-site** ($L=1$) SoV construction: the
covectors $\langle\svx^{[\lambda]}_s| = \langle\svx^{[\lambda]}_0|\,t_{1,s}(\theta)/\prod_{p=1}^s
Q_\theta(\theta-\lambda h+(p-1)h)$, verified symbolically to be left $B$-eigenvectors for
$\lambda=1..5$. There is **no** standalone length-1 SoV covector, length-1 measure/bracket, or any
CG-overlap object in the notebook yet — everything is built directly at $L=2$.

This design covers two related additions, appended as **Section E** (operatorial) and **Section F**
(functional) to `XXX_CG_L2_V2.wb`. Sections A–D are read-only; they are not modified.

## Overview of the two claims

**Part 1 (operatorial, Section E).** The transfer matrices $t_{1,s_1}(\theta_1)$ used to build the
SoV basis act only on the *first* site of the chain. Concretely, the $s_2=0$ slice of the L=2 SoV
basis factorizes as a tensor product of a genuine **length-1** SoV covector on site 1 with the
site-2 pseudovacuum:
$$
\langle x_{s_1,0}| \;=\; c\,\big(\langle x_{s_1}^{[\lambda_1]}| \otimes \langle x_0^{[\lambda_2]}|\big),
$$
where $\langle x_{s_1}^{[\lambda_1]}|$ is built the paper's L=1 way (single-site fused transfer
matrix on the site-1 lowest-weight covector) and $c$ is a nonzero scalar arising from the
$Q_\theta$-vs-$Q_{\theta_1}$ normalization mismatch.

**Part 2 (functional, Section F).** Set up the length-1 FSoV measure/bracket, normalized so that
$\mathrm{br}_1^{[\lambda_1]}[f] = f(\theta_1) + (\text{higher tower})$, and form the 1-site
"$1\times1$ determinant" $\mathrm{br}_1^{[\lambda_1]}\big[k^{u/h}\,Q\big]$ for the already-computed
two-site Q-function $Q = Q_1[\lambda_1,\lambda_2][M,n]$, with $k$ a free parameter. The claim is
that this functional object equals an **operatorial CG overlap**: a global unipotent rotation
$g(k)$ acting on the first site of the two-site chain, applied to the factorized pseudovacuum and
contracted with the two-site eigenvector $Psi$:
$$
\mathrm{br}_1^{[\lambda_1]}\big[k^{u/h}\,Q_1[\lambda_1,\lambda_2][M,n]\big]
\;\overset{?}{=}\;
\big(\langle x_0^{[\lambda_1]}|\,g(k) \;\otimes\; \langle x_0^{[\lambda_2]}|\big)\,\cdot\,
|Psi[\lambda_1,\lambda_2][M,n]\rangle ,
$$
with $g(k)=\exp[\varphi(k)\,E]$, $E\in\{E_{12},E_{21}\}$ a single-site generator. The exact
dictionary — which generator $E$, the map $\varphi(k)$, and whether the two sides agree exactly or
up to a $k$-independent scalar — is **not** assumed; it is discovered by numerical matching and
reported. This is the Clebsch–Gordan overlap the project is ultimately after.

## Decisions (settled during brainstorming)

- **Where:** two new sections (E, F) appended to `XXX_CG_L2_V2.wb`. Sections A–D unchanged.
- **Genuine standalone L=1 block** (not L=2 with a trivialized site 2): build a real single-site
  transfer-matrix hierarchy and covector, so $\langle x_{s_1}^{[\lambda_1]}|$ is a bona fide
  length-1 SoV object matching the paper.
- **Part 1 verification is up-to-scalar with the scalar reported** and cross-checked against the
  predicted normalization factor. Part 1 is a pass/fail check (a known identity).
- **Part 2 discovers the dictionary numerically** (generator, $\varphi(k)$, scalar) starting from
  $\lambda_1=1$ and validating universality; exact-vs-scalar is decided from the data, not asserted.
  Part 2 *reports* the dictionary rather than asserting the correspondence itself.
- **One spec (this document); one implementation plan with two phases** and a hard checkpoint
  between them — Part 1's factorization must verify before Part 2 (which uses the factorized
  $\langle x_0|g\otimes\langle x_0|$) is attempted.

## Part 1 — Operatorial factorization (Section E)

### Build

1. **Single-site fused transfer matrices** `t1[λ][a,s][u]`, the direct L=1 analogue of Section A's
   `t2` block:
   - `t1[λ][1,1][u] := Sum[L[λ][i,j][u-θ1] G[[j,i]], {i,2},{j,2}]` (trace of $G$·single-site Lax);
   - `qdet1[λ][u] := L[λ][1,1][u-θ1].L[λ][2,2][u-θ1-h] - L[λ][2,1][u-θ1].L[λ][1,2][u-θ1-h]`
     (single-site quantum determinant; note both Lax factors are the *same* site's L-matrix at
     $u$ and $u-h$, unlike the two-site `qdetT`);
   - `t1[λ][2,1][u] := χ2 qdet1[λ][u]`; boundary `t1[λ][0,s]=t1[λ][a,0]=Id`;
   - Hirota `t1[λ][1,s>1][u] := t1[λ][1,s-1][u].t1[λ][1,1][u+(s-1)h] - t1[λ][2,1][u+(s-1)h].t1[λ][1,s-2][u]`.
2. **Length-1 SoV covector** with the paper's genuine 1-site normalization ($Q_{\theta_1}(u)=u-\theta_1$):
   `xL1[λ][{s}] := xSingle[λ][{0}] . t1[λ][1,s][θ1] / Product[Qθ1[θ1 - λ h + (p-1) h], {p,1,s}]`,
   with `xL1[λ][{0}] := xSingle[λ][{0}]`.
3. **Internal check (independent of the factorization):** `B1[λ][u] := L[λ][1,1][u-θ1]`, then
   assert `xL1[λ][{s}] . B1[λ][u] == (u-θ1-s h) xL1[λ][{s}]` (symbolically in `u`, `Chop` residual)
   for $\lambda=1..3$, $s=0..\lambda$. This validates the L=1 covectors against the paper before
   they are used.

### Verify the factorization

For $\lambda_1,\lambda_2\in\{1,2,3\}$, $s_1\in 0..\lambda_1$, $s_2=0$:

- `tp = KroneckerProduct[xL1[λ1][{s1}], xSingle[λ2][{0}]]`, `lhs = x2[λ1,λ2][{s1,0}]`;
- over the nonzero components of `tp`, compute the ratio list `lhs[[j]]/tp[[j]]` and assert it is
  constant to the notebook's noise floor (max deviation from the mean $< 10^{-15}$, matching the
  existing hard-assert convention; the working precision is 30 digits so genuine agreement sits at
  $\sim 10^{-20}$ or better);
- with the ratio defined as `c = lhs[[j]]/tp[[j]]`, record that scalar `c[λ1,λ2,s1]` and hard-assert
  it equals the predicted
  `1 / Product[(θ1 - λ1 h + (k-1) h) - θ2, {k,1,s1}]` (the $Q_\theta$-vs-$Q_{\theta_1}$ mismatch:
  `x2` divides by $\prod (u-\theta_1)(u-\theta_2)$ at the site-1 tower points $u=\theta_1-\lambda_1
  h+(k-1)h$, `xL1` only by $\prod(u-\theta_1)$, so `c = ∏ Q_{θ1}/∏ Q_θ = ∏ 1/(u-θ_2)` at those
  points — the **reciprocal** of the product of $(\cdots-\theta_2)$ factors). Note this predicted
  `c` is **independent of $\lambda_2$** (site 2 only contributes its unit pseudovacuum covector);
  computing `c` across $\lambda_2\in\{1,2,3\}$ and confirming it is $\lambda_2$-independent is a free
  extra check.

Two hard asserts fire on failure (constancy of the ratio; match to the predicted scalar), in the
notebook's existing `Message`+`Abort` style. Final report: the table of $c[\lambda_1,\lambda_2,s_1]$
alongside the predicted values, and the max ratio-constancy deviation (digits vanished).

## Part 2 — Functional CG overlap (Section F)

### Build

1. **Length-1 measure and bracket.**
   `mu1[λ1][u] := Nn1[λ1] / Product[(u - θ1 - h k), {k,0,λ1}]` — Section C's `mu[...][1]` stripped of
   all site-2 factors and the cross-site `(1 - Exp[2πI/h (u-θ2)])`.
   `br1[λ1][f] := 2πI · Plus @@ (Residue[f mu1[λ1][u], {u, #}] & /@ Table[θ1 + h k, {k,0,λ1}])`.
   `Nn1[λ1]` is **derived** (not guessed) by requiring the $u=\theta_1$ residue to contribute with
   coefficient exactly 1, so that `br1[λ1][f] = f(θ1) + (higher-tower corrections)`. Internal check
   (**corrected during execution**): the leading node weight `w_0 = 2πI Res[mu1, θ1] == 1` — **not**
   `br1[λ1][1] == 1`. For a single-site rational measure with $\lambda_1\ge1$ the residues sum to
   zero, so `br1[λ1][1] == 0`; the normalization statement is `w_0 == 1`.
2. **LHS (functional) — the 1×1 FSoV determinant, not the bare bracket.**
   `det1[λ1][f] := br1[λ1][f] / (f /. u -> θ1)` — Section C's `det[λ1,λ2][f]` (which divides by `f`
   at each site's node, `1/(f/.u->θ1)/(f/.u->θ2)/(θ2-θ1)`) restricted to one site.
   `FLdet[λ1,λ2][M,n][k] := det1[λ1][ k^(u/h) · Q1[λ1,λ2][M,n][u] ]`. **The `1/f(θ1) = 1/(k^(θ1/h)
   Q1(θ1))` normalization is essential and was the bug in the first execution pass** (which used the
   bare `br1[...]`, leaving a `k^(θ1/h)` prefactor and a state-dependent `Q1(θ1)` factor that
   prevented the dictionary from closing). With it, the `k^(θ1/h)` cancels and `FLdet` is a clean
   degree-$\lambda_1$ **polynomial in `k`** with constant term exactly 1:
   `FLdet = Σ_{j=0}^{λ1} w_j k^j Q1(θ1+jh)/Q1(θ1) = 1 + (higher k)`.
3. **RHS (operatorial).** With `φ` symbolic:
   `g[λ1][φ] := IdentityMatrix[dim] + Sum[φ^m/m! MatrixPower[Ee[λ1][a,b], m], {m,1,λ1}]` for the
   chosen generator $(a,b)\in\{(1,2),(2,1)\}$ — the explicit nilpotent series (exact, polynomial in
   `φ` of degree $\le\lambda_1$). **Note (gotcha found during execution):** do **not** use
   `MatrixExp` here, and do **not** include the `m=0` term inside the `Sum` — `MatrixPower[nilpotent,
   0]` fires `MatrixPower::sing` and stays unevaluated, so the identity term must be split off
   explicitly. `rx[λ1][φ] := xSingle[λ1][{0}] . g[λ1][φ]`;
   `OP[λ1,λ2][M,n][φ] := KroneckerProduct[rx[λ1][φ], xSingle[λ2][{0}]] . Psi[λ1,λ2][M,n]`, a
   polynomial in `φ` of degree $\le\lambda_1$.

### Discover the dictionary

- **Generator choice:** evaluate `OP` for both $E_{12}$ and $E_{21}$; the lowest-weight covector
  $\langle x_0|$ annihilates one of them (trivial φ-dependence). Keep the nontrivial one and record
  which it is.
- **$\varphi(k)$ and scalar:** at $\lambda_1=1$, `OP` is *linear* in `φ`
  (`OP = a0 + a1 φ`); solve `a0 + a1 φ = FL/scalar` for `φ(k)` and the overall factor. This yields a
  candidate `φ(k)` (expected to involve the `k^(θ1/h)` prefactor and/or `log k`; the exact form is
  whatever the match dictates).
- **Universality validation:** substitute the *same* `φ(k)` into `OP` at $\lambda_1=2,3$ and every
  state, across reps $(\lambda_1,\lambda_2)\in\{(1,1),(2,1),(1,2),(2,2)\}$, and compare against `FL`.
  The dictionary must be a single generator + single `φ(k)` independent of $(M,n)$ and rep. Hard-
  assert this internal consistency once discovered; if `φ(k)` comes out state- or rep-dependent,
  that is reported as a finding (not silently patched).
- **Exact-vs-scalar:** report `FL/OP∘φ` as a function of `k` and state — whether it is `1`, a
  `k`-independent constant, or exactly the `k^(θ1/h)` prefactor — per the "decide from the numbers"
  choice.

### Test case and report

Primary: $(\lambda_1,\lambda_2)=(2,1)$, all 6 states of `states21`
(`{{0,0},{1,0},{1,1},{2,0},{2,1},{3,0}}`), consistent with Sections C/D. Universality sweep across
the reps above. Final report: the discovered dictionary
`{generator, φ(k), scalar(k, state)}`, the per-state `FL` and `OP∘φ`, and the residuals — no assert
against a predicted closed form (there is none yet; this discovery is how one would be found).

## Risks

- **`k^(θ1/h)` branch/fractional-power prefactor.** With $\theta_1=1/3$, `k^(θ1/h)=k^{1/3}` is a
  fractional power that the polynomial-in-`φ` RHS cannot carry. Handle by keeping `k` symbolic and
  factoring `k^(θ1/h)` out explicitly (compare `k^(-θ1/h)·FL` against `OP`), letting the discovery
  absorb it into `φ(k)`/scalar cleanly.
- **`Nn1` derivation.** The normalization must be solved so the leading term is exactly `f(θ1)`; a
  wrong `Nn1` silently rescales everything. Guarded by the `br1[λ1][1] == 1` internal check.
- **`Residue` needs simple, distinct poles** at $\theta_1+h k$ — true for the fixed numeric
  $\theta_1=1/3, h=1$; same reasoning already implicit in `FSoV_L2.wb`/Section C.
- **Part 2 is genuinely exploratory.** If the generator/`φ(k)` matching does not close, Phase 2
  reports the mismatch (per-state ratios, best-fit `φ`) rather than asserting a pass — the
  correspondence itself is the object under investigation, not a known identity. Only the *internal*
  consistency (universal `φ(k)` once found) is hard-asserted.
- **Kernel-state hazard (known wolfbook gotcha).** The shared kernel may have lost Section A–D
  definitions; the executor must re-run (not edit) the existing cells to rebuild state before
  Section E/F accessors are trusted, and check `wolfbook_getKernelState` first.

## Out of scope

- Verifying the $s_2\neq 0$ (site-2) factorization or a fully symmetric two-site factorization —
  Part 1 is the $s_2=0$ slice only.
- Writing any of this up in `paper/xxxCG.tex` (still single-site); that is a later step once the
  dictionary is established.
- Comparing `FL`/the CG overlap to a predicted closed-form CG formula — none exists yet.
- Sweeping representations beyond $\{1,2,3\}$ (Part 1) / the four reps above (Part 2).

## Result (executed 2026-07-16)

**Part 1 — factorization confirmed, with the site-2 scalar fully characterized.** The ratio
`x2[{s1,0}] / (xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}])` is bit-exactly constant over all 27
$(\lambda_1,\lambda_2,s_1)$ cases. The site-2 pseudovacuum is a genuine left eigencovector of the
site-2 factor of `t2[1,s1][θ1]` with eigenvalue
$\mu_{s_1}(\lambda_2)=\prod_{k=1}^{s_1}(\theta_1-\theta_2+(k{-}1)h-\lambda_2 h)$, giving
$c(\lambda_1,\lambda_2,s_1)=\prod_{k=1}^{s_1}\frac{\theta_1-\theta_2+(k-1)h-\lambda_2 h}{\theta_1-\theta_2+(k-1)h-\lambda_1 h}$
(=1 at $\lambda_1=\lambda_2$). My original guess that the site-2 factor passes through as the
identity ($c$ = $\lambda_2$-independent) was wrong; it carries this $Q_\theta$-type eigenvalue.
Committed in Section E (cells 53–60), commit `8dbed52`.

**Part 2 — the CG-overlap dictionary closes exactly.** With the corrected 1×1 determinant `det1`
(the bare-bracket `FL` was the bug):
$$
\det{}_1^{[\lambda_1]}\!\big[k^{u/h}\,Q_1[\lambda_1,\lambda_2][M,n]\big]
= \big(\langle x_0^{[\lambda_1]}|\,e^{-k\,E_{21}} \otimes \langle x_0^{[\lambda_2]}|\big)
\cdot Psi[\lambda_1,\lambda_2][M,n].
$$
Dictionary: **generator $E_{21}$, $\varphi(k)=-k$** (α = −1, universal), normalization = the
1-site FSoV determinant. `FLdet[λ1,λ2][M,n][k] − OPg[…][−k]` reduces to the **literal symbolic 0**
across the full rep sweep $\{(1,1),(2,1),(1,2),(2,2),(3,1),(1,3)\}$ (41 states). Section F cells
72–75 (72 `det1`/`FLdet`, 73 λ1=1 universality, 74 full-sweep hard assert, 75 report); cell 71
annotated as resolved.

**Two `Risks`-section notes are now stale** given the above: the "`k^(θ1/h)` prefactor" risk is
resolved *by using `det1`* (the prefactor cancels, no manual factoring needed), and the "`Nn1`…
`br1[λ1][1]==1`" guard should read `w_0==1` per the Build correction above.
