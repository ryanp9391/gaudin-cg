# XXX CG L=3 Notebook Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: none of the standard superpowers execution skills
> apply here — this plan is executed entirely inside a live Wolfram notebook via the project's
> **wolfbook-builder** subagent (see `CLAUDE.md`), not via pytest/git-commit-per-task. Each task
> below specifies notebook cells to insert (via `mcp__wolfbook__wolfbook_insertCells`), run (via
> `mcp__wolfbook__wolfbook_runCell`), and verify (via `mcp__wolfbook__wolfbook_getCellOutput` /
> `wolfbook_getKernelState`), against an exact expected result. Steps use checkbox (`- [ ]`) syntax
> for tracking. There is no separate test file and no per-task git commit — see Global Constraints.

**Goal:** Build `Paul/Mathematica/XXX/Clean/XXX_CG_L3.wb` (currently an empty file) into a
three-site generalization of `XXX_CG_L2_Clean.wb`: the SoV construction (GT basis through
covector basis) and the Baxter TQ-equation solver (`TauEigensystem`/`QSolve`/`Q1`), for
independent gl(2) reps $[\lambda_1,0],[\lambda_2,0],[\lambda_3,0]$ with independent
inhomogeneities $\theta_1,\theta_2,\theta_3$ and a shared companion twist $G$.

**Architecture:** Single notebook, two sections mirroring L2's structure (Section A: GT
basis/generators through SoV covector basis; Section B: Baxter/TQ solver), built directly cell by
cell (no merge from an Experiments/ predecessor this time — see spec).

**Tech Stack:** Wolfram Language, via wolfbook MCP tools operating on the live `.wb` notebook and
kernel.

## Global Constraints

- Numeric parameters: reuse $\theta_1=1/3$, $\theta_2=1/7$, $h=1$, $\kappa_1=N[\mathrm{Exp}[I\,
  \mathrm{Zeta}[3]],30]$, $\kappa_2=1/\kappa_1$ from `XXX_CG_L2_Clean.wb` Cell 1 verbatim. New
  $\theta_3=4/5$ (chosen distinct from $\theta_1,\theta_2$ modulo integers: $4/5-1/3=7/15$,
  $4/5-1/7=23/35$, neither an integer).
  - Numeric parameter reuse comes directly from the approved design; **not** from independent
    verification within this plan — Task 1's own step re-confirms the non-integer-difference
    check mechanically before anything downstream depends on it.
- No dropped validation: every check present in `XXX_CG_L2_Clean.wb` gets an L3 counterpart in
  this notebook (there is no separate `Experiments/` predecessor this time to hold them instead —
  per the approved design, this Clean notebook carries its own full validation suite).
- Hard-assert discipline: every structural check aborts (`Abort[]`) on failure, exactly as in L2.
  If any check in this plan fails when actually run, **stop and investigate — do not weaken the
  tolerance or delete the check to make it pass.** In particular, Task 5's asymptotic cross-check
  and Task 12's end-to-end TQ residual check are the two independent confirmations that this
  plan's re-derived L3 Baxter structure (see Task 10 preamble) is actually correct; a failure
  there means the physics derivation needs revisiting, not the code.
- `wolfbook-builder` stays out of `paper/*.tex` and git commits/pushes (per `CLAUDE.md`) — this
  plan produces no git commits; saving the notebook (`wolfbook_fileOps`) is the only persistence
  step.
- Comment convention: `(* ... *)` cells explain *why*, not what — matches the rest of the
  notebook's cells and the project convention in `CLAUDE.md`.
- Memoized routines (`ps`, `dim`, `T`, `x3`, `t3Coeffs`, `TauEigensystem`, `QSolve`) must have
  their defining cell (and everything downstream) re-run after any fix — stale memoized values
  will not self-clear.

---

## Key derivation this plan encodes (read before starting)

At L=2, the transfer matrix $t_{1,1}(u)$ is quadratic in $u$: two genuinely operator-valued
coefficients ($C_0,C_1$) plus a trivial leading term ($C_2\propto\mathbb{1}$). `TauEigensystem`
reads the weight label $M$ off $t_1$ (the $C_1$-expectation) via a closed formula, and the single
remaining number $J=t_0$ is read directly per state.

At L=3, $t_{1,1}(u)$ is **cubic**: $C_3\propto\mathbb{1}$ (trivial), $C_2$ is a pure
Cartan/total-generator quantity (same style of formula as L2's $C_1$, one order down since the
polynomial is one degree higher — this is proven by the Task 5 asymptotic cross-check before
anything downstream trusts it), but $C_1$ and $C_0$ **do not** reduce to closed forms: expanding
the monodromy $T(u)=L_1(u)L_2(u)L_3(u)$ shows the order-$h^2$ term (which feeds $C_1$) mixes
different sites' raising/lowering operators pairwise ($E^{(k)}E^{(l)}$, $k\ne l$), which is not a
function of the total weight $M$ alone. So L3 needs **two** free per-state numbers, $(J_1,J_0)$,
generalizing L2's single $J$: $\tau(u)=\tau_0(u)+J_1 u+J_0$, with $\tau_0(u)$ containing only the
$u^3,u^2$ terms (by construction, exactly as L2's $\tau_0$ has zero $u^0$ term by construction).

This changes three routines relative to L2:
- `TauEigensystem` jointly diagonalizes a **two**-parameter combination $C_0+r_1C_1+r_2C_2$
  (L2 needed one parameter), extracts $(t_0,t_1,t_2)$, reads $M$ from $t_2$, and records
  $(J_1,J_0)=(t_1,t_0)$ directly (no formula — same "read it straight off the eigenvector" pattern
  L2 used for $J$).
- `BaxterMatrix` subtracts only $\tau_0(u)$ (which has zero $u^1,u^0$ coefficients) when building
  the state-independent linear operator. It still needs exactly **two** high-degree vanishing
  checks, same count as L2 — but shifted up one degree ($M{+}2,M{+}3$, not L2's $M{+}1,M{+}2$)
  since $a(u),Q_\theta(u),\tau_0(u)$ are now degree 3. The row directly above that, degree
  $M{+}1$, is **not** a pure per-column identity at L3 (unlike L2): writing
  $O_0[u^k]=O_{\mathrm{full}}[u^k]+J_1u^{k+1}+J_0u^k$ (relating the state-independent $O_0$, built
  from $\tau_0$ only, to the true $O_{\mathrm{full}}$, built from the full $\tau=\tau_0+J_1u+J_0$),
  the $J_1$ term contaminates row $p=k+1$ — which reaches exactly row $M{+}1$ from the top column
  $k=M$. Rows $M{+}2,M{+}3$ are never reached this way ($k{+}1$ tops out at $M{+}1$ for $k\le M$),
  so they stay clean identities. Row $M{+}1$'s physical content is exactly what `QSolve`'s
  $J_1\cdot(\text{shiftDown})$ correction supplies later — it must not be asserted to vanish in
  `BaxterMatrix` — and its correctness is confirmed instead by Task 12's end-to-end TQ residual
  check (which validates every degree, including $M{+}1$, numerically at sample points). This
  generalizes cleanly: for any chain length $L$, exactly the **top two** extra-degree
  coefficients ($M{+}(L{-}1)$ and $M{+}L$) are always pure per-column identities, regardless of
  $L$; any lower extra-degree row is contaminated by one of the $L{-}1$ free state-dependent
  parameters and belongs to `QSolve`, not `BaxterMatrix`.
- `QSolve` subtracts **both** $J_1\cdot(\text{subdiagonal shift matrix})$ and $J_0\cdot\mathbb{1}$
  from `BaxterMatrix` before the nullspace solve (L2 only needed $J\cdot\mathbb{1}$), because
  $J_1 u^{k+1}$ (from $\tau(u)u^k$) lands one row below the diagonal relative to column $k$.

Everything else (Section A in full; `qdetT`/fusion hierarchy; `dPred`; `Q1`) generalizes
mechanically — the algebraic form doesn't reference $L$ explicitly, only the domain
($\lambda_1,\lambda_2,\lambda_3$ instead of $\lambda_1,\lambda_2$) changes.

---

### Task 1: Notebook scaffold and numeric parameters

**Files:**
- Modify: `Paul/Mathematica/XXX/Clean/XXX_CG_L3.wb` (currently empty — this task creates its
  first two cells)

**Interfaces:**
- Produces: global numerics `θ1,θ2,θ3,h,κ1,κ2,χ1,χ2` used by every later task.

- [ ] **Step 1: Insert Cell 1 (numeric parameters)**

Insert this cell as the notebook's first cell:

```wolfram
(*Numeric parameters (shared across the whole notebook). NOTE: this assigns global numeric values
  to θ1,θ2,θ3,h,κ1,κ2. χ1,χ2 are the companion twist G's eigenvalues (Section A's G={{χ1,-χ2},
  {1,0}} has char. poly λ²-χ1 λ+χ2=(λ-κ1)(λ-κ2), so χ1=κ1+κ2, χ2=κ1κ2=1 here since κ2=1/κ1) —
  fixed here, numerically, once for the whole notebook. θ1,θ2,h,κ1,κ2 are reused verbatim from
  XXX_CG_L2_Clean.wb. θ3 is a new rational inhomogeneity for the third site, chosen distinct from
  θ1,θ2 modulo integers (checked below) to avoid an accidental Baxter-equation degeneracy.*)
θ1 = 1/3;
θ2 = 1/7;
θ3 = 4/5;
h = 1;
κ1 = N[Exp[I Zeta[3]], 30];
κ2 = 1/κ1;
χ1 = κ1 + κ2;
χ2 = 1;
{θ1, θ2, θ3, h, κ1, κ2, χ1, χ2}
```

- [ ] **Step 2: Run Cell 1 and verify**

Run the cell. Expected: no error; output shows the 8-element list with θ3=4/5 and κ1 a 30-digit
complex number on the unit circle.

- [ ] **Step 3: Insert and run a scratch check that θ3 doesn't differ from θ1,θ2 by an integer**

```wolfram
{IntegerQ[θ3 - θ1], IntegerQ[θ3 - θ2]}
```

Expected: `{False, False}`. If either is `True`, stop — pick a different θ3 and redo Step 1.

- [ ] **Step 4: Insert Cell 2 (chain description comment)**

```wolfram
(*Three-site XXX_{su(2)} chain (L=3). Sites carry independent gl(2) representations [λ1,0],
  [λ2,0], [λ3,0], independent inhomogeneities θ1, θ2, θ3, and a SHARED companion twist G=(χ1,χ2)
  for the whole chain. Physical space V_λ1⊗V_λ2⊗V_λ3. Generator/GT-basis machinery (ps,dim,J0,
  Jp,Jm,JJ) and Ee, L are ported verbatim from XXX_su2.wb — they are single-site building blocks,
  independent of L. Structural difference from the L=2 notebook: t3[1,1](u) is CUBIC in u (three
  Lax factors, not two), which changes the Baxter/TauEigensystem construction in Section B — see
  the design note at the top of this notebook's implementation plan for the derivation.*)
```

Run it (comment-only cell, no output expected beyond `Null`).

- [ ] **Step 5: Save the notebook**

---

### Task 2: GT basis, single-site generators, single-site Lax matrix

**Files:**
- Modify: `Paul/Mathematica/XXX/Clean/XXX_CG_L3.wb` (Cells 3-5)

**Interfaces:**
- Consumes: nothing beyond Task 1's params (this block is representation-theory machinery,
  independent of chain length).
- Produces: `ps[rep]`, `dim[rep]`, `J0`/`Jp`/`Jm`/`JJ`, `Ee[λ][i,j]`, `L[λ][i,j][u]` — used by
  every later task.

- [ ] **Step 1: Insert Cell 3 (GT basis / generators), verbatim from `XXX_CG_L2_Clean.wb` Cell 3**

```wolfram
ps[rep_] :=
  ps[rep] =
   Block[{res = (Prepend[rep, #] & /@ (Tuples[
          Table[Range[rep[[1, i + 1]], rep[[1, i]]], {i, 1,
            Length[rep[[1]]] - 1}]]))},
    If[Length[res[[1, 1]]] == 1, Reverse@Sort@res,
     Reverse@Sort@Flatten[ps /@ res, 1]]];
dim[rep_] := dim[rep] = Length[ps[{rep}]];
J0[k_][pat_] :=
  J0[k][pat] =
   Sum[pat[[k, i]], {i, 1, k}] - Sum[pat[[k - 1, i]], {i, 1, k - 1}];
Jp[k_][pat_] :=
  Jp[k][pat] = -Sum[
     Product[pat[[k, i]] - i - pat[[k + 1, j]] + j, {j,
         Range[k + 1]}]/
       Product[pat[[k, i]] - i - pat[[k, j]] + j, {j,
         Complement[Range[k], {i}]}] v[
       ReplacePart[pat, {k, i} -> pat[[k, i]] + 1]], {i, 1, k}];
Jm[k_][pat_] :=
  Jm[k][pat] =
   Sum[Product[
       pat[[k, i]] - i - pat[[k - 1, j]] + j, {j, Range[k - 1]}]/
      Product[pat[[k, i]] - i - pat[[k, j]] + j, {j,
        Complement[Range[k], {i}]}] v[
      ReplacePart[pat, {k, i} -> pat[[k, i]] - 1]], {i, 1, k}];
JJ[k_, l_, rep_] :=
  JJ[k, l, rep] =
   Which[k == l, DiagonalMatrix[J0[k] /@ ps[{rep}]], l == k + 1,
    Transpose[
     Table[Coefficient[Jp[k][ps[{rep}][[i]]], v[#]] & /@
       ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l == k - 1,
    Transpose[
     Table[Coefficient[Jm[l][ps[{rep}][[i]]], v[#]] & /@
       ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l > k + 1,
    JJ[k, k + 1, rep] . JJ[k + 1, l, rep] -
     JJ[k + 1, l, rep] . JJ[k, k + 1, rep], l < k - 1,
    JJ[k, k - 1, rep] . JJ[k - 1, l, rep] -
     JJ[k - 1, l, rep] . JJ[k, k - 1, rep]];
```

- [ ] **Step 2: Run Cell 3, then verify with a scratch check**

```wolfram
{dim[{1, 0}], dim[{2, 0}], dim[{3, 0}]}
```

Expected: `{2, 3, 4}`.

- [ ] **Step 3: Insert and run Cell 4 (`Ee`)**

```wolfram
(*spin-s generators*)
ClearAll[Ee];
Ee[λ_][i_,j_]:=JJ[i,j,{λ,0}]
```

- [ ] **Step 4: Insert and run Cell 5 (single-site Lax matrix)**

```wolfram
(*Lax matrix for spin-s representation (no θ baked in — shift applied when a site is instantiated)*)
ClearAll[L];
L[λ_][i_,j_][u_]:=KroneckerDelta[i,j]*u*IdentityMatrix[dim[{λ,0}]]-h Ee[λ][j,i];
```

Verify with a scratch check:

```wolfram
Dimensions[L[2][1, 1][u]]
```

Expected: `{3, 3}` (matching `dim[{2,0}]=3`).

- [ ] **Step 5: Save the notebook**

---

### Task 3: Companion twist, three-site Lax operators, monodromy

**Files:**
- Modify: `Paul/Mathematica/XXX/Clean/XXX_CG_L3.wb` (Cells 6-10)

**Interfaces:**
- Consumes: `dim`, `L` from Task 2; `θ1,θ2,θ3,h,χ1,χ2` from Task 1.
- Produces: `G`, `L1`/`L2`/`L3`, `T[λ1,λ2,λ3][i,j][u]`.

- [ ] **Step 1: Insert and run Cell 6 (companion twist)**

```wolfram
(*companion twist matrix — shared by the whole chain*)
G={{χ1,-χ2},{1,0}};
```

- [ ] **Step 2: Insert Cell 7 (comment) and Cell 8 (three-site Lax operators)**

```wolfram
(*Three-site Lax operators on V_λ1⊗V_λ2⊗V_λ3, built as KroneckerProduct of the single-site Lax
  matrix (evaluated at its own site's inhomogeneity shift) with the identity on the other two
  factors. Signature L1[λ1,λ2,λ3][i,j][u], etc. — all three λ's must be supplied to each routine
  even though only one representation enters the shifted Lax matrix itself, because the size of
  the identity factors depends on the OTHER two sites' representations.*)
```

```wolfram
ClearAll[L1,L2,L3];
L1[λ1_,λ2_,λ3_][i_,j_][u_]:=KroneckerProduct[L[λ1][i,j][u-θ1],IdentityMatrix[dim[{λ2,0}]],IdentityMatrix[dim[{λ3,0}]]];
L2[λ1_,λ2_,λ3_][i_,j_][u_]:=KroneckerProduct[IdentityMatrix[dim[{λ1,0}]],L[λ2][i,j][u-θ2],IdentityMatrix[dim[{λ3,0}]]];
L3[λ1_,λ2_,λ3_][i_,j_][u_]:=KroneckerProduct[IdentityMatrix[dim[{λ1,0}]],IdentityMatrix[dim[{λ2,0}]],L[λ3][i,j][u-θ3]];
```

Run both cells. Verify with a scratch check:

```wolfram
Dimensions[L1[1, 1, 1][1, 1][u]]
```

Expected: `{8, 8}` (dim{1,0}=2 for each of 3 sites, 2×2×2=8).

- [ ] **Step 3: Insert Cell 9 (comment) and Cell 10 (monodromy)**

```wolfram
(*Monodromy matrix T[λ1,λ2,λ3][i,j][u] = product of the three site Lax operators in auxiliary
  space, T(u)=L1(u).L2(u).L3(u). All three factors carry the SAME spectral parameter u — each
  site's inhomogeneity shift is already built into L1/L2/L3 via θ1,θ2,θ3 respectively.*)
```

```wolfram
ClearAll[T];
T[λ1_,λ2_,λ3_][i_,j_][u_]:=T[λ1,λ2,λ3][i,j][u]=Sum[L1[λ1,λ2,λ3][i,k][u].L2[λ1,λ2,λ3][k,l][u].L3[λ1,λ2,λ3][l,j][u],{k,2},{l,2}];
```

Run it. Verify:

```wolfram
Dimensions[T[1, 1, 1][1, 1][u]]
```

Expected: `{8, 8}`.

- [ ] **Step 4: Save the notebook**

---

### Task 4: Transfer matrix and commutativity check

**Files:**
- Modify: `Paul/Mathematica/XXX/Clean/XXX_CG_L3.wb` (Cells 11-12)

**Interfaces:**
- Consumes: `T`, `G` from Task 3.
- Produces: `t3[λ1,λ2,λ3][1,1][u]`.

- [ ] **Step 1: Insert and run Cell 11 (transfer matrix)**

```wolfram
(*transfer matrix in fundamental representation*)
ClearAll[t3];
t3[λ1_,λ2_,λ3_][1,1][u_]:=Sum[T[λ1,λ2,λ3][i,j][u]G[[j,i]],{i,2},{j,2}];
```

- [ ] **Step 2: Insert and run Cell 12 (commutativity check)**

```wolfram
(*transfer matrices commute at different values of spectral parameter*)
(*so coefficients generate a complete family of integrals of motion. Range kept to {1,2} (not
  L2's {1,...,4}) to bound compute cost — dimension grows as (λ+1)^3 at three sites.*)
Table[t3[λ1,λ2,λ3][1,1][u].t3[λ1,λ2,λ3][1,1][v]-t3[λ1,λ2,λ3][1,1][v].t3[λ1,λ2,λ3][1,1][u],{λ1,1,2},{λ2,1,2},{λ3,1,2}]//Simplify//Flatten//Union
```

Expected: `{0}`. If not, stop — this means the `T`/`t3` construction from Task 3 has a bug; do not
proceed to Task 5 until this passes.

- [ ] **Step 3: Save the notebook**

---

### Task 5: Diagonal-twist transfer matrix, total generators, central charge, asymptotic cross-check

**Files:**
- Modify: `Paul/Mathematica/XXX/Clean/XXX_CG_L3.wb` (Cells 13-17)

**Interfaces:**
- Consumes: `T`, `Ee`, `dim` from earlier tasks.
- Produces: `t3diag`, `EE[λ1,λ2,λ3][i,j]`, `t3diagasympt` — and, most importantly, **numerical
  confirmation that $\tau_0(u)$'s $u^2$ coefficient formula (used in Task 8 onward) is correct.**
  This is the load-bearing validation for the whole Baxter re-derivation described at the top of
  this plan.

- [ ] **Step 1: Insert and run Cell 13 (diagonal-twist transfer matrix)**

```wolfram
(*transfer matrix in fundamental representation with diagonal twist*)
t3diag[λ1_,λ2_,λ3_][1,1][u_]:=Sum[T[λ1,λ2,λ3][i,i][u]κ[i],{i,2}];
```

Note: `κ[1]`,`κ[2]` here are formal symbolic indexed twist eigenvalues, **distinct** from the
numeric plain `κ1`,`κ2` fixed in Task 1 — same convention/gotcha as `XXX_CG_L2_Clean.wb` (see its
"Known gotchas" section). They only appear in this asymptotics cross-check.

- [ ] **Step 2: Insert and run Cell 14 (global/total generators)**

```wolfram
(*global Lie algebra generators (triple coproduct)*)
ClearAll[EE]
EE[λ1_,λ2_,λ3_][i_,j_]:=KroneckerProduct[Ee[λ1][i,j],IdentityMatrix[dim[{λ2,0}]],IdentityMatrix[dim[{λ3,0}]]]+KroneckerProduct[IdentityMatrix[dim[{λ1,0}]],Ee[λ2][i,j],IdentityMatrix[dim[{λ3,0}]]]+KroneckerProduct[IdentityMatrix[dim[{λ1,0}]],IdentityMatrix[dim[{λ2,0}]],Ee[λ3][i,j]];
```

- [ ] **Step 3: Insert and run Cell 15 (central charge check)**

```wolfram
(*central charge*)
Table[EE[λ1,λ2,λ3][1,1]+EE[λ1,λ2,λ3][2,2]-(λ1+λ2+λ3)IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]dim[{λ3,0}]],{λ1,1,3},{λ2,1,3},{λ3,1,3}]//Flatten//Union
```

Expected: `{0}`. If not, stop — `EE` is wrong; fix before Step 4.

- [ ] **Step 4: Insert and run Cell 16 (predicted asymptotics)**

```wolfram
(*Predicted leading-order asymptotics of the diagonal-twist transfer matrix, generalizing L2's
  tdiagasympt: the leading (u^3) term is the trivial (κ1+κ2)Id-times-Π(u-θi) coefficient, and the
  u^2 coefficient is a pure Cartan/total-generator quantity -h[κ(1)E11tot+κ(2)E22tot]. This is the
  ONLY coefficient below leading that reduces to a c-number formula times a global generator —
  lower orders genuinely mix different sites' raising/lowering operators and are NOT predicted
  here (see Task 10's TauEigensystem, where they become the free per-state parameters J1,J0
  instead of a closed form).*)
t3diagasympt[λ1_,λ2_,λ3_][u_]:=(κ[1]+κ[2])u^3 IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]dim[{λ3,0}]] - u^2(κ[1]+κ[2])(θ1+θ2+θ3) IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]dim[{λ3,0}]] - u^2 h(κ[1] EE[λ1,λ2,λ3][1,1]+κ[2] EE[λ1,λ2,λ3][2,2]);
```

- [ ] **Step 5: Insert and run Cell 17 (the load-bearing cross-check)**

```wolfram
(*checking asymptotics matches at u^2 order (analog of L2's u^1 check, one order higher since
  L3's τ(u) is cubic not quadratic). THIS IS THE KEY VALIDATION for the Baxter re-derivation: if
  this does not collapse to zero, the τ0 formula used from Task 8 onward is wrong and must be
  re-derived before proceeding — do not patch downstream code around a failure here.*)
t3diag[1,1,1][1,1][u]-t3diagasympt[1,1,1][u];
Coefficient[%,u,2]//Simplify
```

Expected: the zero matrix (all entries `0`). **If this is nonzero, stop and re-derive — this
validates the central mathematical claim this whole plan depends on.**

- [ ] **Step 6: Save the notebook**

---

### Task 6: Quantum determinant and fusion/Hirota hierarchy

**Files:**
- Modify: `Paul/Mathematica/XXX/Clean/XXX_CG_L3.wb` (Cells 18-19)

**Interfaces:**
- Consumes: `T` from Task 3, `χ2` from Task 1.
- Produces: `qdetT`, `t3[λ1,λ2,λ3][a,s][u]` for general `a,s` (fused transfer matrices).

- [ ] **Step 1: Insert and run Cell 18 (quantum determinant)**

```wolfram
(*quantum determinant*)
ClearAll[qdetT];
qdetT[λ1_,λ2_,λ3_][u_]:=T[λ1,λ2,λ3][1,1][u].T[λ1,λ2,λ3][2,2][u-h]-T[λ1,λ2,λ3][2,1][u].T[λ1,λ2,λ3][1,2][u-h];
```

- [ ] **Step 2: Insert and run Cell 19 (fusion/Hirota hierarchy)**

```wolfram
(*Fusion/Hirota hierarchy for the three-site transfer matrices t3[λ1,λ2,λ3][1,s][u] — algebraic
  form IDENTICAL to the L2 notebook (references T(u)/t3(u) only, not L explicitly); only the
  domain (λ1,λ2,λ3 a triple, not a pair) changes. No standalone check cell here — correctness is
  implicit in the SoV covector eigenvalue check in Task 7, matching the L2 notebook's convention.*)
t3[λ1_,λ2_,λ3_][0,s_][u_]:=IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]dim[{λ3,0}]];
t3[λ1_,λ2_,λ3_][a_,0][u_]:=IdentityMatrix[dim[{λ1,0}]dim[{λ2,0}]dim[{λ3,0}]];
t3[λ1_,λ2_,λ3_][2,1][u_]:=χ2 qdetT[λ1,λ2,λ3][u];
t3[λ1_,λ2_,λ3_][1,s_][u_]:=t3[λ1,λ2,λ3][1,s][u]=t3[λ1,λ2,λ3][1,s-1][u].t3[λ1,λ2,λ3][1,1][u+(s-1)h]-t3[λ1,λ2,λ3][2,1][u+(s-1)h].t3[λ1,λ2,λ3][1,s-2][u];
```

Run it (definitions only, no output expected beyond `Null`).

- [ ] **Step 3: Save the notebook**

---

### Task 7: SoV B-operator, pseudovacuum/covector basis, eigenvalue check (closes Section A)

**Files:**
- Modify: `Paul/Mathematica/XXX/Clean/XXX_CG_L3.wb` (Cells 20-25)

**Interfaces:**
- Consumes: `T`, `t3`, `dim` from earlier tasks.
- Produces: `Qθ[u]`, `B3[λ1,λ2,λ3][u]`, `x3[λ1,λ2,λ3][{s1,s2,s3}]`.

- [ ] **Step 1: Insert and run Cell 20 (`Qθ`)**

```wolfram
(*Qθ factor for the three-site chain. (L2 also kept a single-site Qθ1 factor for parity with its
  source file even though it was never used downstream — dropped here as unused/YAGNI.)*)
ClearAll[Qθ];
Qθ[u_]:=(u-θ1)(u-θ2)(u-θ3);
```

- [ ] **Step 2: Insert Cell 21 (comment) and Cell 22 (`B3`)**

```wolfram
(*SoV B-operator for L=3: the bare (1,1) monodromy element (no G, no trace) — same convention as
  B2[λ1,λ2][u]=T[λ1,λ2][1,1][u] at L=2.*)
```

```wolfram
ClearAll[B3];
B3[λ1_,λ2_,λ3_][u_]:=T[λ1,λ2,λ3][1,1][u];
```

Run both cells.

- [ ] **Step 3: Insert Cell 23 (comment) and Cell 24 (pseudovacuum/covector basis)**

```wolfram
(*SoV pseudovacuum for L=3 as the tensor product of the three single-site lowest-weight
  covectors.*)
```

```wolfram
ClearAll[xSingle,x3];
xSingle[λ_][{0}]:={Table[KroneckerDelta[ii,dim[{λ,0}]],{ii,1,dim[{λ,0}]}]};
x3[λ1_,λ2_,λ3_][{0,0,0}]:=KroneckerProduct[xSingle[λ1][{0}],xSingle[λ2][{0}],xSingle[λ3][{0}]];

x3[λ1_,λ2_,λ3_][{s1_,s2_,s3_}] := x3[λ1,λ2,λ3][{s1,s2,s3}] =
  x3[λ1,λ2,λ3][{0,0,0}].t3[λ1,λ2,λ3][1,s1][θ1].t3[λ1,λ2,λ3][1,s2][θ2].t3[λ1,λ2,λ3][1,s3][θ3]/
  (Product[Qθ[θ1 + (-λ1 + k - 1) h], {k, 1, s1}] Product[Qθ[θ2 + (-λ2 + k - 1) h], {k, 1, s2}] Product[Qθ[θ3 + (-λ3 + k - 1) h], {k, 1, s3}]);
```

Run both cells.

- [ ] **Step 4: Insert and run Cell 25 (eigenvalue check — closes Section A)**

```wolfram
(*checking x3 diagonalizes B3 — λ-range {1,2} (not L2's {1,2,3}) to bound compute cost at three
  sites (max dim (λ+1)^3 grows fast); still covers the general/unequal-λ case.*)
Table[Table[x3[λ1,λ2,λ3][{s1,s2,s3}].B3[λ1,λ2,λ3][u]-(u-θ1-h s1)(u-θ2-h s2)(u-θ3-h s3)x3[λ1,λ2,λ3][{s1,s2,s3}],{s1,0,λ1},{s2,0,λ2},{s3,0,λ3}],{λ1,1,2},{λ2,1,2},{λ3,1,2}]//Simplify//Flatten//Union
```

Expected: `{0}`. If not, stop — this exercises the fusion hierarchy from Task 6 as well as `x3`
itself; investigate before proceeding to Section B.

- [ ] **Step 5: Save the notebook**

---

### Task 8: Baxter ingredients and `BaxterMatrix` (opens Section B)

**Files:**
- Modify: `Paul/Mathematica/XXX/Clean/XXX_CG_L3.wb` (Cells 26-28)

**Interfaces:**
- Consumes: `Qθ` from Task 7; `θ1,θ2,θ3,h,κ1,κ2` from Task 1.
- Produces: `aFun[λ1,λ2,λ3][u]`, `τ0[λ1,λ2,λ3,M][u]`, `qpoly[clist,uu]`, `BaxterMatrix[λ1,λ2,λ3,M]`.

- [ ] **Step 1: Insert Cell 26 (Section B intro comment)**

```wolfram
(*Section B — Baxter (TQ) solver for L=3. Structural difference from L=2, validated in Task 5's
  asymptotic cross-check above: t3[1,1](u) is CUBIC in u, so τ0(u) below carries only its u^3,u^2
  terms (state-independent), and TWO free per-state numbers (J1,J0) are needed instead of L2's
  single J — see TauEigensystem (Task 10) and QSolve (Task 11).*)
```

Run it.

- [ ] **Step 2: Insert and run Cell 27 (Baxter ingredients)**

```wolfram
(*a(u), Qθ(u) (triple product; reuses the Qθ defined in Task 7 — no redefinition needed here since
  Section B was never split from Section A into a separately self-contained file for L3), and
  τ0 = the state-independent part of τ (using ee1=λ1+λ2+λ3-M, ee2=M). τ0 has ONLY u^3,u^2 terms by
  construction — its u^1,u^0 coefficients are deliberately zero; the actual eigenvalue is
  τ=τ0+J1 u+J0 with (J1,J0) read directly per state in TauEigensystem (Task 10).
  qpoly builds a polynomial from an ASCENDING coefficient list {c0,...,cM} (cM=1 for monic q_n);
  formed in a dummy symbol first so substituting uu=0 gives the constant term rather than an
  indeterminate 0^0.*)
ClearAll[aFun, τ0, qpoly];
aFun[λ1_, λ2_, λ3_][u_] := (u - θ1 - h λ1) (u - θ2 - h λ2) (u - θ3 - h λ3);
τ0[λ1_, λ2_, λ3_, M_][u_] := (κ1 + κ2) u^3 - (θ1 + θ2 + θ3) (κ1 + κ2) u^2 - u^2 h (κ1 (λ1 + λ2 + λ3 - M) + κ2 M);
qpoly[clist_, uu_] := Module[{x}, (clist . x^Range[0, Length[clist] - 1]) /. x -> uu];
{aFun[1, 1, 1][u], Qθ[u], τ0[1, 1, 1, 1][u]}
```

Verify: output shows `aFun` as a degree-3 polynomial, `Qθ[u]` as `(u-θ1)(u-θ2)(u-θ3)` (numeric
θ's substituted), and `τ0[1,1,1,1][u]` with only `u^3` and `u^2` terms (no `u^1` or constant term).

- [ ] **Step 3: Insert and run Cell 28 (`BaxterMatrix`)**

```wolfram
(*BaxterMatrix: build O[u^k]=κ1 a(u)(u+h)^k + κ2 Qθ[u](u-h)^k - τ0(u) u^k for k=0..M.
  a(u),Qθ(u),τ0(u) are degree 3 (one more than L2's degree 2), so O[u^k] has nonzero coefficients
  up to degree k+3 — but only the TOP TWO extra-degree coefficients, degree M+2 and M+3, are pure
  per-column identities (same count as L2, just shifted up one degree). Degree M+1 is deliberately
  NOT checked here: it is contaminated by the free parameter J1 (from τ(u)=τ0(u)+J1 u+J0) via its
  u^(k+1) shift landing exactly at the top column k=M — that row's physical content is supplied
  later by QSolve's J1·shiftDown correction, not asserted here, and is confirmed instead by the
  end-to-end TQ residual check (Task 12). HARD ASSERT that the degree M+2, M+3 coefficients vanish
  for every column (this is where a wrong a/Qθ/τ0 or a wrong ee1=λ1+λ2+λ3-M convention would first
  show up). Then columns = degree-0..M coefficients.*)
ClearAll[BaxterMatrix];
BaxterMatrix[λ1_, λ2_, λ3_, M_] := Module[{clists, high, maxhigh, cols},
  clists = Table[
     PadRight[
      CoefficientList[
       Expand[κ1 aFun[λ1, λ2, λ3][u] (u + h)^k + κ2 Qθ[u] (u - h)^k - τ0[λ1, λ2, λ3, M][u] u^k],
       u],
      M + 4],
     {k, 0, M}];
  high = clists[[All, {M + 3, M + 4}]];
  maxhigh = Max[Abs[Flatten[high]]];
  If[maxhigh > 10^-18,
    Print["BaxterMatrix ASSERTION FAILED (λ1,λ2,λ3,M)=", {λ1, λ2, λ3, M},
          "  max|high-degree coeff|=", maxhigh, "  high=", high];
    Abort[]];
  cols = clists[[All, 1 ;; M + 1]];
  Transpose[cols]
];
Dimensions /@ (BaxterMatrix[1, 1, 1, #] & /@ {0, 1, 2, 3})
```

Expected: `{{1,1},{2,2},{3,3},{4,4}}`, no assertion failure printed.

- [ ] **Step 4: Save the notebook**

---

### Task 9: `t3Coeffs` and `dPred`

**Files:**
- Modify: `Paul/Mathematica/XXX/Clean/XXX_CG_L3.wb` (Cells 29-30)

**Interfaces:**
- Consumes: `t3` from Task 4; `κ1,κ2` from Task 1; `dim` from Task 2.
- Produces: `t3Coeffs[λ1,λ2,λ3]` (returns `{C0,C1,C2,C3}`), `dPred[M,λ1,λ2,λ3]`.

- [ ] **Step 1: Insert and run Cell 29 (`t3Coeffs`)**

```wolfram
(* Coefficient matrices C0,C1,C2,C3 with t3[λ1,λ2,λ3][1,1][u] = C0 + C1 u + C2 u^2 + C3 u^3.
   Extracted entrywise via CoefficientList in a dummy variable (NOT by evaluating at numeric
   points). Memoized per (λ1,λ2,λ3) — reusable across all (M,n) of that representation. *)
ClearAll[t3Coeffs];
t3Coeffs[λ1_, λ2_, λ3_] := t3Coeffs[λ1, λ2, λ3] = Module[{uu, mat, cl},
   mat = t3[λ1, λ2, λ3][1, 1][uu];
   cl = Map[PadRight[CoefficientList[#, uu], 4] &, mat, {2}];
   Table[cl[[All, All, k]], {k, 1, 4}]
];
(* check: C3 (leading, index 4 in this 1-indexed list) must equal (κ1+κ2) Id — leading-order
   asymptotics match τ by construction *)
With[{Cs = t3Coeffs[1, 1, 1]}, "max|C3-(κ1+κ2)Id|" -> Max[Abs[Flatten[Cs[[4]] - (κ1 + κ2) IdentityMatrix[8]]]]]
```

Expected: the max-abs-difference is `~0` (at the 30-digit noise floor, e.g. `~10^-27` or smaller).

- [ ] **Step 2: Insert and run Cell 30 (`dPred`)**

```wolfram
(* dPred[M,λ1,λ2,λ3]: predicted dimension of the weight-M subspace of Vλ1⊗Vλ2⊗Vλ3 (number of
   (m1,m2,m3) triples with m1+m2+m3=M, 0<=mi<=λi) — a direct triple-convolution count (no closed
   form; small ranges make the brute double-sum cheap). This is what the GENUINE
   t3-diagonalization below (Task 10) must reproduce exactly for every M. *)
ClearAll[dPred];
dPred[M_, λ1_, λ2_, λ3_] := Sum[If[0 <= M - m1 - m2 <= λ3, 1, 0], {m1, 0, λ1}, {m2, 0, λ2}];
{Sum[dPred[M, 1, 1, 1], {M, 0, 3}], dim[{1, 0}] dim[{1, 0}] dim[{1, 0}]}
```

Expected: `{8, 8}` (the weight-sector counts sum to the total representation dimension).

- [ ] **Step 3: Save the notebook**

---

### Task 10: `TauEigensystem` (the re-derived core)

**Files:**
- Modify: `Paul/Mathematica/XXX/Clean/XXX_CG_L3.wb` (Cell 31)

**Interfaces:**
- Consumes: `t3Coeffs` and `dPred` from Task 9; `θ1,θ2,θ3,h,κ1,κ2` from Task 1.
- Produces: `TauEigensystem[λ1,λ2,λ3]`, a memoized `Association` keyed by `{M,n}` giving
  `<|"J1"->J1_n, "J0"->J0_n, "Psi"->eigenvector|>`.

- [ ] **Step 1: Insert and run Cell 31**

```wolfram
(*TauEigensystem[λ1,λ2,λ3]: diagonalize the GENUINE transfer matrix t3[λ1,λ2,λ3][1,1][u] directly.
  C0,C1,C2,C3 pairwise commute (coefficients of the commuting family t3(u)). C3 is already known to
  be the scalar (κ1+κ2) Id. Unlike L2 (where only C0,C1 carried information, C1 alone determining
  M), here C2 is the analog of L2's C1 — proven in Task 5's asymptotic cross-check to be a pure
  Cartan/total-generator quantity equal to -h[κ1(Λ-M)+κ2 M] plus known θ/κ terms — so M is read off
  ALGEBRAICALLY from t2 (the C2-expectation) exactly as L2 read M from t1. C1 and C0, by contrast,
  genuinely mix different sites' raising/lowering operators (see the plan's derivation note) and
  are NOT predicted by any closed form — they are read off DIRECTLY per eigenvector as (J1,J0), the
  L3 generalization of L2's single J=t0. This is why the diagonalization below uses a
  TWO-parameter generic combination C0+r1 C1+r2 C2 (L2 needed only one parameter, r).
  M-extraction formula (generalizing L2's Mextract by one order, using t2 in place of t1):
    M = [-(θ1+θ2+θ3)(κ1+κ2) - hκ1(λ1+λ2+λ3) - t2] / [h(κ2-κ1)]
  hard-asserted within 10^-10 of an integer. (J1,J0)=(t1,t0) recorded directly, no formula. Each
  M-sector's count is hard-asserted against dPred[M,λ1,λ2,λ3].*)
ClearAll[TauEigensystem];
TauEigensystem::rankdef = "Combination C0+r1 C1+r2 C2 failed to resolve the `1`-dim eigenspace for (λ1,λ2,λ3)=(`2`,`3`,`4`) at every tried (r1,r2) in `5`.";
TauEigensystem::noninteger = "M = `1` is not within tolerance of an integer for eigenvector #`2` at (λ1,λ2,λ3)=(`3`,`4`,`5`) (combo (r1,r2)=`6`) — indicates a bug in the τ0 u^2-coefficient formula or the t3 construction.";
Psi::notgenuine = "Eigenvector #`1` for (λ1,λ2,λ3)=(`2`,`3`,`4`) failed the common-eigenvector residual check against C0/C1/C2 (residual=`5`, combo (r1,r2)=`6`) — the generic-combination diagonalization should make this impossible; investigate before trusting results.";
Psi::degnull = "Weight sector M=`1` for (λ1,λ2,λ3)=(`2`,`3`,`4`) has `5` genuine state(s) but the predicted weight multiplicity dPred is `6`.";
Psi::lastzero = "Genuine eigenvector for (M,n)=(`1`,`2`) at (λ1,λ2,λ3)=(`3`,`4`,`5`) has last component ~0; cannot normalize.";
TauEigensystem[λ1_, λ2_, λ3_] := TauEigensystem[λ1, λ2, λ3] = Module[
   {d, Cs, C0, C1, C2, C3, rCandidates, r1, r2, vals, vecs, recs, byM, table},
   d = (λ1 + 1) (λ2 + 1) (λ3 + 1);
   Cs = t3Coeffs[λ1, λ2, λ3]; {C0, C1, C2, C3} = Cs;
   rCandidates = {{7/3, 11/5}, {N[Pi, 30], N[E, 30]}, {13/4, N[Sqrt[2], 30]}, {17/6, 19/8}, {N[Sqrt[3], 30], 23/9}, {29/11, 31/13}};
   {r1, r2} = SelectFirst[rCandidates,
     Module[{ee, vv, rr1, rr2}, {rr1, rr2} = #;
       {ee, vv} = Eigensystem[C0 + rr1 C1 + rr2 C2];
       MatrixRank[vv] == d && Length[Union[Chop[ee, 10^-12]]] == d] &,
     Missing["NoGenericR"]];
   If[MissingQ[r1], Message[TauEigensystem::rankdef, d, λ1, λ2, λ3, rCandidates]; Abort[]];
   {vals, vecs} = Eigensystem[N[C0 + r1 C1 + r2 C2, 30]];
   recs = Table[
     Module[{v = vecs[[i]], t0, t1, t2v, res0, res1, res2, Mi},
       t0 = (v . (C0 . v))/(v . v);
       t1 = (v . (C1 . v))/(v . v);
       t2v = (v . (C2 . v))/(v . v);
       res0 = Max[Abs[C0 . v - t0 v]];
       res1 = Max[Abs[C1 . v - t1 v]];
       res2 = Max[Abs[C2 . v - t2v v]];
       If[Max[res0, res1, res2] > 10^-15,
         Message[Psi::notgenuine, i, λ1, λ2, λ3, Max[res0, res1, res2], {r1, r2}]; Abort[]];
       Mi = (-(θ1 + θ2 + θ3) (κ1 + κ2) - h κ1 (λ1 + λ2 + λ3) - t2v)/(h (κ2 - κ1));
       If[Abs[Mi - Round[Re[Mi]]] > 10^-10,
         Message[TauEigensystem::noninteger, Mi, i, λ1, λ2, λ3, {r1, r2}]; Abort[]];
       <|"M" -> Round[Re[Mi]], "J1" -> t1, "J0" -> t0, "vec" -> v|>],
     {i, 1, d}];
   byM = GroupBy[recs, #["M"] &];
   table = Association @ Flatten @ KeyValueMap[
     Function[{Mval, group},
       Module[{sorted = SortBy[group, {Re[#["J0"]], Im[#["J0"]], Re[#["J1"]], Im[#["J1"]]} &], predicted = dPred[Mval, λ1, λ2, λ3]},
         If[Length[sorted] != predicted,
           Message[Psi::degnull, Mval, λ1, λ2, λ3, Length[sorted], predicted]; Abort[]];
         MapIndexed[
           Function[{rec, idx}, Module[{n = First[idx] - 1, vv = rec["vec"], last},
             last = vv[[-1]];
             If[Abs[last] < 10^-15 Max[Abs[vv]],
               Message[Psi::lastzero, Mval, n, λ1, λ2, λ3]; Abort[]];
             {Mval, n} -> <|"J1" -> rec["J1"], "J0" -> rec["J0"], "Psi" -> vv/last|>]],
           sorted]]],
     byM];
   If[Length[table] != d, Message[Psi::degnull, "total", λ1, λ2, λ3, Length[table], d]; Abort[]];
   table
];
(*smoke test: (λ1,λ2,λ3)=(1,1,1) table keys and J0 values*)
KeySort[TauEigensystem[1, 1, 1]][[All, "J0"]] // N
```

Expected: no `Abort[]`, no error messages; output is an 8-entry association-slice of numeric `J0`
values (one per state, since `dim{1,0}dim{1,0}dim{1,0}=8`). If any hard assert fires, **stop** —
this means either Task 5's cross-check formula doesn't actually hold for this case (re-examine
Task 5) or there's a bug in this cell's transcription of that formula.

- [ ] **Step 2: Save the notebook**

---

### Task 11: `Psi`/`τ` accessors, `QSolve`, `Q1`

**Files:**
- Modify: `Paul/Mathematica/XXX/Clean/XXX_CG_L3.wb` (Cells 32-34)

**Interfaces:**
- Consumes: `TauEigensystem`, `dPred` from Task 10/9; `BaxterMatrix` from Task 8; `τ0`, `qpoly`
  from Task 8; `κ1,h` from Task 1.
- Produces: `Psi[λ1,λ2,λ3][M,n]`, `τ[λ1,λ2,λ3][M,n][u]`, `QSolve[λ1,λ2,λ3][M,n]`,
  `Q1[λ1,λ2,λ3][M,n][u]`.

- [ ] **Step 1: Insert and run Cell 32 (`Psi`/`τ` accessors)**

```wolfram
(* Psi[λ1,λ2,λ3][M,n]: on-demand accessor over the memoized TauEigensystem cache. Returns the
   fixed (last-component-normalized) eigenvector for a genuine state.
   τ[λ1,λ2,λ3][M,n][u]: the FULL τ eigenvalue polynomial for that state, τ0[λ1,λ2,λ3,M][u] +
   J1_n u + J0_n (generalizing L2's τ0[...][u]+J_n by the extra J1 term — see the plan's
   derivation note). Bounds/failure modes match the L2 accessor style (bnds/spurious). *)
ClearAll[Psi, τ];
Psi::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): need M an integer in 0..λ1+λ2+λ3 and n an integer in 0..M.";
τ::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): need M an integer in 0..λ1+λ2+λ3 and n an integer in 0..M.";
Psi::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): weight sector M=`1` has only `6` genuine state(s) (n=0..`7`), so n=`2` does not exist.";
τ::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): weight sector M=`1` has only `6` genuine state(s) (n=0..`7`), so n=`2` does not exist.";
Psi[λ1_, λ2_, λ3_][M_, n_] := Module[{mult},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 + λ3 && 0 <= n <= M),
    Message[Psi::bnds, M, n, λ1, λ2, λ3]; Return[$Failed]];
  mult = dPred[M, λ1, λ2, λ3];
  If[n > mult - 1,
    Message[Psi::spurious, M, n, λ1, λ2, λ3, mult, mult - 1]; Return[Missing["NoState", {λ1, λ2, λ3, M, n}]]];
  TauEigensystem[λ1, λ2, λ3][{M, n}, "Psi"]];
τ[λ1_, λ2_, λ3_][M_, n_][u_] := Module[{mult, rec},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 + λ3 && 0 <= n <= M),
    Message[τ::bnds, M, n, λ1, λ2, λ3]; Return[$Failed]];
  mult = dPred[M, λ1, λ2, λ3];
  If[n > mult - 1,
    Message[τ::spurious, M, n, λ1, λ2, λ3, mult, mult - 1]; Return[Missing["NoState", {λ1, λ2, λ3, M, n}]]];
  rec = TauEigensystem[λ1, λ2, λ3][{M, n}];
  τ0[λ1, λ2, λ3, M][u] + rec["J1"] u + rec["J0"]];
Print["spurious-guard test (expect a Psi::spurious message and Missing[...]): ", Psi[2, 1, 1][3, 3]];
{Psi[1, 1, 1][1, 0], N[τ[1, 1, 1][1, 0][2], 10]}
```

Expected: a printed `Psi::spurious` message for the deliberate `Psi[2,1,1][3,3]` smoke test
(`dPred[3,2,1,1]=3`, so only `n=0,1,2` exist — `n=3` is gross-bound-valid but exceeds the actual
multiplicity), followed by a `Missing[...]` result; then a valid `{Psi[...], τ value}` pair for
`(M,n)=(1,0)`.

- [ ] **Step 2: Insert and run Cell 33 (`QSolve`)**

```wolfram
(*QSolve[λ1,λ2,λ3][M,n]: linear solve, generalizing L2. τ(u)=τ0(u)+J1 u+J0 now has TWO free
  state-dependent parameters (not L2's single J), so the known-J shift subtracted from
  BaxterMatrix is TWO terms: J0·Identity (same diagonal-shift role as L2's J·Identity) PLUS
  J1·shiftDown, where shiftDown is the (M+1)x(M+1) matrix with 1's on the subdiagonal (row=col+1,
  1-indexed) — this represents multiplying the column's underlying u^k monomial by an extra u^1,
  i.e. "J1 u^(k+1)" landing one row below the diagonal relative to column k. Same
  Chop-before-NullSpace precision-tracked-zero fix as L2 (see L2's QSolve comment for the full
  numerical-precision rationale — unchanged here).*)
ClearAll[QSolve];
QSolve::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): need M an integer in 0..λ1+λ2+λ3 and n an integer in 0..M.";
QSolve::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2,λ3)=(`3`,`4`,`5`): weight sector M=`1` has only `6` genuine state(s) (n=0..`7`), so n=`2` does not exist.";
QSolve::degnull = "Unexpected Q-function nullspace dimension `6` (expected 1) for (M,n)=(`1`,`2`) at (λ1,λ2,λ3)=(`3`,`4`,`5`).";
QSolve::lastzero = "Q-function nullspace vector found but its top (u^M) coefficient is ~0 for (M,n)=(`1`,`2`) at (λ1,λ2,λ3)=(`3`,`4`,`5`); cannot normalize to monic.";
QSolve[λ1_, λ2_, λ3_][M_, n_] := QSolve[λ1, λ2, λ3][M, n] = Module[{mult, rec, J1n, J0n, shiftDown, mat, ns, vec, last},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 + λ3 && 0 <= n <= M),
    Message[QSolve::bnds, M, n, λ1, λ2, λ3]; Return[$Failed]];
  mult = dPred[M, λ1, λ2, λ3];
  If[n > mult - 1,
    Message[QSolve::spurious, M, n, λ1, λ2, λ3, mult, mult - 1]; Return[$Failed]];
  rec = TauEigensystem[λ1, λ2, λ3][{M, n}];
  J1n = rec["J1"]; J0n = rec["J0"];
  shiftDown = Table[If[i == j + 1, 1, 0], {i, 1, M + 1}, {j, 1, M + 1}];
  mat = BaxterMatrix[λ1, λ2, λ3, M] - J1n shiftDown - J0n IdentityMatrix[M + 1];
  ns = NullSpace[Chop[mat, 10^-12], Tolerance -> 10^-10];
  Which[
    Length[ns] != 1, Message[QSolve::degnull, M, n, λ1, λ2, λ3, Length[ns]]; Return[$Failed],
    True,
      vec = ns[[1]]; last = vec[[-1]];
      If[Abs[last] < 10^-15 Max[Abs[vec]],
        Message[QSolve::lastzero, M, n, λ1, λ2, λ3]; Return[$Failed]];
      vec/last]
];
QSolve[1, 1, 1][2, 0]
```

Expected: no error message, output is a length-3 list of numeric coefficients (ascending, monic —
last entry `1.000...+0.×10^-N I`).

- [ ] **Step 3: Insert and run Cell 34 (`Q1`)**

```wolfram
(* Q1[λ1,λ2,λ3][M,n][u] = κ1^(u/h) q_n(u): the actual Baxter Q-function value, built from the
   SAME QSolve cache (mutually consistent with τ by construction). *)
ClearAll[Q1];
Q1[λ1_, λ2_, λ3_][M_, n_][u_] := Module[{q = QSolve[λ1, λ2, λ3][M, n]},
  If[q === $Failed, $Failed, κ1^(u/h) qpoly[q, u]]];
{Q1[1, 1, 1][1, 0][2], N[τ[1, 1, 1][1, 0][2], 10]}
```

Expected: a numeric complex pair, no error.

- [ ] **Step 4: Save the notebook**

---

### Task 12: End-to-end TQ residual regression check, then final save

**Files:**
- Modify: `Paul/Mathematica/XXX/Clean/XXX_CG_L3.wb` (Cell 35)

**Interfaces:**
- Consumes: `aFun`, `Qθ` (Task 8/7), `Q1`, `τ` (Task 11), `dPred` (Task 9).
- Produces: nothing new — this is a pure regression check, the final validation of the whole
  notebook.

- [ ] **Step 1: Insert and run Cell 35**

```wolfram
(*End-to-end Baxter TQ-equation residual check: a(u)Q(u+h) - τ(u)Q(u) + Qθ(u)Q(u-h) = 0, assembled
  purely from Q1 and τ (which trace back to different derivations — Q1 through
  QSolve/BaxterMatrix, τ through TauEigensystem's C0/C1/C2 diagonalization — so a small residual
  here is strong evidence both derivations, including the re-derived L3 M-extraction/J1,J0
  structure, are correct). Sample points avoid roots of a(u) and Qθ(u). Checked for every (M,n) at
  (λ1,λ2,λ3)=(1,1,1) and one asymmetric case (2,1,1), mirroring the L2 notebook's regression
  cases.*)
checkTQ[λ1_, λ2_, λ3_, M_, n_] := Max[Abs[Table[
    aFun[λ1, λ2, λ3][u] Q1[λ1, λ2, λ3][M, n][u + h] - τ[λ1, λ2, λ3][M, n][u] Q1[λ1, λ2, λ3][M, n][u] + Qθ[u] Q1[λ1, λ2, λ3][M, n][u - h],
    {u, {0, 1, 2, -1, 5}}]]];
Table[{M, n, checkTQ[1, 1, 1, M, n]}, {M, 0, 3}, {n, 0, dPred[M, 1, 1, 1] - 1}] // Flatten[#, 1] &
```

Expected: a list of `{M, n, residual}` triples with every `residual` at or below `~10^-18`
(matching the tracked 30-digit numeric precision budget — see L2's `~10^-20` to `~10^-26` range
for comparison). **If any residual is not tiny (e.g. `O(1)`), stop — do not proceed to save a
notebook with a failing end-to-end check.**

- [ ] **Step 2: Insert and run the asymmetric regression case**

```wolfram
Table[{M, n, checkTQ[2, 1, 1, M, n]}, {M, 0, 4}, {n, 0, dPred[M, 2, 1, 1] - 1}] // Flatten[#, 1] &
```

Expected: same as Step 1 — all residuals tiny.

- [ ] **Step 3: Fresh-kernel full run**

Restart the kernel (`wolfbook_kernelControl`) and run every cell top-to-bottom in order. Confirm
zero errors and zero `Abort[]`/hard-assert messages anywhere in the run.

- [ ] **Step 4: Save the notebook**

Save `Paul/Mathematica/XXX/Clean/XXX_CG_L3.wb`. This is the final step — no git commit (per Global
Constraints; `wolfbook-builder` doesn't commit).

---

## Self-Review Notes

- **Spec coverage:** Task 1 covers parameters; Tasks 2-7 cover the spec's "Section A" scope;
  Tasks 8-12 cover "Section B" including the flagged `TauEigensystem`/`BaxterMatrix`/`QSolve`
  re-derivation; Task 12 covers the spec's "no check gets dropped" requirement via the end-to-end
  regression check (which in L2 lived only in the now-nonexistent `Experiments/` predecessor).
  `dim{2,1,1}` and `dim{1,1,1}` test cases match the spec's named representative cases.
- **Placeholder scan:** no TBD/TODO; every step has complete Wolfram code and an exact expected
  result.
- **Type consistency:** `t3Coeffs` returns `{C0,C1,C2,C3}` (checked against Task 10's destructuring
  `{C0,C1,C2,C3}=Cs`); `TauEigensystem`'s returned association keys (`"J1"`,`"J0"`,`"Psi"`) match
  what Task 11's `Psi`/`τ`/`QSolve` read from it.
