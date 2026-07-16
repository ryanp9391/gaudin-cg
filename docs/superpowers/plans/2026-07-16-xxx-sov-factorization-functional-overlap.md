# L=2 SoV-basis factorization and functional CG overlap — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. The executor is the **wolfbook-builder** subagent, working live in the notebook via the wolfbook MCP tools.

**Goal:** Add two sections to `Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb` — Section E verifying the operatorial factorization `x2[λ1,λ2][{s1,0}] = c·(xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}])`, and Section F discovering the dictionary relating the length-1 bracket `br1[k^(u/h) Q1]` to the operatorial CG overlap `(⟨x0|g(k)⊗⟨x0|)·Psi`.

**Architecture:** All new cells append after Section D (last cell ~index 52) of `XXX_CG_L2_V2.wb`. Section E builds a genuine standalone length-1 transfer-matrix hierarchy `t1`/`qdet1` and covector `xL1`, validated by an independent B-eigenvector check, then verifies the factorization up-to-scalar. Section F builds the length-1 measure/bracket `mu1`/`br1`, the functional object `FL`, the operatorial overlap `OP`, and runs an exploratory dictionary-discovery. Sections A–D are read-only.

**Tech Stack:** Wolfram Language in a `.wb` (wolfbook) notebook; wolfbook MCP tools for cell insert/edit/run/save; exact numeric globals (θ1=1/3, θ2=1/7, h=1, κ1=N[Exp[I Zeta[3]],30], κ2=1/κ1, χ1=κ1+κ2, χ2=1) at 30-digit working precision.

## Global Constraints

- **Notebook:** `Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb`. Append only; do **not** edit Section A–D cells (except re-running them to rebuild kernel state — see Task 0).
- **Numeric globals** (already set in cell 0, verbatim): `θ1 = 1/3; θ2 = 1/7; h = 1; κ1 = N[Exp[I Zeta[3]], 30]; κ2 = 1/κ1; χ1 = κ1 + κ2; χ2 = 1;`
- **Existing symbols reused (do not redefine):** `dim[{λ,0}]`, `Ee[λ][i,j]`, `L[λ][i,j][u]` (no θ baked in), `G` (companion twist matrix), `xSingle[λ][{0}]` (1×(λ+1) row covector), `x2[λ1,λ2][{s1,s2}]`, `Qθ1[u]=u-θ1`, `Qθ[u]=(u-θ1)(u-θ2)`, `Psi[λ1,λ2][M,n]` (length-D right eigenvector, last component 1), `Q1[λ1,λ2][M,n][u] = κ1^(u/h) qpoly[q,u]`, `statesOf[λ1,λ2]`, `dPred[M,λ1,λ2]`.
- **Assert style:** match the notebook's existing `Message[name::tag,...]; Abort[]` hard-assert convention (as in `TauEigensystem`, `QSolve`, Section C's `FSoVGram`). Noise floor for genuine-zero at 30-digit precision is ~10^-20..10^-27; the hard-assert threshold convention is **10^-15**.
- **Comments:** every definition cell opens with a `(* ... *)` comment explaining *why*, per repo convention. Multi-statement cells use `;` between statements.
- **Persistence:** cell edits/evaluations do **not** persist until the notebook is explicitly saved (`wolfbook_fileOps` save). Save at the end of every task.
- **Memoization hazard:** `f[x_]:=f[x]=...` caches. After any fix, re-run the defining cell and everything downstream.
- **wolfbook-builder does not commit to git.** Git commits happen at phase boundaries, performed by the orchestrator/user, not the executor. "Save the notebook" replaces "commit" in each task.

## Refinements over the spec (apply these; they are more correct)

1. The LHS exponential base is **`k·κ1`**, not `k` alone: `Q1` already carries `κ1^(u/h)`, so `k^(u/h) Q1 = (k κ1)^(u/h) q_n(u)`, giving an overall `(k κ1)^(θ1/h)` prefactor on the pole tower. The dictionary `φ(k)` absorbs it.
2. The `Nn1` internal check is **leading-node weight = 1** (`w_0 = 2πI Res[mu1, θ1] == 1`), not `br1[λ][1]==1` — for a single-site rational measure with λ≥1 the residues sum to zero, so `br1[λ][1] == 0`. `w_0==1` is exactly the "`br1[f]=f(θ1)+…`" normalization.

---

## Phase 1 — Operatorial factorization (Section E)

### Task 0: Rebuild kernel state and confirm prerequisites

**Files:**
- Inspect: `Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb` (Sections A–D)

**Interfaces:**
- Produces: a live kernel with all Section A–D definitions loaded and verified, ready for Section E/F to build on.

- [ ] **Step 1: Check kernel state**

Call `wolfbook_getKernelState` / `wolfbook_getNotebookContext`. If the notebook is not the active target, `wolfbook_setTarget` to `XXX_CG_L2_V2.wb`.

- [ ] **Step 2: Re-run Sections A–D top to bottom**

Run (do **not** edit) every existing cell from cell 0 through the last Section D cell, in order. This rebuilds `dim`, `Ee`, `L`, `G`, `xSingle`, `x2`, `Psi`, `Q1`, `statesOf`, `dPred`, `Qθ1`, `Qθ`, Section C `br`/`det`/`G21`/`states21`, Section D `PsiLeft`. (Known gotcha: the shared kernel may have lost definitions; earlier sessions had to rebuild this way.)

- [ ] **Step 3: Verify prerequisites resolve**

Run this scratch check (via `wolfbook_evaluateExpression`, not a saved cell):

```wolfram
{
 Dimensions[xSingle[2][{0}]],            (* expect {1,3} *)
 Dimensions[x2[2,1][{1,0}]],             (* expect {1,6} *)
 Length[Psi[2,1][1,0]],                  (* expect 6 *)
 Q1[2,1][1,0][2] =!= $Failed,            (* expect True *)
 Qθ1[u], Qθ[u],                          (* expect u-1/3 ; (u-1/3)(u-1/7) *)
 statesOf[2,1]                           (* expect {{0,0},{1,0},{1,1},{2,0},{2,1},{3,0}} *)
}
```

Expected: dimensions `{1,3}` and `{1,6}`, length `6`, `True`, the two Qθ forms, and the 6-state list. If any differ, stop and report — Section E/F cannot proceed.

- [ ] **Step 4: Save the notebook** (no new cells yet; this just checkpoints the rebuilt outputs).

---

### Task 1: Length-1 transfer-matrix hierarchy and SoV covector `xL1`

**Files:**
- Modify: `Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb` (append cells after Section D)

**Interfaces:**
- Consumes: `L[λ][i,j][u]`, `G`, `dim`, `xSingle[λ][{0}]`, `Qθ1`, globals θ1/h/χ2.
- Produces:
  - `t1[λ][a,s][u]` — length-1 fused transfer matrices ((λ+1)×(λ+1) matrices); `t1[λ][1,1]`, `qdet1[λ]`, `t1[λ][2,1]`, Hirota `t1[λ][1,s]`, boundaries `t1[λ][0,s]=t1[λ][a,0]=Id`.
  - `B1[λ][u]` — length-1 SoV B-operator `L[λ][1,1][u-θ1]`.
  - `xL1[λ][{s}]` — genuine length-1 SoV covector, 1×(λ+1) row, paper (1-site) normalization.

- [ ] **Step 1: Insert the Section E header cell** (markdown, `kind:1`)

```
(*Section E — operatorial factorization of the L=2 SoV basis. Builds a genuine STANDALONE
  length-1 (single-site) fused transfer-matrix hierarchy t1[λ][a,s] and SoV covector xL1[λ][{s}]
  (the paper's L=1 construction), then verifies that the s2=0 slice of the L=2 SoV basis factorizes
  as x2[λ1,λ2][{s1,0}] = c (xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}]) up to a scalar c fixed by the
  Qθ-vs-Qθ1 normalization mismatch. Sections A–D are read-only.
  See docs/superpowers/specs/2026-07-16-xxx-sov-factorization-functional-overlap-design.md.*)
```

- [ ] **Step 2: Insert the `t1`/`qdet1` definition cell** (wolfram, `kind:2`)

```wolfram
(*Single-site (L=1) fused transfer matrices — the L=1 analogue of Section A's two-site t2 block,
  built from the SAME single-site Lax L[λ] (no θ baked in) and shared companion twist G. qdet1 is
  the single-site quantum determinant: both factors are this site's L evaluated at u-θ1 and u-θ1-h
  (the internal -h shift), unlike the two-site qdetT where θ is already in T. t1[λ][1,1](u)=tr_a[G
  L(u-θ1)] matches paper eq. t_{1,1}(u)=tr_a[G L(u-θ)].*)
ClearAll[t1, qdet1];
t1[λ_][0, s_][u_] := IdentityMatrix[dim[{λ, 0}]];
t1[λ_][a_, 0][u_] := IdentityMatrix[dim[{λ, 0}]];
t1[λ_][1, 1][u_] := Sum[L[λ][i, j][u - θ1] G[[j, i]], {i, 2}, {j, 2}];
qdet1[λ_][u_] := L[λ][1, 1][u - θ1] . L[λ][2, 2][u - θ1 - h] - L[λ][2, 1][u - θ1] . L[λ][1, 2][u - θ1 - h];
t1[λ_][2, 1][u_] := χ2 qdet1[λ][u];
t1[λ_][1, s_ /; s > 1][u_] := t1[λ][1, s - 1][u] . t1[λ][1, 1][u + (s - 1) h] - t1[λ][2, 1][u + (s - 1) h] . t1[λ][1, s - 2][u];
```

- [ ] **Step 3: Run the cell; confirm no errors**

Run it. Expected: no output/messages (pure definitions). Then scratch-check the shape:

```wolfram
{Dimensions[t1[2][1,1][u]], Dimensions[t1[2][1,2][u]]}
```

Expected: `{{3,3},{3,3}}`.

- [ ] **Step 4: Insert the `B1`/`xL1` definition cell** (wolfram, `kind:2`)

```wolfram
(*Genuine length-1 SoV covector (paper normalization, 1-site Qθ1(u)=u-θ1). xSingle[λ][{0}] is the
  site lowest-weight covector (Section A) — a 1×(λ+1) row. xL1 stays that shape: (1×n).(n×n)=1×n.
  B1 is the L=1 SoV B-operator.*)
ClearAll[B1, xL1];
B1[λ_][u_] := L[λ][1, 1][u - θ1];
xL1[λ_][{0}] := xSingle[λ][{0}];
xL1[λ_][{s_}] := xL1[λ][{s}] = xSingle[λ][{0}] . t1[λ][1, s][θ1] / Product[Qθ1[θ1 - λ h + (p - 1) h], {p, 1, s}];
```

- [ ] **Step 5: Run the cell; confirm no errors** (pure definitions, no output expected).

- [ ] **Step 6: Insert the B-eigenvector internal-check cell** (wolfram, `kind:2`)

```wolfram
(*Internal check (independent of the factorization): xL1[λ][{s}] must be a LEFT eigencovector of
  B1 with eigenvalue u-θ1-s h — the paper's L=1 SoV property (verified there symbolically for
  λ=1..5). This validates the whole t1/xL1 construction before it is used downstream. Symbolic in
  u; residual Chopped at the 10^-20 noise floor.*)
xL1BCheck = Table[
   xL1[λ][{s}] . B1[λ][u] - (u - θ1 - s h) xL1[λ][{s}],
   {λ, 1, 3}, {s, 0, λ}] // ExpandAll // Chop[#, 10^-20] & // Flatten // Union;
xL1BCheck
```

- [ ] **Step 7: Run the check; verify it collapses to `{0}`**

Expected output: `{0}`. If anything else appears, stop — the L=1 construction is wrong; report the offending `{λ,s}` before proceeding.

- [ ] **Step 8: Save the notebook.**

---

### Task 2: Verify the operatorial factorization

**Files:**
- Modify: `Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb` (append after Task 1 cells)

**Interfaces:**
- Consumes: `x2[λ1,λ2][{s1,0}]`, `xL1[λ1][{s1}]`, `xSingle[λ2][{0}]`, globals θ1/θ2/h.
- Produces: `factResults` — table `{λ1,λ2,s1,c,cPred,spread}` over λ1,λ2∈{1,2,3}, s1∈0..λ1; two hard asserts on failure.

- [ ] **Step 1: Insert the factorization-check cell** (wolfram, `kind:2`)

```wolfram
(*Operatorial factorization: x2[λ1,λ2][{s1,0}] = c (xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}]). Flatten both
  to length-D vectors, take the component ratio over the nonzero components of the tensor product.
  HARD requirement (the actual factorization content — "t2[1,s1][θ1] acts only on site 1"): the
  ratio lhs/tp must be a CONSTANT c (spread below the 10^-15 floor). SEPARATELY, compare c to the
  ideal normalization scalar cPred = 1/Product[(θ1-λ1 h+(k-1)h)-θ2,{k,1,s1}] (reciprocal of the
  product of (tower point − θ2) factors: x2 divides by the 2-site Qθ, xL1 only by the 1-site Qθ1).
  cPred assumes the site-2 pseudovacuum passes through the operator as EXACTLY the identity; the
  ratio μEff = c/cPred is REPORTED (soft-warn if ≠1) rather than hard-asserted, because a nontrivial
  μEff would be a genuine physics finding (a site-2 eigenvalue), not a construction bug. cPred is
  independent of λ2, so c matching across λ2∈{1,2,3} is a free extra check. s1=0 ⇒ cPred=1, ratio 1.*)
SoVFact::notconst = "Ratio lhs/tp not constant for (λ1,λ2,s1)=(`1`,`2`,`3`): spread=`4`.";
SoVFact::mu = "μEff=c/cPred ≠ 1 for (λ1,λ2,s1)=(`1`,`2`,`3`): μEff=`4` (site-2 passthrough not identity — a finding, not a failure).";
factResults = Flatten[Table[
    Module[{lhs, tp, nz, rl, cval, cpred, mueff, spread},
     lhs = Flatten[x2[λ1, λ2][{s1, 0}]];
     tp = Flatten[KroneckerProduct[xL1[λ1][{s1}], xSingle[λ2][{0}]]];
     nz = Position[tp, x_ /; Abs[Chop[x, 10^-20]] > 10^-12];
     rl = Extract[lhs, nz]/Extract[tp, nz];
     cval = Mean[rl];
     spread = Max[Abs[rl - cval]];
     cpred = 1/Product[(θ1 - λ1 h + (k - 1) h) - θ2, {k, 1, s1}];
     mueff = cval/cpred;
     If[spread > 10^-15, Message[SoVFact::notconst, λ1, λ2, s1, N[spread]]; Abort[]];
     If[Abs[N[mueff - 1]] > 10^-12, Message[SoVFact::mu, λ1, λ2, s1, N[mueff]]];  (* soft warn, no Abort *)
     {λ1, λ2, s1, cval, cpred, N[mueff], N[spread]}],
    {λ1, 1, 3}, {λ2, 1, 3}, {s1, 0, λ1}], 2];
Column[factResults]
```

- [ ] **Step 2: Run the cell; verify constancy holds and inspect μEff**

Expected: no `SoVFact::notconst` message, no `Abort` — the ratio is constant in every case (`spread` ~10^-20 or smaller). Output rows are `{λ1,λ2,s1,c,cPred,μEff,spread}`. Inspect `μEff`: if the clean passthrough holds it is `1` everywhere (no `SoVFact::mu` warnings) and `c=cPred`; for `s1=0`, `c=cPred=1`. For fixed `(λ1,s1)`, `c` is the same across `λ2∈{1,2,3}` (λ2-independence). **If `SoVFact::mu` warnings appear, record the μEff values — that is a real result (site-2 eigenvalue), and the constancy still confirms the factorization itself.**

- [ ] **Step 3: Insert a compact summary cell** (wolfram, `kind:2`)

```wolfram
(*Summary: worst ratio-constancy spread across the whole sweep (the hard-asserted factorization
  content), and the range of μEff=c/cPred (should be exactly 1 everywhere if the site-2 passthrough
  is trivial). Columns of factResults: {λ1,λ2,s1,c,cPred,μEff,spread}.*)
<|
  "max spread (ratio non-constancy)" -> Max[factResults[[All, 7]]],
  "max |μEff - 1| over sweep" -> Max[Abs[factResults[[All, 6]] - 1]],
  "num (λ1,λ2,s1) cases checked" -> Length[factResults]
|>
```

- [ ] **Step 4: Run it; confirm the headline numbers**

Expected: max spread ~10^-20 or smaller; `max |μEff-1|` either ~0 (clean passthrough, the expected result) or a definite nonzero value (recorded as a finding); count = sum over λ1∈{1,2,3} of (λ1+1)·3 = (2+3+4)·3 = 27.

- [ ] **Step 5: Save the notebook.**

---

### ⛔ CHECKPOINT — end of Phase 1

Do not start Phase 2 until Phase 1 is reviewed and approved. Phase 2 (functional) *uses* the factorized pseudovacuum `⟨x0|g⊗⟨x0|` and the length-1 machinery; it is only meaningful once the factorization is confirmed. The orchestrator commits Phase 1 (the notebook changes for Sections E) at this point.

---

## Phase 2 — Functional CG overlap (Section F)

### Task 3: Length-1 measure and bracket

**Files:**
- Modify: `Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb` (append after Section E)

**Interfaces:**
- Consumes: globals θ1/h.
- Produces:
  - `Nn1[λ]` — measure normalization (derived, not guessed), so the θ1-node weight is 1.
  - `mu1[λ][u]` — length-1 FSoV measure.
  - `br1[λ][f]` — length-1 bracket: `br1[λ][f] = Σ_{j=0}^{λ} w_j f(θ1+h j)` with `w_0 = 1`.

- [ ] **Step 1: Insert the Section F header cell** (markdown, `kind:1`)

```
(*Section F — functional CG overlap. Sets up the length-1 FSoV measure/bracket (Section C's
  mu[...][1] stripped of all site-2 factors and the cross-site (1-Exp[...]) term), normalized so
  br1[λ1][f] = f(θ1) + (higher-tower). Then builds the functional object FL = br1[λ1][k^(u/h) Q1]
  and the operatorial CG overlap OP = (⟨x0|g(k)⊗⟨x0|)·Psi, and DISCOVERS the dictionary
  {generator, φ(k), scalar} relating them (exploratory — no predicted formula asserted).
  NOTE: Q1 carries κ1^(u/h), so k^(u/h) Q1 = (k κ1)^(u/h) q_n(u): the LHS prefactor is (k κ1)^(θ1/h).
  See docs/superpowers/specs/2026-07-16-xxx-sov-factorization-functional-overlap-design.md.*)
```

- [ ] **Step 2: Insert the `mu1`/`br1`/`Nn1` definition cell** (wolfram, `kind:2`)

```wolfram
(*Length-1 FSoV measure and bracket. mu1 is a single-site rational measure with simple poles at
  the tower θ1+h k, k=0..λ. Nn1 is DERIVED (not guessed) by forcing the θ1-node weight
  w_0 = 2πI Res[mu1, θ1] to equal 1, so that br1[λ][f] = f(θ1) + (higher-tower). br1 is the residue
  sum over the whole tower. (For λ≥1 the residues of mu1 sum to 0, so br1[λ][1]=0 — the correct
  normalization statement is w_0=1, checked next.)*)
ClearAll[Nn1, mu1, br1];
Nn1[λ_] := Nn1[λ] = Module[{uu}, 1/(2 π I Residue[1/Product[(uu - θ1 - h k), {k, 0, λ}], {uu, θ1}])];
mu1[λ_][u_] := Nn1[λ]/Product[(u - θ1 - h k), {k, 0, λ}];
br1[λ_][f_] := 2 π I Plus @@ (Residue[f mu1[λ][u], {u, #}] & /@ Table[θ1 + h k, {k, 0, λ}]);
```

- [ ] **Step 3: Run the cell** (pure definitions, no output expected).

- [ ] **Step 4: Insert the normalization internal-check cell** (wolfram, `kind:2`)

```wolfram
(*Internal check: the leading (θ1-node) weight must be exactly 1 for λ=1,2,3 (this IS the
  br1[f]=f(θ1)+... normalization), and br1[λ][1] must be 0 (rational-measure residues sum to
  zero). Also list the full weight vectors w_j for the record.*)
w1[λ_] := Table[2 π I Residue[mu1[λ][u], {u, θ1 + h j}], {j, 0, λ}];
Nn1Check::badlead = "Leading weight w_0 ≠ 1 for λ=`1`: w_0=`2`.";
Nn1Check::badsum = "br1[λ][1] ≠ 0 for λ=`1`: got `2`.";
Do[
  Module[{w = w1[λ], s},
   If[Abs[N[w[[1]] - 1]] > 10^-15, Message[Nn1Check::badlead, λ, w[[1]]]; Abort[]];
   s = br1[λ][1];
   If[Abs[N[s]] > 10^-15, Message[Nn1Check::badsum, λ, s]; Abort[]]],
  {λ, 1, 3}];
<|"w_0 (λ=1,2,3)" -> Table[w1[λ][[1]], {λ, 1, 3}],
  "full weights w_j" -> Table[λ -> w1[λ], {λ, 1, 3}],
  "br1[λ][1] (λ=1,2,3)" -> Table[br1[λ][1], {λ, 1, 3}]|>
```

- [ ] **Step 5: Run it; verify the asserts pass**

Expected: no `Nn1Check::*` message. `w_0 (λ=1,2,3)` is `{1,1,1}`; `br1[λ][1]` is `{0,0,0}`. If either assert fires, stop and report.

- [ ] **Step 6: Insert a bracket sanity cell** (wolfram, `kind:2`)

```wolfram
(*Sanity: br1 acting on a generic analytic f reproduces the weighted tower sum Σ_j w_j f(θ1+h j),
  and picks out f(θ1) as the leading term. Test with f=g[u] (undefined head, stays symbolic).*)
Clear[g];
{br1[1][g[u]] - Sum[w1[1][[j + 1]] g[θ1 + h j], {j, 0, 1}] // Simplify,
 br1[2][g[u]] - Sum[w1[2][[j + 1]] g[θ1 + h j], {j, 0, 2}] // Simplify}
```

- [ ] **Step 7: Run it; verify `{0,0}`** (confirms `br1[λ][f] = Σ_j w_j f(θ1+h j)` with `w_0=1`).

- [ ] **Step 8: Save the notebook.**

---

### Task 4: LHS functional object `FL` and RHS operatorial overlap `OP`

**Files:**
- Modify: `Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb` (append after Task 3)

**Interfaces:**
- Consumes: `br1`, `Q1[λ1,λ2][M,n][u]`, `Ee[λ][a,b]`, `xSingle[λ][{0}]`, `Psi[λ1,λ2][M,n]`, globals.
- Produces:
  - `FL[λ1,λ2][M,n][k]` — length-1 bracket of `k^(u/h) Q1`, a function of symbolic `k`.
  - `gRot[λ][a,b][φ]` — unipotent rotation `MatrixExp[φ Ee[λ][a,b]]`.
  - `OP[λ1,λ2][a,b][M,n][φ]` — CG overlap `Flatten[⟨x0|g ⊗ ⟨x0|]·Psi`, a polynomial in `φ`.
  - `genNontrivial` — the generator `{a,b}∈{{1,2},{2,1}}` that gives nontrivial φ-dependence.

- [ ] **Step 1: Insert the `FL` (LHS) definition cell** (wolfram, `kind:2`)

```wolfram
(*LHS functional object: FL = br1[λ1][k^(u/h) Q1[λ1,λ2][M,n][u]]. Since Q1=κ1^(u/h) q_n(u), the
  exponential base is (k κ1); on the poles u=θ1+h j this is (k κ1)^(θ1/h) (k κ1)^j, so FL is
  (k κ1)^(θ1/h) times a polynomial in k. Keep k symbolic.*)
ClearAll[FL];
FL[λ1_, λ2_][M_, n_][k_] := br1[λ1][k^(u/h) Q1[λ1, λ2][M, n][u]];
(*smoke: one state — expect a (k κ1)^(θ1/h)-type prefactor times a polynomial in k*)
FL[2, 1][1, 0][k] // Simplify
```

- [ ] **Step 2: Run it; inspect the structure**

Expected: a closed-form expression in `k` (numeric κ1, θ1=1/3, h=1). It should carry an overall fractional power `(k κ1)^(1/3)` (or an equivalent `k^(1/3)` times a κ1 constant) multiplying a low-degree polynomial in `k`. Record the expression for Task 5. (No assert here — this is a build+inspect step.)

- [ ] **Step 3: Insert the `gRot`/`OP` (RHS) definition cell** (wolfram, `kind:2`)

```wolfram
(*RHS operatorial CG overlap: rotate the site-1 pseudovacuum covector by g=exp[φ Ee[λ1][a,b]]
  (Ee[a,b] nilpotent of order ≤λ1+1 on the finite-dim rep). Build g as the explicit truncated
  series Σ_{m=0}^{λ} φ^m/m! Ee[a,b]^m — this is EXACT (higher powers vanish) and guaranteed to be a
  polynomial in φ, avoiding MatrixExp returning an unevaluated symbolic exponential. Tensor with the
  site-2 pseudovacuum covector, contract with the two-site right eigenvector Psi. Shapes:
  (1×n1).(n1×n1)=1×n1 → KroneckerProduct with 1×n2 → Flatten to length D → .Psi = scalar.*)
ClearAll[gRot, OP];
gRot[λ_][a_, b_][φ_] := Sum[φ^m/m! MatrixPower[Ee[λ][a, b], m], {m, 0, λ}];
OP[λ1_, λ2_][a_, b_][M_, n_][φ_] := Flatten[KroneckerProduct[xSingle[λ1][{0}] . gRot[λ1][a, b][φ], xSingle[λ2][{0}]]] . Psi[λ1, λ2][M, n];
```

- [ ] **Step 4: Run it** (pure definitions, no output).

- [ ] **Step 5: Insert the generator-selection cell** (wolfram, `kind:2`)

```wolfram
(*Which single-site generator is nontrivial? The lowest-weight covector ⟨x0| annihilates one of
  E12/E21 (giving φ-independent OP); keep the other. Test on (λ1,λ2)=(2,1), state (1,0).*)
op12 = OP[2, 1][1, 2][1, 0][φ] // Simplify;
op21 = OP[2, 1][2, 1][1, 0][φ] // Simplify;
genNontrivial = Which[
   ! FreeQ[op12, φ] && FreeQ[op21, φ], {1, 2},
   ! FreeQ[op21, φ] && FreeQ[op12, φ], {2, 1},
   True, {1, 2}  (* fallback: if both depend on φ, default to {1,2} and flag in the report *)];
<|"OP(E12)" -> op12, "OP(E21)" -> op21, "nontrivial generator {a,b}" -> genNontrivial,
  "both φ-dependent?" -> (! FreeQ[op12, φ] && ! FreeQ[op21, φ])|>
```

- [ ] **Step 6: Run it; record which generator is nontrivial**

Expected: exactly one of `op12`, `op21` depends on `φ`; `genNontrivial` is set to that `{a,b}`. Record it — Task 5 uses `genNontrivial`. If both depend on φ (the "both φ-dependent?" flag is `True`), note it in the report; the discovery in Task 5 still proceeds with `genNontrivial` but this is a finding to surface.

- [ ] **Step 7: Insert the specialized-`OPg` convenience cell** (wolfram, `kind:2`)

```wolfram
(*OPg fixes the nontrivial generator so downstream code reads OPg[λ1,λ2][M,n][φ]. Degree in φ is ≤λ1.*)
ClearAll[OPg];
OPg[λ1_, λ2_][M_, n_][φ_] := OP[λ1, λ2][genNontrivial[[1]], genNontrivial[[2]]][M, n][φ];
(*smoke: degree in φ for a λ1=1 and a λ1=2 case*)
{Exponent[OPg[1, 1][1, 0][φ] // Simplify, φ], Exponent[OPg[2, 1][1, 0][φ] // Simplify, φ]}
```

- [ ] **Step 8: Run it; verify degrees** — expect `{1, 2}` (degree ≤ λ1: linear for λ1=1, quadratic for λ1=2).

- [ ] **Step 9: Save the notebook.**

---

### Task 5: Discover the dictionary and report

**Files:**
- Modify: `Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb` (append after Task 4)

**Interfaces:**
- Consumes: `FL`, `OPg`, `statesOf`, `genNontrivial`, globals.
- Produces: `dictReport` — the discovered `{generator, φ(k), scalar(k,state)}` and per-state residuals. Exploratory: hard-asserts only the *internal* consistency (a universal `φ(k)` once found).

This task is a guided discovery, not a fixed pass/fail. Build the cells in order; at each stage the executor inspects the output and records what it finds. The only hard assert is universality (Step 5). If the correspondence does not close, the executor **reports the obstruction** (per-state ratios, best-fit φ) rather than forcing a pass — the correspondence itself is the object under investigation.

- [ ] **Step 1: Insert the λ1=1 φ-extraction cell** (wolfram, `kind:2`)

Rationale: at λ1=1, `OPg` is linear in φ (`OPg = c0 + c1 φ`). Tentatively assume exact match (scalar 1) and solve `φ_state(k) = (FL(state,k) − c0)/c1` per state. If the resulting `φ_state(k)` is the **same** function across all states of the rep, the dictionary is `{genNontrivial, that φ(k), scalar 1}` at λ1=1.

```wolfram
(*Discovery step 1 (λ1=1, OPg linear in φ): solve φ(k) per state assuming exact match (scalar 1),
  using rep (λ1,λ2)=(1,2) (dim 6, several states). If φ_state(k) coincides across states, the
  dictionary is exact with that φ(k). Report the per-state φ candidates and their pairwise
  differences.*)
lin12[M_, n_] := Module[{c0, c1}, c0 = OPg[1, 2][M, n][0]; c1 = D[OPg[1, 2][M, n][φ], φ]; {c0, c1}];
phiCand12[M_, n_][k_] := Module[{c0, c1}, {c0, c1} = lin12[M, n]; (FL[1, 2][M, n][k] - c0)/c1];
sts12 = statesOf[1, 2];
phiList12 = Table[{s, phiCand12[s[[1]], s[[2]]][k] // Simplify}, {s, sts12}];
phiDiffs12 = Table[phiList12[[i, 2]] - phiList12[[1, 2]] // Simplify, {i, 2, Length[phiList12]}];
<|"states" -> sts12, "φ candidates per state" -> phiList12,
  "differences from state 1 (expect all 0 if universal at scalar 1)" -> phiDiffs12|>
```

- [ ] **Step 2: Run it; inspect whether φ is universal at scalar 1**

Two outcomes, both informative:
- **All differences 0:** φ(k) is universal at scalar 1 — exact match at λ1=1. Set `φuniversal[k_] := phiList12[[1,2]]` in the next cell. Proceed to Step 4.
- **Nonzero differences:** scalar ≠ 1 for some states. Proceed to Step 3 to introduce a per-state scalar.

- [ ] **Step 3 (only if Step 2 showed nonzero differences): Insert the per-state scalar cell** (wolfram, `kind:2`)

Rationale: with a per-state scalar `c_s`, `FL(s,k)/c_s = c0(s) + c1(s) φ(k)`. Two states share one `φ(k)`; eliminate it to solve the `c_s` ratios, then recover a single universal `φ(k)`.

```wolfram
(*Discovery step 1b: allow a per-state scalar c_s. For each state, FL/c_s = c0 + c1 φ(k) ⇒
  φ(k) = (FL/c_s − c0)/c1. Require φ(k) equal across states: pin c_s (up to one overall scale) by
  matching the k-dependence. Solve using the two lowest states, then read off φ(k).*)
{sa, sb} = sts12[[1 ;; 2]];
{c0a, c1a} = lin12[sa[[1]], sa[[2]]]; {c0b, c1b} = lin12[sb[[1]], sb[[2]]];
(*Set c_sa=1 (overall scale convention); solve c_sb so both give the same φ(k).*)
phiA[k_] := (FL[1, 2][sa[[1]], sa[[2]]][k] - c0a)/c1a;
phiFromB[cb_][k_] := (FL[1, 2][sb[[1]], sb[[2]]][k]/cb - c0b)/c1b;
cbSol = Solve[phiA[k] == phiFromB[cb][k] // Simplify, cb];
<|"c_sb solution(s)" -> cbSol, "φ(k) (from state a)" -> (phiA[k] // Simplify)|>
```

Run it; record the per-state scalar `cbSol` and the resulting `φ(k)`. Then define `φuniversal[k_] := phiA[k]//Simplify` for the next step.

- [ ] **Step 4: Insert the universality-validation cell** (wolfram, `kind:2`)

Rationale: the discovered `φuniversal[k]` (and any per-state scalars) must reproduce `FL` from `OPg` at **higher reps** λ1=2,3 and all their states — this is the real test that the dictionary is a single generator + single `φ(k)`.

```wolfram
(*Discovery step 2: validate the universal φ(k) on the higher reps. For each rep and state, compare
  OPg[λ1,λ2][M,n][φuniversal[k]] against FL[λ1,λ2][M,n][k] (up to the per-state scalar found above;
  if scalar was 1, cScale=1). Report the max residual over all states of each rep. cScale[λ1,λ2][M,n]
  defaults to 1 unless Step 3 introduced per-state scalars — in that case reuse the same convention
  (compare FL/cScale to OPg).*)
cScale[λ1_, λ2_][M_, n_] := 1;  (* replace with the per-state scalars from Step 3 if needed *)
repSweep = {{1, 1}, {2, 1}, {1, 2}, {2, 2}, {3, 1}, {1, 3}};
validation = Table[
   Module[{λ1 = rep[[1]], λ2 = rep[[2]], sts, resids},
    sts = statesOf[λ1, λ2];
    resids = Table[
      Module[{lhs, rhs},
       lhs = FL[λ1, λ2][st[[1]], st[[2]]][k]/cScale[λ1, λ2][st[[1]], st[[2]]];
       rhs = OPg[λ1, λ2][st[[1]], st[[2]]][φuniversal[k]];
       {st, Simplify[lhs - rhs]}],
      {st, sts}];
    {λ1, λ2} -> resids],
   {rep, repSweep}];
validation
```

- [ ] **Step 5: Run it; assess universality**

Inspect every residual. Two outcomes:
- **All residuals 0 (to noise):** the dictionary is confirmed universal. Add a hard assert cell (Step 6) and report success.
- **Some residual nonzero:** φ(k) is *not* universal as constructed — a genuine finding. Report which reps/states fail, with the residuals and the per-state φ they *would* need. Do **not** force a pass. This ends the discovery with a reported obstruction.

- [ ] **Step 6: Insert the final dictionary-report cell** (wolfram, `kind:2`)

```wolfram
(*Final report: the discovered dictionary. genNontrivial = the generator; φuniversal[k] = the
  parameter map; cScale = the (state-dependent?) scalar. Hard-assert internal consistency: the max
  |residual| over the whole rep sweep is at the noise floor. If the discovery ended in an
  obstruction (Step 5), replace this assert with the recorded per-state residual table instead.*)
DictReport::notuniversal = "Dictionary not universal: max residual=`1` (see validation table).";
allResidNums = Flatten[Table[N[validation[[i, 2]][[j, 2]]], {i, Length[validation]}, {j, Length[validation[[i, 2]]]}]];
maxResid = Max[Abs[allResidNums]];
If[maxResid > 10^-12, Message[DictReport::notuniversal, maxResid]; Abort[]];
dictReport = <|
  "generator {a,b}" -> genNontrivial,
  "φ(k)" -> (φuniversal[k] // Simplify),
  "scalar convention" -> "cScale[λ1,λ2][M,n] (see cell; 1 unless Step 3 set per-state)",
  "max |residual| over rep sweep" -> maxResid,
  "digits residual vanished" -> N[-Log10[maxResid + 10^-300]],
  "reps validated" -> repSweep
|>;
dictReport
```

- [ ] **Step 7: Run it; confirm the report**

Expected (success path): no `DictReport::notuniversal` message; `dictReport` shows the generator, a clean closed-form `φ(k)` (likely involving `(k κ1)^(θ1/h)`), max residual ~10^-15 or smaller, and the validated rep list. On the obstruction path: the recorded per-state residual table and best-fit φ, clearly labeled as a non-closure finding.

- [ ] **Step 8: Save the notebook.**

- [ ] **Step 9: Report the dictionary to the orchestrator** — the generator, `φ(k)`, the scalar (exact vs per-state), and whether universality held. This is the headline result of the whole task.

---

### ⛔ CHECKPOINT — end of Phase 2

The orchestrator commits Phase 2 (Section F notebook changes) and reports the discovered dictionary. Consider whether the result warrants: (a) folding the two spec refinements back into the design doc, (b) a `paper/xxxCG.tex` write-up (currently single-site only), (c) a HANDOFF.md update.

---

## Self-review notes

- **Spec coverage:** Part 1 build (t1/qdet1/xL1) → Task 1; Part 1 internal B-check → Task 1 Step 6; Part 1 factorization + predicted scalar → Task 2. Part 2 measure/bracket + Nn1 → Task 3; LHS FL → Task 4 Steps 1–2; RHS OP + generator choice → Task 4 Steps 3–8; discovery + universality + report → Task 5. Kernel-state risk → Task 0. All spec sections mapped.
- **Refinements** (LHS base `k·κ1`; `Nn1` check is `w_0==1` not `br1[1]==1`) are stated up front and applied in Tasks 3–4.
- **Type/shape consistency:** `xSingle`/`xL1` are 1×n rows; `Flatten` before `.Psi`; `Psi` length D; `OP`/`OPg` polynomials in φ; `FL` a function of k. Signatures `FL[λ1,λ2][M,n][k]`, `OP[λ1,λ2][a,b][M,n][φ]`, `OPg[λ1,λ2][M,n][φ]`, `gRot[λ][a,b][φ]`, `t1[λ][a,s][u]`, `xL1[λ][{s}]`, `mu1[λ][u]`, `br1[λ][f]`, `Nn1[λ]` used consistently across tasks.
- **Exploratory honesty:** Task 5 is the only non-deterministic task; it hard-asserts only internal universality and otherwise reports findings, matching the spec's "report the dictionary, don't assert the correspondence."

---

## Execution outcome & corrections (2026-07-16)

**Phase 1 (Section E) — done, committed `8dbed52`.** Factorization confirmed bit-exactly. One
correction to the plan's expectation: the site-2 passthrough is **not** the identity (μ≠1). It is a
`Qθ`-type eigenvalue `μSite2[λ2,s1] = Product[(θ1-θ2+(k-1)h-λ2 h),{k,1,s1}]`, giving
`c = Product[(θ1-θ2+(k-1)h-λ2 h)/(θ1-θ2+(k-1)h-λ1 h),{k,1,s1}]`. Task 2's `cPred` was corrected to
this and the μEff soft-warn became a hard assert; a new cell documents `μSite2`. (Cells 53–60.)

**Task 4 gRot bug (found & fixed during execution):** `Sum[φ^m/m! MatrixPower[Ee,m],{m,0,λ}]` leaves
the `m=0` term unevaluated — `MatrixPower[nilpotent,0]` fires `MatrixPower::sing`. Fixed by splitting
off the identity: `IdentityMatrix[dim] + Sum[..., {m,1,λ}]`. Do not use `MatrixExp` here.

**Task 5 — the plan's object was wrong; the discovery only closed after the fix.** Task 4/5 built the
**bare bracket** `FL = br1[k^(u/h) Q1]`, but the correct object is the **1×1 FSoV determinant**
(Section C's `det` divides by `f` at each node):

```
det1[λ1][f]  := br1[λ1][f] / (f /. u -> θ1)
FLdet[λ1,λ2][M,n][k] := det1[λ1][ k^(u/h) Q1[λ1,λ2][M,n][u] ]
```

The `1/f(θ1) = 1/(k^(θ1/h) Q1(θ1))` normalization cancels the `k^(θ1/h)` prefactor and the
state-dependent `Q1(θ1)`, leaving `FLdet` a degree-λ1 polynomial in k with constant term 1. The
bare-bracket `FL` (no normalization) produced a k-dependent per-state scalar and **did not close** —
that was a plan bug (spec said "the 1×1 determinant collapses to a single bracket", dropping the
normalization), not an executor error.

**Result (cells 72–75, cell 71 annotated):** the dictionary closes **exactly** —
```
det1[λ1][k^(u/h) Q1[λ1,λ2][M,n]] = (⟨x0[λ1]| exp[-k E21] ⊗ ⟨x0[λ2]|) . Psi[λ1,λ2][M,n]
```
generator **E21**, **φ(k) = -k** (α=-1, universal), residual literal-symbolic-`0` across the rep
sweep {(1,1),(2,1),(1,2),(2,2),(3,1),(1,3)}, 41 states. This is the CG-overlap correspondence the
task targeted.
