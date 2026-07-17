# su(3) general gl(3) reps per site — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generalize `Paul/Mathematica/XXX/Experiments/su3_V2_general.wb` from symmetric reps `{S,0,0}` per site to arbitrary gl(3) Young-diagram reps `λ={λ1,λ2,λ3}` (`λ1≥λ2≥λ3≥0`, integers), one per site, preserving V1's diagonal-twist frame, `(H1,H2,n)` labelling, Baxter Q-functions, and hard-assert verification.

**Architecture:** The notebook is edited in place via the **wolfbook-builder** subagent using the wolfbook MCP tools (`wolfbook_editCell`, `wolfbook_runCell`, `wolfbook_getCellOutput`). Most cells change only their signature (`S1,S2 → λ1,λ2`); a small number carry genuine math edits (`ν`, `wtMult`, the degree closed form). Each task edits the relevant defining cell(s), edits its verification cell to run on **non-symmetric** reps, runs it, and confirms zero-residual / expected output before committing.

**Tech Stack:** Wolfram Language (Mathematica) `.wb` notebook; wolfbook MCP tools; git.

## Global Constraints

- **Notebook:** `Paul/Mathematica/XXX/Experiments/su3_V2_general.wb` (40 cells, 0-indexed here; the executor locates cells by their leading `(* … *)` comment / defined symbol, not by a hard index).
- **Interface:** every function takes full 3-vectors `λ1,λ2` per site. No scalar-`S` wrappers.
- **The universal relabel rule** (applied to every "relabel" cell): textual `S1 → λ1`, `S2 → λ2` (signature `S1_→λ1_`, `S2_→λ2_`, and every call `[S1, S2]→[λ1, λ2]`), plus `dim[{S1,0,0}] → dim[λ1]`, `dim[{S2,0,0}] → dim[λ2]`, `{S1, 0, 0} → λ1`, `{S2, 0, 0} → λ2`, and any `S1+S2 → (Total[λ1]+Total[λ2])`. Do **not** otherwise change the body.
- **Yangian weight (the one physics edit), verified from the Lax action on the HWS:** `ν[λ1,λ2][i][u] = (u−θ1−h·λ1[[i]])(u−θ2−h·λ2[[i]])` for `i=1,2,3`.
- **Total boxes:** `N = Total[λ1]+Total[λ2]` (replaces V1's `S1+S2`).
- **Weight↔Cartan inversion (unchanged):** `n3=(N−H1−2H2)/3, n2=n3+H2, n1=n3+H1+H2`; `H1=n1−n2, H2=n2−n3`.
- **Numeric params unchanged:** `θ1=1/3, θ2=1/7, h=1`; twist `z[k]` high-precision, `χ1,χ2,χ3` their symmetric functions. Twist matrix `G` and cell 6 params are rep-independent — **do not touch**.
- **Validation sweep (the six site-pairs `(λ1,λ2)`):**
  ```
  sweep = {{{1,0,0},{1,0,0}}, {{2,0,0},{1,1,0}}, {{1,1,0},{1,1,0}},
           {{2,1,0},{1,0,0}}, {{2,1,0},{1,1,0}}, {{2,2,0},{1,0,0}}};
  ```
- **Cells that stay literally unchanged:** 0 (`ps/dim/J0/Jp/Jm/JJ`), 5 (`G`), 6 (numeric params), 14 (param verify), 35 (`qpoly`).
- **Verification convention:** each check `Abort`s on a real failure; a passing run returns `{0}` / `"…OK"` / `worst→0`. After a passing run, **save the notebook** (`wolfbook_fileOps` save) then commit.
- **Kernel hygiene:** these are memoized defs (`f[x_]:=f[x]=…`). After editing any defining cell, re-run it **and** every downstream cell that may have cached stale values before trusting a check (per `CLAUDE.md`).
- **Commit trailer:** end each commit message with `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`. Work on `main` (repo convention for this project).

---

### Task 1: Section A — embeddings, Lax, monodromy, and the `ν` generalization

**Files:**
- Modify: `su3_V2_general.wb` cell 1 (`Ee/id/idd/EE`), cell 2 (`C1`), cell 4 (`L`), cell 7 (`T/Tg/TG`), cell 8 (`ν/HWS`).
- Verify: cell 3 (`C1` check), cell 9 (HWS check).

**Interfaces:**
- Consumes: `ps/dim/JJ` (cell 0, unchanged); `θ1,θ2,h,z[k],χ` (cells 5–6, unchanged).
- Produces: `Ee[λ][i,j]`, `id[λ]`, `idd[λ1,λ2]`, `EE[λ1,λ2][i,j]`, `C1[λ1,λ2]`, `L[λ][i,j][u]`, `T/Tg/TG[λ1,λ2][i,j][u]`, `ν[λ1,λ2][i][u]`, `HWS[λ1,λ2]`.

- [ ] **Step 1: Rewrite cell 1 (embeddings)** to:

```wolfram
ClearAll[Ee,id,EE,idd];
Ee[λ_][i_,j_]:=JJ[i,j,λ];
id[λ_]:=IdentityMatrix[dim[λ]];
idd[λ1_,λ2_]:=KroneckerProduct[id[λ1],id[λ2]];
EE[λ1_,λ2_][i_,j_]:=KroneckerProduct[Ee[λ1][i,j],id[λ2]]+KroneckerProduct[id[λ1],Ee[λ2][i,j]];
```

- [ ] **Step 2: Rewrite cell 2 (central charge)** to:

```wolfram
(*central charge*)
ClearAll[C1];
C1[λ1_,λ2_]:=Sum[EE[λ1,λ2][i,i],{i,3}];
```

- [ ] **Step 3: Rewrite cell 4 (Lax)** to:

```wolfram
(*Lax operator for a general gl(3) rep λ={λ1,λ2,λ3} at the site*)
ClearAll[L];
L[λ_][i_,j_][u_]:=IdentityMatrix[dim[λ]] u KroneckerDelta[i,j] - h JJ[j,i,λ];
```

- [ ] **Step 4: Rewrite cell 7 (monodromy)** to:

```wolfram
(*monodromy with no twist, diagonal twist, and companion twist*)
ClearAll[T,Tg,TG];
T[λ1_,λ2_][i_,j_][u_]:=Sum[KroneckerProduct[L[λ1][i,ii][u-θ1],L[λ2][ii,j][u-θ2]],{ii,3}];
Tg[λ1_,λ2_][i_,j_][u_]:= z[i] T[λ1,λ2][i,j][u];
TG[λ1_,λ2_][i_,j_][u_]:=Sum[T[λ1,λ2][i,ii][u] G[[ii,j]],{ii,3}];
```

- [ ] **Step 5: Rewrite cell 8 (`ν` — the physics edit) and `HWS`)** to:

```wolfram
(*Yangian weights - eigenvalues of T[i,i] on the highest-weight state.
  General rep: T[i,i]|hw> = (u-θ1-h λ1[[i]])(u-θ2-h λ2[[i]]) |hw>, from L[λ][i,i]|hw>=(u-h λ_i)|hw>.
  Symmetric case λ=(S,0,0) reduces to ν[1]=(u-θ1-hS1)(u-θ2-hS2), ν[2]=ν[3]=(u-θ1)(u-θ2).*)
ClearAll[ν,HWS];
ν[λ1_,λ2_][i_][u_]:=(u-θ1-h λ1[[i]])(u-θ2-h λ2[[i]]);
HWS[λ1_,λ2_]:=Table[{KroneckerDelta[1,ii]},{ii,1,dim[λ1] dim[λ2]}];
```

- [ ] **Step 6: Rewrite the cell-3 `C1` check** (total boxes now) to:

```wolfram
(*checking central charge constraint = total boxes, on general reps*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   C1[r[[1]],r[[2]]]-(Total[r[[1]]]+Total[r[[2]]]) idd[r[[1]],r[[2]]],{r,reps}]]
```

- [ ] **Step 7: Rewrite the cell-9 HWS check** (general reps) to:

```wolfram
(*checking highest-weights against the general ν, on non-symmetric reps*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   T[r[[1]],r[[2]]][i,i][u].HWS[r[[1]],r[[2]]]-ν[r[[1]],r[[2]]][i][u] HWS[r[[1]],r[[2]]],
   {r,reps},{i,3}]]
```

- [ ] **Step 8: Run cells 1,2,4,7,8 (defs), then the cell-3 and cell-9 checks.**
Run each with `wolfbook_runCell`; read output with `wolfbook_getCellOutput`.
Expected: cell 3 → `{0}`; cell 9 → `{0}`. (Any non-`{0}` output means the `ν` formula or a relabel is wrong — stop and diagnose, do not proceed.)

- [ ] **Step 9: Save and commit.**

```bash
git add Paul/Mathematica/XXX/Experiments/su3_V2_general.wb
git commit -m "su(3) general reps: Section A embeddings/Lax/monodromy + general ν (Task 1)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 2: Section A — quantum minors, transfer matrices, qdetT (ν-ordering guard)

**Files:**
- Modify: cell 10 (`qmin/qming/qminG`), cell 11 (`tg/tG`), cell 12 (`qdetT`), cell 15 (`ν1`), cell 16 (`tasympt[1,1]`, exploratory), cell 17 (`Bax` reference form).
- Verify: cell 13 (`qdetT == tg[·][3,1]`).

**Interfaces:**
- Consumes: `T/Tg/TG`, `ν`, `idd` (Task 1).
- Produces: `qmin/qming/qminG[λ1,λ2][II,JJ][u]`, `tg/tG[λ1,λ2][a,1][u]`, `qdetT[λ1,λ2][u]`, `ν1[λ1,λ2][u]`, `Bax[λ1,λ2]`.

- [ ] **Step 1: Relabel cells 10, 11** using the universal relabel rule (`S1→λ1`, `S2→λ2` throughout; these cells contain only signature + call-site occurrences, no `{S,0,0}`).

- [ ] **Step 2: Rewrite cell 12 (`qdetT`)** to:

```wolfram
(*value of quantum determinant (χ3 carried here; general ν)*)
ClearAll[qdetT];
qdetT[λ1_,λ2_][u_]:=χ3 ν[λ1,λ2][1][u-2h] ν[λ1,λ2][2][u-h] ν[λ1,λ2][3][u];
```

- [ ] **Step 3: Rewrite cell 15 (`ν1`)** to:

```wolfram
(*single source of truth for the Yangian weight of T11*)
ClearAll[ν1];
ν1[λ1_,λ2_][u_]:=ν[λ1,λ2][1][u];
```

- [ ] **Step 4: Relabel cell 17 (`Bax`)** using the universal rule (`S1→λ1`, `S2→λ2`; body structure unchanged — it already reads `tg`, `ν1`, `qdetT`).

- [ ] **Step 5: Relabel cell 16 (`tasympt[1,1]`)** using the universal rule. This cell is an exploratory asymptotic expression **not consumed downstream** (`alphaCoeffs` uses `SeriesCoefficient` directly). If it references a symbol that no longer resolves, leave a `(* NB: exploratory, not used by alphaCoeffs *)` note and move on — it is **non-blocking** for this task.

- [ ] **Step 6: Rewrite the cell-13 `qdetT` check** (the ν-ordering guard, general reps) to:

```wolfram
(*qdetT central value == tg[·][3,1]; on non-symmetric reps this is the ν shift-ordering guard*)
Module[{reps={{{2,1,0},{1,0,0}},{{1,1,0},{1,1,0}},{{2,2,0},{1,0,0}},{{2,1,0},{1,1,0}}}},
 Union@Flatten@Table[
   Simplify[qdetT[r[[1]],r[[2]]][u] idd[r[[1]],r[[2]]]-tg[r[[1]],r[[2]]][3,1][u]],{r,reps}]]
```

- [ ] **Step 7: Run cells 10,11,12,15,17 (defs), then cell 13.**
Expected: cell 13 → `{0}`.
**If it is not `{0}`:** the shift assignment in `qdetT` is wrong for non-symmetric reps. Fix by permuting which `ν[i]` receives which `u−(k)h` shift (try the reverse order `ν[3][u−2h] ν[2][u−h] ν[1][u]` or the ascending-`i` order), re-run until `{0}`. Do **not** change the `ν` definition itself (Task 1 already verified it via HWS).

- [ ] **Step 8: Save and commit.**

```bash
git add Paul/Mathematica/XXX/Experiments/su3_V2_general.wb
git commit -m "su(3) general reps: quantum minors + transfer matrices + qdetT (Task 2)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 3: `wtMult` — honest tensor-product weight tally

**Files:**
- Modify: cell 21 (`wtMult`).
- Verify: cell 22 (sum-rule), cell 23 (second-rep guard).

**Interfaces:**
- Consumes: `ps/dim/J0` (cell 0), `Total[λ]` for `N`.
- Produces: `siteWeightMult[λ]` (Association `{w1,w2,w3}→multiplicity`), `wtMult[λ1,λ2][H1,H2]` (Integer ≥ 0).

- [ ] **Step 1: Rewrite cell 21 (`wtMult`)** to the general tensor-product tally:

```wolfram
(*gl(3) weight multiplicities in the tensor product of two general reps λ1⊗λ2.
  siteWeightMult[λ]: tally the gl(3) weight (J0[1],J0[2],J0[3]) of every GT pattern of λ.
  wtMult: recover (n1,n2,n3) from (H1,H2) with N=total boxes, then convolve the two site tallies.*)
ClearAll[siteWeightMult,wtMult];
siteWeightMult[λ_]:=siteWeightMult[λ]=Counts[({J0[1][#],J0[2][#],J0[3][#]}&)/@ps[{λ}]];
wtMult[λ1_,λ2_][H1_,H2_]:=Module[{Ntot=Total[λ1]+Total[λ2],n1,n2,n3,m1,m2},
  n3=(Ntot-H1-2 H2)/3; n2=n3+H2; n1=n3+H1+H2;
  If[!(IntegerQ[n1]&&IntegerQ[n2]&&IntegerQ[n3]),Return[0]];
  m1=siteWeightMult[λ1]; m2=siteWeightMult[λ2];
  Total[KeyValueMap[Function[{w,c}, c Lookup[m2,{n1,n2,n3}-w,0]],m1]]
];
```

- [ ] **Step 2: Rewrite cell 22 (sum-rule check)** to run over several general reps:

```wolfram
(*hard-assert sum-rule: Σ over all (H1,H2) of wtMult == dim[λ1] dim[λ2], for every test rep*)
Module[{reps={{{2,1,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},bad},
  bad=Select[reps,Function[r,
     Module[{Ntot=Total[r[[1]]]+Total[r[[2]]],d,tot},
       d=dim[r[[1]]] dim[r[[2]]];
       tot=Total@Flatten@Table[wtMult[r[[1]],r[[2]]][H1,H2],
           {H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}];
       tot=!=d]]];
  If[bad=!={},Print["SUM-RULE FAIL: ",bad];Abort[]];
  {"all reps sum-rule OK",
   "adjoint⊗fund @(H1,H2)=(0,0) mult"->wtMult[{2,1,0},{1,0,0}][0,0]}
]
```

- [ ] **Step 3: Rewrite cell 23 (second-rep guard)** to a single higher-dim non-symmetric spot check:

```wolfram
(*guard on {2,1,0}⊗{2,1,0} (d=64) — a rep with genuine interior multiplicity on both sites*)
Module[{λ1={2,1,0},λ2={2,1,0},Ntot,d,tot},
  Ntot=Total[λ1]+Total[λ2]; d=dim[λ1] dim[λ2];
  tot=Total@Flatten@Table[wtMult[λ1,λ2][H1,H2],{H1,-3 Ntot,3 Ntot},{H2,-3 Ntot,3 Ntot}];
  If[tot!=d,Print["SUM-RULE FAIL: ",tot," != ",d];Abort[]];
  {"total weights"->tot,"dim"->d}
]
```

- [ ] **Step 4: Run cell 21 (def), then cells 22 and 23.**
Expected: cell 22 → `{"all reps sum-rule OK", "adjoint⊗fund @(H1,H2)=(0,0) mult"->k}` with `k ≥ 1`; cell 23 → `{"total weights"->64, "dim"->64}`.

- [ ] **Step 5: Save and commit.**

```bash
git add Paul/Mathematica/XXX/Experiments/su3_V2_general.wb
git commit -m "su(3) general reps: wtMult tensor-product weight tally (Task 3)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 4: Eigensystem — `tgCoeffs`, `H1op/H2op`, `TauEigensystem3`, `Psi/tau1/tau2`

**Files:**
- Modify: cell 18 (`tgCoeffs`), cell 19 (`H1op/H2op`), cell 24 (`TauEigensystem3`), cell 25 (`Psi/tau1/tau2`).
- Verify: cell 20 (leading-coeff + commuting-family), cell 26 (eigensystem residual spot check).

**Interfaces:**
- Consumes: `tg`, `EE`, `χ1,χ2`, `wtMult` (Tasks 1–3).
- Produces: `tgCoeffs[λ1,λ2]`, `H1op/H2op[λ1,λ2]`, `TauEigensystem3[λ1,λ2]` (Association `{H1,H2,n}→<|"Psi","tau1c","tau2c"|>`, with internal sector-count hard-assert against `wtMult`), `Psi/tau1/tau2[λ1,λ2][H1,H2,n]`.

- [ ] **Step 1: Relabel cell 18 (`tgCoeffs`) and cell 19 (`H1op/H2op`)** using the universal rule (`S1→λ1`, `S2→λ2`).

- [ ] **Step 2: Relabel cell 24 (`TauEigensystem3`)** using the universal rule. In this cell specifically apply: `S1_→λ1_`, `S2_→λ2_`, all `[S1, S2]→[λ1, λ2]` (9 call sites), `{S1, 0, 0}→λ1`, `{S2, 0, 0}→λ2` (so `d=dim[λ1] dim[λ2]`). Leave the `::rankdef/::notgenuine/::noninteger/::degnull` message *text* as-is (cosmetic `S1,S2` in strings are fine), but their `Message[…, S1, S2]` arguments become `…, λ1, λ2`.

- [ ] **Step 3: Relabel cell 25 (`Psi/tau1/tau2`)** using the universal rule (3× signature, 3× `[S1, S2]` calls).

- [ ] **Step 4: Rewrite cell 20 (leading-coeff + commuting check)** for a non-symmetric rep:

```wolfram
(*hard-assert: τ1,τ2 leading coeffs are χ1 Id, χ2 Id; full commuting family pairwise-commutes.
  Tested on the adjoint⊗fundamental (d=24), which has an interior multiplicity-2 weight.*)
Module[{λ1={2,1,0},λ2={1,0,0},d,c1,c2,C2,D4,mats,comm},
  d=dim[λ1] dim[λ2];
  {c1,c2}=tgCoeffs[λ1,λ2];
  C2=c1[[3]]; D4=c2[[5]];
  If[Max[Abs[Flatten[C2-χ1 IdentityMatrix[d]]]]>10^-13,Print["C2 FAIL"];Abort[]];
  If[Max[Abs[Flatten[D4-χ2 IdentityMatrix[d]]]]>10^-13,Print["D4 FAIL"];Abort[]];
  mats=Join[c1,c2,{H1op[λ1,λ2],H2op[λ1,λ2]}];
  comm=Max[Table[Max[Abs[Flatten[mats[[i]].mats[[j]]-mats[[j]].mats[[i]]]]],
      {i,Length[mats]},{j,Length[mats]}]];
  If[comm>10^-12,Print["COMMUTATION FAIL: ",comm];Abort[]];
  {"C2=χ1 Id OK","D4=χ2 Id OK","maxCommutator"->comm}
]
```

- [ ] **Step 5: Relabel cell 26 (eigensystem spot check)** using the universal rule, and change its target rep to the multiplicity-`>1` adjoint pair: set `λ1={2,1,0}, λ2={1,0,0}` (replacing the `S1=2,S2=1` locals). Keep its residual-extraction body intact — it reports the worst common-eigenvector residual. (The per-sector count `== wtMult` assertion lives inside `TauEigensystem3` itself, so running it on this rep exercises the multiplicity-`>1` splitting.)

- [ ] **Step 6: Run cells 18,19,24,25 (defs), then cells 20 and 26.**
Expected: cell 20 → `{"C2=χ1 Id OK","D4=χ2 Id OK","maxCommutator"->c}` with `c` at machine-zero (`< 10^-12`); cell 26 → worst residual at machine-zero and no `::degnull`/`::notgenuine`/`::rankdef` messages.
**If `TauEigensystem3` throws `::rankdef` or `::degnull` on `{2,1,0}⊗{1,0,0}`:** a multiplicity-`>1` sector failed to split. First let its internal r-vector retry run; if it still fails, report the sector `(H1,H2)` and its `wtMult` — this is a genuine finding (two states sharing all conserved charges), not something to silence. Do not weaken the assert.

- [ ] **Step 7: Save and commit.**

```bash
git add Paul/Mathematica/XXX/Experiments/su3_V2_general.wb
git commit -m "su(3) general reps: TauEigensystem3 + Psi/tau accessors (Task 4)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 5: Baxter solver + degree closed-form extraction

**Files:**
- Modify: cell 27 (`BaxOp`), cell 28 (`alphaCoeffs`), cell 29 (`Mdeg`), cell 36 (`QSolve3`), cell 37 (`Qfun`).
- Modify/replace: cells 30–33 (degree tabulations), cell 34 (degree closed-form hard-assert).
- Verify: cell 38 (single-state TQ residual), and the new closed-form assert (cell 34).

**Interfaces:**
- Consumes: `tau1/tau2`, `ν1`, `qdetT`, `TauEigensystem3`, `z[k]`, `χ2,χ3`, `qpoly` (cell 35, unchanged).
- Produces: `BaxOp[λ1,λ2][H1,H2,n][Qf][u]`, `alphaCoeffs[λ1,λ2][H1,H2,n]`, `Mdeg[λ1,λ2][H1,H2,n]` (`{M1,M2,M3}`), `QSolve3[λ1,λ2][H1,H2,n]` (`{q1,q2,q3}` monic coeff lists), `Qfun[λ1,λ2][H1,H2,n][k][u]`, and the **extracted closed form** `Mpred[λ1,λ2][H1,H2]`.

- [ ] **Step 1: Relabel cells 27, 28, 29, 36, 37** using the universal rule (`S1→λ1`, `S2→λ2`; these are call-site/signature only). `Mdeg` (cell 29) is rep-agnostic — its body needs no change beyond the signature.

- [ ] **Step 2: Run cells 27,28,29,36,37 (defs), then cell 38 (single-state TQ) after retargeting it** to `λ1={2,1,0}, λ2={1,0,0}, H1=?,H2=?,n=0`. First evaluate `Keys[TauEigensystem3[{2,1,0},{1,0,0}]]` (via `wolfbook_evaluateExpression`) to pick a valid `(H1,H2,n)` key, then set cell 38's locals to it:

```wolfram
(*single-state TQ-residual check on a non-symmetric rep: each Q_k annihilated by BaxOp*)
Module[{λ1={2,1,0},λ2={1,0,0},H1=<pick>,H2=<pick>,n=0,res},
  res=Table[
     Max@Table[
        Abs@Chop[N[BaxOp[λ1,λ2][H1,H2,n][Function[uu,Qfun[λ1,λ2][H1,H2,n][k][uu]]][u0],30]],
        {u0,{2/3,9/5,4,23/7}}],
     {k,1,3}];
  If[Max[res]>10^-10,Print["TQ RESIDUAL FAIL: ",res];Abort[]];
  {"perQresidual"->res}
]
```

Expected: `{"perQresidual"->{~0,~0,~0}}`. This confirms `QSolve3`/`Mdeg`/`Qfun` work per-state on a general rep.

- [ ] **Step 3: Replace cell 30 with a full-sweep degree-tabulation** that collects `(λ1,λ2,H1,H2,n)→{M1,M2,M3}` as raw data for the fit:

```wolfram
(*collect Mdeg over the whole sweep as {λ1,λ2,H1,H2,n,{M1,M2,M3}} rows for closed-form fitting*)
degData=Flatten[Table[
   Module[{λ1=rep[[1]],λ2=rep[[2]],es},
     es=TauEigensystem3[λ1,λ2];
     KeyValueMap[Function[{key,rec},
        {λ1,λ2,key[[1]],key[[2]],key[[3]],Mdeg[λ1,λ2][key[[1]],key[[2]],key[[3]]]}],es]],
   {rep,{{{1,0,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},
         {{2,1,0},{1,0,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}}}],1];
Column[degData]
```

Run it. **Sanity gate:** every `{M1,M2,M3}` entry must be a triple of non-negative integers (Mdeg already asserts near-integer). Note whether `M_k` depends on `n` at fixed `(λ1,λ2,H1,H2)` — expected **not** to (as in V1).

- [ ] **Step 4: Fit the closed form.** Replace cell 31 with a candidate-testing cell. For each row, recover `(n1,n2,n3)` from `(H1,H2,N)` (`N=Total[λ1]+Total[λ2]`) and test candidate predictors against the collected `{M1,M2,M3}`:

```wolfram
(*test candidate closed forms for M_k against degData. n_k from (H1,H2) with N=total boxes.
  Candidates (each must give M_k in index k, matching z[k]/n_k pairing):
    cand["N-n"]      : M_k = N - n_k
    cand["Lam1-n"]   : M_k = Λ_1 - n_k,   Λ = λ1 + λ2  (pseudovacuum highest weight)
    cand["dualpart"] : M_k = Σ_{j>=k}(Λ_j) ... (partial-sum family, add as needed)
  Report which candidate has zero max-deviation over ALL rows.*)
nvec[H1_,H2_,Ntot_]:=Module[{n3=(Ntot-H1-2 H2)/3},{n3+H1+H2,n3+H2,n3}];
cand["N-n"][λ1_,λ2_,H1_,H2_]:=(Total[λ1]+Total[λ2])-nvec[H1,H2,Total[λ1]+Total[λ2]];
cand["Lam1-n"][λ1_,λ2_,H1_,H2_]:=(λ1+λ2)[[1]]-nvec[H1,H2,Total[λ1]+Total[λ2]];
results=Association@Table[
   name->Max[Abs[Flatten[
      (cand[name][#[[1]],#[[2]],#[[3]],#[[4]]]-#[[6]])&/@degData]]],
   {name,{"N-n","Lam1-n"}}];
results
```

Run it. **The winning candidate has value `0`.** If none is `0`, add further candidates to `cand[...]` (partial-sum families in the `λ` components — e.g. `Σ_{i}λ_i`-based, or `Λ_k`-shifted forms) until one gives exactly `0` across all rows. Record the winning formula in the cell's comment. (This is the "figure it out" deliverable; the numeric `Mdeg` degrees are already trusted, so this step cannot block verification — only the closed form is being pinned down.)

- [ ] **Step 5: Clear cells 32, 33** (the old symmetric-rep per-rep tabulations, now subsumed by Steps 3–4). Replace each with an empty comment cell `(* superseded by degData tabulation + closed-form fit above *)` or delete via `wolfbook_deleteCell`.

- [ ] **Step 6: Rewrite cell 34 as the closed-form hard-assert** using the winning candidate from Step 4. Template (substitute `<WINNER>` = the winning `cand[...]` body):

```wolfram
(*hard-assert: extracted closed form for M_k matches Mdeg exactly (zero residual) over the sweep*)
Mpred[λ1_,λ2_][H1_,H2_]:=<WINNER>;   (* e.g. (Total[λ1]+Total[λ2]) - nvec[H1,H2,Total[λ1]+Total[λ2]] *)
Module[{sweep={{{1,0,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},
               {{2,1,0},{1,0,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}},worst},
  worst=Max@Table[
     Module[{λ1=rep[[1]],λ2=rep[[2]],es},
       es=TauEigensystem3[λ1,λ2];
       Max@KeyValueMap[Function[{key,rec},
          Max[Abs[Mdeg[λ1,λ2][key[[1]],key[[2]],key[[3]]]-Mpred[λ1,λ2][key[[1]],key[[2]]]]]],es]],
     {rep,sweep}];
  If[worst>0,Print["DEGREE-FORMULA MISMATCH: ",worst];Abort[]];
  {"closed-form M_k matches Mdeg over full sweep, worst"->worst}
]
```

Run it. Expected: `{"…worst"->0}`.

- [ ] **Step 7: Save and commit.**

```bash
git add Paul/Mathematica/XXX/Experiments/su3_V2_general.wb
git commit -m "su(3) general reps: Baxter solver + degree closed form (Task 5)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 6: Full-sweep verification (acceptance)

**Files:**
- Modify: cell 39 (full-sweep TQ + Casoratian verification).

**Interfaces:**
- Consumes: everything from Tasks 1–5.
- Produces: the acceptance report (per-rep `worstTQ`, `worstCasoratian`).

- [ ] **Step 1: Rewrite cell 39's `sweep`** to the six general-rep pairs (the rest of the cell body — TQ residual per `Q_k`, 4×3 Casoratian reconstruction of `A1,A2,A3` — is rep-agnostic and stays; apply the universal relabel `S1→λ1, S2→λ2` inside it):

```wolfram
Module[{sweep={{{1,0,0},{1,0,0}},{{2,0,0},{1,1,0}},{{1,1,0},{1,1,0}},
               {{2,1,0},{1,0,0}},{{2,1,0},{1,1,0}},{{2,2,0},{1,0,0}}}, report},
  report=Table[
    Module[{λ1=rep[[1]],λ2=rep[[2]],es,worstTQ=0,worstCas=0},
      es=TauEigensystem3[λ1,λ2];
      KeyValueMap[
       Function[{key,rec},
        Module[{H1=key[[1]],H2=key[[2]],n=key[[3]],Qs,u0,rows,minors,
          A1rec,A2rec,A3rec,A1,A2,A3,N1,N1m,N1mm},
         Do[worstTQ=Max[worstTQ,
             Abs@Chop[N[BaxOp[λ1,λ2][H1,H2,n][
                 Function[uu,Qfun[λ1,λ2][H1,H2,n][k][uu]]][u0],30]]],
           {k,1,3},{u0,{5/4,11/3}}];
         u0=13/5;
         Qs[uu_]:=Table[Qfun[λ1,λ2][H1,H2,n][k][uu],{k,1,3}];
         rows={Qs[u0+h],Qs[u0],Qs[u0-h],Qs[u0-2 h]};
         minors=Table[Det[Delete[rows,i]],{i,1,4}];
         A1rec=minors[[2]]/minors[[1]]; A2rec=minors[[3]]/minors[[1]]; A3rec=minors[[4]]/minors[[1]];
         N1=ν1[λ1,λ2][u0]; N1m=ν1[λ1,λ2][u0-h]; N1mm=ν1[λ1,λ2][u0-2 h];
         A1=tau1[λ1,λ2][H1,H2,n][u0]/N1;
         A2=tau2[λ1,λ2][H1,H2,n][u0]/(N1 N1m);
         A3=qdetT[λ1,λ2][u0]/(N1 N1m N1mm);
         worstCas=Max[worstCas,
            Abs@Chop[N[A1rec-A1,30]],Abs@Chop[N[A2rec-A2,30]],Abs@Chop[N[A3rec-A3,30]]]]],
       es];
      {rep,"nstates"->Length[es],"worstTQ"->worstTQ,"worstCasoratian"->worstCas}],
    {rep,sweep}];
  If[Max[#[[3,2]]&/@report,#[[4,2]]&/@report]>10^-9,
    Print["SWEEP FAIL: ",report];Abort[]];
  Column[report]
]
```

- [ ] **Step 2: Fresh-kernel top-to-bottom run.** Restart the kernel (`wolfbook_kernelControl`), then run all cells 0→39 in order (defs before checks; re-run memoized defs so no stale cache survives). Confirm every check cell (3, 9, 13, 20, 22, 23, 26, 34, 38, 39) returns its zero/OK result and nothing `Abort`s.

- [ ] **Step 3: Read cell 39 output.** Expected: six rows, each `{rep, "nstates"->d, "worstTQ"->~0, "worstCasoratian"->~0}` with both worsts `< 10^-9` (ideally machine-zero). Record the total state count.

- [ ] **Step 4: Save and commit.**

```bash
git add Paul/Mathematica/XXX/Experiments/su3_V2_general.wb
git commit -m "su(3) general reps: full-sweep verification, all zero-residual (Task 6)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

- [ ] **Step 5: Update `HANDOFF.md`** — add a dated state block recording: `su3_V2_general.wb` generalized to arbitrary gl(3) reps, the general `ν` formula, the new `wtMult`, the **extracted degree closed form** (the winner from Task 5), and the full-sweep result (reps, total states, `worstTQ=worstCasoratian=0`). Add an `## Execution outcome` section to this plan's design doc `docs/superpowers/specs/2026-07-17-su3-general-reps-design.md` (mirroring V1's convention). Commit both.

---

## Self-review

**Spec coverage:**
- λ-vectors everywhere → Tasks 1–6 (universal relabel rule). ✓
- General `ν` (physics edit) → Task 1 Step 5, guarded by HWS check (T1) and qdetT check (T2). ✓
- `wtMult` honest tally → Task 3. ✓
- Degree closed-form re-derivation → Task 5 Steps 3–6. ✓
- Mixed-small-dims sweep → Global Constraints + Tasks 3–6. ✓
- Verification suite (TQ, Casoratian, eigenvalue-consistency, sector counts, degree) → Task 4 (sector counts via internal assert), Task 5 (TQ single-state + degree), Task 6 (TQ + Casoratian full sweep). ✓
- Section-A guard checks (C1, qdetT==tg[3,1], leading coeffs) → T1, T2, T4. ✓
- Risk 1 (ν ordering) → Task 2 Step 7 fallback. Risk 2 (sector splitting) → Task 4 Step 6 finding-not-silence. Risk 3 (no clean closed form) → Task 5 Step 4 note. ✓

**Placeholder scan:** the only intentional `<…>` placeholders are (a) the `(H1,H2,n)` key in Task 5 Step 2, resolved by evaluating `Keys[TauEigensystem3[…]]` first, and (b) `<WINNER>` in Task 5 Step 6, resolved by Task 5 Step 4's fit — both are data-dependent by design (the degree closed form is unknown until extracted), with an explicit resolution procedure. No other TBDs.

**Type consistency:** `siteWeightMult`→`wtMult` (T3), `tgCoeffs`→`TauEigensystem3` (T4), `Mdeg`/`Qfun`/`BaxOp`→sweep (T5→T6), `Mpred`/`nvec` defined in T5 and reused in the T5 assert — all names consistent across tasks. Every function uses the `[λ1,λ2]` signature uniformly.
