# su(3) transfer-matrix eigensystem + Baxter Q-functions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. All notebook work is done live via the **wolfbook-builder** subagent / wolfbook MCP tools against the running kernel for `su3_V1.wb`.

**Goal:** Build the gl(3) analog of the su(2) Clean notebook's Section B — an eigenvector builder + transfer-matrix eigenvalues (labeled by `(H1,H2,n)`) + the three Baxter Q-functions — inside `Paul/Mathematica/XXX/Experiments/su3_V1.wb`.

**Architecture:** Direct simultaneous diagonalization of the commuting family {coeffs of `τ1(u)=tg[·][1,1]`, coeffs of `τ2(u)=tg[·][2,1]`, Cartan `H1op,H2op`} in the diagonal-twist frame, numerically at fixed numeric parameters, keeping eigenvalues symbolic in `u`. The 3rd-order Baxter TQ equation is then solved as a linear (nullspace) problem per twist-exponent `z_k`, with degrees `M_k` predicted by asymptotic matching and verified.

**Tech Stack:** Wolfram Language (wolfbook `.wb` notebook), wolfbook MCP tools for live cell insert/run/save.

## Global Constraints

- **Notebook (edit in place):** `Paul/Mathematica/XXX/Experiments/su3_V1.wb`. Append new cells after the existing Section A; do not restructure existing cells except the two reconciliations named in Task 1.
- **Numeric parameters (Task 1, used everywhere):** `θ1 = 1/3`, `θ2 = 1/7`, `h = 1`; `z[1] = N[Exp[I Zeta[3]], 30]`, `z[2] = N[Exp[I Zeta[5]], 30]`, `z[3] = N[Exp[I Zeta[7]], 30]`. `χ1,χ2,χ3` remain the elementary symmetric functions of `z[i]` (already defined in Section A).
- **Twist frame = diagonal** (`tg`, `Tg`, `qming`). Never use the companion-twist `tG`/`qminG` in this work.
- **Diagonal-twist quantum determinant** carries `χ3` inside `qdetT`; the Baxter equation must NOT reintroduce an explicit `χ3` (avoid double-counting).
- **Reps swept in verification:** `{(1,1),(2,1),(1,2),(2,2),(3,1),(1,3)}` as `(S1,S2)`.
- **Notebook conventions (from CLAUDE.md):** multi-statement cells need `;` between statements; after editing a memoized (`f[x_]:=f[x]=…`) definition, re-run the defining cell and everything downstream that used the stale cache; cell edits/evaluations don't persist until the notebook is explicitly saved (`wolfbook_fileOps` save) — save at the end of every task before committing.
- **Check convention:** every routine is followed by a hard-assert check cell (`Abort[]`/`Message` on failure, as in su(2)); a passing run shows the designed success output. Numeric residual tolerance `10^-13` (30-digit arithmetic); integer-rounding tolerance `10^-10`.
- **Commit:** after each task, save the notebook then `git add Paul/Mathematica/XXX/Experiments/su3_V1.wb` and commit. End commit messages with the `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>` trailer.

---

## File Structure

Single file, `su3_V1.wb`. New cells grouped as "Section B" after the existing Section-A cells (the last existing cell is the sketched `Bax`). Logical blocks:

- **Params + reconciliation** (Task 1) — numeric parameter cell; fix `Bax` names; `ν1` alias.
- **Commuting family** (Task 2) — `tgCoeffs`, `H1op`, `H2op`, commutation check.
- **Weight multiplicity** (Task 3) — `wtMult`.
- **Eigensystem** (Task 4) — `TauEigensystem3`, accessors `Psi`, `tau1`, `tau2`.
- **Baxter + degrees** (Task 5) — `BaxOp`, `alphaCoeffs`, `Mdeg`.
- **Q-solver** (Task 6) — `QSolve3`, `Qfun` (accessor) via `qpoly`.
- **Verification** (Task 7) — sweep over reps: per-Q TQ residual + 4×4 Casoratian reconstruction of `A1,A2,A3`.

---

## Task 1: Numeric parameters + Baxter/name reconciliation

**Files:**
- Modify: `Paul/Mathematica/XXX/Experiments/su3_V1.wb` (add one param cell near top-of-Section-B; edit the existing `Bax` cell; add `ν1` alias)

**Interfaces:**
- Produces: global numeric values `θ1,θ2,h,z[1],z[2],z[3]` (hence numeric `χ1,χ2,χ3,G`); `ν1[S1_,S2_][u_]`; a corrected 3rd-order Baxter reference cell.

- [ ] **Step 1: Insert the numeric-parameter cell** (place it immediately after the companion-twist `G` cell so `χ1,χ2,χ3,G` pick up numeric `z`).

```wolfram
(*Numeric parameters (shared across Section B). θ,h copied from the su(2) Clean notebook;
  three DISTINCT unit-modulus twist eigenvalues z[1],z[2],z[3] (generic → non-degenerate;
  no ∏z=1 constraint since gl(3) is not unimodular-normalized). χ1,χ2,χ3 and G were defined
  in Section A as symmetric functions of z[i], so they become numeric automatically.*)
θ1 = 1/3;
θ2 = 1/7;
h = 1;
z[1] = N[Exp[I Zeta[3]], 30];
z[2] = N[Exp[I Zeta[5]], 30];
z[3] = N[Exp[I Zeta[7]], 30];
{θ1, θ2, h, z[1], z[2], z[3], χ1, χ2, χ3}
```

- [ ] **Step 2: Re-run downstream Section-A cells that cached symbolic values.** Re-run (in order) the `G`/`χ` cell, the monodromy cell (`T,Tg,TG`), `ν`/`HWS`, `qmin`/`qming`/`qminG`, `tg`/`tG`, and `qdetT`. (Any memoized `qming`/`qmin` values computed while `z` was symbolic must be cleared — safest is a kernel restart then top-to-bottom run of Section A, then this param cell.)

- [ ] **Step 3: Verify numerics took hold.** Insert-and-run:

```wolfram
{Precision[χ1], N[χ3], Chip = χ1 - (z[1] + z[2] + z[3])} /. Chip -> Chop[Chip]
```
Expected: `χ1` carries ~30-digit precision, `χ3` is a concrete complex number, and `χ1-(z1+z2+z3)` chops to `0`. (If `χ1` is still symbolic `z[1]+z[2]+z[3]` unevaluated, the re-run in Step 2 was skipped.)

- [ ] **Step 4: Add the `ν1` alias cell** (single source of truth for the Yangian weight of `T11`, mirroring su(2)'s `aFun`→`ν1` dedup):

```wolfram
ClearAll[ν1];
ν1[S1_, S2_][u_] := ν[S1, S2][1][u];   (* = (u-θ1-h S1)(u-θ2-h S2) *)
ν1[1, 1][u]
```
Expected output: `(-(1/3) + u) (-(1/7) + u)` shifted by `-h S`, i.e. `(u - 1/3 - 1)(u - 1/7 - 1)` = `(u-4/3)(u-8/7)` for `S1=S2=1`.

- [ ] **Step 5: Reconcile the sketched `Bax` cell.** Replace its body so it uses Section A's real names and drops the explicit `χ3` (kept inside `qdetT`). This cell is documentation of the equation form; the operational operator is built in Task 5.

```wolfram
(*3rd-order Baxter/TQ equation (reference form). Diagonal-twist transfer matrices τ1=tg[·][1,1],
  τ2=tg[·][2,1]; qdetT already carries χ3, so NO explicit χ3 here (would double-count).
  A1=τ1/ν1[u], A2=τ2/(ν1[u]ν1[u-h]), A3=qdetT/(ν1[u]ν1[u-h]ν1[u-2h]).*)
ClearAll[Bax];
Bax[S1_, S2_] :=
  f[u + h]
   - tg[S1, S2][1, 1][u]/ν1[S1, S2][u] f[u]
   + tg[S1, S2][2, 1][u]/(ν1[S1, S2][u] ν1[S1, S2][u - h]) f[u - h]
   - qdetT[S1, S2][u]/(ν1[S1, S2][u] ν1[S1, S2][u - h] ν1[S1, S2][u - 2 h]) f[u - 2 h];
Head[Bax[1, 1]]
```
Expected: `Plus` (evaluates without error to a symbolic combination of `f[u±…]`).

- [ ] **Step 6: Save + commit.**

```bash
git add Paul/Mathematica/XXX/Experiments/su3_V1.wb
git commit -m "su(3) Section B: numeric params + reconcile Baxter/ν1 names"
```

---

## Task 2: Commuting family — coefficient matrices + Cartan operators

**Files:**
- Modify: `su3_V1.wb` (add `tgCoeffs`, `H1op`, `H2op`, commutation check)

**Interfaces:**
- Consumes: `tg[S1,S2][1,1][u]`, `tg[S1,S2][2,1][u]`, `EE[S1,S2][i,j]`, `dim`, `χ1`, `χ2`.
- Produces:
  - `tgCoeffs[S1_,S2_]` → `{c1, c2}` with `c1 = {C0,C1,C2}` (3 matrices, `τ1(u)=Σ C_k u^{k-1}`... i.e. `c1[[k]]` is the `u^{k-1}` coefficient matrix) and `c2 = {D0,…,D4}` (5 matrices).
  - `H1op[S1_,S2_]`, `H2op[S1_,S2_]` → `d×d` matrices, `d = dim[{S1,0,0}] dim[{S2,0,0}]`.

- [ ] **Step 1: Define `tgCoeffs` (memoized), entrywise `CoefficientList` in a dummy var** (mirrors su(2) `t2Coeffs`, extended to both transfer matrices):

```wolfram
ClearAll[tgCoeffs];
tgCoeffs[S1_, S2_] := tgCoeffs[S1, S2] = Module[{uu, m1, m2, cl1, cl2},
   m1 = tg[S1, S2][1, 1][uu];
   m2 = tg[S1, S2][2, 1][uu];
   cl1 = Map[PadRight[CoefficientList[#, uu], 3] &, m1, {2}];
   cl2 = Map[PadRight[CoefficientList[#, uu], 5] &, m2, {2}];
   {Table[cl1[[All, All, k]], {k, 1, 3}], Table[cl2[[All, All, k]], {k, 1, 5}]}
];
```

- [ ] **Step 2: Define the Cartan operators:**

```wolfram
ClearAll[H1op, H2op];
H1op[S1_, S2_] := EE[S1, S2][1, 1] - EE[S1, S2][2, 2];
H2op[S1_, S2_] := EE[S1, S2][2, 2] - EE[S1, S2][3, 3];
```

- [ ] **Step 3: Write the leading-coefficient + commutation check cell** (this is the task's test — hard-assert):

```wolfram
Module[{S1 = 2, S2 = 1, d, c1, c2, C2, D4, mats, comm},
  d = dim[{S1, 0, 0}] dim[{S2, 0, 0}];
  {c1, c2} = tgCoeffs[S1, S2];
  C2 = c1[[3]]; D4 = c2[[5]];
  (* leading coeffs: τ1 ~ χ1 u^2, τ2 ~ χ2 u^4 *)
  If[Max[Abs[Flatten[C2 - χ1 IdentityMatrix[d]]]] > 10^-13, Print["C2 FAIL"]; Abort[]];
  If[Max[Abs[Flatten[D4 - χ2 IdentityMatrix[d]]]] > 10^-13, Print["D4 FAIL"]; Abort[]];
  (* full commuting family pairwise-commutes *)
  mats = Join[c1, c2, {H1op[S1, S2], H2op[S1, S2]}];
  comm = Max[Table[Max[Abs[Flatten[mats[[i]].mats[[j]] - mats[[j]].mats[[i]]]]],
      {i, Length[mats]}, {j, Length[mats]}]];
  If[comm > 10^-12, Print["COMMUTATION FAIL: ", comm]; Abort[]];
  {"C2=χ1 Id OK", "D4=χ2 Id OK", "maxCommutator" -> comm}
]
```
Expected: `{"C2=χ1 Id OK", "D4=χ2 Id OK", "maxCommutator" -> <~1e-13>}` (no `Abort`).

- [ ] **Step 4: Run it; confirm no abort and small commutator.** If `COMMUTATION FAIL`, the transfer matrices were built in inconsistent frames — stop and investigate before proceeding.

- [ ] **Step 5: Save + commit.**

```bash
git add Paul/Mathematica/XXX/Experiments/su3_V1.wb
git commit -m "su(3) Section B: tgCoeffs + Cartan ops + commutation check"
```

---

## Task 3: Weight-multiplicity predictor `wtMult`

**Files:**
- Modify: `su3_V1.wb` (add `wtMult` + check)

**Interfaces:**
- Consumes: `dim`.
- Produces: `wtMult[S1_,S2_][H1_,H2_]` → non-negative integer, the multiplicity of the weight with `H1=n1-n2, H2=n2-n3` in `[S1,0,0]⊗[S2,0,0]`.

- [ ] **Step 1: Define `wtMult`** (recover `(n1,n2,n3)` from `(H1,H2)` and `n1+n2+n3=S1+S2`, then count `(a1,a2,a3)` splittings):

```wolfram
ClearAll[wtMult];
wtMult[S1_, S2_][H1_, H2_] := Module[{Ntot = S1 + S2, n1, n2, n3},
   n3 = (Ntot - H1 - 2 H2)/3; n2 = n3 + H2; n1 = n3 + H1 + H2;
   If[! (IntegerQ[n1] && IntegerQ[n2] && IntegerQ[n3] && n1 >= 0 && n2 >= 0 && n3 >= 0),
     Return[0]];
   Length@Select[
     Flatten[Table[{a1, a2, S1 - a1 - a2}, {a1, 0, S1}, {a2, 0, S1 - a1}], 1],
     (n1 - #[[1]] >= 0 && n2 - #[[2]] >= 0 && n3 - #[[3]] >= 0) &]
];
```

- [ ] **Step 2: Write the sum-rule check cell** (total over all weights = dim of tensor product):

```wolfram
Module[{S1 = 2, S2 = 1, d, tot},
  d = dim[{S1, 0, 0}] dim[{S2, 0, 0}];
  tot = Total@Flatten@Table[
      wtMult[S1, S2][H1, H2], {H1, -3 (S1 + S2), 3 (S1 + S2)}, {H2, -3 (S1 + S2), 3 (S1 + S2)}];
  If[tot != d, Print["SUM-RULE FAIL: ", tot, " != ", d]; Abort[]];
  {"total weights" -> tot, "dim" -> d,
   "spot (2,1)@(H1,H2)=(1,1)" -> wtMult[2, 1][1, 1]}
]
```
Expected: `total weights == dim` (no abort). For `(S1,S2)=(2,1)`, `d = 6·3 = 18`; the printed spot value is a small positive integer.

- [ ] **Step 3: Run; confirm the sum rule holds for one more rep** — re-run the check cell with `S1=2,S2=2` (`d=36`) to guard the combinatorics. Expected: no abort.

- [ ] **Step 4: Save + commit.**

```bash
git add Paul/Mathematica/XXX/Experiments/su3_V1.wb
git commit -m "su(3) Section B: wtMult weight-multiplicity predictor"
```

---

## Task 4: Eigensystem builder `TauEigensystem3` + accessors

**Files:**
- Modify: `su3_V1.wb` (add `TauEigensystem3`, `Psi`, `tau1`, `tau2`)

**Interfaces:**
- Consumes: `tgCoeffs`, `H1op`, `H2op`, `wtMult`, `dim`.
- Produces:
  - `TauEigensystem3[S1_,S2_]` → memoized `Association`, keys `{H1,H2,n}` → `<|"tau1c"->{3 nums}, "tau2c"->{5 nums}, "Psi"->vector|>`.
  - `Psi[S1_,S2_][H1_,H2_,n_]` → eigenvector (or `Missing`/`$Failed` on bad indices).
  - `tau1[S1_,S2_][H1_,H2_,n_][u_]`, `tau2[S1_,S2_][H1_,H2_,n_][u_]` → eigenvalue polynomials in `u`.

- [ ] **Step 1: Define `TauEigensystem3`** (generic-combination diagonalization of the full commuting family, Rayleigh-quotient extraction with residual asserts, weight labels read directly from `H1op,H2op`):

```wolfram
ClearAll[TauEigensystem3];
TauEigensystem3::rankdef = "Generic combination failed to resolve the `1`-dim eigenspace for (S1,S2)=(`2`,`3`) at every tried coefficient set.";
TauEigensystem3::notgenuine = "Eigenvector #`1` for (S1,S2)=(`2`,`3`) failed the common-eigenvector residual check (residual=`4`).";
TauEigensystem3::noninteger = "Weight (H1 or H2)=`1` not near-integer for eigenvector #`2` at (S1,S2)=(`3`,`4`).";
TauEigensystem3::degnull = "Sector (H1,H2)=(`1`,`2`) for (S1,S2)=(`3`,`4`) has `5` state(s) but wtMult predicts `6`.";
TauEigensystem3[S1_, S2_] := TauEigensystem3[S1, S2] = Module[
   {d, c1, c2, mats, nontriv, rSets, rset, Kmix, vals, vecs, rayleigh, recs, byHW, table},
   d = dim[{S1, 0, 0}] dim[{S2, 0, 0}];
   {c1, c2} = tgCoeffs[S1, S2];
   (* full family for extraction; drop the two scalar leading coeffs (c1[[3]],c2[[5]]) from the
      splitting combination since they carry no state info *)
   mats = Join[c1, c2, {H1op[S1, S2], H2op[S1, S2]}];
   nontriv = Join[c1[[1 ;; 2]], c2[[1 ;; 4]], {H1op[S1, S2], H2op[S1, S2]}]; (* 8 matrices *)
   rSets = Table[N[Table[Prime[100 i + 7 j]^(1/3), {j, Length[nontriv]}], 30], {i, 1, 8}];
   rset = SelectFirst[rSets,
     Module[{ee, vv}, {ee, vv} = Eigensystem[#.nontriv];
        MatrixRank[vv] == d && Length[Union[Chop[ee, 10^-12]]] == d] &,
     Missing["NoGenericR"]];
   If[MissingQ[rset], Message[TauEigensystem3::rankdef, d, S1, S2]; Abort[]];
   Kmix = rset . nontriv;
   {vals, vecs} = Eigensystem[N[Kmix, 30]];
   rayleigh[m_, v_] := (v . (m . v))/(v . v);
   recs = Table[
     Module[{v = vecs[[i]], quotients, res, H1v, H2v, t1c, t2c},
       quotients = rayleigh[#, v] & /@ mats;
       res = Max[Abs[#.v - rayleigh[#, v] v]] & /@ mats // Max;
       If[res > 10^-13, Message[TauEigensystem3::notgenuine, i, S1, S2, res]; Abort[]];
       t1c = rayleigh[#, v] & /@ c1;           (* {C0,C1,C2} quotients *)
       t2c = rayleigh[#, v] & /@ c2;           (* {D0..D4} quotients *)
       H1v = rayleigh[H1op[S1, S2], v]; H2v = rayleigh[H2op[S1, S2], v];
       If[Abs[H1v - Round[Re[H1v]]] > 10^-10, Message[TauEigensystem3::noninteger, H1v, i, S1, S2]; Abort[]];
       If[Abs[H2v - Round[Re[H2v]]] > 10^-10, Message[TauEigensystem3::noninteger, H2v, i, S1, S2]; Abort[]];
       <|"H1" -> Round[Re[H1v]], "H2" -> Round[Re[H2v]], "t1c" -> t1c, "t2c" -> t2c, "vec" -> v|>],
     {i, 1, d}];
   byHW = GroupBy[recs, {#["H1"], #["H2"]} &];
   table = Association @ Flatten @ KeyValueMap[
      Function[{hw, group},
        Module[{sorted, predicted, last},
          predicted = wtMult[S1, S2][hw[[1]], hw[[2]]];
          (* order n by (Re,Im) of the C0 Rayleigh quotient (energy-like scalar) *)
          sorted = SortBy[group, {Re[#["t1c"][[1]]], Im[#["t1c"][[1]]]} &];
          If[Length[sorted] != predicted,
            Message[TauEigensystem3::degnull, hw[[1]], hw[[2]], S1, S2, Length[sorted], predicted]; Abort[]];
          MapIndexed[
            Function[{rec, idx}, Module[{n = First[idx] - 1, vv = rec["vec"], lc},
              lc = SelectFirst[vv, Abs[#] > 10^-13 Max[Abs[vv]] &];
              {hw[[1]], hw[[2]], n} -> <|"tau1c" -> rec["t1c"], "tau2c" -> rec["t2c"], "Psi" -> vv/lc|>]],
            sorted]]],
      byHW];
   If[Length[table] != d, Message[TauEigensystem3::degnull, "total", "-", S1, S2, Length[table], d]; Abort[]];
   table
];
```

- [ ] **Step 2: Define the accessors:**

```wolfram
ClearAll[Psi, tau1, tau2];
Psi::bnds = "No state at (H1,H2,n)=(`1`,`2`,`3`) for (S1,S2)=(`4`,`5`).";
Psi[S1_, S2_][H1_, H2_, n_] := Lookup[TauEigensystem3[S1, S2], Key[{H1, H2, n}],
    (Message[Psi::bnds, H1, H2, n, S1, S2]; Missing["NoState", {H1, H2, n}])]["Psi"];
tau1[S1_, S2_][H1_, H2_, n_][u_] := TauEigensystem3[S1, S2][{H1, H2, n}, "tau1c"] . u^Range[0, 2];
tau2[S1_, S2_][H1_, H2_, n_][u_] := TauEigensystem3[S1, S2][{H1, H2, n}, "tau2c"] . u^Range[0, 4];
```

- [ ] **Step 3: Write the eigensystem self-check cell** (residual of the actual transfer matrix against the reconstructed eigenvalue polynomial, per state):

```wolfram
Module[{S1 = 2, S2 = 1, es, worst},
  es = TauEigensystem3[S1, S2];
  worst = Max@KeyValueMap[
     Function[{k, rec},
       Module[{H1 = k[[1]], H2 = k[[2]], n = k[[3]], v = rec["Psi"], uu = 17/5},
         Max[
          Max[Abs[(tg[S1, S2][1, 1][uu] . v) - tau1[S1, S2][H1, H2, n][uu] v]],
          Max[Abs[(tg[S1, S2][2, 1][uu] . v) - tau2[S1, S2][H1, H2, n][uu] v]]]]],
     es];
  If[worst > 10^-11, Print["EIGEN RESIDUAL FAIL: ", worst]; Abort[]];
  {"nstates" -> Length[es], "worstResidual@u=17/5" -> worst}
]
```
Expected: `nstates == 18` for `(2,1)`, `worstResidual` `< 10^-11` (no abort).

- [ ] **Step 4: Run for `(2,1)`; then re-run the check with `S1=S2=1`** (`d=9`) and `S1=S2=2` (`d=36`) to confirm sector counts and residuals across several reps. Expected: no aborts; `nstates` = 9 and 36 respectively.

- [ ] **Step 5: Save + commit.**

```bash
git add Paul/Mathematica/XXX/Experiments/su3_V1.wb
git commit -m "su(3) Section B: TauEigensystem3 eigenvector builder + accessors"
```

---

## Task 5: Baxter operator + degree formula

**Files:**
- Modify: `su3_V1.wb` (add `BaxOp`, `alphaCoeffs`, `Mdeg`)

**Interfaces:**
- Consumes: `tau1`, `tau2`, `qdetT`, `ν1`, `z`, `χ2`, `χ3`, `h`.
- Produces:
  - `BaxOp[S1_,S2_][H1_,H2_,n_][Qf_][u_]` — the TQ residual functional (takes a pure function `Qf`).
  - `alphaCoeffs[S1_,S2_][H1_,H2_,n_]` → `{α1,α2,α3}` (the `1/u` coefficients of `A1,A2,A3`).
  - `Mdeg[S1_,S2_][H1_,H2_,n_]` → `{M1,M2,M3}` predicted polynomial degrees (one per `z[k]`), integers.

- [ ] **Step 1: Define the TQ residual functional** (operator form used by verification; `Qf` is a pure function of `u`):

```wolfram
ClearAll[BaxOp];
BaxOp[S1_, S2_][H1_, H2_, n_][Qf_][u_] :=
  Qf[u + h]
   - tau1[S1, S2][H1, H2, n][u]/ν1[S1, S2][u] Qf[u]
   + tau2[S1, S2][H1, H2, n][u]/(ν1[S1, S2][u] ν1[S1, S2][u - h]) Qf[u - h]
   - qdetT[S1, S2][u]/(ν1[S1, S2][u] ν1[S1, S2][u - h] ν1[S1, S2][u - 2 h]) Qf[u - 2 h];
```

- [ ] **Step 2: Define `alphaCoeffs`** (large-`u` `1/u` coefficients of the rational `A1,A2,A3`):

```wolfram
ClearAll[alphaCoeffs];
alphaCoeffs[S1_, S2_][H1_, H2_, n_] := alphaCoeffs[S1, S2][H1, H2, n] = Module[{u, A1, A2, A3},
   A1 = tau1[S1, S2][H1, H2, n][u]/ν1[S1, S2][u];
   A2 = tau2[S1, S2][H1, H2, n][u]/(ν1[S1, S2][u] ν1[S1, S2][u - h]);
   A3 = qdetT[S1, S2][u]/(ν1[S1, S2][u] ν1[S1, S2][u - h] ν1[S1, S2][u - 2 h]);
   SeriesCoefficient[#, {u, Infinity, 1}] & /@ {A1, A2, A3}
];
```

> **Corrected during execution:** the index must be `+1`, not `-1`. At `x0=Infinity`,
> `SeriesCoefficient[f,{x,Infinity,k}]` returns the coefficient of `x^{-k}`, so the `1/u`
> coefficient is index `+1`; `-1` asks for the (nonexistent) `u^{+1}` term and silently returns
> `0`, which made `Mdeg={0,0,0}` everywhere while still passing the integer guard.

- [ ] **Step 3: Define `Mdeg` from the asymptotic master formula:**

```wolfram
ClearAll[Mdeg];
Mdeg[S1_, S2_][H1_, H2_, n_] := Mdeg[S1, S2][H1, H2, n] = Module[{a, M},
   a = alphaCoeffs[S1, S2][H1, H2, n];
   M = Table[
      (a[[1]] - a[[2]]/z[k] + a[[3]]/z[k]^2)/(h (z[k] - χ2/z[k] + 2 χ3/z[k]^2)),
      {k, 1, 3}];
   If[Max[Abs[M - Round[Re[M]]]] > 10^-9, Print["Mdeg NON-INTEGER: ", M]; Abort[]];
   Round[Re[M]]
];
```

- [ ] **Step 4: Write the degree-tabulation cell** (the "figure out the M_k relation" deliverable — tabulate `M_k` vs `(H1,H2)` and check the sum rule `M1+M2+M3` against the total number of boxes):

```wolfram
Module[{S1 = 2, S2 = 1, es, tab},
  es = TauEigensystem3[S1, S2];
  tab = KeyValueMap[
     Function[{k, rec},
       {"H1H2n" -> k, "M" -> Mdeg[S1, S2][k[[1]], k[[2]], k[[3]]],
        "sumM" -> Total[Mdeg[S1, S2][k[[1]], k[[2]], k[[3]]]]}],
     KeySort[es]];
  Column[tab]
]
```
Expected: every `M` a triple of non-negative integers with no abort; inspect the `(H1,H2) → M` correspondence to extract the closed-form relation. **Record the observed relation** (e.g. `M_k` as a linear map of `(H1,H2,S1,S2)`) as a comment in the notebook and in the plan's Execution outcome.

- [ ] **Step 5: Run; confirm all degrees integer** across `(2,1)`. Then re-run the tabulation for `(1,1)`, `(2,2)`, `(1,3)` and note whether the same closed form holds. If `Mdeg NON-INTEGER` fires, the `alphaCoeffs`/master-formula derivation is off — stop and re-derive before Task 6.

- [ ] **Step 6: Save + commit.**

```bash
git add Paul/Mathematica/XXX/Experiments/su3_V1.wb
git commit -m "su(3) Section B: Baxter operator + asymptotic degree formula Mdeg"
```

---

## Task 6: Q-function solver `QSolve3` + accessor

**Files:**
- Modify: `su3_V1.wb` (add `qpoly`, `QSolve3`, `Qfun`)

**Interfaces:**
- Consumes: `tau1`, `tau2`, `qdetT`, `ν1`, `z`, `h`, `Mdeg`.
- Produces:
  - `qpoly[clist_,uu_]` — polynomial from ascending coeff list.
  - `QSolve3[S1_,S2_][H1_,H2_,n_]` → `{q1,q2,q3}`, each an ascending monic coeff list (length `M_k+1`).
  - `Qfun[S1_,S2_][H1_,H2_,n_][k_][u_]` → `z[k]^(u/h) qpoly[q_k, u]`.

- [ ] **Step 1: Define `qpoly`** (dummy-symbol trick avoids `0^0`, from su(2)):

```wolfram
ClearAll[qpoly];
qpoly[clist_, uu_] := Module[{x}, (clist . x^Range[0, Length[clist] - 1]) /. x -> uu];
```

- [ ] **Step 2: Define `QSolve3`** (per exponent `z[k]`: build the polynomial-identity coefficient matrix at degree `M_k`, nullspace, monic-normalize; `Chop` before `NullSpace` — same precision pitfall as su(2)'s `QSolve`):

```wolfram
ClearAll[QSolve3];
QSolve3::degnull = "Q_`5` nullspace dim `6` (expected 1) at (H1,H2,n)=(`1`,`2`,`3`) for (S1,S2)=(S1,S2).";
QSolve3::lastzero = "Q_`4` top coeff ~0 at (H1,H2,n)=(`1`,`2`,`3`).";
QSolve3[S1_, S2_][H1_, H2_, n_] := QSolve3[S1, S2][H1, H2, n] = Module[
   {Ms, u, T1, T2, QD, N1, N1m, N1mm, sols},
   Ms = Mdeg[S1, S2][H1, H2, n];
   T1 = tau1[S1, S2][H1, H2, n][u]; T2 = tau2[S1, S2][H1, H2, n][u];
   QD = qdetT[S1, S2][u];
   N1 = ν1[S1, S2][u]; N1m = ν1[S1, S2][u - h]; N1mm = ν1[S1, S2][u - 2 h];
   sols = Table[
     Module[{M = Ms[[k]], zk = z[k], cols, mat, ns, vec, last, maxdeg},
       maxdeg = M + 6;
       (* column j (j=0..M): coefficient list of the polynomial identity for q=u^j *)
       cols = Table[
          PadRight[CoefficientList[Expand[
             zk (u + h)^j N1 N1m N1mm
              - T1 u^j N1m N1mm
              + T2/zk (u - h)^j N1mm
              - QD/zk^2 (u - 2 h)^j], u], maxdeg + 1],
          {j, 0, M}];
       mat = Transpose[cols];                 (* (maxdeg+1) x (M+1) *)
       ns = NullSpace[Chop[mat, 10^-12], Tolerance -> 10^-10];
       If[Length[ns] != 1, Message[QSolve3::degnull, H1, H2, n, k, k, Length[ns]]; Abort[]];
       vec = ns[[1]]; last = vec[[-1]];
       If[Abs[last] < 10^-13 Max[Abs[vec]], Message[QSolve3::lastzero, H1, H2, n, k]; Abort[]];
       vec/last],
     {k, 1, 3}];
   sols
];
```

- [ ] **Step 3: Define the Q accessor:**

```wolfram
ClearAll[Qfun];
Qfun[S1_, S2_][H1_, H2_, n_][k_][u_] := z[k]^(u/h) qpoly[QSolve3[S1, S2][H1, H2, n][[k]], u];
```

- [ ] **Step 4: Write the single-state TQ-residual check** (each `Q_k` annihilated by `BaxOp`, tested at several numeric `u`):

```wolfram
Module[{S1 = 2, S2 = 1, H1 = 1, H2 = 1, n = 0, res},
  res = Table[
     Max@Table[
        Abs@Chop[N[BaxOp[S1, S2][H1, H2, n][Function[uu, Qfun[S1, S2][H1, H2, n][k][uu]]][u0], 30]],
        {u0, {2/3, 9/5, 4, 23/7}}],
     {k, 1, 3}];
  If[Max[res] > 10^-10, Print["TQ RESIDUAL FAIL: ", res]; Abort[]];
  {"perQresidual" -> res}
]
```
Expected: `perQresidual` a triple of numbers `< 10^-10` (no abort). (`(H1,H2)=(1,1)` exists for `(2,1)`; if not, pick any key from `TauEigensystem3[2,1]`.)

- [ ] **Step 5: Run; confirm all three Q's solve the TQ equation.** If a residual fails, check the `QSolve3` sign pattern against `BaxOp` and the `Mdeg` value for that `k`.

- [ ] **Step 6: Save + commit.**

```bash
git add Paul/Mathematica/XXX/Experiments/su3_V1.wb
git commit -m "su(3) Section B: QSolve3 three-Q Baxter solver + Qfun accessor"
```

---

## Task 7: Verification suite over the rep sweep

**Files:**
- Modify: `su3_V1.wb` (add the sweep verification cell)

**Interfaces:**
- Consumes: everything above (`TauEigensystem3`, `Qfun`, `BaxOp`, `tau1`, `tau2`, `qdetT`, `ν1`, `Mdeg`, `wtMult`).
- Produces: a printed pass/fail report; no new definitions.

- [ ] **Step 1: Write the sweep verification cell.** For each `(S1,S2)` and each state: (a) per-`Q_k` TQ residual `≈0`; (b) the 4×4 Casoratian reconstruction — the three `Q_k` span the TQ solution space, so the reconstructed `A1,A2,A3` (from 3×3 minors over the four shift rows) must match the coefficients built from `τ1,τ2,qdetT`.

```wolfram
Module[{sweep = {{1, 1}, {2, 1}, {1, 2}, {2, 2}, {3, 1}, {1, 3}}, report},
  report = Table[
    Module[{S1 = rep[[1]], S2 = rep[[2]], es, worstTQ = 0, worstCas = 0},
      es = TauEigensystem3[S1, S2];
      KeyValueMap[
       Function[{key, rec},
        Module[{H1 = key[[1]], H2 = key[[2]], n = key[[3]], Qs, u0, rows, minors,
          A1rec, A2rec, A3rec, A1, A2, A3, N1, N1m, N1mm},
         (* (a) TQ residual per Q_k at sample points *)
         Do[worstTQ = Max[worstTQ,
             Abs@Chop[N[BaxOp[S1, S2][H1, H2, n][
                 Function[uu, Qfun[S1, S2][H1, H2, n][k][uu]]][u0], 30]]],
           {k, 1, 3}, {u0, {5/4, 11/3}}];
         (* (b) Casoratian reconstruction of A1,A2,A3 at a sample u0 *)
         u0 = 13/5;
         Qs[uu_] := Table[Qfun[S1, S2][H1, H2, n][k][uu], {k, 1, 3}];
         rows = {Qs[u0 + h], Qs[u0], Qs[u0 - h], Qs[u0 - 2 h]}; (* 4 x 3 *)
         (* minor_i = det of the 3 rows OTHER than row i (i=1..4) *)
         minors = Table[Det[Delete[rows, i]], {i, 1, 4}];
         (* TQ from Casoratian: Q(u+h) - A1 Q(u) + A2 Q(u-h) - A3 Q(u-2h) = 0
            ⇒ A1 = minor2/minor1, A2 = minor3/minor1, A3 = minor4/minor1 (signs from cofactor) *)
         A1rec = minors[[2]]/minors[[1]];
         A2rec = minors[[3]]/minors[[1]];
         A3rec = minors[[4]]/minors[[1]];
         N1 = ν1[S1, S2][u0]; N1m = ν1[S1, S2][u0 - h]; N1mm = ν1[S1, S2][u0 - 2 h];
         A1 = tau1[S1, S2][H1, H2, n][u0]/N1;
         A2 = tau2[S1, S2][H1, H2, n][u0]/(N1 N1m);
         A3 = qdetT[S1, S2][u0]/(N1 N1m N1mm);
         worstCas = Max[worstCas,
            Abs@Chop[N[A1rec - A1, 30]], Abs@Chop[N[A2rec - A2, 30]], Abs@Chop[N[A3rec - A3, 30]]]]],
       es];
      {rep, "nstates" -> Length[es], "worstTQ" -> worstTQ, "worstCasoratian" -> worstCas}],
    {rep, sweep}];
  If[Max[#[[3, 2]] & /@ report, #[[4, 2]] & /@ report] > 10^-9,
    Print["SWEEP FAIL: ", report]; Abort[]];
  Column[report]
]
```
Expected: one row per rep, each with `worstTQ` and `worstCasoratian` `< 10^-9`; `nstates` equal to `dim[{S1,0,0}]·dim[{S2,0,0}]` (9, 18, 18, 36, 30, 30). No abort.

- [ ] **Step 2: Run the full sweep.** If the Casoratian check fails but the TQ residual passes, the sign/cofactor convention in the `A_i = minor_{i+1}/minor_1` mapping is off — adjust signs so the reconstructed `A_i` match (the physics is already validated by the passing TQ residual). If TQ fails for a specific rep, investigate `Mdeg`/`QSolve3` for that rep.

- [ ] **Step 3: Record the outcome.** Add an `## Execution outcome` section to this plan and to `docs/superpowers/specs/2026-07-17-...-design.md` noting: the observed `M_k(H1,H2,S1,S2)` closed form (from Task 5), the max residuals over the sweep, and the number of states verified.

- [ ] **Step 4: Save + commit.**

```bash
git add Paul/Mathematica/XXX/Experiments/su3_V1.wb docs/superpowers/plans/2026-07-17-su3-transfer-matrix-eigensystem-baxter.md docs/superpowers/specs/2026-07-17-su3-transfer-matrix-eigensystem-baxter-design.md
git commit -m "su(3) Section B: full-sweep verification (TQ residual + Casoratian) + outcome"
```

- [ ] **Step 5: Update `HANDOFF.md`** with the new Section-B state (eigensystem + Q-functions built and verified for gl(3), the `M_k` relation found, what's next — companion-twist SoV / CG-overlap for gl(3)). Commit.

---

## Self-Review

- **Spec coverage:** Component 1 (eigenvector builder + spectrum) → Tasks 2–4. Component 2 (Baxter + degree relation) → Tasks 5–6. Component 3 (verification suite: TQ residual, Wronskian/QQ, eigenvalue consistency, sector counts, degree formula) → sector counts in Task 4, degree formula in Task 5, TQ + Casoratian (which jointly cover the Wronskian/QQ AND eigenvalue-consistency checks) in Tasks 6–7. Numeric parameters + `A3`-normalization decision → Task 1. All spec sections mapped.
- **Placeholder scan:** no TBD/TODO; every code step has complete Wolfram code. The one genuinely open research output — the closed-form `M_k(H1,H2,S1,S2)` — is computed and recorded in Task 5, not assumed.
- **Type consistency:** `TauEigensystem3` returns keys `{H1,H2,n}` with fields `"tau1c"`(len 3)/`"tau2c"`(len 5)/`"Psi"`; `tau1`/`tau2` dot those against `u^Range[0,2]`/`u^Range[0,4]` (consistent lengths). `Mdeg` returns a length-3 integer list; `QSolve3` indexes it `Ms[[k]]` and returns a length-3 list of coeff-lists; `Qfun[...][k]` indexes that `[[k]]` (consistent). `BaxOp` takes a pure function `Qf` and is fed `Function[uu, Qfun[...][k][uu]]` (consistent).

Note on the Casoratian cofactor signs (Task 7 Step 2): the exact `A_i = ± minor_{i+1}/minor_1` signs follow from Laplace expansion of the 4×4 vanishing determinant; the plan flags adjusting them to match, since the TQ residual (Task 6) independently certifies the Q's regardless of that bookkeeping.

---

## Execution outcome (2026-07-17)

All 7 tasks implemented in `Paul/Mathematica/XXX/Experiments/su3_V1.wb` (Section B, cells after
the Section-A `Bax` reference cell), executed via subagent-driven development, each hard-assert
verified live and committed to `main`.

**Result — everything zero-residual:**
- `TauEigensystem3` diagonalizes the diagonal-twist transfer matrices, labels every state by
  `(H1,H2,n)`, sector counts matching `wtMult` exactly.
- The three Baxter Q-functions `Q_k = z_k^{u/h} q_k(u)` are solved by `QSolve3`; each satisfies the
  3rd-order TQ equation with residual `0`.
- **Degree relation (the "figure it out" deliverable):** `deg q_k = M_k = N − n_k`, where
  `N = S1+S2` and `(n1,n2,n3)` are the gl(3) weight components (`n3=(N−H1−2H2)/3, n2=n3+H2,
  n1=n3+H1+H2`); equivalently `M1=(2N−2H1−H2)/3, M2=(2N+H1−H2)/3, M3=(2N+H1+2H2)/3`. `M_k` pairs
  with `z_k` in the same index, depends only on `(H1,H2)` (not the branch index `n`), and
  `M1+M2+M3 = 2N`. Hard-asserted to zero residual over `(2,1),(1,1),(2,2),(1,3)`.
- **Full sweep** `{(1,1),(2,1),(1,2),(2,2),(3,1),(1,3)}`, 141 states total: `worstTQ = 0` and
  `worstCasoratian = 0` for every rep (`nstates` = 9/18/18/36/30/30). The 4×4 Casoratian
  reconstruction of `A1,A2,A3` from the three `Q_k` matches the coefficients built from
  `τ1,τ2,qdetT` — the Wronskian/QQ consistency check.

**Deviations from the as-written plan (both improvements, verified):**
1. `alphaCoeffs` `SeriesCoefficient` index `-1 → +1` (see correction note above) — a real bug in the
   plan; without it every degree came out `0`.
2. `TauEigensystem3` `n`-ordering strengthened from a single C0 quotient to the full tuple of
   `τ1`/`τ2` conserved charges, guaranteeing a total order within a degenerate weight sector.

**Commits:** `b4b68c7` (Task 1) · `93ff66e` (Task 2) · `b783526` (Task 3) · `d09044c` (Task 4) ·
`46d1509` (Task 5) · `dd66326` (Task 6) · `645517d` (Task 7).
