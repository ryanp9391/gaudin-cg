# gl(4) L=2 Twisted Baxter Solver — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a two-site gl(4) XXX spin chain (generic Young-diagram rep per site, inhomogeneities θ, four distinct twist eigenvalues) and solve its 4th-order Baxter/TQ equation for four Q-functions, in `Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb`.

**Architecture:** Direct port of the validated gl(3) pipeline `su3_V2_general.wb` (40 cells) to gl(4). The GT engine `ps/dim/J0/Jp/Jm/JJ` is N-general (gl(4) generators for free); everything above it gains one index/term. Each task appends cells to the single notebook and ends with an in-kernel hard-assert check that must pass (the "test").

**Tech Stack:** Wolfram Language in a wolfbook `.wb` notebook, driven live via the **wolfbook-builder** subagent (wolfbook MCP tools: `insertCells`, `runCell`, `evaluateExpression`, `getCellOutput`, `fileOps` to save). No standalone `.wl` scripts.

## Global Constraints

- **File:** everything lives in one notebook — `Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb` (currently empty). Cells are appended in order; each task's cells come after the previous task's.
- **Numeric params:** `θ1 = 1/3`, `θ2 = 1/7`, `h = 1`, `z[k] = N[Exp[I Zeta[2 k + 1]], 50]` for `k = 1,2,3,4` (Zeta[3], Zeta[5], Zeta[7], Zeta[9]). (Precision raised from 30→50 digits during Task 10 — see deviation note there: 30 digits, matching the gl(3) convention, left the largest rep's Casoratian residual at ~4.5e-8 due to precision cascading through a 60×60 eigendecomposition + 4th-order nullspace solve, not a formula error — TQ residual was exactly 0 everywhere at 30 digits too.)
- **Twist:** four **distinct** unit-modulus eigenvalues; `χ1..χ4` = elementary symmetric functions of `z[1..4]`; **diagonal twist only** — no companion-twist matrix `G`/`TG`/`tG`/`qminG` (dropped 2026-07-19: unused by the gl(3) source pipeline too, which only ever solved the TQ equation off the diagonal-twist branch `Tg`/`tg`/`qming`; kept `TG` there purely for future convention. χ1..χ4 stay as plain elementary-symmetric-function definitions — still needed directly by `qdetT` and the leading-coeff/commutation check).
- **Rep interface:** every function takes full 4-vectors `λ1, λ2` (one per site); `λ = {λ1,λ2,λ3,λ4}`, `λ1≥λ2≥λ3≥λ4≥0` integers. No scalar-`S` wrappers.
- **Degree formula:** verify `M_k = Λ_1 − n_k` (`Λ = λ1+λ2`, `Λ_1 = (λ1+λ2)[[1]]`) holds; do NOT re-derive by candidate-fit unless it fails.
- **Verification convention:** hard-assert with `If[worst > tol, Print[...]; Abort[]]`; all residuals machine-zero (`< 10^-9`). Chop high-precision numerics before rank/nullspace/Union.
- **Sweep (5 rep pairs):** `{1,0,0,0}⊗{1,0,0,0}` (16), `{1,0,0,0}⊗{1,1,1,0}` (16), `{1,1,0,0}⊗{1,0,0,0}` (24), `{1,1,0,0}⊗{1,1,0,0}` (36, interior mult), `{2,1,1,0}⊗{1,0,0,0}` (60, adjoint, zero-weight mult 3).
- **Gotchas:** `;` between statements in multi-statement cells; after fixing a memoized `f[x_]:=f[x]=…` re-run the defining cell AND everything downstream (stale cache); save the notebook explicitly (`fileOps`) before committing; do not trust blank `Out[]` boxes — re-run.
- **Git:** the wolfbook-builder subagent does NOT commit. The controller saves the notebook and commits between tasks.

---

### Task 1: Section A foundations — GT engine, generators, Lax, central charge

**Files:**
- Modify: `Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb` (append cells)

**Interfaces:**
- Produces: `ps[{λ}]`, `dim[λ]`, `J0[k][pat]`, `Jp/Jm[k][pat]`, `JJ[k,l,rep]`, `Ee[λ][i,j]`, `id[λ]`, `idd[λ1,λ2]`, `EE[λ1,λ2][i,j]`, `C1[λ1,λ2]`, `L[λ][i,j][u]`. All indices `i,j,k,l ∈ 1..4`.

- [ ] **Step 1: Insert the GT-engine cell** (copied verbatim from `su3_V2_general.wb` cell 0 — it is N-general; only the top comment changes to `(*su(4) generators*)`):

```wolfram
(*su(4) generators*)
ps[rep_] := ps[rep] = Block[{res = (Prepend[rep, #] & /@ (Tuples[Table[Range[rep[[1, i + 1]], rep[[1, i]]], {i, 1, Length[rep[[1]]] - 1}]]))}, If[Length[res[[1, 1]]] == 1, Reverse@Sort@res, Reverse@Sort@Flatten[ps /@ res, 1]]];
dim[rep_] := dim[rep] = Length[ps[{rep}]];
J0[k_][pat_] := J0[k][pat] = Sum[pat[[k, i]], {i, 1, k}] - Sum[pat[[k - 1, i]], {i, 1, k - 1}];
Jp[k_][pat_] := Jp[k][pat] = -Sum[Product[pat[[k, i]] - i - pat[[k + 1, j]] + j, {j, Range[k + 1]}]/Product[pat[[k, i]] - i - pat[[k, j]] + j, {j, Complement[Range[k], {i}]}] v[ReplacePart[pat, {k, i} -> pat[[k, i]] + 1]], {i, 1, k}];
Jm[k_][pat_] := Jm[k][pat] = Sum[Product[pat[[k, i]] - i - pat[[k - 1, j]] + j, {j, Range[k - 1]}]/Product[pat[[k, i]] - i - pat[[k, j]] + j, {j, Complement[Range[k], {i}]}] v[ReplacePart[pat, {k, i} -> pat[[k, i]] - 1]], {i, 1, k}];
JJ[k_, l_, rep_] := JJ[k, l, rep] = Which[k == l, DiagonalMatrix[J0[k] /@ ps[{rep}]], l == k + 1, Transpose[Table[Coefficient[Jp[k][ps[{rep}][[i]]], v[#]] & /@ ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l == k - 1, Transpose[Table[Coefficient[Jm[l][ps[{rep}][[i]]], v[#]] & /@ ps[{rep}], {i, 1, Length[ps[{rep}]]}]], l > k + 1, JJ[k, k + 1, rep] . JJ[k + 1, l, rep] - JJ[k + 1, l, rep] . JJ[k, k + 1, rep], l < k - 1, JJ[k, k - 1, rep] . JJ[k - 1, l, rep] - JJ[k - 1, l, rep] . JJ[k, k - 1, rep]];
```

- [ ] **Step 2: Insert the generators/embedding cell** (from cell 1, unchanged):

```wolfram
ClearAll[Ee,id,EE,idd];
Ee[λ_][i_,j_]:=JJ[i,j,λ];
id[λ_]:=IdentityMatrix[dim[λ]];
idd[λ1_,λ2_]:=KroneckerProduct[id[λ1],id[λ2]];
EE[λ1_,λ2_][i_,j_]:=KroneckerProduct[Ee[λ1][i,j],id[λ2]]+KroneckerProduct[id[λ1],Ee[λ2][i,j]];
```

- [ ] **Step 3: Insert the central-charge cell** (from cell 2, `{i,3}` → `{i,4}`):

```wolfram
(*central charge*)
ClearAll[C1];
C1[λ1_,λ2_]:=Sum[EE[λ1,λ2][i,i],{i,4}];
```

- [ ] **Step 4: Insert the Lax cell** (from cell 4, comment gl(3)→gl(4), body unchanged):

```wolfram
(*Lax operator for a general gl(4) rep λ={λ1,λ2,λ3,λ4} at the site*)
ClearAll[L];
L[λ_][i_,j_][u_]:=IdentityMatrix[dim[λ]] u KroneckerDelta[i,j] - h JJ[j,i,λ];
```

- [ ] **Step 5: Insert and run the central-charge check** (from cell 3, gl(4) reps). This is the test:

```wolfram
(*checking central charge constraint = total boxes, on general gl(4) reps.*)
Module[{reps={{{1,1,0,0},{1,0,0,0}},{{2,1,1,0},{1,0,0,0}},{{1,1,0,0},{1,1,0,0}}}},
 Union@Flatten@Table[
   C1[r[[1]],r[[2]]]-(Total[r[[1]]]+Total[r[[2]]]) idd[r[[1]],r[[2]]],{r,reps}]]
```

Run via `wolfbook_runCell`. **Expected output: `{0}`** (every entry of `C1 − (total boxes)·Id` is exactly zero across all three test reps). Also sanity-check dims: run `evaluateExpression` on `{dim[{1,0,0,0}],dim[{1,1,0,0}],dim[{2,1,1,0}]}` → **Expected `{4,6,15}`**.

- [ ] **Step 6: Save and commit** (controller):

```bash
# wolfbook-builder saves via fileOps; then controller:
cd /home/paul/git_projects/gaudin-CG
git add Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb
git commit -m "su4 L2: Section A foundations (GT engine, generators, Lax, central charge)"
```

---

### Task 2: Diagonal twist (χ1..χ4) + numeric parameters + monodromies

**Files:**
- Modify: `Baxter_L2_XXX_SU4.wb`

**Interfaces:**
- Consumes: `L`, `dim` (Task 1).
- Produces: `χ1,χ2,χ3,χ4`, `z[1..4]`, `θ1,θ2,h`, `T[λ1,λ2][i,j][u]`, `Tg[λ1,λ2][i,j][u]`. Indices `i,j ∈ 1..4`.

**Deviation from the design spec (confirmed with the user 2026-07-19):** companion twist is dropped entirely for this project — no `G` matrix, no `TG`/`tG`/`qminG`. The gl(3) source pipeline only ever solved the TQ equation off the diagonal-twist branch (`Tg`/`tg`/`qming`); `TG`/`tG`/`qminG` were defined there but unused downstream. `χ1..χ4` stay as plain elementary-symmetric-function definitions of `z[1..4]` — still needed directly by `qdetT` (Task 4) and the leading-coeff/commutation check (Task 5).

- [ ] **Step 1: Insert the χ definitions cell** (elementary symmetric functions of z[1..4]; no companion matrix):

```wolfram
(*χ1..χ4: elementary symmetric functions of the four diagonal-twist eigenvalues z[1..4]*)
χ1=z[1]+z[2]+z[3]+z[4];
χ2=z[1]z[2]+z[1]z[3]+z[1]z[4]+z[2]z[3]+z[2]z[4]+z[3]z[4];
χ3=z[1]z[2]z[3]+z[1]z[2]z[4]+z[1]z[3]z[4]+z[2]z[3]z[4];
χ4=z[1]z[2]z[3]z[4];
```

- [ ] **Step 2: Insert the numeric-parameters cell** (θ,h from gl(3); four distinct unit-modulus twists, adding `z[4]`):

```wolfram
(*Numeric parameters (shared across Section B). θ,h as in the gl(3) notebook;
  four DISTINCT unit-modulus twist eigenvalues z[1..4] (generic → non-degenerate).
  χ1..χ4 were defined above as symmetric functions of z[i], so they become
  numeric automatically once z[i] are assigned.
  NB: 50-digit precision (raised from the gl(3)-matching 30 digits 2026-07-19) — the
  4th-order Baxter equation and larger reps here (up to d=60) chain a 60x60 high-precision
  eigendecomposition through several more linear-algebra stages than gl(3)'s 3rd-order/d<=24
  case, and 30 digits left the worst-case Casoratian residual at ~4.5e-8 (just outside the
  1e-9 hard-assert bar) even though the TQ residual itself was exactly 0 everywhere — i.e. a
  precision-cascade issue, not a formula error. 50 digits gives enough headroom.*)
θ1 = 1/3;
θ2 = 1/7;
h = 1;
z[1] = N[Exp[I Zeta[3]], 50];
z[2] = N[Exp[I Zeta[5]], 50];
z[3] = N[Exp[I Zeta[7]], 50];
z[4] = N[Exp[I Zeta[9]], 50];
{θ1, θ2, h, z[1], z[2], z[3], z[4], χ1, χ2, χ3, χ4}
```

- [ ] **Step 3: Insert the monodromy cell** (from cell 7, `{ii,3}` → `{ii,4}`; diagonal-twist branch only):

```wolfram
(*monodromy with no twist and diagonal twist*)
ClearAll[T,Tg];
T[λ1_,λ2_][i_,j_][u_]:=Sum[KroneckerProduct[L[λ1][i,ii][u-θ1],L[λ2][ii,j][u-θ2]],{ii,4}];
Tg[λ1_,λ2_][i_,j_][u_]:= z[i] T[λ1,λ2][i,j][u];
```

- [ ] **Step 4: Run cells 1→(this) in order, then sanity-check χ1..χ4 and distinctness of the twists** (the test — confirms the z/χ assignments are consistent). Run via `evaluateExpression`:

```wolfram
{Abs[χ1 - Total[{z[1],z[2],z[3],z[4]}]], Min[Abs[Subtract @@@ Subsets[{z[1],z[2],z[3],z[4]}, {2}]]]}
```

**Expected output:** `{r, s}` with `r` a number `< 10^-25` (χ1 numerically equals Σz, confirming z[i] were assigned before χ1's value was captured) and `s` an O(1) number (all four twists well-separated, so the spectrum will be non-degenerate).

- [ ] **Step 5: Save and commit** (controller):

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb
git commit -m "su4 L2: diagonal twist (chi1..chi4), numeric params (z4 added), monodromies"
```

---

### Task 3: Yangian weights ν + highest-weight-state check

**Files:**
- Modify: `Baxter_L2_XXX_SU4.wb`

**Interfaces:**
- Consumes: `T`, `θ1,θ2,h`, `dim` (Tasks 1–2).
- Produces: `ν[λ1,λ2][i][u]` (`i ∈ 1..4`), `HWS[λ1,λ2]`.

- [ ] **Step 1: Insert the Yangian-weights cell** (from cell 8; `ν` body already works for `i=1..4` since it indexes `λ1[[i]],λ2[[i]]`):

```wolfram
(*Yangian weights - eigenvalues of T[i,i] on the highest-weight state.
  General rep: T[i,i]|hw> = (u-θ1-h λ1[[i]])(u-θ2-h λ2[[i]]) |hw>, for i=1..4,
  from L[λ][i,i]|hw>=(u-h λ_i)|hw>.*)
ClearAll[ν,HWS];
ν[λ1_,λ2_][i_][u_]:=(u-θ1-h λ1[[i]])(u-θ2-h λ2[[i]]);
HWS[λ1_,λ2_]:=Table[{KroneckerDelta[1,ii]},{ii,1,dim[λ1] dim[λ2]}];
```

- [ ] **Step 2: Insert and run the HWS check** (from cell 9; gl(4) reps, `{i,3}` → `{i,4}`). This is the test:

```wolfram
(*checking highest-weights against the general ν, on non-symmetric gl(4) reps*)
Module[{reps={{{1,1,0,0},{1,0,0,0}},{{2,1,1,0},{1,0,0,0}},{{1,1,0,0},{1,1,0,0}}}},
 Union@Flatten@Table[
   T[r[[1]],r[[2]]][i,i][u].HWS[r[[1]],r[[2]]]-ν[r[[1]],r[[2]]][i][u] HWS[r[[1]],r[[2]]],
   {r,reps},{i,4}]]
```

**Expected output: `{0}`** (`T[i,i]|hw⟩ − ν[i]|hw⟩ = 0` for all i=1..4 and every test rep).

- [ ] **Step 3: Save and commit** (controller):

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb
git commit -m "su4 L2: Yangian weights nu (i=1..4) + HWS eigenvalue check"
```

---

### Task 4: Quantum minors + transfer matrices + qdetT + shift-ordering guard

**Files:**
- Modify: `Baxter_L2_XXX_SU4.wb`

**Interfaces:**
- Consumes: `T,Tg`, `ν`, `χ4`, `idd`, `h` (Tasks 1–3).
- Produces: `qmin/qming[λ1,λ2][II,JJ][u]`, `tg[λ1,λ2][a,1][u]` (`a ∈ 1,2,3,4`), `qdetT[λ1,λ2][u]`. `τ1=tg[·][1,1]` (deg 2), `τ2=tg[·][2,1]` (deg 4), `τ3=tg[·][3,1]` (deg 6), `tg[·][4,1]` = qdet (deg 8).

**Deviation from the design spec (companion twist dropped, per Task 2):** `qminG`/`tG` are not built — only the diagonal-twist branch (`qmin`/`qming`/`tg`) is needed downstream.

- [ ] **Step 1: Insert the quantum-minors cell** (from cell 10, diagonal-twist branch only):

```wolfram
(*quantum minors (q-determinants) of the monodromy matrix*)
ClearAll[qmin,qming];
qmin[λ1_, λ2_][II_List, JJ_List][u_] := qmin[λ1, λ2][II, JJ][u] = Sum[Signature[σ] Dot @@ Table[T[λ1, λ2][II[[σ[[k]]]], JJ[[k]]][u - (k - 1) h], {k, Length[II]}], {σ, Permutations[Range[Length[II]]]}];
qming[λ1_, λ2_][II_List, JJ_List][u_] := qming[λ1, λ2][II, JJ][u] = Sum[Signature[σ] Dot @@ Table[Tg[λ1, λ2][II[[σ[[k]]]], JJ[[k]]][u - (k - 1) h], {k, Length[II]}], {σ, Permutations[Range[Length[II]]]}];
```

- [ ] **Step 2: Insert the transfer-matrix cell** (from cell 11, subsets of `{1,2,3,4}`, diagonal-twist only):

```wolfram
(*transfer matrices with diagonal twist*)
ClearAll[tg];
tg[λ1_,λ2_][a_,1][u_]:= Total[qming[λ1, λ2][#, #][u] & /@ Subsets[{1, 2, 3, 4}, {a}]]
```

- [ ] **Step 3: Insert the qdetT cell** (gl(4): χ4 carried, 4-fold ν product with shifts 3h,2h,h,0):

```wolfram
(*value of quantum determinant (χ4 carried here; general ν)*)
ClearAll[qdetT];
qdetT[λ1_,λ2_][u_]:=χ4 ν[λ1,λ2][1][u-3h] ν[λ1,λ2][2][u-2h] ν[λ1,λ2][3][u-h] ν[λ1,λ2][4][u];
```

- [ ] **Step 4: Insert and run the shift-ordering guard** (from cell 13; guard is now `qdetT == tg[4,1]`). This is the test:

```wolfram
(*qdetT central value == tg[·][4,1]; on non-symmetric reps this is the ν shift-ordering guard.
  Chop wraps Simplify because z[k] (hence χ4) are high-precision numeric, leaving ~10^-27
  floating "0" coefficients on symbolic-u polynomials that Union never collapses to a literal 0
  even though every entry is numerically exactly 0.*)
Module[{reps={{{1,1,0,0},{1,0,0,0}},{{2,1,1,0},{1,0,0,0}},{{1,1,0,0},{1,1,0,0}}}},
 Union@Flatten@Table[
   Chop[Simplify[qdetT[r[[1]],r[[2]]][u] idd[r[[1]],r[[2]]]-tg[r[[1]],r[[2]]][4,1][u]],10^-15],{r,reps}]]
```

**Expected output: `{0}`** (qdetT matches the top minor `tg[4,1]` on all reps). If a nonzero O(1) residual appears, the ν shift ordering in Step 3 is wrong (the reversed order `ν[4][u-3h]…ν[1][u]` gives an O(1) residual — that is the discriminating signal, per the gl(3) note).

- [ ] **Step 5: Save and commit** (controller):

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb
git commit -m "su4 L2: quantum minors, transfer matrices tg[a,1], qdetT + shift-order guard"
```

---

### Task 5: tgCoeffs + Cartan operators H1,H2,H3 + leading-coeff/commutation check

**Files:**
- Modify: `Baxter_L2_XXX_SU4.wb`

**Interfaces:**
- Consumes: `tg`, `EE`, `χ1,χ2,χ3`, `dim` (Tasks 1–4).
- Produces: `tgCoeffs[λ1,λ2]` → `{c1,c2,c3}` (coeff-matrix lists of τ1,τ2,τ3; lengths 3,5,7; leading coeffs `c1[[3]],c2[[5]],c3[[7]]`), `H1op/H2op/H3op[λ1,λ2]`.

- [ ] **Step 1: Insert the tgCoeffs cell** (from cell 17; add τ3, pad to 3/5/7):

```wolfram
(*coefficient matrices of tau1=tg[.][1,1] (deg 2), tau2=tg[.][2,1] (deg 4), tau3=tg[.][3,1] (deg 6),
  extracted by CoefficientList in a dummy uu. c1={C0,C1,C2}, c2={D0..D4}, c3={G0..G6}.*)
ClearAll[tgCoeffs];
tgCoeffs[λ1_, λ2_] := tgCoeffs[λ1, λ2] = Module[{uu, m1, m2, m3, cl1, cl2, cl3},
   m1 = tg[λ1, λ2][1, 1][uu];
   m2 = tg[λ1, λ2][2, 1][uu];
   m3 = tg[λ1, λ2][3, 1][uu];
   cl1 = Map[PadRight[CoefficientList[#, uu], 3] &, m1, {2}];
   cl2 = Map[PadRight[CoefficientList[#, uu], 5] &, m2, {2}];
   cl3 = Map[PadRight[CoefficientList[#, uu], 7] &, m3, {2}];
   {Table[cl1[[All, All, k]], {k, 1, 3}], Table[cl2[[All, All, k]], {k, 1, 5}], Table[cl3[[All, All, k]], {k, 1, 7}]}
];
```

- [ ] **Step 2: Insert the Cartan cell** (from cell 18; add H3):

```wolfram
(*Cartan operators (gl(4) weight generators H1=E11-E22, H2=E22-E33, H3=E33-E44) on the two-site space*)
ClearAll[H1op, H2op, H3op];
H1op[λ1_, λ2_] := EE[λ1, λ2][1, 1] - EE[λ1, λ2][2, 2];
H2op[λ1_, λ2_] := EE[λ1, λ2][2, 2] - EE[λ1, λ2][3, 3];
H3op[λ1_, λ2_] := EE[λ1, λ2][3, 3] - EE[λ1, λ2][4, 4];
```

- [ ] **Step 3: Insert and run the leading-coeff + commutation check** (from cell 19; τ3 added, H3 added, tested on adjoint⊗fund d=60). This is the test:

```wolfram
(*hard-assert: τ1,τ2,τ3 leading coeffs are χ1 Id, χ2 Id, χ3 Id; commuting family pairwise-commutes.
  Tested on adjoint⊗fundamental (d=60), which has interior weight-multiplicity.*)
Module[{λ1={2,1,1,0},λ2={1,0,0,0},d,c1,c2,c3,C2,D4,G6,mats,comm},
  d=dim[λ1] dim[λ2];
  {c1,c2,c3}=tgCoeffs[λ1,λ2];
  C2=c1[[3]]; D4=c2[[5]]; G6=c3[[7]];
  If[Max[Abs[Flatten[C2-χ1 IdentityMatrix[d]]]]>10^-13,Print["C2 FAIL"];Abort[]];
  If[Max[Abs[Flatten[D4-χ2 IdentityMatrix[d]]]]>10^-13,Print["D4 FAIL"];Abort[]];
  If[Max[Abs[Flatten[G6-χ3 IdentityMatrix[d]]]]>10^-13,Print["G6 FAIL"];Abort[]];
  mats=Join[c1,c2,c3,{H1op[λ1,λ2],H2op[λ1,λ2],H3op[λ1,λ2]}];
  comm=Max[Table[Max[Abs[Flatten[mats[[i]].mats[[j]]-mats[[j]].mats[[i]]]]],
      {i,Length[mats]},{j,Length[mats]}]];
  If[comm>10^-12,Print["COMMUTATION FAIL: ",comm];Abort[]];
  {"C2=χ1 Id OK","D4=χ2 Id OK","G6=χ3 Id OK","maxCommutator"->comm}
]
```

**Expected output:** `{"C2=χ1 Id OK","D4=χ2 Id OK","G6=χ3 Id OK","maxCommutator"->x}` with `x < 10^-12`. No `FAIL`/`Abort`.

- [ ] **Step 4: Save and commit** (controller):

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb
git commit -m "su4 L2: tgCoeffs (tau1..tau3), Cartan H1..H3, leading-coeff+commutation check"
```

---

### Task 6: wtMult (gl(4) tensor-product weight tally) + sum-rule

**Files:**
- Modify: `Baxter_L2_XXX_SU4.wb`

**Interfaces:**
- Consumes: `J0`, `ps`, `dim` (Task 1).
- Produces: `siteWeightMult[λ]` (Association `{w1,w2,w3,w4} → mult`), `wtMult[λ1,λ2][H1,H2,H3]` (integer sector multiplicity).

- [ ] **Step 1: Insert the wtMult cell** (gl(4): 4-component weights; recover `(n1,n2,n3,n4)` from `(H1,H2,H3)` with `n4=(N-H1-2H2-3H3)/4`):

```wolfram
(*gl(4) weight multiplicities in the tensor product λ1⊗λ2.
  siteWeightMult[λ]: tally the gl(4) weight (J0[1..4]) of every GT pattern of λ.
  wtMult: recover (n1,n2,n3,n4) from (H1,H2,H3) with N=total boxes, then convolve the site tallies.
  H1=n1-n2, H2=n2-n3, H3=n3-n4, n1+n2+n3+n4=N ⇒ n4=(N-H1-2H2-3H3)/4.
  NB: wrap the composite weight key in Key[...] — Lookup treats a bare List as several keys.*)
ClearAll[siteWeightMult,wtMult];
siteWeightMult[λ_]:=siteWeightMult[λ]=Counts[({J0[1][#],J0[2][#],J0[3][#],J0[4][#]}&)/@ps[{λ}]];
wtMult[λ1_,λ2_][H1_,H2_,H3_]:=Module[{Ntot=Total[λ1]+Total[λ2],n1,n2,n3,n4,m1,m2},
  n4=(Ntot-H1-2 H2-3 H3)/4; n3=n4+H3; n2=n4+H2+H3; n1=n4+H1+H2+H3;
  If[!(IntegerQ[n1]&&IntegerQ[n2]&&IntegerQ[n3]&&IntegerQ[n4]),Return[0]];
  m1=siteWeightMult[λ1]; m2=siteWeightMult[λ2];
  Total[KeyValueMap[Function[{w,c}, c Lookup[m2,Key[{n1,n2,n3,n4}-w],0]],m1]]
];
```

- [ ] **Step 2: Insert and run the sum-rule check** (from cell 21; gl(4) sweep, extra H3 loop). This is the test:

```wolfram
(*hard-assert sum-rule: Σ over all (H1,H2,H3) of wtMult == dim[λ1] dim[λ2], for every sweep rep*)
Module[{reps={{{1,0,0,0},{1,0,0,0}},{{1,0,0,0},{1,1,1,0}},{{1,1,0,0},{1,0,0,0}},{{1,1,0,0},{1,1,0,0}},{{2,1,1,0},{1,0,0,0}}},bad},
  bad=Select[reps,Function[r,
     Module[{Ntot=Total[r[[1]]]+Total[r[[2]]],d,tot},
       d=dim[r[[1]]] dim[r[[2]]];
       tot=Total@Flatten@Table[wtMult[r[[1]],r[[2]]][H1,H2,H3],
           {H1,-Ntot,Ntot},{H2,-Ntot,Ntot},{H3,-Ntot,Ntot}];
       tot=!=d]]];
  If[bad=!={},Print["SUM-RULE FAIL: ",bad];Abort[]];
  {"all reps sum-rule OK",
   "adjoint⊗fund @(H1,H2,H3)=(1,0,0) mult"->wtMult[{2,1,1,0},{1,0,0,0}][1,0,0]}
]
```

**Expected output:** `{"all reps sum-rule OK", "adjoint⊗fund @(H1,H2,H3)=(1,0,0) mult"->m}` with `m` a positive integer (adjoint⊗fund exercises interior multiplicity). No `Abort`.

- [ ] **Step 3: Save and commit** (controller):

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb
git commit -m "su4 L2: wtMult (gl(4) weight tally) + tensor-dim sum-rule"
```

---

### Task 7: TauEigensystem4 + Psi/tau1/tau2/tau3 accessors + eigen-residual check

**Files:**
- Modify: `Baxter_L2_XXX_SU4.wb`

**Interfaces:**
- Consumes: `tgCoeffs`, `H1op/H2op/H3op`, `wtMult`, `dim`, `tg` (Tasks 4–6).
- Produces: `TauEigensystem4[λ1,λ2]` (Association `{H1,H2,H3,n} → <|tau1c,tau2c,tau3c,Psi|>`), `Psi[λ1,λ2][H1,H2,H3,n]`, `tau1/tau2/tau3[λ1,λ2][H1,H2,H3,n][u]`.

- [ ] **Step 1: Insert the TauEigensystem4 cell** (from cell 23; three commuting families c1,c2,c3, Cartan H1,H2,H3, `(H1,H2,H3,n)` keys, sort by all three tau charges):

```wolfram
ClearAll[TauEigensystem4];
TauEigensystem4::rankdef = "Generic combination failed to resolve the `1`-dim eigenspace for (λ1,λ2)=(`2`,`3`).";
TauEigensystem4::notgenuine = "Eigenvector #`1` for (λ1,λ2)=(`2`,`3`) failed the common-eigenvector residual check (residual=`4`).";
TauEigensystem4::noninteger = "Weight `1` not near-integer for eigenvector #`2` at (λ1,λ2)=(`3`,`4`).";
TauEigensystem4::degnull = "Sector (H1,H2,H3)=(`1`,`2`,`3`) for (λ1,λ2)=(`4`,`5`) has `6` state(s) but wtMult predicts `7`.";
TauEigensystem4[λ1_, λ2_] := TauEigensystem4[λ1, λ2] = Module[
   {d, c1, c2, c3, mats, nontriv, rSets, rset, Kmix, vals, vecs, rayleigh, recs, byHW, table},
   d = dim[λ1] dim[λ2];
   {c1, c2, c3} = tgCoeffs[λ1, λ2];
   (* full family for extraction; drop the scalar leading coeffs c1[[3]],c2[[5]],c3[[7]] from the
      splitting combination since they carry no state info *)
   mats = Join[c1, c2, c3, {H1op[λ1, λ2], H2op[λ1, λ2], H3op[λ1, λ2]}];
   nontriv = Join[c1[[1 ;; 2]], c2[[1 ;; 4]], c3[[1 ;; 6]], {H1op[λ1, λ2], H2op[λ1, λ2], H3op[λ1, λ2]}]; (* 15 matrices *)
   (* NB (fixed 2026-07-19): rSets and Kmix's eigendecomposition were hardcoded to 30-digit
      precision here (a leftover from the gl(3) source's TauEigensystem3), decoupled from the
      z[k]-precision bump in Task 2 (30->70) — this silently capped the whole downstream chain
      (tau-coefficients, Qfun, the Casoratian check) at ~19-25 accurate digits regardless of how
      much precision z carried, which is why raising z from 50->70 digits left the worst-case
      Casoratian residual completely unchanged (bit-identical) instead of shrinking. Both now
      match z's 70-digit precision. *)
   rSets = Table[N[Table[Prime[100 i + 7 j]^(1/3), {j, Length[nontriv]}], 70], {i, 1, 8}];
   rset = SelectFirst[rSets,
     Module[{ee, vv}, {ee, vv} = Eigensystem[#.nontriv];
        MatrixRank[vv] == d && Length[Union[Chop[ee, 10^-12]]] == d] &,
     Missing["NoGenericR"]];
   If[MissingQ[rset], Message[TauEigensystem4::rankdef, d, λ1, λ2]; Abort[]];
   Kmix = rset . nontriv;
   {vals, vecs} = Eigensystem[N[Kmix, 70]];
   rayleigh[m_, v_] := (v . (m . v))/(v . v);
   recs = Table[
     Module[{v = vecs[[i]], res, H1v, H2v, H3v, t1c, t2c, t3c},
       res = Max[Abs[#.v - rayleigh[#, v] v]] & /@ mats // Max;
       If[res > 10^-13, Message[TauEigensystem4::notgenuine, i, λ1, λ2, res]; Abort[]];
       t1c = rayleigh[#, v] & /@ c1;           (* {C0,C1,C2} quotients *)
       t2c = rayleigh[#, v] & /@ c2;           (* {D0..D4} quotients *)
       t3c = rayleigh[#, v] & /@ c3;           (* {G0..G6} quotients *)
       H1v = rayleigh[H1op[λ1, λ2], v]; H2v = rayleigh[H2op[λ1, λ2], v]; H3v = rayleigh[H3op[λ1, λ2], v];
       If[Abs[H1v - Round[Re[H1v]]] > 10^-10, Message[TauEigensystem4::noninteger, H1v, i, λ1, λ2]; Abort[]];
       If[Abs[H2v - Round[Re[H2v]]] > 10^-10, Message[TauEigensystem4::noninteger, H2v, i, λ1, λ2]; Abort[]];
       If[Abs[H3v - Round[Re[H3v]]] > 10^-10, Message[TauEigensystem4::noninteger, H3v, i, λ1, λ2]; Abort[]];
       <|"H1" -> Round[Re[H1v]], "H2" -> Round[Re[H2v]], "H3" -> Round[Re[H3v]], "t1c" -> t1c, "t2c" -> t2c, "t3c" -> t3c, "vec" -> v|>],
     {i, 1, d}];
   byHW = GroupBy[recs, {#["H1"], #["H2"], #["H3"]} &];
   table = Association @ Flatten @ KeyValueMap[
      Function[{hw, group},
        Module[{sorted, predicted},
          predicted = wtMult[λ1, λ2][hw[[1]], hw[[2]], hw[[3]]];
          (* deterministic order within a sector: full real+imag parts of all tau1/tau2/tau3 charges *)
          sorted = SortBy[group, Join[Re[#["t1c"]], Im[#["t1c"]], Re[#["t2c"]], Im[#["t2c"]], Re[#["t3c"]], Im[#["t3c"]]] &];
          If[Length[sorted] != predicted,
            Message[TauEigensystem4::degnull, hw[[1]], hw[[2]], hw[[3]], λ1, λ2, Length[sorted], predicted]; Abort[]];
          MapIndexed[
            Function[{rec, idx}, Module[{n = First[idx] - 1, vv = rec["vec"], lc},
              lc = SelectFirst[vv, Abs[#] > 10^-13 Max[Abs[vv]] &];
              {hw[[1]], hw[[2]], hw[[3]], n} -> <|"tau1c" -> rec["t1c"], "tau2c" -> rec["t2c"], "tau3c" -> rec["t3c"], "Psi" -> vv/lc|>]],
            sorted]]],
      byHW];
   If[Length[table] != d, Message[TauEigensystem4::degnull, "total", "-", "-", λ1, λ2, Length[table], d]; Abort[]];
   table
];
```

- [ ] **Step 2: Insert the Psi/tau accessors cell** (from cell 24; add tau3, keys now 4-tuples, `u^Range` extended):

```wolfram
ClearAll[Psi, tau1, tau2, tau3];
Psi::bnds = "No state at (H1,H2,H3,n)=(`1`,`2`,`3`,`4`) for (λ1,λ2)=(`5`,`6`).";
Psi[λ1_, λ2_][H1_, H2_, H3_, n_] := Lookup[TauEigensystem4[λ1, λ2], Key[{H1, H2, H3, n}],
    (Message[Psi::bnds, H1, H2, H3, n, λ1, λ2]; Missing["NoState", {H1, H2, H3, n}])]["Psi"];
tau1[λ1_, λ2_][H1_, H2_, H3_, n_][u_] := TauEigensystem4[λ1, λ2][{H1, H2, H3, n}, "tau1c"] . u^Range[0, 2];
tau2[λ1_, λ2_][H1_, H2_, H3_, n_][u_] := TauEigensystem4[λ1, λ2][{H1, H2, H3, n}, "tau2c"] . u^Range[0, 4];
tau3[λ1_, λ2_][H1_, H2_, H3_, n_][u_] := TauEigensystem4[λ1, λ2][{H1, H2, H3, n}, "tau3c"] . u^Range[0, 6];
```

- [ ] **Step 3: Insert and run the eigen-residual check** (from cell 25; add tg[3,1]/tau3, adjoint⊗fund d=60). This is the test:

```wolfram
Module[{λ1 = {2, 1, 1, 0}, λ2 = {1, 0, 0, 0}, es, worst},
  es = TauEigensystem4[λ1, λ2];
  worst = Max@KeyValueMap[
     Function[{k, rec},
       Module[{H1 = k[[1]], H2 = k[[2]], H3 = k[[3]], n = k[[4]], v = rec["Psi"], uu = 17/5},
         Max[
          Max[Abs[(tg[λ1, λ2][1, 1][uu] . v) - tau1[λ1, λ2][H1, H2, H3, n][uu] v]],
          Max[Abs[(tg[λ1, λ2][2, 1][uu] . v) - tau2[λ1, λ2][H1, H2, H3, n][uu] v]],
          Max[Abs[(tg[λ1, λ2][3, 1][uu] . v) - tau3[λ1, λ2][H1, H2, H3, n][uu] v]]]]],
     es];
  If[worst > 10^-11, Print["EIGEN RESIDUAL FAIL: ", worst]; Abort[]];
  {"nstates" -> Length[es], "worstResidual@u=17/5" -> worst}
]
```

**Expected output:** `{"nstates" -> 60, "worstResidual@u=17/5" -> w}` with `w < 10^-11`. `nstates` must equal `dim[{2,1,1,0}]·dim[{1,0,0,0}] = 15·4 = 60`. No `Abort`.

- [ ] **Step 4: Save and commit** (controller):

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb
git commit -m "su4 L2: TauEigensystem4 (simultaneous tau1..tau3 diag) + eigen-residual check"
```

---

### Task 8: BaxOp + alphaCoeffs + Mdeg master formula + degree-formula verification

**Files:**
- Modify: `Baxter_L2_XXX_SU4.wb`

**Interfaces:**
- Consumes: `tau1/tau2/tau3`, `qdetT`, `ν`, `z`, `χ2/χ3/χ4`, `h`, `TauEigensystem4` (Tasks 3–7).
- Produces: `ν1[λ1,λ2][u]`, `BaxOp[λ1,λ2][H1,H2,H3,n][Qf][u]`, `alphaCoeffs[λ1,λ2][H1,H2,H3,n]` → `{α1,α2,α3,α4}`, `Mdeg[λ1,λ2][H1,H2,H3,n]` → `{M1,M2,M3,M4}`, `nvec`, `Mpred[λ1,λ2][H1,H2,H3]`.

- [ ] **Step 1: Insert the ν1 single-source cell** (from cell 15):

```wolfram
(*single source of truth for the Yangian weight of T11*)
ClearAll[ν1];
ν1[λ1_,λ2_][u_]:=ν[λ1,λ2][1][u];
```

- [ ] **Step 2: Insert the BaxOp cell** (gl(4): 4th order, add the `+A4 Qf[u-3h]` term):

```wolfram
(*4th-order TQ residual functional (operator form). A1=τ1/ν1, A2=τ2/(ν1 ν1[-h]),
  A3=τ3/(ν1 ν1[-h] ν1[-2h]), A4=qdetT/(ν1 ν1[-h] ν1[-2h] ν1[-3h]); qdetT carries χ4.*)
ClearAll[BaxOp];
BaxOp[λ1_, λ2_][H1_, H2_, H3_, n_][Qf_][u_] :=
  Qf[u + h]
   - tau1[λ1, λ2][H1, H2, H3, n][u]/ν1[λ1, λ2][u] Qf[u]
   + tau2[λ1, λ2][H1, H2, H3, n][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h]) Qf[u - h]
   - tau3[λ1, λ2][H1, H2, H3, n][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h] ν1[λ1, λ2][u - 2 h]) Qf[u - 2 h]
   + qdetT[λ1, λ2][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h] ν1[λ1, λ2][u - 2 h] ν1[λ1, λ2][u - 3 h]) Qf[u - 3 h];
```

- [ ] **Step 3: Insert the alphaCoeffs cell** (from cell 27; four rational combos A1..A4, 1/u coefficient via `{u,Infinity,1}`):

```wolfram
(*large-u asymptotics: alphaCoeffs = {α1,α2,α3,α4}, the 1/u coefficients of A1..A4 in BaxOp.
  NB: for x0=Infinity, SeriesCoefficient[f,{u,Infinity,1}] is the coeff of u^(-1) — use n=+1.*)
ClearAll[alphaCoeffs];
alphaCoeffs[λ1_, λ2_][H1_, H2_, H3_, n_] := alphaCoeffs[λ1, λ2][H1, H2, H3, n] = Module[{u, A1, A2, A3, A4},
   A1 = tau1[λ1, λ2][H1, H2, H3, n][u]/ν1[λ1, λ2][u];
   A2 = tau2[λ1, λ2][H1, H2, H3, n][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h]);
   A3 = tau3[λ1, λ2][H1, H2, H3, n][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h] ν1[λ1, λ2][u - 2 h]);
   A4 = qdetT[λ1, λ2][u]/(ν1[λ1, λ2][u] ν1[λ1, λ2][u - h] ν1[λ1, λ2][u - 2 h] ν1[λ1, λ2][u - 3 h]);
   SeriesCoefficient[#, {u, Infinity, 1}] & /@ {A1, A2, A3, A4}
];
```

- [ ] **Step 4: Insert the Mdeg master-formula cell** (gl(4) derivation — numerator gains `-α4/z^3`, denominator gains `-3χ4/z^3`; k runs 1..4):

```wolfram
(*asymptotic master formula (gl(4)): for each twist eigenvalue z[k] (k=1..4) the large-u balance
  of the 4th-order TQ equation fixes the polynomial degree M[k]. Derived from
    M = -[Σ_j (-1)^j α_j z^{1-j}] / (h Σ_j (-1)^j χ_j (1-j) z^{1-j}), j=0..4,
  giving numerator (α1 - α2/z + α3/z^2 - α4/z^3) and
  denominator h (z - χ2/z + 2χ3/z^2 - 3χ4/z^3), the z-characteristic factor of the companion
  char poly z^4-χ1 z^3+χ2 z^2-χ3 z+χ4 differentiated and divided by z^3.*)
ClearAll[Mdeg];
Mdeg[λ1_, λ2_][H1_, H2_, H3_, n_] := Mdeg[λ1, λ2][H1, H2, H3, n] = Module[{a, M},
   a = alphaCoeffs[λ1, λ2][H1, H2, H3, n];
   M = Table[
      (a[[1]] - a[[2]]/z[k] + a[[3]]/z[k]^2 - a[[4]]/z[k]^3)/(h (z[k] - χ2/z[k] + 2 χ3/z[k]^2 - 3 χ4/z[k]^3)),
      {k, 1, 4}];
   If[Max[Abs[M - Round[Re[M]]]] > 10^-9, Print["Mdeg NON-INTEGER: ", M]; Abort[]];
   Round[Re[M]]
];
```

- [ ] **Step 5: Insert the degree-formula predictor** (gl(4) `nvec` and `Mpred = Λ_1 - n_k`):

```wolfram
(*closed-form degree: M_k = Λ_1 - n_k, Λ=λ1+λ2, Λ_1=(λ1+λ2)[[1]];
  n4=(N-H1-2H2-3H3)/4, n3=n4+H3, n2=n4+H2+H3, n1=n4+H1+H2+H3, N=Total[λ1]+Total[λ2].*)
ClearAll[nvec, Mpred];
nvec[H1_, H2_, H3_, Ntot_] := Module[{n4 = (Ntot - H1 - 2 H2 - 3 H3)/4}, {n4 + H1 + H2 + H3, n4 + H2 + H3, n4 + H3, n4}];
Mpred[λ1_, λ2_][H1_, H2_, H3_] := (λ1 + λ2)[[1]] - nvec[H1, H2, H3, Total[λ1] + Total[λ2]];
```

- [ ] **Step 6: Insert and run the degree-formula verification** (from cell 33; gl(4) sweep, 4-tuple keys). This is the test:

```wolfram
(*hard-assert: closed form M_k=Λ_1-n_k matches the asymptotic Mdeg exactly over the sweep.*)
Module[{sweep={{{1,0,0,0},{1,0,0,0}},{{1,0,0,0},{1,1,1,0}},{{1,1,0,0},{1,0,0,0}},{{1,1,0,0},{1,1,0,0}},{{2,1,1,0},{1,0,0,0}}},worst},
  worst=Max@Table[
     Module[{λ1=rep[[1]],λ2=rep[[2]],es},
       es=TauEigensystem4[λ1,λ2];
       Max@KeyValueMap[Function[{key,rec},
          Max[Abs[Mdeg[λ1,λ2][key[[1]],key[[2]],key[[3]],key[[4]]]-Mpred[λ1,λ2][key[[1]],key[[2]],key[[3]]]]]],es]],
     {rep,sweep}];
  If[worst>0,Print["DEGREE-FORMULA MISMATCH: ",worst];Abort[]];
  {"closed-form M_k=Λ_1-n_k matches Mdeg over full sweep, worst"->worst}
]
```

**Expected output:** `{"closed-form M_k=Λ_1-n_k matches Mdeg over full sweep, worst"->0}`. No `Abort`.

**Fallback (only if this aborts with a mismatch):** the `Λ_1−n_k` form does not carry to gl(4). Insert the gl(3) candidate-fit machinery (`degData` collection + `cand[...]` table, cells 29–30 adapted with the 4-component `nvec` and candidates `{"N-n","Lam1-n","Lamk-n"}`), find the zero-deviation form, replace `Mpred` with it, and re-run this check. Report the discovered form to the user before proceeding.

- [ ] **Step 7: Save and commit** (controller):

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb
git commit -m "su4 L2: BaxOp (4th order), alphaCoeffs, Mdeg master formula, degree-formula check"
```

---

### Task 9: qpoly + QSolve4 + Qfun + single-state TQ check

**Files:**
- Modify: `Baxter_L2_XXX_SU4.wb`

**Interfaces:**
- Consumes: `Mdeg`, `tau1/tau2/tau3`, `qdetT`, `ν1`, `z`, `h`, `TauEigensystem4`, `BaxOp` (Tasks 7–8).
- Produces: `qpoly[clist,uu]`, `QSolve4[λ1,λ2][H1,H2,H3,n]` → list of 4 monic coeff-vectors, `Qfun[λ1,λ2][H1,H2,H3,n][k][u]` = `z[k]^(u/h) q_k(u)`.

- [ ] **Step 1: Insert the qpoly cell** (from cell 34, unchanged):

```wolfram
(*polynomial from ascending coeff list {c0..cM}: qpoly[clist,uu]=Sum ck uu^k. x-substitution avoids 0^0.*)
ClearAll[qpoly];
qpoly[clist_, uu_] := Module[{x}, (clist . x^Range[0, Length[clist] - 1]) /. x -> uu];
```

- [ ] **Step 2: Insert the QSolve4 cell** (gl(4): four Q's; the cleared-denominator column adds the `-T3/z^2 (u-2h)^j N1mmm` and `+QD/z^3 (u-3h)^j` terms; denominator product is 4 ν1 factors; `maxdeg = M+8`):

```wolfram
(*QSolve4 solves the 4th-order Baxter (TQ) equation for the four monic Q-polynomials q_1..q_4
  (one per twist z[k]), degrees M_k from Mdeg. For fixed k, Q_k=z[k]^(u/h) q_k(u); substitute into
  BaxOp and multiply by ν1[u]ν1[u-h]ν1[u-2h]ν1[u-3h] to clear denominators, giving a polynomial
  identity in u linear in the coeffs of q_k. Each q=u^j is one column; NullSpace (Chop first) gives
  the coeff vector (1-dim), monic-normalized by the top coefficient.*)
ClearAll[QSolve4];
QSolve4::degnull = "Q_`5` nullspace dim `6` (expected 1) at (H1,H2,H3,n)=(`1`,`2`,`3`,`4`).";
QSolve4::lastzero = "Q_`5` top coeff ~0 at (H1,H2,H3,n)=(`1`,`2`,`3`,`4`).";
QSolve4[λ1_, λ2_][H1_, H2_, H3_, n_] := QSolve4[λ1, λ2][H1, H2, H3, n] = Module[
   {Ms, u, T1, T2, T3, QD, N1, N1m, N1mm, N1mmm, sols},
   Ms = Mdeg[λ1, λ2][H1, H2, H3, n];
   T1 = tau1[λ1, λ2][H1, H2, H3, n][u]; T2 = tau2[λ1, λ2][H1, H2, H3, n][u]; T3 = tau3[λ1, λ2][H1, H2, H3, n][u];
   QD = qdetT[λ1, λ2][u];
   N1 = ν1[λ1, λ2][u]; N1m = ν1[λ1, λ2][u - h]; N1mm = ν1[λ1, λ2][u - 2 h]; N1mmm = ν1[λ1, λ2][u - 3 h];
   sols = Table[
     Module[{M = Ms[[k]], zk = z[k], cols, mat, ns, vec, last, maxdeg},
       maxdeg = M + 8;
       cols = Table[
          PadRight[CoefficientList[Expand[
             zk (u + h)^j N1 N1m N1mm N1mmm
              - T1 u^j N1m N1mm N1mmm
              + T2/zk (u - h)^j N1mm N1mmm
              - T3/zk^2 (u - 2 h)^j N1mmm
              + QD/zk^3 (u - 3 h)^j], u], maxdeg + 1],
          {j, 0, M}];
       mat = Transpose[cols];
       ns = NullSpace[Chop[mat, 10^-12], Tolerance -> 10^-10];
       If[Length[ns] != 1, Message[QSolve4::degnull, H1, H2, H3, n, k, Length[ns]]; Abort[]];
       vec = ns[[1]]; last = vec[[-1]];
       If[Abs[last] < 10^-13 Max[Abs[vec]], Message[QSolve4::lastzero, H1, H2, H3, n, k]; Abort[]];
       vec/last],
     {k, 1, 4}];
   sols
];
```

- [ ] **Step 3: Insert the Qfun accessor cell** (from cell 36; k up to 4):

```wolfram
(*Q accessor: Q_k(u)=z[k]^(u/h) q_k(u), q_k the monic coeff list from QSolve4.*)
ClearAll[Qfun];
Qfun[λ1_, λ2_][H1_, H2_, H3_, n_][k_][u_] := z[k]^(u/h) qpoly[QSolve4[λ1, λ2][H1, H2, H3, n][[k]], u];
```

- [ ] **Step 4: Insert and run the single-state TQ check** (from cell 38; pick a real key from the eigensystem rather than guessing, k=1..4). This is the test:

```wolfram
(*single-state TQ-residual check on a non-symmetric rep: each Q_k must be annihilated by BaxOp
  at several generic numeric u. State = the first key of TauEigensystem4[{2,1,1,0},{1,0,0,0}].*)
Module[{λ1 = {2, 1, 1, 0}, λ2 = {1, 0, 0, 0}, key, H1, H2, H3, n, res},
  key = First[Keys[TauEigensystem4[λ1, λ2]]];
  {H1, H2, H3, n} = key;
  res = Table[
     Max@Table[
        Abs@Chop[N[BaxOp[λ1, λ2][H1, H2, H3, n][Function[uu, Qfun[λ1, λ2][H1, H2, H3, n][k][uu]]][u0], 30]],
        {u0, {2/3, 9/5, 4, 23/7}}],
     {k, 1, 4}];
  If[Max[res] > 10^-10, Print["TQ RESIDUAL FAIL: ", res]; Abort[]];
  {"state" -> key, "perQresidual" -> res}
]
```

**Expected output:** `{"state" -> {h1,h2,h3,n}, "perQresidual" -> {r1,r2,r3,r4}}` with all `r_k < 10^-10`. No `Abort`.

- [ ] **Step 5: Save and commit** (controller):

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb
git commit -m "su4 L2: QSolve4 (four Q-functions) + Qfun + single-state TQ residual check"
```

---

### Task 10: Full-sweep verification + fresh-kernel acceptance run

**Files:**
- Modify: `Baxter_L2_XXX_SU4.wb`

**Interfaces:**
- Consumes: everything (Tasks 1–9).
- Produces: the acceptance report (per-rep `nstates`, `worstTQ`, `worstCasoratian`).

- [ ] **Step 1: Insert the full-sweep verification cell** (from cell 39; 5 gl(4) reps, 4-tuple keys, k=1..4; Casoratian is a 5×4 shift matrix with `A_j = minor_{j+1}/minor_1`):

```wolfram
(*full sweep verification: over the 5 gl(4) rep pairs and every {H1,H2,H3,n} state, check
  (a) per-Q_k TQ residual and (b) the Casoratian reconstruction — the four Q_k span the 4th-order
  TQ solution space, so A1..A4 rebuilt from 4x4 minors over the five shift rows {u+h,u,u-h,u-2h,u-3h}
  must match tau1/tau2/tau3/qdetT over the appropriate ν1 products. Hard-asserts <10^-9.*)
Module[{sweep = {{{1,0,0,0},{1,0,0,0}},{{1,0,0,0},{1,1,1,0}},{{1,1,0,0},{1,0,0,0}},{{1,1,0,0},{1,1,0,0}},{{2,1,1,0},{1,0,0,0}}}, report},
  report = Table[
    Module[{λ1 = rep[[1]], λ2 = rep[[2]], es, worstTQ = 0, worstCas = 0},
      es = TauEigensystem4[λ1, λ2];
      KeyValueMap[
       Function[{key, rec},
        Module[{H1 = key[[1]], H2 = key[[2]], H3 = key[[3]], n = key[[4]], Qs, u0, rows, minors,
          A1rec, A2rec, A3rec, A4rec, A1, A2, A3, A4, N1, N1m, N1mm, N1mmm},
         Do[worstTQ = Max[worstTQ,
             Abs@Chop[N[BaxOp[λ1, λ2][H1, H2, H3, n][
                 Function[uu, Qfun[λ1, λ2][H1, H2, H3, n][k][uu]]][u0], 30]]],
           {k, 1, 4}, {u0, {5/4, 11/3}}];
         u0 = 13/5;
         Qs[uu_] := Table[Qfun[λ1, λ2][H1, H2, H3, n][k][uu], {k, 1, 4}];
         rows = {Qs[u0 + h], Qs[u0], Qs[u0 - h], Qs[u0 - 2 h], Qs[u0 - 3 h]}; (* 5 x 4 *)
         minors = Table[Det[Delete[rows, i]], {i, 1, 5}];
         A1rec = minors[[2]]/minors[[1]]; A2rec = minors[[3]]/minors[[1]];
         A3rec = minors[[4]]/minors[[1]]; A4rec = minors[[5]]/minors[[1]];
         N1 = ν1[λ1, λ2][u0]; N1m = ν1[λ1, λ2][u0 - h]; N1mm = ν1[λ1, λ2][u0 - 2 h]; N1mmm = ν1[λ1, λ2][u0 - 3 h];
         A1 = tau1[λ1, λ2][H1, H2, H3, n][u0]/N1;
         A2 = tau2[λ1, λ2][H1, H2, H3, n][u0]/(N1 N1m);
         A3 = tau3[λ1, λ2][H1, H2, H3, n][u0]/(N1 N1m N1mm);
         A4 = qdetT[λ1, λ2][u0]/(N1 N1m N1mm N1mmm);
         worstCas = Max[worstCas,
            Abs@Chop[N[A1rec - A1, 30]], Abs@Chop[N[A2rec - A2, 30]],
            Abs@Chop[N[A3rec - A3, 30]], Abs@Chop[N[A4rec - A4, 30]]]]],
       es];
      {rep, "nstates" -> Length[es], "worstTQ" -> worstTQ, "worstCasoratian" -> worstCas}],
    {rep, sweep}];
  If[Max[#[[3, 2]] & /@ report, #[[4, 2]] & /@ report] > 10^-9,
    Print["SWEEP FAIL: ", report]; Abort[]];
  Column[report]
]
```

**Expected output:** a `Column` of 5 rows, one per rep, each `{rep, "nstates"->d, "worstTQ"->0., "worstCasoratian"->0.}` with `nstates` = `{16,16,24,36,60}` respectively and every `worstTQ`/`worstCasoratian` `< 10^-9`. No `Abort`.

- [ ] **Step 2: Fresh-kernel top-to-bottom acceptance run.** Via `wolfbook_kernelControl` restart the kernel, then `runCell` every cell in order from the first to this one. Confirm no cell throws `Abort` and the final sweep cell reproduces the Step-1 output. (This catches stale-memoization / cell-ordering bugs the incremental build could hide.)

**Expected:** clean top-to-bottom run, all hard-assert cells pass in-kernel, final sweep `worstTQ = worstCasoratian = 0` across all 5 reps.

- [ ] **Step 3: Save and commit** (controller):

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb
git commit -m "su4 L2: full-sweep TQ+Casoratian verification (5 reps) + fresh-kernel acceptance run"
```

- [ ] **Step 4: Update `HANDOFF.md`** with the completed state (5 rep pairs verified, worstTQ=worstCasoratian=0, notebook `Baxter_L2_XXX_SU4.wb`), and note the `z→1` untwisted limit as the next thread. Commit.

---

## Notes for the executor

- **Faithful port, not reinvention.** Each cell's gl(3) original is in `su3_V2_general.wb` (cells 0–39). Where a step says "from cell N", diff against that original — the genuine gl(4) changes are exactly the ones written out above; everything else is index/term extension.
- **One notebook, appended in order.** Do not create scratch `.wl` files. Insert with `wolfbook_insertCells`, run with `wolfbook_runCell`, read output with `wolfbook_getCellOutput`, save with `wolfbook_fileOps`.
- **Memoization discipline.** `TauEigensystem4`, `alphaCoeffs`, `Mdeg`, `QSolve4`, `qmin*` all memoize. If any upstream definition is edited after they've been called, re-run the edited cell and clear/re-run these (or restart the kernel) before trusting downstream checks.
- **Blank output boxes.** If a check cell shows a blank `Out[]` after insert/edit, re-run it individually and read the value via `wolfbook_getCellOutput` / `evaluateExpression` — do not treat blank as pass.
- **The only expected hard place** is Task 8 Mdeg (the derived master formula) and Task 9 QSolve4 (the cleared-denominator columns). If Task 10 shows nonzero `worstTQ` but Task 7 eigen-residuals were 0, suspect a sign/shift error in the QSolve4 column expression or the Mdeg denominator before anything upstream.
```
