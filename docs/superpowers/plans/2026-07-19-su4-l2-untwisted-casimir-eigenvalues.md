# gl(4) L=2 Untwisted Transfer-Matrix Eigenvalues via Casimirs — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** In a new section of `Baxter_L2_XXX_SU4.wb`, appended after the existing twisted acceptance
run, find closed-form eigenvalues `Λ_a(u)` of the untwisted (`z[k]=1`) transfer matrices for
`λ1=[S,S,0,0]`, `λ2=[P,P,0,0]`, expressed via su(4) Casimir invariants of the irreps in
`λ1⊗λ2 = ⊕ N_μ μ`.

**Architecture:** Clear the z-dependent memoized caches and set `z[k]=1` (a section-boundary cell);
build the quadratic Casimir `C2op` and a Casimir-resolved eigensystem generalizing the twisted case's
proven "generic linear combination, diagonalize once" trick; extract `Λ_a(u)` per resolved state;
survey for multiplicity; fit a Casimir-based closed form on the multiplicity-free data and report
precisely.

**Tech Stack:** Wolfram Language in the same `.wb` notebook, driven live via the wolfbook MCP tools
(same as the twisted-case build). No standalone `.wl` scripts.

## Global Constraints

- **File:** all new cells appended to `Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb` after
  the existing 34 cells (twisted acceptance run, commit `fcd1349`). Do not modify cells 1–34.
- **Section boundary:** `ClearAll[qmin, qming, tg, tgCoeffs]; z[1]=1; z[2]=1; z[3]=1; z[4]=1;` — the
  only z-dependent *memoized* definitions that need clearing. `T,Tg,ν,ν1,qdetT,H1op,H2op,H3op` are
  not memoized and recompute correctly automatically; do not touch them. **`ClearAll` removes the
  cleared functions' rules entirely, not just their cache** — `qmin,qming,tg,tgCoeffs` must be
  re-defined (identical source, re-inserted) immediately after the boundary cell, or they stay
  undefined (found executing Task 1; see that task for the exact restoration cells).
- **Exact arithmetic simplification (new, not present in the twisted case):** with `z[k]=1` (exact
  integer) and `θ1=1/3,θ2=1/7,h=1` (exact rationals) and `JJ` built from rational GT-pattern
  arithmetic (no square roots in this notebook's convention), every quantity in this section is
  **exact rational** before any `N[...]` conversion. This sidesteps the twisted case's whole
  precision-cascade story (no irrational anchor whose digits need to propagate) — `N[...,30]` in the
  eigensystem step is a generous, not a load-bearing, precision choice here.
- **`qdetT` is central:** it is already a scalar function (same for every state, not per-state) —
  `Λ_4(u) = qdetT[λ1,λ2][u]` requires no per-eigenvector extraction at all.
- **Sweep:** `(S,P) ∈ {(1,1),(2,1),(1,2),(2,2)}`. Dimensions: `dim[{1,1,0,0}]=6`,
  `dim[{2,2,0,0}]=20`, giving tensor dims `36, 120, 120, 400`. `(2,2)` is the most expensive; if it
  proves too slow, descope it and report the finding on the other three (documented per-task, not a
  silent omission).
- **Verification convention:** hard-assert with `If[worst > tol, Print[...]; Abort[]]` where a closed
  form is claimed to hold; explicit, precise reporting (which sector, which coefficient, by how much)
  where it does not — never silently drop a failing case.
- **This is genuine research, not a mechanical port.** Two specific steps (Task 3's top-state formula,
  Task 5's Casimir-affine fit) are **conjectures being tested**, not asserted facts — treat a mismatch
  as a reportable finding, not a bug to force-fix.
- **Gotchas (carried from the twisted-case build):** `;` between statements in multi-statement cells;
  after fixing a memoized `f[x_]:=f[x]=…`, re-run the defining cell and everything downstream; save
  the notebook explicitly before committing; do not trust blank `Out[]` boxes — re-run.
- **Git:** controller commits between tasks (same as the twisted-case build).

---

### Task 1: Section boundary + quadratic Casimir + commutation sanity check

**Files:**
- Modify: `Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb` (append cells after cell 34)

**Interfaces:**
- Consumes: `qmin,qming,tg,tgCoeffs,z,χ1..χ4,EE,H1op,H2op,H3op,dim` (cells 1–34).
- Produces: `C2op[λ1,λ2]` (quadratic Casimir matrix on the two-site space, exact rational).

- [ ] **Step 1: Insert the section-boundary cell**:

```wolfram
(*Section boundary: switch to the untwisted (z=1) regime for the Casimir-eigenvalue exploration
  below. qmin/qming/tg/tgCoeffs are memoized on (λ1,λ2) but read z live inside their definitions —
  without this ClearAll they would silently return stale twisted-z values. T/Tg/ν/ν1/qdetT/H1op/
  H2op/H3op are NOT memoized and need no clearing; they recompute correctly automatically.*)
ClearAll[qmin, qming, tg, tgCoeffs];
z[1] = 1; z[2] = 1; z[3] = 1; z[4] = 1;
{χ1, χ2, χ3, χ4}
```

Run via `wolfbook_runCell`. **Expected output: `{4, 6, 4, 1}`** — the elementary symmetric functions
of four 1's, confirming the reassignment took hold and `χ` recomputed live.

**Gotcha found executing this task (not knowable from the twisted-case build alone): `ClearAll[f]`
wipes `f`'s definitions entirely, not just its memoized cache.** After `ClearAll[qmin,qming,tg,
tgCoeffs]`, those symbols have no rule at all — calling `tgCoeffs[...]` later just returns the
unevaluated symbolic expression (`Set::shape` error when destructured), it does NOT recompute at the
new `z=1`. **Fix:** immediately after the boundary cell, re-insert the *original defining code* for
`qmin`, `qming`, `tg`, `tgCoeffs` (identical source to the twisted-case build's Task 4/5 cells) so the
rules exist again — this naturally reads the freshly-assigned `z=1` on first call, since the memo
cache was correctly wiped. Two extra cells, inserted right after the boundary cell:

```wolfram
(*Restore the qmin/qming/tg/tgCoeffs definitions wiped by the ClearAll above -- ClearAll removes a
  function's rules entirely, not just its memoized cache, so these MUST be redefined (not merely
  re-called) here. Identical source to the twisted-case build; now reads the freshly-assigned z=1.*)
qmin[λ1_, λ2_][II_List, JJ_List][u_] := qmin[λ1, λ2][II, JJ][u] = Sum[Signature[σ] Dot @@ Table[T[λ1, λ2][II[[σ[[k]]]], JJ[[k]]][u - (k - 1) h], {k, Length[II]}], {σ, Permutations[Range[Length[II]]]}];
qming[λ1_, λ2_][II_List, JJ_List][u_] := qming[λ1, λ2][II, JJ][u] = Sum[Signature[σ] Dot @@ Table[Tg[λ1, λ2][II[[σ[[k]]]], JJ[[k]]][u - (k - 1) h], {k, Length[II]}], {σ, Permutations[Range[Length[II]]]}];
tg[λ1_,λ2_][a_,1][u_]:= Total[qming[λ1, λ2][#, #][u] & /@ Subsets[{1, 2, 3, 4}, {a}]]
```

```wolfram
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

- [ ] **Step 2: Insert the quadratic Casimir cell**:

```wolfram
(*Quadratic Casimir of gl(4) on the two-site space: C2 = Sum_ij EE[i,j].EE[j,i]. Exact rational
  (EE is built from JJ, which is rational-valued in this notebook's GT-pattern convention).*)
ClearAll[C2op];
C2op[λ1_, λ2_] := C2op[λ1, λ2] = Sum[EE[λ1, λ2][i, j] . EE[λ1, λ2][j, i], {i, 4}, {j, 4}];
```

- [ ] **Step 3: Insert and run the commutation sanity check** (this is the test — confirms `z=1`
  actually restored full su(4) symmetry, which the twisted case's `z≠1` does NOT have):

```wolfram
(*sanity check: at z=1 the untwisted transfer-matrix coefficients and Cartan operators must commute
  with the quadratic Casimir C2op, confirming full su(4) symmetry restoration (this would FAIL for
  z!=1). Tested on {1,1,0,0}⊗{1,1,0,0} (6⊗6=36, the smallest nontrivial rep pair). tgCoeffs
  recomputes fresh here (cleared in Step 1), now reading z=1.*)
Module[{λ1={1,1,0,0},λ2={1,1,0,0},c1,c2,c3,comm},
  {c1,c2,c3}=tgCoeffs[λ1,λ2];
  comm=Max[Table[Max[Abs[Flatten[m.C2op[λ1,λ2]-C2op[λ1,λ2].m]]],
     {m,Join[c1,c2,c3,{H1op[λ1,λ2],H2op[λ1,λ2],H3op[λ1,λ2]}]}]];
  If[comm>10^-8,Print["C2 COMMUTATION FAIL: ",comm];Abort[]];
  {"C2 commutes with tau-coeffs and Cartan at z=1, maxCommutator"->comm}
]
```

**Expected output:** `{"C2 commutes with tau-coeffs and Cartan at z=1, maxCommutator"->0}` — since
every input is exact rational, `comm` should literally be the exact integer/rational `0`, not merely
numerically small. No `Abort`.

- [ ] **Step 4: Save and commit** (controller):

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb
git commit -m "su4 L2 untwisted: section boundary (z=1) + quadratic Casimir C2op + commutation check"
```

---

### Task 2: Casimir-resolved eigensystem + Λ_a(u) accessor

**Files:**
- Modify: `Baxter_L2_XXX_SU4.wb`

**Interfaces:**
- Consumes: `tgCoeffs,H1op,H2op,H3op,C2op,dim,qdetT` (Task 1, cells 1–34).
- Produces: `TauEigensystemCasimir[λ1,λ2]` (Association `{H1,H2,H3,C2round,n} → <|t1c,t2c,t3c,C2,vec|>`),
  `LambdaU[λ1,λ2][H1,H2,H3,C2round,n][a][u]` (`a=1,2,3,4`).

- [ ] **Step 1: Insert the Casimir-resolved eigensystem cell** (generalizes `TauEigensystem4`'s
  "generic linear combination, diagonalize once" trick — same technique, but mixing in `C2op` instead
  of relying on twist-broken degeneracy; no `wtMult`-style predicted-count check exists here, so the
  routine only asserts the *total* resolved count equals `d`):

```wolfram
ClearAll[TauEigensystemCasimir];
TauEigensystemCasimir::rankdef = "Generic combination failed to resolve the `1`-dim eigenspace for (λ1,λ2)=(`2`,`3`).";
TauEigensystemCasimir::notgenuine = "Eigenvector #`1` for (λ1,λ2)=(`2`,`3`) failed the common-eigenvector residual check (residual=`4`).";
TauEigensystemCasimir::noninteger = "Weight `1` not near-integer for eigenvector #`2` at (λ1,λ2)=(`3`,`4`).";
TauEigensystemCasimir::degnull = "Total resolved states `1` != predicted dim `2` for (λ1,λ2)=(`3`,`4`).";
TauEigensystemCasimir[λ1_, λ2_] := TauEigensystemCasimir[λ1, λ2] = Module[
   {d, c1, c2, c3, mats, nontriv, rSets, rset, Kmix, vals, vecs, rayleigh, recs, byHW, table},
   d = dim[λ1] dim[λ2];
   {c1, c2, c3} = tgCoeffs[λ1, λ2];
   mats = Join[c1, c2, c3, {H1op[λ1, λ2], H2op[λ1, λ2], H3op[λ1, λ2], C2op[λ1, λ2]}];
   (* drop the scalar leading coeffs c1[[3]],c2[[5]],c3[[7]] (=chi_a Id, chi={4,6,4} at z=1) since
      they carry no state info; keep C2op since it's the whole point of this routine *)
   nontriv = Join[c1[[1 ;; 2]], c2[[1 ;; 4]], c3[[1 ;; 6]], {H1op[λ1, λ2], H2op[λ1, λ2], H3op[λ1, λ2], C2op[λ1, λ2]}];
   rSets = Table[N[Table[Prime[100 i + 7 j]^(1/3), {j, Length[nontriv]}], 30], {i, 1, 8}];
   rset = SelectFirst[rSets,
     Module[{ee, vv}, {ee, vv} = Eigensystem[N[#.nontriv, 30]];
        MatrixRank[vv] == d && Length[Union[Chop[ee, 10^-10]]] == d] &,
     Missing["NoGenericR"]];
   If[MissingQ[rset], Message[TauEigensystemCasimir::rankdef, d, λ1, λ2]; Abort[]];
   Kmix = rset . nontriv;
   {vals, vecs} = Eigensystem[N[Kmix, 30]];
   rayleigh[m_, v_] := (v . (m . v))/(v . v);
   recs = Table[
     Module[{v = vecs[[i]], res, H1v, H2v, H3v, C2v},
       res = Max[Abs[#.v - rayleigh[#, v] v]] & /@ mats // Max;
       If[res > 10^-8, Message[TauEigensystemCasimir::notgenuine, i, λ1, λ2, res]; Abort[]];
       H1v = rayleigh[H1op[λ1, λ2], v]; H2v = rayleigh[H2op[λ1, λ2], v]; H3v = rayleigh[H3op[λ1, λ2], v];
       C2v = rayleigh[C2op[λ1, λ2], v];
       If[Abs[H1v - Round[Re[H1v]]] > 10^-6, Message[TauEigensystemCasimir::noninteger, H1v, i, λ1, λ2]; Abort[]];
       If[Abs[H2v - Round[Re[H2v]]] > 10^-6, Message[TauEigensystemCasimir::noninteger, H2v, i, λ1, λ2]; Abort[]];
       If[Abs[H3v - Round[Re[H3v]]] > 10^-6, Message[TauEigensystemCasimir::noninteger, H3v, i, λ1, λ2]; Abort[]];
       <|"H1" -> Round[Re[H1v]], "H2" -> Round[Re[H2v]], "H3" -> Round[Re[H3v]],
         "C2round" -> Round[Re[C2v], 10^-6], "C2" -> Re[C2v],
         "t1c" -> rayleigh[#, v] & /@ c1, "t2c" -> rayleigh[#, v] & /@ c2, "t3c" -> rayleigh[#, v] & /@ c3,
         "vec" -> v|>],
     {i, 1, d}];
   byHW = GroupBy[recs, {#["H1"], #["H2"], #["H3"], #["C2round"]} &];
   table = Association @ Flatten @ KeyValueMap[
      Function[{key, group},
        MapIndexed[Function[{rec, idx}, Append[key, First[idx] - 1] -> rec], group]],
      byHW];
   If[Length[table] != d, Message[TauEigensystemCasimir::degnull, Length[table], d, λ1, λ2]; Abort[]];
   table
];
```

- [ ] **Step 2: Insert the `Λ_a(u)` accessor cell**:

```wolfram
(*Lambda accessor: a=1,2,3 via Rayleigh-quotient coefficients from TauEigensystemCasimir;
  a=4 is qdetT itself (the quantum determinant is central -- same scalar function on EVERY
  state, no per-eigenvector extraction needed).*)
ClearAll[LambdaU];
LambdaU[λ1_, λ2_][H1_, H2_, H3_, C2round_, n_][a_][u_] := Module[
   {rec = TauEigensystemCasimir[λ1, λ2][Key[{H1, H2, H3, C2round, n}]]},
   Which[
     a == 1, rec["t1c"] . u^Range[0, 2],
     a == 2, rec["t2c"] . u^Range[0, 4],
     a == 3, rec["t3c"] . u^Range[0, 6],
     a == 4, qdetT[λ1, λ2][u]
   ]
];
```

- [ ] **Step 3: Insert and run the eigen-residual / consistency check on `(S,P)=(1,1)`** (this is the
  test — mirrors the twisted case's `TauEigensystem4` eigen-residual check):

```wolfram
Module[{λ1 = {1, 1, 0, 0}, λ2 = {1, 1, 0, 0}, es, worst},
  es = TauEigensystemCasimir[λ1, λ2];
  worst = Max@KeyValueMap[
     Function[{k, rec},
       Module[{H1 = k[[1]], H2 = k[[2]], H3 = k[[3]], C2r = k[[4]], n = k[[5]], v = rec["vec"], uu = 17/5},
         Max[
          Max[Abs[(tg[λ1, λ2][1, 1][uu] . v) - LambdaU[λ1, λ2][H1, H2, H3, C2r, n][1][uu] v]],
          Max[Abs[(tg[λ1, λ2][2, 1][uu] . v) - LambdaU[λ1, λ2][H1, H2, H3, C2r, n][2][uu] v]],
          Max[Abs[(tg[λ1, λ2][3, 1][uu] . v) - LambdaU[λ1, λ2][H1, H2, H3, C2r, n][3][uu] v]]]]],
     es];
  If[worst > 10^-8, Print["EIGEN RESIDUAL FAIL: ", worst]; Abort[]];
  {"nstates" -> Length[es], "worstResidual@u=17/5" -> worst}
]
```

**Expected output:** `{"nstates" -> 36, "worstResidual@u=17/5" -> w}` with `w` extremely small
(effectively 0, given the exact-rational inputs). `nstates` must equal
`dim[{1,1,0,0}]·dim[{1,1,0,0}] = 6·6 = 36`. No `Abort`.

**If `TauEigensystemCasimir::rankdef`, `::notgenuine`, or `::degnull` fires:** this is genuinely new,
never-before-run code — debug interactively in the notebook (inspect `rSets`/`rset` selection,
check `MatrixRank[vecs]`, verify `C2op` is actually rational/well-formed) rather than guessing a fix
blind. Report what was found before proceeding.

- [ ] **Step 4: Save and commit** (controller):

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb
git commit -m "su4 L2 untwisted: TauEigensystemCasimir + LambdaU accessor + eigen-residual check"
```

---

### Task 3: Top-state closed-form conjecture test

**Files:**
- Modify: `Baxter_L2_XXX_SU4.wb`

**Interfaces:**
- Consumes: `TauEigensystemCasimir, LambdaU, ν` (Tasks 1–2).

**Context (why this state is special and safe to identify automatically):** the "top" component
`μ = λ1+λ2` — realized by the literal highest-weight vector of the tensor product, weight
`(S+P,S+P,0,0)`, i.e. `(H1,H2,H3)=(0,S+P,0)` — always occurs with multiplicity **exactly 1** in any
tensor product decomposition (a basic representation-theory fact, not specific to this rep family).
This makes it the safest place to test a closed-form conjecture first: no multiplicity ambiguity to
resolve. **Conjecture being tested** (not asserted as fact): since `T[i,i](u)|hw⟩ = ν[i](u)|hw⟩`
already holds on the literal highest-weight vector (established in the twisted-case build, cell 11,
and true independent of `z`), and `Λ_1(u) = tg[1,1](u)` is a sum over rank-1 principal quantum minors
which reduce to `T[i,i](u)` on a state where the off-diagonal `T[i,j](u)` matrix elements vanish, the
conjecture is `Λ_1(u)|_top = ν[λ1,λ2][1][u] + ν[λ1,λ2][2][u] + ν[λ1,λ2][3][u] + ν[λ1,λ2][4][u]`.

- [ ] **Step 1: Insert and run the top-state identification + conjecture test** (this is the test):

```wolfram
(*Test the top-state (mu=lambda1+lambda2) closed-form conjecture: Lambda_1(u) = Sum_i nu[i](u),
  with NO shifts (unlike qdetT's shifted product) since a=1 principal minors are literally T[i,i](u)
  with zero row-shift in the quantum-minor definition. mu=lambda1+lambda2 has weight (S+P,S+P,0,0),
  i.e. (H1,H2,H3)=(0,S+P,0), and occurs with multiplicity exactly 1 (general fact) -- so there is
  exactly one key in TauEigensystemCasimir's output with that (H1,H2,H3).*)
Module[{λ1 = {1, 1, 0, 0}, λ2 = {1, 1, 0, 0}, S = 1, P = 1, es, topKeys, key, H1, H2, H3, C2r, n, conj, actual, worst},
  es = TauEigensystemCasimir[λ1, λ2];
  topKeys = Select[Keys[es], {#[[1]], #[[2]], #[[3]]} == {0, S + P, 0} &];
  If[Length[topKeys] != 1, Print["TOP STATE NOT UNIQUE: ", topKeys]; Abort[]];
  key = First[topKeys];
  {H1, H2, H3, C2r, n} = key;
  conj[u_] := ν[λ1, λ2][1][u] + ν[λ1, λ2][2][u] + ν[λ1, λ2][3][u] + ν[λ1, λ2][4][u];
  worst = Max@Table[Abs[N[LambdaU[λ1, λ2][H1, H2, H3, C2r, n][1][u0] - conj[u0], 30]], {u0, {2/3, 9/5, 4, 23/7}}];
  {"topKey" -> key, "conjectureMatches" -> (worst < 10^-8), "worst" -> worst}
]
```

**Expected outcome (report exactly what happens, do not force a pass):** either
`"conjectureMatches" -> True` with `worst` effectively 0 (confirms the conjecture — proceed to Task 4
using this as the anchor/reference point for the Casimir fit), or `"conjectureMatches" -> False` with
a nonzero `worst` (the conjecture is wrong as stated — report the actual `LambdaU[...][1][u]` value
computed vs. the conjectured one, and flag this for the controller before Task 4 proceeds, since it
changes what baseline Task 5's fit anchors against).

- [ ] **Step 2: Save and commit** (controller) — commit regardless of the outcome, with a message
  reflecting what was actually found:

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb
git commit -m "su4 L2 untwisted: top-state (mu=lambda1+lambda2) closed-form conjecture test"
```

---

### Task 4: Multiplicity survey across the sweep

**Files:**
- Modify: `Baxter_L2_XXX_SU4.wb`

**Interfaces:**
- Consumes: `TauEigensystemCasimir, LambdaU` (Tasks 1–2).
- Produces: `multData` (in-notebook variable: per-`(S,P)` report of which `(H1,H2,H3,C2round)`
  sectors have >1 states and whether those states share `Λ_a(u)` or genuinely differ — i.e. whether
  `N_μ>1` there).

- [ ] **Step 1: Insert and run the multiplicity survey** (this is the test — an honest report, not a
  pass/fail hard-assert, since the outcome is genuinely unknown in advance):

```wolfram
(*For each (S,P) in the sweep, run TauEigensystemCasimir and check every (H1,H2,H3,C2round) sector
  with >1 resolved state: if all states in that sector give the SAME Lambda_a(u) (a=1,2,3), it's a
  harmless numerical near-degeneracy in C2round's rounding, not genuine multiplicity N_mu>1. If they
  give GENUINELY DIFFERENT Lambda_a(u), that sector has real multiplicity -- no pure-Casimir formula
  can apply there (reported, not treated as failure).*)
multData = Table[
  Module[{λ1 = {S, S, 0, 0}, λ2 = {P, P, 0, 0}, es, byKey3, multSectors, genuineMult},
    es = TauEigensystemCasimir[λ1, λ2];
    byKey3 = GroupBy[Keys[es], {#[[1]], #[[2]], #[[3]], #[[4]]} &];
    multSectors = Select[byKey3, Length[#] > 1 &];
    genuineMult = Select[multSectors, Function[keys,
       Module[{lam1s = LambdaU[λ1, λ2][keys[[1, 1]], keys[[1, 2]], keys[[1, 3]], keys[[1, 4]], keys[[1, 5]]][1][17/5]},
         Max[Abs[LambdaU[λ1, λ2][#[[1]], #[[2]], #[[3]], #[[4]], #[[5]]][1][17/5] - lam1s] & /@ keys]] > 10^-6]];
    {S, P, "dim" -> dim[λ1] dim[λ2], "nstates" -> Length[es], "sectorsWithMult" -> Length[multSectors],
     "genuineMultSectors" -> Length[genuineMult], "multFree" -> (Length[genuineMult] == 0)}],
  {SP, {{1, 1}, {2, 1}, {1, 2}, {2, 2}}}, {S = SP[[1]], P = SP[[2]]}];
Column[multData]
```

**Expected:** no `Abort` (this cell reports, it doesn't hard-assert). Read and record the actual
`"multFree"` value for each `(S,P)`. **If `(2,2)` (dim 400) is too slow to complete in reasonable
time, descope it**: rerun the `Table` over just `{{1,1},{2,1},{1,2}}` and note in the commit message
that `(2,2)` was not attempted due to computational cost.

- [ ] **Step 2: Save and commit** (controller), message reflecting the actual multiplicity findings:

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb
git commit -m "su4 L2 untwisted: multiplicity survey across (S,P) sweep"
```

---

### Task 5: Casimir-affine fit on the multiplicity-free data

**Files:**
- Modify: `Baxter_L2_XXX_SU4.wb`

**Interfaces:**
- Consumes: `TauEigensystemCasimir, LambdaU, multData, C2op` (Tasks 1–4).

**Candidate being tested (not asserted as fact):** generalizing the su(2) Heisenberg-chain relation
`H = P − 1 ∝ C2(total) − C2(site1) − C2(site2)`, the conjecture is that each polynomial coefficient of
`Λ_1(u)` (as a function of `u`, degree 2) is an **affine function of `C2(μ)` alone** — i.e. for fixed
`(S,P)`, varying only over the different `μ`'s appearing (the different `(H1,H2,H3,C2)` sectors found
in Task 4 that are multiplicity-free), `coeff_j(Λ_1) = A_j + B_j·C2(μ)` for constants `A_j,B_j`
(`j=0,1,2`, the `u^j` coefficient), fit from **two** anchor points (any two mult-free `μ`'s, e.g. the
top state from Task 3 plus one other) and verified/hard-asserted against **every remaining**
mult-free `μ` in the sweep.

- [ ] **Step 1: Insert and run the data collection + affine fit** (this is the test):

```wolfram
(*Collect (mu-label via (H1,H2,H3,C2), Lambda_1(u) coeffs) rows for every multiplicity-free sector
  across the sweep (using multData's multFree flag; skip sectors flagged genuinely-multiplicity in
  Task 4). Fit each u^j coefficient of Lambda_1 as an affine function of the C2 eigenvalue, using the
  first two distinct-C2 rows as anchors, then hard-assert against every remaining row.*)
casData = Flatten[Table[
   Module[{λ1 = {S, S, 0, 0}, λ2 = {P, P, 0, 0}, es, byKey3, genuineMultKeys},
     es = TauEigensystemCasimir[λ1, λ2];
     byKey3 = GroupBy[Keys[es], {#[[1]], #[[2]], #[[3]], #[[4]]} &];
     genuineMultKeys = Flatten[Values[Select[byKey3, Length[#] > 1 &]], 1]; (* skip these -- ambiguous *)
     KeyValueMap[
       Function[{key, rec}, If[MemberQ[genuineMultKeys, key], Nothing,
         {S, P, key, rec["C2"], CoefficientList[LambdaU[λ1, λ2][key[[1]], key[[2]], key[[3]], key[[4]], key[[5]]][1][uu], uu]}]],
       es]],
   {SP, {{1, 1}, {2, 1}, {1, 2}}}, {S = SP[[1]], P = SP[[2]]}], 2];
(*fit u^j coefficients as affine in C2 from the first two rows with distinct C2, verify on the rest*)
Module[{distinctRows, anchor1, anchor2, coeffFit, worst},
  distinctRows = DeleteDuplicatesBy[casData, #[[4]] &];
  If[Length[distinctRows] < 2, Print["NOT ENOUGH DISTINCT C2 VALUES TO FIT"]; Abort[]];
  {anchor1, anchor2} = distinctRows[[1 ;; 2]];
  coeffFit[j_] := Module[{c1v = anchor1[[5, j + 1]], c2v = anchor2[[5, j + 1]], x1 = anchor1[[4]], x2 = anchor2[[4]], A, B},
     B = (c1v - c2v)/(x1 - x2); A = c1v - B x1; {A, B}];
  worst = Max@Table[
     Module[{A, B, pred},
       {A, B} = coeffFit[j];
       Max[Abs[(A + B #[[4]]) - #[[5, j + 1]]] & /@ casData]],
     {j, 0, 2}];
  {"nDataRows" -> Length[casData], "nDistinctC2" -> Length[distinctRows],
   "coeffFits(A,B)_j=0,1,2" -> Table[coeffFit[j], {j, 0, 2}], "worstDeviation" -> worst}
]
```

**Expected outcome (report precisely, do not force a pass):** if `"worstDeviation"` is effectively 0,
the Casimir-affine conjecture holds exactly across the multiplicity-free sweep — record the fitted
`(A_j,B_j)` as the closed form. If `"worstDeviation"` is large for some `j`, report exactly which
coefficient order fails and by how much, and check (as a documented fallback, not required unless
Step 1 fails) whether a **quadratic** (rather than affine) function of `C2` fits instead, or whether
`C2` alone is insufficient and a second invariant (`H2` itself, or a cubic Casimir) is needed — using
the same "collect data, fit candidates in a named `Association`, report exact match or exact failure"
pattern already established in this notebook's twisted-case degree-formula work (cells 29–33).

- [ ] **Step 2: Save the notebook and commit** (controller), message reflecting the actual result:

```bash
cd /home/paul/git_projects/gaudin-CG
git add Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb
git commit -m "su4 L2 untwisted: Casimir-affine fit for Lambda_1(u) on multiplicity-free sweep data"
```

- [ ] **Step 3: Update `HANDOFF.md`** with the final state of this exploration — what closed form (if
  any) was verified, on what scope (which `(S,P)`, which `μ`'s), and what was found NOT to reduce to
  a pure-Casimir formula (multiplicity sectors from Task 4, any failing coefficients from Task 5).
  Commit.

---

## Notes for the executor

- **This plan is smaller and more exploratory than the twisted-case build.** Tasks 1–2 are
  low-risk (proven pattern, just a new operator family). Tasks 3 and 5 are genuinely testing
  conjectures — a "does not match" result is a valid, useful outcome, not something to iterate on
  blindly. If Task 5's simplest affine ansatz fails, stop and report to the controller with the exact
  data (which coefficient, which sectors, the actual deviation) rather than trying many further
  candidates unprompted — this matches the project's established "ask before ansatz-search loops"
  convention.
- **`TauEigensystemCasimir` is new, untested code.** Debug interactively via the wolfbook tools if it
  doesn't run cleanly on the first attempt (inspect `rSets` selection, `MatrixRank`, `C2op`
  structure) — this is expected for genuinely new research code, unlike the twisted-case port which
  had a validated gl(3) source to check against line-by-line.
- **Exact arithmetic is a feature here, not a risk.** Unlike the twisted case, there is no
  precision-cascade concern in this section (`z=1` is exact, no irrational anchor). If a residual
  check ever shows an unexpectedly large nonzero value, treat it as a genuine bug signal (not a
  precision artifact needing more digits) — the debugging lesson from the twisted case's Task 10 does
  **not** transfer here.
