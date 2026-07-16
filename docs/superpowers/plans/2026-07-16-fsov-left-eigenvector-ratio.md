# Left Eigenvector / <Psi|Psi> Ratio Implementation Plan

> **For agentic workers:** none of the standard superpowers execution skills apply here — this
> plan is executed entirely inside a live Wolfram notebook via the project's **wolfbook-builder**
> subagent (see `CLAUDE.md`), not via pytest/git-commit-per-task. Each task specifies notebook
> cells to insert (via `mcp__wolfbook__wolfbook_insertCells`), run (via
> `mcp__wolfbook__wolfbook_runCell`), and verify (via `mcp__wolfbook__wolfbook_getCellOutput` /
> `wolfbook_getKernelState`), against an explicit expected result. Steps use checkbox (`- [ ]`)
> syntax for tracking. There is no separate test file and no per-task git commit.

**Goal:** Add a new Section D to `Experiments/XXX_CG_L2_V2.wb` that computes left eigenvectors
of the companion-twisted transfer matrix (normalized so the last/bottom component is 1, matching
`Psi`'s convention), pairs each with its same-state right eigenvector to form `<Psi|Psi>`, and
reports the ratio `<Psi|Psi>/det[Q,Q]` for all 6 states of `(λ1,λ2)=(2,1)`.

**Architecture:** Three notebook-cell tasks appended after Section C: (1)
`TauEigensystemLeft`/`PsiLeft` — independent left-eigenvector construction matched to existing
`{M,n}` states via eigenvalue comparison, (2) a biorthogonality hard-assert (the correctness
check on the left-eigenvector construction itself), (3) `braket`/`normRatio` plus a final report
table.

**Tech Stack:** Wolfram Language, via wolfbook MCP tools operating on the live `.wb` notebook and
kernel.

## Global Constraints

- Design doc: `docs/superpowers/specs/2026-07-16-fsov-left-eigenvector-ratio-design.md` — this
  plan implements it in full; do not deviate from its scope (no sweep over other `(λ1,λ2)`, no
  predicted-ratio comparison, no edits to Section B or Section C).
- Test case is fixed: `(λ1,λ2)=(2,1)`, reusing `states21` and `G21` already computed in Section C
  (`{{0,0},{1,0},{1,1},{2,0},{2,1},{3,0}}`, verified numerically diagonal Gram matrix).
- `<Psi|Psi>` is the plain bilinear pairing (no complex conjugation) — `PsiLeft . Psi`, matching
  the notebook's existing non-Hermitian convention (`TauEigensystem`'s `t0=(v.C0.v)/(v.v)`).
- Hard-assert discipline, matching this notebook's existing style (`BaxterMatrix`,
  `TauEigensystem`, `QSolve`, Section C's `FSoVGram`): every structural/numeric check that could
  silently hide a bug aborts (`Abort[]`) with a `Message`, not a soft warning. If a check fails
  when actually run, **stop and investigate — do not weaken the tolerance or delete the check to
  force a pass.**
- Comment convention: `(* ... *)` cells explain *why*, not what.
- `wolfbook-builder` stays out of `paper/*.tex` and git commits/pushes (per `CLAUDE.md`) — this
  plan produces no git commits; saving the notebook (`wolfbook_fileOps`) is the only persistence
  step, done once at the end of Task 3.
- Memoized routines this plan adds (`TauEigensystemLeft`) should be re-run (along with anything
  downstream) if any upstream definition is edited later — same caching caveat as the rest of the
  notebook.
- Do not touch Section B (`TauEigensystem`, `Psi`, `τ`, `QSolve`, `Q1`) or Section C (`mu`, `Nn`,
  `nn`, `br`, `det`, `statesOf`, `GramMatrix`, `states21`, `G21`) — Section D only reads from
  them.

---

### Task 1: `TauEigensystemLeft` + `PsiLeft`

**Files:**
- Modify: `Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb` (append new cells after Section C's
  last cell — the verification `<|...|>` Association cell)

**Interfaces:**
- Consumes: `t2Coeffs[λ1,λ2]` (Section A, memoized `{C0,C1,C2}`), `τ[λ1,λ2][M,n][u]` (Section B),
  `dPred[M,λ1,λ2]` (Section B), `statesOf[λ1,λ2]` (Section C), `κ1,κ2` (cell 0 globals).
- Produces: `TauEigensystemLeft[λ1,λ2]` (memoized `Association` keyed by `{M,n} -> vec` directly
  — no extra `"PsiLeft"` wrapper key, since this cache stores only one field per state), `PsiLeft
  [λ1,λ2][M,n]` (accessor) — both consumed by Task 2 and Task 3.

- [ ] **Step 1: Insert the Section D header comment cell**

Insert as a new cell after Section C's last cell:

```wolfram
(*Section D — left eigenvectors of the companion-twisted transfer matrix, paired with the
  existing right eigenvectors (Psi, Section B) to form <Psi|Psi>, and its ratio against the
  FSoV scalar product det[Q,Q] (Section C). TauEigensystemLeft is an INDEPENDENT construction
  (does not edit TauEigensystem): it diagonalizes Transpose[C0+r C1] for its own generic r, then
  matches each left eigenvector to an existing {M,n} state by comparing eigenvalues (since
  t2(r).v=τ(r) v for the right eigenvector and C2=(κ1+κ2)Id exactly, the target combo-eigenvalue
  for state {M,n} is τ[λ1,λ2][M,n][r]-(κ1+κ2) r^2). Exploratory: no predicted formula for the
  ratio exists yet in paper/xxxCG.tex (still single-site only) — this reports values, it doesn't
  assert against one. See docs/superpowers/specs/2026-07-16-fsov-left-eigenvector-ratio-design.md.*)
```

Run the cell. Expected: comment-only cell, no output, no error.

- [ ] **Step 2: Insert and run the `TauEigensystemLeft` cell**

```wolfram
ClearAll[TauEigensystemLeft];
TauEigensystemLeft::rankdef = "Left combination C0+r C1 failed to resolve the `1`-dim eigenspace for (λ1,λ2)=(`2`,`3`) at every tried r in `4`.";
TauEigensystemLeft::nomatch = "No left eigenvalue matched target for state `1` at (λ1,λ2)=(`2`,`3`) within tolerance (closest distance `4`).";
TauEigensystemLeft::dupmatch = "Left eigenvector index `4` matched more than one state (latest: `1`) at (λ1,λ2)=(`2`,`3`) — combo eigenvalues are not resolving distinctly.";
TauEigensystemLeft::lastzero = "Matched left eigenvector for state `1` at (λ1,λ2)=(`2`,`3`) has last component ~0; cannot normalize.";
TauEigensystemLeft::countmismatch = "TauEigensystemLeft produced `3` states for (λ1,λ2)=(`1`,`2`), expected `4`.";
TauEigensystemLeft[λ1_, λ2_] := TauEigensystemLeft[λ1, λ2] = Module[
   {d, Cs, C0, C1, C2, rCandidatesL, r, valsL, vecsL, states, usedIdx, table},
   d = (λ1 + 1) (λ2 + 1);
   Cs = t2Coeffs[λ1, λ2]; {C0, C1, C2} = Cs;
   rCandidatesL = {17/6, N[Pi, 30] + 1, N[E, 30] + 1, 9/4, 19/6, N[Sqrt[3], 30]};
   r = SelectFirst[rCandidatesL,
     Module[{ee, vv}, {ee, vv} = Eigensystem[Transpose[C0 + # C1]];
       MatrixRank[vv] == d && Length[Union[Chop[ee, 10^-12]]] == d] &,
     Missing["NoGenericR"]];
   If[MissingQ[r], Message[TauEigensystemLeft::rankdef, d, λ1, λ2, rCandidatesL]; Abort[]];
   {valsL, vecsL} = Eigensystem[Transpose[N[C0 + r C1, 30]]];
   states = statesOf[λ1, λ2];
   usedIdx = {};
   table = Association @ Table[
      Module[{Mn = states[[k]], target, dists, idx, w, last},
        target = τ[λ1, λ2][Mn[[1]], Mn[[2]]][r] - (κ1 + κ2) r^2;
        dists = Abs[valsL - target];
        idx = First[Ordering[dists]];
        If[dists[[idx]] > 10^-10,
          Message[TauEigensystemLeft::nomatch, Mn, λ1, λ2, dists[[idx]]]; Abort[]];
        If[MemberQ[usedIdx, idx],
          Message[TauEigensystemLeft::dupmatch, Mn, λ1, λ2, idx]; Abort[]];
        AppendTo[usedIdx, idx];
        w = vecsL[[idx]];
        last = w[[-1]];
        If[Abs[last] < 10^-15 Max[Abs[w]],
          Message[TauEigensystemLeft::lastzero, Mn, λ1, λ2]; Abort[]];
        Mn -> w/last],
      {k, 1, Length[states]}];
   If[Length[table] != d,
     Message[TauEigensystemLeft::countmismatch, λ1, λ2, Length[table], d]; Abort[]];
   table
];
```

Expected: no output (all `:=`/`Set` definitions and `Message` string assignments), no error.

- [ ] **Step 3: Insert and run the `TauEigensystemLeft` smoke test cell**

```wolfram
(*smoke test: (λ1,λ2)=(2,1) must produce the same 6 keys as statesOf[2,1] (already verified in
  Section C), just via the independent left-eigenvector construction.*)
KeySort[TauEigensystemLeft[2, 1]] // Keys
```

Expected: `{{0,0},{1,0},{1,1},{2,0},{2,1},{3,0}}` (same as `states21`, just possibly reordered by
`KeySort`'s default ordering — confirm the *set* of keys matches; if it doesn't, or if
`TauEigensystemLeft::rankdef`/`::nomatch`/`::dupmatch`/`::countmismatch` fired, **stop** — do not
proceed to Step 4 until this matches exactly.

- [ ] **Step 4: Insert and run the `PsiLeft` accessor cell**

```wolfram
(*PsiLeft[λ1,λ2][M,n]: on-demand accessor over the memoized TauEigensystemLeft cache, same
  bounds/spurious-guard pattern as Psi (Section B).*)
ClearAll[PsiLeft];
PsiLeft::bnds = "Invalid indices (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): need M an integer in 0..λ1+λ2 and n an integer in 0..M.";
PsiLeft::spurious = "No genuine state at (M,n)=(`1`,`2`) for (λ1,λ2)=(`3`,`4`): weight sector M=`1` has only `5` genuine state(s) (n=0..`6`), so n=`2` does not exist.";
PsiLeft[λ1_, λ2_][M_, n_] := Module[{tab, mult},
  If[! (IntegerQ[M] && IntegerQ[n] && 0 <= M <= λ1 + λ2 && 0 <= n <= M),
    Message[PsiLeft::bnds, M, n, λ1, λ2]; Return[$Failed]];
  tab = TauEigensystemLeft[λ1, λ2];
  mult = dPred[M, λ1, λ2];
  If[n > mult - 1,
    Message[PsiLeft::spurious, M, n, λ1, λ2, mult, mult - 1]; Return[Missing["NoState", {λ1, λ2, M, n}]]];
  tab[{M, n}]];
(*smoke test: a valid state's left eigenvector, plus confirming its last component is exactly 1*)
{PsiLeft[2, 1][1, 0][[-1]], Length[PsiLeft[2, 1][1, 0]]}
```

Expected: `{1, 6}` (last component exactly `1`, vector length 6 matching `d=(2+1)(1+1)`).

- [ ] **Step 5: Save the notebook**

Use `mcp__wolfbook__wolfbook_fileOps` to save `XXX_CG_L2_V2.wb`.

---

### Task 2: Biorthogonality verification

**Files:**
- Modify: `Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb` (append new cells after Task 1's
  last cell)

**Interfaces:**
- Consumes: `PsiLeft[λ1,λ2][M,n]` (Task 1), `Psi[λ1,λ2][M,n]` (Section B), `states21` (Section
  C).
- Produces: printed verification report (nothing further consumes this programmatically, but it
  is the correctness gate before Task 3 trusts `PsiLeft`).

- [ ] **Step 1: Insert and run the biorthogonality check cell**

```wolfram
(*Biorthogonality: a genuine left/right eigenvector pair of a diagonalizable operator with
  simple spectrum must satisfy PsiLeft[M,n].Psi[M',n']=0 for {M,n}!={M',n'}. This is the
  correctness check on TauEigensystemLeft's construction itself (independent of Section C's
  det[Q,Q] machinery) — the natural companion to Section C's off-diagonal Gram check. Threshold
  10^-15 matches this notebook's existing hard-assert convention.*)
BiorthoCheck::fail = "Off-diagonal <PsiLeft|Psi> pairing does not vanish for (λ1,λ2)=(2,1): max|off-diag|=`1`.";
offBiortho = Flatten[Table[
   If[i != j, PsiLeft[2, 1] @@ states21[[i]] . (Psi[2, 1] @@ states21[[j]]), Nothing],
   {i, 1, 6}, {j, 1, 6}]];
maxOffBiortho = Max[Abs[N[offBiortho]]];
If[maxOffBiortho > 10^-15,
  Message[BiorthoCheck::fail, maxOffBiortho];
  Abort[]];
<|"max|<PsiLeft|Psi>| off-diagonal" -> maxOffBiortho,
  "digits vanished (-Log10)" -> N[-Log10[maxOffBiortho]]|>
```

Expected: no `Abort[]`, no `Message`; the printed `Association` shows `"max|<PsiLeft|Psi>|
off-diagonal"` at or below `~10^-15`. **If the assertion fires, stop and investigate** — per the
Global Constraints, do not raise the threshold; a real failure here means
`TauEigensystemLeft`'s matching (Task 1, Step 2) has a bug (e.g. a mismatched state or a
duplicate-index collision that slipped past the `dupmatch` guard).

- [ ] **Step 2: Save the notebook**

Use `mcp__wolfbook__wolfbook_fileOps` to save `XXX_CG_L2_V2.wb`.

---

### Task 3: `braket`/`normRatio` and final report

**Files:**
- Modify: `Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb` (append new cells after Task 2's
  last cell)

**Interfaces:**
- Consumes: `PsiLeft[λ1,λ2][M,n]` (Task 1), `Psi[λ1,λ2][M,n]` (Section B), `states21`, `G21`
  (Section C).
- Produces: printed report table (terminal deliverable of this plan — nothing further consumes
  this programmatically).

- [ ] **Step 1: Insert and run the `braket`/`normRatio` cell**

```wolfram
(*braket = <Psi|Psi>: plain bilinear pairing (no conjugation), matching this notebook's existing
  non-Hermitian convention. normRatio divides by G21's diagonal entry for the SAME state (G21 is
  Section C's already-verified Gram matrix over states21) — only valid for (λ1,λ2)=(2,1) since
  G21/states21 are fixed globals for that representation, matching this plan's fixed test case.*)
ClearAll[braket, normRatio];
braket[λ1_, λ2_][M_, n_] := PsiLeft[λ1, λ2][M, n] . Psi[λ1, λ2][M, n];
normRatio[λ1_, λ2_][M_, n_] := Module[{i = First[FirstPosition[states21, {M, n}]]},
  braket[λ1, λ2][M, n]/G21[[i, i]]];
(*smoke test*)
{braket[2, 1][1, 0], normRatio[2, 1][1, 0]}
```

Expected: no error; a pair of nonzero numbers (the second is the first divided by `G21`'s
`{1,0}`-state diagonal entry).

- [ ] **Step 2: Insert and run the final report cell**

```wolfram
(*Final report: braket, G21 diagonal entry, and normRatio for all 6 states of (λ1,λ2)=(2,1).
  Hard-assert every braket is nonzero (non-degenerate pairing) before reporting — a zero braket
  here would mean PsiLeft and Psi are somehow orthogonal even on the diagonal, which should be
  impossible given Task 2's biorthogonality check already passed.*)
BraketReport::degenerate = "braket vanished for state `1` at (λ1,λ2)=(2,1) (value `2`) — pairing may be degenerate.";
reportRows = Table[
   Module[{Mn = states21[[i]], b, g, ratio},
     b = braket[2, 1] @@ Mn;
     g = G21[[i, i]];
     ratio = normRatio[2, 1] @@ Mn;
     If[Abs[N[b]] < 10^-15,
       Message[BraketReport::degenerate, Mn, b];
       Abort[]];
     <|"M" -> Mn[[1]], "n" -> Mn[[2]], "braket" -> N[b, 10], "G21 diagonal" -> N[g, 10],
       "normRatio" -> N[ratio, 10]|>],
   {i, 1, 6}];
reportRows // Column
```

Expected: no `Abort[]`, no `Message`; a length-6 `Column` of `Association`s, one per state of
`(2,1)`, each showing `M`, `n`, `braket`, `G21 diagonal`, and `normRatio` to 10 digits. Report the
actual `normRatio` values back — this is the exploratory result the plan exists to produce (is
it constant across states, or state-dependent? no assertion either way, just report what comes
out).

- [ ] **Step 3: Save the notebook**

Use `mcp__wolfbook__wolfbook_fileOps` to save `XXX_CG_L2_V2.wb`. This is the final step of the
plan.
