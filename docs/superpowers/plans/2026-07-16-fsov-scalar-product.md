# FSoV Scalar Product Verification Implementation Plan

> **For agentic workers:** none of the standard superpowers execution skills apply here — this
> plan is executed entirely inside a live Wolfram notebook via the project's **wolfbook-builder**
> subagent (see `CLAUDE.md`), not via pytest/git-commit-per-task. Each task specifies notebook
> cells to insert (via `mcp__wolfbook__wolfbook_insertCells`), run (via
> `mcp__wolfbook__wolfbook_runCell`), and verify (via `mcp__wolfbook__wolfbook_getCellOutput` /
> `wolfbook_getKernelState`), against an explicit expected result. Steps use checkbox (`- [ ]`)
> syntax for tracking. There is no separate test file and no per-task git commit.

**Goal:** Port `Experiments/FSoV_L2.wb`'s FSoV measure/bracket/determinant machinery into a new
Section C of `Experiments/XXX_CG_L2_V2.wb`, then use it to numerically verify that
`det[λ1,λ2][Q_A(u)Q_B(u)]` acts as an SoV scalar product on the notebook's existing `Q1` Baxter
Q-functions: vanishing for distinct states, nonzero on the diagonal, across multiple weight-$M$
sectors including cross-$M$ pairs.

**Architecture:** Three notebook-cell tasks appended after Section B (cell index 33, the `Q1`
definition) of `Experiments/XXX_CG_L2_V2.wb`: (1) port + port-fidelity check of `mu`/`Nn`/`nn`/
`br`/`det`, (2) build the full Gram matrix over one representation's states, (3) verify
diagonal/off-diagonal structure and save.

**Tech Stack:** Wolfram Language, via wolfbook MCP tools operating on the live `.wb` notebook and
kernel.

## Global Constraints

- Design doc: `docs/superpowers/specs/2026-07-16-fsov-scalar-product-design.md` — this plan
  implements it in full; do not deviate from its scope (no predicted-norm comparison, no
  multi-representation sweep).
- Test case is fixed: $(\lambda_1,\lambda_2)=(2,1)$, giving 6 states across $M=0,1,2,3$ with
  `dPred[M,2,1] = 1,2,2,1`.
- Numeric parameters ($\theta_1=1/3,\theta_2=1/7,h=1,\kappa_1=N[\exp(i\zeta(3)),30],\kappa_2=
  1/\kappa_1$) are already set in `XXX_CG_L2_V2.wb` cell 0 — do not redefine them; the ported
  `mu`/`Nn`/`nn`/`br`/`det` cells close over these same globals directly.
- Hard-assert discipline, matching this notebook's existing style (`BaxterMatrix`,
  `TauEigensystem`, `QSolve`): every structural/numeric check that could silently hide a bug
  aborts (`Abort[]`) with a `Message`, not a soft warning. If a check fails when actually run,
  **stop and investigate — do not weaken the tolerance or delete the check to force a pass.**
- Comment convention: `(* ... *)` cells explain *why*, not what.
- `wolfbook-builder` stays out of `paper/*.tex` and git commits/pushes (per `CLAUDE.md`) — this
  plan produces no git commits; saving the notebook (`wolfbook_fileOps`) is the only persistence
  step, done once at the end of Task 3.
- Memoized routines this plan adds (`GramMatrix`) should be re-run (along with anything
  downstream) if any upstream definition is edited later — same caching caveat as the rest of the
  notebook.

---

### Task 1: Port FSoV measures/brackets/determinant + port-fidelity check

**Files:**
- Modify: `Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb` (append new cells after cell index
  33)
- Reference (read-only, source of the ported code): `Paul/Mathematica/XXX/Experiments/FSoV_L2.wb`
  cells 0–4

**Interfaces:**
- Consumes: global `θ1,θ2,h` (from `XXX_CG_L2_V2.wb` cell 0).
- Produces: `mu[λ1,λ2][i][u]`, `Nn[λ1,λ2][i]`, `nn[λ1,λ2][i]`, `br[λ1,λ2][i][f]`,
  `det[λ1,λ2][f]` — all consumed by Task 2.

- [ ] **Step 1: Insert the Section C header comment cell**

Insert as a new cell after cell 33:

```wolfram
(*Section C — FSoV scalar product / orthogonality check. mu/Nn/nn/br/det ported verbatim from
  Experiments/FSoV_L2.wb (that file's own cells 0-4): mu is the FSoV measure per site (supported
  on that site's θ-tower), br is the bracket pairing (residue sum of f·mu over a site's poles),
  det is the 2x2 Wronskian-type determinant built from br[f] and br[u f] at both sites. For
  f=Q_A(u)Q_B(u), det[λ1,λ2][f] is the SoV scalar product <Q_A,Q_B>: it must vanish for A≠B and
  be nonzero for A=B. See docs/superpowers/specs/2026-07-16-fsov-scalar-product-design.md.*)
```

Run the cell. Expected: comment-only cell, no output, no error.

- [ ] **Step 2: Insert and run the `mu` cell**

```wolfram
(*FSoV measures*)
ClearAll[mu];
mu[λ1_,λ2_][1][u_]:=Nn[λ1,λ2][1]/(Product[(u-θ1-h k),{k,0,λ1}]Product[(u-θ2-h k),{k,0,λ2}])(1-Exp[2 π I/h (u-θ2)]);
mu[λ1_,λ2_][2][u_]:=Nn[λ1,λ2][2]/(Product[(u-θ1-h k),{k,0,λ1}]Product[(u-θ2-h k),{k,0,λ2}])(1-Exp[2 π I/h (u-θ1)]);
```

Expected: no output (both are `:=` definitions), no error.

- [ ] **Step 3: Insert and run the `Nn`/`nn` cell**

```wolfram
(*Normalization factors in measure*)
Nn[λ1_,λ2_][1]:=1/(2 π I)/(1-Exp[2 π I/h(θ1-θ2)]) nn[λ1,λ2][1];
Nn[λ1_,λ2_][2]:=1/(2 π I)/(1-Exp[2 π I/h(θ2-θ1)]) nn[λ1,λ2][2];

nn[λ1_,λ2_][1]:= (-1)^(λ1+λ2-1)h^(λ1)Gamma[λ1+1] Product[θ2+k h-θ1,{k,0,λ2}]//Simplify;
nn[λ1_,λ2_][2]:= (-1)^(λ1+λ2-1)h^(λ2)Gamma[λ2+1] Product[θ1+k h-θ2,{k,0,λ1}]//Simplify;
```

Expected: no output, no error.

- [ ] **Step 4: Insert and run the `br` cell**

```wolfram
(*brackets of function f, given by contour integral around all poles in measure*)
ClearAll[br];
br[λ1_,λ2_][1][f_]:=2 π I Plus@@(Residue[f mu[λ1,λ2][1][u],{u,#}]&/@Table[θ1+h k,{k,0,λ1}])//Expand//Simplify;
br[λ1_,λ2_][2][f_]:=2 π I Plus@@(Residue[f mu[λ1,λ2][2][u],{u,#}]&/@Table[θ2+h k,{k,0,λ2}])//Expand//Simplify;
```

Expected: no output, no error.

- [ ] **Step 5: Insert and run the `det` cell**

```wolfram
ClearAll[det];
det[λ1_,λ2_][f_]:=1/(f/.u:>θ1)/(f/.u:>θ2)/(θ2-θ1)({{br[λ1,λ2][1][f],br[λ1,λ2][1][u f]},{br[λ1,λ2][2][f],br[λ1,λ2][2][u f]}}//Det//Expand//Simplify)
```

Expected: no output, no error.

- [ ] **Step 6: Insert and run the port-fidelity check cell**

This transcribes `FSoV_L2.wb`'s own cached test result (its cell 4, `Out[215]`, produced when
`θ1,θ2,h` were still free symbols there) as a literal closed-form expression, and checks that
running the *same* test (`det[1,1][f[u] z^(u/h)]`) in this notebook — where `θ1,θ2,h` are
already bound to `1/3,1/7,1` — agrees with that formula after substituting those same numeric
values in. This catches any transcription slip in Steps 2–5 without needing to hand-simplify
anything:

```wolfram
(*Port-fidelity check: fsovL2CachedFormula is FSoV_L2.wb's own cached Out[215] for
  det[1,1][f[u] z^(u/h)], transcribed verbatim with θ1,θ2,h,z,f still symbolic (as they were in
  that standalone notebook). Here θ1,θ2,h are already bound (cell 0), so comparing the freshly
  computed det[...] against this formula with the SAME substitution applied is a direct
  transcription check, independent of re-deriving the algebra.*)
ClearAll[fsovL2CachedFormula];
fsovL2CachedFormula[θ1_,θ2_,h_,z_,fn_] := -(z^((θ1+θ2)/h-θ1/h-θ2/h) (z fn[h+θ1] (fn[θ2] (h-θ1+θ2)+z (θ1-θ2) fn[h+θ2])+fn[θ1] (fn[θ2] (θ1-θ2)-z fn[h+θ2] (h+θ1-θ2))))/(fn[θ1] fn[θ2] (θ2-θ1));
portCheckDiff = Simplify[det[1,1][f[u] z^(u/h)] - fsovL2CachedFormula[θ1,θ2,h,z,f]]
```

Expected: `portCheckDiff` evaluates to `0`. If it is anything else, **stop** — one of Steps 2–5
has a transcription error relative to `FSoV_L2.wb`; diff the ported cell text against the
original file before continuing.

- [ ] **Step 7: Save the notebook**

Use `mcp__wolfbook__wolfbook_fileOps` to save `XXX_CG_L2_V2.wb`.

---

### Task 2: Build the Gram matrix over one representation's states

**Files:**
- Modify: `Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb` (append new cells after Task 1's
  last cell)

**Interfaces:**
- Consumes: `TauEigensystem[λ1,λ2]` (Section B, cell 30 — memoized `Association` keyed by
  `{M,n}`, each value has `"J"`/`"Psi"`), `Q1[λ1,λ2][M,n][u]` (Section B, cell 33), `det`
  (Task 1).
- Produces: `statesOf[λ1,λ2]` (sorted list of `{M,n}` keys), `GramMatrix[λ1,λ2]` (function
  returning a numeric matrix), `states21`, `G21` — all consumed by Task 3.

- [ ] **Step 1: Insert and run the `statesOf`/`GramMatrix` cell**

```wolfram
(*statesOf: the sorted {M,n} key list for a representation, read off TauEigensystem's memoized
  cache (built once in Section B) rather than re-deriving dPred here.
  GramMatrix: the full Gram matrix G_{(M,n),(M',n')} = det[λ1,λ2][Q1(M,n)(u) Q1(M',n')(u)] over
  ALL states of the representation at once — this covers same-M degeneracy AND cross-M pairs in
  a single object, per the design.*)
ClearAll[statesOf, GramMatrix];
statesOf[λ1_,λ2_] := SortBy[Keys[TauEigensystem[λ1,λ2]], {First,Last}];
GramMatrix[λ1_,λ2_] := Module[{s = statesOf[λ1,λ2]},
  Table[
    det[λ1,λ2][Q1[λ1,λ2][s1[[1]],s1[[2]]][u] Q1[λ1,λ2][s2[[1]],s2[[2]]][u]],
    {s1, s}, {s2, s}]
];
```

Expected: no output, no error.

- [ ] **Step 2: Insert and run the states smoke test**

```wolfram
states21 = statesOf[2, 1]
```

Expected: `{{0,0},{1,0},{1,1},{2,0},{2,1},{3,0}}` — the 6 states of `(λ1,λ2)=(2,1)` (`dPred`
$=1,2,2,1$ for $M=0,1,2,3$), sorted by `M` then `n`. If the count or ordering differs, **stop** —
`TauEigensystem[2,1]` or `dPred` has changed since Section B was built; do not proceed until this
matches.

- [ ] **Step 3: Insert and run the Gram matrix computation**

```wolfram
G21 = GramMatrix[2, 1];
Dimensions[G21]
```

Expected: `{6, 6}`. This may take a little while (36 `det` evaluations, each a `Residue`-based
symbolic-then-numeric computation) — that is expected, not a hang.

- [ ] **Step 4: Save the notebook**

Use `mcp__wolfbook__wolfbook_fileOps` to save `XXX_CG_L2_V2.wb`.

---

### Task 3: Verify orthogonality/non-degeneracy and report precision

**Files:**
- Modify: `Paul/Mathematica/XXX/Experiments/XXX_CG_L2_V2.wb` (append new cell after Task 2's
  last cell)

**Interfaces:**
- Consumes: `G21` (Task 2), `states21` (Task 2).
- Produces: printed verification report (terminal deliverable of this plan — nothing further
  consumes this programmatically).

- [ ] **Step 1: Insert and run the verification cell**

```wolfram
(*Verification: G21 should be numerically diagonal. maxOffDiag is the largest |off-diagonal|
  entry over all 30 off-diagonal pairs (same-M AND cross-M, since states21 spans all of M=0..3);
  minDiagAbs is the smallest |diagonal| entry, confirming non-degeneracy (a genuine scalar
  product, not one that degenerates to 0 on some state). Threshold 10^-15 matches this notebook's
  existing hard-assert convention (e.g. Section B's TauEigensystem residual check, cell 30) — well
  above the ~10^-20..10^-27 noise floor expected at this notebook's 30-digit working precision,
  and far below any genuine O(1)-scale diagonal entry.*)
FSoVGram::notorthogonal = "Off-diagonal Gram entries do not vanish for (λ1,λ2)=(2,1): max|off-diag|=`1`, min|diag|=`2`.";
FSoVGram::degenerate = "A diagonal Gram entry is suspiciously small for (λ1,λ2)=(2,1) (min|diag|=`1`) — scalar product may be degenerate.";
offDiagVals = Flatten[Table[If[i != j, G21[[i, j]], Nothing], {i, 1, 6}, {j, 1, 6}]];
diagVals = Table[G21[[i, i]], {i, 1, 6}];
maxOffDiag = Max[Abs[N[offDiagVals]]];
minDiagAbs = Min[Abs[N[diagVals]]];
If[maxOffDiag > 10^-15,
  Message[FSoVGram::notorthogonal, maxOffDiag, minDiagAbs];
  Abort[]];
If[minDiagAbs < 10^-15,
  Message[FSoVGram::degenerate, minDiagAbs];
  Abort[]];
<|
  "states" -> states21,
  "max|off-diagonal|" -> maxOffDiag,
  "digits vanished (-Log10)" -> N[-Log10[maxOffDiag]],
  "min|diagonal|" -> minDiagAbs,
  "diagonal values" -> diagVals
|>
```

Expected: no `Abort[]`, no `Message`; the printed `Association` shows `"max|off-diagonal|"` at
or below `~10^-15` (report the actual value and its `-Log10` digit count — expect it to land
somewhere in the `~20`–`~30` digit range given the notebook's 30-digit working precision, similar
in spirit to the L3 notebook's Baxter-residual check landing at 18–23 tracked digits), and
`"diagonal values"` a length-6 list of numbers clearly not near that floor (genuine nonzero
norms). **If either assertion fires, stop and investigate** — per the Global Constraints, do not
raise the `10^-15` threshold or remove the check to force a pass; a real failure here means
either the FSoV port (Task 1) or the `Q1` Q-functions themselves need to be revisited.

- [ ] **Step 2: Save the notebook**

Use `mcp__wolfbook__wolfbook_fileOps` to save `XXX_CG_L2_V2.wb`. This is the final step of the
plan.
