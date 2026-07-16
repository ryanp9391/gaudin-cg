# Hand-off

Where the last session left off. Read this first when picking the project back up; it's not a
history log (that's `git log` / `logbook.tex`) — only what's still open or in-flight. Overwrite
freely each session rather than appending.

## Current focus

Computing CG coefficients directly in the XXX spin chain (`Paul/Mathematica/XXX/`). The Gaudin
subproject is parked.

## State as of 2026-07-16 (SoV factorization + functional CG overlap — DONE)

Big result this session: the **Clebsch–Gordan overlap correspondence closed exactly**. Design →
plan → wolfbook-builder execution, in `Experiments/XXX_CG_L2_V2.wb`. Spec:
`docs/superpowers/specs/2026-07-16-xxx-sov-factorization-functional-overlap-design.md`; plan:
`docs/superpowers/plans/2026-07-16-xxx-sov-factorization-functional-overlap.md` (both have a
`Result` / `Execution outcome` section recording what actually happened). Everything below is
**committed**.

### Section E — operatorial factorization (cells 53–60, commit `8dbed52`)

Built a genuine standalone **length-1** fused transfer-matrix hierarchy `t1[λ][a,s]`/`qdet1` and SoV
covector `xL1[λ][{s}]` (validated as a left B-eigencovector, `{0}` for λ=1..3). Proved the s2=0
slice of the L=2 SoV basis factorizes:

    x2[λ1,λ2][{s1,0}] = c (xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}])

ratio bit-exactly constant over all 27 (λ1,λ2,s1). **Key finding (corrected a wrong guess of mine):**
the site-2 pseudovacuum is NOT passed through as the identity — it is a left eigencovector of the
site-2 factor of `t2[1,s1][θ1]` with eigenvalue
`μSite2[λ2,s1] = Product[(θ1-θ2+(k-1)h-λ2 h),{k,1,s1}]`, so
`c = Product[(θ1-θ2+(k-1)h-λ2 h)/(θ1-θ2+(k-1)h-λ1 h),{k,1,s1}]` (=1 at λ1=λ2). All hard-asserted.

### Section F — functional CG overlap (cells 61–75)

The deliverable. **The 1×1 FSoV determinant of the twisted Q-function equals the operatorial overlap
of a site-1 unipotent rotation with the two-site eigenvector:**

    det1[λ1][ k^(u/h) Q1[λ1,λ2][M,n] ] = (⟨x0[λ1]| exp[-k E21] ⊗ ⟨x0[λ2]|) . Psi[λ1,λ2][M,n]

where `det1[λ1][f] := br1[λ1][f]/(f/.u->θ1)` (Section C's `det` restricted to one site). Dictionary:
**generator E21, φ(k) = -k** (α=-1, universal). Residual is the literal symbolic `0` across the rep
sweep {(1,1),(2,1),(1,2),(2,2),(3,1),(1,3)} — 41 states.

**The crux / the thing that nearly went wrong:** the first Phase-2 pass used the bare bracket
`br1[k^(u/h) Q1]` (no determinant normalization) and did NOT close (k-dependent per-state scalar).
The fix was the `1/f(θ1) = 1/(k^(θ1/h) Q1(θ1))` normalization built into the FSoV determinant — it
cancels the `k^(θ1/h)` prefactor and the state-dependent `Q1(θ1)`, leaving `det1[k^(u/h) Q1]` a clean
degree-λ1 polynomial in k starting at 1. Caught by the user asking "was it normalized by
k^θ1/h Q1(θ1)?". The plan's "the 1×1 determinant collapses to a single bracket" wording was the bug;
now corrected in both spec and plan.

**Gotcha found this session:** `MatrixPower[nilpotentMatrix, 0]` fires `MatrixPower::sing` and stays
unevaluated even for a concrete integer power 0. Build `exp[φ E]` for nilpotent E as
`IdentityMatrix[dim] + Sum[φ^m/m! MatrixPower[E,m],{m,1,λ}]` (split off the identity), not `MatrixExp`
and not a Sum starting at m=0.

## Open threads / next steps

- **Nothing queued.** Natural next steps if continuing: (a) write this up in `paper/xxxCG.tex`
  (still single-site only — no L=2 factorization, scalar product, or the CG-overlap dictionary
  written up yet); (b) extract the actual CG coefficients from the overlap generating function (the
  `exp[-k E21]` rotation generates the CG combinations — expanding in k gives per-state coefficients);
  (c) push the factorization to the s2≠0 case / general L; (d) run `logbook-agent` to log today.
- **Docs to promote:** `XXX_CG_L2_V2.wb` is in `Experiments/`; nothing new promoted to `Clean/` this
  session. If Section E/F get promoted, `wolfbook-docs-agent` can regenerate the `.md`/`.tex` refs.
- Older restructuring/Clean-notebook state (L3 build, docs agent, `\mathbbm{1}` LaTeX gotcha, shared
  wolfbook kernel gotcha) is unchanged from the 2026-07-15 session — see `git log` and CLAUDE.md.

## Stable structure

See `CLAUDE.md` for the full layout (XXX vs Gaudin split, Clean/Experiments convention, subagents).
Key for this thread: `XXX_CG_L2_V2.wb` now has Sections A (SoV construction) · B (Baxter/TQ, Psi/Q1) ·
C (2-site FSoV det) · D (left eigenvectors) · E (L=1 factorization) · F (functional CG overlap).
