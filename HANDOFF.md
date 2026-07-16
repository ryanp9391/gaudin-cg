# Hand-off

Where the last session left off. Read this first when picking the project back up; it's not a
history log (that's `git log` / `logbook.tex`) — only what's still open or in-flight. Overwrite
freely each session rather than appending.

## Current focus

Computing CG coefficients directly in the XXX spin chain (`Paul/Mathematica/XXX/`). The Gaudin
subproject is parked.

## State as of 2026-07-16 (SoV factorization + functional CG overlap — DONE, promoted to Clean)

Big result this session: the **Clebsch–Gordan overlap correspondence closed exactly**, then
promoted from `Experiments/` into `Clean/`. Everything below is **committed**. Design/plan docs
(with `Result`/`Execution outcome` sections recording what actually happened):
`docs/superpowers/specs/2026-07-16-xxx-sov-factorization-functional-overlap-design.md`,
`docs/superpowers/plans/2026-07-16-xxx-sov-factorization-functional-overlap.md`.

### The result

**The 1×1 FSoV determinant of the twisted Q-function equals the operatorial overlap of a site-1
unipotent rotation with the two-site eigenvector:**

    det1[λ1][ k^(u/h) Q1[λ1,λ2][M,n] ] = (⟨x0[λ1]| exp[-k E21] ⊗ ⟨x0[λ2]|) . Psi[λ1,λ2][M,n]

where `det1[λ1][f] := br1[λ1][f]/(f/.u->θ1)` (the 1-site FSoV determinant — Section C's `det`
restricted to one site). Dictionary: **generator E21, φ(k) = -k** (α=-1, universal). Residual is
the literal symbolic `0` across the rep sweep {(1,1),(2,1),(1,2),(2,2),(3,1),(1,3)} — 41 states.

This rests on a prerequisite factorization: the s2=0 slice of the L=2 SoV basis,
`x2[λ1,λ2][{s1,0}] = c (xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}])`, where `xL1` is a genuine standalone
length-1 SoV covector (new `t1`/`qdet1` transfer-matrix hierarchy) and the site-2 pseudovacuum
passes through `t2[1,s1][θ1]` with eigenvalue `μSite2[λ2,s1] = Product[(θ1-θ2+(k-1)h-λ2 h),
{k,1,s1}]` — **not** the identity, a wrong guess corrected mid-session — giving
`c = Product[(θ1-θ2+(k-1)h-λ2 h)/(θ1-θ2+(k-1)h-λ1 h),{k,1,s1}]`.

**The crux / what nearly went wrong:** the first functional-discovery pass used the bare bracket
`br1[k^(u/h) Q1]` (no determinant normalization) and did NOT close (k-dependent per-state scalar).
The fix was the `1/f(θ1) = 1/(k^(θ1/h) Q1(θ1))` normalization built into the FSoV determinant — it
cancels the `k^(θ1/h)` prefactor and the state-dependent `Q1(θ1)`, leaving `det1[k^(u/h) Q1]` a clean
degree-λ1 polynomial in k starting at 1. Caught by the user asking "was it normalized by
k^θ1/h Q1(θ1)?". If a functional-SoV object won't universalize, check for a dropped determinant
node-normalization first (see memory `reference_fsov_determinant_normalization`).

**Gotcha found this session:** `MatrixPower[nilpotentMatrix, 0]` fires `MatrixPower::sing` and stays
unevaluated even for a concrete integer power 0. Build `exp[φ E]` for nilpotent E as
`IdentityMatrix[dim] + Sum[φ^m/m! MatrixPower[E,m],{m,1,λ}]` (split off the identity), not `MatrixExp`
and not a Sum starting at m=0.

### Commits (chronological)

`737f6a2` spec → `1bd9bcd` plan → `8dbed52` Section E in `Experiments/XXX_CG_L2_V2.wb` → `9db89eb`
Section F in `Experiments/` + corrected docs → `dcf5c1e` **promoted Sections C,D,E,F into
`Clean/XXX_CG_L2_Clean.wb`** (now a contiguous 71-cell A→F notebook, fresh-kernel top-to-bottom run
verified, every hard-assert exact) → `bfc0ffb` refreshed `Clean/XXX_CG_L2_Clean.md`/`.tex` reference
docs to cover C–F (`.tex` compiles clean with plain `pdflatex`, `\mathbbm{1}` convention).

Section F was promoted in **final form only**: the corrected `det1`/`FLdet` machinery and the closed
dictionary. The exploratory bare-bracket `FL`, the pre-`det1` per-state-scalar cells, and the
`obstructionReport` were dropped from Clean (they remain in `Experiments/XXX_CG_L2_V2.wb` as the
derivation record, same convention as other Clean promotions).

## Open threads / next steps

- **Nothing queued.** Natural next steps if continuing:
  - Write this up in `paper/xxxCG.tex` (still single-site only — no L=2 factorization, scalar
    product, or the CG-overlap dictionary written up yet).
  - Extract the actual CG coefficients from the overlap generating function (the `exp[-k E21]`
    rotation generates the CG combinations — expanding in k gives per-state coefficients).
  - Push the factorization to the s2≠0 case / general L.
  - Run `logbook-agent` to log today.
- `Experiments/XXX_CG_L2_V2.wb` still holds the full discovery narrative (bare-bracket attempt,
  obstruction report, per-state-scalar exploration) — useful if revisiting why the determinant
  normalization matters.
- Older restructuring/Clean-notebook state (L3 build, docs agent, shared wolfbook kernel gotcha) is
  unchanged from the 2026-07-15 session — see `git log` and `CLAUDE.md`.

## Stable structure

See `CLAUDE.md` for the full layout (XXX vs Gaudin split, Clean/Experiments convention, subagents).
`Clean/XXX_CG_L2_Clean.wb` is now the authoritative A–F notebook: A (SoV construction) · B
(Baxter/TQ, Psi/Q1) · C (2-site FSoV scalar product/orthogonality) · D (left eigenvectors,
`<Psi|Psi>=det[Q,Q]`) · E (L=1 factorization + site-2 eigenvalue) · F (functional CG-overlap
dictionary). Companion `.md`/`.tex` reference docs are up to date with it.
