# Hand-off

Where the last session left off. Read this first when picking the project back up; it's not a
history log (that's `git log` / `logbook.tex`) — only what's still open or in-flight. Overwrite
freely each session rather than appending.

## Current focus

Computing CG coefficients directly in the XXX spin chain (`Paul/Mathematica/XXX/`). The Gaudin
subproject is parked.

## State as of 2026-07-19 (gl(4) L=2 untwisted Casimir eigenvalues — DONE, exploratory findings mixed, in Experiments)

Follow-on to the twisted build below: new section appended to `Baxter_L2_XXX_SU4.wb` (cells 35–47),
`z[k]→1` (untwisted/physical point), for the two-row-rectangular rep family `λ1=[S,S,0,0]`,
`λ2=[P,P,0,0]`. Goal: closed-form eigenvalues `Λ_a(u)` of the untwisted transfer matrices via su(4)
Casimir invariants of the irreps `μ` in `λ1⊗λ2=⊕N_μμ`. Design/plan (with gotchas/results recorded
inline as found): `docs/superpowers/specs/2026-07-19-su4-l2-untwisted-casimir-eigenvalues-design.md`,
`docs/superpowers/plans/2026-07-19-su4-l2-untwisted-casimir-eigenvalues.md`. All 5 tasks executed
inline (subagents again lacked wolfbook tools this session).

**Exact-arithmetic simplification:** at `z=1` (vs. the twisted case's irrational `Zeta`-based twists)
every quantity is exact rational — no precision-cascade concerns at all in this section.

**Verified result:** on the "top" state `μ=λ1+λ2` (always multiplicity 1, by general representation
theory — confirmed empirically too: this rep family is **fully multiplicity-free** across the whole
`(S,P)∈{(1,1),(2,1),(1,2),(2,2)}` sweep, up to `dim=400`), `Λ_1(u) = ν[1](u)+ν[2](u)+ν[3](u)+ν[4](u)`
exactly (`worst~0×10⁻²⁵`) — no Casimir needed, follows directly from the pseudovacuum `ν` formula
already built for the twisted case.

**Falsified conjecture (honest negative result, not a bug):** "each `u^j` coefficient of `Λ_1(u)` is
an affine function of `C2(μ)` (quadratic Casimir) alone" — tested two independent ways (2-anchor
extrapolation and least-squares over all 167 multiplicity-free data points across 3 rep pairs) —
**fails** for `j=0` (constant term, deviation ~3–4.7) and `j=1` (subleading term, deviation ~1.25–2.0).
Only `j=2` (leading term) is trivially constant (`=χ1=4`, already known, not a real test). **Stopped
here per convention rather than searching further candidates unprompted** — next step (quadratic-in-C2,
a second invariant like `H2` itself, or a cubic/quartic Casimir) needs a decision before continuing.

**Three real bugs found and fixed while executing** (all documented inline in the plan, useful
precedent for future notebook work): `ClearAll[f]` wipes a function's rules entirely, not just its
memoized cache (had to re-insert the original defining code, not just clear+recall);
`assoc[Key[list]]` single-bracket indexing silently returns `KeyAbsent` even for a genuine key here
(use `assoc[list]`, no wrapper); computing `"field"->expr&/@list` *inline inside* an association
literal silently binds to a scalar instead of a vector (compute to a local variable first, as
`TauEigensystem4` already did). Also: `Table[...,{S=SP[[1]],P=SP[[2]]}]` is invalid iterator syntax
(silently no-ops instead of erroring loudly — `S,P` need to be `Module` locals); `Flatten[list,2]`
over-flattens when list rows have internal list structure (use depth matched to actual nesting).

**Next natural steps:** decide whether to keep searching for a closed form on `j=0,1` (candidates:
quadratic-in-`C2`, a second commuting invariant, cubic/quartic Casimirs `C3op`/`C4op` — not yet
built) or treat the top-state result + falsification as the reportable outcome for this rep family.

## State as of 2026-07-19 (gl(4) L=2 twisted Baxter solver: `Baxter_L2_XXX_SU4.wb` — DONE, in Experiments)

Ported the two-site gl(3) XXX pipeline (`su3_V2_general.wb`) to **gl(4)**: two independent generic
Young-diagram reps `λ={λ1,λ2,λ3,λ4}` per site, four **distinct** diagonal-twist eigenvalues
`z[k]=N[Exp[I Zeta[2k+1]],70]` (k=1..4), 4th-order Baxter/TQ equation, four Q-functions, built fresh
in `Paul/Mathematica/XXX/Experiments/Baxter_L2_XXX_SU4.wb` (was empty). Design/plan (with deviation
notes recording what actually happened): `docs/superpowers/specs/2026-07-19-su4-l2-twisted-baxter-design.md`,
`docs/superpowers/plans/2026-07-19-su4-l2-twisted-baxter.md`. All 10 tasks DONE, executed inline
(not via subagent — the wolfbook MCP tools weren't wired into subagent tool sets this session, see
plan/session notes), controller commits.

**Companion twist dropped entirely** (user direction, confirmed before Task 2): no `G`/`TG`/`tG`/
`qminG` — only the diagonal-twist branch (`Tg`/`tg`/`qming`) is built, matching what the gl(3) source
pipeline actually used downstream anyway. `χ1..χ4` (elementary symmetric functions of `z[1..4]`)
still built directly, needed by `qdetT` and the leading-coeff check.

**Genuine gl(4) math (everything else is pure index/term extension from gl(3)):**
- `ν[λ1,λ2][i][u]=(u−θ1−hλ1[[i]])(u−θ2−hλ2[[i]])` for `i=1..4`; `qdetT=χ4·ν[1][u−3h]ν[2][u−2h]ν[3][u−h]ν[4][u]`.
- 4th-order `BaxOp`: `f[u+h]−A1 f[u]+A2 f[u−h]−A3 f[u−2h]+A4 f[u−3h]`.
- **Degree formula carries over unchanged**: `M_k=Λ_1−n_k` (`Λ=λ1+λ2`, `Λ_1=(λ1+λ2)[[1]]`) — verified
  with **zero deviation** against the derived asymptotic `Mdeg` master formula over the full 152-state
  sweep; no candidate-fit re-derivation was needed (the plan's fallback path was never triggered).
- `TauEigensystem4` simultaneously diagonalizes `τ1,τ2,τ3` (Cartan `H1,H2,H3`); `QSolve4` solves for
  four monic `q_k`, `Qfun`=`z[k]^(u/h) q_k(u)`.

**Full-sweep verification (acceptance):** 5 rep pairs, **152 states total**
(`{1,0,0,0}⊗{1,0,0,0}`=16, `{1,0,0,0}⊗{1,1,1,0}`=16, `{1,1,0,0}⊗{1,0,0,0}`=24,
`{1,1,0,0}⊗{1,1,0,0}`=36, `{2,1,1,0}⊗{1,0,0,0}`=60 adjoint⊗fund), every `Q_k` solves the 4th-order TQ
equation and the 5×4 Casoratian reconstruction of `A1..A4` matches `τ1,τ2,τ3,qdetT` —
**worstTQ = worstCasoratian = 0 exactly**. Verified by a fresh-kernel top-to-bottom run (restart → run
cells 1→34 in order); every hard-assert passed in-kernel with no `Abort`.

**Precision bug found and fixed (the real story of this session):** the plan specified `z[k]` at
30-digit precision (matching gl(3)); the initial full sweep at that precision left the d=60
adjoint⊗fund rep's Casoratian residual at `~4.5×10⁻⁸` — 45x over the `10⁻⁹` bar — even though the TQ
residual was *exactly* 0 everywhere (proof the formulas/implementation were correct; this was a pure
numerics issue). Bumping `z[k]` to 50 then 70 digits barely moved the residual (`~1.6×10⁻⁹`, then
bit-identical again) — the real giveaway that something else was capping precision independent of
`z`. Root cause: `TauEigensystem4`'s `rSets` (generic-combination coefficients for splitting the
degenerate τ-eigenspace) and its `Kmix` eigendecomposition were **hardcoded to 30-digit precision** —
a leftover copied verbatim from the gl(3) source's `TauEigensystem3`, never tied to the `z[k]`
precision bump. This silently capped the entire downstream chain (τ-coefficients, `Qfun`, Casoratian)
at ~19-25 accurate digits no matter how much precision `z` carried. Fix: match both to `z`'s 70-digit
precision — eigen-residual improved from a capped `8.6×10⁻²⁶` to `3.46×10⁻⁶⁴`, and the full sweep then
passed with **exact** zero residuals. Lesson for future gl(N) ports: hardcoded precision literals
copied from a smaller-N source notebook are a silent trap — grep for bare `N[..., <number>]` whenever
raising a global precision parameter.

**Next natural steps:** promote to `Clean/` once satisfied; the deliberate follow-on is the `z→1`
untwisted limit (su(4) irrep multiplet degeneracy, feeds the CG-overlap program) — explicitly out of
scope for this notebook per the design spec.

## State as of 2026-07-17 (gl(3) general reps per site: `su3_V2_general.wb` — DONE, in Experiments)

Generalized the two-site gl(3) XXX notebook from symmetric reps `{S,0,0}` per site to **arbitrary
gl(3) Young-diagram reps `λ={λ1,λ2,λ3}`** (`λ1≥λ2≥λ3≥0`, integers), one independent diagram per
site, in `Paul/Mathematica/XXX/Experiments/su3_V2_general.wb` (a copy of `su3_V1.wb` generalized in
place). Diagonal-twist frame, `(H1,H2,n)` labelling, three Baxter Q-functions, and the hard-assert
zero-residual convention all preserved. Design/plan (with Execution-outcome section):
`docs/superpowers/specs/2026-07-17-su3-general-reps-design.md`,
`docs/superpowers/plans/2026-07-17-su3-general-reps.md`. Tasks 1–6 all DONE; controller commits.

**Interface:** every function is re-signatured `S1,S2 → λ1,λ2` and takes full 3-vectors per site
(no scalar-`S` wrappers). Numeric params unchanged (`θ1=1/3, θ2=1/7, h=1`, high-precision `z[k]`).

**The three genuine math edits (everything else is pure relabeling):**
- **General Yangian weight** `ν[λ1,λ2][i][u] = (u−θ1−h·λ1[[i]])(u−θ2−h·λ2[[i]])` for `i=1,2,3`,
  from the Lax action `L[λ][i,i][u]|hw⟩=(u−h λ_i)|hw⟩` on the two-site HWS. Verified via the HWS
  eigenvalue check (T[i,i]|hw⟩=ν[i]|hw⟩) and the `qdetT == tg[·][3,1]` shift-ordering guard on
  non-symmetric reps. `qdetT[λ1,λ2][u]=χ3 ν[·][1][u−2h] ν[·][2][u−h] ν[·][3][u]` (V1 shift order
  survives).
- **Honest tensor-product `wtMult`:** `siteWeightMult[λ]=Counts` of the gl(3) weight of every GT
  pattern of `λ`; `wtMult[λ1,λ2][H1,H2]` recovers `(n1,n2,n3)` from `(H1,H2)` with `N=Total[λ1]+
  Total[λ2]` and convolves the two single-site tallies — correctly yields interior multiplicity `>1`.
  Sum rule `Σ wtMult == dim[λ1] dim[λ2]` hard-asserted across the sweep; guard on `{2,1,0}⊗{2,1,0}`
  (d=64) passes.
- **Degree closed form (Task 5 headline):** `M_k = Λ_1 − n_k`, with `Λ = λ1+λ2` the pseudovacuum
  highest weight and `Λ_1 = (λ1+λ2)[[1]]`. Sum rule `Σ_k M_k = 3Λ_1 − N`. Reduces to V1's `M_k = N−n_k`
  for symmetric reps (there `Λ_1 = N`). Extracted by fitting candidates `{N−n, Λ_1−n, Λ_k−n}` against
  the rep-agnostic `Mdeg` over the full sweep — only `Λ_1−n` gives zero deviation on the non-symmetric
  reps. Hard-asserted (`Mpred==Mdeg`, worst=0).

**Full-sweep verification (acceptance, Task 6):** six rep pairs, **102 states total**, every `Q_k`
solves the 3rd-order TQ equation and the 4×4 Casoratian reconstruction of `A1,A2,A3` matches
`τ1,τ2,qdetT` — **worstTQ = worstCasoratian = 0** (machine-zero, both `< 10^-9`). Per-rep counts:

| `λ1` | `λ2` | nstates |
|---|---|---|
| `{1,0,0}` | `{1,0,0}` | 9 |
| `{2,0,0}` | `{1,1,0}` | 18 |
| `{1,1,0}` | `{1,1,0}` | 9 |
| `{2,1,0}` | `{1,0,0}` | 24 |
| `{2,1,0}` | `{1,1,0}` | 24 |
| `{2,2,0}` | `{1,0,0}` | 18 |

Verified by a fresh-kernel top-to-bottom run (restart → run cells 1→40 in order); every hard-assert
check passed in-kernel with no `Abort`.

**wolfbook output-persistence glitch — re-run, do not trust blank boxes.** Several check cells created
via `insertCells` / edited during Tasks 3–5 do NOT refresh their on-disk `Out[]` display even after a
correct fresh-kernel re-run (confirmed: re-running an affected cell individually still shows blank).
These cells computed the correct result in-kernel (verified by direct `evaluateExpression`), but their
saved display is empty or stale. **A future reader must re-run the notebook rather than treat these as
failures.** Cells still showing empty/stale output on disk after the full Task-6 run (1-based cell
numbers; correct in-kernel values in parentheses):
- Cell 14 — `qdetT == tg[·][3,1]` guard: **empty** (in-kernel `{0}`).
- Cell 23 — `wtMult` sum-rule: **empty** (in-kernel all reps OK; the informational
  `adjoint⊗fund @(0,0)` spot value is `0` because for N=4 the (0,0) weight is non-integer — not a
  failure).
- Cell 24 — `{2,1,0}⊗{2,1,0}` d=64 guard: **empty** (in-kernel `tot=64=dim`).
- Cell 27 — eigensystem residual spot check on `{2,1,0}⊗{1,0,0}`: **empty** (in-kernel worst≈1.5×10⁻²⁷).
- Cell 31 — `degData` tabulation: **stale** (shows an old `{H1H2n→…,sumM→…}` format; source now emits
  raw `{λ1,λ2,H1,H2,n,{M}}` rows).
- Cell 32 — degree candidate-fit `results` Association: **stale** (shows an old matrix; source now
  emits `{"N-n"→…,"Lam1-n"→0,"Lamk-n"→…}`).
- Cell 40 — full-sweep report: **stale** (shows the old symmetric sweep, 141 states; the current
  λ-vector run gives the 102-state table above with worstTQ=worstCasoratian=0).

Check cells whose on-disk `Out[]` IS correct: cell 4 (`C1`→{0}), cell 10 (HWS→{0}), cell 21
(commuting-family, maxCommutator≈10⁻²⁶), cell 35 (degree assert worst→0; only its result-label string
still reads "N-n_k", cosmetic), cell 39 (single-state TQ perQresidual→{0,0,0}).

**Next natural steps:** promote `su3_V2_general.wb` to `Clean/` once satisfied; then the companion-twist
SoV / FSoV CG-overlap program for general gl(3) reps (mirroring the su(2) Sections C–F), and the
`paper/xxxCG.tex` write-up.

## State as of 2026-07-17 (gl(3) Section B: eigensystem + Baxter Q-functions — DONE, in Experiments)

Built the gl(3) analog of the su(2) Clean notebook's Section B, in
`Paul/Mathematica/XXX/Experiments/su3_V1.wb` (on top of the pre-existing Section A there). All
committed to `main`, all hard-asserts zero-residual. Design/plan (with Execution-outcome sections):
`docs/superpowers/specs/2026-07-17-su3-transfer-matrix-eigensystem-baxter-design.md`,
`docs/superpowers/plans/2026-07-17-su3-transfer-matrix-eigensystem-baxter.md`.

**What's there (Section B cells, diagonal-twist frame, numeric `θ1=1/3,θ2=1/7,h=1`,
`z[k]=N[Exp[I Zeta[3/5/7]],30]`):**
- `tgCoeffs` (coeff matrices of `τ1=tg[·][1,1]`, `τ2=tg[·][2,1]`), Cartan `H1op`/`H2op`, commuting-
  family check (`maxCommutator=0`).
- `wtMult[S1,S2][H1,H2]` — weight-multiplicity predictor (sum rule verified).
- `TauEigensystem3[S1,S2]` — simultaneous diagonalization; every state labeled `(H1,H2,n)`
  (`H1=E11−E22, H2=E22−E33`; `n` sorted by the full `τ1/τ2` charge tuple), sector counts match
  `wtMult`. Accessors `Psi`/`tau1`/`tau2`.
- `BaxOp` (TQ residual functional), `alphaCoeffs`, `Mdeg`, `qpoly`, `QSolve3` (three-Q solver),
  `Qfun` accessor `Q_k = z_k^{u/h} q_k(u)`.

**Key result — the degree relation:** `deg q_k = M_k = N − n_k` (`N=S1+S2`, `(n1,n2,n3)` the gl(3)
weight components), equivalently `M1=(2N−2H1−H2)/3, M2=(2N+H1−H2)/3, M3=(2N+H1+2H2)/3`. `M_k` pairs
with `z_k` in-index, depends only on `(H1,H2)`, `Σ M_k = 2N`. Hard-asserted.

**Verification:** full sweep `{(1,1),(2,1),(1,2),(2,2),(3,1),(1,3)}`, 141 states — every `Q_k`
solves the 3rd-order TQ equation (`worstTQ=0`) and the 4×4 Casoratian reconstruction of `A1,A2,A3`
matches `τ1,τ2,qdetT` (`worstCasoratian=0`). Commits `b4b68c7`…`645517d`.

**Next natural steps:** promote `su3_V1.wb` Section B to `Clean/` once satisfied; build the
companion-twist SoV covectors / FSoV scalar product for gl(3) (the CG-overlap program, mirroring
the su(2) Sections C–F); write up in `paper/xxxCG.tex`.

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
`x2[λ1,λ2][{s1,0}] = xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}]` **exactly** (no leftover scalar), where
`xL1` is a genuine standalone length-1 SoV covector (new `t1`/`qdet1` transfer-matrix hierarchy).

**Normalization correction (same day, post-Clean-promotion):** the first pass of this factorization
check found an apparently nontrivial site-2 "passthrough eigenvalue" `μSite2[λ2,s1] =
Product[(θ1-θ2+(k-1)h-λ2 h),{k,1,s1}]`. **This was wrong** — a bug in `x2`'s own normalization, not
real physics, found by the user. `x2`'s denominator used `Qθ[u]=(u-θ1)(u-θ2)` evaluated at a point
shifted by `-λ1h`, which incorrectly puts that same `-λ1h` shift into the θ2 factor too (which
should carry `-λ2h`, site 2's own weight) — that mismatch produced exactly the spurious `μSite2`.
**Fix:** a new `ν1[λ1,λ2][u] := (u-θ1-hλ1)(u-θ2-hλ2)` (the Yangian weight eigenvalue of T11 on the
HWS; deduped with Section B's `aFun`, same formula) replaces `Qθ` in `x2`'s denominator, argument no
longer carrying the `-λh` shift separately. With the fix, the factorization ratio is **exactly 1**
for all 27 `(λ1,λ2,s1)` cases — proved analytically too: `ν1[λ1,λ2][θ1+(k-1)h]` factors as
`((k-1-λ1)h)·(θ1-θ2+(k-1-λ2)h)`, and that second factor is exactly the old spurious `μSite2`. `xL1`
needed no functional change (`Qθ1` is linear, fix is a numeric no-op) but was rewritten with a
single-site `ν1[λ][u]:=(u-θ1-hλ)` for parity. **Section F confirmed unaffected** (`dictReportDet`
identical before/after). Applied to both `Clean/XXX_CG_L2_Clean.wb` and
`Experiments/XXX_CG_L2_V2.wb`; fresh-kernel Run-All clean on both after the fix (also caught and
fixed one unrelated hidden-dependency bug in V2's `obstructionReport` cell, same class as the
`sts12` issue during Clean promotion — a scratch-only `ABTab` variable that was never a saved cell).

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
`Clean/XXX_CG_L2_Clean.wb` is now the authoritative A–F notebook: A (SoV construction, incl. `ν1`) ·
B (Baxter/TQ, Psi/Q1) · C (2-site FSoV scalar product/orthogonality) · D (left eigenvectors,
`<Psi|Psi>=det[Q,Q]`) · E (exact L=1 factorization, c=1) · F (functional CG-overlap dictionary).
Companion `.md`/`.tex` reference docs are up to date with it (refreshed same day for the Section E
correction; `.tex` compiles clean with plain `pdflatex`).
