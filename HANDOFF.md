# Hand-off

Where the last session left off. Read this first when picking the project back up; it's not a
history log (that's `git log` / `logbook.tex`) — only what's still open or in-flight. Overwrite
freely each session rather than appending.

## Current focus

Computing CG coefficients directly in the XXX spin chain (`Paul/Mathematica/XXX/`). The Gaudin
subproject proper (residue/representation groundwork) is still parked, but the `h→0` scaling limit
of the XXX construction below is an active connection point, not a revival of that parked track.

## State as of 2026-07-22 (C++ CGC reference benchmark set up — `Paul/Cpp/AlexCGC/`, UNCOMMITTED)

Separate thread from the Baxter/Gaudin work below: set up the **state-of-the-art SU(N) CGC code**
(Alex, Kalus, Huckleberry, von Delft, arXiv:1009.0437 — the Gelfand–Tsetlin pattern-calculus C++
program from that paper's appendix) as a **reference benchmark** against Paul's Q-function CGC
generating function (`Paul/Mathematica/Gaudin/Clean/CGGenFun.wb`). New directory `Paul/Cpp/AlexCGC/`
— a C++ cross-check sibling to `Paul/Python/`. **Everything uncommitted** (Paul's choice: place,
don't commit; a compiled binary is machine-specific — rebuild from source elsewhere).

**Toolchain:** needs LAPACK/BLAS. `sudo apt-get install liblapack-dev libblas-dev` was run this
session (installed, 3.10.0). Build all binaries with `bash Paul/Cpp/AlexCGC/build.sh` (plain g++).

**What's there:**
- `ClebschGordan.cpp` → `clebsch` — the authors' original source (downloaded, unmodified); an
  interactive menu calculator (menu 5/6 = CGCs for S⊗S′). `hwcgc.sh` filters its output to just the
  highest-weight state (always the largest index `Q(M″)=dim(S″)`).
- `cgc_hw_bench.cpp` → `cgc_hw_bench` — modified copy that times ONLY the highest-weight-state solve
  (`compute_highest_weight_coeffs`, the raising-op nullspace/SVD) and SKIPS building descendants
  (3 clearly-commented `// added` edits: a `highest_weight_only` ctor flag, a scope timer, a new
  benchmark `main`). Verified `max|hw_only − full| = 0` vs the full run in every case. Modes:
  `cgc_hw_bench N S… S′… S″… [reps] [verify]` (per-irrep timing); `… decompose N S… S′…` (fast μ
  list); `machine` flag = one parseable RESULT line.
- `test.sh S [reps]` — for SU(4), times the highest-weight CGCs of EVERY μ in `[S,S,0,0]⊗[S,S,0,0]`
  (summary table only, no coefficient dumps). `reps` default 51 (report min-over-reps, for the fast
  small cases); use `reps=1` for large S (each solve already runs for seconds).
- `results/` — `hw_timings.txt` (full per-μ tables S=1..8), `summary.csv`, `gen_results.sh`
  (regenerates both). Machine: i5-1345U, 16 GB.

**The benchmark result (the point of all this):** C++ cost is dominated by the SINGLET (0,0,0,0)
— its highest-weight linear system is the largest — and total-over-all-μ hw-solve grows ~×5 per
step in S: S4=0.057, S5=0.355, S6=2.49, S7=15.3, S8=80.2 s. It also has a **memory wall**: the
singlet's dense SVD matrix is ≈ O(dim(ν)³) ≈ tens of GB by S=10 (dim[10,10,0,0]=1716), so the C++
**cannot even run S=10** on a 16 GB laptop. Paul's GF (build-for-all-irreps time) is ~**flat**:
S4=0.063, S6=0.139, S8=0.171, **S10=0.35 s**. ⇒ crossover ~S=4; GF ~470× faster at S=8, and at
S=10 it's "0.35 s vs won't-run" — a scaling win in **both time and memory**, not a constant factor.

**Open / next stage (Paul):** the GF timing is currently **build-only** — the explicit CGCs still
need to be EXTRACTED from the generating function (Paul's stated next implementation stage). That is
the one piece needed to make the head-to-head deliverable identical to the C++ (which outputs
explicit coefficients); then drop GF per-μ times into a second CSV beside the C++ `results/` for a
direct side-by-side. See memory `project_gf_cgc_vs_cpp_benchmark`.

## State as of 2026-07-22 (gl(4) L=2 untwisted τ1,τ2,τ3 ALL DONE exactly via Baxter asymptotics; Gaudin h→0 limit + standalone `GaudinSolver.wb` — DONE)

New notebook `Paul/Mathematica/XXX/Experiments/baxter_SU4_L2_XXX.wb` (Paul's own fresh build,
independent of `Baxter_L2_XXX_SU4.wb`) closes out the τ2 gap left open in the 2026-07-20 entry
below, via a genuinely different method — large-`u` asymptotic expansion of the untwisted Baxter
(TQ) equation, rather than operator-level Casimir fits. **This supersedes that entry's "τ2 —
PARTIAL" status: all three transfer matrices are now closed-form in Cas2(μ), Cas4(μ).**

### The mechanism (Paul's idea, worked through live)

- `n_k = μ_k` on the HWS (mode number = weight component) — verified directly against the main
  notebook's `nvec` inversion.
- Untwisted (`z_k→1`) Q-degree formula: `deg q_k = (S+P+3) − μ̂_k = (S+P)−μ_k+(k−1)`, `μ̂_k=μ_k+4−k`
  the *same* shifted weights that appear in the paper's Casimir generating function (eq 2.13–2.15,
  `G4[u]=Prod_i(u−h(μ_i+4−i))`). The `+(k−1)` staggering is forced, not just plausible: the
  twisted-case `Mdeg` formula's denominator is (up to `z³`) the derivative of the twist
  characteristic polynomial, which has a 4-fold zero at `z=1` — L'Hôpital-resolving that
  degeneracy is exactly what produces the staggered degrees.
- Key structural fact: `Σ_j(−1)^j C(4,j) j^m = 0` for `m=0,1,2,3` (a 4th finite difference) ⇒
  large-`u` orders `1/u,1/u²,1/u³` of the Baxter equation are **state-blind** (independent of the
  Q-data), and since the 4 Q-branches share one transfer matrix but have distinct leading powers
  `p_k`, every power of `p` must vanish separately. Solving that system pins 5 of τ1/τ2/τ3's 7
  unknown `u`-coefficients (`a[1],a[2],b[1],c[2],b[3]`) in terms of `S,P,θ` and one remaining
  scalar `b[2]` — leaving only `b[2]` (τ2's `u²` coeff, carries Cas2) and `d[2]` (τ2's constant,
  carries Cas4) as genuine Casimir-dependent data.
- Pseudovacuum anchor (`μ_top=[S+P,S+P,0,0]`, all τ's = elementary-symmetric-in-ν vacuum
  eigenvalues) gives `b2vac,c2vac,d2vac` independently — **cross-check passed exactly** (vacuum
  `c[2]` matches the asymptotically-predicted `c[2]`, residual 0).

### Result (τ1, τ2, τ3 all closed, all anchored at ΔCas = Cas(μ)−Cas(μ_top))

```
τ1(u) = τ1vac(u) + (h²/2)ΔCas2
τ2(u) = τ2vac(u) + (h²/8)ΔCas2·[τ1vac(u)+τ1vac(u−h)] + h⁴(⅛ΔCas2²−¼ΔCas4) + h⁴(1+2S+2P)/4·ΔCas2
τ3(u) = trivfactor(u)·[ 4u²−4u(θ1+θ2)+4θ1θ2−2hu(S+P+4) + b3vac + (h²/2)ΔCas2 ]
```

**The bug and how it was actually caught:** the first τ2 guess (fusion form
`(h²/8)ΔCas2·[τ1vac(u)+τ1vac(u−h)]` matching u²,u¹ exactly, plus the "obviously right" quartic
`h⁴(⅛ΔCas2²−¼ΔCas4)`) was missing the last `h⁴(1+2S+2P)/4·ΔCas2` term. Found by building the actual
polynomial Baxter solver and testing it: for `μ={2,1,1,0}@(1,1)` (`deg q1=0`), a constant `q1`
solved the equation at `μ_top` but **not** off-top — residual `=5·ν1(u−2h)ν1(u−3h)` exactly (a
pure constant shift), which pinned the missing term via the actual solvability condition, and
cross-checked against 27 irreps symbolically in `θ,h` (all exact) plus 4 more irreps with
`deg q1≥1` to break a degeneracy where `Cas2²` and `Cas4` weren't separable. **Lesson: build the
solver early — it's also the most direct correctness check.**

### Baxter polynomial solver + gauge fixing (in `baxter_SU4_L2_XXX.wb`)

`QsolveG[μ,S,P]` (note: `(μ,S,P)` order in *this* notebook — see naming-convention note below)
solves the untwisted 4th-order Baxter (difference) equation for the four monic Q-polynomials,
**gauge-fixed**: since `deg q1<...<deg q4` strictly, any lower-degree solution can be added to a
higher one (`q_k → q_k+c·q_j`, `j<k`); removed by forcing the coefficient of `u^(deg q_j)` in `q_k`
to vanish for every `j<k`. Verified exact (symbolic in `θ,h`) over full decompositions of every
`(S,P)` up to `(3,2)`, and numerically (rational `θ1=1/3,θ2=1/7,h=1`) over **all 21 irreps at (5,5)**
and **all 45 irreps at (8,8)** — every `q_k` solves the equation exactly, timing flat at
~0.02–0.05s/irrep **independent of irrep dimension** (largest block tested: `dim=69825`, same cost
as `dim=1`) since the solver only ever touches degree-≤19 polynomials, never the representation
itself.

Also built (Schur-polynomial machinery, reused verbatim in `GaudinSolver.wb` below): `irreps[S,P]`
— full `[S,S,0,0]⊗[P,P,0,0]` tensor-product decomposition via `s_λ·s_μ=Σ c^ν_{λμ} s_ν`, peeling off
the lex-highest monomial each step. Confirmed multiplicity-free throughout, dimension sums exact
(`Σdim(μ)=dim([S,S,0,0])²`), box-count (`Total[μ]=2S+2P`) conserved on every irrep as expected.

### Gaudin h→0 limit (derived in `baxter_SU4_L2_XXX.wb`, then built standalone)

Writing `Q=Q0+h·Q1+...`, orders `h⁰..h³` of the cleared Baxter equation vanish identically (same
alternating-binomial identity as above); order `h⁴` gives a **closed 4th-order ODE on Q0 alone**
(no `Q1,Q2` needed) — the gl(4) 2-site Gaudin oper, entirely in `Cas2(μ),Cas4(μ),S,P,θ1,θ2`.
Leading coefficient `8(u−θ1)²(u−θ2)²` (double singular points at the sites, as expected).
Validated: the ODE's polynomial-solution degrees match `degQ[μ,S,P]` exactly (and *only* those
degrees admit a solution) across every irrep of every `(S,P)` tested, symmetric and asymmetric.

**New standalone notebook `Paul/Mathematica/XXX/Experiments/GaudinSolver.wb`** — deliberately
minimal, gl(4) irrep/Casimir/degree machinery + the hardcoded `gaudinODE` only, **no** XXX transfer
matrices or Baxter difference-equation machinery (by design, per Paul's request). Built across
several rounds of naming/convention changes (all re-verified after each):
- Sites renamed `θ1,θ2 → w1,w2`; spectral variable renamed `u → z`.
- **Argument convention, applied everywhere: `(S,P,μ)`** — `degQ`, `GaudinSolve`, `BuildQs`,
  `verifyGaudin`, `CheckAllIrreps` all take `S,P` first, `μ` (or the irrep list) last. (Note: this
  is the *opposite* order from `QsolveG[μ,S,P]` in `baxter_SU4_L2_XXX.wb` — the two notebooks are
  not argument-compatible; don't copy-paste calls between them without swapping.)
- `GaudinSolve[S,P,μ]` returns **monic coefficient lists** (the practical/canonical form); `BuildQs
  [S,P,μ]` is a thin wrapper returning the same four Q's as **explicit polynomials in z** (kept
  purely for display/sharing, never used for internal verification/timing).
- `CheckAllIrreps[S,P,irrList]` takes a **precomputed** irrep list rather than recomputing
  `irreps[S,P]` internally — `irreps[8,8]` alone costs ~1.6s vs ~0.13s to actually solve all 45
  irreps, so precompute-once-reuse is the efficient pattern. Verified at `(8,8)` (45 irreps) and
  `(3,2)` (6 irreps).
- **Verification convention, hard-won:** always reconstruct a coefficient list into a testable
  function via `qf[x_]:=Sum[cf[[j+1]] x^j,{j,0,m}]` (SetDelayed), never
  `Dot[coeffs,z^Range[...]]` wrapped in a pure `Function` — the latter has a genuine WL edge case
  where symbolic `Derivative` of a length-1 (`deg=0`) Dot-built function spuriously returns a
  wrapped `{0}` instead of `0`, causing false verification failures. Cost a full debugging round
  before being isolated.

**Assistant arithmetic error, corrected:** an earlier summary wrongly claimed
`Σdim(μ)=8721²=76,055,441` for the `(8,8)` sweep — `8721` is `dimGL4[{16,16,0,0}]` (the *top*
output irrep), not `dimGL4[{8,8,0,0}]` (the actual input irrep, `=825`). Correct identity:
`Σdim(μ)=825²=680,625` — verified directly, matches the sweep. The sweep/solver results themselves
were never wrong, only that one squared-dimension claim (caught by Paul).

### Next natural steps

- **τ2 is now closed** — this directly unblocks the CG-overlap program the whole `z→1` limit was
  meant to feed (see the 2026-07-19/20 entries below for the original motivation). Natural next
  step: extract CG coefficients using the now-complete Casimir-form transfer matrices.
- Promote `baxter_SU4_L2_XXX.wb`'s Casimir-form τ1,τ2,τ3 + solver to `Clean/` once Paul is
  satisfied (still in `Experiments/`).
- `GaudinSolver.wb` is a clean, minimal, independently-reusable artifact — a natural next step is
  extending it toward the Gaudin-side CG-overlap analog, or documenting it via
  `wolfbook-docs-agent` once/if it's promoted to `Clean/`.
- Kernel-state loss happened repeatedly this session (multiple `KERNEL RESTART` events between
  turns) — always re-run cells top-to-bottom (`?symbolName` to check `DownValues` first) rather
  than trust that a previous turn's definitions are still live.

## State as of 2026-07-20 (gl(4) L=2 untwisted Casimir eigenvalues — τ1,τ3 DONE exactly, τ2 partial, in Experiments)

Follow-on to the twisted build below: new section appended to `Baxter_L2_XXX_SU4.wb` (now 70 cells,
35 onward), `z[k]→1` (untwisted/physical point), for the two-row-rectangular rep family
`λ1=[S,S,0,0]`, `λ2=[P,P,0,0]`. Goal: closed-form eigenvalues `Λ_a(u)` of the untwisted transfer
matrices `τ1,τ2,τ3` via su(4) Casimir invariants of the irreps `μ` in `λ1⊗λ2=⊕N_μμ`. Design/plan:
`docs/superpowers/specs/2026-07-19-su4-l2-untwisted-casimir-eigenvalues-design.md`,
`docs/superpowers/plans/2026-07-19-su4-l2-untwisted-casimir-eigenvalues.md` (the plan's original 5
tasks are DONE; everything below is follow-on work past the plan, done inline/live in the notebook
in direct collaboration with Paul, not via a formal plan). Subagents lacked wolfbook tools this
session (again) — a wolfbook-builder dispatch for the τ2/τ3 investigation partially worked (tools
were attached that one time) and its findings were independently re-verified from primitives.

**Exact-arithmetic simplification:** at `z=1` (vs. the twisted case's irrational `Zeta`-based twists)
every quantity is exact rational — no precision-cascade concerns anywhere in this section.

**Rep family is fully multiplicity-free** across the whole `(S,P)∈{(1,1),(2,1),(1,2),(2,2),(3,1)}`
sweep tested (up to `dim=400`), so every `τ_a` acts as a scalar on each irrep block (Schur's lemma) —
this is what makes closed forms in terms of Casimirs meaningful here at all.

### τ1 — DONE exactly

`Λ_1(u) = 4(u−θ1)(u−θ2) − 2P·h·(u−θ1) − 2S·h·(u−θ2) + (h²/2)(C2(μ) − C2(λ1) − C2(λ2))` (cell 48),
verified over the full sweep, `worst ~ 1e-25`. All `μ`-dependence enters only through the quadratic
Casimir eigenvalue `C2(μ)`, affine with slope `h²/2`. (This *resolved* an earlier false negative: a
naive fit of "each `u^j` coeff of `Λ_1` = `A_j+B_j·C2(μ)` **pooled across all `(S,P)`**" failed for
`j=0,1` — that was a pooling artifact, not a real falsification; the correct per-`(S,P)` formula
above is exact. See cells 45–47 for the superseded pooled-fit attempt, kept for the record with an
updated conclusion cell.)

### τ3 — DONE exactly

`Λ_3(u) = trivfactor[S,P]·( s3scalar[S,P](u)·Id + (h²/2)(C2(μ)−C2(λ1)−C2(λ2)) )` (cell 65), where
`trivfactor[S,P] = (u−θ1)(u−θ2)(u−θ1−h(S+2))(u−θ2−h(P+2))` (Paul's own factorization of τ3's four
trivial zeros, cells 58–63) and `s3scalar[S,P](u)` is an explicit state-independent degree-2
polynomial. **The `h²` term is exactly the same quadratic-Casimir combination as τ1's** — verified
as an exact operator identity on `(1,1),(2,1),(1,2)` including the `S↔P`-asymmetric case. So τ3
carries no higher Casimir, only `C2` — the `Λ³ℂ⁴≅(Λ¹ℂ⁴)*⊗det` self-duality made concrete.

### τ2 — PARTIAL (the self-dual middle case, genuinely harder)

Writing `Λ_2(u) = Σ_k[α_k(S,P)+β_k(S,P)·C2(μ)]u^k`, with `u⁰` additionally carrying
`+⅛C2(μ)²−¼C4(μ)` (universal coefficients, verified): **`α4=6,β4=0` (leading); `α3(S,P)=
−6h(S+P+2)−12(θ1+θ2), β3=0`; `β2=h²` (universal constant); `α2(S,P)` fully resolved (messier
polynomial, exact); `β1(S,P)=−(h³/2)(S+P+2)−h²(θ1+θ2)`** — all verified exactly against `tgCoeffs`
across 5 `(S,P)` points via direct operator-level matching (`M_k = α_k·Id+β_k·C2op`, solved exactly,
not curve-fit) — see cells 49–57 (the original agent-built investigation, confirming `τ2` uniquely
reaches `C4` only in its constant term) plus a from-scratch re-derivation in the live session after a
kernel restart (search for `solveABsym` in the transcript/notebook history).

**Still open: `α0(S,P)` and `α1(S,P)`'s full `(S,P)`-dependence.** These need cubic (`α1`) and
quartic (`α0`) polynomials in `(S,P)` — more unknowns than the 5 available `(S,P)` data points can
resolve. Getting more points is expensive (`(4,1)`→dim 630, `(3,2)`→dim 1000, symbolic-in-`θ,h`
computation at that size is impractical). An operator-level derivation was attempted instead: `T[i,j]`
was expanded exactly in terms of `Ee1,Ee2` building blocks (verified to reproduce `T` and `tg[2,1]`
bit-for-bit — see below) with the idea of reading off `α1,α0` structurally rather than by curve-fitting,
but this did not pay off within the session (a residual-decomposition test against `{X1,X2,Y1,Y2,C2}`
at a single `(S,P)` is underdetermined and not actually informative — would need to be done
simultaneously across multiple `(S,P)` to mean anything). **Promising unexplored lead, added by Paul
independently in cells 58–61**: `t1vac,t2vac,t3vac` — explicit vacuum-eigenvalue formulas built from
`L1=L2=(u−θ1−hS)(u−θ2−hP)`, `L3=L4=(u−θ1)(u−θ2)` — `t2vac` sums all 6 pairs `Li(u)Lj(u−h)` in a
specific descending-index order. An earlier attempt this session to derive the τ2 vacuum eigenvalue
via `Σ_{i<j}ν_i(u)ν_j(u−h)` (ascending-index order) did **not** match `tg[2,1]` on the HWS (nonzero
residual on pairs like `{1,3}`) — Paul's `t2vac` uses the opposite index order, which may be exactly
the fix; reconciling this against the actual `qmin`/`tg[2,1]` construction is the natural next step
and could resolve `α1,α0` cheaply (vacuum eigenvalues need no large-matrix computation at all).

### General Casimir operator formulas (not just eigenvalues — the operators themselves)

Cells 65–68 (note: cell 67, `X1/X2/Y1/Y2` definitions, was a **bugfix** — those operators were
originally only ever created via scratch `evaluateExpression` calls and never actually saved to a
notebook cell; a fresh-kernel run would have failed before the fix). General master formula (any `n`,
verified `n=2,3,4`): `Cn_op[λ1,λ2] = Σ_{pat∈{1,2}^n} patternTerm[pat,λ1,λ2]`, where each of the `2ⁿ`
site-assignment patterns collapses (via `(A⊗B)(C⊗D)=(AC)⊗(BD)` applied repeatedly) to `[ordered
product of site-1 factors]⊗[ordered product of site-2 factors]`, cyclic indices tying both chains
together. Reduced explicit forms with named cross terms:
```
C2op = (Csite[λ1]+Csite[λ2])Id + 2·Omega2                                    [1 cross term]
C3op = (C3site[λ1]+C3site[λ2])Id + 2X1+2X2+Y1+Y2                              [4 cross terms]
C4op = (C4site[λ1]+C4site[λ2])Id + 2Z1+2Z2+2Z3+2Z4+Z5+2Z6+Z7+2Z8               [8 cross terms]
```
Cross-term-type count grows 1→4→8 for `n=2,3,4` (matrix products aren't cyclic, only traces are, so
patterns that look like rotations of each other are usually genuinely different operators) — itself a
notable structural fact. Every formula verified as an **exact** operator identity (`==`, not
tolerance) on asymmetric `λ1≠λ2` pairs (avoiding a `λ1=λ2` coincidence hiding real distinctions), plus
`C3op`'s master formula also checked outside the `[S,S,0,0]` family (adjoint⊗fund) to confirm
generality. Plain operator definition (cell 70, self-contained): `C2=Σ_ab E_ab.E_ba`,
`C3=Σ_abc E_ab.E_bc.E_ca`, `C4=Σ_abcd E_ab.E_bc.E_cd.E_da`, `E_ab:=EE[λ1,λ2][a,b]`.

### General Casimir eigenvalue formulas (cell 69, self-contained, no other cell needed)

Single-site rep `[x,x,0,0]`: `Cₙ(x)=2x(x+2)ⁿ⁻¹` — `C2c[x]=2x(x+2)`, `C3c[x]=2x(x+2)²`,
`C4c[x]=2x(x+2)³`, verified `x=1..5` against the operators. Combined two-site:
`C2[S,P]=2S(S+2)+2P(P+2)`, `C3[S,P]=2S(S+2)²+2P(P+2)²`, `C4[S,P]=2S(S+2)³+2P(P+2)³`.

**Real bugs found and fixed this session** (beyond the earlier three): `X1/X2/Y1/Y2` never persisted
to a cell (see above); `T[i,j]` off-diagonal action on the HWS genuinely does vanish both directions
(`T[i,j][u].HWS=0` for `i≠j`, confirmed identically in symbolic `u`) — useful building block, not a
bug, but non-obvious and worth remembering; a hand-derived vacuum-formula attempt for τ2 had an
index/shift-order bug not yet isolated (see τ2 section above).

**Next natural steps:** reconcile Paul's `t2vac` (cell 61) against `tg[2,1]` on the HWS to get τ2's
vacuum eigenvalue cheaply, then use it (plus the already-known `β1,β0` and the operator-level `T[i,j]`
expansion) to pin down `α1(S,P),α0(S,P)` without more expensive `(S,P)` sweeps. Once τ2 is fully
closed-form, the natural follow-on is the CG-overlap program this whole `z→1` limit was meant to feed
into (see `HANDOFF.md`'s `su4 L2 untwisted` note in the twisted-case section below for the original
motivation).

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
