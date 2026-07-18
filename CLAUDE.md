# gaudin-CG

**Check `HANDOFF.md` first** — it has where the previous session left off (in-flight state, open
threads, next steps). This file covers stable structure/conventions only.

Research repo (Paul Ryan) computing Clebsch-Gordan coefficients for su(4) rectangular
representations. The repo is split into two subprojects:

- **XXX** — computing CG coefficients directly in the rational gl(2)/gl(4) (XXX) spin chain via
  Separation of Variables (SoV). **This is the current focus.**
- **Gaudin** — the Gaudin-model version of the same problem, obtained as a scaling limit of the
  XXX spin chain. Parked for now; only the earlier residue-evaluation and su(4)-representation
  groundwork lives here.

Don't mix the two: a notebook or paper edit for one subproject shouldn't casually pull in the
other's conventions or formulas unless the task is explicitly about the scaling-limit connection
between them.

## Layout

- `paper/xxxCG.tex` — JHEP-format paper for the XXX subproject. Currently covers the Lax matrix,
  companion twist, fusion hierarchy, and SoV construction for the XXX spin chain. Source of truth
  for XXX derivations/formulas; check here before trusting a notebook ansatz.
- `paper/gaudinCG.tex` — JHEP-format paper for the Gaudin subproject. Currently a stub
  (title/author only) — no Gaudin-specific derivations written up yet.
- `logbook.tex` — single, shared dated research logbook (`DD-MM-YYYY` section headings) covering
  both subprojects, research content only (no tooling/admin notes). Maintained via the
  `logbook-agent` subagent, not edited by hand.
- `Paul/Mathematica/XXX/` — XXX subproject notebooks.
  - `Experiments/XXX_su2.wb` — single-site gl(N)/su(2) generator and GT-basis machinery (`ps`,
    `dim`, `J0`, `Jp`, `Jm`, `JJ`, `Ee`, Lax matrix `L`); building blocks reused at multiple sites.
  - `Experiments/su2_XXX_L2.wb` — two-site (L=2) XXX su(2) chain: independent gl(2) reps and
    inhomogeneities per site, shared companion twist. Ports the single-site machinery from
    `XXX_su2.wb` verbatim.
  - `Experiments/Baxter_Solver_L2.wb` — standalone, self-contained numerical solver + validation
    suite for the L=2 su(2) XXX Baxter (TQ) equation. The physics/derivation lives in
    `su2_XXX_L2.wb`; this file is just the reusable solver.
  - `Experiments/su3_V1.wb` — gl(3) analog of the su(2) Clean notebook's Section B: two-site
    eigensystem + three-Q Baxter solver for symmetric reps `{S,0,0}` per site.
  - `Experiments/su3_V2_general.wb` — generalizes `su3_V1.wb` from symmetric reps to arbitrary
    gl(3) Young-diagram reps `λ={λ1,λ2,λ3}`, one independent diagram per site. Full-sweep verified
    (102 states, worstTQ=worstCasoratian=0); next natural step is promotion to `Clean/`.
  - `Experiments/su3_V2_general_Paul.wb` — Paul's personal exploratory fork of
    `su3_V2_general.wb`; diverges from the shared notebook, not a drop-in replacement.
  - `Clean/` — empty so far (nothing promoted yet); validated XXX definitions go here once
    checked, same convention as `Gaudin/Clean/` below.
- `Paul/Mathematica/Gaudin/` — Gaudin subproject notebooks.
  - `Clean/FastQ.wb` — closed-form residue formulas (`res1`/`res0` via Laguerre polynomials) for
    `f(u) = Exp[w u] u^n / ((u-1)^S1 u^S2)`, plus `doRes0`/`doRes1` generalizing to arbitrary
    polynomials.
  - `Clean/Hamiltonians_SU4.wb` — gl(4)/su(4) irrep construction via Gelfand-Tsetlin patterns:
    `ps[rep]` (GT basis), `dim[rep]`, and generator matrix elements `J0`/`Jp`/`Jm`/`JJ`.
  - `Experiments/ResidueBuilder.wb` — exploratory companion to `Clean/FastQ.wb`: validation
    checks, benchmark comparisons of evaluation strategies, scaling/precision experiments that
    motivated what was promoted to `Clean/`.
- `Paul/Python/` — Python cross-check port, orthogonal to the XXX/Gaudin split (not a third
  subproject). Reimplements select Mathematica notebook results independently in Python
  (`sympy`/`uv`-managed) as a correctness cross-check, starting with the gl(3) single-site
  generators ported from `XXX/Experiments/su3_V2_general.wb` (`gl3_reps.py`, tested via
  `test_gl3_reps.py`). Design/plan docs under `docs/superpowers/specs/` and
  `docs/superpowers/plans/` (`*-gl3-python-port-*`).

Each `.wb` notebook has a matching `img/<NotebookName>/` directory (wolfbook eval logs, kernel
debug logs) — tooling artifacts, not research content.

## Conventions

- **Clean vs Experiments** (within each subproject): only promote a notebook/definition to
  `Clean/` once it's validated. When extending a `Clean/` definition, check `Experiments/` first
  for the exploratory work behind it rather than re-deriving from scratch.
- Notebook cells use `(* ... *)` comments explaining *why*, not restating the code.
- `paper/xxxCG.tex` / `paper/gaudinCG.tex` are the authoritative source for formulas in their
  respective subprojects; notebooks implement/check them, not the other way around.
- `Hamiltonians_SU4.wb` (su(4) GT-basis/generators) currently lives under `Gaudin/Clean/` even
  though the representation theory itself isn't Gaudin-specific — that's where it was placed
  during the split, not a statement that it's unusable from XXX work. If XXX work ends up needing
  it, ask before duplicating or relocating it.

## Subagents (`.claude/agents/`)

- **wolfbook-builder** — do live Wolfram computation directly in the `.wb` notebooks via wolfbook
  MCP tools, instead of writing a standalone `.wl` script. Covers everything from routine
  extensions to ports/generalizations and long derivation chains where a subtle error would be
  costly (e.g. L=1 to L=2). Stays out of `paper/*.tex` and git commits/pushes. When a task is
  scoped to one subproject, work only under `Paul/Mathematica/XXX/` or `Paul/Mathematica/Gaudin/`
  as appropriate — don't assume the other subproject's notebooks are relevant.
- **logbook-agent** — invoke with "update the logbook". Reads today's transcripts + git diff,
  writes a research-only entry into `logbook.tex`, compiles, commits, and pushes. The logbook
  stays unified across both subprojects (not split).
- **wolfbook-docs-agent** — invoke with "document the Clean notebooks" or similar. Read-only over
  `Clean/*.wb` (never `Experiments/`): generates/refreshes a companion `<NotebookName>.md` AND
  `<NotebookName>.tex` reference (symbol-by-symbol purpose/signature/dependencies/verification
  status) colocated with each Clean notebook, in both formats. Documents only — equations in
  `Clean/` are already user-verified before promotion, so it doesn't re-check the math. Never
  edits notebooks, never touches `paper/*.tex` (the actual paper, distinct from these reference
  docs), never commits to git — write-only within `Clean/*.md` and `Clean/*.tex`.

## Gotchas worth knowing before touching notebooks

- Multi-statement Wolfram cells need `;` between statements or the front end can silently merge
  lines into one `RuleDelayed`.
- Memoized definitions (`f[x_]:=f[x]=...`) cache bad results — after fixing a definition, re-run
  the defining cell and everything downstream that used the stale cache.
- Cell edits/evaluations don't persist until the notebook is explicitly saved.
