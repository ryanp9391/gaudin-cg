# gaudin-CG

**Check `HANDOFF.md` first** — it has where the previous session left off (in-flight state, open
threads, next steps). This file covers stable structure/conventions only.

Research repo (Paul Ryan). **Project target (since 2026-10-22): the slice engine and the two-site
Gaudin missing-label operators for gl(4)** — the PRV/Kostant slice theorem (string conditions on a
weight space of one site) plus the commuting Laurent coefficients of the two-site gl(4) Gaudin cdet
as the three missing labels, giving the highest-weight Clebsch–Gordan vectors of gl(4) without any
n_W-sized null space (plan: `Paul/Mathematica/Gaudin/gl4-slice-engine-plan.md`, code in
`Paul/Python/gaudin_gl4/`). This revives the repo's original gl(4)/su(4) interest. The gl(3) programme
(2026-10-16..21: CG coefficients of gl(3) from SoV in the two-site Gaudin model — the jet formula, its
general-Λ kernel, the untwisting limit, and finally the gl(3) slice engine + Jacobi chain that beat the
Munich code) is complete and is the template; its write-up is `Slice-engine-notes.*` and
`Gaudin-sl3-CG-notes.*`. The repo is split into two subprojects:

- **Gaudin** — **the current focus.** Plans per round, notebooks, reference data and the working
  notes live under `Paul/Mathematica/Gaudin/` (see Layout); the Python cross-checks and fits under
  `Paul/Python/gaudin_sl3/`.
- **XXX** — CG coefficients directly in the rational gl(3)/gl(4) (XXX) spin chain via SoV. This was
  the scaffold for the Gaudin work (the Gaudin formulas are its h → 0 limit); not actively extended.

Don't mix the two: a notebook or paper edit for one subproject shouldn't casually pull in the
other's conventions or formulas unless the task is explicitly about the scaling-limit connection
between them.

## Layout

- `paper/xxxCG.tex` — JHEP-format paper for the XXX subproject. Currently covers the Lax matrix,
  companion twist, fusion hierarchy, and SoV construction for the XXX spin chain. Source of truth
  for XXX derivations/formulas; check here before trusting a notebook ansatz.
- `paper/gaudinCG.tex` — JHEP-format paper for the Gaudin subproject. Currently a gl(4)-rectangular
  stub (title/author only); the sl(3) Gaudin CG material will replace its body (Paul's call on
  whether to rewrite it or start `paper/gl3GaudinCG.tex`). Until then the authoritative write-up is
  `Paul/Mathematica/Gaudin/Gaudin-sl3-CG-notes.{tex,md,pdf}` (generated from `notes-src/`).
- `Paul/Mathematica/Gaudin/` — the sl(3) Gaudin CG program: `Gaudin-sl3-*plan.md` (one per round,
  with the verdict blocks), `Experiments/Gaudin-*.wb` (chain-side notebooks + their `-spec.md`
  briefs; no `Clean/` yet), `ref/*.txt` (exact machine-readable reference data for the Python
  agents), the notes above. `Paul/Python/gaudin_sl3/` holds the gl(3) Python side: `SPEC-*.md` briefs,
  scripts, `RESULTS-*.md` reports beside each script, and `*.pkl` caches (gitignored, rebuildable).
  `Paul/Python/gaudin_gl4/` is the gl(4) programme (2026-10-22 on): one `SPEC-G<k>-*.md` brief per work
  package G1–G6 of the plan, the scripts (`glN_reps.py` — GT reps of gl(N), general N;
  `slice_engine_gl4.py`; `gaudin_gl4_ops.py` — the cdet operators), `RESULTS-G<k>-*.md` beside each;
  reference data in `Paul/Mathematica/Gaudin/ref/gl4-*.txt`.
  `Paul/Mathematica/old-conventions/old/Gaudin/` is the OLD gl(4)/su(4) Gaudin material, not this.
- `Paul/Mathematica/new-conventions/` — **active** notebooks going forward (2026-10-01 convention
  change). New conventions: Lax ℒ_ij(u) = u δ_ij + h E_ij (was u − h Σ E_ji⊗e_ij), companion
  twist acting on the right (TG = ℒ·G), B/C built from the untwisted T. Reference:
  `Companion-B-C-diagonalize-clean.wb` (gl(3), L=1, adjoint; final formulas + 37 passing checks,
  self-contained, clears the kernel) with its write-up `Companion-B-C-diagonalize-clean.md`/`.tex`/`.pdf`.
  `Companion-B-C-diagonalize.wb` is the working notebook (explorations); `Baxter-QQ-Wronskian.wb`
  runs the Baxter/QQ/Wronskian/orthogonality checks for any rep. L=2 (adj⊗adj, diagonal frame):
  `L2-diagonal-frame-overlap.wb` + `.md`. **Generator convention from end of 2026-10-01: Hermitian
  EH_ij = R E_ij R⁻¹ (EH_ijᵀ = EH_ji)**; `*-EH.wb` notebooks use it; exact zero tests there need
  `PossibleZeroQ` (surds). `Hermitian-structure-abstract.wb` + `.md`: abstract † computation. `paper/xxxCG.tex` still uses the OLD
  conventions until it is updated — don't mix formulas from the two without converting.
- `Paul/Mathematica/old-conventions/` — everything from before the 2026-10-01 convention change,
  structure preserved (paths below are relative to it). Reference only; port deliberately.
- `old-conventions/XXX-CG/` — the XXX notebooks from 2026-09-23 to 2026-09-30 (Building-SoV-bases,
  L1/L2 mathcalA, etc.; HANDOFF.md paths into `Paul/Mathematica/XXX-CG/` now live here).
- `old-conventions/old/` — archive of the pre-2026-09-23 Mathematica tree, kept for reference only
  (not actively maintained, not where new work goes).
  - `old/XXX/` — everything that used to be `Paul/Mathematica/XXX/`, structure preserved as-is:
    `Clean/`, `Experiments/` (the ~90-notebook pile), `rotations/` (self-contained thread on which
    covectors `⟨0|g` are determinants of Q-integrals — gl(3), L=2, adjoint ⊗ fundamental; start at
    its `README.md`).
  - `old/Gaudin/` — everything that used to be `Paul/Mathematica/Gaudin/`, structure preserved:
    `Clean/` (`FastQ.wb` residue formulas, `Hamiltonians_SU4.wb` gl(4)/su(4) GT-basis/generators),
    `Experiments/` (`ResidueBuilder.wb` and others).
- `Paul/Python/` — Python cross-check port, orthogonal to the XXX/Gaudin split (not a third
  subproject). Reimplements select Mathematica notebook results independently in Python
  (`sympy`/`uv`-managed) as a correctness cross-check, starting with the gl(3) single-site
  generators ported from `old/XXX/Experiments/su3_V2_general.wb` (`gl3_reps.py`, tested via
  `test_gl3_reps.py`). Design/plan docs under `docs/superpowers/specs/` and
  `docs/superpowers/plans/` (`*-gl3-python-port-*`).

Each `.wb` notebook has a matching `img/<NotebookName>/` directory (wolfbook eval logs, kernel
debug logs) — tooling artifacts, not research content.

## Conventions

- **Clean vs Experiments** (within `old/XXX/` and `old/Gaudin/`): only promote a notebook/
  definition to `Clean/` once it's validated. When extending a `Clean/` definition, check
  `Experiments/` first for the exploratory work behind it rather than re-deriving from scratch.
  Decide fresh whether `XXX-CG/` needs the same split once it has enough notebooks to warrant one.
- Notebook cells use `(* ... *)` comments explaining *why*, not restating the code.
- `paper/xxxCG.tex` / `paper/gaudinCG.tex` are the authoritative source for formulas in their
  respective subprojects; notebooks implement/check them, not the other way around.
- `Hamiltonians_SU4.wb` (su(4) GT-basis/generators) lives under `old/Gaudin/Clean/` even though the
  representation theory itself isn't Gaudin-specific — that's where it was placed during the
  original split, not a statement that it's unusable from XXX work. If XXX-CG work ends up needing
  it, ask before duplicating or relocating it.

## Subagents (`.claude/agents/`)

- **wolfbook-builder** — do live Wolfram computation directly in the `.wb` notebooks via wolfbook
  MCP tools, instead of writing a standalone `.wl` script. Covers everything from routine
  extensions to ports/generalizations and long derivation chains where a subtle error would be
  costly (e.g. L=1 to L=2). Stays out of `paper/*.tex` and git commits/pushes. Gaudin work goes in
  `Paul/Mathematica/Gaudin/Experiments/` (one `-spec.md` brief per notebook; shared kernel K1); XXX
  work in `Paul/Mathematica/new-conventions/`, never `old-conventions/`. Don't assume the other
  subproject's notebooks are relevant to a task scoped to one. Always dispatch on `model: sonnet`.
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
