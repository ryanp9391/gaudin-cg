# Hand-off

Where the last session left off. Read this first when picking the project back up; it's not a
history log (that's `git log` / `logbook.tex`) — only what's still open or in-flight. Overwrite
freely each session rather than appending.

## Current focus

Computing CG coefficients directly in the XXX spin chain (`Paul/Mathematica/XXX/`). The Gaudin
subproject is parked.

## State as of 2026-07-15

- Repo just got split into two subprojects: `Paul/Mathematica/XXX/` and `Paul/Mathematica/Gaudin/`
  (each with `Clean/`/`Experiments/`), and the paper into `paper/xxxCG.tex` (has the real XXX
  content) and `paper/gaudinCG.tex` (fresh stub). See `CLAUDE.md` for the full layout.
- **Nothing from this restructuring is committed yet** — it's all staged/untracked in the working
  tree. Check `git status` before doing anything else; don't assume HEAD reflects the current
  layout. This covers everything below too (the two Clean promotions, agent consolidation,
  settings change, and the new `XXX_CG_L3.wb` work) — none of it is committed either. The only
  things committed this session are the L3 design spec and implementation plan under
  `docs/superpowers/`.
- Only one wolfbook subagent now: `wolfbook-agent` was removed (was redundant with
  `wolfbook-builder` — both ran on `sonnet`, so the split wasn't buying anything). Use
  `wolfbook-builder` for all live-notebook work, routine or demanding. `CLAUDE.md` and the
  `reference_gaudin_cg_agents` memory are updated to match.
- `.claude/settings.json` permissions broadened: added a `mcp__wolfbook__*` wildcard plus broad
  `Read`/`Grep`/`Glob`, so wolfbook-builder shouldn't need per-call approval anymore. If prompts
  persist, the likely cause is the settings watcher not having picked up the mid-session edit —
  try `/hooks` once or restart before assuming the rule itself is wrong.
- `Paul/Mathematica/XXX/Clean/` now has two built-and-documented notebooks: `XXX_CG_L2_Clean.wb`
  (from the prior session) and `XXX_CG_L3.wb` (new this session — see below).
- New subagent `wolfbook-docs-agent` (`.claude/agents/wolfbook-docs-agent.md`) added: read-only
  over `Clean/*.wb`, writes a companion `<NotebookName>.md` (exhaustive per-symbol reference) AND
  `<NotebookName>.tex` (an honest paper-style LaTeX document — prose organized by mathematical
  topic, real display equations, not a mechanical transliteration of the `.md`). Never touches
  `paper/*.tex` itself, never commits. `CLAUDE.md` and the `reference_gaudin_cg_agents` memory are
  updated to match. Two corrections already folded into the agent definition from the prior
  session: (1) it must not re-verify Clean/ equations — they're user-verified before promotion, so
  it just documents; (2) Greek notebook variable names (`λ`,`θ`,`κ`,`χ`,`σ`) must become LaTeX
  macros (`\lambda` etc.), never literal Unicode glyphs, and the `.tex` must compile with plain
  `pdflatex`, not `xelatex`.
- Gotcha for any future LaTeX in this repo, not just the docs agent: identity matrices/operators
  must use `\mathbbm{1}` (`bbm` package), not plain `amssymb`'s `\mathbb{1}` — the latter compiles
  with zero errors but renders as a garbled substituted glyph, not a clean "1". Caught by visually
  inspecting the compiled PDF, not by the compile log. Already fixed in both `XXX_CG_L2_Clean.tex`
  and `XXX_CG_L3.tex`, and folded into the docs-agent definition/memory.
- Gotcha for anyone compiling `.tex` files by hand in `Clean/`: `pdflatex` leaves `.aux`/`.log`/
  `.out`/`.fls`/`.fdb_latexmk`/`.pdf` byproducts alongside the source — clean these up (or gitignore
  them) before committing; they were left behind once this session and manually removed.
- Gotcha specific to the `wolfbook` VS Code extension: it appears to share **one kernel across all
  open notebooks** in a window, not one kernel per notebook — switching target notebooks may
  require an explicit `wolfbook_setTarget`/`switch` call (`wolfbook_getNotebookContext` can report
  "No active notebook editor" otherwise), and a kernel restart may be needed if stale definitions
  from a different notebook leak in. Also: `wolfbook_editCell`/`runCell`'s persisted output can
  occasionally go stale relative to what was actually evaluated (a cell's *source* update saved
  correctly while its *displayed output* on disk still showed a prior run) — if a result looks
  suspicious, re-run the cell fresh and re-save rather than trusting the last-saved output blindly.

## Open threads / next steps

- Decide whether/when to commit the restructuring (and everything since, including `XXX_CG_L3.wb`
  and its docs).
- `Paul/Mathematica/XXX/Experiments/XXX_su2.wb`, `su2_XXX_L2.wb`, and `Baxter_Solver_L2.wb` are
  untouched and stay in `Experiments/` as the historical record behind `XXX_CG_L2_Clean.wb`
  (confirmed convention — promoting to `Clean/` copies content, it doesn't empty the source).
  `XXX_CG_L3.wb`, by contrast, was built directly in `Clean/` from scratch this session (no
  `Experiments/` predecessor) — a deliberate one-off deviation from that convention, per explicit
  user decision recorded in the design spec below.
- `paper/gaudinCG.tex` is just a stub; no Gaudin-specific derivation has been ported over or
  rewritten since the split.
- Natural next step: start actually building CG coefficients on top of `XXX_CG_L2_Clean.wb` /
  `XXX_CG_L3.wb`, or continue the L-generalization pattern further if needed. Nothing queued.

### Done: `Paul/Mathematica/XXX/Clean/XXX_CG_L3.wb` built, debugged, and documented (completed 2026-07-15)

Three-site (L=3) generalization of `XXX_CG_L2_Clean.wb`, built directly in `Clean/` (no
`Experiments/` merge this time) via `docs/superpowers/specs/2026-07-15-xxx-cg-l3-design.md` →
`docs/superpowers/plans/2026-07-15-xxx-cg-l3-notebook.md` (both committed) → `wolfbook-builder`
executing the plan. 44 cells, fresh-kernel top-to-bottom run clean, end-to-end Baxter TQ residual
check exact zero at 18–23 tracked digits for both $(\lambda_1,\lambda_2,\lambda_3)=(1,1,1)$ and
$(2,1,1)$.

**Key physics difference from L2, worth knowing before touching this again:** $t_{1,1}(u)$ is
*cubic* in $u$ at L=3 (three Lax factors), not quadratic. This means `TauEigensystem` needs a
**two**-parameter joint diagonalization ($C_0+r_1C_1+r_2C_2$, not L2's one-parameter
$C_0+rC_1$) — $M$ is read off $C_2$'s expectation (the pure Cartan/total-generator coefficient,
proven by an asymptotic cross-check), while $(J_1,J_0)$ (not L2's single $J$) are recorded
directly per eigenvector since $C_1,C_0$ genuinely mix different sites' raising/lowering operators
and aren't predictable in closed form. `τ(u)=τ0(u)+J1 u+J0`. See the plan document's "Key
derivation this plan encodes" section for the full argument if this needs extending to L=4+.

**Three bugs found and fixed during the build** (all in the plan's transcription, not the
physics — caught by `wolfbook-builder` stopping and reporting rather than patching around
failures, each independently re-derived/confirmed before applying the fix):
1. `BaxterMatrix`'s hard-assert originally checked 3 high-degree coefficients (M+1,M+2,M+3); only
   the top 2 (M+2,M+3) are pure per-column identities — M+1 is legitimately absorbed by `QSolve`'s
   `J1·shiftDown` term. General rule: for any chain length L, exactly the top two extra-degree
   coefficients are pure identities, regardless of L.
2. `Psi`/`τ`/`QSolve`'s bounds check inherited L2's `n<=M` clause, which relied on an L2-only
   coincidence (`dPred<=M+1` there always). L3's `dPred` is a genuine triple convolution and can
   exceed `M+1` (e.g. `dPred[1,1,1,1]=3`) — dropped to `n>=0`.
3. `checkTQ`'s outer sampling `Table` used bare `u`, colliding with `BaxterMatrix`'s internal
   non-Module-localized `u` (Table's dynamic/Block-style iterator binding clobbers it on any
   first-time uncached `QSolve` call inside the loop) — renamed to `uu`, matching the convention
   the L2 predecessor (`Experiments/Baxter_Solver_L2.wb`) already used for this exact reason.

All three fixes are documented inline in the plan document (kept in sync, committed) as well as in
the notebook's own cell comments.

**One more bug found post-build, by the user manually running the notebook, fixed directly (not
via the plan/agent since it was a small isolated cell):** Cell 32 (a numeric verification cell the
builder added beyond the plan, to work around a benign Simplify-doesn't-collapse-to-`{0}`
precision-display artifact) had one stray extra `]` before `//Chop` — a syntax error, not a math
bug. Fixed, re-run, confirmed `0`. Also found and fixed in the same pass: Cell 18 (the analogous
Task 4 check) had been inserted but never actually run; the end-to-end TQ residual cells (43–44)
had reverted to "not evaluated" (likely a side effect of the user's own manual cell insertion
shifting numbering) — both re-run and reconfirmed clean.

**Docs:** `wolfbook-docs-agent` generated `XXX_CG_L3.md`/`XXX_CG_L3.tex` following the same
conventions as the L2 docs. The agent's first pass (read-only) flagged Cell 32 as showing a stale
error output it couldn't resolve without eval access; that's now fixed in both doc files to
reflect the confirmed-clean state (the agent's flag was correct and appropriately cautious — it
just predated the final fix).

**Next step (if picking this up again):** none queued — this thread is closed unless you want to
extend the pattern to L=4+, or start building CG coefficients on top of these notebooks.

### Done: `XXX_CG_L2_Clean.md`/`.tex` reference docs written (completed 2026-07-15, prior session)

`wolfbook-docs-agent` generated both companion files for
`Paul/Mathematica/XXX/Clean/XXX_CG_L2_Clean.wb`. `Paul/Mathematica/Gaudin/Clean/FastQ.wb` and
`Hamiltonians_SU4.wb` still don't have `.md`/`.tex` docs (not requested; Gaudin is parked anyway,
but the agent works fine on it if asked).

### Done: `Paul/Mathematica/XXX/Clean/XXX_CG_L2_Clean.wb` built (completed 2026-07-15, prior session)

New Clean promotion, assembled by `wolfbook-builder` from two Experiments notebooks (both left
unchanged):
- **Section A** (cells 2–26): all 25 cells of `su2_XXX_L2.wb` copied verbatim, checks included.
- **Section B** (cells 27–34): only the 7 routine-defining cells from `Baxter_Solver_L2.wb`
  (a/d/τ0, BaxterMatrix, t2Coeffs, dPred+TauEigensystem, Psi/τ accessors, QSolve, Q1) — no
  validation/regression cells, and the machinery Baxter had ported verbatim from su2_XXX_L2 was
  deduped down to just the `χ1=κ1+κ2; χ2=1` numeric pinning (moved into the shared params cell),
  reusing Section A's `t2` etc. directly.

Verified: fresh-kernel top-to-bottom run of all 34 cells, zero errors; `Q1`/`QSolve` outputs
spot-checked against the original standalone `Baxter_Solver_L2.wb` and match exactly. Notebook
saved.

**One cosmetic loose end:** cell outputs may show empty on disk (a wolfbook tool-capture glitch
during the merge, confirmed unrelated to correctness via direct kernel-state checks) — a manual
"Run All" in VS Code will populate them if that matters to you; nothing computational is wrong.

### Done: reformulated `Baxter_Solver_L2.wb` (completed 2026-07-15, prior session, superseded above)

Replaced the old Baxter-matrix-diagonalize-then-filter approach with a two-stage flow: genuine
transfer-matrix eigenstates found first (`TauEigensystem`, with a hard residual genuineness check
and M/n read off algebraically from the τ eigenvalue's u¹/u⁰ coefficients), then the TQ equation
solved as a linear `NullSpace` problem (`QSolve`) for the Q-function. Both regression cases pass
(λ1=λ2=1; λ1=2,λ2=1 including (M,n)=(2,0)), residuals ~1e-20 to 1e-26.

One deviation worth remembering if this pattern comes up again: the specific NullSpace tolerance
bug this task was designed around (missed null direction at (λ1,λ2,M,n)=(2,1,2,0)) didn't
reproduce in this system — instead an analogous pitfall showed up in the M=0 sector (a matrix
entry that's exactly zero but with finite tracked precision, which `NullSpace` won't treat as
"confirmed zero" without an explicit `Chop` first). Fixed by `Chop`-ing before `NullSpace`.

This file (`Baxter_Solver_L2.wb`) is now itself superseded by `XXX_CG_L2_Clean.wb` above for
actual reuse, but stays in `Experiments/` unchanged as the derivation record.
