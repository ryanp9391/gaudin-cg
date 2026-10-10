# Execution plan: one-body SoV functionals and CG coefficients

Status: proposed assignments; no research tasks launched by writing this plan.

## 1. Objective and constraints

Complete the adjoint × adjoint, two-site gl(3) construction using one-body covectors and Q-functionals. Start from the established site-1 localization table, especially its singlet-reaching row 8. Deliver explicit normalized formulas, a justified route to CG coefficients, a bounded contour construction, and then the mirror and generalization checks.

Read `HANDOFF-SoV-bases.md` and `main-notes.tex` first. The handoff is a map; the rolling notes specify formulas and conventions. Read relevant notebook context through Wolfbook before proposing code based on its definitions.

Non-negotiable constraints:

- Do not introduce operators on the coupled 64-dimensional space or construct entangled covectors as a new solution method. Use single-site matrices, auxiliary-space identities, Q-data, and highest-weight coefficient relations. Existing recorded coupled-space checks may be cited as prior evidence.
- A nonzero singlet overlap is not a singlet projector, and a generic-twist Bethe eigenstate is not automatically a highest-weight CG state.
- Preserve the distinction between unnormalized transfer operators in the localization table and normalized Q-ratio formulas. Explicitly record all scalar factors.
- Use the true companion-frame SoV basis and the specified untwisted monodromy for B. Do not repeat rejected basis/ratio searches without a new reason.
- A product of two individually localizing operators need not localize. Repeated brackets need an independent derivation.
- Two numerical points support an identity but do not prove symbolic parameter dependence.
- Keep searches finite and declared in advance. Report failure within the tested ansatz; do not claim a general impossibility theorem from it.

## 2. Team and ownership

Use one parent coordinator and at most three active subagents. A–F below are task identities, not six simultaneous workers. Reuse slots as assignments finish.

The parent alone:

1. Owns live-kernel evaluation and notebook edits.
2. Integrates accepted results into `main-notes.tex` and updates the handoff.
3. Resolves convention conflicts and decides which dependent work can proceed.
4. Runs requested computations in short, visible blocks and returns their evidence to the requesting subagent.

Subagents own separate reports under `subagent-results/`: `A-row8.md`, `B-cg.md`, `C-contours.md`, `D-quadratic.md`, `E-mirror.md`, and `F-generalization.md`. They may prepare correspondingly named `.wl` source files, but do not evaluate them themselves or edit shared notebook/notes files. Create these reports only when their assignment starts.

Read-only notebook context requests may be coordinated by the parent to avoid changing the active target. No subagent may call `setPoint`, `computeSpectrum`, kernel restart/abort, or other kernel-mutating tools independently.

Each computation request must contain the intended formula/question, prerequisites, complete bounded code, expected output, and a stopping condition. The parent returns cell IDs, parameter point, output, messages, and execution status.

## 3. Parent preflight

1. Read the Wolfbook skill and applicable repository instructions. Check connection and notebook/kernel binding with `wolfbook_status` and notebook context tools.
2. Read the relevant sections of `L2-localization-tests.wb`; identify the cells defining the localization table, `singletFromQ`, `sovPipeline`, single-site generators, rotation conventions, and existing CG reference data. Record cell IDs, not just positions.
3. Determine current kernel state explicitly. Do not assume the handoff's point-2 state persists. Do not rerun a full spectrum computation unless a scheduled validation needs it.
4. Create a short convention sheet in the parent findings block: GT ordering, covector action, P and g0 directions, normalized versus raw t/R, D definition, Q gauge, and the highest-weight normalization used by existing CG data.
5. Send this convention sheet and relevant formula/cell references to A, C, and D.

Use visible notebook cells, inspect returned output, and explicitly save after each completed block. Use task-specific WL contexts or Module-local variables. Display findings with Grid. Preserve the handoff's multiline Dot precautions.

If the controller is unassociated, reopen the existing notebook with `wolfbook_newNotebook(path=..., overwrite=false, target=true)`, verify controller selection, and retry. Do not use selfTest for this error. Diagnose disconnection before considering any fallback; do not infer license failure from a sandboxed kernel.

## 4. Agent A — explicit row 8 and normalization

Dependencies: parent convention sheet. First-wave priority; its result unlocks B and E.

Target:

`<0| R(theta1) [t21(theta1) - t11(theta1) R(theta2)/(T+1)] = <v8| tensor e8^T`.

Steps:

1. Identify every scalar convention in this raw-operator identity and in the recorded normalized Q-functional. Produce a conversion table.
2. Derive the eight components of v8 using single-site action and auxiliary-space contractions. Use existing recorded localization as evidence, but seek a componentwise symbolic derivation without constructing a coupled operator.
3. Express v8 as a short single-site generator polynomial if possible; the explicit GT row is the required deliverable even if the polynomial is not short.
4. Transport the covector to the diagonal/CG frame, showing every P or g0 factor and inverse. Identify the coefficient responsible for the nonzero singlet overlap in the established conventions.
5. Fix C analytically from the actual SoV normalization. Inspect the local Ryan–Volin source around equation (5.22) and its normalization definitions; give page/equation references. Do not infer Gamma factors from the handoff alone.
6. If analytic C stalls, give a reproducible reference-state prescription for C and distinguish this from a closed formula. Validate its state independence using existing Q/eigenvector data, without introducing new coupled operators.

Acceptance: explicit GT row, frame conversion, normalization formula or precisely labeled remaining normalization gap, and reproduction at both existing points. Symbolic cancellation must list denominator assumptions. If any coefficients are fitted, use a third generic rational point as independent validation.

Stop: after the explicit row and one analytic normalization attempt, report the exact unresolved equation instead of trying unrelated normalization fits. Release usable results to B/E immediately.

Deliverable: `A-row8.md`, computation source if needed, and a proposed append-only notes subsection. Label what is proved, independently validated, conditional, or open.

## 5. Agent B — localized basis to CG data

Dependencies: A's explicit row/conventions; analytic preparation may start earlier after a slot opens. Final rank/extraction claims wait for A.

Steps:

1. Extract or derive the single-site GT rows v1,...,v8 corresponding to all eight localization-table entries. Retain raw and normalized versions separately.
2. Form the 8 × 8 single-site change-of-basis matrix and determine its generic rank and exceptional denominators. Invert it only where justified.
3. Apply the recorded antipode/highest-weight reduction using one-site ladder equations. Construct the map from the six highest-weight channel coefficients, with multiplicity, to the available localized functionals. State its dimensions and basis order explicitly.
4. Compute this map's rank. If it is deficient, give its nullspace and identify exactly which channel or multiplicity combination is missing. Do not equate single-site rank eight with rank six on highest-weight data.
5. Fix the two-octet multiplicity convention by matching existing project CG conventions. If none exists, describe the remaining basis freedom and give a clearly labeled working convention; do not present it as canonical.
6. Identify the bridge from generic-twist Q-data to the intended highest-weight CG states. Reuse a justified project construction if available. If a twist limit is required, specify its path, normalization, and treatment of degeneracies before attempting it; no automatic large limit scan.
7. Give formulas or a finite algorithm for extracting every recoverable coefficient. Check single-site highest-weight relations and normalization against existing recorded CG references.

Acceptance: explicit localized basis and rank evidence, multiplicity convention, and an extraction recipe whose twist/highest-weight assumptions are stated. A rank obstruction with an explicit missing direction is a valid bounded result, but does not complete the overall CG objective.

Deliverable: `B-cg.md`, basis/map tables, and proposed notes subsection.

## 6. Agent C — bounded cross-site contour construction

Dependencies: parent conventions; can begin in parallel with A. Full prefactor validation waits for A's normalization.

Target numerator:

`D1(0,2) Q1(theta2+h) + D1(0,1) Q1(theta2+2h)/(T-1)`.

Track the outer Q1 factor and reference-overlap denominator separately; do not silently absorb state-dependent terms into C.

Steps:

1. Read the existing `L2-det-mathcalA` derivation and its notebook context. Write out its measure, contour, pole locations/orders, and normalization explicitly.
2. Derive a residue table showing what this measure contributes at the required site-1 and site-2 poles. Determine whether the desired two monomials can appear while unwanted monomials cancel.
3. Write the exact finite linear/algebraic constraints on the state-independent measure coefficients. First seek a universal residue identity; if QQ relations are needed, name them and verify that they apply with these shifts.
4. Test the existing measure and, only if needed, one smallest extension justified by a specific missing residue. Declare that extension's poles, numerator degree, and unknown coefficients before solving.
5. Verify coefficient ratio `1/(T-1)` symbolically, including dimensions and large-T behavior. Failure of the earlier factorized ansatz does not rule out all double contours.
6. Only after the undeformed identity passes, introduce a k-deformation with a stated reference value. Derive which one-body covector its coefficients represent. A residue interpolation with no one-body interpretation is an intermediate result.

Acceptance: a contour/residue identity with stated pole prescriptions and parameter restrictions, reproducing the known Q-functional. Numerical checks use available all-state data at the two points and an independent point if fitting occurred.

Stop: if the existing measure and one justified minimal extension fail, return the inconsistent residue equations. Choosing a different FSoV framework or beginning an open-ended ansatz search requires a concrete proposal to Paul first; unaffected assignments continue.

Deliverable: `C-contours.md`, residue table, candidate ledger, and either a derivation or the exact bounded obstruction.

## 7. Agent D — quadratic single-site insertion

Dependencies: parent conventions. Can proceed alongside A and C; use A/B results when available.

Primary target: the correctly conjugated single-site insertion `t11^(1)(theta1)^2` in the existing site-1 family.

Steps:

1. Compute its action with 8 × 8 single-site operators, including g0 conjugation and the family covector. Record its GT components and the coefficient that permits singlet access.
2. Expand the resulting covector in the localized basis when B's basis table is available. Use the established Q-functional for each row to obtain a finite Q-expression; do not insist on a compact bracket before obtaining this expression.
3. Determine how much of the expression lies in row 8 and which other rows are needed. This tests whether the quadratic and localization constructions are two descriptions of the same missing direction.
4. Attempt one compact representation based on the known bracket/QQ relations. If proposing repeated Q12 brackets, derive their composition law explicitly; never substitute a square of eigenvalues for the square of a localized single-site insertion without proof.
5. If the primary insertion fails for a structural reason, examine only the already successful fused `t22^(1)(theta1-h)` alternative and explain the choice.

Acceptance: a Q-functional demonstrably equal to the chosen one-body insertion, with conjugation and scalar factors fixed. A finite sum over the localized basis is acceptable; a compact contour is an optional improvement.

Stop: one primary insertion, at most one justified fallback, and one compactification attempt. Report any dependence on missing row normalizations to A/B.

Deliverable: `D-quadratic.md`, explicit row/basis coefficients, and proposed notes subsection.

## 8. Agent E — site-2 mirror

Dependencies: A's conventions and normalized site-1 identity; B's basis table is useful but not required to begin.

Steps:

1. Derive the mirror transformation from the actual monodromy ordering, shifted tower, and vacuum. Determine whether transpose/duality or a frame change is needed; do not merely replace theta1 with theta2.
2. Derive the site-2 counterpart of row 8 first as `e8^T tensor <v8_mirror|` with an explicit one-site row and nonzero singlet access.
3. Derive its Q-functional using upper/lower-Q and Hodge conventions already recorded. Track every shift and normalization pole.
4. Extend to the remaining seven rows using the same transformation and give the full site-2 table.
5. Compare the mirror singlet functional to the site-1 result after all frame and normalization conversions. Use any discrepancy to identify a specific convention or identity failure.

Acceptance: complete mirror table, explicit rows, and two-point checks; classify symbolic versus numerical evidence separately. If only row 8 can be justified, return that result and the precise obstruction to mirroring the rest.

Deliverable: `E-mirror.md` and proposed notes subsection.

## 9. Agent F — symbolic proof, exceptional points, one further representation

Dependencies: stable A identity and B/C formulas where relevant. First do proof preparation; do not speculate about many representations while waiting.

Steps:

1. Audit the accepted formulas for denominators, singular twist choices, reference overlaps, and branch conventions for Q exponential factors. Separate removable normalization singularities from possible failures of basis rank or localization.
2. Prove remaining table identities with symbolic theta, h, z using one-site/auxiliary-space relations. Start with one representative correction row, then reuse its argument; A remains owner of the row-8 proof.
3. Validate any fitted coefficient at a third generic rational point, chosen away from all listed exceptional loci and recorded before evaluating it. Do not call this a symbolic proof.
4. Select exactly one further representation based on a written criterion. Prefer identical self-dual sl(3) representations if testing the singlet mechanism; the gl(3) one-dimensional target may carry a determinant weight. Check central weights rather than calling it a trivial gl(3) singlet without qualification.
5. A concrete candidate is GT highest weight [4,2,0] at both sites (dimension 27 per site): it retains an sl(3) singlet channel and permits lambda2=2. Confirm conventions first. Work only with single-site matrices and derived functional identities; do not build the 729-dimensional coupled operator problem.
6. Test one prediction: whether the quadratic mechanism becomes available in the original bracket expansion, and whether it has the required one-body singlet-reaching component. Do not promise a full new SoV basis or Bethe spectrum.

Acceptance: a table of genericity conditions with evidence, clearly scoped symbolic proofs, and one interpretable additional-representation test. If the new representation requires unproved formulas, report those prerequisites instead of starting a broad computation.

Deliverable: `F-generalization.md`, exceptional-locus table, and proposed notes subsection.

## 10. Dispatch order and integration gates

1. Parent preflight.
2. Wave 1: A (row 8), C (existing-measure residue analysis), D (single-site quadratic action). These have independent analytic work. Parent schedules kernel requests sequentially, with A first.
3. As A releases explicit formulas, launch B in the first free slot. D's basis expansion waits for B; C's full normalization waits for A. Agents should return intermediate reports and release their slots rather than sit idle on dependencies.
4. Start E once A's conventions are stable and a slot is free. Start F once there is a stable identity to prove or audit. Keep at most three subagents active.
5. Reassign bounded follow-ups to the same task agents when dependencies arrive. A reviews E's normalization; B checks D's basis expansion; C reviews F's pole classification. Reviews are reports, not concurrent edits of each other's files.
6. Parent integrates only claims with supporting formulas and execution evidence. Append accepted results to `main-notes.tex`; preserve earlier findings and flag corrections explicitly. Add notebook Findings cells after each computation block and save.
7. Parent updates the handoff with completed results, exact open questions, and the current kernel state. Compile the notes PDF if the LaTeX source changed and the project build is available.

## 11. Validation and reporting contract

Every agent report must contain:

1. Result and its status: structurally excluded, numerically rejected, pair-fitted, independently validated, conditional, or open; explicitly identify symbolic proofs.
2. Definitions, gauge/frame, normalization, and assumptions needed to use the formula.
3. Derivation and a minimal reproducible computation request or executed cell references.
4. Exact checks where possible; otherwise working precision, absolute and relative residuals, states/points tested, and treatment of vanishing reference overlaps.
5. What remains unresolved, the exact next test, and which other assignment depends on it.

For numerical functional checks, use 60-digit working precision when supported by the underlying data, aim for relative residual below 10^-40, and report actual residuals rather than a pass flag alone. Use absolute residuals and a stated scale when the target is zero. Raise precision only to diagnose a specific precision problem. Do not divide by a vanishing reference overlap; switch to a stated nonzero reference or compare unnormalized identities.

Default computation blocks should finish within about 60 seconds. Split longer jobs into bounded steps; use asynchronous Wolfbook execution and status reporting when needed. A timeout is not evidence of failure or non-execution. Do not interrupt unrelated shared-kernel work.

If a contradiction with the notes appears, preserve both formulas and the minimal counterexample. Pause the affected branch, report the conflict, and seek resolution; continue independent branches. If an open-ended search would be needed, present its precise scope and expected discriminator before asking Paul to choose.

The project is complete only when the one-body construction and its CG extraction are justified. Successful bounded reports that expose an unresolved contour or twist-limit obstruction must remain labeled open; completing all assignments is not itself proof that the mathematical objective has been achieved.
