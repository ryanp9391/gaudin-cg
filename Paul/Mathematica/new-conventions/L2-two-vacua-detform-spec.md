# Spec (Job A): is ⟨Ω_w ⊗ Ω_{w′}|Ψ⟩ a single bracket determinant? — determinant-form search with free per-ladder bra data (2026-10-11)

Delegated to a wolfbook-builder subagent (sonnet); planned and reviewed by the main session. **New notebook:**
`Paul/Mathematica/new-conventions/L2-two-vacua-detform.wb` (rep-generic with `Λ`; [2,1,0] first, then [3,1,0]; build/run as
the previous stages: JSON, `code` CLI, sources in K1, printed results as markdown cells; say so). Paul's question: the
two-vacuum family untwists cleanly but costs 54 products of 2×2 bracket determinants for [3,1,0]; the det₆ (which costs
40 in those units) is ONE 6×6 determinant with 36 bracket entries. Is the two-vacuum overlap, or the CG-complete (∞,∞)
sub-family, also a single determinant of brackets with suitable (non-vacuum) bra data? `Two-rotations-det-overlap.md` §4
(Part B) tested only vacuum-Q data per ladder and failed; this job searches over FREE bra data.

Read first: HANDOFF top section; `Two-rotations-det-overlap.md` §1 (det₆ definition), §4; the Stage B notebook
`L2-two-site-rotation-brackets.wb` cells 2–12 (setup incl. det₆, w-vacuum Q's, Ψ_Q) and 47–49; Stage C summary (last
cell of `L2-210-two-vacua-brackets.wb`). Dumps: `ref-two-site-all.txt`, `ref-rotation-sov-all.txt` (cells 8, 25, 82),
`ref-310-two-rotations-all.txt` (cell 17).

## 0. Objects

Companion frame, h = 1/2, θ = (2/7,−3/5), z = (3/2,5/4,2/3), prec 60 (the per-state fits need ~40 digits). Dual Q's
v(x) = (Q¹,Q²,Q³)(x) of each state. Poles x_{s,n} = θ_s − nh, n = λ₃..λ₁ (3 per ladder for the adjoint, 4 for [3,1,0]);
8-pole weights wgt(p) = 1/Π_{q≠p}(x_p − x_q) (also try per-ladder weights wgt_s — a constant factor per ladder is irrelevant
only if the same on all rows of that ladder; keep both options). The per-ladder determinant det6L: rows (s,(i,j)),
s ∈ {1,2}, (i,j) ∈ {(1,3),(1,2),(2,3)}, entries the six bracket types summed over the poles of ladder s only:
  Σ_{p∈ladder s} wgt(p) {−Q^j(x)F_i(x), −xQ^j(x)F_i(x), Q^j(x)F_i(x+h), xQ^j(x)F_i(x+h), Q^j(x−h)F_i(x), xQ^j(x−h)F_i(x)}.
With F_i = w-vacuum direct Q's on both ladders, det6L ∝ det₆ ∝ ⟨Ωc[w]⊗Ωc[w]|Ψ⟩/⟨x₀|Ψ⟩ · Ψ_Q (verify; the Stage B
setup has det₆ — build det6L from it and check det₆ = const × det6L).
**Free bra data:** F_i^{(s)} are now unknown VALUES at the points where they are evaluated: x_p and x_p + h for the poles
of ladder s, i.e. m = −λ₁..1 (shifted points m = −λ₁−1..0 … list exactly which points occur), per function i ∈ {1,2} and
ladder s. Unknown count: adjoint 2 × 2 × 4 = 16, [3,1,0] 2 × 2 × 5 = 20, plus one overall constant.
**Targets** T(B) := ⟨Φ|Ψ_B⟩/⟨x₀|Ψ_B⟩ for Φ ∈ {Ωc[w]⊗Ωc[w] (control, must solve with F = vacuum Q's up to the gauge
freedom), Ωc[w]⊗Ωc[w′] (w = (5/2,7/3,3/5), w′ = (7/4,−2/5,9/7)), the (∞,∞) covector gaussCov[{a,0,a}]⊗gaussCov[{a′,0,a′}]
at (a,a′) = (5/7,−4/9) (NOT a degenerate vacuum point), gaussCov[g₁]⊗e₁ᵀ generic (control 2: the one-site 5/7-bracket
object — is it a det6L? its Part-A answer was no for restricted forms), and a random covector (control 3: must fail)}.

## 1. Tasks

**Task 0 — setup, det6L, controls.** Rebuild Stage B's setup; define det6L[B, Fdata] with Fdata an explicit table of
values (so that it is a polynomial in the unknowns); check det6L with vacuum data reproduces det₆ (const) and the vacuum
overlap (10⁻⁴⁵); check the gauge freedom: which transformations of F leave det6L/const invariant (row scalings F_i → c_i F_i,
the F₁ ↔ F₂ mixing F₁ → F₁ + κF₂? (rows (1,3),(1,2) use F₁, (2,3) uses F₂: F₂ → F₂ + κF₁ changes row (2,3) by κ×(row
(1,3)) → invariant), overall scale) — count the effective unknowns.

**Task 1 — the polynomial system (exact).** For a target Φ with exact rational SoV coefficients φ_m (from the K tensor),
the identity det6L[F] = c · Σ_m φ_m M¹_{x₁}M²_{x₂} must hold on every state. Approach (a), exact: treat the Q-values
Q^j(θ_s+mh), j = 2,3, at all points that occur as independent symbols; expand det6L (a 6×6 with symbolic entries linear in
the Q-symbols: each entry has ≤ (points)×2 terms; Expand is fine for the adjoint, maybe heavy for [3,1,0] — time it, use
`Det` with `Method -> "CofactorExpansion"` or Laplace along the ladder-1 rows: Σ_{20} det₃ · det₃ and expand each det₃
separately) and match coefficients of the Q-monomials against c·Σ φ_m M_m. This gives polynomial equations in F (degree
≤ 6, bilinear-ish per ladder: degree 3 in ladder-1 unknowns, 3 in ladder-2). Monomials of det6L outside the span of the
M_m must vanish — these give constraints; NB the true identity may hold only modulo the Q-relations (Baxter/QQ) rather
than as a polynomial identity: FIRST check whether the vacuum control holds as a polynomial identity with F = vacuum data
(if not, approach (a) is too strong — say so and rely on (b), but still report which monomials fail). If (a) is
meaningful: solve with `Solve`/`NSolve`/Gröbner over the gauge-fixed unknowns (fix the gauge first), `TimeConstrained`
300 s; report the solution set (dimension, rational points).
Approach (b), numeric on the true states: minimise R(F,c) = Σ_B |det6L[F](B)/Ψ_Q(B) − c·T(B)|² / Σ_B |c·T(B)|² over the
gauge-fixed unknowns and c (eliminate c by least squares), with `FindMinimum` at 40-digit working precision from ≥ 200
random rational starts (and from the vacuum data as a start for the Ωc[w]⊗Ωc[w′] target); report the best residual, the
number of restarts reaching it, and for the control it must reach 10⁻³⁰ or better from a fair fraction of starts (if the
control fails from random starts, the landscape is bad: switch to Newton/`FindRoot` on a square subsystem of states +
verification on the rest, or to approach (a)). A genuine solution has residual at the precision floor; a non-solution
plateaus (report the plateau value). Use all 64 / 225 states; for [3,1,0] time one evaluation of det6L on all states first.

**Task 2 — broader ansätze (only after Task 1, same machinery).** (i) Free j-pattern: rows (s, r) with (i_r, j_r) ∈
{1,2}×{2,3} arbitrary (enumerate the 3-row patterns up to symmetry: ((1,3),(1,2),(2,3)) is det₆'s; try all patterns
with three distinct (i,j) pairs and also repeated i with different j); (ii) six independent bra functions, one per row,
each with its own j_r ∈ {2,3} (unknowns 6 × points); (iii) the columns: in addition to the det₆ six, allow the
"x-polynomial" insertion degree up to 2 (columns Q^jF, xQ^jF, x²Q^jF, Q^jF⁺, xQ^jF⁺, Q^{j−}F, …) with a 6-column
SUBSET chosen from a larger menu — enumerate subsets only if cheap, otherwise take the det₆ columns and the two obvious
alternatives (drop the x-insertions in favour of Q^jF^{++}, i.e. F at x+2h). Report for each ansatz the best residual on
each target. Stop expanding ansätze once one solves the (∞,∞) or the generic two-vacuum target; then verify the solution
on all states at prec 60, try to rationalise the F values (denominators ≤ 8 digits) and, if rational, check at z′ and at
a second (w,w′) / (a,a′) point whether the F's are "vacuum-like" (compare with one-site vacuum Q's of various twists:
fit F_i^{(s)}(m) ∝ y^{−m} q(m) with q polynomial of degree ≤ 2 — the natural guess is the one-site vacuum Q's of twist
(w₁, w₂) on ladder 1 and (w₁′, w₂′) on ladder 2 but with cross-ladder shifted points; report).

**Task 3 — markdown summary** `### determinant-form search for the two-vacuum overlap — result (2026-10-11)`: per rep
and target: whether a det6L-type determinant exists (residual floor vs plateau, restarts), the ansätze tried, the gauge
count, the exact-polynomial diagnostic, any solution found with its bra data and verification, timings. State plainly:
"the two-vacuum overlap IS / IS NOT a 6×6 bracket determinant of the det₆ type with free bra data (within the ansätze
listed)". Save; no HANDOFF.

## 2. Rules and gotchas

As before (K1 shared; first cell wipes; no `*)` in comments; no `D`,`E`,`N`,`I`,`C`,`K` symbols; time one evaluation
before loops; ≤ 3 min per cell without a timed smaller run; `TimeConstrained` on symbolic solves; save after each task;
no other notebooks, no `paper/*.tex`, no git, no HANDOFF). `FindMinimum` with many unknowns: supply the gradient
analytically if slow (det6L is polynomial — `D` is fine but do not name a symbol `D`); use `WorkingPrecision -> 40`,
`PrecisionGoal -> 30`, `MaxIterations -> 2000`. For the exact approach, keep the symbolic Q's as `qq[s,j,m]`.
