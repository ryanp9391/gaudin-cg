# Spec: Stage C + Q-side payoff for the adjoint (∞,∞) two-vacuum family — closed-form measures, u-form, and all six CG data from the bracket Q-number alone (2026-10-11)

Delegated to a wolfbook-builder subagent (sonnet); planned and reviewed by the main session. **New notebook:**
`Paul/Mathematica/new-conventions/L2-210-two-vacua-brackets.wb` (adjoint [2,1,0] only; self-contained; first code cell
`Quiet[ClearAll["Global`*"]]`; build/run as the previous stages: JSON, `code` CLI, sources evaluated in K1, printed results
as markdown cells; say so).

Read first: HANDOFF.md top section (Stages A, B, D and the follow-up), `L2-two-site-rotation-brackets.wb` cells 47–49
(Task 4b: the (∞,∞) family, its blocks, the explicit 12-wedge decomposition at (a,a′) = (2/3,−5/7)) and 55–58,
`L2-two-vacua-untwist.wb` cells 17, 29, 33, 37–40 (closed-form ε = 0 Gauss parameters, K/S, CG extraction, the four
overlaps at small ε), the general notebook's closed-form/u-structure cells (`ref-general-all.txt` cells 31–53: pivot rule,
gauge F_i = 0, F_j = −1, G_j = 0, t-fit over (h,θ) points, gauge-invariant kernels, class tests, explicit F(u), G(u)).
Dumps in the scratchpad: `ref-two-site-all.txt` (Stage B notebook, 55 cells before the follow-up — re-read cells 47–49
from the notebook itself), `ref-two-vacua-untwist-all.txt` (Stage D), `ref-general-all.txt`, `ref-210-untwist-CG-all.txt`.

## 0. The family and the known facts

Companion frame: Φ(a,a′) = e₁ᵀexp(aE₁₂)exp(aE₂₃) ⊗ e₁ᵀexp(a′E₁₂)exp(a′E₂₃) = gaussCov[{a,0,a}]⊗gaussCov[{a′,0,a′}] =
Ωc[(1/a,∞)]⊗Ωc[(1/a′,∞)] (first rows (1,a,a²), (1,a′,a′²)). Diagonal frame: Ω_{(w₁,∞)}⊗Ω_{(w₁′,∞)}, w₁ = 1/a, w₁′ = 1/a′
(Ω_w = e₁ᵀρ(n(z,w)), `oneSiteVac` with w₂ → ∞: r₁ = (w₁−z₁)/(w₁−z₂), r₂ → (w₁−z₁)(z₂−z₃)/((z₁−z₂)(w₁−z₃)) × (−1)·…
— compute the limit, don't guess). Facts: CG-complete at finite twist for a ≠ a′ (Stage D cell 21); untwists cleanly
in the regime z = 1+εφ, w₁ = 1+εω₁, w₁′ = 1+εω₁′ (K/S = O(ε) for all six branches, Stage D cell 29); bracket count 12
(blocks (c,c′) = (0,0),(0,1),(1,0),(1,1) with ranks 3,3,2,2 and P = 4,4,2,2 → wait: Stage B reports P = 4,4,2,2 with
ΣP = 12 and the wedge list "3 terms in (0,0), 3 in (0,1), 2 in (1,0), 2 in (1,1)" — resolve this (the 3+3+2+2 = 10 wedge
count vs ΣP = 12) in Task 1 and state the correct count; lower bound Σr = 10). Rank drops on a = a′ (global rotation,
count 8), a = 0, a′ = 0. The identity to realise in closed form:

  ⟨Φ|Ψ⟩/⟨x₀|Ψ⟩ · M¹_vac M²_vac = Σ_{c,c′∈{0,1}} Q³(θ₁−ch) Q³(θ₂−c′h) Σ_k Bk₂^{(1)}[F_k^{(cc′)}, G_k^{(cc′)}] Bk₂^{(2)}[F′_k^{(cc′)}, G′_k^{(cc′)}],

Bk₂^{(s)}[F,G] = ⟨Q³F⟩_s⟨Q²G⟩_s − ⟨Q²F⟩_s⟨Q³G⟩_s, ⟨Q^jF⟩_s = Σ_{m=−3}^{0} F_m Q^j(θ_s+mh), M^{(s)}_vac = minor_s(−3,−2)·Q³(θ_s−h).
Conventions: h = 1/2, θ = (2/7,−3/5) (t = 62/35), z = (3/2,5/4,2/3), prec 60 for the finite-twist part; 200 digits for the
small-ε part (Stage D's machinery).

## 1. Tasks (each its own cell(s) + markdown result cell with a `Grid`; save after each)

**Task 0 — setup and regression.** Stage B setup (cells 2–12 of its notebook: transfer matrices, states, dual Q's, SoV
basis, monomials, K tensor) + its tools cell; rebuild the blocks of Φ(a,a′) symbolically in (a,a′) (polynomial entries),
confirm ranks (3,3,2,2) and the count at (2/3,−5/7) (reproduce Stage B's numbers exactly), and the overlap identity on all
64 states with explicit bracket products (10⁻⁴⁴).

**Task 1 — closed-form wedge decomposition in (a,a′).** Per block, a rank decomposition A^{(cc′)} = Σ_k ω_k⊗ω′_k over the
function field Q(a,a′) with each ω_k, ω′_k a 2-form on the four ladder points; then each 2-form into wedges F∧G with the
pivot rule and gauge of the general notebook (F normalised at a pivot point, G zero there); report the minimal number of
products you can achieve exactly (try the orientations and a few pivot choices; the lower bound is Σr = 10 — say whether
10, 12 or in between), and print every F, G, F′, G′ as rational functions of (a,a′) (`Factor`). Verify the symbolic identity
(blocks reconstructed exactly) and the numeric identity on all 64 states at three (a,a′) points (10⁻⁴⁴). Structure to
report: which entries are monomials in a (resp. a′) — expectation from the one-site results: a-power law F_m ∝ a^{m+…} on
ladder 1 and a′-power law on ladder 2 (the bra twists w₁^{−(u−θ₁)/h}, w₁′^{−(u−θ₂)/h}); which factors couple a and a′
(the cross-ladder coupling of the two vacua — this is the new information); where a−a′ appears.

**Task 2 — t-dependence and the u-form.** The blocks' coefficients depend on (h,θ) only through t = θ₁₂/h (check: two
(h,θ) with the same t; the K tensor is z-independent and h-homogeneous — verify the homogeneity degree). Fit the t-dependence
of the gauge-fixed F, G entries over ≥ 10 rational (h,θ) points (re-run setup + K per point in a driver cell as the general
notebook's Task 7 driver; ~1 s per point for the adjoint) with the rational-function fitter (`FindRationalFunction`-type or
linear solve for P(t)/Q(t) of bounded degree), validate on 2 held-out points, and report pole positions in t. Then the
u-form: for each wedge, test the classes of the general notebook (SoV measure wgt(m) × twist exponential ρ^m × low-degree
polynomial in m, with cross-ladder poles 1/(θ_s+mh−θ_{s′}+kh) allowed) for F, G, F′, G′; give explicit F(u), G(u) for every
wedge that admits one and say which do not (and what the obstruction is: rank, irrational roots, …). Expectation: ladder-1
kernels ∝ a^{(u−θ₁)/h}, ladder-2 kernels ∝ a′^{(u−θ₂)/h}, with the coupling of the two vacua appearing as a-, a′-dependent
polynomial coefficients. Also the Möbius/natural-measure tests only if cheap.

**Task 3 — the bracket Q-number at small ε and all six CG data from it (the payoff; Stage D's machinery).** Regime z = 1+εφ,
φ = (1/3,−2/7,5/11), a = 1/w₁ = 1/(1+εω₁), a′ = 1/(1+εω₁′); grid ω₁ ∈ {3/4, −5/3}, ω₁′ ∈ {−12/5, 7/2}; ε ∈ {1,2,4}·10⁻⁶;
HW-limit branches (6) at 200 digits, dual Q's of the branches (Baxter residual check), ⟨x₀|Ψ⟩ = 1 normalisation as in
Stage D (diagonal frame) — and the companion↔diagonal dictionary: the bracket formula is in the companion frame with
x₀ = e₁⊗e₁; relate ⟨Φ_c|Ψ_c⟩/⟨x₀_c|Ψ_c⟩ to ⟨Ω_{(w₁,∞)}⊗Ω_{(w₁′,∞)}|Ψ_d⟩/⟨x₀_d|Ψ_d⟩ (the κ scalars of Stage D cell 11; check
the relation numerically at one point). Then: (i) the Q-side number Qnum(B) := Σ_{cc′} Q³Q³ Σ_k Bk₂⁽¹⁾Bk₂⁽²⁾ / (M¹_vac M²_vac)
with the closed-form F, G at (a(ε), a′(ε)) (exact rationals in ε) and the branch Q's — verify Qnum(B) = ⟨Φ_c|Ψ_c⟩/⟨x₀_c|Ψ_c⟩
on all 21 HW-sector states (relative deviation; expect ~10⁻¹⁰⁰, limited by the Q precision); (ii) the ε-exponents of the F, G
entries, of each bracket ⟨Q^jF⟩, of each Bk₂, of each product and of each block sum, per branch — where the cancellations
sit (the Stage D Task 5 table has the exponents of the four overlaps: the mixed one equals the CG-slice exponent); (iii) CG
at leading order from Qnum alone: C^Q = Qnum/(prefactor · [φ_u(0)]_slice) with the ε = 0 undone covector of Stage D cell 17
(a′₀, c′₀, g′₀ at ω₂, ω₂′ → ∞ — take the limit of the closed forms and check it against the direct Gauss decomposition of
n(w)n(w′)⁻¹ at w₂ = w₂′ = ∞), for the four 1-dim slices and the octet 2×2 system from two grid points; grid spread and
error vs the chain slice ∝ ε; Richardson limits Q vs chain (expect ≤ 10⁻¹⁶); compare with the reference numbers
46.67453891, 10.92853892, −5118.538314, octet functional 184.3305612 / −238.5587523 and mixing −1.694 / 1.096 after
matching the normalisation. **Cost statement:** number of bracket products (12 or 10), distinct brackets ⟨Q^jF⟩ actually
needed (count distinct F's up to scale), Q-values needed (Q², Q³ on the four points of each ladder), vs the det₆ (36
bracket entries, one 6×6) and vs dim V² = 64.

**Task 4 — [3,1,0] note (no computation beyond what exists, or one cheap check).** From the Stage B follow-up: the (∞,∞)
family costs 45 (ranks 6,6,6,6; P 12,12,9,12; lower bound 24) vs det₆ 40 and generic 54. One cheap check: is the lower
bound 24 attainable, i.e. do the 2-form factors of a rank decomposition have rank 4 (needing 2 wedges each) in every
orientation, or is there an orientation/pivot with more rank-2 factors? Try the greedy choice (strip a wedge⊗wedge term
whenever some ω_k is decomposable) at one rational point; report the best ΣP found. Nothing else.

**Task 5 — markdown summary** `### adjoint (∞,∞) two-vacuum family: bracket formula and Q-side CG (2026-10-11)`: the
closed-form F, G (or a pointer to the cell with the full list), the count and cost statement, the t-laws and u-forms found,
the Q-number verification and ε-exponent table, the six CG data from the Q side with Richardson limits, the [3,1,0] note,
timings, surprises. Save. Do not write HANDOFF.

## 2. Rules and gotchas

As before. Shared kernel K1 (holds Stage D's state; your first cell wipes it — Stage D's notebook is saved). Exact where
possible; time one point before loops; nothing over ~3 min without a timed smaller run; markdown result cells under code
cells; save after each task; no other notebooks, no `paper/*.tex`, no git, no HANDOFF. Never `*)` in comments; no
`D`,`E`,`N`,`I`,`C`,`K` symbols; purge memoized definitions with literal symbols; `x0Of`-type factorised inverses at small ε;
`SetPrecision` before `Norm` of near-zero vectors; exact zero tests with `===` on rationals only.
