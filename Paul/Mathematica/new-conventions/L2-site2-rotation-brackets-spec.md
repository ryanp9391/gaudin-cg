# Spec: Stage A of `Two-site-rotation-brackets-plan.md` — the site-2 rotation with e₁ᵀ on site 1, the vacuum-family identification, and the meaning of the adjoint quadric (2026-10-11)

Delegated to a wolfbook-builder subagent (sonnet); planned and reviewed by the main session. **New notebook:**
`Paul/Mathematica/new-conventions/L2-site2-rotation-brackets.wb` (self-contained, rep-generic with `Λ` as the only
rep switch; first code cell `Quiet[ClearAll["Global`*"]]`; run [2,1,0] first, then [3,1,0]; [3,2,0] only if everything
else is done and each rep took < 10 min). Create it with `wolfbook_newNotebook` so VS Code has it open and cells
execute natively; if cells of the new file will not execute, run the sources in K1 and store printed results as
markdown cells under each code cell (say so in the summary).

Reference material (read before coding):
- Plan: `Two-site-rotation-brackets-plan.md` §0–§2 (this spec implements §0's "first check" and §2 = Stage A, plus §7's
  side question).
- Write-ups: `L2-210-site1-rotation-brackets.md` (the adjoint site-1 construction), HANDOFF.md top section (the general
  site-1 results: label rule, counts 5/7/11 and 4/5/9, loci, closed-form measures, u-structure).
- Cell sources (scratchpad, UTF-8 text dumps; named in the task message): `ref-general-all.txt` = all 63 cells of
  `L2-site1-rotation-brackets-general.wb` — **copy its machinery**: cells 3–6 (setup: generators `JJ`, `lax1/lax2`,
  exact sparse transfer matrices `tL2`, eigenstates, dual Q's `QU`/`vm`, Gauss factors `Pphys`, companion kets), 9–10
  (B operator `Bop2`, SoV covectors, label rule `labRule`, monomials `siteMono`, wavefunction test `wfTest`), 13 (K map
  `Kfull`, groups `grp`, factors `facL`, `gaussCov`, `ovlExact/ovlK`), 16 (tensors `tensorOf`, Plücker tools `pfaff`,
  `pf4list`, `pf6list`, `decompForm`, brackets `brQ`, `Bk2`, `rankTab`, `supportTab`, `autoGen`), 19–21 (loci:
  `colPf`, `satGen`, `dimGen`, `conicOf`, `ccGuess`), 24–25 (reduced forms, `formData`, `fgList`, 225-state
  verification), 28 (reach table: `L0`, `rhoLow/rhoUp`, `undone`), 31–41 (closed-form measures, t-dependence), 43–53
  (u-structure: gauge-invariant kernels, F(u), G(u)). `ref-rotation-sov-all.txt` = `L2-rotation-sov.wb` (cell 25:
  `Ωc[y_] := uv[1] . PppW[y]` the companion-frame one-site vacuum; cell 98: `oneSiteVac` diagonal-frame Gauss
  parameters r_k(z,w); cells 81–85 Part B).

## 0. Conventions (identical to the general notebook)

Companion frame, L = 2, new conventions ℒ_ij(u) = uδ_ij + hE_ij, plain GT generators, h = 1/2, θ = (2/7, −3/5)
(t = θ₁₂/h = 62/35), ket twist z = (3/2, 5/4, 2/3), z′ = (3/2, 5/4, 3/7) for z-independence, (h,θ) = (1/3, (1/5,−4/9))
(t = 29/15) for t-dependence. `prec` = 60. SoV labels (a|bc) per site, λ₂ ≤ a ≤ λ₁, 0 ≤ c ≤ λ₂, c ≤ b ≤ a; roots
θ_s + {−1−a, −b, −c}h; monomial M^{(s)}_{(a|bc)} = det₂[v(θ_s−(1+a)h), v(θ_s−bh)]_{rows 2,3}·v₃(θ_s−ch), normalised by
the vacuum (λ₁|λ₁λ₂). Rotation g = exp(αE₁₂)exp(βE₁₃)exp(γE₂₃) (`gaussCov`), first row (1, α, β+αγ). Brackets on
ladder s: ⟨Q^jF⟩_s = Σ_m F_m Q^j(θ_s+mh), m = −1−λ₁..0; Bk₂^{(s)}[F,G] = ⟨Q³F⟩_s⟨Q²G⟩_s − ⟨Q²F⟩_s⟨Q³G⟩_s.

## 1. Tasks (each its own code cell(s) + a markdown result cell with a `Grid` of checks; save after each)

**Task 0 — setup + the vacuum-family identification (plan §0, do first, both reps).**
(a) Rebuild the setup (general cells 3–6, 9–10) for Λ; confirm the label rule and the wavefunction monomial
(`wfTest` max deviation) as regression.
(b) Companion-frame one-site vacuum Ωc[{w₁,w₂}] = e₁ᵀ exp(E₁₂/w₁) exp(−E₁₃/(w₁w₂)) exp(((w₁+w₂)/(w₁w₂))E₂₃) (rebuild
from the formula with the rep's generators; `PppW` of the rotation-sov notebook is this product). Check EXACTLY (rational
w): Ωc[{w₁,w₂}] == gaussCov[{1/w₁, −1/(w₁w₂), (w₁+w₂)/(w₁w₂)}], and β+αγ = α² on it. Then the converse: a rational
point (α,γ) of the [3,1,0] best locus β = α²−αγ equals Ωc[{1/α, 1/(γ−α)}] exactly. Confirm Ωc[w]⊗Ωc[w] is an exact
left eigen-covector of the companion-frame transfer matrices with bra twist w (i.e. of t_a(u) twisted by the
companion matrix G(w) — check at one u, residual 0 or 10⁻⁵⁰), so "Ωc" really is the vacuum.
(c) Rewrite the two [3,1,0] closed-form λ₂-group measures of general cells 47–48 with α = 1/w₁, s = γ/α = 1 + w₁/w₂:
F = h⁴α^{(u−θ₁)/h}/(P₃′(u)(u−θ₁)(u−θ₂+h)) etc. become w₁^{−(u−θ₁)/h} × (measure) × 1/(u−θ₂+h) × …; print both forms
and verify them against the tensors at a rational w (this is the "vacuum-Q structure" statement of plan §0 — report
whether the ω(θ₁−h) = 3(s−1)/(3s−2) exception has a w-reading, e.g. 3w₁/(w₂+3w₁)).

**Task 1 — the mirror: Φ = e₁ᵀ ⊗ e₁ᵀρ(g₂), support and K maps (plan §2 (a),(b)).**
For the dim-V GT covectors e₁ᵀ⊗e_νᵀ (and a generic rational g₂): expansion in the L=2 SoV basis (`Kfull`-type solve,
rationalise with error report), the support — which SITE-1 labels occur. Expectation to test: only (λ₁|λ₁ c₁),
c₁ = 0..λ₂, with factors λ₁^{(c₁)} := v₃(θ₁−c₁h)/v₃(θ₁−λ₂h) (k-independence as before). If the site-1 support is
different (the SoV construction is site-asymmetric: B is built from the ordered monodromy, and ⟨x| are left
eigenvectors), record exactly what it is: list of site-1 labels per site-2 label, and whether the factors are still
single Q-ratios. Rank of the K map per site-1 group; z-independence (same exact rationals at z′).

**Task 2 — site-2 tensors and automatic counts (plan §2 (b)).**
Per site-1 group: the column tensor C^{(c₁)}_{(nn′),c′} on the ladder-2 points (minor index (n,n′) over m = −1−λ₁..0,
column c′ = lone point Q³(θ₂−c′h)), populated slots and support per column, 2-form ranks for a generic rotation and a
random covector, automatic count per group and total. Compare with the site-1 numbers (5/7 generic, 4/5 reduced). Verify
the bracket identity on ladder 2 and the full overlap identity ⟨Φ|Ψ⟩/⟨x₀|Ψ⟩ · M^{(2)}_vac = Σ_{c₁} λ₁^{(c₁)} Σ_{c′}
Σ_wedges Bk₂^{(2)}[F,G] Q³(θ₂−c′h) on all states (10⁻⁴⁰ or better).

**Task 3 — loci in g₂ (plan §2 (c)).**
Gauss chart, exact: per group and column the 4×4-Pfaffian ideal, `Factor`, dimension, genuine content, conic constants
cc via `conicOf`/`ccGuess`; the joint loci; the best (fewest-brackets) locus. **Key question:** is the site-2 best locus
again the vacuum family β+αγ = α² (for [3,1,0]) / the quadric 2(t−2)(β+αγ) = (t−3)α² (adjoint), or does the asymmetry
change the constant (try cc(t) with t → −t, t → t−1, t → t+1, t → 2−t, …, and the general ansatz cc = (pt+q)/(rt+s) from
the two (h,θ) points plus one more if needed, (h,θ) = (1/2, (3/7,−1/4)))? Also test directly: do the covectors
Ωc[{w₁,w₂}] (rational w) lie on the site-2 best locus? Reach table at finite twist for the best locus (general cell 28
mirrored: here site 1 = e₁ᵀρ(U₀) is the vacuum and the undone site-2 covector is e₁ᵀρ(g₂L₀D₀) — check the undo
identity on the HW kets exactly as in the site-1 case, note that ρ(U)⊗ρ(U) fixes HW kets so the same argument applies).

**Task 4 — closed-form measures and u-form on the best site-2 locus (plan §2 (d)).**
Mirror general cells 31–36 and 43–48: F, G entries symbolic in (α,γ) (and t via two points + the 10-point fit only if
the form is not obvious); gauge-invariant kernels; explicit F(u), G(u) for the columns that admit them — expectation:
w₁′^{−(u−θ₂)/h} × ladder-2 measure × cross-ladder pole (u−θ₁ …) × low-degree polynomial. Verify to 10⁻⁴⁰ on all
states. Output: the site-2 analogue of the 5-bracket formula, written out.

**Task 5 — the adjoint quadric is a torus image of the vacuum family (plan §7 side question; cheap, exact).**
Analytic fact to verify (main session): for D = diag(d₁,d₂,d₃), e₁ᵀρ(g)ρ(D) ∝ e₁ᵀρ(D⁻¹gD) and D⁻¹gD has
(α,β,γ) → (α d₂/d₁, β d₃/d₁, γ d₃/d₂), so β+αγ scales by d₃/d₁ and α² by (d₂/d₁)². Hence the torus with
d₁d₃/d₂² = cc maps the vacuum family {β+αγ = α²} ONTO the quadric {β+αγ = cc·α²}: **the adjoint quadric is
{Ωc[w]·ρ(D) : D torus, d₁d₃/d₂² = (t−3)/(2(t−2))}**, and the [3,1,0] locus (cc = 1) is the vacuum family modulo the
"unimodular-type" torus d₁d₃ = d₂². Verify exactly on the rep (ρ(D) = exp(Σ log d_i E_ii) on V; Ωc[w]·ρ(D) is on the
quadric for rational w, d with d₁d₃/d₂² = 43/16 at P1). Then: (i) is such a D the resolvent/companion-frame image of a
known object — compare with the §5 product point of `Two-rotations-det-overlap.md` (Ω_{(w₁,y*)}ρ(D) with D =
diag(1/(y−z_k))-type in the diagonal frame: compute d₁d₃/d₂² for that D and see whether it equals cc for some y;
report the relation, don't force it); (ii) the relation cc = (3/4)κ_singlet, κ_singlet = 2(t−3)/(3(t−2)) is then
"d₁d₃/d₂² = (3/4)κ" — note it, no more.

**Task 6 — markdown summary** `### site-2 rotation brackets (Stage A) — result (2026-10-11)`: per rep the vacuum-family
identification (exact), the site-1 support rule for site-2 rotations (as found), K-map ranks, tensor supports,
automatic/reduced counts side by side with the site-1 numbers, the loci and whether the best one is the vacuum family
(with cc(t)), the reach table, the closed-form measures and u-forms, the torus statement for the quadric, timings,
surprises, and what this means for Stage B (two-site K tensor): in particular, whether the two one-site supports are
"vacuum on the other site + one Q-ratio" on both sides, which is what makes the (c,c′) block decomposition of plan §1
well-defined. Save. Do not write HANDOFF (the main session does).

## 2. Rules and gotchas

Exact wherever possible; eigenstates/Q's at `prec`. Time one state/one Gröbner column before loops; nothing over ~3 min
without a timed smaller run; `TimeConstrained` 120 s on symbolic solves, fall back to `Reduce` on a chart or rank
tables at rational points (say so). Markdown result cells under code cells; save after each task. Do not edit other
notebooks, `paper/*.tex`, HANDOFF.md, git. If the site-1 support rule (Task 1) is not of the one-Q-ratio form and you
cannot characterise it within two attempts, record the raw support tables and continue with the groups as found.
- Shared kernel K1: the first cell wipes other notebooks' state. Never write `*)` inside a comment. No `D`, `E`, `N`,
  `I`, `C`, `K` as symbols. Memoized definitions: purge with literal symbols. Semicolons between statements.
- `wedge2`-type symbol clashes with Table variables: use fresh names. Select labels by exact equality, not `Position`.
- `Max[{}]` = −∞; `SetPrecision` before `Norm` of near-zero vectors; exact zero tests with `===` on rationals only.
- `Inverse` of the 225×225 SoV covector matrix at 60 digits is fine (the general notebook did it); report the
  rationalisation error of the K maps (denominators ≤ 6 digits expected).
