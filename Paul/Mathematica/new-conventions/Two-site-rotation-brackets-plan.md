# Plan: two different rotations on the two sites, as ladder brackets — and a different site-2 reference (2026-10-10, to implement later)

Status: plan only, nothing implemented. Builds on `Site1-rotation-brackets-notes.pdf`, `L2-210-site1-rotation-brackets.md`,
`L2-site1-rotation-brackets-general.wb` (63 cells: label rule, K maps, tensors, loci, closed-form measures, u-structure)
and `Two-rotations-det-overlap.md` §1–§4 (companion frame, SoV basis, Part B: Ω_w⊗Ω_{w′}).

## 0. A fact to use from the start: the [3,1,0] locus IS the one-site vacuum family

Companion-frame one-site vacuum (`Two-rotations-det-overlap.md` §1, `Ωc[y]`): Ω_{(y₁,y₂)} = e₁ᵀ exp(E₁₂/y₁) exp(−E₁₃/(y₁y₂))
exp(((y₁+y₂)/(y₁y₂)) E₂₃), the same Gauss ordering as our g = exp(αE₁₂)exp(βE₁₃)exp(γE₂₃). Hence

    α = 1/y₁,   β = −1/(y₁y₂),   γ = (y₁+y₂)/(y₁y₂)   ⇒   β + αγ = (−y₁ + y₁ + y₂)/(y₁²y₂) = 1/y₁² = α².

So the [3,1,0] best locus β + αγ = α² (first row (1, α, α²)) is exactly {Ω_{(w₁,w₂)} : w₁ = 1/α, w₂ = 1/(γ−α)} — the
site-1 covector of the "5-bracket" family is the one-site w-vacuum, with α = 1/w₁ the bra twist. This explains the
exponent law α^{(u+u′)/h} = w₁^{−(u+u′)/h} (vacuum-Q structure w^{−x/h}) and why the λ₂-group measures look like
vacuum data. For the adjoint the best locus has cc = (t−3)/(2(t−2)) ≠ 1 and is NOT the vacuum family (the vacuum
family there costs the generic 5 brackets, the quadric 4). **First check in any implementation:** at a rational
(w₁,w₂), `Ωc[{w1,w2}]` (companion frame) equals e₁ᵀρ(g) with (α,β,γ) as above, exactly; and the (0,1)/(0,0)
closed-form measures F(u), G(u) rewritten with α = 1/w₁ (they should read w₁^{−(u−θ₁)/h} × measure × 1/(u−θ₂+h) …).

Consequence for this plan: "site-1 rotation ⊗ e₁ᵀ on the locus" = ⟨Ω_w ⊗ x₀^{(1)}|. The two obvious moves —
(i) a different site-2 reference, (ii) two different rotations — both land on the family

    ⟨Ω_w ⊗ Ω_{w′}|  (two different one-site vacua, 4 parameters),   and its generalisation  ⟨e₁ᵀρ(g₁) ⊗ e₁ᵀρ(g₂)|  (6 parameters).

Part B (`Two-rotations-det-overlap.md` §4) showed Ω_w⊗Ω_{w′} is not a single det₆ and its SoV coefficient matrix has
full rank 8 (adjoint). Neither matters now: the question is the number of bracket-determinant products.
And the untwisting of Ω_w⊗Ω_{w′} with w, w′ → 1 independently is clean by the exact commutator formula
(`L2-210-localising-span-spec.md` §6): the undone site-1 covector e₁ᵀρ(n(w)n(w′)⁻¹) has all three Gauss parameters
O(1) and generic (four direction parameters ω, ω′ for three Gauss parameters) — the degeneration of the old
det₆×det₂ family (y pinned) does not occur. So the goal is: **a bracket formula for ⟨Ω_w⊗Ω_{w′}|Ψ⟩ (or a cheaper
nearby family) with a count polynomial in the boxes, then its untwist.**

## 1. Objects and notation

- Companion frame, L = 2, rep Λ = (λ₁,λ₂,0); SoV labels (a|bc) per site, roots θ_s + {−1−a, −b, −c}h; wavefunction
  monomials M^{(s)}_{(a|bc)} = det₂[v(θ_s−(1+a)h), v(θ_s−bh)]_{rows 2,3} · Q³(θ_s−ch), normalised by the vacuum
  (verified for both sites in the general notebook, Task 1).
- Per-ladder brackets ⟨Q^jF⟩_s = Σ_m F_m Q^j(θ_s+mh), m = −1−λ₁..0; bracket determinants Bk₂^{(s)}[F,G]; a column
  2-form A on the ladder points; wedges = terms of its rank decomposition.
- **Two-site K tensor.** For GT-basis covectors e_{ν₁}ᵀ⊗e_{ν₂}ᵀ (dim V² of them), the SoV expansion
  e_{ν₁}ᵀ⊗e_{ν₂}ᵀ = Σ_{x₁,x₂} K_{ν₁ν₂;x₁x₂} ⟨x₁x₂|/N_{x₁x₂} in normalised-wavefunction units. Rational, (h,θ)-dependent,
  z-independent (as the one-site K maps were — verify). The one-site K maps are the slices ν₂ = 1 (site 2 = e₁ᵀ).
- For Φ = c₁⊗c₂ (c_s the GT components of e₁ᵀρ(g_s)): φ_{x₁x₂} = Σ c₁_{ν₁}c₂_{ν₂}K_{ν₁ν₂;x₁x₂}, bilinear in (c₁,c₂),
  and ⟨Φ|Ψ⟩/⟨x₀|Ψ⟩ · M^{(1)}_vac M^{(2)}_vac = Σ_{x₁,x₂} φ_{x₁x₂} M^{(1)}_{x₁} M^{(2)}_{x₂}.
- **Blocks.** Group by the lone points (c, c′) (site-1 lone Q³(θ₁−ch), site-2 lone Q³(θ₂−c′h)), c, c′ ∈ {0..λ₂}:
  A^{(c,c′)} = matrix (site-1 minors (m,m′)) × (site-2 minors (n,n′)), i.e. an element of Λ²(C^{λ₁+2}) ⊗ Λ²(C^{λ₁+2}).
  Then ⟨Φ|Ψ⟩·(…) = Σ_{c,c′} Q³(θ₁−ch)Q³(θ₂−c′h) · Σ_{(mm′),(nn′)} A^{(c,c′)}_{(mm′),(nn′)} minor₁(m,m′) minor₂(n,n′).
- **Cost.** Decompose A^{(c,c′)} = Σ_k ω_k ⊗ ω′_k (matrix rank r_{cc′}; ω, ω′ 2-forms), then each 2-form into wedges
  (rank/2 each). Count = Σ_{c,c′} Σ_k (rank ω_k/2)(rank ω′_k/2) products Bk₂^{(1)}·Bk₂^{(2)}. Report also the
  "separable" count r_{cc′} and the alternative decomposition minimising the total (greedy: pick ω_k decomposable
  where possible — the minimal count is a tensor-rank problem; use exact rank bounds, not heuristics, for claims).
- Regression targets with KNOWN Q-forms: (i) ⟨Ω_w⊗Ω_w| = det₆/Ψ_Q — its Laplace expansion along the two ladders
  gives 20 products of per-ladder 3×3 bracket determinants, each 3×3 = Σ_c (2-form)⊗δ_c; so the A^{(c,c′)} of the
  vacuum must have a structure consistent with that (derive the expected ranks before computing: they are an upper
  bound and a check of the bookkeeping); (ii) the §5 product point ⟨Ω_{(w₁,y*)}ρ(D)⊗Ω_{(w₁,y*)}| = det₆×det₂ (a Q-number
  with the ladder-1 insertion y*^{−x/h}: its site-1 measures should carry a second exponential); (iii) the one-site
  results (site 2 = e₁ᵀ): counts 5/7/11 generic, 4/5/9 reduced.

## 2. Stage A — the mirror: a site-2 rotation with e₁ᵀ on site 1 (cheap, do first)

Φ = e₁ᵀ ⊗ e₁ᵀρ(g₂). The SoV construction is not site-symmetric (B is built from the ordered monodromy), so:
(a) which site-1 labels appear in the support (expect labels of the form (λ₁|λ₁ c₁) by analogy, with factors
λ₁^{(c₁)} := Q^k(θ₁−c₁h)/Q^k(θ₁−λ₂h), or something else — record it); (b) the site-2 column tensors and their
automatic counts; (c) the loci in g₂ — is the site-2 best locus again the vacuum family β+αγ = α² (likely, if the
structure is "vacuum on the rotated site"), or does the asymmetry give a different conic (cc(t) with t → −t or
t → t−k?); (d) closed-form measures on the best site-2 locus, their u-form (expect w₁′^{−(u−θ₂)/h}×measure×cross-ladder
pole (u−θ₁ …)). Output: the site-2 analogue of the 5-bracket formula. Compute for [2,1,0] and [3,1,0]; [3,2,0] only
if free.

## 3. Stage B — the two-site K tensor and the general structure (the main computation)

1. Build K_{ν₁ν₂;x₁x₂} exactly (SoV covectors as exact null vectors of B(u) − b_x(u) at a rational u, or at 60 digits
   then rationalised with the error reported, as in the one-site case). [3,1,0]: 225×225, one linear solve per
   covector; time one first. Check: the ν₂ = 1 slice reproduces the one-site K maps; the ν₁ = 1 slice reproduces
   Stage A.
2. Supports: for Φ = c₁⊗c₂ generic, which (x₁,x₂) pairs occur; the structure of the support in terms of the labels
   (does the site-1 label set depend on c₂ and vice versa? — this is the "measure couples the ladders" statement of
   Part B in label language).
3. Blocks A^{(c,c′)}(g₁,g₂) symbolic (bilinear in the two GT covectors, polynomial in (α₁,β₁,γ₁,α₂,β₂,γ₂)); at a generic
   rational point: matrix ranks r_{cc′}, the 2-form ranks of a rank decomposition, the total count. Same for a random
   (non-rotation) pair as control, and for the three regression targets of §1 (their counts and whether the vacuum
   ⟨Ω_w⊗Ω_w| has visibly lower rank than a generic pair — it should, being a single determinant).
4. The family ⟨Ω_w⊗Ω_{w′}| (g_s on the vacuum family, 4 parameters (w₁,w₂,w₁′,w₂′)): counts as functions of the
   parameters (generic; at w′ = w; at common w₁; at w₁′ = w₂ etc.), exact rank drops located by minors (Gröbner in
   4 variables is affordable; 6 variables for the general pair may not be — then slice: fix g₂ generic, vary g₁).
5. The one-site-fixed slices: g₂ on the site-2 best locus, g₁ free: loci in g₁ (3 variables) and the best count;
   and vice versa. Compare with the product of the one-site counts (if the count factorises as (site-1 count)×(site-2
   count) the sites decouple at the bracket level — a strong structural statement; if it is lower, the two sites
   share measures).

## 4. Stage C — explicit measures and the u-form on the best two-site family

For the family with the lowest count among those with generic CG reach (any two distinct generic-enough rotations
reach everything; see Stage D): the wedge decomposition of each block in closed form in the parameters and t
(the one-site machinery: pivot rule, gauge F_i = 0, F_j = −1, G_j = 0; t-dependence by the 10-point fit or by
rebuilding K at symbolic θ if the 225×225 solve is feasible with one symbolic parameter — try `t` symbolic with
h = 1, θ = (t, 0) as in the perturbative notebooks). Then the gauge-invariant kernels A = w₁^{−(u+u′)/h}·a(…) per
site and the class tests (measure × twist exponential × cross-ladder poles × low-degree polynomial; the one-site
result says the λ₂-group columns are natural and the vacuum-group columns are not — here EVERY block has a
non-vacuum label on at least one site when w′ ≠ w, so the expectation is that more blocks become natural; check).
Deliverable: ⟨Ω_w⊗Ω_{w′}|Ψ⟩ = Σ (short) Bk₂^{(1)}[F(u),G(u)] Bk₂^{(2)}[F′(u),G′(u)] Q³(θ₁−ch)Q³(θ₂−c′h) / (Ψ_Q-normalisation)
with explicit functions.

## 5. Stage D — CG reach and the untwist (the payoff)

- Reach at finite twist (exact, as general-notebook Task 6): site 2 = e₁ᵀρ(g₂P⁻¹) = (…)e₁ᵀρ(U₂), U₂ the upper Gauss
  factor of g₂L₀D₀; undone site-1 covector e₁ᵀρ(g₁L₀D₀U₀U₂⁻¹) → weights reached, ranks on double weights, for the
  family of Stage C. Expect complete for generic (w,w′).
- Untwist (`L2-site1-rotation-untwist-spec.md` adapted): z = 1+εφ, w = 1+εω, w′ = 1+εω′ (all → 1 jointly, the
  regime of the det₆×det₂ thread but with w′ free). The ε-scaling question of the one-site plan (Gauss factors of
  P⁻¹ have poles ~ε⁻²) now has the explicit answer from the diagonal-frame picture: in the diagonal frame
  Ω_w⊗Ω_{w′} = e₁ᵀρ(n(z,w))⊗e₁ᵀρ(n(z,w′)) with the rational Gauss parameters r_k(z,w) of `oneSiteVac`, which stay
  O(1) in the joint regime (`L2-210-untwist-CG.md` §1) — so no scaling search is needed for this family: the undone
  covector is e₁ᵀρ(n(w)n(w′)⁻¹), explicitly O(1) and generic at ε = 0. Then: per-weight ε-table (expect all weights
  at ε⁰, rank 2 on doubles), K/S for all HW branches (expect O(ε) for all), CG at leading order with Richardson
  against the chain, the Q-side number from the Stage C bracket formula at ε = 10⁻⁶ (the measures' ε-exponents and
  where the cancellation sits). [2,1,0] first (6 branches, reference numbers in `L2-210-untwist-CG.md` §5), then
  [3,1,0] (10 branches, `L2-310-untwist-CG.md` §4).
- Note the companion-frame ↔ diagonal-frame dictionary needed for the Q side: the Q-functions are frame-independent
  (eigenvalues), ⟨x₀|Ψ⟩ is the companion-frame normalisation; the bracket formula is written in the companion frame
  and the chain-side CG extraction in the diagonal frame — keep both normalisations explicit (the adjoint untwist
  notebook already handles ⟨x₀|Ψ⟩ = 1 with `x0Of`).

## 6. What "success" means and the cost accounting

Success = a family (Ω_w⊗Ω_{w′} or a cheaper relative) with (a) a bracket count polynomial in (λ₁,λ₂) and small in
practice (target ≤ ~25 products for [3,1,0], i.e. ≈ 5×5), (b) explicit measures F(u), G(u) for most blocks,
(c) K/S = O(ε) for all HW branches in the joint regime and all CG data at leading order. Report the count against
dim V² = 225 and against the det₆ Laplace count (20 products of 3×3 bracket determinants ≈ 20×(2+…) in our units).
If the two-site count is simply (one-site count)², the family is "two independent vacua" at the bracket level and the
formula is the product of the two one-site formulas with the site-2 reference replaced — the simplest possible
outcome; check it first at a numeric point before any symbolic work (rank of the 225×225 coefficient matrix φ_{x₁x₂}
vs the product of the one-site supports is NOT the right test (Part B: full rank 8); the block ranks r_{cc′} are).

## 7. Risks, estimates, order of work

- The 6-parameter Gröbner solves may be infeasible; restrict to the vacuum families (4 parameters) and one-site slices.
- Exact 225-dim SoV expansions: the SoV covectors are exact null vectors of sparse rational matrices; if slow, use
  60-digit numerics and rationalise (the one-site K maps rationalised with denominators ≤ 6 digits).
- The det₆ regression is essential: if the vacuum's blocks do not show the expected low rank, the bookkeeping
  (signs of minors, lone-point conventions on site 2) is off — fix before any locus work.
- Order: Stage A ([2,1,0], [3,1,0]) → Stage B items 1–3 at numeric points → the det₆ regression → Stage B item 4
  (vacuum families) → Stage C → Stage D ([2,1,0] then [3,1,0]). Separate notebooks per stage; the sonnet
  wolfbook-builder for all computations; the kernel-wiping drivers last. Each stage ends with a result markdown cell
  and a one-page write-up; HANDOFF after each.
- Side question to settle in Stage A or B, cheaply: the adjoint's quadric cc = (t−3)/(2(t−2)) = (3/4)κ_singlet — is the
  quadric family a known object in the companion frame (e.g. a vacuum of the transposed/dual chain, or the
  ω₂-line Ω₁+κω₂ of §5.2)? Test: does the quadric contain the ω₂-line Ω_{(w₁,z₂)}-type covectors, and is it the
  image of the vacuum family under the ladder-insertion y^{−x/h} (the det₂ factor of the §5 family)?
