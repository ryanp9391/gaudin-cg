# gl(3), L=2, [3,1,0]⊗[3,1,0]: the two-rotation product-state recipe for the next rep

**Date:** 2026-10-09. **Notebook:** `L2-310-two-rotations.wb` (14 cells, self-contained, exact rationals, sparse 225×225; the whole notebook reruns in well under a minute except the mixed-point `Solve`, 29 s). **Builds on:** `Two-rotations-det-overlap.md` (adjoint, §5 and §5.8) and `L2-310-localization.md` (the [3,1,0] localization data). **Question (Paul):** generalize the adjoint construction (v₂ = Ω_{(w₁,y*)}, v₁ = v₂ρ((w₂−Gᵀ)⁻¹), overlap ∝ det₆ × ladder-1 det₂) to [3,1,0].

Conventions: new conventions, plain GT generators (not EH, so everything is exact), Lax (u−θ_s)δ + hE, twist on the right, CBR fusion T_λ(x) = det_{i,j} t_{λ'_j+i−j}(x+(i−1)h) with shapes labelled by row lengths (λ₁,λ₂). Numeric point P1: h = 1/2, θ = (2/7,−3/5), z = (3/2,5/4,2/3), w = (5/2,7/3,3/5); also h = 1/3 and (θ = (1/5,−4/9), z = (3/2,7/6,2/3), w = (7/4,−2/5,9/7)). Diagonal frame; Ω₁ = e₁ᵀe^{r₁E₁₂}e^{r₂E₁₃}e^{r₃E₂₃} with the rep-independent r_k(z,w); Ω⊗Ω is the exact left eigen-covector of t_a(G′) (checked).

**Rotation test (new, exact, linear).** A covector φ of the 15-dim rep is in GL(3)·e₁ᵀ iff the stabilizer of the line [φ] is the 6-dim Borel iff rank{φ, φE_ij (9)} = 4; a generic covector gives 9 (the trace relation removes one), e₁ᵀE₁₂-type (weight (2,2,0)) gives 5, e₁ᵀ(E₁₂+E₂₃)-type 6, e₁ᵀE₁₃-type 7. (`rotRank`, cell 4.) Membership in a given orbit is then settled by an exact `Solve` with a Gauss-decomposed g (cell 9).

## Results

1. **Localization at θ₁−3h** (cell 5): (1,1), (2,2), (3,2) are site-1 local, as known; every other shape with λ₁ ≤ 3 has rank 2 with Ω₁ inside its site-2 span; (3,3) vanishes on ⟨Ω|. The joint site-2 span of all shapes is 2-dimensional, {Ω₁, ω₂}, with ω₂ a rotation and the pencil Ω₁ + κω₂ inside the rotation orbit (cell 6) — the adjoint's structure exactly. ω₂ is not Ω_{(p,q)} for (p,q) among the twist parameters (unlike the adjoint's Ω_{(w₁,z₂)}); it is a torus-rotated vacuum. Shapes with λ₁ = 4 add nothing (cell 8).

2. **Product locus** in the 10-dim span ⟨Ω|(c₀ + Σ_λ c_λ T_λ(θ₁−3h)) (105 quadrics, exact `Solve`, 0.5 s): the localizing span (5 free parameters) and **two isolated points** (cells 6, 7).
   - **Point A.** c = the adjoint's c* = (−2187/1600, −243/80, −27/43, 1, 320/1161) carried verbatim to the shapes shifted by one column: (1,1), (2,1), (2,2), (3,1), (3,2). Site 2 = torus-rotated vacuum with torus invariant y = w₂(θ₁₂−h)/(θ₁₂−2h) (the adjoint's y*), at all three numeric points. Site 1: rotRank 6 — not a rotation, and not in span{e₁ᵀRD, e₁ᵀE₁₂RD, e₁ᵀE₁₂²RD} for the natural R, D (cell 13).
   - **Point B** (all nine shapes and 1; at P1 c = (81/1075, 36/215, 4/39, 1, 860/1053, 4300/9477, −860/729, −92450/85293, −1849000/2302911) on (1, T₁₀, T₁₁, T₂₀, T₂₁, T₂₂, T₃₀, T₃₁, T₃₂)). Site 2 = torus-rotated vacuum with invariant

     **y* = w₂ (θ₁₂−h)/(θ₁₂−3h) = w₂ (θ₁₂−h)/(θ₁₂−λ₁h)**

     at all three numeric points (cell 12; exact rationals −63/43, −203/18, 26/385). Site 1:

     **v₁ ∝ e₁ᵀ(α + βE₁₂) · R(z,(w₁,y*)) · ρ((w₂−Gᵀ)⁻¹)**

     i.e. the adjoint's v₁ = Ω_{(w₁,y*)}ρ(D_{w₂}) with the vacuum covector e₁ᵀ replaced by vacuum + first descendant (cell 13: exactly in that 2-dim span; the other orderings and (w₂,y*)/D_{w₁}, (w₁,w₂) variants fail). β/α = 392/2335 (P1), 294/2255 (h = 1/3), 5400/337687 (third point); no closed form yet. For the adjoint, e₁ᵀ(α+βE₁₂) = α e₁ᵀ e^{(β/α)E₁₂} is itself a rotation because E₁₂² kills e₁ᵀ there — which is why the adjoint's site-1 factor was a rotation; for [3,1,0] the E₁₂² term is missing and v₁ has rotRank 5 (orbit type of e₁ᵀE₁₂; exact orbit test, two solutions, cell 9).

3. **Mixed-point span** (θ₁−3h shapes plus the λ₁ ≤ 2 shapes at θ₁−2h; site-2 span 9-dim, 3549 quadrics, `Solve` 29 s; cells 8, 10): 18 components, no new point with a site-2 factor other than Ω₁. A and B reappear; the only site-1 rotations occur with site 2 = Ω₁ (trivial).

4. **Overlap.** ⟨v₁⊗v₂|Ψ_B⟩ = c·(1, τ_λ(B))(θ₁−3h) ⟨Ω|Ψ_B⟩ automatically (the T's are diagonal on Ψ). The Q-form of this T-polynomial (the [3,1,0] analogue of the ladder-1 contour det₂) is **not done**: it needs the [3,1,0] Q-functions/Baxter machinery and the general-rep vacuum determinant.

5. **Torus.** As for the adjoint (`L2-rotation-sov.wb` §11), y is a torus invariant of the site-2 rotation, so point B gives for every site-2 rotation g exactly one partner e₁ᵀ(α+βE₁₂)RD·ρ(torus); no continuous family.

**Verdict.** "Two distinct rotations of the vacuum with a factorized overlap" does not carry over to [3,1,0] through the θ₁−λ₁h insertions. What carries over, with the same y* pattern (θ₁₂−λ₁h in the denominator), is "rotated vacuum on site 2 × rotated *dressed* vacuum (vacuum + first descendant) on site 1".

## Q-side (cells 15–23; done by an Opus wolfbook-builder subagent from `L2-310-Q-side-spec.md`, reviewed and spot-checked by the main session; everything at P1, 60 digits, 225 states, the whole section reruns in ~10 s)

1. **t's, eigenstates, Q's** (cell 16). t_a(u) of degree (2,4,6) by exact 7-point interpolation; qdet symbolic; simple spectrum (min gap 6.4e-3); eigen-residual 7e-54. All six Q's found for 225/225 states with deg q^j = μ_j, deg q_i = 6 − μ_i (Baxter residual 1.5e-43; 3.4 s).
2. **Vacuum determinant** (cell 17). Vacuum direct q's of degrees (0,4,6). With the rep-generic pole set n = 0..3 on both ladders, insertions 1 and e^{2πix/h}, the same 6×6 columns as the adjoint, and the first-guess reference factor Ψ_Q(x₀) = Π_s [Q²(θ_s−3h)Q³(θ_s−4h) − Q³(θ_s−3h)Q²(θ_s−4h)] Q³(θ_s−h): ⟨Ω⊗Ω|Ψ_B⟩/⟨x₀|Ψ_B⟩ ∝ det₆[Q_B, Q_vac,w]/Ψ_Q, state-independent to 5.1e-45. So the vacuum SoV triple of [3,1,0] is (θ−4h, θ−3h, θ−h), as guessed from the GT-pattern dictionary x₁ = −(m₁₁+1)h, x₂₁ = −m₂₁h, x₂₂ = −m₂₂h.
3. **Ladder-1 determinant** (cells 18–20). With one-site weights on n = 0..3, E = y^{−(x−θ₁)/h}, f_y(B) = det_{j∈{3,2}}[−⟨Q^jE⟩₁, ⟨Q^{j[−]}E⟩₁]/minor₁, minor₁ = Q²(θ₁−3h)Q³(θ₁−4h) − Q³(θ₁−3h)Q²(θ₁−4h): f_y ∈ span{1, τ_λ(θ₁−3h)} for every y (residuals ~1e-47; T₍₃,₃₎(θ₁−3h) ≡ 0 so the basis has 9 functions; controls: random combination 1e-54, a non-member 0.9). With symbolic y the coefficient curve rationalizes exactly and is **graded by the box number**:

   c_λ(y) = a_λ y^{6−|λ|}, |λ| = λ₁+λ₂, with a = (16/9, −2240/387, 78400/15093, 19600/387, −2744000/45279, 96040000/1946997, 2744000/31347, −48020000/407511, 6722800000/52568919) on (1, T₁₀, T₁₁, T₂₀, T₂₁, T₂₂, T₃₀, T₃₁, T₃₂).

   The only nonzero y with c(y) ∥ c_B is **y = −63/43 = y*** (exact Solve); there c(y*)/c(y*)_{T₂₀} = c_B exactly, κ = 34306448400/147008443. No y gives c(y) ∥ c_A: point A (which vanishes on 1, T₁₀, T₂₀, T₃₀) is not on this curve. Review (cell 23): fresh y = 2, −5/3 fit to 1e-46 and reproduce the grading and the a_λ exactly.
4. **End-to-end for point B** (cell 21): with the (w₁,w₂) vacuum Q's,

   ⟨v₁⊗v₂|Ψ_B⟩/⟨x₀|Ψ_B⟩ ∝ det₆N[Q_B, Q_vac,(w₁,w₂)]/Ψ_Q(x₀) × det₂^{(1)}[Q_B; y*]/minor₁, y* = w₂(θ₁₂−h)/(θ₁₂−3h),

   state-independent to 5.1e-45 on all 225 states; ratio4/ratio6 = 1/κ to 4.6e-47; wrong y's (−63/8, 2) fail by 77 and 322.

**So the adjoint formula carries over to [3,1,0] with only rep-generic shifts** (poles n = 0..λ₁, vacuum triple (θ−4h, θ−3h, θ−h), minor₁ at (θ₁−3h, θ₁−4h), y* = w₂(θ₁₂−h)/(θ₁₂−λ₁h)), even though the site-1 factor is a dressed vacuum rather than a rotation. For the adjoint the analogous curve was f_y = 4y⁴ − 70y³τ₁₀ + y²(…τ₁₁ … τ₂₀) + … yτ₂₁: also graded by box number (y^{4−|λ|}) — the [3,1,0] curve is the y^{2λ₁−|λ|} version.

## Open

- Closed form of β/α(h, θ, z, w) (three data points above; the E₁₂ normalization of the GT basis enters).
- Q-side at the other two numeric points (cell 12's); closed form of the a_λ (FactorInteger in cell 20) — compare with the adjoint's 2h²(θ₁₂−h)(θ₁₂−2h)(θ₁₂−3h) y⁴ + … monomials.
- Point A's site-1 structure (rotRank 6).
- A meaning for the dressing: e₁ᵀ(α+βE₁₂) = α e₁ᵀe^{(β/α)E₁₂} − (β²/2α) e₁ᵀE₁₂² is a rotation minus the extremal (1,3,0) covector (itself a rotation), i.e. a sum of two rotations.
- [3,2,0] and [4,2,0] with the same notebook (change Λ in cell 2; the localization data are in `L2-320-localization.md`, `L2-420-localization.md`).

## Log

- Notebook created on a kernel that turned out to be the shared K1: cell 2's `ClearAll` wiped the live state of `L2-rotation-sov.wb` (its file is intact; rebuild = cells 2–21 etc.).
- Cells 2–6: setup, t's (0.2 s, exact sparse), vacuum, rotation test, localization table, product locus (3 components).
- Cell 7: factors of the two points; the adjoint's c* recognized on the shifted shapes.
- Cell 8: enlarged spans S1 (λ₁ ≤ 4) and S2 (mixed points); the cell reported "running" long after it had finished (status bookkeeping), aborts found an idle kernel.
- Cell 9: torus invariants of the site-2 factors (−63/8 and −63/43) and the orbit tests (site 1 of B = e₁ᵀE₁₂-type).
- Cell 10: S2 classification with solutions kept (29 s Solve).
- Cell 11: hypotheses H1–H7 for site 1 of B all negative (pure descendant forms).
- Cell 12: three-point pipeline; first version failed on the `y*)` comment closer (second time today; memory note exists).
- Cell 13: span tests — positive for span{e₁ᵀRD, e₁ᵀE₁₂RD} with (w₁,y*), D = Dy[w₂].
- Cell 14: result cell. Saved; not committed.
- Cells 15–22 (subagent, Opus, ~11 min, 155k tokens): Q-side per `L2-310-Q-side-spec.md`; cell 23 review check by the main session. Saved; not committed.
