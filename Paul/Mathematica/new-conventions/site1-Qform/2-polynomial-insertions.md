# Direction 2: polynomial insertions that kill ladder-2 poles

## Idea (Paul's original guess)

Keep the natural Q-data: the w-vacuum direct Q's for the site-1 vacuum, and E = y₁^{−x/h} for the su(2)-like factor. Insert polynomials π(x) into the bracket integrands. The insertions modify or kill some poles of the measure ρ = ∏_{s,n}(x−θ_s+nh)⁻¹, so that only the "correct" SoV states contribute.
- **Target structure.** For the full site-1 rotation:

  det₆^{π}[Q_B; F_w] × det₄^{π′}[Q³,Q²; E] / (Ψ_Q · ∏minor₂)

- **Why it is genuinely new.** It is not equivalent to changing F, which was Direction 1. In a bracket with insertion π, the columns pick up:
  - −⟨PF⟩ → Σ_p w_p π(x_p) P(x_p) F(x_p);
  - ⟨PF^{[2]}⟩ → Σ_p w_p π(x_p) P(x_p) F(x_p+h);
  - ⟨P^{[−2]}F⟩ → Σ_p w_p π(x_p) P(x_p−h) F(x_p).

  Writing G = πF changes the first and third columns, but the second picks up π(x)/π(x+h). So insertions add freedom that bra data cannot mimic.

## Ansatz family (staged, smallest first)

1. **Per-pole kill masks.**
   - In the all-pole brackets (insertions 1 and e^{2πix/h}), multiply the ladder-2 poles by 0/1 masks per bracket k.
   - The "polynomial" is then any π vanishing on the masked ladder-2 poles. The values at the unmasked poles are free but mostly gauge.
   - Enumerate: 2 brackets × 2³ subsets of the ladder-2 poles {θ₂, θ₂−h, θ₂−2h}, giving 64 combinations, and test each against T. This is cheap.
2. **Low-degree π with free coefficients.**
   - π_k(x) = Σ_{d≤3} a_{k,d} x^d, a separate polynomial per bracket k, possibly a separate one per column type.
   - The det is polynomial in the a's (degree ≤ 6 overall). Fit with the scale fixed: normalize one coefficient to 1 and fix λ = 1.
   - Unknowns: about 8–24. Equations: 64.
3. **Polynomials tied to the ladder-2 roots.** For example π(x) = (x−θ₂)(x−θ₂+h), or factors (x−θ₂+nh) matching the site-2 vacuum roots {θ₂−3h, θ₂−2h, θ₂−h}. These are the natural guesses for "only the correct states survive".

## Two-step strategy

- **(a) The site-1 vacuum first.**
  - Find insertions such that det₆^{π}[Q_B; F_w]/Ψ_Q ∝ T_B, with T_B the Ω_w⊗e₁ᵀ target.
  - Sanity requirement: the same insertion, with F_w → F₀, must still give x₀.
- **(b) The su(2)-like part.**
  - Find π′ such that det₄^{π′}[E]/∏minor₂ equals the site-1 resolvent factor, i.e. the extra factor that turns Ω_w⊗e₁ᵀ into the site-1 rotation.
  - This factor is not ∏(1−z_i/y₁)^{−μ_i}, because the site-1 resolvent is not a global symmetry. So step (b) needs its own target: T^{rot}_B / T_B.

## Safeguards

- Fix the scale. Report ‖det‖/‖T‖.
- The kill masks can make det ≡ 0 identically. Discard those before fitting.
- Positive results must hold at w′ and at a second (h,θ).

## Diagnostics if it fails

- From the exact 12-state SoV expansion: which SoV configurations does each masked det₆ populate? Compare with the target's 12. This tells which poles must be killed or weighted.
- The K table (`../Rotations-SoV-Q.pdf` §8) gives the required weights of the λ₂-shifted states. Matching them may dictate π.

## Cost

- Stage 1: 64 combinations × 64 states × one 6×6 det each. Seconds.
- Stage 2: Gauss–Newton in a few to tens of unknowns. Minutes.
- Stage 3: closed-form guesses. Trivial.

## Outcome (2026-10-02): negative

The code is in `../L2-rotation-sov.wb` §8, cells 47–53.

- **Framework.** General per-pole, per-column weights W[k,f,p] (36 numbers) cover every polynomial insertion, including column-dependent ones. The natural W reproduces Ω⊗Ω and x₀ (10⁻⁵⁰).
- **Stage 1, kill masks.** All 64 masks fail (best ratio spread 3.5).
- **Stage 2, w only.** The residual drifts down to 0.036 but is overfitting: the same W gives 0.996 at w′ and 0.41 for x₀.
- **Stage 2, joint fit** at w and a second real twist w″ = (9/4,4/3,5/7), with λ = 1 at w and a projective comparison at w″. It converges to a firm minimum of **0.121**, with the weights bounded.
- **Conclusion.** There is no w-independent insertion in the 6×6 bracket determinant with natural bra data.
