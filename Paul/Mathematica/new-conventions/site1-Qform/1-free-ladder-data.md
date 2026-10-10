# Direction 1: free bra data on both ladders in the per-ladder 6×6

## Idea

The per-ladder 6×6 has two blocks of rows. Rows (i,j; s=1) see only ladder-1 poles, and rows (i,j; s=2) see only ladder-2 poles. Each block may carry its own bra data F⁽ˢ⁾: the values of F₁, F₂ at θ_s−2h, θ_s−h, θ_s, θ_s+h, which is 8 numbers per ladder.
- **Known cases.** (F_w, F_w) gives Ω_w⊗Ω_w, and (F₀, F₀) gives x₀.
- **Ruled out.** (F_w, anything) does not give Ω_w⊗e₁ᵀ (best residual 0.766).
- **What remains.** Perhaps ladder 1 needs data that is **not** the w-vacuum Q's, e.g. vacuum Q's of a shifted or "half" twist, or Q's with a ladder-2 factor stripped. Direction 1 lets both ladders float.

## Structure (exploited for speed)

- **Laplace expansion.** Expanding along the two row blocks,

  det₆ = Σ_{S⊂{1..6}, |S|=3} ± M⁽¹⁾_S · M⁽²⁾_{S^c} = ⟨ω₁, ω₂⟩.

  Here ω_s ∈ Λ³C⁶ (20-dimensional) is the Plücker vector of the ladder-s row block. It depends on the state B through the dual Q's at ladder-s points, and on F⁽ˢ⁾.
- **Dependence on F.** ω_s is quadratic in F₁⁽ˢ⁾, because two rows use F₁, and linear in F₂⁽ˢ⁾. The coefficient tensors (10 quadratic monomials × 4 F₂ slots × 20 Plücker components, per state and ladder) are precomputed once by polarization. After that, any det₆ is a cheap contraction.
- **Gauges, per ladder.** There are three: F₁ → aF₁ (det scales by a²), F₂ → bF₂ (det scales by b), and F₂ → F₂ + cF₁ (det unchanged). The effective parameter count is 16 − 6 = 10.

## Equations and safeguard

- **Equations.** det₆[B; F⁽¹⁾, F⁽²⁾]/Ψ_Q^B(x₀) = T_B for all 64 states B, with T_B = ⟨Ω_w⊗e₁ᵀ|Ψ_B⟩/⟨x₀|Ψ_B⟩, and **λ fixed to 1**. This is legitimate because the F₁ scales absorb it. The trivial branch det ≡ 0 is then excluded automatically, since T ≠ 0.
- **Count.** 64 equations against 10 effective unknowns, so a solution would be strongly non-generic and therefore meaningful.
- **Report** ‖det‖/‖T‖ and the residual for every solution. Discard anything whose ‖det‖/‖T‖ is not ≈ 1.

## Algorithm

1. Precompute the ω_s coefficient tensors for both ladders and all 64 states (polarization: 10 × 4 = 40 minor evaluations per state and ladder).
2. Validate the contraction:
   - (F_w, F_w) reproduces Ω_w⊗Ω_w;
   - (F₀, F₀) reproduces x₀;
   - (F_w, F₀) gives the known "A⊗e₁ + A′⊗e₂" mismatch.
3. Gauss–Newton with λ = 1 on the 16 values. Gauge directions are handled by the pseudo-inverse. Fixed 60-digit arithmetic (`$MinPrecision = $MaxPrecision = 60`), to avoid the earlier precision collapse.
4. Starts: (F_w, F₀), (F₀, F_w), (F_w, F_w) and (F₀, F₀), plus about 40 random starts.
5. **If solved:**
   - gauge-fix the solution;
   - compare F⁽¹⁾ with natural candidates (the w-vacuum Q's times a polynomial in x, Q's of shifted twists, L=1 vacuum Q's at θ₁) and F⁽²⁾ with F₀;
   - redo at w′ and test whether F⁽²⁾ is universal;
   - then test the full site-1 rotation by multiplying with the su(2)-like 4×4, using E on ladder 1 and its x₀-limit on ladder 2.
6. **If not solved:** record the best residual. Then test the relaxation in which F₂⁽¹⁾⊗F₂⁽²⁾ is replaced by a general 4×4 matrix X. That makes the problem linear in X for fixed F₁'s; if it is solvable, the answer is a sum of determinants, which is a bridge to Direction 3.

## Success criteria

- A solution with residual below 10⁻⁴⁰, ‖det‖/‖T‖ = 1, reproduced at w′, and still valid at a second (h,θ).
- Ideally F⁽¹⁾ and F⁽²⁾ have recognizable closed forms.

## Cost and risks

- **Cost.** The precompute takes about 64 states × 2 ladders × 40 cheap 3×3-minor evaluations, i.e. seconds. Each Gauss–Newton step is a 64×16 Jacobian from tensor contractions, also seconds.
- **Risks:**
  - local minima (mitigated by many starts);
  - near-degenerate gauge directions (mitigated by the pseudo-inverse);
  - a negative result, which is still informative because it motivates Directions 2 and 3.

## Outcome (2026-10-02): negative

The code is in `../L2-rotation-sov.wb` §7, cells 42–46.

- **Validation.** The fast determinant matches det6L (10⁻⁵¹), and the adjugate Jacobian matches finite differences (10⁻³³).
- **Free data on both ladders.** All 16 numbers free, λ = 1 fixed. Levenberg–Marquardt from 44 starts, then 300 extra steps on the best 6, plateaus at a minimal relative residual of **0.273**. Every run gives 0.27–0.28, with ‖det‖/‖T‖ ≈ 0.96, so this is not the trivial branch.
- **X-relaxation.** For fixed F₁'s, a general 4×4 X in place of F₂⁽¹⁾⊗F₂⁽²⁾ reaches only a rank-9 space. The residual is 0.269 at the best F₁'s and 0.61–0.92 at the natural ones.
- **Conclusion.** No per-ladder 6×6, and no sum of them sharing the F₁ data, gives ⟨Ω_w⊗e₁ᵀ|Ψ⟩. The bracket and column structure itself has to change (Direction 2), or the answer is a genuine SoV sum (Direction 4).
