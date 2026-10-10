# Global and local rotations of the SoV vacuum as Q-determinants (gl(3), adjoint, L=1 and L=2)

This is the Markdown companion of `Rotations-SoV-Q.tex`/`.pdf`, which has the full write-up from 2026-10-02.

**Notebooks:**
- `Rotation-e8-SoV.wb` (L=1, Sections 1–10e)
- `L2-rotation-sov.wb` (L=2, Sections 1–5b)

**Numeric points:**
- L=1: h=1/2, θ=2/7, z=(3/2,5/4,2/3), w=(5/2,7/3,3/5).
- L=2: the same h, z, w with θ=(2/7,−3/5), so t=θ₁₂/h=62/35. A second twist w′=(7/4,−2/5,9/7) is used once.
- Caution: z₃=1/z₁ at these points.

## L=1

| # | Result | Status |
|---|---|---|
| 1 | The companion-frame SoV covectors ⟨x_μ\| are constants (z-free, θ-free), listed in Section 1 of the PDF | exact |
| 2 | Σc_μ⟨x_μ\| = s·e₈ᵀρ(g) ⇔ c = c(a,b), bilinear with b·a=0, where a=g⁻¹e₁ and b=e₃ᵀg. Equivalently, X(c) is rank 1 and traceless | exact |
| 3 | A single det N[Q_B,F]/Ψ_Q(x₀) is a rotation only on two components. (I) b∝m(β), with m(x)=(x²,x,1): these are the companion eigen-bras, with one twist eigenvalue α free. (II_k) a⟂m(z_k): degenerate, they kill the μ_k=0 kets. This holds for all 72 working row choices | numeric |
| 4 | Every rotation is a sum of 2 det N's from one companion twist | exact |
| 5 | t₁(u)−(u−θ)χ₁ = hρ_Lie(G(z)ᵀ). The torus ρ(p(Gᵀ)) acts on the eigenkets by ∏p(z_i)^{μ_i}. Every rotation is φ_I·ρ(p(Gᵀ)), with φ_I in component (I) | exact / numeric |
| 6 | ⟨x₀\|M(w)M(y) is a rotation only for y₂=∞ (any y₁), or for y₁=z_k (degenerate) | complete solve |
| 7 | **Resolvent identity**, with s=1/y₁ and T's at θ−2h: ρ((1−sG(z)ᵀ)⁻¹) ∝ 1 + (s/h)T₁₀ + (s²/2h²)(T₁₁+T₂₀) + (s³/3h³)T₂₀T₁₀⁻¹T₁₁ | 10⁻⁶⁰, two z's |
| 8 | (w₁,w₂,y₁) parametrize generic rotations of x₀: x₀ρ(g) ∝ ⟨x₀\|M(w)M(y₁,∞) | Jacobian rank 3; random targets to 10⁻³⁸ |
| 9 | As y₂→∞, F₂ → y₁^{−x/h}/Γ((θ−x)/h) (not a delta) | numeric, O(1/y₂) |
| 10 | **su(2)-like factor:** ∏(1−z_i/y₁)^{−μ_i} ∝ det₂[−⟨Q^jE⟩, ⟨Q^{j[−2]}E⟩]_{j=3,2} / (Q²(θ−2h)Q³(θ−3h)−Q³(θ−2h)Q²(θ−3h)), with E=y₁^{−x/h} | 10⁻⁵⁵ |
| 11 | **Rotation = det₃[Q_B,Q_vac,w]/Ψ_Q × su(2)-like det₂** | follows from 8–10 |

## L=2

| # | Result | Status |
|---|---|---|
| 12 | **Global rotation:** overlap/⟨x₀\|Ψ⟩ ∝ [det₆N(Q_B,Q_vac,w)/Ψ_Q] × [det₄ over (Q³,Q²) × 2 brackets, with columns {−⟨PE⟩,⟨P^{[−2]}E⟩}×{1,x}] / ∏_s minor₂ | 64 states, 10⁻⁵⁰ |
| 13 | The su(2)-like factor is 2L×2L. No single-Q 2×2 works | best residual 0.87 |
| 14 | Site-1 t-products do not localize on a⊗x₀. Only T₁₁ and T₂₁ at θ₁−2h localize, and only on x₀ | exact |
| 15 | The untwisted L=2 B has a simple spectrum and an 8×8 grid of SoV states. The SoV wavefunction is ∏_s det₂[v(x_{s21}),v(x_{s1})]v_k(x_{s22}) normalized at the vacuum | 64×64, 10⁻⁴⁸ |
| 16 | **Site-1 rotation** (e₁ᵀρ(g))⊗e₁ᵀ lives on 12 SoV states: site 2 at the vacuum, or at {2\|20} with factor λ₂=v(θ₂)/v(θ₂−h). Overlap/⟨x₀\|Ψ⟩ = Σ_ν c^{L=1}_ν(g)W_ν, with a rational, g-independent K (table in the PDF) | exact |
| 17 | SoV-vacuum bra data F₀=(F₁=δ_{θ_s−2h}, F₂=δ_{θ_s}) reproduces x₀ | 10⁻⁵¹ |
| 18 | **RETRACTED:** the claimed single per-ladder det₆ for ⟨Ω_w⊗e₁ᵀ\|Ψ⟩. Those fits had λ≈0 and det≡0, the trivial branch. Redone with the scale fixed: the best residual is 0.766 (L=2 vac Q's on ladder 1) or 0.774 (L=1 vac Q's), so the natural ansatz is ruled out | redo: global scan + refinement |
| 19 | Naive F⁽²⁾=F₀ fails: it gives A⊗e₁+A′⊗e₂ | numeric |

## 2026-10-09 continuation

See `Two-rotations-det-overlap.md` (full write-up + log) and `site1-Qform/5-bracket-decomposition.md`, `site1-Qform/6-two-rotations-2603-trick.md`. Main new result: with y* = w₂(θ₁₂−h)/(θ₁₂−2h), v₂ = Ω_{(w₁,y*)}, v₁ = v₂ρ((w₂−G(z)ᵀ)⁻¹) are two distinct rotations and ⟨v₁⊗v₂|Ψ_B⟩/⟨x₀|Ψ_B⟩ ∝ det₆N[Q_B,Q_vac,(w₁,w₂)]/Ψ_Q × ladder-1 su(2)-like det₂[Q_B;y*]/minor₁ (10⁻⁵⁰, 64 states). Item 1 below (site-1 vacuum as a single det) is answered negatively for the natural ladder-1 bracket families.

## Open

1. A Q-form for Ω_w⊗e₁ᵀ. The natural per-ladder det₆ is ruled out (item 18). Try modified ladder-1 data, mixed brackets, or a sum of dets. The 12-state SoV sum (item 16) is exact.

**Lesson:** in det-matching fits, fix the scale, or check \|det\|≠0, before trusting a small residual.
2. The full site-1 rotation: add the resolvent / su(2)-like part to item 18.
3. K in closed form in t (needs a second numeric point).
4. A second (h,θ) point for all L=2 statements, and the explicit normalization constants.
5. e₃=∏z_i^{μ_i} as a natural insertion.
