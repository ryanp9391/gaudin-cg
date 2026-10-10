# Direction 6: two distinct rotations with a contour-integral determinant overlap (the 2603.26617 trick for the adjoint)

**Status:** POSITIVE, 2026-10-09. Code: `../L2-rotation-sov.wb` §10 (cells 86–103). Verified at P1 (h=1/2, θ=(2/7,−3/5), z=(3/2,5/4,2/3), w=(5/2,7/3,3/5)); the structure and y* also at four more (h,θ,z,w) points.

## Result (gl(3), Λ=(2,1,0), L=2)

With Ω_{(p,q)} the companion-frame one-site vacuum e₁ᵀP(p,q), G(z) the companion ket twist, and

    y* = w₂ (θ₁₂ − h)/(θ₁₂ − 2h),

the product covector

    v₂ = Ω_{(w₁,y*)},   v₁ = Ω_{(w₁,y*)} ρ((w₂ − G(z)ᵀ)⁻¹) = Ω_{(w₂,y*)} ρ((w₁ − G(z)ᵀ)⁻¹)

(two distinct exact rotations of the SoV vacuum) has

    ⟨v₁⊗v₂|Ψ_B⟩/⟨x₀|Ψ_B⟩ = c · det₆N[Q_B, Q_vac,(w₁,w₂)]/Ψ_Q(x₀) · det₂[−⟨Q^jE⟩₁, ⟨Q^{j[−]}E⟩₁]_{j=3,2} / (Q²(θ₁−2h)Q³(θ₁−3h) − Q³(θ₁−2h)Q²(θ₁−3h)),

E = y*^{−x/h}, ⟨·⟩₁ = bracket over the ladder-1 poles θ₁−nh (n=0,1,2) with the one-site weights. c is state-independent to 10⁻⁵⁰ on all 64 states. This is the ⟨W_{ω,κ}⊗W_ω| structure of 2603.26617 App. G (global vacuum of twist (w₁,y*), site-1 resolvent rotation at w₂), but the global det carries the (w₁,w₂) vacuum Q's and the ladder det carries y*; the paper-like assignment det₆(w₁,y*)×det₂(w₂) fails (60.8).

## Mechanism

1. T₁₁, T₂₁ at θ₁−2h localize exactly on ⟨Ω_w⊗Ω_w|; T₁₀, T₂₀ do not, but ⟨Ω|T_λ(θ₁−2h) = v_λ⊗Ω₁ + v′_λ⊗ω₂ with ONE extra site-2 direction ω₂ = Ω_{(w₁,z₂)}, and the line Ω₁+κω₂ lies inside the rotation cone (exact).
2. In ⟨Ω|(c₀+c₁T₁₀+c₂T₁₁+c₃T₂₀+c₄T₂₁) the product locus (exact Solve) has one component with c₁,c₃ ≠ 0; on it the site-1 factor is automatically a rotation. It is z-independent: c* = (−3⁷/40², −3⁵/80, −27/43, 1, 320/(27·43)) at (h,θ)=P1.
3. The ladder-1 su(2)-like contour det f_y lies in span{1,τ₁₀,τ₁₁,τ₂₀,τ₂₁}(θ₁−2h) for every y (f_y = 4y⁴ − 70y³τ₁₀ + y²(114τ₁₁ − 181τ₂₀ − 2τ₂₁) + 394yτ₂₁ at P1); its coefficient curve hits c* exactly at y = y*, where Φ_y drops from rank 2 to rank 1.
4. y* decoded from five (h,θ,z,w) points: −y*/w₂ = 27/8, −58/23, 13/32, … = −(t−1)/(t−2), t = θ₁₂/h. v₁ identified by a scan over Ω_{(p,q)}ρ(M) (residual 0).

## Negative side results (same session)

- Enlarged linear family (adding t₂(θ₁), t₂(θ₂+h), t₂(θ₂+h)t₁(θ₂−h)): no continuous family of products of rotations; isolated points only (Ω₁T₂₁⊗Ω₁ with overlap τ₂₁det₆; ⟨Ω|(e₁+e₈t₂(θ₂+h)t₁(θ₂−h)) at a fixed ratio).
- The naive two-twist state Ω_w⊗Ω_{w′} is neither a det₆ with ladder-dependent vacuum data (dev 5–25) nor has a factorized SoV bra wavefunction (dev 1.6).

## Next

- Second numeric point for the determinant identity itself.
- Dependence on w₃ and on a finite second resolvent parameter y₂; the site-2 analogue (resolvent on site 2).
- General reps / L>2, where the localizing set is larger (unified rule T_(k,k), T_(k,d) at θ₁−λ₁h) and a continuous κ-family may reappear.

## Update (2026-10-09, later): beyond the vacuum — `../L2-rotation-sov.wb` §11 (cells 107–112), write-up §5.8

- The ket-twist torus ρ(D)⊗ρ(D) is the weight character χ_μ(D) = Π d_i^{μ_i} on every Ψ_B (exact), so the family comes with its torus orbit for free.
- w₁ is a torus direction and y a torus invariant of the rotation (exact Jacobian ranks 3/4/3; Ω_{(w₁,y)}/Ω_{(w₁′,y)} is a torus character). The family is therefore V = {(g ρ((w₂(y(g))−Gᵀ)⁻¹), g) : g any rotation of x₀}, 3-dimensional: every site-2 rotation (vacuum or not) has exactly one site-1 partner, with overlap ∝ χ_μ(D) det₆[(w₁, w₂(y(g)))] det₂^{(1)}(y(g)). Verified for a random g (10⁻⁵¹, 64 states) and for the family at fresh w₁ (10⁻⁵¹, 10⁻⁵⁰).
- No continuous κ-family even with the torus; the §10 step does not iterate on v₁⊗v₂ (site-1 T's at θ₁−2h have full site-2 column space; product locus = seed only).
- Open: closed form of y(g); site-2 analogue; joint site-1 × site-2 span; general reps.
