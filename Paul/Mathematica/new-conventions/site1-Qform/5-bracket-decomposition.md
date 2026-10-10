# Direction 5: term-by-term bracket decomposition of the 12-state SoV sum

**Status:** done 2026-10-09, negative for single determinants; see Outcome. Code: `../L2-rotation-sov.wb` §8 (cells 59–78).

## Premise

The site-1 rotation Φ = (e₁ᵀρ(g))⊗e₁ᵀ is exactly the 12-state SoV sum of §4,
⟨Φ|Ψ_B⟩/⟨x₀|Ψ_B⟩ = Σ_ν c_ν Σ_x K_{νx} ψ⁽¹⁾_{x₁}(B) λ₂(B)^{[x₂≠vac]}.
The handoff's negative results are for *fitted* determinant ansätze. Here the exact sum is read
term by term instead, looking for bracket structure directly.

## Variables that remove the mess

On ladder 1 the dual Q's satisfy v(θ₁) ∥ v(θ₁−h) and v(θ₁−3h), v(θ₁−2h), v(θ₁−h) span a 2-plane. Per state:
- λ_s = v_k(θ_s)/v_k(θ_s−h), s = 1,2 (k-independent);
- v(θ₁−h) = α v(θ₁−3h) + β v(θ₁−2h).

Then the 8 site-1 wavefunctions (normalized at the vacuum triple (−3,−2,−1)) are
{1, β, −α} × {1, λ₁, λ₁²}, and the overlap is a short polynomial in (α, β, λ₁, λ₂), linear in c_ν.

## Steps

1. Rerun setup cells 2–21 (fresh kernel). Keep the 8 c_ν symbolic (general site-1 covector), impose c = c(a,b) only at the end.
2. Verify the (α,β,λ₁,λ₂) polynomial against the exact overlap for all 64 states (10⁻⁴⁰ or better).
3. Split by site-2 factor: the Q¹(θ₂) group (4 SoV states) and the Q¹(θ₂−h) group (8 states).
4. In each group, split by the third-column point (Q¹(θ₁) vs Q¹(θ₁−h)); the remaining factor is a 2-form on the points {θ₁−3h, θ₁−2h, θ₁−h}. Decompose it as F₁∧F₂ (Plücker coordinates).
5. Compare the F-pairs across sub-terms and groups. Shared bra data ⇒ a determinant; otherwise record which term breaks it.
6. Impose c = c(a,b), map g ↔ (w, y₁), and compare the F's with the w-vacuum Q's and y₁^{−x/h}.
7. If a determinant emerges: test on all 64 states at a second g, then at a second (h, θ).

**Control:** steps 2–5 at L=1 first (same 8 triples, K = identity, answer known: det₃ × det₂).

## Safeguards

- No free overall scale in any comparison; report ‖lhs‖/‖rhs‖.
- Decomposability of a 2-form on 3 points is automatic; only consistency across terms counts as a result.

## Outcome (2026-10-09)

All numbers at P1 (h=1/2, θ=(2/7,−3/5), z, w as in the notebook), 60 digits.

1. **The polynomial (exact, 10⁻⁴⁷ on all 64 states, any site-1 covector c).** With ψ⁽¹⁾ = {1, β, −α, λ₁, λ₁β, λ₁²β, −λ₁α, −λ₁²α} (trip order),

   | site-2 | λ₁ power | coefficient |
   |---|---|---|
   | λ₂⁰ | 1 | c₁ + (8β/27)c₂ + (4α/9)c₃ |
   | λ₂⁰ | λ₁ | −(35c₂+62c₄)/124 − (4β/27)c₅ − α(35c₆+24c₇)/108 |
   | λ₂⁰ | λ₁² | (β/31)c₆ + (3α/62)c₈ |
   | λ₂¹ | 1 | −(3395/3348)c₂ |
   | λ₂¹ | λ₁ | (35/54)c₅ − (35β/124)c₆ + (35α/558)c₈ |

   So the overlap is Σ_x φ⁽¹⁾_x ψ_x + λ₂ Σ_x φ⁽²⁾_x ψ_x with φ⁽ⁱ⁾ = c·K⁽ⁱ⁾, K^vac of rank 8, K^20 of rank 4. (The 31s and 35s are t = 62/35.)

2. **Neither group is a rotation.** The equivariant map X(φ) ∈ sl₃ was rebuilt in this notebook's conventions (X(e₁)=3e₃₁, X(rotation) exactly rank 1). For five generic g, X(c·K^vac) and X(c·K^20) have three nonzero singular values (e.g. 2.9, 0.19, 0.007). So neither site-2 group is an L=1-type determinant with transformed bra data.

3. **The L=1 det₃×det₂ formula is not a polynomial identity.** Fed the L=2 dual Q's on ladder 1, it is not in the span of the 8 site-1 wavefunctions (relative residual 0.11 with either the L=1 or the L=2 w-vacuum Q's). At L=1 every function of the 8 states is in that span, which is why the formula works there. Also the ladder-1 su(2)-type det N_∞(y₁) is not the λ₂⁰ group of the site-1 resolvent rotation (0.93, 0.77).

4. **Locality criterion and the local family (symbolic).** A ladder-1 bracket determinant is a site-1 object (site 2 at the SoV vacuum) only if it reduces to a polynomial in (α,β,λ₁). For the L=1-type 3×3 (rows (1,3),(1,2),(2,3), columns −⟨Q^jF_i⟩, ⟨Q^jF_i^{[+]}⟩, ⟨Q^{j[−]}F_i⟩) the exact condition is **F₂ ∝ F₁ on {θ₁−2h, θ₁−h}** (other components have F₁ vanishing on points). With the gauge F₂ → F₂ + cF₁ (det invariant), the family is F₁ free on four points and F₂ on {θ₁, θ₁+h}; its image in the 8-dim ψ-coefficient space is only a **4-dim cone**, explicitly
   φ = (16f₋₂²p, −32f₋₂f₋₁p, 48f₋₁²p, 8f₋₂f₋₁p − 8f₋₂²q, 16f₋₂f₋₁q, 8f₋₂f₁p − 8f₋₂f₀q, −8f₋₂f₁p − 8(4f₋₁² − f₋₂f₀)q, …),
   with the ideal of 10 quadrics (e.g. φ₂φ₆ = φ₁φ₈, 4φ₃φ₆ = 3φ₂φ₈, φ₆² + φ₆φ₇ = φ₅φ₈). With three independent bra functions the in-span image is still at most 4-dim (32 components, max rank 4).

5. **Exact membership (Gröbner).** Control member: 0. Random vector: 0.41. Group 1 of generic rotations: 0.16, 0.23, 3.9. Group 2: 2.6, 0.18, 0.40. Numeric fits agree (best 0.18 and 0.014 respectively; an exact member fits to 10⁻¹⁷). **No single local 3×3 bracket determinant of this type reproduces either group.** The "Cauchy–Binet product" det₂[⟨Q^jF_i⟩]×⟨Q¹G⟩ (image 6-dim) fails already for L=1 rotations (2.4).

6. **Sums of two local determinants cover everything.** Fitting a sum of two family members: group 1 and group 2 of two rotations, a random 8-vector and a control all fit to 10⁻¹⁵–10⁻²². So every site-1 covector overlap is (D₁+D₂) + λ₂(D₃+D₄) with four local 3×3 bracket determinants, and this says nothing about rotations.

**Reading.** In the 2-plane variables the natural ladder-1 bracket determinants are tiny families (4-dim cones), while the site-1 rotation's two groups are generic in the 8-dim local space. Any single-determinant form must therefore use a different bracket structure: brackets involving v(θ₁+h) (outside the 2-plane), both ladders with insertions that produce λ₂, or larger determinants whose locality is enforced by cancellation. The exact representation remains the 12-state sum, now in the compact form of item 1.

## Follow-up (2026-10-09, later): contour-integral form and two different rotations (`L2-rotation-sov.wb` §9)

Paul's fallback target: Φ = x₀ρ(g₁) ⊗ x₀ρ(g₂) with g₁ ≠ g₂, in genuine contour-integral form (brackets with the full SoV measure). Step 1 = the two vacua Ω_w ⊗ Ω_{w′} (w′ = (7/4, −2/5, 9/7)).

1. **det₆ with ladder-dependent bra data fails.** Full-measure det₆ (insertions 1, e^{2πix/h}) with F = w-vacuum Q's on the ladder-1 poles and w′-vacuum Q's on the ladder-2 poles: max relative deviation of the ratio to the exact overlap 24.9 (L=2 vacuum Q's), 5.8 (one-site vacuum Q's), 6.1 and similar for the mixed choices. Sanity: the w′-vacuum Q's reproduce Ω_{w′}⊗Ω_{w′} to 10⁻⁴⁸. (The per-ladder det₆ is the same determinant up to the constant (e^{2πiθ₂/h} − e^{2πiθ₁/h})³, since the e^{2πix/h} insertion is constant on each ladder.)

2. **The bra wavefunction does not factorize over the sites.** With M_Φ(x) = c_x N_x the SoV coefficients (so overlap/⟨x₀|Ψ⟩ = Σ_x M_Φ(x) ψ_B(x)), a factorized bra wavefunction Φ_{ab}(x) = Φ¹_a(x₁)Φ²_b(x₂) would give M_{ww′}M_{w′w} = M_{ww}M_{w′w′} pointwise, independently of the SoV measure. Measured: relative deviation 1.6; the ratio M_{ww′}/M_{ww} varies by a factor 46 over the site-1 triple. So product bras do not have product wavefunctions in the B-eigenbasis; a "sum over SoV integrals with a factorized bra" does not exist for g₁ ≠ g₂ either. This is the same obstruction as the rank-2 structure of the site-1 rotation (which is the w₂ → ∞ limit).

**Net:** for site-local rotations the only exact representations found are the lattice sums with the explicit coefficient tables (12 states for site 1; 64 states for g₁⊗g₂), i.e. a bracket of a non-factorizing multi-point kernel.
