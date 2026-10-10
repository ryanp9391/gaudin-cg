# gl(3) XXX, L=2, [3,1,0]⊗[3,1,0]: site-localizing transfer matrices and CG reconstruction

These notes cover `L2-310-localization-EH.wb`, from 2026-10-01.
- **What it repeats.** Sections 10 and 12 of `L2-diagonal-frame-overlap-EH.wb` (see `L2-diagonal-frame-overlap.md`), for the rep [3,1,0] (dim 15, 225 states).
- **Running it.** The notebook is self-contained, clears the kernel, and reruns in about 25 s. All 18 checks pass.
- **Conventions.** The new conventions with EH generators. The t's carry the ket twist diag(z), and the bra is ⟨Ω| = ⟨Ω_zw|⊗⟨Ω_zw|.
- **Numeric point.** h = 1/2, θ = (2/7,−3/5), z = (3/2,5/4,2/3), w = (5/2,7/3,3/5).
- **Generic in the rep.** The code takes the rep Λ as a parameter, and [2,1,0] is run as a control. The control reproduces exactly the known objects: T₍₁,₁₎ and T₍₂,₁₎ at θ₁−2h on site 1, and t₂(θ₂+h) and t₂(θ₂+h)t₁(θ₂−h) on site 2.

## Method
- **Locality.** A covector is a d×d matrix C, with rows indexing site 1. Site-1 local means C = vΩ_zwᵀ; site-2 local means C = Ω_zw vᵀ.
- **Search** (machine precision, threshold 10⁻⁸):
  - every CBR shape with λ₁ ≤ 5 and λ₃ = 0, at u = θ_s + ch for c = −8..8;
  - a breadth-first step, multiplying each local covector by t₁ and t₂ at θ_s+ℤh, to depth 2;
  - an unbiased scan of all products of two and of three t's.
- **Why this covers t\*.** The fused t\*'s are CBR shapes of t's at shifted points: t\*_a(v) = t_{3−a}(v+2h)/t₃(v+2h) gives the complementary shape.
- **Confirmation.** Every hit is re-checked exactly (surds included), using covector–Lax products C·(ℒ⊗ℒ) = ℒ⁽¹⁾ᵀCℒ⁽²⁾.

## Results (exact at the numeric point; coefficients fitted at this one point)

**Site 1.** Everything is at **θ₁−3h** (= θ₁−λ₁h; for [2,1,0] it was θ₁−2h).
- **What localizes:** **T₍₁,₁₎ = t₂, T₍₂,₂₎ and T₍₃,₂₎**.
- **What does not:** T₍₃,₁₎ and T₍₂,₁₎ do not localize at any shift c = −8..8.
- **The factorization.** ⟨Ω|T_λ(θ₁−3h) = f_λ ⟨Ω_zw|T^{L=1,θ₁}_λ(θ₁−3h) ⊗ ⟨Ω_zw|, with:
  - f₍₁,₁₎ = (θ₁₂−2h)(θ₁₂−4h);
  - f₍₂,₂₎ = (θ₁₂−h)(θ₁₂−2h)(θ₁₂−3h)(θ₁₂−4h);
  - f₍₃,₂₎ = (θ₁₂−h)²(θ₁₂−2h)(θ₁₂−3h)(θ₁₂−4h).
- **Pure products.** No pure product of t's is local. T₍₂,₂₎(θ₁−3h) = t₂t₂ − t₁t₃ is a genuine combination, since t₃(θ₁−2h) ≠ 0.
- **Local products that add nothing.** T₍₂,₂₎·t₁(θ₁−h) and T₍₃,₂₎·t₂(θ₁+h) are local but do not enlarge the span. The site vector of the first is not of the form f·w + g·Ω.
- **Span.** rank{Ω_zw, v₍₁,₁₎, v₍₂,₂₎, v₍₃,₂₎} = 4.

**Site 2.** The same objects as for [2,1,0].
- **⟨Ω|t₂(θ₂+h)** = (θ₁₂−h)(θ₁₂−4h) · ⟨Ω_zw| ⊗ ⟨Ω_zw|t₂^{L=1,θ₂}(θ₂+h).
- **⟨Ω|t₂(θ₂+h)t₁(θ₂−h)** = (θ₁₂−h)(θ₁₂−4h) · ⟨Ω_zw| ⊗ [−(θ₁₂−2h)·⟨Ω_zw|t₂^{L=1}t₁^{L=1} + q₁(θ₂−h)·⟨Ω_zw|]. Here q₁(θ₂−h) = 6h³χ₃ is the one-site qdet.
- **[2,1,0] for comparison:** the prefactor is (θ₁₂−h)(θ₁₂−3h), and the bracket is −(θ₁₂−h)·(…) + q₁·Ω_zw.
- **Nothing else.** No other CBR shape, pair or triple localizes on site 2.
- **Span:** rank 3.

**Highest-weight reconstruction at small twist.**
- **Setup.**
  - Twists z = e^{εφ} and w = e^{εω}, with 20 values of w.
  - Data rows: the site-1 and site-2 covectors above, plus ⟨x₀| with ⟨x₀|Ψ⟩ = 1.
  - Decomposition: [3,1,0]² = [6,2,0] ⊕ [6,1,1] ⊕ [5,3,0] ⊕ 2[5,2,1] ⊕ [4,4,0] ⊕ 2[4,3,1] ⊕ [4,2,2] ⊕ [3,3,2], so there are 10 HW states.
- **Ranks.** site1+site2 is full rank in every dominant weight space **except (3,3,2): 16 of 18**. This is saturated: 30 values of w give the same ranks.
- **Recovered: 9 of 10 HW states, exactly** (relative error ≤ 10⁻⁸⁶ at ε = 10⁻⁶ and 10⁻⁹). This includes both [5,2,1] copies and both [4,3,1] copies, so the multiplicity mixing is fixed.
- **Not recovered: [3,3,2]** (relative error ≈ 1). Its HW direction lies in the null space of the localized data.
- **Adding the bare overlap row ⟨Ω(w)|Ψ⟩** makes all 10 exact, [3,3,2] included (10⁻¹⁰² at ε = 10⁻⁶, 10⁻⁴⁹ at ε = 10⁻⁹).
- **Contrast with [2,1,0].** There, site1+site2 alone gave everything, singlet scale included.
- **Which state fails.** [3,3,2] is the bottom irrep, with μ₃ = 2 > λ₂ = 1, the kind of state the old ν₃ > λ₂ rule flagged.

## Caveats
- The f's and the bracket coefficients were fitted at one (h, θ) point.
- The reconstruction uses M·Ψ for the Q side. This is an exact identity, but the chain was not rerun with Q functions.

## Traps
- **Nested pure functions.** Inside `Fold[... &]`, write `cbrGen[tApplyX[Λ, #1, #2, #3] &, ...]`. An inner `f[##]&` wrapper picks up the wrong slots.
- **Sort order.** `Sort` puts "t" before "T", so factor lists come out as {{"t",…},{"T",…}}.
- **ClearAll in cell 2.** It prints harmless `ClearAll::wrsym` messages for Wolfbook's protected symbols.

## Open
- **[3,3,2]: resolved (Paul, 2026-10-01).** ⟨Ω|Ψ⟩ itself is acceptable data. So site1 + site2 + ⟨Ω|Ψ⟩ + ⟨x₀|Ψ⟩ = 1 gives **all 10** HW CG states of [3,1,0]² exactly. All three ingredients are needed: dropping any one leaves (3,3,2) rank-deficient (16, 12 or 9 of 18).
- **Fits.** f_λ at a second (h, θ) point.
- **General pattern.** Not yet identified. For [2,1,0] the site-1 shapes are (1,1) and (2,1) at θ₁−2h; for [3,1,0] they are (1,1), (2,2) and (3,2) at θ₁−3h. Run the same notebook (it is generic in Λ) for [3,2,0] and [4,2,0].
