# gl(3) XXX, L=2, [4,2,0]⊗[4,2,0]: site-localizing transfer matrices and CG reconstruction

These notes cover `L2-420-localization-EH.wb`, from 2026-10-01.
- **What it is.** The [4,2,0] copy of the [3,2,0] notebook: dim 27, 729 two-site states.
- **Conventions and numeric point.** As in `L2-310-localization.md`.
- **Precision.** The confirmation step runs at 80 digits rather than exactly, because exact surd arithmetic through T₍₄,₂₎ on 27×27 covectors is too slow.

## Localizing objects (exact to 80 digits at the numeric point; f fitted at this one point)

**Site 1.** Everything is at **θ₁−4h**: **T₍₁,₁₎, T₍₂,₂₎, T₍₃,₂₎ and T₍₄,₂₎**.
- **The factorization.** ⟨Ω|T_λ(θ₁−4h) = f_λ ⟨Ω_zw|T^{L=1,θ₁}_λ ⊗ ⟨Ω_zw|, with:
  - f₍₁,₁₎ = (θ₁₂−2h)(θ₁₂−5h);
  - f₍₂,₂₎ = (θ₁₂−h)(θ₁₂−2h)(θ₁₂−4h)(θ₁₂−5h);
  - f₍₃,₂₎ = (θ₁₂−h)(θ₁₂−2h)²(θ₁₂−4h)(θ₁₂−5h);
  - f₍₄,₂₎ = (θ₁₂−h)²(θ₁₂−2h)²(θ₁₂−4h)(θ₁₂−5h).
- **Local products.** T₍₂,₂₎·t₁(θ₁−2h) and T₍₄,₂₎·t₂(θ₁+h) are also local but add no span.
- **Span:** rank 5.

**Site 2.**
- **t₂(θ₂+h)**, with f = (θ₁₂−2h)(θ₁₂−5h), the same as site-1 f₍₁,₁₎.
- **T₍₂,₂₎(θ₂)**, with f = the site-1 f₍₂,₂₎.
- **T₍₂,₂₎(θ₂)·t₁(θ₂−2h)** is local.
- **Span:** rank 4.
- **No pure pairs or triples** of t's are local on either site.

**Pattern across [2,1,0], [3,1,0], [3,2,0] and [4,2,0]** (observed, not derived):
- **Site 1** uses T₍₁,₁₎ and T₍ₖ, λ₁−λ₂₎ for k = λ₁−λ₂ … λ₁, all at θ₁ − λ₁h.
- **The T₍₁,₁₎ factor** is f₍₁,₁₎ = (θ₁₂ − (λ₁−λ₂)h)(θ₁₂ − (λ₁+1)h).
- **Site 2** always has t₂(θ₂+h), plus T₍₂,₂₎(θ₂) when λ₂ = 2.

## Reconstruction (ε = 10⁻⁶, 20 values of w; ranks saturated, 20 w = 30 w)

**Decomposition.** [4,2,0]² has 13 irreps and 19 HW states: [8,4,0], [8,3,1], [8,2,2], [7,5,0], 2[7,4,1], 2[7,3,2], [6,6,0], 2[6,5,1], **3[6,4,2]**, [6,3,3], [5,5,2], 2[5,4,3], [4,4,4].

**Results with site1 + site2 + Ω + x₀:**
- **16 of 19 exact,** including all three [6,4,2] copies.
- **Not recovered: [4,4,4] and both [5,4,3] copies.** The ranks are (4,4,4) 34/45 and (5,4,3) 32/36.

> **Correction (2026-10-02, `Blindspot-420-check-EH.wb`).** The "structural blind spot" below is wrong. Site-1 T₍₃,₂₎ and T₍₄,₂₎ at θ₁−4h have nonzero overlap with HW_[5,4,3], and T₍₄,₂₎ with HW_[4,4,4]. The overlaps are ≈8×10⁻³ and 8×10⁻⁴ at the finite test twist and scale like ε³ and ε⁴ at small twist. The "1.000" singular values were printed at 3 digits and only show near-containment in the numerical null space, not an exact zero. The reconstruction failure at ε = 10⁻⁶ is real; its cause is open.
>
> **Follow-up (2026-10-02, `Minimal-ingredients-420-EH.wb`, 250 digits).** The seven localized insertions plus ⟨Ω|Ψ⟩ determine **all 19** HW CG states, using clean leading-order combinations (top-killed, then cancelling the structural contamination). That includes [5,4,3]×2 and [4,4,4]. Errors are O(ε) in a fair test: fitted at ε = 10⁻⁶ and applied at a different ε and φ. [6,4,2] (multiplicity 3) needs per-copy clean combinations. The coefficients are fitted, not yet in closed form.

**The blind spot is structural (superseded, see the correction above).**
- The null space of the localized data contains HW_[5,4,3] and HW_[4,4,4] **entirely**: the overlap singular values are exactly 1. So the localized data have zero overlap with those HW vectors at every ε.
- No ε-order, choice of w, or choice of φ can recover them.
- For [4,4,4] (multiplicity 1) only the normalization is lost, since the direction is fixed by the HW condition. For [5,4,3] (multiplicity 2) the mixing itself is lost.

**Non-localized t-rows close it.**
- **Which rows.** ⟨Ω(w)|t_a(u)| at two generic points, plus one product. These are Q-computable, since ⟨Ω|f(t)|Ψ_B⟩ = f(Λ_B)⟨Ω|Ψ_B⟩.
- **Ranks.** With them, every sector is full rank.
- **Reconstruction.** **All 19 HW states are exact** (≤ 4×10⁻³⁷).

## φ-dependence (see `L2-320-localization-EH.wb` Section 7)
- **The limit basis is canonical.** The ε → 0 limit basis in a multiplicity space does not depend on the ket twist direction φ. It is the eigenbasis of the **untwisted** inhomogeneous transfer matrices on HW_μ, because θ₁ ≠ θ₂ already splits the copies.
- **So varying φ is a consistency check.** It gives no new target data, and it cannot repair a blind spot.
