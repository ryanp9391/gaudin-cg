# gl(3) XXX, L=2, [3,2,0]⊗[3,2,0]: site-localizing transfer matrices and CG reconstruction

These notes cover `L2-320-localization-EH.wb`, from 2026-10-01.
- **What it is.** The [3,2,0] copy of `L2-310-localization-EH.wb`. The method, conventions and numeric point are the same; see `L2-310-localization.md`.
- **Running it.** The notebook is self-contained, clears the kernel, and reruns in about 40 s. [2,1,0] is the control.
- **Fits.** One-site fits use all 15 components. A 3-component fit is degenerate for T₍₁,₁₎ here.

## Results (exact at h = 1/2, θ = (2/7,−3/5), z = (3/2,5/4,2/3), w = (5/2,7/3,3/5); coefficients fitted at this one point)

**Site 1.** Everything is at **θ₁−3h**: **T₍₁,₁₎, T₍₂,₁₎ and T₍₃,₁₎**, which matches Paul's guess. For [3,1,0] the shapes were T₍₁,₁₎, T₍₂,₂₎ and T₍₃,₂₎.
- **The factorization.** ⟨Ω|T_λ(θ₁−3h) = f_λ ⟨Ω_zw|T^{L=1,θ₁}_λ ⊗ ⟨Ω_zw|, with:
  - f₍₁,₁₎ = (θ₁₂−h)(θ₁₂−4h);
  - f₍₂,₁₎ = (θ₁₂−h)(θ₁₂−2h)(θ₁₂−4h);
  - f₍₃,₁₎ = (θ₁₂−h)²(θ₁₂−2h)(θ₁₂−4h).
- **Local products.**
  - T₍₁,₁₎·t₁(θ₁−2h) has f = (θ₁₂−h)(θ₁₂−2h)(θ₁₂−4h) and g = f₍₁,₁₎·q₁(θ₁−h), where q₁ = 6h³χ₃ is the one-site qdet.
  - T₍₃,₁₎·t₂(θ₁+h) is also local.
  - Neither adds to the span.
- **Span:** rank 4.

**Site 2.**
- **t₂(θ₂+h)**, with f = (θ₁₂−2h)(θ₁₂−4h).
- **T₍₂,₂₎(θ₂)**, new for this rep, with f = (θ₁₂−h)(θ₁₂−2h)(θ₁₂−3h)(θ₁₂−4h).
- **T₍₂,₂₎(θ₂)·t₁(θ₂−2h)** is local, but its site vector is not of the form f·w + g·Ω.
- **t₂(θ₂+h)t₁(θ₂−h) is NOT local here**, although it was for [2,1,0] and [3,1,0].
- **Span:** rank 4.

**Duality with [3,1,0].** The roles of the two sites swap.
- [3,2,0] site 2 has the same f's as [3,1,0] site 1: T₍₁,₁₎-type and T₍₂,₂₎-type.
- [3,2,0] site 1 matches [3,1,0] site 2: t₂-type with (θ₁₂−h)(θ₁₂−4h), and the product with the same f and g up to sign.
- This is consistent with [3,2,0] ≅ [3,1,0]* ⊗ det³.

**Highest-weight reconstruction at small twist** (ε = 10⁻⁶ and 10⁻⁹, 20 values of w, ⟨x₀|Ψ⟩ = 1).
- **Decomposition.** [3,2,0]² = [6,4,0] ⊕ [6,3,1] ⊕ [6,2,2] ⊕ [5,5,0] ⊕ 2[5,4,1] ⊕ 2[5,3,2] ⊕ [4,4,2] ⊕ [4,3,3], so there are 10 HW states.
- **Locality at other twists.** Exact at ε = 1/3 for 20 values of w.
- **site1 + site2:** 9 of 10 exact (≤ 4×10⁻⁸⁸). This includes both multiplicity-2 pairs.
- **The exception is the bottom irrep [4,3,3]**: its weight space reaches rank 16 of 18.
- **Adding ⟨Ω(w)|Ψ⟩** (allowed data, per Paul) makes **all 10 exact**. The worst is 1.7×10⁻⁴⁸ at ε = 10⁻⁹ and 4.7×10⁻¹⁰² at ε = 10⁻⁶.

## How many ε-orders are needed (Section 6)

**Setup.**
- **Data matrix.** M(ε) has rows site1, site2 and Ω for 8 values of w, plus ε⁻⁴x₀ (‖x₀‖ ~ ε⁴). Its Taylor coefficients M_j come from a Cauchy contour, |ε| = 1/8 with 40 points; resummed, they agree to 3.5×10⁻¹⁸.
- **State.** Ψ_B = ε⁻ᵏ Σ ψ_j εʲ, so the data through relative order N give the block-Toeplitz system Σ_{j≤n} M_j ψ_{n−j} = b_n for n ≤ N.
- **N_min** is the smallest N at which the null space has no ψ₀ component.

| μ | mult | k (Ψ ~ ε⁻ᵏ) | N_min |
|---|---|---|---|
| [6,4,0] | 1 | 4 | 0 |
| [6,3,1], [5,5,0] | 1 | 5 | 1 |
| [5,4,1] | 2 | 6 | 1 |
| [6,2,2] | 1 | 6 | 2 |
| [5,3,2] | 2 | 7 | 2 |
| [4,4,2] | 1 | 8 | 2 |
| [4,3,3] | 1 | 9 | **5** |

- **Reading the table.** Each data function ⟨c|Ψ_B⟩ is needed through ε^{N_min−k}.
- **The HW prior** on ψ₀ makes no difference.
- **The null spaces are exact,** so the cut is clean.
- **Trap.** x0P's generic `Inverse` of P_phys loses all precision at complex ε. Use the factorized inverse `x0S` (U⁻¹D⁻¹L⁻¹).

## φ-dependence (Section 7)
- **Four values of φ give the same limit.** The multiplicity-2 limit vectors in [5,4,1] and [5,3,2] are identical for φ with s = (φ₁−φ₂)/(φ₂−φ₃) = −143/171, 8/57 and 1/2, and for 2φ+(1,1,1). Agreement is 10⁻⁹ = O(ε).
- **They are the untwisted eigenbasis.** The vectors are the eigenvectors of the **untwisted** inhomogeneous transfer matrices restricted to HW_μ, which is an invariant subspace.
- **So the CG basis is canonical,** independent of φ and ω.
- **Consequence for varying φ:** it is a consistency check, not new data. Pooling over φ cannot fix what one φ leaves undetermined, because the per-φ normalizations are unrelated. *(This is an argument, not tested.)*

## Caveats
- One (h, θ) point and one small-twist direction φ were used.
- The reconstruction uses M·Ψ for the Q side. This is an exact identity, but the chain was not rerun with Q functions.
