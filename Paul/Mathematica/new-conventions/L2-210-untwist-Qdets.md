# Untwisting the Q-determinants themselves: ε-structure of det₂^{(1)}/minor₁ and det₆/Ψ_Q in the joint regime (adjoint, 2026-10-10)

Notebook: `L2-210-untwist-Qdets.wb` (45 cells; no native outputs — built from cell sources, run in K1, printed results stored as markdown cells). Spec: `L2-210-untwist-Qdets-spec.md`. Setting as in `L2-210-untwist-CG.md`: P1, joint regime z = 1+εφ, w₁ = 1+εω₁, w₂ = 1+εω₂, y = w₂y₀; ε ∈ {1,2,4,8}·10⁻⁶; 48 direction points (φ, φ′ × 18 ω; two more φ's × 6 ω); 200 digits. Branch order: 27, [4,1,1], [3,3,0], octet 1, octet 2, singlet.

## 0. Summary

1. **det₂^{(1)}/minor₁ is fully untwisted.** [numeric, 10⁻¹⁴⁰ / 10⁻¹⁹⁷] It equals c·F(y)·(1, τ₁₀, τ₁₁, τ₂₀, τ₂₁) with c = 343000/1161 at small ε for every branch and direction; τ_λ(ε) = τ⁰_λ + ετ¹_λ + O(ε²) with τ⁰ the untwisted eigenvalue and τ¹ = Σ_k φ_k a_k the first-order shift (linear in φ, Σ_k a_k = |λ|τ⁰); and τ⁰_λ(x) = P_λ(x)·det[v(x+λ₁h), v(x+(λ₂−1)h), v(x−2h)]/det[v(x), v(x−h), v(x−2h)] with v any basis of the 3-dim polynomial kernel of the untwisted dual Baxter equation (staggered degrees (μ₃, μ₂+1, μ₁+2)), to 10⁻¹⁹⁷ at four values of x. So the leading order is an explicit untwisted-Q formula and the next order needs only the first-order eigenvalue shifts.
2. **det₆/Ψ_Q is not.** [numeric] Its ε-exponents are 10,10,10,9,9,8; the ratios D6_B/D6_top are exact monomials in the direction differences (constants 1, 1, algebraic, −9800/2619) and for [4,1,1], [3,3,0] even exact in ε: D6_[3,3,0]/D6_top = (w₁−z₁)/(w₁−z₂), D6_[4,1,1]/D6_top = (w₁−z₂)(w₂−z₂)/((w₁−z₃)(w₂−z₃)) (deviations ≤ 6·10⁻¹²²), with no Q input. But the old thread's Laplace split (bare block × periodic block) holds only for the top state: for the other five branches the mixed Laplace terms are of the same order, and the untwisted "periodic minor" vanishes identically on non-top states (orthogonality at zero twist). The H column combinations are rational and branch-independent (H1 = c1+c5, H2 = c3−c5, H3 = c4−(5/2)c5−c6 at h = 1/2 — the 5/2 is presumably 5h, untested); the O complement that would make the split work is not rational.
3. **Qnum to O(ε)** is reconstructed to ε² accuracy (10⁻¹² at ε = 3·10⁻⁶, held out) from: the untwisted Wronskian τ⁰, τ¹ (untwisted chain), the exact top-state data, F(y), and the D6 ratio = monomial × (1 + ερ_B) with ρ_B numerical for the octets and singlet (it vanishes for [4,1,1], [3,3,0]).

## 1. Exponent and monomial tables (Task 1)

| factor | ε-exponent per branch | direction dependence |
|---|---|---|
| D2 = det₂/minor₁ | 1, 0, 0, 0, 0, 0 | top: (893025/1024)(ω₂−φ₂)ε; non-top: constants (−1157625/512 for the decuplets, −1157625/128 singlet, algebraic for octets) |
| D6 = det₆/Ψ_Q | 10, 10, 10, 9, 9, 8 | not a monomial (common branch-independent factor); ratios to the top are monomials (§2) |
| D6t (top-kill det₆ at (w₁,y)) | 6 | pure monomial in ω₁: Π_k(ω₁−φ_k)^{4−μ_k} |
| ⟨x₀|Ψ̂⟩ (‖Ψ‖ ~ ε⁻ᵖ) | 2, 3, 3, 4, 4, 6 | |φ₁₂|², |φ₁₂φ₁₃φ₂₃|, |φ₁₂|³, φ₁₂²|φ₁₃φ₂₃| (×2), |φ₁₂φ₁₃φ₂₃|² |
| Qnum | −5, −6, −6, −7, −7, −8 | exact monomials: common φ₁₂⁻²(ω₂−φ₁)⁻²(ω₂−φ₂)⁻¹; [4,1,1]: −70/27·(ω₁−φ₂)/((ω₁−φ₃)(ω₂−φ₃)); [3,3,0]: −70/27·(ω₁−φ₁)/(ω₁−φ₂)·(ω₂−φ₁)⁻²(ω₂−φ₂)⁻²; octets: 6.3383, −8.2030 ·(ω₁−φ₁)/(ω₁−φ₃); singlet: 38.8047·(ω₁−φ₁)²/(ω₁−φ₃)²·(ω₂−φ_k)⁻² for k = 1,2,3 |

## 2. det₆ structure (Task 3)

- Entry exponents p_ab = r_a + s_b, s = (0,−1,0,−1,0,−1); row exponents (bare | periodic): 27 (5,3,3 | 0,−2,−2), [4,1,1] (4,4,2 | 0,0,−2), [3,3,0] (5,2,3 | 0,−3,−2), octets (4,3,2 | 0,−1,−2), singlet (3,3,1 | 0,0,−2).
- Split det₆ = −BO·PH(1+O(ε)): exponents add (12,12,12,11,11,10 = BO + PH), but the ratio tends to 1 only for the top; −66/53, 0.33, 0.72, 0.72, 0.39 (direction dependent) otherwise.
- BO_B/BO_top exact monomials with constants −5/3, −1, 5/2, 5/2, 5/2 (no relation to the old E_λ); PH_B/PH_top not a monomial; PH does not tend to the untwisted minor (which is 13453 for the top, 0 for the rest). The vacuum direct kernel at w = 1 has staggered degrees (0,3,6) and coincides with the top state's.
- D6 octet/singlet ratios carry an O(ε) correction ρ_B that is not a short rational function of the directions.

## 3. Open

1. The det₆ factor: a column basis (possibly state/direction dependent, or with O(ε) admixtures) that makes a bare×periodic term dominant for non-top branches; or a different route (e.g. through the known exact ratios for [4,1,1], [3,3,0], which are pure Gauss-parameter functions and suggest an algebraic rather than SoV-determinant origin).
2. The common non-monomial factor of D6/BO/PH (state-independent, w₃-independent).
3. τ¹ in Q-language (currently from the untwisted chain: ⟨h|T_λ′|h⟩).
