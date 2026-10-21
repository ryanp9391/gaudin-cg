# [3,1,0] point A through the untwist pipeline: all ten HW CG data of [3,1,0]⊗[3,1,0] from the det₆ × det₂^{(0,1,2)} formula (2026-10-10)

gl(3) XXX, L = 2, Λ = (3,1,0) on both sites, new conventions, diagonal frame, plain GT generators (exact sparse). Notebook: `L2-310-untwist-CG.wb` (90 cells; **no native outputs** — VS Code never opened the file, the sonnet agent ran each cell's source in K1 and stored the printed results as markdown cells under each code cell). Spec: `L2-310-untwist-CG-spec.md` (round 1) and the round-2 task message recorded in HANDOFF. The adjoint version of this pipeline is `L2-210-untwist-CG.md`; the point-A data come from `Undo-site2-rotation-notes.md` §5.2–5.4, §6.

**Evidence tags.** [exact], [numeric], [fitted] as in the adjoint write-up. Chain points: P1 (h = 1/2, θ = (2/7,−3/5), t = 62/35), P3 (h = 1/2, θ = (3/7,−1/4), t = 19/14, in no fit set); scans at h = 1, θ = (t,0), 20 values of t (14 fit, 6 held out). Directions and ω grid as for the adjoint; ε ∈ {1,2,4(,8)}·10⁻⁶, 200 digits. One t point ≈ 30 s.

---

## 0. Summary

1. **Regression.** Point A rebuilt from the closed-form dressing: the undo is exact on the 10 HW kets and the identity ⟨v₁⊗v₂|Ψ⟩/⟨x₀|Ψ⟩ ∝ det₆·det₂^{(0,1,2)}[y_A]/(Ψ_Q W(3,4)) holds on all 225 states to 1.3·10⁻¹⁷⁹ (controls y = 2, 3 and pole set {0..3} fail by 10–45). K13 = −189/704 (the notes' +189/704 is a sign typo). [exact/numeric]
2. **Regime (same as the adjoint).** At ε = 0 the dressing cancels E₁₂² exactly: φ_A(0)/δ^Λ = e₁ᵀ(1 + (a′₀/2)E₁₂)e^{c′₀E₂₃}, γ = O(ε), b′/(a′c′) = 1 + ε(t−2)(φ₂−ω₂). Hence **three orders** are needed (two for the adjoint): weights at ε⁰: (3,1,0), (3,0,1), (2,2,0), (2,0,2), one direction of (2,1,1); at ε¹: (1,3,0), one direction of (1,2,1), (1,1,2); at ε²: (0,2,2), (1,0,3), (0,3,1), the second directions of (1,2,1), (1,1,2); at ε³: (0,1,3). [exact]
3. **Clean at leading order** (κ−1 ∝ ε): [6,1,1], [5,3,0], [5,2,1]×2 (u′-functional). **Structural O(ε) constants** (ε-, ω-, φ-independent, closed forms in t by the t-scan; τ_c = the copy's untwisted t₂(θ₁+h)/h⁴ eigenvalue, a Q-number): [4,4,0] κ = (t+1)/2; [4,2,2] κ = (11t−31)/(5(3t−7)); [3,3,2] κ = (t+1)(19t−49)/(40(t−2)); [4,3,1] κ_c = A + Bτ_c with D = t³+3t²−17t+5, A = (−7t⁴+26t³−166t+75)/(9D), B = (7t²−22t+23)/(36D). [fitted, held-out ≤ 10⁻¹⁵]
4. **Second directions.** [5,2,1] pair: R_c = 3/2 − 2t + τ_c/(2(t+4)) [exact to 10⁻¹⁸ on held-out], from the ω₂-slope with signal slopes (t−5)/4, (7−3t)/2 and a remainder affine in τ_c (noise-limited 10⁻⁹ fit). ε² layer: [4,2,2] R = −2(t−3)/(3(3t−7)) [exact], its structural slope constant a rational fit; [4,3,1] R_c = A_R + B_Rτ_c [exact], its structural slope constant only partly identified (β(t) exact, α(t) a degree-8 numerical fit, 10⁻⁹).
5. **FINAL.** All ten branches extracted from Qnum + the top-kill det₆ + the Q-dictionary τ_c + the closed forms, at P1 and P3: errors ≤ 4.7·10⁻⁶ at ε = 10⁻⁶ (O(ε)), Richardson agreement ≤ 1.6·10⁻¹⁸; second directions to 10⁻¹⁵ except [4,3,1] (10⁻⁸, limited by α(t)). [numeric]

---

## 1. Setup and regression (cells 1–7)

Sparse exact machinery of `L2-310-two-rotations.wb` (7 minor families + one twist: 1.0 s), invariant metric G for (3,1,0), G-orthogonal HW projectors from the exact HW kets, factorized x₀, point-A covectors from the closed forms β, γ (Undo notes §5.3), `reg` as in the adjoint, fixed-degree Q solvers (dual q^j: μ_j; vacuum direct: 0, 4, 6), 8-pole det₆, det₂^{(0,1,2)}, W(3,4), Ψ_Q with the vacuum triple (4,3,1). All 225 states at 200 digits: 13.4 s; det₆×det₂ on 225 states: 3.7 s; the 56 states of the 8 HW sectors: 3.1 s per ε. Baxter residuals: dual 3.8·10⁻¹⁸⁴ (regression) / 2.4·10⁻¹⁴⁰ (small ε), vacuum 5.3·10⁻¹⁹⁵.

## 2. ε-structure of φ_A (Task 0, cells 10–11) [exact]

γ(ε) = ε·(−2565/56056) + … at (φ, ω_A); a′ exactly ε-free; b′₀ = a′₀c′₀. Ranks over the 3×6 grid identical for φ and φ′:

| site-1 weight | leading order | rank (ε⁰ / incl. ε¹ / incl. ε²) |
|---|---|---|
| (3,1,0), (3,0,1), (2,2,0), (2,0,2) | ε⁰ | 1 |
| (2,1,1) dim 2 | ε⁰ (u′ = e₁ᵀE₁₂E₂₃ = (2,8)); second at ε¹ (u = e₁ᵀE₁₃ = (−1,8)) | 1 / 2 |
| (1,3,0) | ε¹ | 1 |
| (1,2,1), (1,1,2) dim 2 | ε¹ one direction; second at ε² | 1 / 1 / 2 |
| (0,2,2), (1,0,3), (0,3,1) | ε² | 1 |
| (0,1,3) | ε³ | 0 through ε² |

## 3. Chain side (Task 1, cells 20–22) [numeric]

Undo exact to 4.9·10⁻¹³⁸ at 540 points. ε-exponents of (S, K, ‖Ψ‖): [6,2,0] (−6, 0, −2); [6,1,1], [5,3,0] (−7,−6,−3); [5,2,1]×2 (−8,−7,−4); [4,4,0] (−7,−7,−4); [4,3,1]×2 (−8,−8,−5); [4,2,2] (−9,−9,−6); [3,3,2] (−9,−9,−7). κ at P1: clean branches 1 + O(10⁻⁶); [4,4,0] 1.3857144 (= 97/70), [4,3,1] 1.4105702 / 0.986713, [4,2,2] 1.3661018 (= 1 + 108/295), [3,3,2] 4.6508004; all constant in ε, ω, φ (spreads ≤ 10⁻⁵). Contamination carriers (undone frame): [4,4,0] ← (2,2,0); [4,3,1] ← (2,1,1); [4,2,2] ← (2,0,2); [3,3,2] ← (1,1,2) 85% + (2,0,2) 15%.

## 4. Q side (Task 2, cells 23–24) [numeric]

R_C/R_Q state-independent to 2.5·10⁻⁸², 3.7·10⁻⁸⁶, 8.9·10⁻⁹³ (ε = 1, 2, 4·10⁻⁶) over the 10 branches, 6.4·10⁻⁷⁷ over all 56 HW-sector states (limited by the 104–112-digit precision of R_Q). Leading-order extraction (fac = c′₀, a′₀ [the GT element of e₁ᵀE₁₂ is 2], a′₀c′₀): spread/error 0.85–4.2·10⁻⁶ at ε = 10⁻⁶, ×4 at 4·10⁻⁶; Richardson limits Q = chain to ≤ 1.2·10⁻¹⁹ ([6,1,1] 46.6745…, [5,3,0] −36.8838…, [5,2,1] −1325.40…, −334.33…). The [6,2,0] and [6,1,1] limits coincide with the adjoint's 27 and [4,1,1] values.

## 5. The O(ε) and O(ε²) layers (Tasks 3, A, B; cells 30–34 and round 2)

**Top-kill** T = δ^Λ s⟨v₂⊗v₂|Ψ⟩/⟨x₀|Ψ⟩ equals the (w₁,y_A)-vacuum det₆/Ψ_Q (normalized to the top) to 2.5·10⁻⁸²; T ~ ε⁻⁶ for every branch, relevant only for [6,1,1], [5,3,0] (T/K = 2.30, −1.24).

**Exact untwisted t₂(θ₁+h) on the multiplicity spaces** (h = 1; these feed τ_c): HW_[4,3,1]: trace 8t²+16t−8, det 16t(t−1)(t+1)(t+4); HW_[5,2,1]: trace 8t²+19t−8, det 4(t+4)(t−1)(t+1)(4t+3). t₁(θ₁+h) is degenerate on the [4,3,1] pair.

**O(ε) constants** (t-scan; all held-out residuals ≤ 10⁻¹⁵ unless stated):
- [4,4,0]: κ = (t+1)/2. [4,2,2]: κ = (11t−31)/(5(3t−7)) (pole at t = 7/3, where the chain value is −2.5·10⁶). [3,3,2]: κ = (t+1)(19t−49)/(40(t−2)). [4,3,1]: κ_c = A + Bτ_c, D = t³+3t²−17t+5, A = (−7t⁴+26t³−166t+75)/(9D), B = (7t²−22t+23)/(36D).
- [5,2,1] pair, weight (2,1,1), u′ = (2,8), u = (−1,8), R_c := (u·C)/(u′·C): F₁/F₀ a straight line in ω₂ (≤ 1.3·10⁻⁹), ω₁- and φ-independent slope; exact signal slopes p₁′ = (t−5)/4 (u′-part) and q₁′ = (7−3t)/2 (u-part, times R); remainder slope Δ_c = slope − p₁′ − q₁′R_c = (111+44t−92t²)/105 + (23t−13)τ_c/(105(t+4)) (10⁻⁹, noise-limited fit; the 1/105 coefficients are rationalized, not proven). **R_c = 3/2 − 2t + τ_c/(2(t+4))** (10⁻¹⁸).

**ε² layer** (four-ε Richardson; convention: u′ = (1, x_U) first component normalized, u = (0,1), R = C₂/(C₁ + x_U C₂) on the two GT indices; x_U = −(t+3)/(2(t−3)) at (1,2,1), −3(t−1)/(2(t−3)) at (1,1,2)): the ε² covector coefficients a₂u′ + b₂u have a₂/fac₁, b₂/fac₁ exactly linear in ω₂ and ω₁-independent ((1,2,1): a-slope t(2−t)/(t−3), b-slope 525/3698 at P1; (1,1,2): a-slope −19/35 at P1, b-slope 3(t−1)/(t−3)); F₁/F₀ a straight line in ω₂ (≤ 4·10⁻¹⁵), slope ω₁- and φ-independent; the structural constant kslope = κ₀·(remainder slope):
- [4,2,2]: R = −2(t−3)/(3(3t−7)) [exact]; kslope = (−(28/45)t³ + (2912/675)t² − (448/45)t + 4948/675)/(t−7/3)² (10⁻¹¹ held-out, noise-limited).
- [4,3,1]: R_c = A_R + B_Rτ_c with A_R = (−10 + 28t/3 + 4t² − 4t³ + 2t⁴/3)/D, B_R = (−3/2 + t − t²/6)/D [exact]; kslope_c = α(t) + β(t)τ_c with β = (3331 − 7425t + 7876t² − 2262t³ − 481t⁴ + 247t⁵ − 6t⁶)/(216D²) (10⁻¹¹) and α = P₈(t)/D² a numerical degree-8 fit (10⁻⁹) — **not identified**. The remainder slope itself is not affine in τ_c (because κ_c is).

## 6. Final pure-Q extraction (Task C) [numeric]

Inputs: Qnum, the top-kill det₆, τ_c from t₂(θ₁+h) = χ₃ν₁(θ₁+h)ν₂(θ₁)Q^j(θ₁−h)/Q^j(θ₁) (j-independent), the closed forms above. All 10 branches at P1 and P3: errors ≤ 4.7·10⁻⁶ at ε = 10⁻⁶, Richardson agreement ≤ 1.6·10⁻¹⁸ (P1 and P3). Second directions (|ΔR| / slice error): [5,2,1] ≤ 2.6·10⁻¹⁶ / 1.3·10⁻¹⁵ (P1), ≤ 4.7·10⁻¹⁷ / 2.5·10⁻¹⁶ (P3); [4,2,2] 1.9·10⁻¹⁶ / 3.8·10⁻¹⁶, 1.8·10⁻¹⁷ / 5.0·10⁻¹⁷; [4,3,1] ≤ 2.4·10⁻¹⁰ / 7.4·10⁻⁹, ≤ 2.3·10⁻⁸ / 3.3·10⁻⁷ (limited by α(t)). P1 τ_Q: 13.2820, 37.4788 ([5,2,1]); 9.81505, 35.6315 ([4,3,1]).

## 7. Gotchas

- `Do[DownValues[f] = …]` purges nothing; and a memoized quantity (`qtkC`) that depends on (h, θ) but was keyed only on (φ, ε, ω) survived `setPoint` and made a first ε²-scan t-independent. Purge every (h,θ)-dependent memo with literal symbols.
- Poles of the structural constants at t = 7/3 ([4,2,2]) and t = 2 ([3,3,2]); exclude such t from fits.
- t₁(θ₁+h) does not split the [4,3,1] pair; use t₂(θ₁+h).

## 8. Open

1. Identify α(t) for the [4,3,1] ε² constant (degree-8 numerical fit; needs either higher accuracy — four-ε Richardson at smaller ε — or the perturbative derivation: second-order PT as in `Adjoint-untwist-perturbative.wb`).
2. Prove the fitted forms (κ's, Δ_c, kslope) by first/second-order perturbation theory for [3,1,0].
3. The [5,2,1] remainder and the [4,2,2] kslope are rationalized from 10⁻⁹–10⁻¹¹ fits.
4. [4,1,0] k = 2 (dressing of E₁₃-degree 2: expect four orders).
5. Store outputs (rerun in VS Code).
