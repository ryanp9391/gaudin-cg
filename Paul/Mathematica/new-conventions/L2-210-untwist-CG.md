# Untwisting the two-rotation det₆ × det₂ formula and reading off the CG coefficients of adj ⊗ adj (2026-10-10)

gl(3) XXX, L = 2, Λ = (2,1,0) on both sites, new conventions, diagonal frame, plain GT generators `JJ` (exact rationals). Notebook: `L2-210-untwist-CG.wb` (45 cells, saved; **no stored outputs**, see Tooling). Spec: `L2-210-untwist-CG-spec.md`. Builds on `Two-rotations-det-overlap.md` §5.6b/§5.8 (the family and its det formula), `Undo-site2-rotation-notes.md` §1–§4 (the undo and the per-weight counts) and `Clean-CG-combinations.md` §2–3 (signal/contamination). Implementation by a sonnet wolfbook-builder from the spec; planned and reviewed by the main session.

**Evidence tags.** [exact]: symbolic identity or exact-rational computation; [numeric]: at the stated points with the stated residuals; [fitted]: closed form recognized from a scan and validated on held-out points.

**Chain points.** P1: h = 1/2, θ = (2/7, −3/5), t := θ₁₂/h = 62/35. P2: h = 1/3, θ = (1/5, −4/9), t = 29/15. P3: h = 1/2, θ = (3/7, −1/4), t = 19/14 (in no fit set). Scans in t: h = 1, θ = (t, 0). Directions φ = (1/3, −2/7, 5/11) (control φ′ = (−2/5, 3/8, 7/11)); ω₁ ∈ {3/4, −5/3, 7/2}; ω₂ ∈ {−1/5, 9/7, −12/5, 5/3, 2/9, −7/4}. ε ∈ {10⁻⁶, 2·10⁻⁶, 4·10⁻⁶}, 200 digits.

---

## 0. Summary

1. **Regime.** Untwist everything together: z = 1+εφ, w₁ = 1+εω₁, w₂ = 1+εω₂ and y = w₂·y₀ with y₀ = (θ₁₂−h)/(θ₁₂−2h), so that the product condition w₂(y) = w₂ holds exactly. Then the undone covector φ = v₁ρ(n)⁻¹ has O(1) components at the site-1 weights (2,1,0), (1,2,0), (2,0,1), (1,0,2) and one direction of (1,1,1), and the signal/contamination ratio K/S is O(ε) for [4,1,1], [3,3,0] and both octets. **[exact for the covector, numeric for K/S]**
2. **The E₁₃ direction degenerates by one order:** b′/(a′c′) = 1 + ε(t−2)(φ₂−ω₂) + O(ε²) exactly, so the second octet direction, (0,2,1) and the singlet direction (0,1,2) are O(ε), at the same order as the first contamination. **[exact]**
3. **That contamination is structural.** Singlet: K/S = −t/(3(t−2)), i.e. κ_singlet = 2(t−3)/(3(t−2)), independent of h, φ, ω. Octets: the ω₂-slope of the O(ε) top-killed Q-number has contamination Δ_c = (1−t²) + t·τ_c/(3(t+3)), with τ_c the copy's untwisted eigenvalue of t₂(θ₁+h)/h⁴, a Q-number through the dual-Q dictionary. **[fitted; validated at P2, P3]**
4. **All six HW CG data of adj⊗adj are recovered from the det₆ × det₂ formula alone**, with O(ε) errors at finite ε (≤ 4·10⁻⁶ at ε = 10⁻⁶) and Richardson limits agreeing with the chain side to ≤ 6·10⁻¹⁹, at P1, P2 and P3. Octet mixing R_c agrees to 10⁻¹¹. **[numeric]**
5. **Bug found and fixed:** `Do[DownValues[f] = …, {f, {…}}]` purges nothing (`DownValues` is HoldAll); stale P1 transfer matrices produced a fake singlet law −(62/105)/(t−2). Purge with literal symbols (cell 36).

---

## 1. The regime and the ε = 0 covector (Task 0, cells 10–12) [exact]

With δ_k = 1/(w₂−z_k) = ε⁻¹/(ω₂−φ_k), δ^Λ = δ₁²δ₂ and the Gauss parameters r_k of Ω_{(w₁,y)}:
- a′ = r₁(δ₂−δ₁)/δ₁ = (ω₁−φ₁)(φ₂−φ₁)/((ω₁−φ₂)(ω₂−φ₂)), **exactly ε-free**;
- c′ = r₃(ε)(φ₃−φ₂)/(ω₂−φ₃) = O(1), with r₃'s O(ε) term ω₂-independent;
- b′/(a′c′) = (y−w₂)/(y−z₂) = 1 + ε(t−2)(φ₂−ω₂) + O(ε²) (24/1225 at P1, φ, ω₂ = −1/5);
- φ(0)/δ^Λ = e₁ᵀ e^{a′₀E₁₂} e^{c′₀E₂₃} (b′₀ = a′₀c′₀; the reversed order fails).

Leading ε-orders per site-1 weight and ranks over the 3×6 (ω₁, ω₂) grid (ε⁰ coefficient, ε¹ coefficient, both):

| weight | dim | leading order | ranks (ε⁰, ε¹, both) |
|---|---|---|---|
| (2,1,0) | 1 | ε⁰ | (1,0,1) |
| (1,2,0) | 1 | ε⁰ | (1,0,1) |
| (2,0,1) | 1 | ε⁰ | (1,1,1) |
| (1,1,1) | 2 | ε⁰ | (1,2,2) |
| (1,0,2) | 1 | ε⁰ | (1,1,1) |
| (0,2,1) | 1 | ε¹ | (0,1,1) |
| (0,1,2) | 1 | ε¹ | (0,1,1) |

The (0,1,2) component is 3b′(a′c′−b′)δ^Λ = −3(t−2)(φ₂−ω₂)(a′₀c′₀)²δ^Λ·ε + O(ε²). The (1,1,1) component is b′·u + a′c′·u′ with u = e₁ᵀE₁₃ = (−1,3) and u′ = e₁ᵀE₂₃E₁₂ = (2,0) on the two (1,1,1) GT indices; at ε = 0 only u+u′ = (1,3) is seen.

Why this regime: at fixed y the components are O(ε^depth) (the old table ε⁰,ε¹,ε¹,ε²,ε³,ε³,ε⁴) and the contamination is at or above the signal everywhere. Here the resolvent torus Δ = ε⁻¹diag(1/(ω₂−φ_k)) is a genuine torus element at ε = 0, so the two rotations stay distinct. The E₁₃ degeneracy is the ε = 0 shadow of the old b′ = 0 blind spot: the vacuum's r₂/(r₁r₃) and Δ are both built from φ, and n₀Δ n₀⁻¹ then has b′ = a′c′.

## 2. Chain side: signal vs contamination (Task 1, cells 13–16) [numeric]

Ψ_B normalized by ⟨x₀|Ψ⟩ = 1, P_HW the projector orthogonal for the invariant metric G = diag(1,1,1,2,6,9,4,9) (the plain `JJ` are not Hermitian), S = ⟨v₁⊗v₂|P_HWΨ⟩, K = R_C − S; the exact undo S = s⟨φ⊗e₁ᵀ|P_HWΨ⟩ holds to 10⁻¹⁵² at all 324 points.

| branch | ε-exponent of S, K, ‖Ψ‖ | max |K/S| at ε = 10⁻⁶, 2·10⁻⁶, 4·10⁻⁶ |
|---|---|---|
| [4,2,0] | −5, K ≡ 0, −2 | 0 |
| [4,1,1] | −6, −5, −3 | 9.6e-7, 1.9e-6, 3.8e-6 |
| [3,3,0] | −6, −5, −3 | 9.4e-7, 1.9e-6, 3.7e-6 |
| octet 1, 2 | −7, −6, −4 | 2.8e-6 / 1.1e-6 (×2, ×4) |
| singlet | −8, −8, −6 | 2.58 at every ε and ω (K/S = 31/12) |

Contamination carriers (undone frame): singlet — site-1 weight (1,0,2) only at leading order; octet 1 — (1,1,1), (1,2,0), (2,0,1) at the same order.

## 3. Q side at small ε and leading-order CG (Task 2, cells 17–21) [numeric]

Dual Q's of the branches (deg q^j = μ_j) and vacuum direct Q's of the twist (w₁, w₂) (degrees 0, 2, 4) from fixed-degree least squares; Baxter residuals 8·10⁻¹⁴⁹ and 7·10⁻¹⁹⁰. R_Q(B) = det₆[B, Fvac]/Ψ_Q × det₂^{(1)}[B; y]/minor₁. R_C/R_Q is state-independent over the 6 branches and over all 21 states of the HW sectors to 1.6·10⁻¹¹⁶ (limited by 200 digits). The pure Q-side number is Qnum(B) = R_Q(B)/R_Q(top) × ⟨v₁⊗v₂|e₁e₁⟩/⟨x₀|e₁e₁⟩.

Leading order: C^Q = Qnum/(δ^Λ·fac(ω)) with fac = c′₀ ([4,1,1]), a′₀ ([3,3,0]), a′₀c′₀ (octets, giving the functional (u+u′)·C_c). Spread over the 18 ω's and error vs the chain side are both O(ε): 8.5·10⁻⁷ to 3.5·10⁻⁶ at ε = 10⁻⁶, ×4 at 4·10⁻⁶. Richardson limits (C·ε^p, p = 3,3,4,4) agree Q vs chain to ≤ 10⁻¹⁹.

## 4. The O(ε) layer (Tasks 3–7, cells 22–40)

**Top-killing** (exact): subtract δ^Λ s⟨v₂⊗v₂|Ψ⟩/⟨x₀|Ψ⟩, the plain vacuum overlap of the twist (w₁, y), itself det₆ with the (w₁, y) vacuum Q's over Ψ_Q (identity 1.6·10⁻¹¹⁶).

**Singlet** [fitted, validated]. K/S is constant in ε, ω and φ. Scan over 8 values of t at h = 1 (fresh kernel at t = 5/2 and at P2; φ′ at t = 7/3):

  K/S = −t/(3(t−2)),  κ_singlet = 1 + K/S = 2(t−3)/(3(t−2))  (43/12 at P1), residual 3·10⁻¹¹.

Structural reason (main session): the first-order correction ψ₁ of the untwisted eigenvector comes from the twist derivative of t₁(u), whose off-diagonal part within a weight sector is h²Σ_j E⁽¹⁾_kj E⁽²⁾_jk, θ-independent; the singlet signal carries the explicit factor (t−2) from item 2 of the summary.

**Octets** [fitted, validated]. Let F := Qnum^{tk}/(δ^Λ a′₀c′₀) = F₀ + εF₁ + …, with F₀ ω-independent (4.8·10⁻¹⁷). Then F₁/F₀ is exactly a straight line in ω₂ (residual ≤ 10⁻¹¹), ω₁-independent, with a φ-independent slope (spread 8·10⁻¹²) and a φ-dependent intercept. The signal contributes −(t−2)R_c to the slope, R_c := (u·C_c)/((u+u′)·C_c) the mixing of copy c, and the rest is a structural constant Δ_c = slope_c + (t−2)R_c, independent of ω and φ. Scan over 22 values of t (15 fit, 7 held out, residuals ≤ 4·10⁻¹¹):

  Δ_c = (1−t²) + t·τ_c/(3(t+3)),   τ_c = untwisted eigenvalue of t₂(θ₁+h)/h⁴ of copy c,
  τ₁+τ₂ = 3(2t−1)(t+2),  τ₁τ₂ = 9t(t−1)(t+1)(t+3)  (τ_c are the roots of τ² − 3(2t−1)(t+2)τ + 9t(t−1)(t+1)(t+3) = 0),
  Δ₁+Δ₂ = 3(2−t²)/(t+3),  Δ₁Δ₂ = (t−3)(t−1)(t+1)/(t+3).

τ_c is a Q-number: t₂(θ₁+h) = χ₃ν₁(θ₁+h)ν₂(θ₁)Q^j(θ₁−h)/Q^j(θ₁) (any j, j-spread 0). Hence the mixing from the Q side is **R_c^Q = (Δ_c − slope_c^Q)/(t−2)**.

### 4.1 Derivation (`Adjoint-untwist-perturbative.wb`, 22 cells, 42 exact checks) [exact]

First-order perturbation theory Ψ = h + εψ₁ of the untwisted eigenvector (any generic member of the commuting family; two were used), symbolic in t, φ, ω₁, ω₂, with the twist derivative t′ = Σ_I (∂_ε z_I) m_I built from the untwisted principal minors:
- **Singlet.** K₁/S₁ = −t/(3(t−2)) exactly. The twist perturbation is an adjoint tensor operator, so it connects the singlet only to the two octets: 100% of K₁ comes from the octet pair (27 and decuplets give 0), per copy K₁,c/S₁ = A + Bτ_c with A = (3t⁴+2t³−9t²−11t+6)/(3(t−2)(5t²+4)), B = −(t²−2)/(3(t−2)(5t²+4)); the pair sums to −t/(3(t−2)). S₁ = 3(t−2)(ω₁−φ₁)²(φ₁−φ₃)²(φ₂−φ₃)²/((ω₂−φ₂)(ω₁−φ₃)²(ω₂−φ₃)²).
- **Octets.** t₂(θ₁+h) on HW₍₃,₂,₁₎ is [[3(2t²+5t−1)/2, 3(t−1)/2],[3(t−3)/2, 3(t−1)(2t+3)/2]] (trace 3(2t−1)(t+2), determinant 9t(t−1)(t+1)(t+3)); the discriminant 9(5t²+4) is a conic, rationally parametrized by t = 4m/(5−m²), √(5t²+4) = 2+mt, so both copies are treated symbolically at once. F₁/F₀ is exactly affine in ω₂ with ω₁- and φ-independent slope; the signal part of the slope is exactly −(t−2)R_c; the contamination part is Δ_c = (1−t²) + tτ_c/(3(t+3)) exactly, half from each decuplet ([4,1,1] and [3,3,0]), while the 27 and the other octet contribute only to the φ-dependent intercept.
- **Clean branches** have an O(ε) structural correction too: top-killed, K₁^tk/S₀ = (t−1)(φ₂−φ₃)/4 for [4,1,1] and (t−1)(φ₁−φ₂)/4 for [3,3,0] (ω-independent). Not needed for the extraction.
- Cross-check at t = 8 (τ = 252, 198 rational) with 80-digit twisted operators at ε = 10⁻⁸: Δ_c to 3·10⁻⁸.

## 5. Final: all six HW states from the det₆ × det₂ formula alone (cells 41–45) [numeric]

Per branch: max over 18 ω of |C^Q/C − 1| at ε = 10⁻⁶; Richardson limits of C·ε^p (p = 2,3,3,4,4,6) agree Q vs chain to the stated |Q/chain − 1|.

| branch | P1 | P2 | P3 (out of sample) |
|---|---|---|---|
| 27 | 0 / 2.609467456 / 0 | same | same |
| [4,1,1] | 7.9e-7 / 46.67453891 / 1.2e-19 | 1.1e-6 / 38.57793522 / 6.0e-20 | 7.6e-7 / 100.8170040 / 1.8e-19 |
| [3,3,0] | 9.4e-7 / 10.92853892 / 1.6e-121 | 1.1e-6 / 9.032771962 / 9.7e-121 | 4.3e-7 / 23.60564406 / 3.5e-121 |
| octet 1, (u+u′)·C | 3.4e-6 / 184.3305612 / 3.2e-20 | 3.6e-6 / 122.3424387 / 2.8e-20 | 3.1e-6 / 937.1257171 / 3.8e-19 |
| octet 2, (u+u′)·C | 1.8e-6 / −238.5587523 / 9.5e-20 | 1.2e-6 / −245.4884068 / 1.6e-20 | 3.6e-6 / −254.7138085 / 2.9e-19 |
| singlet | 1.5e-6 / −5118.538314 / 4.5e-19 | 2.0e-6 / −3805.202777 / 6.4e-19 | 1.1e-6 / −20997.70633 / 2.2e-19 |

Octet mixing, Q side vs chain side (|R^Q − R^chain|): P1 4.9e-11, 3.7e-11; P2 8.4e-11, 5.6e-11; P3 1.2e-11, 5.3e-12 (R₁ = −1.694, 1.096 at P1). The differences are at the Richardson O(ε³) level.

## 6. Cell map

| cells | content |
|---|---|
| 3–9 | setup: generators + invariant metric, transfer polynomials `tZof`, factorized x₀, `oneSiteVac`, `ρdiag`/`Dy`, u12 machinery, `reg` (the regime), HW kets/projectors, Q machinery and fixed-degree solvers, `branchesAt` |
| 10–12 | Task 0 (ε-structure, exact) |
| 13–16 | Task 1 (S/K exponents, κ, contamination breakdown) |
| 17–21 | Task 2 (Q side at small ε, Qnum, leading-order CG, Richardson) |
| 22–24 | Task 3 (top-killing, ρ_B diagnostics) |
| 26–28 | Task 4 (singlet κ(t) scan; corrected after the purge bug) |
| 29–32 | Task 5 (octet O(ε) slope/intercept over 8 φ) |
| 34–36 | Task 6 (fresh-kernel leak test; fixed `purge[]`, cell 36) |
| 37–39 | Task 7 (Δ_c, τ_c scan over 22 t; fits) |
| 41–45 | final pure-Q-side extraction at P1, P2, P3; cell 45 = FINAL markdown |

Dependencies for a partial rerun: cells 3–9 always; Task 5 needs 11 and 14; the final cells need 3–9, 11, 14, 28, 36.

## 7. Tooling and gotchas

- **No stored outputs.** VS Code did not execute cells of this file (controller not associated, as for `L2-410-two-rotations.wb`); the agent ran each cell's source in the shared kernel K1. Rerun in VS Code to store outputs (a full run is a few minutes; one t point ≈ 3 s, 21 states at one ε ≈ 1.5 s).
- **`Do[DownValues[f] = …, {f, {…}}]` does nothing** (HoldAll); purge memoized definitions with literal symbols. This produced a fake law before it was caught (the constant was P1's t/3).
- The spec's w₃ = 1+ε·2/9 coincided with w₂ at ω₂ = 2/9 (vacuum Baxter residual 4·10⁻¹⁰); w₃ = 1+ε·5/13 is used (w₃ drops out of the identity).
- The naive Baxter residual printed an exact 0 with Precision ∞; it is evaluated in exact arithmetic on the stored numbers.
- `ClearAll["Global`*"]` had to be wrapped in `Quiet` (a Protected `ClearGlobals` symbol aborts the cell).
- The shared kernel: this notebook's first cell wipes the live state of the other new-conventions notebooks (files intact).

## 8. Open / next

1. **Analytic untwisting of det₆ × det₂:** the leading coefficient of the formula in the regime, written through untwisted (staggered-degree) Q's and the direction parameters (φ, ω); the old thread `old-conventions/old/XXX/Experiments/SU3_L2_Untwisting.{md,tex}` has the Casoratian/degree technology. Not started.
2. **[3,1,0], k = 2 point** (`L2-310-two-rotations.wb` point A) through the same pipeline: the dressing e₁ᵀ(1+βE₁₂+γ(E₁₃+½E₁₂E₂₃)) has its own ε-behaviour in the regime; expect the same pattern (leading-order clean set + structural O(ε) constants for the bottom irreps). Then [4,1,0] k = 2.
3. ~~Derive κ_singlet(t) and Δ_c(t, τ_c)~~ — done, §4.1.
4. Store the outputs (rerun in VS Code).
