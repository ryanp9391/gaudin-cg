# Clean CG functionals from localized transfer matrices (log of 2026-10-02)

gl(3) XXX, L=2, new conventions with EH generators. This file is the Markdown companion of `Clean-CG-combinations.tex`/`.pdf`; both have the same content. The notebooks are listed under [Files](#files).

**Chain points.**
- P1, used for everything unless stated: h = 1/2, θ = (2/7, −3/5), so t ≡ θ₁₂/h = 62/35.
- P2, used for validation: h = 1/3, θ = (1/5, −4/9), t = 29/15.

**Evidence levels.**
- *exact*: a symbolic identity, or agreement to the working precision;
- *validated*: fitted at P1 and then predicted at P2;
- *fitted*: numerical at one point;
- *open*: not settled.

## 1. Setup and the question

**Lax and transfer matrices.**
- Lax: ℒ⁽ˢ⁾_ij(u) = (u−θ_s)δ_ij + h E⁽ˢ⁾_ij, and T(u) = ℒ⁽¹⁾ℒ⁽²⁾ diag(z).
- t_a: sums of principal quantum minors.
- CBR fusion: T_λ(u) = det t_{λ'_j+i−j}(u+(i−1)h).
- ν_n(x) = ∏_s(x−θ_s+hΛ_n).

**Q-functions.**
- Direct: Q_i = z_i^{−u/h} q_i, with deg q_i = 2Λ₁ − μ_i.
- Dual: Q^j = z_j^{u/h} q^j, with deg q^j = μ_j (Λ₃ = 0).

**The bra.** ⟨Ω| = ⟨Ω_zw|⊗⟨Ω_zw|, with Ω_zw = e₁ᵀR and R = e^{r₁E₁₂}e^{r₂E₁₃}e^{r₃E₂₃}. The parameters are
- r₁ = (w₁−z₁)/(w₁−z₂);
- r₂ = (w₁−z₁)(z₁−w₂)z₂₃/(z₁₂(w₁−z₃)(w₂−z₃));
- r₃ = (w₁−z₂)(z₂−w₂)z₁₃/(z₁₂(w₁−z₃)(z₃−w₂)).

They do not depend on the representation, and w₃ drops out.

**Normalization and twists.** ⟨x₀|Ψ⟩ = 1, with x₀ = (e₁ᵀP_phys⁻¹)^{⊗2}. The twists are z = e^{εφ} and w = e^{εω}. HW-limit states behave as Ψ_B → ε^{−k}h_μ.

**Q-side data.** For any insertion X built from t's, ⟨Ω|X|Ψ_B⟩ = X^B⟨Ω|Ψ_B⟩. The factor ⟨Ω|Ψ_B⟩ is the 6×6 determinant formula of `L2-diagonal-frame-overlap.md`, which has been validated for the adjoint only.

**The question.** Can the CG data be read off directly from a Q-expression, without reconstructing the full vector?

## 2. The contamination

### 2.1 The mechanism (adjoint; `Minimal-ingredients-EH.wb`, Sections 5–6)

Undoing R gives ⟨Ω|T_λ(θ₁−2h) = f_λ c_λ(R⊗R), with c_λ = ⟨v_λ|R⁻¹⊗e₁ᵀ. Since (R⊗R)P_HW = P_HW exactly, there is an exact split at any ε:

  T_λ^B⟨Ω|Ψ⟩/f_λ = c_λ·P_HWΨ (**signal** S) + c_λ·(R⊗R)(1−P_HW)Ψ (**contamination** K), with κ = (S+K)/S.

**Why it happens.**
- At z = 1, c_λ ∝ e₁⊗e₁, which is pure top weight. Its lower-weight components are O(ε).
- Write c = a₀ + εa₁ and Ψ = ψ₀ + εψ₁. Then c·Ψ = ε(a₁ψ₀ + a₀ψ₁): the second term is the O(1) top part of c meeting the O(ε) correction ψ₁ raised by R.

**Measured κ values.**
- [3,3,0], and [4,1,1] with (2,1): κ = (t+1)/2 = 97/70.
- Octets: (1.530, 1.072) and (0.642, 0.772).
- Singlet with (2,1): κ ≈ 1.102 = (t+1)(t+3)/12.

For [3,3,0], all of the contamination is the ν = [4,2,0] term: ‖c^{(420)}‖ ~ ε⁰ times ‖[(R⊗R)Ψ⊥]₄₂₀‖ ~ ε⁻². The signal is ε¹·ε⁻³, and the ratio is 27/70 = κ − 1.

### 2.2 Why projectors and w-expansions cannot remove it

These are arguments, plus one numerical counterexample.
- **Weight projectors.** Ψ has a fixed weight and T_λ commutes with the Cartan, so a weight projector on the left-hand side is trivial. With R undone, it would need (R⊗R)⁻¹ to act on Ψ, which is not a Q-quantity.
- **Torus degree.** The torus acts on (r₁,r₂,r₃) with weights (α₁, α₁+α₂, α₂). The component c^{(ν)} has degree top−ν, and [(R⊗R)Ψ]_ν has degree ν−μ. Every term therefore has the same degree, top−μ, so no w-expansion separates the CG term from the [4,2,0] term.
- **Counterexample.** For [3,3,0], the left-hand side is exactly ∝ r₁, a single monomial, and its coefficient is κ × CG.
- **A third independent r adds nothing.** Modulo the torus there is one invariant, ρ = r₂/(r₁r₃) = −(z₁−w₂)z₂₃/((z₂−w₂)z₁₃). It already varies with w₂, so the missing direction is a torus direction, which only rescales by a known factor.
- **The bra only sees the top irrep.** Every coherent bra (e₁ᵀg)^{⊗2} lies in the dual of the top irrep, so all access to other irreps comes from the insertions.

## 3. Clean combinations

### 3.1 Top-killing is exact, and is an eigenvalue difference

For μ ≠ top, e₁e₁·P_HWΨ = 0, so ⟨Ω|Ψ⟩ = e₁e₁·(R⊗R)Ψ⊥. The top-weight contamination of c_λ is therefore exactly (c_λ·e₁e₁)⟨Ω|Ψ⟩. In the bare frame, c_λ·e₁e₁ = ⟨Ω|X|e₁e₁⟩ = X^top. So:

  **top-killed Q-number = (X^B − X^top)·⟨Ω|Ψ_B⟩**

In the bare frame, with W_X = ⟨Ω|X:
- S_X = W_X·P_HWΨ;
- K_X = W_X·Ψ⊥ − X^top⟨Ω|Ψ⟩.

No R⊗R is needed. This holds for any representation.

**The clean-combination recipe.**
1. In each weight sector, take the leading ratios k_B = K/S_ref of every HW-limit copy.
2. The clean combinations are the null space of these k_B's.
3. Least squares over several w then gives P_HWΨ, scale included.

For all three representations studied, the k's are **structural**: independent of φ and ω up to O(ε).

### 3.2 Adjoint [2,1,0]² (`Minimal-ingredients-EH.wb`, Section 7)

**Top-killing alone** gives κ → 1 for [4,1,1], [3,3,0] and both octets.
- The two top-killed T's have rank 2 on the octet HW space.
- So the octet mixing is fixed by site 1 + ⟨Ω|Ψ⟩.

**The singlet.**
- At P1: K_b/S = −18/35 and K_c/S = 19/105.
- This gives β = b/c = h(2t−3)/(t−1) = 19/54 in the f-normalized basis.
- P2 was predicted correctly (−14/15, 13/45, 13/42): **validated**.

**One combination for everything.** v* = (a, β, 1), with a = −3h²z₂z₃(β − hz₃), gives κ → 1 for all five non-top HW states, at P1 and at P2. Section 6 derives this analytically.

### 3.3 [3,1,0]² (`Minimal-ingredients-310-EH.wb`)

**Localizing shapes.** An independent 120-digit scan (Appendix A, shapes with λ₁ ≤ 4, shifts c = −6..2) finds that the only nontrivial site-1-local shapes are T₍₁,₁₎, T₍₂,₂₎ and T₍₃,₂₎, all at θ₁−3h.
- T₍₃,₃₎ and all shapes with λ₁ = 4 vanish on ⟨Ω| there.
- The localization holds to 10⁻¹⁰⁶.

**Site 1 + Ω with sector-wide cleaning** fixes every multiplicity space except [4,3,1]:
- [6,1,1], [5,3,0] and [5,2,1] are clean after top-killing.
- v₁ = 3d₁₁ − 2d₂₂ is clean on [5,3,0], [5,2,1]×2, [4,4,0], [4,3,1]×2 and [4,2,2].
- [5,2,1] has rank 2, so its mixing is fixed.
- **[4,3,1] is not fixed.** Only v₁ is clean on both copies, and v₁'s direction on HW_[4,3,1] does not depend on w (exact rank 1).

**Adding site 2** (A = t₂(θ₂+h) and B = t₂(θ₂+h)t₁(θ₂−h), computed as ⟨Ω|X(R⊗R)⁻¹ so that no prefactors are needed) fixes everything:
- All ratios are structural.
- The HW rank is 2 on both [5,2,1] (second singular value 0.16) and [4,3,1] (0.05).
- **End to end:** all 9 non-top HW states are recovered, scale included. Errors are O(ε): about 10⁻⁵ at ε = 10⁻⁴ and at most about 10⁻⁶ at ε = 10⁻⁵.

**Site 1 only, with per-copy cleaning** (Section 9), also recovers all 9 states. Errors are 10⁻⁴ to 10⁻⁶·⁵ with at least 62 correct digits. Only [4,3,1] needs copy-specific coefficients.

### 3.4 [4,2,0]²: the "blind spot" was wrong (`Blindspot-420-check-EH.wb`)

`L2-420-localization.md` claimed that HW_[5,4,3] and HW_[4,4,4] lie exactly in the null space of the localized data. That is **false**.

**Why the overlaps can't vanish structurally.** A site-1-local covector is u⊗e₁ᵀ, and c·h = u·(1⊗e₁ᵀ)h. The slice (1⊗e₁ᵀ)h is never zero, because the top site-2 weight must occur in h.

**What the localized covectors actually see.**
- T₍₃,₂₎ and T₍₄,₂₎ have components at site-1 weight (1,2,3); T₍₄,₂₎ also has one at (0,2,4).
- The overlaps are about 8×10⁻³ on [5,4,3] and 8×10⁻⁴ on [4,4,4] at the finite twist.
- At small twist they scale exactly as ε³ and ε⁴, while ‖u‖ ≈ 2.5.

**What the old test showed.** Its "overlap singular value = 1.000" was printed to only 3 digits, and only means near-containment in a numerical null space.

### 3.5 [4,2,0]² clean combinations (`Minimal-ingredients-420-EH.wb`, 250 digits)

**Insertions:**
- site 1: T₍₁,₁₎, T₍₂,₂₎, T₍₃,₂₎ and T₍₄,₂₎ at θ₁−4h;
- site 2: t₂(θ₂+h), T₍₂,₂₎(θ₂) and T₍₂,₂₎(θ₂)t₁(θ₂−2h);
- plus ⟨Ω|Ψ⟩.

**Leading order is enough.**
- For every state, each insertion's top-killed contamination has the same ε-power as its signal, or a higher one.
- This includes [5,4,3] (‖Ψ‖ ~ ε⁻¹⁰) and [4,4,4] (‖Ψ‖ ~ ε⁻¹²): the ε^{3,4} suppression of their overlaps is compensated by the growth of Ψ.

**Structural.** The spread over twist directions is 1.4×10⁻³ at ε = 10⁻⁵ and 1.4×10⁻⁴ at ε = 10⁻⁶.

**[6,4,2] (multiplicity 3) is the only problem for sector-wide cleaning.**
- Every other sector's clean space has full HW rank.
- For [6,4,2], the three copies' k's span only 2 dimensions, and the common clean rank is 2 (the third singular value is ∝ ε).
- The raw insertions do have rank 3 on it (1, 0.11, 0.0074).

**Fair end-to-end test.**
- Setup: clean spaces fitted at (ε = 10⁻⁶, φ, ω), then tested at ε = 10⁻⁴ and 10⁻⁵ and at another φ.
- Result: all 18 non-top HW states recovered, scale included, with O(ε) errors and at least 115 correct digits.
- [6,4,2] uses per-copy cleaning.

**Site 1 only** (Section 3) also recovers all 18 states: errors 10⁻⁴·⁶ to 10⁻⁷·¹, at least 118 digits.
- Sector-wide cleaning works for every multiplicity-1 sector, [7,4,1] and [5,4,3].
- [7,3,2], [6,5,1] and [6,4,2] need per-copy cleaning.
- **Site 2 is not essential.** It only reduces the number of copy-specific sectors.

### 3.6 [3,2,0]², [4,1,0]², [4,3,0]² and the λ₁ = 5 reps, site 1 only (`Minimal-ingredients-{320,410,430,510,540,520,530}-EH.wb`)

The two notebooks run the same Λ-generic code: the site-1-only pipeline of 3.3 and 3.5. They differ in two ways.
- The site-1 insertions come from an independent localization scan inside the notebook, not from earlier notebooks.
- The covectors are computed directly as ⟨Ω|T_λ(R⊗R)⁻¹, so no f_λ are needed. The Q-number is T_λ^B⟨Ω|Ψ⟩, with the eigenvalue built by CBR from scalar t-eigenvalues.

Setup: P1, fitted at (ε = 10⁻⁵, φ, ω), tested at (10⁻⁴, φ) and (10⁻⁵, φ′).

**[3,2,0]²** (150 digits, about 3 s):
- **Insertions:** T₍₁,₁₎, T₍₂,₁₎, T₍₃,₁₎ at θ₁−3h, as in `L2-320-localization.md`. The fitted f's match the known ones.
- **All 9 non-top HW states recovered**, scale included. Errors are 10⁻⁴·¹–10⁻⁵·⁷ at 10⁻⁴ and 10⁻⁵–10⁻⁷·¹ at 10⁻⁵, with at least 85 digits.
- **Only [5,3,2] needs per-copy cleaning** (sector clean dimension 1, HW rank 1). This mirrors [4,3,1] in [3,1,0], consistent with [3,2,0] ≅ [3,1,0]*⊗det³.
- [5,4,1] (multiplicity 2) works sector-wide, and so does the bottom irrep [4,3,3].

**[4,1,0]²** (200 digits, about 25 s). No localization had been done for this rep before.
- **Site-1 local** (λ₁ ≤ 4, c = −6..1): exactly **T₍₁,₁₎, T₍₂,₂₎, T₍₃,₃₎, T₍₄,₃₎ at θ₁−4h**.
  - That is four shapes. The pattern "T₍₁,₁₎ + T₍ₖ,λ₁−λ₂₎, k = λ₁−λ₂..λ₁" misses T₍₂,₂₎.
  - Prefactors (fitted at P1, unique within the ansatz ±∏ₖ(θ₁₂−kh)^{nₖ}; checked at P2 to 10⁻¹⁹⁰):
    - f₍₁,₁₎ = (θ₁₂−3h)(θ₁₂−5h);
    - f₍₂,₂₎ = (θ₁₂−2h)(θ₁₂−3h)(θ₁₂−4h)(θ₁₂−5h);
    - f₍₃,₃₎ = (θ₁₂−h)(θ₁₂−2h)(θ₁₂−3h)²(θ₁₂−4h)(θ₁₂−5h);
    - f₍₄,₃₎ = (θ₁₂−h)²(θ₁₂−2h)(θ₁₂−3h)²(θ₁₂−4h)(θ₁₂−5h).
- **All 14 non-top HW states recovered**, scale included. Errors are 10⁻³·⁸–10⁻⁵·² at 10⁻⁴ and 10⁻⁴·⁷–10⁻⁶·³ at 10⁻⁵, with at least 123 digits.
- **Only [6,3,1] needs per-copy cleaning** (sector clean dimension 2, HW rank 1).
- [7,2,1] works sector-wide, but is weakly conditioned: its second singular value is 0.0069, the same at both ε. [5,4,1] and the bottom irrep [4,4,2] also work sector-wide.

**[4,3,0]²** (`Minimal-ingredients-430-EH.wb`, 200 digits, about 25 s). No localization had been done for this rep before.
- **Site-1 local:** exactly **T₍₁,₁₎, T₍₂,₁₎, T₍₃,₁₎, T₍₄,₁₎ at θ₁−4h**. Here the pattern holds.
- **Prefactors** (fitted at P1, unique; checked at P2 to 10⁻¹⁹²):
  - f₍₁,₁₎ = (θ₁₂−h)(θ₁₂−5h);
  - f₍₂,₁₎ = (θ₁₂−h)(θ₁₂−3h)(θ₁₂−5h);
  - f₍₃,₁₎ = (θ₁₂−h)(θ₁₂−2h)(θ₁₂−3h)(θ₁₂−5h);
  - f₍₄,₁₎ = (θ₁₂−h)²(θ₁₂−2h)(θ₁₂−3h)(θ₁₂−5h).
- **All 14 non-top HW states recovered**, scale included. Errors are 10⁻³·⁹–10⁻⁵·³ at 10⁻⁴ and 10⁻⁴·⁹–10⁻⁶·⁶ at 10⁻⁵, with at least 122 digits.
- **Only [7,5,2] needs per-copy cleaning** (sector clean dimension 2, HW rank 1).
  - Under [4,3,0] ≅ [4,1,0]*⊗det⁴, μ ↦ (8−μ₃, 8−μ₂, 8−μ₁) sends [6,3,1] ↦ [7,5,2]. So the per-copy sector is mirrored, as for [3,2,0] ↔ [3,1,0].
  - The conditioning is not mirrored: [7,6,1], the image of [7,2,1], has second singular value 0.13 rather than 0.0069.

**[5,1,0]²** (`Minimal-ingredients-510-EH.wb`, 19 cells, 250 digits; Sections 1–6 about 2.5 min). This is the first rep with λ₁ = 5 (dim 35, 1225 states; 18 HW states).
- **Site-1 local** (λ₁ ≤ 5, c = −7..1): exactly **T₍₁,₁₎, T₍₂,₂₎, T₍₃,₃₎, T₍₄,₄₎, T₍₅,₄₎ at θ₁−5h**.
  - The rule for [λ₁,1,0] is T₍ₖ,ₖ₎ for k = 1..λ₁−1, plus T₍λ₁,λ₁−1₎, all at θ₁−λ₁h. It holds for λ₁ = 2..5.
- **Prefactors** (with t = θ₁₂/h; unique fits at P1, checked at P2 to 10⁻²³⁹):
  - f₁₁ = (t−4)(t−6);
  - f₂₂ = (t−3)(t−4)(t−5)(t−6);
  - f₃₃ = (t−2)(t−3)(t−4)²(t−5)(t−6);
  - f₄₄ = (t−1)(t−2)(t−3)²(t−4)²(t−5)(t−6);
  - f₅₄ = (t−1)f₄₄.
- **All 17 non-top HW states recovered**, scale included.
  - Fit at 10⁻⁵: errors 10⁻³·⁷–10⁻⁵·⁸ at (10⁻⁴, φ) and 10⁻⁴·²–10⁻⁶·³ at (10⁻⁵, φ′), with at least 161 digits.
  - Refit at 10⁻⁶ (Section 6b): errors 10⁻⁵·²–10⁻⁷·³ at (10⁻⁶, φ′). They fall by about a decade, so the stalls at 10⁻⁵ came from the fit point and convergence is O(ε).
- **Only [8,3,1] needs per-copy cleaning** (sector clean dimension 3, HW rank 1).
  - [9,2,1] works sector-wide but is weakly conditioned (second singular value 0.0036, the same at both ε), as [7,2,1] was for [4,1,0].
- **For [λ₁,1,0], the per-copy sector is always [λ₁+3,3,1]:** [4,3,1], [6,3,1], [8,3,1] for λ₁ = 3, 4, 5.

**[5,4,0]²** (`Minimal-ingredients-540-EH.wb`, 250 digits, about 3 min; 18 HW states).
- **Site-1 local:** T₍₁,₁₎, T₍₂,₁₎, T₍₃,₁₎, T₍₄,₁₎, T₍₅,₁₎ at θ₁−5h.
- **Prefactors** (t = θ₁₂/h; checked at P2 to 10⁻²⁴¹):
  - f₁₁ = (t−1)(t−6);
  - f₂₁ = (t−1)(t−4)(t−6);
  - f₃₁ = (t−1)(t−3)(t−4)(t−6);
  - f₄₁ = (t−1)(t−2)(t−3)(t−4)(t−6);
  - f₅₁ = (t−1)f₄₁.
- **All 17 non-top HW states recovered.**
  - Fit at 10⁻⁵: 10⁻³·⁷–10⁻⁵·¹ at (10⁻⁴, φ), 10⁻⁴·⁷–10⁻⁶·⁴ at (10⁻⁵, φ′), with at least 160 digits.
  - Refit at 10⁻⁶: 10⁻⁵·⁷–10⁻⁷·⁴ at (10⁻⁶, φ′).
- **Only [9,7,2] needs per-copy cleaning.** It is the dual image of [5,1,0]'s [8,3,1], as predicted.

**[5,2,0]²** (`Minimal-ingredients-520-EH.wb`, 300 digits, about 9 min; 28 HW states in 18 sectors, two of multiplicity 3).
- **Site-1 local:** T₍₁,₁₎, T₍₂,₂₎, T₍₃,₃₎, T₍₄,₃₎, T₍₅,₃₎ at θ₁−5h.
- **Prefactors** (checked at P2 to 10⁻²⁸⁹):
  - f₁₁ = (t−3)(t−6);
  - f₂₂ = (t−2)(t−3)(t−5)(t−6);
  - f₃₃ = (t−1)(t−2)(t−3)(t−4)(t−5)(t−6);
  - f₄₃ = (t−1)(t−2)²(t−3)(t−4)(t−5)(t−6);
  - f₅₃ = (t−1)²(t−2)²(t−3)(t−4)(t−5)(t−6).
- **All 27 non-top HW states recovered.**
  - Fit at 10⁻⁵: 10⁻³·⁰–10⁻⁵·⁷ at (10⁻⁴, φ), 10⁻³·³–10⁻⁷·¹ at (10⁻⁵, φ′), with at least 199 digits.
  - Refit at 10⁻⁶: 10⁻⁴·³–10⁻⁸·¹ at (10⁻⁶, φ′).
- **Per-copy cleaning required** for [9,3,2], [8,5,1], [8,4,2] (multiplicity 3) and [7,5,2] (multiplicity 3).
- **New: weakly conditioned sectors (Section 6c).**
  - [7,6,1] (σ₂ = 0.025) and [7,4,3] (σ₂ = 0.019) have full sector rank, but their worst copy is amplified by ~1/σ₂: 10⁻⁵·¹ and 10⁻⁴·³ at ε = 10⁻⁶, still O(ε).
  - Per-copy cleaning brings them to 10⁻⁶·² and 10⁻⁵·⁷. With that choice every state is ≤ 10⁻⁵·⁵ at ε = 10⁻⁶.
  - Lesson: when σ_min ≲ 0.05, prefer per-copy cleaning.

**[5,3,0]²** (`Minimal-ingredients-530-EH.wb`, 300 digits, about 8.5 min; 28 HW states in 18 sectors, two of multiplicity 3).
- **Site-1 local:** T₍₁,₁₎, T₍₂,₂₎, T₍₃,₂₎, T₍₄,₂₎, T₍₅,₂₎ at θ₁−5h.
- **Prefactors** (checked at P2 to 10⁻²⁹⁰):
  - f₁₁ = (t−2)(t−6);
  - f₂₂ = (t−1)(t−2)(t−5)(t−6);
  - f₃₂ = (t−1)(t−2)(t−3)(t−5)(t−6);
  - f₄₂ = (t−1)(t−2)²(t−3)(t−5)(t−6);
  - f₅₂ = (t−1)²(t−2)²(t−3)(t−5)(t−6).
- **All 27 non-top HW states recovered.**
  - Fit at 10⁻⁵: 10⁻³·⁷–10⁻⁵·⁸ at (10⁻⁴, φ), 10⁻⁴·⁰–10⁻⁷·⁴ at (10⁻⁵, φ′), with at least 198 digits.
  - Refit at 10⁻⁶: 10⁻⁵·⁰–10⁻⁸·⁴ at (10⁻⁶, φ′).
  - With per-copy cleaning for the weak [7,6,3] (σ₂ = 0.038), every state is ≤ 10⁻⁵·⁵.
- **Per-copy cleaning required** for [9,5,2], [8,7,1], [8,6,2] and [8,5,3]. These are exactly the dual images of [5,2,0]'s four.

### 3.7 Summary across all ten reps (site 1 + ⟨Ω|Ψ⟩ only)

Every λ₁ ≤ 5 rep with λ₁ > λ₂ > 0 has been run: [2,1,0], [3,1,0], [3,2,0], [4,1,0], [4,2,0], [4,3,0], [5,1,0], [5,2,0], [5,3,0], [5,4,0]. **In every rep, all non-top HW CG states are recovered from Q-side numbers, scale included, with O(ε) errors.**

**Site-1 insertions.** Write d = λ₁−λ₂. The site-1-local shapes are:
- T₍ₖ,ₖ₎ for k = 1..d;
- T₍ₖ,d₎ for k = d+1..λ₁;
- all at θ₁−λ₁h.

That is λ₁ insertions in total. Independent scans confirm the rule in all ten reps, and nothing else is site-1 local in the scanned range. The leading prefactor is always f₍₁,₁₎ = (θ₁₂−dh)(θ₁₂−(λ₁+1)h).

**Sectors needing per-copy cleaning** (sector clean space HW-rank deficient):

| Λ | per-copy sectors |
|---|---|
| [2,1,0] | none |
| [3,1,0] / [3,2,0] | [4,3,1] / [5,3,2] |
| [4,1,0] / [4,3,0] | [6,3,1] / [7,5,2] |
| [4,2,0] | [7,3,2], [6,5,1], [6,4,2] |
| [5,1,0] / [5,4,0] | [8,3,1] / [9,7,2] |
| [5,2,0] / [5,3,0] | [9,3,2], [8,5,1], [8,4,2], [7,5,2] / [8,7,1], [9,5,2], [8,6,2], [8,5,3] |

- **Duality.** Under [λ₁,λ₂,0] ↦ [λ₁,λ₁−λ₂,0] (dual ⊗ det^λ₁), the per-copy sectors map exactly, by μ ↦ (2λ₁−μ₃, 2λ₁−μ₂, 2λ₁−μ₁).
- **[λ₁,1,0].** The per-copy sector is always [λ₁+3,3,1].
- **Weak conditioning is not dual-invariant.** It shows up as small σ₂ in some sectors that are still full rank. Examples: [7,2,1] in [4,1,0]; [9,2,1] in [5,1,0]; [7,6,1] and [7,4,3] in [5,2,0]; [7,6,3] in [5,3,0]. In such sectors, per-copy cleaning is the safer choice.

### 3.8 Closed forms of the clean coefficients (`Clean-coeffs-{210,310,410}-EH.wb`)

**No universal combination exists** (`Minimal-ingredients-530-EH.wb` §6d): in [5,3,0] the 23 contaminated states' k's have full rank 5. So the clean coefficients depend on μ, and in per-copy sectors on the copy.

**Method:**
- h = 1, θ = (t, 0);
- f-free one-site covectors, top-killed;
- leading k̂, ŝ, ρ = K_ref/S_ref and λ_c = t₂(u₀) eigenvalue at ε₀, 2ε₀, 4ε₀ (ε₀ = 10⁻⁸), Richardson-extrapolated to ε → 0 (error ≤ 10⁻²³);
- rational fits in t over 20 points, validated on 5 held-out points;
- multiplicity sectors: tested in the form A(t) + B(t)λ_c.

**Control, [2,1,0].** The pipeline reproduces the analytic adjoint result exactly: k̂₁₁ = −(t−1)/(2t−3) and ρ = (2t−3)/3.

**[3,1,0]: everything is in closed form** (validated to ≤ 4×10⁻¹⁹).
- [4,4,0]: k̂ = (−1/6, −1, 1), ρ = t−1.
- [4,2,2]: k̂₁₁ = −(5t−11)/(24(2t−5)), k̂₂₂ = 6k̂₁₁, ρ = (2t−5)/5.
- Bottom irrep [3,3,2]: k̂₁₁ = −(t−1)(3t−7)/(4(7t²−12)), k̂₂₂ = −3(t−1)(3t−2)/(2(7t²−12)), ρ = (7t²−12)/20.
- [4,3,1]: copies labelled by σ_c = ±1:
  - k̂₁₁,c = −[(35t²−104t+53) + σ_c(5t−13)√(1+3t²)] / (12(25t²−65t+24));
  - k̂₂₂ = 6k̂₁₁.
- **Universal on everything but the bottom:** 6d₁₁ − d₂₂, the closed form of the earlier fitted v₁.

**[4,1,0]: everything is in closed form** (one marginal fit).
- **Layer 1** ([6,4,0], [6,3,1]×2, [6,2,2]): exactly k̂₂₂ = 16k̂₁₁ and k̂₃₃ = 72k̂₁₁.
  - So 16d₁₁ − d₂₂ and 72d₁₁ − d₃₃ are t- and copy-independent clean combinations here.
  - [6,3,1] copies: k̂₁₁,c = −[(189t²−767t+445) + σ_c(7t−26)√(25+24t²)] / (360(49t²−182t+75)).
- **Layer 2** ([5,5,0], [5,4,1]×2, [5,3,2]): the common t-independent clean combination is exactly v₂ = 72d₁₁ − 9d₂₂ + d₃₃.
- **Bottom irrep [4,4,2]:** genuinely t-dependent, with P = 11t³+61t−90:
  - k̂ = (−(t−1)(t−2)(2t−7)/(18P), −2(t−1)²(8t−19)/(9P), −2(t−1)(4t²−5t+9)/P, 1);
  - ρ = P/90.

**Multiplicity discriminants** (each copy quantity is P(t) + Q(t)σ_c√D, i.e. linear in the copy's transfer eigenvalue):

| Λ | sector | D |
|---|---|---|
| [3,1,0] | [5,2,1] | (256+105t²)/4 |
| [3,1,0] | [4,3,1] | 16(1+3t²) |
| [4,1,0] | [7,2,1] | 9(100+21t²)/4 |
| [4,1,0] | [6,3,1] | 4(25+24t²) |
| [4,1,0] | [5,4,1] | 25(4+21t²)/4 |

**Structure:**
- Contamination is organized in layers. Each layer's clean combinations are exactly those that annihilate the signal patterns of the shallower layers. Those combinations are t-independent until the bottom irrep, whose coefficients are rational in t.
- The mechanism is tested where it is non-trivial: k̂_B lies in the span of the signal directions of the HW states above B, to 10⁻²⁵ ([6,4,0] and [5,5,0] in [4,1,0]).

**The √D coefficients as a transfer eigenvalue** (`Clean-coeffs-310-EH.wb` §4).
- Let τ_c be the untwisted eigenvalue of t₂(u) at u = θ₁ + h. Normalized per unit signal, the copy-dependent coefficients are affine in τ_c with simple coefficients.
- On [4,3,1]:
  - K₃₂/S = (19t(t−1) − τ_c)/(27t);
  - K₁₁/S = (4t(t−1)(3t−23) − (3t−11)τ_c)/(1296t);
  - ŝ₂₂ = (4t(t−1) − τ_c)/(24t).
- These are exact at all 25 t, held-out ones included, to ≤ 10⁻²².
- **Per-copy clean combination:** v_c = 48[19t(t−1) − τ_c] d₁₁ − [4t(t−1)(3t−23) − (3t−11)τ_c] d₃₂.
- τ_c = 4[t²+2t−1 ± √(1+3t²)].
- k̂₁₁ = K₁₁/K₃₂ is a Möbius function of τ_c, hence its messier form.
- No target is proportional to a single eigenvalue or to a ratio of two, among 39 candidates.
- t₁ never separates the copies.

## 4. The site-1 eigenvalues as single Q-functions

### 4.1 Reduction to 2×2 determinants

Let v(u) = (Q¹,Q²,Q³)(u). The quantum Weyl formula reads

  T_λ(x) = P_λ(x) · det[v(x+λ₁h), v(x+(λ₂−1)h), v(x−2h)] / det[v(x), v(x−h), v(x−2h)],  with P_λ = ∏_{(r,c)∈λ} ν₃(x+(c−r)h).

**Two degeneracies near θ₁.**
- **The 3×3 Casoratian.** The dual Baxter equation gives W(x)/W(x−h) = χ₃ν₁(x)ν₂(x−h)/(ν₃(x)ν₃(x−h)). So det[v(x),v(x−h),v(x−2h)] vanishes at x = θ_s+h, …, θ_s−(Λ₁−1)h.
- **Dual parallelism:** v(θ₁−Λ₂h) ∥ … ∥ v(θ₁), from the dual QQ relation in 4.2.

Together they put v(θ₁−(Λ₁+1)h), …, v(θ₁) in a single 2-plane. Hence

  **T_λ(θ₁−Λ₁h) = P_λ · det₂[v(θ₁+(λ₁−Λ₁)h), v(θ₁+(λ₂−Λ₁−1)h)] / det₂[v(θ₁−Λ₁h), v(θ₁−(Λ₁+1)h)]**,  for λ₁ ≤ Λ₁ and λ₃ = 0,

with any fixed pair of rows. This is the L=2 version of the L=1 "reduction at θ−2h".

**Checked** exactly on 34 states of [3,1,0]² (to 10⁻⁷⁰) and 91 states of [4,2,0]² (to 10⁻¹⁹⁴), for all shapes and every row pair.

### 4.2 The dual QQ relation and single Q's

For [4,2,0], all 91 states satisfy (to at least 226 digits; the shifts 0 and −1 fail):

  **z_b q^a(u)q^b(u+h) − z_a q^b(u)q^a(u+h) = (z_b−z_a) G(u) q_c(u+h),  G(u) = ∏_s ∏_{n=1}^{Λ₂}(u−θ_s+nh).**

**Consequences.**
- The direct vectors (Q₁,Q₂,Q₃) are parallel at θ₁−4h, θ₁−3h and θ₁−2h.
- Each consecutive det₂ becomes a single direct Q.
- The localization factors f_λ cancel all θ₁₂-dependence.

**[4,2,0] at θ₁−4h** (to 10⁻²²⁵, for all c and j):
- T₍₁,₁₎/f₍₁,₁₎ = 10h²χ₃ r₁
- T₍₂,₂₎/f₍₂,₂₎ = 40h⁴χ₃² r₁r₂
- T₍₃,₂₎/f₍₃,₂₎ = −80h⁵χ₃² r₁r₂ρ₁
- T₍₄,₂₎/f₍₄,₂₎ = 80h⁶χ₃² r₁r₂ρ₀

Here r₁ = Q_c(θ₁−3h)/Q_c(θ₁−4h), r₂ = Q_c(θ₁−2h)/Q_c(θ₁−3h), and ρ_k = Q^j(θ₁−kh)/Q^j(θ₁−2h).

**Adjoint at θ₁−2h** (to 10⁻⁵⁰, 18 states):
- T₍₁,₁₎ = 3h²χ₃(θ₁₂−h)(θ₁₂−3h) r
- T₍₂,₁₎ = −3h³χ₃(θ₁₂−h)²(θ₁₂−3h) ρr

Here r = Q_c(θ₁−h)/Q_c(θ₁−2h) and ρ = Q^j(θ₁)/Q^j(θ₁−h). For the top state, r = 1/z₁ and ρ = z₃, which reproduces the known values 3h²z₂z₃ and −3h³z₂z₃².

## 5. Clean combinations in Q-language ([4,2,0], Section 6)

In the basis X_λ = T_λ/f_λ (site 1), normalized to the T₍₄,₂₎ entry, the leading signal ratios are simple rationals:

| state | s̃ over (T₁₁, T₂₂, T₃₂, T₄₂)/f |
|---|---|
| [8,3,1] | (0, 0, −1, 1) |
| [8,2,2] | (0, 0, 0, 1) |
| [7,5,0] | (1, 2, −2, 1) |
| [6,6,0] | (0, 2, −2, 1) |
| [6,3,3] | (0, 0, −1/5, 1) |
| [5,5,2] | (0, 0, −4/5, 1) |
| [4,4,4] | (0, 0, 0, 1) |

**Contamination ratios.**
- **Zero** for [8,3,1] and [7,5,0]. Their only neighbour above is the top, which is killed exactly.
- **[8,2,2]:** k̃ = (t−1)·s̃^[8,3,1], so the clean condition is c₃₂ = c₄₂.
- **[6,6,0]:** k̃ = (t−1)·s̃^[7,5,0], so the clean condition is c₁₁ + 2c₂₂ − 2c₃₂ + c₄₂ = 0.
- **Pattern:** the contamination is the signal pattern of the irrep one simple root above, times (t−1). This holds in normalized form, at each w; the raw coefficient depends on w.
- **Deeper states** ([6,3,3], [5,5,2], [5,4,3], [4,4,4]): k̃ is not a recognizable rational at P1. *Open.*

## 6. The adjoint clean combination, derived analytically (`Adjoint-clean-analytic.wb`)

Method: exact computation in the E basis, with t symbolic and h = 1; h is restored by homogeneity.

**Why first order is enough.** Ω(ε) always lies in the 27, where the untwisted T_λ acts by T_λ^top. So (W_λ − T_λ^top Ω) vanishes at ε⁰ on the weight sector (checked). Writing Ψ = h + εψ₁ + …:
- S_λ = Ω(ε)·T_λ(ε)h;
- K_λ = ε²(W_λ − T_λ^top Ω)₁·ψ₁ + O(ε³).

So only first-order perturbation theory is needed.

**Consequences.**
- [4,1,1], [3,3,0] and the octets have signals of order ε, so they are clean automatically.
- The singlet has S₍₁,₁₎ ≡ 0 and S₍₂,₁₎ = O(ε²). In raw eigenvalues, exactly and for both twist-direction pairs:
  - K₍₁,₁₎/S₍₂,₁₎ = −1/3;
  - K₍₂,₁₎/S₍₂,₁₎ = (2t−3)/3.

**Result.** A single combination is clean for all six HW states:

  **𝒞 = T₍₂,₁₎(θ₁−2h) + h(2θ₁₂−3h)·T₍₁,₁₎(θ₁−2h),  with Q-number (𝒞^B − 𝒞^top)·⟨Ω|Ψ_B⟩.**

This is β = h(2t−3)/(t−1) in the f-normalized basis, i.e. the earlier fit, now derived.

**In Q's.** Dropping the state-independent factor 3h³χ₃(θ₁₂−h)(θ₁₂−3h):

  **{ r[(2θ₁₂−3h) − (θ₁₂−h)ρ] − z₁⁻¹[(2θ₁₂−3h) − (θ₁₂−h)z₃] } · ⟨Ω|Ψ_B⟩.**

**Rescaling the monodromy.** Under T → T/g, the combination changes only through the single number g(θ₁−h). The clean coefficient c in T̃₂₁ + c·T̃₁₁ is then (with h = 1):

| g | c |
|---|---|
| none | 2t−3 |
| ν₁ | (2t−3)/(t+1) |
| ν₂ | singular |
| ν₃ | −(2t−3)/(t−1) |
| f-normalized | (2t−3)/(t−1) |
| ∏_s(u−θ_s−h/2) | **−4/3** |

So rescaling each Lax by its half-shifted factor (u−θ_s−h/2) makes the combination t-free: T̃₂₁ − (4/3)T̃₁₁. Among products of site factors, the half shift on site 2 is forced.

## 7. Numerical traps found today

- **Zeros with low accuracy print as exact 0.** At 120 digits, the [4,4,4] errors at ε ≤ 10⁻⁵ had no correct digits. Always report the Precision. Log10 of such a zero gives Indeterminate.
- **Already-clean insertions.** Their k is O(ε) noise. A relative null-space cut treats it as a constraint and discards the insertion; use an absolute cut (`kTol`).
- **List-minus-vector threading.** Subtracting a vector from a list of vectors of matching length threads element-wise and produces a fake "spread".
- **x₀ at small ε.** `Inverse[P_phys]` is badly conditioned there; use the factorized inverse.
- **Nonzero-entry lookup.** `Position` with a pattern test returned 0, which silently set m₀ = 1. Use `Ordering[Abs[h], -1]`.
- **Symbolic substitution per state.** Substituting symbolic 64×64 operators for each state is far too slow. Substitute once and cache.
- **Shared kernel.** The new-conventions notebooks share one kernel and clear it on their first cell, so rerun a notebook before resuming work in it.

## 8. Open

1. **Closed forms for the deep states** of [3,1,0] and [4,2,0]. They presumably need higher orders (ψ₂, …), and the multiplicity copies above them bring in algebraic data.
2. **A rule for per-copy clean combinations**, for example via the copy's eigenvalues.
3. **Derive the dual QQ factor G** and its +h shift, rather than testing them. Check the single-Q form for [3,1,0].
4. **Validate the 6×6 ⟨Ω|Ψ⟩ formula** beyond the adjoint.
5. **The half-shift normalization:** does it make the [3,1,0]/[4,2,0] clean conditions t-free too?
6. **A second chain point** for [3,1,0] and [4,2,0].

## Files

- `Minimal-ingredients-EH.wb` (23 cells): adjoint. The chain, the contamination, v*, and β at P1 and P2.
- `Minimal-ingredients-310-EH.wb` (33 cells): [3,1,0]. Sections 1–10, plus Appendix A (the localization scan).
- `Blindspot-420-check-EH.wb` (6 cells): the correction of the [4,2,0] blind spot.
- `Minimal-ingredients-420-EH.wb` (28 cells, 250 digits): [4,2,0]. Sections 1–6.
- `Adjoint-clean-analytic.wb` (11 cells): the exact derivation, the Q-form and the rescalings.
- `Minimal-ingredients-{320,410,430,510,540,520,530}-EH.wb`: Λ-generic, site 1 only. Each contains the localization scan, the f fit (checked at P2; not in 320), the structural test and the fair test. 510 onward add the 10⁻⁶ refit (6b); 520 and 530 add the weak-sector comparison (6c).
- `L2-420-localization.md` has a correction note. `HANDOFF.md` has the summary.
