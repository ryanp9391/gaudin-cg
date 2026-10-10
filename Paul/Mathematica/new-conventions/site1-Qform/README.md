# Q-form of the site-1 objects at L=2: research directions

**Goal.** Find a closed Q-function (determinant-type) formula at L=2 for overlaps of states rotated on site 1 only, with site 2 left at the SoV vacuum. There are two targets:
- the site-1 vacuum, ⟨Ω_w⊗e₁ᵀ|Ψ_B⟩/⟨x₀|Ψ_B⟩;
- the full 3-parameter site-1 rotation, ⟨(e₁ᵀρ(g))⊗e₁ᵀ|Ψ_B⟩/⟨x₀|Ψ_B⟩, with g ↔ (w, y₁) via the L=1 resolvent decomposition e₁ᵀρ(g) ∝ Ω_w ρ((y₁−G(z)ᵀ)⁻¹).

The global counterparts (rotation on both sites) are solved:
- ⟨Ω_w⊗Ω_w|Ψ⟩ is the usual 6×6;
- the global rotation is that 6×6 times the su(2)-like 4×4.

See `../Rotations-SoV-Q.{tex,pdf,md}` for those results, and for the retracted "single 6×6 with a numeric ladder-2 F⁽²⁾" claim.

**What is already known** (`../L2-rotation-sov.wb`):
- **Exact SoV representation.** The site-1 rotation lives on 12 SoV states: site 2 at the vacuum, or at {2|20}. Its overlap is Σ_ν c^{L=1}_ν(g) W_ν(B), with a rational, g-independent map K (§3–4).
- **Vacuum data.** F₀ = (F₁ = δ_{θ_s−2h}, F₂ = δ_{θ_s}) on both ladders reproduces x₀ (§5).
- **Per-ladder brackets** reproduce Ω⊗Ω and x₀ (§5).
- **Ruled out:** a per-ladder 6×6 with ladder 1 = the w-vacuum Q's (L=2 or L=1 version) and arbitrary ladder-2 data. The best residual is 0.766 (§6b).

**Mandatory safeguard for every direction.** Fix the overall scale. Never fit |det − λT| with λ free, and report ‖det‖/‖T‖ and λ for every claimed solution. The gauge F₂ → F₂ + cF₁ can make a row duplicate another, giving det ≡ 0 (memory: `feedback_det_fit_trivial_branch`).

| File | Direction | Status |
|---|---|---|
| `1-free-ladder-data.md` | Let **both** ladders carry free bra data in the per-ladder 6×6 | **done: negative** (best residual 0.273; X-relaxation 0.269) |
| `2-polynomial-insertions.md` | Polynomial insertions π(x) in the brackets that kill ladder-2 poles | **done: negative** (joint w/w″ minimum 0.121) |
| `3-sum-of-dets.md` | Allow a short sum of determinants (as at L=1, where generic rotations need 2) | **done: negative** (natural dictionary rank 55, residual 1e-3; two free dets stuck at 0.107) |
| `4-SoV-sum-closed-form.md` | Accept the 12-state SoV sum; get K in closed form in t = θ₁₂/h | planned |
| `5-bracket-decomposition.md` | Read the exact 12-state sum term by term in the 2-plane variables (α,β,λ₁,λ₂); classify local ladder-1 bracket dets | **done 2026-10-09: negative for single dets** (local 3×3 family is a 4-dim cone; groups not in it, exact); compact polynomial form of the sum obtained |
| `6-two-rotations-2603-trick.md` | Two distinct rotations v₁⊗v₂ with overlap = det₆ × ladder-1 contour det (the 2603.26617 App. G trick) | **done 2026-10-09: POSITIVE** — v₂=Ω_{(w₁,y*)}, v₁=v₂ρ((w₂−Gᵀ)⁻¹), y*=w₂(θ₁₂−h)/(θ₁₂−2h); 10⁻⁵⁰ on 64 states |
| `7-undo-site2-rotation-span.md` | Undo the site-2 rotation of the two-rotation family, count site-1 weight-space spans; pole-set chain of product points | **done 2026-10-09 (late): POSITIVE** — adjoint 8/8; [λ,1,0]: product points ↔ pole sets {0..k}, k=λ₂+1..λ₁, y_k=w₂(θ₁₂−h)/(θ₁₂−kh), det₆×det₂^{(0..k)}/W(λ₁,λ₁+1); the k=2 point determines all HW states ([2,1,0],[3,1,0],[4,1,0]); [3,2,0] (d=1) only k=3 → 9/10; τ_λ(θ₁−λ₁h) = one Wronskian/minor. Full notes: `../Undo-site2-rotation-notes.{md,tex,pdf}` |

**Numeric point.** h=1/2, θ=(2/7,−3/5), z=(3/2,5/4,2/3), w=(5/2,7/3,3/5), with a second twist w′=(7/4,−2/5,9/7). Any positive result must then be re-checked at a second (h,θ).
