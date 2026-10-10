# Direction 3: a short sum of determinants

## Motivation

- **L=1 precedent.** A single det N reaches only companion eigen-bras. A generic rotation is exactly a sum of **two** det N's from one companion twist (α, β₁, β₂): Rotation-e8-SoV §4d. Single-determinant forms are therefore not expected for generic objects.
- **L=2 hint.** The site-1 rotation's SoV expansion splits into two pieces (`../L2-rotation-sov.wb` §3–4):
  - 8 states with site 2 at the vacuum;
  - 4 states with site 2 shifted to {2|20}, weighted by λ₂(B) = v(θ₂)/v(θ₂−h).

  The naive mixed det₆ also came out as A⊗e₁ + A′⊗e₂, i.e. two terms.

## Ansätze

1. **Two per-ladder dets with fixed natural data.** T_B ≈ α₁ D[B; F_w, F_a] + α₂ D[B; F_w, F_b], where:
   - D is the per-ladder det₆/Ψ_Q;
   - F_a, F_b are natural ladder-2 data, e.g. F₀ and its {2|20} partner (F₂ shifted by one point);
   - the α's are fitted linearly.

   This is a linear least-squares problem over a small dictionary, with no trivial branch, because the α's multiply nonzero dets.
2. **A dictionary of natural determinants.**
   - **Ladder data:** F_w (L=2 and L=1 versions), F₀, delta-type data supported on 1–2 points, and E-type exponentials.
   - **Dictionary:** all per-ladder det₆ formed from pairs (ladder-1 data, ladder-2 data).
   - **Fit:** T_B ≈ Σ_k α_k D_k(B), using sparse least squares (orthogonal matching pursuit / L1).
   - Look for a 2–3 term representation that holds at w and w′ with the **same structure** and α's that are simple functions of w.
3. **The X-relaxation from Direction 1.** For fixed F₁⁽¹⁾, F₁⁽²⁾ the det is bilinear in (F₂⁽¹⁾, F₂⁽²⁾). If the general 4×4 matrix X replacing F₂⁽¹⁾⊗F₂⁽²⁾ solves the system, then rank(X) is the number of determinants needed.

## Safeguards

- Each α_k multiplies a fixed nonzero det, so there is no λ → 0 branch. Still check the conditioning of the dictionary matrix (64 × #atoms).
- Hold-out test: fit the α's on 48 states and predict the other 16. Then fit at w and predict at w′ with α(w′) from the conjectured closed form.

## What a positive result would look like

⟨Ω_w⊗e₁ᵀ|Ψ⟩/⟨x₀|Ψ⟩ = D[F_w, F₀] + λ₂-type correction, where the second determinant implements the {2|20} shift. This is the Q-form counterpart of the exact SoV split.

## Cost

Building the dictionary: about 20–100 atoms × 64 states × one 6×6 det each, which is seconds to a minute. The linear fits are trivial.

## Outcome (2026-10-02): negative for short sums

The code is in `../L2-rotation-sov.wb` §9, cells 54–57.

- **Natural dictionary** (56 nonzero atoms: F_w, F_w(L=1) and delta-type data on the ladders).
  - The span has rank 55 < 64; the full least-squares residual is 9.6·10⁻⁴.
  - With ladder 1 fixed to F_w, the residual is 0.70.
  - The greedy sparse path is slow (0.67 → 0.44 at 5 atoms).
- **All per-ladder dets.** The polarized atoms span all of C⁶⁴, so a long enough sum always exists. The real question is the minimal number of terms.
- **Two free per-ladder dets** (32 numbers, λ = 1). Stuck at about 0.107, with the two terms each about 6‖T‖ and cancelling. Not a solution.
- **Conclusion.** No compact sum. Within the 6×6 bracket family, Directions 1–3 rule out a compact determinant form for Ω_w⊗e₁ᵀ. Direction 4 (closed-form SoV sum) is the remaining route, unless a different determinant shape is used.
