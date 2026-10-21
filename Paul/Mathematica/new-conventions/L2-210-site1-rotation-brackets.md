# The generic site-1 rotation on the SoV vacuum as ladder-1 Q-brackets: a 2-parameter family of genuine rotations with a four-bracket-determinant overlap (2026-10-10)

Notebook: `L2-210-site1-rotation-brackets.wb` (36 cells, saved; **no native outputs** — VS Code did not execute cells
of the new file, every code cell was run in K1 from its saved source and the printed results are in markdown cells
under the code cells). Spec: `L2-210-site1-rotation-brackets-spec.md`. Built by a sonnet wolfbook-builder, planned
and reviewed by the main session. Setting: gl(3) adjoint, L = 2, companion frame (x₀ = e₁ᵀ⊗e₁ᵀ), new conventions,
h = 1/2, θ = (2/7, −3/5) (t = θ₁₂/h = 62/35), z = (3/2, 5/4, 2/3) (loci re-checked at z′ = (3/2, 5/4, 3/7) and
at two (h,θ) points with h = 1/3). Builds on `Two-rotations-det-overlap.md` §3 (Part A). All statements exact
unless tagged [numeric] (64-state checks at 60 digits).

## 0. Paul's question and the answer

Start from the SoV vacuum, put the most generic rotation on site 1, Φ(g) = e₁ᵀρ(g)⊗e₁ᵀ, and fix the three
parameters so that ⟨Φ(g)|Ψ_B⟩/⟨x₀|Ψ_B⟩ is written only through ladder-1 brackets ⟨Q^jF⟩₁ = Σ_n F_nQ^j(θ₁−nh)
with explicit measures F (no determinant required); pull the θ₂-ladder dependence (λ₂ = Q^k(θ₂)/Q^k(θ₂−h)) out.

**Answer.** Write g = exp(αE₁₂)exp(βE₁₃)exp(γE₂₃) (first row of the 3×3 unipotent: b = (1, α, β+αγ)). On the
quadric

    2(t−2)·(β + αγ) = (t−3)·α²       i.e.   2(t−2) b₁b₃ = (t−3) b₂²,   t = θ₁₂/h   (h-independent coefficient),

a 2-parameter family of genuine rotations (αβγ ≠ 0), the overlap is exactly

    ⟨Φ(g)|Ψ_B⟩/⟨x₀|Ψ_B⟩ · Ψ₁ = Σ_{r∈{−1,0}} ( Bk₂[F⁰_r, G⁰_r] + λ₂ Bk₂[F¹_r, G¹_r] ) · Q³(θ₁ − r h),

    Bk₂[F,G] := ⟨Q³F⟩₁⟨Q²G⟩₁ − ⟨Q²F⟩₁⟨Q³G⟩₁,   Ψ₁ = Q²(θ₁−2h)Q³(θ₁−3h) − Q³(θ₁−2h)Q²(θ₁−3h) · Q³(θ₁−h) (up to sign),

with eight explicit 4-point measures F, G (rational in g; listed in the notebook at the point (α,β,γ) = (1, 11/16, 2))
— four 2×2 bracket determinants, two Q³ point values, and λ₂. [numeric: 10⁻⁵³ on all 64 states.] The λ₂¹ part
(N¹) has this form for EVERY rotation; the quadric is the condition on the λ₂⁰ part (N⁰). Off the quadric N⁰ is
not of this form (R2δ fails: the r = 0 column of the coefficient tensor is not a decomposable 2-form). At finite
twist the undone site-1 covector e₁ᵀρ(g L₀D₀) (site 2 = x₀^{(1)} = e₁ᵀρ(U₀), P_phys⁻¹ = L₀D₀U₀) reaches **all seven
adjoint weights with rank 2 on (1,1,1)** over three points of the quadric — the family is CG-complete at finite twist.

Every more restrictive form contains no genuine rotation: a single product det₂×⟨Q³H⟩ (R1) only for α = β = 0;
a single L=1-type 3×3 bracket determinant (R3) only for β = 0 or on a component at infinity of the Gauss chart;
the λ₂-pencil absorbed into the measures (decomposable for all λ₂) only for β = 0; the composite 3×3 with
λ₂-dependent third column likewise β = 0. So "sum of two bracket determinants per λ₂-power, with point measures
on the Q³ factor" is the sharpest bracket form that genuine rotations admit, and it costs one quadric.

## 1. The exact structure behind it

- **Monomial table.** ψ_x·Ψ₁ = M_x := det₂[v(n′),v(n)]·v_k(r) for the eight site-1 SoV triples x = (n,n′;r):
  (−3,−2;−1), (−3,−1;−1), (−2,−1;−1), (−3,−2;0), (−3,−1;0), (−3,0;0), (−2,−1;0), (−2,0;0) (rows {1,2}/k=1,
  {2,3}/k=3, {1,3}/k=2 all work). Hence N⁰ = Σφ⁰_xM_x, N¹ = Σφ¹_xM_x with φ⁰ = c·Kvac, φ¹ = c·K20 bilinear in
  (a,b); coefficient tensors C⁰ (8 populated slots), C¹ (4 slots: (−3,−2;−1), (−3,−2;0), (−3,−1;0), (−2,−1;0)).
  The minor (−1,0) never occurs (v(θ₁) ∥ v(θ₁−h)).
- **Loci (Gauss chart / whole cone b = a×u):** Pf(C⁰_{·,0}) = −(4/7533)β²(16β+16αγ−43α²) = (4/7533)a₁²(43b₂²−16b₁b₃)
  at t = 62/35; C⁰_{·,−1}∧C⁰_{·,0} = −(32/7533)αβ(β+αγ); C⁰_{·,0}∧C¹_{·,0} = −(280/7533)αβ²(β+αγ);
  Pf(C¹_{·,r}) = 0 and C¹_{·,−1}∧C⁰_{·,0} = 0 identically. 43/16 = c(t) = (t−3)/(2(t−2)) (fitted at six (h,θ)
  points, confirmed at two with h = 1/3: t = 29/15 → 8, t = 31/35 → 37/39).
- **Dictionary.** Bk₂[F,G] = Σ_{n<n′}(F∧G)_{nn′}·minor(n,n′); the L=1-type 3×3 det[⟨Q³F_i⟩,⟨Q²F_i⟩,⟨Q³H_i⟩] =
  Σ_i ±H_i⊗(F_j∧F_k); Pf(A) = A₁₂A₃₄ − A₁₃A₂₄ + A₁₄A₂₃ in the minor order (12),(13),(23),(14),(24),(34).
- **z-independence:** Kvac, K20 are identical exact rational matrices at z and z′; all loci identical.
- **Natural measures** (det₆-entry type F_n = wgt₁(n)·w^{−x_n/h}(p+qx_n)): on the genuine point each column's
  plane contains natural vectors (column −1: two, p = 0, w from (−32±√7990)/54; column 0: three, roots of a cubic),
  so F and G can both be natural but with DIFFERENT w; a common w is impossible (exact). Column −1 is also a product
  of two 2-point bras (supports {−2,−1},{−3,−1},{−3,−2}); column 0 has one 2-point vector ({−3,−2}).
- **Reach table (exact, finite twist, undone covector e₁ᵀρ(gL₀D₀); L₀ = {{1,0,0},{−23/12,1,0},{5/6,−2/3,1}},
  D₀ = diag(1,−1/4,35/72), U₀ = {{1,1,1},{0,1,10/3},{0,0,1}} at z):**

| locus | weights reached | rank on (1,1,1) |
|---|---|---|
| quadric (R2δ), genuine | all 7 | 2 |
| β = 0 (R3 / pencil) | 5 (not (0,1,2), (0,2,1)) | 1 |
| α = β = 0 (R1) | (2,1,0), (2,0,1) | 0 |
| generic rotation (control) | all 7 | 2 |

## 2. Cost

Per state: the dual Q's on the four ladder-1 points (Q², Q³) and Q³ at θ₁−h, θ₁; λ₂ from the two ladder-2 points;
eight fixed 4-point measures. Four 2×2 bracket determinants, no 6×6, no eigenvalue data. The count is set by the
ladder length (λ₁+1 points) and the rank, not by dim V.

## 3. Next (for Paul)

1. **Untwist the quadric family.** Site 2 is the SoV vacuum x₀^{(1)} = e₁ᵀρ(U₀(z)); the undone covector is
   e₁ᵀρ(g L₀(z) D₀(z)) with g on the quadric. In the regime z = 1+εφ the Gauss factors L₀, D₀ of P_phys⁻¹ scale
   with ε; the question is whether (α,β,γ) can be chosen ε-dependent on the quadric so that the undone covector is
   a generic unipotent at ε = 0 (two free parameters + the quadric is t-dependent, not ε-dependent — promising).
   Then K/S for the six branches and CG at leading order, as in `L2-210-untwist-CG.wb` Tasks 1–2.
2. Meaning of the quadric n₁₃ = (t−3)/(2(t−2)) n₁₂²: compare with the §5 product point (y₀ = (t−1)/(t−2)) and the
   singlet constant κ = 2(t−3)/(3(t−2)); is the quadric the image of the det₆×det₂ family's site-1 factor?
3. Site-2 analogue (rotation on site 2 with x₀ on site 1) and both sites rotated (two quadrics?).
4. Tidy: cell 21 ends on a missing check key; cell 25's `pencilCheck` line is wrong (symbol clash `l`, redone in
   the Task 5 cell); the Task 6 driver wipes Global` — run it last. Rerun in VS Code to store outputs.
