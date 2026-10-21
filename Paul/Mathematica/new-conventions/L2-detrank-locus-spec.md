# Spec (Job D): the det₆-rank locus — which covectors have the vacuum's low block ranks, and do any of them read every irrep? (2026-10-11)

Delegated to a wolfbook-builder subagent (sonnet); planned and reviewed by the main session. **New notebook:**
`Paul/Mathematica/new-conventions/L2-detrank-locus.wb` (rep-generic with `Λ`; [2,1,0] then [3,1,0]; build/run as the
previous stages). Stage B: the blocks A^{(c,c′)}(Φ) are LINEAR in the covector Φ ∈ V*⊗V* (A = Σ_ν Φ_ν A_ν with the K
tensor's rows A_ν), and the vacuum Ωc[w]⊗Ωc[w] (= det₆) has block ranks below the generic ones: adjoint (5,3,3,2) vs
(5,3,3,3) (block (1,1): 2 vs 3), [3,1,0] (8,6,6,6) vs (9,6,6,6) (block (0,0): 8 vs 9), and the full matrix φ_{x₁x₂} rank 12
vs 15 for [3,1,0]. The locus L_k := {Φ : rank A^{(k)}(Φ) ≤ r_vac} is a determinantal variety in P(V*⊗V*) containing the
2-parameter vacuum family (and e₁ᵀ⊗anything? check). Paul's question in this form: is L_k bigger than the vacua, and
does it contain covectors that read every irrep of Λ⊗Λ (which the vacua do not: ⟨Ω_w⊗Ω_w|h_ν⟩ = 0 for ν ≠ top)? Such
covectors would be "det₆-cheap and CG-useful".

Read first: HANDOFF top section; Stage B notebook cells 20–25, 37–41, 47–49; Stage D cell 21 (reach machinery; HW kets).
Dumps: `ref-two-site-all.txt`, `ref-two-vacua-untwist-all.txt`.

## 1. Tasks (per rep; each its own cells + markdown result; save after each)

**Task 0 — setup.** Stage B setup, K tensor, blocks as linear maps Φ ↦ A^{(cc′)}(Φ) (store the d₂ matrices A_ν^{(cc′)});
HW kets h_ν of Λ⊗Λ (Stage D); the reach functional "irreps read by Φ" := {ν : ⟨Φ|h_ν⟩ ≠ 0} (exact; for multiplicity-2
sectors the rank of the 1×2 row — note a single covector gives at most one functional per sector; CG-completeness for a
FAMILY needs rank = multiplicity over the family). Regression: vacuum ranks, generic ranks, vacuum reads only the top.

**Task 1 — tangent space of L_k at the vacuum (exact).** For the block k with the vacuum rank drop (adjoint (1,1): r = 2;
[3,1,0] (0,0): r = 8; also the full φ-matrix rank 12 locus for [3,1,0]): the (r+1)×(r+1) minors of A^{(k)}(Φ) vanish on
L_k; their differentials at Φ₀ = Ωc[w]⊗Ωc[w] (rational w) span the conormal space; dim T_{Φ₀}L_k = d₂ − rank(dMinors)
(the Jacobian of all (r+1)-minors, or equivalently of the map Φ ↦ A^{(k)}(Φ) composed with the tangent of the
rank-≤-r variety: T = {δΦ : U_⊥ᵀ A^{(k)}(δΦ) V_⊥ = 0} with U_⊥, V_⊥ bases of the left/right null spaces of A^{(k)}(Φ₀)
— use this linear form, it is exact and cheap: dim T = d₂ − rank of the linear map δΦ ↦ U_⊥ᵀA(δΦ)V_⊥). Compare with the
dimension of the vacuum family (2, plus the scalar) and with the tangent of {products φ₁⊗φ₂} (dim 2d₁ − 1). Repeat at a
second w and at a product-with-e₁ point (Ωc[w]⊗e₁ᵀ: which blocks drop there?). Report: dim T_{Φ₀}L_k for each k; the
tangent directions that are NOT tangent to the vacuum family (a basis, as covectors; are they products? R-dressed
objects? sums of vacua? — test membership of δΦ in span{Ωc[w′]⊗Ωc[w′]} (all w′), in {products}, in the R-dressed site-2
family (e₁ᵀ⊗c)R, in Ω₀-type rotated vacua).

**Task 2 — the full locus near the vacuum, numerically (prec 40) then exact.** From Φ₀ + small random δΦ in T, Newton
iterations on the minors (or alternating projection: project to rank ≤ r in block k, then back to… — Φ is the variable,
A is linear in it: use Gauss–Newton on the residual U_⊥ᵀA(Φ)V_⊥ with U_⊥, V_⊥ updated each step) to land on L_k; collect
≥ 20 points; for each: (i) the other blocks' ranks, (ii) Σr and ΣP, (iii) the irreps read (⟨Φ|h_ν⟩ ≠ 0, numerically with a
10⁻²⁰ threshold), (iv) is it a product covector (rank of the d₁×d₁ matrix Φ), (v) rationalisation attempt. Local dimension
of L_k around Φ₀ from the Newton points (numerical rank of the collected tangent vectors) vs dim T (if equal, L_k is
smooth there of that dimension). Also the intersection of L_k with the product covectors {φ₁⊗φ₂}: restrict Φ = φ₁⊗φ₂
(2d₁ unknowns), minors of A^{(k)} become bilinear-degree polynomials; Newton on this smaller system from random starts
(≥ 50) and classify the solutions found (vacua? Ωc[w]⊗e₁-type? new?) with their irreps read. For [3,1,0] the same for
the rank-12 locus of φ.

**Task 3 — the question answered.** Among all points of L_k found (general and product), the minimum ΣP and the set of
irreps read; whether any point reads every irrep (or at least more than the top); whether a 1-parameter family of such
points exists (then its tangent gives the multiplicity-2 functionals). State plainly: "the det₆-rank locus does / does not
contain CG-useful covectors beyond the vacua (within the search)". If a new CG-reading low-rank covector appears, give it
exactly if rationalisable, its blocks, and its undone/irrep table.

**Task 4 — markdown summary** `### det₆-rank locus — result (2026-10-11)` per rep: tangent dimensions, non-vacuum
tangent directions and what they are, the Newton exploration (counts, ranks, irreps read, product/non-product), the
product-restricted solutions, the answer, timings. Save; no HANDOFF.

## 2. Rules and gotchas

As before (K1 shared; first cell wipes; no `*)` in comments; no `D`,`E`,`N`,`I`,`C`,`K` symbols; time before loops;
`TimeConstrained`; save after each task; no other notebooks, no `paper/*.tex`, no git, no HANDOFF). Exact ranks on
rationals; numerical ranks with explicit tolerances stated. Projective scale: fix Φ·e₁⊗e₁ = 1 (the x₀ component) when
counting dimensions. The K tensor rows: A_ν^{(cc′)} = block of the ν-th unit covector; A(Φ) = Σ_ν Φ_ν A_ν (verify against
`tensorOf`-type code on one Φ).
