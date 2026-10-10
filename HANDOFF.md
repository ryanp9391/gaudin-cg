# Hand-off

Where the last session left off. Read this first when picking the project back up; it's not a
history log (that's `git log`) — only what's still open or in-flight. Overwrite
freely each session rather than appending.

## >>> NEWEST (2026-10-09, late): undoing the site-2 rotation of the two-rotation family — adjoint complete (8/8, singlet); [3,1,0] and [3,2,0]: all HW states reached, 9/10 determined by the det-form point, [3,1,0]'s point A gives 10/10 <<<
**FULL NOTES (read first):** `new-conventions/Undo-site2-rotation-notes.{md,tex,pdf}` (23 pp., sonnet write-up of the whole thread, every result with notebook+cell, conjectures tagged). **Spec + running log:** `new-conventions/site1-Qform/7-undo-site2-rotation-span.md` §1–§16. **Notebooks:** `L2-rotation-sov.wb` §12 (cells 113–123), `L2-310-two-rotations.wb` cells 24–70, `L2-320-two-rotations.wb` (41 cells), `L2-410-two-rotations.wb` (26 cells) — all saved, not committed; sonnet subagents from the spec. Paul's question: for ⟨v₁⊗v₂| with v₂ = ⟨0|R, does ⟨v₁|R⁻¹ projected to the site-1 weight spaces (all coefficients of the free parameters) give a complete basis? The 2026-09-23 families stopped at μ₃ ≤ λ₂ (singlet missed).
- **Undo is exact at finite twist:** in the diagonal frame v₂ = e₁ᵀρ(n), n ∈ N₊ (the vacuum Gauss factors r_k(z,(w₁,y))), and ρ(n)⊗ρ(n) fixes all HW kets, so ⟨v₁⊗v₂|HW⟩ = ⟨v₁ρ(n)⁻¹⊗e₁ᵀ|HW⟩. No untwisting needed for the count.
- **Mechanism (3×3, exact, rep-independent):** φ = v₁ρ(n)⁻¹ = δ^Λ e₁ᵀρ(n′), n′ = Δ⁻¹nΔn⁻¹ = [[1,a′,b′],[0,1,c′],[0,0,1]], Δ = diag(1/(ω−z_k)), ω = resolvent point w₂(y). **b′ = 0 ⟺ ω = y**, and b′/(a′c′) = (y−ω)/(y−z₂). The old families' p₃ relation is exactly b′ = 0 (Δ = diag z); the new family has ω = y(θ₁₂−2h)/(θ₁₂−h) ≠ y — the (t−1)/(t−2) mismatch between the det₆ twist and the ladder insertion is what opens the E₁₃ direction.
- **Adjoint:** per-weight ranks = dims (8/8), exact in y, unchanged under w₁/torus; control ω→y kills (0,1,2),(1,0,2) (old blind spot). HW count: 5 at fixed w₁ (one harmless cross-weight relation), **6 with the torus**; singlet reached. Untwisting leading orders ε⁰,ε¹,ε¹,ε²(rank 2),ε³,ε³,ε⁴. (0,1,2) component = 3b′(a′c′−b′)δ^Λ.
- **[3,1,0] (point-B family over six w₂, β/α re-solved per w₂):** rotation part alone spans 15/15; the actual dressed covector e₁ᵀ(α+βE₁₂)ρ(nΔn⁻¹) spans **12/15** — the fitted dressing cancels (0,1,3), (1,0,3) identically and one direction of (1,1,2) (consistent with the e₁ᵀE₁₂-orbit type). Still **10/10 HW states reached** with the torus ([4,2,2], [3,3,2] included; old 8/10); 9 fully determined, [4,2,2] only up to one combination of its two (1,1,2) components. Control (ω→y*, β=0): 6.
- **Kernel state at end:** K1 holds the [3,1,0] notebook's cells 2–7, 12, 24–33; `L2-rotation-sov.wb` live state wiped (rebuild cells 2–21 + the text of `oneSiteVac`/`w2Of` for §12; 63–111 need their own cells). Private kernels (K310, K2) cannot start in this environment — every notebook ends up on K1.
- **Follow-ups (Paul: "get cheap subagents to start doing these"; spec §9–§10):**
  - **[3,1,0] LW frame** (cells 34–42): exact mirror of the HW frame; [4,2,2] again rank 1 of 2 (at (2,1,1)); frames don't combine. No help.
  - **[3,1,0] point A** (same cells): site-1 covector (rotRank 6) undone exactly → **all 15 directions, all 10 HW states fully determined** (A alone or A∪B). Caveat: A's overlap is c_A·(1,τ_λ(θ₁−3h))⟨Ω|Ψ⟩ (T-polynomial × det₆, a Q-number via T–Q) but not a single ladder det₂.
  - **[3,2,0]** (`new-conventions/L2-320-two-rotations.wb`, 21 cells, saved, not committed): one isolated point only (point-B numbers; **no point A** since T₂₂, T₃₂ vanish), sites swapped (rotation on site 1 with y* = w₂(θ₁₂−h)/(θ₁₂−3h), dressed vacuum e₁ᵀ(1+bE₂₃)R on site 2, b fitted = (6125/5676)w₂/((w₂−z₃)(y*−z₃)) at P1). Undo exact; 12/15 directions; **10/10 HW reached, 9/10 determined, [4,4,2] one direction short**; control blind spot μ₁ ≤ 1 (dual of μ₃ ≥ 2).
- **A1/A2 (Paul's three-avenue plan, spec §11–§14):** **Point A of [3,1,0] IS a contour det₂:** the same ladder-1 det₂ as B but over the adjoint pole set {0,1,2} at y_A = w₂(t−h)/(t−2h), divided by the [3,1,0] minor W(3,4); its site-1 factor is a SECOND-order dressed rotation e₁ᵀ(1+βE₁₂+γ(E₁₃+½E₁₂E₂₃))ρ(n_A)ρ(Dy[w₂]) with closed-form β, γ (spec §13). **Explicit Q-dictionary:** τ_{(λ₁,λ₂)}(θ₁−3h) = C_λ·W(θ₁−(3−λ₁)h, θ₁−(4−λ₂)h)/W(θ₁−3h,θ₁−4h), C_λ explicit in h, t only; t₁, t₂ at any x by Cramer from the dual Baxter equation. So every "T-polynomial × det₆" overlap is explicit in Q's. [3,1,0]²: both points det₆ × det₂^{(P)}, A alone determines 10/10. **A2 negative:** [3,2,0] site-2 insertions give only the localizing span (19-op span, modular Gröbner) → [4,4,2] still one direction short. **A3 DONE — the pole-set chain (spec §14–§16):** product points of the θ₁−λ₁h recipe ↔ pole sets P_k = {0..k}, k = 2..λ₁, at y_k = w₂(t−h)/(t−kh), overlap det₆ × det₂^{(P_k)}[y_k]/W(λ₁,λ₁+1), site-1 dressing of E₁₃-degree λ₁−k. Exact on [3,1,0] (only {0..3}→B, {0,1,2}→A among all 11 subsets) and **[4,1,0]** (`new-conventions/L2-410-two-rotations.wb`, 26 cells, saved; three points k=2,3,4 at −63/8, −63/43, −21/26; 576 states end-to-end 6e-43; vacuum triple (θ−5h,θ−4h,θ−h), minor W(4,5); 14 HW states in 11 irreps). **The k=2 (adjoint pole set) point spans all 24 directions and determines all 11 irreps**; k=3: 21/24, 11/11; k=4 (=B): 16/24, 9/11. [3,2,0]: the k=2 curve collapses into the localizing span (T₂₂,T₃₂ ≡ 0 on ⟨Ω|), so only B exists there ([3,2,0] vacuum triple (θ−4h,θ−3h,θ−2h), minor W(3,4), end-to-end 1e-44). Rep-data pattern: triple (θ−(λ₁+1)h, θ−λ₁h, θ−λ₂h), minor W(λ₁,λ₁+1). NB the [4,1,0] notebook's cells ran through the kernel (controller not associated); outputs not stored — rerun in VS Code (~1 min).
- **Pattern so far:** the pure-rotation partner always spans everything (adjoint: the whole family is pure rotations → complete). Where the recipe produces a *dressed* vacuum ([3,1,0] B, [3,2,0]) the dressing removes exactly the extremal μ-directions and one direction in one dim-2 weight → one bottom irrep ([4,2,2] resp. [4,4,2]) undetermined. [3,1,0] has a second point (A) that fixes it; [3,2,0] does not.
- **Kernel at end of day:** K1 holds `L2-320-two-rotations.wb` (all cells live); `L2-rotation-sov.wb` and `L2-310-two-rotations.wb` live state wiped (files saved: 123 and 42 cells). Private kernels don't start on this machine.
- **NEXT (ask Paul):** (1) reps where the k=2 curve collapses ([3,2,0]-type, λ₂ ≥ 2): find the missing covectors (site-2 analogue was empty; try the dual/LW construction or products); (2) [4,2,0] (dim 27): which k survive? (3) derive the dressing polynomials (closed forms β, γ exist for [3,1,0]) and the general c_k(y) curve from the det₂ residue expansion with the τ-dictionary; (4) the Q-side reading of the CG coefficients: each y-coefficient of det₆·det₂^{(P_k)}[y] is one HW CG functional after the (R⊗R)Ψ→HW cleaning of `Clean-CG-combinations` — put the k=2 point through that pipeline for [3,1,0].

## >>> NEWEST (2026-10-09, evening): [3,1,0]² — the two-rotation recipe gives "rotated vacuum ⊗ rotated DRESSED vacuum", not two rotations <<<
**Notebook:** `new-conventions/L2-310-two-rotations.wb` (14 cells, exact rationals, sparse; saved, not committed). **Write-up:** `new-conventions/L2-310-two-rotations.md`. Paul's request: generalize the adjoint two-rotation result to [3,1,0].
- **Mechanism identical:** at θ₁−3h, (1,1),(2,2),(3,2) localize; all other shapes λ₁ ≤ 3 have rank 2 with ONE common extra site-2 direction ω₂ (a rotation, pencil inside the orbit); exact product locus in the 10-dim span = localizing span + two isolated points. Shapes λ₁ = 4 and the mixed span (+θ₁−2h shapes; 9-dim site-2 span, 29 s Solve) add nothing.
- **Point A** = the adjoint's c* verbatim on the column-shifted shapes (1,1),(2,1),(2,2),(3,1),(3,2); site 2 = rotated vacuum with the adjoint's y = w₂(θ₁₂−h)/(θ₁₂−2h); site 1 rotRank 6, unidentified.
- **Point B** (all shapes): site 2 = rotated vacuum with **y* = w₂(θ₁₂−h)/(θ₁₂−3h)** (three numeric points, exact); **site 1 = e₁ᵀ(α+βE₁₂)·R(z,(w₁,y*))·ρ((w₂−Gᵀ)⁻¹)** — the adjoint's v₁ with the vacuum covector replaced by vacuum + first descendant (β/α = 392/2335 at P1; 294/2255 at h=1/3; 5400/337687 at the third point; no closed form). Not a rotation (rotRank 5 = e₁ᵀE₁₂ orbit type; exact). For the adjoint the same dressing IS a rotation (E₁₂² kills e₁ᵀ), which explains the adjoint result.
- **New tool:** exact linear rotation test for any rep: rank{φ, φE_ij} = 4 ⇔ φ ∈ GL(3)·e₁ᵀ (generic 9; e₁ᵀE₁₂-type 5; e₁ᵀ(E₁₂+E₂₃)-type 6; e₁ᵀE₁₃-type 7).
- **Q-side DONE (cells 15–23, Opus subagent from `L2-310-Q-side-spec.md`, reviewed):** 225 states, all Q's (deg q^j = μ_j, deg q_i = 6−μ_i); vacuum det₆ with poles n = 0..3 and Ψ_Q = Π_s[Q²(θ_s−3h)Q³(θ_s−4h) − Q³(θ_s−3h)Q²(θ_s−4h)]Q³(θ_s−h) (5e-45); ladder-1 det₂ with minor₁ at (θ₁−3h, θ₁−4h) lies in span{1, τ_λ(θ₁−3h)} with an exact box-graded curve c_λ(y) = a_λ y^{6−|λ|}; it hits c_B only at y = −63/43 = y* (exact); END-TO-END: ⟨v₁⊗v₂|Ψ⟩/⟨x₀|Ψ⟩ ∝ det₆(w₁,w₂) × det₂^{(1)}(y*) to 5e-45 (wrong y: 77, 322). Point A is not on the curve.
- **Not done:** Q-side at the other numeric points; closed forms of a_λ and β/α; point A's site 1.
- **Side effect:** the new notebook landed on the shared kernel K1 and its `ClearAll` wiped the live state of `L2-rotation-sov.wb` (file intact; rebuild cells 2–21, 24–27, 60–68, 86–106 if needed). The `y*)`-in-a-comment gotcha struck twice; write y-star in comments.
- **NEXT (ask Paul):** is "rotated dressed vacuum ⊗ rotated vacuum" useful (its overlap is a T-polynomial × ⟨Ω|Ψ⟩, hence a Q-expression)? If so: Q-side for [3,1,0]; [3,2,0]/[4,2,0] with the same notebook (change Λ; localization data in the L2-*-localization.md files); closed form of β/α.

## >>> (2026-10-09, later): beyond the vacuum — the two-rotation family is "any rotation on site 2, one fixed resolvent on site 1" (3-dim); no κ-family; no iteration <<<
**Notebook:** `new-conventions/L2-rotation-sov.wb` §11 (cells 107–112, saved, not committed). **Write-up:** `new-conventions/Two-rotations-det-overlap.md` §5.8 (+ summary item 4, §6 updated); short form appended to `site1-Qform/6-two-rotations-2603-trick.md`. Paul's question: can the §10 construction be generalized to states beyond the vacuum / more general rotations on each site?
- **Torus = character (exact).** ρ(D)⊗ρ(D) (ket-twist Cartan torus) multiplies Ψ_B by Π d_i^{μ_i}; for a resolvent this is the su(2)-like det₄. Every product state with a det overlap comes with its torus orbit.
- **w₁ is a torus direction, y a torus invariant** (exact Jacobian ranks 3 at fixed y / 4 with y free / 3 torus alone, three generic points; Ω_{(w₁,y)}/Ω_{(w₁′,y)} is a torus character). So the (w₁,y) family is (torus, y), and its torus orbit is V = {(g ρ((w₂(y(g))−Gᵀ)⁻¹), g) : g any rotation of x₀}: **every site-2 rotation g has exactly one site-1 partner**, overlap ∝ χ_μ(D) det₆N[Q_B, Q_vac,(w₁,w₂(y(g)))]/Ψ_Q × det₂^{(1)}[Q_B; y(g)]/minor₁. Checked for a random Gauss-decomposed g (y(g) unique from exact Solve at two w₁'s; 10⁻⁵¹ on 64 states) and for the family at fresh w₁ = −7/3, 9/4 (10⁻⁵¹, 10⁻⁵⁰).
- **Negatives:** no continuous κ-family for the adjoint even with the torus (the first version of the §11 intro claimed one — it was wrong, caused by a Jacobian evaluated at the degenerate point y = z₂; rewritten). The §10 step does **not** iterate: on v₁⊗v₂ the site-1 T's at θ₁−2h have 8×8 ranks {1,4,3,8,8}, full site-2 column space, product locus = seed only.
- **Gotchas hit:** Ω_{(w₁,y)} degenerates at y = z_i (0/0 everywhere); the string `y*)` inside a Wolfram comment closes the comment.
- **NEXT (discuss with Paul):** closed form of the invariant y(g); the site-2 analogue (insertions at t₂(θ₂+h), t₂(θ₂+h)t₁(θ₂−h) + non-localizing partners) on ⟨Ω⊗Ω| and on V; joint site-1 × site-2 span; general reps/L>2 (larger localizing set) followed by the torus analysis. Part A says (g₁, x₀) is not a single natural det, so a 6-dim single-det family is unlikely; sums of dets are the fallback.

## >>> (2026-10-09, later): POSITIVE — two distinct rotations with a det₆ × ladder-1 contour-det overlap (2603.26617 App. G trick for the gl(3) adjoint) <<<
**Notebook:** `new-conventions/L2-rotation-sov.wb` §10 (cells 86–103, saved, not committed). **FULL WRITE-UP + LOG (read first):** `new-conventions/Two-rotations-det-overlap.md`; short form `new-conventions/site1-Qform/6-two-rotations-2603-trick.md`.
- **Result.** With y* = w₂(θ₁₂−h)/(θ₁₂−2h): v₂ = Ω_{(w₁,y*)}, v₁ = Ω_{(w₁,y*)}ρ((w₂−G(z)ᵀ)⁻¹) (= Ω_{(w₂,y*)}ρ((w₁−Gᵀ)⁻¹)), both exact rotations, and ⟨v₁⊗v₂|Ψ_B⟩/⟨x₀|Ψ_B⟩ ∝ det₆[Q_B,Q_vac,(w₁,w₂)]/Ψ_Q × ladder-1 su(2)-like det₂[Q_B; y*]/minor₁ (cell 68's det with E = y*^{−x/h}), state-independent to 10⁻⁵⁰ on 64 states (cell 102). y* checked at five (h,θ,z,w) points; the det identity at P1 only.
- **Mechanism.** T₁₁,T₂₁(θ₁−2h) localize exactly on ⟨Ω_w⊗Ω_w|; T₁₀,T₂₀ add one extra site-2 direction ω₂ = Ω_{(w₁,z₂)}, and Ω₁+κω₂ is a line inside the rotation cone. The product locus of ⟨Ω|(c₀+c₁T₁₀+c₂T₁₁+c₃T₂₀+c₄T₂₁) has one nontrivial point, z-independent, = ⟨Ω|𝒦_{y*} with 𝒦_y the operator with eigenvalue the ladder-1 contour det f_y ∈ span{1,τ₁₀,τ₁₁,τ₂₀,τ₂₁}.
- **Negatives:** no continuous family in the enlarged linear span (t₂(θ₁), site-2 localizers); Ω_w⊗Ω_{w′} is not a ladder-dependent det₆ and its SoV bra wavefunction does not factorize.
- **NEXT:** second numeric point for the det identity; w₃ / finite-y₂ dependence; site-2 analogue; general reps and L>2 (bigger localizing set → possibly a continuous κ-family).

## >>> (2026-10-09): site-1 rotation as a Q-determinant — term-by-term reading of the 12-state sum (Direction 5) <<<
**Notebook:** `new-conventions/L2-rotation-sov.wb` §8 (cells 59–78, saved, not committed). **Spec/outcome:** `new-conventions/site1-Qform/5-bracket-decomposition.md`. Paul's idea: don't fit det ansätze; write the exact 12-state sum as a sum over wavefunctions, fix the site-2 factor, and read the ladder-1 terms for bracket structure.
- **2-plane variables.** On ladder 1, v(θ₁)∥v(θ₁−h) and v(θ₁−3h), v(θ₁−2h), v(θ₁−h) span a 2-plane. With λ_s = v(θ_s)/v(θ_s−h) and v(θ₁−h) = αv(θ₁−3h)+βv(θ₁−2h), the 8 site-1 wavefunctions are {1,β,−α}×{1,λ₁,λ₁²} (minus λ₁²·1). **The overlap of any site-1 covector is an explicit short polynomial in (α,β,λ₁,λ₂)** (table in the spec and in cell 61; exact to 10⁻⁴⁷). It splits as Σφ⁽¹⁾ψ + λ₂Σφ⁽²⁾ψ with φ⁽ⁱ⁾ = c·K⁽ⁱ⁾.
- **Negative results (all exact or 60-digit):** neither group is a rotation (X(φ) rank 3); the L=1 det₃×det₂ formula is not a polynomial identity (fed L=2 Q's on ladder 1 it leaves the 8-dim local span, residual 0.11); the L=1-type 3×3 bracket det is site-1-local iff F₂ ∝ F₁ on {θ₁−2h, θ₁−h}, and that family's image is only a 4-dim cone (10 quadrics) which does not contain the groups (Gröbner membership); sums of two such dets cover the whole local space (uninformative).
- **Verdict:** no single ladder-1 bracket determinant of the natural type; the exact representation is the 12-state sum in the compact polynomial form.
- **Also negative (§9, same day):** two different rotations g₁⊗g₂ — the full-measure det₆ with ladder-dependent vacuum Q's fails (deviations 5–25), and the SoV bra wavefunction of Ω_w⊗Ω_{w′} does not factorize over sites (measure-free test M_{ww′}M_{w′w} ≠ M_{ww}M_{w′w′}, dev 1.6). So no 'sum over SoV integrals with a factorized bra' either.
- **Incident:** a symbol clash (`uS`, the tPoly interpolation points, reused in §8) crashed the laptop; fixed (uSym), memory note written.
- **NEXT (ideas, not started):** brackets involving v(θ₁+h) (outside the 2-plane); two-ladder dets with insertions producing λ₂ (locality analysis of the per-ladder 6×6 in these variables); closed-form K(t) (Direction 4) now reduces to the few rational numbers in the polynomial table. Discuss with Paul before another search.

## >>> (2026-10-02, late): rotations of the SoV vacuum as Q-determinants — write-up `Paul/Mathematica/new-conventions/Rotations-SoV-Q.{tex,pdf,md}` <<<
**Notebooks:** `new-conventions/Rotation-e8-SoV.wb` (L=1) and `new-conventions/L2-rotation-sov.wb` (L=2). Both are saved and not committed.
- **L=1.**
  - Every generic rotation of x0 is ⟨x0|M(w)M(y1,∞)⟩. Here ⟨x0|M(w) = Ω_w, and M(y1,∞) = ρ((y1−G(z)ᵀ)⁻¹), the resolvent.
  - Overlap = det3[Q_B, Q_vac,w]/Ψ_Q × an su(2)-like 2×2 in (Q³,Q²) with insertion y1^{−x/h}.
- **L=2, global rotation.** Overlap = det6 (the usual Ω_{z,w} formula) × a 4×4 su(2)-like det, for all 64 states to 10⁻⁵⁰.
- **L=2, rotation on site 1 only.**
  - It is an exact sum over 12 SoV states, with a rational K map.
  - **RETRACTED:** the claim that ⟨Ω_w⊗e1ᵀ|Ψ⟩ is a single per-ladder 6×6 (ladder 1 = w-vacuum Q's, ladder 2 = numeric F^(2)). Those fits converged to the trivial branch, λ≈0 with det≡0.
  - Redone with the scale fixed (L2-rotation-sov.wb §6b): the best residual is 0.766 (L=2 vacuum Q's on ladder 1) or 0.774 (L=1 vacuum Q's). So this ansatz is ruled out.
  - **Lesson:** fix the scale in det fits.
- **NEXT:** a different Q-form for the site-1 vacuum: modified ladder-1 data, mixed brackets, or a sum of dets. The fallback is the exact 12-state SoV sum. Discuss with Paul before searching.

## >>> FULL WRITE-UP of 2026-10-02: `Paul/Mathematica/new-conventions/Clean-CG-combinations.{tex,pdf,md}` <<<
This file logs the whole contamination / clean-combination thread, from the CG chain to the analytic adjoint formula. Read it first; the entries below are the per-step notes.

## >>> NEWEST (2026-10-02, latest): CLOSED FORMS of the clean coefficients ([2,1,0], [3,1,0], [4,1,0]) <<<

**Notebooks:** `Paul/Mathematica/new-conventions/Clean-coeffs-210-EH.wb` (control), `Clean-coeffs-310-EH.wb` and `Clean-coeffs-410-EH.wb`. All are saved and not committed. The write-up is §3.8 of `Clean-CG-combinations`.
- **No universal combination exists.** In [5,3,0] the k's have full rank (`Minimal-ingredients-530-EH.wb` §6d).
- **Method:**
  - h = 1, θ = (t, 0); f-free one-site covectors, top-killed.
  - Richardson extrapolation over ε₀, 2ε₀, 4ε₀ (ε₀ = 10⁻⁸) gives the exact ε → 0 leading data (error ≤ 10⁻²³).
  - Rational fits in t over 20 points, validated on 5 held-out points.
  - Multiplicity sectors: tested in the form A(t) + B(t)λ_c, with λ_c = t₂(u₀) eigenvalue.
- **Results:**
  - The control reproduces β(t) exactly.
  - **[3,1,0] and [4,1,0]: every contamination/signal ratio is in closed form.** Rational in t for multiplicity-1 states, and P(t) + Q(t)σ_c√D(t) in multiplicity sectors, i.e. linear in the copy's eigenvalue. For example, D = 16(1+3t²) for [4,3,1].
  - **Layer structure:** t-independent universal combinations kill everything except the bottom irrep. These are 6d₁₁ − d₂₂ for [3,1,0], and 16d₁₁ − d₂₂, 72d₁₁ − d₃₃ (layer 1) and v₂ = 72d₁₁ − 9d₂₂ + d₃₃ (layer 2) for [4,1,0].
  - The bottom irrep needs t-rational coefficients, with ρ = (2t−3)/3, (7t²−12)/20 and (11t³+61t−90)/90 for λ₁ = 2, 3, 4.
- **Mechanism (tested where non-trivial):** k̂_B lies in the span of the signals of the HW states above B.
- **√D coefficients identified (`Clean-coeffs-310-EH.wb` §4).** Per unit signal they are affine in τ_c = untwisted t₂(θ₁+h) eigenvalue, with simple coefficients. The [4,3,1] clean combination is v_c = 48[19t(t−1) − τ_c] d₁₁ − [4t(t−1)(3t−23) − (3t−11)τ_c] d₃₂, exact to 10⁻²². Not yet tested on [4,1,0] (its [6,3,1] data would need a fresh 18-min run).
- **Next:**
  - Prove or test the layer/annihilator mechanism more sharply. Find the minimal set of states above B.
  - Explain the discriminants and the integers 6, 16, 72, (72, −9, 1).
  - Extend to [4,2,0] and the λ₁ = 5 reps. These need multiplicity-3 sectors (λ² terms) and more fit points; [4,1,0] already took 18 min.

## >>> (2026-10-02, end of day): site-1 clean combinations for ALL λ₁ ≤ 5 reps ([3,2,0] … [5,4,0]) <<<

**Notebooks:** `Paul/Mathematica/new-conventions/Minimal-ingredients-320-EH.wb` (15 cells) and `Minimal-ingredients-410-EH.wb` (17 cells). They share Λ-generic code. Both are saved, ran successfully and are not committed. The write-up is Section 3.6 of `Clean-CG-combinations.{md,tex,pdf}`.
- **Pipeline:** site 1 + ⟨Ω|Ψ⟩ only.
  - The site-1 insertions come from an in-notebook localization scan.
  - The covectors are computed directly, ⟨Ω|T_λ(R⊗R)⁻¹, so no f_λ are needed.
  - The eigenvalues come from CBR of the scalar t-eigenvalues.
  - Fair test: fit at (10⁻⁵, φ, ω), test at (10⁻⁴, φ) and (10⁻⁵, φ′).
- **[3,2,0]²** (150 digits, about 3 s): all 9 non-top HW states recovered, O(ε) errors. Only [5,3,2] needs per-copy cleaning, mirroring [4,3,1] in [3,1,0].
- **[4,1,0]²** (200 digits, about 25 s): this was the first localization for this rep.
  - **Site-1 local:** T₍₁,₁₎, T₍₂,₂₎, T₍₃,₃₎, T₍₄,₃₎ at θ₁−4h. The old pattern missed T₍₂,₂₎.
  - **f's** (fitted at P1, checked at P2): f₁₁ = (θ₁₂−3h)(θ₁₂−5h); f₂₂ = (θ₁₂−2h)(θ₁₂−3h)(θ₁₂−4h)(θ₁₂−5h); f₃₃ = (θ₁₂−h)(θ₁₂−2h)(θ₁₂−3h)²(θ₁₂−4h)(θ₁₂−5h); f₄₃ = (θ₁₂−h)·f₃₃.
  - **Result:** all 14 non-top HW states recovered. Only [6,3,1] needs per-copy cleaning; [7,2,1] works sector-wide but is weakly conditioned (second singular value 0.0069).
- **[4,3,0]²** (`Minimal-ingredients-430-EH.wb`, 17 cells, 200 digits, about 25 s): this was the first localization for this rep.
  - **Site-1 local:** T₍₁,₁₎, T₍₂,₁₎, T₍₃,₁₎, T₍₄,₁₎ at θ₁−4h, so the pattern holds here.
  - **f's** (fitted at P1, checked at P2): f₁₁ = (θ₁₂−h)(θ₁₂−5h); f₂₁ = (θ₁₂−h)(θ₁₂−3h)(θ₁₂−5h); f₃₁ = (θ₁₂−h)(θ₁₂−2h)(θ₁₂−3h)(θ₁₂−5h); f₄₁ = (θ₁₂−h)·f₃₁.
  - **Result:** all 14 non-top HW states recovered. Only [7,5,2] needs per-copy cleaning; it is the dual image of [4,1,0]'s [6,3,1].
- **[5,1,0]²** (`Minimal-ingredients-510-EH.wb`, 19 cells, 250 digits; about 2.5 min plus 1 min for the 6b refit). This is the first λ₁ = 5 rep: dim 35, 1225 states.
  - **Site-1 local:** T₍₁,₁₎, T₍₂,₂₎, T₍₃,₃₎, T₍₄,₄₎, T₍₅,₄₎ at θ₁−5h. The [λ₁,1,0] rule is T₍ₖ,ₖ₎ (k < λ₁) + T₍λ₁,λ₁−1₎.
  - **f's** (t = θ₁₂/h): f₁₁ = (t−4)(t−6); f₂₂ = (t−3)(t−4)(t−5)(t−6); f₃₃ = (t−2)(t−3)(t−4)²(t−5)(t−6); f₄₄ = (t−1)(t−2)(t−3)²(t−4)²(t−5)(t−6); f₅₄ = (t−1)f₄₄.
  - **Result:** all 17 non-top HW states recovered. Only [8,3,1] needs per-copy cleaning; for [λ₁,1,0] the per-copy sector is always [λ₁+3,3,1].
  - **Convergence:** the 10⁻⁶ refit (Section 6b) shows O(ε) convergence, with errors 10⁻⁵·²–10⁻⁷·³ at (10⁻⁶, φ′).
- **[5,4,0]²** (`-540-EH.wb`, 250 digits, about 3 min): T₍ₖ,₁₎ for k = 1..5. Only [9,7,2] needs per-copy cleaning (the dual of [8,3,1], as predicted).
- **[5,2,0]²** (`-520-EH.wb`, 21 cells, 300 digits, about 9 min; 28 HW states, two multiplicity-3 sectors).
  - **Site-1:** T₁₁, T₂₂, T₃₃, T₄₃, T₅₃.
  - **Per-copy required:** [9,3,2], [8,5,1], [8,4,2], [7,5,2].
  - **New Section 6c:** weakly conditioned full-rank sectors ([7,6,1] σ₂ = 0.025, [7,4,3] σ₂ = 0.019) amplify the O(ε) error by about 1/σ₂. Per-copy cleaning fixes them.
- **[5,3,0]²** (`-530-EH.wb`, 21 cells, 300 digits, about 8.5 min): T₁₁, T₂₂, T₃₂, T₄₂, T₅₂. Per-copy needed for [9,5,2], [8,7,1], [8,6,2], [8,5,3], the exact duals of [5,2,0]'s.
- **All ten λ₁ ≤ 5 reps (λ₁ > λ₂ > 0) are now done.** All non-top HW states are recovered with O(ε) errors (≤ 10⁻⁵·⁵ at ε = 10⁻⁶ with the best cleaning choice). The summary is Section 3.7 of `Clean-CG-combinations`.
  - **Unified site-1 rule:** with d = λ₁−λ₂, the insertions are T₍ₖ,ₖ₎ (k ≤ d) and T₍ₖ,d₎ (d < k ≤ λ₁), all at θ₁−λ₁h. The leading prefactor is f₁₁ = (θ₁₂−dh)(θ₁₂−(λ₁+1)h).
- **Not done:** site 2 for λ₁ ≥ 4; closed forms of the clean coefficients; a second chain point for the pipeline.

## >>> (2026-10-02, latest): adjoint clean combination derived analytically <<<

**Notebook:** `Paul/Mathematica/new-conventions/Adjoint-clean-analytic.wb`. Exact, symbolic t = θ₁₂/h, E basis, about 1 s.
- **Mechanism.** Ω lies in the 27, so (W − T^top Ω) vanishes at ε⁰ and top-killing removes the O(ε) contamination identically. The leading contamination is ε²(W − T^top Ω)₁·ψ₁; only first-order perturbation theory is needed.
- **Which states need cleaning.** States with an O(ε) signal ([4,1,1], [3,3,0], octets) are automatically clean. The singlet has an O(ε²) signal; in raw eigenvalues K₁₁/S₂₁ = −1/3 and K₂₁/S₂₁ = (2t−3)/3, independent of φ and ω.
- **The clean combination.** 𝒞 = T₂₁(θ₁−2h) + h(2θ₁₂−3h)T₁₁(θ₁−2h), with Q-number (𝒞^B − 𝒞^top)⟨Ω|Ψ⟩. This is the fitted β = h(2t−3)/(t−1), now derived.
- **In Q's.** T₁₁ = 3h²χ₃(θ₁₂−h)(θ₁₂−3h)·r and T₂₁ = −3h³χ₃(θ₁₂−h)²(θ₁₂−3h)·ρr, with r = Q_c(θ₁−h)/Q_c(θ₁−2h) and ρ = Q^j(θ₁)/Q^j(θ₁−h); checked to 10⁻⁵⁰ on 18 states. Hence the clean functional ∝ {r[(2θ₁₂−3h) − (θ₁₂−h)ρ] − z₁⁻¹[(2θ₁₂−3h) − (θ₁₂−h)z₃]}·⟨Ω|Ψ⟩.
- **Rescaling (Section 3).** Rescaling the monodromy by g(u) changes only the number g(θ₁−h) in the combination.
  - ν₁, ν₂, ν₃ and the f-normalization all leave t-dependence (ν₂ is singular).
  - Rescaling each Lax by (u−θ_s−h/2) gives 𝒞̃ = T̃₂₁ − (4/3)T̃₁₁, which is t-free. The half shift on site 2 is forced.
- **Next:** the same derivation for [3,1,0]/[4,2,0], where deeper states need higher orders (ψ₂, …); also test whether the half-shift normalization makes those clean conditions t-free too.

## >>> (2026-10-02, later): minimal notebook for the contamination; Paul experimenting <<<

**Notebook:** `Paul/Mathematica/new-conventions/Minimal-ingredients-EH.wb` (16 cells, self-contained, clears the kernel, runs in ~10 s, saved, not committed). Built for Paul to experiment with. It contains only the EH generators, `tL2`/`tL1`, T₍₁,₁₎/T₍₂,₁₎ (two-site and one-site), `Rrot`, `Om1`/`Om2`, `x0`, the small-twist HW-limit states, the four localization identities, c_λ, and `Qside` (= T_λ^B⟨Ω|Ψ⟩/f_λ computed from the eigenvector; no Q's).
- **Exact split, at any ε:** `Qside` = c_λ·P_HWΨ (signal) + c_λ·(R⊗R)(1−P_HW)Ψ (contamination), because (R⊗R)P_HW = P_HW.
- **It reproduces every κ in Section 11.** The chain residual is about 10⁻⁴⁶ at 80 digits.
- **`showContam` decomposes the contamination by weight ν.** At [3,3,0] (2,1), the whole contamination is ν = [4,2,0]: c^{(420)} ~ ε⁰ times [(R⊗R)Ψ⊥]₄₂₀ ~ ε⁻², the same order as the signal ε¹·ε⁻³. The octet behaves the same way, and C removes it completely.
- **Discussion conclusions (not in any notebook):**
  - **Expanding in w cannot replace a weight projector.** Every ν term is homogeneous of the same torus degree top−μ in (r₁,r₂,r₃).
  - **A third independent r would add nothing.** The missing direction is a torus direction, and ρ = r₂/(r₁r₃) already depends on w₂.
  - **Only insertions diagonal on Ψ (combinations of T's) can serve as the "RHS projector."**
- **Section 7 (cells 17–23): a + b T₍₁,₁₎ + c T₍₂,₁₎ removes the contamination.** The bare Ω row is the covector e₁⊗e₁, which is pure top weight.
  - **Killing the top weight with a = −3h²z₂z₃(b − h z₃ c) is enough for every state except the singlet.** κ → 1 for [4,1,1], [3,3,0] and both octets, with either T alone.
  - **The two top-killed covectors have rank 2 on the octet HW space**, so the octet mixing is fixed.
  - **The singlet needs b/c = β = h(2t−3)/(t−1).** Its contamination ratios are K_b/S = −(t−1)/(3h) and K_c/S = (2t−3)/3. These were recognized at P1 and **validated at P2**, and are φ- and ω-independent.
  - **The single combination v* = (a, β, 1) gives κ → 1 for all five non-top HW states at P1 and P2.** So the HW CG data come directly from Q-expressions, with no reconstruction.
  - **Not done:** [4,2,0], and deriving β(t).
- **[3,1,0]² (`Minimal-ingredients-310-EH.wb`, 18 cells; generic in Λ, with Λ = {2,1,0} as the control).** Sparse 225-dim matrices at 120 digits; a fresh run takes about 20 s. Saved, not committed.
  - **Ingredients:** site 1 (T₍₁,₁₎, T₍₂,₂₎, T₍₃,₂₎ at θ₁−3h) + ⟨Ω|Ψ⟩ only. Localization holds to 10⁻¹⁰⁶.
  - **Top-killing is exact and works for any rep:** ⟨Ω|Ψ⟩ is the top component of (R⊗R)Ψ⊥. It cleans [6,1,1], [5,3,0] and [5,2,1]×2.
  - **The deep states keep contamination, but it is structural** (φ- and ω-independent).
  - **v₁ = 3d₁₁ − 2d₂₂** is clean on [5,3,0], [5,2,1]×2, [4,4,0], [4,3,1]×2 and [4,2,2].
  - **[5,2,1] mixing is fixed (rank 2). [4,3,1] mixing is not.** The only combination clean on both copies is v₁, and v₁'s direction on HW_[4,3,1] is w-independent (exact rank 1).
  - **Appendix A (independent scan):** the only nontrivial site-1-local CBR shapes (λ₁ ≤ 4, c = −6..2) are T₍₁,₁₎, T₍₂,₂₎, T₍₃,₂₎ at θ₁−3h.
  - **Section 8: adding site-2 A = t₂(θ₂+h) and B = t₂(θ₂+h)t₁(θ₂−h) fixes everything.** Their covectors are computed as ⟨Ω|X(R⊗R)⁻¹ and need no prefactors.
    - All contamination ratios are structural.
    - Clean spaces per weight sector (null space of the copies' k-vectors) have rank 2 on both [5,2,1] and [4,3,1].
    - **End to end:** the Q-side numbers alone recover P_HWΨ, scale included, for all 9 non-top HW states. Errors are O(ε), ≲10⁻⁶ at ε = 10⁻⁵.
    - **Caveat:** the clean combinations are numerical, fitted from the states' own k-vectors. They are structural, but not yet in closed form. One chain point only.
- **[4,2,0] blind spot corrected (`Blindspot-420-check-EH.wb`, 6 cells; one-site objects plus two-site HW vectors).**
  - **It is not blind.** Site-1 T₍₃,₂₎ and T₍₄,₂₎ have nonzero overlap with HW_[5,4,3] (∝ε³), and T₍₄,₂₎ with HW_[4,4,4] (∝ε⁴); the localized u stay O(1).
  - **What the old test showed.** Its "overlap sv = 1.000" (printed at 3 digits) only meant near-containment in the numerical null space. A correction note has been added to `L2-420-localization.md`.
  - **Done, see below.**
- **[4,2,0]² clean combinations (`Minimal-ingredients-420-EH.wb`, 13 cells, 250 digits, about 1 min). Saved, not committed.**
  - **Insertions:** site-1 T₁₁, T₂₂, T₃₂, T₄₂ at θ₁−4h; site-2 t₂(θ₂+h), T₂₂(θ₂), T₂₂(θ₂)t₁(θ₂−2h); plus ⟨Ω|Ψ⟩.
  - **Bare-frame bookkeeping:** S = W·P_HWΨ and K = W·Ψ⊥ − (W·e₁₁)⟨Ω|Ψ⟩, so no R⊗R is needed.
  - **Leading order is enough:** every state has p_S = p_K or better.
  - **Contamination is structural**, with an O(ε) spread.
  - **Sector clean spaces fix every multiplicity space except [6,4,2]** (mult 3; the common clean rank is 2). Per-copy cleaning fixes it.
  - **Fair end-to-end test:** the Q-numbers recover all 18 non-top HW states, scale included, with O(ε) errors (fitted at ε = 10⁻⁶, tested at 10⁻⁴/10⁻⁵ and at another φ). **That includes [5,4,3] and [4,4,4].**
  - **Precision trap:** 120 digits gave fake "0" errors (zero accuracy) for [4,4,4] at ε ≤ 10⁻⁵. Always check Precision.
  - **Open:** closed forms for the clean coefficients (they should be constants in t, like β(t)); a second chain point; a principled way to pick the [6,4,2] per-copy combinations.
- **Site 2 is not essential (2026-10-02).** With site 1 + ⟨Ω|Ψ⟩ only, plus per-copy cleaning where the sector-wide clean space lacks rank, all HW states are recovered (fair test, O(ε) errors).
  - **[2,1,0]²:** sector-wide is enough.
  - **[3,1,0]²** (`Minimal-ingredients-310-EH.wb` Section 9): only [4,3,1] needs per-copy cleaning.
  - **[4,2,0]²** (`Minimal-ingredients-420-EH.wb` Section 3): [7,3,2], [6,5,1] and [6,4,2] need per-copy cleaning; [5,4,3] and [4,4,4] work sector-wide.
  - **What site 2 buys:** fewer sectors with copy-specific coefficients (none for [3,1,0]; only [6,4,2] for [4,2,0]).
  - **Trap:** an already-clean insertion has k = O(ε), so use an absolute cut (`kTol`), not a relative one.
- **Site-1 eigenvalues in dual Q's (2026-10-02; `310` Section 10, `420` Section 4).** Let v(u) = (Q¹,Q²,Q³)(u) be the dual Q's.
  - **Structure near θ₁:** v(θ₁−Λ₂h) ∥ … ∥ v(θ₁), and v(θ₁−(Λ₁+1)h) … v(θ₁) lie in one 2-plane.
  - **Formula:** T_λ(θ₁−Λ₁h) = P_λ · det₂[v(θ₁+(λ₁−Λ₁)h), v(θ₁+(λ₂−Λ₁−1)h)] / det₂[v(θ₁−Λ₁h), v(θ₁−(Λ₁+1)h)], with P_λ = ∏_boxes ν₃(θ₁−Λ₁h+(c−r)h). This is the L=2 version of L=1 row 13.
  - **Validated:** exactly for 34 states of [3,1,0]² and 91 of [4,2,0]², all shapes with λ₁ ≤ Λ₁, every row pair.
  - **[4,2,0] explicit:** T₂₂/f₂₂ = κD₂₃/D₄₅, T₃₂/f₃₂ = −2hκρ₁D₂₃/D₄₅, T₄₂/f₄₂ = 2h²κρ₀D₂₃/D₄₅, with κ = 240h⁴(θ₁₂−3h)(θ₁₂−4h)/((θ₁₂−h)(θ₁₂−2h)); T₁₁/f₁₁ = 20h²(θ₁₂−4h)/(θ₁₂−2h)·D₃₄/D₄₅. Here ρ_k = Q^j(θ₁−kh)/Q^j(θ₁−2h).
  - **Not done:** the 6×6 ⟨Ω|Ψ⟩ formula has not been validated beyond the adjoint.
- **Single-Q form (`420` Section 5, 91 states, ≥225 digits).**
  - **Dual QQ, confirmed:** z_b q^a(u)q^b(u+h) − z_a q^b(u)q^a(u+h) = (z_b−z_a) G(u) q_c(u+h), with G = ∏_s(u−θ_s+h)(u−θ_s+2h). Note the shift +h on the direct q_c.
  - **Direct parallelism:** the direct Q vectors are parallel at θ₁−4h, −3h, −2h.
  - **Insertions:** T₁₁/f₁₁ = 10h²χ₃ r₁; T₂₂/f₂₂ = 40h⁴χ₃² r₁r₂; T₃₂/f₃₂ = −80h⁵χ₃² r₁r₂ρ₁; T₄₂/f₄₂ = 80h⁶χ₃² r₁r₂ρ₀. Here r₁ = Q_c(θ₁−3h)/Q_c(θ₁−4h) and r₂ = Q_c(θ₁−2h)/Q_c(θ₁−3h) (direct, any c), and ρ_k = Q^j(θ₁−kh)/Q^j(θ₁−2h) (dual, any j). The θ₁₂-dependence cancels against f.
  - **Next:** rewrite the fitted clean coefficients in (r₁, r₂, ρ₁, ρ₀); test [3,1,0].
- **Clean combinations, first look (`420` Section 6).**
  - **Top-kill** = eigenvalue difference: (X^B − X^top)⟨Ω|Ψ⟩.
  - **Leading signal ratios s̃** over (T₁₁, T₂₂, T₃₂, T₄₂)/f are simple rationals, e.g. [7,5,0] (1,2,−2,1) and [6,3,3] (0,0,−1/5,1).
  - **Contamination k̃:** 0 for [8,3,1] and [7,5,0] (only the top above them). [8,2,2]: (t−1)·s̃^[8,3,1]. [6,6,0]: (t−1)·s̃^[7,5,0]. In words, it is the signal pattern of the irrep one simple root above, times (t−1); this holds pointwise in w, and the raw factor is w-dependent.
  - **Deeper states:** k̃ is not rational at P1, presumably mixing several irreps above, including the multiplicity copies (algebraic data).
  - **Next:** expand k̃ in the s̃'s of the states above, using extra (non-localized) insertions so the span test is meaningful.

## >>> (2026-10-02): independent polynomial Baxter spectral solver (Section 14) <<<

**Read first:** Section 14 of `Paul/Mathematica/new-conventions/L2-diagonal-frame-overlap.md`. Notebook: `L2-diagonal-frame-overlap-EH.wb`, explicitly saved at revision 675 with 114 cells. Section 13 (cells 90–108) remains the 15/15 audit; new Section 14 is cells 109–114. Cells 110, 112 and 114 were individually executed successfully after fixes. No commit was made.

- **The solver no longer needs physical spectral input.** `specAnsatz` builds monic dual Q polynomials of degrees μⱼ, transfer polynomials t₁,t₂ with weight/twist-fixed leading terms, and known quantum determinant t₃. The ten unknowns (six Q coefficients, four transfer coefficients) satisfy polynomial dual Baxter identities. No physical matrices or eigenvectors enter this solve.
- **Numerical solve.** `specSolve` starts with exact-rational `NSolve`, then `specNewton` performs four scaled overdetermined Newton/least-squares corrections at working precision `pr+80`, seeded only by the Baxter roots. It rounds to the requested precision and checks the stored coefficients by rationalizing them in the equations. Directly increasing `NSolve` to `2pr` was too slow and was aborted; precision comes from polynomial refinement, not seed padding. No kernel operation remains in flight.
- **Validation is downstream.** `specReconstruct` uses spectral solutions, local 8×8 covectors, determinant overlaps, three bra twists, and x₀ normalization. Only afterward does `specValidate` build physical transfer matrices and check eigenvector equations coefficientwise with guard precision. This complete finite-twist test does not use a physical `Eigensystem` as spectral input.
- **Finite-twist result.** At h=1, θ=(5/2,0), z=(2,3,5), all 19 weight sectors yield 64 solutions, matching physical weight dimensions; Q reconstruction succeeds for all 64. The pre-final-refinement run reported max Baxter residual 3.34×10⁻⁹⁴, physical eigenvector residual 4.02×10⁻⁹², and normalization residual 9.28×10⁻⁹³ at 100 digits. The final code also succeeded, but its exact residuals were not recorded; describe it as roughly 90-digit validation rather than assigning new figures. Bra twists: (7/3,11/4,4/7), (9/4,13/5,3/7), (12/5,17/6,5/8).
- **Near-untwisted result.** Along rational paths z=1+εφ and w=1+εω (the same first-order direction as exponentials, but not exponential twists), all 21 states in the five HW-weight sectors (dimensions 1,2,2,6,10) now pass; six HW branches are selected by raising-operator residual. At P1, ε=10⁻⁶/pr=200, max Baxter residual 1.873621426×10⁻¹⁷³, physical residual/norm 1.809978519×10⁻¹³², and CG direction error 1.864923519×10⁻⁷. At P1 ε=10⁻⁹/pr=280 these are 4.779807849×10⁻²⁴², 3.397139631×10⁻¹⁸¹, and 1.864923680×10⁻¹⁰. At P2 ε=10⁻⁹/pr=280 with alternate φ=(-2/5,3/8,7/11), they are 8.105060246×10⁻²⁴⁶, 2.811004503×10⁻¹⁸¹, and 2.638405385×10⁻¹⁰. Primary φ=(1/3,−2/7,5/11); the three ω directions are recorded in Section 13. The small CG-direction error is O(ε), distinct from the much smaller finite-twist residual.
- **Scope.** This supplies a numerical polynomial Baxter solver and validates finite examples; it is not an analytic closed-form or general completeness proof, nor an arbitrary-representation implementation. The untwisted eigensystem is only a verification reference in the limit comparison.

## >>> NEWEST (2026-10-01): return to [2,1,0]² — symbolic factors and full Q reconstruction <<<

**User priority:** Paul asked to return to [2,1,0]², verify the earlier work and improve it. The [4,2,0] notebook was opened read-only; no new [4,2,0] calculation or edits were made.

**Read first:** Section 13 of `Paul/Mathematica/new-conventions/L2-diagonal-frame-overlap.md`. Code: `L2-diagonal-frame-overlap-EH.wb`, cells **90–108**. Rerun **3–14, then 90–108**. The audit ledger has **15 passing checks**.

- **Real correction:** site-2 B=t₂(θ₂+h)t₁(θ₂−h) has scalar coefficient **β=6h⁴χ₃(θ₁₂−h)(θ₁₂−3h)**. The old 3h³ expression only coincided at h=1/2. Cell 81 and the write-up are corrected.
- **All four localization formulas are now symbolic identities.** Undo R, use an arbitrary lower-triangular auxiliary twist (six independent entries), prove entrywise at h=1, θ=(d,0), then restore h and θ by homogeneity/translation. Also checked exactly at three rational chain/twist points.
- **Simpler overlap:** with D_B=det N[B,vac_w], P_B=Ψ_Q^B(x₀), B₀=e₁⊗e₁,
  `ΩΨ_B/(x₀Ψ_B) = D_B P_B₀ / [D_B₀ P_B (z₂−z₁)²]`.
  This cancels y₀, Φ_Q and c₀. Either bracket convention gives the same ratio.
- **The previously missing Q-side test is closed.** Reconstruct using this determinant ratio × transfer eigenvalues recovered from dual Q's, **not M·Ψ**. All 64 baseline states and all six small-twist HW branches pass. At P1/P2, ε=10⁻⁶/10⁻⁹, 240/300 digits, relative errors are below 4×10⁻¹⁶¹ / 3×10⁻¹⁹¹. An independent 160-digit run agrees with the 240-digit result to 7.73×10⁻⁷⁸.
- **Actual CG limit checked:** direction errors fall as O(ε), below 2.64×10⁻¹⁰ at ε=10⁻⁹, including an alternate φ and both octets.
- **Three w values suffice and are minimal for the four operator rows + x₀.** Exact ranks 1,2,2,6,10. Two w yield at most 9 equations for the 10d singlet weight sector. This is not a minimality claim if bare Ω or other rows are added.
- **Exact site-2 HW overlap table** in Section 13e and the `.md`, in terms of the six entries of g⁻ᵀdiag(z)gᵀ.

**Scope / next:** At the time of Section 13, test Q-data still came from Baxter equations fed transfer eigenvalues from the test eigensystem. Section 14 supersedes that limitation for the validated adjoint-square examples with a numerical polynomial Baxter solver; it does not establish general completeness. The historical summary in cells 88–89 includes old rejected approaches; use the new audit ledger. No commit made.

## >>> NEWEST (2026-10-01, final): [4,2,0]², φ-independence, blind spots <<<

**Read first.** `Paul/Mathematica/new-conventions/L2-420-localization.md`; the notebook is `L2-420-localization-EH.wb`, with 20/20 checks passing. Also Section 7 of `L2-320-localization-EH.wb` (φ).

**Results:**
- **[4,2,0], site 1.** T₍₁,₁₎, T₍₂,₂₎, T₍₃,₂₎ and T₍₄,₂₎ at θ₁−4h.
- **[4,2,0], site 2.** t₂(θ₂+h) and T₍₂,₂₎(θ₂).
- **The site-1 pattern across all four reps:** T₍₁,₁₎ plus T₍ₖ,λ₁−λ₂₎ for k = λ₁−λ₂..λ₁, all at θ₁−λ₁h, with f₍₁,₁₎ = (θ₁₂−(λ₁−λ₂)h)(θ₁₂−(λ₁+1)h).
- **[4,2,0] reconstruction.** Localized data + ⟨Ω|Ψ⟩ recover 16 of 19 HW states.
  - **[5,4,3]×2 and [4,4,4] are structurally invisible:** HW_μ lies entirely in the null space, so no ε-order, w or φ helps.
  - **Adding non-localized t-rows** ⟨Ω|t_a(u)| at generic u, which are still Q-computable, recovers all 19.
- **φ-independence.** The ε → 0 limit basis does not depend on φ. It is the untwisted inhomogeneous eigenbasis on HW_μ.
- **Open (Paul to decide).** Are non-localized t's acceptable? If not, a new localized ingredient is needed for [5,4,3].
- **Open, Paul picking this up 2026-10-02: the site-1/site-2 asymmetry.** Site 1 has many more localized objects than site 2 (spans for [2,1,0], [3,1,0], [3,2,0], [4,2,0]: site 1 = 3, 4, 4, 5; site 2 = 3, 3, 4, 4). Paul will look for more site-2 objects.
  - **What the search covered:** CBR shapes λ₁ ≤ 6 at θ₂ + c·h for c = −9..9, which includes the fused t\* objects; products by t₁ and t₂ to depth 2; all pairs and triples of t's at θ₂ + ℤh.
  - **What it did not cover:** non-integer shifts, derivatives in u, and the site-reversed monodromy.

## >>> (2026-10-01, end): [3,2,0]² site localization <<<

**Read first.** `Paul/Mathematica/new-conventions/L2-320-localization.md`. The notebook is `L2-320-localization-EH.wb`, the [3,2,0] copy of the [3,1,0] notebook; it reruns in about 32 s and 16/16 checks pass.

**Results:**
- **Site 1.** **T₍₁,₁₎, T₍₂,₁₎ and T₍₃,₁₎ at θ₁−3h**, which is Paul's guess.
- **Site 2.** t₂(θ₂+h) and **T₍₂,₂₎(θ₂)**.
  - The [2,1,0]/[3,1,0] product t₂(θ₂+h)t₁(θ₂−h) is **not** local here.
- **The f's swap sites relative to [3,1,0]**, consistent with [3,2,0] ≅ [3,1,0]*⊗det³.
- **Reconstruction.** site1 + site2 + ⟨Ω|Ψ⟩ + x₀ recovers **all 10 HW states exactly**.
  - Without ⟨Ω|Ψ⟩, the bottom irrep [4,3,3] is missed, just as [3,3,2] was for [3,1,0].
- **ε-depth (Section 6).** These are the orders beyond leading needed to fix each CG state, using 8 values of w:
  - [6,4,0]: 0;
  - [6,3,1], [5,5,0], [5,4,1]: 1;
  - [6,2,2], [5,3,2], [4,4,2]: 2;
  - [4,3,3]: 5.

  The data are needed through ε^{N−k}, where Ψ ~ ε⁻ᵏ with k = 4..9. This has not been run for [2,1,0] or [3,1,0].

## >>> (2026-10-01, very latest): [3,1,0]² site localization <<<

**Read first.** `Paul/Mathematica/new-conventions/L2-310-localization.md`. The notebook is `L2-310-localization-EH.wb`; it is self-contained and generic in Λ, with [2,1,0] as the control, and reruns in about 25 s. 18/18 checks pass.

**Results:**
1. **Site 1.** Everything is at **θ₁−3h**: **T₍₁,₁₎, T₍₂,₂₎ and T₍₃,₂₎**.
   - **Not T₍₃,₁₎ or T₍₂,₁₎**, at any integer shift.
   - **Factorization.** Each is f_λ × (one-site T_λ at θ₁) ⊗ Ω_zw, with f₍₁,₁₎ = (θ₁₂−2h)(θ₁₂−4h); f₍₂,₂₎ and f₍₃,₂₎ are in the `.md`.
2. **Site 2.** The same as [2,1,0]: t₂(θ₂+h) and t₂(θ₂+h)t₁(θ₂−h), with prefactor (θ₁₂−h)(θ₁₂−4h).
   - The search covered all CBR shapes, t\*, and pairs and triples of t's.
3. **Reconstruction.** site1 + site2 + x₀ reconstruct **9 of the 10 HW states exactly**, including both multiplicity-2 pairs ([5,2,1] and [4,3,1]).
   - **The exception is [3,3,2]**: the (3,3,2) rank is 16 of 18.
   - **Adding the bare ⟨Ω(w)|Ψ⟩ row** makes all 10 exact.

**Decided (Paul):** ⟨Ω|Ψ⟩ is fine as data. So site1 + site2 + ⟨Ω|Ψ⟩ + x₀ gives all 10 HW states of [3,1,0]² exactly.

**Next:**
- Other reps ([3,2,0], [4,2,0]): just change Λ.
- A second (h, θ) point.

## >>> (2026-10-01, latest): site-2 localization closes the octet/singlet gap (Section 12) <<<

**Read first.** "Section 12" of `Paul/Mathematica/new-conventions/L2-diagonal-frame-overlap.md`. The code is in `L2-diagonal-frame-overlap-EH.wb`, cells 79–87, and needs cells 3–14 plus the Section 11 helpers.

**Results:**
1. **t\* reduces to ordinary t's.** t\* (from T⁻¹ = ℒ₂⁻¹ℒ₁⁻¹) satisfies t\*₁(v) ∝ t₂(v+2h) and t\*₂(v) ∝ t₁(v+2h).
2. **The site-2-local objects.** The only ones at θ₂ + ℤh are:
   - **⟨Ω|t₂(θ₂+h)**, which is ∝ t\*₁(θ₂−h);
   - **⟨Ω|t₂(θ₂+h)t₁(θ₂−h)**.
   - **Negative findings.** No CBR T₍₂,₁₎ works. The search covered all CBR shapes and all products of two t's.
3. **Spectator factors (one-point fit).** The site-1 factor is f₍₁,₁₎ = (θ₁₂−h)(θ₁₂−3h), the same as on site 1. The product version has the extra piece α·(one-site t₂t₁) + β·Ω_zw.
4. **What the site-2 data recover.** Site-1 data (Section 10) + site-2 data + ⟨x₀|Ψ⟩ = 1 reconstruct **all six HW CG states exactly at small ε**: both octets and the singlet scale. This closes the Section 11 gap.
   - **Bonus.** Adding the bare ⟨Ω(w)| row alone already fixes the octet, but not the singlet scale.
5. **Caveat.** The reconstruction used M·Ψ in place of the Q side. That is exact by identity, but the full Q chain was not re-run.

**Next (Paul to decide):**
- Closed-form site-2 CG overlaps (the 10d analogue).
- A second (h, θ) point for the fits.
- [3,1,0]².

## >>> (2026-10-01, later): ε→0 of the mixed-twist overlap and the CG chain (Section 11) <<<

**Read first.** "Section 11" of `Paul/Mathematica/new-conventions/L2-diagonal-frame-overlap.md`, the full record. The code is in `L2-diagonal-frame-overlap-EH.wb`, cells 50–78, and needs only cells 3–14 run first. Two chain points: P1 = (1/2, (2/7,−3/5)) and P2 = (1/3, (1/5,−4/9)). Twists are z = e^{εφ}, w = e^{εω}, and t = θ₁₂/h.

**Results:**
1. **The ratio ρ.** ⟨x₀|y₀⟩²ρ_B → ε^# X_B ∏(ω_i−φ_j)^{max(d_j−D_i,0)}, with # = 8, 8, 8, 7, 7, 6 for [4,2,0], [4,1,1], [3,3,0], [3,2,1]±, [2,2,2].
   - **Constants:** X = 1, 1, 1; (2±√(4+5t²))/t² for the octets; −8/(t²−1) for the singlet.
   - **The ε⁸ is normalization only:** ⟨x₀|e₁e₁⟩ = (z₂−z₁)² and ⟨Ω|y₀⟩ ~ ε⁻¹⁰.
2. **The y₀-free ratio** R_B = ρ_B/ρ_[4,2,0]. The [4,1,1] and [3,3,0] rows are exact monomials in (w_i−z_j) at finite twist.
3. **The CG chain.** c_λ·(R⊗R)Ψ = T_λ⟨Ω|Ψ⟩/f_λ holds exactly, but "(R⊗R)Ψ → HW" fails at leading order.
   - **The cause:** c_λ's [4,2,0] component is O(1), and it pairs with the O(ε) correction to Ψ at the same order as the CG term.
   - **The size:** the contamination factors κ are φ- and ω-independent. For [3,3,0] κ = (t+1)/2.
4. **Partial fix.** C = c₂₁ + h z₃ c₁₁ has no [4,2,0] component and gives κ = 1 for [4,1,1] and both octets.
5. **What the allowed data can recover** (T₍₁,₁₎ and T₍₂,₁₎ at θ₁−2h for all w, plus x₀):
   - [4,2,0], [4,1,1] and [3,3,0] are reconstructed exactly.
   - **The octet mixing is undetermined:** there is an exact null vector at every ε, and it tends to a vector in the octet HW space.
   - The singlet's scale is not fixed. Coincident φ and ω don't help.

**Next (Paul to decide).** Fixing the octet needs one more independent ingredient. Paul's constraint is that only site-1-localized T's may be used. **Resolved by site-2-localized t's (Section 12, entry above).**

## >>> (2026-10-01, late): localized T_λ(θ₁−2h) action → highest-weight CG data <<<

**Read first.** The section "Section 10 (EH notebook only)" of `Paul/Mathematica/new-conventions/L2-diagonal-frame-overlap.md`. The code is in `L2-diagonal-frame-overlap-EH.wb`, Sections 10–10f (cells 28–49). Those sections need only cells 3–14 run first.

**Results:**
- **Localization.** Build the t's with the ket's twist diag(z). Then ⟨Ω_zw⊗Ω_zw|T_λ(θ₁−2h) = f_λ ⟨v_λ|⊗⟨Ω_zw| for λ = (1,1), (2,1), with ⟨v_λ| = ⟨Ω_zw|T_λ^{L=1,θ₁}(θ₁−2h).
  - **Exact** at the one numeric point.
  - **Site-2 factors, fitted at that point:** f_(1,1) = (θ₁₂−h)(θ₁₂−3h) and f_(2,1) = (θ₁₂−h)²(θ₁₂−3h).
- **Undoing R.** ⟨v_λ|R⁻¹ = e₁ᵀT_λ^{L=1}[g⁻ᵀ diag(z) gᵀ](θ₁−2h). This is symbolic in z, w and h.
- **Highest-weight overlaps.** For c_λ = ⟨v_λ|R⁻¹⊗e₁ᵀ the overlaps with the HW states of [2,1,0]² are known in closed form; the table is in the `.md`.
  - **(1,1):** misses [4,1,1] and the singlet, as forced by weights.
  - **(2,1):** reaches everything.
- **Small twists.** For z = e^{εφ}, w = e^{εω}, the leading-ε vectors (both λ, ω varied) span every weight sector of e_jᵀ⊗e₁ᵀ.
  - The 2-dimensional sector (3,2,1) needs both λ.
  - HW states are determined by their e_j⊗e₁ components.
  - Hence these vectors generate all the HW CG coefficients (Paul, 2026-10-01).

**In flight elsewhere.** Paul is doing the overlaps → CG step with another agent: assembling ⟨Ω|Ψ_B⟩ formula × T_λ^B(θ₁−2h), then R, then ε → 0.

**Open here:**
1. **Second numeric point.** Localization and the f_λ fit at a second numeric point, or symbolically (needs 64×64 objects).
2. **The pattern.** Test # = max(4−μ₁, μ₃) for the ε-orders.
3. **Other reps.** [3,1,0]² and others.

**Traps:**
- **`MatrixPower[nilpotent, 0]` fails.** Start the exp series from the identity.
- **Arbitrary-precision SVD on rank-deficient matrices** loses precision tracking. Use machine precision with a 10⁻⁹ tolerance for ranks.
- **The quadratic Casimir can't separate [4,1,1] from [3,3,0].**
- **Notebook state.** Sections 5–9 of the EH notebook were not re-run this session, so the summary table is stale.

## >>> CURRENT (2026-10-01, end of day): new conventions; L=2 overlap formula; Hermitian (EH) generators <<<

**Read first, in this order:**
1. `Paul/Mathematica/new-conventions/L2-diagonal-frame-overlap.md`: the L=2 result and the EH variant.
2. `Paul/Mathematica/new-conventions/Hermitian-structure-abstract.md`: whether the t's commute with their conjugates. The answer is no at L=2 with a twist.
3. `Paul/Mathematica/new-conventions/Companion-B-C-diagonalize-clean.md` (with `.tex`/`.pdf`): the L=1 reference.

**Conventions (2026-10-01).** Nothing is committed yet; ask Paul before committing.
- **Lax.** ℒ_ij = u δ_ij + h E_ij (it was u − h Σ E_ji⊗e_ij). At L=1 the twist is the companion G(z), acting on the right: T^G = ℒ·G. At L=2 all work is in the diagonal frame, with T·diag(z). The factors are ν_n = u − θ + hλ_n, multiplied over sites at L=2. `paper/xxxCG.tex` still uses the OLD conventions.
- **Generators (end of day).** From now on, use the **Hermitian generators EH_ij = R E_ij R⁻¹**, which satisfy EH_ijᵀ = EH_ji. R is diagonal, and M = R² solves M⁻¹E_ijᵀM = E_ji with M₁ = 1. For the adjoint, R = diag(1,1,1,√2,√6,3,2,3) in `ps` order.
- **Layout.** Pre-change notebooks live under `Paul/Mathematica/old-conventions/`; new work goes in `Paul/Mathematica/new-conventions/`.

**Files in `new-conventions/`:**
- **L=1 adjoint.**
  - `Companion-B-C-diagonalize-clean.wb` (`.md`/`.tex`/`.pdf`): E basis, 37 checks.
  - `Companion-B-C-diagonalize-clean-EH.wb`: EH basis, 40 checks, all passing. It adds t_a(G)ᵀ = t_a(Gᵀ), and t_a(diag z) is symmetric.
  - `Companion-B-C-diagonalize.wb`: the working and exploration notebook. It hard-codes h=1 in many cells.
  - `Baxter-QQ-Wronskian.wb`: Baxter, QQ, Wronskian and orthogonality checks for any rep.
- **L=2 adj⊗adj, diagonal frame.**
  - `L2-diagonal-frame-overlap.wb`: E basis.
  - `L2-diagonal-frame-overlap-EH.wb`: EH basis.
  - Both end in a summary table whose only 4 FAILs are deliberate records of approaches that don't work (Section 6, plus one non-equivalent row choice).
- **Abstract algebra.** `Hermitian-structure-abstract.wb`: a noncommutative normal-ordering engine, plus t₁ and t₂ for L sites, with † computed abstractly.
- **Other.** `test.wb` is Paul's connectivity test.
- **Spec.** `docs/superpowers/specs/2026-10-01-L2-diagonal-frame-overlap-design.md`.

**Established at L=1** (adjoint [2,1,0]; checked at (h,χ₃) = (1,1), (1,7/2), (1/2,5/4)). Details are in the clean `.md`. ⟨x₀|y₀⟩ = 9h⁴ in the E basis and 3h⁴ in the EH basis; everything else is basis-independent.
1. **Quantum determinant.** t₃ = χ₃ ν₁(u)ν₂(u−h)ν₃(u−2h).
2. **Baxter equations.** The direct and dual equations keep the old forms but **swap roles**: the new dual equation is the old direct form, and vice versa. Equivalently, each is the old one under h → −h with t_a(u) → t_a(u+(a−1)h). The solutions are Q_i = z_i^{−u/h} q_i with deg q_i = λ₁−μ_i, and Q^j = z_j^{+u/h} q^j with **deg q^j = μ_j − λ₃**. The −λ₃ was found at [2,1,1].
3. **QQ relation.** Q_a(u)Q_b(u+h) − Q_b(u)Q_a(u+h) = (1/z_b−1/z_a) χ₃^{−u/h} F_λ Q^c, with F_λ = ∏_{n=λ₂+1}^{λ₁}(u−θ+nh).
4. **Wronskian (quantum Weyl) formula** for every CBR T_λ in terms of the dual Q's: T_λ = ∏_boxes ν₃(x+(c−r)h) · D_{λ−δ}/D_{−δ}. Checked for all 19 λ ⊆ 3×3.
5. **Reduction at θ−2h.** v(θ) ∥ v(θ−h) and the columns v(θ−3h..θ) span a 2-plane, so T_λ(θ−2h) becomes a ratio of 2×2 determinants, or 0. Only (1), (1,1), (2), (2,1) survive, which are exactly the t's that build ⟨x|.
6. **SoV wavefunctions.**
   - ⟨x_μ|Ψ⟩ = N_μ det₂[v(x₂₁),v(x₁)] v_k(x₂₂)/(same at x₀), with N_μ = n_μ h^{k_μ}.
   - ⟨Ψ|y_μ⟩/⟨Ψ|y_low⟩ = N′_μ det₂[u(y₂₁),u(y₂₂)] u_k(y₁)/(same at y_low), using the direct columns u.
7. **Functional orthogonality from Baxter/dual adjointness.** Columns {−⟨P^jQ_i⟩, ⟨P^jQ_i^{[2]}⟩, ⟨P^{j[−2]}Q_i⟩}, bracket Σ_{n=λ₃}^{λ₁} f(θ−nh)/∏_{m≠n}(m−n)h, rows (1,3),(1,2),(2,3). This is the old relation with h → −h. It also holds for [3,1,0], [3,2,0] and [2,1,1].
8. **Normalized scalar product.** ⟨Ψ_A|Ψ_B⟩ = (1/⟨x₀|y₀⟩) · (det N[ket dual, bra direct]/c₀) / (Ψ_Q^B(x₀) Φ_Q^A(y₀)), with c₀ = 1/(4h⁶) read off by expanding det N, and ⟨x₀|y₀⟩ = 9h⁴. The normalizations are ⟨x₀|Ψ⟩ = 1 (x₀ = e₁) and ⟨Ψ|y₀⟩ = 1 (y₀ = y_{2|21}).
9. **The same formula for other twist pairings.**
   - Between **different twists**, ⟨Ψ_{A,w}|Ψ_{B,z}⟩, with each state carrying its own twisted Q's: all 64 entries match.
   - In the **diagonal frame**, Ω_zw·P_phys|Ψ_{B,z}⟩ = (Ω_zw P_phys·y₀) × (the formula with the w-vacuum as the bra).

**Established at L=2** (adj⊗adj, diagonal twist; one numeric point: h=1/2, θ=(2/7,−3/5), z=(3/2,5/4,2/3), w=(5/2,7/3,3/5)):
1. **Ports unchanged.** The quantum determinant, the Baxter equations with ν_n → ∏_s ν_n, the degrees deg q_i = 2λ₁−μ_i and deg q^j = μ_j−2λ₃, and the QQ relation with F = ∏_s(u−θ_s+2h).
2. **The bra.** Ω_zw⊗Ω_zw is the exact vacuum of t(G′), where G′ = P_aux G(w) P_aux⁻¹.
3. **The overlap formula.** Ω·Ψ_B = (Ω·y₀)/(c₀⟨x₀|y₀⟩²Φ_Q^vac(y₀)) · det N[B,vac]/Ψ_Q^B(x₀) · ⟨x₀|Ψ_B⟩, to 50–51 digits for all 64 states.
   - x₀, y₀, Ψ_Q and Φ_Q are the L=1 objects taken as products over sites.
   - The 6×6 det N has rows (i,j;k), with the L=1 pairs (1,3),(1,2),(2,3) and the bracket label k = 1, 2. Its columns are the three L=1 functionals times x^{β−1}.
4. **Brackets (Paul's preferred choice).** Both are residue sums around ALL poles, with insertions 1 and e^{2πiu/h}. Then c₀ = (ε₁−ε₂)³/[16h¹²θ₁₂³(θ₁₂²−h²)⁴(θ₁₂²−4h²)²], with ε_s = e^{2πiθ_s/h}.
   - Per-ladder brackets also work; there (ε₁−ε₂)³ is replaced by −1.
   - **A single bracket over all poles does NOT work.** It passes the orthogonality test but gives a row-dependent det N. See memory `reference_fsov_per_ladder_brackets.md`.

**Hermitian structure** (`Hermitian-structure-abstract.wb`, abstract and exact):
- **L=1.** t_a† = t̄_a (all scalars conjugated), so for real parameters the t's are Hermitian.
- **L=2.** t₁† − t̄₁ = h̄² Σ_{i≠k}(z̄_k−z̄_i) E⁽¹⁾_ik E⁽²⁾_ki, and

  [t₁(u), t₁(v)†] = −h³θ₁₂ Σ(z_i−z_k)² E⁽¹⁾_ik E⁽²⁾_ki + h⁴[C_z, C_z†], with C_z = Σ z_i E⁽¹⁾_ik E⁽²⁾_ki.
- **When it vanishes.** [C_z, C_z†] has cubic words that nothing else can cancel, so the commutator is zero **only for the untwisted chain**. It is not zero for a unitary twist, for partial equality of the z's, or for imaginary parameters.
- **What holds instead.** t_a(u)ᵀ = t_a of the **site-reversed** chain. This was checked numerically in the EH L=2 notebook, Section 1b. The pairing † ∘ site reversal (via the R-matrix that swaps the sites) is unexplored. Paul said to leave it ("that's fine").

**Traps (each cost time today):**
- **Surds with EH generators.** Exact entries contain surds that Mathematica doesn't combine canonically (√6/2 → √(3/2) next to √6). So `Union[...] === {0}` reports false FAILs. Use `PossibleZeroQ` (`zeroQ` in the EH notebooks).
- **Symbolic det N for c₀.** Zero the Q-symbols outside the target monomial *before* `Det`. The coefficients come out the same, and it is about 100× faster than expanding first, which ran for minutes.
- **Orientation of det N.** Pair the **ket's dual** Q's with the **bra's direct** Q's. At equal twists the other orientation looks fine; at different twists it is off by factors of about 100.
- **`Max[{}]`.** It returns `-Infinity` for an exactly-zero polynomial, and `Ordering` sorts that after numbers. Floor at 0, as `maxc` does.
- **Hard-coded, basis-dependent numbers.** For example ⟨x₀|y₀⟩ = 9h⁴. Use the computed values.
- **Wolfbook.**
  - Insertion order: after inserting k cells after N, the next insert goes after N+k (memory `reference_wolfbook_insert_position.md`).
  - Replies are truncated, so read results with single-line `ToString[..., InputForm]` evaluations.
  - Tool replies may name the wrong target notebook. Verify against the saved `.wb` JSON.
- **Shared kernel.** All new-conventions notebooks share one kernel, and each clean or L=2 notebook clears it on its first cell.

**Next steps (Paul has not set the priority — ask):**
1. **Other reps at L=2.** [3,1,0]² first, plus a second numeric point at adj⊗adj.
2. **General L.** Check whether the bracket structure generalizes to L brackets with insertions e^{2πiku/h}, k = 0..L−1.
3. **The paper.** Update `paper/xxxCG.tex` to the new conventions (E → EH, Lax, twist).
4. **The Casoratian.** Find a closed form for D_{1,0,−1} at L=1 (degree |λ|−3λ₃).
5. **The proposed relabelling** (paused by Paul): direct↔dual, upper↔lower, j → 4−j.
6. **Commit** today's work: new-conventions, the old-conventions move, the spec and these notes.

---

## Everything below predates the 2026-10-01 convention change

It refers to the OLD conventions. For paths, read `Paul/Mathematica/XXX-CG/...` as `Paul/Mathematica/old-conventions/XXX-CG/...`, and `Paul/Mathematica/old/...` as `Paul/Mathematica/old-conventions/old/...`.

**2026-09-25 — thread at the time: one-site localization / SoV bases.** Start at
`Paul/Mathematica/XXX-CG/Building-SoV-bases/HANDOFF-SoV-bases.md` (what was tried, what worked and
failed, traps, next steps). Results are recorded in `Building-SoV-bases/main-notes.tex` (a rolling
record: keep appending) and the computations are in `Building-SoV-bases/L2-localization-tests.wb`.

**2026-09-29 update (best result so far):**
- **Adjoint.** In `Building-SoV-bases/transfer-algebra-dim.wb`, a globally rotated vacuum plus R t11(θ1) gives a one-body covector that reaches the adjoint singlet.
- **[3,2,0].** t22(θ1)·t10(θ1+2h) reaches every irrep of [3,2,0]².
- **Still open (at the time).** [3,1,0] and [4,2,0] had gaps; these are now closed, see the later-session entry above.
- **Where it's written up.** `HANDOFF-SoV-bases.md` ("Latest"), the `main-notes.tex` summary section, and the data in `transfer-algebra-dim-results.m`.

**2026-09-29 later session (gaps closed):**
- **The operator.** At L=1, the SoV-basis operators are T_{μ11,0}T_{λ2,0}⁻¹T_{μ21,μ12}(θ1), with fusion by the CBR formula. At L=2 on the twist-w rotated vacuum, the top one, T_{λ1,0}T_{λ2,0}⁻¹T_{λ2,λ2}(θ1), is one-body.
- **What it reaches.** Every irrep of [3,1,0]², [3,2,0]², [4,2,0]² and [4,1,0]².
- **Single-site interpretation.** It acts as the same operator for site 1 alone, times an explicit prefactor from site 2.
- **Where to look.** `HANDOFF-SoV-bases.md` ("Latest"), `main-notes.tex` and `L2-twistw-vacuum-SoV-ops.wb`.

**2026-09-29 end of day (ratios are polynomial):**
- **The identity.** At θ1, T_{n,0}T_{λ2,0}⁻¹T_{λ2,λ2}(θ1) = T_{n,λ2}(θ1) as an exact operator identity, with coefficient 1.
- **Consequence.** O_max = T_{λ1,λ2}(θ1), and R t11 = T_{2,1}(θ1).
- **Where it holds.** Only at θ1, where the two-site t3 vanishes. It fails at generic u.
- **Caveat on "reach".** It means one covector ⟨v_g|O overlaps at least one copy of every irrep. It is **not** "generates all states", and multiplicity coverage was never recorded.
- **Where to look.** `HANDOFF-SoV-bases.md` ("Latest", item 4) and `main-notes.tex`.

**2026-09-23 reorg:** `Paul/Mathematica/XXX/` and `Paul/Mathematica/Gaudin/` were moved wholesale
to `Paul/Mathematica/old/XXX/` and `Paul/Mathematica/old/Gaudin/` (structure preserved, paths below
updated accordingly). New active XXX work goes in the fresh, empty `Paul/Mathematica/XXX-CG/`
instead of adding to the old tree. See `CLAUDE.md` for the current layout.

## !! General reminder for the next hard session (Paul, 2026-09-24) !!

Full write-up: `Paul/Mathematica/XXX-CG/Q-sampling-count.tex`/`.pdf`.
- **The relation.** All three Q_ab = (1/z_b − 1/z_a) F_λ Q^c, with F_λ = ∏_s∏_{n=λ2+1}^{λ1}(u−θ_s−nh).
  Checked for all states: [2,1,0]² at points 1 and 2, [3,1,0]² and [3,2,0]² at point 1 (`L2-mathcalA.wb`
  cells 12–13). The saved outputs are for [3,2,0].
- **[3,1,0]² count.** 2 ratios of q1 (4 with the determinant) and 1 ratio of q^3. Every reached state has
  deg q1 ≤ 2 and deg q^3 ≤ 1. The missed [4,2,2] and [3,3,2] both have deg q^3 = 2.
- **The q^3 count is the ν3 > λ2 rule.** It has λ2 ratios against deg q^3 = Λ3. The q1 count is not binding:
  [4,4,2] in [3,2,0]² is reached with deg q1 = 2 against 1 ratio.
- **[3,2,0]² count.** The only missed state, [4,3,3], is the only one with q^3 undersampled (a cubic at 3 points).
  Every reached state has q^3 covered. This holds across all three reps checked.
- **[4,2,0]² count.** Same result: the 4 missed states ([6,3,3], [5,4,3]×2, [4,4,4]) are exactly those with deg q^3 > 2
  ratios, and all 15 reached states are covered. So in all four reps, a state is missed ⇔ q^3 is undersampled
  (Λ3 > λ2). The Q_ab relation is NOT checked at [4,2,0] (729 states); F = ∏_s(u−θ_s−3h)(u−θ_s−4h) is conjectured.

**Why ⟨k1^{u/h}Q1⟩⟨k2^{u/h}Q12⟩ cannot give all CG states: count the data (Paul's view).**
- **The data.** The coordinate-space (SoV) CG wavefunction is built from Q-functions. To fully
  reconstruct a highest-weight state you need all its Bethe roots, i.e. the polynomials Q1 and Q12.
  In the twisted chain Q1, Q12 fix every state; in the untwisted chain, only up to its multiplet.
- **What the functionals sample.** The k1/k2 expansions sample Q1 and Q12 only at special points,
  short ladders of shifted inhomogeneities: Q1 at θ1+λ2h,…,θ1+λ1h and Q12 at θ1,…,θ1+λ2h. The
  determinant version adds the matching θ2 ladder for Q1. After normalization, only ratios within
  each ladder count.
- **The shortfall.** When a state has more unknown roots than there are sampled ratios, it cannot
  be reconstructed.
- **Count first.** Before building a new functional, count the points of Q1/Q12 it probes and
  compare with the unknown roots in the table below.
- **Way forward.** We need new generating functions that sample the Q's at more points and localize
  to each site separately (site-2 ladders, other shifts). They should supply the missing data.
- **Count the right unknowns.** Q12 = (1/z2 − 1/z1) ν1(u) Q^3(u), no shift, with
  ν1(u) = ∏_s (u−θ_s−λ1h) and q^3 monic. Checked exactly for all 64 states of the twisted
  [2,1,0]², point 1 (`XXX-CG/L2-mathcalA.wb` cell 12). So the unknowns in Q12 are the Λ3 roots of
  q^3, not deg Q12.
  - Not checked: the untwisted version, expected Q12 = (deg q1 − deg q2) h ν1 q^3.
  - Not checked: other reps. For λ1−λ2 > 1 the factor has degree L(λ1−λ2), probably
    ∏_s ∏_{n=λ2+1}^{λ1} (u−θ_s−nh).

Untwisted degrees, [2,1,0]⊗[2,1,0] (L=2, λ=(2,1,0)). The degrees are staggered:
deg q_i = (4−Λ_i)+(i−1) and deg q^j = Λ_j+(3−j) (`old/XXX/Experiments/SU3_L2_Untwisting.md` §1, `untwU`).
So deg Q12 = deg q1 + deg q2 − 1 = 2+Λ3. That is a trivial degree-2 factor times q^3, with deg q^3 = Λ3.
The unknown roots are then (4−Λ1) + Λ3 = M1 + M2, the Bethe-root count.

| irrep Λ | dim | deg q1 | deg q2 | deg Q12 | deg q^3 | unknown roots |
|---|---|---|---|---|---|---|
| [4,2,0] | 27 | 0 | 3 | 2 | 0 | 0 |
| [3,3,0] | 10 | 1 | 2 | 2 | 0 | 1 |
| [4,1,1] | 10 | 0 | 4 | 3 | 1 | 1 |
| [3,2,1] ×2 | 8 | 1 | 3 | 3 | 1 | 2 |
| [2,2,2] | 1 | 2 | 3 | 4 | 2 | 4 |

Also deg Q13 = deg q1 + deg q3 − 1 = 3+Λ2 = (degree-2 factor) × q^2, deg q^2 = Λ2+1.

The adjoint functionals probe Q1 at θ1+h, θ1+2h (1 ratio) and Q12 at θ1, θ1+h (1 ratio). That makes
2 numbers.
- **Linear cases.** One ratio fixes a monic linear q, since (θ1+h−a)/(θ1−a) determines a. So 2
  numbers cover every irrep with ≤2 unknown roots, i.e. all but the singlet.
- **The singlet.** It needs 4 and is exactly the state we missed. The count gives 5/6, matching
  what we reached.
- **Twisted chain.** The same count works per weight μ: deg q1 = 4−μ1, deg Q12 = 2+μ3.
- **Caveats (Claude).** Distinguishing finitely many Bethe solutions needs less than reconstructing
  them, so the count alone is a heuristic. The singlet overlap is also an exact zero (the ν3 > λ2
  rule), not just a degenerate value. A new functional should be checked against both: enough
  sampled points, and non-vanishing at the missed weights (μ3 > λ2).

## Project goal (Paul, 2026-08-27)

**Generalize the localization/CG construction to su(3) GENERIC representations.** Rectangular reps
are more or less understood. For generic (non-rectangular) reps localization *seems to stop
working*, though some hints remain. Closing that gap is the goal of the project, not a side thread.

gl(3) is the minimal setting where the question can even be posed: every gl(2) irrep is rectangular
up to a determinant twist (`{l1,l2} = {l1-l2,0} (x) det^l2`), so su(2) work can exhibit the
mechanism but never the obstruction. In gl(3), `{p,p,0}` and `{p,0,0}` are rectangular; `{l1,l2,l3}`
with `l1>l2>l3` is not.

Directly relevant prior work: the local-rotations thread (2026-07-18, `Experiments/XXX_CG_L2_V2.wb`
Sections G-L, notes in `Experiments/local_rotations_notes.tex`) already established that
localization fails for general reps, ruled out the spectral-parameter-insertion route, and
identified the **product-covector <-> SoV-basis recoupling matrix `A`/`B`** as the object to
generalize. Existing su(3) notebooks: `Paul/Mathematica/old/XXX/Experiments/su3_V1.wb`,
`su3_V2_general.wb` (arbitrary gl(3) Young diagrams, full-sweep verified), `su3_V2_general_Paul.wb`,
`su3_V2_fund_anti_fund.wb`. `Mathematica/su3/` exists but is empty.

## Current focus

Computing CG coefficients directly in the XXX spin chain (`Paul/Mathematica/old/XXX/`). The Gaudin
subproject proper (residue/representation groundwork) is still parked, but the `h→0` scaling limit
of the XXX construction below is an active connection point, not a revival of that parked track.

**Active thread (2026-09-19): twist-overlap reconstruction, section immediately below.** The
2026-09-17 rotations thread is its direct continuation (per-site rotations are the fix).

**Prior thread (2026-09-17): rotations of the SoV reference covector as determinants of
Q-integrals, gl(3) L=2 adjoint ⊗ fundamental. Self-contained in
`Paul/Mathematica/old/XXX/rotations/` — start at its `README.md`.** The gl(3) adjoint thread
(2026-09-16, below) and the gl(4) `N_untw` thread (2026-09-10) are paused, not closed.

## >>> NEWEST (2026-09-23, late): localizing SoV combinations reach the singlet directions (XXX-CG/) <<<

**Read first:** `Paul/Mathematica/XXX-CG/L2-localizing-transfer.tex`/`.pdf`, Sections 4–5 (Results 5–7). Notebook cells 10–16.
All exact at point 1 of [2,1,0]² (h=1, θ=(1/3,−8/57), z=(5/2,7/3,6/35)), finite twist.

- **Product covectors ⟨v⊗w| at θ1 (cells 10–13).** t_{2,1}(θ1)(Ĝ⁻¹⊗1) commutes with the global gl(3), so t_{2,1}(θ1) = R12 G1 exactly.
  - t_{2,1}(θ1) localizes only for v ∝ w, with w a local hw state of ANY Borel (same on both sites); then v' ∝ wĜ.
  - t_{1,1}(θ1) additionally allows a relative E12 rotation, v = e1ᵀe^{aE12}g with w = e1ᵀg. Its localizing subspace is covariant.
  - t_{1,0} and t_{2,0} localize only on the trivial eigen-covectors e1⊗e1 and e8⊗e8.
  - No t_λ(θ1) localizes a relative E13/E23/Weyl rotation. So mixed Borels such as e1⊗e8 never localize.
- **Minimal localizing SoV combinations (cells 14–16).** Basis ⟨x_{g1|g2}| = ⟨vacL⊗vacL|O_{g1}(θ1)O_{g2}(θ2) (per-site L=1 recipe of `SU3_L1_Baxter_Twisted.wb` cell 11; rank 64).
  - Σc_x⟨x| = ⟨v⊗vacL| is an 8-dim space using only 12 labels.
  - The minimal combinations are 4 singles (100|100, 111|100, 200|100, 211|100) and 4 pairs {g|200, g⁺|100} (g⁺ = n11 raised by one), with coefficients −19/9, −19/14, −95/203, −19/28 (19 from θ12 = 9/19).
  - Their site-1 covectors v·g0⁻¹ have rank 8, so they reach (1,0,2) and (0,1,2), which ⟨k1Q1⟩⟨k2Q12⟩ never reaches. ⟨x_{221|100}| − (19/28)⟨x_{211|200}| hits (0,1,2), the singlet component.
- **Next:** the analytic form of the pair coefficients (scan θ, z, p); a Q-functional realization of the pairs (a ⟨k1Q1⟩-type bracket with site-2 poles?); site-2 localization; other reps.
- **Same day, `XXX-CG/xxx-210-sov-limit.wb`** (untwisting limit of the SoV basis, [2,1,0]², z = 1 + φε):
  - Pointwise, all 64 SoV covectors → ∝ vacL⊗vacL (rank 1). The untwisted t's act as scalars on the multiplicity-one [4,2,0].
  - After undoing g0⊗g0 and projecting to total-weight sector ν, every nonzero component starts exactly at ε^{ht(ν)} (ht = simple-root steps below (4,2,0)). The graded limit has rank 64, so any (e1ᵀP)⊗e1ᵀ with P = exp(εaE12 + ε²bE13 + εcE23) is one ε-independent combination.
  - Not checked: that Ψ diverges with the same grading.
  - Cells 5/6 there (the "limit of span" hw-count) are frame-dependent; don't reuse them.
- `xxx-210-cg-test.wb` cell 16: explicit 8-component v_new (finite twist, factored table).
- **Traps:**
  - A multi-line `f[...] := a . b .` ending a line with `.` did not register in wolfbook (use one line).
  - Symbolic 64×64 LinearSolve / `Together` on symbolic inverses hang the kernel, and abort doesn't release it (restart). Work per weight sector (≤10-dim) or with truncated ε-series.

## >>> ACTIVE (2026-09-23): 𝒜 functionals ⟨k1^{u/h}Q1⟩, ⟨k2^{u/h}Q12⟩ at L=1 and L=2 (XXX-CG/) <<<

**Read first:** `Paul/Mathematica/XXX-CG/L2-mathcalA.tex`/`.pdf` (L=2), `L1-mathcalA.tex`/`.pdf` (L=1).
Notebooks `L1-mathcalA.wb`, `L1-orthogonality-diagonal-frame.wb`, `L2-mathcalA.wb` (switches `repT`,
`ptIndexT` in cell 3; ledger `resultsT` in cell 11). Normalization ⟨vacL⊗vacL|Ψ⟩ = 1 = ⟨Ψ|vacR⊗vacR⟩.

L=2, same rep on both sites, checked for [2,1,0], [3,1,0], [3,2,0] at 3 points (≤1e-77):
- 6×6 twisted FSoV det / (product of L=1 c^ref Ψ^ref) = ⟨Ψ_A|Ψ_B⟩ with one constant
  κ2 = (−1)^{λ1+1} h^{6λ1} κ1² [∏_{n=−λ1}^{λ1}(θ12−nh)]³/(e1−e2)³ (exact at points; sign/h-power fitted on λ1=2,3).
- ⟨k1^{u/h}Q1⟩ with the L=1 measure on site 1 only = ⟨v1 ⊗ vacL|Ψ⟩, v1 = the L=1 covector.
- ⟨k2^{u/h}Q12⟩ with site-1 poles only and measure 1/[∏_{k=0}^{λ2}ν3^{(1)}(u−kh) ∏_{n=λ2+1}^{λ1}(u−θ2−nh)]:
  product of both = ⟨v1⊗vacL| Σ_m (−1)^m C(λ2,m) k2^m t^{(1)}_{m,m}(θ1)⊗1 |Ψ⟩ — the L=1 formula exactly.
  Without the site-2 factor it holds only up to θ-dependent c_m; the full ν3(u)ν3(u−h) does not fix it.

**Which transfer matrices localize on a global rotation (2026-09-23)** — `XXX-CG/L2-localizing-transfer.wb`, [2,1,0]², exact,
T = L1 L2 G. ⟨u⊗u|t_{2,1}(θ1) = τ^vac_{2,1}(θ1) ⟨u^G⊗u|, u^G = e1ᵀG⁻¹gG (a,b,c → a z2/z1, b z3/z1, c z3/z2), symbolic in a,b,c:
site 1 only (the regularity picture, t(θ1) = R12 G1). Scan over θ_s+nh, n = −5..5: t_{1,1}(θ1) is also site-1-only, but its w is not a
rotation; t_{1,1}(θ2+3h) acts on site 2 only; t_{2,2}(θ_s±h) kills ⟨u⊗u|; no t_{2,1}(θ2+nh) is site-2-only (ordering). Formula re-checked at point 2 with u = vacL (cell 7).
CG count for ⟨vacL⊗vacL|t_{2,1}(θ1) (`xxx-210-cg-test.wb` cells 11–13): after undoing g0⊗g0, v_t = e1ᵀ G⁻¹g0Gg0⁻¹ = e1ᵀe^{yE23}e^{xE12}
(zero (1,3) entry from the p3 relation), so it is ≡ 0 at (1,0,2),(0,1,2) like v_new. Alone 4/6 hw, combined with old + det still 5/6: nothing new, singlet missed.
**Why the ν3 > λ2 rule — now derived (cells 14–15; tex sections in xxx-210 and a remark in xxx-310):** site-1 μ3 = ν3, and two mechanisms cap it:
(a) t̂ families (old, det, w-vac): rotation part e1ᵀe^{aE12} or e1ᵀ n…n⁻¹ has μ3 = 0, and each quantum minor moves ≤1 unit into slot 3 FOR ANY
TWIST (checked exactly with a general 3×3 twist), with ≤ λ2 insertions. The (1,3) = 0 of G1 is NOT the cause here. (b) t_{2,1} family: e1ᵀe^{yE23}e^{xE12},
zero (1,3) entry, and e1ᵀE23^k = 0 for k > λ2. The "everything else is reached" half of the rule is still observed only. A (1,3)
entry opens E13: hypothetical v^(2) = e1ᵀG⁻²g0G²g0⁻¹ has (1,3) = −p2(z1−z2)(z1−z3)/z1² = O(ε²), hits the singlet, all families + v^(2) = 6/6.
NOT realized by t_{2,1}(θ1)²: that is rank 4, not local (L2-localizing-transfer cell 8). Matches the 2026-09-19 finding (one per-site E31 rotation completes).

**6×6 vs vacuum at twist w (2026-09-23)** — `XXX-CG/L2-vacw-mathcalA.wb` (cell 1 = summary). The 6×6 FSoV with left Q^j
of Ω=e1⊗e1 at twist w is C⟨u⊗u|Ψ_B⟩, a GLOBAL rotation; flag of u: ℓ_k ∝ ℓ0_k (z_k−w2)(z_k−w3), π_k̂ ∝ π0_k̂/(z_k−w3)
((ℓ0,π0) = vacL flag). [2,1,0] pts 1,3 and [3,1,0] pt 1, 3 w each, ≤1e-83; θ-, h-, rep-independent. The plain product
D_B(k1)·S_B(w) is NOT a rotation (σ2/σ1 ~ 0.1–0.45). u(w) is NOT in the vacL(z) family (violates the p3 relation)
but IS a vacL for the Möbius-transformed twist y_k = (αz_k+β)/(z_k−w2) (exact, cell 12).
S_B(w) × the unchanged ⟨k2^{u/h}p Q12⟩ = ⟨V(w)|Σ_m(−1)^m C(λ2,m) k2^m t^{(1)}_{m,m}(θ1)⊗1|Ψ_B⟩: site 1 only again, c'_m = ρ_m,
w-independent (cells 13–14; [3,1,0] pt1, [3,2,0] pt1, [2,1,0] pt3, 3 w each, ≤1e-76). 6×6 flag formula also holds for [3,2,0].
CG count for it (`xxx-310-cg-test.wb` cells 9–10, z and w both untwisted): v''_new = e1ᵀ g_u T g_u⁻¹ reaches 3/6, 3/10, 6/10 —
identical to the det family (same zero pattern, both e1ᵀ n T n⁻¹, n ∈ N+); all three families combined still 5/6, 8/10, 9/10. Asked Paul how the two should be combined (not tried: a joint determinant).

**[4,2,0] check of the k1 det + k2 factor (2026-09-23)** — `L2-det-mathcalA.wb` at [4,2,0] pt1 (dim 729): D = C⟨v1⊗v1|Ψ⟩ to
1e-85 with C the closed form (prediction); c'_m = ρ_m = (1, 67/29, 1608/145), final residual 1e-74. Ledger rows 1–7 re-seeded.
CG count [4,2,0]² (`xxx-310-cg-test.wb`, saved outputs = [4,2,0]): 19 hw; old 15, det 6, w-vac 6, all 15 (finite twist also 15).
Missed [6,3,3], [5,4,3]×2, [4,4,4] = exactly ν3 > λ2 (predicted); v_new ≡ 0 exactly at μ3 ∈ {3,4}. Rule now holds on 4 reps.

**Docs (2026-09-23):** `XXX-CG/L2-det-mathcalA.tex/.pdf` (determinant + second factor), `xxx-210-cg-test.tex/.pdf`
(now with "determinant family", "t_{2,1} family" and "why ν3 ≤ λ2" sections), `xxx-310-cg-test.tex/.pdf` ([3,1,0]², [3,2,0]², the ν3 > λ2
rule + its derivation, reconstruction), `L2-localizing-transfer.tex/.pdf` (one-site action of t_{2,1}(θ1), scan, non-iteration).
Kernel state at end of 2026-09-23: xxx-210 cells 2–14 run for [2,1,0] (this overwrote the [4,2,0] C-suffix variables of xxx-310).
Kernel was restarted after these runs — rerun notebook cells before relying on kernel variables (resultsD, resultsF, jr*).

**2×2 determinant, same rotation on both sites (2026-09-23)** — `XXX-CG/L2-det-mathcalA.wb` (cell 1 = summary,
cell 9 = ledger `resultsD`). det[[⟨k1^{u/h}Q1⟩_α, ⟨k1^{u/h}uQ1⟩_α]]_{α=1,2}, ⟨f⟩_α = ∮ e^{2πi(α−1)u/h} f / ∏_{k=0}^{d} ν1(u+kh)
(two-site ν1, d = λ1−λ2), divided by (e2−e1) Q1(θ1+λ2h) Q1(θ2+λ2h) k1^{(θ1+θ2)/h+2λ2}, equals
C ⟨v1⊗v1|Ψ⟩ with C = (−1)^d / ((d!)² h^{2d} ∏_{n=−d}^{d}(θ12−nh)). 7 runs ([2,1,0] pts 1–3, [3,1,0] pts 1–2, [3,2,0] pts 1,3),
all orders in k1, ≤1e-85; C fitted on 4 runs, predicted on 3. Normalization is at the lowest pole θs+λ2h, not Q1(θs).
Times the unchanged L2-mathcalA second factor ⟨k2^{u/h}p(u)Q12⟩ (site-1 poles, p = site-2 factor): = C⟨v1⊗v1|Σ_m(−1)^m C(λ2,m) k2^m
t^{(1)}_{m,m}(θ1)⊗1|Ψ⟩ — k2 still acts on site 1 only; c'_m with v1⊗v1 = ρ_m, same as with v1⊗vacL (cells 10–11;
[2,1,0] pt1, [3,1,0] pt2, [3,2,0] pt3, ≤1e-80). Open: symmetric version (k2 on site 2 / both sites).
CG count for it (`xxx-210-cg-test.wb` cells 8–10, exact): undo g1 (site-2 rotation of v1, ∈ N+, fixes hw states) → v'_new =
(1−k1z1)² e1ᵀ(1−k2 g1 t̂ g1⁻¹), zero at (1,2,0),(0,2,1),(1,0,2),(0,1,2); reaches 3/6 hw states, combined with old 5/6 —
**singlet still missed**, ⟨v_tot⊗v1|singlet⟩ = 0 exactly. Hypothetical fix: t̂ on BOTH sites, ⟨v1t̂⊗v1t̂|singlet⟩ =
½p2²(z1−z2)²(z1−z3)² ≠ 0 → needs a site-2-pole ⟨k2 p Q12⟩ factor (not built yet; Paul: hold off on proposals).
[3,1,0]² (`xxx-310-cg-test.wb`, switch `lamSwitchC`, reproduces adjoint): 10 hw states; old 8/10, det 3/10, combined 8/10.
Missed [4,2,2], [3,3,2] (site-1 weights (1,1,2), (0,2,2)). v_new ≡ 0 at every site-1 weight with μ3 ≥ 2, as for the adjoint →
in both reps the missed irreps are exactly those with ν3 = 2.
[3,2,0]² (same notebook, cell 3 now has t11 and t22; saved outputs = [3,2,0]): 10 hw; old 9, det 6, combined 9. Missed only
[4,3,3] (site-1 weight (1,1,3)). v_new ≡ 0 exactly at μ3 = 3. Rule over all three reps: v_new ≡ 0 iff μ3 > λ2; missed irreps
= those with ν3 > λ2; det family never adds a state.

**CG test, adjoint (2026-09-23)** — `XXX-CG/xxx-210-cg-test.wb` + `.tex`/`.pdf`, all exact/symbolic.
Undo vacL globally: (v_tot ⊗ vacL)(g0⁻¹⊗g0⁻¹) = v_new ⊗ e1ᵀ, v_new = v_tot g0⁻¹. Untwist z_k = e^{φ_k ε},
leading order per weight space of V. Reached: 6/8 site-1 weights; **5/6 highest-weight states of V⊗V** —
[4,2,0], [3,3,0], [4,1,1], both [3,2,1]; **singlet [2,2,2] missed** (needs site-1 weight (0,1,2), where v_new ≡ 0
even at finite twist; (1,0,2) also ≡ 0). k1, k2 survive the limit, θ, h drop out; in the octet sector the
direction is (2k1−1, 1), so k1 alone separates the two octets. Lower weights come only from the rotated twist
G1 = g0 G g0⁻¹ (upper triangular, O(ε) off-diagonal, (1,3) entry = 0 by the p3 relation).

**Next:** reach the singlet (second t̂ insertion / λ2 ≥ 2 / other functionals); repeat the count for
[3,1,0]², [3,2,0]²; λ2=3 (needs t_{3,3}); different reps per site; derive c_m/κ2 rather than fit.

## >>> PREVIOUS (2026-09-19): twist-overlap reconstruction — the root-step rule <<<

**Read first:** `Paul/Mathematica/old/XXX/Experiments/SU3_Adjoint2_TwistReconstruction_L2.md` (full
reference), notebook `SU3_Adjoint2_TwistReconstruction_L2.wb` (31 cells, saved).

**Question (Paul).** Companion frame, L=2. Knowing `<Om(w)|tau_a(y,u)|Psi(z)>` for all vacuum
twists `w`, all transfer-matrix twists `y`, all `u`, `a = 0,1,2`, how much of `|Psi(z)>` is fixed?

**Answer.** Psi never enters: the solution set is `|Psi_A(z)> + W^perp` with
`W = span{<Om(w)|tau_a(y,u)}`, so the whole question is one rank, the same for all states.

```
W = V_{lamA+lamB}  (+)  everything ONE ROOT STEP away,  with multiplicity
```

| reps | dim | dim W | dim W^perp |
|---|---|---|---|
| `[2,1,0](x)[1,0,0]` | 24 | 24 | 0 |
| `[2,1,0]^2` | 64 | 63 | 1 (the singlet) |
| `[3,0,0]^2` | 100 | 63 | 37 |
| `[3,1,0](x)[2,1,0]` | 120 | 111 | 9 |
| `[3,1,0]^2` | 225 | 171 | 54 |

**Mechanism (derived, cells 29/31).** `qmin_{I,K}` with `|I|=|K|=a` transforms in
`Lam^a(C^3) (x) Lam^a(C^3)*`, which is `adj (+) 1` for a=1 AND a=2 (the `det` factors cancel), and
`1` for a=3 — hence `tau_3 ∝ Id`. Measured: the module generated under `ad` closes at dim 9 with
exactly 1 invariant, the invariant being the untwisted `Sum_i T_ii(u)`. And `<Om(w)|` sits in the
single component `V_{lamA+lamB}`. Wigner-Eckart then gives one root step per insertion.

**Fixes tested.** (i) `k` insertions reach `k` steps — measured 63/90/100 on `[3,0,0]^2`, so k=3
completes it; cost grows linearly with the rep. (ii) A general non-companion twist: NO help (the
module is twist-independent). (iii) The other five polarized vacua: NO help (all are
Cartan-component covectors). (iv) **One per-site rotation `<Om(w)|(exp(c E^(1)_31) (x) 1)`: ranks
63, 64, 64, 64 for 1,2,3,4 values of `c` — ONE nonzero `c` completes `adj^2`.** This is O(1) extra
data, and `rotations/` section 19 showed per-site rotations ARE Q-determinants, so the overlaps are
computable. **This is the bridge to the rotations thread.**

**NOT** the singlet story (`[3,1,0]^2` has no singlet yet loses 54) and **NOT** the
rectangular-vs-generic obstruction (`[3,0,0]` is rectangular and loses 37).

**Next, in order.** (1) Push the one-rotation-parameter fix to `[3,1,0]^2` and `[3,0,0]^2` — does
one `c` still suffice? (2) Does the *cheap* `E_32` rotation suffice, or is `E_21`/`E_31` needed
(sections 17/19 of rotations say those are hard as Q-determinants)? (3) Scaling of the seed ladder
for `[S,S,0]^2` vs SoV's `dim V`. Full list in section 11 of the `.md`.

**Traps.** `Join @@ Table[...,{a},{b}]` does not stack rows (bit 4x); hoist `om.O` out of inner
loops and use `prim` before `MatrixRank` or exact-rational rank will not converge at dim 100+;
`C_2` is NOT injective on partitions of fixed size (`(3,3,0)` and `(4,1,1)` both give 24 at
`|lam|=6`), so `lamOfC2` can silently mislabel.

## >>> CONTINUATION (2026-09-17): which rotations ⟨0|g are Q-determinants? <<<

**Read first:** `Paul/Mathematica/old/XXX/rotations/README.md`, then the `.tex`/`.pdf` there, then the
`.md` for exact numbers and the list of things that already failed.

**Where it stands.** `⟨0| = e₂₄`; the SoV covectors are fixed by `⟨x|Ψ⟩ = Ψ(x)/Ψ(x₀)` (all 24
states, 1e-73). The FSoV scalar product equals the physical one with `e₂₄` on both sides (2e-78).
Pairing it with the vacuum dual Q's at a second twist `w` gives the rotations
`(a,b,c) = (-w₁w₂, w₁, -w₁(w₂+w₃))`, hence the **surface `b = a(c-a)`**, which nothing deforms off
(225 minors, `e^{ku}` dressing, six permuted vacua). Small determinants with a pole-cancelling
polynomial give the `E₃₂` rotations — one site (local) or both (global, `c = -e^{kh}`, measure
`1/(ν₁ν₂) = (u-θ₁)×ladder`), and these are minors of the same scalar-product matrix.

**New (2026-09-18, cells 48-54).** Two of the README's open items are now closed.

- **§16, products of one-site brackets.** The per-site rotation `(e₈e^{c₁E₃₂})⊗(e₃e^{c₂E₃₂})` has
  SoV support on exactly the four `x₁₁` labels with the exact pattern
  `[[1,(47/28)c₂],[c₁-(19/28)c₂,c₁c₂]]`, `det = -19c₂(28c₁-47c₂)/784`. A product of two one-site
  brackets is rank one there, so it is a per-site rotation only on `c₂=0` (= §8) or
  `c₁=((t-2h)/(t-h))c₂` — **independent `(c₁,c₂)` is impossible**, the old conjecture was wrong.
  With `e^{k_su}` dressings `r_s=-e^{k_sh}` and the locus is `k₂-k₁=log(47/28)`. Also: the one-site
  bracket on the **fundamental** site is not a rotation at all (and not a product).
- **§17, Wronskian minors — two no-gos.** (i) Exact/Gröbner: an object built from one site's
  `Q`-values has `c_x` in the 8 site-2-base labels, and the only rotations supported there are the
  local `E₃₂` (`⟨b,a⟩`); no global rotation at all (`⟨c²,b,19a+9c⟩`). So **no local `E₂₁`/`E₃₁` from a
  one-site object**, whatever the measure. (ii) Numerical/complete: of the 9315 minors of the FSoV
  matrix with at least one `i=2` row (sizes 2–5), **none** is an SoV expansion — the best two-twist
  mismatch is `O(10^-1)`. Only the full 6×6 works, as the degree count (`Ψ_A` is degree 6) predicts.

- **§18, free row weights — the surface is rigid.** Writing the 6×6 with free row weights (30
  parameters, base = FSoV measure) and demanding both twist-independence (23 conditions) and
  rotation-ness (`δV ∈ span{V₀ρ(X)}`, 20 conditions) leaves a **14-dimensional** allowed space — twice
  the 7 obvious directions (6 row rescalings + `e^{ku}` dressing), which sit inside it. On all 14,
  `δ[b-a(c-a)] = 0` to 57–65 digits, and the induced `(δa,δb,δc)` has rank exactly **2** = dim of the
  surface. Free weights move the rotation *along* `b = a(c-a)`, reaching the same 2-parameter family as
  varying `w`, never off it. The response is exact, not finite-difference: a determinant is
  multilinear in its rows, so replacing row m by a single-pole row **is** `∂D/∂v^{(m)}_p`.

- **§19, independent per-site rotations — YES, they are Q-determinants.** §16's constraint
  `c₁/c₂ = (t-2h)/(t-h)` was a property of the *product* of two one-site brackets (rank-one coefficient
  matrix), not of bracket objects. A `2×2` **determinant** with free per-site weights `R,S` and a free
  insertion `g(u)` gives `c₁₂ = Sδ`, `c₂₁ = R(1-γ)`, `c₂₂ = RS(δ-γ)` — four knobs against three
  conditions — and hits **every** target `(c₁,c₂)` exactly: `(1,2)`, `(2,½)`, `(-3/7,5/2)`, `(1,0)`,
  `(0,1)`, all twist-independent, `σ₂/σ₁ ~ 1e-94` (genuine products), per-site fit residuals `~1e-94`.
  Even `1 ⊗ e^{c₂E₃₂}`, impossible for a product, works (`R = 19/28`, `S = 47/56`, cubic insertion).
  `γ` stays free: a one-parameter family of realisations per target. Note per-site rotations are **not**
  of the form `⟨0|g`, so they are neither on nor off the surface `b = a(c-a)` — that surface constrains
  global rotations only.

- **§20, Paul's product (global 2x2 det at `k`) x (one-site bracket at `k'`) — nearly right, fails by
  one power of `Q₁`.** Checked on all 24 states. With `X₁ = Q₁(θ₁+2h)/Q₁(θ₁+h)`,
  `X₂ = Q₁(θ₂+h)/Q₁(θ₂)`, the object is *exactly* (72-digit fit residual)
  `D = G_c(1 + c'X₁)`, `G_c = 1 + (9c/28)X₁ + (47c/28)X₂ + c²X₁X₂`, `c = -e^{kh}`, `c' = -e^{k'h}/z₁`,
  while `⟨0|g₁g₂` is `P = 1 + (9c/28 + c')X₁ + (47c/28)X₂ + c(c+c')X₁X₂`. So
  `D - P = (cc'/28)·X₁·(9X₁ + 19X₂ + 28c X₁X₂)`: the constant, `X₁`, `X₂` coefficients are exactly
  right, `X₁X₂` overshoots by `(19/28)cc'`, and `X₁²`, `X₁²X₂` are spurious. `X₁²` is not a
  wave-function ratio — the product puts `Q₁` at a site-1 `x₁₁` point twice (determinant + bracket)
  while `Ψ_A` has one `Q₁(x₁₁)` per site, and `1/(Q₁(θ₁+h)²Q₁(θ₂))` cannot cure that. Measured:
  rank`{c(z),c(z')} = 2` (not an SoV expansion), 14 labels, stabiliser 1, `O(1)` deviation after the
  best constant. **Exact correction:** `⟨0|g₁g₂ ↔ G_c + c'X₁(1 + cX₂)` — the dressing must hit only the
  site-2 leg. That is a sum of two products, i.e. an antisymmetric combination, which is what §19's
  determinant-with-insertion produces.

**Kernel hygiene — read this before running anything.** The kernel inherited at the start of
2026-09-18 held caches from two epochs (the previous session had restarted mid-way), so
`vacAtTwist[w]["qd"]` was inconsistent with the cached §7–§11 results: re-running the §10 object gave a
covector 0.1% off `⟨Ω^L_w|` at 71-digit precision, and recomputing `vacAtTwist` failed outright. After
a restart, **all of cells 1–57 were cold-run in order with no kernel messages**, and every stored
result reproduced exactly (§10 gives `(a,b,c) = (-35/6, 5/2, -263/42)` to 10 digits, residual 1e-70;
§16/§17 exact results bit-identical). Two minor counts in the §17 survey changed (nonzero minors
880 not 884 for shape (2,1), 64 not 73 for (2,2)); the docs carry the cold-kernel numbers. Do not
trust a warm kernel for a base point.

**Next.** Paul's stated want is *different rotations on each site*; §19 delivers that for `E₃₂`. The
natural continuation, and now the most promising item:

1. **Per-site `E₂₁`/`E₃₁`.** Free weights *plus* a free insertion is what unlocked §19; that
   combination has never been tried on an object with `Q₂` rows, which is what moving `x₂₁`, `x₂₂`
   needs. §17's no-gos do **not** cover it (they kill one-site objects and all FSoV minors, not a
   two-site free-weight Wronskian object). Start here.
2. **Finite** weight deformations of the `4×4` shape (3 rows `i=1`, 1 row `i=2`) against the explicit
   degree-4 `E₂₁` target in the §17 result cell — 16 parameters vs 46 conditions, and §18 shows
   first-order rigidity of the global surface, so the structure must conspire. Ask before grinding.
3. Vary the dual state *and* the weights together (each was killed separately, never jointly).
4. A different reference covector, i.e. a different `B`/SoV basis.
5. Derive `b = a(c-a)`; §18 makes this more pressing, since the relation is a property of the
   construction rather than of the measure.

## >>> PAUSED (2026-09-16): gl(3) adjoint L=2 — twisted match + untwisting <<<

**Read first:** `Paul/Mathematica/old/XXX/Experiments/SU3_L2_Untwisting.tex` (5-page PDF, the full
story) and the two companion `.md` files below. `summary.md` at the repo root is the *previous*
session's handoff (superseded by this section).

### Notebooks (Experiments/, all saved, none cold-run end-to-end in one go)

- `SU3_L2_Baxter_Twisted_Final.wb` + `.md`
  - **Twisted L=2 match.** For a unimodular twist `z=(5/2,7/3,6/35)` and `θ=(1/3,−8/57)`, the
    normalized physical scalar product `⟨L_A|R_B⟩/(⟨L_A|0⟩⟨0|R_B⟩)` equals
    `C_vac det N^{AB}/(c^ref_A Ψ^ref_B)` for **all 64×64 pairs, both a/b branches** (residual 1e-76).
    Cells 1–8 were cold-run.
  - **Rows.** `(α,i,j)`, α∈{1,2} (bare / `e^{2πiu/h}`), each bracket one contour around all 6 poles of
    `1/(ν1ν2ν3)`. Columns: `u^β×(−Q^jQ_i, Q^jQ_i^{[−2]}, Q^{j[2]}Q_i)`.
  - **Normalization.** `c^ref` = L=1 coefficient multiplied over θ1, θ2.
  - **Invariance.** All 36 row triples `{(a,c),(a,d),(b,c)}` work, with the normalization following the
    Q-content (`Q_a(θ+1)W(Q_a,Q_b)` and its dual) and a **triple-independent `C_vac`** (cells 9–10).
  - **Bug fixed.** `bra02`/`ket02` needed `Flatten`.
- `SU3_L2_Untwisting.wb` + `.md` + `.tex/.pdf`, with `z_j = e^{εφ_j}`, 200 digits:
  1. **Divergence.** `⟨0|R⟩ ~ ⟨L|0⟩ ~ ε^k`, `k = ht(λ_vac−λ_A) = 0,1,1,2,2,4`.
  2. **Precise relation (Appendix N of 2603.26617 generalized).** In the column basis H (3 surviving
     integrals of motion) / O (3 `O(φ)`), the bare rows are at infinity, and
     `det N_reg = −det[A1]_O · det[A2]_H`. The first factor is `V_λ(z)·E_λ`, with
     `E_λ = −(λ1−λ2+1)(λ2−λ3+1)/[(λ1−λ3+1)!(λ1−λ2)!(λ2−λ3)!]`. The second is the **untwisted**
     periodic 3×3 minor.
  3. **Result.** `lim V_{λvac}/V_{λA} · S(A,A) = B_{λA}/B_{λvac}`, with `B_λ = E_λ det𝒩/(c^ref Ψ^ref)`
     entirely in untwisted q's. `𝒩` is written out explicitly in the tex.
     - **Evidence:** matches the physical limit to 1e-19 at θ P1; `E_λ` independently validated at
       θ=(1/5,46/55) to 1e-13; **not derived** analytically.
  4. **Bare eigenvectors.**
     - `l_0 = (S⊗S) r_0`, `S = diag(1,1,1,2,6,9,4,9)`.
     - With `r_0` fixed by its |hw⟩⊗|extremal⟩ component, `B/B_vac = ‖r_0‖²_S/(N_R N_L)`, with
       `N_R N_L = 1, (1−t²)/4, (1−t²)/4, [(1−t²)(9−t²)/32]²` for 27, 10, 10̄, 1 (`t = θ1−θ2`).
     - Octets: the `B` formula is known in closed form (with √(4+5t²)), but no natural normalization
       has been found yet.
  5. **Separate limits.** `⟨0|R⟩/ε^k` and `⟨L|0⟩/ε^k` are each monomials in `φ_ij`. So
     `|R_untw⟩ = lim ∏(z_i−z_j)^{aR}|R_A⟩ = r_0/b_R` and `⟨L_untw| = lim ∏(z_i−z_j)^{aL}⟨L_A|` exist,
     and the object in (3) equals `⟨L_untw|R_untw⟩`. The split and `b_R`, `b_L` probably depend on
     a13, a23 and the branch.

### Next steps (Paul has not set the priority — ask)
1. **Derive `E_λ`** from the residue at infinity. `O1+O2` cancels at leading order, so the next order is needed.
2. **Octets.** Is `N_R N_L` an eigenvalue / matrix element of the adjoint⊗adjoint `Ř(θ1−θ2)`? That
   would also fix the octet normalization, since the octet block carries √(4+5t²).
3. **`b_R`, `b_L`:** closed form and a13/a23/branch dependence.
4. **Cold run** of `SU3_L2_Untwisting.wb`; then promotion to `Clean/`.

### Traps
- **Kernel state.** `SU3_L2_Untwisting.wb` cells 13/16/17/19 rebuild the chain at other θ.
  `buildAtU` does not rebuild `polesU`/`weightsU`. Run in order.
- **`0^0`.** `qpolyU` must be Horner, because a pole at u = 0 gives `0^0`.
- **Output truncation.** Wolfbook output truncates at `{`.
- **Toy model.** The "untwisted q × exact exponentials" toy gives the wrong coefficient; don't use it.
- **Large-denominator `Rationalize`** results for octets/singlet are unreliable.

## >>> PAUSED: closed form for N_untw — PARTIAL (gl(4), 2026-09-10) <<<

**Full briefing: `Paul/Mathematica/old/XXX/Experiments/TASK_N_untw_closed_form.md`** — read it in full.
Work in `Experiments/SU4_SoV_PerturbativeOverlap_L2.wb`, cells 19-23 (saved). Kernel restart: run
cells 1-15, then 19-22.

**Done.** `N_untw` = determinant of the **untwisted** 8x8 built from the untwisted Baxter
polynomials (`spMatU`/`detU`, cell 20; needed a new untwisted *dual* solver `UntwSolveD`). Verified
exact for `mu = {m,m,m,m}`, where that determinant is nonzero. For the vacua every Q is a pure
exponential, so the matrix is closed-form (`spMatE`, cell 21).

**Empirical, not derived.** For all other `mu` the untwisted determinant vanishes. What we have is an
exact ratio in the inhomogeneities only (`t = (theta1-theta2)/h`, `d = deg Q_12`):

```
N_untw(mu)/N_untw(vac) = c(S,P,d) [(S+t)(P-t)]^3 PROD_{k=1}^{d-1} (S+k+t)(P+k-t)
```

All six ratios come out **exactly rational**, verified at two different `theta2` to 30+ digits. But
this is **not** the Q-function answer Paul asked for: the rewriting via `Qbar_12(theta1+Sh)
Qbar_12(theta2+Ph)` leaves a non-Q factor `[(S+t)(P-t)]^2` and an unexplained `(S,P,d)`-dependent
constant, and `c(S,P,d)` is a fitted table.

**Ruled out:** the whole first-order reduced-minor / Appendix-N (N4) class — `corank(spMatU)` is
strictly less than the vanishing order for every degenerate state, so the degeneracy is not lifted at
first order. Also ruled out: pure Q-products (exhaustively), and both naive exponential dressings.

**Best lead:** for all six states the monic branch polynomials have degree `<= 1` (`q = 1` or `u+b`),
so the twisted Baxter system is **exactly solvable in the twists** — started and works (branch 1 with
`q=1` fixes `tc[3..5]` in closed form in `z_1..z_4`). That makes the confluent limit an analytic
Taylor extraction instead of a numerical fit. Section 5 of the task file has the route.

## Perturbative twist expansion — DONE (2026-09-10)

`Paul/Mathematica/old/XXX/Experiments/SU4_PerturbativeExpansion_L2.wb` (clean, self-contained, cold run
a few seconds) + `SU4_PerturbativeExpansion_L2.md` (reference + the bug list). Scratch predecessor
with the dead ends: `SU4_PerturbativeBaxter_L2.wb`.

**Result.** Two branches alone determine the whole transfer matrix, with a one-order lag: pushing to
`eps^N` with only `Q_1, Q_2` fixes all five `tc` at orders `1..N-1` (verified to N=8 against the
exact pseudovacuum `tau`); only the top order keeps free directions. The dual side `Q^4, Q^3` does
the same, identically.

**Structure.** `Q_k = z_k^(u/h) q_k`, `q_k` monic of FIXED degree `M_k`, exponential kept exact and
expanded by the operator, coefficients Laurent with branch `k` diverging as `eps^-(M_k-M_1)`. Dual:
`Q^j = Exp[-phi_j u] q^j`, `deg q^j = mu_j`. The staggered `Mhat_k = M_k+(k-1)` are the shadow of
this — the untwisted `P_k` are limits of COMBINATIONS, not of individual branches. For `k=2` the
regular combination is the divided difference: `b = -a`, `a ~ 1/eps`, nothing worse.

**Read the .md before touching this.** The expensive bug: accumulated solution rules must be
BACK-SUBSTITUTED, `acc = Join[acc /. new, new]`, or a later order falsely reports INCONSISTENT.
Also: `Solve` treats symbols outside `vars` as parameters and returns `{}` on generic failure, so
undetermined lower-order coefficients must go back into the unknown list.

## NEXT TASK (2026-09-10, Paul): leading small-twist behaviour of the scalar product

**Goal.** Feed the twist expansion we already have into the SoV scalar product and compute its
leading behaviour as the twist is switched off. The hoped-for form of the answer: the result
written in terms of **level-0 (untwisted) Q-functions**, multiplied by overall twist factors whose
powers are fixed by the **quantum numbers of the state** (i.e. by `mu`, equivalently by the degrees
`M_k = (S+P) - mu_k`).

**Start from these two files.**

- `Paul/Mathematica/old/XXX/Clean/SU4_SoV_ScalarProduct_L2.wb` — the finished scalar product, 9 cells,
  self-contained, runs from a cold kernel (~48 s for the orthogonality table). Companion
  `.md` (symbol reference) and `.tex`/`.pdf` (derivation, 6 pp). READ THE `.tex` FIRST: it is much
  cheaper than reconstructing the argument, and it lists the pitfalls.
- `Paul/Mathematica/old/XXX/Experiments/SU4_twisted_Baxter_L2.wb` — Sections 6-10 hold the existing
  twist-expansion machinery (`JointSolve`, the untwisted solver `UntwSolve`/`tc0U`, the
  eps-parametrisation, the divergence bookkeeping). Sections 11-12 are the Hodge duals and the
  scalar product, with the dead ends recorded so they are not retried.

**What is already known about the eps -> 0 structure** (Sections 6-7, verified):

- `q_k = O(eps^{-(M_k - M_1)})`; finer, `[u^j] q_k = O(eps^{-(M_k - max(j,M_1))})`. The individual
  branches diverge, but in the wedge the divergences cancel, leaving one power of eps per colliding
  pair — exactly the Vandermonde `V_I`.
- `W_I = V_I x (vacuum factor) x q~_I`, `V_I = Prod_{a<a'}(z_{i_a}-z_{i_a'})`, so for the triples
  entering the duals `V_Ibar = O(eps^3)`.
- `q_1234 = A(u) A(u-h)` in closed form, state- and twist-independent.
- The regular chain is `{1} subset {1,2} subset {1,2,3} subset {1,2,3,4}` when `M` is increasing
  (state becomes highest weight on untwisting); for other states the flag reverses.

**Why this should work at all.** In the scalar product

    M^{AB}[(alpha,i,j),(r,beta)] = Int_C du E1(u) Exp[2 Pi I u (alpha-1)/h] u^beta [q^j_A q^B_i]_r(u)

the measure `E1 = mu/nu1` and its pole ladder are **completely twist-independent**, and so is the
periodic factor. ALL twist dependence sits in the two Q-functions. So the small-twist limit is
purely a statement about `q^j_A` and `q_i^B`: the `q_i` diverge as above, the duals carry
`V_Ibar = O(eps^3)` and `zeta_Ibar^{(u-h)/h} -> 1`, and the eps-powers should organise into an
overall factor fixed by the `M_k`, leaving untwisted Q-functions inside the brackets.

**TRAP that will cost real time if missed.** The scalar product requires a **unimodular** twist,
`z1 z2 z3 z4 = 1`. The existing eps-parametrisation in Sections 6-10 uses
`psi = {1, 1/3, 1/5, 1/7}` with `z_j = 1 + eps psi_j`, whose product is NOT 1. Those test states
therefore CANNOT be fed into the scalar product as they stand. Use a parametrisation that is
unimodular order by order, e.g. `z_j = Exp[eps psi_j]` with `Sum_j psi_j = 0` (which is exact, not
just perturbative), and re-derive the divergence orders if the psi's change.

**Suggested route.**

1. Pick a unimodular `psi` with `Sum psi_j = 0` and rerun the exact solver at symbolic `eps` for a
   small state (rational, so results are exact) — the `(1,1)` weights `{2,2,0,0}` and `{2,1,1,0}`
   are the cheapest, `(2,2)` `{4,4,0,0}` vs `{4,3,1,0}` the next.
2. Expand each matrix entry of `spMatFSoV` in `eps` and read off the leading power per entry;
   check whether the powers are additive in the row/column labels (they should be, if the answer
   factorises).
3. Factor out the predicted overall power and check that what remains is built from the untwisted
   Q-functions (`UntwSolve`'s `P_k`, staggered degrees `Mhat_k = M_k + (k-1)`).
4. Then do the determinant: the norm `det M^AA` is the object with a chance of a clean closed form.

**Do not** re-derive the measure, the dual normalisation, the contour, or the `T_2` split — all
four are settled and each has a recorded failure mode. If something disagrees, suspect the twist
parametrisation (unimodularity) before suspecting Section 12.

**Practical notes.** Kernel `abort` reliably fails on long `Solve`s in this project — budget for a
`restart` (which loses state; reloading the Clean notebook is ~10 s, the Experiments notebook needs
cells 3-19, 24-28, 78, 82 rerun). Zero-tests must be `Simplify[...] === 0` on exact algebraic
numbers: the determinants are identically zero, so `N[]` merely exhausts `$MaxExtraPrecision` and
warns.

## State as of 2026-09-09 (twisted gl(4) L=2 Baxter solver + Wronskian structure — new notebook)

New notebook `Paul/Mathematica/old/XXX/Experiments/SU4_twisted_Baxter_L2.wb` (**uncommitted**, 37 cells,
sections 1–7). Built step by step at Paul's direction. It replaces item 1 of the 2026-09-02 plan
below: rather than wiring `Baxter_L2_XXX_SU4.wb` into the norm machinery, this is a fresh,
**purely functional** twisted solver — no GT patterns, no monodromy matrices, no eigensystems;
everything is a polynomial identity in `u`. Paul's instruction was explicit: ignore Matteo's
notebooks entirely, and postpone the untwisting/norm thread until the twisted side is understood.

### Conventions

`z[j]` are the primary symbols with `z[j] = Exp[h phi[j]]`, so `Q_j(u) = z_j^(u/h) q_j(u) =
Exp[phi_j u] q_j(u)`, `q_j` monic of degree `M_j = S+P-mu_j` (so `Sum_k M_k = 2(S+P)`). Sites
`theta1, theta2`, shift `h`, site reps `{S,S,0,0} x {P,P,0,0}`. Everything is rational in `z`;
`phi` enters only via a substitution. `mu` labels a gl(4) **weight** (the irrep, for a HWS).

### Established (verified, not conjecture)

- **Transfer-matrix ansatz: 12 free coefficients down to 5.** Two structural facts, both Paul's,
  both then derived and checked symbolically in `S,P,theta,h,z`:
  - `tau3 = nu_1(u-2h) nu_4(u) * (quadratic, leading chi_3)`. Because `nu_1=nu_2` and `nu_3=nu_4`
    for rectangular reps, every term of `tau3 = Sum_{i<j<k} z_i z_j z_k nu_i(u-2h)nu_j(u-h)nu_k(u)`
    shares that factor. Untwisted it reproduces the `trivfactor` of `baxter_SU4_L2_XXX.wb`.
  - **Subleading coefficient, universal:**
    `[u^(2a-1)] tau_a = -a (theta1+theta2+(a-1)h) chi_a - h Xi_a`,
    `Xi_a = Sum_{|I|=a} z_I Sum_{i in I} n_i`, `n = mu` the global Cartan weight.
    Compactly `Xi_a = Sum_i n_i z_i d(chi_a)/d(z_i) =: D_n chi_a` (verified exactly, `a=0..4`), so
    the whole subleading term is a twist prefactor times `chi_a`, plus `h` times the weights acting
    on `chi_a`. Checks: `a=0` gives 0; `a=4` reproduces the `u^7` coefficient of the independently
    known quantum determinant; `a=1,2,3` reproduce the pseudovacuum eigenvalue **identically in u**.
  - Leftover unknowns: `tc[1]` (in `tau1`), `tc[2..4]` (in `tau2`), `tc[5]` (the constant in
    `tau3`'s quadratic).
  - **This formula is not written down anywhere else in the repo** — checked. `paper/xxxCG.tex` is
    gl(2)/single-site only; `Baxter_L2_XXX_SU4.wb`'s
    `alphaCoeffs` are numeric `1/u` coefficients used only to get `Mdeg`; `baxter_SU4_L2_XXX.wb`'s
    `asympSol` is untwisted-only and not in universal form; the structure-constants paper gives
    `tau±` as Q-ratios (C11)–(C12), not a large-`u` subleading. The general-`L` extension
    `[u^(aL-1)]tau_a = -a(Theta + L(a-1)h/2) chi_a - h D_n chi_a`, `Theta = Sum_alpha theta_alpha`,
    follows from the same derivation but is **untested** (only `L=2` verified).
- **The degree formula is derived, not assumed.** For an arbitrary `q` of degree `M`, the top
  coefficient of the Baxter residual is `charPoly[z_k] = 0` automatically, and the next one is
  `h z_k Prod_{j!=k}(z_k - z_j) * (mu_k - S - P + M)` — zero iff `M_k = S+P-mu_k`, given
  `Sum mu = 2(S+P)`. It degenerates exactly when twists collide, which is the untwisting limit.
- **`TwistedSolve[S,P,mu,sub]`** solves the joint bilinear system (5 `tc` + `2(S+P)` q-coefficients;
  `28+2(S+P)` equations, overdetermined by 23). One branch is *not* enough — for `mu_top` at (1,1)
  branch 1 alone leaves a 2-parameter family; all four together are unique.
- **Number of solutions = weight multiplicity**, checked against the product character
  (`wtMultProd`). At (2,2): 1, 2, 3, 8, 12, 22. Expected — the twist breaks su(4), so every state
  of that weight shares the same `M_k`.
- **Validation on the six permuted vacua** `n = pi({S+P,S+P,0,0})` (multiplicity 1, hence rational
  and exactly checkable — Paul's suggestion, and it is enough). At `(1,1),(2,1),(2,2),(3,3),(5,4)`:
  30 states, all unique, Baxter residual **exactly 0** on all four branches, and `tau` matching an
  independently computed pseudovacuum prediction. 0.03–0.2 s each. Repeated at (1,1) with **no
  substitution at all** (`theta1,theta2,h,z_j` symbolic): same 1/0/0, ~3 min.
  - **Shift ordering** for the permuted vacua is *not* the natural index order: the quantum-minor
    shifts follow the Borel adapted to `pi`, `ord = Join[A, complement A]` with
    `A = {i : n_i = S+P}`. Found empirically (searching all 24 permutations against the solver's
    own output selects exactly these). Consistent with `tau4` being ordering-independent, as it
    must be, being central.

### Wronskians (Section 6) — the normalisation `N_I`

With `Q_j = z_j^(u/h) q_j`, the `k x k` Casoratian factorises with the exponentials pulled out:

    Q_I(u) = (Prod_{i in I} z_i)^(u/h) (Prod z_i)^(-(k-1)) W_I(u),
    W_I(u) = det_{a,b} [ z_{i_a}^(k-b) q_{i_a}(u-(b-1)h) ]

`W_I` has degree `Sum_{i in I} M_i` and leading coefficient **exactly the Vandermonde**
`V_I = Prod_{a<a'}(z_{i_a}-z_{i_a'})` (immediate: the `q`'s are monic). Hence

    N_I = V_I / (Prod_{i in I} z_i)^(|I|-1),      q_I = W_I / V_I    (monic, degree Sum M_i)

and `q_I` is regular as `phi -> 0` **exactly when `I` is a lower set for the degrees `M`** (the
`|I|` branches with smallest `M_i`, ties interchangeable). For a state that becomes highest weight
on untwisting `M` is increasing, so the regular chain is `{1} ⊂ {1,2} ⊂ {1,2,3} ⊂ {1,2,3,4}` —
i.e. `q_1, q_12, q_123, q_1234`, as Paul stated. For `I={1,2}`: `N_12 = (z_1-z_2)/(z_1 z_2)`.

- **Divergence of the individual `q_k`:** `q_k = O(eps^(-(M_k - M_1)))`, finer pattern
  `[u^j] q_k = O(eps^(-(M_k - max(j,M_1))))`. Mechanism: the untwisted equation admits the staggered
  degrees `Mhat_k = M_k + (k-1)`, so a monic `q_k` of degree `M_k` must tend to an untwisted
  solution of degree `<= M_k`, and for `k>=2` the only one available is branch 1 of degree `M_1`.
  In the wedge the divergences cancel (the limit of `q_k` is *proportional to* `q_1`), leaving one
  power of `eps` per colliding pair — exactly `V_I`.
- **`q_1234` in closed form, state- and twist-independent:** from `C(u+h) = A_4(u) C(u)` plus
  `nu_2=nu_1, nu_3=nu_4`, `W(u+h)/W(u) = nu_4(u)nu_4(u-h)/(nu_1(u)nu_1(u-h))`, which telescopes:

      q_1234(u) = A(u) A(u-h),   A(u) = Prod_{j=1..S}(u-theta1-jh) Prod_{j=1..P}(u-theta2-jh)

  degree `2(S+P) = Sum M_k`. Verified exactly (residual 0) on 20 states.

### Vacuum factors and the divergence count (Section 7)

Eq. (C6) of the structure-constants paper (`References/2603.26617v1.pdf`) reads
`Q_ij = W(Q_i,Q_j)`, `Q_theta Q_ijk = W(Q_i,Q_j,Q_k)`, `Q_1234 = W(Q_1..Q_4)` with
`Q_theta = Prod_alpha (u-theta_alpha)` of degree `L`. Here `L -> S+P` and `Q_theta -> A`:

    W_12   = V_12 q~_12                deg q~_12   = M1+M2
    W_123  = V_123 * A(u) * q~_123     deg q~_123  = M1+M2+M3 - (S+P)
    W_1234 = V_1234 * A(u) A(u-h)      deg q~_1234 = 0

Verified exactly (`PolynomialRemainder` = 0, quotient monic of exactly that degree) on all six
vacuum weights at (2,1) and (3,1), the pseudovacuum at four `(S,P)`, and a genuine HWS. These
degrees are also `Sum Mhat - k(k-1)/2` with `Mhat_k = M_k+(k-1)` after removing the `A`'s — the
vacuum-factor count and the staggering are the same statement.

**Uniform count for `p_k`.** Divide the `k x k` relation by `V_{1..k-1}`; this strips every
divergence and leaves a `(k-1)`-th order difference operator `O_k` with *regular* coefficients:

    O_k[q_k] = ( Prod_{i<k}(z_i - z_k) ) * (vacuum factor) * q~_{1..k}

with (i) source of `eps`-order `k-1` and degree `Sum_{i<=k} M_i`; (ii) `O_k` at `eps=0` mapping
degree `d` to `Sum_{i<k}M_i + d - (k-1)`, its polynomial kernel being exactly the untwisted
`span{P_1..P_{k-1}}` of degrees `Mhat_1..Mhat_{k-1}`; (iii) the tower `g_0 ~ q_1` (degree `M_1`),
`deg g_j = M_1 + j`. The source enters at `j = p+k-1`, so

    Sum_{i<k}M_i + (M_1+p+k-1) - (k-1) = Sum_{i<=k}M_i   ==>   p_k = M_k - M_1

uniformly in `k`. For `k=2` this is a fully independent first-order/`eps^1` derivation (two ways:
monicity, and the inhomogeneity). Input (ii) is **verified explicitly for `k=3` only** — on state A
the kernel sits at exactly `d = Mhat_1, Mhat_2 = 0, 2` and `deg O_3[f] = M1+M2+d-2` elsewhere; on
the rational vacua `O_3` degenerates to the plain second difference `f(u)-2f(u-h)+f(u-2h)`. For
`k=4` the pattern is asserted, not checked.

### Section 11 (2026-09-09): Hodge-dual `Q^j` — the prefactor `F(u)` is fitted

Convention (asymmetric, Section 6's `Wronsk`; the paper's symmetric (C5) is not used):
`F(u) Q^4(u-h) = det[Q_1(u), Q_2(u-h), Q_3(u-2h)]`, `Q^j = (-1)^j Q_{comp j}`, `F` fixed by
state-independence with the dual polynomial monic. With `I = complement of j`,
`zeta_I = prod z_i`, `V_I` the Vandermonde on `I`, and `A` the Section 7 vacuum factor:

    Q^j(v) = (-1)^j zeta_I^(v/h) q^j(v),  q^j(v) = (W_I/V_I/A)(v+h)  monic, deg q^j = mu_j
    F_j(u) = (-1)^j (V_I/zeta_I) A(u),    A(u) = Prod_{k=1..S}(u-theta1-kh) Prod_{k=1..P}(u-theta2-kh)

So `F` is `A(u)` up to a constant: degree `S+P`, `S` factors with `theta1` and `P` with `theta2`;
the paper's `(u-theta1)(u-theta2)` is the `S=P=1` case. Verified 16/16 on `(2,1) {3,3,0,0}`,
`(1,1) {2,1,1,0}`, `(2,2) {4,3,1,0}` at `h=1`, plus `(2,1) {3,3,0,0}` at `h=3` with different
`theta`. Cells 77-80; symbols `QQfun`, `WasymQ`, `qDual`, `QDual`, `Ffit`, `Gfit`, `numH`.
The pairwise `det[Q^i(u),Q^j(u-h)]` relation holds up to `Ztot^(u/h)` times a pair-dependent,
state-independent constant (normalisation bookkeeping).

Still open: dual Baxter equation `B† Q^j = 0` unchecked; `Q^j` not yet expanded in the twist.
NEW TRAP: a definition with its RHS split over several lines is silently truncated at the first
line in this notebook — keep every RHS on one line.

### Section 12 (2026-09-09): O, O-dagger and the orthogonality identity — VERIFIED

Conventions: `f^[±2n] = f(u ± n h)`; rotated Gamma with `Gam[F(u+h)]/Gam[F(u)] = F(u)`.

- **`O` (12.38), 4th order** = this notebook's Baxter equation divided by `Prod_{m=0..3} nu1(u-mh)`
  with `G(u)=Q(u+h)` and `T_{a,1}=tau_a`. Same operator for all four branches — the `z_k` in
  `BaxRes` are only the stripped exponential. `O q_i^[2] = 0` verified for all `i`.
- **`O†` (12.39), 4th order, `nu3 -> nu4`**: `O† q^j = 0` verified for all `j`, but only for the
  UNCENTRED dual `q^j(u) = QDual[j](u-h)`. Section 11's `QDual` has a `+h` centring baked in.
- **Measure**: `mu(u) = 1/Prod_alpha Prod_{n=0}^{S_alpha-1} (u - theta_alpha - n h)` — `S_alpha`
  one-sided poles per site at `theta_alpha + n h`. `[2,2,0,0]^2`:
  `1/((u-th1)(u-th1-h)(u-th2)(u-th2-h))`.
- **Identity (12.42)** `0 = Contour_C mu q^j_A (O^A - O^B) q_i^{B[2]}` holds exactly — the contour
  must encircle ALL finite poles of the integrand, not just `mu`'s (the `O` coefficients add
  `nu1(u-mh)` denominators: 6 poles at (1,1), 8 at (2,2)). Verified at (1,1) and (2,2), exact
  non-perturbative twisted solutions, unimodular twist `numU`. Both terms vanish separately;
  plain brackets are nonzero.

Cells 81-84; symbols `OpO`, `OpOd`, `qhat`, `qdualhat`, `mMeas`, `contourSum`, `relCheck`, `idTest`.

DEAD ENDS (do not retry): the `theta_alpha +- (2l-1)h/2` ladder with the paper's 8x8 determinant
(zero at (1,1) by luck, fails at (2,1),(2,2)); per-pole functionals with `beta=1..S+P`; fitting
`(a,b)` in `theta-a h, theta, theta+b h` by determinant vanishing — the root-find runs away to
`a~10^3` where the determinant underflows trivially.

NEXT: the determinant/orthogonality packaging on top of the verified identity.

### Section 12.3-12.5 (2026-09-10): all four integrand terms on ONE measure

Trivial factors: `T_1 = t_1`, `T_3 = nu4 nu1^[-4] t_3`, `T_4 = qdet` (cancels). `T_2` has none but
is degree 4; split per state `T_2 = nu1(u+p1) t21 + nu4(u+p2) t22` (deg 2 each; one-parameter
ambiguity `t21 -> t21 + c nu4`, fixed by `t21(0)=0`, cancels in the sum).

Effective measures are contiguous ladders `1/(u-theta_alpha-nh)`:
`E1 = mu/nu1` (n=0..S_alpha), `E2 = mu/(nu1 nu1^[-2])` (n=0..S_alpha+1),
`E3 = mu nu4/(nu1 nu1^[-2])` (n=1..S_alpha+1). Key identity **`E3(u+h) = E1(u)`**.

`p1=-h`, `p2=0` are forced (generically) by demanding the leftover ladder stay contiguous on both
sites. Then `u -> u+h` in `I_22` and `I_3` puts everything on `E1`:

    0 = Int du E1(u) [ -Dt1(u) qA^j(u) qB_i(u) + Dt21(u) qA^j(u) qB_i(u-h)
                       + Dt22^[2](u) qA^j(u+h) qB_i(u) - Dt3^[2](u) qA^j(u+h) qB_i(u-h) ]

The four shift combinations `(0,0),(0,-2),(+2,0),(+2,-2)` are exactly the `a,b=+-1` pattern of the
FSoV determinant, with degree-<=2 coefficient functions. Verified EXACTLY (Simplify) at (1,1) and
(2,2), all four (i,j). Cells 84-86; symbols `splitT2`, `E1meas`, `alignedTerms`.

NEXT: this is the natural starting point for the determinant packaging.

### PROMOTED TO CLEAN (2026-09-10)

`Paul/Mathematica/old/XXX/Clean/SU4_SoV_ScalarProduct_L2.{wb,md,tex,pdf}` — the final scalar product
only, 9 cells, self-contained (verified by running it top to bottom in a FRESH kernel; the
orthogonality table takes ~48 s). `.md` is the symbol reference, `.tex`/`.pdf` (6 pp) the
derivation notes. The Experiments notebook keeps the full history and the dead ends.

### Section 12.6-12.7 (2026-09-10): SoV SCALAR PRODUCT WORKS — orthogonality at (1,1) and (2,2)

Gauge for the T2 split (Paul's): `[u^2]t21 = [u^2]t22 = (1/2)[u^4]T_2 = chi_2/2`. Unique, and since
`chi_2` is state-independent both `u^2` terms cancel in the difference => every `Delta t` is linear,
`Delta t = DeltaI_2 u + DeltaI_1`. (The earlier `t21(0)=0` gauge is WRONG - do not use.)

So the aligned identity is a linear relation among 8 IoM differences, giving an 8x8:

    M^{AB}[(alpha,i,j),(r,beta)] = Int_C du E1(u) Exp[2 Pi I u (alpha-1)/h] u^beta [q^j_A q^B_i]_r(u)

rows `(alpha,i,j)`: alpha=1,2 (alpha=2 inserts the h-periodic `Exp[2 Pi I u/h]`, inert under the
alignment shifts), `i in {1,2}`, `j in {3,4}`; columns `(r,beta)`: the four shift patterns
`(0,0),(0,-h),(+h,0),(+h,-h)` against `beta = 0,1`; poles = the E1 ladder `theta_alpha + n h`,
`n=0..S_alpha`.

RESULT: `det M^{AB} = 0` for A != B, `det M^{AA} != 0`. Verified EXACTLY (Simplify) at (1,1) and
(2,2), both orderings, INCLUDING same-weight pairs (`{2,1,1,0}` #1|#2; `{4,3,1,0}` #1|#2) - the
sharp test. 8/8 cross determinants zero, 5/5 norms nonzero.

Cells 87-88; symbols `ladE1`, `perFac`, `brkE1`, `spMatFSoV`, `orthTest`.

### Test states (for reuse)

`z_j = 1 + eps psi_j`, `psi = {1,1/3,1/5,1/7}`, `theta1=1/3, theta2=1/7, h=1` — rational in `eps`,
so `TwistedSolve` runs with `eps` fully **symbolic** and results expand exactly; no numerics, no
state tracking. States used: **A** = (1,1) `mu={2,1,1,0}`, `M={0,1,1,2}`; **B** = (1,1) lowest
weight `{0,0,2,2}`, `M={2,2,0,0}` (flag reversed — regular chain is `{4}⊂{3,4}⊂...`); **C** = (1,1)
singlet `{1,1,1,1}`, `M={1,1,1,1}` (control: all `M_k` equal, so nothing diverges); **D** = (3,1)
`mu={3,3,1,1}`, `M={1,1,3,3}` (the only one with `M_1>0` *and* real divergence).

Identifying which solution of a degenerate weight sector becomes highest weight: solve the
**untwisted** problem first (`UntwSolve`, staggered degrees `Mhat`, standard gauge fixing) — it is
unique and rational — then pick the twisted branch converging to that `tc0`, validated by
`dist(tc,tc0) ~ eps`. Cheaper and far more robust than degree-spectrum matching on 7 algebraic
branches.

### Traps hit this session (all cost real time)

- `Limit[]` applied to a whole rational expression is **unreliable** on the algebraic
  (multiplicity>1) states — it returned *different branches on repeated runs* of identical code.
  Take the limit coefficient-by-coefficient instead.
- A helper that `Module`-localises the spectral variable returns expressions in a scoped `x$nnn`;
  the caller's `CoefficientList` then silently reads **every coefficient as degree 0**. Pass the
  variable in.
- Detecting admissible untwisted degrees by nullspace **dimension** silently misses a jump when two
  admissible degrees are adjacent (gave `{0,2,4,5}` instead of `{0,1,4,5}` for irrep `{2,2,0,0}`).
  Test **monic solvability** instead; that agrees with `degQ` of `baxter_SU4_L2_XXX.wb`.
- `RootReduce` is exact but far too slow to sweep with; bare `PossibleZeroQ` falls back to an
  undecided numeric test and warns. Hence the rigorous checks are done on the **rational** (vacuum)
  states, where `Expand[Together[...]]` collapses to a literal 0.
- Dividing by `V ~ eps^6` is catastrophic numerically: the low coefficients of `q_1234` read as `0`
  at small `eps` (visible in the accuracy annotations). Confirm such things exactly, not by `N[]`.
- Kernel `abort` again did not take on a long `Solve`; a `restart` was needed (loses state) —
  same as recorded on 2026-09-02.

### Not done

- **The actual perturbative expansion of the polynomial coefficients** — Paul's original step 2. He
  interrupted it to understand the Wronskian first, and that is now understood, so this is the
  immediate next step. The `q_j` stay `Exp[phi_j u] x` polynomial of fixed degree `M_j`
  non-perturbatively; only the coefficients are expanded, and they need Laurent series (orders
  `eps^(-(M_k-M_1))`, now known in advance).
- Nothing is wired to `SU4_norm_untwisting_L2_general.wb` / the FSoV norm yet.
- No scan over reps (item 3 of the 2026-09-02 plan).
- Write-up: the subleading-`tau` formula and the `N_I`/vacuum-factor ledger are worth recording, but
  `paper/xxxCG.tex` is currently gl(2)/single-site only, so adding a gl(4) L=2 section is a scope
  decision Paul has not made. A `Clean/` promotion + `wolfbook-docs-agent` reference doc is the
  other option.

## State as of 2026-09-02 (untwisting the FSoV norm — gl(4) L=2; two notebooks, results below)

Thread: **how to untwist the norm analytically**, started from Appendix N of the new reference
`References/2603.26617v1.pdf` (Bargheer–Bercini–Lefundes–Ryan, "Structure Constants from Q-Systems
and SoV"; Paul is an author). Also new: `References/SU4_twisted_any_L_Qs.nb` — Matteo/Paul's twisted
Q-solver for any L, 6-dim `[1,1,0,0]` site rep.

Two notebooks, both under `Paul/Mathematica/old/XXX/Experiments/`:

- `SU4_norm_untwisting_L2.wb` — LR-symmetric branch (`z1 z4 = z2 z3 = 1`), staged untwisting
  `phi1 -> 0` then `phi2 -> 0`. Ported solver validated against the paper's (P1).
- `SU4_norm_untwisting_L2_general.wb` — general unimodular twists, independent `Q^3, Q^4`,
  high-precision NSolve route. Validated twice: reproduces (P1) exactly, and reproduces the
  first notebook's `B = -0.0274316528219221` to all digits.

### Established (verified, not conjecture)

- `B` (eq. (12)) is **exactly invariant** under `GL(2) x GL(2)` on `(Q1,Q2)` / `(Q^3,Q^4)`,
  rescalings included: the row mixing gives `(det H)^{2L}`, the normalisation (15) gives
  `(det H)^{-2L}`. So the rotation `Q2 -> Q2 + r Q1` is **pure gauge** — it cannot change `B`,
  only whether the limit is computable entrywise. Checked numerically to 20+ digits.
- **Untwisted degrees**: `lambdahat_j = lambda_j + (j-1)`. Verified on two states with different
  predictions: `{0,1,0}` (`lambda={0,1,1,2}`) gives `lambdahat={0,2}` with
  `Q1hat=1, Q2hat=u^2/2+1/3`; `{1,2,1}` (`lambda={1,1,1,1}`) gives `lambdahat={1,2}` with
  `Q1hat=u, Q2hat=u^2+1/3`. Mechanism: subtracting the collision leaves the next order of the
  `Exp[phi u]` expansion, raising the degree by one.
- **The rotation**: fixed by a quadratic with a *double* root. For `{0,1,0}`, coefficients
  `2I, -4I, 2I` i.e. `2I(c-1)^2`, giving `r = 1/phi2`, i.e. subtract the divergent multiple of
  `Q1` (there `p2 = u - Cot[phi2]`, so `Q2 -> -Q1/phi2`).
- **Why the rotation is mandatory**: for `{0,1,0}` the naive entrywise min-cost predicts
  `phi2^(-2)` while the truth is `phi2^12` — 14 orders of cancellation.
- **Everything is a contraction with a universal moment table.**
  `M[(alpha,j,k),(beta,a,b)] = Sum_n c_n mu[n,alpha]`, `mu[n,alpha] = <u^n>_{L,alpha}`.
  For `L=2`: `mu[n,1]` is nonzero **iff n is odd and n >= 2L-1** (parity selection rule, proved by
  closing at infinity: `mu[n,1] = -I Binomial[-L,k] 4^(-k)`, `k=(n-2L+1)/2`); `mu[n,2]` has the
  single zero `mu[1,2]=0`. This explains both previously-unexplained patterns: the `alpha=2`
  vanishing (row `(2,1,4)` of `{0,1,0}` *is* `mu[beta-1,2]`) and the `a=-b` zeros in block `(2,3)`
  (the surviving `u^3` coefficient is `I(a+b)/4`).
- **Overall exponent**: along a ray, `det` (with (12)'s normalisation) vanishes as `eps^nu` with
  `nu = -Sum_{j<k}(lambda_j - lambda_k) = -Sum_j (n+1-2j) lambda_j`. Verified: `{0,1,0}` gives
  `nu=6` (measured ratios 5.963, 5.991, 5.998, 5.999) against `Sum lambda_jk = -6`; `{1,2,1}` gives
  `0` against `0`. This is exactly Appendix L's (L1) prefactor `Prod z_jk^(lambda_jk)`, so
  **(L1)'s `B`, not (12)'s, is the object with a uniform limit** — one determinant, not one per state.
- **Closed results.** `{0,1,0}`: `B -> -32 phi1^2 phi2^4`, where `-32 = -4 * 8` and `8` is the
  reduced minor built from untwisted Qs (`Q12hat = -I u`). `{1,2,1}`: corank 0, so the plain
  untwisted 8x8 determinant already works and gives `B = 7776 = 6^5` exactly, matching the twisted
  limit (7477.7, 7700.4, 7757.0, 7771.3, 7774.8 -> 7776).

### Caveat and what's missing

(L1) makes `B` finite, but when the corank is positive you have `0 * Infinity`: the untwisted
determinant is literally 0 and the prefactor literally infinite. To *evaluate* the answer from
untwisted Qs you still need the reduced minor, whose size is the corank.

**Only two highest-weight states so far** (`{0,1,0}` = `[1,0,1]`, `{1,2,1}` = zero-weight
subspace, 6 solutions = 2+3+1). At `L=2` every state with non-negative Dynkin labels is
automatically LR-symmetric (`6 (x) 6 = 20' + 15 + 1`, all self-conjugate ⇒ `M = Nn`), so the
6-dim site rep cannot supply more data. Hence the plan below.

### Plan agreed for the next session (Paul, 2026-09-02)

> **Status as of 2026-09-09:** item 1 is done, but *differently* — instead of wiring
> `Baxter_L2_XXX_SU4.wb` into the norm machinery, a fresh functional solver was built
> (`SU4_twisted_Baxter_L2.wb`, see the 2026-09-09 section above). Item 2 is the immediate
> next step and is now much better set up (the `eps`-parametrisation makes exact symbolic
> solves possible, and the divergence orders are known in advance). Item 3 is untouched.

1. **Generalise the determinant to `[S,S,0,0] x [P,P,0,0]` site reps.** *Mostly already done*:
   `Experiments/Baxter_L2_XXX_SU4.wb` (2026-07-19 entry below) is already a **twisted** gl(4) L=2
   Baxter solver for two independent **generic Young-diagram** reps per site, with four distinct
   diagonal twists and the rep-dependent vacuum factors
   `nu[lam1,lam2][i][u] = (u - th1 - h lam1[[i]])(u - th2 - h lam2[[i]])`,
   `qdetT = chi4 Prod_i nu[i][u-(4-i)h]`. Full-sweep verified: 152 states, 5 rep pairs,
   worstTQ = worstCasoratian = 0 exactly. Same Q convention as here (`Q_k = z[k]^(u/h) q_k`).
   So `[S,S,0,0] x [P,P,0,0]` is a special case and arbitrary Young diagrams come free — the task
   is wiring it to the norm/bracket machinery in `SU4_norm_untwisting_L2_general.wb`, not deriving
   new taus. (The other notebook, `baxter_SU4_L2_XXX.wb`, has Casimir-form taus for `[S,S,0,0]` but
   **untwisted** — not the one to use here.)
   **Precision trap** recorded in that entry: hardcoded `N[..., 30]` literals inside
   `TauEigensystem4`'s `rSets`/`Kmix` silently capped the whole downstream chain regardless of
   `z[k]`'s precision. Grep for bare precision literals before trusting a scan.
2. **Build a perturbative-in-twist solver.** Every failure this session was NSolve
   ill-conditioning plus state-tracking jumps as `eps -> 0`, because the solutions genuinely
   coalesce there. A series solver in `eps` avoids it, and is structurally right: the leading
   order *is* the degenerate untwisted problem, so this is degenerate perturbation theory and
   `r` should drop out of the order-`eps^0` solvability condition instead of being fitted.
3. **Scan over reps** to deduce the untwisting pattern. Per state record: `lambda_j`;
   `lambdahat_j` (test `lambda_j + (j-1)`); corank and which entries vanish; `nu` against
   `-Sum_{j<k} lambda_jk`. Compute `mu[n,alpha]` once per `L` and check whether the odd-`n`
   rule and `mu[1,2]=0` persist at other `L`.

### Performance notes (matter for the scan)

- Do **not** use `Residue` for the FSoV bracket. Precompute `mu[n,alpha]` symbolically once and
  contract with coefficient lists: ~0.1 s per 8x8 norm instead of minutes.
- Symbolic `Solve` on the general (non-LR) system returns branches with `Sqrt` of perfect squares
  that resist `PowerExpand`; that is why the reference notebook uses `NSolve`. The perturbative
  solver should sidestep this entirely.
- Kernel `abort` did not take on a long `NSolve` sweep; a `restart` was needed (and loses state).

## State as of 2026-08-27 (su(2) SoV localization — teaching notebook, complete; FSoV link pending)

`Mathematica/su2/su2-localization.wb` (uncommitted, along with `Mathematica/su2/su2-localization.md`).
Paul wrote this notebook specifically to teach the localization structure of the SoV basis; Claude
extended it this session with the mirror construction. **Read `Mathematica/su2/su2-localization.md`
first** — it has the conventions, symbol tables, the verified-check table, and the traps.

Established and symbolically verified (L=2 gl(2), companion twist, `Lam[1]=Lam[2]={3,0}`): the SoV
basis built up from the minimal vacuum factorizes on the `mu2=0` edge, `x[{mu,0}] = x1[{mu}] (x)
vac2`, exactly; the mirror basis built down from the maximal `B`-eigenvector factorizes on the other
edge, `xb[{0,nu}] = vacb1 (x) x2b[{nu}]`, also exactly. Both towers build the same basis (`B` has 16
distinct eigenvalues). 7 of the 16 states factorize, not all. Localization is a property of the
vacuum, not of the operator.

Localization is stronger than a vacuum statement: it holds on the whole orbit of the vacuum under
the **global** (diagonal) su(2) rotation, `<0| Exp[r EE[{2,1}]] . t[{1,mu}][w1] == <0| Exp[r
EE[{2,1}]] . (t1[{1,mu}][w1] (x) 1)`. Globality is essential — independent per-site parameters break
it, restored exactly when they agree.

**Next step, from Paul — to be set up explicitly in a later session.** The strategy (sketched, not
yet implemented or verified; full write-up in `su2-localization.md` §7): globally rotated states
`<Omega(r)|` have `<Omega(r)|Psi>` known exactly in Q-functions; localization means a transfer
matrix acting on such a state hits **site 1 only**; certain **combinations of `t`'s sum up to global
rotations**, so such a combination rotates site 1 alone and produces `<Omega(r'), Omega(r)|` with
independent parameters per site; and since `t|Psi> = tau|Psi>` the same overlap is `tau` times the
known one. Net: `<Omega(r'), Omega(r)|Psi>` in Q-functions — the two-parameter family CG
coefficients are read off from. The crux is that the input orbit is one-parameter and the target is
two-parameter, and localization is what breaks the diagonal constraint. Don't start building this
ahead of that session.

Standing constraint restated: **this project only ever deals with two sites.** Unequal reps behave
identically, so the equal-rep choice is convenience, not a special case.

## State as of 2026-08-19 (Matteo's gl(4) Q-space/normalisation notebook imported — read/analyzed only, NOT yet computed on)

New material from Matteo, committed this session (`de97e03`), not yet acted on beyond a structural
read-through. Four files, three of which are the *same* underlying notebook in different formats:

- `Paul/Python/SU(4)_Paul_Normalisation_EXPLICIT.nb` — Matteo's original Mathematica notebook.
- `Paul/Python/SU(4)_Paul_Normalisation_EXPLICIT.wb` — wolfbook copy of the same (open this with
  `wolfbook-builder`/the wolfbook MCP tools to compute live, cell-by-cell, on either machine).
- `Paul/Python/SU4_Paul_Normalisation_EXPLICIT.wl` — a flattened, **directly `Get`-able** script
  version (presentation/output cells stripped; first line says so). Good for `grep`/quick reads or
  a plain `wolframscript`/kernel run; NOT meant to be run top-to-bottom as one script — it's full of
  bare unassigned "example of use" expressions (`Sol=SolveBaxAll[...]` etc.) meant for selective
  cell evaluation, and later sections redefine earlier-named symbols in the 1-site block.
- `Paul/Matteo_Files/SU(4)_Paul_Normalisation_SEND.wl` — **not directly runnable**: this is a raw
  Mathematica "notebook-as-text" export where every cell, *including Input cells*, is wrapped in
  `(* ::Input:: *) (*code*)` comment markers (presumably literally what Matteo emailed/sent). Do not
  `Get[]` it expecting code to run — strip the cell-marker comments first, or just use the `.wb`/`.nb`
  above, which already have live code.

**What the notebook computes** (from reading `SU4_Paul_Normalisation_EXPLICIT.wl` top to bottom, not
yet cross-checked with Paul or numerically re-verified beyond what the notebook's own test cells
already assert): a Gaudin-side two-site (`w1=0,w2=1`) gl(4) FSoV construction —
1. **Ψ-space**: GT-basis generators (`ps/dim/J0/Jp/Jm/JJ`, same machinery as
   `Gaudin/Clean/Hamiltonians_SU4.wb`), Lax connection `LIndex`, explicit `H1..H4` and dual
   `H1†..H3†` Hamiltonians (regulated traces of powers of `L`), and `EigenR/EigenL` — right/left
   joint eigenvectors obtained by matching against the **Q-side** spectrum rather than a bare
   `Eigensystem` call.
2. **Q-space**: `h[rep1,rep2][k,1/2]` (residues at `w1,w2`, closed gl(4)-weight polynomials) plus
   twist elementary-symmetric pieces build the 4th-order Baxter operator `BaxOp`; `SolveBaxAll`
   solves all four Baxter equations for the monic `Q[i]` coefficients **and** the leftover
   Hamiltonian unknowns simultaneously (spectrum + Q-functions in one solve). A dual sector `Q†`
   is built two ways — directly via a dual Baxter equation, and via `Hodge` (Wronskian of the other
   three `Q`'s) — with a "Check Q Hodge" block cross-checking the two agree.
3. **Renormalised scalar product**: `OrthogRen` — a 12×12 residue-moment determinant (poles at
   `z=0,1`) divided by a `Normalisation` built from boundary limits of `Q†_1`, `W(Q†_1,Q†_2)`,
   `Q_4`; the analogous 1-site pieces are `OrthogRen1/OrthogRen2`.
4. **Ψ-vector normalisation (the "EXPLICIT" part)**: `VacR/VacL` build the vacuum state as an
   ordered `MatrixExp` product of the six raising/lowering `EMAll[a,b]` generators with parameters
   `α_ab,β_ab`; `αr[12,24,34]`/`βr[...]` are then **closed-form rational functions of the twists
   `φ`** that fix those parameters (leaving 3 free: `α13,α14,α23`), giving renormalised eigenvectors
   `ReEigenR/ReEigenL` in the same normalisation as the Q-side `OrthogRen`.
5. **`K` vectors** (final ~230 lines): `K22B`/`K22C` — explicit overlap vectors between a 1-site
   and a 2-site factor (SU(2)⊗SU(3)-split), built the same way as `VacR/VacL` but with extra
   `k3,k4`-dependent rational dressing. `TESTK22B/C` + the "Test reconstruction" cells verify, for
   `{1,1,0,0}×{1,1,0,0}` only, that `K·Ψ` (solved for by expanding in `k2,k3,k4`) reproduces the
   `ReEigenR` vector obtained independently via step 4 — i.e. this is a from-scratch construction and
   check of exactly the kind of normalised-eigenvector ↔ Q-function overlap identity that the
   `Gaudin/Clean/Hamiltonians_SU4.wb`/CG-overlap program on the XXX side has been chasing (see
   memories `project_cg_overlap_dictionary`, `reference_fsov_determinant_normalization`) — but done
   independently by Matteo on the Gaudin/FSoV side. **Not yet compared against that program in any
   way; that comparison is the obvious first thing to do before extending this notebook.**

**Two likely bugs, found by inspection, NOT yet fixed (present identically in `.wb`/`.nb`/`.wl`,
so not a conversion artifact) — check before relying on either symbol:**
- `VacR1[rep1_,β12_,...,β34_] := VacR1[rep1,β23,β12,β13,β14,β23,β24,β34] = ...` (1-site vacuum,
  `.wb` line 2235/3829) — the memo-pattern LHS takes 8 arguments (extra leading `β23`, and `β23`
  appears twice) against the 7-arg definition head, so it never matches an actual 7-arg call.
  Effect: memoization is simply dead (recomputes + leaves junk defs behind each call); results are
  still numerically correct, just uncached.
- `K22C` (`.wb` line 4088) — in the first `KroneckerProduct` factor (built on `UPL1[rep1]`), only
  the `E13` exponential uses `EM[...,rep1]`; the `E14`, `E23`, `E24` ones use `EM[...,rep2]` instead.
  Only evaluates without error when `dim[rep1]==dim[rep2]` (true for the notebook's only test case,
  `{1,1,0,0}×{1,1,0,0}`) — will break or silently misbehave for `rep1≠rep2` (e.g. `{p,p,0,0}×
  {q,q,0,0}`, `p≠q`). Compare against the parallel, apparently-correct `K22B` just above it in the
  notebook when fixing.

**To pick this up on a different machine:** `git pull` gets all four files (pushed in `de97e03`).
Needs a Wolfram kernel; use `wolfbook-builder`/wolfbook MCP tools against the `.wb` for live,
checkpointed cell-by-cell computation (matches how `Gaudin/Clean/` work is normally done in this
repo) rather than trying to batch-run the flattened `.wl`. No non-core-WL package dependencies were
seen anywhere in the four files.

**Next steps (none of this has been discussed with Paul yet — treat as an inbox item, not a plan):**
- Fix the two bugs above if/when extending past `{1,1,0,0}×{1,1,0,0}`.
- Work out how this notebook's `OrthogRen`/`ReEigenR`/`K22B`/`K22C` machinery relates to the
  existing XXX-side CG-overlap dictionary and FSoV-determinant-normalisation work — likely the
  actual point of importing it, but unconfirmed.
- This currently sits under `Paul/Python/`, which per `CLAUDE.md` is meant for the sympy Python
  cross-check port, not Mathematica notebooks — a Mathematica `.wb`/`.nb` pair there is a layout
  mismatch. Don't relocate without asking Paul first (same rule already applied to
  `Hamiltonians_SU4.wb`'s Gaudin/XXX placement).

## State as of 2026-07-22 (C++ CGC reference benchmark set up — `Paul/Cpp/AlexCGC/`, UNCOMMITTED)

Separate thread from the Baxter/Gaudin work below: set up the **state-of-the-art SU(N) CGC code**
(Alex, Kalus, Huckleberry, von Delft, arXiv:1009.0437 — the Gelfand–Tsetlin pattern-calculus C++
program from that paper's appendix) as a **reference benchmark** against Paul's Q-function CGC
generating function (`Paul/Mathematica/old/Gaudin/Clean/CGGenFun.wb`). New directory `Paul/Cpp/AlexCGC/`
— a C++ cross-check sibling to `Paul/Python/`. **Everything uncommitted** (Paul's choice: place,
don't commit; a compiled binary is machine-specific — rebuild from source elsewhere).

**Toolchain:** needs LAPACK/BLAS. `sudo apt-get install liblapack-dev libblas-dev` was run this
session (installed, 3.10.0). Build all binaries with `bash Paul/Cpp/AlexCGC/build.sh` (plain g++).

**What's there:**
- `ClebschGordan.cpp` → `clebsch` — the authors' original source (downloaded, unmodified); an
  interactive menu calculator (menu 5/6 = CGCs for S⊗S′). `hwcgc.sh` filters its output to just the
  highest-weight state (always the largest index `Q(M″)=dim(S″)`).
- `cgc_hw_bench.cpp` → `cgc_hw_bench` — modified copy that times ONLY the highest-weight-state solve
  (`compute_highest_weight_coeffs`, the raising-op nullspace/SVD) and SKIPS building descendants
  (3 clearly-commented `// added` edits: a `highest_weight_only` ctor flag, a scope timer, a new
  benchmark `main`). Verified `max|hw_only − full| = 0` vs the full run in every case. Modes:
  `cgc_hw_bench N S… S′… S″… [reps] [verify]` (per-irrep timing); `… decompose N S… S′…` (fast μ
  list); `machine` flag = one parseable RESULT line.
- `test.sh S [reps]` — for SU(4), times the highest-weight CGCs of EVERY μ in `[S,S,0,0]⊗[S,S,0,0]`
  (summary table only, no coefficient dumps). `reps` default 51 (report min-over-reps, for the fast
  small cases); use `reps=1` for large S (each solve already runs for seconds).
- `results/` — `hw_timings.txt` (full per-μ tables S=1..8), `summary.csv`, `gen_results.sh`
  (regenerates both). Machine: i5-1345U, 16 GB.

**The benchmark result (the point of all this):** C++ cost is dominated by the SINGLET (0,0,0,0)
— its highest-weight linear system is the largest — and total-over-all-μ hw-solve grows ~×5 per
step in S: S4=0.057, S5=0.355, S6=2.49, S7=15.3, S8=80.2 s. It also has a **memory wall**: the
singlet's dense SVD matrix is ≈ O(dim(ν)³) ≈ tens of GB by S=10 (dim[10,10,0,0]=1716), so the C++
**cannot even run S=10** on a 16 GB laptop. Paul's GF (build-for-all-irreps time) is ~**flat**:
S4=0.063, S6=0.139, S8=0.171, **S10=0.35 s**. ⇒ crossover ~S=4; GF ~470× faster at S=8, and at
S=10 it's "0.35 s vs won't-run" — a scaling win in **both time and memory**, not a constant factor.

**Open / next stage (Paul):** the GF timing is currently **build-only** — the explicit CGCs still
need to be EXTRACTED from the generating function (Paul's stated next implementation stage). That is
the one piece needed to make the head-to-head deliverable identical to the C++ (which outputs
explicit coefficients); then drop GF per-μ times into a second CSV beside the C++ `results/` for a
direct side-by-side. See memory `project_gf_cgc_vs_cpp_benchmark`.

## State as of 2026-07-22 (gl(4) L=2 untwisted τ1,τ2,τ3 ALL DONE exactly via Baxter asymptotics; Gaudin h→0 limit + standalone `GaudinSolver.wb` — DONE)

New notebook `Paul/Mathematica/old/XXX/Experiments/baxter_SU4_L2_XXX.wb` (Paul's own fresh build,
independent of `Baxter_L2_XXX_SU4.wb`) closes out the τ2 gap left open in the 2026-07-20 entry
below, via a genuinely different method — large-`u` asymptotic expansion of the untwisted Baxter
(TQ) equation, rather than operator-level Casimir fits. **This supersedes that entry's "τ2 —
PARTIAL" status: all three transfer matrices are now closed-form in Cas2(μ), Cas4(μ).**

### The mechanism (Paul's idea, worked through live)

- `n_k = μ_k` on the HWS (mode number = weight component) — verified directly against the main
  notebook's `nvec` inversion.
- Untwisted (`z_k→1`) Q-degree formula: `deg q_k = (S+P+3) − μ̂_k = (S+P)−μ_k+(k−1)`, `μ̂_k=μ_k+4−k`
  the *same* shifted weights that appear in the paper's Casimir generating function (eq 2.13–2.15,
  `G4[u]=Prod_i(u−h(μ_i+4−i))`). The `+(k−1)` staggering is forced, not just plausible: the
  twisted-case `Mdeg` formula's denominator is (up to `z³`) the derivative of the twist
  characteristic polynomial, which has a 4-fold zero at `z=1` — L'Hôpital-resolving that
  degeneracy is exactly what produces the staggered degrees.
- Key structural fact: `Σ_j(−1)^j C(4,j) j^m = 0` for `m=0,1,2,3` (a 4th finite difference) ⇒
  large-`u` orders `1/u,1/u²,1/u³` of the Baxter equation are **state-blind** (independent of the
  Q-data), and since the 4 Q-branches share one transfer matrix but have distinct leading powers
  `p_k`, every power of `p` must vanish separately. Solving that system pins 5 of τ1/τ2/τ3's 7
  unknown `u`-coefficients (`a[1],a[2],b[1],c[2],b[3]`) in terms of `S,P,θ` and one remaining
  scalar `b[2]` — leaving only `b[2]` (τ2's `u²` coeff, carries Cas2) and `d[2]` (τ2's constant,
  carries Cas4) as genuine Casimir-dependent data.
- Pseudovacuum anchor (`μ_top=[S+P,S+P,0,0]`, all τ's = elementary-symmetric-in-ν vacuum
  eigenvalues) gives `b2vac,c2vac,d2vac` independently — **cross-check passed exactly** (vacuum
  `c[2]` matches the asymptotically-predicted `c[2]`, residual 0).

### Result (τ1, τ2, τ3 all closed, all anchored at ΔCas = Cas(μ)−Cas(μ_top))

```
τ1(u) = τ1vac(u) + (h²/2)ΔCas2
τ2(u) = τ2vac(u) + (h²/8)ΔCas2·[τ1vac(u)+τ1vac(u−h)] + h⁴(⅛ΔCas2²−¼ΔCas4) + h⁴(1+2S+2P)/4·ΔCas2
τ3(u) = trivfactor(u)·[ 4u²−4u(θ1+θ2)+4θ1θ2−2hu(S+P+4) + b3vac + (h²/2)ΔCas2 ]
```

**The bug and how it was actually caught:** the first τ2 guess (fusion form
`(h²/8)ΔCas2·[τ1vac(u)+τ1vac(u−h)]` matching u²,u¹ exactly, plus the "obviously right" quartic
`h⁴(⅛ΔCas2²−¼ΔCas4)`) was missing the last `h⁴(1+2S+2P)/4·ΔCas2` term. Found by building the actual
polynomial Baxter solver and testing it: for `μ={2,1,1,0}@(1,1)` (`deg q1=0`), a constant `q1`
solved the equation at `μ_top` but **not** off-top — residual `=5·ν1(u−2h)ν1(u−3h)` exactly (a
pure constant shift), which pinned the missing term via the actual solvability condition, and
cross-checked against 27 irreps symbolically in `θ,h` (all exact) plus 4 more irreps with
`deg q1≥1` to break a degeneracy where `Cas2²` and `Cas4` weren't separable. **Lesson: build the
solver early — it's also the most direct correctness check.**

### Baxter polynomial solver + gauge fixing (in `baxter_SU4_L2_XXX.wb`)

`QsolveG[μ,S,P]` (note: `(μ,S,P)` order in *this* notebook — see naming-convention note below)
solves the untwisted 4th-order Baxter (difference) equation for the four monic Q-polynomials,
**gauge-fixed**: since `deg q1<...<deg q4` strictly, any lower-degree solution can be added to a
higher one (`q_k → q_k+c·q_j`, `j<k`); removed by forcing the coefficient of `u^(deg q_j)` in `q_k`
to vanish for every `j<k`. Verified exact (symbolic in `θ,h`) over full decompositions of every
`(S,P)` up to `(3,2)`, and numerically (rational `θ1=1/3,θ2=1/7,h=1`) over **all 21 irreps at (5,5)**
and **all 45 irreps at (8,8)** — every `q_k` solves the equation exactly, timing flat at
~0.02–0.05s/irrep **independent of irrep dimension** (largest block tested: `dim=69825`, same cost
as `dim=1`) since the solver only ever touches degree-≤19 polynomials, never the representation
itself.

Also built (Schur-polynomial machinery, reused verbatim in `GaudinSolver.wb` below): `irreps[S,P]`
— full `[S,S,0,0]⊗[P,P,0,0]` tensor-product decomposition via `s_λ·s_μ=Σ c^ν_{λμ} s_ν`, peeling off
the lex-highest monomial each step. Confirmed multiplicity-free throughout, dimension sums exact
(`Σdim(μ)=dim([S,S,0,0])²`), box-count (`Total[μ]=2S+2P`) conserved on every irrep as expected.

### Gaudin h→0 limit (derived in `baxter_SU4_L2_XXX.wb`, then built standalone)

Writing `Q=Q0+h·Q1+...`, orders `h⁰..h³` of the cleared Baxter equation vanish identically (same
alternating-binomial identity as above); order `h⁴` gives a **closed 4th-order ODE on Q0 alone**
(no `Q1,Q2` needed) — the gl(4) 2-site Gaudin oper, entirely in `Cas2(μ),Cas4(μ),S,P,θ1,θ2`.
Leading coefficient `8(u−θ1)²(u−θ2)²` (double singular points at the sites, as expected).
Validated: the ODE's polynomial-solution degrees match `degQ[μ,S,P]` exactly (and *only* those
degrees admit a solution) across every irrep of every `(S,P)` tested, symmetric and asymmetric.

**New standalone notebook `Paul/Mathematica/old/XXX/Experiments/GaudinSolver.wb`** — deliberately
minimal, gl(4) irrep/Casimir/degree machinery + the hardcoded `gaudinODE` only, **no** XXX transfer
matrices or Baxter difference-equation machinery (by design, per Paul's request). Built across
several rounds of naming/convention changes (all re-verified after each):
- Sites renamed `θ1,θ2 → w1,w2`; spectral variable renamed `u → z`.
- **Argument convention, applied everywhere: `(S,P,μ)`** — `degQ`, `GaudinSolve`, `BuildQs`,
  `verifyGaudin`, `CheckAllIrreps` all take `S,P` first, `μ` (or the irrep list) last. (Note: this
  is the *opposite* order from `QsolveG[μ,S,P]` in `baxter_SU4_L2_XXX.wb` — the two notebooks are
  not argument-compatible; don't copy-paste calls between them without swapping.)
- `GaudinSolve[S,P,μ]` returns **monic coefficient lists** (the practical/canonical form); `BuildQs
  [S,P,μ]` is a thin wrapper returning the same four Q's as **explicit polynomials in z** (kept
  purely for display/sharing, never used for internal verification/timing).
- `CheckAllIrreps[S,P,irrList]` takes a **precomputed** irrep list rather than recomputing
  `irreps[S,P]` internally — `irreps[8,8]` alone costs ~1.6s vs ~0.13s to actually solve all 45
  irreps, so precompute-once-reuse is the efficient pattern. Verified at `(8,8)` (45 irreps) and
  `(3,2)` (6 irreps).
- **Verification convention, hard-won:** always reconstruct a coefficient list into a testable
  function via `qf[x_]:=Sum[cf[[j+1]] x^j,{j,0,m}]` (SetDelayed), never
  `Dot[coeffs,z^Range[...]]` wrapped in a pure `Function` — the latter has a genuine WL edge case
  where symbolic `Derivative` of a length-1 (`deg=0`) Dot-built function spuriously returns a
  wrapped `{0}` instead of `0`, causing false verification failures. Cost a full debugging round
  before being isolated.

**Assistant arithmetic error, corrected:** an earlier summary wrongly claimed
`Σdim(μ)=8721²=76,055,441` for the `(8,8)` sweep — `8721` is `dimGL4[{16,16,0,0}]` (the *top*
output irrep), not `dimGL4[{8,8,0,0}]` (the actual input irrep, `=825`). Correct identity:
`Σdim(μ)=825²=680,625` — verified directly, matches the sweep. The sweep/solver results themselves
were never wrong, only that one squared-dimension claim (caught by Paul).

### Next natural steps

- **τ2 is now closed** — this directly unblocks the CG-overlap program the whole `z→1` limit was
  meant to feed (see the 2026-07-19/20 entries below for the original motivation). Natural next
  step: extract CG coefficients using the now-complete Casimir-form transfer matrices.
- Promote `baxter_SU4_L2_XXX.wb`'s Casimir-form τ1,τ2,τ3 + solver to `Clean/` once Paul is
  satisfied (still in `Experiments/`).
- `GaudinSolver.wb` is a clean, minimal, independently-reusable artifact — a natural next step is
  extending it toward the Gaudin-side CG-overlap analog, or documenting it via
  `wolfbook-docs-agent` once/if it's promoted to `Clean/`.
- Kernel-state loss happened repeatedly this session (multiple `KERNEL RESTART` events between
  turns) — always re-run cells top-to-bottom (`?symbolName` to check `DownValues` first) rather
  than trust that a previous turn's definitions are still live.

## State as of 2026-07-20 (gl(4) L=2 untwisted Casimir eigenvalues — τ1,τ3 DONE exactly, τ2 partial, in Experiments)

Follow-on to the twisted build below: new section appended to `Baxter_L2_XXX_SU4.wb` (now 70 cells,
35 onward), `z[k]→1` (untwisted/physical point), for the two-row-rectangular rep family
`λ1=[S,S,0,0]`, `λ2=[P,P,0,0]`. Goal: closed-form eigenvalues `Λ_a(u)` of the untwisted transfer
matrices `τ1,τ2,τ3` via su(4) Casimir invariants of the irreps `μ` in `λ1⊗λ2=⊕N_μμ`. Design/plan:
`docs/superpowers/specs/2026-07-19-su4-l2-untwisted-casimir-eigenvalues-design.md`,
`docs/superpowers/plans/2026-07-19-su4-l2-untwisted-casimir-eigenvalues.md` (the plan's original 5
tasks are DONE; everything below is follow-on work past the plan, done inline/live in the notebook
in direct collaboration with Paul, not via a formal plan). Subagents lacked wolfbook tools this
session (again) — a wolfbook-builder dispatch for the τ2/τ3 investigation partially worked (tools
were attached that one time) and its findings were independently re-verified from primitives.

**Exact-arithmetic simplification:** at `z=1` (vs. the twisted case's irrational `Zeta`-based twists)
every quantity is exact rational — no precision-cascade concerns anywhere in this section.

**Rep family is fully multiplicity-free** across the whole `(S,P)∈{(1,1),(2,1),(1,2),(2,2),(3,1)}`
sweep tested (up to `dim=400`), so every `τ_a` acts as a scalar on each irrep block (Schur's lemma) —
this is what makes closed forms in terms of Casimirs meaningful here at all.

### τ1 — DONE exactly

`Λ_1(u) = 4(u−θ1)(u−θ2) − 2P·h·(u−θ1) − 2S·h·(u−θ2) + (h²/2)(C2(μ) − C2(λ1) − C2(λ2))` (cell 48),
verified over the full sweep, `worst ~ 1e-25`. All `μ`-dependence enters only through the quadratic
Casimir eigenvalue `C2(μ)`, affine with slope `h²/2`. (This *resolved* an earlier false negative: a
naive fit of "each `u^j` coeff of `Λ_1` = `A_j+B_j·C2(μ)` **pooled across all `(S,P)`**" failed for
`j=0,1` — that was a pooling artifact, not a real falsification; the correct per-`(S,P)` formula
above is exact. See cells 45–47 for the superseded pooled-fit attempt, kept for the record with an
updated conclusion cell.)

### τ3 — DONE exactly

`Λ_3(u) = trivfactor[S,P]·( s3scalar[S,P](u)·Id + (h²/2)(C2(μ)−C2(λ1)−C2(λ2)) )` (cell 65), where
`trivfactor[S,P] = (u−θ1)(u−θ2)(u−θ1−h(S+2))(u−θ2−h(P+2))` (Paul's own factorization of τ3's four
trivial zeros, cells 58–63) and `s3scalar[S,P](u)` is an explicit state-independent degree-2
polynomial. **The `h²` term is exactly the same quadratic-Casimir combination as τ1's** — verified
as an exact operator identity on `(1,1),(2,1),(1,2)` including the `S↔P`-asymmetric case. So τ3
carries no higher Casimir, only `C2` — the `Λ³ℂ⁴≅(Λ¹ℂ⁴)*⊗det` self-duality made concrete.

### τ2 — PARTIAL (the self-dual middle case, genuinely harder)

Writing `Λ_2(u) = Σ_k[α_k(S,P)+β_k(S,P)·C2(μ)]u^k`, with `u⁰` additionally carrying
`+⅛C2(μ)²−¼C4(μ)` (universal coefficients, verified): **`α4=6,β4=0` (leading); `α3(S,P)=
−6h(S+P+2)−12(θ1+θ2), β3=0`; `β2=h²` (universal constant); `α2(S,P)` fully resolved (messier
polynomial, exact); `β1(S,P)=−(h³/2)(S+P+2)−h²(θ1+θ2)`** — all verified exactly against `tgCoeffs`
across 5 `(S,P)` points via direct operator-level matching (`M_k = α_k·Id+β_k·C2op`, solved exactly,
not curve-fit) — see cells 49–57 (the original agent-built investigation, confirming `τ2` uniquely
reaches `C4` only in its constant term) plus a from-scratch re-derivation in the live session after a
kernel restart (search for `solveABsym` in the transcript/notebook history).

**Still open: `α0(S,P)` and `α1(S,P)`'s full `(S,P)`-dependence.** These need cubic (`α1`) and
quartic (`α0`) polynomials in `(S,P)` — more unknowns than the 5 available `(S,P)` data points can
resolve. Getting more points is expensive (`(4,1)`→dim 630, `(3,2)`→dim 1000, symbolic-in-`θ,h`
computation at that size is impractical). An operator-level derivation was attempted instead: `T[i,j]`
was expanded exactly in terms of `Ee1,Ee2` building blocks (verified to reproduce `T` and `tg[2,1]`
bit-for-bit — see below) with the idea of reading off `α1,α0` structurally rather than by curve-fitting,
but this did not pay off within the session (a residual-decomposition test against `{X1,X2,Y1,Y2,C2}`
at a single `(S,P)` is underdetermined and not actually informative — would need to be done
simultaneously across multiple `(S,P)` to mean anything). **Promising unexplored lead, added by Paul
independently in cells 58–61**: `t1vac,t2vac,t3vac` — explicit vacuum-eigenvalue formulas built from
`L1=L2=(u−θ1−hS)(u−θ2−hP)`, `L3=L4=(u−θ1)(u−θ2)` — `t2vac` sums all 6 pairs `Li(u)Lj(u−h)` in a
specific descending-index order. An earlier attempt this session to derive the τ2 vacuum eigenvalue
via `Σ_{i<j}ν_i(u)ν_j(u−h)` (ascending-index order) did **not** match `tg[2,1]` on the HWS (nonzero
residual on pairs like `{1,3}`) — Paul's `t2vac` uses the opposite index order, which may be exactly
the fix; reconciling this against the actual `qmin`/`tg[2,1]` construction is the natural next step
and could resolve `α1,α0` cheaply (vacuum eigenvalues need no large-matrix computation at all).

### General Casimir operator formulas (not just eigenvalues — the operators themselves)

Cells 65–68 (note: cell 67, `X1/X2/Y1/Y2` definitions, was a **bugfix** — those operators were
originally only ever created via scratch `evaluateExpression` calls and never actually saved to a
notebook cell; a fresh-kernel run would have failed before the fix). General master formula (any `n`,
verified `n=2,3,4`): `Cn_op[λ1,λ2] = Σ_{pat∈{1,2}^n} patternTerm[pat,λ1,λ2]`, where each of the `2ⁿ`
site-assignment patterns collapses (via `(A⊗B)(C⊗D)=(AC)⊗(BD)` applied repeatedly) to `[ordered
product of site-1 factors]⊗[ordered product of site-2 factors]`, cyclic indices tying both chains
together. Reduced explicit forms with named cross terms:
```
C2op = (Csite[λ1]+Csite[λ2])Id + 2·Omega2                                    [1 cross term]
C3op = (C3site[λ1]+C3site[λ2])Id + 2X1+2X2+Y1+Y2                              [4 cross terms]
C4op = (C4site[λ1]+C4site[λ2])Id + 2Z1+2Z2+2Z3+2Z4+Z5+2Z6+Z7+2Z8               [8 cross terms]
```
Cross-term-type count grows 1→4→8 for `n=2,3,4` (matrix products aren't cyclic, only traces are, so
patterns that look like rotations of each other are usually genuinely different operators) — itself a
notable structural fact. Every formula verified as an **exact** operator identity (`==`, not
tolerance) on asymmetric `λ1≠λ2` pairs (avoiding a `λ1=λ2` coincidence hiding real distinctions), plus
`C3op`'s master formula also checked outside the `[S,S,0,0]` family (adjoint⊗fund) to confirm
generality. Plain operator definition (cell 70, self-contained): `C2=Σ_ab E_ab.E_ba`,
`C3=Σ_abc E_ab.E_bc.E_ca`, `C4=Σ_abcd E_ab.E_bc.E_cd.E_da`, `E_ab:=EE[λ1,λ2][a,b]`.

### General Casimir eigenvalue formulas (cell 69, self-contained, no other cell needed)

Single-site rep `[x,x,0,0]`: `Cₙ(x)=2x(x+2)ⁿ⁻¹` — `C2c[x]=2x(x+2)`, `C3c[x]=2x(x+2)²`,
`C4c[x]=2x(x+2)³`, verified `x=1..5` against the operators. Combined two-site:
`C2[S,P]=2S(S+2)+2P(P+2)`, `C3[S,P]=2S(S+2)²+2P(P+2)²`, `C4[S,P]=2S(S+2)³+2P(P+2)³`.

**Real bugs found and fixed this session** (beyond the earlier three): `X1/X2/Y1/Y2` never persisted
to a cell (see above); `T[i,j]` off-diagonal action on the HWS genuinely does vanish both directions
(`T[i,j][u].HWS=0` for `i≠j`, confirmed identically in symbolic `u`) — useful building block, not a
bug, but non-obvious and worth remembering; a hand-derived vacuum-formula attempt for τ2 had an
index/shift-order bug not yet isolated (see τ2 section above).

**Next natural steps:** reconcile Paul's `t2vac` (cell 61) against `tg[2,1]` on the HWS to get τ2's
vacuum eigenvalue cheaply, then use it (plus the already-known `β1,β0` and the operator-level `T[i,j]`
expansion) to pin down `α1(S,P),α0(S,P)` without more expensive `(S,P)` sweeps. Once τ2 is fully
closed-form, the natural follow-on is the CG-overlap program this whole `z→1` limit was meant to feed
into (see `HANDOFF.md`'s `su4 L2 untwisted` note in the twisted-case section below for the original
motivation).

## State as of 2026-07-19 (gl(4) L=2 twisted Baxter solver: `Baxter_L2_XXX_SU4.wb` — DONE, in Experiments)

Ported the two-site gl(3) XXX pipeline (`su3_V2_general.wb`) to **gl(4)**: two independent generic
Young-diagram reps `λ={λ1,λ2,λ3,λ4}` per site, four **distinct** diagonal-twist eigenvalues
`z[k]=N[Exp[I Zeta[2k+1]],70]` (k=1..4), 4th-order Baxter/TQ equation, four Q-functions, built fresh
in `Paul/Mathematica/old/XXX/Experiments/Baxter_L2_XXX_SU4.wb` (was empty). Design/plan (with deviation
notes recording what actually happened): `docs/superpowers/specs/2026-07-19-su4-l2-twisted-baxter-design.md`,
`docs/superpowers/plans/2026-07-19-su4-l2-twisted-baxter.md`. All 10 tasks DONE, executed inline
(not via subagent — the wolfbook MCP tools weren't wired into subagent tool sets this session, see
plan/session notes), controller commits.

**Companion twist dropped entirely** (user direction, confirmed before Task 2): no `G`/`TG`/`tG`/
`qminG` — only the diagonal-twist branch (`Tg`/`tg`/`qming`) is built, matching what the gl(3) source
pipeline actually used downstream anyway. `χ1..χ4` (elementary symmetric functions of `z[1..4]`)
still built directly, needed by `qdetT` and the leading-coeff check.

**Genuine gl(4) math (everything else is pure index/term extension from gl(3)):**
- `ν[λ1,λ2][i][u]=(u−θ1−hλ1[[i]])(u−θ2−hλ2[[i]])` for `i=1..4`; `qdetT=χ4·ν[1][u−3h]ν[2][u−2h]ν[3][u−h]ν[4][u]`.
- 4th-order `BaxOp`: `f[u+h]−A1 f[u]+A2 f[u−h]−A3 f[u−2h]+A4 f[u−3h]`.
- **Degree formula carries over unchanged**: `M_k=Λ_1−n_k` (`Λ=λ1+λ2`, `Λ_1=(λ1+λ2)[[1]]`) — verified
  with **zero deviation** against the derived asymptotic `Mdeg` master formula over the full 152-state
  sweep; no candidate-fit re-derivation was needed (the plan's fallback path was never triggered).
- `TauEigensystem4` simultaneously diagonalizes `τ1,τ2,τ3` (Cartan `H1,H2,H3`); `QSolve4` solves for
  four monic `q_k`, `Qfun`=`z[k]^(u/h) q_k(u)`.

**Full-sweep verification (acceptance):** 5 rep pairs, **152 states total**
(`{1,0,0,0}⊗{1,0,0,0}`=16, `{1,0,0,0}⊗{1,1,1,0}`=16, `{1,1,0,0}⊗{1,0,0,0}`=24,
`{1,1,0,0}⊗{1,1,0,0}`=36, `{2,1,1,0}⊗{1,0,0,0}`=60 adjoint⊗fund), every `Q_k` solves the 4th-order TQ
equation and the 5×4 Casoratian reconstruction of `A1..A4` matches `τ1,τ2,τ3,qdetT` —
**worstTQ = worstCasoratian = 0 exactly**. Verified by a fresh-kernel top-to-bottom run (restart → run
cells 1→34 in order); every hard-assert passed in-kernel with no `Abort`.

**Precision bug found and fixed (the real story of this session):** the plan specified `z[k]` at
30-digit precision (matching gl(3)); the initial full sweep at that precision left the d=60
adjoint⊗fund rep's Casoratian residual at `~4.5×10⁻⁸` — 45x over the `10⁻⁹` bar — even though the TQ
residual was *exactly* 0 everywhere (proof the formulas/implementation were correct; this was a pure
numerics issue). Bumping `z[k]` to 50 then 70 digits barely moved the residual (`~1.6×10⁻⁹`, then
bit-identical again) — the real giveaway that something else was capping precision independent of
`z`. Root cause: `TauEigensystem4`'s `rSets` (generic-combination coefficients for splitting the
degenerate τ-eigenspace) and its `Kmix` eigendecomposition were **hardcoded to 30-digit precision** —
a leftover copied verbatim from the gl(3) source's `TauEigensystem3`, never tied to the `z[k]`
precision bump. This silently capped the entire downstream chain (τ-coefficients, `Qfun`, Casoratian)
at ~19-25 accurate digits no matter how much precision `z` carried. Fix: match both to `z`'s 70-digit
precision — eigen-residual improved from a capped `8.6×10⁻²⁶` to `3.46×10⁻⁶⁴`, and the full sweep then
passed with **exact** zero residuals. Lesson for future gl(N) ports: hardcoded precision literals
copied from a smaller-N source notebook are a silent trap — grep for bare `N[..., <number>]` whenever
raising a global precision parameter.

**Next natural steps:** promote to `Clean/` once satisfied; the deliberate follow-on is the `z→1`
untwisted limit (su(4) irrep multiplet degeneracy, feeds the CG-overlap program) — explicitly out of
scope for this notebook per the design spec.

## State as of 2026-07-17 (gl(3) general reps per site: `su3_V2_general.wb` — DONE, in Experiments)

Generalized the two-site gl(3) XXX notebook from symmetric reps `{S,0,0}` per site to **arbitrary
gl(3) Young-diagram reps `λ={λ1,λ2,λ3}`** (`λ1≥λ2≥λ3≥0`, integers), one independent diagram per
site, in `Paul/Mathematica/old/XXX/Experiments/su3_V2_general.wb` (a copy of `su3_V1.wb` generalized in
place). Diagonal-twist frame, `(H1,H2,n)` labelling, three Baxter Q-functions, and the hard-assert
zero-residual convention all preserved. Design/plan (with Execution-outcome section):
`docs/superpowers/specs/2026-07-17-su3-general-reps-design.md`,
`docs/superpowers/plans/2026-07-17-su3-general-reps.md`. Tasks 1–6 all DONE; controller commits.

**Interface:** every function is re-signatured `S1,S2 → λ1,λ2` and takes full 3-vectors per site
(no scalar-`S` wrappers). Numeric params unchanged (`θ1=1/3, θ2=1/7, h=1`, high-precision `z[k]`).

**The three genuine math edits (everything else is pure relabeling):**
- **General Yangian weight** `ν[λ1,λ2][i][u] = (u−θ1−h·λ1[[i]])(u−θ2−h·λ2[[i]])` for `i=1,2,3`,
  from the Lax action `L[λ][i,i][u]|hw⟩=(u−h λ_i)|hw⟩` on the two-site HWS. Verified via the HWS
  eigenvalue check (T[i,i]|hw⟩=ν[i]|hw⟩) and the `qdetT == tg[·][3,1]` shift-ordering guard on
  non-symmetric reps. `qdetT[λ1,λ2][u]=χ3 ν[·][1][u−2h] ν[·][2][u−h] ν[·][3][u]` (V1 shift order
  survives).
- **Honest tensor-product `wtMult`:** `siteWeightMult[λ]=Counts` of the gl(3) weight of every GT
  pattern of `λ`; `wtMult[λ1,λ2][H1,H2]` recovers `(n1,n2,n3)` from `(H1,H2)` with `N=Total[λ1]+
  Total[λ2]` and convolves the two single-site tallies — correctly yields interior multiplicity `>1`.
  Sum rule `Σ wtMult == dim[λ1] dim[λ2]` hard-asserted across the sweep; guard on `{2,1,0}⊗{2,1,0}`
  (d=64) passes.
- **Degree closed form (Task 5 headline):** `M_k = Λ_1 − n_k`, with `Λ = λ1+λ2` the pseudovacuum
  highest weight and `Λ_1 = (λ1+λ2)[[1]]`. Sum rule `Σ_k M_k = 3Λ_1 − N`. Reduces to V1's `M_k = N−n_k`
  for symmetric reps (there `Λ_1 = N`). Extracted by fitting candidates `{N−n, Λ_1−n, Λ_k−n}` against
  the rep-agnostic `Mdeg` over the full sweep — only `Λ_1−n` gives zero deviation on the non-symmetric
  reps. Hard-asserted (`Mpred==Mdeg`, worst=0).

**Full-sweep verification (acceptance, Task 6):** six rep pairs, **102 states total**, every `Q_k`
solves the 3rd-order TQ equation and the 4×4 Casoratian reconstruction of `A1,A2,A3` matches
`τ1,τ2,qdetT` — **worstTQ = worstCasoratian = 0** (machine-zero, both `< 10^-9`). Per-rep counts:

| `λ1` | `λ2` | nstates |
|---|---|---|
| `{1,0,0}` | `{1,0,0}` | 9 |
| `{2,0,0}` | `{1,1,0}` | 18 |
| `{1,1,0}` | `{1,1,0}` | 9 |
| `{2,1,0}` | `{1,0,0}` | 24 |
| `{2,1,0}` | `{1,1,0}` | 24 |
| `{2,2,0}` | `{1,0,0}` | 18 |

Verified by a fresh-kernel top-to-bottom run (restart → run cells 1→40 in order); every hard-assert
check passed in-kernel with no `Abort`.

**wolfbook output-persistence glitch — re-run, do not trust blank boxes.** Several check cells created
via `insertCells` / edited during Tasks 3–5 do NOT refresh their on-disk `Out[]` display even after a
correct fresh-kernel re-run (confirmed: re-running an affected cell individually still shows blank).
These cells computed the correct result in-kernel (verified by direct `evaluateExpression`), but their
saved display is empty or stale. **A future reader must re-run the notebook rather than treat these as
failures.** Cells still showing empty/stale output on disk after the full Task-6 run (1-based cell
numbers; correct in-kernel values in parentheses):
- Cell 14 — `qdetT == tg[·][3,1]` guard: **empty** (in-kernel `{0}`).
- Cell 23 — `wtMult` sum-rule: **empty** (in-kernel all reps OK; the informational
  `adjoint⊗fund @(0,0)` spot value is `0` because for N=4 the (0,0) weight is non-integer — not a
  failure).
- Cell 24 — `{2,1,0}⊗{2,1,0}` d=64 guard: **empty** (in-kernel `tot=64=dim`).
- Cell 27 — eigensystem residual spot check on `{2,1,0}⊗{1,0,0}`: **empty** (in-kernel worst≈1.5×10⁻²⁷).
- Cell 31 — `degData` tabulation: **stale** (shows an old `{H1H2n→…,sumM→…}` format; source now emits
  raw `{λ1,λ2,H1,H2,n,{M}}` rows).
- Cell 32 — degree candidate-fit `results` Association: **stale** (shows an old matrix; source now
  emits `{"N-n"→…,"Lam1-n"→0,"Lamk-n"→…}`).
- Cell 40 — full-sweep report: **stale** (shows the old symmetric sweep, 141 states; the current
  λ-vector run gives the 102-state table above with worstTQ=worstCasoratian=0).

Check cells whose on-disk `Out[]` IS correct: cell 4 (`C1`→{0}), cell 10 (HWS→{0}), cell 21
(commuting-family, maxCommutator≈10⁻²⁶), cell 35 (degree assert worst→0; only its result-label string
still reads "N-n_k", cosmetic), cell 39 (single-state TQ perQresidual→{0,0,0}).

**Next natural steps:** promote `su3_V2_general.wb` to `Clean/` once satisfied; then the companion-twist
SoV / FSoV CG-overlap program for general gl(3) reps (mirroring the su(2) Sections C–F), and the
`paper/xxxCG.tex` write-up.

## State as of 2026-07-17 (gl(3) Section B: eigensystem + Baxter Q-functions — DONE, in Experiments)

Built the gl(3) analog of the su(2) Clean notebook's Section B, in
`Paul/Mathematica/old/XXX/Experiments/su3_V1.wb` (on top of the pre-existing Section A there). All
committed to `main`, all hard-asserts zero-residual. Design/plan (with Execution-outcome sections):
`docs/superpowers/specs/2026-07-17-su3-transfer-matrix-eigensystem-baxter-design.md`,
`docs/superpowers/plans/2026-07-17-su3-transfer-matrix-eigensystem-baxter.md`.

**What's there (Section B cells, diagonal-twist frame, numeric `θ1=1/3,θ2=1/7,h=1`,
`z[k]=N[Exp[I Zeta[3/5/7]],30]`):**
- `tgCoeffs` (coeff matrices of `τ1=tg[·][1,1]`, `τ2=tg[·][2,1]`), Cartan `H1op`/`H2op`, commuting-
  family check (`maxCommutator=0`).
- `wtMult[S1,S2][H1,H2]` — weight-multiplicity predictor (sum rule verified).
- `TauEigensystem3[S1,S2]` — simultaneous diagonalization; every state labeled `(H1,H2,n)`
  (`H1=E11−E22, H2=E22−E33`; `n` sorted by the full `τ1/τ2` charge tuple), sector counts match
  `wtMult`. Accessors `Psi`/`tau1`/`tau2`.
- `BaxOp` (TQ residual functional), `alphaCoeffs`, `Mdeg`, `qpoly`, `QSolve3` (three-Q solver),
  `Qfun` accessor `Q_k = z_k^{u/h} q_k(u)`.

**Key result — the degree relation:** `deg q_k = M_k = N − n_k` (`N=S1+S2`, `(n1,n2,n3)` the gl(3)
weight components), equivalently `M1=(2N−2H1−H2)/3, M2=(2N+H1−H2)/3, M3=(2N+H1+2H2)/3`. `M_k` pairs
with `z_k` in-index, depends only on `(H1,H2)`, `Σ M_k = 2N`. Hard-asserted.

**Verification:** full sweep `{(1,1),(2,1),(1,2),(2,2),(3,1),(1,3)}`, 141 states — every `Q_k`
solves the 3rd-order TQ equation (`worstTQ=0`) and the 4×4 Casoratian reconstruction of `A1,A2,A3`
matches `τ1,τ2,qdetT` (`worstCasoratian=0`). Commits `b4b68c7`…`645517d`.

**Next natural steps:** promote `su3_V1.wb` Section B to `Clean/` once satisfied; build the
companion-twist SoV covectors / FSoV scalar product for gl(3) (the CG-overlap program, mirroring
the su(2) Sections C–F); write up in `paper/xxxCG.tex`.

## State as of 2026-07-16 (SoV factorization + functional CG overlap — DONE, promoted to Clean)

Big result this session: the **Clebsch–Gordan overlap correspondence closed exactly**, then
promoted from `Experiments/` into `Clean/`. Everything below is **committed**. Design/plan docs
(with `Result`/`Execution outcome` sections recording what actually happened):
`docs/superpowers/specs/2026-07-16-xxx-sov-factorization-functional-overlap-design.md`,
`docs/superpowers/plans/2026-07-16-xxx-sov-factorization-functional-overlap.md`.

### The result

**The 1×1 FSoV determinant of the twisted Q-function equals the operatorial overlap of a site-1
unipotent rotation with the two-site eigenvector:**

    det1[λ1][ k^(u/h) Q1[λ1,λ2][M,n] ] = (⟨x0[λ1]| exp[-k E21] ⊗ ⟨x0[λ2]|) . Psi[λ1,λ2][M,n]

where `det1[λ1][f] := br1[λ1][f]/(f/.u->θ1)` (the 1-site FSoV determinant — Section C's `det`
restricted to one site). Dictionary: **generator E21, φ(k) = -k** (α=-1, universal). Residual is
the literal symbolic `0` across the rep sweep {(1,1),(2,1),(1,2),(2,2),(3,1),(1,3)} — 41 states.

This rests on a prerequisite factorization: the s2=0 slice of the L=2 SoV basis,
`x2[λ1,λ2][{s1,0}] = xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}]` **exactly** (no leftover scalar), where
`xL1` is a genuine standalone length-1 SoV covector (new `t1`/`qdet1` transfer-matrix hierarchy).

**Normalization correction (same day, post-Clean-promotion):** the first pass of this factorization
check found an apparently nontrivial site-2 "passthrough eigenvalue" `μSite2[λ2,s1] =
Product[(θ1-θ2+(k-1)h-λ2 h),{k,1,s1}]`. **This was wrong** — a bug in `x2`'s own normalization, not
real physics, found by the user. `x2`'s denominator used `Qθ[u]=(u-θ1)(u-θ2)` evaluated at a point
shifted by `-λ1h`, which incorrectly puts that same `-λ1h` shift into the θ2 factor too (which
should carry `-λ2h`, site 2's own weight) — that mismatch produced exactly the spurious `μSite2`.
**Fix:** a new `ν1[λ1,λ2][u] := (u-θ1-hλ1)(u-θ2-hλ2)` (the Yangian weight eigenvalue of T11 on the
HWS; deduped with Section B's `aFun`, same formula) replaces `Qθ` in `x2`'s denominator, argument no
longer carrying the `-λh` shift separately. With the fix, the factorization ratio is **exactly 1**
for all 27 `(λ1,λ2,s1)` cases — proved analytically too: `ν1[λ1,λ2][θ1+(k-1)h]` factors as
`((k-1-λ1)h)·(θ1-θ2+(k-1-λ2)h)`, and that second factor is exactly the old spurious `μSite2`. `xL1`
needed no functional change (`Qθ1` is linear, fix is a numeric no-op) but was rewritten with a
single-site `ν1[λ][u]:=(u-θ1-hλ)` for parity. **Section F confirmed unaffected** (`dictReportDet`
identical before/after). Applied to both `Clean/XXX_CG_L2_Clean.wb` and
`Experiments/XXX_CG_L2_V2.wb`; fresh-kernel Run-All clean on both after the fix (also caught and
fixed one unrelated hidden-dependency bug in V2's `obstructionReport` cell, same class as the
`sts12` issue during Clean promotion — a scratch-only `ABTab` variable that was never a saved cell).

**The crux / what nearly went wrong:** the first functional-discovery pass used the bare bracket
`br1[k^(u/h) Q1]` (no determinant normalization) and did NOT close (k-dependent per-state scalar).
The fix was the `1/f(θ1) = 1/(k^(θ1/h) Q1(θ1))` normalization built into the FSoV determinant — it
cancels the `k^(θ1/h)` prefactor and the state-dependent `Q1(θ1)`, leaving `det1[k^(u/h) Q1]` a clean
degree-λ1 polynomial in k starting at 1. Caught by the user asking "was it normalized by
k^θ1/h Q1(θ1)?". If a functional-SoV object won't universalize, check for a dropped determinant
node-normalization first (see memory `reference_fsov_determinant_normalization`).

**Gotcha found this session:** `MatrixPower[nilpotentMatrix, 0]` fires `MatrixPower::sing` and stays
unevaluated even for a concrete integer power 0. Build `exp[φ E]` for nilpotent E as
`IdentityMatrix[dim] + Sum[φ^m/m! MatrixPower[E,m],{m,1,λ}]` (split off the identity), not `MatrixExp`
and not a Sum starting at m=0.

### Commits (chronological)

`737f6a2` spec → `1bd9bcd` plan → `8dbed52` Section E in `Experiments/XXX_CG_L2_V2.wb` → `9db89eb`
Section F in `Experiments/` + corrected docs → `dcf5c1e` **promoted Sections C,D,E,F into
`Clean/XXX_CG_L2_Clean.wb`** (now a contiguous 71-cell A→F notebook, fresh-kernel top-to-bottom run
verified, every hard-assert exact) → `bfc0ffb` refreshed `Clean/XXX_CG_L2_Clean.md`/`.tex` reference
docs to cover C–F (`.tex` compiles clean with plain `pdflatex`, `\mathbbm{1}` convention).

Section F was promoted in **final form only**: the corrected `det1`/`FLdet` machinery and the closed
dictionary. The exploratory bare-bracket `FL`, the pre-`det1` per-state-scalar cells, and the
`obstructionReport` were dropped from Clean (they remain in `Experiments/XXX_CG_L2_V2.wb` as the
derivation record, same convention as other Clean promotions).

## Open threads / next steps

- **Next up (2026-09-09):** build the perturbative-in-twist solver in
  `Experiments/SU4_twisted_Baxter_L2.wb` — keep `Q_j = Exp[phi_j u] x` (polynomial of fixed degree
  `M_j`) non-perturbatively, expand only the polynomial coefficients, as Laurent series with the
  now-known leading orders `eps^(-(M_k-M_1))`. Everything it needs (solver, degrees, `N_I`,
  divergence orders, exact-in-`eps` test states) is in place.
- Then: wire it to the FSoV norm (`SU4_norm_untwisting_L2_general.wb`), and scan over reps
  (item 3 of the 2026-09-02 plan, still untouched).
- Decide where the new gl(4) results get written up — `paper/xxxCG.tex` is gl(2)/single-site only,
  so a gl(4) L=2 section is a scope call; alternatively promote to `Clean/` and run
  `wolfbook-docs-agent`. Also still unwritten: L=2 factorization, scalar product, CG-overlap
  dictionary.
- Older, still open: extract the actual CG coefficients from the overlap generating function (the
  `exp[-k E21]` rotation generates the CG combinations — expanding in k gives per-state
  coefficients); push the factorization to the s2≠0 case / general L.
- `Experiments/XXX_CG_L2_V2.wb` still holds the full discovery narrative (bare-bracket attempt,
  obstruction report, per-state-scalar exploration) — useful if revisiting why the determinant
  normalization matters.
- Older restructuring/Clean-notebook state (L3 build, docs agent, shared wolfbook kernel gotcha) is
  unchanged from the 2026-07-15 session — see `git log` and `CLAUDE.md`.

## Stable structure

See `CLAUDE.md` for the full layout (XXX vs Gaudin split, Clean/Experiments convention, subagents).
`Clean/XXX_CG_L2_Clean.wb` is now the authoritative A–F notebook: A (SoV construction, incl. `ν1`) ·
B (Baxter/TQ, Psi/Q1) · C (2-site FSoV scalar product/orthogonality) · D (left eigenvectors,
`<Psi|Psi>=det[Q,Q]`) · E (exact L=1 factorization, c=1) · F (functional CG-overlap dictionary).
Companion `.md`/`.tex` reference docs are up to date with it (refreshed same day for the Section E
correction; `.tex` compiles clean with plain `pdflatex`).
