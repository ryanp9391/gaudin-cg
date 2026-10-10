# 7. Undoing the site-2 rotation of the two-rotation family: does the site-1 covector span every weight space?

**Date:** 2026-10-09 (evening). **Question (Paul).** The two-rotation family of `Two-rotations-det-overlap.md`
§5.6b/§5.8 gives product covectors ⟨v₁ ⊗ v₂| with a Q-determinant overlap, where v₂ = ⟨0|R is a rotation of
the vacuum. Undo the rotation on site 1, ⟨v₁|R⁻¹, project onto the weight spaces of the site-1 rep V, and
expand in all free parameters. Do these projections span every weight space of V (so that all highest-weight
CG coefficients ⟨μ ⊗ top | ν, HW⟩ are reachable)? The older families (2026-09-23,
`old-conventions/XXX-CG/xxx-210-cg-test.tex`) reached only the site-1 weights with μ₃ ≤ λ₂ and missed the singlet.

**Status:** spec + analytic prediction (Claude, this file). Implementation: `L2-rotation-sov.wb` §12 (adjoint),
then `L2-310-two-rotations.wb` ([3,1,0]).

---

## 1. Why the undo is exact (no untwisting needed for the count)

Diagonal frame (twist diag(z), GT basis, `Ee[i,j] = JJ[i,j,Λ]`). e₁ᵀ is the highest-weight covector of V
(weight Λ = (2,1,0), `J0[k][ps[{Λ}][[1]]]`). Every rotation of the vacuum in the family is of the form
v₂ = s · e₁ᵀ ρ(n) with n ∈ N₊ (upper unipotent 3×3), because `oneSiteVac[z,w]` is literally
e₁ᵀ·MatrixExp[r₁E₁₂]·MatrixExp[r₂E₁₃]·MatrixExp[r₃E₂₃]. For n ∈ N₊, ρ(n)⊗ρ(n) fixes every highest-weight
ket of V⊗V, hence

  ⟨v₁ ⊗ v₂ | ν,HW⟩ = s · ⟨v₁ ρ(n)⁻¹ ⊗ e₁ᵀ | ν,HW⟩ = s · Σ_{m ∈ V_μ} φ_m CG(ν,HW | m, top),  φ := v₁ρ(n)⁻¹, μ = ν − Λ.

So the HW CG coefficients at site-1 weight μ are reachable iff the weight-μ projections of the φ's span V_μ.
This is a finite-twist, exact statement (the old thread's Result "res:nplus"). The untwisting limit only enters
through the eigenstates, not through the covectors; its per-weight leading orders are a secondary check (§4.7).

## 2. The family in the diagonal frame (adjoint)

From `L2-rotation-sov.wb`: `oneSiteVac[zz, ww]` (cell 98), `ρdiag[d]` and `Dy[y1] = 1/(y1 − z_i)` (cell 10),
`w2Of[y] = y(θ₁₂−2h)/(θ₁₂−h)` (cell 109), `x01` (cell 7), `Xdiag` (cell 88). Cell 111 builds the partner as
`g . ρdiag[Dy[w2G]]`, i.e. in the diagonal frame the resolvent ρ((w₂ − Gᵀ)⁻¹) is the diagonal torus element
Δ = ρdiag[1/(w₂ − z_k)].

Family (fixed w₁, free y; w₃ drops out):

  v₂(y) = Ω_{(w₁,y)} = e₁ᵀ ρ(n(y)),  n(y) = e^{r₁E₁₂} e^{r₂E₁₃} e^{r₃E₂₃} = [[1, r₁, r₂ + r₁r₃],[0,1,r₃],[0,0,1]],
  v₁(y) = v₂(y) ρ(Δ(y)),  Δ(y) = diag(δ₁,δ₂,δ₃), δ_k = 1/(ω − z_k), ω := w2Of[y],

with (cell 6 / 98, second vacuum parameter = y)

  r₁ = (w₁−z₁)/(w₁−z₂),  r₂ = (w₁−z₁)(z₁−y)(z₂−z₃)/((z₁−z₂)(w₁−z₃)(y−z₃)),  r₃ = (w₁−z₂)(z₂−y)(z₁−z₃)/((z₁−z₂)(w₁−z₃)(z₃−y)).

Undone covector:

  φ(y) = v₁ ρ(n)⁻¹ = e₁ᵀ ρ(n Δ n⁻¹) = e₁ᵀ ρ(Δ) ρ(n′) = δ^{Λ} · e₁ᵀ ρ(n′),  n′ := Δ⁻¹ n Δ n⁻¹ ∈ N₊,

  n′ = [[1, a′, b′],[0, 1, c′],[0, 0, 1]] = e^{c′E₂₃} e^{b′E₁₃} e^{a′E₁₂}   (matrix entries as exponents, in this order),
  a′ = r₁(δ₂−δ₁)/δ₁,  c′ = r₃(δ₃−δ₂)/δ₂,  b′ = r₁r₃(δ₁−δ₂)/δ₁ − (r₂+r₁r₃)(δ₁−δ₃)/δ₁.

The torus (generic g = Ω_{(w₁,y)}ρ(D) on site 2, §5.8) only rescales the weight components of φ by d^μ/χ(D),
so the per-weight spans of the whole 3-dim family V equal those of the 1-parameter family φ(y) at any fixed w₁.
Varying w₁ must therefore not change any rank (consistency check).

## 3. Analytic prediction (to be verified exactly)

(i) **b′ = 0 ⟺ ω = y.** Using the r's: r₂/(r₁r₃) = (z₃−z₂)(y−z₁)/((z₁−z₃)(y−z₂)), while b′ = 0 ⟺
r₂/(r₁r₃) = (z₃−z₂)(ω−z₁)/((z₁−z₃)(ω−z₂)). Since ω = w2Of[y] = y(θ₁₂−2h)/(θ₁₂−h) ≠ y, b′ ≠ 0 on the family.
(The old blind spot was exactly b′ = 0: `G⁻¹g₀Gg₀⁻¹` has zero (1,3) entry by the p₃ relation; same algebra
with Δ = diag(z).) The "paper-like" assignment (resolvent at the vacuum's own y) would have b′ ≡ 0.

(ii) **b′/(a′c′) = (y − ω)/(y − z₂)** — non-constant in y, w₁-independent.

(iii) Weight components of e₁ᵀρ(n′) for the adjoint (e₁ᵀE₁₂² = e₁ᵀE₂₃² = 0; u's are fixed nonzero covectors):
(2,1,0): 1 · e₁ᵀ; (1,2,0): a′; (2,0,1): c′; (1,1,1) [dim 2]: b′ e₁ᵀE₁₃ + a′c′ e₁ᵀE₂₃E₁₂ → rank 2 by (ii);
(0,2,1): a′²c′ u + a′b′ u′; (1,0,2): b′c′ u″; (0,1,2) [singlet direction]: −b′(b′ + a′c′) · (normalization).
**Prediction: all 8 directions of V are spanned, all 6 HW states of adj⊗adj are reached (rank 6 on the HW span).**

Untwisting z_k = 1 + εφ_k: a′, c′ = O(ε); b′, a′c′ = O(ε²) with ratio → (y−ω)/(y−1) (still non-constant);
(0,1,2) component O(ε⁴). Expected per-weight leading orders ε⁰, ε¹, ε¹, ε² (rank 2), ε³, ε³, ε⁴.

## 4. Implementation — adjoint, `L2-rotation-sov.wb`, new §12 (insert after cell 112)

Rebuild the kernel state first if needed (cells 2–21, then 24–27, 60–68, 86–106; cell 2 does `ClearAll` — make
sure no other notebook shares the kernel). Everything below is exact (rationals) at P1 unless stated. Use
`checks[...]` entries and a final `Grid` result cell as in the rest of the notebook. Write `y-star` (never
`y*)`) inside comments.

4.1 **Setup.** `nMat[y_]` (3×3 N₊ matrix from the r's with w = {w1, y, w3}), `ρN[y_]` = the d1×d1 ρ(n(y)) built as
`MatrixExp[r1 Ee[1,2]].MatrixExp[r2 Ee[1,3]].MatrixExp[r3 Ee[2,3]]`, `v2[y_] = oneSiteVac[zN,{w1,y,wN[[3]]}]`,
`v1[y_] = v2[y].ρdiag[Dy[w2Of[y]]]`, `φ[y_] = v1[y].Inverse[ρN[y]]`. Check v2[y] == e₁ᵀ.ρN[y] exactly, and that
v₁⊗v₂ at y = 7/5 (w₁ = −7/3) reproduces the cell-109 state (same objects). Avoid y ∈ {z_k} and ω ∈ {z_k}.

4.2 **HW kets and exactness of the undo.** Raising generators on kets: `E1[i,j] + E2[i,j]` with i<j (check that
e₁ (ket) is annihilated by Ee[i,j], i<j, in this convention; if instead i>j annihilates e₁, swap and say so).
HW kets of V⊗V = common null space of the total raising generators within each total-weight sector
(6 vectors: (4,2,0),(3,3,0),(4,1,1),(3,2,1)×2,(2,2,2)). Check (ρN[y]⊗ρN[y]).HW == HW exactly for two y's, and
⟨v₁⊗v₂|HW⟩ == s⟨φ⊗e₁ᵀ|HW⟩ with s = (v₂.e₁-component).

4.3 **Weights.** Weight of GT pattern p: `Table[J0[k][p],{k,3}]`. Group the 8 basis indices by weight.

4.4 **Per-weight ranks, numeric sample.** 12 rational y's (e.g. y ∈ {7/5, −2/3, 3, −1/2, 11/9, 4/7, −9/5, 13/6,
−3/11, 8/3, 1/6, −7/4}); for each weight μ, rank of the matrix of φ(y_i)|_μ. Report rank vs dim V_μ (dims
1,1,1,2,1,1,1 → expect all equal).

4.5 **Per-weight ranks, exact in y.** With y symbolic: `Together[φ[y]]`, numerator polynomial coefficient
vectors per weight → rank. Also list which components vanish identically (expect none).

4.6 **Closed forms.** Compute (a′,b′,c′) from n′ = Δ⁻¹nΔn⁻¹ symbolically (y symbolic, rational z, w₁, θ, h) and
verify `Simplify[b′/(a′c′) − (y − ω)/(y − zN[[2]])] == 0`, and that `b′ /. ω -> y` is 0 (substitute the resolvent
point by y, i.e. use Dy[y] instead of Dy[w2Of[y]]: the (0,1,2) component must then vanish identically — the old
blind spot reproduced). Verify the (0,1,2) component of φ(y) ∝ b′(b′ + a′c′) (ratio y-independent).

4.7 **Untwisting leading orders (secondary).** z → {1+ε φ₁, 1+ε φ₂, 1+ε φ₃} with rational φ = {3/7, −5/4, 2/9}
(or symbolic if cheap): everything is rational in z, so substitute into the r's, Dy, ρN, and Series each weight
component of φ(y) in ε to order 5 with y symbolic; report the leading order per weight and the rank of the
leading coefficient vectors over y (expect ε⁰,ε¹,ε¹,ε²(rank 2),ε³,ε³,ε⁴). Do this with the 3×3 (a′,b′,c′) and the
explicit weight formulas of §3(iii) if the 8×8 Series is slow; time a single component first.

4.8 **HW count.** M[i,ν] = ⟨φ(y_i)⊗e₁ᵀ|HW_ν⟩ for the 12 y's and 6 HW kets: rank (expect 6). Per irrep: rank of the
columns of that irrep (octets: 2). Compare with the old families' 5/6.

4.9 **Torus consistency.** Repeat 4.4 with w₁ ∈ {9/4, 5/2} and with the union over the three w₁'s: ranks unchanged.

4.10 **Result cell** (markdown + Grid): the family, the undo, the ranks table, the closed forms, the one-line
mechanism (b′ ≠ 0 because w₂(y) ≠ y), and the comparison with 2026-09-23. Save the notebook.

## 5. Implementation — [3,1,0], `L2-310-two-rotations.wb` (second pass, after the adjoint)

Point B family (handoff 2026-10-09 evening): v₂ = Ω_{(w₁,y*)} with y* = w₂(θ₁₂−h)/(θ₁₂−3h); v₁ =
e₁ᵀ(α + βE₁₂)·R(z,(w₁,y*))·ρ(Δ_{w₂}) with β/α fitted per point (392/2335 at P1). Undo with n = R(z,(w₁,y*))
(the Gauss element of v₂; check it is the same r_k(z,w) formula as the adjoint). Then φ = e₁ᵀ(α+βE₁₂)ρ(nΔn⁻¹) =
(α e₁ᵀ + β δ₂/δ₁ e₁ᵀE₁₂)·δ^{Λ}·ρ(n′). V = [3,1,0] has 15 states, 10 distinct weights, (2,1,1),(1,2,1),(1,1,2) of
dim 2. Prediction: the rotation part alone spans everything by the same (y−ω)/(y−z₂) argument (e₁ᵀE₁₂³ = 0,
e₁ᵀE₂₃² = 0 but e₁ᵀE₁₃³ ≠ 0); the dressing only adds. To get a y-family one must re-run the point-B pipeline at
several w₂ (β/α has no closed form): do 4–5 values of w₂, then the per-weight ranks (need ≤ 2) and the HW count
(10 HW states of [3,1,0]², expect 10; the old families got 8/10, missing [4,2,2],[3,3,2]). Also do point A
(adjoint's c* on the column-shifted shapes) if cheap. Same kernel caveat: this notebook's `ClearAll` wiped the
adjoint notebook's state on 2026-10-09.

---

## 6. Adjoint outcome (2026-10-09, `L2-rotation-sov.wb` §12 = cells 113–123, sonnet wolfbook-builder; saved, not committed)

- **Per-weight spans: full, 8/8** — ranks (1,1,1,2,1,1,1) = dims at w₁ = −7/3 (12 numeric y's and exact in y), unchanged for
  w₁ = 9/4, 5/2 and for the union. No component vanishes identically.
- **Control (resolvent point set to y, i.e. ω = y):** (0,1,2) and (1,0,2) vanish identically and (1,1,1) has rank 1 —
  exactly the old μ₃ > λ₂ blind spot.
- **Closed forms exact:** a′, b′, c′ as in §2; b′/(a′c′) = (y−ω)/(y−z₂); b′ ≡ 0 iff ω = y. The (0,1,2) component is
  3·b′(a′c′ − b′)·δ^Λ (sign differs from §3(iii); the conclusion is unchanged).
- **HW count:** rank 5 at fixed w₁ (one cross-weight linear relation among the functions a′, c′, b′, a′c′ of y — harmless
  for reading CG coefficients state by state, since each state has a definite ν) and **rank 6 with the torus**
  (12 random (y, D) pairs; torus acts as the character Π D_i^{ν_i} on each HW overlap, verified exactly). Per irrep:
  (4,2,0), (4,1,1), (3,3,0) rank 1, octet pair rank 2, singlet nonzero. Control ω = y: rank 4 even with the torus
  (singlet column 0, one octet direction lost).
- **Untwisting leading orders:** ε⁰, ε¹, ε¹, ε² (rank 2), ε³, ε³, ε⁴ at (2,1,0), (1,2,0), (2,0,1), (1,1,1), (0,2,1),
  (1,0,2), (0,1,2) — as predicted.
- Kernel note: only cells 2–21 were rebuilt (plus the text of `oneSiteVac`, `w2Of`); cells 63–111 are not live.

## 7. [3,1,0] implementation details (for the second pass, `L2-310-two-rotations.wb`)

Live helpers there: `rOm[zv,wv]`, `Rrot[zv,wv]` (cell 4), `Om1[zv,wv]`, `ρdiag`, `Dy`, `Ee`, `E1/E2`-type two-site
generators, `rotRank` (cell 4), `pipeline310[zz, ww, hh, th]` (cell 12: redoes the exact product locus for explicit
(h,θ,z,w) and returns, per nontrivial component, the site-2 torus invariant y, rotRanks, c-vector and the
decomposition of site 1 in span{e₁RD, e₁E₁₂RD} with R = Rrot(z,(w₁,y)), D = ρdiag[1/(w₂−z)], i.e. β/α).
Plan: run `pipeline310[zN, {5/2, w2, 3/5}, hN, thN]` for w₂ ∈ {7/3, −8/9, 4/27, 3, −1/2, 11/7} (point B is the
component with y = w₂(θ₁₂−h)/(θ₁₂−3h)); for each, build v₂ = Om1[zN,{w₁,y*,w₃}] = e₁ᵀρ(n), n = Rrot(zN,{w₁,y*,w₃})
(check it is the N₊ element with the same r_k as the adjoint), v₁ = (α e₁ᵀ + β e₁ᵀEe[1,2]).ρ(n).ρdiag[Dy[w₂]], and
φ = v₁.Inverse[ρ(n)]. Then: exactness of the undo on the 10 HW kets of [3,1,0]⊗[3,1,0]; per-weight ranks over the
six w₂'s (10 distinct weights; (2,1,1),(1,2,1),(1,1,2) have dim 2); the same with β = 0 (rotation part only) and
with ω → y* (control); HW count with and without the torus (expect 10); which irreps the old families missed:
[4,2,2], [3,3,2] (site-1 weights (1,1,2), (0,2,2)).

## 8. [3,1,0] outcome (2026-10-09, `L2-310-two-rotations.wb` cells 24–33, sonnet wolfbook-builder; saved, not committed)

Point-B family at w₂ ∈ {7/3, −8/9, 4/27, 3, −1/2, 11/7} (w₁ = 5/2; β/α = 392/2335, 64/731, −32/2409, 8/77, 1/23,
4312/2693, from `pipeline310`); v₁ = (e₁ᵀ + β e₁ᵀE₁₂)ρ(n)ρdiag[Dy[w₂]] reproduces the pipeline's site-1 covector
exactly; undo exact on all 10 HW kets ((6,2,0),(6,1,1),(5,3,0),(5,2,1)×2,(4,4,0),(4,3,1)×2,(4,2,2),(3,3,2)).
[3,1,0] has 12 distinct weights (dims 1 except (2,1,1),(1,2,1),(1,1,2) of dim 2).

| variant | per-weight ranks | HW rank (torus) |
|---|---|---|
| (a) point B as found (dressed) | all full **except** (0,1,3), (1,0,3) ≡ 0 and (1,1,2) rank 1 → 12 of 15 directions | **10/10** |
| (b) rotation part only (β = 0) | all full, 15/15 (the §3 prediction) | 10/10 |
| (c) dressed, resolvent point → y* | as (a) | 10 |
| (c0) β = 0, resolvent point → y* | μ₃ ≥ 2 components vanish (old blind spot) | 6 ([4,2,2],[3,3,2] lost) |

- So the dressing e₁ᵀ(α+βE₁₂) does **not** only add: with the fitted β/α(w₂) it cancels exactly the μ₃ = 3 components
  and one direction of (1,1,2), for every w₂ (also 12 random w₂). Consistent with site 1 lying in the e₁ᵀE₁₂ orbit
  (rotRank 5). Spec §5's "the dressing only adds" was wrong.
- Consequence for HW CG coefficients: all 10 HW states have nonzero overlaps (old families: 8/10, [4,2,2] and
  [3,3,2] missed); 9 of them are fully determined; for [4,2,2] (site-1 weight (1,1,2), dim 2) only one linear
  combination of the two CG components is determined by this family. The missing μ₃ = 3 weights (0,1,3), (1,0,3)
  carry no HW state and matter only for a complete site-1 basis.
- 3×3 identity b′/(a′c′) = (y*−w₂)/(y*−z₂) holds fully symbolically (also symbolic z).
- Fixed-torus HW rank over 12 w₂'s is 7 (cross-irrep relations only); torus-completed 10.
- Kernel: ran on K1; the adjoint notebook's live state is wiped again (file intact, 123 cells).

## 9. [3,1,0]: lowest-weight frame and point A (2026-10-09, `L2-310-two-rotations.wb` cells 34–42, sonnet; saved, not committed)

Paul's follow-up: with point B alone the [4,2,2] HW vector is fixed only up to one ratio (rank 1 of 2 at site-1
weight (1,1,2)); try (1) the lowest-weight frame of the same family, (2) point A.

1. **LW frame.** e_lowᵀ (GT index 15, weight (0,1,3)) is fixed by B₊; v₂ = s·e_lowᵀρ(m), m ∈ N₋ from the exact
   factorization w₀⁻¹n = b₊m; ρ(m)⊗ρ(m) fixes all 10 LW kets; ⟨f1⊗v₂|LW⟩ = s⟨f1ρ(m)⁻¹⊗e_lowᵀ|LW⟩ exact for all six
   w₂. Per-weight ranks of the LW-undone covector, order (0,1,3),(0,2,2),(0,3,1),(1,0,3),(1,1,2),(1,2,1),(1,3,0),
   (2,0,2),(2,1,1),(2,2,0),(3,0,1),(3,1,0): **1,1,1,1,2,2,1,1,1,1,0,0** — the exact mirror of the HW frame
   (0,1,1,0,1,2,1,1,2,1,1,1). [4,2,2]'s LW slice sits at (2,1,1): rank 1 of 2. **Not determined in either frame**;
   the two frames measure different vectors (HW vs LW of the irrep) with unrelated normalizations, so they do not
   combine. Torus-completed LW count 10/10 (all reached). All other nine HW states determined in both frames.
2. **Point A** (site-2 invariant y_A = w₂(θ₁₂−h)/(θ₁₂−2h), site 1 rotRank 6, found at every w₂): v₂_A ∥ e₁ᵀρ(n_A)
   with the plain r_k(z,(w₁,y_A)); undo exact. **Per-weight ranks all full (15/15)**, nothing vanishes; A alone and
   A ∪ B: every HW state fully determined, torus-completed HW count 10. So **point A closes the [4,2,2] gap.**
   Caveat: point A's overlap is c_A·(1, τ_λ(θ₁−3h))·⟨Ω|Ψ_B⟩ — a transfer-matrix polynomial times the vacuum det₆,
   hence a Q-number through the T–Q relations; **superseded by §13: A IS a ladder det₂, over the pole set {0,1,2}.**

Verdict for [3,1,0]²: all 10 HW states reached by B; 9 fully determined by B; the 10th ([4,2,2]) needs A.

## 10. [3,2,0] (2026-10-09, new notebook `L2-320-two-rotations.wb`, 21 cells, sonnet; saved, not committed)

- **Product locus at θ₁−3h** (localizers T₁₁, T₂₁, T₃₁; T₂₂, T₃₂, T₃₃ vanish on ⟨Ω|; joint site-2 span {Ω₁, ω₂} with ω₂
  of rotRank 5): the localizing span plus **one** isolated point with c = (81/1075, 36/215, 4/39, 1, 860/1053,
  −860/729, −92450/85293) on (1,T₁₀,T₁₁,T₂₀,T₂₁,T₃₀,T₃₁) — the [3,1,0] point-B numbers on the shared shapes.
  **No point A** (it would need T₂₂, T₃₂). No continuous component.
- **Sites swapped** ([3,2,0] ≅ [3,1,0]*⊗det³): site 1 = pure rotation e₁ᵀR(z,(w₁,y*))ρ(Dy[w₂]) with
  **y* = w₂(θ₁₂−h)/(θ₁₂−3h)** (same as [3,1,0]); site 2 = dressed vacuum e₁ᵀ(1+bE₂₃)R(z,(w₁,y*)), no torus factor,
  dressing along E₂₃. Fit (18 points, not derived): b = (6125/5676)·w₂/((w₂−z₃)(y*−z₃)) at P1, w₁ = 5/2 — a lead for
  the open closed form of β/α.
- **Undo** (now of the site-1 rotation): φ = f2·ρ(n′)⁻¹; exact on all 10 HW kets; 3×3 identity
  (b″−a″c″)/(a″c″) = (w₂−y*)/(y*−z₂), E₁₃ exponent vanishes iff w₂ = y*.
- **Ranks over six w₂:** dressed point 12/15 ((0,2,3),(0,3,2) ≡ 0, (1,2,2) rank 1 of 2); β = 0: 15/15; control
  resolvent→y*, β = 0: blind spot μ₁ ≤ 1 (the dual of [3,1,0]'s μ₃ ≥ 2). LW frame: mirror ((2,2,1) rank 1).
- **HW states:** 10/10 reached (torus-completed; old site-1+site-2 localization had needed ⟨Ω|Ψ⟩ for [4,3,3]);
  **9/10 fully determined; [4,4,2] one direction short** (site weight (1,2,2)), not repaired by the LW frame; no
  point A to fall back on.
- Tooling: the agent had to associate the notebook controller through the wolfbook stdio bridge; K1 was restarted
  once (all other notebooks' live state gone, files intact). [4,2,0] not attempted.

## 11. Plan (2026-10-09, late): three avenues to complete the determination with "T-polynomial × det₆" overlaps

Paul: T-polynomial × det₆ is acceptable if the polynomial is pinned down generally and the Q-functions are explicit.
Why the adjoint is pure rotations: e₁ᵀE₁₂² = 0 there, so vacuum + first descendant is a rotation; for λ₁−λ₂ ≥ 2 it
is a chord off the orbit. General rotations g⊗g add only the torus character (§5.8), so new covectors must come
from new insertion points/products.

**A1 — Point A: structure, polynomial, explicit Q's** (`L2-310-two-rotations.wb`).
 (a) Is A's site-1 factor an n₊-dressed rotation, f1_A ∈ span{e₁ᵀ, e₁ᵀE₁₂, e₁ᵀE₂₃, e₁ᵀE₁₃}·R(z,(w₁,y_A))·D with
     D ∈ {Dy[w₂], Dy[y_A], 1}? Dressing coefficients vs w₂ (closed form like [3,2,0]'s b).
 (b) A's polynomial: c_A = adjoint c*(h,θ) on the column-shifted shapes, y_A = adjoint y*. Test whether
     f_A(B) = Σ c*_λ τ_{λ+(1,1)}(θ₁−3h) is a ladder-1 det₂ with insertion y_A^{−x/h} over the ADJOINT pole set
     n = 0..2 (minor at (θ₁−2h,θ₁−3h)), or over n = 1..3 (minor at (θ₁−3h,θ₁−4h)), using [3,1,0]'s Q's; also the
     general-(h,θ) adjoint coefficient monomials of `Two-rotations-det-overlap.md` §5.4 against c_A at the three
     numeric points.
 (c) Explicit Q-dictionary: solve the dual Baxter equation (opDual) for t₁(x), t₂(x) as Cramer/Wronskian
     expressions in two dual Q's at x+h, x, x−h, x−2h (t₃ = qdet), hence τ_λ(θ₁−3h) for all shapes via CBR;
     verify on all 225 states; write the simplified forms at x = θ₁−3h. Then any c·(1,τ_λ)·det₆ is explicit.

**A2 — Dual (site-2) recipe** (`L2-320-two-rotations.wb` live; then [3,1,0]). Shapes λ₁ ≤ 3 at θ₂+mh, m ∈ {−1,0,1,2}:
 localization on ⟨Ω⊗Ω| (site-2 local ⟺ rank 1 with site-1 factor Ω₁; known: t₂(θ₂+h), T₂₂(θ₂) for [3,2,0]),
 joint site-1 span of the non-local ones, exact product locus, factor identification (torus invariants, dressing
 type), undo of the rotation side, per-weight ranks, union with the θ₁−3h point → does [3,2,0] reach 10/10
 determined? Prediction by [3,2,0] ≅ [3,1,0]*⊗det³: the dual of [3,1,0]'s point A appears here.

**A3 — General dressing search** (`L2-310-two-rotations.wb`): product locus in enlarged spans — site-1 × site-2
 products T_λ(θ₁−3h)·T_μ(θ₂+mh), and shapes at θ₁−h, θ₁, θ₁+h — catalogue every product point: which side is a
 rotation, dressing type in {e₁ᵀ, e₁ᵀE₁₂, e₁ᵀE₂₃, e₁ᵀE₁₃}·RD, torus invariants, undo ranks; any continuous family?
 Test of the conjecture "every product point is rotation ⊗ n₊-dressed rotation; generic dressing ⇒ complete".

## 12. A2 outcome — [3,2,0] site-2 recipe: EMPTY (`L2-320-two-rotations.wb` cells 22–32, sonnet; saved, not committed)

- Site-2 locals at θ₂+mh (m = −1,0,1,2,−3; all nine shapes): exactly T₍₁,₁₎(θ₂+h) = t₂(θ₂+h) and T₍₂,₂₎(θ₂); rank-2
  operators T₍₂,₀₎(θ₂−h), T₍₁,₀₎(θ₂), T₍₃,₂₎(θ₂); joint site-1 span {Ω₁, ω, ω′} (dim 3; ω a rotation with the 2-plane
  {Ω₁,ω} inside the cone; ω′ rotRank 7).
- Product locus: S₂ (6 ops) → localizing span only, at all six w₂; all shapes at θ₂+h (10 ops) and at θ₂ (10 ops):
  localizing span only; both points together (19 ops): exact Solve times out, GF(1000003) Gröbner basis (436
  elements) shows all non-localizing coefficients nilpotent, affine dimension 3 → localizing span only (modular
  evidence). 22-op span untested. No continuous family.
- So the predicted dual of [3,1,0]'s point A does not appear; [3,2,0] stays at 9/10 determined ([4,4,2] short).
  Likely reason: only two site-2 localizers vs three at θ₁−3h.

## 13. A1 outcome — point A is a det₂ too; explicit Q-dictionary (`L2-310-two-rotations.wb` cells 43–65, sonnet; saved, not committed)

Notation: t = θ₁₂, W(a,b) = Q²(θ₁−ah)Q³(θ₁−bh) − Q³(θ₁−ah)Q²(θ₁−bh) (dual Q's), M = W(3,4) = minor₁.
- **(a) Dressing of A's site 1 is second order:** f1_A ∝ e₁ᵀ(1 + βE₁₂ + γ(E₁₃ + ½E₁₂E₂₃))·ρ(n_A)·ρdiag[Dy[w₂]], site 2 =
  e₁ᵀρ(n_A) (no torus), with exact rational closed forms β = (z₁−z₂)(w₁−z₁)/(2z₂(w₁−z₂))·w₂/(w₂−z₁) and
  γ = ¼(t−h)/(t−2h)·(z₁−z₃)(z₂−z₃)(w₁−z₁)/(z₂z₃(w₁−z₃))·w₂²/((w₂−z₁)(y_A−z₃)) (18 w₂, 10 one-parameter variations,
  28 random points). Not in the first-order span {e₁ᵀ, e₁ᵀE₁₂, e₁ᵀE₂₃, e₁ᵀE₁₃}·RD (rank 5). So the two points are:
  B = first-order dressing (E₁₂) ↔ pole set {0..3}; A = second-order dressing ↔ pole set {0,1,2}.
- **(b) A's polynomial = ladder-1 det₂ over the pole set n = 0,1,2** (one-site weights on that set, insertion
  y_A^{−(x−θ₁)/h}, rows (−⟨Q^jE⟩, ⟨Q^{j[−]}E⟩)), **divided by the [3,1,0] minor M = W(3,4)** (not W(2,3)):
  f_A ∝ det₂^{(0..2)}[y_A]/M, state-independent to 7e-48; end-to-end ⟨Φ_A|Ψ⟩/⟨x₀|Ψ⟩ ∝ det₆/Ψ_Q × det₂^{(0..2)}[y_A]/M to
  5e-45 (controls fail by 5.4, 96). The symbolic-y curve of det₂^{(0..2)}/M is supported on (T₁₁,T₂₁,T₂₂,T₃₁,T₃₂)
  and hits c_A only at y = y_A = −63/8. General form: c_A = (2h²(t−h)(t−2h)(t−3h)y⁴, 2h(t−h)(t−3h)y³, (t−h)y², (t−3h)y²,
  (2/(3h))y) (the adjoint's §5.4 monomials) on the shifted shapes, exact at all three numeric points; symbolically
  f_A(y) = K·det₂^{(0..2)}[y]/M with K = −96h⁸(h−t)(2h−t)(3h−t)²(4h−t). So §9's "A is not on the det₂ curve" was
  true only for the {0..3} pole set.
- **(c) Explicit dictionary.** From the dual Baxter equation (Cramer with any two dual Q's, written for (2,3)):
  t₁(x) = [ν₃(x)W(x+h,x−h) − t₃(x)W(x−2h,x−h)/(ν₃(x−h)ν₃(x−2h))]/W(x,x−h),
  t₂(x) = [ν₃(x)ν₃(x−h)W(x+h,x) − t₃(x)W(x−2h,x)/ν₃(x−2h)]/W(x,x−h), t₃ = χ₃ν₁(x)ν₂(x−h)ν₃(x−2h) (2e-46 on 225
  states). At θ₁−3h: ν₁ = 0 so t₃ = 0. At θ₁: W(0,1) ≡ 0 (Q(θ₁) = κQ(θ₁−h)), t₂(θ₁) = κν₃(θ₁−h)t₁(θ₁). At θ₁+h:
  t₂(θ₁+h) = χ₃ν₁(θ₁+h)ν₂(θ₁)Q^j(θ₁−h)/Q^j(θ₁). **Every fused eigenvalue at the site-1 point is ONE Wronskian:**
  τ_{(λ₁,λ₂)}(θ₁−3h) = C_λ·W(θ₁−(3−λ₁)h, θ₁−(4−λ₂)h)/W(θ₁−3h,θ₁−4h), C_λ = Π_{k=4−λ₁}^{3} kh(kh−t)·Π_{k=5−λ₂}^{4} kh(kh−t)
  (z-independent; verified symbolically in h,t and on all 225 states ≤ 3e-46; (3,3) ↦ W(0,1) = 0). E.g.
  τ₁₁ = 12h²(3h−t)(4h−t)W(2,3)/M, τ₃₂ = 72h⁵(h−t)(2h−t)(3h−t)²(4h−t)W(0,2)/M.
- Status for [3,1,0]²: two points, both det₆ × det₂^{(P)}[y_P]/M with (P, y_P) = ({0..3}, w₂(t−h)/(t−3h)) and
  ({0,1,2}, w₂(t−h)/(t−2h)); A alone fully determines all 10 HW states; all overlaps explicit in Q's.

## 14. A3 (revised) — the pole-set chain conjecture

Conjecture: for [λ₁,λ₂,0] the product points of the θ₁−λ₁h recipe are indexed by pole sets P_k = {0..k},
k = 2..λ₁, with y_k = w₂(θ₁₂−h)/(θ₁₂−kh), overlap det₆ × det₂^{(P_k)}[y_k]/M, dressing of order λ₁−k+1 on the
partner site; k = 1 (y = w₂, the ω = y degeneracy) is excluded; only the k < λ₁ points (higher dressing) give
complete per-weight spans. Tests: (i) [3,1,0]: all pole subsets of {0..3} → which curves hit product points (expect
exactly {0..3}→B, {0,1,2}→A); (ii) [3,2,0]: same with its Q's and minor (expect only {0..3}; see what the {0,1,2}
curve is); (iii) [4,1,0] (dim 24, 576 states): product locus at θ₁−4h, predicted points k = 2,3,4, dressing orders
3,2,1, per-weight spans and HW determination (15 HW states).

## 15. A3 (i),(ii) outcome — the chain holds on [3,1,0]; on [3,2,0] the A-curve is localizing (sonnet; saved, not committed)

- **[3,1,0]** (`L2-310-two-rotations.wb` cells 66–70): all 11 pole subsets P ⊆ {0..3}, |P| ≥ 2, give det₂^{(P)}[y]/W(3,4)
  curves in span{1, τ_λ(θ₁−3h)} (5e-46). Exact rank-1 loci in y (gcd of the 2×2 minors): **only {0..3} → B at
  y = −63/43 and {0,1,2} → A at y = −63/8**; {0,1} is trivially localizing (support T₂₂,T₃₂); all other subsets never
  produce a product covector. So exactly two points, as conjectured.
- **[3,2,0]** (`L2-320-two-rotations.wb` cells 33–41; Q-side ported, 225-state Baxter solve 8 s): **vacuum triple
  (θ−4h, θ−3h, θ−2h)** (Ψ_Q = Π_s W_s(3,4)·Q³(θ_s−2h); the [3,1,0] triple (4,3,1) fails by 22), **minor W(3,4)**
  (W(2,3), W(2,4), W(3,5) fail), end-to-end ⟨Φ_B|Ψ⟩/⟨x₀|Ψ⟩ ∝ det₆/Ψ_Q × det₂^{(0..3)}[y*]/W(3,4) to 1.2e-44. Pole-subset
  scan: {0..3} → B at −63/43 only; **the {0,1,2} curve = the [3,1,0] A-curve with the T₂₂, T₃₂ entries removed (they
  vanish on ⟨Ω|), i.e. supported on {T₁₁,T₂₁,T₃₁} = the localizing span → site 2 = Ω₁, never a nontrivial product.**
  Hence no point A for [3,2,0] (explains §12), 9/10 stays.
- Rep-data pattern so far: vacuum triple (θ−(λ₁+1)h, θ−λ₁h, θ−λ₂h), minor W(λ₁, λ₁+1), det₂ pole sets ⊆ {0..λ₁}.

## 16. A3 (iii) outcome — [4,1,0]: the chain holds, three points (`L2-410-two-rotations.wb`, 26 cells, sonnet; saved, not committed)

Caveat: the notebook controller could not be associated, so the cells were evaluated through the kernel from
their saved sources; outputs are not stored in the file (cell 26 markdown holds the numbers; rerun in VS Code to
regenerate, ~1 min total). 53 Boolean checks pass.
- **Rep data:** localizers at θ₁−4h T₁₁,T₂₂,T₃₃,T₄₃ (T₄₄ ≡ 0 on ⟨Ω|); joint site-2 span {Ω₁, ω₂}, ω₂ a rotation; Q degrees
  deg q^j = μ_j, deg q_i = 8−μ_i (576 states, 13 s); **vacuum triple (θ−5h, θ−4h, θ−h)**, det₆ with poles n = 0..4 (6e-43;
  other triples fail ≥ 1.7); **minor W(4,5)** (next best 0.82). [4,1,0]² = 14 HW states in 11 irreps ([8,2,0],[8,1,1],
  [7,3,0],2[7,2,1],[6,4,0],2[6,3,1],[6,2,2],[5,5,0],2[5,4,1],[5,3,2],[4,4,2]; dims sum to 576), not 15.
- **Product locus** (exact, 276 quadrics, 4.5 s): localizing span + **exactly three isolated points**, site 2 = torus-free
  rotated vacuum e₁ᵀρ(n_k) with **y_k = w₂(θ₁₂−h)/(θ₁₂−kh), k = 2,3,4** (−63/8, −63/43, −21/26); no k = 1 point; no
  continuous component; stable over 24 w₂'s. All 26 pole subsets scanned: only {0,1,2}, {0..3}, {0..4} hit
  product points, each exactly once at y_k. **End-to-end:** ⟨Φ_k|Ψ⟩/⟨x₀|Ψ⟩ ∝ det₆/Ψ_Q × det₂^{(0..k)}[y_k]/W(4,5) on
  all 576 states to 6e-43 for k = 2,3,4 (controls fail by 10–10⁵).
- **Dressing of site 1** = e₁ᵀ·(polynomial in E₁₂,E₁₃,E₂₃)·ρ(n_k)ρ(Δ_{w₂}); minimal monomial sets: k=4 {1,E₁₂,E₁₂²};
  k=3: 7 monomials (E₁₃-degree 1); k=2: 9 monomials (E₁₃-degree 2). So the E₁₃-degree is λ₁−k, not the word length.
- **Undo (exact on all 14 HW kets) and ranks over six w₂** (18 weights, 6 of dim 2):

| k | pole set | directions | irreps determined | reached (torus) |
|---|---|---|---|---|
| 2 | {0,1,2} | **24/24** | **11/11** | 14/14 |
| 3 | {0..3} | 21/24 ((0,1,4),(1,0,4) ≡ 0; (1,1,3) 1/2) | 11/11 | 14/14 |
| 4 | {0..4} (= B) | 16/24 | 9/11 ([5,3,2], [6,2,2] one direction each) | 14/14 |

- **Corrected conjecture:** points ↔ P_k = {0..k}, k = 2..λ₁, y_k = w₂(θ₁₂−h)/(θ₁₂−kh), overlap det₆ × det₂^{(P_k)}[y_k]/
  W(λ₁,λ₁+1), site-1 dressing of E₁₃-degree λ₁−k; the **k = 2 point (adjoint pole set) spans every weight space and
  determines all HW states**; higher k lose directions (k = λ₁ loses the most). Holds on [2,1,0] (k=2 only),
  [3,1,0] (k=2,3), [4,1,0] (k=2,3,4). Fails to produce k < λ₁ points when the shapes they need vanish on ⟨Ω|
  ([3,2,0]: T₂₂,T₃₂ ≡ 0 → only k=3).
