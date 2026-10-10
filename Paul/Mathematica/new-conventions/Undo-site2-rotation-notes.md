# Undoing the site-2 rotation of the two-rotation family, and the pole-set chain of determinant product states

gl(3) XXX, L=2, new conventions. Research notes on the 2026-10-09 (late) thread. Audience: Paul Ryan and future Claude sessions. This file is the Markdown companion of `Undo-site2-rotation-notes.tex`/`.pdf`; both have the same content.

Primary sources: `site1-Qform/7-undo-site2-rotation-span.md` (spec + running log, §1–§16, called "the spec" below), `HANDOFF.md` (entry "NEWEST (2026-10-09, late)"), `Two-rotations-det-overlap.md`, `L2-310-two-rotations.md`, `old-conventions/XXX-CG/xxx-210-cg-test.tex`, and the result cells of the notebooks listed in Appendix A. All paths are under `Paul/Mathematica/new-conventions/` unless stated. Cell numbers are 1-based positions in the `.wb` file (the Python `cells` list index plus one).

**Evidence tags.** Every claim carries one of:
- **[exact]**: a symbolic identity, or an exact-rational computation (ranks, `Solve`), or agreement at 60-digit working precision with a stated residual (the residual is quoted);
- **[numeric]**: true at the numeric point(s) stated, residual in digits quoted, not proved;
- **[conjecture]**: pattern or prediction not (fully) verified; the text says what is verified and what is not;
- **[modular]**: Gröbner evidence over GF(1000003) only.

No notebook, no other `.md`/`.tex`, `HANDOFF.md`, `CLAUDE.md` was edited in writing these notes, nothing was evaluated, and nothing was committed. One independent check was run outside the notebooks while writing (sympy, §3): the 3×3 closed forms of the spec were re-derived symbolically and agree.

---

## 0. Summary

**The question (Paul).** The two-rotation family of `Two-rotations-det-overlap.md` §5.6b/§5.8 gives product covectors ⟨v₁⊗v₂| with a Q-determinant overlap, where v₂ = ⟨0|R is a rotation of the vacuum. Undo the rotation, ⟨v₁|R⁻¹, project on the weight spaces of the site-1 representation V and expand in all free parameters. Do these projections span every weight space of V, so that all highest-weight (HW) CG coefficients ⟨μ ⊗ top | ν, HW⟩ are reachable? The 2026-09-23 families (`old-conventions/XXX-CG/xxx-210-cg-test.tex`) reached only site-1 weights with μ₃ ≤ λ₂, so they missed the adjoint singlet, [4,2,2] and [3,3,2] of [3,1,0]², and so on.

**The mechanism.** In the diagonal frame, v₂ = e₁ᵀρ(n) with n ∈ N₊ (the vacuum Gauss element), and ρ(n)⊗ρ(n) fixes every HW ket. So the undo is exact at finite twist: ⟨v₁⊗v₂|HW⟩ = ⟨v₁ρ(n)⁻¹ ⊗ e₁ᵀ|HW⟩. The undone covector is φ = e₁ᵀ g ρ(nΔn⁻¹) = δ^Λ e₁ᵀ g′ ρ(n′), with Δ = diag(1/(ω−z_k)) the resolvent torus at the det₆ point ω = w₂, g the site-1 dressing (g = 1 for the adjoint) and n′ = Δ⁻¹nΔn⁻¹ ∈ N₊. The (1,3) entry b′ of n′ vanishes iff ω = y (y = the ladder-1 insertion point), and b′/(a′c′) = (y−ω)/(y−z₂). The old families sat exactly at ω = y (the "p₃ relation"), the new family has ω = y (θ₁₂−kh)/(θ₁₂−h) ≠ y for k ≥ 2, which opens the E₁₃ direction and with it all weights with μ₃ > λ₂ **[exact]**.

**Results per representation** (P1 numerics, tags in the body).

| rep | points (pole set P_k, k) | directions spanned | HW states reached (torus-completed) | HW states fully determined |
|---|---|---|---|---|
| [2,1,0] adjoint | k=2, {0,1,2} (pure rotations) | 8/8 | 6/6 incl. the singlet | 6/6 (5/6 at fixed w₁ and fixed torus) |
| [3,1,0] | B: k=3, {0..3} | 12/15 | 10/10 | 9/10 ([4,2,2] short by one direction) |
| [3,1,0] | A: k=2, {0,1,2} | 15/15 | 10/10 | 10/10 |
| [3,2,0] | B only: k=3, {0..3} | 12/15 | 10/10 | 9/10 ([4,4,2] short by one direction) |
| [4,1,0] | k=2 | 24/24 | 14/14 | 11/11 irreps |
| [4,1,0] | k=3 | 21/24 | 14/14 | 11/11 irreps |
| [4,1,0] | k=4 (= B) | 16/24 | 14/14 | 9/11 irreps ([5,3,2], [6,2,2] one direction short each) |

(Reached = nonzero overlap family; determined = the undone covector has full rank on the site-1 weight space ν−Λ carrying the CG vector of ν. See §2.)

**The chain.** For [λ₁,1,0] the product points of the recipe ⟨Ω⊗Ω|(c₀ + Σ c_λ T_λ(θ₁−λ₁h)) are indexed by the initial segments P_k = {0..k}, k = 2..λ₁, at y_k = w₂(θ₁₂−h)/(θ₁₂−kh), with overlap det₆ × det₂^{(P_k)}[y_k]/W(λ₁,λ₁+1), and a site-1 dressing of E₁₃-degree λ₁−k. Exact on [3,1,0] (all 11 subsets scanned) and [4,1,0] (all 26 subsets scanned). On [3,2,0] the k=2 curve collapses into the localizing span, so only k=3 survives. The empirical criterion "a non-trivial point k exists iff λ₂ < k ≤ λ₁" is derived in §9 from three observed rules; it is verified on four representations and is a **[conjecture]** beyond them (predictions for [4,2,0], [λ,1,0]).

**What is explicit now.** Every overlap of a chain point is a product of two Q-determinants (a global det₆ and a ladder-1 det₂), and for [3,1,0] every fused transfer-matrix eigenvalue at the site-1 point is a single Wronskian over a minor (§6).

**What is open.** (i) [3,2,0]-type reps (λ₂ ≥ 2, d = λ₁−λ₂ = 1 in particular): no point with a complete weight span exists in the θ₁−λ₁h recipe, and the site-2 recipe is empty; [4,4,2] stays one direction short. (ii) [4,2,0] not attempted. (iii) Closed forms of the dressing coefficients (known only for [3,1,0] point A, and a fit for [3,2,0]). (iv) The last layer between the covector counts and CG numbers, the untwisting/contamination step, has not been run for this family (§10). (v) The [4,1,0] outputs are not stored in the file (§12).

---

## 1. Setting and conventions

**Conventions (new, 2026-10-01).** Lax ℒ_ij(u) = (u−θ_s)δ_ij + hE_ij; T(u) = ℒ⁽¹⁾ℒ⁽²⁾ with twist M acting on the right (T^M = T·M); t_a = sum of principal quantum minors, k-th column at u−(k−1)h; CBR fusion T_λ(x) = det_{i,j} t_{λ′_j+i−j}(x+(i−1)h) with shapes labelled by row lengths (λ₁,λ₂), t₀ = 1, t_a = 0 outside 0..3 (e.g. T₍₂,₀₎ = t₁t₁^{[+]} − t₂^{[+]}, T₍₂,₁₎ = t₂t₁^{[+]} − t₃^{[+]}). Generators: the plain GT generators `JJ` (not the Hermitian EH ones), so everything is exact rational. The GT basis has the HW vector e₁ first; covectors are rows acted on from the right; the ket e₁ is annihilated by E_ij with i<j (checked in each notebook), so E_ij with i<j raise. Sources: `L2-rotation-sov.wb` cells 2–21, `L2-310-two-rotations.wb` cells 2–7, `L2-310-two-rotations.md`.

**Notation.** t := θ₁₂ = θ₁−θ₂ in this file (so t = 31/35 at P1, θ₁₂−h = 27/70, θ₁₂−2h = −8/70, θ₁₂−3h = −43/70). Caution: `Two-rotations-det-overlap.md` §1 and `Clean-CG-combinations.md` write t := θ₁₂/h = 62/35; the [3,1,0]/[4,1,0] notebooks and spec §13 use t = θ₁₂ in formulas such as kh(kh−t). The two conventions differ by a factor h; formulas below are written with θ₁₂ and h explicit where ambiguity could arise.

**Numeric point P1.** h = 1/2, θ = (2/7, −3/5), ket twist z = (3/2, 5/4, 2/3), bra twist w = (5/2, 7/3, 3/5) so w₁ = 5/2, w₃ = 3/5. Precision 60 digits for Q-side numbers. Other points: h = 1/3 and (θ = (1/5,−4/9), z = (3/2,7/6,2/3), w = (7/4,−2/5,9/7)) (`L2-310-two-rotations.md`, cell 12 of `L2-310-two-rotations.wb`). Caution: at P1 z₃ = 1/z₁ (coincidence flagged in `Rotations-SoV-Q`). The six values w₂ ∈ {7/3, −8/9, 4/27, 3, −1/2, 11/7} are used throughout (with w₁ = 5/2, w₃ = 3/5).

**Diagonal frame, GT basis.** Ket twist diag(z); the SoV vacuum covector on both sites is Ω⊗Ω with the one-site vacuum

Ω_w = e₁ᵀ e^{r₁E₁₂} e^{r₂E₁₃} e^{r₃E₂₃}  (`Om1[z,w]`, `oneSiteVac`; `Rrot[z,w]` = the 3×3 product of exponentials in the rep),

with the representation-independent Gauss parameters (w₃ drops out; the second slot of w is the variable y)

r₁ = (w₁−z₁)/(w₁−z₂),
r₂ = (w₁−z₁)(z₁−y)(z₂−z₃) / ((z₁−z₂)(w₁−z₃)(y−z₃)),
r₃ = (w₁−z₂)(z₂−y)(z₁−z₃) / ((z₁−z₂)(w₁−z₃)(z₃−y)),

with w = (w₁, y, w₃) (spec §2; `rOm[zv,wv]` in `L2-310-two-rotations.wb` cell 4, in `L2-410-two-rotations.wb` cell 4). Ω⊗Ω is the exact left eigencovector of t_a(G′) with G′ = P_aux G(w) P_aux⁻¹ (checked, `L2-320-two-rotations.wb` cell 4). Ω_w = e₁ᵀρ(n) with n ∈ N₊ the 3×3 matrix [[1, r₁, r₂+r₁r₃],[0,1,r₃],[0,0,1]]. **[exact]** (cell 6 / 98 of `L2-rotation-sov.wb`; cell 26 of `L2-310-two-rotations.wb`: v₂ = e₁ᵀρ(n) and ρ(n) = Rrot for all six w₂).

**Resolvent torus.** The resolvent ρ((ω − Gᵀ)⁻¹) rotates the ket-twist frame; in the diagonal frame it is the torus element Δ = ρdiag[1/(ω−z_k)], δ_k = 1/(ω−z_k) (`Dy[ω]`, `ρdiag`, `L2-rotation-sov.wb` cell 10; cell 111 builds the partner as `g . ρdiag[Dy[w2G]]`). Weight components of a covector of weight μ (relative to e₁ᵀ) are multiplied by δ^μ under ρ(Δ). The torus of the ket twist, ρ(D)⊗ρ(D), commutes with the twisted transfer matrices and multiplies an eigenstate of total weight ν by the character χ_ν(D) = Π D_i^{ν_i} (`Two-rotations-det-overlap.md` §5.8, `L2-rotation-sov.wb` cell 108, exact at 60 digits and exact-rational at cell 121).

**The two-rotation family V (adjoint; `Two-rotations-det-overlap.md` §5.6b, §5.8).** With y free, w₂ slaved to it,

ω = w₂(y) = y (θ₁₂−2h)/(θ₁₂−h),
v₂(y) = Ω_{(w₁,y)} = e₁ᵀρ(n(y)),  v₁(y) = v₂(y) ρ(Δ(y)),  Δ = diag(1/(ω−z_k)),

and for every y the overlap is ⟨v₁⊗v₂|Ψ_B⟩/⟨x₀|Ψ_B⟩ ∝ det₆N[Q_B, Q_vac,(w₁,w₂(y))]/Ψ_Q × det₂^{(1)}[Q_B; y]/minor₁ (state independence 10⁻⁵⁰ on 64 states at y = 3, −1/2; at fresh w₁ = −7/3, y = 7/5 and w₁ = 9/4, y = −2/3, cell 109) **[numeric]**. The torus orbit of this family is the 3-dimensional variety

V = {(g ρ((w₂(y(g)) − Gᵀ)⁻¹), g) : g a rotation of x₀},

because w₁ is a torus direction and y is a torus invariant of the rotation (exact Jacobian ranks 3/4/3, cell 108); for a random g (lower unipotent · torus · upper unipotent) the exact `Solve` gives exactly one y(g) at both w₁ tried (y(g) = 1546301/1368484) **[exact]** (cell 111). No continuous κ-family exists (§5.8 of that note). The family V is "any rotation g on site 2, one fixed resolvent partner on site 1, times the torus orbit": torus = character on the overlap.

For general reps the analogous points are products of the same type with y_k = w₂(θ₁₂−h)/(θ₁₂−kh), i.e. ω = w₂ = y_k(θ₁₂−kh)/(θ₁₂−h) (§5, §9).

**Q-side objects used below** (definitions as in `Two-rotations-det-overlap.md` §1 and the Q-side cells of the [3,1,0] notebook): dual Q's Q^j(x) = z_j^{x/h}q^j(x), direct Q_i(x) = z_i^{−x/h}q_i(x); W(a,b) = Q²(θ₁−ah)Q³(θ₁−bh) − Q³(θ₁−ah)Q²(θ₁−bh) (dual Q's); minor₁ = M; det₆ the 6×6 (adjoint) or 8-pole (rep-generic: poles n = 0..λ₁ on both ladders, insertions 1 and e^{2πix/h}) bracket determinant; ladder-1 one-site weights w_n^{P} = 1/Π_{m∈P∖{n}}(m−n)h and

det₂^{(P)}[y] = Σ_{n,m∈P} w_n^{P} w_m^{P} y^{n+m} [Q²(θ₁−nh)Q³(θ₁−(m+1)h) − Q³(θ₁−nh)Q²(θ₁−(m+1)h)],

i.e. det_{j∈{3,2}}[−⟨Q^jE⟩₁, ⟨Q^{j[−]}E⟩₁] with E(θ₁−nh) = y^n (`L2-310-two-rotations.wb` cell 67). For the adjoint P = {0,1,2} and "det₂^{(1)}" is the same object.

---

## 2. The question and why the undo is exact at finite twist

**Setup.** e₁ᵀ is the HW covector of V (weight Λ). Every vacuum rotation in the family has v₂ = s e₁ᵀρ(n), n ∈ N₊ (upper unipotent), because Ω_w is literally e₁ᵀ times the three exponentials. For n ∈ N₊, ρ(n)⊗ρ(n) fixes every HW ket of V⊗V. Hence, with μ = ν − Λ and φ := v₁ρ(n)⁻¹,

⟨v₁⊗v₂|ν,HW⟩ = s ⟨φ ⊗ e₁ᵀ|ν,HW⟩ = s Σ_{m∈V_μ} φ_m CG(ν,HW | m, top).

So the HW CG coefficients at site-1 weight μ are reachable iff the weight-μ projections of the undone covectors φ span V_μ. **[exact]**, a finite-twist statement; it was the "Result `res:nplus`" of the old thread (`xxx-210-cg-test.tex`, Result 5, there for the old families) and is re-checked here: cell 115 of `L2-rotation-sov.wb` (adjoint, all 6 HW kets, (ρ(n)⊗ρ(n))HW = HW and the identity ⟨v₁⊗v₂|HW⟩ = s⟨φ⊗e₁ᵀ|HW⟩), cell 27 of `L2-310-two-rotations.wb` ([3,1,0], 10 HW kets), cell 12 of `L2-320-two-rotations.wb` ([3,2,0], 10 HW kets; here the rotated factor sits on site 1 so φ lands on site 2, "flip symmetry of V⊗V"), cell 12 of `L2-410-two-rotations.wb` ([4,1,0], 14 HW kets, all k, all six w₂). The untwisting limit enters only through the eigenstates, not through the covectors, so no untwisting is needed for the count.

**Reached versus determined.**
- *Reached* (a HW state ν): the overlap ⟨v₁⊗v₂|ν,HW⟩ is not identically zero over the family. The counts "HW rank" are ranks of the matrix M[i,ν] = ⟨φ(p_i)⊗e₁ᵀ|HW_ν⟩ over sample points p_i of the family; rank = number of HW kets (10 for [3,1,0], etc.) means all are reached and independent.
- *Determined* (an irrep ν, "fully determined" in the notebooks): the rank of φ on the weight space V_{ν−Λ} (HW frame) equals dim V_{ν−Λ}. Then the full CG vector ⟨m ⊗ top|ν,HW⟩, m ∈ V_{ν−Λ}, is read off from the functionals φ|_{ν−Λ}. If the rank is smaller, only a subspace of linear combinations of the CG components is accessible.
- For multiplicity-2 irreps the per-irrep column rank is also reported (e.g. "(1,1,1,2,1)" in the adjoint).

**Why the torus is needed for the joint HW rank, but not for per-weight spans.** The torus ρ(D)⊗ρ(D) multiplies the weight-μ components of a covector by d^μ/χ(D) (and column ν of the HW matrix by χ_ν(D)). A per-weight span involves a single weight, so the torus does not change it, and the per-weight spans of the whole 3-dim family V equal those of the 1-parameter slice φ(y) at fixed w₁; varying w₁ must not change any rank (checked, adjoint cell 120–121). The joint HW rank mixes columns of different weights ν, and at fixed D the slice is rank-deficient because of cross-weight linear relations among a′, c′, b′, a′c′ as functions of y (adjoint: rank 3 of 4, cell 121); the torus separates columns of different weight through χ_ν. Per-irrep column ranks are already full at fixed torus (cells 120, 14 of the 320 notebook, 14 of the 410 notebook).

**Which frame.** The undo is stated in the HW frame. A lowest-weight (LW) frame also exists (cells 34–36 of `L2-310-two-rotations.wb`, cell 18 of `L2-320-two-rotations.wb`): v₂ = s e_lowᵀρ(m) with m ∈ N₋ from the exact factorization w₀⁻¹n = b₊m, and ρ(m)⊗ρ(m) fixes all LW kets. The two frames measure different vectors (HW vs LW of each irrep) with unrelated normalizations, so they do not combine; see §5, §7.

---

## 3. The 3×3 mechanism (representation independent)

**Setup.** Let n ∈ N₊ be the vacuum Gauss element above, Δ = diag(δ_k), δ_k = 1/(ω−z_k), where ω is the resolvent point (the second twist parameter w₂ of the det₆ vacuum). Then

φ = v₁ρ(n)⁻¹ = e₁ᵀ g ρ(nΔn⁻¹) = e₁ᵀ g ρ(Δ) ρ(n′),  n′ := Δ⁻¹ n Δ n⁻¹ ∈ N₊,

where g ∈ U(n₊) is the site-1 dressing (g = 1 for the adjoint, where v₁ = v₂ρ(Δ) exactly). Since e₁ᵀ is a torus eigenvector, e₁ᵀρ(Δ) = δ^Λ e₁ᵀ, and the dressing coefficients get rescaled by torus weights (for the first-order dressing, e₁ᵀ(α+βE₁₂)ρ(Δ) = δ^Λ e₁ᵀ(α + β(δ₂/δ₁)E₁₂)).

Write n′ = [[1, a′, b′],[0,1,c′],[0,0,1]] (matrix entries; as a product of exponentials n′ = e^{c′E₂₃}e^{b′E₁₃}e^{a′E₁₂}). With n = [[1,r₁,r₂+r₁r₃],[0,1,r₃],[0,0,1]] (spec §2):

a′ = r₁(δ₂−δ₁)/δ₁,  c′ = r₃(δ₃−δ₂)/δ₂,  b′ = r₁r₃(δ₁−δ₂)/δ₁ − (r₂+r₁r₃)(δ₁−δ₃)/δ₁.

**Results [exact]** (spec §3(i)(ii) predicted; verified symbolically in `L2-rotation-sov.wb` cell 117 for the adjoint with y symbolic; for [3,1,0] in `L2-310-two-rotations.wb` cell 31 with symbolic y, w₁, θ₁₂, h and also symbolic z; for [3,2,0] in `L2-320-two-rotations.wb` cell 17 for n″ = nΔ⁻¹n⁻¹Δ; re-derived by me with sympy while writing, outputs 0 for the lower entries, the formulas for a′, c′, b′ and both identities below):

1. **b′ = 0 ⟺ ω = y.** Using the r's, r₂/(r₁r₃) = (z₃−z₂)(y−z₁)/((z₁−z₃)(y−z₂)), and b′ = 0 ⟺ r₂/(r₁r₃) = (z₃−z₂)(ω−z₁)/((z₁−z₃)(ω−z₂)). Since ω ≠ y on the family, b′ ≠ 0 (the only zero at the numeric point is the isolated y = 0, cell 123).
2. **b′/(a′c′) = (y−ω)/(y−z₂)**, non-constant in y and independent of w₁. This is what makes the two weight-(1,1,1) monomials b′ e₁ᵀE₁₃ and a′c′ e₁ᵀE₂₃E₁₂ independent as functions of y (rank 2 on the 2-dim adjoint octet weight space).
3. For [3,1,0] with ω = w₂ and y = y*: b′/(a′c′) = (y*−w₂)/(y*−z₂), the same identity (cell 31). For [3,2,0] (undone element n″ = nΔ⁻¹n⁻¹Δ, exponents r₁″ = a″, r₃″ = c″, r₂″ = b″−a″c″): (b″−a″c″)/(a″c″) = (w₂−y*)/(y*−z₂), i.e. b″/(a″c″) = (w₂−z₂)/(y*−z₂) (cell 17 prints (5−4W)/(5−4Y) at z₂ = 5/4); the E₁₃ exponent vanishes iff w₂ = y*.

**The old p₃ relation.** In the 2026-09-23 families the covector was e₁ᵀ g₀ … with p₃ = p₂(z₃−z₁)/(p₁(z₂−z₃)); the element G⁻¹g₀Gg₀⁻¹ has zero (1,3) entry by that relation (`xxx-210-cg-test.tex`, Result `mech`, and spec §3(i)). This is the same algebra with Δ = diag(z), i.e. exactly the b′ = 0 degenerate case. The "paper-like" assignment (resolvent point at the vacuum's own y) would give b′ ≡ 0 as well. The new family has ω/y = (θ₁₂−kh)/(θ₁₂−h): the mismatch between the det₆ twist and the ladder insertion is what opens E₁₃.

**The (0,1,2) component (adjoint) [exact].** The weight components of φ/δ^Λ = e₁ᵀρ(n′) are fixed constant covectors times monomials in (a′,b′,c′): (2,1,0): 1; (1,2,0): a′; (2,0,1): c′; (1,1,1) (dim 2): b′ e₁ᵀE₁₃ + a′c′ e₁ᵀE₂₃E₁₂; (0,2,1): a′²c′ u + a′b′ u′; (1,0,2): b′c′ u″; (0,1,2) (singlet direction): **3 b′(a′c′ − b′)**, exact in y (cell 118, fit exact and verified symbolically). Spec §3(iii) predicted −b′(b′+a′c′) up to normalization; that sign/shape was wrong (a′²c′² does not occur); cell 123 item 3 corrects it. The conclusion (nonzero for generic y) is unchanged. The other components agree with spec §3(iii) (cell 123).

**Untwisting leading orders [exact].** With z_k = 1 + εφ_k: a′, c′ = O(ε); b′, a′c′ = O(ε²) with ratio b′/(a′c′) → (y−ω)/(y−1), still non-constant; (0,1,2) component O(ε⁴) (spec §3; computed in cell 119 for the full 8-component φ with φ = (3/7, −5/4, 2/9), y symbolic, series to ε⁷).

---

## 4. Adjoint [2,1,0]² (`L2-rotation-sov.wb` §12, cells 113–123)

**Family.** v₂(y) = e₁ᵀρ(n(y)), v₁(y) = v₂(y)ρ(Δ(y)), w₁ = −7/3 (and 9/4, 5/2 for consistency), ω = w₂(y) = y(θ₁₂−2h)/(θ₁₂−h), 12 rational y's for the numeric ranks: `u12ys` = {7/5, −2/3, 3, −1/2, 11/9, 4/7, −9/5, 13/6, −3/11, 8/3, 1/6, −7/4} (cell 116).

**Per-weight ranks [exact]** (cells 116, 122, 123). V has weights (0,1,2),(0,2,1),(1,0,2),(1,1,1) [dim 2],(1,2,0),(2,0,1),(2,1,0). Control = resolvent point set to y (ω = y).

| weight μ | dim V_μ | rank numeric (w₁=−7/3, 12 y) | rank exact in y | rank union over 3 w₁ | rank control ω=y | untwisting leading order |
|---|---|---|---|---|---|---|
| (0,1,2) | 1 | 1 | 1 | 1 | **0** | ε⁴ |
| (0,2,1) | 1 | 1 | 1 | 1 | 1 | ε³ |
| (1,0,2) | 1 | 1 | 1 | 1 | **0** | ε³ |
| (1,1,1) | 2 | 2 | 2 | 2 | **1** | ε² (y-rank 2) |
| (1,2,0) | 1 | 1 | 1 | 1 | 1 | ε¹ |
| (2,0,1) | 1 | 1 | 1 | 1 | 1 | ε¹ |
| (2,1,0) | 1 | 1 | 1 | 1 | 1 | ε⁰ |

So all 8 directions of V are spanned, no component vanishes identically (**[exact]**, rank-tested exactly in y by numerator coefficient vectors, cell 116), unchanged for w₁ = 9/4, 5/2 and for the union of the three w₁ (cell 120, check 12.9). The control reproduces exactly the old blind spot: (0,1,2) and (1,0,2) vanish identically and the 2-dim weight space is hit only once (rank 1 of 2).

**HW count [exact]** (cells 120, 121, 123). HW kets of adj⊗adj: [4,2,0] (the 27), [4,1,1] and [3,3,0] (the two decuplets), [3,2,1] twice (two octets), [2,2,2] (singlet); 6 kets.

| quantity | value |
|---|---|
| HW rank at fixed w₁ and torus D=1 (w₁ = −7/3, 9/4, 5/2), 12 y's | 5, 5, 5 |
| per-irrep column ranks, fixed (irreps in the order 27, [4,1,1], [3,3,0], octets, singlet) | 1,1,1,2,1 (singlet column nonzero) |
| union over the three w₁ | 6 |
| torus-completed (12 random (y, D) pairs, rows built directly from ⟨v₁ρ(D)⊗v₂ρ(D)|HW⟩) | **6**, per irrep 1,1,1,2,1 |
| torus = character on the HW overlaps | exact residual 0 (cell 121) |
| control ω = y: rank at fixed torus per w₁ | 4,4,4 (per irrep 1,1,1,1,0 at w₁ = −7/3) |
| control ω = y: torus-completed | 4, per irrep 1,1,1,1,0 (singlet missing, one octet direction lost) |
| rank of {a′, c′, b′, a′c′} as functions of y | 3 (one relation) |

The rank-5 at fixed w₁ is a deficiency of the 1-parameter slice, not of the family: the single null vector lives on the columns [4,1,1], [3,3,0] and the two octets (different weights except the octet pair), and corresponds to the linear relation among a′, c′, b′, a′c′. Harmless for reading CG coefficients state by state (each state has a definite ν).

**Comparison with 2026-09-23** (`xxx-210-cg-test.tex`, Results 1–2, 4): the old families reached 6 of the 8 site-1 weight directions ((1,0,2) and (0,1,2) vanish identically in v_new), 5 of 6 HW states (singlet [2,2,2] missed since it needs site-1 weight (0,1,2)), and their leading orders were ε⁰,ε¹,ε¹,ε²,ε³ at (2,1,0),(1,2,0),(2,0,1),(1,1,1),(0,2,1). Here the singlet direction is reached and the family reaches 6/6. The old note explained the misses as the rule μ₃ ≤ λ₂ (Result "why the rule ν₃ ≤ λ₂": each quantum minor moves at most one unit into slot 3 whatever the twist, and the rotation has no E₁₃ component), and said a covector with an E₁₃ component would reach all 6, with (G⁻²g₀G²g₀⁻¹)₁₃ = −p₂(z₁−z₂)(z₁−z₃)/z₁² ≠ 0 as the hypothetical fix (Result "a hypothetical fix"). The new family realizes that fix: b′ ≠ 0 is the E₁₃ exponent.

---

## 5. [3,1,0]² (`L2-310-two-rotations.wb`)

**Rep data.** Λ = (3,1,0), dim 15, 225 states, 12 distinct weights (dim 1 except (2,1,1),(1,2,1),(1,1,2) of dim 2) **[exact]** (cell 33; the spec §7 plan said 10, corrected by §8 and cell 33). HW states of [3,1,0]²: 10 in 8 irreps: [6,2,0], [6,1,1], [5,3,0], 2×[5,2,1], [4,4,0], 2×[4,3,1], [4,2,2], [3,3,2] (cell 27). Site-1 weight carrying the CG vector of each: (3,1,0), (3,0,1), (2,2,0), (2,1,1), (1,3,0), (1,2,1), (1,1,2), (0,2,2) respectively (cell 41).

**Product locus at θ₁−3h [exact]** (`L2-310-two-rotations.md`, cells 5–7): 105 quadrics, exact `Solve` 0.5 s. Site-1 localizers T₍₁,₁₎, T₍₂,₂₎, T₍₃,₂₎; T₍₃,₃₎ vanishes on ⟨Ω|; all other shapes with λ₁ ≤ 3 have rank 2 with Ω₁ in their site-2 span; joint site-2 span {Ω₁, ω₂}, ω₂ a rotation, pencil Ω₁+κω₂ inside the rotation orbit (shapes with λ₁ = 4 add nothing, cell 8). The locus is the localizing span plus two isolated points A and B. A mixed-point span (θ₁−3h shapes plus λ₁ ≤ 2 shapes at θ₁−2h, 3549 quadrics, `Solve` 29 s) gives no new point.

The rotation test (exact, linear; `rotRank`): φ ∈ GL(3)e₁ᵀ iff rank{φ, φE_ij (9)} = 4; generic 9 (trace relation removes one); e₁ᵀE₁₂-type 5; e₁ᵀ(E₁₂+E₂₃)-type 6; e₁ᵀE₁₃-type 7.

### 5.1 Point B (k=3) and its undo

**Factors [exact].** Site 2 = torus-rotated vacuum with torus invariant y* = w₂(θ₁₂−h)/(θ₁₂−3h) at all three numeric points (−63/43 at P1), i.e. v₂ = Ω_{(w₁,y*)} = e₁ᵀρ(n) with ρ(n) = Rrot, checked for all six w₂ (cell 26). Site 1:

v₁ ∝ e₁ᵀ(α + βE₁₂) R(z,(w₁,y*)) ρ(Δ_{w₂}),  Δ = diag(1/(w₂−z_k)),

rotRank (site 1, site 2) = (5, 4) at all six w₂ (cell 25). The dressing e₁ᵀ(α+βE₁₂) = αe₁ᵀe^{(β/α)E₁₂} − (β²/2α)e₁ᵀE₁₂² is a rotation minus the extremal (1,3,0) covector (itself a rotation), i.e. a sum of two rotations (`L2-310-two-rotations.md`, Open). β/α(w₂) (cell 25, 32; **exact rationals, fit not derived; no closed form**):

| w₂ | 7/3 | −8/9 | 4/27 | 3 | −1/2 | 11/7 |
|---|---|---|---|---|---|---|
| y* | −63/43 | 24/43 | −4/43 | −81/43 | 27/86 | −297/301 |
| β/α | 392/2335 | 64/731 | −32/2409 | 8/77 | 1/23 | 4312/2693 |

Other numeric points: β/α = 294/2255 (h = 1/3), 5400/337687 (third point) (`L2-310-two-rotations.md`). c on (1, T₁₀, T₁₁, T₂₀, T₂₁, T₂₂, T₃₀, T₃₁, T₃₂) at P1: (81/1075, 36/215, 4/39, 1, 860/1053, 4300/9477, −860/729, −92450/85293, −1849000/2302911) (normalized c(T₂₀)=1).

**Undo [exact]** (cells 26–27): n = Rrot(z,(w₁,y*,w₃)) the Gauss element of v₂; v₁ = (αe₁ᵀ + βe₁ᵀE₁₂)ρ(n)ρdiag[Dy[w₂]] reproduces the pipeline's site-1 covector exactly; (ρ(n)⊗ρ(n))HW = HW (w₂ = 7/3, −1/2); ⟨v₁⊗v₂|HW⟩ = s⟨φ⊗e₁ᵀ|HW⟩ for all 10 HW kets (s = 1, the e₁-component of v₂). r_k at w₂ = 7/3: (4/5, −1071/605, 1401/484).

**Per-weight ranks over the six w₂ [exact]** (cells 28, 32; "(c0)" = control resolvent point → y* with β = 0; "(c)" = resolvent → y* with β). Weights in this order, dims in brackets:

| weight μ | dim | (a) point B | (b) β=0 | (c) ω=y*, β | (c0) ω=y*, β=0 | LW frame (a) | point A | A∪B |
|---|---|---|---|---|---|---|---|---|
| (0,1,3) | 1 | **0** | 1 | 0 | 0 | 1 | 1 | 1 |
| (0,2,2) | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 1 |
| (0,3,1) | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| (1,0,3) | 1 | **0** | 1 | 0 | 0 | 1 | 1 | 1 |
| (1,1,2) | 2 | **1** | 2 | 1 | 0 | 2 | 2 | 2 |
| (1,2,1) | 2 | 2 | 2 | 2 | 1 | 2 | 2 | 2 |
| (1,3,0) | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| (2,0,2) | 1 | 1 | 1 | 1 | 0 | 1 | 1 | 1 |
| (2,1,1) | 2 | 2 | 2 | 2 | 1 | **1** | 2 | 2 |
| (2,2,0) | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| (3,0,1) | 1 | 1 | 1 | 1 | 1 | **0** | 1 | 1 |
| (3,1,0) | 1 | 1 | 1 | 1 | 1 | **0** | 1 | 1 |
| total directions | 15 | 12 | 15 | 12 | 7 | 12 | 15 | 15 |

(Column totals computed by summing the per-weight ranks above; the notebooks' "total rank (15)" rows are ranks of a 6×15 matrix capped at 6 and are not used. The LW-frame row of this table is from cell 36, point A from cell 39, A∪B from cell 39/41.)

So the rotation part alone spans everything (the spec §3 prediction) **[exact]**, but the actual dressed covector spans only 12 of 15 directions: the fitted dressing cancels the weights (0,1,3), (1,0,3) identically and one direction of (1,1,2), for every w₂ (also 12 random w₂). **Spec §5's "the dressing only adds" is wrong** (corrected in spec §8). Consistent with site 1 lying in the e₁ᵀE₁₂ orbit (rotRank 5). Control (c0) reproduces the old blind spot: all components with μ₃ ≥ 2 vanish ((0,1,3),(0,2,2),(1,0,3),(1,1,2) both,(2,0,2)), plus one direction in each of (1,2,1), (2,1,1). Control with dressing (c) behaves as (a): the dressing alone rescues it at the HW count.

**HW count [exact]** (cells 29, 30, 31, 33).

| variant | fixed torus, six w₂ | fixed torus, 12 random w₂ | torus-completed (12 random (w₂, D)) |
|---|---|---|---|
| (a) point B | 6 (saturates at the 6 rows) | 7 | **10/10**, all irreps incl. [4,2,2], [3,3,2] |
| (b) β=0 | 6 | – | 10 |
| (c) ω=y*, β | 6 | – | 10 |
| (c0) ω=y*, β=0 | – | 5 | **6**: [3,3,2], [4,2,2] zero; [4,3,1] and [5,2,1] one direction short each |

The old families (2026-09-23) reached 8/10, missing [4,2,2] and [3,3,2] (site-1 weights (1,1,2), (0,2,2); `xxx-310-cg-test.tex` as quoted in `old-conventions/XXX-CG/xxx-210-cg-test.tex` §5: "8/10 and 9/10; missed exactly the irreps with ν₃ > λ₂"). Per-irrep determination from B (rank on the weight ν−Λ equal to its dimension): all nine determined except **[4,2,2]**, whose site-1 weight (1,1,2) has dim 2 but rank 1: only one linear combination of its two CG components is fixed by the family. The missing weights (0,1,3), (1,0,3) carry no HW state and matter only for a complete site-1 basis. At fixed torus with 12 random w₂ the rank is 7 (cross-irrep relations only, each irrep separately has full rank).

**Lowest-weight frame [exact]** (cells 34–36). e_lowᵀ is GT index 15 (weight (0,1,3)); v₂ = s e_lowᵀρ(m), m ∈ N₋ (m = (m₂₁, m₃₁, m₃₂), m₃₁ = 11/6, m₃₂ = 22/15 for all w₂; m₂₁ = 2335/1428, 425/324, 825/548, 1925/1164, 575/408, 67325/41916 at the six w₂; s = 925344/73205, 104976/1331, 177552/9317, 754272/62557, 528768/17303, 27161568/1987183); (ρ(m)⊗ρ(m))LW = LW, ⟨f₁⊗v₂|LW⟩ = s⟨φ′⊗e_lowᵀ|LW⟩, 10 LW kets with weights the reverses of the HW weights. Per-weight ranks of the LW-undone covector (table above): 1,1,1,1,2,2,1,1,1,1,0,0, the exact mirror of the HW frame (0,1,1,0,1,2,1,1,2,1,1,1) under μ → reversed. [4,2,2]'s LW slice sits at (2,1,1): rank 1 of 2. **[4,2,2] is not determined in either frame**; the frames measure different vectors, so they do not combine. Torus-completed LW count 10/10.

### 5.2 Point A (k=2)

**Factors.** Site-2 invariant y_A = w₂(θ₁₂−h)/(θ₁₂−2h) (the adjoint's y*; −63/8, 3, −1/2, −81/8, 27/16, −297/56 at the six w₂, found at every w₂), rotRank (site 1, site 2) = (6, 4) (cell 37). Site 2: v₂,A ∥ e₁ᵀρ(n_A) with n_A = Gauss element of rOm(z,(w₁,y_A,w₃)), r(7/3) = (4/5, −630/451, 1095/451); the exact `Solve` returns exactly one solution at each w₂ (cell 38). c_A = the adjoint's c* = (−2187/1600, −243/80, −27/43, 1, 320/1161) carried verbatim to the column-shifted shapes (1,1),(2,1),(2,2),(3,1),(3,2); on (1,T₁₀,T₁₁,T₂₀,T₂₁,T₂₂,T₃₀,T₃₁,T₃₂): (0,0,−2187/1600,0,−243/80,−27/43,0,1,320/1161) (cell 49).

**Undo [exact]** (cells 39–41): n_A the plain Gauss element at y_A; undo exact. Per-weight ranks 15/15 (table above, "point A"), nothing vanishes; A alone and A∪B determine all HW states; torus-completed HW rank 10 for A and for A∪B, per irrep all full. **Point A closes the [4,2,2] gap.** Verdict for [3,1,0]²: all 10 HW states reached by B; 9 fully determined by B; the tenth ([4,2,2]) needs A.

### 5.3 Dressing identifications

**B, first order.** The site-1 factor is the dressed vacuum with β/α the six fitted values above. The dressing order is 1 (a combination of {1, E₁₂}) and the E₁₃-degree is 0. Other orderings and the (w₂,y*)/D_{w₁}/(w₁,w₂) variants fail (`L2-310-two-rotations.md`, result 2, cell 13).

**A, second order [exact]** (cells 44–48, 61–63, 65). f₁,A is not in span{e₁ᵀ, e₁ᵀE₁₂, e₁ᵀE₂₃, e₁ᵀE₁₃}·ρ(R)ρ(D) for any of the 15 subsets and any D ∈ {ρ(Dy[w₂]), ρ(Dy[y_A]), 1} (rank 5 vs 4 at all six w₂). With D = ρ(Dy[w₂]) the undone covector g = f₁D⁻¹R⁻¹ occupies only the weights (3,1,0), (2,2,0), (2,1,1) (no E₂₃ component) and lies in the 4-dim span {e₁ᵀ, e₁ᵀE₁₂, e₁ᵀE₁₃, e₁ᵀE₁₂E₂₃}, in no 3-dim subset (dims of e₁ᵀ·words of length ≤ k: 1, 4, 9, 15). Exactly (rational equalities at 18 w₂ of P1, 10 one-parameter variations of (h, θ, z_k, w₁) and 28 random parameter points):

f₁,A ∝ e₁ᵀ (1 + βE₁₂ + γ(E₁₃ + ½E₁₂E₂₃)) ρ(n_A) ρ((w₂−Gᵀ)⁻¹),  v₂,A = e₁ᵀρ(n_A),

β = (z₁−z₂)(w₁−z₁) / (2z₂(w₁−z₂)) · w₂/(w₂−z₁),
γ = (θ₁₂−h)(z₁−z₃)(z₂−z₃)(w₁−z₁) / (4(θ₁₂−2h) z₂z₃(w₁−z₃)) · w₂² / ((w₂−z₁)(y_A−z₃)),  y_A = w₂(θ₁₂−h)/(θ₁₂−2h).

(β is independent of h, θ, z₃; at P1 K₁₂ = β(w₂−z₁)/w₂ = 2/25 and K₁₃ = γ(w₂−z₁)(y_A−z₃)/w₂² = 189/704, cell 61.) The dressing is a second-order element of U(n₊) with two parameters, E₁₂ and the weight-(2,1,1) vector E₁₃ + ½E₁₂E₂₃, plus the rotation and the resolvent torus; it is not an exponential e^{aE₁₂+bE₁₃} (that would give an E₁₂² component at (1,3,0)). So the two points are: B = first-order dressing ↔ pole set {0..3}; A = second-order dressing ↔ pole set {0,1,2}.

### 5.4 The det₂ identification of point A

**[exact/numeric]** (cells 49–53, 56, 60, 64, 65; "A is not on the det₂ curve" of spec §9 was true only for the pole set {0..3}; corrected by spec §13).

- **State-independence scan** (225 states): (i) n = 0..2 with Q^{j[−]} and minor W(2,3): 22.8 (fails); with M = W(3,4): **7·10⁻⁴⁸ (hit)**; (ii) n = 1..3: 259 (fails); (iii) n = 0..3 with Q^{j[+]}: 4.0 (fails), with xQ^j and M: 7·10⁻⁴⁸ (hit, same function as (i) up to a constant). So the polynomial f_A of point A is the ladder-1 det₂ over the pole set {0,1,2} (one-site weights on that set, insertion y_A^{−(x−θ₁)/h}, second column Q^{j[−]}), divided by the [3,1,0] minor M = W(3,4) (not W(2,3)). **[numeric]**
- **Symbolic y curve:** for P = {0,1,2} the curve lies in span{1, τ_λ(θ₁−3h)} to 10⁻⁴⁷, supported exactly on (T₁₁, T₂₁, T₂₂, T₃₁, T₃₂), with
  c(y) = (9800/5031 y⁴, −171500/5031 y³, 12005000/216333 y², −12005000/135837 y², 3361400000/17522973 y) at P1, and c(y) ∥ c_A only at y = −63/8 = y_A. **[exact]**
- **General form in (h, θ₁₂, y).** c_A = (2h²(θ₁₂−h)(θ₁₂−2h)(θ₁₂−3h) y⁴, 2h(θ₁₂−h)(θ₁₂−3h) y³, (θ₁₂−h) y², (θ₁₂−3h) y², (2/(3h)) y) on (T₁₁, T₂₁, T₂₂, T₃₁, T₃₂) (the adjoint's §5.4 monomials of `Two-rotations-det-overlap.md` on the shifted shapes), exact at all three numeric points (h = 1/2 P1; h = 1/3; θ = (1/5,−4/9), z′, w′). Symbolically f_A(y) = K·det₂^{(0..2)}[y]/M with K = −96h⁸(h−θ₁₂)(2h−θ₁₂)(3h−θ₁₂)²(4h−θ₁₂), identically in (h, θ₁₂, y), where h⁴ det₂ = (y⁴/4)W(2,3) − (y³/2)W(1,3) + (3y²/4)W(1,2) + (y²/4)W(0,3) − (y/2)W(0,2) (W(0,1) = 0). **[exact]**
- **End-to-end:** ⟨Φ_A|Ψ⟩/⟨x₀|Ψ⟩ ∝ det₆/Ψ_Q × det₂^{(0..2)}[y_A]/M to 5·10⁻⁴⁵ over 225 states (controls y = −63/43, 2 fail by 5.4 and 96). **[numeric]**

For B (cells 15–23 of `L2-310-two-rotations.wb`, Q-side done by an Opus subagent from `L2-310-Q-side-spec.md`, reviewed by the main session): with the rep-generic pole set n = 0..3, vacuum triple (θ−4h, θ−3h, θ−h), minor₁ = W(3,4), ⟨v₁⊗v₂|Ψ_B⟩/⟨x₀|Ψ_B⟩ ∝ det₆N[Q_B, Q_vac,(w₁,w₂)]/Ψ_Q(x₀) × det₂^{(0..3)}[Q_B; y*]/minor₁, state-independent to 5.1·10⁻⁴⁵ over 225 states, ratio₄/ratio₆ = 1/κ to 4.6·10⁻⁴⁷ (κ = 34306448400/147008443), wrong y's (−63/8, 2) fail by 77 and 322. The curve f_y ∈ span{1, τ_λ(θ₁−3h)} for every y (residuals ~10⁻⁴⁷; T₃₃(θ₁−3h) ≡ 0 so the basis has 9 functions; controls: random combination 10⁻⁵⁴, non-member 0.9), graded by box number:

c_λ(y) = a_λ y^{6−|λ|}, |λ| = λ₁+λ₂, a = (16/9, −2240/387, 78400/15093, 19600/387, −2744000/45279, 96040000/1946997, 2744000/31347, −48020000/407511, 6722800000/52568919) on (1,T₁₀,T₁₁,T₂₀,T₂₁,T₂₂,T₃₀,T₃₁,T₃₂) (the adjoint curve was y^{4−|λ|}; the pattern is y^{2λ₁−|λ|}). The only nonzero y with c(y) ∥ c_B is y = −63/43 = y* (exact `Solve`); no y gives c(y) ∥ c_A on {0..3}. **[exact/numeric]**

**End-to-end residuals summary** ([3,1,0]): B 5.1·10⁻⁴⁵; A 5·10⁻⁴⁵; all 225 states, 60 digits.

---

## 6. The explicit Q-dictionary (`L2-310-two-rotations.wb` cells 53–60, 65; spec §13(c))

Notation: ν_n(x) = Π_s (x−θ_s+hΛ_n) (so ν₃(θ₁−kh) = kh(kh−θ₁₂) for Λ₃ = 0), t₃ = χ₃ ν₁(x)ν₂(x−h)ν₃(x−2h), with Q's the dual Q's. The dual Baxter equation at x, A Q(x+h) − t₁B Q(x) + t₂C Q(x−h) − t₃Q(x−2h) = 0, A = ν₃(x)ν₃(x−h)ν₃(x−2h), B = ν₃(x−h)ν₃(x−2h), C = ν₃(x−2h), is solved by Cramer using any two dual Q's (written for the pair (2,3), W = W₂₃):

t₁(x) = [ ν₃(x) W(x+h, x−h) − t₃(x) W(x−2h, x−h) / (ν₃(x−h)ν₃(x−2h)) ] / W(x, x−h),
t₂(x) = [ ν₃(x)ν₃(x−h) W(x+h, x) − t₃(x) W(x−2h, x) / ν₃(x−2h) ] / W(x, x−h),
t₃ = χ₃ ν₁(x)ν₂(x−h)ν₃(x−2h)  (2·10⁻⁴⁶ on 225 states).

(Here W(a,b) with arguments x±… means the same bilinear in Q at those points.) Verified on all 225 states: ≤ 2·10⁻⁴⁶ (pair (2,3)), 1.7·10⁻⁴³ (pairs (1,3),(1,2)) at generic x = 3/11, θ₁−3h, θ₁−2h, θ₁−h **[numeric]**.

**Vanishing factors at the special points.**
- x = θ₁−3h: ν₁(x) = 0 (factor x−θ₁+3h), so t₃ = 0 and only the first terms survive; no other factor vanishes at θ₁−2h, θ₁−h.
- x = θ₁: ν₃(x) = 0, ν₂(x−h) = 0, so t₃ = 0, and W_ij(θ₁, θ₁−h) ≡ 0 for all pairs (Q(θ₁) = κQ(θ₁−h), checked 10⁻⁵⁰): Cramer degenerates, only t₂(θ₁) = κν₃(θ₁−h)t₁(θ₁) holds (10⁻⁵⁰); t₁, t₂ separately as limits (10⁻³²).
- x = θ₁+h: ν₃(x−h) = ν₃(θ₁) = 0 (t₁ drops out) and t₂(θ₁+h) = χ₃ν₁(θ₁+h)ν₂(θ₁) Q^j(θ₁−h)/Q^j(θ₁) (10⁻⁵⁰, j = 1,2,3); t₁(θ₁+h) as a limit only.

**Single-Wronskian form [exact]** (symbolic in h, θ₁₂, χ₃; also verified on all 225 states for the 8 nonzero shapes ≤ 3·10⁻⁴⁶ and against the symbolic CBR output at P1): by CBR with W(0,1) = 0 and the Plücker relations of the 2×n matrix (Q²;Q³),

τ_{(λ₁,λ₂)}(θ₁−3h) = C_λ · W(θ₁−(3−λ₁)h, θ₁−(4−λ₂)h) / W(θ₁−3h, θ₁−4h),
C_λ = Π_{k=4−λ₁}^{3} kh(kh−θ₁₂) · Π_{k=5−λ₂}^{4} kh(kh−θ₁₂)

(z- and χ₃-independent; kh(kh−θ₁₂) = ν₃(θ₁−kh)). Special values: (3,3) ↦ W(0,1) = 0, (0,0) ↦ 1. Examples: τ₁₁ = 12h²(3h−θ₁₂)(4h−θ₁₂) W(θ₁−2h, θ₁−3h)/M, τ₃₂ = 72h⁵(h−θ₁₂)(2h−θ₁₂)(3h−θ₁₂)²(4h−θ₁₂) W(θ₁, θ₁−2h)/M. Ladder points needed: θ₁−nh, n ∈ {0..4}; per shape (1,0), (1,1): {2,3,4}; (2,1): {1,3,4} (as single W(1,3)/M); the rest {0..4} before simplification, after it each is one W(a,b) and M.

**Consequence.** Every overlap of the form c·(1, τ_λ(θ₁−3h))·⟨Ω|Ψ⟩ ("T-polynomial × det₆") is explicit in Q's. For [3,1,0], both points are det₆ × det₂^{(P)}[y_P]/M with (P, y_P) = ({0..3}, w₂(θ₁₂−h)/(θ₁₂−3h)) and ({0,1,2}, w₂(θ₁₂−h)/(θ₁₂−2h)). For other reps the analogous single-Wronskian dictionary has **not** been derived (it was tested only indirectly by the det₂ ↔ τ expansions of §7, §8).

**Rep data for the Q-side (all with the poles n = 0..λ₁, insertions 1 and e^{2πix/h}).**

| rep | vacuum triple Ψ_Q = Π_s [Q²(θ_s−n₂₁h)Q³(θ_s−n₁h) − Q³(θ_s−n₂₁h)Q²(θ_s−n₁h)]·Q³(θ_s−n₂₂h) with (n₁, n₂₁, n₂₂) | minor | degrees | det₆ accuracy |
|---|---|---|---|---|
| [2,1,0] | (3,2,1) (Ψ_Q of `Two-rotations-det-overlap.md` §1: Q²(θ_s−2h)Q³(θ_s−3h)−…, Q³(θ_s−h)) | W(2,3) | – | 10⁻⁵⁰ |
| [3,1,0] | (4,3,1) | W(3,4) | deg q^j = μ_j, deg q_i = 6−μ_i; vacuum q degrees (0,4,6) | 5.1·10⁻⁴⁵ |
| [3,2,0] | (4,3,2) ([3,1,0]'s (4,3,1) fails by 22) | W(3,4) (W(2,3) fails by 0.60, W(2,4) 0.98, W(3,5) 0.98; best other (1,3) 0.59) | deg q^j = μ_j, deg q_i = 6−μ_i | 1.2·10⁻⁴⁴ |
| [4,1,0] | (5,4,1) | W(4,5) (next best (1,5) 0.82; W(3,4) 0.69, (3,5) 0.99, (4,6) 0.98) | deg q^j = μ_j, deg q_i = 8−μ_i; vacuum q degrees (0,6,8) (residual 2·10⁻⁴⁵) | 6.1·10⁻⁴³ |

Pattern: triple (θ−(λ₁+1)h, θ−λ₁h, θ−λ₂h), minor W(λ₁, λ₁+1) **[conjecture]** for further reps (verified on the four rows; the [3,1,0] triple was a guess from the GT-pattern dictionary x₁ = −(m₁₁+1)h, x₂₁ = −m₂₁h, x₂₂ = −m₂₂h). For [3,2,0] the other triples scanned fail by (4,2,1): 74, (3,2,1): 86, (4,2,0): 225, (3,2,0): 129, (4,1,0): 3175.

---

## 7. [3,2,0]² (`L2-320-two-rotations.wb`)

**Duality.** [3,2,0] ≅ [3,1,0]*⊗det³ (the weights map μ → (3−μ₃, 3−μ₂, 3−μ₁)). The roles of the two sites are swapped relative to [3,1,0] and the dressing is along E₂₃.

**Rep data.** dim 15, 225 states, 12 distinct weights (dim 2 at (1,2,2), (2,1,2), (2,2,1)); cell 2 reports 36 weight sectors for the two-site space. HW states: [6,4,0], [6,3,1], [6,2,2], [5,5,0], 2×[5,4,1], 2×[5,3,2], [4,4,2], [4,3,3] (10; cell 12). Site weights ν−Λ (Λ = (3,2,0)): (3,2,0), (3,1,1), (3,0,2), (2,3,0), (2,2,1), (2,1,2), (1,2,2), (1,1,3) (cell 14, 19).

**Localization at θ₁−3h [exact]** (cell 5, 21): site-1 localizers T₍₁,₁₎, T₍₂,₁₎, T₍₃,₁₎; T₍₂,₂₎, T₍₃,₂₎, T₍₃,₃₎ vanish on ⟨Ω|; T₍₁,₀₎, T₍₂,₀₎, T₍₃,₀₎ have rank 2 with Ω₁ in their site-2 span (matches `L2-320-localization.md`). Joint site-2 span 2-dim {Ω₁, ω₂} with ω₂ of rotRank 5 (not a rotation, not an Ω_{(p,q)}); the pencil Ω₁+κω₂ has rotRank 5 throughout.

**Single point [exact]** (cells 6–9; 105 quadrics, exact `Solve` 5.5 s, 55 redundant raw components, 53 localizing + 1 isolated nontrivial point). c on (1,T₁₀,T₁₁,T₂₀,T₂₁,T₃₀,T₃₁) = (81/1075, 36/215, 4/39, 1, 860/1053, −860/729, −92450/85293), the same numbers as [3,1,0] point B on the common shapes; c(T₁₀) = c(T₂₀) = c(T₃₀) = 0 on the localizing components. **No analogue of point A** (it lives on T₂₂, T₃₂, which vanish here). No continuous component. κ of the pencil = −50717940/541717.

**Factors [exact]** (cells 8–11, 20, 21).
- Site 1 = e₁ᵀ R(z,(w₁,y*)) ρ(Δ_{w₂}): a torus-rotated vacuum (rotRank 4), y* = w₂(θ₁₂−h)/(θ₁₂−3h) (the k = 3 member of w₂(θ₁₂−h)/(θ₁₂−kh); values for k = 1..4 at w₂ = 7/3: 7/3, −63/8, −63/43, −21/26), torus (2, 20/13, 1) ∝ Dy[w₂].
- Site 2 = e₁ᵀ(1 + bE₂₃)R(z,(w₁,y*)) with **no** torus factor (rotRank 5, e₁ᵀE₂₃ orbit type; E₁₂/E₁₃ dressings and D ∈ {Dy[w₁], Dy[w₂], Dy[w₃], Dy[y*]} all fail).
- **b fit (18 points: six w₂ plus twelve random; fit, not derived):** b = β/α = (6125/5676) w₂ / ((w₂−z₃)(y*−z₃)) at P1, w₁ = 5/2. Values at the six w₂: −343/484, −125/22, 125/308, −1125/2068, −375/286, −128625/113468. (A lead for the open closed form of the [3,1,0] β/α; not derived.)
- Torus invariants of site 1 at the six w₂ (y, d₁, d₂): (−63/43, 2, 20/13), (24/43, 28/43, 8/11), (−4/43, 28/73, 8/17), (−81/43, 14/9, 4/3), (27/86, 7/12, 2/3), (−297/301, 38/3, 76/27).

**Undo [exact]** (cells 12, 17): site 1 = s e₁ᵀρ(n′) with n′ = Δ⁻¹nΔ, so φ = f₂ρ(n′)⁻¹ = (e₁ᵀ+be₁ᵀE₂₃)ρ(nn′⁻¹) sits at site 2. Exact on all 10 HW kets, (ρ(n′)⊗ρ(n′))HW = HW; 3×3 identity as in §3 item 3.

**Per-weight ranks over six w₂ [exact]** (cells 13, 18, 20):

| weight μ | dim | (a) point | (b) β=0 | (c) ω=y*, β | (c0) ω=y*, β=0 | LW frame (a) |
|---|---|---|---|---|---|---|
| (0,2,3) | 1 | **0** | 1 | 0 | 0 | 1 |
| (0,3,2) | 1 | **0** | 1 | 0 | 0 | 1 |
| (1,1,3) | 1 | 1 | 1 | 1 | 0 | 1 |
| (1,2,2) | 2 | **1** | 2 | 1 | 0 | 2 |
| (1,3,1) | 1 | 1 | 1 | 1 | 0 | 1 |
| (2,0,3) | 1 | 1 | 1 | 1 | 1 | 1 |
| (2,1,2) | 2 | 2 | 2 | 2 | 1 | 2 |
| (2,2,1) | 2 | 2 | 2 | 2 | 1 | **1** |
| (2,3,0) | 1 | 1 | 1 | 1 | 1 | **0** |
| (3,0,2) | 1 | 1 | 1 | 1 | 1 | 1 |
| (3,1,1) | 1 | 1 | 1 | 1 | 1 | 1 |
| (3,2,0) | 1 | 1 | 1 | 1 | 1 | **0** |

So (a): 12 of 15 directions ((0,2,3), (0,3,2) ≡ 0; (1,2,2) rank 1 of 2); (b): 15/15; (c) as (a); (c0): the blind spot **μ₁ ≤ 1**, the dual of [3,1,0]'s μ₃ ≥ 2 under μ → (3−μ₃, 3−μ₂, 3−μ₁) (not "μ₃ = 3"); LW frame: mirror, (2,2,1) rank 1, (2,3,0), (3,2,0) ≡ 0. The LW frame does not repair the deficiency.

**HW counts [exact]** (cells 14, 16, 19, 21). Fixed torus six w₂: 6 of 6 rows ((a),(b),(c),(c0): 6,6,6,6); fixed torus with 12 random w₂: 7; **torus-completed (12 random (w₂, D)): 10/10 for (a), (b), (c); 6 for (c0)** (lost entirely: [4,3,3], [4,4,2]; one direction in each of the multiplicity-2 pairs [5,3,2], [5,4,1]). LW frame: 10/10 torus-completed. (The old site-1+site-2 localization needed ⟨Ω|Ψ⟩ for [4,3,3].) **Determination:** 9 of the 10 HW states (7 of 8 irreps) fully determined, all except **[4,4,2]**, whose site-2 weight (1,2,2) has dim 2 but rank 1; the LW frame does not repair it ((2,2,1) rank 1). There is only one family (no point A), hence 9/10, with the dual pattern of missing directions.

**Control blind spot** μ₁ ≤ 1: vanishing at (0,2,3), (0,3,2), (1,1,3), (1,2,2), (1,3,1) (cell 13, c0 column).

### 7.1 A2: the site-2 recipe is empty (cells 22–32; spec §12)

Paul's avenue A2: look for a dual of [3,1,0]'s point A on the site-2 side, using insertions at θ₂+mh. Result: **empty [exact/modular]**.
- Site-2 locals at θ₂+mh (m = −1, 0, 1, 2, −3; all nine shapes λ₁ ≤ 3): exactly T₍₁,₁₎(θ₂+h) = t₂(θ₂+h) and T₍₂,₂₎(θ₂); rank-2 operators T₍₂,₀₎(θ₂−h), T₍₁,₀₎(θ₂), T₍₃,₂₎(θ₂); Ω₁ lies in the site-1 span of every operator with rank ≤ 2; ⟨Ω|T = 0 for T₍₃,₃₎(θ₂−h) and T₍₂,₂₎, T₍₃,₂₎, T₍₃,₃₎ at θ₂−3h. Ranks (shapes (1,0),(1,1),(2,0),(2,1),(2,2),(3,0),(3,1),(3,2),(3,3)): m=−1: 3,3,2,6,6,6,11,12,0; m=0: 2,3,6,8,1,10,12,2,4; m=1: 3,1,6,3,4,10,6,8,8; m=2: 3,3,6,8,6,10,12,12,10; m=−3: 3,3,6,8,0,10,12,0,0.
- S₂ = span{1, T₍₂,₀₎(θ₂−h), T₍₁,₀₎(θ₂), T₍₂,₂₎(θ₂), T₍₃,₂₎(θ₂), T₍₁,₁₎(θ₂+h)}: joint site-1 span 3-dim (Ω₁, ω, ω′) at every w₂; ω a rotation (rotRank 4; the pencil Ω₁+κω has rotRank 4 at 3 values of κ, so {Ω₁, ω} lies in the rotation cone, unlike the θ₁−3h case where the extra direction had rotRank 5); ω′ has rotRank 7. The product locus of S₂ is only the localizing span (3 free coefficients: 1, T₍₂,₂₎(θ₂), T₍₁,₁₎(θ₂+h); the others 0), same at all six w₂. Enlargements: all nine shapes at θ₂+h (10 operators), at θ₂ (10), rank ≤ 4 operators at the five points (13): exact `Solve` with random 3×3 compressions, only the localizing span.
- Both points together (19 independent operators): exact `Solve` times out; Gröbner basis over GF(1000003) of the compressed minors (436 elements): every coefficient except 1, T₍₂,₂₎(θ₂), T₍₁,₁₎(θ₂+h) is nilpotent mod the ideal (order ≤ 4), affine dimension 3, so the locus is the localizing span (**[modular]** evidence, not a symbolic proof). The 22-operator span (adding rank ≤ 4 at m = −1, 2, −3) timed out, untested.
- Consequence: the predicted dual of [3,1,0]'s point A does not appear; [3,2,0] stays at 9/10 determined. Likely reason: only two site-2 localizers vs three at θ₁−3h.

### 7.2 The Q-side and the collapsed k = 2 curve (cells 33–41; spec §15)

- **Q-side port** (cells 34–35): 225 states, six Q's each, deg q^j = μ_j, deg q_i = 6−μ_i, Baxter residual 4.5·10⁻⁴³, eigen-residual 4·10⁻⁶⁰, 8 s for the whole Q-solve.
- **Vacuum triple** (θ−4h, θ−3h, θ−2h) (cell 36): 1.2·10⁻⁴⁴, i.e. Ψ_Q = Π_s W_s(3,4)·Q³(θ_s−2h); other triples fail (table in §6).
- **Minor** W(3,4) (cell 37): det₂^{(0..3)}[y]/W(a,b) lies in span{1, τ₁₀, τ₁₁, τ₂₀, τ₂₁, τ₃₀, τ₃₁} (τ₂₂ = τ₃₂ = τ₃₃ = 0 exactly on the states) only for W(3,4), residual 1.5·10⁻⁴⁵.
- **End-to-end** (cell 39): Φ_P = Σ_k c_{P,k}⟨Ω|O_k is rank 1 (exact); ⟨Φ_P|Ψ⟩/⟨x₀|Ψ⟩ ∝ det₆/Ψ_Q^{(4,3,2)} × det₂^{(0..3)}[y*]/W(3,4) to 1.2·10⁻⁴⁴ over 225 states with y* = −63/43 (controls: wrong y 77, 1692; wrong triples 22, 74, 86). **[numeric]**
- **Pole-subset scan** (cells 38, 40; all 11 subsets in span, residual ≤ 5·10⁻⁴⁶, rationalization ≤ 10⁻⁵³):

| P | support | y-deg | y hits c_P | generic rank Φ | nonzero y with Φ rank 1 |
|---|---|---|---|---|---|
| {0,1} | none (zero curve) | – | – | 0 | – |
| {0,2} | T₁₁, T₃₁ | 2..4 | none | 1 | all y (localizing span) |
| {0,3} | 1, T₂₁, T₃₀ | 3..6 | none | 2 | none |
| {1,2} | T₁₁, T₂₁ | 3..4 | none | 1 | all y (localizing span) |
| {1,3} | 1, T₁₁, T₂₀ | 4..6 | none | 2 | none |
| {2,3} | 1, T₁₀, T₁₁ | 4..6 | none | 2 | none |
| {0,1,2} | T₁₁, T₂₁, T₃₁ | 2..4 | none | 1 | all y (localizing span) |
| {0,1,3} | 1, T₁₁, T₂₀, T₂₁, T₃₀ | 3..6 | none | 2 | none |
| {0,2,3} | 1, T₁₀, T₁₁, T₂₁, T₃₀, T₃₁ | 2..6 | none | 2 | none |
| {1,2,3} | 1, T₁₀, T₁₁, T₂₀, T₂₁ | 3..6 | none | 2 | none |
| {0,1,2,3} | 1, T₁₀, T₁₁, T₂₀, T₂₁, T₃₀, T₃₁ | 2..6 | **−63/43** | 2 | only −63/43 (non-local, ∥ c_P) |

The {0,1,2} curve is (9800/5031 y⁴, −171500/5031 y³, −12005000/135837 y²) on (T₁₁, T₂₁, T₃₁): **the [3,1,0] point-A curve with the T₂₂, T₃₂ entries removed (they vanish on ⟨Ω|), hence entirely inside the localizing span**, a rank-1 curve with site-2 factor Ω₁ (checked), not a new point. So [3,2,0] has exactly one non-trivial product point (B at y*, pole set {0..3}), and the would-be point A collapses (the "k = 2 curve collapses"); this matches §10/§12 of the spec (no point A, 9/10 stays). Cells 8–13 of the notebook were not re-run for the A3 check.

---

## 8. [4,1,0]² (`L2-410-two-rotations.wb`, 26 cells)

**Caveat on provenance.** The notebook controller could not be associated for this new file, so the cells were evaluated in the shared kernel K1 by `ToExpression` of each cell's source (the sources are exactly what ran). **Outputs are not stored in the file**; cell 26 (markdown) holds all numbers, copied from those runs; rerun in VS Code to regenerate (about 1 min total). 53 Boolean checks in `checks`, all True. Timings: T_λ(θ₁−4h) for 14 shapes on 576 dims 11.6 s (exact sparse); exact product locus (276 quadrics, 14 unknowns) 4.5 s, 4 raw components; u21pipe per w₂ 4.7 s; eigenstates (576, 60 digits) 0.6 s, eigenvalue polynomials 8.9 s, all 3456 Q's 13 s, det₆ for all states 1.0 s, 21-pair minor scan 11 s, all 26 pole-subset curves 11 s.

**Rep data.** Λ = (4,1,0), dim 24, 576 states, 18 weights (12 of dim 1, 6 of dim 2). **14 HW states in 11 irreps** (not 15 as the spec §14 plan assumed): [8,2,0], [8,1,1], [7,3,0], 2×[7,2,1], [6,4,0], 2×[6,3,1], [6,2,2], [5,5,0], 2×[5,4,1], [5,3,2], [4,4,2] (dimensions add to 576).

**Localization at θ₁−4h [exact]** (cell 5): site-1 localizers T₍₁,₁₎, T₍₂,₂₎, T₍₃,₃₎, T₍₄,₃₎; T₍₄,₄₎ vanishes on ⟨Ω|; all other shapes rank 2 with Ω₁ in the site-2 span. Joint site-2 span 2-dim {Ω₁, ω₂}, ω₂ a rotation (rotRank 4; pencil Ω₁+κω₂ rotRank 4; not a vacuum Ω_{(p,q)}).

**Product locus [exact]** (cells 6–8): localizing span + **exactly three isolated points** (no continuous component), all with site-2 factor = torus-free vacuum rotation e₁ᵀρ(n_k), n_k from rOm(z,(w₁,y_k,w₃)): y_k = w₂(θ₁₂−h)/(θ₁₂−kh), k = 2, 3, 4: −63/8, −63/43, −21/26 at P1 (θ₁₂ = 31/35), torus d = (1,1); **k = 1 does not occur**. Site-1 rotRank 6, 7, 5 for k = 2, 3, 4 (site-2 rotRank 4). Supports of c: k=2 on T₂₂, T₃₂, T₃₃, T₄₂, T₄₃ (no identity); k=3 on the 9 shapes T₁₁, T₂₁, T₂₂, T₃₁, T₃₂, T₃₃, T₄₁, T₄₂, T₄₃; k=4 on all 14 live entries. The same three points at 6 + 18 further w₂ (stable over 24 w₂).

**Dressing of site 1 [exact]** (cells 8, 25): f₁ = g ρ(n_k) ρ(Δ_{w₂}) with g ∈ span of monomials e₁ᵀE₁₂^a E₁₃^c E₂₃^b; minimal, non-redundant sets (hold at all six w₂): k=4: {1, E₁₂, E₁₂²} (3 monomials); k=3: {1, E₁₃, E₁₂, E₁₂E₂₃, E₁₂E₁₃, E₁₂², E₁₂²E₂₃} (7, E₁₃-degree 1); k=2: those plus E₁₃², E₁₂E₁₃E₂₃ (9, E₁₃-degree 2). The E₁₃-degree is 0, 1, 2 for k = 4, 3, 2, i.e. λ₁−k; the dressing order 1 + E₁₃-degree = λ₁−k+1 = 1, 2, 3 holds only in this sense: the k = λ₁ point already carries E₁₂² (for [3,1,0] only E₁₂), and the minimal word length (words in E₁₂, E₁₃, E₂₃ applied to e₁) is 3, 3, 2 for k = 2, 3, 4, not 3, 2, 1.

**Q-side [numeric]** (cells 18–23): all 576 states, six Q's each; deg q^j = μ_j, deg q_i = 8−μ_i exactly; Baxter residual 3.4·10⁻⁴⁰, eigen-residual 1.4·10⁻⁵⁸; vacuum q degrees (0,6,8) (2·10⁻⁴⁵). det₆ with 10 poles n = 0..4: ⟨Ω|Ψ⟩/⟨x₀|Ψ⟩ ∝ det₆/Ψ_Q for the triple (5,4,1), i.e. (θ−5h, θ−4h, θ−h), 6.1·10⁻⁴³; the other candidates fail: (5,4,2) 850, (5,4,3) 3730, (5,3,1) 1.7, (4,3,1) 4.1; no other triple of the 126 scanned (1 ≤ n₁ ≤ 6) works. Minor: only W(4,5) puts det₂^{(0..4)}[y]/W(a,b) in span{1, τ_λ(θ₁−4h)} (14 functions, τ₄₄ = 0): residual 3.4·10⁻⁴⁴; W(3,4) 0.69, (3,5) 0.99, (4,6) 0.98, best other pair (1,5) 0.82.

**Pole-set curves [exact/numeric]** (cell 22): all 26 subsets P ⊆ {0..4}, |P| ≥ 2, with minor W(4,5): every curve lies in the span (residual ≤ 3.3·10⁻⁴³) and rationalizes (≤ 3·10⁻⁵³). Exactly three subsets hit a product point: {0,1,2} at y = −63/8 = y₂ (c₂), {0..3} at y = −63/43 = y₃ (c₃), {0..4} at y = −21/26 = y₄ (c₄); each is the only nonzero y with Φ(c(y)) rank 1 and non-local site-2 factor. The only other rank-1 curve is P = {0,1} (support T₃₃, T₄₃, site-2 factor Ω₁: inside the localizing span, trivial). All other subsets: generic rank 2 and no rank-1 root.

**End-to-end [numeric]** (cell 23): ⟨Φ_k|Ψ⟩/⟨x₀|Ψ⟩ ∝ det₆/Ψ_Q^{(5,4,1)} × det₂^{(0..k)}[y_k]/W(4,5) over all 576 states, 6.1·10⁻⁴³ for k = 2, 3, 4 (Φ_k rank 1). Controls: wrong y (−63/7, 2): 146/380 (k=2), 888/6062 (k=3), 58196/961 (k=4); wrong pole set {0..j} at y_k: 9.8 to 10⁵; triples (5,4,2), (4,3,1): 850 etc.

**Undo and ranks [exact]** (cells 12–16). Undo exact: (ρ(n_k)⊗ρ(n_k))HW = HW and ⟨f₁⊗f₂|HW⟩ = s⟨φ⊗e₁ᵀ|HW⟩ on all 14 HW kets, all k, all six w₂. Per-weight ranks of φ over six w₂:

| k | pole set | y_k | dressing S_k | directions spanned | weights short | irreps fully determined | HW reached (torus-completed, 18 random (w₂, D)) |
|---|---|---|---|---|---|---|---|
| 2 | {0..2} | w₂(θ₁₂−h)/(θ₁₂−2h) | 9 monomials (E₁₃-degree 2) | **24/24** (all 18 weights full) | none | **11/11** | 14/14 |
| 3 | {0..3} | w₂(θ₁₂−h)/(θ₁₂−3h) | 7 (E₁₃-degree 1) | 21/24 | (0,1,4), (1,0,4) rank 0 of 1; (1,1,3) rank 1 of 2 | **11/11** | 14/14 |
| 4 | {0..4} | w₂(θ₁₂−h)/(θ₁₂−4h) | 3 (E₁₂² only extra) | 16/24 | (0,1,4), (0,2,3), (1,0,4), (2,0,3) rank 0 of 1; (1,1,3) rank 0 of 2; (1,2,2), (2,1,2) rank 1 of 2 | 9/11: [5,3,2] (site weight (1,2,2)) and [6,2,2] ((2,1,2)) each one direction short | 14/14 |

Unions {2,3} and {2,3,4}: all 18 weights full, 11/11 determined. Fixed torus (D = 1) with 6 w₂ saturates at 6 rows; with 18 random w₂ the HW rank is 6, 7, 8 (k = 2, 3, 4), i.e. only cross-irrep relations; each irrep separately has full rank. Each single k reaches all 14 HW states once the torus is included (18 rows needed since 14 HW kets need at least 14 rows).

**Verdict (cell 26).** The chain structure is confirmed. Corrections to the spec §14 conjecture: (i) dressing order λ₁−k+1 holds as 1 + E₁₃-degree, but E₁₂² already appears at k = λ₁ and the span sizes are 3, 7, 9; (ii) "k < λ₁ spans every weight space" is true only for k = 2; k = 3 is 3 directions short of full (but determines all 11 irreps), k = λ₁ = 4 is 8 directions short and misses one direction in each of [5,3,2] and [6,2,2] (two directions short for the HW determination, not one); (iii) the HW count is 14, not 15.

---

## 9. The pole-set chain

**Statement [conjecture; verified instances below].** For Λ = [λ₁,1,0] the product points of the recipe

Φ(c) = ⟨Ω⊗Ω|(c₀ + Σ_λ c_λ T_λ(θ₁−λ₁h))  (a product covector)

(other than the localizing span) are in bijection with the initial segments P_k = {0..k}, k = 2..λ₁ (k = 1, y = w₂, the ω = y degeneracy, is excluded). The point k sits at

y_k = w₂ (θ₁₂−h)/(θ₁₂−kh),  with ω = w₂ = y_k (θ₁₂−kh)/(θ₁₂−h),

has site-2 factor the torus-free vacuum rotation e₁ᵀρ(n_k) (n_k from rOm at (w₁,y_k,w₃)), site-1 factor g ρ(n_k)ρ(Δ_{w₂}) with g in the span of monomials of E₁₃-degree ≤ λ₁−k, and overlap

⟨v₁⊗v₂|Ψ⟩/⟨x₀|Ψ⟩ ∝ det₆[Q_B, Q_vac,(w₁,w₂)]/Ψ_Q × det₂^{(P_k)}[Q_B; y_k] / W(λ₁, λ₁+1),

with the rep-data triple (θ−(λ₁+1)h, θ−λ₁h, θ−λ₂h) in Ψ_Q. The spec §14 version said "dressing order λ₁−k+1", corrected by §16 to "E₁₃-degree λ₁−k" (and the claim "only k < λ₁ points give complete spans" to "only k = 2 does").

**Verified instances.**

| rep | k's present | where verified |
|---|---|---|
| [2,1,0] | k = 2 (pure rotations) | `Two-rotations-det-overlap.md` §5; this note §4 |
| [3,1,0] | k = 2 (A), k = 3 (B); all 11 pole subsets P ⊆ {0..3} scanned: only {0..3}→B (y = −63/43) and {0,1,2}→A (y = −63/8); {0,1} trivially localizing (support T₂₂, T₃₂); all others never produce a product covector; max in-span residual 5·10⁻⁴⁶, rationalization 1.6·10⁻⁵³ | `L2-310-two-rotations.wb` cells 66–70 |
| [4,1,0] | k = 2, 3, 4 at −63/8, −63/43, −21/26; all 26 subsets scanned | `L2-410-two-rotations.wb` cell 22, 26 |
| [3,2,0] | only k = 3; the k = 2 curve collapses into the localizing span | `L2-320-two-rotations.wb` cells 38, 40, 41 |

**The scan table for [3,1,0]** (cell 70; support = shapes with nonzero coefficient, y-degrees = range of powers of y, rank = generic rank of Φ(c(y)); "nonzero y with Φ rank 1" is the exact common-root set of all 2×2 minors of Σ_k c_k(y)A_k; y = 0 is a root for every P but c(0) = 0, so Φ = 0 there):

| P | support | y-degrees | y hits B | y hits A | generic rank | nonzero y with Φ rank 1 |
|---|---|---|---|---|---|---|
| {0,1} | T₂₂, T₃₂ | 1..2 | none | none | 1 (localizing) | all y (trivial: site 2 = Ω₁) |
| {0,2} | T₁₁, T₂₂, T₃₁ | 2..4 | none | none | 2 | none |
| {0,3} | 1, T₂₁, T₃₀ | 3..6 | none | none | 2 | none |
| {1,2} | T₁₁, T₂₁, T₂₂ | 2..4 | none | none | 2 | none |
| {1,3} | 1, T₁₁, T₂₀, T₂₂ | 2..6 | none | none | 2 | none |
| {2,3} | 1, T₁₀, T₁₁ | 4..6 | none | none | 2 | none |
| {0,1,2} | T₁₁, T₂₁, T₂₂, T₃₁, T₃₂ | 1..4 | none | **−63/8** | 2 | only −63/8 (= A) |
| {0,1,3} | 1, T₁₁, T₂₀, T₂₁, T₂₂, T₃₀, T₃₂ | 1..6 | none | none | 2 | none |
| {0,2,3} | 1, T₁₀, T₁₁, T₂₁, T₂₂, T₃₀, T₃₁ | 2..6 | none | none | 2 | none |
| {1,2,3} | 1, T₁₀, T₁₁, T₂₀, T₂₁, T₂₂ | 2..6 | none | none | 2 | none |
| {0,1,2,3} | all nine | 1..6 | **−63/43** | none | 2 | only −63/43 (= B) |

**Three empirical rules and the derived criterion.** (Rules 1 and 3 are statements read off the four representations; rule 2's closed form is my compilation of the localizer lists and is tagged accordingly.) Let d := λ₁−λ₂.

1. **Vanishing rule [conjecture; three reps + adjoint-consistent].** ⟨Ω|T_{(a,b)}(θ₁−λ₁h) = 0 iff b > d. Evidence: [3,1,0] (d=2): T₃₃ vanishes, T₂₂, T₃₂ do not; [3,2,0] (d=1): T₂₂, T₃₂, T₃₃ vanish, shapes with b ≤ 1 survive; [4,1,0] (d=3): T₄₄ vanishes. For the adjoint (d=1) the notebooks only used shapes with b ≤ 1.
2. **Localizers [conjecture; derived from the four lists].** The site-1-local shapes at θ₁−λ₁h are {(a, min(a,d)) : a = 1..λ₁}: adjoint (1,1),(2,1); [3,1,0] (1,1),(2,2),(3,2); [3,2,0] (1,1),(2,1),(3,1); [4,1,0] (1,1),(2,2),(3,3),(4,3). (All four lists are as stated in the notebooks; the closed pattern a ↦ (a, min(a,d)) is an observation made while compiling these notes, not a claim in the sources.)
3. **Support rule [numeric, three instances].** The coefficient curve of point k on [λ₁,1,0] is supported on the support of the adjoint's curve shifted by s = λ₁−k columns, (a,b) ↦ (a+s, b+s): [3,1,0] k=2: (1,T₁₀,T₁₁,T₂₀,T₂₁) ↦ (T₁₁,T₂₁,T₂₂,T₃₁,T₃₂) ✓ (cell 49, 64); [4,1,0] k=2: s=2, adjoint shapes (0,0),(1,0),(1,1),(2,0),(2,1) ↦ T₂₂, T₃₂, T₃₃, T₄₂, T₄₃ ✓; [4,1,0] k=3 (s=1) reproduces the 9 shapes of [3,1,0] B shifted: T₁₁, T₂₁, T₂₂, T₃₁, T₃₂, T₃₃, T₄₁, T₄₂, T₄₃ ✓; k=λ₁ is all live entries (B). On [3,2,0] the k=3 support equals [3,1,0] B minus the vanishing T₂₂, T₃₂, with the same numerical c on the common shapes (cell 21).

**Derived criterion [conjecture beyond the four reps].** A non-trivial point k exists iff its shifted support contains a shape that survives (rule 1) and is not a localizer (rule 2). The entry (1,0) of the adjoint curve is shifted to (1+s, s), with s = λ₁−k; it survives iff s ≤ d and is non-local iff s < d (for s = d, (1+s,s) = (a, min(a,d)) is a localizer; for s > d it vanishes). Hence the point exists iff s < d iff **λ₂ < k ≤ λ₁** (k ≥ 2 always). Check: [3,1,0] k = 2, 3 > 1 ✓; [4,1,0] k = 2, 3, 4 > 1 ✓; [3,2,0] k = 2 is not > λ₂ = 2 (collapses, ✓), k = 3 > 2 ✓.

**The "[3,2,0]-type" class.** d = 1, i.e. λ₂ = λ₁−1: only k = λ₁ survives, so only B-type points exist, with a dressed (non-rotation) factor for λ₁ ≥ 3 and 9/10-type determination deficits. The adjoint (λ₁ = 2, λ₂ = 1) is excepted: there k = λ₁ = 2 is the only point too, but E₁₂² annihilates e₁ᵀ, so the dressing e₁ᵀ(α+βE₁₂) is itself a rotation and the family is pure rotations with full span. (Observed on [2,1,0], [3,1,0] B, [4,1,0] k=4: the dressing at k = λ₁ is the truncated exponential {1, E₁₂, …, E₁₂^{λ₁−2}}: {1} for λ₁=2 (pure rotation), {1,E₁₂} for 3, {1,E₁₂,E₁₂²} for 4; **[conjecture]**, three data points, an observation made while compiling, not stated in the sources.) For [3,2,0] B the dressing is {1, E₂₃} on site 2 (sites swapped).

**Predictions [conjecture, not verified].**
- [4,2,0] (d = 2, dim 27): non-trivial points exist for k = 3, 4 (λ₂ = 2 < k ≤ 4); k = 2 collapses into the localizing span. The k = 4 point is B-type; k = 3 carries one more E₁₃ order. The per-weight spans and whether [4,2,0] has all HW states determined from these points are **open**. `L2-420-localization.md` and `Blindspot-420-check-EH.wb` hold the localization data for [4,2,0] (not read for this note).
- [λ,1,0]: points k = 2..λ, the k = 2 point spans every weight space and determines all HW states (verified λ = 3, 4; adjoint trivially), k = λ loses the most (verified λ = 3, 4).
- General [λ₁,λ₂,0] with d ≥ 2: points k = λ₂+1..λ₁.

(Not verified: any rep with λ₂ ≥ 2 other than [3,2,0]; any rep with three distinct rows (the group here is gl(3), so all weights have three entries, λ₃ = 0 normalization); L > 2.)

---

## 10. What "all CG coefficients" would still need

Three layers separate the present results from CG numbers.

**Layer 1, kinematic completeness (HW ket coverage).** The family reaches/determines the HW kets as in the summary table. This layer is **done** for the adjoint, [3,1,0] (via A), [4,1,0] (via k = 2) and **incomplete** for [3,2,0] ([4,4,2] one direction short) and unattempted for [4,2,0] and beyond. The old 2026-09-23 count (6 of 8 directions, 5 of 6 HW states; the rule μ₃ ≤ λ₂) is superseded by "undo the rotation, count per site-1 weight, reach everything for the k = 2 point".

**Layer 2, explicit Q-form.** The overlap of each chain point is det₆ × det₂^{(P_k)}[y_k]/W(λ₁,λ₁+1) **[numeric, Q-checked]**, and for [3,1,0] each fused T-eigenvalue is a single Wronskian over a minor (§6). What remains: the per-weight projections of φ (the "coefficients of the free parameters") expressed through Q's, that is, each y-coefficient (and each torus D-monomial) of det₆·det₂^{(P_k)}[y] as one HW CG functional. The note `L2-310-two-rotations.md` (Open) and HANDOFF NEXT item (4) say this is the next step: put the k = 2 point through the Clean-CG pipeline for [3,1,0].

**Layer 3, untwisting/contamination.** The overlaps are with the finite-twist eigenstate Ψ_B. The CG coefficients are the HW components of the untwisted state, but (R⊗R)Ψ has components outside the HW subspace ("contamination" K, `Clean-CG-combinations.md` §2): ⟨Ω|T_λΨ⟩/f_λ = c_λ·P_HWΨ (signal S) + c_λ·(R⊗R)(1−P_HW)Ψ (contamination K), κ = (S+K)/S, with κ = (t+1)/2 = 97/70 for [3,3,0] and [4,1,1] etc. (t = θ₁₂/h there). The two statements that matter here: (i) the undo of §2 is exact for HW kets, but the eigenstate Ψ_B at finite twist is not a HW ket, so the leading-order untwisting needs a cleaning step; (ii) in `Clean-CG-combinations` the cleaning is done by top-killed (X^B − X^top)⟨Ω|Ψ_B⟩ and structurally-cleaned insertion combinations (all ten reps [λ₁,λ₂,0], λ₁ ≤ 5, λ₁>λ₂>0, fitted coefficients, closed forms open; adjoint derived analytically in `Adjoint-clean-analytic.wb`, β = h(2t−3)/(t−1) with t = θ₁₂/h). This layer has **not been run with the det₆ × det₂ family**: it would require the untwisting limit of the chain points (does the y_k and the torus survive as free CG-generating parameters? the old note `xxx-210-cg-test.tex` Result 3 says k₁, k₂ survived as generating-function variables, θ₁ and h dropped out) and the cleaning of the E₁₃-direction contributions. Here the leading untwisting orders per weight (adjoint ε⁰,…,ε⁴; §3) are a secondary check only.

**Pointers.** `Clean-CG-combinations.{md,tex,pdf}`, `Adjoint-clean-analytic.wb`, `Minimal-ingredients-{EH,310-EH,320-EH,410-EH,420-EH,...}.wb`, `Clean-coeffs-{210,310,410}-EH.wb` (the cleaned-combination notebooks); `Two-rotations-det-overlap.md` §6 (open questions); the 2026-09-19 twist-reconstruction thread (a single per-site E₃₁ rotation completes adj⊗adj), cited from `xxx-210-cg-test.tex`, Remark after Result "a hypothetical fix".

---

## 11. Open problems and suggested next steps

1. **[3,2,0]-type reps (λ₂ ≥ 2, d = 1 in particular).** The θ₁−λ₁h recipe gives only point B; the site-2 recipe is empty (§7.1). Find the missing covector(s) for [4,4,2] (site weight (1,2,2), one direction): the dual/LW construction was tried (LW frame does not repair it); products of insertions at both sites (spec §11 A3: T_λ(θ₁−3h)·T_μ(θ₂+mh), shapes at θ₁−h, θ₁, θ₁+h) are untested; "every product point is rotation ⊗ n₊-dressed rotation; generic dressing ⇒ complete" is the untested conjecture of A3.
2. **[4,2,0]** (dim 27): which k survive (prediction: k = 3, 4; §9), per-weight spans and HW determination.
3. **Dressing polynomials.** Closed forms of β/α for [3,1,0] point B (six fitted values; the [3,2,0] fit b = (6125/5676) w₂/((w₂−z₃)(y*−z₃)) is a lead), the minimal monomial sets and their coefficients for [4,1,0] (k = 2: 9, k = 3: 7, k = 4: 3 monomials), and the general c_k(y) curve from the det₂ residue expansion with the τ-dictionary of §6 (the τ ↔ Wronskian dictionary is known only for [3,1,0]; for [3,2,0], [4,1,0] the curves were fitted to span{1, τ_λ}).
4. **Q-side reading of CG coefficients** (Layer 2/3): put the k = 2 point of [3,1,0] through the Clean-CG pipeline; the untwisting of the family.
5. **[4,1,0] reproducibility:** rerun `L2-410-two-rotations.wb` in VS Code so the outputs are stored (cell 26 holds the numbers).
6. **Remaining numeric points:** the Q-side checks at the other two numeric points (h = 1/3; θ = (1/5,−4/9)) for [3,1,0] (`L2-310-two-rotations.md`, Open); the determinant identities of the chain at P1 only for [3,2,0], [4,1,0].
7. **The spec §12 22-operator span for A2** timed out (untested); the 19-operator evidence is modular only.

---

## 12. Appendix

### A.1 Notebook/cell map

All notebooks in `Paul/Mathematica/new-conventions/`. Cells are 1-based.

| result | notebook | cells |
|---|---|---|
| adjoint setup, undo exactness, HW kets, per-weight ranks (numeric/exact in y), closed forms (3×3), (0,1,2) component, untwisting, HW count, torus completion, summary | `L2-rotation-sov.wb` | 113 (intro), 114 (12.1), 115 (12.2), 116 (12.3–12.5), 117 (12.6), 118 (12.6b), 119 (12.7), 120 (12.8, 12.9), 121 (12.9b), 122 (result Grid), 123 (markdown summary) |
| family V, torus = character, general rotation | `L2-rotation-sov.wb` | 106 (y free), 108 (torus), 109 (fresh w₁), 110 (no iteration), 111 (generic g), 112 |
| [3,1,0] point B undo, ranks, HW count, 3×3, result | `L2-310-two-rotations.wb` | 24 (intro), 25 (six w₂, β/α), 26 (undo), 27 (HW kets), 28 (per-weight ranks), 29 (HW count), 30 (extra controls), 31 (3×3), 32 (result Grid), 33 (markdown result) |
| [3,1,0] LW frame, point A undo | `L2-310-two-rotations.wb` | 34 (m, s), 35 (LW kets), 36 (LW ranks), 37 (point A records), 38 (v₂,A), 39 (A undo ranks), 40 (A torus HW), 41 (result Grid), 42 (markdown) |
| [3,1,0] A1: dressing of A, det₂ identification, Q-dictionary | `L2-310-two-rotations.wb` | 43 (intro), 44–48 (records, dressing spans), 49–53 (polynomial, minors, curves, general form, Cramer), 54–60 (degeneracies, τ-dictionary, end-to-end 56), 61–63 (β, γ closed forms), 64 (result Grid), 65 (markdown result) |
| [3,1,0] pole-subset scan | `L2-310-two-rotations.wb` | 66–69, 70 (result) |
| [3,1,0] Q-side for B, product points | `L2-310-two-rotations.wb` | 5–13 (localization, product locus, factors), 15–23 (Q-side, vacuum triple, ladder-1 det₂ curve, end-to-end 21); `L2-310-two-rotations.md` |
| [3,2,0] localization, product point, factors, undo, ranks, HW, LW | `L2-320-two-rotations.wb` | 1–21 (cell 21 summary); 5 localization, 6–9 product locus, 11 six w₂, 12 undo, 13 ranks, 14 HW fixed, 15–16 torus, 17 3×3 and b fit, 18–19 LW and per-irrep, 20 result Grid |
| [3,2,0] A2 site-2 recipe | `L2-320-two-rotations.wb` | 22–32 (cell 32) |
| [3,2,0] Q-side, triple, minor, scan | `L2-320-two-rotations.wb` | 33–41 (cell 41); 34–35 Q-side, 36 triple, 37 minor, 38 curves, 39 end-to-end, 40 scan |
| [4,1,0] all | `L2-410-two-rotations.wb` | 1–26: 5 localization, 6–8 locus and factors, 9–16 pipeline, undo, ranks, HW, torus; 17–23 Q-side, triple, minors, pole-set scan (22), end-to-end (23); 24 chain Grid; 25 dressing; 26 (markdown, all numbers; outputs not stored) |

### A.2 Checks summary

- Adjoint `L2-rotation-sov.wb` §12: all boolean `12.x` checks True (cell 122, last row). 
- [3,1,0]: undo exact on 10 HW kets for all six w₂ (cells 26–27, 35); 3×3 identity symbolic including symbolic z; end-to-end 5.1·10⁻⁴⁵ (B), 5·10⁻⁴⁵ (A); Cramer t₁, t₂ ≤ 2·10⁻⁴⁶; τ ↔ Wronskian ≤ 3·10⁻⁴⁶; 11 subsets in span 5·10⁻⁴⁶.
- [3,2,0]: 56 Boolean checks pass (cell 21); Baxter 4.5·10⁻⁴³; end-to-end 1.2·10⁻⁴⁴; 11 subsets in span ≤ 5·10⁻⁴⁶; A2: Gröbner modular only.
- [4,1,0]: 53 Boolean checks pass (cell 26); end-to-end 6.1·10⁻⁴³; 26 subsets in span ≤ 3.3·10⁻⁴³, rationalization ≤ 3·10⁻⁵³.
- My own sympy check (while writing): formulas for a′, b′, c′, the lower entries of n′ vanishing, b′/(a′c′) = (y−ω)/(y−z₂), and b′|_{ω=y} = 0 reproduced symbolically with (z₁,z₂,z₃,w₁,y,ω) symbolic. (Not a notebook cell.)

### A.3 Tooling notes

- **Shared kernel K1.** Private kernels (K310, K2) cannot start on this machine; every notebook ends up on the shared default kernel K1, and the first cell of each notebook does a `ClearAll`, which wipes the other notebooks' live state (files intact). At the end of 2026-10-09: K1 held `L2-320-two-rotations.wb` (all cells live) per HANDOFF; `L2-rotation-sov.wb` and `L2-310-two-rotations.wb` live state wiped (files saved: 123 and 70 cells; HANDOFF's "42" predates the A1/A3 cells). Rebuild for `L2-rotation-sov.wb` §12: cells 2–21 plus the text of `oneSiteVac`/`w2Of`; cells 63–111 need their own cells.
- **Controller association.** For `L2-320-two-rotations.wb` and `L2-410-two-rotations.wb` the agent had to associate the notebook controller through the wolfbook stdio bridge (320: succeeded, K1 restarted once; 410: failed, cells evaluated by `ToExpression` of the saved sources, **outputs not stored**).
- **The `y*)` comment gotcha.** The string `y*)` inside a `(* … *)` comment closes the comment; write `y★` or `yS` in comments. It bit cell 12 of the [3,1,0] notebook twice on 2026-10-09 (memory note exists).
- **Memoized definitions.** After fixing a definition re-run everything downstream; the pipeline cells use fresh definitions (no shared memo) at explicit (h, θ, z, w).
- Output cells of the 410 notebook are empty; the 310 and 320 result Grid cells (32, 41, 64, 70, 21, 20) and the markdown result cells (33, 42, 65, 70, 21, 32, 41, 26) are the stable summaries.

### A.4 Source inconsistencies found while compiling, and how they were resolved

1. **Spec §3(iii) sign** of the (0,1,2) component (−b′(b′+a′c′) predicted vs 3b′(a′c′−b′) found): the notebook (cell 118, exact) wins; spec §6 and cell 123 record the correction.
2. **Spec §5 "the dressing only adds"**: false for [3,1,0] point B (12/15, not 15/15); corrected by spec §8 / cell 33. Same for [3,2,0].
3. **Spec §7 "10 distinct weights"** for [3,1,0]: 12 (spec §8, cell 33).
4. **Spec §9 "A is not on the det₂ curve"**: true only for the pole set {0..3}; A is a det₂ over {0,1,2} with minor W(3,4) (spec §13).
5. **Spec §14 conjecture** (dressing order λ₁−k+1; only k < λ₁ complete; 15 HW states for [4,1,0]): corrected by §16 to E₁₃-degree λ₁−k, only k = 2 complete, 14 HW states.
6. **HANDOFF** says the [3,1,0] notebook has 42 cells ("files saved: 123 and 42 cells"), while the spec cell references run to cell 70 and the file has 70 cells (A1 and A3(i) were added later the same day); the HANDOFF figure refers to an earlier moment.
7. **t = θ₁₂ versus t = θ₁₂/h.** The adjoint notes (`Two-rotations-det-overlap.md` §1, `Clean-CG-combinations.md`) write t := θ₁₂/h = 62/35; the [3,1,0]/[4,1,0] notebooks, spec §13 and the formula for K, C_λ use t = θ₁₂ = 31/35 (consistent with y* = w₂(θ₁₂−h)/(θ₁₂−3h) = −63/43). This note writes θ₁₂ explicitly wherever it matters.
8. **ω rotRank.** The second site-2 direction ω₂ at θ₁−λ₁h has rotRank 5 for [3,2,0] (cell 6, 21) and rotRank 4 for [4,1,0] (cell 26) and is a rotation for the adjoint ([Ω_{(w₁,z₂)}], `Two-rotations-det-overlap.md` §5.2); for [3,1,0] "a rotation, pencil inside the orbit" (`L2-310-two-rotations.md`). Spec §12 reports a different ω (at the site-2 points θ₂+mh) with rotRank 4 and ω′ with rotRank 7. These are different directions at different insertion points, not a contradiction; HANDOFF's "ω₂ of rotRank 5" refers to [3,2,0] at θ₁−3h.
9. **[3,2,0] rotation placement.** Spec §10 says "site 1 = pure rotation" and the undo is "of the site-1 rotation" (so φ lands on site 2); the [3,1,0] convention has the rotation on site 2. Consistent via the flip symmetry of V⊗V; in the tables the site carrying the CG vector is site 2 for [3,2,0] and site 1 for the others.
10. **"Total rank (15)" rows** in the 310/320 per-weight Grids are ranks of a 6×15 matrix (capped at 6); not used.
11. **Not found in the sources** (not guessed): the vacuum q-degrees for [3,2,0]; the outcome of the A2 22-operator span (timed out); β/α(w₂) at the twelve extra random w₂ of the [3,1,0] check (stated only as "also 12 random w₂"); the Q-side checks of [3,1,0] at the two other numeric points; the T-eigenvalue dictionary (single Wronskians, C_λ) for reps other than [3,1,0]; stored outputs of `L2-410-two-rotations.wb`.
