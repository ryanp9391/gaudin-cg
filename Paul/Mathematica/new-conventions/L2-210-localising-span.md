# gl(3) adj⊗adj: what the transfer-matrix algebra reaches on a single SoV vacuum — reachability map of Q-number product covectors (2026-10-10)

Notebook: `L2-210-localising-span.wb` (63 cells, saved; **no native outputs** — VS Code did not execute cells of the
new file, every code cell was run in K1 from its saved source and the printed results sit in a markdown cell under
it). Spec: `L2-210-localising-span-spec.md` (§0–§4 round 1, Addenda A–C rounds 2–4). Built by a sonnet
wolfbook-builder in four rounds, planned and reviewed by the main session. Conventions as in `L2-210-untwist-CG.wb`
(new conventions, GT generators `JJ`, h = 1/2, θ = (2/7, −3/5), diagonal frame). Twist points P1: z = (3/2,5/4,2/3),
w = (5/2,7/3,3/5); P1′: z = (3/2,5/4,3/7), w′ = (7/4,−2/5,9/7). **Every statement below is exact and identical at
P1 and P1′ unless marked [numeric].**

## 0. Question and answer

**Paul's request.** A gl(3) analogue of the gl(4) rectangular result (`L2-gl4-rect-untwist.md`): a covector family
whose overlap with every eigenstate is a Q-number built from a number of SoV brackets / transfer-matrix data that
scales like the boxes of the Young diagram, not like dim V (a determinant is not required), and which untwists
cleanly, i.e. every CG datum of adj⊗adj at leading order in ε with K/S = O(ε) for all five irreps.

**Answer (negative, structural).** Take the base ⟨Ω_w⊗Ω_w| (overlap det₆) and the commutative algebra 𝒜 of
transfer matrices. At L = 2 the degrees of t₁(u), t₂(u) leave only six non-scalar coefficient operators
a₁..a₆ (t₁: 2, t₂: 4; t₃ is the quantum determinant, scalar), so the structured levels are the linear span D₁
(dim 7), the quadratic D₂ (dim 24: four exact quadratic relations), the cubic D₃ (dim 48); D₄ is all 64
functions on the eigenstates, i.e. Lagrange interpolation, the cost Paul excludes. Rational forms P/Q are
non-vacuous for (P,Q) ∈ D₂×D₂ (47 < 64 conditions) and D₃×D₁ (54 < 64) only. Within these levels the product
covectors v⊗Ω₁ (and the two-site variants tried) that are Q-numbers are exactly:

| level | reachable | not reachable |
|---|---|---|
| polynomial D₁ | Ω₁, Ω₁T₁₁^{(1)}(θ₁−2h) | |
| polynomial D₂ = D₃ | + Ω₁T₂₁^{(1)}(θ₁−2h) = −(128/75)·Ω₁ρ(1/z) at P1 (ℓ₁ = ℓ₂ = 3); the §5 product point Ω_{(w₁,y*)}ρ(1/(w₂−z)) ⊗ Ω_{(w₁,y*)} (degree 2, 18 terms) | polynomial torus Ω₁ρ(w₂′−z), ρ(z), resolvent Ω₁ρ(1/(w₂′−z)), second vacuum Ω_{w′}, Ω_{(w₁,y)} |
| rational (2,2) | + the common-w₁ line Ω_{(w₁,y)} on either site (dim 6 / 5 solution spaces), products Ω_{(w₁,y₁)}⊗Ω_{(w₁,y₂)}, the ω₂-lines Ω₁ + κω₂ | polynomial torus, resolvent, different-w₁ vacuum Ω_{w′}, every torus × y-deformed-vacuum target (T1–T4 of spec §6, seven more in the exploratory scan) |
| rational (3,1) | only the polynomial members (Q ∈ D₁ times the same X) — not even the (2,2) y-family | everything else |

So the gl(3) adjoint has, on a single vacuum, **one free rotation-preserving deformation** (the second twist
parameter y at common w₁, which moves only the E₂₃ Gauss parameter) and **no torus with a free point**. The
only tori in the localising span at θ₁−2h are d = 1 and d = 1/z (exact solve; T₁₁^{(1)} never enters), and the
non-localising T₁₀, T₂₀ create no torus line either (effective sum-of-products covector, exact solve).
The gl(4) mechanism needs both ingredients at a free point; the gl(3) adjoint supplies only the first.

## 1. Facts worth keeping

1. **T₂₁(θ₁−2h) is the inverse ket-twist torus on site 1.** ⟨Ω⊗Ω|T₂₁(θ₁−2h) = (Ω₁T₂₁^{(1)})⊗Ω₁ with
   T₂₁^{(1)} diagonal in the GT weight basis, eigenvalue −3h³χ₃²·Π_k z_k^{−μ_k} = −(3/8)Π z_k^{2−μ_k} on every
   weight (symbolic z). The gl(3) adjoint shadow of gl(4)'s t₂(θ₁−h) ∝ ρ₁(z)⁻¹.
2. **T₁₁(θ₁−2h) is not a torus.** T₁₁^{(1)} is diagonal on the six 1-dim weights with eigenvalue (3/4)χ₃/z_m
   (m the index with μ_m = 2) and a non-scalar 2×2 block on (1,1,1): [[3(z₁+z₂)z₃/8, 3(z₁−z₂)z₃/8],
   [(z₁−z₂)z₃/8, z₁z₂/2 + (z₁+z₂)z₃/8]] (trace e₂/2, det (3/16)e₁e₃). Matches none of ê₁(1/z), ρ(1/z)ê₁(z),
   e₂{z^μ}, e₁, e₂ of the complementary multiset {z^{2−μ}}.
3. **b = r₂ + r₁r₃ = (w₁−z₁)/(w₁−z₃) is exactly y-independent** (main session by hand; the notebook sees it as
   the undone covector of (Ω_{(w₁,y)}, Ω_w) having only the weights (2,1,0), (2,0,1)). Vacua of common w₁ differ
   only in the E₂₃ parameter c = r₃(y). Hence the y-family reads only [4,1,1]: it does untwist cleanly
   ([numeric] K/S = 1.1·10⁻⁶, 2.2·10⁻⁶, 4.4·10⁻⁶ at ε = 1,2,4·10⁻⁶; [4,1,1] Richardson limit 46.67453891138
   = the chain, φ′: −12.35023562230), but S ≡ 0 for [3,3,0], both octets and the singlet.
4. **The commutator direction depends on the site-1 y only.** For n_i = n(w₁,y_i), n′ = Δ⁻¹n₁Δn₂⁻¹:
   n′₁₃ − n′₁₂n′₂₃ = (δ₃/δ₁−1)(w₁−z₁)z₂₃/((w₁−z₃)(y₁−z₃)) (main session; the special case y₁ = y₂ is cell 14 of
   the gl(4) notebook). Freeing y₁ would make the undone covector generic at ε = 0 — but
   Ω_{(w₁,y₁)}ρ(1/(w₂−z))⊗Ω_{(w₁,y*)} is not a (2,2) or (3,1) rational Q-number.
5. **ω₂ = Ω₁ρ(exp(γ(E₂₃ − aE₁₃))), a = r₁ = (w₁−z₁)/(w₁−z₂), γ = −r₃(w₂)** (exact; the spec's pure-E₂₃ guess
   was wrong: n⁻¹e₂₃n = e₂₃ − a e₁₃). ⟨Ω⊗Ω|X = v⊗Ω₁ + v′⊗ω₂ exactly for X ∈ {1,T₁₀,T₁₁,T₂₀,T₂₁}(θ₁−2h), v′ ≠ 0
   only for T₁₀, T₂₀; the five effective undone covectors ṽ = v + v′ρ(exp(−γ(E₂₃−aE₁₃))) have rank 5.
6. **Enlarged 8-operator family** (+ t₂(θ₁), t₂(θ₂+h), t₂(θ₂+h)t₁(θ₂−h)): site-2 column space of dimension 5;
   1, T₁₀, T₁₁, T₂₀, T₂₁ stay on span{Ω₁, ω₂}; t₂(θ₁) (rank 2), t₂(θ₂+h) (rank 1), t₂(θ₂+h)t₁(θ₂−h) (rank 1)
   each bring one new site-2 direction, none a common-w₁ vacuum; the two site-2 localisers' directions are
   not rotations at all.
7. [numeric, P1, pr = 200] Every solved target satisfies ⟨Φ|Ψ_B⟩Q^B = P^B⟨Ω⊗Ω|Ψ_B⟩ on the 21 HW-sector states
   to < 10⁻¹⁹²; a 10⁻³ perturbation of one coefficient is detected at 3.5·10⁻³. In the regime the (P,Q)
   coefficients of the y-family stay O(1) while min|Q^B| ~ ε² (the cancellation sits in the denominator).

## 2. Cost statement

For the targets that are Q-numbers: the six eigenvalue coefficients (t₁ at 2 points, t₂ at 4), det₆(w), and one
polynomial (≤ 24 terms) or rational (24 + 24 terms) function of them — the same count for any representation at
L = 2. The failure is not cost: the needed covectors are simply not in the algebra below the vacuous level.

## 3. What is left (for Paul)

1. A different base than a single vacuum: sums of two vacuum determinants ⟨Ω_w⊗Ω_w|P + ⟨Ω_{w′}⊗Ω_{w′}|P′ (two
   det₆'s — within budget) with the target v⊗Ω_w; untested, and the non-vacuity bound is tighter (48 + 48 > 64
   for D₂ × D₂, so only D₁ × D₂-type combinations are structural).
2. Non-product covectors with a short rotation decomposition: the 8-operator family of item 6 has a 5-dim
   site-2 column space of which three directions are not rotations, so the undo is not available there; the
   5-dim θ₁−2h family is undoable (item 5) but its leading-order limit can cover at most 4 of the 5 non-top slice
   directions (dimension count: the family contains the top).
3. The structural reason: for the 6-rep of gl(4) the one-site fused Lax at the localising point is Λ² of a
   group element (t₂ = torus, 1−ωρ rotation-preserving); for the adjoint of gl(3) only T₂₁ is a torus and the
   rotation-preserving deformation moves one nilradical direction. Reps where the localising span contains a
   torus LINE (fundamental-type / rectangular) are the ones where the paper's mechanism works. Worth testing
   next: gl(3) [S,S,0] (rectangular, abelian nilradical — the direct gl(3) analogue of the paper) rather than
   the adjoint, and gl(4) [S,S,0,0] at S = 2 with the same (2,2) scan (settles its open Q-number form).
4. Rerun the notebook in VS Code to store native outputs (round 1–4 ≈ 1 minute of kernel time in total).
