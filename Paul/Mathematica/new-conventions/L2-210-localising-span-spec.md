# Spec: a clean-at-leading-order Q-number family for gl(3) adj⊗adj via the degree-≤2 transfer-matrix span (2026-10-10)

Delegated to a wolfbook-builder subagent (sonnet); planned and reviewed by the main session. **New notebook:**
`Paul/Mathematica/new-conventions/L2-210-localising-span.wb` (create it; self-contained; first code cell
`Quiet[ClearAll["Global`*"]]` as in every new-conventions notebook). Reference cell sources are in the two
scratchpad files named in the task message (`ref-L2-210-untwist-CG-cells1-12.txt`: generators, metric, exact
transfer-matrix coefficient matrices `tBasisC`, `oneSiteVac`, `ρdiag`, HW kets/projectors, `branchesAt`;
`ref-L2-rotation-sov-cells86-103-107-112.txt`: rank scans, product-locus solve, torus analysis). Copy what you
need; do not evaluate cells of other notebooks.

## 0. The problem and the idea (Paul's request)

gl(4) rectangular [S,S,0,0] untwists cleanly at leading order because its covector family has TWO independent
localising deformations (a torus insertion κ and the ω-deformation), both untwisting with z
(`L2-gl4-rect-untwist.md`). In gl(3) adj⊗adj the only known Q-number product family (det₆ × det₂,
`Two-rotations-det-overlap.md` §5) has the site-2 vacuum invariant y slaved to the site-1 resolvent point,
y = w₂(t−1)/(t−2), so y cannot untwist and the E₁₃ direction of the undone covector drops one order
(`Untwist-two-rotation-notes.md` §2.3). Paul does **not** need the overlap to be a determinant any more. He needs:

1. a covector family whose overlap with every eigenstate is a Q-number built from a bounded number of SoV
   brackets / transfer-matrix eigenvalues — the count must scale like the number of boxes of the Young diagram,
   not like dim V (so no Lagrange interpolation over all weights);
2. clean untwisting: every CG datum of adj⊗adj read off at leading order in ε, K/S = O(ε) for all five irreps
   (both octets and the singlet included).

**Idea.** Drop the product/determinant requirement and use the whole commutative algebra. At L = 2,
t₁(u) = Σ_{j=0}^{2} u^j t₁^{(j)} and t₂(u) = Σ_{j=0}^{4} u^j t₂^{(j)} with scalar leading coefficients
(t₁^{(2)} = χ₁·1, t₂^{(4)} = χ₂·1) and t₃(u) scalar (quantum determinant). So the LINEAR span of all
transfer matrices at all points is D₁ = span{1, a₁..a₆}, a = (t₁^{(0)}, t₁^{(1)}, t₂^{(0)}, t₂^{(1)}, t₂^{(2)},
t₂^{(3)}), 7-dimensional, and the QUADRATIC span D₂ = D₁ + span{a_i a_j, i ≤ j} has dimension ≤ 28 —
rep-independent at L = 2 (every fused T_λ(θ_s + mh), λ ≤ (2,2), lies in D₂). The cubic span is ≥ 64 = dim(V⊗V),
i.e. everything, so only D₁, D₂ carry structure. For any X ∈ D₂ the covector ⟨Ω_w⊗Ω_w|X has overlap
X^B·⟨Ω_w⊗Ω_w|Ψ_B⟩ with X^B the same polynomial in the eigenvalue coefficients and ⟨Ω⊗Ω|Ψ_B⟩ = det₆: a Q-number
with at most 6 + 21 terms, whatever the rep. **Questions:**

- **Q1 (finite twist, exact).** Which site-1 covectors v are reachable as ⟨Ω_w⊗Ω_w|X = v ⊗ Ω_w, X ∈ D₂
  (the "site-1 localising subspace" L₁)? Its dimension ℓ₁ ≤ 8. Known members: Ω₁, Ω₁T₁₁^{(1)}(θ₁−2h),
  Ω₁T₂₁^{(1)}(θ₁−2h) (so ℓ₁ ≥ 3). The best case ℓ₁ = 8 means every site-1 covector ⊗ vacuum is a Q-number of
  degree ≤ 2. Same for site 2 (L₂, ℓ₂). The 5-dim span at the single point θ₁−2h gave only the isolated product
  point c* (`Two-rotations-det-overlap.md` §5.3); D₂ mixes all points.
- **Q2.** Are the specific targets in L₁: the polynomial-torus translate Ω₁ρ(diag(w₂′−z)) with w₂′ free
  (the gl(3) analogue of the paper's κ-torus; untwists cleanly, see §3), the plain ket-twist torus Ω₁ρ(diag z),
  the resolvent translate Ω₁ρ(diag 1/(w₂′−z)), and a second vacuum Ω_{w′}?
- **Q3 (untwist).** For the family found, in the joint regime (z, w, w₂′ → 1 together), are all seven site-1
  weights of the undone covector O(1) with full rank, and is K/S = O(ε) for all six HW branches?

Why the torus translate must untwist cleanly if it is a Q-number (main session, from the exact cell-14 identity
of `L2-gl4-rect-untwist.wb`): for n = n(w₁, y) the vacuum's Gauss matrix and Δ = diag(δ), the undone
n′ = Δ⁻¹nΔn⁻¹ has n′₁₃ − n′₁₂n′₂₃ = (δ₃/δ₁ − 1)(w₁−z₁)z₂₃/((w₁−z₃)(y−z₃)). With y = w₂ = 1+εω₂ (the vacuum's
own second twist parameter, untwisting) and δ_k = w₂′−z_k = ε(k−φ_k) every factor is O(ε) in numerator and
denominator: the E₁₃ direction is O(1). With y pinned at y₀ ≠ 1 (the det₂ family) the factor (y−z₃) is O(1)
and the direction is O(ε) — that was the disease.

## 1. Conventions (as in `L2-210-untwist-CG.wb`; copy its cells 3–9 from the reference file)

New conventions, gl(3), Λ = (2,1,0), L = 2, plain GT generators `JJ` (exact rationals), invariant metric `gMet`,
h = 1/2, θ = (2/7, −3/5) (t = 62/35). Diagonal frame, ket twist diag(z) on the right. Finite-twist point P1:
z = (3/2, 5/4, 2/3), w = (5/2, 7/3, 3/5); control point P1′: z = (3/2, 5/4, 3/7), w′ = (7/4, −2/5, 9/7)
(at P1, z₃ = 1/z₁ — do every exact statement at both points). `tBasisC[a, I]` are exact coefficient matrices;
**keep everything exact** (no `N[..]`) in Tasks 1–2: define `tZexact[zz]` = the exact analogue of `tZof`.
One-site vacuum `oneSiteVac[zz, ww]` = e₁ᵀexp(r₁E₁₂)exp(r₂E₁₃)exp(r₃E₂₃); `ρdiag[d]` the torus;
Ω = Ω₁⊗Ω₁ with Ω₁ = oneSiteVac[z, w]. Fused shapes for cross-checks: T₁₀ = t₁, T₁₁ = t₂, T₂₀(v) = t₁(v)t₁(v+h) −
t₂(v+h), T₂₁(v) = t₂(v)t₁(v+h) − t₃(v+h), T₂₂(v) = t₂(v)t₂(v+h) − t₃(v+h)t₁(v) (check this one against
T_λ(x) = det_{i,j} t_{λ′_j+i−j}(x+(i−1)h) before using it).

## 2. Tasks (each its own cell(s) with a `Grid` of checks; `checks` an Association; time one case before loops)

**Task 0 — setup checks.** gl(3) relations, metric, HW kets (6; sectors 420, 411, 330, 321×2, 222), qdet check,
[t₁(u), t₂(v)] = 0 exactly. Build the exact coefficient operators a₁..a₆ at P1 and verify: t₁^{(2)} = χ₁·1,
t₂^{(4)} = χ₂·1, t₃(u) ∝ 1 at two u's; all a_i commute (exact); the 28 operators {1, a_i, a_ia_j} are linearly
independent as 64×64 matrices (rank 28 of the 28×4096 flattening) — report the rank; if < 28, say which
relations. Also report the rank of the cubic span {a_ia_ja_k} ∪ D₂ (expect 64: vacuous, so D₂ is the last
structured level). Verify T₁₁(θ₁−2h), T₂₁(θ₁−2h) ∈ D₂ (exact solve), as a consistency check of the basis.

**Task 1 — the localising subspaces L₁, L₂ (exact, P1 and P1′).** Let Φ_X := ⟨Ω⊗Ω|X reshaped 8×8 (row =
site 1, column = site 2). Φ_X = v ⊗ Ω₁ ⇔ Φ_X·k = 0 for every ket k in the 7-dim annihilator of Ω₁ (NullSpace of
{Ω₁}). With X = c₀ + Σc_i a_i + Σ_{i≤j} c_ij a_ia_j this is a homogeneous linear system (56 equations, 28
unknowns): solve exactly (`NullSpace`), report ℓ₁ = dim of the solution space and a basis of L₁ ⊂ C⁸ (the v's,
v = Φ_X·k₀/(Ω₁·k₀) for a ket k₀ with Ω₁·k₀ ≠ 0). Report: ℓ₁; rank of the D₂-coefficient solution space; which
site-1 weights L₁ covers (per weight (2,1,0),(1,2,0),(2,0,1),(1,1,1),(1,0,2),(0,2,1),(0,1,2): rank of the
projected L₁); whether L₁ contains Ω₁, Ω₁T₁₁^{(1)}, Ω₁T₂₁^{(1)} (controls, must be in) and whether these three
span it; the same within D₁ only (ℓ₁(D₁)). Then the same for site 2 (Φ_X = Ω₁ ⊗ u; ℓ₂, L₂ ⊂ C⁸). Repeat at P1′.
Also record the "column structure": for each basis operator, rank of Φ_X; the dimension of the joint site-2
column space of {Φ_X : X ∈ D₂} (span of all columns) and of its site-1 row space (if both are 8, say so).

**Task 2 — targets (exact, P1 and P1′).** For each target site-1 covector v below test v ∈ L₁ (exact membership:
rank of L₁-basis ∪ {v} vs rank of L₁-basis); if in, solve for the D₂ coefficients of a representative X (the
solution is unique modulo the kernel {X ∈ D₂ : ⟨Ω⊗Ω|X = 0} — report that kernel's dimension too) and print them:
  (a) v = Ω₁ρdiag[w₂′ − z] (polynomial torus) for w₂′ ∈ {3/7, −11/5, 9/4, 7/3, 13/10} — if in L₁ for all,
      fit the coefficients as polynomials in w₂′ (expect degree ≤ 3 = |Λ|, the three operators ê₁, ê₂, ê₃ with
      ê₃ = ρdiag[z]) on four values and check the fifth;
  (b) v = Ω₁ρdiag[z] and Ω₁ρdiag[1/z];
  (c) v = Ω₁ρdiag[1/(w₂′ − z)] (resolvent) for two w₂′ (expect: not in — item 7 of `Rotations-SoV-Q.md`
      needs T₁₀⁻¹ already at L = 1 — but test);
  (d) v = oneSiteVac[z, w′] (a second vacuum at the same z) and oneSiteVac[z, {w₁, y, 1}] for y ∈ {3, −1/2}
      (same w₁, free second parameter);
  (e) the §5 product point: v = Ω_{(w₁,y*)}ρ(diag 1/(w₂−z)) with y* = w₂(t−1)/(t−2) and the vacuum
      Ω_{(w₁,y*)}⊗Ω_{(w₁,y*)} as base — must be in (it is c*·T(θ₁−2h)); regression of the method.
If ℓ₁ = 8 all of (a)–(d) are in automatically; then just print the coefficients and go on. If ℓ₁ < 8 and (a) is
not in L₁: also test the **sum-of-products** version: the subspace M ⊂ ⟨Ω⊗Ω|D₂ whose site-2 column space lies in
span{Ω₁, ω₂}, ω₂ := oneSiteVac[z, {w₁, z₂, 1}] (the second direction of §5.2; check it is a rotation
e₁ᵀexp(r₁E₁₂)exp(r₂E₁₃) with r₃ = 0), dimension of M and of its site-1 projections; then stop after Task 2′ and
report — the main session decides the regime for a sum of two products.

**Task 2′ — Q-number check on the chain (numeric, pr = 200, P1).** For the representative X of target (a) at
one w₂′ (or of the best target found), verify on the 21 states of the HW sectors (`branchesAt` at the exact
P1 twist, no ε): ⟨v⊗Ω₁|Ψ_B⟩ = X^B⟨Ω⊗Ω|Ψ_B⟩ with X^B = c₀ + Σc_iα_i^B + Σc_ijα_i^Bα_j^B, α^B the coefficients of
the state's stored t₁(u), t₂(u) polynomials; report the max relative deviation (expect 10⁻¹⁵⁰ or better).
State in a markdown cell: the overlap is det₆(w) × a polynomial of degree 2 in the six eigenvalue coefficients
(each of which is a Wronskian/Cramer expression in the dual Q's — not re-derived here).

**Task 3 — untwisting the polynomial-torus family (if (a) ∈ L₁).** Family: v₂ = Ω_w, v₁ = Ω_wρdiag[w₂′ − z],
φ_u := v₁ρ(n)⁻¹ (n the vacuum Gauss matrix, `u12rNinv`-pattern with y = w₂). Regime: z = 1+εφ, w = (1+εω₁,
1+εω₂, 1+ε·5/13), w₂′ = 1+εk; φ = (1/3, −2/7, 5/11) (control φ′ = (−2/5, 3/8, 7/11)); grid ω₁ ∈ {3/4, −5/3},
ω₂ ∈ {−1/5, 9/7, 5/3}, k ∈ {−12/5, 7/2, 1/8} (never k = φ_j, k = ω₂, ω_i = φ_j); ε ∈ {1,2,4}·10⁻⁶.
  (i) Exact ε-structure (symbolic ε, as `L2-210-untwist-CG.wb` cells 11–12): per site-1 weight the leading order
      and the rank of the ε⁰ coefficient vectors over the 18-point grid. **Expectation:** all seven weights at
      ε⁰, rank 1 on the 1-dim weights, rank 2 on (1,1,1) (both directions u = e₁ᵀE₁₃ and u′ = e₁ᵀE₂₃E₁₂ seen),
      φ_u(0)/δ^Λ = e₁ᵀρ(n′₀) with n′₀ = Δ₀⁻¹n₀Δ₀n₀⁻¹, Δ₀ = diag(k−φ), n₀ the ε = 0 vacuum Gauss matrix, and
      n′₀,₁₃ − n′₀,₁₂n′₀,₂₃ ≠ 0 (print it). Control: w₂′ = 3/7 fixed (k-regime off) — expect the old depth
      table (ε⁰, ε¹, ε¹, ε², ε³, ε³, ε⁴) since Δ → scalar.
  (ii) Chain side: for each ε and grid point, the 6 HW-limit branches (`branchesAt`), R_C = ⟨v₁⊗v₂|Ψ_B⟩,
      S = ⟨v₁⊗v₂|P_HWΨ_B⟩, K = R_C − S; check S = s⟨φ_u⊗e₁ᵀ|P_HWΨ_B⟩ exactly (undo); ε-exponents of |S|, |K|,
      max|K/S| over the grid at each ε. **Expectation:** K/S = O(ε) for all six branches (singlet and both
      octets included). If a branch has K/S ↛ 0, run the `showContam` breakdown (which site-1 weights carry K)
      and report; do not fit.
  (iii) CG at leading order from the overlap alone: C := R_C/(δ^Λ[φ_u(0)/δ^Λ]_slice) for the 1-dim slices
      ((2,0,1) for [4,1,1], (1,2,0) for [3,3,0], (0,1,2) for the singlet): spread over the grid and error vs the
      chain slice ⟨m_slice⊗e₁|P_HWΨ_B⟩·(normalisation) at two ε's (both ∝ ε), Richardson limits (compare with
      `L2-210-untwist-CG.md` §5: [4,1,1] 46.67453891, [3,3,0] 10.92853892, singlet −5118.538314 at P1).
      Octets: the 2×2 system on two grid points with different k (the (1,1,1) coefficient vectors must be
      independent) → both slice components, error vs chain, and the mixing R_c (P1 chain values −1.694, 1.096).
  (iv) The Q-number at small ε: for one grid point and ε = 10⁻⁶, recompute the exact D₂ coefficients c(ε)
      of X (exact rationals in ε — they will have poles in ε; that is expected) and verify R_C = X^B·⟨Ω⊗Ω|Ψ_B⟩
      on the 6 branches (relative deviation). Report the ε-exponents of the c's (the cancellation depth).
  Then repeat (i)–(ii) for φ′ and, if cheap, for the second-vacuum family v₁ = Ω_{w′} (w′ = 1+εω′,
  ω′ = (−12/5, 7/2, 1/8)) — expect the same cleanliness.

**Task 4 — cost statement (markdown).** Number of distinct transfer-matrix data needed: 6 eigenvalue
coefficients (t₁: 2, t₂: 4), i.e. t₁ at 2 points and t₂ at 4 points, plus det₆ — for ANY rep at L = 2 the
count is the same (the degrees of t_a(u) depend only on L and gl(N)); say explicitly how many nonzero c's the
representative X has. Flag for the follow-ups (not in this notebook): [3,1,0] with the same scan; gl(4)
[S,S,0,0] at S = 2 with D₂ of t₁..t₄ (settles the open Q-number form at S = 2).

## 3. Reporting

End with a markdown cell `### gl(3) adj⊗adj localising-span result (2026-10-10)`: ℓ₁, ℓ₂ at P1 and P1′ with
the weight coverage of L₁; the membership table for targets (a)–(e); the representative coefficients of (a)
and their w₂′-polynomial form; the Task 2′ deviation; the per-weight order/rank table of Task 3(i) and the control;
the K/S table with ε-exponents; the CG-extraction errors and Richardson limits vs the chain; the octet mixing;
every structural surprise. Save (`wolfbook_saveNotebook`). If VS Code does not execute cells of the new file,
run each cell's source in K1 via `wolfbook_evaluateExpression` and store the printed results in a markdown cell
under the code cell (say so in the report). Report back concisely in the same structure, numbers not adjectives.
Do not edit other notebooks, `paper/*.tex`, or git. If Task 1 gives ℓ₁ < 8 and no target of Task 2 is in L₁,
stop after Task 2 (with the sum-of-products diagnostics) and report; the regime for that case is the main
session's call.

## 4. Gotchas (all have bitten before)

- The shared kernel K1 is the only kernel; your first cell wipes the other notebooks' live state (files intact).
- Never write `*)` inside a comment. Do not use `D`, `E`, `N`, `I`, `C`, `K` as symbols (use `KK`, `kk`).
- Memoized functions cache bad results; after fixing a definition, re-run it and everything downstream.
  `Do[DownValues[f] = …, {f, {…}}]` purges nothing (HoldAll) — purge with literal symbols.
- Semicolons between statements in multi-line cells. `Max[{}]` = −∞; floor at 0.
- Exact zero tests on rationals with `===`; numerics: report `Precision`; `SetPrecision` before `Norm` of
  near-zero vectors; never `PossibleZeroQ` on numerics.
- Wolfbook replies truncate at `{`; read long results with `StringReplace[ToString[.., InputForm], "\n" -> " "]`.
- The vacuum degenerates when a w_i equals a z_j (also w₂′ = z_j, y = z_j); check at every new point.
- Exact 64×64 rational linear algebra with the 28-term D₂ basis is fine; exact 4096-column rank computations
  (Task 0) may take a minute — time one before looping. Use `SparseArray` for the Kronecker products.
- `oneSiteVac[zz, ww]` ignores w₃; `reg` in the reference file sets w₃ = 1+2ε/9 — use 1+5ε/13 (2/9 collides with
  ω₂ = 2/9 if you ever add it).

## 5. Addendum A (main session, after round 1): rational Q-numbers P/Q of degree (2,2)

**Round-1 result.** ℓ₁ = ℓ₂ = 3 for D₂ and for D₃ (rank 48); D₄ is the full 64-dim algebra (vacuous). L₁ =
span{Ω₁, Ω₁T₁₁^{(1)}, Ω₁T₂₁^{(1)}} at θ₁−2h — D₂ mixes nothing new in. Exact positive: Ω₁ρdiag[1/z] =
−(128/75)·Ω₁T₂₁^{(1)}(θ₁−2h) at P1 (−3h³χ₃² ρ(1/z) on the HW weight), i.e. the T₂₁ localiser IS the inverse
ket-twist torus on site 1 — the gl(3) adjoint shadow of gl(4)'s t₂(θ₁−h) ∝ ρ₁(z)⁻¹. The polynomial torus
Ω₁ρdiag[w₂′−z], the resolvent, a second vacuum, ρ(z): none in L₁.

**Why a rational form is the natural next level.** At L = 1 the resolvent is ρ((1−sGᵀ)⁻¹) ∝ 1 + (s/h)T₁₀ +
(s²/2h²)(T₁₁+T₂₀) + (s³/3h³)T₂₀T₁₀⁻¹T₁₁ (`Rotations-SoV-Q.md` item 7): a RATIONAL function of the transfer
matrices, numerator and denominator of low degree. The paper's κ-family coefficients κ²/Π_k(1−κz_k) are rational
too. So test, for a target site-1 covector v (so Φ_v := v⊗Ω₁ as a 64-covector): are there P, Q ∈ D₂ with

    Φ_v · Q = ⟨Ω⊗Ω| P        (as covectors; equivalently X^B = P^B/Q^B on every eigenstate).

Counting: (P, Q) ∈ D₂×D₂ has 24+24 coefficients, minus scale = 47 < 64 equations, so a solution is structural
(generic targets have none; verify this with a random v as control). Also try the smaller pairs (D₁, D₂),
(D₂, D₁), (D₁, D₁) for sparser forms (at L = 1 the denominator is T₁₀, one term). Do NOT go to D₃ in either
slot (48 + 24 − 1 ≥ 64: vacuous).

**Task 5 (exact, P1 and P1′).** Build the 64 × 48 exact linear map (P, Q) ↦ Φ_v·Q − ⟨Ω⊗Ω|P over the 24-element
independent bases of D₂ (use the independent 24, not the 28 with relations). For each target report dim of the
solution space, a representative (P, Q) (nonzero-coefficient counts), and whether Q is invertible on the chain
(Q^B ≠ 0 for all 64 states of the exact P1 chain — compute Q^B from `branchesAt`-type eigenvalue coefficients at
pr = 200 on the 21 HW-sector states, and also exactly as the determinant of Q as a 64×64 matrix). Targets:
  (a) Ω₁ρdiag[w₂′ − z], w₂′ ∈ {3/7, −11/5, 9/4} — if solvable, the w₂′-dependence of (P, Q) (polynomial?);
  (c) Ω₁ρdiag[1/(w₂′ − z)], w₂′ ∈ {3/7, 9/4};
  (d) Ω_{w′} (w′ = (7/4, −2/5, 9/7)) and Ω_{(w₁, y, 1)}, y ∈ {3, −1/2};
  (f) the ω₂-line on site 2: Φ = Ω₁ ⊗ (Ω₁ + κω₂), κ ∈ {1, −2}, and on site 1: (Ω₁ + κω₂) ⊗ Ω₁;
  (g) controls: v = Ω₁ρdiag[1/z] (must solve with Q = 1), and a random rational 8-vector (must NOT solve).
If some target solves: verify on the chain (Task 2′ pattern) that ⟨v⊗Ω₁|Ψ_B⟩·Q^B = P^B⟨Ω⊗Ω|Ψ_B⟩ on the 21 states
(rel. dev.), then run **Task 3 of §2 for that family** with the Q-number X^B = P^B/Q^B (in 3(iv) recompute
(P, Q) exactly at the regime point and report the ε-exponents of their coefficients). Priority if several
solve: (d) second vacuum first (cleanest untwisting: two independent limiting rotations n₀(ω′;φ), n₀(ω;φ)),
then (a). If none solves, report the solution-space dimensions (all should be 0 except (g) control 1) and stop.

## 6. Addendum B (main session, after round 2): free the site-1 y of the §5 product point

**Round-2 result.** Rational (2,2) forms exist exactly for the y-deformation of the vacuum at common w₁ —
Ω_{(w₁,y)} on either site (solution dim 6 / 5) — and for nothing else tried. The undone covector of
(Ω_{(w₁,y)}, Ω_w) has only the weights (2,1,0), (2,0,1): because b = r₂ + r₁r₃ = (w₁−z₁)/(w₁−z₃) is exactly
y-independent (main session, by hand from `oneSiteVac`), vacua of common w₁ differ only by the E₂₃ parameter
c = r₃(y). So this family reads [4,1,1] only.

**The fix.** For the §5 product point (v₁, v₂) = (Ω_{(w₁,y*)}ρ(D), Ω_{(w₁,y*)}), D = diag 1/(w₂−z),
y* = w₂(t−1)/(t−2), the undone matrix n′ = D⁻¹n₁Dn₂⁻¹ has, for n₁ = n(w₁,y₁), n₂ = n(w₁,y₂) with a common a, b
and c_i = r₃(y_i):
    n′₁₂ = a(δ₂/δ₁−1),  n′₂₃ = c₁δ₃/δ₂ − c₂,  n′₁₃ − n′₁₂n′₂₃ = (δ₃/δ₁−1)·[b − a c₁ z₂₁/z₃₁] = (δ₃/δ₁−1)(w₁−z₁)z₂₃/((w₁−z₃)(y₁−z₃)).
Only the SITE-1 y enters the commutator direction. In the joint regime (z, w₁, w₂ → 1) it is O(ε) iff y₁ stays
pinned at y* = O(1) ≠ 1 (the disease), and O(1) if y₁ = 1 + εη₁ untwists — while y₂ may stay at y*, since the
site-2 rotation is undone exactly whatever it is. So the candidate family is

    v₂ = Ω_{(w₁, y*)},   v₁ = Ω_{(w₁, y₁)} ρdiag[1/(w₂ − z)],   y₁ free,   y* = w₂(t−1)/(t−2),

a product of two distinct rotations whose undone covector at ε = 0 is a generic unipotent (a′, c′, b′−a′c′ all
O(1) and varying independently over the (ω₁, ω₂, η₁) grid) — every CG slice at ε⁰, including the second
octet direction and the singlet. The y₁-deformation is a (2,2)-rational Q-number on the vacuum and the product
point is polynomial (round-1 (e) regression), so the composite is in 𝒜; the question is its degree.

**Task 6 (exact, P1 and P1′): (P, Q) scan with two-site targets.** Same machinery as Task 5, target Φ a
general 64-covector:
  T1: Φ = Ω_{(w₁,y₁)}ρdiag[1/(w₂−z)] ⊗ Ω_{(w₁,y*)}, y₁ ∈ {3, −1/2, 7/5} (w₂ = the base w₂, y* = w₂(t−1)/(t−2);
      at P1 y* = −63/8). Pairs (D₂,D₂), (D₂,D₁), (D₁,D₂), (D₁,D₁); control y₁ = y* (must solve with Q ∝ 1).
  T2: the solvable locus: Φ = Ω_{(w₁,y₁)}ρdiag[1/(w₂′−z)] ⊗ Ω_{(w₁,y₂)} for (y₁, y₂, w₂′) ∈ {(3, y*, w₂),
      (3, 2, w₂), (3, y*, 9/4), (3, 9/4·(t−1)/(t−2)·… i.e. y₂ = y*(w₂′) with w₂′ = 9/4)} — i.e. is y₂ = y*(w₂′)
      still required, and is w₂′ free of the base w₂?
  T3: T1 with the polynomial torus ρdiag[w₂−z] in place of the resolvent.
  T4: the symmetric variant Φ = Ω_{(w₁,y₁)}ρ(D) ⊗ Ω_{(w₁,y₂)} with y₁ = y₂ = 3 ≠ y* (expect: no solution,
      consistent with §5.6b).
Report per target the solution-space dimension per pair, a representative (nonzero counts), det Q ≠ 0,
min|Q^B| on the 21 states, and the chain identity ⟨Φ|Ψ_B⟩Q^B = P^B⟨Ω⊗Ω|Ψ_B⟩ (relative deviation) for every
solved target. Random-covector control once (must not solve).

**Task 7 (if T1 solves): untwist.** Regime z = 1+εφ, w₁ = 1+εω₁, w₂ = 1+εω₂, w₃ = 1+5ε/13, y* = w₂y₀ with
y₀ = (t−1)/(t−2) (pinned), y₁ = 1+εη₁ (η₁ ∉ {φ_j, ω₂}); grid ω₁ ∈ {3/4, −5/3}, ω₂ ∈ {−1/5, 9/7, 5/3},
η₁ ∈ {−12/5, 7/2, 1/8}; φ and φ′; ε ∈ {1,2,4}·10⁻⁶. Then exactly Task 3 (i)–(iv) of §2 with this family:
(i) symbolic-ε per-weight orders and ranks (expect all seven weights at ε⁰, (1,1,1) rank 2; print a′₀, c′₀,
b′₀ − a′₀c′₀ in closed form in (φ, ω, η₁) if `Factor` gives something short); control η₁-regime off (y₁ = 3
fixed) → expect the (0,1,2) and second-(1,1,1) directions to drop to O(ε) as in `L2-210-untwist-CG.md` §1;
(ii) K/S for the six branches (expect O(ε) for all six, including the singlet and both octets — the decisive
test; if a branch fails, `showContam` breakdown, no fitting); (iii) CG at leading order from the overlap alone
with Richardson against the P1 references ([4,1,1] 46.67453891, [3,3,0] 10.92853892, singlet −5118.538314,
octet (u+u′)-functionals 184.3305612 / −238.5587523, mixing R = −1.694 / 1.096 — mind the normalisations of
`L2-210-untwist-CG.wb` cells 41–45 when comparing: compare the chain-side slice computed in THIS notebook, and
quote the references only as a sanity check), the octet 2×2 system from two grid points with different η₁;
(iv) exact (P,Q) at a regime point at three ε's: solution dimension, ε-exponents of the coefficients and of
min|Q^B|, chain check on the six branches.

**Task 8 (markdown).** Cost: the six eigenvalue coefficients (t₁ at 2 points, t₂ at 4) and det₆(w), one
rational function of degree (2,2) — the same count for any rep at L = 2; number of nonzero (P,Q) coefficients.
Update the final summary cell.

## 7. Addendum C (main session, after round 3): the last non-vacuous scans

**Round-3 result.** No (P,Q) ∈ D₂×D₂ for any target mixing a torus with a y-deformed vacuum; products of two
common-w₁ vacua are (2,2)-rational; the product point is polynomial (degree 2). So within polynomial degree ≤ 3
and rational (2,2) the reachable site-1 covectors on the vacuum are the E₂₃-line, {T₁₁, ρ(1/z)} at θ₁−2h, and
the slaved product point. Two structured levels are still untested and non-vacuous by counting:
(D₃, D₁): 48 + 7 − 1 = 54 < 64 — and this is exactly the shape of the L = 1 resolvent identity
ρ((1−sGᵀ)⁻¹) ∝ 1 + (s/h)T₁₀ + (s²/2h²)(T₁₁+T₂₀) + (s³/3h³)T₂₀T₁₀⁻¹T₁₁ (`Rotations-SoV-Q.md` item 7): numerator
cubic, denominator T₁₀ linear. And the torus content of the localising span itself.

**Task 9 — (D₃, D₁) rational scan (exact, P1 and P1′).** Independent basis of D₃ (rank 48: take the 24 of D₂
plus 24 independent cubic monomials a_ia_ja_k) in the numerator, D₁ (7) in the denominator; also (D₃, 1) i.e.
the polynomial cubic span as a control of the round-1 statement. Targets: (a) Ω₁ρdiag[w₂′−z]⊗Ω₁ (w₂′ ∈ {3/7,
−11/5, 9/4}); (c) Ω₁ρdiag[1/(w₂′−z)]⊗Ω₁ (two w₂′); T1 of §6 (y₁ ∈ {3, −1/2}); (d) Ω_{w′}⊗Ω₁ (different w₁);
controls: the product point (must solve with Q ∝ 1), Ω_{(w₁,3)}⊗Ω₁ (solves already at (2,2); record its
(3,1) dimension), a random 64-covector (must not solve). Report dimensions; for any solved target: a
representative with a SPARSE denominator if one exists (try Q restricted to {1, t₁^{(0)}, t₁^{(1)}} and to
each single a_i), det Q ≠ 0, min|Q^B| on the 21 states, chain identity deviation. (Do not go to (D₃, D₂) or
(D₄, ·): vacuous.)

**Task 10 — torus content of the localising span (exact, P1 and P1′).**
  (i) The one-site operator T₁₁^{(1)} := the site-1 operator with ⟨Ω⊗Ω|T₁₁(θ₁−2h) = (Ω₁T₁₁^{(1)})⊗Ω₁ (you have
      it from round 1): is it diagonal in the GT weight basis? Print its eigenvalues per weight as functions
      of z (compute at symbolic z if cheap — the one-site t₂(u) at u = −2h on the 8-dim adjoint — else at
      three rational z and rationalise), and the same for T₂₁^{(1)} (expect −3h³χ₃²·Π z_k^{−μ_k}). Compare
      T₁₁^{(1)} with the candidates ê₁(1/z) = Σ_kμ_k/z_k, ρ(1/z)ê₁(z), e₂ of the multiset {z_k^{×μ_k}}, and the
      complementary-torus coefficients of Π_k(1−κz_k)^{2−μ_k}; say which (if any) it is, and whether it is
      scalar on the dim-2 weight (1,1,1) (a torus must be).
  (ii) Exact solve: for which (κ₁, κ₂, d₁, d₂, d₃) is Ω₁(1 + κ₁T₁₁^{(1)} + κ₂T₂₁^{(1)}) = Ω₁ρdiag[d] (as
      8-vectors, d up to scale)? Expect only the trivial solutions (κ = 0; κ₁ = 0 with d = 1/z); report any
      nontrivial torus line (that would be the gl(3) κ-family and would send us straight to the untwist).
  (iii) The same question for the 5-dim span at θ₁−2h including the non-localising T₁₀, T₂₀, in the
      sum-of-products form Φ(c) = v(c)⊗Ω₁ + v′(c)⊗ω₂: is there a c with v(c) = Ω₁ρdiag[d] a torus translate
      AND v′(c) = Ω₁ρdiag[d]·(anything)? More simply: the effective undone covector ṽ(c) := v(c) +
      v′(c)ρ(exp(−γE₂₃)) with ω₂ = Ω₁ρ(exp(γE₂₃)) (verify this form of ω₂ exactly and give γ) — the set of c
      with ṽ(c) = Ω₁ρdiag[d] for some d. Exact solve in (c, d).
  (iv) The enlarged 8-dim linear family of `Two-rotations-det-overlap.md` §5.3 ({1, T₁₀, T₁₁, T₂₀, T₂₁ at
      θ₁−2h, t₂(θ₁), t₂(θ₂+h), t₂(θ₂+h)t₁(θ₂−h)}): the site-2 column space of ⟨Ω⊗Ω|X over it (dimension; the
      5 of §5.3 expected) and whether every column direction is a common-w₁ vacuum (lies on the E₂₃-line
      span{Ω₁, ω₂}) or not; if not, which operator brings in which new site-2 direction, and is it a rotation.

**Task 11 — only if Task 9 or 10 gives a free-parameter torus (a (3,1) form for (a) or (c), or a torus line
in (ii)/(iii)):** the untwist of §2 Task 3 for that family (regime z, w, w₂′ → 1 jointly; all seven weights
at ε⁰ expected; K/S for six branches; CG with Richardson; exact (P,Q) at a regime point). Otherwise stop
after Task 10 and update the final summary cell with: the complete map of what is reachable at each level
(poly ≤ 3, rational (2,2), rational (3,1)) and the torus content of the localising span.
