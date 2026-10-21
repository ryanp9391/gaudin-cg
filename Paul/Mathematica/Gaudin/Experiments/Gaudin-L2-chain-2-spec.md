# Spec W2: sl(3) Gaudin chain, L = 2, part 2 — the (∞,∞) two-vacua family, the [3,1,0] fourth localiser, closed forms (2026-10-16)

Delegated to a wolfbook-builder subagent (sonnet); planned by the main session. **New notebook:**
`Paul/Mathematica/Gaudin/Experiments/Gaudin-L2-chain-2.wb` (self-contained; rep-generic with `Λ` as the only switch; first code cell
`Quiet[ClearAll["Global`*"]]`; [2,1,0] first, then [3,1,0] by changing the switch; Task B is [3,1,0]-only, Task C runs several Λ). Read
first: `../Gaudin-sl3-round2-plan.md` §0, §1, §3, §8 and the round-1 plan `../Gaudin-sl3-CG-plan.md` §1 (dictionary, marked ✔ where
verified). **Template:** the round-1 notebook `Gaudin-L2-chain.wb` (48 cells, same directory) — read its JSON (`python -X utf8` script,
print each code cell's `value`) and COPY its setup cells (generators, metric, Lax `Lg`, the symbolic cdet normal ordering and `S[k]`,
Laurent coefficients `Hres`/`S3c`, `gc`/`PauxG`/`ΩD`/`rGauss`/`nMat3`/`nInv3`, `Pphys`-type construction and `x0G`, `hwKets`/`PHWc`,
eigenstates and `branchesAt`, the undo/K/S/CG cells of its Task 4) rather than re-deriving them. Its summary cell (last cell) states the
conventions that are now verified: cdet column-ordered with E_ij; H₁ = −[Σζ_iE1_ii + Σ_{ik}E1_ikE2_ki/θ₁₂] + (|Λ|Σζ + |Λ|²/θ₁₂)·id;
g′ = PauxG[ζ].gc[ω].PauxG[ζ]⁻¹; one-site vacuum rule e₁ρ(U-factor of M^{−T}) with M = left eigenvectors of g′; P_phys^G = (PauxG⁻¹)ᵀ.
Everything exact where possible; 200 digits where eigenstates enter; time one case before any loop; nothing over ~3 min without a timed
smaller run.

## 0. Points and conventions

Λ = {2,1,0} then {3,1,0}; θ = {2/7, −3/5}; P: ζ = {3/2, 5/4, 2/3}, P′: ζ″ = {2, 5/3, −1/2}; (ω₁, ω₁′) = (5/2, 7/4) base, (−2/5, 3)
second; reference finite ω₀ = {5/2, 7/3, 3/5}; regime ζ = εφ, φ = {1/3, −2/7, 5/11}, ω₁ = εω̂₁, ω₁′ = εω̂₁′, ω̂₁ ∈ {3/4, −5/3},
ω̂₁′ ∈ {−12/5, 7/2} (4 grid points; add (1/2, −1/3) and (−7/3, 5/4) if a weight needs more points to reach full rank), ε ∈ {1,2,4}·10⁻⁶,
pr = 200. Never ω₁ = ζ_j, never ω₁ = ω₁′ except in the control.

## 1. Tasks (each: code cell + markdown result cell with a `Grid`; `checks` association; save after each)

**Task A — the Gaudin (∞,∞) family (both reps).**
- A1. `ΩDinf[ζ, ω1]` := e₁ᵀρ(exp(r₁E₁₂)exp(r₂E₁₃)exp(r₃E₂₃)) with r₁ = (ω₁−ζ₁)/(ω₁−ζ₂), r₂ = −(ω₁−ζ₁)(ζ₂−ζ₃)/((ζ₁−ζ₂)(ω₁−ζ₃)),
  r₃ = (ω₁−ζ₂)(ζ₁−ζ₃)/((ζ₁−ζ₂)(ω₁−ζ₃)). Check: equals `Limit[ΩD[ζ, {ω1, ω2, ω3}], ω2 → ∞]` symbolically (ω₃ drops out of ΩD already?
  — state what ΩD depends on) and numerically at ω₂ = 10⁶ to 10⁻⁶ relative. Is ΩDinf⊗ΩDinf a left eigencovector of the Laurent
  coefficients of S₂^{(g′/ω₂)}, S₃^{(g′/ω₂²)} in the limit (rescaled twist)? Compute the exact common left eigencovector of the two-site
  Laurent coefficients with the rescaled, limited twist if a finite limit exists; otherwise just record that the family is defined as a
  limit of eigencovectors.
- A2. Companion frame. `Ωc[ω]` := e₁ᵀρ(upper Gauss factor of `PauxG[ω]^{−T}`) (the one-site rule with M = PauxG[ω], whose rows are the
  left eigenvectors of gc[ω]); verify it is the exact left eigencovector of the one-site S₂, S₃ Laurent coefficients with twist gc[ω]
  (one-site cdet with L¹_ij = gc[ω]_ij + E_ij/(u−θ₁)); its Gauss parameters (α, β, γ) as rational functions of ω (Factor); the ω₂ → ∞
  limit `Ωcinf[ω1]`: is it e₁ᵀρ(e^{αE₁₂}e^{γE₂₃}) with β = 0 and α = γ (the XXX gaussCov{a,0,a} pattern)? Give α(ω₁). Frame tie:
  `Ωc[ω]` ∥ `ΩD[ζ,ω].PphysG[ζ]` (PphysG = rep image of (PauxG[ζ]⁻¹)ᵀ via the Gauss factorisation as in the round-1 `Pphys`): the scalar
  κ^G(ζ,ω) (compute as the ratio of first components; compare with the XXX kappaTie formula ∏_k Det[(n.P₃)_{1..k,1..k}]^{λ_k−λ_{k+1}} with
  n = nMat3[ζ,ω] and P₃ = PauxG⁻¹ᵀ or its inverse — report which, if either, matches) and its ω₂ → ∞ limit.
- A3. Exact undo for Φ∞ := ΩDinf[ζ,ω₁]⊗ΩDinf[ζ,ω₁′]: ⟨Φ∞|h⟩ = ⟨e₁ᵀρ(n∞(ω₁)n∞(ω₁′)⁻¹)⊗e₁ᵀ|h⟩ on all HW kets at P, P′ and the regime point
  (exact). Closed forms of the undone Gauss parameters (a′, c′, g′ = b′ − a′c′) in (φ, ω̂₁, ω̂₁′) (Factor; ε-free — assert); rank of the
  Jacobian of (ω̂₁, ω̂₁′) ↦ (a′, c′, g′) (expect 2); per-weight ranks of φ_u = e₁ᵀρ(n′) over the grid (CG-complete iff rank = multiplicity on
  every weight of V; adjoint has a weight of multiplicity 2, [3,1,0] of multiplicity 2 as well — if a weight falls short, enlarge the grid
  and say so; if it still falls short, the family is NOT CG-complete in Gaudin: report immediately in the result cell, §8 of the plan).
- A4. K/S for every HW branch at the three ε (expect O(ε), doubling), S = O(1); CG at leading order from the overlap alone vs the exact
  G-norm-1 HW-ket slices (Richardson), the 2-dim slices from pairs of grid points, octet mixing R_c vs the exact roots of R² − 5R − 5 = 0;
  controls: ω₁′ = ω₁ (expect a′ = g′ = 0, one irrep read) and ω₁′ = 7/4 fixed (not scaled; expect fine, as in round 1). Both reps.
- A5. Reference tables for the Python agents (markdown): at P with (ω₁, ω₁′) = (5/2, 7/4), for ALL states (G-normalised, sign convention
  of round 1): Qnum(B) = ⟨Φ∞|Ψ_B⟩/⟨x0G|Ψ_B⟩ and the two equal-twist overlaps ⟨ΩD[ζ,ω₀]⊗ΩD[ζ,ω₀]|Ψ_B⟩, ⟨Φ∞|Ψ_B⟩ themselves, 30 digits,
  with the (S₂(u₀), S₃(u₁)) eigenvalue pair of each state to identify it (sorted as in round 1); the same for the second (ω₁, ω₁′) and at
  P′; and at the regime point (ω̂₁, ω̂₁′) = (3/4, −12/5), ε = 10⁻⁶, for the HW branches.

**Task B — the [3,1,0] fourth localiser (Λ = {3,1,0} only; exact at P and P′).**
- B1. With D₁ = {id, W, H₁, S3c[1,2], S3c[1,1], S3c[2,2], S3c[2,1]} and the quadratic monomials a_i.a_j (i ≤ j, over the 6 non-identity
  elements → 21 products), build D₂ (rank 26). The localising subspace 𝓛 := {X ∈ D₂ : ⟨Ω⊗Ω|X ∈ V⊗Ω₁} exactly (NullSpace of the linear
  system Φ_X·k = 0 for k in the annihilator of Ω₁); dim 𝓛; dim(𝓛 ∩ D₁) (expect 3: id, H₁, S3c[1,2]); the quotient dimension (expect 1).
- B2. Sparse representatives of the extra direction: for every support of size 1, 2, 3 among the 21 quadratic monomials (optionally with
  free D₁ coefficients), test whether 𝓛 contains an element with that quadratic support outside D₁ (exact rank test; ≤ 1561 small
  systems — time one, then loop; the 64→225 dimension is 225×225 sparse, fine). Report the minimal supports found and the explicit X.
  Natural candidates to report individually: H₁², H₁.S3c[1,2], S3c[1,2]², S3c[1,1], H₁.S3c[1,1], S3c[1,2].S3c[1,1], S3c[1,1]², W.H₁,
  W.S3c[1,2], and their site-2 mirrors.
- B3. The fourth covector v₄ ∈ L₁ (choose the representative G-orthogonal to the span of Ω₁, Ω₁C_ζ, Ω₁S3^{1-site}[−2]): express it in
  Ω₁·span{one-site operators {1, C_ζ, S2^{1-site}, S3^{1-site} Laurent coefficients at θ₁} and their pair products} — print the
  decomposition with the minimal number of terms; its spectator factor as a Laurent polynomial in θ₁₂ (fit on 6 values, check 3); the
  weights of V reached by the 3-dim L₁ (ranks per weight) vs the 4-dim L₁ (round 1: the 4-dim one reaches every weight with full rank —
  does the 3-dim one fail on exactly one weight? which?). Same for site 2 (ℓ₂ = 4 over D₂ per round 1).
- B4. Compare with the XXX rule for [3,1,0] (T_(k,k) for k ≤ d = 2 and T_(k,d) for d < k ≤ 3 at θ₁ − 3h: shapes (1,1), (2,2), (3,2)):
  three fused shapes plus the identity = 4 = ℓ₁. State which Gaudin element plays the role of the third shape (3,2) — a guess supported by
  the degrees (a product of two Laurent coefficients, as the fused T_(3,2) is cubic in T) is enough; no XXX computation here.

**Task C — closed forms and universality (exact; symbolic ζ, ω, θ where the objects are small).**
- C1. c_H in ⟨Ω⊗Ω|H₁ = [Ω₁(c_H − C_ζ)]⊗Ω₁ (round 1: c_H = |Λ|Σζ + 4/θ₁₂ for [2,1,0], |Λ|Σζ + 6/θ₁₂ for [3,1,0] at P) as a closed form in
  (Λ, ζ, θ₁₂): compute it for Λ ∈ {{1,0,0}, {1,1,0}, {2,0,0}, {2,1,0}, {2,2,0}, {3,1,0}, {3,2,0}, {4,1,0}} (the ⟨Ω⊗Ω|E1_ikE2_ki part is a
  one-site computation on Ω₁E_ik ⊗ Ω₁E_ki; symbolic ζ, ω fine at these dims; use the exact Ω and symbolic θ₁₂) and fit the 1/θ₁₂
  coefficient as a polynomial in (λ₁, λ₂) (conjecture 2λ₁ or |Λ|²−c₂(Λ)-type — report the fit and the check on the two largest reps).
- C2. The spectators f(H₁), f(S3c[1,2]) (round 1: Laurent polynomials in θ₁₂ at fixed ζ, ω: f(H₁) = 6 + 4/θ₁₂ at P, 23/6 + 4/θ₁₂ at P′ for
  [2,1,0]) as closed forms in (ζ, ω₁, ω₂, θ₁₂) for [2,1,0] (symbolic; if the symbolic 64×64 is too slow, fit rational functions of low degree
  on a grid of rational (ζ, ω) points and verify on held-out points); then [3,1,0] (fit).
- C3. The vacuum eigenvalue functions: H_s^{vac}, a_s^{vac}, b_s^{vac} (Laurent data of S₂^{(g′)}, S₃^{(g′)} on Ω⊗Ω) in closed form in
  (ζ, ω, θ) for [2,1,0] (round 1 numbers at P, ω: H₁ = 12539/930, H₂ = 4139/930, a_s = 9/5, b₁ = 2173/310, b₂ = 493/310; e₂(ω) = 131/15,
  e₃(ω) = 7/2) — derive symbolically from the one-site statement (Ω is the HW covector of the Cartan element g′ ⇒ the eigenvalues are
  those of the HW vector of weight Λ for the diagonalised g′: S_k^{vac} should be expressible through ω and ζ by the dictionary of the
  Laurent coefficients; print and verify at P, P′, and the second ω).
- C4. Universality of the untwisted copies: the matrix of S3c[1,2] on the exact HW space of the multiplicity-2 sector, in the basis of
  round 1's `hwKets` (adjoint [3,2,1]: (105/62)[[3,−1],[−1,7]] at P with θ₁₂ = 31/35, i.e. (3/(2θ₁₂))[[3,−1],[−1,7]]?): recompute at P′,
  and at θ = (1/5, −4/9) and (3/7, −1/4) with ζ = P; report whether it equals (3/(2θ₁₂))·[[3,−1],[−1,7]] + s·1 with s scalar (ζ-dependent)
  — then the eigenvectors and the mixing roots of R² − 5R − 5 = 0 are universal; the same for [3,1,0] sectors (5,2,1) and (4,3,1) (for
  (4,3,1) give the mixing as the ratio of the two HW-ket coordinates of each copy in the GT basis rather than through e₁E₁₃, e₁E₂₃E₁₂).

**Task D — summary cell** `### sl(3) Gaudin L = 2 chain, part 2 (2026-10-16)`: A (family status: CG-complete yes/no, K/S, CG, the
companion vacuum's Gauss form and κ^G), B (the fourth localiser), C (closed forms, universality), surprises, what was not done.

## 2. Rules and gotchas

- New-notebook workaround: build/extend the `.wb` as JSON with a Python script (`python -X utf8`; cells `kind` 1 markdown / 2 code,
  `value`, `languageId` "wolfram"/"markdown", `outputs` [], `metadata` {}), evaluate each code cell's exact saved source in K1 via
  `wolfbook_evaluateExpression` (generous `timeoutSeconds`), store the printed `Grid` as a markdown cell under the code cell, save by
  rewriting the JSON (re-read before every rewrite). No native outputs — say so in the report. Kernel K1 is yours alone; the first cell
  wipes it. Do not evaluate or edit any other notebook.
- Semicolons between statements; never `*)` inside a comment; no symbols named D, E, N, I, C, K; memoised definitions cache bad results
  (re-run the definition and everything downstream after a fix); `===` on rationals, `Precision` on numerics, `SetPrecision` before
  `Norm` of near-zero vectors; `Max[{}]` = −∞.
- SparseArray for E1, E2, L_ij and products; dense `N[…, pr]` only per weight sector. [3,1,0]: 225×225 — time one product before loops.
- Report numbers with precision, not adjectives; if an expectation is refuted and two attempts do not fix it, stop that item, keep what
  works, report what you saw.
