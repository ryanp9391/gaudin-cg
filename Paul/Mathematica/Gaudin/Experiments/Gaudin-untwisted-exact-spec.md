# Spec W4: exact chain-side reference for the untwisting limit — HW decompositions and CG tables, untwisted Laurent data via Casimirs, the untwisted oper, perturbation theory in ε (2026-10-16)

Delegated to a wolfbook-builder subagent (SONNET); planned by the main session. **New notebook:** `Paul/Mathematica/Gaudin/Experiments/
Gaudin-untwisted-exact.wb` (self-contained; exact Wolfram algebra; first code cell `Quiet[ClearAll["Global`*"]]`; loops over the representations
Λ ∈ {{2,1,0}, {3,1,0}, {3,2,0}, {2,2,0}, {1,1,0}, {4,1,0}} — the adjoint first and in full, the others as the budget allows, [4,1,0] (dim 24, 576 states) last).
Read first: `../Gaudin-sl3-round4-plan.md` §0 (notation — binding), §2 "W4", §4, §6, and the round-1/2/3 plans' verdict blocks. **Templates:** the
round-1/2/3 notebooks `Gaudin-L2-chain.wb`, `Gaudin-L2-chain-2.wb`, `Gaudin-L2-chain-320.wb` (same directory; read the JSON with a `python -X utf8`
script): generators/metric, Lax, cdet, `S[k]`, Laurent coefficients `Hres`/`S3c`, exact `hwKets`/`PHWc`, the exact HW kets and sign convention, the
universality cells (S3c[1,2] on the HW space), `ΩD`/`ΩDinf`/`x0G`, the regime machinery. **Also write machine-readable data files** (via Bash,
`python -X utf8` or Wolfram `Export`) to `Paul/Mathematica/Gaudin/ref/untwisted-exact-<l1><l2><l3>.txt` (e.g. `untwisted-exact-210.txt`) for the Python
agents U1/U2: plain text, one record per line, Wolfram InputForm of exact rationals (no `Sqrt` — give unnormalised rational HW kets plus the exact
G-norm² separately), with a header describing the GT basis order (state the order explicitly: the notebook's `ps` order with the HW ket first) and
every field. Write the adjoint file as soon as E1–E2 are done for it (the Python agents read it during their run), then the others.

## 0. Conventions

θ = (θ₁, θ₂) SYMBOLIC where feasible (the Laurent data at ζ = 0 are rational functions of θ₁₂ = θ₁ − θ₂ by scaling — verify), numeric θ = (2/7, −3/5) where
not. Untwisted model: ζ = 0 (g = 0): S₁ = |Λ|Σ_s 1/(u−θ_s), S₂ = Σ_s[H_s/(u−θ_s) + c₂/(u−θ_s)²], S₃ = Σ_s[b_s/(u−θ_s) + a_s/(u−θ_s)²], c₂ = λ₂(λ₁+1), with
H₁ + H₂ = 0 at ζ = 0 (round 1: H₁ + H₂ = −W + 2|Λ|Σζ). Regime (E4): ζ = εφ with φ SYMBOLIC (φ₁, φ₂, φ₃), ω̂₁, ω̂₁′ symbolic; HW-limit branch of sector μ =
the eigenstate of the twisted family that tends to the HW vector Ψ₀ (G-norm 1, sign: first nonzero GT component positive; in multiplicity-2 sectors
the copies are the eigenvectors of S3c[1,2] on the HW space — θ-universal, round 2/3). x₀^G = x(ζ)⊗x(ζ), x(ζ) = e₁ᵀP_phys^{G−1} (round-1 construction;
for the adjoint x = (ζ₂−ζ₁, ζ₃−ζ₁, ζ₂−ζ₁, 2ζ₃−ζ₁−ζ₂, 3(ζ₂−ζ₁), 3(ζ₃−ζ₁), 2(ζ₃−ζ₂), 3(ζ₃−ζ₂)) in notebook order — verify, and give x(ζ) exactly for every Λ).
Φ∞ = ΩDinf[ζ,ω₁]⊗ΩDinf[ζ,ω₁′] (round-2 r_k^∞), exactly ε-free in the regime; φ_u = e₁ᵀρ(n′) with the round-2 closed forms (a′, c′, g′).

## 1. Tasks (code cell + markdown result cell per task and rep; `checks` association; save after each; write the data files as you go)

**E1 — exact HW decomposition and CG tables.** For each Λ: the weights of V⊗V, the HW sectors μ with multiplicities (solve E₁₂h = E₂₃h = 0 exactly on
each weight space — `NullSpace` over the rationals), exact HW kets (rational, with the sign convention) and their exact G-norm² (gMet), the universal
copies in multiplicity-2 sectors (eigenvectors of S3c[1,2]|_HW with θ symbolic — confirm θ-independence; they may involve square roots: give the
2×2 matrix and the quadratic exactly, plus the two copies as exact algebraic vectors), and the **CG tables C_ν = ⟨e_ν⊗e₁|Ψ₀⟩/‖Ψ₀‖_G** for every HW
state (site-1 GT label ν, site 2 = HW; exact, as rational×Sqrt[rational]). Write to the data file: sectors, multiplicities, HW kets (unnormalised
rationals), G-norm², copy matrices, CG tables (as unnormalised rationals + the norm), the GT basis order, the weights of the GT basis.

**E2 — untwisted Laurent data via Casimirs.** At ζ = 0 with θ symbolic: the exact values of H₁, H₂, a₁, a₂, b₁, b₂ on every HW state (operators `Hres`,
`S3c` at ζ = 0 applied to the exact HW kets; check they are eigenvectors; for multiplicity-2 sectors the 2×2 matrices and their eigenvalues on the
universal copies). The Casimirs: C₂ = Σ_ikE_ikE_ki and C₃ = Σ_ijkE_ijE_jkE_ki of the total gl(3) on each sector (exact; also C₂(Λ), C₃(Λ) per site). Closed
forms: verify H₁ = (|Λ|² − ½[C₂(μ) − 2C₂(Λ)])/θ₁₂ (derive the identity Σ_ikE¹_ikE²_ki = ½(C₂^{tot} − C₂⁽¹⁾ − C₂⁽²⁾) symbolically) and FIND the closed forms
of a_s and b_s in terms of (Λ, μ, θ₁, θ₂) — fit over all sectors of all reps with an ansatz linear in {1, C₂(μ), C₃(μ), C₂(Λ), C₃(Λ), |Λ|, products} times
{1/θ₁₂, 1/θ₁₂²} (scaling fixes the θ-powers: a_s ∝ 1/θ₁₂? b_s has 1/θ₁₂ and 1/θ₁₂² pieces? — determine); report the sum rules b₁ + b₂ and
a₁ + a₂ + θ₁b₁ + θ₂b₂ per sector (do they vanish? they decide whether the untwisted oper is Fuchsian at ∞). Write the Laurent data and Casimirs to the
data files (per sector, exact rationals at θ = (2/7,−3/5) AND the symbolic forms as strings).

**E3 — the untwisted oper per sector.** D₀ = ∂³ − S₁∂² + S₂∂ − S₃ with the E2 data (θ symbolic). (a) Local exponents at θ₁, θ₂ (expect {0, λ₂+1, λ₁+2}) and at
u = ∞ (substitute u = 1/x; is ∞ a regular singular point? exponents; the Fuchs relation check). (b) Polynomial solutions: for degrees d ≤ 2λ₁ + 4, the
exact solution space of D₀p = 0 among polynomials (dimension per degree; the degrees that occur); compare with the dual degrees μ_j of the sector and
with the direct degrees 2λ₁ − μ_i (the untwisted direct oper = Π^{λ₁−λ₃}∘adjoint∘Π^{−(λ₁−λ₃)} of D₀ — its polynomial solutions too). (c) Hypergeometric
identification: with x = (u−θ₁)/(u−θ₂) (or (u−θ₁)/(θ₂−θ₁)) transform D₀ to a Fuchsian equation on {0, 1, ∞}; if it has three regular singular points it is
a ₃F₂ (generalised hypergeometric) equation up to a gauge factor x^α(1−x)^β — identify the parameters (`DifferentialRootReduce`, or match the
Riemann-scheme exponents), and express the polynomial solutions as `HypergeometricPFQ` (Jacobi-type polynomials); verify by `FullSimplify`. If ∞ is
irregular, say so and characterise the Stokes/rank structure; then instead identify the solutions of D₀ among {polynomials} ∪ {polynomial × log}
(do logs appear? the Frobenius exponents at θ_s differ by integers).

**E4 — degenerate perturbation theory in ε (adjoint; then [3,1,0] if time).** With ζ = εφ symbolic: the HW-limit branches Ψ_B(ε) = Ψ₀ + εΨ₁ + ε²Ψ₂ + … and
their Laurent data H_s(ε), a_s(ε), b_s(ε) to O(ε²) exactly (Rayleigh–Schrödinger on the full 64-dim space, exact rational linear algebra; the unperturbed
family is degenerate within each gl(3) multiplet, so use the exact projector onto the μ-isotypic component and the S3c[1,2] copies to fix the zeroth-order
states; the perturbation is the ζ-dependence of S₂, S₃ — write S_k(ε) = S_k^{(0)} + εS_k^{(1)} + ε²S_k^{(2)} exactly). Report H_s^{(1)}, a_s^{(1)}, b_s^{(1)} as
exact linear forms in φ (and second order if feasible) per branch; check against the regime numerics of round 1 (K/S, branch eigenvalues at ε = 10⁻⁶).
Also the O(ε) correction Ψ₁ (exact vector in φ) per branch and its weight content. Write to the data file.

**E5 — exact leading orders in the regime.** With φ, ω̂₁, ω̂₁′ symbolic and ε symbolic: ⟨x₀^G|Ψ_B(ε)⟩ to the first nonvanishing order (x(εφ) = εx(φ); which
order of Ψ_B first gives a nonzero overlap with x(φ)⊗x(φ)? give the exponent 2 + k_B and the coefficient), ⟨Φ∞|Ψ_B(ε)⟩ = ⟨φ_u⊗e₁|Ψ₀⟩ + ε(…) (the O(1)
term as the exact CG generating function in (ω̂₁, ω̂₁′) — print ĉ·⟨e₁ᵀρ(n′)⊗e₁ᵀ|Ψ₀⟩ expanded as a polynomial in ω̂₁, ω̂₁′ with the C_ν as coefficients, per
branch), and c = ε⁸ĉ. Write ĉ·(generating function) per branch to the data file (adjoint).

**E6 — summary cell** `### Untwisted sl(3) Gaudin: exact reference (2026-10-16)`: per rep the sectors/multiplicities, the closed forms of the untwisted
Laurent data, the oper's nature at ∞ and the hypergeometric verdict, the perturbative first orders, the list of data files written and their formats.

## 2. Rules and gotchas

New-notebook workaround (JSON built by a Python script, `code <path>` once, every code cell's exact saved source evaluated in K1 via
`wolfbook_evaluateExpression`, printed `Grid` stored as a markdown cell under the code cell, re-read before every rewrite, no native outputs — say so).
K1 is yours alone (kernel `k-d8a6a147`, client `VSCode[gaudin-cg]`); the first cell wipes Global`. Semicolons; never `*)` in a comment; no symbols D, E, N, I,
C, K; memoised definitions cache bad results; exact algebra only here (no N[] except for cross-checks); symbolic θ, φ can make 64×64 operators slow —
keep symbols at the level of eigenvalues/Laurent data and HW kets (small), evaluate operators at numeric θ when needed and recover the θ-dependence by
scaling (degree in 1/θ₁₂) or by fitting at 3 θ values; time one case before any loop; nothing over ~3 min without a timed smaller run; [4,1,0] only if
the rest is done. Only this notebook and the `ref/untwisted-exact-*.txt` data files; no other notebooks, `paper/*.tex`, HANDOFF.md, CLAUDE.md,
`Paul/Python/`, or git. Numbers with precision, not adjectives; if an expectation is refuted and two attempts do not fix it, stop that item, keep what
works, report.
