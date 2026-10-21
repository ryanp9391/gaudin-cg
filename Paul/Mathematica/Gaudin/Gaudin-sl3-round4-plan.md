# sl(3) Gaudin, round 4: the general kernel and the analytic untwisting limit — CG coefficients of gl(3) from Gaudin SoV (2026-10-16)

Status: DONE 2026-10-16 — all four SONNET agents (W4, G1, U1, U2) returned and agree; result blocks below §2 (U2, U1, G1, W4). Results: `Experiments/Gaudin-untwisted-exact.wb` + `ref/untwisted-exact-*.txt`, `Paul/Python/gaudin_sl3/RESULTS-general-kernel.md`, `RESULTS-untwisted-Q.md`, `RESULTS-untwist-limit.md`. Next steps in HANDOFF.md.
**Project target (Paul, 2026-10-16): Clebsch–Gordan coefficients of gl(3) from separation of variables in the Gaudin model.** gl(4)
is no longer a target; XXX was the scaffold. Rounds 1–3 (`Gaudin-sl3-CG-plan.md`, `Gaudin-sl3-round2-plan.md`, `Gaudin-sl3-round3-plan.md`)
established everything below. The two questions of this round: **(A)** guess the general-Λ formula for the two-vacua overlap kernel from
the three computed representations; **(B)** take the untwisting limit analytically to get ⟨v ⊗ hws|Ψ₀⟩ in terms of Q-functions, with v a
rotated covector on site 1, hws = e₁ on site 2, Ψ₀ the untwisted HW eigenstate — i.e. the CG coefficients as a generating function in the
rotation parameters, from the Q-functions alone.

## 0. What is established (two-site gl(3) Gaudin, Λ = (λ₁, λ₂, 0) on both sites, θ = (θ₁, θ₂), twist ζ)

**The two-vacua overlap is a Q-number (the jet formula; verified three independent ways for [2,1,0], by fit for [3,1,0], [3,2,0]).**
Family Φ∞(ω₁, ω₁′) = Ω^∞[ζ,ω₁] ⊗ Ω^∞[ζ,ω₁′], Ω^∞[ζ,ω] = e₁ᵀρ(exp(r₁E₁₂)exp(r₂E₁₃)exp(r₃E₂₃)), r₁ = (ω−ζ₁)/(ω−ζ₂), r₂ = −(ω−ζ₁)(ζ₂−ζ₃)/((ζ₁−ζ₂)(ω−ζ₃)),
r₃ = (ω−ζ₂)(ζ₁−ζ₃)/((ζ₁−ζ₂)(ω−ζ₃)) (ω₂ → ∞ vacua; CG-complete iff ω₁ ≠ ω₁′). For every eigenstate Ψ with dual Q-functions P^j = e^{ζ_j u}p_j
(monic, deg p_j = μ_j), jets at θ_s, Plücker coordinates p_kl(θ_s) = P²^{(k)}P³^{(l)} − P³^{(k)}P²^{(l)}:

    c · ⟨Φ∞|Ψ⟩/⟨x₀^G|Ψ⟩ · M(θ₁)M(θ₂) = Σ_{(kl),(k′l′)} K_{(kl),(k′l′)}(ω₁, ω₁′, θ₁₂) · (P³p_kl)(θ₁) · (P³p_k′l′)(θ₂)

- (kl) runs over the NON-FORCED pairs 0 ≤ k < l ≤ λ₁+1 with not both k, l ≤ λ₂ (p_kl ≡ 0 for k < l ≤ λ₂ because P^{(k)}(θ_s) = κ_{k,s}P(θ_s)
  for k ≤ λ₂ — dual local exponents {0, λ₂+1, λ₁+2}); count C(λ₁+2, 2) − C(λ₂+1, 2) = 5, 9, 7 for [2,1,0], [3,1,0], [3,2,0]. The lone factor is
  always P³ (P³^{(c)}p_kl is redundant).
- M = P³·p_{0,λ₂+1} (the first non-forced Plücker coordinate; it is ∝ Π′(θ_s)Q₁(θ_s)·(ζ₃−ζ₂)e^{Σζθ_s} by the Wronskian relation).
- x₀^G = x(ζ)⊗x(ζ), x(ζ) = e₁ᵀP_phys^{G −1}, P_phys^G = ρ((P_aux(ζ)⁻¹)ᵀ); for the adjoint x(ζ) = (ζ₂−ζ₁, ζ₃−ζ₁, ζ₂−ζ₁, 2ζ₃−ζ₁−ζ₂, 3(ζ₂−ζ₁), 3(ζ₃−ζ₁),
  2(ζ₃−ζ₂), 3(ζ₃−ζ₂)) in notebook order (HW first); it is homogeneous of degree 1 in ζ.
- c = (ζ₂−ω₁)^{λ₁−λ₂}(ζ₃−ω₁)^{λ₁}(ζ₂−ω₁′)^{λ₁−λ₂}(ζ₃−ω₁′)^{λ₁}(ζ₁−ζ₂)^{2λ₂} = (ω₁ω₁′)^{2λ₁−λ₂}/(κ^G(ζ,ω₁)κ^G(ζ,ω₁′)), the frame-tie scalar
  (κ^G = Δ₁^{λ₁−λ₂}Δ₂^{λ₂}, Δ₁ = ω²/((ζ₂−ω)(ζ₃−ω)), Δ₂ = ω/((ζ₁−ζ₂)(ζ₃−ω)) at ω₂ = ∞). ζ enters ONLY through c and x₀^G.
- K is unique in this basis, rational, a polynomial in (ω₁, ω₁′) (degree 8 − e_i − e_j for [3,2,0] with e = (0,1,1,2,2,3,4); 10 − w_i − w_j for the
  adjoint with weights (2,3,3,4,5)) and a polynomial in 1/θ₁₂ (degree 4 adjoint, 7 for [3,2,0]); depends on θ only via θ₁₂; site exchange acts as
  K_ij(ω,ω,1/θ₁₂) = K_ji(ω,ω,−1/θ₁₂). All entries known: adjoint (`RESULTS-twovac-limit.md`, `twovac_fit.d2_G`), [3,2,0] (`twovac_fit_320_grid.pkl`,
  base point in `RESULTS-Qside-320.md`), [3,1,0] (base point only).
- **Canonical kernel.** In the Gaudin companion frame (Ψ_c = (P_phys^{G−1})^{⊗2}Ψ, x₀^c = e₁⊗e₁, Φ_c∞ = Ωc∞[ω₁]⊗Ωc∞[ω₁′] with
  Ωc∞[ω] = e₁ᵀρ(e^{aE₁₂}e^{aE₂₃}), a = 1/ω): ⟨Φ_c∞|Ψ_c⟩/⟨x₀^c|Ψ_c⟩·M₁M₂ = Σ Θ_ij m_i(θ₁)m_j(θ₂), m_i = P³p_kl, with **Θ = K/(ω₁ω₁′)^{2λ₁−λ₂}**,
  a POLYNOMIAL in (a, a′) = (1/ω₁, 1/ω₁′) (D3: Θ = G/(4ω³ω′³) for the adjoint, D2's "4c" convention). Θ is the Gaudin K tensor: Φ_c∞ expanded in
  the canonical jet basis, whose covectors ⟨i,j| satisfy ⟨i,j|Ψ_c⟩/⟨x₀^c|Ψ_c⟩ = m_i(θ₁)m_j(θ₂)/(M₁M₂) with the eight one-site monomials
  m = (P³p₀₂, P³′p₀₂, P³p₀₃, P³′p₀₃, P³′p₁₂, P³p₂₃, P³′p₁₃, P³′p₂₃) (adjoint; weights 2,3,3,4,4,5,5,6; P³′p_kl = P³p_{k+1,l}-type by Frobenius); the jet
  matrix Y is exactly ζ-independent and rational, 59 of its 64 covectors couple the sites through 1/θ₁₂. **At 1/θ₁₂ = 0, Θ has rank 1 (the
  sites decouple): Θ → θ(a)⊗θ(a′), a ONE-SITE kernel.**
- Frobenius data: κ_s = P′/P = a_s/c₂ (c₂ = λ₂(λ₁+1) = S₂ double pole), κ₂,_s = P″/P = (a² − Ha + c₂b)/24 for [3,2,0]; direct oper =
  Π^{λ₁−λ₃}∘adjoint∘Π^{−(λ₁−λ₃)}; Wronskians: P^aP^b′ − P^bP^a′ = (ζ_b−ζ_a)e^{Σζu}Π^{λ₂}·Q_c (dual), Q_aQ_b′ − Q_bQ_a′ = (ζ_a−ζ_b)e^{−Σζu}Π^{λ₁−λ₂}P^c
  (direct); W = det[P,P′,P″] ∝ e^{Σζu}Π^{|Λ|}; the Laurent data satisfy three no-log conditions per site (Gaudin "Bethe equations").
- **Chain side of the limit (exact, rounds 1–3).** In the regime ζ = εφ, ω₁ = εω̂₁, ω₁′ = εω̂₁′ the vacua are exactly ε-free and the undo is
  exact: ⟨Φ∞|h⟩ = ⟨φ_u ⊗ e₁ᵀ|h⟩ on every HW ket h, with φ_u = e₁ᵀρ(n′), n′ = exp(a′E₁₂)exp(g′E₁₃)exp(c′E₂₃), a′ = (φ₁−φ₂)(ω̂₁−ω̂₁′)/((φ₂−ω̂₁)(φ₂−ω̂₁′)),
  c′ = (φ₁−φ₃)(φ₂−φ₃)(ω̂₁−ω̂₁′)/((φ₁−φ₂)(φ₃−ω̂₁)(φ₃−ω̂₁′)), g′ = (φ₃−φ₁)(φ₃−φ₂)(ω̂₁−ω̂₁′)/((φ₃−ω̂₁)(φ₂−ω̂₁′)(φ₃−ω̂₁′)) (ε-free). For every HW-limit branch
  Ψ_B = Ψ₀ + εΨ₁ + … (Ψ₀ an exact HW vector of weight μ, G-norm 1): ⟨Φ∞|Ψ_B⟩ = ⟨φ_u⊗e₁ᵀ|Ψ₀⟩·(1 + O(ε)) — **this is the CG generating function**:
  ⟨φ_u(ω̂₁,ω̂₁′)⊗e₁ᵀ|Ψ₀⟩ = Σ_ν [e₁ᵀρ(n′)]_ν C_ν, C_ν = ⟨e_ν⊗e₁|Ψ₀⟩ the CG coefficients in the GT basis. No cleaning is needed on any branch of
  any computed rep. Numerically (D1 F4, adjoint): the Q side reproduces all CG data; the total Σ K N N is O(ε^{N_B}) with N_B = 8,7,7,6,6,4 for
  (4,2,0), (4,1,1), (3,3,0), (3,2,1)×2, (2,2,2) while the individual terms are O(ε²…ε³) — the limit lives in a cancellation of depth up to 6 orders.

## 1. Goal A — the general-Λ kernel (agent G1)

Route: Θ is polynomial in (a, a′, 1/θ₁₂); at 1/θ₁₂ = 0 it is θ(a)⊗θ(a′) with θ the ONE-SITE kernel, which is exactly computable for ANY Λ
(one site: the eigenstates are the GT weight vectors, the Q's are exact rationals, dim V ≤ 24); the 1/θ₁₂ corrections encode the Gaudin
K tensor's site coupling. Steps: (K0) canonical normalisation Θ for the three known reps (polynomiality in a, a′; torus homogeneity per
entry); (K1) the one-site kernel θ_i(a; Λ) exactly for Λ up to [4,2,0], then its closed form in (λ₁, λ₂, k, l) — the conjecture's seed; (K2)
Θ for more reps with the existing fit pipeline: [1,1,0] (dim 3, 2 monomials/site), [2,2,0] (dim 6, 3), [3,1,0] polynomial laws (missing), [4,1,0]
(dim 24, 14 monomials/site, budgeted), and the λ₂ = 0 probes [1,0,0], [2,0,0] (the S₂ double pole is not scalar there — report whether the
structure survives); (K3) structure: Θ(1/θ₁₂ = 0) = θ⊗θ′ verified, the 1/θ₁₂-expansion Θ = Σ_n θ₁₂^{−n}Θ_n (ranks, supports, site-exchange rule
Θ_ij(a,a′,t) = Θ_ji(a′,a,−t)), the **top-state sum rule** Σ Θ_ij m_i^{top}m_j^{top} = M^{top}(θ₁)M^{top}(θ₂) with the top state's Q's in closed
form (P³ = e^{ζ₃u}, P² = e^{ζ₂u}p₂, p₂ the polynomial solution of p₂′ − (ζ₃−ζ₂)p₂ = −(ζ₃−ζ₂)Π^{λ₁−λ₂}·(normalisation) from the Wronskian relation
with Q₁ = e^{−ζ₁u}: p₂ = Σ_k (ζ₃−ζ₂)^{−k}(Π^{λ₁−λ₂})^{(k)} up to sign/normalisation — an exact polynomial identity in (ζ, a, a′, θ) that every
candidate kernel must satisfy); (K4) the conjecture for general Λ (closed form, or a determinantal/factorised structure, or at least the rule
set: degrees, weights, zero pattern, θ⊗θ′ and Θ₁ in closed form), with its verification status per rep.

## 2. Goal B — the analytic untwisting limit (agents U1, U2, W4)

**The difficulty.** At ζ = 0 the three dual Q's e^{ζ_j u}p_j lose their exponentials and their Wronskians vanish ((ζ_b − ζ_a) factors), so the
untwisted limits p_j^{(0)} of the branch Q's are linearly dependent: the ε → 0 limit of the solution space is a ζ-jet (p^{(0)} and ε-derivatives
p^{(1)}, …), and the Fuchs relation shows the ζ = 0 oper is not simply "three polynomials of degrees μ_j". Meanwhile the jet formula's total is
ε^{N_B} smaller than its terms, so the limit involves the ε-expansion of the Plücker monomials to depth up to ~6. The limit is therefore a
statement about the **ε-jet of the Q-functions at the untwisted point**, and the task splits into: the untwisted Q-system and the ε-expansion
of the branch Q's (U1, with W4's exact chain data), and the leading-order extraction of the formula (U2).

**U1 — the untwisted Q-system and the ε-expansion.** (T0) closed-form Laurent data of the untwisted HW states via Casimirs: at ζ = 0, H₂ = −H₁,
H₁ = (|Λ|² − ½[C₂(μ) − 2C₂(Λ)])/θ₁₂ (from Σ_ik E¹_ikE²_ki = ½(C₂^{tot} − C₂⁽¹⁾ − C₂⁽²⁾); C₂(λ) = Σ_i λ_i(λ_i + 4 − 2i)), and a_s, b_s via the cubic Casimir
(derive or fit exactly, verify against W4's tables); the sum rules at ∞ (b₁ + b₂, a₁ + a₂ + θ₁b₁ + θ₂b₂). (T1) the untwisted oper per HW sector:
singularities and exponents (at θ_s {0, λ₂+1, λ₁+2}; at ∞: compute — regular or not?), its polynomial solutions (degrees, dimension), the ε → 0
limits of the branch Q's and of their ε-differences, hence a basis of the limiting 3-dim solution space in closed form (hypergeometric/Jacobi
if Fuchsian with three singular points; W4 attempts the identification in Wolfram). (T2) the ε-expansion P^j(u;ε) = e^{εφ_j u}Σ_n ε^n p_j^{(n)}(u)
of the six adjoint branches: numerically to n ≤ 8 (series fits over ε_k = 10⁻⁶·2^{−k}, k = 0..11, 200 digits), and exactly for n ≤ 2 by
perturbation theory (the branch's Laurent data to O(ε²) — from the chain's degenerate perturbation theory, W4 — then the oper order by order).
(T3) the ε-expansion of the jet monomials N_i(θ_s) = P³p_kl to the orders U2 needs, per branch.

**U2 — the limit of the jet formula.** (L1) ε-order bookkeeping: c = ε^{4λ₁}ĉ(φ,ω̂); K_ij(εω̂₁, εω̂₁′) = ε^{d_ij}K_ij(ω̂₁,ω̂₁′) (homogeneity, d_ij the
total degree); N_i(ε) ~ ε^{n_i(B)}; M₁M₂ ~ ε²; ⟨x₀^G|Ψ_B⟩ ~ ε^{2+k_B} — determine k_B per branch and interpret (⟨x(φ)⊗x(φ)|Ψ₀⟩ vanishes for
non-top μ?); N_B = total exponent. (L2) the normalisation-free limit: ρ_B(ω̂;ω̂⁰) = [ĉ(ω̂⁰)/ĉ(ω̂)]·[Σ K(εω̂)N N]/[Σ K(εω̂⁰)N N] →
⟨φ_u(ω̂)⊗e₁|Ψ₀⟩/⟨φ_u(ω̂⁰)⊗e₁|Ψ₀⟩; the leading polynomial L_B(ω̂₁,ω̂₁′) := Σ_{ij}Σ_{n+m = N_B − d_ij} K_ij(ω̂)N_i^{(n)}N_j^{(m)} from the numerically
extracted N^{(n)} (200 digits), and the identity L_B(ω̂) = κ_B·ĉ(ω̂)·⟨e₁ᵀρ(n′(ω̂))⊗e₁ᵀ|Ψ₀⟩ as polynomials in (ω̂₁, ω̂₁′) (ĉ clears the denominators
of a′, c′, g′) — verified coefficient by coefficient against the exact CG; the depth (which n enter) per branch. (L3) with U1's p^{(n)}:
the explicit formula for the CG generating function of each adjoint branch in terms of the untwisted Q's and their ε-derivatives; solve for
the C_ν; check φ-independence of the C_ν (the generating function depends on φ only through the Gauss parametrisation); then [3,1,0]/[3,2,0]
if time (P4's kernels and regime machinery).

**W4 — exact chain-side reference (Wolfram).** For Λ ∈ {[2,1,0], [3,1,0], [3,2,0], [2,2,0], [1,1,0], [4,1,0]}: (E1) exact HW decomposition of
V_Λ⊗V_Λ (sectors, multiplicities, HW kets as exact rational vectors with their G-norms, the universal copies in multiplicity-2 sectors) and
exact CG tables C_ν = ⟨e_ν⊗e₁|Ψ₀⟩, written to `Paul/Mathematica/Gaudin/ref/untwisted-exact-<Λ>.txt` (machine-readable InputForm) for the
Python agents; (E2) exact Laurent data (H₁, a₁, b₁, H₂, a₂, b₂) at ζ = 0 per HW state with θ₁, θ₂ symbolic, the Casimirs C₂(μ), C₃(μ), and the
closed forms (derive the a_s, b_s law); (E3) the untwisted oper per sector symbolically: exponents at θ₁, θ₂, ∞, polynomial solutions, and
hypergeometric identification (`DifferentialRootReduce`, ₃F₂ matching); (E4) degenerate perturbation theory in ε with symbolic φ for the
HW-limit branches: Laurent data and Ψ_B to O(ε²) exactly (adjoint; the copies in (3,2,1) are fixed by universality); (E5) the exact leading
ε-orders of ⟨x₀^G|Ψ_B⟩ and ⟨Φ∞|Ψ_B⟩ with symbolic φ, ω̂; (E6) summary.

**Result from U2 (returned 2026-10-16, `RESULTS-untwist-limit.md`) — Goal B achieved for the adjoint:** the untwisting limit holds as a polynomial identity in
(ω̂₁, ω̂₁′) for all six branches, κ_B = 1, φ-independent CG coefficients recovered to 1e-77. Structural simplification: with exponents {0,2,4} every jet monomial is
M times a ratio expressible in the Laurent data — ν = (1, κ, π, κπ, ε_s), κ = a_s/3, π = 2H_s − a_s − 2Σζ − 6/δ, ε_s = −(2/3)[H_ta_s/δ − a_sb_s + a_se₂ + 3a_s/δ² − 3a_t/δ²
− 3b_t/δ − 3e₃], b_s = a_s(H_s − a_s)/3 — so **c·⟨Φ∞|Ψ⟩/⟨x₀^G|Ψ⟩ = Σ G_ij ν_i(θ₁)ν_j(θ₂) depends only on the eigenvalue data (H₁, a₁, H₂, a₂)**; the limit is
Σ G(εω̂)νν = ε^{6−k_B}ĉ⟨φ_u⊗e₁|Ψ₀⟩/X_B + O(ε^{7−k_B}), k_B = 0,1,1,2,2,4, needing the ν-series to order 3/5/5/4/4/2. The CG formula therefore needs: G, the ε-series of
(H_s, a_s) per branch (perturbation theory of the Gaudin Hamiltonians / of the no-log equations), the Gauss parameters and ĉ, one constant. U1 and W4 redirected to
the exact ε-series of the Laurent data.

**Result from U1 (returned 2026-10-16, `RESULTS-untwisted-Q.md`):** the Laurent data of the HW-limit branches have EXACT ε-series, polynomial in φ (degree n at
order εⁿ); top: H₁ = 4 + ε(φ₁+2φ₂+3φ₃), a₁ = a₂ = 3εφ₃ exactly; (4,1,1) and (3,3,0): closed algebraic forms with √(1 + ε²(φ_i−φ_j)²/4) (a 2×2 mixing); singlet to n = 6,
octets to n = 5 over Q(√5). Combined with U2's ν-formulas this makes the adjoint CG formula fully explicit (G, the ν-series from these Laurent series, the Gauss
parameters, ĉ). Untwisted facts (five reps, exact): the ζ = 0 oper is Fuchsian at ∞ (H₁+H₂ = 0, b₁+b₂ = 0, a₁+a₂+θ₁b₁+θ₂b₂ = 0 as operator identities); H₁ from C₂,
(a₁−a₂)/2 = (μ₁+2)(μ₂+1)μ₃/(2θ₁₂) (half the product of the Weyl-shifted weight, a C₃ quantity), b₁ = ±√β/θ₁₂² only in multiplicity-2 sectors; **all three untwisted Q's
are polynomials of degrees (μ₃, μ₂+1, μ₁+2)** (plan §2 expected fewer — corrected), not ₃F₂ in general; the monic twisted p_j have poles in ε, the finite objects are the
Plücker monomials and the Laurent data.

**Result from G1 (returned 2026-10-16, `RESULTS-general-kernel.md`) — Goal A, status:** (i) the ONE-SITE kernel in closed form for all Λ: θ_kl(a) = c_kl a^{k+l−λ₂−1},
c_kl = (−1)^{k+l−λ₂−1}C(λ₁−λ₂, l−λ₂−1)C(λ₁+1, k)(l−k)/((λ₂+1)C(l, λ₂+1)), equivalently a dressed Wronskian Σθ_klm_kl ∝ P³·W[(1−a∂)^{λ₁}P², (1−a∂)^{λ₁}P³](θ) (17 reps exact);
(ii) the BLOCK THEOREM: Θ^{[λ₁,λ₂,0]} = [λ₂!(λ₁−λ₂)!/λ₁!]²(aa′)^{−λ₂}·Θ^{[λ₁,0,0]} restricted to the non-forced pairs — one master kernel per λ₁, Θ = θ_i(a)θ_j(a′)R^{λ₁}_ij(at, a′t);
(iii) the master kernel in closed form to first order in t for all pairs and all λ₁, and to ALL orders for the single-index reps [l,l,0] (Γ-function form; l = 5 predicted and
confirmed on 441 states); (iv) the top-state sum rule corrected to Σ Θ m^{top}m^{top} = M M′(1−ζ₂a)^{λ₁−λ₂}(1−ζ₃a)^{λ₁}(same in a′); (v) λ₂ = 0 reps are regular (c₂ = 0) and
are the master kernels. Open: the general four-index master kernel beyond first order.

**Result from W4 (returned 2026-10-16, `Experiments/Gaudin-untwisted-exact.wb`, `ref/untwisted-exact-*.txt`):** exact reference for six reps (45 HW branches): CG tables, the
untwisted Laurent data with the closed form a₁ = −a₂ = [4T(1+T)² − 3(1+T)C₂(μ) + C₃(μ)]/(6θ₁₂) (= (μ₁+2)(μ₂+1)μ₃/(2θ₁₂)), b₁ = ±√disc/θ₁₂² only in multiplicity-2 sectors with
tr = eta + 2c₂, disc = eta² − 4c₂² except [3,2,0] (5,3,2); the untwisted oper is Fuchsian at ∞ with exponents {μ₁+2, μ₂+1, μ₃} on all 45 branches and is a ₃F₂ equation iff some
singular point has a consecutive exponent pair (20/45); exact Rayleigh–Schrödinger series of the adjoint branches to order 7 (agreeing with U1), Ψ₁ irrep content, k_B = 0,1,1,2,2,4
with closed-form leading coefficients ⟨x(φ)⊗x(φ)|Ψ_{k_B}⟩ (products of (φ_i−φ_j)), x(ζ) a homogeneous covector of degree λ₂, and the exact CG generating functions.

## 3. Agents, files, kernel

| agent | type | files (owner) | kernel |
|---|---|---|---|
| **W4** | wolfbook-builder (SONNET) | `Experiments/Gaudin-untwisted-exact.wb` + spec `Gaudin-untwisted-exact-spec.md`; data files `ref/untwisted-exact-*.txt` | K1 (only user) |
| **G1** | general-purpose (SONNET), Python | `Paul/Python/gaudin_sl3/general_kernel.py`, `RESULTS-general-kernel.md`, caches `general_kernel_*.pkl`; spec `SPEC-general-kernel.md` | none |
| **U1** | general-purpose (SONNET), Python | `untwisted_Q.py`, `RESULTS-untwisted-Q.md`, caches `untwisted_Q_*.pkl`; spec `SPEC-untwisted-Q.md` | none |
| **U2** | general-purpose (SONNET), Python | `untwist_limit.py`, `RESULTS-untwist-limit.md`, caches `untwist_limit_*.pkl`; spec `SPEC-untwist-limit.md` | none |

All earlier Python files are read-only (import or copy). Nobody edits `paper/*.tex`, HANDOFF.md, CLAUDE.md, other notebooks, or git.
Cross-agent reads: U2 reads `RESULTS-untwisted-Q.md`/its pkl (U1) and `ref/untwisted-exact-*.txt` (W4) when present; U1 reads W4's files;
G1 is independent. Run heavy stages separately (round 3: a single all-stage run stalled under CPU contention); ≤ 4 Python processes at once.

## 4. Numeric points

θ = (2/7, −3/5); ζ: P = (3/2, 5/4, 2/3), P′ = (2, 5/3, −1/2), extra (7/5, −3/4, 1/6), (−4/3, 2/9, 5/7); (ω₁, ω₁′) = (5/2, 7/4) base, grids as in
rounds 2–3; θ variants (1/5, −4/9), (3/7, −1/4), (1/2, −1/3). Regime: φ = (1/3, −2/7, 5/11), ω̂₁ ∈ {3/4, −5/3, 1/2, −7/3}, ω̂₁′ ∈ {−12/5, 7/2, 1/8, 5/4},
ε_k = 10⁻⁶·2^{−k}, k = 0..11, 200 digits (second direction φ′ = (−2/5, 3/8, 7/11) for the φ-independence check).

## 5. After this round

The CG formula written out (paper section of `paper/gaudinCG.tex`: conventions, the jet formula, the jet basis, the untwisting limit);
general-Λ proof of the kernel conjecture via the K tensor; [4,2,0] (multiplicity 3) as the last cleaning check; Clean promotion + docs agent.

## 6. Risks and fallbacks

- **The untwisted oper is not Fuchsian at ∞ / no hypergeometric closed form.** Then the "untwisted Q's" are defined by the ε-jet (p^{(0)}, p^{(1)}, …)
  computed perturbatively; U1 delivers them as exact rational polynomials (adjoint) and the formula is written in those.
- **The leading order needs very deep jets (n ≳ 6) and the numerics cannot separate orders.** Use more ε values and 300 digits, or the exact
  perturbative p^{(n)} from W4's exact Laurent-data series; report the depth per branch even if the closed form is out of reach.
- **The one-site kernel has no simple closed form in Λ.** Deliver the exact tables for Λ ≤ [4,2,0] and the rule set; the two-site conjecture may
  still be stated relative to the one-site kernel.
- **[4,1,0] too heavy (576 states).** Drop to the base point only, or skip; [2,2,0] and [1,1,0] already add two data points in λ₂.
