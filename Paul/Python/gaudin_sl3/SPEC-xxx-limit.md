# Spec P2: the XXX → Gaudin limit in Python — frames, Hamiltonian dictionary, localised covectors, Baxter → oper (2026-10-15)

Delegated to a general-purpose subagent (sonnet), Python only (the Wolfram kernel is taken by another agent). Read first:
`Paul/Mathematica/Gaudin/Gaudin-sl3-CG-plan.md` §0, §1 (binding conventions; the dictionary table is what you are verifying),
§4, §7. References (UTF-8 cell dumps in `Paul/Mathematica/Gaudin/ref/`): `ref-L2-two-vacua-untwist-all.txt` cells 3, 5, 7, 11
(XXX generators, exact transfer matrices via quantum minors, `rGauss`/`ΩD`/`Pphys`/`x0Of`, `Paux`/`Gc`, the eigencovector check),
`ref-Companion-B-C-diagonalize-clean-all.txt` (L = 1: Baxter equations rows 7–8, the one-site t_a, QQ), `ref-Minimal-ingredients-EH-all.txt`
(the four localisation identities ⟨Ω⊗Ω|T_λ(θ₁−2h) = f_λ ⟨Ω T_λ^{(1)}|⊗⟨Ω| and site 2), and `Paul/Mathematica/new-conventions/
L2-diagonal-frame-overlap.md` Section 13 (symbolic localisation factors). Existing Python: `Paul/Python/gl3_reps.py` (GT generators).

**Files you own** (under `Paul/Python/gaudin_sl3/`): `xxx_chain.py` (XXX chain port), `xxx_limit.py` (Tasks X1–X4),
`RESULTS-xxx-limit.md`. Do not touch `gaudin_chain.py`, `qside*.py`, `RESULTS-Qside.md` (another agent's). `python -X utf8`,
sympy/mpmath/numpy; exact `Fraction` object arrays or sympy `DomainMatrix` for 64×64 operator products (dense sympy `Matrix` is
too slow); symbolic h only at the 3×3 / scalar level. Time one 64×64 exact quantum minor before loops.

## 0. The XXX chain (new conventions) and the Gaudin scaling

XXX: ℒ_ij(u) = u δ_ij + h E_ij; T(u) = ℒ⁽¹⁾(u−θ₁)ℒ⁽²⁾(u−θ₂) (auxiliary matrix product, operator entries E⁽¹⁾ = E⊗1, E⁽²⁾ = 1⊗E);
twist on the right T^G = T·G; quantum minors qmin_{I,J}(u) = Σ_σ sgn σ ∏_k T_{σ(i_k), j_k}(u − (k−1)h) (rows antisymmetrised,
column k at u−(k−1)h); t_a(u) = Σ_{|I|=a} qmin_{I,I}(u); ν_n(u) = ∏_s(u−θ_s+hΛ_n). θ = (2/7, −3/5), Λ = (2,1,0). Scaling:
z = exp(hζ) (so G = diag(z) = exp(h·diag(ζ)) in the diagonal frame), w = exp(hω) for the bra, ζ = (3/2, 5/4, 2/3),
ω = (5/2, 7/3, 3/5), and h ∈ {10⁻², 10⁻³, 10⁻⁴} (rational h, exact arithmetic with z = e^{hζ} replaced by its truncated series
1 + hζ + h²ζ²/2 + h³ζ³/6 + h⁴ζ⁴/24 when exactness matters — say which you used; or mpmath at 60 digits with exact exponentials).
Gaudin side (same definitions as the plan §1; implement your own minimal version): L_g(u) = g + Σ_s E⁽ˢ⁾/(u−θ_s),
cdet(∂ − L_g) = ∂³ − S₁∂² + S₂∂ − S₃ (column-ordered, ∂ to the right).

## 1. Tasks

**X1 — frames and the vacuum.** (a) r_k(z,w) (plan §1 formulas with z, w) at z = e^{hζ}, w = e^{hω}: confirm numerically
  r_k(z,w) → r_k(ζ,ω) and give the O(h) correction (it should be O(h) and simple). (b) XXX P_aux(z) (rows (1, z_i−χ₁, z_i²−χ₁z_i+χ₂))
  and G′(z,w) = P_aux(z)G_c(w)P_aux(z)⁻¹: compute (G′ − 1)/h as h → 0 (does it converge? to which 3×3 matrix?); compare with
  g′ := P_aux^G(ζ)g_c(ω)P_aux^G(ζ)⁻¹ (Lie-algebra companion g_c(ω) = [[e₁,−e₂,e₃],[1,0,0],[0,1,0]], P_aux^G rows (1, ζ_i−e₁, ζ_i²−e₁ζ_i+e₂))
  and with its transpose / inverse variants; report which (if any) matches, or that G′ has no regular limit and what the
  limit of the vacuum eigencovector condition is instead. (c) The XXX x₀(z) = (e₁ᵀP_phys(z)⁻¹)^{⊗2} (reference `x0Of`, factorised
  inverse) at z = e^{hζ}: its leading h-power and leading direction (projectively); compare with x₀^G(ζ) built from P_aux^G(ζ) by the
  same Gauss-factor construction. Report whether they are parallel; if not, say so (the plan expects NOT in general).

**X2 — the Hamiltonian dictionary (the key convention check).** Build exactly, at a few rational u (e.g. u ∈ {3/7, −2/9, 5/3,
  −7/4}) and the three h's, the 64×64 operators A(u,h) := (t₁(u)/ν(u) − 3·1)/h and B(u,h) := (t₂(u)/(ν(u)ν(u−h)) − 3·1)/h with
  ν(u) = (u−θ₁)(u−θ₂) and G = diag(e^{hζ}); extract the h-expansion coefficients A₁, A₂, A₃, B₁, B₂, B₃ by Richardson over the
  three h's (or exactly by using symbolic h at ONE u with the truncated exponential — expensive but exact: try it at one u).
  Build the Gaudin S₁(u), S₂(u), S₃(u) (cdet) at the same u with g = diag(ζ). Find by exact (or 60-digit) linear algebra the
  relations S₂ = α₁A₂ + α₂A₁² + α₃A₁ + α₄·1 (+ B-terms if needed), S₃ = (combination of A₃, B₃, A₁A₂, A₁B₂, A₁³, B₂, …, 1) with
  u-dependent scalar coefficients — the expected structure: S₂ = ½S₁² − A₂ + c(u), S₃ from B₃ and A-products. Report the
  dictionary with the c-numbers as rational functions of u (fit over the u's). If the Gaudin cdet convention of the plan does
  not fit, test the transposed-Lax and row-ordered variants and report which the XXX limit selects.

**X3 — the localised covectors (what W1 must find).** One site at θ₁, Λ = (2,1,0), twist G = diag(z) (ket twist), bra vacuum
  Ω_zw = e₁ᵀρ(n(z,w)) (reference `ΩD`): compute ⟨v_λ| := ⟨Ω_zw|T_λ^{(1)}(θ₁−2h) for the one-site fused shapes T₍₁,₁₎ = t₂ and
  T₍₂,₁₎(v) = t₂(v)t₁(v+h) − t₃(v+h) (one-site t_a: Σ principal quantum minors of ℒ(u−θ₁)G) at z = e^{hζ}, w = e^{hω} for the
  three h's; expand ⟨v_λ| = h^{|λ|}(v₀ + h v₁ + h² v₂ + …) (Richardson) and report v₀ (expect ∝ Ω_{ζω}), v₁ and v₂ as 8-vectors
  and their directions modulo Ω_{ζω}: express v₁ in the basis {Ω·1, Ω·C_ζ, Ω·(one-site Gaudin S₃ Laurent coefficients at θ₁),
  Ω·Σ_{ij}(g−g′)_ij E_ji} (least squares, report residual) where C_ζ = Σ_kζ_kE_kk and the one-site Gaudin S₃ is cdet(∂ − ζ − E/(u−θ₁)).
  Also the spectator factors f₍₁,₁₎ = (θ₁₂−h)(θ₁₂−3h), f₍₂,₁₎ = (θ₁₂−h)²(θ₁₂−3h) → θ₁₂², θ₁₂³. Same for the site-2 pair
  (t₂(θ₂+h), t₂(θ₂+h)t₁(θ₂−h) on the one-site chain at θ₂). Then the two-site statement at small h: verify numerically that
  ⟨Ω_zw⊗Ω_zw|T_λ(θ₁−2h) = f_λ⟨v_λ|⊗⟨Ω_zw| still holds at h = 10⁻³ (sanity of your port; 64×64), and record the leading-order
  two-site operator X₁ := lim (T₍₁,₁₎(θ₁−2h)/(h²θ₁₂²) − c·1)/h (c the scalar leading term) as a 64×64 matrix — express it in
  the Gaudin Laurent coefficients {1, W, H₁, S₃ coefficients at θ₁, θ₂} if it lies in their span (report the rank test).

**X4 — Baxter → oper.** Symbolic in h (sympy, scalars only): the XXX direct equation ν₁(x)ν₁(x+h)ν₁(x+2h)Q(x−h) − t₁(x)ν₁(x+h)
  ν₁(x+2h)Q(x) + t₂(x+h)ν₁(x+2h)Q(x+h) − t₃(x+2h)Q(x+2h) = 0 and the dual ν₃(x)ν₃(x−h)ν₃(x−2h)Q(x+h) − t₁(x)ν₃(x−h)ν₃(x−2h)Q(x) +
  t₂(x)ν₃(x−2h)Q(x−h) − t₃(x)Q(x−2h) = 0 (L = 2: ν_n(x) = ∏_s(x−θ_s+hΛ_n), t₃ = χ₃ν₁(x)ν₂(x−h)ν₃(x−2h)), with t₁(x) = ν(x)[3 + hA₁(x) +
  h²A₂(x) + h³A₃(x)], t₂(x) = ν(x)ν(x−h)[3 + hB₁(x) + h²B₂(x) + h³B₃(x)] as unknown scalar functions, Q(x±kh) Taylor-expanded:
  show orders h⁰, h¹, h² vanish identically (report any condition this imposes on A₁, B₁ — e.g. B₁ = A₁-type relations), and
  print the order-h³ ODE. Substitute X2's dictionary (A₂, A₃, B₂, B₃ in terms of S₁, S₂, S₃) and compare with the Talalaev oper
  ∂³ − S₁∂² + S₂∂ − S₃ and its formal adjoint: which XXX equation (direct/dual) becomes which, and the exponent signs
  (Q_i = z_i^{−x/h}q_i → e^{−ζ_i x}q_i; Q^j → e^{+ζ_j x}q^j). This settles the P1 convention question independently; write the
  conclusion in RESULTS in one unambiguous sentence.

## 2. Reporting

`RESULTS-xxx-limit.md`: §X1–§X4 with the explicit matrices/limits and residuals (with the h-values used and the observed
convergence rate), the dictionary formulas, the v₁ decomposition, the oper/adjoint assignment. Scripts runnable top-to-bottom.
No notebooks, no git, no edits outside `Paul/Python/gaudin_sl3/`. Numbers, not adjectives. If a task's expectation is refuted
and two attempts do not fix it, stop there, keep what works, report what you saw.
