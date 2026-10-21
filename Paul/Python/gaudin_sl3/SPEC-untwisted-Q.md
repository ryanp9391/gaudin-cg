# Spec U1: the untwisted Q-system of the sl(3) Gaudin chain and the ε-expansion of the HW-limit branches' Q-functions (2026-10-16)

Delegated to a general-purpose subagent (SONNET), Python only (the Wolfram kernel belongs to another agent). Read first:
`Paul/Mathematica/Gaudin/Gaudin-sl3-round4-plan.md` §0 (notation), §2 "The difficulty" and "U1" (your goal), §4, §6; then
`Paul/Python/gaudin_sl3/RESULTS-Qside.md` (P1: conventions, the Q solvers, the Q4 regime machinery — branches, their Q's at ε = 10⁻⁶),
`RESULTS-Qside-310.md` and `RESULTS-Qside-320.md` (P3/P4: no-log conditions = algebraic equations on the Laurent data, κ_s = a_s/c₂, κ₂,_s, local
exponents {0, λ₂+1, λ₁+2}), `RESULTS-twovac-fit.md` §F4 (D1: regime orders), the round-1 plan `Gaudin-sl3-CG-plan.md` §1 (dictionary) and the
round-2/3 verdict blocks. Chain-side exact reference data will appear during your run in `Paul/Mathematica/Gaudin/ref/untwisted-exact-<Λ>.txt`
(agent W4: exact HW kets, CG tables, Laurent data at ζ = 0 with symbolic θ, Casimirs, perturbative Laurent data to O(ε²)) — read them when they
exist; do not wait for them.

**Code to reuse (read-only; import or copy, never edit):** `gaudin_chain.py`, `qside.py` (P1: `Chain`, `Operators`, `eigenstates`, Q solvers, `set_regime`,
`raising_residual`, `q4_point(eps, dps)`), `qside_310.py`/`qside_320.py` (rep-generic versions, `frobenius_polys`), `twovac_fit.py` (D1: `jets`, `plucker`,
`StateJets`). **Files you own:** `untwisted_Q.py` (stage flags), caches `untwisted_Q_*.pkl`, `RESULTS-untwisted-Q.md`, and a machine-readable export
`untwisted_Q_series_<Λ>.pkl` (the ε-series of the branch Q's and jet monomials, for agent U2 — document its format in the RESULTS file). Environment:
Python 3.13 (`python -X utf8`), sympy 1.14, mpmath 1.3, numpy; no `uv`. Other agents: W4 (Wolfram — you must NOT use Wolfram), G1 (`general_kernel*`),
U2 (`untwist_limit*`) — do not touch their files. ≤ 2 heavy processes of your own at a time; the 200-digit regime runs are slow — time one ε first.

## 0. Setting

Λ = (2,1,0) first (then (3,1,0), (3,2,0) if time), θ = (2/7, −3/5). Regime ζ = εφ, φ = (1/3, −2/7, 5/11) (second direction φ′ = (−2/5, 3/8, 7/11) for
checks); the ket side does not involve ω. HW sectors of V⊗V for the adjoint: (4,2,0) top, (4,1,1), (3,3,0), (3,2,1)×2 (copies = eigenvectors of the
θ-universal matrix (3/(2θ₁₂))[[3,−1],[−1,7]] on the exact HW space), (2,2,2). The HW-limit branch of sector μ is the eigenstate of the twisted family
with the smallest raising residual (P1's `q4_point` selects it); Ψ_B = Ψ₀ + εΨ₁ + … with Ψ₀ the exact HW vector (G-norm 1, round-1 sign convention).
Dual Q's of the branch P^j(u;ε) = e^{εφ_j u}p_j(u;ε), monic, deg p_j = μ_j; direct Q_i = e^{−εφ_i u}q_i, deg 4 − μ_i.

## 1. Tasks

**T0 — untwisted Laurent data in closed form.** At ζ = 0 the Laurent data of the HW state of weight μ: H₂ = −H₁ (exact: H₁ + H₂ = −W + 2|Λ|Σζ);
derive H₁ = (|Λ|² − ½[C₂(μ) − 2C₂(Λ)])/θ₁₂ from Σ_ik E¹_ikE²_ki = ½(C₂^{tot} − C₂⁽¹⁾ − C₂⁽²⁾), C₂(λ) = Σ_i λ_i(λ_i + 4 − 2i) (gl(3) quadratic Casimir in the
E_ikE_ki normalisation — check the normalisation on the top state where Σ_ik E¹_ikE²_ki = Σλ_i² (W2)). For a_s, b_s (the (u−θ_s)⁻² and (u−θ_s)⁻¹ data
of S₃): compute them EXACTLY at ζ = 0 on exact HW kets (build the exact HW vectors of V⊗V by solving E₁₂, E₂₃ annihilation with Fractions on the
weight spaces; the S₃ Laurent coefficient operators are exact rational matrices; for the multiplicity-2 sector take the universal copies), for the
adjoint, [3,1,0], [3,2,0], [2,2,0], [1,1,0] and all their sectors; then find closed forms in (Λ, μ, θ₁, θ₂) — expect a_s ∝ 1/θ₁₂ (or 1/θ₁₂-free?) and b_s
with 1/θ₁₂ and 1/θ₁₂² pieces; use the cubic Casimir C₃ (compute its eigenvalue on each sector exactly: C₃ = Σ E_ijE_jkE_ki) as a candidate variable;
report the sum rules b₁ + b₂ and a₁ + a₂ + θ₁b₁ + θ₂b₂ (they govern whether the untwisted oper is Fuchsian at ∞: S₃ = O(u⁻³) needs both to vanish).
Cross-check against W4's `untwisted-exact-*.txt` when available.

**T1 — the untwisted oper and its solution space.** For each adjoint HW sector: the ζ = 0 oper D₀ = ∂³ − S₁∂² + S₂∂ − S₃ with S₁ = 3/(u−θ₁) + 3/(u−θ₂),
S₂ = Σ_s[H_s/(u−θ_s) + 3/(u−θ_s)²], S₃ = Σ_s[b_s/(u−θ_s) + a_s/(u−θ_s)²] (e_k(0) = 0). (a) Local exponents at θ₁, θ₂ ({0,2,4} expected) and at ∞
(compute from the u → ∞ behaviour; if the equation is Fuchsian at ∞, the Fuchs relation Σ exponents = 3 must hold — check; if not, characterise the
irregularity). (b) Polynomial solutions: for each degree d ≤ 2λ₁ + 2 solve D₀p = 0 exactly (linear algebra on coefficients): the degrees that admit
solutions and the dimension of the polynomial solution space (expect fewer than 3 independent polynomials, since the three twisted P^j degenerate).
(c) The ε → 0 limits of the branch Q's: compute the regime Q's at ε_k = 10⁻⁶·2^{−k}, k = 0..7 (200 digits) and extrapolate p_j(u;ε) → p_j^{(0)}(u)
coefficient-wise (Richardson); verify the p_j^{(0)} solve D₀ and find their linear dependencies; the limits of the differences (p_a − p_b)/(ε(φ_a − φ_b))
and of ∂_ε(e^{εφ_j u}p_j) — do they solve D₀, or an inhomogeneous equation? Establish a basis of the limiting 3-dimensional solution space
{P^j(u;ε)} as ε → 0 (the "ζ-jet" basis) in closed form (polynomials, possibly of degrees > μ_j). (d) Closed forms: if D₀ is Fuchsian with three
singular points, match it to the generalised hypergeometric ₃F₂ equation after the Möbius map x = (u−θ₁)/(u−θ₂) or x = (u−θ₁)/(θ₂−θ₁) (sympy:
`hyperexpand`, or compare coefficients with the ₃F₂ differential equation; the exponents at 0, 1, ∞ fix the parameters) — give the untwisted Q's
as ₃F₂ polynomials (Jacobi-type) if they are; otherwise say so. Repeat (a)–(c) for one sector each of [3,1,0] and [3,2,0] if time.

**T2 — the ε-expansion of the branch Q's.** (a) Numerically: for each adjoint branch, P^j(u;ε) = e^{εφ_j u}Σ_{n≥0}ε^n p_j^{(n)}(u): fit each polynomial
coefficient of p_j(u;ε) as a polynomial in ε of degree ≤ 10 from the 12 values ε_k = 10⁻⁶·2^{−k}, k = 0..11 at 200 digits (check stability by dropping
points; report the reliable order n_max and the error per order); note the branch selection and sign/normalisation must be continuous in ε (monic
p_j; the branch is identified by sector and copy). (b) Exactly, first orders: the branch's Laurent data H_s(ε), a_s(ε), b_s(ε) to O(ε²) — from W4's
exact perturbation theory if available, else from the no-log conditions (P3/P4 `frobenius_polys`: the Gaudin "Bethe equations" on the six Laurent
data — expand them around the ζ = 0 solution with ζ = εφ and solve order by order; the twist enters S₂, S₃ through e₂(ζ), e₃(ζ) and the W-term:
H₁ + H₂ = −W + 2|Λ|Σζ with W = Σ_kζ_k(E_tot)_kk = ε Σ_k φ_kμ_k on the sector); then the oper order by order: D₀p^{(1)} = −D₁p^{(0)} etc. (D_n the ε^n
part of the oper, including the e^{εφ_j u} factor's contribution: write P = e^{εφ_j u}p and expand) — polynomial linear algebra, exact in φ and θ;
compare p^{(1)}, p^{(2)} with (a). (c) The ε-expansion of the jet monomials: N_i(θ_s;ε) = (P³p_kl)(θ_s) for the five adjoint pairs (02),(12),(03),(13),(23)
and M = P³p₀₂: N_i = Σ_n ε^n N_i^{(n)} to the reliable order, per branch and site, from (a); also p₀₁ ≡ 0 and the Frobenius relation p₁₂ = κp₀₂ with
κ = a_s/3 as a consistency check order by order. Export everything U2 needs in `untwisted_Q_series_210.pkl` (document keys: branch label, site, pair,
list of mp coefficients N^{(n)} as strings at 200 digits; plus p_j^{(n)} polynomial coefficients) and describe the format in RESULTS.

**T3 — closed forms (as far as they go).** Express p^{(0)}, p^{(1)} (and N^{(0)}, N^{(1)}, …) for the adjoint branches in terms of the T1 basis of the
untwisted solution space and the T0 Laurent data; identify which combinations are φ-independent. If T1(d) succeeded, give the untwisted Q's as
explicit hypergeometric polynomials for each sector.

## 2. Reporting

`RESULTS-untwisted-Q.md`: T0 (closed forms, verification per sector/rep, sum rules), T1 (exponents incl. at ∞, polynomial solutions, the ζ-jet basis,
hypergeometric identification or its failure), T2 (series tables per branch — orders, reliability, the exact first orders, the pkl format), T3,
surprises, not done. Exact where possible; numbers with precision. No notebooks, no git, no Wolfram, no edits outside your files. If an expectation
is refuted and two attempts do not fix it, stop that item, keep what works, report what you saw.
