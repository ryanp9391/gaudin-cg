# Spec U2: the analytic untwisting limit of the Gaudin jet formula — ⟨v ⊗ hws|Ψ₀⟩ from Q-functions (2026-10-16)

Delegated to a general-purpose subagent (SONNET), Python only (the Wolfram kernel belongs to another agent). Read first:
`Paul/Mathematica/Gaudin/Gaudin-sl3-round4-plan.md` §0 (the formula, c, M, x₀^G, φ_u closed forms, the regime facts — binding), §2 "The difficulty"
and "U2" (your goal), §4, §6; then `Paul/Python/gaudin_sl3/RESULTS-twovac-fit.md` §F4 (D1: the regime run — branches, ε-exponents, cancellation
depth, the Richardson CG extraction; its code is your starting point), `RESULTS-twovac-limit.md` (D2: the adjoint kernel G as explicit polynomials
in (ω₁, ω₁′, 1/θ₁₂), the ν-form, homogeneity weights (2,3,3,4,5)), `RESULTS-gaudin-B.md` §B4 (D3: Θ = G/(4ω³ω′³), K tensor), and the round-2 notebook
summary cell (last cell of `Paul/Mathematica/Gaudin/Experiments/Gaudin-L2-chain-2.wb`, JSON) for the exact CG values of the adjoint (−1/√2, −1/√2,
0.1178511302, the octet vectors (0.0493116000, 0.2088872897), (0.3379867346, −0.0797878449), mixing roots of R² − 5R − 5 = 0). Exact CG tables and
Laurent data will appear during your run in `Paul/Mathematica/Gaudin/ref/untwisted-exact-<Λ>.txt` (agent W4) and the ε-series of the branch Q's and
jet monomials in `untwisted_Q_series_210.pkl` + `RESULTS-untwisted-Q.md` (agent U1) — use them when they exist, but build your own numerical ε-series
first so you are never blocked.

**Code to reuse (read-only; import or copy, never edit):** `twovac_fit.py` (D1: `d2_G(w1, w1p, thi)` = the 5×5 adjoint kernel as polynomials, `jets`,
`plucker`, `onesite_inf`, `cov_inf`, `x0_onesite`, `cscal`, `kappaG`, the F4 regime code), `qside.py`/`gaudin_chain.py` (P1: `q4_point(eps, dps)`, `set_regime`,
branch selection, Q solvers at 100+ digits), `gl3_reps.py`. **Files you own:** `untwist_limit.py` (stage flags), caches `untwist_limit_*.pkl`,
`RESULTS-untwist-limit.md`. Environment: Python 3.13 (`python -X utf8`), sympy 1.14, mpmath 1.3, numpy; no `uv`. Other agents: W4 (Wolfram — you must
NOT use Wolfram), G1 (`general_kernel*`), U1 (`untwisted_Q*`) — do not touch their files. ≤ 2 heavy processes of your own; 200-digit regime points are
slow — time one ε first and cache per (ε, branch).

## 0. Setting (adjoint; then [3,1,0]/[3,2,0] only if time)

θ = (2/7, −3/5), θ₁₂ = 31/35. Regime ζ = εφ, ω₁ = εω̂₁, ω₁′ = εω̂₁′, φ = (1/3, −2/7, 5/11), grid ω̂₁ ∈ {3/4, −5/3, 1/2, −7/3} × ω̂₁′ ∈ {−12/5, 7/2, 1/8, 5/4}
(16 points; D1 used the first 2×2), ε_k = 10⁻⁶·2^{−k}, k = 0..11, 200 digits. Six HW-limit branches B: (4,2,0), (4,1,1), (3,3,0), (3,2,1) copy 1, (3,2,1)
copy 2, (2,2,2), G-norm 1 and the round-1 sign convention (first nonzero component in GT order positive), copies ordered by the S3c[1,2] eigenvalue.
The jet formula (adjoint): c·Qnum·M₁M₂ = Σ_{ij} G_ij(ω₁,ω₁′,θ₁₂) N_i(θ₁)N_j(θ₂), N = (P³p₀₂, P³′p₀₂, P³p₀₃, P³′p₀₃, P³p₂₃), M = P³p₀₂, c = (ζ₂−ω₁)(ζ₃−ω₁)²·
(ζ₂−ω₁′)(ζ₃−ω₁′)²(ζ₁−ζ₂)², Qnum = ⟨Φ∞|Ψ⟩/⟨x₀^G|Ψ⟩; G_ij is homogeneous in (ω₁, ω₁′, 1/θ₁₂) of degree 10 − w_i − w_j, w = (2,3,3,4,5), and polynomial
in (ω₁, ω₁′) — so G_ij(εω̂₁, εω̂₁′, θ₁₂) = Σ_d ε^d G_ij^{[d]}(ω̂₁, ω̂₁′, θ₁₂) with d the (ω₁,ω₁′)-degree (write G_ij as a sum of (ω-degree d) pieces; the
1/θ₁₂ powers are NOT scaled). Chain side: ⟨Φ∞|Ψ_B⟩ = ⟨φ_u⊗e₁ᵀ|Ψ₀⟩(1 + O(ε)) with φ_u = e₁ᵀρ(exp(a′E₁₂)exp(g′E₁₃)exp(c′E₂₃)) and the ε-free a′, c′, g′
of the plan §0 (functions of φ, ω̂₁, ω̂₁′ with denominators (φ₂−ω̂₁), (φ₂−ω̂₁′), (φ₃−ω̂₁), (φ₃−ω̂₁′), (φ₁−φ₂)); ĉ(φ,ω̂₁,ω̂₁′) := c/ε⁸ = (φ₂−ω̂₁)(φ₃−ω̂₁)²·
(φ₂−ω̂₁′)(φ₃−ω̂₁′)²(φ₁−φ₂)² clears exactly those denominators.

## 1. Tasks

**L1 — ε-order bookkeeping (numbers first).** For each branch and ε_k: the branch state (P1's machinery; check continuity in ε of the branch
selection and sign), ⟨Φ∞|Ψ_B⟩ at the 16 grid points (→ its ε-exponent, expect 0, and the limit ⟨φ_u⊗e₁|Ψ₀⟩ computed directly from the exact HW ket —
D1 did this), ⟨x₀^G|Ψ_B⟩ (exponent 2 + k_B: find k_B per branch; test the interpretation ⟨x(φ)⊗x(φ)|Ψ₀⟩ = 0 for non-top μ by computing it exactly on
the HW kets; if zero, identify which order of Ψ_B first contributes), M₁M₂ (exponent), the five N_i(θ_s) per site (exponents n_i(B,s)), and the total
Σ G N N (exponent N_B; D1: 8,7,7,6,6,4). Tabulate; check N_B = 8 + 0 − (2 + k_B) + (exponent of M₁M₂) consistency with the LHS.

**L2 — the leading polynomial and the normalisation-free limit.** (a) Numerically extract the ε-series of each jet monomial N_i(θ_s;ε) = Σ_n ε^n
N_i^{(n)} (coefficient-wise polynomial fit in ε over the 12 ε_k at 200 digits; report the reliable depth; cross-check with U1's pkl when present).
(b) The leading polynomial L_B(ω̂₁,ω̂₁′) := Σ_{ij} Σ_{d+n+m = N_B} G_ij^{[d]}(ω̂₁,ω̂₁′,θ₁₂)·N_i^{(n)}(θ₁)·N_j^{(m)}(θ₂) — a polynomial in (ω̂₁, ω̂₁′) with 200-digit
coefficients; verify that all lower-order combinations (d+n+m < N_B) cancel to the working precision (this is the analytic content of the
cancellation depth), and list which (d, n, m) triples contribute at order N_B — the "depth" of the limit per branch. (c) The identity to establish:
L_B(ω̂₁,ω̂₁′) = κ_B · ĉ(φ,ω̂₁,ω̂₁′) · ⟨e₁ᵀρ(n′(ω̂₁,ω̂₁′))⊗e₁ᵀ|Ψ₀⟩ · (M-side factor) as polynomials in (ω̂₁, ω̂₁′), where the right side is computed EXACTLY from
the exact HW ket (G-norm 1) and the closed-form Gauss parameters (sympy: expand ĉ·e₁ᵀρ(n′) — ĉ clears the denominators — into a polynomial in ω̂₁, ω̂₁′
whose coefficients are linear in the CG coefficients C_ν = ⟨e_ν⊗e₁|Ψ₀⟩), and the M-side factor is the ε-leading coefficient of M₁M₂/(⟨x₀^G|Ψ_B⟩/ε^{2+k_B})
— state it. Compare coefficient by coefficient (κ_B a single branch constant); report the agreement precision for all six branches, incl. the
two octet copies (the limit covector ⟨φ_u ⊗ e₁| pairs with each copy separately — two generating functions). This is the analytic untwisting of
the formula at the level of the N^{(n)}. (d) Normalisation-free check (D1's route): ρ_B(ω̂;ω̂⁰) = [ĉ(ω̂⁰)/ĉ(ω̂)]·L_B(ω̂)/L_B(ω̂⁰) vs the exact ratio of
⟨φ_u⊗e₁|Ψ₀⟩'s — and the CG coefficients C_ν solved from the coefficients of L_B/(κ_Bĉ·M-factor) as the generating-function coefficients: compare to
the exact C_ν (the "CG from Q's" dictionary made explicit).

**L3 — the formula in terms of Q-functions.** Using U1's results (the untwisted solution-space basis, the exact p^{(0)}, p^{(1)}, …, the closed-form
Laurent data) as far as they exist: express the contributing N_i^{(n)} in terms of the untwisted Q's and their ε-derivatives at θ_s (Plücker
coordinates of the ζ-jet basis), and write the explicit generating function ⟨φ_u(ω̂₁,ω̂₁′)⊗e₁|Ψ₀⟩ ∝ L_B/ĉ for each adjoint branch in those terms; for
the 1-dim slices give the single CG number, for the octets the 2-vectors and the mixing. Check φ-independence: repeat L2 with φ′ = (−2/5, 3/8, 7/11)
(the N^{(n)} change, ĉ changes, the C_ν must not). State clearly which objects the final formula needs: the kernel G, the untwisted Q's to which ε-order,
and the Gauss parametrisation. If U1's exact p^{(n)} are not available, give the formula in terms of the numerically extracted N^{(n)} and say so.

**L4 (only if time) — [3,1,0] or [3,2,0].** Same L1–L2 with P4's kernels (`twovac_fit_320_grid.pkl`; [3,1,0] base-point kernel is insufficient for the
ω-degree split unless G1's polynomial laws exist in `RESULTS-general-kernel.md`) and P3/P4's regime machinery; the exact CG tables from W4's files or
the round-2/3 notebook summaries.

## 2. Reporting

`RESULTS-untwist-limit.md`: L1 (the order table per branch), L2 (the leading polynomial identity: contributing (d,n,m) triples, agreement per branch,
the CG dictionary), L3 (the formula in Q terms, φ-independence), L4 if done, surprises, not done. Numbers with precision; state every normalisation
and sign convention used. No notebooks, no git, no Wolfram, no edits outside your files. If an expectation is refuted and two attempts do not fix it,
stop that item, keep what works, report what you saw.
