# Spec P4: the Q side of the sl(3) Gaudin chain for Λ = (3,2,0), and the [3,2,0] jet formula (2026-10-16)

Delegated to a general-purpose subagent (sonnet), Python only (the Wolfram kernel belongs to another agent). Read first:
`Paul/Mathematica/Gaudin/Gaudin-sl3-round3-plan.md` §0 (the prediction table — rows "S₂ double pole" to "det N" and "jet formula" are yours),
§4, §6; then `Paul/Python/gaudin_sl3/RESULTS-Qside-310.md` (P3: the (3,1,0) run you are repeating for its dual, incl. the gauge-power scan, the
local exponents, the Frobenius no-log conditions, the Lagrange columns, the rank tallies) and `RESULTS-Qside.md` (P1: conventions); for the
jet formula `RESULTS-twovac-fit.md` (D1: the fit pipeline, the hypothesis ladder, κ_s = a_s/3, normalisation by c and M = P³p₀₂, the κ^G scalar)
and `RESULTS-twovac-limit.md` (D2: the formula and its structure), and the round-2 plan `Gaudin-sl3-round2-plan.md` §1.

**Code to reuse (read-only; import or copy, never edit):** `qside_310.py` (P3; `LAM`, `build_states(lam, zeta, dps)`, `vac_data(…, lam1=…)`,
`q1_pscan`, `local_exponents`, `frobenius_polys`, `Nrows(…, gauge)`, `sectionR` — it has `if __name__ == '__main__'`, so importing is safe; copy it
to `qside_320.py` and set LAM = (3,2,0), then fix anything hard-wired to λ₂ = 1), `gaudin_chain.py`, `qside.py`, `gl3_reps.py`, and for the fit
`twovac_fit.py` (D1: `jets`, `plucker`, `onesite_inf`, `cov_inf`, `x0_onesite`, `kappaG`, `cscal`, the `PivQR` high-precision least squares, `design`,
`fit_level`, `rat_recon`). **Files you own:** `qside_320.py`, `twovac_fit_320.py` (both runnable top-to-bottom with `python -X utf8`), caches
`qside_320_*.pkl`, `twovac_fit_320_*.pkl`, `RESULTS-Qside-320.md`. Environment: Python 3.13, sympy 1.14, mpmath 1.3, numpy; no `uv`. (3,2,0):
dim 15, 225 states (P3's (3,1,0) run took ~15 min cold); time first.

## 0. The model

Λ = (3,2,0) on both sites; θ = (2/7, −3/5), ζ = (3/2, 5/4, 2/3), ω = (5/2, 7/3, 3/5), ω′ = (7/4, −2/5, 9/7); find and report the HW index in
gl3_reps order; cdet as before (checks exact); Laurent structure: S₂ double pole (predicted 8 = λ₂(λ₁+1)), S₃ triple pole (predicted none), six
non-scalar data; vacuum eigencovector (ρ via the exponential series; E⁴ = 0 here) and its exact eigenvalue functions (compare with W2's closed
forms: H_s^{vac} = Σ_i(|Λ|−λ_i)ω_i ± k/θ₁₂ with k predicted 12, a_s^{vac} = 5ω₃, b_s^{vac} = ω₃(H_s − 5ω₃)).

## 1. Tasks

**Q1 — opers, degrees, gauge, exponents.** All 225 states; dual degrees μ_j; direct degrees 6 − μ_i; gauge power p ∈ {0..5} (predicted p = λ₁−λ₃ = 3
— [3,2,0] and [3,1,0] have the same λ₁−λ₃, so this does not separate the rules yet; say so); local exponents of the dual oper at θ_s (predicted
{0, 3, 5} = {0, λ₂+1, λ₁+2}) and their consequence: P′(θ_s) = κ_sP(θ_s) AND P″(θ_s) = κ₂,_sP(θ_s) for every solution — derive κ_s and κ₂,_s from the
indicial recursion in terms of the Laurent data (for the adjoint κ_s = a_s/3; find the [3,2,0] formulas), verify on all states, and the rank of the
3×5 jet matrix [P^{j(k)}(θ_s)]_{k≤4} (predicted 2) and 3×6 (predicted 3). No-log (Frobenius) conditions on the Laurent data (P3's
`frobenius_polys` generalised): how many per site, verified on all states. Vacuum Q degrees (direct and dual). Negative controls.

**Q2 — Wronskians, quantum Weyl.** Direct F (predicted Π^{λ₁−λ₂} = Π), dual F (predicted Π), c_ab = ζ_a − ζ_b; W = det[P,P′,P″] ∝ e^{Σζu}Π^{|Λ|} (predicted
Π⁵); residuals with precision.

**Q3 — Lagrange identity, det N, ratio identity.** The six residue columns with Π^{−p} weights (p from Q1), identity residual over ≥ 300 pairs; rank
tallies (A = B: 6; A ≠ B: 5, and 4 exactly when Δμ has one positive and two negative components — P3's rule — verify the tally); det N 6×6 with rows
(s,(i,j)), (i,j) ∈ {(1,3),(1,2),(2,3)}; the x₀-free ratio identity on all 225 states (max relative deviation); the 84-triple scan; controls.

**Q4 (if time) — regime.** As P3's Q4 for the 10 branches (sectors = dual images of [3,1,0]'s; multiplicity 2 in (5,4,1) and (5,3,2)).

**J — the [3,2,0] jet formula (stretch goal; budget ≤ 2 h of compute, stop with diagnostics if not converged).** Generalise D1's pipeline:
(a) objects: Φ∞ = Ω^∞[ζ,ω₁]⊗Ω^∞[ζ,ω₁′] with the rep-independent r_k^∞; x₀^G = (e₁ᵀP_phys^{G −1})^{⊗2} for this rep (D1's `x0_onesite` if rep-generic;
otherwise rebuild from P_aux(ζ) via the Gauss factorisation); Qnum = ⟨Φ∞|Ψ⟩/⟨x₀^G|Ψ⟩ at ζ ∈ {P, P′, +2–4 extra ζ with real spectra}; the frame-tie
scalar c = ω₁^{?}ω₁′^{?}/(κ^Gκ^G′) with κ^G = Δ₁^{λ₁−λ₂}Δ₂^{λ₂−λ₃}Δ₃^{λ₃} = Δ₁Δ₂² (W2's Δ's; adjoint: c = (ζ₂−ω₁)(ζ₃−ω₁)²(ζ₂−ω₁′)(ζ₃−ω₁′)²(ζ₁−ζ₂)² = ω₁³ω₁′³/(κ^Gκ^G′)
— derive the [3,2,0] analogue and verify it removes the ζ-dependence by comparing per-ζ fits); jets of P^j at θ_s to order 5; Plücker p_kl, 0 ≤ k < l ≤ 4
(10; p₀₁ = p₀₂ = p₁₂ = 0 by the two forced jets — zero them explicitly); lone P³^{(c)}, c ≤ 2 (P³′ = κP³, P³″ = κ₂P³ — keep them as separate monomials: the
kernel must be state-independent, so the κ-dependence lives in the monomials). (b) vacuum monomial M_s: first candidate P³(θ_s)p₀₃(θ_s) (the first
nonvanishing Plücker coordinate); test also P³p₀₄, P³′p₀₃ by the state-independence of per-ζ fits. (c) hypothesis ladder as D1's F2: per-site
bases B₁ = {P³^{(c)}p_kl : c ≤ 2, (k,l) ∈ {(0,3),(0,4),(1,3),(1,4),(2,3),(2,4),(3,4)}} (21 per site, 441 unknowns; rows 225 × ≥ 3 ζ = 675 — time the
QR first; D1 did 512 × 576 at 80 digits) → reduce by the ν-variables if needed: ν = (1, κ, κ₂, π₄ := p₀₄/p₀₃, π₃₄ := p₃₄/p₀₃, products…) with target
c·Qnum·(P³p₀₃)₁(P³p₀₃)₂ as a bilinear form in ν-monomials of bounded weight (the adjoint has weights (2,3,3,4,5) and total weight 10; find the
[3,2,0] analogue by trying total weights up to ~14). Report: residual, rank, null space, the number of monomials in the sparse representative
(the adjoint has 5 per site), the kernel at the base point (rational?), and whether it is site-exchange symmetric. If the fit converges, do the
polynomial-law fit over a 2×2 (ω₁,ω₁′) grid and one θ variant as D1's F3 (optional). If time remains after that, run the SAME pipeline for
[3,1,0] (P3's `qside_310.py` states and Q's; lone c ≤ 1 there since λ₂ = 1) and report whether the monomial count differs from [3,2,0].

## 2. Reporting

`RESULTS-Qside-320.md`: §0 (HW index, Laurent structure incl. c₂ and the triple pole, vacuum data vs W2's closed forms), Q1 (gauge, exponents,
κ_s and κ₂,_s formulas, jet ranks, no-log conditions, degree table), Q2, Q3 (tallies, det N, 225-state ratio test, triple scan), Q4 if done, J (the
ladder table, the kernel or the diagnostics, the [3,1,0] comparison if done), surprises, not done. Scripts runnable top-to-bottom; no notebooks, no
git, no Wolfram, no edits outside your files. If an expectation is refuted and two attempts do not fix it, stop there, keep what works, report.
