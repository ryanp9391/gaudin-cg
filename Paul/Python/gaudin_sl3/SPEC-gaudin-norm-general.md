# Spec Q2: the Gaudin norm formula for [3,1,0] and [3,2,0] — exponent rule, constant, sign (round 5 follow-up to WP-D; 2026-10-16)

Delegated to a general-purpose subagent (SONNET), Python only (no Wolfram, no wolfbook tools). Read first: `RESULTS-gaudin-norm.md` (Q1: the
adjoint result you generalise — read §1 (ii) and the D0 conventions carefully), `RESULTS-nu-general.md` §1 (N1: the general-Λ no-log conditions,
3 per site, in `nu_general.py`: `nu_symbolic(l1,l2)`, `frobenius(lam)`, `nu_numeric`), `Paul/Mathematica/Gaudin/Gaudin-sl3-round4-plan.md` §0
(notation — binding), the notes `Paul/Mathematica/Gaudin/Gaudin-sl3-CG-notes.md` §1.4 (frames, x₀^G, G-metric), §2.3. Environment: Python 3.13
(`python -X utf8`), sympy 1.14, mpmath 1.3, numpy; no `uv`; run from `Paul/Python/gaudin_sl3/`. Other agents are running (M1, M2, S1 in Python; W5
on Wolfram): ≤ 2 heavy Python processes; time one unit before loops; cache to `gaudin_norm_gen_*.pkl`.

**Code to reuse (read-only; import or copy, never edit):** `gaudin_norm.py` (Q1: the D0/D1 machinery — Jacobian of the no-log conditions, N(Ψ),
the state-independence test, `--stages`), `nu_general.py` (N1: the general-Λ Frobenius/no-log conditions and the Laurent data per state — N1's
caches `nu_general_res_*.pkl` hold per-state Laurent data for [3,1,0]/[3,2,0] at G1's base points; check their format), `general_kernel.py`
(`get_point2`/`build_point2` — the cached 80-digit points `general_kernel_pt_310_<tag>_*.pkl`, `_320_*` with 225 states each at tags P, Pp, X1; `x_cov`,
`site_monos`), `twovac_fit_320.py`/`qside_320.py`/`qside_310.py` (rep-generic chain, G-metric per rep: find the Gram/metric routine used for
G-normalisation — W4's `ref/untwisted-exact-310.txt` GMET1 record gives the exact one-site diagonal metric in the GT order, use it to cross-check),
`gaudin_chain.py`, `gl3_reps.py`. **Files you own:** `gaudin_norm_gen.py`, `gaudin_norm_gen_*.pkl`, `RESULTS-gaudin-norm-general.md` (write
incrementally). No notebooks, no git, no Wolfram, no edits outside your files.

## 1. Tasks

**N-gen — the formula for general Λ.** Q1 (adjoint): N(Ψ) := ⟨Ψ|Ψ⟩_G/⟨x₀^G|Ψ⟩² = (9/4)·|det J₆|·(ζ₁−ζ₂)^{μ₃−6}(ζ₁−ζ₃)^{μ₂−6}(ζ₂−ζ₃)^{μ₁−6}, J₆ = ∂(F₁..F₆)/∂(H₁,a₁,b₁,H₂,a₂,b₂)
with F the six no-log conditions, μ = the sector (dual-Q degrees). For [3,1,0] and [3,2,0] (225 states each, three twists P, Pp, X1, θ = (2/7,−3/5)):
compute N(Ψ) (G-metric per rep — verify your metric on the HW kets against W4's exact G-norm² records), the Laurent data per state, J₆ from N1's
conditions (general λ₂: condition A involves the other site for λ₂ = 2 — the Jacobian must include those cross-derivatives), and test
N/(|det J₆|·Π_{i<j}(ζ_i−ζ_j)^{p_ij(μ)}) for state-independence with the exponents as unknowns: fit p_ij as affine functions of μ (expect
p₂₃ = μ₁ − c, p₁₃ = μ₂ − c, p₁₂ = μ₃ − c with c = 2|Λ| (= 6 for the adjoint) or c = 2λ₁ + λ₂ or similar — determine), and the constant (adjoint 9/4 =
(c₂/2)²? or (3/2)²? — express it in λ's: compute it for both reps at ≥ 10 ζ points, check θ-independence at 3 θ values). Report the formula and
the residuals (expect ~1e-70 at 80 digits). If the Jacobian normalisation of the conditions matters (the conditions are defined up to constant
factors — N1's polynomials vs U1's), say which normalisation you used and give the constant for it.

**Sign.** det J₆ has no definite sign on the adjoint (32/32). Test two structural hypotheses on all three reps: (a) sign(det J₆) = sign of the
product of the two per-site 3×3 Jacobian blocks (is J₆ block-triangular up to the cross terms? compare det J₆ with det J₃⁽¹⁾·det J₃⁽²⁾ in sign
and in value); (b) sign(det J₆) = (−1)^{(number of states "below" Ψ in a Bethe-type ordering)} — skip (b) if (a) is not informative; just report the
sign pattern vs the sector μ and copy. Do not spend more than a quarter of your effort here.

**Regime (only if time remains).** On N1's regime caches for [3,1,0] (`nu_general_reg_*.pkl`; HW-limit branches, 24 ε's), check lim ε^{p}·N = 1/X_B²
with p = 2(2+k_B)? (adjoint: 4 + 2k_B — find the general p, presumably 2(λ₂+1+k_B) or 2(2λ₂+k_B)… determine from the ε-order of ⟨x₀^G|Ψ_B⟩ = ε^{2λ₂+k_B}X_B?)
and compare with N1's X_B table.

## 2. Reporting

`RESULTS-gaudin-norm-general.md`: §0 conventions (metric, condition normalisation), §1 the general formula with exponents and constant, residuals
per rep/twist, §2 the sign test, §3 the regime check (if done), surprises, not done.
