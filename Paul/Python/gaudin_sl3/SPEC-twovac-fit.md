# Spec D1: the two-vacua Q-number of the sl(3) Gaudin chain by direct fit — jets, Plücker coordinates, kernel laws, Q-side CG (2026-10-16)

Delegated to a general-purpose subagent (sonnet), Python only (the Wolfram kernel belongs to another agent). Read first:
`Paul/Mathematica/Gaudin/Gaudin-sl3-round2-plan.md` §0, §1 (binding: the target, the ladder → jet dictionary, the ansatz, the numeric
points), §4, §8; then `Paul/Python/gaudin_sl3/RESULTS-Qside.md` (P1: conventions, the det N columns, the Q solvers' behaviour) and
`RESULTS-xxx-limit.md` §X1 (x₀^G). Background on the XXX formula you are finding the Gaudin version of:
`Paul/Mathematica/new-conventions/Two-site-rotation-brackets-notes.md` §1 (K tensor, ladders, Bk2), §5 (the ten-product formula, scaling
law, Q-side CG) and `Paul/Mathematica/Gaudin/ref/xxx-210-infinf-formula-cells.txt` (cell 70 = summary). You do not need the XXX chain.

**Code to reuse (read-only, import; never edit):** `gaudin_chain.py` (Chain, Operators, eigenstates, vacuum_covector, gauss_r, Pmat, gc,
gprime), `qside.py` (Q solvers Q1, Qfun objects, vac_data, the det N machinery `Nrows`/`build_N_det`, the regime machinery of Q4:
set_regime, raising_residual, q4_point; HW index 7 in gl3_reps order), `xxx_limit.py` (`mp_x0_onesite` and the P_aux^G Gauss construction
for x₀^G — or re-implement x₀^G = (e₁ᵀP_phys^{G −1})^{⊗2} with P_phys^G the rep image of (P_aux(ζ)⁻¹)ᵀ Gauss-factorised as P₋·diag·P₊;
check against RESULTS-xxx-limit §X1c: one-site direction (1/10, 1/3, 1/10, 17/30, 3/10, 1, 7/15, 7/10) at ζ = (3/2,5/4,2/3) up to the
basis ordering — xxx_limit uses notebook order, gaudin_chain uses gl3_reps order with HW index 7: map carefully and say which order you use).
**Files you own:** `twovac_fit.py` (runnable top-to-bottom: `python -X utf8 twovac_fit.py`, prints all tables; stage flags to skip finished
stages), caches `twovac_fit_*.pkl`, `RESULTS-twovac-fit.md`. Environment: Python 3.13, sympy 1.14, mpmath 1.3, numpy; no `uv`. Time one
64×64 operation and one 60-digit eigen-sector before loops; cache states and Q's (P1's `qside_cache.pkl` may be loaded read-only).

## 0. Objects

Λ = (2,1,0) only. θ = (2/7, −3/5). The (∞,∞) vacuum Ω^∞[ζ,ω₁] = e₁ᵀρ(exp(r₁E₁₂)exp(r₂E₁₃)exp(r₃E₂₃)) with the r_k of the plan §1 (ω₂ → ∞
limits; check once numerically against `vacuum_covector_1site` at ω = (ω₁, 10⁶, ω₃)). Φ∞ := Ω^∞[ζ,ω₁]⊗Ω^∞[ζ,ω₁′]. States Ψ_B: all 64 at
each ζ-point (eigenstates of S₂(3/7) + (π/3)S₃(−2/9), 60 digits; any normalisation — everything below is normalisation-free). Dual Q's
P^j_B = e^{ζ_j u}p_j (monic, deg μ_j) and direct Q_i = e^{−ζ_i u}q_i (monic, deg 4−μ_i) from P1's solvers. The two targets:
Qnum(B) := ⟨Φ∞|Ψ_B⟩/⟨x₀^G|Ψ_B⟩ and T(B) := ⟨Φ∞|Ψ_B⟩·D_B(ω₀)/⟨Ω_D[ζ,ω₀]⊗Ω_D[ζ,ω₀]|Ψ_B⟩ with ω₀ = (5/2, 7/3, 3/5) and D_B(ω₀) = P1's det N
(the vacuum direct Q's of ω₀ in the Π⁻² gauge, the six Lagrange columns, rows (s,(i,j)) ∈ {1,2}×{(1,3),(1,2),(2,3)}). Check: T(B)/Qnum(B)
is state-independent (to precision) — if not, say so and use Qnum only. Jets at θ_s: J^j_k(θ_s) := (d/du)^k P^j(u)|_{θ_s}, k = 0..4
(exact from the polynomial coefficients and e^{ζ_jθ_s}); Plücker p_{kl}(θ_s) := J²_kJ³_l − J³_kJ²_l (0 ≤ k < l ≤ 3); lone P³^{(c)}(θ_s) :=
J³_c; direct jets Q₁^{(k)}(θ_s). Points: the plan §1 (ζ: P, P′, three extra; (ω₁, ω₁′) base (5/2, 7/4); the 3×3 grid; θ variants).

## 1. Tasks

**F0 — setup and sanity.** States and Q's at the five ζ-points (cache). ⟨x₀^G|Ψ_B⟩ ≠ 0 for all states (min |·| after G-normalisation or
max-component normalisation). Φ∞ sanity: per-weight ranks of the undone covector are not your job (W2), but check that ⟨Φ∞|Ψ_B⟩ ≠ 0 for
all 64 states at P (if some vanish, note which weights). T(B)/Qnum(B) state-independence. Print the 64 Qnum(B) at P with (S₂(3/7), S₃(−2/9))
labels to 30 digits (W2 prints the same table from the notebook: cross-check later, same labelling).

**F1 — the jet data.** For each state and site: J^j_k, p_{kl}, Q₁^{(k)}; verify (a) p₀₁(θ_s) = W(P²,P³)(θ_s) = 0 to precision and
p₀₂(θ_s) = W′(θ_s) = (ζ₃−ζ₂)e^{Σζθ_s}Π′(θ_s)Q₁(θ_s) (P1's Wronskian relation differentiated); (b) the Frobenius constraint: P^j′(θ_s) =
κ_s P^j(θ_s) with the SAME κ_s for j = 1,2,3 (report κ_s per state and derive its formula from the oper: with S₁ = e₁ + 3/(u−θ_s) + …,
S₂ = … + H_s/(u−θ_s) + 3/(u−θ_s)², S₃ = … + b_s/(u−θ_s) + a_s/(u−θ_s)² the indicial polynomial is ρ(ρ−2)(ρ−4) and the order-(ρ+1)
recursion for ρ = 0 gives c₁ in terms of the regular parts of S₁, S₂ and of H_s, a_s at θ_s — print the formula and check it on all 64
states); (c) the rank of the 3×4 jet matrix [J^j_k]_{j,k=0..3} (expect 2) and of the 3×5 (k = 0..4; expect 3). Record which jet monomials are
therefore linearly dependent across states.

**F2 — the fit (the main result).** Unknown kernel K, state-independent at fixed (ζ-independent!) (ω₁, ω₁′, θ). Hypothesis ladder, each a
linear least-squares problem over the 64 × 5 = 320 (state, ζ) rows at 80 digits (mpmath; QR or SVD — time a 320×144 case first):
- (a) target Qnum(B)·M¹(B)M²(B) with M_s(B) := Q₁(θ_s)P³(θ_s) [the Gaudin limit of M_vac up to constants], basis per site
  B_s^{dir} = {Q₁^{(k)}(θ_s)P³^{(c)}(θ_s) : k = 0,1,2; c = 0,1} (6) → 36 unknowns; also without the M¹M² factor, and the same two for T(B).
- (b) basis per site B_s^{Pl} = {p_{kl}(θ_s)P³^{(c)}(θ_s) : 0 ≤ k < l ≤ 3, c = 0,1} (12) → 144 unknowns, same four targets.
- (c) enlarge: c ∈ {0,1,2}; jets to order 4 in the Plücker set (10 p_{kl}); P¹ jets; (d) diagnostic only: let K depend linearly on the
  six Laurent data of the state (K = K₀ + Σ_m Δ_m K_m) — if THIS fits and nothing else does, the kernel is itself a Q-number (report).
For each level report: relative residual ‖A·k − y‖/‖y‖, rank of the design matrix A (SVD at 80 digits, tolerance 10⁻⁴⁰), null-space
dimension, and for the first level that fits (residual ≤ 10⁻³⁰): the minimal-norm kernel, a sparse representative (greedy: drop the
smallest entries while the residual stays ≤ 10⁻³⁰; or an L1-type selection), its support, and whether the entries are rational
(continued-fraction reconstruction at 50 digits with denominators ≤ 10¹²; they should be rational at rational (ζ, ω, θ) if the Gaudin K
tensor is rational). Controls that must FAIL: wrong lone index (P² instead of P³), jets at a wrong point (u = 0), swapped sites.

**F3 — the kernel laws.** Repeat the fitting level that worked at the 3×3 (ω₁, ω₁′) grid and at the two θ variants (with ζ = P, P′, +1
extra; 192 rows suffice once the structure is known). Checks: scaling symmetry (u → cu, g → g/c: the data at (θ, ζ, ω) and at (cθ, ζ/c, ω/c)
must give the same kernel up to the jets' rescaling — verify at c = 2 on one point); ζ-independence (fit per ζ-point separately and
compare). Fit every kernel entry as a rational function of low total degree in (ω₁, ω₁′) at fixed θ (exact rational data if F2 found
rationals: use sympy `rational_interpolate` on lines, or solve for numerator/denominator coefficients with a degree bound, and validate
on held-out grid points), then its θ₁₂-dependence. Report the laws in the XXX cell-62 style: {(c,c′), site-1 monomial, site-2 monomial,
κ(ω₁, ω₁′, θ₁₂)}; poles (where, in ω₁−ω₁′? in ω₁−ζ_j — there should be none since K is ζ-independent); the (ω₁, ω₁′) → (εω̂₁, εω̂₁′) scaling of
each entry (its ε-order).

**F4 — the regime and the Q-side CG (only if F2 fits).** ζ = εφ, ω₁ = εω̂₁, ω₁′ = εω̂₁′, grid ω̂₁ ∈ {3/4, −5/3}, ω̂₁′ ∈ {−12/5, 7/2},
ε ∈ {1,2,4}·10⁻⁶, 120 digits: the six HW-limit branches (P1's Q4 machinery: smallest raising residual per HW sector, G-norm 1 and the
round-1 sign convention), their P^j, jets, the formula's prediction vs the chain overlap ⟨Φ∞|Ψ_B⟩ (max relative deviation per ε);
ε-exponents of the jets, of the monomials, of the kernel entries, of the total (cancellation depth in digits, as in the XXX notes §5.6);
CG from the Q side alone: 1-dim slices and the 2-dim octet slices from two grid points, Richardson in ε, vs W1's exact values
−1/√2, −1/√2, 0.1178511302, octets (0.0493116000, 0.2088872897), (0.3379867346, −0.0797878449), mixing roots 5.854101966 / −0.854101966
(the untwisted copies are the eigenvectors of S3c[1,2] on the HW space; the φ_u(0) closed forms come from W2 — if its notebook is not yet
available, compute φ_u = e₁ᵀρ(n∞(εω̂₁)n∞(εω̂₁′)⁻¹) yourself: it is ε-free). Report the normalisation used.

**F5 — cross-checks.** If `RESULTS-twovac-limit.md` (agent D2) exists when you finish F3: compare its predicted Gaudin kernel/monomials with
yours (same point (5/2, 7/4), θ base) and report agreement or the discrepancy. If `Gaudin-L2-chain-2.wb` (agent W2) has its A5 table,
compare the 64 Qnum(B) at P.

## 2. Reporting

`RESULTS-twovac-fit.md`: §0 objects, orders, normalisations, the Frobenius κ_s formula; §F2 the hypothesis ladder table (level, target,
unknowns, rank, null dim, residual) and the kernel at the base point (dense and sparse representatives, rationality); §F3 the kernel laws;
§F4 the regime table and the CG values; §F5 cross-checks; surprises; not done. Numbers with precision, no adjectives. No notebooks, no git,
no edits outside your files. If F2 fails at every level, stop after the diagnostics of §8 of the plan and report the smallest ansatz that
fits (or that none does).
