# Spec P3: the Q side of the sl(3) Gaudin chain for Λ = (3,1,0), and the rank-4 puzzle of (2,1,0) (2026-10-16)

Delegated to a general-purpose subagent (sonnet), Python only (the Wolfram kernel belongs to another agent). Read first:
`Paul/Python/gaudin_sl3/RESULTS-Qside.md` (P1, round 1: the conventions, solvers, the Lagrange-identity columns, the det N and the ratio
identity for (2,1,0) — you are extending exactly this to (3,1,0)), `SPEC-Qside.md` (the original task statements Q1–Q4), and
`Paul/Mathematica/Gaudin/Gaudin-sl3-CG-plan.md` §1 (dictionary, ✔ = verified). Round-1 chain facts for (3,1,0) from the notebook summary
(`Paul/Mathematica/Gaudin/Experiments/Gaudin-L2-chain.wb`, last cell): same cdet, S₃ has no triple pole for (3,1,0) either, S_k(∞) = e_k(ζ),
ℓ₁ = 4; the multiplicity-2 HW sectors are (5,2,1) and (4,3,1).

**Code to reuse (read-only, import; never edit):** `gaudin_chain.py`, `qside.py`, `gl3_reps.py` (in `Paul/Python/`). If a function is
hard-wired to (2,1,0) (degrees 4−μ_i, pole orders, HW index 7, the Π⁻² gauge), copy it into your own module and generalise there; do not
change P1's files (another agent imports them concurrently). **Files you own:** `qside_310.py` (runnable top-to-bottom: `python -X utf8
qside_310.py`), `gaudin_chain_rep.py` (your rep-generic helpers, if needed), caches `qside_310_*.pkl`, `RESULTS-Qside-310.md`. Environment:
Python 3.13, sympy 1.14, mpmath 1.3, numpy; no `uv`. (3,1,0) has dim 15, 225 two-site states; the largest weight sectors are ~20–30-dim:
time one exact 225×225 sparse product and one 60-digit eigen-sector before loops; cache states and Q's.

## 0. The model

Λ = (3,1,0) on both sites, θ = (2/7, −3/5), g = diag(ζ), ζ = (3/2, 5/4, 2/3), bra twists ω = (5/2, 7/3, 3/5), ω′ = (7/4, −2/5, 9/7) — same as
round 1. Find the HW index in gl3_reps order (annihilated by E₁₂, E₁₃, E₂₃; report it). cdet as in round 1 (reuse `symbolic_cdet`); checks:
[S₂(3/7), S₃(−2/9)] = 0 exactly, S₁ scalar = e₁(ζ) + 4/(u−θ₁) + 4/(u−θ₂) (|Λ| = 4), weight conservation. Partial fractions of S₂, S₃
(exact fits, verified at extra points): the double-pole coefficient of S₂ (a scalar c₂(Λ); for (2,1,0) it was 3 — report the value and
whether it equals the quadratic Casimir-type number of Λ), the absence of a triple pole in S₃, the six non-scalar Laurent data
(H₁, H₂, a₁, a₂, b₁, b₂). Vacuum: Ω_{ζω} = e₁ᵀρ(n) with the same r_k (rep-independent formulas; ρ via the finite exponential series — the
nilpotency index is now 3), g′ = P_aux(ζ)g_c(ω)P_aux(ζ)⁻¹; verify Ω⊗Ω is the exact left eigencovector of all Laurent coefficients for ω and
ω′ and record the vacuum eigenvalue functions.

## 1. Tasks

**Q1 (3,1,0) — opers, degrees, the gauge power.** All 225 states (eigenproblem per weight sector, 60 digits). Dual: (∂³ − S₁∂² + S₂∂ − S₃)P = 0
with P^j = e^{+ζ_j u}p_j, expected deg p_j = μ_j. Direct: Q_i = e^{−ζ_i u}q_i with expected deg q_i = 2λ₁ − μ_i = 6 − μ_i, solving
Π^p∘(formal adjoint)∘Π^{−p} for the gauge power p to be FOUND: test p ∈ {0,1,2,3,4} (expect p = 2 or 3; round 1 fixed p = 2 for (2,1,0) by
the degrees and could not decide between λ₁ and λ₁−λ₃ … here both are 3 — if p = 3 works, propose p = λ₁ − λ₃ = λ₁ and say the (2,1,0) result
is consistent with p = λ₁; if p = 2 works, say p = λ₁ − λ₂ + 0? — report the smallest p giving all degrees 6 − μ_i with kernel dimension 1 for
each exponential). Local exponents of the dual oper at θ_s (indicial polynomial from the Laurent data; (2,1,0) had {0,2,4}; expect
{0, ?, ?} — report) and the implied Frobenius constraints on jets (for D1's benefit). Also the vacuum direct Q's for ω and ω′: degrees
(expected (0, 2λ₁−μ₂^{vac}, …) = (0, 4, 6)? — report what you find) and the vacuum dual degrees. Negative controls as in round 1.

**Q2 (3,1,0) — Wronskians and quantum Weyl.** Direct: Q_aQ_b′ − Q_bQ_a′ = c_ab e^{−Σζu} F(u) P^c with F expected ∏_s(u−θ_s)^{λ₁−λ₂} = Π²
(round 1: Π for the adjoint) — find F and c_ab (expect ζ_a − ζ_b with monic normalisation); dual Wronskians likewise; quantum-Weyl: S₁, S₂,
S₃ as ratios of Wronskian-type determinants of the dual Q's (round 1 formulas; W = det[P,P′,P″] = ∏_{i<j}(ζ_j−ζ_i)e^{Σζu}Π^{?} — report the
power). Residuals with precision.

**Q3 (3,1,0) — brackets, orthogonality, det N, the ratio identity.** Derive the Lagrange identity for the gauge power p found in Q1:
ψ = Π^{−p}Q_i (bra A), P = P^j (ket B): ψ·D_B[P] − P·D_B^†[ψ] = C′ ⇒ residues at θ_s of the six Laurent-difference terms give the columns
N[(s,i,j), H_t] = −Res_{θ_s}[ψP′/(u−θ_t)], N[·,a_t] = Res_{θ_s}[ψP/(u−θ_t)²], N[·,b_t] = Res_{θ_s}[ψP/(u−θ_t)] (the same six columns if the
Laurent structure is the same; if S₂'s double pole is still scalar and S₃ has no triple pole they are). Verify Σ_m Δ_m N = 0 over ≥ 300 random
pairs; the stacked 18×6 rank for A = B (6) and A ≠ B (≤ 5; tally the ranks); the 6×6 det N with rows (s,(i,j)), (i,j) ∈ {(1,3),(1,2),(2,3)}:
|det N| normalised, A ≠ B vs A = B; the x₀-free ratio identity R_B(ω)/R_B(ω′) = [D_B(ω)/D_B(ω′)]·[D_{B₀}(ω′)/D_{B₀}(ω)] on ALL 225 states
(max relative deviation, precision); the row-triple scan (84 triples: which pass, which are degenerate, which fail — round 1: the six
permutation triples failed). Controls that must fail as in round 1.

**Q4 (3,1,0) — regime (if time).** ζ = εφ, φ = (1/3, −2/7, 5/11), ω̂ = (3/4, −1/5, 1/2), ω̂′ = (−12/5, 1/8, −1/3), ε = 10⁻⁶ and 2·10⁻⁶ at
100 digits: the 10 HW-limit branches (sectors (6,2,0), (6,1,1), (5,3,0), (5,2,1)×2, (4,4,0), (4,3,1)×2, (4,2,2), (3,3,2); smallest raising
residual, report gaps), their Q degrees, the ratio identity on the branches, ε-exponents of det N.

**R — the rank-4 puzzle of (2,1,0) (cheap; do it first or last, your choice).** Reload P1's (2,1,0) data (`qside_cache.pkl` read-only or
recompute, ~80 s). For all ordered pairs A ≠ B (64·63), compute the rank of the stacked 18×6 matrix (σ-ratio threshold 10⁻⁴⁰ at 60 digits) and
tabulate by (weight of A, weight of B): which sector pairs give rank 4? Is it exactly the pairs where some column is identically zero for
the pair, or where ΔH₁ = ΔH₂ (etc.), or where the weights differ by a specific root? Find the second null vector n₂ (orthogonal to Δ) for
a few rank-4 pairs and test whether it is a fixed linear form in the Laurent data of A and B (e.g. the differences of the SCALAR data e₂,
e₃, c₂, or of the weights) — report the simplest explanation that fits all rank-4 pairs, or that none was found.

## 2. Reporting

`RESULTS-Qside-310.md`: §0 conventions for (3,1,0) (HW index, c₂, Laurent structure, the gauge power p and the rule it suggests), §Q1–Q4
tables (per-sector degrees, residuals with precision, F and c_ab, det N and the 225-state ratio test, the triple scan), §R the rank-4
explanation, surprises, not done. Scripts runnable top-to-bottom; no notebooks, no git, no Wolfram, no edits outside your files. If an
expectation is refuted and two attempts do not fix it, stop there, keep what works, report what you saw.
