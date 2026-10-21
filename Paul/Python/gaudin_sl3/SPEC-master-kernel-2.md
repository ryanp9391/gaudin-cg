# Spec M3: the master kernel, second pass — exact [5,0,0]/[6,0,0] from the operator route, order 3, the one-variable z-family (round 5, WP-A follow-up; 2026-10-16)

Delegated to a general-purpose subagent (SONNET), Python only (no Wolfram, no wolfbook tools). Read first: `RESULTS-master-kernel.md` (M1: the
symbol form Q(ξ;η), the closed orders 0–2, the divisibility Q_n ∝ Π^{λ−n}, the coherent-state/Casimir picture, the translation covariance, the
[5,0,0] test, what was refuted — do not repeat the refuted attempts) and `RESULTS-jet-basis-sym.md` (M2: the operator-derived jet basis
⟨i,j|^c = ⟨x₀^c|ν_i(Ĥ,b̂)ν_j(Ĥ,b̂), the dictionary GT ↔ monomials, Ωc∞ = evaluation at (1,a,a²), the universal A₁, the [5,0,0] coefficients
`jet_basis_sym_Ycoef_ops5.pkl` = sympy matrices Y₀..Y₅ in the GT companion basis), then `Paul/Mathematica/Gaudin/Gaudin-sl3-round5-plan.md` §1
WP-A and `RESULTS-general-kernel.md` K3/K4 as needed. Environment: Python 3.13 (`python -X utf8`), sympy 1.14, mpmath 1.3, numpy; no `uv`; run from
`Paul/Python/gaudin_sl3/`. Other agents are running (S1, Q2 in Python; W5 on Wolfram): ≤ 2 heavy Python processes; time one unit before loops.

**Code to reuse (read-only; import or copy, never edit):** `master_kernel.py`/`master_kernel_bld.py` (M1: the R tables `master_kernel_tab.pkl`, the symbol
machinery, the closed orders as code, the 441-state [5,0,0] points `mk_*.pkl` at 40 digits — find their format), `jet_basis_sym.py` (M2: operator
construction of Y for one-row reps — reuse it to build Y for [6,0,0] (dim 28, 784 covectors) if needed; `t_op4.py`), `general_kernel.py` (`load_rep`,
`pairs_of`, `ckl_closed`, `polys_from_rep`), `gl3_reps.py`, `gaudin_chain.py`. **Files you own:** `master_kernel2.py`, `master_kernel2_*.pkl`,
`RESULTS-master-kernel-2.md` (write incrementally). No notebooks, no git, no Wolfram, no edits outside your files.

## 1. Tasks

**K5 — exact Θ^{[5,0,0]} and the acceptance test.** From M2's Y₀..Y₅: Θ^{[5,0,0]}(a,a′,t) = ev_{a,a′}·Y(t)⁻¹ (ev = the evaluation functional at
((1,a,a²),(1,a′,a′²)) in M2's GT companion dictionary; one 441×441 exact rational solve per rational point (a, a′, t) — time one; use Fraction/
sympy or flint-free exact arithmetic; if a point takes > 10 min use 80-digit mpmath + rational reconstruction). (a) Acceptance test of the
round-5 plan: evaluate this Θ at M1's 441-state base point(s) and print the per-state residual of c·Qnum·M₁M₂ = Σ Θ_ij m_i m_j (expect ~1e-30 at
M1's 40 digits) — this closes WP-A's "predict [5,0,0] on ≥ 200 states" item by the operator route. (b) Reconstruct Θ^{[5,0,0]} as exact
polynomials in (a, a′, t) on enough rational points (M1's support rule bounds the monomials: a-degree ≤ 2λ₁, t-degree ≤ 4(λ₁−1) = 16, with the
homogeneity p+q−n = e_i+e_j — count the unknowns per entry and size the grid; exploit the exchange rule and the translation covariance to cut the
grid), and store it as `general_kernel`-format `T[(i,j)]` in `master_kernel2_rep_500.pkl`. Also R^{5} in M1's table format.

**K6 — order 3 (and 4) closed.** With λ = 2, 3, 4, 5 exact (and λ = 6 via the operator route if the λ-polynomial degree demands it — M2's code builds
Y for [6,0,0]: 784×784, degree 6 in t; time it), extend M1's symbol form: Q₃/Π^{λ−3} as a polynomial in λ (determine its degree from the data —
M1 expects degree 6 in λ at order 3 — and from the structure Q_n = λ²(λ−1)…(λ−n+1)Π^{λ−n}(polynomial in K, A-type invariants); use the invariants M1
found (K, A, e₁−h₁, e₁h₁−e₂−h₂, the third-point values) as the basis of an ansatz and fit the λ-dependence), verify on every λ you have, and state
the recursion/pattern across orders (Q₁ = λ²Π^{λ−1}K and Q₂ suggest Q_n = Π^{λ−n}·(λ)_n-type prefactors × a universal polynomial in K and in the
"connection" invariants — find the generating function if the pattern is visible, e.g. Σ_n Q_n z^n = Π^λ·exp(...) or a hypergeometric-type series).

**K7 — the one-variable z-family and the Casimir completion.** M1: translation covariance + scaling reduce Θ to a matrix family in z = t/(ω−ω′)
only, and Θ(z) on the coherent state v = ℓ_a^λ⊗ℓ_{a′}^λ is NOT a function of the Casimir Ω alone — each irrep component of V_λ⊗V_λ = ⊕_μ V_μ needs a
second vector. Find that second operator: candidates are M2's universal A₁ = (E₂₂+E₃₃)∧E₂₁ + E₂₃∧E₃₁ + E₃₃∧E₃₂ (the first-order coupling of the jet
basis), the sl(2) generators N, F of M1's A3 acting at the third point (the vacuum at ω = ∞), and the "difference" operators E⁽¹⁾_ij − E⁽²⁾_ij.
Test: is Θ(z)v in the span of {Ω^k v, Ω^k A v, …} for a single extra operator A (compute the span dimensions per irrep for λ = 2, 3 — M1 found
component spans 1,2,2 and 1,3,5,3; a correct A must reproduce exactly those)? If yes, write Θ(z)v = Σ_μ P_μ(z)·(projection onto V_μ of a
2-term combination) and read off the z-dependence per irrep (expect rational functions of z with λ-dependent integer coefficients — the
"R-matrix" of the fusion). This is the structural closed form if it exists; report the span test honestly if it does not.

## 2. Reporting

`RESULTS-master-kernel-2.md`: §0 conventions, §1 K5 (the acceptance residuals, the exact [5,0,0] kernel and its sizes), §2 K6 (order 3(,4) closed
forms with verification per λ; the pattern/generating function if found), §3 K7 (span tests, the second operator, the per-irrep z-functions),
surprises, not done. Exact where possible.
