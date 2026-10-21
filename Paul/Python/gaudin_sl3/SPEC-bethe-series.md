# Spec S1: the ε-series of the Laurent data from the Bethe side alone (the degenerate solver), and the end-to-end script `cg_from_Q.py` (round 5, WP-C1 + F3; 2026-10-16)

Delegated to a general-purpose subagent (SONNET), Python only (no Wolfram, no wolfbook tools). Read first: `Paul/Mathematica/Gaudin/
Gaudin-sl3-round5-plan.md` §0 (P3) and §1 WP-C1, WP-F3 (your goals), `Paul/Mathematica/Gaudin/Gaudin-sl3-round4-plan.md` §0 (notation — binding),
the notes `Paul/Mathematica/Gaudin/Gaudin-sl3-CG-notes.md` §2.3 (no-log conditions), §9 (the limit), §10 (the series and how they were obtained —
10.2 describes the rank-5 failure you fix), Appendix D, and `RESULTS-untwisted-Q.md` (U1: `perturb_data`, `perturb_data_deg`, the rules used, where
the degenerate case failed), `RESULTS-untwist-limit.md` (U2: the identity and the CG dictionary — F3 packages it), and whatever exists of
`RESULTS-nu-general.md` (N1: general-Λ no-log conditions) and `Paul/Mathematica/Gaudin/ref/laurent-series-310.txt`/`-320.txt` (W5: exact series
for the other reps — your cross-check targets; if absent, do the adjoint fully and the other reps as far as your own derivation allows).
Environment: Python 3.13 (`python -X utf8`), sympy 1.14, mpmath 1.3, numpy; no `uv`; run from `Paul/Python/gaudin_sl3/`. Other agents may still be
running (M1, M2, N1, Q1 in Python; W5 on Wolfram): ≤ 2 heavy Python processes; cache to `bethe_series_*.pkl`.

**Code to reuse (read-only; import or copy, never edit):** `untwisted_Q.py` (U1: `frobenius_conditions`, `nolog_values`, `zeta0_data`, `perturb_data`,
`perturb_data_deg`, `casimirs_total`, the pkl's `untwisted_Q_laurent_exact.pkl` (n ≤ 5), `untwisted_Q_pert_sym.pkl` (n ≤ 6), `untwisted_Q_octet_fit.pkl`),
`untwist_limit.py` (U2: `nu_funcs`, `nu_from_laurent`, `G_pieces`, `gauss_closed`, `chat_expr`, `gf_polys`, `leading_poly`, `analyse_nu`), `general_kernel.py`
(`load_rep`, `ckl_closed`, `theta1_closed`, `theta_single_closed`), `twovac_fit.py` (`d2_G`, `hw_space`, `gram_diag`), `gl3_reps.py`, `gaudin_chain.py`;
W4's `Paul/Mathematica/Gaudin/ref/untwisted-exact-210.txt` (E4_SERIES, E5_KB = X_B closed forms, E5_GF, HWBASIS/GRAM = exact CG tables). **Files you
own:** `bethe_series.py`, `cg_from_Q.py`, `bethe_series_*.pkl`, `RESULTS-bethe-series.md` (write incrementally). No notebooks, no git, no Wolfram, no
edits outside your files.

## 1. Tasks

**C1 — the degenerate order-by-order solver (sympy, φ symbolic).** Unknowns: the six Laurent data x = (H₁, a₁, b₁, H₂, a₂, b₂) as series in ε with
ζ = εφ, θ₁₂ fixed (= 1 in units; restore by scaling). Equations: the six no-log conditions F(x; ζ, θ) = 0 (U1's `frobenius_conditions` with the
regular parts expressed through the other site's data — three per site; adjoint exponents {0,2,4}) plus the two global rules H₁ + H₂ =
Σ_{i≠j}ζ_iμ_j and b₁ + b₂ = ε²(φ₁φ₂μ₃ + φ₁φ₃μ₂ + φ₂φ₃μ₁) (verified by U1/W4; state which equations are independent). Order by order: J·x_n = r_n(x_0..x_{n−1})
with J = ∂F/∂x at ε = 0 on the branch. U1's solver fails where rank J = 5 (top (4,2,0) and the octets (3,2,1)): implement the Fredholm
alternative — the kernel direction k of J is left undetermined at order n and fixed by the solvability condition ℓ·r_{n+1} = 0 of the next order
(ℓ = the left null vector; if the solvability condition is itself degenerate, go one order further — second-order degenerate perturbation
theory; for the octets the two copies are the two roots of the resulting quadratic — the Galois conjugates over Q(√5)). Deliverable: exact series
for ALL six adjoint branches to n = 8 from these equations alone (no chain data except the ε⁰ untwisted values — and even those should follow
from the equations at order 0: check), agreeing with U1/W4 (n ≤ 7) exactly; the top branch's exact linear form and the sqrt forms of (4,1,1),
(3,3,0) recovered as resummed series. Then [3,1,0] and [3,2,0]: generalise the no-log conditions (exponents {0, λ₂+1, λ₁+2}; N1's `nu_general.py` may
already have them — import if present, else derive with the same Frobenius code) and the global rules (the first is general; the second: derive
from a₁ + a₂ + θ₁b₁ + θ₂b₂ = ? and b₁ + b₂ = ? at finite ζ — W4's untwisted sum rules are the ε⁰ case; find the twisted versions by fitting on the
chain if necessary — `gaudin_chain` is allowed for THIS calibration only), and produce the 10 + 10 series to n = 5; compare with W5's data files.
Report: the mechanism (why rank 5 exactly on the top and the octets — which combination of the six conditions is dependent at ε = 0 and what the
solvability condition says), and whether the Bethe side alone determines everything (count of free parameters at each order).

**F3 — `cg_from_Q.py`.** One file, runs in < 5 min, prints the adj⊗adj CG table and its comparison with the exact values. Inputs, all from
closed forms or the committed data (no chain diagonalisation, no fits): the kernel G^{[2,1,0]} (from `general_kernel_rep_210.pkl` or D2's `d2_G`;
also show the one-site closed form `ckl_closed` + block theorem reproduce its t = 0 part); the ν-formulas (`nu_funcs`); the exact ε-series of
(H_s, a_s) from C1 (fallback: `untwisted_Q_laurent_exact.pkl`) to the depth U2 lists (3, 5, 5, 4, 4, 2); the Gauss parametrisation of φ_u and ĉ
(`gauss_closed`, `chat_expr`); X_B closed forms (W4's E5_KB); U2's identity Σ G(εω̂)νν = ε^{6−k_B}ĉ⟨φ_u⊗e₁|Ψ₀⟩/X_B read as a polynomial identity in
(ω̂₁, ω̂₁′) → the coefficients ⟨e_ν⊗e₁|Ψ₀⟩; the normalisation from `RESULTS-gaudin-norm.md` if Q1 has delivered a Q-side norm (else normalise by
the top coefficient and say so). Output: the six branches' CG slices vs W4's exact table, digits of agreement, a one-paragraph docstring
explaining each step with the equation numbers of the notes. Keep imports to the existing modules (read-only) and the standard libraries.

## 2. Reporting

`RESULTS-bethe-series.md`: §0 the equations and their ranks, §1 the degenerate solver and the adjoint series (n = 8), §2 the other reps, §3 `cg_from_Q.py`
(what it reproduces, timing), surprises, not done.
