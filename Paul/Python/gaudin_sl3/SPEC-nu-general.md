# Spec N1: the ν-formulas for general Λ and the untwisting limit for [3,1,0] and [3,2,0] (round 5, WP-B1 + B2; 2026-10-16)

Delegated to a general-purpose subagent (SONNET), Python only (no Wolfram, no wolfbook tools). Read first: `Paul/Mathematica/Gaudin/
Gaudin-sl3-round5-plan.md` §0 (P1, P3) and §1 WP-B (your goal: B1 and B2; B3 — the EXACT ε-series — belongs to agent W5 on Wolfram, running
concurrently, who writes `Paul/Mathematica/Gaudin/ref/laurent-series-310.txt`/`-320.txt`; use them for a cross-check if they exist when you get
there, never wait for them), `Paul/Mathematica/Gaudin/Gaudin-sl3-round4-plan.md` §0 (notation — binding), the notes `Paul/Mathematica/Gaudin/
Gaudin-sl3-CG-notes.md` §2 (opers, exponents, no-log conditions), §7 (the Frobenius reduction — your B1 generalises 7.1 exactly as 7.2 anticipates),
§9 (the untwisting limit — your B2 repeats 9.2–9.5 for two more reps), §10–11, Appendix D (ε-series, CG tables), and `RESULTS-untwist-limit.md`
(U2: the pipeline you generalise — read it fully), `RESULTS-untwisted-Q.md` (U1), `RESULTS-Qside-310.md`/`RESULTS-Qside-320.md` (P3/P4: rep-generic
Q conventions, Frobenius constants, the D₂ relations), `RESULTS-general-kernel.md` §0 (kernel normalisation). Environment: Python 3.13
(`python -X utf8`), sympy 1.14, mpmath 1.3, numpy; no `uv`; run from `Paul/Python/gaudin_sl3/`. Four other agents run concurrently (M1, M2, Q1 in
Python; W5 on Wolfram): ≤ 2 heavy Python processes of your own at a time; time one unit before loops; cache to `nu_general_*.pkl`.

**Code to reuse (read-only; import or copy, never edit):** `untwist_limit.py` (U2: `build_regime`, `branch_N`, `series_fit`, `gauss_closed`, `chat_expr`,
`gf_polys`, `branch_psi`, `comp_series`, `lead_exp`, `G_pieces`, `N_series`, `leading_poly`, `nu_funcs`, `nu_from_laurent`, `nu_series`, `analyse_nu`),
`untwist_limit_frob.py` (U2's symbolic Frobenius for exponents {0,2,4} — your template for B1), `untwisted_Q.py` (U1: `frobenius_conditions`,
`nolog_values`, `hw_spaces`, `raising_residual`, `oper_coeffs`, `dual_poly`), `twovac_fit.py` (`build_regime`, `hw_space`, `gram_diag`, `gnormalise`,
`S2S3_labels`), `general_kernel.py` (`load_rep(lam)` → exact Θ^{[Λ]} as `T[(i,j)] = {(n,p,q): Fraction}` with Θ_ij = Σ c tⁿ aᵖ a′^q, i,j over `pairs_of(lam)`;
`eval_T`, `e_of`, `S_at`, `site_monos`, `get_point2`), `qside_310.py`/`qside_320.py` (`build_states(lam, zeta, dps)`, `make_S`, `frobenius_polys`,
`frobenius_check`, `local_exponents`, `kappa_checks`, `raising_residual`), `gaudin_chain.py`, `gl3_reps.py`. W4's exact chain-side data:
`Paul/Mathematica/Gaudin/ref/untwisted-exact-310.txt`, `-320.txt` (TAB-separated InputForm; header documents every record: GTBASIS order, HW sectors,
HWBASIS kets, GRAM, LAURMAT, BRANCH, XVEC = x(ζ) as a polynomial covector) and `-210.txt` (also E4/E5 records for the adjoint). **Files you own:**
`nu_general.py`, `nu_general_*.pkl`, `RESULTS-nu-general.md` (write incrementally). No notebooks, no git, no Wolfram, no edits outside your files.
If an expectation is refuted and two attempts do not fix it, stop that item, keep what works, report what you saw.

## 0. What is established and what you must show

Dual oper (∂³ − S₁∂² + S₂∂ − S₃)P = 0, P^j = e^{ζ_j u}p_j. As rational functions of u (no triple poles: the indicial polynomial at θ_s is
ρ(ρ−1)(ρ−2) − |Λ|ρ(ρ−1) + c₂ρ − c₃ with roots {0, λ₂+1, λ₁+2}, forcing c₃ = 0 and c₂ = λ₂(λ₁+1) — verify numerically that S₂, S₃ have no higher
poles for [3,1,0], [3,2,0] with P3/P4's `make_S`):
S₁ = σ₁ + |Λ|Σ_s 1/(u−θ_s), S₂ = e₂(ζ) + Σ_s[H_s/(u−θ_s) + c₂/(u−θ_s)²], S₃ = e₃(ζ) + Σ_s[b_s/(u−θ_s) + a_s/(u−θ_s)²], σ₁, e₂, e₃ = elementary symmetric
functions of ζ. (Verify this partial-fraction form exactly on the states — if the constants are not e₂, e₃, say what they are.) Frobenius at θ_s:
every solution P = Σ c_m x^m (x = u−θ_s) has c_m determined by the recursion for m ∉ {λ₂+1, λ₁+2}; c₀ and c_{λ₂+1} are free; the recursion
equations at the resonant orders are the no-log conditions (3 per site). Hence c_k = c₀u_k + c_{λ₂+1}v_k with u_k, v_k RATIONAL in the Laurent data
(H_s, a_s, b_s; H_t, a_t, b_t through the Taylor coefficients of the other site's poles at θ_s; σ₁, e₂, e₃; δ = θ_s − θ_t), v_k = 0 for k ≤ λ₂,
v_{λ₂+1} = 1, u₀ = 1. Therefore p_kl = P²^{(k)}P³^{(l)} − P³^{(k)}P²^{(l)} = D_s·k!l!(u_kv_l − u_lv_k) with D_s = c₀^{(2)}c_{λ₂+1}^{(3)} − c₀^{(3)}c_{λ₂+1}^{(2)}, M = P³p_{0,λ₂+1} =
P³D_s(λ₂+1)!, and in the jet formula c·Qnum·M₁M₂ = Σ Θ_ij m_i m_j the factor P³D_s cancels: **c·Qnum = Σ_ij Θ_ij(a,a′,t)·(ω₁ω₁′)^{2λ₁−λ₂}·ν_i(θ₁)ν_j(θ₂)
with ν_kl := m_kl/M = k!l!(u_kv_l − u_lv_k)/(λ₂+1)!, ν_{0,λ₂+1} = 1** — a function of the Laurent data only. For the adjoint this is U2's
ν = (1, κ, π, κπ, ε_s). Your B1 makes this explicit for general (λ₁, λ₂) and verifies it; your B2 uses it.

## 1. Tasks

**B1 — ν for general Λ (sympy, exact).** Write `nu_symbolic(lam1, lam2)` → the list of ν_kl (non-forced pairs in `pairs_of` order) as rational
functions of the symbols (H_s, a_s, b_s, H_t, a_t, b_t, σ₁, e₂, e₃, δ), plus the no-log conditions (the 3 resonant-order equations per site, as
polynomials; the adjoint's b_s = a_s(H_s − a_s)/3 must come out of the order-(λ₂+1) one; for [3,1,0] expect b_s = a_s(H_s − a_s)/4 (P3's D₂ relation
S3c[s,1] = ¼S3c[s,2](H_s − S3c[s,2])); for [3,2,0] P4 found an inhomogeneous quadratic — reproduce it) and `nu_numeric(lam, laurent_s, laurent_t,
zeta, theta_s, theta_t)` (mpmath). Verify on actual states: for [3,1,0] and [3,2,0] at G1's cached base points (`general_kernel_pt_310_P_*`, `_320_P_*`;
`site_monos(lam, pt)` gives m_i; the Laurent data per state via `S_at`/P3–P4's `make_S` or by fitting the partial fractions of the S₂, S₃ eigenvalues)
check m_i/M = ν_i to ≥ 50 digits on every state (225 each), and also for the adjoint against U2's `nu_funcs`. Report the general formulas for the
first few ν's (the analogues of κ = a_s/c₂, π, ε_s) in readable form and the no-log conditions in closed form for both reps (and for general λ₂ if
the pattern is clear: the order-(λ₂+1) condition is a polynomial in (H_s, a_s, b_s) only when λ₂ ≤ 1 and involves the other site for λ₂ ≥ 2 — say
exactly what enters).

**B2 — the untwisting limit in ν-form for [3,1,0] and [3,2,0].** Regime as U2: θ = (2/7, −3/5), ζ = εφ with φ = (1/3, −2/7, 5/11) and φ′ = (−2/5, 3/8, 7/11),
ω₁ = εω̂₁, ω₁′ = εω̂₁′ on U2's 16-point grid, ε_k = 10⁻⁶·2^{−k}, k = 0..11, 200 digits. For each ε and φ: the 10 HW-limit branches of the rep (HW
sectors from W4's files; diagonalise PER WEIGHT SPACE — the twisted Hamiltonians commute with the total Cartan, so never diagonalise the full
225×225; identify the branch by the smallest raising residual in the weight-μ space, the two copies of a multiplicity-2 sector by continuity in ε
and by the ε⁰ value of a_s against W4's LAURMAT/BRANCH records — state your copy convention), their Laurent data (H_s, a_s, b_s per site: eigenvalues
of the Laurent-coefficient operators — build them once per ε from `gaudin_chain.Operators`/P3–P4's code), ν via B1, and
Σ_ij G_ij(εω̂₁, εω̂₁′, t)ν_i(θ₁)ν_j(θ₂) with G = (ω₁ω₁′)^{2λ₁−λ₂}Θ^{[Λ]} from `load_rep`. Then: (i) the leading power ε^{N_B} per branch and the leading
coefficient as a polynomial L_B(ω̂₁, ω̂₁′) (fit on the grid; check the lower orders cancel); (ii) the chain side: Ψ_B(ε) G-normalised (sign: first
nonzero GT component > 0), its ε-series Ψ₀ + εΨ₁ + …, x(ζ) from W4's XVEC (degree-λ₂ polynomial covector — parse the record; check x(ζ)⊗x(ζ) = x₀^G
against the code's `x_cov`), k_B := first n with X_B := ⟨x(φ)⊗x(φ)|Ψ_n⟩ ≠ 0, ĉ := c/ε^{4λ₁} (c of the round-4 plan §0), φ_u = e₁ᵀρ(n′) with U2's ε-free
Gauss parameters (`gauss_closed`; ρ for the rep via `gl3_reps`/the chain's generators — nilpotent exponentials are finite sums); (iii) the identity:
is L_B(ω̂) = C_Λ,B·ĉ·⟨φ_u⊗e₁|Ψ₀⟩/X_B as polynomials in (ω̂₁, ω̂₁′), with what constant C_Λ,B (the adjoint had 1 in U2's normalisation; with
ν_{0,λ₂+1} = 1 expect a rep-dependent but branch-INDEPENDENT constant — report it per branch) and what exponent N_B (adjoint: 6 − k_B; look for the
rule, e.g. N_B = 4λ₁ − 2 − k_B or in terms of the weights)? (iv) the CG dictionary: expand ⟨φ_u⊗e₁|Ψ₀⟩ in the Gauss parameters and read off the
CG coefficients C_ν = ⟨e_ν⊗e₁|Ψ₀⟩ for every site-1 GT label ν; compare with W4's exact HW kets (HWBASIS + GRAM) — all 10 + 10 tables to ≥ 30 digits;
(v) φ-independence (φ vs φ′); (vi) the depth of the ν-series needed per branch (as U2's table) — this tells W5/S1 how far the exact series must go.
Cross-check your numerical ε-series of (H_s, a_s) (from `series_fit` over the 12 ε's) against `ref/laurent-series-<Λ>.txt` if W5 has written it.
Budget: build the regime points for [3,1,0] first (both φ), analyse, then [3,2,0]; cache every regime point; if the full grid is too slow, 8 grid
points and k ≤ 9 still determine the leading polynomials (degrees ≤ 2λ₁ per variable — check).

## 2. Reporting

`RESULTS-nu-general.md`: §0 conventions, §1 B1 (the general ν's, the no-log conditions, verification residuals per rep), §2 B2 per rep (branches,
k_B, N_B, the constants, the identity residuals, the 20 CG tables vs exact, φ-independence, depths), surprises, not done. Exact where possible;
numbers with precision. Give the rule for N_B and k_B if you see one, and say clearly whether the ν-reduction holds for both reps.
