# Spec G1: the general-Λ kernel of the Gaudin two-vacua jet formula — one-site kernel, more representations, structure, conjecture (2026-10-16)

Delegated to a general-purpose subagent (SONNET), Python only (the Wolfram kernel belongs to another agent). Read first:
`Paul/Mathematica/Gaudin/Gaudin-sl3-round4-plan.md` §0 (the formula and all notation — binding), §1 (your goal), §4, §6; then
`Paul/Python/gaudin_sl3/RESULTS-twovac-fit.md` (D1: the fit pipeline, normalisation by c and M, κ^G, the adjoint kernel), `RESULTS-Qside-320.md` §J
(P4: the [3,2,0] and [3,1,0] kernels, the non-forced Plücker basis, polynomial laws, site-exchange rule), `RESULTS-twovac-limit.md` (D2: the
adjoint kernel G as explicit polynomials, the ν-form), `RESULTS-gaudin-B.md` (D3: Θ = G/(4ω³ω′³) polynomial in (a, a′), rank 1 at 1/θ₁₂ = 0, the
jet basis), and `RESULTS-Qside-310.md`/`RESULTS-Qside.md` (P3/P1: rep-generic Q conventions, gauge p = λ₁−λ₃, exponents, κ formulas).

**Code to reuse (read-only; import or copy, never edit):** `twovac_fit.py` (D1: `jets`, `plucker`, `onesite_inf`, `cov_inf`, `x0_onesite`, `kappaG`,
`cscal`, `PivQR`, `design`, `fit_level`, `rat_recon`, `d2_G`), `twovac_fit_320.py` and its builders (P4: the rep-generic version — order-(λ₁+1) jets,
non-forced pairs, c for general Λ, `--stages`), `qside_320.py`/`qside_310.py`/`qside.py` (rep-generic Q solvers, `build_states(lam, zeta, dps)`),
`gaudin_chain.py`, `gl3_reps.py`. **Files you own:** `general_kernel.py` (stage flags; run heavy stages separately), caches `general_kernel_*.pkl`,
`RESULTS-general-kernel.md`. Environment: Python 3.13 (`python -X utf8`), sympy 1.14, mpmath 1.3, numpy; no `uv`. Three other agents run
concurrently (W4 on Wolfram — you must NOT use Wolfram; U1, U2 in Python on the untwisting limit — do not touch `untwisted_Q*`, `untwist_limit*`).
Keep to ≤ 2 heavy Python processes of your own at a time; time one eigen-sector and one QR before loops; cache states and Q's per rep.

## 1. Tasks

**K0 — canonical normalisation for the three known reps.** Define Θ := K/(ω₁ω₁′)^{2λ₁−λ₂} where K is the kernel of the diagonal-frame fit with
c = (ζ₂−ω₁)^{λ₁−λ₂}(ζ₃−ω₁)^{λ₁}(ζ₂−ω₁′)^{λ₁−λ₂}(ζ₃−ω₁′)^{λ₁}(ζ₁−ζ₂)^{2λ₂} and M = P³p_{0,λ₂+1} (so that ⟨Φ_c∞|Ψ_c⟩/⟨x₀^c|Ψ_c⟩·M₁M₂ = Σ Θ_ij m_im_j in the
companion frame). Adjoint: Θ = d2_G/(4ω³ω′³) (D3) — check that the "4" is D2's convention (M_vac = 2h²P³p₀₂ per site) and state your convention
once (recommended: M = P³p_{0,λ₂+1} exactly, no 4). [3,2,0]: from P4's grid kernel (`twovac_fit_320_grid.pkl`; rebuild from P4's stages if the pkl
format is unclear). [3,1,0]: redo P4's fit on a (ω₁, ω₁′) grid (P4 did the base point only) and get its polynomial laws. For each rep: verify Θ is a
POLYNOMIAL in (a, a′) = (1/ω₁, 1/ω₁′) and in 1/θ₁₂, the per-entry homogeneity (the XXX torus law: entry (kl),(k′l′) ∝ a^{e_kl}a′^{e_k′l′}·κ(s) with
s = a′/a? — find the exponent pattern), the site-exchange rule Θ_ij(a,a′,t) = Θ_ji(a′,a,−t) with t = 1/θ₁₂, and print the three Θ's in a common
labelling by Plücker pairs.

**K1 — the one-site kernel, exactly, for general Λ.** One site (L = 1, θ₁ = θ), twist ζ = diag(ζ), V_Λ: the eigenstates are the GT weight vectors
e_ν (the one-site Gaudin Hamiltonians S_k are Laurent in 1/(u−θ) with operator coefficients polynomial in the E's; verify they are diagonal in the
GT basis for diagonal ζ, or diagonalise exactly within weight spaces if a weight has multiplicity > 1 — it does not for gl(3) GT). One-site Q's:
P^j_ν = e^{ζ_j u}p_j with deg p_j = ν_j − λ₃ = ν_j (dual), direct Q_i = e^{−ζ_i u}q_i, deg λ₁ − ν_i; the one-site oper has ONE pole: solve exactly
(rational coefficients in ζ, θ) — Frobenius exponents {0, λ₂+1, λ₁+2} at θ as for two sites? (check). One-site vacuum Ω^∞[ζ,ω] (same r_k), one-site
x(ζ) = e₁ᵀP_phys^{G−1}, companion frame Ωc∞[ω] = e₁ᵀρ(e^{aE₁₂}e^{aE₂₃}), a = 1/ω, κ^G one-site = the same Δ₁^{λ₁−λ₂}Δ₂^{λ₂}. One-site jet formula:
⟨Ωc∞[a]|ψ_ν^c⟩/⟨e₁|ψ_ν^c⟩ · M(θ) = Σ_i θ_i(a) m_i(θ) with m_i = P³p_kl over the non-forced pairs and M = P³p_{0,λ₂+1} — LINEAR in the monomials;
fit θ_i(a) exactly (the data are exact rationals: use Fractions/sympy, no mpmath) for Λ ∈ {[1,0,0], [1,1,0], [2,0,0], [2,1,0], [2,2,0], [3,1,0], [3,2,0],
[4,1,0], [4,2,0]} (dims 3, 3, 6, 8, 6, 15, 15, 24, 27 — all cheap one-site); report whether the one-site formula holds at all (it should: D3 found Θ of
rank 1 at 1/θ₁₂ = 0), the θ_i(a) as polynomials in a (degrees, coefficients as functions of ζ? — expect ζ-independent after the κ^G normalisation),
and their closed form in (λ₁, λ₂, k, l): try monomials a^{e_kl}·c_kl(Λ) with c_kl rational functions of λ's; guess and verify on [4,2,0]. Also the
one-site sum rule with the top state e₁ (Q's explicit: P³ = e^{ζ₃u}, P² = e^{ζ₂u}p₂ with p₂ from the Wronskian relation, Q₁ = e^{−ζ₁u}).

**K2 — Θ for more two-site reps.** With the rep-generic fit pipeline (P4's `twovac_fit_320.py` generalised; states and Q's via the rep-generic
solvers), at canonical normalisation, on a (ω₁, ω₁′) grid sufficient for the polynomial laws (degrees ≤ 2·(λ₁+1)? — use P4's degree pattern to
size the grid) and ≥ 2 θ values: [1,1,0] (dim 3, 9 states, 2 monomials/site: pairs (0,2),(1,2)), [2,2,0] (dim 6, 36 states, 3 monomials/site:
(0,3),(1,3),(2,3)), [3,1,0] (if not already in K0), [4,1,0] (dim 24, 576 states, 14 monomials/site, 196 unknowns — budget: states at 50 digits,
rows 576 × 2 ζ; time first, skip the grid if a single fit exceeds 30 min and report the base point only). λ₂ = 0 probes [1,0,0] (dim 3) and [2,0,0]
(dim 6): first check the Laurent structure (P3: the S₂ double pole is NOT scalar for one-row reps), the exponents and which Plücker coordinates
are forced — if the oper/solver conventions break, say exactly where and stop those two. For each rep verify the fit's uniqueness (rank = number
of monomial pairs), rationality, polynomiality in (a, a′, 1/θ₁₂), and the site-exchange rule.

**K3 — structure across reps.** (a) Θ(1/θ₁₂ = 0) = θ(a)⊗θ(a′) with θ from K1 (exact comparison, all reps). (b) Θ = Σ_{n≥0} t^n Θ_n (t = 1/θ₁₂): for each
n the rank, the support pattern in Plücker pairs, homogeneity; Θ_1 in closed form across reps if visible. (c) The top-state sum rule as an exact
polynomial identity: Σ_ij Θ_ij m_i^{top}(θ₁)m_j^{top}(θ₂) = M^{top}(θ₁)M^{top}(θ₂) (companion frame: ⟨Φ_c∞|e₁⊗e₁⟩ = 1 = ⟨x₀^c|e₁⊗e₁⟩), with the top
state's Q's in closed form (P³ = e^{ζ₃u}, P² = e^{ζ₂u}p₂, Q₁ = e^{−ζ₁u}: from P²P³′ − P³P²′ = (ζ₃−ζ₂)e^{Σζu}Π^{λ₂}Q₁ one gets p₂′ − (ζ₃−ζ₂)p₂ =
−(ζ₃−ζ₂)Π^{λ₂}, so p₂ = Σ_{k≥0}(ζ₃−ζ₂)^{−k}(Π^{λ₂})^{(k)} — check the sign/normalisation against the monic convention and the actual top-state Q's).
Verify the identity symbolically for the known Θ's; it is a strong constraint for the conjecture. (d) Any determinantal structure: is Θ the
Laplace-expansion bilinear form of a 6×6 (or 2(λ₁+1)-type) determinant with site-s rows built from the jets? Test by checking Plücker-type
quadratic relations among the entries of Θ viewed as a matrix of minors (rank/quadric tests), and the XXX-style "products of 2×2 brackets"
decomposition count (rank of Θ_n).

**K4 — the conjecture.** Assemble the best general-Λ statement: ideally a closed form Θ_{(kl),(k′l′)}(a, a′, t; λ₁, λ₂); otherwise a structural
theorem-level description (θ(a) closed form, Θ₁ closed form, degrees/weights/zero pattern/rank rules, the sum rules) with explicit verification
on every computed rep and the predictions for the next rep ([4,2,0], [5,1,0]) spelled out so they can be tested later.

## 2. Reporting

`RESULTS-general-kernel.md`: K0 (the three Θ's in common labelling, polynomiality, homogeneity, exchange rule), K1 (the one-site kernels and
their closed form), K2 (new reps: fit diagnostics and Θ's), K3 (factorisation at t = 0, the t-expansion, the sum rule, determinantal tests), K4 (the
conjecture and its status), surprises, not done. Exact where possible; numbers with precision. No notebooks, no git, no Wolfram, no edits outside
your files. If an expectation is refuted and two attempts do not fix it, stop that item, keep what works, report what you saw.
