# Spec Q1: the Gaudin norm from the Q side — ⟨Ψ|Ψ⟩_G/⟨x₀|Ψ⟩² as a function of the Laurent data / jets, and its untwisting limit (round 5, WP-D; 2026-10-16)

Delegated to a general-purpose subagent (SONNET), Python only (no Wolfram, no wolfbook tools). Read first: `Paul/Mathematica/Gaudin/
Gaudin-sl3-round5-plan.md` §0 and §1 WP-D (your goal) and §4 (risk: if none of the candidates is state-independent, flag and stop rather than
invent), `Paul/Mathematica/Gaudin/Gaudin-sl3-round4-plan.md` §0 (notation — binding), the notes `Paul/Mathematica/Gaudin/Gaudin-sl3-CG-notes.md`
§2.3–2.5 (no-log conditions, Wronskians, detN and orthogonality), §6.4–6.5 (the jet basis and the K tensor), §7 (the ν-form), §9.4–9.5, and
`RESULTS-Qside.md` §3 (P1: det N with the Lagrange columns, orthogonality), `RESULTS-gaudin-B.md` (D3: the jet matrix Y, its caches, the remark
"Qnum_eq·M₁M₂/det N_B state-independent"), `RESULTS-untwist-limit.md` (U2: the regime caches and the ν-form), `RESULTS-untwisted-Q.md` (U1: the
no-log conditions as code). Environment: Python 3.13 (`python -X utf8`), sympy 1.14, mpmath 1.3, numpy; no `uv`; run from `Paul/Python/gaudin_sl3/`.
Four other agents run concurrently (M1, M2, N1 in Python; W5 on Wolfram): ≤ 2 heavy Python processes of your own; time one unit before loops;
cache to `gaudin_norm_*.pkl`.

**Code to reuse (read-only; import or copy, never edit):** `twovac_fit.py` (`get_point('P'|'Pp'|'X1')` → the 64 adjoint states with dual Q's and jets at
ζ = P etc.; `gram_diag(ch)`/`gnormalise` (the invariant diagonal metric G per site; two-site metric G⊗G; GT matrices are NOT Hermitian, E_ijᵀG = GE_ji),
`x0_D2`/`x0_onesite`, `detN_for_states(pt, om0, th)`/`detdata`, `S2S3_labels` (Laurent data per state), `onesite_inf`/`cov_inf` (the vacua),
`hw_space`, `build_regime` and U2's caches `untwist_limit_reg_<k>.pkl`/`reg2` (12 regime points, 200 digits)), `untwisted_Q.py` (`frobenius_conditions()`
→ the symbolic no-log conditions in the variables (H, a, b, s₀..s₂, h₀..h₂, t₀..t₂) — the Taylor coefficients of the regular parts; `nolog_values(data,
theta, zeta)` evaluates them on a state), `untwist_limit.py` (`nu_funcs`, `nu_from_laurent`, `G_pieces`), `gaudin_B.py` (D3: `jet_matrix(tag)`,
`comp_jet_matrix`, the pkl's `gaudin_B_jets_P.pkl`, `gaudin_B_Ycomp_P.pkl`, `gaudin_B_cov*.pkl` — see RESULTS-gaudin-B for formats; `omega_comp_sym`,
`d2_G_sym`), `qside.py` (`q2_weyl`, `q2_wronskians`, `Nrows`, `rel_det`: P1's det N machinery), `gaudin_chain.py` (`Chain`, `Operators`, `Pmat`, `gc`,
`gprime`), `gl3_reps.py`. **Files you own:** `gaudin_norm.py`, `gaudin_norm_*.pkl`, `RESULTS-gaudin-norm.md` (write incrementally). No notebooks, no git,
no Wolfram, no edits outside your files. Adjoint only ([2,1,0]²) unless everything closes early.

## 0. The question

The CG generating function from the Q side (U2) is ⟨φ_u⊗e₁|Ψ_B⟩/⟨x₀^G|Ψ_B⟩-type RATIOS; the CG table needs Ψ₀ with G-norm 1, i.e. the functional
N(Ψ) := ⟨Ψ|Ψ⟩_G/⟨x₀^G|Ψ⟩² (⟨Ψ|Ψ⟩_G = Ψᵀ(G⊗G)Ψ; the Hamiltonians are G-symmetric, so the left eigenvector is (G⊗G)Ψ and the overlap is
non-degenerate). Find N(Ψ) from Q-data (Laurent data / jets) with no chain input — the Gaudin analogue of the Sklyanin-measure norm formula /
the Yang–Yang Hessian. Test every candidate on all 64 adjoint eigenstates at ζ = P, Pp, X1 (θ = (2/7, −3/5)); a candidate passes if
N(Ψ)/candidate(Ψ) is state-independent to ≥ 40 digits (a ζ, θ-dependent constant is fine; find it).

## 1. Tasks

**D0 — ground truth.** N(Ψ) for the 64 states at the three twists; also with ⟨Ω₀⊗Ω₀|Ψ⟩ (the twisted vacuum, `cov_inf` at ω = ω₀) in place of x₀^G.
Record the Laurent data, the monomials m_i(θ_s), M_s, and det N_B per state.

**D1 — candidates (in this order; stop at the first that passes, but ALWAYS do (iv)).**
(i) det N_B(ω) with P1's Lagrange columns, and D3's quantity Qnum_eq·M₁M₂/det N_B (read D3 §B4 for what exactly was state-independent: the EQUAL-twist
overlap over det N — that is an overlap identity, not yet a norm; test whether ⟨Ψ|Ψ⟩_G relates to det N_B·(overlaps) — e.g. N(Ψ) ∝ det N_B/(⟨x₀|Ψ⟩-type
products) or ∝ 1/det N_B with the vacuum normalisation).
(ii) The Jacobian of the no-log conditions: F_k(H₁,a₁,b₁,H₂,a₂,b₂; ζ, θ) = 0, k = 1..6 (U1's `frobenius_conditions` evaluated with the regular parts
expressed through the OTHER site's data — write the six conditions as explicit functions of the six Laurent unknowns); candidates: det(∂F/∂(H,a,b)₁,₂)
(6×6), and the reduced 4×4 after eliminating b_s = a_s(H_s − a_s)/3 (then 4 unknowns, 4 remaining conditions); times simple factors (powers of
M_s, of P³(θ_s), of the Wronskian W ∝ e^{Σζu}Π³, of θ₁₂, of ⟨x₀|Ψ⟩). This is the "Yang–Yang Hessian in Laurent variables" hypothesis.
(iii) Wronskian/quantum-Weyl determinants: products Π_s (P³D_s)(θ_s) (D_s = the 2×2 Frobenius determinant = M_s/(P³·2!)), the values of the
direct Q's at θ_s, and the Weyl-type determinants of P1's `q2_weyl`.
(iv) THE SoV MEASURE (exact by construction — the real deliverable): in the companion frame Ψ_c = ⟨x₀^c|Ψ_c⟩ Σ_ij ν_i(θ₁)ν_j(θ₂)|i,j⟩ where |i,j⟩ is
the ket basis dual to D3's jet covectors ⟨i,j| (columns of Y⁻¹) and ν_i = m_i/M (U2's ν = (1,κ,π,κπ,ε_s) extended to all 8 monomials per site —
the 8 one-site monomials m = (P³p₀₂, P³′p₀₂, P³p₀₃, P³′p₀₃, P³′p₁₂, P³p₂₃, P³′p₁₃, P³′p₂₃) of D3 are ALL expressible through the Laurent data by the
Frobenius argument: do it, 8 ν's per site). Then ⟨Ψ|Ψ⟩_G = ⟨x₀|Ψ⟩²·Σ ν_iν_j Γ_{(ij),(i′j′)} ν_{i′}ν_{j′} with Γ = (Y⁻¹)ᵀ·G_c·Y⁻¹, G_c = the metric G⊗G
transported to the companion frame (P_phys^G(ζ)^{⊗2}; D3's `comp_jet_matrix` tells you the frame). Compute Γ exactly (64×64 rational in ζ, θ₁₂ —
rational-reconstruct at the three twists or keep 80-digit numerics), verify the norm identity on all 64 states at all three twists (this MUST hold;
if it does not, your frame/normalisation bookkeeping is off — fix it), and then characterise Γ: ζ-dependence (polynomial? degree?), t = 1/θ₁₂
structure (Γ(t=0) = Γ¹⊗Γ²?), sparsity/selection rules (weights: which (ij),(i′j′) are nonzero), rank, whether Γ factorises as a product of
one-site objects times a coupling (as Θ does), and whether the 5×5 block over U2's five monomials suffices (does the norm need all 8 ν's per site?).
Also test whether Γ-norm equals a candidate from (ii) (state-independent ratio) — that would connect the SoV measure to the no-log Jacobian.

**D2 — the untwisting limit.** On U2's regime caches (six HW-limit branches, 12 ε's, φ and φ′): N(Ψ_B)(ε) ~ ε^{−2(2+k_B)}/X_B² (since ⟨x₀^G|Ψ_B⟩ ~
ε^{2+k_B}X_B and ‖Ψ₀‖_G = 1) — verify, then evaluate the passing Q-side formula (at least (iv)) with the ν-series (U1's exact series in
`untwisted_Q_laurent_exact.pkl`/`untwisted_Q_pert_sym.pkl` or U2's numerical ones) and show it reproduces X_B² (W4's exact closed forms, E5_KB in
`Paul/Mathematica/Gaudin/ref/untwisted-exact-210.txt`) in the limit — the depth of the ν-series needed per branch. Then assemble the FULLY NORMALISED
CG table of adj⊗adj from Q-data only: C_ν = (U2's identity) × (1/√N) — compare with W4's exact CG values (six branches) and report the digits.

## 2. Reporting

`RESULTS-gaudin-norm.md`: §0 ground truth and conventions, §1 the candidates with pass/fail and the state-independent constants, §2 Γ (the SoV
measure: structure, formulas, the identity's residuals), §3 the limit and the normalised CG table, surprises, not done. If NO candidate passes and Γ
has no visible structure, say so plainly (the fallback is the Yang–Yang Hessian in Bethe roots — out of scope; do not start it).
