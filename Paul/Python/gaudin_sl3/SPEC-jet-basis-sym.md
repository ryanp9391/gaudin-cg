# Spec M2: the jet basis of one-row representations in the polynomial realisation — deriving the kernel instead of fitting it (round 5, WP-A2; 2026-10-16)

Delegated to a general-purpose subagent (SONNET), Python only (no Wolfram, no wolfbook tools). Read first: `Paul/Mathematica/Gaudin/
Gaudin-sl3-round5-plan.md` §0 (P2) and §1 WP-A2 (your goal; A1/A3 belong to agent M1, running concurrently — do not fit Pochhammer ansätze),
`Paul/Mathematica/Gaudin/Gaudin-sl3-round4-plan.md` §0 (notation — binding; the "Canonical kernel" paragraph defines the jet basis), the notes
`Paul/Mathematica/Gaudin/Gaudin-sl3-CG-notes.md` §6 (the B operator and the jet basis — read 6.3–6.5 carefully), §8.1–8.3, §8.8–8.9, and
`RESULTS-gaudin-B.md` (D3: how the jet basis Y was computed for the adjoint, its rationality and ζ-independence, the derivation of G from the K
tensor) and `RESULTS-general-kernel.md` K1/K3 (one-site kernel, factorisation at t = 0). Environment: Python 3.13 (`python -X utf8`), sympy 1.14,
mpmath 1.3, numpy; no `uv`; run from `Paul/Python/gaudin_sl3/`. Four other agents run concurrently (M1, N1, Q1 in Python; W5 on Wolfram): ≤ 2
heavy Python processes of your own at a time; time one unit before loops; cache to `jet_basis_sym_*.pkl`.

**Code to reuse (read-only; import or copy, never edit):** `gaudin_B.py` (D3: `jet_matrix`, `comp_jet_matrix`, `adapted_covectors`, `stage_B2`–`B4`,
`omega_comp_sym`, `d2_G_sym`, rational reconstruction helpers `rat_recon`/`mp_to_fr`/`frac_inverse`; its pkl's `gaudin_B_jets_P.pkl`,
`gaudin_B_Ycomp_P.pkl` hold the adjoint jet matrix — see the RESULTS file for formats), `general_kernel.py` (`load_rep(lam)` → exact Θ^{[Λ]} polynomials as
`T[(i,j)] = {(n,p,q): Fraction}` for Θ_ij = Σ c tⁿaᵖa′^q; `pairs_of`, `e_of`, `ckl_closed`; point builders `build_point2`/`get_point2`/`site_monos`/`overlap2`/
`omega_inf`/`x_cov`/`cfac` — the caches `general_kernel_pt_<Λ>_<tag>_<θ>_80.pkl` hold states + dual Q's + jets for Λ = 100, 200, 300, 400 at tags P/Pp/X1
(check which exist with `ls`)), `gaudin_chain.py` (`Chain`, `Operators`, `eigenstates`, `vacuum_covector`, `Pmat`, `gc`), `gl3_reps.py` (`GL3Rep`: GT basis
and generator matrices — find out its basis order and normalisation), `twovac_fit.py` (`PivQR`, `rat_recon`). **Files you own:** `jet_basis_sym.py`,
`jet_basis_sym_*.pkl`, `RESULTS-jet-basis-sym.md` (write incrementally). No notebooks, no git, no Wolfram, no edits outside your files. If an
expectation is refuted and two attempts do not fix it, stop that item, keep what works, report what you saw.

## 0. The idea

For Λ = [λ₁,0,0], V = C[x₁,x₂,x₃]_{λ₁} (homogeneous polynomials of degree λ₁) with E_ij = x_i∂_j; the two-site space is bihomogeneous polynomials
Ψ(x; y) of bidegree (λ₁,λ₁). In the Gaudin companion frame (round-4 plan §0): Ψ_c = (P_phys^{G−1})^{⊗2}Ψ, x₀^c = e₁⊗e₁ (= the coefficient of
x₁^{λ₁}y₁^{λ₁}), Φ_c∞ = Ωc∞[ω₁]⊗Ωc∞[ω₁′] with Ωc∞[ω] = e₁ᵀρ(e^{aE₁₂}e^{aE₂₃}), a = 1/ω. **Hypothesis H0:** in the polynomial realisation
e₁ᵀρ(g) is evaluation at the point e₁·g (the first row of g; or its transpose — fix the convention from the code), so Ωc∞[ω] = evaluation at the
CONIC point (1, a, a²) (the first row of (1+aE₁₂)(1+aE₂₃) = 1 + aE₁₂ + aE₂₃ + a²E₁₃), and the jet formula reads
Ψ_c(1,a,a²; 1,a′,a′²)/Ψ_c(e₁;e₁) · M₁M₂ = Σ_ij Θ_ij(a,a′,t) m_i(θ₁)m_j(θ₂). The canonical jet basis is the set of covectors ⟨i,j| (i,j over the
C(λ₁+2,2) = dim V pairs per site, so (dim V)² covectors) with ⟨i,j|Ψ_c⟩/⟨x₀^c|Ψ_c⟩ = m_i(θ₁)m_j(θ₂)/(M₁M₂) for EVERY eigenstate Ψ_c; D3 found it
ζ-independent and rational (adjoint), and Θ = the expansion coefficients of Φ_c∞ in it. At t = 1/θ₁₂ = 0, Θ = θ(a)⊗θ(a′) with θ_kl(a) = c_kl a^{k+l−1},
c_kl = (−1)^{k+l−1}C(λ₁, l−1)C(λ₁+1, k)(l−k)/l, so the one-site jet covector ⟨kl| is, up to c_kl, the coefficient of a^{k+l−1} in the evaluation
along the conic: a derivative functional at e₁ (a "jet" — hence the name). The 1/θ₁₂ corrections couple the sites; the goal is their closed form.

## 1. Tasks

**J0 — the realisation.** Build the monomial basis x^α (|α| = λ₁) with E_ij = x_i∂_j, and the explicit change of basis to the GT basis used by
`gaudin_chain.Chain`/`gl3_reps.GL3Rep` for [1,0,0], [2,0,0], [3,0,0] (one vector per weight; fix the scalars by matching E_ij matrix elements exactly).
Verify H0: compute e₁ᵀρ(e^{aE₁₂}e^{aE₂₃}) with the code's ρ and conventions (`gaudin_chain.vacuum_covector` at ω₂ → ∞ / `general_kernel.omega_inf`,
and the companion frame `Pmat`/`gc`/`P_phys^G` of the round-1 construction) and compare with evaluation at (1,a,a²) (or at (1,a,a²) transposed /
with other signs) — state the exact dictionary, including the normalisation of x₀^c. If H0 fails, find what Ωc∞ is in the realisation (it must be
evaluation at SOME point of the conic or a derivative thereof) and continue with that.

**J1 — the jet basis exactly for [1,0,0] and [2,0,0].** With the cached base points (tags P, Pp, X1 — three twists; two θ's if cached, else build one
more with `build_point2` at 80 digits — cheap for these dims): Y = W·Ψ⁻¹ where Ψ = the matrix of companion-frame eigenstates (columns) and
W_{(ij),α} = m_i^{(α)}(θ₁)m_j^{(α)}(θ₂)/(M₁^{(α)}M₂^{(α)})·⟨x₀^c|Ψ_c^{(α)}⟩. Rational-reconstruct Y (D3's method), confirm ζ-independence, and write every
covector ⟨i,j| as a linear functional on bihomogeneous polynomials, i.e. as a row over the monomial basis x^αy^β. Then find its form: (a) at
t = 0 the product of one-site jet functionals (verify: ⟨kl| ∝ ∂_a^{k+l−1}|_{a=0} of evaluation at (1,a,a²), i.e. the coefficient extraction above —
check the constants); (b) at t ≠ 0: expand Y = Y₀ + tY₁ + t²Y₂ + … (Y is polynomial or rational in t? D3: 59 of 64 adjoint covectors couple the sites
through 1/θ₁₂ — determine the t-degree) and look for the correction operator: test the hypotheses Y(t) = Y₀·(1 + tA + t²A₂ + …) with A built from the
site-exchange Casimir Ω₁₂ = Σ_ik E⁽¹⁾_ikE⁽²⁾_ki (the Gaudin Hamiltonian's cross term) and/or the diagonal gl(3) generators E⁽¹⁾+E⁽²⁾ (the jet
basis is defined by ζ-independent eigen-data, so A should commute with the total gl(3) or at least with the total Cartan — test which); or
Y(t) = Y₀·exp(tA). Express A as an explicit operator on bihomogeneous polynomials (differential operator in x, y). (c) The monomial side: are the
m_i(θ) exactly the jets of the SoV "wavefunction" in the sense that the left eigencovectors of some commuting family are the ⟨i,j|? (D3's B
hierarchy: b₁, b₂ — if you can cheaply check that ⟨i,j| are joint eigencovectors of D3's b₁ for [1,0,0], do it; otherwise skip.)

**J2 — Θ from Y, and the closed form.** Θ = coefficients of the functional "evaluate at ((1,a,a²),(1,a′,a′²))" in the basis ⟨i,j| — i.e. solve
ev_{a,a′} = Σ Θ_ij ⟨i,j| (a square linear system over Q(a,a′,t); with Y in closed form, Θ follows by a symbolic inverse). Verify against
`general_kernel_rep_100.pkl` and `_200.pkl` exactly (with the normalisation M = P³p_{0,1}, Θ = K/(ω₁ω₁′)^{2λ₁}). If J1 gave A in closed form, derive
Θ^{[λ₁,0,0]} for general λ₁ symbolically (sympy, λ₁ as a parameter, or at least λ₁ = 3, 4 explicitly from the operator formula WITHOUT any
state data) and compare with `general_kernel_rep_300.pkl`/`_400.pkl` — that is the derivation of the master kernel. For λ₁ = 1 the answer must be the
R-matrix form Θ = (1+z)θ⊗θ′ − zθ′⊗θ, z = taa′/(a′−a): derive it.

**J3 — [3,0,0] numerically if J1 did not close.** Y for [3,0,0] (100×100) from the cached 80-digit points (rational reconstruction of 10⁴ entries —
feasible; use denominators bounded by 10¹²; if the reconstruction fails, work with 80-digit numerics), and test the J1 hypotheses there; report the
t-degree per covector and the structure you see.

## 2. Reporting

`RESULTS-jet-basis-sym.md`: §0 the dictionary (realisation ↔ GT basis, H0 verdict with the exact point), §1 Y for [1,0,0] and [2,0,0] (closed forms;
the t-structure; A if found), §2 Θ derived (λ₁ = 1 symbolic derivation; the general formula or the operator rule), §3 [3,0,0] (if done), surprises,
not done. Exact where possible. Keep the RESULTS self-contained: another agent (M1) is fitting the same kernel and the main session reconciles.
