# Hand-off

Where the last session left off. Read this first when picking the project back up; it's not a
history log (that's `git log`) — only what's still open or in-flight. Overwrite
freely each session rather than appending.

## Moved out (2026-10-11): su(4) alternating chain f⊗f̄ thread

The su(4) alternating-chain thread (SoV basis, localisation, free-parameter rotations on few sites) now lives in the
**abjm-c123** repo: `C:\Users\pryan\git-projects\abjm-c123\Paul\Mathematica\su4-alternating\README.md`. Nothing of it
remains here.

## >>> NEWEST (2026-10-16 evening): sl(3) GAUDIN ROUND 5 — ALL EIGHT SONNET agents returned; the per-agent result blocks and the STATE after round 5 are below the "ROUND 5 DISPATCHED" line. Earlier the same day: ROUND 4 DONE — target = gl(3) CG from Gaudin SoV; (A) general-Λ kernel: one-site closed form + block theorem + single-index closed form + first order (G1); (B) the analytic untwisting limit established for the adjoint — the overlap depends only on the Laurent data (U2), whose ε-series are exact/algebraic (U1, W4); all four SONNET agents returned ~17:40 <<<
**Paul (2026-10-16):** "gl(4) is no longer the target, project target is gl(3) CG coefficients from SoV in gaudin model. … can you guess a general formula from the examples … take the untwisting limit analytically to get <v x hws|Psi_0> in terms of Qs where v is some rotated state … SONNET subagents" (a previous session's fable agents burned his token budget — memory `feedback-sonnet-subagents-only`; the target change is memory `project-target-gl3-gaudin-cg`; CLAUDE.md still says su(4) rectangular — outdated, not edited).
**Plan (read first):** `Paul/Mathematica/Gaudin/Gaudin-sl3-round4-plan.md`. §0 restates the established jet formula for general Λ (non-forced Plücker pairs, M = P³p_{0,λ₂+1}, c = (ζ₂−ω₁)^{λ₁−λ₂}(ζ₃−ω₁)^{λ₁}(…)(ζ₁−ζ₂)^{2λ₂}, canonical Θ = K/(ω₁ω₁′)^{2λ₁−λ₂} polynomial in a = 1/ω, rank 1 at 1/θ₁₂ = 0, the jet basis) and the chain-side limit (⟨Φ∞|Ψ_B⟩ = ⟨φ_u⊗e₁|Ψ₀⟩(1+O(ε)) = the CG generating function in (ω̂₁, ω̂₁′) with the ε-free Gauss parameters a′, c′, g′). §1 = Goal A route: one-site kernel exact for general Λ (the θ₁₂ → ∞ factor), Θ for [1,1,0], [2,2,0], [3,1,0] laws, [4,1,0], λ₂ = 0 probes, the top-state sum rule, 1/θ₁₂ expansion, conjecture. §2 = Goal B: the difficulty (at ζ = 0 the dual Q's degenerate — Wronskians ∝ (ζ_b−ζ_a) — so the limit is a statement about the ε-jet of the Q's; the total Σ K N N is ε^{N_B} smaller than its terms, N_B = 8,7,7,6,6,4), split into U1 (untwisted Q-system: Casimir closed forms of the ζ = 0 Laurent data, the untwisted oper's exponents incl. ∞, polynomial solutions, ζ-jet basis, ₃F₂ identification; ε-series of the branch Q's and jet monomials, exported for U2), U2 (the leading polynomial L_B(ω̂) = Σ_{d+n+m=N_B} G^{[d]}N^{(n)}N^{(m)} and the identity with ĉ·⟨e₁ρ(n′)⊗e₁|Ψ₀⟩, the CG dictionary, φ-independence), W4 (exact Wolfram reference: HW decompositions + CG tables + untwisted Laurent data + Casimirs + oper exponents + hypergeometric identification + degenerate perturbation theory in ε with symbolic φ, written to `ref/untwisted-exact-<Λ>.txt` for the Python agents).
**Agents (all SONNET; specs self-contained):** W4 wolfbook-builder on K1 → `Experiments/Gaudin-untwisted-exact-spec.md` / `Gaudin-untwisted-exact.wb` + `ref/untwisted-exact-*.txt`; G1 → `Paul/Python/gaudin_sl3/SPEC-general-kernel.md` / `general_kernel.py`, `RESULTS-general-kernel.md`; U1 → `SPEC-untwisted-Q.md` / `untwisted_Q.py`, `RESULTS-untwisted-Q.md`, `untwisted_Q_series_210.pkl`; U2 → `SPEC-untwist-limit.md` / `untwist_limit.py`, `RESULTS-untwist-limit.md`. Earlier Python files read-only; ≤ 2 heavy processes per agent (round 3: CPU contention stalled a run; a stale round-2 sympy process was killed at ~14:50).
**U2 RETURNED (2026-10-16 ~16:35; `Paul/Python/gaudin_sl3/RESULTS-untwist-limit.md`, `untwist_limit.py` stages R,N,NU; 200 digits, 16-point ω̂ grid, 12 ε-values, both φ and φ′) — THE ANALYTIC UNTWISTING LIMIT IS ESTABLISHED FOR ALL SIX ADJOINT BRANCHES, and the Q-number is a function of the LAURENT DATA ONLY:** (L1) ⟨Φ∞|Ψ_B⟩ ~ ε⁰, M₁M₂ ~ ε², ⟨x₀^G|Ψ_B⟩ ~ ε^{2+k_B} with k_B = 0,1,1,2,2,4 for (4,2,0),(4,1,1),(3,3,0),(3,2,1)×2,(2,2,2), because ⟨x(φ)⊗x(φ)|Ψ₀⟩ = 0 on every non-top HW ket (first contribution from Ψ_{k_B}: X_B := ⟨x(φ)⊗x(φ)|Ψ_{k_B}⟩); total N_B = 8 − k_B = 8,7,7,6,6,4 ✔ D1. (L2) all orders t < N_B of Σ G(εω̂)N N cancel (residuals 1e-51..1e-84); depths: (4,2,0) d = 2..6, n,m ≤ 4; (4,1,1),(3,3,0) d = 1..5, n,m ≤ 5 (odd n only); octets d = 0..4, n,m ≤ 4; singlet six triples only. **Identity (polynomial in ω̂₁, ω̂₁′, coefficient residuals 1e-39..1e-75): L_B(ω̂) = ĉ·⟨φ_u⊗e₁|Ψ₀⟩·m_B/X_B with κ_B = 1 on every branch**, m_B = ε²-coefficient of M₁M₂ (top: m_B = −(φ₃−φ₂)²θ₁₂²); each octet copy pairs with its own HW ket; the CG dictionary recovers the exact slices (−1/√2, −1/√2, octets (0.2088872897, 0.0493116000)/(−0.0797878449, 0.3379867346), singlet 0.1178511302; full 64-component Ψ₀ to 3e-77); **φ-independent** (φ vs φ′ agree to 2e-77 while N^{(n)}, ĉ, X_B, m_B all change). **(L3) THE REDUCTION:** with exponents {0,2,4}, N_i = M·ν_i, ν = (1, κ, π, κπ, ε_s), and EXACTLY (sympy): κ = a_s/3, π = 2H_s − a_s − 2σ₁ − 6/δ, ε_s = −(2/3)[H_t a_s/δ − a_s b_s + a_s e₂ + 3a_s/δ² − 3a_t/δ² − 3b_t/δ − 3e₃], b_s = a_s(H_s − a_s)/3 (no-log condition, 1e-199), δ = θ_s − θ_t, σ₁ = Σζ, e₂, e₃ of ζ ⇒ **c·Qnum = Σ_ij G^s_ij(ω₁,ω₁′,θ₁₂) ν_i(θ₁)ν_j(θ₂)** (G^s = D2's G/4): the two-vacua overlap depends only on (H₁, a₁, H₂, a₂), ζ, θ₁₂, ω — the Q polynomials, their normalisation, M and κ_B never enter (ν form 14+ digits more accurate). **Limit in ν form: Σ G^s_ij(εω̂)ν_iν_j = ε^{6−k_B}·ĉ·⟨φ_u⊗e₁|Ψ₀⟩/X_B + O(ε^{7−k_B})** (1e-55..1e-157), needing the ν-series only to n = 3 (top), 5 ((4,1,1),(3,3,0)), 4 (octets), 2 (singlet). Exact anchors: H₁ + H₂ = Σ_{i≠j}ζ_iμ_j; H₁^{(0)} = [9 − (C₂(μ) − 18)/2]/θ₁₂ = 4, 6, 6, 9, 12 over θ₁₂; octet copies split at ε⁰ by a_s^{(0)} = (4.680852619, −12.25463125) / swapped. **What the CG formula needs: the kernel G split by ω-degree, the ε-series of (H_s, a_s) per site to the depths above (= perturbation theory of the Gaudin Hamiltonians on the HW sector, or the ε-expansion of the no-log 'Bethe equations'), the Gauss parameters (a′, g′, c′) and ĉ, one overall constant.** Surprise: ĉ does not clear the weight-(0,2,1) component's (φ₂−ω̂₁′)⁻² — harmless because that slice vanishes on every HW ket. Not done: L4 ([3,1,0]/[3,2,0]); closed forms of the (H_s, a_s) series beyond zeroth order and of X_B; U1's/W4's files did not exist yet (U2 used its own numerical series and `twovac_fit.hw_space`). → U1 and W4 redirected (~16:40): the closed-form target is now the ε-series of the Laurent data (W4: exact degenerate perturbation theory to O(ε⁵) at numeric φ; U1: perturbative solution of the no-log conditions).
**U1 RETURNED (2026-10-16 ~17:00; `Paul/Python/gaudin_sl3/RESULTS-untwisted-Q.md`, `untwisted_Q.py`; units θ₁₂ = 1: H_s, a_s at order εⁿ scale as θ₁₂^{n−1}, b_s as θ₁₂^{n−2}):** **exact ε-series of the Laurent data for all six adjoint branches** (polynomial in φ, homogeneous of degree n at order εⁿ; rational except √5 at the octets, whose copies are Galois conjugates): **top (4,2,0) exact to all orders: H₁ = 4 + ε(φ₁ + 2φ₂ + 3φ₃), a₁ = a₂ = 3εφ₃**; **(4,1,1): H₁ = 6 + ε(2φ₁+5φ₂+5φ₃)/2 + (√(1+x²/4) − 1), a₁ = 6 + (3/2)ε(φ₂+φ₃) + 3(√(1+x²/4) − 1), x = ε(φ₂−φ₃)**; **(3,3,0): H₁ = 6 + (3/2)ε(φ₁+φ₂+2φ₃) + (√(1+y²/4) − 1), y = ε(φ₁−φ₂), a₁ = a₂ = 3εφ₃** (sqrt forms read off the series to n = 6, not proved); singlet symbolic to n = 6 (H₁^{(2)} = Q/6, a₁^{(3)} = D/72 with Q = Σφ_i² − Σφ_iφ_j, D = (φ₁−2φ₂+φ₃)(φ₁+φ₂−2φ₃)(2φ₁−φ₂−φ₃)); octets exact to n = 5 over Q(√5) (29-point φ fit + PSLQ). Rules used/verified: b_s = a_s(H_s − a_s)/3; H₁ + H₂ = εΣ_{i≠j}φ_iμ_j; b₁ + b₂ = ε²(φ₁φ₂μ₃ + φ₁φ₃μ₂ + φ₂φ₃μ₁); the six no-log conditions hold to 1e-148. The order-by-order solver (no-log conditions + the two rules) works for (4,1,1), (3,3,0), (2,2,2); for the top and the octets the Jacobian has rank 5 (degenerate case not resolved — closed form / fit used instead). U2's ν-formulas independently verified to 1e-185. Exact series in `untwisted_Q_laurent_exact.pkl` (n ≤ 5), `untwisted_Q_pert_sym.pkl` (n ≤ 6), `untwisted_Q_octet_fit.pkl`. **T0 (exact, five reps, all sectors):** H₁ + H₂ = 0, b₁ + b₂ = 0, a₁ + a₂ + θ₁b₁ + θ₂b₂ = 0 as OPERATOR identities on V⊗V at ζ = 0 ⇒ **the untwisted oper is Fuchsian at ∞**; H₁ = (|Λ|² − (C₂(μ) − 2C₂(Λ))/2)/θ₁₂ ✔; **(a₁ − a₂)/2 = (μ₁+2)(μ₂+1)μ₃/(2θ₁₂)** (23 sectors; = C₃/6 − (|Λ|+1)C₂/2 + (2/3)|Λ|(|Λ|+1)² — half the product of the Weyl-shifted weight); a₁ + a₂ = −b₁θ₁₂, b₂ = −b₁, b₁ = 0 except in multiplicity-2 sectors where b₁θ₁₂² = ±√β (β = 45 for [3,2,1]; 105 and 192 for [3,1,0]'s (5,2,1),(4,3,1) and for [3,2,0]'s (5,4,1),(5,3,2); no closed form for β); octet a₁θ₁₂ = (15 ∓ 3√5)/2. W4's `ref/untwisted-exact-210.txt` matches T0 exactly. **T1 (untwisted oper):** exponents {0,2,4} at θ_s; at ∞ the solutions grow as u^k with **k = (μ₃, μ₂+1, μ₁+2) — ALL THREE untwisted solutions are POLYNOMIALS of degrees (μ₃, μ₂+1, μ₁+2)** (the Weyl-shifted weight; Fuchs relation 12 − (|μ|+3) = 3 ✔; Wronskian ∝ Π³; trivial monodromy) — refuting the plan's expectation of fewer than 3 polynomials; top-state basis {1, F, F²}, F = ∫Π; exact echelon bases in x = (u−θ₂)/θ₁₂ in RESULTS/`v0_bases()`; NOT ₃F₂ in general (only (4,1,1), (3,3,0), (2,2,2) match via X = (u−θ₂)/(u−θ₁) with degenerate integer parameters). **The monic branch p_j have NO coefficientwise ε → 0 limit** (coefficient c_m has a pole of order min(μ_j−μ₃, μ_j−m)); the finite objects are the Plücker monomials and the Laurent data; in the untwisted basis the normalised rows share a common first-order direction and differ at ε²..ε⁴ (the ζ-jet structure; jet basis not in closed form). **T2 export:** `untwisted_Q_series_210.pkl` (12 pts, 200 digits, n ≤ 10), `_ext.pkl` (24 pts, 400 digits, n ≤ 20, reliable to 18), `_phi2.pkl` (φ′); first nonzero orders of N_i as in U2's table; p₀₁ = 0 to 1e-195. Not done: degenerate perturbation theory for top/octets (closed form/fit instead; octets n ≤ 5, singlet n ≤ 6), closed forms for the ζ-jet basis and for β, T1/T2 beyond the adjoint.
**G1 RETURNED (2026-10-16 ~17:20; `Paul/Python/gaudin_sl3/RESULTS-general-kernel.md`, `general_kernel.py` stages K1|REP|ZC|BLK|SGL|CORNER|ANALYZE; 80 digits, exact where stated; notation n = λ₁+1, pairs (k,l) with 0 ≤ k < l ≤ n, l ≥ λ₂+1, a = 1/ω₁, a′ = 1/ω₁′, t = 1/θ₁₂, e_kl = k+l−λ₂−1; Θ = K/(ω₁ω₁′)^{2λ₁−λ₂} with M = P³p_{0,λ₂+1} exactly, G = 4K):** **K0:** Θ polynomial in (a, a′, t) with rational coefficients for all 10 two-site reps computed; homogeneity: every monomial a^pa′^qt^m of entry (i,j) has p+q−m = e_i+e_j, 0 ≤ p,q ≤ 2λ₁−λ₂; exchange rule Θ_ij(a,a′,t) = Θ_ji(a′,a,−t) exact everywhere; [3,1,0] laws new (81 entries); [3,2,0] law reproduces P4. **K1 — ONE-SITE KERNEL IN CLOSED FORM (17 reps exact; guessed from Λ ≤ [4,2,0], predicted the other 8):** θ_kl(a) = c_kl a^{e_kl}, **c_kl = (−1)^{k+l−λ₂−1} C(λ₁−λ₂, l−λ₂−1) C(λ₁+1, k)(l−k)/((λ₂+1)C(l, λ₂+1))**, ω- and ζ-independent; equivalently Σ_{k<l}θ^{[λ₁,0,0]}_kl(x₁^kx₂^l − x₁^lx₂^k) = (x₂−x₁)[(1−ax₁)(1−ax₂)]^{λ₁}, i.e. Σθ_klm_kl = ε a^{−λ₂}P³·W[(1−a∂)^{λ₁}P², (1−a∂)^{λ₁}P³](θ) — a DECOMPOSABLE bivector (dressed Wronskian), ε = (−1)^{λ₂}λ₂!(λ₁−λ₂)!/λ₁!. **K2:** exact unique rational kernels for [1,0,0],[1,1,0],[2,0,0],[2,2,0],[3,0,0],[3,3,0],[3,1,0],[3,2,0],[4,0,0],[4,4,0] (+ base-point/block fits [4,1,0],[4,2,0],[4,3,0]); **λ₂ = 0 reps WORK** (exponents {0,1,λ₁+2}, S₂ double pole scalar 0 — P3's 'not scalar' not observed) and are the MASTER kernels. **K3:** Θ(t=0) = θ(a)⊗θ(a′) exact (all reps); rank of Θ_n (t-expansion) 1,2,4,8 for n ≤ 3 then capped; t-degree D_t = 4(λ₁−1) − max(0, 2λ₂−3) (λ₁ ≥ 2); **exact first-order formula for all pairs/all reps:** Θ₁ = tθ_pθ_q[αa + βa′ + γa²/a′ + εa′²/a] with ρ(p) = (n−k)(n−l), σ = k+l−1, f₂ = 2(n−1)−σ, N₀ = n(n−1): α_pq = −(3/n²)ρ(p)(ρ(q)−N₀) − (2/n)ρ(p)σ(q) + (1/n)f₂(p)(ρ(q)−N₀) + f₂(p)σ(q), β_pq = −α_qp, γ = k′l′ρ(p)/n², ε = −klρ(q)/n² (0 mismatches / 1737); extreme m-fold terms (1/m!)(n−k)_m(n−l)_m(k′)_m(l′)_m/((n)_m)² (falling factorials); λ₁ = 1: exact R-matrix form Θ = (1+z)θ(a)⊗θ(a′) − zθ(a′)⊗θ(a), z = t/(ω−ω′). **TOP-STATE SUM RULE CORRECTED:** Σ Θ_ij m_i^{top}m_j^{top} = M M′ R(a)R(a′), R(a) = (1−ζ₂a)^{λ₁−λ₂}(1−ζ₃a)^{λ₁} (the companion top state is not e₁⊗e₁; x_hw = (ζ₂−ζ₁)^{λ₂}); p₂ = Σ_k(ζ₃−ζ₂)^{−k}(Π^{λ₂})^{(k)} confirmed. θ decomposable; Θ₀, Θ₁ Plücker, Θ_{n≥2} and Θ not (Θ is not a compound matrix); det Θ^{[l,l,0]} = c(aa′)^s t^{l(l+1)/2}ν^{T−d}p_l(ν,t), ν = ω−ω′, p₁ = ν+t, p₂ = ν+3t, p₃ = ν²+6tν+3t², p₄ = ν²+10tν+15t², p₅ = ν³+15tν²+45t²ν+15t³; Θ rank-deficient at generic points for [1,0,0],[3,0,0],[3,1,0],[4,0,0],[4,1,0] (not understood). **K4 — THE CONJECTURE (status):** C0 one-site closed form (exact, 17 reps); **C1 BLOCK THEOREM: Θ^{[λ₁,λ₂,0]}_ij = [λ₂!(λ₁−λ₂)!/λ₁!]²(aa′)^{−λ₂}Θ^{[λ₁,0,0]}_ij on the non-forced pairs, i.e. Θ = θ_i(a)θ_j(a′)R^{λ₁}_ij(at, a′t) with R universal (λ₂-independent)** — symbolic for λ₁ ≤ 3 and [4,4,0] vs [4,0,0], direct fits for [4,1,0],[4,2,0],[4,3,0] (dev ≤ 3e-58); corner block (pairs (k, λ₁+1)) = N(aa′)^{λ₁−λ₂} × the single-index form; **C3 CLOSED FORM for [l,l,0] to all orders:** Θ_{kk′} = (−1)^{k+k′}C(l,k)C(l,k′)a^ka′^{k′}Σ_{m≤l−k, m′≤l−k′}C(l−k,m)C(l−k′,m′)F_{mm′}(k,k′)(at)^m(−a′t)^{m′}, F₀₀ = 1, F_{m0} = (k′)_m, F_{0m′} = (k)_{m′}, F_{mm′} = (k+m+1)_{m′−1}(k′+m′+1)_{m−1}(kk′+m′k+mk′) (Γ form: F = Γ(k+m+m′)Γ(k′+m+m′)[(k+m)(k′+m′)−mm′]/(Γ(k+m+1)Γ(k′+m′+1))) — exact l ≤ 4, **l = 5 PREDICTED AND CONFIRMED by direct [5,5,0] fits (441 states)**; every t^n term of Θ̂^{[l,l,0]} divisible by (ξ−η); C4 first order, C5 sum rule. Predictions: [4,2,0] (12 pairs, (1/36)(aa′)^{−2}Θ^{[4,0,0]}|, directly verified), [5,1,0] (20 pairs, (1/25)(aa′)^{−1}Θ^{[5,0,0]}|, a-degree 9, t-degree 16). Surprise: the jet index space carries an auxiliary sl(2) of spin (λ₁+1)/2 (θ^{(a+s)} = Möbius-translate of θ^{(a)}, exact λ₁ ≤ 4) but the two-site kernel is not covariant under a, a′ → a+s, a′+s. **NOT DONE: the closed form of the general four-index master kernel Θ^{[λ₁,0,0]} beyond first order** (single-index all orders, first order all pairs, extreme hops only); full grids for [4,1,0],[4,2,0],[4,3,0]; λ₁ = 5 except [5,5,0]. Helper `general_kernel_bld.py`.
**W4 RETURNED (2026-10-16 ~17:40; `Experiments/Gaudin-untwisted-exact.wb`, 39 cells / 14 code, NO native outputs, one clean run; data files `ref/untwisted-exact-{110,210,220,310,320,410}.txt`, TAB-separated InputForm with a documented header — vectors as rational pairs vA + Sqrt[D0]·vB, exact G-norm², Laurent data as (eta, alpha, beta) with H = eta/θ₁₂, a = alpha/θ₁₂, b = beta/θ₁₂²; the adjoint file also has E4_SERIES (orders 0..7 at φ and φ′), E4_SYMBOLIC (orders 0..3), E4_PSI, E5_KB, E5_SLICE, E5_GF):** **E1:** 45 HW branches over six reps (sectors listed in the report; multiplicity-2 copies = eigenvectors of the θ-independent A₁ = θ₁₂a₁ on the HW space, G-orthogonal); exact adjoint CG: C₁ = 1; (4,1,1) C₂ = −1/√2; (3,3,0) C₃ = −1/√2; (2,2,2) C₈ = 2/√288; octets C₄ = (21−9√5)/√(7560−3240√5), C₅ = (−3+3√5)/√(7560−3240√5) and the conjugate. **E2 (exact on all 37 sectors):** C₂(l) = Σl_i(l_i+4−2i), C₃(l) = l₁³+l₂³+l₃³+4l₁²+l₂²−2l₃²−l₁l₂−l₁l₃−l₂l₃+4l₁−2l₂−2l₃; H₁ = −H₂ = eta/θ₁₂, eta = |Λ|² − ½[C₂(μ) − 2C₂(Λ)] ✔; **a₁ = −a₂ = alpha/θ₁₂, alpha = [4T(1+T)² − 3(1+T)C₂(μ) + C₃(μ)]/6** (T = |Λ|; = U1's (μ₁+2)(μ₂+1)μ₃/2 ✔; alpha = 0 iff μ₃ = 0); b₁ = −b₂ = ±√disc/θ₁₂² on multiplicity-2 sectors only, with tr A₁ = eta + 2c₂ and disc = eta² − 4c₂² on 7 of 8 such sectors — **exception [3,2,0] (5,3,2): tr = 56 = eta + 2c₂ + 18, disc = 192** (no closed form); sum rules b₁ + b₂ = 0 and a₁ + a₂ + θ₁b₁ + θ₂b₂ = 0 as matrix identities ⇒ ∞ regular ✔ U1. **E3 (all 45 branches):** exponents {0, λ₂+1, λ₁+2} at θ_s and **{μ₁+2, μ₂+1, μ₃} at ∞** ✔ U1; three polynomial solutions, θ_s apparent singularities; untwisted DIRECT oper degrees {2λ₁+2−μ₃, 2λ₁+1−μ₂, 2λ₁−μ₁} (43/45; exceptions [1,1,0] (2,2,0), [2,2,0] (4,4,0)); **₃F₂ iff some singular point has a consecutive pair of exponents** (20/45 branches: λ₁ = λ₂ reps, or μ₁ = μ₂ or μ₂ = μ₃; explicit parameters b = (−λ₂, −λ₁−1) etc.; never for μ₁ > μ₂ > μ₃ with λ₁ ≠ λ₂, incl. every multiplicity-2 sector) ✔ U1's adjoint list. **E4 (exact Rayleigh–Schrödinger on the weight-μ space, isolating operator H₁ + a₁/7, octets in Q(√5)):** series to order 7 at φ and φ′, symbolic φ to order 3, Ψ_n to order 4 — AGREES with U1's exact series (top: H₁⁽¹⁾ = φ₁+2φ₂+3φ₃, a₁⁽¹⁾ = 3φ₃, H⁽²⁾ = a⁽²⁾ = 0, b₁⁽²⁾ = φ₃(φ₁+2φ₂), Ψ_B = e₁⊗e₁ exactly; (4,1,1): H₁⁽²⁾ = (θ₁₂/8)(φ₂−φ₃)², a₁⁽²⁾ = (3θ₁₂/8)(φ₂−φ₃)²; (3,3,0): H₁⁽²⁾ = (θ₁₂/8)(φ₁−φ₂)², a₁⁽²⁾ = 0; (2,2,2): H₁⁽²⁾ = (θ₁₂/6)Q, a₁⁽²⁾ = (θ₁₂/4)Q; octet first orders explicit in Q(√5)); H₁ + H₂ = Σ_{i≠j}φ_iμ_j at order 1 only; b_s = a_s(H_s − a_s)/3 order by order; Ψ₁ of (4,1,1)/(3,3,0) lies entirely in the (4,2,0) irrep, Ψ₁ of the singlet entirely in the octet irrep, octet Ψ₁ in all four irreps. **E5:** x(ζ) is a homogeneous polynomial covector of degree λ₂ for all six reps (explicit; `XVEC`); undo identity exact on every adjoint HW ket; **k_B = 0,1,1,2,2,4 ✔ U2**, with exact coef = ⟨x(φ)⊗x(φ)|Ψ_{k_B}⟩: (φ₁−φ₂)²; −(θ₁₂/2)(φ₁−φ₂)(φ₁−φ₃)(φ₂−φ₃); −(θ₁₂/2)(φ₁−φ₂)³; octets −(5766(√5−3)/6125)(φ₁−φ₂)²(φ₁−φ₃)(φ₂−φ₃) and conjugate; singlet (θ₁₂/2)⁴[(φ₁−φ₂)(φ₁−φ₃)(φ₂−φ₃)]² (= closed forms of U2's X_B up to normalisation); the generating function ĉ⟨φ_u⊗e₁|Ψ_B⟩ = Σ_ν W_ν s[ν] with the W_ν polynomials (W₇ has the harmless (φ₂−ω̂₁′)⁻¹ pole U2 saw) and the O(ε) corrections; all written as E5_GF. Not done: E4/E5 beyond the adjoint; closed form for the (5,3,2) exception.
**ROUND-4 RECONCILIATION (zero discrepancies):** U1 ↔ W4 — untwisted Laurent data (incl. the alpha ↔ (μ₁+2)(μ₂+1)μ₃/2 identity), exponents at ∞, polynomial solutions, the ₃F₂ verdict, and the ε-series to order 5–7 all agree exactly; U2 ↔ U1 ↔ W4 — k_B, the ν-formulas (verified independently by U1 to 1e-185), X_B closed forms; G1 ↔ D3 — Θ(t = 0) = θ⊗θ′ with G1's one-site closed form. Plan-§2 expectation corrected: the untwisted oper is Fuchsian at ∞ with THREE polynomial solutions of degrees μ + (2,1,0); the degeneracy of the twisted Q's at ε → 0 sits in the monic normalisation (poles in ε), not in the solution space.
**ROUND 5 DISPATCHED (2026-10-16 ~19:00; Paul: "read it and implement with sonnet subagents"; decisions taken by the session: A and C/D run in parallel, [4,2,0] stays in the backlog, the paper (E2) waits for A/B, CLAUDE.md goal statement updated, `Paul/Python/gaudin_sl3/*.pkl` gitignored):** wave 1 = M1 `SPEC-master-kernel.md` (A1 Pochhammer closed form of R^{λ₁}, [5,0,0] state-residual test, A3 twisted sl(2)/fusion), M2 `SPEC-jet-basis-sym.md` (A2: polynomial realisation, Ωc∞ = evaluation at the conic point (1,a,a²), the jet basis Y for [1,0,0]/[2,0,0] exactly and its t-structure, Θ derived), N1 `SPEC-nu-general.md` (B1 general-Λ ν's by Frobenius with exponents {0,λ₂+1,λ₁+2}, normalised ν_{0,λ₂+1} = 1 so P³D_s cancels trivially; B2 the ν-form untwisting limit for [3,1,0]/[3,2,0], all 20 CG tables), Q1 `SPEC-gaudin-norm.md` (D1 candidates det N / no-log Jacobian / Weyl dets / THE SoV MEASURE Γ = (Y⁻¹)ᵀG_cY⁻¹ in the jet basis; D2 limit → fully normalised CG table), W5 wolfbook-builder on K1, `Experiments/Gaudin-Laurent-series-spec.md` → `Gaudin-Laurent-series.wb` (created by the session via wolfbook_newNotebook so runCell may store native outputs) (B3 exact ε-series for the 10+10 branches of [3,1,0]/[3,2,0] → `ref/laurent-series-310.txt`/`-320.txt`, k_B; C3 the cyclic subspace C_B = span{H₁ⁿΨ₀} and its characteristic polynomial = the closed form of the series; C2 the 2×2 matrices). Wave 2 (after a Python agent returns): S1 `SPEC-bethe-series.md` (C1 degenerate Fredholm solver of the no-log conditions, n = 8 adjoint, then [3,1,0]/[3,2,0]; F3 `cg_from_Q.py` end-to-end script). All specs self-contained on disk — re-dispatch from them if agents are lost.
**N1 RETURNED (2026-10-16 ~19:45; `Paul/Python/gaudin_sl3/RESULTS-nu-general.md`, `nu_general.py` (`nu_symbolic`, `nu_numeric`, `frobenius`), caches `nu_general_reg_*.pkl`/`_res_*.pkl`; 200 digits, 24 ε-points k = 0..23 (12 were too few: deep branches need the series to ε⁹)) — THE ν-REDUCTION HOLDS FOR [3,1,0] AND [3,2,0] AND THE UNTWISTING IDENTITY HAS CONSTANT EXACTLY 1 ON ALL 20 BRANCHES:** B1: with u_m (c₀ = 1, c_{λ₂+1} = 0) and v_m (c₀ = 0, c_{λ₂+1} = 1), ν_kl = k!l!(u_kv_l − u_lv_k)/(λ₂+1)!, ν_{0,λ₂+1} = 1; ν_kl is a POLYNOMIAL in (H_s, a_s, b_s, H_t, a_t, b_t, σ₁, e₂, e₃, 1/δ) (no division by Laurent data); for k ≤ λ₂ < l, ν_kl = k!u_k·ν_{0l} (so ν_{1l} = κν_{0l}, κ = a/c₂); adjoint ν ≡ U2's `nu_funcs` symbolically; [3,1,0]: κ = a/4, π = H − a/2 − σ₁ − 4/δ; [3,2,0]: ν₀₃ = 1, ν₁₃ = a/8, ν₀₄ = 3H − a − 6σ₁ − 30/δ, ν₂₃ = −(Ha − a² − 8b)/24, ν₂₄ = ν₂₃ν₀₄ (ν₃₄, and [3,1,0]'s ν₀₄, ν₂₃, ν₂₄, ν₃₄ long, in RESULTS). No-log conditions = 3 per site: A (order λ₂+1), B1, B2 (order λ₁+2); λ₂ = 1, any λ₁: (λ₁+1)b_s = a_s(H_s − a_s); λ₂ = 2: A is cubic in (H,a,b), linear in b, and involves the other site, e₂, e₃, σ₁, δ (closed form in λ₁ in RESULTS); [3,2,0]'s B2 is P4's inhomogeneous quadratic (involves (H_s,a_s,b_s,H_t) only). Partial-fraction form verified exactly: S₂ constant e₂, S₃ constant e₃, c₂ = λ₂(λ₁+1), no triple pole. Residuals |m_i/M − ν_i| ≤ 1.5e-73 (80-digit data), no-log ≤ 1e-77, all 225 states × 2 sites × all pairs, 3 twists per rep. B2: Σ G ν ν = ε^{N_B}·ĉ·⟨φ_u⊗e₁|Ψ₀⟩/X_B with C = 1 (|C−1| ≤ 2.3e-96 [3,1,0], 4e-117 [3,2,0], both φ; polynomial identity residuals ≤ 1e-95 / 2e-127 at 64 points; lower orders cancel to 1e-103..1e-170). **RULES (26 branches): k_B = height of 2Λ − μ = 2(λ₁−λ₃) − (μ₁−μ₃); N_B = 4λ₁ − 2λ₂ − k_B.** [3,1,0] k_B = 0,1,1,2,2,2,3,3,4,5 for (6,2,0),(6,1,1),(5,3,0),(5,2,1)×2,(4,4,0),(4,3,1)×2,(4,2,2),(3,3,2); [3,2,0] k_B = 0,1,2,1,2,2,3,3,4,5 for (6,4,0),(6,3,1),(6,2,2),(5,5,0),(5,4,1)×2,(5,3,2)×2,(4,4,2),(4,3,3). Copy convention = W4's (copy 1 = smaller a₁; ε⁰ data agree to 2e-197, chain vectors to 1e-200). L_B degree ≤ 2λ₁−λ₂ per variable, total ≤ N_B. x(ζ) from W4's XVEC = `x_cov` exactly. **All 20 CG tables agree with W4's exact HW slices** (0 or ≤ 3e-203 on the 12 rational branches, ≤ 4e-106 / 3e-128 on the multiplicity-2 ones); φ vs φ′ ≤ 1e-106 / 4e-128; duality: five pairs of tables coincide (observation). **Series depth needed (the order of (H_s,a_s[,b_s]) entering the leading polynomial): [3,1,0]: (6,2,0) 5, (6,1,1) 9, (5,3,0) 9, (5,2,1) 8, (4,4,0) 8, (4,3,1) 7, (4,2,2) 6, (3,3,2) 5; [3,2,0]: (6,4,0) 4, (6,3,1) 7, (6,2,2) 6, (5,5,0) 7, (5,4,1) 6, (5,3,2) 5, (4,4,2) 4, (4,3,3) 3** (λ₂ = 1 needs only (H,a); [3,2,0] needs b through ν₂₃, ν₃₄) — W5 told to go to order 9/7 at numeric φ; S1 dispatched (wave 2) with the same target. Not done: closed form of X_B(φ,θ₁₂); [4,2,0]; λ₃ > 0; closed form of A beyond λ₂ ≤ 2; the k_B/N_B rules and C = 1 are patterns, not proofs.
**Q1 RETURNED (2026-10-16 ~20:15; `Paul/Python/gaudin_sl3/RESULTS-gaudin-norm.md`, `gaudin_norm.py --stages D0,D1,G,D2`, `gaudin_norm_output.txt`, caches `gaudin_norm_D0_{P,Pp,X1}.pkl`, `_c0.pkl`, `_G1sym.pkl`; adjoint only) — THE GAUDIN NORM FROM THE Q SIDE: N(Ψ) := ⟨Ψ|Ψ⟩_G/⟨x₀^G|Ψ⟩² = (9/4)·|det J₆|·(ζ₁−ζ₂)^{μ₃−6}(ζ₁−ζ₃)^{μ₂−6}(ζ₂−ζ₃)^{μ₁−6}**, J₆ = the 6×6 Jacobian of the six no-log conditions w.r.t. (H_s, a_s, b_s)_{s=1,2}, μ = the dual-Q degrees (the sector) — state-independent to ~1e-73 on all 64 states at P, Pp, X1, the constant checked at 34 ζ points (≤ 1.7e-73) and 5 θ values (θ-INDEPENDENT); the reduced 4×4 (b eliminated) gives the same since ∂F₁/∂b = −1. The Yang–Yang-Hessian analogue in Laurent variables. **Sign caveat: det J₆ has no definite sign (32/32 at every ζ); no sign rule found (GF(2) searches over Q-values/Laurent data/ν failed) — the formula uses |det J₆|; the overall sign of a CG slice cannot come from the Q side anyway.** Failed candidates: det N_B alone or × M₁M₂; Q-value monomials (log-regression residual 2.8 vs spread 37); D3's equal-twist ratio is an overlap identity, nothing new. **(iv) THE SoV MEASURE Γ = Y⁻ᵀ(G₁⊗G₁)_cY⁻¹ (exact by construction):** all 8 one-site ν from Frobenius, ν = (1, κ, π, κπ, κ², ε, κ²π, κε); N = νᵀΓν to 5e-68..1e-75 (3 twists, 6 θ's); COMPLETENESS Σ_A ν_Aν_Aᵀ/N_A = Γ⁻¹ to 5e-74 (1/N_A is the discrete SoV measure on the spectrum); the 5×5 block over U2's five monomials does NOT suffice (all 8 needed); Γ(t=0) = Γ¹⊗Γ¹ exactly, no factorisation for t ≠ 0 (reshaped singular values 1, .41, .41, .17, …); one-site companion metric G_c1(ζ) rational, S₃-symmetric, denominator Δ⁴, degree −2−(L_i+L_j); Γ-norm ≡ J₆-norm on every state. **SURPRISE (for M2/the derivation): ALL 833 nonzero entries of the adjoint jet matrix Y are pure monomials c·t^k (k = 0..4): Y(t) = t^L·Y(1)·t^{−L} with L the one-site level grading (0,1,1,2,2,3,3,4) summed over sites — the whole θ₁₂ coupling is a GRADING; Y(0) = Y¹⊗Y¹; C(t) = Y(t)Y(0)⁻¹ has log with level shifts 1..7 (not exp of one generator).** D2: the J₆ formula matches the direct norm to 1e-165 in the regime; lim ε^{4+2k_B}N = 1/X_B² (= W4's E5_KB to ~60 digits, sympy-limited); Ψ₀^Q = v_q/√(lim ε^{4+2k_B}N) has G-norm 1 and matches W4's exact Ψ₀ to 3e-148 (top), 1e-65 ((4,1,1),(3,3,0)), 3e-77 (octets), 4e-102 (singlet); φ, φ′ agree; **the fully normalised adj⊗adj CG table from Q-data only** (slices 1; −0.7071067811865 ×2; (0.2088872897342, 0.0493116000130); (−0.0797878448606, 0.3379867346078); 0.1178511301978). Depth: det J₆ needs the Laurent data to ε⁹ (top), ε⁷ ((4,1,1),(3,3,0)), ε⁵ (octets), none (singlet); the Γ-route cancels 12–14 orders (uneconomical in the limit). Not done: the sign of det J₆; closed forms of Γ¹ and C(t); symbolic factorisation of det J₆ (sympy too slow); general Λ (the exponents μ_k − 6 = μ_k − 2|Λ|? — to test on [3,1,0]/[3,2,0]); Yang–Yang in Bethe roots (out of scope).
**M2 RETURNED (2026-10-16 ~20:40; `Paul/Python/gaudin_sl3/RESULTS-jet-basis-sym.md`, `jet_basis_sym.py`, caches `jet_basis_sym_*.pkl` incl. `_Ycoef_ops5.pkl` (operator-derived Y₀..Y₅ for [5,0,0]), `_Ypts4.pkl`; helper `t_op4.py`) — THE JET BASIS IS DERIVED FROM THE GAUDIN OPERATORS, AND THE MASTER KERNELS λ₁ = 3, 4 ARE DERIVED (fits confirmed):** J0: V = C[x]_L, E_ij = x_i∂_j, GT vector v_p = s_p x^{α(p)} (s fixed by matching GL3Rep; HW = x₁^L with s = 1; L=1: s = (2,1,1) for (x₃,x₂,x₁); L=2: (12,6,6,2,2,1)); ρ(g)f(x) = f(xg); **H0 confirmed exactly: Ωc∞[ω] = evaluation at the conic point (1, a, a²), no signs, no transpose; x₀^c = coefficient of x₁^Ly₁^L; frame tie Ωc P_phys⁻¹ = Ω_D/κ with κ = Δ₁^L.** J1: Y ζ-independent, rational (0 reconstruction failures), a polynomial of degree exactly L in t = 1/θ₁₂; every entry a pure monomial in t with t-power = (jet level of the row) − (principal weight α₂+2α₃ of the column ket) (L ≤ 3) ⇒ Y(t) = t^ρY(1)t^{−κ} (Q1's grading, confirmed); Y₀ = y⊗y with the one-site jets = Taylor jets of f/x₁^L at e₁ in (u,v) = (x₂/x₁, x₃/x₁) (L=2: j₀₁ = c₁₁, j₀₂ = −c₁₂/2, j₀₃ = c₂₂ − c₁₃/2, j₁₂ = c₁₃/2, j₁₃ = −c₂₃/2, j₂₃ = c₃₃; Σ c_kl a^{k+l−1}j_kl = f(1,a,a²) incl. the constants); Y(t) = Y₀(1 + tA₁ + … + t^LA_L), T = Y₀⁻¹Y unipotent; **A₁ = (E₂₂+E₃₃)∧E₂₁ + E₂₃∧E₃₁ + E₃₃∧E₃₂ (X∧Y = X⁽¹⁾Y⁽²⁾ − Y⁽¹⁾X⁽²⁾), UNIVERSAL in L** (raises the principal weight by 1, shifts the gl₃ weight by −α₁ or −α₂; commutes with neither the total gl₃ nor the Cartan; not Ω₁₂); L=1: A₁² = 0, E(a,a′)A₁ = (aa′/(a−a′))E(a,a′)(1−P); M(t) = YY₀⁻¹ has small integer entries and is universal in L (one infinite matrix on Λ²U⊗Λ²U; truncations agree L=1,2,3 and L=4,5); operator-Schmidt ranks of M₁, M₂, M₃ = 6, 18, 40; A₂ ≠ A₁²/2, M ≠ exp(tM₁); no closed form of Y(1) or A_{n≥2}. **THE DERIVATION: for one-row reps the no-log condition is a_s = 0 (exponents {0,1,L+2}); every solution's jet lies in span(f₀,f₁) of the Frobenius series, p_kl/p₀₁ = ν_kl = k!l!(γ_kd_l − γ_ld_k) polynomial in (H_s,b_s;H_t,b_t;e_i(ζ);δ) (L=1: ν = (1, H_s, b_s)); hence ⟨i,j|^c = ⟨x₀^c| ν_i(Ĥ₁,b̂₁;Ĥ₂,b̂₂)ν_j(Ĥ₂,b̂₂;Ĥ₁,b̂₁) AS AN OPERATOR IDENTITY (companion-frame Laurent-coefficient operators with twist g_c(ζ))** — exact entry by entry for L = 1,2,3 at P, Pp, X1 and at ζ = 0; explains ζ-independence and t-polynomiality. J2: λ₁ = 1 symbolic: Θ = (1+z)θ⊗θ′ − zθ′⊗θ, z = taa′/(a′−a), derived; λ₁ = 2: Θ from Y = `general_kernel_rep_200` exactly (36/36 polynomials); **λ₁ = 3, 4: Θ from the OPERATOR-ONLY Y (no eigenstates, no Q's) = `rep_300` (100/100), `rep_400` (225/225)** at two random rational points each; [5,0,0]: operator-derived Y(t) exact, degree 5 (Θ^{[5,0,0]} = ev·Y⁻¹ available for M1 — told M1). Not done: the b₁-eigencovector check; ν_kl symbolic in L; closed form of Y(1)/the master kernel beyond first order; Θ^{[5,0,0]} as explicit polynomials.
**M1 RETURNED (2026-10-16 ~21:00; `Paul/Python/gaudin_sl3/RESULTS-master-kernel.md`, `master_kernel.py`, `master_kernel_bld.py` (G1's builder with the SVD tolerance patched to 1e-25 for 40-digit work), `master_kernel_tab.pkl` (R^{λ₁} exact for λ₁ = 1..4: 17, 193, 1231, 5062 entries), `mk_*.pkl`, 441-state [5,0,0] points at 40 digits) — NO FOUR-INDEX CLOSED FORM; ORDERS t⁰, t¹, t² CLOSED FOR ALL λ₁ AND ALL PAIRS, CONFIRMED ON [5,0,0]:** support rule N = m+m′ ≥ 0, −N ≤ m ≤ 2N, m, m′ ≤ 2λ₁, N ≤ 4(λ₁−1); exchange R_ij(x,y) = R_ji(−y,−x). **Symbol form:** Σ Θ_pq(x₁^kx₂^l − x₁^lx₂^k)(y₁^{k′}y₂^{l′} − …) = (x₂−x₁)(y₂−y₁)·Q(ξ;η), F_c(ξ) = Π(1−cξ_i), G_c(η) = Π(1−cη_j), Π = F_aG_{a′}, K = aa′/(a′−a)·(Π − F_{a′}G_a), e, h = elementary symmetric functions of ξ, η: **Q₀ = Π^λ; Q₁ = λ²Π^{λ−1}K; Q₂ = λ²(λ−1)Π^{λ−2}(A + λK²/2), A = −aa′(a−a′)(e₁−h₁)Π + K(K/2 + (a−a′)(aa′(e₁h₁−e₂−h₂) − ae₁ − a′h₁ + 2))**; Q_n divisible by Π^{λ−n} (n ≤ 3, λ ≤ 4) with Q_n/Π^{λ−n} polynomial in λ; order 2 fitted on λ = 2, 3, confirmed on λ = 4 (0 mismatches); order 3 needs λ = 5 data. Structure: the pair space is Sym^λ(C³), θ(a) = the coherent state ℓ_a^λ, ℓ_a = (1, −a, a²); Θ is a polynomial in the values of the monic quadratics P = Π(x−ξ), Q = Π(x−η) at ω = 1/a, ω′ = 1/a′ and at the THIRD point ω = ∞ (the vacuum); coefficients depend on z = t/(ω−ω′) only; **order 1 is exactly the gl(3) Casimir: Q₁ = aa′/(a′−a)·(λ² − Ω)v, v = ℓ_a^λ⊗ℓ_{a′}^λ; Θ(z) is NOT a function of Ω alone on v (λ=2 component spans 1,2,2; λ=3: 1,3,5,3 — each irrep component needs a second vector)**; top order Q_{4(λ−1)} ∝ Π(ξ_i−η_j) (constants 12, 8640 for λ = 2, 3); every tⁿ (n ≥ 1) vanishes on {ξ = η as sets}. [5,0,0] test (441 states, 40 digits, no fit): truncation residuals scale as t¹, t², t³ for orders ≤ 0, 1, 2 (2.0e-5 at t = 10⁻³ for order ≤ 2) ⇒ orders 0–2 correct, first missing order t³. A3: (a) θ(a+s) = exp(sN^{∧2})θ(a) exactly (λ ≤ 5) with Ne_k = −(n−k)e_{k+1}, Fe_k = ke_{k−1}, H = n−2k (an sl(2)); (b) the diagonal shift with t′ = t fails for all λ (works for λ = 1 with t′ = t·aa′/((a+s)(a′+s)); no scalar t′ for λ ≥ 2 — Θ depends on the third point); **what holds exactly (λ ≤ 3): translation of ω with t fixed, [(1−as)(1−a′s)]^{2λ}Θ(a/(1−as), a′/(1−a′s), t) = (e^{sF}⊗e^{sF})Θ(a,a′,t) — with scaling this reduces Θ to a one-variable matrix family in z**; (c) fusion refuted in four forms (naive product of λ=1 exchange factors with the coincident limit after; weighted subset sums; Θ = f(Ω)v from order 2; Toda/Hirota bilinear identities). Not done: orders ≥ 3, [5,1,0], 80-digit λ = 5 builds, the fusion question; M1 did not use M2's operator-derived [5,0,0] Y (message arrived late) — M3 dispatched to combine them.
**S1 RETURNED (2026-10-16 ~21:30; `Paul/Python/gaudin_sl3/RESULTS-bethe-series.md`, `bethe_series.py`, `cg_from_Q.py` (`python -X utf8 cg_from_Q.py [phi2]`, ~3 s), `bethe_series_{210,310,320}_sym.pkl` (exact series polynomial in φ)) — THE BETHE EQUATIONS ALONE DETERMINE EVERY ε-SERIES (26/26 BRANCHES, NO CHAIN INPUT), AND THE END-TO-END SCRIPT GIVES THE ABSOLUTELY NORMALISED adj⊗adj CG TABLE:** equations: 8 in the 6 unknowns (H_s, a_s, b_s) — N1's six no-log conditions (A, B1, B2 per site) + two global rules valid for all (Λ, μ): H₁ + H₂ = Σ_{i≠j}ζ_iμ_j and b₁ + b₂ = ζ₁ζ₂μ₃ + ζ₁ζ₃μ₂ + ζ₂ζ₃μ₁ (derived by matching e^{ζ_ju}u^{μ_j} at ∞; no third rule for a₁+a₂+θ₁b₁+θ₂b₂ at finite ζ — not needed); the no-log conditions alone have rank 3–6 at ε = 0, the sum rules are essential. **At ε = 0 the equations alone have exactly the branch points as solutions (6/6 adjoint, 10/10 [3,1,0]; [3,2,0] Gröbner unfinished, W4's values verified to satisfy them) — the untwisted data are Bethe outputs.** Jacobian (8×6) rank 6 on simple roots, rank 5 exactly on MULTIPLE roots (adjoint: top triple, octets double; 15 rank-5 branches in total). **Mechanism:** at a multiple root 3 dependencies among the 8 linearised conditions (B1/B2 of one site tie to A of the other); the kernel direction = (a_s, b_s) along b_s = a_s(H_s−a_s)/c₂ at fixed H_s with b₁+b₂ = 0 (where the twist acts at first order); Fredholm: x_n = particular + c_n·k, c_n fixed by the COMMON root of all three solvability conditions of a later order — octet-like sectors (lag 1): c₁ from a quadratic at order 2, later c_n linear at order n+1 (adjoint: c₁ = −b₁⁽¹⁾); top-type (lag 2): the order-(n+1) condition is empty in c_n, c₁ from a cubic at order 3 (adjoint top: c₁ = −4φ₃), later c_n at order n+2. U1 failed because it used a single left null vector (a cubic with no rational root); with all three the gcd has a unique rational root. Selection: roots in Q(√D), depth-first; ≤ 2 K-rational candidates for c₁, the loser dies at the next solvability; a UNIQUE solution survives on all 15 degenerate branches. Adjoint: exact to n = 10 (polynomial in φ to n = 8); = U1 exactly (n ≤ 5), = U1's 400-digit numerics to n = 8; **= W5's `ref/laurent-series-210.txt` to n = 10, both φ, 0 mismatches (symbolic orders 0–3: 144 series, 0 mismatches)**; closed forms: top H₁ = 4 + ε(φ₁+2φ₂+3φ₃), a₁ = 3εφ₃, b₁ = εφ₃(4 + ε(φ₁+2φ₂)) (terminates at n = 2); (4,1,1), (3,3,0) the √(1+x²/4) forms; singlet notes-eq. (71) + new H₁⁽⁷⁾ = 0, H₁⁽⁸⁾ ∝ Q(φ); octets eq. (72). **[3,1,0], [3,2,0]: all 10+10 branches solved uniquely, = W5's `laurent-series-310.txt`/`-320.txt` to n = 10 (E4_SERIES, both φ, 0 mismatches; E4_SYMBOLIC orders 0–3: 240+240 series, 0 mismatches)**; N1's depths (9, 7) covered; 0.1–2.4 s per run. **`cg_from_Q.py`:** Bethe series + N1's ν + `general_kernel_rep_210.pkl` + `gf_polys` + Q1's norm N = (9/4)|det J₆|Π|z_ij|^{μ_k−6} (J₆ from `untwisted_Q.frobenius_conditions` on the exact series; det J₆ starts exactly at ε^{8−2k_B}; |X_B| = lim^{−1/2} = W4's closed form (73) to 14 digits, octet rescaled by (31/35)²): k_B = 0,1,1,2,2,4, all lower orders vanish EXACTLY, dictionary consistent in Q(√5) (16 equations, 7 unknowns), **all six slices = W4's exact table to ~1e-91 (90 working digits), φ′ identical; only the overall sign of a slice is aligned to the chain.** Not done: [3,2,0] order-0 Gröbner; proof of unique-survivor selection for all Λ; algebraic closed forms for octet/singlet; `cg_from_Q.py` for [3,1,0]/[3,2,0] (needs Q2's general norm); an independent derivation of the (9/4)|det J₆| norm.
**Q2 RETURNED (2026-10-16 ~22:00; `Paul/Python/gaudin_sl3/RESULTS-gaudin-norm-general.md`, `gaudin_norm_gen.py`, caches `gaudin_norm_gen_tab_*.pkl`, drivers `q2_*.py` + `q2_*_output.txt`) — THE GAUDIN NORM FOR GENERAL Λ (λ₃ = 0), SIGNED, VERIFIED ON 15 REPS:** **N(Ψ) = ⟨Ψ|Ψ⟩_G/⟨x₀^G|Ψ⟩² = (−1)^{μ₂}·C(Λ)·det J₆·(ζ₂−ζ₃)^{μ₁−2|Λ|}(ζ₁−ζ₃)^{μ₂−2|Λ|}(ζ₁−ζ₂)^{μ₃−2|Λ|}, C(Λ) = [λ₂!(λ₁−λ₂)!(λ₁+1)!/((λ₂+1)(λ₁−λ₂+1))]²** (9/4 adjoint, 64 for [3,1,0] and [3,2,0]; NOT (c₂/2)²), with J₆ = the full 6×6 Jacobian (incl. cross-site derivatives) of N1's raw `frobenius(lam)['nolog']` polynomials ordered (A,B1,B2) per site vs (H₁,a₁,b₁,H₂,a₂,b₂) (rescaling a condition by f rescales C by f²; reproduces Q1's adjoint J₆ and 9/4 exactly); G-metric = W4's GMET1 exactly. Exponents pinned by regression (slope exactly 1 in μ at fixed ζ; −8.00000000/−10.00000000 across 13 ζ's); constant fitted on 15 reps [1,0,0]…[5,1,0] (incl. λ₂ = 0, λ₂ = λ₁, the duals) at ≤ 1e-60. Residuals max over 225 states: [3,1,0] 3.4e-73..7.4e-68 (13 ζ), [3,2,0] 1.6e-73..3.1e-68; θ-independent (13 θ's per rep, ≤ 2.4e-68). **SIGN SOLVED: sign(det J₆) = (−1)^{μ₂}·sgn(ζ₂−ζ₃)^{μ₁}sgn(ζ₁−ζ₃)^{μ₂}sgn(ζ₁−ζ₂)^{μ₃}** — depends only on the sector μ, never on the copy/state; 0 violations on [2,1,0],[3,1,0],[3,2,0],[2,2,0],[3,0,0] at 9 ζ's covering every ordering; (−1)^{μ₂} on 3556 states at ordered ζ = P. Hypothesis 'product of the two 3×3 site blocks' is FALSE (≈ half the states; J₆ not block-triangular). Regime (all 20 branches of [3,1,0]/[3,2,0], both φ, 24 ε's, N1's caches): formula = direct norm to ≤ 1e-100 at every ε; **lim ε^p N = 1/X_B² with p = 4λ₂ + 2k_B** (= 4 + 2k_B for λ₂ = 1); = N1's X_B to 1e-99..1e-146 on all 40 branch×φ combinations. Not done: λ₃ > 0; a proof; a closed form of det J₆; a condition normalisation making C = 1; Γ for general Λ. (Note: the shared scratchpad shadowed a helper script between agents — Q2 moved its drivers into the repo directory.)
**W5 RETURNED (2026-10-16 ~22:15; `Experiments/Gaudin-Laurent-series.wb` — NATIVE OUTPUTS STORED (runCell worked on the session-created notebook), saved, final summary cell; data files `ref/laurent-series-{210,310,320}.txt` with records E4_SERIES (both φ, order 10; H, a, b as two site lists of {p,q} pairs), E4_SYMBOLIC (order 3), E5_KB, E5_KBNUM, PSI (Ψ_n to order 4, sector coordinates, ℓ·Ψ_n = 0), GNORM, XVEC; header documents the branch labelling) — ALL OF B3, B3b, C3, C2 AND THE DUALITY DONE, EVERY CHECK TRUE:** B3: exact ε-series of (H₁,H₂,a₁,a₂,b₁,b₂) to ε¹⁰ at φ_A = (1/3,−2/7,5/11) and φ_B = (−2/5,3/8,7/11) for all 6+10+10 branches (fields Q, Q(√5), Q(√105), Q(√192); copy 1 = smaller a₁); symbolic φ to order 3 on all 26 branches (≤ 8 s each); isolating operator H₁ + a₁/7 separates every branch (min gap 1.08). Checks exact on all branches: eigen-equations of all six Laurent operators; ε⁰ = W4; H₁+H₂ = εΣ_{i≠j}φ_iμ_j (order 1 only); **b₁+b₂ = ε²(φ₁φ₂μ₃+φ₁φ₃μ₂+φ₂φ₃μ₁) holds for all three reps**; no-log relations at all orders (b = a(H−a)/c₂ with c₂ = 3, 4; [3,2,0]: W3's quadratic D₂ relation b₁ = 3H² − (7/4)Ha + a²/4 − 3W/θ₁₂ − βH − γa − δ, W = εφ·μ); θ-scaling; symbolic = numeric. B3b: **k_B = 2(λ₁−λ₃) − (μ₁−μ₃) confirmed on all 26 branches** (N1's rule); x(ζ) homogeneous of degree λ₂; exact symbolic X_B = numeric leads. **NEW EXACT OPERATOR STRUCTURE: on every HW weight sector H₁(ζ) = H₁(0) + Σ_i ζ_i(|Λ| − E⁽¹⁾_ii) EXACTLY (no ζ², ζ³ terms), H₁(0) = (|Λ|² + cas₂(Λ) − C₂^{tot}/2)/θ₁₂; a₁ linear, b₁ quadratic in ζ ⇒ χ_B is the characteristic polynomial of a LINEAR PENCIL A + εB.** C3: d_B = dim of the cyclic subspace of Ψ₀ under H₁ (= under {H₁,a₁} = under all six Laurent operators): adjoint 1,2,2,6,6,8; [3,1,0] 1,2,2,6,6,3,10,10,14,18; [3,2,0] 1,2,3,2,6,6,10,10,14,18 — **d_B = the full weight-μ space dimension of V⊗V in every case except the adjoint singlet (8 of 10: the char poly of H₁ there is (linear)² × (irreducible degree 8); the two missing dimensions are a decoupled doubly-degenerate eigenvalue linear in ε, symmetry origin unknown)**; χ_B irreducible over Q(ε) (or Q(ε,√D₀)) and annihilates the branch series to ε¹⁰; adjoint χ's at φ_A: top λ − 140/31 − (260/231)ε; (4,1,1): 51279921λ² − (578966850 + 77474859ε)λ + 1568813400 + 437358075ε + 22237540ε² (/51279921); (3,3,0): 51279921λ² − (578966850 + 147180033ε)λ + 1568813400 + 830855025ε + 100693580ε²; octets degree 6 (same polynomial for both copies = det(λ−H₁) on the weight space); singlet degree 8. **C2 (d = 2 branches, symbolic ζ, θ₁₂; adjoint (4,1,1),(3,3,0); [3,1,0] (6,1,1),(5,3,0); [3,2,0] (6,3,1),(5,5,0)): in the basis (HW ket, G-orthogonal partner = the HW-type state of the NEIGHBOURING irrep ν = μ + e_i − e_j), H₁ = ½Σ_{i≠j}ζ_iμ_j·1 + [[η_μ/θ₁₂, (ζ_j−ζ_i)/2],[(ζ_j−ζ_i)/2, η_ν/θ₁₂]], eigenvalues ½Σ_{i≠j}ζ_iμ_j + (η_μ+η_ν)/(2θ₁₂) ± √(((η_μ−η_ν)/(2θ₁₂))² + (ζ_i−ζ_j)²/4)** (reproduces the sqrt forms; trace/det = χ_B); a₁ = α + βH₁ on C_B with β = c₂ for (4,1,1)-type, 0 for (3,3,0)-type; for d > 2 the pencil entries are linear in ζ. **Duality [3,2,0] ≅ [3,1,0]*⊗det³, exact on all 10 branches to order 10:** μ′ = (6−μ₃, 6−μ₂, 6−μ₁), ζ′ = −(ζ₃,ζ₂,ζ₁), copies swap; H₁′ = H₁ + 6/θ₁₂ − 6Σζ, H₂′ = H₂ − 6/θ₁₂ − 6Σζ, a₁′ = 4H₁ − a₁ − 24/θ₁₂ − 12Σζ, a₂′ = 4H₂ − a₂ + 24/θ₁₂ − 12Σζ (6/θ₁₂ = the dual-rep Casimir shift). Not done: Lanczos/G-orthonormal pencil form for d > 2; a₁, b₁ as polynomials in the H₁-root for d > 2; fully symbolic-ζ χ_B for d > 2; symbolic φ beyond order 3 (X_B symbolic to order k_B). Tooling notes: wolfbook read/search/save need the full absolute notebook path; after an abort the tool status can be stale; long tool messages are truncated (W5 read results via scratch dump files).
**THE STATE OF THE PROJECT (gl(3) CG from Gaudin SoV) AFTER ROUND 5:** the whole chain from Q-data to normalised CG coefficients is now closed for the adjoint and verified for [3,1,0]/[3,2,0]: (1) THE OVERLAP IS A FUNCTION OF THE LAURENT DATA for general Λ: c·⟨Φ∞|Ψ⟩/⟨x₀^G|Ψ⟩ = Σ_ij G_ij ν_i(θ₁)ν_j(θ₂), ν_kl = k!l!(u_kv_l − u_lv_k)/(λ₂+1)! POLYNOMIAL in (H,a,b)₁,₂, σ₁, e₂, e₃, 1/δ (N1), and the jet basis itself is ⟨i,j| = ⟨x₀|ν_i(Ĥ,b̂)ν_j(Ĥ,b̂) as an OPERATOR identity (M2) — the Gaudin SoV basis is derived, not fitted; (2) THE KERNEL: Θ = θ_i(a)θ_j(a′)R^{λ₁}_ij with R derived from the operator jet basis for every λ₁ (one universal matrix N for all λ; exact kernels λ₁ ≤ 6; M2, M3), and in symbol form Q^{(λ)}(t) = N(t)[F_a^λG_{a′}^λ] with N(t) a universal bidifferential operator, normal symbol exp(ΣtⁿD_n), D₁ = the Casimir-type bilinear, tree part by a Hamilton–Jacobi recursion, one-loop terms from n = 4; closed through t⁵ for all λ (M1, M3); the four-index 'Pochhammer' closed form does not exist in that form; (3) THE UNTWISTING LIMIT: Σ G(εω̂)νν = ε^{N_B}ĉ⟨φ_u⊗e₁|Ψ₀⟩/X_B with constant EXACTLY 1 on all 26 branches of three reps, k_B = 2(λ₁−λ₃) − (μ₁−μ₃), N_B = 4λ₁ − 2λ₂ − k_B (N1, W5); (4) THE ε-SERIES OF THE LAURENT DATA are determined by the Bethe side alone (six no-log conditions + two sum rules; degenerate Fredholm solver with unique survivor; S1), exact to ε¹⁰ (W5), and have a closed form as roots of the characteristic polynomial of the LINEAR pencil H₁(ζ) = H₁(0) + Σζ_i(|Λ| − E⁽¹⁾_ii) on the weight-μ space (W5; d = 2 branches fully explicit with the neighbouring-irrep partner); (5) THE NORM from the Q side: N = (−1)^{μ₂}C(Λ)det J₆·Π(ζ_a−ζ_b)^{μ_c−2|Λ|} with J₆ the Jacobian of the no-log conditions (the Yang–Yang Hessian in Laurent variables), signed, 15 reps (Q1, Q2), plus the SoV measure/completeness Σν_Aν_Aᵀ/N_A = Γ⁻¹; (6) `cg_from_Q.py` reproduces the absolutely normalised adj⊗adj table to 1e-91 in 3 s from Q-data only (S1; verified by the session). **Open after round 5:** the master kernel's D_n for all n (D₁..D₅ known; the tree recursion is a Hamilton–Jacobi equation, the one-loop source is unidentified; the per-irrep z-dependence not extracted), λ₃ > 0 (everything is λ₃ = 0), proofs (all statements are exact numerics on ≤ 15 reps), the sign/phase convention of CG slices (not a Q-side quantity), a closed form of det J₆ and of X_B(φ,θ₁₂), the adjoint singlet's two decoupled dimensions, `cg_from_Q.py` for [3,1,0]/[3,2,0] (all ingredients exist), the paper (`paper/gaudinCG.tex` decision), the notes update (`notes-src/` generator — round 5 not yet in the notes), Clean promotion + docs agent, Paul's rerun of the older Gaudin notebooks (the round-5 notebook HAS native outputs), commit of the Gaudin tree (pkl's now gitignored).
**M3 RETURNED (2026-10-16 ~23:00; `Paul/Python/gaudin_sl3/RESULTS-master-kernel-2.md`, `master_kernel2.py`, drivers `mk2_*.py`, `master_kernel2_rep_500.pkl` (exact [5,0,0] kernel, 441 entries, 15875 monomials, t-degree ≤ 16, a-degree ≤ 10 — M1's support rule) and `_rep_600.pkl` (exact [6,0,0], 784×784, t-degree ≤ 20), `_P3coef.pkl`) — THE MASTER KERNEL IS A NORMAL-ORDERED EXPONENTIAL OF ONE UNIVERSAL BIDIFFERENTIAL OPERATOR; CLOSED THROUGH ORDER t⁵; WP-A ACCEPTANCE MET BY THE OPERATOR ROUTE:** K5: M2's Y(t) is graded, Y_rc(t) = t^{ρ_r−κ_c}Y_rc(1) (ρ = jet level, κ = principal weight of the ket; 0 conflicts in Y₀..Y₅) ⇒ ONE exact sparse inverse of the 441×441 Y(1) gives Θ^{[5,0,0]} as exact polynomials (32 s, no grid); reproduces rep_300 (100/100), rep_400 (225/225). **Acceptance on 441 [5,0,0] states (no fit, no truncation): max relative residual 6e-56, median 5e-70 at 80 digits** (40-digit maxima 1e-15..8e-7 are eigenstate conditioning noise). **UNIVERSALITY: M(1) = Y(1)Y(0)⁻¹ and its inverse N are universal in λ₁ (restrictions of N⁽⁵⁾, N⁽⁶⁾ to smaller label sets = the smaller-λ matrices exactly), so ALL kernels at once: Θ^{(λ)}_pq = Σ_rs c_rc_s a^{ρ_r}a′^{ρ_s}t^{ρ_r+ρ_s−ρ_p−ρ_q}N_{rs;pq}, c_kl(λ) = `ckl_closed`** (exact for λ = 2..5 on every monomial). K6: Q₃ = Π^{λ−3}P₃, P₃ = λ²(λ−1)(c₀ + c₁λ + c₂λ² + c₃λ³), c₃ = K³/6 (fitted λ = 2..5, confirmed λ = 1..6 by the operator formula). **Cumulants W_n = Πⁿ[tⁿ]log(ΣQ_mt^m/Π^λ) are polynomials in λ of degree exactly n+1 with no λ⁰, λ¹ term (n ≤ 5), Newton form W_n = λ²Σ_{j<n}(λ−1)_jΓ_{n,j}, Γ_{n,0} = (−1)^{n−1}Kⁿ/n. GENERATING FUNCTION: Q^{(λ)}(t) = N(t)[F_a(ξ)^λG_{a′}(η)^λ] with N(t) a universal bidifferential operator (independent of λ, a, a′) whose NORMAL symbol is exp(S), S = Σ_{n≥1}tⁿD_n: D₁ = Σc_ij p_{E_i}p_{H_j} (c₁₁ = E₁−H₁, c₁₂ = E₂−H₂, c₂₂ = E₂H₁−E₁H₂; = M2's A₁ in symbol space), D₂ = ⟨D₁,D₁⟩₁ (⟨f,g⟩₁ = Σ_k∂f/∂p_k·∂g/∂q_k), D₃ = ½(⟨D₂,D₁⟩+⟨D₁,D₂⟩); the top p-degree (tree) part obeys (n−1)D_n = Σ_{i+j=n}⟨D_i,D_j⟩₁ (n = 2..5) — a Hamilton–Jacobi-type equation tS_t − S = ΣS_pS_q; ONE-LOOP corrections from n = 4 (D₄ gains 3(c₁₂² − c₁₁c₂₂)p_{E₂}²p_{H₂}²; D₅ 42 extra monomials)** ⇒ P_n(λ) = Σ_{i,j≤n}(λ)_i(λ)_j[polynomial]; D₁..D₅ verified against the independent exact kernels (n ≤ 4 at λ = 1..5; n ≤ 5 at λ = 5, 6) ⇒ **orders ≤ 5 closed for all λ**. A_λ(z) diagonal sequence to λ = 4 (1,16,168,1296,7416,30240,79200,100800; top coefficient (2λ−1)!·C(2λ−2,λ−1)). K7: Q₁ = z(λ²−Ω)v and the spans 1,2,2 / 1,3,5,3 reproduced; no operator in {1, E¹_ij, E²_ij, E¹E²} preserves the Θ span except 1 and Ω; the second operator A = D₁ gives Av = z(λ²−Ω)v and the {Ω,A}-cyclic space = the Θ span for λ = 2 but is one vector too big at λ = 3 (1,3,6,3) — **the single-A span test fails; the exact structure is the normal-ordered exponential, whose higher D_n are not polynomials in A and Ω (D₂ is twice what exp(tX) would give).** Not done: a closed form of all D_n (tree recursion + one-loop source unknown; the HJ equation unsolved); the per-irrep z-dependence; [5,1,0]; orders ≥ 6.
**ROUND 5 PLAN WRITTEN, NOT DISPATCHED (2026-10-16 ~18:00; Paul: write the plan, do not implement yet):** `Paul/Mathematica/Gaudin/Gaudin-sl3-round5-plan.md` — WP-A master kernel R^{λ₁} in closed form (Pochhammer ansatz; derivation via the jet basis of Sym^{λ₁}; the auxiliary sl(2)); WP-B general-Λ ν-formulas and the untwisting limit for [3,1,0]/[3,2,0]; WP-C mechanism/closed form of the ε-series of the Laurent data (degenerate solver, the 2×2 mixings, Casimir algebra on isotypic components); WP-D the norm ⟨Ψ₀|Ψ₀⟩ from the Q side (Gaudin norm, no-log Jacobian / det N); WP-E write-up (notes now, paper later — decision on `paper/gaudinCG.tex`); WP-F housekeeping. Suggested sonnet allocation and decision points for Paul in §2–§3. **NOTES WRITTEN (sonnet, 2026-10-16 ~18:30):** `Paul/Mathematica/Gaudin/Gaudin-sl3-CG-notes.{tex,md,pdf}` — 56 pages, 80 equations, 24 tables, 46 footnotes; hep-th conversational style, punchline-first sections 0–13 (overview, setting, Q-functions, localisers, chain-side untwist, jet formula, B operator/jet basis, Laurent-data form, general kernel, untwisting limit, ε-series, untwisted model, 27 surprises, outlook = round-5 plan) + appendices A–F (conventions incl. the kernel-normalisation table and basis orders, all 25 adjoint G entries and the exact 9×9 [3,1,0] K whose 7×7 sub-block is [3,2,0], one-site kernels for 14 reps, ε-series and exact CG tables for all 45 HW branches of six reps, untwisted data, reproducibility). Build sources preserved in `Paul/Mathematica/Gaudin/notes-src/` (build.py regenerates and OVERWRITES the .tex/.md — hand-edit the .tex only if you stop using the generator; README there). Source conflicts resolved in footnotes (later round wins: c₂ = 0 scalar for one-row reps, c₂ω₃/c₂ζ₃/Π^{λ₂} for [3,2,0], direct exponents {0,2,5}, H₁ sign, G = 4K and Θ conventions). The writer's own check: the block theorem in K normalisation, K^{[λ₁,λ₂,0]} = [λ₂!(λ₁−λ₂)!/λ₁!]²K^{[λ₁,0,0]}, exact on [3,1,0]/[3,2,0]/[3,0,0] pairs. **Gaps the notes expose:** all general-Λ claims rest on λ₃ = 0 and three reps; the ν-reduction and the limit are proved for the adjoint only; the ε-series are fits/closed forms, not derived; the master kernel beyond first order; the CG normalisation still needs a chain-side G-norm; no end-to-end `cg_from_Q.py` yet; the XXX-vs-Gaudin octet-mixing comparison and the [3,1,0] (4,3,1) mixing were never done.
**THE STATE OF THE PROJECT (gl(3) CG from Gaudin SoV) after round 4:** the two-vacua overlap is a Q-number depending only on the Laurent data: c·⟨Φ∞|Ψ⟩/⟨x₀^G|Ψ⟩ = Σ G_ij ν_i(θ₁)ν_j(θ₂) with ν rational in (H_s, a_s, b_s = a_s(H_s−a_s)/c₂) (U2, adjoint; expected rep-generic via Frobenius); the kernel G = (ω₁ω₁′)^{2λ₁−λ₂}Θ with Θ = θ_i(a)θ_j(a′)R^{λ₁}_ij(at, a′t), θ in closed form for all Λ, R universal per λ₁ (block theorem), closed to first order in t for all λ₁ and to all orders for λ₂ = λ₁ (G1); the untwisting limit Σ G(εω̂)νν = ε^{6−k_B}ĉ⟨φ_u⊗e₁|Ψ₀⟩/X_B is a polynomial identity giving every CG coefficient of adj⊗adj from the ε-series of (H_s, a_s), which are exact and partly algebraic (U1, W4). **Not done / next (round 5 candidates):** (1) the general four-index master kernel R^{λ₁} beyond first order (the single-index Γ-form and the first-order formula are the data; try a Γ-function/hypergeometric ansatz in the four indices, or derive R from D3's jet basis for [λ₁,0,0]); (2) the ν-formulas for general Λ (Frobenius with exponents {0, λ₂+1, λ₁+2}: ν as rational functions of the Laurent data; then the untwisting limit for [3,1,0]/[3,2,0] — U2's L4); (3) a proof/closed form for the ε-series of the Laurent data (the degenerate cases: top and octets' rank-5 Jacobian; the algebraic sqrt forms suggest 2×2 mixing — find the general mechanism); (4) WRITE THE PAPER SECTION (`paper/gaudinCG.tex`: conventions, the jet formula with the one-site kernel and the block theorem, the jet basis, the Laurent-data form, the untwisting limit and the adjoint CG table) — ask Paul how to handle the existing gl(4) stub; (5) [4,2,0] multiplicity 3; (6) Clean promotion + docs agent; (7) Paul reruns the four Gaudin notebooks in VS Code (no native outputs); gitignore `*.pkl` (≥ 150 MB) before committing; CLAUDE.md's su(4) goal statement is outdated.

## (2026-10-16): sl(3) GAUDIN ROUND 3 DONE — no cleaning needed (verified on [3,2,0]), the Gaudin SoV = jet basis with G derived from first principles (D3), the jet formula generalised to [3,1,0]/[3,2,0] with the count C(λ₁+2,2) − C(λ₂+1,2) (P4), [3,2,0] chain checks (W3). All three agents returned ~15:30
**Paul asked (2026-10-16 afternoon):** can we avoid the per-sector/per-copy cleaning of XXX? what about the Gaudin B operator? can we set up [3,2,0] checks?
**Plan (read first):** `Paul/Mathematica/Gaudin/Gaudin-sl3-round3-plan.md`. §0 answers Q1: the cleaning belonged to the single-vacuum + localised-T approach (`new-conventions/Clean-CG-combinations.md`: O(1) contamination K/S = κ−1 from the top-weight covector meeting the raised O(ε) corrections of the HW-limit state); the two-vacua family removed it already in XXX (Stage D) and in Gaudin rounds 1–2 (K/S = O(ε) every branch, both reps; the jet formula gives CG from the Q's with no cleaning step). Untested in Gaudin: [5,3,2] in [3,2,0]² (this round), multiplicity 3 ([4,2,0]², later), weakly conditioned 2×2 extractions (W3 reports σ₂). §0 Q2: the Gaudin B is not one operator with a simple spectrum — B_XXX = ν³[1 + hb₁ + h²b₂ + …] with b₁ Cartan-like/degenerate; the SoV basis becomes a JET basis (h-expansion of the SoV covectors), the round-2 formula = (Gaudin K tensor) × (jet wavefunction monomials). §0 Q3: the full prediction table for [3,2,0] from the rep-generic rules (c₂ = 8, exponents {0,3,5} ⇒ two forced jets, gauge 3, F = Π, W ∝ Π⁵, ℓ₁ = 4, c_H = 5Σζ + 12/θ₁₂, vacuum H_s = 2ω₁+3ω₂+5ω₃ ± 12/θ₁₂, HW sectors = dual images of [3,1,0] with multiplicity 2 in (5,4,1), (5,3,2), jet formula with p_kl up to l = 4 and lone c ≤ 2, first nonzero p₀₃).
**Agents (all sonnet; specs self-contained):**
- **W3** (wolfbook-builder, K1): `Experiments/Gaudin-L2-chain-320-spec.md` → `Gaudin-L2-chain-320.wb`: A1–A8 = Hamiltonians/Laurent (c₂, triple pole, D₂ relations c_Λ), vacuum/frames/κ^G, eigenstates/HW sectors/universality matrices for (5,4,1),(5,3,2) vs duality, BOTH two-vacua families with K/S for all 10 branches + 2×2 conditioning (the cleaning verdict), localisers (ℓ₁, fourth localiser coefficients vs [3,1,0]), closed forms, A7 reference tables (225 Qnum at P), summary. No native outputs.
- **D3** (Python, exact): `SPEC-gaudin-B.md` → `gaudin_B.py`, `RESULTS-gaudin-B.md`: B0 port/verify the XXX B(u); B1 exact h-expansion (b_k in Gaudin Lax terms, commutators, joint generalised eigenspaces); B2 SoV covectors as rational functions of h → jet basis (classes, filtration, divided differences); B3 Gaudin wavefunction rule in the jet basis (verified on 64 states); B4 Gaudin K tensor for Φ_c∞ and the derivation of D2's G entry by entry (+ equal-twist vs P1 det N); B5 optional single B_G(u).
- **P4** (Python): `SPEC-Qside-320.md` → `qside_320.py`, `twovac_fit_320.py`, `RESULTS-Qside-320.md`: Q1–Q3(–Q4) for [3,2,0] (κ_s, κ₂,_s formulas, no-log conditions, det N on 225 states, Δμ rank rule) and J = the [3,2,0] jet-formula fit (stretch, ≤ 2 h; then [3,1,0] if time).
**P4 RETURNED (2026-10-16 ~14:30; `Paul/Python/gaudin_sl3/RESULTS-Qside-320.md`, `qside_320.py`, `twovac_fit_320.py` per-stage):** [3,2,0]: HW index 14; c₂ = 8 = λ₂(λ₁+1) ✔; no S₃ triple pole ✔; dual exponents {0,3,5} ✔ (**direct exponents {0,2,5}** — duality swaps the patterns: dual {0, λ₂+1, λ₁+2}, direct {0, λ₁−λ₂+1, λ₁+2}); gauge p = 3 ✔ (does not separate λ₁ from λ₁−λ₃); degrees μ_j / 6−μ_i ✔; **κ_s = a_s/c₂** (= a_s/8; adjoint a_s/3 is the same rule), **κ₂,_s = (a² − Ha + 8b)/24**; 3×5 jet rank 2, 3×6 rank 3; p₀₁ = p₀₂ = p₁₂ = 0; 3 no-log conditions per site ✔; vacuum H_s^{vac} = 2ω₁+3ω₂+5ω₃ ± 12/θ₁₂ ✔ (k = 12 = 2e₂(Λ)), b_s^{vac} = ω₃(H_s − 5ω₃) ✔, **CORRECTION a_s^{vac} = c₂·ω₃ = 8ω₃** (not |Λ|ω₃; indistinguishable on the earlier reps where c₂ = |Λ|); vacuum degrees dual (6,4,0), direct (0,2,6); Wronskians: direct F = Π ✔, **CORRECTION dual F′ = Π² = Π^{λ₂}** (= Π^{|Λ|−p}), W ∝ e^{Σζu}Π⁵ ✔; det N with Π⁻³ Lagrange columns: rank tallies obey the Δμ rule (130 rank-4 pairs), ratio identity on all 225 states to 1.6e-51 ✔, 72/84 triples; Q4: 10 branches, multiplicity 2 in (5,4,1), (5,3,2) ✔, det N ~ ε^{16,15,14} (the [3,1,0] pattern by duality). **J — THE JET FORMULA FOR [3,2,0] AND [3,1,0], AND THE GENERAL PATTERN:** c·Qnum·M₁M₂ with M = P³p₀₃ (first non-forced Plücker) and c from κ^G = Δ₁Δ₂² is EXACTLY a bilinear form in P³(θ_s)·p_kl(θ_s) over the seven non-forced Plücker coordinates {p₀₃, p₀₄, p₁₃, p₁₄, p₂₃, p₂₄, p₃₄} (residual 0, rank 49, unique, all 49 entries rational, K₇₇ = 1/9); the lone P³^{(c)} are redundant (P³′p₀ₗ = P³p₁ₗ, P³″p₀ₗ = P³p₂ₗ) — so the adjoint's N₅ = P³·{p₀₂, p₁₂, p₀₃, p₁₃, p₂₃} too. [3,1,0]: P³·p_kl with the nine non-forced pairs (l ≤ 4), unique rational kernel, K₉₉ = 1/9. **Monomials per site = C(λ₁+2, 2) − C(λ₂+1, 2)** (5, 9, 7 for [2,1,0], [3,1,0], [3,2,0]) = Plücker coordinates of the order-(λ₁+1) jet minus the forced zeros. Kernel laws [3,2,0]: polynomial in (ω₁, ω₁′) of total degree 8 − e_i − e_j, e = (0,1,1,2,2,3,4), max degree 4 per variable; depends on θ only through θ₁₂, polynomial of degree exactly 7 in 1/θ₁₂ (adjoint 4); NOT site-symmetric at fixed θ₁₂ but K_ij(ω,ω,1/θ₁₂) = K_ji(ω,ω,−1/θ₁₂) (site exchange = θ₁₂ → −θ₁₂). Not done: the J regime/CG check (W3's chain data), [3,1,0] kernel laws, the Frobenius-system count. Note: the all-stage run stalled under CPU contention from the parallel agents — run stages separately.
**W3 RETURNED (2026-10-16 ~14:40; `Experiments/Gaudin-L2-chain-320.wb`, 75 cells / 32 code, NO native outputs, one clean 8-min run; summary cell 75):** [3,2,0] prediction table: c₂ = 8 ✔; S₃ triple pole 0 ✔; indicial roots {0,3,5} ✔; ℓ₁ = ℓ₂ = 4 ✔ (D₁ 7, D₂ 26, cubic 68; 𝓛∩D₁ = {id, H₁, S3c[1,2]}; the 3-dim L₁ already reaches every weight); **X₄ = H₁² − ½H₁S3c[1,2] − (1/θ₁₂)W − S3c[1,1]** (W coefficient 1 vs 2 for [3,1,0]: λ₁−λ₂ or 2/λ₂ — needs a third rep); v₄ θ₁₂-independent with the same one-site decomposition as [3,1,0]; **the D₂ relation is NOT a scalar c_Λ**: S3c[1,1] = 3H₁² − (7/4)H₁a + ¼a² − (3/θ₁₂)W − βH₁ − γa − δ (a = S3c[1,2], β = 78/θ₁₂ + 15Σζ, γ = −45/(2θ₁₂) − (9/2)Σζ, δ = −504/θ₁₂² − 210Σζ/θ₁₂ − 18(Σζ)² − 3e₂(ζ); exact 225×225 identity at 6 points) + site-2 mirror; c_H = 5Σζ + 12/θ₁₂ ✔; f(H₁) = 5Σζ − ζ·Λ + 12/θ₁₂ ✔ (k = 2e₂(Λ)); **f(S3c[s,2]) = 8ζ₃ = c₂ζ₃** (not |Λ|ζ₃); vacuum H_s ✔, b_s ✔, **a_s^{vac} = 8ω₃ = c₂ω₃** (agrees with P4); Ωc^G Gauss parameters rep-independent ✔, κ^G = Δ₁Δ₂² with P₃ = (PauxG⁻¹)ᵀ ✔; HW sectors = duals of [3,1,0]² ✔ (multiplicity 2 in (5,4,1), (5,3,2)); universality ✔: θ₁₂·S3c[1,2]|_HW θ-independent, no scalar part — (5,4,1): [[13,−6],[−1,22]] (quadratic 6r² + 9r − 1 = 0), (5,3,2): (4/3)[[16,−1],[−2,26]] (r² + 10r − 2 = 0, literally the same as [3,1,0]'s (4,3,1)); duality: spectra = [3,1,0]'s + sector-dependent scalar, mixing R identical (4.3117/−0.8117). **CLEANING VERDICT: no cleaning step is needed anywhere, in either family** — K/S = O(ε) for all 10 branches (doubling 2.000, S exponent 0, K ≡ 0 top); (5,3,2) — the XXX per-copy-cleaning sector — has max|K/S| 1.8e-6 / 1.4e-6 (full family) and 2.4e-6 / 1.3e-6 ((∞,∞)) at ε = 1e-6; CG after Richardson ≤ 1.5e-18 (full) / ≤ 1.9e-18 ((∞,∞); (4,4,2) 1.3e-17); mixing R to 7e-19; 2×2 extraction conditioning σ₂ listed for all pairs — best pairs σ₂ ~ 0.4–2.8, the only weak case is (4,4,2) in the (∞,∞) family (σ₂ ~ 3e-4) and it is still accurate to 1.4e-17 at 200 digits. Controls: ω₁′ = ω₁ reads only top branches; ω₁′ fixed FAILS in the (∞,∞) family (K/S 0.465 at (4,4,2), 0.388 at (4,3,3) — identical to the [3,1,0] duals) but passes in the full family at finite ω₂. A7 reference tables: cells 73–74 (225 Qnum at P for both ω₁ pairs, P′, regime branches). Pitfalls: (4,4,2) has a 2-dim slice with an exactly-zero component (norm-based errors needed); the localising scan cell takes ~6 min cold.
**W3 ↔ P4 reconciliation ([3,2,0]):** zero discrepancies — c₂ = 8, exponents {0,3,5}, HW sectors, the two corrections (a_s^{vac} = c₂ω₃; f(S3c) = c₂ζ₃ matches ⟨Ω⊗Ω|S3c[1,2] = Ω₁S3^{1-site}[−2]⊗Ω₁ with the one-site double pole c₂ζ₃). The rep-generic rule set is now: c₂ = λ₂(λ₁+1); dual exponents {0, λ₂+1, λ₁+2}, direct {0, λ₁−λ₂+1, λ₁+2}; gauge p = λ₁−λ₃; κ_s = a_s/c₂; direct F = Π^{λ₁−λ₂}, dual F = Π^{λ₂}, W ∝ Π^{|Λ|}; c_H = |Λ|Σζ + 2e₂(Λ)/θ₁₂; f(H₁) = |Λ|Σζ − ζ·Λ + 2e₂(Λ)/θ₁₂; f(S3c) = c₂ζ₃; vacuum H_s = Σ(|Λ|−λ_i)ω_i ± 2e₂(Λ)/θ₁₂, a_s = c₂ω₃, b_s = ω₃(H_s − |Λ|ω₃); ℓ₁ = λ₁+1; jet monomials per site = C(λ₁+2,2) − C(λ₂+1,2) with M = P³p_{0,λ₂+1}.
**D3 RETURNED (2026-10-16 ~15:30; `Paul/Python/gaudin_sl3/RESULTS-gaudin-B.md`, `gaudin_B.py`, exact Fractions throughout) — G DERIVED FROM FIRST PRINCIPLES:** B0: XXX B(u) ported, degree 6 in u and in h, [B,B] = 0, simple label-rule spectrum, u-independent covectors. **B1 (the Gaudin B operator):** exact all-order form B/ν³ = (1 + hν′/ν + h²/((u−θ₁)(u−θ₂)))·{t₁₁(u)[t₁₁(u+h)t₂₂(u) − t₂₁(u+h)t₁₂(u)] + t₂₁(u)[t₁₁(u+h)t₂₃(u) − t₂₁(u+h)t₁₃(u)]} with T/ν = 1 + hL + h²C (C_ij = ΣE¹_ikE²_kj/((u−θ₁)(u−θ₂))); b₁ = 2L₁₁ + L₂₂ + ν′/ν (diagonal, eigenvalue Σ_s(1+a+b+c)/(u−θ_s), 25 eigenspaces of dims 1,2,4); b₂ = 2C₁₁ + C₂₂ + L₁₁′ + L₁₁² + 2L₁₁L₂₂ − L₂₁L₁₂ + L₂₁L₂₃ + (ν′/ν)(2L₁₁+L₂₂) + 1/((u−θ₁)(u−θ₂)) (needs the cross-site C terms); [b_i(u),b_j(v)] ≠ 0 individually, Σ_{i+j=n}[b_i,b_j] = 0 for n = 2..8; b₂ compressed to the b₁-eigenspaces resolves all 64 labels; b₂ alone has a Jordan block; the algebra of all u-coefficients is solvable per level block (levels of H_tot = Σ_s(2E^s₁₁+E^s₂₂): sizes 1,4,8,12,14,12,8,4,1) ⇒ simultaneously triangularisable, flag = level flag. **B5: NO single polynomial B_G(u) of degree ≤ 3 in L, L′ triangularises the jet basis**; the resolving operator is the b₁-block compression of b₂ (needs projectors). **B2 (surprise):** in the XXX companion frame the 64 SoV covectors ⟨x|(h) are rational in h of degree ≤ 5/5 and their LEADING terms are already independent (rank 64) — the degeneracy lives in the DIAGONAL frame (P_phys⁻¹(e^{hζ}) = O(h) with singular values h·(1,h,h,h²,h²,h³,h³,h⁴); all 64 diagonal-frame covectors have valuation h² with leading term ∝ x^G⊗x^G); cumulative Taylor ranks 1,5,13,25,39,51,59,63,64; jet basis by echelon saturation over Q[[h]] (divided differences); jet defects per level = coefficients of (1+2t+2t²+2t³+t⁴)² (per site 8 jets of weights 0,1,1,2,2,3,3,4). **B3 — THE GAUDIN SoV WAVEFUNCTION:** one-site monomials m₀ = P³p₀₂, m₁ = P³′p₀₂, m₂ = P³p₀₃, m₃ = P³′p₀₃, m₄ = P³′p₁₂, m₅ = P³p₂₃, m₆ = P³′p₁₃, m₇ = P³′p₂₃ (weights 2,3,3,4,4,5,5,6; 8 = dim V, one per site label); the 64 products m_i(θ₁)m_j(θ₂) are independent on the 64 states; **canonical jet basis ⟨i,j|^c with ⟨i,j|Ψ_c⟩/⟨x₀|Ψ_c⟩ = m_i(θ₁)m_j(θ₂)/(M₁M₂)** in the Gaudin companion frame — its 64×64 matrix Y is exactly ζ-INDEPENDENT (identical rationals at P, P′, X1; 833 nonzeros), only 5 rows are product covectors, the other 59 couple the sites through 1/θ₁₂; verified to 1e-73..1e-76 at three ζ. The XXX-saturation basis is block-triangular w.r.t. this one (filtration = jet-weight filtration). **B4 — G:** the Gaudin K tensor of Φ_c∞ has EXACTLY 25 nonzero entries, all on N₅×N₅, and Θ_ij = G_ij/(4ω₁³ω₁′³) as an exact polynomial identity in (1/ω₁, 1/ω₁′), 25/25, at three θ₁₂; c arises from Ωc∞[1/ω]·P_phys⁻¹(ζ) = Ω_D∞[ζ,ω]/κ^G exactly and the 4 from M_vac = 2h²P³p₀₂ per site; G has rank 1 at 1/θ₁₂ = 0 (sites decouple) and 5 otherwise. Equal twist: all 64 jet coefficients nonzero, Qnum·M₁M₂/det N_B state-independent (12.538… at P). Not done: XXX wavefunction check at h = 1/2; [3,2,0] (the 7/9-monomial kernels of P4 not yet derived this way); θ₁₂ beyond three values.
**ROUND-3 RECONCILIATION:** the K-tensor derivation (D3) = the analytic limit (D2) = the direct fit (D1), 25/25 exact; the one-site monomial set m₀..m₇ (D3) contains N₅ and explains why only 5 enter for Φ∞ (the support of the K tensor); the general count C(λ₁+2,2) − C(λ₂+1,2) with lone factor P³ only (P4) is consistent with the D3 weights (the P³′ monomials equal P³·p_{1l} by the Frobenius relation). W3 = P4 on all [3,2,0] overlaps. Round-3 answers: (Q1) no cleaning anywhere; (Q2) the Gaudin B is a solvable hierarchy, not a single operator — the SoV basis is the ζ-independent jet basis defined by the wavefunction rule; (Q3) all [3,2,0] rules hold with the c₂ corrections.
**Not done / next (round 4 candidates — Paul to choose):** (1) Paul reruns the three Gaudin notebooks in VS Code (no native outputs); (2) general-Λ derivation: the K-tensor route for [3,1,0]/[3,2,0] to DERIVE the 9/7-monomial kernels, and a closed form of the Gaudin K tensor / jet matrix Y for general Λ (Y is ζ-independent and rational — the real structural object); (3) the jet formula as a formula: K_ij polynomial laws for [3,1,0] (grid + θ) and a conjecture for general Λ from the three reps; (4) [4,2,0] (dim 27, multiplicity 3) to finish the cleaning question; (5) gl(4) rectangular Gaudin — the project target (chain template `new-conventions/L2-gl4-rect-untwist.wb`, generators `old/Gaudin/Clean/Hamiltonians_SU4.wb`; Python gl(4) generators needed); (6) the sl(3) Gaudin section of `paper/gaudinCG.tex` (nothing written yet: conventions, localisers, untwist, det N, jet formula, jet basis); (7) Clean promotion + docs agent. Nothing committed; `Paul/Mathematica/Gaudin/` and `Paul/Python/gaudin_sl3/` untracked (~100 MB of .pkl caches — gitignore `*.pkl`).

## (2026-10-16): sl(3) GAUDIN ROUND 2 DONE — layer D: the Q-side JET FORMULA for the CG-complete two-vacua family found (D2, analytic limit) and reproduced independently (D1, direct fit), CG from the Q-functions alone; [3,1,0] closed (W2, P3). All four agents returned ~12:30
**Paul's steer (2026-10-16):** "the goal of the project is Gaudin anyway, not XXX; I assumed XXX would be simpler" — so Gaudin is the deliverable, XXX the scaffold. Round 1 showed Gaudin is simpler at layers A–C (every D₁ Laurent coefficient localises by itself, vacua exactly ε-free, constants Laurent-polynomial in θ₁₂, untwist less fragile) but with its own wrinkles (direct oper = adjoint on Q/Π², state-dependent orthogonality coefficients, [3,1,0] fourth localiser needs two products).
**Plan (read first):** `Paul/Mathematica/Gaudin/Gaudin-sl3-round2-plan.md` — §1 is binding: the Gaudin (∞,∞) family Φ∞ = Ω^∞[ζ,ω₁]⊗Ω^∞[ζ,ω₁′] (ω₂ → ∞ vacua, r_k^∞ formulas there), the normalisation-free targets Qnum(B) = ⟨Φ∞|Ψ⟩/⟨x₀^G|Ψ⟩ and the x₀-free T(B) (P1's det N as anchor), and the ladder → jet dictionary: XXX's four ladder points θ_s + mh become the order-3 jet at θ_s, Bk2's become Plücker coordinates p_kl(θ_s) of the dual Q's, lone Q³(θ_s − ch) → P³, P³′, M_vac → (5/2)h²W′(θ_s)P³(θ_s) with W(θ_s) = 0 (Wronskian zero), so the Gaudin ansatz is a bilinear form in degree-3 jet monomials with a state-independent kernel in (ω₁, ω₁′, θ₁₂); the Frobenius constraint P′(θ_s) = κ_sP(θ_s) (exponents {0,2,4}) makes the kernel unique only modulo a null space. gl(4) rectangular Gaudin = round 3 (plan §7), after the formula's structure is seen on sl(3).
**Agents (all sonnet; specs are self-contained):**
- **W2** (wolfbook-builder, K1): `Experiments/Gaudin-L2-chain-2-spec.md` → `Gaudin-L2-chain-2.wb`. Task A: the (∞,∞) family chain-side in both frames (r_k^∞, companion vacuum Ω_c^G[ω] = e₁ρ(U-factor of PauxG[ω]^{−T}) and its Gauss form, frame tie κ^G, exact undo, per-weight ranks = CG-completeness, K/S, CG, controls; A5 reference tables of Qnum(B) for the Python agents). Task B: the [3,1,0] fourth localiser (exact localising subspace 𝓛 ⊂ D₂, sparse representatives, the fourth covector, XXX (3,2) analogue). Task C: closed forms (c_H in ⟨Ω⊗Ω|H₁ = [Ω₁(c_H − C_ζ)]⊗Ω₁ over 8 reps — conjecture 1/θ₁₂ coefficient = 2λ₁; spectators; vacuum eigenvalue functions; universality of S3c[1,2]|_HW = (3/(2θ₁₂))[[3,−1],[−1,7]] and of the mixing roots R² − 5R − 5 = 0). No native outputs.
- **D1** (Python): `Paul/Python/gaudin_sl3/SPEC-twovac-fit.md` → `twovac_fit.py`, `RESULTS-twovac-fit.md`. Direct Gaudin fit of the jet-bilinear ansatz (hypothesis ladder: Q₁-jet × P³-jet basis 6/site → Plücker × P³ 12/site → larger; with/without M¹M²; five ζ-points for ζ-independence; rank/null space/sparse representative/rationality), kernel laws over a 3×3 (ω₁,ω₁′) grid and θ variants with the scaling-symmetry check, regime at 120 digits with CG from the Q side vs W1's exact values.
- **D2** (Python, exact): `SPEC-twovac-limit.md` → `twovac_limit.py`, `RESULTS-twovac-limit.md`. Formal h → 0 limit of the XXX ten-product formula in the diagonal frame (closed forms parsed from `ref/xxx-210-infinf-formula-cells.txt` = cells 1–70 extract of `L2-210-two-vacua-brackets`; κκ′ and x₀ leading orders; symbolic jets incl. the XXX Q's own h-corrections q₁, q₂; verdict whether q₁ jets enter at the leading order; if not, the explicit Gaudin formula + numerical test; if yes, the q₁ equation from the order-h⁴ Baxter expansion and stop).
- **P3** (Python): `SPEC-Qside-310.md` → `qside_310.py`, `gaudin_chain_rep.py`, `RESULTS-Qside-310.md`. Q1–Q3(–Q4) for [3,1,0] (gauge power of Π, F = Π²?, Lagrange columns, det N ratio identity on 225 states) and the (2,1,0) rank-4 puzzle (R).
**P3 RETURNED (2026-10-16 ~11:30; `Paul/Python/gaudin_sl3/RESULTS-Qside-310.md`, `qside_310.py` ~4 min cached / 15 min cold; 60 digits, 100 in Q4):** (3,1,0): HW index 14 (gl3_reps order), top state 224 = (6,2,0); cdet checks exact; no S₃ triple pole; **S₂ double-pole scalar c₂ = λ₂(λ₁+1)** (= 4 here, 3 for (2,1,0); checked on (1,1,0),(2,2,0),(3,2,0); NOT scalar for (n,0,0); (2,1,1) has an S₃ triple pole); **dual-oper local exponents {0, λ₂+1, λ₁+2}** (= {0,2,5}; {0,2,4} for the adjoint; 450/450); **direct oper = Π^p∘adjoint∘Π^{−p} with p = λ₁−λ₃ = 3** (p = 0,1,2 no solutions to degree 9; p = 4 = Π × the p = 3 solutions; direct Q exponents (0,3,5)); degrees dual μ_j, direct 6−μ_i in all 36 sectors; vacuum direct degrees (0,4,6), S₃ vacuum double poles 12/5 (ω), 36/7 (ω′); Wronskians: direct F = Π^{λ₁−λ₂} = Π², dual F′ = Π, c_ab = ζ_a−ζ_b, W = det[P,P′,P″] = ∏(ζ_j−ζ_i)e^{Σζu}Π^{|Λ|}; det N with the SAME six Lagrange columns (Π⁻³ weights) closes: ratio identity on all 225 states to 6.8e-52, 72/84 triples pass (permutation triples fail), Q4 on the 10 branches 6.6e-52, |det N| ~ ε^{16,15,14}. **No-log (Frobenius) conditions on the Laurent data:** 3 per site, e.g. a² − Ha + 4b = 0 (`frobenius_polys()`), hold on all 225 states × 2 sites — the Gaudin 'Bethe equations' for the Laurent data. **Rank-4 puzzle SOLVED:** rank 4 ⇔ Δμ = μ_B − μ_A has one positive and two negative components; mechanism: for every k with Δμ_k ≤ −1, ψ_k^A P^k_B is rational and decays at ∞, so Σ_s row(s,k,k) = 0 (residue theorem) — each site's rows span a 3-dim space, rank = 6 − dim(S₁∩S₂); the second null vector has no closed form (open, minor). Nilpotency: E^{λ₁−λ₃+1} = 0 (index 4 here).
**W2 RETURNED (2026-10-16 ~11:40; `Experiments/Gaudin-L2-chain-2.wb`, 80 cells / 21 code, NO native outputs; saved Λ = {2,1,0}, every code cell run twice with per-run markdown results `[2,1,0] run`/`[3,1,0] run`; K1 left in the [3,1,0] state):** **Task A — the Gaudin (∞,∞) family IS CG-complete, both reps**: r_k^∞ = exact Limit of the round-1 r_k (ΩD depends on ω₁, ω₂ only); the g′-Laurent coefficients are degree 1 in ω₂ and ΩDinf⊗ΩDinf is the unique common eigencovector of their top coefficients; **companion vacuum Ωc^G[ω] = e₁ρ(U of PauxG[ω]^{−T}) has Gauss parameters (1/ω₁, −1/(ω₁ω₂), (ω₁+ω₂)/(ω₁ω₂)) = the XXX Ωc[w] with w → ω, so Ωc∞ = gaussCov{a,0,a} with a = 1/ω₁ exactly**; frame tie ΩD·PphysG = κ^G·Ωc^G with PphysG = ρ((PauxG⁻¹)ᵀ) and κ^G = XXX kappaTie with P₃ = (PauxG⁻¹)ᵀ (Δ₁ = ω₁²/((ζ₂−ω₁)(ζ₃−ω₁)), Δ₂ = −ω₁ω₂/((ζ₁−ζ₂)(ζ₃−ω₁)(ζ₃−ω₂)), Δ₃ = 1/((ζ₁−ζ₂)(ζ₁−ζ₃)(ζ₃−ζ₂)); finite at ω₂ = ∞); exact undo at 11 points; a′ = (ζ₁−ζ₂)(ω₁−ω₁′)/((ζ₂−ω₁)(ζ₂−ω₁′)), c′ = (ζ₁−ζ₃)(ζ₂−ζ₃)(ω₁−ω₁′)/((ζ₁−ζ₂)(ζ₃−ω₁)(ζ₃−ω₁′)), g′ = (ζ₃−ζ₁)(ζ₃−ζ₂)(ω₁−ω₁′)/((ζ₃−ω₁)(ζ₂−ω₁′)(ζ₃−ω₁′)) (Jacobian rank 2); per-weight ranks = multiplicities; K/S = O(ε) all branches (max 2.5e-6 adjoint, 4.8e-6 [3,1,0] at ε = 1e-6); CG to 1e-18..1e-21; mixing roots R² − 5R − 5 = 0 and [3,1,0] (5,2,1) 4.3117/−0.8117 reproduced. **Control ω₁′ fixed FAILS at ω₂ = ∞** (φ_u rank 1 on multiplicity-2 weights; singlet-type branches K/S = O(1)) — unlike round 1's finite-ω₂ control. A5 reference tables: cells 51 ([2,1,0]) / 52 ([3,1,0]). **Task B — [3,1,0] fourth localiser FOUND:** D₂ has rank 26 because of two exact operator identities S3c[s,1] = ¼S3c[s,2](H_s − S3c[s,2]) (s = 1 with H₁; s = 2 with H₂); dim 𝓛 = 4; **X₄ = H₁² − ½H₁S3c[1,2] − (2/θ₁₂)W − S3c[1,1]** (mod id, H₁, S3c[1,2]; = (H₁−S3c12)(H₁−S3c12/4) − (2/θ₁₂)W), minimal supports {H₁², H₁S3c12} and {H₁², S3c12²}; site 2: X₄′ = H2t² + ½H2t·S3c22 + (2/θ₁₂)W − S3c21 with H2t = H₁+W; coefficients ζ,ω-independent; v₄ = Ω₁·(c₀ + c₁C_ζ + c₂S3[−2] − S3[−1] + C_ζ² + ½C_ζS3[−2]) is θ₁₂-INDEPENDENT (spectator constant); the 3-dim L₁ already reaches every weight (multiplicities ≤ 2) — the fourth is needed for ℓ = λ₁+1, not for reach; degree guess X₄ ↔ T₍₃,₂₎. **Task C — closed forms:** **c_H = |Λ|Σζ + (|Λ|² − Σλ_i²)/θ₁₂ = |Λ|Σζ + 2e₂(Λ)/θ₁₂** (10 reps; the 2λ₁ conjecture is wrong; mechanism Σ_ik ΩE_ik⊗ΩE_ki = (Σλ²)Ω⊗Ω); spectators f(H₁) = |Λ|Σζ − ζ·Λ + k/θ₁₂, f(H₂) = same − k/θ₁₂ (k = 4, 6), f(S3c[s,2]) = |Λ|ζ₃ (ω, θ independent); vacuum Laurent data ζ-independent: H_s^{vac} = Σ_i(|Λ|−λ_i)ω_i ± k/θ₁₂, a_s^{vac} = |Λ|ω₃, b_s^{vac} = ω₃(H_s − |Λ|ω₃); S₂ double pole c₂: W2 quotes |Λ|, P3 λ₂(λ₁+1) — equal on (2,1,0) and (3,1,0), P3 tested five reps, so λ₂(λ₁+1) is the general value; **universality: S3c[1,2]|_HW = (1/θ₁₂)×(θ-independent matrix, no scalar part)**: [3,2,1]: (3/(2θ₁₂))[[3,−1],[−1,7]]; (5,2,1): (1/(3θ₁₂))[[19,−40],[−2,44]]; (4,3,1): (4/(3θ₁₂))[[4,−1],[−2,14]]; coordinate-ratio equations r² + 4r − 1 = 0, 40r² + 25r − 2 = 0, r² + 10r − 2 = 0 — so the mixing (R² − 5R − 5 = 0) is universal. Pitfall: Laurent fits at θ₁ = 3/7 = u₀ hit a pole.
**D2 RETURNED (2026-10-16 ~12:10; `Paul/Python/gaudin_sl3/RESULTS-twovac-limit.md`, `twovac_limit.py` ~1 min) — THE GAUDIN TWO-VACUA FORMULA, case (i), no q₁ obstruction:** the full RHS of the XXX identity (incl. κκ′ ~ h⁻⁸, ⟨x₀^d|ψ⟩ = h²⟨x₀^G|ψ⟩ exactly with x(ζ) = (ζ₂−ζ₁, ζ₃−ζ₁, ζ₂−ζ₁, 2ζ₃−ζ₁−ζ₂, 3(ζ₂−ζ₁), 3(ζ₃−ζ₁), 2(ζ₃−ζ₂), 3(ζ₃−ζ₂)) in notebook order, M_vac ~ h²) has leading order h¹⁰; orders h⁴..h⁹ cancel identically because of ONE relation: the dual Baxter equation at x = θ_s has ν₃ = t₃ = 0, so Q^j(θ_s−h) = ρ_sQ^j(θ_s) with ρ_s j-independent and minor_s(−1,0) = 0 to all orders (generalises W(θ_s) = 0); the h¹⁰ coefficient involves only n = 0 jets (q₁ jets enter at h¹¹). **Formula (verified on all 64 states to 3e-55 at 3 (ω₁,ω₁′), 2 θ, 2 ζ; constant exactly 1):** with jets at θ_s, κ_s = P³′/P³ = P²′/P², p_kl = P²^{(k)}P³^{(l)} − P³^{(k)}P²^{(l)} (p₀₁ = 0, p₁₂ = κp₀₂), the five per-site monomials N = (P³p₀₂, P³′p₀₂, P³p₀₃, P³′p₀₃, P³p₂₃): **4c·Qnum_G·(P³p₀₂)(θ₁)(P³p₀₂)(θ₂) = Σ_ij G_ij(ω₁,ω₁′,θ₁₂) N_i(θ₁)N_j(θ₂)**, Qnum_G = ⟨Φ∞|Ψ⟩/⟨x₀^G|Ψ⟩, c = (ζ₂−ω₁)(ζ₃−ω₁)²(ζ₂−ω₁′)(ζ₃−ω₁′)²(ζ₁−ζ₂)² (= κ^Gκ^G′: ζ enters only here); G 5×5, rank 5, site-exchange symmetric, homogeneous of degree 10 − w_i − w_j (weights 2,3,3,4,5), G[P³p₂₃|P³p₂₃] = 1, G[P³p₀₃|P³p₀₃] = ω₁²ω₁′²; the 64×25 design matrix has full rank 25 ⇒ G unique in this basis. Equivalent forms: Qnum_G = (1/4c)ΣG_ij ν_i(θ₁)ν_j(θ₂), ν = (1, κ, π, κπ, ε), π = p₀₃/p₀₂, ε = p₂₃/p₀₂; twisted jets (e^{−ω_s u}P^j) give a kernel in δ = ω₁′−ω₁ and 1/θ₁₂ only; direct basis p₀₂ = A_sQ₁, p₀₃ = A_s[(Λ_s−κ)Q₁ + 2Q₁′] (A_s = (ζ₃−ζ₂)e^{σθ_s}Π′(θ_s), Λ_s = 2σ + 2/Π′(θ_s)) ⇒ needs only Q₁, Q₁′ × P³ jets with κ_s as a symbol (9 monomials/site). M_vac^G = 2h²W′P³ (not 5/2: the weight-2 minor carries an n = 1 piece fixed by minor(−1,0) = 0). **L4: the naive fixed-a (a = 1/ω₁) t → ∞ limit of the XXX kernels is WRONG** — rank-one U⊗U, state-independent (1−a)³(1−a′)³; the scaling limit a = e^{−hω₁} → 1 is a different limit and all state dependence sits in the cancellation hierarchy. Not done: the h¹¹ layer (q₁), general proof beyond exact random-jet series + symbolic effective model, [3,1,0].
**D1 RETURNED (2026-10-16 ~12:30; `Paul/Python/gaudin_sl3/RESULTS-twovac-fit.md`, `twovac_fit.py` ~2.5 min warm) — INDEPENDENT CONFIRMATION OF D2:** with target 4c·Qnum·(P³p₀₂)₁(P³p₀₂)₂ and the five-monomial basis N₅ the 25-parameter kernel fits at residual < 1e-90 on 320 rows (5 ζ-points), rank 25, unique, all entries rational, **equal to D2's G to 2.8e-74**; per-ζ fits agree with the stacked one to 1e-72 (no free scalar needed once c and M = P³p₀₂ are used); the 64 Qnum at P match W2's A5 table to 3e-31 (printed digits). Kernel laws: every G_ij is a POLYNOMIAL in (ω₁, ω₁′) of total degree ≤ 6 (no poles) and in 1/θ₁₂ of degree ≤ 4, identical to D2's in 25/25 entries; scaling test at c = 2 confirms homogeneity 10 − w_i − w_j with weights (2,3,3,4,5). **κ_s = a_s/3** (the Frobenius constant P′/P at θ_s is one third of the (u−θ_s)⁻² coefficient of S₃; it does not involve H_s or the regular parts). Direct-Q₁ bases do NOT fit with a ζ-independent kernel (p₀₃ = A_s[(Λ_s−κ)Q₁ + 2Q₁′] with Λ_s ∋ σ = Σζ) — the Plücker N₅ form is the ζ-free one. T/Qnum is NOT constant: ⟨x₀^G|Ψ⟩/n_B ∝ Q₁(θ₁)P³(θ₁)Q₁(θ₂)P³(θ₂) exactly, so the x₀-free statement is ⟨Φ∞|Ψ⟩·D_B(ω₀)/⟨Ω₀⊗Ω₀|Ψ⟩ ∝ Σ_ij G_ij N_i(θ₁)N_j(θ₂) with NO vacuum-monomial denominators (P1's det N anchor absorbs them). Controls fail (lone P² or P¹: 0.99; jets at wrong points 0.48; no M: 0.97). Larger bases (144, 576 unknowns) fit with null spaces whose sparse representative is the N₅ kernel. **Regime (120 digits):** the Q-side formula equals the chain Qnum on the six HW branches to 1e-90..1e-98 at ε = 1,2,4e-6 (cancellation 13–31 digits); **CG from the Q side alone** via Richardson ratios ρ_a = Qform_a/Qform₀ reproduces W1's exact values (−1/√2, −1/√2, 0.1178511302, both octet vectors) to the printed digits, Q-side vs chain ρ_a to 12 digits. Caveat: a numerically vanishing feature (p₀₁ ~ 1e-78) must be zeroed explicitly before QR; some ζ triples give a complex S₂ + πS₃ spectrum (skipped).
**ROUND-2 RECONCILIATION (zero discrepancies):** D1 = D2 on the kernel, the constant 4c, the vacuum monomial P³p₀₂ and the kernel laws; W2's A5 Qnum tables = D1's to the printed digits; W2's frame-tie scalar satisfies κ^Gκ^G′ = ω₁³ω₁′³/c (D1) while the XXX κκ′ = h⁻⁸/c (D2) — all consistent; P3's rules (gauge λ₁−λ₃, exponents {0, λ₂+1, λ₁+2}, c₂ = λ₂(λ₁+1)) hold on both reps. The round-2 plan §1 expectations were right in structure and wrong in two details, both corrected there: only FIVE monomials per site (not 12) and M_vac = 2h²W′P³ (not 5/2).
**Not done / next (round 3 candidates — Paul to choose):** (1) Paul reruns BOTH notebooks (`Gaudin-L2-chain.wb`, `Gaudin-L2-chain-2.wb`) in VS Code to store native outputs; (2) the jet formula for [3,1,0] (exponents {0,2,5} ⇒ jets to order 4, Plücker p_kl with l ≤ 4, lone P³ to order 2?; the XXX side had 45 products there — expect the Gaudin monomial count to be fixed by the exponents; D1's fit pipeline is rep-generic once the Q's are, P3's `qside_310.py` has them); (3) closed-form derivation of G (D2's 25 polynomial entries are reconstructed, not derived — look for a determinant/Wronskian structure: G rank 5, weights (2,3,3,4,5), G[P³p₂₃|P³p₂₃] = 1); (4) **gl(4) rectangular Gaudin** (the project target; plan `Gaudin-sl3-round2-plan.md` §7: chain template `new-conventions/L2-gl4-rect-untwist.wb`, generators in `old/Gaudin/Clean/Hamiltonians_SU4.wb`, Python gl(4) generators needed); (5) sl(3) Gaudin section of `paper/gaudinCG.tex` (conventions, localisers, untwist, det N, the jet formula — nothing of this is in the paper yet); (6) promote the two notebooks to Clean/ + docs agent; (7) minor: the h¹¹ layer, the second null vector n₂ of the rank-4 pairs, the λ₃-dependence of the gauge rule (needs a rep with λ₃ > 0). Nothing committed; `Paul/Mathematica/Gaudin/` and `Paul/Python/gaudin_sl3/` untracked (the latter holds ~50 MB of .pkl caches — consider a .gitignore for `*.pkl` before committing).

## (2026-10-15/16): sl(3) GAUDIN ROUND 1 — the XXX "localising T's → CG" program DONE for the two-site gl(3) Gaudin model at layers A–C ([2,1,0] and [3,1,0]); all three agents returned 2026-10-16 and agree exactly
**Paul's ask:** repeat the localising-T / two-vacua-untwist / CG pipeline for the **sl(3) Gaudin model** (not gl(4)); "I don't know which T's to use".
**Plan (read first):** `Paul/Mathematica/Gaudin/Gaudin-sl3-CG-plan.md` — §1 is the binding XXX ↔ Gaudin dictionary (Lax L_g = g + ΣE^{(s)}/(u−θ_s); commuting family = coefficients of the column-ordered cdet(∂−L_g) = ∂³−S₁∂²+S₂∂−S₃; the "T's" = Laurent coefficients of S₂, S₃ at θ₁, θ₂; Lie-algebra companion frame g_c(ζ) since the XXX G_c(z) ≠ exp(h g_c); vacuum Ω_{ζω} = e₁ᵀρ(n) with the SAME r_k formulas (degree-0 homogeneity ⇒ exactly ε-independent in the regime ζ = εφ, ω = εω̂); brackets = residues at θ_s of pole order λ₁−λ₃+1; Q = e^{±ζu}·polynomial; Gaudin = XXX at t = θ₁₂/h → ∞). §2–§4 = the three agent stages, §5 later stages, §7 risks.
**Agents (first dispatched 2026-10-15 ~20:00, lost with that session; re-dispatched 2026-10-16 ~09:20, all sonnet; ALL THREE RETURNED by ~10:10, reports reconciled below and folded into the plan's §1):**
- **W1** (wolfbook-builder, shared kernel K1): `Paul/Mathematica/Gaudin/Experiments/Gaudin-L2-chain-spec.md` → `Gaudin-L2-chain.wb` (Tasks 0–5: cdet Hamiltonians + commutativity, vacuum eigencovector with g′ = P_aux(ζ)g_c(ω)P_aux(ζ)⁻¹ (variants if it fails), eigenstates/HW-limit branches, the localising rank-1 scan in D₂ (answers "which T's"), the two-vacua family Ω_{ζω}⊗Ω_{ζω′}: exact undo, K/S, CG at leading order vs exact HW kets; then [3,1,0]). JSON-built notebook, cells run in K1, no native outputs.
- **P1** (Python): `Paul/Python/gaudin_sl3/SPEC-Qside.md` → `gaudin_chain.py`, `qside.py`, `RESULTS-Qside.md` (opers: which ODE is direct/dual, polynomial degrees; Wronskian/QQ; residue brackets, functional orthogonality, 6×6 det N and the x₀-free overlap-ratio identity on all 64 states).
- **P2** (Python): `Paul/Python/gaudin_sl3/SPEC-xxx-limit.md` → `xxx_chain.py`, `xxx_limit.py`, `RESULTS-xxx-limit.md` (XXX port at z = e^{hζ}: frame/vacuum limits, the t_a ↔ S_k dictionary, h-expansion of the localised covectors ⟨Ω|T_λ^{(1)}(θ₁−λ₁h) = predicted Gaudin localisers, Baxter → oper).
- Reference cell dumps (UTF-8) for all three: `Paul/Mathematica/Gaudin/ref/ref-*.txt` (Stage D untwist notebook, 210-untwist-CG, L=1 clean, rotation-sov, two-site brackets, two-vacua brackets, Minimal-ingredients-EH, diagonal-frame-overlap-EH).
**P1 RETURNED (2026-10-16, Q side closes at 55–67 digits; `Paul/Python/gaudin_sl3/RESULTS-Qside.md`):** cdet column-ordered with E_ij commutes (row-ordered E_ij and column-ordered E_ji do NOT); S3 has NO triple pole; g′ = P g_c(ω) P⁻¹ (first variant) is the one that works; dual oper = (∂³−S1∂²+S2∂−S3)P = 0 with P^j = e^{+ζ_j u}p_j, deg μ_j; **direct oper = Π²∘(formal adjoint)∘Π⁻², Π = (u−θ1)(u−θ2)** (plain adjoint has NO e^{−ζu}×poly solutions) — gauge power 2 fixed for [2,1,0] only; Wronskian Q_aQ_b′−Q_bQ_a′ = (ζ_a−ζ_b)e^{−Σζu}(u−θ1)(u−θ2)P^c; quantum-Weyl S_k = ratios of det[P,P′,P″]-type dets; orthogonality is the Lagrange identity P(ΔS2ψ)′+ΔS3Pψ = C′ with ψ = Π⁻²Q_i, giving six residue columns N[.,H_t], N[.,a_t], N[.,b_t] (state-dependent Δ's — the spec's fixed-coefficient search cannot work); 6×6 det N rows (s,(i,j)) closes, x₀-free ratio identity holds on all 64 states (6e-55) and on the six HW branches at ε = 1e-6 (5e-67); HW ket index is 7 (0-based) in gl3_reps order; open: 12/60 different-sector pairs have rank 4 (unexplained), why no ΔS1/Δe2/Δe3 columns are needed for mixed twists; [3,1,0] not done.
**P2 RETURNED (2026-10-16, XXX→Gaudin limit, exact h-series, 0 failures; `Paul/Python/gaudin_sl3/RESULTS-xxx-limit.md`):** AGREES with P1 on cdet (plan's column-ordered E_ij; row-ordered E_ji = transposed mirror also fits), on g′ = P_aux(ζ)g_c(ω)P_aux(ζ)⁻¹ ((G′−1)/h → g′ at O(h), Jordan block of G_c harmless), and on the oper roles: dual Baxter → Talalaev oper on Q^j = e^{+ζ_j x}q^j, direct Baxter → formal adjoint acting on Q_i/ν², ν = ∏(x−θ_s) (= P1's Π⁻² gauge, found independently). Dictionary: A1 = S1, B1 = 2S1 (forced by Baxter at h¹, h²), S2 = −A2 + S1²/2 + S1′/2, S3 = A3 − B3 + S1³/6 + S1S1′/2 + S1″/2, B2 = A2 + S1²/2 + S1′/2; S_k(∞) = e_k(ζ). **x₀ IS parallel** to x₀^G(ζ) (plan §0.3 expected not; ‖x₀‖ ~ h per site). **Predicted site-1 localiser:** T₍₁,₁₎(θ₁−2h) = 3h²θ₁₂²[1 + h X₁/(…)], X₁ = −(4631/124)·1 + 3H₁ − S3c[1,2] EXACTLY (in span{1,W,H₁,S3 coeffs}, rank 7 of 9), ⟨Ω⊗Ω|X₁ = u⊗Ω with u = Ω·(431/62 − 3C_ζ − S3^{1-site}[−2]); one-site v₁(T₁₁) = Ω·(215/8 − (9/2)C_ζ − S3[−2]), v₁(T₂₁) = Ω·(−215/8 + (9/2)C_ζ) (no S3 term); spectators → θ₁₂², θ₁₂³. Site 2: ⟨Ω|t₂(θ₂+h)t₁(θ₂−h) has valuation h⁴ with leading term Ω·(3C_ζ − 41/4) NOT ∝ Ω. Constants only at the one numerical point (no closed forms in θ₁₂, ζ, ω). xxx_chain re-orders gl3_reps to notebook order (index 0 = HW).
**W1 RETURNED (2026-10-16; `Paul/Mathematica/Gaudin/Experiments/Gaudin-L2-chain.wb`, 48 cells / 16 code, NO native outputs — Paul must rerun in VS Code; saved with Λ = {2,1,0}, [3,1,0] run from a copy and its tables appended to the markdown result cells):** Tasks 0–5 all pass for both reps. cdet convention of the plan works first time (7 commutator pairs exactly 0). **Sign correction:** H₁ = Res_{θ₁}S₂ = −[Σζ_iE1_ii + Σ_{ik}E1_ikE2_ki/θ₁₂] + (|Λ|Σζ + |Λ|²/θ₁₂)·id and H₁+H₂ = −W + 2|Λ|Σζ (plan had +). g′ = PauxG[ζ].gc[ω].PauxG[ζ]⁻¹ exact eigencovector (first variant; all others fail); one-site: e₁ρ(U-factor of M^{−T}) = ΩD[ζ,ω] exactly; P_phys^G = (PauxG⁻¹)ᵀ. **THE ANSWER TO "WHICH T's" ([2,1,0]):** ℓ₁ = ℓ₂ = 3, inside D₁ already — site 1: {id, H₁, S3c[1,2]} EACH SEPARATELY localises (⟨Ω⊗Ω|H₁ = [Ω₁(|Λ|Σζ + 4/θ₁₂ − C_ζ)]⊗Ω₁, ⟨Ω⊗Ω|S3c[1,2] = [Ω₁·S3^{1-site}[−2]]⊗Ω₁ exactly), L₁ = Ω₁·span{1, C_ζ, S3^{1-site}[−2]}, full rank on every weight of V; site 2: {id, W+H₁ = −H₂+scalar, S3c[2,2]}; spectators f(H₁) = 6 + 4/θ₁₂ (P), f(S3c12) = 2 (Laurent polynomials in θ₁₂, point-dependent constants). P2's predicted X₁ = −(4631/124) + 3H₁ − S3c[1,2] reproduced exactly (it is a combination of the three). **[3,1,0]: ℓ₁ = 4 = λ₁+1** (3 over D₁, 4 over D₁+H₁·D₁ and D₂) — the fourth localiser needs ≥ 2 quadratic products, explicit form NOT identified (open). D₂ rank 24 ([2,1,0]) / 26 ([3,1,0]); cubic span 48 / 68 (not vacuous). Multiplicity-2 sectors: S₂(u₀) is scalar on the HW space, S3c[1,2] splits the copies ([3,2,1]: (105/62)[[3,−1],[−1,7]]); the Gaudin untwisted copies = its eigenvectors. **Task 4 payoff: K/S = O(ε) for EVERY HW branch** (max|K/S| at ε = 10⁻⁶: 4.9e-7 … 5.4e-6 for [2,1,0], ≤ 7.7e-6 for [3,1,0]; doubles with ε; S = O(1); K ≡ 0 top branch), exact undo on all HW kets at 9 points, a′,c′,g′ = the XXX leading coefficients exactly, per-weight ranks = multiplicities (CG-complete) both reps; CG from the overlap alone matches the exact G-norm-1 HW slices to 1e-17…1e-19 after Richardson (1-dim slices −1/√2, −1/√2, 0.1178511302; octet mixing R_c = roots of R² − 5R − 5 = 0, i.e. 5.854101966 and −0.854101966, agree to 1e-18; [3,1,0] (5,2,1): 4.311737691 / −0.811737691). **Surprise:** the control "ω′ fixed (not untwisted)" does NOT fail in Gaudin (ΩD at ζ → 0 has a finite direction-dependent limit; K/S still O(ε), ranks full) — unlike XXX. Common ω̂₁′ = ω̂₁ control behaves as expected (a′ = g′ = 0, only (4,2,0),(4,1,1) read). Task 4(f) Q-side numbers in cell c14: equal-twist overlaps ∝ ε (ε² singlet), mixed ⟨Ω_ω⊗Ω_ω′|Ψ⟩ = −⟨Ω_ω′⊗Ω_ω|Ψ⟩ = O(1) for octets, (4,1,1), (3,3,0). Open: [3,1,0] fourth localiser's explicit form; mixing functionals for the [3,1,0] (4,3,1) sector; a both-fixed degenerate control.
**Reconciliation (all three agree, zero discrepancies):** cdet column-ordered with E_ij; S3 no triple pole; g′ first variant; dual oper = Talalaev on e^{+ζ_j u}q^j, direct = adjoint on Q_i/Π² (P1 and P2 independently); vacuum eigenvalue functions and S3 double poles match across W1/P1; X₁ localiser prediction (P2) matches the scan (W1). Only the plan's H₁ sign and its "x₀ not parallel" expectation were wrong (both corrected in plan §0/§1). HW-ket index differs between gl3_reps order (7) and notebook order (0) — P1 uses 7, P2 re-orders, W1 is notebook-native.
**Not done / next:** (1) Paul reruns `Gaudin-L2-chain.wb` in VS Code to store native outputs and reads the summary cell (c47); (2) plan §5: the Q-side formula for the two-vacua overlap ⟨Ω_ω⊗Ω_ω′|Ψ⟩ (P1's det N with the Lagrange-identity columns handles equal twists; mixed twists worked with the same six columns in P1's Q3(c) — unexplained, derive why), via t → ∞ of the XXX 10-product kernels and/or directly; (3) the [3,1,0] fourth localiser (two-product combination) and the (4,3,1) mixing functionals; (4) closed forms in (θ₁₂, ζ, ω) for the point constants (4631/124, 431/62, 1831/124, …) — spectators are already Laurent polynomials in θ₁₂; (5) P1's rank-4 pairs (12/60 different-sector) explained; (6) sl(3) section of `paper/gaudinCG.tex` (not touched this round); (7) promote the notebook to Clean/ + docs agent once Paul is satisfied. Nothing committed; `Paul/Mathematica/Gaudin/` and `Paul/Python/gaudin_sl3/` are untracked.

## >>> (2026-10-11): two-site bracket plan, Stages A, B, D DONE — Ω_w⊗Ω_{w′} untwists cleanly (ALL CG data at leading order, both reps); SoV construction site-asymmetric, cured by the R-matrix; adjoint (∞,∞) two-vacuum family = 10-product bracket formula with all six CG data from the Q side; [3,1,0] stays at 45 <<<
**Plan being executed:** `new-conventions/Two-site-rotation-brackets-plan.md` (Stages A–D, sonnet wolfbook-builders, sequential on K1). **Stage A spec:** `L2-site2-rotation-brackets-spec.md`; **notebook:** `L2-site2-rotation-brackets.wb` (36 cells, [2,1,0], [3,1,0], [3,2,0]; result cell 36; driver cell 32 wipes Global` — run last; **no native outputs**, cells run in K1 from their saved sources, results in markdown cells; file built as JSON by the agent and opened via the `code` CLI). Not committed. **Stage B spec:** `L2-two-site-rotation-brackets-spec.md` (two-site K tensor, (c,c′) blocks, det₆ regression, Ωc[w]⊗Ωc[w′], R-dressing test) — dispatched next.
- **Verified exactly (plan §0):** Ωc[{w₁,w₂}] = e₁ᵀρ(g) with (α,β,γ) = (1/w₁, −1/(w₁w₂), (w₁+w₂)/(w₁w₂)), first row (1, 1/w₁, 1/w₁²), so the [3,1,0] 5-bracket locus β+αγ = α² IS the one-site vacuum family; Ωc[w]⊗Ωc[w] is the exact left eigen-covector of t_a with companion twist G(w₁,w₂,w₃) for ANY w₃. The two [3,1,0] closed-form λ₂-group measures read w₁^{−(u−θ₁)/h} × SoV measure × 1/(u−θ₂+h) × polynomial, amplitudes 48(t−3)(t+1)/(t w₁⁸) and −12(t−3)(w₂+3w₁)/(w₂w₁⁹), end-point exception ω(θ₁−h) = 3w₁/(w₂+3w₁) (cells 13–14).
- **Adjoint quadric = torus image of the vacuum family (exact, cell 30):** e₁ᵀρ(g)ρ(D) ∝ e₁ᵀρ(D⁻¹gD) rescales (α, β+αγ) by (d₂/d₁, d₃/d₁), so {β+αγ = cc·α²} = {Ωc[w]ρ(D) : d₁d₃/d₂² = cc}; cc₁ = (t−3)/(2(t−2)) = (3/4)κ_singlet. The §5 product point's site-1 factor has non-constant b₁b₃/b₂² (on the quadric only along a curve); the diagonal-frame resolvent ratio (y−z₂)²/((y−z₁)(y−z₃)) = cc₁ at no distinguished y — no identification.
- **SURPRISE — site asymmetry (cells 16–24):** e₁ᵀ⊗e_νᵀ (site-2 rotation) expands on ALL dim V site-1 labels (35/64, 103/225, 122/225 SoV states for [2,1,0]/[3,1,0]/[3,2,0]; vacuum label on either site pairs with everything; no weight rule). Per-label automatic counts **17 / 40 / 54** vs the site-1 **5 / 7 / 11**; block form (ladder-1 2-forms ⊗ ladder-2 2-forms, matrix-rank decomposition) ranks 3,3,3,2 → 15 products ([2,1,0]) and 6,6,6,4 → 35 products ([3,1,0]). The vacuum family does NOT reduce the site-2 count; only b₃ = β+αγ = 0 does (13/28/43) and it is not CG-complete (5/7, 7/12, 9/12 weights); CG-complete reducing loci save 1–2 brackets. Adjoint site-2 conic: **cc₂(t) = (t+2)/(2(t+1)) = cc₁(1−t)** (15 (h,θ) points). [3,1,0] site-2 conic constants (t-laws): (2t²+5t−4)/(2(t−1)(t+2)), (t−1)(8t²+13t−7)/(6(t−2)t(t+1)), (t²+3t+6)/(3t(t+1)), (t+2)/(2t), 2(t+3)/(3(t+1)). Reach: generic/vacuum-family/conics all weights, rank 2 on doubles; undo identity ⟨e₁ᵀ⊗g₂|h⟩ = ⟨e₁ᵀ⊗e₁ᵀρ(g₂L₀D₀)|h⟩ exact on HW kets. Closed forms on b₃ = 0: wedges with exponent law gA^{m+m′+k}, kernels rational in t with poles at t ∈ {0,±1,±2}; no natural u-form.
- **Reading (main session, to be tested in Stage B Task 2b):** the reversed-order chain ℒ^{(2)}ℒ^{(1)} = R T R⁻¹ with the Yangian R^{ΛΛ}(θ₁₂) on V⊗V, so t̃ = RtR⁻¹, Ψ̃ = RΨ, identical SoV wavefunctions; mirror symmetry then says the cheap site-2 covectors are the R-dressed ones ⟨(e₁ᵀ⊗g₂)R| with the site-1 counts and constant cc₁(−t); the 1−t of the undressed family would be the R shift. Two-site candidate: ⟨(e₁ᵀ⊗g₂)R(ρ(g₁)⊗1)|.
- **Stage B DONE (`L2-two-site-rotation-brackets-spec.md`, notebook `L2-two-site-rotation-brackets.wb`, 55 cells, summary cell 55, driver cell 52 wipes Global`; no native outputs; not committed).** Two-site K tensor exact (225×225, 2051 nonzero rationals, denominators ≤ 12 digits, z-independent; slices reproduce the one-site maps). **R-dressing CONFIRMED exactly (cells 30–36):** the Yangian R^{ΛΛ}(θ₁₂) on V⊗V (1-dim solution of R·T = T̃·R; commutant dim 8/14 needs two quartic invariants besides 1, C₁₂, D₁₂, D₂₁) satisfies t̃ = RtR⁻¹, B̃R = RB, (e₁⊗e₁)R = e₁⊗e₁; **(e₁ᵀ⊗g₂)·R has support only on site-1 labels (λ₁|λ₁c₁), counts 5/7, adjoint conic cc₁(−t) = (t+3)/(2(t+2)) (167/264 at t = 62/35; checked at 3 more points), [3,1,0] best locus β+αγ = α² with count 5.** Dressing the cheap site-1 rotation makes it expensive (17/40). **Block counts (lower bound Σr / product count ΣP, [2,1,0] / [3,1,0]):** generic pair 14/23, 27/54; **vacuum Ωc[w]⊗Ωc[w] (= det₆) 13/18, 26/40**; **Ωc[w]⊗Ωc[w′] 13/22, 27/54 — generic for [3,1,0]**; R-dressed two-site candidates (e₁⊗Ωc[w′])R(ρ(g₁)⊗1) and the other order: generic 23/54. Sub-families: common w₁′ = w₁ → 18/40 (= the vacuum count); w₂′→∞ → 14/45; both w₂,w₂′→∞ → 8/14 exact. No count is a product of one-site counts; the vacuum φ-matrix has rank 12 < 15 for [3,1,0]. One-site-fixed loci depend on the other site (α = α′, β = 0 …): the sites do not decouple. **Reading:** cheapness correlates with the pair becoming a near-global rotation (common w₁ differs only in E₂₃ ⇒ reads one irrep; w₂ = w₂′ = ∞ is a global rotation ⇒ top irrep only), so no cheap AND CG-complete two-vacuum family exists in this cost metric; the generic count 54 is nevertheless polynomial in the boxes by construction ((λ₂+1)² blocks of size (λ₁+2 choose 2)). The det₆ Laplace bound Σr ≤ 20 is NOT valid (Σr = 26 for the vacuum). Agent's Stage C recommendation (common-w₁ family, 18/40) is overruled by the reach argument — to be confirmed by the finite-twist reach in Stage D.
- **Stage D (chain side) DONE — POSITIVE (`L2-two-vacua-untwist-spec.md`, notebook `L2-two-vacua-untwist.wb`, 40 cells, summary cell 40, saved with Λ = {2,1,0}; no native outputs; not committed).** Family Ω_w⊗Ω_{w′} (diagonal frame, Ω_w = e₁ᵀρ(n(z,w))), regime z = 1+εφ, w = 1+εω, w′ = 1+εω′ all independent. Exact undo ⟨Ω_w⊗Ω_{w′}|h⟩ = ⟨e₁ᵀρ(n(w)n(w′)⁻¹)⊗e₁ᵀ|h⟩ on all HW kets (7 points). **n′ = exp(a′E₁₂)exp(g′E₁₃)exp(c′E₂₃) with all three O(1) at ε = 0:** a′₀ = (φ₁−φ₂)(ω₁−ω₁′)/((ω₁−φ₂)(ω₁′−φ₂)), g′₀ = (φ₁−φ₃)(φ₂−φ₃)(ω₁−ω₁′)(ω₂−ω₁′)/((φ₃−ω₁)(φ₃−ω₂)(φ₃−ω₁′)(ω₁′−φ₂)) (ω₂′-independent), c′₀ = r₃⁰(ω)−r₃⁰(ω′); rank 3. Structure: b = r₂+r₁r₃ = (w₁−z₁)/(w₁−z₃) is w₂-independent, so g′ = (b−b′)−(a−a′)c; common w₁ ⇒ a′ ≡ g′ ≡ 0 (reads one irrep only). Pinning w₂′ (at y₀ or anything) does NOT degrade g′: the old det₆×det₂ disease came from slaving the second vacuum to the first, not from a pinned y per se. **Per weight: everything at ε⁰ with rank = multiplicity (both reps, also for fixed w′). K/S = O(ε) for EVERY HW branch** (adjoint 6 incl. singlet 2.6·10⁻⁶ and both octets; [3,1,0] all 10 incl. [4,4,0], [4,3,1]×2, [4,2,2], [3,3,2], max 2·10⁻⁵ at ε = 10⁻⁶, doubling with ε). **All CG data at leading order from the overlap alone: 6/6, 10/10**, errors ∝ ε, Richardson overlap-vs-chain 10⁻¹⁶–10⁻²⁰; chain limits reproduce 46.67453891, 10.92853892, −5118.538314, 184.3305612/−238.5587523, R = −1.694/1.096 (adjoint) and 46.6745, −36.8838, −1325.402367/−334.3316172 ([3,1,0]); new [3,1,0] slices [4,4,0] −115.8530207, [3,3,2] 23369.80564 (cell 33). **Finite-twist reach (cell 21):** generic full; common w₁ → 2 weights (no); **w₂ = w₂′ = ∞ with w₁ ≠ w₁′ → FULL (CG-complete, 2 parameters)**; w₂ = w₂′ = ∞ with w₁ = w₁′ = global rotation (top only); w₂′ = ∞ only → full. Frame tie Ωc[w] ∝ Ω_w·P_phys up to the scalar κ = Π_k Δ_k(n·P₃)^{λ_k−λ_{k+1}} (Ωc normalised to first component 1). Ω_w⊗Ω_w is the eigen-covector of t_a for the bra twist G′(w) = P_aux G(w) P_aux⁻¹, NOT diag(w). Task 5 (cell 37/39): the four overlaps ⟨ww|, ⟨w′w′|, ⟨ww′|, ⟨w′w| per HW branch with ε-exponents (mixed exponent = CG-slice exponent; ⟨ww′| = ±⟨w′w| at leading order on 1-dim slices, sign = exchange parity), also for the (∞,∞) sub-family; det₆ identity for the diagonal ones 10⁻¹¹⁶ / 10⁻⁷⁷.
- **Stage B follow-up DONE (cells 47–49 of `L2-two-site-rotation-brackets.wb`, now 58 cells; summary cell 58):** the CG-complete (∞,∞) family Ωc[(w₁,∞)]⊗Ωc[(w₁′,∞)] = gaussCov[{a,0,a}]⊗gaussCov[{a′,0,a′}], a ≠ a′: **[2,1,0] ΣP = 12 (ranks 3,3,2,2, Σr = 10; generic 23, det₆ 18) — cheap; [3,1,0] ΣP = 45 (ranks 6,6,6,6, Σr = 24; generic 54, det₆ 40) — NOT cheap.** Rank drops on a = a′ (global rotation, 8/14), a = 0 (8/19), a′ = 0 (trivial). Asymmetric mixes: gaussCov[{a,0,a}]⊗Ωc[w′] 20/54, Ωc[w]⊗gaussCov[{a′,0,a′}] 14/45. **Conclusion for [3,1,0]: no two-vacuum family found that is both CG-complete and near the ~25 target in the product-of-2×2-brackets metric** (all CG-complete ones cost 45–54; the det₆ itself costs 40). Decision point for Paul: accept the generic block formula (polynomial in boxes by construction, 54 products for [3,1,0]), change the cost metric (bracket determinants of size > 2; the det₆ packs 36 bracket entries), or look for a determinant form of ⟨Ω_w⊗Ω_{w′}|Ψ⟩ with non-vacuum bra data now that the exact two-site K tensor is available.
- **Stage C + Q-side payoff for the adjoint (∞,∞) family DONE — POSITIVE (`L2-210-two-vacua-brackets-spec.md`, notebook `L2-210-two-vacua-brackets.wb`, 70 cells; tasks in order 0, 1, 3, 2, 4; drivers cells 59 and 67 wipe Global` — run last; no native outputs; not committed).** **Count 10 = Σr** (the 12 of Stage B was the default-pivot artefact; pivot scan cells 26–30; all 20 ladder 2-forms decomposable, identity 10⁻⁴⁴ on 64 states at three generic points — NB (2/3,−5/7) has w₁ = z₁, a degenerate vacuum with 22 vanishing overlaps). **Torus scaling law (exact, all t):** A^{(cc′)}_{(mm′),(nn′)} = a^{m+m′+e} a′^{n+n′+e′} κ^{(cc′)}(s), s = a′/a, (e,e′) = (7,5),(7,4),(6,5),(6,4); F_m = a^m f_m(s), G_m = a^{m+e} g_m(s) on ladder 1, same with a′ on ladder 2 — the two vacua couple only through s; of 160 vector entries 92 vanish, 33 constant, 31 monomials in s, 4 non-monomial; a−a′ appears as (s−1) factors (rank 3 → 2 at a = a′). Closed forms: cell 29 (s-frame, t = 62/35), 31 (raw in a,a′), 62 (exact in s,t; 40 kernel laws κ(s,t) fitted on 11 t-values, 2 held out exact; kernel poles only t ∈ {0, ±1}; K tensor homogeneity degree 0 in (h,θ), depends on t only). **u-forms: none non-trivial** (cell 64: 18 three-point planes in no tested class; each has one natural vector measure×(m−m*(s,t)); 8 contain a delta). **Q-number (cells 50–52):** companion/diagonal dictionary ⟨Φc|Ψc⟩/⟨x₀c|Ψc⟩ = ⟨Φd|Ψd⟩/(κκ′) (10⁻¹²³; κ, κ′ ~ ε⁻⁴); Qnum = overlap on all 21 HW-sector states to **7.9·10⁻¹¹¹**; ε-exponents: entries −1..1, Bk₂'s 0..2, every product and block sum O(1), total Qnum ε^{6,5,5,4,4,2} for (4,2,0),(4,1,1),(3,3,0),(3,2,1)×2,(2,2,2) — the cancellation is inside the Bk₂'s (up to 26 digits) and in the final sum over the four blocks (11–37 digits). **All six CG data from the Q side alone (cells 54–56):** Richardson Q vs chain ≤ 1.2·10⁻¹⁸: (4,2,0) 2.609467455621, (4,1,1) 46.67453891138, (3,3,0) 10.92853891670, octet functionals 184.3305612334 / −238.5587523317, singlet −5118.538314405, octet mixing −1.6943666102 / 1.0955642150 (vs chain 5·10⁻¹⁹). **Cost:** 10 products = 20 Bk₂; 33 distinct projective ladder vectors → 66 brackets; 16 Q-values (Q², Q³ on 4 points per ladder); vs det₆ 36 bracket entries + one 6×6. **[3,1,0] (cell 68):** full pivot scan gives best 45 (12,12,9,12), a 400-restart greedy elimination 48 (12,12,9,15); lower bound 24 not attained; minimality unproved. Gotcha: a comment containing `m*)` closes early and the cell silently returns Null.
- **Write-up DONE (read first): `Two-site-rotation-brackets-notes.{md,tex,pdf}`** (21 pp., sonnet; spec `Two-site-rotation-brackets-notes-spec.md`; every number cited as [A:n]/[B:n]/[C:n]/[D:n] = site-2 / two-site / Stage C adjoint / untwist notebook, absolute cell positions; App. C lists 13 source inconsistencies — notably: Stage D's result cells misattribute the 8/14 line to w₁ ≠ w₁′ (it is the global rotation), the Stage B follow-up points (2/3,−5/7), (3/2,1/4) are degenerate vacua (w₁ = z₁, z₃) though ranks/counts are unaffected, B and D cite code-cell ordinals while A and C cite absolute positions, "generic" in B means both a generic pair (23/54) and a generic family member (22/54)). Untested readings flagged there: the cc₁(1−t) vs cc₁(−t) "R shift", the triangular support rule, CG-completeness of the type-pair loci (α₁ = 0, β₂ = 0: 40 for [3,1,0]) and of β₂ = 0 with site 1 = Ωc[w] (16/48).
- **Paul (2026-10-11): "try everything" on the [3,1,0] cost question and κ(s,t). Four jobs, sequential on K1 (a second kernel K2 faulted on creation — licence), specs written, dispatched in this order:**
  - **Job A** `L2-two-vacua-detform-spec.md` → `L2-two-vacua-detform.wb`: is ⟨Ω_w⊗Ω_{w′}|Ψ⟩ (or the (∞,∞) family) a single det6L-type 6×6 bracket determinant with FREE per-ladder bra data (16/20 unknowns; exact polynomial-identity diagnostic + numeric fits with restarts; broader row/column ansätze).
  - **Job B** `L2-bracket-count-growth-spec.md` → `L2-bracket-count-growth.wb`: Σr/ΣP for generic, vacuum, Ωc[w]⊗Ωc[w′], (∞,∞) over [2,1,0], [3,1,0], [3,2,0], [4,1,0], [4,2,0] + polynomial fits in (λ₁,λ₂) + reach of (∞,∞) for the new reps.
  - **Job C** `L2-210-principal-torus-spec.md` → `L2-210-principal-torus.wb`: main-session observation Φ(a,a′) ∝ Ω₀·(ρ(D_a)⊗ρ(D_a))·(1⊗ρ(D_s)), Ω₀ = x₀(ρ(g₀)⊗ρ(g₀)), g₀ = exp(E₁₂)exp(E₂₃), D_a = diag(1,a,a²), s = a′/a — so κ(s,t) is a site-2 principal-torus insertion on a globally rotated vacuum; verify, grading of the SoV basis, localisation test (D₁–D₃ and rational (2,2)) of the s-components, vacuum-Q test of the kernels.
  - **Job D** `L2-detrank-locus-spec.md` → `L2-detrank-locus.wb`: the determinantal locus {Φ : rank A^{(k)}(Φ) ≤ r_vac} (adjoint block (1,1) ≤ 2; [3,1,0] block (0,0) ≤ 8 and φ-rank ≤ 12): tangent dimension at the vacuum (exact), Newton exploration, product-covector restriction, which points read more than the top irrep.
  **STOPPED by Paul (2026-10-11, ~12:00) before any result:** Job A was killed ~10 min into its setup evaluation (no notebook file was written; only the four spec files exist); Jobs B–D were never dispatched. The specs are ready to dispatch if Paul wants them later (sequentially on K1; each is self-contained). The kernel evaluation was cancelled.
- **Consequence for the plan:** "two-site count = product of one-site counts" is dead already (5 vs 17 slices). Stage B decides whether Ωc[w]⊗Ωc[w′] (or the R-dressed candidates) has det₆-like block ranks. Stages C/D follow the Stage B recommendation.

## >>> NEWEST (2026-10-10, late night): site-1 rotation brackets generalised to [3,1,0] (and [3,2,0]); closed-form measures; TeX notes <<<
**Read first:** `new-conventions/Site1-rotation-brackets-notes.tex/.pdf` (15 pp., sonnet write-up of the whole bracket thread, evidence tags, cell maps, source ambiguities in App. C). **Notebook:** `L2-site1-rotation-brackets-general.wb` (50 cells, rep-generic, Λ ∈ {210, 310, 320}; **no native outputs**; cells 30–41 = closed-form measures + t-dependence, cell 50 = full summary; Task 7 driver cells 43–48 wipe Global` — run last; **cells 30–50 were inserted by rewriting the .wb JSON with Python because VS Code did not have the file open — open it once in VS Code to confirm it loads**). Spec: `L2-site1-rotation-brackets-general-spec.md`. Not committed.
- **Label rule [exact, all reps]:** SoV labels (a|bc), λ₂ ≤ a ≤ λ₁, 0 ≤ c ≤ λ₂, c ≤ b ≤ a; roots θ_s + {−1−a, −b, −c}h; wavefunction monomial M_{(a|bc)} = det₂[v(θ−(1+a)h), v(θ−bh)]·v_k(θ−ch) (rows {2,3}, k = 3 and two other conventions); site-1 covectors ⊗ e₁ᵀ touch only the site-2 labels (λ₁|λ₁ c₂) with factors λ₂^{(c₂)} = Q^k(θ₂−c₂h)/Q^k(θ₂−λ₂h).
- **Bracket counts:** generic rotation 5 / 7 / 11 for [2,1,0]/[3,1,0]/[3,2,0] (dim V 8/15/15, |Λ| 3/4/5) — automatic for ANY covector (column 2-forms have rank ≤ 4); reduced 4 / 5 / 9 on: the adjoint quadric; the [3,1,0] 2-dim family **β + αγ = α²** (first row of g = (1,α,α²); t-independent) or the curve {243β = 280γ², 351α = 140γ}; the [3,2,0] curve {12493β = 3010γ², 403α = 140γ}. The exact (R2δ) analogue (all vacuum-group columns decomposable) is EMPTY for [3,1,0] and [3,2,0]: "fewest brackets" is the right generalisation. All loci CG-complete at finite twist (all weights, rank 2 on double weights).
- **Closed-form measures [exact]:** on the [3,1,0] locus every F, G entry is a Laurent monomial in α × at most one linear factor over (α−γ); c = 1 columns depend on α only; coefficients are short rational functions of t with poles at t ∈ {0,1,2,3} (cell 41). Natural (wgt·ρ^m·(p+qm)) measures: only isolated vectors on ≤ 3-point supports. Adjoint: c(t) = (t−3)/(2(t−2)) = (3/4)κ_singlet (structural? open).
- **u-structure of the measures [exact, cells 42–54, notebook now 63 cells, JSON-inserted again, backup `backup_general.wb` in the session scratchpad]:** gauge-invariant kernels A^{(c₂,c)} = α^{m+m′+9−c₂−c}·a(s,t), s = γ/α, a affine in s. Explicit F(u), G(u) exist for BOTH λ₂-group columns: F = h⁴α^{(u−θ₁)/h}/(P₃′(u)(u−θ₁)(u−θ₂+h)), G = h³α^{(u−θ₁)/h}/(P₃′(u)(u−θ₂+h)) (lone Q³(θ₁−h)); F = −h⁴α^{(u−θ₁)/h}ω(u)/(P₄′(u)(u−θ₂+h)), G = F·(u−θ₁)/h, ω = 1 except ω(θ₁−h) = 3(s−1)/(3s−2) (lone Q³(θ₁)) — SoV measure × bra twist α^{u/h} × cross-ladder pole × linear polynomial; verified 10⁻⁴⁶ on 225 states. Vacuum-group columns: lone Q³(θ₁−h) = (t−3)α^{…}/((t+m)(t+m′))·v∧(tx+y) with integer v,x,y (v = truncated ladder-1 measure, partner unidentified; Möbius test proves no measure×linear pair); lone Q³(θ₁) rank 4, no class fits. Adjoint: same law with offset 7; λ₂-column (0,1) = δ₋₃∧δ₋₂.
- **KEY IDENTIFICATION (main session, by hand, to verify first next time):** the [3,1,0] locus β+αγ = α² IS the one-site vacuum family: Ω_{(w₁,w₂)} (companion frame, `Ωc`) = e₁ᵀρ(g) with α = 1/w₁, β = −1/(w₁w₂), γ = (w₁+w₂)/(w₁w₂) ⇒ β+αγ = α² identically. So the 5-bracket family is ⟨Ω_w ⊗ x₀^{(1)}| and α^{u/h} = w₁^{−u/h} is the bra twist. (The adjoint quadric, cc ≠ 1, is NOT the vacuum family.)
- **PLAN (Paul: implement another day):** `new-conventions/Two-site-rotation-brackets-plan.md` — two different rotations on both sites / different site-2 reference, headline family ⟨Ω_w⊗Ω_{w′}| (bracket count via the two-site K tensor and (c,c′)-blocks, det₆ as regression, explicit measures, then the untwist with w, w′ → 1 independently, which is clean by the commutator formula and needs no ε-scaling search). Stages A–D with order and risks.
- **NEXT:** untwist of the one-site families — spec `L2-site1-rotation-untwist-spec.md` written (adjoint quadric regression, then [3,1,0]; decisive step: ε-scaling of g on the locus vs the Gauss factors L₀D₀U₀ of P_phys⁻¹(z), 3×3 symbolic), not yet dispatched (Paul to confirm). Then: meaning of the quadric constants; site-2/two-site versions.

## >>> NEWEST (2026-10-10, night): gl(3) adj⊗adj — POSITIVE: the generic site-1 rotation on the SoV vacuum is a four-bracket-determinant Q-expression on a quadric of rotations <<<
**Write-up (read first):** `new-conventions/L2-210-site1-rotation-brackets.md`. **Spec:** `L2-210-site1-rotation-brackets-spec.md`. **Notebook:** `L2-210-site1-rotation-brackets.wb` (36 cells, saved, not committed, **no native outputs** — cells run in K1, results in markdown cells; Task 6 driver cell wipes Global`, run it last). Paul's idea: forget Ω_w, start from x₀ = e₁ᵀ⊗e₁ᵀ (companion frame), generic rotation g on site 1, pull out λ₂ = Q^k(θ₂)/Q^k(θ₂−h), ask for which g the θ₁ part is pure ladder-1 brackets ⟨Q^jF⟩₁ with explicit measures (no det required).
- **Result [exact; 64-state checks 10⁻⁵³]:** with g = exp(αE₁₂)exp(βE₁₃)exp(γE₂₃), on the quadric **2(t−2)(β+αγ) = (t−3)α²** (t = θ₁₂/h, h-independent, z-independent; = 2(t−2)b₁b₃ = (t−3)b₂² for the first row b of g) — a 2-parameter family of genuine rotations — ⟨Φ(g)|Ψ⟩/⟨x₀|Ψ⟩·Ψ₁ = Σ_{r∈{−1,0}}(Bk₂[F⁰_r,G⁰_r] + λ₂Bk₂[F¹_r,G¹_r])·Q³(θ₁−rh), Bk₂[F,G] = ⟨Q³F⟩⟨Q²G⟩ − ⟨Q²F⟩⟨Q³G⟩, eight explicit 4-point measures. The λ₂¹ part has this form for every rotation; the quadric is the condition on the λ₂⁰ part (Plücker decomposability of the r = 0 column of the coefficient tensor). More restrictive forms (single det₂×bracket; single L=1-type 3×3 det; λ₂ absorbed into the measures) contain no genuine rotation (β = 0 or α = β = 0 only).
- **Mechanism:** ψ_x·Ψ₁ = det₂[v(n′),v(n)]·v_k(r) monomials for the eight site-1 SoV triples; N⁰ = Σφ⁰_xM_x (8 slots), N¹ (4 slots, always decomposable); Pf(C⁰_{·,0}) ∝ β²(16β+16αγ−43α²) at t = 62/35.
- **CG reach at finite twist [exact]:** undone covector e₁ᵀρ(gL₀D₀) (P_phys⁻¹ = L₀D₀U₀, site 2 = e₁ᵀρ(U₀)) reaches all 7 weights with rank 2 on (1,1,1) on the quadric (β = 0: 5 weights, rank 1). **Cost:** Q², Q³ on the four ladder-1 points + Q³ at two points + λ₂; four 2×2 bracket dets — fixed by ladder length and rank, not dim V.
- **Natural (det₆-type) measures:** F, G can each be w^{−x/h}·(p+qx)·wgt₁ with different w's per column; no common w (exact).
- **NEXT (ask Paul):** (1) untwist the quadric family (z = 1+εφ; choose (α,β,γ)(ε) on the quadric so the undone e₁ᵀρ(gL₀D₀) is a generic unipotent at ε = 0; then K/S, CG at leading order as in `L2-210-untwist-CG.wb`); (2) meaning of (t−3)/(2(t−2)) vs y₀ = (t−1)/(t−2) and κ_singlet; (3) site-2 analogue / both sites rotated; (4) rerun in VS Code.

## >>> NEWEST (2026-10-10, evening): gl(3) adj⊗adj — reachability map of Q-number product covectors on a single vacuum; NO gl(4)-style clean family below the vacuous level <<<
**Write-up (read first):** `new-conventions/L2-210-localising-span.md`. **Spec:** `L2-210-localising-span-spec.md` (§0–§4 + Addenda A–C = four sonnet rounds). **Notebook:** `L2-210-localising-span.wb` (63 cells, saved, not committed, **no native outputs** — cells run in K1, results in markdown cells; rerun in VS Code ≈ 1 min). Paul's ask: a gl(3) analogue of the gl(4) rect result — Q-number need not be a determinant, but its cost must scale with boxes, not dim V; every CG datum at leading order.
- **Method:** base ⟨Ω_w⊗Ω_w| (overlap det₆); at L = 2 the algebra has six non-scalar coefficient operators (t₁: 2, t₂: 4; t₃ scalar), spans D₁ (7), D₂ (24, four exact quadratic relations), D₃ (48), D₄ = everything (vacuous). Exact scans at two twist points for polynomial membership and rational forms P/Q at the non-vacuous levels (2,2) and (3,1).
- **Reachable as products v⊗Ω₁:** polynomial: only Ω₁, Ω₁T₁₁^{(1)}, Ω₁T₂₁^{(1)} at θ₁−2h (ℓ₁ = ℓ₂ = 3 through D₃); rational (2,2): + the common-w₁ line Ω_{(w₁,y)} (either site, products of two); rational (3,1): nothing new. NOT reachable at any level: polynomial torus at a free point, resolvent, ρ(z), different-w₁ vacuum, every torus × y-deformed vacuum (incl. the y₁-freed product point, which by the exact commutator formula n′₁₃−n′₁₂n′₂₃ = (δ₃/δ₁−1)(w₁−z₁)z₂₃/((w₁−z₃)(y₁−z₃)) WOULD untwist cleanly).
- **Exact facts:** T₂₁(θ₁−2h) on the vacuum = −3h³χ₃²·ρ₁(1/z) (inverse ket-twist torus on site 1; gl(3) shadow of gl(4)'s t₂(θ₁−h) ∝ ρ₁(z)⁻¹). T₁₁^{(1)} is NOT a torus (non-scalar on (1,1,1); eigenvalues (3/4)χ₃/z_m on 1-dim weights). Only tori in the localising span: d = 1, 1/z; T₁₀, T₂₀ add no torus line. b = r₂+r₁r₃ = (w₁−z₁)/(w₁−z₃) is y-independent, so the y-family moves only E₂₃ and reads [4,1,1] alone (clean: K/S = O(ε), Richardson 46.67453891 = chain). ω₂ = Ω₁ρ(exp(γ(E₂₃−aE₁₃))), a = r₁, γ = −r₃(w₂) (not pure E₂₃). Enlarged 8-operator family: site-2 column space 5-dim, three directions off the common-w₁ line, two not rotations.
- **Conclusion:** the gl(3) adjoint has one free rotation-preserving deformation (y) and no free-point torus, so the paper's two-ingredient mechanism has no low-degree realisation here. Suggested next (ask Paul): gl(3) [S,S,0] rectangular (abelian nilradical, the true gl(3) analogue of the paper) and gl(4) [S,S,0,0] S = 2 with the same (2,2) scan; or a two-vacuum base ⟨Ω_w⊗Ω_w|P + ⟨Ω_{w′}⊗Ω_{w′}|P′ (D₁×D₂-type only is non-vacuous).
- Memory saved: `project-qnumber-cost-constraint` (Paul's acceptance criterion).

## >>> NEWEST (2026-10-10, latest): gl(4) rectangular [S,S,0,0] — the 2603.26617 two-rotation family in our conventions, and WHY its untwisting is clean at leading order while the gl(3) family is not <<<
**Write-up (read first):** `new-conventions/L2-gl4-rect-untwist.md`. **Spec:** `L2-gl4-rect-untwist-spec.md` (+ Addendum A — the original §2 family was wrong, see below). **Notebooks:** `L2-gl4-rect-untwist.wb` (S = 1, 34 code cells) and `L2-gl4-rect-untwist-S2.wb` (S = 2, 29 cells), sonnet-built, saved, not committed, **no native outputs** (cells run in K1, results in markdown cells under each code cell; rerun in VS Code). Paul's question: why does untwisting ⟨W_{ω,κ}⊗W_ω| generate everything at leading order in the paper, but the gl(3) det₆×det₂ family needs 2–3 orders?
- **Answer:** (1) abelian nilradical: rotations of the HW covector are e₁ᵀexp(X), X ∈ span{E₁₃,E₁₄,E₂₃,E₂₄}; the undo is entrywise X′_ij = x_ij(δ_j/δ_i−1), no commutator direction (gl(3): b′−a′c′ = a′c′(z₂−w₂)/(y−z₂), exact). (2) Two independent localising insertions at θ₁−Sh (our conventions; the paper's θ₁ is our θ₁−h): the κ-insertion 1+(κ/h)t₁+(κ²/2h²)t₂ = complementary-pair torus = **resolvent** ρ(diag 1/(w₂′−z)) (one-site identities t₁ = −h(z_c′+z_d′), t₂ = 2h²z_c′z_d′ exact at S=1), and the ω-deformation ⟨W_ω| = Ω_zw(1−ωρ(diag z))/(1−ωz₁z₂) (a rotation, paper's (F5)), which moves the vacuum's torus invariant (cross-ratio); nothing ties κ to ω, both untwist with z. In gl(3) the product locus pins y = w₂(t−1)/(t−2), so y cannot untwist and the E₁₃ direction loses an order.
- **Sharp facts:** companion vacuum block x_ij = (−1)^i(z_{i′}−z_j)/(z₁−z₂)·Π_{k=1,2}(w_k−z_i)/(w_k−z_j), cross-ratio = z cross-ratio, w₁,w₂ torus, w₃,w₄ drop; the κ-torus on it (no ω) has det X′ ≡ 0 — blind to the singlet at ANY twist. Holding ω fixed while z→1 reproduces the gl(3) disease (singlet slice O(ε²), K/S constant 0.027 at S=1; at S=2 three irreps contaminated, slices at ε², ε⁴).
- **Joint regime (z,w,ω,w₂′ = 1+ε·dir):** every weight O(1), full rank (S=2 dim-2 weight rank 2), K/S = O(ε) for all irreps (S=1: ≤ 6·10⁻⁷; S=2: ≤ 1.2·10⁻⁶ at ε=10⁻⁶), CG from the overlap alone with O(ε) error and Richardson 10⁻¹². S=1: three-term relation exact, c = κ²/Π(1−κz)·(1, κ/(hf₁), κ²/(2h²f₂)), f₁ = 27/70, f₂ = −54/1225; ⟨W⊗W|Ψ⟩/⟨Ω⊗Ω|Ψ⟩ = [1−ω(E₁+E₂)+ω²χ]/(1−ωz₁z₂)² with explicit τ₂-combinations (36 states, 190 digits) ⇒ the family is a Q-number. **S=2 open:** t₂(θ₁−2h) is not a torus on the 20-rep; ρ(Δ⁻¹) not in the degree-≤2 monomial span at θ₁−2h.
- **NEXT (ask Paul):** Q-number form at S=2 (pole-set products θ₁−kh?); the explicit (G7) Q-det in our conventions; whether gl(3) admits an ω-like second deformation; rerun notebooks in VS Code.

## (2026-10-10, earlier): untwisting the two-rotation det₆ × det₂ formula — ALL six HW CG data of adj⊗adj from the Q-formula alone; structural O(ε) constants κ_singlet(t) and Δ_c(t, τ_c)
**FULL NOTES (read first):** `new-conventions/Untwist-two-rotation-notes.{md,tex,pdf}` (29 pp., sonnet write-up of the whole thread, every result with notebook+cell, evidence tags, 18 source inconsistencies in App. C — note esp. C4: an unresolved mismatch in `L2-310-untwist-CG.wb` cell 61's t = 62/35 sanity line for the [5,2,1] Δ_c, and the two meanings of Δ_c (adjoint: slope + (t−2)R; [3,1,0] [5,2,1]: slope − p₁′ − q₁′R)). Per-topic write-ups: `L2-210-untwist-CG.md` (adjoint, formulas, cell map), **Spec:** `L2-210-untwist-CG-spec.md`. **Notebook:** `L2-210-untwist-CG.wb` (45 cells, saved, **no stored outputs** — VS Code could not execute its cells, the sonnet agent ran every cell's source in K1; rerun in VS Code, a few minutes). Not committed (ask Paul). Also untracked: `img/L2-210-untwist-CG/`; `.vscode/settings.json` has a 12-line change not made by this session.
- **Regime (decided by the main session, sympy-checked, then verified in the notebook):** untwist everything together, z = 1+εφ, w₁ = 1+εω₁, w₂ = 1+εω₂, y = w₂(θ₁₂−h)/(θ₁₂−2h). Then the undone covector has O(1) components at (2,1,0), (1,2,0), (2,0,1), (1,0,2) and one direction (u+u′) of (1,1,1), and K/S = O(ε) for [4,1,1], [3,3,0], both octets (max |K/S| ≤ 3·10⁻⁶ at ε = 10⁻⁶). The E₁₃ direction degenerates by one order: b′/(a′c′) = 1 + ε(t−2)(φ₂−ω₂) exactly (t = θ₁₂/h), so the octet mixing and the singlet sit at O(ε) with the first contamination. At fixed y (old table ε⁰…ε⁴) nothing is clean — don't use it.
- **det₆ × det₂ identity at small ε:** state-independent to 1.6·10⁻¹¹⁶ on all 21 states of the HW sectors (200 digits). Pure Q-side number Qnum(B) = R_Q(B)/R_Q(top)·⟨v₁⊗v₂|e₁e₁⟩/⟨x₀|e₁e₁⟩.
- **Structural O(ε) constants (fitted at h = 1 over 8 resp. 22 values of t, validated on held-out t's and at P2, P3):**
  - singlet: K/S = −t/(3(t−2)), κ = 2(t−3)/(3(t−2)); CG_singlet = Qnum/(κ·φ₍₀,₁,₂₎(ε,ω)). Reason: ψ₁ comes from h²Σ_j E⁽¹⁾_kj E⁽²⁾_jk (θ-independent) while the signal carries (t−2).
  - octets: after top-killing (subtract det₆ of the (w₁,y)-vacuum), F = Qnum^tk/(δ^Λa′₀c′₀) = F₀ + εF₁; F₁/F₀ is a straight line in ω₂ with φ-independent slope; Δ_c := slope_c + (t−2)R_c = (1−t²) + t τ_c/(3(t+3)), τ_c = untwisted t₂(θ₁+h)/h⁴ eigenvalue of copy c (roots of τ² − 3(2t−1)(t+2)τ + 9t(t−1)(t+1)(t+3) = 0; Q-number via t₂(θ₁+h) = χ₃ν₁(θ₁+h)ν₂(θ₁)Q^j(θ₁−h)/Q^j(θ₁)). Mixing R_c = (u·C_c)/((u+u′)·C_c) = (Δ_c − slope_c^Q)/(t−2), u = e₁ᵀE₁₃, u′ = e₁ᵀE₂₃E₁₂.
- **FINAL (cells 41–45):** all six HW states from the Q-formula alone at P1, P2 and P3 = (h=1/2, θ=(3/7,−1/4)) (out of every fit set): errors O(ε) (≤ 3.6·10⁻⁶ at ε = 10⁻⁶), Richardson limits Q vs chain ≤ 6·10⁻¹⁹, octet mixing R^Q vs R^chain ≤ 8·10⁻¹¹.
- **Bug (cost one wrong law):** `Do[DownValues[f] = …, {f, {…}}]` purges nothing (HoldAll); stale P1 transfer matrices gave K/S = −(62/105)/(t−2) with 62/105 = t_P1/3. Fixed in cell 36 (literal symbols). Also: the spec's w₃ = 1+ε·2/9 collided with w₂ at ω₂ = 2/9 (use 5/13); `ClearAll` needs `Quiet` (Protected `ClearGlobals`).
- **[3,1,0] point A through the same pipeline (2026-10-10, later; `L2-310-untwist-CG.wb`, 90 cells after round 2, spec `L2-310-untwist-CG-spec.md`; no native outputs — the agent stored printed results as markdown cells under each code cell; VS Code never opened the file).** Finite-twist regression of the closed-form dressing: undo exact on 10 HW kets, det identity 1.3·10⁻¹⁷⁹ on 225 states; K13 = −189/704 (the notes' +189/704 is a sign typo). **In the joint regime the dressing cancels E₁₂² exactly:** φ_A(0)/δ^Λ = e₁ᵀ(1 + (a′₀/2)E₁₂)e^{c′₀E₂₃}, γ = O(ε). Per-weight orders: ε⁰ at (3,1,0),(3,0,1),(2,2,0),(2,0,2) and one direction of (2,1,1); ε¹ at (1,3,0) and one direction of (1,2,1),(1,1,2); ε² at (0,2,2),(1,0,3),(0,3,1) and the second directions of (1,2,1),(1,1,2); ε³ at (0,1,3). So **three orders** of the det expansion are needed for [3,1,0] (two for the adjoint). Clean at leading order (κ−1 ∝ ε): [6,1,1], [5,3,0], [5,2,1]×2 (u′-functional); contaminated with ε-constant κ: [4,4,0] κ = 97/70 = (t+1)/2 at P1 (the Clean-CG number), [4,3,1] copies 1.41057 / 0.98671, [4,2,2] κ = 1+108/295, [3,3,2] κ = 4.65080 — all ω- and φ-independent (spreads ≤ 10⁻⁵). Q side: identity to 10⁻⁸² (limited by 104-digit R_Q), leading-order extraction O(ε) errors, Richardson Q = chain to 10⁻¹⁹ for the five clean branches. [5,2,1] pair O(ε) layer: F₁/F₀ linear in ω₂, slope φ-independent, Δ_c = −0.341940229, 0.765796900 (u = e₁ᵀE₁₃ = (−1,8), u′ = e₁ᵀE₁₂E₂₃ = (2,8) on the (2,1,1) GT indices, R_c = u·C/u′·C = −0.892189, 1.204070). **Not done:** closed forms in t of these constants (t-scan + τ_c ansatz as for the adjoint), the ε² layer (second directions of [4,3,1]×2, [4,2,2]), the final pure-Q extraction.
- **Adjoint constants DERIVED (2026-10-10, `Adjoint-untwist-perturbative.wb`, 22 cells, ~40 exact checks; spec `Adjoint-untwist-perturbative-spec.md`; no native outputs, VS Code could not open the file).** First-order perturbation theory Ψ = h + εψ₁ (any family member; checked with two), symbolic in t, φ, ω: (i) singlet K₁/S₁ = −t/(3(t−2)) exactly; 100% of K₁ comes from the two octets (the twist perturbation is an adjoint tensor operator: singlet → octets only), per copy K₁,c/S₁ = A + Bτ_c with A = (3t⁴+2t³−9t²−11t+6)/(3(t−2)(5t²+4)), B = −(t²−2)/(3(t−2)(5t²+4)); (ii) octets: t₂(θ₁+h) on HW₍₃,₂,₁₎ is [[3(2t²+5t−1)/2, 3(t−1)/2],[3(t−3)/2, 3(t−1)(2t+3)/2]], discriminant 9(5t²+4), parametrized rationally by t = 4m/(5−m²), √(5t²+4) = 2+mt; Δ_c = (1−t²) + tτ_c/(3(t+3)) exactly on the whole curve, ½ from each decuplet; the 27 and the other octet give ω₂-independent (intercept-only) pieces. Bonus: clean branches' O(ε) contamination, top-killed, is ω-independent: [4,1,1]: K₁^tk/S₀ = (t−1)(φ₂−φ₃)/4; [3,3,0]: (t−1)(φ₁−φ₂)/4. Numerical cross-check at t = 8 (τ = 252, 198 rational) to 3·10⁻⁸ (τ = 252) and 2.4·10⁻⁶ (τ = 198), at ε = 10⁻⁸.
- **[3,1,0] round 2 DONE (`L2-310-untwist-CG.wb`, 90 cells; write-up `L2-310-untwist-CG.md`).** All 10 HW branches from the Q-formula alone at P1 and at t = 19/14 (errors O(ε), Richardson ≤ 1.6·10⁻¹⁸). Closed forms (t-scan, held-out ≤ 10⁻¹⁵): κ[4,4,0] = (t+1)/2, κ[4,2,2] = (11t−31)/(5(3t−7)), κ[3,3,2] = (t+1)(19t−49)/(40(t−2)), κ[4,3,1]_c = A+Bτ_c (D = t³+3t²−17t+5); [5,2,1] mixing R_c = 3/2 − 2t + τ_c/(2(t+4)); exact untwisted t₂(θ₁+h) on HW_[4,3,1] (tr 8t²+16t−8, det 16t(t−1)(t+1)(t+4)) and HW_[5,2,1] (tr 8t²+19t−8, det 4(t+4)(t−1)(t+1)(4t+3)); ε² layer: [4,2,2] R = −2(t−3)/(3(3t−7)), [4,3,1] R_c = A_R+B_Rτ_c exact, their structural slope constants fitted (the [4,3,1] α(t) only numerically, 10⁻⁹ → second directions of [4,3,1] to 10⁻⁸ only). Open: α(t); perturbative proofs for [3,1,0]; [4,1,0].
- **Q-dets untwisting, stage C1 DONE (`L2-210-untwist-Qdets.wb`, 45 cells; write-up `L2-210-untwist-Qdets.md`; spec `L2-210-untwist-Qdets-spec.md`).** det₂^{(1)}/minor₁ fully untwisted: = (343000/1161)·F(y)·(1,τ_λ) at small ε, τ = τ⁰ + ετ¹ (τ¹ = Σφ_k a_k, first-order shifts), τ⁰ = untwisted quantum-Wronskian formula with the staggered-degree polynomial kernel (10⁻¹⁹⁷). det₆/Ψ_Q NOT untwisted: exponents 10,10,10,9,9,8; D6_B/D6_top are exact monomials in the directions, and for [4,1,1], [3,3,0] exact in ε with no Q input ((w₁−z₂)(w₂−z₂)/((w₁−z₃)(w₂−z₃)) and (w₁−z₁)/(w₁−z₂)); the old Laplace split (bare × periodic) holds only for the top (mixed Laplace terms same order; untwisted periodic minor vanishes on non-top states); H columns rational (c1+c5, c3−c5, c4−5h·c5−c6 with "5h" untested), O not. Qnum reconstructed to O(ε²) with ρ_B numerical for octets/singlet. **Open: the det₆ column basis for non-top branches; the common non-monomial factor; τ¹ in Q-language.**
- **All notebooks of 2026-10-10 lack native outputs** (VS Code did not execute cells of new files for any agent): `L2-210-untwist-CG.wb`, `L2-310-untwist-CG.wb`, `Adjoint-untwist-perturbative.wb`, `L2-210-untwist-Qdets.wb`. Rerun each in VS Code to store outputs. Nothing committed.
- **NEXT (ask Paul):** (1) analytic untwisting of det₆ × det₂ in the regime (leading coefficient through untwisted staggered-degree Q's; old technology in `old-conventions/old/XXX/Experiments/SU3_L2_Untwisting.{md,tex}`); (2) the [3,1,0] k = 2 point (`L2-310-two-rotations.wb` point A) through the same pipeline, then [4,1,0] k = 2; (3) derive κ_singlet(t), Δ_c(t,τ_c) by first-order perturbation theory (as `Adjoint-clean-analytic.wb` did for β(t)); (4) rerun the notebook in VS Code to store outputs.

## >>> NEWEST (2026-10-09, late): undoing the site-2 rotation of the two-rotation family — adjoint complete (8/8, singlet); [3,1,0] and [3,2,0]: all HW states reached, 9/10 determined by the det-form point, [3,1,0]'s point A gives 10/10 <<<
**FULL NOTES (read first):** `new-conventions/Undo-site2-rotation-notes.{md,tex,pdf}` (23 pp., sonnet write-up of the whole thread, every result with notebook+cell, conjectures tagged). **Spec + running log:** `new-conventions/site1-Qform/7-undo-site2-rotation-span.md` §1–§16. **Notebooks:** `L2-rotation-sov.wb` §12 (cells 113–123), `L2-310-two-rotations.wb` cells 24–70, `L2-320-two-rotations.wb` (41 cells), `L2-410-two-rotations.wb` (26 cells) — all saved, not committed; sonnet subagents from the spec. Paul's question: for ⟨v₁⊗v₂| with v₂ = ⟨0|R, does ⟨v₁|R⁻¹ projected to the site-1 weight spaces (all coefficients of the free parameters) give a complete basis? The 2026-09-23 families stopped at μ₃ ≤ λ₂ (singlet missed).
- **Undo is exact at finite twist:** in the diagonal frame v₂ = e₁ᵀρ(n), n ∈ N₊ (the vacuum Gauss factors r_k(z,(w₁,y))), and ρ(n)⊗ρ(n) fixes all HW kets, so ⟨v₁⊗v₂|HW⟩ = ⟨v₁ρ(n)⁻¹⊗e₁ᵀ|HW⟩. No untwisting needed for the count.
- **Mechanism (3×3, exact, rep-independent):** φ = v₁ρ(n)⁻¹ = δ^Λ e₁ᵀρ(n′), n′ = Δ⁻¹nΔn⁻¹ = [[1,a′,b′],[0,1,c′],[0,0,1]], Δ = diag(1/(ω−z_k)), ω = resolvent point w₂(y). **b′ = 0 ⟺ ω = y**, and b′/(a′c′) = (y−ω)/(y−z₂). The old families' p₃ relation is exactly b′ = 0 (Δ = diag z); the new family has ω = y(θ₁₂−2h)/(θ₁₂−h) ≠ y — the (t−1)/(t−2) mismatch between the det₆ twist and the ladder insertion is what opens the E₁₃ direction.
- **Adjoint:** per-weight ranks = dims (8/8), exact in y, unchanged under w₁/torus; control ω→y kills (0,1,2),(1,0,2) (old blind spot). HW count: 5 at fixed w₁ (one harmless cross-weight relation), **6 with the torus**; singlet reached. Untwisting leading orders ε⁰,ε¹,ε¹,ε²(rank 2),ε³,ε³,ε⁴. (0,1,2) component = 3b′(a′c′−b′)δ^Λ.
- **[3,1,0] (point-B family over six w₂, β/α re-solved per w₂):** rotation part alone spans 15/15; the actual dressed covector e₁ᵀ(α+βE₁₂)ρ(nΔn⁻¹) spans **12/15** — the fitted dressing cancels (0,1,3), (1,0,3) identically and one direction of (1,1,2) (consistent with the e₁ᵀE₁₂-orbit type). Still **10/10 HW states reached** with the torus ([4,2,2], [3,3,2] included; old 8/10); 9 fully determined, [4,2,2] only up to one combination of its two (1,1,2) components. Control (ω→y*, β=0): 6.
- **Kernel state at end:** K1 holds the [3,1,0] notebook's cells 2–7, 12, 24–33; `L2-rotation-sov.wb` live state wiped (rebuild cells 2–21 + the text of `oneSiteVac`/`w2Of` for §12; 63–111 need their own cells). Private kernels (K310, K2) cannot start in this environment — every notebook ends up on K1.
- **Follow-ups (Paul: "get cheap subagents to start doing these"; spec §9–§10):**
  - **[3,1,0] LW frame** (cells 34–42): exact mirror of the HW frame; [4,2,2] again rank 1 of 2 (at (2,1,1)); frames don't combine. No help.
  - **[3,1,0] point A** (same cells): site-1 covector (rotRank 6) undone exactly → **all 15 directions, all 10 HW states fully determined** (A alone or A∪B). Caveat: A's overlap is c_A·(1,τ_λ(θ₁−3h))⟨Ω|Ψ⟩ (T-polynomial × det₆, a Q-number via T–Q) but not a single ladder det₂.
  - **[3,2,0]** (`new-conventions/L2-320-two-rotations.wb`, 21 cells, saved, not committed): one isolated point only (point-B numbers; **no point A** since T₂₂, T₃₂ vanish), sites swapped (rotation on site 1 with y* = w₂(θ₁₂−h)/(θ₁₂−3h), dressed vacuum e₁ᵀ(1+bE₂₃)R on site 2, b fitted = (6125/5676)w₂/((w₂−z₃)(y*−z₃)) at P1). Undo exact; 12/15 directions; **10/10 HW reached, 9/10 determined, [4,4,2] one direction short**; control blind spot μ₁ ≤ 1 (dual of μ₃ ≥ 2).
- **A1/A2 (Paul's three-avenue plan, spec §11–§14):** **Point A of [3,1,0] IS a contour det₂:** the same ladder-1 det₂ as B but over the adjoint pole set {0,1,2} at y_A = w₂(t−h)/(t−2h), divided by the [3,1,0] minor W(3,4); its site-1 factor is a SECOND-order dressed rotation e₁ᵀ(1+βE₁₂+γ(E₁₃+½E₁₂E₂₃))ρ(n_A)ρ(Dy[w₂]) with closed-form β, γ (spec §13). **Explicit Q-dictionary:** τ_{(λ₁,λ₂)}(θ₁−3h) = C_λ·W(θ₁−(3−λ₁)h, θ₁−(4−λ₂)h)/W(θ₁−3h,θ₁−4h), C_λ explicit in h, t only; t₁, t₂ at any x by Cramer from the dual Baxter equation. So every "T-polynomial × det₆" overlap is explicit in Q's. [3,1,0]²: both points det₆ × det₂^{(P)}, A alone determines 10/10. **A2 negative:** [3,2,0] site-2 insertions give only the localizing span (19-op span, modular Gröbner) → [4,4,2] still one direction short. **A3 DONE — the pole-set chain (spec §14–§16):** product points of the θ₁−λ₁h recipe ↔ pole sets P_k = {0..k}, k = 2..λ₁, at y_k = w₂(t−h)/(t−kh), overlap det₆ × det₂^{(P_k)}[y_k]/W(λ₁,λ₁+1), site-1 dressing of E₁₃-degree λ₁−k. Exact on [3,1,0] (only {0..3}→B, {0,1,2}→A among all 11 subsets) and **[4,1,0]** (`new-conventions/L2-410-two-rotations.wb`, 26 cells, saved; three points k=2,3,4 at −63/8, −63/43, −21/26; 576 states end-to-end 6e-43; vacuum triple (θ−5h,θ−4h,θ−h), minor W(4,5); 14 HW states in 11 irreps). **The k=2 (adjoint pole set) point spans all 24 directions and determines all 11 irreps**; k=3: 21/24, 11/11; k=4 (=B): 16/24, 9/11. [3,2,0]: the k=2 curve collapses into the localizing span (T₂₂,T₃₂ ≡ 0 on ⟨Ω|), so only B exists there ([3,2,0] vacuum triple (θ−4h,θ−3h,θ−2h), minor W(3,4), end-to-end 1e-44). Rep-data pattern: triple (θ−(λ₁+1)h, θ−λ₁h, θ−λ₂h), minor W(λ₁,λ₁+1). NB the [4,1,0] notebook's cells ran through the kernel (controller not associated); outputs not stored — rerun in VS Code (~1 min).
- **Pattern so far:** the pure-rotation partner always spans everything (adjoint: the whole family is pure rotations → complete). Where the recipe produces a *dressed* vacuum ([3,1,0] B, [3,2,0]) the dressing removes exactly the extremal μ-directions and one direction in one dim-2 weight → one bottom irrep ([4,2,2] resp. [4,4,2]) undetermined. [3,1,0] has a second point (A) that fixes it; [3,2,0] does not.
- **Kernel at end of day:** K1 holds `L2-320-two-rotations.wb` (all cells live); `L2-rotation-sov.wb` and `L2-310-two-rotations.wb` live state wiped (files saved: 123 and 42 cells). Private kernels don't start on this machine.
- **NEXT (ask Paul):** (1) reps where the k=2 curve collapses ([3,2,0]-type, λ₂ ≥ 2): find the missing covectors (site-2 analogue was empty; try the dual/LW construction or products); (2) [4,2,0] (dim 27): which k survive? (3) derive the dressing polynomials (closed forms β, γ exist for [3,1,0]) and the general c_k(y) curve from the det₂ residue expansion with the τ-dictionary; (4) the Q-side reading of the CG coefficients: each y-coefficient of det₆·det₂^{(P_k)}[y] is one HW CG functional after the (R⊗R)Ψ→HW cleaning of `Clean-CG-combinations` — put the k=2 point through that pipeline for [3,1,0].

## >>> NEWEST (2026-10-09, evening): [3,1,0]² — the two-rotation recipe gives "rotated vacuum ⊗ rotated DRESSED vacuum", not two rotations <<<
**Notebook:** `new-conventions/L2-310-two-rotations.wb` (14 cells, exact rationals, sparse; saved, not committed). **Write-up:** `new-conventions/L2-310-two-rotations.md`. Paul's request: generalize the adjoint two-rotation result to [3,1,0].
- **Mechanism identical:** at θ₁−3h, (1,1),(2,2),(3,2) localize; all other shapes λ₁ ≤ 3 have rank 2 with ONE common extra site-2 direction ω₂ (a rotation, pencil inside the orbit); exact product locus in the 10-dim span = localizing span + two isolated points. Shapes λ₁ = 4 and the mixed span (+θ₁−2h shapes; 9-dim site-2 span, 29 s Solve) add nothing.
- **Point A** = the adjoint's c* verbatim on the column-shifted shapes (1,1),(2,1),(2,2),(3,1),(3,2); site 2 = rotated vacuum with the adjoint's y = w₂(θ₁₂−h)/(θ₁₂−2h); site 1 rotRank 6, unidentified.
- **Point B** (all shapes): site 2 = rotated vacuum with **y* = w₂(θ₁₂−h)/(θ₁₂−3h)** (three numeric points, exact); **site 1 = e₁ᵀ(α+βE₁₂)·R(z,(w₁,y*))·ρ((w₂−Gᵀ)⁻¹)** — the adjoint's v₁ with the vacuum covector replaced by vacuum + first descendant (β/α = 392/2335 at P1; 294/2255 at h=1/3; 5400/337687 at the third point; no closed form). Not a rotation (rotRank 5 = e₁ᵀE₁₂ orbit type; exact). For the adjoint the same dressing IS a rotation (E₁₂² kills e₁ᵀ), which explains the adjoint result.
- **New tool:** exact linear rotation test for any rep: rank{φ, φE_ij} = 4 ⇔ φ ∈ GL(3)·e₁ᵀ (generic 9; e₁ᵀE₁₂-type 5; e₁ᵀ(E₁₂+E₂₃)-type 6; e₁ᵀE₁₃-type 7).
- **Q-side DONE (cells 15–23, Opus subagent from `L2-310-Q-side-spec.md`, reviewed):** 225 states, all Q's (deg q^j = μ_j, deg q_i = 6−μ_i); vacuum det₆ with poles n = 0..3 and Ψ_Q = Π_s[Q²(θ_s−3h)Q³(θ_s−4h) − Q³(θ_s−3h)Q²(θ_s−4h)]Q³(θ_s−h) (5e-45); ladder-1 det₂ with minor₁ at (θ₁−3h, θ₁−4h) lies in span{1, τ_λ(θ₁−3h)} with an exact box-graded curve c_λ(y) = a_λ y^{6−|λ|}; it hits c_B only at y = −63/43 = y* (exact); END-TO-END: ⟨v₁⊗v₂|Ψ⟩/⟨x₀|Ψ⟩ ∝ det₆(w₁,w₂) × det₂^{(1)}(y*) to 5e-45 (wrong y: 77, 322). Point A is not on the curve.
- **Not done:** Q-side at the other numeric points; closed forms of a_λ and β/α; point A's site 1.
- **Side effect:** the new notebook landed on the shared kernel K1 and its `ClearAll` wiped the live state of `L2-rotation-sov.wb` (file intact; rebuild cells 2–21, 24–27, 60–68, 86–106 if needed). The `y*)`-in-a-comment gotcha struck twice; write y-star in comments.
- **NEXT (ask Paul):** is "rotated dressed vacuum ⊗ rotated vacuum" useful (its overlap is a T-polynomial × ⟨Ω|Ψ⟩, hence a Q-expression)? If so: Q-side for [3,1,0]; [3,2,0]/[4,2,0] with the same notebook (change Λ; localization data in the L2-*-localization.md files); closed form of β/α.

## >>> (2026-10-09, later): beyond the vacuum — the two-rotation family is "any rotation on site 2, one fixed resolvent on site 1" (3-dim); no κ-family; no iteration <<<
**Notebook:** `new-conventions/L2-rotation-sov.wb` §11 (cells 107–112, saved, not committed). **Write-up:** `new-conventions/Two-rotations-det-overlap.md` §5.8 (+ summary item 4, §6 updated); short form appended to `site1-Qform/6-two-rotations-2603-trick.md`. Paul's question: can the §10 construction be generalized to states beyond the vacuum / more general rotations on each site?
- **Torus = character (exact).** ρ(D)⊗ρ(D) (ket-twist Cartan torus) multiplies Ψ_B by Π d_i^{μ_i}; for a resolvent this is the su(2)-like det₄. Every product state with a det overlap comes with its torus orbit.
- **w₁ is a torus direction, y a torus invariant** (exact Jacobian ranks 3 at fixed y / 4 with y free / 3 torus alone, three generic points; Ω_{(w₁,y)}/Ω_{(w₁′,y)} is a torus character). So the (w₁,y) family is (torus, y), and its torus orbit is V = {(g ρ((w₂(y(g))−Gᵀ)⁻¹), g) : g any rotation of x₀}: **every site-2 rotation g has exactly one site-1 partner**, overlap ∝ χ_μ(D) det₆N[Q_B, Q_vac,(w₁,w₂(y(g)))]/Ψ_Q × det₂^{(1)}[Q_B; y(g)]/minor₁. Checked for a random Gauss-decomposed g (y(g) unique from exact Solve at two w₁'s; 10⁻⁵¹ on 64 states) and for the family at fresh w₁ = −7/3, 9/4 (10⁻⁵¹, 10⁻⁵⁰).
- **Negatives:** no continuous κ-family for the adjoint even with the torus (the first version of the §11 intro claimed one — it was wrong, caused by a Jacobian evaluated at the degenerate point y = z₂; rewritten). The §10 step does **not** iterate: on v₁⊗v₂ the site-1 T's at θ₁−2h have 8×8 ranks {1,4,3,8,8}, full site-2 column space, product locus = seed only.
- **Gotchas hit:** Ω_{(w₁,y)} degenerates at y = z_i (0/0 everywhere); the string `y*)` inside a Wolfram comment closes the comment.
- **NEXT (discuss with Paul):** closed form of the invariant y(g); the site-2 analogue (insertions at t₂(θ₂+h), t₂(θ₂+h)t₁(θ₂−h) + non-localizing partners) on ⟨Ω⊗Ω| and on V; joint site-1 × site-2 span; general reps/L>2 (larger localizing set) followed by the torus analysis. Part A says (g₁, x₀) is not a single natural det, so a 6-dim single-det family is unlikely; sums of dets are the fallback.

## >>> (2026-10-09, later): POSITIVE — two distinct rotations with a det₆ × ladder-1 contour-det overlap (2603.26617 App. G trick for the gl(3) adjoint) <<<
**Notebook:** `new-conventions/L2-rotation-sov.wb` §10 (cells 86–103, saved, not committed). **FULL WRITE-UP + LOG (read first):** `new-conventions/Two-rotations-det-overlap.md`; short form `new-conventions/site1-Qform/6-two-rotations-2603-trick.md`.
- **Result.** With y* = w₂(θ₁₂−h)/(θ₁₂−2h): v₂ = Ω_{(w₁,y*)}, v₁ = Ω_{(w₁,y*)}ρ((w₂−G(z)ᵀ)⁻¹) (= Ω_{(w₂,y*)}ρ((w₁−Gᵀ)⁻¹)), both exact rotations, and ⟨v₁⊗v₂|Ψ_B⟩/⟨x₀|Ψ_B⟩ ∝ det₆[Q_B,Q_vac,(w₁,w₂)]/Ψ_Q × ladder-1 su(2)-like det₂[Q_B; y*]/minor₁ (cell 68's det with E = y*^{−x/h}), state-independent to 10⁻⁵⁰ on 64 states (cell 102). y* checked at five (h,θ,z,w) points; the det identity at P1 only.
- **Mechanism.** T₁₁,T₂₁(θ₁−2h) localize exactly on ⟨Ω_w⊗Ω_w|; T₁₀,T₂₀ add one extra site-2 direction ω₂ = Ω_{(w₁,z₂)}, and Ω₁+κω₂ is a line inside the rotation cone. The product locus of ⟨Ω|(c₀+c₁T₁₀+c₂T₁₁+c₃T₂₀+c₄T₂₁) has one nontrivial point, z-independent, = ⟨Ω|𝒦_{y*} with 𝒦_y the operator with eigenvalue the ladder-1 contour det f_y ∈ span{1,τ₁₀,τ₁₁,τ₂₀,τ₂₁}.
- **Negatives:** no continuous family in the enlarged linear span (t₂(θ₁), site-2 localizers); Ω_w⊗Ω_{w′} is not a ladder-dependent det₆ and its SoV bra wavefunction does not factorize.
- **NEXT:** second numeric point for the det identity; w₃ / finite-y₂ dependence; site-2 analogue; general reps and L>2 (bigger localizing set → possibly a continuous κ-family).

## >>> (2026-10-09): site-1 rotation as a Q-determinant — term-by-term reading of the 12-state sum (Direction 5) <<<
**Notebook:** `new-conventions/L2-rotation-sov.wb` §8 (cells 59–78, saved, not committed). **Spec/outcome:** `new-conventions/site1-Qform/5-bracket-decomposition.md`. Paul's idea: don't fit det ansätze; write the exact 12-state sum as a sum over wavefunctions, fix the site-2 factor, and read the ladder-1 terms for bracket structure.
- **2-plane variables.** On ladder 1, v(θ₁)∥v(θ₁−h) and v(θ₁−3h), v(θ₁−2h), v(θ₁−h) span a 2-plane. With λ_s = v(θ_s)/v(θ_s−h) and v(θ₁−h) = αv(θ₁−3h)+βv(θ₁−2h), the 8 site-1 wavefunctions are {1,β,−α}×{1,λ₁,λ₁²} (minus λ₁²·1). **The overlap of any site-1 covector is an explicit short polynomial in (α,β,λ₁,λ₂)** (table in the spec and in cell 61; exact to 10⁻⁴⁷). It splits as Σφ⁽¹⁾ψ + λ₂Σφ⁽²⁾ψ with φ⁽ⁱ⁾ = c·K⁽ⁱ⁾.
- **Negative results (all exact or 60-digit):** neither group is a rotation (X(φ) rank 3); the L=1 det₃×det₂ formula is not a polynomial identity (fed L=2 Q's on ladder 1 it leaves the 8-dim local span, residual 0.11); the L=1-type 3×3 bracket det is site-1-local iff F₂ ∝ F₁ on {θ₁−2h, θ₁−h}, and that family's image is only a 4-dim cone (10 quadrics) which does not contain the groups (Gröbner membership); sums of two such dets cover the whole local space (uninformative).
- **Verdict:** no single ladder-1 bracket determinant of the natural type; the exact representation is the 12-state sum in the compact polynomial form.
- **Also negative (§9, same day):** two different rotations g₁⊗g₂ — the full-measure det₆ with ladder-dependent vacuum Q's fails (deviations 5–25), and the SoV bra wavefunction of Ω_w⊗Ω_{w′} does not factorize over sites (measure-free test M_{ww′}M_{w′w} ≠ M_{ww}M_{w′w′}, dev 1.6). So no 'sum over SoV integrals with a factorized bra' either.
- **Incident:** a symbol clash (`uS`, the tPoly interpolation points, reused in §8) crashed the laptop; fixed (uSym), memory note written.
- **NEXT (ideas, not started):** brackets involving v(θ₁+h) (outside the 2-plane); two-ladder dets with insertions producing λ₂ (locality analysis of the per-ladder 6×6 in these variables); closed-form K(t) (Direction 4) now reduces to the few rational numbers in the polynomial table. Discuss with Paul before another search.

## >>> (2026-10-02, late): rotations of the SoV vacuum as Q-determinants — write-up `Paul/Mathematica/new-conventions/Rotations-SoV-Q.{tex,pdf,md}` <<<
**Notebooks:** `new-conventions/Rotation-e8-SoV.wb` (L=1) and `new-conventions/L2-rotation-sov.wb` (L=2). Both are saved and not committed.
- **L=1.**
  - Every generic rotation of x0 is ⟨x0|M(w)M(y1,∞)⟩. Here ⟨x0|M(w) = Ω_w, and M(y1,∞) = ρ((y1−G(z)ᵀ)⁻¹), the resolvent.
  - Overlap = det3[Q_B, Q_vac,w]/Ψ_Q × an su(2)-like 2×2 in (Q³,Q²) with insertion y1^{−x/h}.
- **L=2, global rotation.** Overlap = det6 (the usual Ω_{z,w} formula) × a 4×4 su(2)-like det, for all 64 states to 10⁻⁵⁰.
- **L=2, rotation on site 1 only.**
  - It is an exact sum over 12 SoV states, with a rational K map.
  - **RETRACTED:** the claim that ⟨Ω_w⊗e1ᵀ|Ψ⟩ is a single per-ladder 6×6 (ladder 1 = w-vacuum Q's, ladder 2 = numeric F^(2)). Those fits converged to the trivial branch, λ≈0 with det≡0.
  - Redone with the scale fixed (L2-rotation-sov.wb §6b): the best residual is 0.766 (L=2 vacuum Q's on ladder 1) or 0.774 (L=1 vacuum Q's). So this ansatz is ruled out.
  - **Lesson:** fix the scale in det fits.
- **NEXT:** a different Q-form for the site-1 vacuum: modified ladder-1 data, mixed brackets, or a sum of dets. The fallback is the exact 12-state SoV sum. Discuss with Paul before searching.

## >>> FULL WRITE-UP of 2026-10-02: `Paul/Mathematica/new-conventions/Clean-CG-combinations.{tex,pdf,md}` <<<
This file logs the whole contamination / clean-combination thread, from the CG chain to the analytic adjoint formula. Read it first; the entries below are the per-step notes.

## >>> NEWEST (2026-10-02, latest): CLOSED FORMS of the clean coefficients ([2,1,0], [3,1,0], [4,1,0]) <<<

**Notebooks:** `Paul/Mathematica/new-conventions/Clean-coeffs-210-EH.wb` (control), `Clean-coeffs-310-EH.wb` and `Clean-coeffs-410-EH.wb`. All are saved and not committed. The write-up is §3.8 of `Clean-CG-combinations`.
- **No universal combination exists.** In [5,3,0] the k's have full rank (`Minimal-ingredients-530-EH.wb` §6d).
- **Method:**
  - h = 1, θ = (t, 0); f-free one-site covectors, top-killed.
  - Richardson extrapolation over ε₀, 2ε₀, 4ε₀ (ε₀ = 10⁻⁸) gives the exact ε → 0 leading data (error ≤ 10⁻²³).
  - Rational fits in t over 20 points, validated on 5 held-out points.
  - Multiplicity sectors: tested in the form A(t) + B(t)λ_c, with λ_c = t₂(u₀) eigenvalue.
- **Results:**
  - The control reproduces β(t) exactly.
  - **[3,1,0] and [4,1,0]: every contamination/signal ratio is in closed form.** Rational in t for multiplicity-1 states, and P(t) + Q(t)σ_c√D(t) in multiplicity sectors, i.e. linear in the copy's eigenvalue. For example, D = 16(1+3t²) for [4,3,1].
  - **Layer structure:** t-independent universal combinations kill everything except the bottom irrep. These are 6d₁₁ − d₂₂ for [3,1,0], and 16d₁₁ − d₂₂, 72d₁₁ − d₃₃ (layer 1) and v₂ = 72d₁₁ − 9d₂₂ + d₃₃ (layer 2) for [4,1,0].
  - The bottom irrep needs t-rational coefficients, with ρ = (2t−3)/3, (7t²−12)/20 and (11t³+61t−90)/90 for λ₁ = 2, 3, 4.
- **Mechanism (tested where non-trivial):** k̂_B lies in the span of the signals of the HW states above B.
- **√D coefficients identified (`Clean-coeffs-310-EH.wb` §4).** Per unit signal they are affine in τ_c = untwisted t₂(θ₁+h) eigenvalue, with simple coefficients. The [4,3,1] clean combination is v_c = 48[19t(t−1) − τ_c] d₁₁ − [4t(t−1)(3t−23) − (3t−11)τ_c] d₃₂, exact to 10⁻²². Not yet tested on [4,1,0] (its [6,3,1] data would need a fresh 18-min run).
- **Next:**
  - Prove or test the layer/annihilator mechanism more sharply. Find the minimal set of states above B.
  - Explain the discriminants and the integers 6, 16, 72, (72, −9, 1).
  - Extend to [4,2,0] and the λ₁ = 5 reps. These need multiplicity-3 sectors (λ² terms) and more fit points; [4,1,0] already took 18 min.

## >>> (2026-10-02, end of day): site-1 clean combinations for ALL λ₁ ≤ 5 reps ([3,2,0] … [5,4,0]) <<<

**Notebooks:** `Paul/Mathematica/new-conventions/Minimal-ingredients-320-EH.wb` (15 cells) and `Minimal-ingredients-410-EH.wb` (17 cells). They share Λ-generic code. Both are saved, ran successfully and are not committed. The write-up is Section 3.6 of `Clean-CG-combinations.{md,tex,pdf}`.
- **Pipeline:** site 1 + ⟨Ω|Ψ⟩ only.
  - The site-1 insertions come from an in-notebook localization scan.
  - The covectors are computed directly, ⟨Ω|T_λ(R⊗R)⁻¹, so no f_λ are needed.
  - The eigenvalues come from CBR of the scalar t-eigenvalues.
  - Fair test: fit at (10⁻⁵, φ, ω), test at (10⁻⁴, φ) and (10⁻⁵, φ′).
- **[3,2,0]²** (150 digits, about 3 s): all 9 non-top HW states recovered, O(ε) errors. Only [5,3,2] needs per-copy cleaning, mirroring [4,3,1] in [3,1,0].
- **[4,1,0]²** (200 digits, about 25 s): this was the first localization for this rep.
  - **Site-1 local:** T₍₁,₁₎, T₍₂,₂₎, T₍₃,₃₎, T₍₄,₃₎ at θ₁−4h. The old pattern missed T₍₂,₂₎.
  - **f's** (fitted at P1, checked at P2): f₁₁ = (θ₁₂−3h)(θ₁₂−5h); f₂₂ = (θ₁₂−2h)(θ₁₂−3h)(θ₁₂−4h)(θ₁₂−5h); f₃₃ = (θ₁₂−h)(θ₁₂−2h)(θ₁₂−3h)²(θ₁₂−4h)(θ₁₂−5h); f₄₃ = (θ₁₂−h)·f₃₃.
  - **Result:** all 14 non-top HW states recovered. Only [6,3,1] needs per-copy cleaning; [7,2,1] works sector-wide but is weakly conditioned (second singular value 0.0069).
- **[4,3,0]²** (`Minimal-ingredients-430-EH.wb`, 17 cells, 200 digits, about 25 s): this was the first localization for this rep.
  - **Site-1 local:** T₍₁,₁₎, T₍₂,₁₎, T₍₃,₁₎, T₍₄,₁₎ at θ₁−4h, so the pattern holds here.
  - **f's** (fitted at P1, checked at P2): f₁₁ = (θ₁₂−h)(θ₁₂−5h); f₂₁ = (θ₁₂−h)(θ₁₂−3h)(θ₁₂−5h); f₃₁ = (θ₁₂−h)(θ₁₂−2h)(θ₁₂−3h)(θ₁₂−5h); f₄₁ = (θ₁₂−h)·f₃₁.
  - **Result:** all 14 non-top HW states recovered. Only [7,5,2] needs per-copy cleaning; it is the dual image of [4,1,0]'s [6,3,1].
- **[5,1,0]²** (`Minimal-ingredients-510-EH.wb`, 19 cells, 250 digits; about 2.5 min plus 1 min for the 6b refit). This is the first λ₁ = 5 rep: dim 35, 1225 states.
  - **Site-1 local:** T₍₁,₁₎, T₍₂,₂₎, T₍₃,₃₎, T₍₄,₄₎, T₍₅,₄₎ at θ₁−5h. The [λ₁,1,0] rule is T₍ₖ,ₖ₎ (k < λ₁) + T₍λ₁,λ₁−1₎.
  - **f's** (t = θ₁₂/h): f₁₁ = (t−4)(t−6); f₂₂ = (t−3)(t−4)(t−5)(t−6); f₃₃ = (t−2)(t−3)(t−4)²(t−5)(t−6); f₄₄ = (t−1)(t−2)(t−3)²(t−4)²(t−5)(t−6); f₅₄ = (t−1)f₄₄.
  - **Result:** all 17 non-top HW states recovered. Only [8,3,1] needs per-copy cleaning; for [λ₁,1,0] the per-copy sector is always [λ₁+3,3,1].
  - **Convergence:** the 10⁻⁶ refit (Section 6b) shows O(ε) convergence, with errors 10⁻⁵·²–10⁻⁷·³ at (10⁻⁶, φ′).
- **[5,4,0]²** (`-540-EH.wb`, 250 digits, about 3 min): T₍ₖ,₁₎ for k = 1..5. Only [9,7,2] needs per-copy cleaning (the dual of [8,3,1], as predicted).
- **[5,2,0]²** (`-520-EH.wb`, 21 cells, 300 digits, about 9 min; 28 HW states, two multiplicity-3 sectors).
  - **Site-1:** T₁₁, T₂₂, T₃₃, T₄₃, T₅₃.
  - **Per-copy required:** [9,3,2], [8,5,1], [8,4,2], [7,5,2].
  - **New Section 6c:** weakly conditioned full-rank sectors ([7,6,1] σ₂ = 0.025, [7,4,3] σ₂ = 0.019) amplify the O(ε) error by about 1/σ₂. Per-copy cleaning fixes them.
- **[5,3,0]²** (`-530-EH.wb`, 21 cells, 300 digits, about 8.5 min): T₁₁, T₂₂, T₃₂, T₄₂, T₅₂. Per-copy needed for [9,5,2], [8,7,1], [8,6,2], [8,5,3], the exact duals of [5,2,0]'s.
- **All ten λ₁ ≤ 5 reps (λ₁ > λ₂ > 0) are now done.** All non-top HW states are recovered with O(ε) errors (≤ 10⁻⁵·⁵ at ε = 10⁻⁶ with the best cleaning choice). The summary is Section 3.7 of `Clean-CG-combinations`.
  - **Unified site-1 rule:** with d = λ₁−λ₂, the insertions are T₍ₖ,ₖ₎ (k ≤ d) and T₍ₖ,d₎ (d < k ≤ λ₁), all at θ₁−λ₁h. The leading prefactor is f₁₁ = (θ₁₂−dh)(θ₁₂−(λ₁+1)h).
- **Not done:** site 2 for λ₁ ≥ 4; closed forms of the clean coefficients; a second chain point for the pipeline.

## >>> (2026-10-02, latest): adjoint clean combination derived analytically <<<

**Notebook:** `Paul/Mathematica/new-conventions/Adjoint-clean-analytic.wb`. Exact, symbolic t = θ₁₂/h, E basis, about 1 s.
- **Mechanism.** Ω lies in the 27, so (W − T^top Ω) vanishes at ε⁰ and top-killing removes the O(ε) contamination identically. The leading contamination is ε²(W − T^top Ω)₁·ψ₁; only first-order perturbation theory is needed.
- **Which states need cleaning.** States with an O(ε) signal ([4,1,1], [3,3,0], octets) are automatically clean. The singlet has an O(ε²) signal; in raw eigenvalues K₁₁/S₂₁ = −1/3 and K₂₁/S₂₁ = (2t−3)/3, independent of φ and ω.
- **The clean combination.** 𝒞 = T₂₁(θ₁−2h) + h(2θ₁₂−3h)T₁₁(θ₁−2h), with Q-number (𝒞^B − 𝒞^top)⟨Ω|Ψ⟩. This is the fitted β = h(2t−3)/(t−1), now derived.
- **In Q's.** T₁₁ = 3h²χ₃(θ₁₂−h)(θ₁₂−3h)·r and T₂₁ = −3h³χ₃(θ₁₂−h)²(θ₁₂−3h)·ρr, with r = Q_c(θ₁−h)/Q_c(θ₁−2h) and ρ = Q^j(θ₁)/Q^j(θ₁−h); checked to 10⁻⁵⁰ on 18 states. Hence the clean functional ∝ {r[(2θ₁₂−3h) − (θ₁₂−h)ρ] − z₁⁻¹[(2θ₁₂−3h) − (θ₁₂−h)z₃]}·⟨Ω|Ψ⟩.
- **Rescaling (Section 3).** Rescaling the monodromy by g(u) changes only the number g(θ₁−h) in the combination.
  - ν₁, ν₂, ν₃ and the f-normalization all leave t-dependence (ν₂ is singular).
  - Rescaling each Lax by (u−θ_s−h/2) gives 𝒞̃ = T̃₂₁ − (4/3)T̃₁₁, which is t-free. The half shift on site 2 is forced.
- **Next:** the same derivation for [3,1,0]/[4,2,0], where deeper states need higher orders (ψ₂, …); also test whether the half-shift normalization makes those clean conditions t-free too.

## >>> (2026-10-02, later): minimal notebook for the contamination; Paul experimenting <<<

**Notebook:** `Paul/Mathematica/new-conventions/Minimal-ingredients-EH.wb` (16 cells, self-contained, clears the kernel, runs in ~10 s, saved, not committed). Built for Paul to experiment with. It contains only the EH generators, `tL2`/`tL1`, T₍₁,₁₎/T₍₂,₁₎ (two-site and one-site), `Rrot`, `Om1`/`Om2`, `x0`, the small-twist HW-limit states, the four localization identities, c_λ, and `Qside` (= T_λ^B⟨Ω|Ψ⟩/f_λ computed from the eigenvector; no Q's).
- **Exact split, at any ε:** `Qside` = c_λ·P_HWΨ (signal) + c_λ·(R⊗R)(1−P_HW)Ψ (contamination), because (R⊗R)P_HW = P_HW.
- **It reproduces every κ in Section 11.** The chain residual is about 10⁻⁴⁶ at 80 digits.
- **`showContam` decomposes the contamination by weight ν.** At [3,3,0] (2,1), the whole contamination is ν = [4,2,0]: c^{(420)} ~ ε⁰ times [(R⊗R)Ψ⊥]₄₂₀ ~ ε⁻², the same order as the signal ε¹·ε⁻³. The octet behaves the same way, and C removes it completely.
- **Discussion conclusions (not in any notebook):**
  - **Expanding in w cannot replace a weight projector.** Every ν term is homogeneous of the same torus degree top−μ in (r₁,r₂,r₃).
  - **A third independent r would add nothing.** The missing direction is a torus direction, and ρ = r₂/(r₁r₃) already depends on w₂.
  - **Only insertions diagonal on Ψ (combinations of T's) can serve as the "RHS projector."**
- **Section 7 (cells 17–23): a + b T₍₁,₁₎ + c T₍₂,₁₎ removes the contamination.** The bare Ω row is the covector e₁⊗e₁, which is pure top weight.
  - **Killing the top weight with a = −3h²z₂z₃(b − h z₃ c) is enough for every state except the singlet.** κ → 1 for [4,1,1], [3,3,0] and both octets, with either T alone.
  - **The two top-killed covectors have rank 2 on the octet HW space**, so the octet mixing is fixed.
  - **The singlet needs b/c = β = h(2t−3)/(t−1).** Its contamination ratios are K_b/S = −(t−1)/(3h) and K_c/S = (2t−3)/3. These were recognized at P1 and **validated at P2**, and are φ- and ω-independent.
  - **The single combination v* = (a, β, 1) gives κ → 1 for all five non-top HW states at P1 and P2.** So the HW CG data come directly from Q-expressions, with no reconstruction.
  - **Not done:** [4,2,0], and deriving β(t).
- **[3,1,0]² (`Minimal-ingredients-310-EH.wb`, 18 cells; generic in Λ, with Λ = {2,1,0} as the control).** Sparse 225-dim matrices at 120 digits; a fresh run takes about 20 s. Saved, not committed.
  - **Ingredients:** site 1 (T₍₁,₁₎, T₍₂,₂₎, T₍₃,₂₎ at θ₁−3h) + ⟨Ω|Ψ⟩ only. Localization holds to 10⁻¹⁰⁶.
  - **Top-killing is exact and works for any rep:** ⟨Ω|Ψ⟩ is the top component of (R⊗R)Ψ⊥. It cleans [6,1,1], [5,3,0] and [5,2,1]×2.
  - **The deep states keep contamination, but it is structural** (φ- and ω-independent).
  - **v₁ = 3d₁₁ − 2d₂₂** is clean on [5,3,0], [5,2,1]×2, [4,4,0], [4,3,1]×2 and [4,2,2].
  - **[5,2,1] mixing is fixed (rank 2). [4,3,1] mixing is not.** The only combination clean on both copies is v₁, and v₁'s direction on HW_[4,3,1] is w-independent (exact rank 1).
  - **Appendix A (independent scan):** the only nontrivial site-1-local CBR shapes (λ₁ ≤ 4, c = −6..2) are T₍₁,₁₎, T₍₂,₂₎, T₍₃,₂₎ at θ₁−3h.
  - **Section 8: adding site-2 A = t₂(θ₂+h) and B = t₂(θ₂+h)t₁(θ₂−h) fixes everything.** Their covectors are computed as ⟨Ω|X(R⊗R)⁻¹ and need no prefactors.
    - All contamination ratios are structural.
    - Clean spaces per weight sector (null space of the copies' k-vectors) have rank 2 on both [5,2,1] and [4,3,1].
    - **End to end:** the Q-side numbers alone recover P_HWΨ, scale included, for all 9 non-top HW states. Errors are O(ε), ≲10⁻⁶ at ε = 10⁻⁵.
    - **Caveat:** the clean combinations are numerical, fitted from the states' own k-vectors. They are structural, but not yet in closed form. One chain point only.
- **[4,2,0] blind spot corrected (`Blindspot-420-check-EH.wb`, 6 cells; one-site objects plus two-site HW vectors).**
  - **It is not blind.** Site-1 T₍₃,₂₎ and T₍₄,₂₎ have nonzero overlap with HW_[5,4,3] (∝ε³), and T₍₄,₂₎ with HW_[4,4,4] (∝ε⁴); the localized u stay O(1).
  - **What the old test showed.** Its "overlap sv = 1.000" (printed at 3 digits) only meant near-containment in the numerical null space. A correction note has been added to `L2-420-localization.md`.
  - **Done, see below.**
- **[4,2,0]² clean combinations (`Minimal-ingredients-420-EH.wb`, 13 cells, 250 digits, about 1 min). Saved, not committed.**
  - **Insertions:** site-1 T₁₁, T₂₂, T₃₂, T₄₂ at θ₁−4h; site-2 t₂(θ₂+h), T₂₂(θ₂), T₂₂(θ₂)t₁(θ₂−2h); plus ⟨Ω|Ψ⟩.
  - **Bare-frame bookkeeping:** S = W·P_HWΨ and K = W·Ψ⊥ − (W·e₁₁)⟨Ω|Ψ⟩, so no R⊗R is needed.
  - **Leading order is enough:** every state has p_S = p_K or better.
  - **Contamination is structural**, with an O(ε) spread.
  - **Sector clean spaces fix every multiplicity space except [6,4,2]** (mult 3; the common clean rank is 2). Per-copy cleaning fixes it.
  - **Fair end-to-end test:** the Q-numbers recover all 18 non-top HW states, scale included, with O(ε) errors (fitted at ε = 10⁻⁶, tested at 10⁻⁴/10⁻⁵ and at another φ). **That includes [5,4,3] and [4,4,4].**
  - **Precision trap:** 120 digits gave fake "0" errors (zero accuracy) for [4,4,4] at ε ≤ 10⁻⁵. Always check Precision.
  - **Open:** closed forms for the clean coefficients (they should be constants in t, like β(t)); a second chain point; a principled way to pick the [6,4,2] per-copy combinations.
- **Site 2 is not essential (2026-10-02).** With site 1 + ⟨Ω|Ψ⟩ only, plus per-copy cleaning where the sector-wide clean space lacks rank, all HW states are recovered (fair test, O(ε) errors).
  - **[2,1,0]²:** sector-wide is enough.
  - **[3,1,0]²** (`Minimal-ingredients-310-EH.wb` Section 9): only [4,3,1] needs per-copy cleaning.
  - **[4,2,0]²** (`Minimal-ingredients-420-EH.wb` Section 3): [7,3,2], [6,5,1] and [6,4,2] need per-copy cleaning; [5,4,3] and [4,4,4] work sector-wide.
  - **What site 2 buys:** fewer sectors with copy-specific coefficients (none for [3,1,0]; only [6,4,2] for [4,2,0]).
  - **Trap:** an already-clean insertion has k = O(ε), so use an absolute cut (`kTol`), not a relative one.
- **Site-1 eigenvalues in dual Q's (2026-10-02; `310` Section 10, `420` Section 4).** Let v(u) = (Q¹,Q²,Q³)(u) be the dual Q's.
  - **Structure near θ₁:** v(θ₁−Λ₂h) ∥ … ∥ v(θ₁), and v(θ₁−(Λ₁+1)h) … v(θ₁) lie in one 2-plane.
  - **Formula:** T_λ(θ₁−Λ₁h) = P_λ · det₂[v(θ₁+(λ₁−Λ₁)h), v(θ₁+(λ₂−Λ₁−1)h)] / det₂[v(θ₁−Λ₁h), v(θ₁−(Λ₁+1)h)], with P_λ = ∏_boxes ν₃(θ₁−Λ₁h+(c−r)h). This is the L=2 version of L=1 row 13.
  - **Validated:** exactly for 34 states of [3,1,0]² and 91 of [4,2,0]², all shapes with λ₁ ≤ Λ₁, every row pair.
  - **[4,2,0] explicit:** T₂₂/f₂₂ = κD₂₃/D₄₅, T₃₂/f₃₂ = −2hκρ₁D₂₃/D₄₅, T₄₂/f₄₂ = 2h²κρ₀D₂₃/D₄₅, with κ = 240h⁴(θ₁₂−3h)(θ₁₂−4h)/((θ₁₂−h)(θ₁₂−2h)); T₁₁/f₁₁ = 20h²(θ₁₂−4h)/(θ₁₂−2h)·D₃₄/D₄₅. Here ρ_k = Q^j(θ₁−kh)/Q^j(θ₁−2h).
  - **Not done:** the 6×6 ⟨Ω|Ψ⟩ formula has not been validated beyond the adjoint.
- **Single-Q form (`420` Section 5, 91 states, ≥225 digits).**
  - **Dual QQ, confirmed:** z_b q^a(u)q^b(u+h) − z_a q^b(u)q^a(u+h) = (z_b−z_a) G(u) q_c(u+h), with G = ∏_s(u−θ_s+h)(u−θ_s+2h). Note the shift +h on the direct q_c.
  - **Direct parallelism:** the direct Q vectors are parallel at θ₁−4h, −3h, −2h.
  - **Insertions:** T₁₁/f₁₁ = 10h²χ₃ r₁; T₂₂/f₂₂ = 40h⁴χ₃² r₁r₂; T₃₂/f₃₂ = −80h⁵χ₃² r₁r₂ρ₁; T₄₂/f₄₂ = 80h⁶χ₃² r₁r₂ρ₀. Here r₁ = Q_c(θ₁−3h)/Q_c(θ₁−4h) and r₂ = Q_c(θ₁−2h)/Q_c(θ₁−3h) (direct, any c), and ρ_k = Q^j(θ₁−kh)/Q^j(θ₁−2h) (dual, any j). The θ₁₂-dependence cancels against f.
  - **Next:** rewrite the fitted clean coefficients in (r₁, r₂, ρ₁, ρ₀); test [3,1,0].
- **Clean combinations, first look (`420` Section 6).**
  - **Top-kill** = eigenvalue difference: (X^B − X^top)⟨Ω|Ψ⟩.
  - **Leading signal ratios s̃** over (T₁₁, T₂₂, T₃₂, T₄₂)/f are simple rationals, e.g. [7,5,0] (1,2,−2,1) and [6,3,3] (0,0,−1/5,1).
  - **Contamination k̃:** 0 for [8,3,1] and [7,5,0] (only the top above them). [8,2,2]: (t−1)·s̃^[8,3,1]. [6,6,0]: (t−1)·s̃^[7,5,0]. In words, it is the signal pattern of the irrep one simple root above, times (t−1); this holds pointwise in w, and the raw factor is w-dependent.
  - **Deeper states:** k̃ is not rational at P1, presumably mixing several irreps above, including the multiplicity copies (algebraic data).
  - **Next:** expand k̃ in the s̃'s of the states above, using extra (non-localized) insertions so the span test is meaningful.

## >>> (2026-10-02): independent polynomial Baxter spectral solver (Section 14) <<<

**Read first:** Section 14 of `Paul/Mathematica/new-conventions/L2-diagonal-frame-overlap.md`. Notebook: `L2-diagonal-frame-overlap-EH.wb`, explicitly saved at revision 675 with 114 cells. Section 13 (cells 90–108) remains the 15/15 audit; new Section 14 is cells 109–114. Cells 110, 112 and 114 were individually executed successfully after fixes. No commit was made.

- **The solver no longer needs physical spectral input.** `specAnsatz` builds monic dual Q polynomials of degrees μⱼ, transfer polynomials t₁,t₂ with weight/twist-fixed leading terms, and known quantum determinant t₃. The ten unknowns (six Q coefficients, four transfer coefficients) satisfy polynomial dual Baxter identities. No physical matrices or eigenvectors enter this solve.
- **Numerical solve.** `specSolve` starts with exact-rational `NSolve`, then `specNewton` performs four scaled overdetermined Newton/least-squares corrections at working precision `pr+80`, seeded only by the Baxter roots. It rounds to the requested precision and checks the stored coefficients by rationalizing them in the equations. Directly increasing `NSolve` to `2pr` was too slow and was aborted; precision comes from polynomial refinement, not seed padding. No kernel operation remains in flight.
- **Validation is downstream.** `specReconstruct` uses spectral solutions, local 8×8 covectors, determinant overlaps, three bra twists, and x₀ normalization. Only afterward does `specValidate` build physical transfer matrices and check eigenvector equations coefficientwise with guard precision. This complete finite-twist test does not use a physical `Eigensystem` as spectral input.
- **Finite-twist result.** At h=1, θ=(5/2,0), z=(2,3,5), all 19 weight sectors yield 64 solutions, matching physical weight dimensions; Q reconstruction succeeds for all 64. The pre-final-refinement run reported max Baxter residual 3.34×10⁻⁹⁴, physical eigenvector residual 4.02×10⁻⁹², and normalization residual 9.28×10⁻⁹³ at 100 digits. The final code also succeeded, but its exact residuals were not recorded; describe it as roughly 90-digit validation rather than assigning new figures. Bra twists: (7/3,11/4,4/7), (9/4,13/5,3/7), (12/5,17/6,5/8).
- **Near-untwisted result.** Along rational paths z=1+εφ and w=1+εω (the same first-order direction as exponentials, but not exponential twists), all 21 states in the five HW-weight sectors (dimensions 1,2,2,6,10) now pass; six HW branches are selected by raising-operator residual. At P1, ε=10⁻⁶/pr=200, max Baxter residual 1.873621426×10⁻¹⁷³, physical residual/norm 1.809978519×10⁻¹³², and CG direction error 1.864923519×10⁻⁷. At P1 ε=10⁻⁹/pr=280 these are 4.779807849×10⁻²⁴², 3.397139631×10⁻¹⁸¹, and 1.864923680×10⁻¹⁰. At P2 ε=10⁻⁹/pr=280 with alternate φ=(-2/5,3/8,7/11), they are 8.105060246×10⁻²⁴⁶, 2.811004503×10⁻¹⁸¹, and 2.638405385×10⁻¹⁰. Primary φ=(1/3,−2/7,5/11); the three ω directions are recorded in Section 13. The small CG-direction error is O(ε), distinct from the much smaller finite-twist residual.
- **Scope.** This supplies a numerical polynomial Baxter solver and validates finite examples; it is not an analytic closed-form or general completeness proof, nor an arbitrary-representation implementation. The untwisted eigensystem is only a verification reference in the limit comparison.

## >>> NEWEST (2026-10-01): return to [2,1,0]² — symbolic factors and full Q reconstruction <<<

**User priority:** Paul asked to return to [2,1,0]², verify the earlier work and improve it. The [4,2,0] notebook was opened read-only; no new [4,2,0] calculation or edits were made.

**Read first:** Section 13 of `Paul/Mathematica/new-conventions/L2-diagonal-frame-overlap.md`. Code: `L2-diagonal-frame-overlap-EH.wb`, cells **90–108**. Rerun **3–14, then 90–108**. The audit ledger has **15 passing checks**.

- **Real correction:** site-2 B=t₂(θ₂+h)t₁(θ₂−h) has scalar coefficient **β=6h⁴χ₃(θ₁₂−h)(θ₁₂−3h)**. The old 3h³ expression only coincided at h=1/2. Cell 81 and the write-up are corrected.
- **All four localization formulas are now symbolic identities.** Undo R, use an arbitrary lower-triangular auxiliary twist (six independent entries), prove entrywise at h=1, θ=(d,0), then restore h and θ by homogeneity/translation. Also checked exactly at three rational chain/twist points.
- **Simpler overlap:** with D_B=det N[B,vac_w], P_B=Ψ_Q^B(x₀), B₀=e₁⊗e₁,
  `ΩΨ_B/(x₀Ψ_B) = D_B P_B₀ / [D_B₀ P_B (z₂−z₁)²]`.
  This cancels y₀, Φ_Q and c₀. Either bracket convention gives the same ratio.
- **The previously missing Q-side test is closed.** Reconstruct using this determinant ratio × transfer eigenvalues recovered from dual Q's, **not M·Ψ**. All 64 baseline states and all six small-twist HW branches pass. At P1/P2, ε=10⁻⁶/10⁻⁹, 240/300 digits, relative errors are below 4×10⁻¹⁶¹ / 3×10⁻¹⁹¹. An independent 160-digit run agrees with the 240-digit result to 7.73×10⁻⁷⁸.
- **Actual CG limit checked:** direction errors fall as O(ε), below 2.64×10⁻¹⁰ at ε=10⁻⁹, including an alternate φ and both octets.
- **Three w values suffice and are minimal for the four operator rows + x₀.** Exact ranks 1,2,2,6,10. Two w yield at most 9 equations for the 10d singlet weight sector. This is not a minimality claim if bare Ω or other rows are added.
- **Exact site-2 HW overlap table** in Section 13e and the `.md`, in terms of the six entries of g⁻ᵀdiag(z)gᵀ.

**Scope / next:** At the time of Section 13, test Q-data still came from Baxter equations fed transfer eigenvalues from the test eigensystem. Section 14 supersedes that limitation for the validated adjoint-square examples with a numerical polynomial Baxter solver; it does not establish general completeness. The historical summary in cells 88–89 includes old rejected approaches; use the new audit ledger. No commit made.

## >>> NEWEST (2026-10-01, final): [4,2,0]², φ-independence, blind spots <<<

**Read first.** `Paul/Mathematica/new-conventions/L2-420-localization.md`; the notebook is `L2-420-localization-EH.wb`, with 20/20 checks passing. Also Section 7 of `L2-320-localization-EH.wb` (φ).

**Results:**
- **[4,2,0], site 1.** T₍₁,₁₎, T₍₂,₂₎, T₍₃,₂₎ and T₍₄,₂₎ at θ₁−4h.
- **[4,2,0], site 2.** t₂(θ₂+h) and T₍₂,₂₎(θ₂).
- **The site-1 pattern across all four reps:** T₍₁,₁₎ plus T₍ₖ,λ₁−λ₂₎ for k = λ₁−λ₂..λ₁, all at θ₁−λ₁h, with f₍₁,₁₎ = (θ₁₂−(λ₁−λ₂)h)(θ₁₂−(λ₁+1)h).
- **[4,2,0] reconstruction.** Localized data + ⟨Ω|Ψ⟩ recover 16 of 19 HW states.
  - **[5,4,3]×2 and [4,4,4] are structurally invisible:** HW_μ lies entirely in the null space, so no ε-order, w or φ helps.
  - **Adding non-localized t-rows** ⟨Ω|t_a(u)| at generic u, which are still Q-computable, recovers all 19.
- **φ-independence.** The ε → 0 limit basis does not depend on φ. It is the untwisted inhomogeneous eigenbasis on HW_μ.
- **Open (Paul to decide).** Are non-localized t's acceptable? If not, a new localized ingredient is needed for [5,4,3].
- **Open, Paul picking this up 2026-10-02: the site-1/site-2 asymmetry.** Site 1 has many more localized objects than site 2 (spans for [2,1,0], [3,1,0], [3,2,0], [4,2,0]: site 1 = 3, 4, 4, 5; site 2 = 3, 3, 4, 4). Paul will look for more site-2 objects.
  - **What the search covered:** CBR shapes λ₁ ≤ 6 at θ₂ + c·h for c = −9..9, which includes the fused t\* objects; products by t₁ and t₂ to depth 2; all pairs and triples of t's at θ₂ + ℤh.
  - **What it did not cover:** non-integer shifts, derivatives in u, and the site-reversed monodromy.

## >>> (2026-10-01, end): [3,2,0]² site localization <<<

**Read first.** `Paul/Mathematica/new-conventions/L2-320-localization.md`. The notebook is `L2-320-localization-EH.wb`, the [3,2,0] copy of the [3,1,0] notebook; it reruns in about 32 s and 16/16 checks pass.

**Results:**
- **Site 1.** **T₍₁,₁₎, T₍₂,₁₎ and T₍₃,₁₎ at θ₁−3h**, which is Paul's guess.
- **Site 2.** t₂(θ₂+h) and **T₍₂,₂₎(θ₂)**.
  - The [2,1,0]/[3,1,0] product t₂(θ₂+h)t₁(θ₂−h) is **not** local here.
- **The f's swap sites relative to [3,1,0]**, consistent with [3,2,0] ≅ [3,1,0]*⊗det³.
- **Reconstruction.** site1 + site2 + ⟨Ω|Ψ⟩ + x₀ recovers **all 10 HW states exactly**.
  - Without ⟨Ω|Ψ⟩, the bottom irrep [4,3,3] is missed, just as [3,3,2] was for [3,1,0].
- **ε-depth (Section 6).** These are the orders beyond leading needed to fix each CG state, using 8 values of w:
  - [6,4,0]: 0;
  - [6,3,1], [5,5,0], [5,4,1]: 1;
  - [6,2,2], [5,3,2], [4,4,2]: 2;
  - [4,3,3]: 5.

  The data are needed through ε^{N−k}, where Ψ ~ ε⁻ᵏ with k = 4..9. This has not been run for [2,1,0] or [3,1,0].

## >>> (2026-10-01, very latest): [3,1,0]² site localization <<<

**Read first.** `Paul/Mathematica/new-conventions/L2-310-localization.md`. The notebook is `L2-310-localization-EH.wb`; it is self-contained and generic in Λ, with [2,1,0] as the control, and reruns in about 25 s. 18/18 checks pass.

**Results:**
1. **Site 1.** Everything is at **θ₁−3h**: **T₍₁,₁₎, T₍₂,₂₎ and T₍₃,₂₎**.
   - **Not T₍₃,₁₎ or T₍₂,₁₎**, at any integer shift.
   - **Factorization.** Each is f_λ × (one-site T_λ at θ₁) ⊗ Ω_zw, with f₍₁,₁₎ = (θ₁₂−2h)(θ₁₂−4h); f₍₂,₂₎ and f₍₃,₂₎ are in the `.md`.
2. **Site 2.** The same as [2,1,0]: t₂(θ₂+h) and t₂(θ₂+h)t₁(θ₂−h), with prefactor (θ₁₂−h)(θ₁₂−4h).
   - The search covered all CBR shapes, t\*, and pairs and triples of t's.
3. **Reconstruction.** site1 + site2 + x₀ reconstruct **9 of the 10 HW states exactly**, including both multiplicity-2 pairs ([5,2,1] and [4,3,1]).
   - **The exception is [3,3,2]**: the (3,3,2) rank is 16 of 18.
   - **Adding the bare ⟨Ω(w)|Ψ⟩ row** makes all 10 exact.

**Decided (Paul):** ⟨Ω|Ψ⟩ is fine as data. So site1 + site2 + ⟨Ω|Ψ⟩ + x₀ gives all 10 HW states of [3,1,0]² exactly.

**Next:**
- Other reps ([3,2,0], [4,2,0]): just change Λ.
- A second (h, θ) point.

## >>> (2026-10-01, latest): site-2 localization closes the octet/singlet gap (Section 12) <<<

**Read first.** "Section 12" of `Paul/Mathematica/new-conventions/L2-diagonal-frame-overlap.md`. The code is in `L2-diagonal-frame-overlap-EH.wb`, cells 79–87, and needs cells 3–14 plus the Section 11 helpers.

**Results:**
1. **t\* reduces to ordinary t's.** t\* (from T⁻¹ = ℒ₂⁻¹ℒ₁⁻¹) satisfies t\*₁(v) ∝ t₂(v+2h) and t\*₂(v) ∝ t₁(v+2h).
2. **The site-2-local objects.** The only ones at θ₂ + ℤh are:
   - **⟨Ω|t₂(θ₂+h)**, which is ∝ t\*₁(θ₂−h);
   - **⟨Ω|t₂(θ₂+h)t₁(θ₂−h)**.
   - **Negative findings.** No CBR T₍₂,₁₎ works. The search covered all CBR shapes and all products of two t's.
3. **Spectator factors (one-point fit).** The site-1 factor is f₍₁,₁₎ = (θ₁₂−h)(θ₁₂−3h), the same as on site 1. The product version has the extra piece α·(one-site t₂t₁) + β·Ω_zw.
4. **What the site-2 data recover.** Site-1 data (Section 10) + site-2 data + ⟨x₀|Ψ⟩ = 1 reconstruct **all six HW CG states exactly at small ε**: both octets and the singlet scale. This closes the Section 11 gap.
   - **Bonus.** Adding the bare ⟨Ω(w)| row alone already fixes the octet, but not the singlet scale.
5. **Caveat.** The reconstruction used M·Ψ in place of the Q side. That is exact by identity, but the full Q chain was not re-run.

**Next (Paul to decide):**
- Closed-form site-2 CG overlaps (the 10d analogue).
- A second (h, θ) point for the fits.
- [3,1,0]².

## >>> (2026-10-01, later): ε→0 of the mixed-twist overlap and the CG chain (Section 11) <<<

**Read first.** "Section 11" of `Paul/Mathematica/new-conventions/L2-diagonal-frame-overlap.md`, the full record. The code is in `L2-diagonal-frame-overlap-EH.wb`, cells 50–78, and needs only cells 3–14 run first. Two chain points: P1 = (1/2, (2/7,−3/5)) and P2 = (1/3, (1/5,−4/9)). Twists are z = e^{εφ}, w = e^{εω}, and t = θ₁₂/h.

**Results:**
1. **The ratio ρ.** ⟨x₀|y₀⟩²ρ_B → ε^# X_B ∏(ω_i−φ_j)^{max(d_j−D_i,0)}, with # = 8, 8, 8, 7, 7, 6 for [4,2,0], [4,1,1], [3,3,0], [3,2,1]±, [2,2,2].
   - **Constants:** X = 1, 1, 1; (2±√(4+5t²))/t² for the octets; −8/(t²−1) for the singlet.
   - **The ε⁸ is normalization only:** ⟨x₀|e₁e₁⟩ = (z₂−z₁)² and ⟨Ω|y₀⟩ ~ ε⁻¹⁰.
2. **The y₀-free ratio** R_B = ρ_B/ρ_[4,2,0]. The [4,1,1] and [3,3,0] rows are exact monomials in (w_i−z_j) at finite twist.
3. **The CG chain.** c_λ·(R⊗R)Ψ = T_λ⟨Ω|Ψ⟩/f_λ holds exactly, but "(R⊗R)Ψ → HW" fails at leading order.
   - **The cause:** c_λ's [4,2,0] component is O(1), and it pairs with the O(ε) correction to Ψ at the same order as the CG term.
   - **The size:** the contamination factors κ are φ- and ω-independent. For [3,3,0] κ = (t+1)/2.
4. **Partial fix.** C = c₂₁ + h z₃ c₁₁ has no [4,2,0] component and gives κ = 1 for [4,1,1] and both octets.
5. **What the allowed data can recover** (T₍₁,₁₎ and T₍₂,₁₎ at θ₁−2h for all w, plus x₀):
   - [4,2,0], [4,1,1] and [3,3,0] are reconstructed exactly.
   - **The octet mixing is undetermined:** there is an exact null vector at every ε, and it tends to a vector in the octet HW space.
   - The singlet's scale is not fixed. Coincident φ and ω don't help.

**Next (Paul to decide).** Fixing the octet needs one more independent ingredient. Paul's constraint is that only site-1-localized T's may be used. **Resolved by site-2-localized t's (Section 12, entry above).**

## >>> (2026-10-01, late): localized T_λ(θ₁−2h) action → highest-weight CG data <<<

**Read first.** The section "Section 10 (EH notebook only)" of `Paul/Mathematica/new-conventions/L2-diagonal-frame-overlap.md`. The code is in `L2-diagonal-frame-overlap-EH.wb`, Sections 10–10f (cells 28–49). Those sections need only cells 3–14 run first.

**Results:**
- **Localization.** Build the t's with the ket's twist diag(z). Then ⟨Ω_zw⊗Ω_zw|T_λ(θ₁−2h) = f_λ ⟨v_λ|⊗⟨Ω_zw| for λ = (1,1), (2,1), with ⟨v_λ| = ⟨Ω_zw|T_λ^{L=1,θ₁}(θ₁−2h).
  - **Exact** at the one numeric point.
  - **Site-2 factors, fitted at that point:** f_(1,1) = (θ₁₂−h)(θ₁₂−3h) and f_(2,1) = (θ₁₂−h)²(θ₁₂−3h).
- **Undoing R.** ⟨v_λ|R⁻¹ = e₁ᵀT_λ^{L=1}[g⁻ᵀ diag(z) gᵀ](θ₁−2h). This is symbolic in z, w and h.
- **Highest-weight overlaps.** For c_λ = ⟨v_λ|R⁻¹⊗e₁ᵀ the overlaps with the HW states of [2,1,0]² are known in closed form; the table is in the `.md`.
  - **(1,1):** misses [4,1,1] and the singlet, as forced by weights.
  - **(2,1):** reaches everything.
- **Small twists.** For z = e^{εφ}, w = e^{εω}, the leading-ε vectors (both λ, ω varied) span every weight sector of e_jᵀ⊗e₁ᵀ.
  - The 2-dimensional sector (3,2,1) needs both λ.
  - HW states are determined by their e_j⊗e₁ components.
  - Hence these vectors generate all the HW CG coefficients (Paul, 2026-10-01).

**In flight elsewhere.** Paul is doing the overlaps → CG step with another agent: assembling ⟨Ω|Ψ_B⟩ formula × T_λ^B(θ₁−2h), then R, then ε → 0.

**Open here:**
1. **Second numeric point.** Localization and the f_λ fit at a second numeric point, or symbolically (needs 64×64 objects).
2. **The pattern.** Test # = max(4−μ₁, μ₃) for the ε-orders.
3. **Other reps.** [3,1,0]² and others.

**Traps:**
- **`MatrixPower[nilpotent, 0]` fails.** Start the exp series from the identity.
- **Arbitrary-precision SVD on rank-deficient matrices** loses precision tracking. Use machine precision with a 10⁻⁹ tolerance for ranks.
- **The quadratic Casimir can't separate [4,1,1] from [3,3,0].**
- **Notebook state.** Sections 5–9 of the EH notebook were not re-run this session, so the summary table is stale.

## >>> CURRENT (2026-10-01, end of day): new conventions; L=2 overlap formula; Hermitian (EH) generators <<<

**Read first, in this order:**
1. `Paul/Mathematica/new-conventions/L2-diagonal-frame-overlap.md`: the L=2 result and the EH variant.
2. `Paul/Mathematica/new-conventions/Hermitian-structure-abstract.md`: whether the t's commute with their conjugates. The answer is no at L=2 with a twist.
3. `Paul/Mathematica/new-conventions/Companion-B-C-diagonalize-clean.md` (with `.tex`/`.pdf`): the L=1 reference.

**Conventions (2026-10-01).** Nothing is committed yet; ask Paul before committing.
- **Lax.** ℒ_ij = u δ_ij + h E_ij (it was u − h Σ E_ji⊗e_ij). At L=1 the twist is the companion G(z), acting on the right: T^G = ℒ·G. At L=2 all work is in the diagonal frame, with T·diag(z). The factors are ν_n = u − θ + hλ_n, multiplied over sites at L=2. `paper/xxxCG.tex` still uses the OLD conventions.
- **Generators (end of day).** From now on, use the **Hermitian generators EH_ij = R E_ij R⁻¹**, which satisfy EH_ijᵀ = EH_ji. R is diagonal, and M = R² solves M⁻¹E_ijᵀM = E_ji with M₁ = 1. For the adjoint, R = diag(1,1,1,√2,√6,3,2,3) in `ps` order.
- **Layout.** Pre-change notebooks live under `Paul/Mathematica/old-conventions/`; new work goes in `Paul/Mathematica/new-conventions/`.

**Files in `new-conventions/`:**
- **L=1 adjoint.**
  - `Companion-B-C-diagonalize-clean.wb` (`.md`/`.tex`/`.pdf`): E basis, 37 checks.
  - `Companion-B-C-diagonalize-clean-EH.wb`: EH basis, 40 checks, all passing. It adds t_a(G)ᵀ = t_a(Gᵀ), and t_a(diag z) is symmetric.
  - `Companion-B-C-diagonalize.wb`: the working and exploration notebook. It hard-codes h=1 in many cells.
  - `Baxter-QQ-Wronskian.wb`: Baxter, QQ, Wronskian and orthogonality checks for any rep.
- **L=2 adj⊗adj, diagonal frame.**
  - `L2-diagonal-frame-overlap.wb`: E basis.
  - `L2-diagonal-frame-overlap-EH.wb`: EH basis.
  - Both end in a summary table whose only 4 FAILs are deliberate records of approaches that don't work (Section 6, plus one non-equivalent row choice).
- **Abstract algebra.** `Hermitian-structure-abstract.wb`: a noncommutative normal-ordering engine, plus t₁ and t₂ for L sites, with † computed abstractly.
- **Other.** `test.wb` is Paul's connectivity test.
- **Spec.** `docs/superpowers/specs/2026-10-01-L2-diagonal-frame-overlap-design.md`.

**Established at L=1** (adjoint [2,1,0]; checked at (h,χ₃) = (1,1), (1,7/2), (1/2,5/4)). Details are in the clean `.md`. ⟨x₀|y₀⟩ = 9h⁴ in the E basis and 3h⁴ in the EH basis; everything else is basis-independent.
1. **Quantum determinant.** t₃ = χ₃ ν₁(u)ν₂(u−h)ν₃(u−2h).
2. **Baxter equations.** The direct and dual equations keep the old forms but **swap roles**: the new dual equation is the old direct form, and vice versa. Equivalently, each is the old one under h → −h with t_a(u) → t_a(u+(a−1)h). The solutions are Q_i = z_i^{−u/h} q_i with deg q_i = λ₁−μ_i, and Q^j = z_j^{+u/h} q^j with **deg q^j = μ_j − λ₃**. The −λ₃ was found at [2,1,1].
3. **QQ relation.** Q_a(u)Q_b(u+h) − Q_b(u)Q_a(u+h) = (1/z_b−1/z_a) χ₃^{−u/h} F_λ Q^c, with F_λ = ∏_{n=λ₂+1}^{λ₁}(u−θ+nh).
4. **Wronskian (quantum Weyl) formula** for every CBR T_λ in terms of the dual Q's: T_λ = ∏_boxes ν₃(x+(c−r)h) · D_{λ−δ}/D_{−δ}. Checked for all 19 λ ⊆ 3×3.
5. **Reduction at θ−2h.** v(θ) ∥ v(θ−h) and the columns v(θ−3h..θ) span a 2-plane, so T_λ(θ−2h) becomes a ratio of 2×2 determinants, or 0. Only (1), (1,1), (2), (2,1) survive, which are exactly the t's that build ⟨x|.
6. **SoV wavefunctions.**
   - ⟨x_μ|Ψ⟩ = N_μ det₂[v(x₂₁),v(x₁)] v_k(x₂₂)/(same at x₀), with N_μ = n_μ h^{k_μ}.
   - ⟨Ψ|y_μ⟩/⟨Ψ|y_low⟩ = N′_μ det₂[u(y₂₁),u(y₂₂)] u_k(y₁)/(same at y_low), using the direct columns u.
7. **Functional orthogonality from Baxter/dual adjointness.** Columns {−⟨P^jQ_i⟩, ⟨P^jQ_i^{[2]}⟩, ⟨P^{j[−2]}Q_i⟩}, bracket Σ_{n=λ₃}^{λ₁} f(θ−nh)/∏_{m≠n}(m−n)h, rows (1,3),(1,2),(2,3). This is the old relation with h → −h. It also holds for [3,1,0], [3,2,0] and [2,1,1].
8. **Normalized scalar product.** ⟨Ψ_A|Ψ_B⟩ = (1/⟨x₀|y₀⟩) · (det N[ket dual, bra direct]/c₀) / (Ψ_Q^B(x₀) Φ_Q^A(y₀)), with c₀ = 1/(4h⁶) read off by expanding det N, and ⟨x₀|y₀⟩ = 9h⁴. The normalizations are ⟨x₀|Ψ⟩ = 1 (x₀ = e₁) and ⟨Ψ|y₀⟩ = 1 (y₀ = y_{2|21}).
9. **The same formula for other twist pairings.**
   - Between **different twists**, ⟨Ψ_{A,w}|Ψ_{B,z}⟩, with each state carrying its own twisted Q's: all 64 entries match.
   - In the **diagonal frame**, Ω_zw·P_phys|Ψ_{B,z}⟩ = (Ω_zw P_phys·y₀) × (the formula with the w-vacuum as the bra).

**Established at L=2** (adj⊗adj, diagonal twist; one numeric point: h=1/2, θ=(2/7,−3/5), z=(3/2,5/4,2/3), w=(5/2,7/3,3/5)):
1. **Ports unchanged.** The quantum determinant, the Baxter equations with ν_n → ∏_s ν_n, the degrees deg q_i = 2λ₁−μ_i and deg q^j = μ_j−2λ₃, and the QQ relation with F = ∏_s(u−θ_s+2h).
2. **The bra.** Ω_zw⊗Ω_zw is the exact vacuum of t(G′), where G′ = P_aux G(w) P_aux⁻¹.
3. **The overlap formula.** Ω·Ψ_B = (Ω·y₀)/(c₀⟨x₀|y₀⟩²Φ_Q^vac(y₀)) · det N[B,vac]/Ψ_Q^B(x₀) · ⟨x₀|Ψ_B⟩, to 50–51 digits for all 64 states.
   - x₀, y₀, Ψ_Q and Φ_Q are the L=1 objects taken as products over sites.
   - The 6×6 det N has rows (i,j;k), with the L=1 pairs (1,3),(1,2),(2,3) and the bracket label k = 1, 2. Its columns are the three L=1 functionals times x^{β−1}.
4. **Brackets (Paul's preferred choice).** Both are residue sums around ALL poles, with insertions 1 and e^{2πiu/h}. Then c₀ = (ε₁−ε₂)³/[16h¹²θ₁₂³(θ₁₂²−h²)⁴(θ₁₂²−4h²)²], with ε_s = e^{2πiθ_s/h}.
   - Per-ladder brackets also work; there (ε₁−ε₂)³ is replaced by −1.
   - **A single bracket over all poles does NOT work.** It passes the orthogonality test but gives a row-dependent det N. See memory `reference_fsov_per_ladder_brackets.md`.

**Hermitian structure** (`Hermitian-structure-abstract.wb`, abstract and exact):
- **L=1.** t_a† = t̄_a (all scalars conjugated), so for real parameters the t's are Hermitian.
- **L=2.** t₁† − t̄₁ = h̄² Σ_{i≠k}(z̄_k−z̄_i) E⁽¹⁾_ik E⁽²⁾_ki, and

  [t₁(u), t₁(v)†] = −h³θ₁₂ Σ(z_i−z_k)² E⁽¹⁾_ik E⁽²⁾_ki + h⁴[C_z, C_z†], with C_z = Σ z_i E⁽¹⁾_ik E⁽²⁾_ki.
- **When it vanishes.** [C_z, C_z†] has cubic words that nothing else can cancel, so the commutator is zero **only for the untwisted chain**. It is not zero for a unitary twist, for partial equality of the z's, or for imaginary parameters.
- **What holds instead.** t_a(u)ᵀ = t_a of the **site-reversed** chain. This was checked numerically in the EH L=2 notebook, Section 1b. The pairing † ∘ site reversal (via the R-matrix that swaps the sites) is unexplored. Paul said to leave it ("that's fine").

**Traps (each cost time today):**
- **Surds with EH generators.** Exact entries contain surds that Mathematica doesn't combine canonically (√6/2 → √(3/2) next to √6). So `Union[...] === {0}` reports false FAILs. Use `PossibleZeroQ` (`zeroQ` in the EH notebooks).
- **Symbolic det N for c₀.** Zero the Q-symbols outside the target monomial *before* `Det`. The coefficients come out the same, and it is about 100× faster than expanding first, which ran for minutes.
- **Orientation of det N.** Pair the **ket's dual** Q's with the **bra's direct** Q's. At equal twists the other orientation looks fine; at different twists it is off by factors of about 100.
- **`Max[{}]`.** It returns `-Infinity` for an exactly-zero polynomial, and `Ordering` sorts that after numbers. Floor at 0, as `maxc` does.
- **Hard-coded, basis-dependent numbers.** For example ⟨x₀|y₀⟩ = 9h⁴. Use the computed values.
- **Wolfbook.**
  - Insertion order: after inserting k cells after N, the next insert goes after N+k (memory `reference_wolfbook_insert_position.md`).
  - Replies are truncated, so read results with single-line `ToString[..., InputForm]` evaluations.
  - Tool replies may name the wrong target notebook. Verify against the saved `.wb` JSON.
- **Shared kernel.** All new-conventions notebooks share one kernel, and each clean or L=2 notebook clears it on its first cell.

**Next steps (Paul has not set the priority — ask):**
1. **Other reps at L=2.** [3,1,0]² first, plus a second numeric point at adj⊗adj.
2. **General L.** Check whether the bracket structure generalizes to L brackets with insertions e^{2πiku/h}, k = 0..L−1.
3. **The paper.** Update `paper/xxxCG.tex` to the new conventions (E → EH, Lax, twist).
4. **The Casoratian.** Find a closed form for D_{1,0,−1} at L=1 (degree |λ|−3λ₃).
5. **The proposed relabelling** (paused by Paul): direct↔dual, upper↔lower, j → 4−j.
6. **Commit** today's work: new-conventions, the old-conventions move, the spec and these notes.

---

## Everything below predates the 2026-10-01 convention change

It refers to the OLD conventions. For paths, read `Paul/Mathematica/XXX-CG/...` as `Paul/Mathematica/old-conventions/XXX-CG/...`, and `Paul/Mathematica/old/...` as `Paul/Mathematica/old-conventions/old/...`.

**2026-09-25 — thread at the time: one-site localization / SoV bases.** Start at
`Paul/Mathematica/XXX-CG/Building-SoV-bases/HANDOFF-SoV-bases.md` (what was tried, what worked and
failed, traps, next steps). Results are recorded in `Building-SoV-bases/main-notes.tex` (a rolling
record: keep appending) and the computations are in `Building-SoV-bases/L2-localization-tests.wb`.

**2026-09-29 update (best result so far):**
- **Adjoint.** In `Building-SoV-bases/transfer-algebra-dim.wb`, a globally rotated vacuum plus R t11(θ1) gives a one-body covector that reaches the adjoint singlet.
- **[3,2,0].** t22(θ1)·t10(θ1+2h) reaches every irrep of [3,2,0]².
- **Still open (at the time).** [3,1,0] and [4,2,0] had gaps; these are now closed, see the later-session entry above.
- **Where it's written up.** `HANDOFF-SoV-bases.md` ("Latest"), the `main-notes.tex` summary section, and the data in `transfer-algebra-dim-results.m`.

**2026-09-29 later session (gaps closed):**
- **The operator.** At L=1, the SoV-basis operators are T_{μ11,0}T_{λ2,0}⁻¹T_{μ21,μ12}(θ1), with fusion by the CBR formula. At L=2 on the twist-w rotated vacuum, the top one, T_{λ1,0}T_{λ2,0}⁻¹T_{λ2,λ2}(θ1), is one-body.
- **What it reaches.** Every irrep of [3,1,0]², [3,2,0]², [4,2,0]² and [4,1,0]².
- **Single-site interpretation.** It acts as the same operator for site 1 alone, times an explicit prefactor from site 2.
- **Where to look.** `HANDOFF-SoV-bases.md` ("Latest"), `main-notes.tex` and `L2-twistw-vacuum-SoV-ops.wb`.

**2026-09-29 end of day (ratios are polynomial):**
- **The identity.** At θ1, T_{n,0}T_{λ2,0}⁻¹T_{λ2,λ2}(θ1) = T_{n,λ2}(θ1) as an exact operator identity, with coefficient 1.
- **Consequence.** O_max = T_{λ1,λ2}(θ1), and R t11 = T_{2,1}(θ1).
- **Where it holds.** Only at θ1, where the two-site t3 vanishes. It fails at generic u.
- **Caveat on "reach".** It means one covector ⟨v_g|O overlaps at least one copy of every irrep. It is **not** "generates all states", and multiplicity coverage was never recorded.
- **Where to look.** `HANDOFF-SoV-bases.md` ("Latest", item 4) and `main-notes.tex`.

**2026-09-23 reorg:** `Paul/Mathematica/XXX/` and `Paul/Mathematica/Gaudin/` were moved wholesale
to `Paul/Mathematica/old/XXX/` and `Paul/Mathematica/old/Gaudin/` (structure preserved, paths below
updated accordingly). New active XXX work goes in the fresh, empty `Paul/Mathematica/XXX-CG/`
instead of adding to the old tree. See `CLAUDE.md` for the current layout.

## !! General reminder for the next hard session (Paul, 2026-09-24) !!

Full write-up: `Paul/Mathematica/XXX-CG/Q-sampling-count.tex`/`.pdf`.
- **The relation.** All three Q_ab = (1/z_b − 1/z_a) F_λ Q^c, with F_λ = ∏_s∏_{n=λ2+1}^{λ1}(u−θ_s−nh).
  Checked for all states: [2,1,0]² at points 1 and 2, [3,1,0]² and [3,2,0]² at point 1 (`L2-mathcalA.wb`
  cells 12–13). The saved outputs are for [3,2,0].
- **[3,1,0]² count.** 2 ratios of q1 (4 with the determinant) and 1 ratio of q^3. Every reached state has
  deg q1 ≤ 2 and deg q^3 ≤ 1. The missed [4,2,2] and [3,3,2] both have deg q^3 = 2.
- **The q^3 count is the ν3 > λ2 rule.** It has λ2 ratios against deg q^3 = Λ3. The q1 count is not binding:
  [4,4,2] in [3,2,0]² is reached with deg q1 = 2 against 1 ratio.
- **[3,2,0]² count.** The only missed state, [4,3,3], is the only one with q^3 undersampled (a cubic at 3 points).
  Every reached state has q^3 covered. This holds across all three reps checked.
- **[4,2,0]² count.** Same result: the 4 missed states ([6,3,3], [5,4,3]×2, [4,4,4]) are exactly those with deg q^3 > 2
  ratios, and all 15 reached states are covered. So in all four reps, a state is missed ⇔ q^3 is undersampled
  (Λ3 > λ2). The Q_ab relation is NOT checked at [4,2,0] (729 states); F = ∏_s(u−θ_s−3h)(u−θ_s−4h) is conjectured.

**Why ⟨k1^{u/h}Q1⟩⟨k2^{u/h}Q12⟩ cannot give all CG states: count the data (Paul's view).**
- **The data.** The coordinate-space (SoV) CG wavefunction is built from Q-functions. To fully
  reconstruct a highest-weight state you need all its Bethe roots, i.e. the polynomials Q1 and Q12.
  In the twisted chain Q1, Q12 fix every state; in the untwisted chain, only up to its multiplet.
- **What the functionals sample.** The k1/k2 expansions sample Q1 and Q12 only at special points,
  short ladders of shifted inhomogeneities: Q1 at θ1+λ2h,…,θ1+λ1h and Q12 at θ1,…,θ1+λ2h. The
  determinant version adds the matching θ2 ladder for Q1. After normalization, only ratios within
  each ladder count.
- **The shortfall.** When a state has more unknown roots than there are sampled ratios, it cannot
  be reconstructed.
- **Count first.** Before building a new functional, count the points of Q1/Q12 it probes and
  compare with the unknown roots in the table below.
- **Way forward.** We need new generating functions that sample the Q's at more points and localize
  to each site separately (site-2 ladders, other shifts). They should supply the missing data.
- **Count the right unknowns.** Q12 = (1/z2 − 1/z1) ν1(u) Q^3(u), no shift, with
  ν1(u) = ∏_s (u−θ_s−λ1h) and q^3 monic. Checked exactly for all 64 states of the twisted
  [2,1,0]², point 1 (`XXX-CG/L2-mathcalA.wb` cell 12). So the unknowns in Q12 are the Λ3 roots of
  q^3, not deg Q12.
  - Not checked: the untwisted version, expected Q12 = (deg q1 − deg q2) h ν1 q^3.
  - Not checked: other reps. For λ1−λ2 > 1 the factor has degree L(λ1−λ2), probably
    ∏_s ∏_{n=λ2+1}^{λ1} (u−θ_s−nh).

Untwisted degrees, [2,1,0]⊗[2,1,0] (L=2, λ=(2,1,0)). The degrees are staggered:
deg q_i = (4−Λ_i)+(i−1) and deg q^j = Λ_j+(3−j) (`old/XXX/Experiments/SU3_L2_Untwisting.md` §1, `untwU`).
So deg Q12 = deg q1 + deg q2 − 1 = 2+Λ3. That is a trivial degree-2 factor times q^3, with deg q^3 = Λ3.
The unknown roots are then (4−Λ1) + Λ3 = M1 + M2, the Bethe-root count.

| irrep Λ | dim | deg q1 | deg q2 | deg Q12 | deg q^3 | unknown roots |
|---|---|---|---|---|---|---|
| [4,2,0] | 27 | 0 | 3 | 2 | 0 | 0 |
| [3,3,0] | 10 | 1 | 2 | 2 | 0 | 1 |
| [4,1,1] | 10 | 0 | 4 | 3 | 1 | 1 |
| [3,2,1] ×2 | 8 | 1 | 3 | 3 | 1 | 2 |
| [2,2,2] | 1 | 2 | 3 | 4 | 2 | 4 |

Also deg Q13 = deg q1 + deg q3 − 1 = 3+Λ2 = (degree-2 factor) × q^2, deg q^2 = Λ2+1.

The adjoint functionals probe Q1 at θ1+h, θ1+2h (1 ratio) and Q12 at θ1, θ1+h (1 ratio). That makes
2 numbers.
- **Linear cases.** One ratio fixes a monic linear q, since (θ1+h−a)/(θ1−a) determines a. So 2
  numbers cover every irrep with ≤2 unknown roots, i.e. all but the singlet.
- **The singlet.** It needs 4 and is exactly the state we missed. The count gives 5/6, matching
  what we reached.
- **Twisted chain.** The same count works per weight μ: deg q1 = 4−μ1, deg Q12 = 2+μ3.
- **Caveats (Claude).** Distinguishing finitely many Bethe solutions needs less than reconstructing
  them, so the count alone is a heuristic. The singlet overlap is also an exact zero (the ν3 > λ2
  rule), not just a degenerate value. A new functional should be checked against both: enough
  sampled points, and non-vanishing at the missed weights (μ3 > λ2).

## Project goal (Paul, 2026-08-27)

**Generalize the localization/CG construction to su(3) GENERIC representations.** Rectangular reps
are more or less understood. For generic (non-rectangular) reps localization *seems to stop
working*, though some hints remain. Closing that gap is the goal of the project, not a side thread.

gl(3) is the minimal setting where the question can even be posed: every gl(2) irrep is rectangular
up to a determinant twist (`{l1,l2} = {l1-l2,0} (x) det^l2`), so su(2) work can exhibit the
mechanism but never the obstruction. In gl(3), `{p,p,0}` and `{p,0,0}` are rectangular; `{l1,l2,l3}`
with `l1>l2>l3` is not.

Directly relevant prior work: the local-rotations thread (2026-07-18, `Experiments/XXX_CG_L2_V2.wb`
Sections G-L, notes in `Experiments/local_rotations_notes.tex`) already established that
localization fails for general reps, ruled out the spectral-parameter-insertion route, and
identified the **product-covector <-> SoV-basis recoupling matrix `A`/`B`** as the object to
generalize. Existing su(3) notebooks: `Paul/Mathematica/old/XXX/Experiments/su3_V1.wb`,
`su3_V2_general.wb` (arbitrary gl(3) Young diagrams, full-sweep verified), `su3_V2_general_Paul.wb`,
`su3_V2_fund_anti_fund.wb`. `Mathematica/su3/` exists but is empty.

## Current focus

Computing CG coefficients directly in the XXX spin chain (`Paul/Mathematica/old/XXX/`). The Gaudin
subproject proper (residue/representation groundwork) is still parked, but the `h→0` scaling limit
of the XXX construction below is an active connection point, not a revival of that parked track.

**Active thread (2026-09-19): twist-overlap reconstruction, section immediately below.** The
2026-09-17 rotations thread is its direct continuation (per-site rotations are the fix).

**Prior thread (2026-09-17): rotations of the SoV reference covector as determinants of
Q-integrals, gl(3) L=2 adjoint ⊗ fundamental. Self-contained in
`Paul/Mathematica/old/XXX/rotations/` — start at its `README.md`.** The gl(3) adjoint thread
(2026-09-16, below) and the gl(4) `N_untw` thread (2026-09-10) are paused, not closed.

## >>> NEWEST (2026-09-23, late): localizing SoV combinations reach the singlet directions (XXX-CG/) <<<

**Read first:** `Paul/Mathematica/XXX-CG/L2-localizing-transfer.tex`/`.pdf`, Sections 4–5 (Results 5–7). Notebook cells 10–16.
All exact at point 1 of [2,1,0]² (h=1, θ=(1/3,−8/57), z=(5/2,7/3,6/35)), finite twist.

- **Product covectors ⟨v⊗w| at θ1 (cells 10–13).** t_{2,1}(θ1)(Ĝ⁻¹⊗1) commutes with the global gl(3), so t_{2,1}(θ1) = R12 G1 exactly.
  - t_{2,1}(θ1) localizes only for v ∝ w, with w a local hw state of ANY Borel (same on both sites); then v' ∝ wĜ.
  - t_{1,1}(θ1) additionally allows a relative E12 rotation, v = e1ᵀe^{aE12}g with w = e1ᵀg. Its localizing subspace is covariant.
  - t_{1,0} and t_{2,0} localize only on the trivial eigen-covectors e1⊗e1 and e8⊗e8.
  - No t_λ(θ1) localizes a relative E13/E23/Weyl rotation. So mixed Borels such as e1⊗e8 never localize.
- **Minimal localizing SoV combinations (cells 14–16).** Basis ⟨x_{g1|g2}| = ⟨vacL⊗vacL|O_{g1}(θ1)O_{g2}(θ2) (per-site L=1 recipe of `SU3_L1_Baxter_Twisted.wb` cell 11; rank 64).
  - Σc_x⟨x| = ⟨v⊗vacL| is an 8-dim space using only 12 labels.
  - The minimal combinations are 4 singles (100|100, 111|100, 200|100, 211|100) and 4 pairs {g|200, g⁺|100} (g⁺ = n11 raised by one), with coefficients −19/9, −19/14, −95/203, −19/28 (19 from θ12 = 9/19).
  - Their site-1 covectors v·g0⁻¹ have rank 8, so they reach (1,0,2) and (0,1,2), which ⟨k1Q1⟩⟨k2Q12⟩ never reaches. ⟨x_{221|100}| − (19/28)⟨x_{211|200}| hits (0,1,2), the singlet component.
- **Next:** the analytic form of the pair coefficients (scan θ, z, p); a Q-functional realization of the pairs (a ⟨k1Q1⟩-type bracket with site-2 poles?); site-2 localization; other reps.
- **Same day, `XXX-CG/xxx-210-sov-limit.wb`** (untwisting limit of the SoV basis, [2,1,0]², z = 1 + φε):
  - Pointwise, all 64 SoV covectors → ∝ vacL⊗vacL (rank 1). The untwisted t's act as scalars on the multiplicity-one [4,2,0].
  - After undoing g0⊗g0 and projecting to total-weight sector ν, every nonzero component starts exactly at ε^{ht(ν)} (ht = simple-root steps below (4,2,0)). The graded limit has rank 64, so any (e1ᵀP)⊗e1ᵀ with P = exp(εaE12 + ε²bE13 + εcE23) is one ε-independent combination.
  - Not checked: that Ψ diverges with the same grading.
  - Cells 5/6 there (the "limit of span" hw-count) are frame-dependent; don't reuse them.
- `xxx-210-cg-test.wb` cell 16: explicit 8-component v_new (finite twist, factored table).
- **Traps:**
  - A multi-line `f[...] := a . b .` ending a line with `.` did not register in wolfbook (use one line).
  - Symbolic 64×64 LinearSolve / `Together` on symbolic inverses hang the kernel, and abort doesn't release it (restart). Work per weight sector (≤10-dim) or with truncated ε-series.

## >>> ACTIVE (2026-09-23): 𝒜 functionals ⟨k1^{u/h}Q1⟩, ⟨k2^{u/h}Q12⟩ at L=1 and L=2 (XXX-CG/) <<<

**Read first:** `Paul/Mathematica/XXX-CG/L2-mathcalA.tex`/`.pdf` (L=2), `L1-mathcalA.tex`/`.pdf` (L=1).
Notebooks `L1-mathcalA.wb`, `L1-orthogonality-diagonal-frame.wb`, `L2-mathcalA.wb` (switches `repT`,
`ptIndexT` in cell 3; ledger `resultsT` in cell 11). Normalization ⟨vacL⊗vacL|Ψ⟩ = 1 = ⟨Ψ|vacR⊗vacR⟩.

L=2, same rep on both sites, checked for [2,1,0], [3,1,0], [3,2,0] at 3 points (≤1e-77):
- 6×6 twisted FSoV det / (product of L=1 c^ref Ψ^ref) = ⟨Ψ_A|Ψ_B⟩ with one constant
  κ2 = (−1)^{λ1+1} h^{6λ1} κ1² [∏_{n=−λ1}^{λ1}(θ12−nh)]³/(e1−e2)³ (exact at points; sign/h-power fitted on λ1=2,3).
- ⟨k1^{u/h}Q1⟩ with the L=1 measure on site 1 only = ⟨v1 ⊗ vacL|Ψ⟩, v1 = the L=1 covector.
- ⟨k2^{u/h}Q12⟩ with site-1 poles only and measure 1/[∏_{k=0}^{λ2}ν3^{(1)}(u−kh) ∏_{n=λ2+1}^{λ1}(u−θ2−nh)]:
  product of both = ⟨v1⊗vacL| Σ_m (−1)^m C(λ2,m) k2^m t^{(1)}_{m,m}(θ1)⊗1 |Ψ⟩ — the L=1 formula exactly.
  Without the site-2 factor it holds only up to θ-dependent c_m; the full ν3(u)ν3(u−h) does not fix it.

**Which transfer matrices localize on a global rotation (2026-09-23)** — `XXX-CG/L2-localizing-transfer.wb`, [2,1,0]², exact,
T = L1 L2 G. ⟨u⊗u|t_{2,1}(θ1) = τ^vac_{2,1}(θ1) ⟨u^G⊗u|, u^G = e1ᵀG⁻¹gG (a,b,c → a z2/z1, b z3/z1, c z3/z2), symbolic in a,b,c:
site 1 only (the regularity picture, t(θ1) = R12 G1). Scan over θ_s+nh, n = −5..5: t_{1,1}(θ1) is also site-1-only, but its w is not a
rotation; t_{1,1}(θ2+3h) acts on site 2 only; t_{2,2}(θ_s±h) kills ⟨u⊗u|; no t_{2,1}(θ2+nh) is site-2-only (ordering). Formula re-checked at point 2 with u = vacL (cell 7).
CG count for ⟨vacL⊗vacL|t_{2,1}(θ1) (`xxx-210-cg-test.wb` cells 11–13): after undoing g0⊗g0, v_t = e1ᵀ G⁻¹g0Gg0⁻¹ = e1ᵀe^{yE23}e^{xE12}
(zero (1,3) entry from the p3 relation), so it is ≡ 0 at (1,0,2),(0,1,2) like v_new. Alone 4/6 hw, combined with old + det still 5/6: nothing new, singlet missed.
**Why the ν3 > λ2 rule — now derived (cells 14–15; tex sections in xxx-210 and a remark in xxx-310):** site-1 μ3 = ν3, and two mechanisms cap it:
(a) t̂ families (old, det, w-vac): rotation part e1ᵀe^{aE12} or e1ᵀ n…n⁻¹ has μ3 = 0, and each quantum minor moves ≤1 unit into slot 3 FOR ANY
TWIST (checked exactly with a general 3×3 twist), with ≤ λ2 insertions. The (1,3) = 0 of G1 is NOT the cause here. (b) t_{2,1} family: e1ᵀe^{yE23}e^{xE12},
zero (1,3) entry, and e1ᵀE23^k = 0 for k > λ2. The "everything else is reached" half of the rule is still observed only. A (1,3)
entry opens E13: hypothetical v^(2) = e1ᵀG⁻²g0G²g0⁻¹ has (1,3) = −p2(z1−z2)(z1−z3)/z1² = O(ε²), hits the singlet, all families + v^(2) = 6/6.
NOT realized by t_{2,1}(θ1)²: that is rank 4, not local (L2-localizing-transfer cell 8). Matches the 2026-09-19 finding (one per-site E31 rotation completes).

**6×6 vs vacuum at twist w (2026-09-23)** — `XXX-CG/L2-vacw-mathcalA.wb` (cell 1 = summary). The 6×6 FSoV with left Q^j
of Ω=e1⊗e1 at twist w is C⟨u⊗u|Ψ_B⟩, a GLOBAL rotation; flag of u: ℓ_k ∝ ℓ0_k (z_k−w2)(z_k−w3), π_k̂ ∝ π0_k̂/(z_k−w3)
((ℓ0,π0) = vacL flag). [2,1,0] pts 1,3 and [3,1,0] pt 1, 3 w each, ≤1e-83; θ-, h-, rep-independent. The plain product
D_B(k1)·S_B(w) is NOT a rotation (σ2/σ1 ~ 0.1–0.45). u(w) is NOT in the vacL(z) family (violates the p3 relation)
but IS a vacL for the Möbius-transformed twist y_k = (αz_k+β)/(z_k−w2) (exact, cell 12).
S_B(w) × the unchanged ⟨k2^{u/h}p Q12⟩ = ⟨V(w)|Σ_m(−1)^m C(λ2,m) k2^m t^{(1)}_{m,m}(θ1)⊗1|Ψ_B⟩: site 1 only again, c'_m = ρ_m,
w-independent (cells 13–14; [3,1,0] pt1, [3,2,0] pt1, [2,1,0] pt3, 3 w each, ≤1e-76). 6×6 flag formula also holds for [3,2,0].
CG count for it (`xxx-310-cg-test.wb` cells 9–10, z and w both untwisted): v''_new = e1ᵀ g_u T g_u⁻¹ reaches 3/6, 3/10, 6/10 —
identical to the det family (same zero pattern, both e1ᵀ n T n⁻¹, n ∈ N+); all three families combined still 5/6, 8/10, 9/10. Asked Paul how the two should be combined (not tried: a joint determinant).

**[4,2,0] check of the k1 det + k2 factor (2026-09-23)** — `L2-det-mathcalA.wb` at [4,2,0] pt1 (dim 729): D = C⟨v1⊗v1|Ψ⟩ to
1e-85 with C the closed form (prediction); c'_m = ρ_m = (1, 67/29, 1608/145), final residual 1e-74. Ledger rows 1–7 re-seeded.
CG count [4,2,0]² (`xxx-310-cg-test.wb`, saved outputs = [4,2,0]): 19 hw; old 15, det 6, w-vac 6, all 15 (finite twist also 15).
Missed [6,3,3], [5,4,3]×2, [4,4,4] = exactly ν3 > λ2 (predicted); v_new ≡ 0 exactly at μ3 ∈ {3,4}. Rule now holds on 4 reps.

**Docs (2026-09-23):** `XXX-CG/L2-det-mathcalA.tex/.pdf` (determinant + second factor), `xxx-210-cg-test.tex/.pdf`
(now with "determinant family", "t_{2,1} family" and "why ν3 ≤ λ2" sections), `xxx-310-cg-test.tex/.pdf` ([3,1,0]², [3,2,0]², the ν3 > λ2
rule + its derivation, reconstruction), `L2-localizing-transfer.tex/.pdf` (one-site action of t_{2,1}(θ1), scan, non-iteration).
Kernel state at end of 2026-09-23: xxx-210 cells 2–14 run for [2,1,0] (this overwrote the [4,2,0] C-suffix variables of xxx-310).
Kernel was restarted after these runs — rerun notebook cells before relying on kernel variables (resultsD, resultsF, jr*).

**2×2 determinant, same rotation on both sites (2026-09-23)** — `XXX-CG/L2-det-mathcalA.wb` (cell 1 = summary,
cell 9 = ledger `resultsD`). det[[⟨k1^{u/h}Q1⟩_α, ⟨k1^{u/h}uQ1⟩_α]]_{α=1,2}, ⟨f⟩_α = ∮ e^{2πi(α−1)u/h} f / ∏_{k=0}^{d} ν1(u+kh)
(two-site ν1, d = λ1−λ2), divided by (e2−e1) Q1(θ1+λ2h) Q1(θ2+λ2h) k1^{(θ1+θ2)/h+2λ2}, equals
C ⟨v1⊗v1|Ψ⟩ with C = (−1)^d / ((d!)² h^{2d} ∏_{n=−d}^{d}(θ12−nh)). 7 runs ([2,1,0] pts 1–3, [3,1,0] pts 1–2, [3,2,0] pts 1,3),
all orders in k1, ≤1e-85; C fitted on 4 runs, predicted on 3. Normalization is at the lowest pole θs+λ2h, not Q1(θs).
Times the unchanged L2-mathcalA second factor ⟨k2^{u/h}p(u)Q12⟩ (site-1 poles, p = site-2 factor): = C⟨v1⊗v1|Σ_m(−1)^m C(λ2,m) k2^m
t^{(1)}_{m,m}(θ1)⊗1|Ψ⟩ — k2 still acts on site 1 only; c'_m with v1⊗v1 = ρ_m, same as with v1⊗vacL (cells 10–11;
[2,1,0] pt1, [3,1,0] pt2, [3,2,0] pt3, ≤1e-80). Open: symmetric version (k2 on site 2 / both sites).
CG count for it (`xxx-210-cg-test.wb` cells 8–10, exact): undo g1 (site-2 rotation of v1, ∈ N+, fixes hw states) → v'_new =
(1−k1z1)² e1ᵀ(1−k2 g1 t̂ g1⁻¹), zero at (1,2,0),(0,2,1),(1,0,2),(0,1,2); reaches 3/6 hw states, combined with old 5/6 —
**singlet still missed**, ⟨v_tot⊗v1|singlet⟩ = 0 exactly. Hypothetical fix: t̂ on BOTH sites, ⟨v1t̂⊗v1t̂|singlet⟩ =
½p2²(z1−z2)²(z1−z3)² ≠ 0 → needs a site-2-pole ⟨k2 p Q12⟩ factor (not built yet; Paul: hold off on proposals).
[3,1,0]² (`xxx-310-cg-test.wb`, switch `lamSwitchC`, reproduces adjoint): 10 hw states; old 8/10, det 3/10, combined 8/10.
Missed [4,2,2], [3,3,2] (site-1 weights (1,1,2), (0,2,2)). v_new ≡ 0 at every site-1 weight with μ3 ≥ 2, as for the adjoint →
in both reps the missed irreps are exactly those with ν3 = 2.
[3,2,0]² (same notebook, cell 3 now has t11 and t22; saved outputs = [3,2,0]): 10 hw; old 9, det 6, combined 9. Missed only
[4,3,3] (site-1 weight (1,1,3)). v_new ≡ 0 exactly at μ3 = 3. Rule over all three reps: v_new ≡ 0 iff μ3 > λ2; missed irreps
= those with ν3 > λ2; det family never adds a state.

**CG test, adjoint (2026-09-23)** — `XXX-CG/xxx-210-cg-test.wb` + `.tex`/`.pdf`, all exact/symbolic.
Undo vacL globally: (v_tot ⊗ vacL)(g0⁻¹⊗g0⁻¹) = v_new ⊗ e1ᵀ, v_new = v_tot g0⁻¹. Untwist z_k = e^{φ_k ε},
leading order per weight space of V. Reached: 6/8 site-1 weights; **5/6 highest-weight states of V⊗V** —
[4,2,0], [3,3,0], [4,1,1], both [3,2,1]; **singlet [2,2,2] missed** (needs site-1 weight (0,1,2), where v_new ≡ 0
even at finite twist; (1,0,2) also ≡ 0). k1, k2 survive the limit, θ, h drop out; in the octet sector the
direction is (2k1−1, 1), so k1 alone separates the two octets. Lower weights come only from the rotated twist
G1 = g0 G g0⁻¹ (upper triangular, O(ε) off-diagonal, (1,3) entry = 0 by the p3 relation).

**Next:** reach the singlet (second t̂ insertion / λ2 ≥ 2 / other functionals); repeat the count for
[3,1,0]², [3,2,0]²; λ2=3 (needs t_{3,3}); different reps per site; derive c_m/κ2 rather than fit.

## >>> PREVIOUS (2026-09-19): twist-overlap reconstruction — the root-step rule <<<

**Read first:** `Paul/Mathematica/old/XXX/Experiments/SU3_Adjoint2_TwistReconstruction_L2.md` (full
reference), notebook `SU3_Adjoint2_TwistReconstruction_L2.wb` (31 cells, saved).

**Question (Paul).** Companion frame, L=2. Knowing `<Om(w)|tau_a(y,u)|Psi(z)>` for all vacuum
twists `w`, all transfer-matrix twists `y`, all `u`, `a = 0,1,2`, how much of `|Psi(z)>` is fixed?

**Answer.** Psi never enters: the solution set is `|Psi_A(z)> + W^perp` with
`W = span{<Om(w)|tau_a(y,u)}`, so the whole question is one rank, the same for all states.

```
W = V_{lamA+lamB}  (+)  everything ONE ROOT STEP away,  with multiplicity
```

| reps | dim | dim W | dim W^perp |
|---|---|---|---|
| `[2,1,0](x)[1,0,0]` | 24 | 24 | 0 |
| `[2,1,0]^2` | 64 | 63 | 1 (the singlet) |
| `[3,0,0]^2` | 100 | 63 | 37 |
| `[3,1,0](x)[2,1,0]` | 120 | 111 | 9 |
| `[3,1,0]^2` | 225 | 171 | 54 |

**Mechanism (derived, cells 29/31).** `qmin_{I,K}` with `|I|=|K|=a` transforms in
`Lam^a(C^3) (x) Lam^a(C^3)*`, which is `adj (+) 1` for a=1 AND a=2 (the `det` factors cancel), and
`1` for a=3 — hence `tau_3 ∝ Id`. Measured: the module generated under `ad` closes at dim 9 with
exactly 1 invariant, the invariant being the untwisted `Sum_i T_ii(u)`. And `<Om(w)|` sits in the
single component `V_{lamA+lamB}`. Wigner-Eckart then gives one root step per insertion.

**Fixes tested.** (i) `k` insertions reach `k` steps — measured 63/90/100 on `[3,0,0]^2`, so k=3
completes it; cost grows linearly with the rep. (ii) A general non-companion twist: NO help (the
module is twist-independent). (iii) The other five polarized vacua: NO help (all are
Cartan-component covectors). (iv) **One per-site rotation `<Om(w)|(exp(c E^(1)_31) (x) 1)`: ranks
63, 64, 64, 64 for 1,2,3,4 values of `c` — ONE nonzero `c` completes `adj^2`.** This is O(1) extra
data, and `rotations/` section 19 showed per-site rotations ARE Q-determinants, so the overlaps are
computable. **This is the bridge to the rotations thread.**

**NOT** the singlet story (`[3,1,0]^2` has no singlet yet loses 54) and **NOT** the
rectangular-vs-generic obstruction (`[3,0,0]` is rectangular and loses 37).

**Next, in order.** (1) Push the one-rotation-parameter fix to `[3,1,0]^2` and `[3,0,0]^2` — does
one `c` still suffice? (2) Does the *cheap* `E_32` rotation suffice, or is `E_21`/`E_31` needed
(sections 17/19 of rotations say those are hard as Q-determinants)? (3) Scaling of the seed ladder
for `[S,S,0]^2` vs SoV's `dim V`. Full list in section 11 of the `.md`.

**Traps.** `Join @@ Table[...,{a},{b}]` does not stack rows (bit 4x); hoist `om.O` out of inner
loops and use `prim` before `MatrixRank` or exact-rational rank will not converge at dim 100+;
`C_2` is NOT injective on partitions of fixed size (`(3,3,0)` and `(4,1,1)` both give 24 at
`|lam|=6`), so `lamOfC2` can silently mislabel.

## >>> CONTINUATION (2026-09-17): which rotations ⟨0|g are Q-determinants? <<<

**Read first:** `Paul/Mathematica/old/XXX/rotations/README.md`, then the `.tex`/`.pdf` there, then the
`.md` for exact numbers and the list of things that already failed.

**Where it stands.** `⟨0| = e₂₄`; the SoV covectors are fixed by `⟨x|Ψ⟩ = Ψ(x)/Ψ(x₀)` (all 24
states, 1e-73). The FSoV scalar product equals the physical one with `e₂₄` on both sides (2e-78).
Pairing it with the vacuum dual Q's at a second twist `w` gives the rotations
`(a,b,c) = (-w₁w₂, w₁, -w₁(w₂+w₃))`, hence the **surface `b = a(c-a)`**, which nothing deforms off
(225 minors, `e^{ku}` dressing, six permuted vacua). Small determinants with a pole-cancelling
polynomial give the `E₃₂` rotations — one site (local) or both (global, `c = -e^{kh}`, measure
`1/(ν₁ν₂) = (u-θ₁)×ladder`), and these are minors of the same scalar-product matrix.

**New (2026-09-18, cells 48-54).** Two of the README's open items are now closed.

- **§16, products of one-site brackets.** The per-site rotation `(e₈e^{c₁E₃₂})⊗(e₃e^{c₂E₃₂})` has
  SoV support on exactly the four `x₁₁` labels with the exact pattern
  `[[1,(47/28)c₂],[c₁-(19/28)c₂,c₁c₂]]`, `det = -19c₂(28c₁-47c₂)/784`. A product of two one-site
  brackets is rank one there, so it is a per-site rotation only on `c₂=0` (= §8) or
  `c₁=((t-2h)/(t-h))c₂` — **independent `(c₁,c₂)` is impossible**, the old conjecture was wrong.
  With `e^{k_su}` dressings `r_s=-e^{k_sh}` and the locus is `k₂-k₁=log(47/28)`. Also: the one-site
  bracket on the **fundamental** site is not a rotation at all (and not a product).
- **§17, Wronskian minors — two no-gos.** (i) Exact/Gröbner: an object built from one site's
  `Q`-values has `c_x` in the 8 site-2-base labels, and the only rotations supported there are the
  local `E₃₂` (`⟨b,a⟩`); no global rotation at all (`⟨c²,b,19a+9c⟩`). So **no local `E₂₁`/`E₃₁` from a
  one-site object**, whatever the measure. (ii) Numerical/complete: of the 9315 minors of the FSoV
  matrix with at least one `i=2` row (sizes 2–5), **none** is an SoV expansion — the best two-twist
  mismatch is `O(10^-1)`. Only the full 6×6 works, as the degree count (`Ψ_A` is degree 6) predicts.

- **§18, free row weights — the surface is rigid.** Writing the 6×6 with free row weights (30
  parameters, base = FSoV measure) and demanding both twist-independence (23 conditions) and
  rotation-ness (`δV ∈ span{V₀ρ(X)}`, 20 conditions) leaves a **14-dimensional** allowed space — twice
  the 7 obvious directions (6 row rescalings + `e^{ku}` dressing), which sit inside it. On all 14,
  `δ[b-a(c-a)] = 0` to 57–65 digits, and the induced `(δa,δb,δc)` has rank exactly **2** = dim of the
  surface. Free weights move the rotation *along* `b = a(c-a)`, reaching the same 2-parameter family as
  varying `w`, never off it. The response is exact, not finite-difference: a determinant is
  multilinear in its rows, so replacing row m by a single-pole row **is** `∂D/∂v^{(m)}_p`.

- **§19, independent per-site rotations — YES, they are Q-determinants.** §16's constraint
  `c₁/c₂ = (t-2h)/(t-h)` was a property of the *product* of two one-site brackets (rank-one coefficient
  matrix), not of bracket objects. A `2×2` **determinant** with free per-site weights `R,S` and a free
  insertion `g(u)` gives `c₁₂ = Sδ`, `c₂₁ = R(1-γ)`, `c₂₂ = RS(δ-γ)` — four knobs against three
  conditions — and hits **every** target `(c₁,c₂)` exactly: `(1,2)`, `(2,½)`, `(-3/7,5/2)`, `(1,0)`,
  `(0,1)`, all twist-independent, `σ₂/σ₁ ~ 1e-94` (genuine products), per-site fit residuals `~1e-94`.
  Even `1 ⊗ e^{c₂E₃₂}`, impossible for a product, works (`R = 19/28`, `S = 47/56`, cubic insertion).
  `γ` stays free: a one-parameter family of realisations per target. Note per-site rotations are **not**
  of the form `⟨0|g`, so they are neither on nor off the surface `b = a(c-a)` — that surface constrains
  global rotations only.

- **§20, Paul's product (global 2x2 det at `k`) x (one-site bracket at `k'`) — nearly right, fails by
  one power of `Q₁`.** Checked on all 24 states. With `X₁ = Q₁(θ₁+2h)/Q₁(θ₁+h)`,
  `X₂ = Q₁(θ₂+h)/Q₁(θ₂)`, the object is *exactly* (72-digit fit residual)
  `D = G_c(1 + c'X₁)`, `G_c = 1 + (9c/28)X₁ + (47c/28)X₂ + c²X₁X₂`, `c = -e^{kh}`, `c' = -e^{k'h}/z₁`,
  while `⟨0|g₁g₂` is `P = 1 + (9c/28 + c')X₁ + (47c/28)X₂ + c(c+c')X₁X₂`. So
  `D - P = (cc'/28)·X₁·(9X₁ + 19X₂ + 28c X₁X₂)`: the constant, `X₁`, `X₂` coefficients are exactly
  right, `X₁X₂` overshoots by `(19/28)cc'`, and `X₁²`, `X₁²X₂` are spurious. `X₁²` is not a
  wave-function ratio — the product puts `Q₁` at a site-1 `x₁₁` point twice (determinant + bracket)
  while `Ψ_A` has one `Q₁(x₁₁)` per site, and `1/(Q₁(θ₁+h)²Q₁(θ₂))` cannot cure that. Measured:
  rank`{c(z),c(z')} = 2` (not an SoV expansion), 14 labels, stabiliser 1, `O(1)` deviation after the
  best constant. **Exact correction:** `⟨0|g₁g₂ ↔ G_c + c'X₁(1 + cX₂)` — the dressing must hit only the
  site-2 leg. That is a sum of two products, i.e. an antisymmetric combination, which is what §19's
  determinant-with-insertion produces.

**Kernel hygiene — read this before running anything.** The kernel inherited at the start of
2026-09-18 held caches from two epochs (the previous session had restarted mid-way), so
`vacAtTwist[w]["qd"]` was inconsistent with the cached §7–§11 results: re-running the §10 object gave a
covector 0.1% off `⟨Ω^L_w|` at 71-digit precision, and recomputing `vacAtTwist` failed outright. After
a restart, **all of cells 1–57 were cold-run in order with no kernel messages**, and every stored
result reproduced exactly (§10 gives `(a,b,c) = (-35/6, 5/2, -263/42)` to 10 digits, residual 1e-70;
§16/§17 exact results bit-identical). Two minor counts in the §17 survey changed (nonzero minors
880 not 884 for shape (2,1), 64 not 73 for (2,2)); the docs carry the cold-kernel numbers. Do not
trust a warm kernel for a base point.

**Next.** Paul's stated want is *different rotations on each site*; §19 delivers that for `E₃₂`. The
natural continuation, and now the most promising item:

1. **Per-site `E₂₁`/`E₃₁`.** Free weights *plus* a free insertion is what unlocked §19; that
   combination has never been tried on an object with `Q₂` rows, which is what moving `x₂₁`, `x₂₂`
   needs. §17's no-gos do **not** cover it (they kill one-site objects and all FSoV minors, not a
   two-site free-weight Wronskian object). Start here.
2. **Finite** weight deformations of the `4×4` shape (3 rows `i=1`, 1 row `i=2`) against the explicit
   degree-4 `E₂₁` target in the §17 result cell — 16 parameters vs 46 conditions, and §18 shows
   first-order rigidity of the global surface, so the structure must conspire. Ask before grinding.
3. Vary the dual state *and* the weights together (each was killed separately, never jointly).
4. A different reference covector, i.e. a different `B`/SoV basis.
5. Derive `b = a(c-a)`; §18 makes this more pressing, since the relation is a property of the
   construction rather than of the measure.

## >>> PAUSED (2026-09-16): gl(3) adjoint L=2 — twisted match + untwisting <<<

**Read first:** `Paul/Mathematica/old/XXX/Experiments/SU3_L2_Untwisting.tex` (5-page PDF, the full
story) and the two companion `.md` files below. `summary.md` at the repo root is the *previous*
session's handoff (superseded by this section).

### Notebooks (Experiments/, all saved, none cold-run end-to-end in one go)

- `SU3_L2_Baxter_Twisted_Final.wb` + `.md`
  - **Twisted L=2 match.** For a unimodular twist `z=(5/2,7/3,6/35)` and `θ=(1/3,−8/57)`, the
    normalized physical scalar product `⟨L_A|R_B⟩/(⟨L_A|0⟩⟨0|R_B⟩)` equals
    `C_vac det N^{AB}/(c^ref_A Ψ^ref_B)` for **all 64×64 pairs, both a/b branches** (residual 1e-76).
    Cells 1–8 were cold-run.
  - **Rows.** `(α,i,j)`, α∈{1,2} (bare / `e^{2πiu/h}`), each bracket one contour around all 6 poles of
    `1/(ν1ν2ν3)`. Columns: `u^β×(−Q^jQ_i, Q^jQ_i^{[−2]}, Q^{j[2]}Q_i)`.
  - **Normalization.** `c^ref` = L=1 coefficient multiplied over θ1, θ2.
  - **Invariance.** All 36 row triples `{(a,c),(a,d),(b,c)}` work, with the normalization following the
    Q-content (`Q_a(θ+1)W(Q_a,Q_b)` and its dual) and a **triple-independent `C_vac`** (cells 9–10).
  - **Bug fixed.** `bra02`/`ket02` needed `Flatten`.
- `SU3_L2_Untwisting.wb` + `.md` + `.tex/.pdf`, with `z_j = e^{εφ_j}`, 200 digits:
  1. **Divergence.** `⟨0|R⟩ ~ ⟨L|0⟩ ~ ε^k`, `k = ht(λ_vac−λ_A) = 0,1,1,2,2,4`.
  2. **Precise relation (Appendix N of 2603.26617 generalized).** In the column basis H (3 surviving
     integrals of motion) / O (3 `O(φ)`), the bare rows are at infinity, and
     `det N_reg = −det[A1]_O · det[A2]_H`. The first factor is `V_λ(z)·E_λ`, with
     `E_λ = −(λ1−λ2+1)(λ2−λ3+1)/[(λ1−λ3+1)!(λ1−λ2)!(λ2−λ3)!]`. The second is the **untwisted**
     periodic 3×3 minor.
  3. **Result.** `lim V_{λvac}/V_{λA} · S(A,A) = B_{λA}/B_{λvac}`, with `B_λ = E_λ det𝒩/(c^ref Ψ^ref)`
     entirely in untwisted q's. `𝒩` is written out explicitly in the tex.
     - **Evidence:** matches the physical limit to 1e-19 at θ P1; `E_λ` independently validated at
       θ=(1/5,46/55) to 1e-13; **not derived** analytically.
  4. **Bare eigenvectors.**
     - `l_0 = (S⊗S) r_0`, `S = diag(1,1,1,2,6,9,4,9)`.
     - With `r_0` fixed by its |hw⟩⊗|extremal⟩ component, `B/B_vac = ‖r_0‖²_S/(N_R N_L)`, with
       `N_R N_L = 1, (1−t²)/4, (1−t²)/4, [(1−t²)(9−t²)/32]²` for 27, 10, 10̄, 1 (`t = θ1−θ2`).
     - Octets: the `B` formula is known in closed form (with √(4+5t²)), but no natural normalization
       has been found yet.
  5. **Separate limits.** `⟨0|R⟩/ε^k` and `⟨L|0⟩/ε^k` are each monomials in `φ_ij`. So
     `|R_untw⟩ = lim ∏(z_i−z_j)^{aR}|R_A⟩ = r_0/b_R` and `⟨L_untw| = lim ∏(z_i−z_j)^{aL}⟨L_A|` exist,
     and the object in (3) equals `⟨L_untw|R_untw⟩`. The split and `b_R`, `b_L` probably depend on
     a13, a23 and the branch.

### Next steps (Paul has not set the priority — ask)
1. **Derive `E_λ`** from the residue at infinity. `O1+O2` cancels at leading order, so the next order is needed.
2. **Octets.** Is `N_R N_L` an eigenvalue / matrix element of the adjoint⊗adjoint `Ř(θ1−θ2)`? That
   would also fix the octet normalization, since the octet block carries √(4+5t²).
3. **`b_R`, `b_L`:** closed form and a13/a23/branch dependence.
4. **Cold run** of `SU3_L2_Untwisting.wb`; then promotion to `Clean/`.

### Traps
- **Kernel state.** `SU3_L2_Untwisting.wb` cells 13/16/17/19 rebuild the chain at other θ.
  `buildAtU` does not rebuild `polesU`/`weightsU`. Run in order.
- **`0^0`.** `qpolyU` must be Horner, because a pole at u = 0 gives `0^0`.
- **Output truncation.** Wolfbook output truncates at `{`.
- **Toy model.** The "untwisted q × exact exponentials" toy gives the wrong coefficient; don't use it.
- **Large-denominator `Rationalize`** results for octets/singlet are unreliable.

## >>> PAUSED: closed form for N_untw — PARTIAL (gl(4), 2026-09-10) <<<

**Full briefing: `Paul/Mathematica/old/XXX/Experiments/TASK_N_untw_closed_form.md`** — read it in full.
Work in `Experiments/SU4_SoV_PerturbativeOverlap_L2.wb`, cells 19-23 (saved). Kernel restart: run
cells 1-15, then 19-22.

**Done.** `N_untw` = determinant of the **untwisted** 8x8 built from the untwisted Baxter
polynomials (`spMatU`/`detU`, cell 20; needed a new untwisted *dual* solver `UntwSolveD`). Verified
exact for `mu = {m,m,m,m}`, where that determinant is nonzero. For the vacua every Q is a pure
exponential, so the matrix is closed-form (`spMatE`, cell 21).

**Empirical, not derived.** For all other `mu` the untwisted determinant vanishes. What we have is an
exact ratio in the inhomogeneities only (`t = (theta1-theta2)/h`, `d = deg Q_12`):

```
N_untw(mu)/N_untw(vac) = c(S,P,d) [(S+t)(P-t)]^3 PROD_{k=1}^{d-1} (S+k+t)(P+k-t)
```

All six ratios come out **exactly rational**, verified at two different `theta2` to 30+ digits. But
this is **not** the Q-function answer Paul asked for: the rewriting via `Qbar_12(theta1+Sh)
Qbar_12(theta2+Ph)` leaves a non-Q factor `[(S+t)(P-t)]^2` and an unexplained `(S,P,d)`-dependent
constant, and `c(S,P,d)` is a fitted table.

**Ruled out:** the whole first-order reduced-minor / Appendix-N (N4) class — `corank(spMatU)` is
strictly less than the vanishing order for every degenerate state, so the degeneracy is not lifted at
first order. Also ruled out: pure Q-products (exhaustively), and both naive exponential dressings.

**Best lead:** for all six states the monic branch polynomials have degree `<= 1` (`q = 1` or `u+b`),
so the twisted Baxter system is **exactly solvable in the twists** — started and works (branch 1 with
`q=1` fixes `tc[3..5]` in closed form in `z_1..z_4`). That makes the confluent limit an analytic
Taylor extraction instead of a numerical fit. Section 5 of the task file has the route.

## Perturbative twist expansion — DONE (2026-09-10)

`Paul/Mathematica/old/XXX/Experiments/SU4_PerturbativeExpansion_L2.wb` (clean, self-contained, cold run
a few seconds) + `SU4_PerturbativeExpansion_L2.md` (reference + the bug list). Scratch predecessor
with the dead ends: `SU4_PerturbativeBaxter_L2.wb`.

**Result.** Two branches alone determine the whole transfer matrix, with a one-order lag: pushing to
`eps^N` with only `Q_1, Q_2` fixes all five `tc` at orders `1..N-1` (verified to N=8 against the
exact pseudovacuum `tau`); only the top order keeps free directions. The dual side `Q^4, Q^3` does
the same, identically.

**Structure.** `Q_k = z_k^(u/h) q_k`, `q_k` monic of FIXED degree `M_k`, exponential kept exact and
expanded by the operator, coefficients Laurent with branch `k` diverging as `eps^-(M_k-M_1)`. Dual:
`Q^j = Exp[-phi_j u] q^j`, `deg q^j = mu_j`. The staggered `Mhat_k = M_k+(k-1)` are the shadow of
this — the untwisted `P_k` are limits of COMBINATIONS, not of individual branches. For `k=2` the
regular combination is the divided difference: `b = -a`, `a ~ 1/eps`, nothing worse.

**Read the .md before touching this.** The expensive bug: accumulated solution rules must be
BACK-SUBSTITUTED, `acc = Join[acc /. new, new]`, or a later order falsely reports INCONSISTENT.
Also: `Solve` treats symbols outside `vars` as parameters and returns `{}` on generic failure, so
undetermined lower-order coefficients must go back into the unknown list.

## NEXT TASK (2026-09-10, Paul): leading small-twist behaviour of the scalar product

**Goal.** Feed the twist expansion we already have into the SoV scalar product and compute its
leading behaviour as the twist is switched off. The hoped-for form of the answer: the result
written in terms of **level-0 (untwisted) Q-functions**, multiplied by overall twist factors whose
powers are fixed by the **quantum numbers of the state** (i.e. by `mu`, equivalently by the degrees
`M_k = (S+P) - mu_k`).

**Start from these two files.**

- `Paul/Mathematica/old/XXX/Clean/SU4_SoV_ScalarProduct_L2.wb` — the finished scalar product, 9 cells,
  self-contained, runs from a cold kernel (~48 s for the orthogonality table). Companion
  `.md` (symbol reference) and `.tex`/`.pdf` (derivation, 6 pp). READ THE `.tex` FIRST: it is much
  cheaper than reconstructing the argument, and it lists the pitfalls.
- `Paul/Mathematica/old/XXX/Experiments/SU4_twisted_Baxter_L2.wb` — Sections 6-10 hold the existing
  twist-expansion machinery (`JointSolve`, the untwisted solver `UntwSolve`/`tc0U`, the
  eps-parametrisation, the divergence bookkeeping). Sections 11-12 are the Hodge duals and the
  scalar product, with the dead ends recorded so they are not retried.

**What is already known about the eps -> 0 structure** (Sections 6-7, verified):

- `q_k = O(eps^{-(M_k - M_1)})`; finer, `[u^j] q_k = O(eps^{-(M_k - max(j,M_1))})`. The individual
  branches diverge, but in the wedge the divergences cancel, leaving one power of eps per colliding
  pair — exactly the Vandermonde `V_I`.
- `W_I = V_I x (vacuum factor) x q~_I`, `V_I = Prod_{a<a'}(z_{i_a}-z_{i_a'})`, so for the triples
  entering the duals `V_Ibar = O(eps^3)`.
- `q_1234 = A(u) A(u-h)` in closed form, state- and twist-independent.
- The regular chain is `{1} subset {1,2} subset {1,2,3} subset {1,2,3,4}` when `M` is increasing
  (state becomes highest weight on untwisting); for other states the flag reverses.

**Why this should work at all.** In the scalar product

    M^{AB}[(alpha,i,j),(r,beta)] = Int_C du E1(u) Exp[2 Pi I u (alpha-1)/h] u^beta [q^j_A q^B_i]_r(u)

the measure `E1 = mu/nu1` and its pole ladder are **completely twist-independent**, and so is the
periodic factor. ALL twist dependence sits in the two Q-functions. So the small-twist limit is
purely a statement about `q^j_A` and `q_i^B`: the `q_i` diverge as above, the duals carry
`V_Ibar = O(eps^3)` and `zeta_Ibar^{(u-h)/h} -> 1`, and the eps-powers should organise into an
overall factor fixed by the `M_k`, leaving untwisted Q-functions inside the brackets.

**TRAP that will cost real time if missed.** The scalar product requires a **unimodular** twist,
`z1 z2 z3 z4 = 1`. The existing eps-parametrisation in Sections 6-10 uses
`psi = {1, 1/3, 1/5, 1/7}` with `z_j = 1 + eps psi_j`, whose product is NOT 1. Those test states
therefore CANNOT be fed into the scalar product as they stand. Use a parametrisation that is
unimodular order by order, e.g. `z_j = Exp[eps psi_j]` with `Sum_j psi_j = 0` (which is exact, not
just perturbative), and re-derive the divergence orders if the psi's change.

**Suggested route.**

1. Pick a unimodular `psi` with `Sum psi_j = 0` and rerun the exact solver at symbolic `eps` for a
   small state (rational, so results are exact) — the `(1,1)` weights `{2,2,0,0}` and `{2,1,1,0}`
   are the cheapest, `(2,2)` `{4,4,0,0}` vs `{4,3,1,0}` the next.
2. Expand each matrix entry of `spMatFSoV` in `eps` and read off the leading power per entry;
   check whether the powers are additive in the row/column labels (they should be, if the answer
   factorises).
3. Factor out the predicted overall power and check that what remains is built from the untwisted
   Q-functions (`UntwSolve`'s `P_k`, staggered degrees `Mhat_k = M_k + (k-1)`).
4. Then do the determinant: the norm `det M^AA` is the object with a chance of a clean closed form.

**Do not** re-derive the measure, the dual normalisation, the contour, or the `T_2` split — all
four are settled and each has a recorded failure mode. If something disagrees, suspect the twist
parametrisation (unimodularity) before suspecting Section 12.

**Practical notes.** Kernel `abort` reliably fails on long `Solve`s in this project — budget for a
`restart` (which loses state; reloading the Clean notebook is ~10 s, the Experiments notebook needs
cells 3-19, 24-28, 78, 82 rerun). Zero-tests must be `Simplify[...] === 0` on exact algebraic
numbers: the determinants are identically zero, so `N[]` merely exhausts `$MaxExtraPrecision` and
warns.

## State as of 2026-09-09 (twisted gl(4) L=2 Baxter solver + Wronskian structure — new notebook)

New notebook `Paul/Mathematica/old/XXX/Experiments/SU4_twisted_Baxter_L2.wb` (**uncommitted**, 37 cells,
sections 1–7). Built step by step at Paul's direction. It replaces item 1 of the 2026-09-02 plan
below: rather than wiring `Baxter_L2_XXX_SU4.wb` into the norm machinery, this is a fresh,
**purely functional** twisted solver — no GT patterns, no monodromy matrices, no eigensystems;
everything is a polynomial identity in `u`. Paul's instruction was explicit: ignore Matteo's
notebooks entirely, and postpone the untwisting/norm thread until the twisted side is understood.

### Conventions

`z[j]` are the primary symbols with `z[j] = Exp[h phi[j]]`, so `Q_j(u) = z_j^(u/h) q_j(u) =
Exp[phi_j u] q_j(u)`, `q_j` monic of degree `M_j = S+P-mu_j` (so `Sum_k M_k = 2(S+P)`). Sites
`theta1, theta2`, shift `h`, site reps `{S,S,0,0} x {P,P,0,0}`. Everything is rational in `z`;
`phi` enters only via a substitution. `mu` labels a gl(4) **weight** (the irrep, for a HWS).

### Established (verified, not conjecture)

- **Transfer-matrix ansatz: 12 free coefficients down to 5.** Two structural facts, both Paul's,
  both then derived and checked symbolically in `S,P,theta,h,z`:
  - `tau3 = nu_1(u-2h) nu_4(u) * (quadratic, leading chi_3)`. Because `nu_1=nu_2` and `nu_3=nu_4`
    for rectangular reps, every term of `tau3 = Sum_{i<j<k} z_i z_j z_k nu_i(u-2h)nu_j(u-h)nu_k(u)`
    shares that factor. Untwisted it reproduces the `trivfactor` of `baxter_SU4_L2_XXX.wb`.
  - **Subleading coefficient, universal:**
    `[u^(2a-1)] tau_a = -a (theta1+theta2+(a-1)h) chi_a - h Xi_a`,
    `Xi_a = Sum_{|I|=a} z_I Sum_{i in I} n_i`, `n = mu` the global Cartan weight.
    Compactly `Xi_a = Sum_i n_i z_i d(chi_a)/d(z_i) =: D_n chi_a` (verified exactly, `a=0..4`), so
    the whole subleading term is a twist prefactor times `chi_a`, plus `h` times the weights acting
    on `chi_a`. Checks: `a=0` gives 0; `a=4` reproduces the `u^7` coefficient of the independently
    known quantum determinant; `a=1,2,3` reproduce the pseudovacuum eigenvalue **identically in u**.
  - Leftover unknowns: `tc[1]` (in `tau1`), `tc[2..4]` (in `tau2`), `tc[5]` (the constant in
    `tau3`'s quadratic).
  - **This formula is not written down anywhere else in the repo** — checked. `paper/xxxCG.tex` is
    gl(2)/single-site only; `Baxter_L2_XXX_SU4.wb`'s
    `alphaCoeffs` are numeric `1/u` coefficients used only to get `Mdeg`; `baxter_SU4_L2_XXX.wb`'s
    `asympSol` is untwisted-only and not in universal form; the structure-constants paper gives
    `tau±` as Q-ratios (C11)–(C12), not a large-`u` subleading. The general-`L` extension
    `[u^(aL-1)]tau_a = -a(Theta + L(a-1)h/2) chi_a - h D_n chi_a`, `Theta = Sum_alpha theta_alpha`,
    follows from the same derivation but is **untested** (only `L=2` verified).
- **The degree formula is derived, not assumed.** For an arbitrary `q` of degree `M`, the top
  coefficient of the Baxter residual is `charPoly[z_k] = 0` automatically, and the next one is
  `h z_k Prod_{j!=k}(z_k - z_j) * (mu_k - S - P + M)` — zero iff `M_k = S+P-mu_k`, given
  `Sum mu = 2(S+P)`. It degenerates exactly when twists collide, which is the untwisting limit.
- **`TwistedSolve[S,P,mu,sub]`** solves the joint bilinear system (5 `tc` + `2(S+P)` q-coefficients;
  `28+2(S+P)` equations, overdetermined by 23). One branch is *not* enough — for `mu_top` at (1,1)
  branch 1 alone leaves a 2-parameter family; all four together are unique.
- **Number of solutions = weight multiplicity**, checked against the product character
  (`wtMultProd`). At (2,2): 1, 2, 3, 8, 12, 22. Expected — the twist breaks su(4), so every state
  of that weight shares the same `M_k`.
- **Validation on the six permuted vacua** `n = pi({S+P,S+P,0,0})` (multiplicity 1, hence rational
  and exactly checkable — Paul's suggestion, and it is enough). At `(1,1),(2,1),(2,2),(3,3),(5,4)`:
  30 states, all unique, Baxter residual **exactly 0** on all four branches, and `tau` matching an
  independently computed pseudovacuum prediction. 0.03–0.2 s each. Repeated at (1,1) with **no
  substitution at all** (`theta1,theta2,h,z_j` symbolic): same 1/0/0, ~3 min.
  - **Shift ordering** for the permuted vacua is *not* the natural index order: the quantum-minor
    shifts follow the Borel adapted to `pi`, `ord = Join[A, complement A]` with
    `A = {i : n_i = S+P}`. Found empirically (searching all 24 permutations against the solver's
    own output selects exactly these). Consistent with `tau4` being ordering-independent, as it
    must be, being central.

### Wronskians (Section 6) — the normalisation `N_I`

With `Q_j = z_j^(u/h) q_j`, the `k x k` Casoratian factorises with the exponentials pulled out:

    Q_I(u) = (Prod_{i in I} z_i)^(u/h) (Prod z_i)^(-(k-1)) W_I(u),
    W_I(u) = det_{a,b} [ z_{i_a}^(k-b) q_{i_a}(u-(b-1)h) ]

`W_I` has degree `Sum_{i in I} M_i` and leading coefficient **exactly the Vandermonde**
`V_I = Prod_{a<a'}(z_{i_a}-z_{i_a'})` (immediate: the `q`'s are monic). Hence

    N_I = V_I / (Prod_{i in I} z_i)^(|I|-1),      q_I = W_I / V_I    (monic, degree Sum M_i)

and `q_I` is regular as `phi -> 0` **exactly when `I` is a lower set for the degrees `M`** (the
`|I|` branches with smallest `M_i`, ties interchangeable). For a state that becomes highest weight
on untwisting `M` is increasing, so the regular chain is `{1} ⊂ {1,2} ⊂ {1,2,3} ⊂ {1,2,3,4}` —
i.e. `q_1, q_12, q_123, q_1234`, as Paul stated. For `I={1,2}`: `N_12 = (z_1-z_2)/(z_1 z_2)`.

- **Divergence of the individual `q_k`:** `q_k = O(eps^(-(M_k - M_1)))`, finer pattern
  `[u^j] q_k = O(eps^(-(M_k - max(j,M_1))))`. Mechanism: the untwisted equation admits the staggered
  degrees `Mhat_k = M_k + (k-1)`, so a monic `q_k` of degree `M_k` must tend to an untwisted
  solution of degree `<= M_k`, and for `k>=2` the only one available is branch 1 of degree `M_1`.
  In the wedge the divergences cancel (the limit of `q_k` is *proportional to* `q_1`), leaving one
  power of `eps` per colliding pair — exactly `V_I`.
- **`q_1234` in closed form, state- and twist-independent:** from `C(u+h) = A_4(u) C(u)` plus
  `nu_2=nu_1, nu_3=nu_4`, `W(u+h)/W(u) = nu_4(u)nu_4(u-h)/(nu_1(u)nu_1(u-h))`, which telescopes:

      q_1234(u) = A(u) A(u-h),   A(u) = Prod_{j=1..S}(u-theta1-jh) Prod_{j=1..P}(u-theta2-jh)

  degree `2(S+P) = Sum M_k`. Verified exactly (residual 0) on 20 states.

### Vacuum factors and the divergence count (Section 7)

Eq. (C6) of the structure-constants paper (`References/2603.26617v1.pdf`) reads
`Q_ij = W(Q_i,Q_j)`, `Q_theta Q_ijk = W(Q_i,Q_j,Q_k)`, `Q_1234 = W(Q_1..Q_4)` with
`Q_theta = Prod_alpha (u-theta_alpha)` of degree `L`. Here `L -> S+P` and `Q_theta -> A`:

    W_12   = V_12 q~_12                deg q~_12   = M1+M2
    W_123  = V_123 * A(u) * q~_123     deg q~_123  = M1+M2+M3 - (S+P)
    W_1234 = V_1234 * A(u) A(u-h)      deg q~_1234 = 0

Verified exactly (`PolynomialRemainder` = 0, quotient monic of exactly that degree) on all six
vacuum weights at (2,1) and (3,1), the pseudovacuum at four `(S,P)`, and a genuine HWS. These
degrees are also `Sum Mhat - k(k-1)/2` with `Mhat_k = M_k+(k-1)` after removing the `A`'s — the
vacuum-factor count and the staggering are the same statement.

**Uniform count for `p_k`.** Divide the `k x k` relation by `V_{1..k-1}`; this strips every
divergence and leaves a `(k-1)`-th order difference operator `O_k` with *regular* coefficients:

    O_k[q_k] = ( Prod_{i<k}(z_i - z_k) ) * (vacuum factor) * q~_{1..k}

with (i) source of `eps`-order `k-1` and degree `Sum_{i<=k} M_i`; (ii) `O_k` at `eps=0` mapping
degree `d` to `Sum_{i<k}M_i + d - (k-1)`, its polynomial kernel being exactly the untwisted
`span{P_1..P_{k-1}}` of degrees `Mhat_1..Mhat_{k-1}`; (iii) the tower `g_0 ~ q_1` (degree `M_1`),
`deg g_j = M_1 + j`. The source enters at `j = p+k-1`, so

    Sum_{i<k}M_i + (M_1+p+k-1) - (k-1) = Sum_{i<=k}M_i   ==>   p_k = M_k - M_1

uniformly in `k`. For `k=2` this is a fully independent first-order/`eps^1` derivation (two ways:
monicity, and the inhomogeneity). Input (ii) is **verified explicitly for `k=3` only** — on state A
the kernel sits at exactly `d = Mhat_1, Mhat_2 = 0, 2` and `deg O_3[f] = M1+M2+d-2` elsewhere; on
the rational vacua `O_3` degenerates to the plain second difference `f(u)-2f(u-h)+f(u-2h)`. For
`k=4` the pattern is asserted, not checked.

### Section 11 (2026-09-09): Hodge-dual `Q^j` — the prefactor `F(u)` is fitted

Convention (asymmetric, Section 6's `Wronsk`; the paper's symmetric (C5) is not used):
`F(u) Q^4(u-h) = det[Q_1(u), Q_2(u-h), Q_3(u-2h)]`, `Q^j = (-1)^j Q_{comp j}`, `F` fixed by
state-independence with the dual polynomial monic. With `I = complement of j`,
`zeta_I = prod z_i`, `V_I` the Vandermonde on `I`, and `A` the Section 7 vacuum factor:

    Q^j(v) = (-1)^j zeta_I^(v/h) q^j(v),  q^j(v) = (W_I/V_I/A)(v+h)  monic, deg q^j = mu_j
    F_j(u) = (-1)^j (V_I/zeta_I) A(u),    A(u) = Prod_{k=1..S}(u-theta1-kh) Prod_{k=1..P}(u-theta2-kh)

So `F` is `A(u)` up to a constant: degree `S+P`, `S` factors with `theta1` and `P` with `theta2`;
the paper's `(u-theta1)(u-theta2)` is the `S=P=1` case. Verified 16/16 on `(2,1) {3,3,0,0}`,
`(1,1) {2,1,1,0}`, `(2,2) {4,3,1,0}` at `h=1`, plus `(2,1) {3,3,0,0}` at `h=3` with different
`theta`. Cells 77-80; symbols `QQfun`, `WasymQ`, `qDual`, `QDual`, `Ffit`, `Gfit`, `numH`.
The pairwise `det[Q^i(u),Q^j(u-h)]` relation holds up to `Ztot^(u/h)` times a pair-dependent,
state-independent constant (normalisation bookkeeping).

Still open: dual Baxter equation `B† Q^j = 0` unchecked; `Q^j` not yet expanded in the twist.
NEW TRAP: a definition with its RHS split over several lines is silently truncated at the first
line in this notebook — keep every RHS on one line.

### Section 12 (2026-09-09): O, O-dagger and the orthogonality identity — VERIFIED

Conventions: `f^[±2n] = f(u ± n h)`; rotated Gamma with `Gam[F(u+h)]/Gam[F(u)] = F(u)`.

- **`O` (12.38), 4th order** = this notebook's Baxter equation divided by `Prod_{m=0..3} nu1(u-mh)`
  with `G(u)=Q(u+h)` and `T_{a,1}=tau_a`. Same operator for all four branches — the `z_k` in
  `BaxRes` are only the stripped exponential. `O q_i^[2] = 0` verified for all `i`.
- **`O†` (12.39), 4th order, `nu3 -> nu4`**: `O† q^j = 0` verified for all `j`, but only for the
  UNCENTRED dual `q^j(u) = QDual[j](u-h)`. Section 11's `QDual` has a `+h` centring baked in.
- **Measure**: `mu(u) = 1/Prod_alpha Prod_{n=0}^{S_alpha-1} (u - theta_alpha - n h)` — `S_alpha`
  one-sided poles per site at `theta_alpha + n h`. `[2,2,0,0]^2`:
  `1/((u-th1)(u-th1-h)(u-th2)(u-th2-h))`.
- **Identity (12.42)** `0 = Contour_C mu q^j_A (O^A - O^B) q_i^{B[2]}` holds exactly — the contour
  must encircle ALL finite poles of the integrand, not just `mu`'s (the `O` coefficients add
  `nu1(u-mh)` denominators: 6 poles at (1,1), 8 at (2,2)). Verified at (1,1) and (2,2), exact
  non-perturbative twisted solutions, unimodular twist `numU`. Both terms vanish separately;
  plain brackets are nonzero.

Cells 81-84; symbols `OpO`, `OpOd`, `qhat`, `qdualhat`, `mMeas`, `contourSum`, `relCheck`, `idTest`.

DEAD ENDS (do not retry): the `theta_alpha +- (2l-1)h/2` ladder with the paper's 8x8 determinant
(zero at (1,1) by luck, fails at (2,1),(2,2)); per-pole functionals with `beta=1..S+P`; fitting
`(a,b)` in `theta-a h, theta, theta+b h` by determinant vanishing — the root-find runs away to
`a~10^3` where the determinant underflows trivially.

NEXT: the determinant/orthogonality packaging on top of the verified identity.

### Section 12.3-12.5 (2026-09-10): all four integrand terms on ONE measure

Trivial factors: `T_1 = t_1`, `T_3 = nu4 nu1^[-4] t_3`, `T_4 = qdet` (cancels). `T_2` has none but
is degree 4; split per state `T_2 = nu1(u+p1) t21 + nu4(u+p2) t22` (deg 2 each; one-parameter
ambiguity `t21 -> t21 + c nu4`, fixed by `t21(0)=0`, cancels in the sum).

Effective measures are contiguous ladders `1/(u-theta_alpha-nh)`:
`E1 = mu/nu1` (n=0..S_alpha), `E2 = mu/(nu1 nu1^[-2])` (n=0..S_alpha+1),
`E3 = mu nu4/(nu1 nu1^[-2])` (n=1..S_alpha+1). Key identity **`E3(u+h) = E1(u)`**.

`p1=-h`, `p2=0` are forced (generically) by demanding the leftover ladder stay contiguous on both
sites. Then `u -> u+h` in `I_22` and `I_3` puts everything on `E1`:

    0 = Int du E1(u) [ -Dt1(u) qA^j(u) qB_i(u) + Dt21(u) qA^j(u) qB_i(u-h)
                       + Dt22^[2](u) qA^j(u+h) qB_i(u) - Dt3^[2](u) qA^j(u+h) qB_i(u-h) ]

The four shift combinations `(0,0),(0,-2),(+2,0),(+2,-2)` are exactly the `a,b=+-1` pattern of the
FSoV determinant, with degree-<=2 coefficient functions. Verified EXACTLY (Simplify) at (1,1) and
(2,2), all four (i,j). Cells 84-86; symbols `splitT2`, `E1meas`, `alignedTerms`.

NEXT: this is the natural starting point for the determinant packaging.

### PROMOTED TO CLEAN (2026-09-10)

`Paul/Mathematica/old/XXX/Clean/SU4_SoV_ScalarProduct_L2.{wb,md,tex,pdf}` — the final scalar product
only, 9 cells, self-contained (verified by running it top to bottom in a FRESH kernel; the
orthogonality table takes ~48 s). `.md` is the symbol reference, `.tex`/`.pdf` (6 pp) the
derivation notes. The Experiments notebook keeps the full history and the dead ends.

### Section 12.6-12.7 (2026-09-10): SoV SCALAR PRODUCT WORKS — orthogonality at (1,1) and (2,2)

Gauge for the T2 split (Paul's): `[u^2]t21 = [u^2]t22 = (1/2)[u^4]T_2 = chi_2/2`. Unique, and since
`chi_2` is state-independent both `u^2` terms cancel in the difference => every `Delta t` is linear,
`Delta t = DeltaI_2 u + DeltaI_1`. (The earlier `t21(0)=0` gauge is WRONG - do not use.)

So the aligned identity is a linear relation among 8 IoM differences, giving an 8x8:

    M^{AB}[(alpha,i,j),(r,beta)] = Int_C du E1(u) Exp[2 Pi I u (alpha-1)/h] u^beta [q^j_A q^B_i]_r(u)

rows `(alpha,i,j)`: alpha=1,2 (alpha=2 inserts the h-periodic `Exp[2 Pi I u/h]`, inert under the
alignment shifts), `i in {1,2}`, `j in {3,4}`; columns `(r,beta)`: the four shift patterns
`(0,0),(0,-h),(+h,0),(+h,-h)` against `beta = 0,1`; poles = the E1 ladder `theta_alpha + n h`,
`n=0..S_alpha`.

RESULT: `det M^{AB} = 0` for A != B, `det M^{AA} != 0`. Verified EXACTLY (Simplify) at (1,1) and
(2,2), both orderings, INCLUDING same-weight pairs (`{2,1,1,0}` #1|#2; `{4,3,1,0}` #1|#2) - the
sharp test. 8/8 cross determinants zero, 5/5 norms nonzero.

Cells 87-88; symbols `ladE1`, `perFac`, `brkE1`, `spMatFSoV`, `orthTest`.

### Test states (for reuse)

`z_j = 1 + eps psi_j`, `psi = {1,1/3,1/5,1/7}`, `theta1=1/3, theta2=1/7, h=1` — rational in `eps`,
so `TwistedSolve` runs with `eps` fully **symbolic** and results expand exactly; no numerics, no
state tracking. States used: **A** = (1,1) `mu={2,1,1,0}`, `M={0,1,1,2}`; **B** = (1,1) lowest
weight `{0,0,2,2}`, `M={2,2,0,0}` (flag reversed — regular chain is `{4}⊂{3,4}⊂...`); **C** = (1,1)
singlet `{1,1,1,1}`, `M={1,1,1,1}` (control: all `M_k` equal, so nothing diverges); **D** = (3,1)
`mu={3,3,1,1}`, `M={1,1,3,3}` (the only one with `M_1>0` *and* real divergence).

Identifying which solution of a degenerate weight sector becomes highest weight: solve the
**untwisted** problem first (`UntwSolve`, staggered degrees `Mhat`, standard gauge fixing) — it is
unique and rational — then pick the twisted branch converging to that `tc0`, validated by
`dist(tc,tc0) ~ eps`. Cheaper and far more robust than degree-spectrum matching on 7 algebraic
branches.

### Traps hit this session (all cost real time)

- `Limit[]` applied to a whole rational expression is **unreliable** on the algebraic
  (multiplicity>1) states — it returned *different branches on repeated runs* of identical code.
  Take the limit coefficient-by-coefficient instead.
- A helper that `Module`-localises the spectral variable returns expressions in a scoped `x$nnn`;
  the caller's `CoefficientList` then silently reads **every coefficient as degree 0**. Pass the
  variable in.
- Detecting admissible untwisted degrees by nullspace **dimension** silently misses a jump when two
  admissible degrees are adjacent (gave `{0,2,4,5}` instead of `{0,1,4,5}` for irrep `{2,2,0,0}`).
  Test **monic solvability** instead; that agrees with `degQ` of `baxter_SU4_L2_XXX.wb`.
- `RootReduce` is exact but far too slow to sweep with; bare `PossibleZeroQ` falls back to an
  undecided numeric test and warns. Hence the rigorous checks are done on the **rational** (vacuum)
  states, where `Expand[Together[...]]` collapses to a literal 0.
- Dividing by `V ~ eps^6` is catastrophic numerically: the low coefficients of `q_1234` read as `0`
  at small `eps` (visible in the accuracy annotations). Confirm such things exactly, not by `N[]`.
- Kernel `abort` again did not take on a long `Solve`; a `restart` was needed (loses state) —
  same as recorded on 2026-09-02.

### Not done

- **The actual perturbative expansion of the polynomial coefficients** — Paul's original step 2. He
  interrupted it to understand the Wronskian first, and that is now understood, so this is the
  immediate next step. The `q_j` stay `Exp[phi_j u] x` polynomial of fixed degree `M_j`
  non-perturbatively; only the coefficients are expanded, and they need Laurent series (orders
  `eps^(-(M_k-M_1))`, now known in advance).
- Nothing is wired to `SU4_norm_untwisting_L2_general.wb` / the FSoV norm yet.
- No scan over reps (item 3 of the 2026-09-02 plan).
- Write-up: the subleading-`tau` formula and the `N_I`/vacuum-factor ledger are worth recording, but
  `paper/xxxCG.tex` is currently gl(2)/single-site only, so adding a gl(4) L=2 section is a scope
  decision Paul has not made. A `Clean/` promotion + `wolfbook-docs-agent` reference doc is the
  other option.

## State as of 2026-09-02 (untwisting the FSoV norm — gl(4) L=2; two notebooks, results below)

Thread: **how to untwist the norm analytically**, started from Appendix N of the new reference
`References/2603.26617v1.pdf` (Bargheer–Bercini–Lefundes–Ryan, "Structure Constants from Q-Systems
and SoV"; Paul is an author). Also new: `References/SU4_twisted_any_L_Qs.nb` — Matteo/Paul's twisted
Q-solver for any L, 6-dim `[1,1,0,0]` site rep.

Two notebooks, both under `Paul/Mathematica/old/XXX/Experiments/`:

- `SU4_norm_untwisting_L2.wb` — LR-symmetric branch (`z1 z4 = z2 z3 = 1`), staged untwisting
  `phi1 -> 0` then `phi2 -> 0`. Ported solver validated against the paper's (P1).
- `SU4_norm_untwisting_L2_general.wb` — general unimodular twists, independent `Q^3, Q^4`,
  high-precision NSolve route. Validated twice: reproduces (P1) exactly, and reproduces the
  first notebook's `B = -0.0274316528219221` to all digits.

### Established (verified, not conjecture)

- `B` (eq. (12)) is **exactly invariant** under `GL(2) x GL(2)` on `(Q1,Q2)` / `(Q^3,Q^4)`,
  rescalings included: the row mixing gives `(det H)^{2L}`, the normalisation (15) gives
  `(det H)^{-2L}`. So the rotation `Q2 -> Q2 + r Q1` is **pure gauge** — it cannot change `B`,
  only whether the limit is computable entrywise. Checked numerically to 20+ digits.
- **Untwisted degrees**: `lambdahat_j = lambda_j + (j-1)`. Verified on two states with different
  predictions: `{0,1,0}` (`lambda={0,1,1,2}`) gives `lambdahat={0,2}` with
  `Q1hat=1, Q2hat=u^2/2+1/3`; `{1,2,1}` (`lambda={1,1,1,1}`) gives `lambdahat={1,2}` with
  `Q1hat=u, Q2hat=u^2+1/3`. Mechanism: subtracting the collision leaves the next order of the
  `Exp[phi u]` expansion, raising the degree by one.
- **The rotation**: fixed by a quadratic with a *double* root. For `{0,1,0}`, coefficients
  `2I, -4I, 2I` i.e. `2I(c-1)^2`, giving `r = 1/phi2`, i.e. subtract the divergent multiple of
  `Q1` (there `p2 = u - Cot[phi2]`, so `Q2 -> -Q1/phi2`).
- **Why the rotation is mandatory**: for `{0,1,0}` the naive entrywise min-cost predicts
  `phi2^(-2)` while the truth is `phi2^12` — 14 orders of cancellation.
- **Everything is a contraction with a universal moment table.**
  `M[(alpha,j,k),(beta,a,b)] = Sum_n c_n mu[n,alpha]`, `mu[n,alpha] = <u^n>_{L,alpha}`.
  For `L=2`: `mu[n,1]` is nonzero **iff n is odd and n >= 2L-1** (parity selection rule, proved by
  closing at infinity: `mu[n,1] = -I Binomial[-L,k] 4^(-k)`, `k=(n-2L+1)/2`); `mu[n,2]` has the
  single zero `mu[1,2]=0`. This explains both previously-unexplained patterns: the `alpha=2`
  vanishing (row `(2,1,4)` of `{0,1,0}` *is* `mu[beta-1,2]`) and the `a=-b` zeros in block `(2,3)`
  (the surviving `u^3` coefficient is `I(a+b)/4`).
- **Overall exponent**: along a ray, `det` (with (12)'s normalisation) vanishes as `eps^nu` with
  `nu = -Sum_{j<k}(lambda_j - lambda_k) = -Sum_j (n+1-2j) lambda_j`. Verified: `{0,1,0}` gives
  `nu=6` (measured ratios 5.963, 5.991, 5.998, 5.999) against `Sum lambda_jk = -6`; `{1,2,1}` gives
  `0` against `0`. This is exactly Appendix L's (L1) prefactor `Prod z_jk^(lambda_jk)`, so
  **(L1)'s `B`, not (12)'s, is the object with a uniform limit** — one determinant, not one per state.
- **Closed results.** `{0,1,0}`: `B -> -32 phi1^2 phi2^4`, where `-32 = -4 * 8` and `8` is the
  reduced minor built from untwisted Qs (`Q12hat = -I u`). `{1,2,1}`: corank 0, so the plain
  untwisted 8x8 determinant already works and gives `B = 7776 = 6^5` exactly, matching the twisted
  limit (7477.7, 7700.4, 7757.0, 7771.3, 7774.8 -> 7776).

### Caveat and what's missing

(L1) makes `B` finite, but when the corank is positive you have `0 * Infinity`: the untwisted
determinant is literally 0 and the prefactor literally infinite. To *evaluate* the answer from
untwisted Qs you still need the reduced minor, whose size is the corank.

**Only two highest-weight states so far** (`{0,1,0}` = `[1,0,1]`, `{1,2,1}` = zero-weight
subspace, 6 solutions = 2+3+1). At `L=2` every state with non-negative Dynkin labels is
automatically LR-symmetric (`6 (x) 6 = 20' + 15 + 1`, all self-conjugate ⇒ `M = Nn`), so the
6-dim site rep cannot supply more data. Hence the plan below.

### Plan agreed for the next session (Paul, 2026-09-02)

> **Status as of 2026-09-09:** item 1 is done, but *differently* — instead of wiring
> `Baxter_L2_XXX_SU4.wb` into the norm machinery, a fresh functional solver was built
> (`SU4_twisted_Baxter_L2.wb`, see the 2026-09-09 section above). Item 2 is the immediate
> next step and is now much better set up (the `eps`-parametrisation makes exact symbolic
> solves possible, and the divergence orders are known in advance). Item 3 is untouched.

1. **Generalise the determinant to `[S,S,0,0] x [P,P,0,0]` site reps.** *Mostly already done*:
   `Experiments/Baxter_L2_XXX_SU4.wb` (2026-07-19 entry below) is already a **twisted** gl(4) L=2
   Baxter solver for two independent **generic Young-diagram** reps per site, with four distinct
   diagonal twists and the rep-dependent vacuum factors
   `nu[lam1,lam2][i][u] = (u - th1 - h lam1[[i]])(u - th2 - h lam2[[i]])`,
   `qdetT = chi4 Prod_i nu[i][u-(4-i)h]`. Full-sweep verified: 152 states, 5 rep pairs,
   worstTQ = worstCasoratian = 0 exactly. Same Q convention as here (`Q_k = z[k]^(u/h) q_k`).
   So `[S,S,0,0] x [P,P,0,0]` is a special case and arbitrary Young diagrams come free — the task
   is wiring it to the norm/bracket machinery in `SU4_norm_untwisting_L2_general.wb`, not deriving
   new taus. (The other notebook, `baxter_SU4_L2_XXX.wb`, has Casimir-form taus for `[S,S,0,0]` but
   **untwisted** — not the one to use here.)
   **Precision trap** recorded in that entry: hardcoded `N[..., 30]` literals inside
   `TauEigensystem4`'s `rSets`/`Kmix` silently capped the whole downstream chain regardless of
   `z[k]`'s precision. Grep for bare precision literals before trusting a scan.
2. **Build a perturbative-in-twist solver.** Every failure this session was NSolve
   ill-conditioning plus state-tracking jumps as `eps -> 0`, because the solutions genuinely
   coalesce there. A series solver in `eps` avoids it, and is structurally right: the leading
   order *is* the degenerate untwisted problem, so this is degenerate perturbation theory and
   `r` should drop out of the order-`eps^0` solvability condition instead of being fitted.
3. **Scan over reps** to deduce the untwisting pattern. Per state record: `lambda_j`;
   `lambdahat_j` (test `lambda_j + (j-1)`); corank and which entries vanish; `nu` against
   `-Sum_{j<k} lambda_jk`. Compute `mu[n,alpha]` once per `L` and check whether the odd-`n`
   rule and `mu[1,2]=0` persist at other `L`.

### Performance notes (matter for the scan)

- Do **not** use `Residue` for the FSoV bracket. Precompute `mu[n,alpha]` symbolically once and
  contract with coefficient lists: ~0.1 s per 8x8 norm instead of minutes.
- Symbolic `Solve` on the general (non-LR) system returns branches with `Sqrt` of perfect squares
  that resist `PowerExpand`; that is why the reference notebook uses `NSolve`. The perturbative
  solver should sidestep this entirely.
- Kernel `abort` did not take on a long `NSolve` sweep; a `restart` was needed (and loses state).

## State as of 2026-08-27 (su(2) SoV localization — teaching notebook, complete; FSoV link pending)

`Mathematica/su2/su2-localization.wb` (uncommitted, along with `Mathematica/su2/su2-localization.md`).
Paul wrote this notebook specifically to teach the localization structure of the SoV basis; Claude
extended it this session with the mirror construction. **Read `Mathematica/su2/su2-localization.md`
first** — it has the conventions, symbol tables, the verified-check table, and the traps.

Established and symbolically verified (L=2 gl(2), companion twist, `Lam[1]=Lam[2]={3,0}`): the SoV
basis built up from the minimal vacuum factorizes on the `mu2=0` edge, `x[{mu,0}] = x1[{mu}] (x)
vac2`, exactly; the mirror basis built down from the maximal `B`-eigenvector factorizes on the other
edge, `xb[{0,nu}] = vacb1 (x) x2b[{nu}]`, also exactly. Both towers build the same basis (`B` has 16
distinct eigenvalues). 7 of the 16 states factorize, not all. Localization is a property of the
vacuum, not of the operator.

Localization is stronger than a vacuum statement: it holds on the whole orbit of the vacuum under
the **global** (diagonal) su(2) rotation, `<0| Exp[r EE[{2,1}]] . t[{1,mu}][w1] == <0| Exp[r
EE[{2,1}]] . (t1[{1,mu}][w1] (x) 1)`. Globality is essential — independent per-site parameters break
it, restored exactly when they agree.

**Next step, from Paul — to be set up explicitly in a later session.** The strategy (sketched, not
yet implemented or verified; full write-up in `su2-localization.md` §7): globally rotated states
`<Omega(r)|` have `<Omega(r)|Psi>` known exactly in Q-functions; localization means a transfer
matrix acting on such a state hits **site 1 only**; certain **combinations of `t`'s sum up to global
rotations**, so such a combination rotates site 1 alone and produces `<Omega(r'), Omega(r)|` with
independent parameters per site; and since `t|Psi> = tau|Psi>` the same overlap is `tau` times the
known one. Net: `<Omega(r'), Omega(r)|Psi>` in Q-functions — the two-parameter family CG
coefficients are read off from. The crux is that the input orbit is one-parameter and the target is
two-parameter, and localization is what breaks the diagonal constraint. Don't start building this
ahead of that session.

Standing constraint restated: **this project only ever deals with two sites.** Unequal reps behave
identically, so the equal-rep choice is convenience, not a special case.

## State as of 2026-08-19 (Matteo's gl(4) Q-space/normalisation notebook imported — read/analyzed only, NOT yet computed on)

New material from Matteo, committed this session (`de97e03`), not yet acted on beyond a structural
read-through. Four files, three of which are the *same* underlying notebook in different formats:

- `Paul/Python/SU(4)_Paul_Normalisation_EXPLICIT.nb` — Matteo's original Mathematica notebook.
- `Paul/Python/SU(4)_Paul_Normalisation_EXPLICIT.wb` — wolfbook copy of the same (open this with
  `wolfbook-builder`/the wolfbook MCP tools to compute live, cell-by-cell, on either machine).
- `Paul/Python/SU4_Paul_Normalisation_EXPLICIT.wl` — a flattened, **directly `Get`-able** script
  version (presentation/output cells stripped; first line says so). Good for `grep`/quick reads or
  a plain `wolframscript`/kernel run; NOT meant to be run top-to-bottom as one script — it's full of
  bare unassigned "example of use" expressions (`Sol=SolveBaxAll[...]` etc.) meant for selective
  cell evaluation, and later sections redefine earlier-named symbols in the 1-site block.
- `Paul/Matteo_Files/SU(4)_Paul_Normalisation_SEND.wl` — **not directly runnable**: this is a raw
  Mathematica "notebook-as-text" export where every cell, *including Input cells*, is wrapped in
  `(* ::Input:: *) (*code*)` comment markers (presumably literally what Matteo emailed/sent). Do not
  `Get[]` it expecting code to run — strip the cell-marker comments first, or just use the `.wb`/`.nb`
  above, which already have live code.

**What the notebook computes** (from reading `SU4_Paul_Normalisation_EXPLICIT.wl` top to bottom, not
yet cross-checked with Paul or numerically re-verified beyond what the notebook's own test cells
already assert): a Gaudin-side two-site (`w1=0,w2=1`) gl(4) FSoV construction —
1. **Ψ-space**: GT-basis generators (`ps/dim/J0/Jp/Jm/JJ`, same machinery as
   `Gaudin/Clean/Hamiltonians_SU4.wb`), Lax connection `LIndex`, explicit `H1..H4` and dual
   `H1†..H3†` Hamiltonians (regulated traces of powers of `L`), and `EigenR/EigenL` — right/left
   joint eigenvectors obtained by matching against the **Q-side** spectrum rather than a bare
   `Eigensystem` call.
2. **Q-space**: `h[rep1,rep2][k,1/2]` (residues at `w1,w2`, closed gl(4)-weight polynomials) plus
   twist elementary-symmetric pieces build the 4th-order Baxter operator `BaxOp`; `SolveBaxAll`
   solves all four Baxter equations for the monic `Q[i]` coefficients **and** the leftover
   Hamiltonian unknowns simultaneously (spectrum + Q-functions in one solve). A dual sector `Q†`
   is built two ways — directly via a dual Baxter equation, and via `Hodge` (Wronskian of the other
   three `Q`'s) — with a "Check Q Hodge" block cross-checking the two agree.
3. **Renormalised scalar product**: `OrthogRen` — a 12×12 residue-moment determinant (poles at
   `z=0,1`) divided by a `Normalisation` built from boundary limits of `Q†_1`, `W(Q†_1,Q†_2)`,
   `Q_4`; the analogous 1-site pieces are `OrthogRen1/OrthogRen2`.
4. **Ψ-vector normalisation (the "EXPLICIT" part)**: `VacR/VacL` build the vacuum state as an
   ordered `MatrixExp` product of the six raising/lowering `EMAll[a,b]` generators with parameters
   `α_ab,β_ab`; `αr[12,24,34]`/`βr[...]` are then **closed-form rational functions of the twists
   `φ`** that fix those parameters (leaving 3 free: `α13,α14,α23`), giving renormalised eigenvectors
   `ReEigenR/ReEigenL` in the same normalisation as the Q-side `OrthogRen`.
5. **`K` vectors** (final ~230 lines): `K22B`/`K22C` — explicit overlap vectors between a 1-site
   and a 2-site factor (SU(2)⊗SU(3)-split), built the same way as `VacR/VacL` but with extra
   `k3,k4`-dependent rational dressing. `TESTK22B/C` + the "Test reconstruction" cells verify, for
   `{1,1,0,0}×{1,1,0,0}` only, that `K·Ψ` (solved for by expanding in `k2,k3,k4`) reproduces the
   `ReEigenR` vector obtained independently via step 4 — i.e. this is a from-scratch construction and
   check of exactly the kind of normalised-eigenvector ↔ Q-function overlap identity that the
   `Gaudin/Clean/Hamiltonians_SU4.wb`/CG-overlap program on the XXX side has been chasing (see
   memories `project_cg_overlap_dictionary`, `reference_fsov_determinant_normalization`) — but done
   independently by Matteo on the Gaudin/FSoV side. **Not yet compared against that program in any
   way; that comparison is the obvious first thing to do before extending this notebook.**

**Two likely bugs, found by inspection, NOT yet fixed (present identically in `.wb`/`.nb`/`.wl`,
so not a conversion artifact) — check before relying on either symbol:**
- `VacR1[rep1_,β12_,...,β34_] := VacR1[rep1,β23,β12,β13,β14,β23,β24,β34] = ...` (1-site vacuum,
  `.wb` line 2235/3829) — the memo-pattern LHS takes 8 arguments (extra leading `β23`, and `β23`
  appears twice) against the 7-arg definition head, so it never matches an actual 7-arg call.
  Effect: memoization is simply dead (recomputes + leaves junk defs behind each call); results are
  still numerically correct, just uncached.
- `K22C` (`.wb` line 4088) — in the first `KroneckerProduct` factor (built on `UPL1[rep1]`), only
  the `E13` exponential uses `EM[...,rep1]`; the `E14`, `E23`, `E24` ones use `EM[...,rep2]` instead.
  Only evaluates without error when `dim[rep1]==dim[rep2]` (true for the notebook's only test case,
  `{1,1,0,0}×{1,1,0,0}`) — will break or silently misbehave for `rep1≠rep2` (e.g. `{p,p,0,0}×
  {q,q,0,0}`, `p≠q`). Compare against the parallel, apparently-correct `K22B` just above it in the
  notebook when fixing.

**To pick this up on a different machine:** `git pull` gets all four files (pushed in `de97e03`).
Needs a Wolfram kernel; use `wolfbook-builder`/wolfbook MCP tools against the `.wb` for live,
checkpointed cell-by-cell computation (matches how `Gaudin/Clean/` work is normally done in this
repo) rather than trying to batch-run the flattened `.wl`. No non-core-WL package dependencies were
seen anywhere in the four files.

**Next steps (none of this has been discussed with Paul yet — treat as an inbox item, not a plan):**
- Fix the two bugs above if/when extending past `{1,1,0,0}×{1,1,0,0}`.
- Work out how this notebook's `OrthogRen`/`ReEigenR`/`K22B`/`K22C` machinery relates to the
  existing XXX-side CG-overlap dictionary and FSoV-determinant-normalisation work — likely the
  actual point of importing it, but unconfirmed.
- This currently sits under `Paul/Python/`, which per `CLAUDE.md` is meant for the sympy Python
  cross-check port, not Mathematica notebooks — a Mathematica `.wb`/`.nb` pair there is a layout
  mismatch. Don't relocate without asking Paul first (same rule already applied to
  `Hamiltonians_SU4.wb`'s Gaudin/XXX placement).

## State as of 2026-07-22 (C++ CGC reference benchmark set up — `Paul/Cpp/AlexCGC/`, UNCOMMITTED)

Separate thread from the Baxter/Gaudin work below: set up the **state-of-the-art SU(N) CGC code**
(Alex, Kalus, Huckleberry, von Delft, arXiv:1009.0437 — the Gelfand–Tsetlin pattern-calculus C++
program from that paper's appendix) as a **reference benchmark** against Paul's Q-function CGC
generating function (`Paul/Mathematica/old/Gaudin/Clean/CGGenFun.wb`). New directory `Paul/Cpp/AlexCGC/`
— a C++ cross-check sibling to `Paul/Python/`. **Everything uncommitted** (Paul's choice: place,
don't commit; a compiled binary is machine-specific — rebuild from source elsewhere).

**Toolchain:** needs LAPACK/BLAS. `sudo apt-get install liblapack-dev libblas-dev` was run this
session (installed, 3.10.0). Build all binaries with `bash Paul/Cpp/AlexCGC/build.sh` (plain g++).

**What's there:**
- `ClebschGordan.cpp` → `clebsch` — the authors' original source (downloaded, unmodified); an
  interactive menu calculator (menu 5/6 = CGCs for S⊗S′). `hwcgc.sh` filters its output to just the
  highest-weight state (always the largest index `Q(M″)=dim(S″)`).
- `cgc_hw_bench.cpp` → `cgc_hw_bench` — modified copy that times ONLY the highest-weight-state solve
  (`compute_highest_weight_coeffs`, the raising-op nullspace/SVD) and SKIPS building descendants
  (3 clearly-commented `// added` edits: a `highest_weight_only` ctor flag, a scope timer, a new
  benchmark `main`). Verified `max|hw_only − full| = 0` vs the full run in every case. Modes:
  `cgc_hw_bench N S… S′… S″… [reps] [verify]` (per-irrep timing); `… decompose N S… S′…` (fast μ
  list); `machine` flag = one parseable RESULT line.
- `test.sh S [reps]` — for SU(4), times the highest-weight CGCs of EVERY μ in `[S,S,0,0]⊗[S,S,0,0]`
  (summary table only, no coefficient dumps). `reps` default 51 (report min-over-reps, for the fast
  small cases); use `reps=1` for large S (each solve already runs for seconds).
- `results/` — `hw_timings.txt` (full per-μ tables S=1..8), `summary.csv`, `gen_results.sh`
  (regenerates both). Machine: i5-1345U, 16 GB.

**The benchmark result (the point of all this):** C++ cost is dominated by the SINGLET (0,0,0,0)
— its highest-weight linear system is the largest — and total-over-all-μ hw-solve grows ~×5 per
step in S: S4=0.057, S5=0.355, S6=2.49, S7=15.3, S8=80.2 s. It also has a **memory wall**: the
singlet's dense SVD matrix is ≈ O(dim(ν)³) ≈ tens of GB by S=10 (dim[10,10,0,0]=1716), so the C++
**cannot even run S=10** on a 16 GB laptop. Paul's GF (build-for-all-irreps time) is ~**flat**:
S4=0.063, S6=0.139, S8=0.171, **S10=0.35 s**. ⇒ crossover ~S=4; GF ~470× faster at S=8, and at
S=10 it's "0.35 s vs won't-run" — a scaling win in **both time and memory**, not a constant factor.

**Open / next stage (Paul):** the GF timing is currently **build-only** — the explicit CGCs still
need to be EXTRACTED from the generating function (Paul's stated next implementation stage). That is
the one piece needed to make the head-to-head deliverable identical to the C++ (which outputs
explicit coefficients); then drop GF per-μ times into a second CSV beside the C++ `results/` for a
direct side-by-side. See memory `project_gf_cgc_vs_cpp_benchmark`.

## State as of 2026-07-22 (gl(4) L=2 untwisted τ1,τ2,τ3 ALL DONE exactly via Baxter asymptotics; Gaudin h→0 limit + standalone `GaudinSolver.wb` — DONE)

New notebook `Paul/Mathematica/old/XXX/Experiments/baxter_SU4_L2_XXX.wb` (Paul's own fresh build,
independent of `Baxter_L2_XXX_SU4.wb`) closes out the τ2 gap left open in the 2026-07-20 entry
below, via a genuinely different method — large-`u` asymptotic expansion of the untwisted Baxter
(TQ) equation, rather than operator-level Casimir fits. **This supersedes that entry's "τ2 —
PARTIAL" status: all three transfer matrices are now closed-form in Cas2(μ), Cas4(μ).**

### The mechanism (Paul's idea, worked through live)

- `n_k = μ_k` on the HWS (mode number = weight component) — verified directly against the main
  notebook's `nvec` inversion.
- Untwisted (`z_k→1`) Q-degree formula: `deg q_k = (S+P+3) − μ̂_k = (S+P)−μ_k+(k−1)`, `μ̂_k=μ_k+4−k`
  the *same* shifted weights that appear in the paper's Casimir generating function (eq 2.13–2.15,
  `G4[u]=Prod_i(u−h(μ_i+4−i))`). The `+(k−1)` staggering is forced, not just plausible: the
  twisted-case `Mdeg` formula's denominator is (up to `z³`) the derivative of the twist
  characteristic polynomial, which has a 4-fold zero at `z=1` — L'Hôpital-resolving that
  degeneracy is exactly what produces the staggered degrees.
- Key structural fact: `Σ_j(−1)^j C(4,j) j^m = 0` for `m=0,1,2,3` (a 4th finite difference) ⇒
  large-`u` orders `1/u,1/u²,1/u³` of the Baxter equation are **state-blind** (independent of the
  Q-data), and since the 4 Q-branches share one transfer matrix but have distinct leading powers
  `p_k`, every power of `p` must vanish separately. Solving that system pins 5 of τ1/τ2/τ3's 7
  unknown `u`-coefficients (`a[1],a[2],b[1],c[2],b[3]`) in terms of `S,P,θ` and one remaining
  scalar `b[2]` — leaving only `b[2]` (τ2's `u²` coeff, carries Cas2) and `d[2]` (τ2's constant,
  carries Cas4) as genuine Casimir-dependent data.
- Pseudovacuum anchor (`μ_top=[S+P,S+P,0,0]`, all τ's = elementary-symmetric-in-ν vacuum
  eigenvalues) gives `b2vac,c2vac,d2vac` independently — **cross-check passed exactly** (vacuum
  `c[2]` matches the asymptotically-predicted `c[2]`, residual 0).

### Result (τ1, τ2, τ3 all closed, all anchored at ΔCas = Cas(μ)−Cas(μ_top))

```
τ1(u) = τ1vac(u) + (h²/2)ΔCas2
τ2(u) = τ2vac(u) + (h²/8)ΔCas2·[τ1vac(u)+τ1vac(u−h)] + h⁴(⅛ΔCas2²−¼ΔCas4) + h⁴(1+2S+2P)/4·ΔCas2
τ3(u) = trivfactor(u)·[ 4u²−4u(θ1+θ2)+4θ1θ2−2hu(S+P+4) + b3vac + (h²/2)ΔCas2 ]
```

**The bug and how it was actually caught:** the first τ2 guess (fusion form
`(h²/8)ΔCas2·[τ1vac(u)+τ1vac(u−h)]` matching u²,u¹ exactly, plus the "obviously right" quartic
`h⁴(⅛ΔCas2²−¼ΔCas4)`) was missing the last `h⁴(1+2S+2P)/4·ΔCas2` term. Found by building the actual
polynomial Baxter solver and testing it: for `μ={2,1,1,0}@(1,1)` (`deg q1=0`), a constant `q1`
solved the equation at `μ_top` but **not** off-top — residual `=5·ν1(u−2h)ν1(u−3h)` exactly (a
pure constant shift), which pinned the missing term via the actual solvability condition, and
cross-checked against 27 irreps symbolically in `θ,h` (all exact) plus 4 more irreps with
`deg q1≥1` to break a degeneracy where `Cas2²` and `Cas4` weren't separable. **Lesson: build the
solver early — it's also the most direct correctness check.**

### Baxter polynomial solver + gauge fixing (in `baxter_SU4_L2_XXX.wb`)

`QsolveG[μ,S,P]` (note: `(μ,S,P)` order in *this* notebook — see naming-convention note below)
solves the untwisted 4th-order Baxter (difference) equation for the four monic Q-polynomials,
**gauge-fixed**: since `deg q1<...<deg q4` strictly, any lower-degree solution can be added to a
higher one (`q_k → q_k+c·q_j`, `j<k`); removed by forcing the coefficient of `u^(deg q_j)` in `q_k`
to vanish for every `j<k`. Verified exact (symbolic in `θ,h`) over full decompositions of every
`(S,P)` up to `(3,2)`, and numerically (rational `θ1=1/3,θ2=1/7,h=1`) over **all 21 irreps at (5,5)**
and **all 45 irreps at (8,8)** — every `q_k` solves the equation exactly, timing flat at
~0.02–0.05s/irrep **independent of irrep dimension** (largest block tested: `dim=69825`, same cost
as `dim=1`) since the solver only ever touches degree-≤19 polynomials, never the representation
itself.

Also built (Schur-polynomial machinery, reused verbatim in `GaudinSolver.wb` below): `irreps[S,P]`
— full `[S,S,0,0]⊗[P,P,0,0]` tensor-product decomposition via `s_λ·s_μ=Σ c^ν_{λμ} s_ν`, peeling off
the lex-highest monomial each step. Confirmed multiplicity-free throughout, dimension sums exact
(`Σdim(μ)=dim([S,S,0,0])²`), box-count (`Total[μ]=2S+2P`) conserved on every irrep as expected.

### Gaudin h→0 limit (derived in `baxter_SU4_L2_XXX.wb`, then built standalone)

Writing `Q=Q0+h·Q1+...`, orders `h⁰..h³` of the cleared Baxter equation vanish identically (same
alternating-binomial identity as above); order `h⁴` gives a **closed 4th-order ODE on Q0 alone**
(no `Q1,Q2` needed) — the gl(4) 2-site Gaudin oper, entirely in `Cas2(μ),Cas4(μ),S,P,θ1,θ2`.
Leading coefficient `8(u−θ1)²(u−θ2)²` (double singular points at the sites, as expected).
Validated: the ODE's polynomial-solution degrees match `degQ[μ,S,P]` exactly (and *only* those
degrees admit a solution) across every irrep of every `(S,P)` tested, symmetric and asymmetric.

**New standalone notebook `Paul/Mathematica/old/XXX/Experiments/GaudinSolver.wb`** — deliberately
minimal, gl(4) irrep/Casimir/degree machinery + the hardcoded `gaudinODE` only, **no** XXX transfer
matrices or Baxter difference-equation machinery (by design, per Paul's request). Built across
several rounds of naming/convention changes (all re-verified after each):
- Sites renamed `θ1,θ2 → w1,w2`; spectral variable renamed `u → z`.
- **Argument convention, applied everywhere: `(S,P,μ)`** — `degQ`, `GaudinSolve`, `BuildQs`,
  `verifyGaudin`, `CheckAllIrreps` all take `S,P` first, `μ` (or the irrep list) last. (Note: this
  is the *opposite* order from `QsolveG[μ,S,P]` in `baxter_SU4_L2_XXX.wb` — the two notebooks are
  not argument-compatible; don't copy-paste calls between them without swapping.)
- `GaudinSolve[S,P,μ]` returns **monic coefficient lists** (the practical/canonical form); `BuildQs
  [S,P,μ]` is a thin wrapper returning the same four Q's as **explicit polynomials in z** (kept
  purely for display/sharing, never used for internal verification/timing).
- `CheckAllIrreps[S,P,irrList]` takes a **precomputed** irrep list rather than recomputing
  `irreps[S,P]` internally — `irreps[8,8]` alone costs ~1.6s vs ~0.13s to actually solve all 45
  irreps, so precompute-once-reuse is the efficient pattern. Verified at `(8,8)` (45 irreps) and
  `(3,2)` (6 irreps).
- **Verification convention, hard-won:** always reconstruct a coefficient list into a testable
  function via `qf[x_]:=Sum[cf[[j+1]] x^j,{j,0,m}]` (SetDelayed), never
  `Dot[coeffs,z^Range[...]]` wrapped in a pure `Function` — the latter has a genuine WL edge case
  where symbolic `Derivative` of a length-1 (`deg=0`) Dot-built function spuriously returns a
  wrapped `{0}` instead of `0`, causing false verification failures. Cost a full debugging round
  before being isolated.

**Assistant arithmetic error, corrected:** an earlier summary wrongly claimed
`Σdim(μ)=8721²=76,055,441` for the `(8,8)` sweep — `8721` is `dimGL4[{16,16,0,0}]` (the *top*
output irrep), not `dimGL4[{8,8,0,0}]` (the actual input irrep, `=825`). Correct identity:
`Σdim(μ)=825²=680,625` — verified directly, matches the sweep. The sweep/solver results themselves
were never wrong, only that one squared-dimension claim (caught by Paul).

### Next natural steps

- **τ2 is now closed** — this directly unblocks the CG-overlap program the whole `z→1` limit was
  meant to feed (see the 2026-07-19/20 entries below for the original motivation). Natural next
  step: extract CG coefficients using the now-complete Casimir-form transfer matrices.
- Promote `baxter_SU4_L2_XXX.wb`'s Casimir-form τ1,τ2,τ3 + solver to `Clean/` once Paul is
  satisfied (still in `Experiments/`).
- `GaudinSolver.wb` is a clean, minimal, independently-reusable artifact — a natural next step is
  extending it toward the Gaudin-side CG-overlap analog, or documenting it via
  `wolfbook-docs-agent` once/if it's promoted to `Clean/`.
- Kernel-state loss happened repeatedly this session (multiple `KERNEL RESTART` events between
  turns) — always re-run cells top-to-bottom (`?symbolName` to check `DownValues` first) rather
  than trust that a previous turn's definitions are still live.

## State as of 2026-07-20 (gl(4) L=2 untwisted Casimir eigenvalues — τ1,τ3 DONE exactly, τ2 partial, in Experiments)

Follow-on to the twisted build below: new section appended to `Baxter_L2_XXX_SU4.wb` (now 70 cells,
35 onward), `z[k]→1` (untwisted/physical point), for the two-row-rectangular rep family
`λ1=[S,S,0,0]`, `λ2=[P,P,0,0]`. Goal: closed-form eigenvalues `Λ_a(u)` of the untwisted transfer
matrices `τ1,τ2,τ3` via su(4) Casimir invariants of the irreps `μ` in `λ1⊗λ2=⊕N_μμ`. Design/plan:
`docs/superpowers/specs/2026-07-19-su4-l2-untwisted-casimir-eigenvalues-design.md`,
`docs/superpowers/plans/2026-07-19-su4-l2-untwisted-casimir-eigenvalues.md` (the plan's original 5
tasks are DONE; everything below is follow-on work past the plan, done inline/live in the notebook
in direct collaboration with Paul, not via a formal plan). Subagents lacked wolfbook tools this
session (again) — a wolfbook-builder dispatch for the τ2/τ3 investigation partially worked (tools
were attached that one time) and its findings were independently re-verified from primitives.

**Exact-arithmetic simplification:** at `z=1` (vs. the twisted case's irrational `Zeta`-based twists)
every quantity is exact rational — no precision-cascade concerns anywhere in this section.

**Rep family is fully multiplicity-free** across the whole `(S,P)∈{(1,1),(2,1),(1,2),(2,2),(3,1)}`
sweep tested (up to `dim=400`), so every `τ_a` acts as a scalar on each irrep block (Schur's lemma) —
this is what makes closed forms in terms of Casimirs meaningful here at all.

### τ1 — DONE exactly

`Λ_1(u) = 4(u−θ1)(u−θ2) − 2P·h·(u−θ1) − 2S·h·(u−θ2) + (h²/2)(C2(μ) − C2(λ1) − C2(λ2))` (cell 48),
verified over the full sweep, `worst ~ 1e-25`. All `μ`-dependence enters only through the quadratic
Casimir eigenvalue `C2(μ)`, affine with slope `h²/2`. (This *resolved* an earlier false negative: a
naive fit of "each `u^j` coeff of `Λ_1` = `A_j+B_j·C2(μ)` **pooled across all `(S,P)`**" failed for
`j=0,1` — that was a pooling artifact, not a real falsification; the correct per-`(S,P)` formula
above is exact. See cells 45–47 for the superseded pooled-fit attempt, kept for the record with an
updated conclusion cell.)

### τ3 — DONE exactly

`Λ_3(u) = trivfactor[S,P]·( s3scalar[S,P](u)·Id + (h²/2)(C2(μ)−C2(λ1)−C2(λ2)) )` (cell 65), where
`trivfactor[S,P] = (u−θ1)(u−θ2)(u−θ1−h(S+2))(u−θ2−h(P+2))` (Paul's own factorization of τ3's four
trivial zeros, cells 58–63) and `s3scalar[S,P](u)` is an explicit state-independent degree-2
polynomial. **The `h²` term is exactly the same quadratic-Casimir combination as τ1's** — verified
as an exact operator identity on `(1,1),(2,1),(1,2)` including the `S↔P`-asymmetric case. So τ3
carries no higher Casimir, only `C2` — the `Λ³ℂ⁴≅(Λ¹ℂ⁴)*⊗det` self-duality made concrete.

### τ2 — PARTIAL (the self-dual middle case, genuinely harder)

Writing `Λ_2(u) = Σ_k[α_k(S,P)+β_k(S,P)·C2(μ)]u^k`, with `u⁰` additionally carrying
`+⅛C2(μ)²−¼C4(μ)` (universal coefficients, verified): **`α4=6,β4=0` (leading); `α3(S,P)=
−6h(S+P+2)−12(θ1+θ2), β3=0`; `β2=h²` (universal constant); `α2(S,P)` fully resolved (messier
polynomial, exact); `β1(S,P)=−(h³/2)(S+P+2)−h²(θ1+θ2)`** — all verified exactly against `tgCoeffs`
across 5 `(S,P)` points via direct operator-level matching (`M_k = α_k·Id+β_k·C2op`, solved exactly,
not curve-fit) — see cells 49–57 (the original agent-built investigation, confirming `τ2` uniquely
reaches `C4` only in its constant term) plus a from-scratch re-derivation in the live session after a
kernel restart (search for `solveABsym` in the transcript/notebook history).

**Still open: `α0(S,P)` and `α1(S,P)`'s full `(S,P)`-dependence.** These need cubic (`α1`) and
quartic (`α0`) polynomials in `(S,P)` — more unknowns than the 5 available `(S,P)` data points can
resolve. Getting more points is expensive (`(4,1)`→dim 630, `(3,2)`→dim 1000, symbolic-in-`θ,h`
computation at that size is impractical). An operator-level derivation was attempted instead: `T[i,j]`
was expanded exactly in terms of `Ee1,Ee2` building blocks (verified to reproduce `T` and `tg[2,1]`
bit-for-bit — see below) with the idea of reading off `α1,α0` structurally rather than by curve-fitting,
but this did not pay off within the session (a residual-decomposition test against `{X1,X2,Y1,Y2,C2}`
at a single `(S,P)` is underdetermined and not actually informative — would need to be done
simultaneously across multiple `(S,P)` to mean anything). **Promising unexplored lead, added by Paul
independently in cells 58–61**: `t1vac,t2vac,t3vac` — explicit vacuum-eigenvalue formulas built from
`L1=L2=(u−θ1−hS)(u−θ2−hP)`, `L3=L4=(u−θ1)(u−θ2)` — `t2vac` sums all 6 pairs `Li(u)Lj(u−h)` in a
specific descending-index order. An earlier attempt this session to derive the τ2 vacuum eigenvalue
via `Σ_{i<j}ν_i(u)ν_j(u−h)` (ascending-index order) did **not** match `tg[2,1]` on the HWS (nonzero
residual on pairs like `{1,3}`) — Paul's `t2vac` uses the opposite index order, which may be exactly
the fix; reconciling this against the actual `qmin`/`tg[2,1]` construction is the natural next step
and could resolve `α1,α0` cheaply (vacuum eigenvalues need no large-matrix computation at all).

### General Casimir operator formulas (not just eigenvalues — the operators themselves)

Cells 65–68 (note: cell 67, `X1/X2/Y1/Y2` definitions, was a **bugfix** — those operators were
originally only ever created via scratch `evaluateExpression` calls and never actually saved to a
notebook cell; a fresh-kernel run would have failed before the fix). General master formula (any `n`,
verified `n=2,3,4`): `Cn_op[λ1,λ2] = Σ_{pat∈{1,2}^n} patternTerm[pat,λ1,λ2]`, where each of the `2ⁿ`
site-assignment patterns collapses (via `(A⊗B)(C⊗D)=(AC)⊗(BD)` applied repeatedly) to `[ordered
product of site-1 factors]⊗[ordered product of site-2 factors]`, cyclic indices tying both chains
together. Reduced explicit forms with named cross terms:
```
C2op = (Csite[λ1]+Csite[λ2])Id + 2·Omega2                                    [1 cross term]
C3op = (C3site[λ1]+C3site[λ2])Id + 2X1+2X2+Y1+Y2                              [4 cross terms]
C4op = (C4site[λ1]+C4site[λ2])Id + 2Z1+2Z2+2Z3+2Z4+Z5+2Z6+Z7+2Z8               [8 cross terms]
```
Cross-term-type count grows 1→4→8 for `n=2,3,4` (matrix products aren't cyclic, only traces are, so
patterns that look like rotations of each other are usually genuinely different operators) — itself a
notable structural fact. Every formula verified as an **exact** operator identity (`==`, not
tolerance) on asymmetric `λ1≠λ2` pairs (avoiding a `λ1=λ2` coincidence hiding real distinctions), plus
`C3op`'s master formula also checked outside the `[S,S,0,0]` family (adjoint⊗fund) to confirm
generality. Plain operator definition (cell 70, self-contained): `C2=Σ_ab E_ab.E_ba`,
`C3=Σ_abc E_ab.E_bc.E_ca`, `C4=Σ_abcd E_ab.E_bc.E_cd.E_da`, `E_ab:=EE[λ1,λ2][a,b]`.

### General Casimir eigenvalue formulas (cell 69, self-contained, no other cell needed)

Single-site rep `[x,x,0,0]`: `Cₙ(x)=2x(x+2)ⁿ⁻¹` — `C2c[x]=2x(x+2)`, `C3c[x]=2x(x+2)²`,
`C4c[x]=2x(x+2)³`, verified `x=1..5` against the operators. Combined two-site:
`C2[S,P]=2S(S+2)+2P(P+2)`, `C3[S,P]=2S(S+2)²+2P(P+2)²`, `C4[S,P]=2S(S+2)³+2P(P+2)³`.

**Real bugs found and fixed this session** (beyond the earlier three): `X1/X2/Y1/Y2` never persisted
to a cell (see above); `T[i,j]` off-diagonal action on the HWS genuinely does vanish both directions
(`T[i,j][u].HWS=0` for `i≠j`, confirmed identically in symbolic `u`) — useful building block, not a
bug, but non-obvious and worth remembering; a hand-derived vacuum-formula attempt for τ2 had an
index/shift-order bug not yet isolated (see τ2 section above).

**Next natural steps:** reconcile Paul's `t2vac` (cell 61) against `tg[2,1]` on the HWS to get τ2's
vacuum eigenvalue cheaply, then use it (plus the already-known `β1,β0` and the operator-level `T[i,j]`
expansion) to pin down `α1(S,P),α0(S,P)` without more expensive `(S,P)` sweeps. Once τ2 is fully
closed-form, the natural follow-on is the CG-overlap program this whole `z→1` limit was meant to feed
into (see `HANDOFF.md`'s `su4 L2 untwisted` note in the twisted-case section below for the original
motivation).

## State as of 2026-07-19 (gl(4) L=2 twisted Baxter solver: `Baxter_L2_XXX_SU4.wb` — DONE, in Experiments)

Ported the two-site gl(3) XXX pipeline (`su3_V2_general.wb`) to **gl(4)**: two independent generic
Young-diagram reps `λ={λ1,λ2,λ3,λ4}` per site, four **distinct** diagonal-twist eigenvalues
`z[k]=N[Exp[I Zeta[2k+1]],70]` (k=1..4), 4th-order Baxter/TQ equation, four Q-functions, built fresh
in `Paul/Mathematica/old/XXX/Experiments/Baxter_L2_XXX_SU4.wb` (was empty). Design/plan (with deviation
notes recording what actually happened): `docs/superpowers/specs/2026-07-19-su4-l2-twisted-baxter-design.md`,
`docs/superpowers/plans/2026-07-19-su4-l2-twisted-baxter.md`. All 10 tasks DONE, executed inline
(not via subagent — the wolfbook MCP tools weren't wired into subagent tool sets this session, see
plan/session notes), controller commits.

**Companion twist dropped entirely** (user direction, confirmed before Task 2): no `G`/`TG`/`tG`/
`qminG` — only the diagonal-twist branch (`Tg`/`tg`/`qming`) is built, matching what the gl(3) source
pipeline actually used downstream anyway. `χ1..χ4` (elementary symmetric functions of `z[1..4]`)
still built directly, needed by `qdetT` and the leading-coeff check.

**Genuine gl(4) math (everything else is pure index/term extension from gl(3)):**
- `ν[λ1,λ2][i][u]=(u−θ1−hλ1[[i]])(u−θ2−hλ2[[i]])` for `i=1..4`; `qdetT=χ4·ν[1][u−3h]ν[2][u−2h]ν[3][u−h]ν[4][u]`.
- 4th-order `BaxOp`: `f[u+h]−A1 f[u]+A2 f[u−h]−A3 f[u−2h]+A4 f[u−3h]`.
- **Degree formula carries over unchanged**: `M_k=Λ_1−n_k` (`Λ=λ1+λ2`, `Λ_1=(λ1+λ2)[[1]]`) — verified
  with **zero deviation** against the derived asymptotic `Mdeg` master formula over the full 152-state
  sweep; no candidate-fit re-derivation was needed (the plan's fallback path was never triggered).
- `TauEigensystem4` simultaneously diagonalizes `τ1,τ2,τ3` (Cartan `H1,H2,H3`); `QSolve4` solves for
  four monic `q_k`, `Qfun`=`z[k]^(u/h) q_k(u)`.

**Full-sweep verification (acceptance):** 5 rep pairs, **152 states total**
(`{1,0,0,0}⊗{1,0,0,0}`=16, `{1,0,0,0}⊗{1,1,1,0}`=16, `{1,1,0,0}⊗{1,0,0,0}`=24,
`{1,1,0,0}⊗{1,1,0,0}`=36, `{2,1,1,0}⊗{1,0,0,0}`=60 adjoint⊗fund), every `Q_k` solves the 4th-order TQ
equation and the 5×4 Casoratian reconstruction of `A1..A4` matches `τ1,τ2,τ3,qdetT` —
**worstTQ = worstCasoratian = 0 exactly**. Verified by a fresh-kernel top-to-bottom run (restart → run
cells 1→34 in order); every hard-assert passed in-kernel with no `Abort`.

**Precision bug found and fixed (the real story of this session):** the plan specified `z[k]` at
30-digit precision (matching gl(3)); the initial full sweep at that precision left the d=60
adjoint⊗fund rep's Casoratian residual at `~4.5×10⁻⁸` — 45x over the `10⁻⁹` bar — even though the TQ
residual was *exactly* 0 everywhere (proof the formulas/implementation were correct; this was a pure
numerics issue). Bumping `z[k]` to 50 then 70 digits barely moved the residual (`~1.6×10⁻⁹`, then
bit-identical again) — the real giveaway that something else was capping precision independent of
`z`. Root cause: `TauEigensystem4`'s `rSets` (generic-combination coefficients for splitting the
degenerate τ-eigenspace) and its `Kmix` eigendecomposition were **hardcoded to 30-digit precision** —
a leftover copied verbatim from the gl(3) source's `TauEigensystem3`, never tied to the `z[k]`
precision bump. This silently capped the entire downstream chain (τ-coefficients, `Qfun`, Casoratian)
at ~19-25 accurate digits no matter how much precision `z` carried. Fix: match both to `z`'s 70-digit
precision — eigen-residual improved from a capped `8.6×10⁻²⁶` to `3.46×10⁻⁶⁴`, and the full sweep then
passed with **exact** zero residuals. Lesson for future gl(N) ports: hardcoded precision literals
copied from a smaller-N source notebook are a silent trap — grep for bare `N[..., <number>]` whenever
raising a global precision parameter.

**Next natural steps:** promote to `Clean/` once satisfied; the deliberate follow-on is the `z→1`
untwisted limit (su(4) irrep multiplet degeneracy, feeds the CG-overlap program) — explicitly out of
scope for this notebook per the design spec.

## State as of 2026-07-17 (gl(3) general reps per site: `su3_V2_general.wb` — DONE, in Experiments)

Generalized the two-site gl(3) XXX notebook from symmetric reps `{S,0,0}` per site to **arbitrary
gl(3) Young-diagram reps `λ={λ1,λ2,λ3}`** (`λ1≥λ2≥λ3≥0`, integers), one independent diagram per
site, in `Paul/Mathematica/old/XXX/Experiments/su3_V2_general.wb` (a copy of `su3_V1.wb` generalized in
place). Diagonal-twist frame, `(H1,H2,n)` labelling, three Baxter Q-functions, and the hard-assert
zero-residual convention all preserved. Design/plan (with Execution-outcome section):
`docs/superpowers/specs/2026-07-17-su3-general-reps-design.md`,
`docs/superpowers/plans/2026-07-17-su3-general-reps.md`. Tasks 1–6 all DONE; controller commits.

**Interface:** every function is re-signatured `S1,S2 → λ1,λ2` and takes full 3-vectors per site
(no scalar-`S` wrappers). Numeric params unchanged (`θ1=1/3, θ2=1/7, h=1`, high-precision `z[k]`).

**The three genuine math edits (everything else is pure relabeling):**
- **General Yangian weight** `ν[λ1,λ2][i][u] = (u−θ1−h·λ1[[i]])(u−θ2−h·λ2[[i]])` for `i=1,2,3`,
  from the Lax action `L[λ][i,i][u]|hw⟩=(u−h λ_i)|hw⟩` on the two-site HWS. Verified via the HWS
  eigenvalue check (T[i,i]|hw⟩=ν[i]|hw⟩) and the `qdetT == tg[·][3,1]` shift-ordering guard on
  non-symmetric reps. `qdetT[λ1,λ2][u]=χ3 ν[·][1][u−2h] ν[·][2][u−h] ν[·][3][u]` (V1 shift order
  survives).
- **Honest tensor-product `wtMult`:** `siteWeightMult[λ]=Counts` of the gl(3) weight of every GT
  pattern of `λ`; `wtMult[λ1,λ2][H1,H2]` recovers `(n1,n2,n3)` from `(H1,H2)` with `N=Total[λ1]+
  Total[λ2]` and convolves the two single-site tallies — correctly yields interior multiplicity `>1`.
  Sum rule `Σ wtMult == dim[λ1] dim[λ2]` hard-asserted across the sweep; guard on `{2,1,0}⊗{2,1,0}`
  (d=64) passes.
- **Degree closed form (Task 5 headline):** `M_k = Λ_1 − n_k`, with `Λ = λ1+λ2` the pseudovacuum
  highest weight and `Λ_1 = (λ1+λ2)[[1]]`. Sum rule `Σ_k M_k = 3Λ_1 − N`. Reduces to V1's `M_k = N−n_k`
  for symmetric reps (there `Λ_1 = N`). Extracted by fitting candidates `{N−n, Λ_1−n, Λ_k−n}` against
  the rep-agnostic `Mdeg` over the full sweep — only `Λ_1−n` gives zero deviation on the non-symmetric
  reps. Hard-asserted (`Mpred==Mdeg`, worst=0).

**Full-sweep verification (acceptance, Task 6):** six rep pairs, **102 states total**, every `Q_k`
solves the 3rd-order TQ equation and the 4×4 Casoratian reconstruction of `A1,A2,A3` matches
`τ1,τ2,qdetT` — **worstTQ = worstCasoratian = 0** (machine-zero, both `< 10^-9`). Per-rep counts:

| `λ1` | `λ2` | nstates |
|---|---|---|
| `{1,0,0}` | `{1,0,0}` | 9 |
| `{2,0,0}` | `{1,1,0}` | 18 |
| `{1,1,0}` | `{1,1,0}` | 9 |
| `{2,1,0}` | `{1,0,0}` | 24 |
| `{2,1,0}` | `{1,1,0}` | 24 |
| `{2,2,0}` | `{1,0,0}` | 18 |

Verified by a fresh-kernel top-to-bottom run (restart → run cells 1→40 in order); every hard-assert
check passed in-kernel with no `Abort`.

**wolfbook output-persistence glitch — re-run, do not trust blank boxes.** Several check cells created
via `insertCells` / edited during Tasks 3–5 do NOT refresh their on-disk `Out[]` display even after a
correct fresh-kernel re-run (confirmed: re-running an affected cell individually still shows blank).
These cells computed the correct result in-kernel (verified by direct `evaluateExpression`), but their
saved display is empty or stale. **A future reader must re-run the notebook rather than treat these as
failures.** Cells still showing empty/stale output on disk after the full Task-6 run (1-based cell
numbers; correct in-kernel values in parentheses):
- Cell 14 — `qdetT == tg[·][3,1]` guard: **empty** (in-kernel `{0}`).
- Cell 23 — `wtMult` sum-rule: **empty** (in-kernel all reps OK; the informational
  `adjoint⊗fund @(0,0)` spot value is `0` because for N=4 the (0,0) weight is non-integer — not a
  failure).
- Cell 24 — `{2,1,0}⊗{2,1,0}` d=64 guard: **empty** (in-kernel `tot=64=dim`).
- Cell 27 — eigensystem residual spot check on `{2,1,0}⊗{1,0,0}`: **empty** (in-kernel worst≈1.5×10⁻²⁷).
- Cell 31 — `degData` tabulation: **stale** (shows an old `{H1H2n→…,sumM→…}` format; source now emits
  raw `{λ1,λ2,H1,H2,n,{M}}` rows).
- Cell 32 — degree candidate-fit `results` Association: **stale** (shows an old matrix; source now
  emits `{"N-n"→…,"Lam1-n"→0,"Lamk-n"→…}`).
- Cell 40 — full-sweep report: **stale** (shows the old symmetric sweep, 141 states; the current
  λ-vector run gives the 102-state table above with worstTQ=worstCasoratian=0).

Check cells whose on-disk `Out[]` IS correct: cell 4 (`C1`→{0}), cell 10 (HWS→{0}), cell 21
(commuting-family, maxCommutator≈10⁻²⁶), cell 35 (degree assert worst→0; only its result-label string
still reads "N-n_k", cosmetic), cell 39 (single-state TQ perQresidual→{0,0,0}).

**Next natural steps:** promote `su3_V2_general.wb` to `Clean/` once satisfied; then the companion-twist
SoV / FSoV CG-overlap program for general gl(3) reps (mirroring the su(2) Sections C–F), and the
`paper/xxxCG.tex` write-up.

## State as of 2026-07-17 (gl(3) Section B: eigensystem + Baxter Q-functions — DONE, in Experiments)

Built the gl(3) analog of the su(2) Clean notebook's Section B, in
`Paul/Mathematica/old/XXX/Experiments/su3_V1.wb` (on top of the pre-existing Section A there). All
committed to `main`, all hard-asserts zero-residual. Design/plan (with Execution-outcome sections):
`docs/superpowers/specs/2026-07-17-su3-transfer-matrix-eigensystem-baxter-design.md`,
`docs/superpowers/plans/2026-07-17-su3-transfer-matrix-eigensystem-baxter.md`.

**What's there (Section B cells, diagonal-twist frame, numeric `θ1=1/3,θ2=1/7,h=1`,
`z[k]=N[Exp[I Zeta[3/5/7]],30]`):**
- `tgCoeffs` (coeff matrices of `τ1=tg[·][1,1]`, `τ2=tg[·][2,1]`), Cartan `H1op`/`H2op`, commuting-
  family check (`maxCommutator=0`).
- `wtMult[S1,S2][H1,H2]` — weight-multiplicity predictor (sum rule verified).
- `TauEigensystem3[S1,S2]` — simultaneous diagonalization; every state labeled `(H1,H2,n)`
  (`H1=E11−E22, H2=E22−E33`; `n` sorted by the full `τ1/τ2` charge tuple), sector counts match
  `wtMult`. Accessors `Psi`/`tau1`/`tau2`.
- `BaxOp` (TQ residual functional), `alphaCoeffs`, `Mdeg`, `qpoly`, `QSolve3` (three-Q solver),
  `Qfun` accessor `Q_k = z_k^{u/h} q_k(u)`.

**Key result — the degree relation:** `deg q_k = M_k = N − n_k` (`N=S1+S2`, `(n1,n2,n3)` the gl(3)
weight components), equivalently `M1=(2N−2H1−H2)/3, M2=(2N+H1−H2)/3, M3=(2N+H1+2H2)/3`. `M_k` pairs
with `z_k` in-index, depends only on `(H1,H2)`, `Σ M_k = 2N`. Hard-asserted.

**Verification:** full sweep `{(1,1),(2,1),(1,2),(2,2),(3,1),(1,3)}`, 141 states — every `Q_k`
solves the 3rd-order TQ equation (`worstTQ=0`) and the 4×4 Casoratian reconstruction of `A1,A2,A3`
matches `τ1,τ2,qdetT` (`worstCasoratian=0`). Commits `b4b68c7`…`645517d`.

**Next natural steps:** promote `su3_V1.wb` Section B to `Clean/` once satisfied; build the
companion-twist SoV covectors / FSoV scalar product for gl(3) (the CG-overlap program, mirroring
the su(2) Sections C–F); write up in `paper/xxxCG.tex`.

## State as of 2026-07-16 (SoV factorization + functional CG overlap — DONE, promoted to Clean)

Big result this session: the **Clebsch–Gordan overlap correspondence closed exactly**, then
promoted from `Experiments/` into `Clean/`. Everything below is **committed**. Design/plan docs
(with `Result`/`Execution outcome` sections recording what actually happened):
`docs/superpowers/specs/2026-07-16-xxx-sov-factorization-functional-overlap-design.md`,
`docs/superpowers/plans/2026-07-16-xxx-sov-factorization-functional-overlap.md`.

### The result

**The 1×1 FSoV determinant of the twisted Q-function equals the operatorial overlap of a site-1
unipotent rotation with the two-site eigenvector:**

    det1[λ1][ k^(u/h) Q1[λ1,λ2][M,n] ] = (⟨x0[λ1]| exp[-k E21] ⊗ ⟨x0[λ2]|) . Psi[λ1,λ2][M,n]

where `det1[λ1][f] := br1[λ1][f]/(f/.u->θ1)` (the 1-site FSoV determinant — Section C's `det`
restricted to one site). Dictionary: **generator E21, φ(k) = -k** (α=-1, universal). Residual is
the literal symbolic `0` across the rep sweep {(1,1),(2,1),(1,2),(2,2),(3,1),(1,3)} — 41 states.

This rests on a prerequisite factorization: the s2=0 slice of the L=2 SoV basis,
`x2[λ1,λ2][{s1,0}] = xL1[λ1][{s1}] ⊗ xSingle[λ2][{0}]` **exactly** (no leftover scalar), where
`xL1` is a genuine standalone length-1 SoV covector (new `t1`/`qdet1` transfer-matrix hierarchy).

**Normalization correction (same day, post-Clean-promotion):** the first pass of this factorization
check found an apparently nontrivial site-2 "passthrough eigenvalue" `μSite2[λ2,s1] =
Product[(θ1-θ2+(k-1)h-λ2 h),{k,1,s1}]`. **This was wrong** — a bug in `x2`'s own normalization, not
real physics, found by the user. `x2`'s denominator used `Qθ[u]=(u-θ1)(u-θ2)` evaluated at a point
shifted by `-λ1h`, which incorrectly puts that same `-λ1h` shift into the θ2 factor too (which
should carry `-λ2h`, site 2's own weight) — that mismatch produced exactly the spurious `μSite2`.
**Fix:** a new `ν1[λ1,λ2][u] := (u-θ1-hλ1)(u-θ2-hλ2)` (the Yangian weight eigenvalue of T11 on the
HWS; deduped with Section B's `aFun`, same formula) replaces `Qθ` in `x2`'s denominator, argument no
longer carrying the `-λh` shift separately. With the fix, the factorization ratio is **exactly 1**
for all 27 `(λ1,λ2,s1)` cases — proved analytically too: `ν1[λ1,λ2][θ1+(k-1)h]` factors as
`((k-1-λ1)h)·(θ1-θ2+(k-1-λ2)h)`, and that second factor is exactly the old spurious `μSite2`. `xL1`
needed no functional change (`Qθ1` is linear, fix is a numeric no-op) but was rewritten with a
single-site `ν1[λ][u]:=(u-θ1-hλ)` for parity. **Section F confirmed unaffected** (`dictReportDet`
identical before/after). Applied to both `Clean/XXX_CG_L2_Clean.wb` and
`Experiments/XXX_CG_L2_V2.wb`; fresh-kernel Run-All clean on both after the fix (also caught and
fixed one unrelated hidden-dependency bug in V2's `obstructionReport` cell, same class as the
`sts12` issue during Clean promotion — a scratch-only `ABTab` variable that was never a saved cell).

**The crux / what nearly went wrong:** the first functional-discovery pass used the bare bracket
`br1[k^(u/h) Q1]` (no determinant normalization) and did NOT close (k-dependent per-state scalar).
The fix was the `1/f(θ1) = 1/(k^(θ1/h) Q1(θ1))` normalization built into the FSoV determinant — it
cancels the `k^(θ1/h)` prefactor and the state-dependent `Q1(θ1)`, leaving `det1[k^(u/h) Q1]` a clean
degree-λ1 polynomial in k starting at 1. Caught by the user asking "was it normalized by
k^θ1/h Q1(θ1)?". If a functional-SoV object won't universalize, check for a dropped determinant
node-normalization first (see memory `reference_fsov_determinant_normalization`).

**Gotcha found this session:** `MatrixPower[nilpotentMatrix, 0]` fires `MatrixPower::sing` and stays
unevaluated even for a concrete integer power 0. Build `exp[φ E]` for nilpotent E as
`IdentityMatrix[dim] + Sum[φ^m/m! MatrixPower[E,m],{m,1,λ}]` (split off the identity), not `MatrixExp`
and not a Sum starting at m=0.

### Commits (chronological)

`737f6a2` spec → `1bd9bcd` plan → `8dbed52` Section E in `Experiments/XXX_CG_L2_V2.wb` → `9db89eb`
Section F in `Experiments/` + corrected docs → `dcf5c1e` **promoted Sections C,D,E,F into
`Clean/XXX_CG_L2_Clean.wb`** (now a contiguous 71-cell A→F notebook, fresh-kernel top-to-bottom run
verified, every hard-assert exact) → `bfc0ffb` refreshed `Clean/XXX_CG_L2_Clean.md`/`.tex` reference
docs to cover C–F (`.tex` compiles clean with plain `pdflatex`, `\mathbbm{1}` convention).

Section F was promoted in **final form only**: the corrected `det1`/`FLdet` machinery and the closed
dictionary. The exploratory bare-bracket `FL`, the pre-`det1` per-state-scalar cells, and the
`obstructionReport` were dropped from Clean (they remain in `Experiments/XXX_CG_L2_V2.wb` as the
derivation record, same convention as other Clean promotions).

## Open threads / next steps

- **Next up (2026-09-09):** build the perturbative-in-twist solver in
  `Experiments/SU4_twisted_Baxter_L2.wb` — keep `Q_j = Exp[phi_j u] x` (polynomial of fixed degree
  `M_j`) non-perturbatively, expand only the polynomial coefficients, as Laurent series with the
  now-known leading orders `eps^(-(M_k-M_1))`. Everything it needs (solver, degrees, `N_I`,
  divergence orders, exact-in-`eps` test states) is in place.
- Then: wire it to the FSoV norm (`SU4_norm_untwisting_L2_general.wb`), and scan over reps
  (item 3 of the 2026-09-02 plan, still untouched).
- Decide where the new gl(4) results get written up — `paper/xxxCG.tex` is gl(2)/single-site only,
  so a gl(4) L=2 section is a scope call; alternatively promote to `Clean/` and run
  `wolfbook-docs-agent`. Also still unwritten: L=2 factorization, scalar product, CG-overlap
  dictionary.
- Older, still open: extract the actual CG coefficients from the overlap generating function (the
  `exp[-k E21]` rotation generates the CG combinations — expanding in k gives per-state
  coefficients); push the factorization to the s2≠0 case / general L.
- `Experiments/XXX_CG_L2_V2.wb` still holds the full discovery narrative (bare-bracket attempt,
  obstruction report, per-state-scalar exploration) — useful if revisiting why the determinant
  normalization matters.
- Older restructuring/Clean-notebook state (L3 build, docs agent, shared wolfbook kernel gotcha) is
  unchanged from the 2026-07-15 session — see `git log` and `CLAUDE.md`.

## Stable structure

See `CLAUDE.md` for the full layout (XXX vs Gaudin split, Clean/Experiments convention, subagents).
`Clean/XXX_CG_L2_Clean.wb` is now the authoritative A–F notebook: A (SoV construction, incl. `ν1`) ·
B (Baxter/TQ, Psi/Q1) · C (2-site FSoV scalar product/orthogonality) · D (left eigenvectors,
`<Psi|Psi>=det[Q,Q]`) · E (exact L=1 factorization, c=1) · F (functional CG-overlap dictionary).
Companion `.md`/`.tex` reference docs are up to date with it (refreshed same day for the Section E
correction; `.tex` compiles clean with plain `pdflatex`).
