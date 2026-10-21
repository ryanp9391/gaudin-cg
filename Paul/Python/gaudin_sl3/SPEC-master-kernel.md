# Spec M1: the master kernel R^{λ₁} in closed form — Pochhammer ansatz, [5,0,0] validation, the twisted sl(2) (round 5, WP-A1 + A3; 2026-10-16)

Delegated to a general-purpose subagent (SONNET), Python only (the Wolfram kernel belongs to agent W5 — you must NOT use Wolfram or any
wolfbook tool). Read first: `Paul/Mathematica/Gaudin/Gaudin-sl3-round5-plan.md` §0 (P2) and §1 WP-A (your goal: A1 and A3; A2 is agent M2's,
running concurrently — do not duplicate it), `Paul/Mathematica/Gaudin/Gaudin-sl3-round4-plan.md` §0 (notation — binding), then the
consolidated notes `Paul/Mathematica/Gaudin/Gaudin-sl3-CG-notes.md` §8 (the general kernel: conventions, one-site kernel, block theorem,
structure, first order, the single-index family, sum rule, the auxiliary sl₂) and Appendix B (all kernels), and `RESULTS-general-kernel.md`
(G1; K3/K4 hold every structural fact you build on). Environment: Python 3.13 (`python -X utf8`), sympy 1.14, mpmath 1.3, numpy; no `uv`;
run from `Paul/Python/gaudin_sl3/`. Four other agents run concurrently (M2, N1, Q1 in Python; W5 on Wolfram): keep to ≤ 2 heavy Python
processes of your own at a time, time one unit before any loop, cache to `master_kernel_*.pkl`.

**Code to reuse (read-only; import or copy, never edit):** `general_kernel.py` — `load_rep(lam)` → `(T, rep)` with `T[(i,j)] = {(n,p,q): Fraction}`
meaning Θ_ij = Σ c·tⁿ·aᵖ·a′^q (i, j index `pairs_of(lam)`; `e_of(lam,(k,l)) = k+l−λ₂−1`; `sym_entry`, `eval_T`), `ckl_closed(lam,k,l)` (one-site
closed form), `theta1_closed` (first order, all pairs), `F_mm`/`theta_single_closed` (the [l,l,0] Γ-form), `poch_fr`, `block_check`, `single_check`,
`corner_check`, `top_monos_exact`/`top_sum_rule_check` (sum rule), and the point builders `build_point2(lam, zeta, th, dps)`/`get_point2`/`site_monos`/
`overlap2`/`cfac`/`omega_inf`/`x_cov` (states + dual Q's + monomials at a base point; `general_kernel_bld.py` builds one point in a subprocess);
`twovac_fit.py` (`PivQR`, `rat_recon`, `mpf`), `qside.py`, `gaudin_chain.py`, `gl3_reps.py`. Exact kernels on disk: `general_kernel_rep_<λ₁λ₂λ₃>.pkl`
for 100, 110, 200, 210, 220, 300, 310, 320, 330, 400, 440 (polynomials in (a, a′, t), rational coefficients). **Files you own:** `master_kernel.py`
(stage flags; heavy stages separately), `master_kernel_*.pkl`, `RESULTS-master-kernel.md` (write it incrementally — a half-written report beats none).
No notebooks, no git, no Wolfram, no edits outside your files. If an expectation is refuted and two attempts do not fix it, stop that item, keep what
works, report what you saw.

## 0. The object (binding)

n := λ₁+1. Pairs (k,l), 0 ≤ k < l ≤ n, l ≥ λ₂+1. Block theorem (G1, exact λ₁ ≤ 3, fits λ₁ = 4): Θ^{[λ₁,λ₂,0]}_ij = [λ₂!(λ₁−λ₂)!/λ₁!]²(aa′)^{−λ₂}Θ^{[λ₁,0,0]}_ij on
the non-forced pairs, so everything is in the one-row kernels Θ^{[λ₁,0,0]} (exact for λ₁ = 1..4 in the pkl's). One-site kernel θ_kl(a) = c_kl a^{k+l−1} with
(λ₂ = 0) c_kl = (−1)^{k+l−1}C(λ₁, l−1)C(λ₁+1, k)(l−k)/l. **Define the master kernel R^{λ₁}_{(kl),(k′l′)} := Θ^{[λ₁,0,0]}_{(kl),(k′l′)}/(θ_kl(a)θ_k′l′(a′)).** By G1's
homogeneity (every monomial aᵖa′^q tⁿ of entry (i,j) has p+q−n = e_i+e_j) R is a sum over t-order N of terms x^m y^{m′}·F with x = at, y = a′t,
m+m′ = N, where m or m′ MAY BE NEGATIVE (G1's first order has a²/a′ and a′²/a: Θ₁ = tθ_pθ_q[αa + βa′ + γa²/a′ + εa′²/a], i.e. (m,m′) ∈ {(1,0),(0,1),(2,−1),(−1,2)}),
bounded below by p ≥ 0. First determine the exact support {(m,m′)} of R for each pair of pairs from the data (λ₁ ≤ 4) and state the rule. Known:
R(0,0) = 1 (Θ(t=0) = θ⊗θ′); the first-order coefficients (G1's α,β,γ,ε formulas with ρ(p) = (n−k)(n−l), σ = k+l−1, f₂ = 2(n−1)−σ, N₀ = n(n−1)); the
"extreme hops" (1/m!)(n−k)_m(n−l)_m(k′)_m(l′)_m/((n)_m)² (falling factorials; read G1's K3 for exactly which coefficients these are); the exchange rule
Θ_ij(a,a′,t) = Θ_ji(a′,a,−t); the rank pattern 1,2,4,8 of the t-expansion; the t-degree 4(λ₁−1) for λ₂ = 0; the single-index block (pairs (k, n)):
R_{(kn),(k′n)} = Σ C(l−k,m)C(l−k′,m′)F_{mm′}(k,k′)(x)^m(−y)^{m′} with F in Γ-form (G1 C3; the normalisation relative to the block theorem — check
G1's corner-block statement "= N(aa′)^{λ₁−λ₂} × the single-index form" — pin it down exactly first); the top-state sum rule Σ Θ_ij m_i^{top}m_j^{top} =
MM′R(a)R(a′), R(a) = (1−ζ₂a)^{λ₁−λ₂}(1−ζ₃a)^{λ₁}.

## 1. Tasks

**A1 — the four-index closed form.** (a) Tabulate R^{λ₁} exactly for λ₁ = 1..4 from the pkl's as `{(k,l,k′,l′,m,m′): Fraction}`; check the support rule,
the exchange rule in these variables, and that the (k,n),(k′,n) block equals the Γ-form. (b) Ansatz: F_{mm′}(k,l;k′,l′) = (sign)·(product of binomials/
Pochhammers in the six integers n,k,l,k′,l′,m,m′)/(product of the same type) — generalise the Γ-form (which has (k+m+1)_{m′−1}(k′+m′+1)_{m−1}(kk′+m′k+mk′) —
note the inhomogeneous factor (kk′+m′k+mk′): expect for general l,l′ a factor like a 2×2 determinant or a low-degree polynomial in the indices
times a Pochhammer ratio). Strategy that worked for G1: fix (m,m′), fit the dependence on (k,l,k′,l′) as a rational function of low degree with
the λ₁-dependence through n, using ALL four λ₁ at once (the same formula must hold for every n); check the λ₁-independence of the fitted rational
function first on the pairs common to several n. Prioritise the order-2 terms (m+m′ = 2 — complete them first, including the negative-m ones), then
order 3, then the general pattern. (c) Use the constraints as checks and as fitting equations: the sum rule with the top monomials (G1's
`top_monos_exact`), the exchange rule, the rank pattern (Θ_N has rank 2^N for N ≤ 3 then capped — if the N-th order is a sum of 2^N rank-1 terms,
find them: that is a product structure, likely the key), and the fact that det Θ^{[l,l,0]} has the explicit form in G1 K3. (d) Report the formula,
or the best partial one (which orders are closed, which remain tabulated).

**A1-val — predict [5,0,0] and test it on states.** dim V = 21, 441 two-site states, 21 pairs; a single base point ζ = P (G1's tag), θ = (2/7, −3/5),
ω = (5/2, 7/4) suffices because the test needs NO FIT: for each state the jet identity c·Qnum·M₁M₂ = Σ_ij Θ_ij m_i(θ₁)m_j(θ₂) is one scalar
equation, so evaluate Θ^{pred} from your formula (plus the block theorem's normalisation) at the base point and print the residual per state
(relative, ≥ 30 digits expected at 40-digit working precision). Build with `build_point2`-style code at 40 digits (G1 found 80 digits too slow for
λ₁ = 5): time ONE state's eigen-solve + dual-Q solve first; if the full 441 is more than ~1 h with 2 processes, use a random subset of ≥ 100
states — any subset is a valid test. If the formula is only partial, test the closed orders by subtracting the known part and fitting only the
rest (then you need the fit: 441 unknowns, use ≥ 2 ζ tags for overdetermination). Also test [5,1,0] (dim 35, 1225 states) ONLY if [5,0,0] passed
and time permits (block theorem prediction (1/25)(aa′)^{−1}Θ^{[5,0,0]}| on 20 pairs — a subset of ~60 states suffices).

**A3 — the twisted sl(2) and the fusion question.** (a) Make G1's observation precise: the vector θ(a) = (θ_kl(a))_{k<l} ∈ Λ²(C^{n+1}) (jet indices
0..n carry spin n/2 under an sl(2) with raising N: e_k → (something)e_{k+1}); find the matrix N (nilpotent, n+1 square) such that θ(a+s) =
exp(sN^{∧2})θ(a) exactly (check λ₁ ≤ 4; give N in closed form — presumably N_{k+1,k} = (n−k) or k+1 up to the c_kl normalisation). (b) The
two-site kernel: test Θ(a+s, a′+s, t′) = (exp(sN^{∧2})⊗exp(sN^{∧2}))Θ(a, a′, t) for t′ = t (G1: fails), t′ = t·f(s,a,a′) for the natural candidates
(the Möbius map on the conic: a → a+s is the unipotent e^{s(E₁₂+2E₂₃)}-type translation of the point (1,a,a²); conjugating the Gaudin model by
a GLOBAL unipotent g acts on θ₁₂? no — but it acts on the twist, so compare Θ at the SAME t and ask whether the mismatch is a known operator:
compute D(s) := Θ(a+s,a′+s,t) − (e^{sN}⊗e^{sN})Θ(a,a′,t) for λ₁ = 1, 2 symbolically and see whether it is of the form t·(commutator with a fixed
matrix) etc.). (c) Fusion: Θ^{[λ₁,0,0]} is the kernel of evaluating a bihomogeneous polynomial (V_{[λ₁,0,0]} = degree-λ₁ polynomials, E_ij = x_i∂_j,
Ωc∞[ω] = evaluation at the conic point (1, a, a²) up to normalisation — M2 verifies this; assume it) — equivalently the restriction of the
2λ₁-site FUNDAMENTAL Gaudin model with λ₁ coincident sites at θ₁ and λ₁ at θ₂ to Sym^{λ₁}⊗Sym^{λ₁}. For λ₁ = 1 the kernel is the R-matrix form
Θ = (1+z)θ⊗θ′ − zθ′⊗θ, z = t/(ω−ω′) = taa′/(a′−a). Test the hypothesis that R^{λ₁} is built from the λ₁ = 1 form by symmetrisation: write
the candidate as the image of ⊗^{λ₁}[(1+z)θ⊗θ′ − zθ′⊗θ] under the Sym projection in a basis where the pair indices (k,l) of [λ₁,0,0] arise from
λ₁ fundamental pair indices (0,1),(0,2),(1,2) (G1 refuted one "box-exchange" version — read exactly which, and try the symmetrised tensor
power with the coincident-site limit taken AFTER the product, i.e. with distinct points θ₁+δ_i first and δ → 0: the coincident limit of the
multi-site kernel involves derivatives in θ, which the naive product misses). Do λ₁ = 2 fully (36 entries) before generalising.

## 2. Reporting

`RESULTS-master-kernel.md`: §0 conventions (R's definition and normalisation, support rule), §1 the closed form (or the closed orders + tables of the
rest), §2 the [5,0,0] test (states used, residuals, timing), §3 A3 (N, the twisted-covariance verdict, the fusion verdict), surprises, not done.
Exact where possible; numbers with precision. Keep the pkl's small (coefficient tables only; the 441-state point is fine to cache).
