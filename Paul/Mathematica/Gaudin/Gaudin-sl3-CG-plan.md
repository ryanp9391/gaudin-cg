# Plan: the sl(3) Gaudin model through the XXX "localising T's → CG" pipeline (2026-10-15)

Status: plan + three agent specs, dispatched 2026-10-15, re-dispatched 2026-10-16 after the first session was lost (W1 on the shared
kernel K1, P1 and P2 in Python so they do not compete for the kernel). **P1 and P2 returned 2026-10-16** (`Paul/Python/gaudin_sl3/
RESULTS-Qside.md`, `RESULTS-xxx-limit.md`); §1 below is updated with what they established (marked ✔ = verified by both, or by the
one agent named). **W1 returned 2026-10-16 as well** (`Paul/Mathematica/Gaudin/Experiments/Gaudin-L2-chain.wb`, 48 cells, no native outputs); all three agree
exactly on every overlapping item; the W1-specific results are in the rows marked (W1) and in §5. Nothing implemented before this date for sl(3) Gaudin; the only Gaudin material in the repo is
gl(4) rectangular (`old-conventions/old/Gaudin/Clean/*`, `paper/gaudinCG.tex`, Matteo's `Paul/Python/SU(4)_Paul_Normalisation_EXPLICIT.wb`).

## 0. Goal and what "repeat the XXX program" means

The XXX program (gl(3), L = 2, new conventions, `Paul/Mathematica/new-conventions/`) has four layers:

| layer | XXX result (where) |
|---|---|
| A. chain + vacuum + Q-functions | Lax u+hE, companion/diagonal frames, bra vacuum Ω_zw = e₁ᵀρ(n(z,w)), Baxter/dual Baxter, QQ, Wronskians, det N overlap (`Companion-B-C-diagonalize-clean.md`, `L2-diagonal-frame-overlap.md`) |
| B. localising T's | ⟨Ω⊗Ω\|T_λ(θ₁−λ₁h) = f_λ ⟨Ω T_λ^{(1)}\|⊗⟨Ω\| for λ = (1,1),(2,1) (adjoint); site 2: t₂(θ₂+h), t₂(θ₂+h)t₁(θ₂−h); the localising span L₁ = {Ω, ΩT₁₁, ΩT₂₁} inside the degree-≤2 transfer-matrix algebra D₂ (`L2-210-localising-span.md`) |
| C. CG-complete covector family + untwist | ⟨Ω_w⊗Ω_{w′}\| with z, w, w′ → 1 jointly: exact undo e₁ᵀρ(n(w)n(w′)⁻¹)⊗e₁ᵀ, all slices O(1), K/S = O(ε) for every HW branch, all CG data at leading order (`L2-two-vacua-untwist.wb`, `Two-site-rotation-brackets-notes.md` §4) |
| D. Q-side formula for the family | adjoint (∞,∞) family = 10 products of per-ladder 2×2 bracket determinants, kernels rational in t = θ₁₂/h (`L2-210-two-vacua-brackets.wb`) |

The Gaudin model is the h → 0 scaling limit of the same chain with twists z = e^{hζ}, w = e^{hω}. Every object above has
a Gaudin counterpart (§1). Layers A–C are reproduced in full for sl(3) Gaudin by this plan (adjoint [2,1,0]⊗[2,1,0] first,
then [3,1,0]); layer D is set up as the limit t → ∞ of the XXX closed forms plus a direct Gaudin test (§5).

**Guiding facts (established in the XXX work, used here):**
1. The vacuum Gauss parameters r_k(z,w) are homogeneous of degree 0 in (z−1, w−1), so the Gaudin vacuum is the SAME
   rational function r_k(ζ,ω), and in the untwisting regime ζ = εφ, ω = εω̂ the Gaudin vacuum covectors are EXACTLY
   ε-independent (in XXX they were O(1) with O(ε) corrections). The undone covector of the two-vacua family is likewise
   exactly ε-independent in Gaudin. Only the eigenstates depend on ε.
2. Gaudin has the scaling symmetry (u−θ) → (u−θ)/c, g → c·g. Hence the untwisting ε → 0 at fixed θ is the same as
   θ₁₂ → ∞ at fixed twist, i.e. the Gaudin model is "XXX at t = θ₁₂/h → ∞". Every structural XXX constant that is a
   rational function of t (κ_singlet = 2(t−3)/(3(t−2)) → 2/3, Δ_c, cc(t) …) has a Gaudin value = its t → ∞ limit.
3. The XXX companion matrix G_c(z) does NOT have the form exp(h g_c): G_c(1,1,1) is a Jordan block. The Gaudin
   companion frame is defined by the Lie-algebra companion matrix g_c(ζ) (§1); its P_aux(ζ) has the same closed form as the
   XXX P_aux(z) with z → ζ. **Corrected by P2 (X1):** the Jordan block is harmless — (G′−1)/h → g′ = P_aux(ζ)g_c(ω)P_aux(ζ)⁻¹
   at O(h), and the XXX x₀(z) at z = e^{hζ} IS parallel (to 10⁻³¹) to the Gaudin x₀^G(ζ) built from P_aux(ζ) by the same Gauss
   construction (‖x₀‖ ~ h per site). The x₀-free ratio form (§3) is still what the Q side uses, but x₀^G is a valid normalisation.

## 1. Conventions and the XXX ↔ Gaudin dictionary (binding for all three agents)

gl(3), L = 2 sites, rep Λ = (λ₁, λ₂, 0) on both sites; GT generators `JJ` (plain, exact rationals, non-Hermitian; the
invariant diagonal metric `gMet`), E⁽¹⁾ = E⊗1, E⁽²⁾ = 1⊗E in `KroneckerProduct` order; covectors are rows acted on from
the right; HW ket e₁; e₁ᵀE_ij ≠ 0 only for i ≤ j. Sites θ = (2/7, −3/5), θ₁₂ = 31/35. **No h anywhere in the Gaudin model.**

| object | XXX (new conventions) | Gaudin (this plan) |
|---|---|---|
| Lax | ℒ_ij(u) = u δ_ij + h E_ij, T = ℒ⁽¹⁾(u−θ₁)ℒ⁽²⁾(u−θ₂), twist on the right T·G | L_g(u)_ij = g_ij·1 + Σ_s E⁽ˢ⁾_ij/(u−θ_s) (T·G = ν(u)[1 + h L_g + O(h²)] with G = e^{hg}) |
| twist, diagonal frame | diag(z) | g = diag(ζ) |
| twist, companion frame | G_c(z) = [[χ₁,−χ₂,χ₃],[1,0,0],[0,1,0]] | g_c(ζ) = [[e₁,−e₂,e₃],[1,0,0],[0,1,0]], e_k = elementary symmetric polynomials of ζ |
| P_aux | rows (1, z_i−χ₁, z_i²−χ₁z_i+χ₂), P_aux G_c P_aux⁻¹ = diag(z) | rows (1, ζ_i−e₁, ζ_i²−e₁ζ_i+e₂), P_aux g_c P_aux⁻¹ = diag(ζ) |
| commuting family | t_a(u) = Σ_{\|I\|=a} qmin_{I,I}(u), rows antisymmetrised, column k at u−(k−1)h | ✔ cdet(∂_u − L_g(u)) = Σ_σ sgn σ · M_{σ(1)1}M_{σ(2)2}M_{σ(3)3}, M = ∂_u·1 − L_g, ∂ moved to the right (∂f = f∂ + f′) =: ∂³ − S₁(u)∂² + S₂(u)∂ − S₃(u). Explicitly (ordered products, l = L_ij, l′ = ∂_u l): S₂ = Σ_{i<j}(l_ii l_jj − l_ij l_ji) − l′₂₂ − 2l′₃₃; S₃ = l″₃₃ − l₁₁l′₃₃ − l₂₂l′₃₃ − l′₂₂l₃₃ + l₃₁l′₁₃ + l₃₂l′₂₃ + l′₃₂l₂₃ + (the six ordered 3×3-determinant terms). Row-ordered with E_ij and column-ordered with E_ji do NOT commute; row-ordered with E_ji is the transposed mirror. **Dictionary (P2, exact):** t₁/ν = 3 + hA₁ + h²A₂ + h³A₃, t₂/(ν(u)ν(u−h)) = 3 + hB₁ + h²B₂ + h³B₃ with A₁ = S₁, B₁ = 2S₁ (forced by the Baxter equations), S₂ = −A₂ + S₁²/2 + S₁′/2, S₃ = A₃ − B₃ + S₁³/6 + S₁S₁′/2 + S₁″/2, B₂ = A₂ + S₁²/2 + S₁′/2; S_k(u→∞) = e_k(ζ). |
| scalar members | t₃ = χ₃ ν₁ν₂ν₃ | ✔ S₁ = tr L_g = Σζ_k + \|Λ\|Σ_s 1/(u−θ_s) (scalar); **S₃ has NO triple pole** (its (u−θ_s)⁻³ coefficient is exactly 0 for [2,1,0]); S₂ = e₂(ζ) + Σ_s[H_s/(u−θ_s) + 3/(u−θ_s)²] (double pole the scalar 3), S₃ = e₃(ζ) + Σ_s[b_s/(u−θ_s) + a_s/(u−θ_s)²] |
| non-scalar coefficient operators ("the T's") | t₁: 2, t₂: 4 (u-coefficients) | S₂: H₁ (residue at θ₁; H₂ = −H₁ + scalar); S₃: the (u−θ_s)⁻² and (u−θ_s)⁻¹ coefficients, s = 1,2 (one linear relation from the residue theorem) — plus W = Σ_k ζ_k (E_tot)_kk (scalar per weight sector) |
| explicit forms ✔ (W1, exact; sign CORRECTED) | — | S₂(u) = e₂(ζ) + Σ_s [H_s/(u−θ_s) + 3/(u−θ_s)²] for [2,1,0]; **H₁ = −[Σ_i ζ_i E⁽¹⁾_ii + Σ_{ik}E⁽¹⁾_ikE⁽²⁾_ki/θ₁₂] + (\|Λ\|Σζ + \|Λ\|²/θ₁₂)·1**, H₁ + H₂ = −W + 2\|Λ\|Σζ·1 (the first version of this row had the opposite sign; the residue of Σ_{i<j}(l_iil_jj − l_ijl_ji) gives −C_ζ⁽¹⁾ + \|Λ\|Σζ and (\|Λ\|² − Σ_{ik}E⁽¹⁾_ikE⁽²⁾_ki)/θ₁₂) |
| bra vacuum | Ω_zw = e₁ᵀρ(n(z,w)), eigencovector of t_a with twist G′ = P_aux(z)G_c(w)P_aux(z)⁻¹ | Ω_{ζω} = e₁ᵀρ(n(ζ,ω)), n = exp(r₁E₁₂)exp(r₂E₁₃)exp(r₃E₂₃), r₁ = (ω₁−ζ₁)/(ω₁−ζ₂), r₂ = (ω₁−ζ₁)(ζ₁−ω₂)(ζ₂−ζ₃)/((ζ₁−ζ₂)(ω₁−ζ₃)(ω₂−ζ₃)), r₃ = (ω₁−ζ₂)(ζ₂−ω₂)(ζ₁−ζ₃)/((ζ₁−ζ₂)(ω₁−ζ₃)(ζ₃−ω₂)); ✔ exact left eigencovector of S_k with g′ = P_aux(ζ)g_c(ω)P_aux(ζ)⁻¹ (this orientation and no other of the six variants; P1 and P2 independently); vacuum eigenvalue functions for ω: S₂ = 131/15 + (12539/930)/(u−θ₁) + 3/(u−θ₁)² + (4139/930)/(u−θ₂) + 3/(u−θ₂)², S₃ double pole 9/5 (ω), 27/7 (ω′) |
| Q-functions | Q_i = z_i^{−u/h}q_i (direct, deg q_i = 2λ₁−μ_i); Q^j = z_j^{u/h}q^j (dual, deg q^j = μ_j−2λ₃) | Q_i = e^{−ζ_i u}q_i, Q^j = e^{ζ_j u}q^j, same degrees |
| Baxter equations | 3rd-order difference equations (clean .md rows 7, 8) | ✔ **dual** = (∂³ − S₁∂² + S₂∂ − S₃)Q^j = 0, Q^j = e^{+ζ_j u}q^j, deg q^j = μ_j (P1 Q1; P2 X4: the dual Baxter equation becomes exactly this at order h³). **direct** = the formal adjoint (−∂)³ − (−∂)²S₁ + (−∂)S₂ − S₃ acting NOT on Q_i but on Q_i/Π², Π = (u−θ₁)(u−θ₂) (i.e. Q_i solves Π²∘A∘Π⁻²); the plain adjoint has NO e^{−ζu}×polynomial solutions; gauge power 2 is the only one giving deg q_i = 2λ₁ − μ_i (shown for [2,1,0]; P1 and P2 independently). Quantum-Weyl (P1): with W = det[P,P′,P″] over the dual Q's, S₁ = det[P,P′,P‴]/W, S₂ = −det[P,P‴,P″]/W, S₃ = det[P‴,P′,P″]/W, W = ∏_{i<j}(ζ_j−ζ_i)e^{Σζu}Π³ |
| QQ relation | Q_aQ_b^{[2]} − Q_bQ_a^{[2]} = (1/z_b−1/z_a)χ₃^{−u/h}F Q^c, F = ∏_s(u−θ_s+2h) (adjoint) | ✔ (P1) Q_aQ_b′ − Q_bQ_a′ = (ζ_a−ζ_b)·e^{−Σζu}·(u−θ₁)(u−θ₂)·Q^c with monic q, p (constant exactly ζ_a−ζ_b); dual Wronskian: same F, constant (ζ_b−ζ_a) |
| brackets | ladder s: ⟨f⟩_s = Σ_{n=λ₃}^{λ₁} f(θ_s−nh)/∏_{m≠n}(m−n)h | ⟨f⟩_s = Res_{u=θ_s} f(u)/(u−θ_s)^{λ₁−λ₃+1} = f^{(λ₁−λ₃)}(θ_s)/(λ₁−λ₃)! (the divided-difference limit) |
| det N overlap | ⟨Ω(w)\|Ψ_B⟩/⟨x₀\|Ψ_B⟩ = D_B(w)P_{B₀}/(D_{B₀}(w)P_B(z₂−z₁)²), D_B = 6×6 det of per-ladder functionals in (ket dual Q^j, bra-vacuum direct Q_i) | ✔ (P1) x₀-free form holds on all 64 states (6·10⁻⁵⁵) and on the HW branches at ε = 10⁻⁶ (5·10⁻⁶⁷): rows (s,(i,j)), (i,j) ∈ {(1,3),(1,2),(2,3)}; columns from the Lagrange identity P(ΔS₂ψ)′ + ΔS₃Pψ = C′ with P = P^j (ket dual), ψ = Π⁻²Q_i (bra direct): N[·,H_t] = −Res_{θ_s}[ψP′/(u−θ_t)], N[·,a_t] = Res_{θ_s}[ψP/(u−θ_t)²], N[·,b_t] = Res_{θ_s}[ψP/(u−θ_t)] (t = 1,2) — six columns, Σ_m Δ_m N[·,m] = 0 with the state-DEPENDENT Laurent differences Δ (so no fixed-coefficient orthogonality exists; the naive per-site {⟨PQ⟩,⟨PQ′⟩,⟨P′Q⟩} columns do not close). Open: 12/60 different-sector pairs have rank 4 instead of 5 |
| localising T's | T₁₁, T₂₁ at θ₁−2h (adjoint); unified rule T_(k,k) (k ≤ d), T_(k,d) (d < k ≤ λ₁) at θ₁−λ₁h, d = λ₁−λ₂ | the Laurent coefficients of S₂, S₃ at θ₁ (and θ₂) and their degree-2 products: found by the rank-1 scan (W1 Task G3). **P2's prediction (X3, exact at the numeric point):** T₍₁,₁₎(θ₁−2h) = 3h²θ₁₂²·1 + h³X₁ + …, X₁ = −(4631/124)·1 + 3H₁ − S3c[1,2] ∈ D₁ (no products needed), ⟨Ω⊗Ω\|X₁ = u⊗Ω₁ with u = Ω₁·(431/62 − 3C_ζ − S3^{1-site}[−2]); one-site v₁(T₁₁) = Ω·(215/8 − (9/2)C_ζ − S3[−2]), v₁(T₂₁) = Ω·(−215/8 + (9/2)C_ζ) (no S3 term); so L₁ ⊇ Ω₁·span{1, C_ζ, S3^{1-site}[−2]} expected. Site 2: ⟨Ω\|t₂(θ₂+h)t₁(θ₂−h) has valuation h⁴ with leading term Ω·(3C_ζ − 41/4) NOT ∝ Ω. **✔ CONFIRMED by W1's exact scan (the answer to "which T's"):** for [2,1,0], ℓ₁ = ℓ₂ = 3 already inside D₁ — the site-1 localisers are **{1, H₁, S3c[1,2]} each separately** (⟨Ω⊗Ω\|H₁ = [Ω₁(\|Λ\|Σζ + 4/θ₁₂ − C_ζ)]⊗Ω₁, ⟨Ω⊗Ω\|S3c[1,2] = [Ω₁·S3^{1-site}[−2]]⊗Ω₁ exactly), L₁ = Ω₁·span{1, C_ζ, S3^{1-site}[−2]} (full rank on every weight of V); site 2: {1, W+H₁ = −H₂ + scalar, S3c[2,2]} with +C_ζ; spectators (v normalised to first component 1) f(H₁) = 6 + 4/θ₁₂ at P (23/6 + 4/θ₁₂ at P′), f(S3c[1,2]) = 2 (−3/2 at P′) — Laurent polynomials in θ₁₂ with ζ,ω-dependent constants; P2's X₁ is a combination of the three. **[3,1,0]: ℓ₁ = 4 = λ₁ + 1** (3 over D₁, 4 over D₁ + H₁·D₁ and over D₂ = rank 26); the fourth localiser needs ≥ 2 quadratic products and its explicit form is open; f(H₁) = 95/12 + 6/θ₁₂, f(S3c[1,2]) = 8/3 at P. Unlike XXX (where only the specific fused shapes T₁₁, T₂₁ localise), in Gaudin every Laurent coefficient in D₁ at θ₁ localises on its own |
| untwisting regime | z = 1+εφ, w = 1+εω̂, w′ = 1+εω̂′ | ζ = εφ, ω = εω̂, ω′ = εω̂′ (vacua exactly ε-independent) |

Numeric points: finite twist ζ = (3/2, 5/4, 2/3), ω = (5/2, 7/3, 3/5), ω′ = (7/4, −2/5, 9/7) (the XXX z, w, w′ numbers reused
as Lie-algebra eigenvalues; all distinct, no ω_i = ζ_j). Regime directions φ = (1/3, −2/7, 5/11) (control φ′ = (−2/5, 3/8, 7/11)),
ω̂ ∈ {3/4, −5/3} × {−1/5, 9/7, 5/3}, ω̂′ ∈ {−12/5, 7/2} × {1/8, −3/4} (the Stage D grid, 24 points), ε ∈ {1, 2, 4}·10⁻⁶,
200 digits on the chain. Generic points for operator identities: u₀ = 3/7, u₁ = −2/9.

## 2. Stage G (chain side, wolfbook, agent W1) — `Paul/Mathematica/Gaudin/Experiments/Gaudin-L2-chain.wb`

Spec: `Paul/Mathematica/Gaudin/Experiments/Gaudin-L2-chain-spec.md`. Rep-generic with Λ as the switch; [2,1,0] in full, then [3,1,0].

- **G0 setup.** Generators, metric, L_g(u) at numeric u (exact rationals), S₁, S₂, S₃ from cdet with the symbolic normal
  ordering done once on noncommuting symbols {L_ij, L′_ij, L″_ij}; checks: [S_k(u), S_l(v)] = 0 exactly; S₁ scalar; the
  expected explicit forms of S₂; weight conservation for diagonal g.
- **G1 vacuum and frames.** g_c(ζ), P_aux(ζ), g′(ζ,ω); Ω_{ζω} from the r_k; VERIFY Ω⊗Ω is the exact left eigencovector of
  S₂^{(g′)}(u), S₃^{(g′)}(u) (which orientation of g′ works); the one-site statement Ω_{ζω} = HW covector of the Cartan
  element g′ (e₁ᵀρ(M) with M g′ M⁻¹ = diag(ω)); Gaudin companion-frame x₀^G := (e₁ᵀP_phys(ζ)⁻¹)^{⊗2} with P_phys the
  rep image of P_aux(ζ) Gauss-factorised as in the XXX code (record it; optional for everything else).
- **G2 eigenstates and HW-limit branches.** Per weight sector, Eigensystem of S₂(u₀) + π S₃(u₁); at ζ = εφ the HW-limit
  branches (smallest raising residual); normalisation: G-norm 1 with sign convention; record ⟨Ω_{ω₀}⊗Ω_{ω₀}|Ψ_B⟩ for a
  reference ω₀ as the alternative normalisation for the Q side. Exact HW kets and G-orthogonal projectors as in Stage D.
- **G3 which T's localise (the answer to "which T's to use").** D₁ = span{1, W, H₁, S₃^{[1,2]}, S₃^{[1,1]}, S₃^{[2,2]}, S₃^{[2,1]}}
  (independent basis), D₂ = D₁ + quadratic monomials (rank ≤ 28). Exact rank-1 scan: X ∈ D₂ with ⟨Ω⊗Ω|X = v⊗Ω₁ → L₁ (and
  the site-2 L₂), at two finite-twist points. Identify the site-1 covectors: test L₁ ⊆ Ω₁·span{1, C_ζ := Σ_kζ_kE_kk, the
  one-site S₃^{1-site}(u) Laurent coefficients at θ₁ for twist diag(ζ), (g−g′)-type elements}. Report: ℓ₁, ℓ₂, which
  combinations X localise (these are the Gaudin "T_λ(θ₁−λ₁h)"), the spectator factors (analogue of f_λ), weights reached.
  Cross-check against P2's limit prediction (§4) when available.
- **G4 the two-vacua family and its untwist (the payoff).** Φ = Ω_{ζω}⊗Ω_{ζω′}; exact undo on all HW kets; regime; K/S for
  every HW branch at three ε's (expect O(ε)); CG at leading order from the overlap alone (1-dim slices and the 2×2
  systems for multiplicity-2 sectors), Richardson vs the exact HW-ket slices; the multiplicity mixing R_c with the
  Gaudin untwisted copies. Controls: ω′ fixed (not untwisted); common ω̂₁ = ω̂₁′.
- **G5 summary cell**, then [3,1,0] by changing Λ.

## 3. Stage Q (Q side, Python, agent P1) — `Paul/Python/gaudin_sl3/`

Spec: `Paul/Python/gaudin_sl3/SPEC-Qside.md`. Independent Python implementation of the Gaudin chain (reusing `gl3_reps.py`
for the GT generators; exact `Fraction` arithmetic for operators, `mpmath` for eigenstates), then:

- **Q1 opers.** Build S_k eigenvalue functions on every eigenstate; test which oper has polynomial×exponential solutions
  with the expected degrees (dual: e^{ζ_j u}, deg μ_j − 2λ₃; direct: e^{−ζ_i u}, deg 2λ₁ − μ_i); gauge-fixed polynomial
  solver (the `GaudinSolver.wb` pattern); fix the sign/role conventions and record them in the results file.
- **Q2 Wronskian/QQ and quantum-Weyl.** Q_aQ_b′ − Q_bQ_a′ ∝ F(u)Q^c (find F, the constant); S₂, S₃ eigenvalues as ratios of
  Wronskian-type determinants of the dual Q's (the Gaudin form of clean .md rows 10–11).
- **Q3 brackets, orthogonality, det N.** With ⟨f⟩_s = Res_{θ_s} f/(u−θ_s)^{λ₁+1}: find the functional orthogonality
  relation (which combinations of ⟨P^jQ_i^{(k)}⟩, ⟨P^{j(k)}Q_i⟩ vanish for A ≠ B — rank test), then the 6×6 det N and the
  x₀-free overlap-ratio identity on all 64 states with the vacuum direct Q's of twists ω, ω′ (state-independence to high
  precision). Optional: derive the same det N as the h → 0 limit of the XXX form and compare.
- **Q4 (if time).** The vacuum-overlap ratios at the regime point for the HW branches; the Q-side reading of the two-vacua
  overlap for the pieces that are det N ratios.

## 4. Stage X (the limit, Python, agent P2) — `Paul/Python/gaudin_sl3/`

Spec: `Paul/Python/gaudin_sl3/SPEC-xxx-limit.md`. A Python port of the XXX chain (one and two sites, exact rationals) at
z = e^{hζ}, w = e^{hω} for small h, to pin every convention by the limit rather than by assumption:

- **X1 frames and vacuum.** r_k(z,w) → r_k(ζ,ω); does P_aux(z)G_c(w)P_aux(z)⁻¹ = 1 + h g′ + O(h²) with g′ =
  P_aux(ζ)g_c(ω)P_aux(ζ)⁻¹ (or which variant)? The projective limit of the XXX x₀(z) vs x₀^G(ζ).
- **X2 Hamiltonian dictionary.** t₁(u)/ν(u) = 3 + hA₁ + h²A₂ + …, t₂/(ν(u)ν(u−h)) = 3 + hB₁ + h²B₂ + h³B₃ + …: express
  S₂, S₃ (cdet) as combinations of the A_k, B_k, products and c-numbers (exact linear algebra at a few u's). This fixes
  the cdet convention against the XXX quantum minors.
- **X3 localised covectors.** ⟨Ω_zw|T_λ^{(1)}(θ₁−λ₁h) (one site, cheap) at small h: expansion v₀, v₁, v₂ in h; the
  directions v₁ (first nontrivial) are the predicted Gaudin site-1 localised covectors — to be matched by W1's L₁. Same
  for the site-2 pair. Also the h-expansion of f_λ.
- **X4 Baxter → oper.** With X2's dictionary, expand the XXX direct and dual Baxter equations in h symbolically (orders
  h⁰..h² must cancel; order h³ is the oper); compare with P1's opers (sign/role of direct vs dual).

## 5. Later stages (not dispatched now)

- **Q-side bracket formula for ⟨Ω_ω⊗Ω_{ω′}|Ψ⟩** (layer D): (a) take t → ∞ of the XXX closed forms (`L2-210-two-vacua-brackets.wb`
  cell 62: 40 kernel laws κ(s,t), poles only at t ∈ {0,±1}) and translate ladder brackets into residue brackets; (b) direct
  Gaudin version via the Gaudin SoV basis (the B-operator limit; needs the Gaudin companion frame, G1's x₀^G) and the
  two-site K tensor (`L2-two-site-rotation-brackets.wb` pattern). Only after G4 and Q3 pass.
- **[3,1,0]**: G0–G5 DONE by W1 (K/S = O(ε) all branches, CG at leading order, (5,2,1) mixing 4.311737691/−0.811737691); OPEN: the
  explicit fourth localiser (≥ 2 quadratic products), the (4,3,1) mixing functionals (e₁E₁₃, e₁E₂₃E₁₂ do not suffice), and Q1–Q3 on
  the Python side (P1 fixed the Π² gauge for [2,1,0] only — check the gauge power for [3,1,0]).
- **Gaudin-specific surprises to understand** (W1 Task 4): the "ω′ fixed" control does NOT break the untwist in Gaudin (ΩD at ζ → 0
  has a finite direction-dependent limit, K/S stays O(ε), ranks full) — in XXX it did; octet mixing R_c = roots of R² − 5R − 5 = 0
  for [2,1,0] (compare the XXX value at finite t). Closed forms in (θ₁₂, ζ, ω) for the point constants (4631/124, 431/62, 1831/124).
- **P1 open points**: 12/60 different-sector pairs have det N rank 4 (not 5); why the same six equal-twist columns also close the
  mixed-twist overlaps (no ΔS₁/Δe₂/Δe₃ columns needed).
- **Paper:** `paper/gaudinCG.tex` is gl(4)-rectangular only; an sl(3) Gaudin section (conventions of §1, the localisers of
  G3, the untwist result of G4, the det N of Q3) once W1/P1 pass. Do not touch the paper in this round.
- **Promote** `Gaudin-L2-chain.wb` to a Clean/ notebook and run the docs agent once Paul is satisfied.

## 6. Agents, kernel, files

- **W1** = `wolfbook-builder` (sonnet) on K1, the only kernel. New-notebook workaround (memory `feedback-wolfbook-new-notebook-workaround`):
  build the `.wb` as JSON, open with the `code` CLI, evaluate each code cell's saved source in K1 via
  `wolfbook_evaluateExpression`, store printed results as markdown cells below the code cells; no native outputs (Paul
  reruns in VS Code to store them). Its first cell wipes `Global``.
- **P1, P2** = `general-purpose` (sonnet), Python 3.13 (`python -X utf8`), sympy 1.14, mpmath 1.3, numpy; no `uv` on PATH
  (plain `python`); scripts in `Paul/Python/gaudin_sl3/`, results in `RESULTS-*.md` there. File ownership: P1 owns
  `gaudin_chain.py`, `qside*.py`, `RESULTS-Qside.md`; P2 owns `xxx_chain.py`, `xxx_limit*.py`, `RESULTS-xxx-limit.md`;
  nobody edits the other's files.
- Reference dumps for all three (UTF-8, scratchpad of this session; copied to `Paul/Mathematica/Gaudin/ref/` for durability):
  `ref-L2-two-vacua-untwist-all.txt` (Stage D notebook: generators, transfer matrices, vacua, HW kets, branches, K/S, CG),
  `ref-L2-210-untwist-CG-all.txt`, `ref-Companion-B-C-diagonalize-clean-all.txt` (L = 1 Baxter/QQ/brackets/det N),
  `ref-L2-rotation-sov-all.txt`, `ref-L2-two-site-rotation-brackets-all.txt`, `ref-L2-210-two-vacua-brackets-all.txt`.
- No agent edits `paper/*.tex`, HANDOFF.md, CLAUDE.md, other notebooks, or git.

## 7. Risks and what to do if something fails

- **cdet convention wrong** (S_k do not commute, or Ω⊗Ω is not an eigencovector): W1 tries the transposed Lax
  (E_ji in place of E_ij) and the row-ordered determinant before anything else; P2's X2 settles it independently.
- **Bra-twist orientation**: four variants of g′ (listed in §1); if none works, compute the exact eigencovector of
  S₂^{(g′)} for g′ = P_aux(ζ)g_c(ω)P_aux(ζ)⁻¹ directly and compare with Ω_{ζω} — the eigencovector is what the later
  stages need; the closed form is secondary.
- **No localiser in D₂** (ℓ₁ = 1): report; G4 does not depend on G3 and still runs.
- **K/S not O(ε) for some branch in G4**: `showContam`-type breakdown, no fitting; this would contradict the XXX result at
  finite t and the t → ∞ argument, so first suspect the branch selection or precision.
- **det N (Q3) does not close**: try the derivative-column variants listed in the spec, then fall back to the h → 0 limit
  of the XXX functional orthogonality (P2's X4 machinery) before giving up.
- Exact 64×64 rational algebra is fine; 225×225 ([3,1,0]) needs sparse arrays and timing of one case before loops.
