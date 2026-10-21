# sl(3) Gaudin, round 2: the Q-side formula for the CG-complete two-vacua family, and closure of [3,1,0] (2026-10-16)

Status: DONE 2026-10-16 — all four agents (W2 on K1; D1, D2, P3 in Python) returned and agree; the formula is in §1 (block added by D2), confirmed independently by D1 (block added by D1). Results: `Experiments/Gaudin-L2-chain-2.wb`, `Paul/Python/gaudin_sl3/RESULTS-twovac-fit.md`, `RESULTS-twovac-limit.md`, `RESULTS-Qside-310.md`. Next steps in HANDOFF.md. Round 1
(`Gaudin-sl3-CG-plan.md`, all three agents returned 2026-10-16) reproduced the XXX layers A–C for sl(3) Gaudin: Hamiltonians,
vacuum, Q-functions, localisers (every Laurent coefficient in D₁ localises; ℓ₁ = λ₁+1), and the two-vacua untwist with K/S = O(ε)
and all CG data at leading order, for [2,1,0] and [3,1,0]. **The project's goal is the Gaudin model; XXX was the scaffold.** What is
missing for Gaudin is layer D — the Q-side (bracket / jet) formula for the CG-complete two-vacua overlap — plus the loose ends of
[3,1,0]. gl(4) rectangular Gaudin (the actual target representation theory) is round 3: the structure of the formula (which jets,
which denominators, how the two sites couple) has to be seen once on sl(3) before it is worth building the gl(4) chain, and the
Python gl(4) generators do not exist yet (Mathematica has them in `old-conventions/old/Gaudin/Clean/Hamiltonians_SU4.wb`; XXX chain
templates in `new-conventions/L2-gl4-rect-untwist.wb`).

## 0. Round-1 facts this round builds on (all verified by ≥ 2 agents unless marked)

- cdet(∂ − L_g) column-ordered with L_ij = g_ij + Σ_s E^{(s)}_ij/(u−θ_s): S₂ = e₂(ζ) + Σ_s[H_s/(u−θ_s) + c₂/(u−θ_s)²] (c₂ = 3 for
  [2,1,0]), S₃ = e₃(ζ) + Σ_s[b_s/(u−θ_s) + a_s/(u−θ_s)²], no triple pole. Six non-scalar Laurent data per state.
- Dual oper (∂³ − S₁∂² + S₂∂ − S₃)P^j = 0, P^j = e^{+ζ_j u}p_j, deg p_j = μ_j (monic). Direct: Q_i = e^{−ζ_i u}q_i, deg 2λ₁−μ_i, with
  Q_i/Π² solving the formal adjoint, Π = (u−θ₁)(u−θ₂). Local exponents of the dual oper at θ_s: {0, 2, 4} ⇒ for EVERY solution
  P′(θ_s) = κ_s P(θ_s) with κ_s fixed by the Laurent data, and the order-3 jet (P, P′, P″, P‴)(θ_s) has only two free parameters.
- Wronskians: P^aP^b′ − P^bP^a′ = (ζ_b−ζ_a)e^{Σζu}Π(u)Q_c; so W(P²,P³)(θ_s) = 0 and W′(θ_s) = (ζ₃−ζ₂)e^{Σζθ_s}Π′(θ_s)Q₁(θ_s).
- Equal-twist overlap: ⟨Ω_ω⊗Ω_ω|Ψ_B⟩/n_B ∝ D_B(ω), the 6×6 det N of residue columns N[·,H_t], N[·,a_t], N[·,b_t] built from Π⁻²Q_i^{vac,ω}
  and P^j_B (P1), rows (s,(i,j)), (i,j) ∈ {(1,3),(1,2),(2,3)}; ratio identity on all 64 states to 6·10⁻⁵⁵.
- Frames: g′ = P_aux(ζ)g_c(ω)P_aux(ζ)⁻¹; P_phys^G(ζ) = (P_aux(ζ)⁻¹)ᵀ (W1); x₀^G = (e₁ᵀP_phys^{G −1})^{⊗2} is the h → 0 limit of the
  XXX x₀ (P2; ‖x₀‖ ~ h per site). Diagonal-frame vacuum Ω_D[ζ,ω] = e₁ᵀρ(n(ζ,ω)), r_k as in round 1, exactly ε-independent in the regime.
- Chain side of the family: Ω_D[ζ,ω]⊗Ω_D[ζ,ω′] untwists cleanly for every HW branch (W1 Task 4), octet mixing R² − 5R − 5 = 0.
- Open from round 1: [3,1,0] fourth localiser (needs ≥ 2 quadratic products), (4,3,1) mixing functionals, closed forms of the
  point constants, P1's rank-4 pairs, the gauge power of Π for [3,1,0].

**Added by P3 (returned 2026-10-16, `RESULTS-Qside-310.md`), rep-generic rules now verified on (2,1,0) and (3,1,0):** direct oper =
Π^{λ₁−λ₃}∘(formal adjoint)∘Π^{−(λ₁−λ₃)}; dual-oper local exponents at θ_s are {0, λ₂+1, λ₁+2}; S₂ double-pole scalar c₂ = λ₂(λ₁+1) (not scalar
for one-row reps); Wronskian F = Π^{λ₁−λ₂} (direct), Π (dual), W = det[P,P′,P″] ∝ e^{Σζu}Π^{|Λ|}; det N closes for (3,1,0) with the same six
Lagrange columns (ratio identity on all 225 states, 6.8·10⁻⁵²); the Laurent data of every eigenstate satisfy three no-log conditions per site
(e.g. a² − Ha + 4b = 0 for (3,1,0)) — the Gaudin analogue of Bethe equations; rank-4 pairs explained by the residue theorem (Δμ with one
positive and two negative components). Consequence for D1: the jets at θ_s of the dual Q's of the adjoint span the exponent-{0,2} Frobenius
solutions only up to order 3 (exponent 4 enters at order 4), as stated below.

**Added by W2 (returned 2026-10-16, `Experiments/Gaudin-L2-chain-2.wb`):** the (∞,∞) family IS CG-complete in Gaudin for both reps (K/S =
O(ε) all branches, CG to 1e-18..1e-21, universal mixing) — no fallback. **The Gaudin companion frame is regular and literally the XXX one
with w → ω:** Ωc^G[ω] = gaussCov{1/ω₁, −1/(ω₁ω₂), (ω₁+ω₂)/(ω₁ω₂)}, so Ωc∞ = gaussCov{a,0,a} with a = 1/ω₁; frame tie ΩD[ζ,ω]·PphysG[ζ] =
κ^G(ζ,ω)Ωc^G[ω], PphysG = ρ((PauxG⁻¹)ᵀ), κ^G = kappaTie with P₃ = (PauxG⁻¹)ᵀ. Consequence for the fit: the kernel is ζ-independent for the
companion-frame Q-number; the diagonal-frame Qnum carries the ζ-dependent scalar κ^G(ζ,ω₁)κ^G(ζ,ω₁′) (D1 informed). Expect the Gaudin
kernel to be the t → ∞ limit of the XXX kernels at fixed s = ω₁/ω₁′ with ladder → jet (D2 informed, optional L4). Control ω₁′ fixed FAILS
at ω₂ = ∞ (unlike round 1 at finite ω₂). [3,1,0] fourth localiser X₄ = H₁² − ½H₁S3c[1,2] − (2/θ₁₂)W − S3c[1,1]; D₂ has two exact relations
S3c[s,1] = ¼S3c[s,2](H_s − S3c[s,2]). Closed forms: c_H = |Λ|Σζ + 2e₂(Λ)/θ₁₂; f(S3c[s,2]) = |Λ|ζ₃; vacuum Laurent data ζ-independent
(H_s^{vac} = Σ_i(|Λ|−λ_i)ω_i ± k/θ₁₂, a_s^{vac} = |Λ|ω₃, b_s^{vac} = ω₃(H_s − |Λ|ω₃)); S3c[1,2]|_HW ∝ 1/θ₁₂ × a fixed matrix ⇒ mixing universal.

**Added by D2 (returned 2026-10-16, `RESULTS-twovac-limit.md`) — THE FORMULA (adjoint, L = 2), verified on all 64 states to 3·10⁻⁵⁵:**
4c·Qnum_G·(P³p₀₂)(θ₁)(P³p₀₂)(θ₂) = Σ_ij G_ij(ω₁,ω₁′,θ₁₂)N_i(θ₁)N_j(θ₂), N = (P³p₀₂, P³′p₀₂, P³p₀₃, P³′p₀₃, P³p₂₃) (five of the twelve
{P³^{(c)}p_kl} monomials), c = (ζ₂−ω₁)(ζ₃−ω₁)²(ζ₂−ω₁′)(ζ₃−ω₁′)²(ζ₁−ζ₂)² = κ^Gκ^G′, G rank 5, symmetric under site exchange, homogeneous in
(ω₁, ω₁′, 1/θ₁₂) of degree 10 − w_i − w_j with weights (2,3,3,4,5), unique (design rank 25). The limit is h¹⁰ (not h⁰: M_vac ~ h², κκ′/⟨x₀⟩ ~ h⁻⁶)
and all lower orders cancel by the single relation minor_s(−1,0) = 0 ⇔ Q^j(θ_s−h) = ρ_sQ^j(θ_s) (Baxter at θ_s with ν₃ = t₃ = 0). The XXX
Q's' own O(h) corrections do NOT enter (first at h¹¹). The naive fixed-a, t → ∞ limit of the kernels is wrong (L4). M_vac^G = 2h²W′P³.
This section's expectations are therefore confirmed in structure (jet-bilinear, state-independent kernel, M_vac ∝ P³Q₁) with the
corrections: only five monomials per site, kernel in (ω₁, ω₁′, 1/θ₁₂) with ζ only in the prefactor c. D1's independent fit is the cross-check.

**Added by D1 (returned 2026-10-16, `RESULTS-twovac-fit.md`):** the direct fit finds the same 5×5 kernel (unique in N₅, rational, equal to
the D2 kernel to 2.8·10⁻⁷⁴; every entry a polynomial of degree ≤ 6 in (ω₁, ω₁′) and ≤ 4 in 1/θ₁₂, no poles; homogeneity confirmed by a scaling
test); κ_s = a_s/3; the direct-Q₁ bases need a ζ-dependent kernel (σ = Σζ enters p₀₃ through Λ_s), so the Plücker form is canonical; the x₀-free
form is ⟨Φ∞|Ψ⟩·D_B(ω₀)/⟨Ω₀⊗Ω₀|Ψ⟩ ∝ Σ G_ij N_i N_j with no vacuum-monomial denominators; in the regime the Q-side formula reproduces all CG data
of adj⊗adj from the Q's alone (Richardson ratios, 12 digits vs the chain). Layer D is closed for the adjoint.

## 1. The target and the expected structure (binding for D1, D2, W2)

**Family.** The Gaudin (∞,∞) sub-family Φ∞(ω₁,ω₁′) := Ω_D^∞[ζ,ω₁] ⊗ Ω_D^∞[ζ,ω₁′] with Ω_D^∞[ζ,ω₁] := lim_{ω₂→∞} Ω_D[ζ,ω] =
e₁ᵀρ(n), n = exp(r₁E₁₂)exp(r₂E₁₃)exp(r₃E₂₃),

    r₁ = (ω₁−ζ₁)/(ω₁−ζ₂),   r₂ = −(ω₁−ζ₁)(ζ₂−ζ₃)/((ζ₁−ζ₂)(ω₁−ζ₃)),   r₃ = (ω₁−ζ₂)(ζ₁−ζ₃)/((ζ₁−ζ₂)(ω₁−ζ₃))

(the ω₂ → ∞ limits of the round-1 r_k; W2 verifies). In XXX the analogous w₂ = w₂′ = ∞ family was CG-complete at finite twist and in
the regime, with K/S = O(ε) for every branch, provided w₁ ≠ w₁′ (`Two-site-rotation-brackets-notes.md` §4.6, §5); its Q-side formula
(adjoint) is the ten-product bracket formula of §5 there. The Gaudin family depends on two parameters only; in the regime ω₁ = εω̂₁,
ω₁′ = εω̂₁′ the vacua are exactly ε-independent. W2 (Task A) verifies the chain side in Gaudin for both reps before anything is built on it.

**Normalisation-free targets** (pure functions of the Q's of the state B):
- (i) Qnum(B) := ⟨Φ∞|Ψ_B⟩/⟨x₀^G|Ψ_B⟩ (companion/x₀ normalisation, as in XXX);
- (ii) x₀-free: T(B) := ⟨Φ∞|Ψ_B⟩ · D_B(ω₀)/⟨Ω_D[ζ,ω₀]⊗Ω_D[ζ,ω₀]|Ψ_B⟩ with a fixed finite reference twist ω₀ and P1's det N D_B(ω₀).
  (T(B) = Qnum(B) × a state-independent constant × [D_B(ω₀)/the equal-twist Q-number] — the last factor is state-independent if P1's
  det N is the equal-twist Q-number, which the ratio identity says it is.)

**Ladder → jet dictionary (what the XXX formula becomes).** In XXX, ladder s has the four points θ_s + mh, m = −3..0, brackets
⟨Q^jF⟩_s = Σ_m F_mQ^j(θ_s+mh), 2×2 bracket determinants Bk2^{(s)}[F,G] = ⟨Q³F⟩⟨Q²G⟩ − ⟨Q²F⟩⟨Q³G⟩ = Σ_{m<m′}(F∧G)_{mm′}minor_s(m,m′),
the lone factors Q³(θ_s − ch), c ∈ {0,1}, the vacuum monomial M^{(s)}_vac = minor_s(−3,−2)Q³(θ_s−h), and the identity
⟨Φ_c|Ψ_c⟩/⟨x₀^c|Ψ_c⟩·M¹_vac M²_vac = Σ_{cc′}Q³(θ₁−ch)Q³(θ₂−c′h)Σ_{k=1}^{10}Bk2^{(1)}[F_k,G_k]Bk2^{(2)}[F′_k,G′_k] with F_m = a^m f_m(s,t),
G_m = a^{m+e}g_m(s,t), a = 1/w₁, s = a′/a, t = θ₁₂/h, kernels rational in (s,t) with poles only at t ∈ {0,±1}. As h → 0 (t → ∞, and
s = e^{−h(ω₁′−ω₁)} → 1, a → 1: the family degenerates to the global-rotation line at leading order, so the ω-dependence sits at O(h)):
- Q^j(θ_s+mh) = Σ_k (mh)^k P^{j(k)}(θ_s)/k! + (h-corrections of the XXX Q itself): the four ladder values ↔ the order-3 jet;
- minor_s(m,m′) ↔ Plücker coordinates p_{kl}(θ_s) := P²^{(k)}P³^{(l)} − P³^{(k)}P²^{(l)} (0 ≤ k < l ≤ 3), six per site, with h^{k+l};
  p₀₁ = W(θ_s) = 0, p₀₂ = W′(θ_s) ∝ Q₁(θ_s);
- M^{(s)}_vac = (5/2)h²·W′(θ_s)·P³(θ_s) + O(h³) (direct expansion: the O(h) term is −W(θ_s) = 0);
- Q³(θ_s − ch) ↔ P³(θ_s), P³′(θ_s) (c = 0,1).
So the Gaudin ansatz is: Qnum(B)·M¹M² (or T(B)·M¹M², or without denominators — D1 tests all) = bilinear form in degree-3 jet monomials,
per site {P³^{(c)}(θ_s)·p_{kl}(θ_s)} (12) or the smaller direct-side basis {Q₁^{(k)}(θ_s)P³^{(c)}(θ_s)} (6; because W-jets are Π·Q₁-jets),
with a state-independent kernel depending on (ω₁, ω₁′, θ₁₂) — by the scaling symmetry u → cu, g → g/c only through (ω₁θ₁₂, ω₁′θ₁₂),
and by the torus covariance of the (∞,∞) family presumably through still fewer combinations. The Frobenius constraint P′(θ_s) = κ_sP(θ_s)
and the oper recursion make the jet monomials linearly dependent in a state-dependent way, so the kernel is determined modulo a null
space: agents report rank, null space and a canonical (minimal-norm and a sparse) representative, as the XXX notebook did with pivot gauges.

**Two routes that must agree:** D1 fits the kernel directly from Gaudin data (states, Q's, overlaps — all in Python, independent of
XXX); D2 derives it as the formal h → 0 limit of the XXX ten-product formula in the diagonal frame (the companion frame degenerates:
P_aux(z) → rank 1 at h = 0). W2 supplies the chain-side verification of the family and reference overlap tables.

**Numeric points.** θ = (2/7, −3/5) (θ₁₂ = 31/35); ζ: P = (3/2, 5/4, 2/3), P′ = (2, 5/3, −1/2), extra ζ for the ζ-independence of
kernels: (7/5, −3/4, 1/6), (−4/3, 2/9, 5/7), (1/3, −2/7, 5/11); (ω₁, ω₁′) = (5/2, 7/4) base, (−2/5, 3) second, and a 3×3 grid
ω₁ ∈ {5/2, −2/5, 9/7} × ω₁′ ∈ {7/4, 3, −5/3} for the kernel laws (never ω₁ = ζ_j or ω₁ = ω₁′); θ variants (1/5, −4/9), (3/7, −1/4);
reference ω₀ = (5/2, 7/3, 3/5) for T(B). Regime: ζ = εφ, φ = (1/3, −2/7, 5/11), ω̂₁ ∈ {3/4, −5/3}, ω̂₁′ ∈ {−12/5, 7/2}, ε ∈ {1,2,4}·10⁻⁶;
chain 200 digits (W2), Q side 100–150 digits (D1).

## 2. Agents, files, kernel

| agent | type | files (owner) | kernel |
|---|---|---|---|
| **W2** | wolfbook-builder (sonnet) | `Paul/Mathematica/Gaudin/Experiments/Gaudin-L2-chain-2.wb` + spec `Gaudin-L2-chain-2-spec.md` | K1 (only user) |
| **D1** | general-purpose (sonnet), Python | `Paul/Python/gaudin_sl3/twovac_fit.py`, `RESULTS-twovac-fit.md`, caches `twovac_fit_*.pkl`; spec `SPEC-twovac-fit.md` | none |
| **D2** | general-purpose (sonnet), Python | `Paul/Python/gaudin_sl3/twovac_limit.py`, `RESULTS-twovac-limit.md`; spec `SPEC-twovac-limit.md`; reference `Paul/Mathematica/Gaudin/ref/xxx-210-infinf-formula-cells.txt` | none |
| **P3** | general-purpose (sonnet), Python | `Paul/Python/gaudin_sl3/qside_310.py`, `gaudin_chain_rep.py` (if a rep-generic copy is needed), `RESULTS-Qside-310.md`; spec `SPEC-Qside-310.md` | none |

Round-1 files (`gaudin_chain.py`, `qside.py`, `xxx_chain.py`, `xxx_limit.py`) are read-only for everyone: import, never edit. Nobody
touches `paper/*.tex`, HANDOFF.md, CLAUDE.md, other notebooks, or git. W2 uses the new-notebook workaround (JSON, evaluate in K1,
results as markdown cells, no native outputs).

## 3. W2 — chain side (notebook `Gaudin-L2-chain-2.wb`)

- **A. The Gaudin (∞,∞) family.** r_k^∞ closed forms (verify against ω₂ = 10⁶); the companion-frame vacuum Ω_c^G[ω] = e₁ᵀρ(upper Gauss
  factor of P_aux(ω)^{−T}) (W1's one-site rule applied to g_c(ω)), its eigencovector check for the one-site S_k^{(g_c(ω))}, its
  ω₂ → ∞ limit and Gauss parameters (is it e₁ᵀρ(e^{aE₁₂}e^{aE₂₃})-like, with which a(ω₁)?); the frame tie Ω_c^G[ω] ∥ Ω_D[ζ,ω]P_phys^G(ζ)
  with the scalar κ^G(ζ,ω) (XXX: kappaTie). Exact undo of Φ∞, closed forms of (a′, c′, g′) in (φ, ω̂₁, ω̂₁′), rank of the map, per-weight
  ranks (CG-complete?), K/S for all branches, CG at leading order, octet mixing — both reps; controls ω₁′ = ω₁ and ω₁′ fixed finite.
  Reference tables for D1/D2: Qnum(B) = ⟨Φ∞|Ψ_B⟩/⟨x₀^G|Ψ_B⟩ for all states at P (30 digits), and at the regime point for the branches.
- **B. [3,1,0] fourth localiser.** The exact localising subspace 𝓛 ⊂ D₂, its dimension, sparsest representatives outside D₁ (supports of
  size 2 and 3 among the quadratic monomials), natural candidates (H₁², H₁S3c[1,2], S3c[1,2]², S3c[1,1], H₁S3c[1,1], …), the fourth
  covector in terms of one-site operators and pair products, spectators, weights reached by the 3-dim vs the 4-dim L₁; site 2.
- **C. Closed forms / universality** (symbolic ζ, ω, θ): c_H(Λ, ζ, θ₁₂) in ⟨Ω⊗Ω|H₁ = [Ω₁(c_H − C_ζ)]⊗Ω₁ over several Λ (conjecture: the
  1/θ₁₂ coefficient is 2λ₁); the spectators f(H₁), f(S3c[1,2]) as functions of (ζ, ω, θ₁₂); the vacuum eigenvalue functions
  S₂^{vac}, S₃^{vac} in closed form; universality of the untwisted-copy matrix S3c[1,2]|_{HW} = (3/(2θ₁₂))[[3,−1],[−1,7]] (+ scalar?)
  across P, P′, two more θ, and the [3,1,0] (5,2,1), (4,3,1) analogues (hence whether R² − 5R − 5 = 0 is universal).

## 4. D1 — direct Gaudin fit of the two-vacua Q-number (Python)

F0 setup (states, Q's, x₀^G, Φ∞, Qnum and T); F1 jets at θ_s (Taylor coefficients to order 4, Plücker coordinates, Q₁ jets, the
Frobenius constraint P′ = κ_sP and κ_s vs the Laurent data, jet-matrix rank 2); F2 the hypothesis ladder (direct-side basis 6×6 → Plücker
basis 12×12 → larger; targets Qnum, T, with/without M¹M²; ≥ 5 ζ-points; rank, null space, residual, canonical representatives); F3 the
kernel laws in (ω₁, ω₁′, θ₁₂) (grid, scaling-symmetry check, rational reconstruction); F4 the regime and the Q-side CG (branches, exponents,
cancellation depth, Richardson CG and octet mixing vs W1's exact values); F5 compare with D2 if its results file exists.

## 5. D2 — the formal h → 0 limit of the XXX ten-product formula (Python, exact)

L0 parse the closed forms (f, g, f′, g′ per product, (e, e′) per block) and reproduce the printed kernel laws exactly; L1 assemble the XXX
identity in the diagonal frame (κκ′, x₀ normalisation, ladder points, a = e^{−hω₁}, s, t) and list every h-dependent ingredient including
the XXX Q's own h-corrections q^j = q^j_0 + hq^j_1 + …; L2 the formal Laurent expansion in h with symbolic jets: leading order, its
coefficient, verdict on whether q₁ jets enter; if not, the explicit Gaudin formula (Plücker and Q₁-jet forms, simplified with W(θ_s) = 0);
if yes, the equation for q₁ from the order-h⁴ Baxter expansion, and stop; L3 numerical test on the 64 Gaudin states (imports from P1/P2
code), and comparison with D1's kernel if available.

## 6. P3 — [3,1,0] Q side and the rank-4 puzzle (Python)

Q1–Q3 of round 1 for [3,1,0] (225 states): degrees 6−μ_i / μ_j, the gauge power of Π (0..4), Laurent structure (triple pole? c₂ value),
Wronskian F = Π^{λ₁−λ₂} = Π²?, Lagrange-identity columns and the det N ratio identity on all 225 states; Q4-type regime check on the 10
branches if time. Plus for [2,1,0]: classify P1's rank-4 pairs (which sectors), find and interpret the second null vector.

## 7. After this round

- **Round 3 = gl(4) rectangular Gaudin** (the project target): chain side from the XXX template `L2-gl4-rect-untwist.wb` with the Gaudin
  dictionary (cdet of a 4×4, S₁..S₄, Lie-algebra companion frame, vacuum r_k from the gl(4) Gauss chart), the localising scan, the
  two-vacua untwist (memory note `project-gl4-rect-untwist-why-clean`: abelian nilradical, κ/ω deformations), then the jet formula found
  here. Python needs gl(4) GT generators (port from `Hamiltonians_SU4.wb` or Matteo's `SU4_Paul_Normalisation_EXPLICIT.wl`).
- Paper: sl(3) Gaudin section of `paper/gaudinCG.tex` (conventions of round 1 §1, the localisers, the untwist, det N, the jet formula).
- Promote `Gaudin-L2-chain.wb` (+ `-2`) to Clean/ and run the docs agent once Paul is satisfied.

## 8. Risks and fallbacks

- **The jet ansatz does not fit (D1 F2 residual stays O(1) at every level).** Then the Gaudin overlap is not bilinear in jets of the ket's
  dual Q's at θ₁, θ₂ with a state-independent kernel. Diagnostics in order: allow the kernel to depend linearly on the six Laurent data of
  the state (if that fits, the "kernel" is itself a Q-number — rewrite it); include P¹ jets; include jets at a third point (u → ∞
  asymptotics of the p_j). Report the smallest ansatz that fits.
- **D2's leading order needs q₁ jets.** Expected possible (the family's ω-dependence is O(h) in the companion frame). Then D2 stops at the
  equation for q₁ and D1's direct fit is the result; D2's partial output still fixes which monomials and which h-orders carry the kernel.
- **The (∞,∞) family is not CG-complete in Gaudin (W2 A).** Then D1/D2 switch to the full Ω_D[ζ,ω]⊗Ω_D[ζ,ω′] (two-parameter ω̂, ω̂′
  grid as in W1 Task 4) with the same jet ansatz — the kernel then depends on four twist parameters; W2 must report this early (its Task A
  is first) and D1 must check W2's notebook results cell before F3.
- **[3,1,0] in P3: the gauge power of Π is not 2.** Expected (λ₁−λ₃ = 3 or λ₁ = 3); the Lagrange-identity columns change accordingly
  (Π^{−p} weights); P3 derives them for the power found rather than reusing P1's.
- Precision: the fit at finite twist is fine at 60–80 digits; the regime needs 100–150 digits on the Q side (the XXX formula lost ~37 digits
  to cancellations at ε = 10⁻⁶).
