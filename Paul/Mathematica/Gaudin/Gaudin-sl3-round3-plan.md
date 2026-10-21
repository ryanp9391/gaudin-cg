# sl(3) Gaudin, round 3: no cleaning, the Gaudin B operator / jet basis, and [3,2,0] checks (2026-10-16)

Status: DONE 2026-10-16 — W3, P4, D3 returned and agree. Verdict blocks below §0 (P4, W3, D3). Results: `Experiments/Gaudin-L2-chain-320.wb`, `Paul/Python/gaudin_sl3/RESULTS-Qside-320.md`, `RESULTS-gaudin-B.md`. Next steps in HANDOFF.md. Rounds 1–2
(`Gaudin-sl3-CG-plan.md`, `Gaudin-sl3-round2-plan.md`) established the conventions, the localisers, the two-vacua untwist, and the
jet formula for the CG-complete two-vacua overlap of the adjoint. Paul's three questions this round:

## 0. The three questions and what is already known

**Q1 — can we avoid the per-sector / per-copy cleaning of XXX?** Yes, and it is already avoided. The cleaning (`new-conventions/
Clean-CG-combinations.md`) belonged to the single-vacuum + localised-T approach ⟨Ω|T_λ(θ₁−λ₁h)|Ψ⟩: the covector c_λ is pure top weight at
ε = 0, so its O(ε) lower-weight part meets the O(1/ε^k) raised corrections of the HW-limit state and gives O(1) contamination K/S = κ−1
(e.g. κ = (t+1)/2 for [3,3,0]); "sector-wide cleaning" = one combination of insertions removing K for a whole sector, "per-copy cleaning" =
copy-specific combinations where the sector clean space lacked HW rank ([4,3,1] in [3,1,0]², [5,3,2] in [3,2,0]², [6,4,2] in [4,2,0]² …).
The two-vacua family Ω_ω⊗Ω_ω′ removed the need for any of it already in XXX (Stage D: K/S = O(ε) for every HW branch, all CG data at
leading order), and in Gaudin the same holds with exactly ε-free vacua (W1 Task 4, W2 Task A: both reps, all branches incl. the
multiplicity-2 sectors), and the jet formula gives the CG data from the Q's alone with no cleaning step (D1 F4: Q-side vs chain to 12
digits, octets included). Not yet tested in Gaudin: the sectors that were worst in XXX — [5,3,2] in [3,2,0]² (this round, W3 Task A4),
multiplicity-3 sectors ([6,4,2] in [4,2,0]², needs dim 27 — later), and the weakly conditioned full-rank sectors (σ_min ≲ 0.05 in XXX; in
Gaudin the 2×2 systems of the CG extraction are solved from two grid points — W3 reports their conditioning for [3,2,0]).

**Q2 — the Gaudin B operator.** Not yet constructed. The XXX untwisted B(u) = T₁₁(u)·qmin_{12,12}(u+h) + T₂₁(u)·qmin_{12,13}(u+h) (gl(3), new
conventions; `ref/ref-L2-rotation-sov-all.txt` line 231 `Bop2`) has a simple spectrum ∏_s∏_k(u − θ_s − m_kh), labels (a|bc) per site with
roots θ_s + {−1−a, −b, −c}h, λ₂ ≤ a ≤ λ₁, 0 ≤ c ≤ λ₂, c ≤ b ≤ a, and its left eigenbasis is the SoV basis in which ⟨x|Ψ_c⟩/⟨x₀|Ψ_c⟩ =
N_x∏_s M^{(s)}_x/M^{(s)}_vac, M_{(a|bc)} = minor_s(−1−a,−b)Q³(θ_s−ch). As h → 0 all roots collapse to θ_s: B(u) = ν(u)³[1 + hb₁(u) + h²b₂(u) + …]
with b₁ Cartan-like (eigenvalue (1+a+b+c)/(u−θ_s) per site — degenerate), so the Gaudin "B operator" is not one operator with a simple
spectrum but a hierarchy whose joint generalised eigenbasis is a JET basis: the h-expansion of the SoV covectors ⟨x|(h) (rational in h) and
their divided differences. The jet formula of round 2 (five monomials N per site, Plücker p_kl and P³ jets) is what the SoV wavefunction
M_x becomes in that basis; the kernel G should be the Gaudin K tensor (Φ∞ expanded in the jet basis) times the wavefunction monomials.
Round 3 builds this (D3): B(u) exactly in h, the b_k, the SoV covectors as functions of h, the jet basis, the Gaudin wavefunction rule,
the Gaudin K tensor, and the derivation of G from it. Whether a single polynomial B_G(u) in the Gaudin Lax entries and derivatives has the
jet basis as generalised eigenbasis is the last, optional item.

**Q3 — [3,2,0] checks.** [3,2,0] (dim 15, 225 states) is the dual of [3,1,0] ([3,2,0] ≅ [3,1,0]*⊗det³), with λ₂ = 2: it tests the λ₂-
dependence of every rep-generic rule found so far. Predictions to check (all from rounds 1–2, nothing assumed beyond them):

| object | rule | [3,2,0] prediction |
|---|---|---|
| S₂ double pole (scalar) | λ₂(λ₁+1) (P3) | 8 |
| S₃ triple pole | none for λ₃ = 0 (P1, P3, W1) | none |
| dual-oper local exponents at θ_s | {0, λ₂+1, λ₁+2} (P3) | {0, 3, 5} ⇒ P′ = κP AND P″ = κ₂P for every solution; 5-jet matrix rank 2 |
| direct oper gauge | Π^{λ₁−λ₃}∘adjoint∘Π^{−(λ₁−λ₃)} (P1, P3, P2) | p = 3 |
| degrees | dual μ_j, direct 2λ₁ − μ_i | μ_j, 6 − μ_i; vacuum direct (0, 2λ₁−μ₂^{vac}, …) → report |
| Wronskians | direct F = Π^{λ₁−λ₂}, dual F = Π, W = det[P,P′,P″] ∝ e^{Σζu}Π^{\|Λ\|} | F = Π, Π, Π⁵ |
| det N | six Lagrange columns with Π^{−p} weights, rows (s,(i,j)) | closes on all 225 states |
| localisers | ℓ₁ = λ₁+1; {1, H₁, S3c[1,2]} ⊂ D₁ always; fourth from two quadratic products (W2) | ℓ₁ = 4; is X₄ = H₁² − ½H₁S3c[1,2] − (2/θ₁₂)W − S3c[1,1] again, or do the coefficients depend on Λ? D₂ relations S3c[s,1] = c_Λ·S3c[s,2](H_s − S3c[s,2]) with c_Λ = ¼ for [3,1,0] — find c_Λ |
| c_H | \|Λ\|Σζ + 2e₂(Λ)/θ₁₂ (W2, 10 reps incl. [3,2,0] — already verified) | 5Σζ + 12/θ₁₂ |
| spectators | f(H₁) = \|Λ\|Σζ − ζ·Λ + k/θ₁₂, f(S3c[s,2]) = \|Λ\|ζ₃; k = 4, 6 for [2,1,0], [3,1,0] ⇒ conjecture k = 2e₂(Λ) | k = 12; f(S3c) = 5ζ₃ |
| vacuum Laurent data | H_s^{vac} = Σ_i(\|Λ\|−λ_i)ω_i ± k/θ₁₂, a_s^{vac} = \|Λ\|ω₃, b_s^{vac} = ω₃(H_s − \|Λ\|ω₃) | 2ω₁ + 3ω₂ + 5ω₃ ± 12/θ₁₂, 5ω₃ |
| HW sectors of V⊗V | dual images of [3,1,0]²: μ* = (6−μ₃, 6−μ₂, 6−μ₁) | 10 HW states; multiplicity 2 in (5,4,1) and (5,3,2) |
| untwist | K/S = O(ε) every branch, both families; CG at leading order; mixing universal (S3c[1,2]\|_HW ∝ 1/θ₁₂ × fixed matrix) | holds; (5,3,2) — the XXX per-copy-cleaning sector — is the key test |
| jet formula | adjoint: 5 monomials/site from Plücker p_kl (l ≤ 3) × P³^{(c)} (c ≤ 1); M_vac = P³p₀₂ | ladder λ₁+2 = 5 points ⇒ jets to order 4, p_kl with l ≤ 4; lone c ≤ λ₂ = 2; p₀₁ = p₀₂ = p₁₂ = 0 ⇒ first nonzero p₀₃; M_vac = P³p₀₃?; monomial count to be found by the fit (P4, stretch) |

Numeric points as in rounds 1–2 (θ = (2/7, −3/5); ζ P, P′; ω, ω′; (ω₁, ω₁′) = (5/2, 7/4); regime φ, ω̂ grids; ε ∈ {1,2,4}·10⁻⁶).

**Verdicts from P4 (returned 2026-10-16, `RESULTS-Qside-320.md`) on the Q-side rows:** c₂ = 8 ✔; no triple pole ✔; dual exponents {0,3,5} ✔ and
direct exponents {0,2,5} (general: dual {0, λ₂+1, λ₁+2}, direct {0, λ₁−λ₂+1, λ₁+2}); gauge p = 3 ✔; degrees ✔; κ_s = a_s/c₂, κ₂,_s = (a² − Ha + 8b)/24;
direct F = Π ✔, dual F = Π^{λ₂} = Π² (table said Π — corrected), W ∝ Π⁵ ✔; det N ✔ on 225 states; vacuum H_s ✔ (k = 12), b_s ✔, a_s^{vac} = c₂ω₃ = 8ω₃
(table said |Λ|ω₃ — corrected; equal on the earlier reps). **Jet formula ✔ and generalised:** c·Qnum·(P³p₀₃)₁(P³p₀₃)₂ = Σ K (P³p_kl)(θ₁)(P³p_k′l′)(θ₂)
over the C(λ₁+2,2) − C(λ₂+1,2) non-forced Plücker coordinates of the order-(λ₁+1) jet (7 for [3,2,0], 9 for [3,1,0], 5 for the adjoint); the lone
factor is always just P³ (higher lone jets are redundant); M_vac = P³·p_{0,λ₂+1}; K unique, rational, polynomial in (ω₁, ω₁′) and in 1/θ₁₂ (degree 7 for
[3,2,0]); site exchange acts as θ₁₂ → −θ₁₂. Chain-side rows await W3; the B-operator/K-tensor derivation awaits D3.

**Verdicts from W3 (returned 2026-10-16, `Experiments/Gaudin-L2-chain-320.wb`) on the chain-side rows:** c₂ = 8 ✔; no triple pole ✔; ℓ₁ = ℓ₂ = 4 ✔ with
X₄ = H₁² − ½H₁S3c[1,2] − (1/θ₁₂)W − S3c[1,1] (W coefficient 1 here, 2 for [3,1,0]); the D₂ relation is an inhomogeneous quadratic, not c_Λ·S3c(H − S3c) (table
row corrected); c_H ✔; f(H₁) ✔ with k = 2e₂(Λ) = 12; f(S3c[s,2]) = c₂ζ₃ = 8ζ₃ (table said |Λ|ζ₃ — corrected, same correction as a_s^{vac} = c₂ω₃); vacuum H_s, b_s ✔;
κ^G = Δ₁Δ₂² ✔; HW sectors = duals of [3,1,0] ✔; universality ✔ ((5,3,2) has literally the same quadratic r² + 10r − 2 as [3,1,0]'s (4,3,1); mixing R identical).
**Q1 of this round answered for [3,2,0]: no cleaning anywhere** — K/S = O(ε) on all 10 branches in both families, (5,3,2) included (max|K/S| ≈ 2·10⁻⁶ at ε = 10⁻⁶),
CG to 10⁻¹⁸ after Richardson, mixing to 7·10⁻¹⁹; the only weakly conditioned 2×2 extraction ((4,4,2), (∞,∞) family, σ₂ ~ 3·10⁻⁴) is still accurate to 10⁻¹⁷ at
200 digits. The ω₁′-fixed control fails in the (∞,∞) family (as for [3,1,0]) and passes in the full family.

**Verdicts from D3 (returned 2026-10-16, `RESULTS-gaudin-B.md`) on Q2:** the Gaudin B is the hierarchy B/ν³ = 1 + hb₁ + h²b₂ + … with b₁ = 2L₁₁ + L₂₂ + ν′/ν
(diagonal, degenerate) and b₂ = 2C₁₁ + C₂₂ + L₁₁′ + L₁₁² + 2L₁₁L₂₂ − L₂₁L₁₂ + L₂₁L₂₃ + (ν′/ν)(2L₁₁+L₂₂) + 1/Π (needs the cross-site C_ij); the b_k do not commute
individually but are simultaneously triangularisable per level block (solvable algebra); no single polynomial B_G(u) of degree ≤ 3 in L, L′ does the job. The SoV
basis is the canonical jet basis defined by the wavefunction rule ⟨i,j|Ψ_c⟩/⟨x₀|Ψ_c⟩ = m_i(θ₁)m_j(θ₂)/(M₁M₂) with the eight one-site monomials m = (P³p₀₂, P³′p₀₂,
P³p₀₃, P³′p₀₃, P³′p₁₂, P³p₂₃, P³′p₁₃, P³′p₂₃) (one per site label, weights 2..6); its 64×64 matrix is exactly ζ-independent and rational; 59 of its 64 covectors couple
the sites through 1/θ₁₂. The leading-order degeneracy sits in the diagonal frame (P_phys⁻¹ = O(h)), not in B. The Gaudin K tensor of Φ∞ has exactly 25 nonzero
entries on N₅×N₅ and equals G/(4ω₁³ω₁′³) as an exact polynomial identity — G derived, D1 = D2 = D3.

## 1. Agents, files, kernel

| agent | type | files (owner) | kernel |
|---|---|---|---|
| **W3** | wolfbook-builder (sonnet) | `Experiments/Gaudin-L2-chain-320.wb` + spec `Gaudin-L2-chain-320-spec.md` | K1 (only user) |
| **D3** | general-purpose (sonnet), Python exact | `Paul/Python/gaudin_sl3/gaudin_B.py`, `RESULTS-gaudin-B.md`, caches `gaudin_B_*.pkl`; spec `SPEC-gaudin-B.md` | none |
| **P4** | general-purpose (sonnet), Python | `Paul/Python/gaudin_sl3/qside_320.py`, `twovac_fit_320.py`, `RESULTS-Qside-320.md`, caches `qside_320_*.pkl`, `twovac_fit_320_*.pkl`; spec `SPEC-Qside-320.md` | none |

All earlier Python files (`gaudin_chain.py`, `qside.py`, `xxx_chain.py`, `xxx_limit.py`, `qside_310.py`, `twovac_fit.py`, `twovac_limit.py`)
are read-only for everyone: import or copy, never edit. Nobody touches `paper/*.tex`, HANDOFF.md, CLAUDE.md, other notebooks, or git.

## 2. W3 — [3,2,0] chain side (notebook `Gaudin-L2-chain-320.wb`)

Copy the machinery of `Gaudin-L2-chain.wb` (round 1) and `Gaudin-L2-chain-2.wb` (round 2) with Λ = {3,2,0}: A1 Hamiltonians and Laurent
structure (c₂, triple pole); A2 vacuum, frames, κ^G; A3 eigenstates, HW sectors (compare with the dual images of [3,1,0]), branches, the
multiplicity-2 matrices S3c[1,2]|_HW at two θ (universality, mixing polynomials; compare with the [3,1,0] ones under duality); A4 the
two-vacua untwist for BOTH families (full ω̂, ω̂′ grid and the (∞,∞) family): exact undo, per-weight ranks, K/S for all 10 branches, CG at
leading order, the conditioning (second singular value) of every 2×2 extraction, controls — the cleaning question is answered by the
(5,3,2) and (5,4,1) rows; A5 localising scan: ℓ₁, ℓ₂, D₂ rank and its relations, the fourth localiser (test W2's X₄ formula; fit c_Λ and
the W coefficient), spectators; A6 closed forms: c_H, f(H₁), f(S3c), the vacuum Laurent data, against the table above; A7 reference tables
(225 Qnum(B) at P for the (∞,∞) family with eigenvalue labels, 30 digits; the HW branches at the regime point); A8 summary.

## 3. D3 — the Gaudin B operator, the jet basis, the Gaudin K tensor, and G from first principles (Python, exact)

B0 port and verify the XXX B(u) (adjoint, exact at h = 1/2: commutativity, degree 6, the 64 eigenvalues with the label rule); B1 the exact
h-expansion B = ν³[1 + hb₁ + h²b₂ + h³b₃ + …], the b_k in Gaudin Lax terms, their commutation relations, degeneracies and joint generalised
eigenspaces; B2 the SoV covectors ⟨x|(h) as exact rational functions of h (nullspaces at many rational h + rational interpolation, or
degenerate perturbation theory), their expansion at h = 0, the jet basis by divided differences, its structure (which labels merge, which
derivatives); B3 the Gaudin wavefunction rule: ⟨x_jet|Ψ_c⟩ as products of Plücker/P³ jets at θ₁, θ₂, verified on all 64 states with P1's
Q's; B4 the Gaudin K tensor (Φ∞ in the jet basis) and the derivation of G: compare with D2's 25 entries; B5 (optional) a single Gaudin
B_G(u) with the jet basis as generalised eigenbasis.

## 4. P4 — [3,2,0] Q side and the [3,2,0] jet formula (Python)

Q1–Q3 (Q4 if time) of P3's round for [3,2,0] (copy `qside_310.py` → `qside_320.py`, LAM = (3,2,0)); the Frobenius structure with two
forced jets; then the jet-formula fit (D1's `twovac_fit.py` machinery generalised: jets to order 4, Plücker k<l ≤ 4, lone c ≤ 2, M_vac
candidates P³p₀₃ …, κ^G for Λ = (3,2,0) = Δ₁^{1}Δ₂^{2}, x₀^G from P_aux^G for this rep) with a time budget; report the monomial count and the
kernel if found, or how far the hypothesis ladder got.

## 5. After this round

[3,1,0] jet formula (same λ₁ as [3,2,0]; if P4's fit works, run it for [3,1,0] too), the closed-form derivation of G from D3's K tensor
for general Λ, [4,2,0] (multiplicity 3) for the cleaning question, gl(4) rectangular Gaudin, the paper section, Clean promotion.

## 6. Risks and fallbacks

- **D3: the SoV covectors are not low-degree rational in h.** Then use degenerate perturbation theory order by order on the b_k, or Richardson
  numerics at 100 digits; the jet basis is defined by the span of the h-expansion coefficients, which exists regardless.
- **D3: the jet basis is not unique (unipotent gauge).** Expected; report the canonical choice used (divided differences in the label order,
  or the generalised eigenvectors of the b_k) and verify G is gauge-independent (it is a bilinear form on jets, not on the basis).
- **P4: the [3,2,0] fit is too large.** Reduce via the Frobenius relations (p₀₁ = p₀₂ = p₁₂ = 0, P³′ = κP³, P³″ = κ₂P³) to the ν-variables
  (Frobenius recursion coefficients), or fit per ζ-point after dividing by c(ζ); if still too large, stop at the rank/null-space diagnostics.
- **W3: [3,2,0] expectations refuted.** Report which; the duality with [3,1,0] is the first thing to check (a refuted rule that holds for
  [3,1,0] but not its dual points to a λ₂- or λ₃-dependence).
