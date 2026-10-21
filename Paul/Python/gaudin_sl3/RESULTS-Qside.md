# RESULTS-Qside: the Q side of the sl(3) Gaudin chain (gl(3), L = 2, Lambda = (2,1,0) on both sites) — agent P1, 2026-10-16

Independent Python implementation (no Wolfram). Files: `gaudin_chain.py` (model, exact sparse Fraction algebra, symbolic normal ordering,
eigenstates at 60/100 digits, bra vacuum), `qside.py` (Q1–Q4; `python -X utf8 qside.py`, about 75 s, prints everything), `qside_output.txt`
(the full printed output of the last run; all numbers below are copied from it), `qside_cache.pkl` (Q1 results cache; delete to recompute).
Precision: 60 digits for the finite-twist states (eigen-residual 1.3e-57), 100 digits in the regime (Q4); vacuum data exact rationals.
Not done: [3,1,0], Q3(d) (no x0^G available), the "moment" column variants (not needed, see 3).

## 0. Conventions fixed

* Rep: GL3Rep (2,1,0), dim 8; two-site E^(1) = E⊗1, E^(2) = 1⊗E (Kronecker order); matrices act on kets as columns, covectors from the right.
  `gl3_reps` only provides the adjacent generators; E_13 = [E_12,E_23], E_31 = [E_32,E_21] (checked [E_13,E_31] = E_11−E_33).
  **HW ket index (0-based) = 7**, GT pattern ((2),(2,1),(2,1,0)) — NOT index 0. It is annihilated by E_12, E_13, E_23 (checked).
  Weight sectors: 19, sizes 1,2,6,... up to 10 (the (2,2,2) sector).
* theta = (2/7, −3/5), g = diag(zeta), zeta = (3/2, 5/4, 2/3), omega = (5/2, 7/3, 3/5), omega' = (7/4, −2/5, 9/7).
* **cdet that works** (exact commutation): column-ordered product `cdet(d − L) = Σ_σ sgn σ M_{σ(1)1} M_{σ(2)2} M_{σ(3)3}`, M = d·1 − L, L_ij = g_ij + Σ_s E^(s)_ij/(u−θ_s)
  (E_ij, not E_ji), d moved right. The row-ordered product with E_ij and the column-ordered product with E_ji do NOT commute (274 nonzero entries in [S2,S3]);
  the doubly-flipped (row-ordered + E_ji) variant is the transposed chain and commutes. So the plan's placement is right. Symbolic result (l = L_ij(u), l′ = dL/du, l″ = d²L/du², ordered products):
  * S1 = l_11 + l_22 + l_33  (= tr g + 3/(u−θ1) + 3/(u−θ2), scalar, checked)
  * S2 = Σ_{i<j}(l_ii l_jj − l_ij l_ji) − l′_22 − 2 l′_33
    (literally: −l′_22 −2 l′_33 + l_11 l_22 + l_11 l_33 − l_21 l_12 + l_22 l_33 − l_31 l_13 − l_32 l_23)
  * S3 = + l″_33 − l_11 l′_33 − l_22 l′_33 − l′_22 l_33 + l_31 l′_13 + l_32 l′_23 + l′_32 l_23 + l_11 l_22 l_33 − l_11 l_32 l_23 − l_21 l_12 l_33 + l_21 l_32 l_13 + l_31 l_12 l_23 − l_31 l_22 l_13
* Checks (all exact rationals): [S2(3/7), S3(−2/9)] = 0; [S_k(u),S_l(v)] = 0 on a 3×3 grid of rational u,v; [S2,E_kk^tot] = [S3,E_kk^tot] = 0.
* Partial fractions (exact fit on 5 / 7 points, verified at 4 / 2 further points): S2 = e2(g) + Σ_s [H_s/x_s + 3/x_s²·1], S3 = e3(g) + Σ_s [b_s/x_s + a_s/x_s²], x_s = u−θ_s.
  So **S3 has NO triple pole** (the coefficient operators of (u−θ_s)^{-3} vanish identically); the double pole of S2 is the scalar 3 (Casimir-type), constants are e2(ζ) = 89/24, e3(ζ) = 5/4 (scalars).
  Verified explicit form: H_s = Σ_{i<j}[g_ii E_jj + g_jj E_ii − g_ij E_ji − g_ji E_ij]^(s) + (1/(θ_s−θ_s′)) Σ_{i<j}[E^(s)_ii E^(s′)_jj + E^(s)_jj E^(s′)_ii − E^(s)_ij E^(s′)_ji − E^(s)_ji E^(s′)_ij] (exact, both s).
  The non-scalar eigenvalue data of a state are thus (H_1, H_2, a_1, a_2, b_1, b_2): 6 numbers.
* Bra vacuum: Ω_{ζω} = e_1^T ρ(n), n = exp(r1 E12) exp(r2 E13) exp(r3 E23), r = (4/5, −7/11, 65/44) at (ζ,ω) (plan formulas). **g′ = P_aux(ζ) g_c(ω) P_aux(ζ)^{-1} (variant 0, first one tried) works**:
  Ω⊗Ω is an exact left eigencovector of all S2, S3 coefficient operators for both ω and ω′; the other five variants (P⁻¹g_cP, transposes, P^T g_c P^-T, P^-T g_c P^T) fail. The vacuum (exact) eigenvalue functions:
  * ω: S1 = 163/30 + 3/x1 + 3/x2; S2 = 131/15 + (12539/930)/x1 + 3/x1² + (4139/930)/x2 + 3/x2²; S3 = 7/2 + (2173/310)/x1 + (9/5)/x1² + (493/310)/x2 + (9/5)/x2²  (e2(ω) = 131/15, e3(ω) = 7/2)
  * ω′: S1 = 369/140 + 3/x1 + 3/x2; S2 = 29/28 + (40463/4340)/x1 + 3/x1² + (1263/4340)/x2 + 3/x2²; S3 = −9/10 + (30501/4340)/x1 + (27/7)/x1² + (−19899/4340)/x2 + (27/7)/x2²
  (the double-pole coefficient of S3 is state-dependent: 9/5 vs 27/7.)
* **Which ODE is direct/dual (THE convention).**
  * dual: (D) `[∂³ − S1 ∂² + S2 ∂ − S3] P = 0`, solutions **P^j = e^{+ζ_j u} p_j, deg p_j = μ_j** (μ = weight of the state; 2λ3 = 0).
  * direct: the formal adjoint (A) `−ψ‴ − (S1ψ)″ − (S2ψ)′ − S3ψ = 0` has NO e^{−ζu}×polynomial solutions (tested gauges Π⁰, Π⁻¹: none up to degree 8). The correct statement is that its solutions are ψ_i = Π(u)^{−2} Q_i
    with **Q_i = e^{−ζ_i u} q_i, deg q_i = 4 − μ_i = 2λ1 − μ_i**, Π(u) = (u−θ1)(u−θ2). I.e. the direct oper is Π² ∘ A ∘ Π⁻²; the gauge power 2 (= λ1 − λ3 here; for other reps unchecked) is forced by the degrees: Π^{-3} gives degrees raised by 2 per exponential
    ((2,4,6) at μ = (0,2,4)); Π^{-2} is the unique gauge with all degrees 4−μ_i. The plain adjoint solutions W_jk/W = Wronskian/(e^{Σζu}Π³) are not polynomial×exponential, which is why the spec's expectation "(A) itself" is refuted.
  * Exponent sign: dual +ζ, direct −ζ (wrong-sign controls: no solutions to degree 8). The same holds for the vacuum with ζ → ω (exact): direct Q_i^{vac,ω} = e^{−ω_i u} q_i, degrees **(0,2,4)**, q_1 = 1, q_2 = u² − (409/35)u + 2448/35, ...; dual degrees (4,2,0).
  * Normalisation: all p_j, q_i monic. Kernel dimension is exactly 1 for each of the 6 exponentials of each of the 64 states (no gauge freedom between different exponentials; degrees need not be strictly ordered).

## 1. Q1 — opers and Q-functions (60 digits)

All 64 states (19 sectors, eigenproblem S2(3/7) + (π/3) S3(−2/9) per sector, residual of all coefficient operators 1.3e-57): degrees of the six Q's as expected in every sector; kernel dimension 1;
max ODE residual (u = 0.37) 3.9e-55; largest smallest-singular-value at the degree found 2.6e-56 (smallest singular value at one degree lower is O(1e-3..1): clear gap).

| sector μ | #states | deg p_j (j=1,2,3) | deg q_i (i=1,2,3) = 4−μ_i | max residual |
|---|---|---|---|---|
| (0,2,4) | 1 | (0,2,4) | (4,2,0) | 5.4e-57 |
| (0,3,3) | 2 | (0,3,3) | (4,1,1) | 2.6e-56 |
| (0,4,2) | 1 | (0,4,2) | (4,0,2) | 1.1e-56 |
| (1,1,4) | 2 | (1,1,4) | (3,3,0) | 2.3e-56 |
| (1,2,3) | 6 | (1,2,3) | (3,2,1) | 4.3e-56 |
| (1,3,2) | 6 | (1,3,2) | (3,1,2) | 9.4e-56 |
| (1,4,1) | 2 | (1,4,1) | (3,0,3) | 2.1e-56 |
| (2,0,4) | 1 | (2,0,4) | (2,4,0) | 6.6e-57 |
| (2,1,3) | 6 | (2,1,3) | (2,3,1) | 9.1e-56 |
| (2,2,2) | 10 | (2,2,2) | (2,2,2) | 1.6e-55 |
| (2,3,1) | 6 | (2,3,1) | (2,1,3) | 1.5e-55 |
| (2,4,0) | 1 | (2,4,0) | (2,0,4) | 8.8e-57 |
| (3,0,3) | 2 | (3,0,3) | (1,4,1) | 3.9e-55 |
| (3,1,2) | 6 | (3,1,2) | (1,3,2) | 1.0e-55 |
| (3,2,1) | 6 | (3,2,1) | (1,2,3) | 1.0e-55 |
| (3,3,0) | 2 | (3,3,0) | (1,1,4) | 1.2e-56 |
| (4,0,2) | 1 | (4,0,2) | (0,4,2) | 6.4e-57 |
| (4,1,1) | 2 | (4,1,1) | (0,3,3) | 7.5e-57 |
| (4,2,0) | 1 | (4,2,0) | (0,2,4) | 6.2e-57 |

Negative controls (one state per sector): dual ODE with e^{−ζu}: none; plain adjoint (gauge 0 or Π⁻¹) with e^{−ζu}: none; Π⁻² gauge with e^{+ζu}: none (all to degree 8).
Regime caveat found in Q4: the leading coefficient of q can be ε-suppressed, so degrees must be read as "smallest d with a kernel" and not by trimming small leading coefficients (fixed).

## 2. Q2 — Wronskian (QQ) relation and quantum-Weyl formulas

(a) **Direct Wronskian**: Q_aQ_b′ − Q_bQ_a′ = (ζ_a − ζ_b) · e^{−(ζ_1+ζ_2+ζ_3)u} · F(u) · P^c (c ∉ {a,b}), all q, p monic, with
**F(u) = (u−θ1)(u−θ2)** (= ∏_s (u−θ_s)^{λ1−λ2}, state-independent). Over 64 states × 3 pairs: relative remainder of the polynomial division ≤ 3.6e-54, |F − (u−θ1)(u−θ2)| ≤ 3.4e-54, c_ab/(ζ_a−ζ_b) = 1 exactly (monic normalisation).
Wronskian degree − deg p_c = 2 always. (XXX analogue: (1/z_b − 1/z_a) χ3^{−u/h} F with z = e^{hζ}: (1/z_b − 1/z_a) ≈ h(ζ_a − ζ_b), χ3^{−u/h} → e^{−Σζ u}: consistent.)
**Dual Wronskian**: P^aP^b′ − P^bP^a′ = (ζ_b − ζ_a) e^{+(ζ_1+ζ_2+ζ_3)u} F(u) Q_c with the same F, remainders ≤ 3.0e-54, |F′−Π| ≤ 2.9e-54.
(b) **Quantum-Weyl**: S1, S2, S3 are exactly the coefficients of the monic operator annihilating P^1,P^2,P^3: with W = det[P^j,P^j′,P^j″],
S1 = det[P,P′,P‴]/W = (ln W)′, S2 = −det[P,P‴,P″]/W, S3 = det[P‴,P′,P″]/W (rows j). Max relative deviation from the eigenvalue functions (64 states × u ∈ {0.37, 1.9, −1.3}): S1 3.5e-51, S2 3.7e-51, S3 3.7e-51.
W = c_W e^{(ζ1+ζ2+ζ3)u} Π(u)³ with c_W = ∏_{i<j}(ζ_j−ζ_i) = −35/288 (state-independent, dev 2.7e-54).
Direct side: ψ_i = Π⁻²Q_i span the kernel of the monic operator ∂³ − S1^ψ ∂² + ... with (−S1, 2S1′+S2, −(S1″+S2′+S3)): deviations 2.6e-48, 6.3e-48, 1.3e-46 (worse only because of 3 numerical derivatives).
The XXX rows 10–11 (ν-factor × determinants of shifts) therefore go over to ratios of Wronskian-type determinants of derivatives, as stated.

## 3. Q3 — residue brackets, orthogonality, det N, x0-free ratio

**(a) Orthogonality.** The spec's search among ⟨P^{(a)}Q^{(b)}⟩ with fixed coefficients cannot work as posed (the relation involves the state-dependent differences of eigenvalues, as δI₁ in XXX), so I derived the exact relation and tested it.
Lagrange identity: ψ·D_B[P] − P·D_B†[ψ] = C(P,ψ)′ with P = P^j (ket B), ψ = Π⁻²Q_i (bra A, same ζ). Since D_B†ψ = ΔD†ψ = −(ΔS2 ψ)′ − ΔS3 ψ (ΔS_k = S_k^B − S_k^A; S1 and all constants cancel), taking Res at each θ_s separately (residue of a derivative vanishes):

  Σ_m Δ_m N[(s,i,j),m] = 0  for s = 1,2 and all i,j,  Δ = (ΔH1, ΔH2, Δa1, Δa2, Δb1, Δb2)  (the six non-scalar Laurent data of ΔS2, ΔS3),
  N[(s,i,j),H_t] = −Res_{θ_s}[ Π⁻²Q_i · P^j′ / (u−θ_t) ],  N[·,a_t] = Res_{θ_s}[ Π⁻²Q_i P^j/(u−θ_t)² ],  N[·,b_t] = Res_{θ_s}[ Π⁻²Q_i P^j/(u−θ_t) ]  (t = 1,2).

(Residues are Taylor coefficients of Q_iP^{j(a)}·(other-site factors); order ≤ 3.) Results: identity residual max |ΣΔ_mN|/Σ|Δ_mN| = 6.3e-57 over 300 random ordered pairs A≠B × 18 rows.
**Rank test** (stacked 18×6 matrix, 9 (i,j) × 2 sites, columns normalised): A = B (64 states): σ6/σ1 ≥ 1.7e-7, rank 6. A≠B: always σ6/σ1 ≤ 4.3e-59, i.e. rank ≤ 5; 40 random same-sector pairs: rank exactly 5; 60 random different-sector pairs: rank 5 in 48, rank 4 in 12 (σ5/σ1 ~ 1e-60; a second relation, apparently of weight/asymptotic origin, not investigated). So rank 5 is the generic value (same-sector pairs and most others), never 6.
Negative result: the naive per-site columns {⟨PQ⟩_s, ⟨PQ′⟩_s, ⟨P′Q⟩_s}, ⟨f⟩_s = f″(θ_s)/2, with the raw direct Q (no Π⁻² weight), s = 1,2: for 16 of 20 pairs A≠B the 18×6 matrix has full rank (σ6/σ1 median 1.7e-3); the other 4 pairs are rank-deficient (not analysed). So these columns do NOT give the orthogonality; the Π⁻² weights and the six (H,a,b)×(site) columns are needed.
**(b) det N**: 6×6, rows (s,(i,j)) with (i,j) ∈ {(1,3),(1,2),(2,3)}, s = 1,2; columns (H1,H2,a1,a2,b1,b2) as above. |det N|/∏(row norms): A≠B max over 100 pairs 1.6e-65; A = B min over 64 states 1.75e-23. So det N vanishes exactly for A≠B (same ζ) and not for A = B.
**(c) x0-free overlap ratio**, bra vacua Ω_{ζω}⊗Ω_{ζω}, ket Ψ_B (right eigenvector, any normalisation), D_B(ω) = det N[ket dual P^j, vacuum direct Q_i^{vac,ω}], B0 = e1⊗e1 (top state, index 63, weight (4,2,0)):
R_B(ω)/R_B(ω′) = [D_B(ω)/D_B(ω′)]·[D_{B0}(ω′)/D_{B0}(ω)] holds on **all 64 states, max relative deviation 6.4e-55** (60-digit arithmetic; worst states (2,1,3), (2,2,2), (1,2,3)); e.g. B = (0,2,4) state: 2.45930790198217 both sides.
Row-triple scan (same triple on both sites): of the 84 triples of (i,j) pairs, 72 pass at < 1e-30 (the XXX triple among them), 6 are degenerate (all rows with the same i or the same j: some det = 0) and 6 fail (deviation ≈ 1): the triples forming a permutation (i_k,j_k all distinct in both indices).
Controls (should fail, do): replacing the a2 column by Res[Π⁻²Q P/(u−θ2)³]: deviation 2.0; swapped orientation (vacuum dual in the ket slot, state direct in the bra slot): deviation 1.6.
The overlap data (normalisation of Ψ_B and the monic Q's) cancels in the ratio. (d) not done.

## 4. Q4 — untwisting regime, ε = 1e-6 and 2e-6 (100 digits)

φ = (1/3, −2/7, 5/11), ω̂ = (3/4, −1/5, 1/2), ω̂′ = (−12/5, 1/8, −1/3) (third components not fixed by the plan; arbitrary). r_k(ζ,ω) is ε-independent (checked: same value at both ε; the vacua are exact rationals).
Branch selection (smallest raising residual per HW sector, number m = 1,1,1,2,1 for (4,2,0),(4,1,1),(3,3,0),(3,2,1),(2,2,2)); residual/gap at ε = 1e-6: 0 | 2.3e-7 vs 1.41 | 1.9e-7 vs 1.41 | 2.2e-7, 2.2e-7 vs 2.05 | 3.5e-7 vs 1.37 — clean gaps of 6–7 orders.
Q's of the six branches: degrees as at finite twist, kernel dimension 1. Ratio test on the 6 branches: max deviation 5.4e-67 (ε = 1e-6), 2.2e-69 .. 1.8e-69 (ε = 2e-6) at 100 digits.
Cancellation (det N with monic Q's, ε-exponents p from ε = 1e-6, 2e-6): |det N(ω)| ~ ε^12 for (4,2,0), (4,1,1), (3,3,0); ε^11 for both (3,2,1) branches; ε^10 for (2,2,2); same for ω′; max |N entry| ~ ε^{−2..−3};
det/∏(column norms) ~ ε^{22,24,28,23,23,20}; σ6/σ1 of the column-normalised N ~ ε^{7,7,8,7,7,6}; |R(ω)| = |⟨Ω⊗Ω|Ψ_B⟩| (Ψ_B max-component normalised) ~ ε^{0,1,1,1,1,2}. So D_B is not O(1): the ratio D_B(ω)/D_B(ω′) is a ratio of two ε^p quantities with equal p; the 100-digit precision suffices (det/∏cols ≈ 1e-185 for (3,3,0) is a column-scale effect, the 6x6 conditioning σ6/σ1 ≥ 1e-55). Record for the later two-vacua stage.

## 5. Surprises and open points

1. The plain formal adjoint of the dual oper is NOT the direct oper: direct Q_i = Π² ψ_i. With the ν-factors of the XXX Baxter equations this is the Gaudin counterpart of ν_1(x+h)ν_1(x+2h)-type weights. The gauge power 2 was fixed by degrees for [2,1,0] only; for [3,1,0] expect λ1−λ3 = 3 or λ1 = 3 (untested).
2. S3 has no triple pole at all (the plan allowed ≤ 3); the double poles are state dependent for S3 (9/5 vs 27/7 for the vacua) and scalar (3) for S2.
3. The spec's naive columns {⟨PQ⟩,⟨PQ′⟩,⟨P′Q⟩} per site fail the generic rank test (see 3a); the replacement columns are residues with the extra weights Π⁻² (from the gauge of the direct oper) and 1/(u−θ_t)^m (the Laurent structure of ΔS2, ΔS3). A moment-type version (poles of order λ1+1..λ1+3 with insertions) was not tried since the derived columns already close.
4. For A≠B in different weight sectors the 18×6 matrix has rank 4 in about 20% of the pairs (a second, unexplained null vector); det N for the XXX row triple vanishes regardless (max 1.6e-65 over 100 pairs).
5. For the mixed-twist overlap (vacuum twist ω ≠ ζ) the same functionals (not the 9-parameter extended Lagrange identity with ΔS1, Δe2, Δe3) give the state-independent ratio at 1e-55; why no extra columns are needed was not derived.
6. 72/84 row triples work; the 6 permutation triples fail. The XXX triple (1,3),(1,2),(2,3) is fine.
7. HW index 7 (not 0) in the `gl3_reps` ordering; any covector/ket indexing in the notebooks' `ps` order must be mapped.
8. Q4 pitfall: with ε-suppressed leading coefficients a threshold-based degree reading gives deg 3 instead of 4 for the (4,2,0) dual Q; the solver now reports the smallest d admitting a kernel.
