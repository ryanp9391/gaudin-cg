# Rotations of the SoV vacuum at L=2: term-by-term reading of the site-1 rotation, and a product of two distinct rotations with a contour-integral determinant overlap

**Date:** 2026-10-09. **Notebook:** `L2-rotation-sov.wb` §8–§10 (cells 59–103). **Companion files:** `site1-Qform/5-bracket-decomposition.md`, `site1-Qform/6-two-rotations-2603-trick.md`, `Rotations-SoV-Q.{md,tex,pdf}` (the 2026-10-02 results this builds on). **Reference:** Bargheer–Bercini–Lefundes–Ryan, arXiv:2603.26617, Appendix G (`References/2603.26617v1.pdf`).

Everything below is gl(3), site representation Λ = (2,1,0) (adjoint), two sites, new conventions (Lax u + hE, twist on the right, companion-frame B/C).

---

## 0. Summary

1. **Site-1 rotation ⊗ SoV vacuum (Direction 5): no single bracket determinant, but a compact exact form.** In the 2-plane variables (α, β, λ₁, λ₂) of §2, the overlap of *any* site-1 covector Σ_ν c_ν x_ν ⊗ e₁ᵀ is a short polynomial (table in §3.1). Ladder-1 bracket determinants of the natural type are site-1 objects only if F₂ ∝ F₁ on {θ₁−2h, θ₁−h}; that family is a 4-dimensional cone and the rotation groups lie outside it (exact Gröbner test). Sums of two such determinants cover everything and are therefore empty statements.

2. **Two different global twists, Ω_w ⊗ Ω_{w′} (§4): negative.** Not a det₆ with ladder-dependent vacuum data (deviations 5–25), and its SoV bra wavefunction does not factorize over the two sites (measure-free test, deviation 1.6).

3. **The 2603.26617 trick works for the adjoint (§5), with one twist of roles.** With

   y* = w₂ (θ₁₂ − h)/(θ₁₂ − 2h),

   the covectors v₂ = Ω_{(w₁, y*)} and v₁ = Ω_{(w₁, y*)} ρ((w₂ − G(z)ᵀ)⁻¹) are two distinct exact rotations of the SoV vacuum, and

   ⟨v₁ ⊗ v₂|Ψ_B⟩ / ⟨x₀|Ψ_B⟩ = c · det₆N[Q_B, Q_vac,(w₁,w₂)] / Ψ_Q^B(x₀) · det₂^{(1)}[Q_B; y*] / minor₁(B),

   with det₂^{(1)} the ladder-1 su(2)-like contour-bracket determinant with insertion y*^{−x/h}. The constant c is state-independent to 10⁻⁵⁰ on all 64 states. The state is the global vacuum of twist (w₁, y*) with a site-1 resolvent rotation at w₂, i.e. the ⟨W_{ω,κ} ⊗ W_ω| structure of the paper with ω ↔ (w₁, y*), κ ↔ 1/w₂. Unlike the rectangular case there is no continuous κ-family: the product condition fixes the insertion point to y*.

4. **Beyond the vacuum (§5.8, added later the same day): any rotation on site 2, one fixed partner on site 1; no κ-family; no iteration.** The ket-twist torus ρ(D)⊗ρ(D) acts on Ψ_B by the character Π d_i^{μ_i}, so the family comes with its torus orbit for free. But w₁ is itself a torus direction and y is a torus *invariant* of the rotation (exact Jacobian ranks 3/4/3). Hence the family is V = {(g ρ((w₂(y(g))−Gᵀ)⁻¹), g) : g any rotation of x₀}, 3-dimensional: every site-2 rotation g (vacuum or not) has exactly one site-1 partner with overlap ∝ χ_μ(D) det₆N[Q_B, Q_vac,(w₁,w₂(y(g)))] det₂^{(1)}[Q_B; y(g)] (10⁻⁵¹ on 64 states for a random g). The paper's continuous κ-family does not exist for the adjoint even with the torus, and the §10 step does not iterate (on v₁⊗v₂ the site-1 T's have full site-2 column space; product locus = the seed only).

---

## 1. Setup and conventions (as in the notebook)

**Numeric point P1.** h = 1/2, θ = (2/7, −3/5) (so t := θ₁₂/h = 62/35), ket twist z = (3/2, 5/4, 2/3), bra twist w = (5/2, 7/3, 3/5), second bra twist w′ = (7/4, −2/5, 9/7). Precision 60 digits unless stated. Caution: z₃ = 1/z₁ at P1 (the §5 structure was re-checked at z = (3/2,5/4,3/7), (2,5/4,2/3), (3/2,7/6,2/3)).

**Chain.** ℒ_ij(u) = u δ_ij + h E_ij on the 8-dimensional adjoint (GT basis, generators `JJ`). Monodromy T(v) = ℒ^{(1)}(v−θ₁) ℒ^{(2)}(v−θ₂), twist on the right T^M = T·M, quantum minors with column k at v−(k−1)h, transfer matrices t_a(u) = Σ_{|I|=a} qmin_{I,I}(u), a = 1,2,3. Fused (CBR) transfer matrices at v: T_{10} = t₁(v), T_{11} = t₂(v), T_{20}(v) = t₁(v)t₁(v+h) − t₂(v+h), T_{21}(v) = t₂(v)t₁(v+h) − t₃(v+h). The site-1 point is v₀ = θ₁ − 2h = θ₁ − λ₁h.

**Frames.** Diagonal frame: twist diag(z); eigenstates Ψ_B (64, simple spectrum of t₁ + π t₂ at u = 3/7). Companion frame: ket Pinv2 = P_phys⁻¹ ⊗ P_phys⁻¹ with P_phys = P₋ P_d P₊ (cell 7), so that the SoV vacuum is x₀ = e₁ᵀ ⊗ e₁ᵀ. The companion matrix of a twist y is G(y) = [[χ₁, −χ₂, χ₃],[1,0,0],[0,1,0]] with χ the elementary symmetric functions.

**One-site vacua.** Companion frame: Ω_{(y₁,y₂)} = e₁ᵀ exp(E₁₂/y₁) exp(−E₁₃/(y₁y₂)) exp(((y₁+y₂)/(y₁y₂)) E₂₃) (`Ωc[y]`, cell 25); it depends on y₁, y₂ only. Diagonal frame: Ω₁ = e₁ᵀ exp(r₁E₁₂) exp(r₂E₁₃) exp(r₃E₂₃) with the rational r_k(z, w) of cell 6 (`oneSiteVac[z,w]`, cell 98). Ω_w ⊗ Ω_w is the exact left eigenvector of t_a(G′_w). Rotation parameters: under the equivariant map X (below), X(Ω_{(p,q)}) = a bᵀ with a = (1, −(p+q), pq), b ∝ (1, 1/p, 1/p²).

**Resolvent rotation.** ρ(M) := exp(Σ_{ik} (log M)_{ik} E_{ik}) for M ∈ GL(3) (`ρG`, cell 69). At L=1, every generic rotation of x₀ is ⟨x₀|M(w)M(y₁,∞) = Ω_w ρ((y₁ − G(z)ᵀ)⁻¹) (Rotations-SoV-Q items 6–8).

**Equivariant map X.** The unique (up to scale) linear map from 8-component covectors (GT basis, right action φ ↦ φ·E_ab) to 3×3 matrices with X(φ·(E_ab − δ_ab)) = [X(φ), e_ab] (cell 63; variant with a 1-dim null space). X(e₁ᵀ) = 3e₃₁. A covector is a rotation of x₀ iff X is rank 1 (traceless automatically). Diagonal-frame covectors are mapped by X(v·P_phys).

**Q-functions.** Dual Q's Q^j(x) = z_j^{x/h} q^j(x) (v(x) := (Q¹,Q²,Q³)(x)), direct Q_i(x) = z_i^{−x/h} q_i(x); vacuum direct Q's Q_i^{vac,w}(x) = w_i^{−x/h} q_i^{vac}(x) from the Baxter equation with the vacuum eigenvalues (degrees 0,2,4). Ψ_Q^B(x₀) = Π_s [Q²(θ_s−2h)Q³(θ_s−3h) − Q³(θ_s−2h)Q²(θ_s−3h)] Q³(θ_s−h).

**Brackets and det₆** (cell 8, Paul's preferred brackets). Poles x_p = θ_s − nh, s = 1,2, n = 0,1,2, weights wgt(p) = 1/Π_{q≠p}(x_p − x_q) over all six poles, insertions ins₁ = 1, ins₂ = e^{2πix/h}. Rows (k, (i,j)) with k ∈ {1,2}, (i,j) ∈ {(1,3),(1,2),(2,3)}; columns Σ_p wgt(p) ins_k(x_p) · {−Q^j(x)F_i(x), −xQ^j(x)F_i(x), Q^j(x)F_i(x+h), xQ^j(x)F_i(x+h), Q^j(x−h)F_i(x), xQ^j(x−h)F_i(x)}. With F = Q^{vac,w}: ⟨Ω_w⊗Ω_w|Ψ_B⟩/⟨x₀|Ψ_B⟩ ∝ det₆/Ψ_Q (10⁻⁵⁰). Since e^{2πix/h} is constant on each ladder, det₆ equals (e^{2πiθ₂/h} − e^{2πiθ₁/h})³ times the per-ladder det₆L (rows (s,(i,j)) with the sum restricted to ladder s).

**Ladder-1 brackets** (cell 68). Poles θ₁ − nh, n = 0,1,2, one-site weights wgt₁(n) = 1/Π_{m≠n}(m−n)h, ⟨f⟩₁ = Σ_n wgt₁(n) f(θ₁−nh). The ladder-1 su(2)-like determinant with insertion E(x) = y^{−(x−θ₁)/h}:

det₂^{(1)}[Q_B; y] = det_{j∈{3,2}} [ −⟨Q^j E⟩₁ , ⟨Q^{j[−]} E⟩₁ ],  Q^{j[−]}(x) = Q^j(x−h),

minor₁(B) = Q²(θ₁−2h)Q³(θ₁−3h) − Q³(θ₁−2h)Q²(θ₁−3h). This is the gl(3) analogue of (G7) of the paper. (The L=1-type 3×3 `detN1` with rows (1,3),(1,2),(2,3) and columns −⟨Q^jF_i⟩, ⟨Q^jF_i^{[+]}⟩, ⟨Q^{j[−]}F_i⟩ is also defined there.)

**SoV basis and 12-state sum** (§3–§4 of the notebook, 2026-10-02). The untwisted B(u) = T₁₁(u) qmin_{12,12}(u+h) + T₂₁(u) qmin_{12,13}(u+h) has a simple spectrum; its 64 left eigenvectors ⟨x| are labelled by (site-1 triple, site-2 triple) of roots θ_s + m h with triples among (−3,−2,−1), (−3,−1,−1), (−2,−1,−1), (−3,−2,0), (−3,−1,0), (−3,0,0), (−2,−1,0), (−2,0,0) (this `trip` order = `gtpats` order {2|21},{2|11},{1|11},{2|20},{2|10},{2|00},{1|10},{1|00}). Wavefunction: ⟨x|Ψ_B⟩/⟨x₀|Ψ_B⟩ = N_x Π_s det₂[v(x_{s,21}), v(x_{s,1})]_{rows 1,2} v₁(x_{s,22}) / (same at (−3,−2,−1)). A site-1 covector x_ν ⊗ e₁ᵀ expands on 12 SoV states (site 2 at vacuum or at {2|20}, the latter with factor λ₂ = v_k(θ₂)/v_k(θ₂−h)); the rational 8×12 map K is tabulated in `Rotations-SoV-Q.md` (the 31s and 35s are t = 62/35).

---

## 2. The 2-plane variables

On ladder 1 the dual Q's satisfy v(θ₁) ∥ v(θ₁−h) and v(θ₁−3h), v(θ₁−2h), v(θ₁−h) span a 2-plane (checked: residual 0 at 60 digits). Define per state B

- λ_s = v_k(θ_s)/v_k(θ_s−h), s = 1,2 (k-independent);
- v(θ₁−h) = α v(θ₁−3h) + β v(θ₁−2h).

Then the eight site-1 wavefunctions, in `trip` order and normalized at the vacuum triple, are exactly

ψ⁽¹⁾ = {1, β, −α, λ₁, λ₁β, λ₁²β, −λ₁α, −λ₁²α}

(max relative deviation 10⁻⁴⁷ over 64 states × 8 triples). By the dual Baxter equation at x = θ₁, θ₁−h, the numbers α, β, λ₁, λ₂ are rational functions of the transfer-matrix eigenvalues at θ₁, θ₁−h, θ₂: e.g. λ₁ = t₂(θ₁)/(t₁(θ₁) ν₃(θ₁−h)).

---

## 3. Part A — the site-1 rotation read term by term (Direction 5; cells 59–80)

### 3.1 The exact polynomial

For a general site-1 covector Σ_ν c_ν x_ν ⊗ e₁ᵀ (c free), the overlap over ⟨x₀|Ψ_B⟩ is, exactly (10⁻⁴⁷ on all 64 states, for a generic rotation c and for a random c):

| site-2 factor | λ₁ power | coefficient |
|---|---|---|
| λ₂⁰ | 1 | c₁ + (8β/27) c₂ + (4α/9) c₃ |
| λ₂⁰ | λ₁ | −(35c₂ + 62c₄)/124 − (4β/27) c₅ − α(35c₆ + 24c₇)/108 |
| λ₂⁰ | λ₁² | (β/31) c₆ + (3α/62) c₈ |
| λ₂¹ | 1 | −(3395/3348) c₂ |
| λ₂¹ | λ₁ | (35/54) c₅ − (35β/124) c₆ + (35α/558) c₈ |

Equivalently Σ_x φ⁽¹⁾_x ψ_x + λ₂ Σ_x φ⁽²⁾_x ψ_x with φ⁽ⁱ⁾ = c·K⁽ⁱ⁾, K^vac of rank 8 and K^20 of rank 4. (Compare L=1: c₁ + βc₂ − αc₃ + λ₁(c₄ + βc₅ − αc₇) + λ₁²(βc₆ − αc₈).) For a rotation parametrized by X = a bᵀ (b·a = 0), the bra coefficients on the 12 states are explicit bilinear forms in (a,b) (cell 80), e.g. φ_{(−3,−2,−1),vac} = a₃b₁/3, φ_{(−3,−1,−1),vac} = 16a₃b₂/81, φ_{(−2,−1,−1),vac} = 8(a₁b₁ + a₂b₂)/27, φ_{(−3,−2,0),vac} = (62a₂b₁ + 27a₃b₂)/186, φ_{(−3,−1,0),vac} = −8(a₁b₁ − a₂b₂)/81, ….

### 3.2 Neither group is a rotation

For five generic g, X(c·K^vac) and X(c·K^20) have three nonzero singular values (e.g. 2.93, 0.19, 0.0070 relative), i.e. generic traceless 3×3, not rank 1. So no site-2 group is an L=1-type determinant with transformed bra data.

### 3.3 The L=1 det₃ × det₂ is not a polynomial identity

Fed the L=2 dual Q's on ladder 1, the L=1 formula det₃N[Q_B, Q_vac,w] · det₂^{(1)}[Q_B; y₁]/(Ψ_Q^{(1)} minor₁) is **not** in the 8-dim span of the site-1 wavefunctions: relative residual 0.112 (L=1 vacuum Q's of w at site 1) and 0.113 (L=2 vacuum Q's), and the fitted coefficient vectors are not rotations (0.19, 0.18). At L=1 every function of the 8 states is in that span, which is why the formula holds there. Likewise the ladder-1 resolvent determinant det N_∞(y) is not the λ₂⁰ group of the site-1 resolvent rotation (0.93, 0.77, 0.63 for three y).

### 3.4 Locality classifies the usable determinants (symbolic)

Write v(−3) = u, v(−2) = w, v(−1) = αu + βw, v(0) = λ₁(αu + βw). Every ladder-1 bracket is A u_j + B w_j with A = μ(−3) + α m, B = μ(−2) + β m, m = μ(−1) + λ₁ μ(0) (μ the measure × bra data × insertion). A bracket determinant over Ψ_Q^{(1)} = (w₂u₃ − w₃u₂)(αu₃ + βw₃) is a function of (α,β,λ₁) alone — a *local* object, i.e. a covector supported on SoV states with site 2 at the vacuum — iff its numerator is divisible by Ψ_Q^{(1)} as a polynomial in u, w. For the L=1-type 3×3 (rows (1,3),(1,2),(2,3)) the monomials u₃²u₂ and w₃²w₂ drop out because rows 1,2 share F₁, and the remaining condition is (exact `Solve`, cell 71): **F₂ ∝ F₁ on {θ₁−2h, θ₁−h}** (other components have F₁ vanishing on points). With the gauge F₂ → F₂ + cF₁ (row 3 → row 3 + c·row 1), the family is F₁ free on the four points m = −2..1 and F₂ supported on {θ₁, θ₁+h}:

φ(F) = (16f₋₂²p, −32f₋₂f₋₁p, 48f₋₁²p, 8f₋₂f₋₁p − 8f₋₂²q, 16f₋₂f₋₁q, 8f₋₂f₁p − 8f₋₂f₀q, −8f₋₂f₁p − 8(4f₋₁² − f₋₂f₀)q, …)

(ψ-basis coefficients; f_m = F₁(θ₁+mh), p = F₂(θ₁), q = F₂(θ₁+h)). Its image is a **4-dimensional cone** (Jacobian rank 4) with an ideal of 10 quadrics, e.g. φ₂φ₆ = φ₁φ₈, 4φ₃φ₆ = 3φ₂φ₈, φ₆² + φ₆φ₇ = φ₅φ₈. Three independent bra functions do not enlarge it (32 components, max in-span rank 4). The "Cauchy–Binet product" det₂[⟨Q^jF_i⟩] × ⟨Q¹G⟩ has a 6-dim image and already fails for L=1 rotations (2.4).

### 3.5 Exact membership and sums

Gröbner membership of the targets: control member 0, random vector 0.41, group 1 of rotations 0.16, 0.23, 3.9, group 2: 2.6, 0.18, 0.40. Numeric fits agree (best 0.18 and 0.014; an exact member fits to 10⁻¹⁷). **Sums of two local 3×3 determinants fit every target, including a random 8-vector, to 10⁻¹⁵**: a "sum of four determinants" form exists for every site-1 covector and says nothing.

**Conclusion of Part A.** The site-1 rotation is a bracket of a short, explicit, non-factorizing three-point kernel on the ladder-1 lattice (times {1, λ₂}); not a single determinant of the natural type.

---

## 4. Part B — two different global twists (cells 81–85)

Target Ω_w ⊗ Ω_{w′}. The w′-vacuum Q's were built as for w (eigen-bra of t(G′_{w′}) exact; degrees 0,2,4; sanity Ω_{w′}⊗Ω_{w′} ∝ det₆[F_vac,w′] to 10⁻⁴⁸).

- det₆ with ladder-dependent bra data (w-vacuum Q's on the ladder-1 poles, w′-vacuum Q's on the ladder-2 poles): max relative deviation of the ratio 24.9 (L=2 vacuum Q's), 5.8 (one-site vacuum Q's), 6.1 and 6.1 (mixed); the per-ladder det₆L gives the same numbers (it is the same determinant up to a constant).
- Measure-free factorization test: with M_Φ(x) = c_x N_x the SoV coefficients (overlap/⟨x₀|Ψ⟩ = Σ_x M_Φ(x) ψ_B(x)), a product bra wavefunction Φ¹(x₁)Φ²(x₂) would force M_{ww′}M_{w′w} = M_{ww}M_{w′w′} pointwise whatever the SoV measure. Measured: 1.6; the ratio M_{ww′}/M_{ww} varies by a factor 46 over the site-1 triple. All four 8×8 coefficient matrices have full rank 8 (including M_{ww}, so the measure itself couples the ladders).

So product bras do not have product SoV wavefunctions in the B-eigenbasis; the site-1 rotation ⊗ e₁ᵀ is the w₂ → ∞ limit of this.

---

## 5. Part C — the 2603.26617 trick for the adjoint (cells 86–103)

### 5.1 What the paper does

⟨W^ℓ_{ω,κ} ⊗ W^{L−ℓ}_ω|Ψ⟩ = K_κ(τ) ⟨W^L_ω|Ψ⟩, with ⟨W_ω| the global twisted vacuum (a contour-bracket determinant), and 1 − κ t₊(θ₁) + κ² t_{[1,1]}(θ₁) a combination of transfer matrices at the site-1 point which (i) localizes on site 1 ((G2), proven for the rectangular 6-rep) and (ii) is diagonal on |Ψ⟩ with eigenvalue a ladder-1 contour-bracket determinant (G7). The deformed vacuum is ⟨W_{ω,κ}| with C_ij = (1−κz_j)/(1−κz_i) R_ij (G14); the κ-family exists because the localizing span meets the (5-dim) rotation cone of the 6-rep in a conic.

### 5.2 Localization on the twisted vacuum, and the extra direction ω₂

Diagonal frame, Ω = Ω₁ ⊗ Ω₁ (w-vacuum), fused T_λ at v₀ = θ₁−2h:

| λ | localizes on ⟨Ω| ? | prefactor f_λ | rank of ⟨Ω|T_λ as 8×8 |
|---|---|---|---|
| (1,0) | no (residual 0.615) | — | 2 |
| (1,1) | yes, exact | −1161/4900 | 1 |
| (2,0) | no (0.724) | — | 2 |
| (2,1) | yes, exact | −31347/343000 | 1 |

(Rank scan over shapes and points θ₁ + mh, m = −3..1: rank ≤ 2 only at m = −2, plus T₁₁ = t₂ at θ₁ (rank 2); t₃ is a scalar. Site-2 localizers t₂(θ₂+h) and t₂(θ₂+h)t₁(θ₂−h) have rank 1 on ⟨Ω|, but joint products with T₁₁, T₂₁ have rank 3, 3, 3, 8.)

Exact decomposition: **⟨Ω|T_λ(θ₁−2h) = v_λ ⊗ Ω₁ + v′_λ ⊗ ω₂** for all four shapes, with a single second site-2 direction ω₂ (v′ = 0 for (1,1), (2,1)). ω₂ is a rotation, **ω₂ = Ω_{(w₁, z₂)}** (a(ω₂) = (1, −(w₁+z₂), w₁z₂) at all data points), and **the whole line Ω₁ + κω₂ lies inside the rotation cone** (minors vanish identically in κ). The seven vectors {Ω₁, v₁₀, v₁₁, v₂₀, v₂₁, v′₁₀, v′₂₀} span 6 dimensions; v′₁₀ ∦ v′₂₀.

### 5.3 The product locus and the point Φ*

In the family Φ(c) = ⟨Ω|(c₀ + c₁T₁₀ + c₂T₁₁ + c₃T₂₀ + c₄T₂₁) = v(c)⊗Ω₁ + v′(c)⊗ω₂, product states are v′(c) ∧ v(c) = 0 (exact `Solve`, 7 components). All but one reduce to the localizing span (rotations only at Ω₁ and Ω₁T₂₁⁽¹⁾). The remaining component,

c* ∝ (−2187/1600, −243/80, −27/43, 1, 320/1161) at (h, θ) = P1,

has **no further rotation constraint**: Φ* = v* ⊗ (Ω₁ + κ*ω₂) with κ* = −1029/1606, both factors exact rotations, distinct from each other and from Ω₁, and ⟨Φ*|Ψ_B⟩ = f*(B)⟨Ω|Ψ_B⟩ with f* = c*·(1, τ₁₀, τ₁₁, τ₂₀, τ₂₁) (10⁻⁵⁶). c* is **independent of z** (checked at three z) and depends on (h, θ) (values at h = 1/3: (390224/419175, 6728/3105, −29/6, 1, −345/116); at θ = (1/5,−4/9): (−107653/4665600, −1183/8640, −13/77, 1, 11520/7007)).

In the enlarged linear family {1, T₁₀, T₁₁, T₂₀, T₂₁, t₂(θ₁), t₂(θ₂+h), t₂(θ₂+h)t₁(θ₂−h)} (site-2 row space 5-dim, 12 product components) there is no continuous family of products of rotations; the only new isolated point is ⟨Ω|(e₁ + e₈ t₂(θ₂+h)t₁(θ₂−h)) with e₁/e₈ = 337851/2195200.

### 5.4 Φ* is the ladder-1 contour determinant at y*

The ladder-1 su(2)-like determinant f_y(B) = det₂^{(1)}[Q_B; y]/minor₁(B) lies, for every y, in span{1, τ₁₀, τ₁₁, τ₂₀, τ₂₁}(θ₁−2h) (residual 10⁻⁶⁰). At P1, exactly,

f_y = 4y⁴ − 70 y³ τ₁₀ + y² ( (4900/43) τ₁₁ − (4900/27) τ₂₀ ) + (1372000/3483) y τ₂₁.

**These coefficients are exact, not fit artefacts** (added later the same day, cell 104): the rationalized values reproduce the 60-digit data to 10⁻⁴⁹, and from the three (h, θ) points (using c* ∝ c(y*)) each coefficient is a unique monomial in h and (θ₁₂ − kh):

det₂^{(1)}[Q_B; y]/minor₁ ∝ 2h²(θ₁₂−h)(θ₁₂−2h)(θ₁₂−3h) y⁴ + 2h(θ₁₂−h)(θ₁₂−3h) y³ τ₁₀ + (θ₁₂−h) y² τ₁₁ + (θ₁₂−3h) y² τ₂₀ + (2/(3h)) y τ₂₁

(all τ at θ₁−2h; in the notebook's normalization the coefficient of y²τ₂₀ is 1/(h³(θ₁₂−h)(θ₁₂−2h))). At P1, θ₁₂−h = 27/70, θ₁₂−2h = −8/70, θ₁₂−3h = −43/70, which produces the 70s, 27s and 43s. The localization prefactors are f₁₁ = (θ₁₂−h)(θ₁₂−3h) = −1161/4900 and f₂₁ = (θ₁₂−h)²(θ₁₂−3h) = −31347/343000 at P1, so with τ̂_λ = τ_λ/f_λ the τ₁₁, τ₂₁ terms are −27 τ̂₁₁ y² − 36 τ̂₂₁ y there.

Its coefficient curve c(y) is proportional to c* exactly at **y* = −63/8** (all 2×2 minors vanish to 10⁻⁴⁶), where Φ_y drops from rank 2 (second singular value 0.06–0.14 for generic y) to rank 1. Consistency: c₀/c₃ = −(108/4900)y*², c₁/c₃ = (27/70)y*, c₂/c₃ = −27/43, c₄/c₃ = −7560/(3483 y*) reproduce c*.

### 5.5 Decoding y* and the two factors

Five data points (self-contained pipeline `pipelineZH[z, w, h, θ]`, cells 98–100, which never touches the memoized Q's):

| h | θ | z | w | t = θ₁₂/h | y* | −y*/w₂ |
|---|---|---|---|---|---|---|
| 1/2 | (2/7,−3/5) | (3/2,5/4,2/3) | w | 62/35 | −63/8 | 27/8 |
| 1/2 | same | (3/2,5/4,3/7), (2,5/4,2/3), (3/2,7/6,2/3) | w | 62/35 | −63/8 | 27/8 |
| 1/2 | same | (3/2,5/4,2/3) | w′ | 62/35 | 27/20 | 27/8 |
| 1/3 | same | same | w | 93/35 | 406/69 | −58/23 |
| 1/2 | (1/5,−4/9) | same | w | 58/45 | −91/96 | 13/32 |
| 1/5 | (1/3,−1/7) | (3/2,7/6,2/3) | w′ | 50/21 | −29/20 | 29/8 |

All satisfy −y*/w₂ = −(t−1)/(t−2), i.e.

**y* = w₂ (θ₁₂ − h)/(θ₁₂ − 2h).**

Site-2 factor: a(Ω₁ + κ*ω₂) = (1, −(w₁ + y*), w₁y*), b ∝ (1, 1/w₁, 1/w₁²): **v₂ = Ω_{(w₁, y*)}** (verified exactly against `oneSiteVac[z, (w₁, y*, w₃)]` at all points).

Site-1 factor: a(v*) = (1, −345/62, 453/62) at P1 (irrational roots, changes with z), not a vacuum. Scan over Ω_{(p,q)} ρ(M) with p, q ∈ {w₁, w₂, y*} and resolvent-type M (cell 101): exact matches (residual 0)

**v₁ = Ω_{(w₁, y*)} ρ((w₂ − G(z)ᵀ)⁻¹) = Ω_{(w₂, y*)} ρ((w₁ − G(z)ᵀ)⁻¹)**,

while Ω_{(w₁,w₂)} ρ((y* − Gᵀ)⁻¹) and the (y* − G), (1 − Gᵀ/y*)⁻¹ variants fail (0.20–0.73). The paper-style rescaling of the Gauss parameters r_k by (1 − z_a/y*)/(1 − z_b/y*) factors is not it either (best 0.096).

### 5.6 End-to-end verification from the closed form (cell 102)

Build v₂ = Ω_{(w₁,y*)} and v₁ = v₂ ρ((w₂ − Gᵀ)⁻¹) in the companion frame from the formulas alone (both descriptions of v₁ agree to 10⁻⁴⁰; both are exact rotations). Then

⟨v₁⊗v₂|Ψ_B⟩/⟨x₀|Ψ_B⟩ ÷ [ det₆N[Q_B, Q_vac,(w₁,w₂)]/Ψ_Q · det₂^{(1)}[Q_B; y*]/minor₁ ] = c, state-independent to **10⁻⁵⁰ on all 64 states** (c ≈ (1.18 + 3.16 i)·10⁻¹⁰ in the notebook's normalizations).

The "paper-like" assignment, det₆ with the (w₁, y*)-vacuum Q's times det₂^{(1)} at w₂, fails (60.8), although Ω_{(w₁,y*)}⊗Ω_{(w₁,y*)} ∝ det₆[F_vac,(w₁,y*)] holds (10⁻⁵⁰).

### 5.6b Is y a free parameter? Yes, with the det₆ twist slaved to it (cell 106)

Read the other way round, y is free and the second twist parameter of the global determinant is fixed by it: w₂(y) = y (θ₁₂−2h)/(θ₁₂−h). For any y, with v₂ = Ω_{(w₁,y)} and v₁ = v₂ ρ((w₂(y) − G(z)ᵀ)⁻¹),

⟨v₁⊗v₂|Ψ_B⟩/⟨x₀|Ψ_B⟩ ∝ det₆N[Q_B, Q_vac,(w₁,w₂(y))]/Ψ_Q · det₂^{(1)}[Q_B; y]/minor₁.

Checked at y = 3 (w₂ = −8/9) and y = −1/2 (w₂ = 4/27) with freshly solved vacuum Q's: state-independence 10⁻⁵⁰, both factors exact rotations, distinct. So the family of such product states is two-dimensional, (w₁, y) (w₃ drops out), inside the eight-dimensional space of pairs of rotations; the ladder-1 insertion y and the det₆ vacuum's w₂ cannot be chosen independently. For fixed w the expansion f_y ∈ span{1, τ₁₀, τ₁₁, τ₂₀, τ₂₁} holds for all y, but ⟨Ω_w|𝒦_y is a product only at y = y*(w).

### 5.7 Interpretation and comparison with the paper

- Structure: v₁⊗v₂ = (Ω_{(w₁,y*)}⊗Ω_{(w₁,y*)}) · (ρ((w₂−Gᵀ)⁻¹) ⊗ 1): the global vacuum of twist (w₁, y*) with a site-1 resolvent rotation at w₂. This is ⟨W_{ω,κ}⊗W_ω| with ω ↔ (w₁, y*) and κ ↔ 1/w₂.
- The overlap is (global det with the (w₁, w₂) vacuum Q's) × (ladder-1 det with insertion y*^{−x/h}). The map (w₁,w₂) ↔ (w₁,y*) is the rescaling of the second twist parameter by (t−1)/(t−2); the global det "sees" w₂ and the ladder det "sees" y*, the opposite of the naive reading of (G18).
- Why no κ-family: the localizing span {1, T₁₁, T₂₁} meets the 5-dim rotation cone in C⁸ only in points (Ω₁, Ω₁T₂₁); the continuous ingredient has to come from the non-localizing T₁₀, T₂₀ through the line Ω₁ + κω₂, and the product condition then fixes the insertion to y*. In the rectangular case the localizing span alone meets the cone in a conic.
- The mechanism depends on two structural facts that should be tested for other reps and L: ⟨Ω|T_λ(θ₁−λ₁h) has site-2 column space of dimension 2 for all shapes, and that second direction is a rotation Ω_{(w₁,z₂)} with the connecting line inside the cone.

### 5.8 Beyond the vacuum: the ket-twist torus, any rotation on site 2, and why there is still no κ-family (cells 107–112, added 2026-10-09 later)

**Question (Paul).** Can the construction be generalized to states beyond the vacuum — more general rotations on each site?

**Torus = character (cell 108).** The Cartan torus of the ket twist, ρ(D)⊗ρ(D) with D = diag(d₁,d₂,d₃) (companion frame: ρ(f(G)) for f a function of G(z)), commutes with every twisted transfer matrix and multiplies an eigenstate Ψ_B of total weight μ by χ_μ(D) = Π_i d_i^{μ_i}; for the resolvent D = (y₁−Gᵀ)⁻¹ this character is the su(2)-like det₄ of `Rotations-SoV-Q` item 12. Checked on random covectors: ⟨φρ(D)⊗ψρ(D)|Ψ_B⟩ = χ_μ(D)⟨φ⊗ψ|Ψ_B⟩ exactly (60 digits). So every product state with a determinant overlap comes with its whole torus orbit.

**w₁ is a torus direction; y is a torus invariant (cell 108, exact).** Symbolic Jacobian of (w₁, y, D) ↦ Ω_{(w₁,y)}ρ(D) (diagonal frame, `oneSiteVac`·`ρdiag`) at three generic rational points: rank 3 at fixed y (= torus 2 + scale), rank 4 with y free (= the full 3-dim orbit, 1-dim fibres), rank 3 for the torus alone; the fibre tangent has zero y-component. Directly: Ω_{(w₁,y)}/Ω_{(w₁′,y)} is componentwise a torus character (residual 0). So the "(w₁, y)" family of §5.6b is really (torus direction, y), and y is an invariant of the rotation modulo the torus. (Caution: Ω_{(w₁,y)} degenerates at y = z_i; the first run sat at y = 5/4 = z₂ and gave spurious numbers.)

**The family at fresh w₁ (cell 109).** All §10 determinant checks had w₁ = 5/2. At (w₁, y) = (−7/3, 7/5) and (9/4, −2/3), with v₂ = Ω_{(w₁,y)}, v₁ = v₂ρ((w₂(y)−Gᵀ)⁻¹), w₂(y) = y(θ₁₂−2h)/(θ₁₂−h): ⟨v₁⊗v₂|Ψ⟩/⟨x₀|Ψ⟩ ∝ det₆N[Q_vac,(w₁,w₂(y))]/Ψ_Q × det₂^{(1)}(y)/minor₁ with state-independence 10⁻⁵¹ and 10⁻⁵⁰.

**A generic rotation g on site 2 (cell 111).** g = x₀ρ(lower unipotent · torus · upper unipotent) with random rationals (exact). Exact `Solve` of Ω_{(w₁,y)}ρ(D) = s·g at w₁ = −7/3 and at w₁ = 9/4: exactly one solution each, with the same y(g) = 1546301/1368484 (D determined up to scale). The partner g ρ((w₂(y(g))−Gᵀ)⁻¹) is a rotation distinct from g, and

⟨g ρ((w₂(y(g))−Gᵀ)⁻¹) ⊗ g|Ψ_B⟩/⟨x₀|Ψ_B⟩ ∝ χ_μ(D) · det₆N[Q_B, Q_vac,(w₁,w₂(y(g)))]/Ψ_Q · det₂^{(1)}[Q_B; y(g)]/minor₁,

state-independent to 10⁻⁵¹ on all 64 states (vacuum Q's solved fresh for the twist (w₁, w₂(y(g)))).

**Consequence.** The torus orbit of the §10 family is the 3-dimensional variety

V = {(g ρ((w₂(y(g))−Gᵀ)⁻¹), g) : g a rotation of x₀},

i.e. "any rotation on site 2, one fixed resolvent on site 1". Vacua Ω_w are not special: every rotation is Ω_{(w₁,y)}ρ(D). Since y(g) is unique, there is no continuous κ-family for the adjoint even after allowing the torus — the negative conclusion of §5.7 stands. The general pair (g₁, g₂) is 6-dimensional; V has codimension 3 (a free resolvent parameter and the second unipotent part of g₁ relative to g₂ are missing).

**The mechanism does not iterate (cell 110).** With Φ₀ = v₁⊗v₂ (the §10 state, exact rationals in the diagonal frame) as the seed, the 8×8 matrices of ⟨Φ₀|{1, T₁₀, T₁₁, T₂₀, T₂₁}(θ₁−2h) have ranks {1, 4, 3, 8, 8}; the site-2 column space is all of C⁸; the exact product locus in the 5-dim span (784 quadrics, `Solve` 9 s) is the seed alone. T₁₁ and T₂₁ do not localize on the deformed state, so a second site-1 resolvent cannot be generated at θ₁−2h.

**Remaining routes within single determinants:** the site-2 analogue (insertions at the site-2 points) on the vacuum and on V; products of T's at different points or the joint site-1 × site-2 span. Part A (§3) already shows that (g₁, x₀) with generic g₁ is not a single natural determinant, so a 6-dimensional single-determinant family is unlikely; sums of determinants are the fallback.

---

## 6. Open questions / next steps

1. Second numeric (h, θ) point for the determinant identity of §5.6 (the structure and y* are verified at five points; the det at P1 only — now also at fresh w₁ and at a generic g, §5.8, but still at P1's (h, θ)).
2. Closed form of the torus invariant y(g) of a rotation (a projective invariant of (a_i b_i)_i in the eigenbasis of G); then the formula of §5.8 is fully explicit for any g.
3. The site-2 analogue (insertions at t₂(θ₂+h), t₂(θ₂+h)t₁(θ₂−h) and their non-localizing partners) on ⟨Ω⊗Ω| and on the states of V; and the joint site-1 × site-2 span. These are the only remaining single-determinant routes to a free resolvent or a second unipotent parameter (§5.8: no κ-family, no iteration at θ₁−2h).
4. General reps and L > 2 — **[3,1,0] done 2026-10-09 evening, see `L2-310-two-rotations.md`**: same mechanism (one extra site-2 direction, two isolated product points), y* = w₂(θ₁₂−h)/(θ₁₂−λ₁h), but the site-1 factor is a rotated *dressed* vacuum e₁ᵀ(α+βE₁₂)RD rather than a rotation. General reps and L > 2 otherwise: the localizing set is larger (T_{(k,k)}, T_{(k,d)} at θ₁−λ₁h), the site-2 column space of the non-localizing T's may be larger than 2, and a continuous family may reappear. The recipe: localization on ⟨Ω⊗Ω|, column space of the non-localizing T's, exact product locus in the linear span, match to the ladder contour det's coefficient curve — and then the torus analysis of §5.8 to see what the family really is modulo the torus.
5. Closed form of c*(h, θ) (it is z-independent); it is the τ-expansion of det₂^{(1)} at y*, so it follows from the residue expansion once the T-eigenvalue dictionary at θ₁−2h is written down.
6. For Part A: a single-determinant form for the pure site-1 rotation would need brackets involving v(θ₁+h) (outside the 2-plane) or two-ladder insertions producing λ₂; nothing of the sort was found.

---

## 7. Log (2026-10-09)

Times are approximate; cell numbers refer to `L2-rotation-sov.wb` as saved at the end of the day (103 cells).

- **Start.** Read HANDOFF (last state 2026-10-02: site-1 vacuum not a per-ladder det₆, Directions 1–3 negative, Direction 4 planned). Paul's request: read the exact 12-state sum term by term, fix the site-2 factor, look for bracket structure; target the generic 3-parameter site-1 rotation. Spec written: `site1-Qform/5-bracket-decomposition.md`.
- **Kernel rebuild.** Fresh kernel; cells 2–21 rerun (~1 min). Setup checks reconfirmed (12-state sum exact; K exact; wavefunction exact).
- **§8, cells 59–61.** 2-plane variables; wavefunction reduction; symbolic polynomial with 8 free c_ν (table §3.1); checks 10⁻⁴⁷.
- **Cells 62–64.** Equivariant map X rebuilt (first attempt applied X to the x_ν-coefficient vector instead of the GT covector — fixed). Orbit test: groups not rotations.
- **Cells 65–66.** "det₂ bracket × bracket" product-form criterion: fails even at L=1 — discarded as the wrong notion.
- **Cells 67–69.** L=1 det₃×det₂ on ladder 1 of L=2 states: not in the local span (0.11). Rotation g = M(w)M(y₁,∞) built and confirmed rank 1 (frame conventions OK).
- **Cells 70–71.** Symbolic locality analysis (first version clobbered the global `uS`; see incident). Local family F₂ ∝ F₁ on {−2,−1}.
- **Cells 72–77.** Image map φ(F), resolvent member, numeric fits (first version discarded any start with a warning and reported Infinity — fixed; then a name clash `res2` with the span-fit residual — fixed, renamed `resTwo`), Gröbner ideal (10 generators), exact membership, sums of two dets.
- **Cells 78–80.** §8 conclusions; bra wavefunction bilinear in (a,b).
- **Paul:** "so it's not a det but still a sum of brackets?" — answered: yes, with a non-factorizing 3-point kernel. **Paul:** wants contour-integral form; fallback: two different rotations on the two sites.
- **§9, cells 81–85.** Second twist w′, det₆ with ladder-dependent data (fails), factorization test (fails), ranks.
- **Incident.** The §9 cell was first run with `uS` (the `tPoly` interpolation points) clobbered by cell 70's symbolic `uS = {u1,u2,u3}`; the transfer matrices were rebuilt on symbolic data, the kernel ran away with memory and **the laptop crashed**. After restart: VS Code/wolfbook reconnected, notebook on disk intact (82 cells), the rename edit had been lost and was re-applied (`uSym/wSym`), kernel state rebuilt (cells 2–21, 24–27, 60–68). Memory note written (`feedback_symbol_clash_setup_cells`). Lesson: grep the setup cells before introducing short global names; time cells that rebuild `tPoly`.
- **Paul:** "try everything; find ⟨v₁⊗v₂| with two distinct rotations such that it is a sum of SoV integrals; read 2603 — there we did it with tricks for rectangular reps." Read the paper (main text + App. D–H, K).
- **§10, cells 86–88.** Exp. 1: det₆ × ladder-1 det₂ covector is rank 2, factors rotations; f_y ∈ span{1,τ₁₀,τ₁₁,τ₂₀,τ₂₁}. Exp. 2: localization of T₁₁, T₂₁ on ⟨Ω|; linear-span rotations only Ω₁, Ω₁T₂₁.
- **Cells 89–92.** Rank-2 structure: Ω₁ ∈ column space, second direction ω₂ y-independent and a rotation; rank scan; joint localization negative; exact decomposition in {Ω₁, ω₂}; line inside the cone; product locus → component with free rotation (Φ*).
- **Cells 93–96.** Φ* verified; f_y expansion in y; y* = −63/8 exact; enlarged family (no continuous family); classification with both factors.
- **Cells 97–100.** Pipeline as a function of the vacuum / of (z, w) / of (z, w, h, θ): y* data at five points → y* = w₂(θ₁₂−h)/(θ₁₂−2h); site-2 factor = Ω_{(w₁,y*)}; ω₂ = Ω_{(w₁,z₂)}. (First `pipelineZ` failed silently because a Module local `ps` shadowed the global GT-pattern function — fixed.)
- **Cell 101.** Identification of v₁ by a candidate scan: Ω_{(w₁,y*)}ρ((w₂−Gᵀ)⁻¹) = Ω_{(w₂,y*)}ρ((w₁−Gᵀ)⁻¹), residual 0.
- **Cell 102.** End-to-end check from the closed form: first run forgot the ⟨x₀|Ψ⟩ normalization (deviation 99); fixed → 10⁻⁵⁰. Paper-like relabelling fails (60.8).
- **Cell 103.** §10 result cell.
- **Cell 106** (Paul: is y free?): y = 3 and y = −1/2 with w₂ = y(θ₁₂−2h)/(θ₁₂−h); identity 10⁻⁵⁰.
- **Cells 104–105** (after Paul asked whether the τ-coefficients are exact): 60-digit check and monomial identification in h, θ₁₂−kh; exact general form recorded. Files: `site1-Qform/6-two-rotations-2603-trick.md`, README row, HANDOFF, memory. Notebook saved; nothing committed.

- **Paul (later):** "Can the construction be generalized to states beyond the vacuum — more general rotations on each of the two sites?" §11, cells 107–112 (`Two-rotations-det-overlap.md` §5.8).
- **Cell 108.** Torus = character on random covectors (exact). Dominance Jacobian: first run at y = 5/4 = z₂ (degenerate — Ω_{(w₁,y)} has (y−z_i) denominators) gave 0/0 and misleading ranks; redone at three generic points: ranks {3, 4, 3}, fibre tangent with zero y-component; Ω_{(w₁,y)}/Ω_{(w₁′,y)} a torus character exactly. So w₁ is a torus direction and y a torus invariant — the "κ-family via the torus" expectation written in the first version of cell 107 was wrong and was rewritten.
- **Cell 109.** (w₁,y) family at fresh w₁ = −7/3 and 9/4: 10⁻⁵¹, 10⁻⁵⁰.
- **Cell 110.** Iteration with Φ₀ = v₁⊗v₂ as seed: ranks {1,4,3,8,8}, column space C⁸, product locus = seed only (784 quadrics, 9 s). (First attempt failed on syntax: the string "y*)" inside a comment closes the comment.)
- **Cell 111.** Random Gauss-decomposed rotation g: exact Solve gives one solution per w₁ with the same y(g); partner formula with χ_μ(D) state-independent to 10⁻⁵¹.
- **Cells 107, 112.** Intro rewritten to the corrected statement; result cell. Notebook saved; nothing committed.

**Checks summary (all in `checks` of the notebook):** 2-plane residual 0; wavefunction reduction 10⁻⁴⁷; Rsym = exact overlap 10⁻⁴⁷; X(cL1) rank 1 exactly; local branch residual identically 0; local f in span exactly; image dimension 4; Gröbner control 0; sum-of-two control 10⁻¹⁷; Ω_{w′}⊗Ω_{w′} ∝ det₆[F_vac,w′] 10⁻⁴⁸; M(Ω_w⊗Ω_{w′}) reproduces the overlap 10⁻⁴⁴; control det₆×det₄ is a product of rotations exactly; T₁₁, T₂₁ localization on ⟨Ω| exact; decomposition in {Ω₁,ω₂} exact; Φ* rank 1 exact; Φ* overlap identity 10⁻⁵⁶; y* minors 10⁻⁴⁶; y* formula at five points exact; v₁ identification 0; MAIN identity 10⁻⁵⁰ on 64 states; §5.8: torus character exact; Jacobian ranks {3,4,3} at three points exact; torus-character residual 0; family at fresh w₁ 10⁻⁵¹/10⁻⁵⁰; generic g partner identity 10⁻⁵¹; iteration product locus = seed only (exact).
