# gl(4) rectangular reps [S,S,0,0], L = 2: the 2603.26617 two-rotation family, its untwisting, and why it is clean where the gl(3) family is not (2026-10-10)

Notebooks: `L2-gl4-rect-untwist.wb` (S = 1, 34 code cells) and `L2-gl4-rect-untwist-S2.wb` (S = 2, 29 code cells),
built by a sonnet wolfbook-builder from `L2-gl4-rect-untwist-spec.md` (+ Addendum A), reviewed by the main session.
**No native outputs** (VS Code did not execute cells of the new files); every cell was run in K1 and the printed
results are stored in markdown cells under the code cells. Main-session sympy cross-checks: the vacuum block and its
cross-ratio, the ω-deformed rotation, the gl(3) commutator identity, and the one-site identities at θ₁−h.

Conventions: new conventions (Lax u δ + hE, twist diag(z) on the right), GT basis (`JJ`, rank-generic code), h = 1/2,
θ = (2/7, −3/5), pr = 200. Finite-twist point z = (3/2,5/4,2/3,4/7), w = (5/2,7/3,3/5,9/8). Evidence tags as in
`Untwist-two-rotation-notes.md`: [exact], [numeric], [conjecture].

## 0. Paul's question and the answer

The paper's ⟨W_{ω,κ} ⊗ W_ω| (two distinct rotated vacua; (G13)–(G18)) gives every CG coefficient of the 6-rep at
**leading order** in the untwisting limit z, κ, ω → 1. The gl(3) det₆ × det₂ family needs two ([2,1,0]) or three
([3,1,0]) orders because several CG slices of the undone covector sit at O(ε), O(ε²), where the first contamination
K = ⟨v₁⊗v₂|(1−P_HW)Ψ⟩ also sits. Why?

**Answer (established below).** Two ingredients, both absent in gl(3):

1. **Abelian nilradical.** For [S,S,0,0] every rotation of the HW covector is e₁ᵀexp(X), X = Σ_{i≤2<j} x_ij E_ij
   with commuting E's. The resolvent/torus undo of the family acts *entrywise*, X′_ij = x_ij(δ_j/δ_i − 1), with no
   commutator direction. In gl(3) the Heisenberg nilradical produces the direction E₁₃ whose undone coefficient is
   the commutator combination b′ − a′c′ = a′c′ (z₂−w₂)/(y−z₂) [exact, sympy + cell 14].
2. **Two independent localising deformations at one point, none pinned.** At the site-1 point θ₁ − S h the fused
   t₁, t₂ localise; at S = 1 they are exactly the complementary-pair torus, so the κ-insertion
   1 + (κ/h)t₁ + (κ²/2h²)t₂ is the resolvent torus ρ(diag 1/(w₂′−z)), and the ω-deformation (1 − ωρ(diag z)) is a
   second, independent, rotation-preserving insertion. Both untwist with z. The vacuum invariant (the cross-ratio
   of X) is moved by ω; the torus point by κ; nothing ties them. In gl(3) the product locus is one isolated point
   c* that ties the ladder insertion y to the resolvent point by y = w₂ (t−1)/(t−2): the resolvent torus is
   non-scalar at ε = 0 only if y → y₀ ≠ 1, and then (z₂−w₂)/(y−z₂) = O(ε) kills the E₁₃ direction.

A sharp corollary: in gl(4) itself, the κ-torus family on the plain companion vacuum (ω = 0) has det X′ ≡ 0 at
**any** twist — it never sees the singlet. And holding ω fixed while z → 1 reproduces the gl(3) disease exactly
(singlet slice O(ε²), K/S a constant ≠ 0). The cleanliness is a property of the *joint* limit (4) of the paper.

## 1. The vacuum, the ω-deformation, the family

**Companion-frame vacuum [exact].** Ω_zw = HW left eigenvector of the one-site t_a(u; G′), G′ = P_aux(z) G_c(w)
P_aux(z)⁻¹ (P_aux rows = left eigenvectors of the companion matrix of z). It is e₁ᵀexp(X) with

    x_ij = (−1)^i (z_{i′} − z_j)/(z₁ − z₂) · Π_{k=1,2} (w_k − z_i)/(w_k − z_j),   i′ = the other index of {1,2},

w₃, w₄ drop out, w₁, w₂ are pure torus directions, the torus invariant x₁₃x₂₄/(x₁₄x₂₃) equals the cross-ratio
(z₁−z₄)(z₂−z₃)/((z₁−z₃)(z₂−z₄)) of the z's (91/95 at P1), and the block is representation-independent (same at S = 2).

**The paper's W_ω [exact, S = 1].** In the diagonal frame ⟨W_ω| = Ω_zw (1 − ω ρ(diag z))/(1 − ωz₁z₂) is again a
rotation e₁ᵀexp(R), R_ij = x_ij (1 − ω z_a z_b)/(1 − ω z₁z₂) with (a,b) the weight of the e₁ᵀE_ij component — the
paper's (F5) — and CR(R) = CR(z)·(1−ωz₂z₃)(1−ωz₁z₄)/((1−ωz₂z₄)(1−ωz₁z₃)). The linear form 1 − ωρ(z) is **not** a
rotation at S = 2; there W_ω is defined as the group element e₁ᵀexp(R) with the same ratios.

**The family.** v₂ = W_ω, v₁ = W_ω ρ(diag 1/(w₂′ − z_k)) (resolvent torus; κ = 1/w₂′). Undo [exact on all HW kets]:
φ_u = v₁ρ(exp R)⁻¹ = δ^Λ e₁ᵀexp(R′), R′_ij = R_ij(δ_j/δ_i − 1), det R′ ∝ CR(R) − CR(z). The polynomial torus
ρ(diag(w₂′−z)) was run as a second family; structurally identical results, numbers quoted where they differ.

**Blind spot [exact, 3 points × 3 w₂′, both S]:** with v₂ = Ω_zw (no ω) the (0,0,S,S) slice of φ_u is exactly 0.

## 2. Q-number property (Task 1, finite twist, exact rationals)

- Localisation: ⟨Ω⊗Ω|t_a(θ₁ + mh) has rank 1 (= site-1 covector ⊗ Ω) only at **m = −S** for a = 1, 2 (S = 1: t₁, t₂;
  S = 2: all five shapes (1,0),(1,1),(2,0),(2,1),(2,2)). The paper's "θ₁" is our θ₁ − h (their ⟨Z| is our dual).
  Site 2 localises at t₁(θ₂) and t₂(θ₂+h).
- **One-site identities at θ₁−h, S = 1 [exact, symbolic z]:** on the weight e_c + e_d with complement {c′,d′},
  t₁^{(1)} = −h(z_{c′}+z_{d′}), t₂^{(1)} = 2h² z_{c′}z_{d′} = 2h²χ₄ ρ(diag z)⁻¹, t₃ = t₄ = 0. Hence
  (1−κz_{c′})(1−κz_{d′}) = 1 + (κ/h)t₁ + (κ²/2h²)t₂, and the paper's (G5) insertion is the **resolvent** torus
  in our conventions ((G14)'s ratio is the i↔j of the transposed rep).
- **Three-term relation [exact, 8 (ω, w₂′) pairs + 1 held out]:**
  ⟨v₁⊗v₂| = ⟨W_ω⊗W_ω| (c₀ + c₁t₁ + c₂t₂)(θ₁−h), c = [κ²/Π_k(1−κz_k)]·(1, κ/(h f₁), κ²/(2h²f₂)), with the site-2
  localisation prefactors f₁ = 27/70, f₂ = −54/1225 (ω-independent, same for W and Ω). Same over ⟨Ω⊗Ω|.
  The polynomial torus satisfies no three-term relation (needs t₁(θ₁−h) + four t₂'s).
- **The ω-deformation [exact]:** ⟨Ω⊗Ω|ρ₁⁻¹ = ⟨Ω⊗Ω|t₂(θ₁−h)/(2h²χ₄f₂); ⟨Ω⊗Ω|ρ₂ = −(1225/27)⟨Ω⊗Ω|t₂(θ₂+h);
  ⟨Ω⊗Ω|ρ₁ = a₀⟨Ω⊗Ω| + Σ_{m=−3}^{0} b_m⟨Ω⊗Ω|t₂(θ₁+mh), (a₀,b) = (62449/10780, −3/11, −8/33, 437/22, 2216/297);
  ρ₁ρ₂ = the character χ = Πz_k^{μ_k}. On all 36 eigenstates (accuracy 190 digits):
  ⟨W_ω⊗W_ω|Ψ⟩/⟨Ω⊗Ω|Ψ⟩ = [1 − ω(E₁+E₂) + ω²χ]/(1−ωz₁z₂)², E₁ = a₀ + Σ b_m τ₂(θ₁+mh), E₂ = −(1225/27)τ₂(θ₂+h).
  So the whole family's overlap is explicitly a function of transfer-matrix eigenvalues — a Q-number — at S = 1.
- **S = 2 [negative]:** the one-site t₂(θ₁−2h) is a degree-2 polynomial in z per weight (e.g. (3/2)z₃z₄ on
  (2,2,0,0), (3/4)(z₂+z₃)z₄ on (2,1,1,0), non-diagonal on the dim-2 weight (1,1,1,1)), not a torus; ρ(Δ⁻¹) on the
  20-rep is not in the span of the degree-≤2 monomials in t₁..t₄ at θ₁−2h (rank 6 → 7). The Q-number property of
  the family at S = 2 is **open** (candidates: pole-set products at θ₁−kh, k = 1..3, as in the gl(3) [4,1,0] chain).

## 3. Untwisting (Tasks 2–4)

Regime: z = 1+εφ, w = 1+εω_w, ω = 1+εo, w₂′ = 1+εk; φ = (1/3,−2/7,5/11,−4/9), φ′ = (−2/5,3/8,7/11,−1/6),
ω_w = (3/4,−1/5,2/9,−7/4), o ∈ {−12/5, 7/2, 1/8, 5/3}, k ∈ {−1/5, 9/7, 5/3}; ε ∈ {1,2,4}·10⁻⁶.

**Undone covector [exact, symbolic ε].** R′(0)_ij = x⁰_ij (o+φ_{i′}+φ_j)/(o+φ₁+φ₂) · (φ_j−φ_i)/(k−φ_j) (resolvent;
(φ_i−φ_j)/(k−φ_i) for the polynomial torus), det R′(0) ≠ 0 on the whole grid, φ_u(0)/δ^Λ = e₁ᵀexp(R′(0)).
S = 1: all 6 weights O(1), rank 1/1. S = 2: all 19 weights O(1), the dim-2 weight (1,1,1,1) rank 2/2. Identical for φ′
and a second ω_w. **Control, ω = 3/7 fixed:** S = 1: singlet slice O(ε²); S = 2: (0,0,2,2) O(ε⁴), (1,0,2,1) and the
other three (0,1,1,2)-type weights O(ε²), one of the two (1,1,1,1) functionals O(ε²).

**Signal vs contamination [numeric, pr = 200, undo exact to 155–198 digits].** Joint regime, max over the grid:

| S | branch | \|K/S\| at ε = 1, 2, 4·10⁻⁶ | ε-exponents S, K |
|---|---|---|---|
| 1 | [2,2,0,0] | K ≡ 0 | −2, — |
| 1 | [2,1,1,0] | 6.1·10⁻⁷, 1.22·10⁻⁶, 2.44·10⁻⁶ | −2, −1 |
| 1 | [1,1,1,1] | 1.89·10⁻⁷, 3.78·10⁻⁷, 7.55·10⁻⁷ | −2, −1 |
| 2 | [4,4,0,0] | K ≡ 0 | −4, — |
| 2 | [4,3,1,0] / [4,2,2,0] | 9.0·10⁻⁸ / 1.22·10⁻⁶ | −4, −3 |
| 2 | [3,3,1,1] / [3,2,2,1] / [2,2,2,2] | 6.8·10⁻⁷ / 4.0·10⁻⁷ / 3.8·10⁻⁷ | −4, −3 |

(φ′: 3.6·10⁻⁸, 2.7·10⁻⁷ at S = 1.) **Control, ω fixed:** K/S constant in ε: S = 1 singlet 0.0274 (resolvent) /
0.0162 (polynomial); S = 2: [2,2,2,2] 0.069, [3,2,2,1] 0.093, [3,3,1,1] 0.112 (resolvent) — the three irreps whose
slices dropped below O(1); the others unchanged.

**CG from the overlap alone [numeric].** C := R_C/(δ^Λ [e₁ᵀexp R′(0)]_slice): grid spread and error against the
chain slice both linear in ε (S = 1: 1.1·10⁻⁶ singlet, 5.9·10⁻⁷ [2,1,1,0] at ε = 10⁻⁶; S = 2: 0.8–2.3·10⁻⁶ for the
five non-top irreps), Richardson agreement 10⁻¹²–10⁻¹³. The dim-2 slice (1,1,1,1) of [3,3,1,1] is determined by an
exact 2×2 system on two grid points with **different o** (k only rescales the (1,1,1,1) row; same-o pairs are
singular): errors 1.3·10⁻⁶ … 4.4·10⁻⁶ at ε = 10⁻⁶, condition numbers 400–1000.

## 4. The gl(3) contrast (cell 14, symbolic)

n = [[1,a,b],[0,1,c],[0,0,1]], Δ = diag δ, n′ = Δ⁻¹nΔn⁻¹: n′₁₃ − n′₁₂n′₂₃ = (δ₃/δ₁−1)(b − ac z₂₁/z₃₁) for
δ_k = 1/(w₂−z_k); with the gl(3) vacuum n(w₁,y): b′ − a′c′ = a′c′ (z₂−w₂)/(y−z₂) [exact]. gl(4): n = 1 + X,
n′ = 1 + X′, X′_ij = x_ij(δ_j/δ_i−1), (X′)² = 0 [exact]. In the gl(3) regime y → y₀ = (t−1)/(t−2) ≠ 1 (forced by
the product condition w₂(y) = y(t−2)/(t−1) and the need for a non-scalar Δ), so (z₂−w₂)/(y−z₂) = O(ε).

## 5. Open

1. Q-number form of the family at S = 2 (which insertion realises ρ(Δ⁻¹) on the 20-rep; pole-set products).
2. The explicit Q-determinant ((G7)/(G18)) in our conventions, S = 1, as a check of the paper's formula.
3. gl(3): is there a second localising deformation (an ω-analogue) for [2,1,0] that moves the vacuum invariant
   y(g) independently of the resolvent point? The 2026-10-09 §11 analysis found the product locus rigid.
4. Rerun both notebooks in VS Code to store native outputs.
