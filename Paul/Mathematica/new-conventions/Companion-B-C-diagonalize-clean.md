# SoV for gl(3) XXX at L=1 in the new conventions — reference

These are notes for `Companion-B-C-diagonalize-clean.wb`, from work on 2026-10-01. The full write-up, with all formulas typeset, is `Companion-B-C-diagonalize-clean.tex`/`.pdf`; this file is the quick lookup. Everything here is at one site (L=1). The complete construction is for the adjoint λ=(2,1,0); results that also hold for other reps are in "Beyond the adjoint" below.

## Files

| File | Role |
|---|---|
| `Companion-B-C-diagonalize-clean.wb` | Final formulas and checks only. Self-contained (cell 3 runs `ClearAll["Global`*"]`), general h, one numeric point; ends in a 37-check pass/fail table, all passing. |
| `Companion-B-C-diagonalize.wb` | The working notebook, including the exploration (failed operator forms, ansatz scans, fixed bugs). |
| `Baxter-QQ-Wronskian.wb` | Baxter, QQ, Wronskian and orthogonality checks for any rep (`gbLam`), run for [2,1,0], [3,1,0], [3,2,0] and [2,1,1]. |
| `old-conventions/old/Gaudin/Companion-B-C-diagonalize-convention-change.nb` | The original notebook that was ported. |

**Numeric point of the clean notebook:** h=1/2, θ=2/7, z=(3/2,5/4,2/3) so χ₃=5/4, second twist w=(5/2,7/3,3/5), 60 digits. The working notebook also checked the same-twist results at h=1, θ=1/3 with z=(5/2,7/3,6/35) and z=(5/2,7/3,3/5).

## Conventions

- **Lax.** ℒ_ij(u) = u δ_ij + h E_ij, with T_ij(u) = ℒ_ij(u−θ). The old conventions were u − h Σ E_ji⊗e_ij with a diagonal twist.
- **Twist.** Companion G(z) = [[χ₁, −χ₂, χ₃], [1,0,0], [0,1,0]], acting on the right: T^G = ℒ·G. Its eigenvalues are z₁, z₂, z₃.
- **Quantum minors.** Rows antisymmetrized, column k at u−(k−1)h. Transfer matrices: t_a = Σ_{|I|=a} qmin_{I,I}[T^G].
- **Factors.** ν_n(u) = u − θ + hλ_n. They were u−θ−hλ_n in the old conventions.
- **Fused (CBR) transfer matrices.** T_λ(u) = det_{i,j≤λ₁} t_{λ'_j+i−j}(u+(i−1)h), with t₀ = 1 and t_a = 0 for a<0 or a>3.
- **Q-functions.** Direct Q_i = z_i^{−u/h} q_i; dual Q^j = z_j^{+u/h} q^j; all q's monic. Shift notation: Q^{[2k]}(u) = Q(u+kh).
- **Column vectors.** v(p) = (Q¹,Q²,Q³)(p) (dual) and u(p) = (Q₁,Q₂,Q₃)(p) (direct).
- **Weights.** μ is defined in the twist-diagonal frame and read off from t₁(θ) = h Σ z_i μ_i.

## Results

| # | Statement | Status |
|---|---|---|
| 1 | **Quantum determinant.** t₃ = χ₃ ν₁(u) ν₂(u−h) ν₃(u−2h)·𝟙 | exact, symbolic; [2,1,0] and [3,1,0] |
| 2 | **B and C.** Built from the *untwisted* T: B(u) = T₁₁(u) qmin₁₂,₁₂(u+h) + T₂₁(u) qmin₁₂,₁₃(u+h), and C(u) is the same with the minors at u | — |
| 3 | **B-eigencovectors.** ⟨x_μ| are built from ⟨x₀| = e₁ᵀ (pattern {2|21}) by T(θ−2h). B-roots x₁ = θ−(μ₁₁+1)h, x₂₁ = θ−μ₂₁h, x₂₂ = θ−μ₂₂h | exact |
| 4 | **C-eigenvectors.** \|y_μ⟩ are built from \|y_low⟩ = e₈ (pattern {1\|00}) by shifted t's at θ. C-roots y₁ = θ−μ₁₁h, y₂₁ = θ−μ₂₁h, y₂₂ = θ−(μ₂₂−1)h | exact |
| 5 | **Overlap matrix.** ⟨x_μ\|y_ν⟩ is triangular and z-independent, with ⟨x₀\|y₀⟩ = 9h⁴ (y₀ = y_{2\|21}). Its off-diagonal entries are 3h⁴ at (2,4) and (6,7) | exact |
| 6 | **Frames.** P_aux G(z) P_aux⁻¹ = diag(z); P_phys t₁(z) P_phys⁻¹ = Σ z_i ℒ_ii; Ω_L = e₁ᵀ P_pp is the companion vacuum; Ω_zw is the vacuum of `tGPzw` = P_phys t(w) P_phys⁻¹ | exact |
| 7 | **Direct Baxter equation.** ν₁(x)ν₁(x+h)ν₁(x+2h) Q(x−h) − t₁(x) ν₁(x+h)ν₁(x+2h) Q(x) + t₂(x+h) ν₁(x+2h) Q(x+h) − t₃(x+2h) Q(x+2h) = 0, solved by Q_i = z_i^{−u/h} q_i with deg q_i = λ₁−μ_i | 3 points; 4 reps |
| 8 | **Dual Baxter equation** (uniform form). ν₃(x)ν₃(x−h)ν₃(x−2h) Q(x+h) − t₁(x) ν₃(x−h)ν₃(x−2h) Q(x) + t₂(x) ν₃(x−2h) Q(x−h) − t₃(x) Q(x−2h) = 0, solved by Q^j = z_j^{+u/h} q^j with deg q^j = μ_j − λ₃ | 3 points; 4 reps |
| 9 | **QQ relation.** Q_a(u)Q_b(u+h) − Q_b(u)Q_a(u+h) = (1/z_b − 1/z_a) χ₃^{−u/h} F_λ(u) Q^c(u), with F_λ = ∏_{n=λ₂+1}^{λ₁}(u−θ+nh) | 3 points; 4 reps |
| 10 | **Wronskians.** With D_{k₁k₂k₃} = det Q^j(x+k_c h): t₁ = ν₃(x) D_{1,−1,−2}/D_{0,−1,−2}, t₂ = ν₃(x)ν₃(x−h) D_{1,0,−2}/D_{0,−1,−2}, t₃ = ν₃(x)ν₃(x−h)ν₃(x−2h) D_{1,0,−1}/D_{0,−1,−2} | adjoint; 4 reps |
| 11 | **Quantum Weyl formula.** T_λ(x) = ∏_{(r,c)∈λ} ν₃(x+(c−r)h) · det Q^j(x+(λ_k−k+1)h) / det Q^j(x+(1−k)h) | all 19 λ ⊆ 3×3, adjoint |
| 12 | **Degenerate columns.** v(θ) ∥ v(θ−h), and v(θ−3h), …, v(θ) span a 2-plane. Likewise u(θ−2h) ∥ u(θ−h), and u(θ−2h), …, u(θ+h) span a 2-plane | adjoint |
| 13 | **Reduction at θ−2h.** T_λ(θ−2h) = P_λ · det₂[v(θ+(λ₁−2)h), v(θ+(λ₂−3)h)] / det₂[v(θ−2h), v(θ−3h)] if λ₃ = 0 and λ₁ ≤ 2, and 0 otherwise. Only (1), (1,1), (2), (2,1) survive | adjoint, any row pair |
| 14 | **x-side wavefunction.** ⟨x_μ\|Ψ⟩ = N_μ · det₂[v(x₂₁), v(x₁)] · v_k(x₂₂) / (the same at x₀), with N_μ = n_μ h^{k_μ}, n = (1,−2,6,−1,2,−2,−6,6) and k_μ = 5−μ₁₁−μ₂₁−μ₂₂ | adjoint, h = 1 and 1/2 |
| 15 | **y-side wavefunction.** ⟨Ψ\|y_μ⟩/⟨Ψ\|y_low⟩ = N′_μ · det₂[u(y₂₁), u(y₂₂)] · u_k(y₁) / (the same at y_low) | adjoint |
| 16 | **Functional orthogonality.** For each (i,j): −δI₁⟨P^jQ_i⟩ + a⟨P^jQ_i^{[2]}⟩ + b⟨P^{j[−2]}Q_i⟩ = 0, with ⟨f⟩ = Σ_{n=λ₃}^{λ₁} f(θ−nh)/∏_{m≠n}(m−n)h. So det N with rows (1,3), (1,2), (2,3) vanishes for A≠B at the same twist | 4 reps |
| 17 | **Norm formula.** ⟨Ψ_A\|Ψ_B⟩ = (1/⟨x₀\|y₀⟩) · (det N[ket B dual, bra A direct]/c₀) / (Ψ_Q^B(x₀) Φ_Q^A(y₀)), with c₀ = 1/(4h⁶) | (h,χ₃) = (1,1), (1,7/2), (1/2,5/4) |
| 18 | **Mixed twists.** The same formula holds for ⟨Ψ_{A,w}\|Ψ_{B,z}⟩, each state with its own twisted Q's | all 64 entries, ~10⁻⁵¹ |
| 19 | **Diagonal frame.** Ω_zw·P_phys\|Ψ_{B,z}⟩ = (Ω_zw P_phys·y₀) · det N[Q^j_{B,z}, Q_{i,vac,w}] / (c₀⟨x₀\|y₀⟩ Ψ_Q^{B,z}(x₀) Φ_Q^{vac,w}(y₀)) | all 8, ~10⁻⁵³ |

**Reference factors** (rows 17–19):
- Ψ_Q(x₀) = [Q²(θ−2h)Q³(θ−3h) − Q³(θ−2h)Q²(θ−3h)] · Q³(θ−h) — from the ket's dual Q's;
- Φ_Q(y₀) = [Q₁(θ−2h)Q₂(θ) − Q₂(θ−2h)Q₁(θ)] · Q₁(θ−2h) — from the bra's direct Q's.

These are the wavefunctions of rows 14 and 15 evaluated at x₀ and y₀. Of the 81 possible combinations of rows and single factors, they are the only one that gives a state-independent ratio. c₀ is the common coefficient of the four (x₀,y₀) monomials when det N is expanded with the Q's as symbols. The leftover constant 1/⟨x₀|y₀⟩ is the (y₀,x₀) entry of the inverse of the overlap matrix in row 5.

## How it fits together

- **Adjointness.** The bracket is the sum of residues of μ·(…) over all poles of the rational measure μ. Because the Q's are entire, this sum is invariant under shifts. Shifting turns ⟨P·O[Q]⟩ into ⟨φ·O†[P]·Q⟩, with φ = ν₁(x+h)ν₁(x+2h)/(ν₃(x−h)ν₃(x−2h)) and μ = ∏_{n=λ₃}^{λ₁+2}(x−θ+nh)⁻¹.
- **Splitting δt₂.** Write δt₂(y) = a ν₁(y) + b ν₃(y−h). Keeping the a part on the direct side and moving the b part to the dual side gives every term the same measure, ρ = ∏_{n=λ₃}^{λ₁}(x−θ+nh)⁻¹.
- **Orientation of det N.** Pair the **ket's dual** Q's with the **bra's direct** Q's. This is forced by the SoV expansion Σ⟨Ψ_A|y⟩M⟨x|Ψ_B⟩. At equal twists the other orientation agrees on the diagonal; at different twists it fails.
- **Row choices.** 77 of the 84 choices of three rows (i,j) work. The 7 that fail have the same i in every row, the same j in every row, or are the diagonal choice.

## Relation to the old conventions

- **The equations swap roles.** The new dual equation is the old direct one (with ν₁ → ν₃), and the new direct equation, after x → y+h, is the old dual one (with ν₃ → ν₁). Equivalently, each new equation is the old one under h → −h followed by t_a(u) → t_a(u+(a−1)h). The asymptotics swap along with them.
- **Why, structurally (not separately verified).** u + hE_ij = u − hE′_ji with E′ = −Eᵀ, which generates the dual rep λ* = (−λ₃,−λ₂,−λ₁). So the new conventions for λ are the old ones for λ*, up to the twist gauge. For example deg q^j = μ_j − λ₃ = λ*₁ − μ*_j.
- **Proposed relabelling (not applied).** Swap direct ↔ dual and upper ↔ lower indices, optionally also reversing j → 4−j (with z_j → z_{4−j}). With the reversal, the reference factors take the old Q₁·Q₁₂ and Q³·(Q³,Q²) patterns.

## Beyond the adjoint (`Baxter-QQ-Wronskian.wb`)

These were run for [2,1,0], [3,1,0], [3,2,0] and [2,1,1] at h=1, θ=1/3, z=(5/2,7/3,3/5).
- **Holds for all four.** Rows 1, 7–10 and 16.
- **Orthogonality numbers.** For A≠B the 9×3 matrix has rank 2: its smallest singular value is ≲10⁻³⁹, against ≳1 for A=B. The det N off-diagonal/diagonal ratio is ≲10⁻⁴⁵.
- **QQ factor.** F = (u−θ+2h)(u−θ+3h) for [3,1,0], and u−θ+3h for [3,2,0].
- **λ₃ correction.** The −λ₃ in the dual degree only showed up at [2,1,1].
- **Not done beyond the adjoint.** The SoV bases, the wavefunctions and the normalized norm formula.

## Pitfalls

- **`Ordering` and `-Infinity`.** An exactly-zero polynomial has an empty coefficient list, so its `Max` is `-Infinity`, and `Ordering` sorts `-Infinity` after numbers. The clean notebook's `maxc` floors at 0.
- **Restoring stripped exponentials.** Stripping z_j^{θ/h} from the rows needs a consistent restoring factor. An early κ value showed a spurious z₃^{2θ/h} for this reason.
- **h = 1 in the working notebook.** Several cells there hard-code h=1; the clean notebook does not.
- **Labelling of det N.** The working notebook's `detNN[A,B]` paired the dual Q's of A with the direct Q's of B. That was harmless at equal twists; read det N as (ket dual, bra direct).

## Open

- **Other reps.** The SoV bases and the normalized scalar product for [3,1,0] and beyond, including how c₀, ⟨x₀|y₀⟩ and the reference points depend on λ.
- **The Casoratian.** A closed form for D_{1,0,−1}. Its roots in x−θ are {0,−1,−2} for [2,1,0], {0..−3} for [3,1,0], {0,−1,−1,−2,−3} for [3,2,0] and {−2} for [2,1,1]; its degree is |λ|−3λ₃.
- **The paper.** Updating `paper/xxxCG.tex`, which still uses the old conventions.
